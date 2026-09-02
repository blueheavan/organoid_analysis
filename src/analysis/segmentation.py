"""Classical 3D instance segmentation, with parameters in physical units."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu
from skimage.morphology import h_maxima
from skimage.segmentation import watershed

from .labels import compact_instance_labels


@dataclass
class SegmentationResult:
    labels: np.ndarray
    threshold: float | None
    foreground_fraction: float
    removed_small_instances: int
    flags: list[str]
    original_ids: dict[int, int]


QC_COLUMNS = [
    "reference_method", "primary_instances", "reference_instances",
    "primary_median_z_extent_um", "reference_median_z_extent_um",
    "object_count_difference_fraction", "z_extent_ratio", "requires_review",
    "qc_flags",
]


def physical_ball(radius_um: float, spacing: tuple) -> np.ndarray:
    if radius_um <= 0:
        return np.ones((1, 1, 1), bool)
    half = np.ceil(radius_um / np.asarray(spacing)).astype(int)
    if np.prod(2 * half + 1) > 5_000_000:
        raise ValueError("Morphological footprint is too large; reduce the physical radius")
    grids = np.ogrid[tuple(slice(-h, h + 1) for h in half)]
    squared = sum((g * s) ** 2 for g, s in zip(grids, spacing))
    return squared <= radius_um ** 2 + 1e-9


def _drop_small(labels: np.ndarray, min_voxels: int) -> tuple[np.ndarray, int]:
    """Remove small dense labels and compact survivors with one LUT pass.

    ``ndi.label`` and watershed produce IDs bounded by ``labels.max()``. A
    direct lookup table avoids ``np.unique(..., return_inverse=True)``, whose
    full-volume int64 inverse array is expensive for large 3D stacks.
    """
    counts = np.bincount(labels.ravel())
    small = (counts > 0) & (counts < min_voxels)
    if len(small):
        small[0] = False
    removed = int(small.sum())
    keep = np.flatnonzero((counts >= min_voxels) & (np.arange(len(counts)) != 0))
    lookup = np.zeros(len(counts), dtype=np.uint32)
    lookup[keep] = np.arange(1, len(keep) + 1, dtype=np.uint32)
    return lookup[labels], removed


def watershed_instances(mask: np.ndarray, spacing: tuple, cfg: dict) -> np.ndarray:
    components, count = ndi.label(mask)
    if not count or not cfg["split_touching"]:
        return components.astype(np.uint32)
    distance = ndi.distance_transform_edt(np.pad(mask, 1), sampling=spacing)[1:-1, 1:-1, 1:-1]
    # h-maxima prominence is a distance in micrometres, not a voxel count.
    peaks = h_maxima(distance, cfg["seed_h_um"]) & mask
    peak_labels, _ = ndi.label(peaks)
    coords = []
    for peak_id, sl in enumerate(ndi.find_objects(peak_labels), start=1):
        if sl is None:
            continue
        values = np.where(peak_labels[sl] == peak_id, distance[sl], -1)
        point = np.array(np.unravel_index(values.argmax(), values.shape)) + [s.start for s in sl]
        coords.append(point)
    # Suppress nearby seeds in the same connected component using physical distance.
    selected: dict[int, list[np.ndarray]] = {}
    for point in sorted(coords, key=lambda p: float(distance[tuple(p)]), reverse=True):
        component = int(components[tuple(point)])
        previous = selected.setdefault(component, [])
        if all(np.linalg.norm((point - other) * spacing) >= cfg["seed_min_distance_um"] for other in previous):
            previous.append(point)
    component_boxes = ndi.find_objects(components)
    for component in range(1, count + 1):
        if not selected.get(component):
            sl = component_boxes[component - 1]
            values = np.where(components[sl] == component, distance[sl], -1)
            point = np.array(np.unravel_index(values.argmax(), values.shape)) + [s.start for s in sl]
            selected[component] = [point]
    markers = np.zeros(mask.shape, np.int32)
    marker = 0
    for component in sorted(selected):
        for point in selected[component]:
            marker += 1
            markers[tuple(point)] = marker
    return watershed(-distance, markers=markers, mask=mask, connectivity=1).astype(np.uint32)


def instance_qc_summary(labels: np.ndarray, spacing: tuple) -> dict:
    """Summarize instance count and Z extent for a mask comparison.

    These are deliberately coarse, ground-truth-free checks. They can reveal a
    severe split/fragmentation discrepancy between two methods, but they do not
    establish which method is biologically correct.
    """
    labels = np.asarray(labels)
    if (labels.ndim != 3 or not np.issubdtype(labels.dtype, np.integer)
            or (np.issubdtype(labels.dtype, np.signedinteger) and np.any(labels < 0))):
        raise ValueError("Instance QC requires a 3D nonnegative integer label mask")
    spacing = np.asarray(spacing, dtype=float)
    if spacing.shape != (3,) or not np.isfinite(spacing).all() or np.any(spacing <= 0):
        raise ValueError("Instance QC spacing must contain three positive finite values")
    ids = np.unique(labels)
    ids = ids[ids != 0]
    if not len(ids):
        return {"instances": 0, "median_z_extent_um": 0.0}

    # Common pipeline masks are compact. Avoid creating a volume-sized inverse
    # array in that case; for sparse, very large IDs compact once so
    # ndi.find_objects cannot allocate a list up to an untrusted max label.
    if int(ids[-1]) <= max(100_000, 4 * len(ids)):
        boxes = ndi.find_objects(labels)
        box_by_id = {object_id: boxes[int(object_id) - 1] for object_id in ids}
        work = labels
    else:
        work, source_id_by_compact_id = compact_instance_labels(labels)
        compact_id_by_source_id = {source_id: compact_id for compact_id, source_id in source_id_by_compact_id.items()}
        boxes = ndi.find_objects(work)
        box_by_id = {object_id: boxes[compact_id_by_source_id[int(object_id)] - 1] for object_id in ids}
    z_extents = []
    for object_id in ids:
        bbox = box_by_id[object_id]
        if bbox is None:
            continue
        local_id = int(object_id) if work is labels else compact_id_by_source_id[int(object_id)]
        occupied_z = np.any(work[bbox] == local_id, axis=(1, 2))
        z_extents.append(float(occupied_z.sum() * spacing[0]))
    return {
        "instances": int(len(ids)),
        "median_z_extent_um": float(np.median(z_extents)) if z_extents else 0.0,
    }


def compare_instance_qc(primary: np.ndarray, reference: np.ndarray, spacing: tuple,
                        count_difference_threshold: float, min_z_extent_ratio: float,
                        *, primary_summary: dict | None = None,
                        reference_summary: dict | None = None) -> dict:
    """Compare primary and reference masks without selecting either automatically.

    ``reference`` is an independent segmentation of the same structural stack.
    Large count disagreement and a short primary median Z extent are common
    symptoms of over-splitting or Z-fragmentation. Both are review prompts, not
    accuracy measurements; only a registered annotated mask can establish that.
    """
    primary_summary = primary_summary or instance_qc_summary(primary, spacing)
    reference_summary = reference_summary or instance_qc_summary(reference, spacing)
    count_difference = abs(primary_summary["instances"] - reference_summary["instances"]) / max(
        reference_summary["instances"], 1
    )
    if reference_summary["median_z_extent_um"] > 0:
        z_extent_ratio = primary_summary["median_z_extent_um"] / reference_summary["median_z_extent_um"]
    else:
        z_extent_ratio = np.nan
    flags = []
    if count_difference > count_difference_threshold:
        flags.append("object_count_disagreement_review")
    if np.isfinite(z_extent_ratio) and z_extent_ratio < min_z_extent_ratio:
        flags.append("primary_z_fragmentation_review")
    if primary_summary["instances"] and not reference_summary["instances"]:
        flags.append("reference_no_instances_review")
    return {
        "primary_instances": primary_summary["instances"],
        "reference_instances": reference_summary["instances"],
        "primary_median_z_extent_um": primary_summary["median_z_extent_um"],
        "reference_median_z_extent_um": reference_summary["median_z_extent_um"],
        "object_count_difference_fraction": float(count_difference),
        "z_extent_ratio": float(z_extent_ratio),
        "requires_review": bool(flags),
        "qc_flags": ";".join(flags),
    }


def segment(structure: np.ndarray, spacing: tuple, cfg: dict,
            probability: np.ndarray | None = None, imported_labels: np.ndarray | None = None) -> SegmentationResult:
    method = cfg["method"]
    threshold = None
    flags = []
    if method == "labels":
        if imported_labels is None or not np.issubdtype(imported_labels.dtype, np.integer):
            raise ValueError("Imported labels must be a 3D nonnegative integer instance-label TIFF")
        if imported_labels.shape != structure.shape or (imported_labels < 0).any():
            raise ValueError("Imported labels must match the image shape and contain nonnegative integers")
        labels, mapping = compact_instance_labels(imported_labels)
        for object_id, sl in enumerate(ndi.find_objects(labels), start=1):
            if sl is not None and ndi.label(labels[sl] == object_id)[1] != 1:
                raise ValueError(f"Imported label {mapping[object_id]} has disconnected components; fix the instance mask")
        removed = 0
    else:
        if method == "probability":
            if probability is None or probability.shape != structure.shape:
                raise ValueError("A foreground probability volume matching the image is required")
            if not np.isfinite(probability).all() or probability.min() < 0 or probability.max() > 1:
                raise ValueError("Foreground probabilities must be finite values in [0,1], not logits or 0-255 intensities")
            threshold = float(cfg["probability_threshold"])
            mask = probability >= threshold
        else:
            image = structure.astype(np.float32)
            if cfg["background_sigma_um"] > 0:
                background = ndi.gaussian_filter(image, np.asarray([cfg["background_sigma_um"]] * 3) / spacing)
                image = image - background
            if cfg["polarity"] == "dark":
                image = -image
                flags.append("exploratory_dark_foreground_segmentation")
            if cfg["gaussian_sigma_um"] > 0:
                image = ndi.gaussian_filter(image, np.asarray([cfg["gaussian_sigma_um"]] * 3) / spacing)
            if float(image.max()) == float(image.min()):
                mask = np.zeros(image.shape, bool)
                flags.append("uniform_structure_image")
            else:
                threshold = float(threshold_otsu(image)) if cfg["threshold"] == "otsu" else float(cfg["threshold"])
                mask = image > threshold
        if cfg["closing_radius_um"] > 0 and mask.any():
            footprint = physical_ball(cfg["closing_radius_um"], spacing)
            pads = tuple((n // 2 + 1, n // 2 + 1) for n in footprint.shape)
            extended = np.pad(mask, pads)
            extended = ndi.binary_closing(extended, structure=footprint)
            mask = extended[tuple(slice(a, -b) for a, b in pads)]
        if cfg["fill_enclosed_holes"]:
            mask = ndi.binary_fill_holes(np.pad(mask, 1))[1:-1, 1:-1, 1:-1]
        components, _ = ndi.label(mask)
        min_voxels = max(1, int(np.ceil(cfg["min_volume_um3"] / np.prod(spacing))))
        components, removed_before = _drop_small(components, min_voxels)
        labels = watershed_instances(components > 0, spacing, cfg)
        labels, removed_after = _drop_small(labels, min_voxels)
        removed = removed_before + removed_after
        mapping = {int(i): int(i) for i in np.unique(labels) if i}
    fraction = float(np.count_nonzero(labels) / labels.size)
    if fraction > cfg["max_foreground_fraction"]:
        flags.append("excess_foreground_review_segmentation")
    if not labels.any():
        flags.append("no_organoids_detected")
    if max(spacing) / min(spacing) > 5:
        flags.append("high_voxel_anisotropy")
    return SegmentationResult(labels, threshold, fraction, removed, flags, mapping)
