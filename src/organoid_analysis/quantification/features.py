"""Physical morphology and raw fluorescence readouts; no inferred cell fractions."""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
from skimage.measure import marching_cubes, mesh_surface_area

from organoid_analysis.microscopy_io.tiff_contract import Sample
from organoid_analysis.segmentation.watershed_instances import SegmentationResult

from . import surface_crofton
from .labels import bbox_touches_volume_boundary

# Consistency constant rescaling MAD into a sigma-equivalent robust scale
# estimator under a Gaussian assumption (Rousseeuw & Croux, 1993, J. Am. Stat.
# Assoc. 88(424):1273-1283, https://doi.org/10.1080/01621459.1993.10476408).
# See docs/ALGORITHM_DECISIONS.md D5. `multilevel3d/qc.py` uses the reciprocal
# form (0.67448975 x ...) of the same statistical concept, at a different
# rounding precision -- 0.67448975 is Phi^-1(0.75) to 8 decimals, not bit-
# identical to 1/1.4826 (they differ by ~1.5ppm, immaterial in practice).
MAD_TO_SIGMA = 1.4826

# Versioned identity of the estimator that defines every exported
# ``surface_area_um2`` (and therefore ``sphericity``). Any change to the
# weights, the stencil rule, the domain constants or the transition counting
# must bump ``method_version`` and repeat the surface V&V, so historical values
# are never silently reinterpreted.
#
# Qualified 2026-09-12 on an untouched confirmation set of 96 in-domain cases:
# worst absolute relative area error 0.951% against the <5% criterion of
# docs/SCIENTIFIC_SPEC.md section 9 (evidence:
# docs/evidence/2026-09-12-surface-crofton-v3). The qualification is
# DOMAIN-RESTRICTED; see ``domain`` below and ``surface_crofton.DOMAIN_SCOPE``.
# Outside that domain -- in particular for surfaces with dihedral creases, at
# any resolution -- this estimator is NOT QUALIFIED and the section 9 criterion
# is not met. ``geometry()`` reports the domain variables and an
# ``in_domain`` flag per object rather than silently gating.
SURFACE_AREA_METHOD: dict[str, str | int | float] = {
    "estimator": "crofton_minimax_sym",
    "method_version": surface_crofton.METHOD_NAME,
    "implementation_version": "padded_edt_v1",
    "implementation": ("organoid_analysis.quantification.surface_crofton: weighted lattice "
                       "transition counts (discrete Cauchy/Crofton); weights from the "
                       "orbit-symmetry-reduced minimax LP over plane orientations"),
    "input": "binary measurement support, zero-padded by the stencil radius",
    "spacing": "native physical ZYX voxel spacing",
    "stencil_radius_rule": "per object, minimises plane-response deviation + grazing deficit",
    "stencil_radius_candidates": "1,2,3,4,5",
    "smoothing": "none",
    "domain": (f"rho_in >= {surface_crofton.DOMAIN_RHO_IN_MIN} and "
               f"anisotropy <= {surface_crofton.DOMAIN_ANISO_MAX}, "
               f"{surface_crofton.DOMAIN_SCOPE}"),
    "not_qualified": surface_crofton.DOMAIN_NOT_QUALIFIED,
    "evidence": "docs/evidence/2026-09-12-surface-crofton-v3",
    "implementation_evidence": "docs/evidence/2026-09-13-morphology-architecture/revalidation.json",
}

# The superseded estimator, preserved under its own identity together with the
# evidence that it FAILS the section 9 criterion: it errs by a median of 12.5%
# and a worst of 21.9% on the very cases the qualified estimator handles to
# 0.33% median, and its error does not shrink with resolution. Retained so
# historical values remain reproducible and comparable, and reachable through
# ``legacy_surface_area()``; it defines no exported column.
LEGACY_SURFACE_AREA_METHOD: dict[str, str | int | float] = {
    "estimator": "marching_cubes_binary",
    "method_version": "marching_cubes_binary_lewiner_v1",
    "implementation": "skimage.measure.marching_cubes (Lewiner) + skimage.measure.mesh_surface_area",
    "input": "binary measurement support, zero-padded by one voxel",
    "isosurface_level": 0.5,
    "spacing": "native physical ZYX voxel spacing",
    "step_size": 1,
    "smoothing": "none",
    "status": "SUPERSEDED 2026-09-12; FAILS the docs/SCIENTIFIC_SPEC.md section 9 area criterion",
    "evidence": "docs/evidence/2026-09-11-measurement-vv",
}

# Sphericity cannot exceed 1 for a continuous solid; discretized values above
# this limit are flagged for review, never clipped (SCIENTIFIC_SPEC section 10,
# item 3). Shared by the classical and Web/mask-feature routes. Heuristic
# review limit; see docs/PARAMETERS.md.
SPHERICITY_REVIEW_LIMIT = 1.05

# Shared output contract: adapters must preserve these fields together. The
# boolean gate covers resolution/anisotropy only, never segmentation validity.
SURFACE_METADATA_COLUMNS = [
    "surface_rho_in", "surface_anisotropy", "surface_stencil_radius",
    "surface_apriori_rel_bound", "surface_in_qualified_domain", "surface_domain_flags",
    "surface_estimator", "surface_method_version", "surface_implementation_version",
    "surface_weights_origin", "surface_weights_evidence_bearing", "surface_evidence", "surface_implementation_evidence", "surface_qualification_scope",
    "sphericity_in_qualified_domain", "surface_to_volume_in_qualified_domain",
]
TOPOLOGY_COLUMNS = ["filled_void_voxels", "filled_void_components", "open_cavity_suspected", "topology_flags"]
DERIVED_GEOMETRY_COLUMNS = ["elongation", "prolate_ratio", "oblate_ratio", "surface_to_volume_ratio_um_inv"]

GEOMETRY_COLUMNS = ["organoid_id", "original_label_id", "segmented_voxels", "envelope_voxels",
                    "segmented_volume_um3", "volume_um3", "surface_area_um2", "sphericity",
                    "equivalent_diameter_um", "enclosed_void_fraction", *TOPOLOGY_COLUMNS, "centroid_z_um", "centroid_y_um",
                    "centroid_x_um", "extent_z_um", "extent_y_um", "extent_x_um", "principal_axis_major_um",
                    "principal_axis_intermediate_um", "principal_axis_minor_um", "axis_ratio_minor_to_major", "n_z_slices",
                    "measurement_basis", *DERIVED_GEOMETRY_COLUMNS, *SURFACE_METADATA_COLUMNS,
                    "touches_border", "morphology_eligible", "morphology_flags"]
MARKER_COLUMNS = ["background_voxels"] + [f"{marker}_{field}" for marker in ["calcein", "pi"]
                  for field in ["mean_raw", "background_median", "background_noise_mad", "mean_bg_corrected",
                                "integrated_bg_corrected", "saturated_fraction"]] + ["viability_measurement_eligible", "viability_measurement_flags"]


def outer_envelope(mask: np.ndarray) -> np.ndarray:
    # binary_fill_holes/slicing on an ndarray input always return an ndarray;
    # np.asarray only fixes scipy's untyped (Any) stub, no value/dtype change.
    return np.asarray(ndi.binary_fill_holes(np.pad(mask, 1))[1:-1, 1:-1, 1:-1])


def surface_mesh(
    mask: np.ndarray,
    spacing: tuple[float, float, float],
    origin_zyx: tuple[float, float, float] = (0, 0, 0),
    step_size: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    # Marching cubes: Lorensen & Cline (1987), ACM SIGGRAPH Computer Graphics
    # 21(4):163-169, https://doi.org/10.1145/37401.37422. See
    # docs/ALGORITHM_DECISIONS.md D3 for the level/spacing/padding rationale.
    # Padding closes the mesh at the crop boundary; clipped image objects remain QC-excluded.
    vertices, faces, _, _ = marching_cubes(np.pad(mask.astype(np.float32), 1), level=0.5,
                                          spacing=spacing, step_size=step_size, allow_degenerate=False)
    vertices += np.asarray(origin_zyx) * spacing - np.asarray(spacing)
    return vertices, faces


def legacy_surface_area(
    mask: np.ndarray,
    spacing: tuple[float, float, float],
    origin_zyx: tuple[float, float, float] = (0, 0, 0),
) -> float:
    """Surface area under the superseded ``marching_cubes_binary_lewiner_v1``.

    Preserved so historical ``surface_area_um2`` values stay reproducible and
    so the evidence that this estimator fails the section 9 criterion can be
    re-executed. It defines no exported column; ``geometry()`` uses the
    qualified estimator. See ``LEGACY_SURFACE_AREA_METHOD``.
    """
    return float(mesh_surface_area(*surface_mesh(mask, spacing, origin_zyx)))


def geometry(
    mask: np.ndarray,
    spacing: tuple[float, float, float],
    origin_zyx: tuple[float, float, float] = (0, 0, 0),
    *,
    fill_holes: bool = True,
) -> tuple[dict, tuple]:
    """Return physical geometry for a nonempty 3D mask.

    The established organoid pipeline measures the filled outer envelope by
    default.  Instance-level multilevel analysis can opt out so its reported
    voxel count, volume, surface area and sphericity describe the same raw
    label voxels.
    """
    if mask.ndim != 3 or not mask.any():
        raise ValueError("Geometry needs a nonempty 3D instance mask")
    # A distinct name from the `spacing` parameter (rather than reassigning
    # it to a different type) -- same values, only the static type changes.
    spacing_arr = np.asarray(spacing, float)
    if spacing_arr.shape != (3,) or not np.isfinite(spacing_arr).all() or (spacing_arr <= 0).any():
        raise ValueError("Geometry needs three positive finite spacings")
    envelope = outer_envelope(mask) if fill_holes else mask.astype(bool, copy=False)
    segmented = int(mask.sum())
    count = int(envelope.sum())
    padded = np.pad(np.asarray(mask, bool), 1)
    background_labels, _ = ndi.label(~padded, structure=ndi.generate_binary_structure(3, 1))
    border_labels = np.unique(np.concatenate([
        background_labels[0].ravel(), background_labels[-1].ravel(),
        background_labels[:, 0, :].ravel(), background_labels[:, -1, :].ravel(),
        background_labels[:, :, 0].ravel(), background_labels[:, :, -1].ravel(),
    ]))
    foreground_coords = np.argwhere(mask)
    bbox_slices = tuple(
        slice(int(axis.min()), int(axis.max()) + 1) for axis in foreground_coords.T
    )
    interior_background = background_labels[1:-1, 1:-1, 1:-1][bbox_slices]
    local_mask = np.asarray(mask, bool)[bbox_slices]
    interior_labels = np.unique(interior_background[~local_mask])
    interior_labels = interior_labels[interior_labels > 0]
    open_cavity = bool(np.any(np.isin(interior_labels, border_labels)))
    filled_void_components = int(np.count_nonzero(~np.isin(interior_labels, border_labels)))
    topology_flags = ["open_cavity_suspected"] if open_cavity else []
    volume = float(count * np.prod(spacing_arr))
    # Surface area comes from the qualified Crofton estimator, not from the
    # mesh: the mesh is still built and returned for visualisation and PLY
    # export, but a marching-cubes mesh of a binarised field measures the area
    # of a staircase, which is the defect SURFACE_AREA_METHOD records.
    vertices, faces = surface_mesh(envelope, tuple(spacing_arr), origin_zyx)
    surface = surface_crofton.measure(envelope, tuple(spacing_arr))
    area = float(surface["surface_area"])
    qualified, domain_flags = surface_crofton.in_domain(surface)
    sphericity = float(np.cbrt(np.pi) * (6.0 * volume) ** (2 / 3) / area)
    points = np.argwhere(envelope)
    center = (points.mean(axis=0) + origin_zyx) * spacing_arr
    extent = (points.max(axis=0) - points.min(axis=0) + 1) * spacing_arr
    centered = (points-points.mean(axis=0))*spacing_arr
    # Include each voxel's intrinsic second moment; axes describe a moment-equivalent ellipsoid.
    covariance = centered.T @ centered / len(points) + np.diag(spacing_arr**2/12)
    lengths = 2*np.sqrt(5*np.linalg.eigvalsh(covariance)[::-1])
    values = {"segmented_voxels": segmented, "envelope_voxels": count,
              "segmented_volume_um3": float(segmented * np.prod(spacing_arr)), "volume_um3": volume,
              "surface_area_um2": area, "sphericity": sphericity,
              "equivalent_diameter_um": float(np.cbrt(6 * volume / np.pi)),
              "enclosed_void_fraction": float((count - segmented) / count),
              "filled_void_voxels": int(count - segmented),
              "filled_void_components": filled_void_components,
              "open_cavity_suspected": open_cavity,
              "topology_flags": ";".join(topology_flags),
              "principal_axis_major_um": float(lengths[0]), "principal_axis_intermediate_um": float(lengths[1]),
              "principal_axis_minor_um": float(lengths[2]), "axis_ratio_minor_to_major": float(lengths[2]/lengths[0]),
              "elongation": float(1.0 - lengths[2] / lengths[0]),
              "prolate_ratio": float(lengths[0] / lengths[1]),
              "oblate_ratio": float(lengths[1] / lengths[2]),
              "surface_to_volume_ratio_um_inv": area / volume,
              "measurement_basis": "filled_envelope" if fill_holes else "raw_label",
              "n_z_slices": int(np.any(envelope, axis=(1, 2)).sum()),
              # Applicability-domain variables of the surface estimator, exported
              # per object so a consumer can tell whether a given
              # surface_area_um2/sphericity pair is inside the qualified domain.
              # ``surface_in_qualified_domain`` covers the machine-checkable gate
              # only; smoothness (SURFACE_AREA_METHOD["domain"]) cannot be
              # checked from a mask and remains the caller's responsibility.
              "surface_rho_in": float(surface["rho_in"]),
              "surface_anisotropy": float(surface["anisotropy"]),
              "surface_stencil_radius": int(surface["stencil_m"]),
              "surface_apriori_rel_bound": float(surface["apriori_bound"]),
              "surface_in_qualified_domain": bool(qualified),
              "surface_domain_flags": ";".join(domain_flags),
              "surface_estimator": SURFACE_AREA_METHOD["estimator"],
              "surface_method_version": surface["method"],
              "surface_implementation_version": SURFACE_AREA_METHOD["implementation_version"],
              "surface_weights_origin": surface["weights_origin"],
              "surface_weights_evidence_bearing": bool(surface["weights_evidence_bearing"]),
              "surface_evidence": SURFACE_AREA_METHOD["evidence"],
              "surface_implementation_evidence": SURFACE_AREA_METHOD["implementation_evidence"],
              "surface_qualification_scope": "numerical_resolution_and_anisotropy_only",
              "sphericity_in_qualified_domain": bool(qualified),
              "surface_to_volume_in_qualified_domain": bool(qualified)}
    for axis, value, width in zip("zyx", center, extent):
        values[f"centroid_{axis}_um"] = float(value)
        values[f"extent_{axis}_um"] = float(width)
    return values, (vertices, faces)


def write_ply(path: Path, vertices_zyx: np.ndarray, faces: np.ndarray) -> None:
    with path.open("w", encoding="ascii") as handle:
        handle.write("ply\nformat ascii 1.0\ncomment coordinates XYZ in micrometres\n")
        handle.write(f"element vertex {len(vertices_zyx)}\nproperty float x\nproperty float y\nproperty float z\n")
        handle.write(f"element face {len(faces)}\nproperty list uchar int vertex_indices\nend_header\n")
        np.savetxt(handle, vertices_zyx[:, ::-1], fmt="%.6f")
        np.savetxt(handle, np.column_stack([np.full(len(faces), 3), faces]), fmt="%d")


def marker_measurements(labels: np.ndarray, object_id: int, bbox: tuple,
                        channels: dict, spacing: tuple, cfg: dict) -> dict:
    fields = ["mean_raw", "background_median", "background_noise_mad", "mean_bg_corrected",
              "integrated_bg_corrected", "saturated_fraction"]
    # Genuinely heterogeneous: most entries are float measurements (nan until
    # filled in), but background_voxels is int, the eligibility flag is bool,
    # and the flags summary is str -- this dict becomes one output row, so
    # its real value type spans all four, not just float.
    result: dict[str, float | int | bool | str] = {
        f"{marker}_{field}": np.nan
        for marker in channels
        for field in fields
    }
    missing = [marker for marker, channel in channels.items() if channel is None]
    # Structure-only runs do not need a per-object distance transform. Preserve
    # the same ineligible/missing-channel semantics while avoiding the dominant
    # marker-measurement cost entirely.
    if len(missing) == len(channels):
        result.update(background_voxels=0, viability_measurement_eligible=False,
                      viability_measurement_flags=";".join(f"missing_{m}_channel" for m in missing))
        return result

    pad = np.ceil(cfg["background_outer_um"] / np.asarray(spacing)).astype(int) + 1
    sl = tuple(slice(max(0, box.start - p), min(length, box.stop + p)) for box, p, length in zip(bbox, pad, labels.shape))
    local_labels = labels[sl]
    roi = outer_envelope(local_labels == object_id)
    outside_distance = ndi.distance_transform_edt(~roi, sampling=spacing)
    shell = (outside_distance > cfg["background_inner_um"]) & (outside_distance <= cfg["background_outer_um"]) & (local_labels == 0)
    background_n = int(shell.sum())
    result["background_voxels"] = background_n
    flags = []
    if background_n < cfg["min_background_voxels"]:
        flags.append("insufficient_local_background")
    for marker, channel in channels.items():
        if channel is None:
            flags.append(f"missing_{marker}_channel")
            continue
        foreground_values = channel[sl][roi].astype(np.float64)
        result[f"{marker}_mean_raw"] = float(foreground_values.mean())
        limit = cfg[f"{marker}_saturation_value"]
        if limit is None and np.issubdtype(channel.dtype, np.integer):
            limit = float(np.iinfo(channel.dtype).max)
        if limit is None:
            flags.append(f"{marker}_saturation_limit_unknown")
        else:
            saturated = float(np.mean(foreground_values >= limit))
            result[f"{marker}_saturated_fraction"] = saturated
            if saturated > cfg["max_saturated_fraction"]:
                flags.append(f"{marker}_saturated")
        if background_n >= cfg["min_background_voxels"]:
            background = channel[sl][shell].astype(np.float64)
            median = float(np.median(background))
            noise = float(MAD_TO_SIGMA * np.median(np.abs(background - median)))
            # Keep negative corrected means. Per-voxel clipping would bias dim objects upward.
            corrected = float(foreground_values.mean() - median)
            result[f"{marker}_background_median"] = median
            result[f"{marker}_background_noise_mad"] = noise
            result[f"{marker}_mean_bg_corrected"] = corrected
            result[f"{marker}_integrated_bg_corrected"] = float(corrected * roi.sum())
    result["viability_measurement_eligible"] = not flags
    result["viability_measurement_flags"] = ";".join(flags)
    return result


def measure_instances(
    sample: Sample, segmentation: SegmentationResult, cfg: dict, mesh_dir: Path | None = None
) -> tuple[list[dict], list[tuple]]:
    labels = segmentation.labels
    rows: list[dict] = []
    preview_meshes: list[tuple[int, np.ndarray, np.ndarray]] = []
    if mesh_dir:
        mesh_dir.mkdir(parents=True, exist_ok=True)
    for object_id, bbox in enumerate(ndi.find_objects(labels), start=1):
        if bbox is None:
            continue
        mask = labels[bbox] == object_id
        origin = tuple(s.start for s in bbox)
        metrics, mesh = geometry(mask, sample.spacing, origin)
        touches = bbox_touches_volume_boundary(bbox, labels.shape)
        flags = []
        if touches:
            flags.append("border_truncated")
        if metrics["n_z_slices"] < cfg["quality"]["min_z_slices"]:
            flags.append("too_few_z_slices")
        if metrics["volume_um3"] < cfg["segmentation"]["min_volume_um3"]:
            flags.append("below_minimum_volume")
        maximum = cfg["quality"]["max_volume_um3"]
        if maximum is not None and metrics["volume_um3"] > maximum:
            flags.append("above_maximum_volume")
        if metrics["sphericity"] > SPHERICITY_REVIEW_LIMIT:
            flags.append("sphericity_above_geometric_range")
        if "excess_foreground_review_segmentation" in segmentation.flags:
            flags.append("excess_foreground_review_segmentation")
        envelope = outer_envelope(mask)
        if ((labels[bbox] != object_id) & (labels[bbox] != 0) & envelope).any():
            flags.append("encloses_other_instance")
        excluded = [flag for flag in flags if flag != "border_truncated" or cfg["quality"]["exclude_border"]]
        metrics.update({"organoid_id": object_id, "original_label_id": segmentation.original_ids[object_id],
                        "touches_border": touches, "morphology_eligible": not excluded,
                        "morphology_flags": ";".join(flags)})
        metrics.update(marker_measurements(labels, object_id, bbox, {"calcein": sample.calcein, "pi": sample.pi}, sample.spacing, cfg["quality"]))
        rows.append(metrics)
        if mesh_dir:
            write_ply(mesh_dir / f"organoid_{object_id:04d}.ply", *mesh)
        if len(preview_meshes) < cfg["report"]["max_meshes_in_preview"]:
            # Only the preview uses a coarser mesh. Measurements and PLYs use step_size=1.
            preview_meshes.append((object_id, *surface_mesh(envelope, sample.spacing, origin, step_size=2)))
    return rows, preview_meshes
