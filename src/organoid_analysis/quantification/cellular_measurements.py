"""Cell-level N:C pairing QC, morphology, and spatial features."""
from __future__ import annotations

from dataclasses import dataclass
from typing import TypedDict

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching
from scipy.spatial import cKDTree

from .features import DERIVED_GEOMETRY_COLUMNS, SURFACE_METADATA_COLUMNS, geometry, outer_envelope
from .labels import (
    bbox_touches_volume_boundary,
    compact_instance_labels,
    relabel_from_source_ids,
)
from .measurement_policy import measurement_policy

PAIR_COLUMNS = [
    "cell_id",
    "original_cell_id",
    "nucleus_id",
    "original_nucleus_id",
    "candidate_nucleus_id",
    "kept",
    "reject_reason",
    "cell_volume_um3",
    "nucleus_volume_um3",
    "overlap_volume_um3",
    "nucleus_containment_fraction",
    "nc_ratio",
]

GEOMETRY_COLUMNS = [
    "cell_id", "cell_volume_um3", "cell_envelope_volume_um3", "surface_area_um2",
    "sphericity", "equivalent_diameter_um", "centroid_z_um", "centroid_y_um",
    "centroid_x_um", "principal_axis_major_um", "principal_axis_intermediate_um",
    "principal_axis_minor_um", "axis_ratio_minor_to_major", "n_z_slices",
    "enclosed_void_fraction", "touches_border", "nucleus_volume_um3", "nc_ratio",
    "nucleus_centroid_z_um", "nucleus_centroid_y_um", "nucleus_centroid_x_um",
    "cell_to_nucleus_centroid_um", "nucleus_centroid_to_cell_border_um",
    "volume_um3", "measurement_basis", *DERIVED_GEOMETRY_COLUMNS, *SURFACE_METADATA_COLUMNS,
]


@dataclass
class CellAnalysisResult:
    cell_labels: np.ndarray
    nucleus_labels: np.ndarray
    features: pd.DataFrame
    pairing: pd.DataFrame
    summary: dict


def _validate_labels(name: str, labels: np.ndarray) -> None:
    if not np.issubdtype(labels.dtype, np.integer) or (labels < 0).any():
        raise ValueError(f"{name} must contain nonnegative integer instance labels")


def _validate_spacing(spacing: tuple) -> tuple[float, float, float]:
    spacing = tuple(float(value) for value in spacing)
    if len(spacing) != 3 or not np.isfinite(spacing).all() or min(spacing) <= 0:
        raise ValueError("spacing must contain three positive finite ZYX values")
    return spacing


def _validate_inputs(cell_masks: np.ndarray, nuclei_masks: np.ndarray, spacing: tuple) -> tuple[float, float, float]:
    if cell_masks.shape != nuclei_masks.shape or cell_masks.ndim != 3:
        raise ValueError("Cell and nucleus masks must be matching 3D ZYX volumes")
    _validate_labels("cell_masks", cell_masks)
    _validate_labels("nuclei_masks", nuclei_masks)
    return _validate_spacing(spacing)


def _validate_single(name: str, labels: np.ndarray, spacing: tuple) -> tuple[float, float, float]:
    if labels.ndim != 3:
        raise ValueError(f"{name} must be a 3D ZYX volume")
    _validate_labels(name, labels)
    return _validate_spacing(spacing)


def _validate_qc(min_cell_volume_um3: float, min_nucleus_volume_um3: float,
                 max_nc_ratio: float, min_nucleus_containment: float) -> None:
    values = (min_cell_volume_um3, min_nucleus_volume_um3, max_nc_ratio, min_nucleus_containment)
    if not np.isfinite(values).all():
        raise ValueError("Cell QC thresholds must be finite")
    if min_cell_volume_um3 < 0 or min_nucleus_volume_um3 < 0 or max_nc_ratio <= 0:
        raise ValueError("Volume thresholds must be nonnegative and max_nc_ratio must be positive")
    if not 0 <= min_nucleus_containment <= 1:
        raise ValueError("min_nucleus_containment must be in [0, 1]")


def _distance_to_mask_surface(mask: np.ndarray, point_zyx: np.ndarray, spacing: tuple) -> float:
    """Exact point distance to the exposed faces of an axis-aligned voxel mask."""
    point_um = point_zyx * spacing
    best_squared = np.inf
    for axis in range(3):
        for direction in (-1, 1):
            neighbor = np.zeros_like(mask)
            source = [slice(None)] * 3
            target = [slice(None)] * 3
            if direction < 0:
                source[axis] = slice(0, -1)
                target[axis] = slice(1, None)
            else:
                source[axis] = slice(1, None)
                target[axis] = slice(0, -1)
            neighbor[tuple(target)] = mask[tuple(source)]
            faces = np.argwhere(mask & ~neighbor)
            if not len(faces):
                continue
            plane = (faces[:, axis] + direction * 0.5) * spacing[axis]
            squared = (point_um[axis] - plane) ** 2
            for other in range(3):
                if other == axis:
                    continue
                lower = (faces[:, other] - 0.5) * spacing[other]
                upper = (faces[:, other] + 0.5) * spacing[other]
                delta = np.maximum(lower - point_um[other], 0) + np.minimum(upper - point_um[other], 0)
                squared += delta ** 2
            best_squared = min(best_squared, float(squared.min()))
    return float(np.sqrt(best_squared))


class _PairingCandidate(TypedDict):
    """The fixed schema of one cell/nucleus pairing candidate record, used
    throughout ``pair_and_filter_cells`` (built, filtered, ranked, and read
    back for the output rows). A plain ``dict[str, int | float | str]`` gives
    every key that same widened union on read, which breaks the `<`/`>`/`*`
    operations this function does on `containment`/`nc_ratio`/`overlap`
    elsewhere; a TypedDict keeps each key's own real type.
    """

    cell: int
    nucleus: int
    overlap: int
    containment: float
    nc_ratio: float
    failure: str


def pair_and_filter_cells(
    cell_masks: np.ndarray,
    nuclei_masks: np.ndarray,
    spacing: tuple[float, float, float],
    *,
    min_cell_volume_um3: float = 200.0,
    min_nucleus_volume_um3: float = 25.0,
    require_nucleus: bool = True,
    max_nc_ratio: float = 1.0,
    min_nucleus_containment: float = 0.5,
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame, dict]:
    """Pair cell and nucleus instances, apply QC, and compact surviving IDs.

    Sparse weighted matching first maximizes the number of one-to-one pairs,
    then maximizes total overlap among those solutions. Unlike the source
    implementation, nucleus volume and N:C ratio use the complete nucleus,
    not only the overlapping part.

    The underlying algorithm is the same assignment-problem theory as the
    Hungarian matching in ``evaluation.py``: Kuhn, H. W. (1955). Naval
    Research Logistics Quarterly, 2(1-2), 83-97.
    https://doi.org/10.1002/nav.3800020109. See docs/ALGORITHM_DECISIONS.md
    D10 for this specific application (bipartite cell-nucleus pairing) and
    docs/PARAMETERS.md for the QC-gate threshold provenance.
    """
    spacing = _validate_inputs(cell_masks, nuclei_masks, spacing)
    _validate_qc(min_cell_volume_um3, min_nucleus_volume_um3, max_nc_ratio, min_nucleus_containment)
    voxel_volume = float(np.prod(spacing))

    cell_ids, cell_counts = np.unique(cell_masks[cell_masks > 0], return_counts=True)
    nucleus_ids, nucleus_counts = np.unique(nuclei_masks[nuclei_masks > 0], return_counts=True)
    cell_sizes = {int(label): int(count) for label, count in zip(cell_ids, cell_counts)}
    nucleus_sizes = {int(label): int(count) for label, count in zip(nucleus_ids, nucleus_counts)}

    overlap_voxels: dict[tuple[int, int], int] = {}
    overlapping = (cell_masks > 0) & (nuclei_masks > 0)
    if overlapping.any():
        pairs, counts = np.unique(
            # A structured pair preserves each integer dtype. Stacking int64
            # and uint64 promotes to float64 and merges IDs above 2**53.
            np.rec.fromarrays((cell_masks[overlapping], nuclei_masks[overlapping])),
            return_counts=True,
        )
        overlap_voxels = {
            (int(cell_id), int(nucleus_id)): int(count)
            for (cell_id, nucleus_id), count in zip(pairs, counts)
        }

    candidates: dict[int, list[_PairingCandidate]] = {int(cell_id): [] for cell_id in cell_ids}
    for (cell_id, nucleus_id), overlap in overlap_voxels.items():
        cell_count = cell_sizes[cell_id]
        nucleus_count = nucleus_sizes[nucleus_id]
        containment = overlap / nucleus_count
        nc_ratio = nucleus_count / cell_count
        if nucleus_count * voxel_volume < min_nucleus_volume_um3:
            failure = "too_small_nucleus"
        elif containment < min_nucleus_containment:
            failure = "nucleus_not_contained"
        elif nc_ratio > max_nc_ratio:
            failure = "nc_ratio_too_high"
        else:
            failure = ""
        record: _PairingCandidate = {
            "cell": cell_id,
            "nucleus": nucleus_id,
            "overlap": overlap,
            "containment": containment,
            "nc_ratio": nc_ratio,
            "failure": failure,
        }
        candidates[cell_id].append(record)

    eligible_cells = {
        cell_id for cell_id, count in cell_sizes.items()
        if count * voxel_volume >= min_cell_volume_um3
    }
    adjacency = {
        cell_id: sorted(
            (candidate for candidate in candidates[cell_id] if not candidate["failure"]),
            key=lambda item: (-item["overlap"], item["nucleus"]),
        )
        for cell_id in eligible_cells
    }
    assignments: dict[int, _PairingCandidate] = {}
    matching_cells = sorted(eligible_cells)
    matching_nuclei = sorted({item["nucleus"] for group in adjacency.values() for item in group})
    if matching_cells:
        cell_index = {cell_id: index for index, cell_id in enumerate(matching_cells)}
        nucleus_index = {nucleus_id: index for index, nucleus_id in enumerate(matching_nuclei)}
        candidate_lookup = {
            (item["cell"], item["nucleus"]): item
            for group in adjacency.values() for item in group
        }
        overlap_total = sum(item["overlap"] for item in candidate_lookup.values())
        cardinality_bonus = overlap_total + 1
        row_indices, column_indices, weights = [], [], []
        # Named `entry`, not `candidate`, to avoid colliding with the
        # `candidate: _PairingCandidate | None` variable used later in this
        # function -- same dict shape, but a different, non-Optional role.
        for (cell_id, nucleus_id), entry in candidate_lookup.items():
            row_indices.append(cell_index[cell_id])
            column_indices.append(nucleus_index[nucleus_id])
            weights.append(cardinality_bonus + entry["overlap"])
        # Every cell has a unique dummy column, guaranteeing a full matching.
        for row_index in range(len(matching_cells)):
            row_indices.append(row_index)
            column_indices.append(len(matching_nuclei) + row_index)
            weights.append(1)
        graph = csr_matrix(
            (np.asarray(weights, float), (row_indices, column_indices)),
            shape=(len(matching_cells), len(matching_nuclei) + len(matching_cells)),
        )
        matched_rows, matched_columns = min_weight_full_bipartite_matching(graph, maximize=True)
        for row_index, column_index in zip(matched_rows, matched_columns):
            if column_index < len(matching_nuclei):
                cell_id = matching_cells[row_index]
                nucleus_id = matching_nuclei[column_index]
                assignments[cell_id] = candidate_lookup[(cell_id, nucleus_id)]

    decisions: dict[int, tuple[bool, str, _PairingCandidate | None]] = {}
    removed: dict[str, int] = {}
    for cell_id in sorted(cell_sizes):
        # Declared once so every branch below is checked against the same
        # real type -- `candidate` genuinely is a dict on some paths (a
        # matched/dominant candidate record) and None on others (no
        # candidate at all), matching `decisions`'s own declared value type.
        kept: bool
        reason: str
        candidate: _PairingCandidate | None
        if cell_sizes[cell_id] * voxel_volume < min_cell_volume_um3:
            kept, reason, candidate = False, "too_small_cell", None
        elif cell_id in assignments:
            kept, reason, candidate = True, "", assignments[cell_id]
        elif not candidates[cell_id]:
            kept, reason, candidate = not require_nucleus, "" if not require_nucleus else "no_nucleus", None
        else:
            dominant = max(candidates[cell_id], key=lambda item: (item["overlap"], -item["nucleus"]))
            if dominant["failure"]:
                kept, reason, candidate = False, dominant["failure"], dominant
            else:
                kept = not require_nucleus
                reason = "" if kept else "nucleus_already_paired"
                candidate = dominant
        decisions[cell_id] = kept, reason, candidate
        if not kept:
            removed[reason] = removed.get(reason, 0) + 1

    rows = []
    next_id = 0
    cell_new_ids: dict[int, int] = {}
    nucleus_new_ids: dict[int, int] = {}
    for cell_id in sorted(decisions):
        kept, reason, candidate = decisions[cell_id]
        new_id = 0
        nucleus_id = 0
        if kept:
            next_id += 1
            new_id = next_id
            cell_new_ids[cell_id] = new_id
            if cell_id in assignments:
                nucleus_id = assignments[cell_id]["nucleus"]
                nucleus_new_ids[nucleus_id] = new_id
        cell_volume = cell_sizes[cell_id] * voxel_volume
        audit_nucleus_id = candidate["nucleus"] if candidate is not None else 0
        rows.append({
            "cell_id": new_id,
            "original_cell_id": cell_id,
            "nucleus_id": new_id if nucleus_id else 0,
            "original_nucleus_id": nucleus_id,
            "candidate_nucleus_id": audit_nucleus_id,
            "kept": kept,
            "reject_reason": reason,
            "cell_volume_um3": cell_volume,
            "nucleus_volume_um3": nucleus_sizes[audit_nucleus_id] * voxel_volume if audit_nucleus_id else np.nan,
            "overlap_volume_um3": candidate["overlap"] * voxel_volume if candidate is not None else np.nan,
            "nucleus_containment_fraction": candidate["containment"] if candidate is not None else np.nan,
            "nc_ratio": candidate["nc_ratio"] if candidate is not None else np.nan,
        })

    # Both label volumes use the same sparse-safe remapping rule: only paired,
    # surviving source IDs receive their compact output ID.
    filtered_cells = relabel_from_source_ids(cell_masks, cell_new_ids)
    filtered_nuclei = relabel_from_source_ids(nuclei_masks, nucleus_new_ids)

    pairing = pd.DataFrame(rows, columns=PAIR_COLUMNS)
    paired = int((pairing["nucleus_id"] > 0).sum()) if len(pairing) else 0
    # A distinct name from the per-row `kept: bool` above -- this is a count,
    # not a flag, and reusing the name made mypy see a bool/int conflict.
    kept_count = int(pairing["kept"].sum()) if len(pairing) else 0
    summary = {
        "cells_before": len(cell_sizes),
        "nuclei_before": len(nucleus_sizes),
        "cells_after": kept_count,
        "nuclei_after": paired,
        "paired_cells": paired,
        "unpaired_cells": kept_count - paired,
        "removed": removed,
    }
    return filtered_cells, filtered_nuclei, pairing, summary


def cell_geometry(cell_labels: np.ndarray, nucleus_labels: np.ndarray | None,
                  spacing: tuple[float, float, float], *, fill_holes: bool = True) -> pd.DataFrame:
    """Measure cells; preserve the legacy envelope default, allow explicit raw support.

    cell_volume_um3 always counts raw cell voxels. The other shape measurements
    use measurement_basis; cell_envelope_volume_um3 retains its envelope meaning.
    """
    if not isinstance(fill_holes, bool):
        raise ValueError("fill_holes must be a boolean")
    if nucleus_labels is None:
        nucleus_labels = np.zeros_like(cell_labels)
    spacing = _validate_inputs(cell_labels, nucleus_labels, spacing)
    voxel_volume = float(np.prod(spacing))
    compact_cells, cell_mapping = compact_instance_labels(cell_labels)
    compact_nuclei, nucleus_mapping = compact_instance_labels(nucleus_labels)
    present_nuclei, nucleus_counts = np.unique(nucleus_labels[nucleus_labels > 0], return_counts=True)
    nucleus_sizes = {int(label): int(count) for label, count in zip(present_nuclei, nucleus_counts)}
    compact_nucleus_ids = np.array(sorted(nucleus_mapping), dtype=int)
    nucleus_centers = dict(zip(
        (nucleus_mapping[label] for label in compact_nucleus_ids),
        ndi.center_of_mass(
            np.ones(compact_nuclei.shape, dtype=np.uint8),
            compact_nuclei,
            compact_nucleus_ids,
        ) if len(compact_nucleus_ids) else [],
    ))

    rows = []
    for compact_id, bbox in enumerate(ndi.find_objects(compact_cells), start=1):
        if bbox is None:
            continue
        cell_id = cell_mapping[compact_id]
        cell_mask = compact_cells[bbox] == compact_id
        origin = tuple(item.start for item in bbox)
        measured, _ = geometry(cell_mask, spacing, origin, fill_holes=fill_holes)
        cell_volume = measured["segmented_volume_um3"]
        nucleus_volume = float(nucleus_sizes.get(cell_id, 0) * voxel_volume)
        row = {
            "cell_id": cell_id,
            "cell_volume_um3": cell_volume,
            "cell_envelope_volume_um3": (measured["volume_um3"] if fill_holes
                                         else float(outer_envelope(cell_mask).sum() * voxel_volume)),
            "surface_area_um2": measured["surface_area_um2"],
            "sphericity": measured["sphericity"],
            "equivalent_diameter_um": measured["equivalent_diameter_um"],
            "centroid_z_um": measured["centroid_z_um"],
            "centroid_y_um": measured["centroid_y_um"],
            "centroid_x_um": measured["centroid_x_um"],
            "principal_axis_major_um": measured["principal_axis_major_um"],
            "principal_axis_intermediate_um": measured["principal_axis_intermediate_um"],
            "principal_axis_minor_um": measured["principal_axis_minor_um"],
            "axis_ratio_minor_to_major": measured["axis_ratio_minor_to_major"],
            "n_z_slices": measured["n_z_slices"],
            "enclosed_void_fraction": measured["enclosed_void_fraction"],
            "touches_border": bbox_touches_volume_boundary(bbox, cell_labels.shape),
            "nucleus_volume_um3": nucleus_volume,
            "nc_ratio": nucleus_volume / cell_volume if nucleus_volume and cell_volume else np.nan,
            **{name: measured[name] for name in ("volume_um3", "measurement_basis", *DERIVED_GEOMETRY_COLUMNS, *SURFACE_METADATA_COLUMNS)},
        }
        if cell_id in nucleus_centers:
            center_voxels = np.asarray(nucleus_centers[cell_id])
            center_um = center_voxels * spacing
            for axis, value in zip("zyx", center_um):
                row[f"nucleus_centroid_{axis}_um"] = float(value)
            cell_center = np.array([row[f"centroid_{axis}_um"] for axis in "zyx"])
            row["cell_to_nucleus_centroid_um"] = float(np.linalg.norm(center_um - cell_center))
            local_center = center_voxels - np.asarray(origin)
            row["nucleus_centroid_to_cell_border_um"] = _distance_to_mask_surface(
                cell_mask, local_center, spacing
            )
        rows.append(row)
    frame = pd.DataFrame(rows, columns=GEOMETRY_COLUMNS)
    frame.attrs["measurement_policy"] = measurement_policy("cell", "filled_envelope" if fill_holes else "raw_label")
    return frame


def cell_neighborhood(cell_labels: np.ndarray, spacing: tuple[float, float, float],
                      radius_um: float = 25.0) -> pd.DataFrame:
    """Measure physical-radius cell neighborhoods with a cKDTree."""
    spacing = _validate_single("cell_labels", cell_labels, spacing)
    if not np.isfinite(radius_um) or radius_um <= 0:
        raise ValueError("radius_um must be positive and finite")
    compact, mapping = compact_instance_labels(cell_labels)
    ids, centers = [], []
    for compact_id, bbox in enumerate(ndi.find_objects(compact), start=1):
        if bbox is None:
            continue
        cell_id = mapping[compact_id]
        local = np.argwhere(compact[bbox] == compact_id)
        centers.append((local.mean(axis=0) + np.array([item.start for item in bbox])) * spacing)
        ids.append(cell_id)
    columns = ["cell_id", "n_neighbors", "nearest_neighbor_distance_um",
               "mean_neighbor_distance_um", "max_neighbor_distance_um",
               "local_density_per_mm3", "neighborhood_complete"]
    if not ids:
        return pd.DataFrame(columns=columns)

    points = np.asarray(centers)
    tree = cKDTree(points)
    neighborhoods = tree.query_ball_tree(tree, radius_um)
    sphere_mm3 = (4 / 3 * np.pi * radius_um ** 3) / 1e9
    image_extent = np.asarray(cell_labels.shape) * spacing
    rows = []
    for index, cell_id in enumerate(ids):
        neighbors = [other for other in neighborhoods[index] if other != index]
        distances = np.linalg.norm(points[neighbors] - points[index], axis=1) if neighbors else np.array([])
        center_to_face = np.minimum(points[index] + np.asarray(spacing) / 2,
                                    image_extent - points[index] - np.asarray(spacing) / 2)
        rows.append({
            "cell_id": cell_id,
            "n_neighbors": len(neighbors),
            "nearest_neighbor_distance_um": float(distances.min()) if len(distances) else np.nan,
            "mean_neighbor_distance_um": float(distances.mean()) if len(distances) else np.nan,
            "max_neighbor_distance_um": float(distances.max()) if len(distances) else np.nan,
            "local_density_per_mm3": len(neighbors) / sphere_mm3 if center_to_face.min() >= radius_um else np.nan,
            "neighborhood_complete": bool(center_to_face.min() >= radius_um),
        })
    return pd.DataFrame(rows, columns=columns)


class _PairingQCConfig(TypedDict):
    """The exact keyword-only QC schema ``pair_and_filter_cells`` accepts.

    A plain dict here mixes float thresholds with one bool
    (``require_nucleus``), so mypy would widen every value to a single type
    and could not check the ``**qc`` call below per-key; a TypedDict keeps
    each key's real, distinct type and lets mypy verify the unpacking.
    """

    min_cell_volume_um3: float
    min_nucleus_volume_um3: float
    require_nucleus: bool
    max_nc_ratio: float
    min_nucleus_containment: float


def analyze_cells(cell_masks: np.ndarray, nuclei_masks: np.ndarray,
                  spacing: tuple[float, float, float], *, neighbor_radius_um: float = 25.0,
                  min_cell_volume_um3: float = 200.0,
                  min_nucleus_volume_um3: float = 25.0,
                  require_nucleus: bool = True, max_nc_ratio: float = 1.0,
                  min_nucleus_containment: float = 0.5) -> CellAnalysisResult:
    """Run N:C pairing QC, geometry, and neighborhood analysis."""
    qc: _PairingQCConfig = {
        "min_cell_volume_um3": min_cell_volume_um3,
        "min_nucleus_volume_um3": min_nucleus_volume_um3,
        "require_nucleus": require_nucleus,
        "max_nc_ratio": max_nc_ratio,
        "min_nucleus_containment": min_nucleus_containment,
    }
    cells, nuclei, pairing, summary = pair_and_filter_cells(cell_masks, nuclei_masks, spacing, **qc)
    measured = cell_geometry(cells, nuclei, spacing)
    kept_pairs = pairing.loc[pairing["kept"].eq(True), [
        "cell_id", "original_cell_id", "original_nucleus_id", "nucleus_containment_fraction",
    ]]
    measured = measured.merge(kept_pairs, on="cell_id", how="left")
    measured = measured.merge(cell_neighborhood(cells, spacing, neighbor_radius_um), on="cell_id", how="left")
    summary["neighbor_radius_um"] = float(neighbor_radius_um)
    summary["qc"] = qc
    summary["measurement_policy"] = measurement_policy("cell", "filled_envelope")
    return CellAnalysisResult(cells, nuclei, measured, pairing, summary)
