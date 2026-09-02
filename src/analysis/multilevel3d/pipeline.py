"""Orchestration of hierarchy, 3D morphology, topology, spatial features, QC."""
from __future__ import annotations

from dataclasses import dataclass
import time

import numpy as np
import pandas as pd

from .config import Multilevel3DConfig
from .hierarchy import assign_parents_by_overlap
from .morphology import measure_instances
from .qc import add_qc_flags
from .spatial import cell_spatial_features
from .topology import direct_contact_edges, topology_per_cell
from .validation import validate_inputs


@dataclass
class Multilevel3DResult:
    """In-memory feature database produced from one registered label volume."""

    organoid_features: pd.DataFrame
    cell_features: pd.DataFrame
    nucleus_features: pd.DataFrame
    cell_topology_edges: pd.DataFrame
    qc_flags: pd.DataFrame
    summary: dict


def _add_metadata(frame: pd.DataFrame, metadata: dict[str, str] | None) -> pd.DataFrame:
    """Prepend only metadata explicitly supplied by a caller; never invent it."""
    if not metadata:
        return frame
    values = {key: value for key, value in metadata.items() if value not in (None, "")}
    if not values:
        return frame
    result = frame.copy()
    for key, value in reversed(list(values.items())):
        result.insert(0, key, value)
    return result


def _add_nuclear_intensity(nuclei: pd.DataFrame, nucleus_labels: np.ndarray,
                           intensity: np.ndarray | None) -> pd.DataFrame:
    """Measure raw nuclear-channel intensity once per labelled nucleus."""
    if intensity is None:
        return nuclei
    if nuclei.empty:
        return nuclei.copy()
    ids = np.sort(nuclei.nucleus_id.to_numpy(dtype=int))
    foreground = nucleus_labels > 0
    label_indices = np.searchsorted(ids, nucleus_labels[foreground])
    values = intensity[foreground].astype(float, copy=False)
    counts = np.bincount(label_indices, minlength=len(ids)).astype(float)
    sums = np.bincount(label_indices, weights=values, minlength=len(ids))
    sum_squares = np.bincount(label_indices, weights=values * values, minlength=len(ids))
    minima = np.full(len(ids), np.inf)
    maxima = np.full(len(ids), -np.inf)
    np.minimum.at(minima, label_indices, values)
    np.maximum.at(maxima, label_indices, values)
    means = sums / counts
    variances = np.maximum(0.0, sum_squares / counts - means * means)
    stds = np.sqrt(variances)
    measurements = pd.DataFrame({
        "nucleus_id": ids,
        "mean_nuclear_intensity": means,
        "max_nuclear_intensity": maxima,
        "min_nuclear_intensity": minima,
        "std_nuclear_intensity": stds,
        "integrated_nuclear_intensity": sums,
        "CV_chromatin": np.divide(stds, means, out=np.full_like(stds, np.nan), where=means != 0),
    })
    return nuclei.merge(measurements, on="nucleus_id", how="left", validate="one_to_one")


def _cell_nucleus_features(cells: pd.DataFrame, nuclei: pd.DataFrame) -> pd.DataFrame:
    grouped = nuclei.loc[nuclei.cell_id > 0].groupby("cell_id", sort=False)
    counts = grouped.size().rename("nucleus_count")
    volumes = grouped.volume_um3.agg([("nucleus_volume_total_um3", "sum"), ("mean_nucleus_volume_um3", "mean"),
                                      ("_largest_nucleus_volume_um3", "max")])
    result = cells.merge(pd.concat([counts, volumes], axis=1), on="cell_id", how="left")
    result["nucleus_count"] = result["nucleus_count"].fillna(0).astype(int)
    for column in ("nucleus_volume_total_um3", "mean_nucleus_volume_um3", "_largest_nucleus_volume_um3"):
        result[column] = result[column].fillna(0.0)
    result["nucleus_to_cell_volume_ratio"] = result["nucleus_volume_total_um3"] / result["volume_um3"]
    result["anucleate"] = result.nucleus_count.eq(0)
    result["multinucleated"] = result.nucleus_count.gt(1)
    result["largest_nucleus_fraction"] = np.where(result.nucleus_volume_total_um3 > 0,
                                                    result._largest_nucleus_volume_um3 / result.nucleus_volume_total_um3,
                                                    np.nan)
    result = result.drop(columns="_largest_nucleus_volume_um3")
    intensity_columns = {
        "mean_nuclear_intensity", "max_nuclear_intensity", "min_nuclear_intensity",
        "std_nuclear_intensity", "integrated_nuclear_intensity", "CV_chromatin",
    }
    if not intensity_columns.issubset(nuclei.columns):
        return result

    nuclear = nuclei.loc[nuclei.cell_id > 0].copy()
    voxel_weights = nuclear["voxel_count"].astype(float)
    nuclear["_intensity_sum_squares"] = voxel_weights * (
        nuclear["std_nuclear_intensity"] ** 2 + nuclear["mean_nuclear_intensity"] ** 2
    )
    intensity = nuclear.groupby("cell_id", sort=False).agg(
        _nuclear_intensity_voxels=("voxel_count", "sum"),
        integrated_nuclear_intensity=("integrated_nuclear_intensity", "sum"),
        min_nuclear_intensity=("min_nuclear_intensity", "min"),
        max_nuclear_intensity=("max_nuclear_intensity", "max"),
        _intensity_sum_squares=("_intensity_sum_squares", "sum"),
    )
    intensity["mean_nuclear_intensity"] = (
        intensity["integrated_nuclear_intensity"] / intensity["_nuclear_intensity_voxels"]
    )
    intensity["std_nuclear_intensity"] = np.sqrt(np.maximum(
        0.0,
        intensity["_intensity_sum_squares"] / intensity["_nuclear_intensity_voxels"]
        - intensity["mean_nuclear_intensity"] ** 2,
    ))
    intensity["CV_chromatin"] = np.divide(
        intensity["std_nuclear_intensity"], intensity["mean_nuclear_intensity"],
        out=np.full(len(intensity), np.nan), where=intensity["mean_nuclear_intensity"].to_numpy() != 0,
    )
    intensity = intensity.drop(columns=["_nuclear_intensity_voxels", "_intensity_sum_squares"])
    return result.merge(intensity, on="cell_id", how="left", validate="one_to_one")


def _organoid_aggregates(organoids: pd.DataFrame, cells: pd.DataFrame, nuclei: pd.DataFrame,
                         config: Multilevel3DConfig) -> pd.DataFrame:
    rows = []
    for organoid in organoids.itertuples(index=False):
        organoid_id = int(organoid.organoid_id)
        child_cells = cells.loc[cells.organoid_id.eq(organoid_id)]
        child_nuclei = nuclei.loc[nuclei.organoid_id.eq(organoid_id)]
        cell_volumes = child_cells.volume_um3.to_numpy(float)
        nucleus_volumes = child_nuclei.volume_um3.to_numpy(float)
        rows.append({"organoid_id": organoid_id,
                     "cell_count": int(len(child_cells)), "nucleus_count": int(len(child_nuclei)),
                     "cell_density_per_um3": float(len(child_cells) / organoid.volume_um3),
                     "nuclear_density_per_um3": float(len(child_nuclei) / organoid.volume_um3),
                     "mean_cell_volume_um3": float(cell_volumes.mean()) if len(cell_volumes) else np.nan,
                     "median_cell_volume_um3": float(np.median(cell_volumes)) if len(cell_volumes) else np.nan,
                     "cell_volume_CV": float(cell_volumes.std() / cell_volumes.mean()) if len(cell_volumes) and cell_volumes.mean() else np.nan,
                     "mean_nucleus_volume_um3": float(nucleus_volumes.mean()) if len(nucleus_volumes) else np.nan,
                     "nucleus_volume_CV": float(nucleus_volumes.std() / nucleus_volumes.mean()) if len(nucleus_volumes) and nucleus_volumes.mean() else np.nan,
                     "anucleate_cell_fraction": float(child_cells.anucleate.mean()) if len(child_cells) else np.nan,
                     "multinucleated_cell_fraction": float(child_cells.multinucleated.mean()) if len(child_cells) else np.nan,
                     "mean_neighbor_count": float(child_cells.neighbor_count.mean()) if len(child_cells) else np.nan,
                     "mean_contact_area_um2": float(child_cells.mean_contact_area_um2.mean()) if len(child_cells) else np.nan,
                     "core_cell_fraction": float((child_cells.normalized_radial_position_equivalent_radius <= config.core_max_normalized_radial_position).mean()) if len(child_cells) else np.nan,
                     "peripheral_cell_fraction": float((child_cells.normalized_radial_position_equivalent_radius >= config.peripheral_min_normalized_radial_position).mean()) if len(child_cells) else np.nan})
    aggregation_columns = ["organoid_id", "cell_count", "nucleus_count", "cell_density_per_um3", "nuclear_density_per_um3",
                           "mean_cell_volume_um3", "median_cell_volume_um3", "cell_volume_CV", "mean_nucleus_volume_um3",
                           "nucleus_volume_CV", "anucleate_cell_fraction", "multinucleated_cell_fraction", "mean_neighbor_count",
                           "mean_contact_area_um2", "core_cell_fraction", "peripheral_cell_fraction"]
    return organoids.merge(pd.DataFrame(rows, columns=aggregation_columns), on="organoid_id", how="left")


def analyze_multilevel_3d(
    organoid_labels: np.ndarray,
    cell_labels: np.ndarray,
    nucleus_labels: np.ndarray,
    spacing_zyx_um: tuple[float, float, float],
    *,
    config: Multilevel3DConfig | None = None,
    metadata: dict[str, str] | None = None,
    nucleus_intensity: np.ndarray | None = None,
) -> Multilevel3DResult:
    """Analyze a registered organoid → cell → nucleus instance hierarchy.

    Parent assignment is maximum voxel overlap.  All scalar morphology,
    topology contact areas, and distances use the supplied physical spacing.
    """
    started = time.perf_counter()
    cfg = config or Multilevel3DConfig()
    cfg.validate()
    spacing = validate_inputs(organoid_labels, cell_labels, nucleus_labels, spacing_zyx_um, nucleus_intensity)
    organoids = measure_instances(organoid_labels, spacing, object_type="organoid")
    cells = measure_instances(cell_labels, spacing, object_type="cell")
    nuclei = measure_instances(nucleus_labels, spacing, object_type="nucleus")
    cell_parent = assign_parents_by_overlap(organoid_labels, cell_labels, parent_name="organoid", child_name="cell")
    nucleus_cell_parent = assign_parents_by_overlap(cell_labels, nucleus_labels, parent_name="cell", child_name="nucleus")
    nucleus_organoid_parent = assign_parents_by_overlap(organoid_labels, nucleus_labels, parent_name="organoid", child_name="nucleus")
    cells = cells.merge(cell_parent, on="cell_id", how="left")
    direct_nucleus_organoid = nucleus_organoid_parent.rename(columns={
        "organoid_id": "direct_organoid_id",
        "organoid_parent_overlap_fraction": "direct_organoid_parent_overlap_fraction",
        "organoid_parent_candidate_count": "direct_organoid_parent_candidate_count",
        "crosses_multiple_organoids": "crosses_multiple_direct_organoids",
    })
    nuclei = nuclei.merge(nucleus_cell_parent, on="nucleus_id", how="left").merge(
        direct_nucleus_organoid, on="nucleus_id", how="left", validate="one_to_one")
    # A nucleus inherits its organoid exclusively through its assigned cell.
    # Direct nucleus-to-organoid overlap is retained as an audit measurement,
    # but never allowed to break the Organoid -> Cell -> Nucleus hierarchy.
    cell_organoid_lookup = cells.set_index("cell_id")["organoid_id"]
    nuclei["organoid_id"] = nuclei["cell_id"].map(cell_organoid_lookup).fillna(0).astype(int)
    # Stable child-oriented aliases are the public hierarchy contract; retain
    # the explicit parent-named fields as useful audit detail.
    cells["cell_parent_overlap_fraction"] = cells["organoid_parent_overlap_fraction"]
    nuclei["nucleus_parent_overlap_fraction"] = nuclei["cell_parent_overlap_fraction"]
    cells["parent_assignment_failed"] = cells.organoid_id.eq(0)
    cells["low_parent_overlap"] = cells.organoid_parent_overlap_fraction.lt(cfg.low_parent_overlap_fraction) & ~cells.parent_assignment_failed
    cells["crosses_multiple_parents"] = cells.crosses_multiple_organoids
    nuclei["parent_assignment_failed"] = nuclei.cell_id.eq(0)
    nuclei["low_parent_overlap"] = nuclei.cell_parent_overlap_fraction.lt(cfg.low_parent_overlap_fraction) & ~nuclei.parent_assignment_failed
    nuclei["crosses_multiple_parents"] = nuclei.crosses_multiple_cells
    nuclei["direct_organoid_parent_mismatch"] = (
        nuclei["direct_organoid_id"].gt(0)
        & nuclei["organoid_id"].gt(0)
        & nuclei["direct_organoid_id"].ne(nuclei["organoid_id"])
    )
    nuclei = _add_nuclear_intensity(nuclei, nucleus_labels, nucleus_intensity)
    cells = _cell_nucleus_features(cells, nuclei)
    parent_lookup = dict(zip(cells.cell_id.astype(int), cells.organoid_id.astype(int)))
    edges = direct_contact_edges(cell_labels, spacing, parent_lookup)
    cells = cells.merge(topology_per_cell(cells.cell_id, edges), on="cell_id", how="left")
    cells = cells.merge(cell_spatial_features(cells, organoids, organoid_labels, spacing), on="cell_id", how="left")
    organoids = _organoid_aggregates(organoids, cells, nuclei, cfg)
    organoids, organoid_qc = add_qc_flags(organoids, object_type="organoid", config=cfg)
    cells, cell_qc = add_qc_flags(cells, object_type="cell", config=cfg)
    nuclei, nucleus_qc = add_qc_flags(nuclei, object_type="nucleus", config=cfg)
    qc_flags = pd.concat([organoid_qc, cell_qc, nucleus_qc], ignore_index=True)
    organoids, cells, nuclei, edges, qc_flags = (_add_metadata(frame, metadata)
                                                  for frame in (organoids, cells, nuclei, edges, qc_flags))
    elapsed = time.perf_counter() - started
    summary = {"status": "complete", "input_shape_zyx": list(organoid_labels.shape), "spacing_zyx_um": list(spacing),
               "organoid_count": int(len(organoids)), "cell_count": int(len(cells)), "nucleus_count": int(len(nuclei)),
               "topology_edge_count": int(len(edges)), "runtime_seconds": elapsed, "config": cfg.as_dict()}
    if metadata:
        summary["metadata"] = {key: value for key, value in metadata.items() if value not in (None, "")}
    return Multilevel3DResult(organoids, cells, nuclei, edges, qc_flags, summary)
