"""Cell positions relative to parent organoids in physical 3D space."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import ndimage as ndi


def cell_spatial_features(cell_features: pd.DataFrame, organoid_features: pd.DataFrame,
                          organoid_labels: np.ndarray, spacing_zyx_um: tuple[float, float, float]) -> pd.DataFrame:
    """Measure centroid radius and EDT depth for cells with an assigned organoid.

    Depth is the physical Euclidean distance transform value at the nearest
    voxel to the cell centroid.  Padding supplies an exterior at all image
    borders, including for border-touching organoids.
    """
    organoids = organoid_features.set_index("organoid_id")
    spacing = np.asarray(spacing_zyx_um, dtype=float)
    # Build one physical depth map per parent organoid, rather than repeating
    # an O(volume) EDT for every child cell.
    parent_cache: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray, float]] = {}
    for organoid_id, organoid in organoids.iterrows():
        mask = organoid_labels == int(organoid_id)
        if not mask.any():
            continue
        centroid = organoid[["centroid_z_um", "centroid_y_um", "centroid_x_um"]].to_numpy(float)
        radius = float((3.0 * float(organoid.volume_um3) / (4.0 * np.pi)) ** (1.0 / 3.0))
        depth_map = ndi.distance_transform_edt(np.pad(mask, 1), sampling=spacing)[1:-1, 1:-1, 1:-1]
        parent_cache[int(organoid_id)] = (mask, depth_map, centroid, radius)

    rows = []
    for cell in cell_features.itertuples(index=False):
        organoid_id = int(cell.organoid_id)
        base = {"cell_id": int(cell.cell_id), "distance_to_organoid_centroid_um": np.nan,
                "equivalent_organoid_radius_um": np.nan, "normalized_radial_position_equivalent_radius": np.nan,
                "distance_to_organoid_surface_um": np.nan, "normalized_depth_equivalent_radius": np.nan}
        if organoid_id <= 0 or organoid_id not in organoids.index:
            rows.append(base)
            continue
        mask, depth_map, parent_centroid, radius = parent_cache.get(organoid_id, (None, None, None, np.nan))
        if mask is None:
            rows.append(base)
            continue
        centroid = np.array([cell.centroid_z_um, cell.centroid_y_um, cell.centroid_x_um], dtype=float)
        radial = float(np.linalg.norm(centroid - parent_centroid))
        point = np.clip(np.rint(centroid / spacing).astype(int), 0, np.asarray(mask.shape) - 1)
        depth = float(depth_map[tuple(point)]) if mask[tuple(point)] else 0.0
        base.update({"distance_to_organoid_centroid_um": radial, "equivalent_organoid_radius_um": radius,
                     "normalized_radial_position_equivalent_radius": radial / radius if radius else np.nan,
                     "distance_to_organoid_surface_um": depth,
                     "normalized_depth_equivalent_radius": depth / radius if radius else np.nan})
        rows.append(base)
    return pd.DataFrame(rows, columns=["cell_id", "distance_to_organoid_centroid_um", "equivalent_organoid_radius_um",
                                       "normalized_radial_position_equivalent_radius", "distance_to_organoid_surface_um",
                                       "normalized_depth_equivalent_radius"])
