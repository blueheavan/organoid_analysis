"""Face-adjacent 3D cell-contact graph construction."""
from __future__ import annotations

import numpy as np
import pandas as pd

EDGE_COLUMNS = ["organoid_id", "cell_id_1", "cell_id_2", "contact_area_um2", "inter_organoid_contact"]


def direct_contact_edges(cell_labels: np.ndarray, spacing_zyx_um: tuple[float, float, float], cell_to_organoid: dict[int, int]) -> pd.DataFrame:
    """Return one edge per pair of labels with voxel-face adjacency.

    A shared Z-face has Y*X area, a shared Y-face Z*X area, and a shared
    X-face Z*Y area; this preserves anisotropic physical contact areas.
    """
    spacing = np.asarray(spacing_zyx_um, dtype=float)
    contact_by_pair: dict[tuple[int, int], float] = {}
    for axis in range(3):
        first = [slice(None)] * 3
        second = [slice(None)] * 3
        first[axis], second[axis] = slice(None, -1), slice(1, None)
        left, right = cell_labels[tuple(first)], cell_labels[tuple(second)]
        valid = (left > 0) & (right > 0) & (left != right)
        if not valid.any():
            continue
        pairs, counts = np.unique(
            np.column_stack((np.minimum(left[valid], right[valid]), np.maximum(left[valid], right[valid]))),
            axis=0, return_counts=True,
        )
        face_area = float(np.prod(np.delete(spacing, axis)))
        for (cell_1, cell_2), count in zip(pairs, counts):
            key = (int(cell_1), int(cell_2))
            contact_by_pair[key] = contact_by_pair.get(key, 0.0) + int(count) * face_area
    rows = []
    for (cell_1, cell_2), area in sorted(contact_by_pair.items()):
        parent_1, parent_2 = cell_to_organoid.get(cell_1, 0), cell_to_organoid.get(cell_2, 0)
        same_parent = parent_1 > 0 and parent_1 == parent_2
        rows.append({"organoid_id": parent_1 if same_parent else 0, "cell_id_1": cell_1, "cell_id_2": cell_2,
                     "contact_area_um2": area, "inter_organoid_contact": not same_parent})
    return pd.DataFrame(rows, columns=EDGE_COLUMNS)


def topology_per_cell(cell_ids: pd.Series, edges: pd.DataFrame) -> pd.DataFrame:
    """Aggregate direct-contact graph degree and physical contact area by cell."""
    values = [int(value) for value in cell_ids]
    rows = []
    for cell_id in values:
        if edges.empty:
            areas = np.array([], dtype=float)
        else:
            areas = edges.loc[(edges.cell_id_1 == cell_id) | (edges.cell_id_2 == cell_id), "contact_area_um2"].to_numpy(float)
        rows.append({"cell_id": cell_id, "neighbor_count": len(areas), "degree": len(areas),
                     "total_contact_area_um2": float(areas.sum()),
                     "mean_contact_area_um2": float(areas.mean()) if len(areas) else 0.0})
    return pd.DataFrame(rows, columns=["cell_id", "neighbor_count", "degree", "total_contact_area_um2", "mean_contact_area_um2"])
