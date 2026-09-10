"""Maximum-overlap parent assignment for label-instance hierarchies."""
from __future__ import annotations

import numpy as np
import pandas as pd


def assign_parents_by_overlap(parent_labels: np.ndarray, child_labels: np.ndarray, *, parent_name: str,
                              child_name: str) -> pd.DataFrame:
    """Assign every child to its maximum-overlap parent instance.

    Ties are deterministic (smallest parent ID).  ``parent_overlap_fraction``
    is overlap voxels / all voxels belonging to that child, so background and
    cross-parent spillover are both visible to QC.
    """
    child_ids, child_counts = np.unique(child_labels[child_labels > 0], return_counts=True)
    totals = {int(identifier): int(count) for identifier, count in zip(child_ids, child_counts)}
    by_child: dict[int, list[tuple[int, int]]] = {identifier: [] for identifier in totals}
    overlap = (parent_labels > 0) & (child_labels > 0)
    if overlap.any():
        pairs, counts = np.unique(
            # Mixed int64/uint64 would promote to float64 and merge IDs >2**53.
            np.column_stack((child_labels[overlap].astype(np.uint64),
                             parent_labels[overlap].astype(np.uint64))), axis=0, return_counts=True
        )
        for (child_id, parent_id), count in zip(pairs, counts):
            by_child[int(child_id)].append((int(parent_id), int(count)))

    child_column = f"{child_name}_id"
    rows = []
    for child_id in sorted(totals):
        candidates = by_child[child_id]
        if candidates:
            parent_id, overlap_count = min(candidates, key=lambda item: (-item[1], item[0]))
        else:
            parent_id, overlap_count = 0, 0
        rows.append({
            child_column: child_id,
            f"{parent_name}_id": parent_id,
            f"{parent_name}_parent_overlap_fraction": overlap_count / totals[child_id],
            f"{parent_name}_parent_candidate_count": len(candidates),
            f"crosses_multiple_{parent_name}s": len(candidates) > 1,
        })
    return pd.DataFrame(rows, columns=[child_column, f"{parent_name}_id", f"{parent_name}_parent_overlap_fraction",
                                       f"{parent_name}_parent_candidate_count", f"crosses_multiple_{parent_name}s"])
