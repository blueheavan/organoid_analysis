"""Small, memory-safe helpers for 3D instance-label volumes.

Instance labels often originate in external tools, where IDs may be sparse or
very large.  These helpers deliberately avoid arrays sized by ``labels.max()``
so importing such masks cannot cause an unexpected large allocation.
"""
from __future__ import annotations

import numpy as np


def compact_instance_labels(labels: np.ndarray) -> tuple[np.ndarray, dict[int, int]]:
    """Return consecutive positive labels and a ``compact_id -> source_id`` map.

    ``scipy.ndimage.find_objects`` indexes its result by label value.  Compact
    labels make that operation safe for sparse source IDs while retaining an
    explicit mapping for measurements and audit output.
    """
    source_ids, inverse = np.unique(labels, return_inverse=True)
    positive = source_ids != 0
    compact_values = np.zeros(len(source_ids), dtype=np.uint32)
    compact_values[positive] = np.arange(1, positive.sum() + 1, dtype=np.uint32)
    source_id_by_compact_id = {
        compact_id: int(source_id)
        for compact_id, source_id in enumerate(source_ids[positive], start=1)
    }
    return compact_values[inverse].reshape(labels.shape), source_id_by_compact_id


def relabel_from_source_ids(labels: np.ndarray, source_to_new_id: dict[int, int]) -> np.ndarray:
    """Apply a sparse source-ID mapping, assigning unmapped IDs to background.

    A direct lookup table would be fast for dense labels but unsafe when an
    input contains one very large ID.  Looking up only the IDs that occur in
    the volume keeps memory proportional to the number of distinct labels.
    """
    source_ids, inverse = np.unique(labels, return_inverse=True)
    new_ids = np.fromiter(
        (source_to_new_id.get(int(source_id), 0) for source_id in source_ids),
        dtype=np.uint32,
        count=len(source_ids),
    )
    return new_ids[inverse].reshape(labels.shape)


def bbox_touches_volume_boundary(
    bbox: tuple[slice, slice, slice], shape: tuple[int, int, int]
) -> bool:
    """Whether a label bounding box reaches any face of its source volume."""
    return any(part.start == 0 or part.stop == length for part, length in zip(bbox, shape))
