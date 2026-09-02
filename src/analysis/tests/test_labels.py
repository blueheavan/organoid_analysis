"""Tests for shared sparse-label helpers used by multiple analysis routes."""
from __future__ import annotations

import numpy as np

from analysis.labels import (
    bbox_touches_volume_boundary,
    compact_instance_labels,
    relabel_from_source_ids,
)


def test_compact_instance_labels_preserves_sparse_source_ids():
    labels = np.array([[[0, 7], [1_000_000_000, 7]]], dtype=np.uint32)

    compact, source_ids = compact_instance_labels(labels)

    assert compact.tolist() == [[[0, 1], [2, 1]]]
    assert source_ids == {1: 7, 2: 1_000_000_000}


def test_relabel_from_source_ids_drops_unmapped_labels_without_dense_lookup():
    labels = np.array([[[0, 7], [1_000_000_000, 7]]], dtype=np.uint32)

    relabelled = relabel_from_source_ids(labels, {7: 3})

    assert relabelled.dtype == np.uint32
    assert relabelled.tolist() == [[[0, 3], [0, 3]]]


def test_bbox_boundary_helper_checks_all_three_axes():
    shape = (5, 6, 7)

    assert bbox_touches_volume_boundary((slice(1, 4), slice(2, 5), slice(0, 3)), shape)
    assert not bbox_touches_volume_boundary((slice(1, 4), slice(2, 5), slice(2, 5)), shape)
