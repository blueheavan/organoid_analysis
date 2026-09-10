from __future__ import annotations

import numpy as np
import pytest
import tifffile

from organoid_analysis.quantification.multilevel_relationships.hierarchy import (
    assign_parents_by_overlap,
)
from organoid_analysis.web_interface.multilevel_results import read_registered_organoid_labels
from organoid_analysis.workflows.multilevel_measurement_workflow import analyze_multilevel_3d


def test_mixed_integer_dtypes_preserve_large_parent_child_ids():
    parents = np.full((3, 4, 4), 17, np.int64)
    children = np.zeros_like(parents, dtype=np.uint64)
    children[0] = 2**53
    children[1] = 2**53 + 1
    result = assign_parents_by_overlap(parents, children, parent_name='organoid', child_name='cell')
    assert result.cell_id.tolist() == [2**53, 2**53 + 1]
    assert result.organoid_id.tolist() == [17, 17]
    assert result.organoid_parent_overlap_fraction.tolist() == [1., 1.]


def test_nuclear_and_pooled_variance_are_invariant_to_large_additive_offset():
    shape = (6, 8, 10)
    organoids = np.ones(shape, np.uint32)
    cells = np.ones(shape, np.uint32)
    nuclei = np.zeros(shape, np.uint32)
    nuclei[1:3, 2:4, 3:5] = 1
    nuclei[3:5, 2:4, 3:5] = 2
    intensity = np.zeros(shape)
    intensity[nuclei == 1] = np.arange(8)
    intensity[nuclei == 2] = np.arange(8) + 4
    base = analyze_multilevel_3d(organoids, cells, nuclei, (1, 1, 1), nucleus_intensity=intensity)
    shifted = analyze_multilevel_3d(organoids, cells, nuclei, (1, 1, 1), nucleus_intensity=intensity + 1e9)
    for result in (base, shifted):
        np.testing.assert_allclose(result.nucleus_features.std_nuclear_intensity, np.sqrt(5.25), rtol=1e-12)
        np.testing.assert_allclose(result.cell_features.std_nuclear_intensity, np.sqrt(9.25), rtol=1e-12)


def test_complex_intensity_is_rejected_not_silently_cast():
    labels = np.ones((3, 4, 5), np.uint32)
    with pytest.raises(ValueError, match='real numeric'):
        analyze_multilevel_3d(labels, labels, labels, (1, 1, 1), nucleus_intensity=np.ones(labels.shape, complex))


def test_web_organoid_grid_mismatch_is_rejected(tmp_path):
    path = tmp_path / 'mask.ome.tif'
    labels = np.ones((3, 4, 5), np.uint32)
    tifffile.imwrite(path, labels, ome=True, photometric='minisblack',
                     metadata={'axes': 'ZYX', 'PhysicalSizeX': 2., 'PhysicalSizeY': 2., 'PhysicalSizeZ': 4.})
    with pytest.raises(ValueError, match='spacing conflicts'):
        read_registered_organoid_labels(path, labels.shape, (2, 1, 1))
    np.testing.assert_array_equal(read_registered_organoid_labels(path, labels.shape, (4, 2, 2)), labels)
