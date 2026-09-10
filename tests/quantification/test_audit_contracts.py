"""Deterministic input-domain and label-identity audit regressions."""
from dataclasses import replace

import numpy as np
import pytest

from organoid_analysis.quantification.cellular_measurements import pair_and_filter_cells
from organoid_analysis.quantification.multilevel_relationships.config import Multilevel3DConfig


@pytest.mark.parametrize('field', ['minimum_voxels', 'mad_z_threshold'])
@pytest.mark.parametrize('value', [np.nan, np.inf, -np.inf, True])
def test_multilevel_qc_rejects_invalid_numeric_values(field, value):
    with pytest.raises(ValueError):
        replace(Multilevel3DConfig(), **{field: value}).validate()


def test_multilevel_voxel_count_requires_integer():
    with pytest.raises(ValueError):
        Multilevel3DConfig(minimum_voxels=1.5).validate()


def test_pairing_preserves_mixed_signed_unsigned_large_ids():
    cells = np.ones((3, 4, 4), dtype=np.int64)
    cells[1:] = 2
    nuclei = np.full(cells.shape, 2**53, dtype=np.uint64)
    nuclei[1:] = 2**53 + 1
    _, _, pairs, summary = pair_and_filter_cells(
        cells, nuclei, (1, 1, 1), min_cell_volume_um3=0, min_nucleus_volume_um3=0,
    )
    assert summary['paired_cells'] == 2
    assert pairs.original_nucleus_id.tolist() == [2**53, 2**53 + 1]
    assert pairs.nucleus_containment_fraction.tolist() == [1., 1.]


def test_mask_feature_sparse_ids_are_compacted_before_regionprops(monkeypatch):
    from skimage import measure

    from organoid_analysis.quantification.mask_features import extract_mask_features
    original = measure.regionprops

    def bounded_regionprops(labels, **kwargs):
        # Prevent the pre-fix implementation allocating billions of boxes.
        assert int(labels.max()) <= 2
        return original(labels, **kwargs)

    monkeypatch.setattr(measure, 'regionprops', bounded_regionprops)
    mask = np.zeros((8, 8, 8), np.uint32)
    mask[1:3, 1:3, 1:3] = 7
    mask[4:6, 4:6, 4:6] = np.iinfo(np.uint32).max
    result = extract_mask_features(mask, spacing_um=(2, 3, 4))
    assert result.index.tolist() == [7, np.iinfo(np.uint32).max]
    assert result.volume_um3.tolist() == [192., 192.]


def test_boolean_mask_contract_and_empty_summary():
    from organoid_analysis.quantification.mask_features import extract_mask_features, summarize_mask
    mask = np.zeros((5, 5, 5), bool)
    mask[1:4, 1:4, 1:4] = True
    assert extract_mask_features(mask).volume_um3.tolist() == [27.]
    empty = summarize_mask(np.zeros_like(mask))
    assert empty.n_objects == 0
    assert empty.total_volume_um3 == 0
    assert np.isnan(empty.mean_volume_um3)
