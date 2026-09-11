"""Hand-counted geometry and failure semantics for Web mask measurements."""
import io

import numpy as np
import pandas as pd
import pytest

from organoid_analysis.quantification.mask_features import FEATURE_COLUMNS, extract_mask_features
from organoid_analysis.quantification.multilevel_relationships.morphology import measure_instances


def asymmetric_hollow_mask():
    mask = np.zeros((9, 10, 11), np.uint32)
    mask[1:8, 1:9, 1:10] = 71
    mask[2:4, 2:5, 2:6] = 0
    return mask


def test_raw_geometry_counts_and_centroid_use_same_voxels():
    mask = asymmetric_hollow_mask()
    result = extract_mask_features(mask, (.2, .3, .7)).loc[71]
    count = 7 * 8 * 9 - 2 * 3 * 4
    assert result.voxel_count == count
    assert result.volume_um3 == pytest.approx(count * .2 * .3 * .7, rel=1e-12)
    expected = [(504 * 4 - 24 * 2.5) / 480 * .7,
                (504 * 4.5 - 24 * 3) / 480 * .3,
                (504 * 5 - 24 * 3.5) / 480 * .2]
    np.testing.assert_allclose(result[[f"centroid_{axis}_um" for axis in "zyx"]].astype(float), expected, rtol=1e-12)
    assert result.measurement_basis == "raw_label"
    assert result.filled_voxel_count == 0
    assert result.solidity == pytest.approx(480 / 504)
    assert result.qc_flags == ""


def test_explicit_envelope_uses_filled_centroid_and_solidity():
    result = extract_mask_features(asymmetric_hollow_mask(), (.2, .3, .7), fill_holes=True).loc[71]
    assert result.volume_um3 == pytest.approx(504 * .2 * .3 * .7)
    assert result.segmented_volume_um3 == pytest.approx(480 * .2 * .3 * .7)
    np.testing.assert_allclose(result[[f"centroid_{axis}_um" for axis in "zyx"]].astype(float), [2.8, 1.35, 1.])
    assert result.solidity == 1
    assert result.measurement_basis == "filled_envelope"
    assert result.filled_voxel_count == 24


def test_raw_web_and_multilevel_measure_same_geometry():
    mask = asymmetric_hollow_mask()
    web = extract_mask_features(mask, (.2, .3, .7)).loc[71]
    multilevel = measure_instances(mask, (.7, .3, .2), object_type="nucleus").iloc[0]
    # Integration consistency is separate from the hand-counted oracle above.
    for name in ("volume_um3", "surface_area_um2", "sphericity", "centroid_z_um", "centroid_y_um", "centroid_x_um"):
        assert web[name] == multilevel[name]


def test_small_volumes_keep_precision_and_missing_solidity_has_qc():
    mask = np.zeros((3, 3, 3), np.uint32)
    mask[1, 1, 1] = 3
    result = extract_mask_features(mask, (.01, .01, .01)).loc[3]
    assert result.volume_um3 == pytest.approx(1e-6, rel=1e-12)
    assert np.isnan(result.solidity)
    # A single voxel's discretized sphericity (~2.79) also exceeds the shared
    # review limit (feature schema 2.1); the value itself is not clipped.
    assert result.qc_flags == "solidity_not_estimable;sphericity_above_geometric_range"
    assert result.qc_status == "review"


def test_border_and_fragment_flags_preserve_objects_and_raw_mask():
    mask = np.zeros((12, 12, 12), np.uint32)
    mask[:3, :3, :3] = 9001
    mask[6:9, 6:9, 6:9] = 9001
    before = mask.copy()
    result = extract_mask_features(mask)
    assert list(result.index) == [9001]
    assert result.loc[9001, "volume_um3"] == 54
    assert result.loc[9001, "connected_component_count"] == 2
    assert result.loc[9001, "qc_flags"] == "border_truncated;fragmented_object"
    np.testing.assert_array_equal(mask, before)


def test_empty_result_retains_download_schema():
    result = extract_mask_features(np.zeros((3, 3, 3), np.uint8))
    assert result.empty and list(result.columns) == FEATURE_COLUMNS
    assert result.index.name == "Label"
    assert pd.read_csv(io.StringIO(result.to_csv())).columns.tolist() == ["Label", *FEATURE_COLUMNS]


def test_discretized_sphericity_above_review_limit_is_flagged_not_clipped():
    # A single voxel has sphericity ~2.79 (> 1 is impossible for a continuous
    # solid). The Web route must flag it like the classical route, not clip it.
    mask = np.zeros((5, 5, 5), np.uint8)
    mask[2, 2, 2] = 1
    result = extract_mask_features(mask).loc[1]
    assert result["sphericity"] > 2.5
    assert "sphericity_above_geometric_range" in result["qc_flags"].split(";")
    assert result["qc_status"] == "review"
