"""Regression tests for physical voxel-spacing propagation.

Guards against the specific failure mode this module was added to close: a
TIFF's real X/Y/Z pixel size silently getting dropped in favour of the
hardcoded ``SegmentationConfig`` fallback (0.414 um XY, 2.9 anisotropy), which
would make every downstream volume/area/distance figure wrong without any
visible error.
"""

from __future__ import annotations

import numpy as np
import pytest

from organoid_analysis.microscopy_io import isotropic_xy_size_um, resolve_spacing_source
from organoid_analysis.quantification.mask_features import extract_mask_features
from organoid_analysis.segmentation.cellpose_inference import SegmentationConfig
from organoid_analysis.web_interface.segmentation_workspace import config_spacing

# --- resolve_spacing_source ------------------------------------------------


def test_matches_metadata_is_classified_as_metadata():
    assert resolve_spacing_source(0.65, 0.65, 0.414) == "metadata"


def test_no_metadata_and_matches_fallback_is_classified_as_default():
    assert resolve_spacing_source(0.414, None, 0.414) == "default"


def test_differs_from_metadata_is_classified_as_user_override():
    assert resolve_spacing_source(0.70, 0.65, 0.414) == "user_override"


def test_no_metadata_and_differs_from_fallback_is_classified_as_user_override():
    assert resolve_spacing_source(0.50, None, 0.414) == "user_override"


# --- SegmentationConfig provenance defaults --------------------------------


def test_bare_config_reports_default_source_and_no_metadata_value():
    config = SegmentationConfig()
    assert config.xy_spacing_source == "default"
    assert config.anisotropy_source == "default"
    assert config.metadata_xy_spacing_um is None
    assert config.metadata_anisotropy is None


def _solid_block_mask(shape_zyx: tuple[int, int, int]) -> np.ndarray:
    mask = np.zeros(shape_zyx, dtype=np.uint16)
    mask[:, :, :] = 1
    return mask


# --- Case 1: real metadata spacing propagates to volume --------------------
# OME-TIFF reports X=0.65um, Y=0.65um, Z=2.0um -- downstream volume must use
# 0.65 * 0.65 * 2.0 per voxel, never the 0.414 um default.


def test_voxel_volume_uses_real_metadata_spacing_not_the_0_414_default():
    x_um = y_um = 0.65
    z_um = 2.0
    anisotropy = z_um / isotropic_xy_size_um(x_um, y_um)

    config = SegmentationConfig(
        xy_spacing_um=x_um,
        xy_spacing_source="metadata",
        metadata_xy_spacing_um=x_um,
        anisotropy=anisotropy,
        anisotropy_source="metadata",
        metadata_anisotropy=anisotropy,
    )
    spacing_xyz = config_spacing(config)
    assert spacing_xyz == pytest.approx((0.65, 0.65, 2.0))

    shape_zyx = (10, 6, 6)
    mask = _solid_block_mask(shape_zyx)
    features = extract_mask_features(mask, spacing_um=spacing_xyz)

    expected_voxel_volume_um3 = 0.65 * 0.65 * 2.0
    expected_total_um3 = expected_voxel_volume_um3 * mask.size
    wrong_default_voxel_volume_um3 = 0.414 * 0.414 * (0.414 * 2.9)

    assert features.loc[1, "volume_um3"] == pytest.approx(expected_total_um3, rel=1e-3)
    # The two per-voxel volumes must not be anywhere near each other -- this
    # is the concrete number a silent fallback-to-default bug would have
    # produced instead.
    assert expected_voxel_volume_um3 == pytest.approx(0.845, rel=1e-3)
    assert wrong_default_voxel_volume_um3 == pytest.approx(0.2054, rel=1e-2)
    assert abs(expected_voxel_volume_um3 - wrong_default_voxel_volume_um3) > 0.5


# --- Case 2: missing metadata falls back explicitly -------------------------


def test_no_metadata_falls_back_to_default_with_explicit_source():
    backend_defaults = SegmentationConfig()
    metadata_xy_spacing_um = None  # e.g. no OME-XML / ImageJ spacing found
    source = resolve_spacing_source(
        backend_defaults.xy_spacing_um, metadata_xy_spacing_um, backend_defaults.xy_spacing_um
    )
    config = SegmentationConfig(
        xy_spacing_um=backend_defaults.xy_spacing_um,
        xy_spacing_source=source,
        metadata_xy_spacing_um=metadata_xy_spacing_um,
    )
    assert config.xy_spacing_source == "default"
    assert config.metadata_xy_spacing_um is None


# --- Case 3: user override is recorded, metadata original preserved --------


def test_override_differing_from_metadata_keeps_original_metadata_value():
    metadata_xy_spacing_um = 0.65
    user_typed_value = 0.80  # user edits the sidebar widget after auto-detect
    source = resolve_spacing_source(
        user_typed_value, metadata_xy_spacing_um, SegmentationConfig().xy_spacing_um
    )
    config = SegmentationConfig(
        xy_spacing_um=user_typed_value,
        xy_spacing_source=source,
        metadata_xy_spacing_um=metadata_xy_spacing_um,
    )
    assert config.xy_spacing_source == "user_override"
    # The metadata's original value must survive the override for audit, not
    # be silently discarded once the user disagrees with it.
    assert config.metadata_xy_spacing_um == pytest.approx(0.65)
    assert config.xy_spacing_um == pytest.approx(0.80)


# --- Case 4: unequal X/Y is rejected, never silently averaged --------------


def test_unequal_xy_pixel_size_raises_instead_of_averaging():
    with pytest.raises(ValueError, match="equal X/Y pixel sizes"):
        isotropic_xy_size_um(0.5, 0.6)
