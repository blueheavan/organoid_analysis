"""Analytical coordinate oracles for center-grid interpolation."""
import numpy as np
import pytest

from organoid_analysis.microscopy_io.resampling import (
    isotropic_xy_scale,
    resampled_spacing_xyz,
    xy_downsample_shape,
)


@pytest.mark.parametrize("shape,factor,target", [
    ((3, 16, 16), .5, (3, 8, 8)),
    ((3, 17, 17), .5, (3, 8, 8)),
    ((3, 15, 15), .25, (3, 4, 4)),
    ((3, 15, 19), 1., (3, 15, 19)),
])
def test_grid_preserves_center_spans(shape, factor, target):
    assert xy_downsample_shape(shape, factor) == target
    spacing = (0.3, 0.3, 1.7)
    resized = resampled_spacing_xyz(shape, target, spacing)
    np.testing.assert_allclose((np.array(shape[::-1]) - 1) * spacing,
                               (np.array(target[::-1]) - 1) * resized, rtol=1e-12)
    assert isotropic_xy_scale(shape, target) == (target[2] - 1) / (shape[2] - 1)


def test_anisotropic_preview_can_represent_separate_xy_scales():
    assert resampled_spacing_xyz((3, 11, 17), (3, 4, 6), (2, 3, 4)) == (6.4, 10., 4.)
    with pytest.raises(ValueError, match="unequal X/Y"):
        isotropic_xy_scale((3, 11, 17), (3, 4, 6))


@pytest.mark.parametrize("factor", [0, -1, 2, np.nan, np.inf, True, .001])
def test_invalid_or_collapsed_grid_rejects(factor):
    with pytest.raises(ValueError):
        xy_downsample_shape((3, 8, 8), factor)


def test_singleton_z_can_be_retained_but_not_rescaled():
    assert resampled_spacing_xyz((1, 5, 5), (1, 3, 3), (1, 1, 2)) == (2, 2, 2)
    with pytest.raises(ValueError, match="singleton"):
        resampled_spacing_xyz((1, 5, 5), (2, 3, 3), (1, 1, 2))
