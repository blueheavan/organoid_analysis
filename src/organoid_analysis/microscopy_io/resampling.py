"""Physical grids for center-aligned resampling of ZYX arrays."""
from __future__ import annotations

from collections.abc import Sequence
from numbers import Integral, Real

import numpy as np

from .voxel_spacing import validate_voxel_spacing_xyz


def _shape_zyx(shape: Sequence[int]) -> tuple[int, int, int]:
    if len(shape) != 3 or any(isinstance(n, bool) or not isinstance(n, Integral) or n < 1 for n in shape):
        raise ValueError("Resampling requires a nonempty ZYX shape")
    return int(shape[0]), int(shape[1]), int(shape[2])


def xy_downsample_shape(shape: Sequence[int], factor: float) -> tuple[int, int, int]:
    """Return the actual rounded shape used by scipy.ndimage.zoom."""
    z, y, x = _shape_zyx(shape)
    if isinstance(factor, bool) or not isinstance(factor, Real) or not np.isfinite(factor) or not 0 < factor <= 1:
        raise ValueError("xy_downsample must be finite and in (0, 1]")
    target = (z, round(y * factor), round(x * factor))
    if min(target[1:]) < 2:
        raise ValueError("xy_downsample would leave fewer than two XY pixels")
    return target


def resampled_spacing_xyz(
    source_shape: Sequence[int], target_shape: Sequence[int], spacing_xyz: Sequence[float],
) -> tuple[float, float, float]:
    """Preserve first/last voxel-center coordinates for ``grid_mode=False``.

    A length-N axis spans N-1 center intervals, not N voxel widths. Singleton
    axes can only be retained unchanged because their span cannot define a scale.
    """
    source = _shape_zyx(source_shape)
    target = _shape_zyx(target_shape)
    spacing = validate_voxel_spacing_xyz(spacing_xyz)
    ratios = []
    for original, resized in zip(source, target):
        if original == resized:
            ratios.append(1.0)
        elif min(original, resized) < 2:
            raise ValueError("Cannot rescale a singleton axis while preserving voxel-center extent")
        else:
            ratios.append((original - 1) / (resized - 1))
    return spacing[0] * ratios[2], spacing[1] * ratios[1], spacing[2] * ratios[0]


def isotropic_xy_scale(source_shape: Sequence[int], target_shape: Sequence[int]) -> float:
    """Effective pixel scale for a model accepting only one XY pixel size.

    Compare rational grid scales exactly; a scalar anisotropy cannot encode two
    different XY spacings. Full-resolution inference always has unit scale.
    """
    z, y, x = _shape_zyx(source_shape)
    rz, ry, rx = _shape_zyx(target_shape)
    if z != rz or min(y, x, ry, rx) < 2:
        raise ValueError("Cellpose XY resampling must preserve Z and at least two XY pixels")
    if (ry - 1) * (x - 1) != (rx - 1) * (y - 1):
        raise ValueError(
            "Rounded downsampling produces unequal X/Y pixel spacings that Cellpose's "
            "scalar anisotropy cannot represent. Use full XY resolution (xy_downsample=1)."
        )
    return (rx - 1) / (x - 1)
