"""Validation and conversion helpers for physical voxel spacing."""

from __future__ import annotations

import math
from collections.abc import Sequence


def validate_voxel_spacing_xyz(
    spacing_um: Sequence[float],
) -> tuple[float, float, float]:
    """Return validated ``(x, y, z)`` spacing in micrometres."""
    if spacing_um is None or len(spacing_um) != 3:
        raise ValueError(
            f"spacing_um must have exactly 3 values (x, y, z), got {spacing_um!r}"
        )
    spacing = tuple(float(value) for value in spacing_um)
    if any(not math.isfinite(value) or value <= 0 for value in spacing):
        raise ValueError(f"spacing_um entries must be finite and > 0: {spacing_um!r}")
    return spacing  # type: ignore[return-value]


def isotropic_xy_size_um(x_um: float, y_um: float) -> float:
    """Return the shared XY size, rejecting unsupported rectangular pixels."""
    x, y, _ = validate_voxel_spacing_xyz((x_um, y_um, 1.0))
    if not math.isclose(x, y, rel_tol=1e-6, abs_tol=1e-9):
        raise ValueError(
            "This segmentation route requires equal X/Y pixel sizes; "
            f"metadata reports x={x:g} um and y={y:g} um. Resample explicitly "
            "or use a workflow that preserves independent X/Y spacing."
        )
    return x
