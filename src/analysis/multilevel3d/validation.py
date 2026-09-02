"""Data-contract validation for multilevel instance masks."""
from __future__ import annotations

import numpy as np


def validate_spacing(spacing_zyx_um: tuple[float, float, float]) -> tuple[float, float, float]:
    """Return finite positive spacing in canonical ``(Z, Y, X)`` order."""
    spacing = tuple(float(value) for value in spacing_zyx_um)
    if len(spacing) != 3 or not np.isfinite(spacing).all() or min(spacing) <= 0:
        raise ValueError("spacing_zyx_um must contain three positive finite values")
    return spacing


def validate_labels(name: str, labels: np.ndarray, expected_shape: tuple[int, int, int] | None = None) -> None:
    """Validate the shared label-mask contract without relabeling instances."""
    if not isinstance(labels, np.ndarray) or labels.ndim != 3:
        raise ValueError(f"{name} must be a 3D (Z, Y, X) numpy label mask")
    if expected_shape is not None and labels.shape != expected_shape:
        raise ValueError(f"{name} shape {labels.shape} does not match {expected_shape}")
    if labels.dtype == np.bool_ or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError(f"{name} must contain nonnegative integer instance IDs")
    if np.issubdtype(labels.dtype, np.signedinteger) and np.any(labels < 0):
        raise ValueError(f"{name} must contain nonnegative integer instance IDs")


def validate_inputs(
    organoid_labels: np.ndarray,
    cell_labels: np.ndarray,
    nucleus_labels: np.ndarray,
    spacing_zyx_um: tuple[float, float, float],
    nucleus_intensity: np.ndarray | None = None,
) -> tuple[float, float, float]:
    """Validate registered multilevel masks and optional raw nucleus channel."""
    validate_labels("organoid_labels", organoid_labels)
    validate_labels("cell_labels", cell_labels, organoid_labels.shape)
    validate_labels("nucleus_labels", nucleus_labels, organoid_labels.shape)
    if nucleus_intensity is not None:
        if not isinstance(nucleus_intensity, np.ndarray) or nucleus_intensity.ndim != 3:
            raise ValueError("nucleus_intensity must be a 3D (Z, Y, X) image")
        if nucleus_intensity.shape != organoid_labels.shape:
            raise ValueError("nucleus_intensity must match label-mask shape")
        if not np.issubdtype(nucleus_intensity.dtype, np.number) or not np.isfinite(nucleus_intensity).all():
            raise ValueError("nucleus_intensity must be finite numeric data")
    return validate_spacing(spacing_zyx_um)
