"""Direct grayscale projections of Z-stacks.

The projection modes operate directly on source pixels rather than through a
VTK colour/opacity transfer function. They therefore retain the black-to-white
appearance of a 2D microscope image instead of producing a heatmap.

Conventions
-----------
* The data array is (Z, Y, X) (single channel) or (C, Z, Y, X).
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from segmentation.io import Spacing

CAMERA_PRESETS: dict[str, str] = {
    "iso": "Isometric",
    "xy": "XY (top)",
    "xz": "XZ (front)",
    "yz": "YZ (side)",
}


def normalize_to_gray(
    volume: np.ndarray,
    low: float | None = None,
    high: float | None = None,
    low_pct: float = 1.0,
    high_pct: float = 99.0,
) -> np.ndarray:
    """Map intensities into a photo-like grayscale ``uint8`` in ``[0, 255]``.

    Clips to ``[low, high]`` (default the ``low_pct``..``high_pct`` percentiles
    of the data) and linearly rescales to 0..255. This reproduces a bright-field /
    fluorescence photo look for direct volume rendering: dark background -> 0,
    bright structures -> 255. Values outside the window saturate at the ends.
    """
    arr = np.asarray(volume, dtype=np.float32)
    if low is None or high is None:
        p_low, p_high = np.percentile(arr, (low_pct, high_pct))
        if low is None:
            low = p_low
        if high is None:
            high = p_high
    if not (np.isfinite(low) and np.isfinite(high)) or high <= low:
        low = float(arr.min())
        high = float(arr.max())
    if high <= low:
        return np.zeros(arr.shape, dtype=np.uint8)
    scale = 255.0 / (high - low)
    out = (arr - low) * scale
    out = np.clip(out, 0.0, 255.0)
    return out.astype(np.uint8)


@dataclass
class RenderResult:
    """A rendered offscreen image and its provenance."""

    image: np.ndarray          # (H, W, 3) uint8
    width: int
    height: int
    mode: str
    preset: str
    channel: int | None


def _project(volume: np.ndarray, preset: str, reducer: str) -> np.ndarray:
    """Project a ``(Z, Y, X)`` volume along the axis for a named view."""
    if preset in ("iso", "xy"):
        axis = 0
    elif preset == "xz":
        axis = 1
    elif preset == "yz":
        axis = 2
    else:
        raise ValueError(f"Unknown camera preset: {preset}")
    if reducer == "mean":
        return volume.mean(axis=axis)
    return volume.max(axis=axis)


def _as_grayscale_image(projection: np.ndarray,
                        norm: tuple[float | None, float | None]) -> np.ndarray:
    gray = normalize_to_gray(projection, low=norm[0], high=norm[1])
    return np.repeat(gray[..., np.newaxis], 3, axis=2)


def render_volume_image(
    volume: np.ndarray,
    spacing: Spacing | None = None,
    preset: str = "iso",
    channel: int | None = None,
    norm: tuple[float | None, float | None] = (None, None),
) -> RenderResult:
    """Render a direct grayscale density projection.

    ``volume`` may be (Z, Y, X) (single channel) or (C, Z, Y, X)
    (multichannel); when multichannel, ``channel`` selects which slice to draw.
    Voxel intensities are averaged along the selected viewing axis, then
    converted to grayscale. This is a 2D density image, not VTK volume rendering.
    """
    if volume.ndim == 4:
        if channel is None:
            raise ValueError("channel is required for a multichannel volume")
        volume = volume[channel]
    image = _as_grayscale_image(_project(np.asarray(volume), preset, "mean"), norm)
    return RenderResult(image, image.shape[1], image.shape[0], "volume", preset, channel)


def render_mip_image(
    volume: np.ndarray,
    spacing: Spacing | None = None,
    preset: str = "iso",
    channel: int | None = None,
    norm: tuple[float | None, float | None] = (None, None),
) -> RenderResult:
    """Render a direct grayscale maximum-intensity projection.

    Unlike VTK maximum blending, this calculates the maximum from the source
    pixels and returns it as an RGB grayscale image. It therefore has the same
    intensity appearance as a 2D microscope image, without a colour transfer
    function or volume-rendering artefacts.
    """
    if volume.ndim == 4:
        if channel is None:
            raise ValueError("channel is required for a multichannel volume")
        volume = volume[channel]
    image = _as_grayscale_image(_project(np.asarray(volume), preset, "max"), norm)
    return RenderResult(image, image.shape[1], image.shape[0], "mip", preset, channel)
