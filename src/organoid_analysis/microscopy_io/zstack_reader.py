"""Z-stack TIFF reader with explicit dimension handling.

Reads TIFF / multi-page / OME-TIFF and normalises the data into our canonical
layouts:

    * no channel:  volume.shape == (Z, Y, X)
    * channel:     volume.shape == (C, Z, Y, X)

We never assume "page order == Z". Axis semantics are resolved from the
TIFF series ``axes`` string / OME-XML DimensionOrder / ImageJ hyperstack
metadata. Time (T>1) is not spliced into the volume; it is carried for the UI
to pick a single time point. RGB data (S/Q axes or samples-per-pixel > 1) is
rejected rather than silently treated as grey.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import tifffile

from .metadata import Spacing, classify_axes, reject_rgb, resolve_spacing


@dataclass
class ZStack:
    """Normalised Z-stack volume and its declared metadata."""

    volume: np.ndarray            # (Z, Y, X) or (C, Z, Y, X)
    spacing: Spacing              # µm / voxel; fields may be None if unknown
    dtype: np.dtype
    axes: str                     # original tifffile series axes string
    channels: int
    frames: int                   # number of T frames (>=1)
    is_multichannel: bool
    meta: dict[str, Any] = field(default_factory=dict)

    @property
    def shape(self) -> tuple[int, ...]:
        return self.volume.shape


def read_axes(path: str | Path) -> str:
    """Return the axes string of the first available TIFF series."""
    with tifffile.TiffFile(path) as tf:
        if not tf.series:
            raise ValueError("TIFF contains no image series.")
        return tf.series[0].axes


def load_zstack(path: str | Path) -> ZStack:
    """Load a Z-stack TIFF and return a normalised :class:`ZStack`.

    Raises:
        ValueError if the axes cannot be resolved, if the image is RGB, or if
        spatial axes Y/X/Z cannot be located.
    """
    path = Path(path)
    series = None
    ome_metadata = None
    imagej_metadata = None
    tiff_tags = {}

    with tifffile.TiffFile(path) as tf:
        if not tf.series:
            raise ValueError("TIFF contains no image series.")
        series = tf.series[0]
        axes = series.axes or ""
        ome_metadata = tf.ome_metadata
        imagej_metadata = getattr(tf, "imagej_metadata", None)
        page = None
        for series_page in series.levels[0]:
            page = series_page
            break
        if page is not None:
            for key in ("XResolution", "YResolution", "ResolutionUnit"):
                tag = page.tags.get(key)
                if tag is not None:
                    tiff_tags[key] = tag.value
            samples = getattr(page, "samplesperpixel", None) or getattr(
                page, "samples_per_pixel", None
            )
        else:
            samples = None

        spacing = resolve_spacing(ome_metadata, imagej_metadata, tiff_tags)

        if reject_rgb(axes, samples):
            raise ValueError(
                "RGB/RGBA image detected (samples-per-pixel > 1). This is not "
                "quantitative Z-stack intensity data and cannot be rendered."
            )

        axis_map = classify_axes(axes)
        # Only spatial + C/T axes are allowed to participate; unknown filler
        # axes are rejected above for S/Q but we require Z/Y/X to exist.
        z_ax = axis_map["Z"]
        y_ax = axis_map["Y"]
        x_ax = axis_map["X"]
        if z_ax is None or y_ax is None or x_ax is None:
            raise ValueError(
                f"Cannot identify Z/Y/X axes from series axes {axes!r}. "
                "Refusing to guess dimension order."
            )

        data = series.asarray()
        data = _reorder_to_czyx(data, axes)

        # After reorder, dims are (T, C, Z, Y, X) minimal form.
        frames = int(data.shape[0])
        channels = int(data.shape[1])
        volume = data[0]  # (C, Z, Y, X)

        is_multichannel = channels > 1
        if not is_multichannel:
            volume = volume[0]  # (Z, Y, X)

    meta = {
        "ome_metadata": ome_metadata,
        "imagej_metadata": imagej_metadata,
        "tiff_tags": tiff_tags,
        "source": str(path),
    }
    return ZStack(
        volume=volume,
        spacing=spacing,
        dtype=volume.dtype,
        axes=axes,
        channels=channels,
        frames=frames,
        is_multichannel=is_multichannel,
        meta=meta,
    )


def _reorder_to_czyx(array: np.ndarray, axes: str) -> np.ndarray:
    """Permute a generic series array into canonical (T, C, Z, Y, X).

    Missing canonical axes are inserted as size-1 singletons. Unknown /
    unplaceable axes (e.g. RGB S/Q, or legacy filler dims) raise, because they
    cannot map to a quantitative Z-stack slot. Fully explicit, no guessing.
    """
    am = classify_axes(axes)
    present = {lbl: idx for lbl, idx in am.items() if idx is not None and lbl in "TCZYX"}

    unknown = set(range(array.ndim)) - set(present.values())
    if unknown:
        raise ValueError(
            f"Cannot map series axes {axes!r} to canonical (T,C,Z,Y,X); "
            "legacy/filler/RGB dimensions are not supported."
        )

    target = ["T", "C", "Z", "Y", "X"]
    # Step 1: bring present axes up front in their original relative order.
    src_order = [idx for idx in range(array.ndim) if idx in present.values()]
    arr = np.transpose(array, src_order)
    present_labels_in_src = [lbl for lbl, idx in present.items()]
    present_labels_in_src.sort(key=lambda lbl: present[lbl])
    pos = {lbl: k for k, lbl in enumerate(present_labels_in_src)}
    # Step 2: permute into target order among the present axes.
    perm = [pos[lbl] for lbl in target if lbl in present]
    arr = np.transpose(arr, perm) if len(perm) > 1 else arr
    # Step 3: insert size-1 for missing axes at their target slots.
    for i, lbl in enumerate(target):
        if lbl not in present:
            arr = np.expand_dims(arr, i)
    assert arr.ndim == 5, arr.shape
    return arr
