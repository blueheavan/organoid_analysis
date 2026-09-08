"""Streamlit + vtk.js GPU volume rendering for Z-stack fluorescence data.

Quick start
-----------
    streamlit run src/ui/vtk_viewer/demo_app.py

Public API
----------
- :func:`st_volume_viewer`        — render a list of per-channel volumes (or a
                                   single matrix) as a browser GPU volume.
- :class:`ChannelConfig`          — per-channel LUT / opacity / window.
- :func:`build_viewer_payload`    — low-level: build the renderable spec.
- :func:`render_viewer_html`      — low-level: produce the HTML string.
"""

from __future__ import annotations

import numpy as np

from .viewer_payload import (  # noqa: F401
    ChannelConfig,
    MAX_CHANNELS,
    ViewerSpec,
    build_mask_overlay_payload,
    build_surface_payload,
    build_viewer_payload,
    pack_channel_stack,
    render_viewer_html,
    spec_to_dict,
    validate_spacing,
)
from .rendering_presets import PRESETS  # noqa: F401

__all__ = [
    "ChannelConfig",
    "MAX_CHANNELS",
    "PRESETS",
    "ViewerSpec",
    "build_mask_overlay_payload",
    "build_surface_payload",
    "build_viewer_payload",
    "pack_channel_stack",
    "render_viewer_html",
    "spec_to_dict",
    "st_volume_viewer",
    "validate_spacing",
]


def st_volume_viewer(
    volumes,
    spacing_um=(1.0, 1.0, 1.0),
    channels=None,
    *,
    render_mode: str = "volume",
    segmentation_mask=None,
    surface_level=None,
    mask_overlay=None,
    mask_overlay_alpha: float = 0.5,
    height: int = 700,
    normalize: bool = True,
    compress: bool = False,
    use_st_iframe: bool = True,
):
    """Render a Z-stack as a browser GPU volume in Streamlit.

    Parameters
    ----------
    volumes : np.ndarray or sequence of np.ndarray
        Either a single ``(Z, Y, X)`` volume, or a list of per-channel volumes
        (each ``(Z, Y, X)``, up to ``MAX_CHANNELS = 4``).
    spacing_um : (x, y, z)
        Physical voxel spacing in microns, in X/Y/Z order.
    channels : sequence of ChannelConfig
        Per-channel LUT / opacity / window. Defaults to cyan/magenta/...
    render_mode : "volume" | "mip" | "surface"
        Defaults to ``"volume"`` (composite ray casting).
    segmentation_mask : np.ndarray, optional
        Required when ``render_mode == "surface"``.
    height : int
        Viewer pixel height passed to the iframe.
    """
    import streamlit as st

    if isinstance(volumes, np.ndarray):
        if volumes.ndim == 3:
            vols = [volumes]
        else:
            raise ValueError(
                "a numpy volume must be (Z, Y, X); pass a list for multi-channel"
            )
    else:
        vols = list(volumes)

    spec = build_viewer_payload(
        vols,
        spacing_um,
        channels,
        render_mode=render_mode,
        segmentation_mask=segmentation_mask,
        surface_level=surface_level,
        mask_overlay=mask_overlay,
        mask_overlay_alpha=mask_overlay_alpha,
        height_px=height,
        normalize=normalize,
        compress=compress,
    )
    html = render_viewer_html(spec)
    if use_st_iframe:
        try:
            st.iframe(html, height=height)
            return
        except Exception:
            pass
    from streamlit.components.v1 import html as _html

    _html(html, height=height, scrolling=False)
