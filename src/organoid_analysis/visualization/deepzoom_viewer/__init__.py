"""Full-resolution deep-zoom review of Z planes (OpenSeadragon, display only).

Public API
----------
- :func:`st_deepzoom_viewer`      -- render planes + contour overlay in Streamlit.
- :func:`build_deepzoom_payload`  -- tile a stack into per-plane pyramids.
- :func:`render_deepzoom_html`    -- standalone HTML string for a payload.
"""
from __future__ import annotations

from .deepzoom_payload import (
    DEFAULT_BUDGET_BYTES,
    TILE_SIZE,
    DeepZoomSpec,
    build_deepzoom_payload,
    level_shape,
    max_level_for,
    render_deepzoom_html,
    st_deepzoom_viewer,
)

__all__ = [
    "DEFAULT_BUDGET_BYTES",
    "DeepZoomSpec",
    "TILE_SIZE",
    "build_deepzoom_payload",
    "level_shape",
    "max_level_for",
    "render_deepzoom_html",
    "st_deepzoom_viewer",
]
