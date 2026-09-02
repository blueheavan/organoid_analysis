"""Workaround for stpyvista failing to import under Streamlit >= 1.61.

stpyvista 0.2.1 declares its frontend assets as file paths
(``js="v2/stpyvista.js"``, ``css="v2/style.css"``) in
``streamlit.components.v2.component(...)``. Since Streamlit 1.61 file-backed
v2 component assets must be declared in the package ``pyproject.toml`` with an
``asset_dir``, and stpyvista's own manifest does not declare one, ``import
stpyvista`` raises ``StreamlitAPIException``.

This module patches ``streamlit.components.v2.component`` *before* stpyvista is
imported so that file-backed ``js``/``css`` are eagerly read from the installed
stpyvista package and passed as inline content. Multiline inline strings are
always treated as inline by Streamlit, so nothing else changes.
"""

from __future__ import annotations

from pathlib import Path

import streamlit.components.v2 as _scv

_ORIGINAL_COMPONENT = _scv.component


def _stpv_asset_dir() -> Path:
    import stpyvista

    return Path(stpyvista.__file__).parent / "backends"


def _load_asset(rel: str) -> str:
    """Read a stpyvista frontend asset identified by a package-relative path."""
    base = _stpv_asset_dir()
    path = base / rel
    if not path.exists():
        raise FileNotFoundError(f"stpyvista asset not found: {path}")
    return path.read_text(encoding="utf-8")


def _looks_like_file_path(value: str) -> bool:
    s = value.strip()
    if "\n" in s or "\r" in s:
        return False
    return ".js" in s or ".css" in s or "/" in s or "\\" in s or s.startswith(".")


def _patched_component(*args, **kwargs):
    """Resolve file-backed js/css to inline before calling the real component."""
    for key in ("js", "css"):
        value = kwargs.get(key)
        if isinstance(value, str) and _looks_like_file_path(value):
            try:
                kwargs[key] = _load_asset(value)
            except FileNotFoundError:
                # Not one of ours (another package's asset); leave untouched.
                pass
    return _ORIGINAL_COMPONENT(*args, **kwargs)


def apply_stpv_patch() -> None:
    """Install the patch. Idempotent."""
    if _scv.component is _patched_component:
        return
    _scv.component = _patched_component


def stpyvista_surface(surface_dicts: list[dict]) -> None:
    """Render interactive surfaces with stpyvista (best-effort).

    ``surface_dicts`` is a list of ``surface_to_dict(...)`` payloads. If
    stpyvista cannot be satisfied at runtime, an ImportError is raised by the
    caller and the offscreen fallback should be used instead.
    """
    import pyvista as pv
    from stpyvista import stpyvista

    plotter = pv.Plotter(window_size=(768, 512))
    try:
        for sdict in surface_dicts:
            poly = sdict["polydata"]
            # Rebuild PolyData from the serialized points/faces.
            pts = _unpack_points(poly["points"], poly["n_points"])
            faces = _unpack_faces(poly["faces"])
            mesh = pv.PolyData(pts, faces)
            color = sdict["color"]
            plotter.add_mesh(
                mesh,
                color=tuple(c / 255.0 for c in color) if color else "red",
                opacity=sdict["opacity"],
                smooth_shading=True,
            )
        plotter.enable_parallel_projection()
        # Panel backend serializes PolyData faithfully (verified in the spike)
        # and, unlike trame, does not require launching an asyncio server.
        stpyvista(plotter, backend="panel", key="stpv_surface")
    finally:
        plotter.close()


def _unpack_points(flat: list[float], n: int):
    import numpy as np

    return np.asarray(flat, dtype=np.float32).reshape(n, 3)


def _unpack_faces(flat: list[int]):
    import numpy as np

    tri = np.asarray(flat, dtype=np.int64).reshape(-1, 3)
    faces = np.empty((tri.shape[0], 4), dtype=np.int64)
    faces[:, 0] = 3
    faces[:, 1:] = tri
    return faces
