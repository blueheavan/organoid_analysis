"""Surface extraction for 3D preview (binary / label volumes).

The interactive path serializes the surface (``surface_to_dict``) for the
browser-side vtk.js viewer (``visualization.volume_viewer``); the fallback
renders the same polydata offscreen into a static image. Both produce the
geometry from marching cubes, so they stay visually consistent.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pyvista as pv

from organoid_analysis.microscopy_io import Spacing


@dataclass
class Surface:
    """An extracted surface mesh and its display metadata."""

    polydata: pv.PolyData
    label: int | None
    color: tuple[int, int, int] | None
    opacity: float = 1.0


def _poly_repr(pd: pv.PolyData) -> dict:
    pts = np.asarray(pd.points)
    tris = np.asarray(pd.faces).reshape(-1, 4)[:, 1:] if pd.n_faces else np.zeros((0, 3))
    return {
        "points": pts.astype(np.float32).ravel(order="C").tolist(),
        "faces": tris.astype(np.uint32).ravel(order="C").tolist(),
        "n_points": int(pts.shape[0]),
        "n_faces": int(tris.shape[0]),
        "bounds": list(float(v) for v in pd.bounds),
        "has_scalars": bool("label" in pd.array_names),
        "scalars": (
            np.asarray(pd["label"]).astype(np.uint32).tolist()
            if "label" in pd.array_names
            else None
        ),
    }


def build_surface(
    volume: np.ndarray,
    spacing: Spacing | None = None,
    threshold: float | None = None,
    labels: Sequence[int] | None = None,
) -> list[Surface]:
    """Extract surfaces from a (Z, Y, X) volume.

    * ``labels`` given -> treat ``volume`` as a label map and contour each
      requested label (each becomes its own mesh, coloured by label).
    * ``labels`` None + ``threshold`` given -> binary iso-surface at the level.
    * ``labels`` None and ``threshold`` None -> contour at ``volume.max()/2``.
    """
    spacing = spacing or Spacing(None, None, None)
    dx = spacing.x or 1.0
    dy = spacing.y or 1.0
    dz = spacing.z or 1.0
    z, y, x = volume.shape
    grid = pv.ImageData(dimensions=(x, y, z), spacing=(dx, dy, dz))
    grid.point_data["scalar"] = np.ascontiguousarray(
        volume.ravel(order="C"), dtype=volume.dtype
    )

    surfaces: list[Surface] = []
    if labels is not None:
        colors = _label_colors(len(labels))
        for i, label in enumerate(labels):
            mask = np.asarray(volume == label, dtype=volume.dtype)
            lab_grid = pv.ImageData(dimensions=(x, y, z), spacing=(dx, dy, dz))
            lab_grid.point_data["scalar"] = np.ascontiguousarray(
                mask.ravel(order="C"), dtype=volume.dtype
            )
            try:
                mesh = lab_grid.contour([0.5], scalars="scalar")
            except Exception:
                continue
            mesh["label"] = np.full(mesh.n_points, label, dtype=np.uint32)
            surfaces.append(
                Surface(mesh, int(label), colors[i % len(colors)], opacity=0.85)
            )
    else:
        level = threshold if threshold is not None else (float(volume.max()) / 2 if volume.max() > 0 else 0.5)
        try:
            mesh = grid.contour([level], scalars="scalar")
        except Exception:
            mesh = pv.PolyData()
        surfaces.append(Surface(mesh, None, None, opacity=1.0))
    return surfaces


def _label_colors(n: int) -> list[tuple[int, int, int]]:
    palette = [
        (231, 76, 60), (46, 204, 113), (52, 152, 219),
        (241, 196, 15), (155, 89, 182), (26, 188, 156),
        (243, 156, 18), (127, 140, 141),
    ]
    return [palette[i % len(palette)] for i in range(n)]


def surface_to_dict(surface: Surface) -> dict:
    """Serialize a Surface into a plain dict consumable by the vtk.js viewer."""
    return {
        "polydata": _poly_repr(surface.polydata),
        "label": surface.label,
        "color": surface.color,
        "opacity": surface.opacity,
    }


def render_surface_image(
    surfaces: Sequence[Surface],
    preset: str = "iso",
    size: tuple[int, int] = (768, 768),
    color: tuple[int, int, int] = (231, 76, 60),
) -> np.ndarray:
    """Render extracted surfaces offscreen into a static RGBA image."""
    pv.OFF_SCREEN = True
    pl = pv.Plotter(off_screen=True, window_size=size)
    try:
        for surf in surfaces:
            color_use = surf.color or color
            pl.add_mesh(
                surf.polydata,
                color=tuple(c / 255.0 for c in color_use),
                opacity=surf.opacity,
                smooth_shading=True,
            )
        if preset == "xz":
            pl.camera_position = "xz"
        elif preset == "yz":
            pl.camera_position = "yz"
        elif preset == "xy":
            pl.camera_position = "xy"
        else:
            pl.camera_position = "iso"
        pl.enable_parallel_projection()
        pl.render()
        return pl.screenshot(return_img=True)
    finally:
        pl.close()
