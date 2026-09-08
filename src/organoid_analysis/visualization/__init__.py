"""3D preview rendering for Z-stack TIFFs.

Implements the "Route B (full)" design:

* Volume / MIP -> server-side offscreen VTK rendering into a numpy image
  (``projection_rendering``), shown with ``st.image`` and repositioned only
  via camera presets (no free rotation).
* Surface -> interactive stpyvista widget (``surface_rendering``) for free 3D
  rotation.
"""

from .capabilities import probe_gpu_mapper
from .projection_rendering import (
    CAMERA_PRESETS,
    normalize_to_gray,
    render_mip_image,
    render_volume_image,
)
from .surface_rendering import build_surface, render_surface_image, surface_to_dict

__all__ = [
    "probe_gpu_mapper",
    "render_volume_image",
    "render_mip_image",
    "normalize_to_gray",
    "CAMERA_PRESETS",
    "build_surface",
    "render_surface_image",
    "surface_to_dict",
]
