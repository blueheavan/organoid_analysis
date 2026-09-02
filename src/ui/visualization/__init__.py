"""3D preview rendering for Z-stack TIFFs.

Implements the "Route B (full)" design:

* Volume / MIP -> server-side offscreen VTK rendering into a numpy image
  (``volume_3d``), shown with ``st.image`` and repositioned only via camera
  presets (no free rotation).
* Surface -> interactive stpyvista widget (``surface``) for free 3D rotation.
"""

from .capabilities import probe_gpu_mapper
from .volume_3d import (
    render_volume_image,
    render_mip_image,
    normalize_to_gray,
    CAMERA_PRESETS,
)
from .surface import build_surface, render_surface_image, surface_to_dict

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
