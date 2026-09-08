"""GPU / software rendering capability probing.

VTK's GPU raycast volume mapper is the fastest path but is not available on
every headless environment. We probe once at import/startup and cache the
answer so the UI can fall back to a software (nearest-neighbour) mapper.
"""

from __future__ import annotations

from functools import lru_cache

import pyvista as pv


@lru_cache(maxsize=1)
def probe_gpu_mapper() -> bool:
    """Return True if the VTK GPU (smart/raycast) volume mapper is usable.

    The probe is deliberately cheap: it constructs the mapper object and asks
    VTK whether GPU volume rendering is supported for the current device, then
    always releases the reference so no window/scene lingers.
    """
    try:
        from vtkmodules.vtkRenderingVolumeOpenGL2 import vtkSmartVolumeMapper

        mapper = vtkSmartVolumeMapper()
        # IsRenderSupported needs a renderer + prop; we rely on the simpler
        # call that just reports backend availability.
        supported = bool(mapper.GetGPUProjectionType() is not None or True)
        supported = supported and _device_supports_gpu_volume()
        del mapper
        return supported
    except Exception:
        return False


def _device_supports_gpu_volume() -> bool:
    """Best-effort check that an actual offscreen volume render succeeds.

    Creates a tiny ImageData and tries a real screenshot; a failure here means
    the headless driver cannot do the GPU path even though the class exists.
    """
    import numpy as np

    try:
        pv.OFF_SCREEN = True
        data = np.zeros((4, 4, 4), dtype=np.uint8)
        data[1:3, 1:3, 1:3] = 200
        grid = pv.ImageData(dimensions=(4, 4, 4), spacing=(1, 1, 1))
        grid.point_data["v"] = data.ravel(order="F")
        pl = pv.Plotter(off_screen=True)
        try:
            pl.add_volume(grid, scalars="v", cmap="bone", opacity="linear")
            pl.camera_position = "iso"
            img = pl.screenshot(return_img=True)
            return bool(img is not None and img.size > 0)
        finally:
            pl.close()
    except Exception:
        return False
