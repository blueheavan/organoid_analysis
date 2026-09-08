"""3D preview page for Z-stack TIFFs (Route B full version).

Run with:  pixi run preview
(or)        streamlit run src/ui/preview_2d.py

Rendering modes:
* Density -> direct mean-intensity projection of source pixels -> static image
* MIP     -> direct maximum-intensity projection of source pixels -> static image
* Surface -> interactive stpyvista widget (with import workaround); falls back
             to an offscreen static render when the widget is unavailable.

Volume/MIP views are repositioned via named camera presets (XY / XZ / YZ /
isometric); free rotation is only available in the interactive Surface mode.
"""

from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import tempfile

import numpy as np
import streamlit as st

from organoid_analysis.microscopy_io import load_zstack
from organoid_analysis.microscopy_io.metadata import Spacing
from organoid_analysis.visualization import build_surface, render_mip_image, render_surface_image, render_volume_image
from organoid_analysis.visualization._stpv import apply_stpv_patch

st.set_page_config(page_title="3D Preview", layout="wide")


@st.cache_data(show_spinner="Loading Z-stack...")
def _load(path: str) -> dict:
    zs = load_zstack(path)
    return {
        "volume": zs.volume,
        "spacing": asdict(zs.spacing),
        "dtype": str(zs.dtype),
        "axes": zs.axes,
        "channels": zs.channels,
        "frames": zs.frames,
        "is_multichannel": zs.is_multichannel,
    }


def _stack_from_upload(upload) -> Path:
    suffix = Path(upload.name).suffix or ".tif"
    with tempfile.NamedTemporaryFile(prefix="3dpreview_", suffix=suffix, delete=False) as tmp:
        tmp.write(upload.getbuffer())
        return Path(tmp.name)


def _render_single(volume, spacing: Spacing, mode: str, opts: dict):
    if mode in ("density", "mip"):
        fn = render_volume_image if mode == "density" else render_mip_image
        return fn(
            volume,
            spacing=spacing,
            preset=opts["preset"],
            channel=opts.get("channel"),
            norm=opts.get("norm", (None, None)),
        ).image


def main() -> None:
    st.title("3D Z-stack Preview")
    st.caption(
        "Density and MIP show original-pixel grayscale projections; "
        "Surface is interactive 3D."
    )

    upload = st.file_uploader(
        "Z-stack TIFF (single-channel or multichannel)", type=["tif", "tiff"]
    )
    if upload is None:
        st.info("Upload a Z-stack TIFF to preview it in 3D.")
        return

    data = _load(str(_stack_from_upload(upload)))
    volume = data["volume"]
    spacing = Spacing(**data["spacing"])
    is_multi = data["is_multichannel"]

    col_meta, col_controls = st.columns([1, 2])
    with col_meta:
        st.subheader("Metadata")
        st.write(f"Original axes: `{data['axes']}`")
        st.write(f"Dtype: `{data['dtype']}`")
        x, y, z = spacing.x, spacing.y, spacing.z
        st.write(
            f"Spacing (µm/voxel): X={x if x is not None else '?'}  "
            f"Y={y if y is not None else '?'}  Z={z if z is not None else '?'}"
        )
        st.write(f"Channels: {data['channels']}  |  Time frames: {data['frames']}")
        if is_multi:
            st.write(f"Multichannel volume shape: {volume.shape}")

    with col_controls:
        st.subheader("Display controls")
        mode = st.selectbox("Render mode", ["Density Projection", "MIP", "Surface"])
        preset = st.selectbox(
            "Projection direction", ["xy", "xz", "yz"],
            format_func={"xy": "XY (top)", "xz": "XZ (front)", "yz": "YZ (side)"}.get,
        ) if mode in ("Density Projection", "MIP") else "iso"
        channel = None
        if is_multi:
            channel = st.slider("Channel", 0, data["channels"] - 1, 0)
        downsample = st.select_slider(
            "XY-downsample (speed)", options=[1.0, 0.5, 0.25], value=1.0
        )
        st.caption("Density and MIP use a fixed black-to-white grayscale palette.")
        lo, hi = int(volume.min()), int(volume.max())
        p1, p99 = int(np.percentile(volume, 1)), int(np.percentile(volume, 99))
        if hi > lo:
            i_lo, i_hi = st.slider(
                "Intensity range (min / max)",
                lo, hi, (p1, p99 if p99 > p1 else hi),
            )
        else:
            i_lo, i_hi = lo, hi

    # Apply downsample consistently with spacing.
    work_volume = volume
    work_spacing = spacing
    if downsample < 1.0 and not is_multi:
        from scipy import ndimage

        down = downsample
        factors = (1.0, down, down) if volume.ndim == 3 else (1.0, 1.0, down, down)
        work = ndimage.zoom(volume.astype(np.float32), factors, order=1)
        work_volume = work.astype(np.uint16) if volume.dtype.itemsize >= 2 else work.astype(np.uint8)
        # XY spacing shrinks when XY is downsampled; Z unchanged.
        work_spacing = Spacing(
            spacing.x / down if spacing.x else None,
            spacing.y / down if spacing.y else None,
            spacing.z,
        )

    st.divider()

    if mode == "Surface":
        st.subheader("Surface mesh")
        thresh = st.slider(
            "Contour (iso) threshold",
            int(volume.min()), int(max(volume.max(), volume.min() + 1)),
            int((volume.min() + volume.max()) / 2),
        )
        surfs = build_surface(
            work_volume if work_volume.ndim == 3 else work_volume[channel],
            work_spacing, threshold=thresh,
        )
        if not surfs:
            st.warning("No surface found at this threshold.")
            return
        try:
            apply_stpv_patch()
            from organoid_analysis.visualization._stpv import stpyvista_surface
            from organoid_analysis.visualization import surface_to_dict

            stpyvista_surface([surface_to_dict(s) for s in surfs])
        except Exception as exc:  # interactive unavailable -> offscreen fallback
            st.caption("Interactive surface unavailable; showing static render.")
            img = render_surface_image(surfs, preset=preset)
            st.image(img, width="stretch")
    else:
        opts = {
            "preset": preset,
            "channel": channel if is_multi else None,
            "norm": (i_lo, i_hi),
        }
        try:
            img = _render_single(
                work_volume,
                work_spacing,
                "density" if mode == "Density Projection" else "mip",
                opts,
            )
            st.image(img, width="stretch")
        except Exception as exc:
            st.error(f"Render failed: {exc}")


main()
