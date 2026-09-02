"""Streamlit demo for the vtk.js GPU volume viewer.

Run:
    streamlit run src/ui/vtk_viewer/demo_app.py

Renders a synthetic 3-channel Z-stack (DAPI/marker/cytoplasm) as a browser GPU
volume by default. Optionally upload a real Z-stack TIFF.
"""

from __future__ import annotations

import numpy as np
import streamlit as st

from ui.vtk_viewer import ChannelConfig, st_volume_viewer

st.set_page_config(page_title="vtk.js Volume Viewer", layout="wide")


def synthetic_organoid(shape=(96, 128, 128), seed: int = 0) -> np.ndarray:
    """A soft spheroid of organoid-like intensity on a dark background."""
    rng = np.random.default_rng(seed)
    z, y, x = shape
    zz, yy, xx = np.mgrid[0:z, 0:y, 0:x]
    cz, cy, cx = z / 2, y / 2, x / 2
    dist = np.sqrt(((zz - cz) / (z * 0.45)) ** 2 + ((yy - cy) / (y * 0.45)) ** 2
                   + ((xx - cx) / (x * 0.45)) ** 2)
    shell = np.clip(1.0 - dist, 0, 1) ** 2
    noise = rng.normal(0, 0.03, shape)
    hollow = np.clip(1.0 - dist, 0, 1) ** 2
    # carve a central lumen so rays see a shell, like a real organoid
    lumen = (np.sqrt((zz - cz) ** 2 + (yy - cy) ** 2 + (xx - cx) ** 2) / (z * 0.20))
    lumen = np.clip(1 - lumen, 0, 1)
    vol = shell * (0.25 + 0.75 * lumen)
    vol = np.clip(vol + noise, 0, 1)
    return (vol * 2000).astype(np.uint16)


def main() -> None:
    st.title("vtk.js GPU Volume Viewer")
    st.caption(
        "Browser-side GPU direct volume ray-casting of a Z-stack. "
        "Rotate/zoom/pan with the mouse; toggle modes and channel opacity in the panel."
    )

    source = st.radio(
        "Data source",
        ["Synthetic 3-channel", "Upload Z-stack TIFF"],
        horizontal=True,
    )

    if source == "Upload Z-stack TIFF":
        upload = st.file_uploader("Z-stack TIFF", type=["tif", "tiff"])
        if upload is None:
            st.info("Upload a TIFF to render it.")
            return
        import tifffile

        zs = tifffile.imread(upload)
        if zs.ndim == 3:
            vols = [zs.astype(np.float32)]
        elif zs.ndim == 4 and zs.shape[-1] in (3, 4):
            raise ValueError("RGB(A) TIFF is not quantitative; refusing to render.")
        else:
            raise ValueError(f"Unsupported TIFF shape {zs.shape}; need (Z, Y, X).")
        z_spacing = st.slider("Z voxel spacing (µm)", 0.5, 5.0, 1.5, 0.1)
        spacing = (0.414, 0.414, z_spacing)
        channels = [ChannelConfig(lut="cyan", opacity=0.9, name="TIFF")]
    else:
        vols = [
            synthetic_organoid(seed=1).astype(np.float32),
            (synthetic_organoid(seed=2) * 0.6).astype(np.float32),
            (synthetic_organoid(seed=3) * 0.35).astype(np.float32),
        ]
        spacing = (0.65, 0.65, 2.0)
        channels = [
            ChannelConfig(lut="cyan", opacity=0.9, name="DAPI"),
            ChannelConfig(lut="magenta", opacity=0.7, name="Marker"),
            ChannelConfig(lut="yellow", opacity=0.5, name="Cytoplasm"),
        ]

    st.write(
        f"**Volume**: {vols[0].shape} ({vols[0].dtype}) | "
        f"**Channels**: {len(vols)} | **Spacing (µm)**: x={spacing[0]}, y={spacing[1]}, z={spacing[2]}"
    )

    st_volume_viewer(
        vols,
        spacing_um=spacing,
        channels=channels,
        render_mode="volume",
        height=700,
    )


main()
