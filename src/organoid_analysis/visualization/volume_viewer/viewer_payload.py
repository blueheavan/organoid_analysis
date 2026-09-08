"""Core data packing and HTML/JS generation for the vtk.js volume viewer.

Architecture
------------
Streamlit
  -> Python: validate + pack a Z-stack as x-fastest ``(Z,Y,X,C)`` float32 wire
     data, base64-encode it, and generate a standalone HTML string served
     through ``st.iframe``.
  -> browser: vanilla JS decodes the payload into one multi-component
     ``vtkImageData`` with ``dimensions = (X, Y, Z)`` and independent per-channel
     transfer functions, then shows the same data in 2D and GPU 3D viewports.

The vendored ``vtk.js`` UMD bundle (``static/vtk.js``) is inlined into the
generated HTML at render time so the iframe has no external / CDN dependency.
This keeps the "vendor locally, no npm build" constraint while remaining robust
in an iframe (an HTML-string iframe cannot resolve relative ``<script src>``
paths to files on the Streamlit server).

All render controls live inside the viewer HTML (client-side vanilla JS) so
that adjusting a LUT / opacity / window never triggers a Streamlit rerun
(which would otherwise re-stream the whole inline payload).
"""

from __future__ import annotations

import base64
import json
import math
import zlib
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np

from organoid_analysis.microscopy_io import validate_voxel_spacing_xyz

from .rendering_presets import DEFAULT_LUTS, PRESETS, Lut

RenderMode = Literal["volume", "mip", "surface"]

MAX_CHANNELS = 4
_ALLOWED_DTYPES = (np.uint8, np.uint16, np.float32)

_VENDOR = Path(__file__).parent / "static" / "vtk.js"


@dataclass
class ChannelConfig:
    """Per-channel rendering settings.

    ``opacity`` is the peak opacity (0..1). ``lut`` is a PRESETS key or an
    explicit ``Lut`` (list of ``(pos, r, g, b)`` control points).
    """

    lut: str | Lut = "cyan"
    opacity: float = 0.9
    name: str = ""
    window: tuple[float, float] | None = None  # (low, high) raw-intensity window

    def resolve_lut(self) -> Lut:
        if isinstance(self.lut, str):
            if self.lut not in PRESETS:
                raise ValueError(
                    f"Unknown LUT preset {self.lut!r}; available: {sorted(PRESETS)}"
                )
            return PRESETS[self.lut]
        return self.lut

    def validate(self) -> None:
        if not (0.0 <= self.opacity <= 1.0):
            raise ValueError(f"channel opacity must be in [0, 1], got {self.opacity}")
        if self.window is not None and not (
            len(self.window) == 2 and self.window[1] > self.window[0]
        ):
            raise ValueError(f"window must be (low, high) with high > low: {self.window}")


@dataclass
class ViewerSpec:
    """Fully packed, ready-to-serve viewer payload."""

    dims: tuple[int, int, int]             # (X, Y, Z)
    spacing: tuple[float, float, float]    # (x, y, z) microns/voxel
    n_channels: int
    channels: list[dict]
    data_b64: str                          # base64 of (optionally zlib) interleaved float32
    compressed: bool
    mode: RenderMode
    surface: dict | None = None            # packed surface mesh (when mode == "surface")
    mask: dict | None = None               # optional mask overlay (labels + surface + palette)
    height_px: int = 700


validate_spacing = validate_voxel_spacing_xyz


def _channel_stack(volumes: Sequence[np.ndarray]) -> np.ndarray:
    """Validate a list of per-channel volumes and stack into (C, Z, Y, X)."""
    if not volumes:
        raise ValueError("no volumes provided")
    if len(volumes) > MAX_CHANNELS:
        raise ValueError(f"up to {MAX_CHANNELS} channels are supported, got {len(volumes)}")
    stack: list[np.ndarray] = []
    first = None
    for i, v in enumerate(volumes):
        arr = np.asarray(v)
        if arr.ndim == 4:
            if arr.shape[-1] != 1:
                raise ValueError(f"channel {i}: 4D volume must be (Z, Y, X, 1), got {arr.shape}")
            arr = arr[..., 0]
        if arr.ndim not in (2, 3):
            raise ValueError(
                f"channel {i}: expected (Z, Y, X) or (Y, X), got shape {arr.shape}"
            )
        if arr.ndim == 2:
            arr = arr[np.newaxis, ...]
        if arr.dtype not in _ALLOWED_DTYPES:
            raise ValueError(
                f"channel {i}: dtype {arr.dtype} not supported; use uint8/uint16/float32"
            )
        if first is None:
            first = arr.shape
        elif arr.shape != first:
            raise ValueError(
                f"all channels must share shape {first}, channel {i} is {arr.shape}"
            )
        stack.append(arr)
    return np.stack([np.ascontiguousarray(a, dtype=np.float32) for a in stack], axis=0)


def _voxel_major_zyxc(stack_czyx: np.ndarray) -> np.ndarray:
    """Interleave ``(C,Z,Y,X)`` as VTK's x-fastest ``(Z,Y,X,C)`` wire array."""
    if stack_czyx.ndim != 4:
        raise ValueError(f"expected (C, Z, Y, X), got {stack_czyx.shape}")
    return np.ascontiguousarray(np.moveaxis(stack_czyx, 0, -1), dtype=np.float32)


def pack_channel_stack(
    volumes: Sequence[np.ndarray],
    channels: Sequence[ChannelConfig],
    *,
    percentile_window: bool = True,
    low_pct: float = 1.0,
    high_pct: float = 99.0,
    compress: bool = False,
) -> tuple[np.ndarray, list[dict], bytes]:
    """Convert channels to x-fastest ``(Z,Y,X,C)`` float32 wire data + metadata.

    uint16/uint8 are converted to float32 losslessly (uint16 max 65535 < 2^24,
    so float32 represents them exactly). Raw scalar values are preserved in the
    payload. Percentiles only choose the initial display window; browser-side
    window changes never overwrite the scalar array.

    NumPy ``(Z,Y,X)`` C-order already matches vtkImageData dimensions
    ``(X,Y,Z)`` because X is the fastest spatial index. The returned array has
    shape ``(Z,Y,X,C)``; only dimensions metadata is reversed later.
    """
    stack = _channel_stack(volumes)
    c, z, y, x = stack.shape
    if len(channels) != c:
        raise ValueError(f"Expected {c} channel configurations, got {len(channels)}")
    chans: list[dict] = []
    for i in range(c):
        cfg = channels[i]
        lut = cfg.resolve_lut()
        flat = stack[i]
        if not np.isfinite(flat).all():
            np.nan_to_num(flat, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
        actual_min, actual_max = float(flat.min()), float(flat.max())
        lo = hi = None
        if percentile_window:
            p_lo, p_hi = np.percentile(flat, (low_pct, high_pct))
            lo, hi = float(p_lo), float(p_hi)
        if cfg.window is not None:
            lo, hi = float(cfg.window[0]), float(cfg.window[1])
        if lo is None or hi is None or hi <= lo:
            lo, hi = actual_min, actual_max
        if hi <= lo:
            hi = lo + 1.0
        # Opacity estimation does not require every voxel. A deterministic
        # sample avoids another full-volume float array while preserving the
        # raw scalar payload and percentile display window.
        sample = flat.ravel()
        if sample.size > 250_000:
            sample = sample[::int(math.ceil(sample.size / 250_000))]
        preview = np.clip((sample - lo) / (hi - lo), 0.0, 1.0)
        opacity_pts = _opacity_control_points(preview, peak_opacity=cfg.opacity)
        chans.append(
            {
                "color": _lut_to_list(lut),
                "lut_key": cfg.lut if isinstance(cfg.lut, str) else None,
                "opacity": float(cfg.opacity),
                "opacity_pts": opacity_pts,
                "name": cfg.name,
                "window": [lo, hi],
                "actual_min": actual_min,
                "actual_max": actual_max,
                "suggested_lo": lo,
                "suggested_hi": hi,
            }
        )
    interleaved = _voxel_major_zyxc(stack)
    raw = interleaved.tobytes()
    if compress:
        raw = zlib.compress(raw, level=9)
    return interleaved, chans, raw


def _lut_to_list(lut: Lut) -> list[list[float]]:
    return [[float(p), float(r), float(g), float(b)] for p, r, g, b in lut]


def _opacity_control_points(flat: np.ndarray, peak_opacity: float) -> list[list[float]]:
    """Opacity TF control points scaled to this channel's value distribution.

    ``flat`` is the already-normalised (0..1) channel. Rather than assuming the
    data fills ``[0, 1]`` (the old TF only became visible above ~0.65, which
    hides signal in z-stack data that is heavily concentrated near zero), the
    ramp is anchored to percentiles so the visible structure maps to a high
    opacity: anything below the 70th percentile stays transparent (suppressing
    background), signal around the 70th..98th percentile becomes bright, and
    the top of the range reaches full opacity.
    """
    if not np.any(flat):
        return [[0.0, 0.0], [1.0, float(peak_opacity)]]
    p_lo = float(np.percentile(flat, 70))
    p_hi = float(np.percentile(flat, 98))
    if p_hi <= p_lo or p_hi <= 0.0:
        p_lo, p_hi = 0.25, 0.6
    p_lo = float(np.clip(p_lo, 0.01, 0.99))
    p_hi = float(np.clip(p_hi, p_lo + 0.01, 1.0))
    peak = float(np.clip(peak_opacity, 0.0, 1.0))
    return [
        [0.0, 0.0],
        [0.7 * p_lo, 0.05 * peak],
        [p_lo, 0.4 * peak],
        [p_hi, peak],
        [1.0, peak],
    ]


def build_viewer_payload(
    volumes: Sequence[np.ndarray],
    spacing_um: Sequence[float],
    channels: Sequence[ChannelConfig] | None = None,
    *,
    render_mode: RenderMode = "volume",
    segmentation_mask: np.ndarray | None = None,
    surface_level: float | None = None,
    mask_overlay: np.ndarray | None = None,
    mask_overlay_alpha: float = 0.5,
    height_px: int = 700,
    normalize: bool = True,
    compress: bool = False,
) -> ViewerSpec:
    """Build a ready-to-serve ``ViewerSpec`` from volumes / config.

    Default ``render_mode`` is ``"volume"`` (composite ray casting); surface is
    produced only when explicitly requested.
    """
    spacing = validate_spacing(spacing_um)
    if channels is None:
        channels = [
            ChannelConfig(lut=DEFAULT_LUTS[i % len(DEFAULT_LUTS)], name=f"Channel {i}")
            for i in range(len(volumes))
        ]
    for cfg in channels:
        cfg.validate()

    if render_mode == "surface":
        if segmentation_mask is None:
            raise ValueError("segmentation_mask is required for render_mode='surface'")
        mask = np.asarray(segmentation_mask)
        if mask.shape != np.asarray(volumes[0]).shape:
            raise ValueError(
                f"segmentation_mask shape {mask.shape} must match volume "
                f"{np.asarray(volumes[0]).shape}"
            )
        spacing_zyx = (spacing[2], spacing[1], spacing[0])  # (z, y, x) for skimage
        surface = build_surface_payload(mask, spacing_zyx, level=surface_level)
        return ViewerSpec(
            dims=(0, 0, 0),
            spacing=spacing,
            n_channels=0,
            channels=[],
            data_b64="",
            compressed=False,
            mode="surface",
            surface=surface,
            height_px=height_px,
        )

    interleaved, chans, raw = pack_channel_stack(
        volumes, channels, percentile_window=normalize, compress=compress
    )
    dims = (interleaved.shape[2], interleaved.shape[1], interleaved.shape[0])
    data_b64 = base64.b64encode(raw).decode("ascii")
    mode: RenderMode = render_mode if render_mode in ("volume", "mip") else "volume"
    mask = None
    if mask_overlay is not None:
        if np.asarray(mask_overlay).shape != np.asarray(volumes[0]).shape:
            raise ValueError(
                f"mask_overlay shape {np.asarray(mask_overlay).shape} must match "
                f"volume {np.asarray(volumes[0]).shape}"
            )
        mask = build_mask_overlay_payload(
            mask_overlay, spacing, alpha=mask_overlay_alpha
        )
    return ViewerSpec(
        dims=dims,
        spacing=spacing,
        n_channels=len(chans),
        channels=chans,
        data_b64=data_b64,
        compressed=compress,
        mode=mode,
        mask=mask,
        height_px=height_px,
    )


def build_surface_payload(
    mask: np.ndarray,
    spacing_zyx: Sequence[float],
    level: float | None = None,
    origin_zyx: Sequence[int] = (0, 0, 0),
) -> dict:
    """Run Marching Cubes and pack the mesh for vtk.js.

    ``spacing_zyx`` must already be in (z, y, x) order because
    ``skimage.measure.marching_cubes`` follows NumPy axis order. The returned
    vertices are transposed to (x, y, z) to match VTK's (X, Y, Z) convention.
    """
    from skimage import measure

    level = 0.5 if level is None else float(level)
    spacing = np.asarray(tuple(float(v) for v in spacing_zyx))
    origin = np.asarray(tuple(int(v) for v in origin_zyx))
    if spacing.shape != (3,) or origin.shape != (3,):
        raise ValueError("spacing_zyx and origin_zyx must each contain three values")
    binary = np.ascontiguousarray(np.asarray(mask) > 0, dtype=np.uint8)
    if binary.ndim != 3 or not np.any(binary):
        raise ValueError("Surface extraction requires a nonempty 3D mask")
    # Padding closes border-touching and all-foreground objects. Correct the
    # one-voxel shift and add the crop origin in physical coordinates.
    verts, faces, _n, _v = measure.marching_cubes(
        np.pad(binary, 1),
        level=level,
        spacing=tuple(spacing),
    )
    verts += origin * spacing - spacing
    verts_xyz = np.ascontiguousarray(verts[:, ::-1], dtype=np.float32)
    # vtkCellArray uses the legacy format: flat [npts, i0, i1, ..., npts, ...].
    # marching_cubes returns (n, 3) triangle indices -> insert a leading 3 per cell.
    n_cells = faces.shape[0]
    legacy = np.empty(n_cells * 4, dtype=np.uint32)
    legacy[0::4] = 3
    legacy[1::4] = faces[:, 0]
    legacy[2::4] = faces[:, 1]
    legacy[3::4] = faces[:, 2]
    return {
        "verts_b64": base64.b64encode(verts_xyz.tobytes()).decode("ascii"),
        "n_points": int(verts_xyz.shape[0]),
        "faces_b64": base64.b64encode(legacy.tobytes()).decode("ascii"),
        "n_faces": n_cells,
        "level": level,
    }


_MASK_COLORS = (
    (0.95, 0.20, 0.35),
    (0.25, 0.55, 1.00),
    (1.00, 0.80, 0.20),
    (0.25, 0.85, 0.45),
    (0.90, 0.40, 0.90),
    (0.30, 0.85, 0.85),
)

_MAX_PER_LABEL_SURFACES = 40


def _mask_palette(label_values: np.ndarray) -> list[list[float]]:
    """Deterministic RGB per label, bounded by the shared palette length."""
    return [
        list(_MASK_COLORS[(int(v) - 1) % len(_MASK_COLORS)])
        for v in label_values
    ]


def _surface_vertex_colors(surface: dict, labels: np.ndarray, spacing_zyx: Sequence[float]) -> str:
    """Per-vertex RGB (uint8, base64) for a single merged multi-label surface.

    Used only when there are too many objects for one actor per label (see
    ``_MAX_PER_LABEL_SURFACES``): every vertex is mapped back to a voxel and
    coloured by that voxel's nearest instance label, so the fallback single
    surface still shows individual objects instead of one flat colour.
    """
    from scipy import ndimage as ndi

    verts_xyz = np.frombuffer(base64.b64decode(surface["verts_b64"]), dtype=np.float32).reshape(-1, 3)
    spacing = np.asarray(tuple(spacing_zyx), dtype=np.float32)
    voxel_zyx = np.round(verts_xyz[:, ::-1] / spacing).astype(np.int64)
    voxel_zyx = np.clip(voxel_zyx, 0, np.array(labels.shape) - 1)
    # Marching-cubes vertices sit on the foreground/background boundary, so a
    # vertex can round to a background voxel; fill every background voxel
    # with its nearest foreground label first so every vertex still resolves
    # to a real object rather than "no label".
    nearest_index = ndi.distance_transform_edt(labels == 0, return_distances=False, return_indices=True)
    filled = np.where(labels > 0, labels, labels[tuple(nearest_index)])
    vertex_labels = filled[voxel_zyx[:, 0], voxel_zyx[:, 1], voxel_zyx[:, 2]].astype(np.int64)
    palette = np.asarray(_MASK_COLORS, dtype=np.float32)
    rgb = palette[(vertex_labels - 1) % len(_MASK_COLORS)]
    rgb_u8 = np.clip(np.rint(rgb * 255.0), 0, 255).astype(np.uint8)
    return base64.b64encode(np.ascontiguousarray(rgb_u8).tobytes()).decode("ascii")


def build_mask_overlay_payload(
    mask: np.ndarray,
    spacing_um: Sequence[float],
    alpha: float = 0.5,
) -> dict | None:
    """Pack a labeled ``(Z, Y, X)`` mask as an overlay payload.

    Returns ``None`` when ``mask`` is empty (no labels). The payload carries the
    raw uint32 labels (zlib-compressed) so the browser can draw 2D contours, a
    marching-cubes surface for the 3D overlay, and a per-label color palette.

    When the label count is at most ``_MAX_PER_LABEL_SURFACES``, individual
    marching-cubes surfaces are also returned so the 3D viewer can colour each
    object distinctly instead of rendering a single monochrome surface for all
    labels.
    """

    labels = np.asarray(mask)
    if labels.ndim != 3:
        raise ValueError(f"mask must be (Z, Y, X), got shape {labels.shape}")
    if labels.dtype != np.bool_ and not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("mask must contain integer instance labels")
    if np.issubdtype(labels.dtype, np.signedinteger) and np.any(labels < 0):
        raise ValueError("mask labels must be nonnegative")
    if labels.size and int(labels.max()) > np.iinfo(np.uint32).max:
        raise ValueError("mask labels exceed the supported uint32 range")
    if not 0.0 <= float(alpha) <= 1.0:
        raise ValueError("mask overlay alpha must be in [0, 1]")
    label_values = np.unique(labels)
    label_values = label_values[label_values != 0]
    if label_values.size == 0:
        return None
    u32 = np.ascontiguousarray(labels, dtype=np.uint32)
    sx, sy, sz = validate_spacing(spacing_um)
    surface = build_surface_payload(
        (u32 > 0).astype(np.uint8),
        spacing_zyx=(sz, sy, sx),
        level=0.5,
    )
    # Per-label coloured surfaces (limit MC calls for large label sets).
    label_surfaces: list[dict] = []
    if label_values.size <= _MAX_PER_LABEL_SURFACES:
        from scipy import ndimage as ndi

        if int(label_values[-1]) <= max(100_000, 4 * len(label_values)):
            local_labels = u32
            boxes = ndi.find_objects(local_labels)
            entries = [(int(lv), int(lv), boxes[int(lv) - 1]) for lv in label_values]
        else:
            # Avoid ndi.find_objects allocating up to a sparse, untrusted max ID.
            local_labels = (np.searchsorted(label_values, u32) + 1).astype(np.uint16)
            local_labels[u32 == 0] = 0
            boxes = ndi.find_objects(local_labels)
            entries = [(int(lv), index, boxes[index - 1])
                       for index, lv in enumerate(label_values, start=1)]
        for original_id, local_id, bbox in entries:
            if bbox is None:
                continue
            binary = np.ascontiguousarray(local_labels[bbox] == local_id, dtype=np.uint8)
            origin = tuple(item.start for item in bbox)
            ls = build_surface_payload(binary, spacing_zyx=(sz, sy, sx), level=0.5,
                                       origin_zyx=origin)
            ls["color"] = list(
                _MASK_COLORS[(original_id - 1) % len(_MASK_COLORS)]
            )
            label_surfaces.append(ls)
    else:
        # Too many objects for one actor each (a real field can have hundreds
        # of nuclei): colour the single merged surface per-vertex by nearest
        # label instead of falling back to one flat colour for everything,
        # which otherwise renders as an unreadable solid blob.
        surface["vertex_colors_b64"] = _surface_vertex_colors(surface, u32, (sz, sy, sx))
    raw = u32.tobytes()
    return {
        "shape": list(labels.shape),  # (Z, Y, X)
        "labels_b64": base64.b64encode(raw).decode("ascii"),
        "labels_compressed": False,
        "surface": surface,
        "palette": _mask_palette(label_values),
        "label_surfaces": label_surfaces,
        "alpha": float(alpha),
    }


def spec_to_dict(spec: ViewerSpec) -> dict:
    return {
        "dims": list(spec.dims),
        "spacing": list(spec.spacing),
        "n_channels": spec.n_channels,
        "channels": spec.channels,
        "data_b64": spec.data_b64,
        "compressed": spec.compressed,
        "mode": spec.mode,
        "surface": spec.surface,
        "mask": spec.mask,
        "height_px": spec.height_px,
        "luts": {name: _lut_to_list(lut) for name, lut in PRESETS.items()},
    }


def render_viewer_html(
    spec: ViewerSpec,
    *,
    inline_vtkjs: bool = True,
) -> str:
    """Produce a complete standalone HTML page string for ``spec``.

    With ``inline_vtkjs=True`` (default) the vendored vtk.js UMD bundle is
    embedded directly into the page so there is no external request.
    """
    spec_json = json.dumps(spec_to_dict(spec), separators=(",", ":"))
    if inline_vtkjs:
        # Merge the vendored vtk.js bundle and the viewer script into a single
        # <script> block. Streamlit's srcdoc serializer drops the opening tag of
        # a *second* consecutive <script> element, which would leave the viewer
        # script unparsed. Concatenating avoids that entirely.
        combined = (_VENDOR.read_text(encoding="utf-8") + "\n" + viewer_js(spec_json))
        script = "<script>" + combined + "</script>"
    else:
        script = (
            '<script src="/visualization/volume_viewer/static/vtk.js"></script>\n'
            + "<script>"
            + viewer_js(spec_json)
            + "</script>"
        )
    return _html_template(script)


def viewer_js(spec_json: str) -> str:
    """The standalone viewer script with ``__SPEC_JSON__`` substituted."""
    return _VIEWER_JS.replace("__SPEC_JSON__", spec_json)


def _html_template(script: str) -> str:
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>vtk.js Volume Viewer</title>
<style>
  html,body{{margin:0;padding:0;height:100%;background:#000;color:#eee;
        font-family:-apple-system,'Segoe UI',Roboto,sans-serif;font-size:13px;overflow:hidden}}
  #workspace{{position:absolute;inset:0;display:grid;grid-template-columns:minmax(280px,42%) 1fr;gap:1px;background:#223}}
  #slice-pane,#viewer-pane{{position:relative;min-width:0;min-height:0;background:#000;overflow:hidden}}
  #slice-pane{{display:flex;flex-direction:column;align-items:stretch}}
  .pane-title{{height:34px;display:flex;align-items:center;justify-content:space-between;padding:0 12px;
        color:#bfe7ff;background:#10141c;border-bottom:1px solid #293342;font-weight:600;box-sizing:border-box}}
  #slice-stage{{position:relative;flex:1;min-height:0;display:flex;align-items:center;justify-content:center;padding:12px;box-sizing:border-box}}
  #slice-canvas{{display:block;max-width:100%;max-height:100%;width:auto;height:auto;background:#000}}
  #slice-footer{{height:46px;display:flex;align-items:center;gap:9px;padding:0 12px;background:#10141c;
        border-top:1px solid #293342;box-sizing:border-box}}
  #slice-slider{{flex:1}}
  #viewer{{position:absolute;inset:0}}
  #controls{{position:absolute;top:8px;left:8px;z-index:20;background:rgba(24,24,30,.9);
        border:1px solid #333;border-radius:6px;padding:8px 10px;max-width:300px}}
  #controls h4{{margin:2px 0 6px;font-size:12px}}
  .row{{margin:5px 0;display:flex;align-items:center;gap:6px}}
  .row label{{width:auto;min-width:40px}}
  input[type=range]{{width:110px;vertical-align:middle}}
  button{{margin:2px;background:#223;border:1px solid #446;color:#def;border-radius:4px;padding:3px 8px;cursor:pointer}}
  button.active{{background:#46c}}
  #loading{{position:absolute;inset:0;display:flex;align-items:center;justify-content:center;
        z-index:30;background:#000}}
  #msg{{position:absolute;bottom:8px;left:8px;z-index:20;color:#ff6;font-size:12px}}
  #spacing{{position:absolute;bottom:8px;right:8px;z-index:20;color:#9ad;font-size:11px}}
  select{{background:#182033;border:1px solid #456;color:#def;border-radius:4px;padding:2px 5px}}
  @media(max-width:760px){{#workspace{{grid-template-columns:1fr;grid-template-rows:45% 55%}}}}
</style>
</head>
<body>
<div id="workspace">
  <section id="slice-pane">
    <div class="pane-title"><span>2D Z slice</span><span id="slice-index"></span></div>
    <div id="slice-stage"><canvas id="slice-canvas"></canvas></div>
    <div id="slice-footer"><label for="slice-slider">Z</label><input id="slice-slider" type="range"/><span id="slice-count"></span></div>
  </section>
  <section id="viewer-pane">
    <div id="viewer"></div>
    <div id="controls"></div>
    <div id="msg"></div>
    <div id="spacing"></div>
  </section>
</div>
<div id="loading"><div style="color:#9cf">Initializing volume renderer…</div></div>
{script}
</body>
</html>
"""


_VIEWER_JS = r"""
(function(){
  const SPEC = __SPEC_JSON__;
  const viewerDiv = document.getElementById('viewer');
  const loadingDiv = document.getElementById('loading');
  const msgDiv = document.getElementById('msg');
  const controlsDiv = document.getElementById('controls');
  const spacingDiv = document.getElementById('spacing');
  const sliceCanvas = document.getElementById('slice-canvas');
  const sliceSlider = document.getElementById('slice-slider');
  const sliceIndex = document.getElementById('slice-index');
  const sliceCount = document.getElementById('slice-count');

  function report(err){ if(msgDiv) msgDiv.textContent='Error: '+err;
    if(window.console) console.error('[viewer]', err); }
  function log(m){ if(window.console) console.log('[viewer]', m); }

  if(!(window.vtk || globalThis.vtk)){
    loadingDiv.style.display='none'; report('vtk.js failed to load'); return; }
  const vtk = window.vtk || globalThis.vtk;
  const RC = vtk.Rendering.Core, DD = vtk.Rendering.OpenGL,
        CD = vtk.Common && vtk.Common.DataModel,
        CC = vtk.Common && vtk.Common.Core,
        Misc = vtk.Rendering && vtk.Rendering.Misc;
  // ── compat resolution across vtk.js versions ────────────────────────────
  const GenericRenderWindow = Misc.vtkGenericRenderWindow || RC.vtkGenericRenderWindow;
  const vtkImageData  = CD.vtkImageData;
  const vtkDataArray  = CC.vtkDataArray || (CD && CD.vtkDataArray);
  const vtkBase64     = CC.vtkBase64;
  const vtkVolumeMapper   = RC.vtkVolumeMapper;
  const vtkVolumeProperty = RC.vtkVolumeProperty;
  const vtkVolume         = RC.vtkVolume;
  const vtkColorTransferFunction = RC.vtkColorTransferFunction;
  const vtkPiecewiseFunction     = CD.vtkPiecewiseFunction;
  const vtkPolyData       = CD.vtkPolyData;
  const vtkPoints         = CC.vtkPoints;
  const vtkCellArray      = CC.vtkCellArray;
  const vtkPolyDataMapper = RC.vtkMapper || DD.vtkPolyDataMapper || RC.vtkActor;
  const vtkActor          = RC.vtkActor;
  const need = {GenericRenderWindow, vtkImageData, vtkDataArray, vtkVolumeMapper,
    vtkVolumeProperty, vtkVolume, vtkColorTransferFunction, vtkPiecewiseFunction};
  for(const k in need){ if(!need[k]){ loadingDiv.style.display='none';
    report('missing vtk.namespace: '+k); return; } }

  const gw = GenericRenderWindow.newInstance();
  gw.setContainer(viewerDiv);
  const renderer = gw.getRenderer();
  const renderWindow = gw.getRenderWindow();
  const interactor = gw.getInteractor();
  renderer.setBackground(0,0,0);
  spacingDiv.textContent = 'spacing xyz = ' + SPEC.spacing.join(', ') + ' µm';
  function sizeToContainer(){
    try{
      const w = viewerDiv.clientWidth || viewerDiv.offsetWidth;
      const h = viewerDiv.clientHeight || viewerDiv.offsetHeight;
      if(w && h && typeof gw.resize === 'function'){
        gw.resize(w, h);
        // GenericRenderWindow normally applies DPR itself. Keep an explicit
        // fallback for versions/backends that only allocate CSS-pixel buffers.
        const canvas=viewerDiv.querySelector('canvas'), dpr=window.devicePixelRatio||1;
        const expectedW=Math.round(w*dpr), expectedH=Math.round(h*dpr);
        if(canvas && (Math.abs(canvas.width-expectedW)>1 || Math.abs(canvas.height-expectedH)>1)){
          const api=gw.getApiSpecificRenderWindow();
          if(api && typeof api.setSize==='function') api.setSize(expectedW,expectedH);
        }
      }
    }catch(e){}
  }
  if(typeof gw.onResize === 'function') gw.onResize(sizeToContainer);
  window.addEventListener('resize', sizeToContainer);
  sizeToContainer();

  function decodeB64(b64){
    if(vtkBase64 && vtkBase64.toArrayBuffer) return vtkBase64.toArrayBuffer(b64);
    const bin = atob(b64); const u8 = new Uint8Array(bin.length);
    for(let i=0;i<bin.length;i++) u8[i]=bin.charCodeAt(i);
    return u8.buffer;
  }

  let decodedData = null;
  function getData(){
    if(decodedData) return decodedData;
    const buf = decodeB64(SPEC.data_b64);
    decodedData = new Float32Array(buf.buffer || buf, buf.byteOffset||0, buf.byteLength/4);
    const [dx,dy,dz] = SPEC.dims, expected=dx*dy*dz*SPEC.n_channels;
    if(decodedData.length !== expected){
      throw new Error('payload length mismatch: got '+decodedData.length+', expected '+expected);
    }
    return decodedData;
  }

  function makeCTF(color, win){
    const ctf = vtkColorTransferFunction.newInstance();
    const lo=win?win.lo:0, hi=win?win.hi:1, span=Math.max(hi-lo, 0.001);
    for(const p of color) ctf.addRGBPoint(lo+p[0]*span, p[1], p[2], p[3]);
    return ctf;
  }
  // Fixed fallback opacity ramp, used only when the caller has not supplied
  // percentile-anchored `opacity_pts` (see Python-side
  // `_opacity_control_points` in this module, which is the normal path and
  // is documented there). This heuristic, display-only default is not tied
  // to any channel's actual value distribution; it does not affect any
  // quantitative measurement, only what intensity range renders as visible
  // "signal" in the 3D preview.
  function makeOp(peak){
    const pf = vtkPiecewiseFunction.newInstance();
    pf.addPoint(0.0, 0.0);
    pf.addPoint(0.30, 0.02*peak);
    pf.addPoint(0.65, 0.45*peak);
    pf.addPoint(1.0, peak);
    return pf;
  }
  function makeOpFromPts(pts){
    const pf = vtkPiecewiseFunction.newInstance();
    for(const p of (pts||[])) pf.addPoint(p[0], p[1]);
    return pf;
  }
  function scaleOp(c, peak, win){
    // keep the percentile-shaped opacity ramp, re-notch overall gain to `peak`
    const pts = (c.opacity_pts || []) ;
    const pf = vtkPiecewiseFunction.newInstance();
    const orig = c.opacity || peak;
    const g = orig>0 ? (peak/orig) : peak;
    const lo=win?win.lo:0, hi=win?win.hi:1, span=Math.max(hi-lo, 0.001);
    for(const p of pts){ pf.addPoint(lo+p[0]*span, p[1]*g); }
    return pf;
  }

  const channelState = SPEC.channels.map((c,i)=>({
    lut: c.lut_key || Object.keys(SPEC.luts)[i % Object.keys(SPEC.luts).length],
    lo: c.suggested_lo, hi: c.suggested_hi,
    min: c.actual_min, max: c.actual_max>c.actual_min ? c.actual_max : c.actual_min+1,
    opacity: c.opacity, visible: true
  }));
  const QUALITY_SAMPLE_DISTANCE={smooth:0.5,balanced:0.2,sharp:0.1};
  let qualityMode='balanced';
  let activeChannel=0, currentSlice=Math.floor(SPEC.dims[2]/2);

  function lutColor(key, t){
    const pts=SPEC.luts[key] || SPEC.channels[activeChannel].color;
    t=Math.max(0,Math.min(1,t));
    for(let i=1;i<pts.length;i++){
      if(t<=pts[i][0]){
        const a=pts[i-1], b=pts[i], f=(t-a[0])/Math.max(b[0]-a[0],1e-6);
        return [a[1]+f*(b[1]-a[1]),a[2]+f*(b[2]-a[2]),a[3]+f*(b[3]-a[3])];
      }
    }
    const p=pts[pts.length-1]; return [p[1],p[2],p[3]];
  }

  function renderSlice(){
    if(!sliceCanvas || SPEC.mode==='surface') return;
    const data=getData(), [dx,dy,dz]=SPEC.dims, nc=SPEC.n_channels;
    currentSlice=Math.max(0,Math.min(dz-1,currentSlice));
    sliceCanvas.width=dx; sliceCanvas.height=dy;
    const ctx=sliceCanvas.getContext('2d'), image=ctx.createImageData(dx,dy);
    const s=channelState[activeChannel], span=Math.max(s.hi-s.lo,0.001);
    for(let y=0;y<dy;y++) for(let x=0;x<dx;x++){
      const v=data[((currentSlice*dy+y)*dx+x)*nc+activeChannel];
      const rgb=lutColor(s.lut,(v-s.lo)/span), q=(y*dx+x)*4;
      image.data[q]=Math.round(rgb[0]*255); image.data[q+1]=Math.round(rgb[1]*255);
      image.data[q+2]=Math.round(rgb[2]*255); image.data[q+3]=255;
    }
    ctx.putImageData(image,0,0);
    if(SPEC.mask){
      getMaskLabels();
      drawMaskContours(image, dx, dy);
      ctx.putImageData(image,0,0);
    }
    sliceIndex.textContent='z = '+currentSlice;
    sliceCount.textContent=(currentSlice+1)+' / '+dz;
  }

  let maskLabels=null;
  function getMaskLabels(){
    if(!SPEC.mask || maskLabels) return maskLabels;
    const buf=decodeB64(SPEC.mask.labels_b64);
    maskLabels=new Uint32Array(buf.buffer||buf, buf.byteOffset||0, buf.byteLength/4);
    return maskLabels;
  }
  function drawMaskContours(image, dx, dy){
    if(!maskLabels) return;
    const [mz,my,mx]=SPEC.mask.shape;
    const pal=SPEC.mask.palette, alpha=SPEC.mask.alpha, a=(alpha===undefined?0.5:alpha);
    for(let y=0;y<dy;y++) for(let x=0;x<dx;x++){
      const idx=(currentSlice*my+y)*mx+x;
      const lab=maskLabels[idx];
      if(lab===0) continue;
      const isBoundary = x===0||y===0||x===mx-1||y===my-1||
        maskLabels[idx-1]!==lab || maskLabels[idx+1]!==lab ||
        maskLabels[idx-mx]!==lab || maskLabels[idx+mx]!==lab;
      if(!isBoundary) continue;
      const col=pal[(lab-1)%pal.length];
      const q=(y*dx+x)*4;
      image.data[q]  =Math.round(col[0]*255);
      image.data[q+1]=Math.round(col[1]*255);
      image.data[q+2]=Math.round(col[2]*255);
      image.data[q+3]=Math.round(255*a+image.data[q+3]*(1-a));
    }
  }

  function buildVolume(){
    if(SPEC.compressed){ console.warn('[viewer] zlib compressed path not supported in sync mode'); }
    const arr = getData();
    const [dx,dy,dz] = SPEC.dims, nc = SPEC.n_channels;
    const colors = SPEC.channels;
    const da = vtkDataArray.newInstance({ values: arr, numberOfComponents: nc });
    const imageData = vtkImageData.newInstance();
    imageData.setDimensions(dx,dy,dz);
    imageData.setSpacing(SPEC.spacing[0],SPEC.spacing[1],SPEC.spacing[2]);
    imageData.getPointData().setScalars(da);
    const mapper = vtkVolumeMapper.newInstance();
    mapper.setInputData(imageData);
    mapper.setAutoAdjustSampleDistances(false);
    mapper.setSampleDistance(QUALITY_SAMPLE_DISTANCE[qualityMode]);
    if(typeof mapper.setMaximumSamplesPerRay==='function') mapper.setMaximumSamplesPerRay(8192);
    // vtk.js WebGL always jitters ray starts internally. Some vtk.js variants
    // expose an explicit toggle; enable it where available.
    if(typeof mapper.setJittering==='function') mapper.setJittering(true);
    if(SPEC.mode === 'mip') mapper.setBlendModeToMaximumIntensity();
    else mapper.setBlendModeToComposite();
    const prop = vtkVolumeProperty.newInstance();
    prop.setIndependentComponents(true);
    prop.setInterpolationTypeToLinear();
    prop.setShade(false);
    const unitDistance=Math.min(...SPEC.spacing);
    for(let i=0;i<nc;i++){
      const c = colors[i];
      const cs=channelState[i], color=SPEC.luts[cs.lut] || c.color;
      prop.setRGBTransferFunction(i,makeCTF(color,cs));
      prop.setScalarOpacity(i,c.opacity_pts?scaleOp(c,cs.opacity,cs):makeOp(c.opacity));
      prop.setScalarOpacityUnitDistance(i,unitDistance);
      if(typeof prop.setComponentWeight==='function') prop.setComponentWeight(i,1);
    }
    const volume=vtkVolume.newInstance();
    volume.setMapper(mapper);volume.setProperty(prop);renderer.addVolume(volume);
    return {mode:SPEC.mode,implementation:'multi-component-single-volume',
      volumes:[volume],mappers:[mapper],props:[prop],volume,mapper,prop,imageData,
      n_channels:nc,quality:qualityMode,jittering:
        typeof mapper.setJittering==='function'?'explicit':'vtkjs-webgl-built-in'};
  }

  function applyChannelState(i){
    if(!state || !state.prop) return;
    const c=SPEC.channels[i], cs=channelState[i], prop=state.prop;
    prop.setRGBTransferFunction(i,makeCTF(SPEC.luts[cs.lut]||c.color,cs));
    prop.setScalarOpacity(i,scaleOp(c,cs.visible?cs.opacity:0,cs));
    if(typeof prop.setComponentWeight==='function') prop.setComponentWeight(i,cs.visible?1:0);
    renderWindow.render(); renderSlice();
  }

  function applyQuality(mode){
    if(!QUALITY_SAMPLE_DISTANCE[mode]) return;
    qualityMode=mode;
    if(state && state.mapper){
      state.mapper.setAutoAdjustSampleDistances(false);
      state.mapper.setSampleDistance(QUALITY_SAMPLE_DISTANCE[mode]);
      state.quality=mode; renderWindow.render();
    }
  }

  function buildSurface(){
    return buildSurfaceFrom(SPEC.surface);
  }

  let maskActor = null;
  function buildMaskOverlay(){
    if(!SPEC.mask || !SPEC.mask.surface) return;
    // Per-label coloured surfaces (best case: individual organoids distinguishable).
    const ls = SPEC.mask.label_surfaces;
    if(ls && ls.length > 0){
      ls.forEach(function(entry){
        const actor = buildSurfaceFrom(entry);
        actor.getProperty().setColor(entry.color[0], entry.color[1], entry.color[2]);
        actor.getProperty().setOpacity(SPEC.mask.alpha === undefined ? 0.5 : SPEC.mask.alpha);
        renderer.addActor(actor);
      });
      return;
    }
    // Fallback: single merged surface. Too many objects for one actor each
    // (see _MAX_PER_LABEL_SURFACES) -- coloured per-vertex by nearest label
    // when available, so a crowded field still shows individual objects
    // instead of one flat colour.
    const vcolors = SPEC.mask.surface.vertex_colors_b64;
    const actor = buildSurfaceFrom(SPEC.mask.surface, vcolors);
    if(!vcolors){
      const pal = SPEC.mask.palette || [[0.95,0.2,0.35]];
      const col = pal[0];
      actor.getProperty().setColor(col[0], col[1], col[2]);
    }
    actor.getProperty().setOpacity(SPEC.mask.alpha === undefined ? 0.5 : SPEC.mask.alpha);
    renderer.addActor(actor);
    maskActor = actor;
  }
  function buildSurfaceFrom(s, vertexColorsB64){
    const vertsBuf = decodeB64(s.verts_b64), facesBuf = decodeB64(s.faces_b64);
    const pts = vtkPoints.newInstance();
    pts.setData(new Float32Array(vertsBuf), 3);
    const polys = vtkCellArray.newInstance();
    polys.setData(new Uint32Array(facesBuf));
    const poly = vtkPolyData.newInstance();
    poly.setPoints(pts); poly.setPolys(polys);
    const colorsB64 = vertexColorsB64 || s.vertex_colors_b64;
    if(colorsB64){
      const colorBuf = decodeB64(colorsB64);
      const colors = vtkDataArray.newInstance({numberOfComponents:3, values:new Uint8Array(colorBuf)});
      poly.getPointData().setScalars(colors);
    }
    const mapper = vtkPolyDataMapper.newInstance();
    mapper.setInputData(poly);
    if(colorsB64){
      mapper.setScalarModeToUsePointData();
      mapper.setColorModeToDirectScalars();
    }
    const actor = vtkActor.newInstance();
    actor.setMapper(mapper);
    return actor;
  }

  let state = null;
  function resetCameraForMode(){
    renderer.resetCamera();
    const cam=renderer.getActiveCamera();
    if(SPEC.mode==='mip'){
      const b=renderer.computeVisiblePropBounds();
      const cx=(b[0]+b[1])/2, cy=(b[2]+b[3])/2, cz=(b[4]+b[5])/2;
      const w=viewerDiv.clientWidth||1, h=viewerDiv.clientHeight||1, aspect=w/h;
      const sx=b[1]-b[0], sy=b[3]-b[2], sz=b[5]-b[4];
      cam.setParallelProjection(true); cam.setFocalPoint(cx,cy,cz);
      cam.setPosition(cx,cy,cz+Math.max(sx,sy,sz)*2.5);
      cam.setViewUp(0,1,0);
      cam.setParallelScale(0.51*Math.max(sy,sx/Math.max(aspect,0.01)));
      renderer.resetCameraClippingRange();
    }else{
      cam.setParallelProjection(false);
    }
  }
  function install(){
    if(SPEC.mode === 'surface'){
      const a = buildSurface(); renderer.addActor(a);
    } else {
      state = buildVolume();
      getMaskLabels();          // decode mask labels if present
      buildMaskOverlay();       // 3D surface overlay on top of the intensity volume
    }
    resetCameraForMode();
    renderWindow.render();
    try { interactor.initialize(); interactor.setDesiredUpdateRate(60); } catch(e){}
    loadingDiv.style.display='none';
    log('rendered mode=' + SPEC.mode + ' dims=' + JSON.stringify(SPEC.dims));
    if(spacingDiv) spacingDiv.style.display='block';
    if(SPEC.mask && !SPEC.channels.length){
      const m = document.getElementById('slice-slider');
      if(m){ m.min=0; m.max=Math.max(0,SPEC.dims[2]-1); }
    }
    try { window.__viewer_state = state; } catch(e){}
    renderSlice();
  }

  function buildControls(){
    const h = document.createElement('h4'); h.textContent = 'Render controls';
    controlsDiv.appendChild(h);
    if(!SPEC.channels.length){
      const row=document.createElement('div'); row.className='row'; row.textContent='Surface';
      controlsDiv.appendChild(row); return;
    }
    const shared=document.createElement('div'); shared.className='row';
    const channelSelect=document.createElement('select'); channelSelect.id='channel-select';
    SPEC.channels.forEach((c,i)=>{ const o=document.createElement('option'); o.value=i;
      o.textContent=c.name||('C'+i); channelSelect.appendChild(o); });
    const lutSelect=document.createElement('select'); lutSelect.id='lut-select';
    ['cyan','magenta','yellow','green','gray'].forEach(k=>{ if(!SPEC.luts[k]) return;
      const o=document.createElement('option'); o.value=k; o.textContent=k; lutSelect.appendChild(o); });
    lutSelect.value=channelState[0].lut;
    shared.appendChild(channelSelect); shared.appendChild(lutSelect); controlsDiv.appendChild(shared);
    const wrow=document.createElement('div'); wrow.className='row';
    const wlab=document.createElement('label'); wlab.textContent='Window';
    const low=document.createElement('input'); low.id='window-low'; low.type='range';
    const high=document.createElement('input'); high.id='window-high'; high.type='range';
    low.style.width=high.style.width='72px'; wrow.appendChild(wlab);wrow.appendChild(low);wrow.appendChild(high);controlsDiv.appendChild(wrow);
    function syncShared(){ const s=channelState[activeChannel], step=Math.max((s.max-s.min)/1000,1e-6);
      lutSelect.value=s.lut;low.min=s.min;low.max=s.max-step;low.step=step;low.value=s.lo;
      high.min=s.min+step;high.max=s.max;high.step=step;high.value=s.hi;renderSlice(); }
    channelSelect.onchange=()=>{activeChannel=+channelSelect.value;syncShared();};
    lutSelect.onchange=()=>{channelState[activeChannel].lut=lutSelect.value;applyChannelState(activeChannel);};
    low.oninput=()=>{const s=channelState[activeChannel],step=+low.step;s.lo=Math.min(+low.value,s.hi-step);low.value=s.lo;applyChannelState(activeChannel);};
    high.oninput=()=>{const s=channelState[activeChannel],step=+high.step;s.hi=Math.max(+high.value,s.lo+step);high.value=s.hi;applyChannelState(activeChannel);};
    syncShared();
    const qrow=document.createElement('div');qrow.className='row';
    const qlab=document.createElement('label');qlab.textContent='Quality';
    const quality=document.createElement('select');quality.id='quality-select';
    [['smooth','Smooth · 0.5 µm'],['balanced','Balanced · 0.2 µm'],['sharp','Sharp · 0.1 µm']].forEach(([v,t])=>{
      const o=document.createElement('option');o.value=v;o.textContent=t;quality.appendChild(o);});
    quality.value=qualityMode;quality.onchange=()=>applyQuality(quality.value);
    qrow.appendChild(qlab);qrow.appendChild(quality);controlsDiv.appendChild(qrow);
    const mrow = document.createElement('div'); mrow.className='row';
    mrow.innerHTML = '<label>Mode</label>';
    ['volume','mip','surface'].forEach(m=>{
      const b=document.createElement('button'); b.textContent=m;
      b.className = (m===SPEC.mode)?'active':'';
      if(m==='surface'){ // surface mode: no toggle unless surface payload bundled
        b.disabled = !SPEC.surface; }
      b.onclick = ()=>{
        if(m!=='volume' && m!=='mip') return;
        SPEC.mode=m;
        if(state && state.volumes) state.volumes.forEach(v=>renderer.removeVolume(v));
        state=buildVolume(); resetCameraForMode(); renderWindow.render();
        controlsDiv.querySelectorAll('button').forEach(x=>x.classList.remove('active'));
        b.classList.add('active');
      };
      mrow.appendChild(b);
    });
    controlsDiv.appendChild(mrow);
    // per-channel opacity + visibility
    SPEC.channels.forEach((c,i)=>{
      const vis=document.createElement('div'); vis.className='row';
      const cb=document.createElement('input'); cb.type='checkbox'; cb.checked=true;
      const lab=document.createElement('label'); lab.textContent=(c.name||('ch'+i))+' ';
      const sl=document.createElement('input'); sl.type='range'; sl.min=0; sl.max=1; sl.step=0.05;
      sl.value=c.opacity;
      cb.onchange=()=>{channelState[i].visible=cb.checked;applyChannelState(i);};
      sl.oninput=()=>{channelState[i].opacity=+sl.value;if(cb.checked)applyChannelState(i);};
      vis.appendChild(cb); vis.appendChild(lab); vis.appendChild(sl);
      controlsDiv.appendChild(vis);
    });
    sliceSlider.min=0;sliceSlider.max=Math.max(0,SPEC.dims[2]-1);sliceSlider.step=1;sliceSlider.value=currentSlice;
    sliceSlider.oninput=()=>{currentSlice=+sliceSlider.value;renderSlice();};
  }

  try {
    buildControls();
    install();
    window.__viewer_ui_state=()=>({activeChannel,currentSlice,qualityMode,
      channels:channelState.map(s=>({...s}))});
  } catch(err){
    loadingDiv.style.display='none';
    report('init failed: ' + ((err && err.stack) || err));
  }

  // ── on-page director diagnostic (read-back of the real GL framebuffer) ──
  function runDiagnostic(readPixels){
    try{
      const dom = viewerDiv.querySelector('canvas');
      const info = {hasCanvas: !!dom,
                    domSize: dom ? [dom.width, dom.height] : null,
                    errGl: null, non0: null, avgLum: null, center: null,
                    camPos: null, camFocal: null, camDist: null, camParallel: null,
                    camScale: null, bounds: null, dataLength: decodedData?decodedData.length:null,
                    dpr:window.devicePixelRatio||1, cssSize: null,
                    implementation:state?state.implementation:null, volumeCount:state?state.volumes.length:null,
                    componentCount:state&&state.imageData?state.imageData.getPointData().getScalars().getNumberOfComponents():null,
                    independentComponents:state&&state.prop?state.prop.getIndependentComponents():null,
                    autoAdjust:state&&state.mapper?state.mapper.getAutoAdjustSampleDistances():null,
                    sampleDistance:state&&state.mapper?state.mapper.getSampleDistance():null,
                    shade:state&&state.prop?state.prop.getShade():null};
      if(!dom){ log('[diag] no canvas'); return info; }
      const rect=dom.getBoundingClientRect(); info.cssSize=[rect.width,rect.height];
      const cam = renderer.getActiveCamera();
      try{ info.camPos = cam.getPosition().slice(); }catch(_){}
      try{ info.camFocal = cam.getFocalPoint().slice(); }catch(_){}
      try{ info.camDist = cam.getDistance(); }catch(_){}
      try{ info.camParallel = cam.getParallelProjection(); }catch(_){}
      try{ info.camScale = cam.getParallelScale(); }catch(_){}
      try{
        const b = renderer.computeVisiblePropBounds(); info.bounds = Array.from(b); }catch(_){}
      // read back from the OpenGL render window's drawing buffer
      const ogl = gw.getApiSpecificRenderWindow();
      const gl = (ogl && ogl.getContext) ? ogl.getContext() : null;
      if(!gl){ log('[diag] no GL context via OpenGL render window'); }
      let p = null;
      if(gl && readPixels){
        info.errGl = gl.getError();
        gl.readPixels(0,0,1,1,gl.RGBA,gl.UNSIGNED_BYTE,(p=new Uint8Array(4)));
        info.corner = Array.from(p);
        const W=dom.width||1024, H=dom.height||1024, step=64;
        let non0=0,total=0,sum=0;
        const tmp=new Uint8Array(4);
        for(let y=0;y<H;y+=step) for(let x=0;x<W;x+=step){
          gl.readPixels(x,y,1,1,gl.RGBA,gl.UNSIGNED_BYTE,tmp); total++;
          if(tmp[3]>0) non0++;
          sum+=(tmp[0]+tmp[1]+tmp[2])/3;
        }
        info.non0=non0; info.total=total; info.avgLum=+(sum/total).toFixed(3);
        gl.readPixels(W>>1,H>>1,1,1,gl.RGBA,gl.UNSIGNED_BYTE,tmp); info.center=Array.from(tmp);
      }
      log('[diag] ' + JSON.stringify(info));
      return info;
    }catch(err){ log('[diag] error: ' + ((err&&err.stack)||err)); return null; }
  }
  window.__viewer_diag = runDiagnostic;
  setTimeout(()=>runDiagnostic(false), 1500);
})();

window.addEventListener('error', function(e){
  try{
    const m = document.getElementById('msg');
    if(m) m.textContent = 'uncaught: ' + ((e && e.message) || e);
  }catch(_){}
});
"""
