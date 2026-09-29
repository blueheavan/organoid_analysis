"""Full-resolution 2D review of Z planes with OpenSeadragon (display only).

The vtk.js volume viewer downsamples large stacks to fit a browser payload;
reviewing a segmentation needs the native pixels. This module cuts each
selected Z plane into a deep-zoom image pyramid (256 px tiles, levels halved
by area averaging) and, when labels are given, a matching pyramid of
instance-contour tiles, then serves them with the vendored OpenSeadragon
bundle (``static/openseadragon.min.js``, BSD-3-Clause, v6.1.1) inlined into a
standalone HTML page. Tiles are PNG/JPEG data URLs, so there is no tile
server, CDN request or file access from the iframe.

Coordinates follow the pipeline's voxel-centre convention: voxel index
``(z, y, x)`` has its centre at ``(z*sz, y*sy, x*sx)`` µm. In OpenSeadragon
image coordinates pixel column ``c`` spans ``[c, c+1)``, so its centre is at
``c + 0.5``; the readout reports the voxel under the pointer and the
physical position of that voxel's centre.

Display only. The page shows an 8-bit windowed rendering of the intensity
(one window for the whole stack, so planes stay comparable); it never feeds a
measurement and is not an annotation tool for segmentation validation (SG-3
needs blinded native-resolution 3D annotation under a written SOP).
"""
from __future__ import annotations

import base64
import io
import json
import math
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

import numpy as np
from scipy import ndimage as ndi  # type: ignore[import-untyped]

TILE_SIZE = 256
DEFAULT_BUDGET_BYTES = 40 * 2**20
MAX_LABEL_MARKERS = 400
ImageFormat = Literal["png", "jpeg"]

_VENDOR = Path(__file__).parent / "static" / "openseadragon.min.js"
# Tableau-like qualitative palette: distinct neighbouring labels, readable on grey.
_PALETTE = np.array([
    (255, 215, 0), (0, 200, 255), (255, 80, 80), (80, 230, 120), (220, 120, 255), (255, 150, 40),
    (120, 160, 255), (255, 110, 200), (170, 255, 60), (60, 255, 220), (255, 255, 140), (200, 170, 120),
], np.uint8)


@dataclass
class DeepZoomSpec:
    """Ready-to-serve deep-zoom payload for a subset of Z planes."""

    width: int
    height: int
    n_planes: int
    max_level: int
    tile_size: int
    spacing_zyx_um: tuple[float, float, float]
    window: tuple[float, float]
    image_format: str
    initial_plane: int
    planes: dict[int, dict] = field(default_factory=dict)   # z -> {"image": [...], "overlay": [...]|None, "ids": [...]}
    has_overlay: bool = False
    payload_bytes: int = 0
    budget_bytes: int = DEFAULT_BUDGET_BYTES
    height_px: int = 640
    title: str = ""

    @property
    def plane_indices(self) -> list[int]:
        return sorted(self.planes)


def max_level_for(width: int, height: int) -> int:
    """OpenSeadragon/DZI level of the full-resolution image (level 0 is 1 x 1)."""
    return max(0, math.ceil(math.log2(max(width, height, 1))))


def level_shape(width: int, height: int, level: int, max_level: int) -> tuple[int, int]:
    """``(rows, cols)`` of ``level``: each level halves the next one, rounding up."""
    scale = 2 ** (max_level - level)
    return math.ceil(height / scale), math.ceil(width / scale)


def _halve(image: np.ndarray, reducer: str) -> np.ndarray:
    rows, cols = image.shape[:2]
    pad = ((0, rows % 2), (0, cols % 2))
    if reducer == "mean":
        padded = np.pad(image, pad, mode="edge")
        return padded.reshape(padded.shape[0] // 2, 2, padded.shape[1] // 2, 2).mean(axis=(1, 3))
    padded = np.pad(image, pad, mode="constant")
    return padded.reshape(padded.shape[0] // 2, 2, padded.shape[1] // 2, 2).max(axis=(1, 3))


def pyramid(base: np.ndarray, max_level: int, reducer: Literal["mean", "max"]) -> list[np.ndarray]:
    """Images for levels ``0..max_level``; ``mean`` for intensity, ``max`` for contour labels."""
    levels = [base]
    for _ in range(max_level):
        levels.append(_halve(levels[-1], reducer))
    return levels[::-1]


def intensity_window(stack: np.ndarray, low_pct: float = 0.5, high_pct: float = 99.8,
                     max_samples: int = 4_000_000) -> tuple[float, float]:
    """One display window for the whole stack (planes stay comparable)."""
    values = np.asarray(stack, dtype=np.float64).ravel()
    values = values[:: max(1, values.size // max_samples)]
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 0.0, 1.0
    low, high = (float(v) for v in np.percentile(values, [low_pct, high_pct]))
    if not high > low:
        high = low + 1.0
    return low, high


def to_uint8(plane: np.ndarray, window: tuple[float, float]) -> np.ndarray:
    low, high = window
    scaled = (np.nan_to_num(np.asarray(plane, np.float64), nan=low) - low) / (high - low)
    return np.clip(np.rint(scaled * 255.0), 0, 255).astype(np.uint8)


def contour_labels(labels_plane: np.ndarray) -> np.ndarray:
    """Label ID on each object's inner boundary pixels, 0 elsewhere."""
    from skimage.segmentation import find_boundaries

    return np.where(find_boundaries(labels_plane, mode="inner"), labels_plane, 0)


def contour_rgba(contours: np.ndarray) -> np.ndarray:
    rgba = np.zeros((*contours.shape, 4), np.uint8)
    on = contours > 0
    rgba[on, :3] = _PALETTE[(contours[on].astype(np.int64) - 1) % len(_PALETTE)]
    rgba[on, 3] = 255
    return rgba


def _data_url(array: np.ndarray, image_format: str) -> str:
    from PIL import Image

    buffer = io.BytesIO()
    image = Image.fromarray(array)
    if image_format == "jpeg" and array.ndim == 2:
        image.save(buffer, format="JPEG", quality=92)
        mime = "image/jpeg"
    else:
        image.save(buffer, format="PNG", optimize=False, compress_level=6)
        mime = "image/png"
    return f"data:{mime};base64," + base64.b64encode(buffer.getvalue()).decode("ascii")


def _tiles(levels: list[np.ndarray], tile_size: int, encode: Callable[[np.ndarray], str], skip_empty: bool) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for image in levels:
        tiles: dict[str, str] = {}
        for row in range(math.ceil(image.shape[0] / tile_size)):
            for col in range(math.ceil(image.shape[1] / tile_size)):
                tile = image[row * tile_size:(row + 1) * tile_size, col * tile_size:(col + 1) * tile_size]
                if skip_empty and not tile.any():
                    continue
                tiles[f"{col}_{row}"] = encode(tile)
        out.append(tiles)
    return out


def _label_markers(labels_plane: np.ndarray) -> list[list[float]]:
    """``[id, x, y]`` per label in this plane, at its in-plane centroid in image coordinates."""
    ids = np.unique(labels_plane[labels_plane > 0])[:MAX_LABEL_MARKERS]
    if ids.size == 0:
        return []
    centres = ndi.center_of_mass(np.ones_like(labels_plane, np.uint8), labels_plane, ids)
    return [[int(i), round(float(c[1]) + 0.5, 2), round(float(c[0]) + 0.5, 2)] for i, c in zip(ids, centres, strict=True)]


def plane_payload(plane: np.ndarray, labels_plane: np.ndarray | None, window: tuple[float, float],
                  image_format: ImageFormat = "png", tile_size: int = TILE_SIZE) -> tuple[dict, int]:
    """Tiles for one plane and their size in bytes."""
    height, width = plane.shape
    top = max_level_for(width, height)
    image_levels = pyramid(np.asarray(plane, np.float64), top, "mean")
    image_tiles = _tiles([to_uint8(level, window) for level in image_levels], tile_size,
                         lambda tile: _data_url(tile, image_format), skip_empty=False)
    payload: dict = {"image": image_tiles, "overlay": None, "ids": []}
    if labels_plane is not None:
        contour_levels = pyramid(contour_labels(labels_plane), top, "max")
        payload["overlay"] = _tiles(contour_levels, tile_size,
                                    lambda tile: _data_url(contour_rgba(tile), "png"), skip_empty=True)
        payload["ids"] = _label_markers(labels_plane)
    size = sum(len(url) for level in image_tiles for url in level.values())
    if payload["overlay"] is not None:
        size += sum(len(url) for level in payload["overlay"] for url in level.values())
    return payload, size


def _plane_order(n_planes: int, centre: int) -> list[int]:
    order = [centre]
    for step in range(1, n_planes):
        order.extend(z for z in (centre + step, centre - step) if 0 <= z < n_planes)
    return order


def build_deepzoom_payload(image: np.ndarray, labels: np.ndarray | None = None,
                           spacing_zyx_um: tuple[float, float, float] = (1.0, 1.0, 1.0), *,
                           centre_plane: int | None = None, budget_bytes: int = DEFAULT_BUDGET_BYTES,
                           max_planes: int | None = None, window: tuple[float, float] | None = None,
                           image_format: ImageFormat = "png", tile_size: int = TILE_SIZE,
                           height_px: int = 640, title: str = "") -> DeepZoomSpec:
    """Tile the planes nearest ``centre_plane`` until ``budget_bytes`` is reached.

    The centre plane is always included; neighbours are added alternately
    above and below while the encoded payload stays within the budget.
    """
    stack = np.asarray(image)
    if stack.ndim == 2:
        stack = stack[np.newaxis]
    if stack.ndim != 3 or 0 in stack.shape:
        raise ValueError(f"image must be (Z, Y, X) or (Y, X), got shape {np.shape(image)}")
    if not (np.issubdtype(stack.dtype, np.integer) or np.issubdtype(stack.dtype, np.floating) or stack.dtype == bool):
        raise ValueError(f"image dtype {stack.dtype} is not numeric")
    label_stack = None
    if labels is not None:
        label_stack = np.asarray(labels)
        if label_stack.ndim == 2:
            label_stack = label_stack[np.newaxis]
        if label_stack.shape != stack.shape:
            raise ValueError(f"labels shape {np.shape(labels)} must match image shape {np.shape(image)}")
        if not np.issubdtype(label_stack.dtype, np.integer) or (label_stack < 0).any():
            raise ValueError("labels must be nonnegative integer instance labels")
    spacing = tuple(float(v) for v in spacing_zyx_um)
    if len(spacing) != 3 or not all(math.isfinite(v) and v > 0 for v in spacing):
        raise ValueError(f"spacing_zyx_um must be three positive finite values, got {spacing_zyx_um}")
    if image_format not in ("png", "jpeg"):
        raise ValueError("image_format must be 'png' or 'jpeg'")
    if tile_size < 16:
        raise ValueError("tile_size must be at least 16 px")
    n_planes, height, width = stack.shape
    centre = n_planes // 2 if centre_plane is None else int(centre_plane)
    if not 0 <= centre < n_planes:
        raise ValueError(f"centre_plane {centre_plane} outside 0..{n_planes - 1}")
    display_window = window if window is not None else intensity_window(stack)
    if not display_window[1] > display_window[0]:
        raise ValueError(f"window must be (low, high) with high > low, got {display_window}")
    spec = DeepZoomSpec(width=width, height=height, n_planes=n_planes, max_level=max_level_for(width, height),
                        tile_size=tile_size, spacing_zyx_um=spacing,  # type: ignore[arg-type]
                        window=(float(display_window[0]), float(display_window[1])), image_format=image_format,
                        initial_plane=centre, has_overlay=label_stack is not None, budget_bytes=int(budget_bytes),
                        height_px=int(height_px), title=title)
    limit = n_planes if max_planes is None else max(1, int(max_planes))
    for z in _plane_order(n_planes, centre)[:limit]:
        payload, size = plane_payload(stack[z], None if label_stack is None else label_stack[z],
                                      spec.window, image_format, tile_size)
        if spec.planes and spec.payload_bytes + size > budget_bytes:
            break
        spec.planes[z] = payload
        spec.payload_bytes += size
    return spec


def _empty_tile() -> str:
    return _data_url(np.zeros((1, 1, 4), np.uint8), "png")


def spec_to_dict(spec: DeepZoomSpec) -> dict:
    return {"width": spec.width, "height": spec.height, "n_planes": spec.n_planes, "max_level": spec.max_level,
            "tile_size": spec.tile_size, "spacing_zyx_um": list(spec.spacing_zyx_um), "window": list(spec.window),
            "initial_plane": spec.initial_plane, "has_overlay": spec.has_overlay, "title": spec.title,
            "plane_indices": spec.plane_indices, "planes": {str(z): spec.planes[z] for z in spec.plane_indices},
            "empty_tile": _empty_tile()}


def _vendor_js() -> str:
    lines = _VENDOR.read_text(encoding="utf-8").splitlines()
    return "\n".join(line for line in lines if not line.startswith("//# sourceMappingURL="))


def render_deepzoom_html(spec: DeepZoomSpec) -> str:
    """Standalone HTML page; vendored OpenSeadragon and the viewer share one ``<script>``.

    (Streamlit's srcdoc serializer drops the opening tag of a second
    consecutive ``<script>``; see the volume viewer.)
    """
    spec_json = json.dumps(spec_to_dict(spec), separators=(",", ":")).replace("</", "<\\/")
    script = "<script>" + _vendor_js() + "\n" + _VIEWER_JS.replace("__SPEC_JSON__", spec_json) + "</script>"
    return _TEMPLATE.replace("__HEIGHT__", str(spec.height_px)).replace("__SCRIPT__", script)


def st_deepzoom_viewer(image: np.ndarray, labels: np.ndarray | None = None,
                       spacing_zyx_um: tuple[float, float, float] = (1.0, 1.0, 1.0), *,
                       height: int = 640, use_st_iframe: bool = True, **kwargs: Any) -> DeepZoomSpec:
    """Render :func:`build_deepzoom_payload` in Streamlit and return the spec (for captions)."""
    import streamlit as st

    spec = build_deepzoom_payload(image, labels, spacing_zyx_um, height_px=height, **kwargs)
    html = render_deepzoom_html(spec)
    if use_st_iframe:
        try:
            st.iframe(html, height=height)
            return spec
        except Exception:
            pass
    from streamlit.components.v1 import html as _html

    _html(html, height=height, scrolling=False)
    return spec


_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>Deep-zoom plane review</title>
<style>
  :root{--bg:#0b0d12;--panel:rgba(16,20,28,.92);--line:#2a3444;--text:#e6edf5;--muted:#9fb0c4;--accent:#5aa9ff}
  html,body{margin:0;padding:0;height:100%;background:var(--bg);color:var(--text);
    font:13px/1.35 -apple-system,'Segoe UI',Roboto,sans-serif;overflow:hidden}
  #wrap{position:relative;height:__HEIGHT__px;max-height:100vh}
  #osd{position:absolute;inset:44px 0 0 0;background:#000}
  #bar{position:absolute;top:0;left:0;right:0;height:44px;display:flex;align-items:center;gap:10px;
    padding:0 10px;background:var(--panel);border-bottom:1px solid var(--line);box-sizing:border-box;overflow-x:auto}
  #bar button{background:#1c2432;border:1px solid #3a475c;color:var(--text);border-radius:4px;
    min-width:30px;height:28px;cursor:pointer;font-size:13px}
  #bar button:hover{border-color:var(--accent)}
  #bar label{color:var(--muted);white-space:nowrap;display:flex;align-items:center;gap:5px}
  #zslider{width:180px}
  #zlabel{min-width:150px;font-variant-numeric:tabular-nums;white-space:nowrap}
  #readout{position:absolute;left:8px;bottom:8px;z-index:5;background:var(--panel);border:1px solid var(--line);
    border-radius:4px;padding:4px 8px;font-variant-numeric:tabular-nums;pointer-events:none;white-space:pre}
  #scale{position:absolute;right:12px;bottom:10px;z-index:5;text-align:center;pointer-events:none;
    text-shadow:0 0 3px #000}
  #scalebar{height:4px;background:#fff;border:1px solid #000;margin:0 auto 2px}
  #note{position:absolute;right:8px;top:52px;z-index:5;background:var(--panel);border:1px solid var(--line);
    border-radius:4px;padding:3px 8px;color:var(--muted);max-width:45%}
  .idtag{color:#fff;font:600 11px/1 monospace;text-shadow:0 0 2px #000,0 0 2px #000;pointer-events:none;white-space:nowrap}
  #err{position:absolute;inset:44px 0 0 0;display:none;align-items:center;justify-content:center;color:#ff8080;padding:20px}
</style>
</head>
<body>
<div id="wrap">
  <div id="bar">
    <button id="zoomout" title="Zoom out">&minus;</button>
    <button id="zoomin" title="Zoom in">+</button>
    <button id="home" title="Fit to view">Fit</button>
    <button id="native" title="One image pixel per screen pixel">1:1</button>
    <label>Z <input id="zslider" type="range" min="0" max="0" step="1" value="0"/></label>
    <span id="zlabel"></span>
    <label id="ovwrap"><input id="ov" type="checkbox" checked/> contours
      <input id="ovop" type="range" min="0" max="1" step="0.05" value="0.9" style="width:70px"/></label>
    <label id="idwrap"><input id="ids" type="checkbox"/> IDs</label>
  </div>
  <div id="osd"></div>
  <div id="readout">pointer: –</div>
  <div id="scale"><div id="scalebar"></div><span id="scaletext"></span></div>
  <div id="note"></div>
  <div id="err"></div>
</div>
__SCRIPT__
</body>
</html>
"""

_VIEWER_JS = r"""
(function () {
  "use strict";
  const SPEC = __SPEC_JSON__;
  const SZ = SPEC.spacing_zyx_um[0], SY = SPEC.spacing_zyx_um[1], SX = SPEC.spacing_zyx_um[2];
  const zs = SPEC.plane_indices;
  const $ = (id) => document.getElementById(id);
  let zi = Math.max(0, zs.indexOf(SPEC.initial_plane));
  let overlayOn = SPEC.has_overlay, overlayOpacity = 0.9, idsOn = false, ready = false;

  function fmt(v) { return (Math.abs(v) >= 100 ? v.toFixed(1) : v.toFixed(2)); }
  function source(z, key) {
    const levels = SPEC.planes[String(z)][key];
    return {
      width: SPEC.width, height: SPEC.height, tileSize: SPEC.tile_size, tileOverlap: 0,
      minLevel: 0, maxLevel: SPEC.max_level,
      getTileUrl: function (level, x, y) { return levels[level][x + "_" + y] || SPEC.empty_tile; }
    };
  }

  if (typeof OpenSeadragon === "undefined") {
    $("err").style.display = "flex"; $("err").textContent = "OpenSeadragon failed to load."; return;
  }
  const viewer = OpenSeadragon({
    element: $("osd"), tileSources: source(zs[zi], "image"), showNavigationControl: false,
    drawer: "canvas", imageSmoothingEnabled: false, crossOriginPolicy: false, maxZoomPixelRatio: 32,
    visibilityRatio: 0.5, minZoomImageRatio: 0.5, immediateRender: true, blendTime: 0,
    gestureSettingsMouse: { clickToZoom: false, dblClickToZoom: true }, preserveViewport: true,
    showNavigator: false
  });
  window.__deepzoom = { viewer: viewer, spec: SPEC };

  function overlayItem() { return viewer.world.getItemCount() > 1 ? viewer.world.getItemAt(1) : null; }
  function applyOverlayOpacity() { const it = overlayItem(); if (it) it.setOpacity(overlayOn ? overlayOpacity : 0); }

  function markers() {
    viewer.clearOverlays();
    if (!idsOn || !SPEC.has_overlay) return;
    const item = viewer.world.getItemAt(0);
    if (!item) return;
    for (const m of SPEC.planes[String(zs[zi])].ids) {
      const el = document.createElement("div");
      el.className = "idtag"; el.textContent = String(m[0]);
      viewer.addOverlay({ element: el, location: item.imageToViewportCoordinates(m[1], m[2]),
                          placement: OpenSeadragon.Placement.CENTER, checkResize: false });
    }
  }

  function zText() {
    const z = zs[zi];
    $("zlabel").textContent = "z " + z + " / " + (SPEC.n_planes - 1) + "  (" + fmt(z * SZ) + " µm)";
  }

  function setPlane(index) {
    zi = Math.max(0, Math.min(zs.length - 1, index));
    $("zslider").value = String(zi);
    zText();
    if (!ready) return;
    const z = zs[zi];
    viewer.addTiledImage({ tileSource: source(z, "image"), index: 0, replace: true, x: 0, y: 0, width: 1,
                           success: function () { markers(); } });
    if (SPEC.has_overlay) {
      viewer.addTiledImage({ tileSource: source(z, "overlay"), index: 1, replace: true, x: 0, y: 0, width: 1,
                             opacity: overlayOn ? overlayOpacity : 0 });
    }
  }

  viewer.addOnceHandler("open", function () {
    ready = true;
    if (SPEC.has_overlay) {
      viewer.addTiledImage({ tileSource: source(zs[zi], "overlay"), x: 0, y: 0, width: 1,
                             opacity: overlayOn ? overlayOpacity : 0 });
    }
    markers(); updateScale();
  });
  viewer.addHandler("open-failed", function (e) {
    $("err").style.display = "flex"; $("err").textContent = "Could not open the image: " + (e.message || "");
  });

  function updateScale() {
    const item = viewer.world.getItemAt(0);
    if (!item) return;
    const screenPerImage = item.viewportToImageZoom(viewer.viewport.getZoom(true));
    if (!(screenPerImage > 0)) return;
    const umPerScreen = SX / screenPerImage;
    const target = 110 * umPerScreen;
    const p = Math.pow(10, Math.floor(Math.log10(target)));
    const nice = [1, 2, 5, 10].map(m => m * p).filter(v => v <= target).pop() || p;
    $("scalebar").style.width = Math.max(1, nice / umPerScreen) + "px";
    $("scaletext").textContent = (nice >= 1 ? String(+nice.toPrecision(3)) : nice.toPrecision(2)) + " µm" +
      "   (" + (screenPerImage >= 1 ? screenPerImage.toFixed(1) + " screen px / px" : (1 / screenPerImage).toFixed(1) + " px / screen px") + ")";
  }
  viewer.addHandler("animation", updateScale);
  viewer.addHandler("resize", updateScale);

  new OpenSeadragon.MouseTracker({
    element: viewer.container,
    moveHandler: function (e) {
      const item = viewer.world.getItemAt(0);
      if (!item) return;
      const p = item.viewerElementToImageCoordinates(e.position);
      const col = Math.floor(p.x), row = Math.floor(p.y), z = zs[zi];
      if (col < 0 || row < 0 || col >= SPEC.width || row >= SPEC.height) { $("readout").textContent = "pointer: outside image"; return; }
      $("readout").textContent =
        "voxel (z, y, x) = (" + z + ", " + row + ", " + col + ")\n" +
        "centre (µm)  z " + fmt(z * SZ) + "  y " + fmt(row * SY) + "  x " + fmt(col * SX);
    },
    leaveHandler: function () { $("readout").textContent = "pointer: –"; }
  });

  $("zoomin").onclick = () => { viewer.viewport.zoomBy(1.6); viewer.viewport.applyConstraints(); };
  $("zoomout").onclick = () => { viewer.viewport.zoomBy(1 / 1.6); viewer.viewport.applyConstraints(); };
  $("home").onclick = () => viewer.viewport.goHome();
  $("native").onclick = () => {
    const item = viewer.world.getItemAt(0);
    if (item) viewer.viewport.zoomTo(item.imageToViewportZoom(1), null, false);
  };
  $("zslider").max = String(zs.length - 1);
  $("zslider").disabled = zs.length < 2;
  $("zslider").oninput = (e) => setPlane(parseInt(e.target.value, 10));
  document.addEventListener("keydown", (e) => {
    if (e.key === "PageUp" || e.key === "]") { setPlane(zi + 1); e.preventDefault(); }
    if (e.key === "PageDown" || e.key === "[") { setPlane(zi - 1); e.preventDefault(); }
  });
  if (SPEC.has_overlay) {
    $("ov").onchange = (e) => { overlayOn = e.target.checked; applyOverlayOpacity(); };
    $("ovop").oninput = (e) => { overlayOpacity = parseFloat(e.target.value); applyOverlayOpacity(); };
    $("ids").onchange = (e) => { idsOn = e.target.checked; markers(); };
  } else {
    $("ovwrap").style.display = "none"; $("idwrap").style.display = "none";
  }
  const notes = [];
  if (zs.length < SPEC.n_planes) notes.push("planes " + zs[0] + "–" + zs[zs.length - 1] + " of 0–" + (SPEC.n_planes - 1) + " loaded");
  if (Math.abs(SY - SX) > 1e-9 * Math.max(SY, SX)) notes.push("pixels are not square (y " + SY + ", x " + SX + " µm); display is not aspect-corrected");
  notes.push("display only • window " + fmt(SPEC.window[0]) + "–" + fmt(SPEC.window[1]));
  $("note").textContent = notes.join(" • ");
  setPlane(zi);
})();
"""
