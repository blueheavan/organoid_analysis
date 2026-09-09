"""Automatic inference of segmentation parameters from a Z-stack.

Two parameters previously set by hand can often be derived from the data:

* **Anisotropy** (``anisotropy = z_step / xy_pixel``) — read from the TIFF's
  OME-XML or ImageJ metadata via :mod:`organoid_analysis.microscopy_io.metadata`. This is a physical
  imaging property that cannot (and should not) be guessed from pixel content;
  when the metadata is present we trust it, otherwise we report it unknown and
  fall back to a manual value.

* **Nuclei / cell diameter** (pixels) — estimated from a fast low-resolution
  2D pre-segmentation: we run the Cellpose network on a middle Z-slice at
  reduced resolution, then measure the median object diameter (the same
  statistic Cellpose uses for its own auto-rescaling, ``cellpose.utils.
  diameters``). The estimate is scaled back to full resolution.

All functions are pure-ish and free of UI/printing so they are easy to unit
test, mirroring the style of ``organoid_analysis.segmentation.cellpose_inference``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile

from organoid_analysis.microscopy_io.metadata import resolve_spacing

# Slices used for the fast diameter pre-estimate (one Z-spaced through the stack).
_ESTIMATE_SLICES = 6
# Downsample factor for the fast 2D pre-segmentation (diameter is not
# resolution-critical, so a coarse pass is plenty and much faster).
_ESTIMATE_DOWNSAMPLE = 0.5
# Coverage quantile used to pick the "typical" object size; 0.5 (median)
# is chosen for robustness to background noise objects.
_DIAMETER_QUANTILE = 0.5

# Cellpose call parameters for the fast diameter *pre-estimate* only (never
# for the final segmentation, whose parameters the user controls directly via
# SegmentationConfig). flow_threshold=0.0/cellprob_threshold=0.0 are
# permissive Cellpose-conventional defaults chosen to avoid under-detecting
# candidate objects at this coarse pre-pass; min_size=5 (px) discards
# single-pixel noise detections; batch_size=8 is an engineering
# throughput/memory choice with no scientific effect. See docs/PARAMETERS.md.
_PRESEGMENTATION_FLOW_THRESHOLD = 0.0
_PRESEGMENTATION_CELLPROB_THRESHOLD = 0.0
_PRESEGMENTATION_MIN_SIZE_PX = 5
_PRESEGMENTATION_BATCH_SIZE = 8

# Plausible full-resolution diameter range (px) for a nucleus; pre-estimate
# diameters outside this range are treated as noise blobs and discarded
# before taking the quantile. Heuristic, not independently calibrated for a
# specific instrument. See docs/PARAMETERS.md.
_PLAUSIBLE_DIAMETER_PX_RANGE = (5.0, 250.0)


def auto_anisotropy(path: str | Path) -> float | None:
    """Return the Z/XY anisotropy from TIFF metadata, or ``None`` if unknown.

    Uses the same OME-XML -> ImageJ -> TIFF-resolution resolution order as the
    reader. Returns ``None`` (never a guessed default) when the physical Z step
    or the XY pixel size cannot be resolved.
    """
    with tifffile.TiffFile(path) as tf:
        ome_metadata = tf.ome_metadata
        imagej_metadata = getattr(tf, "imagej_metadata", None)
        tiff_tags = {}
        for series in tf.series:
            for page in series.levels[0]:
                for key in ("XResolution", "YResolution", "ResolutionUnit"):
                    tag = page.tags.get(key)
                    if tag is not None:
                        tiff_tags[key] = tag.value
                break
            break
    spacing = resolve_spacing(ome_metadata, imagej_metadata, tiff_tags)
    if not spacing.complete:
        return None
    return spacing.anisotropy


def _downsample_2d(slice_image: np.ndarray, factor: float) -> tuple[np.ndarray, float]:
    from scipy import ndimage

    factor = float(factor)
    if factor >= 1.0:
        return slice_image, 1.0
    new_shape = (max(1, int(round(slice_image.shape[0] * factor))),
                 max(1, int(round(slice_image.shape[1] * factor))))
    zoomed = ndimage.zoom(
        slice_image.astype(np.float32),
        (new_shape[0] / slice_image.shape[0],
         new_shape[1] / slice_image.shape[1]),
        order=1,
    )
    return zoomed, factor


def estimate_diameter_from_stack(
    model,
    stack: np.ndarray,
    *,
    downsample: float = _ESTIMATE_DOWNSAMPLE,
    max_slices: int = _ESTIMATE_SLICES,
) -> float | None:
    """Estimate the median object diameter (pixels, full resolution).

    Runs a fast 2D segmentation on a handful of Z-slices (uniformly spaced
    through the stack), measures the median equivalent diameter of the objects,
    and scales the result back up to the full-resolution pixel scale.

    Returns ``None`` when no objects are found.
    """
    from cellpose import utils

    if stack.ndim != 3:
        raise ValueError(f"Expected a (Z, Y, X) stack, got shape {stack.shape}.")
    z = stack.shape[0]
    if max_slices <= 0:
        max_slices = 1
    n_slices = min(max_slices, z)
    indices = np.linspace(0, z - 1, n_slices).round().astype(int)
    indices = np.unique(indices)

    all_diams: list[float] = []
    for zi in indices:
        slice_image = stack[zi]
        work, actual = _downsample_2d(slice_image, downsample)
        work = work[..., np.newaxis] if work.ndim == 2 else work
        try:
            masks = model.eval(
                work, diameter=None,
                flow_threshold=_PRESEGMENTATION_FLOW_THRESHOLD,
                cellprob_threshold=_PRESEGMENTATION_CELLPROB_THRESHOLD,
                min_size=_PRESEGMENTATION_MIN_SIZE_PX, do_3D=False,
                batch_size=_PRESEGMENTATION_BATCH_SIZE,
            )[0]
        except Exception:  # noqa: BLE001
            continue
        if masks is None:
            continue
        masks = masks.astype(np.int32)
        if masks.max() < 1:
            continue
        try:
            _md, diams = utils.diameters(masks)
        except Exception:  # noqa: BLE001
            continue
        if diams is None or diams.size == 0:
            continue
        # diams are in downsampled pixels; scale back to full resolution.
        full = diams / actual
        all_diams.extend(float(d) for d in full)

    if not all_diams:
        return None
    arr = np.asarray(all_diams)
    # Filter to the plausible range for a nucleus so noise blobs don't skew it.
    low, high = _PLAUSIBLE_DIAMETER_PX_RANGE
    arr = arr[(arr >= low) & (arr <= high)]
    if arr.size == 0:
        return None
    return float(np.quantile(arr, _DIAMETER_QUANTILE))


def suggest_config(
    model,
    stack: np.ndarray,
    *,
    anisotropy: float | None = None,
):
    """Return a dict of auto-suggested parameters.

    Brutal but safe: ``anisotropy`` (physical) is returned verbatim when given
    or when inferable from the file; ``nuclei_diameter`` is estimated from the
    stack. ``cell_diameter`` defaults to the nucleus estimate (a cell is roughly
    one nucleus plus cytoplasm), which the user can still override.
    """
    est = estimate_diameter_from_stack(model, stack)
    result = {"nuclei_diameter": est, "cell_diameter": est}
    if anisotropy is not None:
        result["anisotropy"] = anisotropy
    return result
