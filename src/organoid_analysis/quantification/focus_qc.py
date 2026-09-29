"""Z-coverage quality control from focus and signal profiles along Z.

A 3D measurement is only complete if the acquired Z range contains the whole
object. An object cut by the first or last plane is already flagged
``border_truncated`` when its *mask* reaches that plane; this module catches
the cases the mask can miss: a dim object whose segmented boundary stops short
of the stack edge while its signal does not, and widefield / bright-field
stacks whose sharpest plane lies at the edge of the scanned range (Kim et al.
2026, Cell Rep. Methods 6:101601, pick the in-focus bright-field plane with a
7x7 Laplacian contrast metric and note that later-stage organoids exceed a
single plane's depth of field).

Per object, two profiles over Z are computed on the raw structure channel:

* focus  -- variance of the Laplacian-of-Gaussian response (sigma 1.5 px,
  support comparable to a 7x7 Laplacian kernel), a standard sharpness measure,
  within the XY projection of the mask dilated by a margin (edges carry most
  of the sharpness);
* signal -- median intensity minus the stack's background level, within the
  undilated XY projection.

Other objects are excluded plane by plane, dilated in XY by the same margin
(the LoG response of a bright neighbour spills into the margin), so an
organoid lying above or below in the same XY position does not contaminate
this object's profiles.

Flags (review only; never excludes an object or changes a measurement):

* ``z_focus_peak_at_edge``    -- the sharpest plane is the first or last plane;
* ``z_signal_at_edge``        -- an edge plane still carries >= EDGE_SIGNAL_RATIO
  of the peak background-corrected signal;
* ``z_coverage_suspect``      -- either of the above, or the mask touches the
  first/last plane.

EDGE_SIGNAL_RATIO and the focus scale are HEURISTIC review thresholds, not
validated acceptance criteria.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi  # type: ignore[import-untyped]

FOCUS_SIGMA_PX = 1.5
FOOTPRINT_MARGIN_PX = 5
EDGE_SIGNAL_RATIO = 0.5
MIN_PLANES = 3
Z_QC_COLUMNS = ["z_focus_peak_plane", "z_focus_peak_at_edge", "z_edge_signal_ratio", "z_signal_at_edge",
                "z_coverage_suspect", "z_coverage_flags"]
_PLANAR = ndi.generate_binary_structure(2, 1)[np.newaxis]  # XY-only dilation, no Z coupling
STACK_QC_COLUMNS = ["stack_focus_peak_plane", "stack_focus_peak_at_edge", "stack_edge_signal_ratio",
                    "stack_z_coverage_suspect"]


def _region(footprint: np.ndarray | None, exclude: np.ndarray | None, z: int, shape: tuple[int, ...]) -> np.ndarray:
    region = np.ones(shape, bool) if footprint is None else footprint.copy()
    if exclude is not None:
        region &= ~exclude[z]
    return region


def focus_profile(stack: np.ndarray, footprint: np.ndarray | None = None,
                  exclude: np.ndarray | None = None) -> np.ndarray:
    """Variance of the LoG response per Z plane inside a 2D footprint, minus excluded voxels.

    ``exclude`` (ZYX bool) removes other objects, so an organoid lying above or
    below in the same XY position does not contaminate this object's profile.
    """
    values = np.empty(stack.shape[0])
    for z in range(stack.shape[0]):
        response = ndi.gaussian_laplace(np.asarray(stack[z], np.float32), FOCUS_SIGMA_PX)
        region = response[_region(footprint, exclude, z, response.shape)]
        values[z] = float(region.var()) if region.size else 0.0
    return values


def signal_profile(stack: np.ndarray, footprint: np.ndarray | None, background: float,
                   exclude: np.ndarray | None = None) -> np.ndarray:
    values = np.empty(stack.shape[0])
    for z in range(stack.shape[0]):
        region = stack[z][_region(footprint, exclude, z, stack.shape[1:])]
        values[z] = float(np.median(region)) - background if region.size else 0.0
    return values


def _edge_ratio(signal: np.ndarray) -> float:
    peak = float(signal.max())
    if peak <= 0:
        return float("nan")
    return float(max(signal[0], signal[-1], 0.0) / peak)


def _classify(focus: np.ndarray, signal: np.ndarray, touches_z: bool) -> dict[str, object]:
    peak = int(np.argmax(focus))
    at_edge = bool(peak in (0, len(focus) - 1) and np.ptp(focus) > 0)
    ratio = _edge_ratio(signal)
    signal_edge = bool(np.isfinite(ratio) and ratio >= EDGE_SIGNAL_RATIO)
    flags = [name for name, on in (("z_focus_peak_at_edge", at_edge), ("z_signal_at_edge", signal_edge),
                                   ("mask_touches_z_boundary", touches_z)) if on]
    return {"z_focus_peak_plane": peak, "z_focus_peak_at_edge": at_edge, "z_edge_signal_ratio": ratio,
            "z_signal_at_edge": signal_edge, "z_coverage_suspect": bool(flags), "z_coverage_flags": ";".join(flags)}


def stack_background(stack: np.ndarray, labels: np.ndarray | None) -> float:
    """Median of voxels outside every object (or of the whole stack when no labels are given)."""
    if labels is not None and (labels == 0).any():
        return float(np.median(stack[labels == 0]))
    return float(np.median(stack))


def object_z_coverage(stack: np.ndarray, labels: np.ndarray) -> dict[int, dict[str, object]]:
    """Z-coverage QC for every label ID in ``labels`` (same ZYX grid as ``stack``)."""
    if stack.shape != labels.shape or stack.ndim != 3:
        raise ValueError("structure stack and labels must share one 3D ZYX grid")
    results: dict[int, dict[str, object]] = {}
    if stack.shape[0] < MIN_PLANES:
        return {int(i): {"z_focus_peak_plane": -1, "z_focus_peak_at_edge": False, "z_edge_signal_ratio": float("nan"),
                         "z_signal_at_edge": False, "z_coverage_suspect": True, "z_coverage_flags": "too_few_planes"}
                for i in np.unique(labels[labels > 0])}
    background = stack_background(stack, labels)
    for label_id, bbox in enumerate(ndi.find_objects(labels), start=1):
        if bbox is None:
            continue
        ys = slice(max(bbox[1].start - FOOTPRINT_MARGIN_PX, 0), min(bbox[1].stop + FOOTPRINT_MARGIN_PX, stack.shape[1]))
        xs = slice(max(bbox[2].start - FOOTPRINT_MARGIN_PX, 0), min(bbox[2].stop + FOOTPRINT_MARGIN_PX, stack.shape[2]))
        crop = stack[:, ys, xs]
        crop_labels = labels[:, ys, xs]
        own = crop_labels == label_id
        projection = own.any(axis=0)
        footprint = ndi.binary_dilation(projection, iterations=FOOTPRINT_MARGIN_PX)
        others = (crop_labels != 0) & ~own
        if others.any():
            others = ndi.binary_dilation(others, _PLANAR, iterations=FOOTPRINT_MARGIN_PX) & ~own
        touches_z = bool(bbox[0].start == 0 or bbox[0].stop == stack.shape[0])
        results[label_id] = _classify(focus_profile(crop, footprint, others),
                                      signal_profile(crop, projection, background, others), touches_z)
    return results


def stack_z_coverage(stack: np.ndarray, labels: np.ndarray | None = None) -> dict[str, object]:
    """Whole-field Z coverage: does the stack's sharpest plane or signal sit at an edge?"""
    if stack.ndim != 3 or stack.shape[0] < MIN_PLANES:
        return {"stack_focus_peak_plane": -1, "stack_focus_peak_at_edge": False,
                "stack_edge_signal_ratio": float("nan"), "stack_z_coverage_suspect": True}
    footprint = (labels > 0).any(axis=0) if labels is not None and (labels > 0).any() else None
    result = _classify(focus_profile(stack, footprint),
                       signal_profile(stack, footprint, stack_background(stack, labels)), False)
    return {"stack_focus_peak_plane": result["z_focus_peak_plane"],
            "stack_focus_peak_at_edge": result["z_focus_peak_at_edge"],
            "stack_edge_signal_ratio": result["z_edge_signal_ratio"],
            "stack_z_coverage_suspect": result["z_coverage_suspect"]}
