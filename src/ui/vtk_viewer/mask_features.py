"""Per-object feature extraction from a labeled 3D segmentation mask.

Bridge between segmentation (Cellpose masks) and the statistical analysis
layer (``analysis``). ``mask_features`` converts a labeled ``(Z, Y, X)`` mask
plus physical voxel spacing into a per-object feature DataFrame, so the same
generic statistical helpers can operate on data measured directly from a
user's own segmentation rather than only pre-exported Excel tables.

Features are computed with physical spacing so magnitudes are in real units.
Geometry (volume, surface area, sphericity, principal-axis diameters) is
delegated to ``analysis.features.geometry`` so the UI metrics are bit-for-bit
consistent with the canonical analysis pipeline instead of a duplicated
implementation:
  * volume_um3       — object volume in cubic microns
  * surface_area_um2 — isosurface area in square microns (marching cubes)
  * equivalent_disk  — diameter (µm) of the sphere of equal volume
  * sphericity       — ratio to the minimal possible surface for a volume
                      (not clamped; may slightly exceed 1 on discretized voxels)
  * major_axis_um / minor_axis_um / least_axis_um — principal-axis diameters
  * elongation       — 1 - least_axis / major_axis (0 = sphere, up to ~1 = rod)
  * solidity         — volume / convex-hull volume
  * centroid_z / centroid_y / centroid_x — physical positions (µm)

All functions are pure (data in -> DataFrame / dict out) and free of printing
so they are easy to unit test.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .core import validate_spacing


@dataclass(frozen=True)
class MaskSummary:
    """Aggregate statistics for one labeled mask volume."""
    n_objects: int
    total_volume_um3: float
    mean_volume_um3: float
    median_volume_um3: float
    min_volume_um3: float
    max_volume_um3: float
    mean_sphericity: float
    mean_solidity: float


def _physical_spacing_zyx(spacing_um) -> tuple[float, float, float]:
    """Return (z, y, x) element spacing from the (x, y, z) API convention."""
    sx, sy, sz = validate_spacing(spacing_um)
    return (sz, sy, sx)


def validate_mask(mask: np.ndarray) -> np.ndarray:
    """Validate a nonnegative integer-label ``(Z, Y, X)`` mask.

    Fractional floating masks must not be silently truncated into plausible
    instance IDs. Valid integer arrays are returned without an int64 copy,
    which matters for large uint16/uint32 segmentation volumes.
    """
    if not isinstance(mask, np.ndarray):
        mask = np.asarray(mask)
    if mask.ndim != 3:
        raise ValueError(f"Expected a 3D labeled mask (Z, Y, X), got shape {mask.shape}.")
    if min(mask.shape) < 2:
        raise ValueError(f"Mask is too small for feature extraction: {mask.shape}.")
    if mask.dtype == np.bool_:
        return mask
    if not np.issubdtype(mask.dtype, np.integer):
        raise ValueError("Mask must contain integer instance labels, not probabilities or intensities.")
    if np.issubdtype(mask.dtype, np.signedinteger) and np.any(mask < 0):
        raise ValueError("Mask labels must be nonnegative.")
    if mask.size and int(mask.max()) > np.iinfo(np.uint32).max:
        raise ValueError("Mask labels exceed the supported uint32 range.")
    return mask


def count_mask_objects(mask: np.ndarray) -> int:
    """Count labels that actually occur; sparse IDs are not object counts."""
    labels = validate_mask(mask)
    values = np.unique(labels)
    return int(np.count_nonzero(values))


def extract_mask_features(mask: np.ndarray, spacing_um=(1.0, 1.0, 1.0)) -> pd.DataFrame:
    """Return a per-object feature DataFrame for the labeled ``mask``.

    The background label ``0`` is ignored. Output rows are indexed by object
    label and only include labels that appear in the mask.
    """
    from skimage import measure

    from analysis.features import geometry as _geometry

    labels = validate_mask(mask)
    spacing_zyx = _physical_spacing_zyx(spacing_um)

    props = measure.regionprops(labels, spacing=spacing_zyx)
    rows = []
    for p in props:
        # Geometry (volume/surface/sphericity/axes) is computed by the canonical
        # analysis pipeline on the filled envelope so the UI and the pipeline
        # agree exactly; centroids come from regionprops (full-image coords).
        # Crop to the object's bounding box first: geometry() takes an origin
        # offset precisely so callers don't have to run marching cubes etc. on
        # a full-volume array per object (this was previously O(n_objects x
        # volume) and made "Object features" slow to populate on real masks).
        binary = np.ascontiguousarray(labels[p.slice] == p.label, dtype=np.uint8)
        origin_zyx = tuple(s.start for s in p.slice)
        g, _ = _geometry(binary, spacing_zyx, origin_zyx=origin_zyx)
        volume = g["volume_um3"]
        surface = g["surface_area_um2"]
        major = g["principal_axis_major_um"]
        minor = g["principal_axis_intermediate_um"]
        least = g["principal_axis_minor_um"]
        sphericity = g["sphericity"]
        # regionprops.centroid is ordered (z, y, x).
        cz, cy, cx = (float(v) for v in p.centroid)
        solidity = float(p.solidity) if p.solidity is not None else 0.0
        eq_disk = (6.0 * volume / np.pi) ** (1.0 / 3.0)
        elongation = 1.0 - (least / major) if major > 0 else 0.0
        rows.append(
            {
                "Label": int(p.label),
                "volume_um3": round(volume, 4),
                "surface_area_um2": round(surface, 4),
                "equivalent_disk_um": round(eq_disk, 4),
                "sphericity": round(sphericity, 4),
                "solidity": round(solidity, 4),
                "major_axis_um": round(major, 4),
                "minor_axis_um": round(minor, 4),
                "least_axis_um": round(least, 4),
                "elongation": round(elongation, 4),
                "centroid_z_um": round(cz, 4),
                "centroid_y_um": round(cy, 4),
                "centroid_x_um": round(cx, 4),
            }
        )
    feature_df = pd.DataFrame(rows)
    if feature_df.empty:
        return feature_df
    return feature_df.set_index("Label")


def summarize_features(features: pd.DataFrame) -> MaskSummary:
    """Aggregate an existing feature table without repeating 3D geometry."""
    if features.empty:
        return MaskSummary(
            n_objects=0,
            total_volume_um3=0.0,
            mean_volume_um3=0.0,
            median_volume_um3=0.0,
            min_volume_um3=0.0,
            max_volume_um3=0.0,
            mean_sphericity=0.0,
            mean_solidity=0.0,
        )
    return MaskSummary(
        n_objects=int(features.shape[0]),
        total_volume_um3=float(features["volume_um3"].sum()),
        mean_volume_um3=float(features["volume_um3"].mean()),
        median_volume_um3=float(features["volume_um3"].median()),
        min_volume_um3=float(features["volume_um3"].min()),
        max_volume_um3=float(features["volume_um3"].max()),
        mean_sphericity=float(features["sphericity"].mean()),
        mean_solidity=float(features["solidity"].mean()),
    )


def summarize_mask(mask: np.ndarray, spacing_um=(1.0, 1.0, 1.0)) -> MaskSummary:
    """Extract and aggregate features in one convenience call."""
    return summarize_features(extract_mask_features(mask, spacing_um=spacing_um))


def mask_to_binary(mask: np.ndarray) -> np.ndarray:
    """Return a boolean array of the interior (label > 0) of ``mask``."""
    return validate_mask(mask) > 0
