"""Per-object feature extraction from a labeled 3D segmentation mask.

Bridge between segmentation (Cellpose masks) and the statistical analysis
layer (``analysis``). ``mask_features`` converts a labeled ``(Z, Y, X)`` mask
plus physical voxel spacing into a per-object feature DataFrame, so the same
generic statistical helpers can operate on data measured directly from a
user's own segmentation rather than only pre-exported Excel tables.

Features are computed with physical spacing so magnitudes are in real units.
Geometry is delegated to ``quantification.features.geometry``. The default
support is the raw label, as in multilevel analysis; optional hole filling
uses the filled support consistently for every geometry feature. Exported
values retain floating-point precision; rounding belongs in presentation:
  * volume_um3       — object volume in cubic microns
  * surface_area_um2 — surface area in square microns (weighted lattice
    intersection counting, ``surface_crofton``); accuracy is claimed only
    inside the declared domain, reported per object by
    surface_in_qualified_domain / surface_rho_in / surface_anisotropy
  * equivalent_sphere_diameter_um — diameter of the sphere of equal volume
    (equivalent_disk_um is retained as a legacy alias)
  * sphericity       — ratio to the minimal possible surface for a volume
                      (not clamped; may exceed 1 on discretized voxels; values
                      above 1.05 carry a sphericity_above_geometric_range flag)
  * major_axis_um / minor_axis_um / least_axis_um — principal-axis diameters
  * elongation       — 1 - least_axis / major_axis (0 = sphere, up to ~1 = rod)
  * solidity         — volume / convex-hull volume
  * centroid_z / centroid_y / centroid_x — physical positions (µm)

All functions are pure (data in -> DataFrame / dict out) and free of printing
so they are easy to unit test.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from organoid_analysis.microscopy_io import validate_voxel_spacing_xyz
from organoid_analysis.quantification.labels import (
    bbox_touches_volume_boundary,
    compact_instance_labels,
)

# 2.1: qc_flags may contain sphericity_above_geometric_range (same limit as the
# classical route). Columns and numeric definitions are unchanged from 2.0.
FEATURE_SCHEMA_VERSION = "2.1"
FEATURE_COLUMNS = [
    "volume_um3", "surface_area_um2", "surface_rho_in", "surface_anisotropy",
    "surface_in_qualified_domain", "equivalent_disk_um", "equivalent_sphere_diameter_um",
    "sphericity", "solidity", "major_axis_um", "minor_axis_um", "least_axis_um", "elongation",
    "centroid_z_um", "centroid_y_um", "centroid_x_um", "voxel_count", "segmented_volume_um3",
    "filled_voxel_count", "measurement_basis", "touches_image_border", "fragmented_object",
    "connected_component_count", "qc_status", "qc_flags",
]


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


def _physical_spacing_zyx(spacing_um: Sequence[float]) -> tuple[float, float, float]:
    """Return (z, y, x) element spacing from the (x, y, z) API convention."""
    sx, sy, sz = validate_voxel_spacing_xyz(spacing_um)
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


def extract_mask_features(
    mask: np.ndarray, spacing_um: Sequence[float] = (1.0, 1.0, 1.0), *, fill_holes: bool = False,
) -> pd.DataFrame:
    """Return a per-object feature DataFrame for the labeled ``mask``.

    The background label ``0`` is ignored. Output rows are indexed by object
    label and only include labels that appear in the mask. ``fill_holes=True``
    explicitly selects the outer-envelope estimand. QC flags retain all rows.
    """
    from scipy import ndimage as ndi
    from skimage import measure

    from organoid_analysis.quantification.features import SPHERICITY_REVIEW_LIMIT, outer_envelope
    from organoid_analysis.quantification.features import geometry as _geometry

    if not isinstance(fill_holes, bool):
        raise ValueError("fill_holes must be a boolean")
    labels = validate_mask(mask)
    spacing_zyx = _physical_spacing_zyx(spacing_um)

    compact, source_ids = compact_instance_labels(labels)
    props = measure.regionprops(compact, spacing=spacing_zyx)
    rows = []
    for p in props:
        # Crop to the object's bounding box first: geometry() takes an origin
        # offset precisely so callers don't have to run marching cubes etc. on
        # a full-volume array per object (this was previously O(n_objects x
        # volume) and made "Object features" slow to populate on real masks).
        binary = np.ascontiguousarray(compact[p.slice] == p.label, dtype=np.uint8)
        support = outer_envelope(binary) if fill_holes else binary.astype(bool, copy=False)
        origin_zyx = tuple(s.start for s in p.slice)
        g, _ = _geometry(support, spacing_zyx, origin_zyx=origin_zyx, fill_holes=False)
        volume = g["volume_um3"]
        surface = g["surface_area_um2"]
        major = g["principal_axis_major_um"]
        minor = g["principal_axis_intermediate_um"]
        least = g["principal_axis_minor_um"]
        sphericity = g["sphericity"]
        flags = []
        touches = bbox_touches_volume_boundary(p.slice, labels.shape)
        components = int(ndi.label(binary, structure=ndi.generate_binary_structure(3, 1))[1])
        if touches:
            flags.append("border_truncated")
        if components > 1:
            flags.append("fragmented_object")
        if fill_holes and np.any(support & (compact[p.slice] != 0) & (compact[p.slice] != p.label)):
            flags.append("encloses_other_instance")
        # skimage's 3D hull needs non-coplanar voxel centers. Keep undefined
        # solidity missing with QC instead of exporting infinity or zero.
        points = np.argwhere(support)
        if len(points) < 4 or np.linalg.matrix_rank(points - points[0]) < 3:
            solidity = np.nan
        else:
            shape_prop = measure.regionprops(support.astype(np.uint8))[0] if fill_holes else p
            solidity = float(shape_prop.solidity)
        if not np.isfinite(solidity):
            flags.append("solidity_not_estimable")
        # Same review limit as the classical route; the value is kept unclipped.
        if sphericity > SPHERICITY_REVIEW_LIMIT:
            flags.append("sphericity_above_geometric_range")
        # The surface estimator's accuracy claim is domain-restricted, so an
        # object outside the domain is flagged for review rather than silently
        # carrying an unqualified area. The machine-checkable gate only; the
        # smoothness part of the domain is the caller's responsibility.
        if not g["surface_in_qualified_domain"]:
            flags.append("surface_outside_qualified_domain")
        eq_disk = g["equivalent_diameter_um"]
        elongation = 1.0 - least / major
        voxel_count = int(np.count_nonzero(binary))
        rows.append(
            {
                "Label": source_ids[int(p.label)],
                "volume_um3": volume,
                "surface_area_um2": surface,
                "surface_rho_in": g["surface_rho_in"],
                "surface_anisotropy": g["surface_anisotropy"],
                "surface_in_qualified_domain": g["surface_in_qualified_domain"],
                "equivalent_disk_um": eq_disk,
                "equivalent_sphere_diameter_um": eq_disk,
                "sphericity": sphericity,
                "solidity": solidity,
                "major_axis_um": major,
                "minor_axis_um": minor,
                "least_axis_um": least,
                "elongation": elongation,
                "centroid_z_um": g["centroid_z_um"],
                "centroid_y_um": g["centroid_y_um"],
                "centroid_x_um": g["centroid_x_um"],
                "voxel_count": voxel_count,
                "segmented_volume_um3": float(voxel_count * np.prod(spacing_zyx)),
                "filled_voxel_count": int(np.count_nonzero(support)) - voxel_count,
                "measurement_basis": "filled_envelope" if fill_holes else "raw_label",
                "touches_image_border": touches,
                "fragmented_object": components > 1,
                "connected_component_count": components,
                "qc_status": "review" if flags else "not_flagged",
                "qc_flags": ";".join(flags),
            }
        )
    feature_df = pd.DataFrame(rows, columns=["Label", *FEATURE_COLUMNS]).set_index("Label")
    feature_df.attrs["feature_schema_version"] = FEATURE_SCHEMA_VERSION
    feature_df.attrs["spacing_xyz_um"] = spacing_zyx[::-1]
    feature_df.attrs["measurement_basis"] = "filled_envelope" if fill_holes else "raw_label"
    return feature_df


def summarize_features(features: pd.DataFrame) -> MaskSummary:
    """Aggregate an existing feature table without repeating 3D geometry."""
    if features.empty:
        return MaskSummary(
            n_objects=0,
            total_volume_um3=0.0,
            mean_volume_um3=np.nan,
            median_volume_um3=np.nan,
            min_volume_um3=np.nan,
            max_volume_um3=np.nan,
            mean_sphericity=np.nan,
            mean_solidity=np.nan,
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


def summarize_mask(
    mask: np.ndarray, spacing_um: Sequence[float] = (1.0, 1.0, 1.0), *, fill_holes: bool = False,
) -> MaskSummary:
    """Extract and aggregate features in one convenience call."""
    return summarize_features(extract_mask_features(mask, spacing_um=spacing_um, fill_holes=fill_holes))


def mask_to_binary(mask: np.ndarray) -> np.ndarray:
    """Return a boolean array of the interior (label > 0) of ``mask``."""
    return validate_mask(mask) > 0
