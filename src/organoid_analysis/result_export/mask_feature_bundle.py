"""In-memory, traceable CSV exports for a selected nucleus/cell mask."""
from __future__ import annotations

import hashlib
import importlib.metadata
import io
import json
import zipfile
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from organoid_analysis import __version__
from organoid_analysis.microscopy_io.tiff_contract import git_commit_hash, source_code_hashes
from organoid_analysis.microscopy_io.voxel_spacing import validate_voxel_spacing_xyz
from organoid_analysis.quantification.mask_features import FEATURE_SCHEMA_VERSION, validate_mask


def build_mask_feature_bundle(
    mask: np.ndarray,
    features: pd.DataFrame,
    *,
    spacing_um: Sequence[float],
    object_type: str,
    fill_holes: bool,
    segmentation_config: Mapping[str, object],
    segmentation_provenance: Mapping[str, object] | None = None,
) -> bytes:
    """Package full-precision measurements with explicit coordinate/QC semantics.

    No files are written; callers can cache the returned ZIP until the selected
    mask, measurement definition or run changes. Unknown lineage stays unknown.
    """
    labels = validate_mask(mask)
    spacing = validate_voxel_spacing_xyz(spacing_um)
    if not isinstance(fill_holes, bool):
        raise ValueError("fill_holes must be a boolean")
    if object_type not in {"nucleus", "cell"}:
        raise ValueError("Feature exports require an explicit nucleus or cell object_type")
    basis = "filled_envelope" if fill_holes else "raw_label"
    if tuple(features.attrs.get("spacing_xyz_um", ())) != spacing:
        raise ValueError("Feature calibration does not match the requested voxel spacing")
    present = {int(value) for value in np.unique(labels) if value != 0}
    if features.index.has_duplicates or set(features.index) != present:
        raise ValueError("Feature rows do not match the selected mask's instance IDs")
    if features.attrs.get("measurement_basis") != basis or (not features.empty and not features["measurement_basis"].eq(basis).all()):
        raise ValueError("Feature geometry does not match the requested measurement basis")
    csv_data = features.to_csv(index=True, index_label="Label", float_format="%.17g").encode("utf-8")
    mask_header = {"shape_zyx": list(labels.shape), "dtype": labels.dtype.str, "order": "C"}
    digest = hashlib.sha256(json.dumps(mask_header, sort_keys=True).encode("ascii"))
    digest.update(memoryview(np.ascontiguousarray(labels)).cast("B"))
    provenance = {
        "feature_schema_version": FEATURE_SCHEMA_VERSION,
        "pipeline_version": __version__,
        "git_commit": git_commit_hash(),
        "measurement_source_code_sha256": source_code_hashes(),
        "packages": {name: importlib.metadata.version(name) for name in ("numpy", "scipy", "scikit-image", "pandas")},
        "object_type": object_type,
        "measurement_basis": basis,
        "spacing_xyz_um": list(spacing),
        "calibration_sources": {name: segmentation_config.get(name, "unknown") for name in ("xy_spacing_source", "anisotropy_source")},
        "calibration_status": (
            "provided_not_independently_verified"
            if all(segmentation_config.get(name) in {"metadata", "user_override"} for name in ("xy_spacing_source", "anisotropy_source"))
            else "assumed_or_unknown"
        ),
        "coordinates": "array-local voxel centers; origin ZYX=(0,0,0); micrometres",
        "surface_method": "skimage Lewiner marching cubes; level=0.5; native spacing; step_size=1; zero padding",
        "surface_accuracy_status": "NOT ASSESSED for this mask; existing analytical <5% criterion fails",
        "axis_columns": {"major_axis_um": "largest moment-equivalent diameter", "minor_axis_um": "intermediate diameter (legacy name)", "least_axis_um": "smallest diameter"},
        "legacy_aliases": {"equivalent_disk_um": "equivalent_sphere_diameter_um"},
        "qc_policy": "All objects retained; review flags do not imply biological abnormality or validated accuracy",
        "n_objects": len(features),
        "n_review_objects": int(features["qc_status"].eq("review").sum()),
        "statistical_scope": "One field; descriptive object measurements; experimental independence not assessed",
        "selected_mask": {**mask_header, "array_sha256": digest.hexdigest(), "hash_encoding": "SHA256(sorted-key JSON header followed by C-order array bytes)"},
        "files": {"features.csv": {"sha256": hashlib.sha256(csv_data).hexdigest()}},
        "segmentation_config": dict(segmentation_config),
        "segmentation_provenance": dict(segmentation_provenance) if segmentation_provenance is not None else None,
    }
    metadata = (json.dumps(provenance, indent=2, allow_nan=False) + "\n").encode("utf-8")
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        for name, content in (("features.csv", csv_data), ("measurement_provenance.json", metadata)):
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(info, content)
    return buffer.getvalue()
