"""Formal object-level QC flags using robust MAD volume outlier detection."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .config import Multilevel3DConfig

# Phi^-1(0.75), the 0.75-quantile of the standard normal distribution -- the
# z-form of the same MAD-to-sigma consistency-constant concept documented in
# Rousseeuw & Croux (1993), J. Am. Stat. Assoc. 88(424):1273-1283,
# https://doi.org/10.1080/01621459.1993.10476408. See docs/ALGORITHM_DECISIONS.md
# D5. `analysis/features.py` and `analysis/viability.py` use the direct form
# (... x 1.4826) of the same concept; 0.67448975 (8-decimal Phi^-1(0.75)) is
# NOT bit-identical to 1/1.4826 (5-sig-fig rounding) -- they differ by
# ~1.5ppm, immaterial in practice, but kept as this module's own literal
# rather than rewritten to `1 / MAD_TO_SIGMA` to avoid a floating-point
# change to already-validated QC flag output.
NORMAL_MAD_CONSISTENCY_Z = 0.67448975


def _volume_outliers(volumes: pd.Series, threshold: float) -> pd.Series:
    values = volumes.to_numpy(float)
    median = np.median(values) if len(values) else np.nan
    mad = np.median(np.abs(values - median)) if len(values) else np.nan
    if not np.isfinite(mad) or mad == 0:
        # A common, highly regular population can have MAD == 0 while still
        # containing a clear extreme value.  With no usable MAD scale, retain
        # the robust median centre and flag finite values that differ from it.
        return pd.Series(np.isfinite(values) & ~np.isclose(values, median), index=volumes.index)
    robust_z = NORMAL_MAD_CONSISTENCY_Z * (values - median) / mad
    return pd.Series(np.abs(robust_z) > threshold, index=volumes.index)


def add_qc_flags(features: pd.DataFrame, *, object_type: str, config: Multilevel3DConfig) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Attach compact QC status/flags and return a normalized QC flag table."""
    frame = features.copy()
    if frame.empty:
        frame["too_small"] = pd.Series(dtype=bool)
        frame["volume_outlier"] = pd.Series(dtype=bool)
        frame["qc_status"] = pd.Series(dtype=str)
        frame["qc_flags"] = pd.Series(dtype=str)
        return frame, pd.DataFrame(columns=["object_type", "object_id", "flag"])
    frame["too_small"] = frame["voxel_count"] < config.minimum_voxels
    frame["volume_outlier"] = _volume_outliers(frame["volume_um3"], config.mad_z_threshold)
    id_column = f"{object_type}_id"
    flags_by_row: list[list[str]] = []
    for row in frame.itertuples(index=False):
        flags = []
        for field, flag_name in (("touches_image_border", "touches_image_border"), ("fragmented_object", "fragmented_object"),
                                 ("too_small", "too_small"), ("volume_outlier", "volume_outlier")):
            if bool(getattr(row, field)):
                flags.append(flag_name)
        if object_type in {"cell", "nucleus"}:
            if int(getattr(row, f"{object_type}_id")) and int(getattr(row, "parent_assignment_failed")):
                flags.append("parent_assignment_failed")
            if bool(getattr(row, "low_parent_overlap")):
                flags.append("low_parent_overlap")
            if bool(getattr(row, "crosses_multiple_parents")):
                flags.append("crosses_multiple_parents")
        if object_type == "nucleus":
            if bool(getattr(row, "crosses_multiple_direct_organoids", False)):
                flags.append("crosses_multiple_direct_organoids")
            if bool(getattr(row, "direct_organoid_parent_mismatch", False)):
                flags.append("direct_organoid_parent_mismatch")
        if object_type == "cell":
            if bool(getattr(row, "anucleate")):
                flags.append("anucleate")
            if bool(getattr(row, "multinucleated")):
                flags.append("multinucleated")
        flags_by_row.append(flags)
    frame["qc_status"] = ["pass" if not flags else "flagged" for flags in flags_by_row]
    frame["qc_flags"] = [";".join(flags) for flags in flags_by_row]
    qc_rows = [{"object_type": object_type, "object_id": int(object_id), "flag": flag}
               for object_id, flags in zip(frame[id_column], flags_by_row) for flag in flags]
    return frame, pd.DataFrame(qc_rows, columns=["object_type", "object_id", "flag"])
