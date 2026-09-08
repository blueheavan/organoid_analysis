"""Batch-specific control scaling and conservative organoid-level signal states."""
from __future__ import annotations

import numpy as np
import pandas as pd

# Consistency constant rescaling MAD into a sigma-equivalent robust scale
# estimator under a Gaussian assumption (Rousseeuw & Croux, 1993, J. Am. Stat.
# Assoc. 88(424):1273-1283, https://doi.org/10.1080/01621459.1993.10476408).
# See docs/ALGORITHM_DECISIONS.md D5. `multilevel3d/qc.py` uses the reciprocal
# form (0.67448975 x ...) of the same statistical concept, at a different
# rounding precision -- 0.67448975 is Phi^-1(0.75) to 8 decimals, not bit-
# identical to 1/1.4826 (they differ by ~1.5ppm, immaterial in practice).
MAD_TO_SIGMA = 1.4826

STATES = ["viable_like", "mixed_signal", "compromised_like", "indeterminate"]
STATE_COLORS = {"viable_like": "#159D82", "mixed_signal": "#E8AE3B",
                "compromised_like": "#C44E52", "indeterminate": "#9098A3"}
CALIBRATION_COLUMNS = ["batch_id", "status", "reason", "live_control_replicates", "dead_control_replicates",
                       "live_control_units", "dead_control_units", "calcein_low", "calcein_high",
                       "pi_low", "pi_high", "calcein_separation_snr", "pi_separation_snr"]


def _mad(values) -> float:
    values = np.asarray(values, dtype=float)
    return float(MAD_TO_SIGMA * np.median(np.abs(values - np.median(values))))


def calibrate(objects: pd.DataFrame, batch_ids: list[str], cfg: dict) -> pd.DataFrame:
    rows = []
    for batch in batch_ids:
        row = {key: np.nan for key in CALIBRATION_COLUMNS}
        row.update(batch_id=batch, status="unavailable", reason="viability_mode_uncalibrated",
                   live_control_replicates=0, dead_control_replicates=0,
                   live_control_units=0, dead_control_units=0)
        if cfg["mode"] != "controls":
            rows.append(row)
            continue
        valid = objects[(objects.batch_id == batch) & objects.control.isin(["live", "dead"]) &
                        objects.morphology_eligible.astype(bool) & objects.viability_measurement_eligible.astype(bool)]
        if valid.empty:
            row["reason"] = "no_QC_eligible_controls"
            rows.append(row)
            continue
        # A biological_replicate label is only unique *within* a condition (see
        # stats.py's docstring): if this batch's live or dead controls span more
        # than one condition, grouping by (control, biological_replicate) alone
        # below could silently pool unrelated replicates that happen to share a
        # label (e.g. two condition arms each with their own "R1" live control).
        # Reject rather than guess which grouping was intended.
        ambiguous_controls = [name for name in ("live", "dead")
                              if valid.loc[valid.control == name, "condition"].nunique() > 1]
        if ambiguous_controls:
            row["reason"] = "controls_span_multiple_conditions"
            rows.append(row)
            continue
        fields = ["calcein_mean_bg_corrected", "pi_mean_bg_corrected",
                  "calcein_background_noise_mad", "pi_background_noise_mad"]
        units = valid.groupby(["control", "biological_replicate", "unit_id"], as_index=False)[fields].median()
        replicates = units.groupby(["control", "biological_replicate"], as_index=False)[fields].median()
        groups = {name: replicates[replicates.control == name] for name in ["live", "dead"]}
        for name in groups:
            row[f"{name}_control_replicates"] = len(groups[name])
            row[f"{name}_control_units"] = int((units.control == name).sum())
        if min(len(x) for x in groups.values()) < cfg["min_control_replicates"]:
            row["reason"] = "insufficient_independent_control_replicates"
            rows.append(row)
            continue
        reasons = []
        for marker, low_name, high_name in [("calcein", "dead", "live"), ("pi", "live", "dead")]:
            low_values = groups[low_name][f"{marker}_mean_bg_corrected"].to_numpy()
            high_values = groups[high_name][f"{marker}_mean_bg_corrected"].to_numpy()
            low, high = float(np.median(low_values)), float(np.median(high_values))
            # Conservative noise floor: between-replicate variability or single-voxel background noise.
            noise = max(_mad(low_values), _mad(high_values),
                        float(replicates[f"{marker}_background_noise_mad"].median()), 1e-9)
            separation = (high - low) / noise
            row.update({f"{marker}_low": low, f"{marker}_high": high,
                        f"{marker}_separation_snr": separation})
            if not np.isfinite(separation) or separation < cfg["min_control_separation_snr"]:
                reasons.append(f"{marker}_controls_not_separated")
        row["status"] = "available" if not reasons else "unavailable"
        row["reason"] = ";".join(reasons) if reasons else "control_scaled_rules_require_external_validation"
        rows.append(row)
    return pd.DataFrame(rows, columns=CALIBRATION_COLUMNS)


def classify(objects: pd.DataFrame, calibration: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    result = objects.copy()
    result["calcein_control_scaled"] = np.nan
    result["pi_control_scaled"] = np.nan
    result["viability_state"] = "indeterminate"
    result["viability_reason"] = ""
    result["viability_method"] = "unclassified"
    lookup = {str(row["batch_id"]): row for row in calibration.to_dict("records")}
    for index, row in result.iterrows():
        if not bool(row.morphology_eligible):
            result.at[index, "viability_reason"] = "morphology_QC_failed"
            continue
        if not bool(row.viability_measurement_eligible):
            result.at[index, "viability_reason"] = row.viability_measurement_flags
            continue
        controls = lookup.get(str(row.batch_id), {})
        if controls.get("status") != "available":
            result.at[index, "viability_reason"] = controls.get("reason", "missing_batch_controls")
            continue
        c = (row.calcein_mean_bg_corrected - controls["calcein_low"]) / (controls["calcein_high"] - controls["calcein_low"])
        p = (row.pi_mean_bg_corrected - controls["pi_low"]) / (controls["pi_high"] - controls["pi_low"])
        result.at[index, "calcein_control_scaled"] = c
        result.at[index, "pi_control_scaled"] = p
        result.at[index, "viability_method"] = "batch_control_scaled_rules"
        if not np.isfinite([c, p]).all():
            result.at[index, "viability_reason"] = "nonfinite_scaled_signal"
        elif c < cfg["low_gate"] and p < cfg["low_gate"]:
            result.at[index, "viability_reason"] = "both_markers_low"
        elif c >= cfg["high_gate"] and p <= cfg["low_gate"]:
            result.at[index, "viability_state"] = "viable_like"
            result.at[index, "viability_reason"] = "high_calcein_low_PI"
        elif p >= cfg["high_gate"] and c <= cfg["low_gate"]:
            result.at[index, "viability_state"] = "compromised_like"
            result.at[index, "viability_reason"] = "high_PI_low_calcein"
        else:
            result.at[index, "viability_state"] = "mixed_signal"
            result.at[index, "viability_reason"] = "intermediate_or_discordant_marker_signals"
        if c < 0 or c > 1 or p < 0 or p > 1:
            result.at[index, "viability_reason"] += ";outside_control_range"
    return result
