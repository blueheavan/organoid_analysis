"""Fields -> wells -> biological replicates, without treating organoids as independent n."""
from __future__ import annotations

import numpy as np
import pandas as pd

from organoid_analysis.phenotyping.viability import STATES

MORPHOLOGY = ["volume_um3", "surface_area_um2", "sphericity", "equivalent_diameter_um"]
META = ["condition", "biological_replicate", "unit_id", "batch_id", "control"]

# Minimum biological replicates with size data required before a bootstrap CI
# is computed. Coincides with (but is a separate, non-configurable literal
# from) stats.min_replicates_per_condition -- see docs/PARAMETERS.md.
_MIN_REPLICATES_FOR_BOOTSTRAP_CI = 3


def _describe(objects: pd.DataFrame) -> dict:
    eligible = objects[objects.morphology_eligible.astype(bool)]
    # Counts are int; medians, totals and fractions are float (nan when empty).
    result: dict[str, int | float] = {"n_detected": len(objects), "n_included": len(eligible),
                                      "n_excluded": len(objects) - len(eligible)}
    for metric in MORPHOLOGY:
        result[f"median_{metric}"] = float(eligible[metric].median()) if len(eligible) else np.nan
    result["total_volume_um3"] = float(eligible.volume_um3.sum())
    for state in STATES:
        result[f"n_{state}"] = int((eligible.viability_state == state).sum())
        result[f"fraction_{state}"] = result[f"n_{state}"] / len(eligible) if len(eligible) else np.nan
    return result


def make_summaries(objects: pd.DataFrame, samples: pd.DataFrame, cfg: dict) -> tuple:
    sample_rows = []
    for meta in samples.to_dict("records"):
        selected = objects[objects.sample_id == meta["sample_id"]]
        sample_rows.append({**meta, **_describe(selected)})
    sample_summary = pd.DataFrame(sample_rows)
    unit_rows = []
    for keys, fields in sample_summary.groupby(["batch_id", "unit_id"], sort=False, dropna=False):
        meta = fields.iloc[0]
        selected = objects[(objects.batch_id == keys[0]) & (objects.unit_id == keys[1])]
        successful = fields.status.eq("ok")
        unit_rows.append({**{key: meta[key] for key in META}, "n_fields": len(fields),
                          "n_successful_fields": int(successful.sum()),
                          "unit_complete": bool(successful.all()),
                          **_describe(selected)})
    units = pd.DataFrame(unit_rows)
    replicate_rows = []
    # biological_replicate is GLOBAL across batches: repeated acquisition does not create a new independent n.
    for keys, wells in units.groupby(["condition", "biological_replicate", "control"], sort=False, dropna=False):
        complete = wells[wells.unit_complete]
        nonempty = complete[complete.n_included > 0]
        row = dict(zip(["condition", "biological_replicate", "control"], keys))
        row.update(n_units=len(wells), n_complete_units=len(complete), n_nonempty_units=len(nonempty),
                   n_included=int(complete.n_included.sum()), n_detected=int(complete.n_detected.sum()),
                   n_fields=int(complete.n_fields.sum()),
                   mean_detected_count_per_field=float(complete.n_detected.sum() / complete.n_fields.sum()) if len(complete) else np.nan)
        for metric in MORPHOLOGY:
            # Each well receives equal weight, independent of its number of organoids or fields.
            row[f"median_of_unit_medians_{metric}"] = float(nonempty[f"median_{metric}"].median()) if len(nonempty) else np.nan
        for state in STATES:
            row[f"mean_unit_fraction_{state}"] = float(nonempty[f"fraction_{state}"].mean()) if len(nonempty) else np.nan
        replicate_rows.append(row)
    replicates = pd.DataFrame(replicate_rows)
    rng = np.random.default_rng(cfg["seed"])
    condition_rows = []
    for keys, reps in replicates.groupby(["condition", "control"], sort=False, dropna=False):
        available = reps[reps.n_nonempty_units > 0]
        row = {"condition": keys[0], "control": keys[1], "n_biological_replicates_total": len(reps),
               "n_biological_replicates_with_size_data": len(available),
               "n_included_organoids_complete_units": int(reps.n_included.sum()),
               "n_complete_units": int(reps.n_complete_units.sum())}
        for metric in MORPHOLOGY:
            values = available[f"median_of_unit_medians_{metric}"].dropna().to_numpy(float)
            row[f"mean_replicate_median_{metric}"] = float(values.mean()) if len(values) else np.nan
            low, high = np.nan, np.nan
            if len(values) >= _MIN_REPLICATES_FOR_BOOTSTRAP_CI:
                bootstrap = rng.choice(values, (cfg["bootstrap_iterations"], len(values)), replace=True).mean(axis=1)
                low, high = np.quantile(bootstrap, [.025, .975])
            row[f"ci95_low_{metric}"] = low
            row[f"ci95_high_{metric}"] = high
        for state in STATES:
            row[f"mean_replicate_fraction_{state}"] = float(available[f"mean_unit_fraction_{state}"].mean()) if len(available) else np.nan
        condition_rows.append(row)
    return sample_summary, units, replicates, pd.DataFrame(condition_rows)
