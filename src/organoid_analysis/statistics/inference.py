"""Cross-condition significance testing on top of the descriptive summaries in
``aggregation.py``.

For each continuous morphology feature the frozen method of
``statistics.small_sample`` (``METHOD_ID``; owner decisions D-9, D-11, D-16) is
applied: a random-intercept model with the biological replicate as the
independent experimental unit, REML variance components, Satterthwaite degrees
of freedom for every pairwise contrast, and the CR2 cluster-robust fallback with
Satterthwaite degrees of freedom when the replicate variance is at the REML
boundary. Contrasts are Benjamini-Hochberg corrected within a feature. The
omnibus F uses between-within denominator degrees of freedom and is a screen.

Each contrast carries ``interval_qualification``: ``qualified`` only when the
design lies inside the envelope where the SG-5 coverage study met its
predeclared tolerances (``small_sample.design_qualification``), otherwise
``provisional - coverage unqualified``.

``biological_replicate`` labels are only unique *within* a condition, so the
grouping key is the pair (condition, biological_replicate).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from organoid_analysis.statistics import small_sample

# Features log10-transformed before fitting: volume is right-skewed on a
# multiplicative scale. Heuristic (engineering judgment); see
# docs/PARAMETERS.md and docs/ALGORITHM_DECISIONS.md D11.
LOG10_FEATURES = {"volume_um3"}


def fit_model(d: pd.DataFrame, feature: str, group_col: str) -> tuple[small_sample.RandomInterceptFit, str]:
    """Fit the frozen random-intercept method to object rows; returns (fit, branch label)."""
    conditions = list(d["condition"].cat.categories)
    codes = d["condition"].cat.codes.to_numpy()
    data = small_sample.cluster_data(d[feature].to_numpy(float), d[group_col].to_numpy(), codes)
    if data.n_conditions != len(conditions):
        raise ValueError("every retained condition must contribute replicates")
    fit = small_sample.fit(data)
    return fit, fit.branch


def _pairwise_contrasts(fit: small_sample.RandomInterceptFit, conditions: list[str]) -> list[dict]:
    rows = []
    for i in range(len(conditions)):
        for j in range(i + 1, len(conditions)):
            c = np.zeros(len(conditions))
            c[j], c[i] = 1.0, -1.0
            result = small_sample.contrast(fit, c)
            qualified, reason = small_sample.design_qualification(fit.data, (i, j))
            rows.append({"contrast": f"{conditions[j]} vs {conditions[i]}", "estimate": result.estimate,
                         "standard_error": result.standard_error, "ci95_low": result.ci_low,
                         "ci95_high": result.ci_high, "df_denom": result.df, "p_raw": result.p_value,
                         "branch": result.branch,
                         "interval_qualification": "qualified" if qualified else
                         "provisional - coverage unqualified",
                         "qualification_reason": reason})
    return rows


def condition_pairwise_tests(
    objects: pd.DataFrame,
    features: tuple[str, ...] = ("volume_um3", "sphericity"),
    min_replicates_per_condition: int = 3,
) -> tuple[pd.DataFrame, dict]:
    """Omnibus + BH-FDR pairwise condition contrasts for each feature.

    Returns ``(pairwise, omnibus)``: a tidy DataFrame with one row per
    (feature, contrast) and a ``{feature: {"omnibus_p", "model", "n_conditions"}}``
    dict. A feature/condition combination with too few biological replicates
    to fit is silently excluded rather than raising -- callers should treat an
    empty result as "not enough data to test", not an error.
    """
    eligible = objects[objects.morphology_eligible.astype(bool)].copy()
    # Factorize the pair without delimiter collisions (e.g. A::B/C vs A/B::C).
    eligible["_replicate_key"] = eligible.groupby(
        ["condition", "biological_replicate"], sort=False, dropna=False
    ).ngroup()

    pairwise_frames = []
    omnibus: dict[str, dict] = {}
    for feature in features:
        if feature not in eligible.columns:
            continue
        d = eligible[[feature, "condition", "_replicate_key"]].dropna().copy()
        if not np.isfinite(d[feature].to_numpy(dtype=float)).all():
            raise ValueError(f"{feature} contains nonfinite eligible measurements")
        if feature in LOG10_FEATURES and (d[feature] <= 0).any():
            raise ValueError(f"{feature} must be positive for log10 inference; no clipping is permitted")
        replicate_counts = d.groupby("condition")["_replicate_key"].nunique()
        keep_conditions = replicate_counts[replicate_counts >= min_replicates_per_condition].index
        d = d[d.condition.isin(keep_conditions)]
        conditions = sorted(d.condition.unique())
        if len(conditions) < 2:
            continue
        if feature in LOG10_FEATURES:
            d[feature] = np.log10(d[feature])
        d["condition"] = pd.Categorical(d.condition, categories=conditions)

        try:
            fit, model_used = fit_model(d, feature, "_replicate_key")
        except ValueError:
            # No residual degrees of freedom or no variation: not estimable, not an error.
            continue
        statistic, df_num, df_den, omnibus_p = small_sample.omnibus(fit)
        rows = _pairwise_contrasts(fit, conditions)
        if not np.isfinite([row["p_raw"] for row in rows]).all():
            # Zero replicate-level variation within conditions: no test is defined.
            continue
        _, p_adj, _, _ = multipletests([row["p_raw"] for row in rows], method="fdr_bh")
        is_log10 = feature in LOG10_FEATURES
        frame = pd.DataFrame(rows)
        frame.insert(0, "feature", feature)
        frame.insert(3, "estimate_scale", "log10_difference" if is_log10 else "raw_difference")
        # Back-transformed ratio of geometric means; NaN for raw-scale features so
        # it is never mistaken for a ratio that was actually computed.
        frame.insert(4, "geometric_mean_ratio", [10 ** e if is_log10 else np.nan for e in frame.estimate])
        frame["padj"] = p_adj
        frame["significant_fdr05"] = frame.padj < 0.05
        frame["multiplicity_family"] = f"within-feature pairwise contrasts ({feature})"
        frame["method_id"] = small_sample.METHOD_ID
        pairwise_frames.append(frame)
        omnibus[feature] = {"omnibus_p": omnibus_p, "model": model_used, "n_conditions": len(conditions),
                            "omnibus_F": statistic, "omnibus_df": [df_num, df_den],
                            "omnibus_p_small_sample_corrected": False,
                            "omnibus_note": "Wald F with between-within denominator df (G - K); screening only, "
                                            "anti-conservative when G < 10 (owner decision D-9)",
                            "method_id": small_sample.METHOD_ID,
                            "sigma_b2": fit.sigma_b2, "sigma_e2": fit.sigma_e2}

    pairwise = (pd.concat(pairwise_frames, ignore_index=True) if pairwise_frames
                else pd.DataFrame(columns=[
                    "feature", "contrast", "estimate", "estimate_scale", "geometric_mean_ratio",
                    "standard_error", "ci95_low", "ci95_high", "df_denom", "p_raw", "branch",
                    "interval_qualification", "qualification_reason", "padj", "significant_fdr05",
                    "multiplicity_family", "method_id",
                ]))
    return pairwise, omnibus
