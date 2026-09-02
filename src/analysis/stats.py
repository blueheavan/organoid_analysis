"""Cross-condition significance testing, on top of the descriptive summaries in
``summary.py`` (which report replicate-weighted medians and bootstrap CIs but
never a p-value).

For each continuous morphology feature: a linear mixed-effects model (condition
as a fixed effect, biological replicate as a random intercept) gives an omnibus
Wald test plus all-pairwise contrasts between conditions, BH-FDR corrected. The
LMM falls back to OLS with replicate-clustered standard errors when the random
effect is singular or the fit does not converge (small biological-replicate
counts make this common).

``biological_replicate`` labels are only unique *within* a condition (the same
"R1" label is reused by every condition for the well acquired in that batch
position, see ``summary.py``'s ``META`` grouping) -- so the LMM/clustering group
key here is the composite ``condition::biological_replicate``, not the bare
``biological_replicate`` column.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

LOG10_FEATURES = {"volume_um3"}


def fit_model(d: pd.DataFrame, feature: str, group_col: str):
    """LMM with a replicate random intercept; fall back to OLS with
    replicate-clustered SEs when the random-effects fit is singular (variance
    collapses to ~0) or fails to converge."""
    formula = f"{feature} ~ condition"
    for method in ("lbfgs", "cg"):
        try:
            fit = smf.mixedlm(formula, d, groups=d[group_col]).fit(reml=True, method=method)
            re_var = float(np.diag(fit.cov_re).max()) if fit.cov_re.size else 0.0
            degenerate = re_var < 1e-6 * float(fit.scale)
            if fit.converged and np.isfinite(fit.llf) and not degenerate:
                return fit, f"LMM({method})"
        except (np.linalg.LinAlgError, ValueError):
            continue
    fit = smf.ols(formula, d).fit(cov_type="cluster", cov_kwds={"groups": d[group_col]})
    return fit, "OLS(clustered SE)"


def _pairwise_contrasts(fit, conditions: list[str]) -> tuple[list[str], list[float], list[float]]:
    fe_names = list(fit.fe_params.index) if hasattr(fit, "fe_params") else list(fit.params.index)
    reference = conditions[0]
    pairs, estimates, p_values = [], [], []
    for i in range(len(conditions)):
        for j in range(i + 1, len(conditions)):
            ci, cj = conditions[i], conditions[j]
            coefficients = np.zeros(len(fe_names))
            if cj != reference:
                coefficients[fe_names.index(f"condition[T.{cj}]")] = 1
            if ci != reference:
                coefficients[fe_names.index(f"condition[T.{ci}]")] = -1
            test = fit.t_test(coefficients.reshape(1, -1))
            pairs.append(f"{cj} vs {ci}")
            estimates.append(float(np.atleast_1d(test.effect)[0]))
            p_values.append(float(np.atleast_1d(test.pvalue)[0]))
    return pairs, estimates, p_values


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
    eligible["_replicate_key"] = eligible.condition.astype(str) + "::" + eligible.biological_replicate.astype(str)

    pairwise_frames = []
    omnibus: dict[str, dict] = {}
    for feature in features:
        if feature not in eligible.columns:
            continue
        d = eligible[[feature, "condition", "_replicate_key"]].dropna().copy()
        replicate_counts = d.groupby("condition")["_replicate_key"].nunique()
        keep_conditions = replicate_counts[replicate_counts >= min_replicates_per_condition].index
        d = d[d.condition.isin(keep_conditions)]
        conditions = sorted(d.condition.unique())
        if len(conditions) < 2:
            continue
        if feature in LOG10_FEATURES:
            d[feature] = np.log10(d[feature].clip(lower=np.finfo(float).tiny))
        d["condition"] = pd.Categorical(d.condition, categories=conditions)

        fit, model_used = fit_model(d, feature, "_replicate_key")
        wald = fit.wald_test_terms(skip_single=False)
        condition_row = [i for i, name in enumerate(wald.table.index) if "condition" in str(name)][0]
        omnibus_p = float(wald.table.iloc[condition_row]["pvalue"])

        pairs, estimates, p_raw = _pairwise_contrasts(fit, conditions)
        rejected, p_adj, _, _ = multipletests(p_raw, method="fdr_bh")
        pairwise_frames.append(pd.DataFrame({
            "feature": feature, "contrast": pairs, "estimate": estimates,
            "p_raw": p_raw, "padj": p_adj, "significant_fdr05": rejected,
        }))
        omnibus[feature] = {"omnibus_p": omnibus_p, "model": model_used, "n_conditions": len(conditions)}

    pairwise = (pd.concat(pairwise_frames, ignore_index=True) if pairwise_frames
                else pd.DataFrame(columns=["feature", "contrast", "estimate", "p_raw", "padj", "significant_fdr05"]))
    return pairwise, omnibus
