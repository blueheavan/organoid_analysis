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
import scipy.stats as scipy_stats
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

# Features log10-transformed before fitting: volume is right-skewed on a
# multiplicative scale, which the LMM/OLS normal-residual assumption does not
# fit well untransformed. sphericity is bounded near [0,1] and is left as-is.
# Heuristic (engineering judgment); see docs/PARAMETERS.md and
# docs/ALGORITHM_DECISIONS.md D11.
LOG10_FEATURES = {"volume_um3"}

# Below this fraction of the residual scale, the LMM's random-intercept
# variance is treated as numerically collapsed to zero rather than trusted;
# heuristic numerical-stability threshold, not independently validated. See
# docs/PARAMETERS.md.
_DEGENERATE_RANDOM_EFFECT_VARIANCE_RATIO = 1e-6


def fit_model(d: pd.DataFrame, feature: str, group_col: str):
    """LMM with a replicate random intercept; fall back to OLS with
    replicate-clustered SEs when the random-effects fit is singular (variance
    collapses to ~0) or fails to converge.

    See docs/ALGORITHM_DECISIONS.md D11 for the method rationale and citation
    (Benjamini & Hochberg, 1995, J. R. Stat. Soc. B 57(1):289-300, applied to
    the pairwise contrasts in ``condition_pairwise_tests`` below).

    The OLS-fallback branch uses ``use_t=True``, which statsmodels resolves to
    a cluster-robust t(G-1) reference (G = number of replicate clusters) --
    the standard small-cluster correction (Cameron & Miller, 2015, J. Human
    Resources 50(2):317-372) instead of the anti-conservative asymptotic-normal
    default. The LMM branch has no equivalent built-in correction (statsmodels'
    ``MixedLM`` always reports z/chi2-referenced p-values); ``_pairwise_contrasts``
    below applies the same t(G-1) reference manually to that branch's contrasts.
    The LMM branch's *omnibus* Wald test (``condition_pairwise_tests``'s
    ``omnibus_p``) is NOT corrected this way and remains asymptotic -- treat it
    as a rough screening result, not a confirmatory one; the small-sample-
    corrected pairwise contrasts are the primary output.
    """
    formula = f"{feature} ~ condition"
    for method in ("lbfgs", "cg"):
        try:
            fit = smf.mixedlm(formula, d, groups=d[group_col]).fit(reml=True, method=method)
            re_var = float(np.diag(fit.cov_re).max()) if fit.cov_re.size else 0.0
            degenerate = re_var < _DEGENERATE_RANDOM_EFFECT_VARIANCE_RATIO * float(fit.scale)
            if fit.converged and np.isfinite(fit.llf) and not degenerate:
                return fit, f"LMM({method})"
        except (np.linalg.LinAlgError, ValueError):
            continue
    fit = smf.ols(formula, d).fit(cov_type="cluster", cov_kwds={"groups": d[group_col]}, use_t=True)
    return fit, "OLS(clustered SE)"


def _pairwise_contrasts(
    fit, conditions: list[str], n_clusters: int
) -> tuple[list[str], list[float], list[float], list[float], list[int]]:
    """Pairwise condition contrasts with a small-cluster-corrected p-value.

    ``fit.t_test(...)`` gives a correctly t(G-1)-referenced p-value already
    for the OLS-fallback branch (``use_t=True`` in ``fit_model``), but always
    uses an asymptotic z reference for the LMM branch (statsmodels has no
    Satterthwaite/Kenward-Roger correction for ``MixedLM``). Since the
    ``condition`` fixed effect varies only *between* replicate clusters, a
    t(G-1) reference (G = number of replicate clusters; the same convention
    statsmodels' own cluster-robust ``use_t=True`` resolves to, verified
    empirically) is a standard, defensible small-sample correction for the
    LMM branch too -- not full Satterthwaite/Kenward-Roger, but no worse than
    what the OLS-fallback branch already does. See docs/ALGORITHM_DECISIONS.md
    D11 and docs/PARAMETERS.md for the residual limitation this does not
    address (the omnibus test, and non-Satterthwaite df for LMM).
    """
    fe_names = list(fit.fe_params.index) if hasattr(fit, "fe_params") else list(fit.params.index)
    is_lmm = hasattr(fit, "cov_re")
    df_denom = max(n_clusters - 1, 1)
    reference = conditions[0]
    pairs, estimates, standard_errors, p_values, dof = [], [], [], [], []
    for i in range(len(conditions)):
        for j in range(i + 1, len(conditions)):
            ci, cj = conditions[i], conditions[j]
            coefficients = np.zeros(len(fe_names))
            if cj != reference:
                coefficients[fe_names.index(f"condition[T.{cj}]")] = 1
            if ci != reference:
                coefficients[fe_names.index(f"condition[T.{ci}]")] = -1
            test = fit.t_test(coefficients.reshape(1, -1))
            if is_lmm:
                t_stat = float(np.ravel(test.tvalue)[0])
                p_value = float(2 * scipy_stats.t.sf(abs(t_stat), df_denom))
            else:
                p_value = float(np.atleast_1d(test.pvalue)[0])
            pairs.append(f"{cj} vs {ci}")
            estimates.append(float(np.atleast_1d(test.effect)[0]))
            # test.sd is the (model-appropriate) standard error of the
            # contrast; combined with the same t(G-1) reference used for
            # p_value above (not test.conf_int(), which would use the LMM
            # branch's uncorrected asymptotic z reference) to build a CI
            # consistent with the small-sample-corrected p-value.
            standard_errors.append(float(np.ravel(test.sd)[0]))
            p_values.append(p_value)
            dof.append(df_denom)
    return pairs, estimates, standard_errors, p_values, dof


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
        n_clusters = int(d["_replicate_key"].nunique())
        wald = fit.wald_test_terms(skip_single=False)
        condition_row = [i for i, name in enumerate(wald.table.index) if "condition" in str(name)][0]
        omnibus_p = float(wald.table.iloc[condition_row]["pvalue"])

        pairs, estimates, standard_errors, p_raw, dof = _pairwise_contrasts(fit, conditions, n_clusters)
        rejected, p_adj, _, _ = multipletests(p_raw, method="fdr_bh")
        # ``estimate`` is on the model's fitted scale, which is log10 for
        # LOG10_FEATURES (volume_um3) and raw units otherwise (P2-2 audit
        # finding: the CSV previously left this scale implicit). t(G-1) is
        # the same reference used for p_raw above, not statsmodels' own
        # asymptotic-z conf_int() (see _pairwise_contrasts).
        is_log10 = feature in LOG10_FEATURES
        estimate_scale = "log10_difference" if is_log10 else "raw_difference"
        t_crit = float(scipy_stats.t.ppf(0.975, dof[0])) if dof else float("nan")
        ci_low = [e - t_crit * se for e, se in zip(estimates, standard_errors)]
        ci_high = [e + t_crit * se for e, se in zip(estimates, standard_errors)]
        pairwise_frames.append(pd.DataFrame({
            "feature": feature, "contrast": pairs, "estimate": estimates,
            "estimate_scale": estimate_scale,
            # Back-transformed ratio of geometric means between the two
            # conditions; only meaningful when the model was fit on a log10
            # scale. NaN (not 1.0 or the raw difference) for raw-scale
            # features so it is never mistaken for a ratio that was actually
            # computed.
            "geometric_mean_ratio": [10 ** e if is_log10 else np.nan for e in estimates],
            "standard_error": standard_errors,
            "ci95_low": ci_low, "ci95_high": ci_high,
            "df_denom": dof, "p_raw": p_raw, "padj": p_adj, "significant_fdr05": rejected,
            # BH-FDR is applied within one feature's own pairwise contrasts
            # only (P2-2 audit finding) -- e.g. volume_um3 and sphericity each
            # get their own independent FDR correction, not one shared
            # correction across every feature x contrast in this table.
            "multiplicity_family": f"within-feature pairwise contrasts ({feature})",
        }))
        # omnibus_p is NOT small-sample corrected for the LMM branch (see
        # fit_model's docstring); pairwise p_raw/padj above are, for both branches.
        omnibus[feature] = {"omnibus_p": omnibus_p, "model": model_used, "n_conditions": len(conditions),
                            "omnibus_p_small_sample_corrected": not model_used.startswith("LMM")}

    pairwise = (pd.concat(pairwise_frames, ignore_index=True) if pairwise_frames
                else pd.DataFrame(columns=[
                    "feature", "contrast", "estimate", "estimate_scale", "geometric_mean_ratio",
                    "standard_error", "ci95_low", "ci95_high", "df_denom", "p_raw", "padj",
                    "significant_fdr05", "multiplicity_family",
                ]))
    return pairwise, omnibus
