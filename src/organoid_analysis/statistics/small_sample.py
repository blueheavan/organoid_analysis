"""Frozen small-sample inference for condition contrasts (SG-5; owner decisions D-9, D-11, D-16).

Model
-----
``y_ij = mu_{k(i)} + b_i + e_ij`` with replicate random intercept
``b_i ~ N(0, sigma_b^2)`` and object residual ``e_ij ~ N(0, sigma_e^2)``. The
condition varies only between biological replicates, which are the independent
experimental unit; objects are nested measurements. Because every fixed effect
is constant within a replicate, the REML likelihood depends on the data only
through the per-replicate object counts ``n_i``, means ``ybar_i`` and the pooled
within-replicate sum of squares, so the fit below works on those sufficient
statistics exactly (no approximation).

Primary branch (D-9): REML fit; each contrast of condition means is tested with
Satterthwaite degrees of freedom computed as in lmerTest (Kuznetsova, Brockhoff
& Christensen 2017, J. Stat. Softw. 82(13), doi:10.18637/jss.v082.i13):
``df = 2 V^2 / (grad V' A grad V)`` with ``V`` the contrast variance as a
function of the variance parameters and ``A = 2 H^-1`` their asymptotic
covariance from the Hessian ``H`` of the REML deviance.

Fallback branch (D-11): when the REML random-intercept variance is at the
boundary (``sigma_b^2 < 1e-6 sigma_e^2``) the object-level OLS fit is used with
the CR2 cluster-robust covariance and its Satterthwaite degrees of freedom under
a working model of independent homoskedastic errors (Bell & McCaffrey 2002;
Pustejovsky & Tipton 2018, J. Bus. Econ. Stat. 36(4):672-683,
doi:10.1080/07350015.2016.1247004), as implemented in the R package
clubSandwich.

Omnibus (D-9): Wald F on all condition means with between-within denominator
degrees of freedom ``G - K`` (screening only; anti-conservative at small G is
stated with every result).

Verification against lmerTest and clubSandwich golden values:
``tests/statistics/test_small_sample_reference.py``; coverage qualification:
``docs/evidence/2026-09-26-statistical-qualification-protocol``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize  # type: ignore[import-untyped]
from scipy import stats as scipy_stats

METHOD_ID = "reml_satterthwaite_cr2_fallback/1"
BOUNDARY_RATIO = 1e-6  # sigma_b^2 / sigma_e^2 below which the random intercept is treated as absent
LMM_BRANCH = "LMM-REML(Satterthwaite)"

# Design envelope inside which the SG-5 coverage study met its predeclared
# tolerances (docs/evidence/2026-09-26-statistical-qualification-protocol).
# Outside it every interval is labelled provisional. Set from the development
# simulation before the confirmation run and never lowered afterwards (D-10).
QUALIFIED_MIN_REPLICATES = 6                # per compared condition
QUALIFIED_MAX_REPLICATES = 21                 # largest simulated per-condition count
QUALIFIED_MAX_CONDITIONS = 3                  # largest simulated number of conditions
QUALIFIED_MEAN_OBJECTS = (5.0, 150.0)         # simulated range of mean objects per replicate
QUALIFIED_MIN_OBJECTS = 2                     # smallest simulated objects in a replicate
FALLBACK_BRANCH = "OLS-CR2(Satterthwaite)"


@dataclass(frozen=True)
class ClusterData:
    """Sufficient statistics of a nested design: one entry per biological replicate."""
    n: np.ndarray        # objects per replicate
    mean: np.ndarray     # replicate mean of the analysed feature
    group: np.ndarray    # condition index 0..K-1 of each replicate
    ssw: float           # pooled within-replicate sum of squares
    n_conditions: int

    @property
    def n_objects(self) -> int:
        return int(self.n.sum())

    @property
    def n_replicates(self) -> int:
        return int(self.n.size)


def cluster_data(values: np.ndarray, replicate: np.ndarray, condition: np.ndarray) -> ClusterData:
    """Reduce object rows to replicate sufficient statistics (condition index must be constant per replicate)."""
    values = np.asarray(values, float)
    replicate_codes, replicate_index = np.unique(np.asarray(replicate), return_inverse=True)
    condition = np.asarray(condition)
    n = np.bincount(replicate_index).astype(float)
    sums = np.bincount(replicate_index, weights=values)
    means = sums / n
    ssw = float(((values - means[replicate_index]) ** 2).sum())
    group = np.empty(replicate_codes.size, int)
    for i in range(replicate_codes.size):
        levels = np.unique(condition[replicate_index == i])
        if levels.size != 1:
            raise ValueError("a biological replicate spans more than one condition")
        group[i] = int(levels[0])
    return ClusterData(n=n, mean=means, group=group, ssw=ssw, n_conditions=int(group.max()) + 1)


@dataclass(frozen=True)
class RandomInterceptFit:
    data: ClusterData
    sigma_e2: float
    sigma_b2: float
    branch: str
    condition_means: np.ndarray   # GLS means (LMM) or object-weighted means (fallback)

    @property
    def icc(self) -> float:
        return self.sigma_b2 / (self.sigma_b2 + self.sigma_e2)


# ------------------------------------------------------------------- REML
def _profile(rho: float, d: ClusterData) -> tuple[float, float, np.ndarray]:
    """REML deviance profiled over sigma_e^2 at intraclass ratio lambda = rho / (1 - rho)."""
    lam = rho / (1.0 - rho)
    s = lam + 1.0 / d.n
    w = 1.0 / s
    weight_sum = np.bincount(d.group, weights=w, minlength=d.n_conditions)
    mu = np.bincount(d.group, weights=w * d.mean, minlength=d.n_conditions) / weight_sum
    q = float((w * (d.mean - mu[d.group]) ** 2).sum())
    a = (d.ssw + q) / (d.n_objects - d.n_conditions)
    deviance = (d.n_objects - d.n_conditions) * np.log(a) + np.log(s).sum() + np.log(weight_sum).sum()
    return float(deviance), float(a), mu


def fit(d: ClusterData) -> RandomInterceptFit:
    """REML fit of the random-intercept model; boundary fits switch to the CR2 fallback branch."""
    if d.n_objects - d.n_conditions <= 0 or d.ssw + np.ptp(d.mean) <= 0:
        raise ValueError("the design has no residual degrees of freedom or no variation")
    if d.n_objects == d.n_replicates:
        # One object per replicate: sigma_b^2 and sigma_e^2 are not separately
        # identifiable (the REML profile is flat), so the fallback branch is used.
        spread = float(((d.mean - _object_means(d)[d.group]) ** 2).sum())
        return _fallback(d, spread / max(d.n_objects - d.n_conditions, 1), 0.0)
    grid = np.concatenate([[0.0], np.linspace(0.0, 1.0, 81)[1:-1] ** 2, [1.0 - 1e-9]])
    values = [_profile(r, d)[0] for r in grid]
    best = int(np.argmin(values))
    lo, hi = grid[max(best - 1, 0)], grid[min(best + 1, grid.size - 1)]
    rho = float(grid[best])
    if hi > lo:
        found = optimize.minimize_scalar(lambda r: _profile(r, d)[0], bounds=(lo, hi), method="bounded",
                                         options={"xatol": 1e-13, "maxiter": 500})
        if found.fun <= values[best]:
            rho = float(found.x)
    _, a, mu = _profile(rho, d)
    lam = rho / (1.0 - rho)
    b = lam * a
    if lam < BOUNDARY_RATIO:
        return _fallback(d, a, b)
    return RandomInterceptFit(d, a, b, LMM_BRANCH, mu)


def _object_means(d: ClusterData) -> np.ndarray:
    n_k = np.bincount(d.group, weights=d.n, minlength=d.n_conditions)
    return np.bincount(d.group, weights=d.n * d.mean, minlength=d.n_conditions) / n_k


def _fallback(d: ClusterData, a: float, b: float) -> RandomInterceptFit:
    return RandomInterceptFit(d, a, b, FALLBACK_BRANCH, _object_means(d))


def _tau(d: ClusterData, a: float, b: float) -> np.ndarray:
    return b + a / d.n


def _weight_sums(d: ClusterData, a: float, b: float) -> np.ndarray:
    return np.bincount(d.group, weights=1.0 / _tau(d, a, b), minlength=d.n_conditions)


def reml_gradient(d: ClusterData, a: float, b: float) -> np.ndarray:
    """Gradient of the REML deviance (-2 log L) with respect to (sigma_e^2, sigma_b^2)."""
    tau = _tau(d, a, b)
    w = 1.0 / tau
    weight_sum = np.bincount(d.group, weights=w, minlength=d.n_conditions)
    mu = np.bincount(d.group, weights=w * d.mean, minlength=d.n_conditions) / weight_sum
    r = d.mean - mu[d.group]
    common = w - w ** 2 / weight_sum[d.group] - r ** 2 * w ** 2
    g_a = (d.n_objects - d.n_replicates) / a - d.ssw / a ** 2 + float((common / d.n).sum())
    g_b = float(common.sum())
    return np.array([g_a, g_b])


def reml_hessian(d: ClusterData, a: float, b: float) -> np.ndarray:
    """Central finite differences of the analytic gradient (relative step 1e-5)."""
    theta = np.array([a, b])
    hessian = np.empty((2, 2))
    for j in range(2):
        step = 1e-5 * max(abs(theta[j]), 1e-8 * a)
        up, down = theta.copy(), theta.copy()
        up[j] += step
        down[j] -= step
        hessian[:, j] = (reml_gradient(d, *up) - reml_gradient(d, *down)) / (2 * step)
    return 0.5 * (hessian + hessian.T)


# --------------------------------------------------------------- contrasts
@dataclass(frozen=True)
class ContrastResult:
    estimate: float
    standard_error: float
    df: float
    t: float
    p_value: float
    ci_low: float
    ci_high: float
    branch: str


def _finish(estimate: float, variance: float, df: float, branch: str) -> ContrastResult:
    se = float(np.sqrt(variance))
    t = estimate / se if se > 0 else np.nan
    p = float(2 * scipy_stats.t.sf(abs(t), df)) if np.isfinite(t) and df > 0 else np.nan
    half = float(scipy_stats.t.ppf(0.975, df)) * se if df > 0 else np.nan
    return ContrastResult(float(estimate), se, float(df), float(t), p, estimate - half, estimate + half, branch)


def _satterthwaite(model: RandomInterceptFit, c: np.ndarray) -> ContrastResult:
    d, a, b = model.data, model.sigma_e2, model.sigma_b2
    tau = _tau(d, a, b)
    weight_sum = _weight_sums(d, a, b)
    variance = float((c ** 2 / weight_sum).sum())
    factor = (c ** 2 / weight_sum ** 2)[d.group] / tau ** 2
    grad = np.array([float((factor / d.n).sum()), float(factor.sum())])
    covariance = 2.0 * np.linalg.inv(reml_hessian(d, a, b))
    df = 2.0 * variance ** 2 / float(grad @ covariance @ grad)
    return _finish(float(c @ model.condition_means), variance, df, LMM_BRANCH)


def _cr2_terms(model: RandomInterceptFit, c: np.ndarray) -> tuple[float, float]:
    """CR2 variance of c'mu and its Satterthwaite df (working model: independent, homoskedastic)."""
    d = model.data
    n_k = np.bincount(d.group, weights=d.n, minlength=d.n_conditions)
    h = d.n / n_k[d.group]
    if np.any(h >= 1.0):
        raise ValueError("CR2 needs at least two replicates per compared condition")
    residual_sum = d.n * (d.mean - model.condition_means[d.group])
    coefficient = c[d.group] / n_k[d.group]
    variance = float((coefficient ** 2 * residual_sum ** 2 / (1.0 - h)).sum())
    s = coefficient / np.sqrt(1.0 - h)
    same = d.group[:, None] == d.group[None, :]
    cross = np.where(same, (d.n[:, None] * d.n[None, :]) / n_k[d.group][:, None], 0.0)
    gram = s[:, None] * s[None, :] * (np.diag(d.n) - cross)
    df = float(np.trace(gram) ** 2 / (gram ** 2).sum())
    return variance, df


def contrast(model: RandomInterceptFit, c: np.ndarray) -> ContrastResult:
    """Test and 95% interval for the contrast c' mu of condition means."""
    c = np.asarray(c, float)
    if model.branch == LMM_BRANCH:
        return _satterthwaite(model, c)
    variance, df = _cr2_terms(model, c)
    return _finish(float(c @ model.condition_means), variance, df, FALLBACK_BRANCH)


def omnibus(model: RandomInterceptFit) -> tuple[float, float, float, float]:
    """Wald F for equality of all condition means: (F, df_num, df_den, p). df_den = G - K."""
    d = model.data
    k = d.n_conditions
    if model.branch == LMM_BRANCH:
        covariance = np.diag(1.0 / _weight_sums(d, model.sigma_e2, model.sigma_b2))
    else:
        n_k = np.bincount(d.group, weights=d.n, minlength=k)
        h = d.n / n_k[d.group]
        residual_sum = d.n * (d.mean - model.condition_means[d.group])
        covariance = np.diag(np.bincount(d.group, weights=residual_sum ** 2 / (1.0 - h), minlength=k) / n_k ** 2)
    contrast_matrix = np.hstack([-np.ones((k - 1, 1)), np.eye(k - 1)])
    delta = contrast_matrix @ model.condition_means
    try:
        wald = float(delta @ np.linalg.solve(contrast_matrix @ covariance @ contrast_matrix.T, delta))
    except np.linalg.LinAlgError:
        return np.nan, float(k - 1), float(d.n_replicates - k), np.nan  # no replicate-level variation
    df_num, df_den = float(k - 1), float(d.n_replicates - k)
    statistic = wald / df_num
    p = float(scipy_stats.f.sf(statistic, df_num, df_den)) if df_den > 0 else np.nan
    return statistic, df_num, df_den, p


def design_qualification(d: ClusterData, pair: tuple[int, int]) -> tuple[bool, str]:
    """Whether a contrast's design lies inside the qualified envelope, with the first reason if not."""
    counts = np.bincount(d.group, minlength=d.n_conditions)
    mean_objects = d.n_objects / d.n_replicates
    checks = [
        (d.n_conditions <= QUALIFIED_MAX_CONDITIONS, f"{d.n_conditions} conditions > {QUALIFIED_MAX_CONDITIONS}"),
        (min(counts[list(pair)]) >= QUALIFIED_MIN_REPLICATES,
         f"fewer than {QUALIFIED_MIN_REPLICATES} replicates in a compared condition"),
        (max(counts) <= QUALIFIED_MAX_REPLICATES, f"more than {QUALIFIED_MAX_REPLICATES} replicates in a condition"),
        (QUALIFIED_MEAN_OBJECTS[0] <= mean_objects <= QUALIFIED_MEAN_OBJECTS[1],
         f"mean objects per replicate {mean_objects:.1f} outside {QUALIFIED_MEAN_OBJECTS}"),
        (int(d.n.min()) >= QUALIFIED_MIN_OBJECTS, f"a replicate has fewer than {QUALIFIED_MIN_OBJECTS} objects"),
    ]
    for ok, reason in checks:
        if not ok:
            return False, reason
    return True, "inside the SG-5 qualified design envelope"
