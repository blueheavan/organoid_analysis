"""SG-5 engineering verification: the frozen method against published implementations (D-11).

Golden values were computed by R lmerTest (Satterthwaite) and clubSandwich (CR2,
Satterthwaite) on deterministic datasets; see
docs/evidence/2026-09-26-statistical-qualification-protocol/reference/.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from organoid_analysis.statistics import small_sample as ss

REFERENCE = Path(__file__).resolve().parents[2] / "docs/evidence/2026-09-26-statistical-qualification-protocol/reference"
GOLDEN = json.loads((REFERENCE / "reference_values.json").read_text(encoding="utf-8"))
DATA = pd.read_csv(REFERENCE / "reference_datasets.csv", float_precision="round_trip")


def _data(name: str) -> ss.ClusterData:
    d = DATA[DATA.dataset == name]
    codes = {c: i for i, c in enumerate(sorted(d.condition.unique()))}
    return ss.cluster_data(d.y.to_numpy(), d.replicate.to_numpy(), d.condition.map(codes).to_numpy())


def _contrast(k: int, pair: dict) -> np.ndarray:
    c = np.zeros(k)
    c[pair["j"]], c[pair["i"]] = 1.0, -1.0
    return c


@pytest.mark.parametrize("name", sorted(GOLDEN["datasets"]))
def test_boundary_decision_matches_lme4(name):
    fit = ss.fit(_data(name))
    assert (fit.branch == ss.FALLBACK_BRANCH) == GOLDEN["datasets"][name]["singular"]


@pytest.mark.parametrize("name", sorted(n for n, v in GOLDEN["datasets"].items() if not v["singular"]))
def test_reml_and_satterthwaite_match_lmertest(name):
    golden = GOLDEN["datasets"][name]
    fit = ss.fit(_data(name))
    assert fit.sigma_e2 == pytest.approx(golden["sigma_e2"], rel=1e-6)
    assert fit.sigma_b2 == pytest.approx(golden["sigma_b2"], rel=1e-6)
    assert fit.condition_means == pytest.approx(np.array(golden["fixef"]), rel=1e-7, abs=1e-9)
    for pair in golden["pairs"]:
        result = ss.contrast(fit, _contrast(fit.data.n_conditions, pair))
        reference = pair["lmer"]
        assert result.estimate == pytest.approx(reference["estimate"], rel=1e-6, abs=1e-9)
        assert result.standard_error == pytest.approx(reference["se"], rel=1e-6)
        assert result.df == pytest.approx(reference["df"], rel=1e-5)
        assert result.p_value == pytest.approx(reference["p"], abs=1e-6)


@pytest.mark.parametrize("name", sorted(GOLDEN["datasets"]))
def test_cr2_and_its_satterthwaite_df_match_clubsandwich(name):
    golden = GOLDEN["datasets"][name]
    fit = ss.fit(_data(name))
    ols = ss.RandomInterceptFit(fit.data, fit.sigma_e2, 0.0, ss.FALLBACK_BRANCH, np.array(golden["ols"]))
    for pair in golden["pairs"]:
        result = ss.contrast(ols, _contrast(fit.data.n_conditions, pair))
        reference = pair["cr2"]
        assert result.estimate == pytest.approx(reference["estimate"], rel=1e-10, abs=1e-12)
        assert result.standard_error ** 2 == pytest.approx(reference["variance"], rel=1e-9)
        assert result.df == pytest.approx(reference["df"], rel=1e-9)
        assert result.ci_low == pytest.approx(reference["ci_low"], rel=1e-8, abs=1e-10)
        assert result.ci_high == pytest.approx(reference["ci_high"], rel=1e-8, abs=1e-10)


def test_balanced_satterthwaite_df_is_the_between_within_df():
    """Analytic case: balanced random-intercept design, between-replicate contrast => df = G - K."""
    rng = np.random.default_rng(0)
    for k, g, m in [(2, 3, 5), (3, 4, 12), (2, 10, 30)]:
        group = np.repeat(np.arange(k), g)
        data = ss.ClusterData(n=np.full(k * g, float(m)), mean=rng.normal(0, 1, k * g), group=group,
                              ssw=float(rng.chisquare(k * g * (m - 1))) * 0.5, n_conditions=k)
        fit = ss.fit(data)
        if fit.branch == ss.LMM_BRANCH:
            c = np.zeros(k)
            c[1], c[0] = 1, -1
            assert ss.contrast(fit, c).df == pytest.approx(k * g - k, rel=1e-6)


def test_object_rows_and_sufficient_statistics_give_the_same_fit():
    d = DATA[DATA.dataset == "k3_g345_unbal_icc20"]
    codes = d.condition.map({c: i for i, c in enumerate(sorted(d.condition.unique()))}).to_numpy()
    data = ss.cluster_data(d.y.to_numpy(), d.replicate.to_numpy(), codes)
    assert data.n_objects == len(d)
    assert data.ssw == pytest.approx(float(d.groupby("replicate").y.apply(lambda v: ((v - v.mean()) ** 2).sum()).sum()))
