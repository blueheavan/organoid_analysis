import numpy as np
import pandas as pd
import pytest
import scipy.stats as st
import tifffile
import yaml

from organoid_analysis.config import load_config
from organoid_analysis.statistics import small_sample
from organoid_analysis.statistics.inference import (
    _pairwise_contrasts,
    condition_pairwise_tests,
    fit_model,
)
from organoid_analysis.workflows.organoid_measurement_workflow import analyze


def _objects(condition_values: dict, n_replicates: int = 4, feature: str = "volume_um3") -> pd.DataFrame:
    """Build a minimal objects table: one row per (condition, replicate)."""
    rng = np.random.default_rng(0)
    rows = []
    for condition, (mean, noise) in condition_values.items():
        for replicate in range(n_replicates):
            rows.append({
                "condition": condition,
                "biological_replicate": f"R{replicate}",
                "morphology_eligible": True,
                feature: float(mean + rng.normal(0, noise)),
            })
    return pd.DataFrame(rows)


def test_real_effect_is_significant_and_null_feature_is_not():
    rng = np.random.default_rng(1)
    rows = []
    for condition, volume_mean in (("A", 1000.0), ("B", 3000.0)):
        for replicate in range(4):
            rows.append({
                "condition": condition, "biological_replicate": f"R{replicate}",
                "morphology_eligible": True,
                "volume_um3": float(volume_mean + rng.normal(0, 50)),
                "sphericity": float(0.8 + rng.normal(0, 0.01)),  # no real difference
            })
    objects = pd.DataFrame(rows)
    pairwise, omnibus = condition_pairwise_tests(objects, features=("volume_um3", "sphericity"))

    assert omnibus["volume_um3"]["omnibus_p"] < 0.01
    assert omnibus["sphericity"]["omnibus_p"] > 0.05
    volume_contrast = pairwise[pairwise.feature == "volume_um3"].iloc[0]
    assert volume_contrast.padj < 0.05
    assert bool(volume_contrast.significant_fdr05)


def test_bh_fdr_padj_never_below_raw_p():
    objects = _objects({"A": (1000.0, 100.0), "B": (1200.0, 100.0), "C": (2000.0, 100.0)})
    pairwise, _ = condition_pairwise_tests(objects, features=("volume_um3",))
    assert (pairwise.padj >= pairwise.p_raw - 1e-12).all()


def test_pairwise_output_is_self_describing_and_scale_correct():
    """P2-2 audit finding: pairwise_contrasts.csv's ``estimate`` for
    volume_um3 is a log10-scale difference but the column name alone did not
    say so, and the CSV gave no indication that BH-FDR is applied within one
    feature's contrasts only, not globally across every feature."""
    objects = _objects(
        {"A": (1000.0, 50.0), "B": (2000.0, 50.0)}, feature="volume_um3",
    )
    pairwise, _ = condition_pairwise_tests(objects, features=("volume_um3",))
    row = pairwise.iloc[0]

    assert row.estimate_scale == "log10_difference"
    assert row.multiplicity_family == "within-feature pairwise contrasts (volume_um3)"
    # geometric_mean_ratio must back-transform the log10-scale estimate, not
    # just restate the difference or a hardcoded 1.0.
    assert row.geometric_mean_ratio == pytest.approx(10 ** row.estimate)
    assert row.ci95_low < row.estimate < row.ci95_high
    assert row.standard_error > 0
    assert row.ci95_high - row.ci95_low == pytest.approx(
        2 * st.t.ppf(0.975, row.df_denom) * row.standard_error
    )


def test_raw_scale_feature_reports_no_geometric_mean_ratio():
    """sphericity is not log10-transformed; a geometric-mean ratio computed
    from its raw-scale difference would be scientifically meaningless, so it
    must be reported as NaN rather than a plausible-looking number."""
    objects = _objects(
        {"A": (0.80, 0.01), "B": (0.85, 0.01)}, feature="sphericity",
    )
    pairwise, _ = condition_pairwise_tests(objects, features=("sphericity",))
    row = pairwise.iloc[0]
    assert row.estimate_scale == "raw_difference"
    assert np.isnan(row.geometric_mean_ratio)


def _nested(seed: int = 7, offset_sd: float = 25.0, conditions=("A", "B"), replicates: int = 3,
            objects: int = 6) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for index, condition in enumerate(conditions):
        for replicate in range(replicates):
            cluster_offset = rng.normal(0, offset_sd)
            for _ in range(objects):
                rows.append({"condition": condition, "grp": f"{condition}::R{replicate}",
                             "volume_um3": float(1000.0 + 15.0 * index + cluster_offset + rng.normal(0, 15))})
    d = pd.DataFrame(rows)
    d["condition"] = pd.Categorical(d.condition, categories=list(conditions))
    return d


def test_lmm_branch_uses_satterthwaite_df_not_asymptotic_z():
    """D-9: the LMM contrast is referenced to t with Satterthwaite df, never to z.

    For a balanced design with a condition that varies only between replicates,
    the Satterthwaite df equals the between-within df G - K exactly.
    """
    d = _nested()
    fit, model_used = fit_model(d, "volume_um3", "grp")
    assert model_used == small_sample.LMM_BRANCH
    row = _pairwise_contrasts(fit, ["A", "B"])[0]
    assert row["df_denom"] == pytest.approx(d["grp"].nunique() - 2, rel=1e-6)
    t_stat = row["estimate"] / row["standard_error"]
    assert row["p_raw"] == pytest.approx(2 * st.t.sf(abs(t_stat), row["df_denom"]))
    assert row["p_raw"] > 2 * st.norm.sf(abs(t_stat))


def test_boundary_random_effect_switches_to_cr2_fallback():
    d = _nested(offset_sd=0.0, seed=3)
    # Replicate means exactly equal within a condition: the REML replicate variance is 0.
    d["volume_um3"] = d.volume_um3 - d.groupby("grp", observed=True).volume_um3.transform("mean") \
        + d.groupby("condition", observed=True).volume_um3.transform("mean")
    fit, model_used = fit_model(d, "volume_um3", "grp")
    assert model_used == small_sample.FALLBACK_BRANCH
    row = _pairwise_contrasts(fit, ["A", "B"])[0]
    assert row["branch"] == small_sample.FALLBACK_BRANCH
    assert 0 < row["df_denom"] <= d["grp"].nunique() - 1


def test_constant_feature_is_not_estimable_and_skipped():
    rows = [{"condition": c, "biological_replicate": f"R{r}", "morphology_eligible": True, "sphericity": 0.9}
            for c in ("A", "B") for r in range(3)]
    pairwise, omnibus = condition_pairwise_tests(pd.DataFrame(rows), features=("sphericity",))
    assert pairwise.empty and omnibus == {}


def test_designs_outside_the_qualified_envelope_are_labelled_provisional():
    small = _nested(replicates=2)
    fit, _ = fit_model(small, "volume_um3", "grp")
    row = _pairwise_contrasts(fit, ["A", "B"])[0]
    assert row["interval_qualification"] == "provisional - coverage unqualified"
    many = _nested(conditions=("A", "B", "C", "D"), replicates=small_sample.QUALIFIED_MIN_REPLICATES)
    fit, _ = fit_model(many, "volume_um3", "grp")
    assert all(r["interval_qualification"] != "qualified" for r in _pairwise_contrasts(fit, ["A", "B", "C", "D"]))
    inside = _nested(replicates=small_sample.QUALIFIED_MIN_REPLICATES)
    fit, _ = fit_model(inside, "volume_um3", "grp")
    assert _pairwise_contrasts(fit, ["A", "B"])[0]["interval_qualification"] == "qualified"


def test_condition_with_too_few_replicates_is_excluded():
    objects = _objects({"A": (1000.0, 50.0), "B": (2000.0, 50.0)}, n_replicates=4)
    rare = _objects({"C": (1500.0, 50.0)}, n_replicates=2)
    objects = pd.concat([objects, rare], ignore_index=True)
    pairwise, omnibus = condition_pairwise_tests(objects, features=("volume_um3",), min_replicates_per_condition=3)
    assert omnibus["volume_um3"]["n_conditions"] == 2
    contrasts = set(pairwise.contrast)
    assert not any("C" in c for c in contrasts)


def test_single_condition_is_skipped_without_error():
    objects = _objects({"A": (1000.0, 50.0)})
    pairwise, omnibus = condition_pairwise_tests(objects, features=("volume_um3",))
    assert omnibus == {}
    assert pairwise.empty


def _write_field(tmp_path, name: str, mask_scale: float) -> pd.Series:
    z, y, x = np.indices((20, 36, 36))
    mask = ((z - 10) * 2) ** 2 + (y - 18) ** 2 + (x - 18) ** 2 < (10 * mask_scale) ** 2
    data = np.full((3, 20, 36, 36), 100, np.uint16)
    data[0, mask] = 1500
    data[1, mask] = 1000
    data[2, mask] = 130
    image = tmp_path / f"{name}.ome.tif"
    tifffile.imwrite(image, data, ome=True, photometric="minisblack",
                      metadata={"axes": "CZYX", "PhysicalSizeZ": 2.0, "PhysicalSizeY": 1.0, "PhysicalSizeX": 1.0})
    return image.name


def test_pipeline_writes_pairwise_contrasts_end_to_end(tmp_path):
    rows = []
    for condition, scale in (("Vehicle", 1.0), ("Treated", 1.6)):
        for replicate in range(3):
            name = f"{condition}_{replicate}"
            # Replicates differ slightly: identical replicates carry no replicate-level
            # variation, and the frozen method correctly declines to test them.
            image_name = _write_field(tmp_path, name, scale * (1.0 - 0.04 * replicate))
            rows.append({"sample_id": name, "condition": condition,
                        "biological_replicate": f"R{replicate}", "image_path": image_name})
    manifest = tmp_path / "manifest.csv"
    pd.DataFrame(rows).to_csv(manifest, index=False)

    cfg = load_config()
    cfg["channels"] = {"structure": 0, "calcein": 1, "pi": 2}
    cfg["report"]["save_meshes"] = False
    cfg["stats"]["min_replicates_per_condition"] = 3
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump(cfg))

    out = tmp_path / "result"
    result = analyze(manifest, out, config)
    assert result["status"] == "complete"

    pairwise = pd.read_csv(out / "pairwise_contrasts.csv")
    assert set(pairwise.feature) <= {"volume_um3", "sphericity"}
    assert (out / "stats_results.json").stat().st_size > 2
    volume_row = pairwise[pairwise.feature == "volume_um3"].iloc[0]
    assert volume_row.contrast == "Vehicle vs Treated" or volume_row.contrast == "Treated vs Vehicle"
    assert "df_denom" in pairwise.columns
    # Regression check: the computed statistical-test results must be surfaced
    # in the human-facing report, not only saved to disk (P2 audit finding).
    report_html = (out / "report.html").read_text(encoding="utf-8")
    assert "pairwise_contrasts.csv" in report_html
    assert "Cross-condition statistical testing" in report_html


def test_pipeline_skips_stats_cleanly_when_disabled(tmp_path):
    z, y, x = np.indices((20, 36, 36))
    mask = ((z - 10) * 2) ** 2 + (y - 18) ** 2 + (x - 18) ** 2 < 10 ** 2
    data = np.full((3, 20, 36, 36), 100, np.uint16)
    data[0, mask] = 1500
    data[1, mask] = 1000
    data[2, mask] = 130
    image = tmp_path / "stack.ome.tif"
    tifffile.imwrite(image, data, ome=True, photometric="minisblack",
                      metadata={"axes": "CZYX", "PhysicalSizeZ": 2.0, "PhysicalSizeY": 1.0, "PhysicalSizeX": 1.0})
    manifest = tmp_path / "manifest.csv"
    pd.DataFrame([{"sample_id": "F1", "condition": "Vehicle", "biological_replicate": "R1",
                  "image_path": image.name}]).to_csv(manifest, index=False)

    cfg = load_config()
    cfg["channels"] = {"structure": 0, "calcein": 1, "pi": 2}
    cfg["report"]["save_meshes"] = False
    cfg["stats"]["enabled"] = False
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump(cfg))

    out = tmp_path / "result"
    result = analyze(manifest, out, config)
    assert result["status"] == "complete"
    assert not (out / "pairwise_contrasts.csv").exists()
    assert not (out / "stats_results.json").exists()
