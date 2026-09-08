import numpy as np
import pandas as pd
import pytest
import scipy.stats as st
import tifffile
import yaml

from organoid_analysis.config import load_config
from organoid_analysis.workflows.organoid_measurement_workflow import analyze
from organoid_analysis.experimental_statistics.stats import _pairwise_contrasts, condition_pairwise_tests, fit_model


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


def test_lmm_pairwise_uses_small_cluster_t_reference_not_asymptotic_z():
    """Regression test for the P1 anti-conservative-p-value audit finding.

    statsmodels' MixedLM always reports Wald p-values against an asymptotic z
    reference, which is anti-conservative (overstates significance) with few
    replicate clusters. ``_pairwise_contrasts`` must instead use a t(G-1)
    reference (G = number of replicate clusters) for the LMM branch, matching
    the small-cluster correction statsmodels' own cluster-robust OLS fallback
    applies automatically via ``use_t=True``.
    """
    rng = np.random.default_rng(7)
    rows = []
    for condition in ("A", "B"):
        base = 1000.0 if condition == "A" else 1015.0
        for replicate in range(3):
            cluster_offset = rng.normal(0, 25)  # genuine between-replicate variance
            for _ in range(6):
                rows.append({"condition": condition, "grp": f"{condition}::R{replicate}",
                            "volume_um3": float(base + cluster_offset + rng.normal(0, 15))})
    d = pd.DataFrame(rows)
    d["condition"] = pd.Categorical(d.condition, categories=["A", "B"])
    fit, model_used = fit_model(d, "volume_um3", "grp")
    assert model_used.startswith("LMM")  # this test only exercises the LMM branch
    n_clusters = d["grp"].nunique()
    _, _, p_values, dof = _pairwise_contrasts(fit, ["A", "B"], n_clusters)
    assert dof == [n_clusters - 1]
    t_stat = float(np.ravel(fit.t_test(np.array([[0, 1]])).tvalue)[0])
    expected_t_p = float(2 * st.t.sf(abs(t_stat), n_clusters - 1))
    z_p = float(2 * st.norm.sf(abs(t_stat)))
    assert p_values[0] == pytest.approx(expected_t_p)
    assert p_values[0] > z_p  # t(G-1) reference must be more conservative than z


def test_ols_fallback_uses_cluster_robust_t_reference():
    """OLS-fallback p-values must come from use_t=True's t(G-1) reference."""
    rows = []
    for condition in ("A", "B"):
        for replicate in range(3):
            rows.append({"condition": condition, "grp": f"{condition}::R{replicate}", "value": 10.0})
    d = pd.DataFrame(rows)
    d["condition"] = pd.Categorical(d.condition, categories=["A", "B"])
    fit, model_used = fit_model(d, "value", "grp")
    assert model_used == "OLS(clustered SE)"
    assert fit.use_t is True


def test_ols_fallback_on_degenerate_random_effect():
    rows = []
    for condition in ("A", "B"):
        for replicate in range(3):
            rows.append({"condition": condition, "grp": f"{condition}::R{replicate}", "value": 10.0})
    d = pd.DataFrame(rows)
    d["condition"] = pd.Categorical(d.condition, categories=["A", "B"])
    _, model_used = fit_model(d, "value", "grp")
    assert model_used == "OLS(clustered SE)"


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
            image_name = _write_field(tmp_path, name, scale)
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
