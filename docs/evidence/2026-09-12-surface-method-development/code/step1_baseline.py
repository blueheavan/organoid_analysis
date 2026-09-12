"""Step 1 -- re-measure the legacy baseline E0 and lock its identity.

E0 is measured *through the production call path*
``organoid_analysis.quantification.features.geometry`` so the baseline cannot
drift from what the pipeline actually exports.  The reproduced per-case areas
are compared against the frozen record CSV.
"""
from __future__ import annotations

import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import skimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (CRITERION, FROZEN_RECORD, REPO, VOLUME_CRITERION, case_id,  # noqa: E402
                    cases, oracle_self_check, rasterize, rho_stratum, truth)

sys.path.insert(0, str(REPO / "src"))
from organoid_analysis.quantification.features import (SURFACE_AREA_METHOD,  # noqa: E402
                                                       geometry)

OUT = Path(__file__).resolve().parent / "out"
OUT.mkdir(exist_ok=True)


def e0_production(mask: np.ndarray, spacing) -> float:
    """The production estimator, called exactly as the pipeline calls it."""
    return float(geometry(mask, spacing, fill_holes=False)[0]["surface_area_um2"])


def main() -> None:
    oracle = oracle_self_check()

    records = []
    for idx, (kind, p, spacing, variant, rot, offset, min_axis) in enumerate(cases("dev")):
        mask = rasterize(kind, p, spacing, rot, offset)
        a_true, v_true = truth(kind, p)
        a_est = e0_production(mask, spacing)
        v_meas = float(mask.sum() * np.prod(spacing))
        rho = min_axis / max(spacing)
        records.append({
            "case_id": case_id("dev", idx, kind, p, spacing, variant),
            "shape": kind, "rho": rho, "rho_stratum": rho_stratum(rho),
            "spacing": "x".join(f"{s:g}" for s in spacing),
            "area_true_um2": a_true, "area_est_um2": a_est, "area_rel_err": a_est / a_true - 1,
            "volume_true_um3": v_true, "volume_voxel_um3": v_meas,
            "volume_rel_err": v_meas / v_true - 1,
        })
    df = pd.DataFrame.from_records(records)

    # ---- compare against the frozen record -------------------------------
    frozen = pd.read_csv(FROZEN_RECORD / "surface_vv_dev.csv")
    frozen = frozen[frozen.estimator == "E0_marching_cubes_binary_v1"][
        ["case_id", "area_est_um2", "area_rel_err", "volume_voxel_um3"]]
    merged = df.merge(frozen, on="case_id", suffixes=("", "_frozen"), how="outer", indicator=True)
    assert (merged["_merge"] == "both").all(), merged.loc[merged["_merge"] != "both", "case_id"].tolist()
    area_dev = float(np.max(np.abs(merged.area_est_um2 / merged.area_est_um2_frozen - 1)))
    vol_dev = float(np.max(np.abs(merged.volume_voxel_um3 / merged.volume_voxel_um3_frozen - 1)))
    merged.drop(columns="_merge").to_csv(OUT / "baseline_reproduction.csv", index=False,
                                         float_format="%.12g")

    src = (REPO / "src/organoid_analysis/quantification/features.py").read_bytes()
    baseline = {
        "estimator_id": "E0_marching_cubes_binary_v1",
        "role": "legacy baseline -- preserved, not modified, not deleted",
        "method_identity": dict(SURFACE_AREA_METHOD),
        "production_call_path": ("organoid_analysis.quantification.features.geometry"
                                 " -> surface_mesh -> skimage.measure.mesh_surface_area"),
        "production_source_sha256": hashlib.sha256(src).hexdigest(),
        "volume_semantics": "envelope voxel count x prod(spacing); unchanged by this study",
        "sphericity_semantics": "pi^(1/3) (6V)^(2/3) / A; unclipped; flag above 1.05",
        "acceptance_criteria_applied": {
            "surface_abs_rel_error_per_case": CRITERION,
            "volume_abs_rel_error_per_case": VOLUME_CRITERION,
        },
        "current_validation_status": "FAIL (surface); criterion NOT relaxed by this study",
        "fail_evidence_dev_grid": {
            "n_cases": int(len(df)), "n_nonestimable": int(df.area_rel_err.isna().sum()),
            "mean_signed": float(df.area_rel_err.mean()),
            "median_signed": float(df.area_rel_err.median()),
            "mean_abs": float(df.area_rel_err.abs().mean()),
            "median_abs": float(df.area_rel_err.abs().median()),
            "max_abs": float(df.area_rel_err.abs().max()),
            "worst_case": df.loc[df.area_rel_err.abs().idxmax(), "case_id"],
            "n_cases_failing_5pct": int((df.area_rel_err.abs() >= CRITERION).sum()),
            "by_shape_mean_signed": {k: float(g.area_rel_err.mean())
                                     for k, g in df.groupby("shape")},
            "by_rho_stratum_mean_signed": {k: float(g.area_rel_err.mean())
                                           for k, g in df.groupby("rho_stratum")},
        },
        "reproduction_check_vs_frozen_record": {
            "frozen_csv": str((FROZEN_RECORD / "surface_vv_dev.csv").relative_to(REPO)),
            "n_cases_matched": int(len(merged)),
            "max_abs_rel_deviation_area": area_dev,
            "max_abs_rel_deviation_volume": vol_dev,
            "tolerance": 1e-9,
            "status": "REPRODUCED" if max(area_dev, vol_dev) < 1e-9 else "DEVIATION",
        },
        "oracle_self_check_max_rel_diff": max(c["rel_diff"] for c in oracle),
        "environment": {"python": sys.version.split()[0], "platform": platform.platform(),
                        "numpy": np.__version__, "scipy": scipy.__version__,
                        "scikit-image": skimage.__version__, "pandas": pd.__version__},
    }
    (OUT / "E0_legacy_baseline.json").write_text(json.dumps(baseline, indent=2))
    print(json.dumps({"reproduction": baseline["reproduction_check_vs_frozen_record"],
                      "fail_evidence": baseline["fail_evidence_dev_grid"]}, indent=2))


if __name__ == "__main__":
    main()
