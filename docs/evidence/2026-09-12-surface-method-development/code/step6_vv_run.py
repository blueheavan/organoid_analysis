"""Frozen V&V harness -- runs the SAME code on the dev grid and the
confirmation grid, applies the SAME pre-declared acceptance rules, and prints
a verdict it has no way to soften.

    python step6_vv_run.py --grid dev
    python step6_vv_run.py --grid confirm     # run exactly once, after freeze

Acceptance rules (frozen; see FREEZE_RECORD.md). For every case inside the
qualified applicability domain:

    A1  |area error|        <  5 %      (the project §9 criterion, unchanged)
    A2  |volume error|      <  1 %      (existing §9 criterion, re-checked)
    A3  |sphericity error|  <  5 %      (required because area changed)

Out-of-domain cases are measured and reported, and take no part in the
verdict. The verdict is PASS only if A1, A2 and A3 hold for EVERY in-domain
case; a single violation is FAIL. Nothing here aggregates, averages or
trims.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

import common
from surface_area_v2 import METHOD_NAME, measure, sphericity

sys.path.insert(0, str(common.REPO / "src"))
from organoid_analysis.quantification.features import geometry  # noqa: E402

OUT = Path(__file__).resolve().parent / "out"
AREA_CRITERION = 0.05
VOLUME_CRITERION = 0.01
SPHERICITY_CRITERION = 0.05
FIELDS = ["case_id", "grid", "shape", "spacing", "anisotropy", "rho_frozen", "rho_eff",
          "stencil_m", "n_directions", "apriori_bound", "in_domain", "domain_reasons",
          "area_true", "area_v2", "area_rel_err", "area_E0", "area_E0_rel_err",
          "volume_true", "volume_est", "volume_rel_err",
          "sphericity_true", "sphericity_v2", "sphericity_rel_err"]


def run(grid: str) -> dict:
    rows = []
    for case in common.case_records(grid):
        geo = (case["shape"], case["params"], case["spacing"], case["rot"], case["offset"])
        mask = common.rasterize(*geo)
        a_true, v_true = common.truth(case["shape"], case["params"])
        r = measure(mask, case["spacing"])
        g_e0, _ = geometry(mask, case["spacing"], fill_holes=False)
        s_true = sphericity(v_true, a_true)
        s_v2 = sphericity(r["volume"], r["surface_area"])
        rows.append({
            "case_id": case["case_id"], "grid": grid, "shape": case["shape"],
            "spacing": "x".join(str(s) for s in case["spacing"]),
            "anisotropy": round(r["anisotropy"], 4), "rho_frozen": case["rho"],
            "rho_eff": r["rho_eff"], "stencil_m": r["stencil_m"],
            "n_directions": r["n_directions"], "apriori_bound": r["apriori_bound"],
            "in_domain": r["in_domain"], "domain_reasons": ";".join(r["domain_reasons"]),
            "area_true": a_true, "area_v2": r["surface_area"],
            "area_rel_err": r["surface_area"] / a_true - 1,
            "area_E0": float(g_e0["surface_area_um2"]),
            "area_E0_rel_err": float(g_e0["surface_area_um2"]) / a_true - 1,
            "volume_true": v_true, "volume_est": r["volume"],
            "volume_rel_err": r["volume"] / v_true - 1,
            "sphericity_true": s_true, "sphericity_v2": s_v2,
            "sphericity_rel_err": s_v2 / s_true - 1})

    with open(OUT / f"vv_{grid}.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=FIELDS)
        wr.writeheader()
        wr.writerows(rows)

    ind = [r for r in rows if r["in_domain"]]
    out = [r for r in rows if not r["in_domain"]]

    def worst(rs, key):
        return max(rs, key=lambda r: abs(r[key])) if rs else None

    checks = {}
    for name, key, crit in (("A1_area", "area_rel_err", AREA_CRITERION),
                            ("A2_volume", "volume_rel_err", VOLUME_CRITERION),
                            ("A3_sphericity", "sphericity_rel_err", SPHERICITY_CRITERION)):
        viol = [r for r in ind if abs(r[key]) >= crit]
        w = worst(ind, key)
        checks[name] = {"criterion": crit, "n_violations": len(viol),
                        "status": "PASS" if not viol else "FAIL",
                        "worst_case": w["case_id"] if w else None,
                        "worst_abs": abs(w[key]) if w else None,
                        "violating_cases": [r["case_id"] for r in viol][:20]}

    e0_viol = [r for r in ind if abs(r["area_E0_rel_err"]) >= AREA_CRITERION]
    verdict = "PASS" if all(c["status"] == "PASS" for c in checks.values()) else "FAIL"
    summary = {
        "grid": grid, "method": METHOD_NAME, "n_cases": len(rows),
        "n_in_domain": len(ind), "n_out_of_domain": len(out),
        "checks": checks, "verdict": verdict,
        "legacy_E0_in_domain_violations": len(e0_viol),
        "legacy_E0_worst_abs": abs(worst(ind, "area_E0_rel_err")["area_E0_rel_err"]) if ind else None,
        "out_of_domain_worst_area_abs": abs(worst(out, "area_rel_err")["area_rel_err"]) if out else None,
        "out_of_domain_area_violations": sum(
            1 for r in out if abs(r["area_rel_err"]) >= AREA_CRITERION),
    }
    (OUT / f"vv_{grid}_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", required=True, choices=["dev", "confirm"])
    run(ap.parse_args().grid)
