"""Round-2 confirmation: execute confirm2 ONCE and adjudicate by the frozen rules.

Acceptance rules A1-A3 and validity rules V1-V3 are as written in
out/FREEZE_RECORD_V3.md section 6, before this grid was ever executed.
Nothing in this file may be changed after the freeze hashes are taken.
"""
import json
import pandas as pd

import surface_area_v3 as v3
from step8_dev2_run import run_grid   # identical code path to development

SMOOTH = ("sphere", "ellipsoid", "capsule", "torus")
CREASED = ("cylinder", "box")
A1, A2, A3 = 0.05, 0.01, 0.05          # area, volume, sphericity


def adjudicate(csv="out/vv_confirm2.csv"):
    d = pd.read_csv(csv)
    d["in_scope"] = d["shape"].isin(SMOOTH)
    d["gate"] = (d.rho_in >= v3.DOMAIN_RHO_IN_MIN) & \
                (d.anisotropy <= v3.DOMAIN_ANISO_MAX)
    S = d[d.in_scope & d.gate]

    v1 = (len(S) > 0) and set(S["shape"]) == set(SMOOTH)
    v2 = bool(d.area_v3.notna().all() and d.stencil_m.notna().all())

    a1 = float(S.area_rel_err.abs().max())
    a2 = float(S.volume_rel_err.abs().max())
    a3 = float(S.sphericity_rel_err.abs().max())
    passes = {"A1_area": a1 < A1, "A2_volume": a2 < A2, "A3_sphericity": a3 < A3}
    valid = {"V1_set_nonempty_all_classes": v1, "V2_all_measured": v2}
    verdict = "QUALIFIED" if (all(passes.values()) and all(valid.values())) \
        else ("VOID" if not all(valid.values()) else "FAIL / NOT QUALIFIED")

    oos = d[~(d.in_scope & d.gate)]
    rep = {
        "verdict": verdict,
        "criteria": {"A1_area": A1, "A2_volume": A2, "A3_sphericity": A3},
        "domain": {"rho_in_min": v3.DOMAIN_RHO_IN_MIN,
                   "aniso_max": v3.DOMAIN_ANISO_MAX,
                   "scope": v3.DOMAIN_SCOPE,
                   "not_qualified": v3.DOMAIN_NOT_QUALIFIED},
        "n_cases_executed": int(len(d)),
        "n_in_scope_and_gated": int(len(S)),
        "worst": {"area_rel_err": a1, "volume_rel_err": a2,
                  "sphericity_rel_err": a3},
        "rules": {**passes, **valid},
        "violations": {
            "A1": S.loc[S.area_rel_err.abs() >= A1, "case_id"].tolist(),
            "A2": S.loc[S.volume_rel_err.abs() >= A2, "case_id"].tolist(),
            "A3": S.loc[S.sphericity_rel_err.abs() >= A3, "case_id"].tolist()},
        "per_class_in_scope": {
            k: {"n": int((S["shape"] == k).sum()),
                "worst_area": float(S.loc[S["shape"] == k, "area_rel_err"].abs().max()),
                "worst_volume": float(S.loc[S["shape"] == k, "volume_rel_err"].abs().max()),
                "worst_sphericity": float(S.loc[S["shape"] == k, "sphericity_rel_err"].abs().max())}
            for k in SMOOTH},
        "reported_not_adjudicated": {
            "n_out_of_scope_or_below_gate": int(len(oos)),
            "creased_worst_area": float(
                d.loc[d["shape"].isin(CREASED), "area_rel_err"].abs().max()),
            "creased_in_gate_worst_area": float(
                d.loc[d["shape"].isin(CREASED) & d.gate, "area_rel_err"].abs().max()),
            "legacy_E0_worst_area_on_S": float(S.area_E0_rel_err.abs().max()),
            "legacy_E0_violations_of_A1_on_S": int((S.area_E0_rel_err.abs() >= A1).sum()),
            "legacy_E0_median_area_err_on_S": float(S.area_E0_rel_err.abs().median())},
    }
    json.dump(rep, open("out/confirm2_verdict.json", "w"), indent=1)
    return rep


if __name__ == "__main__":
    run_grid("confirm2")
    r = adjudicate()
    print(json.dumps({k: r[k] for k in
                      ("verdict", "n_cases_executed", "n_in_scope_and_gated",
                       "worst", "rules")}, indent=1))
