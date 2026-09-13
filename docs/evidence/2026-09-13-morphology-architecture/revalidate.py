"""Reexecute frozen geometry evidence with current production code, without resealing it.

Run: pixi run python docs/evidence/2026-09-13-morphology-architecture/revalidate.py
This is regression evidence for a dirty working tree, not a canonical release record
or a new independent confirmation study. Historical records are never overwritten.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from organoid_analysis.quantification import surface_crofton as sc
from organoid_analysis.quantification.features import SURFACE_AREA_METHOD
from organoid_analysis.validation.analytical_geometry_evidence import reproduce_production

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
HISTORY = ROOT / "docs/evidence/2026-09-12-surface-method-development"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    old = ROOT / "docs/evidence/2026-09-12-surface-crofton-v3"
    problems, reproduction = reproduce_production(old / "surface_vv.py", old / "surface_vv_dev.csv")
    print(json.dumps(reproduction), flush=True)
    sys.path.insert(0, str(HISTORY / "code"))
    import common2 as phantoms

    phantoms.oracle_self_check2()
    with (HISTORY / "vv_confirm2.csv").open() as handle:
        recorded = {r["case_id"]: r for r in csv.DictReader(handle)
                    if r["shape"] in {"sphere", "ellipsoid", "capsule", "torus"}
                    and float(r["rho_in"]) >= 10.0 and float(r["anisotropy"]) <= 4.0}
    assert len(recorded) == 96, "Frozen confirmation subset changed"
    output = []
    for r in phantoms.case_records2("confirm2"):
        if r["case_id"] not in recorded:
            continue
        mask = phantoms.rasterize2(r["shape"], r["params"], r["spacing"], r["rot"], r["offset"])
        area_true, volume_true = phantoms.truth2(r["shape"], r["params"])
        current = sc.measure(mask, r["spacing"])
        stored = recorded[r["case_id"]]
        for key, column in (("surface_area", "area_v3"), ("volume", "volume_est"),
                            ("rho_in", "rho_in"), ("stencil_m", "stencil_m")):
            if not np.isclose(current[key], float(stored[column]), rtol=1e-9, atol=0):
                problems.append(f"{r['case_id']}: {key} differs from the frozen result")
        assert sc.in_domain(current)[0]
        output.append({"case_id": r["case_id"], "shape": r["shape"],
                       "area_rel_err": current["surface_area"] / area_true - 1,
                       "volume_rel_err": current["volume"] / volume_true - 1,
                       "sphericity_rel_err": sc.sphericity(current["volume"], current["surface_area"])
                           / sc.sphericity(volume_true, area_true) - 1,
                       "weights_origin": current["weights_origin"]})
        if len(output) % 12 == 0:
            print(f"Confirmation subset reexecuted: {len(output)}/96", flush=True)
    worst = {key: max(abs(row[key]) for row in output)
             for key in ("area_rel_err", "volume_rel_err", "sphericity_rel_err")}
    for key, threshold in (("area_rel_err", .05), ("volume_rel_err", .01), ("sphericity_rel_err", .05)):
        if worst[key] >= threshold:
            problems.append(f"{key}: {worst[key]} exceeds the unchanged criterion {threshold}")
    with (OUT / "confirmation-reexecution.csv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    source_paths = sorted((ROOT / "src/organoid_analysis").rglob("*.py"))
    reference_paths = [HISTORY / "vv_confirm2.csv", HISTORY / "code/common.py", HISTORY / "code/common2.py",
                       old / "surface_vv.py", old / "surface_vv_dev.csv", ROOT / "pixi.lock"]
    report = {"scope": "numerical regression; no segmentation or biological validation",
              "method": SURFACE_AREA_METHOD, "production_grid": reproduction,
              "confirmation_cases": len(output), "confirmation_worst_absolute_relative_errors": worst,
              "problems": problems, "status": "PASS" if not problems else "FAIL",
              "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
              "source_state": "dirty working tree; not a canonical validation record",
              "source_sha256": {str(p.relative_to(ROOT)): digest(p) for p in source_paths},
              "reference_sha256": {str(p.relative_to(ROOT)): digest(p) for p in reference_paths}}
    (OUT / "revalidation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("status", "confirmation_cases", "confirmation_worst_absolute_relative_errors", "problems")}), flush=True)
    return bool(problems)


if __name__ == "__main__":
    raise SystemExit(main())
