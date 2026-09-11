"""Scientific validation gate, kept separate from the engineering `ci` gate.

    pixi run science-gate

`pixi run ci` (lint, tests, notebook check, typecheck) is the ENGINEERING
regression gate: it shows that the software behaves as its tests specify. It
does not show that measurements are accurate or that biological or statistical
claims are supported. This gate reports the scientific acceptance items and
exits 0 only when every item is PASS with existing evidence.

The analytical surface/volume items take their status from the frozen surface
V&V evidence. The gate first re-measures the reference sphere with the
production estimator; if the value no longer matches the evidence, the evidence
is stale and the gate fails until the V&V is repeated. No item can report PASS
without an evidence file, and this script never relaxes a criterion.
"""
from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from organoid_analysis.quantification.features import SURFACE_AREA_METHOD, geometry

ROOT = Path(__file__).resolve().parents[1]
SURFACE_EVIDENCE = ROOT / "docs/evidence/2026-09-11-measurement-vv"
SURFACE_CRITERION = 0.05  # SCIENTIFIC_SPEC section 9, analytical shapes
VOLUME_CRITERION = 0.01  # SCIENTIFIC_SPEC section 9, analytical shapes
# Recorded production area for the radius-18 um sphere at ZYX spacing (2,1,1) um
# (docs/evidence/2026-09-11/analytical-geometry.json).
REFERENCE_SPHERE_AREA_UM2 = 4549.28466796875
EVALUATED_METHOD_VERSION = "marching_cubes_binary_lewiner_v1"
VALID_STATUSES = {"PASS", "PARTIAL", "FAIL", "NOT ASSESSED", "INSUFFICIENT EVIDENCE", "NOT APPLICABLE"}


@dataclass(frozen=True)
class GateItem:
    item_id: str
    requirement: str
    criterion: str
    status: str
    basis: str
    evidence: tuple[str, ...]


def _reference_sphere_area() -> float:
    z, y, x = np.indices((25, 49, 49))
    mask = ((z - 12) * 2) ** 2 + (y - 24) ** 2 + (x - 24) ** 2 <= 18.0 ** 2
    return float(geometry(mask, (2.0, 1.0, 1.0))[0]["surface_area_um2"])


def _analytical_items() -> list[GateItem]:
    summary_path = SURFACE_EVIDENCE / "surface_vv_dev_summary.json"
    evidence = (str(summary_path.relative_to(ROOT)), str((SURFACE_EVIDENCE / "SURFACE_SELECTION.md").relative_to(ROOT)))
    stale = []
    if SURFACE_AREA_METHOD["method_version"] != EVALUATED_METHOD_VERSION:
        stale.append(f"production estimator is {SURFACE_AREA_METHOD['method_version']}, "
                     f"evidence covers {EVALUATED_METHOD_VERSION}")
    area = _reference_sphere_area()
    if not np.isclose(area, REFERENCE_SPHERE_AREA_UM2, rtol=1e-6, atol=0):
        stale.append(f"reference sphere area {area!r} no longer matches the evidence value")
    if not summary_path.is_file():
        stale.append("surface V&V summary missing")
    if stale:
        reason = "EVIDENCE STALE: " + "; ".join(stale) + ". Repeat the surface V&V."
        return [GateItem("SG-1", "Analytical surface area", "every analytical case |error| < 5%", "FAIL", reason, evidence),
                GateItem("SG-2", "Analytical voxel volume", "every analytical case |error| < 1%", "FAIL", reason, evidence)]
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    production = [r for r in summary["area"]["overall"] if r["estimator"].startswith("E0_")][0]
    area_status = "PASS" if production["max_abs"] < SURFACE_CRITERION else "FAIL"
    volume_rows = summary["volume"]["by_rho_stratum"]
    volume_status = "PASS" if all(r["max_abs"] < VOLUME_CRITERION for r in volume_rows) else "FAIL"
    volume_failures = ", ".join(f"{r['rho_stratum']} max {100 * r['max_abs']:.2f}%" for r in volume_rows
                                if r["max_abs"] >= VOLUME_CRITERION)
    return [
        GateItem("SG-1", "Analytical surface area (production estimator)", "every analytical case |error| < 5%",
                 area_status, f"{production['n']} cases; max |error| {100 * production['max_abs']:.2f}%, "
                 f"median signed {100 * production['median_signed']:+.2f}%", evidence),
        GateItem("SG-2", "Analytical voxel-count volume", "every analytical case |error| < 1%", volume_status,
                 f"fails for {volume_failures}" if volume_failures else "all strata < 1%", evidence),
    ]


# Items with no computable oracle in this repository. They can only become PASS
# when a qualified evidence file is added here together with its criterion.
STATIC_ITEMS = [
    GateItem("SG-3", "Segmentation accuracy against independent annotation (detection, mask, downstream "
             "measurement bias)", "predeclared per docs/evidence/2026-09-11-measurement-vv/SEGMENTATION_VALIDATION_PROTOCOL.md",
             "INSUFFICIENT EVIDENCE", "no qualified independent 3D annotation set exists", ()),
    GateItem("SG-4", "Calcein/PI viability assay validity", "orthogonal viability reference; control-based calibration "
             "evaluated on held-out batches", "NOT ASSESSED", "no orthogonal assay or held-out control data", ()),
    GateItem("SG-5", "Statistical and study-design validity", "predeclared design with an independent experimental "
             "unit; type-I error and CI coverage for the shipped models", "INSUFFICIENT EVIDENCE",
             "generic tools only; no study design or calibration evidence", ()),
    GateItem("SG-6", "Acquisition calibration and channel registration", "independently verified voxel size and "
             "registration for the intended instruments", "INSUFFICIENT EVIDENCE",
             "metadata is read and traced but not independently verified", ()),
]


def evaluate() -> list[GateItem]:
    items = _analytical_items() + STATIC_ITEMS
    for item in items:
        if item.status not in VALID_STATUSES:
            raise ValueError(f"{item.item_id}: unknown evidence status {item.status!r}")
        if item.status == "PASS" and not (item.evidence and all((ROOT / path).is_file() for path in item.evidence)):
            raise ValueError(f"{item.item_id}: PASS requires existing evidence files")
    return items


def main() -> int:
    items = evaluate()
    print("SCIENTIFIC VALIDATION GATE (separate from the engineering `ci` gate)\n")
    for item in items:
        print(f"[{item.status}] {item.item_id} {item.requirement}\n    criterion: {item.criterion}\n    basis: {item.basis}")
        for path in item.evidence:
            print(f"    evidence: {path}")
    failed = [item.item_id for item in items if item.status != "PASS"]
    print(f"\nGATE: {'PASS' if not failed else 'NOT PASSED'} ({len(items) - len(failed)}/{len(items)} items PASS)")
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
