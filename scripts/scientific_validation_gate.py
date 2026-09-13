"""Scientific validation gate, kept separate from the engineering regression gate.

    pixi run science-gate

`pixi run regression` (lint, tests, notebook check) is the ENGINEERING
regression gate: it shows that the software behaves as its tests specify. It
does not show that measurements are accurate or that biological or statistical
claims are supported. This gate reports the scientific acceptance items and
exits 0 only when every item is PASS backed by a verified, canonical
validation record.

SG-1/SG-2 take their status from the analytical-geometry validation record
named in docs/evidence/analytical_geometry_current_record.json, and only after
that record verifies (analytical_geometry_evidence.verify_record): every file
of the declared dependency scope, the frozen plan and harness, pixi.lock,
package versions, raw and derived artifacts, estimator identity and criteria
must be unchanged; the recorded hashes must be the content of the named source
commit; the statuses must re-derive from the raw results; and the production
estimator is re-executed on every phantom. Any discrepancy rejects the evidence
and both items FAIL until the V&V is repeated (`pixi run validation-record`).
No item can report PASS without such a record, and this script never relaxes a
criterion. See docs/VALIDATION_RECORDS.md.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

from organoid_analysis.validation import analytical_geometry_evidence as analytical
from organoid_analysis.validation.evidence_manifest import EvidenceError

ROOT = Path(__file__).resolve().parents[1]
VALID_STATUSES = {"PASS", "PARTIAL", "FAIL", "NOT ASSESSED", "INSUFFICIENT EVIDENCE", "NOT APPLICABLE"}
SURFACE_CRITERION_TEXT = "every analytical case |error| < 5%"
VOLUME_CRITERION_TEXT = "every analytical case |error| < 1%"


@dataclass(frozen=True)
class GateItem:
    item_id: str
    requirement: str
    criterion: str
    status: str
    basis: str
    evidence: tuple[str, ...]
    # Repository-relative evidence manifest backing this item; None when the
    # item has no verified validation record (and therefore cannot PASS).
    record: str | None = None
    counts_toward_gate: bool = True


Reports = dict[str, analytical.VerificationReport]


def _verify(root: Path, manifest: str, reports: Reports, check_git: bool, reproduce: bool) -> analytical.VerificationReport:
    if manifest not in reports:
        reports[manifest] = analytical.verify_record(root, manifest, check_git=check_git, reproduce=reproduce)
    return reports[manifest]


def _percent(value: float | None, signed: bool = False) -> str:
    return "non-estimable" if value is None else f"{100 * value:{'+' if signed else ''}.2f}%"


def _analytical_items(root: Path, reports: Reports, check_git: bool, reproduce: bool) -> list[GateItem]:
    def rejected(reason: str, evidence: tuple[str, ...] = ()) -> list[GateItem]:
        return [GateItem("SG-1a", "Surface analytical qualification", SURFACE_CRITERION_TEXT,
                         "FAIL", reason, evidence),
                GateItem("SG-1b", "Surface unrestricted characterization", "descriptive only",
                         "NOT APPLICABLE", reason, evidence, counts_toward_gate=False),
                GateItem("SG-2a", "Volume analytical qualification", VOLUME_CRITERION_TEXT,
                         "FAIL", reason, evidence),
                GateItem("SG-2b", "Volume unrestricted characterization", "descriptive only",
                         "NOT APPLICABLE", reason, evidence, counts_toward_gate=False)]

    try:
        manifest_path = analytical.current_manifest_relpath(root)
    except EvidenceError as error:
        return rejected(f"EVIDENCE REJECTED: {error}. Run the analytical V&V (pixi run validation-record).")
    run_dir = str(PurePosixPath(manifest_path).parent)
    evidence = (manifest_path, f"{run_dir}/{analytical.DERIVED_NAMES[1]}", analytical.PLAN_PATH)
    report = _verify(root, manifest_path, reports, check_git, reproduce)
    if not report.valid or report.manifest is None:
        shown = "; ".join(report.problems[:4]) + (f"; ... {len(report.problems) - 4} more" if len(report.problems) > 4 else "")
        return rejected(f"EVIDENCE REJECTED ({len(report.problems)} problems): {shown}. "
                        "Repeat the analytical V&V (pixi run validation-record).", evidence)

    manifest = report.manifest
    record = (f" [record {manifest['validation_run_id']}, {manifest['record_class']}, "
              f"source {manifest['source']['validated_source_commit'][:12]}]")
    sg1, sg2 = manifest["results"]["SG-1"], manifest["results"]["SG-2"]
    area_basis = (f"{sg1['n_cases']} cases; max |error| {_percent(sg1['max_abs_rel_error'])} ({sg1['worst_case']}), "
                  f"median signed {_percent(sg1['median_signed_rel_error'], signed=True)}, "
                  f"non-estimable {sg1['n_nonestimable']}")
    failing = [f"{row['rho_stratum']} max {_percent(row['max_abs_rel_error'])}"
               for row in sg2["by_rho_stratum"] if row["status"] != "PASS"]
    volume_basis = f"fails for {', '.join(failing)}" if failing else "all strata < 1%"
    return [
        GateItem("SG-1a", "Surface analytical qualification", SURFACE_CRITERION_TEXT,
                 "INSUFFICIENT EVIDENCE",
                 "no domain-restricted canonical qualification record; unrestricted record is characterization only",
                 evidence),
        GateItem("SG-1b", "Surface unrestricted characterization", "descriptive only",
                 "NOT APPLICABLE", area_basis + record, evidence, manifest_path, counts_toward_gate=False),
        GateItem("SG-2a", "Volume analytical qualification", VOLUME_CRITERION_TEXT,
                 "INSUFFICIENT EVIDENCE",
                 "no volume-specific domain or canonical qualification record; unrestricted record is characterization only",
                 evidence),
        GateItem("SG-2b", "Volume unrestricted characterization", "descriptive only",
                 "NOT APPLICABLE", volume_basis + record, evidence, manifest_path, counts_toward_gate=False),
    ]


# Items with no computable oracle in this repository. They carry no validation
# record, so evaluate() refuses PASS for them until one exists.
STATIC_ITEMS = [
    GateItem("SG-3A", "Brightfield segmentation accuracy against independent annotation (detection, mask, downstream "
             "measurement bias)", "predeclared per docs/evidence/2026-09-11-measurement-vv/SEGMENTATION_VALIDATION_PROTOCOL.md",
             "INSUFFICIENT EVIDENCE", "no qualified independent 3D annotation set exists", ()),
    GateItem("SG-3B", "Membrane-fluorescence segmentation accuracy against independent annotation (detection, mask, downstream "
             "measurement bias)", "predeclared per docs/evidence/2026-09-11-measurement-vv/SEGMENTATION_VALIDATION_PROTOCOL.md",
             "INSUFFICIENT EVIDENCE", "no qualified independent 3D annotation set exists", ()),
    GateItem("SG-4A", "Calcein/PI four-state analytical classification", "known-rule state and refusal contract with versioned denominator",
             "PARTIAL", "state logic and classifiable-denominator contract are unit-tested; biological identity remains unvalidated", ()),
    GateItem("SG-4B", "Calcein/PI biological validity", "orthogonal object-level reference; control-based calibration "
             "evaluated on held-out batches", "NOT ASSESSED", "no orthogonal assay or held-out control data", ()),
    GateItem("SG-5", "Statistical and study-design validity", "predeclared design with an independent experimental "
             "unit; type-I error and CI coverage for the shipped models", "INSUFFICIENT EVIDENCE",
             "generic tools only; no study design or calibration evidence", ()),
    GateItem("SG-6A", "Relative spacing-ratio verification", "independently verified X:Y and Z:XY ratios for domain stratification",
             "NOT ASSESSED", "no ratio-standard study", ()),
    GateItem("SG-6B", "Absolute acquisition calibration and channel registration", "independently verified voxel size and "
             "registration for the intended instruments", "INSUFFICIENT EVIDENCE",
             "metadata is read and traced but not independently verified", ()),
]


def _require_qualifying_record(root: Path, item: GateItem, reports: Reports, check_git: bool, reproduce: bool) -> None:
    if item.record is None:
        raise ValueError(f"{item.item_id}: PASS requires a verified validation record")
    if not item.evidence or not all((root / path).is_file() for path in item.evidence):
        raise ValueError(f"{item.item_id}: PASS requires existing evidence files")
    report = _verify(root, item.record, reports, check_git, reproduce)
    if not report.valid or report.manifest is None:
        raise ValueError(f"{item.item_id}: PASS rejected; its validation record does not verify: {report.problems}")
    if report.manifest.get("record_class") != "canonical":
        raise ValueError(f"{item.item_id}: PASS requires a canonical record from a clean source commit, "
                         f"not {report.manifest.get('record_class')!r}")
    result = report.manifest.get("results", {}).get(item.item_id)
    if not isinstance(result, dict) or result.get("status") != "PASS":
        raise ValueError(f"{item.item_id}: PASS is not what its validation record reports")


def evaluate(root: Path = ROOT, *, check_git: bool = True, reproduce: bool = True) -> list[GateItem]:
    reports: Reports = {}
    items = _analytical_items(root, reports, check_git, reproduce) + STATIC_ITEMS
    for item in items:
        if item.status not in VALID_STATUSES:
            raise ValueError(f"{item.item_id}: unknown evidence status {item.status!r}")
        if item.status == "PASS":
            _require_qualifying_record(root, item, reports, check_git, reproduce)
    return items


def main() -> int:
    items = evaluate()
    print("SCIENTIFIC VALIDATION GATE (separate from the engineering regression gate)\n")
    for item in items:
        print(f"[{item.status}] {item.item_id} {item.requirement}\n    criterion: {item.criterion}\n    basis: {item.basis}")
        for path in item.evidence:
            print(f"    evidence: {path}")
    counted = [item for item in items if item.counts_toward_gate]
    failed = [item.item_id for item in counted if item.status != "PASS"]
    passed = len(counted) - len(failed)
    print(f"\nGATE: {'PASS' if not failed else 'NOT PASSED'} ({passed}/{len(counted)} qualification items PASS; "
          f"{len(items) - len(counted)} characterization items reported)")
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
