"""Scientific validation gate, kept separate from the engineering regression gate.

    pixi run science-gate

`pixi run regression` (lint, tests, notebook check) is the ENGINEERING
regression gate: it shows that the software behaves as its tests specify. It
does not show that measurements are accurate or that biological or statistical
claims are supported. This gate reports the scientific acceptance items and
exits 0 only when every item is PASS backed by a verified, canonical
validation record.

Qualification items read item-keyed records (owner decision D-14):
SG-1a/SG-2a from the analytical-qualification record
(docs/evidence/analytical_qualification_current_record.json) and SG-4A from the
viability-rule record (docs/evidence/viability_rule_current_record.json). Each
is PASS only when its record verifies (hashes, frozen protocol, Git provenance,
re-derivation from raw results and re-execution of production code), is
canonical, and reports PASS for that item.

The characterization items SG-1b/SG-2b take their numbers from the analytical-geometry validation record
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
from typing import Any

from organoid_analysis.validation import analytical_geometry_evidence as analytical
from organoid_analysis.validation import analytical_qualification_evidence as qualification
from organoid_analysis.validation import statistical_qualification_evidence as statistical
from organoid_analysis.validation import viability_rule_evidence as viability
from organoid_analysis.validation.evidence_manifest import EvidenceError, load_manifest
from organoid_analysis.validation.record_contract import current_manifest_relpath

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


Reports = dict[str, "analytical.VerificationReport | qualification.VerificationReport"]


def _contract_of(root: Path, manifest: str) -> str | None:
    try:
        return load_manifest(root / manifest).get("contract_id")
    except EvidenceError:
        return None


def _verify(root: Path, manifest: str, reports: Reports, check_git: bool, reproduce: bool) -> Any:
    """Verify a record with the verifier of the contract it declares (never a looser one)."""
    if manifest not in reports:
        contract = _contract_of(root, manifest)
        verifier = {qualification.CONTRACT_ID: qualification.verify_record,
                    viability.CONTRACT_ID: viability.verify_record,
                    statistical.CONTRACT_ID: statistical.verify_record}.get(contract or "", analytical.verify_record)
        reports[manifest] = verifier(root, manifest, check_git=check_git, reproduce=reproduce)
    return reports[manifest]


def _percent(value: float | None, signed: bool = False) -> str:
    return "non-estimable" if value is None else f"{100 * value:{'+' if signed else ''}.2f}%"


def _characterization_items(root: Path, reports: Reports, check_git: bool, reproduce: bool) -> list[GateItem]:
    """SG-1b/SG-2b: unrestricted characterization; never PASS/FAIL and never counted."""
    def rejected(reason: str, evidence: tuple[str, ...] = ()) -> list[GateItem]:
        return [GateItem("SG-1b", "Surface unrestricted characterization", "descriptive only",
                         "NOT APPLICABLE", reason, evidence, counts_toward_gate=False),
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
        return rejected(f"EVIDENCE REJECTED ({len(report.problems)} problems): {_shown(report.problems)}. "
                        "Repeat the analytical V&V (pixi run validation-record).", evidence)

    manifest = report.manifest
    record = (f" [record {manifest['validation_run_id']}, {manifest['record_class']}, "
              f"source {manifest['source']['validated_source_commit'][:12]}]")
    sg1, sg2 = manifest["results"]["SG-1"], manifest["results"]["SG-2"]
    area_basis = (f"{sg1['n_cases']} cases; max |error| {_percent(sg1['max_abs_rel_error'])} ({sg1['worst_case']}), "
                  f"median signed {_percent(sg1['median_signed_rel_error'], signed=True)}, "
                  f"non-estimable {sg1['n_nonestimable']}")
    strata = [f"{row['rho_stratum']} max {_percent(row['max_abs_rel_error'])}" for row in sg2["by_rho_stratum"]]
    volume_basis = "per rho_in stratum: " + ", ".join(strata)
    return [
        GateItem("SG-1b", "Surface unrestricted characterization", "descriptive only",
                 "NOT APPLICABLE", area_basis + record, evidence, manifest_path, counts_toward_gate=False),
        GateItem("SG-2b", "Volume unrestricted characterization", "descriptive only",
                 "NOT APPLICABLE", volume_basis + record, evidence, manifest_path, counts_toward_gate=False),
    ]


def _shown(problems: tuple[str, ...]) -> str:
    return "; ".join(problems[:4]) + (f"; ... {len(problems) - 4} more" if len(problems) > 4 else "")


@dataclass(frozen=True)
class QualifiedItem:
    item_id: str
    requirement: str
    criterion: str
    no_record_status: str
    no_record_basis: str


def _qualification_basis(item_id: str, result: dict[str, Any]) -> str:
    if item_id == "SG-1a":
        parts = [f"{grid}: n={row['n']}, max |error| {_percent(row['max_abs_rel_error'])} ({row['worst_case']})"
                 for grid, row in result["by_grid"].items()]
        return "; ".join(parts) + f"; domain rho_in>={result['domain']['rho_in_min']}, aniso<={result['domain']['anisotropy_max']}, smooth"
    if item_id == "SG-2a":
        row = result["confirm3"]
        return (f"confirm3: n={row['n']}, max |error| {_percent(row['max_abs_rel_error'])} ({row['worst_case']}); "
                f"domain rho_vol>={result['domain']['rho_vol_min']}, rho_in>={result['domain']['rho_in_min']}, "
                f"aniso<={result['domain']['anisotropy_max']}, smooth")
    if item_id == "SG-5":
        contrast, interval = result["contrast"], result["condition_interval"]
        env = result["envelope"]
        return (f"{contrast['n_cells']} contrast cells: type-I {contrast['type1_range'][0]:.4f}-"
                f"{contrast['type1_range'][1]:.4f}, coverage {contrast['coverage_range'][0]:.4f}-"
                f"{contrast['coverage_range'][1]:.4f}; {interval['n_cells']} condition-interval cells: coverage "
                f"{interval['coverage_range'][0]:.4f}-{interval['coverage_range'][1]:.4f}; reference verification "
                f"{'passed' if result['reference_verification']['passed'] else 'FAILED'}; envelope "
                f"{env['min_replicates']}-{env['max_replicates']} replicates/condition, <= {env['max_conditions']} "
                f"conditions")
    return (f"{result['n_cases']} known-rule cases, {result['n_comparisons']} comparisons, "
            f"{result['n_mismatched_cases']} mismatched; denominator {result['fraction_definition']}; "
            f"{result['claim_boundary']}")


def _qualified_items(root: Path, module: Any, definitions: list[QualifiedItem], reports: Reports,
                     check_git: bool, reproduce: bool) -> list[GateItem]:
    """Items read from a qualification record: PASS only if it verifies, is canonical and reports PASS."""
    try:
        manifest_path = current_manifest_relpath(root, module.SPEC)
    except EvidenceError as error:
        return [GateItem(d.item_id, d.requirement, d.criterion, d.no_record_status, f"{d.no_record_basis} ({error})", ())
                for d in definitions]
    report = _verify(root, manifest_path, reports, check_git, reproduce)
    if not report.valid or report.manifest is None:
        reason = f"EVIDENCE REJECTED ({len(report.problems)} problems): {_shown(report.problems)}"
        return [GateItem(d.item_id, d.requirement, d.criterion, "FAIL", reason, (manifest_path,)) for d in definitions]
    manifest = report.manifest
    items = []
    for d in definitions:
        result = manifest["results"][d.item_id]
        record = (f" [record {manifest['validation_run_id']}, {manifest['record_class']}, "
                  f"source {manifest['source']['validated_source_commit'][:12]}]")
        basis = _qualification_basis(d.item_id, result) + record
        if result["status"] == "PASS" and manifest["record_class"] == "canonical":
            items.append(GateItem(d.item_id, d.requirement, d.criterion, "PASS", basis, (manifest_path,), manifest_path))
        elif result["status"] == "PASS":
            items.append(GateItem(d.item_id, d.requirement, d.criterion, "INSUFFICIENT EVIDENCE",
                                  "record is not canonical; " + basis, (manifest_path,)))
        elif result["status"] == "VOID":
            items.append(GateItem(d.item_id, d.requirement, d.criterion, "INSUFFICIENT EVIDENCE",
                                  f"confirmation run VOID: {result.get('validity_problems')}; " + basis, (manifest_path,)))
        else:
            items.append(GateItem(d.item_id, d.requirement, d.criterion, "FAIL", basis, (manifest_path,)))
    return items


ANALYTICAL_QUALIFICATION = [
    QualifiedItem("SG-1a", "Surface analytical qualification (declared domain)", SURFACE_CRITERION_TEXT,
                  "INSUFFICIENT EVIDENCE", "no domain-restricted canonical qualification record"),
    QualifiedItem("SG-2a", "Volume analytical qualification (declared volume domain)", VOLUME_CRITERION_TEXT,
                  "INSUFFICIENT EVIDENCE", "no volume-domain canonical qualification record"),
]
STATISTICAL_QUALIFICATION = [
    QualifiedItem("SG-5", "Statistical and study-design validity",
                  "predeclared design with an independent experimental unit; type-I error and CI coverage for the "
                  "shipped models", "INSUFFICIENT EVIDENCE",
                  "no canonical statistical qualification record"),
]
VIABILITY_QUALIFICATION = [
    QualifiedItem("SG-4A", "Calcein/PI four-state analytical classification",
                  "known-rule state and refusal contract with versioned denominator", "PARTIAL",
                  "state logic and denominator contract are unit-tested; no canonical analytical record"),
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
    GateItem("SG-4B", "Calcein/PI biological validity", "orthogonal object-level reference; control-based calibration "
             "evaluated on held-out batches", "NOT ASSESSED", "no orthogonal assay or held-out control data", ()),
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
    items = (_qualified_items(root, qualification, ANALYTICAL_QUALIFICATION, reports, check_git, reproduce)
             + _characterization_items(root, reports, check_git, reproduce)
             + _qualified_items(root, viability, VIABILITY_QUALIFICATION, reports, check_git, reproduce)
             + _qualified_items(root, statistical, STATISTICAL_QUALIFICATION, reports, check_git, reproduce)
             + STATIC_ITEMS)
    order = ["SG-1a", "SG-1b", "SG-2a", "SG-2b", "SG-3A", "SG-3B", "SG-4A", "SG-4B", "SG-5", "SG-6A", "SG-6B"]
    items.sort(key=lambda item: order.index(item.item_id) if item.item_id in order else len(order))
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
