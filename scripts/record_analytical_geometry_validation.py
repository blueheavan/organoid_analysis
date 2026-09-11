"""Run the frozen analytical-geometry V&V and write a sealed validation record.

    pixi run validation-record [--run-id RUN_ID]

Source commit versus validation record
--------------------------------------
The source under validation is the current HEAD commit. The evidence-critical
files (the declared dependency scope, the frozen plan, freeze record and
harness, and pixi.lock) must be unmodified relative to HEAD, otherwise the run
is refused: a record may only name a commit that contains exactly the code it
validated. The record itself is written afterwards into
``docs/evidence/<RUN_ID>/`` and is meant to be committed separately; it names
``validated_source_commit`` but never the commit that stores it.

``record_class`` is ``canonical`` when the whole working tree was clean at HEAD
and ``scope-clean`` when only files outside the evidence-critical set (for
example uncommitted validation tooling, whose hashes are recorded) were dirty.
The scientific gate reports FAIL/INSUFFICIENT EVIDENCE from either class but
grants PASS only from a canonical record.
"""
from __future__ import annotations

import argparse
import datetime
import json
import platform
import shutil
import subprocess
import sys
from pathlib import Path

from organoid_analysis.validation import analytical_geometry_evidence as contract
from organoid_analysis.validation.evidence_manifest import (
    SCHEMA_VERSION,
    dumps_strict,
    file_record,
    git_source_state,
    installed_versions,
    sanitize_non_finite,
    seal,
    sha256_file,
)

ROOT = Path(__file__).resolve().parents[1]


def _utc() -> str:
    return datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def _fail(message: str) -> int:
    print(f"VALIDATION RECORD NOT WRITTEN: {message}", file=sys.stderr)
    return 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run-id", default=f"{datetime.date.today().isoformat()}-analytical-geometry-record")
    args = parser.parse_args(argv)
    run_rel = f"docs/evidence/{args.run_id}"
    run_dir = ROOT / run_rel
    if run_dir.exists():
        return _fail(f"{run_rel} exists; evidence is never overwritten (choose another --run-id)")

    state = git_source_state(ROOT)
    scope_paths = [entry.path for entry in contract.DEPENDENCY_SCOPE]
    critical = [*scope_paths, contract.PLAN_PATH, contract.FREEZE_PATH, contract.FROZEN_HARNESS_PATH, contract.LOCKFILE]
    dirty_critical = sorted(set(critical) & set(state.dirty_paths))
    if dirty_critical:
        return _fail("evidence-critical files differ from HEAD; commit them first so the record names a commit "
                     f"that contains exactly the validated code: {', '.join(dirty_critical)}")
    frozen = contract.read_freeze(ROOT)
    for path in (contract.PLAN_PATH, contract.FROZEN_HARNESS_PATH):
        if frozen.get(path) != sha256_file(ROOT / path):
            return _fail(f"{path} is not the version frozen in {contract.FREEZE_PATH}")
    before = {path: sha256_file(ROOT / path) for path in critical}

    run_dir.mkdir(parents=True)
    shutil.copy2(ROOT / contract.FROZEN_HARNESS_PATH, run_dir / contract.HARNESS_NAME)
    harness_rel = f"{run_rel}/{contract.HARNESS_NAME}"
    command = ["python", harness_rel, "--grid", contract.GRID]
    started = _utc()
    with (run_dir / contract.RAW_NAMES[1]).open("w", encoding="utf-8") as log:
        done = subprocess.run([sys.executable, *command[1:]], cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=False)
        log.write(f"exit={done.returncode}\n")
    finished = _utc()
    if done.returncode != 0:
        return _fail(f"the V&V harness exited {done.returncode}; see {run_rel}/{contract.RAW_NAMES[1]}")

    # The frozen harness writes Python's non-strict JSON (NaN/Infinity for
    # non-estimable candidates); re-serialize it losslessly as strict JSON.
    summary_path = run_dir / contract.DERIVED_NAMES[0]
    summary, replaced = sanitize_non_finite(json.loads(summary_path.read_text(encoding="utf-8")))
    summary_path.write_text(dumps_strict(summary), encoding="utf-8")

    if {path: sha256_file(ROOT / path) for path in critical} != before:
        return _fail("evidence-critical files changed while the V&V was running")
    raw_csv = run_dir / contract.RAW_NAMES[0]
    results = contract.derive_results(raw_csv)
    reproduction_problems, reproduction = contract.reproduce_production(run_dir / contract.HARNESS_NAME, raw_csv)
    if reproduction_problems:
        return _fail("; ".join(reproduction_problems))

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "contract_id": contract.CONTRACT_ID,
        "validation_run_id": args.run_id,
        "validation_scope": contract.VALIDATION_SCOPE,
        "record_class": "canonical" if not state.dirty_paths else "scope-clean",
        "generated_utc": _utc(),
        "source": {
            "validated_source_commit": state.head,
            "validated_source_tree": state.tree,
            "whole_tree_clean_at_execution": not state.dirty_paths,
            "dirty_paths_at_execution": list(state.dirty_paths),
            "scope_clean_at_execution": True,
            "record_storage": contract.RECORD_STORAGE_NOTE,
        },
        "dependency_scope": {
            "description": "Files whose content can change SG-1/SG-2 values or their exported meaning. Any change "
                           "rejects this record.",
            "files": [file_record(ROOT, e.path, role=e.role, rationale=e.rationale) for e in contract.DEPENDENCY_SCOPE],
            "excluded": [{"paths": paths, "rationale": why} for paths, why in contract.EXCLUDED_FROM_SCOPE],
        },
        "environment": {
            "lockfile": file_record(ROOT, contract.LOCKFILE),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "packages": installed_versions(contract.CRITICAL_PACKAGES),
        },
        "vv_plan": [
            file_record(ROOT, contract.PLAN_PATH, role="frozen V&V plan"),
            file_record(ROOT, contract.FREEZE_PATH, role="freeze record: plan and harness SHA-256 before any result"),
        ],
        "vv_implementation": [
            file_record(ROOT, harness_rel, role="V&V harness executed", byte_identical_to=contract.FROZEN_HARNESS_PATH),
            *(file_record(ROOT, path, role="derivation, verification and recording tooling")
              for path in contract.VV_TOOLING),
        ],
        "raw_artifacts": [file_record(ROOT, f"{run_rel}/{name}") for name in contract.RAW_NAMES],
        "derived_artifacts": [
            file_record(ROOT, f"{run_rel}/{contract.DERIVED_NAMES[0]}",
                        strict_json_replacements=[{"pointer": p, "original": o} for p, o in replaced]),
            file_record(ROOT, f"{run_rel}/{contract.DERIVED_NAMES[1]}"),
        ],
        "estimator": contract.production_estimator_identity(),
        "acceptance_criteria": contract.ACCEPTANCE_CRITERIA,
        "commands": [
            {"command": " ".join(command), "cwd": ".", "exit_code": done.returncode,
             "started_utc": started, "finished_utc": finished},
            {"command": "organoid_analysis.validation.analytical_geometry_evidence.reproduce_production", "exit_code": 0},
        ],
        "invocation": " ".join(["pixi run validation-record", *(argv if argv is not None else sys.argv[1:])]).strip(),
        "results": results,
        "resulting_status": {item: result["status"] for item, result in results.items()},
        "reproduction": reproduction,
    }
    manifest_rel = f"{run_rel}/{contract.MANIFEST_NAME}"
    (ROOT / manifest_rel).write_text(dumps_strict(seal(manifest)), encoding="utf-8")
    (ROOT / contract.CURRENT_RECORD_POINTER).write_text(dumps_strict({
        "contract_id": contract.CONTRACT_ID, "manifest": manifest_rel,
        "note": "Selects the record the scientific gate verifies. Not hashed by any record; it can only choose "
                "among records that verify against the current repository."}), encoding="utf-8")

    report = contract.verify_record(ROOT, manifest_rel)
    if not report.valid:
        return _fail("the new record does not verify: " + "; ".join(report.problems))
    print(f"validation record: {manifest_rel}")
    print(f"record class: {manifest['record_class']}; validated source commit {state.head}")
    for item, result in results.items():
        print(f"{item}: {result['status']} (max |error| {result['max_abs_rel_error']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
