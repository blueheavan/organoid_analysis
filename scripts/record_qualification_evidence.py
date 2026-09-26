"""Run a frozen qualification V&V and write a sealed, item-keyed validation record.

    pixi run qualification-record --contract analytical [--run-id RUN_ID] [--processes N]
    pixi run qualification-record --contract viability  [--run-id RUN_ID]

``analytical`` executes the frozen SG-1a/SG-2a confirmation grids
(docs/evidence/2026-09-26-analytical-qualification-protocol) through the
production geometry estimator; ``viability`` executes the SG-4A known-rule
cases. The source under validation is HEAD. The run is refused when any
evidence-critical file (dependency scope, frozen protocol, tooling, pixi.lock)
differs from HEAD, when the run directory exists, or when the protocol differs
from its freeze record. ``record_class`` is ``canonical`` only when the whole
tree was clean; only a canonical record can make a gate item PASS.
"""
from __future__ import annotations

import argparse
import datetime
import sys
from pathlib import Path

from organoid_analysis.validation import analytical_qualification_evidence as analytical
from organoid_analysis.validation import viability_rule_evidence as viability
from organoid_analysis.validation.evidence_manifest import dumps_strict, sha256_file
from organoid_analysis.validation.record_contract import (
    build_manifest,
    refuse_dirty_critical,
    utc_now,
)

ROOT = Path(__file__).resolve().parents[1]


def _fail(message: str) -> int:
    print(f"VALIDATION RECORD NOT WRITTEN: {message}", file=sys.stderr)
    return 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--contract", choices=("analytical", "viability"), required=True)
    parser.add_argument("--run-id")
    parser.add_argument("--processes", type=int, default=4)
    args = parser.parse_args(argv)
    module = analytical if args.contract == "analytical" else viability
    spec = module.SPEC
    run_id = args.run_id or f"{datetime.date.today().isoformat()}-{args.contract}-qualification-record"
    run_rel = f"docs/evidence/{run_id}"
    run_dir = ROOT / run_rel
    if run_dir.exists():
        return _fail(f"{run_rel} exists; evidence is never overwritten (choose another --run-id)")
    state, dirty = refuse_dirty_critical(ROOT, spec)
    if dirty:
        return _fail("evidence-critical files differ from HEAD; commit them first: " + ", ".join(dirty))
    critical = [e.path for e in spec.dependency_scope] + list(spec.protocol) + list(spec.tooling)
    before = {path: sha256_file(ROOT / path) for path in critical}
    invocation = " ".join(["pixi run qualification-record", *(argv if argv is not None else sys.argv[1:])])

    run_dir.mkdir(parents=True)
    started = utc_now()
    extra: dict[str, object]
    if module is analytical:
        frozen = analytical.frozen_hashes(ROOT)
        for path in (analytical.PROTOCOL_PATH, analytical.HARNESS_PATH, *analytical.GENERATOR_FILES):
            if frozen.get(path) != sha256_file(ROOT / path):
                return _fail(f"{path} is not the version frozen before execution")
        done = analytical.run_harness(ROOT, run_dir, args.processes)
        if done.returncode != 0:
            return _fail(f"the harness exited {done.returncode}; see {run_rel}/{analytical.RAW_NAMES[1]}")
        command = f"python {analytical.HARNESS_PATH} --out {run_rel}/{analytical.RAW_NAMES[0]} --processes {args.processes}"
        raw = run_dir / analytical.RAW_NAMES[0]
        results = analytical.derive_results(raw)
        problems, reproduction = analytical.reproduce_subset(ROOT / analytical.HARNESS_PATH, raw)
        if problems:
            return _fail("; ".join(problems))
        extra = {"domain": analytical.domain_declaration(), "reproduction": reproduction}
    else:
        rows = viability.run_cases()
        (run_dir / viability.RAW_NAME).write_text(viability.to_csv(rows), encoding="utf-8")
        command = "organoid_analysis.validation.viability_rule_evidence.run_cases"
        results = viability.derive_results(rows)
        extra = {"specification": {"source": "docs/INTENDED_USE_AND_ESTIMANDS.md sections 5.1-5.3",
                                   "conventions": viability.SPEC_CONVENTIONS}}
    finished = utc_now()
    if {path: sha256_file(ROOT / path) for path in critical} != before:
        return _fail("evidence-critical files changed while the V&V was running")

    manifest = build_manifest(ROOT, spec, run_id=run_id, state=state, run_rel=run_rel, results=results,
                              commands=[{"command": command, "cwd": ".", "exit_code": 0,
                                         "started_utc": started, "finished_utc": finished}],
                              invocation=invocation, extra=extra)
    manifest_rel = f"{run_rel}/evidence_manifest.json"
    (ROOT / manifest_rel).write_text(dumps_strict(manifest), encoding="utf-8")
    (ROOT / spec.pointer).write_text(dumps_strict({
        "contract_id": spec.contract_id, "manifest": manifest_rel,
        "note": "Selects the record the scientific gate verifies. Not hashed by any record; it can only choose "
                "among records that verify against the current repository."}), encoding="utf-8")
    report = module.verify_record(ROOT, manifest_rel)
    if not report.valid:
        return _fail("the new record does not verify: " + "; ".join(report.problems))
    print(f"validation record: {manifest_rel} ({manifest['record_class']}; source {state.head})")
    for item, result in results.items():
        print(f"{item}: {result['status']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
