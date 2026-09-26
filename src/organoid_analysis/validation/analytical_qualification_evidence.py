"""SG-1a / SG-2a: domain-restricted analytical qualification record (owner decision D-14).

The characterization record (``analytical_geometry_evidence``, SG-1b/SG-2b)
measures the estimators everywhere on the 2026-09-11 development grid and
carries no verdict. This contract answers the qualification questions: is the
production surface area within 5 % and the voxel-count volume within 1 % of the
closed-form truth on every case of a *confirmation* set that lies inside the
domain declared for that quantity before the set was executed?

* SG-1a adjudicates the surface domain frozen 2026-09-12 (``rho_in >= 10``,
  anisotropy <= 4, smooth closed) on two confirmation grids: ``confirm2``
  (frozen 2026-09-12) and ``confirm3`` (frozen in the protocol of this contract
  before its first execution; adds lattice-symmetric placements).
* SG-2a adjudicates the volume-specific domain (``surface_crofton.in_volume_domain``)
  derived by the 2026-09-26 development study, on ``confirm3`` only: confirm2's
  volume errors were known when that domain was derived, so using them would be
  post hoc. They are reported as supplementary characterization.

Both items read confirmation evidence only; neither reads a development result
(GATE_SEMANTICS_ANALYSIS.md section 7.3 invariant 3). The frozen generator
files are hash-checked against the 2026-09-12 freeze record and the harness and
protocol against this contract's own freeze file.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import math
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

from organoid_analysis.quantification import surface_crofton
from organoid_analysis.validation.analytical_geometry_evidence import (
    DEPENDENCY_SCOPE as GEOMETRY_SCOPE,
)
from organoid_analysis.validation.evidence_manifest import (
    EvidenceError,
    dumps_strict,
    loads_strict,
    sha256_file,
)
from organoid_analysis.validation.record_contract import (
    ContractSpec,
    ScopeFile,
    VerificationReport,
    current_manifest_relpath,
    verify_common,
    worst,
)

CONTRACT_ID = "analytical-qualification/1"
PROTOCOL_DIR = "docs/evidence/2026-09-26-analytical-qualification-protocol"
PROTOCOL_PATH = f"{PROTOCOL_DIR}/PROTOCOL.md"
HARNESS_PATH = f"{PROTOCOL_DIR}/qualification_harness.py"
FREEZE_PATH = f"{PROTOCOL_DIR}/freeze.json"
GENERATOR_DIR = "docs/evidence/2026-09-12-surface-method-development"
GENERATOR_FREEZE = f"{GENERATOR_DIR}/freeze_record_v3.json"
GENERATOR_FILES = (f"{GENERATOR_DIR}/code/common.py", f"{GENERATOR_DIR}/code/common2.py")
RAW_NAMES = ("qualification_cases.csv", "qualification_run.log")
POINTER = "docs/evidence/analytical_qualification_current_record.json"

SPEC = ContractSpec(
    contract_id=CONTRACT_ID,
    validation_scope=("SG-1a surface-area and SG-2a voxel-count-volume qualification of the production geometry "
                      "estimator inside their declared domains, on prospectively frozen confirmation grids"),
    items=("SG-1a", "SG-2a"),
    dependency_scope=tuple(ScopeFile(e.path, e.role, e.rationale) for e in GEOMETRY_SCOPE),
    protocol=(PROTOCOL_PATH, HARNESS_PATH, FREEZE_PATH, GENERATOR_FREEZE, *GENERATOR_FILES),
    tooling=("src/organoid_analysis/validation/analytical_qualification_evidence.py",
             "src/organoid_analysis/validation/record_contract.py",
             "src/organoid_analysis/validation/evidence_manifest.py",
             "src/organoid_analysis/validation/analytical_geometry_evidence.py",
             "scripts/record_qualification_evidence.py"),
    raw_names=RAW_NAMES,
    pointer=POINTER,
)

SURFACE_THRESHOLD = 0.05
VOLUME_THRESHOLD = 0.01
SMOOTH = ("sphere", "ellipsoid", "capsule", "torus")
REPRODUCTION_RTOL = 1e-9
REPRODUCE_WORST_PER_ITEM = 3
REPRODUCE_EVERY = 120


def frozen_hashes(root: Path) -> dict[str, str]:
    """Path -> SHA-256 fixed before execution: this protocol's freeze file and the generator's 2026-09-12 freeze."""
    hashes: dict[str, str] = {}
    freeze = root / FREEZE_PATH
    if freeze.is_file():
        document = loads_strict(freeze.read_text(encoding="utf-8"))
        hashes.update({str(k): str(v) for k, v in document.get("sha256", {}).items()})
    generator = root / GENERATOR_FREEZE
    if generator.is_file():
        document = json.loads(generator.read_text(encoding="utf-8"))
        for path in GENERATOR_FILES:
            name = Path(path).name
            if name in document:
                hashes[path] = str(document[name])
    return hashes


def domain_declaration() -> dict[str, Any]:
    return {"surface": {"rho_in_min": surface_crofton.DOMAIN_RHO_IN_MIN,
                        "anisotropy_max": surface_crofton.DOMAIN_ANISO_MAX,
                        "scope": surface_crofton.DOMAIN_SCOPE, "weights": "evidence-bearing packaged tables only"},
            "volume": {"rho_in_min": surface_crofton.DOMAIN_RHO_IN_MIN,
                       "anisotropy_max": surface_crofton.DOMAIN_ANISO_MAX,
                       "rho_vol_min": surface_crofton.DOMAIN_RHO_VOL_MIN,
                       "rho_vol_definition": "inscribed radius / cbrt(voxel volume)",
                       "scope": surface_crofton.DOMAIN_SCOPE}}


# ------------------------------------------------------------------ derivation
def read_rows(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
    except (OSError, csv.Error) as error:
        raise EvidenceError(f"raw results unreadable: {error}") from error
    ids = [row.get("case_id") for row in rows]
    if not rows or len(set(ids)) != len(ids):
        raise EvidenceError("raw results are empty or contain duplicate cases")
    return rows


def _flag(row: dict[str, str], column: str) -> bool:
    return row.get(column) == "True"


def _err(row: dict[str, str], column: str) -> float:
    try:
        value = float(row[column])
    except (KeyError, TypeError, ValueError) as error:
        raise EvidenceError(f"case {row.get('case_id')!r}: no numeric {column!r}") from error
    return value


def _summary(rows: list[dict[str, str]], column: str) -> dict[str, Any]:
    case, value = worst((row["case_id"], _err(row, column)) for row in rows)
    return {"n": len(rows), "max_abs_rel_error": value, "worst_case": case,
            "classes": sorted({row["shape"] for row in rows}),
            "pose_types": sorted({"symmetric" if row["pose"].startswith("sym") else "random" for row in rows})}


def _status(sets: list[list[dict[str, str]]], column: str, threshold: float, validity: list[str]) -> str:
    if validity:
        return "VOID"
    values = [abs(_err(row, column)) for rows in sets for row in rows]
    return "PASS" if values and all(math.isfinite(v) and v < threshold for v in values) else "FAIL"


def derive_results(csv_path: Path) -> dict[str, Any]:
    rows = read_rows(csv_path)
    grids = {name: [row for row in rows if row["grid"] == name] for name in ("confirm2", "confirm3")}
    errors = [row["case_id"] for row in rows if row.get("error")]

    # ---- SG-1a
    surface_sets = {}
    validity_1a = []
    excluded_weights = 0
    for name, grid_rows in grids.items():
        in_scope = [r for r in grid_rows if _flag(r, "smooth") and _flag(r, "surface_in_qualified_domain")]
        excluded_weights += sum(not _flag(r, "weights_evidence_bearing") for r in in_scope)
        surface_sets[name] = [r for r in in_scope if _flag(r, "weights_evidence_bearing")]
        if {r["shape"] for r in surface_sets[name]} != set(SMOOTH):
            validity_1a.append(f"{name}: the adjudicated set lacks a smooth class")
    if {r["pose"].startswith("sym") for r in surface_sets["confirm3"]} != {True, False}:
        validity_1a.append("confirm3: the adjudicated set lacks symmetric or random placements")
    if errors:
        validity_1a.append(f"{len(errors)} cases failed to measure")
    sg1a = {"status": _status(list(surface_sets.values()), "area_rel_err", SURFACE_THRESHOLD, validity_1a),
            "criterion": f"every adjudicated case |area_rel_err| < {SURFACE_THRESHOLD}",
            "domain": domain_declaration()["surface"], "validity_problems": validity_1a,
            "by_grid": {name: _summary(s, "area_rel_err") for name, s in surface_sets.items()},
            "excluded_non_evidence_bearing_weights": excluded_weights,
            "not_adjudicated": {"creased_confirm2": _summary(
                [r for r in grids["confirm2"] if not _flag(r, "smooth") and not r.get("error")], "area_rel_err")}}

    # ---- SG-2a (confirm3 only)
    c3 = grids["confirm3"]
    volume_set = [r for r in c3 if _flag(r, "smooth") and _flag(r, "volume_in_qualified_domain")]
    validity_2a = []
    if {r["shape"] for r in volume_set} != set(SMOOTH):
        validity_2a.append("confirm3: the adjudicated volume set lacks a smooth class")
    if {r["pose"].startswith("sym") for r in volume_set} != {True, False}:
        validity_2a.append("confirm3: the adjudicated volume set lacks symmetric or random placements")
    if not any(abs(_err(r, "anisotropy") - 1.0) < 1e-12 for r in volume_set):
        validity_2a.append("confirm3: no isotropic case in the adjudicated volume set")
    if errors:
        validity_2a.append(f"{len(errors)} cases failed to measure")
    d2_only = [r for r in c3 if _flag(r, "smooth") and _flag(r, "surface_in_qualified_domain")
               and not _flag(r, "volume_in_qualified_domain")]
    c2_post_hoc = [r for r in grids["confirm2"] if _flag(r, "smooth") and _flag(r, "volume_in_qualified_domain")]
    sg2a = {"status": _status([volume_set], "volume_rel_err", VOLUME_THRESHOLD, validity_2a),
            "criterion": f"every adjudicated case |volume_rel_err| < {VOLUME_THRESHOLD}",
            "domain": domain_declaration()["volume"], "validity_problems": validity_2a,
            "confirm3": _summary(volume_set, "volume_rel_err"),
            "not_adjudicated": {
                "confirm3_in_surface_domain_below_volume_domain": _summary(d2_only, "volume_rel_err"),
                "confirm2_in_volume_domain_post_hoc": _summary(c2_post_hoc, "volume_rel_err")}}
    return {"SG-1a": sg1a, "SG-2a": sg2a}


# --------------------------------------------------------------- re-execution
def load_harness(path: Path) -> ModuleType:
    name = f"_frozen_qualification_{sha256_file(path)[:16]}"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise EvidenceError(f"cannot load harness {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        sys.modules.pop(name, None)
    return module


def reproduction_subset(rows: list[dict[str, str]]) -> list[str]:
    """Deterministic subset: the worst adjudicated cases of each item plus every Nth row."""
    derived = {"area_rel_err": [r for r in rows if _flag(r, "smooth") and _flag(r, "surface_in_qualified_domain")],
               "volume_rel_err": [r for r in rows if r["grid"] == "confirm3" and _flag(r, "smooth")
                                  and _flag(r, "volume_in_qualified_domain")]}
    chosen: list[str] = []
    for column, subset in derived.items():
        ranked = sorted(subset, key=lambda r: -abs(_err(r, column)))
        chosen += [r["case_id"] for r in ranked[:REPRODUCE_WORST_PER_ITEM]]
    chosen += [r["case_id"] for r in rows[::REPRODUCE_EVERY]]
    return list(dict.fromkeys(chosen))


def reproduce_subset(harness_path: Path, csv_path: Path) -> tuple[list[str], dict[str, Any]]:
    rows = {row["case_id"]: row for row in read_rows(csv_path)}
    harness = load_harness(harness_path)
    cases = {case["case_id"]: case for case in harness.all_cases()}
    problems: list[str] = []
    missing = sorted(set(cases) - set(rows))
    extra = sorted(set(rows) - set(cases))
    if missing or extra:
        problems.append(f"re-execution: raw results do not cover the frozen grids exactly "
                        f"({len(missing)} missing, {len(extra)} extra)")
    subset = [case_id for case_id in reproduction_subset(list(rows.values())) if case_id in cases]
    largest = 0.0
    mismatches = []
    for case_id in subset:
        fresh = harness.measure_case(cases[case_id])
        recorded = rows[case_id]
        for column in ("area_est", "volume_est", "rho_in", "rho_vol", "area_true", "volume_true"):
            a, b = float(recorded[column]), float(fresh[column])
            difference = abs(a - b) / max(abs(a), abs(b)) if a != b else 0.0
            largest = max(largest, difference)
            if difference > REPRODUCTION_RTOL:
                mismatches.append(f"{case_id}: {column} recorded {a!r}, re-executed {b!r}")
        for column in ("n_voxels", "surface_in_qualified_domain", "volume_in_qualified_domain",
                       "weights_evidence_bearing"):
            if str(fresh[column]) != recorded[column]:
                mismatches.append(f"{case_id}: {column} recorded {recorded[column]}, re-executed {fresh[column]}")
    if mismatches:
        problems.append(f"re-execution: {len(mismatches)} values not reproduced: " + "; ".join(mismatches[:5]))
    return problems, {"method": "re-executed a deterministic subset (worst adjudicated cases of each item and "
                                f"every {REPRODUCE_EVERY}th case) through the production estimator",
                      "cases_reexecuted": len(subset), "rtol": REPRODUCTION_RTOL,
                      "max_relative_difference": largest, "mismatches": len(mismatches)}


def run_harness(root: Path, run_dir: Path, processes: int) -> subprocess.CompletedProcess[bytes]:
    with (run_dir / RAW_NAMES[1]).open("w", encoding="utf-8") as log:
        done = subprocess.run([sys.executable, str(root / HARNESS_PATH), "--out", str(run_dir / RAW_NAMES[0]),
                               "--processes", str(processes)], cwd=root, stdout=log, stderr=subprocess.STDOUT,
                              check=False)
        log.write(f"exit={done.returncode}\n")
    return done


# --------------------------------------------------------------- verification
def verify_record(root: Path, manifest_relpath: str | None = None, *, check_git: bool = True,
                  reproduce: bool = True) -> VerificationReport:
    try:
        relpath = manifest_relpath or current_manifest_relpath(root, SPEC)
    except EvidenceError as error:
        return VerificationReport(manifest_relpath, (str(error),))
    manifest, problems = verify_common(root, SPEC, relpath, check_git=check_git, frozen_hashes=frozen_hashes)
    if manifest is None:
        return VerificationReport(relpath, tuple(problems))
    frozen = frozen_hashes(root)
    for path in (PROTOCOL_PATH, HARNESS_PATH, *GENERATOR_FILES):
        if path not in frozen:
            problems.append(f"frozen protocol: {path} has no hash in a freeze record")
    if manifest.get("domain") != domain_declaration():
        problems.append("domain: the recorded domain constants differ from the production constants")
    raw = root / Path(relpath).parent / RAW_NAMES[0]
    reproduction = None
    if raw.is_file():
        try:
            derived = derive_results(raw)
        except EvidenceError as error:
            problems.append(str(error))
        else:
            if loads_strict(dumps_strict(derived)) != manifest.get("results"):
                problems.append("results: recorded statuses/values are not those derived from the raw results")
            if reproduce:
                try:
                    extra, reproduction = reproduce_subset(root / HARNESS_PATH, raw)
                except EvidenceError as error:
                    extra = [f"re-execution: {error}"]
                problems += extra
    return VerificationReport(relpath, tuple(problems), manifest, reproduction)
