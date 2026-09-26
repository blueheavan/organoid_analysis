"""SG-5: statistical qualification record for the frozen contrast method.

The record binds (1) engineering verification of ``statistics.small_sample``
against lmerTest and clubSandwich golden values and (2) the one-time
confirmation simulation of type-I error and 95 % coverage frozen in
docs/evidence/2026-09-26-statistical-qualification-protocol. Verification
re-derives the item status from the raw per-cell counts, re-runs the reference
comparison, and re-executes the worst contrast cell with its recorded seed,
which must reproduce the recorded counts exactly.
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

import numpy as np

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
)

CONTRACT_ID = "statistical-qualification/1"
ITEM = "SG-5"
PROTOCOL_DIR = "docs/evidence/2026-09-26-statistical-qualification-protocol"
PROTOCOL_PATH = f"{PROTOCOL_DIR}/PROTOCOL.md"
HARNESS_PATH = f"{PROTOCOL_DIR}/coverage_confirmation.py"
FREEZE_PATH = f"{PROTOCOL_DIR}/freeze.json"
REFERENCE_FILES = (f"{PROTOCOL_DIR}/reference/reference_datasets.csv",
                   f"{PROTOCOL_DIR}/reference/reference_values.json",
                   f"{PROTOCOL_DIR}/reference/reference.R",
                   f"{PROTOCOL_DIR}/reference/make_reference_datasets.py")
RAW_NAMES = ("sg5_contrast_coverage.csv", "sg5_condition_interval.csv", "sg5_run.log")
POINTER = "docs/evidence/statistical_qualification_current_record.json"
SEED = 20260928
N_SIM = 20000
N_CONTRAST_CELLS = 144
N_INTERVAL_CELLS = 9
TYPE1_RANGE = (0.04, 0.06)
COVERAGE_RANGE = (0.94, 0.96)
# Engineering-verification tolerances (relative), fixed with the protocol.
REFERENCE_RTOL = {"cr2": 1e-9, "lmm_estimate_se": 1e-6, "lmm_df": 1e-5, "lmm_p_abs": 1e-6}

SPEC = ContractSpec(
    contract_id=CONTRACT_ID,
    validation_scope=("SG-5 type-I error and 95% coverage of the frozen condition-contrast method and the "
                      "condition-level interval inside the qualified design envelope, with engineering "
                      "verification against lmerTest and clubSandwich"),
    items=(ITEM,),
    dependency_scope=(
        ScopeFile("src/organoid_analysis/statistics/small_sample.py", "method_core",
                  "REML fit, Satterthwaite df, CR2 fallback, omnibus, qualified design envelope."),
        ScopeFile("src/organoid_analysis/statistics/inference.py", "method_route",
                  "Object rows to replicate sufficient statistics, grouping key, log10 transform, labels."),
        ScopeFile("src/organoid_analysis/statistics/aggregation.py", "interval_core",
                  "condition_interval(): replicate-level t interval, log scale for size metrics."),
    ),
    protocol=(PROTOCOL_PATH, HARNESS_PATH, FREEZE_PATH, *REFERENCE_FILES),
    tooling=("src/organoid_analysis/validation/statistical_qualification_evidence.py",
             "src/organoid_analysis/validation/record_contract.py",
             "src/organoid_analysis/validation/evidence_manifest.py",
             "scripts/record_qualification_evidence.py"),
    raw_names=RAW_NAMES,
    pointer=POINTER,
)


def frozen_hashes(root: Path) -> dict[str, str]:
    path = root / FREEZE_PATH
    if not path.is_file():
        return {}
    document = loads_strict(path.read_text(encoding="utf-8"))
    return {str(k): str(v) for k, v in document.get("sha256", {}).items()}


# ------------------------------------------------------------ verification
def reference_check(root: Path) -> dict[str, Any]:
    """Compare the production method with the lmerTest / clubSandwich golden values."""
    import pandas as pd  # type: ignore[import-untyped]

    from organoid_analysis.statistics import small_sample as ss

    golden = json.loads((root / REFERENCE_FILES[1]).read_text(encoding="utf-8"))
    data = pd.read_csv(root / REFERENCE_FILES[0], float_precision="round_trip")
    worst = {"cr2": 0.0, "lmm_estimate_se": 0.0, "lmm_df": 0.0, "lmm_p_abs": 0.0}
    boundary_mismatch = []
    for name, reference in golden["datasets"].items():
        d = data[data.dataset == name]
        codes = d.condition.map({c: i for i, c in enumerate(sorted(d.condition.unique()))}).to_numpy()
        cd = ss.cluster_data(d.y.to_numpy(), d.replicate.to_numpy(), codes)
        fit = ss.fit(cd)
        if (fit.branch == ss.FALLBACK_BRANCH) != bool(reference["singular"]):
            boundary_mismatch.append(name)
        ols = ss.RandomInterceptFit(cd, fit.sigma_e2, 0.0, ss.FALLBACK_BRANCH, np.array(reference["ols"]))
        for pair in reference["pairs"]:
            c = np.zeros(cd.n_conditions)
            c[pair["j"]], c[pair["i"]] = 1.0, -1.0
            cr2 = ss.contrast(ols, c)
            worst["cr2"] = max(worst["cr2"], abs(cr2.standard_error ** 2 / pair["cr2"]["variance"] - 1),
                               abs(cr2.df / pair["cr2"]["df"] - 1))
            if fit.branch == ss.LMM_BRANCH:
                lmm = ss.contrast(fit, c)
                worst["lmm_estimate_se"] = max(worst["lmm_estimate_se"],
                                               abs(lmm.standard_error / pair["lmer"]["se"] - 1),
                                               abs(lmm.estimate - pair["lmer"]["estimate"]))
                worst["lmm_df"] = max(worst["lmm_df"], abs(lmm.df / pair["lmer"]["df"] - 1))
                worst["lmm_p_abs"] = max(worst["lmm_p_abs"], abs(lmm.p_value - pair["lmer"]["p"]))
    passed = not boundary_mismatch and all(worst[k] <= REFERENCE_RTOL[k] for k in worst)
    return {"passed": passed, "worst": worst, "tolerance": REFERENCE_RTOL,
            "boundary_mismatch": boundary_mismatch, "versions": golden["versions"],
            "n_datasets": len(golden["datasets"])}


def _read(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))
    except (OSError, csv.Error) as error:
        raise EvidenceError(f"raw results unreadable: {error}") from error


def _inside(value: float, bounds: tuple[float, float]) -> bool:
    return math.isfinite(value) and bounds[0] <= value <= bounds[1]


def derive_results(contrast_csv: Path, interval_csv: Path, reference: dict[str, Any]) -> dict[str, Any]:
    cells = _read(contrast_csv)
    intervals = _read(interval_csv)
    validity = []
    if len(cells) != N_CONTRAST_CELLS or len(intervals) != N_INTERVAL_CELLS:
        validity.append(f"expected {N_CONTRAST_CELLS}+{N_INTERVAL_CELLS} cells, found {len(cells)}+{len(intervals)}")
    if any(int(c["n_sim"]) != N_SIM or int(c["seed"]) != SEED for c in cells + intervals):
        validity.append("a cell was not run at the frozen n_sim and seed")
    if any(int(c["failures"]) for c in cells):
        validity.append("fit failures occurred")
    if not reference["passed"]:
        validity.append("reference verification failed")
    t1 = [(c, float(c["type1"])) for c in cells]
    cov = [(c, float(c["coverage"])) for c in cells]
    ci = [(c, float(c["coverage"])) for c in intervals]
    violations = ([_label(c) + f" type1 {v:.4f}" for c, v in t1 if not _inside(v, TYPE1_RANGE)]
                  + [_label(c) + f" coverage {v:.4f}" for c, v in cov if not _inside(v, COVERAGE_RANGE)]
                  + [f"interval G={c['G']} {c['dist']} coverage {v:.4f}" for c, v in ci
                     if not _inside(v, COVERAGE_RANGE)])
    status = "VOID" if validity else ("PASS" if not violations else "FAIL")
    omnibus = [float(c["omnibus_type1"]) for c in cells]
    return {ITEM: {
        "status": status, "criterion": "every cell: type-I in [0.04, 0.06], 95% coverage in [0.94, 0.96]; "
                                       "condition interval coverage in [0.94, 0.96]",
        "validity_problems": validity, "violations": violations[:40], "n_violations": len(violations),
        "contrast": {"n_cells": len(cells), "type1_range": [min(v for _, v in t1), max(v for _, v in t1)] if t1 else None,
                     "coverage_range": [min(v for _, v in cov), max(v for _, v in cov)] if cov else None,
                     "max_fallback_fraction": max(float(c["fallback_fraction"]) for c in cells) if cells else None},
        "condition_interval": {"n_cells": len(intervals),
                               "coverage_range": [min(v for _, v in ci), max(v for _, v in ci)] if ci else None},
        "omnibus_screen_type1_range": [min(omnibus), max(omnibus)] if omnibus else None,
        "reference_verification": {key: reference[key] for key in
                                   ("passed", "tolerance", "boundary_mismatch", "versions", "n_datasets")},
        "envelope": envelope(),
        "claim_boundary": ("method-level qualification inside the design envelope under approximately normal "
                           "replicate effects on the analysis scale; does not establish that a given experiment "
                           "is correctly designed or pre-registered")}}


def _label(cell: dict[str, str]) -> str:
    return f"K={cell['K']} G={cell['G']} m={cell['m']} {cell['balance']} icc={cell['icc']}"


def envelope() -> dict[str, Any]:
    from organoid_analysis.statistics import small_sample as ss
    return {"method_id": ss.METHOD_ID, "min_replicates": ss.QUALIFIED_MIN_REPLICATES,
            "max_replicates": ss.QUALIFIED_MAX_REPLICATES, "max_conditions": ss.QUALIFIED_MAX_CONDITIONS,
            "mean_objects": list(ss.QUALIFIED_MEAN_OBJECTS), "min_objects": ss.QUALIFIED_MIN_OBJECTS}


def load_harness(path: Path) -> ModuleType:
    name = f"_frozen_sg5_{sha256_file(path)[:16]}"
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


def reproduce_worst(root: Path, contrast_csv: Path) -> tuple[list[str], dict[str, Any]]:
    """Re-execute the contrast cell farthest from nominal; counts must match exactly."""
    cells = _read(contrast_csv)
    if not cells:
        return ["re-execution: no cells"], {}
    worst = max(cells, key=lambda c: max(abs(float(c["type1"]) - 0.05), abs(float(c["coverage"]) - 0.95)))
    harness = load_harness(root / HARNESS_PATH)
    cell = {"K": int(worst["K"]), "G": int(worst["G"]), "m": int(worst["m"]), "balance": worst["balance"],
            "icc": float(worst["icc"]), "dist": worst["dist"]}
    fresh = harness.run_cell((cell, int(worst["n_sim"]), int(worst["seed"])))
    keys = ("null_reject", "alt_cover", "omnibus_reject", "null_fallback", "failures")
    mismatches = [f"{k}: recorded {worst[k]}, re-executed {fresh[k]}" for k in keys if str(fresh[k]) != worst[k]]
    problems = [f"re-execution of {_label(worst)} does not reproduce: " + "; ".join(mismatches)] if mismatches else []
    return problems, {"cell": _label(worst), "keys": list(keys), "mismatches": len(mismatches)}


def run_harness(root: Path, run_dir: Path, processes: int) -> subprocess.CompletedProcess[bytes]:
    with (run_dir / RAW_NAMES[2]).open("w", encoding="utf-8") as log:
        done = subprocess.run([sys.executable, str(root / HARNESS_PATH), "--seed", str(SEED), "--n-sim", str(N_SIM),
                               "--out", str(run_dir / RAW_NAMES[0]), "--interval-out", str(run_dir / RAW_NAMES[1]),
                               "--processes", str(processes)], cwd=root, stdout=log, stderr=subprocess.STDOUT,
                              check=False)
        log.write(f"exit={done.returncode}\n")
    return done


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
    for path in (PROTOCOL_PATH, HARNESS_PATH, *REFERENCE_FILES):
        if path not in frozen:
            problems.append(f"frozen protocol: {path} has no hash in {FREEZE_PATH}")
    run_dir = root / Path(relpath).parent
    reproduction = None
    if (run_dir / RAW_NAMES[0]).is_file() and (run_dir / RAW_NAMES[1]).is_file():
        try:
            derived = derive_results(run_dir / RAW_NAMES[0], run_dir / RAW_NAMES[1], reference_check(root))
        except (EvidenceError, KeyError, ValueError) as error:
            problems.append(f"results: {error}")
        else:
            if loads_strict(dumps_strict(derived)) != manifest.get("results"):
                problems.append("results: recorded status/values are not those derived from the raw results")
            if reproduce:
                extra, reproduction = reproduce_worst(root, run_dir / RAW_NAMES[0])
                problems += extra
    return VerificationReport(relpath, tuple(problems), manifest, reproduction)
