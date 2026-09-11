"""Evidence contract for the analytical-geometry scientific gate items.

SG-1 (surface area) and SG-2 (voxel-count volume) are measured on the phantoms
of the frozen surface V&V plan (docs/evidence/2026-09-11-measurement-vv). This
module declares which files can change those results (``DEPENDENCY_SCOPE``),
derives both item statuses from the raw per-case results, re-executes the
production estimator on every phantom, and verifies a recorded validation run
against the current repository (``verify_record``).

The criteria are those of docs/SCIENTIFIC_SPEC.md section 9 and are fixed here;
a record carrying any other criterion is rejected, never honored.
"""
from __future__ import annotations

import csv
import importlib.util
import math
import platform
import statistics
import sys
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from types import ModuleType
from typing import Any

from organoid_analysis.quantification.features import SURFACE_AREA_METHOD
from organoid_analysis.validation.evidence_manifest import (
    SCHEMA_VERSION,
    EvidenceError,
    dumps_strict,
    installed_versions,
    load_manifest,
    loads_strict,
    module_origin,
    recorded_paths,
    sha256_file,
    verify_file_records,
    verify_git_provenance,
    verify_seal,
)

CONTRACT_ID = "analytical-geometry/1"
VALIDATION_SCOPE = ("SG-1 analytical surface area and SG-2 analytical voxel-count volume of the production "
                    "geometry estimator on the frozen development phantom grid")
RECORD_CLASSES = ("canonical", "scope-clean")
RECORD_STORAGE_NOTE = (
    "This record is stored in a commit made after validated_source_commit. The storing commit is "
    "deliberately not named inside the record (a commit cannot contain its own hash); find it with "
    "`git log --format=%H -1 -- <manifest path>`. validated_source_commit identifies the source under "
    "validation, never the commit that stores this record.")

HISTORICAL_VV_DIR = "docs/evidence/2026-09-11-measurement-vv"
PLAN_PATH = f"{HISTORICAL_VV_DIR}/SURFACE_VV_PLAN.md"
FREEZE_PATH = f"{HISTORICAL_VV_DIR}/surface_vv_freeze.txt"
FROZEN_HARNESS_PATH = f"{HISTORICAL_VV_DIR}/surface_vv.py"
CURRENT_RECORD_POINTER = "docs/evidence/analytical_geometry_current_record.json"
MANIFEST_NAME = "evidence_manifest.json"
HARNESS_NAME = "surface_vv.py"
GRID = "dev"
RAW_NAMES = ("surface_vv_dev.csv", "surface_vv_dev.log")
DERIVED_NAMES = ("surface_vv_dev_summary.json", "surface_vv_dev_summary.md")
LOCKFILE = "pixi.lock"
# Harness label of the production estimator (quantification.features.geometry).
PRODUCTION_ESTIMATOR_ID = "E0_marching_cubes_binary_v1"
# The raw CSV stores 10 significant digits (float_format "%.10g"); re-executed
# values must agree to that precision.
REPRODUCTION_RTOL = 1e-9
CRITICAL_PACKAGES = ("numpy", "scipy", "scikit-image", "pandas")
# Code that derives or records the results; a change requires a new record.
VV_TOOLING = (
    "src/organoid_analysis/validation/analytical_geometry_evidence.py",
    "src/organoid_analysis/validation/evidence_manifest.py",
    "scripts/record_analytical_geometry_validation.py",
)

SURFACE_AREA_THRESHOLD = 0.05
VOLUME_THRESHOLD = 0.01
ACCEPTANCE_CRITERIA: dict[str, dict[str, Any]] = {
    "SG-1": {"requirement": "Analytical surface area of the production estimator",
             "metric": "max over all phantom cases of |A_est / A_true - 1|; a non-estimable case fails",
             "comparison": "<", "threshold": SURFACE_AREA_THRESHOLD, "source": "docs/SCIENTIFIC_SPEC.md section 9"},
    "SG-2": {"requirement": "Analytical voxel-count volume",
             "metric": "max over all phantom cases of |V_voxel / V_true - 1|; a non-estimable case fails",
             "comparison": "<", "threshold": VOLUME_THRESHOLD, "source": "docs/SCIENTIFIC_SPEC.md section 9"},
}


@dataclass(frozen=True)
class ScopeEntry:
    path: str
    role: str
    rationale: str


# Every production file whose content can change a validated SG-1/SG-2 value or
# what that value means in an exported measurement.
DEPENDENCY_SCOPE = (
    ScopeEntry("src/organoid_analysis/quantification/features.py", "estimator_core",
               "geometry(): voxel-count volume, marching-cubes surface area and sphericity; surface_mesh(): "
               "isosurface level, padding and physical spacing; outer_envelope(): filled support; "
               "SURFACE_AREA_METHOD identity. Executed for every phantom."),
    ScopeEntry("src/organoid_analysis/quantification/mask_features.py", "support_and_spacing_semantics",
               "Web/mask-feature route: raw-label default versus optional filled support handed to geometry(), "
               "and the (x, y, z) to (z, y, x) spacing conversion."),
    ScopeEntry("src/organoid_analysis/quantification/cellular_measurements.py", "support_semantics",
               "cell_geometry() calls geometry() with the filled-envelope default for cell objects."),
    ScopeEntry("src/organoid_analysis/quantification/multilevel_relationships/morphology.py", "support_semantics",
               "Multilevel morphology calls geometry() on raw labels (fill_holes=False)."),
    ScopeEntry("src/organoid_analysis/microscopy_io/voxel_spacing.py", "spacing_interpretation",
               "validate_voxel_spacing_xyz(): axis order, finiteness and positivity of the spacing that the "
               "mask-feature route passes to geometry()."),
)
EXCLUDED_FROM_SCOPE = (
    ("src/organoid_analysis/microscopy_io/tiff_contract.py, metadata.py, zstack_reader.py",
     "Acquisition spacing is parsed upstream; whether it is correct is gate item SG-6, not an analytical-phantom "
     "result. The V&V passes spacing explicitly."),
    ("src/organoid_analysis/workflows/*",
     "Orchestration; reaches geometry() only through the scoped modules (classical route: "
     "features.measure_instances)."),
    ("src/organoid_analysis/segmentation/watershed_instances.py, src/organoid_analysis/quantification/labels.py",
     "Imported for types, crops and border QC flags; they do not enter area or volume."),
    ("numpy, scipy, scikit-image",
     "Bound by the pixi.lock SHA-256 and recorded package versions, and exercised by re-execution."),
    ("package __init__ modules", "No numerical code (organoid_analysis/__init__.py holds the version string)."),
    ("scripts/scientific_validation_gate.py",
     "Presents verified results and applies the PASS rule; it does not derive any value."),
)


def production_estimator_identity() -> dict[str, Any]:
    return {"harness_label": PRODUCTION_ESTIMATOR_ID,
            "surface_area_method": dict(SURFACE_AREA_METHOD),
            "volume_method": "voxel count x product of ZYX spacing (quantification.features.geometry)"}


def read_freeze(root: Path) -> dict[str, str]:
    """Map path -> SHA-256 from the freeze record written before any result existed."""
    frozen: dict[str, str] = {}
    path = root / FREEZE_PATH
    if not path.is_file():
        return frozen
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) == 2 and len(parts[0]) == 64:
            frozen[parts[1]] = parts[0]
    return frozen


def current_manifest_relpath(root: Path) -> str:
    pointer = root / CURRENT_RECORD_POINTER
    if not pointer.is_file():
        raise EvidenceError(f"no analytical-geometry validation record ({CURRENT_RECORD_POINTER} missing)")
    document = loads_strict(pointer.read_text(encoding="utf-8"))
    if not isinstance(document, dict) or document.get("contract_id") != CONTRACT_ID:
        raise EvidenceError(f"{CURRENT_RECORD_POINTER} does not name a {CONTRACT_ID} record")
    manifest = document.get("manifest")
    if not isinstance(manifest, str) or PurePosixPath(manifest).name != MANIFEST_NAME:
        raise EvidenceError(f"{CURRENT_RECORD_POINTER}: invalid manifest path {manifest!r}")
    return manifest


# ------------------------------------------------------- results derivation
def _production_rows(csv_path: Path) -> list[dict[str, str]]:
    try:
        with csv_path.open(newline="", encoding="utf-8") as handle:
            rows = [row for row in csv.DictReader(handle) if row.get("estimator") == PRODUCTION_ESTIMATOR_ID]
    except (OSError, csv.Error) as error:
        raise EvidenceError(f"raw results unreadable: {error}") from error
    if not rows:
        raise EvidenceError(f"raw results contain no {PRODUCTION_ESTIMATOR_ID} rows")
    case_ids = [row.get("case_id") for row in rows]
    if len(set(case_ids)) != len(case_ids):
        raise EvidenceError("raw results contain duplicate production cases")
    return rows


def _number(row: dict[str, str], column: str) -> float:
    try:
        return float(row[column])
    except (KeyError, TypeError, ValueError) as error:
        raise EvidenceError(f"raw results: case {row.get('case_id')!r} has no numeric {column!r}") from error


def _worst(errors: list[tuple[str, float]], threshold: float) -> dict[str, Any]:
    finite = [(case, value) for case, value in errors if math.isfinite(value)]
    nonestimable = sorted(case for case, value in errors if not math.isfinite(value))
    worst = max(finite, key=lambda item: abs(item[1])) if finite else None
    passed = not nonestimable and worst is not None and abs(worst[1]) < threshold
    return {"status": "PASS" if passed else "FAIL", "n_cases": len(errors), "n_nonestimable": len(nonestimable),
            "nonestimable_cases": nonestimable,
            "max_abs_rel_error": abs(worst[1]) if worst else None, "worst_case": worst[0] if worst else None,
            "median_signed_rel_error": statistics.median(v for _, v in finite) if finite else None}


def derive_results(csv_path: Path) -> dict[str, Any]:
    """Derive SG-1 and SG-2 from the raw per-case CSV alone (every case counts)."""
    rows = _production_rows(csv_path)
    area = [(row["case_id"], _number(row, "area_rel_err")) for row in rows]
    volume = [(row["case_id"], _number(row, "volume_rel_err")) for row in rows]
    strata: dict[str, list[tuple[str, float]]] = {}
    for row, value in zip(rows, volume):
        strata.setdefault(row.get("rho_stratum") or "", []).append(value)
    by_stratum = [{"rho_stratum": name, **_worst(values, VOLUME_THRESHOLD)} for name, values in sorted(strata.items())]
    sg2 = _worst(volume, VOLUME_THRESHOLD)
    sg2["by_rho_stratum"] = by_stratum
    return {"SG-1": _worst(area, SURFACE_AREA_THRESHOLD), "SG-2": sg2}


# ------------------------------------------------------------- re-execution
def load_harness(path: Path) -> ModuleType:
    name = f"_frozen_surface_vv_{sha256_file(path)[:16]}"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise EvidenceError(f"cannot load V&V harness {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        sys.modules.pop(name, None)
    return module


def _relative_difference(a: float, b: float) -> float:
    if a == b or (math.isnan(a) and math.isnan(b)):
        return 0.0
    scale = max(abs(a), abs(b))
    return abs(a - b) / scale if math.isfinite(scale) and scale > 0 else math.inf


def reproduce_production(harness_path: Path, csv_path: Path) -> tuple[list[str], dict[str, Any]]:
    """Re-run the production estimator and the closed-form oracles on every phantom.

    Recomputes voxel count, voxel volume, true and estimated area, both relative
    errors and the rho stratum, and compares them with the raw CSV. A forged or
    stale raw result, or production code that no longer produces the recorded
    values, is reported even when every hash has been updated consistently.
    """
    rows = {row["case_id"]: row for row in _production_rows(csv_path)}
    harness = load_harness(harness_path)
    grid = harness.GRIDS[GRID]
    mismatches: list[str] = []
    seen: set[str] = set()
    largest = 0.0
    for index, (kind, params, spacing, variant, rotation, offset, min_axis) in enumerate(harness.cases(grid)):
        size = params["axes"][0] if kind == "ellipsoid" else params["r"]
        case_id = f"{GRID}-{index:03d}-{kind}-size{size:g}-sp{'x'.join(f'{s:g}' for s in spacing)}-{variant}"
        seen.add(case_id)
        row = rows.get(case_id)
        if row is None:
            mismatches.append(f"{case_id}: missing from the raw results")
            continue
        mask = harness.rasterize(kind, params, spacing, rotation, offset)
        voxel_volume = float(math.prod(spacing))
        area_true, volume_true = (float(value) for value in harness.truth(kind, params))
        area_est = float(harness.e0_production(mask, spacing))
        n_voxels = int(mask.sum())
        recomputed = {"area_true_um2": area_true, "area_est_um2": area_est, "area_rel_err": area_est / area_true - 1,
                      "volume_true_um3": volume_true, "volume_voxel_um3": n_voxels * voxel_volume,
                      "volume_rel_err": n_voxels * voxel_volume / volume_true - 1}
        if str(n_voxels) != row.get("n_voxels"):
            mismatches.append(f"{case_id}: n_voxels recorded {row.get('n_voxels')}, re-executed {n_voxels}")
        stratum = str(harness.rho_stratum(min_axis / max(spacing)))
        if stratum != row.get("rho_stratum"):
            mismatches.append(f"{case_id}: rho_stratum recorded {row.get('rho_stratum')}, re-executed {stratum}")
        for column, value in recomputed.items():
            difference = _relative_difference(_number(row, column), value)
            largest = max(largest, difference)
            if difference > REPRODUCTION_RTOL:
                mismatches.append(f"{case_id}: {column} recorded {row.get(column)}, re-executed {value!r}")
    extra = sorted(set(rows) - seen)
    mismatches += [f"{case_id}: not a case of the frozen {GRID} grid" for case_id in extra]
    problems = []
    if mismatches:
        shown = "; ".join(mismatches[:5]) + (f"; ... {len(mismatches) - 5} more" if len(mismatches) > 5 else "")
        problems.append(f"re-execution: {len(mismatches)} raw values are not reproduced by the production "
                        f"estimator and closed-form oracles: {shown}")
    summary = {"method": "re-executed the production estimator and closed-form oracles on every frozen-grid "
                         "phantom and compared them with the raw CSV",
               "cases_reexecuted": len(seen), "rtol": REPRODUCTION_RTOL,
               "max_relative_difference": largest, "mismatches": len(mismatches)}
    return problems, summary


# --------------------------------------------------------------- verification
@dataclass(frozen=True)
class VerificationReport:
    manifest_path: str | None
    problems: tuple[str, ...]
    manifest: dict[str, Any] | None = None
    reproduction: dict[str, Any] | None = field(default=None)

    @property
    def valid(self) -> bool:
        return self.manifest is not None and not self.problems


def referenced_paths(manifest: dict[str, Any], manifest_relpath: str) -> set[str]:
    """Every repository file the record depends on (used to copy a record for testing)."""
    paths = {manifest_relpath, PLAN_PATH, FREEZE_PATH, FROZEN_HARNESS_PATH, CURRENT_RECORD_POINTER}
    for group in ("vv_plan", "vv_implementation", "raw_artifacts", "derived_artifacts"):
        paths |= recorded_paths(manifest.get(group))
    paths |= recorded_paths(manifest.get("dependency_scope", {}).get("files"))
    paths |= recorded_paths([manifest.get("environment", {}).get("lockfile")])
    return paths


def _as_dict(value: object) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _check_required(records: object, required: set[str], label: str) -> list[str]:
    return [f"{label}: required artifact not recorded: {path}" for path in sorted(required - recorded_paths(records))]


def _cross_check_summary(root: Path, summary_relpath: str, results: dict[str, Any]) -> list[str]:
    path = root / summary_relpath
    if not path.is_file():
        return []
    try:
        summary = loads_strict(path.read_text(encoding="utf-8"))
        production = [row for row in summary["area"]["overall"] if row["estimator"] == PRODUCTION_ESTIMATOR_ID][0]
        strata = {row["rho_stratum"]: row["status"] for row in summary["volume"]["by_rho_stratum"]}
    except (EvidenceError, KeyError, IndexError, TypeError, ValueError) as error:
        return [f"derived summary {summary_relpath} is unreadable or not strict JSON: {error}"]
    problems = []
    recorded = results.get("SG-1", {}).get("max_abs_rel_error")
    if production.get("status") != results.get("SG-1", {}).get("status") or not (
            isinstance(recorded, float) and _relative_difference(float(production["max_abs"]), recorded) < 1e-8):
        problems.append(f"derived summary {summary_relpath} disagrees with the raw-result SG-1 derivation")
    derived = {row["rho_stratum"]: row["status"] for row in results.get("SG-2", {}).get("by_rho_stratum", [])}
    if strata != derived:
        problems.append(f"derived summary {summary_relpath} disagrees with the raw-result SG-2 strata")
    return problems


def verify_record(root: Path, manifest_relpath: str | None = None, *, check_git: bool = True,
                  reproduce: bool = True) -> VerificationReport:
    """Return every reason the recorded evidence no longer supports its results."""
    try:
        relpath = manifest_relpath or current_manifest_relpath(root)
        manifest = load_manifest(root / relpath)
    except EvidenceError as error:
        return VerificationReport(manifest_relpath, (str(error),))
    run_dir = str(PurePosixPath(relpath).parent)
    problems = verify_seal(manifest)
    if manifest.get("schema_version") != SCHEMA_VERSION:
        problems.append(f"manifest: schema {manifest.get('schema_version')!r} is not {SCHEMA_VERSION}")
    if manifest.get("contract_id") != CONTRACT_ID:
        problems.append(f"manifest: contract {manifest.get('contract_id')!r} is not {CONTRACT_ID}")
    if manifest.get("record_class") not in RECORD_CLASSES:
        problems.append(f"manifest: unknown record class {manifest.get('record_class')!r}")
    source = manifest.get("source")
    if not isinstance(source, dict) or source.get("scope_clean_at_execution") is not True:
        problems.append("source: the dependency scope was not clean at the validated commit when the V&V ran")

    scope = manifest.get("dependency_scope", {}).get("files") if isinstance(manifest.get("dependency_scope"), dict) else None
    problems += _check_required(scope, {entry.path for entry in DEPENDENCY_SCOPE}, "dependency scope")
    problems += verify_file_records(root, scope, "dependency scope")
    groups = {
        "vv_plan": {PLAN_PATH, FREEZE_PATH},
        "vv_implementation": {f"{run_dir}/{HARNESS_NAME}", *VV_TOOLING},
        "raw_artifacts": {f"{run_dir}/{name}" for name in RAW_NAMES},
        "derived_artifacts": {f"{run_dir}/{name}" for name in DERIVED_NAMES},
    }
    for group, required in groups.items():
        problems += _check_required(manifest.get(group), required, group)
        problems += verify_file_records(root, manifest.get(group), group)
    environment = _as_dict(manifest.get("environment"))
    lockfile = environment.get("lockfile")
    problems += _check_required([lockfile], {LOCKFILE}, "environment")
    problems += verify_file_records(root, [lockfile], "environment")

    frozen = read_freeze(root)
    for frozen_path, current_path in ((PLAN_PATH, PLAN_PATH), (FROZEN_HARNESS_PATH, f"{run_dir}/{HARNESS_NAME}")):
        target = root / current_path
        if frozen.get(frozen_path) is None or not target.is_file() or sha256_file(target) != frozen[frozen_path]:
            problems.append(f"frozen protocol: {current_path} is not the version frozen in {FREEZE_PATH}")

    if manifest.get("estimator") != production_estimator_identity():
        problems.append("estimator: the production estimator identity differs from the one evaluated")
    if manifest.get("acceptance_criteria") != ACCEPTANCE_CRITERIA:
        problems.append("criteria: recorded acceptance criteria differ from SCIENTIFIC_SPEC section 9 criteria")
    packages = _as_dict(environment.get("packages"))
    for name, version in installed_versions(CRITICAL_PACKAGES).items():
        if packages.get(name) != version:
            problems.append(f"environment: {name} {version} is installed, the record used {packages.get(name)}")
    if environment.get("python") != platform.python_version():
        problems.append(f"environment: Python {platform.python_version()} differs from recorded {environment.get('python')}")
    commands = manifest.get("commands")
    if not isinstance(commands, list) or not commands or any(
            not isinstance(c, dict) or c.get("exit_code") != 0 for c in commands):
        problems.append("commands: validation commands missing or not all successful")

    scope_hashes = {str(r.get("path")): r.get("sha256") for r in scope or [] if isinstance(r, dict)}
    for entry in DEPENDENCY_SCOPE:
        origin = module_origin(entry.path)
        if origin is not None and origin.is_file() and sha256_file(origin) != scope_hashes.get(entry.path):
            problems.append(f"runtime: imported module {origin} is not the recorded {entry.path}")

    raw_csv = root / run_dir / RAW_NAMES[0]
    results = manifest.get("results")
    if raw_csv.is_file():
        try:
            derived = derive_results(raw_csv)
        except EvidenceError as error:
            problems.append(str(error))
        else:
            if loads_strict(dumps_strict(derived)) != results:
                problems.append("results: recorded statuses/values are not those derived from the raw results")
            if manifest.get("resulting_status") != {item: derived[item]["status"] for item in derived}:
                problems.append("results: recorded resulting_status differs from the raw-result derivation")
            problems += _cross_check_summary(root, f"{run_dir}/{DERIVED_NAMES[0]}", derived)

    if check_git:
        problems += verify_git_provenance(root, source, scope)
    reproduction = None
    harness = root / run_dir / HARNESS_NAME
    if reproduce and raw_csv.is_file() and harness.is_file():
        try:
            reproduction_problems, reproduction = reproduce_production(harness, raw_csv)
        except EvidenceError as error:
            reproduction_problems = [f"re-execution: {error}"]
        problems += reproduction_problems
    return VerificationReport(relpath, tuple(problems), manifest, reproduction)
