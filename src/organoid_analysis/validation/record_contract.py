"""Shared build/verify steps for qualification records keyed by gate item.

``analytical_geometry_evidence`` predates this module and keeps its own
characterization contract (SG-1b/SG-2b). The qualification contracts
(``analytical_qualification_evidence`` for SG-1a/SG-2a and
``viability_rule_evidence`` for SG-4A) share the checks below: seal, schema,
record class, source provenance, dependency-scope and tooling hashes,
environment, commands, and a ``results`` map keyed by gate item id
(owner decision D-14, docs/OWNER_DECISIONS.md).

A record is only as good as what its ``verify`` re-derives. Each contract
therefore adds its own re-derivation from raw artifacts and re-execution of the
production code; the checks here only establish that the record still refers
to the repository as it is.
"""
from __future__ import annotations

import datetime
import platform
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from organoid_analysis.validation.evidence_manifest import (
    SCHEMA_VERSION,
    EvidenceError,
    file_record,
    git_source_state,
    installed_versions,
    load_manifest,
    loads_strict,
    module_origin,
    recorded_paths,
    seal,
    sha256_file,
    verify_file_records,
    verify_git_provenance,
    verify_seal,
)

RECORD_CLASSES = ("canonical", "scope-clean")
LOCKFILE = "pixi.lock"
CRITICAL_PACKAGES = ("numpy", "scipy", "scikit-image", "pandas")
RECORD_STORAGE_NOTE = (
    "This record is stored in a commit made after validated_source_commit. The storing commit is "
    "deliberately not named inside the record (a commit cannot contain its own hash); find it with "
    "`git log --format=%H -1 -- <manifest path>`. validated_source_commit identifies the source under "
    "validation, never the commit that stores this record.")
QUALIFICATION_STATUSES = ("PASS", "FAIL", "VOID")


@dataclass(frozen=True)
class ScopeFile:
    path: str
    role: str
    rationale: str


@dataclass(frozen=True)
class ContractSpec:
    """What a qualification record of one contract must bind."""
    contract_id: str
    validation_scope: str
    items: tuple[str, ...]
    dependency_scope: tuple[ScopeFile, ...]
    # Frozen protocol files (plan, freeze record, harness, generators): hashed,
    # and required to equal the hashes the freeze record wrote before any result.
    protocol: tuple[str, ...]
    tooling: tuple[str, ...]
    raw_names: tuple[str, ...]
    pointer: str


@dataclass(frozen=True)
class VerificationReport:
    manifest_path: str | None
    problems: tuple[str, ...]
    manifest: dict[str, Any] | None = None
    reproduction: dict[str, Any] | None = field(default=None)

    @property
    def valid(self) -> bool:
        return self.manifest is not None and not self.problems


def utc_now() -> str:
    return datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def current_manifest_relpath(root: Path, spec: ContractSpec) -> str:
    pointer = root / spec.pointer
    if not pointer.is_file():
        raise EvidenceError(f"no {spec.contract_id} record ({spec.pointer} missing)")
    document = loads_strict(pointer.read_text(encoding="utf-8"))
    if not isinstance(document, dict) or document.get("contract_id") != spec.contract_id:
        raise EvidenceError(f"{spec.pointer} does not name a {spec.contract_id} record")
    manifest = document.get("manifest")
    if not isinstance(manifest, str) or not manifest.endswith("/evidence_manifest.json"):
        raise EvidenceError(f"{spec.pointer}: invalid manifest path {manifest!r}")
    return manifest


def refuse_dirty_critical(root: Path, spec: ContractSpec) -> tuple[Any, list[str]]:
    """Return the Git state and the evidence-critical paths that differ from HEAD."""
    state = git_source_state(root)
    critical = {entry.path for entry in spec.dependency_scope} | set(spec.protocol) | set(spec.tooling) | {LOCKFILE}
    return state, sorted(critical & set(state.dirty_paths))


def build_manifest(root: Path, spec: ContractSpec, *, run_id: str, state: Any, run_rel: str,
                   results: dict[str, Any], commands: list[dict[str, Any]], invocation: str,
                   extra: dict[str, Any], derived: Sequence[str] = ()) -> dict[str, Any]:
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "contract_id": spec.contract_id,
        "validation_run_id": run_id,
        "validation_scope": spec.validation_scope,
        "record_class": "canonical" if not state.dirty_paths else "scope-clean",
        "generated_utc": utc_now(),
        "source": {
            "validated_source_commit": state.head,
            "validated_source_tree": state.tree,
            "whole_tree_clean_at_execution": not state.dirty_paths,
            "dirty_paths_at_execution": list(state.dirty_paths),
            "scope_clean_at_execution": True,
            "record_storage": RECORD_STORAGE_NOTE,
        },
        "dependency_scope": {"files": [file_record(root, e.path, role=e.role, rationale=e.rationale)
                                       for e in spec.dependency_scope]},
        "environment": {
            "lockfile": file_record(root, LOCKFILE),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "packages": installed_versions(CRITICAL_PACKAGES),
        },
        "protocol": [file_record(root, path) for path in spec.protocol],
        "vv_implementation": [file_record(root, path) for path in spec.tooling],
        "raw_artifacts": [file_record(root, f"{run_rel}/{name}") for name in spec.raw_names],
        "derived_artifacts": [file_record(root, f"{run_rel}/{name}") for name in derived],
        "commands": commands,
        "invocation": invocation,
        "results": results,
        "resulting_status": {item: results[item]["status"] for item in spec.items},
        **extra,
    }
    return seal(manifest)


def verify_common(root: Path, spec: ContractSpec, relpath: str, *, check_git: bool,
                  frozen_hashes: Callable[[Path], dict[str, str]] | None = None
                  ) -> tuple[dict[str, Any] | None, list[str]]:
    """Checks every contract shares. Returns the manifest (None if unreadable) and the problems."""
    try:
        manifest = load_manifest(root / relpath)
    except EvidenceError as error:
        return None, [str(error)]
    run_dir = str(Path(relpath).parent)
    problems = verify_seal(manifest)
    if manifest.get("schema_version") != SCHEMA_VERSION:
        problems.append(f"manifest: schema {manifest.get('schema_version')!r} is not {SCHEMA_VERSION}")
    if manifest.get("contract_id") != spec.contract_id:
        problems.append(f"manifest: contract {manifest.get('contract_id')!r} is not {spec.contract_id}")
    if manifest.get("record_class") not in RECORD_CLASSES:
        problems.append(f"manifest: unknown record class {manifest.get('record_class')!r}")
    source = manifest.get("source")
    if not isinstance(source, dict) or source.get("scope_clean_at_execution") is not True:
        problems.append("source: the dependency scope was not clean at the validated commit when the V&V ran")
    if manifest.get("record_class") == "canonical" and (
            not isinstance(source, dict) or source.get("whole_tree_clean_at_execution") is not True):
        problems.append("source: a canonical record requires a clean whole tree at execution")

    scope = manifest.get("dependency_scope", {}).get("files") if isinstance(manifest.get("dependency_scope"), dict) else None
    groups: list[tuple[str, object, set[str]]] = [
        ("dependency scope", scope, {entry.path for entry in spec.dependency_scope}),
        ("protocol", manifest.get("protocol"), set(spec.protocol)),
        ("vv_implementation", manifest.get("vv_implementation"), set(spec.tooling)),
        ("raw_artifacts", manifest.get("raw_artifacts"), {f"{run_dir}/{name}" for name in spec.raw_names}),
    ]
    for label, records, required in groups:
        problems += [f"{label}: required artifact not recorded: {path}" for path in sorted(required - recorded_paths(records))]
        if required or records:
            problems += verify_file_records(root, records, label)
    derived = manifest.get("derived_artifacts")
    if isinstance(derived, list) and derived:
        problems += verify_file_records(root, derived, "derived_artifacts")
    raw_environment = manifest.get("environment")
    environment: dict[str, Any] = raw_environment if isinstance(raw_environment, dict) else {}
    problems += verify_file_records(root, [environment.get("lockfile")], "environment")
    if recorded_paths([environment.get("lockfile")]) != {LOCKFILE}:
        problems.append(f"environment: lockfile {LOCKFILE} not recorded")
    raw_packages = environment.get("packages")
    packages: dict[str, Any] = raw_packages if isinstance(raw_packages, dict) else {}
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
    for entry in spec.dependency_scope:
        origin = module_origin(entry.path)
        if origin is not None and origin.is_file() and sha256_file(origin) != scope_hashes.get(entry.path):
            problems.append(f"runtime: imported module {origin} is not the recorded {entry.path}")
    if frozen_hashes is not None:
        frozen = frozen_hashes(root)
        for path in spec.protocol:
            target = root / path
            if path in frozen and (not target.is_file() or sha256_file(target) != frozen[path]):
                problems.append(f"frozen protocol: {path} is not the version frozen before execution")
    results = manifest.get("results")
    if not isinstance(results, dict) or set(results) != set(spec.items):
        problems.append(f"results: the record must report exactly {list(spec.items)}")
    elif any(not isinstance(results[item], dict) or results[item].get("status") not in QUALIFICATION_STATUSES
             for item in spec.items):
        problems.append("results: an item status is not one of " + ", ".join(QUALIFICATION_STATUSES))
    if check_git:
        problems += verify_git_provenance(root, source, list(scope or []) + [
            r for r in manifest.get("protocol") or [] if isinstance(r, dict)] + [
            r for r in manifest.get("vv_implementation") or [] if isinstance(r, dict)])
    return manifest, problems


def verify_referenced(manifest: dict[str, Any], relpath: str, spec: ContractSpec) -> set[str]:
    """Every repository file a record depends on (used to copy a record for testing)."""
    paths = {relpath, spec.pointer}
    for group in ("protocol", "vv_implementation", "raw_artifacts", "derived_artifacts"):
        paths |= recorded_paths(manifest.get(group))
    paths |= recorded_paths(manifest.get("dependency_scope", {}).get("files"))
    paths |= recorded_paths([manifest.get("environment", {}).get("lockfile")])
    return paths


def worst(errors: Iterable[tuple[str, float]]) -> tuple[str | None, float | None]:
    pairs = list(errors)
    if not pairs:
        return None, None
    case, value = max(pairs, key=lambda item: abs(item[1]))
    return case, abs(value)
