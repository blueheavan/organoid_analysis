"""Adversarial tests: stale or manipulated evidence must never remain valid.

Each test tampers with a byte-identical copy of the current analytical-geometry
record (fixture ``evidence_copy``) and asserts that verification names the
tampered input. They test the evidence mechanism, not today's FAIL result.
"""
from __future__ import annotations

import csv
import importlib.util
import subprocess
from collections.abc import Callable
from pathlib import Path, PurePosixPath
from types import ModuleType
from typing import Any

import numpy as np
import pytest

from organoid_analysis.quantification.features import geometry
from organoid_analysis.validation import analytical_geometry_evidence as contract
from organoid_analysis.validation.evidence_manifest import (
    dumps_strict,
    load_manifest,
    loads_strict,
    seal,
    sha256_file,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
FEATURES = "src/organoid_analysis/quantification/features.py"
GATE_REFERENCE_SPHERE_AREA_UM2 = 4048.66364  # the former single-number sentinel, re-measured
LEGACY_GATE_REFERENCE_SPHERE_AREA_UM2 = 4549.28466796875  # its value under the superseded estimator


def problems(root: Path, relpath: str, **options: bool) -> tuple[str, ...]:
    options.setdefault("check_git", False)
    options.setdefault("reproduce", False)
    return contract.verify_record(root, relpath, **options).problems


def replace_once(path: Path, old: str, new: str) -> None:
    text = path.read_text(encoding="utf-8")
    assert text.count(old) == 1, f"mutation anchor not unique in {path}: {old!r}"
    path.write_text(text.replace(old, new), encoding="utf-8")


def edit_manifest(root: Path, relpath: str, change: Callable[[dict[str, Any]], None], *, reseal: bool = True) -> None:
    path = root / relpath
    manifest = load_manifest(path)
    change(manifest)
    path.write_text(dumps_strict(seal(manifest) if reseal else manifest), encoding="utf-8")


def rehash_everything(root: Path) -> Callable[[dict[str, Any]], None]:
    """Manifest edit that re-records every hash: a forger with full write access."""
    def change(manifest: dict[str, Any]) -> None:
        groups = [manifest[g] for g in ("vv_plan", "vv_implementation", "raw_artifacts", "derived_artifacts")]
        groups += [manifest["dependency_scope"]["files"], [manifest["environment"]["lockfile"]]]
        for records in groups:
            for record in records:
                record["sha256"] = sha256_file(root / record["path"])
    return change


def run_file(relpath: str, name: str) -> str:
    return f"{PurePosixPath(relpath).parent}/{name}"


def reference_sphere(inner_radius: float | None = None) -> np.ndarray:
    z, y, x = np.indices((25, 49, 49))
    distance2 = ((z - 12) * 2) ** 2 + (y - 24) ** 2 + (x - 24) ** 2
    mask = distance2 <= 18.0 ** 2
    return mask & (distance2 > inner_radius ** 2) if inner_radius else mask


def load_module_copy(path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location("organoid_analysis.quantification._mutant_features", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ------------------------------------------------------------------ baseline
def test_current_record_is_rejected_after_source_contract_changes(evidence_copy):
    root, relpath = evidence_copy
    # Make the change here rather than inheriting it from the repository's tree:
    # asserting rejection of an untouched copy would test today's tree state, not
    # the mechanism this name claims.
    with (root / FEATURES).open("a", encoding="utf-8") as handle:
        handle.write("\n# edit\n")
    report = problems(root, relpath)
    assert report
    assert any("changed after the validation run" in problem for problem in report)


def test_scope_verdict_names_exactly_the_differing_files(evidence_copy):
    """The scope verdict must name exactly the dependency-scope files whose content
    differs from the record, and none when the tree matches.

    The live repository currently matches its record, so the naming branch is
    exercised on a tampered copy inside this test. It is deliberately NOT asserted
    that the repository is stale: that is a transient property of the working tree
    and was exactly what made the previous version of this test fail once the
    record was regenerated.
    """
    # Branch 1 (naming): a tampered copy must be rejected, by name.
    root, relpath = evidence_copy
    tampered = contract.DEPENDENCY_SCOPE[0].path
    with (root / tampered).open("a", encoding="utf-8") as handle:
        handle.write("\n# edit\n")
    named = [p for p in problems(root, relpath) if "changed after the validation run" in p]
    assert len(named) == 1
    assert tampered in named[0]

    # Branch 2 (live): the same invariant against the repository's own record.
    if not (PROJECT_ROOT / ".git").exists():
        pytest.skip("Git provenance needs a Git checkout")
    manifest = load_manifest(PROJECT_ROOT / contract.current_manifest_relpath(PROJECT_ROOT))
    recorded = {entry["path"]: entry["sha256"] for entry in manifest["dependency_scope"]["files"]}
    differing = sorted(path for path, digest in recorded.items()
                       if sha256_file(PROJECT_ROOT / path) != digest)
    stale = [problem for problem in contract.verify_record(PROJECT_ROOT).problems
             if problem.startswith("dependency scope:") and "changed after the validation run" in problem]
    assert len(stale) == len(differing)
    for path in differing:
        assert any(path in problem for problem in stale)
    # Non-vacuous in either state: a clean tree names nothing, a stale tree must name.
    assert (differing == []) == (stale == [])


# ------------------------------------------------- 1, 2, 9: production code
def test_surface_change_without_method_version_bump_is_rejected(evidence_copy):
    root, relpath = evidence_copy
    # A change to the area actually exported, with the method_version string
    # left untouched -- the record must reject it on content, not on identity.
    replace_once(root / FEATURES, 'area = float(surface["surface_area"])',
                 'area = float(surface["surface_area"]) * 1.0001')
    assert '"method_version": surface_crofton.METHOD_NAME' in (root / FEATURES).read_text(encoding="utf-8")
    assert any(FEATURES in p and "changed" in p for p in problems(root, relpath))


def test_volume_change_is_rejected(evidence_copy):
    root, relpath = evidence_copy
    replace_once(root / FEATURES, "volume = float(count * np.prod(spacing_arr))",
                 "volume = float(count * np.prod(spacing_arr) * 1.0001)")
    assert any(FEATURES in p for p in problems(root, relpath))


@pytest.mark.parametrize("path", [entry.path for entry in contract.DEPENDENCY_SCOPE])
def test_any_change_inside_the_dependency_scope_is_rejected(evidence_copy, path):
    root, relpath = evidence_copy
    with (root / path).open("a", encoding="utf-8") as handle:
        handle.write("\n# edit\n")
    assert any(path in p for p in problems(root, relpath))


def test_change_invisible_to_the_old_sphere_sentinel_is_still_rejected(evidence_copy):
    root, relpath = evidence_copy
    replace_once(root / FEATURES, "envelope = outer_envelope(mask) if fill_holes else mask.astype(bool, copy=False)",
                 "envelope = outer_envelope(mask)")
    mutant = load_module_copy(root / FEATURES)
    spacing = (2.0, 1.0, 1.0)
    # The old gate re-measured only this sphere: the mutant reproduces it exactly...
    mutant_area = mutant.geometry(reference_sphere(), spacing)[0]["surface_area_um2"]
    assert mutant_area == geometry(reference_sphere(), spacing)[0]["surface_area_um2"]
    assert mutant_area == pytest.approx(GATE_REFERENCE_SPHERE_AREA_UM2, rel=1e-6)
    # ...although it silently replaced raw-label support by the filled envelope.
    shell = reference_sphere(inner_radius=10.0)
    assert (mutant.geometry(shell, spacing, fill_holes=False)[0]["volume_um3"]
            > geometry(shell, spacing, fill_holes=False)[0]["volume_um3"])
    assert any(FEATURES in p for p in problems(root, relpath))


def test_estimator_version_bump_does_not_inherit_old_evidence(evidence_copy, monkeypatch):
    root, relpath = evidence_copy
    bumped = {**contract.SURFACE_AREA_METHOD, "method_version": "marching_cubes_binary_lewiner_v2"}
    monkeypatch.setattr(contract, "SURFACE_AREA_METHOD", bumped)
    assert any(p.startswith("estimator:") for p in problems(root, relpath))


# --------------------------------------------------------- 3, 4: plan/script
def test_plan_change_is_rejected(evidence_copy):
    root, relpath = evidence_copy
    with (root / contract.PLAN_PATH).open("a", encoding="utf-8") as handle:
        handle.write("\nAcceptance: 20% is fine.\n")
    found = problems(root, relpath)
    assert any(contract.PLAN_PATH in p and "changed" in p for p in found)
    assert any(p.startswith("frozen protocol:") and contract.PLAN_PATH in p for p in found)


def test_vv_harness_change_is_rejected(evidence_copy):
    root, relpath = evidence_copy
    harness = run_file(relpath, contract.HARNESS_NAME)
    replace_once(root / harness, "CRITERION = 0.05", "CRITERION = 0.20")
    found = problems(root, relpath)
    assert any(harness in p and "changed" in p for p in found)
    assert any(p.startswith("frozen protocol:") for p in found)


@pytest.mark.parametrize("path", contract.VV_TOOLING)
def test_derivation_tooling_change_is_rejected(evidence_copy, path):
    root, relpath = evidence_copy
    with (root / path).open("a", encoding="utf-8") as handle:
        handle.write("\n# edit\n")
    assert any(path in p for p in problems(root, relpath))


# ---------------------------------------------------- 5, 6, 7: artifacts
def test_raw_result_change_is_rejected(evidence_copy):
    root, relpath = evidence_copy
    raw = root / run_file(relpath, contract.RAW_NAMES[0])
    with raw.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    first = next(row for row in rows if row["estimator"] == contract.PRODUCTION_ESTIMATOR_ID)
    first["area_est_um2"] = str(float(first["area_est_um2"]) * 1.01)
    with raw.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    assert any(contract.RAW_NAMES[0] in p and "changed" in p for p in problems(root, relpath))


def test_summary_only_change_is_rejected(evidence_copy):
    root, relpath = evidence_copy
    summary = root / run_file(relpath, contract.DERIVED_NAMES[0])
    document = loads_strict(summary.read_text(encoding="utf-8"))
    production = next(r for r in document["area"]["overall"] if r["estimator"] == contract.PRODUCTION_ESTIMATOR_ID)
    production["status"], production["max_abs"] = "PASS", 0.01
    summary.write_text(dumps_strict(document), encoding="utf-8")
    found = problems(root, relpath)
    assert any(contract.DERIVED_NAMES[0] in p and "changed" in p for p in found)
    assert any("disagrees with the raw-result SG-1" in p for p in found)


@pytest.mark.parametrize("name", [*contract.RAW_NAMES, *contract.DERIVED_NAMES, contract.HARNESS_NAME])
def test_removed_run_artifact_is_rejected(evidence_copy, name):
    root, relpath = evidence_copy
    (root / run_file(relpath, name)).unlink()
    assert any(name in p for p in problems(root, relpath))


@pytest.mark.parametrize("path", [contract.PLAN_PATH, contract.FREEZE_PATH, contract.LOCKFILE])
def test_removed_protocol_or_environment_file_is_rejected(evidence_copy, path):
    root, relpath = evidence_copy
    (root / path).unlink()
    assert any(path in p for p in problems(root, relpath))


def test_artifact_dropped_from_the_manifest_is_rejected(evidence_copy):
    root, relpath = evidence_copy
    edit_manifest(root, relpath, lambda m: m["raw_artifacts"].pop(0))
    assert any("required artifact not recorded" in p and contract.RAW_NAMES[0] in p for p in problems(root, relpath))


def test_shrinking_the_dependency_scope_is_rejected(evidence_copy):
    root, relpath = evidence_copy
    edit_manifest(root, relpath, lambda m: m["dependency_scope"].update(
        files=[r for r in m["dependency_scope"]["files"] if r["path"] != FEATURES]))
    assert any("required artifact not recorded" in p and FEATURES in p for p in problems(root, relpath))


def test_missing_record_pointer_means_no_evidence(evidence_copy):
    root, _ = evidence_copy
    (root / contract.CURRENT_RECORD_POINTER).unlink()
    report = contract.verify_record(root, check_git=False, reproduce=False)
    assert not report.valid and "no analytical-geometry validation record" in report.problems[0]


# ------------------------------------------------------- 8: fabricated PASS
def test_hand_edited_pass_breaks_the_seal_and_the_raw_derivation(evidence_copy):
    root, relpath = evidence_copy

    def fabricate(manifest: dict[str, Any]) -> None:
        manifest["results"]["SG-1"]["status"] = "PASS"
        manifest["resulting_status"]["SG-1"] = "PASS"
    edit_manifest(root, relpath, fabricate, reseal=False)
    found = problems(root, relpath)
    assert "manifest: content was edited after sealing (content seal mismatch)" in found
    assert any(p.startswith("results:") for p in found)
    edit_manifest(root, relpath, lambda m: None)  # re-sealed by the forger
    assert any(p.startswith("results:") for p in problems(root, relpath))


def test_relaxed_criterion_in_the_record_is_not_honored(evidence_copy):
    root, relpath = evidence_copy
    edit_manifest(root, relpath, lambda m: m["acceptance_criteria"]["SG-1"].update(threshold=0.2))
    assert any(p.startswith("criteria:") for p in problems(root, relpath))


def test_consistently_forged_pass_is_caught_by_reexecution(evidence_copy):
    root, relpath = evidence_copy
    raw = root / run_file(relpath, contract.RAW_NAMES[0])
    with raw.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        if row["estimator"] == contract.PRODUCTION_ESTIMATOR_ID:
            row["area_est_um2"], row["area_rel_err"] = row["area_true_um2"], "0"
    with raw.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = root / run_file(relpath, contract.DERIVED_NAMES[0])
    document = loads_strict(summary.read_text(encoding="utf-8"))
    production = next(r for r in document["area"]["overall"] if r["estimator"] == contract.PRODUCTION_ESTIMATOR_ID)
    production["status"], production["max_abs"] = "PASS", 0.0
    summary.write_text(dumps_strict(document), encoding="utf-8")
    derived = contract.derive_results(raw)
    assert derived["SG-1"]["status"] == "PASS"

    def forge(manifest: dict[str, Any]) -> None:
        rehash_everything(root)(manifest)
        manifest["results"] = loads_strict(dumps_strict(derived))
        manifest["resulting_status"] = {item: result["status"] for item, result in derived.items()}
    edit_manifest(root, relpath, forge)
    # Hashes, seal and derivation are all consistent, so they alone cannot tell...
    assert problems(root, relpath) == ()
    # ...but re-executing the production estimator on the phantoms does.
    assert any(p.startswith("re-execution:") for p in problems(root, relpath, reproduce=True))


# ------------------------------------------------------ environment, JSON
def test_lockfile_change_is_rejected(evidence_copy):
    root, relpath = evidence_copy
    with (root / contract.LOCKFILE).open("a", encoding="utf-8") as handle:
        handle.write("\n")
    assert any(contract.LOCKFILE in p and "changed" in p for p in problems(root, relpath))


def test_different_numerical_package_versions_are_rejected(evidence_copy, monkeypatch):
    root, relpath = evidence_copy
    monkeypatch.setattr(contract, "installed_versions", lambda names: {name: "0.0.0" for name in names})
    assert any(p.startswith("environment: numpy 0.0.0") for p in problems(root, relpath))


def test_failed_validation_command_is_rejected(evidence_copy):
    root, relpath = evidence_copy
    edit_manifest(root, relpath, lambda m: m["commands"][0].update(exit_code=1))
    assert any(p.startswith("commands:") for p in problems(root, relpath))


def test_record_from_a_dirty_scope_is_rejected(evidence_copy):
    root, relpath = evidence_copy
    edit_manifest(root, relpath, lambda m: m["source"].update(scope_clean_at_execution=False))
    assert any(p.startswith("source:") for p in problems(root, relpath))


def test_non_strict_json_manifest_is_rejected(evidence_copy):
    root, relpath = evidence_copy
    replace_once(root / relpath, '"rtol": 1e-09', '"rtol": NaN')
    found = problems(root, relpath)
    assert len(found) == 1 and "non-strict JSON constant 'NaN'" in found[0]


# ---------------------------------------------------------- git provenance
def _git(root: Path, *args: str) -> str:
    done = subprocess.run(["git", "-c", "user.name=evidence-test", "-c", "user.email=evidence-test@example.invalid",
                           *args], cwd=root, check=True, capture_output=True, text=True)
    return done.stdout.strip()


@pytest.fixture
def committed_copy(evidence_copy):
    root, relpath = evidence_copy
    edit_manifest(root, relpath, rehash_everything(root))
    _git(root, "init", "-q")
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", "source under validation")
    head, tree = _git(root, "rev-parse", "HEAD"), _git(root, "rev-parse", "HEAD^{tree}")
    edit_manifest(root, relpath, lambda m: m["source"].update(validated_source_commit=head, validated_source_tree=tree))
    return root, relpath, head


def test_hashes_matching_the_validated_commit_pass_git_provenance(committed_copy):
    root, relpath, _ = committed_copy
    assert problems(root, relpath, check_git=True) == ()


def test_nonexistent_source_commit_is_rejected(committed_copy):
    root, relpath, _ = committed_copy
    edit_manifest(root, relpath, lambda m: m["source"].update(validated_source_commit="0" * 40))
    assert any("is not in this repository" in p for p in problems(root, relpath, check_git=True))


def test_wrong_source_tree_is_rejected(committed_copy):
    root, relpath, head = committed_copy
    edit_manifest(root, relpath, lambda m: m["source"].update(validated_source_tree=head))
    assert any("is not the tree of commit" in p for p in problems(root, relpath, check_git=True))


def test_rehashed_source_that_is_not_the_validated_commit_is_rejected(committed_copy):
    root, relpath, head = committed_copy
    replace_once(root / FEATURES, "level=0.5,", "level=0.45,")
    _git(root, "commit", "-q", "-am", "later change")

    def rehash_scope(manifest: dict[str, Any]) -> None:
        for record in manifest["dependency_scope"]["files"]:
            record["sha256"] = sha256_file(root / record["path"])
    edit_manifest(root, relpath, rehash_scope)  # validated_source_commit still names `head`
    found = problems(root, relpath, check_git=True)
    assert any(f"{FEATURES} is not its content at validated_source_commit {head[:12]}" in p for p in found)
