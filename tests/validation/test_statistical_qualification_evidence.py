"""Adversarial tests for the SG-5 statistical qualification record."""
from __future__ import annotations

import csv
import io
from pathlib import Path

import pytest

from organoid_analysis.validation import statistical_qualification_evidence as sg5
from organoid_analysis.validation.evidence_manifest import (
    dumps_strict,
    load_manifest,
    seal,
    sha256_file,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
REFERENCE_OK = {"passed": True, "worst": {}, "tolerance": sg5.REFERENCE_RTOL, "boundary_mismatch": [],
                "versions": {}, "n_datasets": 10}


def _write(path: Path, rows: list[dict]) -> Path:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return path


def _cells(tmp_path: Path, type1: float = 0.05, coverage: float = 0.95, n: int = sg5.N_CONTRAST_CELLS):
    rows = [{"K": 2, "G": 6, "m": 5, "balance": "balanced", "icc": 0.0, "dist": "normal", "n_sim": sg5.N_SIM,
             "seed": sg5.SEED, "failures": 0, "type1": 0.05, "coverage": 0.95, "omnibus_type1": 0.05,
             "fallback_fraction": 0.1} for _ in range(n)]
    rows[0]["type1"], rows[0]["coverage"] = type1, coverage
    intervals = [{"G": 6, "dist": "normal", "n_sim": sg5.N_SIM, "seed": sg5.SEED, "covered": 19000,
                  "coverage": 0.95} for _ in range(sg5.N_INTERVAL_CELLS)]
    return _write(tmp_path / "c.csv", rows), _write(tmp_path / "i.csv", intervals)


def test_nominal_cells_pass(tmp_path):
    assert sg5.derive_results(*_cells(tmp_path), REFERENCE_OK)["SG-5"]["status"] == "PASS"


@pytest.mark.parametrize("type1,coverage", [(0.0605, 0.95), (0.039, 0.95), (0.05, 0.9395), (0.05, 0.9605)])
def test_one_cell_outside_tolerance_fails(tmp_path, type1, coverage):
    assert sg5.derive_results(*_cells(tmp_path, type1, coverage), REFERENCE_OK)["SG-5"]["status"] == "FAIL"


def test_missing_cells_or_failed_reference_is_void(tmp_path):
    assert sg5.derive_results(*_cells(tmp_path, n=143), REFERENCE_OK)["SG-5"]["status"] == "VOID"
    assert sg5.derive_results(*_cells(tmp_path), {**REFERENCE_OK, "passed": False})["SG-5"]["status"] == "VOID"


def test_reference_verification_passes_on_the_current_implementation():
    report = sg5.reference_check(PROJECT_ROOT)
    assert report["passed"], report


def _problems(root: Path, relpath: str, **options: bool) -> tuple[str, ...]:
    options.setdefault("check_git", False)
    options.setdefault("reproduce", False)
    return sg5.verify_record(root, relpath, **options).problems


def test_untouched_statistical_copy_verifies(statistical_copy):
    root, relpath = statistical_copy
    assert _problems(root, relpath) == ()


def test_edited_cell_is_rejected(statistical_copy):
    root, relpath = statistical_copy
    raw = root / Path(relpath).parent / sg5.RAW_NAMES[0]
    rows = list(csv.DictReader(raw.open(newline="")))
    rows[0]["type1"] = "0.07"
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    raw.write_text(buffer.getvalue(), encoding="utf-8")
    problems = _problems(root, relpath)
    assert any("changed after the validation run" in p for p in problems)
    assert any("not those derived" in p for p in problems)


@pytest.mark.parametrize("path", [sg5.PROTOCOL_PATH, sg5.HARNESS_PATH, sg5.REFERENCE_FILES[1]])
def test_protocol_or_reference_change_is_rejected(statistical_copy, path):
    root, relpath = statistical_copy
    with (root / path).open("a", encoding="utf-8") as handle:
        handle.write("\n")
    assert any(path in p for p in _problems(root, relpath))


def test_forged_counts_with_rehashed_record_are_caught_by_reexecution(statistical_copy):
    root, relpath = statistical_copy
    raw = root / Path(relpath).parent / sg5.RAW_NAMES[0]
    rows = list(csv.DictReader(raw.open(newline="")))
    worst = max(rows, key=lambda c: max(abs(float(c["type1"]) - 0.05), abs(float(c["coverage"]) - 0.95)))
    worst["null_reject"] = str(int(worst["null_reject"]) - 1)  # plausible, still inside tolerance
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    raw.write_text(buffer.getvalue(), encoding="utf-8")
    manifest = load_manifest(root / relpath)
    for record in manifest["raw_artifacts"]:
        record["sha256"] = sha256_file(root / record["path"])
    (root / relpath).write_text(dumps_strict(seal(manifest)), encoding="utf-8")
    assert any(p.startswith("re-execution") for p in _problems(root, relpath, reproduce=True))
