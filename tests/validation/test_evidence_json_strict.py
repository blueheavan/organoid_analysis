"""Evidence JSON must be strict RFC 8259 JSON (no NaN, no Infinity)."""
from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import pytest

from organoid_analysis.validation.evidence_manifest import (
    EvidenceError,
    dumps_strict,
    loads_strict,
    sanitize_non_finite,
    sha256_bytes,
    sha256_file,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
EVIDENCE = PROJECT_ROOT / "docs" / "evidence"


@pytest.mark.parametrize("path", sorted(EVIDENCE.rglob("*.json")), ids=lambda p: p.relative_to(EVIDENCE).as_posix())
def test_every_evidence_json_file_is_strict(path):
    loads_strict(path.read_text(encoding="utf-8"))


def test_strict_parser_rejects_python_extensions():
    for token in ("NaN", "Infinity", "-Infinity"):
        with pytest.raises(EvidenceError, match="non-strict JSON constant"):
            loads_strict(f'{{"value": {token}}}')


def test_strict_serializer_refuses_non_finite_numbers():
    with pytest.raises(ValueError):
        dumps_strict({"max_abs": math.inf})


def test_sanitizer_is_lossless_and_explicit():
    clean, replaced = sanitize_non_finite({"max_abs": math.inf, "mean": math.nan, "n": 3,
                                           "rows": [{"err": -math.inf}, 1.5, math.nan]})
    assert clean == {"max_abs": None, "max_abs_non_finite": "Infinity", "mean": None, "mean_non_finite": "NaN",
                     "n": 3, "rows": [{"err": None, "err_non_finite": "-Infinity"}, 1.5, None]}
    assert replaced == [("/max_abs", "Infinity"), ("/mean", "NaN"), ("/rows/0/err", "-Infinity"), ("/rows/2", "NaN")]
    assert loads_strict(dumps_strict(clean)) == clean


def test_sanitizer_refuses_to_overwrite_an_existing_marker():
    with pytest.raises(EvidenceError, match="already exists"):
        sanitize_non_finite({"x": math.nan, "x_non_finite": "kept"})


def test_non_estimable_values_in_the_current_summary_are_explicit():
    pointer = loads_strict((EVIDENCE / "analytical_geometry_current_record.json").read_text(encoding="utf-8"))
    summary_path = PROJECT_ROOT / Path(pointer["manifest"]).parent / "surface_vv_dev_summary.json"
    summary = loads_strict(summary_path.read_text(encoding="utf-8"))
    rows = [row for row in summary["area"]["overall"] if row.get("max_abs_non_finite")]
    assert rows, "the development grid has non-estimable candidate cases"
    assert all(row["max_abs"] is None and row["n_nonestimable"] > 0 and row["status"] == "FAIL" for row in rows)


def test_historical_conversions_are_recorded_and_recoverable():
    record = loads_strict((EVIDENCE / "STRICT_JSON_CONVERSION.json").read_text(encoding="utf-8"))
    assert record["conversions"]
    for entry in record["conversions"]:
        assert sha256_file(PROJECT_ROOT / entry["path"]) == entry["converted_sha256"]
        assert entry["original_matches_that_commit"] is True
        if (PROJECT_ROOT / ".git").exists():
            commit_path = entry["original_retrievable_with"].removeprefix("git show ")
            original = subprocess.run(["git", "show", commit_path], cwd=PROJECT_ROOT, capture_output=True, check=True)
            assert sha256_bytes(original.stdout) == entry["original_sha256"]
            # Same values: strict null + marker exactly where the original had NaN/Infinity.
            converted = loads_strict((PROJECT_ROOT / entry["path"]).read_text(encoding="utf-8"))
            assert sanitize_non_finite(json.loads(original.stdout))[0] == converted
