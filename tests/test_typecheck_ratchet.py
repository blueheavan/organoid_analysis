"""The mypy ratchet must reject new errors and lock in reductions."""
from __future__ import annotations

import importlib.util
import json
import sys
from collections import Counter
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "scripts" / "typecheck_ratchet.py"

A = "src/pkg/a.py:10: error: Function is missing a type annotation  [no-untyped-def]"
A_MOVED = "src/pkg/a.py:57: error: Function is missing a type annotation  [no-untyped-def]"
B = 'src/pkg/b.py:3: error: Library stubs not installed for "pandas"  [import-untyped]'
C = 'src/pkg/c.py:8: error: Need type annotation for "rows"  [var-annotated]'


def _load():
    spec = importlib.util.spec_from_file_location("typecheck_ratchet", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _output(*lines: str) -> str:
    return "\n".join([*lines, f"Found {len(lines)} errors in 3 files (checked 9 source files)"]) + "\n"


def test_signatures_ignore_line_numbers_but_keep_multiplicity():
    ratchet = _load()
    assert ratchet.parse_mypy_output(_output(A, A_MOVED)) == Counter(
        {"src/pkg/a.py | no-untyped-def | Function is missing a type annotation": 2})
    assert ratchet.normalize("x.py", "no-redef", 'Name "f" already defined on line 45') == \
        'x.py | no-redef | Name "f" already defined on line N'


def test_unparseable_output_is_refused():
    ratchet = _load()
    with pytest.raises(ValueError, match="reported 2 errors"):
        ratchet.parse_mypy_output(A + "\nFound 2 errors in 1 file (checked 1 source file)\n")


def test_removing_one_error_and_adding_another_is_not_unchanged_debt():
    ratchet = _load()
    baseline = ratchet.parse_mypy_output(_output(A, B))
    new, fixed = ratchet.compare(baseline, ratchet.parse_mypy_output(_output(A, C)))
    assert sum(new.values()) == 1 and sum(fixed.values()) == 1


def test_an_extra_occurrence_of_a_baselined_signature_is_new():
    ratchet = _load()
    new, _ = ratchet.compare(ratchet.parse_mypy_output(_output(A)), ratchet.parse_mypy_output(_output(A, A_MOVED)))
    assert sum(new.values()) == 1


def test_classification_rules():
    ratchet = _load()
    classify = ratchet.classify
    assert classify('m.py | import-untyped | Library stubs not installed for "scipy"') == ratchet.STUBS_MISSING
    assert classify('m.py | import-untyped | Skipping analyzing "statsmodels": module is installed, but missing '
                    'library stubs or py.typed marker') == ratchet.LIBRARY_BOUNDARY
    assert classify('m.py | union-attr | Item "TiffFrame" of "TiffPage | TiffFrame" has no attribute "tags"') \
        == ratchet.LIBRARY_BOUNDARY
    assert classify("m.py | no-untyped-def | Function is missing a type annotation") == ratchet.INTERNAL


@pytest.fixture
def isolated(tmp_path, monkeypatch):
    ratchet = _load()
    monkeypatch.setattr(ratchet, "BASELINE", tmp_path / "baseline.json")
    monkeypatch.setattr(ratchet, "ROOT", tmp_path)
    ratchet.BASELINE.write_text(json.dumps(ratchet.baseline_document(ratchet.parse_mypy_output(_output(A, B)), "x")))

    def use(*lines: str):
        monkeypatch.setattr(ratchet, "run_mypy", lambda: (_output(*lines), "x"))
        return ratchet
    return use


def test_new_error_fails_and_update_refuses_to_grow_the_baseline(isolated):
    ratchet = isolated(A, B, C)
    before = ratchet.BASELINE.read_text()
    assert ratchet.main([]) == 1
    assert ratchet.main(["--update"]) == 1
    assert ratchet.BASELINE.read_text() == before


def test_unchanged_debt_passes(isolated):
    assert isolated(A_MOVED, B).main([]) == 0


def test_reduction_must_be_locked_in_and_then_cannot_return(isolated):
    ratchet = isolated(A)
    assert ratchet.main([]) == 1  # reduced but not locked in
    assert ratchet.main(["--update"]) == 0
    assert ratchet.load_baseline(ratchet.BASELINE) == ratchet.parse_mypy_output(_output(A))
    assert isolated(A, B).main([]) == 1  # the fixed error may not come back


def test_committed_baseline_is_internally_consistent():
    ratchet = _load()
    document = json.loads(ratchet.BASELINE.read_text(encoding="utf-8"))
    entries = document["signatures"]
    assert document["total_errors"] == sum(entry["count"] for entry in entries)
    assert all(entry["class"] == ratchet.classify(entry["signature"]) for entry in entries)
    assert sum(document["class_counts"].values()) == document["total_errors"]
