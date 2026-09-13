"""The scientific gate must report failures and never grant PASS without qualifying evidence."""
from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

from organoid_analysis.validation import analytical_geometry_evidence as contract
from organoid_analysis.validation.evidence_manifest import dumps_strict, load_manifest, seal

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "scripts" / "scientific_validation_gate.py"


def _load_gate():
    spec = importlib.util.spec_from_file_location("scientific_validation_gate", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    # @dataclass resolves its defining module through sys.modules, so the
    # script must be registered before it executes.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _item(items, item_id):
    return next(item for item in items if item.item_id == item_id)


def test_gate_exits_nonzero_while_scientific_items_lack_pass_evidence():
    completed = subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True, check=False)
    assert completed.returncode == 1
    assert "GATE: NOT PASSED" in completed.stdout
    # Non-estimable values are reported as such, never as nan%/inf%.
    assert "nan%" not in completed.stdout.lower() and "inf%" not in completed.stdout.lower()


def test_current_analytical_evidence_is_accepted_and_backs_the_items():
    gate = _load_gate()
    items = gate.evaluate()
    for item_id in ("SG-1a", "SG-2a"):
        item = _item(items, item_id)
        assert item.status == "FAIL"
        assert item.basis.startswith("EVIDENCE REJECTED")
        assert item.record is None
    for item_id in ("SG-1b", "SG-2b"):
        item = _item(items, item_id)
        assert item.status == "NOT APPLICABLE"
        assert item.record is None


def test_tampered_evidence_turns_the_analytical_items_into_rejections(evidence_copy):
    root, _ = evidence_copy
    with (root / "src/organoid_analysis/quantification/features.py").open("a", encoding="utf-8") as handle:
        handle.write("\n# edit\n")
    items = _load_gate().evaluate(root, check_git=False, reproduce=False)
    for item_id in ("SG-1a", "SG-2a"):
        assert _item(items, item_id).status == "FAIL"
        assert _item(items, item_id).basis.startswith("EVIDENCE REJECTED")
        assert _item(items, item_id).record is None


def test_missing_record_is_a_failure_not_a_pass(evidence_copy):
    root, _ = evidence_copy
    (root / contract.CURRENT_RECORD_POINTER).unlink()
    items = _load_gate().evaluate(root, check_git=False, reproduce=False)
    assert _item(items, "SG-1a").status == "FAIL" and "no analytical-geometry validation record" in _item(items, "SG-1a").basis


def test_pass_without_any_evidence_is_rejected(monkeypatch):
    gate = _load_gate()
    monkeypatch.setattr(gate, "STATIC_ITEMS", [gate.GateItem("SG-X", "claim", "criterion", "PASS", "no data", ())])
    with pytest.raises(ValueError, match="PASS requires a verified validation record"):
        gate.evaluate()


def test_pass_citing_existing_files_but_no_record_is_rejected(monkeypatch):
    gate = _load_gate()
    fabricated = gate.GateItem("SG-3", "claim", "criterion", "PASS", "looks fine", ("README.md",))
    monkeypatch.setattr(gate, "STATIC_ITEMS", [fabricated])
    with pytest.raises(ValueError, match="PASS requires a verified validation record"):
        gate.evaluate()


def test_pass_borrowing_a_record_that_does_not_report_it_is_rejected(evidence_copy, monkeypatch):
    root, relpath = evidence_copy
    manifest = load_manifest(root / relpath)
    manifest["record_class"] = "canonical"
    (root / relpath).write_text(dumps_strict(seal(manifest)), encoding="utf-8")
    gate = _load_gate()
    borrowed = gate.GateItem("SG-3", "claim", "criterion", "PASS", "see record", (relpath,), relpath)
    monkeypatch.setattr(gate, "STATIC_ITEMS", [borrowed])
    report = _passing_report("canonical")
    monkeypatch.setattr(gate.analytical, "verify_record", lambda *args, **kwargs: report)
    with pytest.raises(ValueError, match="PASS is not what its validation record reports"):
        gate.evaluate(root, check_git=False, reproduce=False)


def _passing_report(record_class: str):
    real = contract.verify_record(PROJECT_ROOT, check_git=False, reproduce=False)
    manifest: dict[str, Any] = {**real.manifest, "record_class": record_class}
    manifest["results"] = {item: {**result, "status": "PASS", "by_rho_stratum": [
        {**row, "status": "PASS"} for row in result.get("by_rho_stratum", [])]} for item, result in real.manifest["results"].items()}
    return contract.VerificationReport(real.manifest_path, (), manifest)


def test_gate_can_pass_an_item_backed_by_a_qualifying_canonical_record(monkeypatch):
    gate = _load_gate()
    report = _passing_report("canonical")  # built before verify_record is replaced
    monkeypatch.setattr(gate.analytical, "verify_record", lambda *args, **kwargs: report)
    fabricated = gate.GateItem("SG-1", "claim", "criterion", "PASS", "see record", (contract.CURRENT_RECORD_POINTER,),
                               contract.current_manifest_relpath(PROJECT_ROOT))
    monkeypatch.setattr(gate, "STATIC_ITEMS", [fabricated])
    assert _item(gate.evaluate(), "SG-1").status == "PASS"


def test_pass_from_a_non_canonical_record_is_refused(monkeypatch):
    gate = _load_gate()
    report = _passing_report("scope-clean")
    monkeypatch.setattr(gate.analytical, "verify_record", lambda *args, **kwargs: report)
    fabricated = gate.GateItem("SG-1", "claim", "criterion", "PASS", "see record", (contract.CURRENT_RECORD_POINTER,),
                               contract.current_manifest_relpath(PROJECT_ROOT))
    monkeypatch.setattr(gate, "STATIC_ITEMS", [fabricated])
    with pytest.raises(ValueError, match="PASS requires a canonical record"):
        gate.evaluate()
