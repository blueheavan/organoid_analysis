"""The scientific gate must report failures and never grant PASS without evidence."""
from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

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


def test_gate_exits_nonzero_while_the_surface_criterion_fails():
    completed = subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True, check=False)
    assert completed.returncode == 1
    assert "[FAIL] SG-1" in completed.stdout
    assert "GATE: NOT PASSED" in completed.stdout


def test_gate_evidence_matches_the_current_production_estimator():
    gate = _load_gate()
    sg1 = [item for item in gate.evaluate() if item.item_id == "SG-1"][0]
    # A stale-evidence FAIL would mean the estimator changed without new V&V.
    assert not sg1.basis.startswith("EVIDENCE STALE")


def test_pass_without_evidence_is_rejected(monkeypatch):
    gate = _load_gate()
    fabricated = gate.GateItem("SG-X", "claim", "criterion", "PASS", "no data", ())
    monkeypatch.setattr(gate, "STATIC_ITEMS", [fabricated])
    with pytest.raises(ValueError, match="PASS requires existing evidence"):
        gate.evaluate()
