"""Adversarial tests for the item-keyed qualification records (SG-1a, SG-2a, SG-4A).

Two kinds: rule tests on synthetic raw results (they hold on any checkout), and
tamper tests on byte-identical copies of the committed records (skipped only on
a checkout that predates the records).
"""
from __future__ import annotations

import csv
import io
from pathlib import Path
from typing import Any

import pytest

from organoid_analysis.validation import analytical_qualification_evidence as qualification
from organoid_analysis.validation import viability_rule_evidence as viability
from organoid_analysis.validation.evidence_manifest import (
    dumps_strict,
    load_manifest,
    seal,
    sha256_file,
)

FIELDS = qualification.SMOOTH


# --------------------------------------------------------------- synthetic rules
def _row(case: str, grid: str, shape: str, *, pose: str = "sym-centre", area: float = 0.001, volume: float = 0.001,
         surface: bool = True, vol_domain: bool = True, weights: bool = True, anisotropy: float = 1.0,
         error: str = "") -> dict[str, Any]:
    return {"case_id": case, "grid": grid, "shape": shape, "smooth": shape in FIELDS, "pose": pose,
            "anisotropy": anisotropy, "area_rel_err": area, "volume_rel_err": volume,
            "surface_in_qualified_domain": surface, "volume_in_qualified_domain": vol_domain,
            "weights_evidence_bearing": weights, "error": error}


def _complete() -> list[dict[str, Any]]:
    rows = []
    for grid in ("confirm2", "confirm3"):
        for shape in FIELDS:
            for pose in ("sym-centre", "rand0"):
                rows.append(_row(f"{grid}-{shape}-{pose}", grid, shape, pose=pose))
    return rows


def _derive(tmp_path: Path, rows: list[dict[str, Any]]) -> dict[str, Any]:
    path = tmp_path / "cases.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return qualification.derive_results(path)


def test_complete_in_domain_set_within_criteria_passes(tmp_path):
    results = _derive(tmp_path, _complete())
    assert results["SG-1a"]["status"] == "PASS" and results["SG-2a"]["status"] == "PASS"


def test_one_volume_case_at_one_percent_fails_sg2a(tmp_path):
    rows = _complete()
    rows[-1]["volume_rel_err"] = -0.0100001
    assert _derive(tmp_path, rows)["SG-2a"]["status"] == "FAIL"


def test_confirm2_volume_never_counts_toward_sg2a(tmp_path):
    rows = _complete()
    rows[0]["volume_rel_err"] = 0.5  # confirm2: post hoc for the volume domain
    results = _derive(tmp_path, rows)
    assert results["SG-2a"]["status"] == "PASS"
    assert results["SG-2a"]["not_adjudicated"]["confirm2_in_volume_domain_post_hoc"]["max_abs_rel_error"] == 0.5


def test_out_of_domain_errors_do_not_fail_the_items(tmp_path):
    rows = _complete() + [_row("c3-out", "confirm3", "sphere", area=0.3, volume=0.3, surface=False, vol_domain=False),
                          _row("c3-crease", "confirm2", "box", area=0.3, volume=0.3)]
    results = _derive(tmp_path, rows)
    assert results["SG-1a"]["status"] == "PASS" and results["SG-2a"]["status"] == "PASS"


def test_non_evidence_bearing_weights_do_not_count_toward_sg1a(tmp_path):
    rows = _complete() + [_row("c3-solved", "confirm3", "sphere", area=0.2, weights=False)]
    results = _derive(tmp_path, rows)
    assert results["SG-1a"]["status"] == "PASS"
    assert results["SG-1a"]["excluded_non_evidence_bearing_weights"] == 1


@pytest.mark.parametrize("change", ["drop_class", "drop_random", "drop_isotropic", "failed_case"])
def test_incomplete_confirmation_is_void_never_pass(tmp_path, change):
    rows = _complete()
    if change == "drop_class":
        rows = [r for r in rows if r["shape"] != "torus"]
    elif change == "drop_random":
        rows = [r for r in rows if not (r["grid"] == "confirm3" and r["pose"] == "rand0")]
    elif change == "drop_isotropic":
        for r in rows:
            r["anisotropy"] = 2.0
    else:
        rows.append(_row("broken", "confirm3", "sphere", error="MemoryError"))
    results = _derive(tmp_path, rows)
    statuses = {results["SG-1a"]["status"], results["SG-2a"]["status"]}
    assert "VOID" in statuses and statuses <= {"VOID", "PASS"}


def test_viability_rules_void_when_a_required_reason_is_not_exercised():
    rows = [r for r in viability.run_cases() if "both_markers_low" not in r["reasons"]]
    result = viability.derive_results(rows)["SG-4A"]
    assert result["status"] == "VOID"
    assert "both_markers_low" in result["uncovered_required_reasons"]


def test_viability_oracle_detects_a_shifted_gate(monkeypatch):
    from organoid_analysis.phenotyping import viability as production

    original = production.classify

    def mutant(objects, calibration, cfg):
        return original(objects, calibration, {**cfg, "high_gate": cfg["high_gate"] - 1e-9})

    monkeypatch.setattr(production, "classify", mutant)
    result = viability.derive_results(viability.run_cases())["SG-4A"]
    assert result["status"] == "FAIL" and result["n_mismatched_cases"] > 0


def test_viability_oracle_detects_a_changed_denominator(monkeypatch):
    from organoid_analysis.statistics import aggregation

    monkeypatch.setattr(aggregation, "VIABILITY_FRACTION_DEFINITION", "all-eligible/1")
    assert viability.derive_results(viability.run_cases())["SG-4A"]["status"] == "FAIL"


# ------------------------------------------------------------- tamper: records
def _problems(module, root: Path, relpath: str, **options: bool) -> tuple[str, ...]:
    options.setdefault("check_git", False)
    options.setdefault("reproduce", False)
    return module.verify_record(root, relpath, **options).problems


def _raw(relpath: str, name: str) -> str:
    return f"{Path(relpath).parent}/{name}"


def test_untouched_qualification_copy_verifies(qualification_copy):
    root, relpath = qualification_copy
    assert _problems(qualification, root, relpath) == ()


def test_edited_raw_volume_is_rejected(qualification_copy):
    root, relpath = qualification_copy
    path = root / _raw(relpath, qualification.RAW_NAMES[0])
    rows = list(csv.DictReader(path.open(newline="")))
    target = next(r for r in rows if r["grid"] == "confirm3" and r["volume_in_qualified_domain"] == "True")
    target["volume_rel_err"] = "0.02"
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    path.write_text(buffer.getvalue(), encoding="utf-8")
    problems = _problems(qualification, root, relpath)
    assert any("changed after the validation run" in p for p in problems)
    assert any("not those derived from the raw results" in p for p in problems)


def test_hand_edited_pass_with_consistent_hashes_is_caught_by_derivation(qualification_copy):
    root, relpath = qualification_copy
    manifest = load_manifest(root / relpath)
    manifest["results"]["SG-2a"]["confirm3"]["max_abs_rel_error"] = 0.0
    (root / relpath).write_text(dumps_strict(seal(manifest)), encoding="utf-8")
    assert any("not those derived" in p for p in _problems(qualification, root, relpath))


def test_forged_raw_value_with_rehashed_record_is_caught_by_reexecution(qualification_copy, monkeypatch):
    root, relpath = qualification_copy
    # Only the worst adjudicated case of each item is re-executed here, to keep the test fast.
    monkeypatch.setattr(qualification, "REPRODUCE_WORST_PER_ITEM", 1)
    monkeypatch.setattr(qualification, "REPRODUCE_EVERY", 10**9)
    raw = root / _raw(relpath, qualification.RAW_NAMES[0])
    rows = list(csv.DictReader(raw.open(newline="")))
    worst = max((r for r in rows if r["grid"] == "confirm3" and r["volume_in_qualified_domain"] == "True"
                 and r["smooth"] == "True"), key=lambda r: abs(float(r["volume_rel_err"])))
    worst["volume_est"] = f"{float(worst['volume_est']) * 0.999:.10g}"  # plausible, still within 1%
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    raw.write_text(buffer.getvalue(), encoding="utf-8")
    manifest = load_manifest(root / relpath)
    for record in manifest["raw_artifacts"]:
        record["sha256"] = sha256_file(root / record["path"])
    (root / relpath).write_text(dumps_strict(seal(manifest)), encoding="utf-8")
    problems = _problems(qualification, root, relpath, reproduce=True)
    assert any(p.startswith("re-execution") for p in problems)


@pytest.mark.parametrize("path", [qualification.PROTOCOL_PATH, qualification.HARNESS_PATH,
                                  qualification.GENERATOR_FILES[1]])
def test_protocol_or_generator_change_is_rejected(qualification_copy, path):
    root, relpath = qualification_copy
    with (root / path).open("a", encoding="utf-8") as handle:
        handle.write("\n# edit\n")
    problems = _problems(qualification, root, relpath)
    assert any(path in p for p in problems)


def test_domain_constant_change_is_rejected(qualification_copy):
    root, relpath = qualification_copy
    path = root / "src/organoid_analysis/quantification/surface_crofton.py"
    text = path.read_text(encoding="utf-8")
    assert text.count("DOMAIN_RHO_VOL_MIN = 36.0") == 1
    path.write_text(text.replace("DOMAIN_RHO_VOL_MIN = 36.0", "DOMAIN_RHO_VOL_MIN = 20.0"), encoding="utf-8")
    assert any("surface_crofton.py changed" in p for p in _problems(qualification, root, relpath))


def test_untouched_viability_copy_verifies_and_reexecutes(viability_copy):
    root, relpath = viability_copy
    assert _problems(viability, root, relpath, reproduce=True) == ()


def test_viability_production_change_is_caught_by_reexecution(viability_copy, monkeypatch):
    root, relpath = viability_copy
    from organoid_analysis.statistics import aggregation

    monkeypatch.setattr(aggregation, "VIABILITY_FRACTION_DEFINITION", "all-eligible/1")
    problems = _problems(viability, root, relpath, reproduce=True)
    assert any(p.startswith("re-execution") for p in problems)


def test_viability_raw_edit_is_rejected(viability_copy):
    root, relpath = viability_copy
    raw = root / _raw(relpath, viability.RAW_NAME)
    raw.write_text(raw.read_text(encoding="utf-8").replace(",True,", ",False,", 1), encoding="utf-8")
    problems = _problems(viability, root, relpath)
    assert any("changed after the validation run" in p for p in problems)
    assert any("not those derived" in p for p in problems)


def test_scope_clean_record_cannot_pass_the_gate(qualification_copy):
    import importlib.util
    import sys

    root, relpath = qualification_copy
    manifest = load_manifest(root / relpath)
    manifest["record_class"] = "scope-clean"
    (root / relpath).write_text(dumps_strict(seal(manifest)), encoding="utf-8")
    script = Path(__file__).resolve().parents[2] / "scripts" / "scientific_validation_gate.py"
    spec = importlib.util.spec_from_file_location("gate_under_test", script)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    items = {item.item_id: item for item in module.evaluate(root, check_git=False, reproduce=False)}
    for item_id in ("SG-1a", "SG-2a"):
        assert items[item_id].status != "PASS"
