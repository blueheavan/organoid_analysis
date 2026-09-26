"""SG-4A: analytical qualification of the Calcein/PI four-state rule.

Scope of the claim
------------------
SG-4A asks whether the software implements the frozen rule, its refusals and
its versioned denominators exactly (docs/INTENDED_USE_AND_ESTIMANDS.md 5.1-5.3,
docs/SCIENTIFIC_VALIDATION_MASTER_PLAN.md section 9). It says nothing about
whether a state is biologically correct; that is SG-4B and needs an independent
object-level reference (owner decision D-8). A PASS here is an analytical,
known-rule result only.

Method
------
Every case is executed through the production functions
(``phenotyping.viability.calibrate``/``classify`` and
``statistics.aggregation._describe``) and compared with an oracle written from
the specification text, not from the production code:

* **State grid**: scaled calcein ``c`` and PI ``p`` on a lattice containing both
  gates, one ulp either side of each gate, 0 and 1 and one ulp outside them,
  interior values, far out-of-range values and every non-finite value, under a
  unit calibration (so ``c``/``p`` are exactly the inputs) and a non-trivial one.
* **Refusals**: morphology QC failure, measurement QC failure, a batch absent
  from the calibration table, and every calibration refusal reason.
* **Calibration arithmetic**: control medians through the object -> unit ->
  biological-replicate hierarchy and the separation SNR with its noise floor.
* **Denominators**: known populations, including zero classifiable and zero
  eligible objects, and the exported definition identifier.

Independence is limited: the oracle's author also read the production code
(audit tier D, docs/SCIENTIFIC_GATE_REASSESSMENT_2026-09-26.md). The
specification's two unstated boundary conventions are fixed explicitly in
``SPEC_CONVENTIONS`` below rather than inferred from the implementation.
"""
from __future__ import annotations

import csv
import io
import math
import statistics
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd  # type: ignore[import-untyped]

from organoid_analysis.validation.evidence_manifest import EvidenceError, dumps_strict, loads_strict
from organoid_analysis.validation.record_contract import (
    ContractSpec,
    ScopeFile,
    VerificationReport,
    current_manifest_relpath,
    verify_common,
)

CONTRACT_ID = "viability-analytical/1"
ITEM = "SG-4A"
RAW_NAME = "viability_rule_cases.csv"
POINTER = "docs/evidence/viability_rule_current_record.json"
SPEC = ContractSpec(
    contract_id=CONTRACT_ID,
    validation_scope=("SG-4A analytical conformance of the four-state Calcein/PI rule, its refusal reasons, "
                      "its control calibration arithmetic and its versioned denominators to the frozen "
                      "specification. Not a biological validation (SG-4B)."),
    items=(ITEM,),
    dependency_scope=(
        ScopeFile("src/organoid_analysis/phenotyping/viability.py", "rule_core",
                  "calibrate() and classify(): control medians, separation SNR, refusal reasons, gates and states."),
        ScopeFile("src/organoid_analysis/statistics/aggregation.py", "denominator_core",
                  "_describe(): classifiable and indeterminate denominators and the fraction-definition id."),
        ScopeFile("src/organoid_analysis/config.py", "defaults",
                  "Default gates, replicate floor and separation threshold, and the low<high validation."),
    ),
    protocol=(),
    tooling=("src/organoid_analysis/validation/viability_rule_evidence.py",
             "src/organoid_analysis/validation/record_contract.py",
             "src/organoid_analysis/validation/evidence_manifest.py",
             "scripts/record_qualification_evidence.py"),
    raw_names=(RAW_NAME,),
    pointer=POINTER,
)

# Frozen specification values (INTENDED_USE_AND_ESTIMANDS.md 5.1).
LOW_GATE, HIGH_GATE = 0.30, 0.60
MIN_CONTROL_REPLICATES, MIN_SEPARATION_SNR = 2, 3.0
MAD_TO_SIGMA = 1.4826
FRACTION_DEFINITION = "classifiable-denominator/2"
SPEC_CONVENTIONS = {
    "both_markers_low": "c < low_gate and p < low_gate (strict); equality to the low gate is not 'low' for this "
                        "refusal, while viable_like/compromised_like accept equality (p <= low_gate, c <= low_gate)",
    "outside_control_range": "flag appended when a finite scaled value is < 0 or > 1; classification uses the "
                             "unclipped value, which is equivalent because both gates lie in [0, 1]",
}
REQUIRED_REASONS = {
    "morphology_QC_failed", "missing_batch_controls", "viability_mode_uncalibrated", "no_QC_eligible_controls",
    "controls_span_multiple_conditions", "insufficient_independent_control_replicates",
    "calcein_controls_not_separated", "pi_controls_not_separated", "nonfinite_scaled_signal", "both_markers_low",
    "high_calcein_low_PI", "high_PI_low_calcein", "intermediate_or_discordant_marker_signals",
    "control_scaled_rules_require_external_validation", "outside_control_range", "measurement_flag_passthrough",
}
CFG = {"mode": "controls", "min_control_replicates": MIN_CONTROL_REPLICATES,
       "min_control_separation_snr": MIN_SEPARATION_SNR, "high_gate": HIGH_GATE, "low_gate": LOW_GATE}


# ------------------------------------------------------------------- oracle
def oracle_state(c: float, p: float) -> tuple[str, str]:
    if not (math.isfinite(c) and math.isfinite(p)):
        return "indeterminate", "nonfinite_scaled_signal"
    if c < LOW_GATE and p < LOW_GATE:
        state, reason = "indeterminate", "both_markers_low"
    elif c >= HIGH_GATE and p <= LOW_GATE:
        state, reason = "viable_like", "high_calcein_low_PI"
    elif p >= HIGH_GATE and c <= LOW_GATE:
        state, reason = "compromised_like", "high_PI_low_calcein"
    else:
        state, reason = "mixed_signal", "intermediate_or_discordant_marker_signals"
    if not (0.0 <= c <= 1.0 and 0.0 <= p <= 1.0):
        reason += ";outside_control_range"
    return state, reason


def _oracle_mad(values: list[float]) -> float:
    centre = statistics.median(values)
    return MAD_TO_SIGMA * statistics.median([abs(v - centre) for v in values])


def oracle_calibration(controls: list[dict[str, Any]]) -> dict[str, Any]:
    """Expected calibration for one batch of QC-eligible live/dead control objects, by explicit loops."""
    fields = ("calcein_mean_bg_corrected", "pi_mean_bg_corrected",
              "calcein_background_noise_mad", "pi_background_noise_mad")
    if not controls:
        return {"status": "unavailable", "reason": "no_QC_eligible_controls"}
    for kind in ("live", "dead"):
        if len({o["condition"] for o in controls if o["control"] == kind}) > 1:
            return {"status": "unavailable", "reason": "controls_span_multiple_conditions"}
    units: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    for obj in controls:
        units.setdefault((obj["control"], obj["biological_replicate"], obj["unit_id"]), []).append(obj)
    unit_medians = {key: {f: statistics.median(o[f] for o in objs) for f in fields} for key, objs in units.items()}
    replicates: dict[tuple[str, str], list[dict[str, float]]] = {}
    for (kind, rep, _), med in unit_medians.items():
        replicates.setdefault((kind, rep), []).append(med)
    rep_medians = {key: {f: statistics.median(m[f] for m in meds) for f in fields} for key, meds in replicates.items()}
    by_kind = {kind: [m for (k, _), m in rep_medians.items() if k == kind] for kind in ("live", "dead")}
    if min(len(v) for v in by_kind.values()) < MIN_CONTROL_REPLICATES:
        return {"status": "unavailable", "reason": "insufficient_independent_control_replicates"}
    out: dict[str, Any] = {}
    reasons = []
    for marker, low_kind, high_kind in (("calcein", "dead", "live"), ("pi", "live", "dead")):
        lows = [m[f"{marker}_mean_bg_corrected"] for m in by_kind[low_kind]]
        highs = [m[f"{marker}_mean_bg_corrected"] for m in by_kind[high_kind]]
        noise = max(_oracle_mad(lows), _oracle_mad(highs),
                    statistics.median(m[f"{marker}_background_noise_mad"] for m in rep_medians.values()), 1e-9)
        low, high = statistics.median(lows), statistics.median(highs)
        snr = (high - low) / noise
        out.update({f"{marker}_low": low, f"{marker}_high": high, f"{marker}_separation_snr": snr})
        if not math.isfinite(snr) or snr < MIN_SEPARATION_SNR:
            reasons.append(f"{marker}_controls_not_separated")
    out["status"] = "unavailable" if reasons else "available"
    out["reason"] = ";".join(reasons) if reasons else "control_scaled_rules_require_external_validation"
    return out


def oracle_fractions(states: list[str], eligible: list[bool]) -> dict[str, Any]:
    kept = [s for s, e in zip(states, eligible) if e]
    classifiable = [s for s in kept if s != "indeterminate"]
    out: dict[str, Any] = {"n_included": len(kept), "n_classifiable": len(classifiable),
                           "fraction_indeterminate_denominator": len(kept),
                           "viability_fraction_definition": FRACTION_DEFINITION}
    for state in ("viable_like", "mixed_signal", "compromised_like"):
        n = sum(s == state for s in kept)
        out[f"n_{state}"] = n
        out[f"fraction_{state}"] = n / len(classifiable) if classifiable else math.nan
    n_ind = sum(s == "indeterminate" for s in kept)
    out["n_indeterminate"] = n_ind
    out["fraction_indeterminate"] = n_ind / len(kept) if kept else math.nan
    return out


# -------------------------------------------------------------------- cases
def _gate_lattice() -> list[float]:
    values = {-5.0, -0.5, 0.0, 0.15, 0.45, 0.8, 1.0, 1.5, 7.0, math.nan, math.inf, -math.inf}
    for anchor in (0.0, LOW_GATE, HIGH_GATE, 1.0):
        values |= {anchor, math.nextafter(anchor, -math.inf), math.nextafter(anchor, math.inf)}
    return sorted(values, key=lambda v: (math.isnan(v), v))


def _object(batch: str, index: int, **values: Any) -> dict[str, Any]:
    row: dict[str, Any] = {"batch_id": batch, "organoid_id": index, "condition": "treated", "control": "sample",
                           "biological_replicate": "R1", "unit_id": "W1", "morphology_eligible": True,
                           "viability_measurement_eligible": True, "viability_measurement_flags": "",
                           "calcein_mean_bg_corrected": 0.0, "pi_mean_bg_corrected": 0.0,
                           "calcein_background_noise_mad": 0.0, "pi_background_noise_mad": 0.0}
    row.update(values)
    return row


def _calibration_row(batch: str, calcein: tuple[float, float], pi: tuple[float, float]) -> dict[str, Any]:
    return {"batch_id": batch, "status": "available", "reason": "control_scaled_rules_require_external_validation",
            "calcein_low": calcein[0], "calcein_high": calcein[1], "pi_low": pi[0], "pi_high": pi[1]}


ControlSpec = dict[tuple[str, str], list[list[tuple[float, float, float, float]]]]


def _controls(spec: ControlSpec) -> list[dict[str, Any]]:
    """Control objects: spec maps (kind, replicate) -> list of (calcein, pi, calcein_noise, pi_noise) per unit."""
    rows: list[dict[str, Any]] = []
    for (kind, replicate), units in spec.items():
        for u, objects in enumerate(units):
            for k, (cal, pi, cn, pn) in enumerate(objects):
                rows.append(_object("B", 1000 + len(rows), control=kind, condition=f"{kind}_ctrl",
                                    biological_replicate=replicate, unit_id=f"{kind}{replicate}U{u}",
                                    calcein_mean_bg_corrected=cal, pi_mean_bg_corrected=pi,
                                    calcein_background_noise_mad=cn, pi_background_noise_mad=pn, organoid_id=k))
    return rows


SEPARATED: ControlSpec = {("live", "L1"): [[(900.0, 12.0, 5.0, 4.0), (1000.0, 10.0, 6.0, 3.0)], [(950.0, 11.0, 5.0, 4.0)]],
             ("live", "L2"): [[(1100.0, 14.0, 7.0, 5.0)]],
             ("dead", "D1"): [[(40.0, 700.0, 5.0, 6.0), (60.0, 800.0, 6.0, 7.0), (50.0, 900.0, 4.0, 5.0)]],
             ("dead", "D2"): [[(30.0, 850.0, 5.0, 6.0)], [(70.0, 650.0, 5.0, 6.0)]]}


def calibration_cases() -> Iterator[tuple[str, list[dict[str, Any]], dict[str, Any], dict[str, Any]]]:
    """(case, objects, cfg, oracle) for calibrate(); objects all in batch B."""
    good = _controls(SEPARATED)
    yield "calib-available", good, CFG, oracle_calibration(good)
    yield "calib-uncalibrated-mode", good, {**CFG, "mode": "uncalibrated"}, \
        {"status": "unavailable", "reason": "viability_mode_uncalibrated"}
    ineligible = [{**o, "morphology_eligible": False} for o in good]
    yield "calib-no-eligible-controls", ineligible, CFG, oracle_calibration([])
    measurement_ineligible = [{**o, "viability_measurement_eligible": False} for o in good]
    yield "calib-no-measurement-eligible-controls", measurement_ineligible, CFG, oracle_calibration([])
    split = [dict(o) for o in good]
    split[0]["condition"] = "other_arm"
    yield "calib-controls-span-conditions", split, CFG, oracle_calibration(split)
    single = _controls({key: value for key, value in SEPARATED.items() if key != ("live", "L2")})
    yield "calib-one-live-replicate", single, CFG, oracle_calibration(single)
    noisy_calcein = _controls({**SEPARATED, ("live", "L2"): [[(80.0, 14.0, 7.0, 5.0)]]})
    yield "calib-calcein-not-separated", noisy_calcein, CFG, oracle_calibration(noisy_calcein)
    noisy_pi = _controls({**SEPARATED, ("dead", "D2"): [[(30.0, 20.0, 5.0, 6.0)], [(70.0, 15.0, 5.0, 6.0)]]})
    yield "calib-pi-not-separated", noisy_pi, CFG, oracle_calibration(noisy_pi)
    high_noise = _controls({key: [[(c, p, 1e4, 1e4) for c, p, _, _ in unit] for unit in units]
                            for key, units in SEPARATED.items()})
    yield "calib-both-not-separated-by-noise-floor", high_noise, CFG, oracle_calibration(high_noise)
    with_samples = good + [_object("B", 1, calcein_mean_bg_corrected=1e6), _object("OTHER", 2, control="live")]
    yield "calib-ignores-samples-and-other-batches", with_samples, CFG, oracle_calibration(good)


def _run_calibration(objects: list[dict[str, Any]], cfg: dict[str, Any]) -> dict[str, Any]:
    from organoid_analysis.phenotyping.viability import calibrate
    row: dict[str, Any] = calibrate(pd.DataFrame(objects), ["B"], cfg).iloc[0].to_dict()
    return row


def _same(a: object, b: object) -> bool:
    if isinstance(a, (int, float, np.floating, np.integer)) and isinstance(b, (int, float, np.floating, np.integer)):
        a, b = float(a), float(b)
        return (math.isnan(a) and math.isnan(b)) or a == b or abs(a - b) <= 1e-12 * max(abs(a), abs(b))
    return a == b


def run_cases() -> list[dict[str, str]]:
    """Execute every case through production code; one row per compared quantity."""
    from organoid_analysis.config import DEFAULTS, validate_config
    from organoid_analysis.phenotyping.viability import classify
    from organoid_analysis.statistics.aggregation import _describe

    rows: list[dict[str, str]] = []

    def record(case: str, family: str, quantity: str, expected: object, observed: object, reasons: str = "") -> None:
        rows.append({"case": case, "family": family, "quantity": quantity, "expected": repr(expected),
                     "observed": repr(observed), "match": str(_same(expected, observed)), "reasons": reasons})

    # Defaults: the frozen values are what ships.
    viability: dict[str, Any] = dict(DEFAULTS["viability"])  # type: ignore[call-overload]
    for key, value in (("low_gate", LOW_GATE), ("high_gate", HIGH_GATE),
                       ("min_control_replicates", MIN_CONTROL_REPLICATES),
                       ("min_control_separation_snr", MIN_SEPARATION_SNR)):
        record(f"default-{key}", "defaults", key, value, viability[key])
    import copy
    inverted: dict[str, Any] = copy.deepcopy(DEFAULTS)
    inverted["viability"].update(low_gate=0.6, high_gate=0.3)
    try:
        validate_config(inverted)
        refused = "accepted"
    except ValueError:
        refused = "refused"
    record("config-inverted-gates", "defaults", "validation", "refused", refused)

    # State lattice under a unit and a non-trivial calibration.
    lattice = _gate_lattice()
    calibrations = {"unit": ((0.0, 1.0), (0.0, 1.0)), "affine": ((125.0, 2125.0), (-40.0, 360.0))}
    for name, (cal, pi) in calibrations.items():
        objects: list[dict[str, Any]] = []
        expected: list[tuple[str, float, float, str, str]] = []
        for i, c in enumerate(lattice):
            for j, p in enumerate(lattice):
                x = cal[0] + c * (cal[1] - cal[0])
                y = pi[0] + p * (pi[1] - pi[0])
                objects.append(_object("B", len(objects), calcein_mean_bg_corrected=x, pi_mean_bg_corrected=y))
                # The oracle receives the scaled value exactly as the specification defines it.
                cs = (x - cal[0]) / (cal[1] - cal[0])
                ps = (y - pi[0]) / (pi[1] - pi[0])
                expected.append((f"state-{name}-{i:02d}-{j:02d}", cs, ps, *oracle_state(cs, ps)))
        result = classify(pd.DataFrame(objects), pd.DataFrame([_calibration_row("B", cal, pi)]), CFG)
        for (case, cs, ps, state, reason), (_, row) in zip(expected, result.iterrows()):
            record(case, "state", "viability_state", state, row.viability_state, reason)
            record(case, "state", "viability_reason", reason, row.viability_reason, reason)
            record(case, "state", "calcein_control_scaled", cs, row.calcein_control_scaled)
            record(case, "state", "pi_control_scaled", ps, row.pi_control_scaled)

    # Object-level refusals.
    refusal_objects = [
        ("refuse-morphology", _object("B", 1, morphology_eligible=False, calcein_mean_bg_corrected=1.0),
         "morphology_QC_failed"),
        ("refuse-measurement", _object("B", 2, viability_measurement_eligible=False,
                                       viability_measurement_flags="saturated_calcein;missing_pi_channel"),
         "saturated_calcein;missing_pi_channel"),
        ("refuse-missing-batch", _object("NOBATCH", 3, calcein_mean_bg_corrected=1.0), "missing_batch_controls"),
        ("refuse-unavailable-batch", _object("U", 4, calcein_mean_bg_corrected=1.0), "pi_controls_not_separated"),
    ]
    calibration = pd.DataFrame([_calibration_row("B", (0.0, 1.0), (0.0, 1.0)),
                                {**_calibration_row("U", (0.0, 1.0), (0.0, 1.0)), "status": "unavailable",
                                 "reason": "pi_controls_not_separated"}])
    result = classify(pd.DataFrame([obj for _, obj, _ in refusal_objects]), calibration, CFG)
    for (case, _, reason), (_, row) in zip(refusal_objects, result.iterrows()):
        tag = "measurement_flag_passthrough" if case == "refuse-measurement" else reason
        record(case, "refusal", "viability_state", "indeterminate", row.viability_state, tag)
        record(case, "refusal", "viability_reason", reason, row.viability_reason, tag)
        record(case, "refusal", "viability_method", "unclassified", row.viability_method, tag)

    # Calibration arithmetic and refusals.
    for case, control_objects, cfg, calibration_expected in calibration_cases():
        observed = _run_calibration(control_objects, cfg)
        for key, value in calibration_expected.items():
            record(case, "calibration", key, value, observed.get(key), str(calibration_expected["reason"]))

    # Denominators.
    populations = {
        "den-mixed": (["viable_like"] * 5 + ["mixed_signal"] * 3 + ["compromised_like"] * 2 + ["indeterminate"] * 4
                      + ["indeterminate"] * 3, [True] * 14 + [False] * 3),
        "den-no-classifiable": (["indeterminate"] * 4, [True] * 4),
        "den-no-eligible": (["indeterminate"] * 3, [False] * 3),
        "den-ineligible-carry-states": (["viable_like", "viable_like", "compromised_like", "indeterminate"],
                                        [True, False, True, True]),
    }
    for case, (states, eligible) in populations.items():
        frame = pd.DataFrame({"morphology_eligible": eligible, "viability_state": states,
                              "volume_um3": 1.0, "surface_area_um2": 1.0, "sphericity": 1.0,
                              "equivalent_diameter_um": 1.0})
        observed = _describe(frame)
        for key, value in oracle_fractions(states, eligible).items():
            record(case, "denominator", key, value, observed.get(key))
    return rows


# ------------------------------------------------------------------ results
def to_csv(rows: list[dict[str, str]]) -> str:
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=["case", "family", "quantity", "expected", "observed", "match",
                                                "reasons"], lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def read_csv(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))
    except (OSError, csv.Error) as error:
        raise EvidenceError(f"raw results unreadable: {error}") from error


def derive_results(rows: list[dict[str, str]]) -> dict[str, Any]:
    mismatches = sorted({row["case"] for row in rows if row.get("match") != "True"})
    families = sorted({row["family"] for row in rows})
    observed_reasons: set[str] = set()
    for row in rows:
        for token in row.get("reasons", "").split(";"):
            if token:
                observed_reasons.add(token)
    uncovered = sorted(REQUIRED_REASONS - observed_reasons)
    required_families = ["calibration", "defaults", "denominator", "refusal", "state"]
    valid = families == required_families and not uncovered
    status = "VOID" if not valid else ("PASS" if not mismatches else "FAIL")
    return {ITEM: {"status": status, "criterion": "every known-rule case equals the specification oracle; every "
                                                  "required state and refusal reason is exercised",
                   "n_comparisons": len(rows), "n_cases": len({row["case"] for row in rows}),
                   "n_mismatched_cases": len(mismatches), "mismatched_cases": mismatches[:50],
                   "families": families, "uncovered_required_reasons": uncovered,
                   "fraction_definition": FRACTION_DEFINITION, "spec_conventions": SPEC_CONVENTIONS,
                   "claim_boundary": "analytical rule conformance only; biological validity is SG-4B"}}


def verify_record(root: Path, manifest_relpath: str | None = None, *, check_git: bool = True,
                  reproduce: bool = True) -> VerificationReport:
    try:
        relpath = manifest_relpath or current_manifest_relpath(root, SPEC)
    except EvidenceError as error:
        return VerificationReport(manifest_relpath, (str(error),))
    manifest, problems = verify_common(root, SPEC, relpath, check_git=check_git)
    if manifest is None:
        return VerificationReport(relpath, tuple(problems))
    raw = root / Path(relpath).parent / RAW_NAME
    reproduction = None
    if raw.is_file():
        try:
            rows = read_csv(raw)
            derived = derive_results(rows)
        except EvidenceError as error:
            problems.append(str(error))
        else:
            if loads_strict(dumps_strict(derived)) != manifest.get("results"):
                problems.append("results: recorded status/values are not those derived from the raw results")
            if reproduce:
                fresh = to_csv(run_cases())
                same = fresh == raw.read_text(encoding="utf-8")
                reproduction = {"method": "re-executed every case through the production code", "identical": same}
                if not same:
                    problems.append("re-execution: the production code no longer reproduces the recorded cases")
    return VerificationReport(relpath, tuple(problems), manifest, reproduction)
