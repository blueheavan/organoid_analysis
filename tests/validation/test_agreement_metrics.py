"""Contract tests for the SG-4B agreement metrics.

These tests check the statistics the module computes, not any qualification:
the module has no thresholds and cannot move SG-4B.
"""
from __future__ import annotations

import math

import pytest

from organoid_analysis.validation.agreement_metrics import (
    CLASSIFIABLE_STATES,
    FROZEN_STATES,
    classification_agreement,
    confusion_matrix,
    inter_rater_agreement,
    macro_f1,
    quadratic_weighted_kappa,
    sensitivity_precision,
    wilson_interval,
)
from organoid_analysis.validation.evidence_manifest import dumps_strict


def test_confusion_matrix_is_reference_by_label_with_every_object_counted() -> None:
    reference = ["viable_like", "viable_like", "compromised_like", "indeterminate"]
    labels = ["viable_like", "compromised_like", "compromised_like", "indeterminate"]
    counts = confusion_matrix(labels, reference)
    assert counts.shape == (4, 4)
    assert int(counts.sum()) == len(reference)
    # Row 0 is reference viable_like: one correct, one mislabelled compromised.
    assert counts[0, 0] == 1
    assert counts[0, 2] == 1
    # Row 2 is reference compromised_like: the single object is correct.
    assert counts[2, 2] == 1
    assert entries(counts, "indeterminate", "indeterminate") == 1


def entries(counts, reference_state: str, label_state: str) -> int:
    return int(counts[FROZEN_STATES.index(reference_state), FROZEN_STATES.index(label_state)])


def test_perfect_agreement_kappa_and_rates() -> None:
    states = ["viable_like", "mixed_signal", "compromised_like", "indeterminate"]
    summary = classification_agreement(states, states)
    assert summary["accuracy"] == 1.0
    assert summary["quadratic_weighted_kappa"] == 1.0
    assert summary["macro_f1"] == 1.0


def test_quadratic_weighted_kappa_two_class_hand_computed() -> None:
    counts = [[2, 1], [0, 3]]
    assert quadratic_weighted_kappa(counts) == pytest.approx(2 / 3, abs=1e-12)


def test_kappa_is_nan_when_reference_is_constant() -> None:
    counts = [[3, 0], [0, 0]]
    assert math.isnan(quadratic_weighted_kappa(counts))


def test_macro_f1_is_nan_when_a_class_is_absent() -> None:
    counts = [[2, 1, 0], [0, 3, 0], [0, 0, 0]]
    assert math.isnan(macro_f1(counts))


def test_sensitivity_precision_uses_match_positions() -> None:
    reference = ["viable_like", "viable_like", "compromised_like"]
    labels = ["viable_like", "compromised_like", "viable_like"]
    rates = sensitivity_precision(confusion_matrix(labels, reference))
    assert rates["viable_like"]["sensitivity"] == 0.5
    assert rates["viable_like"]["precision"] == 0.5
    assert rates["compromised_like"]["sensitivity"] == 0.0
    assert rates["compromised_like"]["precision"] == 0.0


def test_absent_class_reports_none_not_a_fabricated_rate() -> None:
    counts = confusion_matrix(["viable_like"], ["viable_like"], CLASSIFIABLE_STATES)
    rates = sensitivity_precision(counts, CLASSIFIABLE_STATES)
    assert rates["compromised_like"]["support"] == 0
    assert rates["compromised_like"]["sensitivity"] is None
    assert rates["compromised_like"]["precision"] is None


def test_wilson_interval_stays_inside_unit_interval() -> None:
    low, high = wilson_interval(0, 10)
    assert low == 0.0
    assert 0.0 < high < 1.0
    low, high = wilson_interval(10, 10)
    assert high == pytest.approx(1.0)
    assert 0.0 < low < 1.0
    zero_low, zero_high = wilson_interval(0, 0)
    assert math.isnan(zero_low) and math.isnan(zero_high)
    with pytest.raises(ValueError):
        wilson_interval(2, 1)


def test_unknown_state_and_length_mismatch_are_refused() -> None:
    with pytest.raises(ValueError, match="unknown states"):
        confusion_matrix(["viable_like"], ["not_a_state"])
    with pytest.raises(ValueError, match="same length"):
        confusion_matrix(["viable_like"], ["viable_like", "mixed_signal"])


def test_summary_is_strict_json_safe_and_deterministic() -> None:
    reference = ["viable_like", "mixed_signal", "compromised_like", "indeterminate"]
    labels = ["viable_like", "mixed_signal", "viable_like", "indeterminate"]
    first = classification_agreement(labels, reference)
    second = classification_agreement(labels, reference)
    assert first == second
    # dumps_strict raises on NaN/Infinity, so this proves the export is clean.
    assert dumps_strict(first)


def test_inter_rater_agreement_reports_every_pair_and_bounds_the_mean() -> None:
    rater_a = ["viable_like", "mixed_signal", "compromised_like"]
    rater_b = ["viable_like", "mixed_signal", "mixed_signal"]
    rater_c = ["viable_like", "viable_like", "compromised_like"]
    summary = inter_rater_agreement([rater_a, rater_b, rater_c])
    assert summary["n_annotators"] == 3
    assert len(summary["pairwise"]) == 3
    assert summary["pairwise"][0]["annotator_a"] == 0
    assert summary["pairwise"][0]["annotator_b"] == 1
    mean = summary["mean_quadratic_weighted_kappa"]
    assert mean is not None and -1.0 <= mean <= 1.0
    with pytest.raises(ValueError):
        inter_rater_agreement([rater_a])


def test_inter_rater_requires_identical_object_order() -> None:
    with pytest.raises(ValueError, match="same objects"):
        inter_rater_agreement([["viable_like"], ["viable_like", "mixed_signal"]])
