"""Per-object classification agreement metrics for the SG-4B viability study.

Infrastructure only. This module computes the statistics a predeclared,
independent SG-4B study would report against an orthogonal reference. It does
not select the reference, does not set acceptance criteria, and cannot by
itself move SG-4B from ``NOT ASSESSED``: the reference must be independent of
the production segmentation branch (see
``docs/INTENDED_USE_AND_ESTIMANDS.md`` section 5.5 Q1).

Why a confusion matrix and not a correlation
--------------------------------------------
The SG-4B claim is a four-state *per-object classification*. A rank correlation
between a continuous reference fraction and the state ordering can be satisfied
while many individual objects are mis-stated, and it yields no per-state error
rate (``docs/reviews/2026-09-13-solution-report-review.md`` M2). This module
therefore reports per-class sensitivity and precision with intervals, macro-F1,
and quadratic-weighted Cohen's kappa. Rank correlation is deliberately not
offered here; it is retained only for the well-level ATP comparison, whose
estimand really is a continuous fraction.

Cohen's kappa citation: Cohen, J. (1968). Weighted kappa: nominal scale
agreement with provision for scaled disagreement or partial credit.
*Psychological Bulletin*, 70(4), 213-220. https://doi.org/10.1037/h0026256.
The quadratic weight ``w_ij = ((i - j) / (k - 1))^2`` is the standard ordinal
weighting; class order is the caller-supplied ``classes`` order.
"""
from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np

# The four frozen viability states, in ordinal order from most viable-like to
# least. ``N_classifiable`` excludes ``indeterminate`` (see
# ``docs/INTENDED_USE_AND_ESTIMANDS.md`` section 5.3); callers that report the
# frozen V2 estimand pass only the three classifiable states.
CLASSIFIABLE_STATES: tuple[str, ...] = ("viable_like", "mixed_signal", "compromised_like")
INDETERMINATE_STATE = "indeterminate"
FROZEN_STATES: tuple[str, ...] = (*CLASSIFIABLE_STATES, INDETERMINATE_STATE)


def _as_index(
    labels: Sequence[str], reference: Sequence[str], classes: Sequence[str]
) -> tuple[np.ndarray, np.ndarray]:
    """Validate two equal-length label sequences and map them to class indices."""
    if len(classes) < 2:
        raise ValueError("agreement metrics require at least two classes")
    if len(set(classes)) != len(classes):
        raise ValueError("classes must be unique")
    if len(labels) != len(reference):
        raise ValueError(
            f"labels and reference must have the same length ({len(labels)} != {len(reference)})"
        )
    if len(labels) == 0:
        raise ValueError("agreement metrics require at least one object")
    index = {name: i for i, name in enumerate(classes)}
    unknown = sorted({value for value in (*labels, *reference) if value not in index})
    if unknown:
        raise ValueError(f"labels/reference contain unknown states: {', '.join(unknown)}")
    return (
        np.array([index[value] for value in labels], dtype=int),
        np.array([index[value] for value in reference], dtype=int),
    )


def confusion_matrix(
    labels: Sequence[str], reference: Sequence[str], classes: Sequence[str] = FROZEN_STATES
) -> np.ndarray:
    """Counts of ``reference`` (rows) versus ``labels`` (columns).

    Every object is counted exactly once; objects whose two sources disagree on
    topology or state appear in the off-diagonal cells rather than being
    dropped. Rows and columns follow ``classes`` order.
    """
    label_index, reference_index = _as_index(labels, reference, classes)
    counts = np.zeros((len(classes), len(classes)), dtype=int)
    np.add.at(counts, (reference_index, label_index), 1)
    return counts


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion.

    ``z`` defaults to the two-sided 95% normal quantile. The Wilson interval is
    used instead of the normal approximation because per-state sensitivity and
    precision are often near 0 or 1, where the normal interval produces bounds
    outside [0, 1]. Returns ``(nan, nan)`` when ``total == 0`` so a caller
    cannot read a 0/0 as a measured rate.
    """
    if total < 0 or successes < 0 or successes > total:
        raise ValueError("require 0 <= successes <= total")
    if total == 0:
        return (math.nan, math.nan)
    p = successes / total
    denominator = 1 + z * z / total
    centre = (p + z * z / (2 * total)) / denominator
    half = (z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total))) / denominator
    return (max(0.0, centre - half), min(1.0, centre + half))


def sensitivity_precision(
    counts: np.ndarray, classes: Sequence[str] = FROZEN_STATES
) -> dict[str, dict[str, float | None]]:
    """Per-class sensitivity, precision and support from a confusion matrix.

    Sensitivity is recall against the reference row; precision is the label
    column. Both carry a Wilson 95% interval. Classes with no reference or no
    predicted objects report ``None`` rather than a fabricated 0 or 1.
    """
    counts = np.asarray(counts)
    if counts.shape != (len(classes), len(classes)):
        raise ValueError(f"confusion matrix must be {len(classes)}x{len(classes)}")
    rows = counts.sum(axis=1)
    columns = counts.sum(axis=0)
    result: dict[str, dict[str, float | None]] = {}
    for i, name in enumerate(classes):
        true_positive = int(counts[i, i])
        sensitivity_ci = wilson_interval(true_positive, int(rows[i]))
        precision_ci = wilson_interval(true_positive, int(columns[i]))
        result[name] = {
            "support": int(rows[i]),
            "predicted": int(columns[i]),
            "sensitivity": (true_positive / rows[i]) if rows[i] else None,
            "sensitivity_ci95_low": sensitivity_ci[0] if rows[i] else None,
            "sensitivity_ci95_high": sensitivity_ci[1] if rows[i] else None,
            "precision": (true_positive / columns[i]) if columns[i] else None,
            "precision_ci95_low": precision_ci[0] if columns[i] else None,
            "precision_ci95_high": precision_ci[1] if columns[i] else None,
        }
    return result


def macro_f1(counts: np.ndarray, classes: Sequence[str] = FROZEN_STATES) -> float:
    """Unweighted mean of the per-class F1 scores; ``nan`` if any class is absent."""
    counts = np.asarray(counts)
    scores: list[float] = []
    for i in range(len(classes)):
        true_positive = int(counts[i, i])
        support = int(counts[i, :].sum())
        predicted = int(counts[:, i].sum())
        if support == 0 or predicted == 0:
            return math.nan
        precision = true_positive / predicted
        recall = true_positive / support
        scores.append(2 * precision * recall / (precision + recall) if precision + recall else 0.0)
    return sum(scores) / len(scores)


def quadratic_weighted_kappa(counts: np.ndarray) -> float:
    """Quadratic-weighted Cohen's kappa for an ordinal classification.

    ``nan`` when the expected disagreement is zero (a degenerate reference that
    is constant, or a single object), because kappa is undefined there rather
    than zero or one.
    """
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    n_classes = counts.shape[0]
    if counts.shape[0] != counts.shape[1]:
        raise ValueError("confusion matrix must be square")
    if n_classes < 2:
        raise ValueError("kappa requires at least two classes")
    if total == 0:
        raise ValueError("kappa requires at least one object")
    indices = np.arange(n_classes)
    weights = ((indices[:, None] - indices[None, :]) / (n_classes - 1)) ** 2
    observed = counts / total
    expected = np.outer(counts.sum(axis=1), counts.sum(axis=0)) / (total * total)
    expected_disagreement = float((weights * expected).sum())
    if expected_disagreement == 0.0:
        return math.nan
    observed_disagreement = float((weights * observed).sum())
    return 1 - observed_disagreement / expected_disagreement


def classification_agreement(
    labels: Sequence[str],
    reference: Sequence[str],
    classes: Sequence[str] = FROZEN_STATES,
) -> dict[str, object]:
    """Complete agreement summary of one production label set against a reference.

    Returns a strict-JSON-safe mapping (no ``NaN``/``Infinity``) suitable for a
    validation record or an export. Missing rates are ``None``.
    """

    def _finite(value: float) -> float | None:
        return None if math.isnan(value) else float(value)

    counts = confusion_matrix(labels, reference, classes)
    n_objects = int(counts.sum())
    correct = int(np.trace(counts))
    per_class = sensitivity_precision(counts, classes)
    return {
        "n_objects": n_objects,
        "n_classes": len(classes),
        "classes": list(classes),
        "accuracy": correct / n_objects,
        "macro_f1": _finite(macro_f1(counts, classes)),
        "quadratic_weighted_kappa": _finite(quadratic_weighted_kappa(counts)),
        "confusion_matrix_reference_by_label": counts.tolist(),
        "per_class": per_class,
    }


def inter_rater_agreement(
    rater_labels: Sequence[Sequence[str]], classes: Sequence[str] = FROZEN_STATES
) -> dict[str, object]:
    """Pairwise quadratic-weighted kappa among annotators, plus the mean.

    Each element of ``rater_labels`` is one annotator's label for the same
    ordered object set. Inter-rater agreement is reported rather than assumed
    away; it is a property of the reference and is not evidence that the
    production path agrees with it.
    """
    if len(rater_labels) < 2:
        raise ValueError("inter-rater agreement requires at least two annotators")
    length = len(rater_labels[0])
    if any(len(column) != length for column in rater_labels):
        raise ValueError("every annotator must label the same objects in the same order")
    pairwise: list[dict[str, object]] = []
    values: list[float] = []
    for i in range(len(rater_labels)):
        for j in range(i + 1, len(rater_labels)):
            kappa = quadratic_weighted_kappa(
                confusion_matrix(rater_labels[i], rater_labels[j], classes)
            )
            pairwise.append(
                {
                    "annotator_a": i,
                    "annotator_b": j,
                    "quadratic_weighted_kappa": None if math.isnan(kappa) else kappa,
                }
            )
            if not math.isnan(kappa):
                values.append(kappa)
    return {
        "n_annotators": len(rater_labels),
        "n_objects": length,
        "pairwise": pairwise,
        "mean_quadratic_weighted_kappa": (sum(values) / len(values)) if values else None,
    }
