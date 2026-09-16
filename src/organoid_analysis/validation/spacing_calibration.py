"""SG-6A relative spacing-ratio verification workflow (infrastructure only).

SG-6A verifies the X:Y and Z:XY ratios of the acquired voxel size, because
those ratios -- not absolute scale -- decide ``rho_in`` and the anisotropy flag
that classify an object into the analytical domain
(``docs/SCIENTIFIC_VALIDATION_MASTER_PLAN.md`` section 7, and
``docs/evidence/2026-09-13-volume-audit/STUDY_PROPOSALS_SG3_SG6.md`` section
2.1.1: ``rho_in`` is invariant to an isotropic scale error).

This module reports the ratios and their between-session stability. It defines
no acceptance tolerance, no required number of sessions, and no reference
standard: those are owner decisions. ``stability_verdict`` returns
``NOT ASSESSED`` unless a tolerance is supplied, and even a ``PASS`` here is a
computation, not the gate result -- SG-6A remains ``NOT ASSESSED`` until a
predeclared, independent, canonical validation record exists.

Absolute physical calibration and channel registration are SG-6B and are a
separate question; this module does not decide them.
"""
from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from statistics import NormalDist
from typing import TypedDict

RATIOS = ("xy", "zxy")

_NORMAL_975 = NormalDist().inv_cdf(0.975)


class Aggregate(TypedDict):
    n_sessions: int
    mean: float
    sd: float | None
    cv: float | None
    normal_ci95_low: float | None
    normal_ci95_high: float | None


@dataclass(frozen=True)
class AxisScale:
    """One imaging session's measured physical voxel size per axis, in micrometres."""

    session: str
    x_um: float
    y_um: float
    z_um: float

    def __post_init__(self) -> None:
        if not self.session:
            raise ValueError("session identifier must not be empty")
        for name, value in (("x_um", self.x_um), ("y_um", self.y_um), ("z_um", self.z_um)):
            if not isinstance(value, (int, float)) or isinstance(value, bool):
                raise ValueError(f"{name} must be a real number")
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive, got {value!r}")

    @property
    def xy_ratio(self) -> float:
        """X:Y scale ratio; 1.0 means square pixels."""
        return self.x_um / self.y_um

    @property
    def zxy_ratio(self) -> float:
        """Z divided by the geometric mean of X and Y; the anisotropy ratio."""
        return self.z_um / math.sqrt(self.x_um * self.y_um)


def _ratio_values(measurements: Sequence[AxisScale], ratio: str) -> list[float]:
    if ratio not in RATIOS:
        raise ValueError(f"ratio must be one of {RATIOS}, got {ratio!r}")
    return [getattr(measurement, f"{ratio}_ratio") for measurement in measurements]


def _aggregate(values: Sequence[float]) -> Aggregate:
    n = len(values)
    if n == 0:
        raise ValueError("at least one session is required")
    mean = sum(values) / n
    if n == 1:
        return {
            "n_sessions": 1,
            "mean": mean,
            "sd": None,
            "cv": None,
            "normal_ci95_low": None,
            "normal_ci95_high": None,
        }
    sd = math.sqrt(sum((value - mean) ** 2 for value in values) / (n - 1))
    half = _NORMAL_975 * sd / math.sqrt(n)
    return {
        "n_sessions": n,
        "mean": mean,
        "sd": sd,
        "cv": (sd / mean) if mean else None,
        "normal_ci95_low": mean - half,
        "normal_ci95_high": mean + half,
    }


def ratio_report(measurements: Sequence[AxisScale]) -> dict[str, object]:
    """Per-session ratios plus between-session mean, SD, CV and a 95% interval.

    Returns a strict-JSON-safe mapping. The interval is a **normal-approximation**
    interval on the session mean, not a Student-t interval: this keeps the
    module free of an incompletely-typed scientific import, and the study record
    must state the interval construction it uses at its own session count. The
    session, not the field of view, is the independent unit
    (``docs/SCIENTIFIC_VALIDATION_MASTER_PLAN.md`` section 7).
    """
    if not measurements:
        raise ValueError("at least one session is required")
    sessions = [
        {
            "session": measurement.session,
            "x_um": measurement.x_um,
            "y_um": measurement.y_um,
            "z_um": measurement.z_um,
            "xy_ratio": measurement.xy_ratio,
            "zxy_ratio": measurement.zxy_ratio,
        }
        for measurement in measurements
    ]
    return {
        "n_sessions": len(measurements),
        "sessions": sessions,
        "xy": _aggregate(_ratio_values(measurements, "xy")),
        "zxy": _aggregate(_ratio_values(measurements, "zxy")),
        "unit_of_uncertainty": "imaging session",
    }


def stability_verdict(
    measurements: Sequence[AxisScale], *, ratio: str, tolerance: float | None
) -> dict[str, object]:
    """Compare the between-session stability of ``ratio`` with ``tolerance``.

    ``tolerance`` is the owner's predeclared maximum acceptable deviation of the
    ratio from its session mean (relative, e.g. ``0.0025`` for 0.25%). It is a
    required argument and has no default: this function cannot invent one. With
    ``tolerance=None`` the verdict is ``NOT ASSESSED``. A supplied tolerance
    yields an engineering computation, not a gate result; SG-6A still requires a
    canonical study record.
    """
    if ratio not in RATIOS:
        raise ValueError(f"ratio must be one of {RATIOS}, got {ratio!r}")
    if tolerance is not None and (not math.isfinite(tolerance) or tolerance < 0):
        raise ValueError("tolerance must be finite and non-negative when supplied")
    aggregate = _aggregate(_ratio_values(measurements, ratio))
    if tolerance is None:
        status = "NOT ASSESSED"
        worst_deviation = None
    else:
        values = _ratio_values(measurements, ratio)
        mean = aggregate["mean"]
        worst_deviation = max(abs(value - mean) / mean for value in values)
        status = "PASS" if worst_deviation <= tolerance else "FAIL"
    return {
        "ratio": ratio,
        "tolerance_relative": tolerance,
        "worst_session_deviation_from_mean": worst_deviation,
        "aggregate": aggregate,
        "status": status,
        "basis": "computation only; SG-6A requires a predeclared canonical study record",
    }
