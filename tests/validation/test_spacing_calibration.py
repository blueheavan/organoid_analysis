"""Contract tests for the SG-6A spacing-ratio/stability workflow."""
from __future__ import annotations

import math

import pytest

from organoid_analysis.validation.evidence_manifest import dumps_strict
from organoid_analysis.validation.spacing_calibration import (
    AxisScale,
    ratio_report,
    stability_verdict,
)


def test_axis_ratios_are_computed_from_physical_scales() -> None:
    measurement = AxisScale(session="d1", x_um=0.5, y_um=0.5, z_um=0.9)
    assert measurement.xy_ratio == pytest.approx(1.0)
    assert measurement.zxy_ratio == pytest.approx(1.8)


def test_invalid_measurements_are_refused() -> None:
    with pytest.raises(ValueError):
        AxisScale(session="", x_um=1, y_um=1, z_um=1)
    with pytest.raises(ValueError):
        AxisScale(session="d1", x_um=0, y_um=1, z_um=1)
    with pytest.raises(ValueError):
        AxisScale(session="d1", x_um=float("nan"), y_um=1, z_um=1)
    with pytest.raises(ValueError):
        AxisScale(session="d1", x_um=True, y_um=1, z_um=1)


def test_report_aggregates_over_sessions_not_fields() -> None:
    report = ratio_report([
        AxisScale(session="d1", x_um=0.5, y_um=0.5, z_um=0.9),
        AxisScale(session="d2", x_um=0.5, y_um=0.5, z_um=0.9),
    ])
    assert report["unit_of_uncertainty"] == "imaging session"
    assert report["n_sessions"] == 2
    assert report["xy"]["mean"] == pytest.approx(1.0)
    assert report["zxy"]["mean"] == pytest.approx(1.8)
    assert report["zxy"]["sd"] == pytest.approx(0.0)
    assert report["zxy"]["normal_ci95_low"] == pytest.approx(1.8)


def test_single_session_reports_no_dispersion_rather_than_zero() -> None:
    report = ratio_report([AxisScale(session="d1", x_um=0.5, y_um=0.5, z_um=0.9)])
    assert report["zxy"]["sd"] is None
    assert report["zxy"]["normal_ci95_low"] is None


def test_verdict_without_a_tolerance_is_not_assessed() -> None:
    verdict = stability_verdict(
        [AxisScale(session="d1", x_um=0.5, y_um=0.5, z_um=0.9)], ratio="zxy", tolerance=None
    )
    assert verdict["status"] == "NOT ASSESSED"
    assert verdict["worst_session_deviation_from_mean"] is None


def test_verdict_applies_only_the_supplied_tolerance() -> None:
    stable = [AxisScale(session=f"d{i}", x_um=0.5, y_um=0.5, z_um=0.9) for i in range(3)]
    assert stability_verdict(stable, ratio="zxy", tolerance=0.0025)["status"] == "PASS"
    drifting = [
        AxisScale(session="d1", x_um=0.5, y_um=0.5, z_um=0.9),
        AxisScale(session="d2", x_um=0.5, y_um=0.5, z_um=1.0),
    ]
    verdict = stability_verdict(drifting, ratio="zxy", tolerance=0.0025)
    assert verdict["status"] == "FAIL"
    assert verdict["worst_session_deviation_from_mean"] > 0.0025


def test_verdict_refuses_invalid_tolerance_and_ratio() -> None:
    sessions = [AxisScale(session="d1", x_um=0.5, y_um=0.5, z_um=0.9)]
    with pytest.raises(ValueError):
        stability_verdict(sessions, ratio="xy", tolerance=-1)
    with pytest.raises(ValueError):
        stability_verdict(sessions, ratio="zx", tolerance=0.1)


def test_report_and_verdict_are_strict_json_safe() -> None:
    sessions = [AxisScale(session="d1", x_um=0.5, y_um=0.51, z_um=0.9)]
    assert dumps_strict(ratio_report(sessions))
    assert dumps_strict(stability_verdict(sessions, ratio="xy", tolerance=0.01))
    assert not math.isnan(ratio_report(sessions)["xy"]["mean"])
