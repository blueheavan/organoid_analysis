"""A non-finite scaled signal is refused with its own reason and no range flag."""
import math

import pandas as pd
import pytest

from organoid_analysis.phenotyping.viability import classify

CFG = {"mode": "controls", "min_control_replicates": 2, "min_control_separation_snr": 3.0,
       "high_gate": 0.6, "low_gate": 0.3}
CALIBRATION = pd.DataFrame([{"batch_id": "B", "status": "available", "reason": "ok",
                             "calcein_low": 0.0, "calcein_high": 1.0, "pi_low": 0.0, "pi_high": 1.0}])


@pytest.mark.parametrize("calcein,pi", [(-math.inf, 0.1), (math.nan, 5.0), (0.9, math.inf)])
def test_nonfinite_signal_carries_only_the_nonfinite_reason(calcein, pi):
    objects = pd.DataFrame([{"batch_id": "B", "morphology_eligible": True, "viability_measurement_eligible": True,
                             "viability_measurement_flags": "", "calcein_mean_bg_corrected": calcein,
                             "pi_mean_bg_corrected": pi}])
    row = classify(objects, CALIBRATION, CFG).iloc[0]
    assert row.viability_state == "indeterminate"
    assert row.viability_reason == "nonfinite_scaled_signal"


def test_finite_out_of_range_signal_is_still_flagged():
    objects = pd.DataFrame([{"batch_id": "B", "morphology_eligible": True, "viability_measurement_eligible": True,
                             "viability_measurement_flags": "", "calcein_mean_bg_corrected": 1.4,
                             "pi_mean_bg_corrected": -0.2}])
    row = classify(objects, CALIBRATION, CFG).iloc[0]
    assert row.viability_state == "viable_like"
    assert row.viability_reason == "high_calcein_low_PI;outside_control_range"
