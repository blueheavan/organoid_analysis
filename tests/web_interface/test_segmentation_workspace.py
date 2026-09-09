"""Regression tests for the Web upload path's auto-suggestion identity guard.

Guards against the P1-1 audit finding: an auto-detected nuclei diameter /
cell diameter / anisotropy / XY spacing computed for one uploaded stack must
never be silently reused for a later, different upload or a different
channel of the same TIFF -- doing so would also mislabel the stale spacing
as "from metadata" for a sample it was never measured on.
"""

from __future__ import annotations

from organoid_analysis.web_interface.segmentation_workspace import resolve_auto_suggestion


def test_no_suggestion_and_no_identity_returns_empty():
    assert resolve_auto_suggestion(None, None) == {}


def test_suggestion_present_but_current_identity_unknown_is_ignored():
    suggestion = {"nuclei_diameter": 30.0, "source_digest": "digestA", "source_channel": 0}
    assert resolve_auto_suggestion(suggestion, None) == {}


def test_matching_identity_is_applied():
    suggestion = {"nuclei_diameter": 30.0, "source_digest": "digestA", "source_channel": 0}
    assert resolve_auto_suggestion(suggestion, ("digestA", 0)) == suggestion


def test_switching_to_a_different_upload_does_not_inherit_the_suggestion():
    """Sample A auto-detected, then sample B uploaded: B must not inherit A's
    xy_spacing_um / anisotropy / diameter suggestion."""
    suggestion_from_a = {
        "nuclei_diameter": 30.0,
        "cell_diameter": 45.0,
        "xy_spacing_um": 0.65,
        "anisotropy": 3.1,
        "source_digest": "digestA",
        "source_channel": 0,
    }
    resolved_for_b = resolve_auto_suggestion(suggestion_from_a, ("digestB", 0))
    assert resolved_for_b == {}
    assert "xy_spacing_um" not in resolved_for_b
    assert "anisotropy" not in resolved_for_b
    assert "nuclei_diameter" not in resolved_for_b
    assert "cell_diameter" not in resolved_for_b


def test_switching_channel_on_the_same_tiff_does_not_inherit_the_suggestion():
    """Same file (same digest), channel 0 auto-detected, then switched to
    channel 1: the old diameter suggestion must not be silently reused."""
    suggestion_from_channel_0 = {
        "nuclei_diameter": 30.0,
        "source_digest": "digestSame",
        "source_channel": 0,
    }
    resolved_for_channel_1 = resolve_auto_suggestion(suggestion_from_channel_0, ("digestSame", 1))
    assert resolved_for_channel_1 == {}
