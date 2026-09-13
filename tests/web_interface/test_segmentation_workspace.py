"""Regression tests for the Web upload path's auto-suggestion identity guard.

Guards against the P1-1 audit finding: an auto-detected nuclei diameter /
cell diameter / anisotropy / XY spacing computed for one uploaded stack must
never be silently reused for a later, different upload or a different
channel of the same TIFF -- doing so would also mislabel the stale spacing
as "from metadata" for a sample it was never measured on.
"""

from __future__ import annotations

import io
import json
import zipfile

import numpy as np
import pytest

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


@pytest.mark.parametrize("budget", [2000, 150])
def test_preview_budget_preserves_center_extent_and_raw_labels(monkeypatch, budget):
    from organoid_analysis.web_interface import segmentation_workspace as workspace

    monkeypatch.setattr(workspace, "_VIEWER_PAYLOAD_BUDGET_BYTES", budget)
    volume = np.broadcast_to(np.arange(71, dtype=np.uint16), (3, 53, 71)).copy()
    mask = np.zeros(volume.shape, np.uint32)
    mask[:, :, 35:] = 2**24 + 7
    originals = volume.copy(), mask.copy()
    images, labels, spacing = workspace._fit_viewer_payload_budget([volume], mask, (.2, .3, 1.))
    actual_bytes = np.prod(images[0].shape) * 4 * 2 * 1.34
    assert actual_bytes <= budget
    assert images[0].shape == labels.shape
    np.testing.assert_allclose((np.array(images[0].shape[::-1]) - 1) * spacing, [14., 15.6, 2.], rtol=1e-12)
    np.testing.assert_allclose(images[0][0, 0], np.linspace(0, 70, images[0].shape[2]), atol=1e-6)
    assert set(np.unique(labels)) == {0, 2**24 + 7}
    np.testing.assert_array_equal(volume, originals[0])
    np.testing.assert_array_equal(mask, originals[1])


def test_impossible_preview_budget_and_misaligned_channels_reject(monkeypatch):
    from organoid_analysis.web_interface import segmentation_workspace as workspace

    monkeypatch.setattr(workspace, "_VIEWER_PAYLOAD_BUDGET_BYTES", 10)
    with pytest.raises(ValueError, match="budget"):
        workspace._fit_viewer_payload_budget([np.ones((3, 8, 8))], None, (1, 1, 1))
    with pytest.raises(ValueError, match="one ZYX grid"):
        workspace._fit_viewer_payload_budget([np.ones((3, 8, 8))], np.ones((3, 7, 8)), (1, 1, 1))


def test_layer_and_geometry_selection_refresh_features_and_download(monkeypatch):
    from streamlit.testing.v1 import AppTest

    from organoid_analysis.web_interface import segmentation_workspace as workspace

    monkeypatch.setattr(workspace, "st_volume_viewer", lambda *args, **kwargs: None)
    app = AppTest.from_string('''
import numpy as np
import streamlit as st
from organoid_analysis.segmentation.cellpose_inference import SegmentationConfig
from organoid_analysis.web_interface.segmentation_workspace import render_results_tab, render_analysis_tab
if "nuclei_masks" not in st.session_state:
    nuclei = np.zeros((9, 9, 9), np.uint32)
    nuclei[2:5, 2:5, 2:5] = 10
    cells = np.zeros_like(nuclei)
    cells[1:8, 1:8, 1:8] = 70
    cells[3, 3, 3] = 0
    st.session_state["nuclei_masks"] = nuclei
    st.session_state["cell_masks"] = cells
render_results_tab(SegmentationConfig(xy_spacing_um=1., anisotropy=1., xy_spacing_source="user_override", anisotropy_source="user_override"))
render_analysis_tab()
''').run(timeout=30)
    assert not app.exception
    assert app.session_state["features"].volume_um3.tolist() == [27.]
    app.selectbox(key="feature_object_layer").select("Cells").run()
    assert not app.exception
    assert app.session_state["features"].volume_um3.tolist() == [342.]
    app.checkbox(key="feature_fill_holes").check().run()
    assert not app.exception
    assert app.session_state["features"].volume_um3.tolist() == [343.]
    with zipfile.ZipFile(io.BytesIO(app.session_state["feature_bundle"])) as archive:
        meta = json.loads(archive.read("measurement_provenance.json"))
    assert meta["object_type"] == "cell" and meta["measurement_basis"] == "filled_envelope"


def test_nuclear_only_ui_retains_primary_metrics_and_explains_surface_status(monkeypatch):
    from streamlit.testing.v1 import AppTest

    from organoid_analysis.web_interface import segmentation_workspace as workspace

    monkeypatch.setattr(workspace, "st_volume_viewer", lambda *args, **kwargs: None)
    app = AppTest.from_string('''
import numpy as np
import streamlit as st
from organoid_analysis.segmentation.cellpose_inference import SegmentationConfig
from organoid_analysis.web_interface.segmentation_workspace import render_results_tab
mask = np.zeros((5, 5, 5), np.uint32)
mask[2, 2, 2] = 71
st.session_state["nuclei_masks"] = mask
render_results_tab(SegmentationConfig(xy_spacing_um=1., anisotropy=2.))
''').run(timeout=30)
    assert not app.exception
    row = app.session_state["features"].loc[71]
    assert row.volume_um3 == 2.
    assert row.equivalent_diameter_um > 0 and row.least_axis_um > 0
    assert not row.surface_in_qualified_domain
    captions = " ".join(element.value for element in app.caption)
    assert "numerical resolution and anisotropy only" in captions
    assert "does not validate segmentation boundaries" in captions
