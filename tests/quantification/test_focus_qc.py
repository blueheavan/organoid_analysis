"""Z-coverage QC: focus/signal profiles flag objects whose axial extent may be cut by the scanned range."""
import numpy as np
import pytest
from scipy import ndimage as ndi

from organoid_analysis.quantification.focus_qc import (
    MIN_PLANES,
    focus_profile,
    object_z_coverage,
    stack_z_coverage,
)

RNG = np.random.default_rng(7)


def _textured_blob(shape, centre_yx, radius, focus_plane, depth_sigma=2.0, amplitude=1000.0):
    """A textured disc whose texture is sharpest at ``focus_plane`` and blurs with |z - focus_plane|."""
    nz, ny, nx = shape
    y, x = np.indices((ny, nx))
    disc = (y - centre_yx[0]) ** 2 + (x - centre_yx[1]) ** 2 < radius ** 2
    texture = RNG.random((ny, nx))
    stack = np.zeros(shape, np.float32)
    for z in range(nz):
        blur = 0.3 + abs(z - focus_plane) * 1.2
        plane = ndi.gaussian_filter(texture, blur) * disc
        weight = np.exp(-0.5 * ((z - focus_plane) / depth_sigma) ** 2)
        stack[z] = 100 + amplitude * weight * plane
    return stack, disc


def _labels_from(disc, planes, shape, label=1):
    labels = np.zeros(shape, np.int32)
    for z in planes:
        labels[z][disc] = label
    return labels


def test_focus_peak_inside_the_stack_is_not_flagged():
    shape = (15, 64, 64)
    stack, disc = _textured_blob(shape, (32, 32), 12, focus_plane=7)
    labels = _labels_from(disc, range(4, 11), shape)
    result = object_z_coverage(stack, labels)[1]
    assert result["z_focus_peak_plane"] == 7
    assert not result["z_focus_peak_at_edge"]
    assert not result["z_signal_at_edge"]
    assert not result["z_coverage_suspect"]
    assert result["z_coverage_flags"] == ""


def test_object_cut_by_the_top_of_the_stack_is_flagged_even_if_its_mask_stops_short():
    shape = (15, 64, 64)
    stack, disc = _textured_blob(shape, (32, 32), 12, focus_plane=0)
    labels = _labels_from(disc, range(2, 5), shape)   # mask does not touch plane 0
    result = object_z_coverage(stack, labels)[1]
    assert result["z_focus_peak_plane"] == 0
    assert result["z_focus_peak_at_edge"]
    assert result["z_signal_at_edge"] and result["z_edge_signal_ratio"] == pytest.approx(1.0)
    assert result["z_coverage_suspect"]
    assert "mask_touches_z_boundary" not in result["z_coverage_flags"]


def test_mask_touching_the_last_plane_is_flagged():
    shape = (12, 48, 48)
    stack, disc = _textured_blob(shape, (24, 24), 10, focus_plane=6)
    labels = _labels_from(disc, range(5, 12), shape)
    result = object_z_coverage(stack, labels)[1]
    assert result["z_coverage_suspect"]
    assert result["z_coverage_flags"] == "mask_touches_z_boundary"


def test_an_organoid_stacked_in_z_does_not_contaminate_the_other():
    """Two objects share one XY footprint at different depths; each keeps its own focus peak."""
    shape = (20, 64, 64)
    upper, disc = _textured_blob(shape, (32, 32), 12, focus_plane=5, depth_sigma=1.5)
    lower, _ = _textured_blob(shape, (32, 32), 12, focus_plane=19, depth_sigma=1.5, amplitude=3000.0)
    stack = upper + lower - 100
    labels = _labels_from(disc, range(3, 8), shape, label=1)
    labels += _labels_from(disc, range(16, 20), shape, label=2)
    results = object_z_coverage(stack, labels)
    assert not results[1]["z_focus_peak_at_edge"]
    assert results[1]["z_focus_peak_plane"] == 5
    assert results[2]["z_focus_peak_plane"] == 19 and results[2]["z_coverage_suspect"]


def test_too_few_planes_is_suspect_without_computing_profiles():
    shape = (MIN_PLANES - 1, 32, 32)
    labels = np.zeros(shape, np.int32)
    labels[:, 10:20, 10:20] = 4
    result = object_z_coverage(np.ones(shape, np.float32), labels)[4]
    assert result["z_coverage_suspect"] and result["z_coverage_flags"] == "too_few_planes"
    assert stack_z_coverage(np.ones(shape, np.float32))["stack_z_coverage_suspect"]


def test_flat_stack_has_no_focus_edge_flag():
    stack = np.full((9, 32, 32), 50.0, np.float32)
    labels = np.zeros(stack.shape, np.int32)
    labels[3:6, 10:20, 10:20] = 1
    assert np.ptp(focus_profile(stack)) == 0
    result = object_z_coverage(stack, labels)[1]
    assert not result["z_focus_peak_at_edge"] and not result["z_coverage_suspect"]


def test_field_level_coverage_uses_the_union_footprint():
    shape = (15, 64, 64)
    stack, disc = _textured_blob(shape, (32, 32), 12, focus_plane=14)
    labels = _labels_from(disc, range(9, 13), shape)
    result = stack_z_coverage(stack, labels)
    assert result["stack_focus_peak_plane"] == 14
    assert result["stack_focus_peak_at_edge"] and result["stack_z_coverage_suspect"]


def test_shape_mismatch_is_rejected():
    with pytest.raises(ValueError, match="share one 3D ZYX grid"):
        object_z_coverage(np.zeros((5, 8, 8)), np.zeros((5, 8, 9), np.int32))
