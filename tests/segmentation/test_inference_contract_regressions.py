"""Adapter verification: fake models test scale and failure contracts, not accuracy."""
from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from organoid_analysis.segmentation.cellpose_inference import SegmentationConfig, segment_stacks


class RecordingModel:
    def __init__(self):
        self.calls = []

    def eval(self, image, **kwargs):
        self.calls.append((image.shape, kwargs))
        return (np.ones(image.shape[:3], np.uint32),)


@pytest.mark.parametrize('factor', [1., .5, .25])
def test_working_diameter_and_anisotropy_use_same_xy_scale(factor):
    image = np.ones((4, 16, 16), np.uint16)
    model = RecordingModel()
    config = SegmentationConfig(nuclei_diameter=32, cell_diameter=48, anisotropy=4, xy_downsample=factor)
    nuclei, cells = segment_stacks(model, image, image, config)
    # SciPy's center-grid span is N-1, so 16 -> 8 pixels scales by 7/15.
    effective_scale = (round(16 * factor) - 1) / 15
    assert model.calls[0][1]['diameter'] == 32 * effective_scale
    assert model.calls[1][1]['diameter'] == 48 * effective_scale
    assert all(call[1]['anisotropy'] == 4 * effective_scale for call in model.calls)
    assert model.calls[0][0] == (4, int(16 * factor), int(16 * factor))
    assert nuclei.shape == cells.shape == image.shape


@pytest.mark.parametrize('image', [np.zeros((4, 4)), np.full((3, 4, 4), np.nan),
                                 np.full((3, 4, 4), np.inf), np.zeros((3, 4, 4), complex)])
def test_invalid_nuclei_only_input_rejected_before_model(image):
    model = RecordingModel()
    with pytest.raises(ValueError):
        segment_stacks(model, image, None, SegmentationConfig())
    assert not model.calls


def test_unequal_rounded_xy_scales_reject_before_inference():
    model = RecordingModel()
    with pytest.raises(ValueError, match="full XY resolution"):
        segment_stacks(model, np.zeros((3, 11, 17)), None, SegmentationConfig(xy_downsample=.5))
    assert not model.calls


def test_interpolation_keeps_fractional_intensities_and_sparse_labels(monkeypatch):
    from organoid_analysis.segmentation import cellpose_inference

    monkeypatch.setattr(cellpose_inference, "_wrap_run_3d", lambda *args: lambda: None)
    ramp = np.broadcast_to(np.arange(17, dtype=np.uint16), (3, 17, 17)).copy()
    before = ramp.copy()

    class RampModel:
        def eval(self, image, **kwargs):
            np.testing.assert_allclose(image[0, 0], np.linspace(0, 16, 8), atol=1e-6)
            assert image.dtype == np.float32
            assert kwargs["diameter"] == 30 * 7 / 16
            result = np.zeros(image.shape, np.uint32)
            result[:, :, 4:] = 2**24 + 19
            return (result,)

    labels, cells = segment_stacks(RampModel(), ramp, None, SegmentationConfig(xy_downsample=.5))
    assert cells is None
    assert labels.shape == ramp.shape
    assert set(np.unique(labels)) == {0, 2**24 + 19}
    assert np.all(labels[:, :, 0] == 0) and np.all(labels[:, :, -1] == 2**24 + 19)
    np.testing.assert_array_equal(ramp, before)


@pytest.mark.parametrize('key,value', [('anisotropy', -1), ('xy_spacing_um', np.nan),
                                      ('nuclei_diameter', 0), ('model_type', 'unknown'),
                                      ('batch_size', 0), ('xy_downsample', .01)])
def test_invalid_configuration_rejected_before_model(key, value):
    model = RecordingModel()
    with pytest.raises(ValueError):
        segment_stacks(model, np.zeros((3, 8, 8)), None, replace(SegmentationConfig(), **{key: value}))
    assert not model.calls
