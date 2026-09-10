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
    assert model.calls[0][1]['diameter'] == 32 * factor
    assert model.calls[1][1]['diameter'] == 48 * factor
    assert all(call[1]['anisotropy'] == 4 * factor for call in model.calls)
    assert model.calls[0][0] == (4, int(16 * factor), int(16 * factor))
    assert nuclei.shape == cells.shape == image.shape


@pytest.mark.parametrize('image', [np.zeros((4, 4)), np.full((3, 4, 4), np.nan),
                                 np.full((3, 4, 4), np.inf), np.zeros((3, 4, 4), complex)])
def test_invalid_nuclei_only_input_rejected_before_model(image):
    model = RecordingModel()
    with pytest.raises(ValueError):
        segment_stacks(model, image, None, SegmentationConfig())
    assert not model.calls


@pytest.mark.parametrize('key,value', [('anisotropy', -1), ('xy_spacing_um', np.nan),
                                      ('nuclei_diameter', 0), ('model_type', 'unknown'),
                                      ('batch_size', 0), ('xy_downsample', .01)])
def test_invalid_configuration_rejected_before_model(key, value):
    model = RecordingModel()
    with pytest.raises(ValueError):
        segment_stacks(model, np.zeros((3, 8, 8)), None, replace(SegmentationConfig(), **{key: value}))
    assert not model.calls
