import numpy as np
import pytest
from scipy import ndimage as ndi
from analysis.features import marker_measurements
from analysis.config import load_config


def labels_and_bbox():
    labels=np.zeros((20,36,36),np.uint32)
    labels[6:14,12:24,12:24]=1
    return labels,ndi.find_objects(labels)[0]


def test_background_subtraction_retains_negative_signal_instead_of_clipping():
    labels,bbox=labels_and_bbox()
    channel=np.full(labels.shape,100,np.uint16)
    channel[labels==1]=90
    result=marker_measurements(labels,1,bbox,{'calcein':channel,'pi':channel},(2.,1.,1.),load_config()['quality'])
    assert result['calcein_mean_bg_corrected']==-10
    assert result['calcein_integrated_bg_corrected']==-10*(labels==1).sum()


def test_twelve_bit_saturation_in_uint16_is_caught_with_camera_limit():
    labels,bbox=labels_and_bbox()
    channel=np.full(labels.shape,100,np.uint16)
    channel[labels==1]=4095
    cfg=load_config()['quality']
    cfg['calcein_saturation_value']=4095
    cfg['pi_saturation_value']=4095
    result=marker_measurements(labels,1,bbox,{'calcein':channel,'pi':channel},(2.,1.,1.),cfg)
    assert result['calcein_saturated_fraction']==1
    assert not result['viability_measurement_eligible']
    assert 'calcein_saturated' in result['viability_measurement_flags']


def test_floating_point_marker_needs_explicit_saturation_limit():
    labels,bbox=labels_and_bbox()
    channel=np.full(labels.shape,100,dtype=np.float32)
    channel[labels==1]=500
    cfg=load_config()['quality']
    result=marker_measurements(labels,1,bbox,{'calcein':channel,'pi':channel},(2.,1.,1.),cfg)
    assert not result['viability_measurement_eligible']
    cfg['calcein_saturation_value']=1000
    cfg['pi_saturation_value']=1000
    result=marker_measurements(labels,1,bbox,{'calcein':channel,'pi':channel},(2.,1.,1.),cfg)
    assert result['viability_measurement_eligible']


def test_structure_only_measurement_skips_background_distance_transform(monkeypatch):
    labels,bbox=labels_and_bbox()
    monkeypatch.setattr(ndi,'distance_transform_edt',lambda *args,**kwargs: (_ for _ in ()).throw(AssertionError("unexpected EDT")))
    result=marker_measurements(labels,1,bbox,{'calcein':None,'pi':None},(2.,1.,1.),load_config()['quality'])
    assert not result['viability_measurement_eligible']
    assert result['background_voxels']==0
    assert 'missing_calcein_channel' in result['viability_measurement_flags']
