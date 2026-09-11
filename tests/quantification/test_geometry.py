import numpy as np
import pytest

from organoid_analysis.config import load_config
from organoid_analysis.microscopy_io.tiff_contract import Sample
from organoid_analysis.quantification.features import (
    SURFACE_AREA_METHOD,
    geometry,
    measure_instances,
)
from organoid_analysis.segmentation.watershed_instances import SegmentationResult

# REGRESSION tolerance, not the scientific acceptance criterion. SCIENTIFIC_SPEC
# section 9 requires <5% area error for analytical shapes, and this production
# estimator FAILS it (docs/evidence/2026-09-11-measurement-vv). The 15% bound
# only detects an unintended change in the estimator's known terracing bias;
# passing it is not evidence of surface accuracy.
SURFACE_REGRESSION_TOLERANCE=.15


def _reference_sphere():
    z,y,x=np.indices((25,49,49))
    return ((z-12)*2)**2+(y-24)**2+(x-24)**2<=18.**2


def test_anisotropic_sphere_agrees_with_analytic_geometry():
    spacing=(2.,1.,1.)
    radius=18.
    measured,_=geometry(_reference_sphere(),spacing)
    assert abs(measured['volume_um3']/(4*np.pi*radius**3/3)-1)<.025
    # Unsmooth voxel marching cubes has a discretization bias; do not assert perfect spheres.
    assert abs(measured['surface_area_um2']/(4*np.pi*radius**2)-1)<SURFACE_REGRESSION_TOLERANCE
    assert .85<measured['sphericity']<1.04


def test_surface_estimator_version_lock():
    """Change detector for the versioned production surface estimator.

    The expected value is the area recorded for this exact phantom in
    docs/evidence/2026-09-11/analytical-geometry.json and reproduced in
    docs/evidence/2026-09-11-measurement-vv. It is NOT an accuracy oracle: the
    analytical area is 4*pi*18**2 = 4071.5 um2 (+11.7% bias). If this fails, the
    estimator changed: bump SURFACE_AREA_METHOD['method_version'] and repeat the
    surface V&V rather than editing this number.
    """
    measured,_=geometry(_reference_sphere(),(2.,1.,1.))
    assert SURFACE_AREA_METHOD['method_version']=='marching_cubes_binary_lewiner_v1'
    assert measured['surface_area_um2']==pytest.approx(4549.28466796875,rel=1e-6)


def test_units_scale_volume_area_and_not_sphericity():
    mask=np.ones((9,13,17),bool)
    first,_=geometry(mask,(2.,.7,.8))
    second,_=geometry(mask,(4.,1.4,1.6))
    assert first['volume_um3']==pytest.approx(mask.sum()*2*.7*.8)
    assert second['volume_um3']/first['volume_um3']==pytest.approx(8)
    assert second['surface_area_um2']/first['surface_area_um2']==pytest.approx(4,rel=1e-6)
    assert second['sphericity']==pytest.approx(first['sphericity'],rel=1e-6)


def test_hollow_structure_reports_outer_envelope_and_void():
    mask=np.ones((15,15,15),bool)
    mask[3:12,3:12,3:12]=False
    measured,_=geometry(mask,(1.,1.,1.))
    assert measured['volume_um3']==3375
    assert measured['segmented_volume_um3']==3375-729
    assert measured['enclosed_void_fraction']==pytest.approx(729/3375)


def test_translation_does_not_change_shape_or_surface():
    mask=np.ones((9,11,13),bool)
    a,(va,fa)=geometry(mask,(3.,1.,2.))
    b,(vb,fb)=geometry(mask,(3.,1.,2.),origin_zyx=(20,30,40))
    assert a['volume_um3']==b['volume_um3']
    assert a['surface_area_um2']==pytest.approx(b['surface_area_um2'],rel=1e-6)
    assert b['centroid_z_um']-a['centroid_z_um']==60
    np.testing.assert_allclose(vb-va,np.broadcast_to([60,30,80],va.shape))


def test_truncated_organoid_is_preserved_but_excluded():
    labels=np.zeros((20,30,30),np.uint32)
    labels[0:10,5:20,5:20]=1
    sample=Sample(labels.astype(np.uint16),None,None,None,None,(2.,1.,1.),{})
    seg=SegmentationResult(labels,None,.1,0,[],{1:1})
    cfg=load_config()
    rows,_=measure_instances(sample,seg,cfg)
    assert len(rows)==1
    assert rows[0]['touches_border']
    assert not rows[0]['morphology_eligible']
    assert 'border_truncated' in rows[0]['morphology_flags']
