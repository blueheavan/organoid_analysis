import numpy as np
import pytest

from organoid_analysis.config import load_config
from organoid_analysis.microscopy_io.tiff_contract import Sample
from organoid_analysis.quantification.features import (
    LEGACY_SURFACE_AREA_METHOD,
    SURFACE_AREA_METHOD,
    geometry,
    legacy_surface_area,
    measure_instances,
)
from organoid_analysis.segmentation.watershed_instances import SegmentationResult

# REGRESSION tolerance, not the scientific acceptance criterion, and deliberately
# NOT tightened to the section 9 <5% bound: it must keep passing for shapes and
# resolutions outside the qualified domain, where no accuracy claim is made. The
# accuracy claim is tested in tests/quantification/test_surface_qualified_domain.py
# against the declared domain; passing 15% here is not evidence of accuracy.
SURFACE_REGRESSION_TOLERANCE=.15


def _reference_sphere():
    z,y,x=np.indices((25,49,49))
    return ((z-12)*2)**2+(y-24)**2+(x-24)**2<=18.**2


def test_anisotropic_sphere_agrees_with_analytic_geometry():
    spacing=(2.,1.,1.)
    radius=18.
    measured,_=geometry(_reference_sphere(),spacing)
    assert abs(measured['volume_um3']/(4*np.pi*radius**3/3)-1)<.025
    # A wide regression bound by design; see SURFACE_REGRESSION_TOLERANCE above.
    assert abs(measured['surface_area_um2']/(4*np.pi*radius**2)-1)<SURFACE_REGRESSION_TOLERANCE
    assert .85<measured['sphericity']<1.04


def test_surface_estimator_version_lock():
    """Change detector for the versioned production surface estimator.

    The expected value is the area this exact phantom takes under
    ``crofton_minimax_sym_v3``, recorded when that estimator was adopted
    (docs/evidence/2026-09-12-surface-crofton-v3). It is NOT an accuracy
    oracle. If this fails, the estimator changed: bump
    SURFACE_AREA_METHOD['method_version'] and repeat the surface V&V rather
    than editing this number.
    """
    measured,_=geometry(_reference_sphere(),(2.,1.,1.))
    assert SURFACE_AREA_METHOD['method_version']=='crofton_minimax_sym_v3'
    assert measured['surface_area_um2']==pytest.approx(4048.66364,rel=1e-6)


def test_superseded_estimator_stays_reproducible_and_still_fails():
    """The legacy estimator is preserved under its own identity, with its defect.

    Historical ``surface_area_um2`` values must stay reproducible, so the
    superseded estimator keeps its recorded value for this phantom -- and keeps
    failing the SCIENTIFIC_SPEC section 9 <5% criterion that the adopted
    estimator meets on the same mask.
    """
    mask=_reference_sphere()
    analytic=4*np.pi*18.**2
    assert LEGACY_SURFACE_AREA_METHOD['method_version']=='marching_cubes_binary_lewiner_v1'
    assert legacy_surface_area(mask,(2.,1.,1.))==pytest.approx(4549.28466796875,rel=1e-6)
    assert abs(legacy_surface_area(mask,(2.,1.,1.))/analytic-1)>.05
    assert abs(geometry(mask,(2.,1.,1.))[0]['surface_area_um2']/analytic-1)<.05


def test_domain_variables_are_exported_and_gate_is_not_silent():
    """A measurement outside the qualified domain is returned, flagged, not hidden.

    This phantom's inscribed radius is 9.01 voxels of the coarsest axis, just
    below the declared rho_in >= 10, so it is out of the qualified domain --
    which the estimator reports per object rather than refusing or silently
    returning an unqualified number.
    """
    measured,_=geometry(_reference_sphere(),(2.,1.,1.))
    assert measured['surface_rho_in']==pytest.approx(9.014,abs=5e-3)
    assert measured['surface_anisotropy']==pytest.approx(2.)
    assert measured['surface_stencil_radius'] in (1,2,3,4,5)
    assert measured['surface_in_qualified_domain'] is False
    assert 'rho_in' in measured['surface_domain_flags']
    assert measured['surface_area_um2']>0


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


def test_geometry_records_closed_and_open_cavity_topology():
    solid = np.ones((9, 9, 9), dtype=bool)
    closed = solid.copy()
    closed[2:7, 2:7, 2:7] = False
    opened = closed.copy()
    opened[:3, 4, 4] = False

    closed_geometry, _ = geometry(closed, (1., 1., 1.), fill_holes=True)
    open_geometry, _ = geometry(opened, (1., 1., 1.), fill_holes=True)
    assert closed_geometry["filled_void_voxels"] == 125
    assert closed_geometry["filled_void_components"] == 1
    assert closed_geometry["open_cavity_suspected"] is False
    assert open_geometry["filled_void_voxels"] == 0
    assert open_geometry["filled_void_components"] == 0
    assert open_geometry["open_cavity_suspected"] is True


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
