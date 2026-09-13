"""Numerical and reporting contracts independent of biological object identity."""
from __future__ import annotations

import numpy as np
import pytest

from organoid_analysis.quantification.cellular_measurements import cell_geometry
from organoid_analysis.quantification.features import SURFACE_METADATA_COLUMNS, geometry
from organoid_analysis.quantification.mask_features import extract_mask_features
from organoid_analysis.quantification.measurement_policy import measurement_policy
from organoid_analysis.quantification.multilevel_relationships.morphology import measure_instances


@pytest.mark.parametrize("shape", [(1, 1, 1), (3, 3, 3), (12, 12, 12)])
@pytest.mark.parametrize("spacing", [(1., 1., 1.), (2., 1., 1.)])
def test_cropping_preserves_geometry_and_inscribed_radius(shape, spacing):
    mask = np.ones(shape, bool)
    cropped, _ = geometry(mask, spacing, origin_zyx=(2, 2, 2), fill_holes=False)
    padded, _ = geometry(np.pad(mask, 2), spacing, fill_holes=False)
    # EDT is distance between foreground/background voxel centres, not an
    # inferred subvoxel boundary radius. A voxel box supplies a direct oracle.
    radius = min(((n + 1) // 2) * s for n, s in zip(shape, spacing))
    assert cropped["surface_rho_in"] == radius / max(spacing)
    for name in cropped:
        if isinstance(cropped[name], float):
            assert cropped[name] == pytest.approx(padded[name], rel=1e-12), name
        else:
            assert cropped[name] == padded[name], name


@pytest.mark.parametrize("object_type", ["organoid", "cell", "nucleus"])
def test_object_type_never_qualifies_an_underresolved_mask(object_type):
    labels = np.pad(np.ones((12, 12, 12), np.uint32), 1)
    row = measure_instances(labels, (1., 1., 1.), object_type=object_type).iloc[0]
    # Before exterior EDT padding, this tight crop incorrectly cleared rho>=10.
    assert row.surface_rho_in == 6
    assert not row.surface_in_qualified_domain
    assert not row.sphericity_in_qualified_domain
    assert not row.surface_to_volume_in_qualified_domain
    assert "rho_in" in row.surface_domain_flags
    assert row.volume_um3 == 12**3


@pytest.mark.parametrize("fill_holes", [False, True])
def test_same_support_has_same_geometry_across_workflows(fill_holes):
    labels = np.zeros((11, 11, 11), np.uint32)
    labels[1:10, 1:10, 1:10] = 71
    labels[2:5, 2:4, 2:6] = 0
    spacing = (2., 1., 1.)
    direct, _ = geometry(labels == 71, spacing, fill_holes=fill_holes)
    web = extract_mask_features(labels, spacing[::-1], fill_holes=fill_holes).loc[71]
    cell = cell_geometry(labels, None, spacing, fill_holes=fill_holes).iloc[0]
    if fill_holes:
        from organoid_analysis.quantification.features import outer_envelope

        support = outer_envelope(labels == 71).astype(np.uint32) * 71
    else:
        support = labels
    for level in ("organoid", "cell", "nucleus"):
        multilevel = measure_instances(support, spacing, object_type=level).iloc[0]
        for name in ("volume_um3", "equivalent_diameter_um", "surface_area_um2", "sphericity",
                     "axis_ratio_minor_to_major", "elongation", "prolate_ratio", "oblate_ratio",
                     "surface_to_volume_ratio_um_inv", "centroid_z_um", *SURFACE_METADATA_COLUMNS):
            assert web[name] == direct[name], name
            assert cell[name] == direct[name], name
            assert multilevel[name] == direct[name], name
        for field, engine, legacy in (
            ("major_axis_um", "principal_axis_major_um", "major_axis_um"),
            ("intermediate_axis_um", "principal_axis_intermediate_um", "minor_axis_um"),
            ("minor_axis_um", "principal_axis_minor_um", "least_axis_um"),
        ):
            assert multilevel[field] == cell[engine] == direct[engine] == web[legacy]


def test_numerical_gate_does_not_assert_surface_scope_or_volume_accuracy():
    # This creased cube passes the resolution/anisotropy gate. Smoothness
    # cannot be certified from it, and there is no segmentation reference.
    g, _ = geometry(np.ones((30, 30, 30), bool), (1., 1., 1.))
    assert g["surface_in_qualified_domain"]
    assert g["sphericity_in_qualified_domain"] == g["surface_in_qualified_domain"]
    assert g["surface_qualification_scope"] == "numerical_resolution_and_anisotropy_only"
    assert g["surface_weights_origin"] == "packaged"
    assert g["surface_weights_evidence_bearing"] is True
    assert g["surface_method_version"] == "crofton_minimax_sym_v3"
    assert g["surface_implementation_version"] == "padded_edt_v1"
    for level in ("organoid", "cell", "nucleus"):
        policy = measurement_policy(level, "raw_label")
        for name in ("surface_scope_conditions_status", "segmentation_validation_status",
                     "biological_measurement_validity_status", "volume_accuracy_status"):
            assert policy[name] == "NOT ASSESSED"


def test_small_nucleus_primary_metrics_have_independent_voxel_oracles():
    labels = np.zeros((5, 5, 5), np.uint32)
    labels[2, 2, 2] = 2**24 + 7
    row = extract_mask_features(labels, (1., 1., 2.)).iloc[0]
    assert row.voxel_count == 1
    assert row.volume_um3 == 2.
    assert row.equivalent_diameter_um == pytest.approx((12 / np.pi)**(1/3))
    # One cuboid voxel has covariance diag(s²/12), hence diameters 2 sqrt(5λ).
    np.testing.assert_allclose(
        row[["major_axis_um", "minor_axis_um", "least_axis_um"]].to_numpy(float),
        np.sqrt(5 / 3) * np.array([2., 1., 1.]), rtol=1e-12,
    )
    assert row.elongation == .5
    assert not row.surface_in_qualified_domain
    assert "surface_outside_qualified_domain" in row.qc_flags
