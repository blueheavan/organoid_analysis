"""The domain-restricted accuracy claim of the adopted surface estimator.

These tests assert the SCIENTIFIC_SPEC section 9 criterion (<5% absolute
relative area error) on analytical phantoms INSIDE the declared applicability
domain of ``crofton_minimax_sym_v3``, and assert that the exclusions are real
rather than decorative. They are a small in-domain subset, executed on every
CI run; they are not the V&V. The qualifying evidence is the untouched
confirmation set of 96 cases in docs/evidence/2026-09-12-surface-crofton-v3.

Placement matters: every phantom here sits at a sub-voxel offset, as the V&V
phantoms did. Exact lattice alignment is a measure-zero degenerate placement
that the V&V never sampled, and it is pinned separately.

Deliberately NOT asserted: that out-of-domain cases pass. They do not, that is
why the domain exists, and a test that quietly widened the domain by passing
everything would defeat the point.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from organoid_analysis.quantification import surface_crofton as sc
from organoid_analysis.quantification.features import SURFACE_AREA_METHOD, geometry

AREA_CRITERION = 0.05      # SCIENTIFIC_SPEC section 9
VOLUME_CRITERION = 0.01    # SCIENTIFIC_SPEC section 9
# The declared domain's smoothness scope, as shape classes of the frozen
# V&V grid: cylinders carry dihedral creases and are excluded by name.
SMOOTH_SHAPES = frozenset({"sphere", "ellipsoid", "capsule"})


# Sub-voxel offsets matching how the development and confirmation sets placed
# their phantoms. Exact lattice alignment is a degenerate placement and is
# tested separately in test_lattice_aligned_placement_is_a_known_exception.
GENERIC_OFFSET = (0.23, -0.41, 0.17)


def _grid(extent_um, spacing, offset=GENERIC_OFFSET):
    """Voxel-centre physical coordinates covering +/- extent_um, offset sub-voxel."""
    spacing = np.asarray(spacing, float)
    half = np.ceil(np.asarray(extent_um, float) / spacing).astype(int) + 3
    axes = [(np.arange(-h, h + 1) + o) * s
            for h, o, s in zip(half, np.asarray(offset, float), spacing)]
    return np.meshgrid(*axes, indexing="ij")


def sphere(radius, spacing, offset=GENERIC_OFFSET):
    z, y, x = _grid((radius,) * 3, spacing, offset)
    return z ** 2 + y ** 2 + x ** 2 <= radius ** 2


def ellipsoid(axes, spacing):
    z, y, x = _grid(axes, spacing)
    a, b, c = axes
    return (z / a) ** 2 + (y / b) ** 2 + (x / c) ** 2 <= 1.0


def capsule(radius, length, spacing):
    """A smooth elongated body: cylinder with hemispherical caps, no creases."""
    z, y, x = _grid((length / 2 + radius, radius, radius), spacing)
    axial = np.clip(z, -length / 2, length / 2)
    return (z - axial) ** 2 + y ** 2 + x ** 2 <= radius ** 2


def capsule_area(radius, length):
    return 4 * np.pi * radius ** 2 + 2 * np.pi * radius * length


# Smooth closed surfaces inside the declared domain: rho_in >= 10 and
# anisotropy <= 4. Each entry is (label, mask, spacing, area_true, volume_true).
def _in_domain_cases():
    cases = []
    for spacing in [(1.0, 1.0, 1.0), (2.0, 1.0, 1.0), (4.0, 1.0, 1.0)]:
        r = 12.0 * max(spacing)
        cases.append((f"sphere r={r:g} sp={spacing}", sphere(r, spacing), spacing,
                      4 * np.pi * r ** 2, 4 * np.pi * r ** 3 / 3))
    spacing = (2.0, 1.0, 1.0)
    axes = (48.0, 36.0, 24.0)  # min semi-axis 24 um = rho_in 12
    # Thomsen's approximation, p = 1.6075; error below 1e-3 for these axis ratios
    p = 1.6075
    area = 4 * np.pi * (((axes[0] * axes[1]) ** p + (axes[0] * axes[2]) ** p
                         + (axes[1] * axes[2]) ** p) / 3) ** (1 / p)
    cases.append((f"ellipsoid {axes} sp={spacing}", ellipsoid(axes, spacing), spacing,
                  area, 4 * np.pi * axes[0] * axes[1] * axes[2] / 3))
    r, length = 22.0, 66.0
    cases.append((f"capsule r={r:g} l={length:g} sp={spacing}", capsule(r, length, spacing), spacing,
                  capsule_area(r, length), 4 * np.pi * r ** 3 / 3 + np.pi * r ** 2 * length))
    return cases


IN_DOMAIN = _in_domain_cases()


@pytest.mark.parametrize("label,mask,spacing,area_true,volume_true",
                         IN_DOMAIN, ids=[case[0] for case in IN_DOMAIN])
def test_in_domain_phantom_meets_the_section_9_criteria(label, mask, spacing, area_true, volume_true):
    measured, _ = geometry(mask, spacing)
    assert measured["surface_in_qualified_domain"] is True, measured["surface_domain_flags"]
    assert abs(measured["surface_area_um2"] / area_true - 1) < AREA_CRITERION
    assert abs(measured["volume_um3"] / volume_true - 1) < VOLUME_CRITERION


def test_measured_error_stays_within_the_apriori_budget():
    """The a-priori budget predicts the error before the area is measured.

    It is computed from the weights and the inscribed radius, not fitted to the
    result. It is a validated predictive budget, NOT a proven upper bound: it
    held on all 1083 development and 96 confirmation cases, and
    test_lattice_aligned_placement_is_a_known_exception pins the degenerate
    placement that exceeds it.
    """
    for label, mask, spacing, area_true, _ in IN_DOMAIN:
        measurement = sc.measure(mask, spacing)
        realised = abs(measurement["surface_area"] / area_true - 1)
        assert realised <= measurement["apriori_bound"], (label, realised, measurement["apriori_bound"])


def test_lattice_aligned_placement_is_a_known_exception():
    """A sphere centred exactly on a voxel centre with integer radius.

    Found at adoption, not during the V&V: the development and confirmation
    sets placed every phantom at a random sub-voxel offset, so this degenerate
    placement -- where discretisation errors add coherently over the whole
    surface instead of averaging -- was never sampled. At rho_in = 12 it
    exceeds its own a-priori budget by a factor of ~1.3 and misses the section
    9 volume criterion, while the section 9 AREA criterion still holds with a
    factor of 2.8 in hand. Volume here is voxel counting, unchanged from the
    superseded pipeline, so this is a property of the project volume
    measurement rather than of the adopted surface estimator.

    Pinned as a test so the exception cannot be quietly forgotten. If a future
    change removes it, that is an improvement to be re-validated, not a
    regression -- update this test and say so in the evidence record.
    """
    radius, spacing = 12.0, (1.0, 1.0, 1.0)
    mask = sphere(radius, spacing, offset=(0.0, 0.0, 0.0))
    measured, _ = geometry(mask, spacing)
    measurement = sc.measure(mask, spacing)
    area_error = abs(measured["surface_area_um2"] / (4 * np.pi * radius ** 2) - 1)
    volume_error = abs(measured["volume_um3"] / (4 * np.pi * radius ** 3 / 3) - 1)
    assert measured["surface_in_qualified_domain"] is True
    assert area_error < AREA_CRITERION           # the area claim survives it
    assert area_error > measurement["apriori_bound"]   # the budget does not
    assert volume_error > VOLUME_CRITERION       # voxel-count volume does not


def test_creased_surface_is_excluded_and_would_fail():
    """The crease exclusion is load-bearing, not a formality.

    An isotropic cube at rho_in = 15 clears the machine-checkable gate, yet
    misses the <5% criterion -- which is why DOMAIN_SCOPE restricts intended
    use to smooth surfaces and why no resolution threshold can replace it.
    """
    side = 30.0
    cube = np.ones((int(side),) * 3, bool)
    measured, _ = geometry(cube, (1.0, 1.0, 1.0))
    assert measured["surface_in_qualified_domain"] is True   # resolution gate passes
    assert abs(measured["surface_area_um2"] / (6 * side ** 2) - 1) > AREA_CRITERION
    assert "crease" in SURFACE_AREA_METHOD["not_qualified"]


def test_under_resolved_object_is_flagged_out_of_domain():
    mask = sphere(4.0, (1.0, 1.0, 1.0))
    measured, _ = geometry(mask, (1.0, 1.0, 1.0))
    assert measured["surface_in_qualified_domain"] is False
    assert "rho_in" in measured["surface_domain_flags"]


def test_anisotropy_beyond_the_declared_limit_is_flagged():
    spacing = (5.0, 1.0, 1.0)
    measured, _ = geometry(sphere(60.0, spacing), spacing)
    assert measured["surface_anisotropy"] == pytest.approx(5.0)
    assert measured["surface_in_qualified_domain"] is False
    assert "anisotropy" in measured["surface_domain_flags"]


def test_declared_domain_constants_are_the_frozen_ones():
    """The domain was frozen before the confirmation run; it is not tunable."""
    assert sc.DOMAIN_RHO_IN_MIN == 10.0
    assert sc.DOMAIN_ANISO_MAX == 4.0
    assert sc.M_CANDIDATES == (1, 2, 3, 4, 5) and sc.M_EXCLUDED == (6,)
    assert SURFACE_AREA_METHOD["method_version"] == sc.METHOD_NAME == "crofton_minimax_sym_v3"


def test_packaged_weights_are_used_not_resolved():
    """The confirmed weight vectors ship with the package.

    The LP optimum is not unique, so a re-solve is a different conforming
    vector. Every spacing of the V&V grids must resolve to a packaged table, or
    the evidence would not re-execute bit-for-bit.
    """
    for spacing in [(1.0, 1.0, 1.0), (2.0, 1.0, 1.0), (3.0, 1.0, 1.0), (2.0, 0.7, 0.7),
                    (1.5, 1.0, 1.0), (2.5, 0.8, 0.8), (4.0, 1.0, 1.0)]:
        for m in sc.M_CANDIDATES:
            _, _, diagnostics = sc.crofton_weights_sym(spacing, m)
            assert diagnostics["weights_origin"] == "packaged", (spacing, m)


def test_plane_response_deviation_decreases_with_stencil_radius():
    """The mechanism behind the stencil rule: larger m buys plane response."""
    deviations = [sc.crofton_weights_sym((2.0, 1.0, 1.0), m)[2]["plane_response_max_abs_dev"]
                  for m in sc.M_CANDIDATES]
    assert deviations == sorted(deviations, reverse=True)


def test_recorded_vv_table_meets_the_criterion_inside_the_declared_domain():
    """The claim is checked against the record, not only against fresh phantoms.

    SG-1 in the validation record is unrestricted and still FAILS: the frozen
    grid contains creased cylinders and objects down to rho_in 1.5, which the
    declared domain excludes. This test reads the same recorded per-case table
    and asserts the criterion on the subset the accuracy claim actually covers.
    It must not be turned into a restatement of SG-1, and the subset must stay
    defined by the frozen domain constants rather than by a threshold chosen
    to make it pass.
    """
    root = Path(__file__).resolve().parents[2]
    pointer = json.loads((root / "docs/evidence/analytical_geometry_current_record.json").read_text())
    manifest = json.loads((root / pointer["manifest"]).read_text())
    assert manifest["estimator"]["surface_area_method"]["method_version"] == sc.METHOD_NAME
    table = root / pointer["manifest"].rsplit("/", 1)[0] / "surface_vv_dev.csv"
    rows = [r for r in csv.DictReader(table.open())
            if r["estimator"] == manifest["estimator"]["harness_label"]
            and r["shape"] in SMOOTH_SHAPES
            and float(r["rho"]) >= sc.DOMAIN_RHO_IN_MIN
            and float(r["anisotropy"]) <= sc.DOMAIN_ANISO_MAX]
    assert len(rows) >= 40, f"domain subset unexpectedly small: {len(rows)}"
    worst = max(abs(float(r["area_rel_err"])) for r in rows)
    assert worst < AREA_CRITERION, f"worst in-domain recorded area error {worst:.4%}"
