"""A binary mask cannot always identify continuous-object volume to 1%."""

import numpy as np

from organoid_analysis.quantification.features import geometry


def test_two_in_domain_spheres_have_identical_masks_but_disjoint_one_percent_intervals():
    axis = np.arange(-12, 13)
    z, y, x = np.meshgrid(axis, axis, axis, indexing="ij")
    squared_distance = z * z + y * y + x * x

    # No integer lattice point has squared radius 111 or 112. Both smooth,
    # closed spheres therefore produce the same centre-sampled binary mask.
    small_radius = np.sqrt(110) + 1e-5
    large_radius = np.sqrt(113) - 1e-5
    small_mask = squared_distance < small_radius**2
    large_mask = squared_distance < large_radius**2
    assert np.array_equal(small_mask, large_mask)
    assert small_radius > 10

    small_truth = 4 * np.pi * small_radius**3 / 3
    large_truth = 4 * np.pi * large_radius**3 / 3
    # An estimate from identical mask and spacing must be identical for both
    # spheres, yet their required <1% acceptance intervals do not intersect.
    assert small_truth * 1.01 < large_truth * 0.99

    measured, _ = geometry(small_mask, (1.0, 1.0, 1.0))
    assert measured["surface_rho_in"] >= 10
    assert measured["surface_in_qualified_domain"] is True
    estimate = measured["volume_um3"]
    assert np.isfinite(estimate)
    assert max(abs(estimate / small_truth - 1), abs(estimate / large_truth - 1)) > 0.01
