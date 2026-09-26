"""The volume-specific qualified domain: variable, gate and export."""
import math

import numpy as np
import pytest

from organoid_analysis.quantification import surface_crofton
from organoid_analysis.quantification.features import (
    SURFACE_METADATA_COLUMNS,
    VOLUME_METHOD,
    geometry,
)


def _ball(radius: float, spacing: tuple[float, float, float]) -> np.ndarray:
    sp = np.asarray(spacing)
    half = np.ceil(radius / sp).astype(int) + 2
    axes = [(np.arange(-h, h + 1) + 0.37) * s for h, s in zip(half, sp)]
    z, y, x = np.meshgrid(*axes, indexing="ij")
    return z * z + y * y + x * x <= radius * radius


def test_rho_vol_is_inscribed_radius_in_geometric_mean_voxel_edges():
    assert surface_crofton.volume_resolution(8.0, (4.0, 1.0, 1.0)) == pytest.approx(8.0 / 4 ** (1 / 3))
    assert surface_crofton.volume_resolution(8.0, (1.0, 1.0, 1.0)) == pytest.approx(8.0)


def test_volume_domain_requires_the_surface_gate_and_volume_resolution():
    ok = {"rho_in": 12.0, "anisotropy": 2.0, "rho_vol": surface_crofton.DOMAIN_RHO_VOL_MIN}
    assert surface_crofton.in_volume_domain(ok) == (True, ())
    low = {**ok, "rho_vol": surface_crofton.DOMAIN_RHO_VOL_MIN - 1e-9}
    assert surface_crofton.in_volume_domain(low)[0] is False
    coarse = {**ok, "rho_in": 9.99}
    inside, flags = surface_crofton.in_volume_domain(coarse)
    assert not inside and any("rho_in" in flag for flag in flags)
    anisotropic = {**ok, "anisotropy": 4.5}
    assert surface_crofton.in_volume_domain(anisotropic)[0] is False


def test_volume_domain_is_never_wider_than_the_surface_domain():
    # rho_vol >= rho_in always (the geometric mean never exceeds the coarsest edge).
    for spacing in [(1.0, 1.0, 1.0), (2.0, 1.0, 1.0), (4.0, 0.5, 0.5)]:
        assert surface_crofton.volume_resolution(10.0, spacing) >= 10.0 / max(spacing) - 1e-12
    assert surface_crofton.DOMAIN_RHO_VOL_MIN >= surface_crofton.DOMAIN_RHO_IN_MIN


def test_geometry_exports_the_volume_domain_with_every_surface_domain_field():
    radius = surface_crofton.DOMAIN_RHO_VOL_MIN + 2.0
    measured, _ = geometry(_ball(radius, (1.0, 1.0, 1.0)), (1.0, 1.0, 1.0))
    for column in ("volume_rho_vol", "volume_in_qualified_domain", "volume_domain_flags", "volume_method"):
        assert column in SURFACE_METADATA_COLUMNS
        assert column in measured
    assert measured["volume_in_qualified_domain"] is True
    assert measured["volume_rho_vol"] == pytest.approx(measured["surface_rho_in"])
    assert measured["volume_method"] == VOLUME_METHOD["method_version"]
    assert abs(measured["volume_um3"] / (4 / 3 * math.pi * radius**3) - 1) < 0.01


def test_small_object_is_surface_qualified_but_not_volume_qualified():
    measured, _ = geometry(_ball(11.0, (1.0, 1.0, 1.0)), (1.0, 1.0, 1.0))
    assert measured["surface_in_qualified_domain"] is True
    assert measured["volume_in_qualified_domain"] is False
