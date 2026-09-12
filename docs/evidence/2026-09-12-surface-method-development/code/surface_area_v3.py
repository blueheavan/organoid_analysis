"""Round-2 candidate estimator: crofton_minimax_sym_v3.

Production-shaped entry point.  Differences from round 2's predecessor
``crofton_minimax_adaptive_v2`` (all declared in out/ROUND2_PROTOCOL.md before
any round-2 phantom was generated):

  * the domain and the stencil-selection budget use the INSCRIBED radius
    r_in = max(distance_transform_edt(mask, sampling=spacing)) rather than the
    volume-equivalent radius R_eff.  r_in is a local feature-size measure and
    is what actually controls both error terms; R_eff overstates resolution
    for elongated or creased objects, which is what caused the round-1
    confirmation failure;
  * weights come from the orbit-symmetry-reduced minimax LP, so stencil radii
    m = 5, 6 are available (round 1 was capped at 4 because its LP did not
    converge there).

Both quantities are computable from the mask and the voxel spacing alone: no
truth value, no shape class, and no fitted constant enters a measurement.

DOMAIN CONSTANTS ARE DELIBERATELY UNSET in this file until the freeze step.
``measure()`` never applies a domain; ``in_domain()`` raises until the
declaration is written in.  This makes it impossible to run a domain-gated
result before the domain has been declared and hashed.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

from estimators_v2 import transition_counts
from estimators_v3 import apriori_bound_v3, crofton_weights_sym

METHOD_NAME = "crofton_minimax_sym_v3"
SCHEME = "lp_sym"
M_CANDIDATES = (1, 2, 3, 4, 5)

# m = 6 is EXCLUDED.  Recorded before the confirmation grid was executed, per
# out/ROUND2_PROTOCOL.md section 3.  Evidence: the symmetry-reduced LP solved
# m = 6 within budget for 10 of the 11 spacing ratios tried (max 250 s), but
# spacing 1.2x1x1 did not converge in 900 s / 80 cutting-plane iterations
# (out/lp_symmetry_verification.csv, out/m6_exclusion.json).  Production sees
# arbitrary spacings, so a radius that can fail to solve for some of them
# cannot be offered for any of them: the failure would appear at measurement
# time on a spacing not tested here.  The cap is a property of the declared
# method, not a compute-budget convenience -- m = 5 solved for every ratio
# tried, in at most 74 s.
M_EXCLUDED = (6,)

# --- qualified applicability domain -----------------------------------------
# Written by step9_freeze_v3.py from the development evidence, then hashed.
DOMAIN_DECLARED = True
DOMAIN_RHO_IN_MIN = 10.0        # r_inscribed / max(spacing)
DOMAIN_RHO_IN_MIN_VOLUME = 10.0
DOMAIN_ANISO_MAX = 4.0          # max(spacing) / min(spacing)

# SCOPE OF INTENDED USE (not machine-checkable; see FREEZE_RECORD_V3.md S4).
# The domain covers SMOOTH closed surfaces -- surfaces without dihedral
# creases.  Development evidence (out/vv_dev2_all.csv, 1083 cases) shows the
# estimator carries a systematic NEGATIVE area bias at dihedral edges that
# decays far too slowly to be removed by resolution: an isotropic cube is
# still -3.5% at rho_in ~ 30.  Surfaces with creases are therefore NOT
# QUALIFIED at any resolution in this round.
#
# The discriminator is creases, not elongation or topology: the capsule
# (aspect 3:1, smooth) reaches 1.12% while the cylinder (aspect 3:1, two
# circular creases) reaches 4-8% on the same grids; the torus (genus 1,
# smooth) reaches 0.74%.
#
# This is a scope statement, not a gate.  measure() cannot detect creases --
# a mask-computable crease indicator was tested and REJECTED (the marching-
# cubes/Crofton disagreement overlaps completely between classes, r = -0.20;
# out/crease_indicator_rejected.md).  Responsibility for staying in scope
# rests with the caller.
DOMAIN_SCOPE = "smooth closed surfaces (no dihedral creases)"
DOMAIN_NOT_QUALIFIED = "surfaces with dihedral creases, at any resolution"


def inscribed_radius(mask: np.ndarray, spacing) -> float:
    """Largest inscribed sphere radius, from the mask alone."""
    return float(ndi.distance_transform_edt(mask, sampling=spacing).max())


def measure(mask: np.ndarray, spacing, *, force_m: int | None = None) -> dict:
    """Surface area, volume and the domain variables.  No domain gating here."""
    mask = np.asarray(mask, bool)
    sp = tuple(float(s) for s in spacing)
    n_vox = int(np.count_nonzero(mask))
    if n_vox == 0:
        return {"method": METHOD_NAME, "surface_area": 0.0, "volume": 0.0,
                "r_in": 0.0, "rho_in": 0.0, "r_eff": 0.0, "rho_eff": 0.0,
                "stencil_m": None, "apriori_bound": None}

    volume = n_vox * float(np.prod(sp))
    r_eff = (3.0 * volume / (4.0 * np.pi)) ** (1.0 / 3.0)
    r_in = inscribed_radius(mask, sp)

    m = int(force_m) if force_m is not None else min(
        M_CANDIDATES, key=lambda mm: apriori_bound_v3(sp, mm, r_in))
    bound = apriori_bound_v3(sp, m, r_in)

    dirs, wa, diag = crofton_weights_sym(sp, m)
    area = float(np.dot(transition_counts(mask, dirs), wa))

    return {"method": METHOD_NAME, "surface_area": area, "volume": volume,
            "r_in": r_in, "rho_in": r_in / max(sp),
            "r_eff": r_eff, "rho_eff": r_eff / max(sp),
            "anisotropy": max(sp) / min(sp),
            "stencil_m": m, "n_directions": diag["n_directions"],
            "n_orbits": diag["n_orbits"],
            "plane_response_max_abs_dev": diag["plane_response_max_abs_dev"],
            "grazing_term": diag["grazing_constant_D_um2"] / r_in ** 2,
            "apriori_bound": bound}


def in_domain(meas: dict) -> tuple[bool, list[str]]:
    """Apply the frozen applicability domain.  Raises until it is declared."""
    if not DOMAIN_DECLARED:
        raise RuntimeError(
            "applicability domain not declared; run step9_freeze_v3.py first. "
            "A domain-gated result may not be produced before the domain is "
            "written and hashed.")
    reasons = []
    if meas["rho_in"] < DOMAIN_RHO_IN_MIN:
        reasons.append(f"rho_in {meas['rho_in']:.2f} < {DOMAIN_RHO_IN_MIN}")
    if meas["anisotropy"] > DOMAIN_ANISO_MAX:
        reasons.append(f"anisotropy {meas['anisotropy']:.2f} > {DOMAIN_ANISO_MAX}")
    return (not reasons), reasons


def in_domain_volume(meas: dict) -> tuple[bool, list[str]]:
    if not DOMAIN_DECLARED:
        raise RuntimeError("applicability domain not declared")
    reasons = []
    if meas["rho_in"] < DOMAIN_RHO_IN_MIN_VOLUME:
        reasons.append(f"rho_in {meas['rho_in']:.2f} < {DOMAIN_RHO_IN_MIN_VOLUME}")
    if meas["anisotropy"] > DOMAIN_ANISO_MAX:
        reasons.append(f"anisotropy {meas['anisotropy']:.2f} > {DOMAIN_ANISO_MAX}")
    return (not reasons), reasons


def sphericity(volume: float, area: float) -> float:
    return float(np.pi ** (1 / 3) * (6 * volume) ** (2 / 3) / area)
