"""Frozen candidate estimator: crofton_minimax_adaptive_v2.

This is the production-shaped entry point that the V&V qualifies, and the one
the proposed patch installs alongside (never in place of) the legacy
``marching_cubes_binary_lewiner_v1``.

Everything below is FROZEN as of the freeze record. No parameter here was
fitted to a phantom:

  * weights          minimax (Chebyshev) solution of the plane-response
                     problem, from lattice geometry and the voxel spacing only
  * stencil radius m argmin over M_CANDIDATES of the a-priori error budget
                     t*(spacing, m) + D(spacing, m) / R_eff^2, where R_eff is
                     the object's own volume-equivalent radius
  * M_CANDIDATES     (1, 2, 3, 4).  The cap at 4 is a declared limitation, not
                     a result: the minimax weight solve at m >= 5 (K >= 577
                     directions) did not converge in practical time with the
                     dense cutting-plane formulation used here.  It is the
                     reason the qualified anisotropy range stops where it does
                     (see FREEZE_RECORD.md).  Extending it needs the
                     orbit-symmetry reduction noted there, not a new parameter.

The function additionally returns whether the measurement falls inside the
qualified applicability domain. It reports, never silently corrects: an
out-of-domain object still gets a number, flagged ``in_domain=False``.
"""
from __future__ import annotations

import numpy as np

from estimators_v2 import (apriori_bound, crofton_weights, transition_counts)

METHOD_NAME = "crofton_minimax_adaptive_v2"
SCHEME = "lp"
M_CANDIDATES = (1, 2, 3, 4)

# --- qualified applicability domain (frozen; basis in FREEZE_RECORD.md) -----
DOMAIN_RHO_EFF_MIN = 8.0      # R_eff / max(spacing)
DOMAIN_ANISO_MAX = 3.0        # max(spacing) / min(spacing)


def measure(mask: np.ndarray, spacing, *, force_m: int | None = None) -> dict:
    """Surface area of a binary label mask.

    Parameters
    ----------
    mask     boolean array, one connected object, True = inside
    spacing  physical voxel size per axis, same order as ``mask.shape``

    Returns a dict with the area, the selected stencil, the a-priori budget,
    and the applicability-domain verdict.
    """
    mask = np.asarray(mask, bool)
    sp = tuple(float(s) for s in spacing)
    voxel_volume = float(np.prod(sp))
    n_vox = int(np.count_nonzero(mask))
    volume = n_vox * voxel_volume
    r_eff = (3.0 * volume / (4.0 * np.pi)) ** (1.0 / 3.0) if n_vox else 0.0

    if n_vox == 0:
        return {"method": METHOD_NAME, "surface_area": 0.0, "volume": 0.0,
                "r_eff": 0.0, "rho_eff": 0.0, "stencil_m": None,
                "apriori_bound": None, "in_domain": False,
                "domain_reasons": ["empty mask"]}

    if force_m is not None:
        m = int(force_m)
    else:
        m = min(M_CANDIDATES, key=lambda mm: apriori_bound(sp, mm, SCHEME, r_eff))
    bound = apriori_bound(sp, m, SCHEME, r_eff)

    dirs, wa, diag = crofton_weights(sp, m, SCHEME)
    area = float(np.dot(transition_counts(mask, dirs), wa))

    rho_eff = r_eff / max(sp)
    aniso = max(sp) / min(sp)
    reasons = []
    if rho_eff < DOMAIN_RHO_EFF_MIN:
        reasons.append(f"rho_eff {rho_eff:.2f} < {DOMAIN_RHO_EFF_MIN}")
    if aniso > DOMAIN_ANISO_MAX:
        reasons.append(f"anisotropy {aniso:.2f} > {DOMAIN_ANISO_MAX}")

    return {"method": METHOD_NAME, "surface_area": area, "volume": volume,
            "r_eff": r_eff, "rho_eff": rho_eff, "anisotropy": aniso,
            "stencil_m": m, "n_directions": diag["n_directions"],
            "plane_response_max_abs_dev": diag["plane_response_max_abs_dev"],
            "grazing_term": diag["grazing_constant_D_um2"] / r_eff ** 2,
            "apriori_bound": bound,
            "in_domain": not reasons, "domain_reasons": reasons}


def sphericity(volume: float, area: float) -> float:
    return float(np.pi ** (1 / 3) * (6 * volume) ** (2 / 3) / area)
