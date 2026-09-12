"""Candidate surface-area estimators, with derivations.

Every numerical constant in this module is either (a) a property of the voxel
lattice and the physical spacing, or (b) the solution of a convex program whose
objective is stated in closed form below. Nothing here is fitted to a phantom
measurement, and nothing here is tuned against an error target. That is a
requirement of the brief, not a stylistic preference: a coefficient fitted to
phantoms would make the phantom study a calibration exercise rather than a
validation.

--------------------------------------------------------------------------
Candidate A -- continuous-field marching cubes
--------------------------------------------------------------------------
Mechanism M2 showed that the binary estimator's asymptotic bias comes from
pinning every mesh vertex to the midpoint of a lattice edge, which is all the
information a sign test leaves.  If instead a continuous field ``phi`` is
available whose iso-level ``L`` is the object boundary, linear interpolation
along the lattice edge places the vertex at its true subvoxel crossing, the
reconstructed facet normals converge to the true normals, and the orientation
bias vanishes to first order.

Two conventions are implemented because they correspond to two different things
a real pipeline might be able to retain:

  A1 ``signed_distance``  phi = signed distance (inside negative), level 0.
  A2 ``occupancy``        phi = partial volume fraction in [0, 1], level 0.5.

These are *upper-bound / reference* candidates for this project: the production
Cellpose wrapper currently discards the continuous probability field and keeps
only a thresholded mask, so A1/A2 cannot be evaluated on production data
without a pipeline change. They are carried through the full V&V anyway, for
two reasons: they bound what is achievable, and if they pass they define the
pipeline change worth making.

--------------------------------------------------------------------------
Candidate B -- generalised Cauchy-Crofton over an extended lattice stencil
--------------------------------------------------------------------------
Derivation (the weights follow from this; they are not chosen).

For a rectifiable surface S and a fixed unit direction u, integral geometry
gives the projection identity

    int over lines L parallel to u of  #(L cap S)  d(cross-sectional area)
        =  int_S |n(x) . u| dA(x)                                        (1)

On the voxel lattice, let k = (kz, ky, kx) be an integer vector with
gcd(|kz|,|ky|,|kx|) = 1 (primitive), and let v_k = k * spacing be the
corresponding physical vector, u_k = v_k / |v_k|.  The discrete lines in
direction k pass through lattice points, one point every |v_k| of arc length.
The lattice has 1 point per voxel volume V, so the number of distinct discrete
lines per unit cross-sectional area is

    1 / a_k       with      a_k = V / |v_k|                              (2)

Counting N_k = the number of adjacent-along-k voxel pairs whose membership
differs gives the total intersection count over all those lines, so from (1):

    N_k * a_k  ~=  int_S |n . u_k| dA                                    (3)

Now take any non-negative weights w_k over a set of K primitive directions and
sum (3):

    Ahat = sum_k w_k N_k a_k  =  int_S [ sum_k w_k |n . u_k| ] dA
                              =  int_S R(n) dA                           (4)

so the estimator is exact for every surface, at every orientation, exactly to
the extent that the **plane response**

    R(n) := sum_k w_k |n . u_k|                                          (5)

equals 1 for all unit n.  R(n) - 1 *is* the orientation bias of a Crofton
estimator, in the same units as the g(n) - 1 measured in M2. This makes the
weight choice a stated approximation-theory problem rather than a tuning knob.

Two weightings are implemented:

  B-vor  Spherical-Voronoi quadrature.  Treat (5) as a quadrature rule for
         int_{S^2} |n . u| dsigma(u) = 2*pi.  Assigning each direction the
         solid angle Omega_k of its spherical Voronoi cell (both hemispheres)
         and setting w_k = Omega_k / (2*pi) makes (5) exact in the limit of a
         dense direction set.  This is the classical choice and, at K = 13, it
         is exactly the estimator the project has already tested and rejected.
         It is reproduced here as a *reference*, not offered as new.

  B-lp   Minimax weights.  R(n) is linear in w, so

             minimise   t
             over       w >= 0,  t >= 0
             subject to |R(n_j) - 1| <= t   for a dense set {n_j} on S^2

         is a linear program, solved once per (spacing, stencil) from lattice
         geometry alone.  Its optimum t* is a *computable, phantom-free upper
         bound on the orientation bias* of the resulting estimator -- known
         before any phantom is measured.  Non-negativity is imposed so the
         estimator stays monotone in the transition counts.

Extending the stencil from the 3x3x3 neighbourhood (13 primitive directions) to
5x5x5 or 7x7x7 adds directions and can only lower t*.  This matters most under
anisotropy: with 3:1:1 spacing the *physical* directions generated by a 3x3x3
stencil are clustered toward z (index (1,1,0) maps to physical (3,1,0)), which
is precisely why the 13-direction estimator degraded on anisotropic cases.  A
7x7x7 stencil reaches physical directions such as (3,3,0)/|.| = 45 degrees.

Residual error sources of Candidate B, stated in advance:
  * (3) counts *digital* transitions, so features thinner than |v_k| along k
    are missed; this degrades at small rho and bounds the small-object domain.
  * N_k is an integer, so there is quantisation noise ~ 1/sqrt(N_k).
  * The edge term c/rho identified in M1 is not addressed by any binary-mask
    estimator, B included.
"""
from __future__ import annotations

import itertools
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.spatial import SphericalVoronoi
from skimage import measure


# --------------------------------------------------------------------------
# Direction stencils
# --------------------------------------------------------------------------
def primitive_directions(m: int) -> np.ndarray:
    """Primitive integer vectors in [-m, m]^3, one per +/- pair, sorted.

    Primitivity (gcd == 1) is required: a non-primitive vector such as (2,0,0)
    describes the same geometric line family as (1,0,0) but samples only every
    other lattice point, so (2) would not hold for it.
    """
    out = []
    for k in itertools.product(range(-m, m + 1), repeat=3):
        if k == (0, 0, 0):
            continue
        if np.gcd.reduce(np.abs(k)) != 1:
            continue
        if tuple(-np.array(k)) in {tuple(x) for x in out}:
            continue
        out.append(np.array(k))
    return np.array(sorted(out, key=lambda v: (np.abs(v).sum(), tuple(v))))


def fibonacci_sphere(n: int) -> np.ndarray:
    """Deterministic near-uniform unit vectors (upper hemisphere is enough:
    R(n) is even in n, so only n and -n together matter)."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(phi), np.sin(phi) * np.sin(theta),
                            np.sin(phi) * np.cos(theta)])  # (z, y, x)


# --------------------------------------------------------------------------
# Crofton weights
# --------------------------------------------------------------------------
WEIGHT_CACHE = Path(__file__).resolve().parent / "out" / "crofton_weights"


def _cache_key(spacing, m, scheme, n_constraint) -> str:
    s = "_".join(f"{float(x):g}" for x in spacing)
    return f"{scheme}_sp{s}_m{m}_nc{n_constraint}"


@lru_cache(maxsize=None)
def crofton_weights(spacing: tuple, m: int, scheme: str,
                    n_constraint: int = 0) -> tuple:
    """Disk-memoised wrapper -- the LP is deterministic, so a solved weight
    vector is reusable across processes and across V&V stages.  Caching only
    avoids re-solving; it can never change a result."""
    WEIGHT_CACHE.mkdir(parents=True, exist_ok=True)
    path = WEIGHT_CACHE / (_cache_key(spacing, m, scheme, n_constraint) + ".json")
    if path.exists():
        d = json.loads(path.read_text())
        return (tuple(map(tuple, d["dirs"])), tuple(d["wa"]), d["diag"])
    dirs, wa, diag = _crofton_weights_solve(spacing, m, scheme, n_constraint)
    path.write_text(json.dumps({"dirs": [list(map(int, k)) for k in dirs],
                                "wa": list(map(float, wa)), "diag": diag}))
    return dirs, wa, diag


def _crofton_weights_solve(spacing: tuple, m: int, scheme: str,
                           n_constraint: int = 0) -> tuple:
    """Return (dirs, w, diagnostics) for a stencil radius ``m``.

    ``scheme`` is 'voronoi' or 'lp'.  Deterministic: depends only on
    (spacing, m, scheme, n_constraint).

    The LP constraint set is sized at >= 20 normals per free weight.  With too
    few constraints the LP would satisfy them exactly while the true supremum
    of |R(n) - 1| over the sphere sat higher -- i.e. the reported bound would
    be optimistic.  ``t_star`` (the LP optimum, on the constraint set) and
    ``plane_response_max_abs_dev`` (on an independent 20011-normal set) are
    both reported so that gap is visible rather than assumed away.
    """
    sp = np.asarray(spacing, float)
    K = primitive_directions(m)
    if n_constraint == 0:
        n_constraint = max(2000, 20 * len(K))
    V = float(np.prod(sp))
    vk = K * sp                               # physical vectors
    Lk = np.linalg.norm(vk, axis=1)
    uk = vk / Lk[:, None]                     # physical unit directions
    a_k = V / Lk                              # cross-sectional area per line

    if scheme == "voronoi":
        pts = np.vstack([uk, -uk])
        sv = SphericalVoronoi(pts, radius=1.0, center=np.zeros(3))
        areas = sv.calculate_areas()
        omega = areas[:len(uk)] + areas[len(uk):]
        w = omega / (2 * np.pi)
    elif scheme == "lp":
        # Semi-infinite LP solved by cutting planes (exchange algorithm):
        # minimise sup_n |R(n)-1| over the whole sphere, not over a fixed
        # sample.  Start from a coarse normal set, solve, find the worst
        # violators on a dense set, add them, repeat.  On exit t_star (the LP
        # optimum) and the dense-set supremum agree to `tol`, so the reported
        # bound is verified rather than assumed.
        nk = len(uk)
        Ndense = fibonacci_sphere(20011)
        Adense = np.abs(Ndense @ uk.T)
        idx = np.arange(0, len(Ndense), max(1, len(Ndense) // max(600, 4 * nk)))
        c = np.zeros(nk + 1); c[-1] = 1.0
        tol, t_star, w = 1e-5, None, None
        for _ in range(12):
            A = Adense[idx]
            one = np.ones((len(A), 1))
            res = linprog(c, A_ub=np.vstack([np.hstack([A, -one]),
                                             np.hstack([-A, -one])]),
                          b_ub=np.concatenate([np.ones(len(A)), -np.ones(len(A))]),
                          bounds=[(0, None)] * (nk + 1), method="highs")
            if not res.success:
                raise RuntimeError(f"LP failed spacing={spacing} m={m}: {res.message}")
            w, t_star = res.x[:nk], float(res.x[-1])
            dev = np.abs(Adense @ w - 1.0)
            if dev.max() <= t_star + tol:
                break
            idx = np.union1d(idx, np.argsort(dev)[-max(200, nk // 2):])
        n_constraint = int(len(idx))
    else:
        raise ValueError(scheme)
    if scheme == "voronoi":
        t_star = None

    # Verify the plane response on an INDEPENDENT, denser normal set.
    Nv = fibonacci_sphere(20011)
    R = np.abs(Nv @ uk.T) @ w
    # A-priori grazing-chord deficit constant (see GRAZING_DERIVATION below).
    # relative deficit = D / R_eff^2, with D a pure lattice+weight constant.
    D = float(np.sum(w * Lk ** 2) / 24.0)

    diag = {"n_directions": int(len(K)), "stencil_m": int(m), "scheme": scheme,
            "n_constraint": int(n_constraint), "lp_optimum_t_star": t_star,
            "weight_sum": float(w.sum()), "grazing_constant_D_um2": D,
            "plane_response_min": float(R.min()), "plane_response_max": float(R.max()),
            "plane_response_max_abs_dev": float(np.abs(R - 1).max()),
            "plane_response_rms_dev": float(np.sqrt(np.mean((R - 1) ** 2))),
            "n_zero_weights": int((w <= 1e-12).sum())}
    return tuple(map(tuple, K)), tuple(w * a_k), diag


GRAZING_DERIVATION = """
Second error term of Candidate B, derived a priori (no phantom involved).

Equation (3) equates the digital transition count with the true intersection
count. They differ: lattice points along a discrete line of direction k are
spaced |v_k| apart, so a chord shorter than |v_k| may contain no lattice point
and contribute 0 transitions instead of 2. Near the silhouette every convex
body has such grazing chords, so the count is biased LOW, and the bias grows
with the stencil radius -- the opposite direction to the orientation bias.

For a chord of length c and a uniformly-placed lattice of period L = |v_k|, the
probability that at least one sample lands inside is min(1, c/L), so the
expected deficit is 2*max(0, 1 - c/L).  For a sphere of radius R, parametrise
lines by impact parameter b with chord c(b) = 2*sqrt(R^2 - b^2), and substitute
t = sqrt(R^2 - b^2) (so b db = -t dt, and c < L means t < L/2):

    missed  =  int_0^{L/2} 2 (1 - 2t/L) * 2 pi t dt
            =  4 pi [ t^2/2 - 2 t^3/(3L) ]_0^{L/2}
            =  4 pi ( L^2/8 - L^2/12 )
            =  pi L^2 / 6

against a true directional integral int_S |n.u_k| dA = 2 pi R^2, i.e. a
relative deficit of L_k^2 / (12 R^2) for direction k.  Weighting by w_k and
using sum_k w_k = 2 (which follows from R(n) == 1 integrated over S^2):

    relative deficit  =  sum_k w_k L_k^2 / (24 R^2)  =:  D / R^2         (6)

D depends only on the spacing, the stencil and the weights -- it is known
before any measurement.  Combining (5) and (6), the a-priori error budget of a
Crofton estimator on a smooth convex object of volume-equivalent radius R_eff:

    |relative area error|  <~  t*(spacing, m)  +  D(spacing, m) / R_eff^2   (7)

The two terms move oppositely in m, so (7) has an interior optimum.  This is
what the frozen stencil rule will minimise -- a computation on lattice
constants and the object's own voxel count, never on a phantom error.

Not covered by (7): the crease/edge term c/rho of M1 (non-smooth boundaries),
and integer quantisation of N_k.  Both are declared as residuals and are why
(7) is used with a safety margin rather than as an equality.
"""


def grazing_constant(spacing: tuple, m: int, scheme: str) -> float:
    """D in equation (6), physical length^2."""
    return crofton_weights(tuple(float(s) for s in spacing), m, scheme)[2][
        "grazing_constant_D_um2"]


def apriori_bound(spacing: tuple, m: int, scheme: str, r_eff: float) -> float:
    """Equation (7): a-priori relative-error budget, phantom-free."""
    d = crofton_weights(tuple(float(s) for s in spacing), m, scheme)[2]
    return d["plane_response_max_abs_dev"] + d["grazing_constant_D_um2"] / r_eff ** 2


def choose_stencil(spacing: tuple, r_eff: float, scheme: str,
                   m_candidates=(1, 2, 3, 4)) -> tuple:
    """Frozen selection rule: minimise the a-priori bound (7) over m.

    Inputs are the voxel spacing (known from the image metadata) and the
    object's volume-equivalent radius (known from its voxel count). No phantom,
    no truth value, and no tuned constant enters.
    """
    best = min(m_candidates, key=lambda m: apriori_bound(spacing, m, scheme, r_eff))
    return best, apriori_bound(spacing, best, scheme, r_eff)


def transition_counts(mask: np.ndarray, dirs) -> np.ndarray:
    """N_k for each direction: adjacent-along-k voxel pairs of differing label.

    Computed once per (mask, stencil) and reused by every weighting scheme and
    every smaller stencil, since primitive_directions(m) is nested in m.
    """
    dirs = np.asarray(dirs, int)
    pad = int(np.abs(dirs).max()) + 1
    mk = np.pad(mask.astype(bool), pad, constant_values=False)
    out = np.empty(len(dirs), float)
    for i, k in enumerate(dirs):
        a = mk[tuple(slice(max(ki, 0), mk.shape[d] + min(ki, 0))
                     for d, ki in enumerate(k))]
        b = mk[tuple(slice(max(-ki, 0), mk.shape[d] + min(-ki, 0))
                     for d, ki in enumerate(k))]
        out[i] = np.count_nonzero(a ^ b)
    return out


def crofton_area_from_counts(counts: np.ndarray, all_dirs, spacing: tuple,
                             m: int, scheme: str) -> float:
    """Equation (4), using counts precomputed on a superset stencil."""
    dirs, wa, _ = crofton_weights(tuple(float(s) for s in spacing), m, scheme)
    index = {tuple(k): i for i, k in enumerate(map(tuple, np.asarray(all_dirs, int)))}
    return float(sum(w * counts[index[k]] for k, w in zip(dirs, wa) if w > 0.0))


def crofton_area(mask: np.ndarray, spacing: tuple, m: int, scheme: str) -> float:
    """Generalised Cauchy-Crofton surface area, equation (4)."""
    dirs, wa, _ = crofton_weights(tuple(float(s) for s in spacing), m, scheme)
    return float(np.dot(transition_counts(mask, dirs), wa))


# --------------------------------------------------------------------------
# Continuous-field marching cubes
# --------------------------------------------------------------------------
def continuous_mc_area(field: np.ndarray, spacing: tuple, level: float,
                       outside_value: float) -> float:
    f = np.pad(np.asarray(field, float), 1, constant_values=outside_value)
    verts, faces, _, _ = measure.marching_cubes(f, level=level, spacing=tuple(spacing))
    return float(measure.mesh_surface_area(verts, faces))


def area_signed_distance(sdf: np.ndarray, spacing: tuple) -> float:
    """A1: inside-negative signed distance field, iso-level 0."""
    return continuous_mc_area(sdf, spacing, 0.0, float(np.max(sdf)) + float(max(spacing)))


def area_occupancy(occ: np.ndarray, spacing: tuple) -> float:
    """A2: partial-volume occupancy field in [0,1], iso-level 0.5."""
    return continuous_mc_area(occ, spacing, 0.5, 0.0)


# --------------------------------------------------------------------------
if __name__ == "__main__":
    rows = []
    for sp in [(1.0, 1.0, 1.0), (2.0, 0.7, 0.7), (3.0, 1.0, 1.0),
               (2.0, 1.0, 1.0), (1.5, 0.5, 0.5)]:
        for m in (1, 2, 3):
            for scheme in ("voronoi", "lp"):
                _, _, d = crofton_weights(sp, m, scheme)
                d["spacing"] = "x".join(str(s) for s in sp)
                rows.append(d)
    print(json.dumps(rows, indent=1))
