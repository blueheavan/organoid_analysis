"""Round-2 estimator support: orbit-symmetry reduction of the minimax LP.

Round 1 capped the stencil radius at m <= 4 because the 577-variable LP at
m = 5 did not converge in practical time.  That cap forced the anisotropy limit
to 3.0 and left the confirmation grid thin (62 of 78 cases out of domain).

The reduction rests on an invariance argument, not an approximation:

  Let G be the group of signed axis permutations g that preserve the physical
  metric diag(spacing).  G maps the primitive direction set onto itself and
  maps unit normals to unit normals, so the objective

      F(w) = sup_n | sum_k w_k |n . u_k| - 1 |

  satisfies F(g.w) = F(w) for all g in G.  F is convex, so for any optimal w
  the G-average (1/|G|) sum_g g.w is also optimal.  Hence a G-symmetric
  optimum always exists, and restricting the search to weights that are
  constant on G-orbits loses nothing.  The reduced problem has one variable per
  orbit instead of one per direction.

Round 1's solver is left untouched in estimators_v2.py so its evidence stays
reproducible.  This module reuses estimators_v2.primitive_directions and
estimators_v2.fibonacci_sphere verbatim, so the direction set and the
verification normal set are identical between rounds.

Round 2 also changes the curvature scale in the a-priori budget from the
volume-equivalent radius R_eff to the inscribed radius r_in; see
out/ROUND2_PROTOCOL.md section 2.  Nothing else in the derivation changes.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from estimators_v2 import fibonacci_sphere, primitive_directions

WEIGHT_CACHE_V3 = Path(__file__).resolve().parent / "out" / "crofton_weights_v3"


# --------------------------------------------------------------------- group
def metric_symmetry_group(spacing) -> list[np.ndarray]:
    """Signed axis permutations preserving the metric diag(spacing).

    A signed permutation g with (g k)_i = s_i k_{p(i)} preserves
    ||spacing * k|| iff spacing[p(i)] == spacing[i] for every i.
    """
    sp = np.asarray(spacing, float)
    group = []
    for p in itertools.permutations(range(3)):
        if not np.allclose(sp[list(p)], sp):
            continue
        P = np.zeros((3, 3), int)
        for i, j in enumerate(p):
            P[i, j] = 1
        for s in itertools.product((1, -1), repeat=3):
            group.append(np.diag(s).astype(int) @ P)
    return group


def _canon(k: np.ndarray) -> tuple:
    """Representative of the +/- pair {k, -k}."""
    a, b = tuple(int(x) for x in k), tuple(int(-x) for x in k)
    return max(a, b)


def direction_orbits(spacing, m: int):
    """Partition primitive_directions(m) into G-orbits.

    Returns (K, labels, n_orbits) with ``labels[i]`` the orbit index of
    direction ``K[i]``.  Order of K is exactly primitive_directions(m), so
    downstream transition counting is unchanged.
    """
    K = primitive_directions(m)
    index = {_canon(k): i for i, k in enumerate(K)}
    group = metric_symmetry_group(spacing)
    labels = -np.ones(len(K), int)
    nxt = 0
    for i, k in enumerate(K):
        if labels[i] >= 0:
            continue
        members = {index[_canon(g @ k)] for g in group}
        for j in members:
            assert labels[j] < 0 or labels[j] == nxt
            labels[j] = nxt
        nxt += 1
    assert (labels >= 0).all()
    return K, labels, nxt


# ----------------------------------------------------------------- solve
def _cache_key(spacing, m: int) -> str:
    sp = np.asarray(spacing, float)
    r = sp / sp.max()                      # the LP depends only on the ratios
    return "sym_" + "_".join(f"{x:.10g}" for x in r) + f"_m{m}"


def crofton_weights_sym(spacing: tuple, m: int, *, time_budget_s: float = 600.0,
                        verbose: bool = False) -> tuple:
    """Symmetry-reduced minimax weights.

    Same contract as estimators_v2.crofton_weights(scheme='lp'): returns
    (dirs, w_times_area, diagnostics) with dirs in primitive_directions order.
    Raises TimeoutError if the cutting-plane loop exceeds ``time_budget_s``;
    the caller records the exclusion rather than silently lowering m.
    """
    import time

    WEIGHT_CACHE_V3.mkdir(parents=True, exist_ok=True)
    sp = np.asarray(spacing, float)
    path = WEIGHT_CACHE_V3 / (_cache_key(spacing, m) + ".json")
    if path.exists():
        d = json.loads(path.read_text())
        K = primitive_directions(m)
        w = np.asarray(d["w"], float)
        vk = K * sp
        Lk = np.linalg.norm(vk, axis=1)
        a_k = float(np.prod(sp)) / Lk
        diag = dict(d["diag"])
        diag["grazing_constant_D_um2"] = float(np.sum(w * Lk ** 2) / 24.0)
        diag["cached"] = True
        return tuple(map(tuple, K)), tuple(w * a_k), diag

    t0 = time.time()
    K, labels, n_orb = direction_orbits(spacing, m)
    vk = K * sp
    Lk = np.linalg.norm(vk, axis=1)
    uk = vk / Lk[:, None]
    a_k = float(np.prod(sp)) / Lk

    Nd = fibonacci_sphere(20011)
    Adense_full = np.abs(Nd @ uk.T)
    # Collapse columns onto orbits: the orbit weight multiplies the SUM of
    # |n.u_k| over the orbit's members.
    Adense = np.zeros((len(Nd), n_orb))
    np.add.at(Adense.T, labels, Adense_full.T)

    idx = np.arange(0, len(Nd), max(1, len(Nd) // max(600, 4 * n_orb)))
    c = np.zeros(n_orb + 1)
    c[-1] = 1.0
    t_star = None
    w_orb = None
    for it in range(120):
        if time.time() - t0 > time_budget_s:
            raise TimeoutError(f"spacing={tuple(sp)} m={m}: exceeded "
                               f"{time_budget_s}s at iteration {it}")
        A = Adense[idx]
        one = np.ones((len(A), 1))
        res = linprog(c, A_ub=np.vstack([np.hstack([A, -one]),
                                         np.hstack([-A, -one])]),
                      b_ub=np.concatenate([np.ones(len(A)), -np.ones(len(A))]),
                      bounds=[(0, None)] * (n_orb + 1), method="highs")
        if not res.success:
            raise RuntimeError(f"LP failed spacing={tuple(sp)} m={m}: {res.message}")
        w_orb, t_star = res.x[:n_orb], float(res.x[-1])
        dev = np.abs(Adense @ w_orb - 1.0)
        if verbose:
            print(f"    it{it:02d} orbits={n_orb} |cuts|={len(idx)} "
                  f"t*={t_star:.3e} sup={dev.max():.3e}", flush=True)
        if dev.max() <= t_star + 1e-5:
            break
        idx = np.union1d(idx, np.argsort(dev)[-max(200, n_orb // 2):])

    w = w_orb[labels]                      # expand back to per-direction
    R = np.abs(Nd @ uk.T) @ w
    D = float(np.sum(w * Lk ** 2) / 24.0)
    diag = {"n_directions": int(len(K)), "n_orbits": int(n_orb),
            "reduction_factor": float(len(K) / n_orb),
            "stencil_m": int(m), "scheme": "lp_sym",
            "group_order": int(len(metric_symmetry_group(spacing))),
            "n_constraint": int(len(idx)), "lp_optimum_t_star": t_star,
            "weight_sum": float(w.sum()), "grazing_constant_D_um2": D,
            "plane_response_max_abs_dev": float(np.abs(R - 1).max()),
            "plane_response_rms_dev": float(np.sqrt(np.mean((R - 1) ** 2))),
            "n_zero_weights": int((w <= 1e-12).sum()),
            "solve_seconds": float(time.time() - t0), "cached": False}
    store = dict(diag)
    store.pop("grazing_constant_D_um2")
    path.write_text(json.dumps({"w": [float(x) for x in w], "diag": store}))
    return tuple(map(tuple, K)), tuple(w * a_k), diag


def apriori_bound_v3(spacing: tuple, m: int, r_scale: float) -> float:
    """Round-2 budget: plane-response deviation + grazing deficit.

    ``r_scale`` is the INSCRIBED radius (round 2), not the volume-equivalent
    radius (round 1).  See ROUND2_PROTOCOL.md section 2.
    """
    d = crofton_weights_sym(tuple(float(s) for s in spacing), m)[2]
    return d["plane_response_max_abs_dev"] + d["grazing_constant_D_um2"] / r_scale ** 2


def choose_stencil_v3(spacing: tuple, r_scale: float, m_candidates) -> tuple:
    best = min(m_candidates, key=lambda m: apriori_bound_v3(spacing, m, r_scale))
    return best, apriori_bound_v3(spacing, best, r_scale)
