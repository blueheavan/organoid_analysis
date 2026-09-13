"""Qualified surface-area estimator ``crofton_minimax_sym_v3``.

Surface area from weighted counts of mask transitions along lattice
directions, i.e. a discrete Cauchy/Crofton formula: for a convex body the
number of intersections of a line family with the boundary is proportional to
the projected area, so a nonnegative combination of directional transition
counts estimates area. The direction weights are not chosen by a heuristic;
they solve the minimax linear program

    minimise   t
    over       w >= 0,  t
    subject to | sum_k w_k |n . u_k| - 1 | <= t   for all unit normals n

whose optimum is the weight vector whose plane response deviates least, in the
worst case over orientations, from the exact value 1. The LP optimum ``t*``
therefore bounds the estimator's error on locally planar surface elements
before any measurement is taken, which is what makes the a-priori budget in
``measure()`` a budget rather than a fitted constant.

Two terms enter that budget: the plane-response deviation ``t*`` above, and a
grazing deficit ``D / r_in**2`` for boundary elements whose normal is nearly
orthogonal to a lattice direction, where ``D = sum_k w_k L_k**2 / 24`` and
``r_in`` is the inscribed radius. The stencil radius ``m`` is selected per
measurement to minimise their sum; larger ``m`` lowers ``t*`` and raises ``D``.

The budget's STATUS: ``t*`` is a proven worst-case bound for a planar element,
but the grazing term is an asymptotic model of how curved and grazing elements
depart from planarity, so their sum is a *validated predictive budget, not a
proven upper bound on the realised error*. It held on all 1083 development and
all 96 confirmation cases. It is exceeded -- by a factor of about 1.3, while
the section 9 area criterion still holds -- by spheres centred exactly on a
voxel centre with integer radius, a measure-zero placement that neither set
sampled because both offset every phantom sub-voxel. See
``tests/quantification/test_surface_qualified_domain.py``, which pins that
exception. Treat ``apriori_bound`` as the expected error scale of a
measurement, not as a guarantee about it.

Solving the LP directly is impractical beyond ``m = 4`` (577 free variables at
``m = 5``). The search is restricted to weights that are constant on the
orbits of the group G of signed axis permutations preserving the physical
metric ``diag(spacing)``. That restriction is exact, not an approximation: G
maps the primitive direction set onto itself and unit normals to unit normals,
so the objective is G-invariant; being convex, it attains its optimum at the
G-average of any optimum, which is G-symmetric. At ``m = 5`` this leaves 98
variables instead of 577.

VALIDATION AND SCOPE
--------------------
Qualified on an untouched confirmation set of 96 in-domain cases: worst
absolute relative area error 0.951% against the <5% criterion of
docs/SCIENTIFIC_SPEC.md section 9. Evidence, freeze record and the
confirmation set are in docs/evidence/2026-09-12-surface-crofton-v3.

The qualified applicability domain has a machine-checkable part, applied by
``in_domain()``, and a scope of intended use that is NOT machine-checkable and
remains the caller's responsibility -- see ``DOMAIN_SCOPE`` below. Surfaces
with dihedral creases are NOT QUALIFIED at any resolution.

Any change to the weights, the stencil rule, the domain constants or the
transition counting must bump ``method_version`` in
``quantification.features.SURFACE_AREA_METHOD`` and repeat the V&V, so
historical values are never silently reinterpreted.
"""
from __future__ import annotations

import itertools
import json
import os
from functools import cache
from pathlib import Path
from typing import Any, cast

import numpy as np
from numpy.typing import ArrayLike
from scipy import ndimage as ndi  # type: ignore[import-untyped]
from scipy.optimize import linprog  # type: ignore[import-untyped]

METHOD_NAME = "crofton_minimax_sym_v3"
SCHEME = "lp_sym"

# Stencil radii offered to the selection rule.
M_CANDIDATES = (1, 2, 3, 4, 5)

# m = 6 is EXCLUDED, and the exclusion is a property of the declared method
# rather than a compute-budget convenience. The symmetry-reduced LP solved
# m = 6 within budget for 10 of the 11 spacing ratios tried (at most 250 s),
# but ratio 1.2:1:1 did not converge in 900 s / 80 cutting-plane iterations.
# Production sees arbitrary spacings, so a radius that can fail to solve for
# some of them cannot be offered for any of them: the failure would surface at
# measurement time on a spacing never tested. m = 5 solved for every ratio
# tried, in at most 74 s.
M_EXCLUDED = (6,)

# --- qualified applicability domain -----------------------------------------
# Declared and hashed from development evidence (1083 phantom cases) BEFORE the
# confirmation set was executed, and before this estimator was applied to any
# phantom of this repository's frozen grid. The constants were therefore not
# selected post hoc from the results they are reported against, which is the
# condition docs/SCIENTIFIC_CODE_AUDIT_2026-09-11.md F2 places on any change of
# support range.
DOMAIN_RHO_IN_MIN = 10.0        # inscribed radius / max(spacing)
DOMAIN_ANISO_MAX = 4.0          # max(spacing) / min(spacing)

# rho_in >= 10 is set by the VOLUME criterion (<1%), not the area criterion:
# area is comfortable from rho_in ~ 5, while voxel-count volume quantisation is
# the binding term on those phantoms. This historical selection rationale
# does not qualify volume accuracy on an arbitrary mask: the lattice-aligned
# counterexample can pass this gate while failing the volume criterion.
DOMAIN_RHO_IN_MIN_VOLUME = 10.0

# SCOPE OF INTENDED USE -- not machine-checkable.
#
# The domain covers SMOOTH closed surfaces: surfaces without dihedral creases.
# Development evidence shows a systematic NEGATIVE area bias at dihedral edges
# that decays far too slowly to be removed by resolution -- an isotropic cube is
# still -3.5% at rho_in ~ 30, so no resolution gate can recover it. The
# discriminator is creases, not elongation and not topology: a capsule (aspect
# 3:1, smooth) reaches 1.12% while a cylinder (aspect 3:1, two circular
# creases) reaches 4-8% on the same grids, and a torus (genus 1, smooth)
# reaches 0.74%.
#
# ``measure()`` cannot detect creases. A mask-computable crease indicator was
# tested and REJECTED before the confirmation run: the marching-cubes/Crofton
# disagreement, the only truth-free indicator available from a single mask,
# overlaps completely between classes (box 1.3-23.0%, sphere 8.8-17.1%) and
# correlates with the error at r = -0.20. Responsibility for staying in scope
# therefore rests with the caller, and ``in_domain()`` does not pretend to
# enforce it.
#
# Object identity (including an organoid envelope) does not establish smoothness
# or a biologically accurate segmentation boundary.
DOMAIN_SCOPE = "smooth closed surfaces (no dihedral creases)"
DOMAIN_NOT_QUALIFIED = "surfaces with dihedral creases, at any resolution"

# Confirmed ranges. Above rho_in 15.1 the domain rests on the monotone
# development trend, not on confirmation evidence; this is recorded so callers
# can see where the verified band ends.
DOMAIN_RHO_IN_CONFIRMED = (10.4, 15.1)
DOMAIN_ANISO_CONFIRMED = (1.2, 4.0)

# The FROZEN weight tables: the exact vectors solved during round-2 development
# and used to produce the confirmation evidence. They are shipped, not
# re-solved, because the LP optimum is NOT UNIQUE -- several weight vectors
# attain the same minimax deviation t*, and an independent solve of the same
# program lands on a different optimal vertex (weights differing by up to 78%
# elementwise, grazing constant D by up to 14%, with t* agreeing). The
# qualified estimator is therefore a specific frozen vector set, not "whatever
# the solver returns"; re-solving would silently change exported areas across
# solver versions without a method_version bump.
#
# These tables cover every spacing ratio of both phantom grids of the surface
# V&V, so the evidence re-executes exactly and without an LP solve.
#
# A spacing ratio absent from the tables is solved on demand and cached under
# ORGANOID_CROFTON_CACHE (default: the user cache directory). Such a solve is
# CONFORMING BUT NOT EVIDENCE-BEARING: it satisfies the same program and so
# carries the same a-priori bound, which ``measure()`` reports per object, but
# the specific vector was not exercised by the confirmation set.
# ``measure()`` returns ``weights_origin`` so this is visible in exports.
PACKAGED_WEIGHTS = Path(__file__).resolve().parent / "crofton_weights"
LP_TIME_BUDGET_S = 600.0

# This packaged table was added after the untouched confirmation set was
# frozen.  It is conforming and reproducible, but it must not be presented as
# confirmation evidence until a later V&V record exercises it.
_POST_CONFIRMATION_PACKAGED_RATIO = (1.0, 2.0 / 3.0, 2.0 / 3.0)
_POST_CONFIRMATION_PACKAGED_M = 5


def _weights_are_evidence_bearing(origin: str, ratios: tuple[float, float, float], m: int) -> bool:
    return not (
        origin == "packaged"
        and m == _POST_CONFIRMATION_PACKAGED_M
        and np.allclose(ratios, _POST_CONFIRMATION_PACKAGED_RATIO, rtol=0.0, atol=1e-10)
    )


# --------------------------------------------------------------- directions
@cache
def primitive_directions(m: int) -> tuple[tuple[int, int, int], ...]:
    """Primitive integer vectors in [-m, m]^3, one per +/- pair, sorted.

    Primitivity (gcd == 1) is required: a non-primitive vector such as
    (2, 0, 0) describes the same geometric line family as (1, 0, 0) but samples
    only every other lattice point, so the transition count along it does not
    carry the same normalisation.
    """
    out: list[tuple[int, int, int]] = []
    seen: set[tuple[int, int, int]] = set()
    for k in itertools.product(range(-m, m + 1), repeat=3):
        if k == (0, 0, 0) or int(np.gcd.reduce(np.abs(k))) != 1:
            continue
        if tuple(-np.asarray(k)) in seen:
            continue
        direction = cast(tuple[int, int, int], k)  # product(..., repeat=3)
        out.append(direction)
        seen.add(direction)
    return tuple(sorted(out, key=lambda v: (sum(abs(x) for x in v), v)))


def fibonacci_sphere(n: int) -> np.ndarray:
    """Deterministic near-uniform unit vectors in (z, y, x).

    The plane response is even in n, so one hemisphere suffices.
    """
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(phi), np.sin(phi) * np.sin(theta),
                            np.sin(phi) * np.cos(theta)])


def transition_counts(mask: np.ndarray, dirs: ArrayLike) -> np.ndarray:
    """N_k for each direction: adjacent-along-k voxel pairs of differing state."""
    directions = np.asarray(dirs, int)
    pad = int(np.abs(directions).max()) + 1
    padded = np.pad(np.asarray(mask, bool), pad, constant_values=False)
    out = np.empty(len(directions), float)
    for i, k in enumerate(directions):
        a = padded[tuple(slice(max(ki, 0), padded.shape[d] + min(ki, 0))
                         for d, ki in enumerate(k))]
        b = padded[tuple(slice(max(-ki, 0), padded.shape[d] + min(-ki, 0))
                         for d, ki in enumerate(k))]
        out[i] = np.count_nonzero(a ^ b)
    return out


# ------------------------------------------------------------------- group
def metric_symmetry_group(ratios: tuple[float, float, float]) -> list[np.ndarray]:
    """Signed axis permutations preserving the metric diag(ratios).

    A signed permutation g with (g k)_i = s_i k_{p(i)} preserves
    ||ratios * k|| iff ratios[p(i)] == ratios[i] for every i.
    """
    spacing = np.asarray(ratios, float)
    group: list[np.ndarray] = []
    for permutation in itertools.permutations(range(3)):
        if not np.allclose(spacing[list(permutation)], spacing):
            continue
        matrix = np.zeros((3, 3), int)
        for i, j in enumerate(permutation):
            matrix[i, j] = 1
        for signs in itertools.product((1, -1), repeat=3):
            group.append(np.diag(signs).astype(int) @ matrix)
    return group


def _canonical(k: np.ndarray) -> tuple[int, int, int]:
    """Representative of the +/- pair {k, -k}."""
    forward = tuple(int(x) for x in k)
    return cast(tuple[int, int, int], max(forward, tuple(-x for x in forward)))


def direction_orbits(ratios: tuple[float, float, float], m: int) -> tuple[np.ndarray, np.ndarray, int]:
    """Partition ``primitive_directions(m)`` into G-orbits.

    Returns ``(K, labels, n_orbits)`` where ``labels[i]`` is the orbit index of
    ``K[i]``. The order of K is exactly ``primitive_directions(m)``, so
    transition counting is unaffected.
    """
    directions = np.asarray(primitive_directions(m), int)
    index = {_canonical(k): i for i, k in enumerate(directions)}
    group = metric_symmetry_group(ratios)
    labels = -np.ones(len(directions), int)
    nxt = 0
    for i, k in enumerate(directions):
        if labels[i] >= 0:
            continue
        for j in {index[_canonical(g @ k)] for g in group}:
            if labels[j] >= 0 and labels[j] != nxt:
                raise AssertionError("orbit partition is inconsistent")
            labels[j] = nxt
        nxt += 1
    if (labels < 0).any():
        raise AssertionError("orbit partition left directions unassigned")
    return directions, labels, nxt


# ------------------------------------------------------------------- solve
def _ratios(spacing: ArrayLike) -> tuple[float, float, float]:
    """The LP depends only on the spacing ratios, so normalise by the maximum."""
    values = np.asarray(spacing, float)
    normalised = values / values.max()
    return tuple(round(float(x), 10) for x in normalised)  # type: ignore[return-value]


def _cache_key(ratios: tuple[float, float, float], m: int) -> str:
    return "sym_" + "_".join(f"{x:.10g}" for x in ratios) + f"_m{m}"


def _cache_dir() -> Path:
    override = os.environ.get("ORGANOID_CROFTON_CACHE")
    if override:
        return Path(override)
    base = os.environ.get("XDG_CACHE_HOME") or (Path.home() / ".cache")
    return Path(base) / "organoid_analysis" / "crofton_weights_v3"


def _solve_orbit_weights(ratios: tuple[float, float, float], m: int,
                         time_budget_s: float) -> tuple[np.ndarray, dict[str, Any]]:
    import time

    started = time.time()
    directions, labels, n_orbits = direction_orbits(ratios, m)
    scaled = directions * np.asarray(ratios, float)
    lengths = np.linalg.norm(scaled, axis=1)
    units = scaled / lengths[:, None]

    normals = fibonacci_sphere(20011)
    response_full = np.abs(normals @ units.T)
    # Collapse columns onto orbits: an orbit weight multiplies the SUM of
    # |n . u_k| over that orbit's members.
    response = np.zeros((len(normals), n_orbits))
    np.add.at(response.T, labels, response_full.T)

    active = np.arange(0, len(normals),
                       max(1, len(normals) // max(600, 4 * n_orbits)))
    objective = np.zeros(n_orbits + 1)
    objective[-1] = 1.0
    weights = None
    optimum = None
    for iteration in range(120):
        if time.time() - started > time_budget_s:
            raise TimeoutError(f"Crofton LP for ratios={ratios} m={m} exceeded "
                               f"{time_budget_s}s at iteration {iteration}")
        block = response[active]
        ones = np.ones((len(block), 1))
        solution = linprog(
            objective,
            A_ub=np.vstack([np.hstack([block, -ones]), np.hstack([-block, -ones])]),
            b_ub=np.concatenate([np.ones(len(block)), -np.ones(len(block))]),
            bounds=[(0, None)] * (n_orbits + 1), method="highs")
        if not solution.success:
            raise RuntimeError(f"Crofton LP failed for ratios={ratios} m={m}: "
                               f"{solution.message}")
        weights, optimum = solution.x[:n_orbits], float(solution.x[-1])
        deviation = np.abs(response @ weights - 1.0)
        if deviation.max() <= optimum + 1e-5:
            break
        active = np.union1d(active, np.argsort(deviation)[-max(200, n_orbits // 2):])

    if weights is None or optimum is None:  # pragma: no cover - loop always runs
        raise RuntimeError("Crofton LP produced no solution")
    expanded = weights[labels]
    achieved = np.abs(normals @ units.T) @ expanded
    diagnostics = {
        "n_directions": int(len(directions)), "n_orbits": int(n_orbits),
        "reduction_factor": float(len(directions) / n_orbits),
        "stencil_m": int(m), "scheme": SCHEME,
        "group_order": int(len(metric_symmetry_group(ratios))),
        "n_constraint": int(len(active)), "lp_optimum_t_star": optimum,
        "weight_sum": float(expanded.sum()),
        "plane_response_max_abs_dev": float(np.abs(achieved - 1).max()),
        "plane_response_rms_dev": float(np.sqrt(np.mean((achieved - 1) ** 2))),
        "n_zero_weights": int((expanded <= 1e-12).sum()),
        "solve_seconds": float(time.time() - started),
    }
    return expanded, diagnostics


@cache
def _weights_for_ratios(ratios: tuple[float, float, float], m: int
                        ) -> tuple[tuple[float, ...], tuple[tuple[str, Any], ...]]:
    """Dimensionless orbit-symmetric weights for a spacing ratio, with caching.

    Looked up first in the packaged tables, then in the on-disk cache, and
    solved only if absent from both. Weights depend on the spacing ratios and
    ``m`` alone, never on the object being measured.
    """
    name = _cache_key(ratios, m) + ".json"
    for directory, origin in ((PACKAGED_WEIGHTS, "packaged"), (_cache_dir(), "cache")):
        path = directory / name
        if path.is_file():
            stored = json.loads(path.read_text(encoding="utf-8"))
            diagnostics = dict(stored["diag"])
            diagnostics["weights_origin"] = origin
            diagnostics["weights_evidence_bearing"] = _weights_are_evidence_bearing(origin, ratios, m)
            return tuple(float(x) for x in stored["w"]), tuple(sorted(diagnostics.items()))

    weights, diagnostics = _solve_orbit_weights(ratios, m, LP_TIME_BUDGET_S)
    target = _cache_dir() / name
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps({"w": [float(x) for x in weights],
                                      "diag": diagnostics}), encoding="utf-8")
    except OSError:
        pass  # an unwritable cache slows later calls; it never changes a result
    diagnostics["weights_origin"] = "solved"
    diagnostics["weights_evidence_bearing"] = False
    return tuple(float(x) for x in weights), tuple(sorted(diagnostics.items()))


def crofton_weights_sym(spacing: ArrayLike, m: int) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Minimax weights for ``spacing`` and stencil radius ``m``.

    Returns ``(directions, w_times_area, diagnostics)``. The returned weights
    already carry the per-direction lattice area element, so that
    ``dot(transition_counts(mask, directions), w_times_area)`` is an area.
    """
    ratios = _ratios(spacing)
    dimensionless, diagnostic_items = _weights_for_ratios(ratios, int(m))
    diagnostics = dict(diagnostic_items)

    values = np.asarray(spacing, float)
    directions = np.asarray(primitive_directions(int(m)), int)
    lengths = np.linalg.norm(directions * values, axis=1)
    area_element = float(np.prod(values)) / lengths
    weights = np.asarray(dimensionless, float)
    # The LP is solved on ratios; rescale the grazing constant to physical units.
    scale = float(values.max())
    diagnostics["grazing_constant_D_um2"] = float(
        np.sum(weights * (lengths / scale) ** 2) / 24.0) * scale ** 2
    return directions, weights * area_element, diagnostics


def apriori_bound(spacing: ArrayLike, m: int, r_in: float) -> float:
    """Plane-response deviation plus grazing deficit, before any measurement."""
    diagnostics = crofton_weights_sym(spacing, m)[2]
    return float(diagnostics["plane_response_max_abs_dev"]
                 + diagnostics["grazing_constant_D_um2"] / r_in ** 2)


# ------------------------------------------------------------------ measure
def inscribed_radius(mask: np.ndarray, spacing: ArrayLike) -> float:
    """Largest inscribed sphere radius, from the mask and spacing alone.

    The local feature size, not the volume-equivalent radius: the latter
    overstates resolution for elongated or creased objects, which is what made
    an earlier volume-equivalent domain declaration fail its confirmation set.
    """
    # Instance adapters crop to the bounding box. The exterior is background,
    # including when a cropped box contains only foreground; EDT otherwise
    # sees no exterior on some faces and overestimates the inscribed radius.
    return float(ndi.distance_transform_edt(np.pad(mask, 1), sampling=spacing).max())


def measure(mask: np.ndarray, spacing: tuple[float, float, float], *, force_m: int | None = None) -> dict[str, Any]:
    """Surface area, volume and the domain variables for a 3D mask.

    Applies no domain gating: a measurement is always returned, and
    ``in_domain()`` decides separately whether it is inside the qualified
    domain. Every quantity here is computable from the mask and the voxel
    spacing alone -- no truth value, no shape class and no fitted constant
    enters a measurement.
    """
    binary = np.asarray(mask, bool)
    values = tuple(float(s) for s in spacing)
    n_voxels = int(np.count_nonzero(binary))
    if n_voxels == 0:
        return {"method": METHOD_NAME, "surface_area": 0.0, "volume": 0.0,
                "r_in": 0.0, "rho_in": 0.0, "anisotropy": max(values) / min(values),
                "stencil_m": None, "apriori_bound": None}

    volume = n_voxels * float(np.prod(values))
    r_in = inscribed_radius(binary, values)
    m = (int(force_m) if force_m is not None
         else min(M_CANDIDATES, key=lambda candidate: apriori_bound(values, candidate, r_in)))
    bound = apriori_bound(values, m, r_in)

    directions, weights, diagnostics = crofton_weights_sym(values, m)
    area = float(np.dot(transition_counts(binary, directions), weights))

    return {"method": METHOD_NAME, "surface_area": area, "volume": volume,
            "r_in": r_in, "rho_in": r_in / max(values),
            "anisotropy": max(values) / min(values),
            "stencil_m": m, "n_directions": diagnostics["n_directions"],
            "n_orbits": diagnostics["n_orbits"],
            "plane_response_max_abs_dev": diagnostics["plane_response_max_abs_dev"],
            "weights_origin": diagnostics.get("weights_origin", "packaged"),
            "weights_evidence_bearing": bool(diagnostics.get("weights_evidence_bearing", False)),
            "apriori_bound": bound}


def in_domain(measurement: dict[str, Any]) -> tuple[bool, tuple[str, ...]]:
    """Apply the machine-checkable part of the frozen applicability domain.

    Does NOT and cannot check ``DOMAIN_SCOPE``: a True result means the object
    clears the resolution and anisotropy gate, not that it is smooth.
    """
    reasons: list[str] = []
    if measurement["rho_in"] < DOMAIN_RHO_IN_MIN:
        reasons.append(f"rho_in {measurement['rho_in']:.2f} < {DOMAIN_RHO_IN_MIN}")
    if measurement["anisotropy"] > DOMAIN_ANISO_MAX:
        reasons.append(f"anisotropy {measurement['anisotropy']:.2f} > {DOMAIN_ANISO_MAX}")
    return (not reasons), tuple(reasons)


def sphericity(volume: float, area: float) -> float:
    return float(np.pi ** (1 / 3) * (6 * volume) ** (2 / 3) / area)
