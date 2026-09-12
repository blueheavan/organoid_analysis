"""Shared phantoms, analytical oracles and case grids for the surface-area V&V.

The oracle expressions, the rasterizer and the two case grids (``dev`` and
``confirm``) are reproduced *verbatim* in semantics from the frozen study
``docs/evidence/2026-09-11-analytical-geometry-record/surface_vv.py`` so that

  * the legacy baseline E0 is re-measured on exactly the same cases, and
  * the ``confirm`` grid -- which that frozen script defines but which has
    never been executed -- remains an untouched confirmation set.

Nothing in this module depends on any estimator.
"""
from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
from scipy import integrate, special
from scipy.spatial.transform import Rotation

REPO = Path("/Users/jlwang/Desktop/organoid_analysis")
FROZEN_RECORD = REPO / "docs/evidence/2026-09-11-analytical-geometry-record"

CRITERION = 0.05          # SCIENTIFIC_SPEC section 9: |A_est/A_true - 1| < 5% per case
VOLUME_CRITERION = 0.01   # SCIENTIFIC_SPEC section 9: |V_est/V_true - 1| < 1% per case

GRIDS = {
    "dev": {
        "seed": 20260911,
        "spacings": [(1.0, 1.0, 1.0), (2.0, 1.0, 1.0), (3.0, 1.0, 1.0), (2.0, 0.7, 0.7)],
        "sphere_radii": [3, 6, 12, 24, 36], "sphere_offsets": 3,
        "ellipsoid_scales": [6, 12, 24, 48], "ellipsoid_ratio": (1.0, 0.75, 0.5),
        "cylinder_radii": [3, 6, 12, 24], "cylinder_height_ratio": 3.0,
        "n_rotations": 3,
    },
    "confirm": {
        "seed": 424242,
        "spacings": [(1.5, 1.0, 1.0), (2.5, 0.8, 0.8), (4.0, 1.0, 1.0)],
        "sphere_radii": [5, 10, 20, 30], "sphere_offsets": 2,
        "ellipsoid_scales": [10, 20, 40], "ellipsoid_ratio": (1.0, 0.6, 0.45),
        "cylinder_radii": [5, 10, 20], "cylinder_height_ratio": 2.5,
        "n_rotations": 2,
    },
}


# ----------------------------------------------------------------- oracles
def ellipsoid_area(a: float, b: float, c: float) -> float:
    """Closed-form ellipsoid area (Legendre incomplete elliptic integrals)."""
    a, b, c = sorted((a, b, c), reverse=True)
    if np.isclose(a, c):
        return 4 * np.pi * a * a
    phi = np.arccos(c / a)
    m = (a * a * (b * b - c * c)) / (b * b * (a * a - c * c))
    e_inc, f_inc = special.ellipeinc(phi, m), special.ellipkinc(phi, m)
    return float(2 * np.pi * c * c + 2 * np.pi * a * b / np.sin(phi)
                 * (e_inc * np.sin(phi) ** 2 + f_inc * np.cos(phi) ** 2))


def ellipsoid_area_quadrature(a: float, b: float, c: float) -> float:
    """Independent check: integrate |r_theta x r_phi| over the parametric surface."""
    def element(theta: float, phi: float) -> float:
        st, ct, sp, cp = np.sin(theta), np.cos(theta), np.sin(phi), np.cos(phi)
        return st * np.sqrt((b * c * st * cp) ** 2 + (a * c * st * sp) ** 2 + (a * b * ct) ** 2)
    value, _ = integrate.dblquad(element, 0, 2 * np.pi, 0, np.pi, epsabs=0, epsrel=1e-12)
    return float(value)


def truth(kind: str, p: dict) -> tuple[float, float]:
    """Return (analytical area, analytical volume) for a phantom specification."""
    if kind == "sphere":
        r = p["r"]
        return 4 * np.pi * r * r, 4 / 3 * np.pi * r ** 3
    if kind == "ellipsoid":
        a, b, c = p["axes"]
        return ellipsoid_area(a, b, c), 4 / 3 * np.pi * a * b * c
    r, h = p["r"], p["h"]
    return 2 * np.pi * r * h + 2 * np.pi * r * r, np.pi * r * r * h


def oracle_self_check() -> list[dict]:
    """Closed form vs numerical quadrature; must agree before any estimator is used."""
    checks = []
    for axes in [(6, 4.5, 3), (48, 36, 24), (10, 6, 4.5), (40, 24, 18), (5, 5, 5), (9, 3, 3)]:
        closed, quad = ellipsoid_area(*axes), ellipsoid_area_quadrature(*axes)
        checks.append({"axes": axes, "closed_form": closed, "quadrature": quad,
                       "rel_diff": abs(closed / quad - 1)})
    assert all(c["rel_diff"] < 1e-8 for c in checks), checks
    assert abs(ellipsoid_area(5, 5, 5) / (4 * np.pi * 25) - 1) < 1e-12
    return checks


# ---------------------------------------------------------------- phantoms
def _implicit(kind: str, p: dict, q: np.ndarray) -> np.ndarray:
    """Inside-negative implicit function, evaluated on object-frame points ``q``.

    Returned in units of length so its zero level set is the true boundary and
    its gradient magnitude is ~1 near that boundary (an approximate signed
    distance).  Used both to rasterize the binary mask (sign test) and to build
    the Candidate-A continuous field, so the two are guaranteed consistent.
    """
    if kind == "sphere":
        return np.linalg.norm(q, axis=-1) - p["r"]
    if kind == "ellipsoid":
        axes = np.asarray(p["axes"], float)
        # Normalized algebraic distance rescaled by the local axis length; exact
        # sign, and a smooth monotone function of the true signed distance.
        f = np.sqrt(((q / axes) ** 2).sum(-1))
        scale = np.sqrt(((q / axes**2) ** 2).sum(-1))
        with np.errstate(divide="ignore", invalid="ignore"):
            d = np.where(scale > 0, (f - 1) / np.where(scale > 0, scale, 1.0), -min(axes))
        return d
    # Finite cylinder: exact signed distance to the capped cylinder.
    r, h = p["r"], p["h"]
    radial = np.linalg.norm(q[..., 1:], axis=-1) - r
    axial = np.abs(q[..., 0]) - h / 2
    outside = np.sqrt(np.maximum(radial, 0) ** 2 + np.maximum(axial, 0) ** 2)
    inside = np.minimum(np.maximum(radial, axial), 0.0)
    return outside + inside


def _bound(kind: str, p: dict) -> float:
    if kind == "sphere":
        return float(p["r"])
    if kind == "ellipsoid":
        return float(max(p["axes"]))
    return float(np.hypot(p["r"], p["h"] / 2))


def voxel_centers(kind: str, p: dict, spacing, offset, pad: int = 4):
    """Physical-frame voxel-centre coordinate grid enclosing the phantom.

    Convention (identical to the frozen study): index ``i`` along an axis has
    physical coordinate ``(i - half - offset) * spacing``, i.e. voxel *centres*
    sample the continuum, and ``offset`` shifts the object by a subvoxel amount
    relative to the lattice.
    """
    sp = np.asarray(spacing, float)
    half = np.ceil(_bound(kind, p) / sp).astype(int) + pad
    shape = 2 * half + 1
    axes = [(np.arange(n) - (h + o)) * s for n, h, o, s in zip(shape, half, offset, sp)]
    grids = np.meshgrid(*axes, indexing="ij")
    return np.stack(grids, axis=-1), tuple(int(n) for n in shape)


def rasterize(kind: str, p: dict, spacing, rot: np.ndarray, offset) -> np.ndarray:
    """Binary mask by inside/outside test at voxel centres (frozen-study semantics)."""
    pts, shape = voxel_centers(kind, p, spacing, offset)
    q = pts.reshape(-1, 3) @ rot          # object frame: q = R^T p for row vectors
    if kind == "sphere":
        inside = (q ** 2).sum(1) <= p["r"] ** 2
    elif kind == "ellipsoid":
        inside = ((q / np.asarray(p["axes"])) ** 2).sum(1) <= 1
    else:
        inside = (np.abs(q[:, 0]) <= p["h"] / 2) & (q[:, 1] ** 2 + q[:, 2] ** 2 <= p["r"] ** 2)
    mask = inside.reshape(shape)
    assert not (mask[0].any() or mask[-1].any() or mask[:, 0].any() or mask[:, -1].any()
                or mask[:, :, 0].any() or mask[:, :, -1].any()), "phantom touches border"
    return mask


def signed_field(kind: str, p: dict, spacing, rot: np.ndarray, offset) -> np.ndarray:
    """Continuous inside-negative field sampled at the same voxel centres.

    This is the *idealised* Candidate-A input: a continuous segmentation score
    that is an exact monotone function of the signed distance to the boundary.
    A real probability map is not guaranteed to be this (see the report).
    """
    pts, shape = voxel_centers(kind, p, spacing, offset)
    q = pts.reshape(-1, 3) @ rot
    return _implicit(kind, p, q).reshape(shape)


def occupancy_field(kind: str, p: dict, spacing, rot: np.ndarray, offset,
                    sub: int = 4) -> np.ndarray:
    """Partial-volume field: fraction of each voxel inside the object.

    Estimated by ``sub**3`` regular sub-samples per voxel.  This is the field a
    perfectly calibrated partial-volume segmentation would output, and is a
    different (non-distance) continuous convention from ``signed_field``.
    """
    sp = np.asarray(spacing, float)
    _, shape = voxel_centers(kind, p, spacing, offset)
    half = (np.asarray(shape) - 1) // 2
    off = (np.arange(sub) + 0.5) / sub - 0.5        # sub-voxel sample offsets
    acc = np.zeros(shape, float)
    for dz in off:
        for dy in off:
            for dx in off:
                d = np.array([dz, dy, dx])
                axes = [(np.arange(n) - (h + o) + t) * s
                        for n, h, o, s, t in zip(shape, half, offset, sp, d)]
                g = np.meshgrid(*axes, indexing="ij")
                q = np.stack([x.ravel() for x in g], axis=1) @ rot
                acc += (_implicit(kind, p, q) <= 0).reshape(shape)
    return acc / sub ** 3


# ------------------------------------------------------------------ grid
def cases(grid_name: str):
    """Yield the frozen case list for a grid (identical ordering to the frozen study)."""
    grid = GRIDS[grid_name]
    rng = np.random.default_rng(grid["seed"])
    rots = [np.eye(3)] + list(
        Rotation.random(grid["n_rotations"], random_state=grid["seed"]).as_matrix())
    for spacing in grid["spacings"]:
        for r in grid["sphere_radii"]:
            for i in range(grid["sphere_offsets"]):
                yield "sphere", {"r": r}, spacing, f"offset{i}", np.eye(3), rng.uniform(-.5, .5, 3), r
        for s in grid["ellipsoid_scales"]:
            axes = tuple(s * f for f in grid["ellipsoid_ratio"])
            for j, rot in enumerate(rots):
                yield "ellipsoid", {"axes": axes}, spacing, f"rot{j}", rot, rng.uniform(-.5, .5, 3), min(axes)
        for r in grid["cylinder_radii"]:
            h = grid["cylinder_height_ratio"] * r
            for j, rot in enumerate(rots):
                yield "cylinder", {"r": r, "h": h}, spacing, f"rot{j}", rot, rng.uniform(-.5, .5, 3), r


def case_records(grid_name: str):
    """``cases()`` as dicts, with the frozen rho convention applied.

    A pure view over ``cases()``: same order, same content, no new randomness.
    rho = (smallest defining length) / (largest voxel spacing), exactly as in
    the frozen study and as re-verified by step1_baseline.py.
    """
    for idx, (kind, p, spacing, variant, rot, offset, min_axis) in enumerate(cases(grid_name)):
        rho = min_axis / max(spacing)
        yield {"idx": idx, "grid": grid_name, "shape": kind, "params": p,
               "spacing": tuple(float(s) for s in spacing), "variant": variant,
               "rot": rot, "offset": offset, "rho": float(rho),
               "stratum": rho_stratum(rho),
               "case_id": case_id(grid_name, idx, kind, p, spacing, variant)}


def case_id(grid_name: str, idx: int, kind: str, p: dict, spacing, variant: str) -> str:
    size = p["axes"][0] if kind == "ellipsoid" else p["r"]
    return (f"{grid_name}-{idx:03d}-{kind}-size{size:g}"
            f"-sp{'x'.join(f'{s:g}' for s in spacing)}-{variant}")


def rho_stratum(rho: float) -> str:
    return ("very_small(<3)" if rho < 3 else "small(3-6)" if rho < 6
            else "medium(6-12)" if rho < 12 else "large(>=12)")


CROFTON_DIRS_13 = np.array([k for k in itertools.product((-1, 0, 1), repeat=3) if k > (0, 0, 0)])
