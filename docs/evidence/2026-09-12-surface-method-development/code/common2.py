"""Round-2 phantoms and grids.

Round 1 failed because the binding shape class (cylinders) had no in-domain
development case spanning the declared threshold, so the threshold was an
extrapolation exactly where it mattered.  Two structural fixes here:

1. **Cases are parameterised by the domain variable itself.**  Every phantom is
   built from a target ``rho_in`` (inscribed radius / max spacing) so that
   *every* shape class is sampled on the *same* dense ladder across the
   boundary region.  By construction r_inscribed == rho * max(spacing) for all
   six classes, so no class can be under-covered near the threshold.

2. **Three new shape classes** attacking specific assumptions -- box
   (maximally creased), capsule (elongated but smooth), torus (non-trivial
   topology, and R_eff / r_inscribed far from the blob value).  See
   out/ROUND2_PROTOCOL.md section 4.

Round-1 common.py is not modified; round-1 evidence stays reproducible.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial.transform import Rotation

from common import CRITERION, VOLUME_CRITERION, ellipsoid_area  # noqa: F401

SPHERICITY_CRITERION = 0.05

# --------------------------------------------------------------------------
# Grids.  dev2 and confirm2 are disjoint in seed, spacing, rho ladder and
# shape aspect ratios: confirm2 interleaves *between* dev2 rho values and uses
# different aspect ratios, so it probes the gaps rather than repeating dev2.
# --------------------------------------------------------------------------
GRIDS2 = {
    "dev2": {
        "seed": 20260912,
        "spacings": [(1.0, 1.0, 1.0), (2.0, 1.0, 1.0), (3.0, 1.0, 1.0),
                     (4.0, 1.0, 1.0), (2.5, 0.8, 0.8)],
        "rho_ladder": [3, 5, 6, 7, 8, 9, 10, 12, 16],
        "n_rotations": 3,          # identity + 2 random
        "n_sphere_offsets": 2,
        "aspect": {"ellipsoid": (2.0, 1.5, 1.0), "cylinder_h": 3.0,
                   "box": (2.0, 1.5, 1.0), "capsule_L": 3.0, "torus_R": 2.0},
    },
    "confirm2": {
        "seed": 777001,
        "spacings": [(1.2, 1.0, 1.0), (2.2, 1.0, 1.0), (3.5, 1.0, 1.0),
                     (5.0, 1.0, 1.0), (2.4, 0.6, 0.6)],
        "rho_ladder": [4, 6.5, 7.5, 8.5, 9.5, 11, 13, 15],
        "n_rotations": 2,          # identity + 1 random
        "n_sphere_offsets": 2,
        "aspect": {"ellipsoid": (2.5, 1.4, 1.0), "cylinder_h": 2.2,
                   "box": (1.3, 1.1, 1.0), "capsule_L": 4.5, "torus_R": 2.6},
    },
}

SHAPES2 = ("sphere", "ellipsoid", "cylinder", "box", "capsule", "torus")

# dev2b: development EXTENSION, added after dev2 was analysed and before any
# domain was declared.  dev2 showed the binding error is the crease deficit on
# box and cylinder, decaying ~1/rho and still ~4% at the top of the dev2
# ladder (rho=16).  Declaring a threshold there would repeat the round-1
# mistake -- a boundary sitting where coverage runs out.  dev2b extends the
# ladder for the two binding classes only, so the threshold is interpolated
# inside measured evidence rather than extrapolated past it.
# This is development data.  confirm2 is untouched.
GRIDS2["dev2b"] = {
    "seed": 20260912 + 1,
    "spacings": [(1.0, 1.0, 1.0), (2.0, 1.0, 1.0), (3.0, 1.0, 1.0),
                 (2.5, 0.8, 0.8)],
    "rho_ladder": [11, 13, 14, 18, 20, 24, 28, 34],
    "n_rotations": 3,
    "n_sphere_offsets": 2,
    "aspect": GRIDS2["dev2"]["aspect"],
    "shapes": ("cylinder", "box"),
}
# dev2c: the near-cube aspect, the worst crease-to-area ratio a box can have.
# dev2's box aspect is (2, 1.5, 1); a cube concentrates more edge length per
# unit area, so the threshold must be set against the cube, not against an
# elongated box.
GRIDS2["dev2c"] = {
    "seed": 20260912 + 2,
    "spacings": [(1.0, 1.0, 1.0), (2.0, 1.0, 1.0), (3.0, 1.0, 1.0)],
    "rho_ladder": [8, 11, 14, 18, 22, 28, 34],
    "n_rotations": 3,
    "n_sphere_offsets": 2,
    "aspect": {**GRIDS2["dev2"]["aspect"], "box": (1.0, 1.0, 1.0),
               "cylinder_h": 2.0},
    "shapes": ("box", "cylinder"),
}


def shape_params(kind: str, f: float, aspect: dict) -> dict:
    """Phantom parameters with inscribed radius exactly ``f``."""
    if kind == "sphere":
        return {"r": f}
    if kind == "ellipsoid":
        return {"axes": tuple(f * a for a in aspect["ellipsoid"])}
    if kind == "cylinder":
        return {"r": f, "h": aspect["cylinder_h"] * f}
    if kind == "box":
        return {"half": tuple(f * a for a in aspect["box"])}
    if kind == "capsule":
        return {"r": f, "L": aspect["capsule_L"] * f}
    if kind == "torus":
        return {"r": f, "R": aspect["torus_R"] * f}
    raise ValueError(kind)


def inscribed_radius_analytic(kind: str, p: dict) -> float:
    """The exact largest-inscribed-sphere radius, for verifying the mask one."""
    if kind == "sphere":
        return float(p["r"])
    if kind == "ellipsoid":
        return float(min(p["axes"]))
    if kind == "cylinder":
        return float(min(p["r"], p["h"] / 2))
    if kind == "box":
        return float(min(p["half"]))
    if kind == "capsule":
        return float(p["r"])
    return float(p["r"])          # torus: tube radius


# ------------------------------------------------------------------ oracles
def truth2(kind: str, p: dict) -> tuple[float, float]:
    """Closed-form (area, volume)."""
    if kind == "sphere":
        r = p["r"]
        return 4 * np.pi * r * r, 4 / 3 * np.pi * r ** 3
    if kind == "ellipsoid":
        a, b, c = p["axes"]
        return ellipsoid_area(a, b, c), 4 / 3 * np.pi * a * b * c
    if kind == "cylinder":
        r, h = p["r"], p["h"]
        return 2 * np.pi * r * h + 2 * np.pi * r * r, np.pi * r * r * h
    if kind == "box":
        a, b, c = p["half"]
        return 8 * (a * b + b * c + c * a), 8 * a * b * c
    if kind == "capsule":
        r, L = p["r"], p["L"]
        return 2 * np.pi * r * L + 4 * np.pi * r * r, \
            np.pi * r * r * L + 4 / 3 * np.pi * r ** 3
    R, r = p["R"], p["r"]
    return 4 * np.pi ** 2 * R * r, 2 * np.pi ** 2 * R * r * r


def oracle_self_check2(n: int = 400001) -> list[dict]:
    """Independent numerical check of every new closed form, before use.

    Capsule and torus are surfaces of revolution about axis 0, so both area and
    volume follow from one-dimensional quadrature of the profile; the box is
    checked against a direct face-by-face sum.  Required agreement: 1e-6.
    """
    checks = []

    # --- box: face sum and prism volume, independent of the collected formula
    for half in [(2.0, 1.5, 1.0), (13.0, 11.0, 10.0), (5.0, 5.0, 5.0)]:
        a, b, c = half
        faces = 2 * (2 * b) * (2 * c) + 2 * (2 * a) * (2 * c) + 2 * (2 * a) * (2 * b)
        vol = (2 * a) * (2 * b) * (2 * c)
        A, V = truth2("box", {"half": half})
        checks.append({"shape": "box", "p": half,
                       "area_rel_diff": abs(A / faces - 1),
                       "vol_rel_diff": abs(V / vol - 1)})

    # --- capsule: surface/volume of revolution, profile y(t) about axis 0
    for r, L in [(1.0, 3.0), (10.0, 45.0), (4.0, 4.0)]:
        t = np.linspace(-(L / 2 + r), (L / 2 + r), n)
        u = np.clip(np.abs(t) - L / 2, 0, None)
        y = np.sqrt(np.clip(r * r - u * u, 0, None))          # profile radius
        dy = np.gradient(y, t)
        A = np.trapezoid(2 * np.pi * y * np.sqrt(1 + dy ** 2), t)
        V = np.trapezoid(np.pi * y ** 2, t)
        Ac, Vc = truth2("capsule", {"r": r, "L": L})
        checks.append({"shape": "capsule", "p": (r, L),
                       "area_rel_diff": abs(Ac / A - 1),
                       "vol_rel_diff": abs(Vc / V - 1)})

    # --- torus: parametric double integral over (theta, phi)
    for R, r in [(2.0, 1.0), (26.0, 10.0), (1.8, 1.0)]:
        th = np.linspace(0, 2 * np.pi, 4001)
        A = np.trapezoid(2 * np.pi * r * (R + r * np.cos(th)), th)
        V = np.trapezoid(np.pi * r * r * (R + r * np.cos(th)), th)
        Ac, Vc = truth2("torus", {"R": R, "r": r})
        checks.append({"shape": "torus", "p": (R, r),
                       "area_rel_diff": abs(Ac / A - 1),
                       "vol_rel_diff": abs(Vc / V - 1)})

    bad = [c for c in checks if c["area_rel_diff"] > 1e-6 or c["vol_rel_diff"] > 1e-6]
    assert not bad, bad
    return checks


# ----------------------------------------------------------------- geometry
def implicit2(kind: str, p: dict, q: np.ndarray) -> np.ndarray:
    """Inside-negative exact signed distance (box/capsule/torus) in the object
    frame.  Axis 0 is the symmetry axis for capsule and torus."""
    if kind == "box":
        d = np.abs(q) - np.asarray(p["half"], float)
        outside = np.linalg.norm(np.maximum(d, 0.0), axis=-1)
        inside = np.minimum(d.max(-1), 0.0)
        return outside + inside
    if kind == "capsule":
        radial = np.linalg.norm(q[..., 1:], axis=-1)
        axial = np.clip(np.abs(q[..., 0]) - p["L"] / 2, 0, None)
        return np.sqrt(radial ** 2 + axial ** 2) - p["r"]
    if kind == "torus":
        radial = np.linalg.norm(q[..., 1:], axis=-1) - p["R"]
        return np.sqrt(radial ** 2 + q[..., 0] ** 2) - p["r"]
    from common import _implicit
    return _implicit(kind, p, q)


def bound2(kind: str, p: dict) -> float:
    if kind == "box":
        return float(np.linalg.norm(p["half"]))
    if kind == "capsule":
        return float(p["L"] / 2 + p["r"])
    if kind == "torus":
        return float(p["R"] + p["r"])
    from common import _bound
    return _bound(kind, p)


def rasterize2(kind: str, p: dict, spacing, rot: np.ndarray, offset,
               pad: int = 4) -> np.ndarray:
    """Binary mask by sign test at voxel centres (same convention as round 1)."""
    sp = np.asarray(spacing, float)
    half = np.ceil(bound2(kind, p) / sp).astype(int) + pad
    shape = tuple(int(x) for x in 2 * half + 1)
    axes = [(np.arange(n) - (h + o)) * s
            for n, h, o, s in zip(shape, half, offset, sp)]
    # Evaluate in slabs along axis 0: the largest confirm2 phantom is 35M
    # voxels, and a monolithic meshgrid of it costs several GB.
    mask = np.empty(shape, bool)
    slab = max(1, int(4e6 // (shape[1] * shape[2])))
    for lo in range(0, shape[0], slab):
        hi = min(lo + slab, shape[0])
        g = np.meshgrid(axes[0][lo:hi], axes[1], axes[2], indexing="ij")
        q = np.stack([x.ravel() for x in g], axis=1) @ rot
        mask[lo:hi] = (implicit2(kind, p, q) <= 0).reshape(hi - lo, *shape[1:])
    assert not (mask[0].any() or mask[-1].any() or mask[:, 0].any()
                or mask[:, -1].any() or mask[:, :, 0].any()
                or mask[:, :, -1].any()), f"phantom touches border: {kind} {p}"
    return mask


# -------------------------------------------------------------------- cases
def cases2(grid_name: str):
    g = GRIDS2[grid_name]
    rng = np.random.default_rng(g["seed"])
    rots = [np.eye(3)] + list(
        Rotation.random(g["n_rotations"] - 1, random_state=g["seed"]).as_matrix())
    for spacing in g["spacings"]:
        for rho in g["rho_ladder"]:
            f = rho * max(spacing)
            for kind in g.get("shapes", SHAPES2):
                p = shape_params(kind, f, g["aspect"])
                if kind == "sphere":
                    for i in range(g["n_sphere_offsets"]):
                        yield kind, p, spacing, f"offset{i}", np.eye(3), \
                            rng.uniform(-.5, .5, 3), rho
                else:
                    for j, rot in enumerate(rots):
                        yield kind, p, spacing, f"rot{j}", rot, \
                            rng.uniform(-.5, .5, 3), rho


def case_id2(grid: str, idx: int, kind: str, rho: float, spacing, variant: str) -> str:
    return (f"{grid}-{idx:04d}-{kind}-rho{rho:g}"
            f"-sp{'x'.join(f'{s:g}' for s in spacing)}-{variant}")


def case_records2(grid_name: str):
    for idx, (kind, p, spacing, variant, rot, offset, rho) in \
            enumerate(cases2(grid_name)):
        yield {"idx": idx, "grid": grid_name, "shape": kind, "params": p,
               "spacing": tuple(float(s) for s in spacing), "variant": variant,
               "rot": rot, "offset": offset, "rho_target": float(rho),
               "anisotropy": float(max(spacing) / min(spacing)),
               "case_id": case_id2(grid_name, idx, kind, rho, spacing, variant)}
