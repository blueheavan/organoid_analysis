"""SG-2a volume-domain DEVELOPMENT study (2026-09-26). Not confirmation evidence.

Question
--------
The D-2 volume domain (``rho_in >= 10``, anisotropy <= 4, smooth closed) was
challenged by an in-domain pair of lattice-centred spheres with identical
binary masks (docs/SCIENTIFIC_GATE_REASSESSMENT_2026-09-26.md). Every earlier
phantom grid drew sub-voxel offsets uniformly at random, so no grid ever
sampled the lattice-symmetric placements where voxel-count error is largest.
This study measures voxel-count volume error with those placements included,
and derives a volume-specific resolution threshold.

Why voxel count is kept (no estimator change)
---------------------------------------------
Under a uniformly random sub-voxel translation the voxel count times the voxel
volume is a design-unbiased volume estimator for any bounded shape (Cavalieri
point counting), so its expected error is zero and only its dispersion depends
on resolution. A smoothing or surface-reconstruction volume would trade this
unbiasedness for a shape-dependent bias, and the counterexample proves that no
mask-only estimator can remove the ambiguity in the worst case anyway. The
estimator is therefore held fixed and only the domain is derived here.

Domain variable (declared before this script was run)
-----------------------------------------------------
``rho_vol = r_in / cbrt(sz * sy * sx)``: inscribed radius in units of the
geometric-mean voxel edge. The lattice-point discrepancy of a count is governed
by how many voxels resolve the object, not by the coarsest axis alone; the
surface variable ``rho_in = r_in / max(spacing)`` over-penalises anisotropic
grids for volume (an exploratory scan of lattice-centred spheres before this
study showed errors at a given rho_in falling several-fold from isotropic to
4:1 spacing). The volume domain keeps the surface gate as well:
``rho_in >= 10`` and anisotropy <= 4 still apply.

Selection rule (declared before this script was run)
----------------------------------------------------
``RHO_VOL_MIN`` = the smallest ladder value L such that every development case
whose MEASURED rho_vol >= L (and which clears the surface gate) has
``|volume_rel_err| <= 0.5 %`` -- a margin of 2x to the 1 % criterion -- AND the
dense lattice-centred sphere scan (5 symmetric placements, isotropic and
anisotropic spacings) has sup |error| <= 0.5 % over rho_vol in [L, 40].
If no ladder value satisfies this, the study reports NO DOMAIN.

Nothing from the confirmation grid of the later record is generated, executed
or consulted here.
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "docs/evidence/2026-09-12-surface-method-development/code"))
sys.path.insert(0, str(ROOT / "src"))

import common2 as c2  # noqa: E402  frozen 2026-09-12 phantom generator and oracles

from organoid_analysis.quantification.surface_crofton import inscribed_radius  # noqa: E402

MARGIN_LIMIT = 0.005
LADDER = (10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36)
SPACINGS = ((1.0, 1.0, 1.0), (1.5, 1.0, 1.0), (2.0, 1.0, 1.0), (3.0, 1.0, 1.0),
            (4.0, 1.0, 1.0), (2.5, 0.8, 0.8))
SHAPES = ("sphere", "ellipsoid", "capsule", "torus")
ASPECT = c2.GRIDS2["dev2"]["aspect"]
# Lattice-symmetric placements of the object centre (fractions of a voxel):
# voxel centre, voxel corner, face centre, edge centre, and the Z-only corner.
SYMMETRIC_OFFSETS = {"centre": (0.0, 0.0, 0.0), "corner": (0.5, 0.5, 0.5), "face": (0.0, 0.5, 0.5),
                     "edge": (0.0, 0.0, 0.5), "zcorner": (0.5, 0.0, 0.0)}
N_RANDOM_POSES = 3
SEED = 20260926


def gm(spacing) -> float:
    return float(np.prod(spacing) ** (1 / 3))


def poses(rng: np.random.Generator):
    for name, offset in SYMMETRIC_OFFSETS.items():
        yield f"sym-{name}", np.eye(3), np.asarray(offset)
    rotations = Rotation.random(N_RANDOM_POSES, random_state=rng.integers(2**31)).as_matrix()
    for index, rotation in enumerate(rotations):
        yield f"rand{index}", rotation, rng.uniform(-0.5, 0.5, 3)


def phantom_grid():
    rng = np.random.default_rng(SEED)
    for spacing in SPACINGS:
        for rho_vol in LADDER:
            f = rho_vol * gm(spacing)
            for kind in SHAPES:
                params = c2.shape_params(kind, f, ASPECT)
                for pose, rotation, offset in poses(rng):
                    yield kind, params, spacing, rho_vol, pose, rotation, offset


def sphere_count(radius: float, spacing, offset) -> int:
    """Exact lattice-point count of a ball centred at ``-offset`` voxels (sign test at voxel centres)."""
    sp = np.asarray(spacing, float)
    n = np.ceil(radius / sp).astype(int) + 2
    sq = [((np.arange(-k, k + 1) + o) * s) ** 2 for k, o, s in zip(n, offset, sp)]
    yx = np.add.outer(sq[1], sq[2]).ravel()
    yx.sort()
    return int(sum(np.searchsorted(yx, radius * radius - zz, side="right") for zz in sq[0] if zz <= radius * radius))


def dense_sphere_scan(out_csv: Path, step: float = 1 / 40, top: float = 40.0) -> list[dict]:
    rows = []
    for spacing in SPACINGS:
        g = gm(spacing)
        for rho_vol in np.arange(LADDER[0], top + 1e-9, step):
            radius = rho_vol * g
            truth = 4 / 3 * np.pi * radius ** 3
            for name, offset in SYMMETRIC_OFFSETS.items():
                count = sphere_count(radius, spacing, offset)
                rows.append({"spacing": "x".join(f"{s:g}" for s in spacing), "rho_vol": round(float(rho_vol), 6),
                             "rho_in": radius / max(spacing), "placement": name,
                             "volume_rel_err": count * float(np.prod(spacing)) / truth - 1})
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def _measure(item):
    index, (kind, params, spacing, rho_vol, pose, rotation, offset) = item
    mask = c2.rasterize2(kind, params, spacing, rotation, offset)
    _, v_true = c2.truth2(kind, params)
    r_in = inscribed_radius(mask, spacing)
    n = int(mask.sum())
    v_est = n * float(np.prod(spacing))
    return {"case": f"voldev-{index:04d}", "shape": kind, "spacing": "x".join(f"{s:g}" for s in spacing),
            "anisotropy": max(spacing) / min(spacing), "rho_vol_target": rho_vol, "pose": pose,
            "r_in": r_in, "rho_in": r_in / max(spacing), "rho_vol": r_in / gm(spacing),
            "n_voxels": n, "volume_true": float(v_true), "volume_est": v_est,
            "volume_rel_err": v_est / float(v_true) - 1}


def run_grid(out_csv: Path, processes: int = 4) -> list[dict]:
    """Cases are generated sequentially (fixed RNG order) and measured in parallel, order preserved."""
    from multiprocessing import Pool
    items = list(enumerate(phantom_grid()))
    rows = []
    t0 = time.time()
    with out_csv.open("w", newline="") as handle, Pool(processes) as pool:
        writer = None
        for row in pool.imap(_measure, items, chunksize=2):
            if writer is None:
                writer = csv.DictWriter(handle, fieldnames=list(row))
                writer.writeheader()
            writer.writerow(row)
            rows.append(row)
            if len(rows) % 100 == 0:
                handle.flush()
                print(f"{len(rows):5d}/{len(items)} {time.time() - t0:7.0f}s", flush=True)
    return rows


def select_threshold(grid_rows: list[dict], scan_rows: list[dict]) -> dict:
    per_level = []
    chosen = None
    for level in LADDER:
        cases = [r for r in grid_rows if r["rho_vol"] >= level and r["rho_in"] >= 10.0 and r["anisotropy"] <= 4.0]
        scan = [r for r in scan_rows if r["rho_vol"] >= level]
        worst_case = max(cases, key=lambda r: abs(r["volume_rel_err"])) if cases else None
        worst_scan = max(scan, key=lambda r: abs(r["volume_rel_err"])) if scan else None
        ok = bool(cases and scan and abs(worst_case["volume_rel_err"]) <= MARGIN_LIMIT
                  and abs(worst_scan["volume_rel_err"]) <= MARGIN_LIMIT)
        per_level.append({"rho_vol_min": level, "n_grid_cases": len(cases),
                          "grid_worst_abs": abs(worst_case["volume_rel_err"]) if worst_case else None,
                          "grid_worst_case": worst_case["case"] if worst_case else None,
                          "scan_worst_abs": abs(worst_scan["volume_rel_err"]) if worst_scan else None,
                          "scan_worst_at": ({k: worst_scan[k] for k in ("spacing", "rho_vol", "placement")}
                                            if worst_scan else None),
                          "satisfies_rule": ok})
        if ok and chosen is None:
            chosen = level
    return {"rule": f"smallest ladder value with every in-gate dev case and the dense symmetric sphere scan "
                    f"<= {MARGIN_LIMIT:.1%}", "rho_vol_min": chosen, "per_level": per_level}


def _read(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return [{k: (v if k in ("case", "shape", "spacing", "pose", "placement") else float(v))
                 for k, v in row.items()} for row in csv.DictReader(handle)]


def summarise(grid: list[dict], scan: list[dict]) -> dict:
    result = select_threshold(grid, scan)
    by_pose: dict[str, float] = {}
    for r in grid:
        if r["rho_in"] >= 10 and r["anisotropy"] <= 4:
            key = "symmetric" if r["pose"].startswith("sym") else "random"
            by_pose[key] = max(by_pose.get(key, 0.0), abs(r["volume_rel_err"]))
    result["worst_in_D2_gate_by_pose_type"] = by_pose
    return result


if __name__ == "__main__":
    # `--select-only` re-derives the selection from the saved CSVs. It was used
    # once, after the first full run wrote both CSVs and then failed to
    # serialise a numpy bool in this JSON; the data and the rule were unchanged.
    if "--select-only" in sys.argv:
        scan = _read(HERE / "dense_sphere_scan.csv")
        grid = _read(HERE / "volume_dev.csv")
    else:
        scan = dense_sphere_scan(HERE / "dense_sphere_scan.csv")
        print(f"dense scan: {len(scan)} rows", flush=True)
        grid = run_grid(HERE / "volume_dev.csv")
    result = summarise(grid, scan)
    (HERE / "volume_dev_selection.json").write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps({k: result[k] for k in ("rho_vol_min", "worst_in_D2_gate_by_pose_type")}, indent=1))
