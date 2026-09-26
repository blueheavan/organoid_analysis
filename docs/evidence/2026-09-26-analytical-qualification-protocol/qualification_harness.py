"""Frozen harness for the SG-1a / SG-2a analytical qualification record.

Executes two confirmation grids through the PRODUCTION measurement function
``organoid_analysis.quantification.features.geometry`` (filled-envelope default,
the organoid route) and writes one row per phantom. It does not adjudicate:
``analytical_qualification_evidence.derive_results`` applies the frozen rules of
PROTOCOL.md to the raw CSV.

* ``confirm2`` -- the surface confirmation grid frozen on 2026-09-12
  (FREEZE_RECORD_V3.md) before it was first executed, regenerated here from the
  byte-identical frozen generator ``common2.py``. It was executed once with the
  development implementation; this run re-executes it with the production code.
* ``confirm3`` -- a new grid, frozen in PROTOCOL.md before its first execution.
  It adds the lattice-symmetric placements (object centre on a voxel centre,
  corner, face centre, edge centre, Z-only corner; identity orientation) that no
  earlier grid sampled and that carry the largest voxel-count error.

Usage: python qualification_harness.py --out <csv> [--processes N]
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[3]
GENERATOR_DIR = ROOT / "docs/evidence/2026-09-12-surface-method-development/code"
sys.path.insert(0, str(GENERATOR_DIR))

import common2 as c2  # noqa: E402  frozen 2026-09-12 generator (hash-checked by the record)

SMOOTH = ("sphere", "ellipsoid", "capsule", "torus")
CONFIRM3 = {
    "seed": 20260927,
    # Absolute spacings are new; their ratios are exactly the ratios whose
    # packaged Crofton weights were exercised by the 2026-09-12 confirmation, so
    # every surface value is produced by evidence-bearing weights.
    "spacings": [(0.4, 0.4, 0.4), (0.48, 0.4, 0.4), (0.8, 0.4, 0.4), (0.88, 0.4, 0.4),
                 (1.4, 0.4, 0.4), (1.6, 0.4, 0.4)],
    # Target inscribed radius in units of the geometric-mean voxel edge (rho_vol).
    "rho_vol_ladder": [11, 15, 20, 27, 33, 37, 42],
    "aspect": {"ellipsoid": (2.2, 1.3, 1.0), "capsule_L": 3.8, "torus_R": 2.3,
               "cylinder_h": 2.0, "box": (1.0, 1.0, 1.0)},
    "symmetric_offsets": {"centre": (0.0, 0.0, 0.0), "corner": (0.5, 0.5, 0.5), "face": (0.0, 0.5, 0.5),
                          "edge": (0.0, 0.0, 0.5), "zcorner": (0.5, 0.0, 0.0)},
    "n_random_poses": 2,
}
FIELDS = ["case_id", "grid", "shape", "smooth", "spacing", "anisotropy", "pose", "rho_target", "n_voxels",
          "area_true", "volume_true", "area_est", "volume_est", "area_rel_err", "volume_rel_err",
          "rho_in", "rho_vol", "stencil_m", "weights_evidence_bearing", "surface_in_qualified_domain",
          "volume_in_qualified_domain", "error"]


def gm(spacing) -> float:
    return float(np.prod(spacing) ** (1 / 3))


def confirm2_cases():
    for record in c2.case_records2("confirm2"):
        yield {"case_id": record["case_id"], "grid": "confirm2", "shape": record["shape"],
               "params": record["params"], "spacing": record["spacing"], "rot": record["rot"],
               "offset": record["offset"], "pose": record["variant"], "rho_target": record["rho_target"]}


def confirm3_cases():
    g = CONFIRM3
    rng = np.random.default_rng(g["seed"])
    index = 0
    for spacing in g["spacings"]:
        for rho in g["rho_vol_ladder"]:
            f = rho * gm(spacing)
            for kind in SMOOTH:
                params = c2.shape_params(kind, f, g["aspect"])
                poses = [(f"sym-{name}", np.eye(3), np.asarray(offset, float))
                         for name, offset in g["symmetric_offsets"].items()]
                rotations = Rotation.random(g["n_random_poses"], random_state=int(rng.integers(2**31))).as_matrix()
                poses += [(f"rand{k}", rotation, rng.uniform(-0.5, 0.5, 3)) for k, rotation in enumerate(rotations)]
                for pose, rot, offset in poses:
                    yield {"case_id": f"confirm3-{index:04d}-{kind}-rhovol{rho:g}-sp{'x'.join(f'{s:g}' for s in spacing)}-{pose}",
                           "grid": "confirm3", "shape": kind, "params": params, "spacing": tuple(spacing),
                           "rot": rot, "offset": offset, "pose": pose, "rho_target": float(rho)}
                    index += 1


def all_cases() -> list[dict]:
    return [*confirm2_cases(), *confirm3_cases()]


def measure_case(case: dict) -> dict:
    from organoid_analysis.quantification.features import geometry

    spacing = tuple(float(s) for s in case["spacing"])
    area_true, volume_true = (float(v) for v in c2.truth2(case["shape"], case["params"]))
    row = {"case_id": case["case_id"], "grid": case["grid"], "shape": case["shape"],
           "smooth": case["shape"] in SMOOTH, "spacing": "x".join(f"{s:g}" for s in spacing),
           "anisotropy": max(spacing) / min(spacing), "pose": case["pose"], "rho_target": case["rho_target"],
           "area_true": area_true, "volume_true": volume_true, "error": ""}
    try:
        mask = c2.rasterize2(case["shape"], case["params"], spacing, case["rot"], case["offset"])
        values, _ = geometry(mask, spacing)
    except Exception as error:  # recorded, never retried: a failed case voids the run (PROTOCOL.md V2)
        row["error"] = f"{type(error).__name__}: {error}"
        return row
    row.update(n_voxels=int(mask.sum()), area_est=values["surface_area_um2"], volume_est=values["volume_um3"],
               area_rel_err=values["surface_area_um2"] / area_true - 1,
               volume_rel_err=values["volume_um3"] / volume_true - 1,
               rho_in=values["surface_rho_in"], rho_vol=values["volume_rho_vol"],
               stencil_m=values["surface_stencil_radius"],
               weights_evidence_bearing=bool(values["surface_weights_evidence_bearing"]),
               surface_in_qualified_domain=bool(values["surface_in_qualified_domain"]),
               volume_in_qualified_domain=bool(values["volume_in_qualified_domain"]))
    return row


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    parser.add_argument("--processes", type=int, default=4)
    args = parser.parse_args(argv)
    c2.oracle_self_check2()
    cases = all_cases()
    print(f"{len(cases)} cases", flush=True)
    t0 = time.time()
    from multiprocessing import Pool
    with Pool(args.processes) as pool, open(args.out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        for k, row in enumerate(pool.imap(measure_case, cases, chunksize=1)):
            writer.writerow({key: (f"{value:.10g}" if isinstance(value, float) else value)
                             for key, value in row.items()})
            if k % 100 == 0:
                handle.flush()
                print(f"{k:5d}/{len(cases)} {time.time() - t0:7.0f}s", flush=True)
    print(f"done {len(cases)} cases in {time.time() - t0:.0f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
