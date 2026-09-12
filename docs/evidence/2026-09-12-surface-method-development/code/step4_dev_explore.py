"""Step 4 -- DEVELOPMENT-ONLY exploration of the candidate estimators.

Runs on the ``dev`` grid only.  The ``confirm`` grid is not imported, not
rasterized and not scored anywhere in this file.  Every design choice that the
freeze record will fix (which candidate family, which stencil radius rule,
which weighting scheme, the qualified applicability domain) is allowed to be
informed by what this script produces -- and nothing after the freeze is.

Estimators evaluated
  E0            production binary marching cubes, through features.geometry
  A1_sdf        continuous-field marching cubes on a signed distance field
  A2_occ        continuous-field marching cubes on an occupancy field
  Bvor_m{1..M}  Crofton, spherical-Voronoi weights  (m=1 is the project's
                already-rejected 13-direction estimator, carried as reference)
  Blp_m{1..M}   Crofton, minimax (LP) weights

Volume is voxel counting for every row (unchanged from production); sphericity
is recomputed per estimator because it inherits the area error.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
import time
from pathlib import Path

import numpy as np

import common
from estimators_v2 import (area_occupancy, area_signed_distance,
                           crofton_area_from_counts, crofton_weights,
                           primitive_directions, transition_counts)

sys.path.insert(0, str(common.REPO / "src"))
from organoid_analysis.quantification.features import geometry  # noqa: E402

OUT = Path(__file__).resolve().parent / "out"
FIELDS = ["case_id", "grid", "shape", "spacing", "aniso", "rho", "stratum", "estimator",
          "area_true", "area_est", "area_rel_err", "volume_true", "volume_est",
          "volume_rel_err", "sphericity_true", "sphericity_est", "sphericity_rel_err"]


def sphericity(volume: float, area: float) -> float:
    return math.pi ** (1 / 3) * (6 * volume) ** (2 / 3) / area


def run(grid: str, m_max: int, out_name: str) -> None:
    all_dirs = primitive_directions(m_max)
    records = list(common.case_records(grid))
    # warm the weight cache so per-case timing reflects estimation, not setup
    for sp in {c["spacing"] for c in records}:
        for m in range(1, m_max + 1):
            crofton_weights(sp, m, "voronoi")
            crofton_weights(sp, m, "lp")

    fh = open(OUT / out_name, "w", newline="")
    wr = csv.DictWriter(fh, fieldnames=FIELDS)
    wr.writeheader()
    t_start = time.time()
    for i, case in enumerate(records):
        sp, kind, prm = case["spacing"], case["shape"], case["params"]
        geo = (kind, prm, sp, case["rot"], case["offset"])
        mask = common.rasterize(*geo)
        a_true, v_true = common.truth(kind, prm)
        g, _ = geometry(mask, sp, fill_holes=False)
        v_est = float(g["volume_um3"])

        ests = {"E0": float(g["surface_area_um2"]),
                "A1_sdf": area_signed_distance(common.signed_field(*geo), sp),
                "A2_occ": area_occupancy(common.occupancy_field(*geo), sp)}
        counts = transition_counts(mask, all_dirs)
        for m in range(1, m_max + 1):
            for tag, scheme in (("Bvor", "voronoi"), ("Blp", "lp")):
                ests[f"{tag}_m{m}"] = crofton_area_from_counts(counts, all_dirs, sp, m, scheme)

        base = {"case_id": case["case_id"], "grid": grid, "shape": case["shape"],
                "spacing": "x".join(str(s) for s in sp),
                "aniso": round(max(sp) / min(sp), 4), "rho": case["rho"],
                "stratum": case["stratum"],
                "area_true": a_true, "volume_true": v_true, "volume_est": v_est,
                "volume_rel_err": v_est / v_true - 1,
                "sphericity_true": sphericity(v_true, a_true)}
        for name, a_est in ests.items():
            s_est = sphericity(v_est, a_est)
            wr.writerow({**base, "estimator": name, "area_est": a_est,
                         "area_rel_err": a_est / a_true - 1, "sphericity_est": s_est,
                         "sphericity_rel_err": s_est / base["sphericity_true"] - 1})
        fh.flush()
        if i % 10 == 0:
            print(f"  {i:4d}  {case['case_id']:<44s} {time.time()-t_start:7.1f}s", flush=True)
    fh.close()
    print(f"done {grid}: {time.time()-t_start:.1f}s -> {out_name}", flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--grid", default="dev")
    p.add_argument("--m-max", type=int, default=5)
    p.add_argument("--out", default=None)
    a = p.parse_args()
    run(a.grid, a.m_max, a.out or f"candidates_{a.grid}.csv")
