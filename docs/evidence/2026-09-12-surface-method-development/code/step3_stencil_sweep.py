"""Step 3 -- plane-response bound of the Crofton candidates vs stencil radius.

Phantom-free: this computes, from lattice geometry alone, the supremum over
orientations of |R(n) - 1| for each (voxel spacing, stencil radius, weighting
scheme).  That supremum is a rigorous upper bound on the orientation bias of
the resulting estimator, available BEFORE any phantom is measured, and it is
what the frozen stencil-selection rule will be stated against.

Writes out/crofton_plane_response.csv incrementally.
"""
import csv
import sys
import time
from pathlib import Path

from estimators_v2 import crofton_weights

OUT = Path(__file__).resolve().parent / "out"
SPACINGS = [(1.0, 1.0, 1.0), (2.0, 1.0, 1.0), (3.0, 1.0, 1.0), (2.0, 0.7, 0.7),
            (1.5, 1.0, 1.0), (2.5, 0.8, 0.8), (4.0, 1.0, 1.0)]
FIELDS = ["spacing", "aniso", "m", "K", "lp_maxdev", "lp_t_star", "lp_nconstraint",
          "lp_nzero", "vor_maxdev", "secs"]

f = open(OUT / "crofton_plane_response.csv", "w", newline="")
wr = csv.DictWriter(f, fieldnames=FIELDS)
wr.writeheader()
for sp in SPACINGS:
    for m in range(1, 5):
        t0 = time.time()
        _, _, dl = crofton_weights(sp, m, "lp")
        _, _, dv = crofton_weights(sp, m, "voronoi")
        r = {"spacing": "x".join(map(str, sp)), "aniso": round(max(sp) / min(sp), 3),
             "m": m, "K": dl["n_directions"], "lp_maxdev": dl["plane_response_max_abs_dev"],
             "lp_t_star": dl["lp_optimum_t_star"], "lp_nconstraint": dl["n_constraint"],
             "lp_nzero": dl["n_zero_weights"],
             "vor_maxdev": dv["plane_response_max_abs_dev"], "secs": round(time.time() - t0, 1)}
        wr.writerow(r); f.flush()
        print(f"{r['spacing']:12s} aniso={r['aniso']:4.2f} m={m} K={r['K']:5d} "
              f"lp={r['lp_maxdev']*100:6.3f}% t*={r['lp_t_star']*100:6.3f}% "
              f"vor={r['vor_maxdev']*100:6.3f}%  {r['secs']:6.1f}s", flush=True)
        # stop refining once the orientation bias is negligible against the 5% budget
        if r["lp_maxdev"] < 0.005 and m >= 3:
            break
f.close()
