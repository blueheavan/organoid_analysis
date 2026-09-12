"""Round-2 step 2: solve symmetry-reduced minimax weights for every spacing
ratio used by dev2/confirm2, m = 1..6, and verify against round 1 at m <= 4.

Declared before execution: any (ratio, m) exceeding TIME_BUDGET_S is EXCLUDED
and recorded here; the cap is never lowered silently.
"""
import csv, time
import numpy as np
from estimators_v3 import crofton_weights_sym
from estimators_v2 import crofton_weights

TIME_BUDGET_S = 600.0
M_MAX = 6
SPACINGS = [(1.0,1.0,1.0),(1.2,1.0,1.0),(2.0,1.0,1.0),(2.2,1.0,1.0),
            (2.0,0.7,0.7),(2.5,0.8,0.8),(3.0,1.0,1.0),(3.5,1.0,1.0),
            (2.4,0.6,0.6),(4.0,1.0,1.0),(5.0,1.0,1.0)]

rows, excluded = [], []
for sp in SPACINGS:
    for m in range(1, M_MAX + 1):
        t0 = time.time()
        try:
            _, _, d3 = crofton_weights_sym(sp, m, time_budget_s=TIME_BUDGET_S)
        except TimeoutError as e:
            excluded.append({"spacing": str(sp), "m": m, "reason": str(e)})
            print(f"EXCLUDED sp={sp} m={m}: {e}", flush=True)
            continue
        r = {"spacing": "x".join(f"{s:g}" for s in sp),
             "aniso": round(max(sp)/min(sp), 4), "m": m,
             "n_directions": d3["n_directions"], "n_orbits": d3["n_orbits"],
             "group_order": d3["group_order"],
             "reduction": round(d3["reduction_factor"], 2),
             "t_star_v3": d3["lp_optimum_t_star"],
             "sup_dev_v3": d3["plane_response_max_abs_dev"],
             "converged_v3": bool(d3["plane_response_max_abs_dev"]
                                  <= d3["lp_optimum_t_star"] + 2e-5),
             "D_v3": d3["grazing_constant_D_um2"],
             "solve_s": round(d3["solve_seconds"], 2) if not d3["cached"] else 0.0}
        if m <= 4:
            try:
                _, _, d2 = crofton_weights(sp, m, "lp")
                r["t_star_v2"] = d2["lp_optimum_t_star"]
                r["sup_dev_v2"] = d2["plane_response_max_abs_dev"]
                r["t_star_rel_diff"] = abs(r["t_star_v3"]/d2["lp_optimum_t_star"] - 1)
                r["v3_no_worse"] = bool(r["sup_dev_v3"] <= r["sup_dev_v2"] * (1 + 1e-9))
            except Exception as e:
                r["t_star_v2"] = None; r["v3_no_worse"] = None
        rows.append(r)
        print(f"sp={sp} m={m} orb={d3['n_orbits']:3d} t*={r['t_star_v3']:.5e} "
              f"sup={r['sup_dev_v3']:.5e} conv={r['converged_v3']} "
              f"{r['solve_s']:.2f}s", flush=True)

keys = sorted({k for r in rows for k in r})
with open("out/lp_symmetry_verification.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)

print("\n=== SUMMARY ===")
print(f"solved {len(rows)}, excluded {len(excluded)}")
print(f"all v3 converged (sup <= t*+2e-5): {all(r['converged_v3'] for r in rows)}")
cmp = [r for r in rows if r.get("t_star_v2") is not None]
print(f"m<=4 comparisons: {len(cmp)}")
print(f"  max |t*_v3/t*_v2 - 1| = {max(r['t_star_rel_diff'] for r in cmp):.3e}")
print(f"  v3 sup_dev never worse than v2: {all(r['v3_no_worse'] for r in cmp)}")
worse = [r for r in cmp if not r["v3_no_worse"]]
print(f"  cases where v3 worse: {[(r['spacing'], r['m']) for r in worse]}")
better = [r for r in cmp if r["sup_dev_v3"] < r["sup_dev_v2"] * 0.999]
print(f"  cases where v3 STRICTLY better: {len(better)} -> "
      f"{[(r['spacing'], r['m']) for r in better]}")
print(f"  round-1 solves that had not converged (sup_v2 > t*_v2 + 2e-5): "
      f"{sum(1 for r in cmp if r['sup_dev_v2'] > r['t_star_v2'] + 2e-5)}/{len(cmp)}")
