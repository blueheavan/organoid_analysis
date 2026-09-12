"""Round-2 step 3+5: execute the dev2 grid, no domain applied.

Produces out/vv_dev2.csv with one row per case: the domain variables, the
candidate's area/volume/sphericity errors, and the legacy estimator E0 for
comparison.  The domain is chosen from this table afterwards (step 6) and is
NOT applied here -- gating before declaration is exactly the mistake round 1
must not repeat.
"""
import csv, sys, time
import numpy as np
from skimage import measure as skmeasure

import common2 as c2
import surface_area_v3 as v3



def e0_area(mask, spacing):
    """Legacy production estimator: marching cubes on the binary mask."""
    v, f, _, _ = skmeasure.marching_cubes(mask.astype(float), level=0.5,
                                          spacing=spacing, method="lewiner")
    return float(skmeasure.mesh_surface_area(v, f))


def run_grid(grid: str, out: str | None = None) -> str:
    """Execute one grid. IDENTICAL code path for development and confirmation."""
    out = out or f"out/vv_{grid}.csv"
    c2.oracle_self_check2()
    recs = list(c2.case_records2(grid))
    print(f"{grid}: {len(recs)} cases", flush=True)

    fields = ["case_id", "grid", "shape", "spacing", "anisotropy", "rho_target",
              "rho_in", "rho_eff", "r_in", "stencil_m", "n_directions", "n_orbits",
              "apriori_bound", "n_voxels",
              "area_true", "area_v3", "area_rel_err", "area_E0", "area_E0_rel_err",
              "volume_true", "volume_est", "volume_rel_err",
              "sphericity_true", "sphericity_v3", "sphericity_rel_err"]
    f = open(out, "w", newline="")
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()

    t0 = time.time()
    for i, r in enumerate(recs):
        mask = c2.rasterize2(r["shape"], r["params"], r["spacing"], r["rot"], r["offset"])
        A_true, V_true = c2.truth2(r["shape"], r["params"])
        m = v3.measure(mask, r["spacing"])
        aE0 = e0_area(mask, r["spacing"])
        sph_t = v3.sphericity(V_true, A_true)
        sph_e = v3.sphericity(m["volume"], m["surface_area"])
        w.writerow({
            "case_id": r["case_id"], "grid": grid, "shape": r["shape"],
            "spacing": "x".join(f"{s:g}" for s in r["spacing"]),
            "anisotropy": round(r["anisotropy"], 4), "rho_target": r["rho_target"],
            "rho_in": m["rho_in"], "rho_eff": m["rho_eff"], "r_in": m["r_in"],
            "stencil_m": m["stencil_m"], "n_directions": m["n_directions"],
            "n_orbits": m["n_orbits"], "apriori_bound": m["apriori_bound"],
            "n_voxels": int(mask.sum()),
            "area_true": A_true, "area_v3": m["surface_area"],
            "area_rel_err": m["surface_area"] / A_true - 1,
            "area_E0": aE0, "area_E0_rel_err": aE0 / A_true - 1,
            "volume_true": V_true, "volume_est": m["volume"],
            "volume_rel_err": m["volume"] / V_true - 1,
            "sphericity_true": sph_t, "sphericity_v3": sph_e,
            "sphericity_rel_err": sph_e / sph_t - 1})
        if i % 25 == 0:
            f.flush()
            print(f"  {i:4d}/{len(recs)}  {time.time()-t0:6.0f}s  "
                  f"{r['shape']:9s} rho={r['rho_target']:g}", flush=True)
    f.close()
    print(f"done {len(recs)} cases in {time.time()-t0:.0f}s -> {out}", flush=True)
    return out


if __name__ == "__main__":
    run_grid(sys.argv[1] if len(sys.argv) > 1 else "dev2")
