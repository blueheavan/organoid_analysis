"""Surface-area measurement V&V (see SURFACE_VV_PLAN.md, frozen before results).

Run from the repository root:
    pixi run python docs/evidence/2026-09-11-measurement-vv/surface_vv.py --grid dev
    pixi run python docs/evidence/2026-09-11-measurement-vv/surface_vv.py --grid confirm

E0 is measured through the production ``quantification.features.geometry``.
Candidates E1-E3 are implemented here, outside production code, so this study
cannot change production behavior. Oracles are closed-form expressions and do
not use any estimator.
"""
from __future__ import annotations

import argparse
import itertools
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import skimage
from scipy import integrate, ndimage as ndi, special
from scipy.spatial import SphericalVoronoi
from scipy.spatial.transform import Rotation
from skimage.measure import marching_cubes, mesh_surface_area

from organoid_analysis.quantification.features import geometry

OUT = Path(__file__).resolve().parent
CRITERION = 0.05
VOLUME_CRITERION = 0.01

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
        "n_rotations": 2,  # identity + 2 random = 3 orientations
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
    if kind == "sphere":
        r = p["r"]
        return 4 * np.pi * r * r, 4 / 3 * np.pi * r ** 3
    if kind == "ellipsoid":
        a, b, c = p["axes"]
        return ellipsoid_area(a, b, c), 4 / 3 * np.pi * a * b * c
    r, h = p["r"], p["h"]
    return 2 * np.pi * r * h + 2 * np.pi * r * r, np.pi * r * r * h


# ---------------------------------------------------------------- phantoms
def rasterize(kind: str, p: dict, spacing: tuple, rot: np.ndarray, offset: np.ndarray) -> np.ndarray:
    sp = np.asarray(spacing, float)
    if kind == "sphere":
        bound = p["r"]
    elif kind == "ellipsoid":
        bound = max(p["axes"])
    else:
        bound = float(np.hypot(p["r"], p["h"] / 2))
    half = np.ceil(bound / sp).astype(int) + 4
    shape = 2 * half + 1
    grids = np.meshgrid(*[(np.arange(n) - (h + o)) * s for n, h, o, s in zip(shape, half, offset, sp)],
                        indexing="ij")
    pts = np.stack([g.ravel() for g in grids], axis=1)
    q = pts @ rot  # object-frame coordinates: q = R^T p for row vectors
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


# -------------------------------------------------------------- estimators
def e0_production(mask: np.ndarray, spacing: tuple) -> float:
    return float(geometry(mask, spacing, fill_holes=False)[0]["surface_area_um2"])


def e1_gaussian(mask: np.ndarray, spacing: tuple, factor: float) -> float:
    sp = np.asarray(spacing, float)
    sigma_vox = factor * sp.max() / sp
    pad = np.ceil(4 * sigma_vox).astype(int) + 2
    field = ndi.gaussian_filter(np.pad(mask.astype(np.float64), [(int(k), int(k)) for k in pad]),
                                sigma_vox, mode="constant", cval=0.0, truncate=4.0)
    if field.max() <= 0.5:
        return float("nan")  # the object vanished below the isosurface level
    v, f, _, _ = marching_cubes(field, 0.5, spacing=tuple(sp), allow_degenerate=False)
    return float(mesh_surface_area(v, f))


def e2_signed_distance(mask: np.ndarray, spacing: tuple) -> float:
    sp = np.asarray(spacing, float)
    m = np.pad(mask, 2)
    field = ndi.distance_transform_edt(~m, sampling=sp) - ndi.distance_transform_edt(m, sampling=sp)
    v, f, _, _ = marching_cubes(field, 0.0, spacing=tuple(sp), allow_degenerate=False)
    return float(mesh_surface_area(v, f))


CROFTON_DIRS = np.array([k for k in itertools.product((-1, 0, 1), repeat=3)
                         if k > (0, 0, 0)])  # 13 primitive directions, one per +/- pair


def crofton_weights(spacing: tuple) -> tuple[np.ndarray, np.ndarray]:
    d = CROFTON_DIRS * np.asarray(spacing, float)
    u = d / np.linalg.norm(d, axis=1)[:, None]
    areas = SphericalVoronoi(np.vstack([u, -u])).calculate_areas()
    w = (areas[:13] + areas[13:]) / (4 * np.pi)
    return w, np.linalg.norm(d, axis=1)


def e3_crofton(mask: np.ndarray, spacing: tuple) -> float:
    m = np.pad(mask, 1)
    w, length = crofton_weights(spacing)
    voxel_volume = float(np.prod(spacing))
    area = 0.0
    for k, wi, li in zip(CROFTON_DIRS, w, length):
        a_sl, b_sl = [], []
        for step, n in zip(k, m.shape):
            if step == 1:
                a_sl.append(slice(0, n - 1)); b_sl.append(slice(1, n))
            elif step == -1:
                a_sl.append(slice(1, n)); b_sl.append(slice(0, n - 1))
            else:
                a_sl.append(slice(None)); b_sl.append(slice(None))
        crossings = int(np.count_nonzero(m[tuple(a_sl)] != m[tuple(b_sl)]))
        area += 2 * wi * crossings * voxel_volume / li
    return float(area)


ESTIMATORS = {
    "E0_marching_cubes_binary_v1": e0_production,
    "E1_gaussian_phys_sigma1.0": lambda m, s: e1_gaussian(m, s, 1.0),
    "E1a_gaussian_phys_sigma0.5_sensitivity": lambda m, s: e1_gaussian(m, s, 0.5),
    "E1b_gaussian_phys_sigma1.5_sensitivity": lambda m, s: e1_gaussian(m, s, 1.5),
    "E2_signed_distance": e2_signed_distance,
    "E3_crofton_13dir_voronoi": e3_crofton,
}


# ------------------------------------------------------------------ grid
def cases(grid: dict):
    rng = np.random.default_rng(grid["seed"])
    rots = [np.eye(3)] + list(Rotation.random(grid["n_rotations"], random_state=grid["seed"]).as_matrix())
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


def md_table(t: pd.DataFrame) -> str:
    """Minimal GitHub-markdown table (avoids an optional tabulate dependency)."""
    head = "| " + " | ".join(map(str, t.columns)) + " |"
    rule = "|" + "---|" * len(t.columns)
    body = ["| " + " | ".join(map(str, row)) + " |" for row in t.itertuples(index=False)]
    return "\n".join([head, rule, *body])


def rho_stratum(rho: float) -> str:
    return "very_small(<3)" if rho < 3 else "small(3-6)" if rho < 6 else "medium(6-12)" if rho < 12 else "large(>=12)"


def summarize(df: pd.DataFrame) -> dict:
    groupings = {"overall": [], "shape": ["shape"], "rho_stratum": ["rho_stratum"],
                 "spacing": ["spacing"], "anisotropy": ["anisotropy"],
                 "shape_x_rho": ["shape", "rho_stratum"], "spacing_x_rho": ["spacing", "rho_stratum"]}
    out: dict = {}
    for name, cols in groupings.items():
        rows = []
        for est, sub in df.groupby("estimator"):
            parts = [((), sub)] if not cols else sub.groupby(cols)
            for key, g in parts:
                err = g["area_rel_err"]
                finite = err.dropna()
                worst = g.loc[err.abs().fillna(np.inf).idxmax()]
                rows.append({"estimator": est, **({c: (key if isinstance(key, tuple) else (key,))[i]
                                                     for i, c in enumerate(cols)}),
                             "n": int(len(g)), "n_nonestimable": int(err.isna().sum()),
                             "mean_signed": float(finite.mean()), "median_signed": float(finite.median()),
                             "mean_abs": float(finite.abs().mean()), "median_abs": float(finite.abs().median()),
                             "max_abs": float(err.abs().fillna(np.inf).max()),
                             "min_signed": float(finite.min()), "max_signed": float(finite.max()),
                             "worst_case": worst["case_id"],
                             "status": "PASS" if err.notna().all() and err.abs().max() < CRITERION else "FAIL"})
        out[name] = rows
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", choices=sorted(GRIDS), required=True)
    ap.add_argument("--estimators", nargs="*", default=None,
                    help="subset of estimator ids (confirmation runs only the selected ones)")
    args = ap.parse_args()
    grid = GRIDS[args.grid]
    estimators = {k: v for k, v in ESTIMATORS.items() if args.estimators is None or k in args.estimators}

    # Oracle self-check before any estimator result is used.
    oracle_checks = []
    for axes in [(6, 4.5, 3), (48, 36, 24), (10, 6, 4.5), (40, 24, 18), (5, 5, 5), (9, 3, 3)]:
        closed, quad = ellipsoid_area(*axes), ellipsoid_area_quadrature(*axes)
        oracle_checks.append({"axes": axes, "closed_form": closed, "quadrature": quad,
                              "rel_diff": abs(closed / quad - 1)})
    assert all(c["rel_diff"] < 1e-8 for c in oracle_checks), oracle_checks
    assert abs(ellipsoid_area(5, 5, 5) / (4 * np.pi * 25) - 1) < 1e-12

    records = []
    for idx, (kind, p, spacing, variant, rot, offset, min_axis) in enumerate(cases(grid)):
        mask = rasterize(kind, p, spacing, rot, offset)
        a_true, v_true = truth(kind, p)
        v_meas = float(mask.sum() * np.prod(spacing))
        size = p.get("r", p.get("axes", [None])[0]) if kind != "ellipsoid" else p["axes"][0]
        case_id = f"{args.grid}-{idx:03d}-{kind}-size{size:g}-sp{'x'.join(f'{s:g}' for s in spacing)}-{variant}"
        rho = min_axis / max(spacing)
        for est, fn in estimators.items():
            t0 = time.perf_counter()
            a_est = fn(mask, spacing)
            dt = time.perf_counter() - t0
            sph_true = np.cbrt(np.pi) * (6 * v_true) ** (2 / 3) / a_true
            sph_meas = np.cbrt(np.pi) * (6 * v_meas) ** (2 / 3) / a_est if np.isfinite(a_est) else np.nan
            records.append({
                "case_id": case_id, "shape": kind, "size_um": size, "params": json.dumps(p),
                "spacing": "x".join(f"{s:g}" for s in spacing), "anisotropy": round(max(spacing) / min(spacing), 3),
                "variant": variant, "offset_vox": json.dumps([round(float(o), 6) for o in offset]),
                "rotation": json.dumps(np.round(rot, 9).tolist()),
                "min_semi_axis_um": min_axis, "rho": rho, "rho_stratum": rho_stratum(rho),
                "n_voxels": int(mask.sum()), "estimator": est,
                "area_true_um2": a_true, "area_est_um2": a_est, "area_rel_err": a_est / a_true - 1,
                "volume_true_um3": v_true, "volume_voxel_um3": v_meas, "volume_rel_err": v_meas / v_true - 1,
                "sphericity_true": sph_true, "sphericity_est": sph_meas,
                "sphericity_rel_err": sph_meas / sph_true - 1, "runtime_s": dt,
            })
        print(f"{case_id} done", flush=True)

    df = pd.DataFrame.from_records(records)
    df.to_csv(OUT / f"surface_vv_{args.grid}.csv", index=False, float_format="%.10g")
    summary = summarize(df)
    vol = df[df.estimator == df.estimator.iloc[0]]
    volume_summary = {
        "by_rho_stratum": [{"rho_stratum": k, "n": int(len(g)), "max_abs": float(g.volume_rel_err.abs().max()),
                            "median_signed": float(g.volume_rel_err.median()),
                            "status": "PASS" if g.volume_rel_err.abs().max() < VOLUME_CRITERION else "FAIL"}
                           for k, g in vol.groupby("rho_stratum")],
        "overall_max_abs": float(vol.volume_rel_err.abs().max()),
    }
    plane = {}
    rng = np.random.default_rng(0)
    normals = rng.normal(size=(20000, 3))
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    for spacing in grid["spacings"]:
        w, length = crofton_weights(spacing)
        u = CROFTON_DIRS * np.asarray(spacing) / length[:, None]
        resp = 2 * np.abs(normals @ u.T) @ w
        plane["x".join(f"{s:g}" for s in spacing)] = {"min": float(resp.min()), "max": float(resp.max()),
                                                      "mean": float(resp.mean())}
    meta = {"grid": args.grid, "grid_definition": grid, "criterion_abs_rel_area_error": CRITERION,
            "volume_criterion": VOLUME_CRITERION, "n_cases": int(df.case_id.nunique()),
            "estimators": list(estimators), "oracle_checks": oracle_checks,
            "crofton_plane_response": plane,
            "versions": {"python": sys.version.split()[0], "platform": platform.platform(),
                         "numpy": np.__version__, "scipy": scipy.__version__,
                         "scikit-image": skimage.__version__, "pandas": pd.__version__}}
    (OUT / f"surface_vv_{args.grid}_summary.json").write_text(
        json.dumps({"meta": meta, "area": summary, "volume": volume_summary}, indent=2, default=str))

    lines = [f"# Surface V&V — {args.grid} grid", "",
             f"{meta['n_cases']} cases × {len(estimators)} estimators. Criterion: every case |area error| < 5%.",
             "Error = A_est/A_true − 1. Generated by `surface_vv.py`; do not hand-edit.", ""]
    for name in ["overall", "shape", "rho_stratum", "anisotropy", "shape_x_rho"]:
        t = pd.DataFrame(summary[name])
        cols = [c for c in ["estimator", "shape", "rho_stratum", "anisotropy"] if c in t] + \
            ["n", "n_nonestimable", "mean_signed", "median_signed", "median_abs", "max_abs", "status", "worst_case"]
        t = t[cols].copy()
        for c in ["mean_signed", "median_signed", "median_abs", "max_abs"]:
            t[c] = (100 * t[c]).map(lambda x: f"{x:+.2f}%" if c != "max_abs" and c != "median_abs" else f"{x:.2f}%")
        lines += [f"## By {name}", "", md_table(t), ""]
    vt = pd.DataFrame(volume_summary["by_rho_stratum"])
    lines += ["## Voxel-count volume vs §9 <1% (independent of surface estimator)", "", md_table(vt), ""]
    (OUT / f"surface_vv_{args.grid}_summary.md").write_text("\n".join(lines))
    print(json.dumps({r["estimator"]: [r["status"], round(100 * r["max_abs"], 2)] for r in summary["overall"]}))


if __name__ == "__main__":
    main()
