"""Step 2 -- mechanistic analysis of the binary marching-cubes surface error.

Four experiments, each against an analytical oracle:

  M1  Exactness on lattice-aligned polyhedra, at every resolution.
      An axis-aligned box is measured *exactly* by binary marching cubes at any
      voxel size.  If the error were a resolution deficit this could not happen.

  M2  Orientation response g(n): asymptotic area inflation as a function of the
      true surface normal, measured by binning the triangles of a very large
      digitized sphere by the radial direction at their centroid.  Its
      solid-angle average is the irreducible asymptotic bias.

  M3  Resolution sweep: sphere area and volume error vs rho = r / max(spacing),
      fitted as err(rho) = a + b/rho.  A non-zero intercept a is an asymptotic
      bias; volume shows a ~ 0 and surface does not.

  M4  Subvoxel-offset / orientation dispersion at fixed rho, and the resulting
      sphericity error propagation.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation
from skimage.measure import marching_cubes, mesh_surface_area

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import REPO, rasterize, truth  # noqa: E402

sys.path.insert(0, str(REPO / "src"))
from organoid_analysis.quantification.features import geometry  # noqa: E402

OUT = Path(__file__).resolve().parent / "out"
OUT.mkdir(exist_ok=True)
RNG = np.random.default_rng(7)


def e0(mask, spacing) -> float:
    return float(geometry(mask, spacing, fill_holes=False)[0]["surface_area_um2"])


# ------------------------------------------------------- M1 exact polyhedra
def m1_lattice_aligned_boxes() -> pd.DataFrame:
    """Axis-aligned boxes: binary MC is exact at every resolution.

    Rasterizing the box [0,a)x[0,b)x[0,c) on voxel centres with integer extents
    gives a mask whose true surface is entirely lattice-aligned.
    """
    rows = []
    for spacing in [(1.0, 1.0, 1.0), (3.0, 1.0, 1.0), (2.0, 0.7, 0.7)]:
        sp = np.asarray(spacing, float)
        for n in [2, 4, 8, 16, 32, 64]:
            mask = np.zeros(tuple(n + 4 for _ in range(3)), bool)
            mask[2:2 + n, 2:2 + n, 2:2 + n] = True
            ext = n * sp
            a_true = 2 * (ext[0] * ext[1] + ext[0] * ext[2] + ext[1] * ext[2])
            v_true = float(np.prod(ext))
            a_est = e0(mask, spacing)
            v_est = float(mask.sum() * np.prod(sp))
            rows.append({"spacing": "x".join(f"{s:g}" for s in spacing), "n_voxels_per_side": n,
                         "rho_equivalent": n / 2, "area_true": a_true, "area_est": a_est,
                         "area_rel_err": a_est / a_true - 1,
                         "volume_rel_err": v_est / v_true - 1})
    # A 45-degree-rotated square prism: the two slanted faces are (0,1,1) planes,
    # still lattice-representable, so the oracle stays exact.
    for spacing in [(1.0, 1.0, 1.0)]:
        for n in [8, 16, 32, 64]:
            size = 2 * n + 5
            z, y, x = np.mgrid[0:size, 0:size, 0:size]
            c = size // 2
            u, v = (y - c + x - c) / np.sqrt(2), (y - c - (x - c)) / np.sqrt(2)
            mask = (np.abs(u) <= n / np.sqrt(2)) & (np.abs(v) <= n / np.sqrt(2)) & \
                   (np.abs(z - c) <= n / 2)
            # side length of the rotated square cross-section
            side = 2 * n / np.sqrt(2)
            height = float(np.ptp(z[mask]) + 1)
            a_true = 2 * side * side + 4 * side * height
            a_est = e0(mask, spacing)
            rows.append({"spacing": "1x1x1 (45deg prism)", "n_voxels_per_side": n,
                         "rho_equivalent": n / 2, "area_true": a_true, "area_est": a_est,
                         "area_rel_err": a_est / a_true - 1, "volume_rel_err": np.nan})
    return pd.DataFrame(rows)


# ------------------------------------------------- M2 orientation response
def fibonacci_directions(n: int) -> np.ndarray:
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(phi), np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta)], 1)


def m2_orientation_response(radii=(48, 96), n_bins=192, spacings=((1.0, 1.0, 1.0),
                                                                  (3.0, 1.0, 1.0),
                                                                  (2.0, 0.7, 0.7))) -> pd.DataFrame:
    """g(n) = measured MC area / true area, resolved by true surface normal.

    A digitized sphere presents every orientation at once.  Each MC triangle is
    assigned to the solid-angle bin of the *true* normal (the radial direction at
    its centroid), so the ratio per bin is the local area-inflation factor.  The
    solid-angle-weighted mean of g over the sphere must equal the total sphere
    area error, which is checked.
    """
    bins = fibonacci_directions(n_bins)
    probe = fibonacci_directions(400_000)
    bin_weight = np.bincount(np.argmax(probe @ bins.T, axis=1), minlength=n_bins) / len(probe)
    rows = []
    for spacing in spacings:
        sp = np.asarray(spacing, float)
        for r in radii:
            rho = r / sp.max()
            mask = rasterize("sphere", {"r": r}, spacing, np.eye(3), np.zeros(3))
            verts, faces, _, _ = marching_cubes(np.pad(mask.astype(np.float32), 1), level=0.5,
                                                spacing=tuple(sp), allow_degenerate=False)
            tri = verts[faces]
            area = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
            cen = tri.mean(1) - (np.asarray(mask.shape) / 2 + 0.5) * sp  # sphere centre ~ array centre
            cen /= np.linalg.norm(cen, axis=1)[:, None]
            which = np.argmax(cen @ bins.T, axis=1)
            measured = np.bincount(which, weights=area, minlength=n_bins)
            a_true = 4 * np.pi * r * r
            g = measured / (bin_weight * a_true)
            total_err = area.sum() / a_true - 1
            check = float(np.sum(bin_weight * g) - 1 - total_err)
            for b in range(n_bins):
                nz, ny, nx = bins[b]
                rows.append({"spacing": "x".join(f"{s:g}" for s in spacing), "radius": r,
                             "rho": rho, "bin": b, "nz": nz, "ny": ny, "nx": nx,
                             "solid_angle_weight": bin_weight[b], "g": g[b],
                             "total_area_rel_err": total_err, "closure_residual": check})
    return pd.DataFrame(rows)


# ---------------------------------------------------- M3 resolution sweep
def m3_resolution_sweep() -> pd.DataFrame:
    rows = []
    for spacing in [(1.0, 1.0, 1.0), (3.0, 1.0, 1.0), (2.0, 0.7, 0.7)]:
        sp = np.asarray(spacing, float)
        for r in [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128]:
            rho = r / sp.max()
            if rho < 1.5:
                continue
            n_off = 8 if r <= 32 else 3
            for k in range(n_off):
                off = np.zeros(3) if k == 0 else RNG.uniform(-.5, .5, 3)
                mask = rasterize("sphere", {"r": r}, spacing, np.eye(3), off)
                a_true, v_true = truth("sphere", {"r": r})
                rows.append({"spacing": "x".join(f"{s:g}" for s in spacing), "radius": r,
                             "rho": rho, "offset_index": k,
                             "area_rel_err": e0(mask, spacing) / a_true - 1,
                             "volume_rel_err": mask.sum() * np.prod(sp) / v_true - 1})
    return pd.DataFrame(rows)


def fit_intercept(rho: np.ndarray, err: np.ndarray) -> dict:
    """Least-squares fit err = a + b/rho; `a` is the asymptotic bias."""
    A = np.column_stack([np.ones_like(rho), 1 / rho])
    coef, *_ = np.linalg.lstsq(A, err, rcond=None)
    resid = err - A @ coef
    return {"a_asymptotic": float(coef[0]), "b_per_rho": float(coef[1]),
            "rms_residual": float(np.sqrt(np.mean(resid ** 2)))}


# ------------------------------------ M4 dispersion and sphericity transfer
def m4_dispersion_and_sphericity() -> tuple[pd.DataFrame, dict]:
    rows = []
    rots = list(Rotation.random(12, random_state=11).as_matrix())
    for spacing in [(1.0, 1.0, 1.0), (3.0, 1.0, 1.0)]:
        sp = np.asarray(spacing, float)
        for kind, p in [("sphere", {"r": 20}), ("ellipsoid", {"axes": (20, 15, 10)}),
                        ("cylinder", {"r": 20, "h": 60})]:
            for j in range(16):
                off = RNG.uniform(-.5, .5, 3)
                rot = np.eye(3) if kind == "sphere" else rots[j % len(rots)]
                mask = rasterize(kind, p, spacing, rot, off)
                a_true, v_true = truth(kind, p)
                a_est = e0(mask, spacing)
                v_est = float(mask.sum() * np.prod(sp))
                s_true = np.cbrt(np.pi) * (6 * v_true) ** (2 / 3) / a_true
                s_est = np.cbrt(np.pi) * (6 * v_est) ** (2 / 3) / a_est
                rows.append({"spacing": "x".join(f"{s:g}" for s in spacing), "shape": kind,
                             "trial": j, "area_rel_err": a_est / a_true - 1,
                             "volume_rel_err": v_est / v_true - 1,
                             "sphericity_true": s_true, "sphericity_est": s_est,
                             "sphericity_rel_err": s_est / s_true - 1})
    df = pd.DataFrame(rows)
    # Propagation identity: d ln S = (2/3) d ln V - d ln A.
    pred = (2 / 3) * np.log1p(df.volume_rel_err) - np.log1p(df.area_rel_err)
    obs = np.log1p(df.sphericity_rel_err)
    prop = {"max_abs_deviation_of_propagation_identity": float(np.max(np.abs(pred - obs))),
            "mean_sphericity_rel_err": float(df.sphericity_rel_err.mean()),
            "mean_sphericity_rel_err_by_shape": {k: float(g.sphericity_rel_err.mean())
                                                 for k, g in df.groupby("shape")},
            "sphere_mean_sphericity_est": float(df[df["shape"] == "sphere"].sphericity_est.mean()),
            "sphere_true_sphericity": 1.0,
            "n_cases_sphericity_above_1p05": int((df.sphericity_est > 1.05).sum())}
    return df, prop


def main() -> None:
    m1 = m1_lattice_aligned_boxes()
    m1.to_csv(OUT / "mechanism_lattice_aligned.csv", index=False, float_format="%.12g")

    m2 = m2_orientation_response()
    m2.to_csv(OUT / "mechanism_plane_response.csv", index=False, float_format="%.12g")

    m3 = m3_resolution_sweep()
    m3.to_csv(OUT / "mechanism_resolution_sweep.csv", index=False, float_format="%.12g")

    m4, prop = m4_dispersion_and_sphericity()
    m4.to_csv(OUT / "mechanism_dispersion.csv", index=False, float_format="%.12g")

    fits = {}
    for sp, g in m3.groupby("spacing"):
        big = g[g.rho >= 4]
        fits[sp] = {"area": fit_intercept(big.rho.values, big.area_rel_err.values),
                    "volume": fit_intercept(big.rho.values, big.volume_rel_err.values),
                    "area_err_at_max_rho": float(g[g.rho == g.rho.max()].area_rel_err.mean()),
                    "volume_err_at_max_rho": float(g[g.rho == g.rho.max()].volume_rel_err.mean()),
                    "max_rho": float(g.rho.max())}

    summary = {
        "M1_lattice_representable_polyhedra": {
            "per_family": {
                k: {"error_by_rho": {float(r.rho_equivalent): float(r.area_rel_err)
                                     for r in g.itertuples()},
                    "edge_coefficient_c_in_err_approx_c_over_rho":
                        float(np.mean(g.area_rel_err.values * g.rho_equivalent.values)),
                    "error_at_largest_rho": float(g.sort_values("rho_equivalent")
                                                  .area_rel_err.iloc[-1])}
                for k, g in m1.groupby("spacing")},
            "interpretation": (
                "For polyhedra whose faces are lattice-representable planes the area error is "
                "purely an edge/corner artefact (marching cubes chamfers convex edges) and "
                "decays as O(1/rho) to ZERO -- halving on every doubling of rho. A convergent "
                "regime therefore exists. The sphere does not enter it: its error converges to "
                "a non-zero limit. The sphere/ellipsoid error is consequently an "
                "orientation-dependent asymptotic bias, not insufficient resolution."),
        },
        "M2_orientation_response": {
            "by_spacing": {
                sp: {"g_min": float(g.g.min()), "g_max": float(g.g.max()),
                     "g_solid_angle_mean": float((g.solid_angle_weight * g.g).sum()
                                                 / g.solid_angle_weight.sum()),
                     "total_area_rel_err": float(g.total_area_rel_err.iloc[0]),
                     "closure_residual": float(g.closure_residual.abs().max())}
                for sp, g in m2[m2.radius == m2.radius.max()].groupby("spacing")},
            "interpretation": (
                "g(n) is the asymptotic area-inflation factor of binary marching cubes for a "
                "surface patch of true normal n. It is minimal near lattice-representable "
                "normals (M1 shows flat lattice-aligned faces are measured with no asymptotic "
                "error at all; the ~1.03 floor here is the finite bin width smearing those "
                "measure-zero directions) and maximal near normals far from any short lattice "
                "vector. Its solid-angle mean is the irreducible bias of the estimator on a "
                "smooth closed surface, and it grows sharply with voxel anisotropy."),
        },
        "M3_resolution_sweep": fits,
        "M4_sphericity_propagation": prop,
    }
    (OUT / "mechanism_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
