# Why the production surface estimator fails, mechanistically

Estimator under analysis: `E0_marching_cubes_binary_v1`
(`skimage.measure.marching_cubes`, level 0.5, one-voxel zero pad, native
anisotropic spacing, no smoothing), reached through the production call path
`organoid_analysis.quantification.features.geometry`.

Baseline re-measurement reproduces the frozen dev-grid record to
`4.5e-10` relative (tolerance `1e-9`), so everything below describes the
estimator the pipeline actually ships.

---

## Summary of the finding

**The error is not a resolution deficit. It is an orientation-dependent
asymptotic bias of the estimator, and it does not go away at any voxel size.**

Binary marching cubes reads only the *sign* of `f - 0.5` at each voxel. Every
reconstructed vertex is therefore pinned to the exact midpoint of a lattice
edge, and the triangulation in each cell is drawn from a fixed finite catalogue
of configurations. For a smooth surface the resulting mesh converges — as the
voxel size goes to zero — not to the surface but to a *lattice-faceted
approximation of it* whose area is larger by an orientation-dependent factor
`g(n) >= 1`.

Three independent experiments establish this.

---

## M1 — a convergent regime exists, and the sphere is not in it

Digitized axis-aligned cubes, and a cube rotated 45° about one axis (so its
slanted faces are `(011)` lattice planes), were measured at
`rho = half-extent / max(spacing)` from 1 to 32.

| rho | cube, 1×1×1 | cube, 3×1×1 | cube, 2×0.7×0.7 | 45° prism |
|----:|------------:|------------:|----------------:|----------:|
| 1   | −32.43 %    | −29.48 %    | −29.66 %        | —         |
| 2   | −15.43 %    | −13.50 %    | −13.62 %        | —         |
| 4   | −7.52 %     | −6.44 %     | −6.51 %         | +5.21 %   |
| 8   | −3.71 %     | −3.14 %     | −3.18 %         | +2.63 %   |
| 16  | −1.84 %     | −1.55 %     | −1.57 %         | +1.32 %   |
| 32  | −0.918 %    | −0.771 %    | −0.781 %        | +0.662 %  |

The error halves on every doubling of `rho`: it is `c/rho` with
`c ≈ −0.30` (isotropic cube), `−0.26` (anisotropic cubes), `+0.21` (45° prism),
and it extrapolates to **zero**. Its origin is purely the *edges*: marching
cubes chamfers a convex right-angle edge, so the deficit scales as
(edge length × voxel size) / area. The flat faces themselves carry no
asymptotic error.

So the estimator does possess a convergent regime. The sphere never enters it:
its error rises with `rho` and settles at a positive constant (figure, panel a).
Whatever is wrong with the sphere is therefore *not* insufficient sampling.

**Consequence for this project.** The `c/rho` edge term is not an artefact of
contrived test shapes — it is exactly what a cylinder's circular rim, or any
real organoid boundary with a crease, contributes. It is the one error term that
*does* respond to resolution, and no estimator that reads a binary mask can
remove it, because the mask does not record where inside the voxel the edge lay.

---

## M2 — the orientation response `g(n)`

A digitized sphere presents every surface orientation simultaneously. Each
marching-cubes triangle of an `r/s = 96` sphere was assigned to the solid-angle
bin of the *true* normal at its centroid (the radial direction), giving the
local area-inflation factor `g(n)` directly. The solid-angle-weighted mean of
`g` reproduces the total sphere area error to `< 1e-7`, which closes the
accounting: nothing else contributes.

| voxel spacing | `g` min | `g` max | solid-angle mean | total sphere area error |
|---|---:|---:|---:|---:|
| 1 × 1 × 1     | 1.0299 | 1.1261 | **1.0875** | +8.75 % |
| 2 × 0.7 × 0.7 | 1.0251 | 1.2808 | **1.1583** | +15.83 % |
| 3 × 1 × 1     | 1.0233 | 1.2869 | **1.1622** | +16.22 % |

`g` is smallest near lattice-representable normals — M1 shows those are measured
with no asymptotic error at all; the ≈1.03 floor here is finite bin width
smearing a measure-zero set of directions — and largest for normals far from any
short lattice vector. Voxel anisotropy makes this dramatically worse: at 3:1:1
the worst orientation is inflated by 29 %, and even the *average* over
orientations is +16 %.

This is the dominant term, and it is a property of the estimator, not of the
object.

---

## M3 — resolution sweep: area plateaus, volume converges

Spheres from `rho = 1.5` to `rho = 128`, multiple subvoxel offsets, fitted as
`err(rho) = a + b/rho` over `rho >= 4`:

| spacing | area `a` (asymptote) | area at max `rho` | volume `a` | volume at max `rho` |
|---|---:|---:|---:|---:|
| 1 × 1 × 1     | +9.11 %  | +8.78 % (ρ=128) | +0.180 % | −0.0058 % |
| 2 × 0.7 × 0.7 | +16.14 % | +15.88 % (ρ=64) | +0.134 % | −0.0008 % |
| 3 × 1 × 1     | +16.61 % | +16.35 % (ρ=42.7) | +0.098 % | −0.0005 % |

**Volume and surface area behave completely differently.** Voxel counting is a
consistent volume estimator — its error decays to `5.8e-5` at `rho = 128` and it
meets the §9 `<1 %` volume criterion across the sweep. Marching-cubes area on
the same masks converges to `+8.8 %` (isotropic) and `+16 %` (3:1:1).

This is why the project can hold a defensible volume claim and simultaneously
fail on area, and it is why "increase the resolution" is not a remedy: at
`rho = 128` — far beyond any achievable organoid imaging resolution — the area
error is still 8.8 %, i.e. still 1.8× the acceptance limit.

---

## M4 — how the bias propagates into sphericity

Sphericity `S = pi^(1/3) (6V)^(2/3) / A` propagates the two errors as

```
d ln S = (2/3) d ln V  -  d ln A
```

which the measurements satisfy to `3.3e-16` (exact, as it must be — it is an
identity, not a model). Because `d ln V ≈ 0` and `d ln A ≈ +0.09…+0.16`,
sphericity inherits essentially the whole area bias with a negative sign:

| shape | mean sphericity error |
|---|---:|
| sphere   | −10.98 % |
| ellipsoid| −10.52 % |
| cylinder | −9.96 % |

A digitized sphere reports sphericity **0.890** instead of 1.000.

Two operational consequences:

1. Any cross-study or literature comparison of sphericity computed this way is
   biased low by ~10 % (isotropic) and more under anisotropy — and the bias
   depends on voxel spacing, so it is not even a constant offset between
   datasets acquired on different microscopes.
2. The existing `SPHERICITY_REVIEW_LIMIT = 1.05` guard **cannot detect this
   failure**: the bias pushes sphericity *down*, and zero of the 96 dispersion
   cases exceeded 1.05. The guard is designed for the opposite sign of error.

---

## What a candidate estimator must therefore fix

The total error decomposes as

```
relative area error  ≈  [ mean over the surface of ( g(n) - 1 ) ]   +   c / rho
                         \__________ orientation bias __________/       \_edge_/
                          asymptotic; NOT removed by resolution       decays as 1/rho
```

* The **orientation term** is the target. It is removable in principle, because
  it comes from discarding the subvoxel position of the boundary, and any
  estimator that either (i) retains a continuous field or (ii) integrates over
  many lattice directions with a correct quadrature weight can suppress it.
* The **edge term** is *not* removable from a binary mask, and it sets a hard
  floor of roughly `0.2–0.3 / rho` for objects with creases or flat-face
  junctions. At `rho = 6` that is already 3–5 % on its own. This is the honest
  reason why a per-case `<5 %` guarantee across cylinder-like shapes at low
  `rho` may be unreachable, and it is declared in advance rather than
  discovered after the confirmation run.

Files: `mechanism_lattice_aligned.csv`, `mechanism_plane_response.csv`,
`mechanism_resolution_sweep.csv`, `mechanism_dispersion.csv`,
`mechanism_summary.json`, `fig_mechanism.png`.
