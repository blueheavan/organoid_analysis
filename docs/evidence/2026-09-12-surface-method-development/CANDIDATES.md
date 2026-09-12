# Candidate estimators and where their constants come from

The brief requires that any new weighting coefficient have a stated
mathematical basis and provenance, and that no parameter be fitted to the
phantoms. This document states, for every number that enters a candidate,
where it comes from. The implementation is `estimators_v2.py`.

Priority order follows the brief: continuous-field / subvoxel first,
digital-geometry / local-configuration second. Mesh smoothing is **not**
implemented, as instructed — it changes the object rather than the
measurement, and its smoothing length would be a free parameter with no
independent basis.

---

## Candidate A — continuous-field marching cubes (subvoxel)

**Basis.** M2 showed the binary estimator's asymptotic bias arises from pinning
every mesh vertex to a lattice-edge midpoint, which is the only information a
sign test leaves. Given a continuous field `phi` whose iso-level `L` is the
boundary, linear interpolation along each lattice edge recovers the true
crossing point to `O(h^2)`, the facet normals converge to the true normals, and
the orientation term disappears to leading order.

**Parameters.** Two, both fixed by the definition of the field, neither tuned:

| variant | field | iso-level | padding |
|---|---|---|---|
| `A1_sdf` | signed distance, inside negative | `0` | 1 voxel at `max(phi)+max(spacing)` |
| `A2_occ` | occupancy (partial volume) in `[0,1]` | `0.5` | 1 voxel at `0` |

The iso-levels are not choices: `0` is the definition of a zero level set, and
`0.5` is the occupancy at which a voxel is half inside, which is where a
straight boundary crosses the voxel centre. No smoothing kernel, no sigma.

**Provenance constraint (important).** The production Cellpose wrapper consumes
only `cellprob_threshold` and keeps the thresholded mask; the continuous
probability field is discarded. So A1/A2 **cannot be evaluated on production
data as the pipeline currently stands.** They are carried through the full V&V
as *reference / conditional* candidates: they bound what is achievable, and if
they pass they identify the pipeline change worth making. They are not eligible
for unconditional qualification in this round, and the frozen selection rule
records that explicitly rather than discovering it later.

---

## Candidate B — generalised Cauchy–Crofton over an extended lattice stencil

### B.1 The estimator

For a primitive integer direction `k` (gcd of components = 1), with physical
vector `v_k = k * spacing`, unit direction `u_k`, and `V` the voxel volume,
the discrete lines in direction `k` have one line per cross-sectional area
`a_k = V / |v_k|`. Counting `N_k`, the number of adjacent-along-`k` voxel pairs
of differing label, integral geometry gives

```
N_k * a_k  ~=  integral over S of |n . u_k| dA                      (3)
```

so for any non-negative weights `w_k`,

```
Ahat = sum_k w_k N_k a_k = integral over S of R(n) dA,
       R(n) := sum_k w_k |n . u_k|                                  (5)
```

The estimator is exact for every surface at every orientation **exactly to the
extent that the plane response `R(n)` equals 1**. `R(n) - 1` is the orientation
bias, in the same units as the `g(n) - 1` measured in M2. This is what turns
the weight choice into a stated approximation problem instead of a knob.

### B.2 Where the weights come from

**`Bvor` — spherical-Voronoi quadrature.** Read (5) as a quadrature rule for
`integral over S^2 of |n.u| dsigma = 2*pi`. Give each direction the solid angle
`Omega_k` of its spherical Voronoi cell (summing both hemispheres) and set
`w_k = Omega_k / (2*pi)`; then `R -> 1` as the direction set densifies.
Provenance: classical integral geometry. **At `m = 1` (13 directions) this is
exactly the estimator the project has already tested and rejected**; it is
reproduced here as a labelled reference, not offered as new.

**`Blp` — minimax weights.** `R(n)` is linear in `w`, so

```
minimise  t   over w >= 0, t >= 0
subject to  |R(n) - 1| <= t  for all unit n
```

is a semi-infinite linear program, solved here by cutting planes: solve on a
coarse normal set, find the worst violators on a 20011-point set, add them,
repeat until the LP optimum `t*` and the dense-set supremum agree to `1e-5`.
The reported bound is therefore verified, not assumed. Provenance: Chebyshev
approximation; the only inputs are the lattice, the spacing and the choice to
require non-negativity (which keeps the estimator monotone in the counts).

Nothing in either scheme sees a phantom.

### B.3 Second error term, derived in closed form

Equation (3) is not exact: lattice points along a discrete line are spaced
`L_k = |v_k|` apart, so a chord shorter than `L_k` may contain no sample and
contribute 0 transitions instead of 2. Every convex body has such grazing
chords near its silhouette, so the count is biased **low**, and the bias
**grows** with stencil radius — opposite to the orientation bias.

For period `L` the probability of at least one sample inside a chord of length
`c` is `min(1, c/L)`, so the expected deficit is `2*max(0, 1 - c/L)`. For a
sphere of radius `R`, substituting `t = sqrt(R^2 - b^2)` in the impact-parameter
integral:

```
missed = int_0^{L/2} 2 (1 - 2t/L) 2 pi t dt = 4 pi (L^2/8 - L^2/12) = pi L^2 / 6
```

against `int_S |n.u_k| dA = 2 pi R^2`, a relative deficit of `L_k^2/(12 R^2)`.
Weighting by `w_k` and using `sum_k w_k = 2` (which follows from `R == 1`
integrated over the sphere):

```
relative deficit = sum_k w_k L_k^2 / (24 R^2)  =:  D / R^2           (6)
```

`D` depends only on spacing, stencil and weights — known before any
measurement. **Verified against measurement** (spheres, Voronoi weights,
`m = 1..4`): predicted vs observed deficit agrees to `0.01–0.05 %` at
`rho >= 12`, `~0.4 %` at `rho = 6`, degrading only in the extreme
`rho = 2` corner. The formula was derived first and then checked, not fitted.

### B.4 The resulting a-priori error budget

```
|relative area error|  <~  t*(spacing, m)  +  D(spacing, m) / R_eff^2   (7)
```

with `R_eff = (3V / 4pi)^(1/3)` the volume-equivalent radius, obtained from the
object's own voxel count. The two terms move oppositely in `m`, so (7) has an
interior optimum, and **the stencil radius is chosen by minimising it** — a
computation on lattice constants and one number read off the mask. No phantom,
no truth value, no fitted constant. This is the frozen stencil rule.

Two residuals are *not* covered by (7) and are declared as such:

* the crease/edge term `c/rho` of M1 (`c ~ 0.2–0.3`), which no binary-mask
  estimator can remove, because the mask does not record where inside the voxel
  the crease lay;
* integer quantisation of `N_k`.

### B.5 Note on `t*` as a worst case

`t*` is the supremum of `|R(n) - 1|` over *single* orientations, so it is the
error of a flat plate. A closed surface integrates `R(n) - 1` against its own
normal distribution, which for a sphere is the uniform distribution and gives a
much smaller number — this is why a 13-direction rule with `t* = 15 %` at 3:1:1
nonetheless measures a sphere to `0.06 %`. Both quantities are reported: `t*`
as the rigorous worst case over shapes, and the RMS deviation of `R(n) - 1` as
the quantity relevant to closed, orientation-averaging bodies. The applicability
domain is stated in terms of shapes, not only resolution, for exactly this
reason: a plate-like object with all normals in a narrow solid angle is *not*
covered by the sphere/ellipsoid evidence.
