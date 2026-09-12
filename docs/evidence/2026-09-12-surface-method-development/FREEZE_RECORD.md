# Freeze record — surface-area estimator qualification

Written **before** the confirmation grid was rasterized, measured or scored.
`freeze_record.json` carries the SHA-256 of every file named here; that hash
list is what makes this document checkable rather than merely asserted.

Everything in this file is fixed. Nothing below may be revised in the light of
confirmation results. If the confirmation run fails, the outcome is
`FAIL / NOT QUALIFIED` and the remedy is another development cycle against a
*new* confirmation set — not an edit to this record.

---

## 1. Frozen estimator

**`crofton_minimax_adaptive_v2`** — `surface_area_v2.py`, entry point
`measure(mask, spacing)`.

| element | value | where it comes from |
|---|---|---|
| family | generalised Cauchy–Crofton on a primitive-direction stencil | integral geometry; `CANDIDATES.md` §B.1 |
| weights | minimax (Chebyshev) solution of `min_w max_n abs(sum_k w_k abs(n·u_k) − 1)`, `w ≥ 0` | semi-infinite LP by cutting planes, verified against a 20011-normal set to `1e-5`; no phantom input |
| stencil radius `m` | `argmin` over `(1,2,3,4)` of `t*(spacing,m) + D(spacing,m)/R_eff²` | a-priori budget, eq. (7); inputs are the voxel spacing and the object's own voxel count |
| `R_eff` | `(3V/4π)^(1/3)`, `V` = voxel count × voxel volume | measurable from the mask alone |
| input | binary label mask + spacing — **exactly what the production pipeline already retains** | no pipeline change required |

No coefficient in this estimator was chosen by looking at a phantom error.
The weights come from a lattice-geometry optimisation; `m` comes from
minimising a budget whose two terms were derived analytically (orientation
bias from the plane response, grazing-chord deficit from an impact-parameter
integral) and whose grazing term was then *verified* against measurement
rather than fitted to it.

**Declared limitation, not a result.** `M_CANDIDATES` stops at 4 because the
minimax solve at `m ≥ 5` (`K ≥ 577` directions) did not converge in practical
time in the dense cutting-plane formulation used here. This is a compute
limitation of this round. It propagates directly into the anisotropy limit in
§3, because raising `m` is the only mechanism the estimator has for
compensating anisotropy. The route to lifting it is orbit-symmetry reduction
of the LP (for spacings of the form `(a,b,b)` the symmetry group has order 16,
cutting 577 variables to roughly 75), which is a solver change, not a new
parameter. It is named here so that the limit is attributable.

## 2. Comparators carried through, and their standing

| tag | what | standing |
|---|---|---|
| `E0` | `marching_cubes_binary_lewiner_v1`, the current production estimator | **preserved, not replaced.** Its identity, parameters and failing evidence are recorded in `E0_legacy_baseline.json`. It is measured on every case of every grid alongside the candidate. |
| `A1_sdf`, `A2_occ` | continuous-field marching cubes | **conditionally eligible only.** The Cellpose wrapper consumes `cellprob_threshold` and discards the continuous field, so these cannot be computed on production data today. Reported, never qualified in this round. |
| `Bvor_m1` | 13-direction spherical-Voronoi Crofton | **the estimator this project already tested and rejected.** Carried as a labelled reference. It is not the frozen candidate and is not presented as new. |

## 3. Qualified applicability domain

A measurement is in-domain iff **both**:

* **`rho_eff = R_eff / max(spacing) ≥ 8`**
* **`anisotropy = max(spacing) / min(spacing) ≤ 3.0`**

No shape-class restriction. Cylinders — which carry the crease term that no
binary-mask estimator can remove — are **inside** the domain and count towards
the verdict. Excluding them would have narrowed the intended use to shapes
organoids do not reliably have, and would have made the claim easier in a way
the evidence does not support.

**Basis for `rho_eff ≥ 8`.** Both residual error terms scale as inverse powers
of `rho_eff`: the derived grazing deficit as `D/R²`, and the crease term
measured in M1 as `c/rho` with `c ≈ 0.2–0.3`. Development evidence: the
candidate's worst in-domain case is `3.49 %` at `rho_eff ≥ 8`, against
`5.14 %` at `rho_eff ≥ 6` and `5.52 %` at `rho_eff ≥ 4` — i.e. the criterion
is already violated on development at the next threshold down. `8` is the
weakest threshold at which development shows margin, not the threshold that
makes the number look best.

**Basis for `anisotropy ≤ 3.0`.** This is the range development spans
(`1.0, 2.0, 2.857, 3.0`), and the mechanism that would carry the estimator
beyond it — increasing `m` to suppress `t*` — is exhausted at the frozen cap
`m = 4` (§1). Development shows the compensation working across the whole
spanned range with no upward trend (`3.49 %` isotropic, `2.67–3.03 %` at
anisotropy 2–3), so the limit is set by where the mechanism runs out, not by
where the errors start growing. Extrapolating past it is not supported.

**Cost of this declaration, stated in advance.** The confirmation grid uses
spacings `1.5×1×1` (aniso 1.5), `2.5×0.8×0.8` (3.125) and `4×1×1` (4.0). Only
the first is inside the declared anisotropy limit, so a substantial part of
the confirmation set will fall out of domain. That weakens the confirmation
and it is the honest consequence of the `m ≤ 4` cap. All out-of-domain cases
are measured and reported anyway; they cannot rescue a failing verdict and
they do not contribute to a passing one.

## 4. Acceptance rules

For **every** case inside the domain — no averaging, no trimming, no
worst-case exclusion:

| id | quantity | criterion |
|---|---|---|
| A1 | `abs(area_est/area_true − 1)` | `< 5 %` |
| A2 | `abs(volume_est/volume_true − 1)` | `< 1 %` |
| A3 | `abs(sphericity_est/sphericity_true − 1)` | `< 5 %` |

A3 is present because the brief requires sphericity to be re-validated after
any surface change; sphericity inherits the area error directly
(`S ∝ V^(2/3)/A`).

**Verdict rule.** `PASS` iff A1, A2 and A3 hold for every in-domain case on the
confirmation grid. Any single violation ⇒ `FAIL / NOT QUALIFIED`. The harness
`step6_vv_run.py` computes this and has no branch that can soften it.

## 5. Selection predicate (applied before the confirmation run)

Among candidates, in this order:

1. **Eligibility.** Computable from what the production pipeline retains today
   (binary mask + spacing). `A1_sdf`/`A2_occ` fail this and are conditional.
2. **Parameters not phantom-selected.** A fixed-`m` variant is eligible only
   with an independent justification for that `m`. `Blp_m2` has the best
   development number of any fixed variant (`2.747 %`) and is **rejected on
   this ground** — `m = 2` would have been chosen because it won on the
   development phantoms, which is exactly the post-hoc parameter fitting the
   brief forbids. The adaptive rule selects `m` from spacing and object size
   with no phantom input and is therefore eligible despite a worse number.
3. **Development screen.** A1, A2, A3 must hold for every in-domain dev case.
4. **Ranking.** Smallest in-domain dev max `abs(area error)`; ties to fewer
   directions.

Applied to development (156 in-domain rows over 78 cases): `Blp_auto` = `3.49 %`,
`Bvor_auto` = `4.43 %`, `Blp_m1` = `7.95 %` (fails), `Bvor_m1` = `11.35 %`
(fails), `E0` = `18.74 %` (fails). Selected: **`Blp_auto` =
`crofton_minimax_adaptive_v2`**.

If the confirmation run fails, the reported outcome is `FAIL / NOT QUALIFIED`
for this estimator on this domain. No fallback candidate is promoted in its
place, because no other candidate was frozen.

## 6. Confirmation protocol

* Grid `confirm` — seed `424242`, defined in `common.py` since the previous
  frozen study and **never executed until now**.
* Run once: `python step6_vv_run.py --grid confirm`.
* Outputs `vv_confirm.csv` and `vv_confirm_summary.json`.
* No re-runs with altered parameters. No case deletions. If a case errors, the
  error is reported as a failure of the method to produce a measurement.
