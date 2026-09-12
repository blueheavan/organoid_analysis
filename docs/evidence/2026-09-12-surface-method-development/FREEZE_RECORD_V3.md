# Round-2 freeze record — `crofton_minimax_sym_v3`

Written and hashed **before** the `confirm2` grid was executed. Everything below
is fixed. No parameter, threshold, criterion, domain boundary or acceptance rule
may change after this point; the confirmation grid is executed exactly once and
its outcome is reported whatever it is.

## 1. Method

| item | value |
|---|---|
| method name | `crofton_minimax_sym_v3` |
| weights | orbit-symmetry-reduced minimax LP over the metric symmetry group of the voxel lattice |
| stencil radii offered | `m ∈ {1,2,3,4,5}` |
| stencil selection | `argmin_m apriori_bound_v3(spacing, m, r_in)` — computable from mask + spacing alone |
| area | `Σ_k w_k · transitions_k(mask)` |
| volume | voxel count × voxel volume (unchanged from legacy) |
| sphericity | `π^(1/3)(6V)^(2/3) / A` |

`m = 6` is **excluded**. It solved for 10 of 11 spacing ratios (≤ 250 s) but did
not converge for spacing `1.2×1×1` in 900 s / 80 cutting-plane iterations.
Production meets spacings that are not in this table and must solve at
measurement time, so a radius that can hang on an untested ratio cannot be
offered for any. Evidence: `m6_exclusion.json`, `lp_symmetry_verification.csv`.
The cap was **not** chosen to save compute and **not** lowered after seeing any
phantom result.

## 2. What changed from round 1, and why

Round 1 failed confirmation because its domain variable was the
volume-equivalent radius `rho_eff`. That overstates resolution for elongated or
creased objects: the three failing confirmation cylinders had `rho_eff` 8.2–10.5,
comfortably inside a declared `rho_eff ≥ 8`, while their inscribed radii were
6.6 — genuinely under-resolved. The domain variable is now the **inscribed
radius** `rho_in = max(distance_transform_edt(mask, sampling)) / max(spacing)`,
which is a local feature-size measure and is computable from the mask alone.
On round-1 data it places those three failures at 6.57–6.63, correctly below any
sensible threshold (`rho_in_scoping.csv`).

## 3. Declared applicability domain

**Computable gate** (applied by `in_domain()`):

    rho_in ≥ 10.0        and        anisotropy = max(spacing)/min(spacing) ≤ 4.0

**Scope of intended use** (not machine-checkable):

    smooth closed surfaces — surfaces without dihedral creases

## 4. Why the scope statement exists, and why it is not a gate

Development introduced a maximally creased class (box, including the exact
cube) specifically to test whether creases break the estimator. They do, and the
effect is a distinct error mode rather than a resolution deficiency:

| class | aspect | smooth? | worst \|area err\|, rho_in ≥ 6 |
|---|---|---|---|
| sphere | 1:1 | yes | 1.33 % |
| ellipsoid | 2:1.5:1 | yes | 0.73 % |
| capsule | 3:1 | yes | 1.12 % |
| torus (genus 1) | R/r = 2 | yes | 0.74 % |
| cylinder | 3:1 | **no** (2 creases) | 4.96 % |
| box / cube | 1:1:1 … 2:1.5:1 | **no** (12 creases) | 7.60 % |

The discriminator is creases, not elongation and not topology: the capsule and
the cylinder have the same 3:1 aspect and differ by two circular creases, and
their errors differ by a factor of four. The bias is negative, independent of
rotation (−4.17 / −3.33 / −3.59 % for the three box orientations), and exceeds
the estimator's own a-priori budget in 89 % of box cases while never exceeding
it in 150 smooth cases.

It decays far too slowly to be removed by resolution. Isotropic cube, worst
\|area err\| by `rho_in` band: 10.98 % (≤8), 7.60 % (8–11), 6.28 % (11–14),
5.00 % (14–18), 4.62 % (18–22), 4.22 % (22–28), 3.54 % (28–40). A universal
threshold covering creased surfaces would sit near `rho_in ≈ 22` with under one
percentage point of margin, and the predeclared `confirm2` ladder stops at 15,
so it could not be tested. **Creased surfaces are therefore declared NOT
QUALIFIED at any resolution in this round.** They are not removed from the
grids: they are executed, reported, and counted in `confirm2` as out-of-scope
cases with their numbers published.

A machine-checkable crease gate was attempted and **rejected**: the
marching-cubes / Crofton disagreement, the only truth-free indicator available
from a single mask, overlaps completely between classes (box 1.3–23.0 %, sphere
8.8–17.1 %) and correlates with the error at r = −0.20
(`crease_indicator_rejected.md`). Since no gate exists, the scope statement is
the caller's responsibility, and it is recorded as a limitation rather than
disguised as a guarantee.

## 5. Basis for the declared thresholds

From `vv_dev2_all.csv`, smooth classes only, 1083 development cases total:

| gate | n | worst area | worst volume | worst sphericity | violations |
|---|---|---|---|---|---|
| rho_in ≥ 6, aniso ≤ 4 | 361 | 1.33 % | 1.25 % | 1.15 % | 3 volume |
| rho_in ≥ 8, aniso ≤ 4 | 247 | 1.33 % | 1.10 % | 1.15 % | 2 volume |
| **rho_in ≥ 10, aniso ≤ 4** | **136** | **1.33 %** | **0.556 %** | **0.97 %** | **none** |

`rho_in ≥ 10` is set by the **volume** criterion, not the area criterion: voxel-count
volume quantisation is the binding term, and 10 is the first ladder value at
which it clears 1 % across every smooth class and every anisotropy tested.
Margins carried into confirmation: area 3.8×, volume 1.8×, sphericity 5.2×.
`aniso ≤ 4` is the largest anisotropy present in development (`4×1×1`); the
confirmation grid's 3.5 interpolates between development's 3.125 and 4, and its
anisotropy-5 spacing falls outside the domain by declaration.

## 6. Acceptance rules — fixed now, applied once

Let **S** = `confirm2` cases that are in-scope (smooth classes) and pass the
computable gate. The candidate is **QUALIFIED** if and only if all three hold on
**every** case in S — worst-case, not mean:

| # | rule |
|---|---|
| A1 | `abs(area_rel_err) < 5 %` for every case in S |
| A2 | `abs(volume_rel_err) < 1 %` for every case in S |
| A3 | `abs(sphericity_rel_err) < 5 %` for every case in S |

and the run is valid only if:

| # | rule |
|---|---|
| V1 | S is non-empty and contains all four smooth classes |
| V2 | every case in S produced a measurement (no solver failure, no exception) |
| V3 | the frozen file hashes below match at execution time |

Any violation of A1–A3 ⇒ **FAIL / NOT QUALIFIED**. Failure of V1–V3 ⇒ the run is
void and reported as such; it is not retried with altered inputs.

Reported but **not** part of the verdict: all out-of-scope cases (creased
classes, and in-scope cases below the gate), and the legacy estimator `E0`
evaluated on the identical masks.

## 7. Pre-registered predictions

Stated before execution so the confirmation set can falsify them:

1. All A1–A3 pass on S, with worst area error below 2.5 % (development worst on
   the same gate was 1.33 %; confirmation probes new ladder positions, aspects
   and spacings, so some degradation is expected).
2. Creased out-of-scope cases will violate A1 at `rho_in` near 10–15, in the
   3–7 % band.
3. `E0` will violate A1 on essentially every in-scope case (development: 361 of
   361, median 13.1 %).

If prediction 1 fails, the verdict is FAIL and the round ends NOT QUALIFIED.

## 8. Deviations from the round-2 protocol, recorded before execution

1. **Pre-declared check restated.** `ROUND2_PROTOCOL.md` required the reduced LP
   to reproduce round-1 weights to solver tolerance. The LP optimum is not
   unique, so the weight vectors need not match and the literal check is not
   meaningful. It was applied instead to the verified plane-response deviation
   `sup|R−1|`, which is the quantity that enters the error budget: v3 agrees
   within 0.94 % everywhere and is strictly better — by up to 17.8 % — on the
   four `m=4` cases where round 1's cutting-plane loop had not converged. Round
   1's reported bounds at those radii were optimistic relative to its own
   achieved response.
2. **Development grid extended after first analysis.** `dev2b` (higher `rho`
   ladder for the two creased classes) and `dev2c` (exact cube) were added after
   `dev2` was analysed, because `dev2`'s ladder stopped at 16 and the crease
   behaviour was still changing there. Declaring a threshold on that evidence
   would have repeated round 1's error of putting the boundary where coverage
   runs out. These are development data; `confirm2` was not touched, not
   executed, and not consulted.
3. **`m ≤ 5` cap** as in §1.

## 9. Frozen artefacts

Hashes in `freeze_record_v3.json`, written by `step9_freeze_v3.py` immediately
before the confirmation run:
`common2.py`, `estimators_v3.py`, `estimators_v2.py`, `common.py`,
`surface_area_v3.py`, `step8_dev2_run.py`, `step10_confirm2_run.py`,
and every weight file under `out/crofton_weights_v3/`.
