# Surface-area measurement: round-2 method development and V&V

## Verdict

**QUALIFIED**, within the declared domain, on an untouched confirmation set.

    method    crofton_minimax_sym_v3
    domain    smooth closed surfaces, rho_in >= 10, anisotropy <= 4
    result    96/96 confirmation cases pass all three criteria
              worst surface-area error 0.951 %  (criterion 5 %)
              worst volume error        0.761 %  (criterion 1 %)
              worst sphericity error    0.949 %  (criterion 5 %)

The requirement in the brief — *surface-area absolute relative error < 5 % for
every case within the predeclared qualified domain on the untouched
confirmation set* — is met with a factor of 5.3 in hand. The confirmation grid
was executed once, under hashed code, after the domain and the acceptance rules
were frozen.

One part of the problem is **not** solved, and is declared so rather than
hidden: surfaces with dihedral creases are NOT QUALIFIED at any resolution.

---

## 1. What was wrong, and what fixed it

The production estimator `E0` is marching cubes on the binary mask. On the
confirmation set it errs by a median of 12.5 % and a worst of 21.9 % on exactly
the cases the new estimator handles to 0.33 % median. `E0` violates the 5 %
criterion on **96 of 96** in-domain cases. Its error does not shrink with
resolution — this is the defect that motivated the work, and it is a bias, not
noise: marching cubes on a binarised field measures the area of a staircase.

The replacement is an integral-geometric estimator: surface area from weighted
counts of mask transitions along lattice directions (Crofton / Cauchy), with
the direction weights obtained from a minimax linear program that minimises the
worst-case deviation of the estimator's response over all plane orientations.
The weights therefore carry a stated mathematical basis and an a-priori error
budget that is computable per measurement, before any comparison with truth. On
all 96 confirmation cases the realised error stayed inside that budget. See the
addendum at the end of this report: integration found a placement that exceeds
it, so the budget is a validated predictive budget rather than a proven bound.

Two changes made round 2 work where round 1 failed.

**Orbit-symmetry reduction.** Round 1 capped the stencil radius at `m <= 4`
because the LP would not solve at larger radii — 577 free variables at `m = 5`.
Constraining the weights to be constant on the orbits of the voxel lattice's
metric symmetry group reduces this to 98 variables and the solve from a hang to
under a second. The cap is now `m <= 5`, set by a genuine non-convergence
(§4), not by compute.

**The domain variable.** This is the substantive fix. Round 1 declared its
domain in the volume-equivalent radius `rho_eff`, which overstates resolution
for elongated or creased objects. Its three failing confirmation cylinders had
`rho_eff` 8.2–10.5, inside the declared `rho_eff >= 8`, but inscribed radii of
6.6 — genuinely under-resolved. The domain is now declared in the **inscribed
radius** `rho_in = max(EDT(mask, sampling)) / max(spacing)`, a local
feature-size measure computable from the mask alone. It places those three
round-1 failures at 6.57–6.63, correctly outside any sensible threshold.

---

## 2. The declared domain, and why it has two parts

**Computable gate**, enforced by `in_domain()`:

    rho_in >= 10        anisotropy = max(spacing)/min(spacing) <= 4

`rho_in >= 10` is set by the **volume** criterion, not the area criterion.
Surface area is comfortable from `rho_in` ~ 5; voxel-count volume quantisation
is the binding term, and 10 is the first ladder value clearing 1 % across every
smooth class and every anisotropy tested.

**Scope of intended use**, not machine-checkable:

    smooth closed surfaces (no dihedral creases)

This second part exists because development found a distinct error mode that no
resolution threshold removes. It is a limitation of the method, recorded as
such.

---

## 3. The crease limitation

A maximally creased class (box, including the exact cube) was added in round 2
specifically to test whether creases break the estimator. They do:

| class | aspect | smooth | worst \|area err\|, rho_in >= 6 |
|---|---|---|---|
| sphere | 1:1 | yes | 1.33 % |
| ellipsoid | 2:1.5:1 | yes | 0.73 % |
| capsule | 3:1 | yes | 1.12 % |
| torus (genus 1) | R/r = 2 | yes | 0.74 % |
| cylinder | 3:1 | **no** | 4.96 % |
| box / cube | 1:1:1 … 2:1.5:1 | **no** | 7.60 % |

The discriminator is creases, not elongation and not topology. The capsule and
the cylinder share a 3:1 aspect and differ only by two circular creases; their
errors differ by a factor of four. The torus is topologically non-trivial and
behaves like the sphere. The bias is negative, rotation-independent (−4.17 /
−3.33 / −3.59 % across three box orientations, so it is not digitisation of
tilted faces), and exceeds the estimator's own a-priori budget in 89 % of box
cases while never exceeding it in 150 smooth cases.

It decays far too slowly to be gated away. Isotropic cube, worst \|area err\| by
`rho_in` band: 10.98 % (≤8), 7.60 % (8–11), 6.28 % (11–14), 5.00 % (14–18),
4.62 % (18–22), 4.22 % (22–28), 3.54 % (28–40). A universal threshold covering
creased surfaces would sit near `rho_in ≈ 22` with under one percentage point
of margin — and the predeclared confirmation ladder stops at 15, so it could
not have been tested this round. Creased surfaces are therefore declared not
qualified rather than qualified on untested extrapolation.

They were **not** removed from the grids. All 160 creased confirmation cases
were executed and are reported: worst 13.77 %, and 8 of the 49 that pass the
computable gate violate the 5 % criterion. That is the declared limitation
behaving exactly as declared.

A mask-computable crease gate was attempted and **rejected** before the
confirmation run — the marching-cubes/Crofton disagreement, the only truth-free
indicator available from a single mask, overlaps completely between classes
(box 1.3–23.0 %, sphere 8.8–17.1 %) and correlates with the error at r = −0.20
(`crease_indicator_rejected.md`). Since no gate exists, staying in scope is the
caller's responsibility, and the code says so rather than implying a guarantee
it cannot enforce.

**For the organoid application this is a weak constraint** — organoid surfaces
are smooth closed blobs, which is the qualified class. It is a real constraint
for anything faceted.

---

## 4. Validity layers, kept separate

| layer | status | evidence |
|---|---|---|
| mathematical basis | weights from a minimax LP over plane orientations; a-priori bound computable per measurement | `estimators_v3.py`, `lp_symmetry_verification.csv` |
| numerical solve | 63 of 64 spacing×radius solves converged; `m = 6` at spacing 1.2×1×1 did not, in 900 s / 80 cutting-plane iterations | `m6_exclusion.json` |
| oracle correctness | new closed-form areas/volumes agree with independent numerical quadrature to 1.5e-10 | `oracle_self_check2()` |
| development | 1083 cases, six shape classes, shape-class-resolved | `vv_dev2_all.csv` |
| confirmation | 480 cases executed once under hashed code; 96 in-domain | `vv_confirm2.csv`, `confirm2_verdict.json` |
| scope | crease exclusion is a stated limitation, not an enforced gate | `crease_indicator_rejected.md` |

`m = 6` is excluded globally, not per-spacing: production meets spacings absent
from this table and must solve at measurement time, so a radius that can hang
on one untested ratio cannot be offered for any. The exclusion was recorded
with its evidence before any confirmation case ran. The confirmation set used
only `m ∈ {4, 5}`.

---

## 5. Pre-registered predictions, checked after the run

| # | prediction | outcome |
|---|---|---|
| 1 | all criteria pass on the in-scope set, worst area < 2.5 % | **confirmed** — 0.951 % |
| 2 | creased out-of-scope cases violate the criterion near `rho_in` 10–15, in a 3–7 % band | **direction confirmed, band too narrow** — 8 violations observed, spanning 0.99–9.10 %, exceeding the predicted ceiling |
| 3 | `E0` violates on essentially every in-scope case | **confirmed** — 96 of 96, median 12.5 % |

Prediction 2 was partly wrong and is recorded as such: the crease error at the
gate reaches 9.10 %, not 7 %. This does not affect the verdict — those cases are
out of scope by declaration — but it means the crease deficit is somewhat worse
at moderate resolution than development led me to expect.

---

## 6. Margins, and where this is fragile

| quantity | worst on confirmation | criterion | margin |
|---|---|---|---|
| surface area | 0.951 % | 5 % | 5.3× |
| sphericity | 0.949 % | 5 % | 5.3× |
| **volume** | **0.761 %** | **1 %** | **1.31×** |

Volume is the thin one, and it is not a property of the new estimator: volume
is still voxel counting, unchanged from the legacy pipeline, and its error is
quantisation that scales like the surface-to-volume ratio. Development at the
same gate gave 0.556 %; confirmation reached 0.761 % on the torus, whose thin
tube is the worst surface-to-volume case in the set. A more demanding shape or
a lower `rho_in` would break the 1 % volume criterion before it broke the 5 %
area criterion. If volume accuracy matters, that is the next target — a
partial-volume or marching-cubes-on-distance-field volume would remove it.

Other boundaries worth stating plainly: `rho_in` was confirmed over 10.4–15.1
only, so the domain is verified in that band and asserted above it by
monotonicity of the development trend; anisotropy was confirmed at 1.2, 2.2,
3.5 and 4.0; and smooth cases below the gate reach 2.58 %, so the gate has real
headroom on area but is correctly placed for volume.

---

## 7. Reproducing this

    cd surfvv
    python step7_sym_solve.py        # weights (cached in out/crofton_weights_v3/)
    python step8_dev2_run.py dev2    # development grids: dev2, dev2b, dev2c
    python step9_freeze_v3.py --verify
    python step10_confirm2_run.py    # confirmation + adjudication

Development and confirmation call the identical code path (`run_grid`); the
refactor that made this so was verified to reproduce stored development rows
bit-for-bit before the hashes were taken. `step9_freeze_v3.py --verify` returned
`V3_hashes_match: true` immediately before the confirmation run.

## 8. Using it

    from surface_area_v3 import measure, in_domain

    m = measure(mask, spacing)          # no gating; returns rho_in, anisotropy
    ok, why = in_domain(m)              # computable gate only

`measure()` never gates and never silently returns a number it cannot stand
behind; the caller applies `in_domain()` and is responsible for the smoothness
scope. The legacy estimator `E0` is preserved unchanged, under its own name,
together with the evidence that it fails.

---

## Addendum: findings from integration into the project (2026-09-12)

The estimator was adopted as the project's default. Three things surfaced
during that work that the confirmation run could not have shown. All are
recorded here because they qualify claims made above; none of them changes
the verdict, which rests on the area criterion alone.

**1. The a-priori budget is predictive, not a proven bound.** The
plane-response term `t*` is a genuine worst-case bound for a planar surface
element, but the grazing term `D / r_in**2` models how curved and grazing
elements depart from planarity. Their sum held on all 1083 development and all
96 confirmation cases, and it is exceeded — by a factor of about 1.3 — by a
sphere centred exactly on a voxel centre with integer radius. Both the
development and confirmation grids placed every phantom at a random sub-voxel
offset, so this measure-zero placement was never sampled. At `rho_in = 12`,
unit isotropic spacing, it gives 1.81 % area error against a budget of 1.35 %.
The section 9 area criterion still holds there with a factor of 2.8 in hand, so
the qualification stands; the word "bound" does not. Report `apriori_bound` as
the expected error scale of a measurement, not as a guarantee about it.

**2. The same placement misses the volume criterion.** That sphere's
voxel-count volume is 1.18 % off, against the 1 % volume criterion, at a
resolution the declared domain admits. Volume is voxel counting and was not
touched by this work, so this is a property of the project's volume
measurement, not of the surface estimator — but it does mean the volume side of
the declared domain is not robust to lattice-aligned placement. The round-2
volume margin was already the thin one (1.31x versus 5.3x for area).

**3. The minimax optimum is not unique.** Re-solving the same LP returns a
different optimal vertex with the same `t*` and slightly different areas. The
qualified method is therefore the specific frozen weight vectors, not the
program that produces them. They ship with the package; an on-demand solve is
conforming but not evidence-bearing, and every measurement reports which tier
its weights came from.

A fourth item is a consequence rather than a limitation. Sphericity for small
spheres now slightly exceeds 1 (1.01–1.04 at radii of 4–6 voxels). The
superseded estimator's +10 % area error was cancelling the positive bias of
voxel-count volume in that ratio; with an accurate area the volume bias is
exposed instead. The project flags such values for review and never clips them,
which is the correct behaviour — the sphericity was previously closer to 1 for
the wrong reason.

Each of items 1, 2 and 4 is pinned as a test in
`tests/quantification/test_surface_qualified_domain.py` and
`tests/quantification/test_mask_features.py` rather than left in prose, so a
future change cannot quietly erase them.
