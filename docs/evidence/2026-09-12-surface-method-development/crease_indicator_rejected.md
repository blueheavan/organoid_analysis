# Rejected: a mask-computable crease indicator

Development showed the estimator's only substantial error mode is a systematic
negative bias at dihedral creases. If creases were detectable from the mask
alone, the applicability domain could gate on them automatically instead of
relying on a scope statement the caller must honour. This records the attempt
and why it was rejected — before the confirmation run, and without consulting
it.

## Requirement

A crease indicator must be computable from **one binary mask plus its voxel
spacing**, with no truth value and no shape label, since that is all a
measurement has.

## Candidate tested

The disagreement between the two estimators already computed on every mask:

    gap = area_E0 / area_v3 - 1

where `E0` is binary-mask marching cubes. The reasoning: marching cubes and
Crofton respond differently to creases, so their disagreement might carry the
crease signal. `gap` needs no truth and costs nothing extra.

## Result — rejected

Across 339 development cases with `rho_in ≥ 6`, anisotropy ≤ 3:

| class | n | `gap` min | median | max |
|---|---|---|---|---|
| box | 63 | 1.29 % | 12.28 % | 22.99 % |
| cylinder | 57 | 5.20 % | 12.10 % | 24.43 % |
| sphere | 36 | 8.82 % | 12.75 % | 17.06 % |
| ellipsoid | 62 | 8.84 % | 10.75 % | 14.68 % |
| capsule | 61 | 6.51 % | 9.27 % | 18.28 % |
| torus | 60 | 8.52 % | 13.26 % | 18.77 % |

The distributions overlap almost completely and the medians are ordered wrongly
— the torus, a smooth surface, has the largest median disagreement of any
class, while the box spans nearly the entire range. No threshold separates
creased from smooth at any `rho_in` floor tested (6, 8, 10): the smooth maximum
always exceeds the creased minimum.

The indicator also fails to predict the error it is meant to stand in for:
across creased cases, `corr(gap, −signed area error) = −0.20`.

Gating on `gap ≤ G` was evaluated directly for `G ∈ {4.0 … 7.0}` and never
produced a clean survivor set: every threshold either kept violating cases
(worst survivor 5.9–8.5 %) or discarded most of the smooth classes along with
the creased ones.

## Consequence

No computable crease gate is declared. The crease limitation is carried as a
scope-of-intended-use statement in `FREEZE_RECORD_V3.md` §3–4, explicitly
flagged as the caller's responsibility rather than an enforced guarantee, and
creased phantoms remain in both grids with their results reported.
