# Measurement validation update — 2026-09-12

A qualified surface-area estimator was developed and adopted as the production
default. This document records what was qualified, over what domain, what the
evidence is, and what did **not** change. The §9 criteria are unchanged
(<1 % volume, <5 % area for analytical shapes), and no gate criterion, phantom
or acceptance rule was modified.

Record: `docs/evidence/2026-09-12-surface-crofton-v3/evidence_manifest.json`
(class `canonical`, validated source commit `b659e2e89378`). Method
development and its confirmation run are reported separately in the round-2
report; this document covers the project's own frozen V&V grid.

## What was adopted

`crofton_minimax_sym_v3` — surface area from weighted counts of mask
transitions along lattice directions within a stencil of radius `m`, with the
direction weights obtained from a minimax linear program over plane
orientations. Two properties matter for validation:

- the weights have a stated mathematical basis, and the program's optimum
  `t*` bounds the estimator's plane response *before* any measurement;
- the restriction to weights constant on the orbits of the signed axis
  permutations preserving `diag(spacing)` is exact, not an approximation
  (the objective is invariant under that group and convex).

The stencil radius is selected per object to minimise a two-term a-priori
budget `t* + D / rho_in**2`. See
`src/organoid_analysis/quantification/surface_crofton.py` for the derivation
and the frozen constants.

## Qualified domain

`rho_in >= 10`, `anisotropy <= 4`, and **smooth closed surfaces**. Surfaces
with dihedral creases are `NOT QUALIFIED` at any resolution — this is a
property of the plane-response argument, not a resolution artefact, and no
amount of resolution removes it. The first two conditions are computed per
object and exported (`surface_in_qualified_domain`, `surface_rho_in`,
`surface_anisotropy`); objects outside them carry the
`surface_outside_qualified_domain` review flag. Smoothness cannot be checked
from a mask and is **not** asserted by the software; it is the caller's
responsibility, and for real organoid masks it is an assumption, not a fact.

## Results on the frozen V&V grid (188 phantoms, unchanged)

Per-case table: `surface_vv_dev.csv` in the record directory. `E0` in that
table is the harness's frozen slot label; the record's `estimator` block names
the method that actually ran.

| Area error | adopted | superseded |
|---|---|---|
| unrestricted worst | 25.81 % (`dev-078-cylinder-size3`, rho_in 1.50) | 18.74 % (`dev-139-cylinder-size24`) |
| unrestricted median signed | −1.24 % | +9.72 % |
| **inside the declared domain (n = 40)** | **worst 0.790 %, median 0.171 %** | worst 17.85 %, median 10.90 % |
| cases above 5 % inside the domain | 0 of 40 | 40 of 40 |

By shape class (worst / median absolute, adopted vs superseded): sphere
11.63 / 0.65 vs 17.45 / 12.27; ellipsoid 8.38 / 0.59 vs 18.33 / 10.98;
cylinder 25.81 / 4.65 vs 18.74 / 7.85.

The adopted estimator is more accurate on 160 of 188 cases and less accurate
on 28 — **every one of the 28 is outside the declared domain** (under-resolved
objects, and creased cylinders where its plane response is not valid). Its
unrestricted worst case is therefore larger than the superseded estimator's,
on a 3-voxel cylinder. This is a real trade: the method is far more accurate
where it claims to be and degrades faster where it does not.

## Gate status: unchanged, and still FAIL

**SG-1 remains FAIL** (25.81 % > 5 %). SG-1 is unrestricted by construction and
the frozen grid contains creased and under-resolved shapes the domain excludes.
The criterion was not relaxed, the grid was not subset, and no phantom was
removed. The domain-restricted claim is enforced separately, in
`tests/quantification/test_surface_qualified_domain.py`, which recomputes the
statistic above from the record's own per-case table and additionally asserts
the criterion on fresh in-domain phantoms.

**SG-2 (voxel-count volume) is bit-identical to the previous record** — worst
18.26 % for rho<3, 2.74 % for 3≤rho<6, 2.57 % for 6≤rho<12, PASS (0.528 %) for
rho≥12. Volume was not touched by this work. It remains the binding constraint
on the declared domain: inside the domain, area has roughly a factor of 6 in
hand against its criterion while volume has well under 2.

Overall gate outcome is unchanged: **NOT PASSED (0/6 items PASS)**.

## Consequences that must be read before using exported numbers

1. **Every exported area and sphericity changes.** Previously exported CSVs are
   not comparable with new ones. The reference sphere (radius 18 µm at
   (2,1,1) µm spacing) moves from +11.735 % to −0.561 % error. The superseded
   estimator stays reachable as `legacy_surface_area()` under its own identity
   so historical values remain reproducible.
2. **Sphericity for small objects now slightly exceeds 1.** The superseded
   estimator's ≈+10 % area error was cancelling the positive bias of
   voxel-count volume in the isoperimetric ratio; with an accurate area that
   bias is exposed instead (1.01–1.04 for spheres of radius 4–6 voxels; 1.61
   rather than 2.79 for a single voxel). Values above
   `SPHERICITY_REVIEW_LIMIT = 1.05` are flagged for review and never clipped.
   Sphericity closer to 1 before was closer for the wrong reason.
3. **The a-priori budget is predictive, not a proven bound.** It held on all
   1083 development and 96 confirmation cases, and is exceeded by a factor
   ≈1.3 by a sphere centred exactly on a voxel centre with integer radius —
   a measure-zero placement neither set sampled, since both offset every
   phantom sub-voxel. At rho_in 12 that case gives 1.81 % area error (criterion
   still met) and 1.18 % volume error (volume criterion missed, at a resolution
   the domain admits). Report `apriori_bound` as an expected error scale, not
   as a guarantee.
4. **The weight tables are the method.** The minimax optimum is not unique: a
   re-solve returns a different optimal vertex with the same `t*` and slightly
   different areas. The evidence-bearing vectors ship with the package; each
   measurement reports `weights_origin`, and anything other than `packaged`
   means an on-demand solve that is conforming but not evidence-bearing. See
   `src/organoid_analysis/quantification/crofton_weights/PROVENANCE.md`.

Items 2, 3 and the crease exclusion are pinned as tests rather than left in
prose, so a future change cannot quietly erase them.

## What this does not establish

Nothing here validates segmentation boundaries, acquisition calibration, or
any biological interpretation. The qualification is for surface area of
binary masks of smooth, adequately resolved, near-isotropic objects against
analytical phantoms. Whether a real organoid mask is inside that domain is a
judgement the exported per-object variables support but do not make.
