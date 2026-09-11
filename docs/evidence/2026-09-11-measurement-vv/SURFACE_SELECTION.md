# Surface estimator selection record — 2026-09-11

Written after the development grid (`surface_vv_dev_*`) and before any production
change, under the frozen [plan](SURFACE_VV_PLAN.md). The plan and script hashes
are in `surface_vv_freeze.txt`. No amendment to the plan was made.

## Outcome under the predeclared rule

| Estimator | Eligible | Overall (188 cases, every case <5%) | Largest absolute error | Strata that pass (max <5%) |
|---|---|---|---|---|
| E0 `marching_cubes_binary_v1` (production) | yes | **FAIL** | 18.74% | none, including large isotropic objects (9.79% at (1,1,1), ρ≥12) |
| E1 Gaussian, σ = 1.0 × max spacing | yes | **FAIL** | non-estimable in 17 cases (object vanished) | large (ρ≥12) only: 4.36% |
| E2 signed distance | yes | **FAIL** | 22.48% | none: consistent +6–9% positive bias for large objects |
| E3 Crofton 13-direction | yes | **FAIL** | 26.51% | sphere and ellipsoid for ρ≥3 (max 3.40%). Cylinders fail in every stratum, up to 10.82% at ρ≥12 |
| E1a σ = 0.5× (sensitivity only) | **no** | FAIL | 71.30% | large and medium |
| E1b σ = 1.5× (sensitivity only) | **no** | FAIL | non-estimable in 36 cases | none |

- **Rule 2 (replace the default only if every development case passes): not
  triggered.** No candidate passes, so the confirmation grid is **not run**.
  Running it could not change the decision, and leaving it unseen keeps it
  available as an independent set for a future predeclared candidate.
- **Rule 3 (recommend a candidate that is better in every ρ stratum): no
  candidate qualifies.** E3 has a lower maximum error than E0 in the large,
  medium and small strata but a higher one in the very-small stratum (26.51% vs
  18.53%). E1 is non-estimable for very small objects. E2 is worse than E0 for
  very small objects.
- **Decision: production `surface_area_um2` keeps its historical definition
  (E0) unchanged.** The SCIENTIFIC_SPEC §9 surface criterion (<5%) stays
  **FAIL**. Rule 4 is implemented: the estimator identity is now exported as a
  versioned provenance record, so any future change cannot silently
  reinterpret historical values.

## Findings (descriptive, not criteria changes)

1. The production bias does not come from small objects alone. E0
   overestimates area by about +8 to +13% (median) in every ρ stratum. The
   worst case is +18.7% at anisotropy 3. This matches the known non-convergence
   of terraced binary marching cubes. Derived sphericity is biased by −15.8% to
   +7.3% (median −8.9%).
2. For E1 the scientific conclusion **depends on the parameter**. With σ = 0.5×
   the large and medium strata pass. With the predeclared σ = 1.0× only the
   large stratum passes. With σ = 1.5× no stratum passes. Because σ has no
   analytical provenance, choosing 0.5× now, after seeing these results, would
   be post-hoc tuning. That is prohibited by the plan and by the task
   instructions.
3. E3's cylinder failures follow the analytical plane response computed before
   the grid (−15% to +10% for anisotropic lattices). Flat faces aligned with
   lattice planes are the adverse case. Real organoids rarely have flat faces,
   but the §9 criterion is stated for analytical shapes without a shape
   restriction, and this round does not introduce one.
4. **Voxel-count volume also misses the §9 <1% volume criterion outside the
   large stratum.** Maximum absolute errors are 0.53% (ρ≥12, PASS), 2.57%
   (6≤ρ<12), 2.74% (3≤ρ<6) and 18.3% (ρ<3). Earlier reports marked volume PASS
   from one large sphere. That PASS is valid only for that object size.

## Options that need a project-owner decision (not taken here)

- **Intended-use domain restriction**, for example a minimum resolution ratio
  of ρ≥12 for any quantitative surface or volume claim. This changes intended
  use and must be decided by the owner. Even with that restriction, E0 still
  fails (18.42% in the large stratum).
- **A new predeclared candidate.** One example is Crofton with a denser
  direction set, which reduces the flat-face plane-response error. It would need
  a new frozen plan, and the untouched confirmation grid, or a fresh grid, as
  its test set. It was motivated by these results, so it is not independent of
  this development set.
- **Surface estimation from grey-level or probability maps.** This would
  recover the sub-voxel information that binarization discards. It needs
  pipeline changes and a separate V&V.
