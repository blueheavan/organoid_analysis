# Analytical qualification protocol — SG-1a and SG-2a (frozen 2026-09-26)

Written and hashed (`freeze.json`) **before** the `confirm3` grid was generated
or executed, and committed in the source commit that the qualification record
validates. Nothing below may change after that commit; a change requires a new
protocol directory, a new confirmation grid and a new record.

Authority: owner decision D-14 (`docs/OWNER_DECISIONS.md`) authorises the
domain-restricted record contract; D-2 fixes the volume domain as the surface
domain *restricted further* only by evidence. This protocol narrows D-2 for
volume (§3); it never widens a domain.

## 1. Estimators under qualification

| Quantity | Production implementation | Identity |
|---|---|---|
| Surface area | `quantification.features.geometry()` → `surface_crofton.measure()` | `crofton_minimax_sym_v3`, packaged evidence-bearing weights |
| Volume | `geometry()`: filled-envelope voxel count × product of ZYX spacing | `voxel_count_v1` |

No estimator is changed by this protocol. Voxel-count volume is retained on
purpose: under a uniformly random sub-voxel translation it is a design-unbiased
volume estimator for any shape (Cavalieri point counting; Gundersen & Jensen,
*J. Microsc.* 147:229–263, 1987,
[doi:10.1111/j.1365-2818.1987.tb02837.x](https://doi.org/10.1111/j.1365-2818.1987.tb02837.x)),
so only its dispersion depends on resolution. The lattice-centred
counterexample (`docs/SCIENTIFIC_GATE_REASSESSMENT_2026-09-26.md`) proves that no
mask-only estimator can remove the worst-case ambiguity; a smoothing estimator
would add a shape-dependent bias without removing it.

## 2. SG-1a — surface domain (unchanged, frozen 2026-09-12)

`rho_in = r_in / max(spacing) ≥ 10`, anisotropy `≤ 4`, smooth closed surfaces,
and a Crofton weight vector that is evidence-bearing
(`surface_weights_evidence_bearing`). Machine-checkable part:
`surface_crofton.in_domain()` → `surface_in_qualified_domain`.

## 3. SG-2a — volume domain (new, narrower than D-2)

`volume_in_qualified_domain` = surface gate **and**
`rho_vol = r_in / cbrt(sz·sy·sx) ≥ 36`, smooth closed.
Derived in `docs/evidence/2026-09-26-volume-domain-development/` by the rule
declared in `dev_volume.py` before it ran (smallest ladder value with every
development case, including lattice-symmetric placements, `≤ 0.5 %`), and
checked above the ladder top by `dev_volume_ext.py`. No confirmation data were
generated or seen during development.

## 4. Confirmation grids

* **confirm2** — the 480-case grid frozen 2026-09-12 (`FREEZE_RECORD_V3.md`),
  regenerated from the byte-identical frozen `common2.py`/`common.py`
  (hash-checked against `freeze_record_v3.json`). Used for **SG-1a only**. Its
  volume errors were public before §3 was derived, so it cannot confirm §3; its
  in-volume-domain cases are reported, not adjudicated.
* **confirm3** — new, defined in `qualification_harness.py` (`CONFIRM3`):
  seed 20260927; the four smooth classes with aspects different from every
  development grid (ellipsoid 2.2:1.3:1, capsule L = 3.8 r, torus R = 2.3 r);
  six spacings with new absolute values whose ratios (1, 1.2, 2, 2.2, 3.5, 4)
  all carry evidence-bearing packaged weights; `rho_vol` ladder
  11, 15, 20, 27, 33, 37, 42; per cell five **lattice-symmetric placements**
  (centre on a voxel centre, corner, face centre, edge centre, Z-only corner;
  identity orientation) and two random placements (random rotation and
  sub-voxel offset). 1,176 cases. Used for SG-1a and SG-2a.

## 5. Acceptance rules — fixed now, applied once

Let `S1` = cases of confirm2 ∪ confirm3 that are smooth,
`surface_in_qualified_domain` and evidence-bearing. Let `S2` = cases of confirm3
that are smooth and `volume_in_qualified_domain`.

| Item | Rule | Validity (else VOID → INSUFFICIENT EVIDENCE) |
|---|---|---|
| SG-1a | every case in `S1`: `|area_rel_err| < 5 %` | each grid's part of `S1` contains all four smooth classes; confirm3's part contains symmetric and random placements; no case failed to measure |
| SG-2a | every case in `S2`: `|volume_rel_err| < 1 %` | `S2` contains all four classes, symmetric and random placements, and an isotropic case; no case failed to measure |

Worst case, not mean. Any violation ⇒ FAIL. The rules are implemented in
`src/organoid_analysis/validation/analytical_qualification_evidence.derive_results`
and re-derived from the raw CSV every time the gate verifies the record.

## 6. Predictions (stated before execution)

1. SG-1a PASS with worst area error below 2 % (development worst 1.33 %;
   confirm2 worst under the development implementation 0.95 %).
2. SG-2a PASS with worst volume error below 0.6 %.
3. confirm3 cases inside the surface gate but below `rho_vol 36` will include
   volume errors above 1 % at lattice-symmetric placements — the reason for §3.

## 7. What a PASS does and does not establish

It establishes that the production estimators meet the SCIENTIFIC_SPEC §9
criteria on analytical phantoms inside the stated domains, including the
placements that maximise lattice error. It does **not** establish segmentation
accuracy (SG-3), calibration of the physical spacing (SG-6), accuracy for
creased or hollow/open-cavity objects (not qualified), or any per-object
guarantee over every possible placement: the claim is an empirical worst case
over the frozen grids.

## 8. Deviations recorded before the freeze

1. **Harness smoke test.** Before `freeze.json` was written, the harness was
   exercised on three cases to check that it runs: `confirm2-0000`,
   `confirm3-0000` (sphere, `rho_vol` 11, outside the volume domain) and
   `confirm3-1175` (torus, `rho_vol` 42, in the volume domain; volume error
   −0.007 %, area error +0.72 %). No parameter, rule, domain constant or grid
   definition was changed afterwards. The three cases remain in the grid and
   are adjudicated like every other case.
2. **Development extension.** `dev_volume.py`'s rule selected the top of its own
   ladder (36). `dev_volume_ext.py` sampled torus and capsule at 36–48 before
   this protocol was written; the threshold was not changed by it.
