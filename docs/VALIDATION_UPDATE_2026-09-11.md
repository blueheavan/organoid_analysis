# Measurement validation and release-readiness update — 2026-09-11

This round focused on measurement, scientific and release validation, not
general refactoring. Every status below comes from commands executed in this
round, with evidence under
[`evidence/2026-09-11-measurement-vv/`](evidence/2026-09-11-measurement-vv/).
Earlier reports, README text and commit messages were not used as evidence.

## Verdicts (kept separate on purpose)

| Verdict | Result | Deciding evidence |
|---|---|---|
| **Engineering release readiness** | **NOT READY** | The engineering gate `pixi run ci` FAILS at `typecheck` (110 errors in 34 files). All executed tests pass (342 passed, 8 skipped, 8 deselected). Browser (vtk.js) rendering is NOT ASSESSED (8 Playwright tests skipped). Clean-environment install and dependency or secret scans were not run in this round. |
| **Publication / scientific-use readiness** | **NOT READY FOR THE SPECIFIED RESEARCH USE** | `pixi run science-gate` passes 0/6 items. Analytical surface area FAILS <5% in every stratum (worst 18.74%). Voxel volume FAILS <1% for ρ<12. Segmentation accuracy, assay validity, study-design validity and acquisition qualification are INSUFFICIENT EVIDENCE or NOT ASSESSED. |

Passing 342 tests is engineering evidence only. It is not scientific validation.

## 1. Source state and environment

- **Start:** branch `main`, HEAD `ed81d74ac9d8c3b7b2c5efe6e9e9549d00bc7856`,
  equal to `origin/main` (0 ahead / 0 behind after `git fetch`), clean working
  tree. `baseline/git-state.txt` shows one porcelain line: that is this
  evidence directory, created a moment before the state was recorded.
- **End:** uncommitted working-tree changes (no commit or push was requested).
  The source patch hash, untracked-file hashes, lockfile hash and versions are
  in [`snapshot.json`](evidence/2026-09-11-measurement-vv/snapshot.json).
- **Environment:** macOS 26.6.2 arm64; Pixi environment from `pixi.lock`;
  Python 3.12.14, numpy 2.5.2, scipy 1.17.1, scikit-image 0.26.0, pandas 2.3.3,
  statsmodels 0.14.6, tifffile 2026.3.3, cellpose 4.2.1.1, torch 2.7.1.
- **Audit separation: Tier D, LIMITED INDEPENDENCE.** The same agent and
  context implemented and reviewed the changes. The analytical oracles are
  closed-form expressions independent of the estimators. No independent person
  or independent agent reviewed this round.

## 2. Baseline before any change (`evidence/.../baseline/`)

| Check | Baseline result |
|---|---|
| `pixi run lint` | PASS (exit 0) |
| `pixi run check` | exit 0; no notebooks, so NOT APPLICABLE |
| `pixi run test` | **336 passed, 8 skipped, 8 deselected** (exit 0) |
| `pixi run typecheck` | **FAIL: 134 errors in 36 files**. By code: 59 import-untyped, 28 no-untyped-def, 10 arg-type, 8 assignment, 7 union-attr, 7 no-any-return, 5 var-annotated, 3 return-value, 2 index, 2 attr-defined, 1 operator, 1 call-overload, 1 call-arg. Top files: `organoid_measurement_workflow.py` 21, `report.py` 16, `inference.py` 13, `exploration.py` 13 |
| `pixi run ci` | **FAIL** (exit 1, at typecheck) |
| Sphere R = 18 µm, ZYX spacing (2,1,1) µm | volume −0.356%; **area +11.735%** (freshly reproduced) |
| <5% surface criterion | Present and unchanged in SCIENTIFIC_SPEC §9; the 15% bound in `test_geometry.py` is a regression tolerance |

## 3. Changes made

| Area | Change | Scientific effect |
|---|---|---|
| Surface V&V | Frozen plan, study script, 188-case development grid, selection record (`SURFACE_VV_PLAN.md`, `surface_vv.py`, `surface_vv_dev.*`, `SURFACE_SELECTION.md`) | None to production. The evidence is new. |
| Estimator versioning | `features.SURFACE_AREA_METHOD` (`marching_cubes_binary_lewiner_v1`, level, padding, spacing, step, smoothing) is exported in classical `provenance.json`, the multilevel summary and the Web feature bundle. A version-lock test is added. | Values unchanged. Future estimator changes cannot silently reinterpret history. |
| Web sphericity QC | Web and mask-feature rows with sphericity > 1.05 get `sphericity_above_geometric_range` (shared `SPHERICITY_REVIEW_LIMIT`). Feature schema 2.0 → 2.1. | A review flag only. No value, exclusion or threshold change. |
| Channel grid traceability | Classical `load_sample` records `channel_grid_verification` per channel: `same_file_as_primary`, `spacing_metadata_matches_primary` or `shape_only_no_spacing_metadata`. | Acceptance behavior unchanged. The unverified grid assumption now appears in provenance. |
| Gates | `pixi run science-gate` (`scripts/scientific_validation_gate.py`) is separate from `ci`. `test_geometry.py` names its 15% bound as a regression tolerance. | Makes the scientific FAIL explicit and machine-checkable. |
| Typing | Minimal fixes for genuine internal defects (§10). | None; behavior-preserving. |
| Docs | SCIENTIFIC_SPEC §15, VALIDATION_PLAN section, ALGORITHM_DECISIONS D14, PARAMETERS rows, README and VALIDATION_REPORT pointers, segmentation protocol. | The criteria are **not** changed. |

## 4. P1 — Surface-area measurement V&V

**Estimand.** The area of the continuous boundary of the solid whose voxel-centre
(Gauss) digitization produced the mask. This is not identifiable from a binary
mask, so every estimator has irreducible error at a finite resolution.

**Protocol.** [The plan](evidence/2026-09-11-measurement-vv/SURFACE_VV_PLAN.md)
was frozen, with SHA-256 recorded in `surface_vv_freeze.txt`, before any
candidate ran on a phantom.

- **Phantoms:** 188 cases. Spheres (R 3–36 µm, 3 sub-voxel offsets),
  ellipsoids (semi-axes s·(1, 0.75, 0.5), s = 6–48 µm), and flat-capped
  cylinders (r 3–24 µm, h = 3r), each at 4 orientations (identity plus 3
  seeded rotations).
- **Spacings (ZYX):** (1,1,1), (2,1,1), (3,1,1) and (2,0.7,0.7).
- **Oracles:** closed-form areas. The ellipsoid formula is cross-checked
  against numerical quadrature to a relative difference below 1e-8.
- **Criterion:** unchanged §9, applied **per case**. A stratum passes only if
  its worst case is below 5%.

| Estimator | Parameter provenance | Worst |err| | Median signed | Mean signed | Passing strata |
|---|---|---|---|---|---|
| E0 production binary marching cubes | level 0.5: analytical; no free parameter | **18.74%** | +9.72% | +9.90% | **none**, including large isotropic objects (9.79%) |
| E1 physical-space Gaussian + MC | σ = max spacing: **fixed engineering choice**, HEURISTIC, not calibrated | non-estimable ×17 | −6.64% | — | large (ρ≥12, max 4.36%) only |
| E2 signed distance + MC | level 0: analytical; no free parameter | 22.48% | +5.36% | +3.79% | none (+6–9% bias for large objects) |
| E3 Crofton, 13 directions | directions and Voronoi weights: analytical / literature-derived | 26.51% | −1.67% | −2.65% | sphere and ellipsoid with ρ≥3 (max 3.40%); **cylinders fail in every stratum** (up to 10.82% at ρ≥12) |
| E1a / E1b (σ 0.5× / 1.5×) | sensitivity only, **not selectable** | 71.3% / non-estimable | — | — | E1a: large and medium; E1b: none |

**Production E0 by stratum (worst |err|).**

- By resolution ratio ρ: very small (ρ<3) 18.53%, small 17.23%, medium
  18.74%, large 18.42%.
- By anisotropy: 1 → 9.79%, 2 → 18.53%, 2.857 → 18.42%, 3 → 18.74%.
- By shape: sphere 17.45%, ellipsoid 18.33%, cylinder 18.74%.

The bias is positive and systematic. Derived sphericity error ranges from
−15.8% to +7.3% (median −8.9%). All stratified tables, including shape × ρ and
spacing × ρ, are in
[`surface_vv_dev_summary.md`](evidence/2026-09-11-measurement-vv/surface_vv_dev_summary.md).

**Decision under the predeclared rule.**

- No candidate passes every development case, so the confirmation grid was not
  run and remains unseen.
- No candidate is better than E0 in every ρ stratum. E3 is worse for very
  small objects; E1 is non-estimable there.
- E1's result **depends on σ**. Picking σ = 0.5× after seeing these results
  would be post-hoc tuning, so it was not done.

**Result: production is unchanged and the surface criterion is FAIL.** The
historical definition is preserved and versioned. No independent calibration
set was needed, because nothing was calibrated.

**Voxel volume (same phantoms, §9 <1%).** The worst error is 0.53% for ρ≥12
(PASS), 2.57% for 6≤ρ<12, 2.74% for 3≤ρ<6 and 18.3% for ρ<3 (all FAIL).
Earlier "volume PASS" labels came from one large sphere. A nucleus of 5 µm
radius imaged with a 2 µm Z step has ρ = 2.5, which is in the worst stratum.

Options that need an owner decision (intended-use restriction, a new
predeclared candidate, grey-level surfaces) are in
[`SURFACE_SELECTION.md`](evidence/2026-09-11-measurement-vv/SURFACE_SELECTION.md).

## 5. P1 — Segmentation validation: INSUFFICIENT EVIDENCE

There is no qualified, independent 3D annotation set. The 5 local TIFFs have no
annotations and unknown acquisition provenance. The synthetic demo truth masks
are project-generated and share the segmenter's signal assumptions. The metric
code is correct but is not accuracy evidence.

Annotation-independence questions from the task:

| Question | Answer |
|---|---|
| Independent manual annotation? | None exists |
| Blind? | Not applicable, since there are no annotations |
| Multiple biological samples or imaging conditions? | None annotated |
| Overlap with model training data? | Not checked; no candidate set |
| Annotation SOP? | None |
| Inter-rater variability? | Not available |

The [validation protocol](evidence/2026-09-11-measurement-vv/SEGMENTATION_VALIDATION_PROTOCOL.md)
defines:
- reference-standard qualification (SOP, blinding, a training-overlap check,
  and a calibration/test split at the biological-sample level);
- inter-rater agreement as the bound on what can be demanded of the model;
- the biological sample as the unit for uncertainty, via a cluster bootstrap;
- object-level metrics (P/R/F1, FN/FP, split/merge);
- mask-level metrics (Dice/IoU, boundary distance, which must be implemented
  and verified first);
- downstream measurement bias for volume, equivalent diameter, surface area
  (same estimator on both masks), sphericity and axes, with Bland–Altman
  analysis;
- an acceptance framework whose numbers the owner sets from δ before any data
  are seen.

Dice alone can never validate morphology.

## 6. P1 — Biological and assay validity

| Output | Software correctness | Analytical validity | Assay validity | Biological validity | Clinical validity |
|---|---|---|---|---|---|
| Volume, equivalent diameter, axes | PASS (tested contracts and exact voxel arithmetic) | Volume **PARTIAL**: PASS for ρ≥12, FAIL below | INSUFFICIENT EVIDENCE (segmentation unvalidated) | NOT ASSESSED | NOT ASSESSED (not an intended use) |
| Surface area, sphericity | PASS (formula, units, flags) | **FAIL** (<5%) | INSUFFICIENT EVIDENCE | NOT ASSESSED | NOT ASSESSED |
| Calcein/PI viability states | PASS (state logic, abstention paths) | PARTIAL: background arithmetic tested; saturation limit defaults to the dtype maximum, not the detector bit depth | **NOT ASSESSED**: no control data or orthogonal assay; 0.30/0.60 gates and separation SNR are HEURISTIC with no calibration evidence | NOT ASSESSED | NOT ASSESSED |
| Nucleus/cell multilevel features | PASS (hand-computed phantoms) | PARTIAL (inherits surface FAIL; contact area exact on phantoms) | INSUFFICIENT EVIDENCE (registration and stain identity unqualified) | NOT ASSESSED | NOT ASSESSED |
| Phenotype classification (states, exploratory classifiers and clusters) | PASS (fold-local preprocessing) | INSUFFICIENT EVIDENCE | INSUFFICIENT EVIDENCE | NOT ASSESSED; object-level split and clustering-derived labels are circular | NOT ASSESSED |
| Drug response | No dose-response model exists; `dose` and `treatment` are pass-through metadata | NOT APPLICABLE | NOT ASSESSED | NOT ASSESSED | NOT ASSESSED |
| Morphology ↔ biological state | No such inference is implemented | NOT APPLICABLE | NOT ASSESSED | NOT ASSESSED | NOT ASSESSED |

Specific checks:
- **Threshold calibration:** none has calibration evidence.
- **Controls:** the rule requires `min_control_replicates` live and dead
  replicates per batch; whether that is adequate is not assessed.
- **Biological replicate:** a user-supplied manifest field, unique within a
  condition and global across batches; it is not validated.
- **Objects as replicates:** the classical inference uses a replicate random
  intercept or cluster-robust SE. The exploratory tools operate on object rows
  and are labeled exploratory.
- **Repeatability and reproducibility:** no evidence.

## 7. P1 — Statistics and study design: INSUFFICIENT EVIDENCE

- **Experimental unit.** The declared `biological_replicate` within a
  condition. The organoid is the measurement unit. Fields and wells are
  aggregated with equal weight per well (`aggregation.py`).
- **Technical vs biological replicates.** Fields → wells → replicates are
  distinguished. Repeated acquisition across batches does not create a new n.
- **Object-level rows.** `inference.py` models objects with a
  `condition::replicate` random intercept, falling back to OLS with
  cluster-robust SE and t(G−1). Wells nested within a replicate are **not**
  modeled. Paired donors and repeated conditions are unsupported.
- **Mixed-model conditions.** A minimum of 3 replicates per condition. The
  small-sample reference is a project heuristic, the omnibus test is
  asymptotic, and there is no type-I error or coverage evidence.
- **Multiplicity.** BH is applied within each feature's pairwise family, not
  across features. The family must be declared per study.
- **Leakage.**
  - Imputer and scaler are fitted within folds, and model selection uses
    training CV. This part is PASS.
  - Classifier splits are at the **object** level, so same-sample objects can
    fall on both sides. The code labels this "biological generalization NOT
    ASSESSED".
  - Cluster-importance on clustering-derived labels is circular.
  - No automated feature selection exists.
  - Viability calibration and classification use the same batch, but no
    accuracy is claimed from them.
- **Conclusion.** These are generic exploratory tools with no study design.
  **No publication-level inferential validity is claimed.**

## 8. P2 — Resampling, calibration and acquisition

| # | Question | Finding | Status |
|---|---|---|---|
| 1 | Spacing traceable? | Classical: manifest or OME with `spacing_source`; conflicts are rejected. Web: `metadata`, `default` or `user_override` recorded in the run config and bundle. Multilevel CLI: rejects mismatched spacings. | PASS (traceability) |
| 2 | Missing calibration | Classical rejects. Web proceeds with an assumed 0.414 µm XY and anisotropy 2.9, shows a warning, and exports `calibration_status: assumed_or_unknown`. | PASS (explicit) |
| 3 | Silent default calibration? | Not silent in any shipped route. The public `extract_mask_features` API still defaults to (1,1,1) µm; all internal callers pass spacing. | P3 finding, API unchanged |
| 4 | Preview/downsample affects measurement? | Preview is display-only; measurements use full resolution. Cellpose `xy_downsample` changes segmentation, restores masks by nearest neighbour, and is recorded. | PASS (tested contracts) |
| 5 | Downsampling sensitivity analysis | None | NOT ASSESSED |
| 6 | Unequal effective XY rejected on all resampling paths? | The Cellpose adapter rejects it. The classical and multilevel paths do not resample and use separate Y/X spacing. Covered by the resampling and inference contract tests in the passing suite. | PASS (tested paths) |
| 7 | Channel registration check | Spacing and shape equality only. Web requires explicit confirmation when grid metadata is partial. The classical path now records shape-only acceptance. There is no image-content registration check anywhere. | INSUFFICIENT EVIDENCE |
| 8 | Multi-channel mismatch → wrong fluorescence? | Possible when channels are misregistered with equal shape and spacing; software cannot detect it. | NOT ASSESSED |
| 9 | Saturation, dynamic range, background | The saturation flag defaults to the dtype maximum; 12-bit data stored as uint16 is never flagged without an explicit value. The local-shell background has no validation evidence. | INSUFFICIENT EVIDENCE (P2) |

Software unit tests are not used as a substitute for acquisition-level
evidence.

## 9. P2 — Large volumes and runtime robustness ([probe](evidence/2026-09-11-measurement-vv/runtime_probe.json))

| Condition | Observed in this round | Status |
|---|---|---|
| Large Z-stack / object count | 64×1024×1024 uint32 (256 MiB), 3000 objects: 6.9 s, **peak RSS 2.59 GB (about 10× the array)** | PARTIAL: one synthetic size; memory amplification noted (P3). Full-volume Cellpose memory NOT ASSESSED. |
| Sparse labels | IDs 1, 2^31+7, 2^32−2 preserved | PASS |
| uint16 vs uint32 | Identical feature tables | PASS |
| Extreme anisotropy (25×) | Accepted **without any warning**; area +26.4%, volume −1.6% | FAIL as out-of-domain behavior (P2). No domain threshold is chosen here. |
| Empty / degenerate objects | 0 rows for an empty mask. Single voxel, plane and line give NaN solidity with flags. Single-voxel sphericity 2.79 was **unflagged**; now flagged (defect D-1). | PASS after the fix |
| Concurrent Streamlit sessions | Inference lock exists; no multi-session test | NOT ASSESSED |
| Cache invalidation | Existing AppTest layer and geometry-change tests pass in the full suite | PASS (tested cases) |
| Export consistency | Identical inputs give byte-identical ZIP bundles | PASS |
| Native VTK rendering | `pixi run test-render`: 8 passed on this host ([log](evidence/2026-09-11-measurement-vv/render.log)) | PASS on this host only |
| Browser vtk.js rendering | 8 tests skipped (Playwright not installed) | NOT ASSESSED |

## 10. P3 — mypy typing debt

| | Errors |
|---|---|
| Baseline | 134 in 36 files |
| Final | 110 errors in 34 files |
| New errors vs baseline | 0 (message-level diff against the baseline log) |

The 24 errors removed were: 7 assignment, 5 arg-type, 4 union-attr, 2
return-value, 2 var-annotated, 1 attr-defined, 1 call-overload, 1 index and 1
operator.

Fixes made, all genuine internal defects:
- Optional control flow: the `organoid_workspace` config narrowing.
- Reused variables with changed types: `watershed_instances`,
  `viewer_payload`.
- A heterogeneous provenance dict: typed locals in the classical workflow.
- Imprecise tuple lengths: `tiff_contract`, `multilevel_results`.
- A wrong Optional return annotation: `classify_axes`.
- An int/float dict: `aggregation`.

No `Any`, `cast`, `# type: ignore` or placeholder Protocol was added, and no
numerical behavior changed (full suite passes).

Remaining, by cause:
- **By code (110):** 59 import-untyped (pandas 21, scipy 19, 15 imports of
  installed packages that ship no `py.typed` marker, including statsmodels,
  scikit-learn and cellpose, PyYAML 3, seaborn 1), 28 no-untyped-def,
  7 no-any-return, 5 arg-type, 3 var-annotated, 3 union-attr, and 1 each of
  return-value, index, call-arg, attr-defined and assignment.
- **Third-party boundary:** all import-untyped (pandas, scipy, PyYAML,
  statsmodels, scikit-learn, cellpose, seaborn), plus tifffile, matplotlib,
  pyvista and vtk signature mismatches. Recorded as a typing boundary; no stubs
  were forged.
- **Internal debt:** untyped function signatures (no-untyped-def) and Any
  returns from untyped libraries.

## 11. CI and release gates

- **Engineering regression gate:** `pixi run ci` (lint → test → check →
  typecheck). It covers tests, lint, determinism, the export contract and
  interface regressions. **FAILS because of typecheck**, so CI is not described
  as PASS.
- **Scientific validation gate:** `pixi run science-gate`. It covers analytical
  surface (SG-1), analytical volume (SG-2), segmentation (SG-3), assay (SG-4),
  statistics and study design (SG-5), and acquisition (SG-6). It exits nonzero
  until every item has PASS evidence. *(Superseded later on 2026-09-11. The original
  single reference-sphere staleness check was replaced by a sealed evidence
  manifest. CI was split into a required regression job and a mypy ratchet.
  See [VALIDATION_RECORDS.md](VALIDATION_RECORDS.md).)* It is intentionally **not**
  part of `ci` or the GitHub workflow.
- The 15% sphere assertion remains as `SURFACE_REGRESSION_TOLERANCE`. The
  <5% criterion lives only in the specification and the science gate.

## 12. Confirmed software defects (this round)

| ID / severity | Evidence → root cause | Fix and regression test | Status |
|---|---|---|---|
| D-1 P3 | Probe: a single-voxel Web row had sphericity 2.79 and no flag, although SPEC §10.3 says values above 1.05 are flagged. The limit existed only in the classical route. | Shared `SPHERICITY_REVIEW_LIMIT`; Web flag; schema 2.1. `test_discretized_sphericity_above_review_limit_is_flagged_not_clipped`; one existing expectation updated to include the new flag. | Fixed |
| D-2 P2 | `load_sample` accepted a separate channel TIFF without spacing metadata by shape equality alone, and left no record distinguishing it from a metadata-checked channel. | `channel_grid_verification` in sample metadata and `provenance.json`; acceptance unchanged. `test_channel_grid_verification_distinguishes_shape_only_from_metadata_checked`. | Fixed (traceability); registration itself is still unqualified |
| D-3 P3 | The surface estimator was identified only as free text in the Web bundle and was absent from classical and multilevel provenance, which risked silent reinterpretation. | `SURFACE_AREA_METHOD` in all three routes; `test_surface_estimator_version_lock`; bundle provenance assertion. | Fixed |
| D-4 P4 | Genuine typing defects (§10) | Behavior-preserving fixes | Fixed; full suite passes |
| — | My own new gate-test loader did not register the module in `sys.modules`, so the first final run had 2 failures | Test harness fixed; the invalid run is preserved in `final-attempt1-invalid/` | Not a product defect |

## 13. Scientific findings and uncertainties (not fixed by code)

| ID | Finding | Severity / status |
|---|---|---|
| F-1 | Analytical surface area fails <5% in every stratum | **P1, FAIL**. Release and publication blocker. |
| F-2 | Voxel volume fails <1% for ρ<12, including typical nuclei at a 2 µm Z step | **P1** for small-object or nucleus volume claims; **FAIL** |
| F-3 | No warning for extreme anisotropy or low ρ | P2, FAIL (out-of-domain behavior) |
| F-4 | Saturation limit defaults to the dtype maximum | P2, INSUFFICIENT EVIDENCE |
| F-5 | Segmentation accuracy unvalidated | P1 claim risk, INSUFFICIENT EVIDENCE |
| F-6 | Viability assay and thresholds unvalidated | P1 claim risk, NOT ASSESSED |
| F-7 | Study design and inferential calibration | P1 claim risk, INSUFFICIENT EVIDENCE |
| F-8 | Mask-feature memory is about 10× the label array | P3 |
| F-9 | `extract_mask_features` API default spacing is (1,1,1) µm | P3 |

## 14. Acceptance criteria, commands and results

| Criterion (unchanged) | Command | Result | Status |
|---|---|---|---|
| Area: each analytical case <5% | `pixi run python docs/evidence/2026-09-11-measurement-vv/surface_vv.py --grid dev` | worst 18.74% (E0) | FAIL |
| Volume: each analytical case <1% | same | worst 18.3% (ρ<3), 0.53% (ρ≥12) | FAIL (PASS for ρ≥12 only) |
| Segmentation: qualified independent reference | — | no reference set | INSUFFICIENT EVIDENCE |
| Assay: orthogonal reference | — | none | NOT ASSESSED |
| Inference: predeclared design, type-I error and coverage | — | none | INSUFFICIENT EVIDENCE |
| Engineering: tests | `pixi run test` | 342 passed, 8 skipped, 8 deselected | PASS |
| Engineering: lint | `pixi run lint` | All checks passed | PASS |
| Engineering: notebooks | `pixi run check` | no notebooks | NOT APPLICABLE |
| Engineering: types | `pixi run typecheck` | 110 errors in 34 files | FAIL |
| Engineering gate | `pixi run ci` | exit 1: lint, test and check pass; stops at typecheck (110 errors) | FAIL |
| Scientific gate | `pixi run science-gate` | 0/6 PASS, exit 1 | FAIL |
| Native render | `pixi run test-render` | 8 passed | PASS (this host) |

## 15. Changes compared with the previous version

- **Surface.** Before: one sphere, FAIL. Now: 188 cases, 4 estimators, FAIL
  with strata. Worst E0 error went from 11.7% (single case) to 18.7% across
  the grid.
- **Volume.** "PASS" is corrected to FAIL for ρ<12.
- **Provenance.** The surface method is versioned and exported in all routes.
  Channel grid verification is recorded.
- **Web features.** The sphericity review flag is added (schema 2.1).
- **Gates.** The scientific gate now exists, separate from the engineering gate.
- **Typing.** 134 → 110 errors.
- **Tests.** 336 → 342 passing.
- **Unchanged.** Algorithms, thresholds, statistical models, criteria and
  intended use.

## 16. Remaining blockers

- **Engineering release:**
  - mypy errors (§10);
  - browser rendering NOT ASSESSED;
  - no multi-session concurrency test;
  - no clean-environment install, dependency or secret scan in this round;
  - large-volume memory not qualified.
- **Publication / scientific use:**
  - F-1 and F-2 (analytical accuracy);
  - F-5 (segmentation);
  - F-6 (assay);
  - F-7 (study design);
  - acquisition qualification (§8 items 7–9);
  - Tier D review only.

## 17. Evidence files (`docs/evidence/2026-09-11-measurement-vv/`)

`baseline/` (git state, lint/check/typecheck/test/ci logs) ·
`SURFACE_VV_PLAN.md`, `surface_vv_freeze.txt`, `surface_vv.py`,
`surface_vv_dev.csv`, `surface_vv_dev_summary.{json,md}`, `surface_vv_dev.log`,
`SURFACE_SELECTION.md` · `SEGMENTATION_VALIDATION_PROTOCOL.md` ·
`runtime_probe.py`, `runtime_probe.{json,log}` · `render.log` ·
`final/` (exits, lint/check/typecheck/test/ci/science-gate logs, `pytest.xml`) ·
`final-attempt1-invalid/` · `snapshot.json`.

## 18. Diff self-review (Tier D)

Checked against the complete diff:
- No scientific definition changed. `surface_area_um2`, volume, sphericity,
  thresholds, statistical models and intended use are unchanged. The only
  output changes are an added QC flag and added provenance fields.
- No criterion leakage. The criteria are unchanged, and the plan was frozen
  and hashed before any results.
- No undeclared tuning. The only tunable candidate parameter was predeclared
  and was not selected afterwards.
- No silent change to historical results. The version-lock test and the
  provenance field protect the estimator definition.
- Measurements can be rebuilt from provenance: the method record, spacing and
  its source, package versions, source hashes and input hashes.
- The documentation matches the implementation, as checked for the files
  listed in §3.
