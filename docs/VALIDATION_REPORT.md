# Full scientific workflow audit and release assessment — 2026-09-10

Latest scoped update: [2026-09-11 workflow optimization](WORKFLOW_OPTIMIZATION_2026-09-11.md)
corrects Web measurement support/precision and rounded-grid resampling, adds QC
and traceable feature downloads, and records fresh engineering checks. The
2026-09-10 execution counts and source-state statements below are historical;
the broader surface-accuracy and biological/inferential evidence limitations
remain unresolved.

**Overall verdict: NOT READY FOR THE SPECIFIED RESEARCH USE.**

Scope: the current local repository's research workflow, **Z-stack input → preview/inspection → segmentation → quantitative features → analysis/statistics → export**, including classical CLI, Cellpose/Web, multilevel labels and exploratory statistics. This supersedes previous local readiness claims. The intended use and original accuracy criteria were retained. A documented core surface-accuracy criterion fails; representative biological and inferential evidence is insufficient. Passing software tests is not a positive scientific release verdict.

## 1. Source state, evidence and audit separation

Baseline HEAD: `ab6ee90db936de2a0888b9600ddf62f1729b4a1e`. The initial worktree already contained 19 modified tracked files, six untracked regression-test files and prior evidence logs. Those were preserved as audit inputs and revalidated, not credited as new work or accepted as proof. The initial run produced **284 passed, 1 failed, 8 skipped, 8 deselected**; the failure was error-message compatibility, not acceptance of nonfinite intensity. Initial tracked changes were captured; original untracked contents were not separately frozen before extensions, so the initial snapshot is not fully reconstructible from its patch alone.

[Reproducibility record](REPRODUCIBILITY.md) and [final snapshot](evidence/2026-09-10/final-snapshot.json) identify final files, source patch, environment, inputs and built artifact. Evidence in `evidence/2026-09-09` is historical only; fresh results are in `evidence/2026-09-10`.

Review separation is **Tier D — LIMITED INDEPENDENCE**: same agent/context performed implementation and review. No independent person, independent agent, external CI or independent biological validation is claimed. Analytical/unit-conversion oracles are independent of implementation formulas where possible; this does not create reviewer independence. Mock-model tests establish adapter contracts only.

## 2. Scientific impact and actual data flow

[SPEC §13](SCIENTIFIC_SPEC.md#13-current-implementation-audit--2026-09-10-governing-qualification) maps exact components to impact and semantics. Maximum impact is **S3**: default-enabled condition LMM/clustered OLS and exploratory classifiers can affect scientific inference. I/O calibration, segmentation, measurement and filtering are S1–S2. S4 decision-critical use was not requested; its decision-risk/regulatory gates are NOT APPLICABLE.

The project has multiple cooperating routes, not one interchangeable implementation:

- Classical `analyze` uses a manifest, structural watershed/probability/imported labels, filled-envelope geometry, raw fluorescence/background correction, eligibility filtering, well/replicate summaries and default condition inference.
- Web Cellpose segments nuclei and optionally cells; the multilevel runner additionally requires independently supplied, registered organoid labels. Preview-normalized data are not reused as quantitative raw intensities.
- `analyze-3d` measures raw registered organoid/cell/nucleus masks, overlap hierarchy, face contacts, EDT depth, intensity moments and QC flags. Its aggregates retain flagged objects.
- `cells` performs a separate one-to-one pairing/filtering operation. Web `mask_features` uses filled size/axes and raw centroids/solidity and does not apply classical geometry exclusions. Exploratory statistical tables must not be mistaken for independent biological units.

## 3. Engineering QA — final executed results

| Check | Observed result | Evidence status / artifact |
|---|---|---|
| Full default suite | **307 passed, 8 skipped, 8 deselected**, 887 warnings | PASS for executed cases; [full-tests.log](evidence/2026-09-10/full-tests.log), [JUnit](evidence/2026-09-10/pytest.xml). Skipped tests are not passed. |
| Native render suite | Sandbox segfault in VTK screenshot; host graphics retry **8 passed** | PASS on tested host, PARTIAL across environments; [render-host.log](evidence/2026-09-10/render-host.log). |
| Browser rendering | 8 tests skipped: Playwright unavailable | NOT ASSESSED. Native VTK does not validate browser vtk.js interaction. |
| Lint | All checks passed | PASS; [lint.log](evidence/2026-09-10/lint.log). |
| Typecheck | **193 errors in 37 files**; freshly reconstructed HEAD had 194 | FAIL; [typecheck.log](evidence/2026-09-10/typecheck.log), [HEAD log](evidence/2026-09-10/typecheck-head.log). Missing stubs/untyped functions and existing type errors remain. Corrected six input-worktree diagnostics without suppressing checks. |
| Format check | **69 files would be reformatted; 33 already formatted** | FAIL; [format-check.log](evidence/2026-09-10/format-check.log). No repository-wide cosmetic rewrite was made. |
| Notebook validation | No notebooks found | NOT APPLICABLE; [notebooks.log](evidence/2026-09-10/notebooks.log). |
| Wheel build and installation | Declared setuptools backend builds; uv installs actual wheel into isolated target; installed source/static bytes match | PASS for package artifact using current dependencies; [artifact-final.log](evidence/2026-09-10/artifact-final.log). Not a fresh dependency bootstrap. |
| Artifact runtime | Import outside checkout; bundled offline vtk.js; analytic volume; CLI version; demo and real-file analysis | PASS for execution; same artifact log and [paths/hash](evidence/2026-09-10/artifact-paths.json). |
| Web startup | Streamlit AppTest: no startup exceptions | PASS for headless startup only; [streamlit-startup.log](evidence/2026-09-10/streamlit-startup.log). |
| Real Cellpose smoke | MPS, 3×256×256 PDAC C1/C2 crop, 29 nuclei and 23 cells | PASS for execution, not accuracy; [model-smoke-host.log](evidence/2026-09-10/model-smoke-host.log). Sandbox attempt stalled at font-cache setup and was interrupted. |
| Installed dependencies | Duplicate VTK 9.6.2 dist-info and egg-info | FAIL metadata consistency; [dependency-check.log](evidence/2026-09-10/dependency-check.log). No evidence of two different runtime versions. |
| Diff review | Complete tracked changes and all added regression modules re-read; whitespace check passes | PASS for local review scope; Tier D limitation applies. |

The `ci` task's component checks were executed individually so a failure could not hide later checks. No lint/type/format/test settings were relaxed. No vulnerability scan, secret scan or complete SBOM was performed (NOT ASSESSED). Existing test warnings include dependency deprecations and intentionally exercised degenerate statistical fits; the count does not establish benignity of every warning in arbitrary runs.

## 4. Confirmed defects and remediation

All rows below are current-code/source-diff observations plus reproducible tests. “Input edit” means the fix was already in the initial worktree and was reviewed/retested. “Audit repair” means this turn implemented or completed it. Scientific algorithms, calibrated cutoffs, statistical model families and research acceptance criteria were not tuned to improve results.

| ID / severity | Evidence and root cause → affected output | Minimal correction / regression evidence | Resolution |
|---|---|---|---|
| R1 P1 | ImageJ spacing ignored units; reference-frame plane coordinates could become µm; absolute Z differences accepted direction reversals | Input edit: explicit unit conversion and monotonic uniform-grid validation. Audit repair: precise optional types. `test_metadata_contract_regressions.py` verifies nm/mm/µm, unknown units, invalid/reversed/descending Z. | PASS tested contracts; registration/calibration still externally unqualified. |
| R2 P1 | Web reader selected frame 0 without explicit multi-time choice | Input edit: explicit `time_index`, selected-frame provenance, reject ambiguous upload. Audit repair: reject negative/bool/fractional series/time selectors and complex quantitative input. Metadata regression tests. | PASS; Web multi-time selector remains unavailable, safely rejects. |
| R3 P2 | Requested XY downsampling changed image scale without changing diameter; invalid inputs could reach the model | Input edit: scale diameters with requested factor, validate model/config/finite input and output masks. `test_inference_contract_regressions.py`. | PASS adapter contract; odd-grid resampling equivalence NOT ASSESSED. |
| R4 P2 | Mixed int64/uint64 overlap pairs promoted to float64 and merged IDs above 2^53; cell pairing reproduced 1 pair instead of 2 | Input edit: hierarchy preserves integers. Audit repair: structured integer pairs in `cellular_measurements.py`. `test_audit_contracts.py` and multilevel numerical regressions. | PASS exact identity without altering pairing objective. |
| R5 P2 | E[x²]−E[x]² cancellation erased nuclear/pooled variance at large additive offsets | Input edit: centered residual sums and within/between-nucleus pooling. Analytical variances 5.25 and 9.25 remain invariant under +1e9. | PASS tested numerical range; extreme overflow remains unassessed. |
| R6 P1 | Invalid log-volume clipping and preprocessing/model-selection leakage could produce misleading inference/performance | Input edits: reject nonpositive/Inf log inputs; collision-free replicate keys; train-fold imputation/scaling; choose model by training CV. Statistical contract tests. | PASS software controls; biological-unit leakage and statistical calibration remain insufficient. |
| R7 P2 | Undefined effect/assumption/cluster statistics were replaced by 0 effects or p=1; zero-mean CV became 0 | Audit repair: missing numeric results, Not estimable/NOT ASSESSED, nullable significance; inadequate Shapiro domain abstains rather than inventing normality. Constant/small/zero-mean regression cases. | PASS; valid-data method selection unchanged. |
| R8 P2 | NaN/Inf/fractional/bool multi-level QC values could disable flags or change count meaning | Audit repair: finite real thresholds and integral voxel count, without changing valid cutoffs. Pre-fix 8/10 probe cases failed; post-fix all pass. | PASS. |
| R9 P2 | uint32 exports could wrap IDs; partial saves could appear complete; RGB-prone label shapes were written with TIFF RGB defaults | Input edits: range validation, incomplete sentinels/completion ordering. Audit repair: shared calibrated OME ZYX writer and singleton-Z restoration. Fault injection, max-ID round trip, (3,4,4) grayscale check, existing saved-run recovery tests. | PASS. First OME change exposed singleton squeezing; fixed before final suite, not hidden by changing the test. |
| R10 P2 | Web mask-feature `regionprops` allocated by largest sparse ID; accepted bool masks then failed in regionprops; empty means falsely became zero | Audit repair: compact only internally and map source IDs back; bool foreground works; empty mean/range/shape summaries are NaN, count/total remain zero. Guarded allocation test avoids a dangerous pre-fix allocation. | PASS. |
| R11 P2 | Native VTK flattened ZYX in Fortran order despite XYZ grid dimensions; uint8 surface labels truncated IDs | Input edit: C-order layout and uint32 IDs. Asymmetric-coordinate/sparse-ID tests and host native rendering. | PASS tested coordinates; renderer not quantitative ground truth. |
| R12 P2 | Incomplete source hashes, unvalidated uploaded organoid grid, and missing packaged vtk.js weakened reproducibility | Input edits: complete relative package hash keys, mask/input hashes, grid checks, wheel package-data. Audit repair: type narrowing. Export/source hash tests and installed-wheel runtime. | PARTIAL: software capture improved; complete assay/model/study lineage still absent. |

Pre-fix reproductions and intermediate failures are retained as `defect-reproduction.log`, `statistics-defect-reproduction.log`, `mask-defect-reproduction.log`, `index-defect-reproduction.log` and `intermediate-full-tests.log`. Do not mistake them for final results. Removed `_CLASSIC_SOURCE_FILES`/`_PACKAGE_ROOT` lists were internal provenance-only implementation, replaced by complete package enumeration; current callers were checked, no scientific/public capability deleted.

## 5. Algorithm selection and scientific basis

**Conclusion: PARTIAL method implementation support; INSUFFICIENT EVIDENCE for project-data suitability.** [Algorithm decisions](ALGORITHM_DECISIONS.md) map actual versions, definitions, alternative methods, assumptions, failure modes and claim-specific sources. Official Otsu/watershed/Cellpose/EDT/marching-cubes documentation supports operations and API semantics, not the selected physical cutoffs or biological targets. Lewiner marching cubes is the installed implementation; 3D flow-error thresholds are inactive. Mathematical moments/contact/volume have analytical support; overlap, core/periphery and QC are operational definitions.

No independent method comparison establishes the selected model or watershed settings as appropriate for the local assay. Nuclear texture does not by itself establish an organoid boundary. Visual inspection of the real spheroid QC showed multiple bright internal structures receiving instance outlines; this is a reason to qualify target identity, not an annotated determination of which objects are correct. No algorithm was replaced on that visual judgment.

## 6. Parameter / threshold scientific basis

**Conclusion: INSUFFICIENT EVIDENCE.** The [113-row registry](PARAMETERS.md) separates Origin, basis strength and evidence status and records missing calibration/independence/sensitivity/uncertainty. Values such as 2000 µm³ minimum organoid volume, 1 µm smoothing, 15 µm seed distance, 0.5 probability/overlap cutoffs, 3.5 MAD score, 0.5/0.8 radial bins, 0.30/0.60 viability gates, two control replicates and three inference replicates are not justified by a software default or a method citation. No project calibration/held-out artifact establishes them. Model normalization, minimum/maximum mask size, iteration count and tiling defaults were inspected in installed source. No calibrated biological sensitivity range was available; sensitivity is NOT ASSESSED. Exact physical-unit scaling tests are metamorphic verification, not threshold calibration.

## 7. Scientific and numerical robustness

| Condition | Expected behavior / observed result | Consequence / status |
|---|---|---|
| Missing/unknown units, nonuniform Z, conflicting spacing, invalid selection | Reject or retain unknown; tested conversion/selection contracts pass | Data integrity PARTIAL: C=0-only plane checking, unqualified origin/direction and actual registration remain limitations. |
| Invalid/overflowing/mixed IDs | No truncation/collision; reject export overflow and preserve internal source IDs | PASS for tested ranges. Huge dense match matrices and arbitrary int64 populations not stress-validated. |
| Large additive intensity baseline | Variance should not collapse; hand-computed variance preserved at +1e9 | PASS in tested range. Catastrophic overflow at extreme non-microscopy magnitudes NOT ASSESSED. |
| Empty/small/constant statistical inputs | No fabricated effects/p-values; fixed paths abstain/missing | PASS for new regressions; no general inferential adequacy claim. |
| Anisotropic/hollow/border geometry | Report proper units and declared raw/envelope meaning; flags preserve excluded/reviewed objects | PARTIAL. World registration and real boundary truncation accuracy unqualified. |
| Analytic radius-18 µm sphere, spacing (2,1,1) µm | Specification: <1% volume, <5% surface error. Measured **0.3562% volume**, **11.7347% area**, sphericity **0.89285** | Volume PASS for this case; surface **FAIL (P1 release blocker)**. Existing 15% test tolerance cannot override the 5% specification. |
| Low signal / saturation / stain mismatch / strong anisotropy | Correct uncertainty/QC behavior needs acquisition and annotation evidence | INSUFFICIENT EVIDENCE. uint16 range is not detector bit depth; plane thinning/noise/bleedthrough/penetration sensitivity not validated. |
| Cellpose downsample / batch / backend | Equivalent scientific output requires reference agreement | NOT ASSESSED beyond factor/shape adapter tests and one MPS crop. |

The sphere result is in [analytical-geometry.json](evidence/2026-09-10/analytical-geometry.json). No smoothing, formula, threshold or acceptance change was made to hide the failure. Discretization bias is a possible explanation, not permission to claim the requirement passed.

## 8. Engineering robustness

**Conclusion: PARTIAL.** Malformed axes/units, invalid label domains, shape/grid mismatches, incomplete output, overflow and save failures have tests. Output directories cannot silently overwrite nonempty results; incomplete sentinels survive injected failures. Native rendering depends on host graphics; browser interaction is unassessed. Full-volume memory/concurrency/cancellation behavior is not qualified. The Cellpose inference lock serializes the principal model run; auto-diameter model calls and external callers are not comprehensively concurrency tested. Multi-session stress, large-object-count EDT/matching memory, all corrupt TIFF codecs and interrupted historical legacy outputs remain NOT ASSESSED. Static QA and duplicate VTK metadata are unresolved, rather than waived.

## 9. Statistics, experimental units and leakage

**Conclusion: INSUFFICIENT EVIDENCE for biological/confirmatory inference.** Complete-unit filtering and global-across-batch replicate grouping are tested, and software preprocessing leakage was corrected. However object-level random-intercept LMM estimates differ from equal-well descriptive summaries; nested wells, paired donors and repeated conditions are not modeled by the shipped inference. The t(G−1) MixedLM contrast reference is a project heuristic with no type-I-error/coverage validation; omnibus remains asymptotic. BH is within each feature's contrast family, not global across all outputs. Missing/insufficient conditions can be omitted from results, so empty output is not evidence of no effect.

The demo deliberately shares replicate-size factors across conditions; its statistical outputs therefore cannot validate the independent-condition inferential design. It supplies no independent biological N, power or attainable-precision evidence. Classifier object splits remain unsuitable evidence of donor/well generalization, and cluster ANOVA/predictability on clustering-derived labels remains circular. No scientific statistical method was changed merely to obtain a pass.

## 10. Controlled truth and representative real-data results

The final installed-artifact demo runs 15 synthetic stacks / 90 designed objects: **90 detected and matched, 0 FP/FN, median Dice 1.0, minimum Dice 0.9996526572**, maximum relative voxelized-volume error **0.0006944444**. This is PASS for this seeded synthetic execution only. It uses demo-specific segmentation settings and project-generated masks/signals; perfect signal-state agreement is not independent viability validation. [Final synthetic outputs](evidence/2026-09-10/final-artifact/demo/synthetic_validation.json).

Five local TIFFs were freshly hashed/loaded. Complete calibration was available for two files, used without guessing in a real classical artifact run: **2 completed samples, 11 detected objects, 5 geometry-eligible**, no failed sample. Colon: 3 detected/0 eligible and high-anisotropy flag; spheroid: 8 detected/5 eligible. Viability remains indeterminate because channels/controls are unavailable. [Final real summary](evidence/2026-09-10/final-artifact/real/sample_summary.csv). IDs such as `file_identity_only_*` are runtime labels, not assertions of experimental independence.

These are runtime observations, not real-data segmentation accuracy or quantitative biological validation. No qualified organoid+cell+nucleus truth set, expert annotation protocol, independent assay, study design, held-out calibration or predeclared precision/N is present. Full multilevel biological E2E and inference validation are **INSUFFICIENT EVIDENCE**. MPS smoke only establishes execution on one crop.

## 11. Traceability matrix and lifecycle closure

| Frozen requirement / gates | Implementation | Oracle/data / final evidence | Status |
|---|---|---|---|
| Discovery/classification/specification (1–3) | Full source/config/entrypoint inventory; SPEC §13 | Current source, input state, explicit route/units/study limitations | PASS for discovery; intended-use suitability separately bounded. |
| Method/parameter design (4) | Existing algorithms retained; ALGORITHM_DECISIONS/PARAMETERS rewritten to actual behavior | Installed source + official method references + 113-row evidence registry | PARTIAL; project calibration missing. |
| Validation predeclaration (5) | Frozen A1–A10 and original specification criteria | `evidence/2026-09-10/PLAN.md`, VALIDATION_PLAN | PASS process; no retrospective relaxation. |
| A1–A2 identity/calibration (6–8) | microscopy_io, hierarchy, pairing | Exact units/time/IDs, fixtures and real file inventory | PARTIAL overall; tested contracts PASS, acquisition qualification missing. |
| A3 segmentation (6–8) | Cellpose/watershed adapters | Fake-model scale tests, installed method inspection, MPS crop | PARTIAL execution; accuracy INSUFFICIENT EVIDENCE. |
| A4 export/recovery (6–8) | OME/Parquet/sentinels/source hashes | Round trips, failure injection, singleton compatibility, artifact run | PASS tested contracts. |
| A5–A6 numerical stats/leakage (6–8) | inference/exploration | Domain oracles and fold-provenance tests | PARTIAL: software safeguards PASS; study-level inference unqualified. |
| A7–A8 QA/artifact (7) | Package/static assets, CLI/Web | Full suite, host render, lint, type/format, dependency/install checks | PARTIAL: explicit engineering failures/skips remain. |
| A9 scientific verification (8) | geometry/moments/topology/intensity | Analytic/hand-counted/metamorphic tests; area criterion probe | FAIL overall surface-accuracy requirement; other tested arithmetic PASS. |
| A10 intended-use validation (8) | Full research flow | Two unannotated runtime files, no qualified biological reference | INSUFFICIENT EVIDENCE. |
| Separation review / repair (9–10) | Complete final diff reviewed against spec/params/tests/outputs | Reproduction logs, R1–R12, current full regression | PARTIAL independence: Tier D LIMITED INDEPENDENCE; tested repairs PASS. |
| Reproducibility (11) | REPRODUCIBILITY + final snapshot | Patch/source/input/lock/model/artifact hashes and dependency inventory | PARTIAL: traceable local state; external data/assay/model transport unqualified. |
| Verdict (12) | This report | Failed core criterion and explicit missing evidence | NOT READY FOR THE SPECIFIED RESEARCH USE. |

## 12. Unresolved findings and release boundary

- **P1, FAIL:** the existing <5% absolute surface-accuracy requirement is not met. A scientifically justified method/criterion decision with independent revalidation is needed; this audit did not invent one.
- **P1 claim risk, INSUFFICIENT EVIDENCE:** broad accurate organoid/cell/nucleus segmentation, biological viability, and study-level inferential claims lack qualified data, calibration and reference standards. No positive claim is authorized by local execution counts.
- **P2, PARTIAL:** incomplete physical registration/acquisition provenance, inconsistent raw/envelope estimands across routes, requested-factor/rounded-grid downsampling uncertainty, assay saturation limits and large-volume/concurrency limits. These are explicit unsupported conditions, not silently repaired by choosing new scientific rules.
- **P3, FAIL:** static type/format debt and duplicate VTK installation metadata. A clean-environment rebuild and typed-interface cleanup remain engineering work; a large unrelated rewrite was excluded from this remediation.
- **P4, NOT ASSESSED:** external CI, browser interaction, vulnerability/secret scans and complete SBOM.

The repository is usable for the **tested controlled calculations and explicitly exploratory execution**, but that narrow observation does not satisfy the stated full research use. Source and tests remain local and uncommitted; no published release or certification was made. Independent biological/reference evidence and resolution of the failed scientific requirement are necessary before a positive research-readiness assessment.
