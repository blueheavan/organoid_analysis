# Validation Plan — Organoid Pipeline (Multilevel 3D Analysis)

2026-09-13 maintenance plan: the predeclared S2 criteria for shared morphology
reporting and cropped-mask EDT correction are in
[MORPHOLOGY_ARCHITECTURE_REVIEW.md](MORPHOLOGY_ARCHITECTURE_REVIEW.md).
All existing surface/volume acceptance criteria remain unchanged.

## Governing qualification framework — frozen 2026-09-13

This section and `SCIENTIFIC_VALIDATION_MASTER_PLAN.md` supersede conflicting
acceptance or sequencing language below. It is a protocol amendment, not a
validation result.

| Gate | Reference / evaluation unit | Required metrics and strata | Acceptance framework | Current status |
|---|---|---|---|---|
| SG-1a | exact analytical surface / phantom case | per-case error; `rho_in`, anisotropy, smooth-scope and exact-weight evidence | every in-domain confirmation case meets the frozen criterion; missing scope or non-evidence-bearing weights fail closed | `SUPPORTED WITH LIMITATIONS` as a scientific claim; canonical gate-record migration pending |
| SG-1b | exact analytical surface / phantom case | unrestricted worst, median, distribution and strata | characterization only; no verdict | characterized |
| SG-2a | exact analytical volume / phantom case | dense phase/orientation/resolution/anisotropy and topology strata | volume-specific domain first; every in-domain confirmation case meets the engineering criterion | `NOT QUALIFIED` |
| SG-2b | exact analytical volume / phantom case | unrestricted worst, median, distribution and strata | characterization only | characterized |
| SG-3A/B | independent blinded annotation / biological sample is inferential unit | detection, Dice/IoU, boundaries, signed/absolute measurement errors, upper tail, catastrophic failure; size/depth/domain/topology strata | population metrics characterize; qualification separately controls systematic bias, tail and failure rate with precision-calculated N | `INSUFFICIENT EVIDENCE` |
| SG-4A | known-rule synthetic tables / object | exact state and refusal behavior; versioned classifiable denominator | exact contract agreement | `PARTIAL` software evidence only |
| SG-4B | independent object reference plus orthogonal well-level assay / biological replicate | confusion matrix, weighted kappa, macro-F1, per-class sensitivity/precision; indeterminate rate; size/depth/radial/batch/stain/acquisition strata | study-specific values and interval precision frozen before outcome review; correlation is characterization | `NOT ASSESSED` |
| SG-5 | simulation truth / biological replicate cluster | branch-specific type-I error, CI coverage and estimand recovery over ICC, replicate N and imbalance | method frozen before implementation; only cells meeting frozen coverage tolerance support inferential CI | `INSUFFICIENT EVIDENCE` |
| SG-6A | traceable ratio standard / imaging session | X:Y, Z:XY, uncertainty, session and field variation | stratum-specific values frozen from domain-misclassification tolerance | `NOT ASSESSED` |
| SG-6B | traceable absolute/registration standards / imaging session | X/Y/Z bias and uncertainty, depth/field/session variation, channel offset | stratum-specific absolute-unit error budget | `NOT ASSESSED` |

SG-3A and SG-3B are separate studies. SG-6A precedes final SG-3 domain
stratification; SG-6B can run in parallel and is required before absolute-unit
claims. The volume `<1%` rule is an engineering estimator target, not a real
segmentation or biological requirement.

For SG-3 through SG-6, a point estimate cannot PASS without a predeclared
uncertainty/precision requirement at the correct independent-unit level. If N
cannot resolve that requirement, status is `INSUFFICIENT EVIDENCE`.

Version: 1.0.0
Date: 2026-09-02
Scope: Scientific V&V plan for the `analyze-3d` multilevel 3D organoid analysis (hierarchy, morphology, topology, spatial, QC). The classical morphology/viability validation is documented in `docs/VALIDATION.md`.

Historical scope: numerical multilevel checks only. Its recorded PASS labels are not evidence for the current repository; fresh results are in VALIDATION_REPORT.md.

This plan predeclares, for each critical requirement, the reference standard, evaluation unit, dataset, independence level, minimum N, metric, and acceptance criteria. Status is filled in during execution; see `docs/VALIDATION_REPORT.md`.

---

## VR-1 Data-contract and identity integrity

- **Requirement:** Registered organoid/cell/nucleus label masks are accepted only when they share a common physical grid and valid integer IDs; identity is never silently relabeled or reshaped.
- **Reference standard:** Defined input contract (`validate_labels`, `validate_inputs`).
- **Evaluation unit:** per input volume.
- **Dataset:** synthetic + any real registered volumes.
- **Independence:** N/A (definitional).
- **Minimum N:** 1.
- **Metric:** error raised on shape mismatch, non-3D, negative or boolean labels, non-finite intensity.
- **Acceptance:** `PASS` if contract is enforced and tests cover each violation type.
- **Status:** `PASS` — enforced by `validation.py`; tests in `tests/workflows/test_multilevel_measurement_workflow.py`.

## VR-2 Hierarchy correctness (maximum-overlap parent assignment)

- **Requirement:** Each cell is assigned to its max-overlap organoid; each nucleus to its max-overlap cell; nucleus organoid membership inherits through its cell; direct-overlap audit is retained.
- **Reference standard:** Hand-computed overlap counts on a controlled phantom with known geometry.
- **Evaluation unit:** per object.
- **Dataset:** `synthetic_labels()` phantom (`tests/workflows/test_multilevel_measurement_workflow.py`).
- **Independence level:** independent author computes expected assignments by hand.
- **Minimum evaluable N:** 5 cells, 5 nuclei spanning multiple parents.
- **Metric:** fraction of objects with correct parent; QC flags `crosses_multiple_parents`, `parent_assignment_failed`, `low_parent_overlap`, `direct_organoid_parent_mismatch` set correctly.
- **Acceptance criterion:** 100% correct assignments; all flagged states correct.
- **Uncertainty/tolerance:** deterministic; no stochastic tolerance.
- **Status:** `PASS`.

## VR-3 Physical volume and centroid

- **Requirement:** `volume_um3` = voxel_count × voxel volume; centroid in physical µm.
- **Reference standard:** analytical voxel volume formula, independent centroid computation.
- **Dataset:** synthetic phantom + analytical sphere.
- **Metric:** volume exact (integer × spacing product); centroid agreement.
- **Acceptance:** volume must equal exact voxel volume; centroid within one voxel of independent computation.
- **Status:** `PASS`.

## VR-4 Surface area and sphericity (numerical correctness)

- **Requirement:** surface area from native-spacing level-0.5 marching cubes; sphericity = `π^(1/3) (6V)^(2/3)/A`; no silent clipping to [0,1].
- **Reference standard:** independent analytical sphere phantom for known bias; correctness of formula.
- **Dataset:** 15-voxel-radius sphere (spacing 1).
- **Metric:** measured vs analytical V and A ratios; sphericity value; behavior when sphericity >1.05.
- **Acceptance:** document discretization bias (measured V/A ratios); confirm sphericity >1.05 is flagged, not clipped. Do NOT require bias-free area (bias is expected and documented).
- **Status:** `PARTIAL — formula/numerical correctness PASS on phantom; absolute surface accuracy on real organoids INSUFFICIENT EVIDENCE (requires independent surface reference).`

## VR-5 Contact topology

- **Requirement:** Cell adjacency = voxel-face contact only; anisotropic physical contact area correct.
- **Reference standard:** hand-computed face counts and face areas on a controlled phantom.
- **Dataset:** synthetic_labels() phantom (101/102 X-contact).
- **Metric:** edge existence, contact_area_um2, degree, total_contact_area_um2.
- **Acceptance:** X-normal 5×5 shared plane → contact_area = 25 × (Z×Y); degree correct.
- **Status:** `PASS`.

## VR-6 Spatial features (radial position, surface depth)

- **Requirement:** centroid-to-organoid-centroid distance in µm; normalized radial position using equivalent-sphere radius; distance-to-surface from EDT depth.
- **Reference standard:** independent computation on a controlled phantom; defined operational semantics.
- **Dataset:** synthetic phantom.
- **Metric:** distance values; normalized radial positions in plausible range; depth ≥0.
- **Acceptance:** deterministic values match independent computation; depth semantics documented (0 at/outside surface).
- **Status:** `PARTIAL — computational correctness PASS on phantom; biological meaning of radial bins NOT ASSESSED.`
- **Known edge condition:** a cell centroid outside its organoid mask reports `distance_to_organoid_surface_um = 0` (on/outside surface). This is a documented operational choice, not a failure.

## VR-7 QC flags

- **Requirement:** Formal QC flags generated as a feature-database table (`qc_flags.parquet`), including MAD volume outliers, border touch, too-small, parent assignment, anucleate/multinucleated, direct-organoid mismatch.
- **Reference standard:** independent implementation of robust z-score and flag logic.
- **Dataset:** synthetic phantom with engineered border/outlier/anucleate/multinucleated/mismatch cases.
- **Metric:** exact flag set per object.
- **Acceptance:** all engineered flags produced; MAD==0 degenerate case handled.
- **Status:** `PASS`.

## VR-8 Representative real-data end-to-end behavior

- **Requirement:** The pipeline runs end-to-end on representative registered real organoid/cell/nucleus label volumes and produces physically plausible feature tables.
- **Reference standard:** expert review of feature distributions; no ground-truth claims.
- **Dataset:** real registered volumes (when available/permitted).
- **Independence:** domain expertise for plausibility review.
- **Minimum N:** TBD by available dataset.
- **Metric:** successful export; QC flag rates; physically plausible volume/area/ranges.
- **Acceptance:** runs without error; outputs nonempty features; plausible ranges.
- **Status:** `INSUFFICIENT EVIDENCE — depends on available real registered data; not fully assessed in this pass.`

## VR-9 Reproducibility / determinism

- **Requirement:** The measurement pipeline is deterministic given identical inputs and config.
- **Reference standard:** identical re-run reproducibility.
- **Dataset:** synthetic phantom.
- **Metric:** identical feature tables across two runs.
- **Acceptance:** byte-identical outputs.
- **Status:** `PASS — pipeline is deterministic (no RNG in measurement).`

## VR-10 Failure modes and out-of-domain input

- **Requirement:** Missing spacing, conflicting metadata, shape mismatch, non-finite intensity, empty labels are handled (error or explicit empty result), never a silently wrong number.
- **Reference standard:** defined failure behavior.
- **Dataset:** synthetic malformed inputs.
- **Metric:** correct error/exception or empty-DataFrame path.
- **Acceptance:** each case behaves as defined.
- **Status:** `PASS for tested cases (see VALIDATION_REPORT.md for full list).`

## VR-11 TIFF axis/spacing reader behavior unification (S1, D12)

- **Requirement:** Both TIFF-reading paths (`tiff_contract.py` classical-CLI / `zstack_reader.py` + `metadata.py` Streamlit) enforce the same consequential reader behavior: an unrecognized axis is rejected regardless of size (never silently squeezed), and a nonuniform Z-plane grid is rejected rather than silently summarized as one Z spacing. The two unit-conversion tables and spacing-tolerance constants are a single shared definition (`metadata.to_um`), so axis/spacing semantics are identical across readers.
- **Reference standard:** defined input contract per docs/ALGORITHM_DECISIONS.md D12 (strictest behavior of the two readers wins); unit strings "µm"/"μ"/"um"/"micrometer" accepted and converted to µm.
- **Evaluation unit:** per uploaded/parsed TIFF volume.
- **Dataset:** synthetic TIFFs exercising ambiguous singleton axis and nonuniform-Z metadata; 5 real sample images identified in `docs/validation_data_manifest.csv` (by filename + SHA256), re-loaded after the change.
- **Independence:** N/A (definitional contract); re-load of real sample volumes provides reader-path smoke coverage.
- **Minimum evaluable N:** 1 per failure type; 5 real samples.
- **Metric:** correct error (or acceptance) rather than a silently wrong axis/unit interpretation; µm conversion correctness.
- **Acceptance criterion:** `PASS` if the regression tests below pass and real samples reload without regression.
- **Aggregated/naming:** `tiff_contract.py` imports `metadata.to_um` rather than redefining a table; the union of both prior tables is retained.
- **Reproducibility caveat (P2-7 audit finding):** `data/images/*` is gitignored and was never tracked in this repository (`git ls-files data/images/` returns only `.gitkeep`). The 5 real sample files exist only on the machine this validation ran on; a fresh clone of this repository cannot reproduce this re-load step on its own. `docs/validation_data_manifest.csv` records each file's exact identity (filename + SHA256 + size) so the specific files this PASS refers to are at least identifiable and re-obtainable outside the repository, but the re-load itself is not repository-reproducible.
- **Status:** `PASS — tests: test_io.py::test_ambiguous_singleton_axis_is_rejected_not_silently_squeezed, test_io.py::test_nonuniform_z_positions_are_rejected, test_3d_preview.py::SpacingMetadataTests::test_nonuniform_z_positions_are_rejected / test_uniform_z_positions_are_accepted (all repository-reproducible). The "5 real sample images re-loaded with no regression" claim is PARTIAL: true on the machine that ran it, identity recorded in docs/validation_data_manifest.csv, but not independently reproducible from a repository clone alone.`

---

## Predeclared acceptance for this pass (multilevel 3D)

Based on the controlled phantoms and the deterministic, well-scoped nature of the measurement, we predeclare:

- VR-1, 2, 3, 5, 7, 9: **PASS** via controlled phantoms.
- VR-4, VR-6: **PARTIAL** — numerical/computational correctness PASS; absolute accuracy / biological meaning INSUFFICIENT EVIDENCE.
- VR-8: **INSUFFICIENT EVIDENCE** — requires representative real registered data.
- VR-11 (S1, D12): **PASS** for the repository-reproducible regression tests; **PARTIAL** for the real-sample-reload claim specifically (identity recorded in `docs/validation_data_manifest.csv`, but the files are not tracked in git and the re-load is not reproducible from a clone alone).

## Independence level

Audit separation for this pass: Tier D (same agent and context as implementation) → final report records `LIMITED INDEPENDENCE`. Oracle independence: controlled phantoms use hand/independent expected values for hierarchy/topology/QC; morphology uses independent analytical sphere constants.

## Traceability

- Plan mapping to code: see `docs/ALGORITHM_DECISIONS.md` traceability table.
- Tests: `tests/workflows/test_multilevel_measurement_workflow.py`.
- Result outputs: `export_results` → Parquet/JSON hierarchy described in README.


## 2026-09-10 full-workflow audit acceptance contract

The [frozen plan](evidence/2026-09-10/PLAN.md) was written before new remediation/final QA and preserves the input-worktree audit criteria. It extends the earlier multilevel scope across actual I/O, preview, segmentation, quantification, statistics and export. New defect probes exercise the same mathematical/domain invariants, not tuned biological thresholds.

The §VR-4 instruction to document bias does **not** waive SCIENTIFIC_SPEC §9's <5% area requirement. Both are evaluated and the conflict reported; no acceptance criterion was relaxed. Likewise test_geometry.py's 15% surface tolerance is a regression criterion, not release accuracy acceptance.

| Requirement | Oracle / evaluation unit / independence | Dataset and necessary conditions | Acceptance / uncertainty / failure |
|---|---|---|---|
| A1–A2 Input identity/calibration | Unit conversion and distinct-valued arrays; deterministic per series/time/channel, not biological N | Metadata fixtures, invalid units, nonuniform/reversed Z, selected time/series, mixed IDs | Exact identity; stated 1e-12 conversion tolerance; unsupported semantic coercion rejects. |
| A3 Segmentation adapter | Captured model inputs and installed implementation semantics; mock is not accuracy oracle | Full/half/quarter resolution, invalid model/input/config, model output masks | Valid call scale and shape; invalid inputs reject. Rounded resampling equivalence and actual segmentation accuracy require further evidence. |
| A4 Export and recovery | Independent TIFF/Parquet reread and fault injection | Large uint32 IDs, overflow, RGB-prone shapes, singleton-Z saved-run compatibility, partial writes | Exact supported label/spacing recovery; no complete output on failed save/export. |
| A5–A6 Statistics | Definitional domain (positive logs, finite values, nonzero variance); training-fold provenance | Zero/nonfinite measurements, missing feature/fold, too few samples, constant groups | No invented effects/p-values; no train/test preprocessing leakage. This does not establish inferential calibration or independence. |
| A7–A8 Engineering | Actual distribution install and runtime; native/browser tests separate | Locked Pixi environment, wheel, CLI, unit/integration suite, lint, mypy, formatter, notebook discovery, GUI startup | Every command status recorded; skips and static failures cannot become PASS. |
| A9 Geometry/topology | Analytic voxel counts/moments; continuous sphere; deterministic invariants | Sphere, ellipsoid, anisotropic box, hollow/border masks, mixed label types, stable intensity moments | Existing exact arithmetic/scale criteria; §9 volume <1%, surface <5%. A failed surface criterion blocks the broad quantitative claim. |
| A10 Real research validity | Expert/orthogonal assay reference independent of development; independent biological unit | Representative normal/low-signal/saturated/aniso/dense/hollow acquisitions; actual donor/well/batch design; qualified labels/assay | Required minimum N and precision are NOT ASSESSED because intended experiment/reference uncertainty is unavailable. Do not declare PASS from object count or one crop. |

Sample size for deterministic mathematical fixtures is the number of explicitly exercised cases, not an estimate of population accuracy; statistical confidence intervals are inapplicable to exact arithmetic. For real biological performance the missing N/precision/reference qualification leads to INSUFFICIENT EVIDENCE. Any future calibration or acceptance change must be predeclared and evaluated on held-out independent units.

## 2026-09-11 scoped workflow optimization

[O1–O7](evidence/2026-09-11/PLAN.md) predeclare raw/envelope geometry, precision,
QC, effective resampling scales, preview budget, feature exports and interface
verification. They do not relax A9's surface accuracy requirement or replace
A10's qualified biological evidence. Historical A3 tests asserting requested
scale are updated to the independently derived center-grid scale; interpolation
coordinates, fractional intensities and unsupported rectangular grids are tested
explicitly. The [scoped report](WORKFLOW_OPTIMIZATION_2026-09-11.md) separates
these new contract checks from the broader unresolved research-readiness items.

## 2026-09-11 measurement validation and release readiness

The [surface V&V plan](evidence/2026-09-11-measurement-vv/SURFACE_VV_PLAN.md) was
frozen before any candidate estimator result. It keeps the §9 criterion
unchanged and evaluates it **per case**, so a stratum passes only if its
worst case is below 5%. It predeclares candidate estimators, parameter
provenance (only σ = max spacing may be selected for E1), the ρ, shape and
anisotropy strata, and a replacement rule that requires an unseen confirmation
grid. The confirmation grid was not run because no candidate passed the
development grid.

The [segmentation validation protocol](evidence/2026-09-11-measurement-vv/SEGMENTATION_VALIDATION_PROTOCOL.md)
defines the reference-standard qualification, the sample-level independent
unit, and the detection, mask and downstream-measurement metrics that A10/SG-3
require. It reports no result, and SG-3 remains INSUFFICIENT EVIDENCE.

Two gates are now separate. The regression tolerance in `test_geometry.py`
(15%) and the new estimator version-lock test belong to the engineering gate
(`pixi run regression`; `pixi run ci` adds the mypy debt ratchet). The scientific acceptance items are
verified against a sealed evidence record (docs/VALIDATION_RECORDS.md) and reported by
`pixi run science-gate`.
