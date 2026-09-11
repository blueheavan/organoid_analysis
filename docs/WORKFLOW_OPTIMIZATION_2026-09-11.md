# Workflow optimization and scoped validation — 2026-09-11

**Overall verdict: NOT READY FOR THE SPECIFIED RESEARCH USE.**

The repaired S2 contracts pass their controlled checks. The broader intended
publication use remains blocked by the existing surface-accuracy failure and
missing qualified biological/inferential evidence. This is a scoped update to
[the full audit](VALIDATION_REPORT.md), not a new claim of biological validation.

## Changes and scientific consequences

| Stage | Confirmed defect or missing workflow capability | Implemented behavior |
|---|---|---|
| Z-stack preview | Requested scaling misrepresented actual center-coordinate extent; the 0.1 factor floor could exceed the array budget | Use actual resized dimensions and separate X/Y spacings; honor the array budget, retain Z and raw arrays; oversized previews report an explicit limitation while measurement remains available. |
| Cellpose adapter | Requested scale differed from SciPy's center-grid scale; integer output rounded interpolated intensity | Use `(M-1)/(N-1)` for diameter and anisotropy, float32 interpolated intensity and unchanged nearest-neighbor label restoration; save the working shape/spacing/scale in provenance. Unequal effective XY scales reject before inference, with full resolution available. |
| Feature geometry | Web filled volume/axes but used raw centroids/solidity; four-decimal rounding could turn a positive volume into zero | Feature schema 2.0 defaults to raw nucleus/cell labels; explicit envelope mode uses one support for every geometric feature. Keep full precision, raw voxel count/volume and filled-voxel count. |
| QC and analysis | Border/fragment status absent from Web feature table; solidity could be infinite for degenerate support | Retain every object with border, fragment, enclosed-instance and non-estimable-solidity review flags. Undefined solidity remains missing. Descriptive summaries explicitly include flagged objects. |
| Result delivery | Web measured only nuclei and did not provide a feature download with measurement lineage | Select nucleus/cell and raw/envelope mode; selection refreshes cached features, overlay and CSV/JSON ZIP. Record spacing origins, mask/CSV/source hashes, model/run provenance when available and explicit unknown calibration. |

The classically filled organoid workflow, multilevel raw measurements, statistical
models, biological thresholds and <5% surface-accuracy criterion were retained.
No automatic object exclusion, biological-replicate assignment or scientific
performance claim was added. Historical Web values must be compared using an
explicit geometry definition; new CSVs identify it in `measurement_basis`.

## Evidence and traceability

The [O1–O7 plan](evidence/2026-09-11/PLAN.md) was frozen before implementation.
The source began at clean `657204d8d8be2e547063272fd139161f59dc62be`.

| Requirement | Oracle and evidence | Result |
|---|---|---|
| O1 raw/envelope semantics | Hand-counted asymmetric hollow box, physical centroid and raw multilevel comparison; `test_mask_feature_semantics.py` | PASS for tested volume/centroid/support and cross-route consistency. |
| O2 precision/schema | Positive 1e-6 µm³ single-voxel volume, 17-digit CSV round trips, sparse IDs and empty headers | PASS; compute-time rounding removed. |
| O3 QC | Hand-constructed border, face-disconnected fragments, degenerate hull | PASS for exercised flags; operational review rules are not biological classifications. |
| O4 resampling | Even/odd square and rectangular grids, coordinate ramps, captured model calls, installed SciPy 1.17.1 source | PASS for coordinate/adapter contracts; actual model equivalence NOT ASSESSED. |
| O5 preview integrity | Forced 2000/150-byte budgets, large sparse IDs, impossible budget, mismatched grids and immutable arrays | PASS for tested array budgets/coordinate identity. HTML/mesh/browser memory is separately unqualified. |
| O6 export/interface | ZIP/CSV/JSON reread, mask/calibration/basis mismatch rejection, Streamlit AppTest layer and envelope changes | PASS; download content follows the selected mask and geometry. |
| O7 engineering | [Full suite](evidence/2026-09-11/full-tests.log): **336 passed, 8 skipped, 8 deselected**, 899 warnings | PASS for executed tests; [JUnit](evidence/2026-09-11/pytest.xml). Skips are not passes. |
| Static checks | [Lint](evidence/2026-09-11/lint-notebooks.log) passes; notebook check finds no notebooks; [mypy](evidence/2026-09-11/typecheck.log) reports **134 errors in 36 files**, versus baseline 135 in 35 | Lint PASS; notebooks NOT APPLICABLE; typecheck FAIL. Three touched-function annotations repaired; two additional missing-stub diagnostics appear in new imports. No check settings were relaxed. |
| Real TIFF reader/preview | Human-Colon-Organoids-C1.tif, metadata calibrated, `(16,512,512)` → `(16,107,107)` within 1 MB array budget | PASS for runtime, unchanged raw bytes and preserved center extent; [record](evidence/2026-09-11/real-preview.json). No accuracy assessment. |
| Saved-run recovery/features | Existing `(32,128,144)` nucleus labels: 6 objects, raw-volume sum **62622 µm³**, equal to direct voxel-count volume | PASS for recovery/measurement/export; [record](evidence/2026-09-11/saved-run.json) and [example bundle](evidence/2026-09-11/saved-run-nucleus-features.zip). Original biological/synthetic domain was not independently qualified. |
| Existing surface criterion | Radius-18 µm sphere, ZYX spacing `(2,1,1)` µm: **0.356234% volume error; 11.734744% surface error** | Volume PASS for this phantom; surface FAIL against unchanged <5% criterion; [fresh probe](evidence/2026-09-11/analytical-geometry.json). |

## Scientific audit conclusions

- **Algorithm scientific basis — PARTIAL.** Coordinate arithmetic and voxel
  measurement semantics have analytical support. Raw/envelope choice is context
  dependent. Marching-cubes surface accuracy still fails the declared criterion.
  Cellpose and actual biological target suitability remain INSUFFICIENT EVIDENCE.
- **Parameter/threshold scientific basis — PARTIAL.** This change derives scales
  from actual grids and makes geometry explicit. It introduces no calibrated
  biological cutoff. Existing model/QC/assay thresholds still lack project
  calibration and held-out sensitivity evidence.
- **Scientific robustness — PARTIAL.** Ambiguous scalar XY resampling rejects;
  default calibration, raw/envelope semantics and QC are explicit. Image-grid
  agreement does not prove registration, acquisition validity or target identity.
- **Numerical robustness — PARTIAL.** Controlled counts, centroids, scales,
  fractional interpolation, sparse IDs, precision and missing solidity pass.
  Surface bias and extreme dynamic-range/large-population behavior remain open.
- **Engineering robustness — PARTIAL.** Selected-layer cache invalidation,
  aligned preview arrays, budget exhaustion and export-metadata mismatch handling
  pass. Multi-session concurrency, full-volume memory, actual browser GPU
  interaction and native rendering were not newly assessed.
- **Statistics/interpretation — INSUFFICIENT EVIDENCE for publication inference.**
  New exports describe objects from one field. Existing S3 model families and
  their study-unit/coverage limitations are unchanged. Object rows are not
  asserted to be independent biological replicates.

Review separation is **Tier D — LIMITED INDEPENDENCE**: the same agent performed
implementation and final source/spec/output review. Analytical oracles are
separate from the implementation where stated; no independent biological
validation, new Cellpose inference accuracy run or hosted CI is claimed.

## Reproducibility and remaining work

[The snapshot](evidence/2026-09-11/snapshot.json) records the commit, dirty state,
tracked patch hash, relevant untracked hashes, source manifest, lock and package
versions. Reproduce the optional local runtime probes with
`pixi run python docs/evidence/2026-09-11/capture_evidence.py`; local fixtures are
identified by hashes and may be absent in a fresh clone. No software distribution
was published, and no commit or push was performed.

Remaining **P1** findings: failed surface-accuracy requirement; insufficient
independent segmentation annotations, assay validation and study-design evidence
for broad publication claims. **P2** limitations include acquisition registration,
resampling sensitivity, default calibration, and large-volume/concurrency
qualification. **P3** type/stub debt remains. Native/browser rendering, new model
backend equivalence and distribution installation are NOT ASSESSED in this
scoped change. The supported result is improved, tested measurement contracts and
traceable exploratory outputs, not publication-readiness certification.
