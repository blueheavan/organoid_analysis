# 3D morphology architecture review — 2026-09-13

## Scope and predeclared validation

S2 analytical maintenance: `features.geometry`, morphology adapters, surface
qualification metadata, exports and interpretation. The existing inferential
pipeline is S3 but its algorithms, exclusions and acceptance criteria are outside
this change. One common geometry engine remains responsible for numerical
measurements; object policies describe primary versus conditional reporting.

Before changes, source inspection identified lost domain fields and equivalent
diameter in multilevel morphology, lost domain fields in legacy cell morphology,
and stale marching-cubes wording in the UI. A runtime probe also confirmed that
the inscribed-radius EDT lacked exterior background: an all-foreground 3³ cube
at unit spacing reported rho_in=4.1231 and area=36.6899, versus rho_in=2 and
area=42.0576 for the same cube with a background margin. This changes stencil
selection and can falsely admit a small object to the numerical domain.

Validation criteria frozen for this repair:

1. Equal binary support, physical spacing and origin must give identical common
   geometry through every object adapter; raw and filled support are distinct
   estimands and retain explicit provenance. Cropping/background padding must
   not change scalar geometry (floating comparison rtol=1e-12).
2. Organoid/cell/nucleus tables and CSV/Parquet exports retain every surface
   domain variable from `geometry`, estimator/version and weight provenance.
3. An out-of-domain single nucleus retains hand-counted voxel volume,
   `(6V/pi)^(1/3)` diameter and moment axes. No surface flag removes an object.
4. Object type does not change qualification; sphericity and area/volume share
   the surface numerical gate. Smoothness, segmentation accuracy and biological
   validity cannot be inferred from that gate.
5. Keep Crofton weights, directions, stencil candidates and domain thresholds
   fixed. Reexecute the existing 188-case analytical harness against its stored
   values (existing rtol=1e-9), the qualified-domain tests and available frozen
   confirmation phantoms. Preserve original acceptance criteria: area <5%,
   volume <1%; surface and volume evidence remain separate.
6. Run locked Pixi lint, raw mypy, typing ratchet, tests, notebook check, local
   CI and science gate. Record failures without relaxing gates. Review the full
   diff after implementation and verification.

Known-truth tests use analytic voxel-box moments and independent volume/diameter
formulas. Cross-route equality and historical reproduction are regression
evidence, not independent biological validation. Audit separation is Tier D,
same agent/context: LIMITED INDEPENDENCE. Real segmentation and biological
validity remain NOT ASSESSED / INSUFFICIENT EVIDENCE without qualified references.

## Final architecture and policy

`geometry()` owns voxel-count volume, equivalent diameter, covariance axes,
axis ratios/elongation and conditional Crofton metrics. `measurement_policy.py`
describes the organoid, cell and nucleus estimands without numerical code.
Existing adapters map IDs/legacy names and add hierarchy, intensity or QC.
Multilevel joins and exporters preserve the shared surface metadata contract.

For all three object types, volume, equivalent diameter, axes, ratios and
elongation are primary; area, sphericity and area/volume are conditional.
Nuclear-only Web analysis requires no cell mask. Multilevel analysis also
retains orphan nuclei when parent tables are empty. Surface qualification never
drops an object. Existing upstream segmentation/pairing filters are unchanged.

Raw labels remain the multilevel/Web default. Legacy organoid/cell filled
envelopes remain the default for compatibility; `measurement_basis` exposes
that difference. A new optional `cell_geometry(fill_holes=False)` selects raw
shape support; legacy raw/envelope volume columns retain their meanings.
Comparisons across workflows must use the same support, spacing and origin.

## Defects and repairs

| Finding | Severity | Minimal repair | Evidence |
|---|---|---|---|
| Multilevel omitted equivalent diameter and all six engine surface-domain fields; legacy cells omitted domain metadata; Web omitted domain flags/budget/stencil | P2 | Shared `SURFACE_METADATA_COLUMNS`; forward values rather than recompute them | cross-route and CSV/Parquet regression tests |
| Cropped support changed EDT radius, stencil and surface area; a 12³ cube wrongly cleared rho>=10 | P2 | One-voxel background exterior for EDT; implementation revision `padded_edt_v1` | cube radius oracle; padded/cropped equality; unchanged frozen V&V |
| Object identity/smoothness and stale marching-cubes wording could overstate or misidentify measurements | P2 | Explicit numerical-only scope and separate NOT ASSESSED statuses; Crofton UI/report descriptions | policy and nuclear-only Streamlit AppTest |
| Two lint errors and 14 unratcheted mypy diagnostics existed in the Crofton module | P4 | Equivalent cache decorators and type annotations; narrow missing-SciPy-stub import annotations | lint and typing ratchet |

No direction, weight, scientific threshold, stencil candidate, acceptance rule,
segmentation algorithm or statistical procedure changed. The inscribed-radius
boundary repair intentionally changes conditional outputs on affected tight
crops. The saved sphere/holed-box comparison in `adapter-regression.json`
shows unchanged primary fields and unchanged sphere values; holed-box adapters
change only surface-related values. This is a bounded regression check, not a
claim that every possible input has been exhaustively compared.

## Output schema and interpretation

- Multilevel gains `equivalent_diameter_um`, `axis_ratio_minor_to_major`,
  `measurement_basis` and the complete surface metadata contract below.
- Legacy cells gain canonical selected-support `volume_um3`,
  `measurement_basis`, elongation/prolate/oblate ratios, area/volume and surface
  metadata. `cell_volume_um3` remains raw; `cell_envelope_volume_um3` remains filled.
- Web mask schema becomes 2.2, preserving legacy aliases and adding canonical
  diameter, axis ratios, area/volume and the complete surface metadata.
- Shared surface fields: `surface_rho_in`, `surface_anisotropy`,
  `surface_stencil_radius`, `surface_apriori_rel_bound`,
  `surface_in_qualified_domain`, `surface_domain_flags`, `surface_estimator`,
  `surface_method_version`, `surface_implementation_version`,
  `surface_weights_origin`, `surface_evidence`, `surface_implementation_evidence`,
  `surface_qualification_scope`, `sphericity_in_qualified_domain`,
  `surface_to_volume_in_qualified_domain`.
- JSON policy blocks record object/support, primary/conditional metrics, shared
  surface-gate applicability, and separate NOT ASSESSED statuses for smoothness
  scope, segmentation validation, biological validity and volume accuracy.
- Multilevel QC gains `surface_outside_qualified_domain` without excluding rows.

The method remains `crofton_minimax_sym_v3`; `padded_edt_v1` identifies the
boundary fix. Original estimator evidence and current implementation regression
evidence have separate locations. Non-packaged solved/cached weights do not
inherit the exact frozen-vector evidence. The numerical boolean does not
certify smoothness, segmentation, biological boundaries or volume accuracy.

## Scientific V&V and audit conclusions

The frozen 188-case production grid reproduced with zero mismatches at the
existing rtol=1e-9 (largest relative difference 4.6212e-10 from CSV precision).
Reexecuting the exact 96 previously qualified confirmation phantoms preserved
area, volume, rho_in and stencil selection. Worst absolute relative errors:
area 0.950779%, volume 0.761337%, sphericity 0.948692%; each meets its original
criterion on those phantoms only. This is reexecution of existing evidence,
not a new independent confirmation study. Original unrestricted SG-1 and SG-2
failures remain; no subset is substituted for the unrestricted gate.

| Requirement / audit dimension | Implementation and test | Status / limit |
|---|---|---|
| One geometry engine, workflow equality | `test_measurement_policy.py` | PASS for tested raw/filled supports and all three object types |
| Domain propagation / data integrity | `test_surface_metadata_exports.py`, existing bundle tests | PASS through joins and CSV/Parquet round trips |
| Small nuclear-only primary metrics | voxel-box oracle and nuclear-only Streamlit AppTest | PASS; no surface-based object deletion |
| Area/sphericity evidence inheritance | shared fields/policy tests | PASS for numerical-gate semantics; accuracy not inferred from object identity |
| Numerical correctness / robustness | crop invariant tests, qualified-domain tests, frozen 188+96 reexecution | PASS within exercised conditions; lattice-aligned/creased exceptions retained |
| Algorithm selection and scientific basis | unchanged minimax Crofton v3; voxel arithmetic/moments; scope explicit | SUPPORTED WITH LIMITATIONS / PARTIAL for actual biological applications |
| Parameter/threshold basis | original parameters untouched; existing frozen evidence reexecuted | PARTIAL: original domain restrictions still apply |
| Scientific robustness | small/outside-domain objects retained; scope not automatically certified | PASS for reporting contract; biological applicability INSUFFICIENT EVIDENCE |
| Engineering robustness | input checks retained; empty-parent workflow; exports/UI exercised | PASS for targeted paths; full-suite results below |
| Statistical validity / real segmentation / biological validity | no algorithm changes or new reference data | NOT ASSESSED in this task; representative validity INSUFFICIENT EVIDENCE |
| Audit separation | full code/docs/tests/output diff reviewed by implementing agent | Tier D, LIMITED INDEPENDENCE |

## Reproducibility and release boundary

Started from clean `b5b2c03213e92f4bd7fb79bcb8e5a23a4cd89e4e` on `main`.
Locked Pixi environment installed successfully. Current source and reference
SHA-256 manifests are in `evidence/2026-09-13-morphology-architecture/revalidation.json`.
The harness in that directory is rerunnable and leaves old records untouched.
No distribution artifact was created; artifact installation is NOT ASSESSED.

The canonical validation-record command refuses evidence-critical files that
differ from HEAD. Accordingly the old canonical record is rejected against this
dirty working tree, and the science gate remains NOT PASSED. No old evidence
hashes, PASS declarations, source commits or acceptance criteria were rewritten.
Engineering execution results and final release verdict follow below.

## Completed engineering checks — continuation on 2026-09-13

The interrupted task's results were recovered, the full diff was reviewed,
and local CI was rerun against the final implementation. All 62 source hashes
and six reference hashes in `revalidation.json` still match the working tree.
No scientific implementation changed during this continuation. Exact commands,
test identities, implementation/test/script hashes and log hashes are recorded
in [engineering-results.json](evidence/2026-09-13-morphology-architecture/engineering-results.json).

| Check | Actual result | Evidence and interpretation |
|---|---|---|
| `pixi install --locked` | PASS | Existing locked default environment installed successfully |
| Lint, executed by `pixi run ci` | PASS | [CI log](evidence/2026-09-13-morphology-architecture/checks/ci.log); source, tests and scripts checked |
| Full default tests, executed by `pixi run ci` | FAIL: 481 passed, 6 failed, 8 skipped, 8 deselected | 109.35 s; the six failures are the stale-record checks below; eight native-render tests are deselected by the default marker |
| `pixi run ci` | FAIL, exit 1 | Stops at the test failures; the separate checks below are not claimed as a successful combined CI run |
| `pixi run typecheck` | FAIL, exit 1 | [110 existing errors in 34 files](evidence/2026-09-13-morphology-architecture/checks/typecheck.log), 62 source files checked |
| `pixi run typecheck-ratchet` | PASS, exit 0 | [Baseline 110, current 110, new 0](evidence/2026-09-13-morphology-architecture/checks/typecheck-ratchet.log) |
| `pixi run check` | NOT APPLICABLE, exit 0 | [No notebooks found](evidence/2026-09-13-morphology-architecture/checks/check.log) |
| Focused geometry/policy/export/UI/qualified-domain tests | PASS: 43 tests | [Prior-phase log](evidence/2026-09-13-morphology-architecture/checks/focused-tests.log); all also pass in the resumed full run |
| `pixi run test-render` | PASS: 8 tests on host graphics | [Prior-phase log](evidence/2026-09-13-morphology-architecture/checks/render.log); not rerun in the restricted graphics environment |
| Frozen geometry reexecution | PASS for the stated regression scope | [188 production cases and 96 confirmation cases](evidence/2026-09-13-morphology-architecture/checks/numerical-revalidation.log), executed earlier in this same task; source/reference hashes verified unchanged |
| `pixi run science-gate` | FAIL, exit 1; 0/6 PASS | [Current gate log](evidence/2026-09-13-morphology-architecture/checks/science-gate.log) |

Two failures are in `test_scientific_validation_gate.py` and four are in
`test_evidence_integrity.py`. They assume the repository's sealed analytical
record describes the current implementation. Source hashes and estimator
metadata changed, so that assumption is false. Even the forged-PASS test
encounters this earlier estimator-identity rejection before its intended
reexecution check. The canonical writer's
[refusal](evidence/2026-09-13-morphology-architecture/checks/canonical-record-refusal.log)
is consistent with its contract: evidence-critical source must first match
a Git commit. These are unresolved evidence-lifecycle failures, not passing
tests or permission to weaken their assertions. No old record, current-record
pointer, integrity check, threshold or acceptance criterion was changed.

## Release assessment and remaining limitations

**NOT READY FOR THE SPECIFIED RESEARCH USE.** The shared-engine reporting
contract is implemented and its exercised numerical, export and UI paths pass.
The current working tree does not pass combined CI or the scientific gate.
A subsequent release needs a committed source snapshot and a newly executed,
sealed analytical record before repeating those checks. Reissuing that record
does not by itself fix the original unrestricted surface/volume failures or
supply independent segmentation, assay, study-design or calibration evidence.

The EDT correction intentionally changes surface area, sphericity and their
domain metadata for affected tight crops. Existing QC decisions that depend on
sphericity, and downstream summaries that consume those decisions, can also
change; such datasets should be reprocessed and compared. Primary geometry is
unchanged in the exercised regressions, but this is not an exhaustive claim
about every possible dataset. No inference algorithm or exclusion threshold
was changed. Real biological accuracy remains INSUFFICIENT EVIDENCE, and
smoothness, segmentation validation and volume accuracy retain their separate
NOT ASSESSED policy statuses.

The final review included all tracked changes, new policy/tests, and validation
artifacts. `git diff --check` passes. Review separation remains Tier D,
LIMITED INDEPENDENCE; this continuation does not constitute independent
validation. Hosted CI and distribution-artifact installation were NOT ASSESSED.
Changes remain uncommitted on `main`; no commit or push was performed.
