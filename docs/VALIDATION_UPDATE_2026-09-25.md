# Scientific audit and repair update — 2026-09-25

## Scope and scientific impact

Reviewed the current 3D measurement path (`geometry`, voxel-spacing validation,
classical workflow and summaries), the reported condition-inference path, the
existing scientific specification, validation plan and evidence gate. Source
baseline: `a34630178a771d94e7b87f7e7355fa3ebf80cb6a` on `main`, clean
before this repair. The two repairs below are S2 analytical input-contract
changes: invalid input could change reported physical measurements. No method,
weight, threshold, acceptance criterion, or intended-use domain was changed.

## Confirmed defects and minimal repairs

| Finding | Reproduction and scientific consequence | Repair | Severity |
|---|---|---|---|
| Binary-mask contract | A 5×5×5 `uint8` mask with 27 foreground voxels valued `2` produced `segmented_voxels=54` and `segmented_volume_um3=54` at 1 µm spacing, while `volume_um3=27`. This yields internally inconsistent measurements if a caller supplies a label-valued array to `geometry()`. Production instance call sites supply boolean masks, so an affected saved run was not established. | Accept boolean and numeric 0/1 masks; reject other values before measurement. | P2 |
| Physical-spacing type | `validate_voxel_spacing_xyz((True, 1, 1))` and `validate_voxel_spacing_xyz(("1", 1, 1))` both returned `(1, 1, 1)` in µm. Direct `geometry()` likewise accepted these values. A malformed flag/text value could thus become a plausible physical length. | Reject flags and text at both measurement entry points; retain real scalar values including NumPy real scalars. | P2 |

The predeclared repair checks were: a 27-voxel 0/1 block with Z spacing 2 µm
must have `segmented_voxels=27` and both volumes 54 µm³; boolean and 0/1
representations must agree; a value-2 mask must raise; boolean/text spacing
must raise while finite real spacing remains accepted. These are input and
arithmetic checks, not evidence of accuracy on biological specimens. The
existing analytical geometry, segmentation, assay and inference acceptance
criteria remain unchanged.

## Verification and evidence boundaries

| Check | Status |
|---|---|
| Focused geometry and spacing regression | PASS: 30 tests in `test_geometry.py` and `test_voxel_spacing.py`, including the six new input-contract cases. |
| Classical end-to-end measurement and condition report | PASS: `test_pipeline_writes_pairwise_contrasts_end_to_end` on controlled synthetic stacks with the existing valid-mask path. This checks execution and output plumbing, not inferential calibration. |
| Full default engineering suite | INSUFFICIENT EVIDENCE: interrupted after 123 passing cases when the next sparse-ID mask-feature case spent several minutes inside SciPy HiGHS solving previously uncached Crofton weights for spacing `(2, 3, 4)`. That case, and a second affected-workflow selection, were also interrupted during the same solve. No full-suite PASS is claimed. |
| Ruff lint | PASS: `ruff check src tests scripts`. |
| Typecheck ratchet | PASS against existing debt baseline: 110 current, 110 baseline, 0 new. The 110 diagnostics remain unresolved. |
| Scientific gate | FAIL as designed: 0/9 qualification items PASS. The current analytical record is rejected because `features.py` and `voxel_spacing.py` differ from its sealed source hashes. Its previous SG-1/SG-2 values cannot be attributed to this source snapshot. |
| Current-record acceptance test | FAIL as expected for this dirty source: `test_current_analytical_evidence_is_accepted_and_backs_the_items` expects an accepted record (`INSUFFICIENT EVIDENCE` for SG-1a/SG-2a), but the gate correctly returns `FAIL` for the stale record. The test was not relaxed. |
| Representative real-data validation for these repairs | NOT ASSESSED. No qualified independent segmentation/physical calibration reference was introduced. |

`pixi run test` could not prepare its cache under the sandbox and a temporary
cache retry could not reach PyPI. Tests above used the existing Pixi Python
directly with `PYTHONPATH=src`; this verifies the current checkout against the
installed environment, not a fresh locked bootstrap. The uncached Crofton LP
is a pre-existing runtime limit for unusual spacing ratios, not evidence that
the binary-mask change altered the estimator. The interruption leaves its
functional result unassessed in this run.

The current analytical record cannot be regenerated while evidence-critical
source differs from HEAD: `record_analytical_geometry_validation.py` checks
that condition before execution. A future source commit, fresh frozen V&V run
and new sealed record are required for an attributable scientific gate result.
The old record remains unchanged. This repair does not establish that SG-1 or
SG-2 meets its scientific accuracy criterion.

## Bounded conclusion

Algorithm selection and parameter basis are unchanged; their project-data
suitability remains `INSUFFICIENT EVIDENCE` where the governing validation plan
requires independent annotations, assay references, spacing calibration and
study-level statistical validation. Invalid-input behavior is covered by
regression checks; biological validity, inferential calibration and real-data
robustness remain unqualified. Overall verdict: **NOT READY FOR THE SPECIFIED
RESEARCH USE**. Review separation is Tier D, `LIMITED INDEPENDENCE`: the same
agent implemented and reviewed this patch; the 27-voxel arithmetic oracle is
independent of the implementation but is only a controlled example.
