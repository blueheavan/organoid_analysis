# Final separation review (Tier D — LIMITED INDEPENDENCE)

Reconstructed the end-to-end paths from current source and cross-checked the full
tracked diff and six incoming regression modules plus the audit additions.
No separate agent, independent reviewer, or hosted CI was used.

- SPEC versus actual routes: separate classical/Cellpose/cell-pair/multilevel
  estimands, metadata readers, local voxel coordinates and study units are explicit.
- Algorithms versus installed source: Lewiner default, Cellpose normalization,
  ignored 3D flow threshold, dynamics size/iteration defaults, EDT semantics and
  train-fold pipelines verified; mathematical formulas are not biological proof.
- Parameters versus outcomes: no numerical default was calibrated or changed to
  improve validation. New input-domain checks and missing-result handling are
  justified by identity/type/arithmetic requirements. Registry evidence levels
  do not equate method defaults with suitability.
- Acceptance: original <5% sphere surface error fails; current tests' 15% allowance
  is not promoted to an accuracy criterion. Original specification is retained.
- Repairs: sparse IDs compact only internally with source mapping, mixed integer
  pairs remain exact, export is calibrated grayscale ZYX, singleton Z recovery
  is constrained by saved shape, NaN effects/p-values remain missing, and the
  input-worktree fixes for units/time/variance/leakage/provenance were revalidated.
- Regression oracles: hand-counted labels, exact voxel units, analytical variance,
  independent file decoding and training-fold provenance. Mock outputs verify
  adapter invocation, never accuracy. Reproduced failure logs retained.
- Diff scope: no unrelated architecture changes, threshold tuning, validation
  tolerance changes, test disabling, type-ignore expansion or broad formatting.
  Private provenance file-list constants were replaced, not scientific behavior.
- Runtime: native render passed only on host; browser cases remain skipped.
  Real model crop and unannotated classical E2E establish execution only.
  Final wheel is installed outside checkout and tested against final source.
- Release: unresolved scientific accuracy and evidence blockers forbid a positive
  research-readiness claim. Static/dependency failures remain separately visible.

After final executable changes, full non-render regression, lint and installed
artifact E2E were run. Type-only cleanup reduced current mypy errors from 199 to
193 (fresh HEAD 194); remaining added diagnostics are missing sklearn stubs.
The last source edit was a docstring correcting an unsupported MixedLM claim;
final full tests and artifact build also include it. No further scientific code
changes followed final QA. Artifact script's manifest path was made persistent
in the evidence directory; it points to the same copied manifest bytes.
