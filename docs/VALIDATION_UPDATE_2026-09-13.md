# Scientific validation plan update — 2026-09-13

## Scope and impact

This S3 implementation pass records the project owner's scientific-validation
plan and applies the contract changes that can be completed without new
experimental data. It changes claim semantics, validation sequencing,
measurement provenance and viability aggregation behavior; it creates no study
result.

## Changes

- froze Option C qualification/characterization gate architecture;
- retained outer-envelope volume as primary and classified analytical `<1%` as
  a conservative engineering target rather than biological necessity;
- classified analytical volume qualification as `NOT QUALIFIED` pending a
  volume-specific domain and untouched confirmation record;
- split SG-3 by brightfield/membrane modality, SG-4 by analytical/biological
  validity, and SG-6 by relative-ratio/absolute calibration;
- removed SG-6B as a hard prerequisite for annotation work;
- replaced circular viability reference wording with an independent
  object-level reference and a secondary well-level ATP reference;
- froze population-characterization versus measurement-grade SG-3 acceptance,
  independent-unit sample-size planning, and the SG-5
  method-freeze/implementation/simulation sequence;
- added a resolvable bibliography for the candidate methods.

- added `surface_weights_evidence_bearing` and bumped the feature schema to
  `2.4`;
- added explicit filled-void, cavity-topology and open-cavity provenance;
- changed viability fractions to use the predeclared classifiable denominator
  and report indeterminate fraction separately;
- changed the analytical gate to report qualification and unrestricted
  characterization as separate items.

The affected analytical record is rejected until the new source state is
committed and a fresh canonical V&V run is generated.

## Verification

- `pixi run check`: PASS (`No notebooks found`; notebook check not applicable).
- `pixi run regression`: PASS — Ruff passed; pytest reported 510 passed,
  8 skipped and 8 deselected; notebook validation was not applicable because
  the repository contains no notebooks. Existing dependency and numerical
  warnings remain warnings, not scientific evidence.
- `pixi run science-gate`: expected nonzero result. The analytical gate now
  reports SG-1a/SG-1b and SG-2a/SG-2b; the current record is rejected until the
  source-changing topology/schema migration is followed by a new canonical V&V
  run.

## Evidence conclusions

| Layer | Conclusion |
|---|---|
| Engineering behavior | focused topology/aggregation/provenance/gate tests pass; full regression pending |
| Analytical surface | `SUPPORTED WITH LIMITATIONS` in the declared domain; legacy gate remains unrestricted/FAIL until Option C record migration |
| Analytical volume | `NOT QUALIFIED` |
| Real-object morphology | `INSUFFICIENT EVIDENCE` |
| Absolute physical units | `INSUFFICIENT EVIDENCE` |
| Nuclear morphology | exploratory only |
| Biological viability | `NOT ASSESSED` |
| Statistical inference | `INSUFFICIENT EVIDENCE` |

Audit separation is Tier D: `LIMITED INDEPENDENCE`. The overall verdict remains
`NOT READY FOR THE SPECIFIED RESEARCH USE`.
