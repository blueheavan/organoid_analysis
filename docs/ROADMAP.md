# Validation roadmap — 2026-09-13

This roadmap implements `SCIENTIFIC_VALIDATION_MASTER_PLAN.md`. It separates
qualification from characterization and keeps missing studies visible. No item
acquires a PASS from this document.

Overall release verdict: **NOT READY FOR THE SPECIFIED RESEARCH USE**.

## 1. Present claim boundary

| Claim | Status |
|---|---|
| Surface analytical measurement in the declared smooth/closed, resolved domain and evidence-bearing weight configuration | `SUPPORTED WITH LIMITATIONS` |
| Surface outside that domain | characterization only |
| Volume analytical qualification | `NOT QUALIFIED` |
| Real-organoid morphology accuracy | `INSUFFICIENT EVIDENCE` |
| Absolute physical dimensions | `INSUFFICIENT EVIDENCE` |
| Nuclear morphology | `EXPLORATORY` |
| Four-state marker computation | software behavior exists; biological interpretation `NOT ASSESSED` |
| Statistical inferential validity | `INSUFFICIENT EVIDENCE` |
| Engineering and provenance contracts covered by current tests | supported only as software contracts |

## 2. P0 — frozen without new experiments

| ID | Action | State after this revision |
|---|---|---|
| P0-1 | Correct `rho_in`, production-error, shape/resolution and per-axis calibration wording | recorded in `evidence/2026-09-13-volume-audit/ERRATA.md` and the master plan |
| P0-2 | Freeze IU-1/IU-2/IU-3 | complete in `INTENDED_USE_AND_ESTIMANDS.md` and the master plan |
| P0-3 | Freeze `V_env` primary, `V_seg` secondary and `f_void` companion | complete |
| P0-4 | Select Option C gate semantics | owner decision recorded; code migration awaits a dedicated record contract |
| P0-5 | Classify volume `<1%` | conservative engineering target, not biological necessity |
| P0-6 | Require independent viability reference | complete in the estimand and validation plan |
| P0-7 | Remove SG-6B as an SG-3 hard prerequisite | complete; SG-6A ratio verification remains required for stratification |
| P0-8 | Freeze SG-3/SG-4/SG-5 acceptance frameworks | complete in `VALIDATION_PLAN.md` and the master plan |
| P0-9 | Resolve external references | bibliography included in the master plan; project thresholds remain study decisions |

## 3. Parallel work tracks

### Track A — volume analytical V&V

1. Add topology diagnostics for closed voids and suspected open cavities.
2. Freeze the open-cavity policy and a volume-specific domain.
3. Generate development cases across shape, resolution, phase, orientation and
   anisotropy.
4. Freeze a separate seed and artifact hashes for confirmation.
5. Execute confirmation once; every in-domain case must meet the frozen
   engineering criterion if `<1%` is retained.
6. Report unrestricted behavior under SG-2b without a verdict.

### Track B — SG-3 segmentation

1. Maintain separate protocols for brightfield (SG-3A) and membrane
   fluorescence (SG-3B).
2. Train at least two blinded annotators; prefer three.
3. Freeze consensus, inter-rater, detection, boundary and measurement-bias
   metrics before data review.
4. Calculate sample size from precision for bias, upper-tail error, failure
   rate and inter-rater variability at the biological-sample level.
5. Record model checkpoint, specimen type, instrument and acquisition settings.

Annotation preparation starts immediately. Final domain stratification requires
SG-6A, not SG-6B.

### Track C — acquisition calibration

- SG-6A verifies X:Y and Z:XY ratios and uncertainty for `rho_in`/anisotropy.
- SG-6B measures per-axis absolute bias, uncertainty, session/field/depth
  variation and channel registration with traceable standards.
- Stability, accuracy and uncertainty budgets remain separate.

### Track D — statistical inference

1. Select and freeze Satterthwaite or between-within as the primary LMM method.
2. Select an external CR2 implementation or specify an independently checked
   local implementation.
3. Implement and verify known cases and cross-package agreement.
4. Simulate type-I error, CI coverage and estimand behavior separately for LMM
   and fallback branches.
5. Permit inferential intervals only in cells meeting the frozen coverage
   tolerance.

## 4. Gate implementation work

| Gate | Required record before PASS | Current state |
|---|---|---|
| SG-1a | canonical domain-restricted confirmation record including `surface_weights_evidence_bearing` | record migration pending; existing scientific claim remains bounded |
| SG-1b | unrestricted analytical grid | characterized; no verdict |
| SG-2a | volume-specific domain and untouched confirmation record | `INSUFFICIENT EVIDENCE`; analytical claim `NOT QUALIFIED` |
| SG-2b | unrestricted volume grid | characterized; no verdict |
| SG-3A | qualified brightfield annotation record | `INSUFFICIENT EVIDENCE` |
| SG-3B | qualified membrane-fluorescence annotation record | `INSUFFICIENT EVIDENCE` |
| SG-4A | versioned classification/denominator contract and adversarial tests | partial software evidence only |
| SG-4B | independent object-level reference plus orthogonal well-level evidence | `NOT ASSESSED` |
| SG-5 | frozen method plus branch-specific simulations | `INSUFFICIENT EVIDENCE` |
| SG-6A | ratio-standard study | `NOT ASSESSED` |
| SG-6B | traceable absolute calibration/registration study | `NOT ASSESSED` |

Characterization items never contribute to the qualification PASS count.
Qualification fails closed on missing domains/flags, empty subsets,
noncanonical records or changed dependencies.

## 5. Dependencies

```text
Volume V&V -----------------------------------+
                                               |
SG-6A spacing-ratio check --> SG-3 ------------+--> real-object morphology
                                |              |
SG-6B absolute calibration --------------------+
                                |
                                +--> SG-4B biological viability

Statistics method freeze --> implementation --> simulation qualification
                                                --> study-level inference
```

Volume V&V and SG-6B are independent. SG-3 preparation and SG-6B run in
parallel. SG-4B follows a stable SG-3 object definition. Statistical
qualification cannot compensate for biased measurement or invalid objects.

## 6. Remaining owner/study decisions

| Decision | Blocks |
|---|---|
| open-cavity behavior: exclude, flag-only, or a newly defined closing estimand | volume domain and confirmation |
| volume-specific applicability domain | SG-2a |
| SG-3 bias, upper-tail and catastrophic-failure criteria | SG-3A/SG-3B |
| indeterminate-rate ceiling and SG-4 class-performance criteria | SG-4B |
| Satterthwaite versus between-within primary LMM method | SG-5 implementation |
| simulation coverage tolerance | SG-5 qualification |
| SG-6A/SG-6B uncertainty and stability acceptance values by instrument stratum | physical-unit claims |
| calibration requalification interval | maintained SG-6 status |

## 7. Evidence sequence for release

`acquisition → segmentation → geometry → biological readout → inference`

Each claim identifies its estimand, input, reference, acceptance rule and
rationale, domain, independent confirmation data, prospective execution,
current-input domain membership and uncovered errors. Missing answers remain
`NOT ASSESSED`, `INSUFFICIENT EVIDENCE` or `NOT QUALIFIED`.
