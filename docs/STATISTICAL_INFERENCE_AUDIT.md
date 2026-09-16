# Statistical inference audit — 2026-09-16

**Scope:** an audit of the inferential layer as implemented. It changes no
algorithm, no threshold, no statistical method and no gate status. SG-5 remains
`INSUFFICIENT EVIDENCE`. Nothing here qualifies an interval or selects a method.

This document exists because the required sequence for SG-5 is
**freeze method → implement → verify → simulate coverage**, and no primary method
is frozen yet (owner decisions D-9/D-10/D-11 in `docs/OWNER_DECISIONS.md`).
Qualifying the shipped branch and then replacing it would qualify a different
object than the one in use, so the audit records what is actually shipped while
the method choice remains open.

## 1. What is implemented (reconstructed from source)

| Component | Actual behaviour | Source |
|---|---|---|
| Unit of analysis | object is the measurement unit; biological replicate is the inference unit; group key is `condition::biological_replicate` (bare replicate labels repeat across conditions) | `statistics/inference.py:1-17, 137-141` |
| Primary model | `feature ~ condition` with a **biological-replicate random intercept**, REML, fit by `mixedlm` with `lbfgs` then `cg` | `statistics/inference.py:61-68` |
| Fallback | OLS with replicate-clustered SEs and `use_t=True` (statsmodels' cluster-robust **t(G−1)** reference) when the random effect is singular/fails | `statistics/inference.py:69-72` |
| Degenerate-RE guard | random-effect variance < `1e-6 × residual scale` treated as collapsed | `statistics/inference.py:33-37, 65-67` |
| Pairwise contrasts | manual **t(G−1)** reference applied to both branches; LMM branch otherwise has no small-sample correction | `statistics/inference.py:75-121` |
| Multiplicity | Benjamini–Hochberg **within one feature's own contrasts**, not across features | `statistics/inference.py:170, 193-198` |
| Omnibus | Wald test; for the LMM branch it is **asymptotic and uncorrected** and is documented as screening-only | `statistics/inference.py:56-59, 165-167, 199-202` |
| Transform | `volume_um3` log10 (heuristic); positive values enforced, no clipping | `statistics/inference.py:31, 151-152, 159-161` |
| Condition interval | percentile bootstrap of the mean of replicate-level median-of-unit-medians, 2000 resamples, pinned seed, **only at ≥ 3 replicates**, else `NaN` | `statistics/aggregation.py:15, 78-89`; `config.py:51-53` |
| Exploratory module | `statistics/exploration.py` — normality/variance checks, effect sizes, supervised classifiers | module level |

## 2. Findings that constrain any SG-5 claim

1. **No method is frozen, so the shipped path is not the qualification object.**
   The manual t(G−1) applied to the LMM contrasts is explicitly a *project
   heuristic*: sharing the clustered-OLS reference does not establish MixedLM
   small-sample coverage or type-I control, and it is not a
   Satterthwaite/Kenward–Roger correction
   (`statistics/inference.py:80-90`). Choosing the primary method is D-9.
2. **The percentile bootstrap is the most likely coverage failure.** A
   percentile bootstrap of a mean over 3 values resamples from three points and
   cannot extend beyond the observed range; the project's own floor is exactly
   `n = 3` (`aggregation.py:15`). It must be evaluated separately from the
   model's intervals. Interval construction is an *outcome* of the coverage
   study, not a predecided switch.
3. **Branch switching is a study variable.** A variance component collapsing
   near zero switches estimator mid-study; coverage and type-I error must be
   reported **per branch** (`LMM(lbfgs)`, `LMM(cg)`, `OLS(clustered SE)`), not
   averaged.
4. **CR2 is not available in this environment.** `statsmodels` provides no CR2
   sandwich estimator; a hand-written one is new code whose correctness must be
   verified against a published implementation on known cases and simulation
   before it can carry a coverage claim (D-11).
5. **`exploration.py` must not appear in reported inference.** Its classifier
   outputs and normality pre-checks are exploratory by name; a reported
   p-value or interval may not be sourced from it.

## 3. Required sequence before SG-5 can move

1. **Freeze** the primary method (Satterthwaite or between-within) and the
   fallback policy (D-9, D-11). Record it in `docs/ALGORITHM_DECISIONS.md` and
   `docs/PARAMETERS.md`.
2. **Implement** the frozen method; if CR2 is added, verify it against a
   published implementation (known cases and cross-package agreement).
3. **Simulate coverage** over the design grid the pipeline meets (ICC, replicate
   counts 2, 3, 4, 6, 10, 20, object-count imbalance, balanced/unbalanced),
   measuring type-I error under the null and 95% CI coverage, **per branch**,
   including the percentile-bootstrap interval evaluated separately.
4. **Predeclare the coverage tolerance** (D-10) and report no qualified interval
   in any cell that misses it. A cell that misses tolerance is a result
   (descriptive only), never a reason to switch methods.

## 4. Interim scope statement (unchanged from the master plan)

Until steps 1–4 exist, the inferential outputs are **descriptive summaries with
uncertainty indications of unverified coverage**. They may describe an
experiment's observed differences and help plan a confirmatory experiment. They
may not be the primary evidence for a biological conclusion, and a p-value from
this path must not be reported without the replicate count it rests on.
`exploration.py` outputs may not appear in a results claim at all.

## 5. What this audit does not do

- It does not run a coverage simulation (none is possible before the method is
  frozen without qualifying a branch scheduled for replacement).
- It does not choose Satterthwaite, between-within, CR2, BCa or any other
  method.
- It does not move SG-5: the gate stays `INSUFFICIENT EVIDENCE`.
