# SG-5 statistical qualification protocol (frozen 2026-09-26)

Written and hashed (`freeze.json`) before the confirmation simulation was
executed, and committed in the source commit that the record validates.
Sequence required by `SCIENTIFIC_VALIDATION_MASTER_PLAN.md` §10: method freeze →
implementation → engineering verification → simulation qualification.

## 1. Frozen method (D-9, D-11, D-16)

`organoid_analysis.statistics.small_sample`, `METHOD_ID =
reml_satterthwaite_cr2_fallback/1`:

| Component | Definition |
|---|---|
| Independent unit | biological replicate (grouping key = condition × replicate) |
| Model | `y = mu_condition + b_replicate + e_object`, REML variance components |
| Primary contrast | condition mean difference, t reference with Satterthwaite df (lmerTest definition) |
| Boundary fallback | `sigma_b^2 < 1e-6 sigma_e^2` (or one object per replicate): object-level OLS, CR2 covariance, Satterthwaite df (clubSandwich definition) |
| Omnibus | Wald F, denominator df `G − K`; screening only, not gating |
| Multiplicity | Benjamini–Hochberg within a feature (not part of the coverage claim) |
| Condition interval | `aggregation.condition_interval`: Student t on replicate summaries; log scale for positive size metrics (estimand: geometric mean) — replaces the percentile bootstrap |

## 2. Engineering verification (completed before this protocol)

`reference/`: 10 deterministic nested datasets (K = 2–4, balanced and
unbalanced, including three singular fits). Golden values from R 4.4.3,
lme4 2.0.6, lmerTest 3.2.1 and clubSandwich 0.5.10 (package from the anaconda
`r` channel, noarch build for R 4.3; pure R). Agreement: CR2 variance, df and
interval to 2.5e-13; REML variance components, SE and Satterthwaite df to
1.5e-7; boundary decision identical to `lme4::isSingular` in all 10. Pinned by
`tests/statistics/test_small_sample_reference.py` and re-checked by the record.

## 3. Development findings that fixed the envelope (seeds 101, 202, 303)

`docs/evidence/2026-09-26-statistical-method-development/`:

* 288-cell grid (K 2–3; G 2, 3, 4, 6, 10, 20; mean objects 5, 40, 150;
  balanced/unbalanced; ICC 0, 0.1, 0.3, 0.7), n_sim 3,000 per metric.
* At G ≥ 4 the contrast type-I error and coverage behave as nominal (z-score
  dispersion 1.0–1.1, mean type-I 0.050). At G = 2–3 the CR2 fallback is
  conservative in unbalanced designs at ICC ≈ 0 (type-I ≈ 0.03).
* Re-estimation at n_sim 20,000: unbalanced ICC 0 is 0.0425–0.044 at G = 4 and
  0.045 at G = 6; nominal at G = 10.
* **Qualified floor: 6 replicates per compared condition** (D-10 permits raising
  the floor, never lowering it; 4 would leave only 0.25 pp margin).
* The between-within omnibus is anti-conservative in unbalanced three-condition
  designs at G = 2 (0.11–0.13): it remains a screen.
* Robustness (non-gating): strongly skewed replicate effects make the method
  conservative (type-I 0.031–0.039), never anti-conservative, in all tested cells.
* The shipped percentile-bootstrap condition interval covered 75 % at G = 3 and
  85 % at G = 6 under exact normality; it is replaced.

## 4. Confirmation (executed once)

`coverage_confirmation.py`, seed **20260928**, **n_sim = 20,000** per metric per
cell (D-10 re-freeze: MCSE of a 5 % rate 0.154 pp, so the unchanged 1 pp
tolerance is 6.5 MCSE; family-wise chance of a spurious failure over 288 checks
< 0.1 % for an exact method).

* Contrast grid: K ∈ {2, 3}; G ∈ {6, 10, 20} (unbalanced: one condition G + 1);
  mean objects 5, 40, 150; balanced / unbalanced (log-normal objects per
  replicate, σ 0.6, minimum 2); ICC 0, 0.1, 0.3, 0.7 — 144 cells.
* Condition interval: G ∈ {6, 10, 20} × {normal, log-normal CV 0.3, 0.6} — 9 cells.

## 5. Acceptance rules

| Rule | Criterion |
|---|---|
| T1 | every contrast cell: type-I error at α = 0.05 in [0.04, 0.06] |
| T2 | every contrast cell: 95 % CI coverage in [0.94, 0.96] |
| T2c | every condition-interval cell: coverage in [0.94, 0.96] |
| V | all 144 + 9 cells present at n_sim 20,000, no fit failure, reference verification passes |

Any violation of T1/T2/T2c ⇒ FAIL. Violation of V ⇒ VOID (INSUFFICIENT EVIDENCE).
Omnibus type-I and per-branch rates are reported, not adjudicated.

## 6. What a PASS establishes

The frozen method delivers nominal type-I error and 95 % coverage for condition
contrasts in nested designs inside the envelope (2–3 conditions, 6–21
replicates per compared condition, 5–150 mean objects per replicate, ≥ 2
objects per replicate) when replicate effects are approximately normal on the
analysis scale. The runtime labels any contrast outside the envelope
`provisional - coverage unqualified`. It does **not** establish that a
particular experiment was designed correctly (replicate independence,
randomisation, pre-registration of the contrast), and it does not qualify the
omnibus screen, exploratory analyses, or effect sizes on non-normal scales.

## 7. Deviations recorded before the freeze

The condition-interval function was smoke-tested on two cells with seed 999
(n_sim 2,000) before this protocol was written; nothing was changed afterwards.
