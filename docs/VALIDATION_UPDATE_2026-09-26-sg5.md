# SG-5 statistical qualification — 2026-09-26 (second pass)

## Outcome

**SG-5 PASS** from canonical record
`evidence/2026-09-26-statistical-qualification-record/` (source `90c0fec`).
Gate after this pass: **4/9** (SG-1a, SG-2a, SG-4A, SG-5). SG-3A/3B, SG-4B and
SG-6A/6B still need laboratory data; see
[DATA_ACQUISITION_PACKAGE_SG3_SG4B_SG6.md](DATA_ACQUISITION_PACKAGE_SG3_SG4B_SG6.md).

What the PASS means, precisely: the frozen contrast method gives nominal
type-I error and 95 % coverage for condition contrasts in nested designs inside
the qualified envelope (2–3 conditions, **6–21 biological replicates per
compared condition**, 5–150 mean objects per replicate, ≥ 2 objects per
replicate) when replicate effects are approximately normal on the analysis
scale. It does not certify that a particular experiment is correctly designed
or pre-registered; every contrast outside the envelope is labelled
`provisional - coverage unqualified` at runtime.

## 1. Defects found in the shipped inference

| Finding | Evidence | Consequence |
|---|---|---|
| The primary method named by D-9/D-11 (Satterthwaite; CR2 fallback) had never been implemented | `statistics/inference.py` applied a heuristic t(G−1) to statsmodels MixedLM Wald statistics and CR0-type clustering in the fallback | No coverage claim could be made for the shipped path |
| Condition-level percentile bootstrap undercovers | Under exact normality: 75.2 % at 3 replicates, 84.7 % at 6, 92.8 % at 20 (4,000 simulations each) | The `ci95_*` columns in condition summaries were materially too narrow |
| Identical replicates produced tests | zero replicate-level variation was fitted and tested | now not estimable, feature skipped |

## 2. Frozen method (D-16)

`statistics/small_sample.py` (`reml_satterthwaite_cr2_fallback/1`):

* random-intercept REML fitted exactly on replicate sufficient statistics;
* each contrast: t reference with Satterthwaite df as in lmerTest
  (Kuznetsova et al. 2017, J. Stat. Softw. 82(13),
  [doi:10.18637/jss.v082.i13](https://doi.org/10.18637/jss.v082.i13));
* boundary fits (`sigma_b^2 < 1e-6 sigma_e^2`): object-level OLS with CR2 and
  Satterthwaite df as in clubSandwich (Pustejovsky & Tipton 2018, J. Bus. Econ.
  Stat. 36(4):672–683,
  [doi:10.1080/07350015.2016.1247004](https://doi.org/10.1080/07350015.2016.1247004));
* omnibus Wald F with `G − K` df — a screen, not gating;
* condition interval: replicate-level Student t; log scale for volume, surface
  area and equivalent diameter (interval for the geometric mean of replicate
  medians, exported with `ci95_centre_*` and `ci95_estimand_*`).

Both references were checked against the installed packages' own metadata
(lmerTest `citation()`; clubSandwich `DESCRIPTION` and its Wald-test vignette),
not through PubMed, which does not index these journals.

## 3. Engineering verification

R 4.4.3, lme4 2.0.6, lmerTest 3.2.1 and clubSandwich 0.5.10 on 10 deterministic
datasets (K = 2–4; balanced, unbalanced, singular). Maximum discrepancy: CR2
variance/df 2.8e-15; REML variance components and SE 5e-8; Satterthwaite df
1.5e-7; boundary decision identical to `lme4::isSingular` in 10/10.
clubSandwich is not on conda-forge and CRAN was unreachable, so the anaconda
`r` channel noarch build (for R 4.3, pure R) was loaded into R 4.4.3; this is
recorded as a provenance limitation.

## 4. Development and confirmation

| Stage | Design | Result |
|---|---|---|
| Development (seed 101) | 288 cells, G 2–20, n_sim 3,000 | nominal from G = 4; conservative at G = 2–3 in unbalanced ICC ≈ 0 designs (type-I ≈ 0.03) |
| Targeted (seed 202) | 9 cells, n_sim 20,000 | unbalanced ICC 0: 0.0425 at G = 4, 0.045 at G = 6, nominal at G = 10 → floor raised to 6 (D-10) |
| Robustness (seed 303, non-gating) | strongly skewed replicate effects | conservative (type-I 0.031–0.039), never anti-conservative |
| **Confirmation (seed 20260928, frozen in `8d59602`…`90c0fec`)** | 144 contrast cells + 9 interval cells, n_sim 20,000 | **type-I 0.0454–0.0530; coverage 0.9468–0.9555; interval coverage 0.9486–0.9514; 0 violations** |

Omnibus screen type-I in the confirmation grid: 0.046–0.079 (anti-conservative
in unbalanced three-condition designs, as D-9 anticipated); it carries no claim.

## 5. What changes for users

* `pairwise_contrasts.csv` gains `branch`, `interval_qualification`,
  `qualification_reason` and `method_id`; `df_denom` is now the Satterthwaite
  df (non-integer).
* `stats_results.json` gains the omnibus F, its df and a note.
* `condition_summary.csv`: `ci95_*` is now a t interval (geometric mean for size
  metrics) with `ci95_centre_*`, `ci95_estimand_*`, `ci95_qualification_*` and
  `ci95_method`. Earlier bootstrap intervals are not comparable.
* To obtain qualified intervals an experiment needs **≥ 6 independent biological
  replicates per compared condition**. With 3–5 the numbers are computed but
  labelled provisional.

## 6. Process note

The shared recording script is hash-bound by every qualification record, so
extending it for SG-5 invalidated the SG-1a/SG-2a record until it was
re-executed (`2026-09-26-analytical-qualification-record-2`, same frozen grids,
unchanged estimator). This is the evidence system working as designed.
