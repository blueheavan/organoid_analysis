# Scientific validation master plan — frozen 2026-09-13

This document governs scientific qualification. Every output must trace to an
estimand, applicability domain, reference standard, acceptance rule and
validation record. Making a gate green is not the objective.

No study result is created here. Existing `NOT ASSESSED`, `INSUFFICIENT
EVIDENCE` and `NOT QUALIFIED` states remain until prospective evidence exists.

## 1. Frozen owner decisions

### A. Gate architecture

Adopt Option C from `GATE_SEMANTICS_ANALYSIS.md`:

- qualification uses `PASS`, `FAIL` or `INSUFFICIENT EVIDENCE` in a declared
  intended-use domain;
- characterization reports distributions, strata, medians and worst cases
  without a PASS/FAIL verdict;
- surface and volume become SG-1a/SG-1b and SG-2a/SG-2b;
- characterization evidence cannot satisfy qualification.

The gate now reports qualification and characterization separately. Qualification
remains fail-closed until a dedicated domain-restricted canonical record exists.

### B. Organoid volume estimand

Primary: outer-envelope volume, `V_env`. Secondary: segmented-material volume,
`V_seg`. Every `V_env` report also carries
`f_void = (V_env - V_seg) / V_env`.

Fragmented or merged labels are segmentation/QC states, not automatically one
biological organoid. Open cavities require an explicit topology flag because
topological hole filling can change the measurement basis discontinuously.

### C. Meaning of the `<1%` volume criterion

`|relative analytical volume error| < 1%` is retained as a conservative
engineering qualification target for the voxel-count estimator. It is not a
biologically required threshold. A real-segmentation requirement must be
justified and frozen separately before confirmation.

## 2. Intended uses and present claim boundary

| Intended use | Outputs | Present claim |
|---|---|---|
| IU-1 organoid 3D morphology | volume, area, equivalent diameter, axes; derived sphericity, surface-to-volume ratio, elongation and axis ratios | analytical surface: `SUPPORTED WITH LIMITATIONS`; real-object accuracy and absolute units: `INSUFFICIENT EVIDENCE`; analytical volume qualification: `NOT QUALIFIED` |
| IU-2 nuclear morphology | volume, equivalent diameter, axes, elongation and axis ratios; conditional area/sphericity/surface-to-volume ratio | `EXPLORATORY`; not a core-gate blocker |
| IU-3 Calcein AM / PI | four marker-signal states and aggregate fractions | marker-state computation exists; biological viability meaning is `NOT ASSESSED` |

Nuclear area-derived values outside the numerical surface domain remain
exploratory descriptors and carry `OUTSIDE QUALIFIED NUMERICAL DOMAIN`. They
are not validated nuclear surface measurements.

## 3. Domain and provenance rules

- `rho_in = r_in / max(spacing)`, not the voxel diagonal.
- Machine-checkable numerical domain and non-machine-checkable assumptions are
  separate. Smoothness is an intended-use assumption.
- Reports separate declared domain, confirmation-covered range and
  development-supported extrapolation.
- Analytical geometry error is one component of total production measurement
  error; it is not a lower bound on production error.
- Within sampled phantom families, resolution was a major observed
  determinant; anisotropy and shape class also had measurable effects.
- `surface_weights_origin` identifies packaged/cache/solved weights, while
  `surface_weights_evidence_bearing` states whether the exact vector contributed
  to frozen confirmation evidence. The schema migration requires a replacement
  canonical V&V run.

## 4. Gate definitions

| Gate | Scope |
|---|---|
| SG-1a | Surface analytical qualification in the declared domain |
| SG-1b | Unrestricted surface characterization; no verdict |
| SG-2a | Volume analytical qualification in a volume-specific domain |
| SG-2b | Unrestricted volume characterization; no verdict |
| SG-3A | Brightfield organoid segmentation validation |
| SG-3B | Membrane-fluorescence organoid segmentation validation |
| SG-4A | Analytical correctness of four-state fluorescence classification |
| SG-4B | Biological validity of the viability interpretation |
| SG-5 | Statistical inference qualification |
| SG-6A | Relative XY and Z:XY spacing-ratio verification |
| SG-6B | Absolute physical calibration and channel registration |

Brightfield and membrane fluorescence do not share a PASS without a
predeclared transport study.

## 5. SG-1 surface work

Do not reopen Crofton-estimator optimization. The bounded analytical claim is
`SUPPORTED WITH LIMITATIONS` for sufficiently resolved smooth closed objects
in the declared domain and evidence-bearing weight configuration. Extend
confirmation around `rho_in = 10`, to higher `rho_in`, and near anisotropy 1
and 4. SG-3 determines applicability to real organoids.

## 6. SG-2 volume V&V

Volume is `NOT QUALIFIED` pending a dedicated protocol. Development covers
sphere, ellipsoid, capsule, torus, cylinder, box, hollow sphere, shell,
irregular smooth, branched, touching-pair and border-truncated objects;
resolution below/boundary/inside/high; dense subvoxel phase; orientation; and
anisotropy 1:1:1 through 4:1:1 plus real acquisition ratios.

Development may alter the candidate domain. Confirmation uses a different
frozen seed, hashes the generated set before execution, is not inspected before
evaluation, and runs once. If `<1%` is retained, every in-domain confirmation
case must meet it. A median cannot qualify a per-case claim.

## 7. SG-6 acquisition calibration

SG-6A is the minimum prerequisite for SG-3 stratification: verify X:Y and Z:XY
ratios and uncertainty so `rho_in` and anisotropy are not misclassified.
Absolute scale is not required for production/reference relative discrepancy
on the same voxel grid because the common volume factor cancels.

SG-6B is required for absolute µm, µm² and µm³ claims. Measure X, Y and Z
separately using traceable lateral/axial standards, and channel offsets using
multispectral beads. Report per-axis bias, uncertainty, between-session and
field-position variation, depth effects and registration.

For volume, `delta V / V ≈ delta x/x + delta y/y + delta z/z`. Surface has no
universal coefficient; propagation depends on shape, orientation and error
direction. Stability, absolute bias and measurement uncertainty have separate
budgets. The previous `0.25%` is not reused for each component.

## 8. SG-3 real-organoid segmentation

Use at least two independent blinded annotators; three are preferred.
Annotators do not see production results or each other's labels. STAPLE-like
fusion may define consensus, but each annotator-versus-consensus and inter-rater
Dice, IoU and boundary agreement remain visible.

Level 1 characterizes detection sensitivity/precision, F1, IoU/Dice
distributions, boundary distance and Hausdorff distance. Medians, quantiles and
intervals are descriptive.

Level 2 supports measurement claims. Per matched organoid, report signed and
absolute errors, median, upper tail, worst case and uncertainty for `V_env`,
`V_seg`, `f_void` and applicable shape measures. Acceptance separately covers
systematic bias, tail/precision and catastrophic failure rate. Values come from
a pilot, biological requirement and precision calculation; they are not copied
from the analytical estimator criterion.

Sample size targets bias, upper quantile, failure rate and inter-rater
variability. The biological sample/experiment is independent; objects describe
within-sample distributions and do not inflate biological replication.

**Recorded (D-5/D-6, 2026-09-16).** δ, the acceptance criteria C1–C6 and the
precision-derived N are predeclared in
`evidence/2026-09-11-measurement-vv/SEGMENTATION_VALIDATION_PROTOCOL.md` §5:
δ = 0.20 relative for volume and area and 0.10 for diameter, sphericity and
principal-axis length; per-case upper-quantile bias criteria; a detection
criterion on the lower confidence bound of recall and precision; and ≥ 47
biological samples per stratum, ≥ 20 per compared condition. **The inherited
`1% / 5%` pair is not SG-3's criterion** — it governs SG-1a/SG-2a, which qualify
the estimator on analytical cases. SG-3's load-bearing criterion is
**differential**: the 95% CI of the difference in median relative bias between
any two compared conditions must exclude ±2.5%, because a bias common to all
objects cancels in a ratio while a bias that differs between arms does not.

## 9. SG-4 viability estimand and reference

Object outputs remain `viable_like`, `mixed_signal`, `compromised_like` and
`indeterminate`; they do not become live/dead identities until SG-4B passes.

Target aggregate estimands are

`f_viable_like = N_viable_like / N_classifiable`

and

`f_indeterminate = N_indeterminate / N_morphology_QC_eligible`,

where `N_classifiable` is viable-like + mixed-signal + compromised-like. An
excessive predeclared indeterminate fraction blocks a biological claim. The
implementation matches this estimand as of the 2026-09-13 schema migration: the
three classifiable-state fractions divide by `N_classifiable`,
`fraction_indeterminate` divides by `N_morphology_QC_eligible`, and that
denominator is exported as `fraction_indeterminate_denominator`
(`statistics/aggregation.py:19, 28-31`; contract test
`tests/phenotyping/test_viability_summary.py:125`). The biological claim still
requires the independent object-level reference below (owner decision D-8).

The primary object reference is an independent blinded manual, semi-manual, or
separately validated nuclear dead-stain measurement path. It cannot depend on
the production segmentation branch. A 3D ATP assay is a secondary well-level
reference and cannot establish object-level accuracy. A graded insult series
tests monotonicity and dynamic range, not identity by itself.

Primary metrics are the confusion matrix, quadratically weighted Cohen's
kappa, macro-F1, and class-specific sensitivity and precision with intervals.
Spearman correlation is trend characterization only. Stratify by organoid size,
imaging depth, radial location, batch, staining duration and acquisition.

**Recorded (D-7, 2026-09-16).** Per-class sensitivity and precision must reach a
**lower 95% CI bound of 0.80**, and macro-F1 a lower 95% bound of 0.80;
**quadratic-weighted kappa and Spearman are reported but do not gate**, because
kappa behaves incoherently under class imbalance and a worse classifier can score
higher. A class with prevalence below 5% is reported `not evaluable`, not FAIL.
The indeterminate ceiling is **`f_indeterminate ≤ 0.20` per condition together
with a between-condition difference no larger than 5 percentage points on its
upper 95% bound**; above either, V2 is not reportable for that condition. The
differential half is load-bearing: `f_viable_like` divides by `N_classifiable`,
so an indeterminate rate that differs between arms biases the contrast even when
both arms are under the ceiling.

**Recorded (D-8, 2026-09-16).** The primary object-level reference is **blinded
human classification of the same objects**, on a predeclared stratified
subsample, by two annotators, with the reference's own disagreement reported.
The ATP assay remains the secondary well-level reference and R-13 stays out of
scope. The annotation set is the same asset SG-3 requires, so the two studies
share one campaign, and the claim V1 licenses is bounded by the reference's own
accuracy.

## 10. SG-5 statistical inference

Sequence: method freeze → implementation → engineering verification →
simulation qualification. Do not qualify the shipped branch and then replace
it.

**Recorded (D-9, 2026-09-16).** Satterthwaite degrees of freedom for the primary
fixed-effect contrasts. **Between-within is used only for the omnibus test**,
which is reported with a between-within df check and an explicit
anti-conservative-at-`G < 10` statement. Kenward–Roger remains a sensitivity
analysis where a verified implementation exists. The primary claim is a
per-contrast claim, and Satterthwaite is a per-contrast approximation whereas
between-within is an omnibus df defined for a balanced nested design.

**Recorded (D-11, 2026-09-16).** The CR2 fallback is **implemented locally and
verified against a published implementation**: golden values from `clubSandwich`
with its version recorded, analytic cases where CR2 must equal CR0, the published
cases and degrees-of-freedom values of Pustejovsky & Tipton (2016), and then the
D-10 simulation run on the frozen implementation only. "Accept a reliable
external implementation" is not available in the pinned environment —
`statsmodels` 0.14.6 provides CR0-style `cov_type="cluster"` and no bias-reduced
linearization — so the choice is between unverified and verified local code, and
only the latter can carry a coverage claim.

Evaluate LMM and fallback branches separately over ICC, replicate count,
object-count imbalance and balanced/unbalanced designs:

- T1: type-I error under no treatment effect;
- T2: 95% CI coverage at replicate counts 2, 3, 4, 6, 10, 20 and added
  design-specific values;
- T3: pooled-object, replicate-weighted and mean-of-replicate-medians
  estimands.

Choose CI methods only from demonstrated coverage. Candidates are percentile,
BCa, bootstrap-t and replicate-level t intervals. Cells missing the frozen
coverage tolerance are descriptive only.

**Recorded (D-10, 2026-09-16).** The frozen tolerance is **1 percentage point
from nominal, against a simulation MCSE of at most 0.33 pp**, which fixes
`n_sim = 5 000` per cell: `MCSE = sqrt(0.95 × 0.05 / 5 000) = 0.31` pp, so the
tolerance is more than three Monte-Carlo standard errors wide and is a statement
about the method rather than about the simulation. Type-I error at nominal 0.05
must lie in [0.04, 0.06]. An interval may be reported only for cells that pass
T2 coverage; the shipped floor of 3 replicates (`statistics/aggregation.py:15`)
is retained as the computation floor but is **provisional**, so any interval
reported at 2–3 replicates is labelled `provisional — coverage unqualified` and
carries no claim. The study may raise the reportable floor; it may not lower it.

## 11. Dependencies

```text
Volume V&V -----------------------------------+
                                               |
SG-6A spacing-ratio check --> SG-3 ------------+--> real-object morphology
                                               |     (relative comparison)
SG-6B absolute calibration --------------------+--> absolute physical units
                                               |
                                               +--> SG-4 biological viability

Statistics method freeze --> implementation --> simulation qualification
                                                --> study-level inference
```

SG-6A is the only calibration dependency of SG-3: the SG-3 downstream-bias
metric is a same-grid ratio, so absolute scale cancels exactly for volume and to
0.007 pp for area even at a 20% axial error (E-6 in
`evidence/2026-09-13-volume-audit/ERRATA.md`). SG-6B governs absolute-unit
claims and runs in parallel; it is not on the SG-3 critical path.

SG-3 annotation preparation starts without waiting for SG-6B. Volume V&V,
annotation design/training/sample-size planning, calibration work and
statistical-method work proceed in parallel. SG-4B follows a stable SG-3 object
definition.

## 12. Release rule

Every claim answers: estimand, input, reference, criterion, criterion rationale,
domain, confirmation independence, one-time prospective execution, current
input's domain membership and uncovered errors. Missing answers retain
`NOT ASSESSED`, `INSUFFICIENT EVIDENCE` or `NOT QUALIFIED`.

## 13. References

These sources justify candidate methods, not project-specific thresholds.

Traceability: the entries below are this repository's reference set, each with a
DOI or stable URL. They are cited as *candidate-method anchors*; none supplies a
threshold (see items 1–8). Full-text-reviewed versus abstract-only status is
**not recorded** for these entries and must not be assumed: where a future study
relies on a specific published figure, that study's design record must state the
access depth next to the citation. The externally supplied solutions report
reviewed in `docs/reviews/2026-09-13-solution-report-review.md` carried more
than forty unresolved bracketed citations (review defect B1); it is **not** filed
as an evidence document in this repository, and no claim here depends on it.

1. Warfield SK, Zou KH, Wells WM. STAPLE. *IEEE TMI*. 2004;23:903-921.
   https://doi.org/10.1109/TMI.2004.828354
2. Cohen J. Weighted kappa. *Psychological Bulletin*. 1968;70:213-220.
   https://doi.org/10.1037/h0026256
3. Reinke A, et al. Metric-related pitfalls in image analysis validation.
   *Nature Methods*. 2024. https://doi.org/10.1038/s41592-023-02150-0
4. Kuznetsova A, Brockhoff PB, Christensen RHB. lmerTest. *JSS*. 2017;82(13).
   https://doi.org/10.18637/jss.v082.i13
5. Bell RM, McCaffrey DF. Bias reduction in standard errors for linear
   regression with multi-stage samples. *Survey Methodology*. 2002;28:169-181.
   https://www150.statcan.gc.ca/n1/pub/12-001-x/2002002/article/9058-eng.pdf
6. MacKinnon JG, Nielsen MO, Webb MD. Cluster-robust inference. *Journal of
   Econometrics*. 2023;232:272-299.
   https://doi.org/10.1016/j.jeconom.2022.04.001
7. Schenker N. Qualms about BCa bootstrap confidence intervals. *Statistics &
   Probability Letters*. 1992;13:381-385.
   https://doi.org/10.1016/0167-7152(92)90288-G
8. Promega. CellTiter-Glo 3D Cell Viability Assay Technical Manual, TM412,
   revised 2023.
   https://www.promega.com/resources/protocols/technical-manuals/101/celltiter-glo-3d-cell-viability-assay-protocol/
