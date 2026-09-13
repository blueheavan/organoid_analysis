# Study proposals: segmentation accuracy (SG-3) and acquisition calibration (SG-6)

**Status: PROPOSALS ONLY. Neither study has been executed, no data have been
collected, and nothing here asserts or anticipates a PASS.** SG-3 remains
**FAIL / INSUFFICIENT EVIDENCE** and SG-6 remains **FAIL / INSUFFICIENT
EVIDENCE** until a predeclared study is executed and recorded. Acceptance
numbers below marked *illustrative* are placeholders the project owner must
replace with predeclared values before any data are seen; they are not criteria.

**Explicit non-claim.** The analytical geometry evidence in
`docs/evidence/2026-09-13-analytical-geometry-record` and
`docs/evidence/2026-09-12-surface-method-development` concerns exact
rasterizations of analytical solids. It is **not** evidence about segmented
objects and is not offered here as any part of the case for segmentation
validity. A qualified estimator applied to a wrong mask gives a wrong answer
with no warning, and the two error sources have not been characterized jointly.

Written 2026-09-13.

## 1. SG-3 — segmentation accuracy: addendum to the existing protocol

A predeclared protocol already exists at
`docs/evidence/2026-09-11-measurement-vv/SEGMENTATION_VALIDATION_PROTOCOL.md`
and is not restated here. It already specifies the independent reference
standard and its independence conditions, blinding, the training-overlap check,
inter-rater variability, strata, the biological sample as the independent unit,
object matching, and the metric hierarchy from detection through mask overlap to
downstream measurement bias — including the warning that overlap agreement alone
never supports a claim that morphology measurements are validated. This addendum
records only what changed since it was written.

**1.1 Per-axis spacing provenance is now available and should be a recorded
covariate.** As of this round each loaded sample records `spacing_source_by_axis`
and, where an external manifest completed axes the file omitted, the mixed
provenance is preserved rather than collapsed. The annotation study should
record, per biological sample, which axes' spacing came from file metadata and
which from the manifest, because a spacing error propagates linearly into volume
and roughly quadratically into area, and would otherwise be indistinguishable
from segmentation bias in the downstream-bias analysis.

**1.2 The domain flag must be reported, not used as an exclusion.** Every
measured object already carries `surface_in_qualified_domain` and
`surface_domain_flags`. The study should report the downstream-bias metrics
**stratified** by that flag, and must not silently drop out-of-domain objects:
the fraction of real objects that fall outside the checkable domain is itself a
primary result, because it bounds how much of the intended use the qualified
estimator can serve. If that fraction is large, a PASS on the in-domain stratum
would be a narrow claim and must be reported as such.

**1.3 The smoothness clause cannot be checked on real objects either.** Per
`docs/GATE_SEMANTICS_ANALYSIS.md` §2, the non-machine-checkable part of the
domain is load-bearing: creased phantoms inside the checkable domain reach 9.10%
area and 3.20% volume error. Real organoid envelopes are not analytically
smooth. The annotation study therefore cannot verify domain membership from the
masks; the most it can do is measure the *total* discrepancy against the
reference standard, which folds segmentation error and out-of-scope estimator
error together. Separating them requires the analytical extension in
`VOLUME_QUALIFICATION_PROTOCOL.md` §3 (creased and hollow classes) and should be
stated as a limitation of the SG-3 result rather than resolved inside it.

**1.4 The reference standard must be measured on the same basis as the
production output.** Added in round 4, following the estimand decision in
`docs/INTENDED_USE_AND_ESTIMANDS.md` §3. The production organoid volume is
`V_env` — the envelope after topological hole filling — while the annotator
produces a mask. If the reference volume is computed as a raw voxel count and
the production volume as an envelope count, the comparison charges the lumen
handling to segmentation error. The protocol must therefore state that reference
masks pass through the *same* `measurement_policy("organoid", "filled_envelope")`
path, and that the downstream-bias analysis reports `V_env`, `V_seg` and
`f_void` separately. The same requirement applies to nuclear and cell-level
reference masks, which are measured on `raw_label` — the basis must match per
level, not globally.

**1.5 Topology diagnostics are a required covariate, not an optional one.** The
open-cavity case (`INTENDED_USE_AND_ESTIMANDS.md` §3.3) makes `V_env` a
discontinuous function of mask topology: an annotator and a segmenter that
differ only in whether a concavity is closed will differ in `V_env` by the whole
lumen volume while agreeing almost everywhere voxel-wise. Overlap metrics will
not show this, and it is not a small effect. The study must record, per object
and per mask source, the filled-void count and whether an adjacent background
component reaches the image border (roadmap R-3), and must report the
downstream-bias metrics stratified by topology agreement. Objects where the two
masks disagree on topology are a separate reporting category; pooling them with
the rest understates both the agreement and the disagreement.

**1.6 Scope of a possible SG-3 PASS.** Whatever the study returns, it can
support a claim only about the segmentation configuration, specimen type,
instrument and acquisition settings it sampled. It is not transferable to
another model checkpoint or another microscope, and the record must name all
four.

## 2. SG-6 — acquisition calibration and channel registration: proposed study

No protocol exists for this item. The gate's current rationale is that metadata
is read and traced but not independently verified — that is, the pipeline trusts
the voxel size the file declares. This proposal specifies what independent
verification would consist of. It is a wet/instrument study, not a software
study.

### 2.1 Question

For each instrument and objective in the intended use: is the voxel size
recorded in the image metadata equal, within a predeclared tolerance, to the
voxel size measured against a traceable physical standard; and is the spatial
offset between channels smaller than a predeclared fraction of the objects being
measured?

This matters quantitatively and is why it cannot be deferred behind the
estimator work: a scale error is applied to *every* object, is perfectly
correlated across objects, and is invisible to every check in this repository.
The propagation, derived rather than asserted — an earlier form of this passage
stated it incorrectly ([ERRATA](ERRATA.md), E-5):

Write the reported spacing as `s_i = (1 + e_i) · s_i^true` for `i ∈ {z, y, x}`.
Voxel-count volume scales with the product of the three spacings, so

    e_V = (1 + e_z)(1 + e_y)(1 + e_x) − 1 ≈ e_z + e_y + e_x   (first order)

Surface area has no single factor: a surface element spanned by axes `i` and `j`
scales by `(1 + e_i)(1 + e_j)`, so the area error depends on the object's
orientation distribution and lies between the smallest and largest such pairwise
product. Concretely, for a **lateral-only** error `e_y = e_x = e`, `e_z = 0`:

| error | first order | exact at e = 3% |
|---|---|---|
| length, lateral | `e` | +3.00% |
| volume | `2e` | +6.09% |
| area, bound | `e` to `2e` | +3.00% to +6.09% |
| area, sphere (exact oblate-spheroid formula) | `(4/3)e` | +4.02% |

For an **isotropic** error in all three axes, volume goes as `3e` (+9.27% at
e = 3%) and area as `2e` (+6.09%). Either way the systematic term at a
3% miscalibration is several times the qualified estimator's in-domain area
budget (worst case 0.95% on the confirmation set), which is why calibration must
be evidenced before an estimator's accuracy claim means anything in physical
units. The sphere figure above is a numerical evaluation for one shape and is
illustrative; the study propagates the measured deviations per object using the
exact pairwise form, not a single factor.

#### 2.1.1 Shape dependence, and the two quantities that are invariant

Sensitivities as logarithmic derivatives at the true scale
(`d ln Q / d ln k`, where `k = 1 + e` is the lateral factor and the polar axis
is exact). Spheroids with equatorial semi-axis `a` and polar `c`, aligned with
z; computed from the exact spheroid area formula:

| body | area | volume | sphericity |
|---|---|---|---|
| oblate, `c/a = 0.5` | 1.599 | 2 | −0.266 |
| sphere, `c/a = 1` | 1.333 (= 4/3, exact) | 2 | 0 (exact) |
| prolate, `c/a = 2` | 1.138 | 2 | +0.195 |
| prolate, `c/a = 4` | 1.046 | 2 | +0.287 |

Two consequences that a naive treatment gets wrong:

1. **Sphericity is exactly invariant to an isotropic scale error** — it is
   dimensionless, and `(6V)^{2/3}/A` scales as `k^2/k^2`. Under a *lateral-only*
   error the cancellation is exact only for a sphere (`4/3 = (2/3)·2`) and
   fails in both directions for elongated objects, up to ≈0.3 in log-derivative
   at the domain's anisotropy limits. Sphericity must therefore not be described
   as calibration-free; it is free of *isotropic* scale error only.
2. **`rho_in` is also invariant to an isotropic scale error.** It is
   `r_in / max(spacing)` with `r_in` an exact Euclidean distance transform in
   physical units (`surface_crofton.py:418`, `sampling=spacing`): scaling all
   three spacings by `k` scales both numerator and denominator by `k`. So
   **absolute** scale error cannot move an object across the domain boundary —
   only an error in the z:xy *ratio* can. This splits SG-6 into two
   requirements with different consequences: the axis-ratio calibration governs
   domain membership and the anisotropy flag, the absolute scale governs the
   physical units of every reported quantity, and a study that establishes one
   does not establish the other.

### 2.2 Reference standards (must be independent of the instrument's own
metadata)

- **Lateral (XY) scale**: a stage micrometer or a calibrated grid slide with a
  certificate of traceability; alternatively a certified bead-spacing standard.
  The instrument's nominal magnification is *not* a reference standard.
- **Axial (Z) scale**: a certified axial standard — e.g. a calibrated
  multi-layer slide or beads embedded at a certified separation. A nominal piezo
  step is not a reference standard; the commanded step and the delivered step
  are the quantities in question.
- **Refractive-index mismatch**: because the axial scale of a confocal stack
  depends on the mounting medium, the standard must be imaged in the same
  medium and at the same depth range as the specimens, or a measured correction
  factor must be recorded with its own uncertainty. A calibration performed in
  a different medium does not transfer.
- **Channel registration**: sub-resolution multi-spectral beads imaged in the
  same channel configuration as the assay.

### 2.3 Design

1. **Declare intended use**: the instrument/objective/medium/step-size
   combinations to be qualified. Each combination is a separate stratum; a
   result for one objective says nothing about another.
2. **Sampling**: for each combination, image the standard on at least three
   separate days (illustrative), spanning the range of Z depth used for
   specimens, at several field positions including the field corners, so that
   day-to-day reproducibility and field-dependent scale error are both
   estimable. The independent unit for uncertainty is the imaging session, not
   the field of view.
3. **Measurements per session**: lateral scale in X and Y separately (a single
   isotropic factor would hide an XY asymmetry), axial scale, field-dependent
   scale variation, and per-channel centroid offsets in X, Y and Z from the bead
   images.
4. **Analysis**: per stratum, report the estimated scale factor per axis with a
   confidence interval, its deviation from the metadata value, the
   between-session component of variance, and the registration offsets in µm
   and in voxels. Propagate the scale deviations into area and volume to state
   the calibration contribution to measurement error directly in the units the
   project reports.
5. **Acceptance rule, predeclared before any standard is imaged.** Illustrative
   only: per-axis scale deviation from metadata < 1% with the upper confidence
   limit also < 1%; between-session standard deviation < 0.5%;
   channel registration offset < half the lateral resolution limit. A single
   stratum failing is a FAIL for that stratum, and the qualified set of
   instruments is exactly those strata that pass. The owner must set these
   numbers from the measurement requirement, not from the first results.
6. **Record**: the study's outputs become a validation record under
   `docs/VALIDATION_RECORDS.md` — manifest, hashes, source commit, instrument
   identifiers, standard certificate numbers and expiry, and the raw
   measurements — so that SG-6's status is read from a record rather than
   asserted in prose.

### 2.3.1 Error budget, and what it implies about the acceptance numbers

The illustrative numbers in 2.3.5 can be replaced by derived ones once the
question is put the right way round: *how large may a calibration error be
before it dominates the estimator error the project has already qualified?*

The qualified in-domain worst cases on the confirmation set are 0.95% for area
and 0.76% for volume. Requiring the calibration contribution to stay below
those, using the isotropic sensitivities (`3e` for volume, `2e` for area):

| requirement | per-axis tolerance |
|---|---|
| calibration volume error < 0.76% | `e < 0.25%` |
| calibration area error < 0.95% | `e < 0.47%` |
| both | **`e < 0.25%` per axis** |

A 0.25% per-axis tolerance is demanding for the lateral axes and, with
refractive-index mismatch and depth-dependent aberration, is not plausibly
achievable for the axial axis in a confocal stack of a mounted organoid. The
honest consequence is not a looser criterion but a narrower claim:

1. **Absolute physical-unit accuracy at the estimator's precision is out of
   reach**, and a claim of the form "this organoid's volume is X µm³ to within
   1%" should not be made. If absolute units are reported, the calibration
   uncertainty must be reported as a separate systematic term alongside them,
   not folded into the estimator's error.
2. **The intended use does not require it.** IU-1 is comparison of conditions
   within one experiment. A scale error common to all objects in the experiment
   cancels in a ratio and shifts a difference of logs not at all; what does not
   cancel is *instability* — between-session and field-position variation, and
   any change of objective, medium or step size between arms.
3. Therefore the acceptance criteria should be **tiered by claim**, each tier
   predeclared and each admitting its own verdict:

| Tier | Claim it licenses | Criterion (owner to fix the numbers, from this budget) |
|---|---|---|
| **T1 — stability** | within-experiment comparison of conditions in physical units, one instrument configuration | between-session and field-position scale variation per axis below the tolerance derived above (`< 0.25%` if the estimator budget is to dominate); all arms of a contrast acquired in one configuration |
| **T2 — accuracy** | reporting absolute physical units | per-axis deviation from a traceable standard, with the upper confidence limit, below a declared tolerance; if the tolerance cannot be met, the claim becomes "absolute units with a stated systematic uncertainty of ±X%" |
| **T3 — transferability** | comparison across instruments or across long time gaps | T1 and T2 per stratum, plus an estimated cross-stratum bias with its interval |

T1 is achievable and is what the current intended use needs. T2 and T3 are
separable, and failing them does not invalidate T1 — but a report must then
state which tier it rests on. The distinction also determines what the axis-ratio
result means: because `rho_in` is invariant to isotropic scale error (§2.1.1),
domain membership is a T1-type question about the z:xy ratio's stability, while
absolute size is a T2 question.

### 2.4 What this study cannot do

- It cannot qualify specimens, only instrument configurations. Drift after the
  qualification date is uncontrolled unless the study is repeated; the record
  should therefore state a re-qualification interval, which is an owner
  decision.
- It cannot detect a specimen-induced axial distortion (mounting medium,
  clearing, depth-dependent aberration) beyond what the standard imaged in the
  same medium captures.
- It says nothing about whether the segmentation is correct (SG-3) or whether
  the estimators are accurate (SG-1/SG-2). All three are required; none
  substitutes for another.

### 2.5 Interaction with the spacing work in this round

The per-axis provenance added this round is a precondition for this study being
interpretable, not a substitute for it. Recording that an axis's spacing came
from a manifest rather than from the file makes it possible to ask whether the
manifest value was ever verified — but recording provenance is not verification,
and a manifest value entered by hand carries exactly the authority of whoever
typed it. SG-6 stays **FAIL / INSUFFICIENT EVIDENCE** on the strength of this
round's work.

## 3. SG-5 — statistical and inferential scope: proposed study

Added in round 4. SG-5 is **NOT ASSESSED**: no study has been performed on the
inferential layer. This section states what would qualify it, and — more
usefully in the short term — what the layer may be used for meanwhile.

### 3.1 What is implemented

- `statistics/aggregation.py` — per-unit descriptives: counts of detected,
  included and excluded objects, medians of the morphology metrics over
  morphology-eligible objects, total volume, and state counts and fractions.
  At condition level, the reported interval is a **percentile bootstrap of the
  mean of replicate-level median-of-unit-medians**, 2000 resamples on a pinned
  seed, computed only when at least 3 biological replicates carry size data and
  `NaN` otherwise (`aggregation.py:15, 78-79`; `config.py:51`).
- `statistics/inference.py` — per feature, `feature ~ condition` with a
  replicate random intercept (`mixedlm`, REML, grouped on `_replicate_key`;
  `inference.py:64, 163`), falling back to OLS with replicate-clustered
  standard errors and `use_t=True` when the random-effects fit is singular or
  fails to converge (`inference.py:71`). Pairwise contrasts carry a manual
  t(G−1) small-cluster reference in both branches; the LMM omnibus Wald test is
  asymptotic and the docstring already marks it as screening only. Contrast
  p-values are Benjamini–Hochberg adjusted across the family
  (`inference.py:170`). The method rationale and citations are recorded in
  `docs/ALGORITHM_DECISIONS.md` D11.
- `statistics/exploration.py` — normality and equal-variance checks, two-group
  comparisons, effect sizes, and supervised classifiers with fixed seeds. This
  module is exploratory by name and must stay out of any reported inference.

### 3.2 The three questions a qualification would have to answer

**Q1. Does the model's uncertainty statement have its nominal coverage under the
project's own replicate structure?** Method: simulation under known truth. Generate
object-level data with a specified condition effect, a between-replicate variance
component, a between-unit component, and realistic per-unit object counts
including the unbalanced and small-n cases the pipeline will actually meet
(1–3 replicates, wells with few objects). Run the implemented path end to end.
Criterion, predeclared: empirical coverage of the nominal 95% interval within a
declared tolerance across the simulated grid, and type-I error at the nominal
level under the null. Report the cells where coverage fails rather than a single
average — REML with two or three replicate levels is where degradation is
expected, and the result is likely to be a *domain of validity in terms of
replicate count*, which is exactly the useful output.

Three specific sub-questions follow from what is implemented, and the design
should target them rather than treating the path as a black box: (i) the
**branch-switching** behaviour — the simulation must record which branch fired
per run and report coverage per branch, because a variance component that
collapses near zero switches estimator mid-study and the two branches need not
have the same coverage; (ii) the **t(G−1) correction's adequacy** at G = 2 and
G = 3, which is the regime the correction exists for and the one least likely to
hold; (iii) the **omnibus test**, already documented as asymptotic — the
simulation should quantify how anti-conservative it is, so the "screening only"
label is backed by a number rather than by a caveat.

The **percentile bootstrap interval is the single most likely thing here to fail
coverage**, and the simulation must evaluate it separately from the model's
intervals. A percentile bootstrap of a mean over 3 values resamples from three
points: the interval cannot extend beyond the observed range, it carries no
small-sample correction, and it is known to undercover in this regime. The
project's own floor is exactly `n = 3`. Two outcomes are acceptable and both are
useful — a declared minimum replicate count below which the interval is not
reported, or a switch to a small-sample interval (t-based, or BCa) — but the
choice must follow the coverage numbers, and the interval must not be presented
as a 95% interval in any regime where the simulation shows it is not one.

**Q2. Is the unit of analysis handled correctly?** The measurement unit is the
object, the inference unit is the biological replicate. A simulation with a
strong within-replicate correlation and no true condition effect must not
produce significant contrasts; if it does, the object-level rows are being
treated as independent somewhere in the path. This is a specific, cheap,
falsifiable check and should be run first.

**Q3. Do the aggregation choices bias the contrast?** Median-of-objects then
mean-over-replicates is not the same estimator as a pooled mean, and the
difference depends on per-unit object count. Simulate skewed object-level
distributions with unequal counts and report the bias of each aggregation
against the known truth. This is characterisation: the output is a statement of
which aggregation the reported contrast refers to, not a PASS.

### 3.3 Scope statement usable before any of that exists

Until Q1–Q3 are answered, the inferential outputs are **descriptive summaries
with uncertainty indications of unverified coverage**. They may be used to
describe an experiment's observed differences and to plan a confirmatory
experiment. They may not be used as the primary evidence for a biological
conclusion, and a p-value from this path must not be reported without the
replicate count it rests on. `exploration.py`'s classifier outputs may not appear
in a results claim at all.

### 3.4 What this study cannot do

It cannot make an underpowered experiment interpretable, and it cannot validate
the measurements it consumes: a correctly covered interval around a biased
measurement is a precise statement of the wrong quantity. SG-5 is independent of
SG-1/SG-2/SG-3/SG-6 in method and cannot substitute for any of them.

## 4. Sequencing note

SG-6 is logically upstream of SG-3: a segmentation study run on uncalibrated
acquisitions cannot separate a scale error from a mask bias, because both appear
as a systematic volume discrepancy against the reference standard. If only one
study can be resourced, calibration first is the cheaper and more informative
order. This is a recommendation, not a result.
