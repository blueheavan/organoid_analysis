# Literature basis for the fourteen open decisions

Status: **decision support, not evidence.** This document creates no claim about
the pipeline, reports no measurement, runs no study, and changes no gate status.
Its only product is the set of recorded decisions in `docs/OWNER_DECISIONS.md`
and the derivations behind them.

Retrieved 2026-09-16 against live OpenAlex and Crossref. Every citation below
resolves to a registered DOI and was checked for retraction flags. **Access
depth: publisher abstract only; no full text was obtained for any entry** — the
access-depth rule this review asks for is stated in §7 and applies to this
document first.

## 1. The two questions behind all fourteen decisions

The register's fourteen items look heterogeneous — a cavity policy, an axial
tolerance, a degrees-of-freedom choice — but they are instances of two questions
that the measurement-literature answers directly, and this repository has been
answering ad hoc.

The first is *where a tolerance comes from*. Metrology's answer is that a
tolerance is derived from the decision the measurement must support, expressed
before the measurement is made, and never from what the instrument happens to
achieve ([Pendrill 2014](https://doi.org/10.1088/0026-1394/51/4/s206), *Metrologia*).
Clinical chemistry operationalised this into acceptance criteria computed from
biological variation rather than from instrument capability, which is why
`<1%` can be honestly labelled an engineering target while a *different* number
governs a biological contrast
([Ricós 1999](https://doi.org/10.1080/00365519950185229);
[CCLM 2024](https://doi.org/10.1515/cclm-2024-0108)). The same source gives the
rule that matters most here: when several error terms each contribute to one
claim, the budget is *shared*, not granted to each term in full.

The second is *what makes an agreement number mean anything*. Image-analysis
validation has converged on the position that a metric is chosen from the
claim and the data distribution, not from convention, and that the pitfalls are
predictable: class imbalance, undefined denominators, and criteria stated at a
weaker level than the claim they feed
([Maier-Hein 2024](https://doi.org/10.1038/s41592-023-02151-z), *Nature Methods*;
[Reinke 2024](https://doi.org/10.1038/s41592-023-02150-0), *Nature Methods*).
For classification agreement specifically, Cohen's kappa behaves incoherently
under imbalance — a worse classifier can score higher — so the confusion matrix
and class-specific sensitivity and precision are the primary report, with
macro-F1 and MCC as summaries
([Delgado & Tibau 2019](https://doi.org/10.1371/journal.pone.0222916);
[Chicco & Jurman 2020](https://doi.org/10.1186/s12864-019-6413-7)). Reporting
standards for such studies exist and are minimal-cost to follow
([Kottner 2011](https://doi.org/10.1016/j.jclinepi.2010.03.002), GRRAS).

Three of the register's items (D-1, D-2, D-8) are not tolerance questions at all
but *definition* questions — what is being measured, and against what. Those are
decided by the frozen estimand in `INTENDED_USE_AND_ESTIMANDS.md` §3 plus the
biology of the object.

## 2. Why a differential criterion, not an absolute one

The most transferable result from the measurement literature here concerns
which error terms actually threaten a claim. For a comparison between
conditions, an error term that is *common* to all objects, instruments and arms
cancels in the ratio and contributes nothing to a difference of log-volumes;
what does not cancel is an error term that is *differential* — one that varies
with condition, stratum, size or depth — and instability, which is a
differential term by construction.

`STUDY_PROPOSALS_SG3_SG6.md` §2.3.1 already applies this argument to
calibration, and the review's M5 measured it for the segmentation metric
(`reviews/spacing_invariance_check.py`: the volume discrepancy is bit-identical
across a common ×2 spacing change). The consequence for the register is that
D-5, D-6, D-7 and D-12 should each be stated twice: once as an absolute bound
that protects single-object statements, and once as a *differential* bound that
protects the between-condition contrast the intended use actually makes. A
criterion that only bounds the common part qualifies a number nobody reports.

## 3. The object: cavities, heterogeneity, penetration

`INTENDED_USE_AND_ESTIMANDS.md` §3.3 already identifies the open-cavity case as
the defect risk: a lumen connected to the border is not filled, so `V_env`
collapses to `V_seg` and the measured quantity changes discontinuously with
topology, with no biological change.

The biology says a lumen is a *structure*, not an artefact. Lumen formation is
the mechanical driver of folding in cerebral organoids
([Karzbrun 2018](https://doi.org/10.1038/s41567-018-0046-7), *Nature Physics*),
and luminal and cystic architectures are normal in the model systems the
intended use names ([Organoids, *Nat Rev Methods Primers* 2022](https://doi.org/10.1038/s43586-022-00174-y)).
Deleting such objects, or silently closing them, therefore changes the
*population* rather than cleaning it — and sample definition is the classic
source of bias in volume estimation
([Gundersen & Jensen 1987](https://doi.org/10.1111/j.1365-2818.1987.tb02837.x), *Journal of Microscopy*).
This points to the third of the three options the register offers: keep every
object in the report and define the estimand per topology, rather than exclude.

Heterogeneity is the other reason not to use object counts as evidence.
Organoid morphometry platforms exist precisely because size and shape diversity
within a culture is otherwise ignored — `OrganoSeg` stratified 5 167 spheroids
and 5 743 organoids to make its point
([Borten 2018](https://doi.org/10.1038/s41598-017-18815-8)), and `OrganoID`
frames single-organoid tracking as the response readout
([Matthews 2022](https://doi.org/10.1371/journal.pcbi.1010584)). A purpose-built homozygous
high-throughput workflow had to advertise *unusually homogeneous* midbrain
organoids to be screening-compatible
([Nickels 2020](https://doi.org/10.7554/elife.52904)) — the exception that
demonstrates the rule.

Penetration sets a stratification that is not optional. Oxygen tension falls
with spheroid diameter and produces a hypoxic core
([McKeown 2017](https://doi.org/10.1098/rsif.2016.0851)), which is why necrotic
cores are a standing caveat for 3D models
([Nunes 2020](https://doi.org/10.3390/pharmaceutics12121186)). Marker
penetration fails with depth for the same reason, so the viability states and
the segmentation bias are both *depth- and size-dependent*: any ceiling on an
indeterminate rate is a statement about the stratified population, not about
the whole.

## 4. Absolute size: the axial axis is a different problem

The register's D-12 and D-13 sit on a result the repository already derived but
has not yet priced: the per-axis tolerance is `< 0.25%`, and the axial axis
cannot meet it in a mounted confocal stack.

The mechanism is documented and quantified. Axial distances in confocal
microscopy are distorted by refractive-index mismatch, and measured axial
scaling factors depend almost linearly on the refractive index of the immersion
and mounting media ([van Elburg 2014](https://doi.org/10.1111/jmi.12194),
*Journal of Microscopy*). Sample-induced spherical aberration is the dominant
depth-dependent term, and the correction is a per-configuration calibration,
not a property of the instrument
([Diel 2020](https://doi.org/10.1038/s41596-020-0360-2), *Nature Protocols*).
Lateral calibration is a different regime: it is traceable with a stage or grid
standard, and the failure modes are the ordinary ones of acquisition
parameters and precision ([Waters 2009](https://doi.org/10.1083/jcb.200903097),
*Journal of Cell Biology*).

Two consequences follow. Lateral and axial tolerances must be declared
separately rather than as a single "per-axis" number, and the tolerance itself
must be split between *stability* and *accuracy*, because those are the two
independent contributions to one budget. The metrology framework for combining
them is the uncertainty-propagation rule — independent contributions combine in
quadrature, not additively, and each stated limit carries a confidence level
([JCGM 100:2008](https://doi.org/10.59161/jcgm100-2008e);
[NIST TN 1297](https://doi.org/10.6028/nist.tn.1297)) — which is exactly the
defect review M1 found: the same `0.25%` was granted to stability and to
accuracy independently, and two terms at `0.25%` do not sum to a `0.25%` total.

## 5. Inference with three to twenty clusters

The shipped path takes REML `mixedlm`, falls back to OLS with
`cov_type="cluster"`, and applies a manual `t(G−1)` reference to pairwise
contrasts (`statistics/inference.py:61-121`). Every element of that choice is
known to be anti-conservative in the regime the intended use occupies —
"≥ 3 replicates per stratum" is the design the review's B3 and the report's
own Gap 4 both flag.

For mixed models, the small-sample problem has a standard correction: Kenward
& Roger's
adjustment to the fixed-effects variance and degrees of freedom
([Kenward & Roger 1997](https://doi.org/10.2307/2533558), *Biometrics*), and
empirical evaluations of significance testing in LMMs recommend the Satterthwaite
or Kenward–Roger approximations over a naive reference distribution
([Luke 2016](https://doi.org/10.3758/s13428-016-0809-y), *Behavior Research Methods*).

For the cluster-robust branch, the problem is the same and the fix is named:
cluster-robust variance estimation is biased downward when the number of
clusters is small, and the remedy is bias-reduced linearization (CR2) with a
Satterthwaite approximation to the t statistic
([Pustejovsky & Tipton 2016](https://doi.org/10.1080/07350015.2016.1247004),
*Journal of Business and Economic Statistics*; reference implementation
[clubSandwich](https://doi.org/10.32614/cran.package.clubsandwich)). A practical
guide to cluster-robust practice reaches the same conclusion for empirical work
([MacKinnon, Nielsen & Webb 2022](https://doi.org/10.1016/j.jeconom.2022.04.001),
*Journal of Econometrics*), and the specific few-cluster case in group-randomised
trials has been evaluated directly
([Huang & Li 2021](https://doi.org/10.3758/s13428-021-01627-0)).

That is why the register cannot be resolved by picking the cheaper of
Satterthwaite and between-within: both are *degrees-of-freedom* choices, and
neither corrects the variance estimate that the `t(G−1)` path leaves biased.
Between-within is defined for the omnibus test, where the design is balanced and
nested; it is not a per-contrast approximation, and this repository's claim is a
per-contrast claim.

Finally, the method that will be frozen cannot be qualified by assertion. A
simulation study is the instrument for that, and its design has its own standard:
define aims, data-generating mechanisms, estimands, methods and performance
measures; report the Monte-Carlo standard error of every performance measure;
and set tolerances that are justified against it
([Morris, White & Crowther 2019](https://doi.org/10.1002/sim.8086), *Statistics in Medicine*).
Since `MCSE = sqrt(p(1−p)/n_sim)`, a tolerance of one percentage point on
nominal 95% coverage requires roughly five thousand replicates per cell before
it is even meaningful.

## 6. Viability: which reference can carry which claim

The four-state readout is a classification, so its reference must be a
classification of the same object by a *measurement path* independent of the
production branch. The chemistry being independent is not sufficient — that is
the circularity B2 identified.

The assay literature sets the ceiling on what the ATP route can do. ATP-type
viability assays measure metabolic activity of a population, are well-level by
construction, and diverge from dye-based viability in 3D culture precisely
because penetration and metabolic state vary with depth
([Knobel 2024](https://doi.org/10.3233/ch-248101)). They are the standard
secondary well-level reference and cannot establish object identity. The
object-level reference therefore has to be human or an independent tool; a
second segmentation tool is only independent if its measurement path does not
use the production nuclear branch, and the only production nuclear segmentation
available is R-13, which is `NOT ASSESSED` and exploratory-only.

Sample-size and agreement rules for such a reference are also standardised:
report reliability and agreement per the reporting guideline
([Kottner 2011](https://doi.org/10.1016/j.jclinepi.2010.03.002)); size the study
from the required precision rather than from object-count precedent
([Shoukri 2004](https://doi.org/10.1191/0962280204sm365ra);
[Rotondi 2013](https://doi.org/10.3758/s13428-013-0415-1)); and for a
two-method agreement claim, size for the confidence interval of the limits of
agreement against a predefined clinical limit
([Lu 2016](https://doi.org/10.1515/ijb-2015-0039)), reporting bias and limits
rather than a correlation
([Giavarina 2015](https://doi.org/10.11613/bm.2015.015)). Where a criterion is
expressed as a proportion with a confidence bound, use an interval method that
behaves near the boundary rather than the Wald interval
([Brown, Cai & DasGupta 2001](https://doi.org/10.1214/ss/1009213286)).

One further point governs D-7. The indeterminate fraction is not missing data to
be imputed, but the literature's warning still applies in its useful direction:
the *proportion* is not itself a validity criterion, so a bare ceiling on it can
be satisfied while the classifiable subset is unrepresentative
([Madley-Dowd 2019](https://doi.org/10.1016/j.jclinepi.2019.02.016)). Since
`f_viable_like` is computed on the classifiable denominator, the criterion that
actually protects the claim is that the indeterminate fraction must not differ
*between the conditions being compared*.

## 7. Access depth, and the rule this review asked for

The review's B1 asked for a bibliography with a stated access-depth rule. The
rule adopted here, and it applies to every entry above:

| Depth | Meaning in this repository | May support |
|---|---|---|
| `full-text` | the complete paper was obtained and read | a numeric value attributed to the paper |
| `abstract` | the publisher abstract was retrieved; full text not obtained | a mechanism, a framework, a method name, a qualitative finding |
| `citation-only` | known only from another source's reference list | nothing; listed and marked unresolved |

Every citation in this document is `abstract`. **No numeric acceptance value in
§8 is taken from these papers as a measured quantity.** The papers supply the
framework, the mechanism and the method; every number is derived inside this
repository from the frozen intended use, the qualified in-domain error budget
already recorded in `evidence/2026-09-13-volume-audit/STUDY_PROPOSALS_SG3_SG6.md`
§2.3.1, and the stated rule that each error term receives a share of one budget.

## 8. Recorded decisions

Each entry states the decision, the derivation, and the condition that reverses
it. All fourteen were recorded 2026-09-16 under the project owner's delegation
to the analysis agent; they are predeclared before any confirmation or
annotation data are inspected, and each may be superseded only by a dated
entry in `docs/OWNER_DECISIONS.md` recorded before the data its study requires
are seen.

### D-1 — Open-cavity behaviour: define per topology, exclude nothing

**Decision: neither "exclude" nor "flag-only" alone — define the estimand for
every object and report the flag as a stratum.**

`V_env` stays the primary estimand for every object. The cavity policy is
carried by the companions, not by object selection:

- `V_seg` (segmented material) and `f_void = (V_env − V_seg)/V_env` are reported
  for every object, as §3 already requires;
- an object whose void is topologically connected to the exterior background
  carries the flag `open_cavity`, and its `f_void` is reported as
  `0 (by construction)` with the flag, never as a measurement of the lumen;
- no object is dropped. `open_cavity` is a **declared stratum**, and every
  qualification or agreement result is reported with and without it.

Rationale: the discontinuity §3.3 identifies is a defect of *reporting* rather
than of the estimand — the estimand is the outer envelope, which the cavity does
not perturb — and removing objects would change the population in a
size-correlated way (§3). The flag makes the topology visible where it changes
what a reader may infer.

Reversal: if a future development study shows the open/closed classification is
itself unstable between time points of the same organoid, the policy narrows to
reporting `V_seg` as primary for that stratum, recorded with the study.

### D-2 — Volume applicability domain: the surface domain, plus the open-cavity stratum

**Decision: the volume domain is the surface domain — smooth/closed, `rho_in ≥ 10`,
anisotropy ≤ 4 (`SCIENTIFIC_SPEC.md` §9) — extended only by the `open_cavity`
stratum of D-1. Creased and hollow objects are out of domain for volume
qualification.**

The confirmation evidence forces this: smooth in-scope cases meet `<1%` (worst
0.76% for volume) while creased in-gate cases do not (3.20% volume,
`INTENDED_USE_AND_ESTIMANDS.md` §2). A second, looser volume domain would need
its own development and confirmation set before it could be used
(`VOLUME_QUALIFICATION_PROTOCOL.md` §5 item 1), and the honest move — which the
errata and §2.3.1 both take — is a narrower claim, not a looser criterion.

Reversal: a project owns the domain it confirms. If a development study
qualifies a creased volume domain with its own confirmation set, the domain
widens by that record; it cannot widen by argument.

### D-3 — Resource analytical volume qualification now

**Decision: resource it now, in parallel with SG-3, as the first item of
Track A.**

The work needs no imaging data (`VOLUME_QUALIFICATION_PROTOCOL.md` §5 item 2),
and the master plan states the dependency the other way round from the register's
framing: analytical qualification is a *precondition* for SG-3's downstream-bias
analysis, not a substitute for it. Sequencing it behind SG-3 would place a
zero-specimen, zero-annotation task on the critical path of the
annotation-limited study.

### D-4 — `<1%` is an inherited engineering target and is retained

**Decision: `<1%` is not a biological requirement; it is retained unchanged as
the SG-2a qualification criterion, and it is explicitly not the criterion for
SG-3.**

Two separate statements, both needed. Loosening it is refused for a specific
reason rather than inertia: the register permits loosening only before a
confirmation set is measured and with a recorded biological justification, and no
biological requirement looser than `<1%` has been derived — the intended use
(IU-1, comparison of conditions within an experiment) does not require absolute
single-object accuracy at all, which is why the *absolute* claim is bounded by
T2 in D-12 and the *contrast* claim is governed by the differential criteria of
D-6. Retaining `<1%` costs nothing that the claim needs and keeps the estimator's
qualification conservative.

### D-5 — SG-3 δ, acceptance and N

**Decision: predeclare δ on the readout scale, state each criterion per case,
and derive N from precision.**

Smallest biologically relevant effect:

| readout | δ |
|---|---|
| `V_env`, `V_seg` (volume) | 0.20 relative |
| surface area | 0.20 relative |
| equivalent diameter | 0.10 relative |
| sphericity | 0.10 relative |
| principal-axis length | 0.10 relative |

δ is a **predeclared convention**, not a measured value: it is the smallest
between-condition change this project treats as biologically actionable, set
conservatively high because organoid heterogeneity within a culture is the
operative limit on resolution (OrganoSeg's 5 167-spheroid stratification;
§3 above). It is not derived from this repository's observations.

Acceptance, applied per matched object with uncertainty at the biological-sample
level (cluster bootstrap):

| ID | criterion |
|---|---|
| C1 | lower 95% CI bound on recall ≥ 0.90 **and** on precision ≥ 0.90, per stratum carrying the claim |
| C2 | upper one-sided 95% bound on the median absolute relative bias < δ/4 (5% for volume) |
| C3 | upper one-sided 95% bound on the 95th percentile of absolute relative bias < δ/2 (10%) |
| C4 | upper 95% bound on the rate of catastrophic error (absolute relative bias > 0.5) < 0.02 |
| C5 | the bias distribution's 95th percentile bound is the limits-of-agreement criterion; bias and LOA are reported with intervals, never as point estimates |

Minimum evaluable N, derived from C1's precision: a lower bound of 0.90 to
±0.05 needs about 140 objects before clustering; with a predeclared design
effect of `1 + (m − 1)·ICC = 1 + 19 × 0.30 = 6.7` at `m = 20` objects per
biological sample, about 940 objects, i.e. **≥ 47 biological samples per
stratum**. At the level of the contrast, **≥ 20 biological samples per compared
condition**, with the bias criterion resolvable to ±2.5% (`δ/8`) at an assumed
between-sample SD of 0.05. If the achieved CI is wider than the criterion width,
the verdict is `INSUFFICIENT EVIDENCE` for that readout — never PASS on a point
estimate, which the protocol already states.

### D-6 — SG-3 downstream-bias threshold, on SG-3's own terms

**Decision: the inherited `1% / 5%` pair is not SG-3's criterion. SG-3 is
governed by D-5's C1–C5, whose bias criterion is differential as well as
absolute.**

C2–C4 are the absolute bounds. The criterion that protects the intended use is
the differential one: **for any two conditions compared, the 95% CI of the
difference in median relative bias must exclude ±δ/8 (±2.5%), and the same
applies to surface area at ±δ/8.** A common bias cancels in a ratio (§2, and
M5's bit-identical measurement); a bias that differs between arms does not, and
no absolute per-case threshold would catch it.

SG-1a/SG-2a keep `1% / 5%`, because those items qualify the *estimator on
analytical cases*, not the production segmentation against annotation. The two
numbers answer different questions, and the register's M6 concern is resolved by
saying so rather than by re-deriving the analytical pair.

### D-7 — SG-4B criteria and the indeterminate ceiling

**Decision: per-class lower confidence bounds, the confusion matrix as primary,
kappa non-gating, and a two-part indeterminate ceiling.**

- Per-class sensitivity and precision: **lower 95% CI bound ≥ 0.80** for each of
  `viable_like`, `mixed_signal`, `compromised_like`, `indeterminate`. Macro-F1
  lower bound ≥ 0.80.
- Quadratic-weighted Cohen's kappa and Spearman correlation are reported and do
  not gate, for the reason §1 gives.
- A class with prevalence < 5% in the evaluated population is reported as
  **not evaluable**, not as FAIL — a bound cannot be estimated from too few
  positives, and the pitfall is documented.
- Indeterminate ceiling: **`f_indeterminate ≤ 0.20` per condition**, and
  **the between-condition difference in `f_indeterminate` must not exceed
  5 percentage points on its upper 95% bound**. Above either, V2 is not
  reportable for that condition.

The differential half is the load-bearing one: `f_viable_like` divides by
`N_classifiable`, so an indeterminate rate that differs between arms biases the
contrast even when both arms are under the ceiling (§6).

### D-8 — SG-4B primary reference: an independent blinded counting path

**Decision: exit (1). The primary object-level reference is a blinded human
classification of the same objects on a predeclared stratified subsample, by two
annotators, with the reference's own disagreement reported. The ATP assay remains
the secondary well-level reference. R-13 stays out of scope.**

Rationale: it is the only exit that removes the circularity without
re-qualifying the production nuclear branch, which is what V1 already depends on
(§6). Exit (3) would have narrowed the achievable claim to V2 and left the
per-object states permanently `INSUFFICIENT EVIDENCE`. Exit (2) doubles the
qualification burden in the direction of the defect.

Consequences to be recorded with the study, not discovered after it:

- the annotation set is the same asset SG-3 requires, so the two studies share
  one annotation campaign — this is the cost argument for the decision;
- the claim V1 licenses becomes *agreement with blinded human classification of
  calcein/PI signal*, which is bounded by the reference's own accuracy. The
  chemistry's specificity is established separately by the graded-insult
  monotonicity series;
- if inter-annotator disagreement exceeds D-7's per-class bounds, the outcome is
  `INSUFFICIENT EVIDENCE` for that class, not FAIL for the pipeline (the
  protocol already forbids demanding better agreement than the raters achieve).

### D-9 — SG-5 primary method: Satterthwaite; between-within for the omnibus only

**Decision: Satterthwaite degrees of freedom for the primary fixed-effect
contrasts; between-within is used only for the omnibus test.**

The primary claim is a per-contrast claim (a specific pairwise condition
difference), and Satterthwaite is a per-contrast approximation, whereas
between-within is an omnibus-test df defined for a balanced nested design
(§5). The omnibus Wald test keeps its asymptotic form but is reported with a
between-within df check and an explicit anti-conservative-at-`G < 10` statement,
because it is not the claim-bearing quantity.

Kenward–Roger remains a sensitivity analysis where a verified implementation is
available; it is not the primary method because the same verification burden
buys nothing the Satterthwaite choice does not already carry for a per-contrast
claim.

### D-10 — SG-5 coverage tolerance and reportable replicate counts

**Decision: tolerance of one percentage point on nominal 95% coverage, against a
simulation MCSE of at most 0.33 percentage points.**

`n_sim = 5 000` per cell gives `MCSE = sqrt(0.95 × 0.05 / 5 000) = 0.31` pp, so
the ±1.0 pp tolerance is more than three Monte-Carlo standard errors wide and is
therefore a statement about the method rather than about the simulation
(§5). The same rule applies to T1 type-I error at nominal 0.05: the estimate must
lie in [0.04, 0.06].

Reportable replicate counts: an interval may be reported only for cells that
pass T2 coverage at the frozen tolerance. The shipped floor of 3 replicates
(`statistics/aggregation.py:15`) is retained as the computation floor but is
**provisional**: at 2 and 3 replicates a 95% interval from 2 000 percentile
resamples has no coverage that this study has yet certified, so any interval
reported there is labelled `provisional — coverage unqualified` and may not
carry a claim. The study may raise the reportable floor; it may not lower it.

### D-11 — SG-5 CR2 fallback: implement locally, verify against a published implementation

**Decision: commission the local implementation and verify it against the
published reference implementation and published known cases.**

"Accept a reliable external implementation" is not available to this project in
its pinned environment: `statsmodels` 0.14.6 provides `cov_type="cluster"`
(CR0-style) and no bias-reduced linearization, and no Python CR2 implementation
is in the lock. The choice is therefore between unverified local code and
verified local code, and the latter is the only one compatible with a coverage
claim.

The verification protocol, recorded as part of the freeze:

1. golden values from `clubSandwich` (recorded with its version) on a fixed set
   of synthetic designs, stored as fixtures with their provenance;
2. analytic known cases where CR2 must equal CR0 or the model-based estimator
   (saturated and balanced designs);
3. the cases, degrees-of-freedom values and intervals published in
   [Pustejovsky & Tipton 2016](https://doi.org/10.1080/07350015.2016.1247004);
4. the D-10 simulation study, run on the frozen implementation only.

### D-12 — SG-6 tolerances: one budget, split by axis and by tier

**Decision: the `0.25%` per-axis figure is the total calibration contribution to
the volume error and is split, not reused.**

The qualified in-domain worst case is 0.76% for volume, and volume responds to
per-axis scale error as the sum of the affected axes
(`STUDY_PROPOSALS_SG3_SG6.md` §2.3.1): `e_x + e_y + e_z ≤ 0.76%`. Allocating
equally gives the `0.25%` per axis already derived. That single number then
splits twice:

| tier | quantity | lateral (x, y) | axial (z) |
|---|---|---|---|
| **T1 stability** | between-session and field-position scale variation, gated on the upper 95% bound | ≤ 0.10% | ratio stability ≤ 0.25% |
| **T2 accuracy** | deviation from a traceable standard, gated on the upper 95% bound | ≤ 0.15% | **not qualified**; reported as a stated systematic uncertainty ±X% |
| **combined per axis** | quadrature sum (`JCGM 100:2008`) | `sqrt(0.10² + 0.15²) ≈ 0.18%` | stability only |
| **T3 transferability** | cross-instrument/cross-gap bias, derived from the contrast budget `δ/8 = 2.5%` | ≤ 0.8% per axis | reported, not gated |

Three points of substance. Lateral and axial are separated because the axial
term is dominated by a mechanism the lateral term does not have (§4) and because
pretending otherwise is what produced the M1 defect. T1 lateral is 0.10% rather
than 0.15% because the stability term is the one that does *not* cancel across
arms and therefore deserves the larger share of the quadrature budget. T2 axial
is refused rather than loosened: a claim of the form "this organoid's volume is
X µm³ to within 1%" is not made, and absolute axial units are reported with the
measured systematic term stated alongside (§2.3.1 items 1–2).

The ratio-stability requirement of 0.25% is not tight in relation to what it
protects: the smallest relative gap in the eleven predeclared z:xy strata is
4.17% (3 → 3.125), so the tolerance has more than sixteenfold margin.

### D-13 — Re-qualification interval: 6 months, then 12, with event triggers

**Decision: an initial interval of 6 months from the qualification date,
extending to 12 months if the first interval shows drift within T1, with
event-triggered re-verification and ongoing control-chart monitoring.**

Interval-setting in metrology is an inference from observed drift and
reliability rather than a fixed convention, with short initial intervals that
lengthen once drift is characterised
([Yang & Tsai 2005](https://doi.org/10.1109/tim.2004.840234), *IEEE TIM*;
[Pendrill 2014](https://doi.org/10.1088/0026-1394/51/4/s206)). A 6-month first
interval is the conservative implementation of that rule for an instrument whose
axial behaviour is configuration-dependent (§4).

Event triggers force re-verification of T1 before further physical-unit claims:
any change of objective, immersion or mounting medium, Z step, camera, filter
cube or acquisition software version; any instrument move; and any point outside
the control limits of the ongoing monitoring. Monitoring is a standard specimen
measured at a predeclared frequency, plotted against the qualified T1 value with
±2·`e_stab` limits, and its purpose is to detect drift between formal
re-qualifications rather than to replace them.

The interval is a property of the *instrument configuration*, not of the
specimens (as §2.4 already states), and the drift observed in the first interval
governs the second.

### D-14 — Option C record contract: authorised

**Decision: the Option-C domain-restricted qualification record contract is
authorised, with the interface fixed as follows, and its implementation is
recorded as blocked only by commit authority.**

The contract, so that implementation is mechanical:

1. the manifest gains a `record_class` of `canonical` or `scope-clean` and a
   per-item `results` map keyed by gate item id (`SG-1a`, `SG-1b`, `SG-2a`,
   `SG-2b`) carrying `status`, `domain` and the confirmation-evidence hashes;
2. `SG-1a` and `SG-2a` results may report PASS only against domain-restricted
   confirmation evidence — the domain constants (`rho_in ≥ 10`, anisotropy ≤ 4,
   the seed and artifact hashes) are hash-bound into the record;
3. `SG-1b` and `SG-2b` carry the unrestricted characterization with
   `counts_toward_gate = false`, so scope-restriction and scope-behaviour are
   both visible and neither is inferable from the other's absence;
4. the change lands as one commit satisfying `GATE_SEMANTICS_ANALYSIS.md` §6:
   the contract, the adversarial tests, and that decision record's citation.

Two things this decision does not do. It does not move any gate status — the
gate keeps its legacy fail-closed behaviour until the commit lands. And it
cannot be executed in this pass, because `VV_TOOLING` hash-binds
`analytical_geometry_evidence.py`, so any edit to the contract invalidates the
current analytical-geometry record until a new one is generated against a clean
committed tree, which requires commit authority this agent does not hold.

## 9. What this document leaves unresolved

- No gate status, estimator, threshold in code, default, domain constant or
  qualification claim is changed by anything in §8. The gate's statuses are
  identical before and after.
- The decisions are predeclared criteria and study designs. They are not
  evidence, and none of the studies they govern has been executed.
- Every one of the D-1…D-14 reversals is available to the project owner before
  the relevant data are seen, by a dated superseding entry in
  `docs/OWNER_DECISIONS.md`.
