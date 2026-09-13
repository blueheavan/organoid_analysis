# Review of `report_organoid_validation_solutions.md`

**Reviewed:** 2026-09-13 · **Subject:** externally supplied solution designs for
R-1/SG-6, R-6/SG-3, R-10/SG-4, R-5+R-9/SG-5 · **Reviewer action:** evaluation
only. Nothing in the report was executed; no study was run, no code changed, no
status moved.

## Verdict

**Accept the structure, do not file as-is.** The report is the right *kind* of
document — it proposes designs and software changes, states plainly that it is
not validation evidence, marks its numeric defaults as owner decisions, and does
not attempt to move any status. Its scope mapping onto the roadmap is correct
and its description of what this repository already specifies is accurate
wherever I could check it against the source (§1).

Four defects block filing it as a study design, and seven weaken the designs
themselves. Two of the blocking or major defects are inherited from *this*
repository's own documents rather than invented by the report — those are errata
candidates here (§6), not criticisms of the author.

The report's most consequential error runs in the opposite direction from the
ones this project has been guarding against: it makes the critical path **longer
than the evidence requires**, by treating calibration as a prerequisite for the
annotation study when that study's own metric is invariant to what calibration
measures (M5). Fixing it moves the expensive study earlier.

## 1. Claims about this repository — checked against source

| Report claim | Verdict |
|---|---|
| scope mapping R-1/SG-6, R-6/SG-3, R-10/SG-4, R-5+R-9/SG-5 | correct (`ROADMAP.md` §2) |
| out of scope R-2/R-3/R-7/R-11/R-12/R-13 | consistent; it silently *includes* R-4 and R-8, which is fine but should be stated |
| 3% lateral → +6.09% volume; isotropic 3% → +9.27% | correct (`STUDY_PROPOSALS_SG3_SG6.md:127,131`) |
| derived per-axis tolerance `< 0.25%` | correct (§2.3.1, line 237) |
| `rho_in` domain membership invariant to isotropic scale error | correct |
| estimator worst in-domain area error 0.95% | correct (0.951%, n = 96 in-scope smooth confirmation cases) |
| viability §5.5 Q1 preference order, Q2 gate grid, Q3, Q4 | faithful, including "characterisation, never gate selection" |
| gate grid `(0.25,0.55) (0.30,0.60) (0.35,0.65)` | correct, and the middle pair is the shipped default (`config.py:48-49`) |
| inference path: REML LMM → OLS-clustered fallback with `use_t`, manual t(G−1) contrasts, BH, percentile bootstrap at ≥ 3 replicates | correct in every particular |
| `segmentation_metrics.py` is Hungarian-matched IoU/Dice against a single reference, no consensus fusion | correct (`linear_sum_assignment`; IoU/Dice/F1/panoptic quality; no fusion) |
| `configs/imported_instances.yaml`, `microscopy_io/tiff_contract.py`, `exploration.py` exist | all present |
| overall gate 0/6; all items FAIL / INSUFFICIENT EVIDENCE / NOT ASSESSED | correct |

No misattribution to this repository was found. That is worth stating
explicitly, because it is the failure mode round 4 spent most of its effort on.

## 2. Blocking defects

**B1 — No bibliography.** The report carries more than forty bracketed citation
keys and no reference list. Every external anchor — the axial focus-shift
figure, the 594-nucleus benchmark, the ~6-cluster type-I result, the
percentile-bootstrap coverage bounds — is unresolvable as written, so I could
verify none of them. A document whose whole contribution is "the literature
supplies the missing anchors" cannot be filed in an evidence directory in this
state; this repository issued an errata this same round for citing its *own*
files inaccurately. The fix is mechanical: attach the reference list with DOIs
and mark which entries were read in full versus abstract-only (the limitations
section claims four full-text reads — name them).

**B2 — The primary reference standard for the viability claim is circular.**
Study V1's reference (a) is a nuclear dead stain "read by the pipeline's own
nuclear segmentation". Nuclear segmentation is R-13: **NOT ASSESSED**,
exploratory-only, and explicitly out of scope in this report. The per-object
claim would then be validated against a measurement path with no evidence and no
scheduled study. The chemistry is independent; the measurement is not. Three
acceptable exits: count dead nuclei by a path independent of production
(annotator counts on a subsample, or a second tool); bring R-13 into scope for
the modality used; or demote (a) and let well-level ATP be primary, accepting
that only V2 gets validated and the per-object states stay INSUFFICIENT
EVIDENCE. **This defect originates in `INTENDED_USE_AND_ESTIMANDS.md` §5.5**
(§6).

**B3 — Acceptance criteria sit at a weaker level than the claims they feed, and
no study has a precision calculation.** Every qualification rule in this
repository is per-case: the surface claim is `<5%` *for every case* in domain,
the volume claim `<1%` *for every case*. Study S1's acceptances are
central-tendency — "median IoU ≥ 0.70", "median relative discrepancy within 1%
volume / 5% area". A median-based PASS cannot support a per-case claim, and
chaining the two yields a statement stronger than its weakest link. Either state
S1's criteria per case (or as an upper quantile / one-sided upper confidence
limit on the discrepancy), or declare in advance that S1 qualifies central
tendency only and record that in the claim register. Separately: all four
studies set thresholds and none computes the sample size needed to *resolve*
them. With the biological replicate as the inference unit and ≥3 replicates per
stratum, S1's effective n per stratum is 3 — the same few-cluster regime the
report's own Gap 4 identifies as unreliable. Sample sizes should follow from the
precision required on the bias estimate, not from precedent object counts.

**B4 — Gap 4's plan invalidates its own first study.** Study T1 (the null
unit-of-inference check) is scheduled first, against the shipped path. The
method upgrades are then to be implemented "before the grid is frozen" —
including replacing the manual t(G−1) fallback with CR2 plus Satterthwaite
degrees of freedom. That is the branch T1 exercises, so T1 would have to be
re-run after the upgrade, and the report does not say so. Also unstated: the
upgrade changes the object of qualification from "the shipped pipeline" to "a
replacement", it needs its own contract tests and a re-freeze, and `statsmodels`
provides no CR2 — a hand-written sandwich estimator is new code whose
correctness must be verified against a published implementation before it can
carry a coverage claim. Either qualify the shipped path as-is and let the
findings justify replacing it, or sequence the upgrade first and drop T1 from
the pre-upgrade slot.

## 3. Methodological defects

**M1 — The calibration error budget is spent twice.** `0.25%` per axis was
derived as the *total* tolerance for the calibration contribution to stay inside
the estimator's in-domain budget. The report applies it as the T1
between-session CV acceptance *and* leaves T2 accuracy at ≤1% lateral. Two
independent terms at 0.25% do not sum to a 0.25% total. Allocate the budget
across stability and accuracy explicitly, or state that 0.25% governs the
combined term and derive each component from it. T1 also loses the "report the
upper confidence limit" requirement the report correctly keeps for T2 — a CV
point estimate from five sessions has a wide upper limit, and it is the upper
limit that must clear the tolerance.

**M2 — `ρ ≥ 0.8` is the wrong statistic for the viability claim.** The claim
under test is a four-state *per-object classification*. A rank correlation
between a continuous reference fraction and a state ordering can be satisfied
while a large share of objects are individually mis-stated, and it yields no
per-state error rate. The right instruments are a confusion matrix against the
reference (per-state sensitivity and specificity with intervals) plus an ordinal
agreement statistic such as quadratically weighted kappa; keep a correlation for
the well-level V2 comparison, where the estimand really is a continuous
fraction. The cited anchor (r > 0.90 in an imaging-versus-ATP comparison) is a
different estimand than per-object state agreement and does not calibrate this
threshold.

**M3 — The bootstrap remedy does not follow from the cited evidence.** The
report's own citation for small-sample under-coverage is described as covering
percentile *and* BCa intervals; the proposed fix is to switch to BCa at n ≥ 5.
If BCa also under-covers in that regime, the fix is unsupported. Consistent with
this project's rules: let T2 measure coverage for each candidate construction
(t-interval on replicate values, BCa, bootstrap-t) at each replicate count, and
predeclare that no interval is reported in any regime where measured coverage
misses the tolerance. Also unpriced: raising the floor from 3 to 5 replicates
removes the interval from most experiments this pipeline has processed.

**M4 — Precedent numbers are transferred across estimands.** The "≥500 nuclei"
default is anchored to a 594-*nucleus* benchmark and the "F1 > 0.5 satisfactory
floor" comes from the same source, but the claims needing support here are
organoid-level (`V_env`, area, sphericity) and nuclear metrics are out of scope
in this very report. The "≤8% mis-segmented" precedent is likewise a *cell*-level
criterion (boundary area deviation > 30%) applied to organoid strata. Keep the
organoid-level anchors (≥30 objects per stratum, ≥3 replicates), derive the rest
from B3's precision calculation, and drop the nuclear numbers or move them to an
R-13 study.

**M5 — The sequencing is wrong, and wrong in the expensive direction.** The
report requires calibration C1 before annotation study S1, on the grounds that a
scale error and a mask bias are indistinguishable in a volume discrepancy. For
the metric S1 actually reports, that is false, and the size of the residual is
measurable rather than arguable. Production and reference masks are measured on
the *same* voxel grid, so with `V = s_z s_y s_x N` the relative discrepancy
`(V_prod − V_ref)/V_ref = (N_prod − N_ref)/N_ref` is exactly independent of the
spacings — for **any** diagonal spacing error, not merely an isotropic one,
because the spacing product is a common factor of numerator and denominator.
Area ratios scale as `λ²` above and below and cancel likewise; sphericity is
dimensionless.

Measured on two concentric spheres standing in for production and reference
(`spacing_invariance_check.py` in this directory):

| spacing error | relative volume discrepancy | relative area discrepancy |
|---|---|---|
| none (baseline) | −0.146929040 | −0.100203671 |
| common +3% | −0.146929040 | −0.100203581 |
| common ×2 | −0.146929040 | −0.100203581 |
| axial-only +3% | −0.146929040 | −0.100205556 |
| axial-only +20% | −0.146929040 | −0.100273283 |

The volume discrepancy is bit-identical across all five. The area discrepancy
moves by 0.007 percentage points under a **20%** axial error — the magnitude the
report itself attributes to uncorrected refractive-index mismatch, and the
residual is still two orders of magnitude below any threshold S1 would set. So
the scale error does not contaminate S1's bias metric at all; what it perturbs
is the **stratification**, since `rho_in` and anisotropy are functions of the
z:xy ratio. The true dependency is therefore narrow — S1 needs the per-axis
*ratio* to be trustworthy enough to assign strata, not absolute traceable scale.
**Consequence:** S1, the long-lead annotator-limited study, need not queue behind
two instrument studies; absolute-unit calibration can run in parallel.
`ROADMAP.md` §3 item 1 states the same overstated dependency (§6).

**M6 — S1's downstream-bias threshold is an inherited placeholder.** The
"1% volume / 5% area" pair is real — it is the claim register's criterion for M1
and M4 — but volume's domain constants were derived for the *area* estimator and
volume has never declared its own. That is roadmap R-3, which this report puts
out of scope. So S1's most important acceptance rule rests on a number the
repository itself flags as not yet derived for the quantity it governs. Either
bring R-3 in scope ahead of S1, or predeclare S1's threshold on its own terms
and record that it is not the register's criterion.

**M7 — C2 does not tie the claim to the measured residual.** The report
correctly refuses to force a 1% axial criterion and asks for a stated systematic
term instead. But the T-tiers exist to bind claims to outcomes, and the report
never says what claim survives an axial residual of, say, 5%: at some magnitude
the physical-unit export stops being reportable at all. Name the ceiling as a
function of the measured residual.

## 4. What the report contributes

Genuine additions, independent of the defects above:

- **The axial axis is named as the dominant term, with a magnitude.** An axial
  scaling factor near 1.2 under refractive-index mismatch, if it holds for this
  configuration, is orders of magnitude above the derived 0.25% tolerance and
  dwarfs every estimator error in this repository. It also supplies the physical
  reason our own derivation could only conclude "not plausibly achievable
  axially".
- **A concrete measurement method per tier**, which §2 lacked: beads in
  refractive-index-matched gel at specimen-like depth for stability; traceable
  lateral grid plus depth-dependent axial scaling for accuracy.
- **Consensus fusion over majority vote** for multi-annotator ground truth, with
  the right reason (majority voting can flatter the algorithm), and per-rater
  agreement reported as the inter-rater result rather than assumed away.
- **Size and depth stratification of the viability separation metrics** on
  dye-penetration grounds — a mechanism-based stratification §5.5 did not
  require, and well aimed: penetration failure is exactly the artefact the
  `both_markers_low` state would otherwise absorb silently.
- **Naming the few-cluster corrections** (Satterthwaite / between-within, CR2)
  and insisting the omnibus anti-conservatism carry a measured number rather
  than a caveat.
- **Correct refusal to let Q3 count as biological validation**, and correct
  identification of the two zero-cost studies (Q3 refusal contract, T1 null
  check) as the things to run first.

## 5. Required changes before this can be filed

1. Attach the bibliography with DOIs; mark full-text versus abstract-only (B1).
2. Resolve the V1 reference-standard circularity by one of the three named exits
   (B2).
3. Restate S1's acceptances per case or as an upper confidence limit, or scope
   the S1 claim explicitly to central tendency (B3).
4. Add a precision calculation per study and derive sample sizes from it (B3).
5. Re-sequence Gap 4 — upgrade then qualify, or qualify then upgrade, but not T1
   against a path already scheduled for replacement; state the CR2
   implementation and verification burden (B4).
6. Allocate the calibration budget across T1 and T2 instead of applying 0.25%
   twice; require the upper confidence limit for T1 (M1).
7. Replace `ρ ≥ 0.8` for the per-object claim with confusion-matrix and ordinal
   agreement criteria; keep correlation for V2 only (M2).
8. Make the interval construction an outcome of the coverage study rather than a
   predecided switch to BCa; price the n ≥ 5 floor (M3).
9. Drop or relocate the nuclear-level precedent numbers (M4).
10. Re-sequence S1 ahead of absolute-scale calibration, carrying the per-axis
    ratio requirement into C1 (M5).
11. Decide whether R-3 returns to scope; if not, predeclare S1's own threshold
    (M6).

## 6. Errata candidates against this repository

Two defects above are ours, and the report inherited them faithfully. Both are
recorded as candidates only: neither is corrected by this review, and both
change owner-facing content.

- **`INTENDED_USE_AND_ESTIMANDS.md` §5.5, Q1 reference (a).** "A cell-resolved
  live/dead stain with nuclear segmentation on the same objects" does not say
  *whose* nuclear segmentation, and the only one available is the pipeline's own
  — unvalidated (R-13, NOT ASSESSED) and exploratory-only by our own policy. The
  preference order should name an independent counting path, or state that (a)
  is unavailable until R-13 is addressed and (b) is therefore primary, which
  narrows the achievable claim to V2.
- **`ROADMAP.md` §3 item 1.** "R-1 (calibration) before R-6 (segmentation) — a
  segmentation study on uncalibrated acquisitions cannot separate a scale error
  from a mask bias: both appear as a systematic volume discrepancy." The
  discrepancy S1 reports is a ratio of two measurements on one grid, so the
  spacing product cancels exactly for volume and to 0.007 pp for area even at a
  20% axial error (M5, measured). The dependency is real only for the per-axis
  *ratio*, and only through stratum assignment. As written the item overstates
  the dependency and defers the longest-lead study for no gain. The correction
  would reorder the critical path, so it is the owner's call, not a
  documentation fix.

## 7. Status after this review

Unchanged. Gate NOT PASSED, 0/6. No study executed, no threshold adopted, no
`[OWNER]` placeholder filled in. This review is a document review; it is not
evidence about the pipeline and creates no claim.
