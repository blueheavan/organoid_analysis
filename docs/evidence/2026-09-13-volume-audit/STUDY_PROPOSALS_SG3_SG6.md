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

**1.4 Scope of a possible SG-3 PASS.** Whatever the study returns, it can
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
estimator work: a 3% lateral scale error is a ~3% area error and a ~9% volume
error, applied to *every* object and invisible to every check in this
repository. That is larger than the qualified estimator's entire in-domain error
budget.

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

## 3. Sequencing note

SG-6 is logically upstream of SG-3: a segmentation study run on uncalibrated
acquisitions cannot separate a scale error from a mask bias, because both appear
as a systematic volume discrepancy against the reference standard. If only one
study can be resourced, calibration first is the cheaper and more informative
order. This is a recommendation, not a result.
