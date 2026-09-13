# Intended use, estimands and claim status — frozen 2026-09-13

Status vocabulary, used with exactly these meanings throughout this document and
the round-4 roadmap:

| Status | Meaning |
|---|---|
| **SUPPORTED** | A defined estimand with a stated operating domain, an independent reference standard, and a validation record on data not used for tuning, meeting a criterion predeclared before the result was seen. |
| **SUPPORTED WITH LIMITATIONS** | As above, but the validation record covers a narrower range than the claim's declared domain, or rests on at least one condition the software cannot check. The limitation is named. |
| **NOT QUALIFIED** | A validation record exists and the predeclared criterion was not met, or the method is disqualified by a demonstrated defect. |
| **INSUFFICIENT EVIDENCE** | The estimand and criterion are defined but the evidence needed to decide is absent or too weak to decide. |
| **NOT ASSESSED** | No study has been performed and none is in progress. Distinct from INSUFFICIENT EVIDENCE: nothing has been attempted. |

A status is a statement about *evidence*, never about intent, code quality, or
test-suite results. Passing unit tests is not evidence for any row below.

This document freezes the claims. It changes no algorithm, threshold, criterion
or domain constant, and it is not a validation record. Where a claim is narrower
than what the software will compute, the software is not being restricted here —
the *claim* is. The pipeline will still emit a number outside its qualified
domain; the flags and policy fields that accompany it are the mechanism by which
that number is marked as unqualified (`surface_in_qualified_domain`,
`surface_qualification_scope`, `measurement_basis`, `morphology_flags`,
`viability_reason`).

---

## 1. The three intended-use areas

### IU-1 — Quantitative organoid morphology (brightfield or membrane fluorescence)

**Claim.** For organoid-level objects segmented from a 3D stack, the pipeline
reports physical size and shape descriptors in µm-based units, suitable for
comparing conditions within one experiment.

**Experimental unit.** The organoid object is the *measurement* unit; the well
(imaging unit) and the biological replicate are the *inference* units
(`SCIENTIFIC_SPEC.md` §5–6). No claim below is a claim about a single organoid.

**Not permitted.** Clinical diagnosis, patient stratification, any go/no-go
decision without independent validation, absolute cross-instrument comparison of
physical units (blocked by SG-6, below), and comparison of shape descriptors
between conditions whose objects fall in different qualification domains without
reporting the stratification.

### IU-2 — Exploratory nuclear morphology

**Claim.** For nuclear objects, the pipeline reports descriptors intended for
hypothesis generation and internal comparison only. No nuclear measurement is
offered as a validated quantity. See §4 for the policy that makes "exploratory"
operational rather than a disclaimer.

### IU-3 — Quantitative live/dead assessment

**Claim.** For organoid-level objects with calcein and PI channels and
batch-matched live and dead controls, the pipeline assigns each object one of
four signal states and reports state fractions per well. The claim is about
*marker-signal states of segmented objects*, not about the proportion of living
cells. See §5.

**Terminology, binding.** The term **live-cell fraction** is not used anywhere in
this repository's outputs or documents, because no output has a cell-count
denominator: the classification unit is a segmented organoid object, not a cell.
A marker-intensity summary must not be renamed into a cell-count quantity, and a
volume-weighted or intensity-weighted proxy must not be reported as a live-cell
fraction without independent biological validation against a cell-resolved
reference. The implemented state names (`viable_like`, `compromised_like`,
`mixed_signal`, `indeterminate`) and the aggregate names
(`fraction_viable_like`, …) already satisfy this and must not be "cleaned up"
into biological language.

---

## 2. Claim register

Every row is an estimand, not a feature. "Domain" splits into what the software
can check from a mask (**machine-checkable**) and what it cannot
(**scope assumption**) — the distinction matters because only the first can be
enforced at run time, and a claim resting on the second is at best SUPPORTED
WITH LIMITATIONS.

| # | Estimand | Domain (machine-checkable) | Domain (scope assumption) | Reference standard | Record | Status |
|---|---|---|---|---|---|---|
| M1 | Total enclosed volume `V_env` of an organoid object, µm³ (§3) | `rho_in ≥ 10`, anisotropy ≤ 4, QC-eligible (not border-truncated, ≥ min Z slices, within volume limits) | solid or enclosed-lumen object; correct voxel spacing (SG-6); correct segmentation (SG-3) | analytical rasterised solids, exact volume | `2026-09-13-analytical-geometry-record` (grid, FAIL unrestricted); `vv_confirm2.csv` (confirmation, in-scope smooth cases ≤ 0.76%) | **SUPPORTED WITH LIMITATIONS** — see §3.5 for exactly what is and is not covered |
| M2 | Segmented material volume `V_seg`, µm³ (§3) | same as M1 | segmentation resolves the inner lumen boundary — **unevidenced** | none | none | **NOT ASSESSED** |
| M3 | Enclosed void fraction `f_void` (§3) | same as M1 | as M2 | none | none | **NOT ASSESSED** |
| M4 | Surface area of an organoid object, µm² | `rho_in ≥ 10`, anisotropy ≤ 4 (`surface_in_qualified_domain`) | smooth closed surface, no crease at the voxel scale — **not machine-checkable** | analytical rasterised solids, exact area | `2026-09-12-surface-method-development` (development), `vv_confirm2.csv` (confirmation, in-scope smooth cases ≤ 0.95%) | **SUPPORTED WITH LIMITATIONS**; unrestricted SG-1 is **NOT QUALIFIED** (25.81% worst case on the grid) |
| M5 | Sphericity, surface-to-volume ratio | `sphericity_in_qualified_domain`, `surface_to_volume_in_qualified_domain` | as M4 | as M4 | as M4 | **SUPPORTED WITH LIMITATIONS** (inherits M1 and M4) |
| M6 | Equivalent diameter, principal axes, elongation, axis ratios | `rho_in ≥ 10`, QC-eligible | as M1 | analytical solids (exact axes) | grid record | **INSUFFICIENT EVIDENCE** — axis-length accuracy has not been extracted from the record as a criterion-bearing item |
| M7 | Physical scale of every µm-valued output | none — spacing is read from metadata and never verified against a standard | metadata spacing equals true voxel size | traceable stage micrometer / axial standard | none | **NOT ASSESSED** (SG-6). Binds *every* row above: a scale error is systematic across all objects and is invisible to all checks in this repository |
| M8 | Segmentation correspondence: does the mask boundary coincide with the biological object boundary | none | modality, marker, model checkpoint, preprocessing all fixed to the validated configuration | expert manual annotation on real specimens | none | **INSUFFICIENT EVIDENCE** (SG-3) |
| N1 | Nuclear volume, axes, intensity summaries, CV of chromatin | `rho_in ≥ 10` etc. for the conditional metrics | nuclear segmentation validity — unevidenced | none | none | **NOT ASSESSED**, exploratory only (§4) |
| V1 | Per-object marker-signal state (§5) | calibration available for the batch; object QC-eligible for intensity | markers report the intended biology; controls are representative | independent viability assay (§5.4) | none | **INSUFFICIENT EVIDENCE** |
| V2 | Well-level state fractions (§5.3) | as V1, ≥ 1 classified object | as V1 | as V1 | none | **INSUFFICIENT EVIDENCE** |
| S1 | Condition-level contrasts (LMM/OLS-clustered, BH-FDR) and the condition-level percentile-bootstrap interval | replicate structure declared in the design table | the replicate hierarchy in the design table matches the experiment | simulation under known truth; independent reanalysis | none | **NOT ASSESSED** (SG-5) |

Three properties of this table are load-bearing and easy to lose:

1. **A declared domain is not a covered range.** M1 and M4 declare
   `rho_in ≥ 10` with *no upper bound*; the confirmation set covers
   `rho_in` 10.43–15.05 and anisotropy {1.2, 2.2, 3.5, 4.0}. Above `rho_in` 15.05
   the claim rests on an extrapolation argument recorded in
   `surface_crofton.py:135`, not on measurement. Any report that quotes the
   qualified domain must also be able to state the covered range.
2. **The scope-assumption column is where the surface claim actually lives.**
   Inside the machine-checkable gate, the creased analytical classes (box,
   cylinder; n = 49) reach 9.10% area and 3.20% volume error against 0.95% and
   0.76% for the smooth in-scope classes. The criterion is met because the
   *scope assumption* excludes them, and no mask-computable crease indicator
   exists (`crease_indicator_rejected.md`). A user cannot be told the software
   checked this.
3. **The evidence covers a discrete set of voxel spacings.** The frozen Crofton
   weight tables exist for 11 z:xy ratios {1, 1.2, 1.5, 2, 2.2, 2.857, 3, 3.125,
   3.5, 4, 5} with equal lateral spacings. Every spacing in both evidence sets is
   in that set. A realistic acquisition at 1 × 0.325 × 0.325 µm is not, and is
   solved at run time; the minimax program has a non-unique optimum
   (`crofton_weights/PROVENANCE.md`), so a run-time solve is a *different* weight
   vector with the same worst-case plane-response deviation. Such measurements
   carry `surface_weights_origin ∈ {cache, solved}` and are **outside the
   evidence-bearing configuration**, whatever the domain flag says. This is a
   coverage limitation of the evidence, not a code defect; it is roadmap item
   R-7.

---

## 3. The organoid volume estimand (Task 2)

The software already computes both candidate quantities and exports them side by
side; what has been missing is a decision about which one the *claim* is about,
and what the other one means. Both are voxel-count estimators — the count is
multiplied by the voxel volume `sz·sy·sx`.

| Symbol | Column | Definition | Operator |
|---|---|---|---|
| `V_env` | `volume_um3` when `measurement_basis = filled_envelope` | volume enclosed by the object's outer boundary, **including enclosed lumens** | `features.outer_envelope` = `binary_fill_holes` on a 1-voxel pad (`features.py:108`) |
| `V_seg` | `segmented_volume_um3` | volume of voxels the segmenter assigned to the object | raw label |
| `f_void` | `enclosed_void_fraction` | `(envelope_voxels − segmented_voxels) / envelope_voxels` (`features.py:190`) | both |

### 3.1 Decision

**The primary organoid size estimand is `V_env`, the total enclosed volume.**
`V_seg` and `f_void` are reported alongside it and are required accompaniments,
not optional extras. The organoid workflow already measures on this basis
(`organoid_measurement_workflow.py:146`,
`measurement_policy("organoid", "filled_envelope")`); the multilevel cell and
nucleus path already measures raw labels
(`multilevel_measurement_workflow.py:238`, `morphology.py:41-42`), which is
correct for those objects and must not be unified with this decision.

**Reasons, in the order of their weight:**

1. `V_env` is the quantity the existing analytical evidence is about. The
   phantoms are exact rasterisations of *solid* analytical bodies, so their
   voxel count is an envelope volume. The accuracy statement in M1 transfers to
   `V_env` directly, and transfers to `V_seg` only for objects with no enclosed
   void. There is **no analytical evidence for the shell geometry** that `V_seg`
   measures on a lumen-containing object.
2. `V_seg` on a hollow object depends on the *inner* boundary, which depends on
   marker penetration, lumen fluid signal, and how the segmenter behaves on an
   internal surface — three unevidenced factors, none covered by SG-3 as
   currently scoped.
3. `V_env` is the quantity comparable to diameter- and area-based growth
   readouts used in brightfield organoid assays, which is the dominant intended
   comparison in IU-1.

### 3.2 What `V_env` does **not** mean

`V_env` is not tissue mass and must not be described as growth of tissue. A
cystic organoid that swells by fluid accumulation increases `V_env` with no
increase in `V_seg`. Therefore:

> **Binding reporting rule.** A comparison of `V_env` between conditions is
> interpretable as a size comparison only. If the conditions differ in `f_void`,
> a `V_env` difference must not be attributed to tissue growth. Any report
> comparing `V_env` across conditions must also report `f_void` per condition,
> and must state whether `f_void` differs.

This rule is enforceable today from the exported columns and is a documentation
and reporting obligation, not a code change.

### 3.3 The four hard cases, and what the operator actually does

`binary_fill_holes` is a **topological** operator: it fills exactly those
background components that are not connected to the image border (after the
1-voxel pad). That has consequences the estimand definition must confront
rather than assume away.

| Case | What `V_env` becomes | Assessment |
|---|---|---|
| **Solid organoid** | segmented volume; `f_void ≈ 0` | covered by M1 within its domain |
| **Enclosed lumen** (single or multiple, fully interior) | lumen included; `f_void > 0` | the intended behaviour, and the reason for the decision in §3.1 |
| **Cystic / lumen open to the exterior** (a cup, or a lumen breached by the crop or by a segmentation gap) | the void is connected to the border, so it is **not** filled, and `V_env` collapses to `V_seg` | **defect risk, unresolved.** The measured quantity changes discontinuously with topology. The same organoid can move between the two behaviours between time points, or through a Z-truncation, producing an apparent volume change of the lumen's whole size with no biological change. |
| **Necrotic core** | dead material is usually still segmented, so it is neither a void nor distinguishable; `f_void ≈ 0` | `V_env` includes the necrotic core. Distinguishing it requires marker evidence, not geometry. **NOT ASSESSED.** |
| **Fragmented object** (one biological organoid split into several labels, or several organoids merged into one) | one label is one object, so fragmentation/merging changes the object set, not the estimand | handled by QC (`morphology_flags`: `encloses_other_instance`, `excess_foreground_review_segmentation`) and by SG-3, not by the volume definition |

### 3.4 Required additions before M1 can move above SUPPORTED WITH LIMITATIONS

These are the open items the estimand decision creates. They are design
requirements, recorded here and carried into the roadmap; none is implemented by
this document.

1. **A topology flag.** `f_void` cannot distinguish "solid" from "cavity open to
   the exterior" — both give `f_void = 0`. The number of filled background
   components per object, and whether any background component adjacent to the
   object touches the image border, must be exported so that the open-cavity case
   is identifiable in analysis rather than silently measured as solid.
   (Roadmap R-3.)
2. **A declared behaviour for the open-cavity case.** Options are: exclude such
   objects from volume claims; report them with an explicit flag and no
   qualification; or define a geodesic/morphological closing envelope. The third
   changes the estimand and would need its own validation; the first two do not.
   The decision is the owner's and must be made before the volume qualification
   study fixes its case list, because the study's phantom families follow from it.
3. **Hollow and open-cavity phantoms in the volume study.** The present evidence
   contains no hollow object of any kind. Until it does, M1's domain carries
   "solid or enclosed-lumen" as a *scope assumption*, exactly parallel to the
   smoothness assumption on M4 — and equally unverifiable from a mask.
4. **Spacing calibration (M7/SG-6).** `V_env` inherits the full three-axis scale
   error: a 3% error in all three spacings is +9.27% volume, and a lateral-only
   3% error is +6.09%. Derivation in
   `evidence/2026-09-13-volume-audit/STUDY_PROPOSALS_SG3_SG6.md` §2.1.

### 3.5 Exactly what M1's evidence covers

So that the status is not read as stronger than it is:

- Covered: exact rasterised solids of the families sphere, ellipsoid, cylinder,
  capsule, torus, box; `rho_in` 10.43–15.05 on the confirmation set; the 11
  packaged voxel-spacing ratios; smooth closed in-scope classes at ≤ 0.76%
  absolute relative volume error (n = 96).
- Not covered: any hollow object; any real segmentation; `rho_in` above 15.05
  and below 10; anisotropy above 4; the open-cavity topology case; necrotic
  cores; physical scale.
- Unrestricted status: **NOT QUALIFIED.** Over the full analytical grid the
  volume criterion (< 1% for every case) fails, with 18.26% worst case at
  `rho_in < 3` and 68.2% of that stratum above 1%. The qualification is
  domain-restricted, and the domain is not machine-checkable in full.

---

## 4. Nuclear exploratory policy (Task 7)

"Exploratory" is only meaningful if something follows from it. This is what
follows.

**Permitted.** Internal comparison within one experiment and one acquisition
configuration; hypothesis generation; ranking; visual review; quality control of
nuclear segmentation.

**Not permitted.** Any absolute claim about nuclear size, shape, or chromatin
organisation; any cross-experiment or cross-instrument comparison; use of a
nuclear metric as a primary or secondary endpoint; inclusion of a nuclear metric
in a statistical test that is reported as a result rather than as an
exploration; inference from nucleus-to-cell ratios, which compound two
unvalidated segmentations.

**Required labelling.** Every exported nuclear metric must carry the
`measurement_policy` object for its level (already implemented:
`measurement_policy("nucleus", "raw_label")`), and every figure or table
containing a nuclear metric must be labelled exploratory in the artifact itself,
not only in accompanying prose.

**Why nuclear metrics cannot be promoted by the existing surface evidence.** The
analytical qualification is about a rasterised analytical body. A nucleus at
typical confocal spacing is close to the resolution limit, the `rho_in ≥ 10`
condition will frequently fail, and the smoothness scope assumption is untested
for chromatin-textured boundaries. The conditional surface metrics already
inherit `surface_in_qualified_domain` at the nuclear level, which is necessary
but not sufficient: passing that flag does not make a nuclear measurement
validated, because nuclear *segmentation* is unevidenced (M8, N1).

**Promotion path.** A nuclear metric moves out of exploratory status only via:
(i) a nuclear segmentation validation record against expert annotation on the
target modality; (ii) evidence that the nuclear objects fall inside the declared
resolution domain in the intended acquisition regime; (iii) a predeclared
accuracy criterion on a reference standard for the specific metric. Nothing less
— in particular, not agreement between two of the pipeline's own metrics.

---

## 5. The viability estimand (Task 8)

### 5.1 What the implementation actually computes

Read from `phenotyping/viability.py` and `statistics/aggregation.py`:

- One row per **segmented organoid object**. There is no cell-level or
  voxel-level viability path in the repository.
- Per object and per marker, the background-corrected mean intensity inside the
  object mask is scaled against that batch's control medians:
  `c = (calcein_mean_bg_corrected − dead_median) / (live_median − dead_median)`,
  and the mirrored expression for PI, each clipped for interpretation but
  flagged (`outside_control_range`) when the raw scaled value leaves [0, 1].
- Control medians are formed through the replicate hierarchy (object → unit →
  biological replicate → batch), with `calcein_low` the dead-control median and
  `calcein_high` the live-control median, mirrored for PI
  (`viability.py:75-78`). Separation is `(high − low) / noise` where the noise
  floor is the larger of the between-replicate MAD and the median background
  noise MAD (`viability.py:80-82`) — so a large live/dead difference does not
  count as separation if the controls themselves are that variable.
- When calibration succeeds, the recorded reason is
  `control_scaled_rules_require_external_validation` (`viability.py:88`). The
  implementation therefore already declares this readout unvalidated; §5.4 is
  the explicit form of that statement, not a new restriction.
- Calibration is **refused** with a named reason when: the mode is uncalibrated; no QC-eligible control objects exist;
  controls span multiple conditions; there are fewer than
  `min_control_replicates = 2` independent control replicates; or a marker's
  live/dead separation SNR is below `min_control_separation_snr = 3.0`.
- Classification is a rule on the two scaled values with gates
  `low_gate = 0.30`, `high_gate = 0.60`:
  `viable_like` (`c ≥ high` and `p ≤ low`), `compromised_like`
  (`p ≥ high` and `c ≤ low`), `indeterminate` (both markers low, or QC-ineligible,
  or uncalibrated, or nonfinite), `mixed_signal` (anything else).

### 5.2 The estimand, stated

> **V1.** For a segmented organoid object that is QC-eligible for intensity
> measurement, in a batch with an available two-marker calibration, the reported
> `viability_state` is the class assigned by the frozen rule above to the pair of
> batch-control-scaled, background-corrected **mean marker intensities over the
> object mask**.

This is a statement about marker signal, not about cells. It is
**object-average**: a shell of live cells around a dead core and a uniformly
intermediate object can produce the same mean and the same state. The estimand
inherits the geometry of the mask — in particular, for lumen-containing objects
the mean intensity is taken over the measurement mask, so the lumen contributes
to the average whenever the basis includes it. Which basis is used for intensity
must be stated in the record alongside the state (roadmap R-8); it is not a
free choice to be made per analysis.

The four-state output, including an explicit `indeterminate` class, is retained
as a requirement: no object may be forced into live or dead. `mixed_signal` and
`indeterminate` are results, not failures, and must be reported, not dropped.

### 5.3 The aggregate estimand

> **V2.** For an imaging unit (well), `fraction_<state>` is the number of
> objects in that state divided by the number of **morphology-QC-eligible
> objects** in that unit (`aggregation.py:19, 28-29`); unit values are then
> averaged over units within a biological replicate and over replicates within a
> condition.

The denominator is objects, not cells, not area, not volume. Three consequences
that must be stated wherever V2 appears. The fraction is **unweighted by object
size**, so one large and one small organoid count equally. Objects failing
morphology QC leave the denominator entirely, so `n_excluded` must be reported
with the fraction (it is exported). And an object that passes morphology QC but
fails *intensity* QC stays in the denominator as `indeterminate` — so
`fraction_viable_like` is diluted by measurement failures rather than being
computed over classifiable objects only, which is why the state breakdown, not
`fraction_viable_like` alone, is the reportable result. A size-weighted variant is a *different* estimand and would
need its own definition and validation — it is not a presentation choice.

### 5.4 Why V1/V2 are INSUFFICIENT EVIDENCE, not SUPPORTED

There is no validation record. Specifically missing:

1. **No reference standard has been applied.** Neither an independent viability
   assay on the same objects nor a cell-resolved ground truth exists in this
   repository.
2. **The gates are unevidenced.** `0.30` and `0.60` are configuration defaults
   with no record deriving them from control distributions or from a target
   error rate. They must not be adjusted to improve any observed result; that
   would be tuning on the outcome.
3. **Marker validity is assumed.** That calcein retention and PI exclusion
   report the intended biology *in this specimen and mounting medium* is a
   biological claim requiring biological evidence.
4. **Object-mean intensity as the feature is a modelling choice**, not a
   measured quantity, and its adequacy for heterogeneous objects is untested.

### 5.5 The validation study (Task 9)

Design, for the owner to schedule. Acceptance numbers marked `[OWNER]` are
placeholders: they must be fixed **before** any data are seen, and they cannot
be chosen from this repository's existing observations.

**Q1. Do the four states separate on an independent reference?**
Reference standard: an independent viability readout on the same objects, in
declared order of preference — (a) a cell-resolved live/dead stain with nuclear
segmentation on the same objects, giving a per-object fraction of dead nuclei;
(b) a destructive bulk assay (e.g. ATP) at the *well* level, which validates V2
only, not V1; (c) a time-course with a known lethal insult as a positive control
and vehicle as negative, which validates *direction* only. Design: ≥ `[OWNER]`
biological replicates × ≥ `[OWNER]` wells, spanning a graded insult series so
that intermediate states are populated rather than only the extremes. The
biological replicate is the unit of inference. Criterion, predeclared:
monotone association between the per-object dead-nucleus fraction and the state
ordering `viable_like < mixed_signal < compromised_like`, with a `[OWNER]`
minimum separation between the `viable_like` and `compromised_like` groups; and
for the well-level readout, a `[OWNER]` minimum rank correlation with
`fraction_viable_like` across the insult series.

**Q2. Are the gates defensible, and how sensitive is the result to them?**
Report the full two-dimensional `(c, p)` distribution with the gate lines drawn,
per condition; recompute all state fractions at gates
`(0.25, 0.55)`, `(0.30, 0.60)`, `(0.35, 0.65)` and report the change in the
primary contrast. This is a **sensitivity characterisation, not a
qualification**: it may not be used to select gates, and the frozen gates stay
frozen. If the primary contrast reverses within that range, the conclusion is
that the readout is gate-sensitive and V2 cannot support the contrast.

**Q3. Does the calibration refusal machinery work as intended?**
Adversarial cases built to trigger each named refusal
(`insufficient_independent_control_replicates`,
`<marker>_controls_not_separated`, `controls_span_multiple_conditions`,
`no_QC_eligible_controls`, `viability_mode_uncalibrated`), verifying that no
state is emitted and the reason is recorded. This part is mechanical and can be
done today with synthetic tables — it is a software contract check, and it must
be reported as such and never as biological validation.

**Q4. What is the indeterminate rate, and is it informative?**
Report the rate of `indeterminate` and `mixed_signal` per condition with the
reason breakdown. A high `both_markers_low` rate in a condition is itself a
finding (loss of signal, over-fixation, penetration failure) and must not be
treated as missing data. Criterion: none — this is characterisation. If the
indeterminate rate exceeds `[OWNER]` in any condition, V2 is not reportable for
that condition.

**Ordering.** Q3 first (cheap, mechanical, no specimens). Q1 and Q2 from the
same experiment. Q1 cannot be interpreted before segmentation validation (M8) at
least for the modality used, since the object set itself is the denominator of
V2 — this dependency must be stated in the study record rather than discovered
afterwards.

---

## 6. What this document does not do

- It does not change any algorithm, threshold, criterion, domain constant or
  gate status.
- It does not create a validation record. Every row in §2 that says
  NOT ASSESSED or INSUFFICIENT EVIDENCE still says so afterwards.
- It does not redefine intended use around observed failures. IU-1 to IU-3 are
  the areas named in the round-4 brief; the narrowing in §2 is a narrowing of
  *claims to what is evidenced*, and every narrowing is traceable to a missing
  study listed in the roadmap.
- It does not convert the analytical geometry evidence into evidence about
  segmented objects, nuclei, or viability.

## 7. Revision

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-09-13 | First issue. Claim register frozen; organoid volume estimand resolved to `V_env` with `V_seg`/`f_void` as required accompaniments; nuclear exploratory policy and viability estimand defined; five round-3 defects corrected under `evidence/2026-09-13-volume-audit/ERRATA.md`. |
