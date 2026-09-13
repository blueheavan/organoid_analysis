# Validation roadmap — 2026-09-13 (round 4)

Companion to `docs/INTENDED_USE_AND_ESTIMANDS.md`, which freezes the claims this
roadmap is about. Status words are the five fixed terms defined there:
**SUPPORTED**, **SUPPORTED WITH LIMITATIONS**, **NOT QUALIFIED**,
**INSUFFICIENT EVIDENCE**, **NOT ASSESSED**.

Nothing in this round produced a new validation record. The overall verdict is
unchanged: **the pipeline is NOT QUALIFIED for release as a validated
measurement system**, and the parts that are defensible are defensible only
within named domains and with named limitations.

## 1. Final intended-use statement

The software is a **research measurement and exploration tool for 3D organoid
image analysis**, for use by the intended users in `SCIENTIFIC_SPEC.md` §2, in
three areas:

- **IU-1, quantitative organoid morphology.** Size and shape of segmented
  organoid-level objects, for comparison of conditions **within one experiment
  and one acquisition configuration**. Physical-unit outputs carry an
  unquantified systematic scale uncertainty until SG-6 is done; relative
  comparisons within one configuration are the supportable use.
- **IU-2, exploratory nuclear morphology.** Hypothesis generation only, under
  the policy in `INTENDED_USE_AND_ESTIMANDS.md` §4. No nuclear metric is a
  validated quantity and none may serve as an endpoint.
- **IU-3, marker-signal state assessment ("live/dead").** A per-object
  four-state classification of calcein/PI signal against batch-matched
  controls, with **INSUFFICIENT EVIDENCE** for any biological interpretation.
  Not a live-cell fraction; see the binding terminology rule.

Not permitted, unchanged: clinical diagnosis or stratification, treatment
decisions, regulatory submission, absolute cross-instrument comparison of
physical units, and any claim that rests on a status of NOT ASSESSED or
INSUFFICIENT EVIDENCE.

## 2. Prioritized roadmap

Priority bands: **P1** — cheap, needs no new specimens, and either unblocks
other work or prevents a wrong claim being made today. **P2** — required for the
core intended use (IU-1). **P3** — required to widen claims beyond IU-1.
**P4** — outside the current intended use.

| ID | P | Scientific question | Current evidence | Missing evidence | Required study / action | Proposed acceptance logic | Claim enabled if successful |
|---|---|---|---|---|---|---|---|
| **R-1** | P1 | Is the reported voxel spacing the true voxel size, and is it stable? | none. Per-axis provenance is now recorded, which is not verification | traceable lateral and axial scale, in the specimen's medium and depth range; between-session stability | `STUDY_PROPOSALS_SG3_SG6.md` §2, tiered per §2.3.1 | T1: per-axis stability below the derived tolerance (`< 0.25%` if the estimator budget is to dominate). T2: per-axis deviation from a traceable standard within a declared tolerance, upper confidence limit included. Per stratum; one stratum failing fails that stratum only | T1 → physical-unit comparison of conditions within one configuration (the core of IU-1). T2 → reporting absolute units with a stated systematic term |
| **R-2** | P1 | Can the open-cavity topology case be identified at all? | `f_void` is 0 for both a solid object and one whose cavity is open to the exterior — the two are indistinguishable today | a per-object topology diagnostic | export `filled_void_count` and `open_cavity_suspected` from the operator that already runs (`features.py:108`); unit tests on synthetic shells. Changes no measurement | contract test: a closed shell reports ≥ 1 filled component; the same shell with an aperture reports 0 and the open-cavity flag | prerequisite for R-3 and for SG-3's topology stratification; removes a silent failure mode from IU-1 today |
| **R-3** | P1 | What is volume's own operating domain, and what is the criterion? | the domain constants were derived for the **area** estimator; volume has never declared one | a declaration in volume's own terms, frozen before any confirmation measurement | `VOLUME_QUALIFICATION_PROTOCOL.md` §4 steps 1–2 and §6.4; an owner decision on the open-cavity behaviour | none — this is a declaration, not a study. It is a precondition for any volume verdict | makes SG-2a constructible; without it a volume PASS would rest on borrowed constants |
| **R-4** | P1 | Does the refusal machinery in the viability path behave as specified? | the reasons exist in code and are read by tests of the surrounding contract | adversarial cases for each named refusal reason | `INTENDED_USE_AND_ESTIMANDS.md` §5.5 Q3, synthetic tables only | every constructed case yields no state and the specified reason string | a software-contract claim only. **Must never be reported as biological validation** |
| **R-5** | P1 | Is the object row treated as the inference unit? | none | a null simulation with strong within-replicate correlation | `STUDY_PROPOSALS_SG3_SG6.md` §3.2 Q2 | type-I error at the nominal level under the null with no true effect | removes the most consequential possible defect in the inferential layer; cheap and falsifiable |
| **R-6** | P2 | Does the mask boundary correspond to the biological object boundary? | none. INSUFFICIENT EVIDENCE | expert annotation on the target modality, per the existing protocol plus the round-4 addenda §1.1–1.6 | `2026-09-11-measurement-vv/SEGMENTATION_VALIDATION_PROTOCOL.md` + addendum | predeclared in that protocol; downstream bias reported stratified by domain flag and by topology agreement, with the out-of-domain fraction as a primary result | the precondition for **every** measurement claim about real specimens. Without it, M1/M4 are claims about rasterised analytical bodies only |
| **R-7** | P2 | Is `V_env` accurate for hollow and open-cavity objects, and how large is the topology discontinuity? | solid bodies only; in-scope smooth cases ≤ 0.76% on the confirmation set | hollow, open-cavity, multi-void and truncated phantoms; the aperture sweep | `VOLUME_QUALIFICATION_PROTOCOL.md` §6 | `\|rel. err\| < 1%` for **every** in-domain closed/solid confirmation case; the aperture series is characterisation with no verdict | SG-2a on a declared volume domain; and the artefactual-volume-jump number a user needs to interpret a `V_env` change |
| **R-8** | P2 | Which measurement basis do viability intensities use, and is it recorded? | the organoid path measures on `filled_envelope`; the intensity basis is not stated in the viability record | an explicit, recorded basis per level | state it in the export contract and the record; no algorithm change | the record names the basis for every intensity summary | makes V1 reproducible and removes a silent analysis-time choice |
| **R-9** | P2 | Does the inferential layer's uncertainty statement have nominal coverage under this project's replicate structure? | none. NOT ASSESSED | simulation under known truth across replicate counts, including the OLS-fallback branch, G = 2–3, and the percentile-bootstrap interval at its own `n = 3` floor | `STUDY_PROPOSALS_SG3_SG6.md` §3.2 Q1, Q3 | coverage of the nominal 95% interval within a declared tolerance, reported per cell, per branch and separately for the bootstrap interval; a *domain of validity in replicate count* is an acceptable outcome | S1 moves from NOT ASSESSED to a scoped claim; the omnibus test's "screening only" label gets a number; the bootstrap interval either gets a minimum replicate count or a small-sample replacement |
| **R-10** | P3 | Do the four marker-signal states correspond to biological viability? | none. INSUFFICIENT EVIDENCE | an independent viability reference on the same objects, across a graded insult series | `INTENDED_USE_AND_ESTIMANDS.md` §5.5 Q1, Q2, Q4 | predeclared: monotone association with the reference and a minimum separation between `viable_like` and `compromised_like`; gate sensitivity reported but **not** used to select gates | IU-3 becomes a scoped biological claim. Depends on R-6 (the object set is V2's denominator) |
| **R-11** | P3 | Is the area estimator qualified at realistic acquisition spacings? | the evidence covers 11 packaged z:xy ratios; a realistic `1 × 0.325 × 0.325` is solved at run time against a non-unique optimum | confirmation evidence at non-packaged spacings, or packaged tables for the acquisition ratios actually used | extend the packaged ratio set to the project's real acquisition settings, or run a confirmation set with `surface_weights_origin = solved` recorded | the existing SG-1a criterion, evaluated on the added spacings | extends M4/M5 from "the packaged configuration" to the configuration the lab actually images in |
| **R-12** | P3 | Should the gate report qualification and characterisation separately? | `GATE_SEMANTICS_ANALYSIS.md` §3 quantifies four options; §7 specifies the recommended one | an owner signature in the decision record | implement Option C per §7, preserving the five invariants | the `a` items fail closed; the `b` items carry no verdict | the gate stops conflating "qualified for its claim" with "what it does outside it". Depends on R-3 and R-7 for SG-2a |
| **R-13** | P4 | Are nuclear metrics measurable? | none. NOT ASSESSED, exploratory only | nuclear segmentation validation, resolution-domain evidence in the real acquisition regime, and a metric-specific criterion | the three-part promotion path in `INTENDED_USE_AND_ESTIMANDS.md` §4 | as declared there | moves nuclear metrics out of exploratory status. Outside IU-1/IU-3; do not start before R-6 |

## 3. Inter-study dependencies

![Study dependency order]({{artifact:art_e6d4dfe5-f868-4a91-add3-e6335dbc81b3}})

(In-repository copy: `docs/figures/roadmap_dependency_graph.png`. Colour encodes
the priority band; arrows read "must precede".)

The order that matters:

1. **R-1 (calibration) before R-6 (segmentation).** A segmentation study on
   uncalibrated acquisitions cannot separate a scale error from a mask bias:
   both appear as a systematic volume discrepancy against the reference
   standard.
2. **R-2 → R-3 → R-7.** The topology diagnostic must exist before the volume
   domain can be declared coherently, and the domain must be declared before a
   confirmation set can be built to it.
3. **R-6 before any claim about real specimens**, including R-10, whose
   denominator is the object set.
4. **R-7 and R-1 are independent of each other**; the volume qualification is
   analytical and needs no imaging data, which is why it can proceed in
   parallel with the instrument work.
5. **R-12 last among the gate items**, since SG-2a cannot be constructed before
   R-3 and R-7.
6. **R-5, R-9 are independent of everything else** in method — but a correctly
   covered interval around a biased measurement is a precise statement of the
   wrong quantity, so they do not substitute for R-1 or R-6.

## 4. What is supportable today

| Claim | Status |
|---|---|
| Surface area of an organoid object, within `rho_in ≥ 10`, anisotropy ≤ 4, **and** the non-machine-checkable smoothness scope, at the packaged voxel-spacing ratios, as a relative quantity | SUPPORTED WITH LIMITATIONS (worst 0.951%, n = 96, confirmation set) |
| Total enclosed volume `V_env` under the same restrictions, for solid or closed-lumen objects | SUPPORTED WITH LIMITATIONS (worst 0.761%, n = 96) |
| Sphericity and surface-to-volume ratio under the same restrictions | SUPPORTED WITH LIMITATIONS (inherited) |
| The pipeline's refusal behaviour — that it declines to emit measurements and calibrations when its preconditions fail, with named reasons | SUPPORTED as a software contract (unit-tested); not a scientific claim |
| Reproducibility of a given analysis from its record | SUPPORTED (manifests, hashes, pinned seeds, `docs/REPRODUCIBILITY.md`) |

## 5. What is not supportable today

| Claim | Status |
|---|---|
| Any measurement claim about **real segmented objects** | INSUFFICIENT EVIDENCE (R-6) |
| Absolute physical units, or any cross-instrument comparison | NOT ASSESSED (R-1) |
| Area or volume for creased, hollow or open-cavity objects | NOT QUALIFIED / NOT ASSESSED (R-7) |
| Area at non-packaged voxel spacings | outside the evidence-bearing configuration (R-11) |
| Volume or area over the full analytical grid, unrestricted | NOT QUALIFIED (25.81% area, 18.26% volume worst case) |
| Any nuclear metric | NOT ASSESSED, exploratory only (R-13) |
| Any biological interpretation of the viability states; any live-cell fraction | INSUFFICIENT EVIDENCE (R-10) |
| Condition-level p-values and intervals as primary evidence | NOT ASSESSED (R-9) |
| Necrotic-core detection; lumen-versus-tissue growth attribution from `V_env` alone | NOT ASSESSED |

## 6. Next experiments, in the order recommended

1. **R-2, R-3, R-4, R-5** — all four need no specimens and no instrument time.
   They remove a silent failure mode, make a volume verdict constructible, and
   test two contracts. This is the cheapest block with the highest ratio of
   consequence to cost.
2. **R-1 (T1 tier)** — instrument time only, no specimens. Unlocks the core
   intended use.
3. **R-7** — compute only; the phantom extension is the bulk of the analytical
   work and can run in parallel with R-1.
4. **R-6** — the expensive one: annotator time on real specimens. Everything
   about real objects waits on it.
5. **R-9**, then **R-10**, then **R-11**, then **R-12**.
6. **R-13** only if nuclear metrics become part of the intended use.

## 7. Owner decisions outstanding

Each is recorded in the document named, and each blocks the item listed.

| Decision | Recorded in | Blocks |
|---|---|---|
| Gate architecture: Option A / B / C / D | `docs/GATE_SEMANTICS_ANALYSIS.md` §6 (unsigned) | R-12 |
| Behaviour for objects whose cavity is open to the exterior: exclude, flag, or define a closing envelope | `INTENDED_USE_AND_ESTIMANDS.md` §3.4 item 2 | R-3, R-7 |
| Volume's own domain and whether the 1% criterion is the biological requirement or inherited | `VOLUME_QUALIFICATION_PROTOCOL.md` §5 | R-3, R-7 |
| Calibration tolerances per tier, and whether T2 (absolute units) is in scope at all | `STUDY_PROPOSALS_SG3_SG6.md` §2.3.1 | R-1 |
| Re-qualification interval for instrument calibration | `STUDY_PROPOSALS_SG3_SG6.md` §2.4 | R-1 |
| The `[OWNER]` acceptance numbers in the viability study, and the indeterminate-rate ceiling | `INTENDED_USE_AND_ESTIMANDS.md` §5.5 | R-10 |
| Whether the analytical qualification proceeds now or only alongside SG-3 | `VOLUME_QUALIFICATION_PROTOCOL.md` §5 item 2 | scheduling of R-7 |
| Coverage tolerance and the replicate-count grid for the inference simulation | `STUDY_PROPOSALS_SG3_SG6.md` §3.2 Q1 | R-9 |

## 8. Known inconsistency: the gate's status words are not this vocabulary

The scientific gate prints per-item statuses using the same words with a looser
mapping than the one frozen this round. Today it reports:

| Gate item | Gate's word | This vocabulary, applied to the gate's own stated basis |
|---|---|---|
| SG-1 surface area | FAIL | NOT QUALIFIED on the unrestricted criterion; the restricted claim is SUPPORTED WITH LIMITATIONS (M4) |
| SG-2 voxel volume | FAIL | as above (M1) |
| SG-3 segmentation | INSUFFICIENT EVIDENCE | consistent — a protocol exists, no data |
| SG-4 viability assay | NOT ASSESSED | consistent for the *biological* claim; the marker-signal estimand V1 is INSUFFICIENT EVIDENCE, which is a stronger statement than the gate's |
| SG-5 statistics | INSUFFICIENT EVIDENCE | **NOT ASSESSED** — the gate's own basis is "generic tools only; no study design", i.e. nothing was attempted |
| SG-6 calibration | INSUFFICIENT EVIDENCE | **NOT ASSESSED** — "metadata is read and traced but not independently verified" is the absence of a study, not weak evidence |

This is recorded rather than fixed. The frozen vocabulary is new in this round
and the gate's strings predate it; changing them is a status-string change to
the gate and belongs with **R-12**, under the same decision record, not to a
documentation pass. Until then, `INTENDED_USE_AND_ESTIMANDS.md` §2 is
authoritative for what may be claimed, and the gate's strings are read as its
own legacy labelling. The gate's **verdict** — NOT PASSED, 0/6 — is not in
dispute under either vocabulary.

## 9. Round-4 record

This round produced specifications, designs and corrections — no measurements
and no new validation record.

- `docs/INTENDED_USE_AND_ESTIMANDS.md` — new. Claims frozen; volume estimand
  resolved; nuclear and viability policies defined.
- `docs/evidence/2026-09-13-volume-audit/ERRATA.md` — five defects in the
  round-3 documents, with the superseded wording quoted.
- `docs/evidence/2026-09-13-volume-audit/VOLUME_QUALIFICATION_PROTOCOL.md` §6 —
  the volume study design that follows from the estimand decision.
- `docs/evidence/2026-09-13-volume-audit/STUDY_PROPOSALS_SG3_SG6.md` §1.4–1.6,
  §2.1.1, §2.3.1, §3 — segmentation addenda, derived error propagation and
  tiered calibration criteria, and the SG-5 design.
- `docs/GATE_SEMANTICS_ANALYSIS.md` §7 — Option C specified item by item.
- `docs/evidence/2026-09-13-volume-audit/consistency_recheck.csv` — the
  re-derived stratum table both sets were checked against.

The scientific validation gate's verdict and every item status are unchanged by
this round.
