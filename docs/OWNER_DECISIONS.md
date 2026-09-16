# Owner decisions — consolidated register

Status: **ALL FOURTEEN DECISIONS RECORDED 2026-09-16.** Every outstanding
`[OWNER]` placeholder and every "owner must set this before data are seen"
criterion scattered across the validation documents is now resolved below, with
its derivation, the document that carries it, and the condition that reverses
it.

Provenance: the fourteen decisions were recorded under the project owner's
explicit delegation to the analysis agent (`wisp-science`) on 2026-09-16, on the
literature basis in `docs/LITERATURE_BASIS_FOR_DECISIONS.md`. They are
**predeclared criteria and study designs**, recorded before any confirmation,
annotation or instrument data are inspected. They adopt no measurement, report
no study result, and move **no gate status**: the gate keeps its legacy
fail-closed behaviour and remains `NOT PASSED (0/9)`.

Reversal: any entry may be superseded by the project owner with a dated entry in
this table **recorded before the data the corresponding study requires are
seen**. Each row states what that reversal would cost.

Consolidated 2026-09-16. Amended by the same pass as before: item **D-6**
(SG-3 downstream threshold) and item **D-2** (SG-6B dependency) reflect the E-6
correction in `docs/evidence/2026-09-13-volume-audit/ERRATA.md`.

## Recorded decisions

| ID | Decision (recorded 2026-09-16) | Derivation | Recorded in | Reversal condition |
|---|---|---|---|---|
| D-1 | **Define the estimand per topology; exclude no object.** `V_env` primary for all objects; `V_seg` and `f_void` reported per object; objects whose void connects to the exterior carry the flag `open_cavity` and report `f_void = 0 (by construction)`, never a lumen measurement. `open_cavity` is a declared stratum, reported with and without. | `INTENDED_USE_AND_ESTIMANDS.md` §3.3 reads the discontinuity as a reporting defect, and the estimand is the outer envelope, which cavity topology does not perturb. Removing objects would change the population in a size-correlated way. | `INTENDED_USE_AND_ESTIMANDS.md` §3.3; `VOLUME_QUALIFICATION_PROTOCOL.md` §5 item 1; `ROADMAP.md` §6 | If a development study shows the open/closed classification is unstable between time points of the same organoid, the policy narrows to reporting `V_seg` as primary for that stratum. |
| D-2 | **The volume domain is the surface domain** — smooth/closed, `rho_in ≥ 10`, anisotropy ≤ 4 — extended only by D-1's `open_cavity` stratum. Creased and hollow objects are out of domain for volume qualification. | The confirmation evidence already separates them: smooth in-scope worst 0.76% volume versus 3.20% for creased in-gate cases. A looser second domain would need its own development and confirmation set. | `VOLUME_QUALIFICATION_PROTOCOL.md` §5 item 1; `ROADMAP.md` §6 | A development study that qualifies a creased volume domain with its own confirmation set. The domain widens by a record, never by argument. |
| D-3 | **Resource the analytical volume qualification now**, in parallel with SG-3, as the first item of Track A. | It needs no imaging data, and analytical qualification is a *precondition* for SG-3's downstream-bias analysis. Sequencing it behind SG-3 puts a zero-specimen task behind an annotation-limited one. | `VOLUME_QUALIFICATION_PROTOCOL.md` §5 item 2; `ROADMAP.md` §3 Track A, §6 | Only a resourcing decision; no scientific consequence either way. |
| D-4 | **`<1%` is an inherited engineering target, retained unchanged** as the SG-2a criterion, and explicitly **not** the criterion for SG-3. | No biological requirement looser than `<1%` has been derived, and the intended use (IU-1, within-experiment comparison) does not require absolute single-object accuracy. Retaining it is conservative; SG-3 is governed by D-5/D-6 instead. | `VOLUME_QUALIFICATION_PROTOCOL.md` §5 item 3; `SCIENTIFIC_VALIDATION_MASTER_PLAN.md` §1C | Loosening is permitted only before a confirmation set is measured, with the biological justification recorded. |
| D-5 | **δ predeclared as 0.20 relative for volume and surface area, 0.10 for diameter, sphericity and principal-axis length; acceptance C1–C5 per matched object with uncertainty at the biological-sample level; N from precision, not precedent.** C1 recall and precision lower 95% bound ≥ 0.90; C2 upper 95% bound on median absolute relative bias < δ/4; C3 upper 95% bound on its 95th percentile < δ/2; C4 catastrophic-error rate (bias > 0.5) upper 95% bound < 0.02; C5 limits of agreement are the 95th-percentile bound, always reported with intervals. N: ≥ 47 biological samples per stratum (≈ 940 objects, design effect 6.7 at `m = 20`, ICC 0.30), ≥ 20 samples per compared condition, bias resolvable to ±δ/8. | δ is a predeclared convention — the smallest between-condition change treated as actionable — set conservatively because within-culture organoid heterogeneity is the operative limit on resolution. N follows from C1's precision with a predeclared design effect. | `SEGMENTATION_VALIDATION_PROTOCOL.md` §5; `SCIENTIFIC_VALIDATION_MASTER_PLAN.md` §8 | δ and the criterion set may be revised before the confirmation set is measured. Any revision after it invalidates the study. |
| D-6 | **The inherited `1% / 5%` pair is not SG-3's criterion.** SG-3 is governed by D-5's C1–C5, **including the differential criterion**: for any two compared conditions, the 95% CI of the difference in median relative bias must exclude ±δ/8 (2.5%), for volume and for surface area. SG-1a/SG-2a keep `1% / 5%` as estimator criteria. | A bias common to all objects cancels in a ratio; only a bias that differs between arms threatens the between-condition contrast, and no absolute per-case threshold detects it. The two numbers answer different questions. | `VALIDATION_PLAN.md`; `SCIENTIFIC_VALIDATION_MASTER_PLAN.md` §8; `ROADMAP.md` §6 | None short of a change in intended use: the differential criterion is what protects the claim. |
| D-7 | **Per-class lower 95% CI bound on sensitivity and precision ≥ 0.80** for each of the four states, and macro-F1 lower bound ≥ 0.80; confusion matrix primary; **quadratic-weighted kappa and Spearman reported but non-gating**; classes with prevalence < 5% reported `not evaluable`, not FAIL. **Indeterminate ceiling `f_indeterminate ≤ 0.20` per condition, and the between-condition difference in `f_indeterminate` ≤ 5 pp on its upper 95% bound.** Above either, V2 is not reportable for that condition. | Kappa is incoherent under class imbalance, so the confusion matrix and class-specific bounds carry the claim. The differential half is load-bearing because `f_viable_like` divides by `N_classifiable`: an indeterminate rate differing between arms biases the contrast even when both arms are under the ceiling. | `INTENDED_USE_AND_ESTIMANDS.md` §5.5 Q4; `SCIENTIFIC_VALIDATION_MASTER_PLAN.md` §9; `ROADMAP.md` §6 | The numeric bounds may be revised before the study's data are seen. |
| D-8 | **The primary object-level reference is a blinded human classification of the same objects**, on a predeclared stratified subsample, by two annotators, with the reference's own disagreement reported. The ATP assay remains the secondary well-level reference. **R-13 stays out of scope.** The claim V1 licenses becomes *agreement with blinded human classification of calcein/PI signal*, bounded by the reference's accuracy, with marker specificity established separately by the graded-insult series. | It is the only exit that removes B2's circularity without re-qualifying the production nuclear branch that V1 already depends on. It reuses SG-3's annotation campaign. Demoting to ATP-primary would leave the per-object states permanently `INSUFFICIENT EVIDENCE`. | `INTENDED_USE_AND_ESTIMANDS.md` §5.5 Q1; `SCIENTIFIC_VALIDATION_MASTER_PLAN.md` §9; `ROADMAP.md` §6 | Requires new evidence that an independent automated path exists; not a preference. |
| D-9 | **Satterthwaite degrees of freedom for the primary fixed-effect contrasts; between-within is used only for the omnibus test**, which is reported with a between-within df check and an explicit anti-conservative-at-`G < 10` statement. Kenward–Roger remains a sensitivity analysis where a verified implementation exists. | The primary claim is a per-contrast claim; Satterthwaite is a per-contrast approximation, while between-within is an omnibus df defined for a balanced nested design. | `SCIENTIFIC_VALIDATION_MASTER_PLAN.md` §10; `ROADMAP.md` §6 | The D-10 coverage study may overturn it cell by cell; the freeze changes only with a recorded re-freeze. |
| D-10 | **Coverage tolerance 1 percentage point from nominal, against a simulation MCSE of at most 0.33 pp** (`n_sim = 5 000` per cell). Type-I error at nominal 0.05 must lie in [0.04, 0.06]. An interval may be reported only for cells passing T2 coverage; the shipped floor of 3 replicates is retained as a computation floor but is **provisional**, so any interval at 2–3 replicates is labelled `provisional — coverage unqualified` and carries no claim. The study may raise the reportable floor, never lower it. | The tolerance must be wide enough to be a statement about the method rather than about the simulation, which is what fixes `n_sim` relative to MCSE. | `SCIENTIFIC_VALIDATION_MASTER_PLAN.md` §10; `ROADMAP.md` §6 | A recorded re-freeze with a recalculation of MCSE at a different `n_sim`. |
| D-11 | **Commission the local CR2 implementation and verify it against the published reference implementation and published known cases.** Verification: golden values from `clubSandwich` with recorded version; analytic cases where CR2 must equal CR0; the published cases and df values of Pustejovsky & Tipton 2016; then the D-10 simulation on the frozen implementation only. | The pinned environment has no CR2 (`statsmodels` 0.14.6 offers only CR0-style `cov_type="cluster"`), so the choice is between unverified and verified local code, and only the latter can carry a coverage claim. | `SCIENTIFIC_VALIDATION_MASTER_PLAN.md` §10 | Not reversible without abandoning the coverage claim for the fallback branch. |
| D-12 | **`0.25%` per axis is the total calibration contribution and is split, not reused.** Volume budget `e_x + e_y + e_z ≤ 0.76%`; equal allocation gives 0.25%. Split by axis and tier: **T1 stability** lateral ≤ 0.10%, axial (ratio stability) ≤ 0.25%, both gated on the upper 95% bound, monitoring ±2·`e_stab`; **T2 accuracy** lateral ≤ 0.15%, axial **not qualified** and instead reported as a stated systematic uncertainty ±X%; combined per axis in quadrature ≈ 0.18% (lateral); **T3 transferability** ≤ 0.8% per axis, derived from the contrast budget δ/8 = 2.5%, axial reported not gated. | The qualified in-domain worst case is 0.76% for volume and volume responds to the sum of the affected axes' scale errors. Independent contributions combine in quadrature. Lateral and axial separate because the axial term is dominated by refractive-index-mismatch aberration, which is configuration-dependent. | `STUDY_PROPOSALS_SG3_SG6.md` §2.3.1; `ROADMAP.md` §6 | The tier values may be revised before the standard is imaged; the split principle may not be collapsed back into a single reused number. |
| D-13 | **Initial re-qualification interval 6 months**, extending to 12 if the first interval shows drift within T1; the observed drift governs the next interval. **Event triggers** force T1 re-verification before further physical-unit claims: change of objective, immersion or mounting medium, Z step, camera, filter cube or acquisition software version; instrument move; any monitoring point outside the control limits. Ongoing monitoring of a standard specimen at predeclared frequency against the qualified T1 value. | Interval-setting is an inference from observed drift and reliability with conservative first intervals, applied to an instrument whose axial behaviour is configuration-dependent. | `STUDY_PROPOSALS_SG3_SG6.md` §2.4; `ROADMAP.md` §6 | The interval is a property of the instrument configuration; a recorded drift history is sufficient to revise it. |
| D-14 | **The Option-C domain-restricted qualification record contract is authorised**, with the interface fixed: `record_class ∈ {canonical, scope-clean}` plus a per-item `results` map keyed by gate item id carrying `status`, `domain` and confirmation-evidence hashes; PASS only against domain-restricted evidence with the domain constants hash-bound; `SG-1b`/`SG-2b` characterization with `counts_toward_gate = false`; one commit carrying contract, adversarial tests and this record's citation. Implementation is **blocked only by commit authority**, because `VV_TOOLING` hash-binds `analytical_geometry_evidence.py` and editing it invalidates the current record until a new one is generated against a clean committed tree. | `GATE_SEMANTICS_ANALYSIS.md` §6 requires the contract, the adversarial tests and the decision citation to land together. | `GATE_SEMANTICS_ANALYSIS.md` §6; `ROADMAP.md` P0-4 | Withdrawal returns the gate to legacy fail-closed behaviour, which is where it already is. |

## Decisions already recorded (context, not open)

These are listed only so the register is not mistaken for the whole decision
history. They are frozen; do not reopen them here.

| Decision | Recorded |
|---|---|
| Gate architecture Option C (qualification vs characterization) | `GATE_SEMANTICS_ANALYSIS.md`, selected 2026-09-13; contract authorised as D-14 |
| Primary organoid volume estimand `V_env`; `V_seg`/`f_void` companions | `INTENDED_USE_AND_ESTIMANDS.md` §3 |
| Viability states are `*_like`; V2 aggregate estimand and denominators | `INTENDED_USE_AND_ESTIMANDS.md` §5 |
| SG-1 surface claim bounded to smooth/closed resolved domain | `SCIENTIFIC_VALIDATION_MASTER_PLAN.md` §5 |
| Analytical volume classified `NOT QUALIFIED` pending its own domain | `VALIDATION_UPDATE_2026-09-13.md` |

## Non-decisions (explicitly not owner thresholds)

- The `_DEGENERATE_RANDOM_EFFECT_VARIANCE_RATIO = 1e-6` numerical-stability
  threshold and the `LOG10_FEATURES = {"volume_um3"}` transform are documented
  in `docs/PARAMETERS.md` as heuristics pending SG-5; they are not acceptance
  criteria and not open owner decisions here.
- `docs/VALIDATION_PLAN.md` VR-8 "Minimum N: TBD by available dataset" is a
  data-availability limitation, not a threshold the owner sets; it is reported
  `INSUFFICIENT EVIDENCE`.

## What recording these decisions does not do

- No gate status moves. The gate reports `NOT PASSED (0/9)` before and after,
  because every item still lacks its predeclared, independent, one-time
  confirmation record.
- No estimator, threshold in code, default, domain constant or qualification
  claim changes. The decisions are criteria and designs for studies that have
  not been run.
- No study becomes executable by being decided. D-5, D-7, D-8 and D-12 each name
  evidence that requires specimens, annotation or instrument time.
