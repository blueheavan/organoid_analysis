# Errata — round-3 documents (issued 2026-09-13, corrected 2026-09-13)

Five defects found in the round-3 documents during the round-4 consistency
review. Each entry states what was claimed, why it was wrong, what the evidence
actually supports, and where the correction now lives. Nothing here changes any
algorithm, criterion, threshold, domain constant or gate status: E-1 to E-3 and
E-5 are claims *about* the evidence, and E-4 is a framing claim. Both analytical
items remain **FAIL / NOT QUALIFIED** on their unrestricted criteria, exactly as
before.

The original text is not overwritten in place silently — the corrected passages
carry back-references to this file, and the superseded wording is quoted here in
full so the record of what was asserted survives.

---

## E-1 — `rho_in` was defined against the wrong quantity

**Claimed** (`VOLUME_QUALIFICATION_PROTOCOL.md` §2):

> `rho_in` is the inscribed radius divided by the voxel diagonal — a resolution
> measure, not a size measure: it already folds in spacing and anisotropy.

**Wrong because** the implementation divides by the largest single spacing, not
the diagonal:

- `src/organoid_analysis/quantification/surface_crofton.py:448` —
  `rho_in = r_in / max(spacing)`
- record harness,
  `docs/evidence/2026-09-13-analytical-geometry-record/surface_vv.py:264` —
  `rho = min_axis / max(spacing)`, recorded per case as `rho` alongside
  `min_semi_axis_um` (line 276)
- method-development harness,
  `docs/evidence/2026-09-12-surface-method-development/code/common.py:229` — the
  same expression; `code/common2.py:10` states the identity explicitly
  (`r_inscribed == rho * max(spacing)`)

The two differ by a factor of √3 at isotropic spacing, and by a
spacing-dependent factor otherwise, so every domain constant expressed in
`rho_in` (`DOMAIN_RHO_IN_MIN = 10`) would be misread by anyone taking the
document's definition. The second half of the sentence is also wrong as a
consequence: a quantity keyed to the coarsest axis alone does **not** fold
anisotropy in — see E-3.

**Corrected in** `VOLUME_QUALIFICATION_PROTOCOL.md` §2, with both source
citations. The numbers in the document were computed by the harness and are
unaffected; only the prose definition was wrong.

---

## E-2 — the shape-class generalisation was overstated

**Claimed** (§2.1, observation 1):

> **Resolution, not shape class, sets the error size.** At matched `rho_in` the
> three families are interleaved […]

**Wrong because** shape class is detectable at matched resolution even inside
the three sampled families. Re-derived from the frozen grid (n = 188,
voxel-count volume, identical across all six recorded estimators — verified):

| model of `log10 |volume error|` | R² |
|---|---|
| `log10 rho_in` | 0.594 |
| `+ anisotropy indicator` | 0.636 |
| `+ shape family` | 0.677 |

The sphere term is significant (p = 0.002); ellipsoid p = 0.10. Median
`|error|` at matched `rho_in` = 3 / 6 / 12: spheres 0.84 / 0.27 / 0.07%,
ellipsoids 0.28 / 0.15 / 0.04%. A mechanism consistent with the sign: `rho_in`
is keyed to the *smallest* semi-axis, so at matched `rho_in` an elongated object
carries more volume per unit of boundary layer and a smaller relative
quantization error.

What the evidence does support: resolution is the **dominant** term among the
sampled families, and no statement at all can be made about families the grid
does not contain (hollow, branched, touching, border-truncated).

**Corrected in** `VOLUME_QUALIFICATION_PROTOCOL.md` §2.1 observation 1.

---

## E-3 — the anisotropy claim had the wrong direction and the wrong reason

**Claimed** (§2.1, observation 3):

> Anisotropy up to 3:1 adds no separate effect once `rho_in` is fixed […] as
> expected from `rho_in` being defined against the voxel diagonal.

**Wrong in three ways.** The reason given depends on E-1's wrong definition. The
supporting comparison did not hold `rho_in` fixed: it paired cases on shape,
size and variant, but since `rho_in` is keyed to the coarsest axis, an
anisotropic voxel *halves* `rho_in` for the same physical object, so the pairs
differed in resolution. And at genuinely matched `rho_in` the effect is present
and sizeable, in the opposite direction to "no effect":

| matched `rho_in` | isotropic median (n) | anisotropic median (n) |
|---|---|---|
| 3 | 0.970% (11) | 0.370% (22) |
| 6 | 0.356% (11) | 0.129% (22) |
| 12 | 0.084% (11) | 0.029% (25) |

With `log10 rho_in` in the model, the anisotropic indicator carries −0.46 in
`log10 |error|`, p < 0.001 — roughly one third of the isotropic error at the
same `rho_in`. Mechanism: at matched `rho_in` the two finer axes of an
anisotropic voxel are better resolved than `rho_in` implies, so the
boundary-voxel fraction is smaller. `rho_in` is therefore a **conservative**
resolution proxy under anisotropy on these families. Sampled anisotropies are
2, 2.857 and 3; nothing here bounds behaviour beyond 3.

**Corrected in** `VOLUME_QUALIFICATION_PROTOCOL.md` §2.1 observation 3. No
domain constant changes: the anisotropy ceiling is retained, now for stated
reasons (nothing bounds anisotropy > 3; the ceiling also constrains the area
estimator on its own evidence) rather than as a claimed non-effect.

---

## E-4 — "lower bound" was a mathematical claim the evidence cannot support

**Claimed** (§4, item 4):

> The analytical volume error is a *lower bound* on production volume error.

**Wrong because** the analytical error and the segmentation error are not
constrained to have the same sign. A segmenter that systematically under-fills a
boundary can cancel a positive quantization error, so production error may be
*smaller* than the analytical error for a given object. "Lower bound" asserts an
inequality that holds for no shown reason.

**Corrected in** §4 item 4: the analytical error is one component of the
production error budget, alongside segmentation boundary error and calibration
scale error, which may add, partially cancel, or dominate. The conclusion the
passage was there to protect — that analytical geometry evidence is not evidence
about segmented objects — is unchanged and is stated without the false
inequality.

---

## E-5 — the calibration propagation figures were wrong

**Claimed** (`STUDY_PROPOSALS_SG3_SG6.md` §2.1 and `ROUND3_REPORT.md`):

> a 3% lateral scale error is a ~3% area error and a ~9% volume error

**Wrong because** a lateral-only error scales two axes, not one and not three.
With `s_i = (1 + e_i) · s_i^true`, voxel-count volume scales as the product of
the three spacings, so for `e_y = e_x = e`, `e_z = 0`, volume goes as
`(1+e)² − 1 ≈ 2e` = **+6.09%** at e = 3%, not 9%. Surface area has no single
factor — a surface element spanned by axes `i` and `j` scales by
`(1+e_i)(1+e_j)` — so the area error lies between `e` and `2e`, i.e. **+3.00%
to +6.09%**, and is +4.04% for a sphere (exact oblate-spheroid formula,
semi-axes `(1+e)R, (1+e)R, R`). The originally quoted ~3% area / ~9% volume pair
corresponds to no single error model: 9% volume requires a 3% error in all three
axes, and that same isotropic error gives +6.09% area, not 3%.

**Corrected in** `STUDY_PROPOSALS_SG3_SG6.md` §2.1, which now derives the
propagation and tabulates both the lateral-only and isotropic cases, and in
`ROUND3_REPORT.md`. The qualitative conclusion is unchanged and if anything
strengthened: at a 3% miscalibration the systematic term is several times the
qualified estimator's worst in-domain area error (0.95% on the confirmation
set), so calibration cannot be deferred behind estimator work.

**Refined on the same date (round 4).** The first correction gave the sphere's
lateral-only area sensitivity as `≈1.35e` (+4.04% at e = 3%). The exact value
from the oblate-spheroid area formula is `(4/3)e`, i.e. **+4.02%** at e = 3%;
`(4/3)` is exact, and `4/3 = (2/3)·2` is why sphericity is exactly invariant for
a sphere. The shape dependence of the sensitivity, the sphericity invariance
results, and the invariance of `rho_in` to isotropic scale error were absent
from both earlier forms and are now derived in `STUDY_PROPOSALS_SG3_SG6.md`
§2.1.1. The acceptance numbers in §2.3 were illustrative; §2.3.1 now derives a
per-axis tolerance from the estimator's own error budget and tiers the criteria
by the claim each licenses.

---

## E-6 — the SG-6B → real-object-morphology dependency was overstated

**Claimed** (`ROADMAP.md` §5 and `SCIENTIFIC_VALIDATION_MASTER_PLAN.md` §11
dependency diagrams):

> `SG-6B absolute calibration ---> real-object morphology`
> `SG-3 --> real-object morphology`, with SG-6 as a general upstream dependency

**Wrong because** the SG-3 downstream-bias metric is a *ratio of two masks
measured on the same voxel grid*. With `V = s_z s_y s_x N`, the relative
discrepancy `(V_prod − V_ref) / V_ref = (N_prod − N_ref) / N_ref` is exactly
independent of the spacings — for any diagonal spacing error, isotropic or not
— because the spacing product is a common factor of numerator and denominator.
Measured on two concentric spheres standing in for production and reference
(`docs/reviews/spacing_invariance_check.py`): the volume discrepancy is
bit-identical across no error, +3% common, ×2 common, +3% axial-only and +20%
axial-only; the area discrepancy moves by 0.007 percentage points even at a 20%
axial error. Absolute scale therefore does not contaminate SG-3's bias metric.
What `rho_in` and the anisotropy flag respond to is the z:xy *ratio*, so the
real dependency is `SG-6A ratio -> SG-3`, and only through stratum assignment
and the domain flag. `SG-6B` (absolute calibration and registration) governs
*absolute* µm, µm² and µm³ claims — it must be evidenced before physical-unit
reporting, but it is not on the SG-3 critical path.

**Corrected in** `ROADMAP.md` §5 and `SCIENTIFIC_VALIDATION_MASTER_PLAN.md`
§11, whose diagrams now show `SG-6A -> SG-3` for stratification and route
`SG-6B` to absolute-unit claims only. `ROADMAP.md` Track B already stated
"Final domain stratification requires SG-6A, not SG-6B"; the diagrams now agree
with it. This is a documentation correction only: SG-3A/SG-3B remain
`INSUFFICIENT EVIDENCE`, SG-6A/6B unchanged, and no threshold moved.

---

## Consistency items checked and found already correct

Recorded so the review is not repeated:

- **Declared domain vs. confirmation coverage.** The declared domain is
  `rho_in ≥ 10` with no upper bound, anisotropy ≤ 4, smooth closed surfaces. The
  confirmation set covers `rho_in` 10.43–15.05 and anisotropy
  {1.2, 2.2, 3.5, 4.0} for the 96 in-scope smooth cases inside the machine gate.
  Everything above `rho_in` 15.05 is declared but not confirmation-covered; the
  documents state this, and `VOLUME_QUALIFICATION_PROTOCOL.md` §4 item 5 records
  that whether real acquisitions reach `rho_in ≥ 10` at all is unevidenced.
- **Machine-checkable domain vs. scope assumption.** `surface_in_qualified_domain`
  is computed from `rho_in` and anisotropy only; the smooth-closed-surface clause
  is carried separately in `surface_qualification_scope` and cannot be checked
  from a mask (`features.py:85-86`, `crease_indicator_rejected.md`). Inside the
  machine-checkable gate the creased classes (box, cylinder; n = 49) reach 3.20%
  volume and 9.10% area error against 0.76% and 0.95% for the smooth in-scope
  classes — the criterion is honoured by the non-machine-checkable clause, not by
  the code.
- **Packaged vs. solved weights.** Every spacing in both evidence sets resolves
  to a packaged table (verified by running the selector on all nine: dev
  1×1×1, 2×0.7×0.7, 2×1×1, 3×1×1; confirmation 1.2×1×1, 2.2×1×1, 2.4×0.6×0.6,
  3.5×1×1, 5×1×1). The packaged set is 11 discrete z:xy ratios
  {1, 1.2, 1.5, 2, 2.2, 2.857, 3, 3.125, 3.5, 4, 5} with equal lateral spacings.
  A realistic acquisition at 1 × 0.325 × 0.325 µm falls outside it and is solved
  at run time, as does any dataset whose coarse axis is not the first
  (e.g. 1×1×2). `PROVENANCE.md` records that the minimax program has a
  non-unique optimum, so a run-time solve returns a different weight vector with
  the same worst-case plane-response deviation. The per-object column
  `surface_weights_origin` (`features.py:214`) already distinguishes
  `packaged` / `cache` / `solved`, and `measurement_policy.py` already states
  that non-packaged weights are not the frozen evidence-bearing vectors. This is
  a **limitation of the evidence's coverage**, not a defect in the code, and is
  carried into the round-4 roadmap as an open item.

## Reproduction

`consistency_recheck.csv` in this directory holds every number quoted above.
Inputs, as named by `docs/evidence/analytical_geometry_current_record.json`:

- `docs/evidence/2026-09-13-analytical-geometry-record/surface_vv_dev.csv` —
  1128 rows = 188 cases × 6 recorded estimators (frozen grid)
- `docs/evidence/2026-09-12-surface-method-development/vv_confirm2.csv` —
  480 rows (confirmation set)

Voxel-count volume is estimator-independent; that was verified on the record
(every case's `volume_rel_err` identical across all six estimators) rather than
assumed, and one estimator's 188 rows were then used for the volume statistics.
