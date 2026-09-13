# Voxel-count volume: evidence audit and qualification protocol

**Status: FAIL / NOT QUALIFIED (SG-2).** This document audits the evidence that
already exists for the volume estimator and specifies what would have to be
executed to qualify it. It declares no domain, changes no algorithm, weakens no
criterion and reports no new measurement. The volume estimator — a voxel count
times the voxel volume — is untouched by this round's work.

Audited 2026-09-13 from two existing evidence sets, both re-derived here:

- the canonical analytical record `docs/evidence/2026-09-13-analytical-geometry-record`,
  production rows only (harness label `E0_marching_cubes_binary_v1`, n = 188) —
  the frozen grid the science gate reads;
- the surface method-development confirmation set
  `docs/evidence/2026-09-12-surface-method-development/vv_confirm2.csv`
  (n = 480), whose acceptance rules were frozen in `FREEZE_RECORD_V3.md`
  **before** execution and which was executed once.

Tabulated in `volume_audit_by_stratum.csv`. Nothing in this document is new
confirmation evidence; both sets already existed.

## 1. Why volume is a separate scientific problem from surface area

Volume and surface area are separate estimators with separate failure
mechanisms, and this round's surface-area work transfers no confidence to
volume:

- Volume is a **voxel count**. Its error is a partial-volume/quantization error
  at the boundary: whole voxels are counted in or out, so the magnitude is set
  by the boundary-voxel fraction and the *sign* by where the object sits
  relative to the voxel lattice.
- Surface area is a boundary-integral estimator whose error is dominated by the
  discrete direction set and by curvature relative to the voxel. The
  2026-09-12 work changed the surface estimator's weights, direction set and
  stencil. **None of that touches a voxel count** — per-case volume values are
  identical across every estimator row in the record.
- The criteria also differ: SG-1 is `|error| < 5%`, SG-2 is `|error| < 1%`.
  Volume is held to the tighter criterion.

## 2. What the existing evidence shows

`rho_in` is the inscribed radius divided by the voxel diagonal — a resolution
measure, not a size measure: it already folds in spacing and anisotropy.

### 2.1 The frozen grid the gate reads (n = 188)

| `rho_in` stratum | n | max abs. rel. error | median signed | cases ≥ 1% | worst case | vs. SG-2 |
|---|---|---|---|---|---|---|
| very_small (<3) | 44 | **18.26%** | +0.74% | 68.2% | dev-078-cylinder-size3-sp2x1x1-rot0 | FAIL |
| small (3–6) | 44 | **2.74%** | −0.07% | 15.9% | dev-000-sphere-size3-sp1x1x1-offset0 | FAIL |
| medium (6–12) | 44 | **2.57%** | −0.01% | 4.5% | dev-035-cylinder-size6-sp1x1x1-rot0 | FAIL |
| large (≥12) | 56 | 0.53% | +0.00% | 0.0% | dev-090-cylinder-size24-sp2x1x1-rot0 | PASS |

The unrestricted SG-2 status is **FAIL**: the criterion is every analytical case
below 1%, and three of four strata contain cases above it. Two mechanism
observations, both descriptive:

1. **Resolution, not shape class, sets the error size.** At matched `rho_in` the
   three families are interleaved; the per-`rho_in` envelope of `|error|` falls
   steeply and near-monotonically (log–log slope of the envelope on this grid
   −1.78). All three families reach 15–18% at `rho_in < 2` and all sit below
   0.6% at `rho_in ≥ 12`.
2. **Grid phase sets the sign.** At `rho_in < 3` one shape gives −18.3% at one
   sub-voxel offset and +7.3% at another; the stratum median is +0.74% while
   individual cases reach ±18%. This is a per-object quantization error, not a
   systematic bias a single correction factor could remove.
3. Anisotropy up to 3:1 adds no separate effect once `rho_in` is fixed
   (max 0.12–0.53% across anisotropy 1, 2, 2.86, 3 at `rho_in ≥ 8`), as expected
   from `rho_in` being defined against the voxel diagonal.

### 2.2 Prospective in-domain evidence that already exists (n = 480)

This matters and must not be understated. `FREEZE_RECORD_V3.md` predeclared
`A2: |volume relative error| < 1%` as an acceptance rule, alongside the area and
sphericity rules, for the declared domain (`rho_in ≥ 10`, anisotropy ≤ 4,
**smooth** closed surfaces). The confirmation set was then executed once:

| subset of the confirmation set | n | worst \|volume err\| | vs. 1% |
|---|---|---|---|
| in scope and inside the gate (sphere, ellipsoid, capsule, torus) | 96 | **0.76%** | met |
| creased classes inside the same gate (box, cylinder) | 49 | **3.20%** | exceeded |

Per class, in scope: sphere 0.11%, ellipsoid 0.08%, capsule 0.27%, torus 0.76%.
So **voxel-count volume does have predeclared, executed-once, in-domain evidence
at better than 1%** for four smooth classes. Any statement that volume rests on
no prospective evidence would be wrong.

![Volume error against resolution, grid phase, and shape class](volume_error_audit.png)

## 3. Why the status is nonetheless FAIL / NOT QUALIFIED

1. **SG-2's criterion is unrestricted.** It asks for every analytical case below
   1%, and the frozen grid's worst is 18.26%. Whether that is the right question
   for the gate to ask is a separate decision, analysed in
   `docs/GATE_SEMANTICS_ANALYSIS.md`; under the rule as it currently stands the
   status is FAIL, and it is reported as FAIL.
2. **The machine-checkable gate does not bound volume to 1%.** This is the
   substantive finding of the audit. Inside `rho_in ≥ 10` and anisotropy ≤ 4 —
   everything software can verify from a mask — creased confirmation cases reach
   3.20%. The 0.76% figure holds only because the *smoothness* clause of the
   domain excludes them, and `in_domain()` cannot check smoothness: a
   mask-computable crease indicator was tested and rejected before the
   confirmation run. So an out-of-scope object can satisfy every checkable
   condition and still receive a volume error three times the criterion. The
   existing per-object `surface_outside_qualified_domain` flag does not catch
   it, because such an object is inside the checkable domain.
3. **Class coverage is incomplete for volume.** Neither set contains a hollow or
   shell object, a branched shape, a touching pair or a border-truncated object.
   Volume's error mechanism depends on the boundary-voxel fraction, which a
   lumen changes directly at fixed `rho_in` — a shell has two boundaries.
4. **A phantom is not a segmentation.** Both sets are exact rasterizations of
   analytical solids. A production mask is a segmenter output whose boundary
   error is unknown and is assessed separately (SG-3, currently
   FAIL / INSUFFICIENT EVIDENCE). The analytical volume error is a *lower bound*
   on production volume error. Analytical geometry evidence is not evidence
   about segmented objects, and this document does not offer it as such.
5. **Whether real acquisitions reach `rho_in ≥ 10` is unevidenced** (SG-6). A
   domain the data cannot enter would qualify nothing.

**One observation deliberately recorded as a hypothesis, not a result.** On the
frozen grid every production case with `rho_in ≥ 8` is below 1% (n = 67, max
0.53%). That threshold was found by scanning the recorded table for the smallest
cut that passes — hypothesis generation on data already seen, with no error
control. It is listed in `volume_audit_by_stratum.csv` under rows labelled
*exploratory* and **must not be quoted as a resolution guarantee.** Per §3.2 the
binding problem is not the threshold anyway: it is that the checkable part of
the domain does not exclude the shapes that fail.

## 4. Qualification protocol (to be executed; nothing here has been run)

Volume can be qualified for a stated domain by executing the following in order.
Every choice in steps 1–3 is frozen, in a commit, before any confirmation
measurement.

**Step 1 — Declare intended use and a candidate domain in volume's own terms.**
The existing domain constants were selected for the surface estimator (the
source notes `rho_in ≥ 10` was set by the volume criterion on *those*
phantoms); reusing them for volume without restating them is the retrospective
redefinition the brief forbids. The declaration must state shape families
including at least one hollow/shell class, the `rho_in` range, the anisotropy
range, whether border-truncated objects are in scope, and — the point of
§3.2 — **which conditions a caller can verify from a mask alone.** If
smoothness remains load-bearing for volume, the declaration must say that the
domain is not machine-checkable and that responsibility rests with the caller,
exactly as the surface domain does.

**Step 2 — Freeze the acceptance rule.** `|relative volume error| < 1%` for
**every** case inside the declared domain, matching SG-2. Predeclare the
handling of non-estimable cases and predeclare that a single in-domain
violation is a FAIL for the whole domain. No margin, no percentile relaxation,
no per-stratum exemption.

**Step 3 — Build two disjoint phantom sets from one seeded generator.**
   - A **development set** for all exploration, extending the existing classes
     with what §3.3 lacks: hollow spheres and shells across wall thickness,
     branched shapes, touching pairs, border-truncated objects, each at several
     sub-voxel offsets and rotations. Offsets must be sampled densely, since
     §2.1 shows they set the sign.
   - A **confirmation set** from the same generator under a different frozen
     seed, sampling the declared domain independently and including its
     boundary, written and hashed **unmeasured**. Executed exactly once; any
     change of estimator, domain or rule afterwards voids it and requires a new
     set.
   - Both produced by a recorded generator with pinned seeds, per
     `docs/VALIDATION_RECORDS.md`.

**Step 4 — Decide, before execution, whether a crease/lumen indicator is
required.** If the domain's checkable conditions do not exclude the failing
classes, either a mask-computable indicator must be found and frozen on the
development set (the area indicator was rejected on the evidence in
`crease_indicator_rejected.md`; volume needs its own assessment, since volume
error is insensitive to creases in a way area error is not and the indicator
that failed for area may not be the right one here), or the domain must be
declared non-machine-checkable with that limitation stated in the export
contract.

**Step 5 — Execute once, record, and report in the fixed vocabulary.** Extend
the analytical record contract with the volume confirmation result (n, worst
`|error|`, median signed error, per-stratum and per-class breakdown, worst case)
so the gate reads SG-2 from a record rather than from a document, and report
PASS only if every in-domain confirmation case is below 1%; otherwise FAIL, with
the domain recorded as not qualified. A domain narrowed *after* seeing the
confirmation result is a new hypothesis needing a new confirmation set, and must
be reported as such rather than as a PASS on a narrower domain.

**Out of scope for this protocol.** Anything about segmented objects. Volume
accuracy on production masks requires the SG-3 study; analytical qualification
is a precondition for that study's downstream-bias analysis, not a substitute
for it.

## 5. Decisions required from the project owner

1. **Is the volume domain the same domain as the surface domain?** §2.2 shows
   smooth in-scope cases meet 1% while creased in-gate cases do not, so the
   scope clause carries volume as well as area. If volume is to be reported for
   creased objects at any resolution, a new development study is needed before
   any confirmation set is built.
2. **Analytical qualification now, or only alongside SG-3?** Step 3's phantom
   extension is the bulk of the work and needs no imaging data.
3. **Is the 1% criterion the biological requirement, or inherited?** A criterion
   stricter than needed will fail a scientifically adequate domain — but it may
   only be loosened **before** the confirmation set is measured, with the
   biological justification recorded.
