# Gate semantics: should SG-1 and SG-2 be unrestricted or domain-restricted?

**Status: OPTION C SELECTED BY PROJECT OWNER ON 2026-09-13.** The selection
freezes the target architecture. The current gate remains fail-closed under its
legacy implementation until the qualification-record contract and adversarial
tests in §4.3 are implemented; this decision creates no PASS.

Written 2026-09-13. Every number is re-derived from existing evidence
(`docs/evidence/2026-09-13-analytical-geometry-record`, n = 188 production rows;
`docs/evidence/2026-09-12-surface-method-development/vv_confirm2.csv`, n = 480
predeclared confirmation cases). No new measurement was made for it.

## 1. The question

`scripts/scientific_validation_gate.py` evaluates SG-1 (surface area,
`|rel. err| < 5%`) and SG-2 (volume, `|rel. err| < 1%`) over **every case in the
canonical analytical record**, with no reference to the estimator's declared
applicability domain. The record's grid deliberately spans `rho_in` from about
1.2 to 40 and includes creased shapes, because it was built to characterize
behaviour across the whole range — including regions the method never claimed.

Under that rule:

| item | scope | n | worst \|error\| | criterion | status |
|---|---|---|---|---|---|
| SG-1 | every recorded case | 188 | **25.81%** | < 5% | FAIL |
| SG-2 | every recorded case | 188 | **18.26%** | < 1% | FAIL |

Both are reported FAIL today, and both statements are true. The question is not
whether they are true but **what proposition the gate is asserting**, because
the same estimator, evaluated on the domain it actually claims, has a
predeclared confirmation result that met every acceptance rule:

| subset | n | worst area | worst volume | worst sphericity |
|---|---|---|---|---|
| confirmation set, in scope, inside the gate | 96 | 0.95% | 0.76% | 0.95% |
| frozen record grid, inside the declared domain | 40 | 0.79% | 0.22% | — |
| confirmation set, creased classes inside the gate | 49 | **9.10%** | **3.20%** | — |
| frozen record grid, outside the declared domain | 148 | **25.81%** | **18.26%** | — |

So "FAIL" and "qualified for its declared domain" are both accurate, about
different propositions. That ambiguity is the defect, not either number.

## 2. What the third row of that table means

It is the reason this is a scientific question and not a presentation
preference. The declared domain has two parts:

- **machine-checkable**: `rho_in >= 10`, anisotropy <= 4 — `in_domain()` verifies
  these from the mask and spacing alone;
- **not machine-checkable**: smooth closed surfaces, no dihedral creases — no
  mask-computable indicator exists. One was tested and rejected before the
  confirmation run (`crease_indicator_rejected.md`: the marching-cubes/Crofton
  disagreement overlaps completely between classes and correlates with the
  error at r = −0.20).

Creased confirmation cases that satisfy **every checkable condition** reach
9.10% area and 3.20% volume error. A domain-restricted gate that reported PASS
would therefore be asserting a bound that the software cannot enforce and that a
caller can violate without any flag firing, because such an object is inside the
checkable domain. Any restricted option must carry that limitation explicitly;
it cannot be silent.

## 3. Options

### Option A — keep the unrestricted criterion unchanged

Both items stay FAIL. The gate asserts: *this estimator is not accurate
everywhere it can be called.* True, and an honest warning.

- *For*: no change, no possibility of over-claiming, and it keeps visible the
  fact that the software will measure any object handed to it.
- *Against*: the gate cannot pass and cannot be made to pass by any method
  development, because the grid contains cases (`rho_in ≈ 1.2`, creased) where
  no voxel estimator can meet 1%. An item that is red by construction stops
  carrying information, and it gives no credit to — and no visibility of — the
  predeclared confirmation result that does exist. It also cannot distinguish
  "not yet validated" from "validated, for a stated domain".

### Option B — restrict both items to the declared domain

SG-1 and SG-2 would read from an in-domain confirmation result and report PASS
(worst 0.95% and 0.76% against 5% and 1%).

- *For*: the gate would assert what the method actually claims, which is how a
  qualification statement is normally scoped.
- *Against*: a bare PASS would hide §2 — that the bound is not enforceable on
  arbitrary masks — and would hide that out-of-domain objects are still measured
  and exported. Losing the unrestricted number also loses the characterization
  of how bad it gets outside the domain, which is what makes the flag
  meaningful to a downstream analyst. **Not recommended in this bare form.**

### Option C — two items per quantity: qualification and characterization (recommended)

Split each quantity into two gate items with different jobs:

- `SG-1a` / `SG-2a` — **qualification, domain-restricted.** Status comes from a
  predeclared confirmation set inside the declared domain, under the frozen
  rule. PASS/FAIL. Its statement must name the domain and must state that the
  smoothness clause is not machine-checkable.
- `SG-1b` / `SG-2b` — **characterization, unrestricted.** Reports the worst
  unrestricted error as a *number with a domain breakdown*, not a PASS/FAIL.
  Its job is to keep §1's 25.81% and 18.26% in plain sight and to show that the
  failures are out of domain (all 33 cases above 5% and all 39 above 1% are
  outside the declared domain on the frozen grid).

- *For*: each item asserts one proposition. A reader sees both "qualified,
  for this domain" and "here is what happens outside it" without either hiding
  the other. It also makes the §2 limitation a required field rather than a
  footnote.
- *Against*: more items to maintain, and the gate's top-line summary needs a
  rule for combining them. It also requires real work, not renaming — see §4.
- *Cost of being wrong*: if the owner later decides qualification should be
  unrestricted, the `a` items are deleted and the `b` items become the gate;
  nothing is lost, because the unrestricted numbers were never dropped.

### Option D — Option C plus refusal to emit out-of-domain measurements

As C, but the measurement code refuses (rather than flags) to emit area/volume
for objects outside the checkable domain.

- *For*: makes a restricted PASS enforceable for the part of the domain that is
  checkable.
- *Against*: a policy change to production behaviour that would break existing
  workflows and drop data analysts may legitimately want with a flag attached;
  it still cannot enforce the smoothness clause, so it buys less than it
  appears to. Recorded for completeness, not recommended now.

## 4. What Option C would require (so the cost is not understated)

1. **The gate reads only the canonical analytical record.** The predeclared
   confirmation evidence lives in `docs/evidence/2026-09-12-surface-method-development/`,
   which is method-development output, not a validation record. A restricted
   item would require carrying that confirmation result into the record
   contract, with its own manifest and hashes, per `docs/VALIDATION_RECORDS.md`.
   **Pointing a gate item at a development directory is not acceptable** — it
   would convert development evidence into a validation record by wiring.
2. **SG-2a has no confirmation set of its own.** The volume result in §1 comes
   from a set whose domain was declared for the surface estimator. Per
   `docs/evidence/2026-09-13-volume-audit/VOLUME_QUALIFICATION_PROTOCOL.md`,
   volume needs its own declared domain and confirmation set before it can have
   a qualification item. Until then `SG-2a` would be
   **FAIL / INSUFFICIENT EVIDENCE**, not PASS — and the honest outcome of
   adopting Option C today is that the volume qualification item is red for a
   *different and more accurate reason* than it is red now.
3. **Adversarial tests.** The evidence-integrity tests must be extended so a
   restricted item cannot silently widen its domain: the domain constants a
   restricted item evaluates against must be hash-bound to the record, and a
   test must fail if a domain threshold changes without a new confirmation
   record. No existing test may be weakened to accommodate the split.
4. **Status vocabulary.** The `b` items report a number and are not PASS/FAIL,
   so the gate's vocabulary needs a `CHARACTERIZED` (or equivalent) state that
   the summary logic does not count as a pass.

## 5. Recommendation

Option C, with the honest immediate consequence stated: adopting it does **not**
turn the gate green. It would yield `SG-1a` PASS (domain-restricted, on existing
predeclared evidence, once that evidence is carried into a record),
`SG-2a` FAIL / INSUFFICIENT EVIDENCE (no volume-specific confirmation set), and
both `b` items reporting the unrestricted numbers unchanged. The overall gate
remains NOT PASSED, because SG-3 through SG-6 are unaffected by any of this.

The reason to prefer it is not the status it produces but that it separates two
questions the current single item conflates: *is the estimator qualified for
what it claims* and *what does it do outside that claim*. Both need to be
visible, and neither answer should be inferable only from the other's absence.

## 6. Decision record

Recorded from the project owner's 2026-09-13 scientific-validation directive.
Until the prerequisites below are implemented, the gate keeps its legacy
fail-closed behavior.

- [ ] Option A — keep unrestricted
- [ ] Option B — restrict to the declared domain
- [x] Option C — split qualification and characterization
- [ ] Option D — Option C plus refusal to emit out-of-domain measurements

Decided by: project owner  Date: 2026-09-13

If C or D is chosen, the implementing change must land as a single commit that
(a) carries the confirmation evidence into a record contract, (b) adds the
adversarial tests of §4.3, and (c) quotes this decision record. A change to gate
semantics without that reference is out of process.

---

## 7. Specification of Option C — 2026-09-13

Written so that signing the decision record above makes the implementation
mechanical rather than a fresh design exercise. The gate now reports the four
analytical items; qualification remains fail-closed until a dedicated
domain-restricted canonical record exists. Statuses use the vocabulary in
`docs/INTENDED_USE_AND_ESTIMANDS.md`.

### 7.1 Item definitions

| Item | Question it answers | Evidence it reads | Assertion | Status if implemented today |
|---|---|---|---|---|
| **SG-1a** qualification | Is the area estimator accurate inside the domain it claims? | confirmation set, cases with `surface_in_qualified_domain` true **and** in the declared smooth scope | every such case `\|area rel. err\| < 5%` | **INSUFFICIENT EVIDENCE** — the historical 96-case result met the frozen `<5%` criterion, but the current source-changing schema migration invalidates the old record and requires a new domain-restricted canonical record |
| **SG-1b** characterisation | What does it do outside that domain? | the full analytical grid | reports worst case, per-`rho_in`-stratum maxima and the fraction above 1% — **no verdict** | reports 25.81% worst case; no PASS/FAIL |
| **SG-2a** qualification | Is voxel-count volume accurate inside a declared volume domain? | confirmation set, in-domain in-scope cases | every such case `\|volume rel. err\| < 1%` | **INSUFFICIENT EVIDENCE** — the predeclared rule A2 was met (worst 0.761%, n = 96, `FREEZE_RECORD_V3.md`), but volume has no domain of its own: it currently inherits the surface constants, and the evidence contains no hollow or open-cavity object. See `evidence/2026-09-13-volume-audit/VOLUME_QUALIFICATION_PROTOCOL.md` §3 |
| **SG-2b** characterisation | What does it do outside? | the full analytical grid | reports worst case, strata, fraction above 1% — no verdict | reports 18.26% worst case; 68.2% of the `rho_in < 3` stratum above 1% |

### 7.2 Why SG-2a is not simply a PASS

The number is met, and the earlier framing that no prospective volume evidence
exists was itself corrected in round 3. But a qualification item must test a
claim against *its own* declared domain, and the volume claim has never declared
one — `DOMAIN_RHO_IN_MIN` and `DOMAIN_ANISO_MAX` were derived for the area
estimator. Constructing SG-2a out of borrowed constants would create the
appearance of an independently qualified volume claim where none has been
defined. The correct order is: volume declares its domain (roadmap R-3), a
volume confirmation set is built to it, then SG-2a tests it. Until then the item
is INSUFFICIENT EVIDENCE, and reporting it as PASS would be exactly the
gate-semantics change this repository refuses to make silently.

The claim register now classifies M1 as `NOT QUALIFIED`: the met A2 rule is
informative evidence, but no volume-specific domain or dedicated qualification
record exists. SG-2a therefore remains `INSUFFICIENT EVIDENCE` in the gate
architecture and must not borrow the surface domain constants.

### 7.3 Invariants any implementation must preserve

1. The `b` items never carry PASS or FAIL. A characterisation item that can fail
   is a qualification item with a hidden criterion.
2. The `a` items must fail closed: a missing domain flag, a missing scope
   declaration, an empty in-domain subset, or a record whose evidence files do
   not match the head commit must produce FAIL, never PASS by absence. The
   adversarial tests in §4.3 cover exactly these.
3. Neither item may read a method-development directory. Development evidence
   becomes a validation record by being wired into a gate, which is the
   conversion this repository prohibits; the `a` items read confirmation
   evidence only.
4. The domain-restricted subset must be recorded as a count alongside the
   verdict. "Every case passed" over an unstated n is not a reportable result.
5. `surface_weights_origin` must be part of the `a` items' evidence contract: a
   confirmation case measured with cached or solved weights is outside the
   evidence-bearing configuration (see `INTENDED_USE_AND_ESTIMANDS.md` §2,
   note 3) and must not count toward a qualification item. All 480 confirmation
   cases resolve to packaged tables — their five spacings (`1.2x1x1`, `2.2x1x1`,
   `2.4x0.6x0.6`, `3.5x1x1`, `5x1x1`) are all in the packaged ratio set, and the
   55 packaged tables cover every stencil radius in
   `M_CANDIDATES = (1, 2, 3, 4, 5)` — so this is a guard against future drift,
   not a present exclusion.
