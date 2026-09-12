# Round 2 protocol — surface-area method development and V&V

Written before any round-2 grid exists. Round 1 ended `FAIL / NOT QUALIFIED`;
this document states what is inherited unchanged, what changes, and why each
change is a response to the round-1 failure rather than an attempt to rescue
its verdict.

## 0. Status of round-1 material

| material | round-2 status |
|---|---|
| round-1 development grid (188 cases) | **development data** |
| round-1 confirmation grid (78 cases) | **SPENT.** Now development data. May never again serve as confirmation. |
| round-1 verdict `FAIL` | **permanent.** Not revised, not superseded, not deleted. |
| `crofton_minimax_adaptive_v2` frozen record + hashes | preserved; round 2 produces a *new* version, it does not edit the old one |

Everything below that is justified by round-1 numbers is justified by
**development** evidence. No round-2 acceptance decision may cite the round-1
confirmation grid.

## 1. Inherited unchanged

These carry over verbatim. They were not implicated in the failure and
re-deriving them would only add risk.

* **Estimator family** — generalised Cauchy–Crofton: area as a weighted count
  of mask transitions along primitive lattice directions.
* **Weight derivation** — minimax plane response, `min_w max_n |Σ_k w_k|n·u_k| − 1|`,
  solved by cutting-plane semi-infinite LP, verified against a 20 011-normal set.
* **Grazing-deficit derivation** — the closed form `D(spacing, m)`, obtained by
  integrating missed-chord probability over impact parameter; verified against
  measurement to 0.01–0.05% in round 1.
* **Acceptance criteria** — area `< 5%`, volume `< 1%`, sphericity `< 5%`,
  **per case, worst case, in-domain**. Unchanged in value and in form. No mean,
  median or quantile may be substituted.
* **One-shot confirmation discipline** — the confirmation grid is specified and
  hashed before it is generated, executed exactly once, and its verdict stands.

## 2. Change 1 — the domain variable

**This is the substantive scientific change and the direct fix for the failure.**

Round 1 declared the domain on `rho_eff = R_eff / max(spacing)`, where `R_eff`
is the volume-equivalent radius from the voxel count. The three cases that
failed had `rho_eff = 8.20, 8.22, 8.28` — comfortably inside a threshold of 8 —
but their actual cylinder radius was only `6.57–6.63` voxels. `R_eff` is a
*global* size measure; the error mechanism is driven by *local* feature size,
because both the crease term and the grazing term scale with the local radius
of curvature. For a compact blob the two coincide; for an elongated or creased
object `R_eff` overstates how well resolved the surface actually is.

Round 2 therefore declares the domain on

```
rho_in = r_inscribed / max(spacing),
r_inscribed = max( distance_transform_edt(mask, sampling=spacing) )
```

the largest inscribed sphere radius, computed from the mask and the voxel
spacing alone. This is **computable at measurement time** — it requires no
knowledge of shape class, which is the property the production code cannot
have. A shape-class-resolved threshold, the naive reading of the round-1
failure, is *not* implementable in production and is explicitly rejected here.

Evidence from round-1 data (both grids, anisotropy ≤ 3), `out/rho_in_scoping.csv`:

| predicate | n retained | worst \|area err\| | violations |
|---|---|---|---|
| `rho_eff ≥ 8` (round 1) | 94 | 5.411% | **1** |
| `rho_in ≥ 8` | 70 | 3.990% | **0** |
| `rho_eff ≥ 10` | 74 | 3.990% | 0 |
| `rho_in ≥ 10` | 66 | 3.990% | 0 |

Both `rho_in ≥ 8` and `rho_eff ≥ 10` clear round-1 data. **They are not
equivalent claims.** `rho_eff ≥ 10` works here only because sphere, ellipsoid
and cylinder happen to share a narrow range of `R_eff / r_inscribed`; that
coincidence is exactly the kind that produced the round-1 failure. `rho_in` is
preferred on mechanistic grounds, and §4 adds shape classes specifically
designed to discriminate the two. The discriminating test is declared **now**,
before those shapes are run.

The a-priori stencil-selection budget also switches from `R_eff` to
`r_inscribed` as its curvature scale, for the same reason and in the same
direction (more conservative).

## 3. Change 2 — stencil radius cap

Round 1 capped `m ≤ 4` because the minimax LP did not converge at `m = 5` in
practical time (577 free variables). That cap forced the anisotropy limit to
3.0, which put 62 of 78 confirmation cases out of domain.

Round 2 constrains the weights to be constant on the orbits of the symmetry
group of the `(a, b, b)` voxel metric — axis permutation of the two equal axes
and sign changes — which is legitimate because the minimax objective is
invariant under that group and convex, so a symmetric optimum exists. This
reduces the variable count by roughly the group order and is expected to make
`m = 5, 6` tractable.

**Declared in advance:** the reduced solve must reproduce the round-1 `m ≤ 4`
weights to solver tolerance before any new radius is used. Any radius that does
not solve within a stated wall-clock budget is excluded, and the exclusion is
recorded before a single phantom is run — the same discipline round 1 applied
to `m = 5`, so that compute cost never silently shapes the science.

## 4. Change 3 — development coverage and new shape classes

The round-1 failure was possible because cylinders — the binding class — had no
in-domain development case between `rho_eff = 8` and `10.47`. The threshold was
an extrapolation for the one class that mattered.

Round 2 requires:

* **Boundary-spanning coverage.** Every shape class must have development cases
  densely spanning the candidate threshold region, so no threshold is ever set
  by extrapolation for any class.
* **Three new shape classes**, chosen to attack specific assumptions:
  * **box** — maximally creased; every surface point is on a face, edge or
    corner. Round 1 showed an irreducible ~4% crease term on cylinder rims, so
    a box is a genuine risk of outright failure. It is included *because* it
    may fail.
  * **capsule** — hemispherically capped cylinder; elongated but entirely
    smooth. Separates "elongated" from "creased", which cylinders confound.
  * **torus** — non-trivial topology, and `R_eff / r_inscribed` far from the
    blob value. This is the case that discriminates `rho_in` from `rho_eff`.
* **Disjoint randomisation.** New seeds for rotations and offsets, disjoint
  from both round-1 grids.

Analytic area and volume for all three are closed-form and will be
cross-checked against quadrature to `1e-6` before any estimator is run.

## 5. Change 4 — volume domain

Round 1 exposed an independent defect: voxel-count volume violates the existing
`< 1%` criterion on flat-faced objects at coarse axial spacing (worst 10.87%).
Unlike marching-cubes area, voxel-count volume **converges**, so a resolution
domain is a legitimate remedy. Round-1 data gives 0 violations in 79 cases at
`rho_eff ≥ 10`. Round 2 declares a volume domain in the same frozen record, on
the same `rho_in` variable, and subjects it to the same one-shot confirmation.

## 6. What is *not* changing

Stated explicitly so that their absence from the change list is deliberate
rather than an oversight:

* The 5% / 1% / 5% criteria and their per-case worst-case form.
* The legacy estimator's identity, parameters and failing evidence.
* The rule that the a-priori budget may select a stencil but may **never**
  serve as a qualification gate — round 1 established it under-predicts on
  21% of in-domain cases.
* The continuous-field candidates remain **conditional and unqualifiable**
  while the segmentation wrapper discards the field. Round 1 additionally
  refuted the cheap substitute: a field reconstructed from the binary mask by
  distance transform fails 72/78 development cases, because it inherits the
  staircase it was meant to remove.
* The six validity layers stay separated in reporting. Round 2 addresses layers
  1–2 only. Segmentation, assay and biological validity remain unestablished
  whatever the verdict.

## 7. Predeclared failure modes

Round 2 is expected to fail if any of these hold; recording them now prevents
them being rationalised later.

1. The box class exceeds 5% at every resolution — the crease term is
   irreducible for the estimator family. Consequence: boxes leave the qualified
   domain, and the domain must then be stated by a *computable* creasedness
   predicate, not by "not a box".
2. The symmetry-reduced LP does not converge at `m ≥ 5`. Consequence: the
   anisotropy limit stays near 3.0 and the qualified domain stays narrow.
3. The torus shows that neither `rho_in` nor `rho_eff` controls the error for
   non-blob topology. Consequence: the domain needs a topology or
   thickness-uniformity guard.
4. `rho_in` clears development but fails confirmation. Consequence: a second
   `FAIL`, reported as such, and the round-1 conclusion that this estimator
   family cannot be qualified on creased shapes at attainable resolution
   becomes the finding.

A second `FAIL / NOT QUALIFIED` remains an acceptable outcome of this round.
