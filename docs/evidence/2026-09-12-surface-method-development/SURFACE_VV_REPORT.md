# Surface-area measurement: method development and V&V

## Verdict

**`FAIL / NOT QUALIFIED`.**

`crofton_minimax_adaptive_v2` was frozen with a declared applicability domain
and three per-case acceptance rules, then run once against a confirmation grid
that had never been executed. Inside the declared domain it violated all three:

| rule | criterion | result | worst case |
|---|---|---|---|
| A1 area | `< 5 %` | **FAIL**, 1/16 | `confirm-021-cylinder-size10-sp1.5x1x1-rot1`, `5.411 %` |
| A2 volume | `< 1 %` | **FAIL**, 1/16 | `confirm-020-cylinder-size10-sp1.5x1x1-rot0`, `2.273 %` |
| A3 sphericity | `< 5 %` | **FAIL**, 3/16 | `confirm-020-cylinder-size10-sp1.5x1x1-rot0`, `5.737 %` |

No qualified surface-area estimator exists for this pipeline. The production
estimator remains `marching_cubes_binary_lewiner_v1`, unchanged, still failing
the §9 criterion, with its failing evidence intact.

The candidate is nonetheless a large improvement over the legacy estimator on
identical cases — worst `5.41 %` vs `11.70 %`, 1 violation vs 14 of 16 — and
its failure is localised and understood. That is a reason to continue, not a
reason to relabel the result.

---

## 1. What was wrong with the production estimator

Binary marching cubes was diagnosed, not just measured. The error separates
into two terms with different behaviour under refinement.

* **A non-convergent orientation bias.** On a smooth surface the reconstructed
  facet normals are drawn from a finite set fixed by the lattice, so the
  measured area converges to a value above the truth. The inflation factor is
  `1.03–1.13` on isotropic voxels depending on surface orientation and does not
  decrease with resolution. Spheres converge to roughly `+8.7 %`.
* **A convergent edge/crease term.** Lattice-representable polyhedra are *not*
  measured exactly — convex right-angle edges are chamfered — but this error
  decays as `c/ρ` and vanishes under refinement.

The practical consequence is that **"use finer voxels" does not fix the
production estimator**, and the existing resolution guard cannot detect the
problem: it fired on 0 of 96 failing cases. This is the finding that justified
replacing the estimator rather than tightening acquisition.

## 2. What the replacement is, and where its numbers come from

`crofton_minimax_adaptive_v2` is a generalised Cauchy–Crofton estimator: area
is a weighted count of mask transitions along a stencil of primitive lattice
directions. Its accuracy is governed entirely by the *plane response*, how
close `Σ_k w_k |n·u_k|` comes to `1` over all surface normals `n`.

Every constant has a stated origin and none was fitted to a phantom:

* **Weights** — the minimax solution of `min_w max_n |Σ_k w_k |n·u_k| − 1|`,
  computed by a cutting-plane semi-infinite LP and verified against a
  20 011-normal set. Pure lattice geometry plus the voxel spacing.
* **Stencil radius `m`** — chosen as the `argmin` of an *a-priori* error budget
  `t*(spacing, m) + D(spacing, m)/R_eff²`, where `t*` is the LP optimum and `D`
  is a grazing-chord deficit derived in closed form by integrating the
  missed-chord probability over a sphere's impact parameter. The two terms move
  oppositely in `m`, giving an interior optimum. Inputs are the voxel spacing
  and the object's own volume-equivalent radius — both known at measurement
  time, neither derived from a phantom.

The derived grazing term was **verified** against measurement (agreement to
`0.01–0.05 %` at usable resolution), not fitted to it. That distinction is the
reason the stencil rule is a computation rather than a tuning knob.

`Bvor_m1`, the 13-direction estimator this project already tested and rejected,
was carried through the whole study as an explicitly labelled reference. It is
not the candidate and was never presented as new. It failed again (`11.35 %`
on development in-domain).

## 3. Why it failed, precisely

The failure is a **domain-declaration error**, and the confirmation set exists
to catch exactly this.

The threshold `ρ_eff ≥ 8` was justified from development evidence pooled across
shapes: worst in-domain error `3.49 %` at `ρ_eff ≥ 8` against `5.14 %` at
`ρ_eff ≥ 6`. But the shape class that carries the crease term — cylinders —
had **no development case between `ρ_eff = 8` and `10.47`**. Spheres reached
down to `8.00` and ellipsoids to `8.64`, so the pooled statistic looked
supported while the binding class did not actually span the boundary.

The confirmation grid placed three cylinders at `ρ_eff = 8.20, 8.22, 8.28` —
directly in that gap. All three produced area errors of `4.0–5.4 %` and
sphericity errors of `5.0–5.7 %`. The threshold was an extrapolation for the
one class where it mattered, and it did not hold.

This was a genuine blind spot in the freeze, not a marginal statistical
excursion. Panel **b** of the figure is that gap.

A second, independent contributor: the **a-priori budget is a predictor, not an
envelope**. It omits the crease term by construction, and it under-predicts the
measured error on `21 %` of in-domain cases — all cylinders (budget `2.16 %`
vs measured `5.41 %` on the failing case). This was identified during
development and is why the budget was frozen as the *stencil-selection* rule
only and explicitly barred from serving as the qualification gate. Had it been
used as a gate, it would have passed the failing case.

## 4. An independent finding about the existing pipeline

The confirmation run exposed a defect **unrelated to the surface estimator**:
voxel-count volume violates the project's existing §9 `<1 %` criterion.

Worst observed `10.87 %`; several cases in the `2–10 %` range. The violations
are confined almost entirely to **unrotated cylinders**, where a flat cap falls
between slices at coarse axial spacing and a half-voxel quantisation of the cap
position propagates to a percent-level volume error. Rotated versions of the
same cylinders are an order of magnitude better.

Volume is a voxel count and is therefore **identical for the legacy estimator
and the candidate**. This is a property of the production volume measurement
that the previous development grid did not expose, and it means the §9 volume
claim is not universally satisfied. It is recorded as an open item; it is not
caused by this work and is not fixed by it.

## 5. Validity layers — what is and is not established

The brief required these be kept apart. They are not interchangeable and only
the first two were examined here.

1. **Numerical implementation correctness — established.** The baseline was
   re-measured through the production `geometry()` call path and reproduced the
   frozen record for all 188 cases to `4.5e-10`. Oracle areas are analytic or
   quadrature-checked. The LP optimum is cross-checked against an independently
   computed supremum.
2. **Analytical geometry validity — FAILED, and bounded.** The candidate fails
   the frozen criterion on the confirmation set. What *is* established: the
   error mechanism is decomposed and quantified; the failure is localised to
   creased shapes near `ρ_eff ≈ 8`; smooth shapes (spheres, ellipsoids) passed
   every in-domain confirmation case with wide margin (worst `0.84 %`).
3. **Segmentation validity — NOT ADDRESSED.** Every case here is an analytical
   phantom rasterized to a perfect mask. Real masks come from Cellpose and
   carry boundary error of their own. Nothing in this study constrains how a
   surface estimator behaves on a segmentation-derived mask. A phantom result
   is an upper bound on real-data accuracy, never a prediction of it.
4. **Assay validity — NOT ADDRESSED.** No statement is made about whether
   surface area or sphericity measures anything the assay intends to capture.
5. **Biological validity — NOT ADDRESSED.** No statement about whether these
   descriptors track a biological property of organoids.
6. **Publication readiness — NO.** With no qualified estimator, surface area
   and sphericity should not be reported as validated quantities. If they
   appear at all, they must carry the `marching_cubes_binary_lewiner_v1`
   identity and the §9 failure.

Note on the two continuous-field candidates (`A1_sdf`, `A2_occ`): they were
declared in the freeze record as **conditional, never qualifiable this round**,
because the Cellpose wrapper consumes `cellprob_threshold` and discards the
continuous field, so they cannot be computed on production data at all. As
declared, they were measured and reported: both passed all 16 in-domain
confirmation cases (worst `2.18 %`). This is **not** a qualification and must
not be cited as one. It is an argument for a specific pipeline change —
retaining the continuous field — whose own validity would then have to be
established at layer 3, where a Cellpose probability field is emphatically not
a signed distance function.

## 6. Shortcuts the brief forbade — none was taken

| forbidden | what happened |
|---|---|
| loosen the `<5 %` criterion | Unchanged throughout. Two further criteria were *added* (volume, sphericity). |
| narrow intended use to manufacture a pass | The domain was declared before the confirmation run and hashed. Cylinders — the class that fails — were deliberately kept **inside** it. |
| post-hoc parameter tuning | `Blp_m2` had the best development number of any variant (`2.747 %`) and was **rejected** in the frozen selection predicate precisely because `m = 2` would have been chosen by looking at phantom errors. The adaptive rule was selected despite a worse number. |
| report a mean instead of worst-case | Every verdict is per-case worst-case. Means appear in this report only as descriptive context, never in an acceptance rule. |
| delete failing phantoms | All 188 development and all 78 confirmation cases are in the delivered CSVs, in-domain and out. |
| tune on the confirmation set | The confirmation grid was executed exactly once, after the hash was written. It was never imported by any development script. |
| pick `ρ_min` after seeing results | `ρ_min = 8` was frozen before the run — and was **wrong**, which is reported rather than revised. |
| overwrite the legacy estimator | `marching_cubes_binary_lewiner_v1` remains the default with identity, parameters and failing evidence intact. |
| assume sphericity unaffected | Re-validated jointly; it failed, and that failure is part of the verdict. |
| unexplained weighting coefficients | Every coefficient traced to a stated optimisation or derivation in `CANDIDATES.md`. |
| repackage the rejected 13-direction estimator | Carried explicitly as `Bvor_m1`, labelled as the project's prior rejected candidate. |

One thing the brief did not forbid but should be flagged as a limitation of
*this round* rather than of the method: `M_CANDIDATES` stops at `m = 4` because
the minimax solve at `m ≥ 5` did not converge in practical time. That cap is
what forced the anisotropy limit to `3.0`, which in turn put `62` of `78`
confirmation cases out of domain and left the confirmation thinner than it
should have been. This was declared in advance, in §1 and §3 of the freeze
record, along with its consequence.

## 7. What a next round should do

1. **Re-declare `ρ_min` from shape-class-resolved evidence**, not pooled
   strata, and require the development grid to span the boundary *for the
   binding class*. On the evidence here that means `ρ_min ≈ 12` for creased
   shapes — but this must be set against a **new** confirmation set, since the
   present one has now been used and is no longer independent.
2. **Add a crease term to the a-priori budget** so it becomes a genuine
   envelope rather than a predictor, or stop presenting it as a bound.
3. **Lift the `m ≤ 4` cap** via orbit-symmetry reduction of the LP (for
   spacings `(a,b,b)` the symmetry group has order 16, cutting `577` variables
   to about `75`), which would extend the qualified anisotropy range.
4. **Decide on the continuous-field route** — retaining the Cellpose
   probability field would make the strongest candidates computable, at the
   cost of a segmentation-layer validation that does not yet exist.
5. **Treat the voxel-count volume defect separately.** It affects the current
   pipeline today, independently of any surface work.

## 8. Evidence

`freeze_record.json` carries SHA-256 for all 11 frozen files; the confirmation
grid was executed after that record was written.

| file | contents |
|---|---|
| `FREEZE_RECORD.md`, `freeze_record.json` | frozen estimator, domain, acceptance rules, selection predicate, hashes |
| `vv_confirm.csv`, `vv_confirm_summary.json` | all 78 confirmation cases and the verdict |
| `vv_dev.csv`, `vv_dev_summary.json` | all 188 development cases under the frozen code |
| `candidates_dev.csv` | development exploration, 11 estimator variants |
| `conditional_confirm.csv` | the two conditional candidates on in-domain confirmation cases |
| `MECHANISM.md` + 4 CSVs | the error decomposition of §1 |
| `CANDIDATES.md` | derivations and provenance of every constant |
| `E0_legacy_baseline.json` | legacy estimator identity and failing evidence, preserved |
| `surface_vv_patch.tar.gz` | modules, evidence directory, harness scripts, `features.py` diff |
