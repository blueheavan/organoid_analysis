# Scientific review and gate update — 2026-09-26

> Later the same day SG-5 was also qualified (gate 4/9); see
> [VALIDATION_UPDATE_2026-09-26-sg5.md](VALIDATION_UPDATE_2026-09-26-sg5.md).

## Outcome

| Item | Before | After | Basis |
|---|---|---|---|
| SG-1a surface qualification | FAIL (stale record) | **PASS** | canonical record `2026-09-26-analytical-qualification-record`; 966 in-domain cases, worst \|area error\| 1.41 % (< 5 %) |
| SG-2a volume qualification | FAIL (stale record) | **PASS** | same record; 325 cases in the new volume domain, worst \|volume error\| 0.24 % (< 1 %) |
| SG-4A four-state analytical rule | PARTIAL | **PASS** | canonical record `2026-09-26-viability-qualification-record`; 991 known-rule cases, 3,987 comparisons, 0 mismatches |
| SG-1b / SG-2b characterization | rejected (stale) | reported | canonical record `2026-09-26-analytical-geometry-record`; values identical to 2026-09-16 |
| SG-3A, SG-3B, SG-4B, SG-5, SG-6A, SG-6B | not passed | **not passed** | need specimens, blinded annotation, orthogonal assays, a study design or instrument standards |
| **Gate** | NOT PASSED 0/9 | **NOT PASSED 3/9** | |

The three PASS items are **analytical** claims about software on phantoms and
known-rule tables. They do not validate segmentation, biology, physical units
or statistical inference. The release verdict stays
`NOT READY FOR THE SPECIFIED RESEARCH USE`.

## 1. Volume (SG-2a): the D-2 domain was too wide, the estimator was not wrong

**Finding.** The 2026-09-26 counterexample (two lattice-centred spheres with
identical masks, `rho_in` ≈ 10.6) was not an isolated curiosity. Every earlier
phantom grid drew sub-voxel offsets at random and so never sampled
lattice-symmetric placements. With them included, voxel-count error inside the
D-2 gate (`rho_in ≥ 10`, anisotropy ≤ 4) reaches 2.7 % on an isotropic grid
(`evidence/2026-09-26-volume-domain-development/dense_sphere_scan.csv`).
Random placements stayed ≤ 0.47 %, which is why the 2026-09-12 confirmation
(96 random-offset cases, worst 0.76 %) did not see it.

**Why the estimator was kept.** Under a uniformly random sub-voxel translation
the voxel count × voxel volume is a design-unbiased volume estimator for any
bounded shape — the Cavalieri point-counting principle (Gundersen & Jensen,
*J. Microsc.* 147:229–263, 1987,
[doi:10.1111/j.1365-2818.1987.tb02837.x](https://doi.org/10.1111/j.1365-2818.1987.tb02837.x),
verified via PubMed). Its error is therefore dispersion, not bias, and it
shrinks with resolution. A smoothing or surface-reconstruction volume would add
a shape-dependent bias, and the counterexample proves no mask-only estimator can
remove the worst-case ambiguity. The mathematical argument for the lattice gap
(Legendre's three-square theorem: only `4^a(8b+7)` is excluded, so consecutive
representable squared radii differ by ≤ 3) is checked numerically in
`tests/quantification/test_volume_identifiability.py`; it is an elementary
result, not a cited one.

**Domain change (narrowing only, decision D-15).** New variable
`rho_vol = r_in / cbrt(sz·sy·sx)` — inscribed radius in geometric-mean voxel
edges. It unifies lattice error across anisotropies where `rho_in` (coarsest
axis) does not. Rule declared in `dev_volume.py` before execution: the smallest
ladder value at which every development case, including five lattice-symmetric
placements per cell, is ≤ 0.5 % (2× margin). 2,304 development cases plus a
36,030-row dense sphere scan selected **36**; the binding class is the
axis-aligned lattice-centred torus (0.59 % at 32). Because 36 was the ladder
top, a 384-case extension to 48 was run before freezing; worst 0.45 %.

Practical consequence: at 0.65 × 0.65 × 3 µm voxels (geometric mean 1.08 µm)
volume is qualified for smooth objects with inscribed radius ≥ 39 µm. Nuclei
and small organoids are outside the volume claim; they still receive a
measurement and `volume_in_qualified_domain = False`.

**Confirmation.** `confirm3`, frozen in
`evidence/2026-09-26-analytical-qualification-protocol/` and committed
(`8d59602`) before execution: new aspects, new absolute spacings, `rho_vol`
11–42, five lattice-symmetric + two random placements per cell, 1,176 cases.
All three pre-registered predictions held: SG-1a worst < 2 % (1.41 %); SG-2a
worst < 0.6 % (0.24 %); lattice-symmetric cases between the surface and volume
gates exceed 1 % (4 cases, up to 1.35 %). Deviation disclosed in the protocol
§8: three cases were executed as a harness smoke test before the freeze; no
parameter changed afterwards.

## 2. Surface (SG-1a): D-14 contract implemented

The domain-restricted record authorised by D-14 now exists. It re-executes the
2026-09-12 `confirm2` grid through the production estimator (areas reproduce
the development implementation to 4.9 × 10⁻¹⁰) and adds `confirm3`, whose
lattice-symmetric placements test the known a-priori-budget exception. Creased
objects remain NOT QUALIFIED (reported: 13.8 % worst on `confirm2`).

## 3. Viability rule (SG-4A): one defect fixed

An oracle written from `INTENDED_USE_AND_ESTIMANDS.md` §5.1–5.3 (not from the
code) found that `classify()` appended `;outside_control_range` to
`nonfinite_scaled_signal` whenever either scaled value was ±∞ or the other was
finite and out of range — reporting an unusable value as a measurable
out-of-range signal. The state was already correct (`indeterminate`); only the
reason string changed. Summaries now also export
`viability_fraction_definition = "classifiable-denominator/2"` so fractions
under the superseded all-eligible denominator cannot be mixed silently.

Two boundary conventions the specification left implicit are now written in
`viability_rule_evidence.SPEC_CONVENTIONS` (strict `<` for "both markers low";
range flag only for finite values). Independence is limited (tier D): the
oracle's author had read the implementation.

## 4. Record contracts and gate

* `validation/record_contract.py` — shared seal, provenance, scope, protocol
  freeze, environment and item-status checks.
* `analytical-qualification/1` (SG-1a, SG-2a) and `viability-analytical/1`
  (SG-4A): results keyed by gate item; re-derivation from raw CSV on every
  verification; re-execution (full for SG-4A; worst adjudicated cases plus every
  120th case for SG-1a/SG-2a).
* The gate grants PASS only from a verified **canonical** record that reports
  PASS for that item; scope-clean or VOID records give INSUFFICIENT EVIDENCE.
* New adversarial tests (`tests/validation/test_qualification_evidence.py`):
  raw edits, hand-edited results, forged-and-rehashed raw values (caught by
  re-execution), protocol/generator/domain-constant edits, shifted gates,
  changed denominators, incomplete confirmation sets (VOID), scope-clean PASS.

Workflow: `pixi run qualification-record --contract analytical|viability`.

## 5. Verification in this session

| Check | Result |
|---|---|
| Environment | Python 3.12.14, numpy 2.5.2, scipy 1.17.1, scikit-image 0.26.0, pandas 2.3.3 (the `pixi.lock` versions) on linux-x86_64; records store the platform |
| `ruff check src tests scripts` | PASS |
| Typecheck ratchet | identical to the pre-change tree (differences only from cellpose/torch absent in this container) |
| Gate + evidence tests | 159 passed |
| Full suite | 591 selected tests run (render tests excluded by default): 10 failed, all `ModuleNotFoundError: cellpose` in `tests/segmentation/*` (cellpose/torch are not installable in this container; segmentation code is untouched); no other failure or error. One test was deselected: `test_mask_feature_sparse_ids_are_compacted_before_regionprops` needs an uncached Crofton LP solve for spacing ratio (1, 0.75, 0.5) that exceeds its 600 s budget on this CPU — the same behaviour recorded on 2026-09-25, on a path this change does not touch. Both need a run in the Pixi environment on Apple Silicon |

## 6. What remains, and why software cannot close it

| Item | Required evidence (unchanged, D-5…D-13) |
|---|---|
| SG-3A / SG-3B | blinded, modality-matched independent 3D annotation; ≥ 47 biological samples per stratum |
| SG-4B | blinded two-annotator object classification (D-8), graded-insult series, held-out batches |
| SG-5 | frozen design with an independent experimental unit; verified CR2 (clubSandwich golden values) and the D-10 coverage/type-I simulations |
| SG-6A / SG-6B | ratio and absolute calibration standards imaged on the intended instruments |

SG-5 is the only one of these that is mostly computational; it still needs the
owner's actual replicate design before a coverage claim can refer to it.
