# Surface-area measurement V&V plan — frozen 2026-09-11

Status: **FROZEN before any candidate estimator was run on a phantom.** Written
at clean `ed81d74ac9d8c3b7b2c5efe6e9e9549d00bc7856` (equal to `origin/main`).
The only phantom result known when this was written is the historical baseline
already in the repository: the production estimator on a radius-18 µm sphere at
ZYX spacing (2,1,1) µm gives +11.73% area error
(`docs/evidence/2026-09-11/analytical-geometry.json`, reproduced in this round).
The only other computation made before freezing is the analytical
plane-orientation response of the Crofton weights (§3, E3). That is a property of
the estimator on planes; no phantom was involved.

Any change to this plan after the development grid has been run is a
recalibration. It must be recorded as a dated amendment below the freeze line and
cannot support a PASS claim without a new, unseen confirmation set.

## 1. Estimand

The estimand is the **area (µm²) of the boundary of the continuous physical
object** whose Gauss digitization (inside test at voxel centres, on the given ZYX
physical lattice) produced the binary mask. Only the binary mask and the voxel
spacing are available to the estimator. A binary mask does not determine the
continuous boundary uniquely. Every estimator in this plan therefore estimates a
quantity that is not identifiable from its input, and its error cannot be driven
to zero at a fixed resolution.

Out of scope: segmentation error, partial-volume and grey-level information
(real masks are not Gauss digitizations of a known solid), and biological
surface definitions such as membrane folding below optical resolution.

## 2. Acceptance criterion (unchanged from SCIENTIFIC_SPEC §9)

- **Scientific acceptance criterion:** absolute relative area error
  `|A_est/A_true − 1| < 0.05` for analytical shapes. This round does not change
  it, re-scope it or reinterpret it.
- **How it is evaluated:** per phantom case. A stratum PASSES only if its
  **maximum** absolute error is below 0.05. An estimator PASSES overall only if
  every case in the predeclared grid passes. Mean and median values are reported
  but can never produce a PASS.
- **Regression tolerance (distinct):** the 15% bound in
  `tests/quantification/test_geometry.py` is a regression tolerance. It detects
  unintended changes to the production estimator. It is not evidence for the
  scientific criterion and is not modified here.
- The same phantoms also report the §9 volume criterion (`<1%`) for the voxel
  count volume. Volume does not depend on the surface estimator.
- If no estimator passes every case, the result is **FAIL**. A sub-domain in
  which an estimator stays below 5% may be described, but only as descriptive
  information. Restricting the intended use to that sub-domain is a change in
  intended use, which the project owner must decide. It is not a PASS.

## 3. Candidate estimators (defined before results)

| ID | Definition | Rationale | Parameters and provenance |
|---|---|---|---|
| E0 `marching_cubes_binary_v1` (production) | `quantification.features.geometry`: binary mask padded by 1 voxel, scikit-image Lewiner marching cubes at level 0.5 with physical spacing, `mesh_surface_area` | Current method, Lorensen & Cline 1987. On binary input every vertex lies at an edge midpoint, so oblique or curved surfaces become a terrace (staircase). The terrace area does not converge to the true area as resolution increases (Windreich, Kiryati & Lohmann 2003, *Pattern Recognit.* 36:2531; Lindblad 2005, *Image Vis. Comput.* 23:111). | level 0.5 is analytical (midpoint between 0 and 1). No tunable parameter. |
| E1 `marching_cubes_gaussian_phys` | Float indicator smoothed by a Gaussian with the **same physical σ on every axis** (σ_vox = σ_µm / spacing), zero-padded by ⌈4σ_vox⌉+2, then marching cubes at level 0.5 | A symmetric kernel places the 0.5 level of a smoothed planar step exactly at the step, which removes terracing for planes. On curved surfaces the level set moves inward by about σ²·H, where H is mean curvature. This causes a predictable negative area bias for small or highly curved objects and rounds sharp edges. | **σ = 1.0 × max(spacing)**. This is a fixed engineering choice, basis HEURISTIC: terraces have the period of the coarsest sampling interval, so smoothing below that scale cannot remove them. It is **not calibrated** against any phantom. Sensitivity variants σ = 0.5× and 1.5× max(spacing) (E1a/E1b) are reported for sensitivity only and **may not be selected**. |
| E2 `marching_cubes_signed_distance` | Signed Euclidean distance field with physical sampling, `EDT(outside) − EDT(inside)` (distances between voxel centres of opposite class, SciPy exact EDT), padded by 2 voxels, marching cubes at level 0 | The field is 1-Lipschitz. Linear interpolation therefore places vertices at orientation-dependent positions instead of snapping every vertex to an edge midpoint. For axis-aligned faces the zero lies at the in/out midpoint, as in E0. No closed-form bias result is known for curved surfaces (basis CONTEXT_DEPENDENT). | Level 0 is analytical. No tunable parameter. |
| E3 `crofton_13dir_voronoi` | Discrete Cauchy–Crofton estimate `A = 2 Σ_i w_i N_i v / |d_i|`, where `N_i` counts in/out transitions along the 13 primitive lattice directions `d_i` in physical space, `v` is voxel volume, and `w_i` is the spherical-Voronoi solid-angle fraction of ±d̂_i | Crofton/Cauchy integral geometry. The discrete lattice version is used by MorphoLibJ (Legland, Arganda-Carreras & Andrey 2016, *Bioinformatics* 32:3532), and solid-angle direction weights follow Ohser & Schladitz (*3D Images of Materials Structures*, 2009). The analytical plane response computed before freezing is mean 1.000 over orientations, with range 0.929–1.023 at (1,1,1), 0.878–1.094 at (2,1,1), 0.848–1.102 at (3,1,1) and 0.855–1.101 at (2,0.7,0.7). Planar faces may therefore fail at some orientations even with perfect crossing counts. | Direction set (13) and Voronoi weights are analytical in the lattice geometry. No tunable parameter. |

Rejected before any results were seen:
- MC area × constant correction factor. Fitting the factor on these phantoms
  would leak the evaluation set, and the factor would still depend on
  orientation.
- Lindblad weighted local configurations. The published weights are optimized
  for isotropic lattices only; anisotropic spacing would need a new calibration.
- Laplacian/Taubin mesh smoothing. Its iteration and λ/μ parameters have no
  analytical provenance.
- Grey-level/probability-map isosurfaces. The feature pipeline receives only
  label masks. This is a documented future option, because binarization
  destroys the sub-voxel information.

## 4. Development grid (used for comparison and selection)

- Digitization: an object centre at the array centre plus a seeded sub-voxel
  offset, uniform in [−0.5, 0.5) voxels per axis (`numpy.random.default_rng(20260911)`).
  Inside test at voxel centres. Background margin of at least 3 voxels.
- ZYX spacings (µm): (1,1,1), (2,1,1), (3,1,1), (2,0.7,0.7).
- Sphere: R ∈ {3, 6, 12, 24, 36} µm with 3 offsets each.
- Ellipsoid: semi-axes s·(1, 0.75, 0.5) with s ∈ {6, 12, 24, 48} µm.
- Finite cylinder with flat caps: radius r ∈ {3, 6, 12, 24} µm, height 3r.
- Orientations (ellipsoid and cylinder): identity plus 3 uniform random rotations
  (`scipy.spatial.transform.Rotation.random(3, random_state=20260911)`). The
  same rotations are used for every size and spacing (paired design).
- Total: 60 sphere + 64 ellipsoid + 64 cylinder = 188 cases per estimator.

Oracles are closed form: sphere `4πR²`; cylinder `2πrh + 2πr²`; ellipsoid via
Legendre incomplete elliptic integrals (`scipy.special.ellipkinc/ellipeinc`).
The ellipsoid formula is checked against independent numerical quadrature of
the parametric area element, which must agree to a relative error below 1e-8
before any estimator result is used.

## 5. Required strata and metrics

Resolution ratio ρ = (smallest semi-axis, or cylinder radius) / max(spacing).
The ρ strata are `very_small` (ρ<3), `small` (3≤ρ<6), `medium` (6≤ρ<12) and
`large` (ρ≥12). Anisotropy is max/min spacing: 1, 2, 2.857 or 3.

For each estimator and each stratum (overall; shape; ρ; spacing/anisotropy;
shape×ρ) the report gives N, signed mean bias, median signed error, mean and
median absolute error, **maximum absolute error with the worst case identified**,
and PASS/FAIL against 0.05. Derived sphericity error and voxel-volume error are
reported for the same cases.

These are deterministic phantoms. N is the number of exercised configurations,
not a population sample, so confidence intervals do not apply. Spread over
offsets and orientations is reported through min/max.

## 6. Decision rule for production

1. Only E0, E1 (σ = 1.0×max spacing), E2 and E3 are eligible. E1a and E1b are
   not eligible.
2. A candidate may replace the production default only if it **passes every
   development case** and then **passes every case of a separate confirmation
   grid** that is run once, after the selection record is written.
   Confirmation grid: spacings (1.5,1,1), (2.5,0.8,0.8) and (4,1,1); sphere R
   {5,10,20,30} with 2 offsets; ellipsoid s·(1,0.6,0.45) with s {10,20,40};
   cylinder r {5,10,20} with h = 2.5r; 3 orientations (identity + 2 random
   rotations) and offsets from seed 424242; 78 cases.
3. If no candidate meets rule 2, **the production default and the historical
   definition of `surface_area_um2` stay unchanged**, and the surface criterion
   is reported as FAIL. A candidate that is better in every ρ stratum can be
   recommended to the project owner as a versioned option, but it is not made
   the default in this round.
4. Whatever the outcome, any production surface value that is exported must be
   traceable to its estimator name, version, parameters, spacing and package
   versions, so that a later method change cannot silently reinterpret
   historical results.

## 7. Evidence produced

- `surface_vv.py`: phantoms, oracles, candidate estimators, summaries.
- `surface_vv_dev.csv`, `surface_vv_dev_summary.json`, `surface_vv_dev_summary.md`.
- `SURFACE_SELECTION.md`: written after the dev grid and before confirmation.
- `surface_vv_confirm.*`: only if rule 2 is triggered.

---
Amendments after freeze: none.
