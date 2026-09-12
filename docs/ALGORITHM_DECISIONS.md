# Algorithm evidence audit — 2026-09-10

This record describes the current implementation; it does not retrospectively justify its selection. No scientific method was replaced or tuned in this audit. `ESTABLISHED` describes a mathematical/method basis, not suitability for these biological images. Project-data suitability remains `INSUFFICIENT EVIDENCE` unless explicitly bounded below. Installed versions and executable defaults are captured in [installed-methods.txt](evidence/2026-09-10/installed-methods.txt).

## D1. Parent assignment by maximum voxel overlap

`quantification/multilevel_relationships/hierarchy.py` counts label co-occurrences, takes the largest count, and resolves ties to the smallest parent ID. This is an exact operational containment rule, not a calibrated biological assignment model. Origin `HEURISTIC`; basis `CONTEXT_DEPENDENT`; evidence `PARTIAL` (hand-counted assignments and large-ID tests pass). Nearest-centroid assignment is a meaningful alternative for different data, but was not benchmarked. Misregistration, fragments, ties and cross-parent objects require review; same shape/spacing cannot establish registration. No qualified real hierarchy annotations are available.

## D2. Nucleus-organoid membership inherited through the cell

`workflows/multilevel_measurement_workflow.py` inherits the organoid through the cell and retains direct nucleus–organoid overlap as an audit field. Origin/basis `HEURISTIC`; evidence `PARTIAL`. Tests verify transitivity and mismatch flags. Direct overlap is an alternative definition, not an empirically rejected inferior method. No biological transitivity/assignment performance claim is supported.

## D3a. Adopted surface estimator: minimax-weighted lattice transition counts (2026-09-12)

`quantification/surface_crofton.py::measure` is the production surface-area
estimator, `crofton_minimax_sym_v3`. Area is a weighted count of mask
transitions along the lattice directions of a radius-`m` stencil (discrete
Cauchy/Crofton). The weights solve a minimax linear program over plane
orientations; its optimum `t*` bounds the estimator's response to a planar
surface element before any measurement is taken. The search is restricted to
weights constant on the orbits of the signed axis permutations preserving
`diag(spacing)`; that restriction is **exact** — the objective is invariant
under the group and convex, so an optimum exists on the invariant subspace —
which is what made radius 5 solvable (577 free variables reduced to 40 orbits
at isotropic spacing). The stencil radius is chosen per object to minimise
`t* + D/rho_in**2`, a two-term a-priori budget: the first term is a proven
worst-case bound for a planar element, the second an asymptotic model of
grazing curved elements, so their sum is a validated predictive budget rather
than a proven bound (it is exceeded by ≈1.3x for exactly lattice-centred
spheres; see VALIDATION_UPDATE_2026-09-12.md).

Origin `DERIVED`; basis `ESTABLISHED` (integral geometry; the weights carry a
derivation, not a fit); implementation evidence `CONFIRMED` on an untouched
confirmation set. Suitability for the specification's <5 % area criterion is
`SUPPORTED` **within the declared domain only**: ρ_in ≥ 10, anisotropy ≤ 4,
smooth closed surfaces — worst 0.951 % on 96 untouched in-domain cases, 0.790 %
on the 40 in-domain cases of the frozen V&V grid. Outside that domain the
claim is `NOT SUPPORTED`: creased surfaces are excluded at any resolution by
the plane-response argument, and the unrestricted grid still FAILS at 25.81 %.
Objects outside the domain are measured and flagged, never refused or silently
passed.

The minimax optimum is not unique, so the method is the specific frozen weight
vectors and they ship with the package; an on-demand solve is conforming but
not evidence-bearing and is reported through `weights_origin`. Changing the
domain constants, the direction set, the stencil candidates or the weight
tables voids the qualification and requires a new confirmation run on an unused
set — editing the version-lock test's number instead is explicitly not the
remedy.

## D3. Native-spacing marching-cubes surface (SUPERSEDED 2026-09-12)

Superseded by D3a as the production estimator; retained unchanged, reachable
through `features.legacy_surface_area()` under
`marching_cubes_binary_lewiner_v1`, so historical `surface_area_um2` values
stay reproducible, together with the evidence that it fails the §9 criterion.
It defines no exported column. The description below is the historical record
of the decision as it stood on 2026-09-11.


`quantification/features.py::surface_mesh` uses scikit-image 0.26.0's **Lewiner** implementation (the installed default), level 0.5, step 1, padded zero exterior and `allow_degenerate=False`. Lorensen–Cline is historical background, not the exact implementation. [Official method documentation](https://scikit-image.org/docs/stable/api/skimage.measure.html#skimage.measure.marching_cubes) identifies the algorithm and axis-ordered spacing. Origin `PUBLISHED_METHOD`; basis `ESTABLISHED`; implementation evidence `PARTIAL`; suitability for the specification's <5% surface error is `NOT SUPPORTED` / `FAIL`: the radius-18 µm sphere at (2,1,1) µm has 11.7347% area error. Sphericity is the dimensionless isoperimetric expression; it is not clipped. Changing smoothing, surface estimator or tolerance requires a new scientific decision and validation, so none was changed.

Classical organoid geometry fills enclosed holes; multilevel geometry uses raw label voxels. Since 2026-09-11 Web feature schema 2.0 uses raw labels by default and offers an explicit envelope option: centroid and solidity now use the same support as volume/surface/axes. The earlier Web mixed-support behavior is corrected, with full numeric precision. Classical `cells` still reports raw cell volume plus envelope geometry. Border padding closes truncated objects computationally, not biologically. Web QC retains all rows and does not apply the classical eligibility filter. Alternative surface estimators and acquisition-dependent bias correction are NOT ASSESSED.

## D4. Face-contact topology

`topology.py` sums exposed shared voxel faces: Z-normal area = sy*sx, Y-normal = sz*sx, X-normal = sz*sy. Origin `HEURISTIC` (chosen contact definition); basis `ESTABLISHED` for discrete arithmetic, `CONTEXT_DEPENDENT` for biological contact; evidence `PARTIAL`. Controlled tests verify counts and units. Label-face adjacency is a geometric proxy; it does not prove membrane contact. kNN and Delaunay measure different notions of neighborhood. Unassigned/inter-organoid contacts remain explicitly flagged.

## D5. MAD volume QC

`qc.py` uses 0.67448975*(V−median)/MAD with cutoff 3.5. Gaussian consistency follows directly from Phi^-1(0.75); it does not justify the cutoff, population pooling, or the MAD=0 rule (`~isclose` to median). Origin `HEURISTIC`, basis `HEURISTIC`, evidence `PARTIAL` for implementation and `INSUFFICIENT EVIDENCE` for biological abnormality. QC flags do not remove objects in multilevel aggregates. `qc_status=pass` means no implemented flag, not validated measurement accuracy.

## D6. Equivalent-radius position and EDT depth

`spatial.py` normalizes centroid radius and nearest-voxel EDT by (3V/4π)^(1/3). EDT is distance to the nearest **background voxel center**, not the continuous isosurface; [SciPy's definition](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.distance_transform_edt.html) is the numerical reference. A rounded centroid outside its assigned parent gives zero. Boundary foreground voxels can have positive depth. Origin `HEURISTIC`; basis `CONTEXT_DEPENDENT`; evidence `PARTIAL`. This differs from `cellular_measurements.py`'s exact centroid-to-exposed-voxel-face distance. Core/periphery cutoffs are uncalibrated. Full-volume EDT storage per parent has unbenchmarked memory scaling.

## D6b. Moment-equivalent principal axes

`features.py` adds diag(spacing²/12) to the covariance of voxel centers. For a uniform interval of width s, Var=s²/12; for a uniform solid ellipsoid with semi-axis a, the matching covariance eigenvalue is a²/5. Thus the full diameter is 2*sqrt(5λ). This corrects the earlier catalog's confused covariance/inertia explanation. Basis `ESTABLISHED` from analytical integration; evidence `PASS` for controlled voxel/moment arithmetic, `INSUFFICIENT EVIDENCE` for biological shape interpretation. It is not a fitted biological boundary ellipsoid.

## D7. Marker background and viability signal states

`features.py` measures the raw filled ROI against a physical local background shell excluding other labels. It preserves negative corrected means. `phenotyping/viability.py` aggregates controls by wells then biological replicate, scales by batch endpoints, and applies fixed 0.30/0.60 rules. Dtype maximum is only a fallback saturation limit; stored uint16 is not proof of a 16-bit detector. The 1e-9 noise floor is an uncalibrated numerical heuristic.

[Manufacturer documentation](https://www.thermofisher.com/order/catalog/product/P21493/faqs) supports the distinct esterase/retention and membrane-integrity probe mechanisms. It does **not** support this organoid-level scaling, gates, minimum control N, or a live-cell fraction. Origin `HEURISTIC` for the pipeline rule and `CONTROL_DERIVED` for fitted endpoints; basis `HEURISTIC`; evidence `INSUFFICIENT EVIDENCE`. Synthetic controls share the generator's assumptions and are not an orthogonal assay. Viable-like/mixed/compromised-like states are signal patterns. An independent assay and dye penetration/bleedthrough/background/saturation qualification are missing.

## D8. Global Otsu and physical-distance watershed

`watershed_instances.py`: optional Gaussian background subtraction → polarity inversion → Gaussian smoothing → Otsu (or numeric threshold) → physical closing → optional hole filling → small-component removal → physical EDT h-maxima/suppressed seeds → 6-connected watershed → small-instance removal. All thresholds remain unchanged. [Otsu API/method reference](https://scikit-image.org/docs/stable/api/skimage.filters.html#skimage.filters.threshold_otsu) and [watershed API](https://scikit-image.org/docs/stable/api/skimage.segmentation.html#skimage.segmentation.watershed) establish the operations, not their suitability or physical defaults. Basis `ESTABLISHED` for the named methods, `HEURISTIC` for their composition; evidence `PARTIAL`. Probability input must be registered finite [0,1]; imported labels must be connected nonnegative integers. Baselines: unsplit threshold components and imported qualified masks. No superiority comparison was performed. Brightfield halos, nonuniform stain, touching/lumen-containing structures, saturated fluorescence and insufficient Z sampling remain consequential limitations.

## D9. One-to-one instance evaluation

`validation/segmentation_metrics.py` uses SciPy assignment to maximize valid match count then overlap at IoU≥0.5. Precision/recall and PQ include unmatched instances; matched-only Dice must not conceal misses. Origin `PUBLISHED_METHOD` for assignment; basis `ESTABLISHED`; threshold origin/basis `HEURISTIC`; evidence `PARTIAL`. Synthetic overlap tests verify arithmetic. A different truth filename alone does not establish independence. Split/merge hints use 10% overlap, also uncalibrated. Dense pair matrices have unbenchmarked large-N resource limits.

## D10. One-to-one cell–nucleus pairing

`cellular_measurements.py` uses sparse maximum-weight matching with dummy columns, prioritizing cardinality and then overlap among QC-eligible pairs. It measures full nucleus volume. Origin `PUBLISHED_METHOD` for the solver; basis `CONTEXT_DEPENDENT` for biological pairing; evidence `PARTIAL`. Large mixed integer IDs now preserve identity. This route assumes one-to-one pairing and may exclude multinucleated/anucleate cells; the multilevel route is a distinct alternative allowing multiple nuclei. No independent annotation supports one approach as preferable for all organoids. Radius-neighborhood density is censored at image borders and is not face-contact density.

## D11. Condition inference and exploratory analysis

`statistics/aggregation.py` gives each complete well equal weight within a biological replicate and bootstraps replicate summaries. `statistics/inference.py` instead fits eligible **object rows**, with a condition fixed effect and condition+replicate random intercept; singular/nonconverged fits fall back to replicate-clustered OLS. The estimands/weighting differ. Replicate IDs must be global across batches within condition, and independent across conditions; paired donors and well-level residual nesting are not modeled. Delimiter collisions were removed in the input worktree.

Volume is log10 transformed; contrasts have effects, SE, t-based CIs and within-feature BH correction. [statsmodels' test API](https://www.statsmodels.org/stable/generated/statsmodels.regression.mixed_linear_model.MixedLMResults.t_test.html) and [multiple-testing API](https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html) document mechanics. They do not establish that manually using t(G−1) for MixedLM provides calibrated small-sample coverage. LMM omnibus remains asymptotic. Origin `HEURISTIC` for these project choices, basis `CONTEXT_DEPENDENT`, evidence `INSUFFICIENT EVIDENCE` for inferential validity. No type-I-error/CI-coverage simulation, independent study design, or power/precision target is available. No replacement statistical method was selected.

`statistics/exploration.py` retains its normality/variance-driven t/Mann–Whitney selection, object-level splits, RF/logistic/XGBoost and KMeans workflows. Undefined effects, assumption tests and insufficient-sample results now abstain. Learned imputation/scaling occurs within training folds; model selection uses training CV. This follows [scikit-learn's leakage guidance](https://scikit-learn.org/stable/common_pitfalls.html#data-leakage). It does not fix biological-unit leakage from object splits. Cluster ANOVA on the features used to form clusters and cluster-prediction accuracy are circular descriptive diagnostics, not independent phenotype discovery. Shapiro requires [at least three observations](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.shapiro.html); that numerical minimum is not a justified biological sample size.

## D12. Reader and preview contracts

Two readers remain: `tiff_contract.py` (manifest/CLI, CZYX, OME spacing or explicit calibration) and `zstack_reader.py` (Web, ImageJ/TIFF/OME spacing, ZYX/CZYX). They are **not behaviorally identical**: the former does not infer ImageJ calibration; the latter has no UI time selector and rejects multi-time uploads. The Cellpose adapter alone allows metadata-free QYX/IYX as a legacy Z-stack interpretation, recorded without spacing. Time/series indices must be explicit valid integers. Unknown units, complex quantitative intensities, nonuniform/nonmonotonic Z, conflicting complete spacings and label export overflow reject. Only C=0 plane positions are checked; incomplete per-plane metadata, origins/directions and independent-channel registration remain NOT ASSESSED.

Preview percentile normalization, uint8 packing, stride downsampling, mesh coarsening and rendering operate on display payloads; raw source arrays feed segmentation/measurement. [SciPy zoom documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.zoom.html) defines center-grid geometry: resampling length N to M changes spacing by (N-1)/(M-1). Since 2026-09-11 the Cellpose adapter scales diameter/anisotropy by its reciprocal, using the actual rounded shape. Unequal rational X/Y scales reject before model execution because Cellpose accepts one XY anisotropy; full resolution remains available. Preview represents separate X/Y spacings and honors its array budget without a hard minimum 0.1 factor. Singleton collapse is rejected when rescaling would erase the center-span definition.

The existing linear-intensity / nearest-label center-grid interpolation was retained; downsampled intensities now use float32 output so integer inputs do not round interpolated values. Changing to full-pixel geometry or an unrecorded crop would alter the coordinate mapping and was not selected. Basis ESTABLISHED for the grid arithmetic, Origin SOFTWARE_OR_MODEL_DEFAULT for SciPy's grid convention; controlled verification is recorded in [the optimization report](WORKFLOW_OPTIMIZATION_2026-09-11.md). Actual model-output equivalence, anti-aliasing/parameter sensitivity and biological accuracy remain NOT ASSESSED. Valid coordinate scaling does not establish segmentation equivalence.

## D13. Actual Cellpose 4.2.1.1 model route

`cellpose_inference.py` invokes `CellposeModel(pretrained_model='cpdino-vitb')`, 3 orthogonal-plane flows and 3D dynamics, not a volumetric network trained/validated by this project. `do_3D=True` ignores flow-error thresholds; [Cellpose's 3D documentation](https://cellpose.readthedocs.io/en/latest/do3d.html) agrees with installed `dynamics.compute_masks`. The disabled UI controls now reflect this. The cell pass combines cytoplasm and nuclei channels; diameter suggestions use a separate coarse 2D model pass.

Origin `SOFTWARE_OR_MODEL_DEFAULT` for library behavior and `HEURISTIC` for project model/parameter selection; basis `CONTEXT_DEPENDENT`; evidence `PARTIAL` for execution and `INSUFFICIENT EVIDENCE` for object accuracy. Cached weights are hashed when found, but alternate model-cache roots/custom model objects can make recorded identity incomplete. Remaining hidden defaults are listed in PARAMETERS. Alternative foundation models/classical watershed are candidates, not validated superior/inferior choices. One 3×256×256 MPS smoke crop is not full-volume, annotation-based, CPU/GPU-equivalence, or repeatability validation.

## D14. Surface-area estimator comparison — 2026-09-11 (no production change)

A systematic surface V&V was frozen before any candidate result was seen
([plan](evidence/2026-09-11-measurement-vv/SURFACE_VV_PLAN.md), hashes in
`surface_vv_freeze.txt`). It covered 188 Gauss-digitized phantoms: sphere,
ellipsoid and flat-capped cylinder; 4 sizes or radii; ZYX spacings (1,1,1),
(2,1,1), (3,1,1) and (2,0.7,0.7); seeded orientations and sub-voxel offsets.
Every case was compared with a closed-form area.

The candidates were:

- **E0**, the production binary marching cubes;
- **E1**, a physical-space Gaussian with the predeclared σ = max spacing, plus
  0.5× and 1.5× variants that are sensitivity-only and not selectable;
- **E2**, a signed Euclidean distance field;
- **E3**, 13-direction discrete Crofton with spherical-Voronoi weights.

The acceptance rule was the unchanged §9 criterion: **every** case must have
|error| < 5%.

| Estimator | Worst |error| | Median signed | Strata passing |
|---|---|---|---|
| E0 production | 18.74% | +9.72% | none (fails even for large isotropic objects: 9.79%) |
| E1 σ = 1× | non-estimable in 17 cases | −6.64% | large (ρ≥12) only |
| E2 | 22.48% | +5.36% | none |
| E3 | 26.51% | −1.67% | sphere and ellipsoid with ρ≥3; cylinders fail in every stratum, matching the analytical plane response |

**Superseded 2026-09-12 — read D3a.** The conclusion below is the honest
outcome of round 1 and is kept as the record of it. A round-2 development
cycle, with the confirmation set still untouched, produced a qualified
estimator that is now the production default; the replacement rule was met
within a declared domain rather than unrestrictedly.

No candidate met the predeclared replacement rule. None was better in every
resolution stratum. E1's outcome depends on σ, which has no analytical
provenance. The production estimator and the historical definition of
`surface_area_um2` are therefore **unchanged**. The confirmation grid was left
unseen. Suitability for the §9 criterion remains `NOT SUPPORTED` / `FAIL`. The
estimator is now identified by `SURFACE_AREA_METHOD`
(`marching_cubes_binary_lewiner_v1`) in the classical, multilevel and Web
feature provenance, and a version-lock test detects silent changes. Options
that need an owner decision are listed in
[the selection record](evidence/2026-09-11-measurement-vv/SURFACE_SELECTION.md).
