# Algorithm evidence audit — 2026-09-10

This record describes the current implementation; it does not retrospectively justify its selection. No scientific method was replaced or tuned in this audit. `ESTABLISHED` describes a mathematical/method basis, not suitability for these biological images. Project-data suitability remains `INSUFFICIENT EVIDENCE` unless explicitly bounded below. Installed versions and executable defaults are captured in [installed-methods.txt](evidence/2026-09-10/installed-methods.txt).

## D1. Parent assignment by maximum voxel overlap

`quantification/multilevel_relationships/hierarchy.py` counts label co-occurrences, takes the largest count, and resolves ties to the smallest parent ID. This is an exact operational containment rule, not a calibrated biological assignment model. Origin `HEURISTIC`; basis `CONTEXT_DEPENDENT`; evidence `PARTIAL` (hand-counted assignments and large-ID tests pass). Nearest-centroid assignment is a meaningful alternative for different data, but was not benchmarked. Misregistration, fragments, ties and cross-parent objects require review; same shape/spacing cannot establish registration. No qualified real hierarchy annotations are available.

## D2. Nucleus-organoid membership inherited through the cell

`workflows/multilevel_measurement_workflow.py` inherits the organoid through the cell and retains direct nucleus–organoid overlap as an audit field. Origin/basis `HEURISTIC`; evidence `PARTIAL`. Tests verify transitivity and mismatch flags. Direct overlap is an alternative definition, not an empirically rejected inferior method. No biological transitivity/assignment performance claim is supported.

## D3. Native-spacing marching-cubes surface

`quantification/features.py::surface_mesh` uses scikit-image 0.26.0's **Lewiner** implementation (the installed default), level 0.5, step 1, padded zero exterior and `allow_degenerate=False`. Lorensen–Cline is historical background, not the exact implementation. [Official method documentation](https://scikit-image.org/docs/stable/api/skimage.measure.html#skimage.measure.marching_cubes) identifies the algorithm and axis-ordered spacing. Origin `PUBLISHED_METHOD`; basis `ESTABLISHED`; implementation evidence `PARTIAL`; suitability for the specification's <5% surface error is `NOT SUPPORTED` / `FAIL`: the radius-18 µm sphere at (2,1,1) µm has 11.7347% area error. Sphericity is the dimensionless isoperimetric expression; it is not clipped. Changing smoothing, surface estimator or tolerance requires a new scientific decision and validation, so none was changed.

Classical organoid and Web object-feature geometry fill enclosed holes; multilevel geometry uses raw label voxels. These are different estimands. Classical `cells` reports raw cell volume plus envelope geometry; Web centroids/solidity use raw regionprops despite envelope size/axes. Those differences must be retained in interpretation. Border padding closes truncated objects computationally, not biologically. `mask_features` does not apply the classical eligibility filter. Alternative surface estimators and acquisition-dependent bias correction are NOT ASSESSED.

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

Preview percentile normalization, uint8 packing, stride downsampling, mesh coarsening and rendering operate on display payloads; raw source arrays feed segmentation/measurement. [SciPy zoom documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.zoom.html) distinguishes center-grid and full-pixel geometry: the current Cellpose downsample adapter scales diameter/anisotropy by the requested factor while rounded array sizes and center-grid resampling can differ. Downsampled scientific equivalence is NOT ASSESSED, especially for small/odd XY sizes. No interpolation method was changed.

## D13. Actual Cellpose 4.2.1.1 model route

`cellpose_inference.py` invokes `CellposeModel(pretrained_model='cpdino-vitb')`, 3 orthogonal-plane flows and 3D dynamics, not a volumetric network trained/validated by this project. `do_3D=True` ignores flow-error thresholds; [Cellpose's 3D documentation](https://cellpose.readthedocs.io/en/latest/do3d.html) agrees with installed `dynamics.compute_masks`. The disabled UI controls now reflect this. The cell pass combines cytoplasm and nuclei channels; diameter suggestions use a separate coarse 2D model pass.

Origin `SOFTWARE_OR_MODEL_DEFAULT` for library behavior and `HEURISTIC` for project model/parameter selection; basis `CONTEXT_DEPENDENT`; evidence `PARTIAL` for execution and `INSUFFICIENT EVIDENCE` for object accuracy. Cached weights are hashed when found, but alternate model-cache roots/custom model objects can make recorded identity incomplete. Remaining hidden defaults are listed in PARAMETERS. Alternative foundation models/classical watershed are candidates, not validated superior/inferior choices. One 3×256×256 MPS smoke crop is not full-volume, annotation-based, CPU/GPU-equivalence, or repeatability validation.
