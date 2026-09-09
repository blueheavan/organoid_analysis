# Algorithm Decisions — Organoid Pipeline

Version: 1.1.0
Date: 2026-09-03
Applies to: multilevel 3D analysis (`analysis.analyze-3d`), classical organoid morphology/viability analysis, segmentation-validation instance matching, cell–nucleus pairing (`analysis cells`), and cross-condition statistical testing.

## Decision provenance

This document records the scientific rationale for consequential algorithmic choices in the pipeline. For each decision: the scientific task, credible alternatives considered, the selected method, the alternatives rejected and why, and the validation status. It is a living record and is updated when a method is added or changed.

---

## D1. Parent assignment by maximum voxel overlap

**Task:** Assign each cell to an organoid, and each nucleus to a cell, within a registered 3D instance hierarchy.

**Alternatives considered:**
- Nearest centroid (Euclidean) assignment.
- Boundary/contact overlap.
- Voxel-count majority.

**Selected:** Maximum voxel overlap (argmax of shared voxel count), with deterministic tie-break to the smallest parent ID.

**Rationale:** Voxel overlap measures physical co-occupancy of the child within each parent. For nested 3D structures with potentially non-convex shapes, overlap reflects where the object's mass actually lies. Centroid distance can mis-assign a large, elongated or fragmented child whose centroid falls outside its true parent (e.g., a nucleus that protrudes across a boundary). Overlap is robust to such geometry and is a well-established approach for label-map ancestry.

**Rejected:**
- Nearest centroid — fails for elongated/border-crossing objects; the centroid may lie in the neighbor's volume.
- Boundary contact — measures surface adjacency, not containment; two overlapping objects may share no boundary face.

**Failure modes:** A child straddling multiple parents gets the parent with plurality overlap. This is surfaced to the user via `crosses_multiple_parents`, `parent_overlap_fraction < low_parent_overlap_fraction`, and `parent_assignment_failed` QC flags rather than being silently dropped.

**Validation:** `PARTIAL — PASS for algorithmic correctness on controlled phantoms (test_multilevel_measurement_workflow.py); real-world assignment accuracy against expert annotation NOT ASSESSED.`

**Strength labeling:** established method (maximum overlap assignment).

---

## D2. Nucleus-organoid membership inherited through the cell

**Task:** Determine the organoid that contains each nucleus.

**Selected:** The nucleus's organoid is the organoid assigned to its parent cell. The nucleus-to-organoid direct overlap is retained as an audit field (`direct_organoid_id`) but never allowed to override the cell-mediated hierarchy.

**Rationale:** The two-level hierarchy Cell→Organoid, Nucleus→Cell defines a deterministic nesting. A nucleus should belong to the organoid of the cell that houses it, even if the nucleus's own voxels straddle an organoid boundary. This enforces a consistent interpretation of "this nucleus belongs to this cell, which belongs to this organoid."

**Rejected:** Direct nucleus→organoid overlap as the governing assignment — this would break the transitive cell-chain consistency and create contradictory memberships.

**Failure modes:** A mismatch between the cell-derived organoid and the direct nucleus overlap is recorded as the `direct_organoid_parent_mismatch` QC flag, not silently resolved.

**Validation:** `PASS — controlled phantom (test_multilevel_measurement_workflow.py::test_nucleus_inherits_organoid_from_assigned_cell_and_retains_direct_overlap_audit).`

**Strength labeling:** heuristic (design choice to enforce a consistent hierarchy).

---

## D3. True-3D morphology via marching cubes level-0.5 surface

**Task:** Compute physical surface area and sphericity from a 3D label mask.

**Selected:** `skimage.measure.marching_cubes` at `level=0.5`, native physical spacing, step size 1, one-voxel zero padding, no mesh smoothing. Volume = voxel count × voxel volume. Sphericity = `π^(1/3) (6V)^(2/3) / A`.

**Rationale:** Marching cubes at native spacing with no smoothing preserves the underlying voxel geometry and is the standard surface estimator for binary masks. Zero padding closes the crop surface so an object fully inside the image reports a closed surface. SPHERICITY uses the classic isoperimetric ratio; perfect sphere ⇒ 1.0.

**Original literature:** Lorensen, W. E., & Cline, H. E. (1987). Marching cubes: A high resolution 3D surface construction algorithm. *ACM SIGGRAPH Computer Graphics*, 21(4), 163–169. https://doi.org/10.1145/37401.37422

**Rejected:**
- MIP / single-slice / density projections — explicitly excluded by project principle (never treat 2D projection as 3D measurement).
- Smoothed or downsampled mesh for measurements — used only for preview, never for reported area. Documented discretization bias is retained rather than hidden.

**Failure modes / known bias:** Voxel discretization inflates area and can push sphericity slightly above 1.0. Values >1.05 are flagged (`sphericity_above_geometric_range`), never silently clipped into [0,1]. Empirical check on a 15-voxel-radius sphere (spacing 1): volume ratio ≈ 1.0007, area ratio ≈ 1.0895, sphericity ≈ 0.9183. This is the documented discretization bias, not a defect.

**Validation:** `PARTIAL — PASS for volume/centroid/axis analytical agreement on phantoms; surface-area absolute accuracy on real biological organoids NOT ASSESSED (requires independent surface measurement reference).`

**Strength labeling:** established method (skimage marching cubes surface).

---

## D4. Face-contact topology (cell adjacency), not centroid/kNN/Delaunay

**Task:** Define which cells are "neighbors" for a contact graph.

**Selected:** Voxel-face adjacency only. An edge exists only when two cell labels share at least one voxel face. Edge weight = physical contact area computed from the anisotropic face area (Z-face → Y×X, Y-face → Z×X, X-face → Z×Y).

**Rationale:** Face contact is a direct statement of physical adjacency/contact. Centroid distance, k-nearest-neighbor, and Delaunay triangulation infer proximity, not contact, and are deliberately not treated as contact.

**Rejected:** centroid-distance, kNN, Delaunay-based adjacency for contact claims.

**Validation:** `PASS — controlled phantom verifies X-face contact area = 25 × Z×Y for a 5×5×5 cell pair with 5×5 shared X-face (test_multilevel_measurement_workflow.py).`

**Strength labeling:** defined operational definition.

---

## D5. Robust MAD volume-outlier detection

**Task:** Flag unusually large/small organoids/cells/nuclei by volume.

**Selected:** Robust z-score: `0.67448975 × (V − median) / MAD` with a threshold of 3.5 (constant `mad_z_threshold`), i.e. equivalent to ~3.5 z-units on the robust scale (`MAD * 1.4826`). With MAD = 0 (common in near-identical synthetic populations), finite values differing from the median are flagged.

**Rationale:** MAD is robust to outliers, unlike mean/std which are inflated by the very outliers being detected. This prevents a single extreme object from masking the next-level outliers.

**Original literature:** The 1.4826 (equivalently, its reciprocal 0.67448975) consistency constant that rescales MAD into a σ-equivalent robust estimator under a Gaussian assumption is documented in Rousseeuw, P. J., & Croux, C. (1993). Alternatives to the Median Absolute Deviation. *Journal of the American Statistical Association*, 88(424), 1273–1283. https://doi.org/10.1080/01621459.1993.10476408. The constant itself is the reciprocal of the 0.75-quantile of the standard normal distribution, Φ⁻¹(0.75). `src/organoid_analysis/quantification/multilevel_relationships/qc.py` uses the reciprocal form (0.67448975 × …); `src/organoid_analysis/quantification/features.py` and `src/organoid_analysis/phenotyping/viability.py` use the direct form (… × 1.4826); both are the same constant, named and cross-referenced at each definition site so the two spellings do not appear to be different numbers.

**Validation:** `PASS — controlled phantoms; real-data outlier calibration NOT ASSESSED.`

**Strength labeling:** established robust-statistics technique (median + MAD).

---

## D6. Radial position normalized to equivalent-sphere radius

**Task:** Locate a cell within its organoid (core vs periphery).

**Selected:** `normalized_radial_position_equivalent_radius = distance(cell_centroid, organoid_centroid) / r_eq`, where `r_eq = (3V/(4π))^(1/3)` is the equivalent-sphere radius of the organoid's volume.

**Rationale:** Normalizing to an equivalent sphere yields a dimensionless 0..~1+ index. Core = ≤0.5, periphery = ≥0.8 (configurable `core_max_normalized_radial_position`, `peripheral_min_normalized_radial_position`).

**Rejected:** Local-radius normalization — the codebase does not estimate a per-direction local organoid radius and does not claim to.

**Validation:** `PASS — computing path is deterministic/gradient-tested; biological meaning NOT ASSESSED.`

**Strength labeling:** heuristic; explicitly documented as equivalent-radius, not local-radius.

---

## D6b. Moment-equivalent principal axes (ellipsoid fit from voxel second moment)

**Task:** Compute principal axis lengths (major, intermediate, minor) for organoid/cell/nucleus envelope geometry.

**Selected:** Voxel-centroid-centered second-moment matrix with intrinsic voxel moment included: `covariance = centered.T @ centered / n + diag(spacing²/12)`, where `centered` are voxel coordinates relative to centroid scaled by spacing. Eigenvalues λ sorted descending; axis lengths = `2√(5λ)`. The `+ spacing²/12` term accounts for each voxel's intrinsic uniform-density second moment, making the axes describe a moment-equivalent ellipsoid rather than the point-cloud covariance alone.

**Rationale:** The uniform-density voxel has a known intrinsic second moment (I = side⁵/12 per axis), which shifts the covariance to represent the solid object's inertia tensor. The factor `2√(5λ)` converts the second-moment eigenvalues to the equivalent ellipsoid's principal axis lengths (for a uniform ellipsoid, `I_xx = M/5 (b²+c²)`, so axis `a = √(5λ_a)` for the major axis from the diagonalized inertia). This is a standard morphometric definition used in 3D shape analysis.

**Rejected:** Raw point-cloud covariance (without intrinsic voxel moment) — systematically underestimates axis lengths for small objects; single-axis extent (max-min) — ignores shape anisotropy.

**Validation:** `PASS — analytical ellipsoid test (tests/quantification/test_cellular_measurements.py:136–147) verifies major/intermediate/minor = 12/8/6 µm for a 6×4×3 µm spacing-2 ellipsoid; values are pinned by the independent analytical expectation.`

**Strength labeling:** established morphometric definition (moment-equivalent ellipsoid).

---

## D7. Classical morphology marker measurements & viability state gating

**Task:** Classical organoid-level Calcein/PI viability-like states.

**Selected:** Marker means measured on raw (unnormalized) arrays within the outer envelope; background = median in a physical shell; corrected mean = foreground mean − background median; saturation detection uses an explicit camera limit or integer dtype max; control-scaled endpoints learned per batch from live/dead controls; marker scaled as `(mean − low)/(high − low)`; states gated at 0.60 (high) and 0.30 (low) with indeterminate between.

**Rationale, rejected alternatives, and caveats are fully documented in docs/METHODS.md.** This is an intensity-based signal-pattern characterization, not a cell-count or a validated cell-viability fraction.

**Condition-scoping fix (2026-09-03 audit fix):** An independent audit found `calibrate()` grouped controls by `(control, biological_replicate)` without `condition`, unlike `stats.py`'s composite-key convention — since `biological_replicate` labels (e.g. "R1") are only unique *within* a condition, a batch whose live or dead controls spanned more than one `condition` value could have been silently pooled as if they were one replicate, corrupting the batch's calibration endpoints without any error or QC flag. Fixed: `calibrate()` now rejects (status `unavailable`, reason `controls_span_multiple_conditions`) a batch whose live or dead controls span more than one condition, rather than guessing which grouping was intended. Regression test: `test_controls_spanning_multiple_conditions_are_rejected` in `test_viability_summary.py`.

**Validation:** `PARTIAL — controlled phantom agreement on simplified marker patterns; the condition-scoping fix above is unit-tested; biological calibration against matched controls and independent viability measurements NOT ASSESSED.`

**Strength labeling:** empirically calibrated / heuristic.

---

## D8. Otsu global thresholding + physical-distance watershed splitting (classical segmentation)

**Task:** Produce a 3D instance-label mask from a single structural fluorescence/brightfield channel (`segmentation.method: watershed`).

**Selected:** Global Otsu thresholding on the (optionally background-subtracted, Gaussian-smoothed) structural volume, followed by a spacing-aware Euclidean distance transform, h-maxima seed detection at a physical prominence, physical minimum-seed-separation suppression, and masked watershed of the negative distance map.

**Rationale:** Otsu thresholding is a standard, parameter-free global binarization method. Distance-transform watershed is a standard technique for splitting touching, roughly convex objects (such as organoids) without a trained model.

**Original literature:**
- Otsu, N. (1979). A threshold selection method from gray-level histograms. *IEEE Transactions on Systems, Man, and Cybernetics*, 9(1), 62–66. https://doi.org/10.1109/TSMC.1979.4310076
- Vincent, L., & Soille, P. (1991). Watersheds in digital spaces: an efficient algorithm based on immersion simulations. *IEEE Transactions on Pattern Analysis and Machine Intelligence*, 13(6), 583–598. https://doi.org/10.1109/34.87344

**Rejected:** A trained/learned segmentation model for the classical route — deliberately out of scope; this route is the parameter-free/classical alternative to Cellpose (see `src/organoid_analysis/segmentation/`).

**Validation:** `PARTIAL — unit/synthetic-phantom tests cover the split/threshold logic (test_segmentation.py); real-image detection accuracy is assessed only via the optional independent-annotation route in docs/METHODS.md, not a general claim.`

**Strength labeling:** established methods (Otsu; distance-transform watershed), combined by heuristic parameter choices (seed height, minimum seed separation — see docs/PARAMETERS.md).

**Code location:** `src/organoid_analysis/segmentation/watershed_instances.py::segment`, `::watershed_instances`.

---

## D9. Hungarian (Kuhn–Munkres) one-to-one instance matching for segmentation validation

**Task:** Match predicted instance labels to an independent annotated truth mask for validation metrics (precision/recall/Dice/IoU/panoptic quality), and match synthetic phantom predictions to their known ground truth.

**Selected:** `scipy.optimize.linear_sum_assignment` (Hungarian/Kuhn–Munkres algorithm) maximizing a reward that first maximizes the count of IoU ≥ threshold matches, then total IoU, at a fixed `iou_threshold` (default 0.5).

**Rationale:** One-to-one optimal assignment prevents a many-to-one or greedy match from hiding split/merge errors that a foreground-only overlap metric (e.g. plain Dice) would miss.

**Original literature:** Kuhn, H. W. (1955). The Hungarian method for the assignment problem. *Naval Research Logistics Quarterly*, 2(1–2), 83–97. https://doi.org/10.1002/nav.3800020109

**Rejected:** Greedy nearest-IoU matching — can produce inconsistent, order-dependent assignments when several predictions compete for the same truth object.

**Validation:** `PASS — controlled phantom validation (test_analysis.py, run_demo synthetic validation); real-data annotation matching depends on the availability of an independent truth mask (INSUFFICIENT EVIDENCE without one).`

**Strength labeling:** established method (Hungarian assignment); `iou_threshold=0.5` is a conventional choice (see docs/PARAMETERS.md), not independently calibrated for this pipeline.

**Code location:** `src/organoid_analysis/validation/segmentation_metrics.py::match_instances`.

---

## D10. Minimum-weight full bipartite matching for cell–nucleus pairing

**Task:** Pair each cell instance with at most one nucleus instance from two independently segmented, registered label masks (`analysis cells` CLI route, distinct from the `analyze-3d` maximum-overlap hierarchy in D1).

**Selected:** `scipy.sparse.csgraph.min_weight_full_bipartite_matching` on a sparse graph weighted by (cardinality-maximizing bonus + voxel overlap), restricted to candidate pairs that already pass the nucleus-containment/size/N:C-ratio QC gates; every cell gets a guaranteed dummy column so the matching is always full.

**Rationale:** Optimal bipartite matching first maximizes the number of valid one-to-one pairs, then total overlap, which avoids a greedy "closest nucleus wins" rule silently mis-pairing cells that compete for the same nucleus.

**Original literature:** This is an application of the same assignment-problem theory as D9 — Kuhn, H. W. (1955). *Naval Research Logistics Quarterly*, 2(1–2), 83–97. https://doi.org/10.1002/nav.3800020109

**Rejected:** Greedy dominant-overlap pairing (assign each cell its single largest-overlap nucleus without considering competing cells) — can double-assign one nucleus to two cells' "best candidate" lists without resolving the conflict optimally.

**Validation:** `PARTIAL — unit-tested on synthetic label pairs (test_cellular.py per docs/PARAMETERS.md gap noted below); real-data pairing accuracy NOT ASSESSED.`

**Strength labeling:** established method (bipartite assignment); the QC gate thresholds (`min_cell_volume_um3`, `min_nucleus_volume_um3`, `max_nc_ratio`, `min_nucleus_containment`) are heuristic (see docs/PARAMETERS.md).

**Code location:** `src/organoid_analysis/quantification/cellular_measurements.py::pair_and_filter_cells`.

---

## D11. Linear mixed-effects model with BH-FDR pairwise contrasts for condition comparison

**Task:** Test whether a continuous morphology feature (volume, sphericity) differs across experimental conditions while accounting for biological-replicate structure (`src/organoid_analysis/statistics/inference.py`, invoked from `pipeline.py` when `cfg["stats"]["enabled"]`).

**Selected:** A linear mixed-effects model (condition as fixed effect, `condition::biological_replicate` as a random intercept) fit by REML, giving an omnibus Wald test and all-pairwise contrasts between conditions. Falls back to OLS with replicate-clustered standard errors when the random-effects fit is singular (variance below `1e-6 × scale`) or does not converge. Pairwise p-values are Benjamini–Hochberg FDR corrected. `volume_um3` is log10-transformed before fitting; `sphericity` is not.

**Rationale:** A random intercept on the biological-replicate unit avoids treating pooled organoids/technical replicates as independent samples (pseudoreplication). BH-FDR correction controls the false discovery rate across the pairwise contrasts performed for each feature. Log-transforming volume addresses its right-skewed, multiplicative-scale distribution before a model that assumes approximately normal, homoscedastic residuals.

**Original literature:** Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*, 57(1), 289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

**Rejected:** An unpaired t-test/ANOVA on pooled per-organoid values — would treat organoids as independent biological replicates (pseudoreplication), overstating significance.

**Documentation note:** `docs/METHODS.md` previously stated "No automatic inferential tests are included," which was inaccurate once this module shipped; that statement has been corrected (see `docs/METHODS.md`, Treatment summaries) to describe this module and its assumptions.

**Failure modes:** A feature/condition combination with fewer than `stats.min_replicates_per_condition` (default 3) biological replicates in a condition is silently excluded from that feature's test rather than raising; callers must treat an empty result as "not enough data to test." (2026-09-03 audit fix: `organoid_measurement_workflow.py` now filters `objects` through `complete_unit_objects()` before calling `condition_pairwise_tests`, matching `aggregation.py`'s own exclusion of organoids from partially-failed acquisition units; regression test `test_complete_unit_objects_excludes_incomplete_units` in `test_organoid_measurement_workflow.py`.)

**Small-sample p-value correction (2026-09-03 audit fix):** An independent scientific-software audit found that both branches originally reported Wald p-values against an asymptotic z reference, which is anti-conservative (overstates significance) at the pipeline's own documented minimum of 3 replicates/condition (empirically confirmed against this repository's statsmodels 0.14.6). Fixed: the OLS-fallback branch now fits with `use_t=True`, which statsmodels resolves to a cluster-robust t(G−1) reference (G = number of replicate clusters) — the standard small-cluster correction (Cameron & Miller, 2015, *Journal of Human Resources*, 50(2), 317–372). The LMM branch has no equivalent built-in correction in statsmodels (no Satterthwaite/Kenward-Roger for `MixedLM`), so `_pairwise_contrasts` now manually applies the same t(G−1) reference to its pairwise contrasts. **Residual limitation:** this is a standard but approximate small-cluster correction, not full Satterthwaite/Kenward-Roger; and the LMM branch's *omnibus* Wald test (`omnibus_p` in `stats_results.json`, flagged via the `omnibus_p_small_sample_corrected` field) remains asymptotic/uncorrected — treat it as a rough screening result, not confirmatory. The pairwise, BH-FDR-corrected contrasts in `pairwise_contrasts.csv` are the primary, corrected output.

**Validation:** `PARTIAL — the small-sample p-value correction above is independently derived and verified (regression tests in test_inference.py assert the corrected p-value against a hand-computed t(G-1) reference, and that it is strictly more conservative than the uncorrected z reference); no independent review against an external Satterthwaite/Kenward-Roger reference implementation (e.g. R's lmerTest) has been performed, and the omnibus test remains NOT ASSESSED for small-sample validity.`

**Strength labeling:** established statistical methods (LMM; BH-FDR), combined by an engineering fallback rule (LMM→OLS singular-fit threshold) that is heuristic and not independently validated.

**Code location:** `src/organoid_analysis/statistics/inference.py::condition_pairwise_tests`, `::fit_model`.

---

## D12. TIFF axis/spacing readers kept as two implementations, with behavior (not code) unified

**Task:** Two TIFF-reading paths exist -- `microscopy_io/tiff_contract.py` (manifest-driven, multi-role: structure/calcein/pi/probability/labels, for the classical CLI) and `microscopy_io/zstack_reader.py` + `metadata.py` (single uploaded file, for Streamlit). Each independently implements axis canonicalization and OME/ImageJ spacing parsing, which is a real duplication risk: the two had drifted to disagree on two consequential behaviors before this decision.

**Found divergences (2026-09-08 audit):**
1. `tiff_contract.ome_spacing()` already rejected a nonuniform Z-plane grid (via `PositionZ` values); `metadata.parse_spacing_ome()` (the Streamlit path) did not check this at all, so a Streamlit-uploaded stack with unevenly spaced Z planes would silently get one (wrong) Z spacing value instead of an error.
2. `tiff_contract.canonical_czyx()` silently squeezed away any unrecognized axis of size 1 (e.g. a stray legacy/filler dimension); `zstack_reader._reorder_to_czyx()` already rejected *any* unrecognized axis regardless of size ("Fully explicit, no guessing").

**Selected:** Fixed both divergences toward the stricter behavior rather than merging the two readers into one implementation:
- Ported the nonuniform-Z-grid check into `metadata.parse_spacing_ome()` (using `ome_types`' `Pixels.planes`, mirroring `ome_spacing()`'s own ElementTree-based check).
- Removed `canonical_czyx()`'s size-1 exception so it always rejects an unrecognized axis, matching `_reorder_to_czyx()`.
- Consolidated the two independent unit-conversion tables and the two independent `SPACING_RTOL`/`SPACING_ATOL_UM` constant definitions into one copy in `metadata.py` (`to_um`), which `tiff_contract.py` now imports rather than redefining -- the table became the union of both (adding the µ/μ codepoint variants and the "micrometer" spelling that only one side previously accepted).

**Rejected (deferred, not abandoned):** A full merge into one canonical reader/one OME-parsing implementation. Investigated switching `tiff_contract.ome_spacing()`'s hand-rolled `xml.etree.ElementTree` parsing to `metadata.py`'s `ome_types`-based parsing, and found `ome_types` performs strict pydantic validation of unit strings that the OME schema's controlled vocabulary does not actually require real-world writers to follow exactly (e.g. `PhysicalSizeZUnit="um"` -- ASCII, not the canonical "µm" -- raises an uncaught `pydantic.ValidationError` under `ome_types`, but parses fine under raw ElementTree). Switching would trade the classical pipeline's current tolerance of such files for the Streamlit path's `ome_types` robustness elsewhere, which is not a clear improvement and was not evaluated against real acquisition files from either code path's actual user base. The manifest-vs-single-file/multi-role-vs-none input shapes also do not collapse into one API without changing either path's public contract. Two implementations are kept, with the specific behaviors above now enforced identically and verified by shared-intent regression tests, rather than forcing a structural merge that would need to pick a winner on an open question (unit-string strictness) with no evidence to decide it.

**Validation:** Both changed behaviors have regression tests: `tests/microscopy_io/test_io.py::test_ambiguous_singleton_axis_is_rejected_not_silently_squeezed` and `tests/visualization/test_3d_preview.py::SpacingMetadataTests::test_nonuniform_z_positions_are_rejected`/`test_uniform_z_positions_are_accepted`. All 5 real sample images in `data/images/` were re-loaded through `load_zstack()` after the change with no regression (none use ambiguous axes or nonuniform Z metadata, so none were expected to be affected, and none were).

**Code location:** `src/organoid_analysis/microscopy_io/tiff_contract.py::canonical_czyx`, `::ome_spacing`; `src/organoid_analysis/microscopy_io/metadata.py::parse_spacing_ome`, `::to_um`.

---

## Rejected overall approaches

- **Replace the existing segmentation model** — preserved by requirement. Segmentation is upstream; this pipeline consumes its output.
- **Use centroid-based hierarchy** — rejected as less robust than overlap (D1).
- **Use 2D projections for any quantitative 3D measurement** — rejected by project principle.

## Traceability

| Decision | Code location |
|---|---|
| D1 | `src/organoid_analysis/quantification/multilevel_relationships/hierarchy.py` |
| D2 | `src/organoid_analysis/workflows/multilevel_measurement_workflow.py` (lines ~193-212) |
| D3 | `src/organoid_analysis/quantification/features.py::geometry`, `::surface_mesh` |
| D4 | `src/organoid_analysis/quantification/multilevel_relationships/topology.py` |
| D5 | `src/organoid_analysis/quantification/multilevel_relationships/qc.py::_volume_outliers` |
| D6 | `src/organoid_analysis/quantification/multilevel_relationships/spatial.py` |
| D7 | `src/organoid_analysis/quantification/features.py::marker_measurements`, `src/organoid_analysis/phenotyping/viability.py`, `src/organoid_analysis/result_export/report.py` |
| D8 | `src/organoid_analysis/segmentation/watershed_instances.py::segment`, `::watershed_instances` |
| D9 | `src/organoid_analysis/validation/segmentation_metrics.py::match_instances` |
| D10 | `src/organoid_analysis/quantification/cellular_measurements.py::pair_and_filter_cells` |
| D11 | `src/organoid_analysis/statistics/inference.py::condition_pairwise_tests`, `::fit_model` |
| D12 | `src/organoid_analysis/microscopy_io/tiff_contract.py::canonical_czyx`, `::ome_spacing`; `src/organoid_analysis/microscopy_io/metadata.py::parse_spacing_ome`, `::to_um` |
