# Parameters Catalog — Organoid Pipeline

Version: 1.1.0
Date: 2026-09-03

This catalog records the provenance of every scientifically consequential parameter, threshold, and default in the pipeline. Properties follow the rigor protocol (section 6). Engineering-only parameters that cannot affect scientific results (e.g. output verbosity, file naming) are not catalogued.

Legend for status: `Validated` / `Partially validated` / `Not validated`; sensitivity column: `Tested` / `Not tested` / `N/A`.

---

## Multilevel 3D analysis parameters

Source of truth: `src/organoid_analysis/quantification/multilevel_relationships/config.py` (dataclass) and CLI flags in `src/organoid_analysis/command_line/organoid_commands.py`.

| Parameter | Value | Unit | Origin | Rationale | Scientific impact | Calibration artifact | Status | Sensitivity | User configurable |
|---|---|---|---|---|---|---|---|---|---|
| `minimum_voxels` | 5 | voxel | heuristic (engineer-defined) | Objects with <5 voxels are too small to be reliably measured / non-degenerate | Flags too-small objects as `too_small` QC; affects QC only, not measurement set | N/A (validation only) | Partially validated | Not tested | Yes (`--minimum-voxels`) |
| `low_parent_overlap_fraction` | 0.5 | fraction (of child voxels) | heuristic | <50% of a child's voxels overlapping its assigned parent suggests a weak/ambiguous assignment | Flags `low_parent_overlap` QC | N/A | Partially validated | Not tested | Yes |
| `mad_z_threshold` | 3.5 | robust z-units | heuristic | ~3.5 σ beyond median/MAD identifies extremes without catastrophic masking | Flags `volume_outlier` QC | N/A | Partially validated | Not tested | Yes |
| `core_max_normalized_radial_position` | 0.5 | dimensionless | heuristic | Equivalent-radius ≤0.5 = core | Defines `core_cell_fraction` | N/A | Not validated (operational definition) | Not tested | Yes |
| `peripheral_min_normalized_radial_position` | 0.8 | dimensionless | heuristic | Equivalent-radius ≥0.8 = periphery | Defines `peripheral_cell_fraction` | N/A | Not validated (operational definition) | Not tested | Yes |

---

## Classical analysis config defaults

Source of truth: `src/organoid_analysis/config.py` (module `DEFAULTS`), validated by `validate_config`.

### Segmentation

| Parameter | Value | Unit | Origin | Rationale / impact | Status | Sensitivity |
|---|---|---|---|---|---|---|
| `method` | `watershed` | – | established (skimage watershed) | Classical threshold+watershed split | Validated (algorithm documented) | N/A |
| `polarity` | `bright` | – | heuristic | Dark-polarity is an exploratory brightfield baseline | Partially validated | N/A |
| `gaussian_sigma_um` | 1.0 | µm | heuristic | Pre-smoothing scale; smooths then segments | Not validated on real data | Not tested |
| `threshold` | `otsu` | – | established (Otsu) | Global automatic threshold | Validated (established statistic) | N/A |
| `probability_threshold` | 0.5 | probability | conventional | Foreground-probability decision boundary | Partially validated | Not tested |
| `closing_radius_um` | 1.0 | µm | heuristic | Binary closing scale | Not validated | Not tested |
| `fill_enclosed_holes` | `true` | – | heuristic | Envelope fills internal lumens (documented trade-off) | Partially validated | N/A |
| `min_volume_um3` | 2000.0 | µm³ | heuristic | Minimum object volume to retain | Not validated on real data | Not tested |
| `seed_h_um` | 2.0 | µm | heuristic | H-maxima seed height | Not validated | Not tested |
| `seed_min_distance_um` | 15.0 | µm | heuristic | Physical seed suppression distance | Not validated | Not tested |
| `split_touching` | `true` | – | heuristic | Enable watershed splitting | Partially validated | N/A |
| `max_foreground_fraction` | 0.70 | fraction | heuristic | Dense foreground ⇒ QC review flag | Partially validated | N/A |
| `qc_reference_method` | `none` | – | heuristic | Separate watershed as agreement check; `none` disables | Partially validated | N/A |
| `qc_count_difference_threshold` | 0.40 | fraction | heuristic | Count disagreement triggers review | Partially validated | N/A |
| `qc_min_z_extent_ratio` | 0.65 | fraction | heuristic | Z-extent disagreement threshold | Partially validated | N/A |
| `_HIGH_ANISOTROPY_RATIO` (segmentation.py, not YAML-configurable) | 5 | ratio (max/min voxel spacing) | heuristic | Raises `high_voxel_anisotropy` QC flag; relates to docs/SCIENTIFIC_SPEC.md's "Isotropic voxels" limitation | Not validated | Not tested |

### Quality

| Parameter | Value | Unit | Origin | Rationale / impact | Status | Sensitivity |
|---|---|---|---|---|---|---|
| `exclude_border` | `true` | – | established | Border-truncated objects unreliable | Partially validated | N/A |
| `min_z_slices` | 5 | slice | heuristic | Too few Z ⇒ insufficient 3D sampling | Partially validated | Not tested |
| `max_volume_um3` | `null` | µm³ | – | Optional upper exclusions | N/A | N/A |
| `background_inner_um` | 2.0 | µm | heuristic | Background shell inner radius | Not validated | Not tested |
| `background_outer_um` | 7.0 | µm | heuristic | Background shell outer radius | Not validated | Not tested |
| `min_background_voxels` | 100 | voxel | heuristic | Minimum local background to trust background stats | Not validated | Not tested |
| `max_saturated_fraction` | 0.01 | fraction | heuristic | Saturation >1% ⇒ invalidate marker | Partially validated | Not tested |

### Viability

| Parameter | Value | Unit | Origin | Rationale / impact | Status | Sensitivity |
|---|---|---|---|---|---|---|
| `mode` | `uncalibrated` | – | heuristic | No controls ⇒ uncalibrated states | Partially validated | N/A |
| `min_control_replicates` | 2 | replicate | statistical (min for separation) | Minimum independent control replicates | Partially validated | N/A |
| `min_control_separation_snr` | 3.0 | – | heuristic | Endpoint separation relative to noise | Partially validated | Not tested |
| `high_gate` | 0.60 | scaled index | heuristic | Start gate for viable-like | Not validated as assay criterion | Not tested |
| `low_gate` | 0.30 | scaled index | heuristic | Start gate for compromised-like | Not validated as assay criterion | Not tested |

### Report / stats

| Parameter | Value | Unit | Origin | Rationale / impact | Status | Sensitivity |
|---|---|---|---|---|---|---|
| `bootstrap_iterations` | 2000 | iteration | conventional | Bootstrap CI resamples | Partially validated | N/A |
| `seed` | 20260831 | int | fixed | Determinism of bootstrap | Partially validated | N/A |
| `save_meshes` | `true` | – | engineering | Persist PLY meshes | N/A | N/A |
| `max_meshes_in_preview` | 20 | count | engineering | Preview mesh cap | N/A | N/A |
| `stats.enabled` | `true` | – | engineering | Enable stats section | Partially validated | N/A |
| `stats.features` | `[volume_um3, sphericity]` | – | heuristic | Which features are summarized/compared | Partially validated | N/A |
| `stats.min_replicates_per_condition` | 3 | replicate | statistical | Minimum biological replicates per condition for `stats.py`'s LMM/OLS-fallback hypothesis test (docs/ALGORITHM_DECISIONS.md D11); a condition below this is silently excluded from that feature's test. Distinct from the bootstrap-CI threshold below (was previously mislabeled as gating the CI here). | Partially validated | N/A |
| `_MIN_REPLICATES_FOR_BOOTSTRAP_CI` (summary.py, not YAML-configurable) | 3 | replicate | statistical/heuristic | Minimum biological replicates with size data before `condition_summary.csv`'s bootstrap 95% CI is computed at all; coincidentally equal to, but a separate literal from, `stats.min_replicates_per_condition` above | Partially validated | N/A |

---

## Cellpose 3D segmentation parameters

Source of truth: `src/organoid_analysis/segmentation/cellpose_inference.py::SegmentationConfig` and `src/organoid_analysis/segmentation/parameter_estimation.py`. Exposed as UI controls in `src/organoid_analysis/web_interface/segmentation_workspace.py`. This route is independent of the classical `src/organoid_analysis/segmentation/watershed_instances.py` route catalogued above and was not previously catalogued here.

| Parameter | Value | Unit | Origin | Rationale / impact | Status | Sensitivity | User configurable |
|---|---|---|---|---|---|---|---|
| `model_type` | `cpdino-vitb` | – | heuristic (engineer-defined default) | Selects which Cellpose v4 foundation model runs; alternatives `cpsam_v2`, `cpdino`, `cpsam` trade accuracy for speed. Not independently benchmarked for this pipeline's organoid images. | Not validated | Not tested | Yes (UI dropdown) |
| `nuclei_diameter` | 30.0 | px (full resolution) | heuristic (engineer-defined default); can be auto-estimated per-image by `auto_config.estimate_diameter_from_stack` | Directly sets the scale Cellpose expects objects at; a wrong diameter is a common cause of over/under-segmentation. | Not validated on real data | Not tested | Yes (UI slider 5-200) |
| `cell_diameter` | 50.0 | px (full resolution) | heuristic (engineer-defined default); can be auto-estimated (defaults to the nucleus estimate) | Same role as `nuclei_diameter` for the cell/cytoplasm pass. | Not validated on real data | Not tested | Yes (UI slider 5-300) |
| `anisotropy` | 2.9 | dimensionless (Z step / XY pixel size) | heuristic fallback; preferentially read from OME/ImageJ TIFF metadata via `auto_config.auto_anisotropy` when available, never guessed from pixel content | Rescales the Z axis so Cellpose's 3D network sees near-isotropic voxels; a wrong value distorts 3D shape. | N/A — should be sourced from acquisition metadata for real data | Not tested | Yes |
| `xy_spacing_um` | 0.414 | µm/pixel | heuristic (default matching one specific microscope/objective configuration used during development) | Only informational for physical-scale display; segmentation itself is pixel-scale. | Not validated | N/A | Yes |
| `nuclei_flow_threshold` | 0.4 | dimensionless (Cellpose flow-error threshold) | conventional (Cellpose's own suggested operating range) | Higher values reject more flow-inconsistent masks (fewer, more confident objects). | Not validated on real data | Not tested | Yes (UI slider) |
| `cell_flow_threshold` | 0.6 | dimensionless | conventional | Same role for the cell pass. | Not validated on real data | Not tested | Yes (UI slider) |
| `nuclei_cellprob_threshold` | 0.0 | dimensionless (Cellpose logit threshold) | conventional (Cellpose default) | Foreground/background decision boundary in logit space. | Not validated | Not tested | Yes (UI slider -6..6) |
| `cell_cellprob_threshold` | 0.0 | dimensionless | conventional (Cellpose default) | Same role for the cell pass. | Not validated | Not tested | Yes |
| `flow3d_smooth` | 1.0 | dimensionless (Gaussian smoothing passed to Cellpose `flow3D_smooth`) | heuristic (engineer-defined default) | Smooths the estimated 3D flow field before instance construction; affects splitting of touching objects. | Not validated | Not tested | Yes (UI slider 0-5) |
| `xy_downsample` | 1.0 | fraction (0,1] | engineering (performance/memory trade-off, not scientific) | `1.0` = full resolution; lower values speed up inference at the cost of XY precision, then upsample masks with nearest-neighbor. | N/A (engineering) | N/A | Yes |
| `batch_size` | 8 | images/batch | engineering (GPU/MPS memory trade-off) | Does not change segmentation results, only throughput. | N/A (engineering) | N/A | Yes |
| `_ESTIMATE_SLICES` (auto_config.py) | 6 | Z-slices | heuristic | Number of Z-slices sampled for the fast diameter pre-estimate; more slices cost more time for marginal stability gain. | Not validated | Not tested | No (internal constant) |
| `_ESTIMATE_DOWNSAMPLE` (auto_config.py) | 0.5 | fraction | heuristic | Downsample factor for the fast 2D pre-segmentation pass used only to estimate diameter, not to segment. | Not validated | Not tested | No |
| `_DIAMETER_QUANTILE` (auto_config.py) | 0.5 (median) | quantile | heuristic | Which quantile of detected pre-segmentation object diameters is reported as "the" diameter estimate; median chosen for robustness to background-noise blobs. | Not validated | Not tested | No |
| plausible diameter filter (auto_config.py:142) | [5.0, 250.0] | px (downsample-corrected) | heuristic | Diameter pre-estimates outside this range are discarded as noise before taking the quantile. | Not validated | Not tested | No |
| quick pre-segmentation call (auto_config.py:116-119) | `flow_threshold=0.0, cellprob_threshold=0.0, min_size=5, batch_size=8` | mixed | heuristic (Cellpose defaults chosen for the diameter pre-estimate only, not the final segmentation) | Only affects the diameter *estimate* fed as a suggestion; does not affect the final segmentation parameters, which the user can override. | N/A (estimate only) | N/A | No |

**Model weight provenance:** `create_model` (`src/organoid_analysis/segmentation/cellpose_inference.py:136`) loads Cellpose's pretrained weights for `model_type` by name via the `cellpose` package; the specific weight file/version is whatever the installed `cellpose==4.2.1.1` package resolves (pinned in `pixi.lock`), not independently hashed or pinned by this repository. See `pyproject.toml` for the pinned `cellpose` package version.

---

## Cell–nucleus pairing parameters (`analysis cells` route)

Source of truth: `src/organoid_analysis/quantification/cellular_measurements.py`, exposed via the `cells` CLI subcommand in `src/organoid_analysis/command_line/organoid_commands.py`. This route is independent of `analyze-3d` (which uses maximum-overlap assignment, D1) and was not previously catalogued here. See `docs/ALGORITHM_DECISIONS.md` D10 for the bipartite-matching algorithm decision.

| Parameter | Value | Unit | Origin | Rationale / impact | Status | Sensitivity | User configurable |
|---|---|---|---|---|---|---|---|
| `min_cell_volume_um3` | 200.0 | µm³ | heuristic | Cells smaller than this are dropped before pairing (`too_small_cell`). | Not validated on real data | Not tested | Yes (`--min-cell-volume-um3`) |
| `min_nucleus_volume_um3` | 25.0 | µm³ | heuristic | Candidate nuclei smaller than this fail pairing (`too_small_nucleus`). | Not validated on real data | Not tested | Yes (`--min-nucleus-volume-um3`) |
| `max_nc_ratio` | 1.0 | fraction (nucleus voxels / cell voxels) | heuristic | Candidate pairs with nucleus volume exceeding the cell volume are rejected (`nc_ratio_too_high`) as biologically implausible. | Not validated | Not tested | Yes (`--max-nc-ratio`) |
| `min_nucleus_containment` | 0.5 | fraction (overlap voxels / nucleus voxels) | heuristic | A candidate nucleus must have at least this fraction of its own volume inside the cell to be paired (`nucleus_not_contained`). | Not validated | Not tested | Yes (`--min-nucleus-containment`) |
| `require_nucleus` | `True` | boolean | engineering/QC policy choice | When true, an unpaired cell is dropped entirely rather than kept anucleate; `--allow-nucleusless` flips this. | N/A (policy choice, not a measurement threshold) | N/A | Yes (`--allow-nucleusless`) |
| `neighbor_radius_um` | 25.0 | µm | heuristic | Physical radius for the cell-neighborhood density/nearest-neighbor query (`cell_neighborhood`); a cell whose search sphere would extend past the image border is flagged `neighborhood_complete=False` rather than silently biased. | Not validated | Not tested | Yes (`--neighbor-radius-um`) |

---

## Segmentation-validation matching parameter

Source of truth: `src/organoid_analysis/validation/segmentation_metrics.py::match_instances`. Stated as a narrative fact in `docs/METHODS.md` and `docs/SCIENTIFIC_SPEC.md` ("Hungarian assignment at IoU ≥ 0.5") but not previously entered in this catalog. See `docs/ALGORITHM_DECISIONS.md` D9.

| Parameter | Value | Unit | Origin | Rationale / impact | Status | Sensitivity | User configurable |
|---|---|---|---|---|---|---|---|
| `iou_threshold` | 0.5 | fraction (IoU) | conventional (common minimum-overlap convention in instance-segmentation benchmarks; not independently calibrated for this pipeline) | A predicted instance below this IoU with its best-matching truth object counts as a false positive/negative rather than a weak match; directly sets precision/recall/Dice/panoptic-quality values. | Not validated (conventional value, not calibrated) | Not tested | Yes (function argument; not currently exposed as a CLI flag) |

---

## Data-contract spacing/position tolerances

Source of truth: `src/organoid_analysis/microscopy_io/tiff_contract.py` (`ome_spacing`, `load_sample`) and `src/organoid_analysis/command_line/organoid_commands.py` (`_run_cells`, `_run_multilevel`). The same literal tolerance is repeated at `io.py:142,182,194,245` and `cli.py:91,99,156,163`.

| Parameter | Value | Unit | Origin | Rationale / impact | Status | Sensitivity | User configurable |
|---|---|---|---|---|---|---|---|
| spacing/Z-position consistency tolerance | `rtol=0.01, atol=1e-5` | relative fraction / µm | heuristic (engineering judgment for "close enough to be the same acquisition metadata", not independently calibrated) | Governs whether the pipeline accepts or rejects (a) OME plane Z-positions as a uniform grid, (b) explicit manifest/CLI spacing against OME metadata, and (c) spacing agreement between paired label volumes in the `cells`/`analyze-3d` routes. Too loose silently accepts mismatched acquisitions; too tight rejects valid metadata rounding. | Not validated | Not tested | No (hardcoded) |

---

## Statistical testing parameters (`stats.py`)

Source of truth: `src/organoid_analysis/statistics/inference.py`. See `docs/ALGORITHM_DECISIONS.md` D11 for the LMM/BH-FDR method decision.

| Parameter | Value | Unit | Origin | Rationale / impact | Status | Sensitivity | User configurable |
|---|---|---|---|---|---|---|---|
| `LOG10_FEATURES` | `{"volume_um3"}` | – | heuristic (engineering judgment: volume is right-skewed/multiplicative, sphericity is not) | Determines which features are log10-transformed before the LMM/OLS fit; changes the scale on which the omnibus test and pairwise contrasts are computed for that feature. | Not validated | Not tested | No (hardcoded set) |
| degenerate random-effect threshold | `1e-6` (× `fit.scale`) | dimensionless (variance ratio) | heuristic (numerical-stability engineering judgment) | Below this ratio the LMM's random-effect variance is treated as collapsed to zero and the fit falls back to clustered-SE OLS instead of trusting a numerically singular mixed model. | Not validated | Not tested | No (hardcoded) |

---

## Exploratory statistics/ML parameters (`src/organoid_analysis/statistics/exploration.py`, Tutorials 2-5)

Source of truth: `src/organoid_analysis/statistics/exploration.py`. These functions are explicitly scoped as reusable exploratory tutorial workflows (module docstring: "centralises the analysis workflows that live in Tutorials 2-5"), not part of the core `analyze`/`analyze-3d` measurement pipeline, and their outputs (classifier accuracy, feature importance, shape category) are exploratory/descriptive, not validated biological classifications.

| Parameter | Value | Unit | Origin | Rationale / impact | Status | Sensitivity | User configurable |
|---|---|---|---|---|---|---|---|
| Rod/Disk/Sphere gating percentiles | 75th / 25th percentile of the loaded dataset's prolate/oblate ratio | percentile of sample distribution | heuristic (descriptive-statistics split, not a biologically validated shape boundary) | Assigns a categorical shape label per object from a distribution-relative cutoff; the same object can receive a different label depending on what other objects are in the loaded dataset. | Not validated as a biological shape classification | Not tested | No |
| Binary-classifier hyperparameters (Tutorial 4) | `RandomForestClassifier(n_estimators=100, max_depth=10)`; `XGBClassifier(n_estimators=100, max_depth=5, learning_rate=0.1)` | – | engineering default (not tuned or cross-validated for this pipeline's data) | Determines the reported "best model" accuracy and feature-importance ranking shown to the user; untuned defaults may understate achievable accuracy or misrank features. | Not validated | Not tested | No |
| Cluster-classifier hyperparameters | `RandomForestClassifier(n_estimators=200, max_depth=15)` | – | engineering default (not tuned) | Determines cluster-membership prediction accuracy and the feature-importance ranking used to characterize clusters. | Not validated | Not tested | No |

---

## Notes

- All values have a machine-readable home in `src/organoid_analysis/config.py` or `src/organoid_analysis/quantification/multilevel_relationships/config.py`; no silent magic numbers gate scientific behavior.
- Thresholds marked `Not validated on real data` require independent calibration to justify a biological claim. The pipeline labels such results as research-use only and does not claim clinical or biological endpoint validity.
- Viability gates (0.30/0.60) are documented as starting gates, not experimentally calibrated assay thresholds.

## Traceability

- Multilevel config: `src/organoid_analysis/quantification/multilevel_relationships/config.py`
- Classical config defaults + validation: `src/organoid_analysis/config.py`
- CLI exposure: `src/organoid_analysis/command_line/organoid_commands.py`
- Cellpose 3D segmentation: `src/organoid_analysis/segmentation/cellpose_inference.py`, `src/organoid_analysis/segmentation/parameter_estimation.py`, UI exposure in `src/organoid_analysis/web_interface/segmentation_workspace.py`
- Cell–nucleus pairing: `src/organoid_analysis/quantification/cellular_measurements.py`, CLI exposure in `src/organoid_analysis/command_line/organoid_commands.py::_run_cells`
- Segmentation-validation matching: `src/organoid_analysis/validation/segmentation_metrics.py`
- Data-contract tolerances: `src/organoid_analysis/microscopy_io/tiff_contract.py`, `src/organoid_analysis/command_line/organoid_commands.py`
- Statistical testing: `src/organoid_analysis/statistics/inference.py`
- Exploratory tutorial statistics/ML: `src/organoid_analysis/statistics/exploration.py`
