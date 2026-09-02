# Parameters Catalog — Organoid Pipeline

Version: 1.0.0
Date: 2026-09-02

This catalog records the provenance of every scientifically consequential parameter, threshold, and default in the pipeline. Properties follow the rigor protocol (section 6). Engineering-only parameters that cannot affect scientific results (e.g. output verbosity, file naming) are not catalogued.

Legend for status: `Validated` / `Partially validated` / `Not validated`; sensitivity column: `Tested` / `Not tested` / `N/A`.

---

## Multilevel 3D analysis parameters

Source of truth: `src/analysis/multilevel3d/config.py` (dataclass) and CLI flags in `src/analysis/cli.py`.

| Parameter | Value | Unit | Origin | Rationale | Scientific impact | Calibration artifact | Status | Sensitivity | User configurable |
|---|---|---|---|---|---|---|---|---|---|
| `minimum_voxels` | 5 | voxel | heuristic (engineer-defined) | Objects with <5 voxels are too small to be reliably measured / non-degenerate | Flags too-small objects as `too_small` QC; affects QC only, not measurement set | N/A (validation only) | Partially validated | Not tested | Yes (`--minimum-voxels`) |
| `low_parent_overlap_fraction` | 0.5 | fraction (of child voxels) | heuristic | <50% of a child's voxels overlapping its assigned parent suggests a weak/ambiguous assignment | Flags `low_parent_overlap` QC | N/A | Partially validated | Not tested | Yes |
| `mad_z_threshold` | 3.5 | robust z-units | heuristic | ~3.5 σ beyond median/MAD identifies extremes without catastrophic masking | Flags `volume_outlier` QC | N/A | Partially validated | Not tested | Yes |
| `core_max_normalized_radial_position` | 0.5 | dimensionless | heuristic | Equivalent-radius ≤0.5 = core | Defines `core_cell_fraction` | N/A | Not validated (operational definition) | Not tested | Yes |
| `peripheral_min_normalized_radial_position` | 0.8 | dimensionless | heuristic | Equivalent-radius ≥0.8 = periphery | Defines `peripheral_cell_fraction` | N/A | Not validated (operational definition) | Not tested | Yes |

---

## Classical analysis config defaults

Source of truth: `src/analysis/config.py` (module `DEFAULTS`), validated by `validate_config`.

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
| `stats.min_replicates_per_condition` | 3 | replicate | statistical | Bootstrap CI requires ≥3 size-evaluable replicates | Partially validated | N/A |

---

## Notes

- All values have a machine-readable home in `src/analysis/config.py` or `src/analysis/multilevel3d/config.py`; no silent magic numbers gate scientific behavior.
- Thresholds marked `Not validated on real data` require independent calibration to justify a biological claim. The pipeline labels such results as research-use only and does not claim clinical or biological endpoint validity.
- Viability gates (0.30/0.60) are documented as starting gates, not experimentally calibrated assay thresholds.

## Traceability

- Multilevel config: `src/analysis/multilevel3d/config.py`
- Classical config defaults + validation: `src/analysis/config.py`
- CLI exposure: `src/analysis/cli.py`
