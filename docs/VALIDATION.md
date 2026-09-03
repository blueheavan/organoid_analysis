# Validation record — 2026-08-31

This record distinguishes software verification from biological validation. No real microscopy accuracy, patient-derived viability performance, clinical response prediction, or macOS installation performance has been established.

## Executed checks

- **37 automated tests passed**, including end-to-end nonempty, empty, and partial-failure runs. The final test run took approximately 7.85 seconds in the provided runtime; this is not a benchmark for the user's Mac. Machine-readable results: `test_results.xml`.
- The final seeded demo processed **15 synthetic three-channel Z-stacks and 90 organoids** without failed fields. Analysis took approximately 33.68 seconds in the same runtime, excluding environment setup and data generation.
- Ground truth was evaluated using one-to-one Hungarian instance assignment at IoU ≥0.5. All 90 objects matched; false positives=0 and false negatives=0. Median matched Dice=1.0000; minimum Dice≈0.99965. Median absolute relative volume error=0; maximum≈0.000694 (0.0694%).
- Synthetic marker-state agreement was 100% on matched objects. This is **not a cross-validated classifier result**: no machine-learning classifier was trained, the phantoms are easy, and the simulator and control-scaling rules deliberately share simplified marker assumptions.
- Treatment and calibration table counts were reconciled. The per-object output has 90 rows and 51 columns; required volume, area and sphericity values were nonmissing. Every treatment group had 3 synthetic biological replicates and 18 measured objects.
- PNG treatment charts and representative orthogonal/3D QC panels were visually inspected. SVG and PDF versions were exported from the same Matplotlib figures. The HTML report embeds its images and does not require an external plotting CDN. Native 3D masks retain physical spacing metadata.
- `pyproject.toml` and `pixi.toml` parse successfully. Existing-output refusal, TIFF axes, physical-spacing conflict detection, nonuniform Z positions, 12-bit saturation, missing controls and biological-replicate weighting were exercised.

The excellent demo overlap is a consequence of simple, high-contrast ellipsoids, not evidence that this baseline will perform similarly on hollow, budding, dense, noisy or brightfield organoids. Dice alone also does not establish surface-area accuracy: unsmooth voxel surfaces have a known discretization bias, explicitly checked with analytical sphere phantoms using a finite tolerance. Full-resolution surface measurements are not replaced by visually smoother preview meshes.

## Validation environment

| Component | Version |
|---|---|
| Python | 3.12.13 |
| Platform | Linux x86_64, glibc 2.39 |
| NumPy | 2.5.2 |
| SciPy | 1.18.1 |
| scikit-image | 0.26.0 |
| tifffile | 2026.8.23 |
| pandas | 2.2.3 |
| Matplotlib | 3.10.8 |
| PyYAML | 6.0.3 |
| pytest | 9.1.1 |

The run emitted 64 upstream deprecation warnings because scikit-image currently assigns NumPy array shapes in places deprecated by NumPy 2.5. No tests failed; warnings were not suppressed. This should be revisited before upgrading those dependencies further. The table above records the tested direct Python packages for this validation run (Python 3.12.13 / Linux x86_64, i.e. not the target osx-arm64 pixi environment); it is not itself a lock file. `pixi.lock` (committed at the repository root, generated and pinned on osx-arm64) is the actual fully resolved, cross-run-reproducible lock for the pixi environment used elsewhere in this repository (`docs/VALIDATION_REPORT.md`'s Gate 11 baseline). There is no separate `requirements-tested.txt` file in this repository; do not rely on such a file existing.

## What still requires real data

Boundary segmentation needs independent annotations across treatments, organoid sizes and imaging batches. Review missed/extra instances, split/merge errors, surface bias and treatment-dependent detection. Confirm microscope voxel spacing, uniform Z steps, optical resolution, channel registration, exposure, bit depth and stain penetration. Calibrate viability gates using matched controls and biological reference measurements. A morphology/brightfield viability predictor additionally requires independent labeled training/test data partitioned at the relevant patient/culture/batch level. Statistical tests require the actual experimental hierarchy and pairing. None of these validations is replaced by the synthetic test set.
