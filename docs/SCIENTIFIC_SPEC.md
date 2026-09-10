# Scientific Specification — Organoid Pipeline

Version: 1.0.0
Date: 2026-09-02
Classification: S3 (Inferential) for the pipeline as a whole, because `src/organoid_analysis/statistics/inference.py`'s condition-comparison hypothesis testing (linear mixed-effects model + Benjamini-Hochberg FDR-corrected pairwise contrasts; `cfg["stats"]["enabled"]` defaults to `true`) is exactly the "hypothesis testing, differential analysis" example category the S3 rubric names. Core measurement/QC (morphology, hierarchy, topology, viability-state gating) is S2 (Analytical) on its own; see docs/ALGORITHM_DECISIONS.md D11 for the S3 component's method, validation status, and known small-sample-inference limitation. The bundled exploratory UI tutorial workflows (`src/organoid_analysis/statistics/exploration.py`: trained classifiers, clustering) are also S3-shaped but are explicitly out of scope of this specification's measurement claims (see docs/PARAMETERS.md, "Exploratory statistics/ML parameters").

## 1. Scientific Objective and Question

**Primary objective:** Quantitative 3D morphological analysis of organoid structures from microscopy Z-stack images.

**Scientific question:** How do organoid size, shape, internal organization, and cellular composition vary across experimental conditions?

**Scope:** This software provides morphological measurements and descriptive statistics, plus, by default, exploratory cross-condition hypothesis testing (linear mixed-effects model omnibus test and BH-FDR-corrected pairwise contrasts on the two summary features in `stats.features`; see docs/ALGORITHM_DECISIONS.md D11). It does not perform causal inference or predictive modeling as part of the core measurement pipeline. The statistical tests it does perform are exploratory (small documented replicate minimums, an approximate rather than full Satterthwaite/Kenward-Roger small-sample correction, no predeclared experimental design) and must not be treated as confirmatory without independent statistical review appropriate to the actual study design. Biological interpretation remains with the researcher.

## 2. Intended Users and Permitted Research Use

**Intended users:**
- Research scientists in biology, pharmacology, and related fields
- Lab technicians performing organoid assays
- Computational biologists analyzing organoid data

**Permitted research use:**
- Characterizing organoid morphology in controlled experiments
- Comparing morphological features across treatment conditions
- Quality control of organoid cultures
- Method development for organoid analysis

**Explicitly NOT permitted:**
- Clinical diagnosis or patient stratification
- Treatment selection or go/no-go decisions without independent validation
- Regulatory submissions without separate regulatory assessment
- Any use where morphological measurements directly drive clinical decisions

## 3. Supported and Unsupported Input Domain

### Supported inputs
- **Image format:** 3D TIFF or OME-TIFF Z-stacks
- **Axes:** CZYX (canonized internally)
- **Spacing:** Uniform Z-step with known physical spacing (µm/voxel)
- **Channels:** Structural fluorescence (required), Calcein/PI viability dyes (optional)
- **Bit depth:** 8-bit, 12-bit, 16-bit integer; 32-bit float
- **Organoid type:** Spherical or near-spherical organoids with clear boundaries
- **Size range:** 2,000 µm³ to ~10⁶ µm³ (configurable)

### Unsupported inputs
- 2D images or single optical sections
- Non-uniform Z-spacing
- Multi-well plate images without per-field processing
- Organoids with complex internal architecture (lumen, branching)
- Brightfield images without probability maps (exploratory only)
- Tissue samples, organoids >1 mm diameter, or non-organoid biological structures
- Images with motion artifacts, significant drift, or poor registration

## 4. Sample/Specimen/Organism/Assay Platform Context

**Biological system:** 3D organoid cultures derived from:
- Stem cells (iPSC, ESC)
- Primary tissue
- Cell lines
- Patient-derived samples

**Assay platforms:**
- Confocal microscopy (preferred)
- Light-sheet microscopy
- Two-photon microscopy
- Spinning disk confocal

**Staining:**
- Nuclear stain (Hoechst, DAPI) — required for Cellpose segmentation
- Viability dyes (Calcein AM, PI) — optional for state classification
- Structural markers — optional for feature extraction

**Acquisition requirements:**
- Isotropic or near-isotropic XY resolution (≤1 µm/pixel recommended)
- Z-step ≤2 µm for reliable 3D reconstruction
- Minimum 5 Z-slices per organoid
- No saturation in intensity channels

## 5. Scientific/Experimental Unit

**Primary unit:** Individual organoid instance

**Secondary units:**
- Individual cell instance (within organoid)
- Individual nucleus instance (within cell)
- Well (group of organoids from same culture well)
- Biological replicate (independent culture preparation)
- Condition (experimental treatment group)

**Independence structure:**
- Organoids within a well are pseudo-replicates (shared culture conditions)
- Wells within a batch are technical replicates
- Biological replicates are independent culture preparations
- Statistical inference should use biological replicate as the unit of analysis

## 6. Biological vs Technical Replicate Structure

**Biological replicates:** Independent culture preparations from the same or different donors. These capture biological variability and should be used for inference.

**Technical replicates:** Multiple images from the same well or culture. These capture measurement variability but not biological variability.

**Common mistake:** Treating organoids or technical replicates as biological replicates (pseudoreplication). The software outputs per-organoid measurements and its classical workflow can perform exploratory mixed-effects tests, but those tests depend on correctly supplied well, batch, and biological-replicate identifiers. Researchers must verify the nested design and must not interpret object-level rows as independent biological replicates.

## 7. Output Semantics

### Direct measurements (from label masks)
- Volume (µm³)
- Surface area (µm²)
- Sphericity (dimensionless, π^(1/3)(6V)^(2/3)/A)
- Principal axes lengths (µm)
- Elongation ratios
- Surface-to-volume ratio
- Centroid position (physical coordinates)

### Derived measurements
- Nucleus-to-cell volume ratio
- Cell-nucleus count
- Radial position (normalized to organoid equivalent radius)
- Distance to organoid centroid/surface (µm)
- Contact area between adjacent cells (µm²)
- Cell degree (number of neighbors)

### Intensity measurements (when raw channel provided)
- Mean/min/max nuclear intensity
- Integrated nuclear intensity
- Coefficient of variation (CV) of chromatin intensity

### Classification outputs
- Viability state (viable-like, compromised-like, indeterminate)
- QC flags (border touch, fragmentation, size outliers, etc.)

### Statistical summaries
- Well-level medians
- Condition-level means and bootstrap CIs
- State fractions

## 8. Downstream Consequence

**Immediate downstream use:**
- Feature tables for statistical analysis
- Quality control reports
- Visualization and exploration

**Indirect downstream use (researcher responsibility):**
- Statistical testing across conditions
- Dose-response modeling
- Correlation with other assays
- Biological interpretation and reporting

**Consequence of error:**
- Morphological mischaracterization could lead to:
  - Incorrect biological conclusions
  - Misleading treatment effects
  - Wasted resources on follow-up experiments

**Mitigation:**
- Explicit QC flags for potentially unreliable measurements
- Transparent reporting of measurement uncertainty
- Clear documentation of assumptions and limitations

## 9. Acceptable Reference Standard

### For segmentation accuracy
- **Reference:** Manually annotated 3D instance masks by domain expert
- **Independence:** Annotations must be independent of algorithm development
- **Matching criterion:** Hungarian assignment at IoU ≥ 0.5
- **Metrics:** Precision, recall, Dice, IoU, panoptic quality

### For morphological measurements
- **Reference:** Analytical phantoms with known geometry (spheres, ellipsoids)
- **Validation:** Compare measured volume/area to analytical values
- **Tolerance:** Volume error <1%, area error <5% for analytical shapes

### For biological validity
- **Reference:** Independent biological measurements (e.g., dry weight, cell count)
- **Status:** NOT ASSESSED — requires separate biological validation study
- **Limitation:** Synthetic validation does not establish biological accuracy

## 10. Explicit Limitations and Unresolved Assumptions

### Known limitations
1. **Segmentation accuracy:** Depends on image quality, staining, and organoid density. No claim of general accuracy without independent validation.
2. **Envelope geometry:** Measures outer envelope, not internal structure. Enclosed lumens are filled by default.
3. **Sphericity bias:** Voxel discretization can produce values >1.0. Values >1.05 are flagged but not corrected.
4. **Intensity measurements:** Raw intensities are not normalized across experiments. Cross-experiment comparison requires separate normalization.
5. **Viability states:** Calcein/PI states describe signal patterns, not cell viability fraction. States are not validated against independent viability assays.
6. **Statistical summaries:** Bootstrap CIs (well/condition-level descriptive summaries) are exploratory with small samples. The default-on hypothesis tests (linear mixed-effects omnibus test, BH-FDR-corrected pairwise contrasts; see §1 and docs/ALGORITHM_DECISIONS.md D11) are likewise exploratory — small documented replicate minimums, an approximate rather than full Satterthwaite/Kenward-Roger small-sample correction, no predeclared design — and must not be treated as confirmatory.
7. **Spatial features:** Contact area is voxel-face based. Centroid distance and kNN are not used for contact definition.

### Unresolved assumptions
1. **Isotropic voxels:** Algorithm assumes near-isotropic spacing. Highly anisotropic data may produce biased measurements.
2. **Single time point:** No temporal tracking or dynamic analysis.
3. **Homogeneous staining:** Assumes uniform dye penetration. Edge effects or incomplete staining are not modeled.
4. **Independent organoids:** Assumes organoids do not touch or overlap. Overlapping organoids may be merged or split incorrectly.

### What this software does NOT claim
- Clinical accuracy or diagnostic utility
- Biological validity without independent validation
- Generalizability to all organoid types or imaging conditions
- Replacement for expert manual analysis
- That its default-on exploratory significance tests (docs/ALGORITHM_DECISIONS.md D11) are confirmatory, independently reviewed, or a substitute for statistical analysis matched to a predeclared experimental design

## 11. Revision History

| Version | Date | Changes |
|---------|------|---------|
| 1.0.0 | 2026-09-02 | Initial scientific specification |
| 1.1.0 | 2026-09-03 | Reclassified S2→S3: `stats.py`'s default-on hypothesis testing (D11) was shipped without updating this document's classification/scope/"does NOT claim" sections; corrected per an independent scientific-software audit. |
| 1.1.1 | 2026-09-08 | §10 item 6 ("Statistical summaries") still read "No inferential tests are provided," left over from before the 1.1.0 S2→S3 reclassification and directly contradicting §1's Scope (which correctly describes the default-on LMM/BH-FDR hypothesis tests). Corrected to describe those tests' exploratory status consistently with §1. |

## 12. References

1. scikit-image measurement API: https://scikit-image.org/docs/stable/api/skimage.measure.html
2. Cellpose: https://cellpose.readthedocs.io/
3. Thermo Fisher viability dyes: https://www.thermofisher.com/uk/en/home/life-science/cell-analysis/cell-viability-and-regulation/cell-viability.html

## 13. Current implementation audit — 2026-09-10 (governing qualification)

The preceding intended research use is retained, **not validated by its presence in this document**. Species/platform, ≤1 µm XY, ≤2 µm Z, minimum five slices, shape and size recommendations above lack project calibration/held-out evidence; treat them as claimed operating assumptions with INSUFFICIENT EVIDENCE, not established acquisition specifications. The numerical requirements in §9 (<1% volume, <5% area) are retained. Surface accuracy currently fails that requirement; see [current validation report](VALIDATION_REPORT.md). Neither defaults nor previous PASS labels establish suitability.

### Reconstructed end-to-end paths and scientific impact

| Stage / scientific impact | Actual implementation and contract | Required validation |
|---|---|---|
| Read (S1–S2: identity and units) | `microscopy_io/tiff_contract.py` reads manifest-driven CZYX with ≥3 planes, OME or explicit calibration. `zstack_reader.py` reads first series into ZYX/CZYX with ImageJ/TIFF/OME calibration; time selection explicit in API, multi-time upload rejected by Web. Cellpose legacy QYX/IYX fallback assumes a grayscale Z-stack and leaves spacing unknown. | Exact identity/units/time/channel tests; acquisition grid qualification; corrupted/unknown metadata rejection. |
| Inspect (S2 if it guides selection) | Web preview uses a separate display payload with percentile contrast and possible uint8/downsampling. Native VTK, browser vtk.js and 2D orthogonal views are distinct rendering paths. | Axis/spacing and raw-data preservation tests; actual renderer checks; browser upload and visual performance assessed separately. |
| Segment (S2) | Classical structural watershed/probability/imported labels and Web Cellpose cell/nucleus masks are separate routes. Cellpose alone does not produce a qualified organoid envelope. | Qualified independent 3D annotations, object-level detection/overlap and acquisition strata; parameter calibration and held-out sensitivity. |
| Quantify (S2) | Classical geometry measures filled organoid envelopes and eligibility-filtered summaries. Web object features measure filled size/axes but raw centroid/solidity and no classical eligibility exclusion. `analyze-3d` accepts registered organoid/cell/nucleus labels and measures raw voxels, with review flags but no automatic aggregate exclusion. | Analytic volumes/moments/contact; surface accuracy; raw/envelope semantics; intensity/registration/saturation qualification. |
| Infer/explore (S3) | `statistics/inference.py` fits object-level LMM or clustered OLS; `aggregation.py` reports equal-well replicate summaries; `exploration.py` supplies row-level tests/classifiers/clustering. | Correct independent units/design, residual assumptions, type-I error/CI coverage, multiplicity family, group-aware validation. |
| Export (S1–S2) | CLI CSV/JSON/OME; multilevel Parquet/JSON/OME; Cellpose ZIP/OME/summary/provenance. Missing data are not evidence of zero effects or successful tests. | ID/unit round trips, complete/incomplete outputs, source/config/model/input traceability. |

### Study-unit and output qualification

An organoid is the **measurement unit**, not necessarily the independent experimental unit. Wells are not inherently technical or biological replicates: their independence depends on randomization, donor/preparation and treatment assignment, which must be documented outside these data tables. For the shipped classical inference, `biological_replicate` identifies the same replicate across batches within a condition; the code treats conditions as independent groups. Paired donor experiments, cross-condition repeated measures, nested well effects and temporal inference are unsupported without a justified design. This audit did not replace that method or invent study metadata.

Coordinates are array-local voxel-center coordinates in micrometres, not microscope stage/world coordinates. Equal shape and spacing are necessary but insufficient to prove registration or matching origins/directions. Descending monotonic Z acquisition is accepted for shape measurements; anatomical orientation is not preserved as a world-coordinate transform. EDT depth is distance to a background voxel center at the rounded child centroid, not a continuous-surface distance.

The earlier list of supported stains does not establish that a nuclear channel outlines whole organoids. Whole-organoid segmentation needs a qualified structural boundary signal or imported annotated/model masks. Without that evidence, a successful pipeline may count substructures instead of organoids. Local images have no independent organoid/cell/nucleus annotation set, documented experimental design or orthogonal viability assay. No biological accuracy, live-cell percentage, calibrated classifier or confirmatory treatment effect is supported by the current validation.

### Reference standard and uncertainty qualification

Analytical voxel solids and hand-counted label relationships qualify numerical arithmetic only. The demo generator is project-owned and shares signal assumptions with viability classification; it is not independent assay validation. The local real files are identified by hashes, but acquisition provenance, stain identity, annotator protocol, blinding, adjudication, inter-annotator agreement and biological sample independence are unavailable. The model smoke test measures executability on one crop, not accuracy.

For biological segmentation/viability/statistical claims, required N, precision target and operating strata must come from the actual experiment and qualified reference uncertainty. No defensible universal N or error target can be inferred here; status is INSUFFICIENT EVIDENCE. Clinical/decision-critical intended use is absent; S4 decision-risk and regulatory conformity assessment are NOT APPLICABLE to this audit.
