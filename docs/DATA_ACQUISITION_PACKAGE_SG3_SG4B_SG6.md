# What the remaining gate items need from the laboratory

SG-3A, SG-3B, SG-4B, SG-6A and SG-6B cannot be passed by changing software.
Each requires physical specimens, blinded human annotation or traceable
instrument standards that do not exist in this repository. Their acceptance
criteria and sample sizes were frozen on 2026-09-16 (D-5 to D-13 in
[OWNER_DECISIONS.md](OWNER_DECISIONS.md)); this page turns them into an
acquisition checklist and names the tooling that already exists to analyse the
data once it is collected.

## Public data checked (2026-09-26)

A PubMed search for independently annotated 3D organoid segmentation data found
no set that meets the frozen SG-3 protocol (blinded two-annotator 3D instance
masks, the intended modalities and voxel sizes, no overlap with Cellpose
training data, ≥ 47 biological samples per stratum):

| Dataset | Why it does not qualify |
|---|---|
| OrgaQuant, Kassis et al. 2019, [doi:10.1038/s41598-019-48874-y](https://doi.org/10.1038/s41598-019-48874-y) | 2D brightfield, bounding boxes, not 3D instance masks |
| COSMOS nuclear training set, Mukashyaka et al. 2023, [doi:10.1101/2023.03.03.531019](https://doi.org/10.1101/2023.03.03.531019) | nuclei only (3,862 cells in 36 organoids, one laboratory, 5 µm Z); annotator independence and agreement not reported in the abstract |
| 3D-Organoid-SwinNet, Sohaib et al. 2025, [doi:10.1109/JBHI.2024.3511422](https://doi.org/10.1109/JBHI.2024.3511422) | semantic (not instance) masks, breast-cancer subtypes only |

Such sets could support an *external characterization*, never a qualification
of this project's intended use.

## SG-6A / SG-6B — acquisition calibration (do first; no biology required)

| Step | Requirement (D-12, D-13; `evidence/2026-09-13-volume-audit/STUDY_PROPOSALS_SG3_SG6.md` §2) |
|---|---|
| Standards | certified stage micrometer or grid slide (XY); certified axial standard or beads at certified separation in the **same mounting medium** as organoids (Z); sub-resolution multi-spectral beads (registration) |
| Configurations | every objective × immersion × Z-step × camera combination used for organoids |
| Sessions | ≥ 3 independent sessions per configuration, fields at centre and edges |
| Acceptance (T1 stability) | lateral ≤ 0.10 %, axial ratio ≤ 0.25 % (upper 95 % bound) |
| Acceptance (T2 accuracy) | lateral ≤ 0.15 %; axial reported as a stated systematic uncertainty |
| Tooling | `validation/spacing_calibration.py` (`ratio_report`, `stability_verdict`) |
| Re-qualification | every 6 months and after any objective, medium, Z-step, camera or software change |

## SG-3A / SG-3B — segmentation accuracy

| Step | Requirement (D-5, D-6; `evidence/2026-09-11-measurement-vv/SEGMENTATION_VALIDATION_PROTOCOL.md`) |
|---|---|
| Specimens | ≥ 47 independent biological samples (donors or lines) **per stratum** (brightfield; membrane fluorescence), ≈ 20 objects each; strata cover voxel size, density, lumen/hollow, border truncation |
| Annotation | manual 3D instance masks at native resolution under a written SOP, blinded to model output and condition; two annotators on ≥ 20 % of fields; third annotator adjudicates |
| Independence | calibration split and test split disjoint by biological sample; test split used once; images not in Cellpose training data |
| Acceptance | C1 recall and precision lower 95 % bound ≥ 0.90; C2–C5 bias bounds relative to δ = 0.20 (volume, area) / 0.10 (diameter, sphericity); 95 % CI of the between-condition difference in median relative bias within ±δ/8 (2.5 %) (D-6) |
| Tooling | `validation/segmentation_metrics.py` (Hungarian IoU matching, P/R/F1, Dice, PQ, split/merge); the report's segmentation-validation section reads `truth_labels_path` |

## SG-4B — Calcein/PI biological validity

| Step | Requirement (D-7, D-8) |
|---|---|
| Reference | blinded human classification of the same objects (viable / mixed / compromised / indeterminate) by two annotators on a predeclared stratified subsample; ATP assay as secondary well-level reference |
| Controls | live and dead controls in every batch; a graded-insult series for monotonicity and dynamic range; held-out batches |
| Acceptance | per-class sensitivity and precision lower 95 % bound ≥ 0.80; macro-F1 lower bound ≥ 0.80; indeterminate ≤ 20 % per condition and between-condition difference ≤ 5 pp |
| Tooling | `validation/agreement_metrics.py` (confusion matrix, Wilson intervals, macro-F1, weighted kappa, inter-rater agreement); the four-state rule itself is analytically qualified (SG-4A) |

## Order of work

1. SG-6 standards (days; no biology) — unblocks physical-unit claims and the
   SG-3 spacing-ratio dependency.
2. One shared annotation campaign serves SG-3 and SG-4B (D-8): annotate
   organoid, cell and nucleus masks and the viability class on the same objects.
3. Each study: freeze its protocol and seed in a committed directory **before**
   data are opened, then generate a canonical record, as was done for SG-1a,
   SG-2a, SG-4A and SG-5.
