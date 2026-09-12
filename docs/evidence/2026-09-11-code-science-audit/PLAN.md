# Current code scientific audit plan

Frozen before new runtime probes on 2026-09-11. Scope: audit only; no production
algorithm, parameter, acceptance criterion, or test-suite changes.
Source baseline: 3639bec910337ec623782b657fe62cf1dafc3800, clean main.

Intended use is reconstructed from current source and SCIENTIFIC_SPEC.md:
registered 3D microscopy label/intensity measurements and exploratory inference.
Organoids/cells/nuclei are measurement units; independent culture/donor units
must be supplied for biological inference. Highest impact is S3 because
statistics/inference.py runs hypothesis tests by default and exploration.py
provides classifiers and cluster analyses. This audit is Tier D, LIMITED
INDEPENDENCE; independently calculated oracles do not make the reviewer independent.

| ID | Requirement / predeclared expected behavior | Oracle and cases | Evidence limit |
|---|---|---|---|
| C1 | Volume error <1%, area error <5% per existing specification; scale covariance holds | Voxelized spheres at radii 10, 20, 40 um, spacings (1,1,1), (2,1,1); analytic 4*pi*r^3/3 and 4*pi*r^2 | Deterministic geometry verification; no biological N or accuracy claim |
| C2 | Known physical calibration conflicts reject even when another axis is missing | Hand-authored OME TIFF with X/Y=1 um and missing Z; manifest X/Y=2 um; compare two partial grids with conflicting X/Y | Exact metadata identity, not image registration validation |
| C3 | Recommended k and displayed silhouette must describe the same feature transformation as fitted KMeans | Seed 20260911; change a feature's measurement unit by 1000; independently compute standardized silhouette | Metamorphic consistency; no proof of biological clusters |
| C4 | Parent IDs are preserved exactly, including unmatched children | Int64 organoid ID 2^53+1, assigned and unassigned nuclei; hand-calculated parent identity and counts | Supported integer-domain boundary, not normal ID prevalence |
| C5 | Impossible geometric metrics have visible object QC in each route | Interior 2x2x2 label, minimum_voxels default 5; compare sphericity to existing 1.05 review limit | Numerical warning consistency, not new threshold calibration |
| C6 | Failed segmentation/QC states remain visible downstream | Trace sample/object flags, aggregation and exported artifacts; controlled disagreement if warranted | Review flags do not automatically require exclusion; distinguish policy from implementation defect |
| C7 | Engineering baseline remains separate from science acceptance | Locked Pixi default non-render suite, lint, typecheck debt gate, science gate; capture exact exits | No rendering, build or model-accuracy claim from these checks |

Statistical study validation, qualified real annotations, viability assay
validation, biological sample-size/precision targets and held-out parameter
calibration are assessed by inventory; missing evidence is INSUFFICIENT EVIDENCE.
Do not invent acceptance targets or infer biological accuracy from synthetic
fixtures, previous reports, existing tests, or gate PASS declarations.

Additional source-confirmed defects may receive analytical counterexamples;
record any added probe and its expected invariant before executing it.

## Extension before follow-up execution

C8 (2026-09-11): the specification lists 12-bit data as supported, while
features.marker_measurements infers detector saturation from uint16 storage.
Use a known 12-bit clipping ceiling of 4095, stored as uint16, with Calcein
ROI entirely at 4095 and unsaturated PI/background. Compare default limits
against the explicit, independently known 4095 ceiling. Expected: default
configuration must not certify this ROI as unsaturated/eligible; unknown
acquisition saturation needs qualification. This checks QC semantics, not
the biological meaning of a viability state. No thresholds are changed.
