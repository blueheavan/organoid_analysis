# Workflow optimization contract — 2026-09-11

Baseline: clean `657204d8d8be2e547063272fd139161f59dc62be`; Pixi lock SHA256
`fdba57d5d030dc0c72eaae1e3f1d950f4b58c44f1ab471aa54ce193308f9ae8c`.
This plan precedes implementation and final validation. Scope is S2:
Cellpose resampling/calibration, preview geometry, and measured mask features
and their export. Existing S3 statistical methods are outside this change.

## Intended behavior and method decisions

- The Web object is a nucleus or cell label, selected explicitly. Its default
  measurements use raw label voxels, consistent with multilevel measurement.
  Optional hole filling measures an envelope; volume, centroid, axes, surface
  and solidity must all use that same geometry. Export the geometry definition,
  raw voxel count, raw volume and QC with full numeric precision. Historical
  filled-envelope results can be requested explicitly; no silent QC exclusion.
- Border contact and face-disconnected fragments are review flags, not evidence
  of biological abnormality. No new biological threshold or size filter.
- Retain SciPy's existing center-grid interpolation (`grid_mode=False`). For
  input length N and actual output length M, spacing changes by (N-1)/(M-1).
  Cellpose diameter and Z/XY anisotropy use the reciprocal effective XY scale.
  Its scalar XY contract cannot represent unequal effective X/Y spacings:
  reject such downsampling before model execution and offer full resolution.
  Nearest-neighbor labels return to the original grid. Fractional interpolated
  intensities must not be rounded to an integer input dtype.
- Preview can represent separate X/Y spacings and must preserve center-coordinate
  extents and source arrays. Honor the array payload budget even below 0.1 scale;
  if two XY samples per axis cannot fit, raise a clear size error.
- Feature downloads include both a CSV and JSON provenance (selected mask,
  physical spacing and origin of calibration, geometry definition, model/run
  identity where available, source and content hashes). One field's rows remain
  descriptive observations, not independent replicates or inferential evidence.

Sources: [SciPy zoom](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.zoom.html)
defines center-grid coordinates and output dtype; [Cellpose 3D](https://cellpose.readthedocs.io/en/latest/do3d.html)
defines anisotropy as Z/XY spacing. These support API/math semantics, not
project-data segmentation accuracy. Basis ESTABLISHED for coordinate arithmetic,
CONTEXT_DEPENDENT for raw versus envelope estimands; project suitability PARTIAL.

## Frozen acceptance

| ID | Reference / fixture / unit | Acceptance |
|---|---|---|
| O1 | Hand-counted asymmetric hollow label, anisotropic voxel spacing | Raw volume/count/centroid exact to 1e-12; envelope metrics use filled support; label IDs unchanged. |
| O2 | Voxel volumes below 1e-4 µm³ and CSV round trip | Positive values preserved; relative error <=1e-12; no computation-time decimal rounding; empty table retains its schema. |
| O3 | Border and face-disconnected labels | Correct explicit flags, no automatic row removal; raw/envelope basis recorded. |
| O4 | Coordinate ramps on even/odd square and asymmetric rectangular grids | Original center extents preserved to 1e-12; model uses actual scalar scale; unequal scales reject before model invocation; interpolated ramp values within 1e-6. |
| O5 | Budget-forced preview, sparse large labels, immutable source arrays | Estimated array payload within budget; exact source arrays and surviving label IDs; all channels/masks aligned; impossible budget rejects. |
| O6 | Saved/unsaved runs, nucleus/cell and raw/envelope choices | Download CSV retains IDs/values/QC; JSON agrees with selected mask/config and includes hashes; changing selection invalidates cached features. |
| O7 | Existing default tests, lint, notebook check, mypy and headless Web startup | Affected/full regression and lint pass; compare mypy to baseline without suppressing diagnostics. Renderer/model tests remain distinct evidence. |

Deterministic fixtures require each listed case, not biological N or confidence
intervals. Real model equivalence, runtime/memory performance on large volumes,
annotated biological accuracy, assay validity and inferential coverage remain
NOT ASSESSED unless new qualified evidence is supplied. The existing <5% surface
accuracy specification remains unchanged and its known failure is not resolved
by these repairs. Tier D review: LIMITED INDEPENDENCE.
