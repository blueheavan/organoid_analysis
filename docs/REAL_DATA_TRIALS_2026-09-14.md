# Public 3D microscopy trial — 2026-09-14

## Purpose and status

This is an end-to-end input and execution-capability trial, not segmentation
validation. It supplies no reference labels, independent annotations,
calibration study, biological control, or held-out confirmation set. Therefore
it cannot change SG-3, SG-4, SG-5, SG-6, or any scientific qualification item
to `PASS`.

## Dataset

- Source: [Zenodo record 20624491](https://zenodo.org/records/20624491), DOI
  `10.5281/zenodo.20624491`.
- Description supplied by the depositor: 4-channel 2-photon fluorescence
  kidney-organoid OME-TIFF, 16-bit CZYX; DAPI is channel C1 and Phalloidin is
  a cell-boundary channel.
- The OME series consists of four companion files. A single companion file is
  not a valid standalone input; all four must be placed in the same directory.
- Downloaded C1 companion SHA-256:
  `8f7183dd99f5ea9741727b6c33d0e9bea3cfc562a1e564e6db363b5224e31ea2`.

## Read-path result

With all companion files present, the repository reader accepted the series as
`CZYX = 4 × 267 × 512 × 512`, `uint16`, finite, and reported physical spacing
`Z/Y/X = 0.498127/0.576023/0.576023 µm`. The inferred Z:XY ratio is 0.865.
This verifies format/axis/metadata handling only.

The first one-file attempt failed safely because `tifffile` could not resolve
the other OME companion files. That failure was caused by incomplete download,
not image quality. The reader now tolerates an OME `AnnotationRef` that does
not affect pixel or spacing semantics; this is covered by a regression test.

## Fixed inference attempts

All attempts used C1, the metadata-derived anisotropy 0.865, native pixels,
the shipped `cpdino-vitb` model, and unmodified default segmentation
thresholds. They were engineering capability attempts only; no parameter was
tuned against the image.

| Attempt | ROI ZYX | Result |
|---|---:|---|
| 1 | 64 × 512 × 512 | cancelled after approximately seven minutes in DINO attention forward pass |
| 2 | 32 × 256 × 256 | cancelled after approximately four minutes in DINO forward pass |
| 3 | 16 × 128 × 128 | cancelled after approximately three minutes in DINO attention forward pass |

The runtime detected `cpu`. A subsequent direct Torch probe emitted the known
duplicate `libomp.dylib` initialization error before it could report MPS status.
The trial therefore establishes neither a successful Cellpose segmentation nor
a scientific failure of the image. It establishes an environment-performance
blocker for this host.

## Conditions before a real-data segmentation claim

1. Run the same locked environment with a verified MPS or CUDA accelerator;
   resolve the duplicate OpenMP runtime rather than setting the unsafe
   `KMP_DUPLICATE_LIB_OK` override.
2. Keep all files of a multi-file OME series together and retain original
   physical-spacing metadata.
3. Select a nuclei channel with through-depth signal and no widespread
   saturation; for cell measurements use a registered membrane/cytoplasm
   channel acquired on the same grid.
4. Record image-level QC: dynamic range/saturation, depth attenuation, object
   truncation, fragmentation, and spacing provenance.
5. For SG-3, use a predeclared independent 3D annotation study. This public
   image has no supplied instance-mask reference and cannot be repurposed as
   ground truth.

## Other public candidate considered

The [Cellos repository](https://github.com/TheJacksonLaboratory/Cellos)
documents a publicly available confocal organoid example with EGFP, mCherry,
and brightfield channels. Its download is approximately 11 GB and requires the
source XML metadata to reconstruct the well/channel/Z layout. It was not
downloaded for this host trial because the present CPU-only Cellpose runtime
cannot complete even the much smaller fixed ROI attempts above. It remains a
potential future engineering/interoperability dataset, not an SG-3 reference
standard unless independent annotations and the complete acquisition protocol
are separately qualified.
