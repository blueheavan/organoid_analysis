# Scientific gate reassessment — 2026-09-26

## Decision from the current evidence

The request to change algorithms and pass the scientific gate authorizes a
method review; it does not supply the independent reference observations that
the existing nine-item gate requires. The gate remains `NOT PASSED`. Its
analytical record is currently rejected because the source fixes from
2026-09-25 differ from that record's sealed hashes. This is an evidence-integrity
result, not a new analytical-accuracy measurement. No gate status, threshold,
Crofton weight, or statistical decision rule was changed in this reassessment.

The mathematical counterexample below also rules out a universal per-object
`<1%` continuous-volume guarantee from a binary mask over the current
`rho_in >= 10` smooth-sphere domain. The production voxel count remains a
well-defined **discrete mask volume**. Calling it an exact continuous-object
volume would exceed what the input identifies.

## Binary-volume non-identifiability, independently checked

At isotropic 1 µm spacing, sample a centred sphere by testing whether each
integer voxel centre is inside it. Let

- `r_small = sqrt(110) + 0.00001` µm;
- `r_large = sqrt(113) - 0.00001` µm.

No 3D integer-grid point has squared distance 111 or 112 from the origin, so
the two spheres have the **same 25×25×25 binary mask** (4,945 foreground
voxels). Both are smooth and closed; both have true radius above 10 µm. The
current `geometry()` also reports `surface_rho_in = 10.6301458` and its
machine-checkable surface-domain flag as true for their shared mask. Their exact
sphere volumes are approximately 4,832.58 and 5,031.59 µm³, a 4.12% ratio.
The intervals within 1% of these two truths do not overlap:

`1.01 × V_small < 0.99 × V_large`.

Any deterministic estimator given only this mask and spacing must return the
same number for both spheres, so it cannot satisfy the `<1%` per-case criterion
for both. The current voxel-count result is 4,945 µm³, with errors above 1%
against each truth. The executable controlled-truth check is
`tests/quantification/test_volume_identifiability.py`. This is an adversarial
**method-limit test**, not a new representative biological validation set or a
post-hoc replacement for the frozen SG-2 confirmation protocol.

The sampling ambiguity agrees with the general non-invertibility of discrete
shape observations discussed in [Fatemi et al., *Shapes From Pixels*](https://arxiv.org/abs/1508.05789).
That source supports the principle; the numerical 3D example and its 1% result
above are derived directly here, not attributed to the paper.

## Algorithm decision

| Option | Consequence | Current decision |
|---|---|---|
| Change the binary-mask volume formula or fit a correction to known phantoms | One mask still represents both spheres; a fitted correction cannot guarantee the existing continuous-volume claim. | Reject as a route to gate PASS. |
| Report voxel count × voxel volume as discrete mask volume | Exact for the stated mask estimand; does not resolve partial-volume or segmentation error. | Retain current algorithm and explicit interpretation. |
| Estimate subvoxel occupancy from raw intensities with an acquisition/PSF model | Potentially adds information absent from a binary mask, but changes the estimand and requires instrument calibration, independent truth, frozen method selection and held-out confirmation. | Design candidate only; no unsupported production switch. |
| Restrict intended use to a better-resolved domain | Could support a bounded claim after a new development study, prospectively frozen domain and independent confirmation; it cannot retroactively convert the existing confirmation set into PASS. | Requires new evidence and a versioned scope decision. |

No published algorithm can substitute for the project-specific validation of
that algorithm. The [3DCellScope study](https://www.nature.com/articles/s41592-025-02685-4.pdf)
illustrates 3D organoid and nuclear segmentation with manually annotated nuclei,
but its complete raw/segmented data are available on request rather than in the
repository, and its reference/task do not establish this project's SG-3A/3B,
SG-4B, or SG-6 calibration outcomes.

## Gate-by-gate evidence still required

| Item | What would be needed for an attributable PASS under the existing plan |
|---|---|
| SG-1a | Source commit and new canonical, domain-restricted analytical surface record on a prospectively frozen confirmation set; exact weight-vector and method provenance. |
| SG-2a | Resolve the continuous-versus-discrete volume claim and the counterexample above; then a prospectively frozen volume-specific qualification record meeting the existing criterion for its stated domain. |
| SG-3A / SG-3B | Separate modality-matched, blinded independent 3D instance annotations, qualified annotator agreement and predeclared sample-level precision; the repository has only synthetic truth masks. |
| SG-4A | A canonical analytical state/refusal record with known-rule cases and denominator checks. Existing unit tests support only `PARTIAL`. |
| SG-4B | Blinded object-level biological reference, orthogonal/graded-insult controls, held-out batches and the predeclared class-specific uncertainty criteria. |
| SG-5 | Freeze the primary statistical method, verify its small-sample implementation independently, then run the predeclared type-I-error and interval-coverage simulations for the actual replicate design. |
| SG-6A / SG-6B | Traceable per-axis spacing-ratio and absolute calibration measurements and channel-registration standards from the intended instruments, with uncertainty and session strata. |

## Executed checks

- Controlled-truth adversarial test plus affected geometry, spacing and
  classical end-to-end inference selection: **32 passed** with the current
  Pixi Python and `PYTHONPATH=src`.
- `ruff check src tests scripts` and `git diff --check`: **PASS**.
- `scientific_validation_gate.py`: **NOT PASSED (0/9)**. SG-1a/SG-2a report
  `FAIL` because the prior sealed record rejects the 2026-09-25 changed source;
  the remaining evidence statuses are as listed above. No old record was
  rewritten and no gate test was relaxed.

The publication evidence can guide a candidate method or study design. It
cannot fill these absent measurements. Until qualified evidence is available,
the overall verdict remains `NOT READY FOR THE SPECIFIED RESEARCH USE`. Audit
separation for this local review is Tier D, `LIMITED INDEPENDENCE`.
