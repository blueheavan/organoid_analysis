# Round 3 report — spacing defect, record re-issue, volume audit, gate semantics

2026-09-13. Three commits, all local (push blocked in this environment).

| Item | Status after this round |
|---|---|
| Per-axis voxel spacing | Defect fixed, demonstrated, 17 regression tests |
| Canonical analytical record | Re-executed at the fixed source; bit-identical results |
| SG-1 surface area | FAIL (unrestricted, worst 25.81%) — unchanged |
| SG-2 volume | FAIL / NOT QUALIFIED — audited, no algorithm change |
| SG-3 segmentation | FAIL / INSUFFICIENT EVIDENCE — protocol addendum written |
| SG-4 viability | NOT ASSESSED — unchanged |
| SG-5 study design | INSUFFICIENT EVIDENCE — unchanged |
| SG-6 acquisition calibration | FAIL / INSUFFICIENT EVIDENCE — study proposed |
| Overall gate | **NOT PASSED (0/6)** |

## 1. `be223ba` — per-axis voxel spacing

`ome_spacing()` returned `None` when any single axis was missing, discarding
the Y/X values the file did state. Because the manifest-conflict check required
complete metadata, it was then skipped entirely: a manifest value silently
overrode a *stated* file axis.

Demonstrated against the previous commit, not merely asserted: a file stating
Z = 2.0 µm with Y/X omitted, loaded with a manifest declaring Z = 1.0 µm.
Pre-fix, the load succeeded and used 1.0 µm — a 2× single-axis error in every
physical measurement, no warning. Post-fix it is refused, naming the axis and
both values.

The fix separates two rules the old shape conflated: **completion** of axes a
file omits, and **comparison** on axes it states. Neither is skipped because
the other axes are incomplete. Mixed provenance is preserved in
`spacing_source_by_axis`. Measurement behaviour is unchanged wherever metadata
was already complete. 17 tests cover both sources complete and agreeing, both
conflicting, partial completed from manifest, partial conflicting on a stated
axis, neither source resolving an axis, invalid stated values, multi-source
merging, and the same rule on companion channels and annotation files.

Recording provenance is not verification — a manifest value carries the
authority of whoever typed it. That is SG-6.

## 2. `10fb983` — record re-issue

The record was stale relative to the source. It was re-executed at the
committed source and is canonical again.

**This is a re-execution, not new independent confirmation.** All 1128 rows are
bit-identical to the superseded record across area, volume, sphericity and
`rho_in`. It restores the binding between record and source under test;
confidence in the estimators is unchanged. The superseded directories are
retained unmodified.

Both gate items are reported exactly as re-derived: SG-1 worst 25.81 %, median
signed −1.24 %; SG-2 per stratum 2.57 % / 2.74 % / 18.26 %.

## 3. `3e975a2` — volume audit, gate semantics, study proposals

Analysis of existing evidence. No algorithm, criterion, domain constant,
phantom, acceptance rule or gate behaviour changed; no new measurement.

**Volume (SG-2) stays FAIL / NOT QUALIFIED.** The audit reports the prospective
evidence that does exist rather than understating it: the 2026-09-12 freeze
record predeclared `|volume err| < 1%` and the confirmation set, executed once,
met it on 96 in-scope smooth cases (worst 0.76 %).

The finding that matters is why that does not qualify volume: **the
machine-checkable part of the domain does not bound volume.** Inside
`rho_in >= 10` and anisotropy <= 4 — everything software can verify from a mask —
creased confirmation cases reach 3.20 %. The 0.76 % holds only because the
*smoothness* clause excludes them, and no mask-computable crease indicator
exists. An object can satisfy every checkable condition, fire no flag, and
carry three times the criterion. On the frozen grid the error is set by
resolution and its sign by sub-voxel grid phase (one shape: −18.3 % at one
offset, +7.3 % at another), so it is a per-object quantization error, not a
correctable bias.

`rho_in >= 8` passes on the frozen grid (n = 67, max 0.53 %) but was found by
scanning data already seen; it is recorded as a post-hoc hypothesis, labelled
exploratory, and is **not** declared a guarantee. Per the above, the threshold
is not the binding problem anyway.

**Gate semantics: analysis and recommendation only, decision record unsigned.**
The gate keeps its current unrestricted behaviour until the owner records a
choice. The ambiguity — not either number — is the defect: "FAIL" (unrestricted)
and "qualified for its declared domain" (predeclared confirmation) are both
true of different propositions. Four options are quantified; the recommendation
is to split each quantity into a domain-restricted qualification item and an
unrestricted characterization item that reports a number rather than PASS/FAIL.
Adopting it **does not turn the gate green**: the volume qualification item
would be FAIL / INSUFFICIENT EVIDENCE for a more accurate reason than it is red
now, and SG-3 to SG-6 are unaffected.

**Study proposals assert nothing.** The SG-3 addendum adds per-axis spacing
provenance as a covariate and requires the domain flag to be reported
stratified rather than used as a silent exclusion. The SG-6 proposal is new:
traceable lateral and axial standards imaged in the specimen's mounting medium,
strata per instrument and objective, imaging session as the independent unit,
scale deviations propagated into area and volume. A 3 % *lateral* scale error is
a +6.09 % volume error and a +3 % to +6.09 % area error depending on object
orientation (+4.02 % for a sphere); a 3 % error in all three axes is +9.27 %
volume and +6.09 % area. Either way it lands on every object, is invisible to
every check in this repository, and exceeds the qualified estimator's whole
in-domain budget. (The figures in the first issue of this report, "~3 % area and
~9 % volume" for a lateral error, were wrong — see
`docs/evidence/2026-09-13-volume-audit/ERRATA.md`, E-5.)
Calibration is recommended before the segmentation study, since a scale error
and a mask bias both appear as a systematic volume discrepancy.

## 4. Checks

`ruff` clean. Type-check ratchet unchanged (110 signatures, 0 new, 0 fixed) —
typing debt was treated as non-blocking and was not lowered. Notebook
validation not applicable. 124 record, gate and domain tests pass. Gate output
unchanged.

Pre-existing failures, unchanged and environmental to this sandbox: several
segmentation and export-contract tests fail on a permission error against
`~/.cellpose`. They may or may not reproduce on your runner.

## 5. Open decisions for the project owner

1. `docs/GATE_SEMANTICS_ANALYSIS.md` §6 — which gate-semantics option, signed
   and dated. No implementation may precede it.
2. Whether the volume domain is the surface domain, given that the scope clause
   carries volume as well as area.
3. Whether analytical volume qualification is resourced now or alongside SG-3.
4. Whether the 1 % volume criterion is the biological requirement or inherited
   — loosening is permissible only *before* a confirmation set is measured,
   with the justification recorded.
