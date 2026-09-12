# Packaged Crofton weight tables (`crofton_minimax_sym_v3`)

One file per (spacing ratio, stencil radius `m`): `w` is the orbit-resolved
weight vector, `diag` the solver diagnostics (`lp_optimum_t_star`,
`plane_response_max_abs_dev`, `grazing_constant_D_um2`).

**Why these files ship instead of being solved on demand.** The minimax linear
program has a non-unique optimum: an independent re-solve of the same program
returns a different optimal vertex with the same worst-case plane-response
deviation `t*`, and therefore slightly different areas. The qualified method is
this specific set of vectors, not "whatever the solver returns". Re-solving
would let exported areas drift across solver versions with no method-version
bump, so these tables are part of the method identity.

## Provenance

| Tables | Origin |
|---|---|
| 54 files: every (ratio, `m`) exercised by the round-2 development and confirmation sets | Solved during round-2 method development, **frozen before** the confirmation run, and the exact vectors that produced the qualifying evidence in `docs/evidence/2026-09-12-surface-crofton-v3`. |
| 1 file: `1_0.6666666667_0.6666666667_m5.json` (ratio 1 : 2/3 : 2/3, `m` = 5) | Solved at adoption, from the same frozen program and code path. This combination is selected by the repository's own confirmation grid at spacing (1.5, 1.0, 1.0) but was not exercised by the round-2 sets, so no table existed for it. It is **conforming but not evidence-bearing**: no case in the qualifying confirmation set used it. |

Shipping it removes an on-demand solve from the production path; leaving it out
would have meant a non-reproducible fallback for a spacing the project actually
validates against. Cases that use it are covered by the section 9 criterion
through the repository's own V&V grid, not by the round-2 confirmation set.

A measurement reports which tier its weights came from in
`weights_origin` (`packaged` / `cache` / `solved`); anything other than
`packaged` in an export means an on-demand solve, and is not evidence-bearing.
