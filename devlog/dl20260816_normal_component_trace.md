# Normal-Component Trace: the Solve Sees h_n, the Plots Don't

**Date:** 2026-08-16
**Purpose:** Record the blind-parallel Codex+Grok trace of the By ≡ 0
observation on thin-shell layer blocks, and the three defects it
registered
**Module:** maxwell, fem/kernel
**Round:** `tmp/ai_exchange/normal_component_trace.md`

## Trigger

Christian, inspecting frame 91 of the tapestack3d quench study: the
field in the ybco layer appeared to collapse to zero at the tape edges;
suspects were the side connectors and the ghost penalty. Claude's
netCDF read of the frame first established the pre-registered facts:
|B| in the layer is textbook critical state (~70 mT at the edges, µT
flux-free core), JJC healthy — but **By and Hy are exactly 0.0 on every
node of every layer block**, and the normal component is physically the
dominant one at a tape edge. His Mac view (Y+ outer sheet at "Bx = 0",
"1.4 T" max) is not in the file (0.089–0.100 T medians on both outer
sheets, global max 0.139 T) — viewer-side candidates recorded, including
Grok's sharp one: 1.4 is this campaign's J/Jc number, so a JJC-colored
view would explain the ×10 "field".

## Verdicts (both voices independently, Claude source-checked the pivots)

- **The solve is intact.** Maxwell residual/tangent and thermal Joule
  both consume b = µ(h_t + h_n) with h_n recovered from the air φ
  (projected: `hn = dot(hn,n)*n`), through `compute_h_ts_edge` /
  `compute_b_ts` / `compute_rho`. The jc(T,|B|,θ) lookup in the solve
  sees the full field. Christian's worst-case scenario — the kernels
  sampling only the ab-plane column — is refuted.
- **Ghost penalty refuted** as a cause: matrix-only, never writes field
  values, and the outermost sheets have no ghost neighbor anyway.
- **Side connectors refuted**: the wall's own residual uses the full
  ht+hb+hn; its postprocessor copies tape values onto wall nodes, never
  the reverse.
- **The gap is the node-field write**: layer H/B ship `E*q` only, and
  PENTA6TS has no through-thickness dof, so the normal row of E is
  identically zero on a flat stack. The recovered h_n — computed and
  used for jc in the same function — is never added to the written H/B.
  Registered as **DR-84** (P2). The fix is nearly free: the projected
  bn exists two lines above the output `combine`.

## The finding that needed the second voice

Codex reported the postproc J/Jc b as "includes the normal component" —
true, but Grok's R2 caught what that hides: the postproc adds the
**unprojected** air average, so the JJC b carries an O(|B_t|)
tangential contamination and the ":791 same-as-solve" comment (added in
0d70a825 with the fold-unification) is false as written. Claude
confirmed against source. Registered as **DR-85** (P2, one-line fix, to
land with DR-69's signed flip so the plot stays the solve's witness).

Third registration: **DR-86** (P3) — Grok's aura-φ lag hypothesis
(phi distributed at step finalize, not per iterate; partition-boundary
shells may assemble h_n from the previous step). Medium confidence,
probe defined.

## DR-69 cross-links (same evening)

The round also confirmed what the git blame had shown earlier: commit
0d70a825 (2026-08-15) switched the postproc angle to the folded
bn_angle, closing the plot-vs-solve *discrepancy* that triggered DR-69
while leaving the ruled signed-dot implementation outstanding — and
removing the output's honest-lookup witness. The DR-69 row now records
this; the eventual implementation flips solve and postprocessor
together, with DR-85's projection fix.

By-catch (not registered): `compute_hn`'s doc comment says "flux
density" but the function returns H; `dof_manager_usage_guide.md`
§h_ts_metal documents a deleted kernel with the unprojected average —
both for the next doc pass.

Status: reviewed, not verified — static trace, no executable gate ran.

## Addendum: DR-85 fixed same evening

Christian's go ("solve DR-85"): the projection line `bn = dot(bn,n)*n`
landed in `compute_superconductor_ts` before `b = mB + bn`, mirroring
`compute_hn`. With it, the ":791 same jc evaluation as the solve"
comment becomes true as written. TU syntax-gated green with the build's
flags.make. Executable gate owed at the next rebuild: one-frame JJC A/B,
differences expected to concentrate at the tape edges. DR-84 (the H/B
node-field write) remains open — that one is Christian's call on whether
it lands before the release.

## Addendum 2: DR-84 fixed same evening

Christian's go: new `compute_conductor_ts` for ThinShellConductor blocks
(bn recovery identical in convention to the superconductor path — air
average, −0.5·µ0, projected onto the shell normal — then `mB += bn`,
`mH += bn/µ`), bound in the dispatch while bulk Conductor keeps the
plain `compute_conductor`. `compute_superconductor_ts` ships the same
addition after the jc lookup (its b was formed before the mutation — no
double count). `copy_seam_fields` propagates the recovered component to
the wall nodes for free. TU gates green in release and debug configs;
Codex verify pass on the diff dispatched. Executable gate for both rows
at the next rebuild: By on layer blocks nonzero, edge-peaked, matching
the air-side By across the interface; JJC one-frame A/B for the DR-85
projection.
