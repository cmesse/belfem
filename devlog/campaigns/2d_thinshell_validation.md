# Campaign: 2-D Thin-Shell Validation

**Status seed:** 2026-08-05 from dl20260710/0804 + todo trackers — every claim
`[seeded — confirm]` until Christian's correction pass.
**Branch:** `sideconnectors` · **Plans:** `todo/2d_thinshell_todo.md`, dl20260710_2d_thinshell_plan_v2

## Current accepted design `[seeded — confirm]`

2-D thin shells as QUAD4TS layer stacks over LINE2 facets, EF_QUAD4TS edge functions
aligned with the validated 3-D convention since `4a42d982`: top-edge dof enters with
**+mS[1]** (shared inter-layer edges keep tangential H continuous), curl parameters are
the true thickness derivatives (−0.5, +0.5), ∇η from actual element geometry (no
axis-aligned assumption), mDetJ = area/4 (quad weight sum 4). Top node-tie uses
master-oriented order (`get_top_nodes` fix, `43474c9f`). Ghost facets only at
material-interface layers; identical-material stacks share edges strongly.

## Last passing reproducer `[seeded — confirm]`

`make check-fast` → `InterfaceOrientation.*` (10 tests): unit circulation, sign sweep,
face activity, Stokes E–C consistency, inter-layer continuity, N=4 Ampère telescope
(the test class that would have caught the alternation bug in minutes), top/bottom
node-tie alignment. Mutation-checked red/green against the historical `-mS[1]` flip
on 2026-08-05. Physics expected values pending Christian sign-off.

## Open P0/P1 `[seeded — confirm]`

- greg2 4-layer smoke rerun (DR-15, Gregory): expect net shell current = I(t),
  same-sign near-uniform layer currents.
- B9: `fix_facet_masters` orientation propagation never runs on all-air 2-D meshes
  (DR-16).
- 2-D edges-on-edges-top path blocked by missing size-1 `to_master_orientation`
  case (DR-14, latent).
- Ghost 12-dof contract test needs a Calculator fixture (DR-46).

## Superseded approaches `[seeded — confirm]`

Pre-4a42d982 EF_QUAD4TS conventions (−mS[1] positional reversal, ±0.5-both curl pars,
axis-aligned ∇η, mSumW=1); facet-order top node tie (pre-43474c9f).

## Dated entries

dl20260710_2d_thinshell_plan_v2 · dl20260710_2d_thinshell_t0_flag_leak ·
dl20260804_2d_thinshell_layer_alternation · (battery: dl20260805_falsification_tooling)
