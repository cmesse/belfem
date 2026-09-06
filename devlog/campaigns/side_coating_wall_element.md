# Campaign: Side-Connector Wall Element (HEX8TB)

**Status seed:** 2026-08-05 from dl20260728/0729/0731/0803/0804 — every claim
`[seeded — confirm]` until Christian's correction pass.
**Branch:** `sideconnectors` · **Plan:** `todo/hex8tb_phase2_fem_wiring.md`

## Current accepted design `[seeded — confirm]`

Degenerated HEX8TB wall element bridging stacked thin-shell layers at the tape sides:
4 longitudinal edge dofs, +curve oriented; exact-cuboid metric (width and layer thickness
are exact block data, node positions partly a drawing proxy); j_t ≡ 0, h_t is a stream
function. Hidden sideset/curve; recovery facet is READ-ONLY (no weak-form terms) with
**facet id = wall id + 1**, master = layer-block element — thickness and ownership
contracts route through the master. Ownership: connector blocks selected per shell in
`Kernel::partition_mesh`. Doc: `src/fem/maxwell/doc/side_coating_wall_element.md`.
(The 2026-06-05 removal verdict applied to the OLD HEX8TS wrap — free binormal-H at the
fold; the revival is the collapsed-wall element, not that geometry.)

## Last passing reproducer `[seeded — confirm]`

None end-to-end yet — construction runs to the deliberate WIP factory stop. Executable
evidence to date: `.bfm` persistence round-trip of the per-shell record (dl20260804);
battery `make check-fast` covers the QUAD4TS/PENTA6TS/HEX8TS conventions the wall couples
to; `Hex8TbUnitCirculation` was DISABLED as the first executable evidence of the R4 gap and
is **re-enabled since 2026-08-10** — the interpolation-factory case landed (DR-47 closed), so
this test is now the wall element's only automated check. Its first green run doubles as the
`EXPECTED` sign-off.

## Open P0/P1

- Christian's numeric smoke run of the construction (DR-18): MeshChecker volumes both
  signs, `compute_orientation`, facet resolution.
- R5 plumbing checklist P1–P9 (`todo/hex8tb_phase2_fem_wiring.md`): the
  `h_side_connector` kernel scaffolding landed 2026-08-05 (Christian) and survived the
  frozen jury round (`review_h_side_connector`; fixes applied same day) — what remains
  is wiring: connector facet linking (id+1 invariant), group selection + FieldList,
  link_to_group dispatch, material assignment, input keys, WIP-stop removal, Newton
  part (DR-19).
- ~~Lagrange interpolation-factory case for HEX8TB (R4 / DR-47) — also blocks the
  disabled `Hex8TbUnitCirculation` test.~~ **Both done** (case landed 2026-08-08, test
  re-enabled 2026-08-10).
- B7 stale-cache detection for connector settings (DR-21); save→reload validation with
  the WIP stop gated (DR-22).

## Superseded approaches `[seeded — confirm]`

HEX8TS side-wrap (removed 2026-06-05, unphysical free binormal-H); hand-rolled facet
node order in the .bfm (loader re-derives from master); Lagrange-multiplier coupling
never considered (static condensation house rule).

## Dated entries

dl20260728_side_connector_revival · dl20260729_side_connector_fixes_hex8tb ·
dl20260731_hex8tb_edge_function · dl20260803_side_connector_recovery_facets ·
dl20260804_sideconnector_bfm_persistence
