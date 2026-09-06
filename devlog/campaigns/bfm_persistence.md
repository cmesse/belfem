# Campaign: .bfm Mesh Persistence

**Status seed:** 2026-08-05 from dl20260622 plan + dl20260804 — every claim
`[seeded — confirm]` until Christian's correction pass.
**Branch:** `sideconnectors` · **Plan:** `todo/closed/meshfile_refactor_plan.md` (base refactor DONE)

## Current accepted design `[seeded — confirm]`

Cache-reload of processed meshes: topology type map = PRE-enrichment snapshot (reload
uses `Topology::run_on_enriched_mesh`, never plain `run()`); `shell_NN` sidesets are
empty husks post-TSF; no connectivity persistence by decision (element→edge/neighbor
reconstruction sped up instead — binary search in `find_index_in_unique_cell`);
edge node-pair keys are NOT unique on enriched meshes (twin edges on cut boundary
curves — ties resolve find-first for HEX8TB slots 2/3, last-wins elsewhere).
Per-shell side-connector record since dl20260804: `layers` (renamed from `blockids`,
deliberate pre-1.0 break), optional `coatings`/`seams`/`widths`; caching enriched
meshes is a hard `BELFEM_ERROR` (cut enrichment is dead).

## Last passing reproducer `[seeded — confirm]`

Base refactor: reload verified serial + 2/4/8 procs (2026-07-01, D13–D25 closed).
Connector record: loader verified by review; the save-path P0s were caught pre-run —
**no executable round-trip of the connector record yet** (blocked behind the WIP stop,
DR-22).

## Open P0/P1 `[seeded — confirm]`

- B7 stale-cache detection (DR-21): checksum hashes only geometry; processing-option
  changes (connector flags, thicknesses) silently load stale meshes.
- Save→reload validation cycle with the WIP stop gated (DR-22).

## Superseded approaches `[seeded — confirm]`

Storing connectivity in the file (rejected on size, 2026-07-17 — do not re-propose);
`blockids` dataset name; legacy `HDF5Writer`/`HDF5Reader` (retirement question O7).

## Dated entries

dl20260622_meshfile_refactor_plan · dl20260804_sideconnector_bfm_persistence
