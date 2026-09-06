# BFM save/load re-audit + plan refresh

**Date:** 2026-06-28
**Purpose:** Re-audit the BFM mesh-file save/load stack against the current tree and refresh
`todo/meshfile_refactor_plan.md`. Audit-only — no source modified.
**Module:** `src/mesh` (`cl_Mesh_BfmFile`, `cl_ProtoMesh`), `src/io` (`cl_HDF5_Dataset`)
**AIs:** Claude (primary), Grok (third voice), Codex (precision + prose). Consensus, no dissent.

## What this was

A standing-back re-audit of `BfmFile`/`ProtoMesh`/`Distributor` plus a structural refresh of the
plan: reconcile R0–R12 / D1–D7 against the code, delete the now-dead legacy bug catalogue
(Appendix A) and its cascade, and correct the "still unimplemented" list.

## Net findings (new defects — all latent; `Mesh::save/load` still errors `"not implemented"`
for `.bfm` at `cl_Mesh.cpp:152-155,352-354`, so `BfmFile` is unwired and the save path is never run)

- **D8 (CRITICAL)** `save_element_data` never advances the vlen row index: `set_size(tCount,n)` at
  `cl_Mesh_BfmFile.cpp:298` with no `++tCount` in `:296-304`. `Dataset::set_size` uses the explicit
  index (`cl_HDF5_Dataset.hpp:86`). ⇒ every element targets row 0 (debug assert on element #2 /
  release clobber+leak). Element topology never round-trips.
- **D9 (CRITICAL)** `save_facet_data` elements/indices loop (`:436-492`, set_size `:457-479`) — same
  missing per-facet `++tCount`; dormant topology branch `:418` too.
- **D10 (CRITICAL)** `save_face_data` — no `++tCount` (`:785-839`) **and** `close()` instead of
  `save()` (`:841-842`) → face elements/indices freed without `H5Dwrite`.
- **D11 (LOW)** `load_face_data` opens `elements`/`indices` without the per-dataset guard the facet
  loader has (`:861-862` vs `:571-572`).
- **D12 (HIGH, data loss)** `save_control_point_data` fills `tTopo` (`:955`) but **never calls
  `tTopo.save()`** before `close_active_group()` (`:966`) → element→control-point incidence never
  written (ids/coords are; points reload without element links). *Grok's catch; I had earlier marked
  this save path "safe" — wrong.* Same class as D10(b).

D8–D10 are the exact row-index / close-vs-save defects recorded **fixed on 2026-06-27**
(Dataset-refactor pass); the current tree does not contain those fixes — apparently lost in a later
file restructure. Correct by contrast (verified): `save_edge_data` (flat index `e`),
`save_node_duplicate_data`, and all five `save_hanging_*` (`tIDs(tCount++)`).

## Status delta recorded in the plan

- Done: D1 (node-dup mechanism retired, replaced by `save/load_node_duplicate_data`), D2
  (`create_sidesets` moved to load), D3 (`reconstruct_edge/face_connectivity`), D4 (single
  `finalize()` on load); R2/R3/R4 (realized as `BfmFile`+`hdf5::Dataset`, not separate writer/reader
  classes). Open: D5–D7, R5–R12.
- Obsolete: **R0** (struck) — legacy `cl_Mesh_HDF5{Writer,Reader}` deleted.
- De-facto: **R1** — an implemented schema exists; §6 is now partly stale.
- Implemented (out of "missing"): node duplicates, hanging entities (= static-condensation/T-matrix
  data, R-row 9), control points (**partial** — D12 breaks incidence save). Still missing:
  vertices, curves+segments, thin shells, periodicity, fields/globals/time, block↔material.
- Element **type** is per-**block** (`blocks/types`), not per-element (`ElementData::mTypes`) — still
  resolves the row-3 ambiguity since each block is single-type (plan corrected).

## Plan-doc edits (`todo/meshfile_refactor_plan.md`)

Added D8–D12 to §4.0 (dated, severity, file:line) + a master-proc-guard note (`BfmFile::save/load`
have no `comm_rank()==master` guard — later pass). Reconciled R0–R5 status. **Deleted Appendix A**
(tombstone left; the one shared item, W6 `get_array_size`, flagged for re-check only if still used).
Cascade: §1 and §3 got "HISTORICAL — dead citations" banners (§3's a/b/c classification kept; only
the "Saved today?" column + HDF5Writer/Reader cites are dead); O5 updated (control-point path now
exists, with D12 + a 2D-`tCoords(2,p)` OOB caveat in `create_control_points`); O7 already
RESOLVED→retire; R10 renamed `ProtoMeshHDF5Writer`→`BfmFile`.

## AI cross-check

Claude found D8–D11 + the status delta; Grok independently confirmed D8–D10 and **found D12**; Codex
confirmed D8–D12, the unimplemented list, and the `Mesh::save/load`-unwired fact, and supplied the
prose/accuracy fixes applied above (element-type-per-block; "have been retired"; R10 naming; O5).
No dissent. One trivial unreconciled citation nit: Codex cites `set_size` at `cl_HDF5_Dataset.hpp:121`,
Claude's direct read puts it at `:86` — same function, behavior identical.

## Next (fix pass, not this session)
D8/D9/D10/D12 are the blockers for any round-trip — fix the row-index increments + the missing/
wrong `save()`, then R5 (base-mesh save→load→save) which would have caught all of them.
