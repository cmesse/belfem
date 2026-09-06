# Devlog 2026-08-04 — Side Connectors in the .bfm: Gap Analysis and Closure

**Date:** 2026-08-04
**Topic:** Three-AI gap analysis of side-connector persistence in the `.bfm` mesh cache,
Christian's per-shell connector record, its audit (two save-path P0s), and the closure of
B1–B6; only the stale-cache detection (B7) stays open, deprioritized.
**AIs involved:** Claude (primary), Codex + Grok (independent analysis and audits, every
round)
**Claude Confidence:** high on everything file:line-cited and 3/3-confirmed; the B5 slot
rule encodes a factory convention in the loader (cross-referenced comments both ends)
**Literature References:** none new

## Summary

The `.bfm` reload path never rebuilds side connectors (`create_thinshells` early-returns
when shells came from the file), so whatever the file does not carry is absent for the
whole run. A three-way gap analysis (full agreement, no finding contradicted) showed the
raw entities mostly survive — HEX8TB elements/blocks incl. domain types, recovery facets
with master/slave/orientation, facet ids (so the id+1 invariant holds), hanging outer
entities with sources and weights — while the shell-side associations and the connector
width did not. Christian implemented the per-shell record; the audit caught two P0s in
the save path before it ever ran; B3/B5/B6 followed the same day.

## The schema (thinshells group)

`blockids` renamed `layers` (deliberate pre-1.0 break, no shim — an old file fails loudly
at thin-shell load; remedy = delete the cached `.bfm`, documented in
`src/mesh/doc/bfm_file_format.md` §6). New optional per-shell datasets, written only when
a shell has connectors and read under a `dataset_exists` guard: `coatings` (connector
block ids), `seams` (recovery sideset ids, index-aligned with coatings), `widths`
(connector block thickness — the only record of it, since block thickness is not part of
the group data and `ThinShell::set_thicknesses` reaches layer blocks only). Sparse rows
work because both `hdf5::Dataset` constructors value-initialize to `{len=0, p=null}`.

## Key results

- **Audit of the record (Codex + Grok + Claude, 3/3):** two P0s in the save loop — the
  freshly allocated output buffer was used as the container index
  (`side_connector_blocks()( tIDs[k] )` → `( k )`), and `side_connector_sidesets()` was
  empty because the factory never registered it (the fresh-path half of the gap). Both
  compile-silent; both fixed. The loader was verified correct from the start.
- **Factory registration:** `create_side_connectors` now pushes block and sideset
  together, one per curve, with the index-alignment contract commented (the bfm pairs
  them by position).
- **B3 — recovery facet edges after reload:** re-derived, no new dataset. The facet's two
  dof edges are exactly the wall's inner-face edge slots, and the slave face index is
  persisted: `Element_HEX8TB::get_edges_of_facet` gained the lateral faces 0 (slots 0,2)
  and 2 (slots 1,3), bottom edge first, and the coating loop rebuilds each facet's QUAD4
  edge container from `slave()->get_edges_of_facet( index_on_slave() )` — pointer-exact
  against the fresh path for both connector signs and both orientation branches.
- **B5 — twin-sheet ambiguity, resolved without persistence:** the wall mixes edge sheets
  by a deterministic SLOT rule (bottom slots 0/1 take the material-interface duplicate
  sheet when it exists, top slots 2/3 always the primary — the tA–tD selection in the
  factory). `reconstruct_edge_connectivity` now applies the same rule: HEX8TB slots 2/3
  resolve via a new find-first lookup (first = primary, because primaries are appended to
  the mesh before duplicates and the stable sort keeps container order), everything else
  keeps the accepted last-wins. No-op when no duplicates exist.
- **B6 — Christian's design call:** rather than teaching the pre-enrichment topology
  snapshot about connectors, the factory refuses to cache enriched meshes:
  `BELFEM_ERROR( ! mUseEnrichment, ... )` before `mMesh->save()`. Cut enrichment is
  currently dead — the implemented bubble functions are not the right space for it.
  The residual type-map leak of connector entities on reload is benign:
  `Left/RightCoating` fall through `default: pass` in `Topology::select_blocks`, and
  no consumer queries those buckets (verified).
- **B7 stays open** (deprioritized): the checksum hashes only the original gmsh geometry,
  so connector-setting changes cannot invalidate a stale cache; wants a processing-options
  stamp and a `/meta/format_version` reader gate.

## Changes made

- `src/mesh/cl_Mesh_BfmFile.cpp` — per-shell connector record (save + load), the two P0
  fixes, facet edge rebuild in the coating loop, redundant layer `set_thickness` commented
  out with a pointer to the surviving write in `ThinShell::set_thicknesses`.
- `src/fem/kernel/cl_ThinShellFactory.cpp` — sideset registration next to the block push.
- `src/mesh/cl_Element_HEX8TB.hpp` — lateral faces 0/2 in `get_edges_of_facet`.
- `src/mesh/cl_ProtoMesh.cpp` — find-first lookup + HEX8TB slot rule in
  `reconstruct_edge_connectivity`.
- `src/fem/maxwell/cl_MaxwellFactory.cpp` — enriched-mesh save guard.
- `src/mesh/cl_Block.hpp` — thickness reverted to a scalar (the 3-slot array from the
  previous session was never needed: the layer thickness routes through the recovery
  facet's master).
- `src/mesh/doc/bfm_file_format.md` — `layers` rename, the three new datasets with their
  alignment contract, and the versioning note.
- `todo/hex8tb_phase2_fem_wiring.md` — B-tracker updated live (B1–B6 ticked, B7 open).
- Thread: `tmp/ai_exchange/sideconnector_bfm_gaps.md` (gap analysis, three audit rounds,
  fix record; ephemeral).

## Open items

B7 (stale-cache detection). Validation for the whole set: a save→reload cycle with the
WIP factory stop temporarily gated — that stop is why compile-silent save-path defects
could not have surfaced in any run yet. Then the r′ wall-term theory session (R4/R5).
