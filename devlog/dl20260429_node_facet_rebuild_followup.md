# Devlog 2026-04-29 — Node-Facet Rebuild Follow-Up

**Date:** 2026-04-29
**Topic:** Read-only audit of `ThinShellFactory::update_node_facet_tables()` and its interaction with `create_facet_table()`
**AIs involved:** Codex
**Codex Audit Confidence:** high for code-path tracing, medium-high for BELFEM-convention interpretation
**Literature References:** N/A

## Summary

Audited the new shell-local node→facet adjacency rebuild added to `ThinShellFactory`. The new helper materially improves `create_facet_table()` for the current connector path: the call-site placement is acceptable, the intended original-node lookup is now direct, and no other pre-finalize consumer on this path appears to depend on full mesh adjacency.

Two issues remain inside the helper:

- duplicate flag-2 reset is wrong on line 2934, so duplicate facet containers are not allocated robustly if a duplicate enters with stale flag state;
- the index backup/restore tail is dead or unfinished, because no node indices are modified and the restore loop ignores the saved values.

## Key Findings

- `update_node_facet_tables()` is inserted after `update_node_edge_tables()` and before `create_side_connectors()`, which is acceptable for the current connector flow because the only later `Node::facet()` reader in this path is `create_facet_table()`. See `src/mesh/cl_ThinShellFactory.cpp:172-174`, `302-315`, `3040-3044`.
- The placement is not a general “full adjacency is current” point when ghost facets exist, because `create_ghost_facets()` runs later. That does not affect the current connector path because `create_side_connectors()` rejects ghost facets. See `src/mesh/cl_ThinShellFactory.cpp:249-276`, `1982-1985`.
- Duplicate bookkeeping is not airtight: inside the duplicate-reset loop the code unflags `tNode`, not `tOrg->duplicate(d)`, so stale duplicate flag 2 can suppress `allocate_facet_container()`. See `src/mesh/cl_ThinShellFactory.cpp:2931-2935`, `2963-2967`.
- The shell-only adjacency has narrow blast radius before mesh finalization. `MaxwellFactory::create_thinshells()` does not appear to read `Node::facet()` after `tFactory.create()` and before `mMesh->finalize()`. Full node→facet adjacency is rebuilt during finalize. See `src/fem/maxwell/cl_MaxwellFactory.cpp:863-903`, `src/mesh/cl_Mesh.cpp:824`, `src/mesh/cl_Mesh_ConnectivityCalculator.cpp:271-303`.
- Giving duplicates the same shell-facet adjacency is reasonable for this helper’s local shell-only view, but it is not a global BELFEM invariant. Full-mesh adjacency can legitimately differ between originals and duplicates after relink/update. See `src/fem/maxwell/doc/thin_shell_facet_orientation.md:195-199`.
- The index backup `tIndices` is unused in practice: no node indices are modified, and the tail loop writes `0..n-1` instead of restoring the saved values. See `src/mesh/cl_ThinShellFactory.cpp:2921-2928`, `2990`, `529-533`.
- The `sideset_id()` filter in `create_facet_table()` is now redundant by construction but still acceptable as a cheap defensive guard. See `src/mesh/cl_ThinShellFactory.cpp:217-224`, `2972-2984`, `3029-3045`.

## Changes Made / Proposed

- No source-code changes made.
- Audit recorded in `todo/ai_exchange.md`.

## Open Questions

- Should `update_node_facet_tables()` be kept as a shell-local helper only, or should it be renamed/commented to make that restricted scope explicit?
- If ghost facets are ever re-enabled together with side connectors, should the rebuild move later or remain intentionally shell-only?

## Files Updated

- `todo/ai_exchange.md`
- `devlog/dl20260429_node_facet_rebuild_followup.md`
