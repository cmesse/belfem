# Devlog 2026-04-29 — Connector Facet Trace Follow-Up

**Date:** 2026-04-29
**Topic:** Read-only trace of `LeftCoating` / `RightCoating` FEM-element connectivity to adjacent shell facets
**AIs involved:** Codex
**Codex Audit Confidence:** high for source-level findings, medium-high for modeling rationale
**Literature References:** Alves et al. 2022b (paper6) §3.1; Monk 2003 §5.5.1

## Summary

Confirmed the user's physical rationale: the shell/connector fold already shares the longitudinal edge DOFs, but the binormal / interface-normal field component is not enforced by that shared-edge topology. The connector therefore needs a valid pointer to the adjacent shell facet so BELFEM can recover or weakly couple that missing component through the facet/sideset normal-calculator infrastructure.

The current code only partially achieves that. `ThinShellFactory::create_side_elements()` populates the intended connector-element -> facet-id map, but the upstream handoff is stale and the downstream Maxwell hookup is unsafe in parallel.

## Key Findings

- `ThinShellFactory::create()` still calls `create_side_connectors(..., tSideConnectorFacetTable)` and then `append_move( mSideConnectorFacetMap, tSideConnectorFacetTable )`, but the declaration/definition of `create_side_connectors()` have no such output argument. This is a source-level mismatch in the current tree. See `src/mesh/cl_ThinShellFactory.cpp:303-318`, `src/mesh/cl_ThinShellFactory.hpp:303-315`, `src/mesh/cl_ThinShellFactory.cpp:1958-1969`.
- The actual per-element mapping is populated directly in `ThinShellFactory::create_side_elements()` via `mSideConnectorFacetMap[ tElement->id() ] = aFacetIDs( k )` in both orientation branches. That is the correct place to build the mapping. See `src/mesh/cl_ThinShellFactory.cpp:2670`, `2709`.
- Reusing one shell facet id for all connector rows at the same curve segment is consistent with existing thin-shell facet linking, where each through-thickness copy of a shell element reuses the same in-plane surface facet. See `src/fem/kernel/cl_FEM_DofMgr_BlockData.cpp:408-425`.
- `ThinShellFactory::create_facet_table()` is plausible but fragile: it can silently leave facet id `0`, it relies on `Node::facet()` adjacency, and it implicitly assumes the side-curve nodes being flagged are originals because it flags `tA`/`tB` directly but checks `tFacet->node(...)->original()->is_flagged(1)`. See `src/mesh/cl_ThinShellFactory.cpp:2940-2966`; stale-adjacency warning in `src/fem/maxwell/doc/thin_shell_facet_orientation.md:193-201`.
- `MaxwellFactory::connect_side_connectors_with_facets()` has a real MPI bug. It first dereferences `mThinShellSideConnectorFacetMap[ elem_id ]` on every rank before the root rank converts and redistributes facet ids. Since `Map::operator[]` default-inserts missing keys, non-root ranks can fabricate facet id `0` and then hit `Mesh::facet(0)`. See `src/fem/maxwell/cl_MaxwellFactory.cpp:2679-2687`, `2720-2758`, `src/containers/cl_Map.hpp:146-150`, `src/mesh/cl_Mesh.hpp:1320-1327`.

## Changes Made / Proposed

- No source-code changes made.
- Proposed direction:
  - Remove or reconcile the stale `tSideConnectorFacetTable` handoff in `ThinShellFactory::create()`.
  - Make `create_facet_table()` validate "found exactly one adjacent shell facet" rather than silently returning `0`.
  - In `connect_side_connectors_with_facets()`, use either the checked serial path or the existing root-distribution path, but do not dereference `mThinShellSideConnectorFacetMap[ ... ]` on non-root ranks before redistribution.

## Open Questions

- Are `aProtoShell->side_curves()` guaranteed to remain original-node curves under all CutFactory paths that will coexist with side connectors? `create_facet_table()` currently depends on that invariant implicitly.
- Will connector recovery use only the facet pointer for postprocessing / `compute_bn`-like reconstruction, or is the longer-term target the `h_penalty()` fold coupling path as well?

## Files Updated

- `todo/ai_exchange.md`
- `devlog/dl20260429_connector_facet_trace_followup.md`
