# Devlog 2026-04-29 — Side Connector Facet Connectivity Audit

**Date:** 2026-04-29
**Topic:** Read-only audit of side-connector FEM element to adjacent shell-facet connectivity
**AIs involved:** Codex
**Codex Audit Confidence:** high for source-level findings, medium-high for modeling rationale
**Literature References:** Monk 2003 §5.5.1; Alves et al. 2022b (paper6) §3.1

## Summary

Audited the new `ThinShellFactory` facet-table path and `MaxwellFactory::connect_side_connectors_with_facets()`. The concept is correct: `LeftCoating`/`RightCoating` FEM elements need a pointer to the adjacent original shell facet so connector recovery can use the shell-side surface geometry/basis context for the perpendicular/binormal H component that is not represented as a standard thin-shell mapped trace. The current implementation is not yet safe.

## Key Findings

- Compile blocker: `ThinShellFactory::create()` calls `create_side_connectors(..., tSideConnectorFacetTable)` and then `append_move(mSideConnectorFacetMap, tSideConnectorFacetTable)`, but the declaration/definition of `create_side_connectors()` do not accept that final table argument, and `mSideConnectorFacetMap` is a `Map<id_t,id_t>`, not a `Cell<id_t>`.
- The actual connector-element-to-facet mapping is already populated directly in `create_side_elements()` via `mSideConnectorFacetMap[tElement->id()] = aFacetIDs(k)`, so the local `tSideConnectorFacetTable` path appears stale/redundant.
- `create_facet_table()` can silently leave a segment mapped to facet id `0` if no matching facet is found. That should be an explicit error before the map is consumed.
- `create_facet_table()` loops through `tA->facet(f)` and does not verify that the matched facet belongs to the passed `aFacets` set. This can select the wrong adjacent facet on non-manifold or overlapping sideset topology.
- `MaxwellFactory::connect_side_connectors_with_facets()` first sets facets directly on every rank using `mThinShellSideConnectorFacetMap[elem_id]`. In parallel, non-root ranks do not own the root-populated map; `operator[]` inserts a default id if missing. This can assign facet 0 or fail before the later root-mediated distribution runs.
- The second half of `connect_side_connectors_with_facets()` is the right parallel pattern: collect connector element IDs, map them to facet IDs on rank 0, distribute facet IDs back, and set `Element::facet()` locally. The first direct loop should be removed or restricted to serial execution with checked map access.

## Changes Made / Proposed

- No source-code changes made.
- Proposed fixes: remove the stale `tSideConnectorFacetTable` handoff; add checked facet-table generation; and make `connect_side_connectors_with_facets()` use only the serial branch or the existing distributed branch depending on `comm_size()`.

## Open Questions

- Should connector elements be linked to the original geometry-only shell facet only, or to a layer/block-specific facet view if duplicated shell interfaces return? Current connector preconditions disable ghost facets, so original facet mapping is currently consistent.

## Files Updated

- `devlog/dl20260429_side_connector_facet_connectivity_audit.md`
- `devlog/README.md`
