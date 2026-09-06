# Devlog 2026-06-09 — Periodic Match Nodes Simplification Audit

**Date:** 2026-06-09
**Topic:** Counter-check of simplified periodic `match_nodes()` after thin-cut duplication
**AIs involved:** Codex
**Codex Audit Confidence:** medium-high
**Literature References:** N/A

## Summary

Performed a read-only audit of the proposed simplification for Step 3 of the periodic thin-cut continuity fix: rely on `PeriodicityFactory::match_nodes()` being facet-pair-local instead of making it duplicate-aware.

## Key Findings

- `match_nodes()` is confirmed facet-pair-local: for each source facet it searches only the single paired target facet (`src/mesh/cl_Mesh_PeriodicityFactory.cpp:617-633`).
- Normal original/duplicate twins should not compete inside a single facet because `relink_element()` overwrites element node slots and `Mesh::update_facet_nodes()` refreshes sideset facets from the relinked master element (`src/homology/cl_CutProcessor.cpp:1188-1195`, `src/mesh/cl_Mesh.cpp:692-699`).
- The simplification still depends on an unproven invariant: `map_facets()` pairs facets by centroid geometry only and does not check cut side, hanging status, or original-vs-duplicate role (`src/mesh/cl_Mesh_PeriodicityFactory.cpp:385-465`).
- The new positional alignment in `match_nodes()` is sound when matching is complete: each match assigns equal source/target indices, fills `aTargetNodes(tCount)`, and sorts source nodes by `opNodeIndex` (`src/mesh/cl_Mesh_PeriodicityFactory.cpp:648-658`, `src/mesh/op_Node_Index.hpp:24-28`).
- The completeness check uses `BELFEM_ASSERT`, so it is debug-only and not an always-active release guard (`src/core/assert.hpp:150-166`).

## Changes Made / Proposed

- Updated `todo/ai_exchange.md` with per-claim verdicts C1-C6.
- Proposed adding an acceptance invariant rather than immediately deleting Step 3 duplicate-awareness work: after `Periodicity::update()`, reject any original-to-cut-duplicate periodic node pair and verify paired periodic facets have matching per-node role patterns.

## Open Questions

- Does the translational periodic reproducer pass the role-consistency invariant after the current simplified matcher?
- Should `map_facets()` filter sidesets by periodic domain type or explicitly exclude generated/non-periodic plane sidesets?
- Should `match_nodes()` promote the `tCount == n` check to an always-active `BELFEM_ERROR` for release safety?

## Files Updated

- `todo/ai_exchange.md`
- `devlog/dl20260609_periodic_match_nodes_simplification_audit.md`
- `devlog/README.md`
