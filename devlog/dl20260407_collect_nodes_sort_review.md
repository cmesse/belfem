# Devlog 2026-04-07 — collect_nodes Sort Review

**Date:** 2026-04-07
**Topic:** Read-only review of sorting selected thin-shell source nodes by ID in `ThinShellFactory::collect_nodes()`
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Reviewed the new `sort( aNodes, opVertexID )` change at the end of `ThinShellFactory::collect_nodes()`.

I do not see a correctness bug introduced by the change. The updated local numbering remains internally consistent with later thin-shell routines that use `node->original()->index()` as their lookup key.

However, the new comment appears too strong. Sorting the selected node container by ID does not, by itself, preserve facet-local edge orientation. Edge objects are still canonicalized later, and local element orientation is recovered through explicit sign logic.

## Key Findings

- [cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp) now sorts the selected node container by `opVertexID` and then reassigns contiguous local indices.
- That change is internally consistent with:
  - node-normal assembly in `process_nodes_*()`
  - layer extrusion in `create_nodes_on_layers()`
  - node-based lookup during element/edge creation via `node->original()->index()`
- The actual temporary edge construction in [cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp) still canonicalizes edges later:
  - unordered edge key
  - final edge-node order by node ID
- The actual facet/element-local orientation still comes from explicit sign recovery in:
  - [cl_Mesh.cpp](/home/christian/codes/belfem/src/mesh/cl_Mesh.cpp)
  - [cl_FEM_Element.cpp](/home/christian/codes/belfem/src/fem/kernel/cl_FEM_Element.cpp)
- So the sort is likely behaviorally neutral for the current first-order edge-orientation question.

## Changes Made / Proposed

- No source-code changes made by Codex.
- Added an audit note to `todo/ai_exchange.md`.

## Open Questions

- Whether this sort is worth keeping as a determinism/readability aid, even if it is not the orientation fix.
- The layered patch-test failure still points more strongly toward the coupling path than toward local shell edge orientation.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260407_collect_nodes_sort_review.md
- devlog/README.md
