# Devlog 2026-04-02 — ThinShell Selective Duplication Follow-Up Review

**Date:** 2026-04-02
**Topic:** Read-only follow-up audit of the latest `ThinShellFactory` selective duplication fixes
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Reviewed the latest fixes to the `ThinShellFactory` selective duplication work.

Most of the earlier build-breaking leftovers are now fixed. Two issues remain: the first-order edge linker still does not advance the layer index, and the new ghost suppression path sizes `GhostFacets` before the `hasDuplicates` early exit, which leaves null ghost entries on shared interfaces.

## Key Findings

- The previous helper-signature mismatch and removed-member references in `create()` are fixed ([src/mesh/cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L210), [src/mesh/cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L282)).
- The first-order branch of `link_elements_with_edges()` still never advances `l`, so all first-order blocks would reuse the same two layer edge containers ([src/mesh/cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L1558)).
- `create_ghost_facets()` now guards on `hasDuplicates`, but it still executes `tFacets.set_size( aFacets.size(), nullptr )` before the guard. That leaves nonzero `GhostFacets` containers containing only null pointers on same-material interfaces ([src/mesh/cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L1678)).
- `create()` later counts those null entries via `GhostFacets.size()` and iterates them unconditionally when constructing the ghost sideset, which can still dereference null pointers ([src/mesh/cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L298), [src/mesh/cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L312)).
- The `hasDuplicates` layer-marking sequence still appears correct for both first- and second-order shells ([src/mesh/cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L186)).

## Changes Made / Proposed

- No source changes made.
- Logged the follow-up audit in `todo/ai_exchange.md`.

## Open Questions

- Should `create_ghost_facets()` simply skip untouched layers entirely, or should it explicitly zero out `GhostFacets` on non-duplicated interfaces for defensive clarity?
- Are first-order shells exercised in the current Maxwell workflow, or is the remaining `link_elements_with_edges()` defect latent for now?

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260402_thinshell_selective_duplication_followup_review.md
