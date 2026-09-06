# Devlog 2026-04-02 — ThinShell Selective Duplication Review

**Date:** 2026-04-02
**Topic:** Read-only audit of the new `ThinShellFactory` selective edge/face duplication logic
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Reviewed the in-progress `ThinShellFactory` changes that introduce per-layer `hasDuplicates` and shared-vs-duplicated edge/face containers keyed off adjacent material labels.

The interface-layer indexing for `hasDuplicates` looks conceptually correct, but the implementation still contains several build-breaking leftovers in `create()`, one first-order layer-linking bug, and one unresolved semantic issue: ghost facets are still created on all interfaces even when the interface is intentionally kept continuous.

## Key Findings

- `ThinShellFactory::create()` still calls `create_edges_on_layers()` / `create_faces_on_layers()` with the removed `aOrder` argument ([src/mesh/cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L211), [src/mesh/cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L233)).
- `ThinShellFactory::create()` still appends removed layer members `EdgesBottom/EdgesMid/EdgesTop` and `FacesBottom/FacesMid/FacesTop` after the `Layer` struct was changed to `Edges/EdgeDuplicates/Faces/FaceDuplicates` ([src/mesh/cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L282), [src/mesh/cl_ThinShellFactory.hpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.hpp#L35)).
- `ThinShellFactory::create()` redeclares `tCount` in the same scope ([src/mesh/cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L189), [src/mesh/cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L261)).
- First-order `link_elements_with_edges()` no longer advances `l`, so multiple blocks would reuse the same two layer edge containers ([src/mesh/cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L1567)).
- `hasDuplicates` is set on the expected interface layers for both first- and second-order shell layouts ([src/mesh/cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L186)).
- `create_ghost_facets()` still creates ghost sidesets across all interfaces, so same-material interfaces do not yet cleanly collapse back to CG behavior ([src/mesh/cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L275), [src/mesh/cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L1665), [src/fem/maxwell/cl_MaxwellFactory.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_MaxwellFactory.cpp#L1712)).

## Changes Made / Proposed

- No source changes made.
- Logged the audit in `todo/ai_exchange.md`.

## Open Questions

- Should same-material interfaces suppress ghost facet creation entirely, or is the intent to retain a stabilization term even when edge/face DOFs are shared?
- If material equality is the trigger, is exact string equality sufficient, or should the decision be based on a resolved material identity / material class?

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260402_thinshell_selective_duplication_review.md
