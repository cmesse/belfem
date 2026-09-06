# Devlog 2026-04-02 — TS / `h_ghost()` Compatibility Review

**Date:** 2026-04-02
**Topic:** Read-only review of Claude's implemented TS integration changes, focused on `h_ghost()` compatibility
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Reviewed the implemented thin-shell facet integration changes described in `todo/ai_exchange.md` and compared them against the current `h_ghost()` path. For the active ghost case (`PENTA6TS`, first-order 3D), the new TS facet-count, orientation-count, and master/slave integration-point logic appear compatible with the current ghost assembly.

## Key Findings

- The `PENTA6TS` ghost path now looks internally consistent:
  - ghost facets use master face `1` and slave face `0` in [src/mesh/cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L1696)
  - `Calculator::link( mesh::Facet * )` indexes the slave integration table as `index_on_slave * 3 + orientation_on_slave - 1` in [src/fem/kernel/cl_FEM_Calculator.cpp](/home/christian/codes/belfem/src/fem/kernel/cl_FEM_Calculator.cpp#L1014)
  - the new TS overrides in [src/mesh/meshtools.cpp](/home/christian/codes/belfem/src/mesh/meshtools.cpp#L563) and [src/mesh/meshtools.cpp](/home/christian/codes/belfem/src/mesh/meshtools.cpp#L803) make that indexing coherent for `PENTA*TS`
  - the resulting parametric points match the active-column expectations in [src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L2124)
- `QUAD*TS` slave orientation is still not fully plumbed through `Calculator::slave_integration_2d()` in [src/fem/kernel/cl_FEM_Calculator.cpp](/home/christian/codes/belfem/src/fem/kernel/cl_FEM_Calculator.cpp#L378), but this does not affect current `h_ghost()`.
- `PENTA18TS` still has an internal facet-count inconsistency between [src/mesh/cl_Element_PENTA18TS.hpp](/home/christian/codes/belfem/src/mesh/cl_Element_PENTA18TS.hpp#L46) and [src/mesh/meshtools.cpp](/home/christian/codes/belfem/src/mesh/meshtools.cpp#L563), but `h_ghost()` currently excludes higher-order TS elements.

## Changes Made / Proposed

- No source changes made.
- Appended a Codex review entry to `todo/ai_exchange.md`.
- Added this devlog entry.

## Open Questions

- Whether the repo wants to finish the analogous `QUAD*TS` slave-orientation plumbing even though it is not needed for the current 3D ghost path.
- Whether `PENTA18TS` should ultimately expose 2 physical facets or keep the third mid-plane facet with a separate TS-specific interpretation.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260402_ts_hghost_compatibility_review.md
