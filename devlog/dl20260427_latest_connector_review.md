# Devlog 2026-04-27 — Latest Connector Review

**Date:** 2026-04-27
**Topic:** Read-only review of the latest weighted side-connector changes
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Re-reviewed the current working tree after the weighted binomial-edge hookup was extended into `ThinShellFactory` and `DofData`. The earlier verdict that the new DOF-manager branch was dead no longer applies, because `create_binomial_edges()` now attaches multi-source edge bases. The current primary risks are still on the mesh / connection side rather than in the postprocessor.

## Key Findings

- `src/mesh/cl_ThinShellFactory.cpp:2913` has a precedence bug in the source-edge predicate: `xor` is evaluated after `&&`, so the code does not test `(flagged xor flagged) && incident` as intended.
- `src/mesh/cl_ThinShellFactory.cpp:2025`, `:2081-2086`, and `:343` now create and append null side-connector `SideSet*` entries into the mesh because side-facet generation is commented out but `aExtraSideSets` is still sized and appended.
- `src/mesh/cl_Mesh.cpp:644-666` and `src/fem/maxwell/cl_MaxwellFactory.cpp:805-814` both dereference mesh sidesets unconditionally, so those null connector sidesets are not benign.
- `src/mesh/cl_ThinShellFactory.cpp:2555` uses `aLayers( j )->Edges` instead of the `mIndices`-selected physical layer; this is a latent wrong-layer binding if selective connector layers return.
- The connector postprocessor path still looks coherent: `LeftCoating` / `RightCoating` blocks are collected and sent through the normal conductor postprocessor path (`src/fem/maxwell/cl_MaxwellFactory.cpp:1894-2040`, `src/fem/maxwell/cl_MaxwellPostprocessor.cpp:209-235`, `:462-476`).

## Changes Made / Proposed

- Added a current audit note to `todo/ai_exchange.md`
- Added this devlog and updated `devlog/README.md`
- Proposed prioritizing the mesh-side fixes before deeper postprocessor debugging

## Open Questions

- Are the connector interface sidesets intentionally disabled for the current experiment, or should the old side-facet generation path be restored?
- Will `mConnectorsForAllLayers` remain hard-wired to true, or should the binomial-edge source lookup be corrected now to preserve the selective-layer mode?

## Files Updated

- `todo/ai_exchange.md`
- `devlog/dl20260427_latest_connector_review.md`
- `devlog/README.md`
