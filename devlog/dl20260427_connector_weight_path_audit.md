# Devlog 2026-04-27 — Connector Weight Path Audit

**Date:** 2026-04-27
**Topic:** Read-only audit of the new weighted side-connector path in `ThinShellFactory` and `cl_FEM_DofMgr_DofData`
**AIs involved:** Codex
**Codex Audit Confidence:** high

## Summary

Investigated whether the bad result after adding weighted connector-edge sources is more likely caused by connection mechanics or by postprocessing.

Verdict:

- The stronger suspect is the **connection mechanics**.
- The current `preprocess_binomial_edges()` implementation computes source-edge weights but does not attach them to the actual connector horizontal edges.
- The new multi-source edge branch in `cl_FEM_DofMgr_DofData.cpp` therefore appears not to be exercised yet for this path.
- The current postprocessor wiring for connector blocks looks consistent with the earlier fix and is not the first place I would blame.

## Key Findings

- `preprocess_binomial_edges()` computes `tSources` and `tWeights` but ends without any `set_sources(...)` call.
  - `src/mesh/cl_ThinShellFactory.cpp:2851-2929`

- The target edges that should carry those sources, `EdgesHorizontal`, are only created later by `create_binomial_edges()`.
  - `src/mesh/cl_ThinShellFactory.cpp:2037-2049`
  - `src/mesh/cl_ThinShellFactory.cpp:2530-2557`

- The source-edge selection predicate in `preprocess_binomial_edges()` has suspicious operator precedence:
  - current code: `tA->is_flagged() xor tB->is_flagged() && (...)`
  - likely intended: `( tA->is_flagged() xor tB->is_flagged() ) && (...)`
  - `src/mesh/cl_ThinShellFactory.cpp:2884`

- Connector blocks are now passed to the postprocessor through the dedicated connector-block collection.
  - `src/fem/maxwell/cl_MaxwellFactory.cpp:1894-2040`

- `MaxwellPostprocessor::Conductor` accepts `LeftCoating` / `RightCoating` and uses `compute_conductor()`.
  - `src/fem/maxwell/cl_MaxwellPostprocessor.cpp:209-235`
  - `src/fem/maxwell/cl_MaxwellPostprocessor.cpp:462-476`

## Changes Made / Proposed

- No source-code changes made.
- Added AI exchange note and this devlog.

## Open Questions

- The next direct check is whether `EdgesHorizontal(k)->number_of_sources()` remains zero after creation. If yes, the weighted DOF code path is currently dead.
- Once sources are actually attached, re-check whether the weighted combination should collapse onto shell-edge DOFs directly or onto already-hanging edge DOFs from duplicated inner edges.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260427_connector_weight_path_audit.md
- devlog/README.md
