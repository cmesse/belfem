# Devlog 2026-05-06 - `map_extra_nodes()` Flag-0 Reuse

**Date:** 2026-05-06
**Topic:** Read-only investigation of `ThinShellFactory::map_extra_nodes()` second-call BFS failure
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Investigated the reported second-call failure in `ThinShellFactory::map_extra_nodes()`. The diagnostic output is consistent with stale element flag 0, not memory corruption: the first BFS marks reached air/ferro elements with flag 0 and the function never clears that flag before returning or before the next call.

## Key Findings

- `ThinShellFactory` constructor initially clears element flag 0 and marks air/ferro elements with flag 1 at `src/fem/kernel/cl_ThinShellFactory.cpp:45-59`.
- `map_extra_nodes()` resets element levels at `src/fem/kernel/cl_ThinShellFactory.cpp:3132-3137`, but does not clear flag 0.
- The BFS uses flag 0 as its visited marker at `src/fem/kernel/cl_ThinShellFactory.cpp:3148-3150` and `src/fem/kernel/cl_ThinShellFactory.cpp:3179-3183`.
- On the next side-curve call, the stale visited flags suppress BFS expansion. This matches the second diagnostic: one seed level exists, `level 1 : 0`, and only nodes inside the seed elements are found before the level-count guard fires.
- Secondary issue: `tNumNodes` is taken from `aCurve->nodes().size()` at `src/fem/kernel/cl_ThinShellFactory.cpp:3203` even though the loop is placing `aNodes`; this is harmless only while both counts match.

## Changes Made / Proposed

- No source changes made.
- Proposed fix: clear element flag 0 at the start of each `map_extra_nodes()` call, before seeding the BFS, or avoid mesh-global flag 0 and track visited state locally.

## Open Questions

- Whether the preferred fix should preserve any pre-existing caller-owned element flag 0 state. If yes, a local `DynamicBitset`/level sentinel is safer than `mMesh->unflag_all_elements(0)`.

## Files Updated

- `devlog/dl20260506_map_extra_nodes_flag0.md`
- `devlog/README.md`
