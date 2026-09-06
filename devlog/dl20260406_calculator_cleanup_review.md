# Devlog 2026-04-06 — Calculator Cleanup Review

**Date:** 2026-04-06
**Topic:** Read-only review of the latest `cl_FEM_Calculator.cpp` cleanup and prism ghost changes
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Reviewed the current `cl_FEM_Calculator.cpp` changes as a code review. Two blockers remain: the file no longer builds under `-Werror` because the A-cleanup is incomplete, and the new `slave_integration_penta()` indexing is wrong for `PENTA6TS`/`PENTA18TS`, which is the active ghost thin-shell path.

## Key Findings

- `allocate_memory()` still sets `tHaveA`, but the new logic no longer uses it, so `cmake --build cmake-build-debug --target hphirun -j4` fails with `-Werror=unused-but-set-variable` on `tHaveA` (`src/fem/kernel/cl_FEM_Calculator.cpp:208`).
- The A-formulation cleanup is only partial: `mFunNedelecDataA`, `nedelec_data_a()`, and the `"edge_a"` / `"face_a"` helpers are still declared and defined in `src/fem/kernel/cl_FEM_Calculator.hpp:268` and `src/fem/kernel/cl_FEM_Calculator.hpp:1692`.
- The new `slave_integration_penta()` helper uses full-prism cumulative offsets `{0,4,8,12,15}` (`src/fem/kernel/cl_FEM_Calculator.cpp:418`).
- Thin-shell prisms use only 2 facets with 3 orientations each (`src/mesh/meshtools.cpp:804`), and `SideSet` stores slave integration tables in facet-major dense order based on those counts (`src/fem/kernel/cl_FEM_SideSet.cpp:515`).
- As a result, `slave_integration_penta()` misindexes `PENTA6TS` / `PENTA18TS` slave integration data, including a possible out-of-range access on facet 1 orientation 3.

## Changes Made / Proposed

- No source changes made.
- Appended the review result to `todo/ai_exchange.md`.
- Proposed immediate fixes:
  - finish removing the dead A-path state or keep using `tHaveA` until the cleanup is complete,
  - split `slave_integration_penta()` into TS vs non-TS handling, or compute the offset from `mesh::number_of_orientations( mGroup->slave_type(), f )` instead of hard-coding full-prism counts.

## Open Questions

- Whether the same “full prism vs thin-shell prism” distinction should also be made explicit in any future slave-face or orientation helper added for `PENTA18TS`.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260406_calculator_cleanup_review.md
