# Devlog 2026-05-05 - ElementMapper Audit

**Date:** 2026-05-05
**Topic:** Read-only audit of `cl_IF_ElementMapper.*pp`
**AIs involved:** Codex
**Codex Audit Confidence:** high for source-level defects, medium for intended API assumptions
**Literature References:** N/A

## Summary

Audited the new `fem::ElementMapper` class in `src/fem/interpolation/cl_IF_ElementMapper.*pp`.
No source changes were made. A syntax-only compile of the `.cpp` with the current interpolation
compile flags passed, but the file is not currently part of the interpolation CMake source list.

## Key Findings

- `cl_IF_ElementMapper.cpp` is not listed in `src/fem/interpolation/CMakeLists.txt`, so the class
  will not be built or linked by the current library target.
- The default `mDim = 3` prevents `link()` from calling `set_dimension()`, leaving `mJ` and `mInvJ`
  unsized before fast-path evaluators write into them.
- Higher-order/general element types never set `mFunEval` to `evaluate_general()`, so `evaluate()`
  can call a null or stale function pointer.
- The QUAD4 path has inconsistent `mA` dimensions/layout and writes outside the allocated matrix;
  the affine and non-affine inverse formulas also use inconsistent coefficient indexing.
- The TRI3 path has a Jacobian entry typo mixing node-1 `y` with node-2 `x`.
- The TET4 path calls `inv2()` on a 3x3 Jacobian instead of `inv3()`.
- The Newton/general path does not recompute `N` after each update, overwrites `aXi` instead of
  applying an incremental correction, and always uses `inv3()` even for 2D mappings.
- `ElementMapper` owns raw pointers but has implicit copy/assignment, which risks double deletion
  if the class is copied.

## Changes Made / Proposed

- No source changes made.
- Proposed fixes: add the file to CMake, initialize/synchronize `mDim` and Jacobian storage in
  `link()`, repair the linear fast-path formulas, set `mFunEval = &ElementMapper::evaluate_general`
  for default element families, and rewrite the Newton loop with dimension-specific linear solves.

## Open Questions

- Should this mapper support embedded 2D elements in 3D coordinates, or only elements whose
  topological dimension matches the physical coordinate dimension?
- Should degenerate/singular element inversion be a debug assertion only, or an always-active
  `BELFEM_ERROR` because the mapper may be used in production search/mapping paths?

## Files Updated

- `devlog/dl20260505_element_mapper_audit.md`
- `devlog/README.md`
- `todo/ai_exchange.md`
