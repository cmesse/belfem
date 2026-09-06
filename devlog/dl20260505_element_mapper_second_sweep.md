# Devlog 2026-05-05 - ElementMapper Second Sweep

**Date:** 2026-05-05
**Topic:** Second read-only audit of `cl_IF_ElementMapper.*pp`
**AIs involved:** Codex
**Codex Audit Confidence:** high for remaining source-level issues, medium for intended embedded-surface scope
**Literature References:** N/A

## Summary

Re-audited `fem::ElementMapper` after the first-pass fixes. The original build integration,
dimension initialization, TET4 inverse, stale/null evaluator, QUAD coefficient layout, and copy
ownership blockers are addressed. A syntax-only compile of `cl_IF_ElementMapper.cpp` with the
current interpolation compile flags passed.

## Key Findings

- `evaluate_quad4()` still indexes `aXi(0)`/`aXi(1)` by reference in the non-affine branch before
  sizing the output vector.
- Non-affine QUAD failure cases can produce NaNs and still return true because neither the
  discriminant nor final parametric coordinates are checked for finiteness.
- `evaluate_general()` returns only the inside-reference check, even if Newton fails to converge
  within 100 iterations.
- The current implementation silently treats 2D elements with `set_dimension(3)` as projected
  x-y mappings; arbitrary embedded 3D surface elements are not correctly inverted.
- `inside_pyra()` still checks a full box in xi/eta instead of the shrinking pyramid cross-section.

## Changes Made / Proposed

- No source changes made.
- Proposed fixes: size `aXi` before all direct writes; reject negative QUAD discriminants and
  non-finite coordinates; require convergence in `evaluate_general()`; clarify or enforce the
  physical-dimension contract; tighten the pyramid inside test.

## Open Questions

- Should 2D elements embedded in 3D be supported by this class, or should `link()` require
  `mDim == aElement->dimension()` unless a caller deliberately asks for x-y projection?

## Files Updated

- `devlog/dl20260505_element_mapper_second_sweep.md`
- `devlog/README.md`
- `todo/ai_exchange.md`
