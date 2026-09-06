# Devlog 2026-04-07 — `compute_bn()` Helper Review

**Date:** 2026-04-07
**Topic:** Read-only review of the new thin-shell `compute_bn()` helper
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Reviewed the new `compute_bn()` helper added to `src/fem/maxwell/matrices/mt_maxwell_h.hpp`.

The helper body is correct for the intended fix: it computes the averaged air-side field from `phi_m` / `phi_s`, projects it onto the facet normal, and returns the projected `bn`.

However, the helper is not yet used by any thin-shell kernel, so it currently has no effect and does not fix the original `bn` bug.

## Key Findings

- `compute_bn()` is declared and defined in [mt_maxwell_h.hpp](../src/fem/maxwell/matrices/mt_maxwell_h.hpp), but the active shell kernels in [mt_maxwell_h.cpp](../src/fem/maxwell/matrices/mt_maxwell_h.cpp) still use the old pattern
  `bn = tCalc->Bm(0) * phi_m + tCalc->Bs(0) * phi_s;`
  followed by `bn *= -0.5 * mu0;`.

- The helper implementation itself is sound:
  - it uses `Bm/Bs` to build full `bm` and `bs`,
  - averages them,
  - projects the result onto `tCalc->normal(k)`.

- The scratch vectors `bm` and `bs` are already provisioned in [cl_IWG_Maxwell.cpp](../src/fem/maxwell/cl_IWG_Maxwell.cpp), so using the helper should not fail for missing calculator storage.

- Minor API concern:
  - `compute_bn( aCalc, k )` suggests per-quadrature evaluation,
  - while the current shell kernels explicitly assume the flat/linear case and compute `bn` once outside the loop.
  - That is not a correctness bug by itself, but the intended usage should be made explicit.

## Changes Made / Proposed

- No source-code changes were made in this review.
- Added a Codex audit entry to `todo/ai_exchange.md`.

## Open Questions

- None about the helper body itself.
- The remaining required step is simply to wire the helper into the shell kernels.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260407_compute_bn_helper_review.md
