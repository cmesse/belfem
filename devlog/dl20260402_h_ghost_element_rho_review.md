# Devlog 2026-04-02 — `h_ghost()` Element-Rho Review

**Date:** 2026-04-02
**Topic:** Read-only review of the new `element_rho` cache path and latest `h_ghost()` logic
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Reviewed the latest `mt_maxwell_h.cpp` changes after the introduction of the `element_rho` mesh field. The cache approach is a better architectural direction than trying to reconstruct full HTS resistivity logic inside `h_ghost()`, but the current `h_ghost()` implementation still contains two high-confidence logic errors.

## Key Findings

- The new interface coefficient in [src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L2080) is mathematically wrong: `(rho_m + rho_s) / (rho_m - rho_s)` is not a harmonic resistivity and becomes singular when `rho_m == rho_s`.
- The current `Dm` / `Ds` assembly in [src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L2132) and [src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L2137) fills both local halves with the same sign, so it no longer represents a through-thickness difference operator.
- The new cache plumbing itself looks coherent: `element_rho` is registered in [src/fem/maxwell/cl_IWG_Maxwell.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_IWG_Maxwell.cpp#L75), resolves to an element field in [src/mesh/fn_entity_type.cpp](/home/christian/codes/belfem/src/mesh/fn_entity_type.cpp#L48), and block assembly precedes sideset assembly in [src/fem/kernel/cl_FEM_DofManager.cpp](/home/christian/codes/belfem/src/fem/kernel/cl_FEM_DofManager.cpp#L631).
- The cached `element_rho` is an element-average proxy, not an exact interface resistivity. That is acceptable if documented as an approximate stabilization coefficient.

## Changes Made / Proposed

- No source changes made.
- Added collaboration notes to `todo/ai_exchange.md`.
- Added this devlog entry to record the review outcome.

## Open Questions

- Whether the same-sign `Dm` / `Ds` fill was an intentional reformulation or an accidental sign error.
- Whether `h_ghost()` should use the standard positive interface coefficient directly once cached `element_rho` values are guaranteed to be positive.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260402_h_ghost_element_rho_review.md
