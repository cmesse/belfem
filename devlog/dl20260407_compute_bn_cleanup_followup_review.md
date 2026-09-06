# Devlog 2026-04-07 — compute_bn Cleanup Follow-up Review

**Date:** 2026-04-07
**Topic:** Read-only review of the latest thin-shell `compute_bn()` cleanup after wiring and build fixes
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Reviewed the current source after the new `compute_bn()` wiring, calculator-vector additions, and `get_thermal_calculator()` cleanup. The normal-field projection fix is now actually compiled into the thin-shell kernels, and I did not find a new correctness bug in this cleanup itself.

## Key Findings

- `compute_bn()` is now called from every thin-shell shell kernel family in [src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp).
- The helper prerequisites are present in [src/fem/maxwell/cl_IWG_Maxwell.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_IWG_Maxwell.cpp): vectors `n`, `bm`, and `bs` are now created.
- The linear-element cache guard is active: `mIsLinear` exists in [src/fem/kernel/cl_FEM_Calculator.hpp](/home/christian/codes/belfem/src/fem/kernel/cl_FEM_Calculator.hpp) and is assigned in [src/fem/kernel/cl_FEM_Calculator.cpp](/home/christian/codes/belfem/src/fem/kernel/cl_FEM_Calculator.cpp).
- `get_thermal_calculator()` is now centralized and the previously dangling local thermal-calculator cleanup problem is gone.
- Residual risk remains in dataflow style only: callers still rely on `compute_bn()` being called before reading `aCalc->vector("n")`.
- This cleanup does not explain the Cu/Ag multilayer patch-test failure; that investigation remains open.

## Changes Made / Proposed

- No source changes made by Codex in this session.
- Added audit notes to `todo/ai_exchange.md`.

## Open Questions

- Why the Cu/Ag multilayer case still sits around `20 dB` despite the normal-field cleanup.
- Whether the remaining patch-test defect is in thin-shell topology / interface attachment, or in the inter-layer coupling algebra.

## Files Updated

- devlog/dl20260407_compute_bn_cleanup_followup_review.md
- devlog/README.md
- todo/ai_exchange.md
