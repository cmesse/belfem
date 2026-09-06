# Devlog 2026-04-20 - HEX8 Nedelec Fix Review

**Date:** 2026-04-20
**Topic:** Read-only review of the user's `cl_EF_HEX8.*pp` fixes
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** Local Nedelec implementation convention in `src/fem/interpolation/doc/nedelec.md` Section 5

## Summary

Reviewed the user's fixes for the standard `HEX8` Nedelec edge function. The prior compile blockers are addressed and targeted syntax checks pass. One formulation defect remains: bottom-face edge basis signs are opposite to BELFEM's `mesh::HEX8` local edge directions.

## Key Findings

- `EF_HEX8` now declares `E()` and `C()`, addressing the abstract-class factory blocker.
- `cl_EF_HEX8.cpp` now includes the HEX8 interpolation template headers, addressing the missing-template blocker.
- The column-11 curl orientation typo is fixed.
- Bottom-face edge factors `F(0..3)` still integrate to `-1` along local `mesh::HEX8` edge directions; they and their derivatives should be sign-flipped.

## Changes Made / Proposed

- Appended the audit to `todo/ai_exchange.md`.
- No source-code edits were made.

## Open Questions

- A dedicated 12x12 reference-edge circulation test should be added once the bottom signs are corrected.

## Files Updated

- `todo/ai_exchange.md`
- `devlog/dl20260420_hex8_nedelec_fix_review.md`
- `devlog/README.md`
