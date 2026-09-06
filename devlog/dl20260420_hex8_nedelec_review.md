# Devlog 2026-04-20 - HEX8 Nedelec Review

**Date:** 2026-04-20
**Topic:** Read-only review of `cl_EF_HEX8.*pp`
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** Local Nedelec implementation convention in `src/fem/interpolation/doc/nedelec.md` Section 5

## Summary

Reviewed the standard `HEX8` Nedelec edge-function implementation against BELFEM's edge-function contract, adjacent `HEX8TS` implementation, and `mesh::HEX8` edge ordering.

## Key Findings

- `EF_HEX8` is currently abstract because `cl_EF_HEX8.hpp` does not declare `E()` or `C()`, while `EdgeFunction::C()` is pure virtual.
- `cl_EF_HEX8.cpp` uses `InterpolationFunctionTemplate<HEX,LAGRANGE,3,8>` without including `lagrange/cl_IF_HEX8.hpp`.
- The implemented scalar factors and `E()` mapping follow the `HEX8TS` face-loop ordering, not the standard `mesh::HEX8` edge ordering.
- `C()` uses `mS[7]` instead of `mS[11]` for the `mEy(:,11)` entries.

## Changes Made / Proposed

- Appended the audit to `todo/ai_exchange.md`.
- No source-code edits were made.

## Open Questions

- After compile fixes, the element still needs an explicit 12x12 unit-circulation test on the reference hexahedron before being trusted in Maxwell assembly.

## Files Updated

- `todo/ai_exchange.md`
- `devlog/dl20260420_hex8_nedelec_review.md`
- `devlog/README.md`
