# Devlog 2026-04-02 — `h_ghost()` Follow-Up Review

**Date:** 2026-04-02
**Topic:** Read-only follow-up review of the latest `h_ghost()` fixes
**AIs involved:** Codex
**Codex Audit Confidence:** medium-high (~85%)
**Literature References:** N/A

## Summary

Reviewed the latest revision of `h_ghost()` in `mt_maxwell_h.cpp` after the fixes to the interface coefficient and `Dm` / `Ds` sign pattern. Within the stated scope of first-order 3D only and a hard-coded `eta`, I do not see a new high-confidence logic defect in the current implementation.

## Key Findings

- The interface coefficient in [src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L2078) is back to the expected `rho_m * rho_s / ( rho_m + rho_s )` form.
- The current `Dm` / `Ds` fill in [src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L2131) now has the expected opposite signs between local bottom/top slots.
- The active-face mapping is consistent with ghost creation in [src/mesh/cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L1696) and with the `PENTA6TS` upper/lower face operator split in [src/fem/interpolation/nedelec/cl_EF_PENTA6TS.cpp](/home/christian/codes/belfem/src/fem/interpolation/nedelec/cl_EF_PENTA6TS.cpp#L282).
- Remaining caveat: the construction still uses an active-face proxy for the local two-face derivative operator. That seems acceptable for the current linear shell element, but a small numerical check would still be prudent.
- Remaining caveat: `element_rho` is an element-average stabilization proxy, not an exact interface resistivity.

## Changes Made / Proposed

- No source changes made.
- Added collaboration notes to `todo/ai_exchange.md`.
- Added this devlog entry.

## Open Questions

- Whether a one-element / one-interface patch test should be added to lock in the `Dm` / `Ds` convention numerically.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260402_h_ghost_followup_review.md
