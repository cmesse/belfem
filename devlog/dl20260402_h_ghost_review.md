# Devlog 2026-04-02 — h_ghost Review

**Date:** 2026-04-02
**Topic:** Read-only review of the latest `h_ghost` edits
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** `todo/thinshell_selective_nitsche_coupling.md`

## Summary

Reviewed the current local changes in `src/fem/maxwell/matrices/mt_maxwell_h.cpp`. The edit removes the `aCalc->E(k)` crash on ghost sidesets, but the replacement `Dm` / `Ds` assembly is still inconsistent with the `PENTA6TS` face-local shell basis layout.

## Key Findings

- The ghost master facet is face `1` and the slave facet is face `0`. See `src/mesh/cl_ThinShellFactory.cpp:1738-1739`.
- For `PENTA6TS`, the top-face basis lives in edge columns `3..5`, while the bottom-face basis lives in columns `0..2`. See `src/fem/interpolation/nedelec/cl_EF_PENTA6TS.cpp:282-306`.
- The new `Dm` uses only `Em(:,0..2)` on the master top face, so it collapses to zero instead of representing a through-thickness derivative. See `src/fem/maxwell/matrices/mt_maxwell_h.cpp:1882-1894`.
- The new `Ds` reuses `Es(:,0..2)` for both local halves and with positive sign in both halves, so it does not match the intended `[-E | +E]/h` structure.
- The ghost coefficient logic still nulls the operator on HTS-adjacent interfaces when `T < T_crit`. See `src/fem/maxwell/matrices/mt_maxwell_h.cpp:1829-1840`.
- The follow-up sign revision
  `Dm(:,0:2)=-Em(:,0:2), Dm(:,3:5)=+Em(:,0:2), Ds(:,0:2)=+Es(:,0:2), Ds(:,3:5)=-Es(:,0:2)`
  fixes the sign duplication issue on `Ds`, but not the underlying data problem: `Em(:,0:2)` is still zero on the master top face, and `Es(:,0:2)` still represents only the slave bottom-face basis.

## Changes Made / Proposed

- Proposed only: none in source.
- Proposed only: re-derive the local derivative operator against the actual `PENTA6TS` column layout before continuing with this `Dm/Ds` path.

## Open Questions

- Whether the consistency term should continue through a local `E`-difference operator or be replaced entirely by a verified `Cm/Cs`-based formulation.
- If the `E`-difference formulation is kept, how best to expose `E_bot` and `E_top` for each shell element at the same interface `(xi, eta)` coordinates.
- The newest active-column version using `Em(:,3..5)` for the lower element and `Es(:,0..2)` for the upper element is plausible for linear `PENTA6TS`, but still worth verifying numerically against one affine element pair.
- The coefficient choice for HTS-adjacent interfaces remains open: `rho(T)` is a valid positive normal-state proxy for YBCO, but not the superconducting operating resistivity envisioned by the selective-DG design.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260402_h_ghost_review.md
