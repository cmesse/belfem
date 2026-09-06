# Devlog 2026-04-06 — Symmetric Nitsche Audit

**Date:** 2026-04-06
**Topic:** Read-only audit of the new symmetric `h_ghost()` assembly and penalty scaling
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Audited the current `maxwell::h_ghost()` implementation after the recent changes intended to make the ghost/Nitsche matrix symmetric and to revise the penalty scaling.

Conclusion: the current implementation is still unsymmetric. The diagonal adjoint-consistency terms were added correctly, but the off-diagonal adjoint terms were inserted into the wrong blocks. Separately, the current `alpha = eta * max( rho_m / hm, rho_s / hs )` penalty scaling is a better practical choice for the HTS/Hastelloy stack than the previous fully harmonic penalty.

## Key Findings

- `IWG_Maxwell` now advertises `SymmetryMode::GeneralSymmetric` in [src/fem/maxwell/cl_IWG_Maxwell.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_IWG_Maxwell.cpp#L40).
- In [src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L2147), `Kmm` and `Kss` include the expected self-adjoint additions `-Dm^T Em` and `+Ds^T Es`.
- The off-diagonal adjoint terms are currently swapped:
  - `Kms` includes `-Ds^T Em`
  - `Ksm` includes `+Dm^T Es`
- Symmetry requires the opposite placement:
  - `Kms` should receive `+Dm^T Es`
  - `Ksm` should receive `-Ds^T Em`
- The current source now uses `alpha = eta * max( rho_m / hm, rho_s / hs )` in [src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L2086), while retaining `rho_harm` in the consistency weighting `rwdS`. For the thin HTS plus thicker Hastelloy stack, that split is more suitable than a fully harmonic penalty because it does not let the continuity penalty collapse to the HTS resistivity.

## Changes Made / Proposed

- No source changes made.
- Logged the audit result in `todo/ai_exchange.md`.
- Added this devlog entry.

## Open Questions

- Whether to patch the off-diagonal adjoint terms immediately and keep `SymmetryMode::GeneralSymmetric`.
- Whether to expose `eta` through `psi()` or a similar runtime parameter now that the penalty scaling is being revisited.
- Whether a small one-interface patch test should be added to check `Kms = Ksm^T` numerically.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260406_symmetric_nitsche_audit.md
- devlog/README.md
