# Devlog 2026-04-06 — MUMPS `-40` Maxwell Matrix Audit

**Date:** 2026-04-06
**Topic:** Read-only audit of Maxwell matrix kernels after the reported MUMPS `-40` failure
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Reviewed the Maxwell matrix kernels in `src/fem/maxwell/matrices/` to look for other suspicious nonsymmetric or non-coercive contributions after the switch to symmetric Nitsche and symmetric MUMPS mode.

Conclusion: I do not see another high-confidence nonsymmetric Maxwell kernel. The current `h_ghost()` block is now algebraically symmetric. The strongest remaining matrix-side concern is that the active ghost penalty still uses the fully harmonic coefficient, which is likely too weak for the HTS/Hastelloy stack and can compromise coercivity of the symmetric Nitsche block.

## Key Findings

- MUMPS error `-40` is decoded locally as “The matrix was indicated to be positive definite but is not” in [src/sparse/cl_SolverMUMPS.cpp](/home/christian/codes/belfem/src/sparse/cl_SolverMUMPS.cpp#L816).
- `IWG_Maxwell` currently declares `SymmetryMode::GeneralSymmetric` in [src/fem/maxwell/cl_IWG_Maxwell.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_IWG_Maxwell.cpp#L40), so a true `-40` still suggests either a runtime mode mismatch or a definiteness issue rather than a plain asymmetry issue.
- The current ghost/Nitsche block in [src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L2142) is symmetric on inspection.
- The active ghost penalty coefficient is still
  `alpha = eta * rho_harm * ( 1./hm + 1./hs )`
  in [src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L2084). For the thin HTS plus thicker Hastelloy stack, this remains the most suspicious matrix-level cause of losing positive definiteness after moving to symmetric Nitsche.
- Volume `phi` and `h` kernels are assembled from symmetric Gram / rank-1 terms in [src/fem/maxwell/matrices/mt_maxwell_phi.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_phi.cpp#L25) and [src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L517).
- Symmetry and background-field kernels are also symmetric in [src/fem/maxwell/matrices/mt_maxwell_symmetry.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_symmetry.cpp#L23) and [src/fem/maxwell/matrices/mt_maxwell_background.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_background.cpp#L25).
- `phi_phi_2d()` still contains a debug `exit(0)` in [src/fem/maxwell/matrices/mt_maxwell_phi_phi.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_phi_phi.cpp#L23), but `InterfaceFerroAir` is deactivated for the `HPhi` path in [src/fem/maxwell/cl_MaxwellFactory.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_MaxwellFactory.cpp#L1197), so it is likely inactive in the present case.
- If ferromagnetic regions are active, the Newton correction assembled in [src/fem/iwg/cl_TimestepMatrices.cpp](/home/christian/codes/belfem/src/fem/iwg/cl_TimestepMatrices.cpp#L131) is symmetric but not guaranteed positive definite because it uses `dMdx_times_x - dMdx_times_h + dKdx_times_x * dt`.

## Changes Made / Proposed

- No source changes made.
- Logged the audit result in `todo/ai_exchange.md`.
- Added this devlog entry.

## Open Questions

- Whether the runtime that produced `-40` is definitely using `SYM=2` and not `SYM=1`.
- Whether the harmonic ghost penalty should now be replaced by `alpha = eta * max( rho_m / hm, rho_s / hs )` to restore coercivity on the HTS/Hastelloy interface.
- Whether a one-interface eigenvalue or patch test should be added to check that the symmetric ghost block is actually positive semidefinite / coercive for the chosen `alpha`.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260406_mumps_minus40_maxwell_audit.md
- devlog/README.md
