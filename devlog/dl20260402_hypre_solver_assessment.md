# Devlog 2026-04-02 — hypre Solver Assessment

**Date:** 2026-04-02
**Topic:** Read-only assessment of whether hypre would be a worthwhile addition to BELFEM's current solver stack
**AIs involved:** Codex
**Codex Audit Confidence:** medium-high (~80%)
**Literature References:** Messe et al. 2023 (paper1) as summarized in repo docs; PETSc `PCHYPRE` manual page; hypre BoomerAMG / AMS / ADS / ILU / Euclid docs

## Summary

Evaluated hypre against BELFEM's present solver architecture and documented solver guidance. Conclusion: hypre is a plausible PETSc-side complement for elliptic subproblems, but not a compelling new primary backend for BELFEM's hard Maxwell/HTS cases, where STRUMPACK and MUMPS remain the better fit for robustness.

## Key Findings

- BELFEM currently positions STRUMPACK as first choice and MUMPS as robustness fallback for difficult production solves, with PETSc reserved for very large / better-conditioned iterative cases.
- The PETSc wrapper does not currently expose `PCHYPRE`; the preconditioner enum lacks a `HYPRE` value and the wrapper sets the PC type from that enum.
- BELFEM's Maxwell systems are mixed H(curl)/H1 systems (`edge_h`, `phi`, and sometimes `face_h`), so the attractive hypre components are the more specialized ones, not just generic BoomerAMG.
- hypre AMS/ADS require auxiliary operators and geometry data that BELFEM does not currently pass through its PETSc wrapper.
- For the stated scale ceiling (~10M DOFs), robustness under ill-conditioning is still likely to dominate over extreme-scale AMG scalability.

## Changes Made / Proposed

- Added a read-only assessment entry to `todo/ai_exchange.md`.
- No source files modified.

## Open Questions

- Are the dominant user workloads memory-bound in STRUMPACK/MUMPS, or primarily convergence/robustness-bound?
- Is there interest in a narrow PETSc+BoomerAMG path for Poisson/thermal/homology only, rather than a general hypre backend?
- Are future Maxwell iterations expected to depend on low-order H(curl) iterative solves strongly enough to justify AMS integration work?

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260402_hypre_solver_assessment.md
