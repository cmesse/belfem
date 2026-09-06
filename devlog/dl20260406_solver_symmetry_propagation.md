# Devlog 2026-04-06 — Solver Symmetry Propagation

**Date:** 2026-04-06
**Topic:** Read-only verification that `IWG_Maxwell` symmetry mode propagates into the solver / MUMPS path
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Checked whether the current `IWG_Maxwell` symmetry flag really reaches the sparse solver backend.

Conclusion: yes, the propagation is correct. `IWG_Maxwell` currently declares `SymmetryMode::Unsymmetric`, `SolverData` forwards that declaration into the `Solver`, and the MUMPS bridge receives `SYM=0`.

## Key Findings

- `IWG_Maxwell` currently uses `SymmetryMode::Unsymmetric` in [src/fem/maxwell/cl_IWG_Maxwell.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_IWG_Maxwell.cpp#L36).
- `MaxwellFactory::configure_solver()` passes solver parameters to the field in [src/fem/maxwell/cl_MaxwellFactory.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_MaxwellFactory.cpp#L2189).
- `DofManager::set_solver()` forwards to `SolverData::set_solver()` in [src/fem/kernel/cl_FEM_DofManager.cpp](/home/christian/codes/belfem/src/fem/kernel/cl_FEM_DofManager.cpp#L1015).
- `SolverData::set_solver()` creates the solver and immediately copies `mParent->iwg()->symmetry_mode()` into it in [src/fem/kernel/cl_FEM_DofMgr_SolverData.cpp](/home/christian/codes/belfem/src/fem/kernel/cl_FEM_DofMgr_SolverData.cpp#L2023).
- On first solve, `Solver::solve()` initializes the backend wrapper with that stored symmetry mode in [src/sparse/cl_Solver.cpp](/home/christian/codes/belfem/src/sparse/cl_Solver.cpp#L132).
- `MUMPS::initialize()` forwards the mode into the MUMPS parameter block in [src/sparse/cl_SolverMUMPS.cpp](/home/christian/codes/belfem/src/sparse/cl_SolverMUMPS.cpp#L127).
- The Fortran bridge sets `tMUMPS%SYM = tSymmetryMode` and documents the mapping `SYM=0 unsymmetric`, `SYM=1 SPD`, `SYM=2 general symmetric` in [src/sparse/mumpstools.f90](/home/christian/codes/belfem/src/sparse/mumpstools.f90#L256).

## Changes Made / Proposed

- No source changes made.
- Logged the propagation check in `todo/ai_exchange.md`.
- Added this devlog entry.

## Open Questions

- None on propagation itself.
- Separate question remains whether the assembled Maxwell matrix is actually symmetric when ghost/Nitsche interfaces are active.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260406_solver_symmetry_propagation.md
- devlog/README.md
