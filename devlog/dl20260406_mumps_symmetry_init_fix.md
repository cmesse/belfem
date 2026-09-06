# Devlog 2026-04-06 — MUMPS Symmetry Init Fix

**Date:** 2026-04-06
**Topic:** Initialize MUMPS with the requested symmetry mode before `JOB=-1`
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

The solver-side symmetry flag propagation was already correct, but the MUMPS wrapper created the solver instance too early. `mumpstools_create_solver()` called `DMUMPS` with `JOB=-1` before `PAR` and `SYM` were assigned, so the instance started in the default unsymmetric mode even when the Maxwell path requested `SYM=1`.

The fix was to pass `WorkingHost` and `SymmetryMode` into the creation routine and assign both before the `JOB=-1` initialization call. `hphirun` rebuilt successfully after the change.

## Key Findings

- `MUMPS::initialize()` now forwards `WorkingHost` and `SymmetryMode` into the solver-creation bridge in [src/sparse/cl_SolverMUMPS.cpp](/home/christian/codes/belfem/src/sparse/cl_SolverMUMPS.cpp#L127).
- The C binding for `mumpstools_create_solver()` now accepts those two inputs in [src/sparse/mumpstools.hpp](/home/christian/codes/belfem/src/sparse/mumpstools.hpp#L33).
- The Fortran bridge now sets `tMUMPS%PAR` and `tMUMPS%SYM` before `tMUMPS%JOB = -1` and `call DMUMPS( tMUMPS )` in [src/sparse/mumpstools.f90](/home/christian/codes/belfem/src/sparse/mumpstools.f90#L58).
- The earlier user-facing symptom, `Symmetry 1` in debug output plus an unsymmetric MUMPS banner, is consistent with the old init order and no longer points to a propagation bug.

## Changes Made / Proposed

- Updated [src/sparse/mumpstools.hpp](/home/christian/codes/belfem/src/sparse/mumpstools.hpp) to extend the `mumpstools_create_solver()` signature.
- Updated [src/sparse/cl_SolverMUMPS.cpp](/home/christian/codes/belfem/src/sparse/cl_SolverMUMPS.cpp) to pass `WorkingHost` and `SymmetryMode` when creating the MUMPS solver instance.
- Updated [src/sparse/mumpstools.f90](/home/christian/codes/belfem/src/sparse/mumpstools.f90) so `PAR` and `SYM` are assigned before the `JOB=-1` initialization call.

## Open Questions

- A fresh Maxwell run is still needed to confirm that MUMPS now prints the symmetric banner instead of the unsymmetric LU banner.
- Temporary debug prints can be removed after runtime confirmation.

## Files Updated

- src/sparse/mumpstools.hpp
- src/sparse/cl_SolverMUMPS.cpp
- src/sparse/mumpstools.f90
- todo/ai_exchange.md
- devlog/dl20260406_mumps_symmetry_init_fix.md
- devlog/README.md
