# Devlog 2026-09-16 — pardisotools: MKL's own interface, no dparm, no determinant

**Date:** 2026-09-16
**Topic:** The two gfortran argument-mismatch warnings in `src/sparse/pardisotools.f90`, and what an explicit interface exposed
**AIs involved:** Claude (Fable 5.1), Codex (gpt-5.6-terra, high), Grok (grok-4.6, high) — blind jury on the proposal before any edit
**Claude Confidence:** high (interface, dparm), medium (ILP64 finding — needs an ilp64 build to observe)
**Codex Audit Confidence:** high
**Literature References:** Intel oneMKL `mkl_pardiso.f90` (vendor header, sixteen-argument `PARDISO_D`), oneMKL Fortran developer reference, `pardiso iparm` parameter table (iparm( 32 )–( 33 ) reserved)
**Verification:** focused regression — `pardisotools.f90` compiles without any argument-mismatch warning against MKL's interface; `test_sparse` `SolverSolve.PARDISO*` 2/2 (tridiagonal solve through phases 11/22/33 and the release phase, determinant throws); `make check` 17/17 on `cmake-build-debug` (debug, MKL lp64), 2026-09-16

## Summary

Christian's test build in `./build` showed two `-Wargument-mismatch` warnings in `pardisotools.f90`: argument 7
of `pardiso` is the rank-1 `aValues` in the factorizations but a scalar dummy in `pardisotools_free`; argument 14 is
the rank-2 `aRHS` in the solve but a scalar elsewhere. Cause: `pardiso` had no explicit interface, so gfortran could
only compare the call sites with each other. gfortran has no `!GCC$ diagnostic push`, so the `stringtools` wrapper
pattern does not transfer; Christian's direction was to do what `mumpstools.f90` does — include the vendor header.

Writing the interface exposed three things the implicit interface had hidden: (1) the phase-33 call passed a
**seventeenth** argument, `gDPARM`, which is the Panua/Schenk PARDISO API and does not exist in MKL's
sixteen-argument `pardiso`; (2) the determinant path (`pardisotools_get_determinant` = `gDPARM( 33 )`) could
therefore never be written by MKL — on the live path the C++ flag is hard-wired off and `get_determinant()` returned
a signaling NaN (Codex corrected my "silent 0.0"); oneMKL documents `iparm( 33 )` as reserved and required zero;
(3) Christian's caution that a Fortran `INTEGER`'s width follows the compiler flags, not the source, turned into a
finding: `USE_MKL_64BIT_API` links MKL ilp64 and makes `int_t` 64-bit, but the module's plain `integer` variables
stayed 32-bit and no `-fdefault-integer-8` is passed — a latent ABI mismatch in that configuration.

## Changes Made

- `src/sparse/pardisotools.f90`: `include 'mkl_pardiso.f90'` and `use mkl_pardiso` (vendor interface, as
  `mumpstools` includes `dmumps_struc.h`); `gMemoryPointers` is `type( MKL_PARDISO_HANDLE )`; one size-1 placeholder
  array per dummy MKL may write (`tPerm`, `tB`, `tX`, plus `tA`, `tIA`, `tJA` in the release phase) instead of one
  shared scalar; `aRHS` is `intent( inout )` because MKL's `B` is; `gDPARM` and `pardisotools_get_determinant`
  removed; `gParameters( 33 )` always 0; `#error` under `BELFEM_INT64` naming the missing `-fdefault-integer-8`.
- `src/sparse/pardisotools.hpp`: `aRHS` non-const, determinant prototype removed.
- `src/sparse/cl_SolverPARDISO.{hpp,cpp}`: `get_determinant()` raises `BELFEM_ERROR( "MKL PARDISO does not compute
  a determinant …" )`; parameter-8 and info-slot-8 comments say "reserved".
- `tests/sparse/test_Solver.cpp`: `SolverSolve.PARDISODeterminantIsUnsupported`.
- Also this morning: `fn_symrationspace.hpp` half-count in integer arithmetic (the last `-Wfloat-conversion`
  warning of the test build); the material citation typos and the γ(T\*) table (73e7734).

## Jury

Both auditors rejected my drafted interface block (mixed kinds, wrong intents for `perm`/`b`/`x`, shared dummies
violating the INTENT( OUT ) aliasing rule, and an `aRHS` intent that would not have compiled) and agreed on the
redesign above, on dropping `gDPARM`, on not claiming ILP64, and on error-plus-deletion for the determinant. Grok's
reminder to keep `-fallow-argument-mismatch` (MUMPS and ARPACK still need it) is followed. Exchange:
`tmp/ai_exchange/review_pardiso_interface.md`.

## Open Questions

- ILP64 + PARDISO is now a compile-time refusal rather than a silent mismatch; making it work means passing
  `-fdefault-integer-8` to the Fortran sources (Intel's documented route) or calling `pardiso_64` with every
  integer 64-bit, including `iparm` — a decision for when someone needs that configuration.

## Files Updated

- src/sparse/pardisotools.f90, src/sparse/pardisotools.hpp
- src/sparse/cl_SolverPARDISO.cpp, src/sparse/cl_SolverPARDISO.hpp
- tests/sparse/test_Solver.cpp
- src/math/tools/fn_symrationspace.hpp
- CHANGELOG.md, devlog/README.md
