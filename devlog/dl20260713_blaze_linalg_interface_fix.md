# Devlog 2026-07-13 — Blaze Linalg Interface Fix

**Date:** 2026-07-13
**Topic:** Fix Blaze matrix/vector wrapper build failures
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Fixed the Blaze backend build failures reported from `cl_TimestepMatrices.cpp`,
`mt_maxwell_h.cpp`, and `cl_Material_UserDefined.cpp`. The root cause was a
pair of wrapper-level issues: overly broad matrix expression conversion in
`Matrix<T>`, ambiguous compound assignment overloads for Blaze matrix
expressions, and an Armadillo-style `set_size(size, 1)` call left in the Blaze
`Vector<T>` wrapper.

The requested `hphirun` target now builds successfully with the Blaze
configuration in `cmake-build-debug`.

## Key Findings

- `Matrix<T>::Matrix(const blaze::Expression<ET>&)` accepted Blaze vector
  expressions and then failed inside `DynamicMatrix` construction. The wrapper
  now accepts Blaze matrix expressions explicitly, and maps dense vector
  expressions to `n x 1` or `1 x n` matrices to preserve existing Armadillo-style
  call-site behavior such as `Matrix<real> Ctj(trans(C) * j)`.
- `Matrix<T>::operator-=(...)` had competing implicit materialization paths for
  Blaze matrix expressions, causing ambiguous overloads at the reported IWG and
  Maxwell sites. Dense Blaze matrix expressions now bind directly to compound
  operators.
- `Vector<T>(const std::vector<T>&)` called `mVector.set_size(aVector.size(), 1)`,
  which is not a Blaze `DynamicVector` API. It now uses `resize(size, false)`.
- The final executable link exposed a separate header ODR issue: two
  `cl_FEM_Calculator.hpp` member definitions were missing `inline`, unlike the
  surrounding header-defined functions.

## Changes Made

- `src/linalg/blaze/cl_BZ_Matrix.hpp`: constrained matrix expression
  construction/assignment to Blaze matrix CRTP types, added dense vector
  construction/assignment as column/row matrices, and added direct Blaze matrix
  compound operators.
- `src/linalg/blaze/cl_BZ_Vector.hpp`: fixed the `std::vector<T>` constructor
  resize call for Blaze.
- `src/fem/kernel/cl_FEM_Calculator.hpp`: added `inline` to
  `compute_lambda_metal()` and `compute_dlambdadT_metal()` definitions to fix
  multiple-definition link errors.

## Verification

- `make -C cmake-build-debug hphirun -j2` completed successfully.

## Open Questions

- No open Blaze wrapper compile blockers remain from this session.
- The worktree contained unrelated in-progress edits before and during the
  session; those were preserved.

## Files Updated

- `src/linalg/blaze/cl_BZ_Matrix.hpp`
- `src/linalg/blaze/cl_BZ_Vector.hpp`
- `src/fem/kernel/cl_FEM_Calculator.hpp`
- `devlog/dl20260713_blaze_linalg_interface_fix.md`
- `devlog/README.md`
