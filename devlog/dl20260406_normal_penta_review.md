# Devlog 2026-04-06 — normal_penta Review

**Date:** 2026-04-06
**Topic:** Read-only review of the new `normal_penta()` / `normal_penta_ts()` implementation
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Reviewed the new prism normal functions in `cl_FEM_Calculator`. The overall direction is correct, but the current patch is not ready: it introduces immediate compile errors, the new functions are not yet connected to the normal-dispatch switch, and the TS/general top-face formulas still contain sign/normalization defects.

## Key Findings

- The current branch does not build because `mFunInvertJ` was changed to `mFunInverJ` and `mDetJ` / `mDetJIndex` were changed to `mDeJ` / `mDeJIndex` in `cl_FEM_Calculator.cpp`, while the declarations in `cl_FEM_Calculator.hpp` were not changed.
- `normal_penta()` and `normal_penta_ts()` are defined, but `Calculator::allocate()` still lacks a `GeometryType::PENTA` dispatch branch, so the original `No normal function assigned` abort path remains.
- `normal_penta_ts()` case `0` and case `1` do not produce opposite normals for the two TS faces, even though `PENTA6TS` facets `0` and `1` are opposite `zeta=-1` / `zeta=+1` faces.
- `normal_penta()` case `0` has a `break` before normalization, leaving `mSurfaceIncrement` unset for that face.
- The z-component sign in `normal_penta()` case `4` and `normal_penta_ts()` case `1` is opposite to the standard `cross( dX/dxi, dX/deta )` expression used elsewhere in the file.

## Changes Made / Proposed

- No source changes made.
- Verified the build failure with `cmake --build cmake-build-debug --target hphirun -j4`.
- Appended the audit to `todo/ai_exchange.md`.

## Open Questions

- Whether the intended dispatch should distinguish `PENTA6TS` / `PENTA18TS` from regular prism volume elements inside the `GeometryType::PENTA` branch.
- Whether the TS facet formulas should be derived directly from the TS canonical facet mapping to avoid another orientation mismatch.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260406_normal_penta_review.md
