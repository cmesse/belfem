# Devlog 2026-07-01 - NLopt Wrapper Design

**Date:** 2026-07-01
**Topic:** Read-only design pass for a general `nlopt` wrapper
**AIs involved:** Codex
**Codex Audit Confidence:** medium-high
**Literature References:** N/A

## Summary

Reviewed BELFEM's numerics, material-property, and CMake patterns to guide a first implementation of a general `nlopt` wrapper for future material-property optimization work. No source files were modified.

## Key Findings

- The top-level CMake already has a commented `USE_NLOPT` option, so the intended build integration point is visible in `CMakeLists.txt`.
- A wrapper belongs under `src/numerics/`, likely `src/numerics/optimization/`, not directly in `src/physics/materials/`. The material model should consume the optimizer rather than own the third-party interface.
- BELFEM currently models stored elastic data as `E` and `nu`; `G(T)` and `K(T)` are algebraically derived in `cl_Material.hpp`. A later K/G-first generator can still emit final `E`/`nu` splines without immediately expanding `MaterialProperty`.
- Existing numerics code favors explicit interfaces and preallocated BELFEM containers (`Vector<real>`, `Cell<T>`) over `std::function`-heavy APIs.
- External-library failure handling should be normalized into BELFEM-style results and `BELFEM_ERROR` for fatal configuration/runtime failures.

## Changes Made / Proposed

- Proposed a small optimizer module:
  - `src/numerics/optimization/cl_OptimizerProblem.hpp`
  - `src/numerics/optimization/cl_NloptWrapper.hpp`
  - `src/numerics/optimization/en_OptimizerAlgorithm.hpp`
  - `src/numerics/optimization/en_OptimizerStatus.hpp`
- Proposed using the nlopt C API internally to avoid exposing `std::vector` copies or exception-style semantics at BELFEM call sites.
- Proposed first tests as Rosenbrock/quadratic bounded problems, plus a compile-time unlinked path when `USE_NLOPT=OFF`.

## Open Questions

- Whether the material generator should remain an offline executable/tool or become a reusable in-library material fitting module.
- Whether final material storage should keep the existing `E`/`nu` spline contract or add first-class `K`/`G` material properties.
- Whether nonlinear constraints should be included in the first wrapper pass or deferred until the basic bounded objective path is stable.

## Files Updated

- `devlog/dl20260701_nlopt_wrapper_design.md`
- `devlog/README.md`
