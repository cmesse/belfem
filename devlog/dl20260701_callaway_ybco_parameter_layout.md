# Devlog 2026-07-01 — Callaway YBCO Parameter Layout

**Date:** 2026-07-01
**Topic:** Callaway thermal-conductivity parameter packing and YBCO regression coverage
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** Existing in-tree references only

## Summary

Fixed the Callaway thermal-conductivity unpacking bug where the Fortran kernel
read the Umklapp temperature-correction parameter `d` from the optical Raman
shift slot. The corrected layout is:

- `b`: C++ `params[11]` / Fortran `params(12)`
- `d`: C++ `params[12]` / Fortran `params(13)`
- `lambda_opt`: C++ `params[18]` / Fortran `params(19)`
- `omega_opt`: C++ `params[19]` / Fortran `params(20)`

## Changes Made

- `src/physics/materials/debye.f90`: changed `d = params(20)` to
  `d = params(13)` and removed the stale FIXME. Added local C++/Fortran index
  comments for the Umklapp and optical-phonon slots.
- `src/physics/materials/debye.hpp`: updated the 20-parameter layout comment so
  every slot lists both the C++ zero-based index and the Fortran one-based index.
- `src/physics/materials/cl_Material_YBCO.cpp`: clarified comments in both
  YBCO Callaway packers; no physical parameter values were changed.
- `src/physics/materials/doc/callaway_thermal_conductivity.md` and
  `src/physics/materials/doc/README.md`: removed stale aliasing warnings and
  documented the corrected layout.
- Added active `tests/physics` with:
  - a direct `callaway_conductivity()` regression that verifies changing
    `params[12]` changes the Umklapp result while changing `params[19]` does not
    affect Umklapp when optical scattering is off;
  - public `YBCO::lambda(T)` regression values at 20 K, 77 K, 100 K, and 300 K;
  - finite/positive/smooth checks over 1-400 K and around `T_crit = 92.5 K`.

## Verification

- `cmake -S . -B cmake-build-debug -DUSE_TEST=ON`
- `cmake --build cmake-build-debug --target material -j2`
- `cmake --build cmake-build-debug --target test_physics -j2`
- `./cmake-build-debug/test/test_physics --gtest_filter=CallawayConductivityParameterLayout.*:YBCOThermalConductivity.*`
- `ctest --test-dir cmake-build-debug -R physics --output-on-failure`
- `./cmake-build-debug/bin/material YBCO -t 20 300 57`
- `./cmake-build-debug/bin/material YBCO -t 88 98 1`
- `./cmake-build-debug/bin/material YBCO -t 1 400 39.9`

Focused builds/tests passed. The default full build was attempted with
`cmake --build cmake-build-debug -j2` and failed in an unrelated existing target:
`src/fem/kernel/staticheat.cpp:92` references `MaterialType::Copper`, which is
not a member of the current `MaterialType` enum.

The corrected YBCO public thermal-conductivity values at the regression points
are unchanged from the old rounded CLI output and the direct old-aliased
comparison was identical to 17 printed digits for 20 K, 77 K, 100 K, and 300 K.
No refitting appears necessary from this verification.

## Files Updated

- `src/physics/materials/debye.f90`
- `src/physics/materials/debye.hpp`
- `src/physics/materials/cl_Material_YBCO.cpp`
- `src/physics/materials/doc/callaway_thermal_conductivity.md`
- `src/physics/materials/doc/README.md`
- `tests/CMakeLists.txt`
- `tests/physics/CMakeLists.txt`
- `tests/physics/test_physics_main.cpp`
- `tests/physics/test_YBCO.cpp`
- `devlog/dl20260701_callaway_ybco_parameter_layout.md`
