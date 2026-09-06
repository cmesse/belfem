# Devlog 2026-07-01 - NLopt CMake Audit

**Date:** 2026-07-01
**Topic:** Audit of initial nlopt CMake configuration
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Reviewed the new nlopt CMake integration and surrounding link-list changes. No source changes were made by Codex. The `USE_NLOPT=OFF` configure path succeeds with HDF5/Exodus/examples disabled for the audit command, but `USE_NLOPT=ON` currently fails during configuration.

## Key Findings

- `CMakeLists.txt` includes `config/numerics/find_nlopt.cmake` and then calls `find_package(NLOPT REQUIRED)`, but the included script is not a CMake package module or config file. Configure with `-DUSE_NLOPT=ON` fails looking for `FindNLOPT.cmake`, `NLOPTConfig.cmake`, or `nlopt-config.cmake`.
- `config/numerics/find_nlopt.cmake` never sets `NLOPT_LIB_DIR`, yet appends `-L${NLOPT_LIB_DIR}` to `BELFEM_OTHER_LIBS`.
- The custom search checks for `nlopt.h` but does not append the discovered include directory to `BELFEM_INCLUDES`, so non-system installs would still fail when wrapper source includes `<nlopt.h>`.
- `USE_NLOPT` remains commented out in the top-level options block, so it is not exposed as a normal project option.
- The broad `BELFEM_IO_LIBS` to `BELFEM_OTHER_LIBS` rename appears internally consistent in the touched build scripts, but it is a cross-cutting change and should stay intentionally documented.

## Changes Made / Proposed

- Proposed either removing the extra `find_package(NLOPT REQUIRED)` call and making `find_nlopt.cmake` fully self-contained, or converting the custom file into a real `FindNLOPT.cmake` on `CMAKE_MODULE_PATH`.
- Proposed finding both header and library paths, appending the include directory to `BELFEM_INCLUDES`, appending a concrete library path or `-lnlopt`, and defining `USE_NLOPT` as a proper option.

## Open Questions

- Whether BELFEM should prefer an SCLS-only lookup first and treat `/usr`/`/usr/local` as fallback, or allow any package-config installation through `find_package`.

## Files Updated

- `devlog/dl20260701_nlopt_cmake_audit.md`
- `devlog/README.md`
