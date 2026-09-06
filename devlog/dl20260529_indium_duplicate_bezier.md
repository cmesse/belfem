# Devlog 2026-05-29 - Indium Duplicate Bezier Symbols

**Date:** 2026-05-29
**Topic:** Diagnosis of duplicate `belfem::Bezier` symbols when linking `bin/indium`
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Investigated the `bin/indium` link failure reporting duplicate `belfem::Bezier`
symbols. The failure is caused by `src/physics/materials/cl_Material_Indium.hpp`
including `cl_Bezier.cpp` instead of `cl_Bezier.hpp`.

## Key Findings

- `src/physics/materials/cl_Material_Indium.hpp:10` includes the implementation
  file `cl_Bezier.cpp`, so every translation unit that includes the indium
  header emits out-of-line `Bezier` function definitions.
- `src/physics/materials/indium.cpp:12` includes `cl_Material_Indium.hpp`, so
  the executable object contains `Bezier` definitions.
- `src/physics/materials/cl_Material_Indium.cpp:5` also includes
  `cl_Material_Indium.hpp`, and this object is archived into
  `libbelfem_materials.a`, so the same `Bezier` definitions appear there too.
- `src/numerics/bezier/CMakeLists.txt:3-4` already compiles `cl_Bezier.cpp` as
  the `bezier` library; the implementation file should not be included by a
  public material header.

## Changes Made / Proposed

- Proposed source fix: change `#include "cl_Bezier.cpp"` to
  `#include "cl_Bezier.hpp"` in `src/physics/materials/cl_Material_Indium.hpp`.
- No source-code edits were made in this read-only investigation.

## Open Questions

- After applying the include fix, rebuild `indium` from a clean or refreshed
  build tree to ensure stale object files are not reused.

## Files Updated

- devlog/dl20260529_indium_duplicate_bezier.md
- devlog/README.md
