# Devlog 2026-07-01 — User Material Backend Boundary

**Date:** 2026-07-01
**Topic:** Read-only investigation of user-defined material compilation and linalg backend coupling
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** medium-high
**Literature References:** N/A

## Summary

Investigated why the standalone user-defined material example under
`cmake-build-debug/MatData` cannot be backend-agnostic as currently structured.
No source code was modified.

## Key Findings

- `cl_Material.hpp` includes `cl_Spline.hpp` directly (`src/physics/materials/cl_Material.hpp:58-59`), stores `Cell<Spline*>` (`:281-283`), exposes `set_spline()` (`:1040-1043`), and has inline methods that dereference `Spline` (`:1689-1707`). A user material that includes `cl_Material.hpp` therefore imports spline at compile time.
- `cl_Material.hpp` also exposes `set_user_defined_polynomial(..., const Vector<real>&)` (`:984-987`), so even the user-material convenience API names the backend-selected `Vector` type.
- `cl_Spline.hpp` includes `cl_Vector.hpp` and `cl_SpMatrix.hpp` (`src/numerics/spline/cl_Spline.hpp:15-20`) and exposes `Vector<real>`, `SpMatrix`, and `Matrix<real>` in its public interface/storage (`:66`, `:72-90`, `:135-141`, `:254-270`, `:301-341`).
- `cl_SpMatrix.hpp` includes `cl_Vector.hpp` (`src/sparse/cl_SpMatrix.hpp:17-22`) and exposes `Vector<real>` in multiplication APIs/operators (`:384-398`, `:598-604`).
- `cl_Vector.hpp` selects a concrete implementation only under `BELFEM_ARMADILLO` or `BELFEM_BLAZE` (`src/linalg/cl_Vector.hpp:15-24`). With neither macro defined, `Vector` and then `Matrix` are undefined.
- The MatData material target currently receives only `-Dmat_EXPORTS` and does not import BELFEM backend definitions or third-party include paths (`cmake-build-debug/MatData/build/CMakeFiles/mat.dir/flags.make`). The template has `BELFEM_CACHE_LOCATIONS` but does not consume it (`cmake-build-debug/MatData/CMakeLists.txt:68-73`).
- A syntax-only compile with complete local BELFEM include paths but no backend macro fails at `cl_Vector.hpp`/operator headers with "`Vector` does not name a type". Adding `-DBELFEM_ARMADILLO` progresses to the next missing dependency (`<armadillo>` include path), confirming the issue is usage-requirement propagation, not material logic.

## Changes Made / Proposed

- No source changes made.
- Recommended immediate direction: make the standalone material build import the exact compile definitions, include directories, and relevant third-party include paths from the BELFEM build it will be loaded into.
- Recommended longer-term direction: split the user-material-facing ABI from spline/linalg-heavy framework internals. The plugin-facing header should expose scalar property registration and backend-neutral coefficient inputs; spline construction and `Vector`/`SpMatrix` ownership should remain inside BELFEM or behind non-inline implementation boundaries.

## Open Questions

- Whether user material libraries should ever be allowed to construct BELFEM splines directly, or whether they should only register scalar functions/data and let BELFEM construct interpolation objects.
- Whether BELFEM should provide an exported CMake package/usage target for out-of-tree plugins, or a smaller generated config file dedicated to user materials.

## Files Updated

- `devlog/dl20260701_user_material_backend_boundary.md`
