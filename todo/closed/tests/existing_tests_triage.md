# Existing Test Suite Triage

**Date:** 2026-03-22
**Purpose:** Assess existing tests in `./test/`, determine what to keep, salvage, or retire before implementing the new test plans (tests_00 through tests_13).

---

## Key Finding: No Framework Migration Needed

All existing tests already use **Google Test**. Every directory has a `test_*_main.cpp` that initializes MPI (`gComm`), the Logger (`gLog`), calls `testing::InitGoogleTest()`, and `RUN_ALL_TESTS()`. The build system links with `-lgtest -lgtest_main` via `config/scripts/Add_Test.cmake`.

---

## Current Test Inventory

| Directory | Files | What It Tests | Quality | Build Status |
|---|---|---|---|---|
| `test/core/` | 1 test | `clean_string()`, `string_to_words()` | Smoke | Builds |
| `test/linalg/` | 2 tests | Vector/Matrix basics (init, min/max, operators) | Moderate | Builds |
| `test/sparse/` | 3 tests | SpMatrix CSR/CSC, indexing, SpMV, UMFPACK/PARDISO/MUMPS/PETSc solve | Good | Builds |
| `test/spline/` | 1 test | Spline interpolation on NASA thermo data (cp, h, s, derivatives) | Thorough | Builds |
| `test/math/tools/` | 2 tests | Cardano (2 cases), cubic Bezier (1 point) | Minimal | Builds |
| `test/math/tensor/` | 5 tests | contract42, contract44, rotate, kelvin_christoffel, invert_symmetric | Very good | Builds |
| `test/fem/` | 9 tests | Normal computation (TRI3/6, QUAD4/9, TET4/10, HEX8/27), interface integration | Data-driven | **Partial** — `cl_IntegrationData_Interface.cpp` cannot compile (missing `FEM_geometry.hpp`) |
| `test/maxwell/` | 4 tests | End-to-end Maxwell h-phi solver (TRI3/6, TET4/10) | Integration | **Disabled** — commented out in `test/CMakeLists.txt` |
| `test/physics/gastables/` | 4 tests | Gas data parsing, N2 reference properties vs NIST | Thorough | Builds |
| `test/physics/gasmodels/` | 15 tests | EoS, caloric, vapor, transport for H2/O2/CH4/air vs NIST/VDI | Excellent | Builds |

---

## Recommended Approach: Keep In Place, Extend Alongside

**Do NOT move existing tests to `test/old/`.** Reasons:

1. The existing tests already use GTest and follow BELFEM conventions — they're not broken legacy code.
2. The gasmodels and gastables tests (19 files) are outside the scope of tests_01–13 and should continue running as-is.
3. The maxwell tests are disabled but represent future integration test targets.
4. Moving files breaks the existing CMake build and git history for no gain.

**Instead:** Add new test files alongside the existing ones, in the same directory structure. Where an existing test partially covers a new plan's scope, extend it rather than replacing it.

---

## Salvage Plan by Module

### Tier 1: Direct Salvage (copy test logic into new test files)

| Existing File | New Plan | Action |
|---|---|---|
| `test/math/tensor/fn_ddot_contract42.cpp` | tests_04 | **Keep as-is.** Already uses the naive-loop-vs-optimized pattern that tests_04 recommends. Add the missing tests (construction, fill, copy/move, mat-to-ten conversion) in a new file alongside. |
| `test/math/tensor/fn_ddot_contract44.cpp` | tests_04 | **Keep as-is.** Same reasoning. |
| `test/math/tensor/fn_rotate.cpp` | tests_04 | **Keep as-is.** Uses C/C-SiC elasticity data — excellent physics-based validation. |
| `test/math/tensor/fn_invert_symmetric.cpp` | tests_04 | **Keep as-is.** C * C^{-1} = I identity check. |
| `test/math/tensor/fn_kelvin_christoffel.cpp` | tests_04 | **Keep as-is.** |
| `test/sparse/cl_SpMatrix_CSR.cpp` | tests_08 | **Keep as-is.** Graph-based construction, indexing, SpMV, solver backends. Add dense-constructor and structural-zero tests in new file. |
| `test/sparse/cl_SpMatrix_CSC.cpp` | tests_08 | **Keep as-is.** |
| `test/sparse/cl_SolverPETSC.cpp` | tests_08 | **Keep as-is** (conditional on PETSc availability). |
| `test/spline/cl_Spline.cpp` | tests_09 | **Keep as-is.** NASA thermo data validation is excellent. Add construction/boundary-condition/edge-case tests in new file. Do NOT add save/load roundtrip until BUG-S1 is fixed. |
| `test/linalg/cl_Vector.cpp` | tests_02 | **Keep as-is.** Add comprehensive tests (copy/move, data layout, BLAS, submatrix) in new file. |
| `test/linalg/cl_Matrix.cpp` | tests_02 | **Keep as-is.** Same reasoning. |

### Tier 2: Extend with Caution

| Existing File | New Plan | Action |
|---|---|---|
| `test/math/tools/fn_cardano.cpp` | tests_05 | **Keep.** Add root-substitution validation, degenerate cases, and the other math functions (rotation_matrix, find_interval, etc.) in new file(s). |
| `test/math/tools/fn_cubic_bezier.cpp` | tests_05 | **Keep.** Only covers 1 test point — extend with endpoint/derivative tests in new file. |
| `test/fem/fn_normal_*.cpp` (8 files) | tests_11, tests_13 | **Keep.** These test normal computation and surface integrals — useful geometry validation not directly covered by the new plans. |
| `test/fem/cl_IntegrationData_Interface.cpp` | tests_13 | **Keep but note it cannot compile.** The test pattern (master/slave physical-space coincidence for all orientations) is exactly what tests_13 targets. Use it as reference for the Layer 1 approach described in the new plan. When the `create_test_mesh` infrastructure is built, this file can be enabled. |

### Tier 3: Leave Alone (out of scope)

| Existing Files | Reason |
|---|---|
| `test/physics/gasmodels/` (15 files) | Not covered by tests_01–13. Excellent quality, keep running. |
| `test/physics/gastables/` (4 files) | Same. |
| `test/maxwell/` (4 files) | Disabled, integration-level. Future work. |
| `test/core/stringtools.cpp` | No core test plan yet. Keep as-is. |

---

## New Test File Placement

For each new test plan, add files in the existing directory structure:

```
test/
├── containers/              # NEW — tests_01
│   ├── CMakeLists.txt
│   ├── test_containers_main.cpp
│   ├── test_Cell.cpp
│   ├── test_DynamicBitset.cpp
│   ├── test_ShiftRegister.cpp
│   ├── test_Bitset.cpp
│   ├── test_Map.cpp
│   ├── test_Set.cpp
│   ├── test_Queue.cpp
│   └── test_StringList.cpp
│
├── linalg/
│   ├── cl_Vector.cpp         # EXISTING — keep
│   ├── cl_Matrix.cpp         # EXISTING — keep
│   ├── test_LinalgExtended.cpp  # NEW — additional tests_02 coverage
│   └── ...
│
├── math/
│   ├── tools/
│   │   ├── fn_cardano.cpp       # EXISTING — keep
│   │   ├── fn_cubic_bezier.cpp  # EXISTING — keep
│   │   ├── test_MathTools.cpp   # NEW — tests_05 remaining functions
│   │   └── ...
│   ├── tensor/
│   │   ├── fn_ddot_contract42.cpp  # EXISTING — keep
│   │   ├── fn_rotate.cpp          # EXISTING — keep
│   │   ├── ...                    # EXISTING — keep all
│   │   ├── test_Tensor.cpp        # NEW — tests_04 class/construction tests
│   │   └── ...
│   └── quaternion/            # NEW — tests_03
│       ├── CMakeLists.txt
│       ├── test_quaternion_main.cpp
│       └── test_Quaternion.cpp
│
├── sparse/
│   ├── cl_SpMatrix_CSR.cpp    # EXISTING — keep
│   ├── cl_SpMatrix_CSC.cpp    # EXISTING — keep
│   ├── cl_SolverPETSC.cpp     # EXISTING — keep
│   ├── test_SpMatrixExtended.cpp  # NEW — tests_08 additions
│   └── ...
│
├── spline/
│   ├── cl_Spline.cpp          # EXISTING — keep
│   ├── test_SplineExtended.cpp  # NEW — tests_09 additions
│   └── ...
│
├── comm/                      # NEW — tests_06
│   ├── CMakeLists.txt
│   └── ...
│
├── mesh/                      # NEW — tests_10
│   ├── CMakeLists.txt
│   └── ...
│
├── graph/                     # NEW — tests_07
│   ├── CMakeLists.txt
│   └── ...
│
├── fem/
│   ├── fn_normal_*.cpp          # EXISTING — keep all 8
│   ├── cl_IntegrationData_Interface.cpp  # EXISTING — keep (non-compilable, future)
│   ├── test_InterpolationFactory.cpp     # NEW — tests_11
│   ├── test_InterpolationLagrange.cpp    # NEW — tests_11
│   ├── test_IntegrationDispatch.cpp      # NEW — tests_12
│   ├── test_IntegrationInvariants.cpp    # NEW — tests_12
│   ├── test_FacetOrientationTet.cpp      # NEW — tests_13
│   ├── test_FacetOrientationHex.cpp      # NEW — tests_13
│   └── ...
│
├── maxwell/                   # EXISTING — keep disabled
├── physics/                   # EXISTING — keep all
└── CMakeLists.txt             # UPDATE — add new subdirectories
```

---

## Test Main Pattern (for new directories)

All existing test mains follow this pattern — new ones should too:

```cpp
#include <gtest/gtest.h>
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"

belfem::Communicator gComm;
belfem::Logger       gLog( 5 );

int main( int argc, char * argv[] )
{
    gComm = belfem::Communicator( argc, argv );
    testing::InitGoogleTest( &argc, argv );
    int tResult = RUN_ALL_TESTS();
    gComm.finalize();
    return tResult;
}
```

For modules that don't need MPI (containers, linalg, math), the Communicator init can be omitted but keeping it is harmless and consistent.

---

## Patterns Worth Emulating from Existing Tests

1. **Naive-loop reference implementations** (`test/math/tensor/`) — the gold standard for testing optimized kernels.
2. **NIST/literature-validated data** (`test/physics/gasmodels/`) — r2 goodness-of-fit against published tables.
3. **Graph-based SpMatrix construction** (`test/sparse/`) — realistic sparse pattern, not just dense-to-sparse conversion.
4. **NASA polynomial data for splines** (`test/spline/`) — uses real thermodynamic coefficients, analytical derivatives.

---

## What NOT to Do

- **Do not move/rename existing test files.** It breaks CMake, git blame, and existing CI.
- **Do not merge existing and new tests into single files.** Keep them separate for clear provenance.
- **Do not add `#define protected public` hacks** (seen in `test/physics/`) — test through public API only, per tests_00_strategy.md.
- **Do not duplicate existing coverage.** If `cl_Matrix.cpp` already tests `operator+`, don't re-test it in the new file. Add only what's missing.

---

## Summary: Action Items

| Priority | Action | Effort |
|---|---|---|
| 1 | Create `test/containers/` directory with CMakeLists.txt and test files (tests_01) | New work |
| 2 | Create `test/math/quaternion/` for tests_03 | New work |
| 3 | Create `test/comm/`, `test/graph/`, `test/mesh/` for tests_06, 07, 10 | New work |
| 4 | Add extended test files alongside existing in `test/linalg/`, `test/sparse/`, `test/spline/`, `test/math/tools/`, `test/math/tensor/` | Extend existing |
| 5 | Add interpolation/integration/orientation tests in `test/fem/` | New work alongside existing |
| 6 | Update top-level `test/CMakeLists.txt` to include new subdirectories | Build config |
| 7 | Fix `cl_IntegrationData_Interface.cpp` compile issues when mesh test infrastructure is ready | Future |
