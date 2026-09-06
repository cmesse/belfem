# `tests/old` Triage: Revive `core`, Delete the Rest

**Date:** 2026-08-05
**Purpose:** Decide, per module, which of the stale `tests/old` suites to revive and which to remove
**Module:** `tests/`

## Starting point

`tests/old/` held 47 files across `core`, `linalg`, `sparse`, `spline`,
`math/{tools,tensor}`, `fem`, `maxwell` and `physics`. None of it was being
built: `tests/CMakeLists.txt` never had an `add_subdirectory( old )`, and
`tests/old/maxwell` was additionally commented out inside
`tests/old/CMakeLists.txt`. The `physics/` subtree (gasmodels, gastables) was
out of scope here — a concurrent session owns it under
`todo/gasmodels_open_source_migration.md`.

## Verdicts

| Old suite | Verdict | Basis |
|---|---|---|
| `core/stringtools.cpp` | **revived** | `clean_string` / `string_to_words` alive and load-bearing for the input parser; no coverage anywhere in the new suite |
| `linalg/cl_{Vector,Matrix}.cpp` | deleted | superseded by `tests/linalg/test_Vector.cpp` (76 cases) + `test_Matrix.cpp` (70) |
| `sparse/cl_SpMatrix_{CSR,CSC}.cpp`, `cl_SolverPETSC.cpp` | deleted | superseded by `tests/sparse/test_SpMatrix.cpp` and `test_Solver.cpp::SolverSolve.PETScSolveTridiagonal` |
| `spline/cl_Spline.cpp` | deleted | superseded by `tests/math/test_Spline.cpp` (33 cases) |
| `math/tools/{fn_cardano,fn_cubic_bezier}.cpp` | deleted | superseded by `tests/math/test_MathTools.cpp` + `test_Polynomials.cpp` (which also covers `ferrari`) |
| `math/tensor/*.cpp` (5) | deleted | superseded by `tests/math/test_Tensor.cpp` + `test_TensorKernels.cpp` (35 cases, same four kernels) |
| `fem/fn_normal_*.cpp` (8), `cl_IntegrationData_Interface.cpp` | deleted, gap recorded | dead API — see below |
| `maxwell/test_{tri3,tri6,tet4,tet10}.cpp` | deleted | dead API — see below |

## The two dead-API suites

Both were unrevivable rather than merely stale, so neither was ported.

**`fem/`** — all nine tests are written against `FEM_geometry.hpp` and a
`belfem::fem::geometry` namespace (`normal_tri`, `normal_quad`, `normal_tet`,
`normal_hex`, `test_pipette`, `test_intpoints`, `create_test_mesh`). That header
no longer exists anywhere in `src/`; the normal computation now lives as member
functions on the calculator (`src/fem/kernel/cl_FEM_Calculator.hpp:1249-1284`).
An earlier triage (`todo/closed/tests/existing_tests_triage.md`) had marked these
"keep" and noticed only that `cl_IntegrationData_Interface.cpp` wouldn't compile;
in fact none of the nine can.

`cl_IntegrationData_Interface.cpp` is genuinely superseded — its concern (mapping
facet integration points through every orientation permutation) is what
`tests/fem/test_FacetIntegrationPoints.cpp` now does. The eight normal tests are
**not** superseded: `Calculator::normal_tri_straight/_curved`,
`normal_quad_straight/_curved`, `normal_hex` and `cl_Pipette` have zero coverage
today. Their reference fixture is API-independent and still valid, so it was
kept and moved to `tests/fem/support/test_database.hdf5`; the rebuild is tracked
in `todo/test_normals_and_pipette_coverage.md`.

**`maxwell/`** — these drove `fem::MaxwellJob::initialize_test` / `run_test`
against per-element `IwgType::MAXWELL_HPHI_{TRI3,TRI6,TET4,TET10}` enum entries.
`MaxwellJob` does not exist in `src/` at all, and the enum has since collapsed to
a single `IwgType::Maxwell` (`src/fem/iwg/en_IWGs.hpp:73`). An end-to-end h-φ
regression test remains worth having, but it is a new design rather than a
restoration; `test_maxwell.hdf5` (281 KB of results keyed to the old job API) was
dropped with it.

## What landed

New `tests/core/` (`test_core_main.cpp`, `test_StringTools.cpp`, `CMakeLists.txt`),
registered in `tests/CMakeLists.txt` and labelled `fast`. `core` is already in
`BELFEM_LIBLIST_BASE`, so no extra `LIBLIST` or include directories were needed.

The old file was one test over two functions. The revived version is 45 cases
over the whole header — `clean_string`, `first_word`, `string_to_words`,
`search_and_replace`, `string_to_{lower,upper,bool}`, `to_real`, `is_integer`,
the four path helpers, `format_with_leading_zeros`, `utf8_character_count`,
`unit_to_si`, `string_to_cell`, `to_pair`, `datatype_string`.

Several expectations pin down current quirks rather than an idealized contract,
and are flagged `QUIRK` in the source so nobody "fixes" the test instead of the
code:

- `clean_string` truncates at `#` but leaves the space in front of it, because
  the trailing-space erase only fires when the loop runs off the end of the
  string.
- `to_real("")` returns `0.0`, not NAN — the guard requires a non-terminating
  character after the failed parse. `to_real("3.5Alpha")` returns `3.5`.
- `filetype()` on an extensionless path returns the whole path: `find_last_of`
  gives `npos`, and `npos + 1` wraps to 0.
- `format_with_leading_zeros()` returns a printf *format string* (`"%02u"`), not
  a formatted number.
- `string_to_bool(" true")` is false — the input is not trimmed.

Two of the hand-derived expectations were wrong and were caught by running the
binary, not by inspection. The second is worth knowing about:

**Unqualified `basename()` does not resolve to `belfem::basename`.** glibc's
`<string.h>` declares `char * basename( const char * )`, and for a string-literal
argument that C overload is the exact match, so `EXPECT_EQ` compared two
pointers and passed nothing meaningful. The test now qualifies `belfem::` on all
four path helpers. This is not a defect in `stringtools` — callers passing a
`std::string` get the BELFEM overload either way — but it is a live footgun for
any translation unit that does `using namespace belfem` and hands a `const
char *` to `basename`.

## Verification

Compiled and linked standalone against the prebuilt `libbelfem_{core,containers,comm}.a`
with the binary written to scratch, so the shared `cmake-build-debug/` tree was
never touched: **45/45 pass**. The suite has not yet been run through the real
CMake target — that build is Christian's to drive.

## Left open

- `tests/old/CMakeLists.txt` survives as the only file under `tests/old/`, still
  listing the seven now-deleted subdirectories. It is inert (nothing includes
  it). Deleting it and the directory waits on the gasmodels session finishing
  with `tests/old/physics`.
- `todo/test_normals_and_pipette_coverage.md` — R1-R5 for the normals/pipette
  rebuild. R1 (confirming the outward-normal sign convention against the stored
  reference columns) is a physics call, not a test-style call.
