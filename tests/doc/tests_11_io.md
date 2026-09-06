# BELFEM I/O Module Tests — Detailed Plan (HDF5 Focus)

**Date:** 2026-03-24 (updated 2026-03-24)
**Purpose:** Method-level test matrix for the I/O module (`src/io/`), focused on HDF5 read/write
**Depends on:** `tests_0_strategy.md`, `tests_1_containers.md` (Cell, Vector, Matrix)
**Confidence:** High on HDF5 class round-trips. Medium on group tree navigation. Low on parallel I/O.
**Status:** 42 tests implemented and passing.

---

## Module Overview

| Component | In Scope | Deferred |
|---|---|---|
| `HDF5` class | Scalar, string, vector, matrix, bool save/load round-trips | — |
| `HDF5_Tools` | Low-level template functions (tested indirectly via `HDF5` class) | Direct template testing |
| Group management | `create_group`, `select_group`, `close_active_group`, nested trees | — |
| File modes | `NEW`, `OPEN_RDONLY`, `OPEN_RDWR` | `OPEN_RDONLY_PARALLEL` |
| `Cell<string>` I/O | Variable-length string arrays | — |
| Error paths | Duplicate dataset, missing dataset, nonexistent file | — |
| `Ascii` / `CsvFile` | — | Future plan |
| `InputFile` / `XML` | — | Future plan |

---

## Key Design Observations

### Everything Is `#ifdef BELFEM_HDF5`

All HDF5 functionality compiles to no-ops without the `BELFEM_HDF5` define. The entire test file should be wrapped in `#ifdef BELFEM_HDF5`.

### HDF5 Class Is a High-Level Wrapper

The `HDF5` class wraps the HDF5 C API. It manages file lifetime (constructor opens, destructor closes), active group tracking, and provides overloaded `save_data`/`load_data` for all supported types. Tests should use this high-level interface, not the `hdf5::` namespace functions directly.

### Matrix Storage Is Row-Major With Copy

`save_matrix_to_file` copies the BELFEM column-major `Matrix<T>` into a row-major `T**` buffer before writing. `load_matrix_from_file` reads into a row-major buffer and copies back element-by-element. This means the round-trip should preserve values exactly, but the internal layout differs from BELFEM's column-major format.

### Save Errors on Duplicate Labels

Every `save_*` function checks `!dataset_exists(aFileID, aLabel)` via `BELFEM_ERROR` (always active). Saving the same label twice throws. This is testable.

### Load Errors on Missing Labels

Every `load_*` function checks `dataset_exists(aFileID, aLabel)` via `BELFEM_ERROR` (always active). Loading a nonexistent label throws. This is testable.

### Test Files Must Be Cleaned Up

Each test should create a temporary file, write data, close, reopen, read, verify, and delete. Use a unique filename per test to avoid collisions.

---

## Floating-Point Comparison

```cpp
namespace
{
    const belfem::real tEps = 1e-12;
}
```

---

## Test File Structure

```
tests/io/
├── test_HDF5.cpp                # HDF5 class round-trip tests
```

---

## 1. Scalar Round-Trips

**File:** `test_HDF5.cpp`

### 1.1 Tests `[semantic]`

| Test Name | Type | What It Verifies |
|---|---|---|
| `ScalarSintRoundTrip` | `sint` | Save -42, load, verify == -42 |
| `ScalarUintRoundTrip` | `uint` | Save 12345, load, verify == 12345 |
| `ScalarLuintRoundTrip` | `luint` | Save large value, load, verify exact |
| `ScalarRealRoundTrip` | `real` | Save 3.14159265358979, load, verify within tEps |
| `ScalarBoolTrueRoundTrip` | `bool` | Save true, load, verify == true |
| `ScalarBoolFalseRoundTrip` | `bool` | Save false, load, verify == false |

---

## 2. String Round-Trips

### 2.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `StringRoundTrip` | Save "Hello BELFEM", load, verify exact match |
| `StringFromCharPointer` | `save_data(label, "literal")` → load as string, verify |

**Known limitation:** Empty string round-trip is not supported. The HDF5 C API rejects
zero-length fixed-size strings (`H5Tset_size` requires size > 0). The implementation in
`HDF5_Tools.cpp:152` passes `aValue.length()` directly to `H5Tset_size` with no guard.
Callers should use a sentinel value or `Cell<string>` for potentially-empty strings.

---

## 3. Vector Round-Trips

### 3.1 Tests `[semantic]`

| Test Name | Type | What It Verifies |
|---|---|---|
| `VectorSintRoundTrip` | `Vector<sint>` | Save {-1, 0, 1, 42}, load, verify element-by-element |
| `VectorUintRoundTrip` | `Vector<uint>` | Save {1, 2, 3, 4, 5}, load, verify |
| `VectorLuintRoundTrip` | `Vector<luint>` | Save large values, load, verify |
| `VectorRealRoundTrip` | `Vector<real>` | Save known values, load, verify within tEps |
| `VectorEmptyRoundTrip` | `Vector<real>` | Save empty vector (length 0), load, verify length == 0 |

---

## 4. Matrix Round-Trips

### 4.1 Tests `[semantic]`

| Test Name | Type | What It Verifies |
|---|---|---|
| `MatrixSintRoundTrip` | `Matrix<sint>` | Save 3×4 matrix, load, verify all elements |
| `MatrixUintRoundTrip` | `Matrix<uint>` | Save 2×3 matrix, load, verify |
| `MatrixRealRoundTrip` | `Matrix<real>` | Save 4×5 matrix with known values, load, verify within tEps |
| `MatrixRealDimensionsPreserved` | `Matrix<real>` | Verify n_rows() and n_cols() match after load |
| `MatrixSingleElement` | `Matrix<real>` | 1×1 matrix round-trip |

---

## 5. Cell of Strings Round-Trip

### 5.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CellStringRoundTrip` | Save {"alpha", "beta", "gamma"}, load, verify size and content |
| `CellStringEmpty` | Save empty Cell, load, verify size == 0 |
| `CellStringWithWhitespace` | (idea from ChatGPT) Save {"hello world", "foo bar"}, load, verify spaces preserved. Exercises the variable-length string HDF5 path which is hand-written and more bug-prone than the templated scalar/vector/matrix paths. |

---

## 5b. Raw Array Round-Trip

### 5b.1 Tests `[semantic]`

These test the `load_data(label, T*, memorySize)` overloads that load directly into caller-allocated buffers.

| Test Name | What It Verifies |
|---|---|
| `RawArrayIndexRoundTrip` | Save `Vector<index_t>`, load into raw `index_t*` buffer via `load_data(label, ptr, size)`, verify values |
| `RawArrayRealRoundTrip` | Same for `real*` |

**Note:** The raw-array integer overload uses `index_t*` (not `uint*`). The header
(`cl_HDF5.hpp`) originally declared `uint*` but the source (`cl_HDF5.cpp`) defined
`index_t*`. The header was corrected to match the source during implementation.

---

## 6. Group Management

### 6.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CreateAndSelectGroup` | `create_group("grp")`, data saved in group, reopen, `select_group("grp")`, data loaded |
| `NestedGroups` | Create "outer", create "inner" inside it, save data at each level, reopen, navigate and verify |
| `MultipleGroupsSameLevel` | Create "grp1" and "grp2", save different data in each, reopen, verify isolation |
| `CloseActiveGroupNavigatesBack` | Open nested group, close_active_group, verify parent is active via `tree()` AND by loading data from the parent group |
| `CloseTreeReturnsToRoot` | (idea from ChatGPT) After navigating multiple levels, `close_tree()` resets state to root; verified via `tree()` AND by loading data from root |
| `TreePathTrackingCorrect` | (idea from ChatGPT/Gemini) `tree()` returns correct path string during nested navigation |

---

## 7. File Modes

### 7.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `FileModeNewCreatesFile` | `HDF5(path, FileMode::NEW)` → file exists after close |
| `FileModeNewTruncatesExisting` | Create file with data, reopen with `NEW`, old data gone |
| `FileModeOpenRdonly` | Open existing file as RDONLY → data readable |
| `FileModeOpenRdwr` | Open existing file as RDWR → can add new datasets |

---

## 8. Error Paths

### 8.1 Tests `[semantic]`

All use `BELFEM_ERROR` (always active) — NOT debug-only.

| Test Name | What It Verifies |
|---|---|
| `EmptyPathThrows` | (idea from ChatGPT) `HDF5("", FileMode::NEW)` → `BELFEM_ERROR` fires |
| `SaveDuplicateLabelThrows` | Save "x" twice → `BELFEM_ERROR` fires |
| `LoadMissingLabelThrows` | Load nonexistent label → `BELFEM_ERROR` fires |
| `OpenNonexistentFileRdonlyThrows` | `HDF5(path, OPEN_RDONLY)` on missing file → `BELFEM_ERROR` |
| `OpenNonexistentFileRdwrThrows` | (idea from ChatGPT) `HDF5(path, OPEN_RDWR)` on missing file → `BELFEM_ERROR` |
| `CreateDuplicateGroupThrows` | `create_group("grp")` twice → `BELFEM_ERROR` |
| `SelectMissingGroupThrows` | (idea from ChatGPT) `select_group("nonexistent")` → `BELFEM_ERROR` |

---

## 9. Mixed Types in Same File

### 9.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `MixedTypesInRootGroup` | Save scalar, vector, matrix, string under root → all loadable |
| `MixedTypesInSubGroup` | Same in a named group → all loadable |

---

## 10. Source Bugs Found and Fixed

Five bugs (BUG-IO-1 through BUG-IO-3c) were found and fixed during test implementation.
All are resolved. Full details in `devlog/dl20260324_test_suite_bugs.md` under "I/O Module (HDF5)".

---

## 11. What We Do NOT Test (Deferred)

- Parallel file I/O (`make_path_parallel`, parallel mode constructor)
- `Ascii` class file read/write
- `CsvFile` class
- `InputFile` / `Input_Section` / `Input_Settings`
- `XML` class
- Mesh-level I/O (GmshReader, HDF5Reader/Writer, ExodusWriter, VtkWriter)
- Cross-platform endianness (HDF5 handles this internally)

---

## 12. Implementation Notes for Claude Code

1. **All tests need `#ifdef BELFEM_HDF5`.** Without HDF5, everything compiles to no-ops.
2. **Use unique temp file per test.** Pattern: `std::string tPath = "test_hdf5_<testname>.hdf5";` Delete at start and end with `std::remove(tPath.c_str());`.
3. **The test pattern is always:** create file (NEW) → save data → close → reopen (OPEN_RDONLY) → load data → verify → close → delete file.
4. **BELFEM_ERROR is compiled in for every build.** Duplicate-label and missing-label tests stay unguarded; `EXPECT_THROW` works in release because each test `main` calls `belfem::assert::set_throw_on_error( true )`. "Always active" describes the *check*, not the reaction — without that call a release build would `MPI_Abort` instead of throwing. Assertion-tier tests are the separate case: guard those with `#if BELFEM_ASSERTIONS_ACTIVE`, exported by `assert.hpp`.
5. **Matrix values are compared element-by-element** with `EXPECT_NEAR` for `real` or `EXPECT_EQ` for integer types.
6. **The `HDF5` destructor calls `close()` automatically.** Scoping with `{}` blocks is a clean pattern for separating write and read phases.
7. **Node IDs and element data from BUG-S1** — the Spline module's save/load field name mismatch (BUG-S1) is NOT an HDF5-level bug. It's a Spline-level bug where the wrong label strings are passed to working HDF5 functions.
8. **Bool save/load uses `hbool_t` internally.** The `bool` round-trip tests verify the cast between C++ `bool` and HDF5 `hbool_t` works correctly.

---

## 13. Codex Audit Checklist

When reviewing Claude Code's test implementation, verify:

- [ ] `#ifdef BELFEM_HDF5` wraps the entire test file
- [ ] Each test creates a unique temp file and cleans up
- [ ] Scalar round-trips cover sint, uint, luint, real, bool
- [ ] Vector round-trips cover all 4 numeric types + empty
- [ ] Matrix round-trips verify both values and dimensions
- [ ] String round-trip tested with non-empty strings (empty string unsupported — HDF5 API limitation)
- [ ] Group navigation tested: create, select, nested, close (with I/O verification, not just path string)
- [ ] File mode tests verify NEW truncation and RDONLY read
- [ ] `BELFEM_ERROR` paths tested unguarded; any `BELFEM_ASSERT` path guarded with `#if BELFEM_ASSERTIONS_ACTIVE`
- [ ] `EXPECT_NEAR` for real comparisons, `EXPECT_EQ` for integers
- [ ] BELFEM naming conventions (`t` prefix for locals)
- [ ] Deferred items documented in test file header
