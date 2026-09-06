# Bugs Found and Fixed During Test Suite Implementation

**Date:** 2026-03-23 / 2026-03-24
**Purpose:** Permanent record of all bugs discovered during the BELFEM test suite implementation campaign
**Module:** All (containers through mesh)

---

## Summary

24 source bugs were identified across 11 modules during test development. All were fixed except BUG-M4 (platform-specific, documented).

---

## Linalg Module

### BUG-L1: `Matrix::submat()` bounds checking
- **Files:** `cl_AR_Matrix.hpp`, `cl_BZ_Matrix.hpp`
- **Problem:** All 4 `BELFEM_ASSERT` calls compared row indices against `n_cols()` and reused `aLastRow` where `aLastCol` was intended. 16 assert lines across 4 overloads.
- **Fix:** Row checks now use `n_rows()`, column checks use `aFirstCol`/`aLastCol`.
- **Regression:** `MatrixDebug.SubmatRowOutOfBoundsThrows`, `MatrixDebug.SubmatColOutOfBoundsThrows`

### BUG-L2: `r2()` matrix version loop bounds
- **File:** `fn_r2.hpp`
- **Problem:** Variable names `m`/`n` swapped, both loops iterated `n` (cols) for row index.
- **Fix:** Corrected variable names and loop bounds.
- **Regression:** `R2.R2MatrixVersion` (3x5 non-square matrix)

---

## Math Tools Module

### BUG-M1: `symratiospace()` midpoint formula
- **File:** `fn_symrationspace.hpp`
- **Problem:** `aX(tN) = 0.5 * (aXmax - aXmin)` missing `aXmin +` prefix. Masked when `aXmin = 0`.
- **Fix:** `aX(tN) = aXmin + 0.5 * (aXmax - aXmin)`

### BUG-M2: `find_interval()` stale output
- **File:** `fn_find_interval.hpp`
- **Problem:** `aIndex`/`aXi` assignments after `break` — never executed.
- **Fix:** Moved assignments before break, using `i`/`k` interval endpoints.

### BUG-M3: `quadratic_gradient()` commented-out assertion
- **File:** `fn_quadratic_gradient.hpp`
- **Problem:** `BELFEM_ASSERT(aX.length() == aF.length())` commented out.
- **Fix:** Uncommented.

### BUG-M4: `rotation_matrix()` Euler-angle Blaze path (NOTED, NOT FIXED)
- **File:** `fn_rotation_matrix.hpp`
- **Problem:** Blaze path assumes 4x4 padding for 3x3 matrix. Platform-specific.
- **Status:** Tests use `operator(i,j)` which abstracts away padding.

### BUG-M5: `cardano()` D=0 branch
- **File:** `fn_cardano.hpp`
- **Problem:** D=0 returned single root `-b/(3a)`, wrong for double-root case.
- **Fix:** Full contract redesign — returns distinct real roots only, sorted. Quadratic D=0 returns 1 root, cubic D=0 returns 1 (triple) or 2 (double+simple) roots with corrected formulas.

---

## Comm Module

### BUG-C1: `comm_tag()` modulo-by-zero
- **File:** `commtools.cpp`
- **Problem:** `% gComm.max_tag()` with `max_tag() == 0` in non-MPI mode.
- **Fix:** Early return `if( gComm.max_tag() == 0 ) return 0;`

---

## Graph Module

### BUG-G1: `graph::sort()` is a no-op
- **File:** `fn_Graph_sort.cpp`
- **Problem:** `std::sort(.begin(), .begin(), ...)` — empty range.
- **Fix:** Second `.begin()` → `.end()`.

### BUG-G3: `metis_partition()` debug counter sized by MPI ranks instead of partitions
- **File:** `fn_Graph_METIS.cpp:297-303`
- **Problem:** Debug block `Vector<index_t> tCounters(gComm.size(), 0)` used `gComm.size()` (1 in serial) to size a counter indexed by `owner()` (0..N-1). With N=2 partitions, `++tCounters(1)` was out-of-bounds on a size-1 vector.
- **Fix:** `gComm.size()` → `aNumPartitions`.
- **Regression:** `GraphMetis.MetisPartition2`

### BUG-G2: `find_pseudo_peripheral_node()` confuses width with depth
- **File:** `fn_Graph_find_pseudo_peripheral_node.cpp`
- **Problem:** `bfs()` returns max width, code used it as max depth.
- **Fix:** Replaced with explicit max-level scan after BFS.
- **Build fix:** File was missing from `src/math/graph/CMakeLists.txt`.

### Empty graph crash in `find_connected_partitions()`
- **File:** `fn_Graph_find_connected_partitions.cpp`
- **Problem:** `.front()` on empty vector for empty graph.
- **Fix:** Early return `if( aGraph.size() == 0 ) return 0;`

---

## Tensor Module

### BUG-T1: Move assignment self-assignment crash
- **File:** `cl_Tensor.hpp`
- **Problem:** `free(mData)` before `mData = aTensor.mData` on self-move.
- **Fix:** `if( this == &aTensor ) return *this;` in both copy and move assignment.

### BUG-T2: `operator==` reads past second tensor's buffer
- **File:** `cl_Tensor.hpp`
- **Problem:** Used `aA.capacity()` without checking `aB.capacity()`.
- **Fix:** `if( aA.capacity() != aB.capacity() ) return false;`

---

## Sparse Module

### mIndexFunction uninitialized in dense/external-arrays constructors
- **File:** `cl_SpMatrix.cpp`
- **Problem:** Dense and external-arrays constructors did not call `set_indexing_base()`.
- **Fix:** Added `this->set_indexing_base(SpMatrixIndexingBase::Cpp);` to both.

### BUG-S1: `multiply()` indexing side-effect
- **File:** `cl_SpMatrix.cpp`
- **Problem:** `set_indexing_base(Fortran)` never restored after BLAS call.
- **Fix:** Save/restore pattern. Added `BELFEM_ASSERT(indexing_base() == 0)` to `operator()`.

### BUG-S2: SolverParameters spelling typos
- **Files:** `cl_SolverParameters.hpp/.cpp`, `cl_SolverPETSC.cpp`, `cl_SolverSTRUMPACK.cpp`
- **Problem:** `mUseInitualGuess`, `set_use_inital_guess`, `use_inital_guess`.
- **Fix:** Renamed to `mUseInitialGuess`, `set_use_initial_guess`, `use_initial_guess`.

### `indexing_base()` null pointer on empty matrix
- **File:** `cl_SpMatrix.hpp`
- **Fix:** `return mPointers != nullptr ? mPointers[0] : 0;`

### `create_coo_indices()` indexing side-effect
- **File:** `cl_SpMatrix.cpp`
- **Fix:** Save/restore pattern matching `multiply()`.

---

## Mesh Module

### BUG-ME1: `Mesh::node(i,j,k)` asserts wrong dimension
- **File:** `cl_Mesh.hpp`
- **Problem:** 3D accessor asserted `number_of_dimensions() == 2`.
- **Fix:** `== 2` → `== 3`, message `"2d"` → `"3d"`.

### BUG-ME2: `QUAD4TS` wrong facet count
- **Files:** `cl_Element_QUAD4TS.hpp`, `cl_Element_Factory.cpp`
- **Problem:** Template `<4,4,2,1,1>` — 1 facet, but top/bottom shell has 2.
- **Fix:** Template `<4,4,2,2,1>`.

---

## Container Module (Source Improvements, Not Bugs)

These were code quality improvements found by Junie during review:
- `Cell::print()` — `(int)` cast replaced with `std::ostringstream` for generic types
- `DynamicBitset` — added `#include <bit>` for C++20 `std::popcount`
- `DynamicBitset` — FNV-1a comment corrected
- `ShiftRegister` — added `is_shift_register_safe<T>` trait
- `Bitset` — added `const data()` overload
- `Map`/`OrderedMap` — `find()` changed to `const Key&`
- `Genome` — added `#include "constants.hpp"`

---

## Spline Module (DEFERRED)

### BUG-S1 (spline): `save()`/`load()` field name mismatch (FIXED)
- **File:** `cl_Spline.cpp`
- **Problem:** `save(hid_t&)` wrote fields as `"min"`, `"max"`, `"n"`, `"data"` but `load(hid_t&)` read `"xmin"`, `"xmax"`, `"npoints"`, `"coeffs"`. Additionally, `load()` transposed the data but `save()` did not.
- **Fix:** Updated `save()` to use `load()`'s field names (`"xmin"`, `"xmax"`, `"npoints"`, `"coeffs"`) and added `trans(mData)` before saving. This aligns with `save_to_database()` which already used the correct names and transposition.
- **Note:** Any HDF5 files saved with the old `save()` will NOT load correctly. The `save_to_database()` path was always consistent and is unaffected.

### Entropy mode assert-before-set (FIXED)
- **File:** `cl_Spline.cpp`, `create_entropy()`
- **Problem:** `mExtraMode = Entropy` was set at line 667, AFTER `this->entropy(aXref)` was called at line 660. The `entropy()` method asserts `mExtraMode == Entropy`, so this was a guaranteed assert failure when constructing a spline with entropy mode.
- **Fix:** Moved `mExtraMode = spline::ExtraMode::Entropy` before the `this->entropy(aXref)` call.
- **Discovery:** Found during spline test suite implementation when `EntropyModeConstruction` test failed.

---

## I/O Module (HDF5)

### BUG-IO-1: Header/source prototype mismatch for raw-array load
- **File:** `cl_HDF5.hpp:295`, `cl_HDF5.cpp:847`
- **Problem:** Header declared `load_data(label, uint*, memorySize)` but source defined `load_data(label, index_t*, memorySize)`. With `BELFEM_INT64`, `index_t` is `uint64_t` while `uint` is `uint32_t` — linker error.
- **Fix:** Changed header declaration to `index_t*` to match the source definition.
- **Regression:** `RawArrayIndexRoundTrip`

### BUG-IO-2: `create_group()` closed parent handle and lost tree state
- **File:** `cl_HDF5.cpp:158-201`
- **Problem:** Two interrelated bugs: (1) `create_group()` called `H5Gclose(mActiveGroup)` before using `mActiveGroup` as the parent for `H5Gcreate2`, invalidating the parent handle for nested group creation. (2) `create_group()` did not push to `mTree`/`mTreeLabels`, so `close_active_group()` and `tree()` could not track groups created this way.
- **Fix:** Removed the incorrect `H5Gclose` call. Added `mTreeLabels.push(aLabel)`, `mTree.push(mActiveGroup)`, and `mActiveGroupLabel = hdf5::create_tree(...)` to match the state management in `select_group()`.
- **Regression:** `NestedGroups`, `CloseActiveGroupNavigatesBack`, `CloseTreeReturnsToRoot`, `TreePathTrackingCorrect`

### BUG-IO-3: Error message typos (3 instances)
- **Files:** `cl_HDF5.cpp:580`, `HDF5_Tools.cpp:187`, `HDF5_Tools.cpp:245`
- **Problems:** (a) "Error reating dataset" — missing 'c'. (b) Trailing dash in save-string error message. (c) "trying to store" used inside `load_string_from_file`.
- **Fix:** (a) "Error creating dataset". (b) Removed dash. (c) "trying to load".

### Known limitation: empty string round-trip unsupported
- **File:** `HDF5_Tools.cpp:152`
- **Problem:** `H5Tset_size(tDataType, aValue.length())` — HDF5 C API rejects size 0 for fixed-size strings.
- **Status:** Documented in test plan. Not a BELFEM bug — HDF5 C library limitation. Callers should use a sentinel or `Cell<string>` for potentially-empty strings.

---

## ODE Module

### BUG-ODE-1: DOP853 missing from Integrator dispatcher (FIXED BY USER)
- **File:** `cl_ODE_Integrator.cpp:29-47`, `en_ODE_Type.hpp`
- **Problem:** The `Integrator` constructor only handled `Type::RK45`. DOP853 existed as standalone functions but was unreachable via the high-level API.
- **Fix:** User added `Type::DOP853 = 1` to the enum and a `case(Type::DOP853)` branch in the constructor that calls `DOP853_init` and assigns the `DOP853` function pointer.
- **Regression:** `DOP853.ViaIntegratorExponentialDecay`, `DOP853.ViaIntegratorHarmonic`, `Integrator.TypeAccessorDOP853`

### BUG-ODE-2: `Integrator::time()` was a dead property
- **File:** `cl_ODE_Integrator.cpp:53-65`
- **Problem:** `mTime` was stored and had a public accessor, but `step()` never read or wrote it. Users could set `time()` expecting it to reflect integration state, but it remained at whatever value was assigned.
- **Fix:** Added `mTime = aT` after the integration function returns in `step()`. Now `time()` tracks the actual integration time.
- **Regression:** `Integrator.TimeDoesNotDriveStep` (updated to verify sync)

### BUG-ODE-3: Persistent TRAPPED status after step rejection
- **Files:** `fn_ODE_RK45.cpp:92`, `fn_ODE_DOP853.cpp:189`
- **Problem:** If a step was initially clamped to `aTmax` (setting `aStatus = TRAPPED`) but then rejected due to excessive error, the status remained TRAPPED on retry even when the reduced step size no longer reached `aTmax`.
- **Fix:** Added `aStatus = Status::OK` reset at the top of each iteration, before the trap check. Applied to both RK45 and DOP853.
- **Regression:** Existing `StatusTrapped` tests still pass (the fix only affects the rejected-then-retried path).
