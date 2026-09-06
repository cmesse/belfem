# BELFEM Sparse Module Tests — Detailed Plan

**Date:** 2026-03-22
**Purpose:** Method-level test matrix for the sparse module (`src/sparse/`)
**Depends on:** `tests_0_strategy.md` (conventions), `tests_2_linalg.md` (Vector/Matrix), `tests_7_graph.md` (Graph/Vertex)
**Confidence:** High on SpMatrix core. Medium on solver backends (conditional compilation). Low on DistMatrix (complex MPI+permutation interactions).

---

## Module Overview

| Layer | Content | Test Priority |
|---|---|---|
| `SpMatrix` class | CSR/CSC sparse matrix. Manual malloc/free. Dense constructor ("for testing"). Binary-search element access. Function-pointer indexing dispatch. | **Highest** |
| `Solver` class | High-level façade: dispatches to backend wrappers, lazy initialization, parameter sync | High |
| `solver::Wrapper` hierarchy | UMFPACK, MUMPS, PARDISO, PETSc, STRUMPACK — all conditionally compiled | Medium (smoke tests) |
| `SolverParameters` | Configuration: solver type, reordering, compression, Krylov method, tolerances | Medium |
| `DistMatrix` / `DistMatrixCSR` / `DistMatrixAIJ` | Distributed sparse layouts with optional METIS permutation | Low (MPI integration) |
| Utilities | `create_graph_from_matrix`, `compute_permutation`, `fn_rcond`, ARPACK | Low |
| Enums | `SolverType`, `SymmetryMode`, `SpMatrixType`, `Preconditioner`, `KrylovMethod`, etc. | Low |

---

## Key Design Observations

### Dense-Matrix Constructor (Testing Gateway)

`SpMatrix(const Matrix<real>&, SpMatrixType)` is explicitly documented as "constructor for testing purposes." It counts nonzeros in a dense matrix, builds CSR or CSC structure, and populates values. All SpMatrix tests should use this constructor as the primary entry point.

### Operator() Asymmetry

- **Writable** `operator()(row, col)`: asserts (`BELFEM_ASSERT`) if position is a structural zero. Returns reference to value.
- **Const** `operator()(row, col)`: returns `mZero` (member = 0.0) for structural zeros. No assertion.

This means writing to a zero position is a debug-only failure, while reading a zero position silently returns 0.0. Both behaviors need explicit tests.

### Indexing Base Mutation

`multiply()` internally calls `set_indexing_base(Fortran)` before delegating to FSPBLAS/MKL. After multiply returns, the matrix remains in Fortran indexing. If element access via `operator()` is attempted after multiply, the function-pointer dispatch uses the Fortran-based binary search, which expects 1-based indices in the pointer/index arrays. This should be tested: verify element access works correctly both before and after multiply.

### sort_entries()

The external-array constructor and HDF5 loader call `sort_entries()` automatically. This sorts index-value pairs within each row (CSR) or column (CSC) to enable binary search in `index()`. Tests with deliberately unsorted input verify this contract.

---

## Floating-Point Comparison

```cpp
namespace
{
    const belfem::real tEps = 1e-12 ;  // for exact-in-theory results
    const belfem::real tTol = 1e-9 ;   // for solver residuals
}
```

---

## Test File Structure

```
tests/sparse/
├── test_SpMatrix.cpp              # SpMatrix construction, access, multiply, transpose
├── test_Solver.cpp                # Solver façade, backend smoke tests
```

Backend-specific integration tests (MUMPS, STRUMPACK, PETSc) should be in `test_Solver.cpp` wrapped in `#ifdef` guards. DistMatrix tests belong in a future MPI integration file.

---

## 1. SpMatrix Construction

**File:** `test_SpMatrix.cpp`

### 1.1 Dense Constructor `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `ConstructCSRFromDense` | Build from 4×4 dense matrix → `type() == CSR`, `n_rows() == 4`, `n_cols() == 4` |
| `ConstructCSCFromDense` | Same dense matrix → `type() == CSC` |
| `NonzeroCountCorrect` | Known dense matrix with specific zero pattern → `number_of_nonzeros()` matches hand count |
| `DiagonalMatrixCSR` | Identity-like diagonal matrix → nnz = N, each row has exactly 1 entry |
| `FullDenseMatrixCSR` | Fully dense 3×3 matrix → nnz = 9 |
| `EmptyMatrixCSR` | All-zero dense matrix → nnz = 0, no crash |
| `RectangularCSR` | 3×5 dense matrix → correct dimensions and structure |
| `RectangularCSC` | 5×3 dense matrix → correct dimensions and structure |

### 1.2 Dense Round-Trip `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CSRReadBackMatchesDense` | For every (i,j): `const SpMatrix(i,j) == dense(i,j)` including structural zeros |
| `CSCReadBackMatchesDense` | Same for CSC |
| `CSRAndCSCAgree` | Build CSR and CSC from same dense matrix → `operator()(i,j)` returns same values for all (i,j) |

### 1.3 Pointer/Index Structure `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CSRPointersMonotone` | `pointers()[i] <= pointers()[i+1]` for all i |
| `CSRPointersLastEqualsNNZ` | `pointers()[n_rows()] == number_of_nonzeros()` |
| `CSCPointersLastEqualsNNZ` | `pointers()[n_cols()] == number_of_nonzeros()` |
| `CSRIndicesSorted` | Within each row, column indices are strictly increasing |
| `CSCIndicesSorted` | Within each column, row indices are strictly increasing |
| `IndicesReturnsCorrectArray` | For CSC, `indices()` returns row indices; for CSR, column indices |

---

## 2. SpMatrix Indexing and Access

### 2.1 Element Access `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `WritableAccessNonzero` | `tM(i,j) = 5.0` for a structural nonzero → value stored and readable |
| `ConstAccessNonzero` | Const ref reads correct value |
| `ConstAccessStructuralZero` | Const ref at zero position returns 0.0 |
| `FillSetsAllValues` | `fill(3.14)` → all nonzero entries are 3.14, structural zeros still return 0.0 |
| `IndexReturnsNNZForZeroPosition` | `index(i,j)` for a structural zero returns `number_of_nonzeros()` |
| `IndexReturnValidForNonzero` | `index(i,j)` for a nonzero returns value < `number_of_nonzeros()` |

### 2.2 Element Access `[debug]`

| Test Name | What It Verifies |
|---|---|
| `WritableAccessStructuralZeroThrows` | `tM(i,j) = val` at a zero position → assertion |
| `DataIndexOutOfBoundsThrows` | `data(number_of_nonzeros())` → assertion |

### 2.3 Indexing Base Conversion `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CppToFortranRoundTrip` | Convert to Fortran → `indexing_base() == 1` → convert back to C++ → `indexing_base() == 0` |
| `FortranBasePointersShifted` | After Fortran conversion, `pointers()[0] == 1` |
| `ElementAccessAfterFortranConversion` | `operator()(i,j)` still returns correct values in Fortran mode |
| `DoubleConversionIsIdempotent` | Two consecutive C++→C++ conversions don't corrupt state |

---

## 3. SpMatrix Operations

### 3.1 Multiply `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `MultiplyCSRIdentity` | Identity matrix × vector → same vector |
| `MultiplyCSRTridiagonal` | Known tridiagonal × known vector → known result |
| `MultiplyCSCTridiagonal` | Same test in CSC format |
| `MultiplyAlphaBeta` | `y = 2*A*x + 3*y` → verify against hand calculation |
| `MultiplyTransposed` | `y = Aᵀ*x` → verify against dense transpose multiplication |
| `OperatorStarMatchesMultiply` | `A * x` (operator) matches `A.multiply(x, y)` |
| `MultiplyRectangular` | 3×5 sparse matrix × length-5 vector → length-3 result |
| `ElementAccessAfterMultiply` | After multiply, `operator()(i,j)` still returns correct values (tests indexing base mutation recovery) |

### 3.2 Transpose `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `TransposeSwapsType` | CSR → transpose → CSC (and vice versa) |
| `TransposeSwapsDimensions` | 3×5 matrix → transpose → 5×3 |
| `TransposePreservesValues` | `A(i,j) == Aᵀ(j,i)` for all entries |
| `DoubleTransposeRoundTrip` | Transpose twice → original structure and values |

### 3.3 COO Indices `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CreateCOOForCSR` | After `create_coo_indices()`, `rows()` is non-null and `have_coo_indices() == true` |
| `CreateCOOForCSC` | After `create_coo_indices()`, `cols()` is non-null |
| `COOIndicesCorrect` | For each entry k: `(rows()[k], cols()[k])` identifies a valid nonzero with value `data()[k]` |
| `FreeCOOIndices` | After `free_coo_indices()`, the supplementary array is null and `have_coo_indices() == false` |

### 3.4 sort_entries `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `SortEntriesFixesUnsortedInput` | Build SpMatrix from external arrays with deliberately unsorted indices → after construction, `index(i,j)` works correctly via binary search |
| `SortEntriesPreservesValues` | Values move with their indices during sort |

---

## 4. Graph and Permutation Utilities

**File:** `test_SpMatrix.cpp`

### 4.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CreateGraphFromMatrix` | Build SpMatrix → `create_graph_from_matrix` → graph has correct vertex count and adjacency |
| `ComputePermutationIsValid` | After graph reordering → forward permutation is a bijection of 0..N-1 |
| `PermutedSolveMatchesOriginal` | Build permuted SpMatrix from reordered graph → `Q * permute(x) ≈ permute(A * x)` |

---

## 5. Solver Façade

**File:** `test_Solver.cpp`

### 5.1 Construction `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `SolverTypeStored` | `Solver(SolverType::UMFPACK).type() == SolverType::UMFPACK` (if linked) |
| `SolverWrapperNotNull` | `wrapper()` returns non-null pointer after construction |
| `SolverFromParameters` | `Solver(SolverParameters(type))` creates matching type |

### 5.2 Known-Solution Solve Pattern `[semantic]`

For each available backend, use the canonical pattern from `solvertest.cpp`:

```
1. Build dense tridiagonal K (2 on diagonal, -1 off-diagonal)
2. Convert to SpMatrix
3. Choose known x: x(k) = k+1
4. Compute b = K * x (dense multiply)
5. Solve K * x_solved = b
6. Verify x_solved ≈ x within tTol
```

| Test Name | Backend Guard | What It Verifies |
|---|---|---|
| `UMFPACKSolveTridiagonal` | `#ifdef BELFEM_SUITESPARSE` | Full solve round-trip |
| `MUMPSSolveTridiagonal` | `#ifdef BELFEM_MUMPS` | Full solve round-trip |
| `PARDISOSolveTridiagonal` | `#ifdef BELFEM_PARDISO` | Full solve round-trip |
| `PETScSolveTridiagonal` | `#ifdef BELFEM_PETSC` | Full solve round-trip (needs MPI init) |
| `STRUMPACKSolveTridiagonal` | `#ifdef BELFEM_STRUMPACK` | Full solve round-trip (needs MPI init) |

### 5.3 Re-Solve with Changed Values `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `ReSolveReuseFactorization` | Solve once, change matrix values (not structure), solve again → correct answer. Tests MUMPS job-5 vs job-6 distinction. |

### 5.4 Solver Lifecycle `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `LazyInitialization` | `wrapper()->is_initialized() == false` before first solve, `true` after |
| `FreeResetsState` | After `free()`, `is_initialized() == false` |

### 5.5 Solver `[debug]`

| Test Name | What It Verifies |
|---|---|
| `UnlinkedBackendThrows` | Creating a solver for an unlinked backend → `BELFEM_ERROR` fires |

---

## 6. SolverParameters

**File:** `test_Solver.cpp`

### 6.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `DefaultParametersValid` | `SolverParameters(type)` has sensible defaults for each type |
| `ReorderingMethodStored` | `set_reordering_method` → getter returns same value |
| `SymmetryModeStored` | `set_symmetry_mode` → getter returns same value |
| `InitialGuessFlag` | `set_use_inital_guess(true)` → `use_inital_guess() == true` (**NOTE:** API is misspelled in source — `inital` not `initial`. Use actual API names.) |

---

## 7. ARPACK Error Mapping (Low Priority)

### 7.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `ArpackCheckSuccessOnZero` | `arpack::check(0)` → no error |
| `ArpackCheckKnownCodes` | `arpack::check(-1)`, `arpack::check(-3)`, etc. → correct `BELFEM_ERROR` messages |

---

## 8. What We Do NOT Test

- External solver library correctness (MUMPS, PETSc, STRUMPACK, etc.)
- HDF5 save/load (requires HDF5 linkage; separate integration concern)
- DistMatrix distribution and collection (requires MPI; future plan)
- FSPBLAS/MKL internal correctness
- Performance benchmarks
- Large-scale problems (>100 DOFs)

---

## 9. Implementation Notes for Claude Code

1. **Use the dense-matrix constructor** for all SpMatrix tests. It's the designed test entry point.
2. **Build small test matrices by hand.** The tridiagonal Laplacian (2 on diagonal, -1 off-diagonal) is the canonical test matrix for solver round-trips.
3. **CSR and CSC tests should be symmetric.** Every structural test done for CSR should have a CSC counterpart.
4. **The `multiply()` indexing base mutation** is a real concern. Test element access after multiply to confirm the matrix is still usable via `operator()`.
5. **Backend solver tests are smoke tests.** Verify the answer is correct, not that the solver is optimal. Use known solutions.
6. **Backend tests need `#ifdef` guards.** Each backend's tests must be wrapped in the corresponding `BELFEM_*` define.
7. **PETSc and STRUMPACK tests may need MPI.** If so, they belong in a separate MPI test binary (like the comm module's Tier 2).
8. **The re-solve test** is one of the highest-value solver tests. It verifies that BELFEM correctly distinguishes first factorization from value-only update, which is a BELFEM policy decision layered on top of the solver library.
9. **sort_entries test:** Build a SpMatrix from raw arrays with deliberately reversed column indices within a row. After construction (which calls sort_entries automatically), verify that binary search via `index()` works.

---

## 10. Codex Audit Checklist

When reviewing Claude Code's test implementation, verify:

- [ ] All SpMatrix tests use the dense-matrix constructor
- [ ] CSR and CSC both tested for construction, access, and multiply
- [ ] Round-trip test verifies every (i,j) entry including structural zeros
- [ ] Operator() asymmetry tested: writable asserts on zero, const returns 0.0
- [ ] Indexing base C++↔Fortran round-trip tested
- [ ] Element access after multiply tested (indexing base mutation)
- [ ] sort_entries tested with unsorted external input
- [ ] Solver tests use known-solution pattern (b = A*x_known, solve, verify)
- [ ] Each backend wrapped in correct `#ifdef`
- [ ] Re-solve test present for at least one backend
- [ ] Debug tests: writable zero access, data index bounds, unlinked backend
- [ ] `EXPECT_NEAR` for all floating-point comparisons
- [ ] BELFEM naming conventions (`t` prefix for locals)
