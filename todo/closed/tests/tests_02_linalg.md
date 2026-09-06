# BELFEM Linear Algebra Tests — Detailed Plan

**Date:** 2026-03-22
**Purpose:** Method-level test matrix for the `src/linalg/` module
**Depends on:** `tests_0_strategy.md` (read that first for conventions and architecture)
**Confidence:** High on API surface. Medium on exact numerical tolerances (backend-dependent). High on submat bounds bug (see Section 4 notes).

---

## Architecture Notes

The linalg module uses a **compile-time backend selection**: either Armadillo (`BELFEM_ARMADILLO`) or Blaze (`BELFEM_BLAZE`). Both backends expose the same BELFEM API through `cl_Vector.hpp` and `cl_Matrix.hpp`. All tests go through the **unified BELFEM API** — no `#ifdef` for backend selection in test files. The same test binary, compiled against whichever backend is active, validates that backend.

**Implication for CI:** The test suite should be compiled and run under both backend configurations. If only one backend is available in the CI environment, that's fine — the tests are backend-agnostic by design. But if both are available, run both.

**What we test vs. what we don't:**

- We **do** test: BELFEM wrapper contracts (construction, access, sizing, copy/move, `data()` pointer layout), BELFEM-authored logic (bounds checks, `inv2`/`inv3`, `polyfit` rescaling, `eigen` NaN policy, `crossmat`, `r2`), and correct dispatch of free functions (`dot`, `norm`, `det`, etc.) through small deterministic examples.
- We **do not** test: Armadillo/Blaze internal arithmetic correctness, BLAS/LAPACK numerical accuracy, expression template optimization. Those are the backend's responsibility.

---

## Floating-Point Comparison Policy

All floating-point comparisons use `EXPECT_NEAR` with explicit tolerance. Recommended constants:

```cpp
// At the top of each linalg test file
namespace
{
    const belfem::real tEps   = 1e-12 ;  // for exact-in-theory results (integer arithmetic via reals)
    const belfem::real tTol   = 1e-9 ;   // for solver residuals, eigenvalues, polyfit
}
```

For solver tests, prefer **residual checks** (`norm(A*x - b) < tTol`) over comparing `x` to a known solution, since pivoting and round-off can change the solution path without changing correctness.

For `inv2`/`inv3`, verify `A * inv(A) ≈ I` entrywise.

---

## Module Inventory

### Classes

| Class | Backend (Armadillo) | Backend (Blaze) | Notes |
|---|---|---|---|
| `Vector<T>` | `arma::Mat<T>` (Nx1) | `blaze::DynamicVector<T>` | Size via `length()` |
| `Matrix<T>` | `arma::Mat<T>` | `blaze::DynamicMatrix<T>` | Column-major storage |

### Free Functions — Math

| Function | Header | BELFEM logic? | Priority |
|---|---|---|---|
| `dot` | `fn_dot.hpp` | No (thin wrapper) | Phase 1 |
| `cross` | `fn_cross.hpp` | Assert length==3 | Phase 1 |
| `crossmat` | `fn_crossmat.hpp` | Yes (2D/3D, dust removal) | Phase 2 |
| `det` | `fn_det.hpp` | Thin wrapper + BELFEM overload | Phase 1 |
| `inv` | `fn_inv.hpp` | Thin wrapper | Phase 2 |
| `inv2` | `fn_inv2.hpp` | Yes (manual 2x2 formula) | Phase 1 |
| `inv3` | `fn_inv3.hpp` | Yes (manual 3x3 formula) | Phase 1 |
| `trans` | `fn_trans.hpp` | Thin wrapper | Phase 1 |
| `norm` | `fn_norm.hpp` | Thin wrapper (L2) | Phase 1 |
| `sum` | `fn_sum.hpp` | Thin wrapper | Phase 2 |
| `min` / `max` | `fn_min.hpp` / `fn_max.hpp` | Thin wrapper | Phase 2 |
| `linspace` | `fn_linspace.hpp` | Thin wrapper / manual (Blaze) | Phase 2 |
| `sort` | `fn_sort.hpp` | Thin wrapper / manual (Blaze) | Phase 2 |
| `unique` | `fn_unique.hpp` | Thin wrapper / manual (Blaze) | Phase 2 |
| `reverse` | `fn_reverse.hpp` | Thin wrapper | Phase 2 |
| `append` | `fn_append.hpp` | Backend-specific logic | Phase 2 |
| `combine` | `fn_combine.hpp` | Yes (2/3/4-vector variants) | Phase 2 |
| `polyfit` | `fn_polyfit.hpp` | Yes (rescaling logic) | Phase 3 |
| `polyval` | `fn_polyval.hpp` | Horner's method | Phase 3 |
| `dpolyval` | `fn_dpolyval.hpp` | Yes (derivative) | Phase 3 |
| `ddpolyval` | `fn_ddpolyval.hpp` | Yes (second derivative) | Phase 3 |
| `r2` | `fn_r2.hpp` | Yes (Vector and Matrix versions) | Phase 3 |
| `eigen` | `fn_eigen.hpp` | Yes (NaN policy for complex eigenvalues) | Phase 3 |

### Free Functions — LAPACK Wrappers

| Function | Header | Notes | Priority |
|---|---|---|---|
| `gesv` (Vector RHS) | `fn_gesv.hpp` | LU solve, mutates A and x | Phase 1 |
| `gesv` (Matrix RHS) | `fn_gesv.hpp` | Multi-RHS variant | Phase 2 |
| `posv` | `fn_posv.hpp` | Cholesky solve (SPD matrices) | Phase 3 |
| `lapack::gemm` | `fn_LAPACK_gemm.hpp` | Raw BLAS GEMM | Phase 3 |
| `lapack::gesv` | `fn_LAPACK_gesv.hpp` | Raw LAPACK interface | Phase 3 |
| `lapack::getrf` / `getri` | `fn_LAPACK_getrf/i.hpp` | LU factorization / inversion | Phase 3 |

### Operators

| Operator | Operands | Header | Priority |
|---|---|---|---|
| `+` | Vector±Vector, Vector±scalar | `op_VectorPlus.hpp` | Phase 1 |
| `-` | Vector-Vector, Vector-scalar, scalar-Vector | `op_VectorMinus.hpp` | Phase 1 |
| `*` | Vector*scalar, scalar*Vector | `op_VectorTimes.hpp` | Phase 1 |
| `/` | Vector/scalar | `op_VectorDivide.hpp` | Phase 1 |
| `==` | Vector==Vector, Vector==scalar | `op_VectorEqualEqual.hpp` | Phase 2 |
| `%=` | Vector element-wise multiply | (in Vector class) | Phase 2 |
| `+` | Matrix+Matrix, Matrix+scalar | `op_MatrixPlus.hpp` | Phase 1 |
| `-` | Matrix-Matrix, Matrix-scalar | `op_MatrixMinus.hpp` | Phase 1 |
| `*` | Matrix*Matrix, Matrix*Vector, Matrix*scalar | `op_MatrixTimes.hpp` | Phase 1 |

---

## Test File Structure

```
tests/linalg/
├── test_Vector.cpp              # Vector wrapper contracts
├── test_Matrix.cpp              # Matrix wrapper contracts
├── test_LinalgOperators.cpp     # Binary operators (+, -, *, /, ==)
├── test_LinalgFunctions.cpp     # Free functions (dot, cross, norm, det, inv, trans, etc.)
├── test_LinalgSolvers.cpp       # gesv, posv, LAPACK wrappers
├── test_LinalgPolynomials.cpp   # polyfit, polyval, dpolyval, ddpolyval, r2
```

---

## 1. Vector\<T\>

**File:** `test_Vector.cpp`
**Approach:** Typed tests over `real` (= `double`). If `float` is used anywhere in the codebase, add it. No `std::complex` — BELFEM uses real arithmetic.

### 1.1 Construction & Destruction `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `DefaultConstructorEmpty` | `length() == 0` |
| `SizedConstructor` | `Vector<real>(5)` → `length() == 5` |
| `SizedWithFillValue` | `Vector<real>(5, 3.0)` → all elements are `3.0` |
| `InitializerListMultipleElements` | `{1.0, 2.0, 3.0}` → `length() == 3`, values match |
| `InitializerListSingleElement` | `{42.0}` → `length() == 1`, `(0) == 42.0` |
| `CopyConstructorDeepCopies` | Modify copy, original unchanged |
| `MoveConstructor` | Source is in valid moved-from state |
| `CopyAssignment` | Deep copy, independent |
| `MoveAssignment` | Source is moved-from |
| `SelfCopyAssignment` | No corruption |
| `SelfMoveAssignment` | No corruption |
| `AssignFromScalar` | `tVec = 5.0` fills all elements with `5.0` |
| `AssignFromInitializerList` | `tVec = {1.0, 2.0}` resizes and sets values |

### 1.2 Memory & Access `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `DataPointerNonNull` | `data()` is non-null for sized vector |
| `DataPointerMatchesFirstElement` | `data()[0] == (0)` and `data()[i] == (i)` for all i |
| `ColumnMajorLayout` | `data()[i]` equals `operator()(i)` — verifies contiguous access |
| `ConstDataPointer` | `const Vector& → data()` returns same pointer |
| `VectorDataExposesBackend` | `vector_data()` returns reference to internal type |
| `ParenthesisReadWrite` | `(i) = x` then `(i) == x` |
| `ConstParenthesisAccess` | `const Vector& ref → ref(i)` compiles and returns correctly |
| `LengthMatchesSized` | `length()` equals constructor argument |
| `Iteration` | Range-based for loop visits `length()` elements in order |

### 1.3 Mutation `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `FillSetsAllElements` | `fill(7.0)` → all elements are `7.0` |
| `SetSizeChangesLength` | `set_size(10)` → `length() == 10` |
| `SetSizeWithValue` | `set_size(10, 3.0)` → `length() == 10`, all `3.0` |

### 1.4 Compound Operators `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `PlusEqualsScalar` | `tVec += 2.0` adds 2 to every element |
| `PlusEqualsVector` | `tA += tB` adds elementwise |
| `MinusEqualsScalar` | `tVec -= 2.0` subtracts 2 from every element |
| `MinusEqualsVector` | `tA -= tB` subtracts elementwise |
| `TimesEqualsScalar` | `tVec *= 3.0` scales all elements |
| `DivideEqualsScalar` | `tVec /= 2.0` divides all elements |
| `ElementwiseMultiplyEquals` | `tA %= tB` multiplies elementwise, lengths must match |

### 1.5 Access `[debug]`

| Test Name | What It Verifies |
|---|---|
| `OutOfBoundsThrows` | `tVec(length())` triggers assertion |
| `ElementwiseMultiplySizeMismatchThrows` | `tA(3) %= tB(5)` triggers assertion |

### 1.6 Print `[semantic]` (smoke test only)

| Test Name | What It Verifies |
|---|---|
| `PrintProducesOutput` | Capture stdout, verify non-empty |

---

## 2. Matrix\<T\>

**File:** `test_Matrix.cpp`
**Approach:** Tests use `real` type. Focus on BELFEM-specific behavior (bounds checks, views, `data()` layout).

### 2.1 Construction & Destruction `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `DefaultConstructorEmpty` | `n_rows() == 0`, `n_cols() == 0` |
| `SizedConstructor` | `Matrix<real>(3, 4)` → `n_rows() == 3`, `n_cols() == 4` |
| `SizedWithFillValue` | `Matrix<real>(3, 4, 2.0)` → all elements are `2.0` |
| `InitializerListConstructor` | `{{1,2,3},{4,5,6}}` → 2×3, correct values |
| `CopyConstructorDeepCopies` | Modify copy, original unchanged |
| `MoveConstructor` | Source is moved-from |
| `CopyAssignment` | Deep copy, independent |
| `MoveAssignment` | Source is moved-from |
| `SelfCopyAssignment` | No corruption |
| `SelfMoveAssignment` | No corruption |
| `AssignFromScalar` | `tMat = 5.0` fills all elements with `5.0` |

### 2.2 Memory & Access `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `DataPointerNonNull` | `data()` is non-null for sized matrix |
| `ColumnMajorLayout` | `data()[i + j * n_rows()] == (i, j)` for all i, j |
| `ConstDataPointer` | `const Matrix& → data()` returns same pointer |
| `MatrixDataExposesBackend` | `matrix_data()` returns reference to internal type |
| `ParenthesisReadWrite` | `(i, j) = x` then `(i, j) == x` |
| `ConstParenthesisAccess` | const ref access works |
| `CapacityMatchesProduct` | `capacity() == n_rows() * n_cols()` |

### 2.3 Sizing `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `FillSetsAllElements` | `fill(3.0)` → all elements are `3.0` |
| `SetSizeChanges` | `set_size(5, 7)` → dimensions match |
| `SetSizeWithValue` | `set_size(5, 7, 1.0)` → all elements `1.0` |

### 2.4 Row / Column / Submatrix Views `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `RowViewReadsCorrectly` | `row(i)` returns elements of row i |
| `RowViewWriteModifiesParent` | Modify via `row(i)`, verify parent matrix changed |
| `ColViewReadsCorrectly` | `col(j)` returns elements of column j |
| `ColViewWriteModifiesParent` | Modify via `col(j)`, verify parent matrix changed |
| `SubmatReadsCorrectly` | Extract 2×2 submatrix from 4×4, verify values |
| `SetRowFromVector` | `set_row(i, tVec)` → row i matches vector |
| `SetColFromVector` | `set_col(j, tVec)` → column j matches vector |
| `SubmatOnNonSquareMatrix` | Extract submatrix from 3×5 matrix — exercises potential bounds bug (see notes) |

**NOTE on `submat` bounds checking:** The Armadillo backend's `submat()` method has assert conditions that compare row indices against `n_cols()` instead of `n_rows()`, and reuse `aLastRow` where `aLastCol` is intended. The Blaze backend has the same pattern. A test with a non-square matrix (e.g., 3 rows × 7 cols) where row index ≥ 3 but < 7 would pass the buggy assert but should fail. Include this test explicitly so the bug is caught and documented. If the bug is fixed before testing, the test still validates correctness.

### 2.5 Row / Column / Submatrix `[debug]`

| Test Name | What It Verifies |
|---|---|
| `RowOutOfBoundsThrows` | `row(n_rows())` triggers assertion |
| `ColOutOfBoundsThrows` | `col(n_cols())` triggers assertion |
| `ParenthesisRowOutOfBoundsThrows` | `(n_rows(), 0)` triggers assertion |
| `ParenthesisColOutOfBoundsThrows` | `(0, n_cols())` triggers assertion |
| `SetRowLengthMismatchThrows` | `set_row(i, wrongLengthVec)` triggers assertion |
| `SetColLengthMismatchThrows` | `set_col(j, wrongLengthVec)` triggers assertion |

### 2.6 Compound Operators `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `PlusEqualsScalar` | `tMat += 2.0` adds 2 to every element |
| `PlusEqualsMatrix` | `tA += tB` adds elementwise |
| `MinusEqualsScalar` | `tMat -= 2.0` subtracts from every element |
| `MinusEqualsMatrix` | `tA -= tB` subtracts elementwise |
| `TimesEqualsScalar` | `tMat *= 3.0` scales all elements |
| `TimesEqualsMatrix` | `tA *= tB` performs matrix multiplication in-place |
| `DivideEqualsScalar` | `tMat /= 2.0` divides all elements |

### 2.7 Print `[semantic]` (smoke test only)

| Test Name | What It Verifies |
|---|---|
| `PrintProducesOutput` | Capture stdout, verify non-empty |

---

## 3. Binary Operators

**File:** `test_LinalgOperators.cpp`
**Approach:** Small deterministic matrices/vectors with hand-computed expected results.

### 3.1 Vector Operators `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `VectorPlusVector` | `{1,2,3} + {4,5,6}` → `{5,7,9}` |
| `VectorPlusScalar` | `{1,2,3} + 10.0` → `{11,12,13}` |
| `ScalarPlusVector` | `10.0 + {1,2,3}` → `{11,12,13}` |
| `VectorMinusVector` | `{5,7,9} - {4,5,6}` → `{1,2,3}` |
| `VectorMinusScalar` | `{5,7,9} - 1.0` → `{4,6,8}` |
| `ScalarMinusVector` | `10.0 - {1,2,3}` → `{9,8,7}` |
| `VectorTimesScalar` | `{1,2,3} * 2.0` → `{2,4,6}` |
| `ScalarTimesVector` | `2.0 * {1,2,3}` → `{2,4,6}` |
| `VectorDivideScalar` | `{2,4,6} / 2.0` → `{1,2,3}` |
| `VectorEqualityTrue` | `{1,2,3} == {1,2,3}` → `true` |
| `VectorEqualityFalse` | `{1,2,3} == {1,2,4}` → `false` |
| `VectorEqualityScalar` | `{5,5,5} == 5.0` → `true` (if both overloads exist) |

### 3.2 Vector Operators `[debug]`

| Test Name | What It Verifies |
|---|---|
| `VectorPlusSizeMismatchThrows` | `Vector(3) + Vector(5)` triggers assertion |
| `VectorMinusSizeMismatchThrows` | `Vector(3) - Vector(5)` triggers assertion |

### 3.3 Matrix Operators `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `MatrixPlusMatrix` | `{{1,2},{3,4}} + {{5,6},{7,8}}` → `{{6,8},{10,12}}` |
| `MatrixMinusMatrix` | Inverse of addition |
| `MatrixTimesMatrix` | 2×2 * 2×2 with known result |
| `MatrixTimesVector` | 2×3 matrix * length-3 vector → length-2 result, hand-computed |
| `MatrixTimesScalar` | `{{1,2},{3,4}} * 2.0` → `{{2,4},{6,8}}` |
| `ScalarTimesMatrix` | `2.0 * {{1,2},{3,4}}` → `{{2,4},{6,8}}` |

### 3.4 Mixed Expression Smoke Test `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `ChainedMatrixVectorProduct` | `A * B * x` matches `A * (B * x)` within tolerance |
| `LinearCombination` | `2.0 * tA + 3.0 * tB` matches entrywise |

---

## 4. Free Functions — Core Math

**File:** `test_LinalgFunctions.cpp`
**Approach:** Small deterministic inputs with known analytical answers.

### 4.1 Dot Product `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `DotProductBasic` | `dot({1,2,3}, {4,5,6})` → `32.0` |
| `DotProductOrthogonal` | `dot({1,0}, {0,1})` → `0.0` |
| `DotProductSelf` | `dot(v, v) == norm(v)^2` |

### 4.2 Cross Product `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CrossProductBasis` | `cross({1,0,0}, {0,1,0})` → `{0,0,1}` |
| `CrossProductAnticommutative` | `cross(a, b) == -cross(b, a)` |
| `CrossProductSelf` | `cross(a, a)` → `{0,0,0}` |

### 4.3 Cross Product `[debug]`

| Test Name | What It Verifies |
|---|---|
| `CrossProductWrongLengthThrows` | `cross(Vector(2), Vector(2))` triggers assertion |

### 4.4 Crossmat `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `Crossmat2DBasic` | 2D normal × 2-row matrix, verify result |
| `Crossmat2DWithScale` | Verify scale factor is applied |
| `Crossmat3DBasic` | 3D normal × 3-row matrix → 3×N result |
| `Crossmat3DWithScale` | Verify accumulation with scale |
| `Crossmat2DDustRemoval` | Near-zero entries (relative to norm) are zeroed |

### 4.5 Determinant `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `DetIdentity` | `det(I)` → `1.0` for 2×2 and 3×3 |
| `DetKnownMatrix` | `det({{1,2},{3,4}})` → `-2.0` |
| `DetSingular` | `det({{1,2},{2,4}})` → `0.0` (within tolerance) |

### 4.6 Inverse `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `InvIdentity` | `inv(I) ≈ I` |
| `InvRoundTrip` | `A * inv(A) ≈ I` for well-conditioned matrix |

### 4.7 Inv2 and Inv3 (BELFEM-authored) `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `Inv2ReturnsCorrectInverse` | 2×2 known matrix: `A * B ≈ I`, return value is `det(A)` |
| `Inv2DetReturnValue` | Verify the returned determinant matches `det(A)` |
| `Inv3ReturnsCorrectInverse` | 3×3 known matrix: `A * B ≈ I`, return value is `det(A)` |
| `Inv3DetReturnValue` | Verify the returned determinant |

### 4.8 Inv2 and Inv3 `[debug]`

| Test Name | What It Verifies |
|---|---|
| `Inv2SingularThrows` | Singular 2×2 triggers assertion (`BELFEM_EPS` threshold) |
| `Inv3SingularThrows` | Singular 3×3 triggers assertion |

### 4.9 Transpose `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `TransposeSwapsDimensions` | 3×4 → 4×3 |
| `TransposeValues` | `trans(A)(j, i) == A(i, j)` |
| `DoubleTransposeIdentity` | `trans(trans(A)) ≈ A` |

### 4.10 Norm `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `NormUnitVector` | `norm({1,0,0})` → `1.0` |
| `NormKnownVector` | `norm({3,4})` → `5.0` |
| `NormZeroVector` | `norm({0,0,0})` → `0.0` |

### 4.11 Sum, Min, Max `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `SumVector` | `sum({1,2,3,4})` → `10.0` |
| `SumMatrix` | `sum({{1,2},{3,4}})` → `10.0` (if sum is total sum) |
| `MinVector` | `min({3,1,4,1,5})` → `1.0` |
| `MaxVector` | `max({3,1,4,1,5})` → `5.0` |
| `MinMatrix` | `min(2×2 matrix)` → smallest element |
| `MaxMatrix` | `max(2×2 matrix)` → largest element |

### 4.12 Linspace `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `LinspaceEndpoints` | `linspace(0, 1, 5)` → first element `0.0`, last `1.0` |
| `LinspaceLength` | Result has exactly N elements |
| `LinspaceUniformSpacing` | Differences between consecutive elements are equal |
| `LinspaceOutputOverload` | 4-argument version (`linspace(start, end, N, tVec)`) fills tVec |

### 4.13 Sort, Unique, Reverse `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `SortAscending` | Produces non-decreasing order |
| `UniqueRemovesDuplicates` | `{3,1,3,2,1}` → sorted unique values |
| `ReverseFlipsOrder` | `{1,2,3}` → `{3,2,1}` |

### 4.14 Append and Combine `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `AppendConcatenates` | `append({1,2}, {3,4})` → `{1,2,3,4}` in first argument |
| `AppendConstOverload` | `append(tA, const tB)` works |
| `Combine2Vectors` | `combine(A, B, C)` → C has `A.length() + B.length()` elements |
| `Combine3Vectors` | 3-input version produces correct concatenation |
| `Combine4Vectors` | 4-input version produces correct concatenation |

---

## 5. Solvers

**File:** `test_LinalgSolvers.cpp`
**Approach:** Small well-conditioned systems with known solutions. Verify via residual norm, not exact `x` comparison. Explicitly test mutation semantics.

### 5.1 gesv (Vector RHS) `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `GesvIdentitySystem` | `I * x = b` → `x == b` |
| `GesvKnown2x2` | Solve `{{2,1},{1,3}} * x = {5,7}`, verify residual `< tTol` |
| `GesvKnown3x3` | Solve a 3×3 system, verify residual |
| `GesvMutatesMatrix` | After `gesv(A, x, P)`, `A` is overwritten (not equal to original) |
| `GesvMutatesRhs` | After `gesv(A, x, P)`, `x` contains solution (not original RHS) |
| `GesvPivotVectorPopulated` | Pivot vector `P` has non-trivial values |

### 5.2 gesv (Matrix RHS) `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `GesvMultipleRhs` | Solve `A * X = B` with 3×2 RHS, verify all columns |

### 5.3 gesv `[debug]`

| Test Name | What It Verifies |
|---|---|
| `GesvRowMismatchThrows` | `A(3×3)` with `x(4)` triggers assertion |
| `GesvColMismatchThrows` | Non-square A or wrong pivot length triggers assertion |

### 5.4 posv `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `PosvSpdSystem` | Solve SPD 3×3 system, verify residual |

---

## 6. Polynomial and Statistical Functions

**File:** `test_LinalgPolynomials.cpp`
**Approach:** Round-trip tests and analytical derivatives.

### 6.1 Polyval `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `PolyvalConstant` | Coeffs `{5.0}` → `polyval == 5.0` for any x |
| `PolyvalLinear` | Coeffs `{2.0, 3.0}` (= 2x+3) at x=4 → `11.0` |
| `PolyvalQuadratic` | Coeffs `{1.0, 0.0, -1.0}` (= x²-1) at x=3 → `8.0` |
| `PolyvalVectorized` | `polyval(coeffs, xVec, yVec)` → yVec matches scalar polyval for each x |

### 6.2 Dpolyval (first derivative) `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `DpolyvalLinear` | Derivative of `{2.0, 3.0}` → `2.0` for any x |
| `DpolyvalQuadratic` | Derivative of `{1.0, 0.0, -1.0}` at x=3 → `6.0` |
| `DpolyvalConstantIsZero` | Single coefficient → `0.0` |

### 6.3 DDpolyval (second derivative) `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `DDpolyvalQuadratic` | Second derivative of `{1.0, 0.0, -1.0}` → `2.0` |
| `DDpolyvalLinearIsZero` | Two coefficients → `0.0` |
| `DDpolyvalConstantIsZero` | One coefficient → `0.0` |

### 6.4 Polyfit `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `PolyfitLinearRoundTrip` | Fit degree-1 to `y = 2x + 3` samples → coefficients ≈ `{2, 3}` |
| `PolyfitQuadraticRoundTrip` | Fit degree-2 to `y = x²` samples → `polyval(coeffs, x) ≈ x²` within tTol |
| `PolyfitReproducesData` | `polyval(polyfit(x, y, n), x)` ≈ `y` for the training points |

### 6.5 Polyfit `[debug]`

| Test Name | What It Verifies |
|---|---|
| `PolyfitLengthMismatchThrows` | `aX.length() != aY.length()` triggers assertion |
| `PolyfitInsufficientSamplesThrows` | `aX.length() <= aN` triggers assertion |

### 6.6 R² (coefficient of determination) `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `R2PerfectFit` | `r2(y, y)` → `1.0` |
| `R2KnownValue` | Approximation with known residual → verify formula |
| `R2ConstantExact` | If all exact values are constant and approximation matches → `1.0` |
| `R2MatrixVersion` | Matrix overload produces valid result |

**NOTE on `r2` matrix version:** The inner loop uses `j < n` (cols) where it should use `j < m` (rows). For square matrices this is masked. Include a non-square test case (e.g., 3×5) to catch this.

### 6.7 Eigen `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `EigenDiagonalMatrix` | `diag(1,2,3)` → eigenvalues `{1,2,3}` (in some order) |
| `EigenSymmetricMatrix` | Known symmetric matrix → real eigenvalues, no NaN |
| `EigenComplexEigenvaluesReturnNaN` | Non-symmetric matrix with complex eigenvalues → some entries are NaN |
| `EigenResultLength` | Output vector has `n_cols()` elements |

### 6.8 Eigen `[debug]`

| Test Name | What It Verifies |
|---|---|
| `EigenNonSquareThrows` | Non-square matrix triggers assertion |

---

## 7. Known Bugs to Test Against

These tests explicitly target suspicious code patterns found during source review. They should be implemented even if the bugs are fixed before testing, to prevent regressions.

| ID | Location | Issue | Test That Catches It |
|---|---|---|---|
| BUG-L1 | `Matrix::submat()` (both backends) | Row indices compared against `n_cols()` instead of `n_rows()`; `aLastRow` used instead of `aLastCol` for column bounds | `SubmatOnNonSquareMatrix`: use 3×7 matrix, access `submat(0, 0, 2, 5)` — should succeed on correct code, but buggy assert compares row 2 < 7 (n_cols) instead of 2 < 3 (n_rows) |
| BUG-L2 | `r2()` matrix version (`fn_r2.hpp`) | Inner loop uses `j < n` (cols) instead of `j < m` (rows) | `R2MatrixVersion` with non-square matrix: exact and approximate are 3×5, bug would iterate wrong range |

---

## 8. Implementation Notes for Claude Code

1. **Floating-point comparisons:** Use `EXPECT_NEAR(actual, expected, tEps)` for exact-in-theory results (integer arithmetic done in reals). Use `EXPECT_NEAR(actual, expected, tTol)` for solver residuals and fitted polynomials. Never use `EXPECT_EQ` on `real` values.
2. **Solver tests:** Verify via residual norm: `EXPECT_NEAR(norm(A_original * x_solution - b_original), 0.0, tTol)`. This requires saving copies of A and b before calling `gesv` because it mutates inputs.
3. **Backend agnosticism:** Do not use `arma::` or `blaze::` types directly in test code. Use only `belfem::Vector<real>`, `belfem::Matrix<real>`, and `belfem::` free functions. The exception is if you need to construct test data via backend APIs — in that case, wrap in a helper and note the dependency.
4. **Do not test backend correctness:** If `dot({1,0,0}, {0,1,0})` returns `0.0`, that's sufficient. Do not verify BLAS internal behavior or test GEMM with random large matrices.
5. **Expression template results:** When a free operator returns a backend expression type rather than a `Vector<T>` or `Matrix<T>`, assign the result to a `Vector<T>` or `Matrix<T>` before comparing. This forces evaluation and avoids type-deduction issues.
6. **Print tests:** Smoke test only — capture stdout, check non-empty. No golden-file comparison.
7. **Prioritize tests that catch the known bugs** (BUG-L1, BUG-L2) in the first implementation pass.

---

## 9. Codex Audit Checklist

When reviewing Claude Code's test implementation, verify:

- [ ] Every row in the test matrices above has a corresponding `TEST` or `TYPED_TEST`
- [ ] Debug tests wrapped in `#ifndef NDEBUG`
- [ ] BELFEM naming conventions (`t` prefix for locals, BELFEM types)
- [ ] `EXPECT_NEAR` used for all floating-point comparisons, never `EXPECT_EQ`
- [ ] No backend-specific types (`arma::`, `blaze::`) in test assertions
- [ ] Solver tests save copies of inputs before calling mutating functions
- [ ] Solver tests verify via residual norm, not exact solution comparison
- [ ] `inv2`/`inv3` tests verify both the inverse matrix (A*B≈I) and the returned determinant
- [ ] BUG-L1 and BUG-L2 regression tests are present with non-square matrices
- [ ] Expression results assigned to concrete BELFEM types before comparison
