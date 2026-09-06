# Linear Algebra Module Test Suite

**Date:** 2026-03-23

---

## Test Count by File

| File | Suites | Tests |
|------|--------|-------|
| `test_Vector.cpp` | `Vector`, `VectorDebug` | 76 |
| `test_Matrix.cpp` | `Matrix`, `MatrixDebug` | 71 |
| `test_LinalgSolvers.cpp` | `LinalgSolver`, `LinalgSolverDebug`, `Crossmat`, `Polynomial`, `PolynomialDebug`, `R2`, `Eigen`, `EigenDebug` | 38 |
| **Total** | **12 suites** | **185** |

---

## Key Design Decisions

### Floating-Point Comparison
- `tEps = 1e-12` for exact-in-theory results (integer arithmetic via reals)
- `tTol = 1e-9` for solver residuals, eigenvalues, polyfit
- Never use `EXPECT_EQ` on `belfem::real` values

### Expression Templates
- Armadillo returns expression templates from operator overloads, not `Vector<T>`.
- Use `Vector<real> tResult( expr )` (direct initialization), NOT `Vector<real> tResult = expr` (copy initialization fails).

### Backend Agnosticism
- `VectorDataExposesBackend` uses `data()` pointer (works for both Armadillo and Blaze), NOT `vector_data()` with backend-specific `operator()`.
- No `#ifdef BELFEM_ARMADILLO` in test code.

---

## @warning API Gotchas

- `gesv()` and `posv()` mutate all their arguments. Save copies before calling.
- `polyval` uses descending power order: `{a_n, ..., a_1, a_0}`.
- `r2()` returns 1.0 when SSres or SStot is below `BELFEM_EPSILON`.
- `eigen()` aborts by default on a complex eigenvalue; pass `aAbortOnComplex = false` to get
  `BELFEM_QUIET_NAN` entries and a count of them instead. The threshold is `BELFEM_EPSILON`
  on both backends (it was a hard-coded `1e-15` on Blaze until 2026-08-18).
- `eigen_sym()` is the symmetric entry point: real eigenvalues by construction, ascending,
  reading the upper triangle. Symmetry is not checked.
- `cross()` BELFEM wrapper has `BELFEM_ASSERT(length == 3)`.

---

## @todo Future Work

- [ ] Blaze backend compilation verification
- [ ] Additional `gesv` matrix-RHS tests
- [ ] ARPACK eigenvalue tests (requires ARPACK linkage)
