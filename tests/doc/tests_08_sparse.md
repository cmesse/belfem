# Sparse Module Test Suite

**Date:** 2026-03-23

---

## Test Count by Suite

| Suite | Tests | Coverage |
|-------|-------|----------|
| `SpMatrixConstruction` | 7 | Dense CSR/CSC, nonzero count, diagonal, full, rectangular |
| `SpMatrixRoundTrip` | 3 | CSR/CSC readback, full dense |
| `SpMatrixStructure` | 4 | Pointer/index arrays, sizes |
| `SpMatrixAccess` | 4 | Element read, structural zero, const |
| `SpMatrixDebug` | 3 | Writable zero throws, data OOB |
| `SpMatrixIndexing` | 4 | Cpp/Fortran base, Fortran round-trip, double conversion |
| `SpMatrixMultiply` | 6 | CSR/CSC, rectangular, transposed, operator*, indexing restore |
| `SpMatrixTranspose` | 5 | Type swap, dimension swap, value preservation, non-symmetric rectangular, double transpose |
| `SpMatrixCOO` | 3 | CSR/CSC COO creation, free |
| `SolverConstruction` | 3 | UMFPACK creation, facade type |
| `SolverSolve` | 2 | Known-solution, re-solve |
| `SolverLifecycle` | 2 | Lazy init, free/re-init |
| `SolverParameters` | 8 | All parameter accessors |
| `SolverUnlinked` | 4 | Error paths for unlinked backends |
| **Total** | **56** | |

---

## @warning API Gotchas

- **Compilation flags must match library:** The test binary MUST be compiled with the same defines as the library (especially `-DBELFEM_NETLIB`). Mismatched flags cause struct layout differences that produce silent wrong results, not compilation errors.
- `operator()` writable: uses `BELFEM_ASSERT` — asserts if position is structural zero (debug only).
- `operator()` const: returns `mZero` (0.0) for structural zeros — always safe.
- `operator()` now asserts `indexing_base() == 0` — catches Fortran-mode access attempts.
- `multiply()` restores indexing base after call (transparent to caller).
- `transpose()` is void and in-place. Do NOT write `SpMatrix tMT = tM.transpose()`.
- `indexing_base()` returns `int` (0 or 1), NOT enum.
- `index()` returns `number_of_nonzeros()` for structural zeros.
- Node/element IDs in graph constructor start at 1.

---

## @todo Future Work

- [ ] Raw-array constructor + `sort_entries()` tests
- [ ] MUMPS job-5/job-6 re-solve distinction
- [ ] Graph/permutation utility tests
- [ ] PETSc/STRUMPACK MPI tests (separate file with custom main)
- [ ] ARPACK error mapping tests
