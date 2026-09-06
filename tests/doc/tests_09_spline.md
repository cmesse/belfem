# Spline Module Test Suite

**Date:** 2026-03-24
**Test file:** `tests/math/test_Spline.cpp`
---

## Test Count by Suite

| Suite | Tests | Coverage |
|-------|-------|----------|
| `TrussPoly` | 5 | Value interpolation, derivative interpolation, cubic exactness, linear function, non-origin interval |
| `GluePoly` | 1 | All 5 conditions (value, 1st/2nd deriv at both endpoints) |
| `HelpMatrix` | 4 | Size, tridiagonal interior, NoCurvature BC, Tangent BC |
| `HelpMatrixDebug` | 1 | Too few points throws |
| `Spline` | 9 | Constant, linear, cubic clamped, knot interpolation, continuity, natural BC, extrapolation (below/above), accessors |
| `SplineIntegral` | 2 | Linear integral, continuity across knots |
| `SplineIntegralDebug` | 1 | Integrate without create throws |
| `SplineEntropy` | 2 | Entropy construction, dentropy finite-difference consistency |
| `SplineEntropyDebug` | 2 | Entropy/dentropy without create throws |
| `SplineUpdate` | 2 | Coefficient change, grid preservation |
| `SplineVariant` | 1 | Empty container constructor |
| **Total** | **30** | |

---

## Three-Layer Compilation Structure

1. **Always runs:** TrussPoly, GluePoly, HelpMatrix, EmptyContainerConstructor
2. **Debug only (`#ifndef NDEBUG`):** HelpMatrixDebug, SplineIntegralDebug, SplineEntropyDebug
3. **SuperLU required (`#ifdef BELFEM_SUPERLU`):** All Spline, SplineIntegral, SplineEntropy, SplineUpdate

   Corrected 2026-08-10. This layer was gated on `BELFEM_SUITESPARSE`, which compiled every
   spline construction/evaluation test out of the default build — SuiteSparse is off by
   default (legacy, not BSD-3 clean), so spline coverage in the fast gate was zero. The gate
   was wrong on its own terms too: `Spline::update_data` solves with **SuperLU**
   (`cl_Spline.cpp:427-429`), which is ON by default, so these tests were gated on a
   dependency the code under test does not use.

---

## Key Design Decisions

### Constructor Disambiguation
Two Spline constructors are ambiguous with 3 args. All tests pass `0.0` as 4th argument: `Spline(tX, tY, tA, 0.0)`.

### Extrapolation Behavior
`find_col()` clamps interval selection, but `eval()` uses unclamped x. Result: polynomial extrapolation. Tests verify via `eval(x, col)` two-argument overload.

### Natural BCs Don't Preserve Quadratics
Approximate tolerances used for non-linear functions with natural BCs.

---

## API Gotchas

- `create_entropy` is **private** — activate via constructor with `aXref > 0`
- `create_integral` is **public** — call after construction
- `create_helpmatrix` takes `const real&` for size, not `index_t`
- Equidistant spacing required (asserted with 1e-9 relative tolerance)
- `entropy(x)` requires `x > 0` (uses `log(x)`)
- `update_data()` only runs on rank 0 — caller must synchronize

---

## Known Issues (Deferred)

- save/load HDF5 round-trip broken (field name mismatch)
- `find_col()` bounds check has `#HACK` comment
- `update_data()` lacks MPI synchronization

---

## Future Work

- [ ] Fix save/load field name mismatch
- [ ] Parabolic and mixed BC tests
- [ ] Input validation debug tests
- [ ] GluePolyQuadraticInput test
