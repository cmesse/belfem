# BELFEM Spline Module Tests — Detailed Plan

**Date:** 2026-03-22
**Purpose:** Method-level test matrix for the spline module (`src/numerics/spline/`)
**Depends on:** `tests_0_strategy.md`, `tests_2_linalg.md` (Vector/Matrix), `tests_8_sparse.md` (SpMatrix, Solver/UMFPACK)
**Confidence:** High on interpolation and polynomial helpers. Medium on entropy/integral modes. Low on I/O (known bug).

---

## Module Overview

| Component | What It Does |
|---|---|
| `Spline` class | Cubic spline on **equidistant** knots. Stores per-interval polynomial coefficients as `Matrix<real>(4, N)`. Solves tridiagonal system via UMFPACK. Supports 3 boundary conditions, optional entropy/integral augmentation. |
| `create_truss_poly()` | Cubic Hermite polynomial from 2 endpoints with prescribed values and derivatives. Core per-interval building block. |
| `create_glue_poly()` | Quartic polynomial connecting 3 points (value+derivative at endpoints, value at midpoint). |
| `spline::create_helpmatrix()` | Assembles tridiagonal SpMatrix for the derivative system. Boundary-condition dependent. |

**Critical constraint:** X-data must be equidistant (asserted with `1e-9` relative tolerance). Non-uniform knots are not supported.

---

## Known Regression Targets

| ID | Location | Issue | Regression Test |
|---|---|---|---|
| BUG-S1 | `save()` / `load()` | Field name mismatch: `save()` writes `"min"`, `"n"`, `"data"` but `load()` reads `"xmin"`, `"npoints"`, `"coeffs"`. A spline saved with `save()` cannot be loaded with `load()`. | Defer to I/O integration phase — document as known issue. |
| NOTE-S1 | `find_col()` line 417 | Bounds-check assertion is commented out (`#HACK`). Out-of-range x selects the boundary interval but **evaluates the polynomial at the original (out-of-range) x value** — this is polynomial extrapolation, not clamping. The x is NOT clamped before polynomial evaluation. | Test extrapolation behavior explicitly: verify that `eval(xmin - h)` uses the first interval's polynomial at `xmin - h`, not at `xmin`. |
| NOTE-S2 | `update_data()` line 941 | Only executes on rank 0. Does NOT call `synchronize()` afterward. Non-root ranks are not updated. | Document as "caller must synchronize" or test in MPI context. |

---

## Floating-Point Comparison

```cpp
namespace
{
    const belfem::real tEps = 1e-12 ;  // for exact polynomial recovery
    const belfem::real tTol = 1e-9 ;   // for approximate interpolation
}
```

---

## Test File Structure

```
tests/numerics/
├── test_Spline.cpp              # All spline tests in one file
```

---

## 1. Helper Polynomials

### 1.1 create_truss_poly `[semantic]`

The polynomial `f(x) = a·x³ + b·x² + c·x + d` must satisfy four conditions: `f(x₀) = f₀`, `f'(x₀) = f'₀`, `f(x₁) = f₁`, `f'(x₁) = f'₁`.

| Test Name | What It Verifies |
|---|---|
| `TrussPolyInterpolatesValues` | Evaluate at x₀ and x₁ → matches f₀ and f₁ |
| `TrussPolyInterpolatesDerivatives` | Evaluate derivative at x₀ and x₁ → matches f'₀ and f'₁ |
| `TrussPolyLinearFunction` | Input from f(x) = 3x+1 → cubic coefficient ≈ 0, correct linear result |
| `TrussPolyCubicExact` | Input from f(x) = x³ at two points → recovered coefficients match {1, 0, 0, 0} |
| `TrussPolyNonOriginInterval` | Interval [2, 5] → still correct interpolation |

### 1.2 create_glue_poly `[semantic]`

The quartic `f(x) = a·x⁴ + b·x³ + c·x² + d·x + e` must satisfy: `f(X-dX) = f₀`, `f'(X-dX) = df₀`, `f(X) = f₁`, `f(X+dX) = f₂`, `f'(X+dX) = df₂`.

| Test Name | What It Verifies |
|---|---|
| `GluePolyInterpolatesAllFiveConditions` | Evaluate at X-dX, X, X+dX and derivatives at endpoints → all match |
| `GluePolyQuadraticInput` | Input from f(x) = x² → quartic coefficients have a≈0 |

---

## 2. Help Matrix

### 2.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `HelpMatrixSizeCorrect` | `create_helpmatrix(N, dX, A)` → A is N×N square |
| `HelpMatrixTridiagonal` | Inner rows have the pattern: off-diagonal = 1/dX, diagonal = 4/dX |
| `HelpMatrixNoCurvatureBC` | First row: diagonal = 2/dX, off-diag = 1/dX |
| `HelpMatrixParabolicBC` | First row: diagonal = 1/dX, off-diag = 1/dX |
| `HelpMatrixTangentBC` | First row: diagonal = 1, row(0,1) = 0 |

### 2.2 Tests `[debug]`

| Test Name | What It Verifies |
|---|---|
| `HelpMatrixTooFewPointsThrows` | `create_helpmatrix(3, ...)` → assertion ("Need at least four datapoints") |

---

## 3. Spline Construction and Evaluation

### 3.1 Polynomial Exactness `[semantic]`

A cubic spline through data sampled from a polynomial of degree ≤ 3 should reproduce it exactly.

| Test Name | Function | BC | What It Verifies |
|---|---|---|---|
| `SplineConstantFunction` | f(x) = 7 | NoCurvature | `eval(x) ≈ 7` for all x in domain |
| `SplineLinearFunction` | f(x) = 3x + 1 | NoCurvature | `eval(x_mid) ≈ f(x_mid)`, `deval ≈ 3` |
| `SplineQuadraticNatural` | f(x) = x² | NoCurvature | `eval` and `deval` match at knots and midpoints (approximate — natural BC doesn't preserve quadratics exactly) |
| `SplineCubicClamped` | f(x) = x³ | Tangent (f'=3x² at endpoints) | `eval(x) ≈ x³`, `deval(x) ≈ 3x²`, `ddeval(x) ≈ 6x` — exact recovery |
| `SplineQuadraticClamped` | f(x) = 2x²-3x+5 | Tangent (f'=4x-3) | Exact recovery of value, derivative, second derivative |

### 3.2 Knot Interpolation `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `SplineInterpolatesKnots` | For every data point: `eval(x_i) ≈ y_i` |
| `SplineContinuousAcrossKnots` | `eval(x_i - ε) ≈ eval(x_i + ε)` at interior knots |
| `SplineDerivativeContinuous` | `deval(x_i - ε) ≈ deval(x_i + ε)` at interior knots |

### 3.3 Boundary Conditions `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `NoCurvatureBC` | `ddeval(x_min) ≈ 0` and `ddeval(x_max) ≈ 0` (natural spline) |
| `TangentBC` | `deval(x_min) ≈ prescribed_slope_0` and `deval(x_max) ≈ prescribed_slope_1` |
| `ParabolicBC` | First interval cubic coefficient ≈ 0 (quadratic runout) |
| `MixedBCStartTangentEndNatural` | `deval(x_min) ≈ prescribed`, `ddeval(x_max) ≈ 0` |

### 3.4 Extrapolation (NOTE-S1) `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `EvalBelowRangeUsesFirstInterval` | `eval(x_min - 1)` returns first interval polynomial evaluated at x_min (clamped) |
| `EvalAboveRangeUsesLastInterval` | `eval(x_max + 1)` returns last interval polynomial evaluated at x_max (clamped) |
| `DevalBelowRange` | `deval(x_min - 1) ≈ deval(x_min)` |

### 3.5 Accessors `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `SplineNReturnsPointCount` | `n()` matches input vector length |
| `SplineXminXmax` | `x_min()` and `x_max()` match input range |
| `SplineDeltaX` | `delta_x()` matches expected step size |
| `SplineCoefficientsMatrix` | `coefficients()` returns 4×N matrix |

---

## 4. Integral Mode

### 4.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `IntegralModeFlagSet` | After `create_integral()`, integration functions don't assert |
| `IntegralConsistentWithEval` | `integrate(a, b) ≈ numerical_quadrature(eval, a, b)` within tolerance |
| `IntegralOfLinear` | For f(x) = 2x+1: `integrate(0, L) ≈ L² + L` |
| `IntegralReferenceOffset` | `create_integral(xref, yref)` → `integrate(xref) ≈ yref` |
| `IntegralContinuousAcrossKnots` | `integrate(x_i - ε) ≈ integrate(x_i + ε)` |

### 4.2 Tests `[debug]`

| Test Name | What It Verifies |
|---|---|
| `IntegrateWithoutCreateThrows` | `integrate(x)` without `create_integral()` → assertion |

---

## 5. Entropy Mode

### 5.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `EntropyModeFlagSet` | After construction with `aXref > 0`, entropy functions don't assert |
| `EntropyContinuousAcrossKnots` | `entropy(x_i - ε) ≈ entropy(x_i + ε)` |
| `EntropyReferenceOffset` | `entropy(xref) ≈ sref` |
| `DentropyConsistentWithFiniteDifference` | `dentropy(x) ≈ (entropy(x+h) - entropy(x-h)) / (2h)` |

### 5.2 Tests `[debug]`

| Test Name | What It Verifies |
|---|---|
| `EntropyWithoutCreateThrows` | `entropy(x)` on a plain spline → assertion |
| `DentropyWithoutCreateThrows` | `dentropy(x)` on a plain spline → assertion |

---

## 6. update_data

### 6.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `UpdateDataChangesCoefficients` | After `update_data()` with new y-values, `eval` matches new function |
| `UpdateDataPreservesGrid` | `x_min()`, `x_max()`, `delta_x()` unchanged after update |

---

## 7. Construction Variants

### 7.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `EmptyContainerConstructor` | `Spline(N, xmin, xmax)` → `n() == N`, `x_min() == xmin`, `x_max() == xmax` |

### 7.2 Tests `[debug]`

| Test Name | What It Verifies |
|---|---|
| `CheckInputNonSortedXThrows` | Unsorted x-vector → assertion |
| `CheckInputNonEquidistantXThrows` | Non-uniform spacing → assertion |
| `CheckInputLengthMismatchThrows` | `x.length() != y.length()` → assertion |
| `CheckInputTooFewPointsThrows` | 3 or fewer points → assertion |
| `CheckInputMatrixSizeMismatchThrows` | Help matrix size ≠ N → assertion |

---

## 8. What We Do NOT Test (Deferred)

- HDF5 save/load (BUG-S1: known field-name mismatch)
- `save_to_database` (separate database format)
- MPI synchronization (NOTE-S2: requires multi-process)
- Parallel constructors (`Spline(masterProc)`)

---

## 9. Implementation Notes for Claude Code

1. **All x-data must be equidistant.** Use `linspace`-style generation: `x(k) = xmin + k * deltaX`. Do not use arbitrary spacing.
2. **Spline requires UMFPACK.** The constructor calls `Solver(SolverType::UMFPACK)`. Tests will fail if UMFPACK is not linked. Wrap in `#ifdef BELFEM_SUITESPARSE` if needed.
3. **`create_helpmatrix` takes the number of points as `aSize`, not the number of intervals.** The first argument is `aSize` (= N points), the second is `aDeltaX`.
4. **The truss-poly evaluation is Horner form.** `eval(x)` computes `((a*x + b)*x + c)*x + d`. The x-values are raw (not shifted to interval-local coordinates). Coefficients are in global x-coordinates.
5. **Extrapolation is polynomial extension, not constant.** `find_col()` clamps x to `[xmin, xmax]` and evaluates the boundary polynomial at the clamped x. This gives the appearance of clamping but the actual returned value depends on the boundary polynomial evaluated at `xmin` or `xmax`.
6. **`coefficients()` returns a 4×N matrix** where N = number_of_points. Row 0 = cubic coefficient, row 3 = constant. Column k = interval k. The last column duplicates the second-to-last.
7. **For polynomial exactness tests:** use enough points (≥ 10) to avoid boundary effects dominating interior accuracy. The clamped BC (Tangent) gives exact recovery for cubics; natural BC does not.
8. **BUG-S1 is confirmed.** Do not write a save/load round-trip test that uses both `save()` and `load()` — it will fail due to the field-name mismatch.

---

## 10. Codex Audit Checklist

When reviewing Claude Code's test implementation, verify:

- [ ] All x-data is equidistant (computed, not hardcoded with potential rounding)
- [ ] Polynomial exactness tests use clamped BC for cubic recovery
- [ ] `create_truss_poly` tests verify all 4 interpolation conditions
- [ ] `create_glue_poly` tests verify all 5 interpolation conditions
- [ ] Boundary condition tests cover all 3 types
- [ ] Extrapolation behavior tested and documented
- [ ] Integral and entropy modes tested with continuity checks
- [ ] Debug-only tests for all `BELFEM_ASSERT` guards
- [ ] `#ifdef BELFEM_SUITESPARSE` guard if UMFPACK might not be linked
- [ ] `EXPECT_NEAR` for all floating-point comparisons
- [ ] BELFEM naming conventions (`t` prefix for locals)
- [ ] BUG-S1 documented, not tested with broken round-trip