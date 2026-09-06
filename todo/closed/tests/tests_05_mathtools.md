# BELFEM Math Tools Tests — Detailed Plan

**Date:** 2026-03-22
**Purpose:** Method-level test matrix for the math tools module
**Depends on:** `tests_0_strategy.md` (conventions), `tests_2_linalg.md` (Vector/Matrix/gesv/polyval assumed tested)
**Confidence:** High on API surface. High on identified regression targets. Medium on `circle_from_points` divisor formula (flagged for verification).

---

## Module Overview

A collection of standalone mathematical utility functions — no class, no shared state. Each header is self-contained.

| Function | Header | Category | BELFEM Logic |
|---|---|---|---|
| `sign()` | `fn_sign.hpp` | Scalar | Trivial |
| `cardano()` | `fn_cardano.hpp` | Root-finding | Complex, multi-branch |
| `circle_from_points()` | `fn_circle_from_points.hpp` | Geometry | Newton iteration |
| `create_beam_poly()` | `fn_create_beam_poly.hpp` | Interpolation | C¹ cubic Hermite via gesv |
| `create_fifth_order_beam_poly()` | `fn_create_fifth_order_beam_poly.hpp` | Interpolation | C² quintic Hermite via gesv |
| `cubic_bezier()` | `fn_cubic_bezier.hpp` | Curve evaluation | BELFEM ξ∈[-1,1] parametrization |
| `cubic_bezier_derivative()` | `fn_cubic_bezier.hpp` | Curve evaluation | Derivative of above |
| `find_interval()` | `fn_find_interval.hpp` | Search | Binary search |
| `quadratic_gradient()` | `fn_quadratic_gradient.hpp` | Differentiation | 3-point finite difference |
| `rotation_matrix()` (axis-angle) | `fn_rotation_matrix.hpp` | Rotation | Rodrigues formula, backend-specific |
| `rotation_matrix()` (Euler) | `fn_rotation_matrix.hpp` | Rotation | DIN 9300 yaw-pitch-roll, backend-specific |
| `rotation_matrix_strip()` | `fn_rotation_matrix.hpp` | Rotation | First 2 rows of axis-angle matrix |
| `symratiospace()` | `fn_symrationspace.hpp` | Mesh spacing | Symmetric geometric grading |

---

## Floating-Point Comparison

```cpp
namespace
{
    const belfem::real tEps = 1e-12 ;  // for exact-in-theory results
    const belfem::real tTol = 1e-9 ;   // for iterative solver results (circle_from_points)
}
```

For `cardano`: prefer **root substitution** over comparing to approximate values. If cardano returns root x, verify `|a·x³ + b·x² + c·x + d| < tEps` rather than `|x - x_expected| < tEps`. This is more robust and catches sign errors.

---

## Known Regression Targets

These are suspicious code patterns found during source review. Tests must exercise them with inputs that would expose the bug if present.

| ID | Location | Issue | Regression Test |
|---|---|---|---|
| BUG-M1 | `symratiospace()` line 930 | Midpoint set to `0.5 * (aXmax - aXmin)` instead of `aXmin + 0.5 * (aXmax - aXmin)`. For `aXmin = 0` this is masked. | Test with `aXmin = 2.0, aXmax = 10.0` — midpoint should be 6.0, not 4.0 |
| BUG-M2 | `find_interval()` interior loop | `aIndex` and `aXi` are assigned inside the loop body, but the loop can `break` at `k - i == 1` before that assignment executes. Output may be stale from a previous iteration. | Test with value exactly in the middle of a 5-element sorted array; verify `aIndex` and `aXi` are correct |
| BUG-M3 | `quadratic_gradient()` | Length-match assertion is commented out. Functions access `aF(0..2)` or `aF(n-3..n-1)` without checking `n >= 3`. | Document precondition in tests; test with exactly 3 points to verify boundary formulas |
| BUG-M4 | `rotation_matrix()` Euler-angle, Blaze path | Writes `tData[10]` for a 3×3 matrix, skipping indices 3 and 7 (assumes Blaze padding). If padding assumptions are wrong, writes to wrong locations. | Test via `operator()` access rather than `data()` pointer; verify R * Rᵀ ≈ I |

---

## Test File Structure

```
tests/math/
├── test_MathTools.cpp           # All mathtools tests in one file
```

The module is compact enough for a single file with well-named test groups.

---

## 1. Sign Function

### 1.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `SignPositive` | `sign(5.0)` → `1` |
| `SignNegative` | `sign(-3.2)` → `-1` |
| `SignZero` | `sign(0.0)` → `0` |
| `SignNegativeZero` | `sign(-0.0)` → `0` (IEEE 754 negative zero) |
| `SignInfinity` | `sign(+inf)` → `1`, `sign(-inf)` → `-1` |
| `SignInteger` | `sign(7)` → `1`, `sign(-7)` → `-1`, `sign(0)` → `0` (template works with int) |

---

## 2. Cardano (Cubic Equation Solver)

### 2.1 Root Substitution Helper

Claude Code should implement a test-local helper that verifies a root by substitution:

```cpp
void verify_root( const Vector<real> & aA, real aX )
{
    real tVal = aA(0)*aX*aX*aX + aA(1)*aX*aX + aA(2)*aX + aA(3);
    EXPECT_NEAR( tVal, 0.0, tEps );
}
```

### 2.2 Degenerate Cases (a=0 fallback) `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CardanoLinearFallback` | `{0, 0, 5, 7}` → 1 root, x = -7/5. Verify by substitution |
| `CardanoConstantNoRoots` | `{0, 0, 0, 1}` → 0 roots (empty vector) |
| `CardanoQuadraticTwoRoots` | `{0, 1, 3, 2}` → 2 roots (-2, -1). Sorted. Verify by substitution |
| `CardanoQuadraticNoRealRoots` | `{0, 1, 0, 1}` → discriminant < 0, 0 roots |
| `CardanoZeroAllCoeffs` | `{0, 0, 0, 0}` → 0 roots (c coefficient also zero) |

### 2.3 Cubic One Real Root (D > 0) `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CardanoOneRealRoot` | `{1, 0, 1, 1}` → 1 root. Verify by substitution |
| `CardanoOneRealRootKnown` | `{1, -6, 12, -8}` = (x-2)³ → at least 1 root near 2.0. Verify by substitution |

### 2.4 Cubic Three Real Roots (D < 0) `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CardanoThreeRealRoots` | `{1, -6, 11, -6}` = (x-1)(x-2)(x-3). Roots sorted: 1, 2, 3. Verify each by substitution |
| `CardanoThreeRealRootsIrreducible` | `{1, 0, -3, -1}` → 3 roots. All verified by substitution. Sorted. |
| `CardanoThreeRealRootsNegativeLeading` | `-1·x³ + ...` — verify sign of leading coefficient doesn't break the solver |

### 2.5 Cubic D = 0 (Repeated Root) `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CardanoTripleRoot` | `{1, -3, 3, -1}` = (x-1)³. Returns 1 root at 1.0 |
| `CardanoDoubleRoot` | `{1, 0, -3, 2}` = (x-1)²(x+2). Verify returned roots by substitution |

### 2.6 Numerical Robustness `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CardanoSmallLeadingCoeff` | `{1e-14, 1, 0, -1}` → effectively quadratic. Verify result makes sense |
| `CardanoLargeCoeffs` | Coefficients in 1e6 range — verify roots by substitution |
| `CardanoReturnedRootsSorted` | For 3-root case, verify `aX(0) <= aX(1) <= aX(2)` |

---

## 3. Circle from Points

### 3.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CircleFromPointsUnitCircle` | Points on unit circle: (1,0), (0,1), (-1,0) → center (0,0), radius 1.0 |
| `CircleFromPointsTranslatedCircle` | Points on circle with center (5, 3), radius 2 → verify center and radius |
| `CircleFromPointsLargeRadius` | Points on circle with radius 1000 → verify convergence |
| `CircleFromPointsAllPointsEquidistant` | Verify all 3 input points are equidistant from returned center (within `tTol`) |

### 3.2 Tests `[debug]`

| Test Name | What It Verifies |
|---|---|
| `CircleFromPointsWrongXLengthThrows` | `aX.length() != 3` → assertion |
| `CircleFromPointsWrongYLengthThrows` | `aY.length() != 3` → assertion |

### 3.3 Notes

`circle_from_points` uses `BELFEM_ERROR` (always active) for the 1000-iteration guard. Collinear points would cause `tDiv ≈ 0` in the initial guess, producing garbage that the Newton solver may not converge from. A collinear-point test would verify the `BELFEM_ERROR` fires. However, this is hard to trigger reliably — skip unless easy to set up.

---

## 4. Beam Polynomial Interpolation

### 4.1 Cubic Beam Poly `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `BeamPolyInterpolatesValues` | Create poly at x1=0, x2=1 with known f and f'. Evaluate poly at x1 and x2 via `polyval` → matches f1 and f2 |
| `BeamPolyInterpolatesDerivatives` | Evaluate derivative (via `dpolyval`) at x1 and x2 → matches df1 and df2 |
| `BeamPolyKnownCubic` | Input is a known cubic f(x) = x³ at two points → recovered coefficients match |
| `BeamPolyLinearFunction` | f(x) = 2x+3, f'(x) = 2 → poly should reduce to linear |

### 4.2 Fifth-Order Beam Poly `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `FifthOrderInterpolatesValues` | Evaluate at endpoints → matches prescribed f |
| `FifthOrderInterpolatesFirstDerivatives` | Evaluate 1st derivative at endpoints → matches prescribed f' |
| `FifthOrderInterpolatesSecondDerivatives` | Evaluate 2nd derivative at endpoints → matches prescribed f'' |
| `FifthOrderKnownQuintic` | Input is f(x) = x⁵ at two points → recovered coefficients match |

---

## 5. Cubic Bézier

### 5.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `BezierEndpointMinus1` | At ξ=-1, result passes through first control point column |
| `BezierEndpointPlus1` | At ξ=+1, result passes through last control point column |
| `BezierMidpoint` | At ξ=0, result is a known weighted combination of control points |
| `BezierStraightLine` | Collinear control points → Bézier evaluates on the line for any ξ |
| `BezierDerivativeConsistency` | `cubic_bezier_derivative` at ξ matches finite difference `(bezier(ξ+h) - bezier(ξ-h)) / 2h` within tolerance proportional to h² |
| `BezierDustCleanup` | Near-zero values (< BELFEM_EPSILON) are cleaned to exactly 0.0 |

### 5.2 Tests `[debug]`

| Test Name | What It Verifies |
|---|---|
| `BezierWrongPointLengthThrows` | `aPoint.length() != aPoints.n_rows()` → assertion |
| `BezierWrongColCountThrows` | `aPoints.n_cols() != 4` → assertion |
| `BezierWrongWorkLengthThrows` | `aWork.length() != 4` → assertion |

---

## 6. Find Interval (Binary Search)

### 6.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `FindIntervalBelowRange` | Value ≤ second element → `aIndex == 0`, `aXi` computed correctly |
| `FindIntervalAboveRange` | Value ≥ second-to-last element → `aIndex == N-1`, `aXi` correct |
| `FindIntervalExactKnot` | Value exactly on a data point → returns correct interval and `aXi = 0.0` or `1.0` |
| `FindIntervalInterior` | Value between interior knots → `aIndex` identifies correct interval, `aXi ∈ [0,1]` |
| `FindIntervalAllIntervals` | Sweep a value through every pair of adjacent knots in a 10-element array → verify `aIndex` and `aXi` for each |
| `FindIntervalRegressionBugM2` | 5-element sorted array, value at exact midpoint — verify `aIndex` and `aXi` are assigned correctly (BUG-M2 target) |

---

## 7. Quadratic Gradient

### 7.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `QuadraticGradientExactQuadratic` | f(x) = 3x² + 2x + 1 sampled at 5 points → `quadratic_gradient` recovers exact derivative `6x + 2` at interior points |
| `QuadraticGradientLeftBoundary` | Derivative at index 0 matches analytical value for quadratic |
| `QuadraticGradientRightBoundary` | Derivative at last index matches analytical value for quadratic |
| `QuadraticGradientLinearFunction` | f(x) = 5x + 3 → gradient is 5.0 everywhere |
| `QuadraticGradientMinimumThreePoints` | Test with exactly 3 points — all three branches exercised (index 0, 1, 2) |

### 7.2 Notes

The length-match assertion (`aX.length() == aF.length()`) is commented out. Tests should use matching-length vectors and document that the function requires `n >= 3` and matching lengths even though it doesn't assert this.

---

## 8. Rotation Matrices

### 8.1 Axis-Angle `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `RotationMatrixIdentity` | Angle = 0 → R ≈ I |
| `RotationMatrix90AboutZ` | Axis {0,0,1}, angle π/2 → R * {1,0,0} ≈ {0,1,0} |
| `RotationMatrix180AboutX` | Axis {1,0,0}, angle π → R * {0,1,0} ≈ {0,-1,0} |
| `RotationMatrixOrthogonal` | R * Rᵀ ≈ I for arbitrary axis and angle |
| `RotationMatrixDetPlusOne` | det(R) ≈ +1 |
| `RotationMatrixConsistentWithQuaternion` | Same axis and angle → quaternion module's `quaternion_to_rotation_matrix` matches (cross-module test) |

### 8.2 Euler Angles (DIN 9300) `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `EulerIdentity` | Yaw=pitch=roll=0 → R ≈ I |
| `EulerPureYaw` | Yaw = π/2, rest 0 → known matrix entries |
| `EulerPurePitch` | Pitch = π/2, rest 0 → known matrix entries |
| `EulerPureRoll` | Roll = π/2, rest 0 → known matrix entries |
| `EulerOrthogonal` | R * Rᵀ ≈ I for arbitrary angles |
| `EulerDetPlusOne` | det(R) ≈ +1 |
| `EulerRegressionBugM4` | Arbitrary non-trivial angles → verify via `operator()` that matrix entries are correct, not by inspecting `data()` pointer (catches Blaze padding bug if present) |

### 8.3 Strip Version `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `StripMatchesFullMatrixFirstTwoRows` | `rotation_matrix_strip` output matches first 2 rows of `rotation_matrix` for same axis/angle |

### 8.4 Rotation Matrices `[debug]`

| Test Name | What It Verifies |
|---|---|
| `AxisAngleWrongAxisLengthThrows` | `aAxis.length() != 3` → assertion |
| `AxisAngleWrongMatrixSizeThrows` | Matrix not 3×3 → assertion |
| `EulerWrongMatrixSizeThrows` | Matrix not 3×3 → assertion |
| `StripWrongMatrixSizeThrows` | Matrix not 2×3 → assertion |

---

## 9. Symratiospace (Symmetric Geometric Spacing)

### 9.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `SymratioEndpoints` | First element is `aXmin`, last is `aXmax` |
| `SymratioUniformWhenRatioOne` | `aRatio = 1.0` → uniform spacing |
| `SymratioMonotone` | `aX(i) < aX(i+1)` for all i |
| `SymratioSymmetric` | Spacing from left mirrors spacing from right: `aX(k) - aX(k-1) ≈ aX(N-k) - aX(N-k-1)` |
| `SymratioGeometricProgression` | Adjacent spacing ratios match `aRatio` for the graded portion |
| `SymratioRegressionBugM1` | `aXmin=2.0, aXmax=10.0, N=5` → midpoint should be 6.0, not 4.0 (BUG-M1 target) |
| `SymratioRegressionBugM1Symmetric` | Same setup → verify symmetry around midpoint |

### 9.2 Tests `[debug]`

| Test Name | What It Verifies |
|---|---|
| `SymratioEvenNThrows` | `aN = 4` → assertion ("aN must be odd") |

---

## 10. Error Mechanism Summary

| Function | Macro | Active In | Notes |
|---|---|---|---|
| `cardano()` unreachable branch | `BELFEM_ERROR` | Always | Nearly impossible to trigger in practice |
| `circle_from_points()` iteration guard | `BELFEM_ERROR` | Always | Fires if > 1000 iterations |
| `circle_from_points()` input length | `BELFEM_ASSERT` | Debug only | |
| `cubic_bezier()` size checks | `BELFEM_ASSERT` | Debug only | |
| `find_interval()` | None | — | No assertions at all |
| `quadratic_gradient()` index check | `BELFEM_ASSERT` | Debug only | Length check commented out |
| `rotation_matrix()` size checks | `BELFEM_ASSERT` | Debug only | |
| `symratiospace()` odd-N check | `BELFEM_ASSERT` | Debug only | |

---

## 11. Implementation Notes for Claude Code

1. **Root substitution over golden values.** For `cardano`, always verify by plugging roots back into the polynomial. This catches bugs that "expected ≈ -0.6823" comparisons miss.
2. **Beam poly round-trip.** Use `polyval` and `dpolyval` from the linalg module to evaluate the returned coefficients. This tests the beam-poly function through the same infrastructure the framework uses.
3. **BUG-M1 is the highest-priority regression.** The `symratiospace` midpoint formula looks wrong for non-zero `aXmin`. If the test passes, the bug analysis was wrong and the formula is compensated elsewhere. If it fails, you've found a real bug. Either way, the test is valuable.
4. **BUG-M2 requires careful setup.** The `find_interval` interior loop's output assignment issue is subtle. Use a 5+ element sorted array and test with values that land in the middle intervals, not just the edge cases handled before the loop.
5. **Bézier parametrization is ξ ∈ [-1, 1].** Do not assume standard t ∈ [0, 1] Bézier. At ξ=-1 and ξ=+1, evaluate the four basis functions by hand and verify they select the correct control points.
6. **Rotation matrix Blaze padding (BUG-M4).** Test the Euler-angle rotation matrix by reading values through `operator()(i,j)`, not through `data()`. This abstracts away the backend's memory layout.
7. **`find_interval` has no assertions.** It will silently produce wrong results for bad inputs (unsorted data, length < 2). Tests should only use valid inputs and verify outputs.

---

## 12. Codex Audit Checklist

When reviewing Claude Code's test implementation, verify:

- [ ] `cardano` tests use root substitution, not hardcoded approximate root values
- [ ] BUG-M1 regression test uses `aXmin != 0` (e.g., `[2, 10]`)
- [ ] BUG-M2 regression test uses a non-trivial interior interval
- [ ] BUG-M4 regression test uses `operator()` for Euler matrix, not `data()`
- [ ] Beam poly tests verify interpolation conditions via `polyval`/`dpolyval`
- [ ] Bézier endpoint tests verify against the actual ξ∈[-1,1] parametrization
- [ ] `quadratic_gradient` tests document the `n >= 3` precondition
- [ ] Rotation matrix tests check orthogonality and det ≈ +1
- [ ] Debug tests wrapped in `#ifndef NDEBUG` for `BELFEM_ASSERT`, no guard for `BELFEM_ERROR`
- [ ] `EXPECT_NEAR` used for all floating-point comparisons
- [ ] BELFEM naming conventions (`t` prefix for locals)
