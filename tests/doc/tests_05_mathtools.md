# Math Tools Module Test Suite

**Date:** 2026-03-23

---

## Test Count by Suite

| Suite | Tests | Coverage |
|-------|-------|----------|
| `Sign` | 6 | Positive, negative, zero, -0, infinity, integer |
| `Cardano` | 14 | Linear fallback, constant, quadratic (2 roots, no roots, double root), all-zero, D>0, D<0 (3 variants), triple root, double root, sorted output |
| `Circle` | 3 | Unit circle, translated, equidistant |
| `CircleDebug` | 2 | Wrong X/Y length |
| `BeamPoly` | 9 | Cubic (endpoints, derivatives, known coeffs, linear), fifth-order (values, 1st/2nd deriv, non-zero endpoints, known coeffs) |
| `Bezier` | 6 | Endpoints xi=-1/+1, straight line, midpoint, dust cleanup, derivative finite-diff |
| `BezierDebug` | 3 | Wrong point/column/work length |
| `FindInterval` | 4 | Below range, above range, exact knot, interior value |
| `QuadGrad` | 5 | Exact quadratic (interior/left/right), linear, minimum 3 points |
| `QuadGradDebug` | 1 | Out-of-bounds index |
| `Rotation` | 12 | Axis-angle (identity, 90-Z, ortho, det), Euler (identity, roll/yaw/pitch, ortho, det), strip |
| `RotationDebug` | 4 | Wrong axis length, wrong matrix sizes |
| `Symratio` | 7 | Endpoints, uniform, monotonicity, step symmetry, point symmetry, non-zero offset (2) |
| `SymratioDebug` | 1 | Even N |
| **Total** | **77** | |

---

## @warning API Gotchas

- `cardano()` contract: returns distinct real roots only, sorted.
- Bezier uses xi in [-1, 1], NOT standard t in [0, 1]. At xi=-1: selects column 0. At xi=+1: selects column 3.
- `symratiospace` filename is misspelled in source: `fn_symrationspace.hpp`.
- `find_interval()` has NO assertions — will silently produce wrong results for bad inputs.
- `rotation_matrix()` Euler: DIN 9300 yaw-pitch-roll convention.
- Beam poly coefficients in descending power order: `{a_n, ..., a_0}`, matching `polyval`.
