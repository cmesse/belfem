# Spline Module Documentation {#numerics_spline_index}

**Module:** `numerics/spline`

1-parameter cubic spline (`Spline`) on a uniform grid, plus auxiliary polynomial
fitters. Used throughout the property databases (gas tables, materials) for smooth
1-D property curves.

## Documents

| Document | Contents |
|----------|----------|
| [spline_usage_guide.md](spline_usage_guide.md) | How to build/evaluate a `Spline`, boundary conditions, extrapolation, thermo extra modes, pitfalls; and orientation on the related 2-/3-parameter B-spline lookup tables (`physics/database`) |

## Quick reference

| Class / function | Purpose | Source |
|------------------|---------|--------|
| `Spline` | Cubic spline, uniform grid, C² (natural/parabolic/clamped BCs) | `cl_Spline.hpp/.cpp` |
| `spline::create_helpmatrix` | Build the shared help matrix for a grid | `cl_Spline.hpp:27-33` |
| `create_truss_poly` | Cubic through 2 values + 2 slopes | `fn_Create_Truss_Poly.hpp` |
| `create_glue_poly` | Quartic through 5 conditions | `fn_Create_Glue_Poly.hpp` |

## Related modules

- `physics/database` — 2-/3-parameter B-spline lookup tables (`Database`).
- `tests/math/test_Spline.cpp` — unit tests.
