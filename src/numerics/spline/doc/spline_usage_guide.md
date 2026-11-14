# Spline Usage Guide {#numerics_spline_spline_usage_guide}

**Date:** 2026-06-14
**Purpose:** Practical guide to BELFEM's 1-parameter cubic `Spline` and its relation
to the multi-parameter B-spline lookup tables used for property databases.
**Module:** `numerics/spline` (1-D `Spline`); related: `physics/database` (2-D/3-D tables)

---

## 1. What this module provides

| Need | Use | Defined in |
|------|-----|-----------|
| Smooth 1-parameter curve `f(x)` (e.g. `c_p(T)`, `J_c(T)`, a B–H segment) | `Spline` | `cl_Spline.hpp`, `cl_Spline.cpp` |
| Cubic through 2 values + 2 slopes | `create_truss_poly()` | `fn_Create_Truss_Poly.hpp` |
| Quartic through 5 conditions ("glue" between segments) | `create_glue_poly()` | `fn_Create_Glue_Poly.hpp` |
| 2-/3-parameter property table `f(x,y[,z])` | `Database` (separate module) | `physics/database/cl_Database.hpp` — see §6 |

The 1-D `Spline` is a **cubic** spline on a **uniform** abscissa grid. It is the
workhorse for tabulated 1-parameter material/fluid properties.

---

## 2. The 1-D cubic Spline

### Type and storage

- **Cubic (3rd order).** Each interval stores 4 coefficients. Evaluation uses the cubic Horner polynomial `((a·x + b)·x + c)·x + d` (`Spline::eval`, `cl_Spline.hpp`; `mData.set_size( 4, mNumberOfPoints )` in `Spline::create_coeffs`, `cl_Spline.cpp`).
- **Uniform grid required.** The spline uses one step, `mDeltaX`. Equidistance is checked during construction through `BELFEM_ASSERT` (`Spline::check_input`, `cl_Spline.cpp`). The check is therefore **compiled out in a release build**: a non-uniform grid then produces silently wrong interpolation instead of an error. The interval is found arithmetically, not by searching (`Spline::find_col`, `cl_Spline.hpp`), so lookup is O(1).
- **At least four points.** `create_helpmatrix` and the constructors require `X.length() > 3`. This check also uses `BELFEM_ASSERT`, so it is absent in a release build.
- **Continuity: C².** The coefficients come from one global system for the nodal first derivatives (`Spline::initialize`, `cl_Spline.cpp`). This gives the classical global cubic spline with continuous second derivative for natural and clamped end conditions.

### Boundary conditions

`enum class spline::SplineBC` (`Spline_Enums.hpp`):

| Value | Meaning |
|-------|---------|
| `NoCurvature` | **Natural** spline — 2nd derivative = 0 at the endpoint (default) |
| `Parabolic` | Parabolic runout — end polynomial's cubic term = 0 |
| `Tangent` | **Clamped** — prescribed first derivative (`dYdX0` / `dYdX1`) |

The start and end conditions are independent. They appear in two places: the first and last rows of the help matrix (`create_helpmatrix`), and the first and last right-hand-side entries built from the conditions passed to the constructor or `update_data`.

**The two must use the same pair.** If a help matrix built for `NoCurvature` ends is used with a `Tangent` spline, the solve uses the wrong end equations and produces wrong coefficients without reporting an error.

### Construction

The constructor takes a **pre-built help matrix** (`SpMatrix`). The matrix depends only on the number of points, the step, and the boundary-condition pair. Splines that share all three can therefore share one matrix. Each solve factorizes the matrix again and frees the factorization afterward. What is shared is the assembled matrix, not a factorization.

```cpp
#include "cl_Spline.hpp"
using namespace belfem;

// tX must be uniformly spaced, with at least four points; tY are the sampled values
Vector< real > tX = ...;   // abscissae (equidistant)
Vector< real > tY = ...;   // ordinates

// 1) build the help matrix once for this grid and this boundary-condition pair
SpMatrix tA;
spline::create_helpmatrix( tX.length(), tX( 1 ) - tX( 0 ), tA,
                           spline::SplineBC::NoCurvature,
                           spline::SplineBC::NoCurvature );

// 2) construct the spline with the same boundary conditions
Spline tSpline( tX, tY, tA,
                spline::SplineBC::NoCurvature,
                spline::SplineBC::NoCurvature );
```

The data constructors are declared in `cl_Spline.hpp`. One accepts boundary conditions, slopes, `aXref`, `aSref`, and `aMasterProc`; another takes `aXref` directly after the help matrix and uses natural ends.

Other constructors create a spline from an HDF5 file or group (`Spline( aFile, aLabel, aMaster )`, `Spline( aGroup, aMaster )`), receive one through `Spline( aMasterProc )`, or create an empty container with `Spline( aN, aXmin, aXmax )`.

The file constructors read on the master rank and `synchronize` to the others. The `X`/`Y` constructors do not synchronize by default (`aMasterProc = gNoOwner`): each calling rank solves its own coefficients, and nothing is communicated. Pass a master rank to make that rank solve and broadcast. The other ranks then construct with `Spline( aMasterProc )`.

The constructors solve with SuperLU when it is linked, otherwise UMFPACK, otherwise `gDefaultSolver`.

### Re-fitting with `update_data`

To fit new `Y` values on the **same** grid without rebuilding the help matrix, call:

`update_data( tA, tYnew, aStartBC, aEndBC, adYdX0, adYdX1, aXref, aSref )`

(`Spline::update_data`, `cl_Spline.hpp`). Use the same boundary conditions that built the help matrix.

Important details:

- It runs only on rank 0 and does not broadcast. Other ranks retain the old table.
- It always solves with SuperLU, so the build must use `USE_SUPERLU` (ON by default). Without SuperLU, the solver constructor stops with an error.
- It rebuilds the 4-row coefficient table and resets the extra mode to `None` (`Spline::create_coeffs`, `cl_Spline.cpp`). It rebuilds the entropy row only when `aXref > 0`. If the spline had an entropy row, pass `aXref` and `aSref` again. If it had an integral row, call `create_integral()` again afterward. Otherwise, `entropy()`/`integrate()` behave as on a spline without the mode: a `BELFEM_ASSERT` in a debug build and an out-of-bounds read in a release build.

### Evaluation

All evaluation functions are in `cl_Spline.hpp`:

| Call | Returns |
|------|---------|
| `eval(x)` | `f(x)` |
| `deval(x)` | `f'(x)` |
| `ddeval(x)` | `f''(x)` |
| `eval(x, col)`, `deval(x, col)`, `ddeval(x, col)` | the same, with the interval index supplied (no `find_col`) |
| `entropy(x)` / `dentropy(x)`, also with `col` | entropy mode (see below) |
| `integrate(x)` / `integrate(x0,x1)` | antiderivative (integral mode) |

The `col` overloads let callers evaluating several quantities at the same `x` find the interval once (`find_col( x )`) and reuse it. The gas models use this pattern.

```cpp
real f   = tSpline.eval(   1.5 );
real df  = tSpline.deval(  1.5 );   // first derivative — needed for Newton-Raphson
real ddf = tSpline.ddeval( 1.5 );
```

### Extrapolation behavior (important)

Out-of-range queries do **not** throw and are **not** clamped to a constant. The interval *index* is clamped (`find_col`), but the boundary interval's **cubic polynomial is evaluated at the actual `x`** (`eval`). Thus, beyond `[x_min, x_max]`, the spline gives a cubic extension of the edge segment. Guard the input range yourself if you do not want this behavior.

### Thermo extra modes

`enum class spline::ExtraMode { None, Entropy, Integral }` (`Spline_Enums.hpp`) records whether a 5th coefficient row is present (`Spline::add_row_to_data`, `cl_Spline.cpp`). No constructor takes an `ExtraMode`; the row is built in one of the following ways.

- **Entropy** — built by a data constructor or `update_data` called with `aXref > 0`. The row is offset so that `entropy( aXref ) == aSref`:

  ```cpp
  Spline tSpline( tT, tH, tA, tTref, tSref );   // tH: enthalpy samples h(T)
  ```

  `dentropy(x)` equals `deval(x) / x`, and `entropy(x)` contains `log x`. The stored cubic must therefore represent the **enthalpy** `h(T)`, so that `entropy()` integrates `c_p / T`. If you fit `c_p` samples instead, the result integrates `c_p' / T` without complaint. The abscissa must be positive (absolute temperature); nothing checks this.
- **Integral** — `create_integral( aXref, aYref )` (`cl_Spline.hpp`) adds an analytic antiderivative row. It is offset so that `integrate( aXref ) == aYref`; `aXref` defaults to `x_min`. `integrate(x0,x1)` is then O(1).

A spline holds at most one extra row. `create_integral()` requires a 4-row table, so it cannot be called on a spline that already has an entropy row.

`set_extra_mode()` only records the flag. It is intended for callers that write coefficients directly through `matrix_data()`. It does not build a row. Calling `entropy()` or `integrate()` without the matching mode triggers a `BELFEM_ASSERT` in a debug build and reads past the table in a release build.

### Persistence

`save`/`load` for an HDF5 group and `save_to_database( aDatabase, aLabel )` (declared in `cl_Spline.hpp`) both store the extra mode. `save_to_database` also writes an `order` entry, always `0` (`Spline::save_to_database`, `cl_Spline.cpp`); `load` does not read it.

---

## 3. Quick recommendations

- **One help matrix per grid and boundary-condition pair.** Build it once with `create_helpmatrix` and share it across splines that use the same grid and end conditions. Use `update_data` to re-fit, passing the same end conditions.
- **Choose the end condition deliberately.** `NoCurvature` (natural) is the default. Use `Tangent` when the boundary slope is known, such as for a known asymptote.
- **Guard the domain.** Extrapolation is a cubic extension, not a clamp.
- **Need a derivative for a Newton solve?** Use `deval`. It is the analytic derivative of the same cubic and is therefore consistent with `eval` (no finite differencing).

---

## 4. Common pitfalls

| Pitfall | Symptom | Fix |
|---------|---------|-----|
| Non-uniform `X` | `BELFEM_ASSERT` at construction in debug builds; silently wrong values in release (`Spline::check_input`) | Resample onto a uniform grid first |
| Fewer than four points | `BELFEM_ASSERT` in debug builds; not checked in release | Supply at least four samples |
| Help matrix built for other end conditions | Wrong coefficients, no error | Build the help matrix with the same `SplineBC` pair as the spline |
| Relying on out-of-range values | Silent cubic extrapolation, not clamping | Clamp/validate `x` before `eval` |
| Calling `entropy()`/`integrate()` without the mode, including after `update_data()` dropped the extra row | `BELFEM_ASSERT` (debug) / out-of-bounds read (release) | For entropy, construct or re-fit with `aXref > 0`; for integration, call `create_integral` again after `update_data()` |
| Expecting coefficients on all ranks | `X`/`Y` constructors with the default `aMasterProc = gNoOwner` solve per rank; `update_data()` runs on rank 0 only and does not broadcast | Pass `aMasterProc` (other ranks: `Spline( aMasterProc )`); after `update_data()` re-create the spline |

---

## 5. Testing

`tests/math/test_Spline.cpp` exercises construction, evaluation, boundary
conditions, the integral and entropy modes, and the `TrussPoly`/`GluePoly` helpers.

---

## 6. Related: multi-parameter B-spline lookup tables

For **2-parameter** and **3-parameter** properties, BELFEM uses B-spline lookup tables implemented in the **`physics/database`** module (`cl_Database`), not in this spline class. The two are easy to confuse, so this section gives a short orientation.

### Orders supported

BELFEM supports **2nd-order (quadratic)** and **3rd-order (cubic)** B-spline lookup tables. A 1st-order/linear form also exists:

| Order | 2-D element | 3-D element | Evals / integration point (3-D) | Continuity |
|-------|-------------|-------------|---------------------------------|-----------|
| 1 (linear) | QUAD4 | HEX8 | 8 | C⁰ |
| **2 (quadratic)** | **QUAD9** | **HEX27** | **27** | **C¹** |
| 3 (cubic) | QUAD16 | HEX64 | 64 | C² |

The continuity column applies to projected tables, as described in the next section.

### Why this is genuinely a B-spline table (continuity)

The table stores **nodal values** and evaluates them at runtime with tensor-product Lagrange shape functions (`Database::evaluate` and the `eval_quad*`/`eval_hex*` members, `cl_Database.hpp`).

The nodal values come from a **B-spline projection**: the construction (`database::Projector` + `TensorMeshFactory::t_matrix`) samples the fitted B-spline at the element's Lagrange nodes. Within each element, the B-spline is a single degree-`p` polynomial. Lagrange interpolation at `(p+1)` nodes reproduces that polynomial **exactly**. Runtime evaluation therefore reconstructs the B-spline and inherits its **Cᵖ⁻¹ continuity** across elements: C¹ for quadratic and C² for cubic. The uniform grid also preserves cheap O(1) element lookup.

This argument requires projection. `Database( aMesh, aField, aProject, aMaterial )` projects by default (`aProject = true`). With `aProject = false`, the raw field values are used as nodal values, and the table is guaranteed only to be C⁰.

Like `Spline`, `Database::evaluate` clamps the element index but evaluates the edge element at the actual coordinate. Out-of-range queries therefore extrapolate instead of clamp.

Besides `evaluate`, `cl_Database.hpp` provides first derivatives (`evaluate_derivx`/`_derivy`, and `_derivz` in 3-D), the mixed 3-D derivatives `evaluate_derivxy` and `evaluate_derivyz`, and two certified bounds over the face `y = min`: `upper_bound_derivy_on_lower_y_face()` and `lower_bound_on_lower_y_face()`.

### Design guidance: prefer 2nd order (quadratic) at FEM runtime

For nonlinear FEM, **2nd-order tables (HEX27) are usually the better choice than 3rd-order tables (HEX64)**:

- **Cost.** A 3-D evaluation sums over all element nodes at each integration point: **27** shape-function evaluations for HEX27 versus **64** for HEX64. During a nonlinear solve, this occurs at every integration point and every Newton iteration. The 27-vs-64 difference dominates.
- **What the solver actually needs.** Newton–Raphson needs a **continuous first derivative**—a continuous, consistent Jacobian—for good convergence. Quadratic B-splines already provide **C¹** continuity, which satisfies that requirement. The additional **C²** continuity of cubic tables is rarely worth the approximately 2.4× evaluation cost in practice.

**Rule of thumb:** use **2nd-order (C¹)** lookup tables for nonlinear material/fluid properties unless a specific need for continuous curvature (C²) is demonstrated.

> **Test coverage (3rd-order construction).** The cubic T-matrix builders
> (`TensorMeshFactory::compute_tmatrix_quad16` / `_hex64`, `cl_TensorMeshFactory.cpp`)
> evaluate the uniform cubic basis `bspline_cub`, ordered left to right like
> `bspline_quad`. `tests/mesh/test_TensorMeshBSpline.cpp` builds QUAD16 and HEX64
> tensor meshes through `Mesh` and checks the node-to-control-point weights, which are
> rows of the same T-matrix: partition of unity, use of every control point, and exact
> reproduction of x, x² and x³. The test does not run the `Projector` or evaluate a
> `Database` table, so the C² continuity of a projected cubic table follows from the
> argument above but is not tested directly.

For the full description of the database tables, see the `physics/database` module guide, `src/physics/database/doc/database_usage_guide.md`.
