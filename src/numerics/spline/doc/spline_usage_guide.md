# Spline Usage Guide {#numerics_spline_spline_usage_guide}

**Date:** 2026-06-14
**Purpose:** Practical guide to BELFEM's 1-parameter cubic `Spline` and its relation
to the multi-parameter B-spline lookup tables used for property databases.
**Module:** `numerics/spline` (1-D `Spline`); related: `physics/database` (2-D/3-D tables)

---

## 1. What this module provides

| Need | Use | Defined in |
|------|-----|-----------|
| Smooth 1-parameter curve `f(x)` (e.g. `c_p(T)`, `J_c(T)`, a B–H segment) | `Spline` | `cl_Spline.hpp/.cpp` |
| Cubic through 2 values + 2 slopes | `create_truss_poly()` | `fn_Create_Truss_Poly.hpp` |
| Quartic through 5 conditions ("glue" between segments) | `create_glue_poly()` | `fn_Create_Glue_Poly.hpp` |
| 2-/3-parameter property table `f(x,y[,z])` | `Database` (separate module) | `physics/database/cl_Database.hpp` — see §6 |

The 1-D `Spline` is a **cubic** spline on a **uniform** abscissa grid. It is the
workhorse for tabulated 1-parameter material/fluid properties.

---

## 2. The 1-D cubic Spline

### Type and storage

- **Cubic (3rd order).** Each interval stores 4 coefficients; evaluation is a cubic
  Horner polynomial `((a·x + b)·x + c)·x + d` (`cl_Spline.hpp:432-440`,
  `mData.set_size(4, n)` in `cl_Spline.cpp:571`).
- **Uniform grid required.** A single step `mDeltaX` is used. Equidistance is checked at
  construction, but through `BELFEM_ASSERT` (`cl_Spline.cpp:479`), so the check is **compiled
  out in a release build** — a non-uniform grid then produces silently wrong interpolation
  rather than an error. The interval is found by
  arithmetic, not a search (`find_col`, `cl_Spline.hpp:418-426`) — O(1) lookup.
- **Continuity: C².** Coefficients come from solving one global system for the
  nodal first derivatives (`cl_Spline.cpp:303-316`), i.e. the classical global
  cubic spline with continuous second derivative (for the natural/clamped end
  conditions).

### Boundary conditions

`enum class spline::SplineBC` (`Spline_Enums.hpp`):

| Value | Meaning |
|-------|---------|
| `NoCurvature` | **Natural** spline — 2nd derivative = 0 at the endpoint (default) |
| `Parabolic` | Parabolic runout — end polynomial's cubic term = 0 |
| `Tangent` | **Clamped** — prescribed first derivative (`dYdX0` / `dYdX1`) |

Start and end conditions are set independently.

### Construction

The constructor takes a **pre-built help matrix** (`SpMatrix`) so that many splines
sharing the same grid can reuse one factorization:

```cpp
#include "cl_Spline.hpp"
using namespace belfem;

// X must be uniformly spaced; Y are the sampled values
Vector<real> X = ...;   // abscissae (equidistant)
Vector<real> Y = ...;   // ordinates

// 1) build the help matrix once for this grid
SpMatrix A;
spline::create_helpmatrix( X.length(), X(1) - X(0), A,
                           spline::SplineBC::NoCurvature,
                           spline::SplineBC::NoCurvature );   // cl_Spline.hpp:27-33

// 2) construct the spline (natural BCs shown)
Spline tSpline( X, Y, A,
                spline::SplineBC::NoCurvature,
                spline::SplineBC::NoCurvature );              // cl_Spline.hpp:81-90
```

Other constructors: from an HDF5 file/group (`Spline(file, label, master)`,
`cl_Spline.hpp:94-101`), and parallel/empty forms. The file constructors read on
the master rank and `synchronize` to the others. The `X`/`Y` constructors above do
not by default (`aMasterProc = gNoOwner`): each calling rank solves its own
coefficients and nothing is communicated. Pass a master rank to have that rank
solve and broadcast; the other ranks then construct with `Spline( aMasterProc )`.
`update_data()` runs on rank 0 only and does not broadcast.

To re-fit new `Y` on the **same** grid without rebuilding the help matrix, call
`update_data(A, newY, ...)` (`cl_Spline.hpp:254-263`).

### Evaluation

| Call | Returns | Source |
|------|---------|--------|
| `eval(x)` | `f(x)` | `cl_Spline.hpp:432-441` |
| `deval(x)` | `f'(x)` | `:445-453` |
| `ddeval(x)` | `f''(x)` | `:479-485` |
| `eval(x, col)` | `f(x)` with the interval index supplied | `:457-465` |
| `entropy(x)` / `dentropy(x)` | thermo entropy modes (see below) | `:489-513` |
| `integrate(x)` / `integrate(x0,x1)` | antiderivative | `:515-537` |

```cpp
real f   = tSpline.eval(   1.5 );
real df  = tSpline.deval(  1.5 );   // first derivative — needed for Newton-Raphson
real ddf = tSpline.ddeval( 1.5 );
```

### Extrapolation behavior (important)

Out-of-range queries do **not** throw and are **not** clamped to a constant: the
interval *index* is clamped, but the boundary interval's **cubic polynomial is
evaluated at the actual `x`** (`cl_Spline.hpp:425-426` then `:435-440`). So beyond
`[x_min, x_max]` you get a cubic extension of the edge segment. Guard the input
range yourself if that is not what you want.

### Thermo extra modes

`enum spline::ExtraMode { None, Entropy, Integral }` adds a 5th coefficient row
(`add_row_to_data`, `cl_Spline.cpp:852`) for:
- **Entropy** — `entropy()`/`dentropy()` use a `c_p/T`-style form with `log x`
  (`cl_Spline.hpp:489-513`); built for gas heat polynomials.
- **Integral** — `create_integral(...)` (`cl_Spline.hpp:287`) precomputes an
  analytic antiderivative so `integrate(x0,x1)` is O(1).

These are opt-in; calling them without the matching mode trips a `BELFEM_ASSERT`.

### Persistence

`save`/`load` (HDF5 group), `save_to_database(file, label)` (`cl_Spline.hpp:222-240`).
An `order` tag is written, with `0` reserved as a special indicator
(`cl_Spline.cpp:828-829`).

---

## 3. Quick recommendations

- **One help matrix per grid.** Build it once with `create_helpmatrix` and share it
  across all splines on that grid; use `update_data` to re-fit.
- **Pick the end condition deliberately.** `NoCurvature` (natural) is the default;
  use `Tangent` when you know the boundary slope (e.g. a known asymptote).
- **Guard the domain.** Remember extrapolation is a cubic extension, not a clamp.
- **Need a derivative for a Newton solve?** Use `deval` — it is the analytic
  derivative of the same cubic, so it is consistent with `eval` (no finite
  differencing).

---

## 4. Common pitfalls

| Pitfall | Symptom | Fix |
|---------|---------|-----|
| Non-uniform `X` | `BELFEM_ASSERT` at construction in debug builds; silently wrong values in release (`cl_Spline.cpp:479`) | Resample onto a uniform grid first |
| Relying on out-of-range values | Silent cubic extrapolation, not clamping | Clamp/validate `x` before `eval` |
| Calling `entropy()`/`integrate()` without the mode | `BELFEM_ASSERT` (debug) / garbage (release) | Construct with the right `ExtraMode` / call `create_integral` |
| Expecting coefficients on all ranks | `X`/`Y` constructors with the default `aMasterProc = gNoOwner` solve per rank; `update_data()` runs on rank 0 only and does not broadcast | Pass `aMasterProc` (other ranks: `Spline( aMasterProc )`); after `update_data()` re-create the spline |

---

## 5. Testing

`tests/math/test_Spline.cpp` exercises construction, evaluation, boundary
conditions, the integral and entropy modes, and the `TrussPoly`/`GluePoly` helpers.

---

## 6. Related: multi-parameter B-spline lookup tables

For **2-parameter** and **3-parameter** properties, BELFEM uses B-spline lookup
tables implemented in the **`physics/database`** module (`cl_Database`), not in this
spline class. A short orientation, because the two are easy to confuse:

### Orders supported

BELFEM supports **2nd-order (quadratic)** and **3rd-order (cubic)** B-spline lookup
tables (a 1st-order/linear form also exists):

| Order | 2-D element | 3-D element | Evals / integration point (3-D) | Continuity |
|-------|-------------|-------------|---------------------------------|-----------|
| 1 (linear) | QUAD4 | HEX8 | 8 | C⁰ |
| **2 (quadratic)** | **QUAD9** | **HEX27** | **27** | **C¹** |
| 3 (cubic) | QUAD16 | HEX64 | 64 | C² |

### Why this is genuinely a B-spline table (continuity)

The table stores **nodal values** and evaluates them at runtime with tensor-product
Lagrange shape functions (`cl_Database.hpp:108-363`). This still yields B-spline
smoothness because the nodal values are produced by a **B-spline projection**: the
construction (`database::Projector` + `TensorMeshFactory::t_matrix`) samples the
fitted B-spline at the element's Lagrange nodes. Within each element the B-spline is
a single degree-`p` polynomial, and Lagrange interpolation at `(p+1)` nodes
reproduces that polynomial **exactly** — so the runtime evaluation reconstructs the
B-spline, inheriting its **Cᵖ⁻¹ continuity** across elements (C¹ for quadratic, C²
for cubic). The cheap O(1) element lookup of the uniform grid is preserved.

### Design guidance: prefer 2nd order (quadratic) at FEM runtime

For nonlinear FEM, **2nd-order tables (HEX27) are usually the better choice over
3rd-order (HEX64)**:

- **Cost.** A 3-D evaluation sums over all element nodes per integration point —
  **27** shape-function evaluations for HEX27 vs **64** for HEX64. During a
  nonlinear solve this happens at every integration point, every Newton iteration;
  the 27-vs-64 difference dominates.
- **What the solver actually needs.** Newton–Raphson needs a **continuous first
  derivative** (a continuous, consistent Jacobian) for good convergence. Quadratic
  B-splines already provide **C¹** continuity, which satisfies that requirement.
  The extra **C²** continuity of cubic tables is rarely worth the ~2.4× evaluation
  cost in practice.

**Rule of thumb:** use **2nd-order (C¹)** lookup tables for nonlinear material/fluid
properties unless a specific need for continuous curvature (C²) is demonstrated.

> **Status note (3rd-order construction).** A code review of the cubic T-matrix
> builders (`TensorMeshFactory::compute_tmatrix_quad16` / `_hex64`) found they call
> the *quadratic* basis helper `bspline_quad`, while the cubic helper `bspline_cub`
> is defined but never called. This should be verified before relying on 3rd-order
> tables; it is one more reason to prefer 2nd order in the meantime. (See the
> property-core architecture notes for details.)

For the full description of the database tables, see the `physics/database` module
and the B2a property-core architecture notes.
