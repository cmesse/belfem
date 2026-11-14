# Database Module Usage Guide {#physics_database_database_usage_guide}

**Date:** 2026-08-10
**Purpose:** What `Database` and `Projector` are for, how to use them, and why the projection
step exists at all
**Module:** `src/physics/database`

---

## 1. What this module is

A `Database` is a **precomputed lookup table over a structured tensor grid**, evaluated by
shape-function interpolation and persisted to HDF5. It moves an expensive property evaluation
from the nonlinear loop to setup.

Its main user is `physics/materials`: normal-metal resistivity depends on temperature, field
magnitude and field angle through a Kohler fit too costly to evaluate at every quadrature
point in every Newton iteration. The fit is sampled once on a three-dimensional grid,
projected, stored, then read by interpolation.

The module has exactly two classes:

| Class | Role |
|---|---|
| `Projector` | L2-projects an **already sampled** field: node values → B-spline coefficients → projected node values |
| `Database` | holds the finished table, evaluates it, saves and loads it |

**The `Projector` does not sample anything.** It consumes a mesh field the caller has already
filled (`mMesh->field( aField )->data()`); producing those samples is the caller's job. It is
build-time only: `Database`'s mesh constructor creates it, runs it, and destroys it, and
external code should not construct one.

## 2. The life cycle

```
      caller                        Projector                    Database
      (expensive, exact)                                         (cheap, approximate)

      f(T, log10 B, angle)          L2 projection                node values
      sampled onto a       ------>  onto a smooth      ------>   + interpolation
      tensor mesh field             B-spline basis               + HDF5
```

1. The caller creates a tensor mesh and fills a nodal field with sampled values.
2. `Database( Mesh *, field, aProject = true )` runs the `Projector` over that field.
3. The resulting node values live in `mValues` and are evaluated with
   `evaluate( x, y[, z] )`.
4. `save( hid_t )` writes the table to HDF5. The other two constructors load it again, so
   caching the table avoids rebuilding the samples on every run. All three constructors are
   collective: call them on every rank (§5).

## 3. Why there is a projection at all

**The projector is not a fitting convenience. It exists because the derivative accessors must
be usable.** The solve converts a raw nodal field into the closest smooth B-spline field in
the L2 metric; it does not try to improve the property model itself.

`Database` exposes `evaluate_derivx`, `evaluate_derivy` and `evaluate_derivz` as well as
`evaluate`, and material code uses those derivatives in Newton tangents — `drho/dT`, `drho/dB`
and related terms. If sampled values were stored directly and interpolated with the element's
Lagrange shape functions, the table would be only C⁰: values would stay continuous across
element boundaries, but **derivatives would jump at every one of them**. Tangents built from a
saw-tooth derivative degrade Newton convergence, while the values themselves still look
reasonable — which is what makes the failure hard to spot.

The projector removes the derivative jumps by doing a least-squares (L2) projection of the
sampled data onto a **B-spline basis** and storing the node values of *that* smooth
representation. In `compute_element_matrices` (`cl_DatabaseProjector.cpp`):

| Symbol | Meaning |
|---|---|
| `mMel` | element mass matrix in the Lagrange basis, `Σ_k w_k Nᵀ N`, scaled by `½ · element_step` per dimension (the tensor grid's Jacobian) |
| `mT` | the T-matrix from `TensorMeshFactory::t_matrix` — the element's B-spline basis evaluated at its Lagrange nodes, mapping B-spline **control-point** coefficients to Lagrange **nodal** values |
| `mBel` | `Tᵀ · M` — maps sampled nodal values to the B-spline-side right-hand side |
| `mAel` | `Bel · T = Tᵀ · M · T` — the B-spline mass matrix |

`project()` assembles `A` and `Y` over all elements, solves `A · X = Y` for the control-point
coefficients `X`, then writes back **nodal** values `Fel = mT · Xel`. The stored table is the
smooth B-spline field sampled at the grid nodes; at order 2 that representation is C¹, so the
derivative accessors return a continuous field.

If it seems odd that smoothness survives being stored as node values and evaluated with
Lagrange shape functions, note that those node values are not arbitrary samples: on each
element they are exactly `mT · Xel`, so the Lagrange polynomial reconstructs the extracted
B-spline element polynomial rather than approximating it.

Two consequences worth remembering:

- **The stored values are not the sampled values.** They are the projection of them. Expect
  small differences from a direct evaluation of the underlying property function. That
  discrepancy is the cost of the smoothness.
- **Reproducing the table bit-for-bit requires a deterministic solve.** The projection uses a
  direct sparse solver. A multithreaded solver may sum in a different order from one run to
  another, so two rebuilds of the same table can differ in their trailing digits. The module
  does not pin thread counts or request a deterministic mode. Compare rebuilt tables with a
  tolerance and treat such differences as solver noise, not as a defect.

## 4. Parallel build contract

**The build is collective in call and master-only in work, deliberately.**

Every rank needs the finished material, so `Database`'s constructor calls the `Projector` on
**all** ranks (`cl_Database.cpp`), and the non-projecting branch pairs `share` / `receive`
instead. But the *work* is not distributed:

- Rank 0 owns every node of the tensor work grid, assembles the system, and does the
  evaluation.
- The other ranks hold an element-less copy of the grid and join only where the solver
  requires it: `project()` calls `mSolver->solve()` on a worker **only if**
  `mSolver->wrapper()->uses_mpi()`. With a serial solver the workers simply wait at the
  barriers.
- Workers still allocate an empty `SpMatrix`; the MPI-solver path binds a reference to it
  before entering the solve. Keep that allocation.
- The finished node values reach every rank through the `share` / `receive` pair.

**Why it is built this way:** the projector's system is comparatively small, and every solver
it can select — MUMPS, PARDISO, SUPERLU, UMFPACK — assembles on the main process anyway.
Distributing the build would add more communication overhead than it saves in compute.
STRUMPACK is deliberately *not* offered here, for the same reason: the matrices are too small
to interest it.

The solver is chosen at compile time in the `Projector` constructor, in this order of
preference: **MUMPS → PARDISO → SUPERLU → UMFPACK**. If none is available the constructor
raises a hard error.

**Do not read the packing code as evidence of a distributed design.** A producer such as
`Metal::populate_rho_database` packs values densely over the nodes it owns, while `project()`
reads the field positionally by `node->index()`. These look like two index spaces; they are
one. The value vector is a full-length mesh field and `Mesh` assigns `node->set_index(...)`
walking the same container, so under single ownership slot *k* belongs to the node whose index
is *k*. The `owner()` filter in such a producer is a permanent no-op under this contract.

## 5. API

### Constructors

| Constructor | Use |
|---|---|
| `Database( Mesh *, field, aProject = true, material = "" )` | build from a sampled tensor mesh; `aProject = false` skips the projection and shares the raw nodal values, so the table is then only C⁰ |
| `Database( const string & path, const string & material )` | load from an HDF5 file, selecting the material's group |
| `Database( hid_t, const string & label )` | load from an already-open HDF5 group |

**All three constructors are collective.** Rank 0 reads the file or mesh, then
`broadcast`s the grid parameters and `share`s the values. Every other rank must enter the same
constructor to `receive` them (`Database::load`, `cl_Database.cpp`). Constructing a `Database`
on rank 0 only deadlocks. With the `hid_t` constructor, only rank 0's handle is read; other
ranks may pass any value. `save( hid_t )` is different: it returns immediately on every rank
except rank 0.

**A non-tensor mesh raises an error.** On the default `aProject = true` path, the error comes
from the projector through `TensorMeshFactory::create_bsplines` ("Mesh must be a tensor
mesh"). With `aProject = false`, it comes from the `aMesh->is_tensormesh()` check in the
`Database` constructor. This module works only on structured grids; that structure is the basis
of the O(1) element lookup in `evaluate`.

### Evaluation

```cpp
real v  = tDatabase.evaluate( x, y, z );          // 2-D overload also available
real dx = tDatabase.evaluate_derivx( x, y, z );   // ... derivy, derivz
real lo = tDatabase.min( 0 );                     // grid bounds per dimension
real hi = tDatabase.max( 0 );
```

| Method | Returns |
|---|---|
| `evaluate( x, y[, z] )` | interpolated table value |
| `evaluate_derivx` / `evaluate_derivy` (2-D and 3-D), `evaluate_derivz` (3-D) | first physical derivatives |
| `evaluate_derivxy( x, y, z )`, `evaluate_derivyz( x, y, z )` | mixed second derivatives (3-D) |
| `upper_bound_derivy_on_lower_y_face()` | an upper bound of ∂f/∂y over the whole face `y = min( 1 )`: the largest Bernstein coefficient of the face polynomial over all elements (3-D) |
| `lower_bound_on_lower_y_face()` | a lower bound of `f` over the face `y = min( 1 )` (3-D) |
| `dimension()` | 2 or 3 |
| `min( d )`, `max( d )` | grid bounds in dimension `d` |
| `values()` | read-only reference to the stored node values, in the stored scale |

All values and derivatives use the **stored** scale. The materials module stores `log(rho)` for
resistivity and log₁₀ for the jc/n tables, so the inverse transform depends on the table (§6).

`evaluate` locates the element through integer arithmetic on the grid config, maps the query
point to parametric coordinates, evaluates the Lagrange shape functions, and contracts them
with the stored node values. There is no mesh search; the remaining work is a fixed-size
shape-function contraction over the element's nodes.

**Clamping to `min`/`max` is the caller's responsibility, and skipping it fails silently.**
`TensorMeshConfig::element_ijk` (`cl_TensorMeshConfig.hpp`) floors the scaled coordinate,
clamps the result to the grid as a real, and only then converts it to `index_t`. No error is raised in a
release build:

- **Above `max`**, the index clamps to the last element.
- **Below `min`**, regardless of distance, the index clamps to element 0.
- **NaN** also maps to element 0 in a release build; a debug build stops at a `BELFEM_ASSERT`.

In every case, the selected element is evaluated at a parametric coordinate outside `[-1, 1]`;
that is, the table **extrapolates**. A polynomial extrapolated even slightly beyond its element
can diverge quickly, producing a plausible-looking but incorrect value. Every production
consumer clamps first—see the pattern in §6.

### Persistence

`save( hid_t )` writes the table into the given HDF5 group as the datasets `order`, `origin`,
`points`, `step` and `values`. The record contains **no format version or content stamp**, so a
consumer that caches tables on disk must determine for itself whether a file is still current.

The materials module does this with a revision number: `Metal` and `Alloy` write a `revision`
dataset next to the table, and `rho_database_is_current` (`fn_rho_database_is_current.hpp`)
accepts a cached file only if it contains an `RRR` dataset and a single integer `revision` equal
to `gRhoDatabaseRevision`. Any other file is rebuilt. A check that tests only for an expected
dataset could accept a stale table with the right shape.

### `db2exo`: viewing a table

`db2exo <hdf5file>` converts the tables in a database file into an Exodus file with the same
name and the `.exo` extension, with one node field per table for inspection in ParaView. A
table is any group containing all five datasets above. The first 2-D or 3-D table defines the
grid. Tables with a different grid or dimension are skipped with a message. Order-3 tables
are written with corner connectivity only. The target is built only when both `USE_HDF5` and
`USE_EXODUS` are on (`src/physics/database/CMakeLists.txt`).

## 6. Consumer contract: materials as the worked example

The grid for a metal's resistivity comes from `create_database_mesh` in `physics/materials`:
an order-2 tensor mesh with **95 × 35 × 37 nodes**. Note that `Mesh`'s second argument is
`aNumNodes`, not an element count, so at order 2 this gives **47 × 17 × 18 elements**
(`(num_nodes − 1) / order` per axis; the element *step* is `node_step × order`). The node
steps are `{ 4.0, 0.1, 5° }`, from the origin `{ 0.0, -2, 0.0 }`, over the axes of
temperature, log₁₀ field magnitude, and field angle. The covered ranges are:

| Axis | Range |
|---|---|
| temperature | 0 … 376 K |
| log₁₀ B | −2 … 1.4, i.e. **B from 0.01 to ≈ 25 T** |
| angle | 0 … π rad (180°) |

These are the bounds reported by `min( d )` and `max( d )`, and they are the bounds to which a
consumer must clamp. `min` and `max` use the table's own units, so the angle bound is
`max( 2 )` = π, not 180.

The query side (`Metal::rho_table`) demonstrates the three conventions a consumer must follow:

```cpp
real theta  = std::clamp( T, mDatabaseTmin, mDatabaseTmax );              // 1. clamp
real log10B = std::log( std::clamp( B, mDatabaseBmin, mDatabaseBmax ) ) * mInvLog10 ;
real angle  = beta < 0 ? beta + constant::pi                              // 2. fold
                       : beta > constant::pi ? beta - constant::pi : beta ;

return std::exp( mRhoData->evaluate( theta, log10B, angle ) );            // 3. undo log
```

1. **Clamp every coordinate** to the grid bounds before evaluation.
2. **Match the axis transform** used by the sampler—here, log₁₀ for field magnitude and an
   angle in `[0, π]`. The fold is a single shift by π: it maps `[-π, 2π]` into `[0, π]` and
   leaves an angle outside that range outside the grid, where the table extrapolates. The field
   angles computed by `Calculator` (`bj_angle_2d`, `bj_angle_3d`, `bn_angle`) already lie
   within `[0, π]`.
3. **Undo the value transform.** This table stores `log(rho)`, not `rho`. Storing the logarithm
   keeps a quantity spanning orders of magnitude better conditioned for the projection;
   consumers exponentiate on return. Derivatives from such a table are derivatives of
   `log(rho)` and require the matching chain-rule conversion.

## 7. Pitfalls

| Pitfall | Symptom | Avoid by |
|---|---|---|
| Passing a non-tensor mesh | hard error, from `TensorMeshFactory::create_bsplines` when projecting, or from the `Database` constructor's `is_tensormesh()` check when `aProject = false` | build the grid with a `TensorMeshFactory` / `create_database_mesh`-style helper |
| Querying outside the grid | **no error** — silent extrapolation off the nearest edge element, diverging fast | clamp to `min`/`max` first — always the caller's job |
| Forgetting the value transform | resistivity off by orders of magnitude | remember tables may store a transformed quantity (`log` for rho, log₁₀ for jc/n) |
| Assuming stored == sampled | small mismatches against a direct property call | the table is the **projection** of the samples (§3) |
| Expecting bit-identical rebuilds | spurious "regression" on a rebuilt table | compare with a tolerance (§3) |
| Constructing a `Database` on rank 0 only | deadlock in `broadcast` | call every constructor on all ranks (§5) |
| Trying to parallelize the build | more communication, no speed-up | the master-only design is deliberate (§4) |
| Removing the workers' empty `SpMatrix` | null dereference on the MPI-solver worker path | leave it; workers bind it when the selected solver uses MPI |
| Trusting a cached HDF5 file by name alone | stale table from an older law, loaded silently | store and check a content revision, as `rho_database_is_current` does (§5) |

## 8. Related documentation

- `src/physics/materials/doc/materials_contracts_and_invariants.md` — the ownership and
  database contracts on the consumer side, including the master-build invariant.
- `src/physics/materials/doc/materials_usage_guide.md` — how material property lookups are
  wired end to end.
- `doc/coding_philosophy.md` — the MPI and container rules this module follows.
