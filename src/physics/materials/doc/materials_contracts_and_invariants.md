# BELFEM Materials Module - Contracts and Invariants {#physics_materials_materials_contracts_and_invariants}

**Module:** `src/physics/materials`
**Purpose:** The rules whose violation produces memory errors, deadlocks or physically wrong results, stated as contracts with their reasons
**Date:** 2026-08-25 (first version 2026-01-16)
**Revision:** 2.0

---

## Overview

The usage guide explains how to call the module. This document lists what must be true for those calls to be safe. Each contract names the code that enforces it, or states plainly that nothing does.

---

## Ownership and Lifetime

### 1. Property-function ownership transfers on assignment

> `set_jc_function( const JcFunction * )`, `set_n_function( const JcFunction * )`,
> `load_bh_curve( const BhCurve * )` and `set_bh_curve( const BhCurve * )` hand the object to the
> material permanently. The caller must never delete it. Use `load_bh_curve`: it also routes
> `H`, `mu` and `dmudH` to the curve, which `set_bh_curve` does not.

The destructors enforce this: `Material::~Material` deletes the J_c and n functions (`cl_Material.cpp:174`), and `Metal::~Metal` deletes the B-H curve (`cl_Material_Metal.cpp:61`). A second caller-side `delete` is a double free.

```cpp
material::JcFunction * tJc = tFactory.create_jc_function( 1e9, 5.0, 0.5, 2.0 );
tYbco->set_jc_function( tJc );      // tYbco owns tJc from here
delete tYbco ;                      // deletes tJc
```

### 2. The factory owns nothing it returns

> `MaterialFactory::create_material`, `create_bh_curve` and `create_jc_function` return
> caller-owned pointers. The label-based path keeps no registry; the input-deck constructor
> fills a label → pointer map (`mMaterialsMap`) for the kernel's lookup, and the factory
> destructor deletes nothing in either case.

| Call | Returns | Owner until `set_*()` | Owner after |
|---|---|---|---|
| `create_material( … )` | `Material *` | caller | caller — always |
| `create_bh_curve( … )` | `BhCurve *` | caller | the material, after `load_bh_curve` / `set_bh_curve` |
| `create_jc_function( … )` | `JcFunction *` | caller | the material, after `set_jc_function` / `set_n_function` |

An `Alloy` owns the component `Metal` objects the factory built for it.

### 3. Materials outlive everything that references them

> Mesh blocks, FEM kernels, DOF managers and calculators hold raw pointers to materials.
> Destroy materials last — after the kernel and the mesh.

Nothing enforces this; a dangling material pointer surfaces as a crash in assembly.

---

## Physical Contracts

### 4. HTS: `rho()` and `rho_powerlaw()` are different physics

> `rho( T )` is the normal-state resistivity of the matrix. The superconducting E–J power law
> is `rho_powerlaw( normJ, T, normB, angleNxB )`. Choosing one is a modeling decision, not a
> convenience.

The assembly path always calls the full T-bearing overload and lets the material route on its
declared dependencies (usage guide §7.3). All arguments are `real`: an argument order copied
from an older document — `( normJ, normB, angle, T )` — compiles and silently evaluates the
wrong point.

### 5. The thermal mass matrix uses the undeformed-mesh density

> Use `density( gTroom )` or `ref_density()` in the mass matrix, never `density( T )`.

BELFEM computes on the undeformed mesh, and `density( T ) = ρ_ref / l(T)³` already accounts for
the expansion the mesh does not undergo. `calculator::MaxwellData` caches the reference value
once at construction.

### 6. A metal's RRR is fixed before use; an alloy's RRR is fixed at construction

> `Metal::set_RRR` sets ρ₀ and, when tables are enabled, builds the lookup table. Call it once,
> before the material is assigned to a block. `Alloy` accepts its RRR only through the factory
> and raises an error on `set_RRR`.
>
> `set_RRR` also rebuilds the λ spline. Because `set_spline` resets that property's dependency
> flags to T-only, `set_RRR` then re-registers the Kohler flags when `depends( rho, normB )` is
> true. This preserves `depends( lambda, normB )` for every pure metal, so the tool and the
> calculator select `lambda( T, B, beta )` on that flag alone.

Most material constants default to NaN — the exceptions set in the constructor are `mu`
(to μ₀), `q` (to 1) and `density_correction` (`cl_Material.cpp:164-170`).
So a metal without an RRR has no ρ(T) and no λ(T). Nothing asserts that; the NaN propagates.

---

## Lookup-Table Contracts

### 7. The table is a cache, built on rank 0, in lockstep on all ranks

> `Metal::populate_rho_database` runs when `set_RRR` finds `depends( rho, angleBxJ )` and
> tables are enabled. Every rank enters it. Rank 0 evaluates every node of the tensor work mesh;
> the other ranks hold an element-less copy and join the projector's collective solve
> (`Database( Mesh *, … , aProject = true )`, `cl_Database.cpp:57`).
> Later constructions load `<label>_RRR<n>.hdf5` from the run directory after a rank-0 format
> probe is broadcast (`cl_Material_Metal.cpp:891`).
>
> Do **not** "fix" this into a distributed build.

The projector's system is small and every solver it can use assembles on the main process
anyway; distributing the build would add more communication than it saves. A consequence for
anyone reading the packing code: `populate_rho_database( Mesh *, Vector< real > & )` packs
values densely over nodes it owns, while the projector reads them back positionally
(`F( tElement->node( i )->index() )`, `cl_DatabaseProjector.cpp:260`).
These are one index space: `Mesh` assigns `node->set_index( tCount++ )` walking the same
container, so under single ownership slot *k* is node *k*. The `owner() != rank` filter is a
permanent no-op under this contract, not a hook for a distributed future.

### 8. Who needs the table

> `Metal::rho( T, B, beta )` and `lambda( T, B, beta )` evaluate Kohler's rule directly without
> a table and read the table when one exists. `Alloy::rho( T, B, beta )` reads the table and
> asserts its presence (`cl_Material_Alloy.hpp:147`) — an
> alloy queried for field dependence must have been constructed with tables enabled.

Two consequences: a `Metal` built with `aBuildTables = false` is slower per query but complete;
an `Alloy` built that way is field-blind, and the assertion is compiled out in release builds.

### 9. Building the table is a construction-time cost

> The build samples a three-dimensional tensor mesh over (T, log₁₀B, β) and solves a
> projection. Do it once per (label, RRR) — the cache file exists so that it happens once per
> run directory, not once per run.

No wall-clock figures are recorded in the tree; an earlier version of this document quoted
timings that had no source. Measure before planning around a number.

---

## Availability and Dependencies

### 10. Check `have()` before querying a property that may be absent

> `have( MaterialProperty )` is the only availability check. Querying an unset property
> dispatches through a null routing pointer.

Built-in materials set everything they list; user-defined materials have exactly what their
init function sets; `Magnesia` has no resistivity; `HastelloyC276` has no field dependence.

### 11. `depends()` is declarative, not protective

> `depends( property, dependency )` reports what the property varies with. It validates
> nothing.

The FEM calculator uses `depends()` to choose between the `( T )` and `( T, B, beta )`
overloads; a caller that ignores it and passes fewer arguments than the material consumes gets
the field-free value without warning. Declare dependencies truthfully — for a user-defined J_c
the routing trusts them completely (usage guide §7.3).

---

## Immutability and Independence

### 12. Materials are immutable after assignment

> After a material has been handed to a mesh block or a kernel, call no `set_*()` on it.

Materials are shared across blocks. Kernels cache material-dependent quantities, and a
mid-simulation `set_RRR` would also try to rebuild the lookup table collectively.

### 13. Materials do not depend on the mesh

> A property is a function of state (T, B, β, J), never of element size, shape or orientation.
> Geometry belongs in the kernels.

---

## Conventions

### 14. Polynomial coefficients are descending

> `set_user_defined_polynomial` and every `polyval` in the module take coefficients from highest
> power to constant term (MATLAB order); the overloads take `Cell< real >` or `std::vector< real >`.

```cpp
// lambda( T ) = 0.01 T + 50
mat->set_user_defined_polynomial( MaterialProperty::lambda, std::vector< real >{ 0.01, 50.0 } );
```

### 15. The label `buffer` is reserved

> The factory rejects a user-defined (`usermat`) material registered under `buffer`; the
> label names the thin-shell φ-formulation buffer layer and resolves to `Magnesia`. A built-in
> section with that label is not checked and also resolves to `Magnesia`.

### 16. Thread safety is a design intent, not a verified contract

> Construct and configure on the master thread. Read-only queries after construction are
> `const` and touch no mutable state, and are expected to be safe from OpenMP threads — but
> BELFEM's stated policy is that nothing internal is thread-safe (`doc/coding_philosophy.md`),
> and no test exercises concurrent material queries. Protect with `#pragma omp critical` when in
> doubt.

---

## Summary

1. Never delete a J_c, n or B-H object after handing it to a material; delete the material.
2. The factory owns nothing; the caller deletes materials, last of all.
3. HTS: `rho_powerlaw( normJ, T, normB, angleNxB )` for the superconducting state.
4. Mass matrix: `density( gTroom )`, never `density( T )`.
5. RRR before use — through the factory for alloys; a metal without one has NaN resistivity.
6. The lookup table builds on rank 0 in lockstep, is cached as `<label>_RRR<n>.hdf5`, and is not to be distributed.
7. `Metal` evaluates Kohler without a table; `Alloy` needs one.
8. `have()` before querying; `depends()` describes, it does not guard.
9. No `set_*()` after assignment to a block.
10. Descending polynomial coefficients; `buffer` is reserved; thread safety is an intent.
