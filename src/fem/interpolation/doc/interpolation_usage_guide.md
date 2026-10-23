# Interpolation Module Usage Guide {#fem_interpolation_interpolation_usage_guide}

**Date:** 2026-01-20
**Module:** src/fem/interpolation
**Purpose:** Comprehensive guide to BELFEM's shape function and integration system for finite element analysis

**Revision History:**

| Date | Version | Changes |
|------|---------|---------|
| 2026-01-20 | 1.0 | Initial comprehensive guide |

---

## Table of Contents

1. [Common Pitfalls](#common-pitfalls)
2. [Overview](#overview)
3. [Architecture](#architecture)
4. [Mental Model: The Four Pillars](#mental-model-the-four-pillars)
5. [Lifetime & Mutation Contract](#lifetime--mutation-contract)
6. [Usage Examples](#usage-examples)
7. [Integration Data (The Main Workhorse)](#integration-data-the-main-workhorse)
8. [Nédélec Edge Elements](#nédélec-edge-elements)
9. [Factory Patterns and Caching](#factory-patterns-and-caching)
10. [Element Support Matrix](#element-support-matrix)
11. [Integration Point System](#integration-point-system)
12. [Derivative Matrix Conventions](#derivative-matrix-conventions)
13. [Performance Considerations](#performance-considerations)
14. [Thread Safety and MPI](#thread-safety-and-mpi)
15. [Common Patterns](#common-patterns)
16. [Development Notes](#development-notes)

---

## Common Pitfalls

**Place pitfalls first** - these are the most common mistakes that lead to performance degradation or incorrect results.

### 1. Creating Factories/IntegrationData Inside Element Loops

```cpp
// WRONG - Performance disaster!
for (Element* e : elements) {
    InterpolationFunctionFactory factory;                    // Factory per element!
    InterpolationFunction* shape = factory.create_lagrange_function(e->type());

    IntegrationData data(e->type());  // Allocation + precomputation per element!
    data.populate(4);

    for (uint k = 0; k < data.number_of_integration_points(); ++k) {
        // ... assembly
    }
    delete shape;
}
```

**Why wrong:** Each factory creation, shape function allocation, and IntegrationData precomputation is expensive, and none of it depends on the element, only on its type. In this loop it is repeated for every element.

```cpp
// CORRECT - Cache by element type ( aElements: Cell< mesh::Element * > & )
Map< ElementType, IntegrationData * > tCache ;

for ( mesh::Element * tElement : aElements )
{
    const ElementType tType = tElement->type() ;

    // Create once per element type
    if ( ! tCache.key_exists( tType ) )
    {
        IntegrationData * tNew = new IntegrationData( tType );
        tNew->populate( 4, IntegrationScheme::GAUSS );
        tCache[ tType ] = tNew ;
    }

    const IntegrationData * tData = tCache( tType );
    for ( uint k = 0; k < tData->number_of_integration_points(); ++k )
    {
        // Reuse precomputed N, dNdXi
    }
}

// Cleanup
for ( auto & tEntry : tCache )
{
    delete tEntry.second ;
}
tCache.clear() ;
```

`Map` (`containers/cl_Map.hpp`) tests for a key with `key_exists()`. `operator[]` inserts a missing key; `operator()` looks a key up and aborts with `BELFEM_ERROR` if it is missing.

---

### 2. Calling d2NdXi2 on an Unspecialized Template

```cpp
InterpolationFunction* shape = factory.create_lagrange_function(ElementType::TET10);

Matrix<real> d2NdXi2;
shape->d2NdXi2(xi, d2NdXi2);  // fine: every factory-created element implements it
```

**Second derivative support:** every element the factory can create implements `d2NdXi2`; the `BELFEM_ERROR` in `InterpolationFunctionTemplate` is reached only for a template combination without a specialization, which the factory never instantiates.

---

### 3. Ignoring Facet Orientation for Slave Sides

```cpp
// WRONG - Missing orientation
IntegrationData tSlaveData( ElementType::TET4 );
tSlaveData.populate_for_slave(
    tSlaveIndex,   // Which facet
    0,             // Orientation 0 assumed for every facet!
    tOrder,
    tScheme );
```

**Problem:** In contact, DG, and interface problems, the slave integration points must follow the master facet's node ordering. If orientation is ignored, the slave points can be placed at the wrong physical locations, corrupting interface fluxes.

```cpp
// CORRECT - take the orientation from the mesh facet ( aFacet: mesh::Facet * )
// Facet::orientation_on_slave() is 1-based; populate_for_slave() takes a 0-based index
tSlaveData.populate_for_slave( aFacet->index_on_slave(),
                               aFacet->orientation_on_slave() - 1,
                               tOrder,
                               tScheme );
```

`Mesh::finalize()` computes slave orientations for all facets by calling `Mesh::compute_facet_orientations()`. This step is skipped if `set_compute_facet_orientation_flag( false )` was called. In the FEM kernel, `fem::SideSet` creates one `IntegrationData` object for each slave facet and orientation. `fem::Calculator` then selects the table using `orientation_on_slave() - 1` (`fem/kernel/cl_FEM_SideSet.cpp`, `fem/kernel/cl_FEM_Calculator.cpp`, `slave_integration_index_tet`).

For a 2D slave element (TRI, and QUAD other than `QUAD4TS`/`QUAD9TS`), `populate_for_slave()` ignores `aOrientation`. A line facet of a volume TRI/QUAD has one orientation (`mesh::number_of_orientations` returns 1), and the slave points are generated by reversing the facet's point order.

**When orientation matters:**
- Contact mechanics (slave-master pairing)
- Discontinuous Galerkin (interior facet integration)
- Interface coupling (h-φ formulation coupling, thin-shell interfaces)

---

### 4. Confusing Auto Integration Order

```cpp
IntegrationData data(ElementType::TRI6);
data.populate(0, IntegrationScheme::GAUSS);  // Order = 0 means "auto"
```

**Problem:** Users assume `0` means "minimal" or "error". It actually triggers automatic order selection based on element type.

**Auto order selection:** BELFEM does not treat `0` as a minimal rule. It maps the element's interpolation order to a fixed integration order through `auto_integration_order()` (`numerics/integration/fn_intpoints_auto_integration_order.cpp`):

| Element interpolation order | Integration order |
|---|---|
| constant | 1 |
| linear | 4 |
| quadratic, serendipity | 7 |
| cubic | 10 |
| quartic | 13 |
| quintic | 16 |

These values are above `2p + 1`, the degree of a mass-like product with a linearly varying coefficient on an affine element. The source does not document why these particular values were chosen.

This is a default, not an exactness guarantee. For distorted QUAD/HEX elements and for curved higher-order mappings, the Jacobian makes the integrand non-polynomial, so exact integration is not available anyway. Such integrands, and nonlinear ones, can still be under-integrated by the default.

The integer order does carry one guarantee: for `IntegrationScheme::GAUSS`, `intpoints()` returns a rule that integrates every polynomial of total degree up to that order exactly on the reference element, for every geometry it serves (locked by `tests/fem/test_IntegrationExactness.cpp`). Several rules deliver more than requested (the tensor-product rules, and the tetrahedron tables for orders 7 and 9), so the number is a floor, not the exact degree. It says nothing about mapped, non-polynomial integrands.

**Recommendation:** Explicitly specify order for clarity:
```cpp
data.populate(7, IntegrationScheme::GAUSS);  // Explicit order for this TRI6 example
```

---

### 5. Ownership Confusion with IntegrationData

```cpp
InterpolationFunction* shape = factory.create_lagrange_function(ElementType::TET4);

// Constructor 2, aClaimOwnership = false: IntegrationData does NOT own shape
IntegrationData data1(ElementType::TET4, shape, false);
// ... use data1 ...
delete shape;  // Caller must delete!

// Constructor 2, aClaimOwnership = true: IntegrationData OWNS shape2
InterpolationFunction* shape2 = factory.create_lagrange_function(ElementType::TET4);
IntegrationData data2(ElementType::TET4, shape2, true);
// ... use data2 ...
// shape2 is deleted in ~IntegrationData()
```

**Problem:** Memory leaks if caller doesn't understand ownership semantics.

**Safe patterns:**
```cpp
// Pattern A: Let IntegrationData create and own (RECOMMENDED)
IntegrationData data(ElementType::TET4);  // Auto-creates Lagrange shape, owns it
data.populate(4);

// Pattern B: Explicit ownership transfer
InterpolationFunction* shape = factory.create_lagrange_function(ElementType::TET4);
IntegrationData data(ElementType::TET4, shape, true);  // Claim ownership
data.populate(4);
// shape is now owned by data - do NOT delete manually
```

---

### 6. Using Linear Bernstein and Expecting Different Results

```cpp
InterpolationFunction* bern = factory.create_bernstein_function(ElementType::TRI3);
// This actually returns a LAGRANGE function internally!
```

**Issue:** For linear elements (TRI3, QUAD4, LINE2), Bernstein polynomials are mathematically identical to Lagrange polynomials. The factory returns `InterpolationType::LAGRANGE` for these cases.

**Impact:** Low (mathematically correct), but can confuse users expecting separate Bernstein implementation.

**Documentation:** Explicitly state: "Linear Bernstein = Linear Lagrange (factory returns Lagrange for efficiency)."

---

### 7. Reading `det_J()` Before the Shape Function Has Updated It

```cpp
real detJ = edge->det_J();  // stale for some elements: no accessor has run at point k yet
const Matrix<real>& E = edge->E(k);
const Matrix<real>& C = edge->C(k);
```

**Problem:** `det_J()` returns a stored value. Which call updates it depends on the element:

| Element | `det_J()` is set by |
|---|---|
| `EF_TRI3`, `EF_TET4`, `EF_QUAD4TS`, `EF_PENTA6TS`, `EF_HEX8TB` | `link()` (one constant value per element) |
| `EF_TRI6`, `EF_TET10` with straight edges | `link()` (one constant value per element) |
| `EF_TRI6`, `EF_TET10` curved | `E(k)`, `C(k)` and `G(k)`, for point `k` |
| `EF_HEX8`, `EF_HEX8TS` | `E(k)`, `C(k)` and `G(k)`, for point `k` (`update_nabla()`) |

For the last two rows, reading `det_J()`, or calling a helper such as
`dV_ts()` that reads it, before any accessor has run at point `k` returns the
determinant from the previous integration point.

**Safe pattern:**
```cpp
for (uint k = 0; k < npts; ++k) {
    const Matrix<real>& E = edge->E(k);
    const Matrix<real>& C = edge->C(k);  // Refreshes detJ
    real detJ = edge->det_J();

    // Safe: any dV helper that reads edge->det_J() now sees the current point
}
```

**Scope:** This affects the elements that refresh `det_J()` per point: curved
`EF_TRI6` and `EF_TET10`, `EF_HEX8` and `EF_HEX8TS`. If a new element caches
Jacobian data lazily, document which accessor performs the update.

---

## Overview

The **interpolation module** provides the core machinery for evaluating shape functions, their derivatives, and integration points in finite element analysis. It is designed for **zero-overhead abstraction** with extensive precomputation capabilities.

### Purpose

- **Shape function evaluation**: N(ξ), ∂N/∂ξ, ∂²N/∂ξ² for all standard element geometries
- **Integration point management**: Gauss quadrature for volume and facet integration
- **Precomputation**: `IntegrationData` stores shape values/derivatives at integration points to avoid redundant computation in assembly loops
- **Nédélec elements**: Edge-based (H(curl)) shape functions for electromagnetic FEM

### Key Features

1. **Template-based architecture**: `InterpolationFunctionTemplate<Geometry, Type, Dimension, Bases>` provides compile-time specialization
2. **Factory pattern**: `InterpolationFunctionFactory` and `EdgeFunctionFactory` create shape functions by element type
3. **Multiple interpolation families**: Lagrange, Hermite, Bernstein, Bubble functions
4. **Integration schemes**: Gauss quadrature with auto-order selection
5. **Facet integration**: Master/slave facet handling with orientation for contact/DG
6. **Performance optimization**: Precomputed IntegrationData is the **recommended** usage pattern

### Design Philosophy

Following `doc/coding_philosophy.md`:
- **Manual memory management**: Factories return raw pointers; caller or IntegrationData owns
- **Zero abstraction penalty**: Release builds compile to direct array access (no virtual function overhead after precompute)
- **Explicit ownership**: IntegrationData constructor specifies ownership (`aClaimOwnership` flag)
- **Caching-first**: IntegrationData is designed to be cached and reused across elements of the same type

---

## Architecture

### File Organization

```
src/fem/interpolation/
├── cl_IF_InterpolationFunction.hpp        # Abstract base class
├── cl_IF_InterpolationFunctionTemplate.hpp # Template implementation
├── cl_IF_InterpolationFunctionFactory.{hpp,cpp} # Factory for shape functions
├── cl_IF_IntegrationData.{hpp,cpp}        # Precomputed integration bundle
├── cl_EdgeFunctionFactory.{hpp,cpp}       # Factory for Nédélec elements
├── fn_IF_initialize_integration_points.*  # Integration point generation
├── fn_IF_initialize_integration_points_on_facet.* # Facet integration points
├── fn_IF_initialize_shape_function.*      # Shape function initialization
├── lagrange/                              # Lagrange polynomial specializations
├── hermite/                               # Hermite polynomial specializations
├── bernstein/                             # Bernstein polynomial specializations
├── bubble/                                # Bubble function stabilizations
└── nedelec/                               # Nédélec edge element implementations
    ├── cl_EF_TRI3.{hpp,cpp}               # 1st order triangle edge element
    ├── cl_EF_TRI6.{hpp,cpp}               # 2nd order triangle edge element
    ├── cl_EF_TET4.{hpp,cpp}               # 1st order tetrahedron edge element
    ├── cl_EF_TET10.{hpp,cpp}              # 2nd order tetrahedron edge element
    ├── cl_EF_HEX8.{hpp,cpp}               # 1st order hexahedron edge element
    ├── cl_EF_QUAD4TS.{hpp,cpp}            # Thin shell with curl (quad)
    ├── cl_EF_PENTA6TS.{hpp,cpp}           # Thin shell with curl (tri)
    ├── cl_EF_HEX8TS.{hpp,cpp}             # Thin shell (hex)
    ├── cl_EF_HEX8TB.{hpp,cpp}             # Side-connector wall element, 4 DOFs
    ├── cl_EF_LINE3.{hpp,cpp}              # Proof of concept, not dispatched by the factory
    └── fn_num_nedelec_dofs.hpp            # DOF count per element type
```

### Class Hierarchy

```
InterpolationFunction (abstract base)
    ├── InterpolationFunctionTemplate<G,T,D,B>
    │   ├── Lagrange specializations (LINE2-5, TRI3-15, QUAD4-16, TET4-35, ...)
    │   ├── Hermite specializations (LINE2, QUAD4)
    │   ├── Bernstein specializations (LINE, TRI)
    │   └── Bubble specializations (TRI3/6 edges, TET4/10 faces)
    └── (Direct derivations for special cases)

EdgeFunction (abstract base for Nédélec)
    ├── EF_TRI3, EF_TRI6    # 2D electromagnetics
    ├── EF_TET4, EF_TET10, EF_HEX8   # 3D electromagnetics
    ├── EF_QUAD4TS, EF_PENTA6TS, EF_HEX8TS, EF_HEX8TB # Thin shells
    └── EF_LINE3            # Proof of concept (not created by the factory)

IntegrationData (precompute bundle)
    └── Owns or borrows an InterpolationFunction

InterpolationFunctionFactory
    ├── create_lagrange_function(ElementType)
    ├── create_hermite_function(ElementType)
    ├── create_bernstein_function(ElementType)
    └── create_bubble_function(ElementType, uint aFacet)

EdgeFunctionFactory
    └── create_edge_function(ElementType)
```

---

## Mental Model: The Four Pillars

Understanding the interpolation module requires grasping four distinct but interconnected concepts:

### Pillar 1: InterpolationFunction (Shape Function Interface)

**What it is:** Abstract base class defining the contract for evaluating shape functions.

**Key methods:**
```cpp
virtual void N(const Vector<real> & aXi, Matrix<real> & aN) const = 0;
virtual void dNdXi(const Vector<real> & aXi, Matrix<real> & adNdXi) const = 0;
virtual void d2NdXi2(const Vector<real> & aXi, Matrix<real> & ad2NdXi2) const = 0;
virtual void param_coords(Matrix<real> & aXiHat) const = 0;
```

**Matrix dimension contracts:**
- `aN`: (1 × num_bases) - shape function values
- `adNdXi`: (num_dimensions × num_bases) - first derivatives
- `ad2NdXi2`: See [Derivative Matrix Conventions](#derivative-matrix-conventions)
- `aXiHat`: (num_dimensions × num_nodes) - parametric node coordinates. For Lagrange elements the number of nodes equals the number of bases; for Hermite it does not (`LINE2`: 1 × 2 for 4 bases, `QUAD4`: 2 × 4 for 16 bases)

**When to use directly:** Rarely. Most users should use `IntegrationData` (Pillar 3).

**When to derive custom classes:** Implementing non-standard shape functions (e.g., NURBS, hierarchical p-refinement).

---

### Pillar 2: InterpolationFunctionFactory (Creation)

**What it is:** Factory for creating shape function instances.

**Usage pattern:**
```cpp
InterpolationFunctionFactory factory;
InterpolationFunction* shape = factory.create_lagrange_function(ElementType::TET10);

// Use shape...

delete shape;  // Caller owns the pointer
```

**Methods:**
- `create_lagrange_function(ElementType)` - Standard Lagrange polynomials
- `create_hermite_function(ElementType)` - C¹ Hermite polynomials (beams, plates)
- `create_bernstein_function(ElementType)` - Bézier-like basis
- `create_bubble_function(ElementType, uint aFacet)` - Stabilization bubbles (the second argument is a facet index, not a `BubbleType`)

**Ownership:** Factory returns **raw pointer**; caller must delete or transfer to IntegrationData.

---

### Pillar 3: IntegrationData (The Main Workhorse)

**What it is:** Precomputed bundle of integration points, weights, and shape function values/derivatives.

**Why it matters:** This is the **recommended** way to use the interpolation module. It:
1. Generates integration points (Gauss quadrature)
2. Evaluates `N`, `dNdXi`, `d2NdXi2` **once** at each point
3. Stores results for **O(1) access** in assembly loops

**Usage pattern:**
```cpp
IntegrationData data(ElementType::TET10);
data.populate(4, IntegrationScheme::GAUSS);  // Order 4 Gauss

for (uint k = 0; k < data.number_of_integration_points(); ++k) {
    const Matrix<real>& N = data.N(k);      // Precomputed!
    const Matrix<real>& dN = data.dNdXi(k); // Precomputed!
    real w = data.weights()(k);

    // Element assembly: Ke += N^T * D * N * w * detJ
}
```

**Critical for performance:** See [Factory Patterns and Caching](#factory-patterns-and-caching).

---

### Pillar 4: EdgeFunction (Nédélec for Electromagnetics)

**What it is:** Separate branch for vector-valued (H(curl)) shape functions used in electromagnetic FEM (h-φ formulation, Maxwell solvers).

**Key difference from scalar shapes:**
- Vector-valued: `E(ξ)` returns `(ndim × ndofs)` matrix
- Curl operator: `C(k)` returns curl matrix
- Requires linking to an FEM element: `link( fem::Element * )` reads node coordinates and edge directions

**Usage pattern:**
```cpp
EdgeFunctionFactory edgeFactory;
EdgeFunction* edge = edgeFactory.create_edge_function(ElementType::TET4);

edge->precompute(integrationPoints);  // once per set of integration points
edge->link(femElement);               // per element ( fem::Element * )

const Matrix<real>& E = edge->E(integrationPointIndex);  // Interpolation (3 × 6 for TET4)
const Matrix<real>& C = edge->C(integrationPointIndex);  // Curl (3 × 6)
real detJ = edge->det_J();

delete edge;
```

Call `precompute()` first. For `EF_TRI6` and `EF_TET10` with straight edges, `link()` evaluates the Jacobian from the precomputed data. `fem::Calculator` follows the same sequence: it precomputes when it creates the edge function and links once per element (`fem/kernel/cl_FEM_Calculator.cpp`).

**Prerequisites:**
- Mesh must have edges/faces created (`Mesh::create_edges()`, `Mesh::create_faces()`)
- Element must provide edge orientation data
- Edge DOFs are managed by `DofManager` (`fem/kernel/cl_FEM_DofManager.hpp`)

**See:** [Nédélec Edge Elements](#nédélec-edge-elements) for details.

---

## Lifetime & Mutation Contract

Understanding when objects can be mutated and when they should be treated as immutable is critical for safe, performant code.

### InterpolationFunction Lifetime & Mutability

**Lifetime:**
- Created by factory → Owned by caller or IntegrationData
- Must live until all dependent `IntegrationData` objects are destroyed (if borrowed)
- Safe to delete after `IntegrationData` with `aClaimOwnership = true` takes ownership

**Mutability:**
```cpp
InterpolationFunction* shape = factory.create_lagrange_function(ElementType::TET10);

// ✅ IMMUTABLE after creation - safe to call from multiple threads
Vector<real> xi(3);
Matrix<real> N, dNdXi;
shape->N(xi, N);        // Read-only, thread-safe
shape->dNdXi(xi, dNdXi); // Read-only, thread-safe

// ❌ NO mutation methods exist - object is logically const after construction
```

**Contract:** `InterpolationFunction` is **effectively immutable** after factory creation. All methods are `const` and thread-safe.

---

### IntegrationData Lifetime & Mutability

**Lifetime:**
- Lives from construction to explicit `delete` (manual memory management)
- If owns shape function (`aClaimOwnership = true`), deletes it in destructor
- If borrows shape function (`aClaimOwnership = false`), caller must ensure shape outlives `IntegrationData`

**Mutability - Two phases:**

#### Phase 1: Mutable (Setup)
```cpp
IntegrationData data(ElementType::TET10);

// ✅ MUTABLE - can call populate() multiple times
data.populate(2, IntegrationScheme::GAUSS);  // First populate
data.populate(4, IntegrationScheme::GAUSS);  // Re-populate with different order (overwrites)

// OR populate for facet
data.populate_for_master(2, 3);  // Overwrites volume integration
```

**Warning:** Each `populate()` call **overwrites** previous integration data. Don't assume cached values survive re-population.

#### Phase 2: Immutable (After populate)
```cpp
data.populate(4);

// ✅ IMMUTABLE - all access methods are const, thread-safe
for (uint k = 0; k < data.number_of_integration_points(); ++k) {
    const Matrix<real>& N = data.N(k);      // Read-only
    const Matrix<real>& dN = data.dNdXi(k); // Read-only
    real w = data.weights()(k);             // Read-only
}
```

**Contract:** After final `populate()` call, treat `IntegrationData` as **immutable** for performance and thread safety.

---

### Ownership Patterns

#### Pattern 1: IntegrationData Owns Shape (Recommended)
```cpp
IntegrationData* data = new IntegrationData(ElementType::TET10);
data->populate(4);

// ... use data ...

delete data;  // Deletes owned shape function automatically
```

**Lifetime:** `data` owns shape → both deleted together.

---

#### Pattern 2: IntegrationData Borrows Shape
```cpp
InterpolationFunctionFactory factory;
InterpolationFunction* shape = factory.create_lagrange_function(ElementType::TET10);

IntegrationData* data = new IntegrationData(ElementType::TET10, shape, false);  // Borrow
data->populate(4);

// ... use data ...

delete data;   // Does NOT delete shape
delete shape;  // Caller must delete
```

**Lifetime:** Caller must ensure `shape` outlives `data`.

---

#### Pattern 3: Shared Shape Across Multiple IntegrationData
```cpp
InterpolationFunction* shape = factory.create_lagrange_function(ElementType::TET10);

// Multiple IntegrationData instances borrow same shape
IntegrationData* volumeData = new IntegrationData(ElementType::TET10, shape, false);
volumeData->populate(4);

IntegrationData* facetData = new IntegrationData(ElementType::TET10, shape, false);
facetData->populate_for_master(2, 3);

// ... use both ...

delete volumeData;  // Does NOT delete shape
delete facetData;   // Does NOT delete shape
delete shape;       // Caller deletes once
```

**Use case:** When you need volume and facet integration for the same element type.

---

### EdgeFunction Lifetime & Mutability

**Lifetime:**
- Created by factory → Owned by caller
- Must be deleted explicitly after use

**Mutability - Two phases:**

#### Phase 1: Mutable (Setup)
```cpp
EdgeFunction* edge = edgeFactory.create_edge_function(ElementType::TET4);

// ✅ MUTABLE - precompute and link modify internal state
edge->precompute(gaussPoints);  // Reference-element data at the integration points
edge->link(element);            // Element context (node coords, edge orientations); fem::Element *
```

**Contract:** `precompute()` and then `link()` **must** be called before accessing `E()`, `C()` or `G()`.

#### Phase 2: Evaluation (After link)
```cpp
// E(), C() and G() are non-const: they compute into member buffers
for (uint k = 0; k < npts; ++k) {
    const Matrix<real>& E = edge->E(k);  // Overwritten by the next E() call
    const Matrix<real>& C = edge->C(k);  // Overwritten by the next C() call
    real detJ = edge->det_J();
}
```

**Warning:** `EdgeFunction` is **not thread-safe** at any stage. `E()`, `C()` and `G()` write to member buffers, and some elements also update `det_J()` (see Common Pitfall 7). A returned reference remains valid only until the next call to the same accessor.

---

### Summary Table

| Object | Creation | Mutable Phase | Immutable Phase | Thread-Safe After? |
|--------|----------|---------------|-----------------|-------------------|
| **InterpolationFunction** | Factory | Never (effectively const) | Always | ✅ Yes |
| **IntegrationData** | Constructor | During `populate()` calls | After final `populate()` | ✅ Yes |
| **EdgeFunction** | Factory | Always (`precompute()`, `link()`, `E()`, `C()`, `G()`) | Never | ❌ No |

---

### Common Lifetime Bugs

**Bug 1: Deleting borrowed shape too early**
```cpp
InterpolationFunction* shape = factory.create_lagrange_function(ElementType::TET10);
IntegrationData data(ElementType::TET10, shape, false);  // Borrow
data.populate(4);

delete shape;  // ❌ BUG! data still references shape

data.N(0);  // ❌ Undefined behavior (dangling pointer)
```

**Fix:** Ensure `shape` outlives `data`.

---

**Bug 2: Re-populating in parallel region**
```cpp
IntegrationData data(ElementType::TET10);
data.populate(4);

#pragma omp parallel for
for (int i = 0; i < n; ++i) {
    data.populate(2);  // ❌ BUG! Race condition
    // ... assembly ...
}
```

**Fix:** Populate **once** before parallel region.

---

**Bug 3: Reusing EdgeFunction without re-linking**
```cpp
EdgeFunction* edge = edgeFactory.create_edge_function(ElementType::TET4);

edge->precompute(gaussPoints);
for (fem::Element* e : elements) {
    // ❌ BUG! Missing link(e)
    // ... assembly (uses the previous element's context) ...
}
```

**Fix:** Precompute once, then call `link(e)` for each element.

```cpp
edge->precompute(gaussPoints);
for (fem::Element* e : elements) {
    edge->link(e);                  // ✅ Correct
    // ... assembly ...
}
```

---

## Usage Examples

### Example 1: Basic Shape Function Evaluation

```cpp
#include "cl_IF_InterpolationFunctionFactory.hpp"

using namespace belfem::fem;

InterpolationFunctionFactory factory;
InterpolationFunction* shape = factory.create_lagrange_function(ElementType::TRI6);

// Parametric coordinates (evaluate at element center)
Vector<real> xi(2);
xi(0) = 1.0/3.0;
xi(1) = 1.0/3.0;

// Shape function values (1 × 6)
Matrix<real> N;
shape->N(xi, N);

// First derivatives (2 × 6)
Matrix<real> dNdXi;
shape->dNdXi(xi, dNdXi);

// Parametric node coordinates (2 × 6)
Matrix<real> xiHat;
shape->param_coords(xiHat);

// Matrix has no operator<<; print entry by entry
for (uint i = 0; i < shape->number_of_bases(); ++i) {
    std::cout << i << ": N = " << N(0, i)
              << ", dN/dxi = " << dNdXi(0, i)
              << ", dN/deta = " << dNdXi(1, i) << std::endl;
}

delete shape;
```

**Output interpretation:**
- `N(0,i)` = shape value at node `i`
- `dNdXi(0,i)` = ∂N_i/∂ξ
- `dNdXi(1,i)` = ∂N_i/∂η

---

### Example 2: Integration Data for Element Assembly (Recommended)

```cpp
#include "cl_IF_IntegrationData.hpp"

IntegrationData data(ElementType::TET10, InterpolationType::LAGRANGE);
data.populate(4, IntegrationScheme::GAUSS);  // Order 4 Gauss quadrature

uint npts = data.number_of_integration_points();
std::cout << "Number of integration points: " << npts << std::endl;

for (uint k = 0; k < npts; ++k) {
    const Matrix<real>& N_k = data.N(k);          // (1 × 10)
    const Matrix<real>& dN_k = data.dNdXi(k);     // (3 × 10)
    real w_k = data.weights()(k);

    // Access parametric coordinates if needed
    Vector<real> xi_k = data.points().col(k);

    // Element stiffness integration:
    // Ke += dN^T * D * dN * w * detJ
    // (Compute detJ from mesh element coordinates and dN_k)
}
```

**Why this is recommended:**
- Shape functions evaluated **once** during `populate()`
- Access via `N(k)`, `dNdXi(k)` in loops returns a stored matrix; nothing is evaluated

---

### Example 3: Hermite Shape Functions for Beam Elements

```cpp
// Hermite LINE2: 4 bases (value + derivative at each node)
InterpolationFunction* hermite = factory.create_hermite_function(ElementType::LINE2);

Vector<real> xi(1);
xi(0) = -1.0;  // Left end ( node 0 sits at xi = -1, node 1 at xi = +1 )

Matrix<real> N;
hermite->N(xi, N);  // N = (1 × 4)

// N = [1, 0, 0, 0] at the left end
// ( at xi = 0 it is [0.5, 0.25, 0.5, -0.25] )
// Used for Euler-Bernoulli beam (cubic interpolation)

delete hermite;
```

**Hermite basis interpretation** (LINE2, `hermite/cl_IF_BEAM.hpp`):
- N(0) = value at node 0
- N(1) = derivative with respect to ξ at node 0
- N(2) = value at node 1
- N(3) = derivative with respect to ξ at node 1

The slope bases are defined in ξ. Therefore, the coefficient of N(1) and N(3) is dw/dξ = (dx/dξ) · dw/dx, which is L/2 times the physical slope for a straight element of length L.

---

### Example 4: Facet Integration for Boundary Conditions

```cpp
IntegrationData boundaryData(ElementType::TET4);

uint masterFacetIndex = 2;  // Which facet of the tet
boundaryData.populate_for_master(masterFacetIndex, 3, IntegrationScheme::GAUSS);

// Integration points are now on facet 2 of TET4 (a triangle)
for (uint k = 0; k < boundaryData.number_of_integration_points(); ++k) {
    const Matrix<real>& N_k = boundaryData.N(k);  // Shape on facet
    real w_k = boundaryData.weights()(k);

    // Apply Neumann BC: ∫_Γ N^T * q * dS
}
```

**Master vs. Slave:**
- **Master facet**: Own element's facet (for boundary conditions)
- **Slave facet**: Neighboring element's facet (for contact/DG, requires orientation)

---

### Example 5: Slave Facet Integration with Orientation (Contact/DG)

```cpp
IntegrationData slaveData(ElementType::HEX8);

uint slaveFacetIndex = 4;     // Which facet of this hex
uint orientation = 3;         // 0-based: facet->orientation_on_slave() - 1, range 0..3 for a HEX facet

slaveData.populate_for_slave(slaveFacetIndex, orientation, 3, IntegrationScheme::GAUSS);

// Integration points now match the master element's facet node ordering
```

**Orientation encoding:** `Facet::compute_orientation()` finds the first corner node of the slave facet in the master facet's corner-node list. Its 1-based position is `Facet::orientation_on_slave()`. `populate_for_slave()` takes that value minus one.

The number of orientations per facet is given by `mesh::number_of_orientations()`: 3 for a TET facet, 4 for a HEX facet, 4 or 3 for the quad or triangle facets of a PENTA, 1 for the line facets of a volume TRI/QUAD, and 2 for the line facets of `QUAD4TS`/`QUAD9TS`.

---

### Example 6: Nédélec Edge Elements for Electromagnetics

```cpp
#include "cl_EdgeFunctionFactory.hpp"

EdgeFunctionFactory edgeFactory;
EdgeFunction* edge = edgeFactory.create_edge_function(ElementType::TET4);

// Precompute at integration points ( once per point set )
IntegrationData gauss(ElementType::TET4);
gauss.populate(4);
edge->precompute(gauss.points());

// Link to the FEM element ( fem::Element *, provides edge orientations )
edge->link(femElement);

// Access interpolation and curl at integration point 0
const Matrix<real>& E = edge->E(0);  // (3 × 6) for TET4 (6 edges)
const Matrix<real>& C = edge->C(0);  // (3 × 6) curl operator

real detJ = edge->det_J();
real absDetJ = edge->abs_det_J();

// Assembly: ∫ (curl E)^T * μ^-1 * (curl E) dV
// Ke += C^T * (μ^-1 / detJ) * C * w * detJ

delete edge;
```

**Nédélec DOF counts:**
- TET4: 6 DOFs (one per edge)
- TET10: 20 DOFs (2 per edge for 2nd order)
- TRI3: 3 DOFs
- TRI6: 8 DOFs

---

## Integration Data (The Main Workhorse)

`IntegrationData` is the **recommended** entry point for most users. It combines integration point generation with precomputed shape function evaluation.

### Constructors

```cpp
// Constructor 1: Auto-create a shape function of type aType
IntegrationData(const ElementType aElementType,
                const InterpolationType aType = InterpolationType::LAGRANGE,
                InterpolationFunction* aShapeFunction = nullptr);

// Constructor 2: Provide shape function with ownership control
IntegrationData(const ElementType aElementType,
                InterpolationFunction* aShapeFunction,
                const bool aClaimOwnership);
```

**Ownership semantics:**
- Constructor 1 with `nullptr`: Creates and **owns** shape function (deleted in destructor)
- Constructor 1 with a non-null `aShapeFunction`: **Borrows** it (caller must delete); `aType` is then not used
- Constructor 2 with `aClaimOwnership = true`: **Owns** provided shape function
- Constructor 2 with `aClaimOwnership = false`: **Borrows** shape function (caller must delete)

**Recommended pattern:** Use constructor 1 (auto-create and own):
```cpp
IntegrationData data(ElementType::TET10);  // Simple, safe
```

---

### Populate Methods

#### Volume Integration

```cpp
void populate(const uint aIntegrationOrder = 0,
              const IntegrationScheme aScheme = IntegrationScheme::GAUSS);
```

**Parameters:**
- `aIntegrationOrder`: Integration order (0 = auto-select based on element type)
- `aScheme`: one of `IntegrationScheme` (`numerics/integration/en_IntegrationScheme.hpp`):
  - `GAUSS`: the default; tabulated rules for every geometry
  - `GAUSSCLASSIC`: tensor-product Gauss-Legendre rules for QUAD and HEX; for every other geometry `intpoints()` falls back to `GAUSS`
  - `LOBATTO`: Gauss-Lobatto, LINE only; other geometries abort with `BELFEM_ERROR`

**Auto order:** When `aIntegrationOrder = 0`, BELFEM calls `auto_integration_order( aElementType )`. It returns 4 / 7 / 10 / 13 / 16 for linear / quadratic (and serendipity) / cubic / quartic / quintic elements (see Common Pitfall 4).

For distorted QUAD/HEX elements and curved higher-order mappings, the Jacobian makes the integrand non-polynomial, so no fixed order integrates it exactly.

**Usage:**
```cpp
IntegrationData data(ElementType::TET10);
data.populate(7);  // Explicit order for this TET10 example
// or
data.populate();   // Auto-select: order 7 for the quadratic TET10
```

---

#### Master Facet Integration (Boundary Conditions)

```cpp
void populate_for_master(const uint aMasterIndex,
                         const uint aIntegrationOrder = 0,
                         const IntegrationScheme aScheme = IntegrationScheme::GAUSS);
```

**Parameters:**
- `aMasterIndex`: Local facet index (0-based, element-dependent)
  - TET4: 0-3 (4 triangular faces)
  - HEX8: 0-5 (6 quadrilateral faces)
  - etc.

**Use case:** Neumann boundary conditions, surface integrals

```cpp
IntegrationData bcData(ElementType::HEX8);
bcData.populate_for_master(2, 3);  // Face 2, order 3

for (uint k = 0; k < bcData.number_of_integration_points(); ++k) {
    const Matrix<real>& N = bcData.N(k);
    real w = bcData.weights()(k);
    // ∫_Γ N^T * traction * dS
}
```

---

#### Slave Facet Integration (Contact/DG)

```cpp
void populate_for_slave(const uint aSlaveIndex,
                        const uint aOrientation = 0,
                        const uint aIntegrationOrder = 0,
                        const IntegrationScheme aScheme = IntegrationScheme::GAUSS);
```

**Parameters:**
- `aSlaveIndex`: Local facet index on slave element
- `aOrientation`: 0-based orientation index, `Facet::orientation_on_slave() - 1`
  - Range `0 .. mesh::number_of_orientations( type, facet ) - 1` (see Example 5)
  - Ignored for TRI and for QUAD other than `QUAD4TS`/`QUAD9TS` (one orientation per line facet)

**Use case:** Contact mechanics, discontinuous Galerkin, interface coupling (h-φ, thin-shell)

```cpp
// aFacet: mesh::Facet * with master and slave; orientations are set by Mesh::finalize()
slaveData.populate_for_slave(aFacet->index_on_slave(),
                             aFacet->orientation_on_slave() - 1,
                             3);
```

**Orientation example (TRI facet on TET):** if the slave facet's first corner node is at position 2 (1-based) in the master facet's corner list, then `orientation_on_slave()` is 2 and `aOrientation` is 1.

---

### Access Methods

```cpp
uint number_of_integration_points() const;
const Vector<real>& weights() const;      // (npts,)
const Matrix<real>& points() const;       // parametric coords, one column per point
                                          // (ndim × npts; volume TRI/TET rules carry ndim+1 rows)

const Vector<real>& phi(uint k) const;    // (nbases,) as vector
const Matrix<real>& N(uint k) const;      // (1 × nbases) as matrix
const Matrix<real>& Nvector(uint k) const; // For vector fields (ndim × ndim*nbases)
const Matrix<real>& dNdXi(uint k) const;  // (ndim × nbases)
const Matrix<real>& d2NdXi2(uint k) const; // See [Derivative Matrix Conventions]

const Vector<real>& dphidxi(uint k) const; // LINE elements only
```

**Common usage:**
```cpp
for (uint k = 0; k < data.number_of_integration_points(); ++k) {
    const Matrix<real>& N = data.N(k);
    const Matrix<real>& dN = data.dNdXi(k);
    real w = data.weights()(k);

    // Stiffness: Ke += (B^T * D * B) * w * detJ
    // where B = dN * inv(J), J = Jacobian from mesh
}
```

---

## Nédélec Edge Elements

Nédélec (or Whitney) elements are vector-valued shape functions for H(curl)-conforming FEM, primarily used in electromagnetic applications (Maxwell solvers, h-φ formulation).

### What Makes Edge Elements Different

| Property | Scalar (Lagrange) | Edge (Nédélec) |
|----------|-------------------|----------------|
| **DOF location** | Nodes | Edges (tangential component) |
| **Continuity** | C⁰ (value) | Tangential component continuous |
| **Output** | Scalar N | Vector E (tangent-aligned) |
| **Operator** | Gradient (grad N) | Curl (curl E) |
| **FEM space** | H¹ (continuous) | H(curl) (tangential continuity) |

### Supported Elements

| Element | DOFs | Order | Implementation |
|---------|------|-------|----------------|
| **TRI3** | 3 | 1st | `EF_TRI3` |
| **TRI6** | 8 | 2nd | `EF_TRI6` |
| **TET4** | 6 | 1st | `EF_TET4` |
| **TET10** | 20 | 2nd | `EF_TET10` |
| **LINE3** | 2 | 2nd | `EF_LINE3` — ⚠️ proof of concept, and **not dispatched by `EdgeFunctionFactory`** (no `LINE3` case); see `nedelec.md` §6.5 |
| **HEX8** | 12 | 1st | `EF_HEX8` |
| **QUAD4TS** | 2 | Thin shell | `EF_QUAD4TS` (line-in-2D extrusion) |
| **PENTA6TS** | 6 | Thin shell | `EF_PENTA6TS` (triangle-in-3D extrusion) |
| **HEX8TS** | 8 | Thin shell | `EF_HEX8TS` (quad-in-3D extrusion) |
| **HEX8TB** | 4 | Thin wall | `EF_HEX8TB` (side-connector wall element, longitudinal edges only) |
| **QUAD9TS**, **PENTA18TS** | 3, 16 (`num_nedelec_dofs`) | Thin shell | none: `EdgeFunctionFactory` aborts with `BELFEM_ERROR` (quadratic thin-shell edge functions are not implemented) |

**DOF count formula:**
- **1st order**: `num_edges` (one DOF per edge)
  - TRI3: 3 DOFs (3 edges)
  - TET4: 6 DOFs (6 edges)
- **2nd order**: `2 × num_edges + num_faces × face_dofs`
  - TRI6: 8 DOFs (3 edges × 2 = 6 edge DOFs + 1 face × 2 = 2 face DOFs)
  - TET10: 20 DOFs (6 edges × 2 = 12 edge DOFs + 4 faces × 2 = 8 face DOFs)
- **Thin-shell elements** use reduced DOF counts from dimensional reduction — see [nedelec_thinshell.md](nedelec_thinshell.md) for the construction.

For a deeper treatment of the Nédélec theory in BELFEM, see:
- [nedelec.md](nedelec.md) — general framework, `EdgeFunction` base class, volume elements, unit-circulation convention
- [nedelec_thinshell.md](nedelec_thinshell.md) — thin-shell (`QUAD4TS`, `PENTA6TS`) elements

---

### EdgeFunction API

```cpp
namespace belfem::fem {
class EdgeFunction {
public:
    virtual void link(Element* aElement) = 0;  // fem::Element: node coords, edge directions
    virtual void precompute(const Matrix<real>& aXi) = 0;  // Precompute at points

    virtual const Matrix<real>& E(const uint aIndex);          // Interpolation operator
    virtual const Matrix<real>& C(const uint aIndex = 0) = 0;  // Curl operator
    virtual const Matrix<real>& G(const uint aIndex = 0) = 0;  // Gradient operator

    real det_J() const;       // Jacobian determinant
    real abs_det_J() const;   // |det(J)|
    real sum_w() const;       // Sum of the integration weights
    uint ndofs() const;       // Number of DOFs

    virtual void update_nabla(const uint aIndex);  // Base version aborts; overridden by EF_HEX8, EF_HEX8TS
};
}
```

The base-class `E()`, `C()` and `update_nabla()` bodies abort with `BELFEM_ERROR`; every element the factory creates overrides `E()` and `C()`.

**Matrix dimensions:**
- `E(k)`: (ndim × ndofs) - interpolates edge DOFs to physical vector field
- `C(k)`: 3D: (3 × ndofs), vector curl; 2D: (1 × ndofs), scalar curl
- `G(k)`: (ndim² × ndofs) - gradient of each basis function, column-major vectorized (layout documented at `mGrad` in `nedelec/cl_EF_EdgeFunction.hpp`)

---

### Usage Pattern

```cpp
EdgeFunctionFactory factory;
EdgeFunction* edge = factory.create_edge_function(ElementType::TET4);

// Precompute at integration points ( once per point set )
IntegrationData gauss(ElementType::TET4);
gauss.populate(4);
edge->precompute(gauss.points());

// CRITICAL: Link to the FEM element ( fem::Element * ) before evaluating
edge->link(femElement);  // Provides edge orientations, node coordinates

// Assembly loop
for (uint k = 0; k < gauss.number_of_integration_points(); ++k) {
    const Matrix<real>& E = edge->E(k);  // (3 × 6) for TET4
    const Matrix<real>& C = edge->C(k);  // (3 × 6)
    real detJ = edge->det_J();
    real w = gauss.weights()(k);

    // Weak form: ∫ (curl E)^T * μ^-1 * (curl E) dV
    // Ke += C^T * (μ^-1 / detJ) * C * w * detJ
}

delete edge;
```

---

### Prerequisites for Edge Elements

**Mesh must have edges created:**
```cpp
mesh->create_edges();  // Or load from file with edges
```

**Element must provide:**
- Edge connectivity (which nodes form each edge)
- Edge orientation (tangent direction convention)
- Node coordinates (for Jacobian computation)

**Typical use cases:**
- h-φ formulation for superconductor magnets (h-field uses edge elements)
- Time-harmonic Maxwell's equations
- Magnetostatics with multiply-connected domains (via cohomology cuts)

---

### Curved vs. Straight Element Optimization

Higher-order Nédélec elements (TRI6, TET10) select their evaluation routines in `link()` through member-function pointers, depending on whether the mesh element is curved. Abridged from `EF_TRI6`:

```cpp
// cl_EF_TRI6.hpp
void ( EF_TRI6::*mFunInterpolation )( const uint aIndex );
void ( EF_TRI6::*mFunCurl )( const uint aIndex );
void ( EF_TRI6::*mFunGrad )( const uint aIndex );

// EF_TRI6::link()
if( tElement->is_curved() )
{
    mFunInterpolation = & EF_TRI6::E_curved ;
    mFunCurl          = & EF_TRI6::C_curved ;   // Jacobian at every point
    mFunGrad          = & EF_TRI6::G_curved ;
}
else
{
    this->compute_nabla( 0 );                   // constant Jacobian, once per element
    mFunInterpolation = & EF_TRI6::compute_E ;
    mFunCurl          = & EF_TRI6::C_straight ;
    mFunGrad          = & EF_TRI6::G_straight ;
}

// cl_EF_TRI6.hpp
inline const Matrix< real > &
EF_TRI6::C( const uint aIndex )
{
    ( this->*mFunCurl )( aIndex );
    return mC ;
}
```

`EF_TET10` uses the same mechanism with the pointers `mFunInterpolation` (`E_curved` / `E_straight`), `mFunDerivatives` (`E_xi_curved` / `E_xi_straight`, used by `C()`) and `mFunGrad` (`G_curved` / `G_straight`).

**Performance impact:** On a straight element the Jacobian is computed once in `link()` rather than at every integration point.

---

## Factory Patterns and Caching

### The Performance Problem

```cpp
// WRONG - Catastrophic performance
for (Element* e : elements) {
    InterpolationFunctionFactory factory;  // ❌ Per element!
    InterpolationFunction* shape = factory.create_lagrange_function(e->type());

    IntegrationData data(e->type());  // ❌ Per element!
    data.populate(4);

    // ... assembly ...
    delete shape;
}
```

**Why catastrophic:**
- 10,000 elements × (factory creation + shape allocation + IntegrationData precompute)
- Each precompute allocates and evaluates `N`, `dNdXi` and `d2NdXi2` at every integration point (11 points for TET10 at order 4)
- Total: about 10,000 × 11 = 110,000 redundant evaluations of each of the three

---

### The Solution: Cache by Element Type

```cpp
// CORRECT - Cache IntegrationData ( aElements: Cell< mesh::Element * > & )
Map< ElementType, IntegrationData * > tCache ;

for ( mesh::Element * tElement : aElements )
{
    const ElementType tType = tElement->type() ;

    // Create once per element type
    if ( ! tCache.key_exists( tType ) )
    {
        IntegrationData * tNew = new IntegrationData( tType );
        tNew->populate( 4, IntegrationScheme::GAUSS );
        tCache[ tType ] = tNew ;
    }

    const IntegrationData * tData = tCache( tType );

    // Use precomputed data
    for ( uint k = 0; k < tData->number_of_integration_points(); ++k )
    {
        const Matrix< real > & tN  = tData->N( k );
        const Matrix< real > & tdN = tData->dNdXi( k );
        // ... assembly ...
    }
}

// Cleanup after assembly
for ( auto & tEntry : tCache )
{
    delete tEntry.second ;
}
tCache.clear() ;
```

**Performance gain:** the shape functions are evaluated once per element type instead of once per element. No benchmark of the speed-up is kept in the repository.

---

### Advanced: Multi-Order Caching

```cpp
// For problems requiring multiple integration orders (e.g., mass + stiffness)
struct IntegrationKey {
    ElementType type;
    uint order;
    IntegrationScheme scheme;

    bool operator<(const IntegrationKey& other) const {
        if (type != other.type) return type < other.type;
        if (order != other.order) return order < other.order;
        return scheme < other.scheme;
    }
};

std::map<IntegrationKey, IntegrationData*> cache;

IntegrationKey key{ElementType::TET10, 4, IntegrationScheme::GAUSS};
if (!cache.count(key)) {
    cache[key] = new IntegrationData(ElementType::TET10);
    cache[key]->populate(4, IntegrationScheme::GAUSS);
}

IntegrationData* data = cache[key];
```

---

### Block-Level Caching Pattern (BELFEM Best Practice)

```cpp
// Mesh blocks are typically homogeneous (all elements same type)
for (Block* block : mesh->blocks()) {
    ElementType blockType = block->elements()(0)->type();

    // Cache once per block
    IntegrationData data(blockType);
    data.populate(4);

    for (Element* e : block->elements()) {
        // All elements in block use same IntegrationData
        for (uint k = 0; k < data.number_of_integration_points(); ++k) {
            // ... assembly ...
        }
    }
}
```

**When blocks are mixed:** Fall back to element-type caching (Map-based approach).

---

## Element Support Matrix

### Lagrange Elements

| Geometry | Order 1 | Order 2 | Order 3 | Order 4 |
|----------|---------|---------|---------|---------|
| **LINE** | LINE2 | LINE3 | LINE4 | LINE5 |
| **TRI** | TRI3 | TRI6 | TRI10 | TRI15 |
| **QUAD** | QUAD4 | QUAD8, QUAD9 | QUAD16 | — |
| **TET** | TET4 | TET10 | TET20 | TET35 |
| **PENTA** | PENTA6 | PENTA15, PENTA18 | — | — |
| **PYRA** | PYRA5 | PYRA13, PYRA14 | — | — |
| **HEX** | HEX8 | HEX20, HEX27 | HEX64 | — |

**Note:** QUAD8, HEX20 are serendipity elements (no interior nodes).

**Thin-shell types:** `create_lagrange_function()` maps the thin-shell element types onto the volume Lagrange function of the same node layout: `QUAD4TS` → QUAD4, `QUAD9TS` → QUAD9, `PENTA6TS` → PENTA6, `PENTA18TS` → PENTA18, `HEX8TS` and `HEX8TB` → HEX8.

---

### Hermite Elements (C¹ Continuity)

| Element | Bases | Derivatives at Nodes | Use Case |
|---------|-------|----------------------|----------|
| **LINE2** | 4 | Value + slope | Euler-Bernoulli beams (cubic) |
| **QUAD4** | 16 | Value + ∂/∂x + ∂/∂y + ∂²/∂x∂y | Kirchhoff plates (bicubic) |

**Hermite basis interpretation:**
- Each node contributes multiple DOFs (value + derivatives)
- Ensures C¹ continuity across elements (continuous value and slope)

---

### Bernstein Elements

| Element | Order | Bases | Note |
|---------|-------|-------|------|
| LINE2 | 1 | 2 | Returns LAGRANGE (identical) |
| LINE3 | 2 | 3 | True Bernstein basis |
| TRI3 | 1 | 3 | Returns LAGRANGE (identical) |
| TRI6 | 2 | 6 | True Bernstein basis |

**Mathematical note:** Bernstein polynomials = Lagrange polynomials for linear elements. Factory returns Lagrange for efficiency.

---

### Bubble Functions (Stabilization)

Bubble functions are facet-attached enrichment functions (one per edge in 2D, one per face in 3D) used for stabilization in mixed formulations (e.g., inf-sup stability); each vanishes on every facet except its own.

| Element | Bubble Type | Facet functions | Bases per function | Use Case |
|---------|-------------|-----------------|--------------------|----------|
| TRI3 | Edge | 3 (one per edge) | 1 | Stabilize TRI3-based mixed methods |
| TRI6 | Edge | 3 | 2 | Stabilize TRI6-based mixed methods |
| TET4 | Face | 4 (one per face) | 1 | Stabilize TET4-based mixed methods |
| TET10 | Face | 4 | 3 | Stabilize TET10-based mixed methods |

`create_bubble_function( aElementType, aFacet )` returns the function of one facet; its `number_of_bases()` is the "Bases per function" column (`cl_IF_InterpolationFunctionFactory.cpp`).

**Bubble basis properties:**
- Zero on every facet except the one it belongs to
- Non-zero in element interior
- Can be statically condensed (eliminated before global assembly)

---

### Second Derivative Support

Every element the factory can create implements `d2NdXi2` (LINE, TRI, QUAD, TET, PENTA, PYRA, HEX, plus the Hermite, Bernstein and bubble specializations). The `BELFEM_ERROR` in `InterpolationFunctionTemplate` is reached only for a template combination without a specialization, which the factory never instantiates.

---

## Integration Point System

### Volume Integration

```cpp
void initialize_integration_points(
    const ElementType & aElementType,
    Vector<real> & aWeights,
    Matrix<real> & aPoints,
    const uint aIntegrationOrder = 0,
    const IntegrationScheme aIntegrationScheme = IntegrationScheme::GAUSS
);
```

**Output:**
- `aWeights`: (npts,) integration weights
- `aPoints`: parametric coordinates of the points, one column per point. TRI and TET rules carry `ndim + 1` rows: the last row is the remaining barycentric coordinate, 1 − ξ − η (− ζ). Shape functions read only the first `ndim` rows

**Example:**
```cpp
Vector<real> weights;
Matrix<real> points;

initialize_integration_points(ElementType::TET4, weights, points, 2, IntegrationScheme::GAUSS);

// For TET4, order 2:
// npts = 4 (standard 4-point Gauss for tetrahedron)
// points is (4 × 4): rows xi, eta, zeta and 1 - xi - eta - zeta
// weights is (4,)
```

---

### Facet Integration

```cpp
void initialize_integration_points_on_facet(
    const ElementType aElementType,
    const uint aFacetIndex,
    Vector<real> & aWeights,
    Matrix<real> & aPoints,
    const uint aOrder = 0,
    const IntegrationScheme aScheme = IntegrationScheme::GAUSS
);
```

**Parameters:**
- `aFacetIndex`: Local facet index (0-based)
  - TET4: 0-3 (triangular faces)
  - HEX8: 0-5 (quad faces)

**Output:**
- `aPoints`: Parametric coordinates **in the parent element's coordinate system**
- Points lie on the specified facet

**Example:**
```cpp
Vector<real> weights;
Matrix<real> points;

// Face 2 of TET4 (a triangle)
initialize_integration_points_on_facet(ElementType::TET4, 2, weights, points, 3);

// points(0:2, :) are (ξ, η, ζ) coords in TET4 parameter space
// but lie on face 2 (one coordinate fixed)
```

---

### Facet Integration with Orientation (Advanced)

For slave facets in contact/DG:

```cpp
// Internal function (called by IntegrationData::populate_for_slave)
facetintpoints::intpoints_tet(
    const uint aSlaveIndex,
    const uint aOrientation,
    Vector<real> & aWeights,
    Matrix<real> & aPoints,
    const uint aOrder,
    const IntegrationScheme aScheme
);
```

**Orientation:** 0-based index, `Facet::orientation_on_slave() - 1` (see Example 5).

**Supported geometries in `IntegrationData::populate_for_slave()`:**
- TET (via `intpoints_tet`)
- PENTA (via `intpoints_penta`)
- HEX (via `intpoints_hex`)
- QUAD4TS, QUAD9TS (via `intpoints_quad`, with orientation)
- TRI and other QUAD (via `populate_for_slave_tri`, orientation ignored)

PYRA has an oriented slave rule (`facetintpoints::intpoints_pyra`), but `populate_for_slave()` has no PYRA case and aborts with `BELFEM_ERROR`.

---

## Derivative Matrix Conventions

### First Derivatives (dNdXi)

**Dimension:** (ndim × nbases)

**Interpretation:**
```cpp
Matrix<real> dNdXi;  // (ndim × nbases)
shape->dNdXi(xi, dNdXi);

// Row 0: ∂N/∂ξ  for all bases
// Row 1: ∂N/∂η  for all bases (2D/3D)
// Row 2: ∂N/∂ζ  for all bases (3D)
```

**Example (TRI6):**
```cpp
// dNdXi is (2 × 6)
real dN1_dxi = dNdXi(0, 0);  // ∂N₁/∂ξ
real dN1_deta = dNdXi(1, 0); // ∂N₁/∂η
real dN6_dxi = dNdXi(0, 5);  // ∂N₆/∂ξ
real dN6_deta = dNdXi(1, 5); // ∂N₆/∂η
```

---

### Second Derivatives (d2NdXi2)

**Storage format depends on dimensionality:**

#### 1D Elements (LINE)

**Dimension:** (1 × nbases)

```cpp
Matrix<real> d2NdXi2;  // (1 × nbases)
shape->d2NdXi2(xi, d2NdXi2);

// Row 0: ∂²N/∂ξ²
```

---

#### 2D Elements (TRI, QUAD)

**Dimension:** (3 × nbases)

**Row layout:**
```cpp
// d2NdXi2 is (3 × nbases)
// Row 0: ∂²N/∂ξ²
// Row 1: ∂²N/∂η²
// Row 2: ∂²N/(∂ξ∂η)  (mixed derivative)
```

**Example (TRI6):**
```cpp
Matrix<real> d2N;  // (3 × 6)
shape->d2NdXi2(xi, d2N);

real d2N1_dxi2 = d2N(0, 0);      // ∂²N₁/∂ξ²
real d2N1_deta2 = d2N(1, 0);     // ∂²N₁/∂η²
real d2N1_dxideta = d2N(2, 0);   // ∂²N₁/(∂ξ∂η)
```

---

#### 3D Elements (TET, PENTA, PYRA, HEX)

**Dimension:** (6 × nbases)

**Row layout:**
```cpp
// d2NdXi2 is (6 × nbases)
// Row 0: ∂²N/∂ξ²
// Row 1: ∂²N/∂η²
// Row 2: ∂²N/∂ζ²
// Row 3: ∂²N/(∂η∂ζ)
// Row 4: ∂²N/(∂ξ∂ζ)
// Row 5: ∂²N/(∂ξ∂η)
```

**Example (TET10):**
```cpp
Matrix<real> d2N;  // (6 × 10)
shape->d2NdXi2(xi, d2N);

real d2N1_dxi2 = d2N(0, 0);       // ∂²N₁/∂ξ²
real d2N1_deta2 = d2N(1, 0);      // ∂²N₁/∂η²
real d2N1_dzeta2 = d2N(2, 0);     // ∂²N₁/∂ζ²
real d2N1_detadzeta = d2N(3, 0);  // ∂²N₁/(∂η∂ζ)
real d2N1_dxidzeta = d2N(4, 0);   // ∂²N₁/(∂ξ∂ζ)
real d2N1_dxideta = d2N(5, 0);    // ∂²N₁/(∂ξ∂η)
```

---

### Why This Layout?

**Voigt notation compatibility:** The row ordering matches stress/strain tensor storage in solid mechanics:

```cpp
// 2D strain tensor (symmetric 2×2)
ε = [ε_xx  ε_xy]  →  Voigt: [ε_xx, ε_yy, γ_xy]
    [ε_xy  ε_yy]

// d2NdXi2 rows:  [∂²/∂ξ², ∂²/∂η², ∂²/(∂ξ∂η)]  (same order)
```

**Use in plate/shell elements:** Second derivatives directly form curvature matrices.

---

## Performance Considerations

### Direct Evaluation vs. IntegrationData

Creating and evaluating a shape function for every element repeats work that depends only on the element type. With a cached `IntegrationData`, the assembly loop reads stored matrices instead. The repository contains no benchmark comparing the two approaches.

**Recommendation:** use a cached `IntegrationData` in assembly loops.

---

### Memory Footprint

**IntegrationData storage** (TET10, order 4), payload only:
```cpp
// 11 integration points, 10 bases, 3 dimensions; populate() fills every table
mWeights:   11 × sizeof(real)                = 88 bytes
mPoints:    4 × 11 × sizeof(real)            = 352 bytes   // TET rules carry 4 rows
mN:         11 × (1 × 10) × sizeof(real)     = 880 bytes
mPhi:       11 × 10 × sizeof(real)           = 880 bytes
mNvector:   11 × (3 × 30) × sizeof(real)     = 7,920 bytes
mdNdXi:     11 × (3 × 10) × sizeof(real)     = 2,640 bytes
md2NdXi2:   11 × (6 × 10) × sizeof(real)     = 5,280 bytes

Total: ~18 KB per element type
```

The size grows with the number of points and bases, but a cache holds one entry per element type and order, so it stays small compared to mesh data.

---

### When NOT to Precompute

**Adaptive quadrature:** If integration points change per element (e.g., adaptive p-refinement), precomputation loses value.

**Solution:** Use `InterpolationFunction` directly for dynamic integration schemes.

---

## Thread Safety and MPI

### Thread Safety

BELFEM parallelizes with MPI rather than threads, and its internals are deliberately not thread-safe (see `doc/coding_philosophy.md`). This module adds no locks. The classes manage their own state as follows:

**InterpolationFunction:**
- `N`, `dNdXi`, `d2NdXi2`, `param_coords` are `const` and write only into the caller's output matrix

**IntegrationData:**
- `populate*()` reallocates and overwrites every table
- The accessors are `const` and return references to stored tables

**EdgeFunction:**
- `precompute()`, `link()`, `E()`, `C()` and `G()` all write member state (see the Lifetime & Mutation Contract)

**Factories:**
- Both factories have no data members

**Containers:** `Map::operator[]` is non-const and inserts a missing key; use `key_exists()` and `operator()` for lookups. `Cell` is indexed with `operator()`; it has no `operator[]`.

Code that calls into this module from threads must provide its own protection, for example with `#pragma omp critical`.

---

### MPI Compatibility

**All classes are MPI-compatible by design:**
- No global state (each rank creates independent instances)
- No communication required for shape functions/integration
- Works seamlessly with distributed meshes

**Edge elements with distributed meshes:**
```cpp
// Each rank has local mesh elements ( aGroup: fem::Group * )
EdgeFunction* edge = edgeFactory.create_edge_function(aGroup->element_type());
edge->precompute(gaussPoints);
edge->link(localElement);  // Links to rank-local fem::Element
// Edge orientations are part of local mesh data (no communication)
```

---

## Common Patterns

### Pattern 1: Single Element Type Assembly

```cpp
// aElements: Cell< mesh::Element * > &, all of type TET10
IntegrationData tData( ElementType::TET10 );
tData.populate( 4 );

for ( mesh::Element * tElement : aElements )
{
    BELFEM_ASSERT( tElement->type() == ElementType::TET10, "Unexpected element type" );

    for ( uint k = 0; k < tData.number_of_integration_points(); ++k )
    {
        const Matrix< real > & tN  = tData.N( k );
        const Matrix< real > & tdN = tData.dNdXi( k );
        const real tW = tData.weights()( k );

        // Compute Jacobian J from element node coords and dN
        // Assembly: Ke += B^T * D * B * w * detJ
    }
}
```

`BELFEM_ASSERT` is compiled out of release builds, so in a release build an element of another type is not detected here. Use `BELFEM_ERROR` where the check must also run in release builds.

---

### Pattern 2: Multi-Type Mesh with Dispatch

```cpp
// aElements: Cell< mesh::Element * > &
Map< ElementType, IntegrationData * > tCache ;

for ( mesh::Element * tElement : aElements )
{
    const ElementType tType = tElement->type() ;

    if ( ! tCache.key_exists( tType ) )
    {
        IntegrationData * tNew = new IntegrationData( tType );
        tNew->populate( 4 );
        tCache[ tType ] = tNew ;
    }

    const IntegrationData * tData = tCache( tType );

    // ... assembly using tData ...
}

// Cleanup
for ( auto & tEntry : tCache )
{
    delete tEntry.second ;
}
```

---

### Pattern 3: Boundary Condition Application

One `IntegrationData` per local facet index, built before the facet loop. `fem::SideSet` builds its master tables the same way (`fem/kernel/cl_FEM_SideSet.cpp`).

```cpp
// aMesh: Mesh *, aNeumannID: id of the Neumann sideset, master elements of type TET4
const ElementType tMasterType = ElementType::TET4 ;
const uint tNumFacets = mesh::number_of_facets( tMasterType );

Cell< IntegrationData * > tMasterData( tNumFacets, nullptr );
for ( uint f = 0; f < tNumFacets; ++f )
{
    tMasterData( f ) = new IntegrationData( tMasterType );
    tMasterData( f )->populate_for_master( f, 3 );
}

for ( mesh::SideSet * tSideSet : aMesh->sidesets() )
{
    if ( tSideSet->id() != aNeumannID ) continue ;

    for ( mesh::Facet * tFacet : tSideSet->facets() )
    {
        BELFEM_ASSERT( tFacet->master()->type() == tMasterType,
                       "unexpected master element type" );

        const IntegrationData * tData = tMasterData( tFacet->index_on_master() );

        for ( uint k = 0; k < tData->number_of_integration_points(); ++k )
        {
            const Matrix< real > & tN = tData->N( k );
            const real tW = tData->weights()( k );

            // Apply traction: F += N^T * traction * w * dS
        }
    }
}

for ( IntegrationData * tData : tMasterData )
{
    delete tData ;
}
```

---

### Pattern 4: Edge Element Assembly for Electromagnetics

```cpp
// aGroup: fem::Group * whose elements are TET4
IntegrationData tGauss( ElementType::TET4 );
tGauss.populate( 4 );

EdgeFunctionFactory tFactory ;
EdgeFunction * tEdge = tFactory.create_edge_function( ElementType::TET4 );
tEdge->precompute( tGauss.points() );      // once per point set

for ( fem::Element * tElement : aGroup->elements() )
{
    tEdge->link( tElement );               // once per element

    for ( uint k = 0; k < tGauss.number_of_integration_points(); ++k )
    {
        const Matrix< real > & tC = tEdge->C( k );  // Curl
        const real tDetJ = tEdge->det_J() ;
        const real tW = tGauss.weights()( k );

        // Curl-curl term: Ke += C^T * μ^-1 * C * w * detJ
    }
}

delete tEdge ;
```

---

## Development Notes

### Adding New Element Types

**Steps:**

1. **Create shape function class:**
   - Lagrange: Add to `lagrange/` subdirectory
   - Hermite: Add to `hermite/` subdirectory
   - Bernstein: Add to `bernstein/` subdirectory

2. **Specialize template or derive from base:**
   ```cpp
   template <>
   class InterpolationFunctionTemplate<
       GeometryType::NEW_GEO,
       InterpolationType::LAGRANGE,
       3,  // Dimensions
       20  // Number of bases
   > : public InterpolationFunction {
       // Implement N, dNdXi, d2NdXi2, param_coords
   };
   ```

3. **Add to factory:**
   ```cpp
   // In InterpolationFunctionFactory::create_lagrange_function
   case ElementType::NEW_ELEMENT:
       return new InterpolationFunctionTemplate<...>();
   ```

4. **Add integration points:**
   - Update `fn_IF_initialize_integration_points.cpp`
   - Add Gauss quadrature rules for new geometry
   - Update auto-order selection logic

5. **Test:**
   - Partition of unity: ∑ N_i(ξ) = 1
   - Patch test: Constant strain produces constant stress
   - Derivative verification (numerical vs. analytical)

6. **Document:**
   - Update element support matrix in README.md
   - Add to this guide's element list
   - Note any limitations (e.g., d2NdXi2 support)

---

### Adding Nédélec Elements

**Steps:**

1. **Create `EF_<NAME>` class** in `nedelec/` subdirectory:
   ```cpp
   class EF_NEW : public EdgeFunction {
       // mE, mC, mGrad, mJ, mDetJ, mNumDofs are inherited; size mE, mC, mGrad
       // in the constructor: E (ndim × ndofs), C (3 × ndofs in 3D, 1 × ndofs in 2D),
       // G (ndim² × ndofs)
       // ... element-specific Jacobian data ...

       void link(Element* aElement) override;              // fem::Element
       void precompute(const Matrix<real>& aXi) override;
       const Matrix<real>& E(const uint aIndex) override;
       const Matrix<real>& C(const uint aIndex) override;
       const Matrix<real>& G(const uint aIndex) override;
   };
   ```

2. **Implement curl operator:**
   - Reference Nédélec (1980) or Monk (2003) for basis definitions
   - Test against analytical solutions (e.g., constant field)

3. **Add to EdgeFunctionFactory:**
   ```cpp
   case ElementType::NEW_ELEMENT:
       return new EF_NEW();
   ```

4. **Update `num_nedelec_dofs()`:**
   ```cpp
   case ElementType::NEW_ELEMENT:
       return <number_of_edges>;  // Or 2 * num_edges for 2nd order
   ```

5. **Document:**
   - Add to Nédélec element table
   - Specify DOF count and order

---

### Known Limitations

1. **`EdgeFunctionFactory` has no LINE2 or LINE3 case:**
   - `num_nedelec_dofs()` (`nedelec/fn_num_nedelec_dofs.hpp`) returns 1 and 2 for LINE2 and LINE3
   - `EdgeFunctionFactory::create_edge_function()` aborts for both; `EF_LINE3` exists but is not created by the factory

2. **No quadratic thin-shell edge functions:** `num_nedelec_dofs()` returns 3 for `QUAD9TS` and 16 for `PENTA18TS`, but `EdgeFunctionFactory::create_edge_function()` aborts with `BELFEM_ERROR` for both.

3. **Bernstein linear → Lagrange transparency:**
   - Document explicitly that linear Bernstein returns Lagrange
   - Or create separate Bernstein template specializations reporting correct type

---

### Testing Shape Functions

**Partition of unity test:**
```cpp
InterpolationFunction* shape = factory.create_lagrange_function(ElementType::TRI6);

Vector<real> xi(2);
xi(0) = 0.3; xi(1) = 0.4;  // Arbitrary point

Matrix<real> N;
shape->N(xi, N);

real sum = 0.0;
for (uint i = 0; i < shape->number_of_bases(); ++i) {
    sum += N(0, i);
}

BELFEM_ASSERT(std::abs(sum - 1.0) < 1e-12, "Partition of unity failed");
```

**Derivative verification (finite difference):**
```cpp
real h = 1e-7;
Vector<real> xi(2), xi_plus(2);
xi(0) = 0.3; xi(1) = 0.4;
xi_plus = xi; xi_plus(0) += h;

Matrix<real> N, N_plus, dN_analytical;
shape->N(xi, N);
shape->N(xi_plus, N_plus);
shape->dNdXi(xi, dN_analytical);

for (uint i = 0; i < shape->number_of_bases(); ++i) {
    real dN_numerical = (N_plus(0,i) - N(0,i)) / h;
    real error = std::abs(dN_numerical - dN_analytical(0,i));
    BELFEM_ASSERT(error < 1e-5, "Derivative verification failed");
}
```

---

## See Also

### Internal Documentation

- **FEM Kernel**: `../../kernel/doc/` - Uses IntegrationData for assembly
- **Mesh Module**: `../../../mesh/doc/` - Element definitions, facet orientations
- **Homology Module**: `../../../homology/doc/` - Cohomology cuts for Nédélec elements
- **Linear Algebra**: `../../../linalg/doc/` - Matrix/Vector for shape storage

### Literature

**Shape Functions:**
- **Zienkiewicz & Taylor**, "The Finite Element Method" Vol. 1, Ch. 6-7
- **Hughes**, "The Finite Element Method", Ch. 3 (Isoparametric elements)
- **Bathe**, "Finite Element Procedures", Ch. 5

**Nédélec Elements:**
- **Nédélec (1980)**, "Mixed finite elements in ℝ³", Numer. Math. 35, 315-341
- **Monk (2003)**, "Finite Element Methods for Maxwell's Equations", Oxford
- **BELFEM papers**: `literature/papers/fem/messe2023.txt` (h-φ formulation), `literature/papers/fem/arsenault2023.txt` (magnetodynamic coupling)

**Integration:**
- **Stroud (1971)**, "Approximate Calculation of Multiple Integrals"
- **Dunavant (1985)**, "High degree efficient symmetrical Gaussian quadrature rules for the triangle", Int. J. Numer. Methods Eng. 21, 1129-1148

### Project References

- **Project README**: `../../../../README.md`
- **Claude Instructions**: `../../../../CLAUDE.md`
- **Coding Philosophy**: `../../../../doc/coding_philosophy.md`
- **Documentation Guidelines**: `../../../../doc/documentation_guidelines.md`

---

**Contributors:** Based on code by Christian Messe and Gregory Giard

**Prepared by the BELFEM development team.**
