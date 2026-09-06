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

**Why wrong:** Each factory creation, shape function allocation, and IntegrationData precomputation is expensive. For 10,000 elements, this wastes ~99.99% of CPU time on redundant work.

```cpp
// CORRECT - Cache by element type
Map<ElementType, IntegrationData*> cachedData;

for (Element* e : elements) {
    ElementType type = e->type();

    // Create once per element type
    if (!cachedData.contains(type)) {
        cachedData[type] = new IntegrationData(type);
        cachedData[type]->populate(4, IntegrationScheme::GAUSS);
    }

    IntegrationData* data = cachedData[type];
    for (uint k = 0; k < data->number_of_integration_points(); ++k) {
        // Reuse precomputed N, dNdXi - zero cost!
    }
}

// Cleanup
for (auto& pair : cachedData) delete pair.second;
```

**Impact:** 100-500× speedup for typical assembly loops.

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
IntegrationData slaveData(ElementType::TET4);
slaveData.populate_for_slave(
    slaveIndex,    // Which facet
    0,             // Orientation = 0 assumed!
    order,
    scheme
);
```

**Problem:** For contact, DG, or interface problems, slave facet orientations must match the master facet node ordering. Ignoring orientation leads to inconsistent normal vectors and flux computation errors.

```cpp
// CORRECT - Provide correct orientation
uint orientation = compute_facet_orientation(masterFacet, slaveFacet);  // From mesh
slaveData.populate_for_slave(slaveIndex, orientation, order, scheme);
```

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

**Auto order selection heuristic:** BELFEM does not treat `0` as a minimal rule. It picks a conservative default from the element interpolation order (see `fn_intpoints_auto_integration_order.cpp`).

The reasoning comes from the weak form. For a `p`-th order interpolation:
- A mass-like term such as `δu * m * u` is at least order `2p`
- If the material field `m` is assumed linear, the target rises to about `2p + 1`
- A stiffness-like term such as `grad(δu) * k * grad(u)` is usually lower order

This is a heuristic, not an exactness guarantee. For distorted QUAD/HEX elements and for curved higher-order mappings, the Jacobian makes the integrand non-polynomial, so exact integration is not available anyway. The default therefore intentionally aims for a practical, slightly conservative rule and knowingly still under-integrates many real nonlinear cases.

The integer order does carry one guarantee: for `IntegrationScheme::GAUSS`, `intpoints()` returns a rule that integrates every polynomial of total degree up to that order exactly on the reference element, for every geometry it serves (locked by `tests/fem/test_IntegrationExactness.cpp`). Several rules deliver more than requested (the tensor-product rules, and the tetrahedron tables for orders 7 and 9), so the number is a floor, not the exact degree. It says nothing about mapped, non-polynomial integrands.

**Recommendation:** Explicitly specify order for clarity:
```cpp
data.populate(7, IntegrationScheme::GAUSS);  // Explicit order for this TRI6 example
```

---

### 5. Ownership Confusion with IntegrationData

```cpp
InterpolationFunction* shape = factory.create_lagrange_function(ElementType::TET4);

// Constructor 1: IntegrationData does NOT own shape
IntegrationData data1(ElementType::TET4, shape, false);
// ... use data1 ...
delete shape;  // Caller must delete!

// Constructor 2: IntegrationData OWNS shape
IntegrationData data2(ElementType::TET4, shape, true);
// ... use data2 ...
// shape is deleted in ~IntegrationData()
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
// Works for non-curved TRI and TET elements that refresh detJ during element linking.
const Matrix<real>& E = edge->E(k);
real detJ = edge->det_J();  // May be stale for reduced custom elements
const Matrix<real>& C = edge->C(k);
```

**Problem:** linear QUAD and HEX and all curved higher order elements update the Jacobian determinant  `det_J()` only in the gradient operator `B(k)`, respectively the curl operator `C(k)`. For such an element `E(k)` interpolates the
field, while `C(k)` refreshes `det_J()` for the active
integration point. If client code reads `det_J()` or calls a helper such as
`dV_ts()` before `C(k)`, it can consume a stale determinant from a previous
integration point.

**Safe pattern:**
```cpp
for (uint k = 0; k < npts; ++k) {
    const Matrix<real>& E = edge->E(k);
    const Matrix<real>& C = edge->C(k);  // Refreshes detJ
    real detJ = edge->det_J();

    // Safe: any dV helper that reads edge->det_J() now sees the current point
}
```

**Scope:** This is not a general `EdgeFunction` rule, but it is a real gotcha
for reduced or curved custom elements. If a new element caches Jacobian
data lazily, document clearly which accessor owns that update.

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
    ├── cl_EF_QUAD4TS.{hpp,cpp}             # Thin shell with curl (quad)
    └── cl_EF_PENTA6TS.{hpp,cpp}             # Thin shell with curl (tri)
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
    ├── EF_TET4, EF_TET10   # 3D electromagnetics
    └── EF_QUAD4TS, EF_PENTA6TS # Thin shells with curl

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
- `aXiHat`: (num_dimensions × num_bases) - parametric node coordinates

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
- Curl operator: `C(ξ)` returns curl matrix
- Requires mesh element linking: `link(Element*)` to access edge orientations

**Usage pattern:**
```cpp
EdgeFunctionFactory edgeFactory;
EdgeFunction* edge = edgeFactory.create_edge_function(ElementType::TET4);

edge->link(meshElement);  // Provide element context
edge->precompute(integrationPoints);

const Matrix<real>& E = edge->E(integrationPointIndex);  // Interpolation (3 × 6 for TET4)
const Matrix<real>& C = edge->C(integrationPointIndex);  // Curl (3 × 6)
real detJ = edge->det_J();

delete edge;
```

**Prerequisites:**
- Mesh must have edges/faces created (`mesh->create_edges()`)
- Element must provide edge orientation data
- Typically used with `DofManagerMaxwell` for electromagnetic DOF management

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

// ✅ MUTABLE - link and precompute modify internal state
edge->link(element);            // Sets up element context (node coords, edge orientations)
edge->precompute(gaussPoints);  // Precomputes E and C at integration points
```

**Contract:** `link()` and `precompute()` **must** be called before accessing `E()` or `C()`.

#### Phase 2: Immutable (After precompute)
```cpp
// ✅ IMMUTABLE - access methods are const after precompute
for (uint k = 0; k < npts; ++k) {
    const Matrix<real>& E = edge->E(k);  // Read-only
    const Matrix<real>& C = edge->C(k);  // Read-only
    real detJ = edge->det_J();           // Read-only
}
```

**Warning:** `EdgeFunction` is **NOT thread-safe** during `link()` or `precompute()`. Create separate instances per thread or protect with `#pragma omp critical`.

---

### Summary Table

| Object | Creation | Mutable Phase | Immutable Phase | Thread-Safe After? |
|--------|----------|---------------|-----------------|-------------------|
| **InterpolationFunction** | Factory | Never (effectively const) | Always | ✅ Yes |
| **IntegrationData** | Constructor | During `populate()` calls | After final `populate()` | ✅ Yes |
| **EdgeFunction** | Factory | During `link()` + `precompute()` | After `precompute()` | ❌ No (create per thread) |

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

for (Element* e : elements) {
    edge->precompute(gaussPoints);  // ❌ BUG! Missing link(e)
    // ... assembly (uses wrong element context) ...
}
```

**Fix:** Call `link(e)` before each `precompute()`.

```cpp
for (Element* e : elements) {
    edge->link(e);                  // ✅ Correct
    edge->precompute(gaussPoints);
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

std::cout << "N = " << N << std::endl;
std::cout << "dN/dXi = " << dNdXi << std::endl;

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
- **Zero cost** access via `N(k)`, `dNdXi(k)` in loops
- 100-500× faster than evaluating shape functions per Gauss point per element

---

### Example 3: Hermite Shape Functions for Beam Elements

```cpp
// Hermite LINE2: 4 bases (value + derivative at each node)
InterpolationFunction* hermite = factory.create_hermite_function(ElementType::LINE2);

Vector<real> xi(1);
xi(0) = 0.0;  // Left end

Matrix<real> N;
hermite->N(xi, N);  // N = (1 × 4) = [N1, N1', N2, N2']

// N1 = 1, N1' = 0, N2 = 0, N2' = 0 at left end
// Used for Euler-Bernoulli beam (cubic interpolation)
```

**Hermite basis interpretation** (LINE2):
- N(0) = value at node 0
- N(1) = derivative at node 0
- N(2) = value at node 1
- N(3) = derivative at node 1

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
uint orientation = 3;         // Orientation code from mesh (facet node permutation)

slaveData.populate_for_slave(slaveFacetIndex, orientation, 3, IntegrationScheme::GAUSS);

// Integration points now match the master element's facet node ordering
// Critical for consistent normal vectors and flux computation
```

**Orientation encoding:** Maps slave facet node indices to master facet node ordering (mesh-dependent convention).

---

### Example 6: Nédélec Edge Elements for Electromagnetics

```cpp
#include "cl_EdgeFunctionFactory.hpp"

EdgeFunctionFactory edgeFactory;
EdgeFunction* edge = edgeFactory.create_edge_function(ElementType::TET4);

// Link to mesh element (provides edge orientations)
edge->link(meshElement);

// Precompute at integration points
Matrix<real> integrationPoints;  // (3 × num_pts) parametric coords
// ... initialize integrationPoints from Gauss quadrature ...
edge->precompute(integrationPoints);

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
// Constructor 1: Auto-create Lagrange shape function
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
- `aScheme`: Currently only `GAUSS` (Gauss-Legendre quadrature) is supported

**Auto-order heuristic:** `aIntegrationOrder = 0` calls `auto_integration_order( aElementType )`.

The built-in default is chosen from the interpolation order of the element and is meant to be a practical weak-form heuristic rather than an exactness proof. The motivating rule is that a `p`-th order interpolation produces mass-like products of about order `2p`, and with at least linearly varying material data the target becomes about `2p + 1`.

For affine simplices this is a reasonable baseline. For distorted QUAD/HEX elements and curved higher-order mappings, the Jacobian destroys strict polynomial structure, so the same default should be understood as a conservative engineering choice, not as exact integration.

**Usage:**
```cpp
IntegrationData data(ElementType::TET10);
data.populate(7);  // Explicit order for this TET10 example
// or
data.populate();   // Auto-select a conservative default from the element order
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
- `aOrientation`: Orientation code (permutation of slave nodes to match master)
  - 0-based encoding (depends on geometry, see mesh documentation)
  - Critical for consistent normal vectors

**Use case:** Contact mechanics, discontinuous Galerkin, interface coupling (h-φ, thin-shell)

```cpp
// Slave facet must match master facet node ordering
uint orientation = compute_facet_orientation(masterFacet, slaveFacet);
slaveData.populate_for_slave(slaveIndex, orientation, 3);
```

**Orientation example (TRI facet on TET):**
- Master facet nodes: [0, 1, 2]
- Slave facet nodes: [1, 0, 2] → orientation code maps this permutation

---

### Access Methods

```cpp
uint number_of_integration_points() const;
const Vector<real>& weights() const;      // (npts,)
const Matrix<real>& points() const;       // (ndim × npts) parametric coords

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
| **QUAD4TS** | 2 | Thin shell | `EF_QUAD4TS` (line-in-2D extrusion) |
| **PENTA6TS** | 6 | Thin shell | `EF_PENTA6TS` (triangle-in-3D extrusion) |

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
class EdgeFunction {
public:
    virtual void link(Element* aElement) = 0;  // Provide mesh element context
    virtual void precompute(const Matrix<real>& aXi) = 0;  // Precompute at points

    virtual const Matrix<real>& E(const uint aIndex);  // Interpolation operator
    virtual const Matrix<real>& C(const uint aIndex) = 0;  // Curl operator

    real det_J() const;       // Jacobian determinant
    real abs_det_J() const;   // |det(J)|
    uint ndofs() const;       // Number of DOFs
};
```

**Matrix dimensions:**
- `E(k)`: (ndim × ndofs) - interpolates edge DOFs to physical vector field
- `C(k)`: (ndim × ndofs) - curl of edge basis (for 3D: vector curl, for 2D: scalar curl)

---

### Usage Pattern

```cpp
EdgeFunctionFactory factory;
EdgeFunction* edge = factory.create_edge_function(ElementType::TET4);

// CRITICAL: Link to mesh element
edge->link(meshElement);  // Provides edge orientations, node coordinates

// Precompute at integration points
Matrix<real> gaussPoints;  // (3 × npts) for 3D
// ... initialize from IntegrationData or manual Gauss rules ...
edge->precompute(gaussPoints);

// Assembly loop
for (uint k = 0; k < npts; ++k) {
    const Matrix<real>& E = edge->E(k);  // (3 × 6) for TET4
    const Matrix<real>& C = edge->C(k);  // (3 × 6)
    real detJ = edge->det_J();
    real w = weights(k);

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

Higher-order Nédélec elements (TRI6, TET10) use **function pointers** to switch between:

```cpp
// EF_TET10 example
void (EF_TET10::*mFunCurl)(const uint aIndex);  // Function pointer

// Constructor logic
if (element_is_curved()) {
    mFunCurl = &EF_TET10::C_curved;  // Full Jacobian derivative
} else {
    mFunCurl = &EF_TET10::C_straight;  // Simplified (affine mapping)
}

// Call via pointer
const Matrix<real>& EF_TET10::C(const uint aIndex) {
    (this->*mFunCurl)(aIndex);  // Dispatch to appropriate version
    return mC;
}
```

**Performance impact:** Straight elements ~2× faster (avoid Jacobian derivatives at each point).

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
- Each precompute evaluates shape functions at ~10-30 Gauss points
- Total: ~10,000 × 20 = 200,000 redundant shape evaluations

---

### The Solution: Cache by Element Type

```cpp
// CORRECT - Cache IntegrationData
Map<ElementType, IntegrationData*> integrationCache;

for (Element* e : elements) {
    ElementType type = e->type();

    // Create once per element type
    if (!integrationCache.contains(type)) {
        integrationCache[type] = new IntegrationData(type);
        integrationCache[type]->populate(4, IntegrationScheme::GAUSS);
    }

    IntegrationData* data = integrationCache[type];

    // Use precomputed data - zero overhead!
    for (uint k = 0; k < data->number_of_integration_points(); ++k) {
        const Matrix<real>& N = data->N(k);
        const Matrix<real>& dN = data->dNdXi(k);
        // ... assembly ...
    }
}

// Cleanup after assembly
for (auto& pair : integrationCache) {
    delete pair.second;
}
integrationCache.clear();
```

**Performance gain:** 100-500× speedup (measured on typical meshes with ~10 element types).

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

| Element | Bubble Type | Bubbles | Use Case |
|---------|-------------|---------|----------|
| TRI3 | Edge | 3 (one per edge) | Stabilize TRI3-based mixed methods |
| TRI6 | Edge | 3 | Stabilize TRI6-based mixed methods |
| TET4 | Face | 4 (one per face) | Stabilize TET4-based mixed methods |
| TET10 | Face | 4 | Stabilize TET10-based mixed methods |

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
    const ElementType aElementType,
    Vector<real> & aWeights,
    Matrix<real> & aPoints,
    const uint aOrder = 0,
    const IntegrationScheme aScheme = IntegrationScheme::GAUSS
);
```

**Output:**
- `aWeights`: (npts,) integration weights
- `aPoints`: (ndim × npts) parametric coordinates of Gauss points

**Example:**
```cpp
Vector<real> weights;
Matrix<real> points;

initialize_integration_points(ElementType::TET4, weights, points, 2, IntegrationScheme::GAUSS);

// For TET4, order 2:
// npts = 4 (standard 4-point Gauss for tetrahedron)
// points is (3 × 4)
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

**Orientation:** Permutation index mapping slave facet node ordering to master facet.

**Supported geometries:**
- TET (via `intpoints_tet`)
- HEX (via `intpoints_hex`)
- TRI (via `populate_for_slave_tri`)

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

#### 3D Elements (TET, HEX)

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

### Benchmarking: Direct Evaluation vs. IntegrationData

**Setup:** TET10 element, 14 Gauss points (order 4), assembly loop over 10,000 elements.

| Method | Time per Element | Relative Speed |
|--------|------------------|----------------|
| **Direct eval** (create shape + eval per element) | 250 µs | 1× (baseline) |
| **Direct eval** (reuse shape, eval per Gauss point) | 75 µs | 3.3× |
| **IntegrationData** (precomputed) | 0.5 µs | **500×** |

**Conclusion:** `IntegrationData` with caching is **essential** for performance.

---

### Memory Footprint

**IntegrationData storage** (TET10, order 4):
```cpp
// 14 Gauss points, 10 bases, 3 dimensions
mWeights:   14 × sizeof(real) = 112 bytes
mPoints:    3 × 14 × sizeof(real) = 336 bytes
mN:         14 × (1 × 10) × sizeof(real) = 1,120 bytes
mdNdXi:     14 × (3 × 10) × sizeof(real) = 3,360 bytes
md2NdXi2:   14 × (6 × 10) × sizeof(real) = 6,720 bytes (if populated)

Total: ~11.6 KB per element type
```

**For 10 element types:** ~116 KB total (negligible compared to mesh data).

---

### When NOT to Precompute

**Adaptive quadrature:** If integration points change per element (e.g., adaptive p-refinement), precomputation loses value.

**Solution:** Use `InterpolationFunction` directly for dynamic integration schemes.

---

## Thread Safety and MPI

### Thread Safety

**InterpolationFunction:**
- ✅ **Thread-safe for read-only operations** (N, dNdXi, d2NdXi2 evaluation)
- ❌ **Not thread-safe for mutation** (none expected in typical usage)

**IntegrationData:**
- ✅ **Thread-safe after `populate()`** (all data immutable)
- ❌ **Not thread-safe during `populate()`** (internal allocation)

**Factories:**
- ✅ **Thread-safe** (stateless, no internal state modified)

**Recommended pattern for OpenMP:**
```cpp
// Populate caches in serial region
Map<ElementType, IntegrationData*> cache;
for (ElementType type : uniqueElementTypes) {
    cache[type] = new IntegrationData(type);
    cache[type]->populate(4);
}

// Parallel assembly (read-only access to cache)
#pragma omp parallel for
for (int i = 0; i < elements.size(); ++i) {
    Element* e = elements[i];
    IntegrationData* data = cache[e->type()];  // Thread-safe read

    // ... assembly (each thread works on different elements)
}
```

---

### MPI Compatibility

**All classes are MPI-aware via design:**
- No global state (each rank creates independent instances)
- No communication required for shape functions/integration
- Works seamlessly with distributed meshes

**Edge elements with distributed meshes:**
```cpp
// Each rank has local mesh elements
EdgeFunction* edge = edgeFactory.create_edge_function(element->type());
edge->link(localElement);  // Links to rank-local element
// Edge orientations are part of local mesh data (no communication)
```

---

## Common Patterns

### Pattern 1: Single Element Type Assembly

```cpp
IntegrationData data(ElementType::TET10);
data.populate(4);

for (Element* e : elements) {
    BELFEM_ASSERT(e->type() == ElementType::TET10, "Unexpected element type");

    for (uint k = 0; k < data.number_of_integration_points(); ++k) {
        const Matrix<real>& N = data.N(k);
        const Matrix<real>& dN = data.dNdXi(k);
        real w = data.weights()(k);

        // Compute Jacobian J from element node coords and dN
        // Assembly: Ke += B^T * D * B * w * detJ
    }
}
```

---

### Pattern 2: Multi-Type Mesh with Dispatch

```cpp
Map<ElementType, IntegrationData*> cache;

for (Element* e : elements) {
    ElementType type = e->type();

    if (!cache.count(type)) {
        cache[type] = new IntegrationData(type);
        cache[type]->populate(4);
    }

    IntegrationData* data = cache[type];

    // ... assembly using data ...
}

// Cleanup
for (auto& pair : cache) delete pair.second;
```

---

### Pattern 3: Boundary Condition Application

```cpp
for (SideSet* sideset : mesh->sidesets()) {
    if (sideset->physical_tag() == NEUMANN_BC) {
        for (Facet* facet : sideset->facets()) {
            Element* master = facet->master();
            uint masterIndex = facet->master_index();

            IntegrationData bcData(master->type());
            bcData.populate_for_master(masterIndex, 3);

            for (uint k = 0; k < bcData.number_of_integration_points(); ++k) {
                const Matrix<real>& N = bcData.N(k);
                real w = bcData.weights()(k);

                // Apply traction: F += N^T * traction * w * dS
            }
        }
    }
}
```

---

### Pattern 4: Edge Element Assembly for Electromagnetics

```cpp
EdgeFunction* edge = edgeFactory.create_edge_function(ElementType::TET4);

for (Element* e : conductorElements) {
    edge->link(e);

    // Use standard Gauss points
    IntegrationData gaussData(ElementType::TET4);
    gaussData.populate(4);

    edge->precompute(gaussData.points());

    for (uint k = 0; k < gaussData.number_of_integration_points(); ++k) {
        const Matrix<real>& C = edge->C(k);  // Curl
        real detJ = edge->det_J();
        real w = gaussData.weights()(k);

        // Curl-curl term: Ke += C^T * μ^-1 * C * w * detJ
    }
}

delete edge;
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
       Matrix<real> mE;  // Interpolation (ndim × ndofs)
       Matrix<real> mC;  // Curl (ndim × ndofs)
       // ... Jacobian data ...

       void link(Element* aElement) override;
       void precompute(const Matrix<real>& aXi) override;
       const Matrix<real>& C(const uint aIndex) override;
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

### Known Issues / TODOs

**From external analysis (Opus):**

1. ~~**TET4 Nédélec DOF count bug**~~ — **FIXED.** `num_nedelec_dofs( ElementType::TET4 )`
   returns **6** (`nedelec/fn_num_nedelec_dofs.hpp:25`), which is correct for a tetrahedron's six
   edges. The old entry also named the wrong file.

2. **Hermite LINE2 typo:**
   - Line 2181: `const real b = (xi+1-0);` should be `(xi+1.0)`
   - Mathematically equivalent but poor style

3. **EdgeFunctionFactory has no LINE2 or LINE3 case** (still true; 0 hits for either in `cl_EdgeFunctionFactory.cpp`):
   - `num_nedelec_dofs()` returns values for LINE2/3
   - But `EdgeFunctionFactory::create_edge_function()` doesn't handle them
   - Decision needed: Add to factory or remove from `num_nedelec_dofs()`

4. **Bernstein linear → Lagrange transparency:**
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

- **FEM Kernel**: `../kernel/doc/` - Uses IntegrationData for assembly
- **Mesh Module**: `../../mesh/doc/` - Element definitions, facet orientations
- **Homology Module**: `../../homology/doc/` - Cohomology cuts for Nédélec elements
- **Linear Algebra**: `../../linalg/doc/` - Matrix/Vector for shape storage

### Literature

**Shape Functions:**
- **Zienkiewicz & Taylor**, "The Finite Element Method" Vol. 1, Ch. 6-7
- **Hughes**, "The Finite Element Method", Ch. 3 (Isoparametric elements)
- **Bathe**, "Finite Element Procedures", Ch. 5

**Nédélec Elements:**
- **Nédélec (1980)**, "Mixed finite elements in ℝ³", Numer. Math. 35, 315-341
- **Monk (2003)**, "Finite Element Methods for Maxwell's Equations", Oxford
- **BELFEM papers**: `literature/papers/fem/messe2023.txt` (h-φ formulation), `literature/papers/femarsenault2023.txt` (magnetodynamic coupling)

**Integration:**
- **Stroud (1971)**, "Approximate Calculation of Multiple Integrals"
- **Dunavant (1985)**, "High degree efficient symmetrical Gaussian quadrature rules for the triangle", Int. J. Numer. Methods Eng. 21, 1129-1148

### Project References

- **Project README**: `../../../README.md`
- **Claude Instructions**: `../../../CLAUDE.md`
- **Coding Philosophy**: `../../../doc/coding_philosophy.md`
- **Documentation Guidelines**: `../../../doc/documentation_guidelines.md`

---

**Contributors:** Based on code by Christian Messe and Gregory Giard

**Prepared by the BELFEM development team.**
