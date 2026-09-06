# Interpolation Module Documentation {#fem_interpolation_index}

**Module:** src/fem/interpolation
**Purpose:** Index of documentation for BELFEM's shape function and integration module

---

## Overview

The `interpolation` module provides shape functions, integration point management, and Nédélec (edge) elements for finite element analysis. It supports:

- Shape function evaluation (Lagrange, Hermite, Bernstein, Bubble)
- Integration point generation (Gauss quadrature, volume and facet)
- Nédélec edge elements for electromagnetic applications
- Precomputed integration data for efficient assembly

---

## Documentation Files

### User Guides

- **[interpolation_usage_guide.md](interpolation_usage_guide.md)** - Comprehensive usage guide for the interpolation module
  - Shape function evaluation workflow
  - Integration data precomputation
  - Nédélec edge elements for electromagnetics
  - Factory patterns and caching strategies
  - Element support matrix
  - Common patterns and pitfalls

### Nédélec Edge Elements

- **[nedelec_derivation.md](nedelec_derivation.md)** - The mathematical derivation (extracted 2026-08-14 from Christian's pre-BELFEM theory notes)
  - Interpolation operators N/E/B/C in 3D, 2D and axisymmetric form
  - Barycentric coordinates, Lagrange TRI3/TRI6, geometry Jacobian and the J-transpose pitfall
  - Whitney edge functions for TRI3/TRI6/TET4/TET10, curl operators, face-function redundancy and ownership
  - Second-order circulation convention (**TRI6 only**: its parent-edge functions carry circulation 1/2; TET10's are twice those and carry unit circulation)
  - Edge/face generation concept (sort-unique keying)
  - Flags the known `EF_TET4::E()` edge defect (E and C disagree on one edge)

- **[nedelec.md](nedelec.md)** - General Nédélec framework
  - Role in H-phi formulations and why edge elements are needed
  - The `EdgeFunction` base class, state, and lifecycle (link → precompute → E/C)
  - Factory pattern and DOF counts per element type
  - Edge orientation and the `s_k` sign convention
  - Unit-circulation DOF convention (`∫ E_k · dℓ = δ_{jk}`) and why it matters for shared edges
  - Volume elements (TRI3, TRI6, TET4, TET10, LINE3) and the Whitney 1-form construction
  - Implementation checklist and common pitfalls

- **[nedelec_thinshell.md](nedelec_thinshell.md)** - Thin-shell edge elements
  - Motivation for dimensional reduction in REBCO tape modeling
  - Shared infrastructure: pseudo-inverse Jacobian, normal/binormal, layer thickness
  - `QUAD4TS` (line-in-2D) and `PENTA6TS` (triangle-in-3D) surface thin-shell elements

---

## Quick Reference

### Entry Point Classes

| Class | File | Purpose |
|-------|------|---------|
| **`InterpolationFunction`** | cl_IF_InterpolationFunction.hpp | Abstract base for shape functions |
| **`InterpolationFunctionFactory`** | cl_IF_InterpolationFunctionFactory.hpp | Creates shape functions by element/type |
| **`IntegrationData`** | cl_IF_IntegrationData.hpp | **Main workhorse**: precomputed integration bundle |
| **`EdgeFunction`** | (in nedelec/) | Abstract base for Nédélec elements |
| **`EdgeFunctionFactory`** | cl_EdgeFunctionFactory.hpp | Creates edge functions by element type |

### Interpolation Types

| Type | Use Case | Elements |
|------|----------|----------|
| **LAGRANGE** | Standard FEM | All element geometries, orders 1-4 |
| **HERMITE** | C¹ continuity (beams, plates) | LINE2, QUAD4 |
| **BERNSTEIN** | Bézier-like interpolation | LINE, TRI |
| **Bubble** | Stabilization (mixed formulations) | TRI3/6, TET4/10 |

### Common Operations

```cpp
// Create shape function via factory
InterpolationFunctionFactory tFactory;
InterpolationFunction* tFunc = tFactory.create_lagrange_function(ElementType::TET10);

// Evaluate at a point
Vector<real> xi(3);  // Parametric coordinates
Matrix<real> N;      // Shape values (1 x num_bases)
Matrix<real> dNdXi;  // Derivatives (ndim x num_bases)

tFunc->N(xi, N);
tFunc->dNdXi(xi, dNdXi);

delete tFunc;  // Caller owns

// RECOMMENDED: Use IntegrationData for assembly loops
IntegrationData tData(ElementType::TET10);
tData.populate(4, IntegrationScheme::GAUSS);  // Order 4 Gauss quadrature

for(uint k = 0; k < tData.number_of_integration_points(); ++k) {
    const Matrix<real>& N_k = tData.N(k);          // Precomputed!
    const Matrix<real>& dN_k = tData.dNdXi(k);     // Precomputed!
    real w_k = tData.weights()(k);

    // Element integration: Ke += N^T * D * N * w * detJ
}
```

---

## Shape Function API

### InterpolationFunction Interface

```cpp
// Parametric node coordinates (ndim x nbases)
void param_coords(Matrix<real> & aXiHat) const;

// Shape function values (1 x nbases)
void N(const Vector<real> & aXi, Matrix<real> & aN) const;

// First derivatives (ndim x nbases)
void dNdXi(const Vector<real> & aXi, Matrix<real> & adNdXi) const;

// Second derivatives (see derivative matrix conventions below)
void d2NdXi2(const Vector<real> & aXi, Matrix<real> & ad2NdXi2) const;

// Metadata
uint number_of_bases() const;
uint number_of_dimensions() const;
GeometryType geometry_type() const;
ElementType element_type() const;
InterpolationType interpolation_type() const;
```

### Derivative Matrix Conventions

**Second derivatives** (`d2NdXi2`) storage:

| Dimension | Matrix Size | Row Contents |
|-----------|-------------|--------------|
| **1D** | `(1 x nbases)` | Row 0: ∂²N/∂ξ² |
| **2D** | `(3 x nbases)` | Row 0: ∂²N/∂ξ², Row 1: ∂²N/∂η², Row 2: ∂²N/(∂ξ∂η) |
| **3D** | `(6 x nbases)` | Rows: ∂²N/∂ξ², ∂²N/∂η², ∂²N/∂ζ², ∂²N/(∂η∂ζ), ∂²N/(∂ξ∂ζ), ∂²N/(∂ξ∂η) |

---

## IntegrationData Quick Reference

### Constructors

```cpp
// Constructor 1: Auto-create shape function
IntegrationData(ElementType, InterpolationType = LAGRANGE);

// Constructor 2: Provide shape function (optional ownership)
IntegrationData(ElementType, InterpolationFunction*, bool aClaimOwnership);
```

### Populate Methods

```cpp
// Volume integration
void populate(uint aIntegrationOrder = 0,  // 0 = auto
              IntegrationScheme = GAUSS);

// Master facet (for boundary conditions)
void populate_for_master(uint aMasterIndex, uint aOrder = 0, ...);

// Slave facet (with orientation, for contact/DG)
void populate_for_slave(uint aSlaveIndex, uint aOrientation = 0, uint aOrder = 0, ...);
```

### Access Methods

```cpp
uint number_of_integration_points() const;
const Vector<real> & weights() const;           // Integration weights
const Matrix<real> & points() const;            // Parametric points (ndim x npts)

const Vector<real> & phi(uint k) const;         // Shape as vector (nbases x 1)
const Matrix<real> & N(uint k) const;           // Shape as matrix (1 x nbases)
const Matrix<real> & Nvector(uint k) const;     // For vector fields
const Matrix<real> & dNdXi(uint k) const;       // First derivatives
const Matrix<real> & d2NdXi2(uint k) const;     // Second derivatives
```

---

## Nédélec Edge Elements

For electromagnetic applications (H(curl) formulations). See [nedelec.md](nedelec.md) for the general framework and [nedelec_thinshell.md](nedelec_thinshell.md) for the thin-shell variants.

| Element | DOFs | Class | Use Case |
|---------|------|-------|----------|
| LINE3 | 2 | EF_LINE3 | Higher-order 1D edge |
| TRI3 | 3 | EF_TRI3 | 2D electromagnetics |
| TRI6 | 8 | EF_TRI6 | Higher-order 2D |
| TET4 | 6 | EF_TET4 | 3D electromagnetics |
| TET10 | 20 | EF_TET10 | Higher-order 3D |
| QUAD4TS | 2 | EF_QUAD4TS | 2D thin-shell (LINE2 sideset extrusion) |
| PENTA6TS | 6 | EF_PENTA6TS | 3D thin-shell (TRI3 sideset extrusion) |

```cpp
EdgeFunctionFactory tEdgeFactory;
EdgeFunction* tEdge = tEdgeFactory.create_edge_function(ElementType::TET4);

tEdge->link(aElement);  // Provide mesh element context
tEdge->precompute(integrationPoints);

const Matrix<real>& E = tEdge->E(integrationPointIndex);  // Interpolation (ndim x ndofs)
const Matrix<real>& C = tEdge->C(integrationPointIndex);  // Curl operator (ndim x ndofs)
real detJ = tEdge->det_J();
```

---

## Element Support Matrix

### Lagrange Elements (Full Support)

| Geometry | Elements | Order Range |
|----------|----------|-------------|
| LINE | LINE2, LINE3, LINE4, LINE5 | 1-4 |
| TRI | TRI3, TRI6, TRI10, TRI15 | 1-4 |
| QUAD | QUAD4, QUAD8, QUAD9, QUAD16 | 1-3 |
| TET | TET4, TET10, TET20, TET35 | 1-4 |
| PENTA | PENTA6, PENTA15, PENTA18 | 1-2 |
| PYRA | PYRA5, PYRA13, PYRA14 | 1-2 |
| HEX | HEX8, HEX20, HEX27, HEX64 | 1-3 |

### Hermite Elements (C¹ Continuity)

| Element | Bases | Derivatives | Use Case |
|---------|-------|-------------|----------|
| LINE2 | 4 | Value + slope | Euler-Bernoulli beams (cubic) |
| QUAD4 | 16 | Value + ∂/∂x + ∂/∂y + ∂²/∂x∂y | Kirchhoff plates (bicubic) |

---

## Integration Schemes

| Scheme | Description | Use Case |
|--------|-------------|----------|
| **GAUSS** | Gauss-Legendre quadrature | Default, optimal for polynomials |

**Auto order selection** (`aIntegrationOrder = 0`):
- Automatically selects sufficient order based on element type and interpolation order
- Conservative choice (safe but not minimal)

---

## Performance Considerations

### Critical: Cache IntegrationData Outside Loops

```cpp
// BAD - Performance disaster
for (Element* e : elements) {
    IntegrationData data(e->type());  // Allocation + factory per element!
    data.populate(4);
    // ... assembly
}

// GOOD - Cache by element type
Map<ElementType, IntegrationData*> cachedData;
for (Element* e : elements) {
    if (!cachedData.contains(e->type())) {
        cachedData[e->type()] = new IntegrationData(e->type());
        cachedData[e->type()]->populate(4);
    }
    IntegrationData* data = cachedData[e->type()];
    // ... assembly (reuse precomputed values)
}
```

| Operation | Without Cache | With IntegrationData | Speedup |
|-----------|---------------|----------------------|---------|
| Shape eval | O(num_gp × num_basis²) | O(1) per point | ~100× |
| Derivatives | O(num_gp × num_basis² × ndim) | O(1) per point | ~100-500× |

---

## Source Code

**Module location:** `../../`

**Key source files:**
- **Base classes**: `cl_IF_InterpolationFunction.hpp`, `cl_IF_InterpolationFunctionTemplate.hpp`
- **Integration**: `cl_IF_IntegrationData.{hpp,cpp}`, `fn_IF_initialize_integration_points.{hpp,cpp}`
- **Factories**: `cl_IF_InterpolationFunctionFactory.{hpp,cpp}`, `cl_EdgeFunctionFactory.{hpp,cpp}`
- **Lagrange shapes**: `lagrange/` subdirectory
- **Hermite shapes**: `hermite/` subdirectory
- **Bernstein shapes**: `bernstein/` subdirectory
- **Bubble functions**: `bubble/` subdirectory — **dormant by design.** The bubble
  machinery (14 headers, `InterpolationFunctionFactory::create_bubble_function`, four
  call sites) hangs off `MaxwellFactory::mUseEnrichment`, which is never set from the
  deck: it is the residue of a dropped experiment to enrich the phi elements of the
  iron, where the enrichment space was the wrong one — bubble rather than hierarchical
  (Dular et al. 2021). It is deliberately kept, not deleted: the machinery is working
  scaffolding for a future hierarchical enrichment and could be compile-gated (e.g. a
  `BELFEM_BUBBLE` define) should it ever get in the way. Do not wire a consumer to it
  without revisiting the enrichment-space choice.
- **Nédélec elements**: `nedelec/` subdirectory

---

## External References

### Shape Functions and Integration

- **Zienkiewicz & Taylor**, "The Finite Element Method" Vol. 1, Ch. 6-7 (Shape functions and integration)
- **Hughes**, "The Finite Element Method", Ch. 3 (Isoparametric elements and integration)
- **Bathe**, "Finite Element Procedures", Ch. 5 (Element formulation)

### Nédélec Elements

- **Nédélec (1980)**, "Mixed finite elements in ℝ³", Numer. Math. 35, 315-341
- **Monk (2003)**, "Finite Element Methods for Maxwell's Equations"
- See `literature/papers/fem` for BELFEM-specific electromagnetics applications

### Related BELFEM Modules

- **Mesh** (`src/mesh/`): Element definitions, facet orientations
- **FEM Kernel** (`src/fem/kernel/`): Uses IntegrationData for assembly
- **Linear Algebra** (`src/linalg/`): Matrix/Vector for shape function storage
- **Numerics** (`src/numerics/`): Parent module

---

## Common Pitfalls

1. **Creating factories inside element loops** → Cache IntegrationData by element type
2. **Calling d2NdXi2 on an unspecialized template** → every factory-created element implements it; the BELFEM_ERROR in InterpolationFunctionTemplate fires only for a combination without a specialization
3. **Ignoring facet orientation** → Slave facets require correct orientation for consistent normals
4. **Auto order = 0** → Understand what "auto" selects for your element (see usage guide)
5. **Ownership confusion** → IntegrationData can own or borrow the shape function (see constructors)

---

## Development Notes

### Adding New Element Types

When adding support for new element geometries:

1. Create shape function class in appropriate subdirectory (lagrange/, hermite/, etc.)
2. Specialize `InterpolationFunctionTemplate<G,T,D,B>` or derive from `InterpolationFunction`
3. Add case to `InterpolationFunctionFactory::create_*_function()`
4. Add integration point generation to `fn_IF_initialize_integration_points.cpp`
5. Update element support matrix in documentation

### Thread Safety

- **InterpolationFunction**: Thread-safe for read-only operations (N, dNdXi evaluation)
- **IntegrationData**: Thread-safe after `populate()` (read-only access)
- **Factories**: Thread-safe (stateless)

### MPI Compatibility

- Each rank creates independent instances
- IntegrationData is local (no communication needed)
- Edge functions work with distributed meshes (via element linking)

---

## See Also

- **Project README**: `../../../README.md`
- **Claude Instructions**: `../../../CLAUDE.md`
- **Documentation Guidelines**: `../../../doc/documentation_guidelines.md`
- **FEM Kernel**: `../kernel/doc/`
- **Mesh Module**: `../../mesh/doc/`
