# Hanging DOFs and Static Condensation via T-Matrices {#fem_kernel_hanging_dofs_static_condensation}

**Date:** 2026-01-30
**Module:** src/fem/kernel
**Purpose:** Comprehensive documentation on hanging DOF implementation and static condensation in BELFEM

---

## Overview

**Hanging DOFs** are degrees of freedom that are constrained to be linear combinations of other (source/master) DOFs. They are eliminated from the global system via **static condensation** using **T-matrices** (transformation matrices).

**Key terminology:**
- **Hanging DOF**: A DOF that depends linearly on other DOFs (eliminated from system)
- **Hanging edge/node/face**: A mesh entity with hanging DOFs associated with it
- **Source/Master DOF**: Independent DOF that hanging DOFs depend on
- **T-matrix**: Transformation matrix implementing the change of basis

**Current applications in BELFEM:**
1. **H-φ interfaces** - Thin-shell conductor edges → air/ferro nodes (edge-to-node coupling)
2. **H-H interfaces** - Thin-shell conductor edges → volume conductor edges (edge-to-edge coupling)
3. **Periodic boundary conditions** - Slave DOFs constrained to master DOFs on opposite boundary

**Future applications:**
- **Adaptive mesh refinement** - Hanging nodes at refinement boundaries
- **Multi-point constraints** - General linear constraints between DOFs

---

## Table of Contents

1. [Hanging DOF Concept and Use Cases](#hanging-dof-concept-and-use-cases)
2. [Physics Motivation: H-φ Mixed Formulation](#physics-motivation-h-φ-mixed-formulation)
3. [Why Static Condensation Over Lagrange Multipliers](#why-static-condensation-over-lagrange-multipliers)
4. [T-Matrix Mathematical Formulation](#t-matrix-mathematical-formulation)
5. [Interface Coupling Types](#interface-coupling-types)
6. [Implementation Architecture](#implementation-architecture)
7. [Code Walkthrough](#code-walkthrough)
8. [Known Issues and Debugging](#known-issues-and-debugging)
9. [Literature References](#literature-references)

---

## Hanging DOF Concept and Use Cases

### What are Hanging DOFs?

A **hanging DOF** is a degree of freedom whose value is **completely determined** by other DOFs through a linear constraint:

```
u_hanging = T × u_sources
```

Where:
- `u_hanging` is a vector of hanging DOFs (eliminated from the system)
- `u_sources` is a vector of independent source/master DOFs (retained in the system)
- `T` is the transformation matrix (typically sparse)

**Key property:** Hanging DOFs are **NOT part of the global system matrix** - they are eliminated during assembly via static condensation.

### Current Use Cases in BELFEM

#### 1. H-φ Interfaces (Edge-to-Node Coupling)

**Problem:** Incompatible discretizations at conductor-air interfaces

**Setup:**
- **Conductor domain (HTS tapes):** H-formulation with edge DOFs (Nédélec elements)
- **Air/ferro domain:** φ-formulation with node DOFs (Lagrange elements)
- **Interface:** Thin-shell 2D conductor embedded in 3D non-conducting volume

**Constraint:**
```
∫_edge H·dl ≈ (φ₁ - φ₀)    (tangential H continuity)
```

**Hanging relationship:**
- Thin-shell edge DOF (H-field) hangs on volume node DOFs (φ-field)
- 1 edge DOF → 2 node DOFs (LINE2) or 3 node DOFs (LINE3)

**Why this works:**
- Maintains tangential H continuity across interface
- Reduces DOFs significantly (2× in 2D, 3× in 3D using φ instead of H in air)
- Avoids saddle-point system from Lagrange multipliers

**Literature:** Arsenault et al. 2023, Messe et al. 2023 Section 2.2

---

#### 2. H-H Interfaces (Edge-to-Edge Coupling)

**Problem:** Orientation mismatch at thin-shell-to-conductor interfaces

**Setup:**
- **Thin-shell conductor (2D):** H-formulation with edge DOFs
- **Volume conductor (3D):** H-formulation with edge DOFs
- **Interface:** Thin-shell 2D conductor embedded in 3D conducting volume

**Constraint:**
```
∫_shell H·dl = ± ∫_volume H·dl    (sign depends on edge orientation)
```

**Hanging relationship:**
- Thin-shell edge DOF hangs on volume edge DOF
- 1 edge DOF → 1 edge DOF with sign `±1.0`
- Sign determined by comparing node ordering (topology)

**Why this works:**
- Direct 1:1 coupling (simplest T-matrix possible)
- Handles orientation automatically via sign
- Enables HTS tapes in contact with copper stabilizers

**Literature:** Messe et al. 2023 Section 2.5; Monk 2003 Chapter 5 (Nédélec edge elements)

---


#### 4. Adaptive Mesh Refinement (Potential)

**Problem:** Hanging nodes at mesh refinement boundaries

**Setup:**
```
Coarse mesh:  [a]----------[b]
Refined mesh: [a]---[h]----[b]
```

**Hanging relationship:**
- Node `h` is hanging: `u_h = 0.5 * u_a + 0.5 * u_b`
- Maintains C⁰ continuity across refinement boundary

**Benefits:**
- Local mesh refinement without global remeshing
- Standard FEM technique (well-documented)

**Literature:** Zienkiewicz & Taylor Vol. 1, Chapter 9 (adaptive refinement)

---

## Physics Motivation: H-φ Mixed Formulation

### Why We Need Interface Coupling

BELFEM uses a **mixed H-φ formulation** for electromagnetic field analysis of superconducting systems:

**Conducting domains (e.g., HTS tapes, copper stabilizers):**
- **H-formulation**: Solve for magnetic field `H` directly using Nédélec edge elements
- DOF type: **Edge elements** (Nédélec, `H(curl)` space)
- Physics: `∂(μ₀H)/∂t + ∇×(ρ∇×H) = 0` (magnetodynamic)
- DOFs = number_of_edges per element

**Non-conducting domains (e.g., air, ferromagnetic cores):**
- **φ-formulation**: Solve for scalar magnetic potential `φ` with `H = -∇φ`
- DOF type: **Nodal elements** (Lagrange, `H¹` space)
- Physics: `∇·(μ∇φ) = 0` (magnetostatic)
- DOFs = number_of_nodes per element

**Why this matters:**
- DOF count reduction: 2× in 2D, 3× in 3D (scalar potential vs. vector field)
- Better conditioning: Scalar equations are easier to solve
- Physics-accurate: Reflects that `∇×H = J = 0` in non-conducting regions

### Interface Continuity Requirements

At conductor-air interfaces, tangential H must be continuous:
```
H_tangential(conductor) = H_tangential(air)
```

**Problem:** Incompatible discretizations
- Conductor side: H-field on **edges** (Nédélec, vector)
- Air side: φ-field on **nodes** (Lagrange, scalar)

**Solution:** Static condensation via T-matrices
- Eliminate edge DOFs on one side
- Express them as linear combinations of node DOFs on the other side
- Maintain tangential continuity

**Literature:** Messe et al. 2023, Section 2.2; Arsenault et al. 2023, Equations 5-7

---

## Why Static Condensation Over Lagrange Multipliers

### The Lagrange Multiplier Approach

**Standard approach for constraints:**
```
Minimize: E(u) = ½u^T K u - u^T f
Subject to: C u = 0
```

**Lagrange multiplier formulation:**
```
[  K   C^T ] [ u ]   [ f ]
[  C    0  ] [ λ ] = [ 0 ]
```

**Problem:** Zeros on the main diagonal!

This creates a **saddle-point system**:
- Not positive definite (eigenvalues have mixed signs)
- LU decomposition can fail (pivoting issues)
- Iterative solvers converge slowly (indefinite preconditioners)

**Real failure mode in MUMPS/STRUMPACK:**
```
ERROR: Zero or near-zero pivot encountered during factorization
ERROR: Matrix is singular or ill-conditioned
```

**Literature:** Messe et al. 2023, lines 509-517; Bathe 2016, §3.2.3

### The Static Condensation Alternative

**Key idea:** Eliminate dependent DOFs directly via change of basis

**Condensation formulation:**
```
u_local = T × u_global    (change of basis)
K_global = T^T × K_local × T    (transform stiffness)
f_global = T^T × f_local         (transform force)
```

**Result:**
- **No constraint equations** (DOFs eliminated before assembly)
- **No zeros on diagonal** (system remains positive definite)
- **Smaller system** (fewer DOFs to solve)
- **Robust solvers** (standard Cholesky/LU works)

**Trade-off:**
- Static condensation: More complex assembly, simpler solve
- Lagrange multipliers: Simpler assembly, more complex solve

**BELFEM choice:** Static condensation (prioritizes solver robustness)

**Literature:** Bathe 2016, §8.2 (static condensation); Messe et al. 2023, "We prefer static condensation over Lagrange multipliers to avoid zeros on the main diagonal"

---

## T-Matrix Mathematical Formulation

### Change of Basis Representation

**T-matrix represents the constraint:**
```
u_hanging = T × u_master
```

Where:
- `u_hanging` = local/dependent/slave DOFs (eliminated from system)
- `u_master` = global/independent/master DOFs (retained in system)
- `T` = transformation matrix (sparse, typically 1-3 nonzeros per row)

**Dimensions:**
```
T: (n_hanging × n_master)
```

**Example: Hanging node at edge midpoint**
```
Coarse edge:  [a]----------[b]
Refined edge: [a]---[h]----[b]

Constraint: u_h = 0.5 * u_a + 0.5 * u_b

T-matrix:
       [a   b]
  [h] [0.5 0.5]

u_hanging = [u_h]
u_master = [u_a, u_b]
```

### Assembly with T-Matrices

**Element-level assembly (local coordinates):**
```cpp
K_local: (n_hanging × n_hanging)  // Element stiffness
f_local: (n_hanging)               // Element force
```

**Transformation to global coordinates:**
```cpp
K_global = T^T × K_local × T     // Matrix: (n_master × n_master)
f_global = T^T × f_local          // Vector: (n_master)
```

**Implementation in BELFEM:**
```cpp
// src/fem/kernel/cl_FEM_DofManager.cpp:430-629
if( tElement->has_t_matrix() )
{
    mIWG->compute_jacobian( tElement, tJ );          // K_local
    tElement->t_matrix()->project( tJ, tTJT );      // K_global = T^T × K_local × T
    mSolverData->assemble_jacobian( tElement, tTJT);
}
```

### T-Matrix Storage Format

**BELFEM uses CSR (Compressed Sparse Row) format:**

```cpp
// src/fem/kernel/cl_FEM_Tmatrix.hpp:23-31
class Tmatrix
{
    const uint mNumRows ;      // n_hanging (local DOFs)
    const uint mNumCols  ;     // n_master (global DOFs)
    const uint mNumNonZeros  ; // Number of nonzero entries

    uint * mPointers = nullptr ; // Row pointers (size: mNumRows + 1)
    uint * mIndices  = nullptr ; // Column indices (size: mNumNonZeros)
    real * mValues   = nullptr ; // Nonzero values (size: mNumNonZeros)
};
```

**CSR Access Pattern:**
```cpp
// Row i has nonzeros at columns mIndices[j] with values mValues[j]
// for j in range [mPointers[i], mPointers[i+1])

for (uint i = 0; i < mNumRows; ++i) {
    for (uint j = mPointers[i]; j < mPointers[i+1]; ++j) {
        uint col = mIndices[j];
        real val = mValues[j];
        // T(i, col) = val
    }
}
```

### T-Matrix Projection Operations

**Vector projection: `b = T^T × a`**
```cpp
// src/fem/kernel/cl_FEM_Tmatrix.cpp:100-115
void Tmatrix::project( const Vector< real > & aA, Vector< real > & aB ) const
{
    aB.set_size( mNumCols, 0.0 );
    for ( uint i = 0; i < mNumRows; ++i)
    {
        for ( uint j = mPointers[ i ]; j < mPointers[ i + 1 ]; ++j)
        {
            aB( mIndices[ j ] ) += mValues[ j ] * aA( i );
        }
    }
}
```

**Matrix projection: `B = T^T × A × T`**
```cpp
// src/fem/kernel/cl_FEM_Tmatrix.cpp:117-153
void Tmatrix::project( const Matrix< real > & aA, Matrix< real > & aB ) const
{
    aB.set_size( mNumCols, mNumCols, 0.0 );
    for ( uint i = 0; i < mNumRows; ++i)
    {
        for ( uint p = mPointers[ i ]; p < mPointers[ i + 1 ]; ++p )
        {
            uint j = mIndices[ p ];
            real T_ji = mValues[ p ];        // T(i,j) → T^T(j,i)
            for( uint k=0; k<mNumRows; ++k )
            {
                for( uint q=mPointers[ k ]; q<mPointers[ k+1 ]; ++q )
                {
                    uint l = mIndices[ q ];
                    real T_kl = mValues[ q ];  // T(k,l)
                    aB( j, l ) += T_ji * aA( i, k ) * T_kl ;  // B(j,l) += T^T(j,i) A(i,k) T(k,l)
                }
            }
        }
    }
}
```

**Complexity:**
- Vector projection: `O(nnz)` where `nnz` = number of nonzeros in T
- Matrix projection: `O(nnz² × n_hanging)` (can be expensive for dense A)

---

## Interface Coupling Types

BELFEM handles two distinct coupling scenarios at H-φ interfaces:

### 1. Edge-to-Node Coupling (H-φ Interface)

**Use case:** Thin-shell conductor edge → Air/ferro nodes

**Physics:**
- Conductor side: H-field on edge (circulation DOF)
- Air side: φ-field on nodes (scalar potential)
- Coupling: `∫_edge H·dl ≈ (φ₁ - φ₀) / length`

**Geometry:**
```
Thin-shell (conductor, H-formulation):
  [n₀]--------[edge DOF]-------[n₁]

Volume (air/ferro, φ-formulation):
  [φ₀]------------------------[φ₁]
```

**T-matrix structure:**
```
       [φ₀   φ₁]
[edge] [1.0  -1.0]   (for LINE2 edge)
```

**For LINE3 edge (quadratic):**
```
         [φ₀    φ₁   φ₂]
[edge₀]  [1.0   1/3  -4/3]   (first DOF)
[edge₁]  [1/3   1.0  -4/3]   (second DOF)
```

**Implementation:**
- `MaxwellFactory::hang_thinshell_edges_on_nodes_bottom()` (lines 1251-1317)
- `MaxwellFactory::hang_thinshell_edges_on_nodes_top()` (lines 1320-1387)
- `DofData::create_dofwise_t_matrices_master()` (lines 3572-3686)

**Literature:** Arsenault et al. 2023, Equation 8; Alves et al. 2022b, Section 3.2

### 2. Edge-to-Edge Coupling (H-H Interface)

**Use case:** Thin-shell conductor edge → Volume conductor edge

**Physics:**
- Both sides: H-formulation (edge elements)
- Coupling: Direct 1:1 correspondence with orientation
- Sign: `h_shell = ± h_conductor` (depends on edge direction)

**Geometry:**
```
Thin-shell (conductor, 2D):
  [n₀]--------[edge_shell]--------[n₁]

Volume (conductor, 3D):
  [n₀]--------[edge_volume]-------[n₁]
     (or possibly [n₁]------------[n₀] if reversed)
```

**T-matrix structure:**
```
            [edge_volume]
[edge_shell] [+1.0]         if same orientation (n₀→n₁ matches)
[edge_shell] [-1.0]         if reversed orientation (n₀→n₁ opposite)
```

**Sign determination (by node ordering):**
```cpp
// Compare node IDs to determine orientation
if ( edge_shell->node(0)->id() == edge_volume->node(0)->id() &&
     edge_shell->node(1)->id() == edge_volume->node(1)->id() )
{
    sign = +1.0;  // Same orientation
}
else if ( edge_shell->node(0)->id() == edge_volume->node(1)->id() &&
          edge_shell->node(1)->id() == edge_volume->node(0)->id() )
{
    sign = -1.0;  // Reversed orientation
}
else
{
    ERROR("Edge node IDs don't match!");
}
```

**Why sign matters (Nédélec edge elements):**

Edge DOF represents **circulation** (line integral):
```
DOF_edge = ∫_edge H·dl
```

Direction matters:
- Traversing edge from n₀→n₁: `∫_edge H·dl`
- Traversing edge from n₁→n₀: `-∫_edge H·dl`

**Implementation:**
- `MaxwellFactory::hang_thinshell_edges_on_edges_bottom()` (lines 1390-1444)
- `MaxwellFactory::hang_thinshell_edges_on_edges_top()` (lines 1447-1501)
- `DofData::create_dofwise_t_matrices_master()` (lines 3687-3741)

**Literature:** Messe et al. 2023, Section 2.5; Monk 2003, Chapter 5 (Nédélec elements)

---

### 3. Higher-Order Elements: The Mesh→DOF Complexity

**CRITICAL INSIGHT:** Mesh-level coupling (edges hang on edges) ≠ DOF-level mapping!

**The two-phase process:**

**Phase 1 - Mesh level (MaxwellFactory):**
```cpp
edge_shell hangs on edge_volume with sign = -1.0 (reversed orientation)
```
**What you know:** Edge relationship and overall sign
**What you DON'T know:** How many DOFs? How do individual DOFs map?

**Phase 2 - DOF level (DofData):**
Now you must determine per-DOF mapping based on element order:

---

#### Linear Elements (LINE2, 1 DOF per edge)

**Simple case:**
```
edge_volume:  [n₀]---(DOF 0)---[n₁]
edge_shell:   [n₀]---(DOF 0')--[n₁]  (same orientation)

Mapping:
  shell_dof(0) → +1.0 * volume_dof(0)
```

**Reversed:**
```
edge_volume:  [n₀]---(DOF 0)---[n₁]
edge_shell:   [n₁]---(DOF 0')--[n₀]  (reversed)

Mapping:
  shell_dof(0) → -1.0 * volume_dof(0)
```

**Simple:** 1 DOF → mesh sign directly transfers to DOF sign

---

#### Quadratic Elements (LINE3, 2 DOFs per edge)

**Complexity:** DOF ordering depends on traversal direction!

**Same orientation (+1.0):**
```
edge_volume:  [n₀]---(DOF 0)---(DOF 1)---[n₁]
edge_shell:   [n₀]---(DOF 0')---(DOF 1')--[n₁]

Mapping:
  shell_dof(0) → +1.0 * volume_dof(0)  ← DOF 0 to DOF 0
  shell_dof(1) → +1.0 * volume_dof(1)  ← DOF 1 to DOF 1
```

**Reversed orientation (-1.0):**
```
edge_volume:  [n₀]---(DOF 0)---(DOF 1)---[n₁]
edge_shell:   [n₁]---(DOF 0')---(DOF 1')--[n₀]  ← REVERSED TRAVERSAL!

Mapping:
  shell_dof(0) → -1.0 * volume_dof(1)  ← DOF 0' maps to DOF 1 (swapped!)
  shell_dof(1) → -1.0 * volume_dof(0)  ← DOF 1' maps to DOF 0 (swapped!)
```

**Why order swaps:** Nédélec edge DOFs represent circulation in specific direction
- `DOF 0` represents `∫_{n₀→center} H·dl`
- `DOF 1` represents `∫_{center→n₁} H·dl`
- When traversing backward (n₁→n₀), `DOF 0'` corresponds to volume `DOF 1`, not `DOF 0`!

**T-matrix for reversed quadratic edge:**
```
T (2×2):
         [vol_dof_0  vol_dof_1]
[shell_dof_0] [   0       -1.0  ]  ← maps to vol_dof_1 with sign
[shell_dof_1] [ -1.0        0   ]  ← maps to vol_dof_0 with sign
```

**Implementation (current code lines 3841-3866):**
```cpp
if ( tW == 1.0 )
{
    // Same orientation: preserve DOF order
    shell_dof(0)->set_sources( volume_dof(0), +1.0 );
    shell_dof(1)->set_sources( volume_dof(1), +1.0 );
}
else if ( tW == -1.0 )
{
    // Reversed: swap DOF order AND apply sign
    shell_dof(0)->set_sources( volume_dof(1), -1.0 );  ← Swapped!
    shell_dof(1)->set_sources( volume_dof(0), -1.0 );  ← Swapped!
}
```

**Key insight:** **Both sign AND order must change** for reversed higher-order elements!

**Status:** ⚠️ **UNTESTED** - Code has logic for quadratic edges but hasn't been validated

**Literature:** Monk 2003, §5.3 (edge element orientation); Nédélec 1980 (H(curl) finite elements)

---

## Implementation Architecture

### Three-Phase Workflow

**Phase 1: Mesh-level hanging relationships (MaxwellFactory)**
- Identify which thin-shell edges/facets should be "hung"
- Populate `mesh::Basis::mSources` and `mesh::Basis::mWeights` containers
- Store topological relationships (which edges hang on which nodes/edges)

**Phase 2: DOF-level T-matrix construction (DofData)**
- Convert mesh-level relationships to DOF-level dependencies
- Build T-matrix for each hanging DOF
- Handle multi-level dependencies (cascading hanging)

**Phase 3: Element assembly (DofManager)**
- Project element matrices through T-matrix: `K_global = T^T × K_local × T`
- Assemble transformed contributions to global system
- Hanging DOFs never appear in system matrix (eliminated)

### Data Flow Diagram

```
MaxwellFactory::create_hanging_edges_and_facets()
  ↓
  ├─→ hang_thinshell_edges_on_nodes_bottom()
  │     ↓
  │     For each thin-shell edge at interface:
  │       edge->add_source(node0, weight0)  ← mesh::Basis container
  │       edge->add_source(node1, weight1)
  │
  ├─→ hang_thinshell_edges_on_edges_bottom()
  │     ↓
  │     For each thin-shell edge at interface:
  │       edge->add_source(volume_edge, ±1.0)  ← orientation sign!
  │
DofData::create_dofwise_t_matrices_master()
  ↓
  For each hanging DOF:
    ├─→ Read sources from dof->mesh_basis()->source(i)
    ├─→ Read weights from dof->mesh_basis()->weight(i)
    ├─→ Construct T-matrix (CSR format)
    └─→ Store in dof->set_sources(sources, weights)
  ↓
Element::relink_dofs()
  ↓
  Replace hanging DOFs with T-matrix transformation
  ↓
DofManager::compute_jacobian_and_rhs()
  ↓
  For each element:
    if (element->has_t_matrix()):
      IWG->compute_jacobian(element, K_local)
      element->t_matrix()->project(K_local, K_global)  ← T^T × K × T
      SolverData->assemble_jacobian(element, K_global)
```

### Container Hierarchy and the Mesh → DOF Distinction

**CRITICAL:** Hanging relationships exist at **two separate levels** with different information:

#### Level 1: Mesh-Level Relationships (Geometric)

**What you know:**
- Geometric coupling: "This edge hangs on these nodes/edges"
- Entity relationships only

**What you DON'T know yet:**
- How many DOFs does each entity have?
- How do individual DOFs map?
- What are the signs for each DOF?

**Storage:**
```cpp
// src/mesh/cl_Mesh_Basis.hpp
class Basis {
    uint mNumberOfSources;
    Basis ** mSources;  // Pointers to source entities (nodes/edges/faces)
    real * mWeights;    // Entity-level weights (e.g., orientation sign for edges)
};
```

**Example - Edge hanging on edge:**
```cpp
edge_shell->add_source(edge_volume, -1.0);  // Mesh-level: reversed orientation
// But we still don't know:
// - Does edge_volume have 1 DOF (linear) or 2 DOFs (quadratic)?
// - How do individual DOFs map?
```

---

#### Level 2: DOF-Level Mapping (Per-DOF)

**Now you know:**
- Exactly which DOF maps to which DOF
- Signs for each individual DOF
- Complete T-matrix entries

**Storage:**
```cpp
// src/fem/kernel/cl_FEM_Dof.hpp
class Dof {
    uint mNumberOfSources;
    Dof ** mSources;      // Pointers to individual source DOFs
    real * mCoefficients; // Per-DOF T-matrix weights
};
```

**Example - Quadratic edge (2 DOFs per edge):**
```cpp
// Mesh-level told us: edge_shell hangs on edge_volume with sign -1.0

// DOF-level mapping (LINEAR element, 1 DOF per edge):
dof_shell(0) -> -1.0 * dof_volume(0)

// DOF-level mapping (QUADRATIC element, 2 DOFs per edge):
dof_shell(0) -> -1.0 * dof_volume(1)  // Reversed order!
dof_shell(1) -> -1.0 * dof_volume(0)  // Reversed order!
```

**Why the order swap?** Nédélec edge elements traverse edge in opposite direction when reversed.

---

#### Level 3: Element-Level T-Matrix (Assembly)

**Complete transformation for element:**
```cpp
// src/fem/kernel/cl_FEM_Tmatrix.hpp
class Tmatrix {
    uint mNumRows;       // Number of hanging DOFs on this element
    uint mNumCols;       // Number of source DOFs
    uint mNumNonZeros;   // Nonzeros in sparse T-matrix

    uint * mPointers;    // CSR row pointers
    uint * mIndices;     // CSR column indices
    real * mValues;      // CSR nonzero values (individual DOF weights)
};
```

**Example - Element with 1 quadratic hanging edge:**
```
T-matrix (2×2):
       [dof_vol_0  dof_vol_1]
[dof_shell_0] [   0      -1.0  ]    (reversed: maps to dof_vol_1 with sign)
[dof_shell_1] [ -1.0       0   ]    (reversed: maps to dof_vol_0 with sign)
```

---

### The Vertex Container "Hack" (Visualization)

**Important clarification:** The vertex container usage is **intentional, not a bug!**

**Why it exists:**
- At interfaces, nodes are **duplicated and uncoupled** for better Paraview visualization
- The nodes edges link to for reference ≠ the geometric nodes the edge is connected to
- This allows discontinuous fields to be visualized correctly at interfaces

**Two separate node sets:**
1. **Geometric nodes** (`edge->node(i)`): Physical connectivity in mesh
2. **Reference nodes** (`edge->vertex(i)`): Nodes used for visualization/reference

**Why both containers are needed:**
- **Source container** (`mSources`): Geometric coupling (which entity hangs on which)
- **Vertex container** (`mVertices`): Visualization nodes (for field output)

**This is a feature, not a bug** - it enables proper discontinuous field visualization at H-φ and H-H interfaces!

---

## Code Walkthrough

### 1. MaxwellFactory: Mesh-Level Hanging Setup

**File:** `src/fem/maxwell/cl_MaxwellFactory.cpp`

**Function:** `create_hanging_edges_and_facets()`

**Purpose:** Top-level orchestration for interface and thin-shell coupling

The function refuses higher-order meshes
(`BELFEM_ERROR( mMesh->max_element_order() == 1, "Not implemented for higher order" )`),
then walks `mMesh->sidesets()`: `InterfaceCondAir` sidesets (and `InterfaceCondFerro` under
`Formulation::HPhi`) get their master-facet edges hung on the slave-side nodes. Thin-shell
layer facets are dispatched per facet according to the `DomainType` of the neighboring
block — φ-side facets first, then conductor-side ones:

```cpp
// bottom layer of the shell; the top layer is handled the same way
// with the slave block and hang_thinshell_edges_on_*_top()
for ( index_t k : tPhiFacets )          // neighbour block is not a Conductor
{
    this->hang_thinshell_edges_on_nodes_bottom( tWork, tFacets( k ), tElementsBottom( k ) );
}
for ( index_t k : tConductorFacets )    // neighbour block is a Conductor
{
    this->hang_thinshell_edges_on_edges_bottom( tWork, tFacets( k ), tElementsBottom( k ) );
}
```

**Key design choice:**
- φ neighbor (air, ferro) → Edge-to-node (H-φ interface)
- `DomainType::Conductor` neighbor → Edge-to-edge (H-H interface)

---

**Function:** `hang_thinshell_edges_on_nodes_bottom( EdgeWorkData & aWork, mesh::Facet * aFacet, mesh::Element * aElement )`

**Purpose:** Populate edge→node dependencies for one bottom-layer facet

`mesh::get_bottom_nodes()` / `mesh::get_bottom_edges()` collect the shell-side nodes and
edges of `aElement`; the volume nodes come from `aFacet->master()->get_nodes_of_facet()`.
Each unprocessed shell edge gets a source container sized to its node count and one volume
node per edge node:

```cpp
tEdge->allocate_source_container( tEdge->number_of_nodes() );

for ( uint k=0; k<tEdge->number_of_nodes(); ++k )
{
    mesh::Node * tNode = tNodesOnVolume( tEdge->node( k )->index() );

    // note: we determine the weights later in
    // DofData::create_dofwise_t_matrices_master()
    tEdge->add_source( tNode );
}

tEdge->flag();  // tag edge as processed
```

An edge that already hangs on a periodic master edge keeps that tie and is only flagged.

---

**Function:** `hang_thinshell_edges_on_edges_bottom( EdgeWorkData & aWork, mesh::Facet * aFacet, mesh::Element * aElement )`

**Purpose:** Populate edge→edge dependencies for one bottom-layer facet

Shell-side nodes and edges are collected as above; the volume edges come from
`aFacet->master()->get_edges_of_facet()`. Each unprocessed shell edge gets a single source,
the matching volume edge, with the orientation sign as weight:

```cpp
tEdge->allocate_source_container( 1 );

mesh::Edge * tOther = tEdgesOnVolume( tEdge->index() );

if ( tEdge->node( 0 )->index() == tOther->node( 0 )->index() &&
     tEdge->node( 1 )->index() == tOther->node( 1 )->index() )
{
    tWeight = 1.0 ;     // same orientation
}
else if ( tEdge->node( 0 )->index() == tOther->node( 1 )->index() &&
     tEdge->node( 1 )->index() == tOther->node( 0 )->index() )
{
    tWeight = -1.0 ;    // reversed orientation
}
else
{
    BELFEM_ERROR( false, "Could not determine edge orientation" );
}
tEdge->add_source( tOther, tWeight );

tEdge->flag();
```

**Design features:**
- **Sign determination:** Uses **node index** (topology), not geometric tangent vectors
- **Periodic guard:** an edge that already hangs on a periodic master edge keeps that tie

---

### 2. DofData: DOF-Level T-Matrix Construction

**File:** `src/fem/kernel/cl_FEM_DofMgr_DofData.cpp`

**Function:** `create_dofwise_t_matrices_master()` (lines 3486-3784)

**Purpose:** Convert mesh-level hanging to DOF-level T-matrices

**Structure:**
```cpp
void DofData::create_dofwise_t_matrices_master()
{
    Cell< Dof * > tSources;
    Vector< real > tWeights;
    Vector< real > tCoefficients;
    Vector< index_t > tIndices;
    Matrix< real > tNodeWeights;

    // Connect DOFs to mesh entities
    this->connect_dofs_to_mesh();

    // Process each hanging DOF
    for ( Dof * tDof : mHangingDOFs )
    {
        // Branch based on entity type
        if (tDof->entity_type() == EntityType::NODE) {
            // Handle hanging nodes
        }
        else if (tDof->entity_type() == EntityType::EDGE) {
            // Handle hanging edges (TWO CASES!)
            if (source is NODE) {
                // Edge-to-node coupling (H-φ interface)
            }
            else if (source is EDGE) {
                // Edge-to-edge coupling (H-H interface)
            }
        }
        else if (tDof->entity_type() == EntityType::FACE) {
            // Handle hanging faces
        }
    }
}
```

---

**Case 1: Trivial 1:1 Dependency** (lines 3505-3520)

**When:** Single source, same entity type, single DOF

```cpp
// Check for trivial dependency
if ( tDof->number_of_sources() == 1 &&
     tDof->source( 0 )->entity_type() == tDof->entity_type() &&
     tDof->source( 0 )->mesh_basis()->number_of_dofs() == 1 )
{
    BELFEM_ASSERT( ! tDof->source( 0 )->mesh_basis()->is_hanging(),
                  "Source of dof is hanging. This should not happen" );

    // Grab other dof
    Dof * tSource = reinterpret_cast< Dof * >( tDof->source( 0 )->mesh_basis()->dof( 0 ) );

    tSources.set_size( 1, tSource );
    tWeights.set_size( 1, 1.0 );  // ⚠️ BUG! Should use stored weight!

    tDof->set_sources( tSources, tWeights );
    continue;
}
```

**FIXED (historical):** an earlier revision hardcoded the trivial-dependency
weight to `1.0`, so edge-to-edge couplings picked up the wrong sign on
reversed edges. The current code passes the basis weight —
`tDof->set_source( ..., tDof->mesh_basis()->weight( 0 ) )` in
`cl_FEM_DofMgr_DofData.cpp` — which is exactly the correction this section
used to demand.

---

**Case 2: Hanging Nodes** (lines 3523-3570)

**When:** Node depends on other nodes (AMR, mesh refinement)

```cpp
if ( tDof->entity_type() == EntityType::NODE )
{
    tSources.set_size( tDof->mesh_basis()->number_of_sources(), nullptr );
    tWeights.set_size( tDof->mesh_basis()->number_of_sources(), 0.0 );

    uint t = tDof->type_id();

    for ( uint k = 0; k < tDof->mesh_basis()->number_of_sources(); ++k )
    {
        // Find matching DOF type
        for ( uint j = 0; j < tDof->mesh_basis()->source( k )->number_of_dofs(); ++j )
        {
            Dof * tSource = reinterpret_cast< Dof * >(
                tDof->mesh_basis()->source( k )->dof( j ) );

            if ( tSource->type_id() == t )
            {
                tSources( k ) = tSource;
                tWeights( k ) = tDof->mesh_basis()->weight( k );  // ✅ Uses stored weight
                break;
            }
        }
    }

    tDof->set_sources( tSources, tWeights );
}
```

**This case works correctly** - uses stored weights from mesh

---

**Case 3: Hanging Edges on Nodes (H-φ Interface)** (lines 3572-3686)

**When:** Edge DOF hangs on node DOFs (thin-shell H → air/ferro φ)

```cpp
if ( tDof->entity_type() == EntityType::EDGE )
{
    mesh::Edge * tEdge = tDof->edge();

    mesh::Basis * tSource = reinterpret_cast< mesh::Basis * >( tEdge->vertex( 0 ) );
    // Uses vertex container (reference nodes for visualization)

    if ( tSource->entity_type() == EntityType::NODE )
    {
        // Collect all source DOFs from nodes
        for ( uint i = 0; i < tEdge->number_of_nodes(); ++i )
        {
            mesh::Node * tNode = reinterpret_cast< mesh::Node * >( tEdge->vertex( i ) );
            // Uses vertex container (may be different from geometric nodes at interface)

            if ( tNode->is_hanging() )
            {
                // Node itself is hanging - cascade dependencies
                for ( uint j=0; j < tNode->number_of_sources(); ++j )
                {
                    tSources.push( reinterpret_cast< Dof * >( tNode->source( j )->dof( 0 ) ) );
                }
            }
            else
            {
                // Node is independent
                tSources.push( reinterpret_cast< Dof * >( tNode->dof( 0 ) ) );
            }
        }

        unique( tSources );  // Remove duplicates

        // Build node weight matrix
        uint n = tSources.size();
        tNodeWeights.set_size( n, tEdge->number_of_nodes(), 0.0 );

        for ( uint i=0; i<tEdge->number_of_nodes(); ++i )
        {
            mesh::Node * tNode = reinterpret_cast< mesh::Node * >( tEdge->vertex( i ) );
            // Uses vertex container (reference nodes, may be duplicated at interface)

            if ( tNode->is_hanging() )
            {
                for ( uint k=0; k<tNode->number_of_sources(); ++k )
                {
                    tNodeWeights( tNode->source( k )->dof( 0 )->index(), i ) = tNode->weight( k );
                }
            }
            else
            {
                tNodeWeights( tNode->dof( 0 )->index(), i ) = 1.0 ;
            }
        }

        // Compute edge DOF weights from node weights
        if ( tEdge->number_of_nodes() == 2 )
        {
            // LINE2: h_edge = φ₁ - φ₀
            tCoefficients.set_size( 2, 0.0 );
            tCoefficients( 0 ) = 1.0 ;
            tCoefficients( 1 ) = -1.0 ;

            tWeights = tNodeWeights * tCoefficients;
            tDof->set_sources( tSources, tWeights );
        }
        else if ( tEdge->number_of_nodes() == 3 )
        {
            // LINE3: Two edge DOFs (quadratic)
            tCoefficients.set_size( 3, 0.0 );

            // First DOF
            tCoefficients( 0 ) =  1.0 ;
            tCoefficients( 1 ) =  1./3. ;
            tCoefficients( 2 ) = -4./3. ;
            tWeights = tNodeWeights * tCoefficients;
            reinterpret_cast< Dof * >( tEdge->dof( 0 ) )->set_sources( tSources, tWeights );

            // Second DOF
            tCoefficients( 0 ) =  1./3. ;
            tCoefficients( 1 ) =  1.0 ;
            tCoefficients( 2 ) = -4./3. ;
            tWeights = tNodeWeights * tCoefficients;
            reinterpret_cast< Dof * >( tEdge->dof( 1 ) )->set_sources( tSources, tWeights );
        }
    }

    tEdge->reset_vertex_container();  // ✅ Cleanup after use
    tEdge->flag();
}
```

**Design notes:**
1. Uses `vertex()` container for reference nodes (visualization) - requires `reinterpret_cast` by design
2. Vertex container cleaned up after T-matrix construction
3. Quadratic edge formulas (LINE3) are specific to H-φ coupling

---

**Case 4: Hanging Edges on Edges (H-H Interface)** (lines 3687-3741)

**When:** Edge DOF hangs on edge DOF (thin-shell H → conductor H)

```cpp
else if ( tSource->entity_type() == EntityType::EDGE )
{
    mesh::Edge * tOther = reinterpret_cast< mesh::Edge * >( tSource );

    if ( tEdge->number_of_nodes() == 2 )
    {
        Dof * tSourceDof = reinterpret_cast< Dof * >( tOther->dof( 0 ) );

        if ( tSourceDof->is_hanging() )
        {
            // Cascade: source edge is itself hanging
            tSources.set_size( tSourceDof->number_of_sources(), nullptr );
            tWeights.set_size( tSourceDof->number_of_sources(), tEdge->weight( 0 ) );

            for ( uint k=0; k<tSourceDof->number_of_sources(); ++k )
            {
                tSources( k ) = tSourceDof->source( k );
                tWeights( k ) *= tSourceDof->weight( k );  // ✅ Multiply signs!
            }
        }
        else
        {
            // Direct 1:1 coupling
            tSources.set_size( 1, tSourceDof );
            tWeights.set_size( 1, tEdge->weight( 0 ) );  // ✅ Uses stored sign!
        }

        reinterpret_cast< Dof * >(tEdge->dof( 0 ))->set_sources( tSources, tWeights );
    }
    else if ( tEdge->number_of_nodes() == 3 )
    {
        // LINE3 edge: Two DOFs (quadratic)
        // ⚠️ WARNING: this is not tested

        tSources.set_size( 1, nullptr );
        real tW = tEdge->weight( 0 );
        tWeights.set_size( 1, tW );

        if ( tW == 1.0 )
        {
            // Same orientation
            tSources( 0 ) = reinterpret_cast< Dof * >(tOther->dof( 0 ));
            reinterpret_cast< Dof * >(tEdge->dof( 0 ))->set_sources( tSources, tWeights );
            tSources( 0 ) = reinterpret_cast< Dof * >(tOther->dof( 1 ));
            reinterpret_cast< Dof * >(tEdge->dof( 1 ))->set_sources( tSources, tWeights );
        }
        else if ( tW == -1.0 )
        {
            // Reversed orientation - swap DOF order
            tSources( 0 ) = reinterpret_cast< Dof * >(tOther->dof( 1 ));
            reinterpret_cast< Dof * >(tEdge->dof( 0 ))->set_sources( tSources, tWeights );
            tSources( 0 ) = reinterpret_cast< Dof * >(tOther->dof( 0 ));
            reinterpret_cast< Dof * >(tEdge->dof( 1 ))->set_sources( tSources, tWeights );
        }
    }
}
```

**Key features:**
1. ✅ Uses stored sign from `tEdge->weight(0)`
2. ✅ Handles cascading dependencies (edge hangs on hanging edge)
3. ✅ Multiplies signs correctly when cascading
4. ⚠️ Quadratic edge case (LINE3) is untested

---

### 3. DofManager: Element Assembly with T-Matrices

**File:** `src/fem/kernel/cl_FEM_DofManager.cpp`

**Function:** `compute_jacobian()` (lines 430-629)

**Purpose:** Shows how T-matrices are used during assembly

```cpp
void DofManager::compute_jacobian( const bool aReset )
{
    // ... initialization ...

    // Loop over all elements
    for ( Element * tElement : mBlock->elements() )
    {
        // Check if element has T-matrix (hanging DOFs)
        if( tElement->has_t_matrix() )
        {
            // Compute element Jacobian in local coordinates
            mIWG->compute_jacobian( tElement, tJ );

            // Project through T-matrix: J_global = T^T × J_local × T
            tElement->t_matrix()->project( tJ, tTJT );

            // Assemble to global system
            mSolverData->assemble_jacobian( tElement, tTJT);
        }
        else
        {
            // No T-matrix: direct assembly
            mIWG->compute_jacobian( tElement, tJ );
            mSolverData->assemble_jacobian( tElement, tJ );
        }
    }
}
```

**Key point:** Assembly is **transparent** to the IWG (physics layer)
- IWG always computes in local coordinates (doesn't know about hanging DOFs)
- T-matrix projection happens in kernel layer
- Global system never sees hanging DOFs

---

## Periodic DOF Constraints

**Files:** `src/mesh/cl_Mesh_Periodicity.{hpp,cpp}`, `src/fem/kernel/cl_FEM_DofMgr_DofData.cpp`

Periodicity is imposed at mesh level, not by a separate DOF-level pass.
`mesh::Periodicity::set_entity_dependencies()` (called from `MaxwellFactory` before the
hanging setup) adds each slave node/edge/face's master as its single source with weight 1.0;
if the master is itself hanging, the slave takes over the master's sources and weights
instead. Matching is performed by `mesh::Periodicity` (see `src/mesh/doc/periodicity.md`).

`DofData::create_dofwise_t_matrices_master()` then treats those dependencies exactly like
interface hanging entities. A slave edge whose master is itself hanging is flattened onto the
master's sources in the deferred-chain pass at the end of that routine, so the T-matrix
assembly eliminates both interface-hanging and periodic-hanging DOFs.

---

## Known Issues and Debugging

### Critical Bugs

**1. Hardcoded weight (FIXED)**

**Location:** `DofData::create_dofwise_t_matrices_master()` in `src/fem/kernel/cl_FEM_DofMgr_DofData.cpp`

An earlier revision set the trivial-dependency weight to `1.0`, which ignored the
orientation sign of reversed edge-to-edge couplings. The code now passes
`tDof->mesh_basis()->weight( 0 )`. Kept as history.

---

**2. Vertex vs Source Container Usage (CLARIFIED - Not a Bug!)**

**Previous concern:** Using `vertex()` with `reinterpret_cast` seemed like type-unsafe code

**Actual purpose:** **Intentional design for visualization!**

**Why vertex containers are used:**
- Nodes are **duplicated and uncoupled** at interfaces for Paraview visualization
- `edge->vertex(i)` returns **reference nodes** (for visualization/output)
- `edge->node(i)` returns **geometric nodes** (physical mesh connectivity)
- These may be different nodes at interfaces!

**Why reinterpret_cast is needed:**
- `graph::Vertex` is base class (for graph algorithms)
- Storing mesh entities (`Node*`, `Edge*`) requires cast
- This is **by design**, not a bug

**Recommendation:** **Keep as-is** - this enables proper field visualization at discontinuous interfaces

**See also:** Container hierarchy section above for detailed explanation

---

### Debugging Techniques

**1. Enable debug output for hanging edges**

**Add** a temporary print in `DofData::create_dofwise_t_matrices_master()` when a hanging edge is processed:
```cpp
std::cout << "#edge " << tEdge->node(0)->id() << "->" << tEdge->node(1)->id()
          << " " << tEdge->number_of_sources() << " " << tEdge->number_of_vertices()
          << std::endl;
```

**Interpretation:**
- `number_of_sources()` = how many entities this edge hangs on
- `number_of_vertices()` = size of vertex container (should match sources)

---

**2. Print T-matrix for inspection**

```cpp
// After T-matrix creation
if ( tDof->is_hanging() )
{
    std::cout << "Hanging DOF " << tDof->id() << " depends on:" << std::endl;
    for ( uint k = 0; k < tDof->number_of_sources(); ++k )
    {
        std::cout << "  Source DOF " << tDof->source( k )->id()
                  << " with weight " << tDof->weight( k ) << std::endl;
    }
}
```

---

**3. Verify sign correctness**

**For edge-to-edge coupling, check orientation:**
```cpp
mesh::Edge * tShellEdge = ...;
mesh::Edge * tVolumeEdge = ...;

std::cout << "Shell edge: " << tShellEdge->node(0)->id() << " → "
          << tShellEdge->node(1)->id() << std::endl;
std::cout << "Volume edge: " << tVolumeEdge->node(0)->id() << " → "
          << tVolumeEdge->node(1)->id() << std::endl;

real expected_sign = (tShellEdge->node(0)->id() == tVolumeEdge->node(0)->id()) ? 1.0 : -1.0;
real actual_sign = tShellEdge->weight(0);

BELFEM_ERROR( std::abs(expected_sign - actual_sign) < 1e-12,
             "Sign mismatch! Expected %g, got %g", expected_sign, actual_sign );
```

---

**4. Test with simple geometry**

**Create minimal test case:**
- Single thin-shell element (QUAD4)
- Single volume conductor element (HEX8)
- Known current input → verify circulation matches

**Expected result:**
- If orientations match: `∫_shell H·dl = ∫_volume H·dl`
- If orientations opposite: `∫_shell H·dl = -∫_volume H·dl`

---

### MPI Considerations

**Hanging DOFs must be consistent across ranks:**

1. **DOF numbering** - Same global ID on all procs
2. **Source identification** - Same source DOFs on all procs
3. **Weight values** - Identical on all procs (deterministic)

**Why this matters:**
- Assembly tables use global DOF indices
- Inconsistent T-matrices → inconsistent assembly → wrong solution

**How to verify:**
```cpp
// On each rank, print hanging DOF info
if ( comm_rank() == 0 || mHangingDOFs.size() > 0 )
{
    std::cout << "Rank " << comm_rank() << ": "
              << mHangingDOFs.size() << " hanging DOFs" << std::endl;

    for ( Dof * tDof : mHangingDOFs )
    {
        std::cout << "  DOF " << tDof->id() << " → ";
        for ( uint k = 0; k < tDof->number_of_sources(); ++k )
        {
            std::cout << tDof->source(k)->id() << "(" << tDof->weight(k) << ") ";
        }
        std::cout << std::endl;
    }
}
```

---

## Literature References

### Primary BELFEM Papers

**Messe et al. 2023** - BELFEM core reference
- Section 2.2: H-φ formulation motivation
- Section 2.5: Thin-shell interface coupling
- Lines 509-517: Static condensation vs Lagrange multipliers
- DOI: (check literature/papers/fem/messe2023.txt)

**Arsenault et al. 2023** - Magnetodynamic H-φ formulation
- Equation 5-7: Interface conditions
- Equation 8: Edge-to-node coupling formula
- DOI: (check literature/papers/fem/arsenault2023.txt)

**Alves et al. 2022b** - Thin-shell theory
- Section 3.2: H-φ interface implementation
- Equation 10: Cohomology basis (cuts)
- DOI: (check literature/papers/fem/alves2022b.txt)

### General FEM Theory

**Bathe 2016** - Finite Element Procedures
- §3.2.3: Constraint equations and Lagrange multipliers
- §8.2: Static condensation
- Chapter 8: Solvers and numerical stability

**Monk 2003** - Finite Element Methods for Maxwell's Equations
- Chapter 5: Nédélec edge elements
- Section 5.3: Edge element orientation and circulation
- DOI: 10.1093/acprof:oso/9780198508885.001.0001

**Hughes 2000** - The Finite Element Method
- Chapter 3: Weak forms and mixed formulations
- Chapter 4: Element technology

**Brenner & Scott 2008** - Mathematical Theory of FEM
- Chapter 5: Mixed methods
- Chapter 11: Inf-sup stability

**Arnold 2018** - Finite Element Exterior Calculus
- Chapter 4: H(curl) spaces and edge elements
- Section 4.3: Whitney forms

**Boffi et al. 2013** - Mixed Finite Element Methods
- Chapter 2: Saddle-point problems
- Chapter 4: Inf-sup condition
- Section 6.3: Static condensation vs Lagrange multipliers

### Implementation References

**See also:**
- `dof_manager_usage_guide.md` - DofManager comprehensive guide

---

## Quick Reference

### Key Files

| File | Lines | Purpose |
|------|-------|---------|
| `cl_MaxwellFactory.cpp` | 1038-1064 | Top-level hanging orchestration |
| `cl_MaxwellFactory.cpp` | 1251-1317 | Edge→node (bottom) |
| `cl_MaxwellFactory.cpp` | 1320-1387 | Edge→node (top) |
| `cl_MaxwellFactory.cpp` | 1390-1444 | Edge→edge (bottom) |
| `cl_MaxwellFactory.cpp` | 1447-1501 | Edge→edge (top) |
| `cl_FEM_DofMgr_DofData.cpp` | 3486-3784 | DOF-level T-matrix construction |
| `cl_FEM_DofMgr_DofData.cpp` | 3505-3520 | Trivial dependency (weight bug fixed — uses the basis weight) |
| `cl_FEM_DofMgr_DofData.cpp` | 3572-3686 | Edge→node T-matrix |
| `cl_FEM_DofMgr_DofData.cpp` | 3687-3741 | Edge→edge T-matrix |
| `cl_FEM_DofManager.cpp` | 430-629 | Assembly with T-matrices |
| `cl_FEM_Tmatrix.cpp` | 100-153 | T-matrix projection operations |

### Container Access Patterns

| Container | Access | Type | Purpose |
|-----------|--------|------|---------|
| Source (mesh) | `edge->source(i)` | `mesh::Basis*` | Source entities (geometric coupling) |
| Weight (mesh) | `edge->weight(i)` | `real` | Entity-level weight (e.g., edge orientation) |
| Source (DOF) | `dof->source(i)` | `fem::Dof*` | Source DOFs (individual DOF dependencies) |
| Weight (DOF) | `dof->weight(i)` | `real` | Per-DOF T-matrix coefficient |
| Vertex (graph) | `edge->vertex(i)` | `graph::Vertex*` | Reference nodes/edges for visualization |
| Node (mesh) | `edge->node(i)` | `mesh::Node*` | Geometric nodes (physical connectivity) |

### T-Matrix Dimensions

| Coupling Type | Rows (hanging) | Columns (master) | Typical NNZ |
|---------------|----------------|------------------|-------------|
| Edge→2 nodes (LINE2) | 1 | 2 | 2 |
| Edge→3 nodes (LINE3) | 2 | 3 | 6 |
| Edge→1 edge (LINE2) | 1 | 1 | 1 |
| Edge→2 edges (LINE3) | 2 | 2 | 2 |
| Face→nodes | 1 | 3-4 | 3-4 |

---

**Last updated:** 2026-01-30
