# DofManager and FEM Kernel Usage Guide {#fem_kernel_dof_manager_usage_guide}

**Module:** src/fem/kernel
**Version:** 1.0
**Date:** 2026-01-20
**Purpose:** Comprehensive user guide for BELFEM's FEM Kernel and DOF management system

**Revision History:**
- **v1.0** (2026-01-20) - Initial comprehensive guide
  - Documented then-current bug in `number_of_hanging_dofs()`
  - Added Safe Defaults Quick Start
  - Added Glossary
- **v1.1** (2026-08-14) - The `number_of_hanging_dofs()` bug has been fixed
  in the source; all "CONFIRMED BUG" sections rewritten as historical notes
  - Added Ownership and Lifetime section
  - Comprehensive JEDI system documentation
  - MPI distribution patterns
  - Calculator usage patterns
- **v1.2** (2026-08-28) - Sections 19.7 and 20.3-20.7 rewritten against the
  current thin-shell code
  - The per-material `h_ts_*` kernel family documented here no longer exists;
    kernel selection is by solver algorithm, and the thin-shell specifics live
    in the `calculator::MaxwellData` constructor dispatch
  - The normal-field listing showed an unprojected master/slave average; the
    live `compute_hn` projects onto the facet normal
  - Anisotropy section rewritten around `bn_angle` (unfolded to `[0, pi]`),
    with the Jc extrema corrected and `bn_angle` separated from `bj_angle`

---

## Table of Contents

1. [Common Pitfalls (Read This First!)](#1-common-pitfalls-read-this-first)
2. [Safe Defaults Quick Start](#2-safe-defaults-quick-start)
3. [Glossary](#3-glossary)
4. [Overview and Architecture](#4-overview-and-architecture)
5. [The Four Pillars of FEM Kernel](#5-the-four-pillars-of-fem-kernel)
6. [DOF Allocation and Distribution](#6-dof-allocation-and-distribution)
7. [The JEDI Block Matrix System](#7-the-jedi-block-matrix-system)
8. [MPI Parallel Distribution](#8-mpi-parallel-distribution)
9. [Calculator and Element-Level Assembly](#9-calculator-and-element-level-assembly)
10. [Boundary Conditions](#10-boundary-conditions)
11. [Hanging Nodes and T-Matrices](#11-hanging-nodes-and-t-matrices)
12. [Postprocessing (L2 Projection)](#12-postprocessing-l2-projection)
13. [Ownership and Lifetime](#13-ownership-and-lifetime)
14. [Performance Considerations](#14-performance-considerations)
15. [Thread Safety and MPI](#15-thread-safety-and-mpi)
16. [Known Issues and Bugs](#16-known-issues-and-bugs)
17. [Development Notes](#17-development-notes)
18. [Literature References](#18-literature-references)
19. [H-Phi Mixed Formulation (Electromagnetics)](#19-h-phi-mixed-formulation-electromagnetics)
20. [Thin Shell Formulation](#20-thin-shell-formulation)
21. [Maxwell-Specific Implementation](#21-maxwell-specific-implementation)

---

## 1. Common Pitfalls (Read This First!)

⚠️ **Critical bugs and common mistakes to avoid:**

### 1.1 FIXED: `number_of_hanging_dofs()` copy-paste accessor (historical)

An earlier revision of `cl_FEM_DofMgr_DofData.hpp` returned
`mNumberOfFixedDofs` from `number_of_hanging_dofs()`. **The accessor is
fixed** — it returns `mNumberOfHangingDofs` — so both
`number_of_hanging_dofs()` and `my_number_of_hanging_dofs()` are safe to use.
The entry is kept only as a reminder that the count accessors in this header
are copy-paste-prone; when adding one, check the returned member matches the
name.

---

### 1.2 Forgetting to Call `initialize()`

**The Mistake:**
```cpp
DofManager* field = kernel.create_field(iwg);
field->sideset(1)->impose_dirichlet(300.0);

// ❌ Forgot to initialize!
field->compute_jacobian_and_rhs();  // CRASH! Matrices are nullptr
```

**Correct Pattern:**
```cpp
DofManager* field = kernel.create_field(iwg);
field->sideset(1)->impose_dirichlet(300.0);

field->initialize();  // ✅ Allocates DOFs and matrices

field->compute_jacobian_and_rhs();  // Works correctly
```

**Why This Happens:**
- `initialize()` orchestrates DOF creation, graph extraction, matrix allocation
- Without it, all internal pointers are `nullptr`
- Assembly functions will dereference null pointers → segfault

**Detection:**
```cpp
BELFEM_ERROR(field->solver() != nullptr, "Call initialize() first!");
```

---

### 1.3 Accumulating Assembly Without Reset

**The Mistake:**
```cpp
// First assembly
field->compute_jacobian_and_rhs(false);  // aReset = false

// Second assembly (same timestep)
field->compute_jacobian_and_rhs(false);  // VALUES ACCUMULATE!
```

**Result:**
- Matrix and RHS values are **added** to previous values
- System becomes 2× (or N×) the correct magnitude
- Solver produces garbage results

**Correct Pattern:**
```cpp
// Default behavior: reset before assembly
field->compute_jacobian_and_rhs();  // aReset = true (default)

// OR explicitly reset:
field->compute_jacobian_and_rhs(true);
```

**When to Use `aReset=false`:**
- Multi-block assembly where you want to accumulate contributions
- Custom assembly loops where you manually reset at the start

---

### 1.4 Accessing JEDI Matrices Before Initialization

**The Mistake:**
```cpp
DofManager* field = kernel.create_field(iwg);

SpMatrix* J = field->system_matrix();  // nullptr!
J->set_value(0, 0, 1.0);  // CRASH!
```

**Correct Pattern:**
```cpp
DofManager* field = kernel.create_field(iwg);
field->initialize();  // ✅ Allocates matrices

SpMatrix* J = field->system_matrix();  // Valid pointer
BELFEM_ASSERT(J != nullptr, "Jacobian should be allocated");
```

**All JEDI Matrices:**
```cpp
SpMatrix* A = field->system_matrix(); // A: Free → Free
SpMatrix* D = field->dirichlet();     // D: Free → Fixed
SpMatrix* E = field->enforcement();   // E: Fixed → Free
SpMatrix* I = field->imposition();    // I: Fixed → Fixed
```

All are `nullptr` until `initialize()` is called.

---

### 1.5 Modifying DOF Values During Assembly

**The Mistake:**
```cpp
// Inside element loop during assembly
for (Element* elem : block->elements()) {
    iwg->compute_jacobian(elem, K);

    // ❌ WRONG: Modifying DOF values during assembly
    for (Dof* dof : elem->dofs()) {
        dof->set_value(some_new_value);  // Race condition in parallel!
    }

    solver_data->assemble_jacobian(elem, K);
}
```

**Why This Is Wrong:**
- Assembly reads DOF values (for nonlinear terms, material properties)
- Modifying values during assembly creates race conditions
- Results depend on element traversal order (non-deterministic)

**Correct Pattern:**
```cpp
// Modify DOF values BEFORE assembly or AFTER solve
field->init_dof_values();  // Set initial values

field->compute_jacobian_and_rhs();  // Assembly (read-only access to DOF values)
field->solve();  // Solver updates DOF values

// Now safe to read updated values
for (Dof* dof : field->dofs()) {
    real new_value = dof->value();  // Updated by solver
}
```

---

### 1.6 Not Synchronizing Fields in MPI

**The Mistake:**
```cpp
// Parallel run (MPI)
field->solve();  // Master solves, workers have stale data

// ❌ Workers access stale field values!
mesh.field("T")->data()[node_index];  // Wrong on workers!
```

**Correct Pattern:**
```cpp
// After solve, synchronize fields
field->solve();
field->synchronize_fields({"T", "q"});  // Collect then distribute

// Now all procs have updated values
mesh.field("T")->data()[node_index];  // Correct on all procs
```

**Field Synchronization Methods:**
```cpp
field->collect_fields({"T"});      // Workers → Master
field->distribute_fields({"T"});   // Master → Workers
field->synchronize_fields({"T"});  // Round-trip (collect + distribute)
```

---

### 1.7 Deleting IWG Before DofManager

**The Mistake:**
```cpp
IWG* iwg = factory.create_iwg(...);
DofManager* field = kernel.create_field(iwg);

delete iwg;  // ❌ WRONG: DofManager still holds pointer!

field->compute_jacobian_and_rhs();  // Dangling pointer access!
```

**Correct Pattern 1: Let Kernel Manage Lifetime**
```cpp
// Kernel owns both IWG and DofManager
Kernel kernel(&params);
IWG* iwg = kernel.create_equation(...);
DofManager* field = kernel.create_field(iwg);

// ... use field ...

// Kernel destructor cleans up both
```

**Correct Pattern 2: Manual Management**
```cpp
IWG* iwg = factory.create_iwg(...);
DofManager* field = new DofManager(&kernel, 0);
field->set_equation(iwg);

// ... use field ...

delete field;  // Delete DofManager first
delete iwg;    // Then delete IWG
```

**See:** [Section 13: Ownership and Lifetime](#13-ownership-and-lifetime)

---

### 1.8 Mixing DOF ID and DOF Index

**The Mistake:**
```cpp
// DOF ID (global, permanent)
id_t dof_id = dof_data->node_dof_id(node->id(), dof_type);

// ❌ Using ID as index into array!
Vector<real> values(num_dofs);
values(dof_id) = 1.0;  // WRONG! ID can be > num_dofs
```

**Correct Pattern:**
```cpp
// DOF ID (global, permanent, before reordering)
id_t dof_id = dof_data->node_dof_id(node->id(), dof_type);

// Get DOF object from ID
Dof* dof = dof_data->dof(dof_id);

// Use DOF index (reordered position)
index_t dof_index = dof->index();

// Now safe to use as array index
Vector<real> values(num_dofs);
values(dof_index) = 1.0;  // ✅ CORRECT
```

**Key Difference:**
- **ID:** Global unique identifier, permanent, can be sparse
- **Index:** Reordered position in `[0, num_dofs)`, dense, changes after reordering

---

### 1.9 Creating DofManager Inside Element Loops

**The Mistake:**
```cpp
// ❌ EXTREMELY INEFFICIENT:
for (Element* elem : block->elements()) {
    DofManager* field = kernel.create_field(iwg);  // Allocates DOFs, matrices!
    field->initialize();
    // ...
}
```

**Why This Is Wrong:**
- DOF creation is O(num_nodes + num_edges + num_faces + num_elements)
- Matrix allocation is O(num_dofs × avg_connectivity)
- Creating inside loop → O(num_elements × num_dofs²) complexity!

**Correct Pattern:**
```cpp
// ✅ Create once, use many times:
DofManager* field = kernel.create_field(iwg);
field->initialize();  // Once

for (Element* elem : block->elements()) {
    // Use field for assembly
    iwg->compute_jacobian(elem, K);
    solver_data->assemble_jacobian(elem, K);
}
```

---

## 2. Safe Defaults Quick Start

Minimal working example with safe defaults:

```cpp
#include "cl_FEM_Kernel.hpp"
#include "cl_IwgFactory.hpp"

int main()
{
    // 1. Load the mesh first -- KernelParameters is built around it.
    //    There is no default KernelParameters and no Kernel::import_mesh().
    Mesh* mesh = new Mesh("mesh.exo");

    KernelParameters params(mesh);
    Kernel kernel(&params);

    // 2. Create equation (heat conduction example)
    IWG* iwg = kernel.create_equation(
        IwgType::TransientHeatConduction,  // ✅ Well-tested physics
        ModelDimensionality::ThreeD
    );

    // Configure time-stepping (SAFE choices)
    iwg->select_blocks({1, 2, 3});
    iwg->set_timestepping_method(EulerMethod::BackwardDifference2);  // ✅ BDF2: stable, 2nd order
    iwg->set_algorithm(SolverAlgorithm::NewtonRaphson);
    iwg->set_omega(0.9);  // Line search relaxation

    // 3. Create field (DofManager)
    DofManager* field = kernel.create_field(iwg);

    // 4. Set boundary conditions (block 1 = Dirichlet, block 2 = Neumann)
    field->sideset(1)->impose_dirichlet(300.0);  // T = 300K
    field->sideset(2)->impose_neumann(1000.0);   // q = 1000 W/m²

    // 5. Set solver (SAFE choices)
    SolverParameters solver_params( SolverType::MUMPS );  // ✅ Robust, parallel-ready; the type is fixed at construction
    solver_params.set_reordering_method( ReorderingMethod::METIS );  // ✅ Nested dissection
    field->set_solver( solver_params );

    // 6. Initialize (CRITICAL STEP - allocates DOFs and matrices)
    field->initialize();

    // 7. Time-stepping loop
    real dt = 0.01;  // Time step
    real t = 0.0;
    real t_end = 1.0;

    while (t < t_end) {
        // Assembly (aReset=true by default - safe)
        field->compute_jacobian_and_rhs(true);

        // Solve
        field->solve();

        // Postprocess (L2 projection of secondary fields)
        field->postprocess();

        // Synchronize fields for MPI (safe even in serial)
        field->synchronize_fields({"T"});

        // Export results
        if (comm_rank() == 0) {
            mesh->save("output_" + std::to_string(t) + ".exo");
        }

        t += dt;
    }

    // 8. Cleanup handled by Kernel destructor
    return 0;
}
```

**Why These Defaults Are Safe:**

| Choice | Reason |
|--------|--------|
| `BDF2` | A-stable, 2nd order accuracy, well-tested |
| `MUMPS` | Robust, parallel-ready, handles ill-conditioned systems |
| `METIS` | Optimal fill-minimization for sparse factorization |
| `NewtonRaphson` | Full Jacobian, quadratic convergence |
| `aReset=true` | Prevents accumulation bugs |
| `synchronize_fields()` | Safe in both serial and parallel |
| `postprocess()` | Ensures secondary fields are up-to-date |

---

## 3. Glossary

| Term | Full Name / Definition |
|------|------------------------|
| **DOF** | Degree Of Freedom - unknowns in the FEM system (e.g., temperature at each node) |
| **DofManager** | Central orchestrator for DOF allocation, assembly, and solving |
| **DofData** | Low-level DOF container and ID calculation |
| **JEDI** | Jacobian/Enforcement/Dirichlet/Imposition - block matrix system for handling BCs |
| **Free DOF** | Unknown to solve for (not constrained by Dirichlet BC) |
| **Fixed DOF** | Constrained by Dirichlet BC (imposed value) |
| **Hanging DOF** | Constrained by mesh refinement (weighted sum of parent DOFs) |
| **T-matrix** | Transformation matrix for hanging node constraints |
| **Aura** | Halo of elements/nodes on neighboring MPI procs (for ghost data) |
| **Owned** | Entities owned by this MPI proc (master copy) |
| **BC** | Boundary Condition (Dirichlet, Neumann, Robin) |
| **RCM** | Reverse Cuthill-McKee - bandwidth reduction algorithm |
| **METIS** | Graph partitioning library for nested dissection and MPI decomposition |
| **Calculator** | Element-level computation engine (shape functions, gradients, Jacobian) |
| **N** | Shape function matrix (interpolation) |
| **B** | Gradient operator matrix (B = inv(J) * dN/dξ) |
| **J (Jacobian)** | Free → Free coupling matrix (main system to solve) |
| **D (Dirichlet)** | Free → Fixed coupling matrix (Dirichlet BCs) |
| **E (Enforcement)** | Fixed → Free coupling matrix (reaction forces) |
| **I (Imposition)** | Fixed → Fixed coupling matrix (self-coupling of imposed values) |
| **L2 Projection** | Least-squares projection of discontinuous fields to continuous mesh fields |
| **Assembly** | Accumulation of element contributions into global system |
| **Facet** | Mesh face used for boundary conditions (triangle, quadrilateral) |
| **Lambda DOF** | Lagrange multiplier DOF (for interface constraints) |

---

## 4. Overview and Architecture

### 4.1 What Is the FEM Kernel?

The **FEM Kernel** module is the **core orchestration layer** that connects:
- Mesh data structures (`Mesh`, `Node`, `Element`)
- Physics equations (IWG - Integral Weak Form)
- DOF management (allocation, numbering, constraints)
- Sparse linear algebra (matrix assembly, solvers)
- Parallel computing (MPI distribution, field synchronization)

**Central Class:** `DofManager` - orchestrates all DOF operations.

**Top-Level Class:** `Kernel` - owns mesh, materials, boundary conditions, and multiple DofManagers.

---

### 4.2 Class Hierarchy

```
Kernel (Top-level orchestrator)
├── Mesh* (borrowed - externally owned or kernel-owned)
├── Materials (owned)
├── BoundaryConditions (owned)
├── Cell<IWG*> (equations - owned)
└── Cell<DofManager*> (fields - owned)

DofManager (Central orchestrator)
├── Inherits from: DofManagerBase
├── Composition pattern: owns 8 specialized data managers
│   ├── Parameters*
│   ├── DofData*
│   ├── BlockData*
│   ├── SideSetData*
│   ├── BearingData*
│   ├── FieldData*
│   ├── SolverData*       (JEDI block matrix system)
│   └── EigenValues*
└── Cell<Postprocessor*> (L2 projection - owned)
```

---

### 4.3 Composition Pattern

**DofManager uses composition to delegate responsibilities:**

| Manager | File | Responsibility |
|---------|------|----------------|
| **Parameters** | `cl_FEM_DofMgr_Parameters` | Integration orders, scheme settings |
| **DofData** | `cl_FEM_DofMgr_DofData` | DOF creation, numbering, hanging nodes |
| **BlockData** | `cl_FEM_DofMgr_BlockData` | Block management, thin-shell facet linking |
| **SideSetData** | `cl_FEM_DofMgr_SideSetData` | Sideset management, BCs, wetted nodes |
| **BearingData** | `cl_FEM_DofMgr_BearingData` | Bearing element management |
| **FieldData** | `cl_FEM_DofMgr_FieldData` | Field collect/distribute, linear projection |
| **SolverData** | `cl_FEM_DofMgr_SolverData` | Matrix assembly, solver interface, JEDI system |
| **EigenValues** | `cl_FEM_DofMgr_EigenValues` | Eigenvalue problem setup |

**Each manager is:**
- Created in `DofManager` constructor
- Initialized in `DofManager::initialize()`
- Deleted in `DofManager` destructor

**User interacts primarily with `DofManager`**, which forwards calls to appropriate managers.

---

### 4.4 Key Abstractions

**Group → Block / SideSet:**
- `Group` is abstract base class
- `Block` represents volume elements (material domains)
- `SideSet` represents surface elements (boundary conditions)
- Both own a `Calculator` for element-level computations

**Element (FEM wrapper):**
- Wraps `mesh::Element` with FEM-specific data
- Links DOFs to mesh entities (nodes, edges, faces, elements, facets)
- Manages T-matrices for hanging nodes
- Stores edge directions for Nédélec elements

**Dof (Degree of Freedom):**
- Inherits from `graph::Vertex` (for graph algorithms)
- Links to `mesh::Basis` (Node, Edge, Face, Element, Facet)
- Stores value, fixed flag, Dirichlet BC value
- Manages hanging node sources and weights

---

## 5. The Four Pillars of FEM Kernel

Understanding DofManager requires understanding **four core concepts**:

### Pillar 1: DOF Allocation

**Process:**
1. IWG specifies DOF requirements (nodes, edges, faces, cells, lambdas)
2. DofData counts DOFs per entity type
3. DofData computes ID offsets for contiguous ranges
4. DofData allocates DOF objects
5. DofData creates fast lookup map (ID → Dof*)

**Result:** All DOFs have unique global IDs

---

### Pillar 2: Graph Extraction and Reordering

**Process:**
1. SolverData extracts DOF-DOF connectivity graph from mesh
2. DofData splits DOFs into free and fixed based on BCs
3. DofData reorders free DOFs (symrcm; solver-internal nested dissection happens in the solver wrapper)
4. DofData reorders fixed DOFs separately
5. Each DOF gets new index in [0, num_dofs)

**Result:** DOFs are renumbered for optimal sparse solver performance

---

### Pillar 3: JEDI Matrix System

**Process:**
1. SolverData allocates 4 sparse matrices (J, D, E, I)
2. SolverData creates assembly tables for MPI distribution
3. Element loop: IWG computes local matrices
4. SolverData assembles into global JEDI matrices
5. Solver solves `J*x = f - D*y` for free DOFs
6. SolverData computes reaction forces `g = E*x + I*y`

**Result:** Solution vector + reaction forces at constraints

---

### Pillar 4: MPI Field Synchronization

**Process:**
1. Each proc has owned + aura DOFs
2. Assembly: each proc contributes to global system
3. Solve: master solves (or parallel solve with subdomain procs)
4. Distribute: solution scattered to all procs
5. Synchronize: round-trip to ensure consistency

**Result:** All procs have consistent field values

---

These four pillars execute sequentially in `DofManager::initialize()` (Pillar 1-2) and `DofManager::solve()` (Pillar 3-4).

---

## 6. DOF Allocation and Distribution

### 6.1 DOF Creation Workflow

The `DofManager::initialize()` method triggers DOF creation:

```cpp
DofManager::initialize()
├── init_dofs()
│   └── DofData::create_dofs(IWG*)
│       ├── 1. Query IWG for DOF requirements
│       │   ├── IWG::number_of_dofs_per_node()
│       │   ├── IWG::number_of_dofs_per_edge()
│       │   ├── IWG::number_of_dofs_per_face()
│       │   ├── IWG::number_of_dofs_per_cell()
│       │   └── IWG::number_of_lambda_dofs()
│       │
│       ├── 2. Count DOFs per entity type
│       │   ├── count_node_dofs()
│       │   ├── count_edge_dofs()
│       │   ├── count_face_dofs()
│       │   ├── count_cell_dofs()
│       │   └── count_lambda_dofs()
│       │
│       ├── 3. Compute ID offsets
│       │   └── compute_dof_offsets(IWG*)
│       │       ├── mEdgeDofOffset = max_node_id * num_dof_types
│       │       ├── mFaceDofOffset = mEdgeDofOffset + max_edge_id * num_dof_types
│       │       ├── mCellDofOffset = mFaceDofOffset + max_face_id * num_dof_types
│       │       └── mLambdaDofOffset = mCellDofOffset + max_cell_id * num_dof_types
│       │
│       ├── 4. Allocate DOF objects
│       │   ├── Create Dof objects (one per DOF)
│       │   ├── Link to mesh::Basis (Node, Edge, Face, Element, Facet)
│       │   └── Set DOF ID (unique global identifier)
│       │
│       └── 5. Create DOF map
│           └── create_dof_map()
│               └── mDofMap( dof_id ) = dof_pointer
│
├── init_work()
│   └── SolverData::extract_graph_from_mesh()
│       ├── compute_element_dof_connectivity()
│       ├── compute_dof_element_connectivity()
│       └── compute_dof_dof_connectivity()
│
└── init_matrices()
    └── SolverData::allocate_matrices()
        ├── Reorder DOFs (symrcm)
        ├── Allocate sparse matrices (J, D, E, I)
        ├── Create assembly tables
        └── Initialize solver
```

---

### 6.2 DOF ID Calculation

**Contiguous ID ranges with offsets:**

```cpp
// Node DOFs: [0, max_node_id * num_dof_types)
node_dof_id = nodeID * mNumDofTypes + dofType

// Edge DOFs: [edge_offset, edge_offset + max_edge_id * num_dof_types)
edge_dof_id = mEdgeDofOffset + edgeID * mNumDofTypes + dofType

// Face DOFs: [face_offset, face_offset + max_face_id * num_dof_types)
face_dof_id = mFaceDofOffset + faceID * mNumDofTypes + dofType

// Cell DOFs: [cell_offset, cell_offset + max_cell_id * num_dof_types)
cell_dof_id = mCellDofOffset + cellID * mNumDofTypes + dofType

// Lambda DOFs: [lambda_offset, lambda_offset + max_facet_id * num_dof_types)
lambda_dof_id = mLambdaDofOffset + facetID * mNumDofTypes + dofType
```

**Example:** Heat conduction (1 DOF type: temperature)
- Node 5: DOF ID = 5 * 1 + 0 = 5
- Node 10: DOF ID = 10 * 1 + 0 = 10

**Example:** Electromagnetics (3 DOF types: Ax, Ay, Az)
- Edge 100, type 0 (Ax): DOF ID = edge_offset + 100 * 3 + 0
- Edge 100, type 1 (Ay): DOF ID = edge_offset + 100 * 3 + 1
- Edge 100, type 2 (Az): DOF ID = edge_offset + 100 * 3 + 2

**Implementation:** `cl_FEM_DofMgr_DofData.hpp:428-474`

---

### 6.3 DOF Numbering vs. Reordering

**DOF ID** (initial, global, permanent):
- Calculated from mesh entity IDs
- Unique across all procs
- Sparse (gaps in numbering)
- Never changes

**DOF Index** (reordered, dense, local):
- Position in reordered array [0, num_dofs)
- Dense (no gaps)
- Changes after reordering (RCM, METIS)
- Used for sparse matrix indexing

**Workflow:**
```cpp
// 1. Initial numbering (ID-based)
Dof* dof = new Dof();
dof->set_id(node_dof_id);          // ID from mesh entity
dof->set_index(dof->id());         // Initially index = ID

// 2. Reordering (index updated)
Graph free_dofs;  // Create connectivity graph
symrcm(free_dofs);  // bandwidth reduction; nested dissection is solver-internal

for (index_t i = 0; i < free_dofs.size(); ++i) {
    Dof* dof = (Dof*) free_dofs(i);
    dof->set_index(i);  // New index in [0, num_free_dofs)
}

// 3. Matrix assembly uses index, not ID
SpMatrix* J = new SpMatrix(num_free_dofs, num_free_dofs);
J->set_value(dof_i->index(), dof_j->index(), value);  // ✅ Use index
```

---

### 6.4 Free vs. Fixed DOFs

**Split Criteria:**
```cpp
DofData::split_dof_container(free_dofs, fixed_dofs)
{
    for (Dof* dof : mDOFs) {
        if (dof->is_fixed()) {
            fixed_dofs.push(dof);  // Dirichlet BC applied
        } else {
            free_dofs.push(dof);   // Unknown to solve for
        }
    }
}
```

**Free DOFs:**
- Unknowns to solve for
- Form the main system `J*x = f - D*y`
- Numbered [0, num_free_dofs)

**Fixed DOFs:**
- Constrained by Dirichlet BCs
- Known values (imposed)
- Numbered [0, num_fixed_dofs)
- Used to compute reaction forces `g = E*x + I*y`

**Hanging DOFs:**
- NOT free or fixed
- Eliminated via T-matrix transformation
- NOT part of global system

---

### 6.5 Hanging DOF Detection and Elimination

**What Are Hanging Nodes?**

Hanging nodes occur at mesh refinement boundaries:
```
Coarse element:       Refined element:
a-----------b         a-----h-----b
                           ^
                      Hanging node h = 0.5*a + 0.5*b
```

**Data Structure:**
```cpp
class Dof {
    uint mNumberOfSources;      // Number of parent DOFs
    Dof** mSources;             // Parent DOF pointers
    real* mCoefficients;        // Interpolation weights
};

// Example: h = 0.5*a + 0.5*b
h->mNumberOfSources = 2;
h->mSources = [a_dof, b_dof];
h->mCoefficients = [0.5, 0.5];
```

**Workflow:**
```cpp
DofData::collect_hanging_dofs()
├── 1. Traverse all nodes/edges/faces
├── 2. Check if entity is hanging (mesh::Node::is_hanging())
├── 3. Get parent entities and weights
└── 4. Store in mHangingDOFs container

DofData::remove_hanging_dofs_from_container()
├── Remove hanging DOFs from mDOFs
└── Keep in separate mHangingDOFs container

Element::relink_dofs()
├── Replace hanging DOF pointers with T-matrix
└── T-matrix transforms element DOFs to free DOFs
```

**Note:** the historical `number_of_hanging_dofs()` accessor bug (§1.1, fixed)
never affected the hanging-node elimination logic — that uses T-matrices.

---

## 7. The JEDI Block Matrix System

### 7.1 What Is JEDI?

**JEDI = Jacobian / Enforcement / Dirichlet / Imposition**

A block matrix system that **separates free and fixed DOFs** to enable:
1. **Efficient solving** (J matrix is smaller and better conditioned)
2. **Reaction force computation** (E matrix computes constraint forces)
3. **Static condensation** (D and E couple free and fixed DOFs)

**The JEDI system literally computes "the Force"** (reaction forces at constraints)!

---

### 7.2 Block Matrix Structure

**Full system:**
```
[ J  D ] [ x ]   [ f ]
[ E  I ] [ y ] = [ g ]
```

Where:
- `x` = free DOFs (unknowns to solve for)
- `y` = fixed DOFs (imposed Dirichlet values)
- `f` = RHS for free DOFs
- `g` = RHS for fixed DOFs (reaction forces)

**Matrix meanings:**

| Matrix | Dimension | Meaning | Physical Interpretation |
|--------|-----------|---------|-------------------------|
| **J** | (N_free × N_free) | Free → Free coupling | Stiffness/Jacobian for unknowns |
| **D** | (N_free × N_fixed) | Fixed → Free dependency | How imposed BCs affect free DOFs |
| **E** | (N_fixed × N_free) | Free → Fixed enforcement | **Reaction forces** at constraints |
| **I** | (N_fixed × N_fixed) | Fixed → Fixed self-coupling | BC self-interaction |

---

### 7.3 Solving the JEDI System

**Step 1: Solve for free DOFs**
```
J*x = f - D*y
```

- `y` is known (imposed Dirichlet values)
- `D*y` term moves fixed DOF effects to RHS
- Solve for `x` using sparse solver (UMFPACK, MUMPS, STRUMPACK, etc.)

**Step 2: Compute reaction forces (optional)**
```
g = E*x + I*y
```

- `x` is now known (from Step 1)
- `g` contains reaction forces at Dirichlet boundaries
- In structural mechanics, this is literally "the Force" at supports!

**Implementation:**
```cpp
// Step 1: Solve for free DOFs
field->solve();  // Internally solves J*x = f - D*y

// Step 2: Access reaction forces (if needed)
solver_data->use_jedi_force(true);  // Enable reaction force computation
Vector<real>& g = solver_data->rhs_vector();  // Contains reaction forces
```

---

### 7.4 JEDI Assembly Pattern

**Element stiffness matrix K_e maps element DOFs to element DOFs.**

**Split contributions based on free/fixed status:**

```cpp
for each (dof_i, dof_j) in K_e:

    if (dof_i.is_free() && dof_j.is_free()):
        J[i_free, j_free] += K_e[i, j]      // Free → Free

    elif (dof_i.is_free() && dof_j.is_fixed()):
        D[i_free, j_fixed] += K_e[i, j]     // Fixed → Free

    elif (dof_i.is_fixed() && dof_j.is_free()):
        E[i_fixed, j_free] += K_e[i, j]     // Free → Fixed

    elif (dof_i.is_fixed() && dof_j.is_fixed()):
        I[i_fixed, j_fixed] += K_e[i, j]    // Fixed → Fixed
```

**Implementation:** `SolverData::assemble_jacobian()` in `cl_FEM_DofMgr_SolverData.cpp`

---

### 7.5 JEDI vs. Standard FEM

**Standard FEM:**
```
K*u = f

// Apply Dirichlet BCs by:
// 1. Zero out rows/columns for fixed DOFs
// 2. Set diagonal to 1, RHS to imposed value
// 3. Lose reaction force information
```

**JEDI FEM:**
```
[ J  D ] [ x ]   [ f ]
[ E  I ] [ y ] = [ g ]

// Advantages:
// 1. Preserve full system structure
// 2. Compute reaction forces explicitly
// 3. Better conditioning (J is smaller)
// 4. Support static condensation
```

**When to use JEDI:**
- Structural mechanics (need reaction forces at supports)
- Contact problems (enforce constraints, compute contact forces)
- Interface coupling (static condensation)
- Large Dirichlet BC count (J is much smaller than K)

**When standard FEM is sufficient:**
- Pure Neumann problems (no fixed DOFs)
- Don't need reaction forces
- Small Dirichlet BC count (overhead not worth it)

---

### 7.6 Accessing JEDI Matrices

```cpp
// After initialize()
field->initialize();

// Access JEDI matrices
SpMatrix* A = field->system_matrix(); // A: Free → Free
SpMatrix* D = field->dirichlet();     // D: Free → Fixed
SpMatrix* E = field->enforcement();   // E: Fixed → Free
SpMatrix* I = field->imposition();    // I: Fixed → Fixed

// Check dimensions
BELFEM_ASSERT(J->n_rows() == num_free_dofs, "J should be square");
BELFEM_ASSERT(D->n_cols() == num_fixed_dofs, "D columns = fixed DOFs");
BELFEM_ASSERT(E->n_rows() == num_fixed_dofs, "E rows = fixed DOFs");
BELFEM_ASSERT(I->n_rows() == num_fixed_dofs, "I should be square");
```

**All matrices are `nullptr` until `initialize()` is called!**

---

## 8. MPI Parallel Distribution

### 8.1 MPI Architecture

**Two-tier parallelism:**
1. **Domain decomposition** (MPI) - coarse-grained, distributed memory
2. **Sparse solver parallelism** (optional) - fine-grained, shared memory (MUMPS, STRUMPACK)

**Proc roles:**
- **Master (rank 0):** Owns full mesh, assembles global system, runs solver (or coordinates parallel solver)
- **Workers (rank 1-N):** Own submesh, assemble local contributions, send to master

---

### 8.2 Mesh Partitioning

**Partitioning strategy:**
```cpp
Kernel::distribute_mesh()
├── Master (rank 0):
│   ├── Partition mesh using METIS/SCOTCH
│   ├── Create submeshes for each proc
│   ├── Add aura (ghost) elements for inter-proc boundaries
│   └── Send submeshes to workers
│
└── Workers (rank 1-N):
    ├── Receive submesh from master
    └── Store owned + aura elements
```

**Owned vs. Aura:**
- **Owned:** Entities with `proc_owner() == comm_rank()`
- **Aura:** Halo of entities on neighboring procs (for interpolation)

**CommTables:**
- Define communication patterns (which DOFs to send/receive)
- Created during mesh distribution
- Used for field synchronization

---

### 8.3 DOF Distribution

**DOF creation:**
```cpp
// Each proc creates DOFs for owned + aura entities
DofData::create_dofs(IWG*)
├── Count DOFs for owned entities
├── Count DOFs for aura entities
├── Allocate DOF objects
└── DOF numbering is globally consistent
```

**Global consistency:**
- All procs use same DOF ID calculation formula
- DOF IDs are unique and consistent across procs
- Reordering is done globally (master coordinates)

---

### 8.4 Matrix Assembly in Parallel

**Assembly pattern:**
```cpp
// Each proc assembles local contributions
for (Block* block : blocks) {
    for (Element* elem : block->owned_elements()) {
        // Compute local matrices
        iwg->compute_jacobian(elem, K_e);
        iwg->compute_rhs(elem, f_e);

        // Assemble into local matrices
        solver_data->assemble_jacobian(elem, K_e);
        solver_data->assemble_rhs(elem, f_e);
    }
}

// Collect contributions from all procs
solver_data->collect_matrices();
solver_data->collect_rhs_vector();
```

**Assembly tables:**
```cpp
Cell<Vector<index_t>> mSystemTable;  // Maps global DOF indices to local matrix positions

// Example: Element with DOFs [5, 12, 7]
// Global indices: [5, 12, 7]
// Local matrix positions in J: [3, 8, 4] (after reordering)
// mSystemTable stores this mapping for MPI_Send/Recv
```

**Parallel assembly workflow:**
```
Worker:                       Master:
┌─────────────┐              ┌─────────────┐
│ Assemble K_e│              │ Allocate J  │
│ for owned   │              │             │
│ elements    │              │             │
└─────┬───────┘              └─────┬───────┘
      │                            │
      │  MPI_Send(K_local)         │
      ├───────────────────────────>│
      │                            │
      │                      ┌─────▼───────┐
      │                      │ Accumulate  │
      │                      │ contributions│
      │                      │ into J      │
      │                      └─────────────┘
```

**Implementation:** `SolverData::collect_matrices()` in `cl_FEM_DofMgr_SolverData.cpp`

---

### 8.5 Field Synchronization

**Three synchronization modes:**

```cpp
// 1. Collect: Workers → Master
field->collect_fields({"T", "q"});
// Master now has full field data
// Workers still have stale data

// 2. Distribute: Master → Workers
field->distribute_fields({"T", "q"});
// All procs now have updated field data
// Workers have subset for their submesh

// 3. Synchronize: Round-trip (collect + distribute)
field->synchronize_fields({"T", "q"});
// Ensures all procs have consistent, up-to-date data
```

**When to synchronize:**
- **After solve:** Solution computed on master, must distribute to workers
- **Before postprocessing:** Workers need full solution for L2 projection
- **Before output:** All procs need consistent data for VTK/Exodus export

**Implementation:**
```cpp
void DofManager::synchronize_fields(const Cell<string>& labels)
{
    // Step 1: Collect from all procs to master
    collect_fields(labels);

    // Step 2: Distribute from master to all procs
    distribute_fields(labels);
}
```

**⚠️ Common Pitfall:** Forgetting to synchronize after solve → workers have stale data!

---

### 8.6 MPI-Safe Patterns

**Safe pattern:**
```cpp
// Master-only operations
if (comm_rank() == 0) {
    mesh->save("output.exo");  // Only master writes
}

// All-proc operations
field->synchronize_fields({"T"});  // All procs participate

// Worker-only operations (rare)
if (comm_rank() > 0) {
    // Custom worker-specific logic
}
```

**Collective operations (all procs must call):**
- `field->initialize()`
- `field->solve()`
- `field->synchronize_fields()`
- `field->collect_fields()`
- `field->distribute_fields()`

**⚠️ Deadlock if only master calls collective operation!**

---

## 9. Calculator and Element-Level Assembly

### 9.1 What Is the Calculator?

The **Calculator** is the element-level computation engine. It computes:
- **N**: Shape functions (interpolation)
- **B**: Gradient operator (B = inv(J_geom) * dN/dξ)
- **J_geom**: Geometric Jacobian (dX/dξ)
- **dV**: Differential volume (for integration)
- **dS**: Differential surface (for boundary integrals)
- **Normals**: Outward normal vectors

**File:** `cl_FEM_Calculator.{hpp,cpp}` (1742 lines - largest file in kernel!)

---

### 9.2 Calculator Usage Pattern

**Typical assembly loop:**
```cpp
Block* block = field->block(1);
Calculator* calc = block->calculator();

for (Element* elem : block->elements()) {
    // Link calculator to element
    calc->link(elem);

    // Access integration points
    uint num_gauss = elem->number_of_integration_points();

    for (uint g = 0; g < num_gauss; ++g) {
        // Shape functions N (interpolation)
        const Matrix<real>& N = calc->N(g);

        // Gradient operator B
        const Matrix<real>& B = calc->B(g);

        // Differential volume
        real dV = calc->dV(g);

        // Compute element matrix K_e
        // K_e += B^T * D * B * dV  (for diffusion/elasticity)
        K_e += trans(B) * material_D * B * dV;
    }

    // Assemble into global system
    solver_data->assemble_jacobian(elem, K_e);
}
```

---

### 9.3 Function Pointer Pattern (Performance Optimization)

**Why function pointers instead of virtual functions?**

The Calculator is called **millions of times** during assembly. Virtual function calls have ~5-10% overhead. Function pointers avoid this overhead while maintaining runtime polymorphism.

**Pattern:**
```cpp
class Calculator {
    // Function pointers for different formulations
    const Matrix<real>& (Calculator::*mFunN)(uint);
    const Matrix<real>& (Calculator::*mFunB)(uint);
    real (Calculator::*mFundV)(uint);

    void allocate() {
        // Select implementation based on dimensionality, field type and element geometry
        if (dimensionality == 2D) {
            mFunN = &Calculator::N2D;
            mFunB = &Calculator::Bscalar;          // Bplanestress for plane stress
            mFundV = &Calculator::dV_tri6_tet10;   // per element type: dV_quad4ts, dV_axsymmx, dV_axsymmy, ...
        } else {
            mFunN = &Calculator::N3D;
            mFunB = &Calculator::Bscalar;          // Bvoigt for elasticity
            mFundV = &Calculator::dV_hex;          // per element type: dV_tri6_tet10, ...
        }
    }

public:
    // User calls these (no virtual dispatch!)
    const Matrix<real>& N(uint g) { return (this->*mFunN)(g); }
    const Matrix<real>& B(uint g) { return (this->*mFunB)(g); }
    real dV(uint g) { return (this->*mFundV)(g); }
};
```

**Specializations:**
- **N functions:** Scalar, Vec2D, Vec3D (vector field interpolation)
- **B functions:** Gradient, PlaneStress, Voigt (strain-displacement operators)
- **dV functions:** Cartesian, AxiSymmX, AxiSymmY (2D), ThreeD
- **Geometry:** Straight vs. curved elements

---

### 9.4 Shape Functions (N) and Gradients (B)

**Shape Function Matrix N:**

For scalar field (e.g., temperature):
```
N = [ N1  N2  N3  ...  Nn ]    (1 × num_bases)
```

For 2D vector field (e.g., displacement):
```
N = [ N1  0   N2  0   N3  0   ... ]    (2 × 2*num_bases)
    [ 0   N1  0   N2  0   N3  ... ]
```

For 3D vector field:
```
N = [ N1  0   0   N2  0   0   ... ]    (3 × 3*num_bases)
    [ 0   N1  0   0   N2  0   ... ]
    [ 0   0   N1  0   0   N2  ... ]
```

**Gradient Operator B:**

For diffusion (scalar field):
```
B = [ dN1/dx  dN2/dx  ...  dNn/dx ]    (ndim × num_bases)
    [ dN1/dy  dN2/dy  ...  dNn/dy ]
    [ dN1/dz  dN2/dz  ...  dNn/dz ]  (3D only)
```

For elasticity (2D plane stress):
```
B = [ dN1/dx    0     dN2/dx    0     ... ]    (3 × 2*num_bases)
    [   0     dN1/dy    0     dN2/dy  ... ]
    [ dN1/dy  dN1/dx  dN2/dy  dN2/dx  ... ]
```

For elasticity (3D Voigt notation):
```
B = [ dN1/dx    0       0     ... ]    (6 × 3*num_bases)
    [   0     dN1/dy    0     ... ]
    [   0       0     dN1/dz  ... ]
    [ dN1/dy  dN1/dx    0     ... ]
    [   0     dN1/dz  dN1/dy  ... ]
    [ dN1/dz    0     dN1/dx  ... ]
```

**Computation:**
```cpp
// Geometric Jacobian: J = dN/dξ * X
// X = nodal coordinates
// ξ = parametric coordinates

Matrix<real> J_geom = dN_dxi * X_coords;

// Inverse Jacobian
Matrix<real> inv_J = inv(J_geom);

// Gradient operator: B = inv(J) * dN/dξ
Matrix<real> B = inv_J * dN_dxi;
```

---

### 9.5 Integration Point Data

**Differential volume (dV):**
```cpp
// Cartesian 2D/3D
dV = det(J_geom) * weight

// Axisymmetric about x-axis
dV = det(J_geom) * weight * 2*pi * y_coord

// Axisymmetric about y-axis
dV = det(J_geom) * weight * 2*pi * x_coord
```

**Differential surface (dS):**
```cpp
// 2D boundary (line element)
dS = sqrt((dx/dξ)^2 + (dy/dξ)^2) * weight

// 3D boundary (surface element)
dS = ||∂X/∂ξ × ∂X/∂η|| * weight
```

**Outward normal:**
```cpp
// 2D
normal = normalize([dy/dξ, -dx/dξ])

// 3D
normal = normalize(∂X/∂ξ × ∂X/∂η)
```

---

### 9.6 Nédélec Edge Elements (Electromagnetics)

For electromagnetics (H-φ formulation), edge elements use **curl-conforming** basis functions:

**Edge functions (E):**
```cpp
const Matrix<real>& E = calc->E(g);  // Edge basis functions

// E matrix maps edge DOFs to vector field in physical space
```

**Curl operator (C):**
```cpp
const Matrix<real>& C = calc->C(g);  // Curl of edge functions

// C matrix: curl of edge basis functions
// For H-formulation: ∇ × H appears in weak form
```

**Edge directions:**
```cpp
Bitset< 12 > tDirs;
elem->edge_directions( tDirs );

// Edge direction affects sign of edge basis function
// Ensures consistency across elements sharing an edge
```

**Implementation:** `cl_FEM_Calculator.cpp` (specialized functions for Nédélec elements)

---

## 10. Boundary Conditions

### 10.1 BC Types

BELFEM supports three main BC types:

| Type | Name | Equation | Example |
|------|------|----------|---------|
| **Dirichlet** | Essential BC | u = g | T = 300K (fixed temperature) |
| **Neumann** | Natural BC | ∂u/∂n = q | q = 1000 W/m² (heat flux) |
| **Robin** | Mixed BC | α*(u - u_∞) = q | h*(T - T_∞) = q (convection) |

**Implementation:**
```cpp
// Dirichlet: u = g
sideset->impose_dirichlet(300.0);  // T = 300K

// Neumann: q = constant
sideset->impose_neumann(1000.0);   // q = 1000 W/m²

// Robin (Alpha): α*(u - u_∞) = q
sideset->impose_alpha(h, T_inf);   // Convection BC
```

---

### 10.2 Dirichlet BCs (Essential)

**Application:**
```cpp
SideSet* sideset = field->sideset(1);
sideset->impose_dirichlet(300.0);  // All DOFs on sideset = 300
```

**Effect on DOF:**
```cpp
Dof* dof = /* node on sideset 1 */;
dof->is_fixed();      // true
dof->value();         // 300.0
dof->index();         // Index in fixed DOF array
```

**Effect on JEDI system:**
- DOF moved from free to fixed container
- Contributes to D, E, I matrices (not J)
- RHS modified: `f -= D*y` where `y[i] = 300.0`

---

### 10.3 Neumann BCs (Natural)

**Application:**
```cpp
SideSet* sideset = field->sideset(2);
sideset->impose_neumann(1000.0);  // Flux q = 1000
```

**Weak form contribution:**
```cpp
// Surface integral: ∫ N^T * q * dS

for (Element* facet_elem : sideset->elements()) {
    calc->link(facet_elem);

    for (uint g = 0; g < num_gauss; ++g) {
        const Matrix<real>& N = calc->N(g);
        real dS = calc->dS(g);
        real q = sideset->bc_value( 0 );  // 1000.0

        f_e += trans(N) * q * dS;
    }

    solver_data->assemble_rhs(facet_elem, f_e);
}
```

**No modification to Jacobian** - only RHS vector.

---

### 10.4 Robin BCs (Convection)

**Application:**
```cpp
SideSet* sideset = field->sideset(3);
real h = 100.0;       // Convection coefficient (W/m²K)
real T_inf = 298.0;   // Ambient temperature (K)

sideset->impose_alpha(h, T_inf);
```

**Weak form contribution:**
```cpp
// Jacobian: ∫ N^T * h * N * dS
// RHS:      ∫ N^T * h * T_inf * dS

for (Element* facet_elem : sideset->elements()) {
    calc->link(facet_elem);
    const Vector<real>& tTinf = mesh->field_data( "Tinf" );  // written by impose_alpha()

    for (uint g = 0; g < num_gauss; ++g) {
        const Matrix<real>& N = calc->N(g);
        real dS = calc->dS(g);
        real h = sideset->bc_value( 0 );   // alpha, stored by impose_alpha()
        real T_inf = tTinf( facet_elem->element()->node( 0 )->index() );  // impose_alpha() writes T_inf per node into the "Tinf" mesh field

        K_e += trans(N) * h * N * dS;         // Jacobian
        f_e += trans(N) * h * T_inf * dS;     // RHS
    }

    solver_data->assemble_jacobian(facet_elem, K_e);
    solver_data->assemble_rhs(facet_elem, f_e);
}
```

**Modifies both Jacobian and RHS.**

---

### 10.5 BC Imposing Modes

**Enum:** `BoundaryConditionImposing`

```cpp
enum class BoundaryConditionImposing
{
    Free,
    Dirichlet,       // nodal value such as displacement or temperature
    Neumann,         // nodal flux, such as heat load or force
    Alpha,           // convective heat flux
    Lambda,          // Maxwell-specific
    Weak,            // Maxwell-specific
    UNDEFINED
};
```

`SideSet::bc_type( aDimension )` returns the entry a sideset carries for that dof dimension;
Dirichlet values are imposed strongly through the JEDI split.

---

## 11. Hanging Nodes and T-Matrices

### 11.1 What Are Hanging Nodes?

Hanging nodes occur at mesh refinement boundaries where a node on a coarse element lies on an edge/face of refined elements.

**Example (1D):**
```
Coarse element:       Refined element:
a-----------b         a-----h-----b
                           ^
                      Hanging node h
```

**Constraint:** `h = 0.5*a + 0.5*b`

**Example (2D quad refinement):**
```
Coarse:              Refined:
a-------b            a---e---b
|       |            | 1 | 2 |
|       |            f---h---g
d-------c            | 3 | 4 |
                     d---i---c
```

**Hanging nodes:** e, f, g, h, i
**Constraints:**
- `e = 0.5*a + 0.5*b`
- `f = 0.5*a + 0.5*d`
- `g = 0.5*b + 0.5*c`
- `h = 0.25*a + 0.25*b + 0.25*c + 0.25*d`
- `i = 0.5*d + 0.5*c`

---

### 11.2 Hanging DOF Data Structure

```cpp
class Dof {
    uint mNumberOfSources;      // Number of parent DOFs
    Dof** mSources;             // Parent DOF pointers (malloc)
    real* mCoefficients;        // Interpolation weights (malloc)
};

// Example: h = 0.25*a + 0.25*b + 0.25*c + 0.25*d
Dof* h = /* hanging node DOF */;
h->mNumberOfSources = 4;
h->mSources = [a, b, c, d];           // Pointers to parent DOFs
h->mCoefficients = [0.25, 0.25, 0.25, 0.25];

// Check if DOF is hanging
if (h->is_hanging()) {
    // Get parent DOFs and weights
    for (uint i = 0; i < h->number_of_sources(); ++i) {
        Dof* parent = h->source(i);
        real weight = h->weight(i);
        // h->value() = sum(weight * parent->value())
    }
}
```

---

### 11.3 T-Matrix Transformation

**Element-level transformation:**

Original element DOFs include hanging nodes:
```
u_elem = [u1, u2, u3, u_h, ...]  // u_h is hanging
```

Transform to free DOFs only:
```
u_elem = T * u_free

where T is the transformation matrix:
T = [ 1   0   0   ... ]  (row 1: u1 = u_free[1])
    [ 0   1   0   ... ]  (row 2: u2 = u_free[2])
    [ 0   0   1   ... ]  (row 3: u3 = u_free[3])
    [ w1  w2  w3  ... ]  (row 4: u_h = w1*u1 + w2*u2 + w3*u3)
```

**Element matrix transformation:**
```
K_elem = T^T * K_free * T
```

This eliminates hanging DOFs from the element matrix before assembly.

---

### 11.4 T-Matrix Workflow

```cpp
// 1. Collect hanging DOFs
DofData::collect_hanging_dofs()
├── Traverse all nodes/edges/faces
├── Check if entity is hanging (mesh::Node::is_hanging())
├── Get parent entities and weights
└── Store in mHangingDOFs container

// 2. Remove from active DOF list
DofData::remove_hanging_dofs_from_container()
├── Remove hanging DOFs from mDOFs
└── Keep in separate mHangingDOFs container

// 3. Create T-matrices for elements
DofData::create_dofwise_t_matrices_master()
├── For each element with hanging DOFs
├── Identify which DOFs are hanging
├── Create T-matrix transformation
└── Store in Element::mTmatrix

// 4. Relink element DOFs
Element::relink_dofs()
├── Replace hanging DOF pointers with parent DOFs
└── Apply T-matrix during assembly
```

---

### 11.5 Assembly with T-Matrices

**Standard assembly (no hanging nodes):**
```cpp
iwg->compute_jacobian(elem, K_e);
solver_data->assemble_jacobian(elem, K_e);
```

**Assembly with T-matrix:**
```cpp
iwg->compute_jacobian(elem, K_e);  // K_e includes hanging DOFs

// Apply T-matrix transformation
if (elem->has_t_matrix()) {
    elem->t_matrix()->project(K_e, K_cond);  // K_cond = Tᵀ K_e T, only free DOFs
    solver_data->assemble_jacobian(elem, K_cond);
} else {
    solver_data->assemble_jacobian(elem, K_e);
}
```

**Automatic in BELFEM:** the projection is done by `DofManager::compute_jacobian()` /
`compute_jacobian_and_rhs()` before `SolverData::assemble_jacobian()` is called.

---

### 11.6 Hanging-Node Counts vs. Elimination (historical note)

The `number_of_hanging_dofs()` accessor bug (§1.1) is **fixed**; both count
accessors are safe. Even while it was live it never affected the hanging-node
*elimination* — assembly, solving and solution accuracy all go through the
T-matrix mechanism, which computes its maps independently of the count
accessors. Only logging and user queries ever saw the wrong number.

---

## 12. Postprocessing (L2 Projection)

### 12.1 What Is L2 Projection?

**Problem:** IWG computes secondary fields (e.g., heat flux, stress) at integration points. These are **discontinuous** across elements.

**Goal:** Project discontinuous field to **continuous** nodal field for visualization.

**Method:** L2 projection (least-squares fit)

**Formulation:**
```
Minimize: ||q_gauss - N*q_nodes||²

Solve: M*q_nodes = ∫ N^T * q_gauss * dV

where M = ∫ N^T * N * dV  (mass matrix)
```

---

### 12.2 Postprocessor Workflow

```cpp
// 1. Create postprocessors
field->initialize_postprocessors();

// This creates separate DofManager for each secondary field
// E.g., for heat conduction:
//   - Temperature (T) - primary field (already in field)
//   - Heat flux (q) - secondary field (needs postprocessor)

// 2. Compute secondary fields
field->postprocess();

// Workflow:
field->postprocess()
├── For each postprocessor
│   ├── 1. Evaluate secondary field at integration points
│   │   └── iwg->compute_secondary_field(elem, gauss_point, value)
│   │
│   ├── 2. Assemble RHS: f = ∫ N^T * value * dV
│   │
│   ├── 3. Solve: M*q = f
│   │   (M is mass matrix, precomputed during initialize_postprocessors())
│   │
│   └── 4. Store nodal values in mesh field
│       └── mesh->field("q")->set_value(node_index, q_value)
│
└── Synchronize projected fields (MPI)
```

---

### 12.3 Creating Postprocessors

```cpp
// Automatic creation based on IWG requirements
field->initialize_postprocessors();

// IWG specifies secondary fields via:
IWG::create_fields(DofManager*)
{
    // Create mesh fields for secondary quantities
    // (fields are scalar per label — one field per component)
    mesh->create_field("heat_flux_x");
    mesh->create_field("gradient_T_x");
}

// DofManager creates postprocessors for each secondary field
// Each postprocessor is a separate mini-DofManager with:
//   - DOFs for secondary field
//   - Mass matrix M
//   - RHS vector f
//   - Solver (usually direct, since M is well-conditioned)
```

---

### 12.4 Accessing Projected Fields

```cpp
// After postprocessing
field->postprocess();

// Access on mesh
mesh::Field* q_field = mesh->field("heat_flux_x");

// Get nodal values
for (mesh::Node* node : mesh->nodes()) {
    real qx = q_field->value(node->index());
}

// Or direct access to data array
Vector<real>& q_data = q_field->data();
// q_data = [qx0, qy0, qz0, qx1, qy1, qz1, ...]
```

---

### 12.5 Postprocessing Performance

**Cost:** L2 projection requires solving `M*q = f` for each secondary field.

**Optimization:**
- Mass matrix M is SPD (Symmetric Positive Definite)
- Can use Cholesky factorization (faster than LU)
- M is constant (compute once, reuse for all timesteps)
- Can use lumped mass matrix for cheaper approximation

**Lumped mass (faster, less accurate):**
```
M_lumped[i,i] = sum_j(M[i,j])  (row sum)
M_lumped[i,j] = 0  (i ≠ j)

// Solve becomes trivial:
q[i] = f[i] / M_lumped[i,i]
```

**BELFEM uses consistent mass by default** (more accurate).

---

## 13. Ownership and Lifetime

### 13.1 Ownership Hierarchy

**BELFEM uses explicit ownership via raw pointers.**

```
Kernel (top-level owner)
├── borrows → Mesh* (the caller who did `new Mesh` deletes it; `mOwnMesh` is not set by any constructor)
├── owns → Materials
├── owns → BoundaryConditions
├── owns → Cell<IWG*>
└── owns → Cell<DofManager*>

DofManager
├── borrows → Kernel* (parent, NOT owned)
├── borrows → Mesh* (NOT owned)
├── borrows → IWG* (NOT owned)
├── owns → Parameters*
├── owns → DofData*
├── owns → BlockData*
├── owns → SideSetData*
├── owns → BearingData*
├── owns → FieldData*
├── owns → SolverData*
├── owns → EigenValues*
└── owns → Cell<Postprocessor*>

DofData
├── owns → Cell<Dof*>
├── owns → Cell<Dof*> (hanging DOFs)
└── owns → Map<id_t, Dof*> (map does NOT own Dof objects)
```

**Key Rules:**
- `owns →` means destructor calls `delete`
- `borrows →` means destructor does NOT call `delete`
- `Cell<T*>` is non-owning (does NOT call `delete` on elements)
- Ownership is documented in header comments

---

### 13.2 Lifetime Management Patterns

**Pattern 1: Kernel Manages Everything (RECOMMENDED)**

```cpp
int main() {
    Mesh* mesh = new Mesh("mesh.exo");

    // Kernel owns everything downstream of the parameters
    KernelParameters params(mesh);
    Kernel kernel(&params);

    // Kernel creates and owns
    IWG* iwg = kernel.create_equation(...);
    DofManager* field = kernel.create_field(iwg);

    // ... use field ...

    // Kernel destructor deletes:
    //   - DofManager (which deletes 8 data managers)
    //   - IWG
    //   - Materials
    //   - BoundaryConditions

    delete mesh;  // the mesh stays ours

    return 0;
}
```

**Pattern 2: Manual Management (ADVANCED)**

```cpp
// Create kernel
Kernel* kernel = new Kernel(&params);

// Create IWG
IwgFactory factory(mesh);
IWG* iwg = factory.create_iwg(...);

// Create DofManager
DofManager* field = new DofManager(kernel, 0);
field->set_equation(iwg);

// ... use field ...

// Manual cleanup (CORRECT ORDER!)
delete field;   // 1. Delete DofManager first
delete iwg;     // 2. Then delete IWG
delete kernel;  // 3. Finally delete Kernel

// ❌ WRONG ORDER:
// delete iwg;     // DofManager still holds pointer!
// delete field;   // Dangling pointer access!
```

**Pattern 3: DofManager Owns IWG (ALTERNATIVE)**

```cpp
IwgFactory factory(mesh);
IWG* iwg = factory.create_iwg(...);

DofManager* field = kernel.create_field(iwg);

// ⚠️ If DofManager created via factory pattern,
// it may take ownership of IWG!
// Check implementation or header comments.

delete field;  // May delete IWG automatically
// DO NOT: delete iwg;  // Potential double-free!
```

**⚠️ CRITICAL:** Always check header comments for ownership documentation!

---

### 13.3 Common Ownership Mistakes

**Mistake 1: Double-Free**
```cpp
IWG* iwg = kernel.create_equation(...);
DofManager* field = kernel.create_field(iwg);

delete iwg;  // ❌ WRONG: Kernel still owns it!

// Later: kernel destructor tries to delete iwg → CRASH!
```

**Mistake 2: Dangling Pointer**
```cpp
IWG* iwg = factory.create_iwg(...);
DofManager* field = new DofManager(&kernel, 0);
field->set_equation(iwg);

delete iwg;  // ❌ WRONG: DofManager still holds pointer!

field->compute_jacobian_and_rhs();  // Dangling pointer access!
```

**Mistake 3: Memory Leak**
```cpp
for (int i = 0; i < 100; ++i) {
    IWG* iwg = factory.create_iwg(...);
    // ❌ Never deleted → 100 leaked IWG objects!
}
```

**Mistake 4: Deleting Borrowed Pointer**
```cpp
DofManager* field = kernel.create_field(iwg);
Mesh* mesh = field->mesh();  // Borrowed from Kernel

delete mesh;  // ❌ WRONG: Kernel owns mesh!

// Later: kernel destructor tries to delete mesh → CRASH!
```

---

### 13.4 Smart Pointer Policy

**BELFEM uses raw pointers on critical paths** for performance reasons.

**When smart pointers are OK:**
- High-level utilities (I/O, post-processing)
- One-off allocations (setup, initialization)
- Outside hot paths (assembly, solve)

**When raw pointers are required:**
- Inside loops (assembly, element traversal)
- Interfacing with C libraries (MPI, BLAS, Fortran)
- Performance-critical code

**See:** `doc/coding_philosophy.md` for complete rationale.

---

### 13.5 Checking Ownership

**Header comments:**
```cpp
class DofManager {
    Kernel* mParent;          //! the parent object (borrowed, NOT owned)
    Mesh* mMesh;              //! the mesh (borrowed, NOT owned)
    IWG* mIWG;                //! the equation (borrowed, NOT owned)
    DofData* mDofData;        //! owned, deleted in destructor
};
```

**Ownership documentation pattern:**
- `//! owned` - Destructor calls `delete`
- `//! borrowed` - Destructor does NOT call `delete`
- `//! the parent object` - Usually borrowed (child does not outlive parent)

**If unclear, check destructor:**
```cpp
DofManager::~DofManager()
{
    delete mDofData;       // ✅ Owned (deleted)
    delete mBlockData;     // ✅ Owned (deleted)
    // ...

    // mParent NOT deleted   ✅ Borrowed (not deleted)
    // mMesh NOT deleted     ✅ Borrowed (not deleted)
    // mIWG NOT deleted      ✅ Borrowed (not deleted)
}
```

---

## 14. Performance Considerations

### 14.1 Computational Complexity

**DOF Allocation:**
- O(N_entities × N_dof_types) - Linear in number of mesh entities
- Dominated by: Graph extraction O(N_elements × N_dofs_per_elem)

**Graph Reordering:**
- RCM: O(N_dofs + E_graph + N_dofs × log(N_dofs)) for sorting
- METIS: O(N_dofs × log(N_dofs)) expected

**Matrix Assembly:**
- O(N_elements × N_dofs_per_elem² × N_gauss)
- Dominated by: Element matrix computation (IWG::compute_jacobian)

**Sparse Solver:**
- Direct (MUMPS, STRUMPACK): O(N_dofs^1.5 to N_dofs^2) depending on fill-in
- Iterative (PETSc): O(iterations × N_dofs × avg_connectivity)

**Overall:** Assembly is usually bottleneck for large problems.

---

### 14.2 Memory Usage

**DOF Storage:**
```
Per DOF: ~200 bytes (Dof object + graph::Vertex overhead + hanging node data)
Total: N_dofs × 200 bytes

Example: 1M DOFs → ~200 MB
```

**Sparse Matrix Storage (CSR format):**
```
J matrix: (N_free × N_free) with avg_connectivity entries per row
Storage: N_free × avg_connectivity × sizeof(real) + indexing overhead

Example: 1M DOFs, avg 50 connections per row
  Values: 1M × 50 × 8 bytes = 400 MB
  Indices: 1M × 50 × 4 bytes = 200 MB
  Total: ~600 MB
```

**JEDI Matrices:**
```
J: ~600 MB (from above)
D: N_free × N_fixed × sparsity ~ 10-50 MB (typically small)
E: N_fixed × N_free × sparsity ~ 10-50 MB
I: N_fixed × N_fixed × sparsity ~ 1-10 MB (typically very sparse)

Total JEDI: ~700 MB (dominated by J)
```

**Peak Memory:**
- During assembly: Element matrices (small, reused)
- During solve: Factorization (can be 2-10× matrix storage for direct solvers)

**Example: 1M DOF problem**
- DOFs: 200 MB
- JEDI matrices: 700 MB
- Factorization: 2-5 GB (MUMPS with fill-in)
- **Total: ~3-6 GB**

---

### 14.3 Performance Optimization Tips

**1. Reuse DofManager:**
```cpp
// ❌ BAD: Create inside loop
for (int step = 0; step < n_steps; ++step) {
    DofManager* field = kernel.create_field(iwg);  // Expensive!
    field->initialize();
    // ...
    delete field;
}

// ✅ GOOD: Create once, reuse
DofManager* field = kernel.create_field(iwg);
field->initialize();  // Once

for (int step = 0; step < n_steps; ++step) {
    field->compute_jacobian_and_rhs();
    field->solve();
}
```

**2. Use Optimal Reordering:**
```cpp
SolverParameters params( SolverType::MUMPS );

// For direct solvers (MUMPS, STRUMPACK)
params.set_reordering_method( ReorderingMethod::METIS );  // ✅ Nested dissection
```

**3. Choose Right Solver:**
```cpp
// Small-medium problems (< 100k DOFs)
SolverParameters params( SolverType::UMFPACK );  // Fast for small problems

// Large problems, serial
SolverParameters params( SolverType::MUMPS );    // Robust, good memory management

// Large problems, parallel
SolverParameters params( SolverType::STRUMPACK );  // ~2× faster than MUMPS (if available)

// Very large, well-conditioned problems
SolverParameters params( SolverType::PETSc );    // Iterative, scales to millions of DOFs
```

**4. Minimize Field Synchronization:**
```cpp
// ❌ BAD: Synchronize every iteration
for (int iter = 0; iter < max_iter; ++iter) {
    field->compute_jacobian_and_rhs();
    field->solve();
    field->synchronize_fields({"T"});  // Expensive MPI call!
}

// ✅ GOOD: Synchronize only when needed
for (int iter = 0; iter < max_iter; ++iter) {
    field->compute_jacobian_and_rhs();
    field->solve();
}
field->synchronize_fields({"T"});  // Once at end
```

**5. Preallocate Work Arrays:**
```cpp
// Inside IWG implementation
class IWG_Custom : public IWG {
    Matrix<real> mWorkK;  // Preallocated work array
    Vector<real> mWorkF;

    void compute_jacobian(Element* elem, Matrix<real>& K) {
        // Reuse mWorkK instead of allocating inside loop
        mWorkK.fill(0.0);
        // ... compute ...
        K = mWorkK;
    }
};
```

---

### 14.4 Profiling and Benchmarking

**Use BELFEM profiling tools:**
```cpp
#include "cl_Timer.hpp"
#include "cl_Profiler.hpp"

Timer timer;

// Measure initialization
timer.reset();
field->initialize();
message(InfoLevel::Verbose, "Initialize: %u ms", (uint) timer.stop());

// Measure assembly
timer.reset();
field->compute_jacobian_and_rhs();
message(InfoLevel::Verbose, "Assembly: %u ms", (uint) timer.stop());

// Measure solve
timer.reset();
field->solve();
message(InfoLevel::Verbose, "Solve: %u ms", (uint) timer.stop());
```

**Typical time breakdown:**
- **Initialize:** 5-10% (one-time cost)
- **Assembly:** 30-50% (per Newton iteration)
- **Solve:** 40-60% (per Newton iteration)
- **Postprocess:** 5-10% (per output step)

**Optimization Priority:**
1. Reduce Newton iterations (better initial guess, line search)
2. Optimize assembly (vectorize element loops, reduce function calls)
3. Optimize solver (choose best reordering, preconditioner)

---

## 15. Thread Safety and MPI

### 15.1 Thread Safety

**BELFEM is deliberately NOT thread-safe internally.**

**Rationale:**
- MPI for parallelism (distributed memory), not threading
- No internal mutexes (adds latency to every operation)
- Users needing OpenMP: protect calls externally

**Example: OpenMP with BELFEM**
```cpp
// ❌ WRONG: Race condition
#pragma omp parallel for
for (int i = 0; i < n_elements; ++i) {
    field->compute_jacobian(elem[i]);  // Multiple threads accessing same DofManager!
}

// ✅ CORRECT: Protect with critical section
#pragma omp parallel for
for (int i = 0; i < n_elements; ++i) {
    #pragma omp critical
    {
        field->compute_jacobian(elem[i]);
    }
}

// ⚠️ But this defeats parallelism! Better approach:
// Use MPI for parallelism, not OpenMP
```

**Safe for OpenMP:**
- Read-only access to DOF values
- Independent element computations (no shared state)

**NOT safe for OpenMP:**
- Matrix assembly (writes to shared sparse matrix)
- DOF value modifications
- Graph operations (reordering, connectivity)

---

### 15.2 MPI Collective Operations

**Collective operations require ALL procs to participate:**

```cpp
// ✅ CORRECT: All procs call
field->initialize();              // Collective
field->solve();                   // Collective
field->synchronize_fields({"T"}); // Collective
```

**⚠️ DEADLOCK if only master calls:**
```cpp
// ❌ WRONG: Deadlock!
if (comm_rank() == 0) {
    field->initialize();  // Workers waiting for MPI_Bcast!
}
```

**Master-only operations:**
```cpp
// ✅ CORRECT: Only master writes
if (comm_rank() == 0) {
    mesh->save("output.exo");  // File I/O is master-only
}
```

---

### 15.3 MPI Communication Patterns

**Point-to-point:**
```cpp
// Worker → Master (assembly contributions)
MPI_Send(K_local, master_rank, ...);

// Master → Worker (distributed solution)
MPI_Recv(x_local, worker_rank, ...);
```

**Collective:**
```cpp
// Gather: All → Master
MPI_Gather(local_data, global_data, master_rank, ...);

// Scatter: Master → All
MPI_Scatter(global_data, local_data, master_rank, ...);

// Broadcast: Master → All
MPI_Bcast(shared_data, master_rank, ...);

// Reduce: All → Master (sum, max, min, etc.)
MPI_Reduce(local_value, global_value, MPI_SUM, master_rank, ...);
```

**BELFEM abstracts MPI via `Communicator` class:**
```cpp
#include "commtools.hpp"

// Rank and size. The free functions are the API to use; gComm.rank() and
// gComm.size() also exist, but nothing creates or hands out a communicator.
proc_t rank = comm_rank();
proc_t size = comm_size();

// Point-to-point. There is no tag parameter: the tag is derived from the
// rank pair by comm_tag(), so a send and its receive match automatically.
if (comm_rank() == 0) send(data, 1);       // to rank 1
if (comm_rank() == 1) receive(data, 0);    // from rank 0

// Collective
broadcast(data, 0);                        // every rank calls it
comm_barrier();
```

Two things worth knowing before you reach for something that is not there:

- **There is no global-sum helper.** `allreduce()` exists but is fixed to `MPI_MAX`
  (`commtools.hpp:234-250`). A global sum has to be written with `MPI_Allreduce` directly; if
  you need one more than once, add it to `commtools.hpp` rather than open-coding it at the call
  site.
- **`share`/`receive` and `distribute`/`collect` are not collective.** They are asymmetric
  halves that must be paired across ranks, and calling one without its counterpart deadlocks.
  See the [communication module guide](../../../comm/doc/comm_usage_guide.md).

---

### 15.4 MPI-Safe Field Access

```cpp
// After solve, synchronize before accessing field values
field->solve();
field->synchronize_fields({"T"});

// Now safe to access on all procs
for (mesh::Node* node : mesh->owned_nodes()) {
    real T = mesh->field("T")->value(node->index());
    // ... use T ...
}

// ⚠️ Aura nodes may have stale data without synchronization!
for (mesh::Node* node : mesh->aura_nodes()) {
    real T = mesh->field("T")->value(node->index());  // May be stale!
}
```

---

## 16. Known Issues and Bugs

### 16.1 FIXED: `number_of_hanging_dofs()` (reported 2026-01-20)

The accessor once returned `mNumberOfFixedDofs`; it now returns
`mNumberOfHangingDofs` (`cl_FEM_DofMgr_DofData.hpp`). No open issue remains.
Kept as the record of the report; the standing lesson is in §1.1.

---

### 16.2 (retracted)

An earlier revision of this guide reported a typo in the `assemble_rhs()`
method name; the report itself was garbled (it showed the same spelling on
both sides) and the method exists as `assemble_rhs()` in
`cl_FEM_DofMgr_SolverData.hpp`. Nothing to work around.

---

### 16.3 CBorg False Positives (Verified Not Bugs)

During AI-assisted code review, CBorg (having access only to `.hpp` header files, not `.cpp` implementations) flagged **7 issues** that were verified to be **false positives**:

**❌ FALSE POSITIVE #1: DOF offsets not initialized before use**
- **Claim:** `mEdgeDofOffset`, `mFaceDofOffset`, etc. are `gNoID` and never set
- **Reality:** `compute_dof_offsets(aIWG)` is called at the START of `create_dofs()` (before any DOF ID calculations)
- **Evidence:** `compute_dof_offsets()` call at the top of `DofData::create_dofs()` in `cl_FEM_DofMgr_DofData.cpp`

**❌ FALSE POSITIVE #2: `create_dofwise_t_matrices_master()` declared but not defined**
- **Claim:** Function is declared in header but never implemented
- **Reality:** Full ~100-line implementation exists
- **Evidence:** `DofData::create_dofwise_t_matrices_master()` in `cl_FEM_DofMgr_DofData.cpp`

**❌ FALSE POSITIVE #3: Destructor order issues**
- **Claim:** DofManager destructor deletes in wrong order or misses objects
- **Reality:** Destructor correctly deletes all 8 owned data managers in proper order
- **Evidence:** `DofManager::~DofManager()` in `cl_FEM_DofManager.cpp`

**❌ FALSE POSITIVE #4: `node_dof_id()` formula breaks for sparse meshes**
- **Claim:** Formula `nodeID * mNumDofTypes + dofType` causes collisions with sparse node IDs
- **Reality:** Offset system ensures non-overlapping ranges (nodes, edges, faces, cells, lambdas each have separate ranges)
- **Explanation:** May waste ID space but cannot cause collisions. Even node ID 1,000,000 produces valid ID in node range.

**❌ FALSE POSITIVE #5: `EigenValues::mSolver` not initialized**
- **Claim:** `mSolver` pointer is never set before use
- **Reality:** Lazy initialization pattern using `mFirstRun` flag
- **Evidence:** `EigenValues::mFirstRun` in `cl_FEM_DofMgr_EigenValues.cpp` - `mSolver` created on first use

**❌ FALSE POSITIVE #6: `bearing(id)` returns nullptr**
- **Claim:** `mEmptyBearing` is nullptr, causing crash
- **Reality:** `mEmptyBearing = new Bearing(mParent)` is always allocated in `create_bearings()`
- **Evidence:** `BearingData::create_bearings()` in `cl_FEM_DofMgr_BearingData.cpp` - Returns valid (empty) Bearing, not nullptr

**❌ FALSE POSITIVE #7: Calculator function pointers null before allocation**
- **Claim:** Function pointers (`mFunN`, `mFunB`, etc.) can be null when accessed
- **Reality:** Guard exists in `link()` method: `BELFEM_ERROR(mIsAllocated, "calculator is not allocated")`
- **Evidence:** `Calculator::link()` in `cl_FEM_Calculator.cpp` - Must call `link()` before any accessor, `link()` checks allocation

---

**⚠️ PARTIALLY VALID (Defensive Programming Recommendations):**

1. **`solve()` has no explicit solver null check**
   - No explicit null check for `mSolver` before `mSolver->solve()` (line 17893)
   - Normal usage always calls `set_solver()` first, so unlikely to occur
   - **Recommendation:** Add `BELFEM_ASSERT(mSolver != nullptr, "Solver not set")` for safety

2. **Calculator accessors could have explicit allocation checks**
   - Guard exists in `link()` but not in individual accessors (`N()`, `B()`, etc.)
   - Since `link()` must be called before any accessor, current design is safe
   - **Recommendation:** Could add `BELFEM_ASSERT(mIsAllocated, ...)` to individual accessors for extra safety

---

**Lesson Learned:**

Header-only review (`.hpp` files) can identify:
- ✅ Interface design issues
- ✅ Declaration/definition mismatches (when both are in headers)
- ✅ Inconsistent naming

But may **incorrectly assume** missing functionality that actually exists in `.cpp` implementation files:
- ❌ Missing initializations (may be in constructors)
- ❌ Missing function definitions (may be in .cpp)
- ❌ Incorrect execution order (may be correct in implementation)

**Verification Method:** Always check `.cpp` implementation files before concluding that functionality is missing!

---

### 16.4 Debugging Tips

**Enable verbose logging:**
```cpp
gLog.set_info_level( InfoLevel::Verbose );

// Now see detailed initialization messages
field->initialize();
// Output:
//   Creating 10000 node DOFs...
//   Creating 5000 edge DOFs...
//   Reordering with METIS...
//   Allocating Jacobian: 50000 × 50000, 2.5M non-zeros...
```

**Print DOF information:**
```cpp
field->print(0);  // Print for master proc (rank 0)

// Output:
//   DOF Manager 0:
//     Free DOFs: 45000
//     Fixed DOFs: 5000
//     Hanging DOFs: 1000
//     Blocks: 3
//     Sidesets: 5
```

**Check matrix dimensions:**
```cpp
SpMatrix* A = field->system_matrix();
index_t n_free = field->solver_data()->number_of_free_dofs();

BELFEM_ASSERT(A->n_rows() == n_free, "A should be square");
BELFEM_ASSERT(A->n_cols() == n_free, "A should be square");

message(InfoLevel::Info, "J: %lu × %lu, nnz = %lu",
    J->n_rows(), J->n_cols(), J->number_of_nonzeros());
```

**Check DOF connectivity:**
```cpp
Dof* dof = field->dof(dof_id);

message(InfoLevel::Info, "DOF %lu:", dof->id());
message(InfoLevel::Info, "  Index: %lu", dof->index());
message(InfoLevel::Info, "  Fixed: %s", dof->is_fixed() ? "yes" : "no");
message(InfoLevel::Info, "  Hanging: %s", dof->is_hanging() ? "yes" : "no");
message(InfoLevel::Info, "  Connectivity: %u neighbors", dof->number_of_vertices());

for (uint i = 0; i < dof->number_of_vertices(); ++i) {
    Dof* neighbor = (Dof*) dof->vertex(i);
    message(InfoLevel::Info, "    Neighbor %u: DOF %lu", i, neighbor->id());
}
```

**Valgrind for memory issues:**
```bash
valgrind --leak-check=full --track-origins=yes ./my_belfem_app

# Look for:
#   - Invalid reads/writes (use-after-free, out-of-bounds)
#   - Memory leaks (missing delete)
#   - Uninitialized values
```

---

## 17. Development Notes

### 17.1 Adding New BC Types

When implementing new boundary condition types:

1. **Add to enum** (`en_FEM_BoundaryConditionImposing.hpp` — the per-dimension kind a sideset stores):
   ```cpp
   enum class BoundaryConditionImposing {
       Free,
       Dirichlet,
       Neumann,
       Alpha,
       Lambda,
       Weak,
       CustomBC,  // New BC type
       UNDEFINED
   };
   ```

2. **Update SideSet** (`cl_FEM_SideSet.{hpp,cpp}`):
   ```cpp
   // Header
   void impose_custom_bc(const real aValue, const uint aDofType = 0);

   // Implementation
   void SideSet::impose_custom_bc(const real aValue, const uint aDofType)
   {
       mBcValues(aDofType) = aValue;
       mBcTypes(aDofType) = BoundaryConditionImposing::CustomBC;
   }
   ```

3. **Implement assembly** in IWG:
   ```cpp
   void IWG::compute_surface_loads(Element* elem, Vector<real>& f)
   {
       SideSet* sideset = (SideSet*) elem->parent();

       if (sideset->bc_type(0) == BoundaryConditionImposing::CustomBC) {
           // Compute custom BC contribution
           real bc_value = sideset->bc_value(0);
           // ... weak form assembly ...
       }
   }
   ```

4. **Test:**
   - Manufactured solution with known analytical result
   - Mesh convergence study (h-refinement)
   - MPI parallel consistency check

---

### 17.2 Adding New DOF Entity Types

Currently supported: Nodes, Edges, Faces, Cells, Lambdas (facets).

To add new entity type (e.g., "Volumes" for 3D interior DOFs):

1. **Update IWG interface** (`cl_IWG.hpp`):
   ```cpp
   virtual uint number_of_dofs_per_volume( const id_t aBlockID ) const { return 0; }
   ```

2. **Update DofData** (`cl_FEM_DofMgr_DofData.{hpp,cpp}`):
   ```cpp
   // Header
   id_t mVolumeDofOffset = gNoID;

   id_t volume_dof_id(const id_t aVolumeID, const uint aDofType) const;

   // Implementation
   index_t count_volume_dofs(IWG* aIWG, ...);

   void compute_dof_offsets(IWG* aIWG)
   {
       // ... existing offsets ...
       mVolumeDofOffset = mCellDofOffset + max_cell_id * mNumDofTypes;
   }

   id_t DofData::volume_dof_id(const id_t aVolumeID, const uint aDofType) const
   {
       return mVolumeDofOffset + aVolumeID * mNumDofTypes + aDofType;
   }
   ```

3. **Update Element** (`cl_FEM_Element.{hpp,cpp}`):
   ```cpp
   // Link volume DOFs to element
   void link_volume_dofs();
   ```

4. **Test thoroughly:**
   - DOF ID uniqueness (no collisions)
   - Mesh convergence (verify DOF count scaling)
   - Parallel consistency (all procs agree on DOF IDs)

---

### 17.3 Debugging DOF Issues

**Common DOF problems:**

1. **DOF ID collisions:**
   ```cpp
   // Check for duplicates
   Map<id_t, uint> dof_id_count;
   for (Dof* dof : dofs) {
       dof_id_count(dof->id())++;
   }

   for (auto& pair : dof_id_count) {
       BELFEM_ERROR(pair.second == 1,
           "DOF ID %lu has %u instances (should be 1)",
           pair.first, pair.second);
   }
   ```

2. **Missing DOFs:**
   ```cpp
   // Verify expected count
   index_t expected_node_dofs = mesh->number_of_nodes() * num_dof_types;
   index_t actual_node_dofs = /* count DOFs with ID < edge_offset */;

   BELFEM_ERROR(expected_node_dofs == actual_node_dofs,
       "Expected %lu node DOFs, found %lu",
       expected_node_dofs, actual_node_dofs);
   ```

3. **Hanging node constraint errors:**
   ```cpp
   // Verify hanging node constraint
   for (Dof* h : hanging_dofs) {
       real h_value = 0.0;
       for (uint i = 0; i < h->number_of_sources(); ++i) {
           h_value += h->weight(i) * h->source(i)->value();
       }

       real error = std::abs(h_value - h->value());
       BELFEM_ERROR(error < 1e-10,
           "Hanging node constraint violated: error = %e", error);
   }
   ```

---

## 18. Literature References

### 18.1 Finite Element Theory

**Core FEM Textbooks:**

- **Zienkiewicz & Taylor**, "The Finite Element Method" Vol. 1-2
  - Vol. 1: DOF management, assembly, Chapter 2-3
  - Vol. 2: Advanced topics, Chapter 8-10

- **Hughes**, "The Finite Element Method" (2000)
  - Chapter 1-3: FEM fundamentals
  - Chapter 4: Locking, mixed methods
  - Chapter 5: Element technology

- **Bathe**, "Finite Element Procedures" (2014)
  - Chapter 3: Formulation of FEM
  - Chapter 8: Solvers and algorithms
  - Chapter 9: Nonlinear analysis

**Hanging Nodes and Constraints:**

- **Brenner & Scott**, "The Mathematical Theory of Finite Element Methods" (2008)
  - Chapter 5: Hanging node constraints
  - Chapter 9: Multi-level methods

---

### 18.2 Sparse Solvers

**Direct Solvers:**

- **Davis**, "Direct Methods for Sparse Linear Systems" (2006)
  - UMFPACK implementation details
  - Sparse matrix storage formats (CSR, CSC)
  - Fill-minimizing orderings

- **MUMPS User Guide** - [http://mumps.enseeiht.fr/](http://mumps.enseeiht.fr/)
  - Multifrontal method
  - Parallel sparse factorization

- **STRUMPACK Documentation** - [https://portal.nersc.gov/project/sparse/strumpack/](https://portal.nersc.gov/project/sparse/strumpack/)
  - Structured matrix methods
  - HSS compression

**Iterative Solvers:**

- **Saad**, "Iterative Methods for Sparse Linear Systems" (2003)
  - Chapter 6: Krylov subspace methods
  - Chapter 10: Preconditioning

- **PETSc User Manual** - [https://petsc.org/](https://petsc.org/)
  - KSP (Krylov Subspace) solvers
  - PC (Preconditioner) options

---

### 18.3 Graph Partitioning and Reordering

**Algorithms:**

- **Cuthill & McKee** (1969), "Reducing the Bandwidth of Sparse Symmetric Matrices"
  - ACM Conference Proceedings
  - Reverse Cuthill-McKee (RCM) algorithm

- **George** (1973), "Nested Dissection of a Regular Finite Element Mesh"
  - SIAM J. Numer. Anal., Vol. 10, No. 2
  - Fill-minimizing ordering

- **Liu & Sherman** (1976), "Comparative Analysis of the Cuthill-McKee and the Reversed Cuthill-McKee Ordering Algorithms"
  - SIAM J. Numer. Anal., Vol. 13, No. 2

**Software:**

- **METIS** - [http://glaros.dtc.umn.edu/gkhome/metis/metis/overview](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
  - Karypis & Kumar (1998), "A Fast and High Quality Multilevel Scheme for Partitioning Irregular Graphs"
  - Nested dissection, graph partitioning

- **SCOTCH** - [https://www.labri.fr/perso/pelegrin/scotch/](https://www.labri.fr/perso/pelegrin/scotch/)
  - Alternative to METIS
  - PT-SCOTCH for parallel partitioning

---

### 18.4 MPI and Parallel Computing

**MPI Standards:**

- **MPI Forum**, "MPI: A Message-Passing Interface Standard" (2021)
  - Version 4.0
  - [https://www.mpi-forum.org/](https://www.mpi-forum.org/)

**Parallel FEM:**

- **Gropp et al.**, "Using MPI" (3rd Edition, 2014)
  - Practical MPI programming
  - Collective operations, derived datatypes

- **Bathe & Wilson**, "Numerical Methods in Finite Element Analysis" (1976)
  - Chapter 8: Parallel solution strategies
  - Domain decomposition methods

---

### 18.5 BELFEM-Specific References

**Related BELFEM Modules:**

- **IWG Module** (`src/fem/iwg/doc/`): Physics layer, weak form assembly, time-stepping
- **Interpolation Module** (`src/fem/interpolation/doc/`): Shape functions, integration points, Nédélec elements
- **Sparse Module** (`src/sparse/doc/`): Sparse matrix storage, solver interfaces
- **Graph Module** (`src/math/graph/doc/`): Graph algorithms (BFS, DFS, RCM, METIS)
- **Mesh Module** (`src/mesh/doc/`): Mesh data structures, I/O formats

**BELFEM Papers:**

- **messe2023.txt** - BELFEM core paper (Messe) - Architecture, h-φ formulation
- **arsenault2023.txt** - Magnetodynamic h-φ coupling (Arsenault) - JEDI system usage
- **alves2022b.txt** - Thin-shell theory and cohomology cuts - Hanging nodes, constraints

**See:** `doc/literature_references.md` for complete reference list with DOIs.

---

## 19. H-Phi Mixed Formulation (Electromagnetics)

**Note:** This section covers Maxwell-specific applications of the DOF Manager. For general DOF management, see sections 1-18.

### 19.1 Overview

BELFEM implements a sophisticated **H-Phi mixed formulation** (`maxwell::Formulation::HPhi`) for electromagnetic simulations that combines:

1. **H-formulation** in conductors: Magnetic field intensity **H** as primary variable (edge DOFs)
2. **Phi-formulation** in air/ferro: Scalar potential **φ** as primary variable (node DOFs)

**Why Mixed Formulation?**
- **H-formulation** directly captures eddy currents in conductors (∇ × H = J)
- **Phi-formulation** is computationally efficient in non-conducting regions
- **Thin shells** naturally couple both formulations at interfaces

---

### 19.2 Governing Equations

| Domain Type | Primary Variable | Governing Equation | Weak Form |
|-------------|------------------|-------------------|-----------|
| **Conductor** | **H** (edge DOF) | ∇ × (ρ ∇ × **H**) + μ₀ ∂**H**/∂t = 0 | ∫ (ρ ∇ × **w**) · (∇ × **H**) dV + ∫ μ₀ **w** · ∂**H**/∂t dV = 0 |
| **Air/Ferro** | **φ** (node DOF) | ∇ · (μ ∇φ) = 0, **H** = -∇φ | ∫ (μ ∇w) · (∇φ) dV = 0 |

**Where:**
- ρ = electrical resistivity [Ω·m]
- μ = magnetic permeability [H/m] (μ₀ in air, μ(H,T) in ferro)
- **w** = edge test functions (Nédélec)
- w = nodal test functions (Lagrange)

**References:**
- **Messe et al. 2023** - BELFEM core paper (SUST) - Equations 6-8
- **Arsenault et al. 2023** - Magnetodynamic h-φ coupling (IEEE TASC) - Section II, Eq. 5-7.
  Read together with the **2026 erratum** (DOI:10.1109/TASC.2026.3686487), which corrects the
  air-domain strong and weak forms of that paper

---

### 19.3 Domain Types and DOF Assignments

**Domain Type Enum:** `src/mesh/en_DomainType.hpp` (excerpt)

```cpp
enum class DomainType {
    Air                    =  1,   // φ-formulation (node DOFs)
    Buffer                 =  2,
    Ferro                  =  3,   // φ-formulation with μ(H,T)
    Coil                   =  4,   // Simplified H-formulation (source)
    Conductor              =  5,   // H-formulation (edge DOFs)
    ThinShell              =  9,   // Mixed: edge_h + φ interface

    InterfaceCondAir       = 25,   // Conductor-Air (currently disabled)
    InterfaceCondFerro     = 26,   // Conductor-Ferro (H-Phi only)
    InterfaceFerroAir      = 27,   // Ferro-Air (φ-φ coupling)
};
```

**DOF Field Assignments:**

| Domain | DOF Fields | Entity | Physics |
|--------|-----------|--------|---------|
| **Air** | `"phi"` | Node | H = -∇φ |
| **Ferro** | `"phi"` | Node | H = -∇φ, μ = μ(H,T) |
| **Conductor** | `"edge_h"` | Edge | ∇ × H = J |
| (Higher order) | `+ "face_h"` | Face | Second-order edge elements |
| **ThinShell** | `"edge_h"` + `"phi"` | Edge + Node | **Both formulations!** |

**Facet Master Priority** (for interface orientation):
```
Conductor > Ferro > Air > Coil
```

This ensures consistent normal vectors at interfaces.

---

### 19.4 Field List Construction for Thin Shells

**Reference:** `cl_Maxwell_FieldList.cpp:116-156`

Thin shells get DOFs from **both** conductor and air domains:

```cpp
// Conductor DOFs added to ThinShell (line 124)
for (const string & tDof : Conductor) {
    ThinShell.push(tDof);  // Gets "edge_h" (and "face_h" if higher order)
}

// Air DOFs also added to ThinShell (line 151)
for (const string & tDof : Air) {
    ThinShell.push(tDof);  // Gets "phi"
}

// Result: ThinShell = {"edge_h", "phi"} (+ "face_h" if higher order)
```

This is **automatic** - users don't need to manually specify thin shell DOFs.

---

### 19.5 Hanging DOFs in Electromagnetic Context: Cohomology Cuts

In electromagnetic problems, hanging DOFs arise from enforcing **current boundary conditions** via cohomology cuts.

**Problem:** Impose electric current `I` through a conductor in a multiply-connected non-conducting domain.

**Solution:** Create a "cut" surface where the scalar potential φ has a discontinuity:
```
[φ] = φ⁺ - φ⁻ = I
```

**Implementation:**

1. **Duplicate nodes** on one side of the cut surface
2. **Original nodes** → independent DOFs (φ⁻)
3. **Duplicated nodes** → hanging DOFs (φ⁺)
4. **T-matrix** enforces: φ⁺ = φ⁻ + I

**Creation:** `MaxwellFactory::create_hanging_edges_and_facets()` at line 1008 in `cl_MaxwellFactory.cpp`

**Collection:** `DofData::collect_hanging_dofs()` in `cl_FEM_DofMgr_DofData.cpp`

**T-Matrix Setup:** `DofData::create_dofwise_t_matrices_master()` in `cl_FEM_DofMgr_DofData.cpp`

**Workflow:**
```
1. Mesh Creation
   └─> MaxwellFactory::create_hanging_edges_and_facets()
       └─> Creates hanging edges/facets at interfaces

2. DOF Initialization
   └─> DofData::create_dofs()
       └─> DofData::collect_hanging_dofs()
           └─> Identifies all hanging DOFs
           └─> DofData::create_dofwise_t_matrices_master()
               └─> Computes T-matrices for each hanging DOF

3. Matrix Assembly
   └─> DofManager::compute_jacobian()
       └─> For each element:
           └─> if has_t_matrix(): K_condensed = T^T * K_element * T

4. System Solution
   └─> Solve reduced system (independent DOFs only)
   └─> Hanging DOF values computed via: φ_hanging = Σ w_i * φ_source_i
```

**References:**
- **Alves et al. 2022b** - Thick cuts / cohomology basis (Equation 10)
- **Riva et al. 2023** - Cohomology cuts for transport current (IEEE TASC)
- **Schnaubelt et al. 2023** - Thick cuts interpretation (IEEE TASC)

**See also:** [Section 11: Hanging Nodes and T-Matrices](#11-hanging-nodes-and-t-matrices) for general mesh refinement hanging nodes.

---

### 19.6 Static Condensation vs. Lagrange Multipliers

**BELFEM uses static condensation exclusively** for all interface coupling, including:
- Thin shell h-φ interfaces
- Cohomology cuts for transport current
- Hanging DOFs from duplicated nodes

**Why Static Condensation?** (following Messe et al. 2023's recommendation)

| Aspect | Static Condensation | Lagrange Multipliers |
|--------|---------------------|---------------------|
| **Diagonal zeros** | ✅ Avoids zeros on main diagonal | ❌ Introduces zeros (saddle-point) |
| **Solver robustness** | ✅ Better for direct solvers (MUMPS/STRUMPACK) | ⚠️ Requires specialized saddle-point solvers |
| **Global DOF count** | ✅ Reduces DOF count | ❌ Increases DOF count (adds λ DOFs) |
| **Matrix structure** | ✅ Positive definite (or close) | ❌ Indefinite saddle-point |
| **Implementation** | ✅ Element-level (T-matrices) | ⚠️ Requires global constraint assembly |
| **Parallelization** | ✅ Better (local elimination) | ⚠️ More complex (global coupling) |

**Code Evidence:**
- Extensive T-matrix infrastructure (`cl_FEM_Tmatrix.hpp/cpp`) with CSR storage
- Hanging DOF mechanisms (`is_hanging()`, `has_t_matrix()`)
- Assembly using `T^T * K * T` projection (element-level condensation **before** global assembly)
- No saddle-point matrix structure

**Theoretical Foundation:**
- **Felippa §10.2** — Static condensation theory
- **Bathe §8.5** — Implementation in production FEM codes
- **Messe et al. 2023, Section 3** — BELFEM's application to h-φ interfaces

**Historical Note:** Early versions of Messe et al. 2023 mentioned Lagrange multipliers, but the production code has evolved to use static condensation exclusively, which is the **"preferred method"** cited in Messe et al. 2023.

---

### 19.7 Matrix Function Selection

**Reference:** `IWG_Maxwell::link_to_group`, `cl_IWG_Maxwell.cpp` (anchor: the
`switch ( aGroup->domain_type() )`)

Matrix assembly functions are selected by domain type. For the conductor and
thin-shell domains, the selection is by solver algorithm only; see §20.3 for why:

```cpp
void IWG_Maxwell::link_to_group(Group * aGroup) {
    switch (aGroup->domain_type()) {
        case DomainType::Air:
            mFunMKF = &maxwell::phi_tri3;     // or phi_tet4, phi_tri6_tet10
            break;

        case DomainType::Ferro:
            mFunMKF = (this->algorithm() == SolverAlgorithm::NewtonRaphson) ?
                  &maxwell::phi_ferro_newton
                : &maxwell::phi_ferro_picard;  // Nonlinear μ(H,T)
            break;

        case DomainType::Conductor:
        case DomainType::ThinShell:
            // one case for both: the material ( metal vs. HTS, bulk vs.
            // shell ) is resolved inside MaxwellData, not here. Only the
            // mu tangent needs a kernel distinction
            if (this->algorithm() == SolverAlgorithm::NewtonRaphson) {
                mFunMKF = tMat->is_constant(MaterialProperty::mu) ?
                      &maxwell::h_newton_mu0
                    : &maxwell::h_newton_mu;
            } else {
                mFunMKF = &maxwell::h_picard;
            }
            break;

        case DomainType::InterfaceFerroAir:
            // no ferro-air weak form exists in 2d or 3d — the interface nodes
            // are duplicated for visualization only, and the sideset is set
            // Inactive under h-phi before it can be linked
            BELFEM_ERROR( false, "Ferro-Air Interfaces must be disabled!" );
            break;
    }
}
```

**Disabled Interfaces** (currently throw errors):
- `InterfaceCondAir`: Conductor-Air direct interface
- `InterfaceCondFerro`: Conductor-Ferro direct interface
- `InterfaceFerroAir`: Ferro-Air interface — nodes are duplicated for
  visualization only; no weak form exists in 2d or 3d
- `InterfaceAirCoil`: Air-Coil interface
- `InterfaceFerroCoil`: Ferro-Coil interface

**Reason:** Thin shells now handle conductor-to-air/ferro transitions more accurately.

---

## 20. Thin Shell Formulation

### 20.1 Physics and Motivation

**Problem:** Thin conducting regions (tapes, films, shells) where:
- Thickness `δ << λ` (characteristic length scale)
- Normal field variation is negligible
- In-plane currents dominate

**Solution:** 2.5D formulation where:
- **Tangential fields** computed from edge DOFs (∇ × H = J)
- **Normal fields** imposed from the volume elements on either side (−∇φ on an air/ferro side, the conductor's own Nédélec trace on an h-conductor side, §20.4)

**Benefits:**
- Avoids meshing thin dimension (would require prohibitively fine mesh)
- Captures physics accurately (validated against 3D in Alves et al. 2022b, Alves et al. 2024)
- Computationally efficient (reduces problem dimensionality)

---

### 20.2 Mathematical Decomposition

The magnetic flux density is decomposed:

```
B = B_tangent + B_normal
```

**Where:**
- **B_tangent** = f(edge_h) : In-plane component from curl of H
- **B_normal** = f(φ_master, φ_slave) : Through-thickness component from external φ

**Key Insight:** Normal field penetration is an **external boundary condition**, while tangential fields are determined by **local current flow**.

**Thin-Shell Interface Conditions** (Alves et al. 2022b, Appendix A, Eq. 19):
```
[n × h] = j_s  (tangential jump)
[n · b] = 0    (normal continuity)
```

**Magnetodynamic Coupling** (Arsenault et al. 2023, Eq. 5-7; with the 2026 erratum for the air domain):
Interface coupling derived from Faraday's law ensures tangential **E** continuity when `E^CC = E^C`.

---

### 20.3 Where the Thin-Shell Assembly Lives

There is no thin-shell kernel to call. A thin shell is assembled by the same
H-formulation kernels as a bulk conductor, and everything that makes it a shell
is resolved one level down, in the Calculator's material data object. The split
is worth internalizing before reading any of the code:

| Layer | Anchor | Decides |
|-------|--------|---------|
| Kernel selection | `DomainType::ThinShell` in `IWG_Maxwell::link_to_group`, `cl_IWG_Maxwell.cpp` | solver algorithm (Picard vs. Newton) and whether mu is constant |
| Per-point field and material math | `tIsThinShell` in the `calculator::MaxwellData` constructor, `cl_FEM_Calculator.cpp` | edge vs. node field assembly, the normal-field term, the rho family, the temperature source |

**Kernel selection.** `Conductor` and `ThinShell` share one `case` label, and the
choice inside it never looks at the material's shape:

```cpp
case DomainType::Conductor :
case DomainType::ThinShell :
{
    Material * tMat = aGroup->material();

    // material specifics ( rho law, T source, mu source ) are
    // resolved inside MaxwellData; only the mu tangent needs a
    // kernel distinction: constant mu skips the dMdx blocks
    if ( this->algorithm() == SolverAlgorithm::NewtonRaphson )
    {
        mFunMKF = tMat->is_constant( MaterialProperty::mu ) ?
              & maxwell::h_newton_mu0
            : & maxwell::h_newton_mu ;
    }
    else
    {
        mFunMKF = & maxwell::h_picard ;
    }
    break ;
}
```

`mt_maxwell_h.hpp` exports exactly six kernels: the three above, the facet
stabilization kernel `h_ghost`, and the edge-coating wall pair
`h_side_connector` / `h_side_connector_newton`. The per-material `h_ts_*` family
that earlier revisions of this guide tabulated no longer exists.

**What a kernel actually looks like.** `h_picard` is the short one, and it shows
the pattern: ask `aCalc->maxwell()` for the point values, assemble, and never
branch on the material:

```cpp
void
h_picard( Calculator * aCalc, TimestepMatrices * aMatrices )
{
    const Vector< real > & w = aCalc->integration()->weights();
    calculator::MaxwellData * mx = aCalc->maxwell();
    ...
    for ( uint k = 0; k < aCalc->num_intpoints(); ++k )
    {
        const Matrix< real > & E = aCalc->E( k );
        const Matrix< real > & C = aCalc->C( k );

        real mu   = mx->compute_mu( k );
        real rho  = mx->compute_rho( k );

        real wdV = w( k ) * aCalc->dV( k );

        aMatrices->M() += trans( E ) * E * ( mu * wdV );
        aMatrices->K() += trans( C ) * C * ( rho * wdV );
        ...
    }
}
```

`compute_mu` and `compute_rho` are function-pointer calls into `MaxwellData`.
On a thin-shell block those pointers were bound, at construction, to
implementations that know about the normal field and about the tape normal.
That binding is the whole thin-shell formulation.

**The bindings that matter.** From the `tIsThinShell` branch of the
`MaxwellData` constructor:

```cpp
mFunH = & MaxwellData::compute_h_ts_edge ;   // h = hn + ht, not just E * q
mFunB = & MaxwellData::compute_b_ts ;        // keeps bt and bn separately
```

plus a rho family chosen from the material (§20.6), and a temperature source
chosen from whether a thermal kernel exists (§20.7).

There is no edge-versus-node test in that branch, and there does not need to be:
a `ThinShell` block shares the `Conductor` dof table, which is `edge_h`
(`FieldList::collect_block_dofs`), so a shell is always edge-interpolated. Beware
the near-miss: the `FieldList::ThinShell` member in the same file is the
**sideset** table and does carry `phi`; block routing is the one that decides
interpolation. Ferro and Buffer blocks take `compute_h_bulk_node` instead and do
not run thin-shell facet recovery at all (Air never reaches `MaxwellData`).

---

### 20.4 The Normal Field

The tangential part comes from the block's own edge DOFs. The normal part does
not: it is obtained from the volume elements on either side of the facet through
the free function `compute_hn` in `cl_FEM_Calculator.hpp`, which delegates the
per-side trace to `compute_h_trace` in the same header.

```cpp
inline const Vector< real > &
compute_hn( Calculator * aCalc , const uint k )
{
    Vector< real > & hn = aCalc->vector("hn");

    BELFEM_ASSERT( aCalc->element_is_linear(), "compute_hn doesn't work for ..." );

    if ( k == 0 )
    {
        Vector< real > & phi_m = aCalc->vector("phi_m");   // master dofs ( phi-side only )
        Vector< real > & phi_s = aCalc->vector("phi_s");   // slave dofs  ( phi-side only )

        bool tMasterIsConductor ;
        bool tSlaveIsConductor ;
        Calculator * tCalc = aCalc->get_normal_calculator(
                phi_m, phi_s, tMasterIsConductor, tSlaveIsConductor );

        Vector< real > & hk = aCalc->vector("hk");         // per-point scratch

        Vector< real > & hm = aCalc->vector("hm");
        compute_h_trace( tCalc, true,  tMasterIsConductor, phi_m, hk, hm );  // h on master side

        Vector< real > & hs = aCalc->vector("hs");
        compute_h_trace( tCalc, false, tSlaveIsConductor,  phi_s, hk, hs );  // h on slave side

        hn = 0.5 * ( hm + hs ) ;                           // average

        Vector< real > & n = aCalc->vector("normal");
        n = tCalc->normal( k );

        hn = dot( hn, n ) * n ;                            // PROJECT onto n
    }
    return hn;
}
```

`compute_h_trace` dispatches on the kind of the volume block on that side:

| side | test ( `Calculator::volume_is_conductor` ) | trace |
|---|---|---|
| φ-region ( air, ferro, buffer ) | block `DomainType` is not `Conductor` | `-B_σ(0) * phi_σ`, constant on a linear element |
| h-conductor | block `DomainType::Conductor` | integration-weighted mean over the facet rule of `E_σ(k) * q_σ`, with `q_σ` read from the `edge_h` mesh field ( `nedelec_data_master_h` / `nedelec_data_slave_h` ) |

The conductor side has no potential: its nodal `phi` is bookkeeping ( zero on
conductor-only nodes, whatever the neighboring air holds on shared nodes ) and
must not be differentiated. Its own lowest-order Nédélec field is `a + b × r`, so
the normal trace is linear over the facet and the weighted mean is the exact
centroid value in 2D and 3D. The rule rests on `[n · B] = 0` across the sheet with
`μ = μ₀` on the conductor side, which is why `MaxwellFactory::assign_materials`
refuses a magnetic conductor as a volume neighbor of a shell.

Five properties of this routine are load-bearing:

**1. The result is projected, not just averaged.** The last line,
`hn = dot( hn, n ) * n`, is what makes `hn` a normal field. Dropping it leaves
an `O(|H_t|)` tangential contamination from the volume side, which is the same
order as the tape's own field through `[n x h] = j_s`. Any code that recomputes
this average for display, rather than calling `compute_hn`, has to repeat the
projection or it will report a different `|B|` and a different angle than the
solve used.

**2. It is computed once per element, at `k == 0`.** For flat linear shells the
φ-side gradient is constant over the facet and the h-side mean is taken over the
facet rule inside the `k == 0` branch, so the guard is a cache, not an
approximation. The `BELFEM_ASSERT` on `element_is_linear()` is the guardrail: for
higher-order shells the volume integration points are not the surface points this
routine assumes, and the caching is invalid as well.

**3. `get_normal_calculator` is the sanctioned bridge** from a shell block to
the volume elements on either side. It links the tape sideset's calculator to the
facet, reports per side whether the volume is an h-conductor, fills `phi_m` and
`phi_s` on the φ-sides, and returns a calculator that can produce `Bm`, `Bs`,
`Em`, `Es`, and the facet `normal`. Do not reach into the neighbor blocks by hand.

**4. The mesh field is the live DOF storage.** `Calculator::q()` reads
`mMesh->field( dof->field_index() )`, rank 0 writes every DOF into its field after
each linear solve, and `DofManager::solve_from_residual` distributes the fields to
the workers before the next assembly; worker-rank `Dof::value()` is never updated.
Reading `edge_h` by edge index is therefore identical to the master's `q()` and is
the only accessor that exists for an aura master, which carries no DOF objects.

**5. The normal direction has a sign, and it is a modeling input.** `n` points
from the master volume element to the slave one, which is also the layer-stack
direction. A minus prefix on the `thinshell` sidesets key flips both together.
§20.5 depends on this. The projection itself is sign-blind.

`compute_b_ts` then keeps the two parts apart, because the postprocessor and the
anisotropy law both want them:

```cpp
inline const Vector< real > &
calculator::MaxwellData::compute_b_ts( const uint aIndex )
{
    this->compute_h( aIndex );
    real mu = this->compute_mu( aIndex );
    mBt = mu * mHt ;
    mBn = mu * mHn ;
    mB = mBt + mBn ;
    return mB ;
}
```

The material laws are evaluated on the **total** `b`. That matters: it is what
couples the surrounding field into the tape's critical current and into a
metal's magnetoresistance.

---

### 20.5 Superconductor Anisotropy

An HTS tape's critical current depends on where the field points relative to the
tape plane. BELFEM measures that with `Calculator::bn_angle`, the angle between
**b** and the tape normal **n**:

```cpp
inline real
Calculator::bn_angle( const Vector< real > & b, const Vector< real > & n, real & norm_b ) const
{
    norm_b = norm( b ) ;
    BELFEM_ASSERT( std::abs( norm( n ) - 1.0 ) < 1e-6, "tape normal must be normalized" ) ;

    // UNFOLDED angle between field and tape normal, [ 0, pi ]
    return norm_b < 1e-6 ? constant::pi*0.5 :
        std::acos( std::clamp( dot( n, b ) / norm_b, -1.0, 1.0 ) ) ;
}
```

**The angle is unfolded onto `[0, pi]`, deliberately.** `theta < pi/2` means
**b** has a component along `+n`. Measured `jc(theta)` tables are asymmetric
about `pi/2`, so folding the angle into `[0, pi/2]` would average two physically
different lobes into one. Analytic laws that are even in theta (ModifiedKim)
fold internally, by construction, so they are unaffected. Do not re-fold this
angle; the sideset sign (§20.4, item 4) is what selects the lobe.

**Physical reading:**
- `theta = 0` or `pi`: **B** perpendicular to the tape, giving the minimum Jc for
  a typical REBCO tape (easy vortex penetration through the ab-planes)
- `theta = pi/2`: **B** in the tape plane, giving the maximum Jc (intrinsic
  pinning between the CuO layers)

The lobes at 0 and pi are not required to be equal, and for a real tape they are
not.

**`bn_angle` is not `bj_angle`.** The second angle, between **b** and **j**, is
the Kohler angle a normal metal's magnetoresistance needs, and it is what the
Newton tangent block `add_rho_field_tangent` in `mt_maxwell_h.cpp`
differentiates. The two are bound separately in `MaxwellData` for exactly that
reason: for an HTS block `mFundRhodBeta` stays `return_zero`, because binding it
would apply the metal law's `dbeta/dq` rows to a tape angle they do not
describe. In 2-D, `bj_angle` is `pi/2` by convention (the current is
perpendicular to the model plane), while `bn_angle` is a genuine geometric
quantity in both 2-D and 3-D.

---

### 20.6 Material Models

The rho law is selected in the **deck**, not in C++, and it is bound once per
block in the `MaxwellData` constructor. Nothing downstream branches on it.

**Pure metals.** `compute_rho_metal` is evaluated as `rho(T, |B|, beta_bj)`.
Field dependence follows the Kohler magnetoresistance, and this case binds
the `|B|` and `beta` Newton channels (`compute_drhodb_metal`,
`compute_drhodbeta_metal`).

**HTS.** Three laws, chosen with the `resistivity type` key, each in a plain and
a defect flavor:

| `resistivity type` | Thin-shell binding | Notes |
|---|---|---|
| `power-law` (default), `powerlaw` | `compute_rho_powerlaw_ts` | `E = ec (J/Jc)^n`, in parallel with the normal-state channel (Duron 2004) |
| `piecewise` | `compute_rho_piecewise_ts` | Bezier flux-flow blend above roughly 1.4 Jc |
| `riva` | `compute_rho_riva_ts` | same parallel model as `powerlaw`, made total over the whole table: dead defects and spline over/underflow fall back to the fully normal branch, and `n -> 1` near Tc takes a finite ohmic form instead of producing NaN |

A `defect { file : ...; label : ...; }` subsection swaps each of these for its
`_ts_defect` twin, which additionally passes the integration point's coordinates
and the current time, allowing a plugin to modulate Jc in space and time. The bulk
conductor equivalents carry `_bulk` instead of `_ts` in the same naming scheme.

`jc` and `n` come either from an explicit `jc` + `n` pair or from a
`file` holding `Jc(B,theta,T)` and `n(B,theta,T)` tables; the two are
alternatives, not to be combined. `ec` defaults to 1e-4 V/m. The
`resistivity type` value is case sensitive.

For the table contract, the interpolation scheme, and what the tables are
allowed to do at their edges, see `src/physics/materials/doc/`.

---

### 20.7 Thermal Coupling

There is no thermally coupled kernel either. The temperature source is a third
function pointer, bound according to whether the run has a thermal kernel at all:

```cpp
if ( mThermalCalculator != nullptr )
{
    mFunT = & MaxwellData::compute_T_fem ;
}
else
{
    mFunT = & MaxwellData::compute_T_const ;
}
```

`mThermalCalculator` is resolved in the `MaxwellData` constructor by looking up
this block's id in the thermal kernel's DofManager. The `block_exists()` guard
there is not decoration: `block()` returns the EmptyBlock, which carries a live
calculator, for an unknown id, so a block that is not part of the thermal
problem (a side connector, for instance) would silently adopt the wrong peer.

`compute_T_fem` interpolates the thermal block's own solution at the
integration point:

```cpp
real T = dot( mThermalCalculator->Nvec( aIndex ), mThermalCalculator->q() ) ;

mTClamped = ( T < gTmin ) || ( T > mTmax ) ;
mT = mTClamped ? std::clamp( T, gTmin, mTmax ) : T ;
```

**The clamp is part of the contract, and so is the flag.** Transient nonlinear
iterates routinely swing outside a material table's window, and clamping keeps
the evaluation legal. While clamped, every `dT` derivative is forced to zero, so
the tangent stays consistent with the value actually used. A **converged**
solution sitting at a clamp is a modeling error rather than a numerical one, and
the controller can ask for `T_clamped()` to detect it. `compute_T_const` honors
the same contract for the uncoupled case, pinning a user-set `gTbulk` that falls
outside the window.

Coupling direction: the electromagnetic assembly reads temperature, and the
Joule heating that closes the loop is assembled on the thermal side
(`src/fem/thermal/matrices/mt_thermal_h.cpp`).

---

## 21. Maxwell-Specific Implementation

### 21.1 Custom Vectors for Thin Shells

Thin shell functions request work vectors by name:

**Reference:** `IWG_Maxwell::create_custom_vectors_and_matrices()`

```cpp
Vector<real> & bt = aCalc->vector("bt");      // In-plane field
Vector<real> & bn = aCalc->vector("bn");      // Normal field
Vector<real> & b  = aCalc->vector("b");       // Total field
Vector<real> & j  = aCalc->vector("j");       // Current density
Vector<real> & phi_m = aCalc->vector("phi_m");  // Master phi DOFs ( phi-side )
Vector<real> & phi_s = aCalc->vector("phi_s");  // Slave phi DOFs ( phi-side )
Vector<real> & hk    = aCalc->vector("hk");     // per-point scratch of compute_h_trace
```

On a tape-sideset calculator, `"nedelec_h"` ( master-sized ) and `"nedelec_h_s"`
( slave-sized ) hold the edge DOFs of an h-conductor volume side, filled by
`nedelec_data_master_h()` / `nedelec_data_slave_h()` from the `edge_h` field.

These are **allocated once** during IWG initialization, then **reused** during assembly (avoids repeated allocation in hot path).

---

### 21.2 Calculator Interface for Electromagnetics

**Key Methods:**

```cpp
class Calculator {
    // Nédélec Edge Basis Functions
    const Matrix<real> & E(uint k);     // Edge basis (H(curl) space)
    const Matrix<real> & C(uint k);     // Curl operator (∇ ×)

    // Lagrange Node Basis Functions
    const Matrix<real> & N(uint k);     // Node basis (H1 space)
    Vector<real> Nvec(uint k);          // Node basis as vector
    const Matrix<real> & Bm(uint k);    // Gradient on master side
    const Matrix<real> & Bs(uint k);    // Gradient on slave side

    // Integration
    const Vector<real> & integration()->weights();
    uint num_intpoints();
    real dV(uint k);                    // Jacobian at integration point k

    // DOF values
    const Vector<real> & q();           // Current DOF values (edge_h or phi)

    // Geometric
    const Vector<real> & normal();      // Outward normal vector
    const Matrix<real> & X();           // Node coordinates

    // Work vectors (reusable storage)
    Vector<real> & vector(const string & aName);

    // Normal field calculator (for thin shells): links the tape-sideset
    // calculator to the facet, gathers phi on the phi-sides, reports
    // per side whether the volume is an h-conductor
    Calculator * get_normal_calculator(Vector<real> & phi_m, Vector<real> & phi_s,
                                       bool & aMasterIsConductor, bool & aSlaveIsConductor);

    // edge dofs of the facet's master / slave volume, from the edge_h field
    const Vector<real> & nedelec_data_master_h();
    const Vector<real> & nedelec_data_slave_h();
    bool volume_is_conductor(const mesh::Element * aVolume) const;
};
```

---

### 21.3 Mesh Entity Access for Interfaces

Elements can query their master/slave neighbors and facet information:

**Reference:** `cl_FEM_Element.hpp:244-264`

```cpp
class Element {
    Element * master();           // Master element (for sidesets)
    Element * slave();            // Slave element (for sidesets)
    mesh::Facet * facet();        // Associated mesh facet

    bool has_t_matrix() const;    // True if has hanging DOFs
    const Tmatrix * t_matrix();   // T-matrix for static condensation
};
```

**Usage in Thin Shell Functions:** do not gather neighbor DOFs by hand — the
per-side traces go through `get_normal_calculator` and `compute_h_trace`
( §20.4 ). In particular `Dof::value()` is not the live value off rank 0, and an
aura neighbor has no DOF objects at all; the mesh fields are.

---

### 21.4 Matrix Function File Organization

The matrix kernels live in `src/fem/maxwell/matrices/` (they moved there in
the maxwell kernel split; the layout below is the current one, not the
pre-split single-file layout this guide originally described):

| File | Purpose | Functions |
|------|---------|-----------|
| `mt_maxwell_phi.hpp/cpp` | Phi-formulation matrices | see header |
| `mt_maxwell_h.hpp/cpp` | H-formulation matrices | `h_picard`, `h_newton_mu0`, `h_newton_mu`, `h_ghost`, `h_side_connector`, `h_side_connector_newton` |
| `mt_maxwell_background.hpp/cpp`, `mt_maxwell_symmetry.hpp/cpp`, `mt_maxwell_l2_*.hpp/cpp` | Background field, symmetry, L2 projections | see headers |

**Note:** the old `h_ts_*` thin-shell kernel family and the
`mt_maxwell_thinshell.hpp/cpp` stub are gone — thin shells are handled through
the layer-block machinery, and the stub file was removed as this guide once
recommended.

---

### 21.5 Key Line References for Maxwell

#### H-Phi Formulation
- Formulation enum: `en_Maxwell_Formulations.hpp:22-29`
- Domain types: `src/mesh/en_DomainType.hpp`
- Field assignments: `cl_IWG_Maxwell.cpp:51-62`
- Matrix function routing: `cl_IWG_Maxwell.cpp:165-350`

#### Thin Shells
- Field list construction: `cl_Maxwell_FieldList.cpp:116-156`
- H-formulation kernels (post-split): `src/fem/maxwell/matrices/mt_maxwell_h.cpp` — `h_picard`, `h_newton_mu0`, `h_newton_mu`, `h_ghost`
- HTS anisotropy tangent: `mt_maxwell_h.cpp`, the `rho(|B|,beta)` field-derivative block at the top of the file
- Thermal coupling: `src/fem/thermal/matrices/mt_thermal_h.cpp` (the `M += Nᵀ ρ cp N dV` consumer)

#### Maxwell Factory
- Hanging edge creation: `cl_MaxwellFactory.cpp:1008+`
- Cohomology cuts and condensation: See `hanging_dofs_static_condensation.md`
- Thin shell design: See `../../maxwell/doc/README.md` and `../../maxwell/doc/maxwell_usage_guide.md`

---

### 21.6 Literature References for H-Phi Formulation

**Essential Reading Order:**
1. **messe2023.txt** → Main BELFEM document (Equations 6–8, Section 2.7)
2. **alves2022b.txt** → Thin-shell theory (Appendix A, Eq. 19)
3. **arsenault2023.txt** → Magnetodynamic coupling (Section II, Eq. 5-7)
4. **alves2024.txt** → Thin cuts implementation
5. **dular2021.txt** → Stability theory (inf-sup, hierarchical enrichment)

**Quick Reference Table:**

| Topic | Search Terms | Primary Papers |
|-------|--------------|----------------|
| **H-φ Interface Coupling** | "Faraday", "phi", "interface", "static condensation" | Messe et al. 2023, Arsenault et al. 2023, Alves et al. 2024 |
| **Thin-Shell 2D** | "thin-shell", "auxiliary 1-D", "interface conditions" | Alves et al. 2022b, Messe et al. 2023, Alves et al. 2024 |
| **Thin-Shell 3D** | "3-D thin-shell", "surface representation" | Alves et al. 2022a, Schnaubelt et al. 2023 |
| **Transport Current/Cuts** | "thick cuts", "cohomology", "up(phi) - down(phi)" | Schnaubelt et al. 2023, Alves et al. 2024, Riva et al. 2023 |
| **Stability** | "inf-sup", "hierarchical", "oscillations" | Dular et al. 2021, Messe et al. 2023 |
| **Solver Performance** | "MUMPS", "STRUMPACK", "computation time" | Messe et al. 2023, Arsenault et al. 2021 |

**Paper Locations:** `/home/christian/codes/belfem/literature/papers/fem` as `dular2021.txt` through `schnaubelt2023.txt`

**FEM Textbook References** (for foundational theory):

| Topic | Primary Reference | BELFEM Application |
|-------|-------------------|-------------------|
| **Mixed Formulations** | Brenner Ch. 8 (inf-sup/LBB), Bathe §4.4 | H-Phi coupling stability |
| **Edge Elements (H(curl))** | Zienkiewicz Vol. 1, Hughes, Bathe | Nédélec basis for edge_h DOFs |
| **Static Condensation** | Felippa §10.2, Bathe §8.5 | T-matrix implementation |
| **Newton-Raphson Methods** | Bathe §8.6, Zienkiewicz Vol. 2 Ch. 3 | Nonlinear HTS solver |

**See:** [Section 18: Literature References](#18-literature-references) for complete list.

---

### 21.7 Recommended Study Workflow for Maxwell

**For understanding BELFEM's H-Phi implementation:**

1. **Start with general DOF management**:
   - Read Sections 1-18 of this guide
   - Understand JEDI system, hanging nodes, T-matrices

2. **BELFEM-specific formulation**:
   - Read **Messe et al. 2023** (main BELFEM document - Equations 6-8)
   - Read **Arsenault et al. 2023** (magnetodynamic coupling - Section II)
   - Read **Alves et al. 2022b** (thin-shell theory - Appendix A)
   - Cross-reference with Bathe/Hughes for general FEM context

3. **Deep dive into thin shells**:
   - **Theory**: Alves et al. 2022b (equations) → Code: `mt_maxwell_h.cpp`
     (`h_picard`, `h_newton_mu0`, `h_newton_mu`) plus the thin-shell branch of
     the `calculator::MaxwellData` constructor in `cl_FEM_Calculator.cpp`
   - **Material models**: `src/physics/materials/doc/` → Code: `Material::rho()`, `Material::jc()`
   - **Thermal coupling**: Messe et al. 2023 → Code: `MaxwellData::compute_T_fem()`
     and the Joule-heating consumer in `src/fem/thermal/matrices/mt_thermal_h.cpp`

4. **Validation and debugging**:
   - Benchmarks: Messe et al. 2023, Alves et al. 2022b, Alves et al. 2024 (expected results)
   - Common pitfalls: See Section 1 (Common Pitfalls)
   - Convergence issues: Bathe Ch. 8, Messe et al. 2023, Section 2.7

**Cross-referencing strategy:**
- **"Why does BELFEM do X?"** → Start with papers, verify with books
- **"How should I implement Y?"** → Start with books (algorithms), cross-check with papers (BELFEM specifics)
- **"What's the theory behind Z?"** → Books for general theory, papers for HTS application

---

### 21.8 Summary of Maxwell Implementation

The BELFEM framework implements a sophisticated **H-Phi mixed formulation** for electromagnetic simulations:

1. **DOF Manager** provides unified management of different DOF types (node, edge, face, lambda)
2. **Hanging DOFs** enable static condensation for cohomology cuts and constraints via T-matrices
3. **H-Phi coupling** combines edge-based (H) and node-based (φ) formulations in a single system
4. **Thin shells** serve as interface elements, coupling both formulations through field decomposition:
   - Tangential components from edge DOFs (in-plane currents)
   - Normal components from neighboring φ DOFs (external field)

**This design enables:**
- Accurate modeling of thin conducting regions
- Efficient handling of conductor-insulator interfaces
- Support for complex material models (HTS, magnetoresistance, thermal coupling)
- Topologically correct current boundary conditions via cohomology

**Production-ready with full support for:**
- Pure metals, alloys, and high-temperature superconductors
- Thermal coupling with separate thermal kernel
- Defects and spatial variations
- Multiple constitutive laws (power-law, piecewise E-J)

**Key Insight:** Thin shells naturally couple H and φ formulations by recognizing that normal field penetration is an external boundary condition, while tangential fields are determined by local current flow.

---

## End of DofManager Usage Guide

**For additional information:**
- Module README: `README.md` (Quick reference)
- Project Documentation: `../../../../doc/README.md`
- Coding Philosophy: `../../../../doc/coding_philosophy.md`
- CLAUDE.md: `../../../../CLAUDE.md` (AI assistant instructions)

**Report bugs to:** christian.messe@lbl.gov

**Last updated:** 2026-08-28
