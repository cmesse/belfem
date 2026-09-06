# FEM Kernel Module Documentation {#fem_kernel_index}

**Module:** src/fem/kernel
**Purpose:** Index of documentation for BELFEM's FEM kernel and DOF management system

---

## Overview

The **FEM Kernel** module implements the core finite element orchestration layer in BELFEM. It manages:

- **Degree of Freedom (DOF) allocation and distribution** - Global numbering, reordering, hanging nodes
- **JEDI block matrix system** - Jacobian/Enforcement/Dirichlet/Imposition separation
- **Element-level assembly coordination** - Calculator, weak form integration
- **Sparse solver interface** - UMFPACK, MUMPS, STRUMPACK, PARDISO, PETSc
- **MPI parallel distribution** - Domain decomposition, field synchronization
- **Boundary condition application** - Dirichlet, Neumann, Robin, constraints
- **Postprocessing** - L2 projection of secondary fields

The central class is **DofManager**, which orchestrates all DOF operations through composition of 8 specialized data managers.


## Documentation Files

### User Guides

- **[dof_manager_usage_guide.md](dof_manager_usage_guide.md)** - Comprehensive DofManager usage guide
  - Common pitfalls (historical bug reports are marked FIXED where the source has moved on)
  - Safe defaults quick start
  - Mental model (orchestration pattern, JEDI system)
  - DOF allocation and distribution workflow
  - JEDI block matrix system
  - MPI parallel patterns
  - Ownership and lifetime management
  - Element-level assembly with Calculator
  - Hanging node implementation (T-matrices)
  - Known issues and debugging tips

- **[hanging_dofs_static_condensation.md](hanging_dofs_static_condensation.md)** - Hanging DOFs and static condensation (NEW)
  - **Overview:** Hanging DOF concept and current/future use cases
  - **H-φ interfaces:** Thin-shell conductor to air/ferro (edge-to-node coupling)
  - **H-H interfaces:** Thin-shell conductor to volume conductor (edge-to-edge coupling)
  - **Higher-order elements:** Quadratic edge DOF mapping complexity
  - **Two-phase process:** Mesh-level (geometric) vs DOF-level (per-DOF mapping)
  - **Vertex containers:** Intentional design for visualization (not a bug!)
  - Physics motivation: Why mixed H-φ formulation requires interface coupling
  - Why static condensation over Lagrange multipliers (robustness)
  - T-matrix mathematical formulation and CSR storage
  - Implementation architecture (MaxwellFactory, DofData, DofManager)
  - Detailed code walkthrough with line references
  - Historical bug record (the hardcoded-weight defect is fixed in the source)
  - Debugging techniques and MPI considerations
  - Future applications: Periodic BCs, adaptive mesh refinement
  - Literature references (Messe 2023, Arsenault 2023, Bathe, Monk)

### Theory References

- **[nonlinear_controller_theory.md](nonlinear_controller_theory.md)** - The nonlinear Controller
  - Hybrid Picard→Newton strategy with per-scheme relaxation stores and switch latches
  - Adaptive relaxation rule (Messe et al. 2023, Eq. 14, with the arctan sign note) and backtracking line search
  - Safety nets: divergence rules, escalation, stall guard, progress watchdog, retry hygiene
  - Adaptive timestep control and its configuration warnings
  - Complete input-key reference for the nonlinear sections

- **[bearing_gauge_eigenmode.md](bearing_gauge_eigenmode.md)** - The single-node gauge pin and the near-null φ mode
  - Why a point Dirichlet is a measure-zero pin in 3-D and leaves the Jacobian near-singular
  - The "Newton does nothing" signature (amplified constant-φ Δx, frozen residual) and why Picard is immune
  - Why surface gauging is inadmissible with net transport current (multivalued φ)
  - Remedies: pure Picard (recommended), mean-value gauge row, deflation; diagnostics (`#bearing`; the Δx probe is described in §5)

- **[anderson_acceleration_theory.md](anderson_acceleration_theory.md)** - Anderson-accelerated Picard iteration
  - Type-II Anderson mixing formulation (Walker & Ni) and its depth-1 Aitken limit
  - Column-normalized QR least squares with rank/coefficient safeguards
  - Stage/commit handshake with the controller line search; history flush rules
  - Cost, MPI layout (master-side, no new communication), and the `anderson depth` input key

- **[bdf_timestepping_theory.md](bdf_timestepping_theory.md)** - Variable-step BDF time integration (BDF1–BDF5)
  - Variable-step coefficients from Lagrange differentiation, with fixed-step limits
  - Startup order ramp, step-size history, rejection handling, savepoints
  - Stability notes: A(α)-stability, step-ratio bounds, implicit overestimation of runaway modes

---

## Quick Reference

### Entry Point Classes

| Class | File | Purpose |
|-------|------|---------
| **`Kernel`** | cl_FEM_Kernel.{hpp,cpp} | Top-level orchestrator (owns mesh, materials, BCs, DofManagers) |
| **`DofManager`** | cl_FEM_DofManager.{hpp,cpp} | **Main workhorse**: DOF management, assembly coordination, solver interface |
| **`DofData`** | cl_FEM_DofMgr_DofData.{hpp,cpp} | DOF creation, numbering, hanging nodes |
| **`SolverData`** | cl_FEM_DofMgr_SolverData.{hpp,cpp} | JEDI block matrix system, assembly tables, solver interface |
| **`Calculator`** | cl_FEM_Calculator.{hpp,cpp} | Element-level computation (N, B, J, dV, dS) |

### DofManager Composition (8 Specialized Data Managers)

| Manager | File | Responsibility |
|---------|------|----------------|
| **Parameters** | cl_FEM_DofMgr_Parameters.{hpp,cpp} | Integration orders, scheme settings |
| **DofData** | cl_FEM_DofMgr_DofData.{hpp,cpp} | DOF creation, numbering, hanging nodes |
| **BlockData** | cl_FEM_DofMgr_BlockData.{hpp,cpp} | Block management, thin-shell facet linking |
| **SideSetData** | cl_FEM_DofMgr_SideSetData.{hpp,cpp} | Sideset management, BCs, wetted nodes |
| **BearingData** | cl_FEM_DofMgr_BearingData.{hpp,cpp} | Bearing element management |
| **FieldData** | cl_FEM_DofMgr_FieldData.{hpp,cpp} | Field collect/distribute, linear projection |
| **SolverData** | cl_FEM_DofMgr_SolverData.{hpp,cpp} | Matrix assembly, solver interface, JEDI system |
| **EigenValues** | cl_FEM_DofMgr_EigenValues.{hpp,cpp} | Eigenvalue problem setup |

### Domain Representation Classes

| Class | File | Purpose |
|-------|------|---------|
| **`Group`** | cl_FEM_Group.{hpp,cpp} | Abstract base class for blocks and sidesets |
| **`Block`** | cl_FEM_Block.{hpp,cpp} | Volume elements (material domains) |
| **`SideSet`** | cl_FEM_SideSet.{hpp,cpp} | Surface elements (boundary conditions) |
| **`Bearing`** | cl_FEM_Bearing.{hpp,cpp} | Special constraint elements |
| **`Element`** | cl_FEM_Element.{hpp,cpp} | FEM element wrapper (links DOFs to mesh entities) |

### Core FEM Components

| Class | File | Purpose |
|-------|------|---------|
| **`Dof`** | cl_FEM_Dof.{hpp,cpp} | Degree of freedom (inherits from graph::Vertex) |
| **`Calculator`** | cl_FEM_Calculator.{hpp,cpp} | Element-level engine (N, B, J, dV, dS) |
| **`Tmatrix`** | cl_FEM_Tmatrix.{hpp,cpp} | Hanging node transformations |

---

## Core Enumerations

### MatrixType (JEDI System)

```cpp
enum MatrixType
{
    System,       // A: Free → Free (Picard operator)
    Jacobian,     // A + dA/dx: Free → Free (Newton tangent, same pattern as A)
    Enforcement,  // E: Fixed → Free (enforces constraints, computes forces)
    Dirichlet,    // D: Free → Fixed (Dirichlet boundary conditions)
    Imposition,   // I: Fixed → Fixed (self-coupling of imposed values)
    FullMass,     // Full mass matrix including fixed DOFs
    FullStiffness // Full stiffness matrix including fixed DOFs
};
```

### SolverAlgorithm

```cpp
enum class SolverAlgorithm
{
    Direct,         // Linear solve (no iteration)
    NewtonRaphson,  // Full Newton-Raphson with Jacobian
    Picard,         // Picard iteration (linearized)
    UNDEFINED
};
```

### BoundaryConditionType

```cpp
enum class BoundaryConditionType
{
    Neumann,
    Dirichlet,
    Bearing,
    Gauge,
    Current,              // Maxwell specific
    Voltage,              // Maxwell specific
    CircuitCurrent,       // Maxwell specific
    CircuitVoltage,       // Maxwell specific
    Background,           // Maxwell specific, impose a background field weakly
    BackgroundDirichlet,  // Maxwell specific, impose a background field through phi
    UNDEFINED
};
```

`BoundaryConditionImposing` ( Free, Dirichlet, Neumann, Alpha, Lambda, Weak ) is the per-dof
imposition mode used by `SideSet`; the Robin/convection case is `Alpha` there, not a
`BoundaryConditionType`.

---

## Common Operations

### Creating a DofManager

```cpp
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_DofManager.hpp"

// Create kernel. KernelParameters is constructed from the mesh --
// there is no default constructor.
Mesh* mesh = new Mesh("mesh.exo");
KernelParameters params(mesh);
Kernel kernel(&params);

// Create equation (IWG)
IWG* iwg = kernel.create_equation(
    IwgType::TransientHeatConduction,
    ModelDimensionality::ThreeD
);

// Create field (DofManager)
DofManager* field = kernel.create_field(iwg);

// Set boundary conditions
field->sideset(1)->impose_dirichlet(300.0);  // T = 300K
field->sideset(2)->impose_neumann(1000.0);   // q = 1000 W/m²

// Initialize (allocates DOFs, matrices)
field->initialize();

// Cleanup handled by Kernel destructor
```

---

### Assembly and Solve

```cpp
// Compute Jacobian and RHS
field->compute_jacobian_and_rhs();

// Solve linear system (or Newton iteration)
field->solve();

// Access solution
Cell<Dof*>& dofs = field->dofs();
for (Dof* dof : dofs) {
    real value = dof->value();
    // ...
}
```

---

### JEDI Block Matrix Access

```cpp
// Access JEDI matrices
SpMatrix* A = field->system_matrix(); // A: Free → Free
SpMatrix* D = field->dirichlet();     // D: Free → Fixed
SpMatrix* E = field->enforcement();   // E: Fixed → Free (reaction forces)
SpMatrix* I = field->imposition();    // I: Fixed → Fixed

// Access RHS vector
Vector<real>& f = field->rhs_vector();

// The JEDI system computes "the Force":
// J*x = f - D*y        (solve for unknowns)
// E*x + I*y = g        (compute reaction forces!)
```

---

### MPI Field Synchronization

```cpp
// Collect field data from all procs to master
field->collect_fields({"T", "q"});

// Distribute field data from master to all procs
field->distribute_fields({"T", "q"});

// Round-trip synchronization
field->synchronize_fields({"T", "q"});
```

---

### Postprocessing (L2 Projection)

```cpp
// Create postprocessors for secondary fields
field->initialize_postprocessors();

// Compute secondary fields (e.g., heat flux from temperature gradient)
field->postprocess();

// Access projected fields on mesh
mesh.field("heat_flux")->data();
```

---

## DofManager Workflow

The `DofManager::initialize()` method orchestrates the entire setup:

**Workflow:**
1. **Create DOFs** (DofData)
   - Count DOFs per entity type (nodes, edges, faces, cells, lambdas)
   - Compute offsets for ID calculation
   - Allocate DOF objects
   - Create DOF map for fast lookup

2. **Split Free/Fixed** (DofData)
   - Separate into free and fixed DOFs based on BC flags

3. **Extract Graph** (SolverData)
   - Compute element-DOF connectivity
   - Compute DOF-element connectivity
   - Compute DOF-DOF connectivity (adjacency)

4. **Reorder DOFs** (DofData)
   - Apply RCM or graph partitioning (METIS)
   - Assign new indices to DOFs

5. **Allocate Matrices** (SolverData)
   - Populate adjacency for J, D, E, I matrices
   - Allocate sparse matrix structures
   - Create assembly tables for MPI distribution

6. **Initialize Solver** (SolverData)
   - Create Solver object (UMFPACK, MUMPS, STRUMPACK, etc.)

---

## JEDI Block Matrix System

**The JEDI system literally computes "the Force"** - reaction forces at constraints!

**Block form:**
```
[ J  D ] [ x ]   [ f ]
[ E  I ] [ y ] = [ g ]
```

Where:
- `x` = free DOFs (unknowns to solve for)
- `y` = fixed DOFs (imposed Dirichlet values)

**Physical interpretation:**
- **First equation:** `J*x = f - D*y` → Solve for unknowns
- **Second equation:** `E*x + I*y = g` → Compute reaction forces at constraints

**Why separate matrices?**
1. Enable reaction force computation (E matrix computes constraint forces)
2. Optimize solver performance (J is typically smaller and better conditioned)
3. Support static condensation (D and E couple free and fixed DOFs)

**Matrix dimensions:**
```
J: (N_free × N_free)    - Main system to solve
D: (N_free × N_fixed)   - Free DOF dependency on fixed
E: (N_fixed × N_free)   - Fixed DOF dependency on free (reaction forces)
I: (N_fixed × N_fixed)  - Fixed DOF self-coupling
```

---

## DOF Numbering Strategy

**Goal:** Minimize bandwidth for direct solvers, good locality for iterative solvers

**DOF ID Calculation (contiguous ranges with offsets):**
```cpp
node_dof_id   = nodeID * numDofTypes + dofType
edge_dof_id   = edgeDofOffset + edgeID * numDofTypes + dofType
face_dof_id   = faceDofOffset + faceID * numDofTypes + dofType
cell_dof_id   = cellDofOffset + cellID * numDofTypes + dofType
lambda_dof_id = lambdaDofOffset + facetID * numDofTypes + dofType
```

**Reordering Approach:**
1. Create graph where nodes = DOFs, edges = DOF connectivity
2. Apply reordering algorithm (RCM, METIS nested dissection, etc.)
3. Assign new indices based on reordered graph

**Implementation:** `DofData::reorder_dofs()` in `cl_FEM_DofMgr_DofData.cpp`

---

## Hanging Node Implementation

**What are hanging nodes?**
Hanging nodes occur at mesh refinement boundaries where a node on a coarse element lies on an edge/face of refined elements. The hanging node's DOF is a weighted sum of parent DOFs.

**Data Structure:**
```cpp
class Dof {
    uint mNumberOfSources;      // Number of parent DOFs
    Dof** mSources;             // Parent DOF pointers
    real* mCoefficients;        // Interpolation weights
};

// Example: hanging node h = 0.5 * a + 0.5 * b
h->mNumberOfSources = 2;
h->mSources = [a, b];
h->mCoefficients = [0.5, 0.5];
```

**Workflow:**
1. `DofData::collect_hanging_dofs()` - Identify hanging DOFs
2. `DofData::remove_hanging_dofs_from_container()` - Separate from active DOFs
3. `Element::relink_dofs()` - Replace hanging DOFs with T-matrix transformation
4. Hanging DOFs are NOT part of the system matrix (eliminated via T-matrix)


---

## MPI Distribution Strategy

**Mesh Partitioning (Kernel):**
- Master (rank 0) owns full mesh
- Workers receive submesh with owned + aura elements
- CommTables define communication patterns

**DOF Distribution:**
- Each proc creates DOFs for owned + aura entities
- DOF numbering is global and consistent
- Free/fixed split is consistent across all procs

**Matrix Assembly:**
- Each proc assembles local element contributions
- SolverData maintains assembly tables: `Cell<Vector<index_t>>`
- Tables map global DOF indices to local matrix positions
- `collect_matrices()` gathers contributions to master using tables

**Field Synchronization:**
- `collect_fields()` - Gather field data from all procs to master
- `distribute_fields()` - Scatter field data from master to all procs
- `synchronize_fields()` - Collect then distribute (round-trip)

---

## Calculator Function Pointer Pattern

**Why function pointers instead of virtual functions?**

Performance. Virtual function calls have overhead (~5-10% in tight loops). The Calculator is called millions of times during assembly.

**Pattern:**
```cpp
// Function pointers for different formulations
const Vector<real>& (Calculator::*mFunNormal)(uint);
const Matrix<real>& (Calculator::*mFunN)(uint);
const Matrix<real>& (Calculator::*mFunB)(uint);
real (Calculator::*mFundV)(uint);

// Set function pointer based on problem type
if (dimensionality == 2D) {
    mFunN = &Calculator::N2D;
    mFundV = &Calculator::dV_2D;
} else {
    mFunN = &Calculator::N3D;
    mFundV = &Calculator::dV_3D;
}

// Call through pointer (no virtual dispatch)
Matrix<real>& N = (this->*mFunN)(integration_point);
```

---

## Common Pitfalls

1. **Count accessors are copy-paste-prone** → the historical `number_of_hanging_dofs()` bug is fixed; both count accessors are safe. When adding one, check the returned member matches the name.

2. **Forgetting to call `initialize()`** → DOFs and matrices not allocated. Call before assembly.

3. **Calling `compute_jacobian_and_rhs()` without reset** → Previous values accumulate. Use `aReset=true` (default).

4. **Accessing JEDI matrices before initialization** → Matrices are `nullptr`. Call `initialize()` first.

5. **Modifying DOF values during assembly** → Race conditions in parallel. Modify only in `solve()` or after.

6. **Not synchronizing fields after solve** → MPI procs have stale data. Call `synchronize_fields()`.

7. **Deleting IWG before DofManager** → Dangling pointer. Let Kernel manage lifetime.

8. **Mixing DOF ID and DOF index** → ID is global unique, index is reordered position. Use mapping.

---

## Source Code

**Module location:** `../`

**Key source files:**

- **Orchestration**: `cl_FEM_Kernel.{hpp,cpp}`, `cl_FEM_DofManager.{hpp,cpp}`, `cl_FEM_DofManagerBase.{hpp,cpp}`
- **Data Managers**: `cl_FEM_DofMgr_{Parameters,DofData,BlockData,SideSetData,BearingData,FieldData,SolverData,EigenValues}.{hpp,cpp}`
- **Domain Classes**: `cl_FEM_{Group,Block,SideSet,Bearing,Element}.{hpp,cpp}`
- **Core FEM**: `cl_FEM_{Dof,Calculator,Tmatrix}.{hpp,cpp}`
- **Utilities**: `cl_FEM_{Postprocessor,Controller,KernelParameters}.{hpp,cpp}`, `cl_{Pipette,MeshChecker}.{hpp,cpp}`
- **Enums**: `en_FEM_{BoundaryConditionType,BoundaryConditionImposing,GroupActivationMode,SolverAlgorithm}.hpp`

---

## External References

### Finite Element Theory

- **Zienkiewicz & Taylor**, "The Finite Element Method" Vol. 1-2 (DOF management, assembly, solvers)
- **Hughes**, "The Finite Element Method" (FEM fundamentals, Chapter 1-5)
- **Bathe**, "Finite Element Procedures", Ch. 3 (Formulation), Ch. 8 (Solvers)

### Sparse Solvers

- **Davis**, "Direct Methods for Sparse Linear Systems" (2006) - UMFPACK, sparse matrix storage
- **MUMPS**: Multifrontal Massively Parallel Solver - [http://mumps.enseeiht.fr/](http://mumps.enseeiht.fr/)
- **STRUMPACK**: Structured Matrices Package - [https://portal.nersc.gov/project/sparse/strumpack/](https://portal.nersc.gov/project/sparse/strumpack/)

### Graph Partitioning

- **Karypis & Kumar**, "METIS - Unstructured Graph Partitioning" (1998)
- **George**, "Nested Dissection of a Regular Finite Element Mesh" (1973)
- **Cuthill & McKee**, "Reducing the Bandwidth of Sparse Symmetric Matrices" (1969)

### Related BELFEM Modules

- **IWG** (`src/fem/iwg/`): Physics layer, weak form assembly, time-stepping
- **Interpolation** (`src/fem/interpolation/`): Shape functions, integration points
- **Sparse** (`src/sparse/`): Sparse matrix storage, solver interfaces
- **Graph** (`src/math/graph/`): Graph algorithms (BFS, DFS, RCM, METIS)
- **Mesh** (`src/mesh/`): Mesh data structures, element definitions
- **Materials** (`src/physics/materials/`): Material properties

---

## Development Notes

### Adding New BC Types

When implementing new boundary condition types:

1. **Add to enum** (`en_FEM_BoundaryConditionType.hpp`):
   ```cpp
   enum class BoundaryConditionType {
       ...,
       CustomBC,
       UNDEFINED
   };
   ```

2. **Update SideSet** (`cl_FEM_SideSet.cpp`):
   ```cpp
   case BoundaryConditionType::CustomBC:
       impose_custom_bc(value);
       break;
   ```

3. **Implement assembly** in IWG:
   - Modify `compute_rhs()` to handle new BC type
   - Add BC-specific data storage if needed

4. **Test**:
   - Manufactured solution with known analytical result
   - Mesh convergence study
   - MPI parallel consistency check

### Debugging DOF Issues

**Common debugging commands:**

```cpp
// Print DOF information for rank 0
field->print(0);

// Check DOF counts
index_t n_free = field->solver_data()->number_of_free_dofs();
index_t n_fixed = field->solver_data()->number_of_fixed_dofs();


// Verify DOF connectivity
Dof* dof = field->dof(dof_id);
for (uint i = 0; i < dof->number_of_vertices(); ++i) {
    Dof* neighbor = (Dof*) dof->vertex(i);
    // ... check connectivity
}
```

---

## Performance Profiling

Use `cl_Profiler` and `cl_Timer` for benchmarking:

```cpp
Timer tTimer;

field->initialize();
message(InfoLevel::Verbose, "DOF initialization: %u ms", (uint) tTimer.stop());

tTimer.reset();
field->compute_jacobian_and_rhs();
message(InfoLevel::Verbose, "Assembly: %u ms", (uint) tTimer.stop());

tTimer.reset();
field->solve();
message(InfoLevel::Verbose, "Solve: %u ms", (uint) tTimer.stop());
```

---

## See Also

- **Project README**: `../../../../README.md`
- **Claude Instructions**: `../../../../CLAUDE.md`
- **Coding Philosophy**: `../../../../doc/coding_philosophy.md`
- **Documentation Guidelines**: `../../../../doc/documentation_guidelines.md`
- **General Documentation**: `../../../../doc/README.md`
- **IWG Module Documentation**: `../iwg/doc/README.md`
