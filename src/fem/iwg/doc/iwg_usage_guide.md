# IWG Module Usage Guide {#fem_iwg_iwg_usage_guide}

**Date:** 2026-01-20
**Module:** src/fem/iwg
**Purpose:** Comprehensive guide to BELFEM's Integral Weak Form (IWG) physics module for finite element analysis

---

## Table of Contents

1. [Common Pitfalls](#common-pitfalls)
2. [Safe Defaults Quick Start](#safe-defaults-quick-start)
3. [Glossary](#glossary)
4. [Overview](#overview)
5. [Architecture](#architecture)
6. [Ownership and Lifetime](#ownership-and-lifetime)
7. [Mental Model: The Four Pillars](#mental-model-the-four-pillars)
8. [Usage Examples](#usage-examples)
9. [TimestepMatrices (The Main Workhorse)](#timestepmatrices-the-main-workhorse)
10. [Time-Stepping Algorithms](#time-stepping-algorithms)
11. [Weak Form Assembly](#weak-form-assembly)
12. [Nonlinear Solvers](#nonlinear-solvers)
13. [Factory Patterns](#factory-patterns)
14. [Material Property Integration](#material-property-integration)
15. [Performance Considerations](#performance-considerations)
16. [Thread Safety and MPI](#thread-safety-and-mpi)
17. [Common Patterns](#common-patterns)
18. [Known Limitations](#known-limitations)
19. [Development Notes](#development-notes)

All code examples assume `using namespace belfem;` and `using namespace belfem::fem;`.

---

## Common Pitfalls

**Place pitfalls first** - these are the most common mistakes that lead to crashes, incorrect results, or performance issues.

### 1. Calling `set_omega()` Before the IWG Is Linked

`IWG::set_omega()` aborts with `BELFEM_ERROR( mField != nullptr, "equation is not linked to a
field" )` unless the IWG is already linked to a DofManager. The link is made inside
`Kernel::create_field()`. `set_omega()` is also collective: every rank must call it, and the
value from rank 0 is broadcast.

```cpp
// WRONG - aborts: the IWG has no field yet
tIWG->set_omega( 0.9 );
DofManager * tField = tKernel.create_field( tIWG );

// CORRECT - link first, then set the relaxation factor on every rank
DofManager * tField = tKernel.create_field( tIWG );
tIWG->set_omega( 0.9 );
```

---

### 2. Calling `reset_flags()` on the TimestepMatrices

The `MatrixFlag` bits record which matrices an IWG produces. The IWG sets them **once, in its
constructor**. `TimestepMatrices::initialize()` sizes exactly the flagged matrices, `reset()`
zeroes them, and `assemble_dJdx()` adds only the flagged blocks. `reset_flags()` has no
production caller.

```cpp
// WRONG - clears the flags the constructor set; every later reset()
// and assemble_dJdx() then skips those matrices
tIWG->matrices()->reset_flags();
```

Do not clear the flags between elements. To zero the matrices before an element is computed,
call `reset()`.

---

### 3. Calling Element Routines Outside the DofManager

`compute_jacobian()`, `compute_jacobian_and_rhs()` and `compute_mkf()` read the calculator,
material and dof tables of the group the IWG is currently linked to. `DofManager` makes that
link with `link_to_group( block )` before it loops over a block's elements. In normal use you
never call these routines yourself. You call `DofManager::compute_jacobian()` or
`DofManager::compute_jacobian_and_rhs()`. If you do call one directly, first link the group
yourself, as [Example 3](#example-3-inspecting-one-element-matrix) shows.

---

### 4. Shifting the History at the End of a Step

`IWG_Timestep::shift_fields()` copies the current dof fields into the history fields and pushes
the current `delta_time()` onto the step-size history. `FEM::Controller` calls it at the
**start** of every step, **before** it assigns the new `delta_time()`. A driver loop must follow
the same order:

```cpp
for ( uint tStep = 0; tStep < tNumSteps; ++tStep )
{
    tIWG->shift_fields();              // current -> history, previous step size -> history
    tIWG->delta_time() = tDeltaTime;   // size of the step about to be taken
    tField->compute_jacobian_and_rhs();
    tField->solve();
}
```

If a loop shifts only after the solve, step 1 is assembled against history fields that were
never filled from the initial condition, and the startup ramp runs one step late (see
[Startup order ramp](#bdf-coefficient-formulas)).

---

### 5. Treating `IWG_TransientHeatConduction` as a Working Solver

`IWG_TransientHeatConduction::compute_mkf` (`cl_IWG_TransientHeatConduction.cpp`) does not use
the material:

```cpp
// from cl_IWG_TransientHeatConduction.cpp
M += w( k ) * trans( N ) * 1.0 *  N * aCalc->dV( k );
K += w( k ) * trans( B ) * 111.0 * B * aCalc->dV( k );
```

The mass coefficient is the constant `1.0` and the conductivity is the constant `111.0`. The load
vector `f` is zeroed and nothing is added to it. The constructor sets no `MatrixFlag`, and nothing
calls `TimestepMatrices::initialize()` for this class. `M` and `K` are therefore never sized before
`compute_mkf` accumulates into them. Do not use this class for a production solve. See
[Known Limitations](#known-limitations).

---

### 6. `is_maxwell()` Must List Every Maxwell Variant

`is_maxwell()` controls which factory may build a type. `IwgFactory` and
`Kernel::create_equation` refuse Maxwell types, while `MaxwellFactory` accepts them. It covers
`Maxwell` and `MaxwellThermal` (`en_IWGs.hpp`, `is_maxwell`). The comment on the `UNDEFINED`
line of the enum reminds you to update it when adding a Maxwell IWG.

```cpp
inline bool is_maxwell( const IwgType aType )
{
    return aType == IwgType::Maxwell || aType == IwgType::MaxwellThermal ;
}
```

---

### 7. IwgFactory Covers Three Types, Rejects the Rest Explicitly

```cpp
// cl_IwgFactory.cpp, IwgFactory::create_iwg
switch ( aType )
{
    case IwgType::Poisson:
        return new IWG_Poisson( ... );

    case IwgType::StaticHeatConduction:
        return new IWG_StaticHeatConduction( ... );

    case IwgType::TransientHeatConduction:
        return new IWG_TransientHeatConduction( ... );

    // Gradient2D, Gradient3D, SurfaceGradient, PlaneStress and
    // LinearElasticity are enum values with no IWG class behind them

    default:
        BELFEM_ERROR( false, "invalid type" );  // explicit, not silent
}
```

**Impact:** the factory builds the three heat-conduction and Poisson types. It rejects Maxwell
types up front with "Use MaxwellFactory instead". Every other value hits the explicit error
above.

---

## Safe Defaults Quick Start

**For new users - minimal configuration that follows the lifecycle of `src/fem/kernel/poisson.cpp`:**

```cpp
// 1. Load the mesh on every rank (only the master reads the file)
Mesh tMesh( "model.exo" );

// 2. Create the kernel and let it create the IWG (it uses IwgFactory internally)
KernelParameters tParams( tMesh );
Kernel tKernel( &tParams );                       // collective
IWG * tIWG = tKernel.create_equation(
    IwgType::Poisson,
    ModelDimensionality::ThreeD );                // the kernel owns the IWG

// 3. Select the domain BEFORE create_field
tIWG->select_blocks( { 1 } );
tIWG->select_sidesets( { 1, 2 } );

// 4. Create the DofManager (initializes the IWG and links the two)
DofManager * tField = tKernel.create_field( tIWG );  // collective; the kernel owns it
tField->set_solver( SolverParameters( SolverType::MUMPS ) );

// 5. Boundary conditions, then initialize
tField->sideset( 1 )->impose_dirichlet( 0.0 );
tField->sideset( 2 )->impose_dirichlet( 1.0 );
tField->initialize();  // creates the fields, allocates DOFs and matrices

// 6. Assemble and solve
tField->compute_jacobian();
tField->solve();

// ~Kernel deletes the DofManager and the IWG; tMesh outlives the kernel
```

### What's Currently Safe to Use

| Component | Safe Choice | Why |
|-----------|-------------|-----|
| **Factory-built IWGs** | `Poisson`; `StaticHeatConduction` with a material assigned to every block | `TransientHeatConduction` uses placeholder coefficients (see [Known Limitations](#known-limitations)) |
| **Time-stepping** | BDF1–BDF5, through an IWG whose `compute_mkf` fills its flagged matrices | Variable-step coefficients; Crank–Nicolson and Galerkin are disabled |
| **Nonlinear iteration** | `NewtonRaphson`, with `set_omega()` called after `create_field()` | ω is a fixed relaxation factor |
| **Initialization** | Let `Kernel::create_field()` initialize the IWG | It calls `IWG::initialize()` itself |

### What to Avoid

| Avoid | Why | Alternative |
|-------|-----|-------------|
| `set_omega()` before `create_field()` | `BELFEM_ERROR`: not linked to a field | Call it after `create_field()`, on every rank |
| Factory for types without an IWG class | Runtime error | Use one of the three factory types, or write a new IWG |
| `reset_flags()` | Clears the matrix flags the constructor set | Use `reset()` to zero the matrices |
| `IWG_TransientHeatConduction` for real material data | Hard-coded `1.0` / `111.0` coefficients, unsized `M`/`K` | None in this module |
| `select_blocks()` after `create_field()` | The dof tables are already built | Select blocks and sidesets before `create_field()` |

---

## Glossary

**Common acronyms and terms used in the IWG module:**

| Term | Full Name / Definition |
|------|------------------------|
| **IWG** | Integral Weak form of Governing equations - implements physics PDEs |
| **DOF** | Degree Of Freedom - unknowns in the FEM system (e.g., T at each node) |
| **BDF** | Backward Difference Formula - implicit time integration methods (BDF1-5) |
| **MPI** | Message Passing Interface - parallel computing standard |
| **PDE** | Partial Differential Equation |
| **FEM** | Finite Element Method |
| **Galerkin** | Weighted residual method where test functions = trial functions |
| **Gauss points** | Integration points for numerical quadrature (2-27 per element) |
| **Jacobian** | Matrix of derivatives ∂R/∂x for Newton-Raphson |
| **Residual** | R(x) - error in satisfying the PDE |
| **Stiffness matrix (K)** | ∫ B^T λ B dV - spatial derivative terms |
| **Mass matrix (M)** | ∫ N^T ρcp N dV - time derivative terms |
| **Force vector (f)** | ∫ N^T q dV - source/load terms |
| **Shape functions (N)** | Interpolation functions (e.g., linear, quadratic) |
| **B-matrix** | Derivative of shape functions (e.g., ∇N for diffusion) |
| **Weak form** | Variational formulation of PDE (multiply by test function, integrate by parts) |
| **Assembly** | Accumulating element contributions into global system |
| **Newton-Raphson** | Nonlinear solver with quadratic convergence |
| **Picard** | Fixed-point nonlinear iteration (linear convergence) |
| **A-stable** | Unconditionally stable time integration (no timestep limit) |
| **CSR** | Compressed Sparse Row - sparse matrix storage format |
| **UMFPACK, MUMPS, STRUMPACK** | Direct sparse linear solvers |
| **Voigt notation** | 6-component vector for 3D stress/strain tensors |

---

## Overview

The **IWG (Integral Weak form of Governing equations)** module provides the physics layer for BELFEM’s finite element framework. It provides:

- **PDE formulations:** Poisson (the Laplace operator), static heat conduction, and transient heat conduction in this directory. Maxwell IWGs are in `src/fem/maxwell` and `src/fem/thermal`. The enum includes elasticity types, but no IWG class implements them.
- **Time integration:** BDF1–BDF5, static mode, and forward-Euler mode. Crank–Nicolson and Galerkin are present in the `EulerMethod` enum but **disabled**. Selecting either raises `BELFEM_ERROR` (`cl_IWG_Timestep.cpp`, `set_timestepping_method`) because Newton tangent assembly supports only BDF.
- **Nonlinear solvers:** Newton-Raphson, Picard, and direct linear solves.
- **Element assembly:** Mass (M), damping (D), stiffness (K), and force (f) matrices.
- **Jacobian linearization:** Hand-written contracted derivative blocks—∂M/∂x·x, ∂K/∂x·x, and ∂f/∂x—computed by each IWG kernel.

### Purpose

The IWG module connects the mathematical formulation to the discrete FEM implementation:

1. **Mathematical:** Weak form of the governing PDE
2. **IWG:** Element-level integral evaluation
3. **Assembly:** Construction of the global sparse system
4. **Solver:** Iterative or direct solution

### Key Features

1. **Modular physics:** Each PDE has its own IWG class, such as Poisson and HeatConduction.
2. **Automatic time integration:** BDF methods support variable time steps.
3. **Nonlinear iteration:** Newton-Raphson uses a fixed relaxation factor ω and no line search.
4. **Material coupling:** IWGs integrate directly with the materials module.
5. **Dense element matrices:** The implementation is optimized for typical element sizes of 4–27 DOFs.

### Design Philosophy

Following `doc/coding_philosophy.md`:

- **Performance first:** Dense element matrices and flag-based assembly avoid unnecessary recomputation.
- **Zero abstraction penalty:** Virtual-function overhead is limited to the element level, not the Gauss-point level.
- **Manual memory management:** The Kernel owns the IWGs it creates or adopts. Each Group owns its Calculator, while the IWG only borrows it.
- **Explicit ownership:** `TimestepMatrices` is owned by `IWG_Timestep`.

---

## Architecture

### File Organization

```
src/fem/iwg/
├── en_IWGs.hpp                            # Enumerations (IwgType, ModelDimensionality, etc.)
├── en_IWG_SideSetDofLinkMode.hpp         # SideSet DOF linking modes
├── cl_IWG.{hpp,cpp}                       # Abstract base class
├── cl_IWG_Timestep.{hpp,cpp}             # Time-stepping base class
├── cl_TimestepMatrices.{hpp,cpp}         # Element matrix container
├── cl_IwgFactory.{hpp,cpp}               # Factory for IWG creation
├── cl_IWG_Poisson.{hpp,cpp}              # Poisson equation
├── cl_IWG_StaticHeatConduction.{hpp,cpp} # Steady-state thermal
└── cl_IWG_TransientHeatConduction.{hpp,cpp} # Transient thermal
```

### Class Hierarchy

```
IWG (abstract base)
├── IWG_Timestep (time-stepping base)
│   ├── IWG_Poisson
│   ├── IWG_StaticHeatConduction
│   ├── IWG_TransientHeatConduction
│   │   └── IWG_MaxwellThermal   (src/fem/thermal; overrides compute_mkf)
│   └── IWG_Maxwell              (src/fem/maxwell)
└── IWG_MaxwellPostproc          (src/fem/maxwell)

TimestepMatrices (container for element matrices)
    └── Owned by IWG_Timestep

IwgFactory (factory for non-Maxwell IWGs)
```

### Key Enumerations

**IwgType** (Physics selection):
```cpp
enum class IwgType
{
    Poisson,                    // IWG_Poisson: Laplace operator, no volume source
    StaticHeatConduction,       // IWG_StaticHeatConduction: -∇·(λ∇T) = 0
    TransientHeatConduction,    // IWG_TransientHeatConduction (placeholder coefficients)
    Gradient2D,                 // enum value only, no IWG class
    Gradient3D,                 // enum value only, no IWG class
    SurfaceGradient,            // enum value only, no IWG class
    PlaneStress,                // enum value only, no IWG class
    LinearElasticity,           // enum value only, no IWG class
    Maxwell,                    // IWG_Maxwell (h-φ formulation)
    MaxwellThermal,             // IWG_MaxwellThermal: thermal side of the coupled problem
    UNDEFINED
};
```

**ModelDimensionality** (Spatial dimension):
```cpp
enum class ModelDimensionality
{
    TwoD,       // 2D Cartesian (x, y)
    AxSymmX,    // Axisymmetric about X-axis (r, z with x=0)
    AxSymmY,    // Axisymmetric about Y-axis (r, z with y=0)
    ThreeD,     // 3D Cartesian (x, y, z)
    UNDEFINED
};
```

**IwgMode** (Linear vs. nonlinear):
```cpp
enum class IwgMode
{
    Direct,     // Linear solve (assemble once, solve once)
    Iterative,  // Nonlinear solve (Newton-Raphson or Picard)
    UNDEFINED
};
```

`IWG_Poisson` defaults to `IwgMode::Direct`, and `set_algorithm()` rejects a Direct-mode IWG.

`IWG_StaticHeatConduction` and `IWG_TransientHeatConduction` default to `IwgMode::Iterative`. Their default algorithm is `NewtonRaphson`.

**SolverAlgorithm** (Nonlinear solver type):
```cpp
enum class SolverAlgorithm
{
    Direct,         // Linear solve (no iteration)
    NewtonRaphson,  // Full Newton with Jacobian
    Picard,         // Fixed-point iteration (linearized)
    UNDEFINED
};
```

### Factory Support Matrix

⚠️ **Important:** `IwgFactory` builds only 3 of the 10 concrete `IwgType` values:

| IwgType | Factory Support | Alternative |
|---------|----------------|------------|
| **Poisson** | ✅ Yes | - |
| **StaticHeatConduction** | ✅ Yes | - |
| **TransientHeatConduction** | ✅ Yes (placeholder coefficients) | - |
| Gradient2D | ❌ No | None: no IWG class exists |
| Gradient3D | ❌ No | None: no IWG class exists |
| SurfaceGradient | ❌ No | None: no IWG class exists |
| PlaneStress | ❌ No | None: no IWG class exists |
| LinearElasticity | ❌ No | None: no IWG class exists |
| **Maxwell** | ❌ Explicitly blocked | `MaxwellFactory` |
| **MaxwellThermal** | ❌ Explicitly blocked | `ThermalFactory` (`cl_ThermalFactory.cpp`) |

A type outside the factory needs its own IWG class. Construct it with `new`, then hand it to the
kernel with `Kernel::add_equation()` (see [Ownership Pattern 2](#common-ownership-patterns) and
[Adding New IWG Types](#adding-new-iwg-types)).

---

## Ownership and Lifetime

**Understanding memory ownership prevents leaks and double-frees:**

### Ownership Rules

| Object | Owned By | Lifetime | Deletion |
|--------|----------|----------|----------|
| **IWG** | Kernel (via `create_equation`/`add_equation`) or the user (manual `new`, never handed over); ownership passes to the kernel only at `add_equation()` | Until `~Kernel` | Kernel destructor deletes it; the DofManager only borrows it |
| **TimestepMatrices** | IWG_Timestep | Same as IWG | IWG_Timestep destructor frees |
| **Calculator** (wraps the IntegrationData) | Group (Block or SideSet of a DofManager) | Created by `Group::create_calculator()` | Group frees it; the IWG borrows it in `link_to_group()` |
| **DofManager** | Kernel (`create_field`) | Until `~Kernel` | Kernel destructor deletes it |
| **Material** | Kernel when added with `add_material()` | Until `~Kernel` | Kernel destructor deletes it; Groups and IWGs borrow it |
| **Mesh** | User, unless a factory handed it over (`claim_mesh_ownership`) | Must outlive the kernel | User must delete |

### Common Ownership Patterns

**Pattern 1: Kernel owns IWG and DofManager (typical)**
```cpp
KernelParameters tParams( tMesh );
Kernel tKernel( &tParams );
IWG * tIWG = tKernel.create_equation( ... );           // Kernel owns the IWG
DofManager * tField = tKernel.create_field( tIWG );    // Kernel owns the DofManager

// ... use tField ...

// ~Kernel deletes the DofManagers, then the IWGs
// DO NOT: delete tIWG; delete tField;  ❌ Double-free!
```

**Pattern 2: User-built IWG handed to the Kernel (advanced)**
```cpp
IWG * tIWG = new IWG_CustomPhysics( ModelDimensionality::ThreeD );  // your own class
tIWG->select_blocks( { 1 } );
tKernel.add_equation( tIWG );                         // Kernel now owns it
DofManager * tField = tKernel.create_field( tIWG );

// an IWG never handed to the kernel stays the user's: delete it yourself
```

### What Gets Invalidated When

| Call | What Happens | Notes |
|------|--------------|-------|
| `iwg->select_blocks()` | Replaces the block selection and recounts the dofs per block | Call before `create_field()` |
| `iwg->select_sidesets()` | Replaces the sideset selection | Call before `create_field()`; `DofManager::sideset( id )` exists only for selected sidesets |
| `kernel.create_field( iwg )` | Calls `iwg->initialize()`, then `DofManager::set_equation()`, which calls `iwg->set_field()` | Collective |
| `iwg->initialize()` | Builds the dof-type map and counts the dofs per block and sideset | Called by `create_field()`; `DofManager::set_equation()` refuses an uninitialized IWG |
| `iwg->set_field()` | Stores the DofManager and mesh pointers; `IWG_Timestep` also creates the history fields | Called by `create_field()`, not by user code |
| `iwg->shift_fields()` | x_old ← x (time stepping) | - |
| `iwg->reset_fields()` | x ← x_old (reject timestep) | - |
| `dofManager->initialize()` | Creates the fields and allocates the DOF system | Returns at once on a second call |

### Memory Leaks to Avoid

```cpp
// LEAK 1: Creating IWG in loop without deletion
for ( int i = 0; i < 1000; ++i )
{
    IWG* iwg = factory.create_iwg( ... );  // ❌ Leaks 999 IWGs
    // ... use iwg ...
    // Missing: delete iwg  (or kernel.add_equation( iwg ))
}

// LEAK 2: Factory-built IWG never handed to the kernel
IWG* iwg = factory.create_iwg( ... );
DofManager* dofManager = kernel.create_field( iwg );  // links, but does not take ownership
// Missing: kernel.add_equation( iwg )  ❌ Leaks the IWG

// CORRECT: hand it over, ~Kernel deletes it
kernel.add_equation( iwg );  // ✅
```

### Double-Free to Avoid

```cpp
// DOUBLE-FREE: Deleting what the Kernel owns
IWG* iwg = kernel.create_equation( ... );             // Kernel owns it
DofManager* dofManager = kernel.create_field( iwg );  // Kernel owns it

delete dofManager;  // ❌ CRASH at ~Kernel: double-free
delete iwg;         // ❌ same

// CORRECT: let ~Kernel clean up
```

---

## Mental Model: The Four Pillars

Understanding BELFEM's FEM framework requires grasping four interconnected pillars:

### Pillar 1: DofManager (Global System)

**What it is:** Manages the global sparse linear system Ax = b.

**Responsibilities:**
- Creates DOFs (degrees of freedom) for all mesh nodes/edges/faces
- Assembles element contributions into global sparse matrix
- Applies boundary conditions (Dirichlet, Neumann)
- Calls linear solver (UMFPACK, MUMPS, STRUMPACK)

**Interaction with IWG:**
```cpp
DofManager* dofManager = kernel.create_field( iwg );  // Link IWG for physics (field "T" comes from the IWG)
dofManager->set_solver( SolverParameters( SolverType::MUMPS ) );
dofManager->initialize();  // Create global DOFs
```

---

### Pillar 2: Domain (Blocks and SideSets)

**What it is:** Groups of elements sharing material properties or boundary conditions.

**Types:**
- **Block**: Volume elements (e.g., conductor, air, insulator)
- **SideSet**: Surface elements (e.g., Neumann BC, interface)

**Interaction with IWG:**
```cpp
iwg->select_blocks( { 1, 2, 3 } );       // Material regions
iwg->select_sidesets( { 10, 11 } );      // Boundary conditions
iwg->set_blocks( blockIDs, blockTypes ); // Domain types for activation
```

---

### Pillar 3: IWG (Physics)

**What it is:** Implements the weak form of the governing PDE.

**Responsibilities:**
- Define DOF fields (e.g., "T" for temperature, "phi" for potential)
- Compute element matrices M, K, f
- Declare in its constructor, through the matrix flags, which matrices it produces
- Handle material property evaluation

**Key virtual methods:**
```cpp
virtual void compute_jacobian( Element* aElement, Matrix<real>& aJacobian );
virtual void compute_mkf( Element* aElement );  // For time-stepping
```

---

### Pillar 4: IntegrationData (Numerical Integration)

**What it is:** Precomputed shape functions and integration points (from `interpolation` module).

**Responsibilities:**
- Generate Gauss quadrature points and weights
- Evaluate shape functions N(ξ) at all points (once!)
- Evaluate derivatives ∂N/∂ξ at all points (once!)

**Interaction with IWG:**
```cpp
// Inside IWG::compute_mkf()
Calculator* calc = iwg->calc();  // Calculator wraps IntegrationData
const Matrix<real>& N = calc->N( k );      // Shape at point k
const Matrix<real>& B = calc->B( k );      // Derivatives at point k
real w = calc->integration()->weights()( k );             // Weight
real dV = calc->dV( k );                   // Volume element
```

---

### How the Pillars Work Together

**Workflow for a transient nonlinear solve:**

1. **DofManager** creates global DOFs based on **IWG** field definitions
2. **IWG** selects **Blocks/SideSets** (domains)
3. **DofManager** loops over elements in each block/sideset
4. For each element:
   - **IWG** uses **IntegrationData** (via Calculator) to evaluate weak form
   - **IWG** populates element matrices M, K, f
   - **DofManager** assembles element contributions into global system
5. **DofManager** solves global system
6. Repeat for next timestep

---

## Usage Examples

### Example 1: Static Poisson Equation

`IWG_Poisson` assembles the Laplace operator ∫ ∇N·∇N dV and has no volume source term. With
Dirichlet data on two sidesets, it solves ∇²φ = 0.

```cpp
#include "cl_Mesh.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_KernelParameters.hpp"
#include "cl_FEM_DofManager.hpp"

// Load mesh (every rank; only the master reads the file)
Mesh tMesh( "poisson.exo" );

// Create the kernel and the IWG
KernelParameters tParams( tMesh );
Kernel tKernel( &tParams );
IWG * tIWG = tKernel.create_equation(
    IwgType::Poisson,
    ModelDimensionality::ThreeD );  // IwgMode::Direct: linear problem

// Select domain
tIWG->select_blocks( { 1 } );          // All elements in block 1
tIWG->select_sidesets( { 1, 2 } );     // needed for DofManager::sideset( id )

// Create DOF manager (its field list comes from the IWG)
DofManager * tField = tKernel.create_field( tIWG );  // initializes the IWG and links the two
tField->set_solver( SolverParameters( SolverType::MUMPS ) );

// Dirichlet BCs: phi = 0 on sideset 1, phi = 1 on sideset 2
tField->sideset( 1 )->impose_dirichlet( 0.0 );
tField->sideset( 2 )->impose_dirichlet( 1.0 );

tField->initialize();  // creates the fields, allocates DOFs and matrices

// Solve
tField->compute_jacobian();   // Assemble stiffness matrix
tField->solve();

// Extract solution
const Vector< real > & tPhi = tField->field_data( "phi" );

// Cleanup
// ~Kernel deletes the DofManager and the IWG
```

---

### Example 2: Transient Problems

No IWG in `src/fem/iwg` supports a transient run end to end. `IWG_TransientHeatConduction` uses
placeholder coefficients (see [Known Limitations](#known-limitations)). The production transient
path is the coupled Maxwell–thermal solver: `FEM::Controller` driving `IWG_Maxwell` and
`IWG_MaxwellThermal` (see `src/executables/belfem.cpp`). The skeleton below shows the order of
calls that the controller uses for an `IWG_Timestep`-derived equation whose `compute_mkf` fills
its flagged matrices. See [Time-Stepping Workflow](#time-stepping-workflow).

---

### Example 3: Inspecting One Element Matrix

This example uses the same calls that `DofManager::compute_jacobian()` makes for one block. Run it
after `tField->initialize()` from Example 1.

```cpp
Block * tBlock = tField->block( 1 );

// link the IWG to the block: calculator, material and dof tables
tIWG->link_to_group( tBlock );

const uint tN = tIWG->number_of_dofs_per_element( tBlock );
Matrix< real > tJ( tN, tN );

Element * tElement = tBlock->elements()( 0 );
tIWG->compute_jacobian( tElement, tJ );   // IWG_Poisson: ∫ Bᵀ B dV
```

The next `DofManager` assembly re-links the IWG to each group, so this call leaves no lasting
state. `IWG_StaticHeatConduction` does not implement `compute_mkf`, and
`IWG_Timestep::compute_mkf` aborts with "compute_mkf not implemented for this IWG".

---

### Example 4: Nonlinear Steady-State Heat Conduction

`IWG_StaticHeatConduction` evaluates λ(T) at the current temperature, so every block needs a
material. Its kernel (`compute_conduction`) assembles K(T) and forms no dλ/dT block. The Newton
tangent therefore lacks the conductivity derivative. Do not expect quadratic convergence.

```cpp
#include "cl_MaterialFactory.hpp"

Mesh tMesh( "thermal.exo" );
KernelParameters tParams( tMesh );
Kernel tKernel( &tParams );

IWG * tIWG = tKernel.create_equation(
    IwgType::StaticHeatConduction,
    ModelDimensionality::ThreeD );       // IwgMode::Iterative, NewtonRaphson by default

tIWG->select_blocks( { 1 } );
tIWG->select_sidesets( { 1, 2 } );

DofManager * tField = tKernel.create_field( tIWG );
tIWG->set_omega( 0.9 );                  // after create_field, on every rank
tField->set_solver( SolverParameters( SolverType::MUMPS ) );

// material: the kernel owns it, the block borrows it by label
MaterialFactory tMaterialFactory ;
tKernel.add_material( "copper", tMaterialFactory.create_material( "copper" ) );
tField->block( 1 )->set_material( "copper" );

tField->sideset( 1 )->impose_dirichlet( 300.0 );
tField->sideset( 2 )->impose_dirichlet( 350.0 );
tField->initialize();
tField->field_data( "T" ).fill( 300.0 );  // starting guess

const real tTolerance = 1e-6 ;
for ( uint tIter = 0; tIter < 50; ++tIter )
{
    tField->compute_jacobian_and_rhs();
    tField->solve();                      // applies the update scaled by ω

    // relative residual, broadcast from rank 0
    if ( tField->residual( tIter ) < tTolerance )
    {
        break ;
    }
}

// ~Kernel deletes the DofManager, the IWG and the material
```

---

## TimestepMatrices (The Main Workhorse)

`TimestepMatrices` is the **central container** for element-level matrices in transient and nonlinear FEM. It serves as the interface between the physics-aware IWG and the time-stepping scheme.

### Design Decisions

**Why dense matrices?**
- Element matrices are small (typically 4-27 DOFs for common elements)
- Dense storage is faster than sparse for small matrices
- Allows use of optimized BLAS routines (gemm, gemv)

**Why contracted derivatives?**
- Full derivative tensors ∂K/∂x would be (n_e × n_e × n_e) - too large
- IWG computes contracted form directly: (∂K/∂x) · x → (n_e × n_e)
- Same memory footprint as K itself

**Why bitset flags?**
- The IWG declares once, in its constructor, which matrices it produces
- `initialize()` sizes only those, `reset()` zeroes only those
- `assemble_dJdx()` adds only the flagged derivative blocks

### Matrix Flags

```cpp
enum class MatrixFlag : index_t
{
    M            = 0,   // Mass matrix populated
    D            = 1,   // Damping matrix (reserved for 2nd-order systems)
    K            = 2,   // Stiffness matrix populated
    F            = 3,   // Force vector populated
    dMdX_times_x = 4,   // (∂M/∂x) · x for Newton-Raphson
    dMdX_times_h = 5,   // (∂M/∂x) · h (history correction)
    dKdX_times_x = 6,   // (∂K/∂x) · x for Newton-Raphson
    dFdX         = 7    // ∂f/∂x for Newton-Raphson
};
```

### Workflow

This is the sequence `IWG_MaxwellThermal` follows (`cl_IWG_MaxwellThermal.cpp`):

**Step 1: the IWG constructor sets the flags**
```cpp
mTimeStepMatrices->set_flag( MatrixFlag::M ) ;
mTimeStepMatrices->set_flag( MatrixFlag::K ) ;
mTimeStepMatrices->set_flag( MatrixFlag::F ) ;
```

**Step 2: `link_to_group()` sizes the flagged matrices**
```cpp
uint tNumDofs = mGroup->elements()( 0 )->number_of_local_dofs() ;
mTimeStepMatrices->initialize( tNumDofs ) ;   // flagged matrices plus J and dJdx
```

**Step 3: per element, `compute_mkf()` zeroes and fills them**
```cpp
// Inside IWG::compute_mkf()
mTimeStepMatrices->reset();   // zeroes the flagged matrices, J and dJdx

for ( uint k = 0; k < mCalc->num_intpoints(); ++k )
{
    // M += ..., K += ..., f += ... ( see Weak Form Assembly )
}
```

**Step 4: the time-stepper combines them**

`IWG_Timestep::compute_jacobian_and_rhs()` calls `compute_mkf()`, then
`assemble_dJdx( delta_time, alpha )`, then the active scheme (`bdf1`, `bdf2`, ...), which writes
the element Jacobian and right-hand side (see [Time-Stepping Algorithms](#time-stepping-algorithms)).

**Step 5: global assembly**

`DofManager::compute_jacobian_and_rhs()` passes the element Jacobian, `matrices()->dJdx()` and the
right-hand side to its solver data (`assemble_jacobian`, `assemble_newton`, `assemble_rhs`). User
code does not call these.

### API Reference

**Initialization:**
```cpp
void initialize( const index_t aNumDofs );  // Size the flagged matrices, plus J and dJdx, to n_e
void reset();                                // Zero the flagged matrices, J and dJdx; keep flags
void reset_flags();                          // Clear flag bitset (no production caller; do not use)
```

**Flag management:**
```cpp
void set_flag( const MatrixFlag aFlag );             // Declare a matrix (IWG constructor)
bool has_flag( const MatrixFlag aFlag ) const;       // Is the matrix produced?
const Bitset<8>& flags() const;                      // Direct access to all flags
```

**Matrix accessors (const):**
```cpp
const Matrix<real>& M() const;           // Mass matrix
const Matrix<real>& D() const;           // Damping matrix
const Matrix<real>& K() const;           // Stiffness matrix
const Vector<real>& f() const;           // Force vector
const Matrix<real>& dMdx_times_x() const;  // (∂M/∂x) · x
const Matrix<real>& dMdx_times_h() const;  // (∂M/∂x) · h
const Matrix<real>& dKdx_times_x() const;  // (∂K/∂x) · x
const Matrix<real>& dfdx() const;          // ∂f/∂x
const Matrix<real>& J() const;             // Written only by assemble_J
const Matrix<real>& dJdx() const;          // Newton correction term
```

**Matrix accessors (non-const for IWG):**
```cpp
Matrix<real>& M();    // IWG fills this
Matrix<real>& K();    // IWG fills this
Vector<real>& f();    // IWG fills this
// ... (same for all matrices)
```

**Assembly:**
```cpp
void assemble_J( const real adt );                          // J = M + adt·K (flagged terms only)
void assemble_dJdx( const real adt, const real aAlpha = 1.0 );
```

`assemble_J` has no caller in `src/`. The element Jacobian that the DofManager assembles comes from
the scheme functions in `IWG_Timestep`, not from `mJ`.

`assemble_dJdx` forms, from the flagged blocks only:

```
dJdx = α·(∂M/∂x)·x − (∂M/∂x)·h + Δt·(∂K/∂x)·x − Δt·(∂f/∂x)
```

---

## Time-Stepping Algorithms

The `EulerMethod` enum (`src/sparse/en_SolverEnums.hpp`) has 12 concrete values, two of them
disabled. They target transient PDEs of the form:

**M(x) dx/dt + K(x) x = f(x, t)**

Where:
- M = mass matrix (capacity)
- K = stiffness matrix (diffusion/conductivity)
- f = force vector (source term)

### `EulerMethod` values

All rows below exist in the enum. **Crank–Nicolson and Galerkin are disabled** — selecting either
raises `BELFEM_ERROR` (`cl_IWG_Timestep.cpp`, `set_timestepping_method`), because the Newton
tangent assembly supports the BDF family only. They are listed for completeness, not as options.

Every scheme assembles the element system **multiplied by the step size h**. That is the
unscaled form (for example M/h + K for BDF1) multiplied through by h, so the solution is the
same. In the table, q₀ is the previous step and qhist is the β-weighted history described below.

`IWG_Poisson` overrides `compute_jacobian` and `IWG_StaticHeatConduction` overrides
`compute_jacobian_and_rhs`, so neither runs these scheme functions.

| EulerMethod | Order | A-Stability | Element J | Element RHS | Use Case |
|-------------|-------|-------------|-----------|-------------|----------|
| **Static** | 0 | - | h K | h f | Steady-state |
| **ForwardExplicit** | 1 | ❌ Conditionally | M | h f + (M − h K) q₀ | Explicit (dt limited) |
| **BackwardDifference1** (BDF1) | 1 | ✅ A-stable | M + h K | h f + M q₀ | Implicit Euler |
| **CrankNicolson** | 2 | ✅ A-stable | — | — | ❌ **disabled** — raises `BELFEM_ERROR` |
| **Galerkin** | 1 | ✅ A-stable | — | — | ❌ **disabled** — raises `BELFEM_ERROR` |
| **BackwardDifference2** (BDF2) | 2 | ✅ A-stable | α M + h K | h f + M qhist | Variable step BDF2 |
| **BackwardDifference3** (BDF3) | 3 | A(α)-stable (α ≈ 86°) | α M + h K | h f + M qhist | Higher accuracy |
| **BackwardDifference4** (BDF4) | 4 | A(α)-stable (α ≈ 73°) | α M + h K | h f + M qhist | Higher accuracy |
| **BackwardDifference5** (BDF5) | 5 | A(α)-stable (α ≈ 52°) | α M + h K | h f + M qhist | Highest BDF order |
| **StiffnessOnly** | - | - | K | 0 | Stiffness matrix alone (used by the eigenvalue path) |
| **MassOnly** | - | - | M | 0 | Mass matrix alone (used by the eigenvalue path) |
| **Derivative** | - | - | M | −K q₀ | Solves M ẋ = −K q₀ (no load term) |

### BDF Coefficient Formulas

All BDF methods support **variable time-steps**. Coefficients are computed automatically.

**BDF1 (Backward Euler):**
```
J   = M + h K
RHS = h f + M q0

Where:
  h  = current timestep (delta_time)
  q0 = solution at previous timestep
```

**BDF2 (variable step), `compute_bdf_coeffs_2`:**
```
Given:
  h  = upcoming timestep (delta_time)
  h1 = previous timestep (mH(0))

Coefficients:
  α  = (2h + h1) / (h + h1)
  β0 = (h + h1) / h1
  β1 = h² / (h1 (h + h1))

Assembly:
  J     = α M + h K
  RHS   = h f + M qhist
  qhist = β0 q0 − β1 q1          (for h = h1: α = 3/2, qhist = 2 q0 − ½ q1)
```

**BDF3-5:** The same pattern with more history terms; see `compute_bdf_coeffs_3`,
`compute_bdf_coeffs_4` and `compute_bdf_coeffs_5` in `cl_IWG_Timestep.cpp`. `mBeta` stores
positive magnitudes. `collect_qhist()` applies the alternating signs:
qhist = β0 q0 − β1 q1 + β2 q2 − β3 q3 + β4 q4, truncated to the order that is running.

**Coefficient lifecycle:** `compute_bdf_coefficients()` is invoked lazily from
`IWG_Timestep::compute_jacobian_and_rhs()`. `shift_fields()` and `reset_fields()` mark the
coefficients dirty; the first element assembly of a step recomputes them once, using the
`delta_time()` set at that moment. No manual call is needed.

**Startup order ramp:** a BDF-p run needs p history states. The order that runs is
`ramped_order()` (`cl_IWG_Timestep.hpp`): min(p, n), where n is the number of `shift_fields()`
calls so far, and n = 0 counts as 1. The controller calls `shift_fields()` at the start of every
step, so its step k runs BDF-min(k, p): step 1 runs BDF1, step 2 BDF2, and so on until the
requested order is reached. A driver that calls `shift_fields()` only after the solve runs BDF1 on
both step 1 and step 2. A rejected timestep (`reset_fields()`) rolls the counter back by one.
`order_pending()` reports the order the next assembly will run.

### Newton Jacobian for Nonlinear Time-Stepping

For BDF schemes with nonlinear M(x), K(x) and f(x), `assemble_dJdx` combines the flagged
contracted blocks:
- α·(∂M/∂x)·x: current-iterate derivative, scaled by the BDF coefficient α
- −(∂M/∂x)·qhist: history derivative, unscaled; it must be contracted against the same
  `collect_qhist()` vector the residual uses
- +h·(∂K/∂x)·x: stiffness derivative, scaled by the time-step
- −h·(∂f/∂x): load derivative, scaled by the time-step

Assembled in `TimestepMatrices::assemble_dJdx( adt, aAlpha )`; the α factor is passed
from `IWG_Timestep::compute_jacobian_and_rhs()`. For every scheme that scales M by one
(BDF1 and explicit Euler; also Crank–Nicolson and Galerkin, were they enabled), α = 1. The
DofManager assembles `dJdx()` next to the scheme's element Jacobian (`assemble_newton`).

Reference: Hairer & Wanner 1996, II.4 for the BDF Jacobian structure.

### Selecting a Time-Stepping Method

```cpp
iwg->set_timestepping_method( EulerMethod::BackwardDifference2 );

// For problems with stiffness matrix:
iwg->set_timestepping_method( EulerMethod::BackwardDifference2, true );

// For problems with only mass matrix (e.g., pure diffusion with lumped M):
iwg->set_timestepping_method( EulerMethod::BackwardDifference2, false );
```

**With stiffness (`aHaveStiffness = true`):**
- Assembles J = α M + h K
- Standard for diffusion-reaction, elastodynamics

**Without stiffness (`aHaveStiffness = false`):**
- BDF2–BDF5 assemble J = α M; BDF1 and `ForwardExplicit` both map to `bdf1_nok`, which
  assembles J = M. The right-hand side is h f + M qhist (BDF1: h f + M q₀)
- `Derivative` maps to `derivative_no_stiffness`, which adds nothing
- `Static`, `CrankNicolson` and `Galerkin` raise `BELFEM_ERROR`; `StiffnessOnly` and `MassOnly`
  raise "Invalid euler method"

The method may be set before or after `create_field()`. If the IWG is already linked, a changed
order rebuilds the history fields and the calculators' history tables.

### Variable Time-Step Support

All BDF methods automatically adapt to changing time-steps:

```cpp
real tDeltaTime = 0.01;  // Initial timestep

for ( uint tStep = 0; tStep < tNumSteps; ++tStep )
{
    tIWG->shift_fields();              // history; mH(0) <- size of the step just completed

    // Adaptive time-stepping (your controller decides)
    if ( tConvergedQuickly )
    {
        tDeltaTime *= 1.5;
    }
    else if ( tConvergedSlowly )
    {
        tDeltaTime *= 0.5;
    }

    tIWG->delta_time() = tDeltaTime;   // size of the upcoming step

    // BDF coefficients recomputed lazily from the step-size history
    tField->compute_jacobian_and_rhs();
    tField->solve();
}
```

### Time-Stepping Workflow

**Initialization:**
```cpp
tIWG->set_timestepping_method( EulerMethod::BackwardDifference2 );
DofManager * tField = tKernel.create_field( tIWG );
tField->set_solver( SolverParameters( SolverType::MUMPS ) );
tField->initialize();
tField->field_data( "T" ).fill( tTinitial );  // Set initial condition
```

**Time loop (the order `FEM::Controller` uses):**
```cpp
for ( uint tStep = 0; tStep < tNumSteps; ++tStep )
{
    // 1. Shift history: q1 <- q0, q0 <- q; push the previous step size
    tIWG->shift_fields();

    // 2. Set the size of the step about to be taken
    tTime += tDeltaTime;
    tIWG->delta_time() = tDeltaTime;

    // 3. Nonlinear iteration within the step
    for ( uint tIter = 0; tIter < tMaxIter; ++tIter )
    {
        tField->compute_jacobian_and_rhs();   // compute_mkf + scheme, per element
        tField->solve();
        if ( tField->residual( tIter ) < tTolerance ) break;
    }
}
```

**Field shifting:**
- `shift_fields()`: x_old = x (start of a step)
- `reset_fields()`: x = x_old; reverses exactly one `shift_fields()` (reject the step, retry with a smaller dt)
- `make_savepoint()` / `restore_savepoint()`: undo any number of shifts (sub-stepped thermal equation)

---

## Weak Form Assembly

The IWG module implements the **weak form** (variational form) of governing PDEs.

### Mathematical Background

**Strong form** (classical PDE):
```
ρcp ∂T/∂t - ∇·(λ∇T) = q   in Ω
T = T₀                      on ∂Ω_D (Dirichlet)
-λ ∇T·n = h                 on ∂Ω_N (Neumann)
```

**Weak form** (multiply by test function v, integrate by parts):
```
∫_Ω v ρcp ∂T/∂t dV + ∫_Ω ∇v · λ∇T dV = ∫_Ω v q dV + ∫_∂Ω_N v h dS
```

**Discrete weak form** (Galerkin FEM, T ≈ N T_e):
```
M_ij = ∫_Ω N_i ρcp N_j dV           (mass matrix)
K_ij = ∫_Ω ∇N_i · λ∇N_j dV          (stiffness matrix)
f_i  = ∫_Ω N_i q dV + ∫_∂Ω_N N_i h dS  (force vector)

Time-discrete (BDF1):
(M/Δt + K) T^(n+1) = f + M T^n / Δt
(assembled multiplied by Δt: (M + Δt K) T^(n+1) = Δt f + M T^n)
```

### Element-Level Assembly (IWG Implementation)

**Pattern for compute_mkf():** the constructor must have flagged M, K and F, and
`link_to_group()` must have called `mTimeStepMatrices->initialize( n )` (see
[Workflow](#workflow)).

```cpp
void IWG_CustomPhysics::compute_mkf( Element * aElement )
{
    // the element matrices live on the IWG, not in the signature
    Matrix< real > & tM = mTimeStepMatrices->M();
    Matrix< real > & tK = mTimeStepMatrices->K();
    Vector< real > & tF = mTimeStepMatrices->f();

    // 1. Link the group's calculator to this element
    mCalc->link( aElement );

    // 2. Zero the flagged matrices
    mTimeStepMatrices->reset();

    // 3. Nodal temperatures (calculator scratch, no allocation)
    const Vector< real > & tTnodes = mCalc->node_data( "T" );

    // density of the undeformed mesh, once per element
    const real tRho = mMaterial->ref_density();

    // 4. Gauss quadrature loop
    const Vector< real > & tW = mCalc->integration()->weights();

    for ( uint k = 0; k < mCalc->num_intpoints(); ++k )
    {
        const Matrix< real > & N = mCalc->N( k );   // (1 × num_nodes)
        const Matrix< real > & B = mCalc->B( k );   // (ndim × num_nodes)
        const real tdV = tW( k ) * mCalc->dV( k );

        // Interpolate T at the Gauss point
        const real tT = mCalc->node_interp( k, tTnodes );

        // Material properties at tT
        const real tRhoCp  = tRho * mMaterial->cp( tT );
        const real tLambda = mMaterial->lambda( tT );
        const real tQ      = 1000.0;  // heat source

        tM += trans( N ) * N * ( tRhoCp * tdV );
        tK += trans( B ) * B * ( tLambda * tdV );
        tF += trans( N ) * ( tQ * tdV );
    }
}
```

### Nfunction and Bfunction

**Nfunction** determines shape function dimensionality:
```cpp
enum class Nfunction
{
    Scalar,   // N: (1 × num_nodes) for scalar fields (T, φ, p)
    Vec2d,    // N: (2 × 2*num_nodes) for 2D vector fields (u_x, u_y)
    Vec3d,    // N: (3 × 3*num_nodes) for 3D vector fields (u_x, u_y, u_z)
    UNDEFINED
};
```

**Bfunction** determines derivative operator:
```cpp
enum class Bfunction
{
    Gradient,     // B = ∇N for diffusion: -∇·(κ∇T)
    Planestress,  // B for plane stress (εxx, εyy, γxy from u, v)
    Voigt,        // B for 3D elasticity (6 strain components from u, v, w)
    UNDEFINED
};
```

The `Vec2d`/`Vec3d` and `Planestress`/`Voigt` operators exist in the calculator, but no IWG class
in the tree uses them.

**Example - Scalar diffusion:**
```cpp
// N: (1 × 4) for TET4
// B: (3 × 4) for 3D gradient

// Stiffness: ∫ B^T λ B dV
K += w * trans( B ) * lambda * B * dV;
// Result: K is (4 × 4) symmetric
```

---

## Nonlinear Solvers

BELFEM supports two nonlinear iteration methods:

### Newton-Raphson

**Algorithm:**
```
For iteration k = 0, 1, 2, ...
  1. Assemble Jacobian J(x^k) and residual R(x^k)
  2. Solve J(x^k) Δx = -R(x^k)
  3. Update x^(k+1) = x^k + ω Δx  (ω = fixed relaxation factor)
  4. Check convergence: residual( k ) < tol
```

`DofManager::residual()` returns a **relative** residual: the norm of the solver's current
right-hand-side vector divided by a reference norm, computed on rank 0 and broadcast.
`absolute_residual()` returns the unscaled norm.

ω is a fixed relaxation factor, not a line search: `set_omega()` stores it unchanged (it does not
clamp it), and 1 is an undamped step. `FEM::Controller` adapts ω between trials on its own.

**Advantages:**
- Quadratic convergence near solution, when the IWG supplies the full tangent
- Fewer iterations than Picard

**Disadvantages:**
- Requires Jacobian computation (expensive)
- May diverge if initial guess is poor

**BELFEM usage:**
```cpp
tIWG->set_algorithm( SolverAlgorithm::NewtonRaphson );  // Iterative-mode IWGs only

DofManager * tField = tKernel.create_field( tIWG );
tIWG->set_omega( 0.9 );  // after create_field, on every rank
tField->set_solver( SolverParameters( SolverType::MUMPS ) );
tField->initialize();

for ( uint tIter = 0; tIter < tMaxIter; ++tIter )
{
    // Assemble J and R
    tField->compute_jacobian_and_rhs();

    // Solve and update; the update is scaled by ω internally
    tField->solve();

    // Check convergence
    if ( tField->residual( tIter ) < tTolerance ) break;
}
```

**Jacobian assembly for nonlinear K(x):**

With K(x) = ∫ Bᵀ λ(T) B dV and T = N x, the derivative of the residual term K(x)·x has two
parts. The term K·Δx is the stiffness itself. The term (∂K/∂x)·x comes from λ depending on T.
The thermal kernel `T_h_newton` (`src/fem/thermal/matrices/mt_thermal_h.cpp`) assembles that
second part as the mixed-operator block

```
(∂K/∂x)·x = ∫ ( Bᵀ (B x) ) ⊗ ( dλ/dT · N ) dV
```

which is unsymmetric:

```cpp
// from mt_thermal_h.cpp, T_h_newton ( Btg is an n × 1 calculator workspace )
Btg = trans( B ) * ( B * aCalc->q() );
aMatrices->dKdx_times_x() += Btg * N * ( dlambdadT * wdV );
```

The IWG flags `dKdX_times_x` in its constructor. `assemble_dJdx` adds the block, scaled by h, to
`dJdx()`. The DofManager assembles `dJdx()` next to the scheme's element Jacobian (for BDF2,
J = α M + h K).

---

### Picard Iteration

**Algorithm:**
```
For iteration k = 0, 1, 2, ...
  1. Assemble K(x^k) using x^k
  2. Solve (M/dt + K(x^k)) x^(k+1) = f + M x^n / dt
  3. Check convergence: ||x^(k+1) - x^k|| < tol
```

**Advantages:**
- Simpler than Newton (no Jacobian derivatives)
- More robust for poor initial guesses

**Disadvantages:**
- Linear convergence (slower than Newton)
- May require many iterations

**BELFEM usage:**
```cpp
tIWG->set_algorithm( SolverAlgorithm::Picard );   // Iterative-mode IWGs only

for ( uint tIter = 0; tIter < tMaxIter; ++tIter )
{
    // Assemble K(x^k)
    tField->compute_jacobian_and_rhs();

    // Solve (M + h K) x^(k+1) = RHS
    tField->solve();

    // Check convergence
    if ( tField->residual( tIter ) < tTolerance ) break;
}
```

---

### Hybrid Picard-Newton Strategy

BELFEM papers (Messe et al. 2023, Section 2.7) recommend **hybrid approach** for HTS electromagnetics:

1. **Picard stage** (5-10 iterations): Robust, gets close to solution
2. **Quasi-Newton stage** (optional): Approximate Jacobian
3. **Full Newton stage**: Quadratic convergence to tight tolerance (ε < 10⁻¹¹)

```cpp
// Stage 1: Picard
tIWG->set_algorithm( SolverAlgorithm::Picard );
for ( uint k = 0; k < 10; ++k )
{
    tField->compute_jacobian_and_rhs();
    tField->solve();
    if ( tField->residual( k ) < 1e-4 ) break;
}

// Stage 2: Newton
tIWG->set_algorithm( SolverAlgorithm::NewtonRaphson );
tIWG->set_omega( 0.9 );   // the IWG is linked: create_field has run
for ( uint k = 0; k < 50; ++k )
{
    tField->compute_jacobian_and_rhs();
    tField->solve();
    if ( tField->residual( k ) < 1e-11 ) break;  // Tight tolerance!
}
```

**Why tight tolerance (ε < 10⁻¹¹)?**
From Messe et al. 2023: Prevents checkerboarding in HTS superconductor simulations.

---

## Factory Patterns

### IwgFactory Usage

`Kernel::create_equation()` builds an `IwgFactory` internally, selects all blocks when you pass no
block list, and takes ownership of the IWG. Use `IwgFactory` directly only when you need the IWG
before the kernel exists, and then hand it over:

```cpp
#include "cl_IwgFactory.hpp"

Mesh tMesh( "model.exo" );
IwgFactory tFactory( tMesh );

IWG * tIWG = tFactory.create_iwg(
    IwgType::StaticHeatConduction,
    ModelDimensionality::ThreeD );       // caller owns it

// All block IDs of the mesh (gathered on rank 0, broadcast)
tIWG->select_blocks( tFactory.all_block_ids() );

KernelParameters tParams( tMesh );
Kernel tKernel( &tParams );
tKernel.add_equation( tIWG );            // the kernel owns it from here on
DofManager * tField = tKernel.create_field( tIWG );
```

### Supported IWG Types (Factory Limitation)

⚠️ **Only 3 of the 10 concrete `IwgType`s have a factory branch:**

| IwgType | Factory Support | Alternative |
|---------|----------------|------------|
| Poisson | ✅ Yes | - |
| StaticHeatConduction | ✅ Yes | - |
| TransientHeatConduction | ✅ Yes (placeholder coefficients) | - |
| Gradient2D | ❌ No | None: no IWG class exists |
| Gradient3D | ❌ No | None: no IWG class exists |
| SurfaceGradient | ❌ No | None: no IWG class exists |
| PlaneStress | ❌ No | None: no IWG class exists |
| LinearElasticity | ❌ No | None: no IWG class exists |
| Maxwell | ❌ Explicitly blocked | `MaxwellFactory` |
| MaxwellThermal | ❌ Explicitly blocked | `ThermalFactory` |

**Manual IWG construction:** only for a class you wrote yourself (see
[Adding New IWG Types](#adding-new-iwg-types)):
```cpp
IWG * tIWG = new IWG_CustomPhysics( ModelDimensionality::ThreeD );
tIWG->select_blocks( { 1 } );
tKernel.add_equation( tIWG );                         // the kernel owns it
DofManager * tField = tKernel.create_field( tIWG );   // initializes and links it
```

---

## Material Property Integration

IWG implementations should use the materials module, not hardcoded values.

### Problem: Hardcoded Properties

`IWG_TransientHeatConduction` uses the constants `1.0` (mass) and `111.0` (conductivity) instead
of material data; see [Pitfall 5](#5-treating-iwg_transientheatconduction-as-a-working-solver).

### Solution: Material System

`mMaterial` is the material of the group the IWG is linked to. `IWG::link_to_group()` sets it
from `Group::material()`, so it changes with every block. `DofManager::block( id )->set_material(
label )` assigns it.

```cpp
// inside compute_mkf, after mCalc->link( aElement )
const Vector< real > & tTnodes = mCalc->node_data( "T" );

// density of the undeformed mesh ( Material::ref_density ), not rho():
// Material::rho( T ) is the electrical resistivity
const real tRho = mMaterial->ref_density();

for ( uint k = 0; k < mCalc->num_intpoints(); ++k )
{
    // Interpolate temperature at Gauss point
    const real tT = mCalc->node_interp( k, tTnodes );

    // Evaluate material properties at tT
    const real tCp     = mMaterial->cp( tT );       // Specific heat
    const real tLambda = mMaterial->lambda( tT );   // Thermal conductivity

    const real tdV = tW( k ) * mCalc->dV( k );
    tM += trans( mCalc->N( k ) ) * mCalc->N( k ) * ( tRho * tCp * tdV );
    tK += trans( mCalc->B( k ) ) * mCalc->B( k ) * ( tLambda * tdV );
}
```

### Material Derivatives for Newton-Raphson

```cpp
// For nonlinear Newton, the material supplies dλ/dT
const real tDlambdaDT = mMaterial->dlambdadT( tT );

// conductivity tangent: ( Bᵀ·∇T ) ⊗ ( dλ/dT · N ), see Nonlinear Solvers.
// tBtg is an n × 1 workspace sized once per group, not per call
tBtg = trans( B ) * ( B * mCalc->q() );
mTimeStepMatrices->dKdx_times_x() += tBtg * N * ( tDlambdaDT * tdV );
```

### Material Types

See `src/physics/materials/doc/` for material property framework:
- Pure metals (Cu, Al, etc.)
- HTS materials (REBCO, Bi-2223)
- User-defined materials

---

## Performance Considerations

### Critical: Let the DofManager Drive Assembly

The IWG is created once per equation and stays alive for the whole run. `DofManager` links it
to one group at a time (`link_to_group()`), and each group's Calculator holds the integration data
for its element type. Creating or deleting IWGs during assembly gains nothing and breaks the
kernel's ownership.

---

### Element Matrix Assembly Complexity

| Operation | Complexity | Notes |
|-----------|------------|-------|
| **Shape evaluation** | O(n_gp × n_b) | Cached in IntegrationData |
| **Derivative evaluation** | O(n_gp × n_b × n_dim) | Cached in IntegrationData |
| **Mass matrix** | O(n_gp × n_b²) | Dense BLAS gemm |
| **Stiffness matrix** | O(n_gp × n_dim × n_b²) | Dense BLAS gemm |
| **Material property eval** | O(n_gp) | Typically cheap (polynomial) |

Where:
- n_gp = number of Gauss points (4-27)
- n_b = number of bases (4-27 for Lagrange)
- n_dim = spatial dimensions (2 or 3)

**Bottleneck:** Not the IWG itself, but **global sparse assembly** (O(nnz) per element) and **sparse solve** (O(nnz^1.5) for direct solvers).

---

### Flag-Based Assembly

`reset()`, `initialize()` and `assemble_dJdx()` touch only the matrices the IWG flagged in its
constructor. An IWG without a Newton tangent therefore pays nothing for the derivative blocks.
No benchmark of this saving exists in the tree.

---

### Memory Footprint

There is one `TimestepMatrices` per IWG, not per element. It is re-sized for each group in
`link_to_group()`. For TET10 (10 DOFs per element):
```cpp
M:  10 × 10 × sizeof(real) = 800 bytes (dense)
K:  10 × 10 × sizeof(real) = 800 bytes
f:  10 × sizeof(real) = 80 bytes
J, dJdx: 800 bytes each
```

The DofManager assembles element by element, so the only large storage is the **global sparse
Jacobian**.

---

## Thread Safety and MPI

### Thread Safety

**IWG:**
- ❌ **Not thread-safe**. The IWG writes into its own `TimestepMatrices` and into the scratch
  buffers of the group's Calculator.
- ❌ **Per-thread IWGs do not fix this**. `link_to_group()` points every IWG on a group at the
  *same* Calculator (`mCalc = aGroup->calculator()`). A per-thread IWG that skipped
  `link_to_group()` has no calculator at all.

BELFEM's parallelism is MPI (see `doc/parallel_execution.md`). There is no supported threaded
element assembly. **DofManager** is not thread-safe during assembly either.

---

### MPI Compatibility

**Every rank runs the same sequence of calls:**
- `Mesh( path )` is collective: only the master reads the file; the others get an empty mesh
- The `Kernel` constructor is collective. It partitions the mesh on the master and distributes
  the submeshes
- Each rank's kernel creates its own IWG instance
- `Kernel::create_field()`, `IWG::set_omega()` and `DofManager::set_solver()` are collective;
  every rank must call them
- Element assembly needs no communication; DofManager handles the MPI assembly and the
  distributed solve

**Typical MPI workflow:**
```cpp
// every rank: the master reads the file, the Kernel distributes the submeshes
Mesh tMesh( "model.exo" );

KernelParameters tParams( tMesh );
Kernel tKernel( &tParams );                            // collective

// every rank creates its own IWG; the kernel owns it
IWG * tIWG = tKernel.create_equation( IwgType::Poisson, ModelDimensionality::ThreeD );
tIWG->select_blocks( { 1 } );
tIWG->select_sidesets( { 1, 2 } );

DofManager * tField = tKernel.create_field( tIWG );    // collective
tField->set_solver( SolverParameters( SolverType::MUMPS ) );
tField->sideset( 1 )->impose_dirichlet( 0.0 );
tField->sideset( 2 )->impose_dirichlet( 1.0 );
tField->initialize();  // Creates global DOF numbering via MPI

// Each rank assembles its local elements
tField->compute_jacobian();

// Distributed solve (e.g., MUMPS, PETSc)
tField->solve();

// ~Kernel deletes the DofManager and the IWG
```

---

## Common Patterns

### Pattern 1: Static Linear Problem

```cpp
IWG * tIWG = tKernel.create_equation(
    IwgType::Poisson,
    ModelDimensionality::ThreeD );

tIWG->select_blocks( { 1 } );
tIWG->select_sidesets( { 1, 2 } );

DofManager * tField = tKernel.create_field( tIWG );  // initializes the IWG and links the two
tField->set_solver( SolverParameters( SolverType::MUMPS ) );
tField->sideset( 1 )->impose_dirichlet( 0.0 );
tField->sideset( 2 )->impose_dirichlet( 1.0 );
tField->initialize();  // creates the fields, allocates DOFs and matrices

// Assemble and solve once
tField->compute_jacobian();
tField->solve();

// ~Kernel deletes the DofManager and the IWG
```

---

### Pattern 2: Transient Nonlinear Problem

Use the loop in [Time-Stepping Workflow](#time-stepping-workflow) with an IWG whose
`compute_mkf` fills its flagged matrices. `IWG_TransientHeatConduction` does not qualify (see
[Known Limitations](#known-limitations)).

---

### Pattern 3: Custom IWG Implementation

See [Adding New IWG Types](#adding-new-iwg-types) for the class skeleton and
[Element-Level Assembly](#element-level-assembly-iwg-implementation) for the body of
`compute_mkf`. Keep per-call scratch (`Matrix`/`Vector`) out of the quadrature loop. Size it once
as a member, or use a named calculator workspace (`mCalc->matrix( "label" )`, as `T_h_newton`
does).

---

## Known Limitations

### `IWG_TransientHeatConduction` Uses Placeholder Coefficients

**Location:** `cl_IWG_TransientHeatConduction.cpp`, `IWG_TransientHeatConduction::compute_mkf`; search for `TODO(cm): the mass matrix uses a unit coefficient`.

**Current behaviour:**

- The mass matrix uses the constant `1.0` instead of ρ·cp.
- The stiffness matrix uses the constant `111.0` instead of λ(T).
- The load vector `f` is zeroed, and nothing is added to it.
- The constructor sets no `MatrixFlag`, and nothing calls `TimestepMatrices::initialize()` for this class. Therefore, `M` and `K` are not sized before `compute_mkf` accumulates into them.

The class is reachable through `IwgFactory` and `Kernel::create_equation()`. `IWG_MaxwellThermal` derives from it but overrides `compute_mkf`, sets its own flags, and sizes its matrices in `link_to_group()`. The Maxwell thermal path therefore does not call this implementation of `compute_mkf`.

### `IWG_StaticHeatConduction`

- Every selected block needs a material because `compute_conduction` calls `mMaterial->lambda( T )` without a null check.
- The volume right-hand side is zero, so there is no heat source, and no dλ/dT tangent is formed.
- `compute_mkf` is not implemented.
- On a sideset with `SideSet::impose_alpha( aAlpha, aTinf )`, `compute_alpha_boundary_condition` reads the nodal fields `"alpha"` and `"Tinf"`. `SideSet::set_boundary_conditions()` writes `aAlpha` and `aTinf` to those fields on the sideset’s nodes. If `aTinf` is omitted (NaN), it becomes `BELFEM_TREF`; if `aAlpha` is omitted, it remains NaN.

---

## Development Notes

### Adding New IWG Types

**Steps:**

1. **Create derived class:**
   ```cpp
   class IWG_CustomPhysics : public IWG_Timestep
   {
   public:
       IWG_CustomPhysics( const ModelDimensionality aDimensionality,
                          const IwgType aType = IwgType::CustomPhysics,
                          const IwgMode aMode = IwgMode::Iterative )
       : IWG_Timestep( aType, aDimensionality, aMode )
       {
           mNumberOfDofsPerNode = 1;
           mNumberOfSpatialDimensions = 3;
           mNumberOfDerivativeDimensions = 3;
           mDofFields = { "T" };

           this->initialize( IwgType::CustomPhysics );   // skipped by children

           // declare the matrices compute_mkf produces
           mTimeStepMatrices->set_flag( MatrixFlag::M );
           mTimeStepMatrices->set_flag( MatrixFlag::K );
           mTimeStepMatrices->set_flag( MatrixFlag::F );
       }

       void link_to_group( Group * aGroup ) override;   // IWG::link_to_group, then, for a
                                                         // non-empty group,
                                                         // mTimeStepMatrices->initialize( n )
   protected:
       void compute_mkf( Element * aElement ) override; // fills matrices()->M()/K()/f()
   };
   ```

2. **Add to enum** (`en_IWGs.hpp`):
   ```cpp
   enum class IwgType { ..., CustomPhysics, UNDEFINED };
   ```

3. **Update factory** (`cl_IwgFactory.cpp`):
   ```cpp
   case IwgType::CustomPhysics:
       return new IWG_CustomPhysics( aDimensionality );
   ```

4. **Implement weak form** in `compute_mkf()`.

5. **Test:**
   - Manufactured solution
   - Mesh convergence (h-refinement)
   - Time-step convergence (for transient)

---

### Testing Shape Functions

**Partition of unity test:**
```cpp
Calculator* calc = iwg->calc();
calc->link( element );

for ( uint k = 0; k < calc->num_intpoints(); ++k )
{
    const Matrix<real>& N = calc->N( k );

    real sum = 0.0;
    for ( uint i = 0; i < N.n_cols(); ++i )
    {
        sum += N( 0, i );
    }

    BELFEM_ASSERT( std::abs( sum - 1.0 ) < 1e-12, "Partition of unity failed at point %u", k );
}
```

---

## See Also

### Internal Documentation

- **Interpolation**: `../../interpolation/doc/` - Shape functions, IntegrationData
- **DofManager**: `../../kernel/doc/` - Global assembly, DOF management
- **Materials**: `../../../physics/materials/doc/` - Material property framework
- **Mesh**: `../../../mesh/doc/` - Element definitions, blocks, sidesets

### Literature

**Finite Element Theory:**
- **Zienkiewicz & Taylor**, "The Finite Element Method" Vol. 1-2
- **Hughes**, "The Finite Element Method", Ch. 9 (Transient analysis)
- **Bathe**, "Finite Element Procedures", Ch. 8-9

**Time Integration:**
- **Hairer & Wanner 1996**, "Solving Ordinary Differential Equations II" (BDF stability)
- **Gear 1971**, "Numerical Initial Value Problems in Ordinary Differential Equations"

**BELFEM Papers:**
- **Messe et al. 2023** - Section 2.7 (Nonlinear iteration strategy)
- **Arsenault et al. 2023** - Magnetodynamic coupling (read with the **Arsenault et al. 2026** erratum, which corrects the air-domain form)

### Project References

- **Project README**: `../../../../README.md`
- **Claude Instructions**: `../../../../CLAUDE.md`
- **Coding Philosophy**: `../../../../doc/coding_philosophy.md`
- **Documentation Guidelines**: `../../../../doc/documentation_guidelines.md`

---

**Contributors:** Based on code by Christian Messe and Gregory Giard

**Prepared by the BELFEM development team.**
