# IWG Module Usage Guide {#fem_iwg_iwg_usage_guide}

**Date:** 2026-01-20
**Module:** src/fem/iwg
**Purpose:** Comprehensive guide to BELFEM's Integral Weak Form (IWG) physics module for finite element analysis

**Revision History:**

| Date | Version | Changes |
|------|---------|---------|
| 2026-01-20 | 1.1 | Updated to reflect critical bug fixes (BDF5, assemble_J) |
| 2026-01-20 | 1.0 | Initial comprehensive guide |

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
18. [Known Issues and Historical Bug Record](#known-issues-and-historical-bug-record)
19. [Development Notes](#development-notes)

---

## Common Pitfalls

**Place pitfalls first** - these are the most common mistakes that lead to crashes, incorrect results, or performance issues.

### 1. Using BDF5 Time Stepping (FIXED - 2026-01-20)

✅ **This bug has been fixed.** BDF5 is now safe to use.

**What was wrong:** The `mBeta` coefficient array was allocated with 4 elements but BDF5 accessed `mBeta(4)`, causing out-of-bounds access.

**Location (historical):** the pre-refactor `cl_IWG_Timestep.cpp`; in the
current file, search for `mBeta.set_size` — the allocation is 5.

**The fix applied:**
```cpp
// Changed from:
mBeta.set_size( 4, BELFEM_QUIET_NAN );  // Indices 0-3 only

// To:
mBeta.set_size( 5, BELFEM_QUIET_NAN );  // Indices 0-4 for BDF5 ✅
```

**BDF5 is now safe to use for high-order temporal accuracy:**
```cpp
// NOW SAFE - 5th order accuracy
iwg->set_timestepping_method( EulerMethod::BackwardDifference5 );
```

---

### 2. Trusting Assembled Jacobian (FIXED - 2026-01-20)

✅ **This bug has been fixed.** Jacobian assembly now works correctly.

**What was wrong:** The `assemble_J()` function reset `mJ` but then wrote to `mdJdx` instead, leaving `mJ` all zeros.

**Location (historical):** the pre-refactor `cl_TimestepMatrices.cpp`; in the
current file, `TimestepMatrices::assemble_J` correctly assembles `mJ`.

**The fix applied:**
```cpp
void TimestepMatrices::assemble_J( const real adt )
{
    mJ.fill( 0.0 );

    if ( mFlags.test( static_cast< index_t >( MatrixFlag::M )) )
    {
        mJ += mM;  // ✅ FIXED (was: mdJdx += mM)
    }

    if ( mFlags.test( static_cast< index_t >( MatrixFlag::K )) )
    {
        mJ += mK*adt;  // ✅ FIXED (was: mdJdx += mK*adt)
    }
}
```

**Newton-Raphson and transient solvers now work correctly:**
```cpp
// NOW SAFE - Jacobian assembly is correct
TimestepMatrices* matrices = iwg->matrices();
matrices->assemble_J( dt );
const Matrix<real>& J = matrices->J();  // Correctly populated ✅
```

---

### 3. Creating IWG Inside Element Loops

```cpp
// WRONG - Performance disaster!
IwgFactory factory( mesh );

for ( Element* e : elements )
{
    IWG* iwg = factory.create_iwg(        // ❌ Per element!
        IwgType::Poisson,
        ModelDimensionality::ThreeD
    );

    iwg->compute_jacobian( e, J );
    delete iwg;  // ❌ Allocation + deallocation per element
}
```

**Why wrong:** Factory creates new IWG, allocates IntegrationData, computes integration points, all per element.

**Impact:** 100-1000× slower than necessary for typical meshes.

```cpp
// CORRECT - Create once, reuse
IwgFactory factory( mesh );
IWG* iwg = factory.create_iwg(
    IwgType::Poisson,
    ModelDimensionality::ThreeD
);

iwg->select_blocks( { 1, 2, 3 } );
iwg->set_field( dofManager );
iwg->initialize();

for ( Element* e : elements )
{
    iwg->compute_jacobian( e, J );  // Reuses precomputed integration data
}

delete iwg;
```

---

### 4. Forgetting to Reset Flags Before Assembly

```cpp
TimestepMatrices* matrices = iwg->matrices();

for ( Element* e : elements )
{
    // WRONG - Missing reset_flags()
    iwg->compute_mkf( e );  // May use stale flags from previous element
    matrices->assemble_J( dt );
}
```

**Why wrong:** `MatrixFlag` bitset is not cleared, so assembly may include matrices from previous element.

**Impact:** Incorrect assembly; results depend on element ordering.

```cpp
// CORRECT - Reset before each element
TimestepMatrices* matrices = iwg->matrices();

for ( Element* e : elements )
{
    matrices->reset_flags();  // ✅ Clear all flags
    iwg->compute_mkf( e );
    matrices->assemble_J( dt );
}
```

---

### 5. Hardcoding Material Properties

```cpp
// WRONG - from cl_IWG_TransientHeatConduction.cpp (search "todo:: replace by cp and lambda")
aM += w( k ) * trans( N ) * 1.0 *  N * aCalc->dV( k );    // cp = 1.0 hardcoded!
aK += w( k ) * trans( B ) * 111.0 * B * aCalc->dV( k );   // lambda = 111.0 hardcoded!
```

**Why wrong:** Material properties are hardcoded test values, not actual material data.

**Impact:** `IWG_TransientHeatConduction` produces incorrect results for any real material.

```cpp
// CORRECT - Use material system
real tT = aCalc->node_interp( k, tTnodes );  // Temperature at Gauss point
aM += w( k ) * trans( N ) * aCalc->material()->cp( tT ) * N * aCalc->dV( k );
aK += w( k ) * trans( B ) * aCalc->material()->lambda( tT ) * B * aCalc->dV( k );
```

---

### 6. `is_maxwell()` must list every Maxwell variant (FIXED)

`is_maxwell()` gates the factories: `IwgFactory` refuses Maxwell types, `MaxwellFactory`
accepts them. It covers `Maxwell` and `MaxwellThermal` (`en_IWGs.hpp:99-102`); the
`UNDEFINED` line of the enum reminds you to extend it when adding a Maxwell IWG.

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

    // Not in the factory: Gradient2D, Gradient3D, SurfaceGradient,
    //                     PlaneStress, LinearElasticity

    default:
        BELFEM_ERROR( false, "invalid type" );  // explicit, not silent
}
```

**Impact:** the factory covers the three heat-conduction/Poisson types; Maxwell
types are rejected up front ("Use MaxwellFactory instead"), and everything
else hits the explicit error above.

**Workaround:** construct unlisted IWG types manually.

---

## Safe Defaults Quick Start

**For new users - minimal safe configuration to get started:**

```cpp
// 1. Create the kernel and let it create the IWG (it uses IwgFactory internally)
KernelParameters params( &mesh );
Kernel kernel( &params );
IWG* iwg = kernel.create_equation(
    IwgType::TransientHeatConduction,  // ✅ Supported in factory
    ModelDimensionality::ThreeD );     // no mode argument

// 2. Configure time-stepping (SAFE choices)
iwg->select_blocks( { 1, 2 } );
iwg->set_timestepping_method( EulerMethod::BackwardDifference2 );  // ✅ BDF2: stable, 2nd order
iwg->set_algorithm( SolverAlgorithm::NewtonRaphson );
iwg->set_omega( 0.9 );  // Line search relaxation

// 3. Create the DofManager (initializes the IWG and links the two)
DofManager* dofManager = kernel.create_field( iwg );  // initializes the IWG and links the two
dofManager->set_solver( SolverParameters( SolverType::MUMPS ) );
dofManager->initialize();  // creates the fields, allocates DOFs and matrices

// 4. Time loop
real dt = 0.01;
for ( uint step = 0; step < numSteps; ++step )
{
    iwg->delta_time() = dt;
    dofManager->compute_jacobian_and_rhs();
    dofManager->solve();
    iwg->shift_fields();
}

// ~Kernel deletes the DofManager and the IWG
```

### What's Currently Safe to Use

| Component | Safe Choice | Why |
|-----------|-------------|-----|
| **Factory-supported IWGs** | Poisson, StaticHeatConduction, TransientHeatConduction | Only 3/11 types implemented |
| **Time-stepping** | BDF2 or BDF3 | Stable, accurate, variable step |
| **Solver** | NewtonRaphson with ω=0.9 | Quadratic convergence |
| **Initialization** | Call `iwg->initialize()` ONCE after setup | Caches integration data |

### What to Avoid

| Avoid | Why | Alternative |
|-------|-----|-------------|
| Creating IWG in element loops | 1000× slower | Create once, reuse |
| Factory for unsupported types | Runtime error | Manual construction |
| Forgetting `initialize()` | No cached integration data | Always call after `set_field()` |
| Missing `reset_flags()` | Stale matrix flags | Reset before each element |

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

The **IWG (Integral Weak form of Governing equations)** module implements the physics layer for BELFEM's finite element framework. It provides:

- **PDE formulations**: Poisson, heat conduction, elasticity, Maxwell electromagnetics
- **Time integration**: the BDF family (BDF1–BDF5), plus static and forward-Euler modes. Crank–Nicolson and Galerkin are present in the `EulerMethod` enum but **disabled** — selecting either raises `BELFEM_ERROR` (`cl_IWG_Timestep.cpp:88-94`), because the Newton tangent assembly supports BDF only
- **Nonlinear solvers**: Newton-Raphson, Picard, direct linear solves
- **Element assembly**: Mass (M), damping (D), stiffness (K), force (f) matrices
- **Jacobian linearization**: Automatic differentiation for nonlinear problems

### Purpose

The IWG module bridges the gap between mathematical formulation and discrete FEM implementation:

1. **Mathematical**: Weak form of governing PDE
2. **IWG**: Element-level integral evaluation
3. **Assembly**: Global sparse system construction
4. **Solver**: Iterative or direct solution

### Key Features

1. **Modular physics**: Each PDE has its own IWG class (Poisson, HeatConduction, etc.)
2. **Automatic time integration**: BDF methods with variable time-step support
3. **Nonlinear iteration**: Newton-Raphson with line search (ω relaxation)
4. **Material coupling**: Direct integration with materials module
5. **Dense element matrices**: Optimized for typical element sizes (4-27 DOFs)

### Design Philosophy

Following `doc/coding_philosophy.md`:
- **Performance first**: Dense element matrices, flag-based assembly avoids recomputation
- **Zero abstraction penalty**: Virtual function overhead only at element level, not Gauss point level
- **Manual memory management**: IWG owns IntegrationData; caller owns IWG
- **Explicit ownership**: TimestepMatrices owned by IWG_Timestep

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
│   └── IWG_TransientHeatConduction
└── (Maxwell IWGs in separate module: src/fem/maxwell)

TimestepMatrices (container for element matrices)
    └── Owned by IWG_Timestep

IwgFactory (factory for non-Maxwell IWGs)
```

### Key Enumerations

**IwgType** (Physics selection):
```cpp
enum class IwgType
{
    Poisson,                    // Scalar Poisson: ∇²φ = f
    StaticHeatConduction,       // Steady thermal: -∇·(λ∇T) = q
    TransientHeatConduction,    // Transient: ρcp ∂T/∂t - ∇·(λ∇T) = q
    Gradient2D,                 // Gradient operator (2D)
    Gradient3D,                 // Gradient operator (3D)
    SurfaceGradient,            // Surface gradient
    PlaneStress,                // 2D elasticity (plane stress)
    LinearElasticity,           // 3D elasticity
    Maxwell,                    // Electromagnetics (h-φ formulation)
    MaxwellThermal,             // Coupled EM + thermal
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

⚠️ **Important:** Only 3 of the 10 concrete `IwgType`s are currently implemented in `IwgFactory`:

| IwgType | Factory Support | Workaround |
|---------|----------------|------------|
| **Poisson** | ✅ Yes | - |
| **StaticHeatConduction** | ✅ Yes | - |
| **TransientHeatConduction** | ✅ Yes | - |
| Gradient2D | ❌ No | Manual construction |
| Gradient3D | ❌ No | Manual construction |
| SurfaceGradient | ❌ No | Manual construction |
| PlaneStress | ❌ No | Manual construction |
| LinearElasticity | ❌ No | Manual construction |
| **Maxwell** | ❌ Explicitly blocked | Use `MaxwellFactory` |
| **MaxwellThermal** | ❌ No | Use `MaxwellFactory` |

**Use factory-supported types for new projects**. Manual construction example:
```cpp
// For unsupported types:
IWG* iwg = new IWG_PlaneStress( ModelDimensionality::TwoD, ... );
iwg->set_field( dofManager );
iwg->initialize();
```

---

## Ownership and Lifetime

**Understanding memory ownership prevents leaks and double-frees:**

### Ownership Rules

| Object | Owned By | Lifetime | Deletion |
|--------|----------|----------|----------|
| **IWG** | Kernel (via `create_equation`/`add_equation`) or the user (manual `new`, never handed over) | Until `~Kernel` | Kernel destructor deletes it; the DofManager only borrows it |
| **TimestepMatrices** | IWG_Timestep | Same as IWG | IWG destructor frees |
| **IntegrationData** | IWG (Calculator) | Cached after `initialize()` | IWG destructor frees |
| **Calculator** | IWG | Same as IWG | IWG destructor frees |
| **DofManager** | Kernel (`create_field`) | Until `~Kernel` | Kernel destructor deletes it |
| **Mesh** | User | Entire program | User must delete |
| **Material** | Mesh or User | Depends on setup | Check setup |

### Common Ownership Patterns

**Pattern 1: Kernel owns IWG and DofManager (typical)**
```cpp
KernelParameters params( &mesh );
Kernel kernel( &params );
IWG* iwg = kernel.create_equation( ... );             // Kernel owns the IWG
DofManager* dofManager = kernel.create_field( iwg );  // Kernel owns the DofManager

// ... use dofManager ...

// ~Kernel deletes the DofManagers, then the IWGs
// DO NOT: delete iwg; delete dofManager;  ❌ Double-free!
```

**Pattern 2: User-built IWG handed to the Kernel (advanced)**
```cpp
IWG* iwg = new IWG_CustomPhysics( ModelDimensionality::ThreeD );
kernel.add_equation( iwg );                           // Kernel now owns it
DofManager* dofManager = kernel.create_field( iwg );

// an IWG never handed to the kernel stays the user's: delete it yourself
```

**Pattern 3: Per-thread IWGs (OpenMP)**
```cpp
#pragma omp parallel
{
    // Each thread creates its own IWG
    IWG* iwg_thread = factory.create_iwg( ... );
    iwg_thread->set_field( dofManager );
    iwg_thread->initialize();

    #pragma omp for
    for ( int i = 0; i < n; ++i )
    {
        iwg_thread->compute_jacobian( elements[i], J[i] );
    }

    delete iwg_thread;  // ✅ Thread cleans up its own IWG
}
```

### What Gets Invalidated When

| Call | What Happens | Must Re-initialize |
|------|--------------|-------------------|
| `iwg->select_blocks()` | Clears block selection | Call `initialize()` |
| `iwg->set_field()` | Links DofManager, material | Call `initialize()` |
| `iwg->initialize()` | Caches IntegrationData for all element types | - |
| `iwg->shift_fields()` | x_old ← x (time stepping) | - |
| `iwg->reset_fields()` | x ← x_old (reject timestep) | - |
| `dofManager->initialize()` | Creates the fields and allocates the DOF system | - |

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
- Set matrix flags to indicate what was computed
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

Solve ∇²φ = f in domain Ω with φ = 0 on boundary ∂Ω.

```cpp
#include "cl_Mesh.hpp"
#include "cl_FEM_Kernel.hpp"

// Load mesh
Mesh mesh( "poisson.exo" );

// Create the kernel and the IWG
KernelParameters params( &mesh );
Kernel kernel( &params );
IWG* iwg = kernel.create_equation(
    IwgType::Poisson,
    ModelDimensionality::ThreeD  // Linear problem
);

// Select domain
iwg->select_blocks( { 1 } );  // All elements in block 1

// Create DOF manager (its field list comes from the IWG)
DofManager* dofManager = kernel.create_field( iwg );  // initializes the IWG and links the two
dofManager->set_solver( SolverParameters( SolverType::MUMPS ) );
// Set Dirichlet BC: phi = 0 on sideset 1
dofManager->sideset( 1 )->impose_dirichlet( 0.0 );

dofManager->initialize();  // creates the fields, allocates DOFs and matrices

// Solve
dofManager->compute_jacobian();   // Assemble stiffness matrix
dofManager->solve();               // Solve Kφ = f

// Extract solution
const Vector<real>& phi = dofManager->field_data( "phi" );

// Cleanup
// ~Kernel deletes the DofManager and the IWG
```

---

### Example 2: Transient Heat Conduction

Solve ρcp ∂T/∂t - ∇·(λ∇T) = q with BDF2 time-stepping.

```cpp
#include "cl_FEM_Kernel.hpp"

Mesh mesh( "thermal.exo" );

// Create the kernel and the IWG
KernelParameters params( &mesh );
Kernel kernel( &params );
IWG* iwg = kernel.create_equation(
    IwgType::TransientHeatConduction,
    ModelDimensionality::ThreeD  // Nonlinear (temperature-dependent properties)
);

iwg->select_blocks( { 1, 2 } );  // Conductor + insulator
iwg->set_algorithm( SolverAlgorithm::NewtonRaphson );

// Set time-stepping method
iwg->set_timestepping_method( EulerMethod::BackwardDifference2 );

// Create DOF manager (its field list comes from the IWG)
DofManager* dofManager = kernel.create_field( iwg );  // initializes the IWG and links the two
dofManager->set_solver( SolverParameters( SolverType::MUMPS ) );
dofManager->initialize();  // creates the fields, allocates DOFs and matrices

// Set initial condition: T = 300 K
mesh.field_data( "T" ).fill( 300.0 );

// Time-stepping loop
real dt = 0.01;  // Time step
real t = 0.0;
uint maxSteps = 1000;

for ( uint step = 0; step < maxSteps; ++step )
{
    t += dt;
    iwg->delta_time() = dt;

    // Newton-Raphson iteration
    dofManager->compute_jacobian_and_rhs();  // Assembles M/dt + K and RHS
    dofManager->solve();                      // Solve for increment

    // Shift fields: T_old = T for next timestep
    iwg->shift_fields();

    std::cout << "Step " << step << ", t = " << t << std::endl;
}

// ~Kernel deletes the DofManager and the IWG
```

---

### Example 3: Element-Level Matrix Assembly

Directly use IWG for element matrices (advanced use).

```cpp
IWG* iwg = factory.create_iwg( IwgType::StaticHeatConduction, ModelDimensionality::TwoD );
iwg->select_blocks( { 1 } );
iwg->set_field( dofManager );
iwg->initialize();

// Get first element
Element* element = mesh.block( 1 )->elements()( 0 );

// Get element matrix container
TimestepMatrices* matrices = iwg->matrices();

// Reset flags (critical!)
matrices->reset_flags();

// Compute element matrices
iwg->compute_mkf( element );  // Fills M, K, f; sets flags

// Check which matrices were populated
bool hasM = matrices->has_flag( MatrixFlag::M );
bool hasK = matrices->has_flag( MatrixFlag::K );
bool hasF = matrices->has_flag( MatrixFlag::F );

std::cout << "Mass matrix populated: " << hasM << std::endl;
std::cout << "Stiffness matrix populated: " << hasK << std::endl;
std::cout << "Force vector populated: " << hasF << std::endl;

// Access matrices
const Matrix<real>& K_elem = matrices->K();
const Vector<real>& f_elem = matrices->f();

std::cout << "Element stiffness matrix:\n" << K_elem << std::endl;
std::cout << "Element force vector:\n" << f_elem << std::endl;
```

---

### Example 4: Nonlinear Steady-State with Newton-Raphson

Solve nonlinear steady-state problem with Newton-Raphson.

```cpp
IWG* iwg = factory.create_iwg(
    IwgType::StaticHeatConduction,
    ModelDimensionality::ThreeD);

iwg->set_algorithm( SolverAlgorithm::NewtonRaphson );
iwg->set_omega( 0.9 );  // Relaxation parameter for line search

DofManager* dofManager = kernel.create_field( iwg );  // initializes the IWG and links the two
dofManager->set_solver( SolverParameters( SolverType::MUMPS ) );
dofManager->initialize();  // creates the fields, allocates DOFs and matrices

// Newton-Raphson iteration
real tol = 1e-6;
uint maxIter = 50;

for ( uint iter = 0; iter < maxIter; ++iter )
{
    // Assemble Jacobian and RHS
    dofManager->compute_jacobian_and_rhs();

    // Solve for increment: J * ΔT = -R
    dofManager->solve();

    // Check convergence
    real residual = dofManager->residual( iter );
    std::cout << "Iteration " << iter << ", residual = " << residual << std::endl;

    if ( residual < tol )
    {
        std::cout << "Converged in " << iter << " iterations" << std::endl;
        break;
    }
}

// ~Kernel deletes the DofManager and the IWG
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
- Avoids recomputing matrices that haven't changed
- Time-stepper knows which matrices to include in assembly
- Zero-overhead check: single bit test

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

**Step 1: Reset flags**
```cpp
TimestepMatrices* matrices = iwg->matrices();
matrices->reset_flags();  // Clear all 8 flag bits
```

**Step 2: IWG populates matrices**
```cpp
// Inside IWG::compute_mkf()
Calculator* calc = iwg->calc();
Matrix<real>& M = matrices->M();
Matrix<real>& K = matrices->K();
Vector<real>& f = matrices->f();

M.fill( 0.0 );
K.fill( 0.0 );
f.fill( 0.0 );

for ( uint k = 0; k < calc->num_intpoints(); ++k )
{
    const Matrix<real>& N = calc->N( k );
    const Matrix<real>& B = calc->B( k );
    real w = calc->integration()->weights()( k );
    real dV = calc->dV( k );

    // Mass matrix: ∫ ρcp N^T N dV
    M += w * trans( N ) * rho_cp * N * dV;

    // Stiffness matrix: ∫ B^T λ B dV
    K += w * trans( B ) * lambda * B * dV;

    // Force vector: ∫ N^T q dV
    f += w * trans( N ) * heat_source * dV;
}

// Set flags to indicate what was computed
matrices->set_flag( MatrixFlag::M );
matrices->set_flag( MatrixFlag::K );
matrices->set_flag( MatrixFlag::F );
```

**Step 3: Time-stepper assembles Jacobian**
```cpp
real dt = 0.01;
matrices->assemble_J( dt );  // J = M/dt + K (for BDF1)

// Or for BDF2 (variable step):
// J = α*M/h + K
// where α, h computed from timestep history
```

**Step 4: Global assembly**
```cpp
// DofManager accumulates element J into global sparse Jacobian
const Matrix<real>& J_elem = matrices->J();
globalJacobian->add_element_contribution( element, J_elem );
```

### API Reference

**Initialization:**
```cpp
void initialize( const index_t aNumDofs );  // Resize all matrices to n_e × n_e
void reset();                                // Zero all matrices, keep flags
void reset_flags();                          // Clear flag bitset
```

**Flag management:**
```cpp
void set_flag( const MatrixFlag aFlag );             // Set flag after populating matrix
bool has_flag( const MatrixFlag aFlag ) const;       // Check if matrix was populated
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
const Matrix<real>& J() const;             // Assembled Jacobian
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
void assemble_J( const real adt );      // J = f(M, K, dt) based on time method
void assemble_dJdx( const real adt );   // dJ/dx for Newton correction
```

---

## Time-Stepping Algorithms

The `EulerMethod` enum has 12 concrete schemes (two of them disabled) for transient PDEs of the form:

**M(x) dx/dt + K(x) x = f(x, t)**

Where:
- M = mass matrix (capacity)
- K = stiffness matrix (diffusion/conductivity)
- f = force vector (source term)

### `EulerMethod` values

All rows below exist in the enum. **Crank–Nicolson and Galerkin are disabled** — selecting either
raises `BELFEM_ERROR` (`cl_IWG_Timestep.cpp:88-94`), because the Newton tangent assembly supports
the BDF family only. They are listed for completeness, not as options.

| EulerMethod | Order | A-Stability | Formula | Use Case |
|-------------|-------|-------------|---------|----------|
| **Static** | 0 | - | K x = f | Steady-state |
| **ForwardExplicit** | 1 | ❌ Conditionally | M(x^(n+1) - x^n)/h = f^n - K^n x^n | Explicit (dt limited) |
| **BackwardDifference1** (BDF1) | 1 | ✅ A-stable | M(x^(n+1) - x^n)/h + K x^(n+1) = f^(n+1) | Implicit Euler |
| **CrankNicolson** | 2 | ✅ A-stable | M(x^(n+1) - x^n)/h + K(x^(n+1) + x^n)/2 = f^(n+½) | ❌ **disabled** — raises `BELFEM_ERROR` |
| **Galerkin** | 1 | ✅ A-stable | Discontinuous Galerkin variant | ❌ **disabled** — raises `BELFEM_ERROR` |
| **BackwardDifference2** (BDF2) | 2 | ✅ A-stable | M(α x^(n+1) - β₁x^n + β₀x^(n-1))/h + K x^(n+1) = f^(n+1) | Variable step BDF2 |
| **BackwardDifference3** (BDF3) | 3 | A(α)-stable (α ≈ 86°) | Variable step BDF3 | Higher accuracy |
| **BackwardDifference4** (BDF4) | 4 | A(α)-stable (α ≈ 73°) | Variable step BDF4 | Higher accuracy |
| **BackwardDifference5** (BDF5) | 5 | A(α)-stable (α ≈ 52°) | Variable step BDF5 | ✅ **NOW SAFE** (fixed 2026-01-20) |

### BDF Coefficient Formulas

All BDF methods support **variable time-steps**. Coefficients are computed automatically.

**BDF1 (Backward Euler):**
```
J = M/h + K
RHS = f + M x₀/h

Where:
  h = current timestep
  x₀ = solution at previous timestep
```

**BDF2 (variable step):**
```
Given:
  h = current timestep (t^(n+1) - t^n)
  h₀ = previous timestep (t^n - t^(n-1))

Coefficients:
  α = (2h + h₀) / (h + h₀)
  β₁ = (h + h₀) / h₀
  β₀ = h² / (h₀ (h + h₀))

Assembly:
  J = α M/h + K
  RHS = f + M (β₁ x^n − β₀ x^(n-1)) / h      (for h = h₀: 2 x^n − ½ x^(n-1))
```

**BDF3-5:** Similar pattern with more history terms. See `cl_IWG_Timestep.cpp:compute_bdf_coeffs_N()`.

**Coefficient lifecycle (added 2026-07-03):** `compute_bdf_coefficients()` is invoked
lazily from `IWG_Timestep::compute_jacobian_and_rhs()`. `shift_fields()` and
`reset_fields()` mark the coefficients dirty; the first element assembly of a step
(which always runs after `delta_time()` has been set) recomputes them once. No manual
call is needed.

**Startup order ramp (added 2026-07-03):** a BDF-p run needs p history states, but only
the initial condition exists at step 1. The scheme therefore ramps up automatically:
step 1 runs BDF1, step 2 BDF2, ..., until the requested order is reached. A rejected
timestep (`reset_fields()`) rolls the ramp back by one step.

### Newton Jacobian for Nonlinear Time-Stepping

For BDF schemes with nonlinear M(x), the Newton Jacobian includes:
- α·(∂M/∂x)·x^n: Current iterate derivative (scaled by BDF coefficient α)
- (∂M/∂x)·(Σ β_i·x^{n-i}): History term derivatives (not scaled)
- Δt·(∂K/∂x)·x^n: Stiffness derivative (scaled by time-step)
- Δt·(∂f/∂x): Load derivative (scaled by time-step)

Assembled in `TimestepMatrices::assemble_dJdx( adt, aAlpha )`; the α factor is passed
from `IWG_Timestep::compute_jacobian_and_rhs()`. For all schemes that scale M by one
(BDF1 and explicit Euler; also Crank–Nicolson and Galerkin, were they enabled), α = 1.

Reference: Hairer & Wanner (1996) II.4 for the BDF Jacobian structure.

### Selecting a Time-Stepping Method

```cpp
iwg->set_timestepping_method( EulerMethod::BackwardDifference2 );

// For problems with stiffness matrix:
iwg->set_timestepping_method( EulerMethod::BackwardDifference2, true );

// For problems with only mass matrix (e.g., pure diffusion with lumped M):
iwg->set_timestepping_method( EulerMethod::BackwardDifference2, false );
```

**With stiffness (`aHaveStiffness = true`):**
- Assembles J = α M/h + K
- Standard for diffusion-reaction, elastodynamics

**Without stiffness (`aHaveStiffness = false`):**
- Assembles J = M/h only
- Faster for pure transport, explicit-implicit hybrid

### Variable Time-Step Support

All BDF methods automatically adapt to changing time-steps:

```cpp
real dt = 0.01;  // Initial timestep

for ( uint step = 0; step < maxSteps; ++step )
{
    // Adaptive time-stepping
    if ( converged_quickly )
    {
        dt *= 1.5;  // Increase timestep
    }
    else if ( slow_convergence )
    {
        dt *= 0.5;  // Decrease timestep
    }

    iwg->delta_time() = dt;  // Update timestep

    // BDF coefficients automatically recomputed from timestep history
    dofManager->compute_jacobian_and_rhs();
    dofManager->solve();

    iwg->shift_fields();  // x^(n-1) = x^n, x^n = x^(n+1)
}
```

### Time-Stepping Workflow

**Initialization:**
```cpp
iwg->set_timestepping_method( EulerMethod::BackwardDifference2 );
mesh.field_data( "T" ).fill( T_initial );  // Set initial condition
```

**Time loop:**
```cpp
for ( uint step = 0; step < numSteps; ++step )
{
    // 1. Update time
    t += dt;
    iwg->delta_time() = dt;

    // 2. Assemble J and RHS
    //    - Calls compute_mkf() for each element
    //    - Assembles M/dt + K → J
    //    - Assembles RHS with history terms
    dofManager->compute_jacobian_and_rhs();

    // 3. Solve linear system J Δx = RHS
    dofManager->solve();

    // 4. Shift history: x^(n-1) = x^n, x^n = x^(n+1)
    iwg->shift_fields();
}
```

**Field shifting:**
- `shift_fields()`: x_old = x (prepare for next timestep)
- `reset_fields()`: x = x_old (reject timestep, retry with smaller dt)

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
```

### Element-Level Assembly (IWG Implementation)

**Pattern for compute_mkf():**

```cpp
void IWG_CustomPhysics::compute_mkf( Element * aElement )
{
    // the element matrices live on the IWG, not in the signature
    Matrix< real > & tM = this->matrices()->M();
    Matrix< real > & tK = this->matrices()->K();
    Vector< real > & tF = this->matrices()->f();

    // 1. Get calculator (wraps IntegrationData)
    Calculator* calc = this->calc();
    calc->link( aElement );  // Link to current element

    // 2. Collect node data
    Matrix< real > tNodeCoords;
    collect_node_coords( aElement, tNodeCoords );  // X, Y, Z

    Matrix< real > tNodeTemps;
    collect_node_data( aElement, { "T" }, tNodeTemps );  // T at nodes

    // 3. Zero matrices
    tM.fill( 0.0 );
    tK.fill( 0.0 );
    tF.fill( 0.0 );

    // 4. Gauss quadrature loop
    const Vector<real>& w = calc->integration()->weights();
    uint npts = calc->num_intpoints();

    for ( uint k = 0; k < npts; ++k )
    {
        // Shape functions
        const Matrix<real>& N = calc->N( k );      // (1 × num_nodes)
        const Matrix<real>& B = calc->B( k );      // (ndim × num_nodes)
        real dV = calc->dV( k );                   // Volume element

        // Interpolate T at Gauss point
        real T_gp = calc->node_interp( k, tNodeTemps.col(0) );

        // Material properties at T_gp
        const Material* mat = mMaterial;
        real rho_cp = mat->rho() * mat->cp( T_gp );
        real lambda = mat->lambda( T_gp );
        real q = 1000.0;  // Heat source (could be field)

        // Mass matrix: ∫ ρcp N^T N dV
        tM += w(k) * trans( N ) * rho_cp * N * dV;

        // Stiffness matrix: ∫ B^T λ B dV
        tK += w(k) * trans( B ) * lambda * B * dV;

        // Force vector: ∫ N^T q dV
        tF += w(k) * trans( N ) * q * dV;
    }

    // 5. Set flags
    TimestepMatrices* matrices = this->matrices();
    matrices->set_flag( MatrixFlag::M );
    matrices->set_flag( MatrixFlag::K );
    matrices->set_flag( MatrixFlag::F );
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

**Example - Scalar diffusion:**
```cpp
// N: (1 × 4) for TET4
// B: (3 × 4) for 3D gradient

// Stiffness: ∫ B^T λ B dV
K += w * trans( B ) * lambda * B * dV;
// Result: K is (4 × 4) symmetric
```

**Example - 3D elasticity:**
```cpp
// N: (3 × 12) for TET4 (3 DOFs per node)
// B: (6 × 12) for Voigt strain-displacement

// Stiffness: ∫ B^T D B dV
K += w * trans( B ) * D_elasticity * B * dV;
// D_elasticity: (6 × 6) constitutive matrix
// Result: K is (12 × 12) symmetric
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
  3. Update x^(k+1) = x^k + ω Δx  (ω = relaxation parameter)
  4. Check convergence: ||R(x^(k+1))|| < tol
```

**Advantages:**
- Quadratic convergence near solution
- Fewer iterations than Picard

**Disadvantages:**
- Requires Jacobian computation (expensive)
- May diverge if initial guess is poor

**BELFEM usage:**
```cpp
iwg->set_algorithm( SolverAlgorithm::NewtonRaphson );
iwg->set_omega( 0.9 );  // Line search relaxation (0 < ω ≤ 1)

DofManager* dofManager = kernel.create_field( iwg );
dofManager->set_solver( SolverParameters( SolverType::MUMPS ) );
dofManager->initialize();

for ( uint iter = 0; iter < maxIter; ++iter )
{
    // Assemble J and R
    dofManager->compute_jacobian_and_rhs();

    // Solve J Δx = -R
    dofManager->solve();  // x^(k+1) = x^k + ω Δx done internally

    // Check convergence
    real residual = dofManager->residual( iter );
    if ( residual < tol ) break;
}
```

**Jacobian assembly for nonlinear K(x):**

IWG computes **linearization** of K(x):
```cpp
// Stiffness depends on x (e.g., temperature-dependent conductivity)
K(x) = ∫ B^T λ(T(x)) B dV

// Linearization:
dK/dx = ∫ B^T (dλ/dT)(dT/dx) B dV
      = ∫ B^T (dλ/dT) N B dV  (since T = N x_T)

// Contracted form (what IWG computes):
(dK/dx) · x = ∫ B^T (dλ/dT) (N x) B dV
```

IWG populates `dKdX_times_x` matrix:
```cpp
for ( uint k = 0; k < npts; ++k )
{
    real T = calc->node_interp( k, T_nodes );
    real lambda = mat->lambda( T );
    real dLambda_dT = mat->dlambdadT( T );  // Material derivative
    real T_gp = calc->node_interp( k, T_nodes );

    // Standard stiffness
    K += w(k) * trans( B ) * lambda * B * dV;

    // Linearized derivative: (dK/dx) · x
    Matrix<real>& dKdx_times_x = matrices->dKdx_times_x();
    dKdx_times_x += w(k) * trans( B ) * dLambda_dT * T_gp * B * dV;
}

matrices->set_flag( MatrixFlag::K );
matrices->set_flag( MatrixFlag::dKdX_times_x );
```

Time-stepper assembles:
```cpp
J = M/dt + K + dK/dx  // Full Newton Jacobian
```

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
iwg->set_algorithm( SolverAlgorithm::Picard );

for ( uint iter = 0; iter < maxIter; ++iter )
{
    // Assemble K(x^k)
    dofManager->compute_jacobian_and_rhs();

    // Solve (M/dt + K) x^(k+1) = RHS
    dofManager->solve();

    // Check convergence
    real residual = dofManager->residual( iter );
    if ( residual < tol ) break;
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
iwg->set_algorithm( SolverAlgorithm::Picard );
for ( uint k = 0; k < 10; ++k )
{
    dofManager->compute_jacobian_and_rhs();
    dofManager->solve();
    if ( dofManager->residual( k ) < 1e-4 ) break;
}

// Stage 2: Newton
iwg->set_algorithm( SolverAlgorithm::NewtonRaphson );
iwg->set_omega( 0.9 );
for ( uint k = 0; k < 50; ++k )
{
    dofManager->compute_jacobian_and_rhs();
    dofManager->solve();
    if ( dofManager->residual( k ) < 1e-11 ) break;  // Tight tolerance!
}
```

**Why tight tolerance (ε < 10⁻¹¹)?**
From Messe et al. 2023: Prevents checkerboarding in HTS superconductor simulations.

---

## Factory Patterns

### IwgFactory Usage

```cpp
#include "cl_IwgFactory.hpp"

Mesh mesh( "model.exo" );
IwgFactory factory( &mesh );

// Create IWG by type
IWG* iwg = kernel.create_equation(
    IwgType::TransientHeatConduction,
    ModelDimensionality::ThreeD);

// Get all block IDs from mesh
const Vector<id_t>& blockIDs = factory.all_block_ids();

delete iwg;
```

### Supported IWG Types (Factory Limitation)

⚠️ **Only 3 of the 10 concrete `IwgType`s are currently implemented in factory:**

| IwgType | Factory Support | Workaround |
|---------|----------------|------------|
| Poisson | ✅ Yes | - |
| StaticHeatConduction | ✅ Yes | - |
| TransientHeatConduction | ✅ Yes | - |
| Gradient2D | ❌ No | Manual construction |
| Gradient3D | ❌ No | Manual construction |
| SurfaceGradient | ❌ No | Manual construction |
| PlaneStress | ❌ No | Manual construction |
| LinearElasticity | ❌ No | Manual construction |
| Maxwell | ❌ Explicitly blocked | Use MaxwellFactory |
| MaxwellThermal | ❌ No | Use MaxwellFactory |

**Manual IWG construction:**
```cpp
// For types not in factory
IWG* iwg = new IWG_CustomType(
    ModelDimensionality::ThreeD,
    IwgType::CustomType,
    IwgMode::Iterative
);

iwg->set_field( dofManager );
iwg->initialize();
```

---

## Material Property Integration

IWG implementations should use the materials module, not hardcoded values.

### Problem: Hardcoded Properties

```cpp
// WRONG - from cl_IWG_TransientHeatConduction.cpp (search "todo:: replace by cp and lambda")
aM += w( k ) * trans( N ) * 1.0 *  N * aCalc->dV( k );    // cp = 1.0
aK += w( k ) * trans( B ) * 111.0 * B * aCalc->dV( k );   // lambda = 111.0
```

### Solution: Material System

```cpp
// CORRECT - Use material properties
const Material* mat = mMaterial;  // Set via IWG::set_field()

for ( uint k = 0; k < npts; ++k )
{
    // Interpolate temperature at Gauss point
    real T = calc->node_interp( k, T_nodes );

    // Evaluate material properties at T
    real rho = mat->rho();              // Density
    real cp = mat->cp( T );             // Specific heat (temperature-dependent)
    real lambda = mat->lambda( T );     // Thermal conductivity

    // Assembly
    aM += w(k) * trans( N ) * (rho * cp) * N * dV;
    aK += w(k) * trans( B ) * lambda * B * dV;
}
```

### Material Derivatives for Newton-Raphson

```cpp
// For nonlinear Newton, compute dλ/dT
real lambda = mat->lambda( T );
real dLambda_dT = mat->dlambdadT( T );  // Material provides derivative

// Linearized stiffness derivative
Matrix<real>& dKdx_x = matrices->dKdx_times_x();
dKdx_x += w(k) * trans( B ) * dLambda_dT * T * B * dV;
```

### Material Types

See `src/physics/materials/doc/` for material property framework:
- Pure metals (Cu, Al, etc.)
- HTS materials (REBCO, Bi-2223)
- User-defined materials

---

## Performance Considerations

### Critical: Cache IWG, Reuse IntegrationData

**WRONG - Per-element factory:**
```cpp
IwgFactory factory( &mesh );

for ( Element* e : elements )
{
    IWG* iwg = factory.create_iwg( ... );  // ❌ 1000× slower
    iwg->compute_jacobian( e, J );
    delete iwg;
}
```

**CORRECT - Single IWG, cached integration:**
```cpp
IwgFactory factory( &mesh );
IWG* iwg = factory.create_iwg( ... );  // ✅ Once

iwg->select_blocks( blockIDs );
iwg->set_field( dofManager );
iwg->initialize();  // Caches IntegrationData per element type

for ( Element* e : elements )
{
    iwg->compute_jacobian( e, J );  // Reuses cached integration data
}

delete iwg;
```

**Performance gain:** 100-1000× for typical meshes.

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

### Flag-Based Assembly Avoids Recomputation

```cpp
TimestepMatrices* matrices = iwg->matrices();

// First element
matrices->reset_flags();
iwg->compute_mkf( elem1 );  // Sets M, K, F flags

if ( matrices->has_flag( MatrixFlag::M ) )
{
    // M was computed - use it
}

// Time-stepper checks flags before assembly
matrices->assemble_J( dt );  // Only includes M, K if flags set
```

**Benefit:** Avoids zeroing/summing matrices that weren't computed. ~10-20% faster for large element counts.

---

### Memory Footprint

**TimestepMatrices for TET10 (10 nodes, 10 DOFs):**
```cpp
// 10 DOFs per element
M:  10 × 10 × sizeof(real) = 800 bytes (dense)
K:  10 × 10 × sizeof(real) = 800 bytes
f:  10 × sizeof(real) = 80 bytes
Flags: 1 byte (bitset)

Total per element: ~2 KB
```

**For 100,000 elements:** ~200 MB for all element matrices (if stored, which DofManager doesn't do).

**Actual memory:** DofManager assembles on-the-fly, so only **global sparse Jacobian** is stored (~10-50 MB for typical problems).

---

## Thread Safety and MPI

### Thread Safety

**IWG:**
- ❌ **Not thread-safe** (shared Calculator, TimestepMatrices)
- ✅ **Solution:** Create separate IWG per thread

```cpp
// WRONG - Shared IWG
#pragma omp parallel for
for ( int i = 0; i < elements.size(); ++i )
{
    iwg->compute_jacobian( elements[i], J[i] );  // ❌ Race condition
}

// CORRECT - Per-thread IWG
#pragma omp parallel
{
    IWG* iwg_thread = factory.create_iwg( ... );
    iwg_thread->set_field( dofManager );
    iwg_thread->initialize();

    #pragma omp for
    for ( int i = 0; i < elements.size(); ++i )
    {
        iwg_thread->compute_jacobian( elements[i], J[i] );
    }

    delete iwg_thread;
}
```

**DofManager:**
- ❌ **Not thread-safe** during assembly
- Global sparse matrix assembly requires locking or coloring

**Best practice:** Let DofManager handle parallelism via block-based or coloring schemes.

---

### MPI Compatibility

**All IWG classes are MPI-aware:**
- Each rank creates independent IWG instances
- No communication required for element assembly
- DofManager handles MPI assembly (scatter/gather)

**Typical MPI workflow:**
```cpp
// Each rank has local mesh partition
Mesh mesh_local( comm );

// Each rank creates independent IWG
IWG* iwg = factory.create_iwg( ... );
iwg->select_blocks( blockIDs );

// DofManager handles MPI assembly
KernelParameters params( &mesh_local );
Kernel kernel( &params );
kernel.add_equation( iwg );                            // the kernel now owns the IWG
DofManager* dofManager = kernel.create_field( iwg );
dofManager->set_solver( SolverParameters( SolverType::MUMPS ) );
dofManager->initialize();  // Creates global DOF numbering via MPI

// Each rank assembles local elements
dofManager->compute_jacobian();  // Local assembly

// DofManager calls MPI solver (e.g., MUMPS, PETSc)
dofManager->solve();  // Distributed solve

// ~Kernel deletes the DofManager and the IWG
```

---

## Common Patterns

### Pattern 1: Static Linear Problem

```cpp
IWG* iwg = kernel.create_equation(
    IwgType::Poisson,
    ModelDimensionality::ThreeD);

iwg->select_blocks( { 1 } );

DofManager* dofManager = kernel.create_field( iwg );  // initializes the IWG and links the two
dofManager->set_solver( SolverParameters( SolverType::MUMPS ) );
dofManager->initialize();  // creates the fields, allocates DOFs and matrices

// Assemble and solve once
dofManager->compute_jacobian();
dofManager->solve();

// ~Kernel deletes the DofManager and the IWG
```

---

### Pattern 2: Transient Nonlinear Problem

```cpp
IWG* iwg = kernel.create_equation(
    IwgType::TransientHeatConduction,
    ModelDimensionality::ThreeD);

iwg->set_algorithm( SolverAlgorithm::NewtonRaphson );
iwg->set_timestepping_method( EulerMethod::BackwardDifference2 );
iwg->set_omega( 0.9 );

DofManager* dofManager = kernel.create_field( iwg );  // initializes the IWG and links the two
dofManager->set_solver( SolverParameters( SolverType::MUMPS ) );
dofManager->initialize();  // creates the fields, allocates DOFs and matrices

// Time-stepping
real dt = 0.01;
for ( uint step = 0; step < numSteps; ++step )
{
    iwg->delta_time() = dt;

    dofManager->compute_jacobian_and_rhs();
    dofManager->solve();

    iwg->shift_fields();
}

// ~Kernel deletes the DofManager and the IWG
```

---

### Pattern 3: Custom IWG Implementation

```cpp
class IWG_MyPhysics : public IWG_Timestep
{
public:
    IWG_MyPhysics(
        const ModelDimensionality aDimensionality,
        const IwgType aType = IwgType::CustomType,
        const IwgMode aMode = IwgMode::Iterative )
    : IWG_Timestep( aType, aDimensionality, aMode )
    {
        // Set DOF fields
        mDofFields = { "u", "v", "w" };  // 3D displacement
        mNumberOfDofsPerNode = 3;
        mNumberOfSpatialDimensions = 3;
        mNumberOfDerivativeDimensions = 6;  // Voigt notation
    }

protected:
    void compute_mkf( Element * aElement ) override
    {
        Matrix< real > & tM = this->matrices()->M();
        Matrix< real > & tK = this->matrices()->K();
        Vector< real > & tF = this->matrices()->f();

        Calculator* calc = this->calc();
        calc->link( aElement );

        tM.fill( 0.0 );
        tK.fill( 0.0 );
        tF.fill( 0.0 );

        const Vector<real>& w = calc->integration()->weights();
        uint npts = calc->num_intpoints();

        for ( uint k = 0; k < npts; ++k )
        {
            const Matrix<real>& N = calc->N( k );  // (3 × 3*num_nodes)
            const Matrix<real>& B = calc->B( k );  // (6 × 3*num_nodes) Voigt
            real dV = calc->dV( k );

            // Material properties
            real rho = mMaterial->rho();
            Matrix<real> D = elasticity_matrix( mMaterial );  // 6×6

            // Mass matrix (if needed for dynamics)
            tM += w(k) * trans( N ) * rho * N * dV;

            // Stiffness matrix: ∫ B^T D B dV
            tK += w(k) * trans( B ) * D * B * dV;

            // Body force
            Vector<real> f_body = { 0, 0, -9.81 * rho };  // Gravity
            tF += w(k) * trans( N ) * f_body * dV;
        }

        TimestepMatrices* matrices = this->matrices();
        matrices->set_flag( MatrixFlag::M );
        matrices->set_flag( MatrixFlag::K );
        matrices->set_flag( MatrixFlag::F );
    }
};
```

---

## Known Issues and Historical Bug Record

> **Note (2026-08-14):** this section is a dated record. The line numbers in
> the 2026-01-20 entries refer to the pre-refactor sources and are kept as
> written; the timestep machinery has since been rewritten (variable-step BDF
> family). Verified against the tree on 2026-08-14: items 1, 2, 3 and 5 are
> resolved in the current source; item 4 is the only one still open.

### Resolved (2026-01-20 record, re-verified 2026-08-14)

| # | Issue | Status in the current tree |
|---|-------|----------------------------|
| 1 | BDF5 `mBeta` out-of-bounds | Fixed — `mBeta.set_size( 5, BELFEM_QUIET_NAN )` in `cl_IWG_Timestep.cpp` (search `mBeta.set_size`) |
| 2 | `assemble_J()` wrote `mdJdx` instead of `mJ` | Fixed — `TimestepMatrices::assemble_J` assembles `mJ` from `mM` and `mK*adt` |
| 3 | `is_maxwell()` missing `MaxwellThermal` | Fixed — `en_IWGs.hpp`, `is_maxwell` returns true for both `Maxwell` and `MaxwellThermal` |
| 5 | `IwgFactory` switch silently incomplete | Resolved — the `default:` branch is an explicit `BELFEM_ERROR( false, "invalid type" )`, and Maxwell types are rejected up front with "Use MaxwellFactory instead" |

---

### Issue 4 (OPEN): Hardcoded Material Properties in TransientHeatConduction

**Location:** `cl_IWG_TransientHeatConduction.cpp` — search for the
`todo:: replace by cp and lambda` comment; the mass matrix uses `1.0` and the
stiffness `111.0` in place of `material()->cp(T)` / `material()->lambda(T)`.

**Impact:** produces incorrect results for any real material. The Maxwell
thermal path does **not** go through this IWG; this affects the standalone
transient-heat-conduction solver only.

**Fix:** See [Material Property Integration](#material-property-integration).

---

## Development Notes

### Adding New IWG Types

**Steps:**

1. **Create derived class:**
   ```cpp
   class IWG_CustomPhysics : public IWG_Timestep
   {
       void compute_mkf( Element * aElement ) override;  // writes matrices()->M()/K()/f() and sets the flags
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
       // the derived constructor forwards ( IwgType::CustomPhysics, aDimensionality, IwgMode::Iterative ) to IWG_Timestep
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
- **Hairer & Wanner**, "Solving Ordinary Differential Equations II" (BDF stability)
- **Gear (1971)**, "Numerical Initial Value Problems in ODEs"

**BELFEM Papers:**
- **messe2023.txt** - Section 2.7 (Nonlinear iteration strategy)
- **arsenault2023.txt** - Magnetodynamic coupling (read with **arsenault2026.txt**, the erratum correcting the air-domain form)

### Project References

- **Project README**: `../../../../README.md`
- **Claude Instructions**: `../../../../CLAUDE.md`
- **Coding Philosophy**: `../../../../doc/coding_philosophy.md`
- **Documentation Guidelines**: `../../../../doc/documentation_guidelines.md`

---

**Contributors:** Based on code by Christian Messe and Gregory Giard

**Prepared by the BELFEM development team.**
