# IWG Module Documentation {#fem_iwg_index}

**Module:** src/fem/iwg
**Purpose:** Index of documentation for BELFEM's Integral Weak Form (IWG) physics module

---

## Overview

The **IWG (Integral Weak form of Governing equations)** module implements the physics layer for BELFEM's finite element framework. It provides:

- Various PDE formulations (Poisson, heat conduction, elasticity, Maxwell electromagnetics)
- Time-stepping schemes (BDF1-5 live; Crank-Nicolson and Galerkin parse but hard-error — disabled)
- Nonlinear solvers (Newton-Raphson, Picard)
- Element-level matrix assembly (mass, damping, stiffness, force)
- Jacobian linearization for nonlinear problems

---

## Documentation Files

### User Guides

- **[iwg_usage_guide.md](iwg_usage_guide.md)** - Comprehensive usage guide for the IWG module
  - Common pitfalls and performance traps
  - Mental model (The Four Pillars of FEM)
  - Time-stepping algorithms and BDF methods
  - Factory patterns and IWG creation
  - Nonlinear solver workflows
  - Weak form assembly examples
  - Known issues and critical bugs

---

## Quick Reference

### Entry Point Classes

| Class | File | Purpose |
|-------|------|---------|
| **`IWG`** | cl_IWG.{hpp,cpp} | Abstract base class for all physics |
| **`IWG_Timestep`** | cl_IWG_Timestep.{hpp,cpp} | Base class for transient problems |
| **`TimestepMatrices`** | cl_TimestepMatrices.{hpp,cpp} | **Main workhorse**: element matrix container |
| **`IwgFactory`** | cl_IwgFactory.{hpp,cpp} | Factory for creating IWG instances |

### Key IWG Implementations

| IWG Type | Class | Purpose |
|----------|-------|---------|
| **Poisson** | IWG_Poisson | Scalar Poisson equation (∇²φ = f) |
| **StaticHeatConduction** | IWG_StaticHeatConduction | Steady-state thermal diffusion |
| **TransientHeatConduction** | IWG_TransientHeatConduction | Transient thermal diffusion |
| **Maxwell** | (Separate module) | Electromagnetics (h-φ formulation) |

---

## Core Enumerations

### IwgType

```cpp
enum class IwgType
{
    Poisson,                    // Scalar Poisson equation
    StaticHeatConduction,       // Steady-state thermal
    TransientHeatConduction,    // Transient thermal
    Gradient2D,                 // Gradient computation (2D)
    Gradient3D,                 // Gradient computation (3D)
    SurfaceGradient,            // Surface gradient
    PlaneStress,                // Plane stress (2D elasticity)
    LinearElasticity,           // Linear elasticity (3D)
    Maxwell,                    // Electromagnetics
    MaxwellThermal,             // Coupled electromagnetics + thermal
    UNDEFINED
};
```

### ModelDimensionality

```cpp
enum class ModelDimensionality
{
    TwoD,       // 2D Cartesian
    AxSymmX,    // Axisymmetric about X-axis
    AxSymmY,    // Axisymmetric about Y-axis
    ThreeD,     // 3D Cartesian
    UNDEFINED
};
```

### SolverAlgorithm

```cpp
enum class SolverAlgorithm
{
    Direct,         // Linear solve (no iteration)
    NewtonRaphson,  // Full Newton-Raphson
    Picard,         // Picard iteration (linearized)
    UNDEFINED
};
```

---

## Common Operations

### Creating an IWG via Factory

```cpp
#include "cl_IwgFactory.hpp"

IwgFactory factory( mesh );
IWG* iwg = factory.create_iwg(
    IwgType::TransientHeatConduction,
    ModelDimensionality::ThreeD );   // two parameters; there is no mode argument

// Select blocks (material domains)
iwg->select_blocks( { 1, 2, 3 } );

// Set time-stepping method
iwg->set_timestepping_method( EulerMethod::BackwardDifference2 );

// Initialize DOF manager link
iwg->set_field( dofManager );
iwg->initialize();

// Cleanup
delete iwg;
```

---

### Element-Level Assembly

```cpp
// IWG populates element matrices
TimestepMatrices* matrices = iwg->matrices();
iwg->compute_mkf( element );  // Fills M, K, f and sets flags

// Time-stepper assembles Jacobian
real dt = 0.01;  // Time step
matrices->assemble_J( dt );  // J = M/dt + K (for BDF1)

// Access assembled matrices
const Matrix<real>& J = matrices->J();
const Vector<real>& f = matrices->f();
```

---

## TimestepMatrices Workflow

The `TimestepMatrices` class serves as the interface between physics (IWG) and time-stepping:

**Workflow:**
1. Time-stepper calls `matrices()->reset_flags()`
2. IWG populates element matrices and sets flags:
   - `M()` - Mass matrix
   - `K()` - Stiffness matrix
   - `f()` - Force vector
3. Time-stepper assembles J via `assemble_J(dt)`
4. Global assembler accumulates into sparse system

**Matrix flags:**
```cpp
enum class MatrixFlag
{
    M,               // Mass matrix populated
    D,               // Damping matrix (reserved for 2nd-order)
    K,               // Stiffness matrix populated
    F,               // Force vector populated
    dMdX_times_x,    // Mass derivative contracted with x
    dMdX_times_h,    // Mass derivative contracted with history
    dKdX_times_x,    // Stiffness derivative contracted with x
    dFdX             // Force derivative (Jacobian ∂f/∂x)
};
```

---

## Time-Stepping Methods

| EulerMethod | Order | Stability | Formula |
|-------------|-------|-----------|---------|
| **Static** | 0 | - | No time derivative |
| **ForwardExplicit** | 1 | Conditionally stable | `M * (x - x₀)/h = f` |
| **BackwardDifference1** (BDF1) | 1 | A-stable | `M * (x - x₀)/h + K * x = f` — the validated baseline |
| **BackwardDifference2** (BDF2) | 2 | A-stable | Variable step BDF2 |
| **BackwardDifference3-5** (BDF3-5) | 3-5 | A(α)-stable | Variable step BDF |
| ~~CrankNicolson~~ | 2 | — | **Disabled**: parses, then `BELFEM_ERROR` in `cl_IWG_Timestep.cpp` (the Newton correction is exact for the BDF family only) |
| ~~Galerkin~~ | 1 | — | **Disabled**: same hard error as Crank-Nicolson |

**Auto-coefficient computation:** All BDF methods support variable time-steps. Coefficients α, β are computed automatically in `compute_bdf_coeffs_N()`.

---

## Algorithm Selection Guide

| Problem Type | Recommended IWG | Time Method | Solver |
|--------------|-----------------|-------------|--------|
| **Steady diffusion** | StaticHeatConduction | Static | Direct |
| **Transient diffusion** | TransientHeatConduction | BDF2 | Direct |
| **Nonlinear steady** | Custom IWG | Static | NewtonRaphson |
| **Nonlinear transient** | Custom IWG | BDF2 | NewtonRaphson |
| **Electromagnetics** | Maxwell | BDF2 | NewtonRaphson |

---

## Quick Reference: Nfunction and Bfunction

**Nfunction (Shape function dimensionality):**
```cpp
enum class Nfunction
{
    Scalar,     // Scalar field (N: 1 × num_bases)
    Vec2d,      // 2D vector field (N: 2 × 2*num_bases)
    Vec3d,      // 3D vector field (N: 3 × 3*num_bases)
    UNDEFINED
};
```

**Bfunction (Derivative operator type):**
```cpp
enum class Bfunction
{
    Gradient,      // B = ∇N (for diffusion: -∇·(κ∇T))
    Planestress,   // B for plane stress (2D elasticity)
    Voigt,         // B for 3D elasticity (Voigt notation)
    UNDEFINED
};
```

---

## Source Code

**Module location:** `../`

**Key source files:**

- **Base classes**: `cl_IWG.{hpp,cpp}`, `cl_IWG_Timestep.{hpp,cpp}`
- **Timestep container**: `cl_TimestepMatrices.{hpp,cpp}`
- **Factory**: `cl_IwgFactory.{hpp,cpp}`
- **Poisson**: `cl_IWG_Poisson.{hpp,cpp}`
- **Heat conduction**: `cl_IWG_StaticHeatConduction.{hpp,cpp}`, `cl_IWG_TransientHeatConduction.{hpp,cpp}`
- **Enums**: `en_IWGs.hpp`, `en_IWG_SideSetDofLinkMode.hpp`

---

## External References

### Finite Element Theory

- **Zienkiewicz & Taylor**, "The Finite Element Method" Vol. 1-2 (Weak forms, time integration)
- **Hughes**, "The Finite Element Method" (Transient analysis, Chapter 9)
- **Bathe**, "Finite Element Procedures", Ch. 9 (Time integration), Ch. 8 (Nonlinear analysis)

### Time Integration

- **Hairer & Wanner**, "Solving Ordinary Differential Equations II" (BDF stability)
- **Gear (1971)**, "Numerical Initial Value Problems in Ordinary Differential Equations" (BDF methods)

### BELFEM Papers

- **messe2023.txt** - BELFEM h-φ formulation (Section 2.7: nonlinear iteration)
- **arsenault2023.txt** - Magnetodynamic coupling

### Related BELFEM Modules

- **Interpolation** (`src/fem/interpolation/`): Shape functions, integration points
- **DofManager** (`src/fem/kernel/`): DOF management, global assembly
- **Materials** (`src/physics/materials/`): Material properties (cp, λ, μ)
- **Mesh** (`src/mesh/`): Element definitions, block/sideset structure

---

## Development Notes

### Adding New IWG Types

When implementing new physics:

1. **Derive from IWG or IWG_Timestep**:
   ```cpp
   class IWG_CustomPhysics : public IWG_Timestep
   {
       void compute_mkf( Element * aElement ) override;
       // write into matrices()->M(), matrices()->K(), matrices()->f() and set the flags
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

4. **Implement weak form** in `compute_mkf()`:
   - Populate `aM`, `aK`, `aF`
   - Set flags: `matrices()->set_flag( MatrixFlag::M )`

5. **Test**:
   - Manufactured solution with known analytical result
   - Mesh convergence study (h-refinement)
   - Time-step convergence (for transient)

---

## Common Pitfalls

1. **Creating IWG inside element loops** → Cache IWG per block type
2. **Forgetting to reset flags** → Call `matrices()->reset_flags()` before `compute_mkf()`
3. **Selecting Crank-Nicolson or Galerkin** → disabled; raises `BELFEM_ERROR` (`cl_IWG_Timestep.cpp:88-94`). Use a BDF scheme
4. **Hardcoded material properties** → Use `material()->cp(T)`, not constants
5. **Incomplete factory switch** → Only 3/11 IWG types currently handled
6. **Assuming `is_maxwell()` misses `MaxwellThermal`** → it does not; the check covers both `Maxwell` and `MaxwellThermal` (`en_IWGs.hpp:99-102`)

---

## See Also

- **Project README**: `../../../../README.md`
- **Claude Instructions**: `../../../../CLAUDE.md`
- **Coding Philosophy**: `../../../../doc/coding_philosophy.md`
- **Documentation Guidelines**: `../../../../doc/documentation_guidelines.md`
- **General Documentation**: `../../../../doc/README.md`
