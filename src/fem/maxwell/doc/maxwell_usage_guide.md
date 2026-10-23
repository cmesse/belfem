# Maxwell Module Usage Guide {#fem_maxwell_maxwell_usage_guide}

**Module:** src/fem/maxwell
**Version:** 1.2
**Date:** 2026-08-14
**Purpose:** Core usage patterns for BELFEM's electromagnetic Maxwell module

**Revision History:**
- **v1.0** (2026-01-20) - Initial guide covering H-φ formulation patterns
- **v1.1** (2026-08-14) - Corrected against the tree: the input file is
  `input.conf` (the guide previously showed an XML deck that no parser ever
  read), factory/controller/postprocessor snippets now match the real API
  (`hphirun.cpp` is the reference driver), postprocessor types and mesh field
  names verified against `cl_MaxwellPostprocessor.hpp`
- **v1.2** (2026-10-22) - Corrected against the tree: the reference driver is
  `src/executables/belfem.cpp` (`hphirun` is retired and not built); API,
  ordering and distribution examples match the current code

---

## Table of Contents

1. [Common Pitfalls (Read This First!)](#1-common-pitfalls-read-this-first)
2. [Safe Defaults Quick Start](#2-safe-defaults-quick-start)
3. [Mental Model: H-Phi Mixed Formulation](#3-mental-model-h-phi-mixed-formulation)
4. [Factory Pattern Usage](#4-factory-pattern-usage)
5. [DOF Structure and FieldList](#5-dof-structure-and-fieldlist)
6. [Material Model Selection](#6-material-model-selection)
7. [Boundary Conditions and Cuts](#7-boundary-conditions-and-cuts)
8. [Thin-Shell Formulation](#8-thin-shell-formulation)
9. [Field Postprocessing](#9-field-postprocessing)
10. [Performance Considerations](#10-performance-considerations)
11. [Thread Safety and MPI](#11-thread-safety-and-mpi)
12. [Development Notes](#12-development-notes)
13. [Literature References](#13-literature-references)

---

## 1. Common Pitfalls (Read This First!)

⚠️ **Critical mistakes to avoid:**

### 1.1 Missing Edges in Mesh

**The Mistake:** building the magnetic kernel on a mesh with no edges.

**Why it fails:** the H-φ formulation places **edge DOFs** (`edge_h`) on conductor regions. A mesh read from file contains nodes and elements, but no edges. Edges must be created before the DOF manager can attach `edge_h` to them.

**Correct approach:** create the mesh edges before building the kernel and its field. Faces are needed only for higher-order `face_h` DOFs:

```cpp
// ✅ CORRECT
aMesh->create_edges();     // Nédélec edge functions ( edge_h )
aMesh->create_faces();     // only for higher-order face_h DOFs
```

`Mesh` has no `create_elements()` or `create_facets()`; those belong to `ProtoMesh`. The factory performs this step through `MaxwellFactory::create_edges_and_faces_on_mesh()`. It restricts edge creation to conductor groups and creates faces only when `max_element_order() > 1`.

### 1.2 Using L2 Formulations for Solving

**The Mistake:**
```cpp
// ❌ WRONG - L2 formulations are for postprocessing ONLY
IWG_Maxwell * tIWG = new IWG_Maxwell(
    maxwell::Formulation::L2PhiH,  // WRONG!
    ModelDimensionality::ThreeD );
```

**Why it's wrong:**

- `L2PhiH`, `L2PhiB`, and `L2EdgeH` are **postprocessing formulations** for field recovery through L2 projection.
- They do not implement the full magnetodynamic equations.
- `IWG_Maxwell` sets up DOFs only for `HPhi` (`IWG_Maxwell` constructor). The L2 formulations are handled by `IWG_MaxwellPostproc`.

**Correct approach:**
```cpp
// ✅ CORRECT - Use HPhi for solving
IWG_Maxwell * tIWG = new IWG_Maxwell(
    maxwell::Formulation::HPhi,
    ModelDimensionality::ThreeD );

// L2 formulations belong to the postprocessing IWG; the dimensionality is
// required, only the higher-order flag has a default
IWG_MaxwellPostproc * tPostIWG = new IWG_MaxwellPostproc(
    maxwell::Formulation::L2PhiH,
    ModelDimensionality::ThreeD );
```

There is no `Kernel::solve()`. The `Controller` drives the solve (§4.2).

**Impact:** Incomplete physics, wrong results.

### 1.3 Thin-Shell Without Thickness

**The Mistake:**
```cpp
// ❌ WRONG
aBlock->set_domain_type( DomainType::ThinShell );
// Missing thickness specification!
```

**Why it fails:** the thin-shell formulation scales matrices by thickness δ (see Alves et al. 2022b). A block's thickness defaults to `BELFEM_QUIET_NAN` (`cl_Block.hpp`), so a block without an explicit thickness has no usable value.

**Correct approach:** in the factory path, the deck's `layers` section supplies the thicknesses. `ThinShellFactory` then sets the domain type and thickness for each layer block. A hand-built mesh must set both:

```cpp
// ✅ CORRECT
aBlock->set_domain_type( DomainType::ThinShell );
aBlock->set_thickness( 1e-6 );  // 1 μm thick tape
```

**Typical values:**

- **HTS tapes:** 1-10 μm (1e-6 to 1e-5 m)
- **Metal films:** 0.1-1 μm (1e-7 to 1e-6 m)

### 1.4 Loose Convergence Tolerance for HTS

**The Mistake:** running an HTS deck with the class-default nonlinear tolerance, `1e-6` (`mRelativeEpsilonTarget` in `cl_FEM_Controller.hpp`).

**Why it's wrong (from Messe et al. 2023, Section 2.7):**

- HTS has extreme nonlinearity (power-law E-J with n=20-40)
- Loose tolerance → checkerboarding in current density
- Oscillations do not damp out

**Correct approach (from Messe et al. 2023):** the tolerance is a deck key read by `Controller::set_params()`. The solver and controller have no C++ setter for it:

```
solver
{
    nonlinear
    {
        tolerance : 1e-11 ;     // tight tolerance prevents checkerboarding
    }
}
```

**What the iteration does.** `SolverAlgorithm` has three values: `Direct`, `NewtonRaphson`, and `Picard`. There is no separate quasi-Newton stage between them.

The deck key `algorithm` defaults to `Picard`, and the run stays on Picard. With `algorithm : Newton`, each timestep starts with Picard and switches to Newton-Raphson when the residual drops below `tolerance switch` (default `1e-4`, `mEpsilonSwitch`). The code cites Messe et al. 2023, Eq. 13 for this handoff. The keys and their guidance are in the input file reference, §4.2 and §4.5.

**Impact:** Spurious current oscillations, non-physical hot spots, non-convergence.

### 1.5 Postprocessor Type Mismatch

**The Mistake:** constructing a `MaxwellPostprocessor` by hand with a type
that does not match the domain — e.g. `Conductor` (no J/Jc) for an HTS block.

**Why it's wrong:**
- Postprocessor uses wrong material model
- Missing J/Jc computation for HTS
- Incorrect field recovery formulas

**Correct approach:** don't construct postprocessors manually. The factory's
`create_postprocessors()` step derives one postprocessor per domain from the
`topology` domain types in the deck (`cl_MaxwellFactory.cpp`,
`create_postprocessors`), so a `thinshell` over an HTS layer automatically gets
`ThinShellSuperConductor` with J/Jc. If you must construct one, match the real
signature — `MaxwellPostprocessor( Kernel *, Map< id_t, DomainType >,
Map< id_t, string >, MaxwellPostprocessorType, bool )`
(`cl_MaxwellPostprocessor.hpp`) — and match the type to the domain.

**Type mapping:**

| Block DomainType | Postprocessor Type | Has J/Jc? |
|---|---|---|
| Air | Air | No |
| Ferro | Ferro | No |
| Conductor (metal) | Conductor | No |
| Conductor (HTS) | SuperConductor | Yes |
| ThinShell (metal) | ThinShellConductor | No |
| ThinShell (HTS) | ThinShellSuperConductor | Yes |
| Edge-coating wall (metal) | SideConnector | No |

---

### 1.6 Forgetting Static Condensation at Interfaces

**The Pattern (from Messe et al. 2023):**

BELFEM uses **static condensation** (NOT Lagrange multipliers) for interface coupling.

**Why (Messe et al. 2023):**
- Lagrange multipliers → zeros on diagonal → solver instability
- Static condensation → positive-definite system → better performance

**What you need to know:**
- This is automatic in `MaxwellFactory`
- Do NOT manually add Lagrange multiplier DOFs for interfaces
- Use hierarchical enrichment instead (`use_enrichment=true`)

**If implementing custom interfaces:**
```cpp
// ✅ CORRECT pattern (static condensation, Messe et al. 2023)
// 1. Assemble interface element matrices
// 2. Condense interface DOFs via T^T * K * T
// 3. Add to global system (no separate lambda DOFs)

// ❌ WRONG pattern
// Adding explicit Lagrange multiplier rows (creates zero diagonal blocks)
```

---

### 1.7 Distorted QUAD/HEX Elements in the Magnetic Solve

**The rule: quadrilateral and hexahedral elements in a Maxwell problem must be perfectly
rectangular**, meaning axis-aligned or rigidly rotated rectangles/bricks. Trapezoidal,
sheared, or otherwise distorted quads/hexes silently break the edge-element convergence
theory.

**Why (Monk 2003, §6.1 and §8.2–8.3; Falk et al. 2011; Arnold et al. 2001/2002/2005):**
- The Nédélec convergence theory for hexahedra assumes affine (parallelepiped) elements
- A bi-/trilinear map to a distorted element destroys the mapped basis's completeness.
  Lowest-order elements can lose convergence entirely, not just lose one order
- Nothing asserts. The run assembles and solves, but the fields are wrong or refine at a
  broken rate

**What to do:**
- Hex-meshed conductors: use structured rectangular bricks only. Grading is fine; shear is not
- Geometry that will not brick-mesh: use `TET4`. Simplex elements have no such restriction
- Thermal/scalar problems are not affected. This restriction is specific to edge elements

Full mechanism and literature trail: `src/fem/interpolation/doc/nedelec.md`, §6.6.

---

## 2. Safe Defaults Quick Start

BELFEM is driven by a plain-text **`input.conf`** deck, read by the `belfem`
executable from the working directory. The deck decides whether the run is
magnetic-only or coupled h-ɸ/T. The complete key-by-key contract is
[the input file reference](@ref doc_input_file_reference); the shipped decks
under `examples/` are working starting points. A minimal 3-D thin-shell deck
looks like this:

```
mesh
{
    file : tape.msh ;
    unit : mm ;
}

solver
{
    linear
    {
        library : strumpack ;      // first choice; mumps = robust fallback
    }
    nonlinear
    {
        tolerance : 1e-11 ;        // tight, prevents HTS checkerboarding
    }
    timestep
    {
        initial timestep : 0.1 ms ;
        simulation time  : 100 ms ;
    }
}

materials
{
    ybco
    {
        builtin : ybco ;
        jc : 3e10 ;                // A/m^2
        n  : 25 ;
    }
}

layers : tape
{
    ybco : 1 mum ;
}

homology
{
    algorithm : generalized pellikka ;   // automatic cohomology cuts
}

topology
{
    thinshell : tape { sidesets : 1 ; }
    air       { blocks : 1 ; }
}

boundary conditions
{
    current
    {
        input curves : 1 ; output curves : 2 ;
        type : sine ; amplitude : 100 A ; frequency : 50 Hz ;
    }
}
```

Run it from the deck's directory:

```bash
cd myrun/                                   # contains input.conf and the mesh
../cmake-build-debug/bin/belfem             # serial
mpirun -np 4 .../bin/belfem                 # MPI
```

The defaults this leaves in place — Picard iteration, BDF1 time stepping,
Anderson off — are the validated baseline (Messe et al. 2023, §4); treat every
departure as deliberate (see the input file reference, §4.5).

---

## 3. Mental Model: H-Phi Mixed Formulation

### 3.1 The Core Idea

**Problem:** Electromagnetic field simulation in regions with:
- **Conductors** (HTS tapes, metal coils) - eddy currents, nonlinear E-J
- **Air/vacuum** - no currents, simple permeability
- **Ferromagnetics** - nonlinear B-H curves

**Solution:** Mixed formulation
- **H-formulation** in conductors → captures eddy currents directly (∇ × H = J)
- **φ-formulation** in air/ferro → computationally efficient (∇²φ = 0 in air)

### 3.2 Governing Equations

**Conductors (H-formulation, Messe et al. 2023, Eq. 6-8, Arsenault et al. 2023, Eq. 5-7):**

```
∇ × (ρ ∇ × H) + μ₀ ∂H/∂t = 0
```

Where:
- **H** = magnetic field intensity [A/m] - **PRIMARY VARIABLE** (edge DOFs)
- ρ = electrical resistivity [Ω·m] - material property
- μ₀ = permeability of free space = 4π × 10⁻⁷ H/m

The equation above assumes a conductor with constant permeability μ₀. The code also supports a variable-μ branch. Under Newton-Raphson, `IWG_Maxwell::link_to_group` binds `h_newton_mu0` when the conductor material's `mu` is constant and `h_newton_mu` when it is not. Picard uses `h_picard` in both cases.

**Air/Ferro (φ-formulation):**

```
∇ · (μ ∇φ) = 0
H = -∇φ
```

Where:
- **φ** = magnetic scalar potential [A] - **PRIMARY VARIABLE** (node DOFs)
- μ = magnetic permeability [H/m] - μ₀ in air, μ(H,T) in ferro

**Weak form (conductor):**

```
∫ (ρ ∇ × w) · (∇ × H) dV + ∫ μ₀ w · ∂H/∂t dV = 0
```

Where **w** = edge test functions (Nédélec H(curl) basis)

**Weak form (air/ferro):**

```
∫ (μ ∇w) · (∇φ) dV = 0
```

Where w = nodal test functions (Lagrange C⁰ basis)

### 3.3 DOF Placement

| Domain | Primary Variable | DOF Type | Mesh Entity |
|--------|------------------|----------|-------------|
| **Conductor** | **H** (magnetic field) | edge_h (+ face_h when higher order) | Edges (+ faces) |
| **Air** | **φ** (scalar potential) | phi | Nodes |
| **Buffer** | **φ** (scalar potential) | phi (uses the Air DOF table) | Nodes |
| **Ferro** | **φ** (scalar potential) | phi | Nodes |
| **ThinShell** | **H** + **φ** (coupled) | edge_h + phi (from neighbors) | Facets |

**Key insight (Arsenault et al. 2023):**
- Edge DOFs (H) naturally enforce **tangential continuity** across elements
- Node DOFs (φ) naturally enforce **normal continuity** (H = -∇φ)
- Interface coupling via **static condensation** (Messe et al. 2023)

### 3.4 Why This Formulation? (vs alternatives)

| Formulation | BELFEM Choice | Why? |
|-------------|---------------|------|
| **A-φ** (vector potential) | ❌ NOT used | Requires gauge fixing, more DOFs |
| **T-Ω** (current potential) | ❌ NOT used | Limited to simply-connected regions |
| **H-φ** (BELFEM's choice) | ✅ USED | Natural for conductor/air split, no gauge |

**References:**
- **Messe et al. 2023** - BELFEM core (Equations 6-8)
- **Arsenault et al. 2023** - Magnetodynamic coupling (Section II, Eq. 5-7)
- **Monk (2003)** - "Finite Element Methods for Maxwell's Equations"

---

## 4. Factory Pattern Usage

### 4.1 MaxwellFactory: The Orchestrator

`MaxwellFactory` handles the complete workflow:

```
MaxwellFactory
    │
    ├─ Read input file (input.conf)
    ├─ Load mesh (gmsh .msh, or cached .bfm)
    ├─ Create cohomology cuts (if transport current)
    ├─ Create thin-shell structures
    ├─ Assign materials to blocks
    ├─ Create IWG_Maxwell
    ├─ Create Kernel
    ├─ Create DofManager (Kernel::create_field)
    ├─ Configure solver
    └─ Initialize boundary conditions
```

### 4.2 Standard Factory Usage

The built driver is `src/executables/belfem.cpp`. A custom driver follows the same
shape. The excerpt below is its magnetic-only and fully coupled path:

```cpp
#include "cl_InputFile.hpp"
#include "cl_MaxwellFactory.hpp"
#include "cl_ThermalFactory.hpp"
#include "cl_ElectricalCircuitFactory.hpp"
#include "cl_FEM_Controller.hpp"

const string tFile = "input.conf" ;

// the deck selects the coupled h-ɸ/T problem with a `linear thermal` or
// `nonlinear thermal` section in its solver block ( deck_requests_thermal()
// in belfem.cpp ). Every rank parses the file and takes the same branch:
// the factories are collective
bool tHaveThermal ;
{
    InputFile tInput( tFile );
    tHaveThermal = tInput.section_exists( "solver" )
        && ( tInput.section( "solver" )->section_exists( "linear thermal" )
          || tInput.section( "solver" )->section_exists( "nonlinear thermal" ) );
}

// ========================================
// 1. CREATE FACTORY ( reads input.conf, loads the mesh )
// ========================================
MaxwellFactory tMFactory( tFile );

// the material laws need a reference temperature even without a thermal solve
if ( std::isnan( gTbulk ) )
{
    gTbulk = 77.0 ;
}

// ========================================
// 2. CIRCUIT — before the kernel: the circuit factory adds its terminal
//    boundary conditions to the Maxwell factory's list
// ========================================
electronics::ElectricalCircuitFactory tElFactory( tFile, tMFactory.boundary_conditions() );

// ========================================
// 3. CREATE MAGNETIC KERNEL
// ========================================
// This single call creates the cuts, the thin shells, the IWG_Maxwell,
// the kernel and its DofManager ( phi and edge_h; face_h only for
// higher-order elements ), configures the solver, assigns the materials,
// and wires the postprocessors
auto tKernel = tMFactory.create_magnetic_kernel();   // std::shared_ptr< Kernel >
tKernel->mesh()->field_data( "phi" ).fill( 0.0 );

// one PhysicalBoundaryCondition per current condition in the deck
Cell< PhysicalBoundaryCondition * > tCurrentBCs = tMFactory.current_BCs();

// seed with BELFEM_EPS, not zero — avoids a singular first matrix
Vector< real > tI( tCurrentBCs.size(), BELFEM_EPS );
reinterpret_cast< IWG_Maxwell * >( tKernel->dofmgr()->iwg() )->set_currents( tI );

// ========================================
// 4. CREATE CONTROLLER ( time stepping, nonlinear loop, saving )
// ========================================
// timestep sizes, tolerances etc. come from input.conf — there are no
// setter shortcuts on the controller for them
auto tControl = tMFactory.create_controller();       // std::shared_ptr< Controller >

// ========================================
// 5. THERMAL KERNEL — only when the deck asks for it; the controller is
//    never handed a null kernel
// ========================================
std::unique_ptr< ThermalFactory > tTFactory ;
std::shared_ptr< Kernel >         tKernel2 ;

if ( tHaveThermal )
{
    tTFactory.reset( new ThermalFactory( tFile, tKernel.get() ) );
    tKernel2 = tTFactory->create_thermal_kernel();
    tKernel2->mesh()->field_data( "T" ).fill( gTbulk );
    reinterpret_cast< IWG_Timestep * >( tKernel2->dofmgr()->iwg() )
        ->set_timestepping_method( tControl->euler_method() );
    tControl->set_thermal_kernel( tKernel2.get() );
}

tControl->set_circuit( tElFactory.circuit() );

// resumes from a backup if one exists
tControl->load_memdump( "memdump.hdf5" );

string tExodusFile = sprint( "%s.e-s", tMFactory.label().c_str() );

// ========================================
// 6. TIME LOOP ( magnetic-only, or thermal and fully coupled )
// ========================================
while ( tControl->time() < tControl->simulation_time() )
{
    tControl->initialize_timestep();      // also updates BC values

    // sample each current BC and hand the vector to the IWG — WITHOUT
    // this the whole run solves at zero transport current
    uint tCount = 0;
    for ( PhysicalBoundaryCondition * tBC : tCurrentBCs )
    {
        tI( tCount++ ) = tBC->value();
    }
    reinterpret_cast< IWG_Maxwell * >( tKernel->dofmgr()->iwg() )->set_currents( tI );

    tControl->solve_coupled();               // skips the thermal phase without a thermal kernel
    if ( tControl->reset() ) { continue; }   // timestep was cut — retry

    // postprocessing runs inside the controller; the fields land on the mesh
    if ( tControl->save() || tControl->time() >= tControl->simulation_time() )
    {
        tControl->finalize( true );

        // save_IV must come before the Exodus write: it refreshes the
        // terminal I/U globals that the Exodus file records
        tControl->save_IV( "iv_results.csv" );
        tControl->save( tExodusFile );
        tControl->save_memdump( "memdump.hdf5" );
    }
    else
    {
        tControl->finalize( false );
    }
}
// kernels and controller are smart-pointer-owned; no manual cleanup
```

`belfem` runs this loop when there is no thermal kernel or when `tControl->is_fullycoupled()` is true. A thermal run that is not fully coupled uses a segregated loop instead.

That loop calls `initialize_magnetic()` / `solve_magnetic()`. It then iterates `initialize_thermal()` / `solve_thermal()` until the thermal time catches up, and writes `hphi_results.e-s`. See `belfem.cpp` for the full branch.

### 4.3 Setup Without the Factory

There is no short manual equivalent of `create_magnetic_kernel()`. In addition to the kernel and field, it creates the cuts, thin shells, gauge pins, materials, and postprocessors.

`MaxwellFactory::create_magnetic_kernel()` (`cl_MaxwellFactory.cpp`) is the reference for the required order of operations. A custom setup must preserve the same constraints:

| Constraint | Reason in the code |
|---|---|
| Create edges (and faces for `face_h`) on the mesh before the kernel | §1.1; the factory does it in `create_edges_and_faces_on_mesh()` |
| Select the kernel's blocks in `KernelParameters` before constructing the `Kernel` | the constructor partitions the mesh over the selected blocks; the factory selects air, buffer, conductor, ferro and thin-shell blocks, and excludes coils |
| Pass block and sideset types with `IWG::set_blocks( ids, types )` / `IWG::set_sidesets( ids, types )` **before** `Kernel::create_field()` | `create_field()` calls `IWG::initialize()` immediately, and `IWG_Maxwell::initialize()` builds the per-block DOF tables from those types. `select_blocks()` stores IDs only and assigns no types. Factory: `set_block_types_in_magnetic_equation()` |
| Call `DofManager::set_solver( SolverParameters )` before `DofManager::initialize()` | `initialize()` asserts `"set solver first"`. The factory calls `configure_solver()` right after `create_field()` |
| Decide ownership explicitly | `Kernel` borrows its `KernelParameters` and its mesh unless `claim_parameter_ownership()` / `claim_mesh_ownership()` is called; the factory claims both. `create_field()` does not take ownership of the IWG |

---

## 5. DOF Structure and FieldList

### 5.1 FieldList Organization

The `maxwell::FieldList` class organizes DOFs by domain type:

```cpp
// Located at: src/fem/maxwell/cl_Maxwell_FieldList.hpp

class FieldList
{
public:
    // ===== MAIN DOFS (set by IWG) =====
    Cell<string> Conductor;   // edge_h (+ face_h when higher order, Nédélec)
    Cell<string> Air;         // phi (Lagrange node); Buffer blocks use this table
    Cell<string> Ferro;       // phi (Lagrange node)
    Cell<string> Coil;        // empty: coil blocks carry no DOFs
    Cell<string> ThinShell;   // edge_h + phi (from master/slave)
    Cell<string> Cut;         // not filled by IWG_Maxwell
    Cell<string> Ghost;       // copy of the Conductor list (ghost elements)
    Cell<string> Lambda;      // DOF names starting with "lambda"; IWG_Maxwell creates none

    // ===== INTERFACE DOFS =====
    Cell<string> InterfaceCondAir;   // Conductor-Air
    Cell<string> InterfaceCondFm;    // Conductor-Ferro
    Cell<string> InterfaceFmFm;      // Ferro-Ferro
    Cell<string> InterfaceFmAir;     // Ferro-Air

    // ===== SYMMETRY BCS =====
    Cell<string> SymmetryAir;
    Cell<string> SymmetryFerro;
    Cell<string> SymmetryConductor;
    Cell<string> AntiSymmetryAir;
    Cell<string> AntiSymmetryFerro;
    Cell<string> AntiSymmetryConductor;

    // ===== BOUNDARY BCS =====
    Cell<string> BoundaryAir;
    Cell<string> BoundaryFerro;
    Cell<string> BoundaryConductor;

    // ===== NON-DOF FIELDS (computed) =====
    Cell<string> MagneticFieldDensity;  // B (from H or φ)
    Cell<string> CurrentDensity;        // J (from ∇ × H)
    Cell<string> CurrentBC;             // Cut interface currents
};
```

### 5.2 Typical DOF Assignment Example

**Problem setup:**
- Block 1: HTS conductor
- Block 2: Ferromagnetic yoke
- Block 3: Air region
- Sideset 10: HTS-Air interface
- Sideset 11: Ferro-Air interface

**Resulting DOF structure:**

```cpp
Conductor = { "edge_h" }                     // Block 1 (HTS)
Ferro = { "phi" }                            // Block 2
Air = { "phi" }                              // Block 3

InterfaceCondAir = { "edge_h", "phi" }       // Sideset 10 (HTS-Air)
InterfaceFmAir = { "phi" }                   // Sideset 11 (Ferro-Air)

ThinShell = { "edge_h", "phi" }              // If thin-shell used

// Total DOF types: 2 (edge_h for HTS, phi for air/ferro)
```

**DOF counts (example 10k element mesh):**
- Node DOFs (phi): ~2000 nodes → 2000 DOFs
- Edge DOFs (edge_h): ~15000 edges in conductor → 15000 DOFs
- **Total: ~17000 DOFs**

### 5.3 Higher-Order Edge Elements

```cpp
IWG_Maxwell * tIWG = new IWG_Maxwell(
    maxwell::Formulation::HPhi,
    ModelDimensionality::ThreeD,
    true,     // aHigherOrder
    false );  // aUseEnrichment

// Now the conductor DOF list is { "edge_h", "face_h" }, and the
// multiplicity is 2 per edge and 2 per face ( IWG_Maxwell constructor ).
// For the example above:
// - Edge DOFs: ~30000 ( 2 per edge, twice the low-order count )
// - Face DOFs: ~2 per conductor face
```

**Use higher-order when:**
- Need better accuracy for complex geometries
- Willing to pay a higher cost per solve

**Avoid when:**
- Thin-shells (low-order sufficient, Alves et al. 2022b)
- Large meshes (>100k elements)

---

## 6. Material Model Selection

**There is no per-material kernel to call.** The 26-variant `h_metal`/`h_hts`
family this section once tabulated was collapsed: `mt_maxwell_h.hpp` exports
iteration-scheme kernels only (`h_picard`, `h_newton_mu0`, `h_newton_mu`,
`h_ghost`, `h_side_connector*`), and the material behavior — metal vs. HTS,
power-law vs. piecewise resistivity, defects, temperature coupling — is
resolved inside them through the Calculator's material interface
(`compute_rho`, `compute_mu`, `compute_drhodb`, …).

Selection therefore happens in the **deck**, not in C++:

```
materials
{
    ybco
    {
        builtin : ybco ;
        jc : 3e10 ; n : 25 ;              // power-law E-J
        // resistivity type : piecewise ; // piecewise linearization
        // defect { file : lib.so ; label : MyDefect ; }
    }
    copper { builtin : copper ; RRR : 50 ; }
}
```

Temperature coupling follows from giving the deck a thermal solver section,
which makes `belfem` select the coupled h-ɸ/T problem. See `doc/input_file_reference.md` §5 for every material key.
Defects and the piecewise law are also deck keys: `defect { file ; label ; }` and `resistivity type : piecewise`. They are not kernel variants.

---

## 7. Boundary Conditions and Cuts

### 7.1 Cohomology Cuts (Transport Current)

**Problem:** Imposing transport current I_transport through a conductor loop.

**Traditional approach (wrong):**
- Apply current density J uniformly → assumes current path known
- Fails for complex geometries

**BELFEM approach (Alves et al. 2022b, Alves et al. 2024, Schnaubelt et al. 2023):**
- Use cohomology to find **topologically correct cuts**
- Duplicate nodes along cut
- Impose potential difference: φ_up - φ_down = I_transport

**Factory handles automatically** — the deck names the cut algorithm and the
transport current, and the factory generates the cuts:

```
homology
{
    algorithm : generalized pellikka ;
}

boundary conditions
{
    current
    {
        input curves : 1 ; output curves : 2 ;
        type : ramp ; amplitude : 1000 A ; period : 100 ms ; offset : 0 ms ;
    }
}
```

**There is no public manual cut API.** Cut generation runs inside `create_magnetic_kernel()` through the private `MaxwellFactory::create_cuts()`.

On the master rank, it calls `mesh::Topology::run()`, then `synchronize_maps()`. It reads `homology { algorithm }` and passes the mesh, topology, and thin-shell data to a `mesh::CutFactory`. Its `run()` method produces the cut sidesets.

Worker ranks receive the maps through `synchronize_maps()`. `Topology` itself exposes `run()`, the group and map accessors, and `synchronize_maps()`. It has no cohomology or cut methods.

**See:** `src/homology/doc/` for cohomology theory.

### 7.2 Symmetry Planes: Magnetic Wall vs Flux-Normal

When modeling a fraction (half / quarter) of a magnet, each mirror plane carries
one of two conditions. Which one applies is decided by how the **current** mirrors
across the plane:

| Current across plane | Field on plane | Condition | BC needed |
|----------------------|----------------|-----------|-----------|
| **inverts** (anti-symmetric source) | tangential — flux runs *parallel* | **B·n = 0** ("magnetic wall") | **none** |
| **preserved** (symmetric source) | normal — flux crosses *perpendicular* | **B×n = 0** | must be imposed |

**Key point — the magnetic wall (anti-symmetry) needs no boundary condition.**
In the scalar-φ region the weak form `∫ ∇w · μ∇φ` has `B·n = 0` as its *natural*
boundary condition, so an anti-symmetry plane is obtained for free simply by adding
no term (deactivating the sideset). Only the flux-normal `B×n = 0` symmetry plane is
non-natural and must be actively enforced. That is why the module ships only
`mt_maxwell_symmetry.cpp` (`symmetry_phi`, `h_symmetry`) with **no anti-symmetry
counterpart** — none is required.

**Example — cosθ dipole quadrant.** The current is symmetric across the mid-plane
(x-axis) and anti-symmetric across the pole plane (y-axis):

```
topology {
    air symmetry     { sidesets : <x-axis curves> ; }   // B×n = 0, enforced
    // y-axis is a magnetic wall (B·n = 0): natural, needs no entry.
    // Optionally tag it 'air antisymmetry' for documentation — the group
    // is simply deactivated.
}
```

**Implementation notes:**

- `air/ferro/buffer symmetry` is implemented as a hard `impose_dirichlet(φ=0)` (`MaxwellFactory::create_magnetic_kernel()`, the `AirSymmetry` / `BufferSymmetry` / `FerroSymmetry` case of the sideset loop). This is a special case of `B×n=0` that also pins the potential. The general `symmetry_phi` weak form (`B×n=0` with φ free to float) exists in the IWG (`IWG_Maxwell::link_to_group`, which binds `symmetry_phi_2d` / `symmetry_phi_3d`). However, `SideSet::impose_dirichlet()` retypes the sideset to `DomainType::Dirichlet` with geometry-only activation, so that binding does not run.
- Untagged one-sided air/ferro/conductor boundaries are automatically classified as `*AntiSymmetry` by `Topology::detect_sideset_types()` (`cl_Topology.cpp`). They therefore default to the magnetic wall.
- The “natural = B·n = 0” argument is exact only in the scalar-φ region. A plane that cuts a **conductor / H(curl) region**, the **h-φ interface**, or a **cut** is not automatically a magnetic wall and may need explicit treatment.

---

## 8. Thin-Shell Formulation

See Section 20 of `../../kernel/doc/dof_manager_usage_guide.md` for full theory.

**There is no thin-shell kernel.** A shell block is assembled by the same
`h_picard` / `h_newton_mu0` / `h_newton_mu` kernels as a bulk conductor
(`mt_maxwell_h.cpp`); `Conductor` and `ThinShell` share one `case` in
`IWG_Maxwell::link_to_group`. What makes it a shell is the function-pointer
dispatch bound at construction in `calculator::MaxwellData`
(`cl_FEM_Calculator.cpp`, anchor `tIsThinShell`):

```cpp
mFunH = & MaxwellData::compute_h_ts_edge ;   // h = hn + ht
mFunB = & MaxwellData::compute_b_ts ;        // keeps bt and bn apart
```

The normal part comes from the volume elements on either side through `compute_hn`
(`cl_FEM_Calculator.hpp`): on a φ-region ( air, ferro ) the side's trace is
`-grad phi`, on an h-conductor ( `DomainType::Conductor`, e.g. the solder annulus
of corc_solder ) it is the conductor's own Nédélec trace `E * q` averaged over the
facet rule ( `compute_h_trace` ); the two traces are averaged and the average is
**projected onto the facet normal**:

```cpp
hn = 0.5 * ( hm + hs ) ;
n  = tCalc->normal( k );
hn = dot( hn, n ) * n ;      // the projection is not optional
```

Skipping that last line leaves an `O(|H_t|)` tangential contamination, the same
order as the tape's own field. The routine is evaluated once per element, at
`k == 0`, which is exact for flat linear shells and is guarded by an assert on
`element_is_linear()`.

**Critical insight:** The normal field is **imposed** from the surrounding air/ferro region; the tangential field is **solved** from edge DOFs.

### 8.1 Closed rims and terminals: what is accepted

A 3-D thin shell with a closed outer rim (a closed side curve, in the code) and no terminal on the
rim is accepted as of 2026-10-08. For example, a tape can be fed through copper leads, with the
terminals on bulk-conductor faces. The rim nodes remain single-valued, as required by the free-edge
condition. Three limits apply:

- **The closed rim must be isolated.** A closed rim that shares a node with a terminal curve of any
  shell, or with a side curve of another shell, is refused before any node is duplicated. The error
  message names the node and the curve. Keep terminals and other shells off a closed rim.
- **A closed rim that meets a bulk terminal face only along an edge is not detected.** The terminal
  generator flags the edges of that face, so a closed side curve lying on those edges is picked up
  as part of the terminal. There is no check for this case. Keep the rim off terminal faces.
- **Edge coating is not available on a closed rim.** The coating's anchor search needs the rim's
  first node to lie on a terminal curve. On a closed rim, setup stops in
  `ThinShellFactory::compute_binomial_vectors`. Use a bulk-coated rim or an open tape.

The facet alignment, patch relation, and block-side guard in `MaxwellFactory::fix_facet_masters`
run only in 3-D. In 2-D, the per-block element numbering of the mesh provides a uniform master side
along the tape.

---

## 9. Field Postprocessing

### 9.1 MaxwellPostprocessor Usage

Postprocessors are created by the factory, not by the user:
`MaxwellFactory::create_postprocessors()` pushes one `MaxwellPostprocessor`
per domain kind onto the field's postprocessor list, deriving the type from
the `topology` domain types. The controller runs them each saved timestep
(`run()` / `initialize()` — the class interface has no separate
`process()`-style entry points), and the results land as node fields on the
mesh, which the controller writes into the Exodus results file:

```cpp
// Fields the factory creates for the postprocessors (3-D; z-components
// are omitted in 2-D) — see MaxwellFactory::create_postprocessors():
mesh->field_data( "Hx" );   // magnetic field [A/m]
mesh->field_data( "Hy" );
mesh->field_data( "Hz" );
mesh->field_data( "Bx" );   // flux density [T]
mesh->field_data( "By" );
mesh->field_data( "Bz" );
```

The real constructor signature, for the rare case that a custom driver builds
one directly (`cl_MaxwellPostprocessor.hpp`):

```cpp
MaxwellPostprocessor(
    Kernel * aKernel,
    const Map< id_t, DomainType > & aBlockTypes,
    const Map< id_t, string >     & aMaterialMap,
    const MaxwellPostprocessorType aType,
    const bool aCreateElementFields = false );
```

### 9.2 Postprocessor Type Selection

The types, from `cl_MaxwellPostprocessor.hpp`:

```cpp
enum class MaxwellPostprocessorType
{
    Air,                        // B = -μ₀ ∇φ, H = B/μ₀
    Ferro,                      // B = -μ(H) ∇φ, H = B/μ (nonlinear)
    Conductor,                  // H from edge dofs, no J/Jc
    SuperConductor,             // + J/Jc
    ThinShellConductor,         // thin-shell, no J/Jc
    ThinShellSuperConductor,    // thin-shell, + J/Jc
    SideConnector,              // edge-coating walls: full h = ht + hb + hn
    UNDEFINED
};
```

The factory selects them from the deck's `topology` section; a `thinshell`
whose layer stack contains an HTS material gets `ThinShellSuperConductor`
automatically.

### 9.3 Ohmic Dissipation Global (`dotQ`)

Alongside the node fields, the contributing h-domain kernels accumulate each
element's ohmic dissipation into one mesh-wide global. BELFEM writes it to the
Exodus results file as the global variable `dotQ`:

```
dotQ  =  ∫ rho |j|^2 dV        over the whole h-domain        [W]
```

**`dotQ` is a power, not an energy.** The dot is the point of the name: it is a
rate, in watts, at one instant. An AC loss per cycle is its time integral, and
**nothing in BELFEM takes that integral** — integrate the per-timestep values
from the Exodus file downstream.

**Scope.** One scalar for the whole mesh. There is no per-element, per-block or
per-domain breakdown, so a run cannot attribute loss separately to the tape, the
stabilizer or the former.

**Lifecycle.** Accumulated per integration point, added to the global per
element, then summed across ranks:

| Stage | Where |
|---|---|
| created, value 0 | `MaxwellFactory::create_magnetic_kernel()` — `create_global_variable( "dotQ", 0 )` |
| zeroed before every assembly | `Controller::reset_dotQ()` |
| accumulated per element | `save_dotQ()` in `mt_maxwell_h.hpp` |
| summed over MPI ranks | `Controller::collect_dotQ()` — rank 0 collects and sums; no-op on one rank |

`reset_dotQ()` and `collect_dotQ()` bracket each `compute_jacobian_and_rhs()`.
The value written for a timestep therefore comes from that timestep's **last**
assembly; it is not an average over the nonlinear iterations.

**Which kernels contribute.** Five of the six h-kernels:

| Kernel | Contributes | Why |
|---|---|---|
| `h_picard()` | yes | volumetric `rho \|j\|^2` |
| `h_newton_mu0()` / `h_newton_mu()` | yes | same term in the Newton branches |
| `h_side_connector()` / `h_side_connector_newton()` | yes | edge-coating wall elements |
| `h_ghost()` | **no** | its `rho` enters a Nitsche stabilization coefficient, not a dissipation term — a penalty is not a loss |

The phi-domain kernels (`mt_maxwell_phi.cpp`) contribute nothing, correctly:
that domain is non-conducting.

> **Note — the resistivity clamp, and why it does not bite.**
> The kernels obtain `rho` through `MaxwellData::compute_rho()`, which clamps it
> to `[ gRhoMin, gRhoMax ]` and supplies the *same* clamped value to every
> consumer (stiffness, Joule source, element mean). The defaults are `0` and
> `1e10` Ohm*m, set in `Communicator::set_globals()`, whose own comment calls
> them a no-op — and the material layer guarantees it. All three HTS laws put
> the power-law channel **in parallel** with the normal-state channel:
> `rho_powerlaw()` and `rho_riva()` return the parallel combination explicitly,
> and `rho_piecewise()` returns the unbounded power-law branch only below its
> flux-flow knee, `rhon` above it. So `rho <= rhon` always, however far
> over-critical an iterate drives `|j|`. A deep Newton excursion cannot reach
> `gRhoMax` — the parallel combination forbids it.
>
> This is deliberate, not a coincidence of magnitudes. The riva channel is
> evaluated in log10 space behind an explicit overflow early-out, and the source
> says why: *"past this cap 1/ρPL vanishes to machine precision against any
> physical ρn — the caller takes the fully-normal branch"* (`powerlaws.hpp`,
> grep that sentence). Past the cap the helper returns `false` and `rho_riva()`
> returns `rhon` outright. `rho_powerlaw()` reaches the same place by arithmetic:
> if `rhoPL` overflows, `1/(1/rhon + 1/rhoPL)` tends to `rhon`.
>
> The upper bound can therefore only fire if a material's own normal-state
> `rho( T )` exceeds `1e10` Ohm*m — no shipped material does, though a user
> plugin material could — or if an executable narrows the window, which none in
> this tree does; a deck cannot touch it at all. The lower bound sits at zero,
> so it can only trip on a *negative* resistivity, which is a material defect
> worth surfacing rather than a guard doing its job.
>
> The consequence for `dotQ` is a real but latent one: where the clamp does
> fire, dissipation is silently capped and `dotQ` becomes a lower bound. The
> calculator records that it happened (`rho_clamped()`), but the flag is not
> propagated to `dotQ`, so a reported value carries no indication either way.
> With stock materials and stock bounds, treat `dotQ` as unclamped.

---

## 10. Performance Considerations

### 10.1 Computational Cost Breakdown

No measured cost breakdown for the Maxwell module is published here. The shares of factorization, assembly, and postprocessing depend on the deck, so measure them on the deck in question.

The controller's timestep footer reports iteration time and, where measured, postprocessing and eigenvalue-analysis times.

### 10.2 Optimization Strategies

**1. Solver selection (from Messe et al. 2023):** the deck selects the solver. There is no runtime fallback switch:

```
solver
{
    linear
    {
        library : strumpack ;   // FIRST CHOICE — ~2× faster than MUMPS
        // library : mumps ;    // FALLBACK — most robust; switch the deck
        //                      // and rerun if STRUMPACK struggles
    }
}
```

Avoid iterative solvers (GMRES, BiCGSTAB) for HTS work. The E-J nonlinearity ill-conditions the system, so direct solvers are the validated path.

**2. Material evaluation happens inside the calculator.** The kernels never call a material directly. They obtain `rho` through `MaxwellData::compute_rho()`, which is bound per group to one of these paths (`cl_FEM_Calculator.hpp`):

- **Metal conductors:** `compute_rho_metal()` calls `Material::rho( T, |b|, beta )`, the normal-state resistivity with magnetoresistance (Kohler's rule in `Metal`, or its precomputed table when the material was built with tables). The signature is **temperature first**. Every parameter is `real`, so a call written in another order compiles cleanly but returns an incorrect resistivity.
- **HTS:** the `compute_rho_powerlaw_*`, `compute_rho_piecewise_*` and `compute_rho_riva_*` families are selected by the deck's `resistivity type`. They evaluate the E-J law (`rho_powerlaw()` and its siblings), which puts the power-law channel in parallel with the normal-state resistivity; see `src/physics/materials/doc/resistivity_laws.md`.

`Material::rho( T, B, beta )` is **not** the E-J law. On a superconductor such as `YBCO` (a `Metal` subclass), it returns the Kohler normal-state value, not the resistivity assembled by the kernels. `MaxwellData` evaluates these quantities at each integration point and caches intermediate values (`b`, `|j|`, `beta`) through its `is_current()` flags.

**3. Reuse calculators:** a `Calculator` is created once per group (its constructors take a `Group *`), and the element loop relinks it. This is what `IWG_Maxwell::compute_mkf()` does:

```cpp
// per element: relink the group's calculator, then call the bound kernel
mGroup->calculator()->link( aElement );
mTimeStepMatrices->reset() ;
( *mFunMKF )( mGroup->calculator(), mTimeStepMatrices );
```

**4. MPI partitioning:** partitioning is automatic. Nothing in `src/` calls `Mesh::partition()`.

In the factory path, rank 0 reads the mesh (`new Mesh( path, 0, false, false )`), broadcasts the mesh dimension, and each worker creates an empty `new Mesh( dim, 0 )`. The `Kernel` constructor then partitions on the master (`Kernel::partition_mesh()`, which calls `graph::metis_partition`) and runs the `mesh::Distributor` on every rank. See "Mesh Distribution (MPI)" in the mesh usage guide (`src/mesh/doc/mesh_usage_guide.md`).

---

## 11. Thread Safety and MPI

**Maxwell module is NOT thread-safe:**
- No internal mutexes
- Use MPI for parallelism (distributed memory)

**MPI usage:**

```cpp
#include "cl_Communicator.hpp"

// Initialize MPI (belfem.cpp does exactly this)
gComm.init( argc, argv );

// Factory is MPI-aware: the mesh is partitioned during kernel creation,
// DOFs are distributed across ranks, the solver uses distributed matrices
MaxwellFactory factory( "input.conf" );
auto kernel  = factory.create_magnetic_kernel();
auto control = factory.create_controller();

// Each rank solves its partition inside the §4.2 time loop; saving
// collects on rank 0 through the controller's save calls

return gComm.finalize();
```

Run with `mpirun -np <ranks> belfem` — Open MPI only (see
[MPI support](@ref doc_mpi_support)).

---

## 12. Development Notes

### 12.1 Adding New Material Behavior

New material behavior goes into the material classes
(`src/physics/materials/`) and the Calculator's material interface — not into
new kernels; the iteration-scheme kernels in `mt_maxwell_h.hpp` are
material-agnostic. The material side is documented in
`src/physics/materials/doc/materials_usage_guide.md`, and the HTS resistivity laws in
`src/physics/materials/doc/resistivity_laws.md`.

### 12.2 Testing New Implementations

**Use manufactured solutions:**

```cpp
// 1. Choose analytical solution
// Example: H = [y, -x, 0] (uniform rotation)

// 2. Compute source term
// ∇ × (ρ ∇ × H) + μ₀ ∂H/∂t = f_source

// 3. Implement in code as RHS forcing

// 4. Run simulation

// 5. Compare: ||H_fem - H_exact|| / ||H_exact||
//    Should converge as O(h^p) where p = element order
```

---

## 13. Literature References

### 13.1 BELFEM Papers (Essential)

**Primary references** (`literature/papers/fem`):

1. **Messe et al. 2023** - BELFEM core (SUST)
   - H-φ formulation (Eq. 6-8)
   - Static condensation (avoiding zero diagonal)
   - Solver strategy (Section 2.7, ε < 10⁻¹¹)

2. **Arsenault et al. 2023** - Magnetodynamic coupling (IEEE TASC)
   - Interface conditions (Eq. 5-7)
   - Faraday's law coupling

3. **Alves et al. 2022b** - Thin-shell theory
   - Thin-shell equations (Appendix A, Eq. 19)
   - Why N > 1 for stacked tapes: at N = 1 the in-plane field does not fully
     penetrate the tape, so the thin-shell model underestimates AC losses at low
     transport current; raising N converges onto the h-φ reference (§4, figure 12)

4. **Alves et al. 2024**, **Schnaubelt et al. 2023** - Transport current and cuts

### 13.2 Electromagnetics FEM

- **Monk (2003)** - "Finite Element Methods for Maxwell's Equations"
- **Jin (2014)** - "The Finite Element Method in Electromagnetics"
- **Nédélec (1980)** - Original edge element paper
- **Arnold et al. (2001, 2002, 2005)**, **Falk et al. (2011)** - Distorted QUAD/HEX
  elements and the loss of edge-element convergence (see §1.7)

### 13.3 General FEM

- **Bathe** - "Finite Element Procedures" (mixed methods, Newton)
- **Brenner** - "Mathematical Theory of FEM" (inf-sup)

---

## See Also

- **Module README**: [Maxwell Module README](@ref fem_maxwell_index) (quick reference)
- **H-Phi Theory**: `maxwell_weak_forms.md` (detailed equations)
- **Kernel Guide**: `../../kernel/doc/dof_manager_usage_guide.md` (Section 19-21 for Maxwell specifics)
- **Project Documentation**: `../../../../doc/README.md`

**Last updated:** 2026-10-22
