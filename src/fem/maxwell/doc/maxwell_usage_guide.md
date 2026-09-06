# Maxwell Module Usage Guide {#fem_maxwell_maxwell_usage_guide}

**Module:** src/fem/maxwell
**Version:** 1.1
**Date:** 2026-08-14
**Purpose:** Core usage patterns for BELFEM's electromagnetic Maxwell module

**Revision History:**
- **v1.0** (2026-01-20) - Initial guide covering H-φ formulation patterns
- **v1.1** (2026-08-14) - Corrected against the tree: the input file is
  `input.conf` (the guide previously showed an XML deck that no parser ever
  read), factory/controller/postprocessor snippets now match the real API
  (`hphirun.cpp` is the reference driver), postprocessor types and mesh field
  names verified against `cl_MaxwellPostprocessor.hpp`

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

**The Mistake:**
```cpp
// ❌ WRONG
mesh->create_elements();
IWG_Maxwell* iwg = new IWG_Maxwell(...);
iwg->initialize();  // CRASH! No edge DOFs available
```

**Why it fails:** H-φ formulation requires **edge DOFs** for conductor regions. Without edges created, the DOF manager cannot allocate edge_h DOFs.

**Correct approach:**
```cpp
// ✅ CORRECT
mesh->create_elements();
mesh->create_edges();      // Required for Nédélec edge functions
mesh->create_facets();     // Required for thin-shells (if used)
IWG_Maxwell* iwg = new IWG_Maxwell(...);
iwg->initialize();
```

**Impact:** Segfault or assertion failure during initialization.

---

### 1.2 Using L2 Formulations for Solving

**The Mistake:**
```cpp
// ❌ WRONG - L2 formulations are for postprocessing ONLY
IWG_Maxwell* iwg = new IWG_Maxwell(
    maxwell::Formulation::L2PhiH,  // WRONG!
    ModelDimensionality::ThreeD
);
kernel->solve();  // Will not work correctly
```

**Why it's wrong:**
- `L2PhiH`, `L2PhiB`, `L2EdgeH` are **postprocessing formulations** for field recovery (L2 projection)
- They do NOT implement the full magnetodynamic equations
- Only `HPhi` implements the solving formulation

**Correct approach:**
```cpp
// ✅ CORRECT - Use HPhi for solving
IWG_Maxwell* iwg = new IWG_Maxwell(
    maxwell::Formulation::HPhi,  // CORRECT
    ModelDimensionality::ThreeD
);
kernel->solve();

// THEN use L2 formulations for postprocessing
IWG_MaxwellPostproc* postIWG = new IWG_MaxwellPostproc(
    maxwell::Formulation::L2PhiH  // OK for postprocessing
);
```

**Impact:** Incomplete physics, wrong results.

---

### 1.3 Thin-Shell Without Thickness

**The Mistake:**
```cpp
// ❌ WRONG
block->set_type( DomainType::ThinShell );
// Missing thickness specification!
```

**Why it fails:** Thin-shell formulation scales matrices by thickness δ (see Alves et al. 2022b). Without it, scaling is undefined.

**Correct approach:**
```cpp
// ✅ CORRECT
block->set_type( DomainType::ThinShell );
block->set_thickness( 1e-6 );  // 1 μm thick tape
```

**Typical values:**
- **HTS tapes:** 1-10 μm (1e-6 to 1e-5 m)
- **Metal films:** 0.1-1 μm (1e-7 to 1e-6 m)

---

### 1.4 Loose Convergence Tolerance for HTS

**The Mistake:**
```cpp
// ❌ WRONG for HTS
solver->set_tolerance( 1e-6 );  // Too loose!
```

**Why it's wrong (from Messe et al. 2023, Section 2.7):**
- HTS has extreme nonlinearity (power-law E-J with n=20-40)
- Loose tolerance → checkerboarding in current density
- Oscillations do not damp out

**Correct approach (from Messe et al. 2023):**
```cpp
// ✅ CORRECT for HTS
solver->set_tolerance( 1e-11 );  // Tight tolerance prevents checkerboarding

// Hybrid solver strategy (Messe et al. 2023, Section 2.7):
// 1. Picard iterations (5-10 steps) - robust startup
// 2. Quasi-Newton (optional) - transition
// 3. Full Newton - final convergence with ε < 1e-11
```

**Impact:** Spurious current oscillations, non-physical hot spots, non-convergence.

---

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
| **Conductor** | **H** (magnetic field) | edge_h, face_h | Edges, faces |
| **Air** | **φ** (scalar potential) | phi | Nodes |
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
    ├─ Create DofManager
    ├─ Create Kernel
    ├─ Configure solver
    └─ Initialize boundary conditions
```

### 4.2 Standard Factory Usage

The reference consumer is `src/executables/hphirun.cpp` — a custom driver
follows the same shape:

```cpp
#include "cl_MaxwellFactory.hpp"

// ========================================
// 1. CREATE FACTORY (reads input.conf from the working directory)
// ========================================
MaxwellFactory factory( "input.conf" );

// ========================================
// 2. CREATE MAGNETIC KERNEL
// ========================================
// This single call:
// - Creates mesh
// - Creates cuts (if transport current specified)
// - Creates thin-shells (if HTS tapes specified)
// - Creates IWG_Maxwell
// - Creates DofManager with correct DOF types
// - Allocates DOFs (node phi, edge_h, face_h, lambda)
// - Initializes solver (STRUMPACK/MUMPS/etc.)

auto kernel = factory.create_magnetic_kernel();   // std::shared_ptr< Kernel >

// ========================================
// 3. CREATE CONTROLLER (time stepping, nonlinear loop, saving)
// ========================================
auto control = factory.create_controller();       // std::shared_ptr< Controller >
// Timestep sizes, tolerances etc. come from input.conf — there are no
// setter shortcuts on the controller for them.

// ========================================
// 4. CURRENT BOUNDARY CONDITIONS
// ========================================
// One PhysicalBoundaryCondition per current condition in the deck; the
// driver samples each and hands the vector to the IWG every timestep.
Cell< PhysicalBoundaryCondition * > currentBCs = factory.current_BCs();

// ========================================
// 5. TIME LOOP (the pattern hphirun.cpp uses)
// ========================================
// seed with BELFEM_EPS, not zero — avoids a singular first matrix
Vector< real > tI( currentBCs.size(), BELFEM_EPS );
reinterpret_cast< IWG_Maxwell * >( kernel->dofmgr()->iwg() )->set_currents( tI );

while ( control->time() < control->simulation_time() )
{
    control->initialize_timestep();      // also updates BC values

    // sample each current BC and hand the vector to the IWG — WITHOUT
    // this the whole run solves at zero transport current
    uint tCount = 0;
    for ( PhysicalBoundaryCondition * tBC : currentBCs )
    {
        tI( tCount++ ) = tBC->value();
    }
    reinterpret_cast< IWG_Maxwell * >( kernel->dofmgr()->iwg() )->set_currents( tI );

    control->solve_coupled();               // nonlinear iteration ( certified exit )
    if ( control->reset() ) { continue; }   // timestep was cut — retry

    // ========================================
    // 6. SAVE — postprocessing runs inside the controller; the
    //    postprocessors were wired by the factory from the topology
    //    domain types, and the fields land on the mesh
    // ========================================
    if ( control->save() || control->time() >= control->simulation_time() )
    {
        control->finalize( true );
        control->save( "hphi_results.e-s" );
        control->save_memdump( "memdump.hdf5" );
    }
    else
    {
        control->finalize( false );
    }
}
// kernel and controller are shared_ptr-owned; no manual cleanup
```

### 4.3 Manual Setup (Without Factory)

For custom workflows:

```cpp
#include "cl_IWG_Maxwell.hpp"
#include "cl_FEM_Kernel.hpp"

// ========================================
// 1. LOAD MESH
// ========================================
Mesh* mesh = new Mesh( "mesh.msh" );
mesh->create_edges();    // Required!
mesh->create_facets();   // For thin-shells

// ========================================
// 2. CREATE IWG
// ========================================
IWG_Maxwell* iwg = new IWG_Maxwell(
    maxwell::Formulation::HPhi,
    ModelDimensionality::ThreeD,
    false,  // higher_order
    true    // use_enrichment
);

// ========================================
// 3. CREATE KERNEL & DOF MANAGER
// ========================================
KernelParameters* params = new KernelParameters( mesh );
Kernel* kernel = new Kernel( params );   // Kernel takes KernelParameters*, not Mesh*
DofManager* dofMgr = kernel->create_field( iwg );

// ========================================
// 4. SELECT BLOCKS/SIDESETS
// ========================================
// Conductor blocks
iwg->select_blocks( { 1, 2, 3 } );

// Interface sidesets
iwg->select_sidesets( { 10, 11 } );

// ========================================
// 5. INITIALIZE
// ========================================
dofMgr->initialize();

// ========================================
// 6. SOLVE
// ========================================
// The solver lives on the DofManager, configured from SolverParameters —
// the factory path fills these from the deck's solver { linear { } } section
SolverParameters tParams( SolverType::STRUMPACK );  // cl_SolverParameters.hpp
dofMgr->set_solver( tParams );                      // cl_FEM_DofManager.hpp

// Time stepping and the nonlinear loop belong to the Controller in the
// factory path (§4.2); a fully manual loop drives the DofManager directly:
dofMgr->compute_jacobian_and_rhs();
dofMgr->solve();

// Cleanup
delete kernel;
delete iwg;
delete mesh;
```

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
    Cell<string> Conductor;   // edge_h, face_h (Nédélec)
    Cell<string> Air;         // phi (Lagrange node)
    Cell<string> Ferro;       // phi (Lagrange node)
    Cell<string> Coil;        // (user-defined)
    Cell<string> ThinShell;   // edge_h + phi (from master/slave)
    Cell<string> Cut;         // Interface DOFs for cohomology

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
IWG_Maxwell* iwg = new IWG_Maxwell(
    maxwell::Formulation::HPhi,
    ModelDimensionality::ThreeD,
    true,   // higher_order = TRUE
    false
);

// Now conductor has:
Conductor = { "edge_h", "face_h" }  // + face DOFs

// DOF counts increase:
// - Edge DOFs: ~15000 (same)
// - Face DOFs: ~30000 (2 per face in 3D)
// Total: ~47000 DOFs (3× more)
```

**Use higher-order when:**
- Need better accuracy for complex geometries
- Willing to pay 3-5× computational cost

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
- `_defect`: Includes defect model
- `_piecewise`: Piecewise linearization for Newton solver

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

**Manual cut creation:**

```cpp
#include "cl_Topology.hpp"

// Create topology analyzer
Topology* topo = new Topology( mesh );

// Find cohomology generators
topo->compute_cohomology();

// Create cuts
Cell<mesh::Facet*> cuts = topo->create_cuts();

// Apply to DOF manager
for ( mesh::Facet* cut : cuts )
{
    dofMgr->create_cut_interface( cut );
}
```

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
- `air/ferro/buffer symmetry` is currently realized as a hard `impose_dirichlet(φ=0)`
  (`cl_MaxwellFactory.cpp:626`) — a special case of `B×n=0` that also pins the
  potential. The general `symmetry_phi` weak form (`B×n=0` with φ free to float)
  exists in the IWG (`cl_IWG_Maxwell.cpp:419`) but is presently shadowed by that
  Dirichlet.
- Untagged one-sided air/ferro/conductor boundaries are auto-classified as
  `*AntiSymmetry` by `Topology` (`cl_Topology.cpp:221`), i.e. they default to the
  magnetic wall.
- The "natural = B·n = 0" argument is exact only in the scalar-φ region. A plane
  that cuts a **conductor / H(curl) region**, the **h-φ interface**, or a **cut**
  is not automatically a magnetic wall and may need explicit treatment.

---

## 8. Thin-Shell Formulation

See Section 20 of `/src/fem/kernel/doc/dof_manager_usage_guide.md` for full theory.

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

**For typical HTS simulation (100k DOFs, 100 timesteps):**

| Operation | Time % | Scaling | Bottleneck |
|-----------|--------|---------|------------|
| **Sparse factorization** | 60% | O(N^1.5-2.3) | Fill-in, pivoting |
| **Element assembly** | 25% | O(N_elem × n_dof²) | Material evaluation (HTS E-J) |
| **Postprocessing** | 10% | O(N_elem) | L2 projection |
| **MPI communication** | 5% | O(log P) | Sparse matrix assembly |

### 10.2 Optimization Strategies

**1. Solver selection (from Messe et al. 2023):** the deck picks the solver;
there is no runtime fallback switch:

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

Avoid iterative solvers (GMRES, BiCGSTAB) for HTS work — the E-J nonlinearity
ill-conditions the system; the direct solvers are the validated path.

**2. Material property caching:**

The signature is `rho( T, B, beta )` — **temperature first**. Every parameter is `real`, so a
call written in any other order compiles cleanly and returns a wrong resistivity.

```cpp
// ❌ SLOW - Evaluate material every integration point
for ( uint k = 0; k < num_intpoints; ++k )
{
    real rho = material->rho( T, norm(b), angle );  // Expensive!
    ...
}

// ✅ FAST - Cache if field doesn't vary much within element
real b_avg = average_field_in_element();
real rho = material->rho( T, norm(b_avg), angle );  // Once per element

for ( uint k = 0; k < num_intpoints; ++k )
{
    // Use cached rho
    ...
}
```

**3. Reuse calculators:**

```cpp
// ❌ SLOW - Create calculator per element
for ( Element* elem : block->elements() )
{
    Calculator* calc = new Calculator( elem );  // Allocation overhead!
    maxwell::h_picard( calc, matrices );
    delete calc;
}

// ✅ FAST - Reuse calculator
Calculator* calc = block->calculator();  // Allocated once
for ( Element* elem : block->elements() )
{
    calc->link( elem );  // Just relink
    maxwell::h_picard( calc, matrices );
}
```

**4. MPI partitioning:**

```cpp
// Use METIS for load-balanced partitioning
mesh->partition( comm_size(), PartitionMethod::METIS );

// Minimizes:
// - Edge cuts (MPI communication)
// - Load imbalance
```

---

## 11. Thread Safety and MPI

**Maxwell module is NOT thread-safe:**
- No internal mutexes
- Use MPI for parallelism (distributed memory)

**MPI usage:**

```cpp
#include "cl_Communicator.hpp"

// Initialize MPI (hphirun.cpp does exactly this)
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
material-agnostic. See Section 17 in
`/src/fem/kernel/doc/dof_manager_usage_guide.md`.

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

**Last updated:** 2026-01-20
