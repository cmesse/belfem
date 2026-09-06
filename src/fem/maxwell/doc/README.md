# Maxwell Module Documentation {#fem_maxwell_index}

**Module:** src/fem/maxwell
**Purpose:** Index of documentation for BELFEM's electromagnetic (Maxwell) physics module

---

## Overview

The **Maxwell** module implements electromagnetic field simulation using the **H-Phi mixed formulation** for BELFEM's finite element framework. It provides:

- H-φ magnetodynamic formulation for conductors, ferromagnetics, and air
- Nédélec edge/face elements for H(curl) spaces
- Thin-shell formulation for superconducting tapes and films
- Cohomology cuts for transport current boundary conditions
- Material models: pure metals, alloys, high-temperature superconductors (HTS)
- Static condensation at interfaces (avoiding zero-diagonal blocks)
- Field postprocessing (B, H, J, J/Jc)

---

## Documentation Files

### User Guides

- **[maxwell_usage_guide.md](maxwell_usage_guide.md)** - Comprehensive usage guide for the Maxwell module
- **[maxwell_weak_forms.md](maxwell_weak_forms.md)** - The weak-form derivations behind the formulations (extracted 2026-08-14 from Christian's pre-BELFEM theory notes): fundamental lemma, divergence-theorem corollaries, least squares projection, thermal warm-up example, b-conform (a and a-v, derived but not implemented) and the implemented h-conform weak form
- **[thin_shell_virtual_domains.md](thin_shell_virtual_domains.md)** - What a virtual domain is (thin-shell surface extruded by `input.conf` thicknesses, added to the output `.exo` for visualization, intersecting air by design), and the edge-allocation rule: a thin-shell block gets edges only when it is a **conductor** — ferro/void/non-conducting buffer get none. Governs why non-conducting virtual facets must be skipped by edge-consuming passes (e.g. the periodic `collect_edges`).
- **[thin_shell_facet_orientation.md](thin_shell_facet_orientation.md)** - Why facet master/slave orientation matters for thin shells, the two-stage fix algorithm, and common pitfalls
- **[thinshell_postprocessor_node_sharing.md](thinshell_postprocessor_node_sharing.md)** - Known visualization-only artifact: thin-shell layer blocks share interface nodes, so SPR-recovered nodal fields (e.g. `J/Jc`) appear on both sides of a `hasDuplicates` interface even though the H system is decoupled correctly. Includes a future fix sketch and an explanation of why we're deferring it (thermal module dependency).
- **[buffer_cut_topology.md](buffer_cut_topology.md)** - How the cohomology cut graph must be modified when a phi-formulation buffer layer splits the thin-shell conductor into two halves. Single-cut rerouting approach (alternative to Schnaubelt's two-cut method). Includes orientation rule and CutFactory modification checklist.
- **[contact_impedance_theory.md](contact_impedance_theory.md)** - Derivation of the contact-impedance surface bilinear form for a collapsed rint layer. Distinguishes the measurable contact resistance `rho*h` from the `rho/h` that actually enters the matrix, shows the TSA gives the same conditioning but eliminates DOFs and the Nitsche penalty parameter, and gives the matrix structure for an `h_contact_impedance()` kernel.
- **[coulomb_gauge_penalty_theory.md](coulomb_gauge_penalty_theory.md)** - Whether a Coulomb-gauge penalty can condition the h-formulation, and where. Derives the dimensionally consistent form (`gamma = chi*rho*/mu^2`, K-channel), proves the element-local no-go on TET4/TRI3 and the thin-shell kernels (their curl null space is element-wise constant, cohomology modes included), tabulates the basis divergence per element family, defines the `G` gradient operator and its layout contract, and records what the measured quench-deck condition numbers actually are (flat, and anti-correlated with solver difficulty). Deck key: `nonlinear magnetic { coulomb gauge penalty { chi } }`. **Opt-in** since 2026-09-01 (absent block = off; it was on at `chi = 1e-4` from 2026-08-27 to then); switch on with a positive `chi`.
- **[side_coating_wall_element.md](side_coating_wall_element.md)** - The side-connector campaign reference: physics of HTS tape side conductivity (current escape around the insulating buffer through the plated edge wall), why resolving/enriching failed and the old HEX8TS wrap was unphysical, the collapsed-wall impedance r′, the HEX8TB wall element (4 longitudinal edge dofs, fused/hanging traces), its exact-cuboid edge function (constant Jacobian, E/C operators, stream-function picture, j_t ≡ 0), the h_b/h_n recovery design, and implementation status.
- **[ghost_penalty_stabilization.md](ghost_penalty_stabilization.md)** - What the `h_ghost()` symmetric weighted interior-penalty (Nitsche) kernel solves (weak tangential-H continuity across stacked thin-shell layers), how the regularized harmonic-mean penalty works (Ern, Stephansen & Zunino 2009; Burman & Zunino 2006), how the two stabilization constants (`eta=4.0`, `k_reg=1e-3` Ohm) were chosen, and the `nonlinear magnetic { nitsche ghost penalty { } }` deck block that overrides them (since 2026-08-24; absent block = ghost OFF since 2026-09-01 (opt-in; `eta : 0` is also off, `eta > 0` switches it on)).
  - Common pitfalls and known bugs
  - Mental model (H-φ formulation, DOF structure)
  - Factory patterns and problem setup
  - Material model selection (26 variants)
  - Boundary conditions and cuts
  - Thin-shell conductors and superconductors
  - Field postprocessing and recovery
  - Performance considerations

---

## Quick Reference

### Entry Point Classes

| Class | File | Purpose |
|-------|------|---------|
| **`IWG_Maxwell`** | cl_IWG_Maxwell.{hpp,cpp} | Main electromagnetic IWG (H-φ formulation) |
| **`MaxwellFactory`** | cl_MaxwellFactory.{hpp,cpp} | High-level problem setup orchestrator |
| **`FieldList`** | cl_Maxwell_FieldList.{hpp,cpp} | DOF organization (conductor/air/ferro/cuts) |
| **`MaxwellPostprocessor`** | cl_MaxwellPostprocessor.{hpp,cpp} | Field recovery (B, H, J, J/Jc) |
| **`TMatrix`** | cl_Maxwell_TMatrix.{hpp,cpp} | Nédélec DOF transformations for hanging nodes |
| **`IWG_MaxwellPostproc`** | cl_IWG_MaxwellPostproc.{hpp,cpp} | L2 projection IWG for postprocessing |
| **`MaxwellBoundaryConditionFactory`** | cl_MaxwellBoundaryConditionFactory.{hpp,cpp} | Physical BC setup |

### Matrix Implementations

| File | Purpose | Key Functions |
|------|---------|---------------|
| **mt_maxwell_h.{hpp,cpp}** | **Conductor formulation** (iteration-scheme kernels; material resolved inside) | `h_picard()`, `h_newton_mu0()`, `h_newton_mu()`, `h_ghost()`, `h_side_connector()` |
| **mt_maxwell_phi.{hpp,cpp}** | Air/ferro scalar potential | `phi()`, `phi_ferro()`, `phi_tri3()`, `phi_tet4()` |
| **mt_maxwell_l2_h.{hpp,cpp}** | L2 projection: edge_h → H | `l2_h()`, `l2_ah_2d()` |
| **mt_maxwell_l2_b.{hpp,cpp}** | L2 projection: → B | (implementation) |
| **mt_maxwell_l2_phi.{hpp,cpp}** | L2 projection: φ → H or B | `l2_phi()`, `l2_phi_ferro()` |
| **mt_maxwell_symmetry.{hpp,cpp}** | Symmetry/anti-symmetry BCs | (implementation) |
| **mt_maxwell_background.{hpp,cpp}** | Background field coupling | (implementation) |

---

## Core Enumerations

### Formulation

```cpp
enum class maxwell::Formulation
{
    HPhi,       // H-φ magnetodynamic (SOLVING formulation)
    L2PhiH,     // L2 projection φ → H (postprocessing only)
    L2PhiB,     // L2 projection φ → B (postprocessing only)
    L2EdgeH,    // L2 projection edge_h → H (postprocessing only)
    UNDEFINED
};
```

**Note:** Only `HPhi` is used for solving. The L2* formulations are for field recovery in postprocessing.

---

### MaxwellPostprocessorType

```cpp
enum class MaxwellPostprocessorType
{
    Air,                        // Air/vacuum regions
    Ferro,                      // Ferromagnetic materials
    Conductor,                  // Conductors (no J/Jc output)
    SuperConductor,             // Superconductors (has J/Jc)
    ThinShellConductor,         // Thin-shell conductors
    ThinShellSuperConductor,    // Thin-shell superconductors (has J/Jc)
    SideConnector,              // Edge-coating walls: full h = ht + hb + hn
    UNDEFINED
};
```

**Usage:** Determines which fields are computed during postprocessing:
- **All types:** H, B
- **Conductor/SuperConductor:** + J (current density)
- **SuperConductor only:** + J/Jc (normalized current)

---

## Common Operations

### Creating a Maxwell Problem via Factory

```cpp
#include "cl_MaxwellFactory.hpp"

// Create factory from the input.conf deck in the working directory
MaxwellFactory factory( "input.conf" );

// Create magnetic kernel (mesh, cuts, thin shells, DOFs, IWG, solver)
auto kernel = factory.create_magnetic_kernel();     // shared_ptr< Kernel >

// Create controller (time stepping, nonlinear loop, saving)
auto control = factory.create_controller();         // shared_ptr< Controller >

// Current boundary conditions, one per deck condition
Cell< PhysicalBoundaryCondition * > currentBCs = factory.current_BCs();

// Time loop — see maxwell_usage_guide.md §4.2 and hphirun.cpp for the
// full pattern (initialize_timestep / solve_coupled / finalize + save). Postprocessing runs inside the controller; the
// factory wires the postprocessors from the topology domain types.
```

---

### Manual IWG Creation

```cpp
#include "cl_IWG_Maxwell.hpp"

// Create H-φ IWG
IWG_Maxwell* iwg = new IWG_Maxwell(
    maxwell::Formulation::HPhi,
    ModelDimensionality::ThreeD,
    false,  // higher-order edge functions (default: false)
    true    // hierarchical enrichment at interfaces (default: false;
            // explicitly enabled here)
);

// Link to DOF manager
iwg->set_field( dofManager );
iwg->initialize();

// Select blocks (by material type)
iwg->select_blocks( { 1, 2, 3 } );  // Conductor blocks
iwg->select_sidesets( { 10 } );      // Interface sidesets

// Cleanup
delete iwg;
```

---

### Field Postprocessing

Postprocessors are created by the factory (`create_postprocessors()`), one per
domain kind, with the type derived from the deck's `topology` section; the
controller runs them at every saved timestep. The real constructor, for custom
drivers (`cl_MaxwellPostprocessor.hpp`):

```cpp
#include "cl_MaxwellPostprocessor.hpp"

MaxwellPostprocessor( Kernel * aKernel,
                      const Map< id_t, DomainType > & aBlockTypes,
                      const Map< id_t, string >     & aMaterialMap,
                      const MaxwellPostprocessorType aType,
                      const bool aCreateElementFields = false );

// The interface is run() / initialize() (Postprocessor base class);
// results land as node fields on the mesh:
mesh->field_data( "Bx" );   // flux density components (By, Bz in 3-D)
mesh->field_data( "Hx" );   // magnetic field components
```

---

### T-Matrix for Hanging Nodes

```cpp
#include "cl_Maxwell_TMatrix.hpp"

// Create T-matrix processor
TMatrix tmat( mesh );

// Process a facet (e.g., at a hanging node interface)
for ( mesh::Facet* facet : sideset->facets() )
{
    const Vector<real>& result = tmat.process( facet );

    // result is 12-component vector:
    // - result(0..5):  First edge function contribution
    // - result(6..11): Second edge function contribution

    // Use in hanging DOF elimination (see Alves et al. 2022b, cohomology cuts)
}
```

---

## H Formulation Kernels

**The per-material kernel family this section once tabulated (26 `h_metal`/
`h_hts`/`h_ts_*` variants) no longer exists.** Since the maxwell kernel
collapse, `src/fem/maxwell/matrices/mt_maxwell_h.hpp` exports a small set of
*iteration-scheme* kernels, and material behavior is resolved **inside** them
through the Calculator's material interface (`compute_rho`, `compute_mu`,
`compute_drhodb`, …), not through per-material entry points:

| Function | Role |
|----------|------|
| `h_picard()` | Picard (fixed-point) conductor matrices |
| `h_newton_mu0()` / `h_newton_mu()` | Newton tangent, constant-μ / field-dependent-μ |
| `h_ghost()` | Symmetric weighted interior-penalty tangential-H continuity across stacked thin-shell layers |
| `h_side_connector()` / `h_side_connector_newton()` | Edge-coating wall elements |

All of them except `h_ghost()` also accumulate the element's ohmic dissipation
into the mesh-wide `dotQ` global — a **power in watts**, not an energy; see
§9.3 of the [usage guide](maxwell_usage_guide.md).

The IWG registers the kernel per domain type
(`cl_IWG_Maxwell.cpp`, the `mFunMKF` assignments); whether an element behaves
as metal, alloy or HTS (power-law or piecewise, with or without defects,
temperature-coupled or not) is decided by the material object the deck's
`materials` section assigned to the block — see
[Materials](../../../physics/materials/doc/README.md) and Alves et al. 2022b
for the thin-shell theory.

---

## Field List DOF Structure

The `FieldList` class organizes DOFs by domain type:

| Category | Fields | Mesh Entity |
|----------|--------|-------------|
| **Conductor** | `edge_h`, `face_h` (if higher-order) | Block elements (HTS/metal) |
| **Air** | `phi` | Block elements (air/vacuum) |
| **Ferro** | `phi` | Block elements (ferromagnetic) |
| **Coil** | (user-defined) | Block elements (coils) |
| **ThinShell** | `edge_h`, `face_h`, `lambda` | Sideset facets |
| **Cut** | (interface DOFs) | Sideset facets (cohomology cuts) |
| **Interface** | Mix of conductor + air/ferro DOFs | Sideset facets |

**Non-DOF Fields** (computed, not solved for):
- `MagneticFieldDensity` (B)
- `CurrentDensity` (J)
- `CurrentBC` (cut interface currents)

---

## Quick Reference: Nédélec Elements

**Nédélec elements** are H(curl)-conforming basis functions used for edge DOFs:

### Why Nédélec?

| Property | Standard Lagrange | Nédélec |
|----------|-------------------|---------|
| **Continuity** | Value (C⁰) | Tangential component |
| **DOF location** | Nodes | Edges (and faces for higher-order) |
| **Curl operator** | Requires derivatives | Natural (∇ × H well-defined) |
| **Best for** | Scalar fields (φ) | Vector fields (H) |

### BELFEM Implementation

- **Edge DOFs:** Located on mesh edges, 1 DOF per edge (low-order)
- **Face DOFs:** Located on mesh faces, multiple DOFs per face (higher-order)
- **Edge functions:** Polynomial order = element order + 1 (e.g., linear element → quadratic edge functions)

**See:** `src/fem/interpolation/doc/` for shape function details.

**References:**
- **Nédélec (1980)** - "Mixed finite elements in R³", Numerische Mathematik
- **Monk (2003)** - "Finite Element Methods for Maxwell's Equations"
- **Arsenault et al. 2023** - BELFEM magnetodynamic coupling (Section II)

---

## Workflow: Typical Maxwell Simulation

```
1. **Setup** (MaxwellFactory)
   ├─ Read input file (input.conf)
   ├─ Load mesh (gmsh .msh, or cached .bfm)
   ├─ Create cohomology cuts (if transport current)
   ├─ Create thin-shell structures (if HTS tapes)
   ├─ Assign materials to blocks
   └─ Configure solver (STRUMPACK/MUMPS/...)

2. **Initialization** (create_magnetic_kernel)
   ├─ Create IWG_Maxwell
   ├─ Create DofManager
   ├─ Allocate DOFs (node φ, edge_h, face_h, lambda)
   ├─ Link IWG to DOF manager
   └─ Initialize boundary conditions

3. **Time-Stepping** (controller->run)
   ├─ For each timestep:
   │  ├─ Assemble Jacobian (J) and RHS (f)
   │  │  ├─ Loop over elements
   │  │  ├─ Compute element matrices (M, K, f)
   │  │  ├─ Assemble global system
   │  │  └─ Apply boundary conditions
   │  ├─ Solve linear system (J * Δx = f)
   │  ├─ Update solution (x += Δx)
   │  ├─ Check convergence
   │  └─ Write output
   │
4. **Postprocessing**
   ├─ Create MaxwellPostprocessor
   ├─ Compute derived fields (B, H, J)
   ├─ Export to visualization (VTK/Exodus)
   └─ ( integrating those into global quantities -- loss, force -- is left to the caller;
      the postprocessor writes fields, not integrals )
```

---

## Algorithm Selection Guide

| Problem Type | Deck ingredients |
|--------------|------------------|
| **Eddy currents (AC), metal** | `conductor` block with a metal material |
| **HTS AC loss (bulk)** | `conductor`/`superconductor` block, HTS material (power-law or piecewise) |
| **HTS AC loss (tape)** | `thinshell` + `layers` stack |
| **Transport current** | `current` boundary condition + `homology { algorithm }` (automatic cuts) |
| **Multiphysics (thermo-EM)** | Add a `linear thermal` or `nonlinear thermal` solver section. `belfem` then selects the coupled problem automatically. |

**Notes:**
- The nonlinear machinery defaults to `Picard` + `bdf1` — the validated
  baseline (Messe et al. 2023 §4); `algorithm : Newton` and higher BDF orders
  are deliberate departures (input file reference §4.5).
- The iteration kernel and material model follow from the deck; there is no
  per-problem C++ function to pick.

---

## Maxwell Factory Configuration

Key settings come from the `input.conf` deck — the complete contract is
[the input file reference](@ref doc_input_file_reference); a compact deck
covering the same ground:

```
solver
{
    linear    { library : strumpack ; }   // or mumps, pardiso, umfpack
    nonlinear { tolerance : 1e-11 ; }     // tight, for HTS convergence
    timestep
    {
        initial timestep : 0.1 ms ;
        simulation time  : 100 ms ;
        scheme : bdf1 ;                    // validated baseline
    }
}

materials
{
    ybco { builtin : ybco ; jc : 3e10 ; n : 25 ; }
}

layers : tape { ybco : 1 mum ; }

homology  { algorithm : generalized pellikka ; }

topology
{
    thinshell : tape { sidesets : 1 ; }
    air              { blocks : 1 ; }
    ferro : yoke     { blocks : 2 ; material : iron ; }
}
```

---

## Performance Considerations

### Assembly Performance

| Operation | Complexity | Bottleneck |
|-----------|------------|------------|
| **DOF allocation** | O(N_nodes + N_edges) | Graph traversal |
| **Element loop** | O(N_elem × n_dof²) | Material evaluation |
| **Sparse assembly** | O(N_dof × bandwidth) | Cache locality |
| **Solver factorization** | O(N_dof^1.5 to 2.3) | Fill-in, pivoting |
| **Postprocessing** | O(N_elem × n_int) | L2 projection |

**Optimization strategies:**
- **Reuse calculators** - Don't create per element
- **Cache material properties** - Expensive E-J lookups
- **Use STRUMPACK first** - ~2× faster than MUMPS (Messe et al. 2023)
- **MPI distribution** - Partition by METIS for load balance

---

## Common Pitfalls

### 1. **Missing Edges in Mesh**

```cpp
// ❌ WRONG - Mesh without edges
mesh->create_elements();
iwg->initialize();  // CRASH! No edge DOFs available

// ✅ CORRECT
mesh->create_elements();
mesh->create_edges();   // Required for edge_h DOFs
iwg->initialize();
```

**Why:** H-φ formulation requires edge DOFs for conductors.

---

### 2. **Wrong Formulation for Solving**

```cpp
// ❌ WRONG - L2 formulations are for postprocessing only
IWG_Maxwell* iwg = new IWG_Maxwell( maxwell::Formulation::L2PhiH, ... );
kernel->solve();  // WRONG! L2PhiH is not a solving formulation

// ✅ CORRECT
IWG_Maxwell* iwg = new IWG_Maxwell( maxwell::Formulation::HPhi, ... );
kernel->solve();
```

**Note:** `L2PhiH`, `L2PhiB`, `L2EdgeH` are ONLY for field recovery, not solving.

---

### 3. **Thin-Shell Without Thickness**

```cpp
// ❌ WRONG - Thin-shell without thickness specification
block->set_type( DomainType::ThinShell );
// Missing: block->set_thickness( 1e-6 );

// ✅ CORRECT
block->set_type( DomainType::ThinShell );
block->set_thickness( 1e-6 );  // Required for thin-shell formulation
```

**Why:** Thin-shell formulation scales by thickness (see Alves et al. 2022b).

---

### 4. **Loose Convergence Tolerance**

```cpp
// ❌ WRONG - Tolerance too loose for HTS
solver->set_tolerance( 1e-6 );  // May cause checkerboarding

// ✅ CORRECT (from Messe et al. 2023, Section 2.7)
solver->set_tolerance( 1e-11 );  // Tight tolerance prevents oscillations
```

**Why:** HTS nonlinearity requires tight convergence (Messe et al. 2023).

---

### 5. **Postprocessor Type Mismatch**

```cpp
// ❌ WRONG - ThinShell block but Conductor postprocessor
block->set_type( DomainType::ThinShell );
postproc->select_domain_type( MaxwellPostprocessorType::Conductor );

// ✅ CORRECT
block->set_type( DomainType::ThinShell );
postproc->select_domain_type( MaxwellPostprocessorType::ThinShellSuperConductor );
```

**Impact:** Wrong material model used in field computation.

---

## Literature References

### BELFEM Papers

**Essential reading (in `literature/papers/fem`):**

1. **messe2023.txt** - BELFEM core paper (Messe, SUST 2023)
   - H-φ formulation (Equations 6-8)
   - Static condensation at interfaces (avoiding zero diagonal)
   - Nonlinear solver strategy (Section 2.7, ε < 10⁻¹¹)
   - Solver performance (STRUMPACK vs MUMPS)

2. **arsenault2023.txt** - Magnetodynamic h-φ coupling (Arsenault, IEEE TASC 2023)
   - Faraday's law coupling (Section II, Eq. 5-7)
   - Interface conditions (E^CC = E^C)
   - Recommended formulation

3. **alves2022b.txt** - Thin-shell theory and cohomology cuts (Alves)
   - Thin-shell equations (Appendix A, Eq. 19)
   - Interface conditions ([n × h] = j_s, [n · b] = 0)
   - Why N > 1 for stacked tapes: N = 1 leaves the in-plane field short of full
     penetration and underestimates AC losses at low transport current; the solution
     converges towards h-φ as N grows (§4, figure 12; Table 1 for the cost)
   - Cohomology cut generation (Eq. 10)

4. **dular2021.txt** - Mixed formulation stability (Dular)
   - Inf-sup condition (LBB test)
   - Hierarchical enrichment (avoiding spurious oscillations)

5. **alves2024.txt**, **schnaubelt2023.txt** - Transport current and cuts
   - Thick cuts implementation
   - Boundary condition enforcement

---

### Finite Element Theory

**Electromagnetics:**

- **Monk (2003)**, "Finite Element Methods for Maxwell's Equations"
  - Chapter 4: Edge elements and H(curl) spaces
  - Chapter 5: Mixed formulations

- **Nédélec (1980)**, "Mixed finite elements in R³", Numerische Mathematik
  - Original Nédélec element paper

- **Jin (2014)**, "The Finite Element Method in Electromagnetics" (3rd Ed.)
  - Chapter 6: Edge elements
  - Chapter 9: Time-domain analysis

**General FEM (for context):**

- **Brenner & Scott**, "Mathematical Theory of FEM" - Ch. 8 (inf-sup/LBB)
- **Bathe**, "Finite Element Procedures" - §4.4 (mixed methods), §8.6 (Newton-Raphson)
- **Hughes**, "The Finite Element Method" - Ch. 4 (locking remedies)

---

### Related BELFEM Modules

- **IWG** (`src/fem/iwg/doc/`): Time-stepping, weak forms, BDF methods
- **Kernel** (`src/fem/kernel/doc/`): DOF management, JEDI system, hanging nodes
- **Interpolation** (`src/fem/interpolation/doc/`): Nédélec elements, shape functions
- **Homology** (`src/homology/doc/`): Cohomology cut generation algorithms
- **Materials** (`src/physics/materials/doc/`): HTS properties, power-law E-J
- **Sparse** (`src/sparse/doc/`): Solver interfaces (STRUMPACK, MUMPS, PARDISO)

---

## Development Notes

### Adding New Material Models

To add a new material variant to `mt_maxwell_h.cpp`:

1. **Implement weak form**:
   ```cpp
   namespace belfem { namespace fem { namespace maxwell {
       void h_newmaterial( Calculator* aCalc, TimestepMatrices* aMatrices )
       {
           // Get material
           Material* mat = aCalc->group()->material();

           // Element matrices
           Matrix<real>& M = aMatrices->M();
           Matrix<real>& K = aMatrices->K();
           Vector<real>& f = aMatrices->f();

           // Integration loop
           for ( uint k = 0; k < aCalc->num_intpoints(); ++k )
           {
               const Matrix<real>& C = aCalc->C(k);  // Curl operator
               const Matrix<real>& E = aCalc->E(k);  // Edge basis

               // Material properties
               real rho = mat->rho( /* arguments */ );

               // Accumulate
               K += trans(C) * C * (w(k) * rho * aCalc->dV(k));
               M += trans(E) * E * (w(k) * constant::mu0 * aCalc->dV(k));
           }

           // Set flags
           aMatrices->set_flag( MatrixFlag::M );
           aMatrices->set_flag( MatrixFlag::K );
       }
   }}}
   ```

2. **Update IWG function pointer selection** in `cl_IWG_Maxwell.cpp`

3. **Test with manufactured solution**

---

### Extending Postprocessor

To add a new field to `MaxwellPostprocessor`:

1. **Add field to mesh** in `process()`:
   ```cpp
   mesh->create_field( "MyNewField", EntityType::NODE );
   ```

2. **Compute at each node/element**:
   ```cpp
   for ( Node* node : mesh->nodes() )
   {
       real value = compute_my_field( node );
       node->field("MyNewField") = value;
   }
   ```

3. **Export in output**

---

## See Also

- **Project README**: `../../../../README.md`
- **Claude Instructions**: `../../../../CLAUDE.md`
- **Coding Philosophy**: `../../../../doc/coding_philosophy.md`
- **Documentation Guidelines**: `../../../../doc/documentation_guidelines.md`
- **General Documentation**: `../../../../doc/README.md`
- **Literature Routing**: `../../../../literature/README.md`
