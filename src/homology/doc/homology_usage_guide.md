# Homology Module Usage Guide {#homology_homology_usage_guide}

**Date:** 2026-01-16
**Module:** src/homology
**Purpose:** Comprehensive usage guide for BELFEM's computational topology module

**Revision History:**

| Date | Changes |
|------|---------|
| 2026-01-16 | Initial documentation |
| 2026-01-16 | Applied external review feedback: Added Critical Constraints section (mesh mutation, MPI, memory ownership, thread safety, topological vs FEM correctness, mesh requirements, scratch-field usage); Added Quick-Start section; Added comprehensive Glossary; Moved Common Pitfalls to top with 10 concrete examples |
| 2026-08-31 | Marked the manifold-filtering description as a design of record with no implementation: `manifold_filter_3d()` and `check_surface()` exist nowhere in the tree. Added a design note under "Manifold Cleanup"; corrected the glossary, the Internal Workflow step, and the troubleshooting entry that told readers to tune thresholds that do not exist. Named the shipped path, `Cohomology::clean_spfa()` → `remove_cut_pockets()`. Content status remains Gregory Giard's to confirm |

---

## Module Contracts

Before diving into usage details, understand these fundamental design contracts:

1. **Mesh Dependency**: All topology operations require a valid, connected Mesh object
2. **Integer Arithmetic**: Homology computations use exact integer arithmetic (no floating-point errors)
3. **Algorithm Independence**: Results are topologically equivalent across algorithms; no performance comparison is recorded
4. **Manifold Requirement**: Thin cuts must form valid 2-manifolds for FEM assembly
5. **MPI Awareness**: Cohomology computation is serial: it runs on the master rank, which holds the whole, undistributed mesh. `CutFactory::run()` is nevertheless collective

---

## Critical Constraints

**⚠️ READ THIS SECTION BEFORE USING THE HOMOLOGY MODULE ⚠️**

### Mesh Mutation Contract

The homology module **directly modifies mesh topology** in invasive, irreversible ways:

**Operations performed:**
- **Duplicates nodes** on cut surfaces (creates new Node objects with new IDs)
- **Relinks elements** to point to duplicated nodes (modifies Element connectivity arrays)
- **Creates new SideSets** for cut surfaces (adds to mesh's sideset collection)
- **Creates abstract nodes** not belonging to any element (new DOF carriers)
- **Calls mesh->unfinalize() and mesh->finalize()** multiple times during processing
- **Renumbers the nodes** with reverse Cuthill-McKee (`symrcm`), which reorders `Mesh::nodes()`
  and rewrites every node index, then updates the node, edge and face indices
- **Writes the node field `"phi"`**: an existing field of that name is overwritten with the new
  node indices; otherwise the field is created
- **Removes the edges and faces** at the end of `run()` (`reset_edges()`, `reset_faces()`), so
  code that needs them must create them again

**Consequences:**
- **NOT REVERSIBLE**: Cannot undo node duplication or element relinking
- **Invalidates indices**: Any code relying on original node indices will break
- **Invalidates fields**: Custom mesh fields indexed by original nodes become inconsistent
- **Changes mesh size**: Node count, element connectivity, and sideset count all increase

**Execution order requirement:** run the cut factory before anything that holds node
references or node indices — in particular before the FEM kernel and its `DofManager` are
built. `MaxwellFactory::create_magnetic_kernel()` follows this order: it calls `create_cuts()`
first, hands `Mesh::abstract_nodes()` to the equation (`IWG::set_abstract_nodes`), and builds
the kernel afterwards.

**Cannot run twice:**
```cpp
// WRONG
tFactory.run();  // First run succeeds
tFactory.run();  // ❌ UNDEFINED BEHAVIOR - mesh already mutated
```

---

### MPI Serial Limitation

**CRITICAL:** Cohomology computation is **serial-only** and must run on an **undistributed mesh**.

**What breaks if you ignore this:**
- Running on partitioned mesh → **incorrect simplicial complexes** (missing connectivity across ranks)
- Running on distributed mesh → **incomplete generators** (topological holes span process boundaries)
- Running on multiple ranks → **rank-dependent results** (nondeterministic cut generation)

**Correct MPI workflow.** This follows `MaxwellFactory` (`create_cuts()`,
`create_cuts_sub_master()`, `create_cuts_sub_slave()`). Every call below runs on **every**
rank. Several calls are collective; calling them only on rank 0 deadlocks:

- the file-reading `Mesh( aPath )` constructor is collective: the master reads the file, and
  every other rank gets an empty mesh that knows the dimension count;
- `Topology::run()` works on rank 0 only and does nothing elsewhere;
- `Topology::synchronize_maps()` is collective. It broadcasts the domain-type map, which every
  rank needs before it constructs a `CutFactory`, because the constructor reads
  `groups( DomainType::ThinShell )`;
- `CutFactory::run()` does its work on the master rank but calls `comm_barrier()`, so every
  rank must call it.

```cpp
// every rank: the master reads the file, the others get an empty mesh
Mesh * tMesh = new Mesh( "problem.exo" );

// master only: block domain types decide what is phi and non-phi
if ( comm_rank() == 0 )
{
    tMesh->block( 1 )->set_domain_type( DomainType::Air );
    tMesh->block( 2 )->set_domain_type( DomainType::Conductor );
}

mesh::Topology tTopology( tMesh );
tTopology.run();                // works on rank 0, no-op elsewhere
tTopology.synchronize_maps();   // collective

Cell< Protoshell * > tProtoshells ;   // empty: no thin shells

mesh::CutFactory tFactory( tMesh, &tTopology, tProtoshells,
                           mesh::CutAlgorithm::PellikkaGeneralized,
                           false );
tFactory.set_periodicity( nullptr );  // optional: nullptr is the default
tFactory.run();                       // collective; the work happens on the master
```

After `run()` the cut mesh is still whole on the master. Do not call `Mesh::partition()` for
distribution: it does nothing off the master, and nothing in `src/` calls it. The FEM
`Kernel` constructor partitions the mesh on the master (`Kernel::partition_mesh()` →
`graph::metis_partition`) and then distributes it on every rank. See
[Mesh Distribution (MPI)](../../mesh/doc/mesh_usage_guide.md) in the mesh guide.

The cohomology computation itself is not parallelized.

---

### Memory Ownership

**Critical ownership rules:**

| Object | Created By | Owned By | Caller Action |
|--------|-----------|----------|---------------|
| `SideSet*` (cuts) | CutFactory | **Mesh** | **Do NOT delete** |
| `Node*` (abstract) | CutFactory | **Mesh** | **Do NOT delete** |
| `Node*` (duplicates) | CutFactory | **Mesh** | **Do NOT delete** |
| `Node*` (orphaned) | — (existing mesh nodes) | **Mesh** | **Do NOT delete** |
| `Cochain*` (generators) | Cohomology | **Cohomology** | **Do NOT delete** |
| `Chain*` (generators) | Homology | **Homology** | **Do NOT delete** |
| `Protoshell*` | **User** | **User** | **Caller must delete** |

**Key rule:** CutFactory transfers all created mesh entities (nodes, sidesets) to the mesh. The mesh takes ownership and will delete them in its destructor.

The cut sidesets are reached through `Mesh::sidesets()` (see "Result Retrieval"). The
abstract and orphaned nodes are in `Mesh::abstract_nodes()` and `Mesh::orphaned_nodes()`,
which are non-owning views onto nodes that `Mesh::nodes()` owns.

```cpp
// CORRECT
tFactory.run();
mesh::SideSet * tCut = tMesh->sideset( tCutID );   // reference only
delete tMesh;                                     // ~Mesh deletes cuts and nodes

// WRONG
delete tMesh->sideset( tCutID );   // ❌ DOUBLE DELETE when ~Mesh runs
```

---

### Thread Safety

The homology module is **NOT thread-safe**:

- `SimplicialComplex::reduce_*()` modifies internal state (chain/cochain maps)
- `CutFactory::run()` modifies mesh topology (node duplication, element relinking)
- `smithForm()` and friends modify their input matrix in place (`aMat` is taken by reference)

**OpenMP:** Do not call homology functions from parallel regions

**MPI:** The computation runs on the master rank only, but `run()` must be called on every rank (see MPI Serial Limitation above)

---

### Topological vs FEM Correctness

**Critical distinction:** The homology module provides two separate guarantees:

#### Topological Guarantees

✅ **What homology guarantees:**
- Generators span correct cohomology groups H^k
- Algorithms produce topologically equivalent results (different bases, same groups)
- Betti numbers are exact (integer arithmetic, no rounding errors)
- Orientation is consistent after cleanup

#### FEM Guarantees

⚠️ **What homology does NOT automatically guarantee:**
- Thin cuts form valid 2-manifolds → **Requires manifold filtering (may fail)**
- Node duplication yields independent DOFs → **Requires proper element relinking (usually works)**
- No element has mixed φ⁺/φ⁻ nodes → **Requires correct cut-side detection (usually works)**
- Static condensation produces correct constraints → **Requires correct FEM implementation (user responsibility)**

**Implication:** A topologically correct cohomology generator may still produce a non-manifold thin cut that breaks FEM assembly. Always inspect the cut sidesets (debug output, ParaView) — there is no `is_manifold()` query on `SideSet`.

---

### Mesh Requirements

**Preconditions for correctness:**

✅ **Required mesh properties:**
- Mesh must be **connected** (single component, or separate components with distinct domain types)
- Mesh must be **orientable** (no Möbius strips)
- Elements must form a **manifold** (no hanging faces, no non-manifold edges)
- Facet orientations must be **computable** (consistent normal directions)

⚠️ **Topology requirements:**
- Phi/non-phi domains must be **correctly tagged** (see `Topology::run()`)
- Thin-shell interfaces must be **properly labeled** as sidesets
- Terminals (if used) must **exist as mesh entities**: block ids (2D bulk conductor), sideset
  ids (3D bulk conductor) or curve ids (thin shell); see `set_terminals()`
- Block domain types must be set before `Topology::run()`: `Air`, `Buffer` and `Ferro` blocks
  are phi blocks, `Conductor` and `Coil` blocks are non-phi blocks, and a block of any other
  type (including `Default`) is in neither list

❌ **Unsupported configurations:**
- Non-orientable manifolds (Möbius strip, Klein bottle)
- Meshes with dangling tetrahedra (elements connected only at nodes/edges)

**Validation pattern** (master rank, where the mesh and the topology lists are populated):
```cpp
// Check basic mesh properties
BELFEM_ERROR(tMesh->number_of_elements() > 0, "Empty mesh");
BELFEM_ERROR(tMesh->number_of_blocks() > 0, "No blocks defined");

// Check topology analysis succeeded
mesh::Topology tTopology(tMesh);
tTopology.run();
BELFEM_ERROR(tTopology.phi_block_ids().length() > 0,
             "No phi-domain blocks identified");
```

---

### Scratch-Field Usage

**Warning:** the module uses the mesh entities' flag bits (`flag_*()` / `unflag_*()`) as
scratch and, in the SPFA infeasibility path, writes `Element::level`
(`cl_Cohomology.cpp`, `set_level`). It does not write `owner`. It does write `index`: the RCM
renumbering in `CutFactory::compute_rcm_problem()` sets every node index, and
`compute_cohomologies()` calls `update_node_indices()`, `update_edge_indices()` and
`update_face_indices()`. The same RCM step overwrites, or creates, the node field `"phi"`
with the new node indices.

**Consequences:**
- Do not rely on entity flags, `Element::level`, entity indices or a field named `"phi"` from
  before `CutFactory::run()`
- Run homology **before** the mesh is distributed, because cohomology needs the undistributed
  mesh (see "MPI Serial Limitation"), not because of `owner`. Distribution happens in the FEM
  `Kernel` constructor, so build the cuts before the kernel

---

## Table of Contents

1. [Quick-Start](#quick-start)
2. [Glossary](#glossary)
3. [Common Pitfalls](#common-pitfalls)
4. [Introduction](#introduction)
5. [Module Architecture](#module-architecture)
6. [CutFactory - Main Entry Point](#cutfactory---main-entry-point)
7. [Cohomology Computation](#cohomology-computation)
8. [Homology Computation](#homology-computation)
9. [Cut Algorithm Selection](#cut-algorithm-selection)
10. [Topology Analysis](#topology-analysis)
11. [Thin Shell Handling](#thin-shell-handling)
12. [Cut Processing](#cut-processing)
13. [Data Structures](#data-structures)
14. [Advanced Topics](#advanced-topics)
15. [Common Patterns](#common-patterns)
16. [Troubleshooting](#troubleshooting)
17. [Usage Examples](#usage-examples)

---

## Quick-Start

**Minimal example to generate cuts for a toroidal mesh:**

```cpp
#include "cl_CutFactory.hpp"
#include "cl_Topology.hpp"

// Every rank: the master reads the file, the others get an empty mesh
Mesh * tMesh = new Mesh( "toroid.exo" );

// Master: set the block domain types (here block 1 is air, block 2 the conductor)
if ( comm_rank() == 0 )
{
    tMesh->block( 1 )->set_domain_type( DomainType::Air );
    tMesh->block( 2 )->set_domain_type( DomainType::Conductor );
}

// Analyze topology
mesh::Topology tTopology( tMesh );
tTopology.run();
tTopology.synchronize_maps();

// Generate cuts
Cell< Protoshell * > tProtoshells ;   // no thin shells; must be an lvalue
mesh::CutFactory tFactory( tMesh, &tTopology, tProtoshells,
                           mesh::CutAlgorithm::PellikkaGeneralized,
                           false );   // aUseEnrichment
tFactory.set_periodicity( nullptr );
tFactory.run();

// Results (master): one abstract node in tMesh->abstract_nodes(), and one sideset
// labelled cut_<n>, per H^1 generator
```

The cut sidesets are stored in the mesh and labelled `cut_<n>` in `Mesh::sidesets()`. The
abstract nodes in `Mesh::abstract_nodes()` carry the jump DOFs (see "Result Retrieval").

**Leave the cut sidesets at `DomainType::Default`.** These are thick cuts: the jump is carried by
the abstract nodes, not by dofs on the cut surface. The Maxwell field list builds a thin-cut dof
table for every sideset typed `DomainType::Cut` (`cl_Maxwell_FieldList.cpp`), so retyping a
`cut_<n>` sideset would add dofs the thick-cut formulation does not use. `MaxwellFactory` does not
retype them, so the Maxwell solver does not select them.

---

## Glossary

**Core Concepts:**

| Term | Definition |
|------|------------|
| **Chain** | Formal sum of k-simplices with integer coefficients (homology side): `c = Σ aᵢ σᵢ` where aᵢ ∈ ℤ |
| **Cochain** | Dual to chain; element of cochain complex (cohomology side): `φ : C_k → ℤ` |
| **Thick Cut** | Cohomology generator as a set of directed edges with ±1 coefficients (topological object) |
| **Thin Cut** | FEM-compatible surface where φ is discontinuous (minimal 2-manifold derived from thick cut). Worked example: [thick_thin_cuts_and_conjugate_edges.md](thick_thin_cuts_and_conjugate_edges.md) |
| **Conjugate Edges** | Face-loop edges bounding the thin cut, dual to the apex edges the thick cut intersects (see [thick_thin_cuts_and_conjugate_edges.md](thick_thin_cuts_and_conjugate_edges.md)) |
| **Abstract Node** | Node not belonging to any element; carries DOF for jump condition [φ] = φ⁺ - φ⁻ |
| **Duplicate Node** | Copy of original node created during cut processing; allows φ⁺ ≠ φ⁻ on cut surface |
| **Simplicial Complex** | Collection of simplices (nodes, edges, faces, elements) with incidence relations |
| **Boundary Operator** | Map ∂_k : C_k → C_{k-1} (takes k-chain to its (k-1)-boundary) |
| **Coboundary Operator** | Map δ^k : C^k → C^{k+1} (dual to boundary operator) |
| **Smith Normal Form** | Diagonal matrix decomposition Q⁻¹AR = D revealing Betti numbers and torsion (`smithForm` returns Q and Q⁻¹ as `Q`, `Q_`) |
| **Betti Number** | Rank of k-th homology/cohomology group; counts k-dimensional "holes" |
| **Homology Group** | H_k = Ker(∂_k) / Im(∂_{k+1}) (k-cycles modulo k-boundaries) |
| **Cohomology Group** | H^k = Ker(δ^k) / Im(δ^{k-1}) (k-cocycles modulo k-coboundaries) |
| **Terminal** | Electrical connection for current injection/extraction, given by mesh entity ids: block ids (2D bulk), sideset ids (3D bulk) or curve ids (thin shell) |
| **Protoshell** | Configuration object for thin shell structures (sidesets, terminals, material, thickness) |
| **Phi-domain** | Region where magnetic scalar potential φ is solved (non-conducting) |
| **Non-phi-domain** | Conducting region or external boundary (not part of φ solution domain) |
| **Manifold** | Surface where every edge belongs to exactly 2 faces (topologically valid) |
| **Non-manifold** | Surface with edges shared by >2 faces or other defects (invalid for FEM) |
| **Double Pocket** | Topological defect: small cluster of faces attached via single articulation point |

**Algorithm Terms:**

| Term | Definition |
|------|------------|
| **Elementary Collapse** | Removal of k-simplex and its free (k-1)-face (preserves homology) |
| **Interior Face Reduction** | Merging two k-simplices sharing a common face |
| **Generalized Combine** | Enhanced reduction updating boundaries after removing shared faces |
| **Pellikka Algorithm** | Original reduction via elementary collapses and combines |
| **PellikkaGeneralized** | Enhanced Giard algorithm combining CCR with Pellikka |
| **CCR** | Chain Complex Reduction (Kaczynski); classical acyclic matching |
| **BeltedTree** | Spanning tree approach with belt fasteners for cohomology |
| **Manifold Filtering** | Post-processing to remove non-manifold artifacts (double pockets, etc.). Design of record — **not implemented**; see the design note under "Manifold Cleanup" |
| **Tarjan's Algorithm** | DFS-based articulation point detection, part of the manifold-cleanup design. No implementation exists in the tree |

**FEM Integration:**

| Term | Definition |
|------|------------|
| **Jump Condition** | Discontinuity enforced on cut: [φ] = φ⁺ - φ⁻ = I (Ampere current) |
| **Static Condensation** | FEM technique to eliminate internal DOFs via transformation matrix T |
| **Poincaré-Lefschetz Duality** | Theorem allowing thick-to-thin cut conversion (edge set → surface) |
| **h-φ Formulation** | Magnetostatic FEM using magnetic field H = -∇φ in non-conducting regions |
| **Ampere's Law** | ∮ H·dl = I for any loop encircling conductor with current I |

---

## Common Pitfalls

**⚠️ Read this section to avoid frequent mistakes.**

### 1. Building the CutFactory Before the Topology Has Run

**Problem:** Topology not analyzed before creating CutFactory.

```cpp
Cell< Protoshell * > tProtoshells ;

// ❌ WRONG
mesh::Topology tTopology( aMesh );
mesh::CutFactory tFactory( aMesh, &tTopology, tProtoshells,
                           mesh::CutAlgorithm::PellikkaGeneralized, false );
// ❌ Topology::run() and synchronize_maps() were never called

// ✅ CORRECT
mesh::Topology tTopology( aMesh );
tTopology.run();                // Analyze FIRST
tTopology.synchronize_maps();
mesh::CutFactory tFactory( aMesh, &tTopology, tProtoshells,
                           mesh::CutAlgorithm::PellikkaGeneralized, false );
tFactory.set_periodicity( nullptr );
tFactory.run();
```

**Consequence:** the `CutFactory` constructor reads `groups( DomainType::ThinShell )`. Without
the domain-type map, the lookup stops with `BELFEM_ERROR` ("Key ... not found in map"). On a
rank other than 0, the map exists only after `synchronize_maps()`.

---

### 2. Running Homology After Mesh Partitioning

**Problem:** Distributing the mesh before homology breaks the simplicial complex.

The FEM `Kernel` constructor distributes the mesh: it calls `Kernel::partition_mesh()` on the
master and then distributes the mesh on every rank. Therefore, call `CutFactory::run()` before
building the kernel, as `MaxwellFactory::create_magnetic_kernel()` does. `Mesh::partition()` does
not distribute anything; it only assigns ownerships on the master.

**Consequence:** Incomplete generators, missing cuts, or crashes due to missing connectivity.

---

### 3. Deleting Cuts or Abstract Nodes

**Problem:** Manually deleting mesh entities owned by the mesh.

```cpp
// ❌ WRONG
delete tMesh->sideset( tCutID );        // ❌ DOUBLE DELETE - mesh owns this!
delete tMesh->abstract_nodes()( 0 );    // ❌ same for abstract nodes

// ✅ CORRECT
mesh::SideSet * tCut = tMesh->sideset( tCutID );
// Just use it - ~Mesh deletes it
```

**Consequence:** Double delete → crash or memory corruption.

---

### 4. Forgetting Edge/Face Creation

**Problem:** SimplicialComplex needs edges and faces, but mesh doesn't have them.

```cpp
// ❌ WRONG
Mesh* tMesh = new Mesh("mesh.msh");
mesh::SimplicialComplex tComplex(tMesh);  // ❌ No edges!

// ✅ CORRECT
Mesh* tMesh = new Mesh("mesh.msh");
tMesh->create_edges();  // Create edges first
tMesh->create_faces();  // Create faces (for 3D)
mesh::SimplicialComplex tComplex(tMesh);
```

**Consequence:** `number_of_ksimplices(1) == 0` → cohomology computation fails.

---

### 5. Mixing Chain and Cochain

**Problem:** Using Chain where Cochain expected (or vice versa).

```cpp
// CutProcessor takes the whole Cohomology object in its constructor
// (cl_CutProcessor.hpp:94-101); it is never fed individual Chain or Cochain
// pointers. Homology generators are Chain*, cohomology generators are
// Cochain* — do not pass one where the other is expected, e.g.
// Cohomology::updatekGeneratorsFromHomology() takes Cell<Chain*>&.
Chain*   tChain   = tHomology.get_Generators()( 1 )(0);    // homology
Cochain* tCochain = tCohomology.get_Generators()( 1 )(0);  // cohomology
```

**Consequence:** Compile error or incorrect cut generation.

---

### 6. Not Checking Manifold Status

**Problem:** Non-manifold cuts break FEM assembly.

```cpp
// ⚠️ RISKY
tFactory.run();
// Assume all cuts are manifold

// ✅ SAFER (pseudo-code: SideSet has no edges() accessor, and run() removes the mesh
// edges, so walk the facets' edges yourself; tCutSideSets holds the cut_<n> sidesets
// collected from tMesh->sidesets())
tFactory.run();
for (SideSet* cut : tCutSideSets) {
    // Check manifold status
    uint nonManifoldEdges = 0;
    for (Edge* edge : cut->edges()) {
        if (edge->number_of_facets() != 2) {
            ++nonManifoldEdges;
        }
    }
    BELFEM_ERROR(nonManifoldEdges == 0,
        "Cut %lu has %u non-manifold edges",
        (luint)cut->id(), nonManifoldEdges);
}
```

**Consequence:** FEM assembly fails due to edges shared by >2 faces.

---

### 7. Using Wrong Algorithm for Mesh Size

**Problem:** Using CCR (the theory/debugging path) for production meshes.

```cpp
// ❌ Theory/debugging path
mesh::CutFactory tFactory( aMesh, &tTopology, tProtoshells, mesh::CutAlgorithm::CCR, false );
tFactory.run();

// ✅ Production default
mesh::CutFactory tFactory( aMesh, &tTopology, tProtoshells, mesh::CutAlgorithm::PellikkaGeneralized, false );
tFactory.run();
```

**Consequence:** No timing comparison between the algorithms is recorded (see README, "Cut Algorithms").

---

### 8. Running CutFactory Twice

**Problem:** Calling `run()` multiple times on same factory.

```cpp
// ❌ WRONG
tFactory.run();  // First run - mutates mesh
tFactory.run();  // ❌ UNDEFINED BEHAVIOR - mesh already modified

// ✅ CORRECT
tFactory.run();  // Run exactly once
// Create new factory if you need different cuts
```

**Consequence:** Crash, incorrect cuts, or memory corruption.

---

### 9. Ignoring Orphaned Nodes

**Problem:** Orphaned nodes are nodes of `Mesh::nodes()` that no block element uses and that
are not abstract nodes (`CutFactory::collect_orphan_nodes()`). They stay in the mesh, which
owns them.

```cpp
// ⚠️ CHECK (master)
tFactory.run();
if (tFactory.orphaned_nodes().size() > 0) {
    message(InfoLevel::Default, "Warning: %lu orphaned nodes",
            (long unsigned int) tFactory.orphaned_nodes().size());
    // Do NOT delete them: ~Mesh deletes every node in Mesh::nodes()
}
```

**Consequence:** Deleting an orphaned node is a double free when `~Mesh` runs. Downstream,
`MaxwellFactory` hands the list to the equation (`set_orphaned_nodes`).

---

### 10. Incorrect Terminal-to-Shell Mapping

**Problem:** The second argument of `set_terminals()` is mistaken for a list of shell numbers.

`aThinShellIndices` lists, in ascending order, the **positions in `aTerminals`** whose entries
belong to a thin shell. It does not map terminals to protoshells.
`Homology::suggest_Homology()` walks `aTerminals` and compares each position `i` with the next
unused entry of `aThinShellIndices`, so the list must be ascending.

```cpp
// ❌ WRONG: one entry per terminal, read as "shell number"
Cell<Cell<id_t>> tTerminals = {{101, 102}, {7, 8}};
Cell<id_t> tShellIndices = {0, 0};   // marks position 0 as a thin shell, twice

// ✅ CORRECT: entry 0 is a bulk conductor (sideset ids in 3D),
// entry 1 is a thin shell (curve ids); only position 1 is listed
Cell<Cell<id_t>> tTerminals = {{101, 102}, {7, 8}};
Cell<id_t> tShellIndices = {1};
```

**Consequence:** A bulk entry is read as curve ids (or the reverse), so the suggested homology
is built from the wrong entities.

---

## Introduction

The `homology` module provides computational topology tools for analyzing and processing finite element meshes, with a focus on electromagnetic field simulations in multiply-connected domains.

**Primary use case:** Generating topological cuts for enforcing current constraints in magnetostatic h-φ formulations.

### What is a Cohomology Cut?

In multiply-connected domains (e.g., toroidal magnets), the magnetic scalar potential φ must satisfy Ampere's law:

```
∮_C H · dl = I
```

Where C is any closed loop encircling a conductor carrying current I. To enforce this in FEM:

1. **Compute cohomology generator**: Identifies topological "hole" in the mesh
2. **Convert to thin cut**: Minimal surface where φ is discontinuous
3. **Duplicate nodes**: Create separate DOFs for φ⁺ and φ⁻ on cut surface
4. **Enforce jump condition**: [φ] = φ⁺ - φ⁻ = I via static condensation

**Location:** `src/homology/`
**CMake Target:** `belfem_homology`

---

## Module Architecture

### File Organization

```
src/homology/
├── Entry Point
│   └── cl_CutFactory.{hpp,cpp}
│
├── Computational Core
│   ├── cl_Cohomology.{hpp,cpp}
│   ├── cl_Homology.{hpp,cpp}
│   ├── cl_SimplicialComplex.{hpp,cpp}
│   └── fn_Smith.{hpp,cpp}
│
├── Data Structures
│   ├── cl_Chain.{hpp,cpp}
│   ├── cl_Cochain.{hpp,cpp}
│   ├── cl_CutData.{hpp,cpp}
│   └── cl_CutSet.{hpp,cpp}
│
├── Processing
│   ├── cl_CutProcessor.{hpp,cpp}
│   ├── cl_CutProcessorManual.{hpp,cpp}
│   └── cl_InterfaceProcessor.{hpp,cpp}
│
├── Configuration
│   ├── cl_Topology.{hpp,cpp}
│   └── en_CutAlgorithm.hpp
│
└── Algorithms
    └── cl_BeltedTree.{hpp,cpp}
```

### Dependencies

- **mesh**: Mesh, Node, Element, Facet, SideSet classes
- **linalg**: Matrix, Vector for boundary/coboundary matrices
- **containers**: Cell, OrderedMap, Set for data structures
- **sparse**: Integer matrix solvers for Smith Normal Form
- **comm**: MPI support (computations currently serial)

---

## CutFactory - Main Entry Point

### Purpose

`CutFactory` is the **primary user-facing class** that orchestrates the entire cut generation pipeline.

**File:** `cl_CutFactory.{hpp,cpp}`

### Constructor

```cpp
CutFactory::CutFactory(
    Mesh*               aMesh,
    Topology*           aTopology,
    Cell<Protoshell*>&  aProtoshells,
    const CutAlgorithm  aAlgorithm,
    const bool          aUseEnrichment
);
```

**Parameters:**
- `aMesh`: Finite element mesh to process
- `aTopology`: Pre-analyzed mesh topology (identifies domains, interfaces)
- `aProtoshells`: Configuration for thin shell structures
- `aAlgorithm`: Cut generation algorithm (see [Cut Algorithm Selection](#cut-algorithm-selection))
- `aUseEnrichment`: controls one step of `run()`. When it is true,
  `restore_thin_shell_sidesets()` moves the temporary facets of each thin-shell sideset into a
  new hidden sideset of domain type `EnrichedInterface`, instead of deleting them

Neither of the last two arguments has a default; both are passed explicitly. `aProtoshells`
is a non-const reference, so it must be an lvalue (`{}` does not compile); pass an empty
`Cell< Protoshell * >` when there are no thin shells.

**When the mesh is periodic, call `set_periodicity()` before `run()`.** The member
`mPeriodicFactory` defaults to `nullptr`, which `run()` reads on the master as "no periodicity".
When it computes the cuts, `MaxwellFactory` calls `set_periodicity()` on the master, passing its
own factory pointer.

**⚠️ Ownership:** CutFactory does **not** take ownership of aMesh or aTopology. Caller must ensure they outlive the factory.

### Key Methods

#### Workflow Methods

```cpp
void set_terminals(
    const Cell<Cell<id_t>>& aTerminals,
    const Cell<id_t>&       aThinShellIndices
);
```

Sets up electrical terminal definitions for current injection/extraction points. The lists
are copied; they are consumed by `Homology::suggest_Homology()` during `run()`.

- `aTerminals`: one entry per terminal condition. `MaxwellFactory::create_terminal_list()`
  builds one entry per current, circuit-current and voltage boundary condition, from its
  `domains()`. Each entry holds mesh entity **ids**, not node ids: block ids for a 2D bulk
  conductor, sideset ids for a 3D bulk conductor, curve ids for a thin shell. The first half of
  an entry is the input terminal; in 3D, the second half is the output terminal. In 2D, the
  second half of a bulk entry is not read.
- `aThinShellIndices`: the **positions** in `aTerminals` whose entries belong to a thin shell,
  in ascending order. It is not a list of protoshell numbers.

```cpp
void set_periodicity( PeriodicityFactory * aPeriodicFactory );
```

Stores the periodicity factory, or `nullptr` when periodicity is not used (the default). Call
it before `run()` (see above).
On the master, `run()` calls `create_periodicity()` on a non-null factory and installs the
result on the mesh.

```cpp
bool create_thin_shell_cuts();
```

Thin-shell preprocessing. `run()` does **not** call it. `MaxwellFactory` calls it on the master
after `set_periodicity()` and before `set_terminals()`/`run()`. It returns `false` and does
nothing when the topology has no `ThinShell` sidesets. Otherwise it creates the thin-shell
curves, duplicates and relinks the nodes and facets on the shell sidesets, closes the terminal
loops, and returns `true`. Afterwards `thin_shell_master_nodes()` and `thin_shell_slave_nodes()`
hold the shell node pairs.

```cpp
void run();
```

**Main execution pipeline.** On the master rank it:
1. Rebuilds the mesh connectivity, creates edges (and faces in 3D) and installs the
   periodicity, if any
2. Renumbers the nodes with reverse Cuthill-McKee and writes the `"phi"` node field
3. Computes the cohomology: the suggested homology from the terminals, the complex reduction
   for the selected algorithm, `Cohomology`, then recombination with the suggested homology and
   a second `clean()`
4. Converts thick cuts to thin cuts (`CutProcessor`): node duplication, abstract nodes, cut
   sidesets; then duplicates interface nodes (`InterfaceProcessor`) and collects orphaned nodes
5. Restores the thin-shell sidesets and computes the element adjacencies
6. Removes the edges and faces and finalizes the mesh

All of this work runs on the master. `compute_rcm_problem()` and the end of `run()` each call
`comm_barrier()`, so every rank must call `run()`. Generator cleaning, the step that runs
instead of a manifold filter, is part of step 3 (see "Manifold Cleanup").

**⚠️ Call exactly once.** Calling `run()` multiple times is undefined behavior.

#### Result Retrieval

```cpp
Cell<Node*>& abstract_nodes();
```

Returns the abstract nodes (new degrees of freedom), one per H^1 generator, which carry the
cut jumps. This is a reference to the mesh's own list, `Mesh::abstract_nodes()`; the nodes
themselves are owned by `Mesh::nodes()`.

**Usage:** `MaxwellFactory` passes `Mesh::abstract_nodes()` to the equation
(`IWG::set_abstract_nodes`) before the kernel is built; later
`DofManager::extract_abstract_dofs_from_mesh()` reads the abstract DOF data. There is no
per-node "add abstract DOF" call.

**Cut sidesets.** `CutFactory` has no accessor for the cuts. `run()` adds them to the mesh
through `CutData::add_thin_cut_sidesets_to_mesh()`, one for each H^1 generator, with the label `cut_`
followed by a zero-padded number and domain type `Default`. Find them through
`Mesh::sidesets()` by label, and leave the domain type at `Default` (see Quick-Start).

```cpp
Cell<Node*>& orphaned_nodes();
```

Returns the nodes of `Mesh::nodes()` that no block element uses and that are not abstract
nodes, collected after node duplication. This is a reference to `Mesh::orphaned_nodes()`, a
non-owning view.

**Usage:** Read only. Do **not** delete these nodes: `~Mesh` deletes them.

```cpp
Vector<id_t>& thin_shell_boundaries();
```

Returns IDs of boundaries created for thin shell sidesets.

```cpp
Cell<Facet*>& thin_shell_facets(const id_t aID);
```

Returns the facets of the temporary sideset at **position** `aID` in
`mTemporaryThinShellSidesets`. This is a positional index, not a sideset id. The list is filled
by `create_thin_shell_cuts()`. `run()` restores the thin-shell sidesets, deletes the temporary
sidesets and clears the list. This accessor is valid only between `create_thin_shell_cuts()`
and `run()`; do not call it after `run()`.

#### Debugging

```cpp
void save_curve_debug_meshes();
```

The only public debug writer on `CutFactory`: it writes one `curve_<id>.vtk` per mesh curve.
There is no `CutFactory::save_debug_meshes()`. The per-cut meshes `cut_<index>.vtk` come from
`CutProcessor::save_debug_meshes()`. `CutFactory` calls it inside
`compute_thin_cuts_and_duplicate_interface_nodes()` only in debug builds
(`#if !defined( NDEBUG ) || defined( DEBUG )`). No code writes separate thick-cut, thin-cut or
manifold files.

### Usage Pattern

```cpp
#include "cl_CutFactory.hpp"
#include "cl_Topology.hpp"

// Every rank (see "MPI Serial Limitation"); aMesh holds the whole mesh on the master

// 1. Analyze mesh topology
mesh::Topology tTopology( aMesh );
tTopology.run();
tTopology.synchronize_maps();

// 2. Create CutFactory
mesh::CutFactory tFactory( aMesh, &tTopology, mProtoshells,
                           mesh::CutAlgorithm::PellikkaGeneralized,
                           false );   // no enrichment
tFactory.set_periodicity( nullptr );

// 3. Thin shells and terminals (master; MaxwellFactory order)
if ( comm_rank() == 0 )
{
    tFactory.create_thin_shell_cuts();

    // 3D: entry 0 = input sideset 101, output sideset 102 (a bulk conductor)
    Cell< Cell< id_t > > tTerminals = { { 101, 102 } };
    Cell< id_t > tThinShellIndices ;   // no thin-shell entries
    tFactory.set_terminals( tTerminals, tThinShellIndices );
}

// 4. Run pipeline
tFactory.run();

// 5. Results (master): cut sidesets cut_<n> in aMesh->sidesets(),
//    jump DOFs in aMesh->abstract_nodes()
```

---

## Cohomology Computation

### Purpose

`Cohomology` computes cohomology groups H^k = Ker(δ^k) / Im(δ^{k-1}) from a simplicial complex.

**File:** `cl_Cohomology.{hpp,cpp}`

### Constructor

```cpp
// From simplicial complex
Cohomology::Cohomology(SimplicialComplex* aComplex, Mesh* aFullMesh);

// With belted tree
Cohomology::Cohomology(SimplicialComplex* aComplex, Mesh* aFullMesh, BeltedTree* aTree);
```

**Parameters:**
- `aComplex`: Pre-constructed (and reduced) simplicial complex
- `aFullMesh`: Original full mesh (for node/element lookups)
- `aTree`: Belted tree for alternative algorithm (optional)

Both constructors leave the object with computed, cleaned generators, ready
to use. The `SimplicialComplex` constructor computes all H^k groups; the
BeltedTree constructor takes its H^1 generators from the tree.

### Key Methods

```cpp
void cohomologyGroupOfChainComplex();
```

**Main computation method.** Performs:

1. Kernel/image decomposition of the coboundary maps
2. Quotient group computation (Smith Normal Form)

Complex reduction happens beforehand on the `SimplicialComplex`; generator
extraction is `generatorsOfCohomology()` and cleaning is `clean()`, both of
which the constructor runs afterward.

**⚠️ Called internally** by the `SimplicialComplex`-based constructor. The
BeltedTree constructor copies its generators from the tree instead. Do not
call it yourself.

```cpp
void quotientGroup();
```

Computes quotient space H^k = Ker(δ^k) / Im(δ^{k-1}) via Smith Normal Form.

**Note:** Typically called internally by `cohomologyGroupOfChainComplex()`.

```cpp
void generatorsOfCohomology();
```

Extracts cohomology generators from Smith decomposition.

**Output:** Populates `mGenerators` with Cochain objects representing basis elements.

```cpp
void clean();
```

Cleans generator coefficients so every coefficient is in {-1, 0, 1}.

**Method:** Dispatches to `clean_spfa()`: SPFA feasibility certificate, greedy
coboundary rectification, then `remove_cut_pockets()`.

**Note:** Every constructor already calls `clean()`, so a constructed
`Cohomology` always has unit-coefficient generators. `CutFactory` calls it a
**second** time after `updatekGeneratorsFromHomology()`, because integer
recombination of unit generators is generally non-unit again. Both runs are
required: the second run is the only clean of the final generators, and the
downstream `CutData` machinery can only represent coefficients in {-1, 0, 1}.
See `cohomology_algorithms.md`, "Cleaning runs twice on the factory path".

```cpp
Cell<Cell<Cochain*>>& get_Generators();
```

Returns cohomology generators organized by dimension k.

**Structure:** `mGenerators[k]` contains generators for H^k.

```cpp
Cell<Matrix<int>>& get_CoboundaryMatrix();
```

Returns coboundary operator matrices δ^k.

**Structure:** `mCoboundaryMatrix[k]` is the δ^k : C^k → C^{k+1} matrix.

### Visualization

```cpp
void create_kGeneratorsField(uint k, Mesh* aMesh, string aLabel);
```

Creates mesh fields to visualize cohomology generators.

**Parameters:**
- `k`: Dimension (0 for nodes, 1 for edges, 2 for faces)
- `aMesh`: Mesh to attach fields to
- `aLabel`: Field name prefix

**Output:** Creates fields showing generator support with ±1 coefficients.

### Usage Pattern

```cpp
#include "cl_Cohomology.hpp"

// Build the simplicial complex on the flagged mesh entities and reduce it
mesh::SimplicialComplex tComplex(aFullMesh, true);
tComplex.coreduce_complexPellikkaGeneralized();

// The constructor computes the cohomology groups and cleans the generators
mesh::Cohomology tCohomology(&tComplex, aFullMesh);

// Extract H^1 generators (most common for electromagnetic cuts)
Cell<Cochain*>& generators = tCohomology.get_Generators()( 1 );

message(InfoLevel::Default, "Found %d cohomology generators", generators.size());

// Visualize generators
tCohomology.create_kGeneratorsField(1, aFullMesh, "H1_gen");

// Access individual generator
Cochain* gen = generators(0);
for (auto& [edge_id, coeff] : gen->get_simplices_map()) {
    message(InfoLevel::Verbose, "Edge %lu: coefficient %d", edge_id, coeff);
}
```

---

## Homology Computation

### Purpose

`Homology` computes homology groups H_k = Ker(∂_k) / Im(∂_{k+1}), the dual counterpart to cohomology.

**File:** `cl_Homology.{hpp,cpp}`

**Use cases:**
- Relative homology for terminal definitions
- Orientation computation for current flow direction
- Topological analysis of conductor geometries

### Constructor

```cpp
// With terminal suggestion
Homology::Homology(
    Mesh*               aMesh,
    Cell<Cell<id_t>>    aTerminals,
    Cell<id_t>          aThinShells
);

// From simplicial complex
Homology::Homology(SimplicialComplex* aComplex, Mesh* aFullMesh);
```

The two constructors do different things:

- The **terminal constructor** only calls `suggest_Homology( aTerminals, aThinShells )`. It
  fills `get_Generators()( 1 )` with suggested 1-chains: the boundary of each input terminal,
  and in 3D also of each output terminal. It builds no simplicial complex and no boundary
  matrices, so `homologyGroupOfChainComplex()` and `get_BoundaryMatrix()` have nothing to work
  on. `aTerminals` and `aThinShells` follow the `CutFactory::set_terminals()` convention:
  entity ids per entry, and ascending positions of the thin-shell entries. `CutFactory` uses
  this constructor, then `reorient_generators()`.
- The **complex constructor** builds the boundary matrices from `aComplex`, computes the
  groups (`homologyGroupOfChainComplex()`) and extracts the generators.

### Key Methods

```cpp
void homologyGroupOfChainComplex();
```

Main computation method: kernel/image decomposition of the boundary matrices, then the
quotient groups. It needs the boundary matrices that only the complex constructor builds, and
that constructor already calls it. Generator extraction is `generatorsOfHomology()`.

```cpp
void suggest_Homology(
    Cell<Cell<id_t>> aTerminals,
    Cell<id_t>       aThinShellIndices
);
```

Suggests homology generators based on terminal connectivity.

**Use case:** Pre-seeds homology computation with physically meaningful cycles.

```cpp
Cell<int> generators_orientation(Cell<Vector<real>>& aVectors);
```

Computes orientation of generators relative to specified direction vectors.

**Returns:** Cell of ±1 values indicating generator orientation.

**Use case:** Ensures current flow direction matches physical conventions.

```cpp
void reorient_generators();
```

Flips generator orientation based on computed directions.

```cpp
Cell<Cell<Chain*>>& get_Generators();
```

Returns homology generators organized by dimension k.

```cpp
Cell<Matrix<int>>& get_BoundaryMatrix();
```

Returns boundary operator matrices ∂_k.

### Usage Pattern

```cpp
#include "cl_Homology.hpp"

// Terminal hints (3D): one bulk conductor, input sideset 101, output sideset 102
Cell< Cell< id_t > > tTerminals = { { 101, 102 } };
Cell< id_t > tThinShellIndices ;   // no thin-shell entries

// The constructor builds the suggested generators; it computes no groups
mesh::Homology tHomology( aMesh, tTerminals, tThinShellIndices );

// Suggested 1-chains, one per terminal
Cell< mesh::Chain * > & tGenerators = tHomology.get_Generators()( 1 );

// Fix the sense of each loop, as CutFactory does
tHomology.reorient_generators();
```

---

## Cut Algorithm Selection

### Available Algorithms

Four algorithms are available via the `CutAlgorithm` enum:

**File:** `en_CutAlgorithm.hpp`

```cpp
enum class CutAlgorithm
{
    Pellikka            = 0,  // Original Pellikka (2013)
    CCR                 = 1,  // Chain Complex Reduction (Kaczynski 2004)
    BeltedTree          = 2,  // Spanning tree approach
    PellikkaGeneralized = 3,  // Enhanced Giard - RECOMMENDED
    UNDEFINED           = 4
};
```

### Algorithm Comparison

| Algorithm | Complexity | Performance | Manifold Issues | Recommendation |
|-----------|-----------|-------------|-----------------|----------------|
| **PellikkaGeneralized** | unmeasured | unmeasured | Rare | ✅ **Use this** (production default) |
| Pellikka | unmeasured | unmeasured | Occasional | Legacy support |
| CCR | unmeasured | unmeasured | Frequent | Theory/debugging |
| BeltedTree | unmeasured | unmeasured | Rare | Alternative |

### Implementation Details

**PellikkaGeneralized** (Giard et al., in preparation):
- **Location:** `cl_SimplicialComplex.cpp` — `pGeneralizedCombine()`, `reduce_complexPellikkaGeneralized()`, `pGeneralizedCocombine()`, `coreduce_complexPellikkaGeneralized()`
- **Method:**
  - Downward pass: `pReduce(k)` for k from d to 1 (elementary collapses)
  - Second pass (again from d down to 1): interleaves `pGeneralizedCombine(k)` with `pReduce(k-1)`; the cochain variant runs both passes upward, k = 0 … d−1
  - For cohomology: Reversed dimensions (k from 0 to d-1)
- **Intended advantages** (no complexity or timing measurement is recorded):
  - Smaller Smith Normal Form computation
  - Compact generator representations

**Pellikka** (Pellikka et al., 2013):
- **Location:** `cl_SimplicialComplex.cpp` — `pReduce()`, `pCombine()`, `reduceOmit()`
- **Method:**
  - `pReduce()`: Elementary collapses (remove k-simplex and free (k-1)-face)
  - `pCombine()`: Interior face reductions (merge neighboring simplices)
  - `reduceOmit()`: Initialization pass
- **Advantages:**
  - Well-tested in electromagnetic community
  - Widely used in GetDP, Sparselizard
  - Stable for moderate mesh sizes

**CCR** (Kaczynski et al., 2004):
- **Location:** `cl_SimplicialComplex.cpp` — `reduce_complexCCR()`, `coreduce_complexCCR()`
- **Method:** Classical chain complex reduction via acyclic matching
- **Advantages:**
  - Mathematically pure approach
  - Good for theoretical understanding
  - Reliable for small meshes
- **Disadvantages:**
  - May produce non-manifold cuts requiring extensive cleanup

**BeltedTree**:
- **Location:** `cl_BeltedTree.{hpp,cpp}`
- **Method:** Spanning tree construction with belt identification
- **Use case:** Alternative approach when reduction methods struggle

### Selection Guidance

```cpp
// Production code - use generalized Pellikka
mesh::CutFactory tFactory( aMesh, aTopology, aProtoshells,
                           mesh::CutAlgorithm::PellikkaGeneralized, false );

// Legacy compatibility:  mesh::CutAlgorithm::Pellikka
// Debugging/theory:      mesh::CutAlgorithm::CCR
// Requires the suggested homology (on by default): mesh::CutAlgorithm::BeltedTree
```

In a deck, `MaxwellFactory` reads the algorithm from the `algorithm` key of the `homology`
section (`to_cut_algorithm`).

**Benchmark results:** none recorded — see README, "Cut Algorithms".

---

## Topology Analysis

### Purpose

`Topology` analyzes mesh structure to identify domain types, interfaces, and boundary conditions.

**File:** `cl_Topology.{hpp,cpp}`

### Constructor

```cpp
Topology::Topology(Mesh* aMesh);
```

### Key Methods

```cpp
void run();                      // and run_on_enriched_mesh() for a reloaded .bfm mesh
```

Performs mesh topology analysis **on rank 0 only** (it does nothing on other ranks):
1. Collects the block and sideset domain types into a map of `DomainType` to id lists
2. Detects the sideset types from the neighboring blocks (`sideset_type()`)
3. Selects the phi and non-phi blocks
4. Selects the phi interface, phi boundary and phi periodic sidesets

The block domain types must already be set: `run()` reads them, it does not infer them.

```cpp
void synchronize_maps();
```

Collective. Broadcasts the domain-type map from rank 0 to every rank, so that `groups()`
works everywhere. It must be called before a `CutFactory` is constructed on a rank other than 0.

```cpp
const Map< id_t, DomainType > & block_map() const;
const Vector< id_t > & groups( const DomainType aType ) const;
```

`block_map()` maps a block id to its domain type. `groups( aType )` returns the ids recorded
for one domain type; an absent type stops with `BELFEM_ERROR`.

**Domain types** (`en_DomainType.hpp`, namespace `belfem`). The block types that matter here:
- `DomainType::Air`, `DomainType::Buffer`, `DomainType::Ferro`: phi blocks (φ is solved)
- `DomainType::Conductor`, `DomainType::Coil`: non-phi blocks
- `DomainType::ThinShell`: thin-shell sidesets
- `DomainType::Cut`: thin-cut sidesets; the field list builds a thin-cut dof table on them. The
  thick `cut_<n>` sidesets from `CutFactory` stay `Default` (see Quick-Start)

The sideset interface types include `InterfaceCondAir`, `InterfaceCondFerro`,
`InterfaceFerroAir`, `InterfaceAirCoil` and `InterfaceFerroCoil`.

```cpp
const Vector< id_t > & phi_block_ids() const;
const Vector< id_t > & non_phi_block_ids() const;
```

Returns lists of block IDs categorized by domain type.

```cpp
const Vector< id_t > & phi_interface_ids() const;
const Vector< id_t > & phi_boundary_ids() const;
const Vector< id_t > & phi_periodic_ids() const;
```

Returns sideset IDs. The interface list holds the `InterfaceCondAir`, `InterfaceCondFerro`,
`InterfaceAirCoil`, `InterfaceFerroCoil` and `ThinShell` sidesets. The boundary list holds the
air, buffer and ferro symmetry and antisymmetry sidesets and `BackgroundField`. The periodic
list holds the periodic sidesets.

### Usage Pattern

```cpp
#include "cl_Topology.hpp"

mesh::Topology tTopology( aMesh );
tTopology.run();

// on rank 0
for ( id_t tID : tTopology.phi_block_ids() )
{
    message( InfoLevel::Default, "Block %lu: phi domain\n", ( long unsigned int ) tID );
}

for ( id_t tID : tTopology.phi_interface_ids() )
{
    mesh::SideSet * tInterface = aMesh->sideset( tID );
    // Apply interface conditions
}

if ( tTopology.block_map()( 5 ) == DomainType::Conductor )
{
    // Apply source current density
}
```

---

## Thin Shell Handling

### Purpose

`Protoshell` configures thin shell structures for layered conductors (e.g., superconducting tapes, cables).

**File:** `../mesh/cl_Protoshell.hpp` (header-only, namespace `belfem`)

### Structure

```cpp
class Protoshell
{
    const id_t            mID ;
    string                mLabel ;
    Vector< id_t >        mSideSetIDs ;      // sidesets forming the shell surface
    Vector< id_t >        mTerminals ;       // terminal ids
    Vector< real >        mThicknesses ;     // layer thicknesses in m
    Cell< string >        mMaterials ;       // material name per layer
    Cell< mesh::Curve * > mSideCurves ;      // deleted by the mesh
    Cell< mesh::Curve * > mTerminalCurves ;  // deleted by the mesh
    bool                  mEdgeCoating ;     // 3D: create side connector walls
    real                  mEdgeCoatingWidth ;
    Vector< id_t >        mFlippedSideSets ; // sidesets whose facets are flipped
};
```

### Key Methods

```cpp
Protoshell( const id_t aID );   // the only constructor

Vector< id_t > & sidesets();
Vector< id_t > & terminals();
Vector< real > & thicknesses();
Cell< string > & materials();
Cell< mesh::Curve * > & side_curves();
Cell< mesh::Curve * > & terminal_curves();
string & label();
bool & edge_coating();
real & edge_coating_width();
Vector< id_t > & flipped_sidesets();
id_t id() const;
```

There are no `add_*` methods. Each accessor returns a mutable reference to the container,
which the caller fills directly. `MaxwellFactory::read_thin_shell_data()` builds the
protoshells from the input deck. `CutFactory` creates the side and terminal curves.

### Usage Pattern

```cpp
#include "cl_Protoshell.hpp"

Protoshell * tShell = new Protoshell( 1 );

tShell->sidesets()    = { 101, 102 };   // shell surface
tShell->terminals()   = { 201, 202 };
tShell->materials().push( "ybco" );
tShell->thicknesses() = { 1e-6 };        // 1 micron

// Pass to CutFactory (a non-const reference: the Cell must outlive the factory)
Cell< Protoshell * > tProtoshells = { tShell };
mesh::CutFactory tFactory( aMesh, &tTopology, tProtoshells,
                           mesh::CutAlgorithm::PellikkaGeneralized, false );
```

**⚠️ Memory Management:** CutFactory does **not** take ownership of Protoshell objects. Caller
must delete them after use (`MaxwellFactory` deletes its protoshells in its destructor).

---

## Cut Processing

### Purpose

`CutProcessor` converts cohomology generators (thick cuts) to FEM-compatible thin cuts.

**File:** `cl_CutProcessor.{hpp,cpp}`

**Key operations:**
1. Converts thick cuts (edge sets) to thin cuts (face sets)
2. Duplicates nodes on cut surface
3. Relinks elements to create separate DOFs for φ⁺ and φ⁻
4. Creates abstract nodes for jump conditions
5. Performs manifold cleanup

### Constructor

```cpp
CutProcessor::CutProcessor(
    Mesh*               aMesh,
    Cohomology*         aCohomology,
    const Vector<id_t>& aPhiBlocks,
    const Vector<id_t>& aNonPhiBlocks,
    const Vector<id_t>& aPhiInterfaces,
    const Vector<id_t>& aPhiBoundaries,
    const Vector<id_t>& aPhiPeriodic
);
```

**Parameters:**
- `aMesh`: Mesh to process
- `aCohomology`: Cohomology computation results
- `aPhiBlocks`: Block IDs for phi domains
- `aNonPhiBlocks`: Block IDs for non-phi domains
- `aPhiInterfaces`: Interface sideset IDs
- `aPhiBoundaries`: phi boundary sideset IDs
- `aPhiPeriodic`: phi periodic sideset IDs

`CutFactory::compute_thin_cuts_and_duplicate_interface_nodes()` passes the five
`Topology` lists in this order (`phi_block_ids()`, `non_phi_block_ids()`,
`phi_interface_ids()`, `phi_boundary_ids()`, `phi_periodic_ids()`).

### Key Methods

```cpp
Cell<Node*>& abstract_nodes();
```

Returns abstract nodes created for cut DOFs.

```cpp
id_t max_node_id();
id_t max_sideset_id();
```

Returns maximum IDs after node duplication and cut creation.

**Use case:** Ensures unique IDs when creating additional mesh entities.

```cpp
void save_debug_meshes();
```

Exports one debug mesh per cut for visualization, as `cut_<index>.vtk`
(`CutProcessor::save_debug_meshes` in `cl_CutProcessor.cpp`).

### Internal Workflow

The cut processing pipeline (called internally by CutFactory):

1. **Thick cut extraction**: Extract edges from cohomology generators
2. **Face identification**: Find faces incident to thick cut edges
3. **Manifold filtering**: Remove non-manifold artifacts (double pockets, etc.) — design of record; the step that actually runs is `Cohomology::remove_cut_pockets()`, see "Manifold Cleanup" below
4. **Node duplication**: Duplicate nodes on the positive side of the cut. Each cut set assigns
   copy ids while walking its originals in ascending id order, with the non-periodic originals
   first and the periodic pairs in the order encountered by that walk. The sets are processed
   in sorted order of their pattern hex strings. `collect_duplicates` then appends the abstract
   nodes and the copies in ascending copy-id order, so a copy's position never depends on the
   hash order of a duplicate map
5. **Element relinking**: Update element connectivity to duplicated nodes
6. **Abstract node creation**: Create DOFs for [φ] = φ⁺ - φ⁻
7. **SideSet creation**: Package cut as SideSet for FEM assembly

### Manifold Cleanup

> **Design note — design of record, not wired in.** `manifold_filter_3d()` and `check_surface()`
> do **not** exist anywhere in the tree, and no shipped code applies the three-phase scheme or the
> metric thresholds below. No Tarjan articulation-point implementation exists in the tree.
>
> **What runs today is `Cohomology::clean_spfa()` → `remove_cut_pockets()`**
> (`cl_Cohomology.hpp:135,145`), applied after SPFA rectification. Both are **private** — they run
> internally during cohomology computation and are not part of the caller-facing API. There is no
> `CutProcessor` method of either name.
>
> The description is kept because it may be the design of record; it is marked so that no reader
> goes looking for the functions or tunes the thresholds. **Its status is Gregory Giard's to
> confirm** (`doc/ai_collaboration_protocol.md` §7.1). The fuller note is in
> [cohomology_theory_and_implementation.md](cohomology_theory_and_implementation.md),
> "Post-Processing Cleanup".

**Problem:** Cohomology generators may produce non-manifold surfaces (edges shared by >2 faces).

**Solution as designed:** Multi-phase cleanup:
- **Phase 1**: Tarjan's algorithm to detect articulation points (necks)
- **Phase 2**: Region growing to detect embedded pockets
- **Phase 3**: BFS layering with metric validation

**Metrics as designed:**
- Cycle density > 0.3
- Compactness < 2.0
- Size < 20 faces

**Implementation as shipped:** `Cohomology::clean_spfa()` → `remove_cut_pockets( const bool aFireTierA )`
(`cl_Cohomology.hpp:135,145`). The designed filter above is not among the code paths that run.

---

## Data Structures

### Chain

**Purpose:** A chain is a formal sum of k-simplices with integer coefficients, keyed by
the mesh entity's *index* (`Node::index()`, `Edge::index()`, `Face::index()` for the
triangles of a 3D mesh, `Element::index()` for the top-dimensional cells), never by its id. It can also carry its own boundary and coboundary
chains, which the simplicial-complex builder fills alongside it.

**File:** `cl_Chain.{hpp,cpp}`

**Structure:**
```cpp
class Chain
{
    int                        mDim ;           // 0 nodes, 1 edges, 2 faces, 3 elements
    Mesh *                     mMesh ;
    OrderedMap< index_t, int > mSimplicesMap ;  // simplex index -> coefficient
    Chain *                    mBoundary ;      // (k-1)-chain, owned, or nullptr
    Chain *                    mCoboundary ;    // (k+1)-chain, owned, or nullptr
};
```

The constructor `Chain( aDim, aMesh, aIsBound, aIsCobound )` allocates `mBoundary` unless
`aDim == 0` or `aIsBound` is set, and allocates `mCoboundary` unless `aDim == 3` or
`aIsCobound` is set. The two flags mark a chain that *is* someone's boundary or
coboundary. A chain built with both flags false therefore owns a boundary (unless `aDim == 0`)
and a coboundary (unless `aDim == 3`). Each of those owns at most one chain in the opposite
direction (the boundary's coboundary and the coboundary's boundary), and those own nothing.
The destructor deletes both members.

**Key operations:** `operator+`, `operator-`, and `operator*` modify `this` in place.
`operator()` changes neither chain and returns the pairing. None of the operators returns a
new chain.
```cpp
void operator+( Chain * aChain ) ;        // this += aChain
void operator-( Chain * aChain ) ;        // this -= aChain
void operator*( const int aMult ) ;       // every coefficient *= aMult
int  operator()( Chain * aChain ) ;       // sum over shared indices of coeff * coeff

void add_simplex_to_chain( const index_t aIndex, const int aCoeff ) ;  // adds, erases a zero
void remove_simplex_from_chain( const index_t aIndex ) ;
void add_chain_to_chain( Chain * aChain, int aCoeff ) ;                // this += aCoeff * aChain
void remove_chain_from_chain( Chain * aChain ) ;

OrderedMap< index_t, int > & get_simplices_map() ;
int  get_coefficient( const index_t aIndex ) ;    // 0 when the index is absent
void set_coefficient( const index_t aIndex, const int aCoeff ) ;  // 0 erases the entry
int  get_dim() ;

Chain * get_boundary() ;
Chain * get_coboundary() ;
void add_simplex_to_boundary( index_t aIndex, const int aCoeff ) ;    // forwards to mBoundary
void add_simplex_to_coboundary( index_t aIndex, const int aCoeff ) ;  // forwards to mCoboundary
void print() ;
```

`add_chain_to_chain` refuses a dimension mismatch (it prints and returns). Otherwise, it also
propagates into `mBoundary` and `mCoboundary` wherever this chain owns them, so a generator
assembled from other chains keeps its boundary consistent. It adds the stored boundary and
coboundary of `aChain`; it does not recompute them from the coefficients. Therefore, `aChain`
must own every member that this chain owns: a chain built with a flag set has a null member in
that direction.
`remove_chain_from_chain` does not propagate.

**Usage:**
```cpp
Chain tC1( 1, tMesh, false, false ) ;
tC1.add_simplex_to_chain( tEdgeA->index(),  1 ) ;
tC1.add_simplex_to_chain( tEdgeB->index(), -1 ) ;

Chain tC2( 1, tMesh, false, false ) ;
tC2.add_simplex_to_chain( tEdgeA->index(), 2 ) ;

tC1 + &tC2 ;                       // tC1 is now { A: 3, B: -1 }
int tPairing = tC1( &tC2 ) ;       // 3 * 2 = 6

for ( const auto & [ tIndex, tCoeff ] : tC1.get_simplices_map() )
{
    // tIndex is Edge::index(); look the edge up through the mesh, not by id
}
```

---

### Cochain

**Purpose:** A `Cochain` is dual to `Chain`: a formal sum of k-cosimplices keyed by the
same simplex indices. The generators that `Cohomology` hands out are `Cochain *`.

**File:** `cl_Cochain.{hpp,cpp}`

**Structure:** The structure is identical to `Chain`, but uses `Cochain *` boundary and
coboundary members. The constructor and destructor follow the same rules.

**Key operations:** The in-place API corresponds to `Chain`, with `cochain` in the method
names, with one difference: `add_simplex_to_cochain` and `add_cochain_to_cochain` return
immediately for a zero coefficient. A zero-coefficient `add_cochain_to_cochain` therefore
neither checks the dimension nor propagates. The `Chain` versions do not return early.
```cpp
void operator+( Cochain * aCochain ) ;
void operator-( Cochain * aCochain ) ;
int  operator()( Chain * aChain ) ;       // evaluation <phi, c>: sum over shared indices

void add_simplex_to_cochain( const index_t aIndex, const int aCoeff ) ;
void remove_simplex_from_cochain( const index_t aIndex ) ;
void add_cochain_to_cochain( Cochain * aCochain, const int aCoeff ) ;
void remove_cochain_from_cochain( Cochain * aCochain ) ;

OrderedMap< index_t, int > & get_simplices_map() ;
int  get_coefficient( const index_t aIndex ) ;
void set_coefficient( const index_t aIndex, const int aCoeff ) ;
int  get_dim() ;

Cochain * get_coboundary() ;
Cochain * get_boundary() ;
void add_simplex_to_coboundary( index_t aIndex, const int aCoeff ) ;
void add_simplex_to_boundary( index_t aIndex, const int aCoeff ) ;
bool isBound() ;      // no callers in the tree
bool isCobound() ;
void print() ;
```

There is no scalar `operator*` on `Cochain`. `add_cochain_to_cochain` accumulates
(`this += aCoeff * aCochain`), so to scale a cochain, add it with the factor to a new, empty
cochain of the same dimension. `Cohomology` recombines generators this way into fresh
cochains.

**Usage:**
```cpp
Cochain tPhi( 1, tMesh, false, false ) ;
tPhi.add_simplex_to_cochain( tEdgeA->index(),  1 ) ;
tPhi.add_simplex_to_cochain( tEdgeB->index(), -1 ) ;

Chain * tC = ... ;
int tValue = tPhi( tC ) ;          // cochain evaluated on the chain
```

---

### CutData

**Purpose:** Stores topology and metadata for a single cut (one H^1 generator). `CutData` is
internal to `CutProcessor`, which creates one per generator.

**File:** `cl_CutData.{hpp,cpp}`

**Structure:**
```cpp
class CutData
{
    Mesh *            mMesh ;
    Cochain *         mCochain ;            // the generator
    const index_t     mIndex ;
    const id_t        mID ;
    const uint        mNumberOfCuts ;
    Cell< Edge * >    mCohomologyEdges ;
    DynamicBitset *   mCohomologyPlus ;     // edge index -> coefficient +1
    DynamicBitset *   mCohomologyMinus ;    // edge index -> coefficient -1
    Cell< Element * > mCutElements ;        // elements that contribute to this cut
    Vector< int >     mCutCases ;
    Map< id_t, int >  mCutCaseMap ;
    Cell< Edge * >    mThinCutEdges ;       // thin cut, 2D
    Cell< Face * >    mThinCutFaces ;       // thin cut, 3D
};
```

**Key methods:**
```cpp
CutData( Mesh * aMesh, Cohomology * aCohomology, const index_t aCohomologyIndex );

int  weight( const Edge * aEdge );                      // +1, -1 or 0, by edge index
int  cut_case( const index_t aLocalElementIndex );      // case of the i-th cut element
void collect_coefficients( const index_t aNumEdges );
void collect_elements( const Vector< id_t > & aNonPhiDomains );
void collect_edges();
void add_thin_cut_sidesets_to_mesh( id_t & aMaxSideSetID, id_t & aMaxElementID );

Cell< Element * > & elements();
Cell< Edge * > &    cohomology_edges();
Cell< Edge * > &    thin_cut_edges();
Cell< Face * > &    thin_cut_faces();
index_t             index() const;
Mesh *              create_debug_mesh();
void flag_nodes();
void flag_all_elements();
void unflag_all_elements();
```

The per-element cut case is computed by the private `determine_cut_case_2d()` and
`determine_cut_case_3d()`. There is no public node or edge membership test.

---

### SimplicialComplex

**Purpose:** Builds and reduces simplicial chain/cochain complexes.

**File:** `cl_SimplicialComplex.{hpp,cpp}`

**Structure:**
```cpp
class SimplicialComplex
{
    Cell< Map< index_t, Chain * > >   mChainsMap;     // k-chains
    Cell< Map< index_t, Cochain * > > mCochainsMap;   // k-cochains
    Mesh *                            mMesh;          // the mesh the complex is built from
};
```

**Key methods:**
```cpp
void create_complex(Mesh* aMesh, bool aPeriodicity);
```

Builds chain/cochain complex from flagged mesh entities.

**Flags:**
- Nodes: 0-simplices
- Edges: 1-simplices
- Faces: 2-simplices
- Elements: 3-simplices

```cpp
void reduce_complexPellikkaGeneralized();
void reduce_complexPellikka();
void reduce_complexCCR();
```

Reduction algorithms (see [Cut Algorithm Selection](#cut-algorithm-selection)).

```cpp
Cell<Matrix<int>> createMatrixFromBoundaryMap();
Cell<Matrix<int>> createMatrixFromCoboundaryMap();
```

Exports boundary/coboundary operators as integer matrices.

**Structure:** `∂_k : C_k → C_{k-1}`, `δ^k : C^k → C^{k+1}`

---

### Cut Sideset Creation

**Purpose:** Cohomology generators become mesh sidesets (the thin cuts).

**Live path:** `CutFactory::run()` calls
`compute_thin_cuts_and_duplicate_interface_nodes()`, which constructs a
`CutProcessor`. The `CutProcessor` constructor calls
`create_thin_cut_sidesets()`, which delegates to
`CutData::add_thin_cut_sidesets_to_mesh()`. Users do not call any of this
directly; it runs inside `CutFactory::run()`.

Each cut becomes one new sideset, added to the mesh with `Mesh::add_sideset()`. It is labelled
`cut_` followed by a zero-padded number, and it keeps domain type `Default`.
`CutProcessor`/`CutData` is the only sideset path.

---

### InterfaceProcessor

**Purpose:** Handles node duplication and element relinking at domain interfaces (e.g., air-conductor boundaries).

**File:** `cl_InterfaceProcessor.{hpp,cpp}`

**Structure:**
```cpp
class InterfaceProcessor
{
    Mesh* mMesh;
    Topology* mTopology;
    Cell<Node*>& mAbstractNodes;          // Output: abstract nodes
    Map<id_t, InterfaceSet*> mSetsMap;    // Interface sets by sideset ID
    DynamicBitset* mIsAirBlock;           // Air block flags
    DynamicBitset* mIsFerroBlock;         // Ferromagnetic block flags
};
```

**Key methods:**
```cpp
InterfaceProcessor(
    Mesh* aMesh,
    Topology* aTopology,
    Cell<Node*>& aAbstractNodes,
    const uint aNumOriginalSideSets,
    const id_t aMaxNodeID
);
```

Constructor performs full interface processing pipeline:
1. Creates interface sets (groups of related interfaces)
2. Connects sidesets to sets
3. Duplicates nodes at interfaces
4. Relinks elements to duplicated nodes. A node is copied once per interface set, except on the
   conductor side. Each tie-weight set whose master block is a conductor takes its copy from one
   shared map keyed by the node's own id. A final sweep installs that copy in every element of
   every conductor block that still holds the original. Two conductor blocks that meet the air
   along a common line therefore reference the same copy, and the conductor edges leaving that
   line along their common face are created once. Without this sharing, each block had its own
   copy, so the edges were doubled and untied, and tangential h could jump across the air-adjacent
   row of the conductor|conductor face: a lossless sheet. Ferro, buffer and air-master sets keep
   private per-set copies; a ferro|ferro sideset stays split, as its material jump requires.
   On a periodic mesh, the copies of a periodic node pair are also tied to each other, either
   within one interface set or across two; this pairing runs at the end of node duplication,
   before the relink. Across sets, a pair is tied only when both sets have the same treatment
   and, for tie-weight sets, the same master domain type (`pair_within`). When no compatible set
   holds the partner's copy, because it lies only on a buffer|air, ferro|air, or air-master set
   against a conductor|air set, or only on a decoupled coil set against a tie-weight set, setup
   refuses the pair with a message naming both nodes and both domain types or treatments. The
   reason is the tie itself: it identifies like unknowns, h on conductor edges and φ elsewhere
   (Messe et al. 2023, §1), and cannot represent a change of formulation across the plane. The
   rule is applied to the set's master domain type, so it also refuses ferro against buffer,
   although both are φ blocks. A periodic plane that coincides with a change of domain type, or
   that a coil touches on one side only, can be avoided by shifting the plane to a section where
   the interface has one domain type and the same treatment on both sides. The plain
   "not duplicated in any interface set" refusal remains for a partner that lies on no interface
   set at all.
   Two rules come first. First, every interface sideset, that is, every sideset with a facet
   this processor duplicates at (one side air or ferro, not both air), must have one master
   domain type on all of its facets that have a slave; a sideset without such a facet is left
   alone. A sideset is admitted to an interface set, and the set is classed as shared or
   private and paired across a periodic plane, by its first facet, as in 0.9.2: a native Gmsh
   mesh builds one sideset per geometric entity, so its facets cannot differ in kind.
   `refuse_mixed_master_types` checks, after `CutProcessor` has run, every sideset that
   existed before it (generated cut sidesets are not examined, and thin-shell sidesets have
   no slave by then), before any sideset is admitted and before any interface copy exists. If
   it finds a mismatch, it reports the sideset id, the two domain types (name and enumerator
   value, since several types print under one name), the two facets with their master
   elements, and the remedy: split the sideset so that each piece has one master domain type.
   A conductor|air facet beside an air|air facet is refused the same way, since the air|air
   facet's master is air. A sideset over several blocks of one domain type stays accepted. A
   reused `.bfm` (one whose checksum and configuration tag match the deck) skips all of this:
   it replays the enrichment it was written with.
   Second, interface-copy numbering is a function of the mesh alone. Copies are allocated one
   set at a time, in set order (the sorted hex strings of the block-index bitsets of each
   set's block pair, a fixed function of the mesh's block order), and within each set in
   ascending order of the original id. A shared conductor copy is allocated by the first set
   that reaches the original. Copies are appended to the mesh in ascending copy-id order,
   which is the allocation order. Periodic pairs are backed up set by set: first, all pairs
   tied within one set, in ascending order of the original id; then the next set; finally,
   pairs tied across sets, again set by set and in ascending order. The maps that hold
   originals and copies are hash maps, but their iteration order never reaches the mesh.
   Therefore, interface-copy ids and positions are the same on every platform for the same
   incoming node order. The cut copies, which precede the interface copies in the mesh, are
   appended the same way: their ids come from a sorted walk of each cut set's originals (see
   "CutSet"), and `CutProcessor::collect_duplicates` sorts the kept copies by id after the
   abstract nodes before appending them. Therefore, no hash-map iteration order reaches a
   copy's position. Cut copies and interface copies occupy positions determined by the node
   order entering the cut processor. That is not yet a platform-independent node numbering.
   `CutFactory` renumbers the mesh with reverse Cuthill-McKee before the cohomology, and
   `MaxwellFactory` renumbers it again after enrichment, in the order that reaches the `.bfm`
   and the dof manager. `symrcm` sorts a node's neighbors by degree with `std::sort`, so equal
   degrees follow whatever order the C++ standard library's sort happens to leave them in. That
   behavior is reproducible on one toolchain but is not guaranteed across toolchains until the
   graph sort has a tie-break. The thin-shell copies were already appended in ascending index
   order of their originals. This does not cover a reused `.bfm`, which replays the numbering
   with which it was written, the dof manager's own ordering, or partitioned runs.
5. Ties each duplicate to its original (`InterfaceTreatment::TieWeight1`) or decouples it when a coil touches the interface (`InterfaceTreatment::Decouple`)

**Helper class: InterfaceSet**
```cpp
class InterfaceSet
{
    Mesh* mMesh;
    DynamicBitset* mBitset;               // Interface pattern
    Cell<Element*> mElements;              // Elements on interface
    Map<id_t, Node*> mOriginals;          // Original nodes
    Map<id_t, Node*> mDuplicates;         // Duplicated nodes
    Cell<SideSet*> mSideSets;             // Sidesets in this set
};
```

**Physical context:** Duplicates the nodes of every sideset that separates a φ block (Air/Buffer/Ferro) from a non-φ block or Ferro (conductor-air, conductor-ferro, ferro-air), so the φ side and the other side carry distinct DOFs; coil-touching interfaces are decoupled (`InterfaceTreatment::Decouple`). It does not create abstract nodes; it receives the factory's list only to keep them out of the duplicate scan.

**Usage:**
```cpp
// Internal use by CutFactory
Cell<Node*> tAbstractNodes;
mesh::InterfaceProcessor tProcessor(
    aMesh, &tTopology, tAbstractNodes,
    originalSidesetCount, maxNodeID
);
// tProcessor constructor performs all processing
// tAbstractNodes is only read (excluded from the duplicate scan)
```

**Note:** Runs unconditionally from `CutFactory::compute_thin_cuts_and_duplicate_interface_nodes()`, thin shells or not.

---

### CutSet

**Purpose:** Manages node duplication for a single cohomology cut.

**File:** `cl_CutSet.{hpp,cpp}`

**Structure:**
```cpp
class CutSet
{
    Mesh* mMesh;
    Cell<Node*>& mNodeOriginals;          // Candidate nodes for duplication
    DynamicBitset* mBitset;               // cut pattern: bit c set ⇔ this set's nodes are split by cut c
    DynamicBitset* mNodeBitset;           // Node switches (which side of cut)
    Map<id_t, Node*> mNodeDuplicates;     // Original → duplicate mapping
};
```

**Key methods:**
```cpp
CutSet(
    Mesh* aMesh,
    Cell<Node*>& aNodeOriginals,
    const string& aHexString,          // Bitset pattern as hex
    const index_t aNumberOfCuts
);
```

Constructor parses hex-encoded bitset describing cut topology.

```cpp
bool test(const Node* aNode) const;
```

Returns `mBitset->test( aNode->index() )`. `mBitset` is the set's cut-pattern bitset, one bit
per cut (`aNumberOfCuts` bits), so the call reads bit number *node index* of the pattern. It
is not a node-on-cut membership test. No code in `src/` calls it.

```cpp
void create_duplicates(id_t& aMaxNodeID, Cell<Node*>& aAbstractNodes);
```

Creates duplicate nodes for positive side of cut.

**Algorithm:**
1. Iterate through nodes on the cut in ascending id (`CutProcessor::collect_nodes` sorts the
   originals): the non-periodic originals first, then the periodic pairs in the order
   encountered by the walk
2. Create a new Node with a unique ID in that order
3. Store it in the `mNodeDuplicates` map (a hash map; its iteration order never reaches the
   mesh, because `CutProcessor::collect_duplicates` sorts the copies by id before appending
   them)
4. Tie each duplicate to its original plus the abstract nodes of the cuts in this set's pattern (weight 1 each; the abstract nodes are created by `CutProcessor::create_abstract_nodes()`)

```cpp
Node* duplicate(Node* aNode);
```

Retrieves previously created duplicate for a node.

**Returns:** Duplicate node pointer (must exist in map).

```cpp
DynamicBitset* node_bitset();
Map<id_t, Node*>& duplicate_map();
```

Accessors for internal data structures.

**Usage:**
```cpp
// Internal use by CutProcessor
mesh::CutSet tCutSet(
    aMesh, nodeOriginals, hexBitset, numCuts
);

// Create duplicates
tCutSet.create_duplicates(maxNodeID, abstractNodes);

// Access duplicate
mesh::Node* dup = tCutSet.duplicate(originalNode);
```

**Bitset encoding:** Uses DynamicBitset with hex string serialization for efficient cut pattern storage and transmission.

**Note:** This is an internal class used by CutProcessor during thick-to-thin conversion.

---

## Algorithm Classes and Functions

### BeltedTree

**Purpose:** Alternative algorithm for cohomology computation using spanning tree construction with belt identification.

**File:** `cl_BeltedTree.{hpp,cpp}`

**Structure:**
```cpp
class BeltedTree
{
    Cell<Chain*> m1HomologyGenerators;       // H_1 generators (input)
    Cell<Cochain*> m1CohomologyGenerators;   // H^1 generators (output)
    Cell<Edge*> mBeltFasteners;              // Edges connecting tree to belts
    Cell<index_t> mTree;                     // Spanning tree edge indices
    SimplicialComplex* mSimplicialComplex;   // Reduced complex
    Mesh* mMesh;                             // Full mesh
};
```

**Algorithm concept:**

The belted tree method constructs cohomology generators by:
1. Building a spanning tree of the simplicial complex
2. Identifying "belts" (cycles not in the tree) corresponding to homology generators
3. Selecting "belt fasteners" (edges connecting tree to belts)
4. Computing dual cochains from the fastener configuration

**Key methods:**
```cpp
BeltedTree(
    Mesh* aMesh,
    SimplicialComplex* aSimplicialComplex,
    Cell<Chain*> a1HomologyGenerators
);
```

Constructor requires pre-computed H_1 homology generators.

**Precondition:** Homology must be computed first.

```cpp
void select_belt_fasteners();
```

Identifies edges that connect the spanning tree to independent cycles.

**Algorithm:** For each homology generator (1-cycle), find minimal edge set that completes the cycle when added to the tree.

```cpp
void create_tree();
```

Constructs spanning tree of the 1-skeleton (edge graph).

**Method:** Breadth-first search or depth-first search from arbitrary root node.

```cpp
void compute_cohomology();
```

Computes H^1 cochains from belt fastener configuration.

**Duality:** Each belt fastener corresponds to a cohomology generator via Poincaré duality.

```cpp
Cell<Cochain*>& get_cohomology();
```

Returns computed cohomology generators.

**Visualization:**
```cpp
void create_tree_field( Mesh * aEdgeMesh, const string aFieldName );
void create_cohomology_field( Mesh * aEdgeMesh );
```

Creates mesh fields to visualize spanning tree and cohomology generators.

**Usage pattern:**
```cpp
// tComplex is a SimplicialComplex* built on the flagged mesh
// 1. Compute homology first (the constructor computes the groups)
mesh::Homology tHomology(tComplex, aFullMesh);
Cell<Chain*>& h1_gens = tHomology.get_Generators()( 1 );

// 2. Create BeltedTree
mesh::BeltedTree tTree(aMesh, tComplex, h1_gens);

// 3. Execute algorithm
tTree.select_belt_fasteners();
tTree.create_tree();
tTree.compute_cohomology();

// 4. Retrieve cohomology
Cell<Cochain*>& cohom_gens = tTree.get_cohomology();
```

**When to use:**
- Alternative to reduction-based methods (Pellikka, CCR)
- Useful when homology generators have specific physical meaning
- Educational: Explicit geometric interpretation of duality

**Comparison to other algorithms:**
- **Advantage:** Direct geometric construction (easy to visualize)
- **Disadvantage:** Requires pre-computed homology (extra step)
- **Performance:** unmeasured (no timing comparison recorded)

**Selection in CutFactory:**
```cpp
mesh::CutFactory tFactory( aMesh, &tTopology, tProtoshells,
                           mesh::CutAlgorithm::BeltedTree, false );
```

On the factory path, the belted tree takes the suggested H_1 generators from the terminal
`Homology` (`set_terminals()`), not a homology computed from a complex. Without the suggested
homology, `compute_cohomologies()` stops with `BELFEM_ERROR`.

---

### fn_Smith - Smith Normal Form Functions

**Purpose:** Integer matrix decomposition for computing homology/cohomology groups.

**File:** `fn_Smith.{hpp,cpp}`

**Mathematical background:**

The Smith Normal Form decomposes an integer matrix A into:
```
Q⁻¹ * A * R = D        equivalently  A = Q * D * R⁻¹
```
Where (Kaczynski et al. 2004; `smithForm` returns Q, Q⁻¹, R, R⁻¹ as `Q`, `Q_`, `R`, `R_`):
- Q, R are unimodular matrices (det = ±1)
- D is diagonal matrix with d_1 | d_2 | ... | d_r (divisibility chain)
- Diagonal entries d_i are **invariant factors**
- Number of nonzero entries r is the **rank**

**Homology application:**

Given boundary operators ∂_k : C_k → C_{k-1} (rank r_k) and ∂_{k+1} (rank r_{k+1}):
- dim Ker(∂_k) = n_k - r_k (n_k = number of k-simplices = columns of ∂_k)
- dim Im(∂_{k+1}) = r_{k+1}
- H_k rank (Betti number) β_k = n_k - r_k - r_{k+1}
- Torsion of H_k = the invariant factors d_i > 1 of ∂_{k+1}

**Core functions:**

```cpp
std::tuple<Matrix<int>, Matrix<int>, uint>
rowEchelon(Matrix<int>& aMat);
```

Computes row echelon form of integer matrix.

**Returns:**
- Q: Row transformation matrix
- Q_: Inverse of Q
- rank: Number of nonzero rows

**Algorithm:** Euclidean-style elimination in exact integer arithmetic: pivot on the smallest non-zero entry, subtract integer quotients of rows (`floor(a/b)`), swap; no fractions ever appear.

```cpp
std::tuple<Matrix<int>, Matrix<int>>
kernelImage(Matrix<int>& aMat);
```

Computes kernel and image of integer matrix.

**Returns:**
- Kernel basis (columns span Ker(A))
- Image basis (columns span Im(A))

**Method:** Row echelon form followed by basis extraction.

```cpp
std::tuple<Matrix<int>, Matrix<int>, Matrix<int>, Matrix<int>, uint, uint>
smithForm(Matrix<int>& aMat);
```

**Main Smith Normal Form decomposition.**

`aMat` is overwritten with the diagonal matrix D.

**Returns** `( Q, Q_, R, R_, s, t )`, with Q⁻¹ · A · R = D:
- Q: Left transformation matrix
- Q_: Inverse of Q
- R: Right transformation matrix
- R_: Inverse of R
- s: number of diagonal entries equal to 1 (unit invariant factors)
- t: number of nonzero diagonal entries, i.e. the rank of A

The torsion coefficients are the diagonal entries `D(i,i)` for `s <= i < t`. `Homology` and
`Cohomology` read the generators from columns `s` onward of their basis matrix.

**Algorithm:**
1. Row reduction to echelon form
2. Column reduction to eliminate above-diagonal entries
3. Iterative divisibility checking (ensure d_i | d_{i+1})
4. Minimize diagonal entries via elementary operations

**Helper functions:**

```cpp
std::pair<uint, uint> minNonzero(Matrix<int>& aMat, const uint k);
```

Finds smallest nonzero entry in submatrix (optimization for faster reduction).

```cpp
void moveMinNonzero(Matrix<int>& aMat, Matrix<int>& aQ, Matrix<int>& aQ_,
                    Matrix<int>& aR, Matrix<int>& aR_, const uint k);
```

Swaps rows/columns to move minimal element to pivot position.

```cpp
std::tuple<bool, uint, uint, int>
checkForDivisibility(Matrix<int>& aMat, const uint k);
```

Checks whether the pivot `B[k,k]` divides every entry of the trailing submatrix `B[k+1:end, k+1:end]` (1-based; Kaczynski et al. 2004, Smith-form step); returns the first offending `(row, col, quotient)`.

**Returns:** (divisible, row_index, col_index, quotient): the 1-based position of the first
entry the pivot does not divide, and the truncated integer quotient of that entry by the
pivot. A divisible matrix returns `(true, 0, 0, 0)`.

```cpp
void partSmithForm(Matrix<int>& aMat, Matrix<int>& aQ, Matrix<int>& aQ_,
                   Matrix<int>& aR, Matrix<int>& aR_, const uint k);
```

Performs partial Smith form reduction for diagonal entry k.

**Method:** Euclidean algorithm to eliminate non-divisible entries.

```cpp
Matrix<int> SolveInt(Matrix<int> aMat, Matrix<int>& aVec);
```

Solves integer linear system A*x = b using Smith decomposition.

**Returns:** Integer solution vector (if exists).

**Usage in homology:**

```cpp
// Compute boundary matrix
Cell<Matrix<int>> boundary = tComplex->createMatrixFromBoundaryMap();

// Smith decomposition of ∂_k (overwrites boundary( k ) with D)
auto [ tQ, tQ_, tR, tR_, tUnitPivots, tRank ] = smithForm( boundary( k ) );

// Ker(∂_k) = span of the columns tRank ... n_k-1 of tR
// Betti_k  = ( n_k - tRank ) - rank(∂_{k+1})
```

**Performance notes:**

- **Complexity:** O(n³) for dense matrices (n = dimension)
- **Optimization:** Uses smallest-entry pivoting for faster convergence
- **Integer overflow:** No protection for very large matrices (TODO: arbitrary-precision)
- **Sparsity:** Does not exploit sparsity (future improvement)

**Elementary operations:**

The following operations preserve equivalence and update transformation matrices. Each takes
the matrix as its first argument, and the row and column indices are **1-based**:
- `rowExchange(i, j)`: Swap rows i and j
- `columnExchange(i, j)`: Swap columns i and j
- `rowMultiply(i)`: Multiply row i by -1
- `columnMultiply(i)`: Multiply column i by -1
- `rowAdd(i, j, q)`: Add q times row j to row i
- `columnAdd(i, j, q)`: Add q times column i to column j

**Corresponding operations** (with transformation tracking):
- `rowExchangeOperation()`: Updates Q, Q_, and matrix
- `columnExchangeOperation()`: Updates R, R_, and matrix
- `rowAddOperation()`: Updates Q, Q_, and matrix
- `columnAddOperation()`: Updates R, R_, and matrix

**Example:**

```cpp
#include "fn_Smith.hpp"

// aDk  = ∂_k     : rows = (k-1)-simplices, columns = k-simplices
// aDk1 = ∂_{k+1} : rows = k-simplices,     columns = (k+1)-simplices
void betti_number( Matrix< int > & aDk, Matrix< int > & aDk1 )
{
    const uint tNk = aDk.n_cols();   // number of k-simplices

    // both calls overwrite their input with D
    auto [ tQ,  tQ_,  tR,  tR_,  tUnitPivots,  tRankK  ] = smithForm( aDk );
    auto [ tQ1, tQ1_, tR1, tR1_, tUnitPivots1, tRankK1 ] = smithForm( aDk1 );

    const uint tBetti   = tNk - tRankK - tRankK1;    // free rank of H_k
    const uint tTorsion = tRankK1 - tUnitPivots1;    // invariant factors > 1

    message( InfoLevel::Default, "beta_k = %u, torsion factors: %u\n", tBetti, tTorsion );

    // Ker(∂_k): the last tNk - tRankK columns of tR
    Matrix< int > tKernel( tR.n_rows(), tNk - tRankK );
    for ( uint j = tRankK; j < tNk; ++j )
    {
        for ( uint i = 0; i < tR.n_rows(); ++i )
        {
            tKernel( i, j - tRankK ) = tR( i, j );
        }
    }
}
```

**Thread safety:** ⚠️ **Not thread-safe on shared inputs:** modifies the input matrix in place.

---

## Common Patterns

### Pattern 1: Basic Cut Generation

```cpp
#include "cl_CutFactory.hpp"
#include "cl_Topology.hpp"

void generate_cuts(Mesh* aMesh)
{
    // 1. Analyze topology
    mesh::Topology tTopology(aMesh);
    tTopology.run();

    tTopology.synchronize_maps();   // collective

    // 2. Create empty protoshell list (no thin shells)
    Cell< Protoshell * > tProtoshells;

    // 3. Create factory
    mesh::CutFactory tFactory( aMesh, &tTopology, tProtoshells,
                               mesh::CutAlgorithm::PellikkaGeneralized,
                               false );
    tFactory.set_periodicity( nullptr );

    // 4. Run pipeline (every rank)
    tFactory.run();

    // 5. Process results (master): the cut sidesets are labelled cut_<n>
    if ( comm_rank() == 0 )
    {
        for ( mesh::SideSet * tSideSet : aMesh->sidesets() )
        {
            if ( tSideSet->label().rfind( "cut_", 0 ) == 0 )
            {
                // leave the domain type at Default (see Quick-Start)
                message( InfoLevel::Default, "Cut %lu: %lu facets\n",
                         ( long unsigned int ) tSideSet->id(),
                         ( long unsigned int ) tSideSet->number_of_facets() );
            }
        }
    }
}
```

---

### Pattern 2: Terminal-Driven Cuts

```cpp
void generate_terminal_cuts( Mesh * aMesh )
{
    // Topology analysis
    mesh::Topology tTopology( aMesh );
    tTopology.run();
    tTopology.synchronize_maps();

    Cell< Protoshell * > tProtoshells ;
    mesh::CutFactory tFactory( aMesh, &tTopology, tProtoshells,
                               mesh::CutAlgorithm::PellikkaGeneralized,
                               false );
    tFactory.set_periodicity( nullptr );

    if ( comm_rank() == 0 )
    {
        // 3D, two bulk conductors: entry = { input sideset, output sideset }
        Cell< Cell< id_t > > tTerminals = { { 101, 102 }, { 201, 202 } };

        // positions of thin-shell entries in tTerminals: none
        Cell< id_t > tThinShellIndices ;

        tFactory.set_terminals( tTerminals, tThinShellIndices );
    }

    tFactory.run();

    if ( comm_rank() == 0 )
    {
        message( InfoLevel::Default, "Generated %lu cuts\n",
                 ( long unsigned int ) aMesh->abstract_nodes().size() );
    }
}
```

---

### Pattern 3: Direct Cohomology Computation

```cpp
void compute_cohomology_groups(Mesh* aMesh)
{
    // Build and reduce the simplicial complex on the flagged entities
    mesh::SimplicialComplex tComplex(aMesh, true);
    tComplex.coreduce_complexPellikkaGeneralized();

    // Compute cohomology (the constructor computes and cleans the generators)
    mesh::Cohomology tCohomology(&tComplex, aMesh);

    // Extract H^1 generators
    Cell<Cochain*>& generators = tCohomology.get_Generators()( 1 );

    // Visualize
    tCohomology.create_kGeneratorsField(1, aMesh, "H1");

    // Analyze generators
    for (uint i = 0; i < generators.size(); ++i) {
        Cochain* gen = generators(i);
        message(InfoLevel::Default, "Generator %u: %lu edges",
                i, gen->get_simplices_map().size());
    }
}
```

---

### Pattern 4: Algorithm Performance Comparison

`run()` mutates the mesh (see "Mesh Mutation Contract"). The algorithms therefore cannot be
compared by running several factories in a loop over one mesh: every algorithm after the first
would see an already cut mesh. For a fair comparison, load a fresh mesh and build a fresh
`Topology` for each algorithm, then time `run()` with a `Timer`. `run()` itself logs the time of its
Smith normal form and thin-cut steps at `InfoLevel::Detailed`. No comparison of this kind is
recorded (see README, "Cut Algorithms").

---

### Pattern 5: Thin Shell with Terminals

The sidesets 101 and 102 must be `ThinShell` sidesets in the topology (their domain type
decides that, not the protoshell). The ids in the terminal list must match the curves
that `create_thin_shell_cuts()` builds; `MaxwellFactory` takes both from the deck. The
example shows the call order and the container types.

```cpp
void process_thin_shell_coil( Mesh * aMesh )
{
    mesh::Topology tTopology( aMesh );
    tTopology.run();
    tTopology.synchronize_maps();

    // Protoshell configuration (filled through its accessors)
    Protoshell * tShell = new Protoshell( 1 );
    tShell->sidesets()    = { 101, 102 };
    tShell->terminals()   = { 201, 202 };
    tShell->materials().push( "ybco" );
    tShell->thicknesses() = { 1e-6 };

    Cell< Protoshell * > tProtoshells = { tShell };

    mesh::CutFactory tFactory( aMesh, &tTopology, tProtoshells,
                               mesh::CutAlgorithm::PellikkaGeneralized,
                               false );
    tFactory.set_periodicity( nullptr );

    if ( comm_rank() == 0 )
    {
        // thin-shell preprocessing: run() does not do it
        tFactory.create_thin_shell_cuts();

        // entry 0: a bulk conductor (sideset ids); entry 1: the thin shell (curve ids)
        Cell< Cell< id_t > > tTerminals = { { 301, 302 }, { 7, 8 } };
        Cell< id_t > tThinShellIndices = { 1 };   // position of the thin-shell entry
        tFactory.set_terminals( tTerminals, tThinShellIndices );
    }

    tFactory.run();

    // the factory holds a reference to tProtoshells: delete after the factory is done
    delete tShell;
}
```

---

## Troubleshooting

### Problem: Non-Manifold Cut Surfaces

**Symptom:** Error message "Cut surface is non-manifold" or FEM assembly fails.

**Cause:** Cohomology generator contains "double pockets" or other topological defects.

**Solution:**
1. **Try different algorithm:**
   ```cpp
   // PellikkaGeneralized may create complex cuts
   // Try original Pellikka for cleaner results
   mesh::CutFactory tFactory( aMesh, &tTopology, tProtoshells,
                              mesh::CutAlgorithm::Pellikka, false );
   ```

2. **Enable debug output:**
   ```cpp
   tFactory.run();
   // save_debug_meshes() is a public member of CutProcessor, not of CutFactory.
   // CutFactory calls it in debug builds only, inside
   // compute_thin_cuts_and_duplicate_interface_nodes() (compiled out under NDEBUG);
   // it writes one file per cut: cut_<index>.vtk
   ```
   For curve meshes, `CutFactory::save_curve_debug_meshes()` is public and writes
   `curve_<id>.vtk`. Inspect either family in ParaView.

3. **Inspect the cleanup that actually runs:**
   - The shipped path is `Cohomology::clean_spfa()` → `remove_cut_pockets()`
     (`cl_Cohomology.hpp:135,145`), not the metric-threshold filter described under
     "Manifold Cleanup" — that scheme is a design of record with no implementation, so there
     are no cycle-density or size thresholds to tune. See the design note there before
     going looking for them.

---

### Problem: Incorrect Number of Cuts

**Symptom:** Expected N cuts, got M cuts (M ≠ N).

**Cause:** Mesh topology doesn't match expected domain connectivity.

**Diagnosis:**
```cpp
mesh::Topology tTopology(aMesh);
tTopology.run();

// on rank 0
message( InfoLevel::Default, "Phi blocks: %lu\n",
         ( long unsigned int ) tTopology.phi_block_ids().length() );
message( InfoLevel::Default, "Interfaces: %lu\n",
         ( long unsigned int ) tTopology.phi_interface_ids().length() );
```

**Solution:**
- **Check block domain types:** Ensure phi/non-phi classification is correct
- **Verify mesh connectivity:** Use mesh viewer to confirm expected topology
- **Terminal hints:** Provide explicit terminals via `set_terminals()` to guide computation

---

### Problem: Slow Cohomology Computation

**Symptom:** `run()` takes minutes/hours for moderate mesh sizes.

**Cause:** Using suboptimal algorithm (CCR or original Pellikka).

**Solution:**
```cpp
// Switch to generalized Pellikka
mesh::CutFactory tFactory( aMesh, &tTopology, tProtoshells,
                           mesh::CutAlgorithm::PellikkaGeneralized, false );
```

If still slow:
- **Profile:** Use `Profiler` to identify hotspot
- **Mesh quality:** Simplify mesh if possible (remove small features)
- **MPI:** The cohomology runs on the master rank only; more ranks do not speed it up

---

### Problem: Abstract Nodes Not Created

**Symptom:** `tFactory.abstract_nodes().size() == 0` after `run()`.

**Cause:** No cuts were generated (topology is simply-connected).

**Diagnosis:** count the `cut_<n>` sidesets in `Mesh::sidesets()` on the master (see "Result
Retrieval").

**Solution:**
- **Check mesh:** Verify domain is multiply-connected (e.g., toroid, solenoid)
- **Topology:** Ensure phi/non-phi domains are correctly identified
- **Terminals:** Cuts may not be needed if terminals span all connections

---

### Problem: Jump Condition Not Satisfied

**Symptom:** FEM solution violates Ampere's law (∮ H·dl ≠ I).

**Cause:** Cut orientation or DOF assignment incorrect.

**Diagnosis:**
1. **Visualize cut:**
   ```cpp
   tCohomology.create_kGeneratorsField(1, aMesh, "H1_gen");
   // Check edge orientations in ParaView
   ```

2. **Check abstract node count:** `CutProcessor` creates one abstract node and one `cut_<n>`
   sideset per H^1 generator. Compare `tFactory.abstract_nodes().size()` with the number of
   `cut_<n>` sidesets in `Mesh::sidesets()`.

**Solution:**
- **Orientation:** Use `Homology::generators_orientation()` to ensure consistent direction
- **Static condensation:** Verify transformation matrix T incorporates [φ] = I correctly
- **Gauge:** Confirm φ is gauged (set to zero at reference point)

---

## Usage Examples

### Example 1: Toroidal Magnet

```cpp
#include "cl_CutFactory.hpp"
#include "cl_Topology.hpp"

// called on every rank, after MPI initialization
void toroidal_magnet()
{
    // every rank: the master reads the file, the others get an empty mesh
    Mesh * tMesh = new Mesh( "toroid.exo" );

    if ( comm_rank() == 0 )
    {
        tMesh->block( 1 )->set_domain_type( DomainType::Air );
        tMesh->block( 2 )->set_domain_type( DomainType::Conductor );
    }

    mesh::Topology tTopology( tMesh );
    tTopology.run();
    tTopology.synchronize_maps();

    // Create factory (no thin shells)
    Cell< Protoshell * > tProtoshells ;
    mesh::CutFactory tFactory( tMesh, &tTopology, tProtoshells,
                               mesh::CutAlgorithm::PellikkaGeneralized,
                               false );
    tFactory.set_periodicity( nullptr );

    // Generate cuts
    Timer tTimer;
    tFactory.run();
    const uint64_t tElapsed = tTimer.stop();

    if ( comm_rank() == 0 )
    {
        message( InfoLevel::Default, "Generated %lu cuts in %lu ms\n",
                 ( long unsigned int ) tMesh->abstract_nodes().size(),
                 ( long unsigned int ) tElapsed );
    }

    // FEM: hand tMesh->abstract_nodes() to the equation (IWG::set_abstract_nodes),
    // then build the Kernel, which partitions and distributes the mesh.
    // MaxwellFactory::create_magnetic_kernel() is the reference for this sequence.

    delete tMesh;
}
```

---

### Example 2: Superconducting Cable

```cpp
#include "cl_CutFactory.hpp"
#include "cl_Protoshell.hpp"
#include "cl_Topology.hpp"

// called on every rank; block domain types are assumed set on the master,
// and sideset 101 is a ThinShell sideset
void superconducting_cable( Mesh * aMesh )
{
    mesh::Topology tTopology( aMesh );
    tTopology.run();
    tTopology.synchronize_maps();

    // Protoshell for thin shell approximation
    Protoshell * tShell = new Protoshell( 1 );
    tShell->sidesets()    = { 101 };      // cable surface
    tShell->terminals()   = { 201, 202 };
    tShell->materials().push( "bscco" );
    tShell->thicknesses() = { 0.2e-3 };   // 0.2 mm

    Cell< Protoshell * > tProtoshells = { tShell };

    mesh::CutFactory tFactory( aMesh, &tTopology, tProtoshells,
                               mesh::CutAlgorithm::PellikkaGeneralized,
                               false );
    tFactory.set_periodicity( nullptr );

    if ( comm_rank() == 0 )
    {
        tFactory.create_thin_shell_cuts();

        // one thin-shell entry: input curve 7, output curve 8
        Cell< Cell< id_t > > tTerminals = { { 7, 8 } };
        Cell< id_t > tThinShellIndices = { 0 };
        tFactory.set_terminals( tTerminals, tThinShellIndices );
    }

    tFactory.run();

    if ( comm_rank() == 0 )
    {
        message( InfoLevel::Default, "Abstract DOFs (= cuts): %lu\n",
                 ( long unsigned int ) aMesh->abstract_nodes().size() );
    }

    delete tShell;
}
```

---

### Example 3: Performance Analysis

```cpp
#include "cl_Cohomology.hpp"
#include "cl_SimplicialComplex.hpp"

void analyze_cohomology_performance(Mesh* aMesh)
{
    // Build simplicial complex on the flagged entities of aMesh
    // (the constructor builds the complex itself)
    Timer t1;
    mesh::SimplicialComplex* tComplex = new mesh::SimplicialComplex(aMesh, /*aPeriodicity=*/false);
    uint64_t time_build = t1.stop();

    message(InfoLevel::Default, "Complex construction: %lu ms", time_build);
    message(InfoLevel::Default, "  0-simplices: %u", tComplex->number_of_ksimplices(0));
    message(InfoLevel::Default, "  1-simplices: %u", tComplex->number_of_ksimplices(1));
    message(InfoLevel::Default, "  2-simplices: %u", tComplex->number_of_ksimplices(2));
    message(InfoLevel::Default, "  3-simplices: %u", tComplex->number_of_ksimplices(3));

    // Reduction
    Timer t2;
    tComplex->reduce_complexPellikkaGeneralized();
    uint64_t time_reduce = t2.stop();

    message(InfoLevel::Default, "Reduction: %lu ms", time_reduce);
    message(InfoLevel::Default, "  0-simplices: %u", tComplex->number_of_ksimplices(0));
    message(InfoLevel::Default, "  1-simplices: %u", tComplex->number_of_ksimplices(1));
    message(InfoLevel::Default, "  2-simplices: %u", tComplex->number_of_ksimplices(2));
    message(InfoLevel::Default, "  3-simplices: %u", tComplex->number_of_ksimplices(3));

    // Cohomology
    Timer t3;
    mesh::Cohomology tCohomology(tComplex, aMesh);   // computes the groups
    uint64_t time_cohomology = t3.stop();

    message(InfoLevel::Default, "Cohomology: %lu ms", time_cohomology);
    message(InfoLevel::Default, "  H^0: %lu", tCohomology.get_Generators()( 0 ).size());
    message(InfoLevel::Default, "  H^1: %lu", tCohomology.get_Generators()( 1 ).size());
    message(InfoLevel::Default, "  H^2: %lu", tCohomology.get_Generators()( 2 ).size());

    delete tComplex;
}
```

---

## Summary

The `homology` module provides a complete computational topology toolkit for finite element electromagnetics:

| Component | Purpose | Entry Point |
|-----------|---------|-------------|
| **CutFactory** | Orchestrates cut generation | Main user class |
| **Cohomology** | Computes H^k groups | Direct access for advanced use |
| **Homology** | Computes H_k groups | Dual to cohomology |
| **SimplicialComplex** | Reduction engine | Internal use |
| **Topology** | Mesh analysis | Prerequisite for cuts |
| **Protoshell** | Thin shell config | Specialized geometries |

### Recommended Workflow

1. **Load mesh** and ensure proper block/sideset labeling
2. **Analyze topology** with `Topology::run()`
3. **Create CutFactory** with `PellikkaGeneralized` algorithm
4. **Set terminals** if using current constraints
5. **Run pipeline** via `CutFactory::run()`
6. **Extract results** (cuts, abstract nodes)
7. **Integrate with FEM** (DOF assignment, jump conditions)

### Performance Best Practices

- Use `PellikkaGeneralized` (the production default; no timing comparison is recorded)
- Debug builds write one `cut_<index>.vtk` per cut; release builds do not
- Pre-filter mesh to remove unnecessary entities before cohomology computation

### Further Reading

- **Theory:** `cohomology_theory_and_implementation.md`
- **Algorithms:** `cohomology_algorithms.md`
- **FEM Integration:** `src/fem/maxwell/doc/` (when available)
- **Literature:** See references in algorithm documentation

---

**Document version:** 1.0
**Last updated:** 2026-10-22
**Author:** Claude Code (claude.ai/code)
