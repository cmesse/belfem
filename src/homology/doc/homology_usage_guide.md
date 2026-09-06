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
5. **MPI Awareness**: Cohomology computation is currently serial (mesh must be on a single process)

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

**Consequences:**
- **NOT REVERSIBLE**: Cannot undo node duplication or element relinking
- **Invalidates indices**: Any code relying on original node indices will break
- **Invalidates fields**: Custom mesh fields indexed by original nodes become inconsistent
- **Changes mesh size**: Node count, element connectivity, and sideset count all increase

**Execution order requirement:**
```cpp
// CORRECT order
Mesh* tMesh = new Mesh("problem.exo");
mesh::CutFactory tFactory(tMesh, ...);
tFactory.run();                          // Mutates mesh
fem::DofManager* tDofManager = new fem::DofManager(tMesh);  // AFTER cuts

// WRONG order
fem::DofManager* tDofManager = new fem::DofManager(tMesh);  // BEFORE cuts
tFactory.run();  // ❌ DOF manager now has invalid node references
```

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

**Correct MPI workflow:**
```cpp
// 1. Load mesh on master rank ONLY
Mesh* tMesh = nullptr;
if (comm_rank() == 0) {
    tMesh = new Mesh("problem.h5");

    // 2. Run homology on master (serial)
    mesh::CutFactory tFactory(tMesh, ...);
    tFactory.run();
}

// 3. Partition AFTER homology
if (comm_rank() == 0) {
    tMesh->partition(comm_size());
}

// 4. Distribute to all ranks
if (comm_rank() == 0) {
    tMesh->distribute();
} else {
    tMesh = new Mesh();
    tMesh->receive(0);
}
```

**Future work:** Parallel cohomology is planned but not yet implemented.

---

### Memory Ownership

**Critical ownership rules:**

| Object | Created By | Owned By | Caller Action |
|--------|-----------|----------|---------------|
| `SideSet*` (cuts) | CutFactory | **Mesh** | **Do NOT delete** |
| `Node*` (abstract) | CutFactory | **Mesh** | **Do NOT delete** |
| `Node*` (duplicates) | CutFactory | **Mesh** | **Do NOT delete** |
| `Cochain*` (generators) | Cohomology | **Cohomology** | **Do NOT delete** |
| `Chain*` (generators) | Homology | **Homology** | **Do NOT delete** |
| `Protoshell*` | **User** | **User** | **Caller must delete** |

**Key rule:** CutFactory transfers all created mesh entities (nodes, sidesets) to the mesh. The mesh takes ownership and will delete them in its destructor.

```cpp
// CORRECT
CutFactory tFactory(tMesh, ...);
tFactory.run();
Cell<SideSet*>& cuts = tFactory.cuts();  // Reference only
// cuts are owned by tMesh - do NOT delete

delete tMesh;  // Mesh destructor deletes cuts and nodes

// WRONG
Cell<SideSet*>& cuts = tFactory.cuts();
delete cuts(0);  // ❌ DOUBLE DELETE when mesh destructor runs
```

---

### Thread Safety

The homology module is **NOT thread-safe**:

- `SimplicialComplex::reduce_*()` modifies internal state (chain/cochain maps)
- `CutFactory::run()` modifies mesh topology (node duplication, element relinking)
- `smithForm()` and friends modify their input matrix in place (`aMat` is taken by reference)

**OpenMP:** Do not call homology functions from parallel regions

**MPI:** Must run on single rank (see MPI Serial Limitation above)

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
- Terminals (if used) must **exist as nodes** in the mesh

❌ **Unsupported configurations:**
- Non-orientable manifolds (Möbius strip, Klein bottle)
- Meshes with dangling tetrahedra (elements connected only at nodes/edges)
- Multiply-connected phi-domains where cuts span >10^7 faces (memory limit)

**Validation pattern:**
```cpp
// Check basic mesh properties
BELFEM_ERROR(tMesh->number_of_elements() > 0, "Empty mesh");
BELFEM_ERROR(tMesh->number_of_blocks() > 0, "No blocks defined");

// Check topology analysis succeeded
mesh::Topology tTopology(tMesh);
tTopology.run();
BELFEM_ERROR(tTopology.phi_block_ids().size() > 0,
             "No phi-domain blocks identified");
```

---

### Scratch-Field Usage

**Warning:** the module uses the mesh entities' flag bits (`flag_*()` / `unflag_*()`) as
scratch and, in the SPFA infeasibility path, writes `Element::level`
(`cl_Cohomology.cpp:510-528`). It does not write `owner` or `index`.

**Consequences:**
- Do not rely on entity flags or `Element::level` after `CutFactory::run()`
- Run homology **before** partitioning because cohomology needs the undistributed mesh
  (see "MPI Serial Limitation"), not because of `owner`

```cpp
// WRONG order
tMesh->partition();  // Distributes the mesh
tFactory.run();      // ❌ Cohomology needs the undistributed mesh

// CORRECT order
tFactory.run();      // Serial mesh
tMesh->partition();  // Now distribute
```

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

// Load mesh
Mesh* tMesh = new Mesh("toroid.exo");

// Analyze topology
mesh::Topology tTopology(tMesh);
tTopology.run();

// Generate cuts
mesh::CutFactory tFactory(
    tMesh, &tTopology, {},  // Empty protoshells
    mesh::CutAlgorithm::PellikkaGeneralized
);
tFactory.run();

// Use results
for (mesh::SideSet* cut : tFactory.cuts()) {
    cut->set_domain_type(mesh::DomainType::Cut);
}
```

**That's it!** The `tFactory.cuts()` are now ready for FEM assembly with jump conditions.

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
| **Smith Normal Form** | Diagonal matrix decomposition QAR = D revealing Betti numbers and torsion |
| **Betti Number** | Rank of k-th homology/cohomology group; counts k-dimensional "holes" |
| **Homology Group** | H_k = Ker(∂_k) / Im(∂_{k+1}) (k-cycles modulo k-boundaries) |
| **Cohomology Group** | H^k = Ker(δ^k) / Im(δ^{k-1}) (k-cocycles modulo k-coboundaries) |
| **Terminal** | Electrical connection point (set of nodes) for current injection/extraction |
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
| **Tarjan's Algorithm** | DFS-based articulation point detection for manifold cleanup. Exists only as an archived prototype, `archive/graph/fn_Graph_tarjan.{hpp,cpp}`, which is outside the build |

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
// ❌ WRONG
mesh::Topology tTopology(aMesh);
mesh::CutFactory tFactory(aMesh, &tTopology, {});
tFactory.run();  // ❌ Topology::run() was never called!

// ✅ CORRECT
mesh::Topology tTopology(aMesh);
tTopology.run();  // Analyze FIRST
mesh::CutFactory tFactory(aMesh, &tTopology, {});
tFactory.run();
```

**Consequence:** `phi_block_ids()` is empty → no cohomology generators → no cuts.

---

### 2. Running Homology After Mesh Partitioning

**Problem:** Partitioning before homology breaks simplicial complex.

```cpp
// ❌ WRONG
tMesh->partition(comm_size());  // Distributes across ranks
tFactory.run();  // ❌ Mesh is partitioned - cohomology sees incomplete topology

// ✅ CORRECT
tFactory.run();  // Run on full mesh FIRST
tMesh->partition(comm_size());  // Partition AFTER
```

**Consequence:** Incomplete generators, missing cuts, or crashes due to missing connectivity.

---

### 3. Deleting Cuts or Abstract Nodes

**Problem:** Manually deleting mesh entities owned by the mesh.

```cpp
// ❌ WRONG
Cell<SideSet*>& cuts = tFactory.cuts();
delete cuts(0);  // ❌ DOUBLE DELETE - mesh owns this!

// ✅ CORRECT
Cell<SideSet*>& cuts = tFactory.cuts();
// Just use them - mesh will delete when tMesh destructor runs
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

// ✅ SAFER (pseudo-code: SideSet has no edges() accessor; walk the facets' edges yourself)
tFactory.run();
for (SideSet* cut : tFactory.cuts()) {
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
mesh::CutFactory tFactory(..., mesh::CutAlgorithm::CCR);
tFactory.run();

// ✅ Production default
mesh::CutFactory tFactory(..., mesh::CutAlgorithm::PellikkaGeneralized);
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

**Problem:** Orphaned nodes increase memory usage and may leak.

```cpp
// ⚠️ CHECK
tFactory.run();
if (tFactory.orphaned_nodes().size() > 0) {
    message(InfoLevel::Default, "Warning: %lu orphaned nodes",
            tFactory.orphaned_nodes().size());
    // Optionally: delete orphaned nodes if not needed
}
```

**Consequence:** Memory leak or unexpected behavior if orphaned nodes are referenced elsewhere.

---

### 10. Incorrect Terminal-to-Shell Mapping

**Problem:** Terminal indices don't match protoshell indices.

```cpp
// ❌ WRONG
Cell<Cell<id_t>> tTerminals = {{101}, {201}};
Cell<id_t> tShellIndices = {1, 0};  // Terminal 0 on shell 1? Mismatch!

// ✅ CORRECT
Cell<Cell<id_t>> tTerminals = {{101}, {201}};
Cell<id_t> tShellIndices = {0, 0};  // Both terminals on shell 0
```

**Consequence:** Terminals assigned to wrong thin shell → incorrect current distribution.

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
- `aUseEnrichment`: Use XFEM enrichment (experimental)

Neither argument has a default; both are passed explicitly.

**⚠️ Ownership:** CutFactory does **not** take ownership of aMesh or aTopology. Caller must ensure they outlive the factory.

### Key Methods

#### Workflow Methods

```cpp
void set_terminals(
    const Cell<Cell<id_t>>& aTerminals,
    const Cell<id_t>&       aThinShellIndices
);
```

Sets up electrical terminal definitions for current injection/extraction points.

- `aTerminals`: Cell of node ID lists, one per terminal
- `aThinShellIndices`: Maps terminals to thin shell indices

```cpp
void run();
```

**Main execution pipeline.** Performs:
1. Thin shell preprocessing
2. Cohomology computation
3. Thick-to-thin cut conversion
4. Node duplication
5. Manifold cleanup

**⚠️ Call exactly once.** Calling `run()` multiple times is undefined behavior.

#### Result Retrieval

```cpp
Cell<Node*>& abstract_nodes();
```

Returns abstract nodes (new degrees of freedom) created for cut discontinuities.

**Usage:** Assign these to FEM DOF manager for φ⁺ and φ⁻ values.

```cpp
Cell<SideSet*>& cuts();
```

Returns generated cut surfaces as SideSet objects.

**Usage:** Mark as `DomainType::Cut` and apply discontinuity constraints in FEM assembly.

```cpp
Cell<Node*>& orphaned_nodes();
```

Returns nodes not belonging to any element after node duplication.

**Usage:** Typically freed after FEM setup, but may be used for debugging.

```cpp
Vector<id_t>& thin_shell_boundaries();
```

Returns IDs of boundaries created for thin shell sidesets.

```cpp
Cell<Facet*>& thin_shell_facets(const id_t aID);
```

Returns facets for a specific thin shell ID.

#### Debugging

```cpp
void save_debug_meshes();
```

Saves one debug mesh per cut, named `cut_<index>.vtk` (`cl_CutProcessor.cpp:132-140`).
It does **not** write separate thick-cut, thin-cut or manifold files — those filenames
appear in older documentation but no code emits them.

### Usage Pattern

```cpp
#include "cl_CutFactory.hpp"
#include "cl_Topology.hpp"

// 1. Analyze mesh topology
mesh::Topology tTopology(aMesh);
tTopology.run();

// 2. Create CutFactory
mesh::CutFactory tFactory(
    aMesh,
    &tTopology,
    mProtoshells,
    mesh::CutAlgorithm::PellikkaGeneralized,
    false   // No enrichment
);

// 3. Set terminals (if needed)
Cell<Cell<id_t>> tTerminals = {{101, 102, 103}, {201, 202, 203}};
Cell<id_t> tShellIndices = {0, 0};
tFactory.set_terminals(tTerminals, tShellIndices);

// 4. Run pipeline
tFactory.run();

// 5. Extract results
for (mesh::SideSet* cut : tFactory.cuts()) {
    cut->set_domain_type(mesh::DomainType::Cut);
}

for (mesh::Node* node : tFactory.abstract_nodes()) {
    aDofManager->add_abstract_dof(node);
}
```

---

## Cohomology Computation

### Purpose

`Cohomology` computes cohomology groups H^k = Ker(δ^k) / Im(δ^{k+1}) from a simplicial complex.

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

Computes quotient space H^k = Ker(δ^k) / Im(δ^{k+1}) via Smith Normal Form.

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
for (auto& [edge_id, coeff] : gen->getSimplicesMap()) {
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

### Key Methods

```cpp
void homologyGroupOfChainComplex();
```

Main computation method. Performs:
1. Simplicial complex reduction
2. Quotient group computation (Ker/Im decomposition)
3. Generator extraction

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

// Create with terminal hints
Cell<Cell<id_t>> tTerminals = {{101, 102}, {201, 202}};
Cell<id_t> tShellIndices = {0, 0};

mesh::Homology tHomology(aMesh, tTerminals, tShellIndices);

// Compute homology
tHomology.homologyGroupOfChainComplex();

// Get H_1 generators (1-cycles)
Cell<Chain*>& generators = tHomology.get_Generators()( 1 );

// Compute orientations
Cell<Vector<real>> tDirections = {{1, 0, 0}, {0, 1, 0}};
Cell<int> orientations = tHomology.generators_orientation(tDirections);

// Reorient if needed
tHomology.reorient_generators();
```

---

## Cut Algorithm Selection

### Available Algorithms

Four algorithms are available via the `CutAlgorithm` enum:

**File:** `en_CutAlgorithm.hpp`

```cpp
enum class CutAlgorithm : uint
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
- **Advantages:**
  - Near-linear complexity for large meshes
  - Minimal Smith Normal Form computation
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
  - Slowest for large meshes
  - May produce non-manifold cuts requiring extensive cleanup

**BeltedTree**:
- **Location:** `cl_BeltedTree.{hpp,cpp}`
- **Method:** Spanning tree construction with belt identification
- **Use case:** Alternative approach when reduction methods struggle

### Selection Guidance

```cpp
// Production code - use generalized Pellikka
tFactory = new CutFactory(
    aMesh, aTopology, aProtoshells,
    CutAlgorithm::PellikkaGeneralized
);

// Legacy compatibility - original Pellikka
tFactory = new CutFactory(
    aMesh, aTopology, aProtoshells,
    CutAlgorithm::Pellikka
);

// Debugging/theory - CCR
tFactory = new CutFactory(
    aMesh, aTopology, aProtoshells,
    CutAlgorithm::CCR
);
```

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

Performs mesh topology analysis:
1. Categorizes blocks by domain type (phi/non-phi)
2. Identifies interfaces between domains
3. Classifies sidesets by type
4. Builds connectivity maps

```cpp
DomainType block_domain_type(const id_t aBlockID) const;
```

Returns domain type for a given block ID.

**Domain types:**
- `DomainType::Phi`: Region where scalar potential φ is solved
- `DomainType::Conductor`: Conducting region (source of magnetic field)
- `DomainType::Air`: Air/vacuum region
- `DomainType::Boundary`: Domain boundary
- `DomainType::Interface`: Interface between domains
- `DomainType::Cut`: Cohomology cut surface

```cpp
const Cell<id_t>& phi_block_ids() const;
const Cell<id_t>& nonphi_block_ids() const;
```

Returns lists of block IDs categorized by domain type.

```cpp
const Cell<id_t>& interface_ids() const;
```

Returns IDs of interface sidesets.

### Usage Pattern

```cpp
#include "cl_Topology.hpp"

mesh::Topology tTopology(aMesh);
tTopology.run();

// Check domain types
for (id_t blockID : tTopology.phi_block_ids()) {
    message(InfoLevel::Default, "Block %lu: Phi domain", blockID);
}

// Identify interfaces
for (id_t ssID : tTopology.interface_ids()) {
    SideSet* interface = aMesh->sideset(ssID);
    // Apply interface conditions
}

// Query specific block
if (tTopology.block_domain_type(5) == DomainType::Conductor) {
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
    Cell<id_t>   mSideSetIDs;       // Sidesets defining shell surfaces
    Cell<id_t>   mTerminalIDs;      // Terminal connection points
    Vector<real> mThicknesses;      // Shell thicknesses per layer
    Cell<string> mMaterials;        // Material names
    Cell<id_t>   mCurveIDs;         // Geometric curves

    // ... methods ...
};
```

### Key Methods

```cpp
void add_sideset(const id_t aID);
void add_terminal(const id_t aID);
void add_thickness(const real aValue);
void add_material(const string& aName);
void add_curve(const id_t aID);
```

Adds components to protoshell configuration.

```cpp
const Cell<id_t>& sideset_ids() const;
const Cell<id_t>& terminal_ids() const;
const Vector<real>& thicknesses() const;
```

Accessors for configuration data.

### Usage Pattern

```cpp
#include "cl_Protoshell.hpp"

// Create protoshell configuration
mesh::Protoshell* tShell = new mesh::Protoshell();

// Add sidesets defining shell surfaces
tShell->add_sideset(101);  // Top surface
tShell->add_sideset(102);  // Bottom surface

// Add terminal connection points
tShell->add_terminal(201);
tShell->add_terminal(202);

// Add material and thickness
tShell->add_material("YBCO");
tShell->add_thickness(1e-6);  // 1 micron

// Pass to CutFactory
Cell<Protoshell*> tProtoshells = {tShell};
CutFactory tFactory(aMesh, &tTopology, tProtoshells, ...);
```

**⚠️ Memory Management:** CutFactory does **not** take ownership of Protoshell objects. Caller must delete them after use.

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
    const Vector<id_t>& aShellBlocks
);
```

**Parameters:**
- `aMesh`: Mesh to process
- `aCohomology`: Cohomology computation results
- `aPhiBlocks`: Block IDs for phi domains
- `aNonPhiBlocks`: Block IDs for non-phi domains
- `aPhiInterfaces`: Interface sideset IDs
- `aShellBlocks`: Thin shell block IDs

### Key Methods

```cpp
Cell<Node*>& abstract_nodes();
```

Returns abstract nodes created for cut DOFs.

```cpp
id_t max_node_id() const;
id_t max_sideset_id() const;
```

Returns maximum IDs after node duplication and cut creation.

**Use case:** Ensures unique IDs when creating additional mesh entities.

```cpp
void save_debug_meshes();
```

Exports one debug mesh per cut for visualization, as `cut_<index>.vtk`
(`cl_CutProcessor.cpp:139-147`).

### Internal Workflow

The cut processing pipeline (called internally by CutFactory):

1. **Thick cut extraction**: Extract edges from cohomology generators
2. **Face identification**: Find faces incident to thick cut edges
3. **Manifold filtering**: Remove non-manifold artifacts (double pockets, etc.) — design of record; the step that actually runs is `Cohomology::remove_cut_pockets()`, see "Manifold Cleanup" below
4. **Node duplication**: Duplicate nodes on positive side of cut
5. **Element relinking**: Update element connectivity to duplicated nodes
6. **Abstract node creation**: Create DOFs for [φ] = φ⁺ - φ⁻
7. **SideSet creation**: Package cut as SideSet for FEM assembly

### Manifold Cleanup

> **Design note — design of record, not wired in.** `manifold_filter_3d()` and `check_surface()`
> do **not** exist anywhere in the tree, and no shipped code applies the three-phase scheme or the
> metric thresholds below. Tarjan articulation detection exists only as an archived prototype,
> `archive/graph/fn_Graph_tarjan.{hpp,cpp}`, which is outside the build — nothing in
> `CMakeLists.txt` references `archive/`.
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

**Purpose:** Represents formal sums of k-simplices with integer coefficients.

**File:** `cl_Chain.{hpp,cpp}`

**Structure:**
```cpp
class Chain
{
    OrderedMap<id_t, int> mData;  // Simplex ID → coefficient
    uint mDimension;               // 0 (nodes), 1 (edges), 2 (faces), 3 (elements)
};
```

**Key operations:**
```cpp
Chain operator+(const Chain& aOther) const;  // Addition
Chain operator-(const Chain& aOther) const;  // Subtraction
Chain operator*(const int aValue) const;     // Scalar multiplication

int inner_product(const Chain& aOther) const;  // <c1, c2>

const OrderedMap<id_t, int>& data() const;  // Access coefficients
```

**Usage:**
```cpp
Chain c1;
c1.add(101, 1);   // Add edge 101 with coefficient +1
c1.add(102, -1);  // Add edge 102 with coefficient -1

Chain c2;
c2.add(101, 2);

Chain c3 = c1 + c2;  // {101: 3, 102: -1}
```

---

### Cochain

**Purpose:** Dual to Chain; represents cochains in the cochain complex.

**File:** `cl_Cochain.{hpp,cpp}`

**Structure:**
```cpp
class Cochain
{
    OrderedMap<index_t, int> mData;  // Simplex index → coefficient
    uint mDimension;
};
```

**Key operations:**
```cpp
Cochain operator+(const Cochain& aOther) const;
Cochain operator-(const Cochain& aOther) const;
Cochain operator*(const int aValue) const;

int evaluate(const Chain& aChain) const;  // <φ, c>

const OrderedMap<index_t, int>& data() const;
```

**Usage:**
```cpp
Cochain phi;
phi.add(0, 1);   // Add edge at index 0 with coefficient +1
phi.add(1, -1);  // Add edge at index 1 with coefficient -1

Chain c = ...;
int value = phi.evaluate(c);  // Cochain evaluation on chain
```

---

### CutData

**Purpose:** Stores topology and metadata for a single cut.

**File:** `cl_CutData.{hpp,cpp}`

**Structure:**
```cpp
class CutData
{
    Cell<id_t>           mEdgeIDs;        // Edges in thick cut
    Cell<id_t>           mElementIDs;     // Elements adjacent to cut
    Cell<Facet*>         mFacets;         // Facets forming thin cut
    DynamicBitset        mEdgeBitset;     // Fast edge membership test
    DynamicBitset        mNodeBitset;     // Fast node membership test
    DomainType           mType;           // Cut classification
};
```

**Key methods:**
```cpp
void determine_cut_case(Element* aElement);
```

Determines how cut intersects an element (which faces are on cut).

**Use case:** Guides node duplication and element relinking.

```cpp
bool node_is_on_cut(const id_t aNodeID) const;
bool edge_is_on_cut(const id_t aEdgeID) const;
```

Fast membership tests using bitsets.

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

**Note:** The earlier `SideSetFactory` class and the
`CutFactory::create_sidesets_2d/3d` methods that built sidesets directly from
generator support were retired in 2026. `CutProcessor`/`CutData` is now the
only sideset path.

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
4. Relinks elements to duplicated nodes
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

Tests if node is on the cut surface.

**Returns:** `true` if node index is set in cut bitset.

```cpp
void create_duplicates(id_t& aMaxNodeID, Cell<Node*>& aAbstractNodes);
```

Creates duplicate nodes for positive side of cut.

**Algorithm:**
1. Iterate through nodes on cut
2. Create new Node with unique ID
3. Store in `mNodeDuplicates` map
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

// Test membership
if (tCutSet.test(someNode)) {
    // Node is on cut
}

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
void create_TreeField(Mesh* tEdgeMesh, string tFieldName);
void create_cohomologyField(Mesh* tMeshEdge);
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
mesh::CutFactory tFactory(
    aMesh, &tTopology, tProtoshells,
    mesh::CutAlgorithm::BeltedTree  // Use belted tree algorithm
);
```

---

### fn_Smith - Smith Normal Form Functions

**Purpose:** Integer matrix decomposition for computing homology/cohomology groups.

**File:** `fn_Smith.{hpp,cpp}`

**Mathematical background:**

The Smith Normal Form decomposes an integer matrix A into:
```
Q * A * R = D
```
Where:
- Q, R are unimodular matrices (det = ±1)
- D is diagonal matrix with d_1 | d_2 | ... | d_r (divisibility chain)
- Diagonal entries d_i are **invariant factors**
- Number of nonzero entries r is the **rank**

**Homology application:**

Given boundary operator ∂_k : C_k → C_{k-1} with matrix representation A:
- Ker(∂_k) dimension = n - r (n = number of k-chains)
- Im(∂_{k+1}) dimension = r
- H_k rank (Betti number) = dim(Ker) - dim(Im)
- Torsion subgroup = {coefficients d_i where d_i > 1}

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

**Returns:**
- Q: Left transformation matrix
- Q_: Inverse of Q
- R: Right transformation matrix
- R_: Inverse of R
- rank: Number of nonzero diagonal entries
- Number of rows/columns

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

Checks whether the pivot `B[k,k]` divides every entry of the trailing submatrix `B[k+1:end, k+1:end]` (Kaczynski et al., Smith-form step); returns the first offending `(row, col, quotient)`.

**Returns:** (divisible, row_index, col_index, violating_entry)

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

// Smith decomposition of ∂_k
auto [Q, Q_, R, R_, rank, size] = smithForm(boundary[k]);

// Extract homology generators
// Ker(∂_k) = columns of R corresponding to zero diagonal entries
// Betti_k = number of zero diagonal entries - rank(∂_{k+1})
```

**Performance notes:**

- **Complexity:** O(n³) for dense matrices (n = dimension)
- **Optimization:** Uses smallest-entry pivoting for faster convergence
- **Integer overflow:** No protection for very large matrices (TODO: arbitrary-precision)
- **Sparsity:** Does not exploit sparsity (future improvement)

**Elementary operations:**

The following operations preserve equivalence and update transformation matrices:
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

// Boundary matrix for 2-dimensional complex
Matrix<int> boundary2(4, 6);
// ... populate with boundary coefficients ...

// Compute Smith Normal Form
auto [Q, Q_, R, R_, rank, dim] = smithForm(boundary2);

// Betti number: number of zero diagonal entries
uint betti = dim - rank;

message(InfoLevel::Default, "H_2 Betti number: %u", betti);

// Extract generators
Matrix<int> generators(R.n_rows(), betti);
uint col = 0;
for (uint i = rank; i < dim; ++i) {
    for (uint row = 0; row < R.n_rows(); ++row) {
        generators(row, col) = R(row, i);
    }
    col++;
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

    // 2. Create empty protoshell list (no thin shells)
    Cell<mesh::Protoshell*> tProtoshells;

    // 3. Create factory
    mesh::CutFactory tFactory(
        aMesh,
        &tTopology,
        tProtoshells,
        mesh::CutAlgorithm::PellikkaGeneralized
    );

    // 4. Run pipeline
    tFactory.run();

    // 5. Process results
    for (mesh::SideSet* cut : tFactory.cuts()) {
        cut->set_domain_type(mesh::DomainType::Cut);
        message(InfoLevel::Default, "Cut %lu: %u facets",
                cut->id(), cut->number_of_facets());
    }
}
```

---

### Pattern 2: Terminal-Driven Cuts

```cpp
void generate_terminal_cuts(Mesh* aMesh)
{
    // Topology analysis
    mesh::Topology tTopology(aMesh);
    tTopology.run();

    // Define terminals (node IDs)
    Cell<Cell<id_t>> tTerminals;
    tTerminals(0) = {101, 102, 103};  // Terminal 1
    tTerminals(1) = {201, 202, 203};  // Terminal 2

    // Map to thin shell indices (0 = not on thin shell)
    Cell<id_t> tShellIndices = {0, 0};

    // Create factory
    mesh::CutFactory tFactory(
        aMesh, &tTopology, {},
        mesh::CutAlgorithm::PellikkaGeneralized
    );

    // Set terminals
    tFactory.set_terminals(tTerminals, tShellIndices);

    // Run
    tFactory.run();

    // Results
    message(InfoLevel::Default, "Generated %lu cuts",
            tFactory.cuts().size());
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
                i, gen->getSimplicesMap().size());
    }
}
```

---

### Pattern 4: Algorithm Performance Comparison

```cpp
void benchmark_algorithms(Mesh* aMesh)
{
    mesh::Topology tTopology(aMesh);
    tTopology.run();

    Cell<mesh::Protoshell*> tProtoshells;

    CutAlgorithm algorithms[] = {
        CutAlgorithm::PellikkaGeneralized,
        CutAlgorithm::Pellikka,
        CutAlgorithm::CCR,
        CutAlgorithm::BeltedTree
    };

    for (auto alg : algorithms) {
        Timer timer;

        mesh::CutFactory tFactory(aMesh, &tTopology, tProtoshells, alg);
        tFactory.run();

        uint64_t elapsed = timer.stop();
        message(InfoLevel::Default, "Algorithm %d: %lu ms, %lu cuts",
                static_cast<uint>(alg), elapsed, tFactory.cuts().size());
    }
}
```

---

### Pattern 5: Thin Shell with Terminals

```cpp
void process_thin_shell_coil(Mesh* aMesh)
{
    // Topology
    mesh::Topology tTopology(aMesh);
    tTopology.run();

    // Protoshell configuration
    mesh::Protoshell* tShell = new mesh::Protoshell();
    tShell->add_sideset(101);  // Top surface
    tShell->add_sideset(102);  // Bottom surface
    tShell->add_terminal(201);  // Input terminal
    tShell->add_terminal(202);  // Output terminal
    tShell->add_material("YBCO");
    tShell->add_thickness(1e-6);

    Cell<mesh::Protoshell*> tProtoshells = {tShell};

    // Terminals
    Cell<Cell<id_t>> tTerminals = {{301, 302}, {401, 402}};
    Cell<id_t> tShellIndices = {0, 0};  // Both on shell 0

    // Factory
    mesh::CutFactory tFactory(
        aMesh, &tTopology, tProtoshells,
        mesh::CutAlgorithm::PellikkaGeneralized
    );

    tFactory.set_terminals(tTerminals, tShellIndices);
    tFactory.run();

    // Cleanup
    delete tShell;

    // Use results
    for (mesh::SideSet* cut : tFactory.cuts()) {
        // Apply jump conditions in FEM
    }
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
   CutFactory tFactory(..., CutAlgorithm::Pellikka);
   ```

2. **Enable debug output:**
   ```cpp
   tFactory.run();
   // save_debug_meshes() is a member of CutProcessor (public,
   // cl_CutProcessor.hpp:116), not of CutFactory. CutFactory calls it
   // automatically in debug builds (cl_CutFactory.cpp:477-479, compiled out
   // under NDEBUG); it writes one file per cut: cut_<index>.vtk
   ```
   For curve meshes, `CutFactory::save_curve_debug_meshes()` is public
   (`cl_CutFactory.hpp:184`) and writes `curve_<id>.vtk`. Inspect either family in ParaView.

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

message(InfoLevel::Default, "Phi blocks: %lu",
        tTopology.phi_block_ids().size());
message(InfoLevel::Default, "Interfaces: %lu",
        tTopology.interface_ids().size());
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
CutFactory tFactory(..., CutAlgorithm::PellikkaGeneralized);
```

If still slow:
- **Profile:** Use `Profiler` to identify hotspot
- **Mesh quality:** Simplify mesh if possible (remove small features)
- **MPI:** Current implementation is serial (future parallelization planned)

---

### Problem: Abstract Nodes Not Created

**Symptom:** `tFactory.abstract_nodes().size() == 0` after `run()`.

**Cause:** No cuts were generated (topology is simply-connected).

**Diagnosis:**
```cpp
message(InfoLevel::Default, "Num cuts: %lu", tFactory.cuts().size());
```

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

2. **Check abstract node count:**
   ```cpp
   // Should match number of cuts
   BELFEM_ERROR(tFactory.abstract_nodes().size() == tFactory.cuts().size(),
                "Abstract node count mismatch");
   ```

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

int main()
{
    // Load toroidal mesh
    Mesh* tMesh = new Mesh("toroid.exo");

    // Analyze topology
    mesh::Topology tTopology(tMesh);
    tTopology.run();

    // Create factory (no thin shells)
    mesh::CutFactory tFactory(
        tMesh,
        &tTopology,
        {},  // Empty protoshells
        mesh::CutAlgorithm::PellikkaGeneralized
    );

    // Generate cuts
    Timer timer;
    tFactory.run();
    uint64_t elapsed = timer.stop();

    message(InfoLevel::Default,
            "Generated %lu cuts in %lu ms",
            tFactory.cuts().size(), elapsed);

    // Apply to FEM
    for (mesh::SideSet* cut : tFactory.cuts()) {
        cut->set_domain_type(mesh::DomainType::Cut);
    }

    // Create DOF manager
    fem::DofManager* tDofManager = new fem::DofManager(tMesh);
    for (mesh::Node* node : tFactory.abstract_nodes()) {
        tDofManager->add_abstract_dof(node, "phi");
    }

    // ... FEM assembly and solve ...

    delete tDofManager;
    delete tMesh;
    return 0;
}
```

---

### Example 2: Superconducting Cable

```cpp
#include "cl_CutFactory.hpp"
#include "cl_Protoshell.hpp"
#include "cl_Topology.hpp"

int main()
{
    Mesh* tMesh = new Mesh("cable.h5");

    // Topology
    mesh::Topology tTopology(tMesh);
    tTopology.run();

    // Protoshell for thin shell approximation
    mesh::Protoshell* tShell = new mesh::Protoshell();
    tShell->add_sideset(101);  // Cable surface
    tShell->add_terminal(201);  // Input
    tShell->add_terminal(202);  // Output
    tShell->add_material("BSCCO");
    tShell->add_thickness(0.2e-3);  // 0.2 mm

    Cell<mesh::Protoshell*> tProtoshells = {tShell};

    // Terminals
    Cell<Cell<id_t>> tTerminals;
    tTerminals(0) = {1001, 1002, 1003};  // Input terminal nodes
    tTerminals(1) = {2001, 2002, 2003};  // Output terminal nodes
    Cell<id_t> tShellIndices = {0, 0};

    // Factory
    mesh::CutFactory tFactory(
        tMesh, &tTopology, tProtoshells,
        mesh::CutAlgorithm::PellikkaGeneralized
    );
    tFactory.set_terminals(tTerminals, tShellIndices);
    tFactory.run();

    // Results
    message(InfoLevel::Default, "Cuts: %lu", tFactory.cuts().size());
    message(InfoLevel::Default, "Abstract DOFs: %lu",
            tFactory.abstract_nodes().size());

    // Cleanup
    delete tShell;
    delete tMesh;
    return 0;
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
    delete tFlagged;
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
- Enable debug output only when troubleshooting (overhead ~10%)
- Pre-filter mesh to remove unnecessary entities before cohomology computation
- For very large meshes (>10^7 elements), consider domain decomposition (future feature)

### Further Reading

- **Theory:** `cohomology_theory_and_implementation.md`
- **Algorithms:** `cohomology_algorithms.md`
- **FEM Integration:** `src/fem/maxwell/doc/` (when available)
- **Literature:** See references in algorithm documentation

---

**Document version:** 1.0
**Last updated:** 2026-08-31
**Author:** Claude Code (claude.ai/code)
