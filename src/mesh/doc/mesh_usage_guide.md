# BELFEM Mesh Module - Usage Guide {#mesh_mesh_usage_guide}

**Module:** `src/mesh`
**Purpose:** Comprehensive guide to BELFEM's mesh data structures, I/O, and topology management

**Date:** 2026-01-16
**Revision:** 1.0

---

## Revision History

| Date | Version | Changes |
|------|---------|---------|
| 2026-01-16 | 1.0 | Initial documentation |
| 2026-01-16 | 1.1 | Applied external review feedback: Added Quick-Start section; Added comprehensive Glossary (Core Concepts, Element Types, Topology, Partitioning, I/O, Fields); Updated METIS link to GitHub |
| 2026-01-16 | 1.2 | Added Critical Contracts and Invariants section (finalize() MPI semantics, unfinalize() cache invalidation, ownership vs derived containers, connectivity invalidation rules, ID vs index stability) |
| 2026-01-16 | 1.3 | Added complete enum reference tables: ElementType (29 types), Connectivity (36 types), InterpolationOrder (8 types), InterpolationType (11 types) with detailed descriptions and use cases |
| 2026-01-16 | 1.4 | Added ThinShell class documentation with layered structure details, usage examples, and physical context for electromagnetic simulations |

---

## Quick-Start

**Minimal example to load, finalize, and query a mesh:**

```cpp
#include "cl_Mesh.hpp"

// Load mesh from file (the reader finalizes the mesh itself)
Mesh* tMesh = new Mesh("mesh.msh");

// finalize() is only needed after building or modifying a mesh by hand
tMesh->finalize();

// Query basic properties
message(InfoLevel::Default, "%lu nodes, %lu elements, %u blocks",
    (luint)tMesh->number_of_nodes(),
    (luint)tMesh->number_of_elements(),
    tMesh->number_of_blocks());

// Clean up
delete tMesh;
```

**That's it!** Most operations require finalization - without it, you'll get crashes or incorrect results.

---

## Glossary

**Core Concepts:**

| Term | Definition |
|------|------------|
| **Node** | Spatial point with (x,y,z) coordinates; inherits from Vertex (graph properties) and Basis (DOF connectivity) |
| **Element** | Volumetric (3D) or surface (2D) finite element; connects nodes, has type (TRI, QUAD, TET, HEX, etc.) |
| **Edge** | Topological 1D entity connecting two nodes; used for H(curl) DOF allocation in Nédélec elements |
| **Face** | Topological 2D entity (triangle/quad face of 3D element); used for H(div) DOF allocation |
| **Facet** | Lightweight reference to element's face/edge; stores master element, local index, orientation |
| **Block** | Group of elements of same type, typically one material or subdomain |
| **SideSet** | Group of boundary facets (surface elements on domain boundaries) |
| **ID (id_t)** | Permanent entity identifier from mesh file (may be sparse: 1, 5, 100, ...) |
| **Index (index_t)** | 0-based position in container after reordering (continuous: 0, 1, 2, ...) |
| **Owner (proc_t)** | MPI rank owning this entity (relevant only in parallel) |
| **Ghost** | Entity owned by another MPI rank but stored locally for shared boundary |
| **Duplicate Node** | Copy of original node created for cohomology cuts; allows φ⁺ ≠ φ⁻ on cut surfaces |
| **Finalization** | Critical step that computes indices, connectivities, maps, and orientations |

**Element Types:**

| Term | Definition |
|------|------------|
| **ElementType** | Enum identifying element geometry: LINE2, TRI3, QUAD4, TET4, HEX8, etc. Not every enumerator is constructible — see the element-type table below |
| **Order** | Polynomial degree of element: 1 (linear), 2 (quadratic), 3 (cubic), 4 (quartic), 5 (quintic) |
| **Corner Nodes** | Vertices of element (e.g., 4 for TET4, 8 for HEX8) |
| **Mid-Edge Nodes** | Nodes on element edges (quadratic+ elements only) |
| **Face Center Nodes** | Nodes on element faces (cubic+ elements only) |
| **Volume Center Node** | Node at element centroid (HEX27, HEX64, etc.) |
| **Curved Element** | Element with non-linear geometry (requires higher-order integration) |
| **Thin-Shell Element** | Special element variant for layered structures (cohomology-aware node ordering) |

**Topology and Connectivity:**

| Term | Definition |
|------|------------|
| **Connectivity** | Relationship between entities (NodeToElement, ElementToElement, etc.) - 36 types total |
| **Master Element** | Element that owns a facet (for boundary facets, master is the interior element) |
| **Slave Element** | Neighbor element sharing a facet (null on domain boundary) |
| **Neighbor** | Element sharing a facet with another element |
| **Manifold** | Mesh where each edge/face belongs to correct number of elements (1 for boundary, 2 for interior) |
| **Hanging Node** | Node on element face/edge but not a vertex of neighboring element (p-refinement artifact) |

**Partitioning and MPI:**

| Term | Definition |
|------|------------|
| **Partition** | Decomposition of mesh into subdomains for MPI parallelization |
| **Distributor** | Class that sends mesh partitions to MPI ranks |
| **Ghost Layer** | Shell of elements owned by neighbors but stored locally for shared boundaries |
| **METIS** | Graph partitioning library for load balancing (default) |
| **SCOTCH** | Graph partitioning library wrapped in src/math/graph; not used by the mesh Partitioner, which is METIS-only |
| **Continuous Partition** | Partition where each subdomain is a single connected component (no isolated regions) |

**I/O Formats:**

| Term | Definition |
|------|------------|
| **Gmsh (.msh)** | Open-source mesh format (industry standard for generation) |
| **HDF5 (.hdf5, .bfm)** | BELFEM native format (.bfm, HDF5-based), written and read on rank 0 |
| **Exodus II (.exo)** | Sandia format for time series (ParaView-compatible) |
| **VTK (.vtk)** | Visualization format (ParaView-compatible, simple ASCII/binary) |

**Fields and Data:**

| Term | Definition |
|------|------------|
| **Field** | Data array associated with mesh entities (nodes, elements, edges, faces) |
| **EntityType** | Enum specifying which entities field is defined on: NODE, ELEMENT, EDGE, FACE |
| **Global Variable** | Scalar mesh-level data (e.g., total current, total energy) |
| **Time Step** | Index of current time in transient simulation |
| **Time Stamp** | Physical time value for current time step |

---

## Common Pitfalls

**Read this section first to avoid frequent mistakes:**

### 1. Forgetting to Finalize the Mesh

```cpp
// File-loaded meshes are finalized by the reader itself:
Mesh* tMesh = new Mesh("mesh.msh");   // already finalized

// WRONG: Using a hand-assembled mesh without finalization
// Access indices, connectivities - WILL FAIL

// CORRECT: Finalize after assembling or modifying a mesh by hand
tMesh->finalize();  // Sets up connectivities, indices, etc.
```

**Why:** `finalize()` computes element indices, sets up node-element connectivities, orients facets, and performs essential topology setup. Without it, many mesh operations will fail or produce incorrect results.

### 2. Accessing Entities by ID vs Index

```cpp
// ID (permanent identifier from mesh file)
mesh::Node* node = tMesh->node(12345);  // Access by ID

// Index (position in container, 0-based)
mesh::Node* node = tMesh->nodes()(42);  // Direct container access

// PITFALL: Mixing ID and index
id_t nodeID = 100;
// tMesh->nodes()(nodeID);  // WRONG if nodeID != index
tMesh->node(nodeID);        // CORRECT
```

**Why:** IDs are permanent identifiers from the mesh file and may be sparse (e.g., 1, 5, 100, 1000). Indices are continuous 0-based positions after reordering. Use `node(id_t)` for ID lookup (uses internal map), `nodes()(index_t)` for direct container access.

### 3. Memory Management of Mesh Entities

```cpp
// WRONG: Manual deletion of mesh entities
mesh::Node* node = tMesh->nodes()(0);
delete node;  // DON'T DO THIS

// CORRECT: Mesh owns all entities
delete tMesh;  // Mesh destructor handles cleanup
```

**Why:** The `Mesh` class owns all entities (nodes, elements, blocks, sidesets, etc.) and deallocates them in its destructor. Manual deletion causes double-free errors.

### 4. Modifying Mesh After FEM Kernel Initialization

```cpp
// WRONG: Creating edges after DofManager setup
tKernel.create_fields();  // DofManager initialized
tMesh->create_edges();    // TOO LATE - DOFs already allocated

// CORRECT: Complete mesh topology before FEM setup
tMesh->create_edges();     // Create edges first
tMesh->create_faces();     // Create faces if needed
tKernel.create_fields();   // Then initialize FEM
```

**Why:** The FEM kernel allocates DOFs based on mesh topology. Modifying topology afterward invalidates DOF assignments.

### 5. Parallel Mesh Ownership Confusion

```cpp
// In MPI context
Mesh* tMesh = new Mesh("mesh.msh");
tMesh->partition(comm_size());  // Partition mesh

// PITFALL: Assuming partition() shrinks the containers
for (mesh::Node* node : tMesh->nodes()) {
    // still every node of the full mesh on the master; other ranks are empty
}
```

**Why:** `partition()` only assigns `owner()` on the master's full mesh; `nodes()` / `elements()` still hold everything. Only after `Distributor::run()` does a rank hold a partial mesh.

### 6. Element Type Confusion

```cpp
// PITFALL: Assuming homogeneous element types
mesh::Block* block = tMesh->block(1);
ElementType type = block->element_type();  // Type of FIRST element

// BETTER: Check if heterogeneous (rare)
ElementType firstType = block->element(0)->type();
for (mesh::Element* elem : block->elements()) {
    if (elem->type() != firstType) {
        // Handle heterogeneous block
    }
}
```

**Why:** Blocks typically contain one element type, but BELFEM allows heterogeneous blocks. Most FEM operations assume homogeneity within blocks.

### 7. Node Coordinates Modification

```cpp
// Modifying coordinates after finalization
mesh::Node* node = tMesh->node(100);
node->set_coords(1.0, 2.0, 3.0);  // OK, but...

// PITFALL: Forgetting to update element geometry
// Curved elements cache jacobians - invalidated by coord changes!

// SAFER: Scale entire mesh
tMesh->scale_mesh(0.001);  // Convert mm to m (affects all nodes)
```

**Why:** Coordinate changes invalidate cached geometric data in curved elements. Use `scale_mesh()` for uniform scaling, or manually update affected elements.

### 8. Sideset vs Block Confusion

```cpp
// Blocks contain volumetric elements (or 2D surface elements in 2D mesh)
mesh::Block* block = tMesh->block(1);
Cell<mesh::Element*>& elements = block->elements();

// Sidesets contain boundary facets (surface elements on 3D volume boundaries)
mesh::SideSet* sideset = tMesh->sideset(2);
Cell<mesh::Facet*>& facets = sideset->facets();

// PITFALL: Confusing facets with elements
// Facets are lightweight wrappers around elements
mesh::Facet* facet = facets(0);
mesh::Element* master = facet->master();  // The actual element
```

**Why:** Facets are references to element faces/edges, not standalone elements. They store master element, local facet index, and orientation.

### 9. Edge/Face Creation Order (3D Elements)

```cpp
// WRONG: Creating faces before edges (3D elements)
tMesh->create_faces();     // May fail - needs edges first
tMesh->create_edges();

// CORRECT: Edges before faces
tMesh->finalize();         // Finalize main mesh first
tMesh->create_edges();     // Create edges
tMesh->finalize_edges();   // Finalize edge connectivity
tMesh->create_faces();     // Then faces (depend on edges)
tMesh->finalize_faces();   // Finalize face connectivity
```

**Why:** Face DOFs often reference edge DOFs in Nédélec H(curl) and H(div) elements. Edge containers must exist before face creation.

### 10. Storing Indices Across Operations

```cpp
// WRONG: Storing indices across finalize/partition
Cell<index_t> nodeIndices;
for (mesh::Node* node : tMesh->nodes()) {
    if (node->x() > 0.0) {
        nodeIndices.push(node->index());
    }
}
tMesh->partition(comm_size());  // REORDERS INDICES!
// nodeIndices now points to wrong nodes

// CORRECT: Store IDs (permanent)
Cell<id_t> nodeIDs;
for (mesh::Node* node : tMesh->nodes()) {
    if (node->x() > 0.0) {
        nodeIDs.push(node->id());  // ID is stable
    }
}
tMesh->partition(comm_size());
// nodeIDs still valid - access via tMesh->node(id)
```

**Why:** Indices are reordered by `finalize()`, `partition()`, and `Distributor::run()`. IDs are permanent identifiers that remain stable across all operations.

---

## Critical Contracts and Invariants

**Understanding these rules is essential for correct mesh usage.**

### finalize() Behavior and MPI Semantics

The `finalize()` method has **different behavior on master vs non-master MPI ranks**:

```cpp
// Master proc (rank 0)
if (comm_rank() == 0) {
    Mesh* tMesh = new Mesh("mesh.msh");
    tMesh->finalize();  // Full finalization: indices, connectivity, orientation
}

// Non-master proc
else {
    Mesh* tMesh = new Mesh("mesh.msh");
    tMesh->finalize();  // Limited finalization: basic setup only
}
```

**Master proc finalization:**
- Computes all element and node indices
- Builds ID→entity maps
- Computes facet orientations
- Populates connectivity structures

**Non-master proc finalization (before distribution):**
- Minimal setup (entity counts, basic properties)
- Full finalization occurs AFTER the distributor has run and `partial_mesh()` has handed the rank its local mesh

**Rule:** call `finalize()` after mesh loading. After distribution, only the **root** needs it — `Distributor::run()` finalizes each worker's mesh itself (`cl_Mesh_Distributor.cpp:248`).

### unfinalize() Is Cache Invalidation, Not Rollback

The `unfinalize()` method is **NOT** a symmetric undo operation:

```cpp
tMesh->finalize();      // Computes indices, connectivity, maps
tMesh->unfinalize();    // Clears caches, resets flags
// Mesh is now in PARTIALLY INCONSISTENT state
```

**What unfinalize() DOES:**
- Resets the Node*/Edge*/Face* connectivity flags (ElementToNode is kept)
- Resets `mIsFinalized` flag
- Invalidates cached data

**What unfinalize() DOES NOT DO:**
- Restore original indices (indices remain modified)
- Reset the element-to-element and facet-to-facet links (kept because they are expensive to recompute; see the commented-out block in `Mesh::unfinalize()`). Node/edge/face adjacency arrays ARE freed.
- Undo facet orientation (orientations remain)

**Use case:** Modifying mesh topology (adding elements, repartitioning) requires `unfinalize()` → modify → `finalize()`.

**Warning:** Avoid relying on unfinalize() to restore exact pre-finalization state. It's a cache invalidation tool, not a rollback mechanism.

### Ownership vs Derived Containers

**Primary ownership containers (allocated explicitly):**
- `mNodes` - Mesh owns these nodes
- `mBlocks` - Mesh owns these blocks
- `mSideSets` - Mesh owns these sidesets

**Derived containers (views into primary containers):**
- `mElements` - References to elements owned by blocks
- `mFacets` - References to facets owned by sidesets

**Critical distinction:**

```cpp
// Elements are owned by blocks
mesh::Block* block = tMesh->block(1);
Cell<mesh::Element*>& blockElems = block->elements();  // Owner

// Mesh mElements is a VIEW into all block elements
Cell<mesh::Element*>& meshElems = tMesh->elements();   // Derived view

// WRONG: Deleting from derived view
for (mesh::Element* elem : meshElems) {
    delete elem;  // DON'T DO THIS - blocks own these
}

// CORRECT: Mesh destructor deletes blocks, which delete elements
delete tMesh;  // Proper cleanup chain
```

**Facets work similarly:**
- Sidesets OWN facets
- `Mesh::mFacets` is a derived view
- Facet master/slave elements are REFERENCES, not owned

**Rule:** Never manually delete entities from derived containers. Delete the mesh or the owning container (block/sideset).

### Connectivity Invalidation Rules

**Connectivities are invalidated by topology changes:**

| Operation | Invalidates | Reason |
|-----------|-------------|--------|
| `Block::insert_element()` + re-finalize | NodeToElement, ElementToElement | New element neighbors |
| `partition()` | All connectivities | Element ownership changes |
| `Distributor::run()` | All connectivities | Mesh decomposition |
| `create_edges()` | Edge-related connectivities | New edge entities created |
| `create_faces()` | Face-related connectivities | New face entities created |
| `Node::add_duplicate()` | NodeToElement | Node connectivity changes |
| `unfinalize()` | Node*/Edge*/Face* connectivities (selective; ElementToNode kept) | Cache invalidation |

**Safe operations (do NOT invalidate connectivity):**
- `node->set_coords()` - Geometric change only
- `scale_mesh()` - Coordinate scaling
- `field_data()` access - Pure data operations

**Pattern for topology modifications:**

```cpp
tMesh->finalize();              // Initial finalization
// ... use mesh ...

// Modify topology
tMesh->unfinalize();            // Invalidate caches
tMesh->block(id)->insert_element(newElem);  // Topology change
tMesh->finalize();              // Recompute everything

// Connectivities are now valid again
```

**Rule:** Never store cached connectivity references across topology-modifying operations. Always re-access after `finalize()`.

### ID vs Index Stability

**IDs are permanent, indices are NOT:**

```cpp
// Load mesh
Mesh* tMesh = new Mesh("mesh.msh");
mesh::Node* node = tMesh->node(12345);  // ID = 12345

id_t id = node->id();           // 12345 (STABLE)
index_t idx1 = node->index();   // e.g., 42 (UNSTABLE)

// Finalize changes indices
tMesh->finalize();
index_t idx2 = node->index();   // NOW DIFFERENT (e.g., 5042)

// Partition changes indices again
tMesh->partition(comm_size());
index_t idx3 = node->index();   // DIFFERENT AGAIN (e.g., 17)

// But ID remains stable
BELFEM_ASSERT(node->id() == 12345, "ID changed!");  // Always passes
```

**Operations that reorder indices:**
1. `finalize()` - assigns continuous 0-based indices in container order (does not sort)
2. `partition()` - Reorders by MPI rank assignment
3. `Distributor::run()` - Reorders within local partition
4. Any operation involving `unfinalize()` → modify → `finalize()`

**Hard rule:** **NEVER** store indices across these operations. Store IDs instead.

**Example pitfall:**

```cpp
// WRONG: Storing indices
Cell<index_t> boundaryNodeIndices;
for (mesh::Node* node : tMesh->nodes()) {
    if (node->x() < 0.01) {
        boundaryNodeIndices.push(node->index());
    }
}

tMesh->partition(comm_size());  // INDICES NOW INVALID

// CORRECT: Storing IDs
Cell<id_t> boundaryNodeIDs;
for (mesh::Node* node : tMesh->nodes()) {
    if (node->x() < 0.01) {
        boundaryNodeIDs.push(node->id());
    }
}

tMesh->partition(comm_size());  // IDs still valid

// Access later via ID lookup
for (id_t id : boundaryNodeIDs) {
    mesh::Node* node = tMesh->node(id);  // Works correctly
}
```

**Performance note:** ID lookup uses `Map<id_t, Node*>` with O(1) average time. Index access via `nodes()(index)` is direct array access but requires index stability.

**Rule:** Use IDs for persistent references, indices only for temporary iteration within a single finalized state.

---

## Architecture Overview

The mesh module provides data structures and I/O for finite element meshes with support for:
- Multiple element types (1D/2D/3D, linear to 5th order)
- Parallel mesh partitioning via METIS
- Multiple file formats (Gmsh, HDF5, Exodus, VTK)
- Higher-order topological entities (edges, faces for Nédélec elements)
- Tensor product meshes for structured grids
- Thin-shell elements for layered structures
- Node duplication for cohomology cuts

**Key design principles:**
- **Single ownership:** Mesh owns all entities (nodes, elements, etc.)
- **ID-based lookup:** Fast access via `Map<id_t, T*>` for sparse IDs
- **Index-based iteration:** Continuous 0-based indices for containers
- **Eager connectivity:** computed in `finalize()` when `Connectivity::Compute` is set; the bitset records which ones exist
- **Factory pattern:** Element creation via `Element_Factory`
- **Format-agnostic:** Common interface across I/O backends

---

## Core Classes

### Mesh Container

**File:** `cl_Mesh.{hpp,cpp}`

The `Mesh` class is the top-level container for all mesh entities:

```cpp
class Mesh {
    // Entity containers
    Cell<mesh::Node*>     mNodes;
    Cell<mesh::Element*>  mElements;
    Cell<mesh::Edge*>     mEdges;
    Cell<mesh::Face*>     mFaces;
    Cell<mesh::Facet*>    mFacets;

    // Grouping containers
    Cell<mesh::Block*>    mBlocks;
    Cell<mesh::SideSet*>  mSideSets;

    // Data containers
    Cell<mesh::Field*>    mFields;
    Cell<mesh::GlobalVariable*> mGlobalVariables;

    // ID-based maps for fast lookup
    Map<id_t, mesh::Node*>    mNodeMap;
    Map<id_t, mesh::Element*> mElementMap;
    Map<id_t, mesh::Block*>   mBlockMap;
    Map<id_t, mesh::SideSet*> mSideSetMap;
    // ... (additional maps for edges, faces, facets)
};
```

**Key properties:**
- `mMasterProc`: MPI rank that owns this mesh (default: 0)
- `mNumberOfDimensions`: Spatial dimensionality (1, 2, or 3)
- `mNumberOfPartitions`: MPI partitions (1 for serial)
- `mIsFinalized`: Flag indicating finalization status
- `mConnectivities`: Bitset tracking computed connectivities

### Node Class

**File:** `cl_Node.{hpp,cpp}`

Represents a spatial point with 3D coordinates:

```cpp
class Node : public Vertex {
    real mCoords[3];             // x, y, z coordinates
    Node** mDuplicates;          // For cohomology cuts
    int mNumberOfDuplicates;     // Negative if this is duplicate
};
```

**Inheritance:** `Node` inherits from `Vertex`, which inherits from `Basis`, which inherits from `graph::Vertex`:
- `graph::Vertex`: id, index, owner, level and the 8-slot flag bitset
- `Basis`: adds the hanging-source (weights) table and the DOF container
- `Vertex`: adds node/edge/face/facet/element connectivity containers
- `Node`: adds spatial coordinates and duplication

**Key methods:**
```cpp
real x() const;                         // X-coordinate
real y() const;                         // Y-coordinate
real z() const;                         // Z-coordinate
Vector<real> coords() const;            // {x, y, z}
void set_coords(real x, real y, real z);

// Duplication (for cohomology cuts)
uint number_of_duplicates() const;
Node* duplicate(uint aIndex);
Node* original();                       // Returns this if not duplicate
bool is_duplicate() const;
```

### Element Class

**File:** `cl_Element.{hpp,cpp}`

Base class for all element types (LINE, TRI, QUAD, TET, HEX, PENTA, PYRA):

```cpp
class Element : public Basis {
protected:
    bool mCurvedFlag;                // Geometry curvature flag
    uint16_t mElementTags[ 2 ];      // [0] geometry tag (= block id), [1] physical tag, both from gmsh

    Element** mElements;             // Connected elements
    Element** mNeighbors;            // Face neighbors (may contain nullptrs)
    Facet** mFacets;                 // Boundary facets
    ControlPoint** mControlPoints;   // For spline-based elements
};
```

**Key methods (virtual, implemented in derived classes):**
```cpp
virtual uint number_of_nodes() const;
virtual uint number_of_corner_nodes() const;
virtual uint number_of_edges() const;
virtual uint number_of_faces() const;
virtual uint number_of_facets() const;
virtual ElementType type() const;

virtual Node* node(uint aIndex);
virtual Edge* edge(uint aIndex);
virtual Face* face(uint aIndex);

void insert_node(Node* aNode, uint aIndex);
void insert_edge(Edge* aEdge, uint aIndex);
void insert_face(Face* aFace, uint aIndex);
```

**Element hierarchy:**
```
Element (base)
├── ElementTemplate< 2, 2, 1, 0, 0 >     // LINE2
├── ElementTemplate< 3, 3, 3, 3, 1 >     // TRI3
├── ElementTemplate< 4, 4, 4, 4, 1 >     // QUAD4
├── ElementTemplate< 4, 4, 6, 4, 4 >     // TET4
├── ElementTemplate< 8, 8, 12, 6, 6 >    // HEX8
└── ... (one instantiation per ElementType)
```
Template parameters are ( nodes, corner nodes, edges, facets, faces ); the `ElementType` is
recovered from the count tuple by the `type()` specialization in each `cl_Element_<TYPE>.hpp`.

### Block Class

**File:** `cl_Block.{hpp,cpp}`

Groups elements of the same type (typically one material/subdomain):

```cpp
class Block {
    id_t mID;                      // Block identifier
    Cell<Element*> mElements;      // Elements in this block
    DomainType mDomainType;        // Default, ThinShell, etc.
    string mLabel;                 // User-assigned name
    bool mHasEdges;                // Nédélec edge flags
    bool mHasFaces;                // Nédélec face flags
    real mThickness;               // For thin-shell blocks
};
```

**Usage:**
```cpp
mesh::Block* block = tMesh->block(blockID);
Cell<mesh::Element*>& elems = block->elements();
ElementType type = block->element_type();  // Type of first element
```

### SideSet Class

**File:** `cl_SideSet.{hpp,cpp}`

Groups boundary facets (surface elements):

```cpp
class SideSet {
    id_t mID;                      // SideSet identifier
    Cell<Facet*> mFacets;          // Boundary facets
    Cell<Node*> mNodes;            // Nodes on boundary (computed)
    DomainType mDomainType;        // Default, Boundary, Cut, etc.
    string mLabel;                 // User-assigned name
    bool mIsHidden;                // Hide from Exodus output
};
```

**Usage:**
```cpp
mesh::SideSet* sideset = tMesh->sideset(sidesetID);
Cell<mesh::Facet*>& facets = sideset->facets();
sideset->collect_nodes();  // Populate mNodes from facets
Cell<mesh::Node*>& nodes = sideset->nodes();
```

### Edge Class

**File:** `cl_Edge.{hpp,cpp}`

Topological edge entity for Nédélec-type (edge-based) finite elements:

```cpp
class Edge : public Vertex {
    // Inherits connectivity from Vertex
    // Used for H(curl) DOF allocation
};
```

**Key methods:**
```cpp
EntityType entity_type() const;  // Returns EntityType::EDGE
size_t memory() const;           // Memory footprint
```

**Usage:**
- Created by `Mesh::create_edges()` for blocks/sidesets requiring edge-based DOFs
- Used in H(curl) formulations (electromagnetic fields, Nédélec elements)
- Each edge has unique global ID and connectivity to owning elements

**Example:**
```cpp
tMesh->create_edges();
Cell<mesh::Edge*>& edges = tMesh->edges();
mesh::Edge* edge = tMesh->edge(edgeID);

// Edge connectivity
mesh::Element* elem = tMesh->element(elemID);
uint nEdges = elem->number_of_edges();
for (uint i = 0; i < nEdges; ++i) {
    mesh::Edge* edge = elem->edge(i);
    bool isPlus = elem->edge_direction(i);  // Tangent orientation
}
```

### Face Class

**File:** `cl_Face.{hpp,cpp}`

Topological face entity for Nédélec-type (face-based) finite elements:

```cpp
class Face : public Vertex {
    // Inherits connectivity from Vertex
    // Used for H(div) DOF allocation
};
```

**Key methods:**
```cpp
EntityType entity_type() const;  // Returns EntityType::FACE
size_t memory() const;           // Memory footprint
```

**Usage:**
- Created by `Mesh::create_faces()` for 3D elements requiring face-based DOFs
- Used in H(div) formulations (fluid flow, Raviart-Thomas elements)
- Each face has unique global ID and connectivity to owning elements

**Example:**
```cpp
tMesh->create_faces();
Cell<mesh::Face*>& faces = tMesh->faces();
mesh::Face* face = tMesh->face(faceID);

// Face connectivity
mesh::Element* elem = tMesh->element(elemID);
uint nFaces = elem->number_of_faces();
for (uint i = 0; i < nFaces; ++i) {
    mesh::Face* face = elem->face(i);
}
```

### Facet Class

**File:** `cl_Facet.{hpp,cpp}`

Wraps an owned lower-dimensional element and links it to the master/slave elements it sits on:

```cpp
class Facet : public Vertex {
    Element* mElement;         // owned facet element (its own node list)
    Element* mMaster;          // element this facet sits on
    Element* mSlave;           // neighbor element (nullptr on boundary)
    suint mMasterFaceID;       // local facet index on master
    suint mSlaveFaceID;        // local facet index on slave
};
```

**Key properties:**
- **Master element:** The element this facet references
- **Slave element:** Neighbor across this facet (null on domain boundary)
- **Local indices:** Facet numbering within master/slave elements
- **Orientation:** Computed to ensure consistent normal directions

**Important distinction:**
- `Edge` and `Face` are **topological entities** with global IDs (for DOF allocation)
- `Facet` is a **geometric reference** to element faces/edges (for boundary conditions)

### ThinShell Class

**File:** `cl_ThinShell.{hpp,cpp}`

Data container linking thin-shell element blocks to boundary sidesets for layered conductor modeling:

```cpp
class ThinShell {
    SideSet* mSideSet;              // Boundary sideset (facets)
    Cell<Block*> mBlocks;            // Element blocks (one per layer)
    Vector<real> mThicknesses;       // Thickness per layer [m]
    Cell<string> mMaterials;         // Material name per layer
};
```

**Key properties:**
- **Layered structure:** Each block represents one physical layer
- **Sideset coupling:** Facets link thin-shell elements to adjacent "air" elements
- **Physical parameters:** Thickness and material per layer
- **Cohomology-aware:** Node ordering supports discontinuous scalar potential

**Key methods:**
```cpp
Cell<Facet*>& facets();                      // Boundary facets
Cell<Block*>& blocks();                      // Element blocks (layers)
const Vector<real>& thicknesses() const;     // Layer thicknesses [m]
void set_thicknesses(const Vector<real>&);   // Set thicknesses
void set_materials(Cell<string>&);           // Set material names
const Cell<string>& materials() const;       // Get material names
ElementType element_type() const;            // Element type (QUAD4TS, PENTA6TS, etc.)
id_t id() const;                             // ID (from sideset)
const string& label() const;                 // Label (from sideset)
```

**Usage example:**
```cpp
// Create thin-shell configuration
SideSet* sideset = tMesh->sideset(sidesetID);
// the constructor takes both sidesets; the ghost may be nullptr
ThinShell* shell = new ThinShell(sideset, ghostSideset);

// Set layer thicknesses (e.g., 3 layers: 100 µm, 50 µm, 100 µm)
Vector<real> thicknesses = {100e-6, 50e-6, 100e-6};
shell->set_thicknesses(thicknesses);

// Set material names
Cell<string> materials = {"Copper", "Insulation", "Copper"};
shell->set_materials(materials);

// Access blocks and facets
Cell<Block*>& layers = shell->blocks();
Cell<Facet*>& facets = shell->facets();

// Add to mesh
tMesh->thin_shells().push(shell);
```

**Typical workflow:**
1. Define sideset representing thin-shell surface
2. Create `ThinShell` object from sideset
3. Set layer thicknesses and materials
4. FEM kernel creates element blocks for each layer
5. Cohomology module creates cuts if needed for multiply-connected shells
6. Facets couple thin-shell elements to surrounding domain

**Physical context:** Used for modeling superconducting tapes, cables, and layered conductors in electromagnetic simulations. Thickness is typically orders of magnitude smaller than in-plane dimensions, justifying 2D shell approximation.

**See also:**
- Thin-shell element types: `QUAD4TS`, `QUAD9TS`, `PENTA6TS`, `PENTA18TS`
- Advanced Features → Thin-Shell Elements section (line 1462)
- Homology module documentation (cohomology cuts for shells)

---

## Mesh Construction and I/O

### Reading Mesh Files

**Gmsh format (.msh):**
```cpp
#include "cl_Mesh.hpp"

// Read Gmsh file (ASCII MSH 2.2 or 4.1); the reader finalizes the mesh
Mesh* tMesh = new Mesh("mesh.msh");
```

**HDF5 format (.hdf5):**
```cpp
// Read HDF5 mesh (BELFEM native format); the reader finalizes the mesh
Mesh* tMesh = new Mesh("mesh.hdf5");
```

**Parallel mesh loading:**
```cpp
#ifdef BELFEM_MPI
// Only the master proc reads the file; the other ranks receive the
// dimension count and stay empty until Distributor::run()
Mesh* tMesh = new Mesh("mesh.msh",
    0,       // master proc
    true,    // compute connectivities
    true);   // parallel mode

// Partition and distribute
tMesh->partition(comm_size());
tMesh->finalize();
#endif
```

### Writing Mesh Files

**HDF5 (BELFEM native, read/write):**
```cpp
// Save mesh topology and enrichment (no field data; fields go to the
// restart file via Mesh::save_fields(hid_t))
tMesh->save("output.hdf5");

// Save with element connectivity (larger file)
tMesh->save("output.hdf5");   // one argument; there is no field-saving flag
```

**Exodus II (Sandia format, write-only):**
```cpp
#include "cl_Mesh_ExodusWriter.hpp"

mesh::ExodusWriter writer(tMesh);
writer.save("output.exo");

// Time series: ExodusWriter::save() always recreates the file, so write
// one numbered file per step via the "e-s" extension
tMesh->set_time_step(1);
tMesh->time_stamp() = 0.0;
tMesh->save("output.e-s");   // -> output.e-s.00001

tMesh->set_time_step(2);
tMesh->time_stamp() = 0.01;
tMesh->save("output.e-s");   // -> output.e-s.00002
```

**VTK (ParaView format, write-only):**
```cpp
#include "cl_Mesh_VtkWriter.hpp"

mesh::VtkWriter writer("output.vtk", tMesh);  // writes on construction
// or simply: tMesh->save("output.vtk");
```

### Creating Meshes Programmatically

**Empty mesh:**
```cpp
// Create empty 3D mesh
Mesh* tMesh = new Mesh(3);  // 3 dimensions

// Manually populate nodes and elements
// (advanced, see tensor mesh example)
```

**Tensor product mesh:**
```cpp
// 2D structured mesh: 11x11 nodes, quadratic (order 2)
uint order = 2;
Vector<index_t> numNodes = {11, 11};
Vector<real> step = {0.1, 0.1};
Vector<real> origin = {0.0, 0.0};

Mesh* tMesh = new Mesh(order, numNodes, step, origin);
tMesh->finalize();

// Access nodes by (i, j) index
mesh::Node* node = tMesh->node(5, 5);

// Access elements by (i, j) index
mesh::Element* elem = tMesh->element(2, 3);
```

---

## Mesh Finalization

**Critical step after mesh construction:**

```cpp
tMesh->finalize();
```

**What finalize() does:**
1. **Compute element indices:** Assigns continuous 0-based indices to all elements
2. **Update node indices:** Assigns continuous 0-based indices to all nodes
3. **Set block/sideset IDs:** Propagates block IDs to elements, sideset IDs to facets
4. **Compute facet orientations:** Ensures consistent normals on sidesets
5. **Build connectivity maps:** Creates ID→entity maps for fast lookup

`checksum()` is computed lazily on first call; hanging entities are collected by `collect_hanging_basis()`, not by `finalize()`.

**Advanced finalization control:**
```cpp
// Disable facet orientation (if you handle it manually)
tMesh->set_compute_facet_orientation_flag(false);
tMesh->finalize();
```

**Undoing finalization (rare):**
```cpp
tMesh->unfinalize();  // Resets selected connectivity flags; indices are untouched
// Modify mesh topology...
tMesh->finalize();    // Re-finalize
```

---

## Element Types and Topology

### Complete ElementType Enum Reference

**File:** `Mesh_Enums.hpp`

**Complete table of all supported element types:**

| Enum Value | Nodes | Geometry | Order | Description |
|------------|-------|----------|-------|-------------|
| `EMPTY` | 0 | — | — | Empty/placeholder element |
| `VERTEX` | 1 | Point | 0 | Single point (constant) |
| **1D Elements** | | | | |
| `LINE2` | 2 | Line | 1 | Linear line segment |
| `LINE3` | 3 | Line | 2 | Quadratic line |
| `LINE4` | 4 | Line | 3 | Cubic line |
| `LINE5` | 5 | Line | 4 | Quartic line |
| `LINE6` | 6 | Line | 5 | Quintic line — **enum only, not built by `ElementFactory`** |
| **2D Triangles** | | | | |
| `TRI3` | 3 | Triangle | 1 | Linear triangle |
| `TRI6` | 6 | Triangle | 2 | Quadratic triangle |
| `TRI10` | 10 | Triangle | 3 | Cubic triangle |
| `TRI15` | 15 | Triangle | 4 | Quartic triangle |
| `TRI21` | 21 | Triangle | 5 | Quintic triangle — **enum only, not built by `ElementFactory`** |
| **2D Quadrilaterals** | | | | |
| `QUAD4` | 4 | Quad | 1 | Bilinear quadrilateral |
| `QUAD8` | 8 | Quad | 2 | Serendipity quadratic quad |
| `QUAD9` | 9 | Quad | 2 | Biquadratic quad (with center node) |
| `QUAD16` | 16 | Quad | 3 | Bicubic quadrilateral |
| `QUAD4TS` | 4 | Quad | 1 | Thin-shell bilinear quad |
| `QUAD9TS` | 9 | Quad | 2 | Thin-shell biquadratic quad |
| **3D Tetrahedra** | | | | |
| `TET4` | 4 | Tet | 1 | Linear tetrahedron |
| `TET10` | 10 | Tet | 2 | Quadratic tetrahedron |
| `TET20` | 20 | Tet | 3 | Cubic tetrahedron |
| `TET35` | 35 | Tet | 4 | Quartic tetrahedron |
| **3D Hexahedra** | | | | |
| `HEX8` | 8 | Hex | 1 | Trilinear hexahedron |
| `HEX20` | 20 | Hex | 2 | Serendipity triquadratic hex |
| `HEX27` | 27 | Hex | 2 | Triquadratic hex (with center node) |
| `HEX64` | 64 | Hex | 3 | Tricubic hexahedron |
| **3D Prisms (Wedges)** | | | | |
| `PENTA6` | 6 | Prism | 1 | Linear prism |
| `PENTA15` | 15 | Prism | 2 | Serendipity quadratic prism |
| `PENTA18` | 18 | Prism | 2 | Quadratic prism (with face centers) |
| `PENTA6TS` | 6 | Prism | 1 | Thin-shell linear prism |
| `PENTA18TS` | 18 | Prism | 2 | Thin-shell quadratic prism |
| **3D Pyramids** | | | | |
| `PYRA5` | 5 | Pyramid | 1 | Linear pyramid |
| `PYRA13` | 13 | Pyramid | 2 | Serendipity quadratic pyramid |
| `PYRA14` | 14 | Pyramid | 2 | Quadratic pyramid (with base center) |
| `UNDEFINED` | — | — | — | Undefined/unknown element type |

**Notes:**

1. **Thin-Shell Elements (TS suffix):** Special variants for cohomology-aware layered structures. Nodes ordered with top/bottom pairs for discontinuous DOF support.
2. **Serendipity vs Full:** Serendipity elements (QUAD8, HEX20, PENTA15, PYRA13) omit interior nodes for same polynomial order.
3. **Gmsh Compatibility:** Most element types follow Gmsh numbering (enum values 1-30, 92); QUAD16 = 32 deviates from Gmsh's 36. Thin-shell/beam types (103, 105, 106, 110, 118, 125) are BELFEM extensions.
4. **Node Ordering:** All elements follow Gmsh node ordering conventions (see Gmsh documentation for detailed node numbering).

### InterpolationOrder Enum

**File:** `Mesh_Enums.hpp`

| Enum Value | Description | Typical Elements |
|------------|-------------|------------------|
| `CONSTANT` | Order 0 (piecewise constant) | VERTEX |
| `LINEAR` | Order 1 (piecewise linear) | LINE2, TRI3, QUAD4, TET4, HEX8, PENTA6, PYRA5 |
| `QUADRATIC` | Order 2 (quadratic, includes face/volume centers) | LINE3, TRI6, QUAD9, TET10, HEX27, PENTA18, PYRA14 |
| `SERENDIPITY` | Order 2 (quadratic, edge nodes only) | QUAD8, HEX20, PENTA15, PYRA13 |
| `CUBIC` | Order 3 (cubic) | LINE4, TRI10, QUAD16, TET20, HEX64 |
| `QUARTIC` | Order 4 (quartic) | LINE5, TRI15, TET35 |
| `QUINTIC` | Order 5 (quintic) | LINE6, TRI21 |
| `UNDEFINED` | Unknown/invalid order | — |

**Key distinction:** `SERENDIPITY` elements have the same polynomial order as `QUADRATIC` but omit interior nodes (face/volume centers), reducing DOF count while maintaining boundary accuracy.

### InterpolationType Enum

**File:** `Mesh_Enums.hpp`

| Enum Value | Description | Use Case |
|------------|-------------|----------|
| `LAGRANGE` | Standard Lagrange polynomials | Default for most FEM elements |
| `HERMITE` | Hermite polynomials (C¹ continuous) | Beam elements, plate/shell elements |
| `BERNSTEIN` | Bernstein basis (numerically stable) | High-order elements, isogeometric analysis |
| `BubbleEdge0` | Edge bubble function (edge 0) | Stabilization, enrichment |
| `BubbleEdge1` | Edge bubble function (edge 1) | Stabilization, enrichment |
| `BubbleEdge2` | Edge bubble function (edge 2) | Stabilization, enrichment |
| `BubbleFace0` | Face bubble function (face 0) | Stabilization, enrichment |
| `BubbleFace1` | Face bubble function (face 1) | Stabilization, enrichment |
| `BubbleFace2` | Face bubble function (face 2) | Stabilization, enrichment |
| `BubbleFace3` | Face bubble function (face 3) | Stabilization, enrichment |
| `UNEFINED` | Unknown/invalid type | — |

**Note:** Bubble functions vanish on element boundaries and are used for local enrichment in mixed FEM formulations (e.g., MINI element for incompressible flow).

### Supported Element Types (Examples)

**1D Elements:**
```cpp
ElementType::LINE2   // 2-node linear line
ElementType::LINE3   // 3-node quadratic line
ElementType::LINE4   // 4-node cubic line
ElementType::LINE5   // 5-node quartic line
```

**2D Elements:**
```cpp
// Triangles
ElementType::TRI3    // 3-node linear
ElementType::TRI6    // 6-node quadratic
ElementType::TRI10   // 10-node cubic
ElementType::TRI15   // 15-node quartic

// Quadrilaterals
ElementType::QUAD4   // 4-node bilinear
ElementType::QUAD8   // 8-node serendipity quadratic
ElementType::QUAD9   // 9-node biquadratic
ElementType::QUAD16  // 16-node bicubic
```

**3D Elements:**
```cpp
// Tetrahedra
ElementType::TET4    // 4-node linear
ElementType::TET10   // 10-node quadratic
ElementType::TET20   // 20-node cubic
ElementType::TET35   // 35-node quartic

// Hexahedra
ElementType::HEX8    // 8-node trilinear
ElementType::HEX20   // 20-node serendipity triquadratic
ElementType::HEX27   // 27-node triquadratic
ElementType::HEX64   // 64-node tricubic

// Prisms (wedges)
ElementType::PENTA6  // 6-node linear prism
ElementType::PENTA15 // 15-node quadratic prism
ElementType::PENTA18 // 18-node quadratic prism

// Pyramids
ElementType::PYRA5   // 5-node linear pyramid
ElementType::PYRA13  // 13-node quadratic pyramid
ElementType::PYRA14  // 14-node quadratic pyramid
```

**Thin-Shell Elements (special variants):**
```cpp
ElementType::QUAD4TS   // Thin-shell quad (cohomology-aware)
ElementType::QUAD9TS   // Thin-shell quad, quadratic
ElementType::PENTA6TS  // Thin-shell prism
ElementType::PENTA18TS // Thin-shell prism, quadratic
```

### Element Topology Queries

```cpp
mesh::Element* elem = tMesh->element(elemID);

// Type and geometry
ElementType type = elem->type();
uint dim = elem->dimension();  // 1, 2, or 3

// Node connectivity
uint nNodes = elem->number_of_nodes();
uint nCorners = elem->number_of_corner_nodes();
for (uint i = 0; i < nNodes; ++i) {
    mesh::Node* node = elem->node(i);
    real x = node->x();
    real y = node->y();
    real z = node->z();
}

// Facet connectivity
uint nFacets = elem->number_of_facets();
for (uint i = 0; i < nFacets; ++i) {
    Cell<mesh::Node*> facetNodes;
    elem->get_nodes_of_facet(i, facetNodes);
}

// Edge/face connectivity (if created)
if (elem->has_edges()) {
    uint nEdges = elem->number_of_edges();
    for (uint i = 0; i < nEdges; ++i) {
        mesh::Edge* edge = elem->edge(i);
    }
}

if (elem->has_faces()) {
    uint nFaces = elem->number_of_faces();
    for (uint i = 0; i < nFaces; ++i) {
        mesh::Face* face = elem->face(i);
    }
}
```

### Node Ordering Conventions

BELFEM follows **Gmsh node ordering** for all elements:

**Example: QUAD9 (9-node biquadratic quad)**
```
3---6---2
|       |
7   8   5
|       |
0---4---1
```
- Nodes 0-3: Corner nodes (CCW from origin)
- Nodes 4-7: Edge midpoints
- Node 8: Face center

**Example: HEX20 (20-node serendipity triquadratic hex)**
```
      7------18------6
     /|             /|
   19 |           17 |
   /  15          /  14
  4------16------5   |
  |   |          |   |
  |   3------10--|---2
 12  /          13  /
  | 11           | 9
  |/             |/
  0-------8------1
```
- Nodes 0-7: Corner nodes
- Nodes 8-19: Edge midpoints

See Gmsh documentation for complete node ordering tables.

---

## Mesh Connectivity

### Connectivity Types

BELFEM supports 36 connectivity types defined in the `Connectivity` enum (`Mesh_Enums.hpp`):

**Complete Connectivity Enum Table:**

| Enum Value                     | Integer | Description |
|--------------------------------|---------|-------------|
| `Compute`                      | 0 | Trigger automatic connectivity computation |
| **Node Connectivities**        | | |
| `NodeToVertex`                 | 1 | Nodes to graph vertices |
| `NodeToNode`                   | 2 | Node neighbors via shared elements |
| `NodeToEdge`                   | 3 | Nodes to edges containing them |
| `NodeToFace`                   | 4 | Nodes to faces containing them |
| `NodeToFacet`                  | 5 | Nodes to boundary facets containing them |
| `NodeToElement`                | 6 | Nodes to elements containing them |
| **Edge Connectivities**        | | |
| `EdgeToVertex`                 | 7 | Edges to graph vertices |
| `EdgeToNode`                   | 8 | Edges to their endpoint nodes |
| `EdgeToEdge`                   | 9 | Edge neighbors via shared elements |
| `EdgeToFace`                   | 10 | Edges to faces containing them |
| `EdgeToFacet`                  | 11 | Edges to facets containing them |
| `EdgeToElement`                | 12 | Edges to elements containing them |
| **Face Connectivities**        | | |
| `FaceToVertex`                 | 13 | Faces to graph vertices |
| `FaceToNode`                   | 14 | Faces to their corner/edge nodes |
| `FaceToEdge`                   | 15 | Faces to edges bounding them |
| `FaceToFace`                   | 16 | Face neighbors via shared elements |
| `FaceToFacet`                  | 17 | Faces to facets |
| `FaceToElement`                | 18 | Faces to elements containing them |
| **Facet Connectivities**       | | |
| `FacetToVertex`                | 19 | Facets to graph vertices |
| `FacetToNode`                  | 20 | Facets to their nodes |
| `FacetToEdge`                  | 21 | Facets to edges on them |
| `FacetToFace`                  | 22 | Facets to topological faces |
| `FacetToFacet`                 | 23 | Facet neighbors |
| `FacetToElement`               | 24 | Facets to master/slave elements |
| **Element Connectivities**     | | |
| `ElementToVertex`              | 25 | Elements to graph vertices |
| `ElementToNode`                | 26 | Elements to their nodes (intrinsic) |
| `ElementToEdge`                | 27 | Elements to edges on them |
| `ElementToFace`                | 28 | Elements to faces on them |
| `ElementToFacet`               | 29 | Elements to boundary facets |
| `ElementToElement`             | 30 | Element neighbors via shared entities |
| **Specialized Connectivities** | | |
| `TsElementToTsElement`         | 31 | Thin-shell element neighbors |
| `ShellToShell`                 | 32 | Shell neighbors via shared edges |
| `ControlPointToControlPoint`   | 33 | Spline control point neighbors |
| `ElementToControlPoint`        | 34 | Elements to control points (B-spline) |
| `ControlPointToElement`        | 35 | Control points to elements |
| `UNDEFINED`                    | 36 | Undefined connectivity type |

**Intrinsic vs Computed:**
- **Intrinsic:** `ElementToNode` (stored directly in element)
- **Computed:** All others (built on-demand, cached)

**Most commonly used connectivities:**
- `NodeToElement`: Which elements contain each node
- `ElementToElement`: Element neighbors via shared facets
- `NodeToNode`: Node neighbors via shared elements
- `EdgeToElement`: Which elements contain each edge (for Nédélec)
- `FaceToElement`: Which elements contain each face (for Raviart-Thomas)

### Computing Connectivities

**Automatic (via finalize):**
```cpp
Mesh* tMesh = new Mesh("mesh.msh",
    0,       // master proc
    true);   // compute connectivities = true (default)
tMesh->finalize();
```

**Manual (on-demand):**
```cpp
Mesh* tMesh = new Mesh("mesh.msh",
    0,       // master proc
    false);  // compute connectivities = false

tMesh->finalize();

// Later, check whether a connectivity was built (finalize() builds them
// eagerly; nothing is computed on access)
if (!tMesh->test_connectivity(Connectivity::NodeToElement)) {
    // NodeToElement was not built; node->elements() is not populated
}
```

### Using Connectivities

**Node-to-Element:**
```cpp
mesh::Node* node = tMesh->node(nodeID);

// Access connected elements
uint nElems = node->number_of_elements();
for (uint i = 0; i < nElems; ++i) {
    mesh::Element* elem = node->element(i);
}
```

**Element-to-Element:**
```cpp
mesh::Element* elem = tMesh->element(elemID);

// Access all connected elements (via shared nodes/edges/faces)
uint nElems = elem->number_of_elements();
for (uint i = 0; i < nElems; ++i) {
    mesh::Element* neighbor = elem->element(i);
}
```

**Element neighbors (via facets):**
```cpp
mesh::Element* elem = tMesh->element(elemID);

// First, populate neighbors
tMesh->populate_element_neighbors();

// Access facet neighbors
uint nFacets = elem->number_of_facets();
for (uint i = 0; i < nFacets; ++i) {
    mesh::Element* neighbor = elem->neighbor(i);
    if (neighbor == nullptr) {
        // This facet is on domain boundary
    } else {
        // This facet is shared with neighbor
    }
}
```

---

## Mesh Partitioning (MPI Parallelization)

### Graph-Based Partitioning

**Using METIS (default):**
```cpp
#ifdef BELFEM_METIS
// Partition mesh into comm_size() subdomains
tMesh->partition(comm_size());

// Each element now has owner() == MPI rank
for (mesh::Element* elem : tMesh->elements()) {
    proc_t owner = elem->owner();
}

// Each node has owner() == owning proc
for (mesh::Node* node : tMesh->nodes()) {
    proc_t owner = node->owner();
}
#endif
```

**Partitioning options:**
```cpp
// Partition with selected blocks only
Vector<id_t> selectedBlocks = {1, 2, 3};
tMesh->partition(comm_size(), selectedBlocks);

// Partition with blocks and sidesets
Vector<id_t> selectedSideSets = {10, 20};
tMesh->partition(comm_size(), selectedBlocks, selectedSideSets);

// Force continuous partitions (eliminate small disconnected regions)
tMesh->partition(comm_size(),
    selectedBlocks,
    selectedSideSets,
    true,   // set proc owners
    true);  // force continuous partitions
```

### Mesh Distribution (MPI)

**After partitioning, distribute mesh to processors:**
```cpp
#ifdef BELFEM_MPI
#include "cl_Mesh_Distributor.hpp"

// Master proc (rank 0) has full mesh, others empty
Mesh* tMesh = nullptr;
if (comm_rank() == 0) {
    tMesh = new Mesh("mesh.msh");
    tMesh->partition(comm_size());
}

// Distribute: Each proc gets its partition + ghost layers
mesh::Distributor distributor(tMesh);
distributor.run();                       // every rank

Mesh* localMesh = tMesh;                 // root keeps the mesh it already has
if (comm_rank() != 0) {
    // workers only -- partial_mesh() is a BELFEM_ERROR on root
    localMesh = distributor.partial_mesh();
}

// Now each proc has localMesh with:
// - Elements owned by this proc
// - Ghost elements shared with neighbors
// - Nodes owned or shared

// run() has already finalized the worker meshes, so only the root's own
// mesh still needs it ( cl_Mesh_Distributor.cpp:248 ).
if (comm_rank() == 0) {
    localMesh->finalize();
}

// Do NOT delete tMesh on root: localMesh IS tMesh there, and it is still in use.
#endif
```

### Ghost Layers and Ownership

**After distribution:**
- Each proc stores elements with `owner() == comm_rank()`
- Plus ghost elements: `owner() != comm_rank()` (shared with neighbors)
- Nodes are owned by lowest-rank proc containing them

**Checking ownership:**
```cpp
for (mesh::Element* elem : localMesh->elements()) {
    if (elem->owner() == comm_rank()) {
        // This element is owned by this proc
    } else {
        // This is a ghost element
    }
}

for (mesh::Node* node : localMesh->nodes()) {
    if (node->owner() == comm_rank()) {
        // This node is owned by this proc
    } else {
        // This is a ghost node
    }
}
```

---

## Fields and Data Management

### Creating Fields

**Node-based fields (most common):**
```cpp
// Create node field (automatically sized to number of nodes)
Vector<real>& temperature = tMesh->create_field(
    "Temperature",           // Field name
    EntityType::NODE,        // Entity type
    0);                      // Field ID (0 = auto-assign)

// Initialize field data
temperature.fill(300.0);

// Access field data later
Vector<real>& T = tMesh->field_data("Temperature");
```

**Element-based fields:**
```cpp
Vector<real>& stress = tMesh->create_field(
    "Stress",
    EntityType::ELEMENT,
    0);

// Size = number of elements
for (index_t i = 0; i < tMesh->number_of_elements(); ++i) {
    stress(i) = compute_element_stress(tMesh->elements()(i));
}
```

**Edge/Face fields (for Nédélec elements):**
```cpp
// Create edge field (for H(curl) DOFs)
tMesh->create_edges();  // First create edges
Vector<real>& edgeField = tMesh->create_field("EdgeFlux", EntityType::EDGE);

// Create face field (for H(div) DOFs)
tMesh->create_faces();  // First create faces
Vector<real>& faceField = tMesh->create_field("FaceFlux", EntityType::FACE);
```

### Global Variables

**Scalar mesh-level data:**
```cpp
// Create global variable
real& current = tMesh->create_global_variable("TotalCurrent", 1000.0);

// Access later
real& I = tMesh->global_variable_data("TotalCurrent");
I = 1500.0;

// Check existence
if (tMesh->global_variable_exists("TotalCurrent")) {
    // ...
}
```

### Field I/O

**Fields are not part of the mesh file.** `Mesh::save("*.hdf5")` writes topology and enrichment only; fields, global variables and the time cursor go to the separate restart file through `Mesh::save_fields(hid_t)` / `load_fields(hid_t)` (see `bfm_file_format.md`).
```cpp
// Create fields
tMesh->create_field("Temperature", EntityType::NODE);
tMesh->create_field("Pressure", EntityType::NODE);

// Save mesh (no fields)
tMesh->save("output.hdf5");

// A reloaded mesh has no fields until load_fields() is called
Mesh* loadedMesh = new Mesh("output.hdf5");
```

**Exodus time series:**
```cpp
#include "cl_Mesh_ExodusWriter.hpp"

// one numbered file per step; ExodusWriter::save() always recreates its file
for (uint step = 1; step <= 100; ++step) {
    tMesh->set_time_step(step);
    tMesh->time_stamp() = step * 0.01;

    // Update field data
    Vector<real>& T = tMesh->field_data("Temperature");
    // ... compute T ...

    tMesh->save("output.e-s");   // -> output.e-s.00001, .00002, ...
}
```

---

## Advanced Features

### Edges and Faces (Nédélec Elements)

**For H(curl) and H(div) finite elements:**

```cpp
// Create edges (for edge-based DOFs)
tMesh->create_edges(
    true,                 // print statistics
    {1, 2},               // Nédélec blocks (require edges)
    {10, 20},             // Nédélec sidesets (require edges)
    true);                // create edges on all sidesets

// Access edges
Cell<mesh::Edge*>& edges = tMesh->edges();
mesh::Edge* edge = tMesh->edge(edgeID);

// Edge connectivity
mesh::Element* elem = tMesh->element(elemID);
uint nEdges = elem->number_of_edges();
for (uint i = 0; i < nEdges; ++i) {
    mesh::Edge* edge = elem->edge(i);

    // Edge direction (tangent orientation)
    bool isPlus = elem->edge_direction(i);
}

// Create faces (for face-based DOFs)
tMesh->create_faces(
    true,                 // print statistics
    {1, 2},               // Nédélec blocks
    {10, 20});            // Nédélec sidesets

// Access faces
Cell<mesh::Face*>& faces = tMesh->faces();
uint nFaces = elem->number_of_faces();
for (uint i = 0; i < nFaces; ++i) {
    mesh::Face* face = elem->face(i);
}
```

**Finalize edges/faces separately:**
```cpp
tMesh->finalize();        // Finalize main mesh
tMesh->create_edges();    // Create edges
tMesh->finalize_edges();  // Finalize edge connectivity
```

### Thin-Shell Elements

**For layered conductor modeling (cohomology):**

```cpp
#include "cl_ThinShell.hpp"

// Define thin-shell configuration: the constructor takes both sidesets;
// the ghost may be nullptr. Layers are attached afterwards.
mesh::ThinShell* shell = new mesh::ThinShell(sideset, ghostSideset);
shell->blocks().push(layerBlock);          // one entry per layer
shell->set_thicknesses(thicknesses);       // Vector<real>, one per layer [m]

// Add to mesh
tMesh->thin_shells().push(shell);

// Set block domain types
for (id_t blockID : {1, 2, 3}) {
    mesh::Block* block = tMesh->block(blockID);
    block->set_domain_type(DomainType::ThinShell);
    block->set_thickness(0.001);
}
```

**Thin-shell elements have special node ordering for cohomology:**
- Nodes on "top" surface (positive normal)
- Nodes on "bottom" surface (negative normal)
- Allows discontinuous scalar potential across shell

### Tensor Product Meshes

**Structured Cartesian grids:**

```cpp
// 3D tensor mesh: 21x21x11 nodes, cubic elements (order 3)
uint order = 3;
Vector<index_t> numNodes = {21, 21, 11};
Vector<real> step = {0.05, 0.05, 0.1};     // Grid spacing
Vector<real> origin = {0.0, 0.0, 0.0};     // Origin

Mesh* tMesh = new Mesh(order, numNodes, step, origin);
tMesh->finalize();

// Access via structured indices
mesh::Node* node = tMesh->node(10, 10, 5);  // i, j, k
mesh::Element* elem = tMesh->element(5, 5, 2);

// Access tensor configuration
const TensorMeshConfig* config = tMesh->tensorconf();
index_t numElemsX = config->num_elements(0);
```

**Tensor mesh properties:**
- Nodes and elements accessible by (i, j) or (i, j, k) indices
- Automatically creates blocks for each element type
- Fast structured-grid operations
- Ideal for image-based meshing, regular domains

### Node Duplication (Cohomology Cuts)

**For multiply-connected domain topological cuts:**

```cpp
// Duplicate node for cohomology cut
mesh::Node* original = tMesh->node(nodeID);
original->allocate_duplicate_container(2);  // Prepare for 2 duplicates

mesh::Node* dup1 = new mesh::Node(newID1, original->x(), original->y(), original->z());
mesh::Node* dup2 = new mesh::Node(newID2, original->x(), original->y(), original->z());

dup1->set_original(original);
dup2->set_original(original);

original->add_duplicate(dup1);
original->add_duplicate(dup2);

// Access duplicates
uint nDups = original->number_of_duplicates();
for (uint i = 0; i < nDups; ++i) {
    mesh::Node* dup = original->duplicate(i);
}

// Check if node is duplicate
if (dup1->is_duplicate()) {
    mesh::Node* orig = dup1->original();  // Points back to original
}
```

**Typical use:** Homology module creates duplicates for cohomology cuts, allowing discontinuous DOFs across cut surfaces.

### Curved Elements

**Flag curved elements (for geometry-exact integration):**

```cpp
tMesh->flag_curved_elements();

for (mesh::Element* elem : tMesh->elements()) {
    if (elem->is_curved()) {
        // Use higher-order integration
    } else {
        // Linear geometry, use standard integration
    }
}
```

**Manually set curved flag:**
```cpp
mesh::Element* elem = tMesh->element(elemID);
elem->set_curved_flag();    // Mark as curved
elem->unset_curved_flag();  // Mark as straight-sided
```

---

## Common Usage Patterns

### Pattern 1: Reading and Querying Mesh

```cpp
#include "cl_Mesh.hpp"

// Load mesh
Mesh* tMesh = new Mesh("mesh.msh");
tMesh->finalize();

// Query mesh properties
uint nDims = tMesh->number_of_dimensions();
index_t nNodes = tMesh->number_of_nodes();
index_t nElems = tMesh->number_of_elements();
uint nBlocks = tMesh->number_of_blocks();
uint nSideSets = tMesh->number_of_sidesets();

message(InfoLevel::Default, "Mesh has %lu nodes, %lu elements, %u blocks, %u sidesets",
    (luint)nNodes, (luint)nElems, nBlocks, nSideSets);

// Iterate over blocks
for (mesh::Block* block : tMesh->blocks()) {
    message(InfoLevel::Default, "Block %lu: %lu elements of type %s",
        (luint)block->id(),
        (luint)block->number_of_elements(),
        to_string(block->element_type()).c_str());
}

// Cleanup
delete tMesh;
```

### Pattern 2: Creating Node-Based Field

```cpp
// Create temperature field
Vector<real>& T = tMesh->create_field("Temperature", EntityType::NODE);

// Initialize based on node coordinates
Cell<mesh::Node*>& nodes = tMesh->nodes();
for (index_t i = 0; i < nodes.size(); ++i) {
    mesh::Node* node = nodes(i);
    T(i) = 300.0 + 100.0 * node->x();  // Linear temperature gradient
}

// Save to HDF5
tMesh->save("result.hdf5");
```

### Pattern 3: Boundary Condition Application

```cpp
// Get boundary sideset
mesh::SideSet* boundary = tMesh->sideset(boundaryID);

// Collect nodes on boundary
boundary->collect_nodes();
Cell<mesh::Node*>& boundaryNodes = boundary->nodes();

// Apply Dirichlet BC
for (mesh::Node* node : boundaryNodes) {
    node->flag();  // Mark for BC application
}

// Later, check flags
for (mesh::Node* node : tMesh->nodes()) {
    if (node->is_flagged()) {
        // Apply BC to this node
    }
}

// Unflag after use
tMesh->unflag_all_nodes();
```

### Pattern 4: Element Integration Loop

```cpp
// Integrate over specific block
mesh::Block* block = tMesh->block(blockID);
real integral = 0.0;

for (mesh::Element* elem : block->elements()) {
    // Get element nodes
    uint nNodes = elem->number_of_nodes();
    Matrix<real> X(nNodes, 3);  // Node coordinates
    for (uint i = 0; i < nNodes; ++i) {
        mesh::Node* node = elem->node(i);
        X(i, 0) = node->x();
        X(i, 1) = node->y();
        X(i, 2) = node->z();
    }

    // Compute element contribution (simplified)
    real elemVolume = compute_element_volume(elem, X);
    integral += elemVolume;
}

message(InfoLevel::Default, "Total volume: %.6e", integral);
```

### Pattern 5: Parallel Mesh Loading and Distribution

```cpp
#ifdef BELFEM_MPI
#include "cl_Mesh_Distributor.hpp"

Mesh* localMesh = nullptr;

if (comm_rank() == 0) {
    // Master proc loads and partitions the full mesh
    Mesh* globalMesh = new Mesh("mesh.msh");
    globalMesh->partition(comm_size());

    mesh::Distributor distributor(globalMesh);
    distributor.run();

    // The root does NOT take a partition -- partial_mesh() rejects it
    // ( BELFEM_ERROR, cl_Mesh_Distributor.cpp:2663 ). It keeps the mesh
    // it already has, so globalMesh must NOT be deleted here.
    localMesh = globalMesh;
} else {
    // Worker procs receive their partition
    mesh::Distributor distributor(nullptr);
    distributor.run();
    localMesh = distributor.partial_mesh();
}

// run() has already finalized the worker meshes; finalize the root's own
if (comm_rank() == 0) {
    localMesh->finalize();
}

// Work with local mesh
index_t localNodes = localMesh->number_of_nodes();
index_t localElems = localMesh->number_of_elements();

message(InfoLevel::Default, "Rank %d: %lu nodes, %lu elements",
    comm_rank(), (luint)localNodes, (luint)localElems);

delete localMesh;
#endif
```

### Pattern 6: Mesh Scaling and Coordinate Modification

```cpp
// Scale entire mesh (e.g., convert mm to m)
tMesh->scale_mesh(0.001);

// Or manually modify specific nodes
mesh::Node* node = tMesh->node(nodeID);
real x = node->x();
real y = node->y();
real z = node->z();

// Move node
node->set_coords(x + 0.1, y, z);
```

### Pattern 7: Mesh Statistics and Memory Usage

```cpp
// Compute memory footprint
size_t memBytes = tMesh->memory();
message(InfoLevel::Default, "Mesh memory: %.2f MB", memBytes / 1e6);

// Compute checksum (for verification)
std::size_t checksum = tMesh->checksum();
message(InfoLevel::Default, "Mesh checksum: %zu", checksum);

// Max element order
uint maxOrder = tMesh->max_element_order();
message(InfoLevel::Default, "Max element order: %u", maxOrder);
```

---

## Performance Considerations

### Memory Management

**Mesh entity allocation:**
- Nodes: Stored in `Cell<Node*>`, manually allocated via `new`
- Elements: heap-allocated via `ElementFactory::create_element()`, owned and deleted by their `Block`
- Blocks/SideSets: Allocated via `new`
- Maps: BELFEM `Map<K,V>`, a wrapper around `std::unordered_map`

**Memory ownership:**
- Mesh owns all entities - DO NOT manually delete nodes/elements
- Mesh destructor handles cleanup
- `abstract_nodes()` is a non-owning view; the nodes themselves are appended to `nodes()` by `set_abstract_nodes()` and freed by `~Mesh`

**Memory profiling:**
```cpp
size_t meshMem = tMesh->memory();
message(InfoLevel::Default, "Total mesh memory: %.2f MB", meshMem / 1e6);

// Per-entity memory
for (mesh::Node* node : tMesh->nodes()) {
    meshMem += node->memory();  // Includes connectivity arrays
}
```

### Connectivity Computation

**Eager connectivity:**
- Connectivities computed in `finalize()` when `Connectivity::Compute` is set
- Cached in bitset `mConnectivities`
- Expensive: `O(N * avg_degree)` for node-element, element-element
- Cheap: Intrinsic connectivities (element-node, node coords)

**Optimize connectivity:**
```cpp
// Only compute needed connectivities
Mesh* tMesh = new Mesh("mesh.msh",
    0,       // master
    false);  // don't compute all connectivities

// Manually trigger specific connectivity
// (triggered implicitly by accessing node->elements(), etc.)
```

### Partitioning Performance

**METIS partitioning:**
- Time: `O(N log N)` for N elements
- Memory: Temporary graph structures `~O(N * avg_degree)`

**Distribution overhead:**
- MPI communication: `O(N / P * P)` for N elements, P procs
- Ghost layer: Typically 1-2 element layers (~5-10% overhead)
- Minimize: Use continuous partitions (fewer ghost elements)

### I/O Performance

**File format comparison:**

| Format | Read | Write | Parallel | Compression | Metadata |
|--------|------|-------|----------|-------------|----------|
| Gmsh   | Fast | N/A   | Serial   | ASCII/Binary | Partial  |
| HDF5   | Fast | Fast  | Serial   | No          | Full   |
| Exodus | N/A  | Medium | Serial   | No          | Time series |
| VTK    | N/A  | Medium | Serial   | No          | Visualization |

**Recommendations:**
- **Gmsh** for mesh generation (industry standard)
- **HDF5** for production runs (restart files)
- **Exodus** for time series visualization (ParaView)
- **VTK** for quick visualization checks

**HDF5 options:**
```cpp
// Enable HDF5 compression (slower write, smaller file)
// (requires HDF5 compiled with zlib support)
// TODO: Document compression API once exposed
```

---

## Debugging and Visualization

### Mesh Integrity Checks

**Check mesh validity:**
```cpp
// Verify all entities exist
BELFEM_ASSERT(tMesh->number_of_nodes() > 0, "Mesh has no nodes");
BELFEM_ASSERT(tMesh->number_of_elements() > 0, "Mesh has no elements");

// Verify finalization
BELFEM_ASSERT(tMesh->is_finalized(), "Mesh not finalized");

// Check element connectivity
for (mesh::Element* elem : tMesh->elements()) {
    uint nNodes = elem->number_of_nodes();
    for (uint i = 0; i < nNodes; ++i) {
        mesh::Node* node = elem->node(i);
        BELFEM_ASSERT(node != nullptr,
            "Element %lu has null node at index %u",
            (luint)elem->id(), i);
    }
}

// Verify checksum
std::size_t expected = /* stored checksum */;
std::size_t actual = tMesh->checksum();
BELFEM_ASSERT(expected == actual, "Mesh checksum mismatch");
```

### Quick Visualization

**Export to VTK for ParaView:**
```cpp
#include "cl_Mesh_VtkWriter.hpp"

mesh::VtkWriter writer("debug_mesh.vtk", tMesh);  // writes on construction
// or simply: tMesh->save("debug_mesh.vtk");
```

**Then in ParaView:**
```bash
paraview debug_mesh.vtk
```

### Debugging Node/Element Lookup

**Print node info:**
```cpp
mesh::Node* node = tMesh->node(nodeID);
message(InfoLevel::Verbose, "Node %lu: (%.6e, %.6e, %.6e), owner=%d, index=%lu",
    (luint)node->id(), node->x(), node->y(), node->z(),
    node->owner(), (luint)node->index());
```

**Print element info:**
```cpp
mesh::Element* elem = tMesh->element(elemID);
elem->print();  // Prints type, nodes, connectivity
```

### Flagging System for Debugging

**Mesh entities have flag bits for temporary marking:**
```cpp
// Flag specific nodes
mesh::Node* node = tMesh->node(nodeID);
node->flag();

// Check if flagged
if (node->is_flagged()) {
    message(InfoLevel::Verbose, "Node %lu is flagged", (luint)node->id());
}

// Unflag
node->unflag();

// Bulk operations
tMesh->unflag_all_nodes();
tMesh->unflag_all_elements();

// Flag all nodes in block
mesh::Block* block = tMesh->block(blockID);
block->flag_nodes();
```

---

## See Also

- **Project README:** `../../../README.md`
- **Claude Instructions:** `../../../CLAUDE.md`
- **Documentation Guidelines:** `../../../doc/documentation_guidelines.md`
- **Coding Philosophy:** `../../../doc/coding_philosophy.md`
- **General Documentation:** `../../../doc/README.md`

**Related Modules:**
- **Containers** (`src/containers/`): `Cell`, `Map` used throughout mesh
- **Core** (`src/core/`): `typedefs.hpp` (id_t, index_t, proc_t), `Logger`
- **Communication** (`src/comm/`): MPI utilities for parallel mesh distribution
- **Graph** (`src/math/graph/`): Graph algorithms for partitioning (METIS, SCOTCH)
- **Homology** (`src/homology/`): Cohomology cuts, node duplication, topological analysis
- **FEM Kernel** (`src/fem/kernel/`): Mesh consumers (DofManager, Domain, etc.)

**External Tools:**
- **Gmsh:** Open-source mesh generator (https://gmsh.info/)
- **ParaView:** Visualization for VTK/Exodus files (https://www.paraview.org/)
- **METIS:** Graph partitioning library (https://github.com/KarypisLab/METIS)

---

**End of Mesh Usage Guide**
