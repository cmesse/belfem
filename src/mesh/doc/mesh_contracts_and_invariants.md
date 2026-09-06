# BELFEM Mesh Module - Contracts and Invariants {#mesh_mesh_contracts_and_invariants}

**Module:** `src/mesh`
**Purpose:** Critical contracts, invariants, and rules for safe mesh usage

**Date:** 2026-01-16
**Revision:** 1.0

---

## Overview

This document specifies **hard contracts** that must be respected when using the mesh module. Violating these contracts leads to undefined behavior, data corruption, or crashes.

---

## Finalization Contracts

### 1. MPI Finalization Semantics

> **Contract:**
>
> * **Before distribution:** Only the master rank (`mMasterProc`, default 0) owns valid topology after `finalize()`.
> * **After `Distributor::run()` / `partial_mesh()`:** Each rank owns a *local mesh* that must be finalized independently.
> * **Non-master ranks before distribution:** Accessing full topology on non-master ranks before distribution is **undefined behavior**.

**Implementation detail (from `cl_Mesh.cpp`):**
```cpp
void Mesh::finalize() {
    if (comm_rank() == mMasterProc) {
        // Only master builds topology
        this->update_element_indices();
        this->update_node_indices();
        // ...
    }
    mIsFinalized = true;  // All ranks set flag
    this->set_block_ids();              // runs on every rank
    this->compute_facet_orientations(); // runs on every rank
}
```

**Safe usage:**
```cpp
#ifdef BELFEM_MPI
// Master loads and finalizes
Mesh* globalMesh = nullptr;
if (comm_rank() == 0) {
    globalMesh = new Mesh("mesh.msh");
    globalMesh->finalize();  // Only master has valid topology
}

// Distribute to all ranks
mesh::Distributor dist(globalMesh);
dist.run();                              // every rank
Mesh* localMesh = globalMesh;            // root keeps the full mesh
if (comm_rank() != 0) {
    localMesh = dist.partial_mesh();     // workers only -- root aborts here
}

// run() already finalized the worker meshes ( cl_Mesh_Distributor.cpp );
// the root finalizes the mesh it kept
if (comm_rank() == 0) { localMesh->finalize(); }  // workers were finalized by run()
#endif
```

### 2. `unfinalize()` is NOT the Inverse of `finalize()`

> **Contract:**
>
> `unfinalize()` is **cache invalidation**, not rollback. It:
> * Clears `mElements` and `mFacets` containers (derived lists)
> * Selectively resets connectivity bitset
> * Does **not** restore mesh to pre-finalize state
> * Assumes elements will not change topology

**Implementation detail (from `cl_Mesh.cpp`):**
```cpp
void Mesh::unfinalize() {
    mElements.clear();  // Derived container
    mFacets.clear();    // Derived container
    reset_connectivity(Connectivity::NodeToElement); // ElementToNode is kept
    // ... selective reset, NOT full reconstruction
    mIsFinalized = false;
}
```

> **Warning:**
> `unfinalize()` is intended for **advanced workflows only** (e.g., mesh modification before re-finalization). Do NOT use casually.

**Dangerous example:**
```cpp
// WRONG: Assuming unfinalize() restores pristine state
tMesh->finalize();
// ... use mesh ...
tMesh->unfinalize();  // Does NOT restore original state!
tMesh->partition(comm_size());   // May fail - expects specific state
```

**Safe example:**
```cpp
// Correct: unfinalize() for cache invalidation before topology change
tMesh->finalize();
tMesh->unfinalize();       // Invalidate caches
tMesh->create_edges();     // Modify topology
tMesh->finalize();         // Rebuild caches
```

---

## Ownership and Container Contracts

### 3. Ownership vs. Derived Containers

> **Contract:**
>
> * **Mesh owns all entity objects** (Nodes, Elements, Edges, Faces, Facets).
> * **Some containers are derived views** (e.g., `mElements`, `mFacets`) and may be cleared/rebuilt during `finalize()`/`unfinalize()`.
> * **Do not store raw pointers to containers** across finalize/unfinalize boundaries.

**Primary ownership (from `cl_Mesh.hpp`):**
```cpp
// Primary owner containers (entities allocated here)
Cell<mesh::Node*>    mNodes;      // Mesh owns nodes
Cell<mesh::Edge*>    mEdges;      // Mesh owns edges
Cell<mesh::Face*>    mFaces;      // Mesh owns faces
Cell<mesh::Block*>   mBlocks;     // Mesh owns blocks
Cell<mesh::SideSet*> mSideSets;   // Mesh owns sidesets
```

**Derived containers (rebuilt from blocks/sidesets):**
```cpp
// Derived views (rebuilt by finalize())
Cell<mesh::Element*> mElements;   // Collected from blocks
Cell<mesh::Facet*>   mFacets;     // Collected from sidesets
```

**Safe pattern:**
```cpp
// Access container reference each time
Cell<mesh::Element*>& elems = tMesh->elements();
for (mesh::Element* elem : elems) {
    // Use element
}

// WRONG: Store container pointer across finalize()
Cell<mesh::Element*>* pElems = &tMesh->elements();
tMesh->finalize();     // Invalidates container
// (*pElems)[0];       // UNDEFINED - container may be reallocated
```

---

## Index Stability Contract

### 4. Indices are NOT Stable Across Operations

> **Rule:**
> Never store indices across `finalize()`, `partition()`, or `Distributor::run()`.
> IDs are stable; indices are not.

**Why:** Indices are reordered by:
- `finalize()` — assigns continuous 0-based indices
- `partition()` — reorders by MPI ownership
- `Distributor::run()` — creates new local ordering
- Graph algorithms (RCM, METIS) — bandwidth/partitioning reordering

**Safe pattern:**
```cpp
// Store IDs (permanent)
Vector<id_t> importantNodeIDs = {100, 200, 300};

tMesh->partition(comm_size());  // Reorders indices
tMesh->finalize();               // Rebuilds indices

// Access by ID (stable)
for (id_t nodeID : importantNodeIDs) {
    mesh::Node* node = tMesh->node(nodeID);  // ID lookup always works
}
```

**Dangerous pattern:**
```cpp
// WRONG: Store indices
Cell<index_t> nodeIndices;
for (mesh::Node* node : tMesh->nodes()) {
    if (node->x() > 0.0) {
        nodeIndices.push(node->index());  // Index at this moment
    }
}

tMesh->partition(comm_size());  // INVALIDATES INDICES

// These indices now point to wrong nodes!
for (index_t idx : nodeIndices) {
    // mesh::Node* node = tMesh->nodes()(idx);  // WRONG NODE!
}
```

---

## Connectivity Invalidation Contract

### 5. Connectivity Invalidation Rules

> **Contract:**
> Cached connectivity becomes **invalid** after certain operations. Access after invalidation is undefined behavior.

| Operation | Invalidates Connectivity |
|-----------|--------------------------|
| `create_edges()` | Node↔Edge, Edge↔Element, Element↔Edge |
| `create_faces()` | Node↔Face, Face↔Element, Element↔Face |
| `scale_mesh()` | nothing in the mesh; only node/control-point coordinates change |
| Node duplication | Node↔Element, Node↔Node, Element↔Node |
| `partition()` | **All** ownership & adjacency |
| `Distributor::run()` | **All** (creates new mesh) |
| `unfinalize()` | Node*/Edge*/Face* connectivities (selective; ElementToNode kept) |

**Safe pattern:**
```cpp
// Compute connectivity
tMesh->finalize();  // Triggers connectivity computation

// Use connectivity
mesh::Node* node = tMesh->node(nodeID);
uint nElems = node->number_of_elements();  // Accesses NodeToElement

// Invalidate connectivity
tMesh->create_edges();  // Adds new entities

// Recompute connectivity
tMesh->finalize_edges();  // Rebuild edge connectivities

// Safe to access again
uint nEdges = node->number_of_edges();
```

**Dangerous pattern:**
```cpp
// Cache connectivity result
mesh::Node* node = tMesh->node(nodeID);
uint nElems = node->number_of_elements();

// Invalidate connectivity
tMesh->create_edges();

// WRONG: Use cached result from before invalidation
// for (uint i = 0; i < nElems; ++i) {
//     mesh::Element* elem = node->element(i);  // UNDEFINED
// }

// CORRECT: Re-query after invalidation
uint nElemsNew = node->number_of_elements();  // Recomputes
```

---

## Type Homogeneity Contract

### 6. Block Homogeneity Assumption

> **Assumption:**
> Most FEM kernels assume blocks contain a **single element type**.
> While BELFEM allows heterogeneous blocks, you must handle element-type dispatch yourself.

**Why this matters:**
```cpp
// Typical FEM code assumes homogeneity
mesh::Block* block = tMesh->block(blockID);
ElementType type = block->element_type();  // Returns FIRST element type

// Allocates shape functions for this type
ShapeFunction* shape = new ShapeFunction(type);

// Loops over all elements (assumes all same type)
for (mesh::Element* elem : block->elements()) {
    // shape->evaluate(...);  // WRONG if elem->type() != type
}
```

**Safe handling:**
```cpp
// Check for heterogeneity
mesh::Block* block = tMesh->block(blockID);
bool isHomogeneous = true;
ElementType firstType = block->element(0)->type();

for (mesh::Element* elem : block->elements()) {
    if (elem->type() != firstType) {
        isHomogeneous = false;
        break;
    }
}

if (!isHomogeneous) {
    // Handle heterogeneous block with element-specific dispatch
    for (mesh::Element* elem : block->elements()) {
        ElementType type = elem->type();
        // Create type-specific shape function for each element
    }
}
```

---

## Facet Lifetime Contract

### 7. Facets are Views, Not Standalone Geometry

> **Invariant:**
> Facets wrap an owned lower-dimensional element and link to the master/slave elements they sit on.

**From `cl_Facet.hpp`:**
```cpp
class Facet : public Vertex {
    Element* mElement;    // owned facet element (its own node list)
    Element* mMaster;     // element this facet sits on
    Element* mSlave;      // neighbor (nullptr on boundary)
    suint mMasterFaceID;  // local facet index on master
    suint mSlaveFaceID;   // local facet index on slave
};
```

**Implications:**
- Deleting a master element invalidates all facets referencing it
- Facet nodes are available directly via `facet->node(k)`; the same nodes can be read from the master via `master->get_nodes_of_facet(facet->index_on_master(), nodes)`
- Facet orientation is relative to master element

**Safe pattern:**
```cpp
mesh::SideSet* sideset = tMesh->sideset(sidesetID);
for (mesh::Facet* facet : sideset->facets()) {
    mesh::Element* master = facet->master();  // check has_master() on slave-only facets
    mesh::Element* slave = facet->slave();    // nullptr on boundary

    // Access facet nodes via master
    uint masterFacetIdx = facet->index_on_master();
    Cell<mesh::Node*> facetNodes;
    master->get_nodes_of_facet(masterFacetIdx, facetNodes);
}
```

**Both routes give the same nodes:**
```cpp
mesh::Facet* facet = sideset->facets()(0);
mesh::Node* node = facet->node(0);          // facet's own element

mesh::Element* master = facet->master();
uint facetIdx = facet->index_on_master();
Cell<mesh::Node*> nodes;
master->get_nodes_of_facet(facetIdx, nodes);
```

---

## Edge/Face Creation Order Contract

### 8. Edges Must Be Created Before Faces (3D Elements)

> **Rule:**
> In 3D, `create_edges()` must be called before `create_faces()`.

**Why:** Face DOFs often reference edge DOFs in Nédélec H(curl) and H(div) elements. Edge containers must exist before face creation.

**Correct order:**
```cpp
tMesh->finalize();        // 1. Finalize main mesh
tMesh->create_edges();    // 2. Create edges first
tMesh->finalize_edges();  // 3. Finalize edge connectivity
tMesh->create_faces();    // 4. Create faces (depends on edges)
tMesh->finalize_faces();  // 5. Finalize face connectivity
```

**Wrong order:**
```cpp
tMesh->finalize();
tMesh->create_faces();    // WRONG: Faces need edges
tMesh->create_edges();    // TOO LATE
```

---

## Parallel Ghost Layer Guarantees

### 9. Ghost Layer Completeness After Distribution

> **`partial_mesh()` is worker-only.** It carries
> `BELFEM_ERROR( mCommRank > 0, ... )` (`cl_Mesh_Distributor.cpp:2663`), so the root rank must
> not call it — the root keeps the mesh it started with and takes the communication tables
> instead (`cl_FEM_Kernel.cpp:688-708`).
>
> **Guarantees after `Distributor::run()` / `partial_mesh()`:**
>
> * All elements have **valid node pointers** (owned or ghost)
> * All ghost elements have **complete connectivity**
> * Nodes are **owned by lowest-rank proc** containing them
> * Sidesets may be **incomplete** unless explicitly redistributed

**Checking ownership:**
```cpp
for (mesh::Element* elem : localMesh->elements()) {
    if (elem->owner() == comm_rank()) {
        // Owned element - compute on this proc
    } else {
        // Ghost element - only for connectivity
    }
}

for (mesh::Node* node : localMesh->nodes()) {
    if (node->owner() == comm_rank()) {
        // Owned node - this proc writes DOF values
    } else {
        // Ghost node - receive values from owner
    }
}
```

---

## Performance-Critical Contracts

### 10. ID Lookup vs. Index Access Performance

> **Performance note:**
> Access via `node(id_t)` uses a hash map (`O(1)` average, `O(N)` worst case).
> Prefer index-based iteration in tight loops.

**Slow (ID-based access in loop):**
```cpp
// O(N) hash lookups
for (id_t nodeID : nodeIDs) {
    mesh::Node* node = tMesh->node(nodeID);  // Map lookup each iteration
    // compute something...
}
```

**Fast (index-based iteration):**
```cpp
// O(1) direct container access
Cell<mesh::Node*>& nodes = tMesh->nodes();
for (index_t i = 0; i < nodes.size(); ++i) {
    mesh::Node* node = nodes(i);  // Direct pointer access
    // compute something...
}
```

---

## Coordinate Modification Invalidation

### 11. What Breaks When Coordinates Change

> **Contract:**
> Modifying node coordinates invalidates cached geometric data:
> * Jacobians (element transformation matrices)
> * Normals (facet orientations)
> * Integration weights (quadrature points)

**Safe pattern:**
```cpp
// Uniform scaling of all coordinates; the mesh keeps no geometry caches,
// downstream kernels must recompute theirs
tMesh->scale_mesh(0.001);

// Manual coordinate change (advanced)
mesh::Node* node = tMesh->node(nodeID);
node->set_coords(newX, newY, newZ);

// Must manually invalidate affected elements
// (typically by re-computing element jacobians before next use)
```

---

## Summary: Critical Rules for Claude Code

**Paste this checklist into Claude Code guidance:**

1. ✓ Treat `finalize()` as establishing mesh invariants; do not access topology before it.
2. ✓ Do not assume non-master MPI ranks have valid topology before distribution.
3. ✓ Treat `unfinalize()` as cache invalidation, not rollback.
4. ✓ Never store indices across `finalize()`, `partition()`, or `Distributor::run()`.
5. ✓ Assume blocks may be heterogeneous unless explicitly checked.
6. ✓ Treat facets as views into elements, never standalone entities.
7. ✓ Assume connectivities can be invalidated by topology changes.
8. ✓ Prefer index-based iteration for performance-critical loops.
9. ✓ Create edges before faces in 3D.
10. ✓ Check `owner()` for all entities in MPI context.

---

**End of Mesh Contracts and Invariants**
