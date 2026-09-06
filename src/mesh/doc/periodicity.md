# Periodic Boundary Conditions {#mesh_periodicity}

**Date:** 2026-03-03
**Module:** src/mesh
**Purpose:** Documentation of periodic boundary condition infrastructure in BELFEM


> **⚠️ API references in this document are stale and are being reviewed.** Several function names
> used below do not exist in the tree (searched `src/**/*.{hpp,cpp}`):
> `create_dofwise_periodicities_master()`, `collect_nodes_from_flags_12()`,
> `flag_periodic_entities_12()` and `match_faces()`. `Periodicity`'s actual surface is
> `master_nodes()` / `slave_nodes()` / `master_edges()` / … plus `update()` and
> `set_entity_dependencies()` (`cl_Mesh_Periodicity.hpp:87-125`), and the DOF-side entry point is
> `DofData::create_dofwise_t_matrices_master()`.
>
> The live path is `Periodicity::update()` + `set_entity_dependencies()`, driven from
> `MaxwellFactory` (`cl_MaxwellFactory.cpp:684-689`, `cl_Mesh_Periodicity.cpp:49-75`) and consumed
> by `create_dofwise_t_matrices_master()`.
>
> **What the described behavior is worth:** the *shape* of the algorithm — how master/slave pairs
> are matched and how the constraints reach the DOFs — is the reason to read this page, and no
> part of it has been shown wrong. But it has not been re-derived against that live path either,
> so treat it as a sketch to check rather than a specification. **Do not copy the snippets.**

---

## Overview

Periodic boundary conditions constrain DOFs on one boundary plane to match their counterparts on an opposite boundary plane. BELFEM implements this via a three-stage pipeline:

1. **PeriodicityFactory** - Detects boundary planes, pairs facets via k-d tree, matches nodes per facet pair
2. **Periodicity** - Stores matched entity pairs (nodes, edges, faces), derives edge/face pairs from node mapping
3. **Distributor** - Preserves periodic pairing across MPI ranks

After thin shells and cuts are created, the periodicity is rebuilt to include all duplicate nodes, and DOF-level constraints are applied from `cl_FEM_DofMgr_DofData.cpp`. The real entry point there is `DofData::create_dofwise_t_matrices_master()` (`cl_FEM_DofMgr_DofData.hpp:266`, defined at `:3478`); the mesh-side rebuild is driven from `MaxwellFactory` (`cl_MaxwellFactory.cpp:684-690`).

---

## Table of Contents

1. [PeriodicityFactory: Node Matching](#periodicityfactory-node-matching)
2. [Periodicity: Entity Pair Storage](#periodicity-entity-pair-storage)
3. [Self-Paired Entities](#self-paired-entities)
4. [MPI Distribution](#mpi-distribution)
5. [Periodicity Propagation Through Node Creation Stages](#periodicity-propagation-through-node-creation-stages)
6. [Key Files](#key-files)

---

## PeriodicityFactory: Node Matching

**File:** `cl_Mesh_PeriodicityFactory.{hpp,cpp}`

### Plane Definition

Each periodic boundary is a plane defined by three `.geo` Point IDs (A, B, C). The factory computes an orthonormal basis `[P; Q; N]` and plane offset `D`:

```
P = normalized(B - A)          in-plane direction 1
N = normalized(cross(P, C-A))  surface normal
Q = cross(N, P)                in-plane direction 2
T = [P; Q; N]                  transformation matrix
D = dot(A, N)                  plane offset
```

The transform `Y = T * X` maps 3D coordinates to in-plane coordinates `(Y(0), Y(1))` and signed distance from the plane `Y(2)`.

### Sideset Detection

`select_sidesets()` iterates all mesh sidesets. A sideset is on the plane if **all** its nodes satisfy `|dot(X, N) - D| < epsilon`. The matching sidesets are used by `map_facets()` to pair facets; they are not stored on the `Periodicity` object.

### K-d Tree Matching

1. The facet centroids of the slave-side sidesets are projected to in-plane coordinates and inserted into a 2D k-d tree (median-based, alternating x/y splits).
2. Each master-side facet centroid is projected with the master transform, then the closest slave centroid is found via nearest-neighbor search (Euclidean distance, pruning). This pairs the facets.
3. Within each facet pair, `match_nodes()` compares the projected node coordinates exhaustively. The result is a 1:1 master-slave node pairing stored in `Periodicity::master_nodes()` and `slave_nodes()`.

### Usage

```cpp
mesh::PeriodicityFactory tFactory( tMesh );
tFactory.set_master_plane( pointA_id, pointB_id, pointC_id );
tFactory.set_slave_plane(  pointD_id, pointE_id, pointF_id );
tFactory.create_periodicity();
// tMesh->periodicity() now holds the matched pairs
```

The three points defining each plane must produce the same in-plane coordinate system (same P and Q directions). This is ensured by choosing corresponding geometry points (e.g., A maps to D, B to E, C to F).

---

## Periodicity: Entity Pair Storage

**File:** `cl_Mesh_Periodicity.{hpp,cpp}`

### Data Layout

The `Periodicity` class stores the master and slave entities in separate arrays:

```
mMasterNodes / mSlaveNodes
mMasterEdges / mSlaveEdges
mMasterFaces / mSlaveFaces
```

These arrays are positionally paired once the factory has run: `match_edges()` compacts the edge containers to tied pairs, and `crosslink()` / `to_proto()` / `from_proto()` rely on that alignment. The bidirectional `periodic()` pointer set by `match_edges()` / `match_facets_and_faces()` carries the same pairing and is what `set_entity_dependencies()` follows.

Node pairs are populated externally (by `PeriodicityFactory` from the deterministic node-pair backup, or `collect_nodes_from_flags_12()`). Edge and face pairs are derived by `PeriodicityFactory::update_periodicity()`.

### The `update()` Pipeline

`Periodicity::update()` (`cl_Mesh_Periodicity.cpp:38-46`) rebuilds a `PeriodicityFactory` from the stored master/slave planes and calls `PeriodicityFactory::update_periodicity( this )`; the edge/face matching (`collect_edges()`, `match_edges()`, `match_facets_and_faces()`, `fix_face_slaves()`, `crosslink()`) lives in the factory, not on `Periodicity`.

### Node Index Space

`update_node_indices()` assigns periodic indices 0..N-1 to both master and slave nodes:

```
master(0) -> index 0    slave(0) -> index 0
master(1) -> index 1    slave(1) -> index 1
...
```

This shared index space allows edge and face matching: an edge with corner periodic indices `(A, B)` on the master side matches the edge with the same indices `(A, B)` on the slave side.

### Edge Matching

`match_edges()` uses a `Map<key_t, Edge*>` with sorted pair keys:

```
key = max(A,B) * N + min(A,B)
```

where A, B are periodic indices of the edge's corner nodes and N is the number of periodic nodes. The indices are taken from `node->original()->index()`, not `node->index()`: a cohomology cut or thin-shell duplicate on the seam keys to its geometric original, so the master and slave sides resolve to the same shared index. The duplicate's own periodicity is then induced through its sources (hanging chain), so it never needs to be a matched periodic node itself.

**Pass 1:** Flag master nodes, collect master edges (both corners flagged), insert into map.
**Pass 2:** Flag slave nodes, find slave edges (both corners flagged), look up master via same key.

Slave edge nodes are reordered to match master edge orientation using `insert_node()`.

### Face Matching

`match_faces()` uses a `Map<key128_t, Face*>` with sorted triple keys:

```
key = (N * sorted(2) + sorted(1)) * N + sorted(0)
```

where sorted(0..2) are the three smallest periodic indices of the face's corner nodes.

### Slave Face Orientation

`fix_face_slaves()` reassigns each slave face from master to slave role with the correct orientation. For each slave face:

1. Find which facet of its volume element matches the face's nodes
2. Determine the rotation index by comparing the first corner node of the element facet against the face's nodes (using `original()->id()` for thin-shell duplicate nodes)
3. Call `set_slave(element, facet_index, orientation)` to reassign

### Bidirectional Pointers

`set_periodic_entities()` sets `periodic()` pointers in both directions:

```cpp
tMaster->set_periodic( tSlave );
tSlave->set_periodic( tMaster );
```

This applies to nodes, edges, and faces. Self-paired entities (on the rotation axis) point to themselves.

---

## Self-Paired Entities

Nodes that lie on both periodic planes (e.g., on the rotation axis of a cylindrical geometry) are paired with themselves: `node->periodic() == node`. This is necessary for edge and face matching to work correctly, since matching keys are built from periodic indices which require every boundary node to have a valid index.

Self-paired entities are skipped in `create_dofwise_periodicities_master()` (since constraining a DOF to itself is a no-op).

---

## MPI Distribution

**File:** `cl_Mesh_Distributor.cpp`

### Pairing Preservation

The distributor does not track periodic roles with flags and does not transmit a `Periodicity` object. Periodic slaves are hanging entities whose sources are their masters (`Periodicity::set_entity_dependencies()`), so `select_sources()` ghosts the master with the slave (`cl_Mesh_Distributor.cpp:326-330`), and the constraints themselves travel as the T-matrices (`send_t_matrices()` / `receive_t_matrices()`). The partner table is serialized only on the `.bfm` path, by `PeriodicityFactory::to_proto()` / `from_proto()` via `ProtoMesh::create_periodicitiy()`.

### Face Reconstruction

On receiving procs, `ProtoMesh::create_faces()` reconstructs faces from transmitted data. Slave periodic faces (which have no master element, only a slave) are handled by the Face constructor's slave-only branch, using `to_master_orientation()` with the stored `mOrientationOnSlave` to reorder nodes correctly.

---

## Periodicity Propagation Through Node Creation Stages

After the initial `PeriodicityFactory::create_periodicity()` call, several operations create new nodes that must also be paired periodically:

1. **CutFactory::duplicate_nodes_on_face_sidesets()** - Thin-shell duplicate nodes
2. **CutFactory::link_node_duplicates_and_originals()** - Cut duplicate nodes
3. **ThinShellFactory::create_nodes_on_layers()** - Extruded layer nodes

Each stage propagates the `periodic()` pointer and flags 1/2 from the original node to the new node's periodic counterpart.

### Flag Propagation Design

The key insight is that `unflag_all_nodes()` only clears flag 0, so flags 1 and 2 survive through mesh operations. Each node creation stage checks `is_periodic()` on the source node and, if true:

1. Finds the periodic counterpart's corresponding duplicate
2. Sets bidirectional `periodic()` pointers between the new nodes
3. Copies flags 1/2 from the original pair to the new pair

### Rebuilding After All Stages

In `MaxwellFactory::create_magnetic_kernel()`, after all node creation stages are complete:

```cpp
// 1. Set flags 1/2 on original periodic entities
tPeriodicity->flag_periodic_entities_12();

// 2. CutFactory and ThinShellFactory propagate flags to new nodes

// 3. Rebuild node lists from all flagged nodes
tPeriodicity->collect_nodes_from_flags_12();

// 4. Derive edge and face pairs from the expanded node lists
tPeriodicity->update();
```

`collect_nodes_from_flags_12()` walks all mesh nodes, collects those with flag 1 into `mMasterNodes` and their `periodic()` partners into `mSlaveNodes`. This captures the original boundary nodes plus all thin-shell duplicates, cut duplicates, and layer nodes.

---

## Key Files

| File | Purpose |
|------|---------|
| `cl_Mesh_PeriodicityFactory.{hpp,cpp}` | Plane detection, k-d tree matching, node pairing |
| `cl_Mesh_Periodicity.{hpp,cpp}` | Entity pair storage, edge/face derivation, update pipeline |
| `cl_Mesh_Distributor.cpp` | MPI distribution with periodic pairing preservation |
| `cl_ProtoMesh.cpp` | Face reconstruction on receiving procs |
| `fn_to_master_orientation.{hpp,cpp}` | Slave-to-master node reordering for faces |
| `src/homology/cl_CutFactory.cpp` | Periodicity propagation through thin-shell and cut duplicates |
| `src/fem/kernel/cl_ThinShellFactory.cpp` (namespace `mesh`) | Periodicity propagation through layer nodes |
| `src/fem/maxwell/cl_MaxwellFactory.cpp` | Orchestrates rebuild via `collect_nodes_from_flags_12()` + `update()` |
| `src/fem/kernel/cl_FEM_DofMgr_DofData.cpp` | DOF-level periodic constraints (`create_dofwise_periodicities_master()`) |
