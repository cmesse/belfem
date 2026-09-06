# BELFEM Mesh Module Tests — Detailed Plan

**Date:** 2026-03-22
**Purpose:** Method-level test matrix for the mesh module (`src/mesh/`)
**Depends on:** `tests_0_strategy.md`, `tests_1_containers.md` (Cell, Map), `tests_7_graph.md` (Graph/Vertex)
**Confidence:** High on entity contracts and element catalog. Medium on connectivity workflows (complex interactions). Low on I/O, distribution, periodicity (deferred).

---

## Module Overview

The mesh module is BELFEM's geometric backbone (~49K lines, ~100 source files). For this test plan, we focus on the testable core and defer I/O, parallel distribution, periodicity, and thin shells to future phases.

| Layer | In Scope | Deferred                                             |
|---|---|------------------------------------------------------|
| Entity hierarchy | `Node`, `Element`, `ElementFactory`, `Basis`, `Vertex` | `ControlPoint`, `Segment`                            |
| Element catalog | All 30+ element types via `ElementFactory` | —                                                    |
| Mesh class | Tensor mesh constructor, accessors, maps | File-based constructors                              |
| Connectivity | `ConnectivityCalculator` on small meshes | —                                                    |
| Containers | `Block`, `SideSet` basics | `Curve`, `ThinShell`                                 |
| Processing | `scale_mesh`, `finalize` | `OrderConverter`, `Distributor`, `Partitioner`       |
| I/O | — | `GmshReader`, `BfmFile`, `ExodusWriter`, `VtkWriter` |

---

## Testing Entry Points

### 1. ElementFactory (Data-Driven Catalog Test)

`ElementFactory::create_element(ElementType, id)` instantiates any element type as `ElementTemplate<N,C,E,T,F>`. The template parameters encode the topology. A single data-driven test can verify the entire catalog.

### 2. Tensor Mesh Constructor (Integration Gateway)

`Mesh(order, numNodes, step, origin)` creates structured grids programmatically — no file I/O needed. Supported:

| Order | 2D Element | 3D Element |
|---|---|---|
| 2 | QUAD4 | HEX8 |
| 3 | QUAD9 | HEX27 |
| 4 | QUAD16 | HEX64 |

This constructor returns a complete mesh with nodes, elements, blocks, and correct coordinates.

---

## Test File Structure

```
tests/mesh/
├── test_MeshEntities.cpp         # Node, Element, factory, topology catalog
├── test_Mesh.cpp                 # Tensor mesh, accessors, maps, connectivity
```

---

## 1. Element Factory & Topology Catalog

**File:** `test_MeshEntities.cpp`

### 1.1 Element Catalog (Data-Driven) `[semantic]`

A single parameterized test verifies every element type. The expected topology for each type:

| ElementType | N | C | E | T | F | Dim |
|---|---|---|---|---|---|---|
| VERTEX | 1 | 1 | 0 | 0 | 0 | 0 |
| LINE2 | 2 | 2 | 1 | 0 | 0 | 1 |
| LINE3 | 3 | 2 | 1 | 0 | 0 | 1 |
| LINE4 | 4 | 2 | 1 | 0 | 0 | 1 |
| LINE5 | 5 | 2 | 1 | 0 | 0 | 1 |
| TRI3 | 3 | 3 | 3 | 3 | 1 | 2 |
| TRI6 | 6 | 3 | 3 | 3 | 1 | 2 |
| TRI10 | 10 | 3 | 3 | 3 | 1 | 2 |
| TRI15 | 15 | 3 | 3 | 3 | 1 | 2 |
| QUAD4 | 4 | 4 | 4 | 4 | 1 | 2 |
| QUAD8 | 8 | 4 | 4 | 4 | 1 | 2 |
| QUAD9 | 9 | 4 | 4 | 4 | 1 | 2 |
| QUAD16 | 16 | 4 | 4 | 4 | 1 | 2 |
| TET4 | 4 | 4 | 6 | 4 | 4 | 3 |
| TET10 | 10 | 4 | 6 | 4 | 4 | 3 |
| PENTA6 | 6 | 6 | 9 | 5 | 5 | 3 |
| PENTA15 | 15 | 6 | 9 | 5 | 5 | 3 |
| PENTA18 | 18 | 6 | 9 | 5 | 5 | 3 |
| PYRA5 | 5 | 5 | 8 | 5 | 5 | 3 |
| PYRA13 | 13 | 5 | 8 | 5 | 5 | 3 |
| PYRA14 | 14 | 5 | 8 | 5 | 5 | 3 |
| HEX8 | 8 | 8 | 12 | 6 | 6 | 3 |
| HEX20 | 20 | 8 | 12 | 6 | 6 | 3 |
| HEX27 | 27 | 8 | 12 | 6 | 6 | 3 |
| HEX64 | 64 | 8 | 12 | 6 | 6 | 3 |

For each row: create via `ElementFactory::create_element(type, 1)`, verify `number_of_nodes()`, `number_of_corner_nodes()`, `number_of_edges()`, `number_of_facets()`, `number_of_faces()`, `dimension()`, `type()`.

**VERIFIED 2026-03-22:** All 25 rows confirmed correct against source (`cl_Element_Factory.cpp`, `cl_ElementTemplate.hpp`).

**VERIFIED 2026-03-30:** `TET20` and `TET35` now have dedicated mesh headers providing `type()`, `dimension()`, facet extraction, and edge extraction. Keep them in the normal catalog tests.

**Also not in catalog:** Thin-shell variants (QUAD4TS, QUAD9TS, PENTA6TS, PENTA18TS), EMPTY, TRI21, LINE6, UNDEFINED — either specialized or incomplete.

**BUG-ME1 (VERIFIED 2026-03-23):** `Mesh::node(i,j,k)` at `cl_Mesh.hpp:1205` asserts `this->number_of_dimensions() == 2` instead of `== 3`. The 3-argument accessor is broken for 3D tensor meshes. Add a regression test.

**Implementation:** Use `TEST_P` with `testing::Values(...)` or a simple loop.

### 1.2 Element Factory Unknown Type `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `FactoryUnknownTypeThrows` | `create_element(ElementType::UNDEFINED, 1)` → `BELFEM_ERROR` fires |

---

## 2. Node Entity

**File:** `test_MeshEntities.cpp`

### 2.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `NodeConstruction` | `Node(1, 2.0, 3.0, 4.0)` → `id()==1`, `x()==2.0`, `y()==3.0`, `z()==4.0` |
| `NodeDefaultCoords` | `Node(1)` → `x()==0`, `y()==0`, `z()==0` |
| `NodeSetCoords3D` | `set_coords(1.0, 2.0, 3.0)` → coordinates updated |
| `NodeSetCoordsVector` | `set_coords(Vector{5.0, 6.0, 7.0})` → coordinates updated |
| `NodeEntityType` | `entity_type() == EntityType::NODE` |
| `NodeCoordsVector` | `coords()` returns length-3 vector matching x(), y(), z() |

---

## 3. Element Entity Contracts

**File:** `test_MeshEntities.cpp`

### 3.1 Node Insertion and Access `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `ElementInsertAndAccessNodes` | Create TRI3 via factory, insert 3 nodes, verify `node(0..2)` returns correct pointers |
| `ElementFlagNodes` | Insert 3 nodes into TRI3, call `flag_nodes()` → all nodes `is_flagged()` |
| `ElementUnflagNodes` | After flagging, `unflag_nodes()` → all nodes `!is_flagged()` |
| `ElementIdStored` | `Element(42)` → `id() == 42` |

### 3.2 Facet Topology (Representative Types) `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `TRI3FacetNodes` | For each facet 0–2: `get_nodes_of_facet` returns correct pair of nodes |
| `QUAD4FacetNodes` | For each facet 0–3: correct node pairs |
| `TET4FacetNodes` | For each facet 0–3: correct triangle of nodes |
| `HEX8FacetNodes` | For each facet 0–5: correct quad of nodes |

### 3.3 Element Base Class Guard `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `BaseElementNodeThrows` | Calling `node()` on base `Element` (not via factory) → `BELFEM_ERROR` fires (always active) |
| `BaseElementTypeThrows` | `type()` on base → `BELFEM_ERROR` |

---

## 4. Tensor Mesh Construction

**File:** `test_Mesh.cpp`

### 4.1 2D Mesh `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `TensorMesh2DQuad4` | `Mesh(2, {3,3}, {1.0,1.0})` → 9 nodes, 4 QUAD4 elements, 1 block |
| `TensorMesh2DQuad4NodeCoords` | Corner nodes at (0,0), (2,0), (0,2), (2,2) |
| `TensorMesh2DQuad4Dimensions` | `number_of_dimensions() == 2` |
| `TensorMesh2DQuad9` | `Mesh(3, {5,5}, {1.0,1.0})` → 25 nodes, 4 QUAD9 elements |
| `TensorMesh2DNodeCount` | For order=2: `(Nx)*(Ny)` nodes. For order=3: node count matches grid formula |

### 4.2 3D Mesh `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `TensorMesh3DHex8` | `Mesh(2, {3,3,3}, {1.0,1.0,1.0})` → 27 nodes, 8 HEX8 elements |
| `TensorMesh3DHex27` | `Mesh(3, {5,5,5}, {1.0,1.0,1.0})` → 125 nodes, 8 HEX27 elements |
| `TensorMesh3DNodeCoords` | Corner nodes at expected positions |

### 4.3 Mesh with Origin Offset `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `TensorMesh2DWithOrigin` | `Mesh(2, {3,3}, {1.0,1.0}, {5.0, 10.0})` → node at (5,10) exists |

### 4.4 Mesh Accessors `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `MeshNodeCount` | `number_of_nodes()` matches expected |
| `MeshElementCount` | `number_of_elements()` matches expected |
| `MeshBlockCount` | `number_of_blocks()` ≥ 1 |
| `MeshIsTensormesh` | `is_tensormesh() == true` |
| `MeshNodeByID` | `node(id)` returns correct node |
| `MeshElementByID` | `element(id)` returns correct element |
| `MeshNodeByIJ` | `node(i,j)` returns node at grid position (tensor mesh only) |

---

## 5. Map and Index Consistency

**File:** `test_Mesh.cpp`

### 5.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `NodeMapConsistent` | For every node in `nodes()`: `node(node->id()) == node` |
| `ElementMapConsistent` | For every element: `element(element->id()) == element` |
| `BlockMapConsistent` | For every block: `block(block->id()) == block` |
| `NodeIndicesContiguous` | After finalize: node indices are 0..N-1 |
| `ElementIndicesContiguous` | After finalize: element indices are 0..N-1 |

---

## 6. Connectivity on Tensor Mesh

**File:** `test_Mesh.cpp`

### 6.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `NodeToElementConnectivity` | After finalize: each interior node connected to 4 elements (2D QUAD4 grid) |
| `CornerNodeConnectivity` | Corner node of 2D grid connected to exactly 1 element |
| `EdgeNodeConnectivity` | Edge node (not corner) connected to exactly 2 elements |
| `ElementToElementConnectivity` | Interior element has 4 neighbors (QUAD4 2D), boundary element has fewer |
| `EdgesExistAfterFinalize` | `edges_exist() == true` after finalize |
| `EdgeCountCorrect` | Number of edges matches expected for structured grid |

---

## 7. Scale Mesh

**File:** `test_Mesh.cpp`

### 7.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `ScaleMeshFactor` | `scale_mesh(0.001)` → all node coordinates multiplied by 0.001 |
| `ScaleMeshPreservesTopology` | After scaling: same number of nodes, elements; connectivity unchanged |

---

## 8. Block and SideSet Basics

**File:** `test_Mesh.cpp`

### 8.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `BlockElementType` | Block in tensor mesh reports correct element type |
| `BlockElementCount` | Block element count matches mesh element count (single-block mesh) |
| `BlockElementAccess` | `block->element(0)` returns valid element pointer |

---

## 9. Memory and Cleanup

### 9.1 Tests `[semantic]` `[valgrind]`

| Test Name | What It Verifies |
|---|---|
| `MeshDestructorCleansUp` | Create tensor mesh in scope, let it destruct → no leaks (Valgrind target) |
| `ElementFactoryOwnership` | Factory creates element, caller deletes → no double-free |

---

## 10. Implementation Notes for Claude Code

1. **Use `ElementFactory` for all element creation.** Never instantiate `ElementTemplate` directly in tests.
2. **Use tensor mesh constructor** for all mesh-level tests. It's the `SpMatrix(dense)` equivalent — a complete, valid mesh without file I/O.
3. **Memory ownership:** `Mesh` owns and deletes all its entities. Elements created via factory for standalone tests must be deleted by the test.
4. **The element catalog test is the single highest-value test.** It locks down the topology constants for all 25+ element types in one sweep.
5. **Connectivity tests require `finalize()`.** The tensor mesh constructor may or may not call finalize internally — verify before writing connectivity tests.
6. **Node IDs in tensor meshes are 1-based** (typical for mesh formats). Don't assume 0-based.
7. **Facet topology tests** should create elements via factory, insert real `Node` objects, then call `get_nodes_of_facet` and verify the correct nodes are returned. Remember to clean up.
8. **2D meshes have z=0 for all nodes.** Verify this explicitly.
9. **The `Element` base class virtual methods fire `BELFEM_ERROR` (always active)**, not `BELFEM_ASSERT`. This means these guards work in release builds too.
10. **Deferred items** (I/O, distribution, periodicity, thin shells, curves, order conversion) should be noted as "future" in the test file header, not silently omitted.

---

## 11. Codex Audit Checklist

When reviewing Claude Code's test implementation, verify:

- [ ] Element catalog test covers at least 20 element types with correct (N,C,E,T,F,dim) values
- [ ] Factory is used for all element creation, never direct `ElementTemplate` instantiation
- [ ] Tensor mesh tests verify node count, element count, and representative coordinates
- [ ] Map consistency tests verify node/element/block maps agree with containers
- [ ] Connectivity tests use a small tensor mesh and verify node-to-element counts
- [ ] Facet topology tested for at least TRI3, QUAD4, TET4, HEX8
- [ ] All standalone elements created with factory are deleted by the test
- [ ] Mesh destructor cleanup is a Valgrind target
- [ ] BELFEM naming conventions (`t` prefix for locals)
- [ ] Deferred items documented in test file header
