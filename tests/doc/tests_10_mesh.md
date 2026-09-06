# Mesh Module Test Suite

**Date:** 2026-03-23

---

## Test Count by Suite

| Suite | Tests | Coverage |
|-------|-------|----------|
| `ElementCatalog` | 2 | 29-type data-driven catalog (N, C, E, T, F, dim), topology values |
| `NodeEntity` | 6 | Construction, default coords, set_coords scalar, entity type, coords vector, set_coords vector |
| `ElementEntity` | 4 | Node insertion/access, flags, type, dimension |
| `FacetTopology` | 4 | TRI3, QUAD4, TET4, HEX8 facet node ordering (Exodus convention) |
| `ElementBaseGuard` | 2 | Base class `BELFEM_ERROR` guards (always active, NOT debug-only) |
| `TensorMesh2D` | 6 | Construction, node/element count, coords, grid access, element type |
| `TensorMesh3D` | 3 | Construction, node count, 3D grid node access |
| `MeshAccessors` | 3 | Node/element by 1-based ID |
| `MeshMaps` | 3 | Node/element map consistency, block map |
| `MeshConnectivity` | 4 | Node-to-element counts, edges NOT auto-created |
| `MeshEdges` | 3 | create_edges(), edge count, edges_exist() |
| `MeshScale` | 2 | Scale factor, coordinate verification |
| `MeshBlock` | 3 | Block count, type, element access |
| `MeshMemory` | 1 | Destructor cleanup |
| **Total** | **46** | |

---

## Key Design Decisions

### Tensor Mesh Order → Element Type Mapping

@warning The test plan's table was WRONG. Corrected from source:

| Order | 2D Element | 3D Element |
|-------|-----------|-----------|
| 1     | QUAD4     | HEX8      |
| 2     | QUAD9     | HEX27     |
| 3     | QUAD16    | HEX64     |

### Connectivity Setup Ritual

Fresh tensor mesh does NOT have connectivity. To test node-to-element counts:

```cpp
tMesh.unfinalize();
tMesh.set_connectivity( belfem::Connectivity::Compute );
tMesh.finalize();
```

### Edges Are Lazy

`finalize()` does NOT create edges. `edges_exist()` is false after finalize. Call `create_edges()` explicitly to materialize them.

### 1-Based IDs

Node and element IDs start at 1, not 0. `mesh.node(0)` crashes. `mesh.node(1)` returns the first node.

### Compilation

Mesh tests require C++17 (`-std=gnu++17`) due to `std::clamp` in `TensorMeshConfig`. Also requires many libraries: mesh, interpolation, integration, spline, sparse, graph, io, containers, comm, core, plus external libs (petsc, armadillo, lapack, blas, umfpack, suitesparse, hdf5, exodus, netcdf, metis, scotch, gfortran, gomp).

---

## Element Catalog (29 Types Tested)

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
| QUAD4TS | 4 | 4 | 2 | **2** | 1 | 2 |
| QUAD9TS | 9 | 4 | 3 | 2 | 1 | 2 |
| PENTA6TS | 6 | 6 | 6 | 2 | 2 | 3 |
| PENTA18TS | 18 | 6 | 9 | 3 | 3 | 3 |

@note TET20 and TET35 are excluded — they lack `type()`/`dimension()` specializations.

---

## @todo Future Work

- [ ] I/O tests: GmshReader, HDF5Reader/Writer, ExodusWriter, VtkWriter
- [ ] Parallel: Distributor, Partitioner, periodicity, ghost workflows
- [ ] Advanced: ThinShell, Curve, OrderConverter
- [ ] Node/element index contiguity checks
- [ ] Element-to-element connectivity tests
- [ ] Reference mesh file (.msh) integration tests
