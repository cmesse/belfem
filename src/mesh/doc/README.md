# Mesh Module Documentation {#mesh_index}

**Module:** src/mesh
**Purpose:** Index of documentation for BELFEM's mesh data structures and I/O module

---

## Overview

The `mesh` module provides data structures and I/O capabilities for finite element meshes, supporting:

- Multiple element types (1D/2D/3D, linear to 5th order)
- Parallel mesh partitioning (METIS) and MPI distribution
- Multiple file formats (Gmsh, HDF5, Exodus, VTK)
- Higher-order topological entities (edges, faces for Nédélec elements)
- Tensor product meshes for structured grids
- Thin-shell elements for layered structures
- Node duplication for cohomology cuts
- Field and global variable management
- Mesh connectivity computation and caching

---

## Documentation Files

### User Guides

- **[mesh_usage_guide.md](mesh_usage_guide.md)** - Comprehensive usage guide for the mesh module
  - Common pitfalls (MUST READ FIRST!)
  - Core classes (Mesh, Node, Element, Block, SideSet)
  - Mesh construction and I/O (Gmsh, HDF5, Exodus, VTK)
  - Element types and topology
  - Mesh connectivity
  - Mesh partitioning (MPI parallelization)
  - Fields and data management
  - Advanced features (edges, faces, thin shells, tensor meshes)
  - Common usage patterns
  - Performance considerations
  - Debugging and visualization

- **[thin_shell_geometry_and_periodicity.md](thin_shell_geometry_and_periodicity.md)** - How `ThinShellFactory` (namespace `mesh`, in `src/fem/kernel/cl_ThinShellFactory.cpp`) builds the thin-shell virtual domain (protoshell → surface S → extruded layers L) and propagates periodicity from the seam nodes P to the extruded layer nodes Q. Key subtlety: Q inherits P's node periodicity by *link* (flags 4/5), not geometric position — Q need not lie on the periodic plane — and the layer wrapper facets get no edges — `PeriodicityFactory::collect_edges()` reconstructs their edge set from the master volume element instead.
- **[periodicity.md](periodicity.md)** - Periodic boundary condition infrastructure
  - PeriodicityFactory (plane detection, k-d tree facet matching, per-facet node matching)
  - Periodicity class (entity pair storage, edge/face derivation)
  - Self-paired entities (rotation axis)
  - MPI distribution with pairing preservation
  - Periodicity propagation through CutFactory and ThinShellFactory

- **[bfm_file_format.md](bfm_file_format.md)** - **The `.bfm` (HDF5) enriched-mesh file format**
  - Design philosophy (enrichment cache; input file stays source of truth)
  - Encoding conventions (ID-keyed references, vlen datasets, transposed coordinate storage)
  - Group-by-group dataset reference (`/meta` … `/curves`)
  - Positional membership rules (element→block, facet→sideset)
  - The reload contract: what is recomputed vs restored vs synthesized
  - Verified by end-to-end restart (serial + 2/4/8 MPI procs, 2026-07-01)

- **[mesh_contracts_and_invariants.md](mesh_contracts_and_invariants.md)** - **Critical contracts for safe usage**
  - MPI finalization semantics
  - `unfinalize()` behavior (NOT inverse of `finalize()`)
  - Ownership vs derived containers
  - Index stability rules
  - Connectivity invalidation rules
  - Block homogeneity assumptions
  - Facet lifetime contracts
  - Edge/face creation order
  - Parallel ghost layer guarantees
  - Performance-critical contracts

---

## Quick Reference

### Main Classes

| Class | File | Purpose |
|-------|------|---------|
| **`Mesh`** | cl_Mesh.{hpp,cpp} | **Top-level container** for all mesh entities |
| `Node` | cl_Node.{hpp,cpp} | Spatial point with 3D coordinates |
| `Element` | cl_Element.{hpp,cpp} | Base class for all element types |
| `Block` | cl_Block.{hpp,cpp} | Collection of elements (subdomain) |
| `SideSet` | cl_SideSet.{hpp,cpp} | Collection of boundary facets |
| `Facet` | cl_Facet.{hpp,cpp} | Reference to element face/edge |
| `Edge` | cl_Edge.{hpp,cpp} | Topological edge (for Nédélec elements) |
| `Face` | cl_Face.{hpp,cpp} | Topological face (for Nédélec elements) |
| `Field` | cl_Mesh_Field.{hpp,cpp} | Data field on mesh entities |
| `GlobalVariable` | cl_Mesh_GlobalVariable.{hpp,cpp} | Scalar mesh-level data |

### I/O Classes

| Class | File | Purpose | Read | Write |
|-------|------|---------|------|-------|
| `GmshReader` | cl_Mesh_GmshReader.{hpp,cpp} | Gmsh format (.msh) | ✓ | |
| `BfmFile` | cl_Mesh_BfmFile.{hpp,cpp} | BELFEM mesh format (.bfm, HDF5-based) — see [bfm_file_format.md](bfm_file_format.md) | ✓ | ✓ |
| `ExodusWriter` | cl_Mesh_ExodusWriter.{hpp,cpp} | Exodus II format (.exo) | | ✓ |
| `VtkWriter` | cl_Mesh_VtkWriter.{hpp,cpp} | VTK format (.vtk) | | ✓ |

### Utility Classes

| Class | File | Purpose |
|-------|------|---------|
| `Partitioner` | cl_Mesh_Partitioner.{hpp,cpp} | Graph-based mesh partitioning (METIS) |
| `Distributor` | cl_Mesh_Distributor.{hpp,cpp} | MPI mesh distribution with ghost layers |
| `OrderConverter` | cl_Mesh_OrderConverter.{hpp,cpp} | Upgrades an all-linear mesh to second order (LINE3, TRI6, QUAD9, TET10, PENTA18, HEX27) |
| `ConnectivityCalculator` | cl_Mesh_ConnectivityCalculator.{hpp,cpp} | Compute entity-to-entity connectivities |
| `ThinShell` | cl_ThinShell.{hpp,cpp} | Thin-shell configuration for layered structures |
| `Periodicity` | cl_Mesh_Periodicity.{hpp,cpp} | Matched periodic entity pairs (nodes, edges, faces) |
| `PeriodicityFactory` | cl_Mesh_PeriodicityFactory.{hpp,cpp} | Plane detection, k-d tree facet matching, per-facet node matching |

---

## Typical Usage Pattern

```cpp
#include "cl_Mesh.hpp"

// 1. Load mesh from file (the reader finalizes the mesh itself)
Mesh* tMesh = new Mesh("mesh.msh");

// 2. Finalize: only needed after building or modifying a mesh by hand;
//    a second call on a file-loaded mesh is a no-op
tMesh->finalize();

// 3. Access mesh entities
index_t nNodes = tMesh->number_of_nodes();
index_t nElements = tMesh->number_of_elements();

mesh::Node* node = tMesh->node(nodeID);         // By ID
mesh::Element* elem = tMesh->element(elemID);   // By ID

Cell<mesh::Block*>& blocks = tMesh->blocks();
Cell<mesh::SideSet*>& sidesets = tMesh->sidesets();

// 4. Create fields
Vector<real>& temperature = tMesh->create_field(
    "Temperature",
    EntityType::NODE);

// 5. Save results
tMesh->save("output.hdf5");

// 6. Cleanup
delete tMesh;
```

---

## Supported Element Types

### 1D Elements
- **LINE2**, **LINE3**, **LINE4**, **LINE5** (linear to quartic)

### 2D Elements
- **Triangles:** TRI3, TRI6, TRI10, TRI15 (linear to quartic). `TRI21` exists in the enum but `ElementFactory::create_element` has no case for it, so a quintic triangle cannot be built
- **Quadrilaterals:** QUAD4, QUAD8, QUAD9, QUAD16 (bilinear to bicubic)

### 3D Elements
- **Tetrahedra:** TET4, TET10, TET20, TET35 (linear to quartic)
- **Hexahedra:** HEX8, HEX20, HEX27, HEX64 (trilinear to tricubic)
- **Prisms:** PENTA6, PENTA15, PENTA18 (linear to quadratic)
- **Pyramids:** PYRA5, PYRA13, PYRA14 (linear to quadratic)

### Special Elements
- **Thin-Shell:** QUAD4TS, QUAD9TS, PENTA6TS, PENTA18TS (cohomology-aware)

---

## File Format Support

| Format | Extension | Read | Write | Parallel I/O | Use Case |
|--------|-----------|------|-------|--------------|----------|
| **Gmsh** | .msh | ✓ | | Serial | Mesh generation (industry standard) |
| **HDF5** | .hdf5, .bfm | ✓ | ✓ | Serial (rank 0) | Production runs, restart files |
| **Exodus II** | .exo | | ✓ | Serial | Time series visualization (ParaView) |
| **VTK** | .vtk | | ✓ | Serial | Quick visualization checks |

---

## Mesh Connectivity Types

**Common connectivities** (36 total types):

| Connectivity | Description |
|--------------|-------------|
| **NodeToElement** | Which elements contain each node |
| **ElementToNode** | Which nodes belong to each element |
| **ElementToElement** | Element neighbors via shared entities |
| **NodeToNode** | Node neighbors via shared elements |
| **EdgeToElement** | Which elements contain each edge |
| **FaceToElement** | Which elements contain each face |
| **FacetToElement** | Which element owns each facet |

See `Connectivity` enum in `Mesh_Enums.hpp` for complete list.

---

## Parallel Mesh Workflow

```cpp
#ifdef BELFEM_MPI
#include "cl_Mesh_Distributor.hpp"

// 1. Master proc loads mesh
Mesh* globalMesh = nullptr;
if (comm_rank() == 0) {
    globalMesh = new Mesh("mesh.msh");
    globalMesh->partition(comm_size());  // METIS partitioning
}

// 2. Distribute to all procs
mesh::Distributor distributor(globalMesh);
distributor.run();                       // every rank
Mesh* localMesh = globalMesh;            // root keeps the full mesh
if (comm_rank() != 0) {
    localMesh = distributor.partial_mesh();   // workers only -- root aborts here
}

// 3. Finalize. run() has already finalized the worker meshes
//    ( cl_Mesh_Distributor.cpp:248 ), so this is the root's own.
if (comm_rank() == 0) {
    localMesh->finalize();
}

// 4. Work with local partition
// Each proc has elements with owner() == comm_rank()
// Plus ghost elements for communication

// 5. Cleanup
if (comm_rank() == 0) {
    delete globalMesh;
}
delete localMesh;
#endif
```

---

## Common Pitfalls Summary

**ALWAYS** finalize a mesh built or modified by hand (the file readers finalize on load):
```cpp
tMesh->finalize();  // REQUIRED after hand assembly or topology changes
```

**NEVER** manually delete mesh entities:
```cpp
delete tMesh->node(0);  // WRONG - causes double-free
```

**DISTINGUISH** between ID (sparse, from file) and index (continuous, 0-based):
```cpp
mesh::Node* node = tMesh->node(12345);      // By ID (uses map)
mesh::Node* node = tMesh->nodes()(42);      // By index (direct access)
```

**CREATE** edges/faces BEFORE FEM kernel initialization:
```cpp
tMesh->create_edges();  // BEFORE create_fields()
```

**USE** correct entity types for fields:
```cpp
// Node field (most common)
tMesh->create_field("Temperature", EntityType::NODE);

// Edge field (for Nédélec H(curl))
tMesh->create_edges();  // Create edges first!
tMesh->create_field("EdgeFlux", EntityType::EDGE);
```

See **mesh_usage_guide.md** for detailed pitfall explanations.

---

## Advanced Features

### Nédélec Elements (Edge/Face DOFs)
```cpp
tMesh->create_edges();      // For H(curl) elements
tMesh->finalize_edges();

tMesh->create_faces();      // For H(div) elements
tMesh->finalize_faces();
```

### Tensor Product Meshes
```cpp
Mesh* tMesh = new Mesh(
    2,                      // order (quadratic)
    {11, 11},               // 11x11 nodes
    {0.1, 0.1},             // grid spacing
    {0.0, 0.0});            // origin

mesh::Node* node = tMesh->node(5, 5);  // Structured access
```

### Thin-Shell Elements
```cpp
mesh::Block* block = tMesh->block(shellBlockID);
block->set_domain_type(DomainType::ThinShell);
block->set_thickness(0.001);  // 1 mm
```

### Node Duplication (Cohomology)
```cpp
mesh::Node* original = tMesh->node(nodeID);
original->allocate_duplicate_container(2);

// Create duplicates for cohomology cuts
// (typically handled by Homology module)
```

---

## Development Notes

### Adding New Element Types

1. Define element type in `Mesh_Enums.hpp`
2. Create element class `cl_Element_<TYPE>.hpp` inheriting from `ElementTemplate`
3. Implement virtual methods (number_of_nodes, node ordering, etc.)
4. Add to `Element_Factory` switch statement
5. Update `to_string()` helper in `Mesh_Enums.hpp`
6. Test with Gmsh-generated meshes

### Adding New I/O Format

1. Create reader/writer classes: `cl_Mesh_<FORMAT>Reader.{hpp,cpp}`
2. Implement mesh entity parsing (nodes, elements, blocks, sidesets)
3. Handle element type mapping to BELFEM types
4. Add to `Mesh::save()` or constructor (format auto-detection via file extension)
5. Document limitations (e.g., Exodus is write-only)

### Performance Profiling

Use `cl_Profiler` and `cl_Timer`:
```cpp
Timer tTimer;

tMesh->partition(comm_size());
message(InfoLevel::Verbose, "Partitioning: %u ms", (uint)tTimer.stop());

tMesh->finalize();
message(InfoLevel::Verbose, "Finalization: %u ms", (uint)tTimer.stop());
```

---

## Source Code

**Module location:** `../../`

**Key source files:**
- **Main container:** `cl_Mesh.{hpp,cpp}`
- **Entities:** `cl_Node.{hpp,cpp}`, `cl_Element.{hpp,cpp}`, `cl_Block.{hpp,cpp}`, `cl_SideSet.{hpp,cpp}`
- **Topology:** `cl_Edge.{hpp,cpp}`, `cl_Face.{hpp,cpp}`, `cl_Facet.{hpp,cpp}`
- **I/O:** `cl_Mesh_GmshReader.{hpp,cpp}`, `cl_Mesh_BfmFile.{hpp,cpp}`, `cl_Mesh_ExodusWriter.{hpp,cpp}`, `cl_Mesh_VtkWriter.{hpp,cpp}`
- **Parallel:** `cl_Mesh_Partitioner.{hpp,cpp}`, `cl_Mesh_Distributor.{hpp,cpp}`
- **Elements:** `cl_Element_<TYPE>.hpp` (50+ element classes)
- **Enums:** `Mesh_Enums.hpp`
- **Utilities:** `meshtools.hpp`, `cl_Mesh_ConnectivityCalculator.{hpp,cpp}`

---

## External References

### Mesh Generators

- **Gmsh:** Open-source mesh generator with GUI and scripting (https://gmsh.info/)
- **Cubit/Trelis:** Commercial mesh generator (supports Exodus format)

### File Formats

- **Gmsh format:** MSH file format specification (https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format)
- **Exodus II:** Sandia finite element data format (https://gsjaardema.github.io/seacas/)
- **HDF5:** Hierarchical Data Format (https://www.hdfgroup.org/)

### Partitioning Libraries

- **METIS:** Graph partitioning library (http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
- **ParMETIS:** Parallel graph partitioning (http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview)
- **SCOTCH:** Alternative partitioning library (https://www.labri.fr/perso/pelegrin/scotch/)

### Visualization

- **ParaView:** Open-source visualization (https://www.paraview.org/)
- **VisIt:** DOE visualization tool (https://visit-dav.github.io/visit-website/)

---

## Related BELFEM Modules

- **Containers** (`src/containers/`): `Cell`, `Map` containers for mesh entities
- **Core** (`src/core/`): `typedefs.hpp` (id_t, index_t, proc_t), `Logger`, `Timer`
- **Communication** (`src/comm/`): MPI utilities for parallel mesh operations
- **Graph** (`src/math/graph/`): Graph algorithms for partitioning (METIS, SCOTCH)
- **Homology** (`src/homology/`): Cohomology cuts, node duplication, topological analysis
- **FEM Kernel** (`src/fem/kernel/`): DOF management, domain integration (mesh consumers)
- **FEM Interpolation** (`src/fem/interpolation/`): Shape functions for element types

---

## See Also

- **Project README:** `../../../README.md`
- **Claude Instructions:** `../../../CLAUDE.md`
- **Documentation Guidelines:** `../../../doc/documentation_guidelines.md`
- **Coding Philosophy:** `../../../doc/coding_philosophy.md`
- **General Documentation:** `../../../doc/README.md`
