# Graph Module Documentation {#math_graph_index}

**Module:** src/math/graph
**Purpose:** Index of documentation for BELFEM's graph algorithms and partitioning module

---

## Overview

The `graph` module provides graph algorithms and external library integrations for:

- Mesh reordering (bandwidth reduction, nested dissection)
- Domain decomposition (graph partitioning for MPI parallelization)
- Connectivity analysis (connected components, BFS/DFS traversal)

---

## Documentation Files

### Module Documentation

- **[graph_usage_guide.md](graph_usage_guide.md)** - Comprehensive guide to the graph module
  - Architecture and file organization
  - Vertex class (adjacency lists, manual memory management)
  - Graph algorithms (BFS, DFS, RCM, pseudo-peripheral finding, connected partitions)
  - External library integrations (METIS, ParMETIS, SCOTCH, PT-SCOTCH)
  - CSR adjacency format for external libraries
  - Usage examples and performance considerations
  - Thread safety and MPI awareness
  - Common pitfalls
  - Literature references

---

## Quick Reference

### Key Classes

| Class | File | Purpose |
|-------|------|------------|
| `Vertex` | cl_Graph_Vertex.{hpp,cpp} | Graph node with adjacency list |

### Graph Type

```cpp
typedef Cell<graph::Vertex*> Graph;
```

### Key Algorithms

| Algorithm | File | Purpose | Complexity |
|-----------|------|---------|------------|
| **BFS** | fn_Graph_bfs.{hpp,cpp} | Breadth-first search, level computation | O(V + E) |
| **DFS** | fn_Graph_dfs.{hpp,cpp} | Depth-first search, connected components | O(V + E) |
| **RCM** | fn_Graph_symrcm.{hpp,cpp} | Reverse Cuthill-McKee bandwidth reduction | O(V + E) traversal + O(V log Δ) neighbor sort + O(V log V) final reorder; up to O(V²) with many components |
| **Pseudo-peripheral** | fn_Graph_find_pseudo_peripheral_vertex.{hpp,cpp} | Find vertex with large eccentricity | O(k(V + E)), k ≈ 3 |
| **Connected partitions** | fn_Graph_find_connected_partitions.{hpp,cpp} | Find and sort components by size | O(V + E + P log P) |

### External Library Integrations

| Library | File | Purpose | Compilation Flag |
|---------|------|---------|------------------|
| **METIS** | fn_Graph_METIS.{hpp,cpp} | Nested dissection, partitioning | `BELFEM_METIS` |
| **ParMETIS** | fn_Graph_ParMETIS.{hpp,cpp} | Parallel nested dissection (ordering) | `BELFEM_PARMETIS` |
| **SCOTCH** | fn_Graph_SCOTCH.{hpp,cpp} | Alternative partitioning | `BELFEM_SCOTCH` |
| **PT-SCOTCH** | fn_Graph_PTSCOTCH.{hpp,cpp} | Parallel SCOTCH | `BELFEM_PTSCOTCH` |

### Common Operations

```cpp
// Bandwidth reduction via RCM
Graph tGraph = /* mesh connectivity */;
symrcm(tGraph);  // Reorders in-place

// Nested dissection via METIS (fill-minimizing for direct solvers)
#ifdef BELFEM_METIS
    metis_nd(tGraph);  // Optimal for sparse factorization
#endif

// Graph partitioning for MPI domain decomposition
#ifdef BELFEM_METIS
    metis_partition(tGraph, comm_size(), true);
    // Each vertex->owner() now holds MPI rank
#endif

// Connected component analysis
index_t tMainDomainSize = find_connected_partitions(tGraph);
// Largest component has owner() == 0

// BFS traversal
Vertex* tStart = tGraph(0);
index_t tMaxWidth = bfs(tGraph, tStart);
// Each vertex->level() holds distance from start

// DFS to find components
proc_t tNumComponents = dfs(tGraph);
// Each vertex->owner() holds component ID
```

---

## Algorithm Selection Guide

| Objective | Recommended Algorithm | Notes |
|-----------|----------------------|-------|
| **Bandwidth reduction** | `symrcm()` | For banded direct solvers |
| **Fill-in minimization** | `metis_nd()` | For sparse factorization (MUMPS, STRUMPACK) |
| **Parallel factorization** | `metis_ndp(tGraph, comm_size())` | Top-level partitions before nested dissection |
| **MPI domain decomposition** | `metis_partition()` | Balanced partitioning with minimal edge cuts |
| **Connectivity check** | `find_connected_partitions()` | Identifies disconnected components |
| **Distance computation** | `bfs()` | Level = shortest path distance |

---

## Vertex Class Quick Reference

### Properties

```cpp
id_t     id() const;                    // Permanent identifier
index_t  index() const;                 // Reordering position
proc_t   owner() const;                 // MPI rank or component ID
index_t  level() const;                 // BFS/DFS depth
bool     is_flagged() const;            // Visited flag
uint     number_of_vertices() const;    // Degree (number of neighbors)
Vertex*  vertex(uint k);                // k-th neighbor
```

### Construction Pattern

```cpp
// Two-pass construction (count neighbors, then allocate)
Vertex* v = new Vertex();

// Pass 1: Count
v->increment_vertex_counter();  // For each neighbor

// Pass 2: Allocate and populate
v->init_vertex_container();     // malloc based on counter
v->insert_vertex(neighbor);     // Add each neighbor

// Cleanup
delete v;  // Calls ~Vertex(), which frees adjacency list
```

---

## CSR Adjacency Format

External libraries (METIS, SCOTCH) require Compressed Sparse Row (CSR) format:

```cpp
#include "graphtools.hpp"

Vector<metis_t> tVertices;  // Row pointers: size V+1
Vector<metis_t> tEdges;     // Column indices: size E

build_graph_adjacency(tGraph, tVertices, tEdges);

// tVertices[i] = start index in tEdges for vertex i
// tEdges[tVertices[i] : tVertices[i+1]] = neighbors of vertex i
```

---

## Comparison Operators

For sorting graphs by vertex properties:

```cpp
sort(tGraph, opVertexIndex);   // Sort by index()
sort(tGraph, opVertexID);      // Sort by id()
sort(tGraph, opVertexOwner);   // Sort by owner()
sort(tGraph, opVertexDegree);  // Sort by number_of_vertices()
sort(tGraph, opVertexLevel);   // Sort by level()
```

---

## Source Code

**Module location:** `../../`

**Key source files:**

- **Vertex class**: `cl_Graph_Vertex.{hpp,cpp}`
- **Type definitions**: `graph_typedefs.hpp`, `graphtools.hpp`
- **Traversal**: `fn_Graph_bfs.{hpp,cpp}`, `fn_Graph_dfs.{hpp,cpp}`
- **Reordering**: `fn_Graph_symrcm.{hpp,cpp}`, `fn_Graph_find_pseudo_peripheral_vertex.{hpp,cpp}`
- **Components**: `fn_Graph_find_connected_partitions.{hpp,cpp}`
- **METIS**: `fn_Graph_METIS.{hpp,cpp}`, `fn_Graph_ParMETIS.{hpp,cpp}`
- **SCOTCH**: `fn_Graph_SCOTCH.{hpp,cpp}`, `fn_Graph_PTSCOTCH.{hpp,cpp}`
- **Utilities**: `fn_Graph_sort.{hpp,cpp}`, `fn_Graph_clear.{hpp,cpp}`
- **Operators**: `op_Graph_Vertex_{Index,ID,Owner,Level,Degree}.hpp`

---

## External References

### Classic Graph Algorithms

**Cuthill-McKee:**
- Cuthill & McKee (1969): "Reducing the Bandwidth of Sparse Symmetric Matrices", ACM Conference
- Liu & Sherman (1976): "Comparative Analysis of the Cuthill-McKee and the Reversed Cuthill-McKee Ordering Algorithms", SIAM J. Numer. Anal.
- Referenced in **Bathe** (§8.2.3) and **Hughes** (FEM textbooks)

**Nested Dissection:**
- George (1973): "Nested Dissection of a Regular Finite Element Mesh", SIAM J. Numer. Anal.
- METIS implements multilevel nested dissection (Karypis & Kumar, 1998)

### External Libraries

- **METIS**: [http://glaros.dtc.umn.edu/gkhome/metis/metis/overview](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
- **ParMETIS**: [http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview)
- **SCOTCH**: [https://www.labri.fr/perso/pelegrin/scotch/](https://www.labri.fr/perso/pelegrin/scotch/)

### Related BELFEM Modules

- **Containers** (`src/containers/`): `Cell`, `Queue`, `DynamicBitset` used in graph algorithms
- **Core** (`src/core/`): `typedefs.hpp` (index_t, id_t, proc_t), `Logger`, `Timer`
- **Communication** (`src/comm/`): MPI utilities for ParMETIS/PT-SCOTCH

---

## Development Notes

### Adding New Graph Algorithms

When implementing new graph algorithms:

1. Follow BELFEM naming conventions (see `doc/coding_philosophy.md`)
2. Use `fn_Graph_<algorithm>.{hpp,cpp}` naming pattern
3. Operate on `Graph` (Cell of Vertex pointers)
4. Document time/space complexity
5. Handle disconnected graphs if applicable
6. Use vertex flags/level for visited tracking (reset at start)

### Integrating New Partitioning Libraries

To add new external library (e.g., KaHIP, Zoltan):

1. Add type definition in `graph_typedefs.hpp` with `#ifdef`
2. Create `fn_Graph_<LIBRARY>.{hpp,cpp}` wrapper
3. Use `build_graph_adjacency<T>()` template for CSR conversion
4. Add fallback to RCM if library unavailable
5. Document compilation flag in README

### Performance Profiling

Use `cl_Profiler` and `cl_Timer` for algorithm benchmarking:

```cpp
Timer tTimer;

symrcm(tGraph);

message(InfoLevel::Verbose, "RCM reordering: %u ms", (uint) tTimer.stop());
```

---

## Common Pitfalls

1. **Forgetting continuous indices**: Most algorithms require `vertex->index()` in `[0, N)`. Set before calling algorithms.

2. **Mixing ID and index**: `id()` is permanent identifier, `index()` is reordering position. Use mapping `id → index` after reordering.

3. **Memory leaks**: Call `clear(tGraph)` or manually `delete` each vertex to free adjacency lists.

4. **Thread safety**: Graph algorithms modify vertex properties (not thread-safe). Use separate graphs per thread.

5. **METIS availability**: Check `#ifdef BELFEM_METIS` before calling METIS functions. Provide RCM fallback.

6. **Disconnected graphs**: Use `find_connected_partitions()` to check connectivity before partitioning.

---

## See Also

- **Project README**: `../../../README.md`
- **Claude Instructions**: `../../../CLAUDE.md`
- **Documentation Guidelines**: `../../../doc/documentation_guidelines.md`
- **Coding Philosophy**: `../../../doc/coding_philosophy.md`
- **General Documentation**: `../../../doc/README.md`
