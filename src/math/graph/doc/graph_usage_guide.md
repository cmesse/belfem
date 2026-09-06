# Graph Module Documentation {#math_graph_graph_usage_guide}

**Module:** `src/math/graph`
**Purpose:** Graph algorithms and external partitioning library integrations for FEM mesh ordering and domain decomposition
**Date:** January 16, 2026
**Last Updated:** January 16, 2026

---

## Table of Contents

0. [Glossary](#glossary)

1. [Overview](#overview)
2. [Architecture](#architecture)
3. [Vertex Class](#vertex-class)
4. [Graph Type Definition](#graph-type-definition)
5. [Graph Algorithms](#graph-algorithms)
   - [Breadth-First Search (BFS)](#breadth-first-search-bfs)
   - [Depth-First Search (DFS)](#depth-first-search-dfs)
   - [Reverse Cuthill-McKee (RCM)](#reverse-cuthill-mckee-rcm)
   - [Pseudo-Peripheral Vertex Finding](#pseudo-peripheral-vertex-finding)
   - [Connected Partitions](#connected-partitions)
6. [External Library Integrations](#external-library-integrations)
   - [METIS](#metis-integration)
   - [ParMETIS](#parmetis-integration)
   - [SCOTCH and PT-SCOTCH](#scotch-integration)
7. [Utility Functions](#utility-functions)
8. [CSR Adjacency Format](#csr-adjacency-format)
9. [Usage Examples](#usage-examples)
10. [Thread Safety and MPI](#thread-safety-and-mpi)
11. [Common Pitfalls](#common-pitfalls)
12. [Literature References](#literature-references)
13. [Performance Considerations](#performance-considerations)

---

## Glossary

### Core Types

| Type | Definition | Notes |
|------|------------|-------|
| `Graph` | `typedef Cell<graph::Vertex*> Graph;` | Dynamic array of vertex pointers |
| `Vertex` | `graph::Vertex` | Graph node with adjacency list |
| `Cell<T>` | BELFEM dynamic array | See `src/containers/doc/` |
| `Queue<T>` | BELFEM FIFO queue | Used in BFS traversal |
| `DynamicBitset` | BELFEM bitset | Visited tracking in algorithms |
| `Vector<T>` | BELFEM vector | CSR arrays, partitions |

### Primitive Types (from `core/typedefs.hpp`)

| Type | Definition | Purpose |
|------|------------|---------|
| `index_t` | `uint32_t` or `uint64_t` | Array indexing, ordering positions |
| `id_t` | `unsigned int` (`typedefs.hpp:41`) | Permanent unique identifiers |
| `proc_t` | `int` | MPI rank, partition/component ID |
| `uint` | `unsigned int` | Counts, degrees |
| `real` | `double` | Floating-point (not used in graph module) |

### Sentinel Values

| Constant | Type | Value | Meaning |
|----------|------|-------|---------|
| `gNoIndex` | `index_t` | `std::numeric_limits<index_t>::max()` | Invalid/unset index |
| `gNoID` | `id_t` | `std::numeric_limits<id_t>::max()` | Invalid/unset ID |
| `gNoOwner` | `proc_t` | `numeric_limits<proc_t>::max()` (`typedefs.hpp:59`) | No owner assigned |

### External Library Types

| Type | Definition | Purpose |
|------|------------|---------|
| `metis_t` | `idx_t` (from METIS) | METIS integer type (32 or 64-bit) |
| `scotch_t` | `SCOTCH_Num` (from SCOTCH) | SCOTCH integer type |

### Naming Conventions (from `doc/coding_philosophy.md`)

| Prefix | Scope | Example |
|--------|-------|---------|
| `a` | Argument/parameter | `aGraph`, `aStart` |
| `t` | Temporary/local variable | `tVertex`, `tCount` |
| `m` | Member variable | `mVertices`, `mID`, `mIndex` |
| `g` | Global variable | `gNoIndex`, `gNoID` |

---

## Overview

The `graph` module provides graph algorithms and partitioning tools for finite element mesh operations. The module serves three primary purposes:

1. **Bandwidth Reduction** — Reverse Cuthill-McKee (RCM) ordering to reduce matrix bandwidth for direct solvers
2. **Nested Dissection** — METIS/SCOTCH integration for fill-minimizing orderings optimal for sparse factorization
3. **Domain Decomposition** — Graph partitioning for parallel computing (MPI-based mesh distribution)

### Key Features

- Classic graph traversal algorithms (BFS, DFS)
- Bandwidth reduction via RCM ordering
- Integration with METIS, ParMETIS, SCOTCH, and PT-SCOTCH
- Connected component analysis
- CSR (Compressed Sparse Row) adjacency structure building
- Manual memory management for cache-efficient vertex storage

### Design Philosophy

The graph module follows BELFEM's coding philosophy (see `doc/coding_philosophy.md`):

- **Manual memory management**: Vertex adjacency lists use raw pointers (`malloc`/`free`) for contiguous storage
- **Zero-overhead abstraction**: Graph algorithms operate directly on `Cell<Vertex*>` with minimal indirection
- **External library isolation**: Conditional compilation (`#ifdef BELFEM_METIS`) allows building without external dependencies
- **Prefix naming**: Temporary variables `t`, arguments `a`, members `m` (Vertex class)

---

## Architecture

### Module Structure

```
src/math/graph/
├── cl_Graph_Vertex.{hpp,cpp}              # Core vertex class
├── graph_typedefs.hpp                     # METIS/SCOTCH type definitions
├── graphtools.hpp                         # CSR adjacency building templates
│
├── fn_Graph_bfs.{hpp,cpp}                 # Breadth-first search
├── fn_Graph_dfs.{hpp,cpp}                 # Depth-first search
├── fn_Graph_symrcm.{hpp,cpp}              # Reverse Cuthill-McKee ordering
├── fn_Graph_find_pseudo_peripheral_vertex.{hpp,cpp}
├── fn_Graph_find_pseudo_peripheral_node.{hpp,cpp}
├── fn_Graph_find_connected_partitions.{hpp,cpp}
├── fn_Graph_sort.{hpp,cpp}                # Graph sorting by vertex properties
│
├── fn_Graph_METIS.{hpp,cpp}               # METIS nested dissection & partitioning
├── fn_Graph_ParMETIS.{hpp,cpp}            # Parallel METIS integration
├── fn_Graph_SCOTCH.{hpp,cpp}              # SCOTCH partitioning
├── fn_Graph_PTSCOTCH.{hpp,cpp}            # Parallel SCOTCH
│
├── fn_Graph_clear.{hpp,cpp}               # Graph cleanup
└── op_Graph_Vertex_{Index,ID,Owner,Level,Degree}.hpp  # Comparison operators
```

### Dependencies

**BELFEM Internal:**
- `core/typedefs.hpp` — `index_t`, `id_t`, `proc_t`, sentinel values
- `containers/cl_Cell.hpp` — Dynamic array for graph storage
- `containers/cl_Queue.hpp` — BFS/DFS traversal
- `containers/cl_DynamicBitset.hpp` — Visited tracking
- `linalg/cl_Vector.hpp` — CSR adjacency arrays
- `comm/commtools.hpp` — MPI utilities (`comm_size()`, `comm_rank()`)

**External Libraries (Optional):**
- **METIS** — Serial graph partitioning and nested dissection (`BELFEM_METIS`)
- **ParMETIS** — Parallel nested dissection (ordering) (`BELFEM_PARMETIS`)
- **SCOTCH** — Alternative serial partitioning (`BELFEM_SCOTCH`)
- **PT-SCOTCH** — Parallel SCOTCH (`BELFEM_PTSCOTCH`)

---

## Vertex Class

### Declaration

`src/math/graph/cl_Graph_Vertex.{hpp,cpp}`

```cpp
namespace belfem::graph {
    class Vertex {
        id_t     mID;            // Vertex identifier
        index_t  mIndex;         // Position/reordering index
        proc_t   mOwner;         // MPI rank owner (partitioning)
        index_t  mLevel;         // BFS/DFS level (traversal depth)
        uint8_t  mFlags;         // 8 marking flags (bitset), see flag()/unflag()/is_flagged()

        uint32_t mVertexCounter; // Number of adjacent vertices (overflow asserted)
        Vertex** mVertices;      // Adjacency list (raw pointer array)
    };
}
```

### Properties

| Property | Type | Purpose | Sentinel Value |
|----------|------|---------|----------------|
| `mID` | `id_t` | Permanent vertex identifier | `gNoID` |
| `mIndex` | `index_t` | Current position/ordering | `gNoIndex` |
| `mOwner` | `proc_t` | MPI rank (partitioning) | `gNoOwner` (`cl_Graph_Vertex.hpp:42`) — see the note below on the `comm_size()` convention |
| `mLevel` | `index_t` | Traversal depth (BFS/DFS) | `0` (`cl_Graph_Vertex.hpp:45`) |
| `mFlags` | `uint8_t` | Bitset of marking flags (`cl_Graph_Vertex.hpp:48`) | `0` |
| `mVertexCounter` | `uint` | Degree (number of neighbors) | 0 |
| `mVertices` | `Vertex**` | Adjacency list | `nullptr` |

> **`mOwner` initializes to `gNoOwner`** (`cl_Graph_Vertex.hpp:42`), which is
> `numeric_limits<proc_t>::max()`. Until 2026-08-31 it was initialized to `gNoID` — a different
> type and a different value (`numeric_limits<id_t>::max()` over `unsigned int`) — which narrowed
> to an out-of-range `proc_t` and made `owner() == gNoOwner` false on a fresh vertex.
>
> That the corrected value is the intended one is visible in the algorithms rather than argued
> from naming: `fn_Graph_dfs.cpp:32` resets with `set_owner( gNoOwner )`, and the `std::min`
> ownership sweeps in `cl_Mesh_Partitioner.cpp:253-260` and `cl_FEM_Kernel.cpp:427` only work if
> the sentinel behaves as **+∞** — under the old value it was an absorbing element and every
> facet would have ended unassigned.
>
> **`gNoOwner` is not `comm_size()`.** A separate convention marks a deliberately unassigned
> entity with `owner() == comm_size()` (`cl_FEM_Kernel.cpp:195,203`), and that is what the guard
> at `graphtools.hpp:142` tests. The two mean different things: "never set" versus "set, and
> known to have no home". Do not conflate them.

### Key Methods

```cpp
// Property access
void     set_id(id_t aID);
id_t     id() const;

void     set_index(index_t aIndex);
index_t  index() const;

void     set_owner(proc_t aOwner);
proc_t   owner() const;

void     set_level(index_t aLevel);
index_t  level() const;

void     flag();
void     unflag();
bool     is_flagged() const;

// Adjacency list management
void     increment_vertex_counter();                // Count neighbors first
void     init_vertex_container();                   // Allocate based on counter
void     init_vertex_container(uint aSize);         // Allocate with explicit size
void     reset_vertex_container();                  // Free memory

void     insert_vertex(Vertex* aVertex);            // Add neighbor
uint     number_of_vertices() const;                // Get degree
Vertex*  vertex(uint aIndex);                       // Access neighbor
void     sort_vertices();                           // Sort by index
```

### Memory Management

The vertex adjacency list uses **manual memory management** for performance:

```cpp
void Vertex::init_vertex_container() {
    if (mVertices != nullptr) free(mVertices);  // Free existing

    if (mVertexCounter > 0) {
        mVertices = (Vertex**) malloc(mVertexCounter * sizeof(Vertex*));
    } else {
        mVertices = nullptr;
    }
    mVertexCounter = 0;  // Reset for insertion
}

Vertex::~Vertex() {
    this->reset_vertex_container();  // Cleanup
}
```

**Rationale** (see `doc/coding_philosophy.md`):
- **Contiguous allocation**: Single `malloc` for all neighbor pointers → cache-friendly
- **Exact sizing**: No over-allocation like `std::vector` growth strategy
- **Zero overhead**: No smart pointer bookkeeping or atomic reference counting

### Typical Construction Pattern

```cpp
// Two-pass construction (common in mesh algorithms)
Graph tGraph(num_vertices, nullptr);

// Pass 1: Count neighbors
for (Vertex* v : tGraph) {
    for (/* each edge incident to v */) {
        v->increment_vertex_counter();
    }
}

// Pass 2: Allocate and populate
for (Vertex* v : tGraph) {
    v->init_vertex_container();  // malloc based on counter
    for (Vertex* neighbor : /* neighbors of v */) {
        v->insert_vertex(neighbor);
    }
}
```

---

## Graph Type Definition

```cpp
typedef Cell<graph::Vertex*> Graph;
```

A graph is a `Cell` (BELFEM's dynamic array, see `src/containers/doc/`) containing vertex pointers. This allows:

- **Random access**: `aGraph(i)` returns the i-th vertex
- **Iteration**: `for (Vertex* v : aGraph)`
- **Dynamic resizing**: `aGraph.push(vertex)`
- **Sorting**: `sort(aGraph, opVertexIndex)`

---

## Graph Algorithms

### Breadth-First Search (BFS)

**Files:** `fn_Graph_bfs.{hpp,cpp}`

**Purpose:** Level-by-level graph traversal, computes distances from start vertex, used in pseudo-peripheral vertex finding and RCM.

```cpp
namespace belfem::graph {
    // BFS from specific start vertex
    index_t bfs(Graph& aGraph, Vertex* aStart);   // no default — use the one-argument overload for all components

    // BFS handling disconnected graphs (all components)
    index_t bfs(Graph& aGraph);
}
```

**Return Value:** Maximum level width (largest number of vertices at any single level).

**Side Effects:**
- Sets `vertex->level()` to distance from start vertex (0 = start, 1 = neighbors, etc.)
- Unvisited vertices have `level() == gNoIndex`
- Sets `vertex->index()` to continuous values `[0, N)`

**Algorithm:**

1. Initialize all vertices to `level = gNoIndex`
2. Set start vertex to `level = 0`
3. Use two queues (`tCurrentLevel`, `tNextLevel`) for level-synchronous traversal
4. For each vertex at current level, visit all unvisited neighbors → add to next level
5. Track maximum width across all levels
6. For disconnected variant: Repeat for each unvisited component

**Example:**

```cpp
Graph tGraph = /* construct graph */;
Vertex* tStart = tGraph(0);

index_t tMaxWidth = bfs(tGraph, tStart);

for (Vertex* v : tGraph) {
    if (v->level() != gNoIndex) {
        std::cout << "Vertex " << v->id() << " at distance " << v->level() << std::endl;
    } else {
        std::cout << "Vertex " << v->id() << " unreachable from start" << std::endl;
    }
}
```

**Time Complexity:** O(V + E) where V = vertices, E = edges

**Space Complexity:** O(V) for visited bitset and queues

**Literature:** Classic algorithm (Cormen et al., "Introduction to Algorithms")

---

### Depth-First Search (DFS)

**Files:** `fn_Graph_dfs.{hpp,cpp}`

**Purpose:** Recursive-style graph traversal (implemented with explicit stack), finds connected components.

```cpp
namespace belfem::graph {
    // Find all connected components, assign owner IDs
    proc_t dfs(Graph& aGraph);

    // DFS from specific start vertex (helper)
    void dfs_from_start(Graph& aGraph, Vertex* aStart);
}
```

**Return Value:** `dfs(aGraph)` returns the number of connected components.

**Side Effects:**
- Sets `vertex->owner()` to component ID (0, 1, 2, ...)
- Sets `vertex->level()` to depth in DFS tree from component root
- Flags all visited vertices (`is_flagged() == true`)

**Algorithm:**

1. Reset all vertices: `owner = gNoOwner`, `level = gNoIndex`, `unflag()`
2. For each unvisited vertex:
   - Assign new component ID
   - Run DFS from this vertex using stack
3. DFS traversal:
   - Push start vertex onto stack with `level = 0`
   - Pop vertex, visit all unvisited neighbors
   - Set neighbor `level = parent_level + 1`, `owner = start_owner`
   - Push neighbors onto stack

**Example:**

```cpp
Graph tGraph = /* construct graph */;

proc_t tNumComponents = dfs(tGraph);

std::cout << "Found " << tNumComponents << " connected components" << std::endl;

for (proc_t i = 0; i < tNumComponents; ++i) {
    index_t count = 0;
    for (Vertex* v : tGraph) {
        if (v->owner() == i) ++count;
    }
    std::cout << "Component " << i << " has " << count << " vertices" << std::endl;
}
```

**Time Complexity:** O(V + E)

**Space Complexity:** O(V) for stack and flags

**Use Cases:**
- Connected component detection in mesh connectivity graphs
- Checking mesh partition connectivity
- Finding isolated subdomains

---

### Reverse Cuthill-McKee (RCM)

**Files:** `fn_Graph_symrcm.{hpp,cpp}`

**Purpose:** Reorder graph vertices to reduce bandwidth of adjacency matrix, improving cache performance of banded direct solvers.

```cpp
namespace belfem::graph {
    void symrcm(Graph& aGraph, Vertex* aStart = nullptr);
}
```

**Parameters:**
- `aGraph` — Graph to reorder (modified in-place)
- `aStart` — Starting vertex (if `nullptr`, uses pseudo-peripheral vertex)

**Side Effects:**
- Reorders `aGraph` in-place (vertices sorted by new RCM indices)
- Sets `vertex->index()` to new position in RCM ordering

**Algorithm:**

1. Find starting vertex:
   - If `aStart == nullptr`: find pseudo-peripheral vertex (see next section)
   - Heuristic: if exactly one vertex has minimum degree, the pseudo-peripheral search is seeded from it; otherwise it is seeded from `aGraph(0)`

2. Cuthill-McKee traversal:
   - BFS from start vertex
   - At each level, **sort neighbors by degree (ascending)**
   - This prioritizes low-degree vertices → reduces bandwidth

3. **Reverse the permutation** (key step):
   - RCM ordering is Cuthill-McKee in reverse
   - Empirically reduces bandwidth more than forward CM

4. Apply permutation:
   - Update `vertex->index()` to RCM positions
   - Sort graph by new indices

5. Handle disconnected components:
   - Process remaining unvisited vertices with new CM runs
   - Maintains component separation in final ordering

**Example:**

```cpp
Graph tGraph = /* construct mesh connectivity graph */;

// Reorder for bandwidth reduction
symrcm(tGraph);

// Now tGraph is sorted by RCM ordering: tGraph(i)->index() == i
// tGraph(0) is the last vertex visited by Cuthill-McKee (RCM index 0)
// tGraph(N-1) is the start vertex (RCM index N-1)

// Use for matrix assembly with reduced bandwidth
for (index_t i = 0; i < tGraph.size(); ++i) {
    index_t tNewIndex = tGraph(i)->index();
    // Assemble matrix using tNewIndex
}
```

**Time Complexity:** three terms, and the last one dominates on a well-connected mesh.
O(V + E) for the BFS traversal; **O(V log Δ)** for the neighbor sort — each vertex is sorted into
exactly one list, of its unvisited neighbors (`fn_Graph_symrcm.cpp:126-129`); and **O(V log V)** for
the final `sort( aGraph, opVertexIndex )` that reorders the whole graph (`:210`). The
disconnected-component scan is a nested search over remaining vertices (`:138-150`), so a graph of
many small components degrades toward O(V²).

**Space Complexity:** O(V) for visited bitset, queue, permutation

**Literature References:**
- **Cuthill & McKee** (1969): "Reducing the Bandwidth of Sparse Symmetric Matrices", ACM Conference
- **Liu & Sherman** (1976): "Comparative Analysis of the Cuthill-McKee and the Reversed Cuthill-McKee Ordering Algorithms for Sparse Matrices", SIAM J. Numer. Anal.
- **Bathe**: Bandwidth reduction discussion in matrix storage (FEM textbook)
- **Hughes**: References to RCM in sparse matrix algorithms

**Use Cases:**
- Preprocessing for banded direct solvers (LAPACK)
- Reducing cache misses in iterative solvers
- Improving matrix-vector product locality
- Fallback when METIS nested dissection unavailable

---

### Pseudo-Peripheral Vertex Finding

**Files:** `fn_Graph_find_pseudo_peripheral_vertex.{hpp,cpp}`

**Purpose:** Find a vertex with large eccentricity (distance to farthest vertex), used as RCM starting point.

```cpp
namespace belfem::graph {
    Vertex* find_pseudo_peripheral_vertex(
        Graph& aGraph,
        Vertex* aStart = nullptr
    );
}
```

**Parameters:**
- `aGraph` — Input graph
- `aStart` — Initial guess (if `nullptr`, uses first vertex)

**Return Value:** Pointer to a pseudo-peripheral vertex.

**Algorithm:**

1. Start with initial vertex (argument or `aGraph(0)`)
2. Loop until no improvement:
   - BFS from current vertex to find farthest vertices
   - Compute eccentricity (maximum level reached)
   - If eccentricity ≤ previous: **stop, return current**
   - Otherwise: Among farthest vertices, choose one with **minimum degree**
   - Repeat BFS from this new candidate

**Heuristic Rationale:**
- Peripheral vertices tend to have large eccentricity
- Minimum degree among farthest vertices avoids high-degree "hubs"
- Usually converges in 2-4 iterations for mesh graphs

**Example:**

```cpp
Graph tGraph = /* construct graph */;

Vertex* tPeripheral = find_pseudo_peripheral_vertex(tGraph);

std::cout << "Pseudo-peripheral vertex: ID = " << tPeripheral->id()
          << ", degree = " << tPeripheral->number_of_vertices() << std::endl;

// Use for RCM starting point
symrcm(tGraph, tPeripheral);
```

**Time Complexity:** O(k(V + E)) where k = number of iterations (typically k ≈ 3)

**Literature:** Gibbs, Poole, Stockmeyer (1976) — standard heuristic for RCM start vertex

---

### Connected Partitions

**Files:** `fn_Graph_find_connected_partitions.{hpp,cpp}`

**Purpose:** Find connected components and **sort by size** (largest first), useful for identifying main mesh domain vs. disconnected islands.

```cpp
namespace belfem::graph {
    index_t find_connected_partitions(Graph& aGraph);
}
```

**Return Value:** Number of vertices in the **largest** component.

**Side Effects:**
- Finds connected components via DFS
- Sorts components by size (descending)
- Sets `vertex->owner()` to component ID (0 = largest, 1 = second-largest, ...)
- **Reorders `aGraph`** to group vertices by component

**Algorithm:**

1. Ensure continuous indices: set `vertex->index()` to `[0, N)`
2. Run DFS to find components (sets `vertex->owner()`)
3. Count vertices per component
4. **Sort components by size** (descending):
   - Create pairs `(size, component_id)`
   - Sort using `std::sort` with `std::greater`
5. Remap owner IDs: largest component → owner 0, etc.
6. **Sort graph** by owner using `opVertexOwner` comparator
7. Return size of largest component

**Example:**

```cpp
Graph tGraph = /* construct mesh graph with potential islands */;

index_t tMainDomainSize = find_connected_partitions(tGraph);

// Count components
proc_t tNumComponents = 0;
for (Vertex* v : tGraph) {
    if (v->owner() + 1 > tNumComponents) {
        tNumComponents = v->owner() + 1;
    }
}

std::cout << "Main domain: " << tMainDomainSize << " vertices" << std::endl;
std::cout << "Total components: " << tNumComponents << std::endl;

// Process main component (owner == 0)
for (Vertex* v : tGraph) {
    if (v->owner() == 0) {
        // Main domain vertex
    } else {
        // Disconnected island
        std::cout << "Warning: Vertex " << v->id()
                  << " in disconnected component " << v->owner() << std::endl;
    }
}
```

**Time Complexity:** O(V + E + P log P) where P = number of components

**Use Cases:**
- Detecting mesh quality issues (disconnected elements)
- Separating main domain from boundary layer islands
- Validating mesh connectivity before partitioning

---

## External Library Integrations

### METIS Integration

**Files:** `fn_Graph_METIS.{hpp,cpp}`

**Compilation Flag:** `BELFEM_METIS`, set by `USE_METIS` (default **ON**, `CMakeLists.txt:106`). METIS is found directly by `config/linalg/config_metis.cmake` — it is **not** linked through SuiteSparse, whose own option defaults OFF.

**Type Definition:** `graph_typedefs.hpp` defines `metis_t` as METIS's `idx_t` integer type

#### Nested Dissection (NodeND)

**Purpose:** Reorder graph for **fill-minimizing factorization** of sparse matrices, optimal for direct solvers like MUMPS, STRUMPACK.

```cpp
namespace belfem::graph {
    void metis_nd(Graph& aGraph);
}
```

**Side Effects:**
- Reorders `aGraph` in-place using METIS nested dissection
- Sets `vertex->index()` to new position in ND ordering
- Sorts `aGraph` by new indices

**METIS Options Used:**
- `METIS_OPTION_NUMBERING = 0` (zero-based indexing)
- `METIS_OPTION_COMPRESS = 0` (no graph compression)
- `METIS_OPTION_CTYPE = METIS_CTYPE_SHEM` (sorted heavy-edge matching coarsening)
- `METIS_OPTION_RTYPE = METIS_RTYPE_SEP1SIDED` (one-sided refinement for separators)

**Algorithm:** METIS_NodeND (external library):
1. Recursively partition graph into two subgraphs + separator
2. Order separator vertices last (high elimination order)
3. Recursively order each subgraph
4. Produces shallow elimination tree → less fill-in

**Example:**

```cpp
#ifdef BELFEM_METIS
    Graph tGraph = /* mesh connectivity graph */;

    metis_nd(tGraph);

    // Use for sparse factorization
    // Vertices ordered to minimize fill-in during Cholesky/LU
#else
    // Fallback if METIS unavailable
    symrcm(tGraph);
#endif
```

**Time Complexity:** O(E) average case (METIS internal, highly optimized)

**Literature:**
- Karypis & Kumar (1998): "METIS: A Software Package for Partitioning Unstructured Graphs"
- STRUMPACK documentation recommends NodeND or NodeNDP for optimal performance

#### Nested Dissection with Partitioning (NodeNDP)

**Purpose:** Variant of NodeND that creates **top-level partitions** before nested dissection, often better for parallel factorization.

```cpp
namespace belfem::graph {
    void metis_ndp(Graph& aGraph, uint aNumPartitions);
}
```

**Parameters:**
- `aNumPartitions` — Number of top-level domains before nested dissection

**STRUMPACK Note:** Documentation suggests `--sp_enable_METIS_NodeNDP` may work better than NodeND for certain problem types.

**Algorithm:** METIS_NodeNDP:
1. Partition graph into `aNumPartitions` subgraphs
2. Apply nested dissection within each partition
3. Produces balanced elimination tree suitable for parallel factorization

**Example:**

```cpp
#ifdef BELFEM_METIS
    Graph tGraph = /* mesh graph */;

    // Use number of MPI ranks for top-level partitions
    uint tNumPartitions = comm_size();

    metis_ndp(tGraph, tNumPartitions);
#endif
```

#### Graph Partitioning (PartGraphKway)

**Purpose:** Partition graph into **balanced subdomains** with minimal edge cuts, for MPI domain decomposition.

```cpp
namespace belfem::graph {
    void metis_partition(
        Graph& aGraph,
        uint aNumPartitions,
        bool aForceContinuousPartitions = true,
        Vector<proc_t>* aPartitions = nullptr
    );
}
```

**Parameters:**
- `aNumPartitions` — Number of partitions (typically `comm_size()`)
- `aForceContinuousPartitions` — If true, each partition is connected (slower but better quality)
- `aPartitions` — If `!= nullptr`, store partition IDs here instead of `vertex->owner()`

**Side Effects:**
- Sets `vertex->owner()` to partition ID (0, 1, ..., aNumPartitions-1)
- If `aPartitions` provided, stores partition IDs there instead

**METIS Options:**
- `METIS_OPTION_CONTIG` — Forces continuous partitions
- `METIS_OPTION_OBJTYPE = METIS_OBJTYPE_VOL` — Minimize communication volume

**Example:**

```cpp
#ifdef BELFEM_METIS
    Graph tGraph = /* mesh graph */;
    uint tNumProcs = comm_size();

    metis_partition(tGraph, tNumProcs, true);

    // Each vertex now has owner() set to MPI rank
    for (Vertex* v : tGraph) {
        proc_t tRank = v->owner();
        // Distribute vertex to rank tRank
    }
#endif
```

**Time Complexity:** O(E) with K-way refinement

**Use Cases:**
- MPI mesh distribution
- Load balancing for parallel FEM
- Minimizing inter-processor communication

---

### ParMETIS Integration

**Files:** `fn_Graph_ParMETIS.{hpp,cpp}`

**Compilation Flag:** `BELFEM_PARMETIS`

**Purpose:** parallel nested-dissection ordering through ParMETIS.

**Input model — note this is not the usual ParMETIS contract.** BELFEM's wrapper does *not* take
an already-distributed graph. The **complete** graph must be on rank 0: `build_pargraph_adjacency()`
asserts `comm_rank() == 0` (`graphtools.hpp:171`), builds every rank's CSR slice there, and
distributes them from within `parmetis_nd()` (`fn_Graph_ParMETIS.cpp:57-99`, root branch). The public entry
point is `parmetis_nd( Graph & )`, and it is **collective**: every rank calls it, non-root ranks
with an empty `Graph`.

Every vertex must carry an owner in `[0, comm_size())` before the call. The builder maps the two
unowned sentinels (`comm_size()` — the Kernel's marker — and the `gNoOwner` default) to rank 0 and
refuses anything else with an always-on error. For a graph nobody has partitioned (the sparse-matrix
graph `DistMatrix` builds), `graph::block_distribution( aGraph, comm_size() )` assigns owners in
contiguous blocks of graph order — `N/P` vertices on the first `P − (N mod P)` ranks, one more on
the rest, the same split `DistMatrix` uses for its rows. That is a *working* distribution for the
ordering, not a partition.

ParMETIS refuses a rank with no vertex. `parmetis_nd()` checks the distribution on root, broadcasts
the verdict to every rank before any slice is sent, and on a failed verdict prints a warning and
falls back to serial METIS on root while the other ranks return; no rank reaches `ParMETIS_V3_NodeND`
with an empty slice. PT-Scotch accepts empty ranks, so `ptscotch_nd()` has no such guard.

What wiring does **not** change: root still holds the complete graph. The parallel wrappers
parallelize the ordering *compute* after a root-side scatter; they do not reduce rank-0 memory.

**Key Functions:**
- Uses ParMETIS for parallel nested-dissection ordering (`ParMETIS_V3_NodeND`)
- Requires MPI communicator setup
- More complex adjacency structures (per-process CSR arrays)

**Template Function:**

```cpp
template<typename T>
void build_pargraph_adjacency(
    Graph& aGraph,
    Vector<T>& aDistribution,
    Cell<Vector<T>>& aVertices,
    Cell<Vector<T>>& aEdges
);
```

**Purpose:** Build distributed CSR adjacency for ParMETIS (called by root process).

**Algorithm:**
1. Reorder graph by owner (MPI rank)
2. Compute distribution array: `aDistribution(p+1) - aDistribution(p)` = vertices on rank p
3. Build per-process CSR arrays (vertices and edges)
4. Validate all indices are set (debug mode)

**Use Cases:**
- Ordering time on large systems when the ranks are already there (the PETSc path at
  `comm_size() > 1` with `reordering scheme : parmetis` or `ptscotch`)
- Not a memory remedy: the complete graph is still built and held on rank 0

---

### SCOTCH Integration

**Files:** `fn_Graph_SCOTCH.{hpp,cpp}`, `fn_Graph_PTSCOTCH.{hpp,cpp}`

**Compilation Flags:** `BELFEM_SCOTCH` (serial), `BELFEM_PTSCOTCH` (parallel)

**Type Definition:** `scotch_t` as `SCOTCH_Num` (SCOTCH integer type)

**Purpose:** Alternative to METIS for graph partitioning and ordering, often used in European HPC codes.

**Features:**
- `fn_Graph_SCOTCH.cpp` — Serial SCOTCH partitioning
- `fn_Graph_PTSCOTCH.cpp` — Parallel PT-SCOTCH

**Differences from METIS:**
- Different algorithm (spectral methods + greedy refinement)
- Sometimes produces better partitions for specific graph types
- Used when METIS license restrictions apply

**Typical Usage:** Similar API to METIS functions (implementation details in `.cpp` files).

---

## Utility Functions

### Graph Sorting

**Files:** `fn_Graph_sort.{hpp,cpp}`, `op_Graph_Vertex_*.hpp`

**Purpose:** Sort graph vertices by various properties.

```cpp
namespace belfem::graph {
    void sort(Graph& aGraph);  // Default sort
}
```

**Comparison Operators:**

```cpp
// src/math/graph/op_Graph_Vertex_Index.hpp
struct {
    bool operator()(const Vertex* a, const Vertex* b) {
        return a->index() < b->index();
    }
} opVertexIndex;

// src/math/graph/op_Graph_Vertex_Degree.hpp
struct {
    bool operator()(const Vertex* a, const Vertex* b) {
        return a->number_of_vertices() < b->number_of_vertices();
    }
} opVertexDegree;

// Similar: opVertexID, opVertexOwner, opVertexLevel
```

**Example:**

```cpp
// Sort by index (after RCM reordering)
sort(tGraph, opVertexIndex);

// Sort by degree (for greedy coloring)
sort(tGraph, opVertexDegree);

// Sort by owner (group by MPI rank)
sort(tGraph, opVertexOwner);
```

---

### Graph Clearing

**Files:** `fn_Graph_clear.{hpp,cpp}`

**Purpose:** Deallocate vertex objects and clean up graph.

```cpp
namespace belfem::graph {
    void clear(Cell<Vertex*>& aGraph);
}
```

**Implementation:**

```cpp
void clear(Cell<Vertex*>& aGraph) {
    for (Vertex* v : aGraph) {
        delete v;  // Calls ~Vertex(), which frees adjacency list
    }
    aGraph.clear();  // Clear the Cell
}
```

**When to Use:**
- After converting graph back to mesh data structures
- When graph is temporary (e.g., connectivity analysis)

**Warning:** Do **not** call if vertices are managed elsewhere (e.g., mesh owns them).

---

## CSR Adjacency Format

### Purpose

External libraries (METIS, SCOTCH) require **Compressed Sparse Row (CSR)** format for graph adjacency:

```
Vertices: [0, 2, 5, 8, 10]
Edges:    [1, 3, 0, 2, 4, 1, 3, 5, 2, 4]
```

Interpretation:
- Vertex 0 has neighbors `[Edges[0:2]] = [1, 3]`
- Vertex 1 has neighbors `[Edges[2:5]] = [0, 2, 4]`
- Vertex 2 has neighbors `[Edges[5:8]] = [1, 3, 5]`
- Etc.

### Template Function

```cpp
template<typename T>
void build_graph_adjacency(
    Graph& aGraph,
    Vector<T>& aVertices,
    Vector<T>& aEdges
);
```

**Template Parameter:** `T = metis_t` or `T = scotch_t` (library-specific integer type)

**Algorithm:**

```cpp
// First pass: ensure continuous indices
T tCount = 0;
for (Vertex* v : aGraph) {
    v->set_index(tCount++);
}

// Second pass: count edges, self-loops excluded
T tNumEdges = 0;
for (Vertex* v : aGraph) {
    for (uint k = 0; k < v->number_of_vertices(); ++k) {
        if (v->vertex(k) != v) ++tNumEdges;
    }
}

// Allocate CSR arrays
aVertices.set_size(tNumVertices + 1, 0);
aEdges.set_size(max(tNumEdges, 1), 0);  // At least 1 to avoid NULL

// Third pass: build CSR structure, same predicate
tCount = 0;
for (Vertex* v : aGraph) {
    aVertices(v->index()) = tCount;
    for (uint k = 0; k < v->number_of_vertices(); ++k) {
        if (v->vertex(k) != v) {
            aEdges(tCount++) = static_cast<T>(v->vertex(k)->index());
        }
    }
}
aVertices(tNumVertices) = tCount;  // Sentinel
```

**Self-loops are dropped from the CSR, not from the `Graph`.** A graph built from a sparse matrix
carries one loop per diagonal entry (`create_graph_from_matrix` inserts every column, the diagonal
included). METIS and SCOTCH define their input as loop-free, and `METIS_NodeNDP` corrupts its heap
on a looped graph; ParMETIS passes loops through unchecked. So both builders skip `v->vertex(k) == v`
when counting and filling, while the `Graph` keeps the loop: `DistMatrix` builds the permuted
`SpMatrix` from the `Graph` and needs the diagonal. A vertex whose only neighbour was itself ends up
with an empty CSR row, and a graph or slice with no edges left keeps the one-element placeholder
buffer — the wrappers take the logical edge count from the CSR terminal, never from that buffer.

**Example Usage:**

```cpp
Graph tGraph = /* construct graph */;

Vector<metis_t> tVertices;
Vector<metis_t> tEdges;

build_graph_adjacency(tGraph, tVertices, tEdges);

// Now call METIS_NodeND, METIS_PartGraphKway, etc.
```

---

## Usage Examples

### Example 1: Bandwidth Reduction for Direct Solver

```cpp
#include "graphtools.hpp"
#include "fn_Graph_symrcm.hpp"

// Build mesh connectivity graph
Graph tGraph(num_nodes, nullptr);

// ... populate graph with element adjacencies ...

// Apply RCM reordering
symrcm(tGraph);

// Use reordered indices for matrix assembly
for (Vertex* v : tGraph) {
    index_t tNewIndex = v->index();
    index_t tOldID = v->id();

    // Map old node ID to new index
    mNodeReordering(tOldID) = tNewIndex;
}

// Assemble matrix with reduced bandwidth
// aMatrix(mNodeReordering(i), mNodeReordering(j)) += contribution;
```

---

### Example 2: METIS Nested Dissection

```cpp
#ifdef BELFEM_METIS
    #include "fn_Graph_METIS.hpp"

    Graph tGraph = /* mesh graph */;

    // Apply METIS nested dissection ordering
    metis_nd(tGraph);

    // Extract permutation. metis_nd() does not touch id(); this assumes the
    // vertices were created with id() == original index
    Cell<index_t> tPermutation(tGraph.size());
    for (index_t i = 0; i < tGraph.size(); ++i) {
        index_t tOldIndex = tGraph(i)->id();
        index_t tNewIndex = tGraph(i)->index();
        tPermutation(tOldIndex) = tNewIndex;
    }

    // Apply to DOF manager
    mDofManager->apply_permutation(tPermutation);
#else
    // Fallback to RCM if METIS unavailable
    symrcm(tGraph);
#endif
```

---

### Example 3: MPI Domain Decomposition

```cpp
#include "fn_Graph_METIS.hpp"
#include "commtools.hpp"

Graph tGraph = /* global mesh graph on rank 0 */;
proc_t tNumProcs = comm_size();

if (comm_rank() == 0) {
    // Partition graph
    metis_partition(tGraph, tNumProcs, true);

    // Count vertices per rank
    Vector<index_t> tCounts(tNumProcs, 0);
    for (Vertex* v : tGraph) {
        ++tCounts(v->owner());
    }

    // Distribute vertices to ranks
    for (proc_t p = 0; p < tNumProcs; ++p) {
        Cell<id_t> tLocalIDs(tCounts(p));
        index_t tCount = 0;

        for (Vertex* v : tGraph) {
            if (v->owner() == p) {
                tLocalIDs(tCount++) = v->id();
            }
        }

        // Send to rank p (MPI communication)
        if (p == 0) {
            mLocalNodes = tLocalIDs;
        } else {
            send(tLocalIDs, p, TAG_NODE_DISTRIBUTION);
        }
    }
} else {
    // Receive local node IDs
    receive(mLocalNodes, 0, TAG_NODE_DISTRIBUTION);
}
```

---

### Example 4: Connected Component Analysis

```cpp
#include "fn_Graph_find_connected_partitions.hpp"

Graph tGraph = /* mesh connectivity */;

index_t tMainDomainSize = find_connected_partitions(tGraph);

// Identify isolated vertices/elements
Cell<Vertex*> tIsolatedVertices;

for (Vertex* v : tGraph) {
    if (v->owner() != 0) {  // Not in main component
        tIsolatedVertices.push(v);
    }
}

if (tIsolatedVertices.size() > 0) {
    message(InfoLevel::Warning,
            "Found %u vertices in %d disconnected components",
            (uint) tIsolatedVertices.size(),
            (int) (/* max owner ID */ + 1));

    // Option 1: Remove isolated vertices
    // Option 2: Add connectivity constraints
    // Option 3: Warn user about mesh quality
}
```

---

### Example 5: Custom Graph Algorithm

```cpp
#include "cl_Graph_Vertex.hpp"
#include "cl_Queue.hpp"
#include "cl_DynamicBitset.hpp"

// Example: Find graph diameter (longest shortest path)
index_t find_graph_diameter(Graph& aGraph) {
    index_t tDiameter = 0;

    for (Vertex* tStart : aGraph) {
        // BFS from each vertex
        DynamicBitset tVisited(aGraph.size());
        Queue<Vertex*> tQueue;

        tStart->set_level(0);
        tQueue.push(tStart);
        tVisited.set(tStart->index());

        index_t tMaxDist = 0;

        while (!tQueue.empty()) {
            Vertex* v = tQueue.pop();
            tMaxDist = std::max(tMaxDist, v->level());

            for (uint k = 0; k < v->number_of_vertices(); ++k) {
                Vertex* neighbor = v->vertex(k);
                if (!tVisited.test(neighbor->index())) {
                    neighbor->set_level(v->level() + 1);
                    tQueue.push(neighbor);
                    tVisited.set(neighbor->index());
                }
            }
        }

        tDiameter = std::max(tDiameter, tMaxDist);
    }

    return tDiameter;
}
```

---

## Thread Safety and MPI

### Thread Safety

| Component | Thread Safe? | Notes |
|-----------|--------------|-------|
| **Vertex class** | ❌ No | Adjacency list modification not atomic |
| **BFS/DFS** | ❌ No | Uses vertex flags and level (shared state) |
| **RCM** | ❌ No | Modifies vertex indices in-place |
| **METIS** | ⚠️ Partial | METIS_NodeND thread-safe, PartGraph uses OpenMP internally |
| **Graph sorting** | ❌ No | Reorders Cell in-place |
| **CSR building** | ✅ Yes | If each thread builds separate graph |

**General Rule:** Graph algorithms are **not thread-safe** because they modify vertex properties (`level`, `index`, `owner`, `flag`) as side effects.

**Safe Patterns:**

```cpp
// ✅ GOOD: Each thread works on separate graph
#pragma omp parallel
{
    Graph tLocalGraph = /* thread-local construction */;
    symrcm(tLocalGraph);
}

// ❌ BAD: Multiple threads modifying shared graph
#pragma omp parallel for
for (index_t i = 0; i < tGraph.size(); ++i) {
    bfs(tGraph, tGraph(i));  // Race condition on level/flag
}
```

### MPI Awareness

| Function | MPI-Aware? | Notes |
|----------|------------|-------|
| **bfs/dfs** | ❌ No | Operates on local graphs only |
| **symrcm** | ❌ No | Local reordering (no communication) |
| **metis_nd** | ❌ No | Serial METIS on local graph |
| **ParMETIS** | ✅ Yes | Parallel nested dissection (ordering) |
| **PT-SCOTCH** | ✅ Yes | Parallel nested-dissection *ordering* (`ptscotch_nd`) |
| **build_pargraph_adjacency** | ✅ Yes | Requires `comm_size()`, `comm_rank()` |

**Key Points:**
- Most graph algorithms assume **replicated graph** (each rank has full graph or its portion)
- **ParMETIS/PT-SCOTCH** wrappers take the full graph on rank 0 and distribute CSR slices themselves — they do not accept an already-distributed graph
- Use **gather/scatter** patterns when mixing local and distributed algorithms

**Example: Distributed Partitioning**

```cpp
// Rank 0 holds the complete graph; every other rank passes an empty one.
Graph tGraph;

if ( comm_rank() == 0 )
{
    tGraph = /* full graph, e.g. sparse::create_graph_from_matrix */;

    // every vertex needs an owner before the collective call; for a graph
    // nobody has partitioned, contiguous blocks of graph order are enough
    graph::block_distribution( tGraph, comm_size() );
}

// collective — OUTSIDE the rank guard. Root builds and scatters the per-rank
// CSR slices, every rank calls ParMETIS_V3_NodeND, root gathers the
// permutation and applies it: afterwards index() on root is the new order.
graph::parmetis_nd( tGraph );      // or graph::ptscotch_nd( tGraph )
```

This is the shape `DistMatrix::create_matrix` uses (`cl_SolverDistMatrix.cpp`); the one-line
edit that puts `parmetis_nd()` inside the rank guard is a hang, because the non-root ranks never
enter the barrier.

---

## Common Pitfalls

### 1. Forgetting to Set Continuous Indices

**Problem:** Many algorithms assume `vertex->index()` is in range `[0, N)` and continuous.

❌ **Bad:**

```cpp
Graph tGraph = /* vertices with arbitrary IDs */;
symrcm(tGraph);  // May fail if indices not set
```

✅ **Good:**

```cpp
Graph tGraph = /* construct graph */;

// Ensure continuous indices before algorithms
index_t tCount = 0;
for (Vertex* v : tGraph) {
    v->set_index(tCount++);
}

symrcm(tGraph);  // Now safe
```

**Why:** Bitsets, arrays indexed by `vertex->index()` require dense `[0, N)` range.

---

### 2. Memory Leaks from Adjacency Lists

**Problem:** Forgetting to call `reset_vertex_container()` or `delete vertex`.

❌ **Bad:**

```cpp
Graph tGraph(N, nullptr);
for (index_t i = 0; i < N; ++i) {
    tGraph(i) = new Vertex();
    tGraph(i)->init_vertex_container(degree);
    // ... populate ...
}
// Graph goes out of scope → vertices leaked
```

✅ **Good:**

```cpp
{
    Graph tGraph(N, nullptr);
    // ... use graph ...

    // Cleanup
    for (Vertex* v : tGraph) {
        delete v;  // ~Vertex() calls reset_vertex_container()
    }
    tGraph.clear();
}

// Or use helper:
clear(tGraph);
```

---

### 3. Mixing ID and Index

**Problem:** Confusing `vertex->id()` (permanent identifier) with `vertex->index()` (reordering position).

❌ **Bad:**

```cpp
symrcm(tGraph);

for (Vertex* v : tGraph) {
    // BUG: id() may not match index() after reordering
    aMatrix(v->id(), v->id()) = value;
}
```

✅ **Good:**

```cpp
symrcm(tGraph);

// Create mapping: original ID → new index
Cell<index_t> tIDtoIndex(maxID + 1, gNoIndex);
for (Vertex* v : tGraph) {
    tIDtoIndex(v->id()) = v->index();
}

// Use mapping for matrix assembly
for (Element* e : mesh.elements()) {
    for (Node* n : e->nodes()) {
        index_t tNewIndex = tIDtoIndex(n->id());
        aMatrix(tNewIndex, tNewIndex) += ...;
    }
}
```

---

### 4. Assuming Connected Graphs

**Problem:** Mesh graphs may have disconnected components (isolated elements, floating nodes).

❌ **Bad:**

```cpp
Vertex* tStart = find_pseudo_peripheral_vertex(tGraph);
symrcm(tGraph, tStart);  // Only reorders component containing tStart
```

✅ **Good:**

```cpp
// Check connectivity first
index_t tMainSize = find_connected_partitions(tGraph);

if (tMainSize < tGraph.size()) {
    message(InfoLevel::Warning,
            "Mesh has disconnected components: %u / %u vertices in main domain",
            (uint) tMainSize, (uint) tGraph.size());
}

// RCM handles disconnected graphs automatically
symrcm(tGraph);  // Processes all components
```

---

### 5. Using METIS Without Checking Availability

**Problem:** Calling METIS functions when library not linked.

❌ **Bad:**

```cpp
metis_nd(tGraph);  // Crashes if BELFEM_METIS not defined
```

✅ **Good:**

```cpp
#ifdef BELFEM_METIS
    metis_nd(tGraph);
#else
    message(InfoLevel::Warning, "METIS unavailable, using RCM fallback");
    symrcm(tGraph);
#endif
```

**Or use error handling** — only where a failed `BELFEM_ERROR` throws, i.e. in debug builds or
after `set_throw_on_error( true )`; a release build aborts before the `catch`. Prefer the
`#ifdef BELFEM_METIS` guard.

```cpp
try {
    metis_nd(tGraph);  // BELFEM_ERROR if unavailable: throws in debug, aborts in release
} catch (...) {
    message(InfoLevel::Warning, "METIS failed, falling back to RCM");
    symrcm(tGraph);
}
```

---

### 6. Incorrect CSR Index Base

**Problem:** METIS/SCOTCH may use 0-based or 1-based indexing depending on options.

❌ **Bad:**

```cpp
build_graph_adjacency(tGraph, tVertices, tEdges);
// Assumes METIS_OPTION_NUMBERING = 0, but may be wrong
```

✅ **Good:**

```cpp
build_graph_adjacency(tGraph, tVertices, tEdges);

// Explicitly set METIS options
Cell<metis_t> tOptions(METIS_NOPTIONS, 0);
METIS_SetDefaultOptions(tOptions.data());
tOptions(METIS_OPTION_NUMBERING) = 0;  // 0-based (match CSR builder)

METIS_NodeND(&tNumVertices, tVertices.data(), tEdges.data(),
             nullptr, tOptions.data(), tPerm.data(), tIPerm.data());
```

---

### 7. Ignoring Vertex Degree in RCM

**Problem:** Starting RCM from high-degree vertex produces poor bandwidth.

❌ **Bad:**

```cpp
Vertex* tStart = tGraph(0);  // May be high-degree hub
symrcm(tGraph, tStart);
```

✅ **Good:**

```cpp
// Let symrcm find pseudo-peripheral vertex automatically
symrcm(tGraph);  // Seeds the pseudo-peripheral search heuristically
```

**Or manually:**

```cpp
// Find minimum degree vertex as start candidate
Vertex* tStart = tGraph(0);
uint tMinDegree = tStart->number_of_vertices();

for (Vertex* v : tGraph) {
    if (v->number_of_vertices() < tMinDegree) {
        tMinDegree = v->number_of_vertices();
        tStart = v;
    }
}

// Then find pseudo-peripheral from this candidate
tStart = find_pseudo_peripheral_vertex(tGraph, tStart);
symrcm(tGraph, tStart);
```

---

## Literature References

### Classic Graph Algorithms

**Cuthill-McKee Ordering:**
- **Cuthill, E. & McKee, J.** (1969): "Reducing the Bandwidth of Sparse Symmetric Matrices", *Proceedings of 24th National Conference, ACM*, pp. 157-172.
- **Liu, W.-H. & Sherman, A. H.** (1976): "Comparative Analysis of the Cuthill-McKee and the Reversed Cuthill-McKee Ordering Algorithms for Sparse Matrices", *SIAM Journal on Numerical Analysis*, Vol. 13, pp. 198-213.
- **Gibbs, N. E., Poole, W. G., & Stockmeyer, P. K.** (1976): "An Algorithm for Reducing the Bandwidth and Profile of a Sparse Matrix", *SIAM Journal on Numerical Analysis*, Vol. 13, pp. 236-250.

**Referenced in BELFEM literature:**
- **Bathe** (FEM textbook): Bandwidth reduction discussion (Section 8.2.3)
- **Hughes** (FEM textbook): RCM references in sparse matrix algorithms

### External Partitioning Libraries

**METIS:**
- **Karypis, G. & Kumar, V.** (1998): "A Fast and High Quality Multilevel Scheme for Partitioning Irregular Graphs", *SIAM Journal on Scientific Computing*, Vol. 20, No. 1, pp. 359-392.
- **METIS Manual**: http://glaros.dtc.umn.edu/gkhome/metis/metis/overview

**ParMETIS:**
- **Karypis, G., Schloegel, K., & Kumar, V.** (2003): "ParMETIS: Parallel Graph Partitioning and Sparse Matrix Ordering Library", University of Minnesota.

**SCOTCH:**
- **Pellegrini, F. & Roman, J.** (1996): "SCOTCH: A Software Package for Static Mapping by Dual Recursive Bipartitioning of Process and Architecture Graphs", *HPCN Europe*, LNCS 1067, pp. 493-498.
- **PT-SCOTCH**: https://www.labri.fr/perso/pelegrin/scotch/

### General Graph Theory

- **Cormen, T. H., Leiserson, C. E., Rivest, R. L., & Stein, C.** (2009): *Introduction to Algorithms* (3rd ed.), MIT Press. [BFS/DFS algorithms]

---

## Performance Considerations

### Algorithm Selection Guide

| Objective | Recommended Algorithm | Fallback | Notes |
|-----------|----------------------|----------|-------|
| **Bandwidth reduction** | `symrcm()` | Manual ordering | Fast, works well for banded solvers |
| **Fill-in minimization** | `metis_nd()` | `symrcm()` | Best for direct solvers (MUMPS, STRUMPACK) |
| **Parallel factorization** | `metis_ndp()` | `metis_nd()` | Use `aNumPartitions = comm_size()` |
| **Domain decomposition** | `metis_partition()` | `scotch_partition()` | Balances load, minimizes edge cuts |
| **Parallel nested-dissection ordering** | `parmetis_nd()` / `ptscotch_nd()` | `metis_nd()` | Full graph on rank 0; the wrapper distributes the CSR slices |
| **Connected components** | `find_connected_partitions()` | `dfs()` | Also sorts by size |

### Time Complexity Summary

| Algorithm | Time | Space | Notes |
|-----------|------|-------|-------|
| BFS | O(V + E) | O(V) | Level-synchronous |
| DFS | O(V + E) | O(V) | Stack-based |
| RCM | O(V + E + V log V) | O(V) | V log V from the final whole-graph sort; neighbor sorts are O(V log Δ) |
| Pseudo-peripheral | O(k(V + E)) | O(V) | k ≈ 3 iterations typical |
| METIS NodeND | O(E) average | O(V + E) | Highly optimized (coarsening) |
| METIS Partition | O(E) average | O(V + E) | K-way refinement |

### Memory Usage

**Vertex Adjacency:**
- Per vertex: `sizeof(Vertex*) * degree` bytes for adjacency list
- Typical mesh: average degree ≈ 6-8 for 3D elements
- Example: 1M vertices, degree 8 → 8MB for adjacency lists (64-bit pointers)

**CSR Adjacency:**
- `(V + 1) * sizeof(metis_t)` for vertex array
- `E * sizeof(metis_t)` for edge array
- Example: 1M vertices, 8M edges → (1M + 1) * 4 + 8M * 4 ≈ 36 MB (32-bit `metis_t`)

**Algorithm Temporaries:**
- BFS/DFS: `V/8` bytes for bitset + queue overhead
- RCM: `V * sizeof(index_t)` for permutation array

### Optimization Tips

1. **Minimize graph construction overhead:**
   - Use two-pass pattern (count, then allocate) to avoid reallocation
   - Pre-sort vertices by index before CSR building

2. **Cache efficiency:**
   - RCM improves locality for subsequent matrix operations
   - Sort vertices by index after reordering

3. **Parallel partitioning:**
   - Use ParMETIS for > 5M elements on distributed memory
   - Use METIS with OpenMP for < 5M elements on shared memory

4. **Reuse graphs:**
   - If mesh topology unchanged, reuse partitioning across time steps
   - Store permutation separately rather than rebuilding graph

---

**End of Graph Module Documentation**

---

## See Also

- **Project README**: `../../../README.md`
- **Claude Instructions**: `../../../CLAUDE.md`
- **Documentation Guidelines**: `../../../doc/documentation_guidelines.md`
- **Coding Philosophy**: `../../../doc/coding_philosophy.md`
- **General Documentation**: `../../../doc/README.md`
- **Containers Module**: `../../../src/containers/doc/` (Cell, Queue, Bitset)
- **Core Module**: `../../../src/core/doc/` (typedefs, Logger, Timer)
