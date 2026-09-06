# BELFEM Graph Tests — Detailed Plan

**Date:** 2026-03-22
**Purpose:** Method-level test matrix for the graph module (`src/math/graph/`)
**Depends on:** `tests_0_strategy.md` (conventions), `tests_1_containers.md` (Cell, DynamicBitset, Queue)
**Confidence:** High on Vertex class and core algorithms. Medium on external library wrappers (METIS/SCOTCH — conditional compilation).

---

## Module Overview

A graph algorithm library for mesh reordering and partitioning, built around `graph::Vertex` with manual adjacency list management.

| Layer | Content |
|---|---|
| `Vertex` class | Node with ID, index, owner, level, 8-bit flags, manual `malloc`/`free` adjacency buffer |
| Core algorithms | `bfs`, `dfs`, `find_connected_partitions`, `find_pseudo_peripheral_node`, `find_pseudo_peripheral_vertex`, `symrcm`, `clear`, `sort` |
| External wrappers | `metis_nd`, `metis_ndp`, `metis_partition` (METIS); ParMETIS; SCOTCH; PT-SCOTCH |
| Utilities | `build_graph_adjacency` (CSR export), `build_pargraph_adjacency`, `apply_graph_permutation`, comparison operators (`opVertexDegree`, `opVertexID`, `opVertexIndex`, `opVertexLevel`, `opVertexOwner`) |

**Key type:** `Graph` is a typedef for `Cell<graph::Vertex*>`. The graph owns vertex pointers; `graph::clear()` deletes them and clears the Cell.

---

## Known Regression Target

| ID | Location | Issue | Regression Test |
|---|---|---|---|
| BUG-G1 | `fn_Graph_sort.cpp` line 2430–2432 | `graph::sort()` passes `.begin()` as both first and last iterator to `std::sort`, making it a no-op. Should be `.begin()` and `.end()`. | Test: call `graph::sort()` on a shuffled graph → verify indices are *not* sorted (confirms bug), or verify they *are* sorted (confirms fix). |

---

## Test Topology Helpers

Claude Code should implement a test-local graph builder function. All test topologies use this pattern:

```cpp
// Build a graph from an adjacency list.
// aAdjacency[i] is a list of neighbor indices for vertex i.
Graph build_test_graph( const Cell< Cell< index_t > > & aAdjacency )
{
    index_t tN = aAdjacency.size();
    Graph tGraph( tN, nullptr );

    // create vertices
    for( index_t i = 0; i < tN; ++i )
    {
        tGraph( i ) = new graph::Vertex();
        tGraph( i )->set_id( i );
        tGraph( i )->set_index( i );
    }

    // count neighbors
    for( index_t i = 0; i < tN; ++i )
    {
        for( index_t j = 0; j < aAdjacency( i ).size(); ++j )
        {
            tGraph( i )->increment_vertex_counter();
        }
    }

    // allocate and fill
    for( index_t i = 0; i < tN; ++i )
    {
        tGraph( i )->init_vertex_container();
        for( index_t j = 0; j < aAdjacency( i ).size(); ++j )
        {
            tGraph( i )->insert_vertex( tGraph( aAdjacency( i )( j ) ) );
        }
    }

    return tGraph;
}
```

### Standard Test Topologies

| Name | Vertices | Edges | Properties |
|---|---|---|---|
| **Path₅** | 5 | 0—1—2—3—4 | Bandwidth 1, endpoints are pseudo-peripheral |
| **Cycle₅** | 5 | 0—1—2—3—4—0 | All degree 2, BFS depth depends on start |
| **Star₅** | 5 | 0—{1,2,3,4} | Center degree 4, leaves degree 1, BFS depth 1 from center |
| **K₄** | 4 | Complete graph | All degree 3, BFS depth 1 from any vertex |
| **BinaryTree₇** | 7 | 0→{1,2}, 1→{3,4}, 2→{5,6} (symmetric edges) | Depth 2, pseudo-peripheral = any leaf |
| **Disconnected** | 6 | Triangle {0,1,2} + Triangle {3,4,5} | 2 connected components |
| **SingleVertex** | 1 | None | Edge case |
| **EmptyGraph** | 0 | None | Edge case |

All adjacency lists should be symmetric (undirected): if vertex i lists j as a neighbor, vertex j must list i.

---

## Test File Structure

```
tests/graph/
├── test_GraphVertex.cpp          # Vertex class semantics
├── test_GraphAlgorithms.cpp      # BFS, DFS, partitions, pseudo-peripheral, symrcm
├── test_GraphTools.cpp           # CSR export, permutation, sort operators
```

---

## 1. Vertex Class

**File:** `test_GraphVertex.cpp`

### 1.1 Construction & Properties `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `DefaultConstruction` | `id() == gNoID`, `index() == gNoIndex`, `owner() == gNoID`, `level() == 0` |
| `SetGetId` | `set_id(42)` → `id() == 42` |
| `SetGetIndex` | `set_index(7)` → `index() == 7` |
| `SetGetOwner` | `set_owner(3)` → `owner() == 3` |
| `SetGetLevel` | `set_level(5)` → `level() == 5` |

### 1.2 Flag System `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `FlagDefault` | `is_flagged(0)` through `is_flagged(7)` all return false initially |
| `FlagAndTest` | `flag(0)` → `is_flagged(0) == true`, other flags still false |
| `UnflagAndTest` | `flag(3)`, `unflag(3)` → `is_flagged(3) == false` |
| `MultipleFlagsIndependent` | `flag(0)`, `flag(5)` → both true, others false |
| `AllEightFlags` | Flag and test all indices 0–7 |

### 1.3 Flag System `[debug]`

| Test Name | What It Verifies |
|---|---|
| `FlagIndexOutOfBoundsThrows` | `flag(8)` → assertion |
| `UnflagIndexOutOfBoundsThrows` | `unflag(8)` → assertion |
| `IsFlaggedIndexOutOfBoundsThrows` | `is_flagged(8)` → assertion |

### 1.4 Adjacency Container Lifecycle `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CountAllocFillPattern` | `increment_vertex_counter()` × 3, `init_vertex_container()`, `insert_vertex()` × 3 → `number_of_vertices() == 3`, all pointers accessible |
| `InitWithExplicitSize` | `init_vertex_container(5)` → can insert up to 5 vertices |
| `ResetClearsContainer` | After populating, `reset_vertex_container()` → `number_of_vertices() == 0` |
| `ReInitAfterFill` | Fill with 3, then re-init with 2 → works correctly (old buffer freed) |
| `EmptyContainer` | `init_vertex_container()` with counter=0 → `number_of_vertices() == 0`, no crash |

### 1.5 Adjacency Access `[debug]`

| Test Name | What It Verifies |
|---|---|
| `VertexAccessOutOfBoundsThrows` | `vertex(3)` when only 2 neighbors → assertion |

### 1.6 Sort and Reverse `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `SortVerticesByIndex` | Insert neighbors with indices {5, 1, 3}; `sort_vertices()` → order is {1, 3, 5} |
| `ReverseVertices` | Insert {A, B, C}; `reverse_vertices()` → order is {C, B, A} |
| `ReverseEmpty` | `reverse_vertices()` on empty container → no crash |

### 1.7 Element Container (Disabled) `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `InitElementContainerThrows` | `init_element_container()` → `BELFEM_ERROR` fires (always active) |
| `ResetElementContainerThrows` | `reset_element_container()` → `BELFEM_ERROR` fires |

---

## 2. BFS

**File:** `test_GraphAlgorithms.cpp`

### 2.1 BFS with Start Vertex `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `BfsPath5FromEnd` | Path₅, start=vertex 0 → max width = 1, levels are 0,1,2,3,4 |
| `BfsPath5FromMiddle` | Start=vertex 2 → levels are 2,1,0,1,2; max width = 2 |
| `BfsStar5FromCenter` | Start=center → all leaves at level 1, max width = 4 |
| `BfsStar5FromLeaf` | Start=leaf → center at level 1, other leaves at level 2 |
| `BfsK4` | Complete K₄, any start → all others at level 1, max width = 3 |
| `BfsBinaryTree` | Start=root → levels correct, max width = 4 (leaf level) |
| `BfsSingleVertex` | 1 vertex → max width = 1, level = 0 |
| `BfsEmptyGraph` | 0 vertices → returns 0 |

### 2.2 BFS Disconnected Overload `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `BfsDisconnectedGraph` | Two triangles → all vertices get levels, max width across components |
| `BfsConnectedViaDisconnectedOverload` | Path₅ through disconnected overload → same result as single-start version |

---

## 3. DFS

**File:** `test_GraphAlgorithms.cpp`

### 3.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `DfsConnectedReturnsOne` | Path₅ → returns 1 component, all owners equal |
| `DfsDisconnectedReturnsTwoComponents` | Two triangles → returns 2, vertices in each triangle share owner, owners are distinct |
| `DfsK4ReturnsOne` | K₄ → 1 component |
| `DfsSingleVertex` | 1 vertex → 1 component |
| `DfsEmptyGraph` | 0 vertices → returns 0 |
| `DfsAllVerticesFlagged` | After DFS, every vertex `is_flagged()` returns true |
| `DfsLevelsSet` | All vertices have `level() != gNoIndex` after DFS |

---

## 4. Connected Partitions

**File:** `test_GraphAlgorithms.cpp`

### 4.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `FindConnectedPartitionsOne` | Path₅ → returns 1 |
| `FindConnectedPartitionsTwo` | Two triangles → returns 2 |
| `FindConnectedPartitionsThree` | Three isolated vertices → returns 3 |
| `FindConnectedPartitionsEmpty` | Empty graph → returns 0 |

---

## 5. Pseudo-Peripheral Vertex

**File:** `test_GraphAlgorithms.cpp`

### 5.1 `find_pseudo_peripheral_vertex` `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `PseudoPeripheralPath5` | Path₅ → returns an endpoint (vertex 0 or vertex 4) |
| `PseudoPeripheralBinaryTree` | Returns a leaf (maximum eccentricity) |
| `PseudoPeripheralK4` | All vertices equivalent → returns any vertex, no crash |
| `PseudoPeripheralSingleVertex` | Returns the single vertex |
| `PseudoPeripheralEmptyGraph` | Returns `nullptr` |
| `PseudoPeripheralEccentricity` | Run BFS from returned vertex → level depth ≥ BFS from arbitrary start |

### 5.2 `find_pseudo_peripheral_node` `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `PseudoPeripheralNodePath5` | Returns an endpoint |
| `PseudoPeripheralNodeMatchesVertex` | Both functions return vertices with similar eccentricity on the same graph |

---

## 6. SymRCM (Reverse Cuthill-McKee)

**File:** `test_GraphAlgorithms.cpp`

### 6.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `SymrcmPath5` | Path₅ → bandwidth does not increase (already optimal, bandwidth = 1) |
| `SymrcmBinaryTree` | BinaryTree₇ → bandwidth after RCM ≤ bandwidth before |
| `SymrcmStar5` | Star₅ → reordering reduces or maintains bandwidth |
| `SymrcmDisconnected` | Two triangles → handles both components, no crash |
| `SymrcmSingleVertex` | No crash |
| `SymrcmEmptyGraph` | No crash |
| `SymrcmPermutationIsValid` | After symrcm, every index 0..N-1 appears exactly once |
| `SymrcmPreservesAdjacency` | After symrcm, each vertex's neighbors are the same set of vertices (by ID) as before |

### 6.2 Bandwidth Helper

Claude Code should implement a test-local bandwidth calculator:

```cpp
index_t compute_bandwidth( const Graph & aGraph )
{
    index_t tBandwidth = 0;
    for( index_t i = 0; i < aGraph.size(); ++i )
    {
        for( uint k = 0; k < aGraph( i )->number_of_vertices(); ++k )
        {
            index_t tDiff = std::abs(
                (long) aGraph( i )->index() - (long) aGraph( i )->vertex( k )->index() );
            tBandwidth = std::max( tBandwidth, tDiff );
        }
    }
    return tBandwidth;
}
```

---

## 7. Graph Utilities

**File:** `test_GraphTools.cpp`

### 7.1 CSR Export `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `BuildAdjacencyPath5` | Path₅ → `aVertices` has 6 entries, total edges = 8 (4 edges × 2 directions), CSR structure is valid |
| `BuildAdjacencyK4` | K₄ → 12 edges (6 × 2), every vertex has 3 neighbors |
| `BuildAdjacencyEmpty` | Empty graph → `aVertices` has 1 entry (value 0) |
| `BuildAdjacencySingleVertex` | 1 vertex, 0 edges → `aVertices` = {0, 0} |
| `BuildAdjacencyVerticesMonotone` | `aVertices(i) <= aVertices(i+1)` for all i |
| `BuildAdjacencyDegreeConsistent` | `aVertices(i+1) - aVertices(i) == vertex(i)->number_of_vertices()` |

### 7.2 Graph Sort (BUG-G1) `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `GraphSortRegressionBugG1` | Create a graph with shuffled indices, call `graph::sort()`. If the bug is present (`.begin(), .begin()`), graph is NOT sorted. If fixed, graph IS sorted. Document either outcome. |

### 7.3 Apply Graph Permutation `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `ApplyPermutationIdentity` | Identity permutation → graph unchanged |
| `ApplyPermutationReverse` | Reverse permutation → indices flipped, adjacency preserved |
| `ApplyPermutationPreservesAdjacency` | After permutation, each vertex's neighbors (by ID) are the same set |

### 7.4 Comparison Operators `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `OpVertexDegree` | Sorts by `number_of_vertices()` ascending |
| `OpVertexID` | Sorts by `id()` ascending |
| `OpVertexIndex` | Sorts by `index()` ascending |
| `OpVertexLevel` | Sorts by `level()` ascending |
| `OpVertexOwner` | Sorts by `owner()` ascending |

---

## 8. External Library Wrappers (Conditional)

### 8.1 METIS (only if `BELFEM_METIS` defined) `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `MetisNdPath5` | `metis_nd()` produces a valid permutation (all indices unique) |
| `MetisNdBinaryTree` | Valid permutation, adjacency preserved |
| `MetisPartition2` | `metis_partition(graph, 2)` → all owners in {0, 1} |
| `MetisPartition4K4` | K₄ partitioned into 4 → each vertex gets distinct owner |

### 8.2 Notes

METIS, ParMETIS, SCOTCH, and PT-SCOTCH tests should be wrapped in `#ifdef BELFEM_METIS` / `#ifdef BELFEM_SCOTCH` etc. ParMETIS and PT-SCOTCH additionally require MPI (Tier 2 tests). Smoke-test level only — verify valid output, not optimality.

---

## 9. Memory Management

### 9.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `GraphClearDeletesVertices` | `graph::clear(tGraph)` → `tGraph.size() == 0`, no memory leak (Valgrind target) |
| `VertexDestructorFreesAdjacency` | Create vertex, fill adjacency, delete → no leak (Valgrind target) |

---

## 10. Implementation Notes for Claude Code

1. **Use the graph builder helper** for all tests. Never manually `malloc` adjacency buffers in individual tests — use the `build_test_graph` pattern.
2. **Always call `graph::clear()` at the end of each test** that creates vertices with `new`. Otherwise every test leaks. Alternatively, use a test fixture with cleanup in `TearDown()`.
3. **BUG-G1 is a real bug.** The `graph::sort()` function passes `.begin()` twice to `std::sort`. The test should document this: either it confirms the bug exists (sort is no-op) or it confirms the fix works.
4. **Adjacency must be symmetric.** When building test graphs for undirected algorithms, ensure that if vertex i lists j, vertex j lists i. The `build_test_graph` helper takes care of this if the input adjacency list is symmetric.
5. **BFS returns max width, not depth.** The return value of `bfs()` is the maximum number of vertices at any single level, not the number of levels. Test accordingly.
6. **DFS uses vertex flags.** After `dfs()`, all reachable vertices are flagged. Tests that call DFS on a graph and then call other algorithms should clear flags first.
7. **`find_pseudo_peripheral_node` vs `find_pseudo_peripheral_vertex`:** Two different implementations with different heuristics. The `_node` version picks max-degree from farthest level; the `_vertex` version picks min-degree. Both should return high-eccentricity vertices, but not necessarily the same one.
8. **Symrcm modifies the graph in-place.** After `symrcm()`, vertex indices are changed and the graph Cell is reordered. Tests that need the original graph for comparison should save adjacency (by vertex ID) before calling symrcm.
9. **External library tests are smoke tests.** Don't test METIS/SCOTCH optimality — just verify valid output (unique indices, valid owners, no crashes).

---

## 11. Codex Audit Checklist

When reviewing Claude Code's test implementation, verify:

- [ ] Graph builder helper exists and produces symmetric adjacency
- [ ] Every test that creates vertices with `new` calls `graph::clear()` or has fixture cleanup
- [ ] BFS tests check max width (return value) AND vertex levels
- [ ] DFS tests verify component count, owner consistency, and flag state
- [ ] Symrcm tests compute bandwidth before and after, verify bandwidth ≤ original
- [ ] Symrcm tests verify permutation validity (all indices 0..N-1 appear exactly once)
- [ ] BUG-G1 test is present and documents the `.begin(), .begin()` issue
- [ ] CSR export tests verify `aVertices` is monotone and consistent with vertex degrees
- [ ] External library tests wrapped in appropriate `#ifdef`
- [ ] Memory tests are `[valgrind]` targets
- [ ] BELFEM naming conventions (`t` prefix for locals)
