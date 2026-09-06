# Graph Module Test Suite 

**Date:** 2026-03-23

Source bugs found and fixed during test development are documented in devlog/dl20260324_test_suite_bugs.md.

---

## Test Count by Suite

| Suite | Tests | Coverage |
|-------|-------|----------|
| `GraphVertex` | 5 | Construction, ID, index, owner, level |
| `GraphVertexFlag` | 5 | Default, flag/unflag, multiple, all 8 flags |
| `GraphVertexDebug` | 4 | Flag/unflag/is_flagged OOB, vertex access OOB |
| `GraphVertexAdj` | 5 | Count/alloc/fill pattern, explicit size, reset, re-init, empty |
| `GraphVertexSort` | 3 | Sort by index, reverse, reverse empty |
| `GraphVertexElement` | 2 | init/reset element container (BELFEM_ERROR — always active) |
| `GraphBfs` | 9 | Path5 (end/middle), Star5 (center/leaf), K4, BinaryTree, single, empty, disconnected |
| `GraphDfs` | 6 | Connected, disconnected, K4, single, empty, all flagged, levels set |
| `GraphPartitions` | 4 | One component, two components, three isolated, empty |
| `GraphPseudoPeripheral` | 6 | Path5 (vertex), BinaryTree, K4, single, empty, Path5 (node endpoint check) |
| `GraphSymrcm` | 7 | Path5 (bandwidth + adjacency), Star5, BinaryTree, disconnected, single, empty, permutation validity |
| `GraphMemory` | 1 | clear() deletes vertices |
| `GraphTools` | 18 | Adjacency (path5, K4, empty, single, monotone, degree), self-loop drop (serial + parallel CSR, loops kept on the Graph, all-loop placeholder, `metis_ndp` on a looped graph under `BELFEM_METIS`), parallel CSR owner sentinels (unowned → root, `comm_size()` marker → root, `-1` and `comm_size()+1` throw), sort regression, permutation (identity, reverse+adjacency) |
| `GraphOperators` | 5 | opVertexDegree, ID, Index, Level, Owner |
| `GraphMetis` | 2 | METIS ND, METIS partition (conditional on `BELFEM_METIS`) |
| `GraphSpfa` | 16 | Difference-constraint feasibility: degenerate/structural (no arcs, isolated vertices, chain, tight unit chain, components), self loops (non-negative feasible, negative infeasible), cycles (zero-weight and balanced **feasible**, negative triangle, greedy-hang c=2, buried cycle), randomized (200 feasible-by-construction, 100 with injected negative cycle), stress (2000-vertex chain, 500-vertex negative return arc) |
| **Total** | **88** | |

### `GraphSpfa` — what these tests actually assert

The suite verifies the *answer*, not the return value. A feasible verdict is
checked by testing **every** constraint against the returned theta, so a solver
returning `true` with a wrong solution cannot pass. An infeasible verdict is
checked by verifying the returned certificate is a genuine closed walk whose
weights sum **strictly** negative.

Two cases are regression guards for defects found during development, and are
worth keeping for that reason alone:

- **Zero-weight and balanced cycles must be FEASIBLE.** They only force the
  potentials around the cycle to be equal. An early solver accepted a zero-weight
  predecessor cycle as proof of infeasibility.
- **Large feasible random instances.** The classic SPFA update-count trigger
  (`cnt > n`) can fire *before* the predecessor graph contains a cycle; an early
  version treated that spurious trigger as an internal error and aborted on a
  feasible system.

The randomized instances are feasible by construction — potentials are drawn
first and every arc weight is set no smaller than the difference it must admit —
so any infeasible verdict there is a false negative, not a fixture bug. The PRNG
is a fixed-seed xorshift so a failure is reproducible.

---

## @warning API Gotchas

- `Graph` is `typedef Cell<graph::Vertex*>` — the graph owns vertex pointers allocated with `new`. Every test MUST call `graph::clear(tGraph)` at the end.
- `bfs()` return value is **max width** (largest level size), NOT max depth. To get eccentricity, scan vertex levels.
- `bfs()` reindexes vertices at start: `tVertex->set_index(tCount++)`. After BFS, `tGraph(i)->index() == i`.
- `find_connected_partitions()` returns **size of largest partition**, NOT the number of partitions.
- `find_pseudo_peripheral_node()` picks **max degree** among farthest-level candidates. `find_pseudo_peripheral_vertex()` picks **min degree** (works correctly, different implementation).
- `graph::sort(Graph&)` is the no-argument overload (was buggy, now fixed). `sort(Cell<T>, Comp)` is the two-argument Cell sort (always worked).
- Element container is DISABLED: `init_element_container()` fires `BELFEM_ERROR` (always active, NOT debug-only).
- Flag bounds checks use `BELFEM_ASSERT` (debug only).

---
