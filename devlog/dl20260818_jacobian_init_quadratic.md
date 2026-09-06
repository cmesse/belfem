# Jacobian Initialization Is Quadratic in the DOF Count

**Date:** 2026-08-18
**Purpose:** Read-only three-AI investigation into why `DofManager::init_matrices()` takes ~200 s
             for a 2.1e6 DOF problem on 5 ranks
**Module:** `fem/kernel` (DofManager / SolverData), `containers` (DynamicBitset)
**Status:** Investigation only. No source was modified. No fix has been applied or tested.

## Trigger

A production run (2,100,000 DOFs, 5 MPI ranks, Intel i9-10900X, 64 GB) is healthy but spends
~200 s in the phase logged as " Initialize Jacobian ... " / "... time for initializing Jacobian"
(`src/fem/kernel/cl_FEM_DofManager.cpp:399-424`). Question: is there an obvious algorithmic
cause in the DofManager and kernel initialization routines.

Method: Claude pre-registered hypotheses C1-C5 in `tmp/ai_exchange/jacobian_init_perf.md`,
then Codex and Grok audited in parallel blind jury mode with an explicit instruction to
refute C1. Claude additionally built standalone microbenchmarks to replace estimates with
measurements. Everything below is static reading plus those microbenchmarks — **nothing was
measured inside BELFEM itself**, and no fix was built or run.

## Root cause: `DynamicBitset::where()` scans the whole bitset, once per DOF

`SolverData::populate_graph` (`src/fem/kernel/cl_FEM_DofMgr_SolverData.cpp:789`) allocates one
`DynamicBitset` sized to the whole local DOF space (`:810`), marks each row's neighbors, then
extracts them with

```cpp
tBitset.where( tIndices );          // cl_FEM_DofMgr_SolverData.cpp:845
```

`where()` defaults to `aAssumeSparse = true` (`src/containers/cl_DynamicBitset.hpp:206`,
`:559-573`) and dispatches to `where_sparse` (`src/containers/cl_DynamicBitset.cpp:367-399`),
whose loop is

```cpp
for (index_t b = 0; b < mMemorySize; ++b)   // mMemorySize = ceil(N/64)
```

with no early exit and no population shortcut. `where_dense` (`:326-364`) has the same
full-width loop plus a `count()` pass. **Neither function has a sparse path** — the flag only
selects `push()` versus preallocated writes. So one `where()` costs O(N/64) no matter how few
bits are set, and it is called once per row DOF.

On rank 0 the graph is the global one: `reorder_dofs` sets
`mMyNumberOfFreeDofs = mNumberOfFreeDofs` (`src/fem/kernel/cl_FEM_DofMgr_DofData.cpp:3295-3296`).
At N = 2.1e6 the bitset is 32,813 words, so a full pass is ~6.9e10 word inspections.

The author was aware full-width bitset ops are expensive: `tBitset.reset()` is commented out at
`:825` and individual bits are cleared instead at `:872`. `where()` silently puts the full-width
scan back.

### Measurement

Standalone benchmark reproducing the pattern exactly (mark ~30 neighbor bits, `where_sparse`
full-width scan, clear only the listed bits), `c++ -O2 -std=c++17`, Apple silicon dev machine —
**not** the i9-10900X that produced the 200 s:

| N (DOFs) | bitset words | one full-width pass |
|---------:|-------------:|--------------------:|
|  100,000 |        1,563 |   0.116 s |
|  200,000 |        3,125 |   0.393 s |
|  400,000 |        6,250 |   1.527 s |
| 2,100,000|       32,813 |  37.806 s |

3.4x / 3.9x per doubling — clean O(N^2). The 262 kB bitset is L2-resident on the i9 (1 MB
L2/core), so this is instruction-bound, not bandwidth-bound; the earlier DRAM-streaming
rationale in the pre-registration was the wrong machine model even though it landed in the
right decade.

### How many passes actually pay it

`where()` sits behind `if ( tA->is_fixed() == tRowIsFixed || tUseFullMatrix )` (`:823`), so a
call whose row type is "fixed" only scans fixed rows. Grok caught this; the pre-registration's
"5 x 37.8 s = 189 s accounts for the 200 s in full" was wrong and is withdrawn.

| Call site | Rows hitting `where()` | Full-width pass? |
|---|---|---|
| `reorder_dofs` Jacobian (`DofData.cpp:3299`) | all free | yes |
| `reorder_dofs` Imposition (`DofData.cpp:3315`) | fixed only | no, unless n_fixed ~ N |
| `allocate_matrices` FullMass (`SolverData.cpp:282`) | all DOFs | yes, iff `mUseFullForce` |
| `allocate_matrices` Dirichlet (`SolverData.cpp:351`) | all free | yes (bitset nearly empty, scan cost identical) |
| `allocate_matrices` System (`SolverData.cpp:388`) | all free | yes |
| `allocate_matrices` Enforcement / Imposition (`:362`, `:373`) | fixed | not called; `mUseJediForce = false` (`SolverData.hpp:181`) |

`MaxwellFactory.cpp:686` calls `use_full_force( true )` unconditionally for the magnetic field,
so an `hphirun` job pays **four** full-width passes; a non-Maxwell job pays three.

4 x 37.8 s = 151 s of `where()` scanning, plus the hash-lookup cost below, on a machine that is
not the one that reported 200 s. Consistent with the observation; not a closed identity.

## Secondary costs (measured or verified, ranked)

| # | Cost | Site | Magnitude |
|---|---|---|---|
| 2 | `mDofData->dof( id )` is an `std::unordered_map` lookup, called once per DOF and once per adjacency entry | `SolverData.cpp:816`, `:830`; `DofData.hpp:406-410`; `cl_Map.hpp:71` | **measured 5.6 s** for 6.3e7 lookups over a 2.1e6-entry map; ~22 s over four passes |
| 3 | `compute_dof_dof_connectivity` and `unite_dofs` each build `Cell< Vector<id_t> >` of size N and `unique()` every row; under Armadillo `unique()` allocates a fresh vector per row (`fn_unique.hpp:33`) | `SolverData.cpp:1297-1342`, `:1424-1473` | ~2 x 2.1e6 small allocations + 2 x 2.1e6 small sorts, est. 8-20 s |
| 4 | `reset_vertex_container` / `init_vertex_container` free+malloc per DOF per pass, including rows that are then skipped | `SolverData.cpp:817`, `:859-860`; `cl_Graph_Vertex.cpp:49-77` | ~N malloc/free pairs per pass, est. 2-6 s |
| 5 | `create_assembly_tables` ships the System COO twice — `Jacobian` is a child of `System` sharing indices (`SolverData.cpp:405`) but both are in the `tNumMatrices = 6` loop | `SolverData.cpp:605-663` | duplicated GB-scale collect; the `index()` arithmetic itself is **measured 0.26 s** and is not a concern |

`SpMatrix::index()`'s binary search was pre-registered as suspect C4 and is refuted by
measurement: 6.3e7 lookups take 0.26 s.

## Two further instances of the same bug, outside this timer

1. **`Kernel::partition_mesh`** (`src/fem/kernel/cl_FEM_Kernel.cpp:255-301`) builds the
   element-element graph with the identical pattern, and worse: per element it calls both
   `tBitset->reset()` (full-width memset) **and** `tBitset->where()` (full-width scan) over a
   bitset sized to the element count. That is O(2 * N_elem^2 / 64) on rank 0, in the
   "Partitionig mesh" timer rather than the Jacobian one. Codex surfaced this.

2. **`graph::symrcm`'s disconnected-component fallback**
   (`src/math/graph/fn_Graph_symrcm.cpp:139-158`) rescans all remaining vertices for each new
   component: O(C * V). Harmless on the connected free-DOF graph. The **fixed**-DOF graph is
   built with `aLinkToSelf = false` and can be highly fragmented — an isolated Dirichlet DOF is
   its own component — so `symrcm( aFixedDofs )` (`DofData.cpp:3319`) can approach O(n_fixed^2).
   Both auditors flagged it; the pre-registration had waved it off. Whether it fires depends on
   n_fixed and boundary-patch connectivity, neither of which is known for this run.

## Direction for a fix (not implemented, not tested)

Keep the bitset for O(1) membership; delete only the extraction scan. Record indices during
the marking loop via test-and-set into a scratch `Cell`, then sort. Benchmarked on the same
workload and machine: **0.654 s versus 37.806 s, 57.8x**.

Invariants any such rewrite must preserve — the first was verified in-tree, the rest are
Grok's and are unverified:

1. **Sorted CSR columns.** `SpMatrix` is constructed with `aSortGraph = false`
   (`SolverData.cpp:353`, `:396`), `create_csr_indices` copies `vertex(k)->index()` in insertion
   order (`cl_SpMatrix.cpp:546-555`), and `index_csr_zero_based` binary-searches it
   (`cl_SpMatrix.cpp:1515`). Today sortedness is a side effect of `where()`'s ascending scan.
   **Verified safe:** `reorder_dofs` sorts both graphs by index and then assigns `my_index` by
   position (`DofData.cpp:3423-3439`), so ascending `my_index` implies ascending `index` on every
   rank. Sorting the scratch list by `my_index + tOff` preserves it.
2. **Sort key is `my_index + tOff`, not `index()`.** FullMass's `populate_graph` runs *before*
   the fixed-DOF index shift (`SolverData.cpp:282` then `:284-290`).
3. **Uniqueness.** `set()` is not test-and-set today; the bit is what dedups overlapping
   element contributions. A scratch-list path must test before pushing.
4. **`aLinkToSelf = false`.** Self is currently cleared *after* the mark loop (`:842-843`); a
   push-during-mark path must skip or erase it.
5. **Do not switch to `where_dense`** — it adds a full `count()` scan on top.

Secondary fixes in priority order: rewrite the neighbor IDs in `aGraphData` to `my_index` once
so `populate_graph` stops hashing (kills item 2); hoist the per-DOF `Vector` allocations out of
`compute_dof_dof_connectivity` / `unite_dofs`; stop shipping the System COO twice in
`create_assembly_tables`.

## Reconciliation

| Hypothesis | Codex | Grok | Outcome |
|---|---|---|---|
| C1 `where()` full-width scan is the dominant cost | CONFIRMED | mechanism CONFIRMED, "5 passes / accounts for 200 s in full" REFUTED | mechanism stands; pass count corrected 5 -> 4 (Maxwell) or 3; Grok correct |
| C2 unordered_map per adjacency entry | CONFIRMED, line cites corrected | CONFIRMED, line cites corrected | stands; pre-registration cited `:1391`/`:1424`, correct sites are `:1409`/`:1448` |
| C3 rank-0 serialization in `extract_graph_from_mesh` | PARTIALLY CORRECT | CONFIRMED structurally, plus: rank 0 already holds the complete graph, so `unite_dofs` is redundant serial work | stands, second-order; Grok's addition is new |
| C4 `create_assembly_tables` binary searches | CONFIRMED | CONFIRMED as structure, not as seconds | **refuted by measurement** (0.26 s); the real cost there is the duplicated System/Jacobian COO transfer |
| C5 symrcm is fine | PARTIALLY CORRECT | PARTIALLY CORRECT | both auditors right, Claude wrong to wave it off: disconnected fallback is O(C*V) |

Claude's contributions the auditors did not have: the O(N^2) scaling table, the 57.8x
fix benchmark, the 5.6 s hash measurement, the 0.26 s refutation of C4, and the in-tree
verification of the sorted-column invariant.
Grok's decisive correction: the `tRowIsFixed` guard at `:823`.
Codex's decisive contribution: the same bug in `Kernel::partition_mesh`.

## Open questions for the next session

- Which executable produced the 200 s, and is `mUseFullForce` on? Determines 3 versus 4 passes.
- Was the binary built with `USE_DEBUG=ON` (the tree default)? The assertion set inflates all of
  the above and none of the numbers here account for it.
- What is n_fixed / N for this run? Decides whether the `symrcm` fallback matters.
- What does the "Partitionig mesh" timer report on the same run? Predicted to be large.
- **No fix has been written or benchmarked in-tree.** All timings above are from standalone
  microbenchmarks on a different machine.

## Artifacts

- `tmp/ai_exchange/jacobian_init_perf.md` — pre-registration, measured evidence, both audits
  (AI-only tier, ephemeral)
