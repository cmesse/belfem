# Wire ParMETIS / PT-Scotch Parallel Nested Dissection

**Date:** 2026-08-31
**Purpose:** Bring the two parallel nested-dissection entry points into use. They are fully written
and reachable from nothing; wiring them requires fixing a latent out-of-bounds write in the shared
adjacency builder first, because wiring is exactly what makes it reachable.
**Module:** `src/math/graph`
**AIs involved:** Claude (survey + plan)
**Status:** CLOSED 2026-09-04 — substance landed and verified: D1–D5 fixed, `reordering scheme :
parmetis | ptscotch` wired through `DistMatrix::order_graph`, `make check` 19/19 with the new
Tier-2 `sparsempi` suite green at np 2 and np 4 (12/12 each). **Sole remnant is a run gate, W-R5,
the benchmark — tracked as DR-156 in `todo/debt_register.md`, not here.** Also open, not gates:
`compute_permutation` is dead (Christian's call whether to retire it); the new guide prose has had
no Codex language sweep. Six jury rounds in `tmp/ai_exchange/review_parmetis_d1.md`; devlog
`dl20260904_parmetis_d1_owner_guard.md`.

> **Scope guards:**
> - D1 below must land **before** the first caller, not after. It is harmless only while dead.
> - Out of scope: changing the partitioner or the ordering strategy itself.

---

## 1. Current state

| Symbol | Definition | Declaration | Callers |
|---|---|---|---|
| `parmetis_nd` | `fn_Graph_ParMETIS.cpp:31` | `fn_Graph_ParMETIS.hpp:23` | **none** |
| `ptscotch_nd` | `fn_Graph_PTSCOTCH.cpp:28` | `fn_Graph_PTSCOTCH.hpp:23` | **none** |
| `build_pargraph_adjacency` | `graphtools.hpp:126` | — | the two above, and the four R1 tests in `tests/math/test_GraphTools.cpp` |

Verified with `command grep` across `src/`, `tests/`, `archive/` and `nonfree/`. (The `ptscotch_nd`
hits under `tmp/STRUMPACK/` are a vendored tree and a different symbol,
`sep_tree_from_ptscotch_nd_tree`.)

Build gates already exist: `BELFEM_PARMETIS` (`config/linalg/config_metis.cmake:7-9`) and
`BELFEM_PTSCOTCH` (`config/linalg/config_scotch.cmake:6-8`), so the guards a caller needs are in
place.

## 2. Defects that block wiring

- [x] **D1 — CRITICAL once reachable. Out-of-bounds write in `build_pargraph_adjacency`.** *(fixed and verified 2026-09-04; the cites below describe the pre-fix source)*
      `graphtools.hpp:136` sizes the counter to the rank count:
      ```cpp
      Vector< T > tCount( tCommSize, 0 );
      ```
      `:142-146` guards exactly **one** unassigned marker before using the owner as an index:
      ```cpp
      if ( tVertex->owner() == tCommSize )   // the comm_size() marker only
      { tVertex->set_owner( 0 ); }
      ++tCount( tVertex->owner() );          // unguarded for gNoOwner
      ```
      The guard catches the deliberate `owner() == comm_size()` marker set at
      `cl_FEM_Kernel.cpp:195,203`. It does **not** catch `gNoOwner`, which is
      `std::numeric_limits<proc_t>::max()` with `proc_t = int` (`src/core/typedefs.hpp:42,59`) —
      i.e. **2147483647**. A vertex still holding it indexes `tCount` at 2147483647 against a
      length of `comm_size()`.
      Debug builds assert (`Vector` is bounds-checked). **Release builds perform an out-of-bounds
      increment — a write, not a read.**
      *Fix:* remap both sentinels to root, then `BELFEM_ERROR( tOwner >= 0 && tOwner < tCommSize, … )`
      — always-on, because the loop is setup-rate and the release failure mode is a heap write;
      `>= 0` because `proc_t` is signed. See §6.1.

      **Note on failure loudness.** Before the `mOwner` default was corrected from `gNoID` to
      `gNoOwner` (`cl_Graph_Vertex.hpp:42`, fixed 2026-08-31), the stale value was `-1`, which under
      `Vector`'s `size_t` subscript became `2^64-1`, an offset that on this platform faulted in
      practice (the `gNoID` `uint`→`int` conversion is implementation-defined in C++17; the
      `-1`→`size_t` subscript conversion is well-defined `SIZE_MAX`; the out-of-bounds access
      itself guarantees nothing — so read this as "did fault here", not "must fault"). `2^31-1` is a
      plausible offset and so **less** certain to fault. The `mOwner` fix is correct and must not be
      reverted — it is required by the `std::min` ownership sweeps at
      `cl_Mesh_Partitioner.cpp:253-260` and `cl_FEM_Kernel.cpp:427`, which need the sentinel to
      behave as +∞ — but it moved this site from loud to quiet, which is why D1 is a blocker rather
      than a cleanup.

- [x] **D2 — ParMETIS separator-size buffer under-allocated** (Codex, plan jury 2026-09-04;
      adjudicated and fixed 2026-09-04: `tMySizes( 2 * tCommSize, 0 )`). *Pre-fix:* the
      allocation (then at `fn_Graph_ParMETIS.cpp:91`) had `2 * floor( log2( P ) )` entries and was
      passed as the `sizes` output of `ParMETIS_V3_NodeND`. The library documents `sizes` as "the
      2*nparts array" (`KarypisLab/ParMETIS` `libparmetis/ometis.c`); `sizes[0] = 2*npes-1` is a
      cursor and the pre-decremented separator writes reach index `2·npes−2`, the first `npes`
      slots hold subdomain sizes — the API length is `2·npes`. At P = 2 that is 2 entries against 4 required — an
      out-of-bounds write by the library, independent of D1. *Fix:* allocate `2 * tCommSize`.
      Blocks R3, not R1.
- [x] **D3 — PT-Scotch receives a phantom edge on an empty rank** (Codex, plan jury 2026-09-04;
      adjudicated and fixed 2026-09-04: `edgelocnbr` = CSR terminal, `edgelocsiz` = buffer length). *Pre-fix:* `graphtools.hpp:212-216` allocates a
      one-element placeholder edge buffer for an owner with zero vertices (to keep ParMETIS off a
      null pointer) while its CSR terminal stays 0 (`:236-239`); the wrapper (then at
      `fn_Graph_PTSCOTCH.cpp:88`) took `edgelocnbr = tEdges.length()` — 1, not 0 — into
      `SCOTCH_dgraphBuild`.
      *Fix:* use the CSR terminal `tVertices( tLocalNumVerts )` for `edgelocnbr`, keep the
      placeholder as storage only. Blocks R3, not R1.

- [ ] **D4 — ParMETIS refuses a rank with zero vertices, and the builder produces such slices**
      (Codex, code jury round 3, 2026-09-04; confirmed from the library source
      `libparmetis/weird.c`, `CheckInputsNodeND`: `vtxdist[mype+1]-vtxdist[mype] < 1` →
      "Poor initial vertex distribution. Processor %d has no vertices assigned to it!" → error
      status). `graphtools.hpp:214-224` hands any owner with `tCount( p ) == 0` a zero-vertex slice,
      so `parmetis_nd` on such a distribution fails at the `BELFEM_ERROR( tStatus == METIS_OK )`
      after the call — loudly, but rank-locally, which in a debug build strands the peers. The
      production trigger is **`N < P`** under a block split (`tDiv = 0`), not a hand-made
      all-on-root assignment. Not a D2/D3 patch: the remedy is a
      wiring decision tied to O1 — validate the distribution **collectively** before the call (every
      rank ≥ 1 vertex, else all ranks take the same fallback, e.g. the serial `metis_ndp` path) or
      guarantee the distribution upstream. PT-Scotch has no such restriction (D3 makes the empty
      rank legal there). Blocks R3.

- [x] **D5 — self-loops in the matrix graph corrupt METIS_NodeNDP's heap** (found by the first
      real run of `sparsempi`, 2026-09-04; Grok had flagged the self-loops in rounds 4 and 5 for
      ParMETIS — it was plain METIS that broke first). `create_graph_from_matrix` inserts the
      diagonal (`fn_create_graph_from_matrix.cpp:50-57`), so the `DistMatrix` graph carries one loop
      per row; `metis_ndp` on it dies in `FM_2WayNodeRefine1Sided → rpqDestroy → free()` with
      `munmap_chunk(): invalid pointer`, at 2 and at 4 parts. **This is today's production path**
      (PETSc at `comm_size() > 1` with `reordering scheme : metis`), latent because no example deck
      uses that combination. Reproduced with a serial scratchpad probe (64-row Poisson graph: crash
      with loops, clean bijection without); **fixed** by dropping self-loops in
      `build_graph_adjacency` and `build_pargraph_adjacency` — the CSR the graph libraries see,
      not the `Graph` itself, which `DistMatrix` still needs with its diagonal to build the permuted
      `SpMatrix`. Probe against the fixed header: bijection at 2 and 4 parts with the loops present.
      Two serial regression tests added (`BuildAdjacencyDropsSelfLoops`,
      `BuildParAdjacencyDropsSelfLoops`). Code jury round 6.

Adjacent, pre-existing, noted (Grok, round 3), not filed: a rank that *has* vertices but no
edges gets a length-0 edge buffer from `graphtools.hpp:216-224` (the placeholder is only for
zero-vertex owners), so `data()` may be null there — degenerate for nested dissection, not
reachable from a connected FEM graph. The serial `scotch_nd` reads `tNumEdges = tEdges.length()`
against the serial one-element placeholder (`fn_Graph_SCOTCH.cpp:48`, `graphtools.hpp:99-107`) —
different API without `edgelocsiz`, only for a graph with no edges at all; do not copy D3's split
there.

Sibling pattern, out of scope, noted for whoever next touches it: `cl_FEM_Kernel.cpp:375-381`
guards an owner subscript with `tNeighborOwner < mCommSize` alone — same signed-`proc_t` gap as
D1 (skips `gNoOwner`, admits `-1`). Not reachable with today's owner values.

## 3. Ordered steps

- [x] **R1** — Fix D1. Standalone, independently reviewable, no behaviour change while dead.
      Implemented 2026-09-04 per §6; diff jury in `tmp/ai_exchange/review_parmetis_d1.md`.
- [x] **R2** *(after R1)* — decided (O1 ruling) and implemented in §7, 2026-09-04. Decide and record the call site: which solver path requests a parallel
      nested-dissection ordering, and on what condition (see O1).
- [x] **R3** *(after R2 and D4)* — done as W-R2 / W-O2, 2026-09-04. Wire it behind the existing `BELFEM_PARMETIS` / `BELFEM_PTSCOTCH`
      guards. ~~including the fallback when neither is compiled in~~ *(superseded by W1: a deck
      value whose library is not compiled in is a collective `BELFEM_ERROR`, never a silent
      fallback; the only fallback is W3's D4 path, and it announces itself)*.
- [x] **R4** *(after R3)* — done as W-R4, green 2026-09-04. Test. **Nothing in `make check` exercises this path today**, and
      `make check` runs MPI at np 2/4 only for `tests/comm` — so a green suite is not a gate for
      this code. A new test must actually run multi-rank (`TESTRANKS`, `Add_Test.cmake:141-172`),
      one per compiled backend, and must (a) include an owner with **zero** vertices so D3's path
      is exercised, and (b) check on rank 0 that the result is a bijective permutation that
      preserves adjacency — "enters the code" alone would not expose D3.
- [ ] **R5** *(after R4)* — Benchmark against the current ordering on a deck large enough to matter;
      an ordering that is not faster is not worth the dependency.

## 4. Open design questions

- [x] **O1 — RESOLVED 2026-09-04 (Christian's ruling: contiguous row blocks; see §7).** Who calls it, and is the graph in the right place?
      `build_pargraph_adjacency` asserts it runs on **rank 0** with the complete graph
      (`graphtools.hpp:133`) and then redistributes. Confirm that matches how the caller holds the
      graph, or the assert fires at the first real use. **And:** the only serial ND callers today
      build **ownerless** graphs — `sparse::create_graph_from_matrix` never calls `set_owner`
      (`fn_create_graph_from_matrix.cpp:39-41`), and `fn_compute_permutation.cpp:38` /
      `cl_SolverDistMatrix.cpp:87` hand those to `metis_ndp`. If R2 wires `parmetis_nd` there,
      D1's remap sends every vertex to rank 0 and the parallel ND degenerates to a serial ND
      with extra MPI. D1 is a crash fix, not a partition: the caller must supply a real
      distribution first. Design question for Christian.
- [x] **O2 — RESOLVED 2026-09-04 (W1: explicit deck values `parmetis` / `ptscotch`, no automatic promotion, so no default is needed).** ~~Needs a default and a deck override.~~
- [ ] **O3 — Documentation.** `src/math/graph/doc/graph_usage_guide.md:796-810,1234,1257` already
      documents `build_pargraph_adjacency()` with a worked call example, written as though it were
      in use. When R3 lands that becomes true; until then it overstates. Update it with the wiring.
      The example at `:1258` also names `ParMETIS_V3_PartKway` (partitioning) where the wrapper
      calls `ParMETIS_V3_NodeND` (ordering) — fix so the wiring pass does not copy it.

## 5. Definition of done

- [x] D1 fixed, with the always-on `BELFEM_ERROR` range check that prevents recurrence (R1)
- [x] D2 and D3 adjudicated by Christian and fixed (code jury round 3 in the exchange)
- [x] D4 resolved as part of the wiring (collective verdict in `parmetis_nd`, tested both ways)
- [x] Call site wired behind the existing build guards (`DistMatrix::order_graph`; the only fallback is D4's, announced)
- [x] A test that runs at more than one rank and actually enters the code (`sparsempi_np2/_np4`, green)
- [ ] Benchmarked against current ordering
- [x] `graph_usage_guide.md` reconciled with reality (O3; Codex language sweep not done)

---

## 6. R1 implementation plan (pre-registered 2026-09-04, Claude)

Scope: D1 only. No behaviour change for any current caller (there are none). Reviewed by the
plan-jury before coding; the diff goes to a second jury before commit.

### 6.1 The edit — `graphtools.hpp:136-146`

Normalise the owner once into a local, keep the existing "unowned → root" policy for the **two**
sentinels actually in use, then refuse anything else with an always-on check:

```cpp
Vector< T > tCount( tCommSize, 0 );
for ( Vertex * tVertex : aGraph )
{
    tVertex->unflag();

    // unowned vertices (e.g. circuit dofs) go to root. Two sentinels are in
    // use: the Kernel's comm_size() marker and the gNoOwner default.
    proc_t tOwner = tVertex->owner();
    if ( tOwner == tCommSize || tOwner == gNoOwner )
    {
        tOwner = 0;
        tVertex->set_owner( tOwner );
    }
    BELFEM_ERROR( tOwner >= 0 && tOwner < tCommSize,
        "Vertex %lu has owner %d outside [0, %d)",
        ( long unsigned int ) tVertex->id(), ( int ) tOwner, ( int ) tCommSize );

    ++tCount( tOwner );
}
```

Design decisions, each open to the jury:

1. **`BELFEM_ERROR`, not `BELFEM_ASSERT`.** The loop runs once per ordering (setup, by the
   call-rate rule in `CLAUDE.md` §Error Handling) and the release failure mode is a heap write.
   An assert would leave the release path exactly as quiet as it is now.
2. **Two named sentinels, not `>= tCommSize`.** *(Rationale corrected by the plan jury.)* The
   never-assigned case **is** `gNoOwner`, and both policies send it to rank 0 — the happy-path
   test below asserts exactly that. The two policies differ only on the open interval
   `( tCommSize, gNoOwner )`: a value there is garbage from a bug, not a sentinel anyone set,
   and naming the sentinels makes it die instead of quietly landing on root. Whether an
   all-unowned graph *should* reach this builder at all is O1, not R1.
3. **`tOwner >= 0` is load-bearing.** `proc_t` is a signed `int` (`typedefs.hpp:42`).
   `tOwner < tCommSize` alone passes `-1`, which the `Vector` subscript then turns into
   `2^64-1`. The pre-2026-08-31 `gNoID` default shows `-1` has already occurred in this code.
4. **Normalise once.** The later owner-indexed loops — the reorder at `:158-161`, `tNNZ` at
   `:197-200`, `tXcount/tAcount` at `:223-231`, and the permutation gather in the wrappers after
   return (`fn_Graph_ParMETIS.cpp:126-129`, `fn_Graph_PTSCOTCH.cpp:148-151`) — read
   `tVertex->owner()` after this loop, so writing the remapped value back with `set_owner`
   covers them. No second guard is added downstream. (`:149-153` is a prefix sum over `p`, not
   an owner subscript.)

Not changed: the `comm_size()` / `tCommSize` mixing in the distribution loop (`:148-152`) is
cosmetic and out of scope for a defect fix (minimal-edit rule).

### 6.2 The regression tests — `tests/math/test_GraphTools.cpp`

Four serial tests next to the existing `BuildAdjacency*` group, using the same
`build_test_graph` fixture, whose vertices carry the **default** owner `gNoOwner`
(`cl_Graph_Vertex.hpp:42`) because the fixture never calls `set_owner`. The first is the D1
happy path:

```cpp
TEST( GraphTools, BuildParAdjacencyUnownedGoesToRoot )
{
    // fixture vertices keep the gNoOwner default: this is the D1 path
    belfem::Graph tG = make_path5();
    belfem::Vector< int > tDistribution;
    belfem::Cell< belfem::Vector< int > > tVertices;
    belfem::Cell< belfem::Vector< int > > tEdges;

    belfem::graph::build_pargraph_adjacency( tG, tDistribution, tVertices, tEdges );

    const belfem::proc_t tSize = belfem::comm_size();
    ASSERT_EQ( tDistribution.length(), static_cast< size_t >( tSize + 1 ) );
    EXPECT_EQ( tDistribution( 0 ), 0 );
    EXPECT_EQ( tDistribution( tSize ), 5 );
    EXPECT_EQ( tVertices( 0 ).length(), 6u );  // all 5 on root, N+1 offsets
    EXPECT_EQ( tVertices( 0 )( 5 ), 8 );        // 4 undirected edges × 2
    for ( belfem::graph::Vertex * tV : tG ) { EXPECT_EQ( tV->owner(), 0 ); }

    belfem::graph::clear( tG );
}
```

`test_math_main.cpp` initialises `gComm`, so `comm_size()` is 1 under `make check` and the
`comm_rank() == 0` assert at `graphtools.hpp:133` holds. The test enters the `gNoOwner` branch
in both debug and release. It also locks a **degenerate** distribution (everything on root) as
the golden path; that is the correct D1 regression and must not be read as "an ownerless graph
is prepared for parallel ND" (see O1).

The other two exercise the range check. `test_math_main.cpp:26` calls
`set_throw_on_error( true )`, so `BELFEM_ERROR` throws `std::runtime_error` even under `NDEBUG`
and `EXPECT_THROW` is the project pattern (278 uses under `tests/`, no death tests). Both
`clear()` the graph after the throw, which happens in the first loop, before the builder has
moved any vertex:

```cpp
TEST( GraphTools, BuildParAdjacencyRejectsNegativeOwner )   // owner -1  → throws
TEST( GraphTools, BuildParAdjacencyRejectsOwnerAboveCommSize ) // owner comm_size()+1 → throws
```

`comm_size()+1` is 2 at one rank, which is neither sentinel, so it lands in the open interval
that decision 2 says must die.

The fourth, added after the code jury (Codex and Grok both flagged the gap), sets every owner to
`comm_size()` — the Kernel's own marker — and asserts the same all-on-root result as the happy
path. Without it, dropping the `== tCommSize` clause would have left the suite green.

```cpp
TEST( GraphTools, BuildParAdjacencyCommSizeMarkerGoesToRoot )
```

### 6.3 Gate

`make check-fast` (the `math` tests carry the `fast` label). Success = the four new tests pass and
the `BuildAdjacency*` group is unchanged. **Passed** — Christian ran the full `make check`,
2026-09-04. This is R1's executable gate; it does not stand in for
R4.

---

## 7. R2/R3 wiring plan (pre-registered 2026-09-04, Claude; ruling by Christian: distribute by contiguous row blocks)

### 7.1 Findings that shape the plan

- **The only live nested-dissection call site is `DistMatrix::create_matrix`**
  (`cl_SolverDistMatrix.cpp:55-110`). With the permutation switch on it is reached from
  `SolverPETSC` at `comm_size() > 1` (`cl_SolverPETSC.cpp:407-409`, `PETScAIJ` derives from
  `DistMatrix`); `DistMatrix` is also constructed by `SolverSTRUMPACK` (`cl_SolverSTRUMPACK.cpp:146,151`)
  and `DofMgr_EigenValues` (`cl_FEM_DofMgr_EigenValues.cpp:713-724`, which pins `NATURAL`), where the
  switch is off. `set_permutation_switch` is **not STRUMPACK and `reordering scheme : metis`**
  (`cl_SolverDistMatrix.cpp:41-52`) — the STRUMPACK exclusion is load-bearing (it runs its own ND
  and redistributes) and **stays** when the condition grows. `sparse::compute_permutation` has **no callers** in `src/` or `nonfree/` (only the
  usage guide mentions it). So the parallel-ND call site is inherently multi-rank, and it is PETSc's.
- **Today the ordering runs on rank 0 alone**: non-root ranks return from `create_matrix` right after
  `broadcast( mNumRows )` (`:57-63`), and `metis_ndp` runs on root only (`:87`). A collective
  `parmetis_nd` cannot be dropped in at `:87`; the non-root early return has to move below the call.
- **Root holds the complete matrix and graph** (`cl_SolverDistMatrix.cpp:25-26`, `:59-63`), which is
  exactly the input contract `parmetis_nd` / `ptscotch_nd` already implement (root builds every
  rank's slice, `graphtools.hpp:133`; non-root ranks pass an empty graph and `receive`). O1 is
  therefore answered by the existing shape: **the caller supplies owners, the wrapper distributes.**
- **The graph is ownerless** (`fn_create_graph_from_matrix.cpp:39-41`). Without an ownership step
  D1's remap sends everything to rank 0 and, with D4, ParMETIS then refuses the empty ranks.
- **`DistMatrix` already splits rows into contiguous blocks** in `determine_sizes_and_offsets`
  (`:170-190`: `tDiv` rows on the first `tSplit = P − (N mod P)` ranks, `tDiv + 1` on the rest) —
  but on the *permuted* matrix, after ordering. ParMETIS needs a distribution *before* ordering, and
  it only needs to be a working distribution: its output permutation is global and `DistMatrix`'s
  own post-permutation block split is unchanged. Using the same split formula on the **original**
  row index is the ruling's "contiguous row blocks".
- **`SolverParameters` is synchronised positionally** (`synchronize()`, `cl_SolverParameters.cpp:365-420`,
  `mReorderingMethod` as a `uint` slot), so new enumerators must be **appended** before `UNDEFINED`;
  the string parser loops `k < UNDEFINED` (`en_SolverEnums.cpp`), so appended values parse for free
  once `to_string` knows them.
- **Every switch on `ReorderingMethod`** that must learn the new values: `strumpacktools.cpp:83-116`
  (its `default: // pass` would silently keep STRUMPACK's default), **`cl_SolverMUMPS.cpp:115-140`**
  (its `default` writes `AUTOMATIC` — a silent misconfiguration; missed in the first draft, caught
  by both auditors), `cl_SolverPETSC.cpp:727-730`, `cl_SolverDistMatrix.cpp:48`,
  `en_SolverEnums.cpp` `to_string`. Note the **naming overload**: for MUMPS and STRUMPACK the
  existing `metis` / `scotch` values *already* engage ParMETIS / PT-Scotch at `comm_size() > 1`
  (`cl_SolverMUMPS.cpp:117-131`, `strumpacktools.cpp:88-110`); the new values must behave
  identically there or the deck value is a trap.
- **Where the parameters become uniform across ranks:** `SolverParameters( const input::Section * )`
  parses on whichever ranks the factory runs (`cl_MaxwellFactory.cpp:2698-2701`);
  `Solver::create_wrapper()` calls `synchronize()` first on every rank (`cl_Solver.cpp:43-46`). A
  check that must abort collectively belongs **after** that `synchronize()`.
- **Latent, not introduced here:** `create_graph_from_matrix` inserts the diagonal, so the graph
  carries self-loops (`fn_create_graph_from_matrix.cpp:50-57`) although `build_graph_adjacency`
  says "no self-loops expected". METIS tolerates it today; if ParMETIS's input check objects, that is
  what W6 tests 2 and 5 will see first — read a ParMETIS input error there as this, not as a
  `DistMatrix` bug.

### 7.2 Design decisions (open to the jury)

- **W1 — Explicit deck values, no automatic promotion.** `reordering scheme : parmetis | ptscotch`,
  appended as `PARMETIS = 4`, `PTSCOTCH = 5`, `UNDEFINED = 6`. `metis` keeps today's behaviour
  (serial `metis_ndp` on root) so no existing deck changes behaviour; `automatic` unchanged. A
  value whose library is not compiled in is a **`BELFEM_ERROR` in `Solver::create_wrapper()`, right
  after `synchronize()`** (setup-time, unsupported configuration, and collective: every rank holds
  the same enum and every rank errors; a parse-time check on rank 0 alone would be a debug-throw
  strand) — not a silent fallback. Answers O2: no "both available" default is needed because
  neither is chosen unless the deck asks. The input contract states in one sentence that for MUMPS
  and STRUMPACK `metis` / `scotch` already mean ParMETIS / PT-Scotch in parallel, and that
  `parmetis` / `ptscotch` select the same thing there and the BELFEM-side parallel ND for PETSc.
- **W2 — The distribution step is a graph helper called by the owner of the graph.**
  `graph::block_distribution( Graph & aGraph, const proc_t aCommSize )` in `graphtools.hpp`,
  declared **`inline`** (the header holds only templates today; a plain function there would be an
  ODR violation): `owner( k ) = p` for the k-th vertex in graph order, with the
  `determine_sizes_and_offsets` split written as its **two loops** (`tDiv` on the first
  `P − (N mod P)` ranks, `tDiv + 1` on the rest, remainder on the *last* ranks) — never
  `k / tDiv`, which divides by zero at `N < P`. Root calls it in `create_matrix` before the
  collective call. Explicit over implicit: the wrapper does not guess a distribution for an
  ownerless graph.
- **W3 — D4 is closed inside `parmetis_nd`, collectively.** After `build_pargraph_adjacency` root
  checks every `tDistribution( p+1 ) − tDistribution( p ) ≥ 1` into an **`int` verdict** (not
  `bool`; a 0/1 is clearer and avoids the `MPI_CXX_BOOL` fuss). The verdict is broadcast **in both
  rank branches**, right after the existing `comm_barrier()` and **before** `broadcast( tSeed )` and
  the `share` / `receive` pairs, which are not collective:
  ```
  root                                   non-root
  build_pargraph_adjacency               —
  comm_barrier()                         comm_barrier()
  broadcast( tOk )                       broadcast( tOk )
  if !tOk: message(); metis_ndp(); return   if !tOk: return
  broadcast( tSeed ); share; distribute  broadcast( tSeed ); receive ×3
  ```
  The fallback **announces itself** (`message( InfoLevel::Warning, … )`): a user who asked for
  `parmetis` gets a serial ordering and is told so — a configuration surprise, not an algorithmic
  retry. The graph root hands to `metis_ndp` has been reordered by owner and re-indexed by the
  builder; `metis_ndp` and `DistMatrix` depend only on `id()` and the final `index()` bijection.
  No rank ever reaches `ParMETIS_V3_NodeND` with an empty slice, and the reaction is identical on
  all ranks. PT-Scotch accepts empty ranks (D3), so `ptscotch_nd` gets no such guard.
- **W4 — `create_matrix` restructure.** `mReorderingMethod` becomes a `const` member set in the
  initializer list by `set_reordering_method( aParams )` — root reads, `broadcast`, same pattern as
  `set_permutation_switch`, whose condition grows to
  `type() != STRUMPACK && ( METIS || PARMETIS || PTSCOTCH )` — **the STRUMPACK exclusion stays.**
  The inviting one-line edit (`metis_ndp` → `parmetis_nd` at today's `:87`) is a **hang**: everything
  after `broadcast( mNumRows )` is root-only today, and `parmetis_nd` is collective. Non-negotiable
  control flow:
  ```
  all : comm_barrier ; broadcast( mNumRows )                        // unchanged
  root: aMatrix->set_indexing_base( Cpp )                            // unchanged, root-only
  if ( ! mPermutationSwitch ) { root: mMatrix = aMatrix ; return ; } // both ranks return, as today
  root: create_graph_from_matrix( *aMatrix, tGraph )                 // non-root: tGraph stays empty
  all : this->order_graph( tGraph )                                  // NEW private method — every rank enters it
  if ( ! root ) return
  root: permutations, permuted SpMatrix, index permutation, delete   // unchanged
  ```
  `order_graph`: `PARMETIS` → root `block_distribution`, **all ranks** `parmetis_nd`; `PTSCOTCH` →
  root `block_distribution`, **all ranks** `ptscotch_nd`; anything else → `metis_ndp` (a no-op on
  the non-root empty graph, `fn_Graph_METIS.cpp:104-107`, so today's root-only line keeps its
  meaning). `aMatrix` is still never touched off root; root still owns and deletes the vertices.
- **W5 — MUMPS and STRUMPACK map the new values onto what `metis` / `scotch` already do there.**
  *(Revised after the plan jury: an explicit error would make `parmetis` fail on the two solvers
  where `metis` already means ParMETIS.)* `strumpacktools.cpp`: `PARMETIS` → the `METIS` branch's
  strategy pair (PARMETIS at `comm_size() > 1`, METIS serially), `PTSCOTCH` → the `SCOTCH` pair —
  identical code paths, zero new behaviour. `cl_SolverMUMPS.cpp`: `PARMETIS` → the `METIS` pair
  (serial METIS / parallel PARMETIS), `PTSCOTCH` → the `SCOTCH` pair. PETSc's own-ordering switch
  (`cl_SolverPETSC.cpp:727-730`) adds both to the `NATURAL` group — the BELFEM-side permutation is
  already applied upstream, exactly as for `metis`. No `default:` swallows either value anywhere.
- **W6 — Tier-2 test binary `sparsempi`** in `tests/sparse` with `TESTRANKS 2 4`, following
  `tests/comm`'s two-binary layout. **Two** things are lifted into a shared header
  (`tests/common/tier2_launcher_sentinel.hpp`, included by both Tier-2 mains; `Add_Test.cmake` gets
  `include_directories( ${CMAKE_SOURCE_DIR}/tests/common )` since it adds no `tests/` path today):
  the launcher sentinel test, with its rank floor **parameterised** rather than the comm-specific
  text copied, and the **`allreduce` fold of `RUN_ALL_TESTS()`** (`test_commmpi_main.cpp:122-134`)
  — measured 2026-08-30: without it a rank-3-only failure is green at np 2. The suite's
  `CMakeLists.txt` registers `sparsempi` only under `if( USE_METIS OR USE_SCOTCH )`, because
  `Add_Test.cmake` gates on `USE_MPI` and `TESTRANKS` alone and an `#ifdef` inside the sources would
  leave a backend-less binary green; no `fast` label (it must not leak into `check-fast`). Tests
  (each `#ifdef` on its library; test 5 additionally on `BELFEM_PETSC`):
  1. `block_distribution` + `build_pargraph_adjacency` at np 2/4 on a 2-D grid graph: slice sizes
     match the split formula, terminals equal per-rank degree sums, every global index appears once.
  2. `parmetis_nd` on the same graph: on root the index set is a bijection of `0..n−1`, ids
     untouched, adjacency preserved; non-root graphs stay empty.
  3. `parmetis_nd` on the D4 path, **two ways**: (a) `block_distribution` with `N < P` (3 vertices
     at np 4 — the production trigger, `tDiv = 0`), and (b) all owners hand-set to rank 0. In both
     the verdict fires on every rank, the fallback runs, the result is still a bijection — no hang
     (the test itself is the gate).
  4. `ptscotch_nd` on the grid, and with an empty rank (D3's path): bijection, no hang.
  5. End to end: `SolverPETSC` at np 2 on a 2-D Poisson matrix with `reordering scheme : parmetis`
     (and `ptscotch`) against the `metis` solution, **with the Krylov method and preconditioner
     pinned** (`preonly` + LU, or CG with a stated PC) — under the inherited ASM/ILU default the
     ordering changes the preconditioner and "same to tolerance" is a flaky gate. This is the line
     that proves the call site, not just the wrapper.
- **W7 — R5 benchmark is Christian's run**: tapestack3d (or the cheapest deck that reaches PETSc at
  np ≥ 4) with `reordering scheme : metis` vs `parmetis` vs `ptscotch`; record ordering time and
  total solve wall time. An ordering that is not faster is not worth the dependency (R5 as written).

### 7.3 Ordered steps

- [x] **W-O2** — *(implemented 2026-09-04, code jury round 5)* enum values, `to_string`, the collective library check in `create_wrapper()`, the
      **five** switches (W5), both input-contract artifacts (`doc/input_file_reference.md:213`,
      `doc/input_schema.yaml:368-373`).
- [x] **W-R2** — *(implemented 2026-09-04, code jury round 5; syntax-checked Debug tree, not run)*
      `graph::block_distribution`; D4 guard in `parmetis_nd` (W3); `DistMatrix`
      restructure (W4). Serial gate: existing `make check` unchanged (no serial path changes
      behaviour) — **owed**.
- [x] **W-R4** — **passed 2026-09-04: `make check` 19/19, `sparsempi_np2` and `_np4` green, 12/12
      tests each, after D5.** Tier-2 `sparsempi` with the shared sentinel; tests 1–5 above. *(written
      2026-09-04, `tests/sparse/test_SparseMPI.cpp`, 12 tests behind their `#ifdef`s; compiles;
      first `make check` 2026-09-04: both lines **Failed** because the binary was never built —
      `make check` builds only the suites named in the `foreach( TESTDIR … )` list at
      `CMakeLists.txt:473`, and `sparsempi` was not in it; added. The box closes when the two
      lines pass)*
- [x] **W-O3** — *(applied 2026-09-04; Codex language sweep not done)* `graph_usage_guide.md`: **replace** the example at `:1244-1260` (it puts the
      ParMETIS call inside `if ( comm_rank() == 0 )` and names `PartKway`; the real shape is root
      builds ownership, then `parmetis_nd()` collectively outside the guard); strike the "very large
      meshes / serial METIS memory too large" use case at `:826-829` — root still holds the complete
      graph (`graphtools.hpp:133`), so wiring parallelises the ordering *compute* after a root-side
      scatter and does not reduce rank-0 memory; describe `block_distribution` and the announced D4
      fallback. `sparse_usage_guide.md:1144-1153` (enum dump) and `:532-534` (setter examples) learn
      the two values; `:694-705` marks `compute_permutation` as uncalled or the function is retired
      (Christian's call — not in this plan's scope). Both input-contract artifacts carry the
      naming-overload sentence from W1.
- [ ] **W-R5** — benchmark (W7), Christian. A wall-time gate, **not** a memory-scaling gate (see W-O3).

Jury for this plan and for its diff: the xhigh row (MPI collectives are a safety boundary).
Plan jury round 4 (2026-09-04, Codex terra/xhigh + Grok 4.6/xhigh): design spine accepted by both;
every constraint above marked "revised" or "non-negotiable" was written in from their findings.
