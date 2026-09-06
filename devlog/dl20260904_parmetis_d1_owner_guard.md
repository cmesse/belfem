# ParMETIS/PT-Scotch D1: owner sentinel guard in `build_pargraph_adjacency`

**Date:** 2026-09-04
**Purpose:** Close D1 / R1 of `todo/parmetis_ptscotch_wiring.md` — the out-of-bounds owner
subscript in the shared parallel-adjacency builder — through a plan jury, then a diff jury.
**Module:** `src/math/graph`, `tests/math`
**AIs involved:** Claude (Fable 5.1; survey, plan, implementation, verification), Codex
(gpt-5.6-terra / high, both rounds), Grok (grok-4.6 / high, both rounds)
**Exchange:** `tmp/ai_exchange/review_parmetis_d1.md` (both rounds, one thread)

## Context

`parmetis_nd` and `ptscotch_nd` are fully written and called from nowhere. Their shared helper
`build_pargraph_adjacency` (`graphtools.hpp`) sizes a per-rank counter to `comm_size()`, remaps
only the `owner() == comm_size()` marker to root, and then subscripts the counter by owner. The
vertex default `gNoOwner` is `INT_MAX` (`typedefs.hpp:42,59`), so a graph that reaches the
builder without a full ownership pass performs an unguarded increment at offset 2^31−1 in release
builds. Christian's 2026-08-31 ruling was to wire the entry points rather than delete them, which
makes this reachable; the todo made D1 a blocker to land first.

## What was done

**Plan (todo §6) → jury round 1.** Pre-registered the edit shape and a single serial regression
test, then dispatched both auditors on the todo file. Verification against the tree, all cites
re-read:

- Edit shape confirmed by all three voices: read the owner once, remap both sentinels
  (`comm_size()` marker and `gNoOwner`) to root with write-back, then
  `BELFEM_ERROR( tOwner >= 0 && tOwner < tCommSize, … )`. Always-on because the loop is setup-rate
  and the release failure is a heap write (`coding_philosophy.md`, frequency rule); `>= 0` because
  `proc_t` is signed and a `-1` becomes `SIZE_MAX` in the subscript.
- **My rationale for "two named sentinels over a blanket `>= tCommSize`" was wrong** (Grok,
  confirmed): the never-assigned case *is* `gNoOwner`, and both policies send it to root — the
  happy-path test asserts exactly that degenerate distribution. The two differ only on the open
  interval `( tCommSize, gNoOwner )`. Policy kept, rationale rewritten in the todo.
- **"The range check cannot be tested without death tests" was a false excuse** (Codex and Grok
  independently): `test_math_main.cpp:26` sets `throw_on_error`, `BELFEM_ERROR` throws
  `std::runtime_error`, and `EXPECT_THROW` has 278 uses under `tests/`. Two throw tests added to R1.
- Plan-text defects fixed: §2 and §5 still said "assert"; my owner-loop citation `:153-157` pointed
  at a prefix sum over `p`; the "always faulted" history was overstated (implementation-defined
  conversion, no guaranteed fault).
- **Two new defects filed by Codex, single-raiser, awaiting Christian's adjudication, both blocking
  R3 rather than R1:** **D2** — `fn_Graph_ParMETIS.cpp:91` allocates the `sizes` output with
  `2·floor(log2 P)` entries where the ParMETIS manual specifies `2·P` (2 vs 4 at two ranks; the
  header on this machine declares the call without documenting the length, so the contract is cited,
  not read). **D3** — for an owner with zero vertices the builder allocates a one-element placeholder
  edge buffer while the CSR terminal stays 0, and `fn_Graph_PTSCOTCH.cpp:88` takes `edgelocnbr` from
  the buffer length, so an empty rank advertises one edge.
- Refuted: a METIS part could equal `comm_size()` and collide with the marker — parts are `0..k−1`.
- O1 sharpened (Grok): the only serial ND callers today build **ownerless** graphs
  (`fn_create_graph_from_matrix.cpp:39-41` → `fn_compute_permutation.cpp:38`,
  `cl_SolverDistMatrix.cpp:87`). Wiring `parmetis_nd` there without a distribution step turns the
  parallel ordering into a serial one with extra MPI. D1 is a crash fix, not a partition.
- Sibling one-sided check noted, not touched: `cl_FEM_Kernel.cpp:375-381`.

**Implementation.** `graphtools.hpp` count loop as above. `tests/math/test_GraphTools.cpp`: four
tests before the `graph::sort` group — `BuildParAdjacencyUnownedGoesToRoot` (fixture owners at the
`gNoOwner` default → distribution `{0, 5}`, CSR 6/8, every owner 0),
`BuildParAdjacencyRejectsNegativeOwner` (`-1` → throws), `BuildParAdjacencyRejectsOwnerAboveCommSize`
(`comm_size()+1` → throws, guarded `!= gNoOwner`), `BuildParAdjacencyCommSizeMarkerGoesToRoot`
(the Kernel's marker → root, added after round 2). Behaviour for the existing `comm_size()` marker is
unchanged; there are no live callers.

**Executable evidence.** `-fsyntax-only` of the test TU (which instantiates the template) with the
real `compile_commands.json` flags of both trees — `cmake-build-debug` (Debug) and `build` (Release,
MKL) — rc 0. Then Christian ran the full **`make check`: passes**, the four new tests included —
R1 is at the regression rung.

**Diff → jury round 2.** Both auditors: no logic defect in the hunk; guard, signedness check,
write-back and tier hold. Codex raised a **debug multi-rank deadlock** (rank 0 throws in
`BELFEM_ERROR`, peers wait at the barrier) — **refuted as a defect of this diff**: the reaction is
the documented policy (`assert.cpp:312` initialises throw-on-error to the assertion state;
`coding_philosophy.md:633`: a debug run throws in parallel *deliberately*), and the identical shape
already sits two lines away at `graphtools.hpp:194`. Two findings confirmed by both voices and
applied after the round: the `owner() == comm_size()` marker was untested (dropping that clause
would have left the suite green) — fourth test `BuildParAdjacencyCommSizeMarkerGoesToRoot` added;
and the throw-test comment "before any vertex is moved" overclaimed — earlier vertices are already
unflagged and remapped when a later one throws, only the reordered graph is unbuilt — comment
rewritten. Grok's wording point kept: the throw tests prove the check fires under the test hook, not
that production aborts. No third round — one comment and one same-shape test changed; the next gate
is executable. Re-syntax-checked in both trees after the edits, rc 0.

## D2 and D3 — same day, after adjudication

Christian accepted D1, asked for D2 to be checked against the ParMETIS source, and accepted the D3
fix. `KarypisLab/ParMETIS` `libparmetis/ometis.c` settles D2: the `sizes` argument is documented in
the source as "the 2*nparts array", `MultilevelOrder` sets `sizes[0] = 2*npes-1` and
`LabelSeparators` writes `sizes[--sizes[0]]`, so the largest index written is `2·npes−1`. The old
`2·floor(log2 npes)` allocation was short at every rank count from two upward.

- `fn_Graph_ParMETIS.cpp`: `tMySizes( 2 * tCommSize, 0 )`, with the source citation in a comment.
- `fn_Graph_PTSCOTCH.cpp`: `edgelocnbr` now comes from the CSR terminal `tVertices( tLocalNumVerts )`,
  `edgelocsiz` from the buffer length. On a rank with edges the two are equal, so nothing changes
  there; on an empty rank the phantom edge disappears while the one-element placeholder from the
  builder stays as legal slack storage.

Both trees define `BELFEM_PARMETIS` and `BELFEM_PTSCOTCH`, so the syntax check is not vacuous:
`-fsyntax-only` with the real flags, rc 0 for both files in both trees. Still dead code, still no
run. **Code jury round 3:** both auditors confirm the two calculations. Codex raised a
**blocking** point that is real and pre-existing: ParMETIS's own `CheckInputsNodeND`
(`libparmetis/weird.c`, read from the repository) rejects any processor with zero vertices —
"Poor initial vertex distribution" — and the builder hands an empty owner exactly such a slice. The
failure would be loud at the status check but rank-local. Not patched here; filed as **D4**, a wiring
decision (collective distribution validation or an upstream guarantee, tied to O1), blocking R3.
Grok caught that my D2 comment overstated the write indices (`sizes[0]` is a cursor; the API length
is `2·npes` regardless) — rewritten to state the documented length only — and a dead `<cmath>`
include, removed. Two pre-existing adjacent cases noted in the todo, not filed: a vertex-bearing
rank with zero edges gets a length-0 buffer, and the serial `scotch_nd` reads its edge count from
the buffer length too. Re-syntax-checked, rc 0 both trees. No fourth round.

## Wiring (R2/R3/R4), same day, after Christian's O1 ruling

Ruling: distribute by contiguous row blocks of the matrix layout. Survey first: the only live
nested-dissection call site is `DistMatrix::create_matrix`, reached by PETSc at more than one rank
(STRUMPACK is excluded from the pre-permutation, `EigenValues` pins `NATURAL`);
`sparse::compute_permutation` has no callers. Today the ordering runs on rank 0 alone because the
non-root ranks leave `create_matrix` right after the row-count broadcast — the inviting one-line
swap of `metis_ndp` for `parmetis_nd` is a hang.

**Plan jury (round 4, both xhigh).** Spine accepted: explicit deck values `parmetis` / `ptscotch`;
a caller-side `block_distribution`; D4 closed inside `parmetis_nd` collectively; all ranks through
a new `order_graph`; a Tier-2 `sparsempi` binary. Constraints written in from the findings: the
MUMPS switch I had missed (its `default` would have turned the new values into "automatic" while
`metis` there already means ParMETIS); R3's "fallback when not compiled" struck in favour of a
collective `BELFEM_ERROR` after `synchronize()`; the STRUMPACK exclusion kept explicit; the verdict
broadcast in both rank branches as an `int` before any `share`/`receive`; the two-loop split so
`N < P` cannot divide by zero; the Tier-2 `allreduce` fold lifted with the sentinel; a CMake-level
library gate; pinned CG + Jacobi for the end-to-end test; the usage guide's memory-scaling claim
struck (root still holds the whole graph). W5 revised from "error on STRUMPACK" to "map onto what
`metis`/`scotch` already do on MUMPS and STRUMPACK" — Grok's naming-overload point. One refutation:
my "implementation-defined" sentence was about the `uint`→`int` conversion, not `-1`→`size_t`.

**Implementation.** `ReorderingMethod` gains `PARMETIS = 4`, `PTSCOTCH = 5` (appended; the value
is synchronised positionally as a `uint`). Five switches learn them. `Solver::create_wrapper()`
refuses a value whose library is not linked, on every rank. `graph::block_distribution` in
`graphtools.hpp`. `parmetis_nd` computes the D4 verdict on root and broadcasts it in both branches;
on failure root warns and runs serial METIS, everyone returns. `DistMatrix` carries a broadcast
`mReorderingMethod`; `create_matrix` keeps every rank through `order_graph`. Tests:
`tests/common/tier2_launcher_sentinel.hpp` (sentinel with a parameterised floor, verdict fold — the
comm main now includes it), `sparsempi` under `TESTRANKS 2 4` with eleven tests: builder slices and
the `N < P` placeholder, `parmetis_nd` on an 8×8 grid plus both D4 triggers, `ptscotch_nd` with and
without an empty rank, and a PETSc Poisson solve under `metis` (control), `parmetis`, `ptscotch`.
Both input-contract artifacts and both usage guides updated; the guide example that put the
collective call inside the rank guard is replaced.

**Executable evidence.** `-fsyntax-only` with both trees' real flags after the final edits: rc 0
for eight source TUs and four test TUs. `check_doc_claims`
38/38. **Nothing has run.** `make check` now registers `sparsempi_np2` and `_np4`; those two lines
are W-R4's gate, and `make check` as a whole is W-R2's serial gate.

**Code jury (round 5, both xhigh).** MPI pairing, the `create_matrix` rewrite, the enum append and the
Tier-2 lift accepted by both. Revised after the round: the library check was scoped to PETSc
(Codex: a STRUMPACK tree without `USE_METIS` would have been refused `parmetis` although STRUMPACK
maps it itself — `config_strumpack.cmake` requires SCOTCH, not METIS); the PETSc test now mirrors
the production contract, matrix and vectors on rank 0 only, empty objects elsewhere; the METIS
control test gained its `#ifdef`; a one-row PETSc solve exercises the D4 fallback through the whole
`DistMatrix` tail; two comments and the guide's line cites corrected. Two library-contract questions
settled from the ParMETIS repository: a non-power-of-two rank count is truncated internally
(`npes = 1<<log2Int(npes)`), not refused — which also explains the original `2·floor(log2 P)`
sizing — and self-loops are neither rejected nor stripped at input, so the end-to-end PETSc test
is the only thing that will tell whether ParMETIS's coarsening tolerates the diagonal the matrix
graph carries. The library's input check folds its verdict over the communicator, so the pre-D4
failure was collective in release; the rank-local strand I had claimed holds only for the debug
throw. Re-syntax-checked in both trees after the edits. No sixth round; the gate is `make check`.

**First `make check` (Christian): `sparsempi_np2` and `_np4` Failed.** Not the code: no
`test_sparsempi` binary existed in either tree and its object directory was empty. `make check`
does not build every registered test — it depends only on the suites named in the hand-maintained
`foreach( TESTDIR … )` list at `CMakeLists.txt:473`, where `commmpi` had been added by hand when it
appeared. `sparsempi` was registered with ctest and absent from that list, so ctest reported a
missing executable. Added to the list with a comment naming the failure mode.

**Second `make check`: the tests ran. Ten of twelve pass at np 2 and np 4** — the launcher
sentinel, both builder tests, ParMETIS on the block-distributed grid, both D4 fallbacks (the
warning prints, serial METIS takes over, every rank returns), PT-Scotch with and without an empty
rank. The PETSc section crashed — in the **`metis` control test**, before either new value was
tried: `munmap_chunk(): invalid pointer` inside `METIS_NodeNDP`
(`FM_2WayNodeRefine1Sided → rpqDestroy → free`), called from `DistMatrix::order_graph`'s default
branch. That is today's production path for PETSc at more than one rank with
`reordering scheme : metis`, latent because no example deck uses the combination. Grok had flagged
the self-loops the matrix graph carries (one per diagonal entry) in rounds 4 and 5 as a ParMETIS
risk; plain METIS broke on them first. A serial scratchpad probe (64-row Poisson graph through
`create_graph_from_matrix` and `metis_ndp`) reproduced the crash at 2 and 4 parts and ran clean
with the loops stripped. **D5** fixed in the two adjacency builders — they drop self-loops when
counting and filling the CSR the graph libraries see, while the `Graph` itself keeps its diagonal,
which `DistMatrix` needs to build the permuted matrix. The probe rebuilt against the fixed header
orders the looped graph cleanly at 2 and 4 parts; two serial regression tests added to the graph
suite; syntax-clean in both trees.

**Third `make check`: 19/19, `sparsempi_np2` and `_np4` green, 12/12 each.** The three PETSc
end-to-end solves — `metis` control, `parmetis`, `ptscotch` — and the one-row `parmetis` fallback
through the whole `DistMatrix` tail all pass at both rank counts. R2, R3, R4 and D4, D5 are closed
at the regression rung.

**Code jury round 6 (D5 diff).** Both auditors accept the layer — drop loops in the CSR the
libraries see, keep them on the `Graph` — and both find the same two holes, fixed the same hour:
the parallel builder chose its one-element placeholder by *vertex* count, so a slice whose vertices
had only self-loops would now get a zero-length buffer and a possibly null pointer (branch moved to
the edge count); and the serial SCOTCH wrappers still took their edge count from the buffer length
(now the CSR terminal, as PT-Scotch already did). Grok's test points taken: the serial test pins
that the `Graph` keeps its loops, the parallel test scans its CSR for self-indices, an all-loop
fixture checks the placeholder in both builders, and `metis_ndp` itself runs on a looped graph under
`BELFEM_METIS` — the crash's own shape in the serial suite. The usage guide's pseudocode, which
still described the unfiltered count, and the graph test catalogue are updated. `make check` owed
once more for these test additions; nothing on the ParMETIS/PT-Scotch path changed after the green
run.

**Disposition.** The plan moved to `todo/closed/` with its one remnant, the W-R5 benchmark, filed as
**DR-156** in the debt register (`[RUN] [P]`); the register's `[P]` count went 7 → 8 and
`check_doc_claims` holds at 38/38. Open beside it: `compute_permutation` is dead code awaiting
Christian's call, and the new guide prose has not had its Codex language sweep.

## Tripwires

L-11 (serial correctness says nothing about ownership) applies by subject. R1 changes no ownership
semantics — it guards an index — and its tests are serial by design; the multi-rank gate is R4,
whose wording now requires an empty-owner partition and a bijection check (Codex F4).

## State of the todo

R1 ticked. R2–R5 open: R2 (call site) is blocked on O1, which is now a design question for Christian;
D2, D3 fixed; D4 closed by the verdict in `parmetis_nd`; O1, O2 resolved by ruling; R2/R3 implemented; R4 written. Open: `make check` (serial + Tier 2), W-R5 benchmark, `compute_permutation` disposition, Codex language sweep of the guide text.

Status: **verified** — D1 through D5 and the wiring (R1–R4) are behind `make check` 19/19 with
`sparsempi` at np 2 and np 4. Open: DR-156 (W-R5 benchmark, Christian's run, a wall-time question);
`compute_permutation` disposition (dead, Christian's call); Codex language sweep of the new guide
text; one more `make check` for the round-6 test additions. Six jury rounds in the exchange file.
Plan closed.
