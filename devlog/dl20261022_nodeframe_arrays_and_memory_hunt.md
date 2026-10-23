# Yufan's 0.9.3 memory report: the Vertex counters are not it, the postprocessor's per-node `NodeFrame` was a `Vector` misuse and is now `real[3]`

**Date:** 2026-10-22
**Purpose:** Yufan reported that a 3.5 M-dof 3D run needs much more memory under 0.9.3 than 0.9.2 and attributed it to the `mesh::Vertex` counters widening from `uint8_t` to `uint16_t`, proposing to narrow them again behind an overflow flag. Christian asked for an independent opinion, then a jury on whether anything else in v0.9.2..v0.9.3 explains the increase; the hunt found the recovery postprocessor's `NodeFrame`, which Christian diagnosed as a coding-guideline violation (`Vector< real >` for 3-component geometry) and asked to fix through the standing loop.
**Module:** `src/mesh` (analysis only), `src/fem/kernel` (`cl_FEM_Postprocessor`), `tests/fem`
**AIs involved:** Claude Fable 5.1 (primary: analysis, probes, plan, implementation, sandbox gates); memory jury: Codex `gpt-6-astra`/high **failed on vendor quota**, Grok `grok-4.7`/high (one-auditor round); plan jury and code jury on the `NodeFrame` change: Codex `gpt-6-astra`/high and Grok `grok-4.7`/high, both seats both rounds
**Status:** COMPLETE — the `NodeFrame` change landed and is verified; the memory question is answered as far as source and probes reach, with the deck-level measurement handed to Yufan
**Verification:** executed, in Claude's private sandbox trees (detached worktree of `90d2e07` + two build trees in the scratchpad, GNU 16.1.0 from `/opt/scls`, `USE_DEBUG=ON`, `USE_MAXWELL=ON`; the shared `build/` and `cmake-build-debug/` untouched), serial runs of `./test/test_fem`. Positive control, pre-change sources + new tests, Blaze: `PostprocessorRecovery.*` 15/15. Change, Blaze: 15/15, zero compiler warnings. Change, Armadillo (full build): 15/15, zero warnings; whole `test_fem` 337/337 in 59 s. Mutation (`T1 := Normal` in the 2D branch, Blaze): `ThicknessJumpSurvivesCurvature2D`, `ConstantReproducedOnCurvedShell2D`, `TangentLinearFieldReproducedOnFlatShell2D` red; restore 15/15. Evidence level: regression + mutation. Not run: the other `check-fast` suites (nothing outside `fem` changed), Release builds, MPI ranks (the recovery cases skip there by design), any deck-level memory measurement. The memory-hunt findings themselves are source trace plus two standalone `sizeof` probes — reviewed, not verified.
**Literature References:** N/A (Zienkiewicz & Zhu 1992 patch recovery untouched)
**Exchange:** `tmp/ai_exchange/review_memory_092_093.md` (hunt), `tmp/ai_exchange/review_nodeframe_arrays.md` (plan round, code round, gate)

## Summary

The counter widening costs 8 bytes per `Node`, `Edge`, `Face`, `Facet` and `Segment` — tens of MB
on Yufan's mesh, not GB — and there is no sound 0.9.2 baseline for a mesh that needs counters above
255, because a 0.9.2 release build wraps the `uint8_t` and writes past the adjacency array. The
largest code-side memory delta between the tags is the recovery postprocessor's per-node
`NodeFrame` (commit `5ab4d41`): four `Vector< real >` for 2/3-component geometry, 784 B per node
under Armadillo, heap-allocated and held for the run, one map per postprocessor instance. That is
now four `real[3]` held by value — 112 B, no heap — with four new 2D recovery tests, verified on
both backends. What the code cannot explain, Yufan's deck can: the factor fill-in after the
interface-copy renumbering and RCM tie-break change, a stale `.bfm`, or the mesh's own valence.

## The memory question (hunt round, one auditor)

**Counter widening (high).** `b26b660` + `ef70431` turned ten `uint8_t` counters/capacities in
`cl_Vertex.hpp:49-66` into `uint16_t`. Probe against the v0.9.2 and current headers (debug-tree
flags): `mesh::Vertex` 136→144, `Node` 184→192, `Edge` 144→152, `Face`/`Facet` 176→184; `Element`
128→128 (derives from `Basis`, `cl_Element.hpp:34`); `Segment` 152→160 by the same layout (Grok).
The pointer arrays stay sized by count. With glibc's 16-byte chunk rounding the +8 lands on `Node`
and `Segment` (+16 each) and vanishes on `Edge`, `Face`, `Facet` — Grok had the Node and Edge rows
swapped; conclusion unchanged. ~10–100 MB on a 3.5 M-dof mesh, concentrated on rank 0.

**No valid baseline (high; both voices independently).** Callers increment then call
`allocate_*_container()` with no argument (`cl_Mesh_ConnectivityCalculator.cpp:74`), so a 0.9.2
release build wraps the counter mod 256 (the increment guard is `BELFEM_ASSERT`), allocates the
wrapped count, passes the capacity `BELFEM_ERROR`, and `add_*` (assert-only bound,
`cl_Vertex.cpp:90-91`) overruns the heap. Yufan's 0.9.2 number is from an unsound or a different
run. Grok adds that the two binaries disagree on the linear system even below 255 (interface
copies, RCM order, gauge pin, thin-shell orientation). The increment path is still assert-only, so
a release build wraps at 65536 the same way — a design question for Christian, not a delta.

**`NodeFrame` (Claude, single raiser in that round — Grok read the changelog, not the diff).**
`cl_FEM_Postprocessor.hpp:80-90` as of 0.9.3: `Vector< real > Origin, Normal, T1, T2`, allocated
`new` per selected node (`:693-694`), read at every recovery (`:1122-1125`), freed from the
destructor. `Vector` wraps `arma::Mat` 1:1 (`cl_AR_Vector.hpp:38`), 192 B with its
`mem_local[ 16 ]`; probe: `sizeof( NodeFrame )` = 784 B Armadillo, 112 B Blaze (which instead
heap-allocates each sized member: 3/5/6 chunks per frame). One instance per postprocessor type —
Air, Ferro, Conductor, SuperConductor, ThinShellConductor, ThinShellSuperConductor, SideConnector
(`cl_MaxwellFactory.cpp:2751-2767`, `:2836-2845`, `:2899-2915`). Estimate ~0.5–1 GB on Yufan's
run (medium). Not the dominant term; the direct factorization is (the 0.9.3 release gate measured
the pancake deck at 38 GB resident, `dl20261010_release_gate_0_9_3_decks.md:31-32`).

**Refuted.** Grok's low-confidence candidate that the MUMPS `ICNTL(23)` memory budget is new in
0.9.3: `git show v0.9.2:src/sparse/cl_SolverMUMPS.cpp` already has `memory_budget_mb`; the tag diff
of that file is the warning-bit decoder only. Grok's "`ef70431` unattested": it is an ancestor of
v0.9.3, the rewritten SHA of `168b861`. Any change to the distributor or partitioner: the tag diff
of both files is empty. No solver default changed; `cl_SolverParameters.cpp` moved validation into
setters.

**Confirmed but unverifiable from source (Grok).** Interface copies are numbered by sorted
original id and RCM tie-breaks via `std::sort` (`CHANGELOG.md:74-89`, commits `cfaf7e8`,
`e30bd98`): the factor ordering changed and fill may have moved by GB in either direction; shared
conductor-side copies *reduce* the dof count (`:61-73`); a `.bfm` written by 0.9.2 is still reused
and replays the old enrichment (`:90-93`). These are the questions for Yufan's deck: delete the
`.bfm`, log MUMPS `INFOG(16)`/`INFOG(36)` and peak RSS per rank, compare the max node and edge
valence of the mesh.

**On the proposal to narrow the counters again** (dynamic widening on overflow): the saving is the
8 B/entity above, ~tens of MB; a widen-on-overflow scheme adds a branch and an indirection to the
connectivity path for that. Not pursued. The `NodeFrame` fix below saves roughly ten times as much
per node.

## The `NodeFrame` change (plan round, code round, gate)

**Plan** (`todo/closed/postprocessor_nodeframe_arrays.md`): `struct NodeFrame { real Origin[ 3 ],
Normal[ 3 ], T1[ 3 ], T2[ 3 ]; real H, Tau; }` zero-initialised, `static_assert(
is_trivially_copyable )`, `Map< mesh::Node *, NodeFrame >` by value; every `Vector` operation on the
frame becomes a loop to `mNumDimensions`; `mOffset`, `mShellDirection` stay `Vector` (member scratch,
`compute_poly` consumes a `Vector`); `mNodeMatrices` out of scope (O4). Precedent:
`Node::mCoords` is `real[3]` (`cl_Node.hpp:32`).

**Plan jury, both seats, no refutation of the direction.** Codex and Grok both caught that the
plan as first written would not compile (the `*tFrame` dereference at the `add_node_gram` call
`:698`, `Origin( )` at `:959` and `:1075` not listed), and that the gate had no 2D case — all three
test fixtures are `Mesh( 3, … )`, the 2D tangent branch `:999-1003` had never run, and the closed
thin-shell test plan had left its 2D case (T5) unbuilt. Codex: `norm` → BLAS `dnrm2` is an
algorithm change (scaled summation), not summation order, and a round-off perturbation can flip the
SVD rank cut in `recovery_block_inverse` (`fn_FEM_recovery_inverse.hpp:108,119`) — "last bit only"
was too strong; the hot `dot` was already Armadillo's inline loop for n ≤ 32
(`op_dot_meat.hpp:156`), so "removes a BLAS call" and "cannot be slower" were false — the
justification is footprint and allocation count, runtime unmeasured; a 4×4 inverse fits
Armadillo's `mem_local[ 16 ]`, so O4's "192 + 128" double-counted; and **D1**: `initialize()`
nulls both selection bitsets (`:105-106`, `:128-129`) before `compute_node_matrices()`, so a second
`initialize()` after a throw dereferences the null `mElementBitset` at `:228` — the `:685` comment
describes a path that crashes earlier. Pre-existing, logged, not fixed. Grok: Blaze allocates
3/5/6 chunks per frame, not 2; `Vector` wraps `Mat`, not `Col`; `NodeFrame` is `private`; six of
my line citations were off; "a loop to 3 would be safe" is false for the `Vector` operands (sized
to `mNumDimensions`); T1's zero third component is load-bearing (the 3D branch writes `tA`, `tB`
only); add `<type_traits>`. All folded in before implementation.

**Implementation.** R1–R5 as planned, syntax-gated under the tree's `-Wall -Werror`. R8b: fixtures
`Shell2D( aCurved )` — two QUAD4TS layers (blocks 1, 2) on an arc of radius 1 or a flat tape of
length 0.2, node rows at `R − t, R, R + t`, LINE2 tape facets on the middle row, a real
`mesh::ThinShell`, node order (bottom0, bottom1, top1, top0) as `ThinShellFactory` (`:1464-1467`)
— and `FarTriPatch` (8 TRI3, Air, 1e3 from the origin); cases `ThicknessJumpSurvivesCurvature2D`,
`ConstantReproducedOnCurvedShell2D`, `TangentLinearFieldReproducedOnFlatShell2D`,
`LinearFieldFarFromOrigin2D`. **Christian's ruling:** the fixture helper `create_field` takes no
dimensionality parameter (I had added one; nothing in `src/` was touched by it); it now derives
`ModelDimensionality` from `mMesh->number_of_dimensions()`.

**Code jury, both seats: no blocking defect.** The one refutation of the session is Grok's, of my
own fixture comment: I wrote that the counterclockwise station order "keeps the mesh checker from
reordering the bottom edge"; `MeshChecker` swaps only on a negative `Pipette` measure
(`cl_MeshChecker.cpp:207-210`) and `measure_quad4ts` returns `std::abs(...)`
(`cl_Pipette.cpp:322-335`), so the swap never runs for QUAD4TS — the node order alone protects the
bottom edge, and the winding is immaterial (`dV_quad4ts` returns `|det J|`,
`cl_FEM_Calculator.hpp:2155-2186`; the patch weights are `_Volumes` = facet length × thickness).
Grok also cut two overclaims from the test comments (the flat-tape case accepts any `T1` with an x
component since `η` already spans `y`, and `s = −Δx` for the normal `+e_y`; the curved cases have no
in-plane dependence) and two from the header (lifetime; "every loop to `mNumDimensions`" — the
zeroing loop and the axis search correctly run to 3). Codex and Grok, from different angles, logged
**O5**: both shell fixtures are axis-aligned, so the least-aligned-axis branches of the 3D tangent
construction and a *tilted* 2D `T1` are not pinned. Comments corrected after Christian's approval.

**Gate** (table in `**Verification:**`). Two things the gate taught: the 2D thin-shell-only h-φ
mesh initialises and recovers through the dof manager on the first try (positive control green),
and the mutation shows the curved 2D cases *do* discriminate a `T1` parallel to the normal — the
`s` column becomes collinear with `η`, the Gram matrix goes near-singular and
`EXPECT_EQ( tTruncated, 0 )` fails (6 patches each) — which both Grok and I had called "does not
read the tangent". The flat-tape case remains the only one that would catch a tilted `T1`; O5
stands. The `ThicknessJumpSurvivesCurvature2D` comment now says so.

## Open

- **Yufan's deck:** `.bfm` deleted, `INFOG(16)`/`INFOG(36)`, peak RSS per rank, max valence — the
  measurement that settles the GB-scale question. Not a code item.
- **O4** `mNodeMatrices` (~192 B header + two allocations per node under Armadillo), **O5**
  recovery-frame coverage, **R8c** deck-level VmHWM — in the closed plan's §5.
- **D1** re-`initialize()` after a throw crashes on a null bitset — pre-existing, Christian's call
  whether it gets a debt-register row.
- The release-build `increment_*_counter()` is still assert-only and wraps at 65536 as 0.9.2 did at
  255 (Grok) — a tier choice, Christian's.

## Files Updated

`src/fem/kernel/cl_FEM_Postprocessor.hpp`, `src/fem/kernel/cl_FEM_Postprocessor.cpp`,
`tests/fem/test_PostprocessorRecovery.cpp`, `todo/closed/postprocessor_nodeframe_arrays.md` (new,
retired on completion), `todo/README.md`, `devlog/README.md`, this file.
