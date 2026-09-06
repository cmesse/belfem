# SPFA Rectify-or-Certify: clean_spfa() Implementation

**Date:** 2026-07-15 (evening session, same day as the diagnosis devlog
`dl20260715_corc_clean_divergence.md`)
**Purpose:** Implement the SPFA difference-constraint replacement for the greedy
`Cohomology::clean()` (plan T1/T4/T5/T6), Christian-approved.
**Modules:** src/math/graph, src/homology
**Exchange thread:** `tmp/ai_exchange/corc_clean_divergence.md` (audits appended there)

## What landed

**Christian's scaffold (pre-session):** `Cohomology::clean()` now dispatches through
the `mFunClean` member-function pointer; old body preserved as `clean_greedy()`;
`clean_spfa()` is the default in all three constructors — switching algorithms is a
one-line pointer change.

**1. `src/math/graph/fn_Graph_spfa.{hpp,cpp}`** (new; registered in the graph
CMakeLists): `graph::spfa_difference_constraints( nV, tails, heads, weights, theta,
cycle )` — a generic feasibility solver for `θ(head) − θ(tail) ≤ w` systems.
Queue-based Bellman–Ford–Moore with a virtual super source (all θ = 0), CSR arc
buckets, preallocated ring buffer sized to the active-vertex count (filled directly by
`DynamicBitset::where()` — Christian's side-agent refinement), `int64_t` distances,
negative-self-arc fast path (length-1 certificate for folded periodic edges),
update-count trigger with **verified** predecessor-graph cycle extraction, and a hard
op-cap `BELFEM_ERROR` failsafe. Per Christian's Q3 the solver lives in `math/graph`;
`Vertex::mLevel` was deliberately not used (potentials are signed 64-bit; `mLevel` is
`index_t`; external index-keyed arrays per the audited plan).

**2. `Cohomology::clean_spfa()`** (`cl_Cohomology.cpp`): per generator, two arcs per
`original_edges()` edge `e = (i, j)` with `u = rep(i)`, `v = rep(j)` (flagged slaves
fold onto masters): `(u→v, 1+c)` and `(v→u, 1−c)`, i.e. `|c + θ(u) − θ(v)| ≤ 1`.
- **Feasible:** fire every representative node `u` with weight `−θ(u)` via the static
  helper `cohomology::fire_node_coboundary()`, which replicates `clean_greedy()`'s
  firing pattern (master∪slave quotient adjacency, cochain sign convention, manual
  boundary bookkeeping) scaled by the weight — net `Δc = θ(rep(i)) − θ(rep(j))` per
  edge, so every coefficient lands in {−1, 0, 1}; zero entries auto-erase.
- **Infeasible:** release-active `BELFEM_ERROR` reporting the negative-cycle
  certificate as mesh edge IDs (≤ 32 printed) with a targeted-refinement
  recommendation — replacing the greedy stall/plateau failure mode.
- **Guards:** null-complex `BELFEM_ERROR`; int-overflow guard on firing weights;
  release-active post-write-back `|c| ≤ 1` sweep over the whole map (closedness is
  preserved by construction — only node coboundaries are added, `dd = 0`).
- One deliberate deviation from greedy: the boundary node fold uses `is_flagged()` in
  **both** adjacency passes (greedy's second pass used `is_periodic()`, which would
  mis-fold a master endpoint of another periodic pair onto its slave).

## Verification

- Both files compile clean under the real `flags.make` flags
  (`-Wall -Werror -pedantic-errors`, Blaze backend, mpicxx).
- Scratch test `test_spfa.cpp` (session scratchpad): 18 checks — feasible chains
  rectify to unit, the greedy-hang 3-cycle (all c = 2) returns a chained
  negative-weight certificate, self-loop and multi-component cases, 400 randomized
  instances with independently verified constraints/certificates, n = 2000 stress —
  ALL PASS; valgrind clean.
- **The test caught a real soundness bug before it reached the repo:** the classic
  SPFA update-count trigger (`cnt > n`) can fire **spuriously**, before the
  predecessor graph contains a cycle. First version treated failed extraction as an
  internal error → aborted on a feasible instance. Fix: only a predecessor-graph
  cycle (provably negative-weight) is accepted as an infeasibility proof; spurious
  triggers reset the counter and keep relaxing (Cherkassky–Goldberg style
  verify-or-continue), with the op-cap as the honest last-resort failure.

## Audit round (same session) — verdicts and applied fixes

Codex + Grok audited the implementation (thread: exchange file). **A1 sign chain:
PASS (both, high). A4 solver: no false-feasible path (both).** All residual findings
were applied the same evening and re-verified (syntax + 18-check suite + valgrind):

- **Zero-weight parent cycle is not an infeasibility proof** (Grok): the extractor
  now accepts only strictly negative cycle sums; zero cycles are dead-end walks.
- **Op-cap raised** to (V+1)·(E+1)+1000, past SPFA's O(V·E) worst case (Codex).
- **Self-paired periodic nodes** (`periodic() == this` in `slave_nodes()`, Codex):
  write-back now fires them (skip requires `periodic() != this`), and the helper's
  second pass guards against double-firing their own adjacency.
- **Zero-generator early return** no longer skips `unflag_everything()` (Codex).
- Accepted as documented limits: generic solver doesn't overflow-guard arbitrary
  int64 weights (production weights are 1±c, int c); certificate prints edge IDs
  without per-arc orientation (location certificate).

## T7 (same session, adapted)

`CutData::weight()` reads sign bitsets only, so a non-unit magnitude can never reach
`determine_cut_case_3d()` — the T5 guard in `clean_spfa()` is the real split. Landed
instead: the inadmissible-pattern `default:` no longer says "mesh might need a finer
resolution" (post-certificate era) and names a generator/closedness defect with the
pattern code; the two sign-coherence errors now print the per-edge coefficient values.
Also removed a stray block-scope `collect_thin_cut_edges` declaration pasted
mid-function in `determine_cut_case_3d()` (committed in `6ee74493`; legal C++, zero
effect, clearly accidental).

## Open

- **T3 verdict pending Christian's run:** build, then run the double-layer corc.
  Outcome A — all seven generators rectify (gen 6 was Regime 1) and the thin-cut
  pipeline proceeds; outcome B — a certificate names the throat edges to refine
  (Regime 2). Either outcome is a success vs. the greedy stall. Then single-layer
  corc as the periodic regression (T9).
- T8 port of the scratch cases into the repo test suite; periodic quotient case
  through the real `clean_spfa()` path needs a mesh fixture.

## Files touched

- `src/math/graph/fn_Graph_spfa.hpp` (new), `src/math/graph/fn_Graph_spfa.cpp` (new),
  `src/math/graph/CMakeLists.txt` (+1 line)
- `src/homology/cl_Cohomology.cpp` (`clean_spfa()` body, `cohomology::fire_node_coboundary()`,
  two includes) — on top of Christian's scaffold (`mFunClean`, `clean_greedy()` split,
  `cl_Cohomology.hpp` declarations)
- `src/homology/cl_CutData.cpp` (T7 message split; stray declaration removed)
- `todo/thin_cut_nonunit_rectification_implementation.md` (T1/T2/T4/T5/T6/T7 ticked
  with notes, Status updated)
