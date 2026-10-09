# DR-113 struck on run evidence that was never recorded; DR-122's LINE3 terminal-loop closure fixed behind a static helper, the rest of the order-2 residue stated in the module doc

**Date:** 2026-10-08
**Purpose:** Investigate whether DR-113 (comm fabric drain check, gate (a) owed) and DR-122
(order-2 cut-pipeline residue, watch row) can be cleared within the standing policies, and clear
them where they can.
**Module:** `src/comm` (read only), `src/homology` (`cl_CutFactory`), `tests/homology`,
`src/homology/doc`
**AIs:** Claude Fable 5.1 (DrClean session: investigation, plan, implementation, juries),
Codex `gpt-5.6-terra`/high and Grok `grok-4.7`/high (plan jury and code jury), Codex
`gpt-5.6-luna`/medium (language sweep of the doc section); peer sessions Yufan (local, DR-194/197
in flight, holds `cl_CutFactory.*` and the register files tonight) and Server (release audit,
nightly state)
**Status:** COMPLETE — DR-113 struck in the Yufan session's `f681f98` on this session's run
evidence; DR-122 fixed, tested and struck in this session's commit; the DR-103 guard found never
landed, recorded. Verified in a private debug MKL Blaze tree on `f681f98`: `test_homology` 49/49,
both mutation checks red where predicted, `make check` 17/17 (fast 8, mpi 4, 43 s real).
Christian's `make check` in a shared tree owed
**Exchange:** `tmp/ai_exchange/review_dr122_line3_terminal_loop.md`,
`tmp/ai_exchange/dr113_dr122_register_text.md`

## Summary

- **Coordination first.** Both peers answered before any work started. Yufan: landing DR-194 and
  DR-197 tonight, `cl_CutFactory.cpp/.hpp`, `cl_ThinShellFactory.cpp`, `cl_Vertex.hpp`, two new
  tests and the register, index and devlog-index files are theirs until their commit; three
  uncommitted todo files of theirs sit in the shared checkout and ride their commit. Server: main
  pinned at v0.9.2 and never pushed, devel is the release candidate, the nightly (4-way matrix,
  `make -j2 check`) has run `commmpi` at 2 and 4 ranks green but that is gate (b) again, not a
  chokepoint at more than one rank; every homology commit on devel costs another blind
  Codex + Grok xhigh audit of the final tip; Codex and Server both judged neither row Gregory's.
- **DR-113.** The row tied gate (a) to the Tier 2 FEM fixture, but the gate as filed on
  2026-08-31 was a parallel deck in a debug build with zero drain aborts. That run happened in the
  DR-105 weekend campaign (2026-10-02/03) and was never written against this row. Re-read at the
  source in `tmp/dr105_weekend`, not from its devlog: the campaign's debug library carries the
  `comm_drain_check` error string and the release library does not; every parallel log's banner
  says `DEBUG flags are on`; nine runs at 2, 4 and 8 ranks on the tapestack3d quench deck with
  cohomology cuts; all nineteen logs grep clean for a stray or drain message, with the positive
  control that the per-step lines are found in all nineteen; `p2_np8_newton` completed the full
  0.3 s at 8 ranks (15826 steps, 304 saves, exit 0), so the five chokepoints were each crossed
  thousands of times with the check live; `p1_np4` crossed 1661 steps before an unrelated
  floor-retry abort; `p1_np2` ran clean to its cap. That is the top rung of the register's closure
  ladder. Christian's ruling: strike (the DR-42/49 pattern); T-62 stays in
  `todo/test_gap_backlog.md` behind F1 as coverage work. Strike text in
  `tmp/ai_exchange/dr113_dr122_register_text.md`, folded into the Yufan session's commit by
  agreement.
- **DR-122.** Four items re-read in the tree. (4) is real: `close_terminal_loops` built the return
  pass of a closed terminal loop with `insert_node( tNode, n-1-k )`, which for a LINE3
  `[e0, e1, mid]` yields `[mid, e1, e0]` while the same function's arclength extension and
  node-list rebuild read node 2 as the midpoint; `MaxwellFactory::create_cuts` runs before the
  order-1 guard, so an order-2 thin-shell deck executes the wrong closure and then dies at the
  guard. (2) `CutSet::create_duplicates` sources every duplicate from its original without the
  hanging cascade `InterfaceProcessor` applies — confirmed, but its order-2-only trigger is
  DR-103's audit claim, not re-derived here. (3) no TRI6/TET10 fixture — confirmed. (1) the doc
  states the diagonal coverage is not verified face by face — as recorded. The row's own gate is
  "once the order-1 guard is lifted", a scope decision ruled out for 1.0. Christian's ruling: fix
  what is gateable (the LINE3 reversal, with a test), state the rest in the module doc, strike.
  The per-curve body of the closure became the static `CutFactory::close_terminal_loop( Curve *,
  const Cell< Node * > & aDuplicates, id_t & aMaxID )` (the precedent is the static
  `refuse_closed_side_curve_sharing_nodes`; a test cannot own a `CutFactory`, whose constructor
  dereferences a typed `Topology` map), with the slot formula `k < 2 ? 1 - k : k` and a bounds
  assert on the duplicate lookup. `tests/homology/test_CutFactory_terminal_loop.cpp` pins the
  LINE3 closure (node order, node list without the closing node, every arclength value, closed
  flag, id advance), the `gNoIndex` branch, LINE3 after `Curve::reverse`, the live LINE2 closure,
  the empty-curve no-op and three refusals. The order-2 boundary is a new "Element Order"
  subsection in `src/homology/doc/README.md`, four items.
- **The plan jury found a second defect, on the live LINE2 path.** Codex and Grok, blind and
  independently: the closure's arclength extension appended segment lengths for `k` ascending
  while the return pass is stored in reverse segment order, so every interior return-node
  parameter is wrong for unequal lengths; total and monotonicity are preserved, which is exactly
  what my pre-registered oracle checked. No solver reads `Curve::arclength()`; the BFM curve
  save/load is the only other consumer, so a `.bfm` written before the fix carries the old
  interior samples. Fixed in the helper (both branches walk `k = n-1 .. 0`), pinned with
  unequal lengths 1 and 3 and every value asserted; the ascending-order mutant fails on LINE3
  entries 5–7 and LINE2 entry 3.
- **Grok found that DR-103's fix never landed.** The closed DR-103 row says the `check_midside`
  `is_flagged` guard landed in `0db03e32` with nine call sites updated. `cl_CutProcessor.hpp:329-345`
  has no such guard, `git log -S` finds no change since the history re-init, and the 2026-08-31
  DR-109 devlog records the hunk as left in `stash@{0}`, "stashed and unclaimed"; the stash list
  is empty today. Christian's ruling: a correction note on the closed row and a fourth item in
  the order-2 boundary; re-implementation is future work behind the order-2 scope.
- **The code jury asked for one correction, applied as a separate approved step.** Codex: the
  duplicate-lookup error fired after a segment was allocated and the id counter advanced, so the
  throwing policy of the test main would leak the segment and mutate `aMaxID`. Grok: the public
  helper did not state or check its arclength-size precondition, and a null duplicate in range
  would be stored. One preflight block before the first allocation closes all three; three
  refusal tests pin that the curve and the counter survive a refusal.

## Juries

Plan round (Codex `gpt-5.6-terra`/high, Grok `grok-4.7`/high, blind): both REVISE on the
arclength order; Grok's F3 (DR-103 guard absent) confirmed in the tree; Codex's contract-before-
allocation adopted for the element type, the closed-input guard declined on Grok's L-16 argument
and the upstream refusal in the Maxwell factory. Code round (same depth): Codex one correction
(leak on throw), Grok no blocking defect, two open risks; all closed by one preflight block.
Predictions: plan round P1 raised, P2 not, P3 refuted (`assign_edges` is a consumer), P4 not;
code round P1–P3 not raised, P4 held.

## Key Findings

- A by-catch run can be the strongest evidence a row ever gets, and it is invisible unless
  someone re-reads other campaigns' logs with the row in mind. The binary-string check
  (`grep -a -c` for the error text in the shared library) is a cheap discriminator for "was the
  debug-only check compiled in", stronger than the banner line.
- The clang diagnostics the IDE raises on a new test file without the project's include paths
  are noise; the syntax gate with the build tree's own `flags.make` is the one that counts.

## Changes Made

- `src/homology/cl_CutFactory.hpp/.cpp`: static `close_terminal_loop`, `close_terminal_loops`
  reduced to the two loops and the call; the slot formula fixed; the return-pass arclength
  walked in reversed segment order; type, duplicate and arclength-size preconditions as
  `BELFEM_ERROR` before any allocation.
- `tests/homology/test_CutFactory_terminal_loop.cpp` (new), `tests/homology/CMakeLists.txt`.
- `src/homology/doc/README.md`: "Element Order" under Common Pitfalls.
- `todo/closed/dr122_line3_terminal_loop.md` (the plan, retired), `todo/README.md`,
  `todo/debt_register.md` (DR-122 removed, freeze lens recounted) and
  `todo/debt_register_closed.md` (DR-122 archived with its closing cell; correction note on
  DR-103; DR-113 archived by the Yufan commit `f681f98` on this session's text), this entry,
  `devlog/README.md`, `devlog/campaigns/release_1.0.md`.

## Open

- Christian's `make check` in a shared tree (the private-tree gate is recorded above).
- The DR-103 guard: re-implementation from the row's description, future work behind the order-2
  scope; no live row (ruling).
- Live register after this session: 5 rows (DR-29, DR-71, DR-177 homology; DR-193, DR-196).
  No row filed this session.
