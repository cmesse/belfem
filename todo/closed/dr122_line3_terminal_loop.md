# DR-122: the LINE3 terminal-loop closure puts the midside node in a corner slot

**Date:** 2026-10-08
**Purpose:** `CutFactory::close_terminal_loops` builds the return pass of a closed terminal loop by
reversing every segment with `insert_node( tNode, n-1-k )`. For a LINE3 segment `[e0, e1, mid]`
that yields `[mid, e1, e0]`: the midside lands in corner slot 0 while every consumer of the
closed loop reads node 2 as the midpoint. Fix the reversal so corners swap and the midside stays
at slot 2, as `Curve::reverse` and the curve factory already do, and pin the closure of a LINE2
and a LINE3 loop with a test. The order-2 residue that cannot be gated (the missing hanging-source
cascade in `CutSet::create_duplicates`, the unverified diagonal-case midside coverage, the absent
TRI6/TET10 fixture) moves into the module documentation as the stated order-2 boundary, and the
register row is to be struck on Christian's ruling of 2026-10-08 once the gate has run.
**Module:** `src/homology` (`cl_CutFactory`), `tests/homology`, `src/homology/doc`
**AIs involved:** Claude Fable 5.1 (investigation + plan), Codex `gpt-5.6-terra`/high and Grok
`grok-4.7`/high (plan jury and code jury)
**Status:** ✅ COMPLETE 2026-10-08 — landed on devel: static `CutFactory::close_terminal_loop` (slot formula `k < 2 ? 1 - k : k`, return-pass arclength walked `k = n-1 .. 0`, type/duplicate/arclength-size preconditions as `BELFEM_ERROR` before any allocation), `tests/homology/test_CutFactory_terminal_loop.cpp` with eight cases, the order-2 boundary (four items, incl. the DR-103 guard that never landed) in `src/homology/doc/README.md`. Plan jury and code jury Codex terra/high + Grok 4.7/high, reconciled. Verified in a private debug MKL Blaze tree on `f681f98`: `test_homology` 49/49, both mutation checks red where predicted, `make check` 17/17 (fast 8, mpi 4, 43 s). Residual: Christian's `make check` in a shared tree; order-2 h-φ itself stays behind the guard (not 1.0).

> **Scope guards:**
> - Order-2 h-φ stays behind the deliberate order-1 guard in
>   `MaxwellFactory::create_hanging_edges_and_facets` (`cl_MaxwellFactory.cpp:1966`). Lifting it is
>   a scope decision (DR-103 ruling: not a 1.0 feature) and is OUT of scope here.
> - `CutSet::create_duplicates` is not touched: its missing cascade has no constructible trigger
>   without the order-2 Maxwell pipeline, so a fix there could not be gated. It is documented.
> - The LINE2 path is live in every 2-D tape deck. Its behaviour must not change; the test pins it.
> - `cl_CutFactory.cpp/.hpp` are held by the Yufan session tonight (DR-194/DR-197 destructor and
>   restore). Implementation runs in a private worktree and lands only after that commit.
> - The cohomology core (protocol §7.1) is not touched; `cl_CutFactory` is outside the ban.

---

## 1. Current Behaviour and How It Fails

`CutFactory::close_terminal_loops` (`cl_CutFactory.cpp:1576-1706`), per terminal curve of every
protoshell: for each original segment it creates a duplicate element of the same type, fills it
with the thin-shell duplicate of each node (or the node itself when its index is `gNoIndex`)
at slot `n-1-k` (`:1596-1605`), stores the duplicates in reverse order, extends the arclength
vector (`:1611-1662`), rebuilds the node list from `node( 0 )` of every segment and, for LINE3,
`node( 2 )` as "the center node" (`:1668-1697`), and sets the closed flag.

| Failure | Mechanism | Evidence |
|---|---|---|
| LINE3 duplicate segment has the midside at slot 0 | `n = 3`, `k = 2` is the midside, `n-1-k = 0` | `cl_CutFactory.cpp:1602`; a segment is built `[corner, corner, mid]` by `CurveFactory::create_segments` (`cl_CurveFactory.cpp:1009-1026`) |
| Node list of the closed loop lists a corner as the center | `tNodes( tCount++ ) = tSegment->node( 2 )` on the reversed half | `cl_CutFactory.cpp:1686` |
| Arclength bookkeeping assumes node 2 is the midpoint | midpoint arclength written at `s + 0.5 l` for the reversed half | `cl_CutFactory.cpp:1648-1656` |

Reachability: `MaxwellFactory::create_cuts` (`cl_MaxwellFactory.cpp:675`) runs before the order-1
guard (`:1966`); `CutFactory::create_thin_shell_cuts` calls `close_terminal_loops` (`:1166`). On an
order-2 thin-shell deck the wrong loop is built and the run then aborts at the guard, so no
completed h-φ solve observes the defect today; `CutFactory::run` does walk the wrong segments in
`assign_edges` (`cl_CutFactory.cpp:135-138`, `cl_Curve.cpp:166-174`) before the guard, and in 2-D
`orient_terminal_curves_2D` may `reverse()` them (`:1909`). The two in-tree reversals that handle LINE3 keep the
midside at slot 2: `Curve::reverse` (`cl_Curve.cpp:136-141`) and the curve factory's end-node sort
(`cl_CurveFactory.cpp:1218-1221`, "only corners 0 and 1 swap: a LINE3 midpoint stays node 2").
Confidence high: all three sites re-read tonight.

**Bottom line:** one slot formula is wrong for `n = 3`; everything downstream in the same function
assumes the correct one.

## 2. Mechanism

Replace the slot formula with the corner swap the rest of the mesh layer uses:

```cpp
// corners swap, a LINE3 midside stays node 2 ( Curve::reverse )
tDup->element()->insert_node( tNode, k < 2 ? 1 - k : k );
```

To make the closure testable without a `CutFactory` instance — its constructor dereferences
`aTopology->groups( DomainType::ThinShell )` (`cl_CutFactory.cpp:51`), and `Map::operator()`
raises `BELFEM_ERROR` on a missing key (`cl_Map.hpp:227-233`), so a test would need a populated
`Topology` — the per-curve body of the loop becomes a static member, following the precedent of
`refuse_closed_side_curve_sharing_nodes` (`cl_CutFactory.hpp:213`):

```cpp
/**
 * Closes one terminal curve into a loop: every segment is duplicated in
 * reverse, its nodes replaced by aDuplicates( index ) where the node has an
 * index, the arclength vector extended, the node list rebuilt and the closed
 * flag set. aMaxID is advanced by one per new segment element. LINE2 and
 * LINE3 only; a LINE3 midside stays node 2 of the reversed segment.
 */
static void
close_terminal_loop( Curve * aCurve, const Cell< Node * > & aDuplicates, id_t & aMaxID );
```

`close_terminal_loops` keeps its loop over protoshells and terminal curves and calls the helper
with `mThinShellDuplicates` and `mMaxID`. The body moves verbatim except for the slot formula.

Rejected: fixing the formula in place with no test (the register's own gate for this row is an
order-2 deck behind the guard, so an untested fix would be reviewed, never verified); building a
`Topology` in the test (the type map is filled by `Topology::run` from typed sidesets, which the
hand-built mesh does not have — larger fixture, same assertion).

## 3. Gap Table

| # | State | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| 1 | LINE3 reversed segment node order | any order-2 terminal loop | no | (c) | `cl_CutFactory.cpp:1602` |
| 2 | LINE2 reversed segment node order | every 2-D tape deck | yes, `[e1, e0]` | (a) | `n = 2`: `n-1-k` and `1-k` coincide |
| 3 | Arclength extension, both types | loop parameterization | **no** — lengths appended for `k = 0..n-1` while the return pass is stored reversed (`:1605`); wrong interior values for unequal lengths, total and monotonicity preserved | (c) | `:1624-1628`, `:1645-1652`; found by Codex and Grok, plan round |
| 4 | Node list rebuild | `orient_terminal_curves_2D` and consumers | yes, given item 1 | (a) | `:1668-1697`, unchanged |
| 5 | Ownership of the new segments | no leak in test and production | yes | (a) | `Curve::~Curve` deletes segments (`cl_Curve.cpp:45-51`), `Segment::~Segment` deletes its element (`cl_Segment.cpp:26-29`) |
| 6 | Hanging-source cascade in `CutSet::create_duplicates` | order-2 cuts through hanging originals | no | (b) → doc | `cl_CutSet.cpp:80`; precedent `InterfaceProcessor` `cl_InterfaceProcessor.cpp:612-630`; no constructible trigger behind the guard |
| 7 | Diagonal thick-cut midside coverage | order-2 thick cuts | not verified face by face | (b) → doc | `src/homology/doc/thick_thin_cuts_and_conjugate_edges.md:173` |
| 8 | TRI6/TET10 homology fixture | order-2 cut pipeline | none (five test files in `tests/homology/`) | (b) → doc | order-2 scope |
| 9 | DR-103's `check_midside` `is_flagged` guard | diagonal TET10/TRI6 thick cuts at order 2 | **no** — `cl_CutProcessor.hpp:329-345` has no guard; the closed DR-103 row says it landed, the 2026-08-31 devlog says it stayed in `stash@{0}` | (b) → doc + register correction, **Christian's ruling** | Grok F3, confirmed in tree |

### 3.1 Cross-cutting

The helper takes `aDuplicates` by index. In production `mThinShellDuplicates` is sized to the
thin-shell node count (`cl_CutFactory.cpp:1316`) and a terminal node with an index is a thin-shell
node, so the lookup is in range; nodes with `gNoIndex` bypass it. The test exercises both branches.

## 4. Ordered Steps

- [x] **R1** (done 2026-10-08; reconciled in the exchange; Christian approved the helper, then the widened scope and the DR-103 correction) Plan jury (Codex terra/high + Grok 4.7/high), reconcile, Christian approves plan and edits.
- [x] **R2** (done 2026-10-08 in a private worktree on `0b25a66`, to be rebased on Yufan's tip; helper approved by Christian before the jury returned) Private worktree on that tip. Add the
      static `close_terminal_loop` to `cl_CutFactory.hpp` (public, beside the other static helper)
      and move the per-curve body into it in `cl_CutFactory.cpp` with the corrected slot formula;
      `close_terminal_loops` becomes the two loops plus the call.
- [x] **R2b** (approved and done 2026-10-08; worktree rebased on `f681f98`) Arclength extension iterates `k = n-1 .. 0` in both branches so the appended lengths follow the reversed return pass; the element-type `BELFEM_ERROR` moves ahead of the first allocation; the duplicate-lookup check becomes `BELFEM_ERROR`.
- [x] **R3** (done 2026-10-08; four cases, 44/44 in the private tree, mutation check run by hand: old formula red on both LINE3 cases, LINE2 and empty-curve green) `tests/homology/test_CutFactory_terminal_loop.cpp`, registered in
      `tests/homology/CMakeLists.txt`. Fixture: nodes owned by the test, a `Curve` (owns its
      segments) built like `CurveFactory::create_segments` does, open, with arclength filled.
      Cases:
      - LINE3, two segments `[1,2,m12]`, `[2,3,m23]`, nodes indexed 0..4 with a duplicates cell of
        five fresh nodes: after the call 4 segments; segment 2 is `[d3, d2, dm23]`, segment 3 is
        `[d2, d1, dm12]`; node list `[1, m12, 2, m23, d3, dm23, d2, dm12]`; arclength 9 entries,
        strictly increasing, last equal to twice the original length; closed flag set; `aMaxID`
        advanced by 2.
      - LINE3 with every node at `gNoIndex` and an empty duplicates cell: the reversed segments
        reuse the originals (`[3, 2, m23]`, `[2, 1, m12]`).
      - LINE2, two segments: 4 segments, reversed pairs, node list of 4, arclength 5 entries.
      - (R2b) Unequal segment lengths (1 and 3) for LINE2 and LINE3, every arclength value asserted; a `Curve::reverse()` after the close keeps the LINE3 midside at node 2; an empty curve is left untouched.
      - Mutation check (run once by hand, not committed): the old formula must turn the LINE3
        case red on the node-2 assertion and leave LINE2 green.
- [x] **R4** (code jury reconciled 2026-10-08: preflight fix applied as a separate approved step, re-gated) (syntax gate green on both files with the build tree's `flags.make`; own-files diff at `tmp/ai_exchange/review_dr122_line3_terminal_loop.diff`; code jury dispatched 2026-10-08 22:50, terra/high + grok-4.7/high) Syntax gate with the real flag set; `git diff HEAD -- <own files>` to
      `tmp/ai_exchange/review_dr122_line3_terminal_loop.diff`; code jury, same slug; verification
      pass; reconciliation.
- [x] **R5** (private tree on `f681f98`, debug MKL Blaze: `make check` 17/17, fast 8, mpi 4, 43 s real; `test_homology` 46/46; both mutation checks red where predicted; Christian's shared-tree `make check` still owed) `make check` in my private scratchpad
      tree (homology suite plus the rest); report N/N. Christian's `make check` in a shared tree
      is the second gate.
- [x] **R6** (done 2026-10-08, Codex luna/medium sweep applied, no suspected errors) Text half, no jury: the order-2 boundary (items 6–8) stated in
      `src/homology/doc/README.md` "Development Notes" and cross-referenced from
      `thick_thin_cuts_and_conjugate_edges.md`, facts inline, no register IDs; Codex language
      sweep (luna/medium) over the new prose.
- [x] **R7** (done 2026-10-08) Close out: register row DR-122 struck with the status cell naming the fix, the test,
      the documented boundary and the gate that ran; DR-113's status cell rewritten on the DR-105
      run evidence (separate ruling, same session); `devlog/dl20261008_dr113_dr122_clearance.md`,
      `devlog/README.md`, campaign page, this plan to `todo/closed/`. Register and index edits go
      through the Yufan session's commit, which holds those files tonight. Server receives the
      SHA and file list for the release audit.

## 5. Open Design Questions

- **O1** Static helper with the whole per-curve body (recommended) versus the one-token fix with a
  helper-only slot test. The former pins the live LINE2 closure as a regression test for free;
  the latter is a smaller diff. Jury input wanted; decision Christian's.
- **O2** RESOLVED 2026-10-08 → no closed-input guard (Grok: deck terminals are refused upstream, `cl_MaxwellFactory.cpp:3643-3648`; a new guard is a new failure mode, L-16); the element-type error moves before allocation (Codex). Whether the helper should refuse a closed input curve (`is_closed()` already true) with a
  `BELFEM_ERROR`. Today the caller skips only empty segment lists (`:1586`). Setup code, once per
  run, so the error tier is allowed; not part of the defect.

## 6. Interface

See §2. No input key, no deck, no `input_schema.yaml` change.

## 7. Definition-of-Done

- [x] R1–R7 ticked, Status line rewritten with the gate that ran and the suite count.
- [x] DR-122 struck; DR-113 rewritten and struck (DR-113 via `f681f98`); both cells name their evidence and residuals.
- [x] Homology doc states the order-2 boundary; language sweep applied.
- [x] Server told the SHA and files (after the push).

- **O3** The DR-103 guard (gap row 9): recover it from `stash@{0}` if it still exists, or re-file it; and correct the closed DR-103 row. Not this plan's code; Christian's ruling.

## 7a. DR-113 evidence cited for the ruling

`tmp/dr105_weekend/bin/dbg/libbelfem.so.0.9.3` carries the `comm_drain_check` error string (`grep -a -c`: 1; `bin/rel`: 0); `*/belfem.log` banners `DEBUG flags are on`; nine runs at 2, 4, 8 ranks; zero logs with a stray or drain message of nineteen; `p2_np8_newton/belfem.log` 15826 steps, 304 saves, `exit.txt` 0. Christian's ruling 2026-10-08: strike, T-62 stays in the backlog.

## 8. Audit Trail

- `tmp/ai_exchange/review_dr122_line3_terminal_loop.md` — plan jury, code jury, reconciliation.
