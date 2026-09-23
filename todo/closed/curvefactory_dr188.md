# CurveFactory: Quadratic Midpoints, 3-D Pairs, Mixed Order, Closed Curves (DR-188)

**Date:** 2026-09-22
**Purpose:** Close DR-188's four items in `mesh::CurveFactory`. The mechanism: build every curve's
segment pairs and LINE3 midpoints from facet edges instead of the node-clique graph. Then order the
segments by traversal instead of re-walking an adjacency, so that closed curves close.
**Module:** `src/mesh` (`cl_CurveFactory.{hpp,cpp}`), `tests/mesh/test_CurveFactory.cpp`
**AIs involved:** Claude (Opus 5.5, plan and implementation); Codex gpt-6-astra + Grok 4.7 (plan
jury, `tmp/ai_exchange/review_curvefactory_dr188.md`)
**Status:** CLOSED 2026-09-23 — landed in the combined commit with the mesh bundle and DR-181. Plan jury, three code rounds; `CurveFactory.*` 24/24 (14 new cases red on the old library), the mandatory A/B 10/10 decks identical, `ThinShellTerminalCurves.*` 6/6 (3 red without the guard). DR-188 stays open as `[W]` for its residue (a)–(f). Devlog `dl20260922_debt_register_oneliners_dr181_dr188.md`.

> **Baseline — read before citing a line.** Every `cl_CurveFactory.*:N` below is from the **working
> tree**, not `HEAD`. `src/mesh/cl_CurveFactory.{hpp,cpp}` and `tests/mesh/test_CurveFactory.cpp`
> carry the uncommitted mesh one-liner bundle (`todo/mesh_oneliners_dr161_169_176_179_180.md`,
> jury-closed, `test_mesh` green on Christian's release tree, commit pending). The bundle's helper
> `create_pairs_from_facets` (`cl_CurveFactory.cpp:890-986`), the `unique_ptr` in
> `from_2d_sidesets` (`:700`), the `element()->insert_node( …, 2 )` midpoint write (`:1014-1022`)
> and the three quadratic/corner tests exist **only** in the working tree. This plan must land
> after that commit, or be rebased onto it. Re-read every anchor before editing (L-05, L-08).
>
> **Scope guards.**
> - In scope: DR-188 (D), (Z1), (Z2), (Z3) as filed and widened on 2026-09-22.
> - Out of scope: the cohomology core (`cl_Homology.cpp` is read here as a consumer only); the
>   downstream closed-curve support in `CutFactory` / `ThinShellFactory` (see O3); the
>   `sort_end_nodes` throw-mode index residue already recorded in the bundle (§4.2 G2).
> - Compatibility promise, narrowed by the plan jury: an open curve comes out bit-identical — the
>   same pairs in the same key order, the same segment element ids, the same node order, the same
>   arclength — **when the edges its sidesets share form one path**. The old clique rule can also
>   bridge two shared-edge components with a co-element chord; such a curve works today and is
>   refused as "not one connected chain" after R2. Whether a shipped deck does that is what the
>   mandatory A/B gate (§4.1) settles.

---

## 1. Current Behaviour and How It Fails

| # | Failure | Mechanism | Evidence (working tree) |
|---|---|---|---|
| D | **Silent midpoint drift** on a LINE3 curve | `collect_midnodes( aPairs, aMidNodes )` writes `aMidNodes( tCount++ )` only inside the match (`cl_CurveFactory.cpp:856-885`, write at `:877`). One unmatched pair shifts every later midpoint one slot early and leaves the last slot `nullptr`. `create_segments` then assigns them positionally (`:1014-1022`). Callers: `thin_shell_side_curves` (`:122-126`) and the 3-id `intersect` (`:565-570`). | source trace (Grok round 2 and Claude, `review_mesh_oneliners…md` "Open Risks"). Not executed |
| D′ | **The 3-D path has the clique defect too.** The pairs for `intersect` come from `create_pairs` (`:792-839`). It pairs every two end nodes that are graph neighbours, and `connect_nodes_to_nodes` links every node of every incident element (`cl_Mesh_ConnectivityCalculator.cpp:343-351`). Wherever one element holds 3+ end nodes, the clique adds a non-edge chord. | (i) **Curve corner inside one hex.** A = top face, B = two adjacent side faces. The top face diagonal `5–7` is added, and its three nodes form a triangle that `sort_end_nodes` **accepts as a closed loop** (all degrees 2, `:1048-1094`). On LINE3 the diagonal is not a hex edge, so `collect_midnodes` misses it and the drift in D hands `create_segments` a null midpoint. (ii) **Ear facet.** Any tet or prism face with three end nodes (a tet, or a sideset-A triangle whose three corners all lie on B) adds its third edge. On a longer curve that is a degree-3 node, and `sort_end_nodes` refuses it as a "branch" (`:1076`): loud. On a 3-node curve it is again a silent triangle loop. | source trace (Claude), fixtures in R2 |
| Z1 | **Mixed order is silent or out of range** | `from_2d_sidesets` concatenates every sideset's facets (`:681-697`) but takes the curve order from `aSideSets(0)` alone (`:700` → `Curve` ctor → `SideSet::element_type()`, which reads facet 0 only, `cl_SideSet.hpp:271-282`). `create_pairs_from_facets` takes `tNodes( 2 )` whenever `size() > 2` (`:932`). Consequences: a LINE4 edge supplies its ξ = −1/3 node as a "midpoint"; a LINE3 sideset after a LINE2 one silently loses its midsides; a non-line facet reaches the pair key with a stale node index and indexes `aEndNodes` out of range. The reverse order, LINE3 first and then LINE2, is already refused (`:947-954`). | source trace (astra + Grok, round 3, confirmed) |
| Z2 | **A closed curve is mis-closed** (LINE2 and LINE3, 2-D and 3-D) | Traced in Appendix A. `orient_segments` (`:1195-1212`) flips every segment to low→high **including the closing one**, so no ordering of the segments is a head-to-tail cycle. `sort_segments` (`:1216-1241`) starts at the walk-index-0 node, from its adjacency slot 0 (`:1219`); on a closed curve that slot holds either of the node's two segments, depending on the order of the segment cell. Its first `next()` step then picks the smaller-id candidate between the two ends (`:1397-1400`), so the traversal direction is arbitrary. `collect_nodes` then trips its closure `BELFEM_ASSERT` (`:1445-1451`; it throws in debug under `throw_on_error`, `assert.hpp:190-193`). In release it returns a wrong node list, e.g. `1,2,2,3` for a triangle. | source trace, plus a model check: every one of 36 + 576 + 14 400 enumerated 3-, 4- and 5-loops fails (Appendix A). `SingleLoopIsAccepted` (`test_CurveFactory.cpp:354-381`) stops at `sort_end_nodes` |
| Z3 | **Leak in an uncalled overload** | `intersect( aThinShellSideSet, aBoundarySideSet, DynamicBitset & )` (`:597-674`) does `new Curve` (`:601-604`) and leaks it on the `nullptr` return (`:622-625`) and on any refusal. | source; no caller (§3.1 C5) |

**Bottom line:** every defect in DR-188 traces to one of two design faults. First, curve pairs and
midpoints are recovered from **mesh-wide node and element adjacency** rather than from the facets
the curve is made of. Second, the segment order is recomputed by a **second, independent walk**
rather than taken from the node walk `sort_end_nodes` has already done. The bundle removed the first
fault for the 2-D path. This plan removes it for the 3-D and thin-shell paths and removes the second
fault.

## 2. Architecture: Why Facet-Built Pairs and Traversal-Keyed Segments Are the Right Spine

**Pairs.** `create_pairs_from_facets` (`:890-986`) is jury-closed and green in `test_mesh`. It emits
one pair per facet edge whose two corners are end nodes, keyed and ordered exactly like
`create_pairs` (`:817-819` vs `:930`), and gathers each pair's midpoint in the same pass, refusing
a missing or conflicting midpoint. Routing the 3-D `intersect` through it gives D and D′ one
mechanism. For an intersection, the pair set is **the edges that a facet of A and a facet of B
share**, not merely "edges of A facets with both corners on B" as the brief proposed. The tet ear
(R2 fixture T2) refutes the A-only rule: A's triangle `0–2–1` has all three corners on B, so the
A-only rule keeps the chord `0–2`, which is an edge of no B facet. On a conforming mesh, the edges
common to A and B are exactly the discrete curve. Confidence **high** (source trace plus the fixture
enumeration in R2).

*Rejected — option (a), a `BELFEM_ERROR` on an unmatched pair in the old `collect_midnodes`.* It
turns the drift into a refusal but leaves D′ in place. The corner case would then be refused on
LINE3 and silently closed on LINE2 (a triangle), and the ear would still be a branch refusal. The
refusal is kept, inside the facet-built route (the helper already refuses a missing midpoint).

**Segment order.** After `sort_end_nodes`, every end node carries its walk index 0…n−1
(`:1175-1181`), consecutive along the curve (`SingleChainIsMonotone`, `SingleLoopIsAccepted`).
Orient each segment **along the walk**: low→high, except that on a closed curve the closing segment
{0, n−1} goes n−1→0. Then segment *k* is the one whose `node( 0 )` has index *k*, and
`sort_segments` becomes "set each segment's index to its first node's index, then sort". There is no
second walk, no first-step ambiguity, and no segment flags. For open chains this is the same result
the old walk gives: on an open chain the start is a tip, so the walk is forced. Appendix A checks
2000 random open chains in the model, all identical.

*Rejected — patch only the closing segment's orientation (the DR row's first suggestion).*
Appendix A shows that this alone is not enough: `sort_segments`' first step can still turn toward
the closing segment, so the traversal can run against the orientation. *Rejected — "let
`collect_nodes` close on `tSegments( 0 )->node( 0 )`" (the DR row's second suggestion).* That would
hide the misclosure in the node list and leave the segments non-consecutive and mis-oriented, and
`Curve::reverse`, `EdgeFactory::link_curves_to_edges` and the terminal orientation (§3.1 C4) all
read the segments.

## 3. Gap Table

| # | State / behaviour | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| G1 | 3-D pairs = edges shared by an A facet and a B facet | D′, D | no — clique | (c) → R2 | `:559`, `:792-839`; `cl_Mesh_ConnectivityCalculator.cpp:343-351` |
| G2 | 3-D LINE3 midpoint per pair | D | drift-prone | (c) → R2 | `:565-570`, `:856-885` |
| G3 | Thin-shell side-curve LINE3 midpoint per pair | D | drift-prone, and depends on node→element links | (c) → R3 | `:113-126`; pairs are already facet edges (`:114-120`, edges from `:195-250`) |
| G4 | `create_pairs`, old `collect_midnodes` | — | callers vanish after R2/R3 | (a) delete → R3 | private, `hpp:116-121` |
| G5 | 2-D curve order uniform and LINE2/LINE3 | Z1 | no | (c) → R5 | `:681-700`, `:932` |
| G6 | Empty curve (no facets / no shared edge) | robustness | two failures: with **no end nodes**, the degree check passes (`0 == 0`, `:1076`) and `sort_end_nodes` indexes `aEndNodes( 0 )` (`:1099`); with **end nodes but no pair**, every degree is 0 and the branch refusal fires (`:1076`), loud but misnamed. An empty id vector indexes `aSideSets( 0 )` | (b) → O6 | pre-existing |
| G7 | Closing segment oriented along the walk | Z2 | no | (c) → R4 | `:1200-1211` |
| G8 | Segment order from the node walk | Z2 | no — a second walk | (c) → R4 | `:1216-1241`, `:1321-1402` |
| G9 | Closed-curve node list repeats the start (`n+1` entries), arclength `n+1` | Z2 | intended by `collect_nodes` (`:1444-1451`), never achieved | (a) kept | consumers already tolerate either convention (`cl_ThinShellFactory.cpp:2152-2158`) |
| G10 | LINE3 midpoint under a flip | Z2 | yes — only nodes 0/1 swap (`:1204-1209`); `compute_coordinates`' `|aξ+b|` is even under x0↔x1 over symmetric Gauss points (`:1532-1547`) | (a) | Appendix A |
| G11 | Bitset `intersect` overload | Z3 | leaks, no caller | (a) delete → R1 | `:597-674`, `hpp:60-65` |
| G12 | `intersect` throw safety | D, R2 | raw `new Curve` (`:546-549`) | (c) → R2 | same `unique_ptr` pattern as `:700` |
| G13 | Downstream consumers of *closed* CurveFactory curves | Z2 follow-through | never exercised | (b) → O3 | §3.1 C4 |

### 3.1 Cross-cutting findings

- **C1 — the DR row's text is stale on one point.** It names "the facet-keyed `collect_midnodes(
  aFacets, aNumEndNodes, … )` overload the DR-180 fix pass added". Fix pass 2 (G1) removed that
  overload; its successor is `create_pairs_from_facets` (`hpp:123-138`). The row should be updated
  at close-out.
- **C2 — the brief's A-only pair rule is insufficient** (§2, fixture T2). The helper needs a
  "shared with B" filter. The change is additive and does not touch the 2-D path. See O1.
- **C3 — shipped decks.** `examples/` holds **14** decks with a `curves` section, enumerated with
  an awk scan of every `examples/*/input.conf`.
  - 2-D (list syntax → `from_2d_sidesets`): `tapestack2d_layered`, `tapestack2d_racetrack`,
    `undulator2d`, `gantry`.
  - 3-D (`@` → 3-id `intersect`): `rutherford` (untracked, new), `tapestack3d`,
    `tapestack3d_quench`, `pancake`, `tapestack_circuit`, `corc_periodic_bc`,
    `tape_quench_usermat`, `corc_twolayer`, `corc_garber_effect`, `racetrack_usermat`.

  **None is second order.** `rg -i "ElementOrder|order\s*2|-order|SetOrder" examples` returns no
  mesh-order hit (only prose: Runge-Kutta, tape-order). That scan covers the `.geo` files and the
  three Python generators, and Gmsh's default order is 1. Confidence **high**. **None is closed.**
  Each 3-D curve is a tape surface meeting a boundary or terminal plane, which gives an open arc or
  segment (the deck comments say so: "tape 1, left", "inner end"). Each 2-D curve is a
  straight tape trace (bundle §4.2 consumer check), and thin-shell side curves of tapes whose ends lie on
  domain-boundary sidesets are open: `select_edges` (`:305-353`) keeps the whole rim when no
  boundary unflags its nodes (Grok plan-jury finding 8), and R3 now refuses a closed rim.
  Confidence **medium**: geometric reading, not a run. **Consequence for D:** a deck that runs today has a clique pair set that passes
  `sort_end_nodes`. That set contains the A∩B edge set; the two are equal unless the clique adds a
  tip-to-tip chord (which closes the curve) or a chord that bridges two A∩B components (which the
  plan jury found; the facet rule then refuses the curve). The A/B gate (§4.1) decides. With the same set, `collect_end_nodes_from_intersection`
  (unchanged) gives the same slots, so the keys, the element ids (`++mMaxElementID` in pair order,
  `:1003`) and the output are unchanged.
- **C4 — consumers of node order, closure and segment orientation.** Enumerated with `rg` over
  `src/` for `nodes()`, `segments()`, `edges()`, `arclength()`, `is_closed()`,
  `set_closed_flag` and `reverse()` on curves.

  | Consumer | Reads | Open-curve impact | Closed-curve impact |
  |---|---|---|---|
  | `Curve::reverse` (`cl_Curve.cpp:126-141`) | segment order + orientation | none | now reverses a consistent cycle |
  | `EdgeFactory::link_curves_to_edges` (`cl_EdgeFactory.cpp:371-393`) | node pair key | none (order-free) | none |
  | `Homology` terminal chain (`cl_Homology.cpp:267-285`, read-only core) | `segment->edge_direction()` | none | consistent directions instead of mixed |
  | `CutFactory` tip count (`cl_CutFactory.cpp:1221-1248`) | node set | none | the repeated start node is counted twice. Harmless when it has 2 flagged facets; see O3 |
  | `CutFactory` side-curve open/closed (`:1256-1305`) | node set; **overwrites** `is_closed` | none | none |
  | `CutFactory::close_terminal_loops` (`:1576-1700`) | segments, arclength | none | would double a curve that is already closed. Its curves are the deck's `curves`, so it was reachable before in release (mis-closed) and aborted in debug; residue R-3 |
  | `CutFactory::orient_terminal_curves` 3-D/2-D (`:1705-1760`, `:1810-1911`) | `segments()( 0 )` tangent | none | first segment now leaves node 0 along the walk |
  | `ThinShellFactory` side nodes, binormals (`cl_ThinShellFactory.cpp:327-349`, `:489-575`, `:2137-2230`) | node order, `nodes()( i )` vs edges | none | already documents "may or may not repeat the start" (`:2152-2158`) and refuses closed-loop tapes (`:2194-2201`) |
  | `Mesh_BfmFile` save/load (`cl_Mesh_BfmFile.cpp:2459-2530`, `:2580-2625`) | nodes, segments, closed flag, arclength | none | round-trips the corrected list |
  | `EdgeCutter` (`cl_EdgeCutter.cpp:106-300`) | node order | none | none — no reference outside its own file (`rg -n EdgeCutter src` → only its definition) |

  **Segment flags:** the old `sort_segments` leaves every visited segment flagged (`:1221`,
  `:1235`) and nobody clears them; the new one touches no flag. No consumer reads a `Segment` flag
  (`rg "tSegment->(is_flagged|unflag|flag)\("` outside `cl_CurveFactory.cpp` → no hits).
- **C5 — Z3 caller scope.** `rg -uu -n "intersect\s*\("` across the repository, excluding the
  build trees, `tmp/`, `devlog/`, `todo/` and `literature/`, hits the two definitions, the two
  declarations and `cl_MaxwellFactory.cpp:613` (the 3-id overload). **`nonfree/` is not present
  in this checkout**, so a `nonfree` caller cannot be excluded; see O5. Git history begins at the
  release re-init (`f81dd84`), so it cannot show a former caller.
- **C6 — by-catch, outside DR-188 (report, do not fix here).**
  `CutFactory::close_terminal_loops` reverses a duplicate segment with
  `tDup->element()->insert_node( tNode, n-1-k )` (`cl_CutFactory.cpp:1596-1600`). For LINE3
  (n = 3) that puts the **midpoint in slot 0** and the original node 0 in slot 2, where a LINE3
  flip must swap slots 0/1 and keep slot 2. LINE2 is correct. It is unreachable in shipped decks
  (first order). `CutFactory` is not one of the six sealed units. See O4.

## 4. Ordered Steps

Every step lands with its tests in `tests/mesh/test_CurveFactory.cpp`, which is already wired in
`tests/mesh/CMakeLists.txt:16`, so no CMake change is needed. The binary is `test/test_mesh` in the
build tree. The mesh main forces `throw_on_error` (`tests/mesh/test_mesh_main.cpp:22`), so a
`BELFEM_ERROR` or a debug `BELFEM_ASSERT` is an observable `std::runtime_error`, which gtest records
as a failure rather than a crash. "Clean red" below means the test fails by expectation or by
exception on the current working tree. "Crash red" means a SEGV, which must run under its own
filter.

### R1 — Z3: delete the bitset `intersect` overload

- [x] Landed 2026-09-22 (O5 ruled: `nonfree/` has no caller).

- **Diff:** the overload and its declaration are gone. The `create_pairs_from_facets` contract no
  longer names it; the sentence about `intersect()` resetting every node index to `gNoIndex` first
  (`collect_end_nodes_from_intersection`) landed together with R2, which made it true.
- **Why deletion, not repair:** no caller (C5). Its one distinct idea, a thin-shell ∩ boundary
  curve through a node bitset, is served in production by `thin_shell_side_curves`
  (`cl_CutFactory.cpp:1961-1972`).
- **Consumer impact:** none in `src/`, `tests/` or `examples/`.
- **Gate:** compile; `CurveFactory.*` green. No red is constructible for a deletion.

### R2 — D + D′: the 3-id `intersect` builds its pairs from shared facet edges (after: R1)

- [x] Landed 2026-09-22 (O1 ruled: the facet-built rule).

- **Diff, helper:** `create_pairs_from_facets` takes a trailing
  `const Cell< Facet * > * aSharedWith = nullptr`. When it is given, a first pass over its facet
  edges between two end nodes fills a `Map< key_t, Node * >` of their midpoints, and **checks every
  insertion**: a key seen twice with different midpoints, or a quadratic edge with none, is refused
  with the same messages as on the A side, after `mMesh->update_node_indices()`. The main pass
  skips every key that map lacks, and refuses an A-side midpoint that differs from B's. The two
  refusals are two small lambdas shared by both passes. The 2-D call is unchanged (no
  `aSharedWith`).
- **Diff, caller:** `std::unique_ptr< Curve > tOwner`, local `Curve * tCurve`,
  `unflag_all_nodes()` on A moved above the helper (the helper reads no flag), the helper called
  with `& tCurve->sideset_b()->facets()`, Steps 3–4 (`create_pairs`, `collect_midnodes`) gone, and
  `return tOwner.release()`.
- **O6:** zero pairs are refused after `update_node_indices()`: "Curve %lu: sidesets %lu and %lu
  share no facet edge, so the curve would be empty". This covers both G6 failures on this path.
- **Precondition:** `collect_end_nodes_from_intersection` sets every mesh node to `gNoIndex` before
  it assigns the end-node slots, so every non-end facet corner is `gNoIndex`. 3-D facet elements are
  QUAD4, TRI3, QUAD8 (TRI6, QUAD9 also fit): `get_nodes_of_edge` returns `( c0, c1[, mid] )`.
- **Open-curve identity:** holds when the shared edges form one path (see the scope guard).
  The bridge case of the plan jury (A∩B = two components, joined in the clique by a co-element
  chord) is refused after R2 where it passed before. The A/B gate (§4.1) is **mandatory** before
  acceptance.
- **Tests** (`Fixture` in `test_CurveFactory.cpp`: `Mesh( 3, 0, true )`, hand-pushed nodes, one
  block, facet elements QUAD4 / TRI3 / QUAD8 linked with `set_master( element, k )`, `finalize()`;
  a LINE facet cannot hold a 3-D face):
  - [x] **T1 `ThreeDCornerInOneLinearHexGivesTwoSegments`** — HEX8, A = {5}, B = {0, 3}; expects
    the open curve 5–4–7. Red on the pre-fix tree: clean (the clique triangle `5–7` closes it;
    debug throws from `collect_nodes`, release fails the counts).
  - [x] **T2 `ThreeDEarFacetOfOneTetGivesTwoSegments`** — TET4, A = {3}, B = {0, 1}; expects
    0–1–2. Red: clean (clique triangle). Discriminates O1: the A-only rule also gives the triangle.
  - [x] **T3 `ThreeDCornerInOneQuadraticHexGivesTwoMidpoints`** — T1 on HEX20 with QUAD8 facets;
    required; checks midpoints 16 and 19 and the arclength. Red: **crash** (null midpoint from
    the drift, dereferenced in `collect_nodes`). Own filter.
  - [x] **`ThreeDSharedEdgeMidpointsMustAgreeInEitherOrder`** — two HEX20 sharing a face's corners
    but not its midpoints; B holds that face from both elements, in both insertion orders
    (sidesets 2 and 3). Both refuse "different midpoints". The second order is the one where a
    keep-last map would agree with A and hide the conflict. Red: clean (today's intersect returns
    a curve).
  - [x] **`ThreeDDisjointSidesetsAreRefused`** — HEX8, A = top, B = bottom; refuses "Curve 7",
    "share no facet edge". Red: debug throws the wrong message; release **crash**
    (`aEndNodes( 0 )` on an empty cell). Own filter.
- **What could break:** a shipped deck whose curve the clique bridged (A/B gate); a mesh whose A
  facets carry no nodes (never seen: the readers link them); 3-D mixed order (residue R-2).
- **Gate:** the tests above; `CurveFactory.*`; `make check-fast`; the A/B gate.

### R3 — D: thin-shell side curves take their midpoints from their own facet; delete the clique helpers (after: R2)

- [x] Landed 2026-09-22.

- **Diff:** a private `midnode_of_facet_edge( aFacet, aA, aB, aWork )` returns node 2 of the facet
  element's edge `{ aA, aB }`, or `nullptr` (no such edge, or no midpoint). The scratch cell is the
  caller's, hoisted out of the per-edge loop. `thin_shell_side_curves` gathers every protocurve's
  type and midpoints **before any `Curve` exists**; a missing midpoint deletes the protocurves,
  releases the scratch state and refuses, naming both node ids. The scratch release (vertex
  containers, the temporary edges, `update_node_indices()`) is one local lambda, shared with the
  success path and the closure refusal below.
- **Closed side curves (ruling on Z2):** refused. After the first `identify_subchains` and before
  `create_protocurves`, each component's node count and degree sum are compared; degree sum ≥
  2 × node count means the component holds a loop, and its closing edge would keep a stale index
  that `create_protocurves` subtracts one from and uses to index past `tNumChains`. The refusal
  releases the scratch state first and names a node: "The thin-shell side curve through node %lu
  is closed; closed side curves are not supported". It also refuses a loop with a branch.
- **Deleted:** `create_pairs` and `collect_midnodes( aPairs, aMidNodes )`, no caller left.
- **Tests:**
  - [x] **`ThinShellQuadraticSideCurveTakesItsFacetsMidpoints`** — HEX20 on a
    `Mesh( 3, 0, false )` (asserts that node 5 has no element), thin shell = top face, domain
    boundary = right face, so the side curve is 5–4–7–6; checks every midpoint (16, 19, 18), the
    node order and the arclength. Red: **crash** (`collect_midnodes` finds no element, writes no
    midpoint, `collect_nodes` dereferences the null). Own filter.
  - [x] **`ThinShellClosedSideCurveIsRefused`** — HEX8, thin shell = top, no boundary: refuses
    "is closed", "not supported"; then every mesh node has its position as index and the rim nodes
    have no vertex container. Red: undefined (out-of-range write in `create_protocurves`, or the
    old Z2 path). Own filter.
- **Consumer impact:** LINE3 midpoints only on the quadratic path; no shipped deck is quadratic.
  The closure refusal is new behavior for any deck with a closed thin-shell rim (today: undefined
  behavior; `CutFactory` passes `collect_boundary_sidesets`, so an open rim needs the tape ends on
  domain boundaries — Grok plan-jury finding 8, medium).
- **Gate:** the tests above; `CurveFactory.*`; `make check-fast`.

### R4 — Z2: orient along the walk, order segments by their first node (after: R2 for T6; independent otherwise)

- [x] Landed 2026-09-22.

- **Diff, `orient_segments( const bool aClosed, aEndNodes, aSegments )`:** `tLast = n − 1`; for
  each segment `tI`, `tJ` are its corners' walk indices; the closing segment
  (`aClosed` and `{ tI, tJ } = { 0, tLast }`) is flipped when `tI == 0`, every other one when
  `tI > tJ`. Only slots 0 and 1 swap. `aClosed` is needed because an open two-node chain also has a
  segment `{ 0, tLast }`.
- **Diff, `sort_segments( aSegments )`:** set each segment's index to its node 0's walk index,
  then `sort( aSegments, opSegmentIndex )`. The private `next( …, Segment * )` is deleted. No
  segment flag is touched any more (none was read, C4).
- **Call sites:** the three pipelines call `orient_segments( tCurve->is_closed(), … )` and
  `sort_segments( tCurve->segments() )`. The header's pipeline contract says what each sorter
  now does.
- **Tests:**
  - [x] **T4 `TwoDClosedLoopLinear` / T5 `TwoDClosedLoopQuadratic`** — QUAD4 / QUAD8, one
    `from_2d_sidesets` sideset on all four edges. Expects closed, `P0 P3 P2 P1 P0` (LINE3 with
    `M30 M23 M12 M01` interleaved), segments head to tail with the closing segment `P1 → P0`,
    arclength `0 … 4` in steps of 1 or 0.5. Red: clean (debug `collect_nodes` ASSERT; release
    chain checks).
  - [x] **T6 `ThreeDTopRimOfOneLinearHexIsOneClosedLoop` / `…QuadraticHex…`** — HEX8 / HEX20, A =
    {5}, B = {0, 1, 2, 3}: `P4 P7 P6 P5 P4`, midpoints 19, 18, 17, 16. Red: clean (the clique adds
    both top diagonals, degree 3, branch refusal) on HEX8; crash possible on HEX20 (the drift).
  - The thin-shell closed case is refused in R3, not ordered (ruling).
- **Consumer impact:** open curves unchanged (both orders are forced on a path). Closed curves
  from `from_2d_sidesets` / `intersect` now satisfy `collect_nodes`' contract instead of a debug
  abort or a wrong release list. **Not guarded downstream:** see residue R-3.
- **Gate:** T4–T6 green; `CurveFactory.*` on a **debug** tree too (closure ASSERT active).

### R5 — Z1: `from_2d_sidesets` refuses mixed or non-line sidesets (independent)

- [x] Landed 2026-09-22.

- **Diff:** before anything is allocated or re-indexed: an empty id vector is refused ("no sideset
  given"), zero facets are refused ("its sidesets have no facets, so the curve would be empty"),
  and in the facet concatenation loop every facet is checked to be LINE2 or LINE3 and of the type
  of sideset `aSideSets( 0 )`'s first facet, both with `to_string( ElementType )`:
  "all sidesets of one curve must have the same element type".
- **Tests** (`QuadEight` fixture: one QUAD8 with a LINE3 sideset 2 on edge 1; further sidesets
  are registered after `finalize()` with `Mesh::add_sideset`, facets without master, nodes
  inserted by hand with `tFacet->element()->insert_node` — `finalize()` would overwrite a linked
  facet from its master):
  - [x] **T8 `TwoDMixedOrderSidesetsAreRefused`** — a hand-noded LINE2 sideset 1 on edge 0:
    `{ 1, 2 }` and `{ 2, 1 }` both refuse "Curve 7", "same element type". Red: clean (`{1,2}`
    returns a LINE2 curve; `{2,1}` refuses with the old "has no midpoint" text).
  - [x] **T9 `TwoDNonLineSidesetsAreRefused`** — node-less LINE4 and QUAD4 facets: refuse with
    "line4" / "quad4" and "LINE2 or LINE3". Red: **crash** (null node). Own filter.
  - [x] **`TwoDEmptyCurvesAreRefused`** (O6) — an empty sideset and an empty id vector. Red:
    undefined or crash. Own filter.
- **Consumer impact:** the four 2-D decks pass LINE2 in every sideset (C3).
- **Gate:** the tests above; `CurveFactory.*`; `make check-fast`.

### R6 — Bundle gate (after: R1–R5; Christian runs it)

- [ ] `test/test_mesh --gtest_filter='CurveFactory.*'` on the release tree, and again on a debug
      tree (the asserts in `collect_nodes` and `Vertex::insert_node` are active only there).
- [ ] The red-on-pre-fix record, if wanted: the clean-red tests in one filter; each crash-red test
      (T3, `ThreeDDisjointSidesetsAreRefused`, `ThinShellQuadraticSideCurveTakesItsFacetsMidpoints`,
      `ThinShellClosedSideCurveIsRefused`, T9, `TwoDEmptyCurvesAreRefused`) under its own filter.
      The new test file calls only members the pre-fix header also has, so it builds against a
      tree with the old `src/mesh/cl_CurveFactory.*` restored (e.g. `scratchpad/dr188_ab/baseline/`).
- [ ] `make check-fast` (`mesh` is labeled fast).
- [ ] **The A/B gate (§4.1): mandatory before acceptance.**
- [ ] The code jury.

### 4.1 The A/B gate: curve creation on the 10 shipped 3-D decks, old vs new

**Question:** does the facet-built rule (R2) change, or refuse, any curve that a shipped 3-D deck
builds today? Compared per curve: segment count, ordered node ids, closed flag (and, as
information, the element type and the segment element ids). Curve creation only: no solve, no
`.bfm`, nothing written into the checkout.

**Decks (3-D, `@` syntax → `intersect`):** `corc_garber_effect`, `corc_periodic_bc`,
`corc_twolayer`, `pancake`, `racetrack_usermat`, `rutherford`, `tape_quench_usermat`,
`tapestack3d`, `tapestack3d_quench`, `tapestack_circuit`. The four 2-D decks (`gantry`,
`tapestack2d_layered`, `tapestack2d_racetrack`, `undulator2d`) can ride along with `ALL=1`; they
exercise R4/R5 only.

**Mechanism:** a standalone probe in the session scratchpad,
`/tmp/claude-1000/-home-christian-codes-belfem/a22f4264-8919-418a-a910-1b3ebff13dc2/scratchpad/dr188_ab/`:

| File | Role |
|---|---|
| `dr188_curve_probe.cpp` | loads the deck's `.msh` as `MaxwellFactory::read_mesh` does on a fresh build (read without connectivity, `MeshChecker`, `unfinalize`, `Connectivity::Compute`, `finalize`), then builds every curve of `topology { curves }` in key order with **one** factory, as `create_curves` does, and writes a canonical dump. Mode `old` or `new`. A refusal is caught (throw mode) and dumped as `REFUSED: <message>` |
| `cl_CurveFactoryOld.{hpp,cpp}` | the pre-DR-188 working-tree `cl_CurveFactory.*`, class renamed to `CurveFactoryOld` and nothing else (`baseline/` holds the unrenamed snapshot) — compiled into the probe, so one binary runs both and no stash or second library is needed |
| `decks.sh`, `mesh_decks.sh` | copy each deck to `work/<deck>` and mesh it there (gmsh `-3`/`-2` on the deck's `.geo`, or `python3 python/main.py` for `corc_twolayer`, `pancake`, `rutherford`); the checkout is not touched |
| `build_probe.sh` | compiles both sources with the `belfem` executable's `flags.make` and links with its `link.txt` (object and output swapped) against `build/lib/libbelfem.a` |
| `run_ab.sh` | per deck, runs `old` and `new` in separate processes, diffs the dumps, prints `SAME` / `DIFF` / `NO DUMP` |

**Run (Christian):**

```bash
cd /home/christian/codes/belfem/build && make belfem        # libbelfem.a with DR-188
P=/tmp/claude-1000/-home-christian-codes-belfem/a22f4264-8919-418a-a910-1b3ebff13dc2/scratchpad/dr188_ab
$P/mesh_decks.sh          # ten 3-D decks into $P/work ( ALL=1 adds the 2-D four )
$P/build_probe.sh         # BUILD=... to use another build tree
$P/run_ab.sh              # per deck: SAME / DIFF; dumps and diffs in $P/out
```

The `new` mode calls the new `sort_segments( Cell< Segment * > & )`, so a `libbelfem.a` built
before DR-188 fails to link rather than comparing old with old.

**Pass:** every 3-D deck prints `SAME`. A `DIFF` whose new side reads `REFUSED: … not one
connected chain` is the plan jury's bridge case: the clique joined two shared-edge components
through a co-element chord. That needs a ruling before acceptance (refuse, or widen the rule), not
a waiver. If the scratchpad is gone, the four sources are small enough to recreate from this
section and `baseline/` = the pre-DR-188 `src/mesh/cl_CurveFactory.*`.

### 4.0 Implementation Progress

- 2026-09-22: R1–R5 landed together in `src/mesh/cl_CurveFactory.{hpp,cpp}` and
  `tests/mesh/test_CurveFactory.cpp` (15 new tests). Gate reached: `-fsyntax-only` against
  `build/`'s `flags.make` (release defines, and again with `-UNDEBUG -DDEBUG`), plus
  `cl_CutFactory.cpp` and `cl_MaxwellFactory.cpp` against the new header. **Below the compile
  rung:** nothing linked, no test run. Red-on-pre-fix statuses above are predicted from the
  source, not observed.

## 5. Design Questions and Rulings

Christian's rulings of 2026-09-22 on the plan jury's reconciliation:

- **O1 — 3-D pair rule:** the facet-built rule (edges shared by an A facet and a B facet), with a
  B-side midpoint conflict check and a both-orders test. **Mandatory A/B gate** (§4.1) before
  acceptance.
- **O2 — run gate:** superseded by the mandatory A/B gate on all ten 3-D decks.
- **O3 — closed curves:** orientation fixed for `from_2d_sidesets` and the 3-id `intersect`, with
  full-pipeline closed-loop tests on LINE2 and LINE3. `thin_shell_side_curves` refuses a closed
  side curve at setup; full closed thin-shell support is residue R-4.
- **O4 — the C6 by-catch:** residue R-1, not fixed here.
- **O5 — `nonfree/` callers of the bitset overload:** none (Christian); deleted.
- **O6 — empty curves:** refused with named errors.
- D: facet-built midpoints on both 3-D callers; `create_pairs` and the old `collect_midnodes`
  deleted; required quadratic intersect (HEX20/QUAD8) and quadratic thin-shell
  (`Mesh( 3, 0, false )`) fixtures.
- Z1: per-facet type refusal with `to_string( ElementType )`. Z3: delete.

### 5.1 Residue (not fixed here)

- **R-1 — `CutFactory::close_terminal_loops` reverses a LINE3 duplicate wrongly.**
  `tDup->element()->insert_node( tNode, n-1-k )` (`src/homology/cl_CutFactory.cpp:~1598-1602`)
  puts the LINE3 midpoint into slot 0 and corner 0 into slot 2; a LINE3 flip swaps slots 0/1 and
  keeps slot 2. Reachable only from order-2 decks. `src/homology` is Gregory's module. **Open
  question for Gregory (Christian):** could cohomology run on the linear skeleton, with the order-2
  lift confined to one place?
- **R-2 — 3-D mixed order silently becomes LINE2.** A quadratic A with a linear B gives a LINE2
  curve (`Curve` ctor, `cl_Curve.cpp:27-31`), and `create_pairs_from_facets` then skips the
  midpoint checks. No shipped deck is second order.
- **R-3 — closed `curves` reach `close_terminal_loops` unguarded.** Terminal curves are the deck's
  `curves` (`MaxwellFactory`, `cl_MaxwellFactory.cpp:~3233-3260`), not thin-shell side curves, so
  R3's refusal does **not** cover this consumer. `close_terminal_loops` duplicates every non-empty
  terminal curve without reading `is_closed()` (`cl_CutFactory.cpp:1576-1700`). Before R4 a closed
  terminal curve aborted in debug at `collect_nodes` and went on, mis-closed, in release; after R4
  it goes on, correctly closed, in both. No shipped deck has a closed curve (C3, medium). A guard
  belongs in `MaxwellFactory`'s terminal association or in `close_terminal_loops`; both are outside
  this plan's files.
- **R-4 — closed thin-shell side curves are refused, not supported.** Support needs
  `identify_subchains` to give every edge its component (the closing edge keeps a stale index) and
  a downstream audit of `CutFactory`'s side-curve consumers.
- **R-5 — `sort_end_nodes` throw-mode index residue** (bundle §4.2 G2), unchanged.

## 6. Interface Changes (all inside `mesh::CurveFactory`)

| Member | Before | After |
|---|---|---|
| `intersect( id, id, DynamicBitset & )` | public | **deleted** (R1) |
| `create_pairs_from_facets(…)` | private | gains `const Cell< Facet * > * aSharedWith = nullptr` (R2); contract text rewritten |
| `create_pairs(…)` | private | **deleted** (R3) |
| `collect_midnodes( aPairs, aMidNodes )` | private | **deleted** (R3) |
| `midnode_of_facet_edge( aFacet, aA, aB, aWork )` | — | **new**, private (R3) |
| `orient_segments(…)` | public | gains leading `const bool aClosed` (R4) |
| `sort_segments(…)` | public | becomes `sort_segments( Cell< Segment * > & )` (R4) |
| `next( …, Segment * )` | private | **deleted** (R4) |

`orient_segments`, `sort_segments` and `sort_end_nodes` are public only for the tests. In this tree
nothing outside `cl_CurveFactory.cpp` calls the first two; a `nonfree` caller would fail to compile
(loud). `sort_end_nodes`' `aAdjacency` output is now read only by the tests.

## 7. Definition-of-Done Checklist

- [x] Every gap-table row mapped to a step, a ruling or a residue (G1–G13 → R1–R5, §5, §5.1).
- [x] Each claimed gap backed by a citation (§1, §3).
- [x] Ordered steps with dependencies (§4).
- [x] Open questions ruled (§5).
- [ ] The R6 bundle gate: `CurveFactory.*` release and debug; red-on-pre-fix record if wanted;
      `make check-fast`.
- [ ] The A/B gate (§4.1) passes, or its DIFFs are ruled on.
- [ ] Code jury (diff: `tmp/ai_exchange/review_curvefactory_dr188.diff`).
- [ ] DR-188 row updated: the stale overload name (C1), and the evidence and gates per item.
- [ ] Plan registered in `todo/README.md`.

## 8. Audit Trail

- Source of the findings: `tmp/ai_exchange/review_mesh_oneliners_dr161_169_176_179_180.md`.
  D comes from code round 2 (Grok, "Open Risks"); Z1, Z2 and Z3 from code round 3 (astra, and the
  Grok rerun). The bundle plan `todo/mesh_oneliners_dr161_169_176_179_180.md` §4.1–4.3 records G1
  (`create_pairs_from_facets`), G2 (`unique_ptr`) and fix pass 3 (`element()->insert_node`).
- This plan's own additions, each re-traced by Claude against the working tree and none executed:
  - D′: the clique defect on the 3-D path;
  - C2: the A-only rule refuted by the tet ear;
  - Z2 widened: `sort_segments`' first step is also ambiguous, so orienting the closing segment
    alone is not enough;
  - C6: the `close_terminal_loops` LINE3 reversal;
  - O6: empty curves.
- Plan jury (2026-09-22, `tmp/ai_exchange/review_curvefactory_dr188.md`): Codex gpt-6-astra and
  Grok 4.7. Accepted into this plan: the closed thin-shell failure in `identify_subchains` /
  `create_protocurves` (astra), the unguarded `close_terminal_loops` (both), the false open-curve
  identity lemma (both), the missing D tests (both), the B-side midpoint conflict (both), cleanup
  on the new refusal path (astra), and Grok's plan-text corrections (T2/T3 labels, the R4
  dependency, G6, the `sort_segments` start, the premature R1 contract sentence, the `set_master`
  citation, checkboxes as list items) and naming.
- Implementation finding (Claude, 2026-09-22): the thin-shell closure refusal does not cover
  `close_terminal_loops`, whose curves come from the deck's `curves` (residue R-3).

## Appendix A — Z2 trace and model check

**Trace (closed 4-loop, walk indices 0-1-2-3).** In `sort_end_nodes`, the start is the smallest
id and `tEnd = tStart` (`:1096-1108`). The walk assigns 0…3 along the loop and breaks when both
neighbours are flagged (`:1123-1135`), and the start's index is forced back to 0 (`:1140`).
`orient_segments` then gives `0→1, 1→2, 2→3` **and `0→3`**. In `sort_segments`, the first segment
is whichever segment sits in row 0, slot 0 of the adjacency (`:1219`), which depends on the order
of the segment cell. `next()` from `0→1` sees unflagged candidates at both ends, `0→3` at node 0
and `1→2` at node 1, and returns the smaller element id (`:1397-1400`). The traversal can therefore
start in either direction. In both directions the segment list is not head-to-tail, because `0→3`
points against the walk. `collect_nodes` appends `node( 1 )` of each segment
(`:1419-1422`, `:1462-1466`) and asserts `first->node(0) == last->node(1)` (`:1447`).

**Model check.** A Python transcription of `sort_end_nodes`, `next()` on nodes and segments,
`orient_segments`, `sort_segments` and `collect_nodes` (LINE2) was run over every node-id
permutation × segment-cell order for 3-, 4- and 5-loops (36, 576 and 14 400 cases). The script is
in the session scratchpad and was not kept in the tree.
- Today, every case fails both the closure test and the head-to-tail test. A triangle gives
  `1, 2, 2, 3`.
- The R4 rule (closing segment n−1→0, segment index = first node's index) passes every case, with
  `1, 4, 3, 2, 1` for the 4-loop, which is the order `SingleLoopIsAccepted` expects.
- On 2000 random open chains, with random segment-cell order and random segment orientation, the
  R4 rule and today's code produce identical output.

Evidence rung: a numeric probe of a **model** of the code. It sits below a regression on the real
code, and T4–T6 are that regression.

**LINE3 under a flip.** `orient_segments` swaps element nodes 0 and 1 only (`:1204-1209`), and
node 2 is the midside (LINE3 numbering `(0)-(2)-(1)`, `:1429-1432`). `compute_coordinates`
integrates `|a ξ + b|` with `a = x0 + x1 − 2 x2` and `b = (x1 − x0)/2` (`:1530-1547`). A swap of
x0 and x1 keeps `a` and negates `b`, and on symmetric Gauss points `Σ w |aξ − b| = Σ w |aξ + b|`.
The midpoint and the segment length are therefore invariant. Confidence **high**.
