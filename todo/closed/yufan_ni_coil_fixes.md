# Yufan's NI-Coil Fixes onto devel: Shell Alignment Guard and Closed Side Curves

**Date:** 2026-10-06
**Purpose:** Land Yufan's branch `yufan` (one commit, `f29206f`, on v0.9.2) on `devel` with the
two defects the code jury found closed and the closed-side-curve policy narrowed as Christian ruled.
Mechanism: her `align_thin_shell_facets` per-patch vote stays as written; a block-consistent
relation between a sideset's patches and a block-set refusal after it make the block-by-block cut
rule explicit; the DR-188 closed-rim refusal moves from `CurveFactory` to `CutFactory`, where the
terminal curves are known, and refuses only a closed rim that shares a node with a terminal curve
of its shell.
**Module:** `src/fem/maxwell` (`cl_MaxwellFactory`), `src/mesh` (`cl_CurveFactory`),
`src/homology` (`cl_CutFactory`, outside the closed core), `tests/mesh`, `tests/homology`, `tests/fem`
**AIs involved:** Claude Fable 5.1 (merge, pre-registration, plan, reconciliation); Codex
gpt-6-astra/high and Grok grok-4.7/high (four jury rounds: code diff, policy question, plan, R2b);
Claude Opus 5.5 (probe builds and runs, the cut diagnosis, the generator fix, the kernel diagnosis,
implementation)
**Status:** ✅ COMPLETE 2026-10-08 — landed on devel in one commit (Christian's rulings of
2026-10-06/08). Landed: Yufan's `align_thin_shell_facets` and `identify_subchains` fixes unchanged;
R2b block-consistent patch relation and R2 `check_thin_shell_block_sides` (one-sided facet and
block-on-both-sides refusals); R3 closed rims built and accepted when isolated
(`CutFactory::refuse_closed_side_curve_sharing_nodes`, before duplication, per-shell mesh ownership,
shrink-before-delete on the throw), degree-3 branch refusal in `thin_shell_side_curves`; R4 tests
(`test_mesh` 131/131, `test_homology` 32/32, `test_fem` 319/319; every refusal red without its
guard); R5 DR-194/195/196/197; R7 `examples/pancake_tape` with the `SolderAir` sideset. Verified:
`make check` 17/17 three times in Claude's private debug tree (MKL, Armadillo); the 20-turn
pancake positive control (unpatched devel red in `assign_edges`, patched green through
`CutFactory::run`). Juries: four plan-level rounds, a two-seat code round, Grok xhigh on fix round
1, Codex astra/high on fix round 2; no open P0/P1. **Owed:** the copper-lead end-to-end gate (R6),
blocked on the kernel φ-dof gap DR-196, which also keeps `pancake_tape` from running to a solution;
untested residues: the parity-conflict refusal, the two flag-clear throw paths, the later-shell
build-throw path (DR-195, DR-188 (b)). Residue: `thin_shell_facet_orientation.md` still cites
`todo/orientation_fix_block_ordering_bug.md` three times (pre-existing).

> **Scope guards.**
> - In scope: Yufan's two source changes and her records; the two P1s of the code jury; the
>   narrowed closed-rim policy with its tests; the DR rows the round produced; `pancake_tape` as a
>   shipped example (R7).
> - Out of scope: the `update_node_facet_tables` allocate/unflag inconsistency (DR row, not fixed
>   here); the kernel dof-source assertion behind the SEGV (own diagnosis, O5); 2-D, which never
>   runs the same-type walk; the cohomology core (`cl_Homology`, `cl_Cohomology`,
>   `cl_SimplicialComplex`, `cl_Chain`, `cl_Cochain`, `fn_Smith`), read only, Gregory's.
> - Compatibility promise (narrowed by the plan jury): a 3-D deck whose thin-shell sidesets are
>   each **one connected patch**, uniformly wound, whose side curves are open, and whose sidesets
>   each have disjoint master-block and slave-block sets behaves bit-identically. A uniformly wound
>   sheet with one block on both of its sides is refused by R2 from now on; a sideset in several
>   patches may have patches flipped by R2b. Evidence: the alignment pass flipped nothing on the
>   200-turn `pancake_tape` deck without an interface sideset ("no facet flipped", not output
>   equivalence); with the interface sideset it flipped 14 + 22 facets and removed a devel abort
>   (see §1, last row).
> - Not shipped: `examples/helix/Copper_RRR50.hdf5` (a 0.9.2 rho cache; devel's caches carry a
>   revision and would rebuild it on the first run).

---

## 1. Current Behaviour and How It Fails

Yufan's branch, read against devel `8565205` (the trial merge is conflict-free in source; only
`CHANGELOG.md` and `devlog/README.md` needed merging):

| Failure | Mechanism | Evidence |
|---|---|---|
| **Mixed master side on a 3-D shell sideset** (her trigger, NI pancake with copper between the turns) | the same-type walk in `fix_facet_masters` passes through faces that meet the tape at a T-junction, so one sideset is reached from two directions with opposite windings; stock 0.9.2 then corrupted the heap in `ThinShellFactory::update_node_facet_tables` | `cl_MaxwellFactory.cpp:1512-1535` (walk); her devlog `dl20261003_thin_shell_alignment.md`. The devlog's mechanism (allocate over `tNode`'s duplicates, fill over `tOrg`'s) is **unproven**: the function is called on `mMasterNodes`, which `reset_node_indices` asserts are non-duplicates (`cl_ThinShellFactory.cpp:701`; the nodes arrive via `collect_nodes` `:851` and the call at `:174`), and for a non-duplicate `original()` returns `this` (`cl_Node.hpp:317-321`), so the two loops alias today |
| **Her fix leaves a new hole: a split cross-type vote** | `align_thin_shell_facets` lets the cross-type facets decide a patch and flips the minority, so where a conductor sits on the top in one region and the bottom in another, the lower domain type becomes master on the flipped facets; only a `message` follows | `cl_MaxwellFactory.cpp:1645-1663`, `:1691-1699` |
| **…which the cut pipeline cannot take** | `relink_slave_elements_with_duplicate_nodes` unions the *slave blocks* of the whole sideset, flags their elements, and relinks every flagged element that touches a shell node; a block that is slave on one part and master on another is relinked on both sides | `cl_CutFactory.cpp:1437-1446`, `:1448-1479` |
| **Closed side curve (her second trigger, copper leads)** | the rim of a shell that reaches no domain boundary is one loop; the pre-DR-188 `identify_subchains` left the closing edge without a chain index and `create_protocurves` indexed out of range. Her hunk stamps every reached edge. On devel the hunk is inert: DR-188 refuses any loop before `create_protocurves` | `cl_CurveFactory.cpp:528-530` (hunk), `:112-141` (refusal) |
| **Records contradict devel** | her CHANGELOG line "no longer aborts setup" and devlog `dl20261003_closed_side_curve.md` describe 0.9.2; the comment at `cl_CurveFactory.cpp:112-115` ("cannot place the edge that closes it") is false once the hunk is in | `CHANGELOG.md` `[Unreleased] ### Fixed` |
| **Contract drift** | `fix_facet_masters` promises domain-type order of masters; alignment now permits exceptions on shell sidesets; `thin_shell_facet_orientation.md:46` is titled "Two-Stage Algorithm" and `:48` cites a stale range `cl_MaxwellFactory.cpp:1459-1536` | `cl_MaxwellFactory.hpp:228-231` |
| **Why DR-188 refused closed rims, and what the literature says** | the refusal was the crash mechanism plus a scope fence (residue R-4/O3), not physics. Alves et al. 2022a §II-B duplicate the shell "except nodes and edges on its lateral boundaries ∂Γs"; their §III-B disc has a wholly interior rim and needs no cut. Schnaubelt et al. 2023 model the NI contact layer "in-between conducting materials" fed by conductor leads. A rim that **carries a thin-shell terminal** is the one configuration both rules forbid at once (free edge = single-valued, terminal = duplicated), and in BELFEM it cancels the terminal chain: the rim's nodes are kept out of the duplication bitset, `close_terminal_loops` then builds the return segments on the same nodes, `link_curves_to_edges` / `assign_edges` key both on the unordered node pair so they land on one mesh edge with opposite `edge_direction`, and `Chain::addSimplexToChain` erases the zero sum | `todo/closed/curvefactory_dr188.md:249-254`, `:398-400`; `cl_CutFactory.cpp:1299-1306`, `:1601`; `cl_EdgeFactory.cpp:382-392`; `cl_Curve.cpp:166-176`; `cl_Segment.hpp:170-173`; `cl_Homology.cpp:381-384`; `cl_Chain.cpp:115-118`; Alves et al. 2022a, 10.1109/TASC.2022.3143076; Schnaubelt et al. 2023, 10.1109/TASC.2023.3258905 (the §-quotations are Claude's reading of the local extractions; no auditor could check them) |
| **Her alignment fix is positive on a shipped-generator deck (executed)** | 20-turn `pancake_tape` with a solder/air interface sideset: unpatched devel, debug, aborts in `Curve::assign_edges` from `CutFactory::run` ("Edge not found for segment 138168", `ASSERT tEdge != nullptr`) because sidesets 83 and 88 keep a mixed master side and the duplication misroutes a curve segment; patched devel, debug, logs "aligned 14 of 200 facets of sideset 83", "aligned 22 of 480 facets of sideset 88" and completes `CutFactory::run` (cut creation 0.435 s) | gdb logs `scratchpad/pancake_gen/run20_dbg{,_base}/gdb.log`; `cl_CutFactory.cpp:135-138` |

**Bottom line:** her per-patch alignment is correct, leaves uniform decks alone, and removes a real
devel abort on our own mesh; what is missing is the refusal of the one geometry it cannot make
consistent with the cut (a block on both sides of a sideset) and a well-defined relation between a
sideset's disconnected patches. Her curve hunk is correct and necessary, and the policy it was
written for is sound in the literature as long as the closed rim carries no terminal.

## 2. Architecture

**Patch relation (R2b).** A sideset's patches are disconnected pieces (the BFS walks same-sideset
neighbors only, `cl_MaxwellFactory.cpp:1605-1627`; `connect_facets_to_facets` links every pair of
facets sharing an edge, `cl_Mesh_ConnectivityCalculator.cpp:505-518`). The only well-defined
relation between two pieces is through their **blocks**: pieces that have the same two blocks on
their sides are aligned so the same block is master on both; pieces with distinct block pairs have
no combinatorial relation and are left as voted. Rejected (focused jury, `review_r2b_patch_propagation.md`,
Codex and Grok independently): an edge-parity walk through carrier surfaces. On an edge shared by
three facets the shell facet is edge-consistent with exactly one of the two carrier facets and
inverted with the other (`fn_check_facet_orientation.hpp:97-106`), so no edge-parity rule defines a
shell↔carrier relation, and over the full adjacency graph the result depends on facet-index order.
That is the ambiguity that produced the original bug. Also rejected: refusing anchorless
multi-patch sidesets (Christian, 2026-10-06).

**Alignment guard (R2).** After R2b, for every sideset of every protoshell, gather the set of
master block indices and the set of slave block indices over its facets (two `DynamicBitset`s over
`Block::index()`, which `Mesh::finalize` writes densely, `cl_Mesh.cpp:799-800`; the idiom of
`cl_CutFactory.cpp:1432-1440`; `mesh::Block` carries no flag slots); if they intersect,
`BELFEM_ERROR` naming the shell, the sideset and one offending block **id**, with the remedy in the
message: split the sideset into pieces that each have one block per side and give each piece its
own deck sign. After R2b this catches a split **inside** a patch — (i) a split cross-type vote on a
single block pair, where the minority flip puts that pair's higher block into the slave set, and
(ii) a same-type mesh that breaks the rule with no cross-type facet at all, which the vote never
sees. It does **not** catch a split vote across two distinct block pairs (masters `{H1, L2}`, slaves
`{L1, H2}`, Codex's counter-example): that configuration is accepted, and the existing split-vote
`message` stays as its diagnostic. Rejected: turning the split-vote message alone into an error
(misses (ii)); excluding cross-type facets from the flip set (leaves the mixed winding her deck
crashed on).

**Closed-rim policy (R3).** `CurveFactory::thin_shell_side_curves` keeps the loop detection
(degree sum ≥ 2·nodes) but stops refusing; it builds the closed `Curve` through the DR-188 closure
path (`sort_end_nodes` sets the flag at `:1110`, `orient_segments` runs the closer n−1 → 0 at
`:1222-1232`, `sort_segments` keys by walk index at `:1240-1247`, `collect_nodes` repeats the start
node). The refusal moves to `CutFactory::create_side_curves_for_thinshells_3d`, the first place that
holds a protoshell's side curves and its terminal curves together: a closed side curve that **shares
a node** with any terminal curve of its shell is refused, before `duplicate_nodes_on_face_sidesets`
(`cl_CutFactory.cpp:1146` vs `:1157`) mutates anything. Node sharing, not edge sharing, is the
criterion on purpose: a node cannot be both a free-edge node and a duplicated terminal node, so a
terminal that only touches the rim at a corner is refused too (the release note says so). A branched
loop still dies in `sort_end_nodes` (loud; its throw-mode scratch leak is DR-188 residue (b),
`todo/debt_register_closed.md:189`). Rejected: keeping the refusal in `CurveFactory` (it does not
know the terminals, so it can only refuse everything or nothing); lifting without the terminal check
(the silent cancellation above).

**What "closed" means downstream, qualified (G10).** The free-edge exclusion in `CutFactory`
is keyed on the flag the 3-D branch **overwrites** at `:1258-1284` ("shares no node with another
shell's side curve"), not on `CurveFactory`'s geometric flag. On a single shell every side curve,
open or closed, is excluded from duplication; a closed rim that shares a node with another shell's
rim is marked open and duplicated. This plan does not change that rule. The supported case is the
isolated closed rim; two terminal-free closed rims sharing a node is **unsupported and untested**,
recorded in R5's register row, not refused here.

## 3. Gap Table

| # | State | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| G1 | master/slave block sets of each shell sideset disjoint | the block-by-block relink | no (Python generator only) | (c) → R2 | `cl_CutFactory.cpp:1437-1479` |
| G2 | split cross-type vote | diagnosis of accepted two-pair splits; refusal of one-pair splits | `message` only | (c) → R2: one-pair case refused by G1, two-pair case keeps the message | `cl_MaxwellFactory.cpp:1648-1651`, `:1691-1699` |
| G3 | disconnected patches of one sideset | relative sign of the pieces | majority winding per patch, warning that does not say so | (c) → R2b | `:1595-1603`, `:1645-1646`, `:1685-1691` |
| G4 | `fix_facet_masters` contract | readers | false | (c) → R1 | `cl_MaxwellFactory.hpp:228-231` |
| G5 | orientation doc, three stages, stale range | readers | two stages, `:48` range stale | (c) → R1 | `thin_shell_facet_orientation.md:46`, `:48`, `:132` |
| G6 | closed rim: build the curve | copper-lead decks | refused | (c) → R3 | `cl_CurveFactory.cpp:112-141` |
| G7 | closed rim sharing a node with a thin-shell terminal | not cancelling the terminal chain | not detected | (c) → R3 | `cl_CutFactory.cpp:1299-1306`, `:1601`; `cl_Chain.cpp:115-118` |
| G8 | closed rim sharing an **edge** with a bulk terminal sideset | the bulk-terminal generator adds any curve with one flagged edge | not detected | (b) → O2, ruled: register row `[W]` | `cl_Homology.cpp:240`, `:267-285`; side curves are mesh curves, `cl_CutFactory.cpp:1985-1991` |
| G9 | edge coating / side connectors on a closed rim | `compute_binomial_vectors` anchor | hard error, message says "closed-loop tapes have no terminals" (wrong for a copper-lead tape) | (c) → R3 (message) | `cl_ThinShellFactory.cpp:2199-2205`, reached from `:491` only under edge coating |
| G10 | rim nodes single-valued | physics (Alves ∂Γs) | yes for an isolated rim; the predicate is the overwritten flag | (a) qualified, see §2 | `cl_CutFactory.cpp:1258-1284`, `:1299-1306` |
| G11 | `.bfm` reload skips `create_cuts`, keeps old windings | restart from an earlier processed mesh | by design | (a) documented in R1 (CHANGELOG) | `cl_MaxwellFactory.cpp:377`, `:445`, `:672-676` (verified by Grok) |
| G12 | `update_node_facet_tables` allocates/unflags over `tNode`, increments/fills over `tOrg` | the crash behind her trigger | alias under the non-duplicate assert; `Vertex::add_facet` ASSERT-only | (c) → R5 (DR row, mechanism unproven) | `cl_ThinShellFactory.cpp:2754`, `:2770`, `:2788`, `:2802`; `cl_Vertex.cpp:250-253` |
| G13 | records: CHANGELOG, comment `:112-115`, her devlogs' "no longer aborts" | truth of the release text | false on devel | (c) → R1 | — |
| G14 | side curves built but not yet owned by the mesh when R3 refuses | no leak on the throwing path, no double free on success | n/a today | (c) → R3 | `cl_Protoshell.hpp:55-57`; `cl_Mesh.cpp:218-220`; `cl_Curve.cpp:45-50`; late push loop `cl_CutFactory.cpp:1985-1991` |
| G15 | node flag slot 0 dirty when R3 runs | correctness of the terminal-node test | `sort_end_nodes` flags and never clears (`:1126-1147`); `unflag_all_nodes` only after the side-curve function returns (`cl_CutFactory.cpp:1916-1922`) | (c) → R3 | — |
| G16 | a shell facet without a master or a slave | R2 must not dereference null; the cut relinks dereference `slave()` unconditionally | `align_thin_shell_facets` refuses a slave-less facet only when it flips it (`:1657-1661`) | (c) → R2 | `cl_Facet.hpp:546-549`; `cl_CutFactory.cpp:1437-1440`, `:1518-1522` |
| G17 | `pancake_tape` ships without a conductor/air interface sideset, so its thin cut collapses | a runnable pancake example; the copper-lead gate | generator fix exists as a diff in a private copy | (c) → R7 | `cl_CutProcessor.cpp:395-398`, `:503-546`; `cl_Topology.cpp:217-245`; probe numbers in O3 |

### 3.1 Cross-cutting findings

- **One rule, two guards.** G1 and G7 are the same kind of statement: the cut opens a shell block
  by block, and it duplicates terminal nodes; any mesh that asks for both of two incompatible things
  at one entity (a block on both sides; a rim node that is both free and a terminal) must be refused
  at setup, with the sideset or curve named. Neither existed before; both are `BELFEM_ERROR`,
  setup-tier, rank 0, and both leave flags and scratch clean before throwing, because the test
  harness catches the exception and keeps the mesh (`test_ThinShellTerminalCurves.cpp:293-303`).
- **The per-patch vote is kept verbatim.** Nothing in three readings and two executed decks found a
  defect inside it. In-function edits: two message texts and the R2b patch relation, which runs
  after all patches are internally consistent and before the block-set check.
- **Testability drives the shape of R3.** The refusal is a static helper on `CutFactory` that takes
  a `Mesh *` and a `Protoshell *`, so a unit test can build a `Protoshell` by hand and call it
  without a `Topology` or a cohomology run.

## 4. Ordered Steps

- [x] **R1 — Records and contracts (no behavior).** Rewrite her CHANGELOG `[Unreleased] ### Fixed`
  entries: the alignment entry gains the patch relation and the block-set refusal and names the
  configuration that is now refused (one block on both sides of a sideset, including a uniformly
  wound one); the closed-curve entry becomes "a 3-D thin shell whose outer edge is a closed loop is
  accepted when no terminal curve of that shell touches it; a closed rim that shares a node with a
  terminal curve is refused", and notes that a processed `.bfm` from an earlier build keeps its
  windings. Replace the `cl_CurveFactory.cpp:112-115` comment with the loop-detection intent. Fix
  the `fix_facet_masters` Doxygen (`cl_MaxwellFactory.hpp:228-231`): domain-type order, then
  per-patch alignment, then the block-consistent patch relation, then the block-set refusal. Add the
  third stage to `thin_shell_facet_orientation.md`, replacing the stale range at `:48`, and run the
  Codex language sweep on that doc only. Leave her two devlogs as written (records); add a
  `> **Superseded on devel**` note at the top of `dl20261003_closed_side_curve.md` pointing at this
  plan's devlog.
- [x] **R2 — Block-set refusal** (after: R2b). New private
  `MaxwellFactory::check_thin_shell_block_sides()` called from `fix_facet_masters` after
  `align_thin_shell_facets`, 3-D only, rank 0. Two `DynamicBitset`s sized `mMesh->number_of_blocks()`
  constructed once, reset per sideset. Per facet: `BELFEM_ERROR` if `! has_master() || ! has_slave()`
  (the cut's relinks dereference `slave()` unconditionally, `cl_CutFactory.cpp:1437-1440`,
  `:1518-1522`; the implementer may cite an earlier always-active check instead); set the master and
  the slave block index. A block index set in both → `BELFEM_ERROR` with shell label, sideset id,
  `Block::id()` and the split-the-sideset remedy. **Keep** the split-vote message (`:1691-1699`); it
  is the only signal for the accepted two-pair case. Reword the patch message (`:1685-1691`) to what
  the code does: each patch takes the **majority** winding of its own facets (the seed's only on a
  tie); R2b then aligned the patches that share a block pair to one master block, or left them
  unrelated because their pairs differ (say which); a deck sign flips the whole sideset.
- [x] **R2b — Block-consistent alignment of a sideset's patches** (after: nothing; Christian's
  third O1 ruling). Inside `align_thin_shell_facets`, once every patch of a sideset is internally
  consistent and before the per-sideset unflag at `:1702-1706`: record per patch its block pair
  {master block index, slave block index}, read from the patch's first facet after the vote's flips
  (every facet of a consistent patch has the same pair; `BELFEM_ASSERT` it over the patch in debug).
  Group the patches by unordered pair. Within a group, the reference patch is the first one whose
  vote had cross-type voters (anchored), else the one holding the lowest facet index; every other
  patch in the group whose master block differs from the reference's is flipped whole (`Facet::flip`
  per facet; a slave-less facet is refused as at `:1657-1661`). Patches whose unordered pair is
  unique in the sideset are left as voted. Storage: one `Cell< Facet * >` per patch already exists
  as `tPatch`; keep the patches of the current sideset in a `Cell< Cell< Facet * > >` sized to the
  patch count, plus a `Cell` of block-index pairs; no new scratch flag. The multi-patch message
  reports the number of patches, how many were aligned to a reference and how many were left
  unrelated. Rejected: the edge-parity carrier walk (§2) and refusing anchorless multi-patch
  sidesets.
- [x] **R3 — Closed rim** (after: nothing; shares the devlog with R2).
  `thin_shell_side_curves`: keep the per-component E ≥ N detection as a comment-level invariant,
  drop the `BELFEM_ERROR`, keep `tReleaseScratch` on the remaining refusal paths.
  New `static void CutFactory::refuse_closed_side_curve_with_terminal( Mesh * aMesh, Protoshell * aShell )`,
  called in `create_side_curves_for_thinshells_3d` right after `thin_shell_side_curves` returns for
  that protoshell: `aMesh->unflag_all_nodes()` first (G15); flag the nodes of `terminal_curves()`;
  for each side curve with `is_closed()` (the geometric flag, still untouched at this point) and a
  flagged node: delete this protoshell's side curves (each `Curve` deletes its segments), clear
  `side_curves()`, `unflag_all_nodes()`, then `BELFEM_ERROR` naming the shell, the side curve's
  first node id and the terminal curve id, with the remedy "end the tape on a domain boundary, or
  keep its terminal curves off a closed rim". On success `unflag_all_nodes()` and push this
  protoshell's curves into `mMesh->curves()` at once, **deleting the late push loop at
  `cl_CutFactory.cpp:1985-1991`** (G14: otherwise every accepted curve is pushed twice and
  `Mesh::~Mesh` double-frees it; the labels loop `:1977-1983` stays and runs before the pushes). Fix
  the `compute_binomial_vectors` message (`:2203-2204`): the anchor search needs the rim's first node
  on a terminal curve, and closed rims are not supported under edge coating.
- [x] **R4 — Tests** (after: R2, R2b, R3).
  - `tests/mesh/test_CurveFactory.cpp`: `ThinShellClosedSideCurveIsRefused` (`:1205-1237`)
    becomes `ThinShellClosedSideCurveIsBuilt`: HEX8, shell = top face, no boundary → one curve,
    `is_closed()`, four segments head-to-tail, the closer oriented n−1 → 0, five nodes with the
    start repeated, monotone arclength, every mesh node back at its position index, rim nodes
    without vertex containers. Red on devel.
  - `tests/homology/test_CutFactory_closed_rim.cpp` (new, helper-level): a hand-built `Protoshell`
    whose `side_curves()` holds the closed HEX8 rim from `CurveFactory` and whose
    `terminal_curves()` holds (i) a curve on that rim → refused, message contains the shell label
    and "terminal"; (ii) a curve sharing exactly one corner node → refused (the corner-touch case);
    (iii) a curve disjoint from the rim → passes, all node flags clear afterwards, the side curve
    still intact. Red on devel: the helper does not exist.
  - `tests/fem/test_ThinShellAlignment.cpp` (new fixture, modelled on the terminal-curve harness
    but with a `homology` section, two volume blocks and the materials a thin shell needs, edge
    coating off): (a) **positive mixed-winding regression** — a hand-written Gmsh 2.2 mesh whose
    shell sideset the same-type walk reaches from two sides (a T-junction face), so
    `align_thin_shell_facets` flips facets; after `create_cuts` every shell facet of the sideset has
    the same master block (the oracle for a conductor|conductor shell, where no domain type can be
    asserted), the block sets are disjoint, and `create_thinshells` completes; (b) the same-type
    intersection case: one block on both sides of a sideset, uniformly wound → "block … on both
    sides" refusal; (c) the uniform control from (a) with a non-split sideset → passes with no
    alignment message; (d) for R2b: a sideset of two disconnected conductor|conductor pieces with
    the same two blocks on their sides, wound anti-parallel in the mesh → after `create_cuts` both
    pieces have the same master block and the log names the alignment; (e) two pieces whose block
    pairs differ → both keep their voted winding and the message names them as unrelated. The
    fixture's cost (a cohomology run on a ~100-element mesh) is the risk; if it cannot be made to
    pass through `create_cuts` on a hand mesh in one session, (a)/(d)/(e) fall back to a
    `MaxwellFactory`-level unit of `fix_facet_masters` on a mesh with hand-set masters, and that
    fallback is recorded in §4.0, not silently taken. **The 20-turn `pancake_tape` run with the
    interface sideset is the executed positive control for (a)** until the fixture exists
    (§1, last row).
- [x] **R5 — Register** (after: nothing). File **DR-194 `[CODE] [W]`**: `update_node_facet_tables`
  allocates (`:2770`) and unflags (`:2802`) over `tNode`'s duplicates, increments and fills
  (`:2754`, `:2788`) over `tOrg`'s; benign while `reset_node_indices` asserts master nodes are not
  duplicates (`:701`), a heap overflow the moment that alias breaks, with `Vertex::add_facet`
  ASSERT-guarded; remedy: walk `tOrg` in all four loops and a `BELFEM_ERROR` capacity check in the
  setup function, not in the container primitive. Yufan's devlog names this as the 0.9.2 crash
  mechanism; that attribution is unproven (medium, no debugger trace). File **DR-195 `[RULING] [W]`**
  (DR-188 is struck, `todo/debt_register_closed.md:189`): closed thin-shell rims — the narrowed
  policy, the G8 bulk-terminal edge-contact residue (O2), the shared-rim case of §2, and the
  edge-coating refusal. File **DR-196 `[CODE] [P]`** for the kernel dof-source assertion (O5) once
  the diagnosis names the mechanism. Run `scripts/check_doc_claims.py` after editing the register.
- [ ] **R6 — Gates.** `make check` in Claude's private debug tree (SCLS debug, Armadillo), count
  recorded. The 200-turn pancake probe on the patched release build without the interface sideset
  must still print no alignment line and abort at the same `CutProcessor` assertion as unpatched
  devel. The 20-turn run with the interface sideset must still log the two "aligned" lines and pass
  `CutFactory::run`. The copper-lead end-to-end gate, **owed, blocked on O5**: a deck whose tape ends
  inside a conductor, through `create_thinshells` and a solve; the current oracle is
  **terminal-in = terminal-out** (continuity through the single current path) and, separately,
  the azimuthal current I_C2 reported by the solver equals the integral of the tape plus solder
  current through one radial half-plane within 1 % — two numbers, because radial bypass through the
  contact layer makes source and azimuthal current differ (Schnaubelt et al. 2023 §IV;
  `src/fem/maxwell/doc/buffer_cut_topology.md:18-21`).
- [x] **R7 — `pancake_tape` becomes a shipped example** (after: R1–R5 landed; Christian's ruling).
  Move `cmake-build-debug/pancake_tape/{python,input.conf,README.md}` into `examples/pancake_tape/`
  with the interface-sideset diff applied (`scratchpad/pancake_gen/interface_sidesets.diff`: a
  `SolderAir` sideset over every solder/air contact triangle, appended last so existing ids hold,
  undeclared in the deck because `Topology` types it `InterfaceCondAir`, `cl_Topology.cpp:221-237`,
  `:342-348`, `:501`; a `validate` check that each such triangle sits between one solder and one air
  tet and that the contact is fully covered). Add `Allrun`/`Allclean` wrappers like the other decks,
  regenerate `pancake_topology.conf`, keep the 200-turn default and document `main.py 20`. Register
  in `examples/README.md` (the deck count there changes); `todo/pancake_inner_leadout.md` points at
  the new location. The example is **not** runnable to a solve until O5 is fixed; its README says so
  in its closing status paragraph, and the CHANGELOG `### Added` entry says the same.

### 4.0 Implementation Progress

Started 2026-10-06 after approval: R1–R5 to an Opus agent in the scratchpad devel worktree, R7 as a
second run.

**R1–R5 landed 2026-10-06 (Opus agent, unstaged on top of Yufan's staged patch; not committed).**
Debug MKL tree with tests, zero warnings; `test_mesh` 130/130, `test_homology` 29/29 (25 + 4 new),
`test_fem` 316/316 (310 + 6 new). The code jury and `make check` are still owed.

- **R2b** (`align_thin_shell_facets`): each patch is kept with its anchored flag, smallest facet id
  and block pair; patches are grouped by unordered pair after the conflict check; the reference is
  the first anchored patch, else the smallest facet id; others flipped whole. **Deviation:** the
  plan's `BELFEM_ASSERT` that every facet of a patch has the same pair was not written — one
  connected sheet can lie between block 1 and blocks 2 and 3, which is valid input. Such a patch
  gets the pair `gNoIndex` and is counted as unrelated (left as voted); R2 still checks its sideset.
  The multi-patch message now reports patch count, aligned (and flipped) count, unrelated count.
- **R2** `check_thin_shell_block_sides()`: two `DynamicBitset`s over `Block::index()`, a
  `BELFEM_ERROR` for a facet without master or slave ("facet %lu of sideset %lu needs an element on
  both sides") and for a block in both sets ("block %lu lies on both sides of sideset %lu; split the
  sideset …"). The split-vote message is unchanged.
- **R3** the `CurveFactory` loop refusal and its now-unused `identify_subchains` call are gone (a
  comment states the loop intent); `CutFactory::refuse_closed_side_curve_with_terminal` is public
  static, called per protoshell, and the protoshell's curves are pushed into `mesh->curves()` right
  after it; the late push loop is deleted, the label loop is unchanged (labels are assigned after
  the push; the pointers are the same objects). `compute_binomial_vectors` message reworded.
- **R4** all listed cases plus one: `tests/fem/test_ThinShellAlignment.cpp` runs the full
  `create_magnetic_kernel()` (no fallback needed): a Kuhn-tetrahedra box 5–6 × 5 × 2 cells, every
  conductor|air face and every outer face in a sideset (without them the air cut collapses, the O3
  mechanism; with a single z layer the gauge autopin finds no free air node), a minimal
  `solver { linear { library : mumps ; } }` so the kernel completes, and the shell facets read from
  the thin-shell factory's `GeometryOnly` sideset. Extra case: a shell on the domain boundary,
  refused by R2's one-sided-facet check. Red-without-guard checked by temporary local edits, each
  restored by `cp` from a backup and verified with `cmp`: R2 block check off → (b) red; R2b flip
  off → (d) red (R2 refuses it); the vote flip off → (a) red (R2 refuses it); R2 one-sided check off
  → the boundary case SEGVs; R3 guard off → helper cases (i) and (ii) red; the staged
  `cl_CurveFactory.cpp` → `ThinShellClosedSideCurveIsBuilt` red.
- **R5** DR-194 `[CODE] [W]` and DR-195 `[RULING] [W]` filed, header `[W]` count 3 → 5;
  `check_doc_claims.py` 38/38. DR-196 not filed (O5 diagnosis pending).
- **Found on the way, not fixed:** `CutFactory::~CutFactory` deletes the facets of
  `mTemporaryThinShellSidesets` and then the sidesets, whose own destructor deletes the same
  facets again: any exception between `duplicate_and_relink_facets` and
  `restore_thin_shell_sidesets` (for example "Thin cut 0 is empty") turns into a double free on
  unwinding (seen as a SEGV in the first fixture draft). A register row candidate.
- **R1 open item:** the Codex language sweep of `thin_shell_facet_orientation.md` (parent session).
- **R7 landed 2026-10-06 (Opus agent, unstaged):** `examples/pancake_tape/` (python with the interface-sideset diff, `diff -r` identical to `scratchpad/pancake_gen/python`; `input.conf` 200-turn ids, `SolderAir` = 813 noted; `Allrun`/`Allclean`; README with the `SolderAir` section and a closing status paragraph naming the dof defect in words, no DR id), `examples/README.md` 20 → 21 decks, CHANGELOG `### Added`, `todo/pancake_inner_leadout.md` pointer; `main.py 20` validates (SolderAir 6276 facets), `check_doc_claims.py` rerun. Codex language sweep of the README owed.

**Fix round 1, 2026-10-08 (Opus agent, unstaged; jury record `review_yufan_ni_coil_code.md`).**
Christian's rulings: fix everything; a closed rim is accepted only when isolated. Debug MKL tree,
zero warnings; `test_mesh` 131/131, `test_homology` 32/32, `test_fem` 318/318;
`check_doc_claims.py` 38/38.

- **Closed-rim scope (Grok 7).** `CutFactory::refuse_closed_side_curve_sharing_nodes( Mesh *,
  Cell< Protoshell * > & )` replaces the per-shell helper: a geometrically closed side curve is
  refused when it shares a node with a terminal curve of any shell (pass 0) or with a side curve
  of another shell (pass 1). `create_side_curves_for_thinshells_3d` now builds every shell's side
  curves, runs the check once, then pushes all of them; a refusal deletes every shell's side
  curves. Residue: a throw inside `thin_shell_side_curves` for a later shell leaves the earlier
  shells' unpushed curves unowned (throw mode only; release aborts). Three helper tests added
  (other shell's terminal, other shell's side curve, two disjoint rims pass and are pushed).
- **Branch check (Codex 2).** `thin_shell_side_curves` refuses a node with more than two rim
  edges right after `connect_edges_to_nodes`, after `tReleaseScratch()`;
  `ThinShellBranchedSideCurveIsRefused` (two HEX8 top faces touching at one node).
- **Messages (Grok 1–4).** Block-set remedy now says to split the volume block, not the sideset;
  the split-vote line states only the fact; the multi-patch line counts unique-pair and mixed-pair
  patches separately; the `compute_binomial_vectors` message names both causes of a missing anchor.
- **Flag hygiene (Grok 6).** `align_thin_shell_facets` clears slots 1–3 of the current sideset
  before the missing-slave and the not-orientable refusals; neither path has a test (no
  Kuhn-box mesh reaches them). `check_thin_shell_block_sides` uses no flags.
- **Group once (Codex 3).** One pass assigns each patch a group by its unordered pair and updates
  the group's reference incrementally; same reference rule.
- **Tests (Grok 8).** `SplitVoteIsRefusedByTheBlockSets` (f) and `MixedPairPatchIsAccepted` (g;
  one patch prints no patch line). The parity conflict stays untested.
- Each new or changed guard shown red with the guard removed (local edit, rebuild, restore by
  `cp`, `cmp`). Doc ranges fixed (`1473-1560`, `1565-1918`) and the Stage 2 text matched to the
  new messages; CHANGELOG and DR-195 ( 2 ) updated; README sweep applied except two edits.

**Fix round 2, 2026-10-08 (Opus agent, unstaged; Grok xhigh leg G2–G6).** Debug MKL tree;
`test_mesh` 131/131, `test_homology` 32/32, `test_fem` 319/319; `check_doc_claims.py` 38/38.

- **Ownership (G2).** `create_side_curves_for_thinshells_3d` records `mMesh->curves().size()`,
  pushes each shell's side curves onto the mesh right after its `thin_shell_side_curves` call,
  then runs the sharing check. The helper gained `const index_t aNumMeshCurves`; its refusal
  first shrinks `aMesh->curves()` back to that count (`Cell::set_size`, i.e. `std::vector::resize`,
  which keeps the prefix), then deletes the side curves, empties each `side_curves()`, clears the
  flags and throws. A throw while building a later shell leaves the earlier curves owned by the
  mesh. Correction to G2: devel also built every shell before pushing any, so the leak was the
  pre-existing residue round 1 named, not a regression of round 1; it is closed now either way.
  The closed-rim tests push the curves first, behind one curve the mesh already owns, and check
  that a refusal leaves exactly that curve and that survivors sit in the mesh once. Red without
  the shrink: count 2 vs 1, then SEGV from the double free in `~Mesh`. Not driven end to end:
  the function is private and needs a `Topology`.
- **Grouping (G3).** `Map< luint, index_t >` from the ordered pair ( min, max ), keyed
  `min * number_of_blocks + max`, to the group; same reference rule; O(P).
- **README (G4), DR-195 ( 4 ) (G5).** The generator cuts the solder; the guard is 3-D only.
- **Mixed-pair count (G6).** `MixedPairPieceIsCountedSeparately`: two pieces, one mixed-pair,
  one unique-pair, line reports 1 and 1. Red with the mixed counter not incremented. A swap of
  the two counters is caught by `PatchesWithDistinctBlockPairsAreUnrelated` (2 vs 0), not by the
  new 1/1 case.

## 5. Open Design Questions

- **O1 — Disconnected patches of one sideset.** Three rulings on 2026-10-06. First: propagate
  across the T-junction (after the plan jury showed "warn" was unreachable for one block pair).
  Second: implement it here → R2b v1, an edge-parity carrier walk. **R2b v1 refuted by its focused
  jury** (Codex and Grok independently): on a 3-facet edge the shell facet is edge-consistent with
  one carrier and inverted with the other, so no edge-parity rule relates shell to carrier; the
  full-graph walk is facet-index-order dependent; the conflict refusal would false-refuse every
  T-edge with two carriers; propagation did not run to completion; all anchored patches were seeded
  at one parity. **RESOLVED 2026-10-06 (third ruling) → block-consistent alignment**, R2b as written
  above: the only well-defined relation between pieces is through their blocks. Audited by the code
  jury together with the code (Christian, approval ruling).
- **O2 — Closed rim sharing an edge with a bulk terminal sideset (G8).** RESOLVED 2026-10-06 →
  register row only (DR-195), no check; the bulk branch flags *edges* (`cl_Homology.cpp:240`), so
  the residue text says "shares a flagged edge". Changing the bulk branch would touch the closed core.
- **O3 — Pancake thin-cut collapse on devel.** RESOLVED 2026-10-06 (probe, verified) → **not a
  regression, and not the cohomology generator.** `pancake.exo` of 09-24 is `msh2exo` output, the
  deck never ran, the guard predates it. The dangling-edge probe on the 20-turn reproducer (first
  pruning pass, cut 0): 3606 emitted faces, 5634 distinct edges, 349 exempt boundary edges,
  **207 dangling edges, all 207 on the air/conductor surface, 0 interior**; 192 of the 207 emitted
  faces at those edges lie in air, 15 on the interface, 0 in the conductor. The sheet lives in the
  air and ends on the solder surface, which carries no sideset in `pancake_tape`, so its edges are
  never exempt (`cl_CutProcessor.cpp:395-398` exempts `mPhiBoundariesAndPeriodic` only, and
  `Topology` types existing sidesets, `cl_Topology.cpp:217-245`, it creates none). The generator fix
  (R7) takes the 20-turn case through cut creation (0.14 s release / 0.435 s debug, no pruning
  line; the pristine generator on identical nodes still collapses). The cut supports changed with
  identical nodes (1867/7822 → 1329/7311), so the sideset also changes the cohomology input — in
  Gregory's note. Gregory's three questions are in `tmp/gregory_note_cut_questions.md`.
- **O4 — Where the closed-rim refusal lives.** RESOLVED 2026-10-06 → `CutFactory`, before
  duplication. Gregory is told; `cl_CutFactory.cpp` is outside the six closed units.
- **O5 — The kernel dof-source assertion behind the SEGV.** With the interface sideset, the
  20-turn run passes the cuts and then dies: release SEGVs at binary-dependent sites; patched debug
  gives `DofData::create_dofwise_t_matrices_master`, `cl_FEM_DofMgr_DofData.cpp:3112`,
  `ASSERT tNode->source(j)->number_of_dofs() == 1`, "couldn't figure out how to assign dependencies
  of dof 1275943 (edge 561783)", rethrown from `Kernel::create_field` (`cl_FEM_Kernel.cpp:620`).
  **DIAGNOSED 2026-10-06 (Opus, read-only, gdb; citations spot-checked by Fable):** not the thin
  shell, not Yufan's patch, not the interface sideset as such. `create_hanging_edges_and_facets`
  (`cl_MaxwellFactory.cpp:~1755-1830`) hangs the conductor-side edges of every `InterfaceCondAir`
  facet on the air element's facet nodes, which after the cut relink can be cut duplicates;
  `CutSet::create_duplicates` (`cl_CutSet.cpp:70-80`) makes each duplicate depend on its original
  plus the λ node of every generator in its pattern; the edge branch (`:3094-3116`) assumes every
  such source carries exactly one dof (φ) and takes `dof(0)`. φ goes only to air, abstract and
  orphan nodes (`count_node_dofs` `:870-990`; `collect_orphan_nodes` `cl_CutFactory.cpp:2312-2342`
  adopts block-less nodes only). Failing source: node 6276 (r = 4.12 mm, inner lead bar, block 42,
  interior node of `SolderAir`, on no tape sideset): both cut sheets co-terminate on the interface
  edge 6571→6276 and their plus sides take all 8 of its air tets into duplicates (24175: λ₀,
  25930: λ₁, 26164: λ₀+λ₁); the original keeps 12 conductor tets, gets no φ, is no orphan →
  `ndofs = 0`. Release reads a null `mDofs` in the same function. Silent second effect: the node
  branch (`:3031-3086`) un-hangs a node dof whose source lacks the dof, dropping λ jumps there.
  Pre-existing devel gap, triggered by two or more cut sheets co-terminating on an interface edge;
  a single sheet always leaves the original a minus side; not R13's ±1-weight refusal. Alignment
  flips as cause: unlikely (node on neither 83 nor 88), untestable on unpatched devel (earlier
  abort). 0.9.3 blocker: probably not (no shipped deck known to have co-terminating sheets; none
  run; medium). Fix direction, own plan + jury, not a patch: an interface original with cut
  duplicates and no air element still needs a φ dof (adopt it like an orphan). Gregory asked whether
  two representatives should share a boundary edge on ∂Ωc at all (note, question 4). → **DR-196**.
- **O6 — By-catch of R4 (Opus, reading + a gdb trace on an early fixture draft):**
  `CutFactory::~CutFactory` (`cl_CutFactory.cpp:63-72`) deletes each facet of
  `mTemporaryThinShellSidesets` and then deletes the sideset, whose destructor deletes the same
  facets again; any exception between `duplicate_and_relink_facets` and
  `restore_thin_shell_sidesets` becomes a double free while unwinding (debug/throw mode). Not fixed
  here; → **DR-197 `[CODE] [W]`**, filed with R5's rows.

## 6. Interface

No input key changes. Contract strings for R4: "block %lu … on both sides of sideset %lu"; "side
curve … of thin shell %s shares node %lu with terminal curve %lu"; the reworded multi-patch
message (patch count, aligned count, unrelated count). `pancake_tape` adds one undeclared sideset
to its mesh; the deck's `topology` section is unchanged.

## 7. Definition-of-Done

- [ ] Every gap row maps to a step or an O.
- [ ] R2 and R3 refusals are `BELFEM_ERROR`, rank 0, setup tier, flags and scratch clean before
  the throw, each named in a test that is red without it.
- [ ] R2b: fixture (d) shows two same-pair anti-parallel pieces aligned to one master block, (e)
  shows two distinct-pair pieces left alone with the message naming it; no new flag slots.
- [ ] `ThinShellClosedSideCurveIsBuilt`, the three helper-level cases, and the alignment fixture's
  (a)–(e) green; (a) proves the flipping branch repairs a mesh, or its fallback is recorded.
- [x] `make check` 17/17 in the private debug tree (MKL flavor, Armadillo, debug): 2026-10-06 on R1–R5 + Yufan's patch (39 s), 2026-10-08 after fix round 1 (43.5 s) and after fix round 2 (45.3 s).
- [ ] Pancake probes: 200-turn without interface → no alignment line, same abort as devel; 20-turn
  with interface → two "aligned" lines, `CutFactory::run` passes.
- [ ] Records: CHANGELOG, comment, header, doc, DR-194, DR-195, DR-196, devlog, `todo/README.md`,
  `examples/README.md`, R7's README with its closing status paragraph.
- [ ] Gregory informed of the cut-side change (isolated closed rims now reach
  `duplicate_nodes_on_face_sidesets`) and of O3's questions.
- [ ] The Thickcuts session receives the devel head SHA once the commit is on local devel (it holds
  the devel→thickcuts merge until then; devel itself was pushed at `8565205` on 2026-10-06).

## 8. Audit Trail

- `tmp/ai_exchange/review_yufan_ni_coil_fixes.md` — code jury on her diff: both auditors found the
  split-vote → slave-block relink hazard independently.
- `tmp/ai_exchange/review_closed_side_curve_policy.md` — policy round: keep the hunk, do not lift
  wholesale; both traced the terminal-cancellation chain; Claude added the Alves/Schnaubelt reading.
- `tmp/ai_exchange/review_yufan_ni_coil_plan.md` — plan jury: Codex — split vote ≠ block
  intersection, G10's predicate is the overwritten flag, the refusal-before-ownership leak, no
  positive alignment regression, the current gate undefined, the compatibility promise too wide.
  Grok — O1(a) contradicted R2, the message reword overstated the seed rule, the tNode/tOrg alias,
  the dirty node flags on R3's entry, the missing `has_slave` guard, the harness never reaching the
  refusals, node- vs edge-contact for G7/G8, G11 verified.
- `tmp/ai_exchange/review_r2b_patch_propagation.md` — focused R2b jury: the carrier walk refuted
  (3-facet edge argument, incomplete propagation, anchor seeding), R3's double push, the missing
  slave guard's downstream dereference, citation fixes, the test (a) oracle. Every citation was
  re-read against the staged tree before it entered this version.
