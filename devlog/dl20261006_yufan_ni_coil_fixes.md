# Yufan's NI-coil fixes: merge onto devel, three jury rounds, closed-rim policy, pancake cut collapse

**Date:** 2026-10-06
**Purpose:** record how Yufan's branch `yufan` (`f29206f`, on v0.9.2) was read against devel, what
three jury rounds found, the rulings Christian made, and what the two pancake probes showed
**Module:** `fem/maxwell`, `mesh` (CurveFactory), `homology` (CutFactory, CutProcessor; the
cohomology core was read, not touched)
**AIs:** Claude Fable 5.1 (merge, pre-registration, plans, reconciliation); Codex gpt-6-astra/high
and Grok grok-4.7/high (jury seats, all rounds); Claude Opus 5.5 (probe builds, diagnosis and the
dangling-edge probe, the generator fix)
**Status:** CLOSED 2026-10-08 — committed on local devel as the commit that carries this devlog (see the
close-out line at the end); plan `todo/closed/yufan_ni_coil_fixes.md`; owed: the copper-lead end-to-end gate, blocked
on DR-196.

## What Yufan changed (`f29206f`, "revisions for ni coils", 8 files)

- `MaxwellFactory::align_thin_shell_facets()`, called at the end of the 3-D same-type walk in
  `fix_facet_masters`: per thin-shell sideset, a BFS over same-sideset neighbors records each facet's
  winding relative to its patch seed; cross-type facets vote when the patch has any, else all
  facets; the minority is flipped; a non-orientable patch is a `BELFEM_ERROR`. Trigger: an NI
  pancake deck with copper between the turns, where the walk reached one sideset from two sides and
  stock 0.9.2 corrupted the heap in `ThinShellFactory::update_node_facet_tables`.
- `CurveFactory::identify_subchains`: every edge reached from a chain gets the chain's index, so
  the edge that closes a loop is no longer left at `gNoIndex`. Trigger: a tape fed by copper leads,
  whose side curve is one closed loop; 0.9.2 crashed with an Armadillo bounds error.
- CHANGELOG `[Unreleased]`, two devlogs dated 2026-10-03, and `examples/helix/Copper_RRR50.hdf5`
  (a 1 MB rho cache, committed by accident; devel's caches carry a revision since `77b4567` and
  would rebuild it).

## Reading it against devel (133 commits ahead)

Trial merge (`git merge-tree`): source merges cleanly, `CHANGELOG.md` and `devlog/README.md`
conflict trivially. The clean patch (her commit cherry-picked onto `8565205` without the `.hdf5`,
conflicts resolved) is staged in a scratchpad worktree; `tmp/ai_exchange/yufan_ni_coil_fixes.diff`
and `_src.diff` hold it.

Two collisions with devel decisions:
- **DR-188** (`0243908`, 2026-09-23) refuses any closed thin-shell side curve before
  `create_protocurves` (`cl_CurveFactory.cpp:112-141`), so her index hunk is inert on devel and her
  CHANGELOG line "no longer aborts setup" is false there. The ruling's stated mechanism was exactly
  the stale closing-edge index her hunk fixes; the rest was a scope fence (residue R-4/O3).
- The thin-shell μ0/linear-neighbor checks that 0.9.2 had commented out are back on devel
  (`8243c5a`, `77b4567`); her copper-neighbor decks are the first real test of the relative band.

## Jury round 1 — her diff (`review_yufan_ni_coil_fixes.md`)

Pre-registered: no P0; call placement, parity, conflict count, tie rule, flag hygiene and error
tier sound; B1 = the hunk is dead code under DR-188 (P1, records). Both auditors, independently:
**a split cross-type vote flips the minority so the lower domain type becomes master on part of a
sideset, and only a `message` follows (`cl_MaxwellFactory.cpp:1691-1699`); downstream,
`relink_slave_elements_with_duplicate_nodes` (`cl_CutFactory.cpp:1437-1479`) unions the slave
blocks of the whole sideset and relinks every flagged element, so a block that is slave on one part
and master on another is torn on both sides.** P1, confirmed in source. P2s: the `fix_facet_masters`
contract (`.hpp:228-231`), anchorless multi-patch sidesets, `thin_shell_facet_orientation.md`
still two-stage, stale comment at `cl_CurveFactory.cpp:112-115`. `git status` clean after both legs.

## Jury round 2 — the closed-loop policy (`review_closed_side_curve_policy.md`)

Christian asked why DR-188 refused closed loops and whether Yufan's fix makes the policy worth
changing. Literature read (local extractions): **Alves et al. 2022a §II-B** duplicate the shell
"except nodes and edges on its lateral boundaries ∂Γs", and their §III-B disc has a wholly interior
rim and needs no cut; **Schnaubelt et al. 2023** put the NI contact layer "in-between conducting
materials" fed by conductor leads. A closed free rim is published and validated; BELFEM's 3-D
exclusion of side-curve nodes from duplication (`cl_CutFactory.cpp:1299-1306`) is that rule.
Both auditors traced the same hazard: a thin-shell terminal curve lying on the excluded rim makes
`close_terminal_loops` (`:1601`) build the return segments on undoubled nodes; both segments land
on one mesh edge (`cl_EdgeFactory.cpp:382-392`, `cl_Curve.cpp:166-176`) with opposite
`edge_direction`, and `Chain::addSimplexToChain` (`cl_Chain.cpp:115-118`) erases the zero sum — the
current constraint vanishes silently. Reconciled: keep the hunk; do not lift wholesale; **narrow**
the refusal to a closed rim that shares a node with a terminal curve of its shell (a free edge and a
terminal are contradictory at one node). Christian ruled: narrow it. The refusal moves to
`CutFactory::create_side_curves_for_thinshells_3d` (`:1146`, before `duplicate_nodes_on_face_sidesets`
at `:1157`), the first place that holds side curves and terminal curves together.

## Plan and jury round 3 (`todo/yufan_ni_coil_fixes.md`, `review_yufan_ni_coil_plan.md`)

Christian's rulings (question form): block-set refusal after alignment; narrowed closed-rim
refusal; O2 (rim on a bulk terminal face) register row `[W]` only; O4 refusal in `CutFactory`;
Fable plans, plan+audit → switch to devel → code+audit; probe the cut collapse first.

Plan v1 → v2 after the jury; every finding re-read against the tree:
- split vote ≠ block intersection: two patches with distinct block pairs give disjoint sets and
  only the message signals it → keep the message (Codex, Grok);
- with no cross-type voters the **majority** wins, the seed only on a tie (`:1645-1646`);
- the exclusion predicate is the flag the 3-D branch **overwrites** (`:1258-1284`, "shares no node
  with another shell's side curve"), not the geometric one → G10 qualified to isolated rims;
- the relocated refusal would leak the side curves (`Protoshell` borrows, `Mesh` deletes,
  `cl_Protoshell.hpp:55-57`, `cl_Mesh.cpp:218-220`) → own them across the throw;
- node flag slot 0 is dirty when the check runs (`sort_end_nodes` flags and never clears,
  `cl_CurveFactory.cpp:1126-1147`) → `unflag_all_nodes()` first;
- `has_slave()`/`has_master()` guards in the block check (`align_thin_shell_facets` itself treats a
  slave-less facet as reachable, `:1657-1661`);
- the terminal-curve harness omits the homology section and stops at `cl_MaxwellFactory.cpp:1142`
  → the new refusals need a new fixture; the closed-rim refusal becomes a static helper on a
  `Protoshell *` so a unit test can call it;
- no positive mixed-winding regression existed → added (a hand mesh the walk reaches from two sides);
- "tape current = lead current" is not an NI invariant under radial bypass (Schnaubelt §IV) → the
  gate is terminal-in = terminal-out plus I_C2 against the integrated azimuthal current;
- **the `update_node_facet_tables` tNode/tOrg mismatch is an alias today**: master nodes are
  asserted non-duplicate (`cl_ThinShellFactory.cpp:701`), `original()` returns `this`
  (`cl_Node.hpp:317-321`); her devlog's crash mechanism is unproven → DR-194 filed as `[W]`;
- anti-parallel patches of **one** block pair are already refused by the block-set guard, so
  O1 "warn" was unreachable for that case; the propagation can only matter for distinct block pairs.
  Christian, re-asked with that fact: **implement the propagation here** → R2b (carrier BFS through
  non-shell facets with a parity-conflict refusal), plan v3, focused jury dispatched.

## Pancake probes (executed, private trees, release MKL, np = 4)

- Patched build on `cmake-build-debug/pancake_tape` (200 turns, conductor|conductor tape pieces,
  self-touching solder): `align_thin_shell_facets` printed nothing — every sideset already uniform,
  the pass a no-op. Then `CutProcessor::collect_facets` aborted: "Thin cut 0 collapsed … pruning
  removed 7860 of 7965 emitted faces in 32 passes". **Unpatched devel fails identically.**
- Diagnosis (Opus, read-only): `pancake.exo` of 09-24 holds only `NodeID`/`NodeOwner` — it is
  `msh2exo` output, **the deck never ran**; the guard predates it; the four devel commits since
  09-20 on the cut path are refusals/ownership only. 20-turn smoke case reproduces (3553/3606).
- Dangling-edge probe (approved, one file in a private worktree, reverted): first pruning pass,
  cut 0, 20 turns — 3606 emitted faces, 207 dangling edges, **207 on the air/conductor surface, 0
  interior**; 192 of the emitted faces there in air, 15 on the interface, 0 in the conductor. The
  sheet lives in the air and ends on the solder surface, which has **no sideset** in this generator;
  `cl_CutProcessor.cpp:395-398` exempts `mPhiBoundariesAndPeriodic` only and `Topology`
  (`cl_Topology.cpp:217-245`) types existing sidesets, creates none. Not the cohomology generator.
  Christian: fix the generator now in parallel (Opus, private copy, emit the solder/air contact as a
  sideset as `examples/pancake` does with `WallLeft_k/Right_k`); Gregory's three questions drafted
  in `tmp/gregory_note_cut_questions.md` for Christian to send.

## Generator fix and the kernel assertion behind the SEGV (executed, private trees)

Opus added a `SolderAir` sideset to a private copy of the `pancake_tape` generator (every
solder/air contact triangle, 6276 at 20 turns; air-facing tape facets and the terminals excluded;
left undeclared because `Topology` types a conductor|air sideset `InterfaceCondAir` on its own,
`cl_Topology.cpp:221-237`, `:342-348`, `:501`; a `validate` check covers the contact). 20-turn case:
cut creation completes, 0.14 s, no pruning line; the pristine generator on identical nodes still
collapses (3553/3606). Cut supports changed with identical nodes (1867/7822 → 1329/7311), so the
sideset also changes the cohomology input — noted for Gregory. Diff:
`scratchpad/pancake_gen/interface_sidesets.diff` (+71/−5, not yet in the tree).

The run then SEGVs in release at binary-dependent sites: unpatched devel right after RCM, patched
devel after `computing element connectivities ( thin shells )`. Patched debug build under gdb:
`DofData::create_dofwise_t_matrices_master`, `cl_FEM_DofMgr_DofData.cpp:3112`,
`ASSERT tNode->source(j)->number_of_dofs() == 1`, "couldn't figure out how to assign dependencies of
dof 1275943 (edge 561783)", rethrown from `Kernel::create_field` (`cl_FEM_Kernel.cpp:620`). A
debug-only guard; release proceeds into undefined behaviour, hence the moving SEGV. This is the
hanging-edge dependency path of `todo/test_gap_followups.md` R13, a kernel defect outside the
plan. Same debug run: **"thin shell tape: aligned 14 of 200 facets of sideset 83", "aligned 22 of
480 facets of sideset 88"** — Yufan's alignment pass fired on a shipped-generator deck, so her
T-junction case is real on our own meshes, not only on hers.

**Unpatched devel, debug build, same mesh, gdb:** a different and earlier abort —
`Curve::assign_edges` from `CutFactory::run` (`cl_CutFactory.cpp:135-138`), `ASSERT tEdge != nullptr`,
"Edge not found for segment 138168". With sidesets 83 and 88 mixed, the node duplication misroutes a
curve segment and no mesh edge matches it; the release SEGV after RCM is that path without the
guard. **Yufan's alignment pass removes this abort** (patched debug passes `CutFactory::run`, cut
creation 0.435 s), so her fix is verified positive on a shipped-generator deck, with unpatched
devel as the red control. The dof assertion in `create_dofwise_t_matrices_master` is what remains,
and it is a separate defect.

## The kernel assertion, diagnosed (Opus, read-only, gdb on the patched debug build; cites spot-checked)

Not the thin shell, not Yufan's patch, not the interface sideset as such. Mechanism:
`create_hanging_edges_and_facets` (`cl_MaxwellFactory.cpp:~1755-1830`) hangs the conductor-side
edges of every `InterfaceCondAir` facet on the air element's facet nodes, which after the cut relink
can be cut duplicates; `CutSet::create_duplicates` (`cl_CutSet.cpp:70-80`) makes each duplicate
depend on its original plus the λ node of every generator in its pattern; the edge branch of
`DofData::create_dofwise_t_matrices_master` (`cl_FEM_DofMgr_DofData.cpp:3094-3116`) then assumes
every such source carries exactly one dof, the scalar φ, and takes `dof(0)`. φ goes only to air,
abstract and orphan nodes (`count_node_dofs`, `:870-990`; `CutFactory::collect_orphan_nodes`,
`cl_CutFactory.cpp:2312-2342`, adopts block-less nodes only). The failing source is node **6276**
(r = 4.12 mm, on the LeadIn bar, block 42, interior node of the `SolderAir` sideset, on no tape
sideset): both cut sheets co-terminate on the interface edge 6571→6276, their plus sides cover all 8
of its air tets, so every air tet moved to a duplicate (24175: λ₀, 25930: λ₁, 26164: λ₀+λ₁), the
original keeps only its 12 conductor tets, gets no φ and is not an orphan — `ndofs = 0`. Debug
asserts; release reads a null `mDofs` and SEGVs in the same function. A second, silent effect: the
node branch (`:3031-3086`) un-hangs a node dof whose source has no matching dof, so in a release
build that happens not to crash the λ jumps at such a node are dropped. Pre-existing devel gap on
the cut / h-φ path, triggered by two or more cut sheets that co-terminate on an interface edge; a
single sheet always leaves the original a minus side. Not the ±1-weight refusal of
`todo/test_gap_followups.md` R13. Alignment flips as the cause: unlikely (the node is on neither
83 nor 88) but not excluded, and not testable on unpatched devel, which aborts earlier in
`assign_edges` on this mesh. 0.9.3 blocker: probably not — no shipped deck is known to have
co-terminating sheets (medium; none was run). Fix direction, for a plan and jury, not a patch: an
interface original with cut duplicates and no air element still needs a φ dof (adopt it like an
orphan). → DR-196.

## Implementation R1–R5 (Opus, scratchpad devel worktree, unstaged on top of her staged patch)

15 files, +382/−84 over her patch: R2b block-consistent alignment inside `align_thin_shell_facets`
(patches grouped by unordered block pair, reference = first anchored else smallest facet id,
whole-patch flip; a patch whose facets do not share one pair is left alone — the plan's debug assert
was dropped because a sheet between block 1 and blocks {2,3} is legal input); R2
`check_thin_shell_block_sides` (both-sides refusal, bitset intersection refusal); R3
`CutFactory::refuse_closed_side_curve_with_terminal` (static, public for the unit test; unflag
first, delete the shell's curves on the throw, immediate push into the mesh, late push loop
deleted); the `CurveFactory` refusal removed; messages, header contract, orientation doc (three
stages, stale range replaced; Codex sweep still owed), CHANGELOG, DR-194/195 filed; tests:
`ThinShellClosedSideCurveIsBuilt`, `test_CutFactory_closed_rim.cpp` (4), `test_ThinShellAlignment.cpp`
(6, full `MaxwellFactory` through `create_magnetic_kernel` on a Kuhn-tet box; the fallback was not
needed). Every refusal shown red with its guard removed. Debug build zero warnings; `test_mesh`
130/130, `test_homology` 29/29, `test_fem` 316/316; **`make check` 17/17** in the private debug tree
(MKL, Armadillo). By-catch: `~CutFactory` double-deletes temporary sideset facets on unwind →
DR-197. DR-196 filed from the kernel diagnosis. Code jury dispatched on `yufan_ni_coil_code.diff`
(astra/high + grok-4.7/high), R2b's design included.

## R7 and the doc sweeps (Opus, records only)

`examples/pancake_tape/` added: the generator with the `SolderAir` sideset (diff-identical to the
verified private copy), `Allrun`/`Allclean` wrappers, the 200-turn deck, a README whose closing
status says the example does not yet run to a solution and describes the kernel defect in words
(no register id in a doc); `examples/README.md` (twenty-one decks), CHANGELOG `### Added`,
`todo/pancake_inner_leadout.md` repointed; generator smoke run 6.8 s, validation line "SolderAir:
6276 facets, each between one solder and one air tet". By-catch, not fixed: the shared launcher's
hint `( cd ./python && python3 main.py )` writes the mesh into `python/` for every Python-generated
deck (`examples/scripts/Allrun`); `*.exo` and `*_RRR*.hdf5` are not git-ignored anywhere.
`thin_shell_facet_orientation.md`: Codex luna/medium sweep, 36 of 40 edits applied (4 skipped for
changing meaning), and six pre-existing factual errors fixed against the source (`DomainType`
values, the element-ID master rule applies to newly connected facets only, same-type flagging is
3-D only, `update_facet_nodes` skips thin-shell sidesets, conflict detection is in
`align_thin_shell_facets`, stale `Node::facet()` pointers last until `connect_nodes_to_facets`).
Residue: three pre-existing `todo/orientation_fix_block_ordering_bug.md` citations in that doc
break the docs-cite-no-todos rule; left for a records pass. The pancake README sweep is running.

## Codex quota exhausted mid-round

At 17:18 PDT Codex returned "You've hit your usage limit … try again at 8:30 PM". The pancake README
sweep died (0-byte output) and the code jury's Codex leg (`review_yufan_ni_coil_code`, astra/high)
died silently the same way — `_cross_review.log` carries the error, the leg's file is empty while
Grok's leg keeps running. A detached runner (`scratchpad/codex_rerun_2031.sh`, `setsid nohup`) fires
at 20:31 PDT and re-runs both: the jury leg with the same brief `cross_review.sh` builds, appended
to `review_yufan_ni_coil_code.md` as a `# CODEX` entry, then the README sweep. Until it lands, the
code round is Grok + Fable only, and the commit waits on Christian's ruling whether to hold for the
Codex leg or to take it post-commit.

## Code jury verdicts and the fix round (2026-10-06 → 10-08)

Grok 4.7/high (10-06 17:28): the vote and the block-set guard sound; **P1** the one-block remedy
text ("split the sideset") cannot fix a shell inside one volume block; **P1** the closed-rim
refusal checks only the shell's own terminals, so a rim touching another shell's terminal keeps a
node single-valued where a duplicate is needed, and a rim touching another shell's side curve is
rewritten open and torn by the duplication pass; P2s: split-vote message overstates, "unrelated"
count conflates two causes, edge-coating message claims closure, two doc ranges off by one,
flag slots 1–3 left set on the new throw paths, no split-vote / mixed-pair tests, 2-D outside the
check (by design). Codex astra/high (rerun 10-06 20:32 after the quota reset, on the pre-fix diff):
the diff file lacked the two untracked test files (packaging, not code); **P2** a branched cyclic
rim now reaches `sort_end_nodes` and throws past `tReleaseScratch` (leak in throw mode); P2 the
O(P²) patch grouping. Every finding verified against the source by Fable; no P0. Christian ruled
(10-06): fix everything now; a closed rim is accepted only when isolated (no node shared with any
shell's terminal or another shell's side curve); re-audit = the Codex rerun plus Grok 4.7/xhigh on
the fix hunks. Fix round dispatched 10-08 12:50 (Opus): refusal scope across all protoshells with
ownership handled at the throw, a degree-3 branch check before any curve is built, remedy and
message texts, flag hygiene before the throws, group-once, the two extra fixture cases, doc
ranges, and the pancake README sweep (33 edits, no technical errors flagged).

## Fix round 1 landed (Opus, 2026-10-08)

Closed-rim refusal rescoped to "isolated": `CutFactory::refuse_closed_side_curve_sharing_nodes(
Mesh *, Cell< Protoshell * > & )` flags every shell's terminal nodes (pass 0) and every other
shell's side-curve nodes (pass 1); all side curves are built first, checked once, then pushed; a
refusal deletes every shell's still-unowned side curves before the throw. Branch refusal in
`thin_shell_side_curves` before any curve exists ("branches at node %lu"), scratch released
first. Remedy text: split the volume block, not the sideset. Split-vote and multi-patch messages
state facts only (four counts: aligned, flipped, unique-pair, mixed-pair). Flag slots 1–3
cleared before both refusals in `align_thin_shell_facets`. Group-once patch grouping. Edge-coating
anchor message no longer claims closure. Tests: closed-rim suite rebuilt on two cubes / two shells
(other shell's terminal → refused; other shell's side curve → refused; disjoint closed rims pass),
`ThinShellBranchedSideCurveIsRefused`, `SplitVoteIsRefusedByTheBlockSets`,
`MixedPairPatchIsAccepted`; `test_mesh` 131/131, `test_homology` 32/32, `test_fem` 318/318, zero
warnings, each new refusal red without its guard; **`make check` 17/17** again in the private debug
tree. README sweep applied (29 of 33; three skipped for changing a fact, one no-op). Residue stated
by the implementer: a `thin_shell_side_curves` throw on a later shell leaks earlier shells' unpushed
curves in throw mode (same class as DR-188 (b)); the parity-conflict refusal has no fixture. Grok
grok-4.7/xhigh is auditing the fix-round diff; the Codex seat is the 10-06 rerun on the pre-fix diff.

## Grok xhigh on fix round 1, and fix round 2 (2026-10-08)

Grok grok-4.7/xhigh verified all round-1 fixes (reference rule shown equivalent to the old one;
same-shell contact of a closed rim is caught by the degree-3 check; pass 1 correctly skips the
shell's own side curves) and found **one regression from the round**: building every shell's side
curves before pushing any means a throw in `thin_shell_side_curves` for a later shell leaks the
earlier shells' curves in throw mode, where the per-shell push had given the mesh ownership before.
P3s: the unique-pair grouping scan still Θ(P²); a sweep edit made the README say BELFEM cuts the
solder; DR-195 silent on the 3-D-only scope; no test prints a positive mixed-pair count. Fix round 2
dispatched (Opus): per-shell push restored, the sharing refusal shrinks `mMesh->curves()` back
before deleting, a Map by ordered pair for the grouping, the README subject, the DR-195 clause, and
a two-piece mixed-pair fixture if the Kuhn box can build it. The Codex seat (astra/high) audits the
round-2 diff before the commit.

## Fix round 2 landed (Opus, 2026-10-08)

Per-shell push restored in `create_side_curves_for_thinshells_3d`; `refuse_closed_side_curve_sharing_nodes`
takes the pre-push curve count and, on a refusal, shrinks `mMesh->curves()` back (`Cell::set_size`
keeps the prefix) before deleting the side curves, so no curve is ever deleted while the mesh holds
it; red without the shrink: count 2 ≠ 1, then a double free in `~Mesh`. Correction to Grok's xhigh
finding, checked against `git show HEAD:…`: devel also built all shells before pushing any, so the
leak on a later shell's build throw was pre-existing, not a round-1 regression; closed either way.
Grouping by `Map< luint, index_t >` keyed on the ordered block pair; README subject restored
(the generator cuts the solder); DR-195 names the 3-D-only scope; `MixedPairPieceIsCountedSeparately`
asserts a 1 unique / 1 mixed line. `test_mesh` 131/131, `test_homology` 32/32, `test_fem` 319/319,
`check_doc_claims` 38/38, **`make check` 17/17** (third run). Codex astra/high audits the round-2
diff; then Christian's sign-off and the commit.

## Evidence ladder

Executed: the trial merge; three probe builds and five deck runs (200-turn patched/unpatched,
20-turn with the interface sideset on patched release, unpatched release, patched debug and
unpatched debug under gdb); the dangling-edge probe; the kernel diagnosis under gdb; the three test
binaries after each round; `make check` 17/17 three times in Claude's private debug tree (MKL
flavor, Armadillo); every new refusal red with its guard removed. Reviewed only: the parity-conflict
refusal, the two flag-clear throw paths, the later-shell build-throw path, the MPI debug throw
window, the Alves/Schnaubelt §-quotations (local extractions, no auditor could check them).
Not run: a release-tree `make check`, Blaze, valgrind on the new paths, the copper-lead end-to-end
gate (blocked on DR-196).

## Records touched this session

`todo/yufan_ni_coil_fixes.md` (new, v3), `todo/README.md` (registered), this devlog,
`tmp/gregory_note_cut_questions.md`, `tmp/ai_exchange/{yufan_ni_coil_fixes*.diff,
review_yufan_ni_coil_fixes.md, review_closed_side_curve_policy.md, closed_side_curve_policy_question.md,
review_yufan_ni_coil_plan.md, review_r2b_patch_propagation.md}`. Register rows DR-194/DR-195 are
R5 of the plan, not yet filed.

## Close-out

Committed on local devel 2026-10-08 as the single commit that carries this devlog, the plan in
`todo/closed/`, the source, the tests and `examples/pancake_tape` (36 files, +4908/−107; the first
commit id `c7472c4` was amended to add this paragraph). Not pushed: the Thickcuts session holds the
devel push and the devel→thickcuts merge for this SHA. Gregory's note is at
`tmp/gregory_note_cut_questions.md`, for Christian to send. Owed and recorded in the register: the
copper-lead end-to-end gate (DR-196), DR-194/195/197.
