# One Conductor-Side Interface Duplicate per Node (Yufan's Seam Short)

**Date:** 2026-10-08
**Purpose:** Close the lossless current paths Yufan found on the NI pancake (`bug_seam.md`,
2026-10-03; root cause in `devlog/dl20261008_seam_short_interface_duplicates.md`). Mechanism:
`InterfaceProcessor` creates the conductor-side copy of every conductor|air interface node once
per (conductor block, air block) *set*; a node where two conductor blocks meet the air gets one
copy per block, the conductor edges leaving it along the conductor|conductor face are created
twice (edge keys are node indices) and nothing ties them, so tangential h may jump across the
air-adjacent row of every conductor|conductor face — a lossless sheet. The fix makes the duplicate
a property of the node: one copy shared by every tie-weight-1 set, so adjacent conductor blocks
share it and the edges are created once.
**Module:** `src/homology` (`cl_InterfaceProcessor`, outside the closed cohomology core), `tests/homology`
**AIs involved:** Claude Fable 5.1 (diagnosis, plan, pre-registration, reconciliation); Codex gpt-6-astra/high and Grok grok-4.7/high (plan and code juries); Claude Opus 5.5 (implementation in a scratchpad worktree)
**Status:** ✅ COMPLETE 2026-10-08 — landed on devel in one commit (Christian's rulings of 2026-10-08). Landed: one conductor-side interface copy per node shared by every conductor-master set, keyed by the node's own id, plus a sweep over all conductor elements (R1, R4); periodic copies paired once with identity checks (R2); copies appended once (R3); `tests/homology/test_InterfaceProcessor.cpp`, eight tests, five red on HEAD (R5); gates G0–G5 as recorded under R6 — `make check` 17/17 zero warnings, F mesh strips and frames gone, F deck seam jumps 50 A → 1 A with the Ampère sum closing and I_2 charging, H deck unchanged to 0.1 % (R6); devlog, guide paragraph (Codex luna sweep), DR-198 struck (R7). Juries: plan round Codex astra/high ×2 + Grok 4.7/high; code round astra/high + Grok 4.7/high with one approved fix round; Codex terra/medium quick pass clean. Residue: no periodic deck run (the hand-built periodic fixture is the only unit gate for R2); O2 orphan copies and the Decouple source inheritance stay as noted; Yufan's note delivered beside her report.

> **Scope guards.**
> - In scope: the duplicate sharing in `InterfaceProcessor`, the per-facet master-block relink,
>   the periodic pairing that follows from sharing, a unit test with a positive control, and the
>   executable gates on Yufan's F and H decks.
> - Out of scope: the thin-shell cut relinking (`CutFactory`, confirmed correct at the seams);
>   the orphan copies the set machinery leaves behind (O2); `Decouple` (coil) semantics; the
>   cohomology core; 2-D (the same code runs, nothing 2-D-specific changes).
> - Promise kept: the visualization split between a conductor and the air stays (one conductor
>   copy per node, the air keeps the original); `TieWeight1` sourcing unchanged.

---

## 1. Current Behaviour and How It Fails

`InterfaceProcessor::create_interface_sets` (`cl_InterfaceProcessor.cpp:325-397`) keys an
`InterfaceSet` by the hex of a two-bit block bitset {master block, slave block} of each sideset's
**first** facet; sidesets between two air blocks, or with neither side air/ferro, are skipped.
`InterfaceSet::collect_nodes_and_elements` (`:66-146`) flags the nodes of the set's facets together
with their `original()` and all its duplicates, records the master block of each sideset's **first
facet** only (`:103-110`), and collects every element of those blocks that touches a flagged node.
`InterfaceSet::duplicate_nodes` (`:156-213`) allocates one new node per flagged original, **per
set**. `relink_elements` (`:229-252`, `insert_node` at `:243`) rewires the set's elements to the set's copies. The air side
(slave) keeps the originals; the copies hang on their originals with weight 1
(`add_duplicate_nodes_to_mesh`, `:516-575`), a visualization split by design (INC-453).

| Failure | Mechanism | Evidence |
|---|---|---|
| Lossless sheet along every conductor\|conductor face where it meets the air | node n on the line C1∩C2∩Air lies in set (C1,A) and set (C2,A); each makes its own copy; C1's elements get copy₁, C2's copy₂; the edges (m, copy₁) and (m, copy₂) to an interior node m of the C1\|C2 face are two edges (`cl_EdgeFactory.cpp:287-310`, index-keyed); the hanging-edge pass sources only edges **on** interface facets (`cl_MaxwellFactory.cpp:1994-2040`) | run F exodus: five copies at (r 14.750, θ 0°, z 1.50): blocks 10, 11, 12, 13, air 18; 6 non-conforming faces per radial sub-face at 0°, all touching the top/bottom edge; run H: 814 non-conforming frame faces on the mid-gap surface |
| Short between tape sheets | the strip connects two points of the shell: seam to seam along the 0°/180° radial faces (A, B, D, F), tape end to tape end via the liner faces and the mid-gap frame (I) | Yufan's §3 table; run G (seam without block face) is clean; run H (no cu\|cu face touches the tape) is clean |
| Same defect without `Mesh.SaveAll` | one air-interface sideset spanning two conductor blocks: only the first facet's block is relinked, the second block keeps the originals — doubled edges again, and the visualization split is missing on the second block | `:103-110` (first facet only); untested on Yufan's meshes, which name every surface |

**Bottom line:** the duplicate is keyed by (set, node) where the design needs it keyed by node;
every conductor block touching an interface node must reference the same conductor-side copy.

## 2. Architecture: Why Sharing the Copy Is the Right Spine

The duplication between domains of different material, or between a conductor and a void, is deliberate (Christian, 2026-10-08): h is continuous there but j is not, so the nodal j field must not be averaged across such a surface — the conductor needs *a* copy distinct from the air's original, not one per block. Two conductor blocks of the same material share their common face today (no interface set is made for a conductor|conductor sideset, `:336-341`); the fix only removes the accidental per-block split on the line where that face meets the air. A processor-level
map original-id → copy, consulted by every `TieWeight1` set before allocating, restores exactly
that, leaves `relink_elements`, `add_duplicates`, `unify_duplicates` and the hanging-edge sourcing
(`original()`-keyed, `cl_MaxwellFactory.cpp:2010-2014`) untouched, and makes the periodic pairing
simpler (dup(A) ↔ dup(B), one object each).

Rejected: (i) tying the doubled edges afterwards (an edge-level hanging relation between two
conductor edges) — a new mechanism for a defect the node level causes; (ii) keying sets by air block
only — still splits a conductor at the line where two *air* blocks meet it; (iii) dropping the
duplication for conductors — removes the visualization split Christian confirmed as intended.

## 3. Gap Table

| # | State | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| 1 | one conductor-side copy per interface node | shared edges along cu\|cu faces | no — per set | (c) R1 | `:156-170` |
| 2 | every conductor block of a sideset relinked | copy used by all blocks of a multi-block sideset | no — first facet's block | (c) R4 | `:103-110` |
| 3 | periodic pairing of copies | periodic decks (helix, corc) | per set + cross-set pass | (c) R2 | `:181-225`, `:452-505` |
| 4 | each copy appended to the mesh once, sourced once | mesh integrity | per set | (c) R3 | `:516-575` |
| 5 | `Decouple` (coil) copies private, unsourced | coil decoupling | yes | (a) unchanged | `:19-35`, `:255-268` |
| 6 | thin-shell factory tolerant of fewer copies per original | NI decks | unknown by reading: `collect_nodes` keys by `original()->index()` (`:826-838`); the duplicate walks are the normal assembly (`:1102-1113`) and `update_node_facet_tables`; the `create_nodes_on_layers` link never fires on this path (flags cleared at `:164-167`, Grok) | (b) O1 → gate G3/G4 | see O1 |
| 7 | `.bfm` cache | every rerun | the `.bfm` is the **processed** mesh, written after the hanging edges (`cl_MaxwellFactory.cpp:720-740`) and reused whenever the base-mesh checksum and the deck tag match (`:400-455`); a code-only change moves neither key, so a sibling `.bfm` replays the old enrichment | (c) G2–G4 delete it first; Yufan's note says so (Grok) | — |
| 8 | parallel | all | built on the mesh master (`cl_CutFactory.cpp:115`), distributed later by `Kernel::distribute_mesh` (`cl_FEM_Kernel.cpp:111`); no collective in this code | (a) untouched | — |

### 3.1 Cross-cutting findings

- The set key stays as the grouping for *treatment* and periodic bookkeeping; only the **copy
  identity** becomes global. A node in both a `TieWeight1` set and a `Decouple` set keeps a
  private decoupled copy for the coil side (O3).
- `collect_nodes_and_elements` runs for all sets before any relink (`populate_interface_sets`,
  `:424-431`), so every set sees the originals; sharing does not change the collection.

## 4. Ordered Steps

- [x] **R1** One copy per node for the **Conductor** master class (v4: Ferro and every other master type keep private per-set copies as today — ferro carries node dofs, so the doubled-edge mechanism cannot occur there, and a ferro|ferro sideset is an admitted set whose material split is deliberate; a refusal of other master types would abort Buffer|Air and Air|Default sets of shipped decks). `InterfaceProcessor` owns one shared map, keyed by the collected node's **own id** (`tPair.first`, never `original()->id()`: `mOriginals` holds root nodes and pre-existing cut/thin-shell duplicates side by side, `cl_InterfaceProcessor.cpp:75-91`, and canonicalizing through `original()` would collapse distinct cut branches — Codex P2). A `TieWeight1` set's class is the domain type of its master block (Conductor for conductor|air and conductor|ferro sets, Ferro for ferro|air sets; `en_DomainType.hpp:23-30`, master = higher type). `InterfaceSet::duplicate_nodes( id_t & aMaxNodeID, Map< id_t, Node * > & aShared )` looks the node up in its class map first, allocates and registers on a miss; on a hit it still writes `mDuplicates[ id ]` **before** the periodic loop (`relink_elements` `:243` and the cross-set pass `:465` index it unconditionally — Grok); only the miss path allocates or bumps `aMaxNodeID`; the maps borrow mesh-owned nodes (the mesh owns them after `append`, `:573`; no destructor deletes map values). `Decouple` sets allocate privately as today. Conductor and ferro therefore keep separate copies on a conductor∩ferro∩air junction line, exactly as today (Codex P2: equal treatment is not equal visualization identity). The set's `mDuplicates` still maps node id → copy so `relink_elements` is unchanged.
- [x] **R2** (after R1) periodic pairing: in `duplicate_nodes`, guard on `Node::is_periodic` (not on flag slot 6, which is cleared per set, `:191-224`): pair C ↔ D only when neither is periodic yet (a shared copy may already be paired by an earlier set) so `Periodicity::add_node_pair_to_backup` — which appends without deduplication, and whose restore consumes every entry (`cl_Mesh_Periodicity.cpp:312-326`, `:275-290`) — is called **once per pair**; in `pair_cross_set_periodic_duplicates`, candidates are searched only among sets of the **same treatment and class**, a partner copy found in several such sets is ambiguous only when the objects differ, and the reachable failure — C not periodic, D periodic with another object (`:469` skips a periodic C before `:494`, Grok) — stays a `BELFEM_ERROR` judged within one treatment and class; a copy met already periodic is checked, not skipped: its partner must be the expected D or the configure refuses (`BELFEM_ERROR`, setup code — Codex, second entry). A `Decouple` copy pairs with a `Decouple` copy only (Codex P2). `test_CohomologyPeriodic` builds the complexes directly and never runs this path, so R5 (e) is the only unit gate for R2.
- [x] **R3** (after R1) `add_duplicate_nodes_to_mesh`: collect each copy once — iterate `mSharedDuplicates` for the tie-weight copies and the `Decouple` sets' maps for the private ones; `tOriginalMap` as today. `number_of_nodes()` is no longer a valid count of new nodes — size from the collected cell.
- [x] **R4** Class-level relink sweep (after the per-set relinking, which stays and installs the same shared copies): v4: conductors only — for every element of every block of that domain type (`Topology::groups( Conductor )`, `groups( Ferro )`), replace any node found in the class map by its copy. This reaches a conductor block that touches the interface node without owning an interface facet (Codex P1: C2 shares the conductor face through n but meets the air only at a node or an edge) and a sideset that spans two conductor blocks (today only the first facet's block is relinked, `:103-110`). `collect_nodes_and_elements` additionally takes the master block of every facet, not of the first facet only (`tBlockIDs` sized by the facet count before `unique`, not by `mSideSets.size()` as at `:96` — Grok), and asserts that every master block of one set has the same domain type, so the per-set collection is complete too.
- [x] **R5** (after R1–R4) unit test `tests/homology/test_InterfaceProcessor.cpp`, registered in the explicit source list of `tests/homology/CMakeLists.txt:9-14`: a hand-built tet mesh with two conductor blocks sharing a face and one air block above both, (a) one air sideset per block, (b) one air sideset over both, (c) a third conductor block that touches the triple-line node only through an edge of the conductor|conductor face and owns no air facet (Codex P1), (d) a ferro block beside a conductor block under the air (conductor and ferro copies distinct, both distinct from the air original — Codex P2), (e) a periodic variant built the way `test_CohomologyPeriodic.cpp:160-176` builds its `Periodicity` by hand, with a periodic pair whose partners lie in different sets: one backup entry per pair, no ambiguity refusal, and a genuine conflict still refused. Assert for a node on the triple line: exactly one conductor-side copy, referenced by elements of every conductor block touching it, the air elements on the original, the copy hanging on the original with weight 1; assert the face shared by two conductor blocks references identical node ids in both elements. Positive control: (a), (b), (c) red on HEAD.
- [x] **R6** gates (Claude's private debug worktree, Christian's ruling 2026-10-08). **Results 2026-10-08, fixed binary (plan v4 diff), fresh directories:**
  - G1 `make check` 17/17 (debug, MKL, Armadillo), zero warnings; `test_homology` 40/40 — on the v4 diff and again on the final diff after the fix round.
  - G2 F deck first save: zero non-conforming faces between copper blocks other than the 600 tape facets per turn pair (HEAD: 6 per radial sub-face pair at 0°, 200–208 per mid-gap pair); nodes 20463 → 19632. H deck: mid-gap frame 816 → 0; nodes 19223 → 18467.
  - G3 F deck at 5 s: seam jumps 1.05 / 0.16 / 0.09 / −1.07 A (HEAD at 5 s: 50.2 / 8.8 / −6.6 / −51.2 A, I_2 = −43.4 A; Yufan's 0.9.200 run at 5 s: 59.1 / 9.5 / −6.6 / −60.3 A, I_2 = −43.2 A, with the same windows); net tape current 113.0–113.2 A at 90°/175°/185°/270° against I_2 = 112.6 A (HEAD: 73–75 A against I_2 = 41.7 A — the mid-gap frames carried azimuthal current the tape sum cannot see); I_2 rises −28 → −113 A over 1–5 s where HEAD sits at −42 A. The 0° windows are not used for the Ampère sum: the lead arcs overlap turns 0 and 4 there in this diagnostic.
  - G4 H deck: I_2 within 0.1 % at every save (−37.55 vs −37.52 A at 5 s); sheet currents at 270° within 0.7 A, and the Ampère sum now closes (37.52 vs 37.52 A) where HEAD was 2.3 A short — the frame loop had carried that; the 0.5 A bar of v3 assumed HEAD was right.
  - G5 `test_homology` 40/40 incl. `CohomologyPeriodic`; no periodic deck run (residue).
  Original gate text:
  - **G1** `make check`.
  - **G0** every gate run starts in a directory without `ni_coil_test.bfm` (gap 7); the HEAD baselines of 2026-10-08 ran on fresh directories.
  - **G2** mesh-level: `tmp/seam_check/conformity.py` (which now excludes the faces of the thin-shell layer block) on the first exodus save of the F deck — on HEAD the (10,12), (11,13), … 6-face strips and the (10,11), … frames are present (positive control, measured: 6 faces per radial sub-face pair, 200–208 per mid-gap pair); with the fix, zero non-conforming faces between two conductor blocks.
  - **G3** physics, F deck (`geo_thick/ni_coil_flip_no180.py`, 5 s ramp as shipped), measured with `tmp/seam_check/sheet_currents.py` (±3° windows at 355° and 5°): on Yufan's own run F output at 5 s the seam jumps are 59.1 / 9.5 / −6.6 / −60.3 A and the tape's net current near 0° collapses to ≈ 2 A against ≈ 82 A at 90°; her ±6° measure gave 112 / 14 / 8.5 / 113 A (`bug_seam.md` §3, in the gitignored archive `cmake-build-debug/tmp_20261008_for_Christian/archive_20261005/`). HEAD must reproduce that picture (positive control); the fix must show jumps ≤ 1 A at every seam (run H: 0.54 / 0.05 / 0.08 / −0.31 A) and a net current near 0° equal to the other angles within 2 %.
  - **G4** H deck on HEAD and on the fix: I_2 at every save within 1 % and every sheet current at 90°/270° within 0.5 A (the fix removes the dead mid-gap frame; bit identity is not the bar — the dof set changes). HEAD baseline at 4 s: 29.36 / −2.08 / −2.63 / −4.36 / 4.25 A at 90°, I_2 = −25.99 A, identical to Yufan's 0.9.200 run.
  - **G5** `test_homology` with the new suite, and `test_CohomologyPeriodic` unchanged; the periodic decks (`examples/helix`, corc) are not run here — they are the residue if G5's hand-built periodic variant is judged too weak.
- [x] **R7** docs (devlog, usage-guide paragraph swept by Codex luna; the register strike is owed at the commit): `src/homology/doc/homology_usage_guide.md` (interface processing paragraph) states the one-copy-per-node contract; register row filed and struck in the same session if the gates pass; note to Yufan.

## 5. Open Design Questions (not silently decided)

- **O1** Does any consumer rely on the rim nodes of a thin shell being *distinct* copies on the master and the slave block (today true by this accident)? Open by reading: the duplicate walks that matter are the normal assembly (`cl_ThinShellFactory.cpp:1102-1113` and the line2/tri6 twins) and `update_node_facet_tables`; Grok puts a ~40 % chance on rim normals changing, which a uniform facet list on every duplicate would renormalize away. Decided by G3/G4, not by reading.
- **O2** The set machinery also copies the *original* of a thin-shell duplicate even when no element of the set's block references it (the seam-line "()" copies in run F); they become orphaned nodes with a φ dof. Pre-existing, out of scope; note for a later tidy.
- **O3** A node in both a `Decouple` and a `TieWeight1` set: keep the decoupled copy private. Proposed; confirm in the jury. **Codex P2 (confirmed against `:156-170`, `:255-268`, `:538`):** "decoupled copies stay unsourced" is not an established invariant — `duplicate_nodes` copies an already-hanging original's sources into the new copy regardless of treatment, and the later `Decouple` exclusion only skips the weight-1 sourcing. Pre-existing and conditional (it needs a hanging node on a coil interface); recorded as a register row, not fixed here — a coil-side behavior change is outside this plan's scope.
- **O5** (Grok) R5's fixture must put the conductor element on `master()` and the air element on `slave()` by hand — the test does not run `fix_facet_masters`; layout (b) is red on HEAD only when `facets()(0)` belongs to the block the assertion treats as first.
- **O6** If a conductor|conductor interface between *different* materials is ever to be split for the same j reason, the class rule must become per material, not per domain type. Out of scope; noted so the class keying is not mistaken for a material rule.
- **O4** RESOLVED 2026-10-08 (v4) → the only shared class is Conductor. The v3 claim that no admitted set has a master other than Conductor or Ferro was wrong: `mIsAirBlock` is seeded from `groups( Air )` only, so a Buffer|Air sideset is admitted with the buffer as master, and an Air|Default sideset with air as master; both keep today's private per-set copies, not a refusal. Ferro sets keep private copies too (ferro|ferro is admitted and its split is deliberate). Only the mixed-type refusal stays: the master blocks of one set must share one domain type.

## 6. Interface

`InterfaceSet::duplicate_nodes( id_t & aMaxNodeID, Map< id_t, Node * > & aSharedDuplicates )`;
`InterfaceProcessor::mSharedDuplicates`. No deck key, no input-contract change.

## 7. Definition-of-Done Checklist

- [ ] Every gap-table row mapped to a step or an open question.
- [ ] Unit test red on HEAD, green with the fix, for both sideset layouts.
- [ ] G1–G4 run and recorded with numbers.
- [ ] Devlog, register row, Yufan's note.

## 8. Audit Trail

- Exchange thread: `tmp/ai_exchange/review_interface_duplicates_shared.md`.
- Diagnosis devlog: `devlog/dl20261008_seam_short_interface_duplicates.md`.
