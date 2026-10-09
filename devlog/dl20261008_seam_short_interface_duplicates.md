# Devlog 2026-10-08 — Lossless current paths at conductor|conductor faces: the interface duplicates

**Date:** 2026-10-08
**Topic:** Root cause of Yufan's open defect "lossless current paths where a conductor|conductor block boundary touches a thin shell" (`bug_seam.md`, 2026-10-03)
**AIs involved:** Claude Fable 5.1 (diagnosis, plan, pre-registration, reconciliation, gates); Codex gpt-6-astra/high (plan round twice — the jury's own retry after a usage-limit error and a re-dispatch — and code round); Grok grok-4.7/high (plan and code rounds); Claude Opus 5.5 (implementation and tests in Claude's private worktree)
**Claude Confidence:** high — the mechanism was measured on her output files, reproduced on devel HEAD, and the fix removes it on the same mesh (gates below)
**Literature References:** N/A (mesh plumbing, not physics)

## Summary

Yufan's NI pancake runs A–I showed a thin-shell current that transfers between tape sheets without passing through copper wherever a copper|copper block face ends on a tape seam or a tape end. Her report ruled out node, facet and element relinking of the thin-shell cut and expected the defect at the edge/dof level. The defect is one level earlier than the thin shell: the **interface duplicates** made for the conductor|air surfaces by `InterfaceProcessor` are created **per (conductor block, air block) pair**, so a node on the line where two conductor blocks meet the air gets one conductor-side copy per block. The conductor edges that leave such a node into the interior along the conductor|conductor face are then created twice (edge keys are node indices, `cl_EdgeFactory.cpp:287-310`), and nothing ties the two copies together: the hanging-edge pass only sources the edges that lie *on* interface facets. Tangential h is therefore free to jump across the strip of conductor|conductor faces adjacent to the air, which is a lossless sheet, exactly the device the thin shell itself uses for its sheet current.

Where that strip joins two parts of the tape, the sheet current takes it:

| run | strip | effect |
|---|---|---|
| A, B, D | top and bottom rows of the radial copper faces at 0° and 180°, from seam P[k] to P[k+1] | sheets exchange current radially along the split planes (sheet 1 → 5 at 180°) |
| F | radial faces at 0° only | 0° jumps at P[1]..P[4]; the middle turns are bypassed |
| G | tape seam at 180° without a block face | no strip at 180°, no short — the seam itself is innocent |
| H | the only cu\|cu face is the mid-gap surface; its frame (top, bottom, both ends at 0°) is a strip but touches no tape | correct NI charge |
| I | as H plus two Out\|In faces at 0° that end on the tape's free ends P[0], P[5]; their top/bottom rows connect P[0] → mid-gap frame → P[5] | shortcut between the tape ends |

This is the trigger rule of her §1 stated mechanically: not "a block boundary touching the tape" but "a conductor|conductor face that reaches the air and whose air-adjacent row connects two points of the shell".

## Evidence

Measured on the enriched meshes BELFEM wrote (`ni_coil_test.e-s.00001` of runs F, H, I under
`cmake-build-debug/tmp_20261008_for_Christian/archive_20261005/`), with `tmp/seam_check/conformity.py`
(groups coincident node positions, then flags every geometric tet face whose two volume elements
reference different node ids) and `conformity2.py` (classifies the duplicated positions and lists
the faces of a named block pair with every copy at each position).

- Every air/copper interface position is multi-copy in all three runs (0 single-copy of 1448): the
  air keeps the original, each conductor block its own copy — the per-set relinking of
  `InterfaceSet::relink_elements` (`cl_InterfaceProcessor.cpp:229-252`).
- Run F, the mid-gap edge at the top of the 0° plane, position (r 14.750 mm, θ 0°, z 1.50 mm):
  five copies, one each for blocks 10, 11, 12, 13 and air 18. The interior node below it
  (z 1.00 mm) is one node shared by 10, 11, 12, 13. The radial face 10|12 (gap 1 | gap 2, inner
  sub-layer) has 6 non-conforming faces, all with a corner on the top or bottom edge; 11|13 the
  same; likewise 12|14, 13|15, 14|16, 15|17. Those 6-face strips are the only non-conforming
  faces on the 0° plane.
- Run F, the seam line P[1] at (14.500, 0°, 1.00): three copies — bore air, block 11 (inside turn 2),
  and one shared by blocks 10 and 12 (outside turn 1 and outside turn 2). The thin-shell relinking
  is correct at the seam; Yufan's §4 is confirmed.
- Every run: the mid-gap surfaces carry a non-conforming frame (F: 208 faces per gap on 10|11 etc.;
  H: the 814 faces of 31|32 with fewer than three duplicated positions). In H it is a dead loop
  spanning a surface with negligible flux, which is why H passes.
- In H the top and bottom faces of each conductor block against air are also non-conforming
  (30|31: 824, 30|32: 1352), by design: those edges are hanging edges sourced on the air nodes, and
  the conductor copies hang on their originals with weight 1 (`InterfaceTreatment::TieWeight1`).

Design intent of the duplication (Christian, 2026-10-08): the split between domains of different
material, or between a conductor and a void, is deliberate — h is continuous there but j is not, so
the nodal j must not be averaged across the surface (INC-453 had recorded it as visualization
only). Nothing in that intent needs one copy per conductor block of the *same* material: those
blocks share their common face everywhere except, by this accident, on the line where it meets
the air.

## Implementation and gates (same session, after Christian's approval)

Plan `todo/interface_duplicates_shared.md` (v1 → v4 across the rounds). Landed in Claude's private
worktree, not yet committed: `InterfaceProcessor` keeps **one conductor-side copy per interface
node**, shared by every tie-weight set whose master block is a conductor and keyed by the node's own
id; a sweep after the per-set relinking puts that copy into every element of every conductor block
that still holds the original, so a block that meets the air only at an edge is covered; periodic
copies are paired once per pair, with identity checks instead of silent skips; copies are appended
to the mesh once. Ferro, buffer, air-master and decoupled (coil) sets keep private per-set copies as
before — ferro carries node dofs, so the doubled-edge mechanism cannot occur there, and a ferro|ferro
split is a deliberate material split. New `tests/homology/test_InterfaceProcessor.cpp` (eight
tests on hand-built tet junctions, incl. a periodic pair across sets, a ferro|ferro split and a
buffer master): (a)–(e) red on HEAD, all green with the fix.

Juries: plan round Codex astra/high ×2 + Grok 4.7/high (two P1s folded in: share per master class,
reach blocks without an interface facet; Grok's `.bfm` finding below); code round Codex astra/high +
Grok 4.7/high (one P2 → the per-facet master walk reverted in favour of the sweep; Grok's P0 "the
periodic tests do not compile" refuted by the build and an 8/8 run). Thread:
`tmp/ai_exchange/review_interface_duplicates_shared.md`.

Gates, fixed binary, debug, MKL, Armadillo, fresh directories (the `.bfm` is the *processed* mesh
and is reused whenever the base-mesh checksum and the deck tag match — a changed binary moves
neither, so a stale `.bfm` replays the old enrichment; delete it before any rerun):

| gate | HEAD | fix |
|---|---|---|
| `make check` | — | 17/17; `test_homology` 40/40; zero warnings |
| F mesh, non-conforming copper\|copper faces off the tape | 6 per radial sub-face pair at 0°, 200–208 per mid-gap pair | 0 (nodes 20463 → 19632) |
| F deck 5 s, seam jumps P[1]..P[4] | 50.2 / 8.8 / −6.6 / −51.2 A (Yufan's 0.9.200 run: 59 / 9.5 / −6.6 / −60 A) | 1.05 / 0.16 / 0.09 / −1.07 A |
| F deck 5 s, Σ sheets at 90°–270° vs I_2 | 83 A vs 43.4 A | 113.0–113.2 A vs 112.6 A |
| F deck, I_2(t) | stalls near −43 A | −28 → −113 A over 1–5 s, charging |
| H mesh, mid-gap frame faces | 816 | 0 (nodes 19223 → 18467) |
| H deck 5 s, I_2 | −37.55 A | −37.52 A; Σ sheets = I_2 to 0.01 A (HEAD 2.3 A short: the frame loop) |

The 0° windows of the sheet-current diagnostic are not Ampère-valid (the lead arcs overlap turns 0
and 4 there); the seam jumps use windows clear of the arcs. Positive control: devel HEAD on the
regenerated F mesh reproduces her run.

## Proposed fix (as first written, before the plan rounds)

In `InterfaceProcessor`, make the duplicate a property of the original node, not of the set:
one processor-level map original → duplicate shared by every `TieWeight1` set, consulted by
`InterfaceSet::duplicate_nodes` before allocating. Adjacent conductor blocks then share the
conductor-side copy, the edges along their common face are created once, and the strip closes.
`Decouple` (coil) sets keep private copies. `pair_cross_set_periodic_duplicates` simplifies to
dup(A) ↔ dup(B) and its "duplicated in more than one set" refusal becomes moot. Second, let
`InterfaceSet::collect_nodes_and_elements` take the master block of every facet, not of the first
facet only, so a single air-interface sideset spanning two conductor blocks (a mesh written without
`Mesh.SaveAll`) relinks both; today the second block keeps the originals and the same doubled edges
appear.

Gates for the fix: `make check`; a doubled-edge count on the conductor blocks (edges whose nodes
map to the same original pair, excluding thin-shell and cut edges) that must drop to zero on
Yufan's F mesh; the F deck (regenerated with `geo_thick/ni_coil_flip_no180.py`) showing no 0° jump
in the sheet currents; the H deck unchanged. Her §1 checks (Ampère sum vs I_2, I_2 → N·I_1) are the
physics acceptance.

Open: whether the thin-shell factory relies anywhere on the rim nodes being distinct copies on the
master and slave blocks (today they are, by this accident; after the fix the rim copy is shared,
as the side-curve design intends). Read `update_node_facet_tables` before implementing.

## Files

- `devlog/dl20261008_seam_short_interface_duplicates.md` (this file)
- `tmp/seam_check/conformity.py`, `tmp/seam_check/conformity2.py` — the checkers (ephemeral; a
  conformity test "no non-conforming faces outside thin shells and cuts" would be the durable form)
- `src/homology/cl_InterfaceProcessor.{hpp,cpp}`, `tests/homology/test_InterfaceProcessor.cpp`, `tests/homology/CMakeLists.txt` — landed on devel (byte-identical to the juried diff `tmp/ai_exchange/interface_duplicates_shared.diff`), Christian's ruling 2026-10-08
- `todo/interface_duplicates_shared.md`, `todo/debt_register.md` (DR-198), `tmp/seam_check/sheet_currents.py`
- note to Yufan: `cmake-build-debug/tmp_20261008_for_Christian/archive_20261005/note_to_yufan.md`

## Addendum, same evening: `examples/helix/Copper_RRR50.hdf5` removed from history

Christian's ruling. The 1 MB copper-resistivity cache had entered the tree in Yufan's commit
f29206f (2026-10-05) and was an ancestor of `main`, `devel` and `yufan` on GitLab, though absent
from every current tree; it was also sitting staged in the shared checkout. Done: unstaged and
deleted; full bundle backup of all refs (`belfem_pre_rewrite.bundle`, 66 MB, in Claude's scratchpad);
`git filter-repo --invert-paths --path examples/helix/Copper_RRR50.hdf5 --force`; verified that no
ref and no object holds the blob; `origin` re-added; `main`, `devel`, `yufan` force-pushed with
leases on the old SHAs (`fb543a6` → `f8abe20`, `c2adb3c`/`c182960` → `f9b9421`, `ca331d2` → `b26b660`;
`origin/thickcuts` did not contain the commit and is unchanged). Every SHA cited in devlogs and
memory for commits after 2026-10-05 is stale; the bundle maps them. Every other clone — Yufan's,
Gregory's, the CI runner's, the two sibling worktrees — must hard-reset or re-clone; GitLab keeps
the old objects until its housekeeping runs. The uncommitted fix of this session was backed up
file by file before the rewrite and restored byte-identical afterwards.

## Addendum, late evening: DR-194, DR-197, DR-195 and DR-113 closed; the release gate started

On Christian's autopilot ruling the two latent rows of the NI-coil campaign were fixed through the
standing loop (plan `todo/closed/dr194_dr197_small_fixes.md`, exchange
`tmp/ai_exchange/review_dr194_dr197.md`, juries astra/high + grok-4.7/high on plan and code, no
P0; three text-level points of the code round applied: guard message, header contract sentence,
license block). DR-194: `ThinShellFactory::update_node_facet_tables` walks the duplicates of the
original in all five loops and refuses, in setup code, to add a facet to a full container
(`Vertex::facet_capacity()` new); the function is `public static` so the production body is
tested on hand-built nodes — the duplicate-valued case throws on HEAD (the release fault is a store
through a null container, Grok's correction of the row's "overflow"). DR-197: `~CutFactory`
deletes only the owning temporary sidesets and the restore nulls each slot after its delete; the
unwind test (`duplicate_and_relink_facets()` then destruction without `run()`) is a SEGV and four
valgrind contexts on HEAD, 0 errors with the fix in the `/opt/scls/debug` tree, whole suite and
baseline 0. Gates: `test_fem` 321/321, `test_homology` 41/41, `make check` 17/17 zero warnings,
NI F deck smoke unchanged. DR-195 struck by documentation (usage guide §8.1). DR-113 struck on
the DrClean session's run evidence (its text, folded into this commit's register edit on its
request). DR-193 planned (`todo/dr193_usermat_table_extension.md`), DR-196 waits for Gregory.

Coordination: `main` is pinned at v0.9.2 by the Server session (it mirrors to public GitHub at
02:00), devel only; the Server session's blind xhigh audit of c1d2bf7 raised three points
(cross-class periodic refusal message and test, first-facet master type, append order from an
unordered map) that await Christian's ruling — no code change. The Server session also asked, on
Christian's word, for the shipped decks whose connectivity c1d2bf7 changes to be run against a
f9b9421 reference (rlc_circuit, tapestack2d_racetrack, tapestack3d, tapestack3d_quench, pancake,
helix, racetrack_usermat, pancake_tape); that gate is running in Claude's scratchpad as this is
written and is reported in the next entry. Christian's new close-out rule applies from this
session: say which register rows filed in the session are still open.
