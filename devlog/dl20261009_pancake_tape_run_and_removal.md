# The 20-turn pancake_tape case observed on devel and on the commit before DR-198: three φ null pivots on element-less originals where the two cut sheets coincide; the deck removed from the tree

**Date:** 2026-10-09
**Purpose:** Christian and Gregory doubt that DR-196 is real. Its row rests on one debugger
diagnosis of 2026-10-06 that nobody reproduced, and the release gate of 2026-10-09 showed the
symptom change across c1d2bf7. Task 1, relayed by the Server session: run the 20-turn
`pancake_tape` smoke case on the devel tip and report observations for the jury that decides
the row. Task 2: remove `examples/pancake_tape` from devel, since it has never produced a
solution and must not ship in 0.9.3, and fix every reference. Folded in: the wording pass on
the three periodic-refusal passages from the blind audit of 00b7972.
**Module:** `examples` (removal), `src/homology/doc`, `doc/`, `CHANGELOG.md`, `todo/`; the run
touched `src/fem/kernel` and `src/sparse` only through scratchpad probes that never reached
the tree
**AIs:** Claude Fable 5.1 (DrClean session); the Server session (brief, relay to Christian and
to the DR-196 jury record); Codex `gpt-6-astra`/xhigh and Grok `grok-4.7`/xhigh (the Server's
blind audit of 00b7972, whose doc findings are applied here)
**Status:** COMPLETE — findings sent to the Server; deck removed; no register row touched
(DR-196 stays with the jury and Christian); no `src/` change in the commit
**Artifacts:** scratchpad only (gone with the session): `pancake/run20/` (log, `cut_0.vtk`,
`cut_1.vtk`), `run20_pre/` (f9b9421 log), `run20_probe/` (`probe_dofmap.txt`,
`probe_nodes.txt`, `probe_nullpivots.txt`), `probes_scratchpad_only.diff`, the three analysis
scripts. The deck itself is recoverable from git at 00b7972

## Summary

- **Setup.** Debug binary at 00b7972 (SCLS `debug`, gcc, Armadillo, non-MKL), serial, fresh
  directory, no `.bfm`, mesh from `python3 python/main.py 20` (20445 nodes, 117556 tets, 44
  blocks, 50 sidesets, SolderAir = 94). Setup completes: two 1-cochains, 129952 free, 3 fixed,
  72368 condensed dofs. No assertion fires. The first solve soft-fails with MUMPS −10 at
  Δt = 10 ms and at every halving down to 0.078 ms, eight in a row, then the controller's own
  abort (exit 134).
- **Before c1d2bf7.** The same mesh and deck on f9b9421 abort in setup at
  `cl_FEM_DofMgr_DofData.cpp:3112` ("couldn't figure out how do assign dependencies of dof
  1275943 ( edge 561783 )"). c1d2bf7 is the commit that turned the DR-196 assertion into a
  singular first solve.
- **No node without a φ dof.** Scratchpad probes in the dof manager: of 5023 originals that
  carry duplicates, none has zero dofs; the node branch's silent un-hang fired zero times; the
  edge-branch assertion did not fire. Node 6276 (the row's node) has one φ dof, seven
  duplicates and no element of its own. 37 originals have no element and still carry a φ dof;
  33 of them lie on both cut sheets.
- **The sheets.** The run writes `cut_0.vtk` (3246 triangles) and `cut_1.vtk` (655). They share
  10 boundary edges, all on SolderAir: 6276–6281 on the inner lead bar and nine on the sliver
  at z = −1 and −2; 17 boundary nodes; and 66 identical interior triangles in the inner-lead
  and bore region (r 0.8–7.6 mm, z −9.3 to −0.5 mm). Sheet 0 ends on SolderAir (280 edges) and
  the domain boundary (133); sheet 1 ends on SolderAir only. The row's premise, two sheets on
  one conductor-surface edge, holds, and goes further: the sheets coincide over a patch.
- **The singularity.** With MUMPS null-pivot detection (`ICNTL(24) = 1`, probe), every attempt
  reports `INFOG(28) = 3`. The three null pivots are φ node dofs and the first three dofs in
  RCM order (lowest graph degree): node 14913 and node 17427, pure air nodes on both sheets
  where the sheets carry the same five triangles, whose two cut duplicates hold all 24 and all
  20 tets between them (11 + 13, 5 + 15), both condensed on original plus λ, the originals left
  with no element and a free φ dof; and node 5119, a tape-rim node on tape sideset 83 and
  SolderAir, on no sheet, whose one conductor-side copy holds 20 elements and no dof and whose
  original has no element and a free φ dof. Its neighbors 5117 and 5118 have the same shape
  and are not null pivots. With the three pivots replaced, Newton converges to −156 dB on the
  first step and the run proceeds: apart from a three-dimensional null space the system is
  consistent. MUMPS reports −10, not −6, so the rows are not structurally empty.
- **Terminal conditions.** The three fixed dofs are φ-type on nodes 47338 and 23908 (beyond
  the mesh, not duplicates: the two generators' λ nodes, fixed by the single current condition
  in 89 / out 90) and 16291 (an air node at r ≈ 31.7 mm, z = 8.5 mm: the autopin). One φ pin,
  both λ fixed, no current dof left free. Whether one current pair on two generators whose
  sheets overlap is the right posing for a coil whose solder touches itself is Gregory's
  question; the run does not answer it.
- **Reading for the row.** The mechanism DR-196 describes (interface original without φ dof,
  assertion at :3112) reproduces only on the pre-c1d2bf7 binary. On devel the failure is three
  φ null pivots on element-less originals, two of them air nodes where the two cut sheets
  coincide. The "no φ dof" claim is stale; the co-termination premise is confirmed and extends
  to a shared patch. No theory beyond that is offered here; the jury and Christian decide.

## Removal

- `examples/pancake_tape/` removed (13 files: README, input.conf, Allrun, Allclean, the
  generator `python/main.py` and its eight-module package). `examples/README.md`: the `.geo` sentence, the deck count (twenty-one), the table
  row, the "does not yet run" paragraph. `CHANGELOG.md`: the "New example deck" entry deleted;
  the rutherford entry's count corrected to twenty-one, which is what the tree holds without
  the deck (19 in 0.9.2 plus `tape_with_contact_and_defect`, which has no entry of its own, plus
  rutherford). `todo/pancake_inner_leadout.md` parked with the git pointer; its index hook
  updated. Devlogs, `todo/closed/`, the DR-196 row and the code the deck exercised (SolderAir
  typing, closed-rim acceptance) untouched. No `CMake`, test or script referenced the deck.
- The three periodic-refusal passages (CHANGELOG `### Refusals`, the `periodic` section of the
  input reference, guide item 4) now say the refusal applies only when no compatible set holds
  the partner's copy, that the pairing runs at the end of node duplication before the relink,
  that 0.9.2 tied such copies regardless, and that the plane "can be avoided", not "is always
  avoidable" — the Server's reconciled list from the blind audit of 00b7972, no `src/` change.

## DR-196 struck (same day, Christian's ruling)

The jury split: Codex voted keep-as-unverified pending exactly the inventory above (every
required source carries a φ dof on the tip), Grok voted strike. Christian ruled "Strike
DR-196". The row moved to `todo/debt_register_closed.md` with the facts of this run in its
status cell: the assertion was real on the build of 2026-10-06 and reproduces on `f9b9421`;
it is gone since `c1d2bf7`, whose conductor sweep and orphan adoption give every such
original a φ dof; the stale line ranges and the two unsupported sentences ("a single sheet
always leaves a minus side"; the later MUMPS −10 read as this mechanism) are retracted in
the row; the three-null-pivot observation is recorded there and not re-filed; Codex's
residual hazard (the sweep visits conductor blocks only, so an original retained by an
edge-only coil block is not covered) stands in the closing note. Live rows 4, archive 189.

## Open Questions

- For Christian's list, from the Server's audit, not fixed here: `to_string( DomainType )`
  prints both coatings as "Coating", so a LeftCoating|Air against RightCoating|Air refusal
  would show equal names; the reverse directions of each mismatch and several excluded sets at
  once are untested.
- `tape_with_contact_and_defect` ships since 0.9.2 without a CHANGELOG entry.
