# Stoyan's enclosed-slit tape: the closed-rim abort, a devel run, and the deck error behind "isolation not shown"

**Date:** 2026-10-09
**Purpose:** Stoyan (FNAL, binary 0.9.2 install) reported that a REBCO tape with a fully enclosed
longitudinal slit aborts BELFEM in `CurveFactory::thin_shell_side_curves` with
`Mat::operator(): index out of bounds`, while a slit open to one end face runs. Christian asked
whether the problem was already fixed, then to check the deck in `cmake-build-debug/stoyan`, then to
run the enclosed geometry on devel and to mail the findings.
**Module:** `src/mesh` (`cl_CurveFactory`), `src/homology` (`cl_CutFactory`), `src/fem/maxwell`
(`cl_MaxwellFactory` terminal association); no source change
**AIs involved:** Claude Fable 5.1 (investigation, run, analysis); Codex `gpt-5.6-luna`/medium
(language sweep of the outgoing mail only, no audit); the Server session (mail relay)
**Status:** COMPLETE — mail sent; no code owed for this case
**Verification:** executed — devel 337a318, release build on `/opt/scls/mkl`, scratch worktree;
`belfem` on the coarsened enclosed-slit deck: setup passes, 18 BDF2 steps to 100 µs, exit 0;
post-processing of the Exodus output and `iv_results.csv` in Python
**Exchange:** `tmp/ai_exchange/stoyan_slit_mail_sweep.md` — the Codex sweep thread, nothing else

## Summary

The abort is the closed-rim `identify_subchains` defect recorded in
`devlog/dl20261003_closed_side_curve.md`, fixed on devel in 9f569c7 (2026-10-08), in no release:
an enclosed slit makes the shell's free edge a closed loop, and 0.9.2 indexes out of range on the
edge that closes it. On devel the slit rim is accepted as an isolated closed side curve, the
cohomology finds the hole generator on its own, and the run completes. The colleague's
"isolation not shown" result comes from his deck, not from the slit: the two unfed curves of the
second ribbon are free terminals, which shorts that ribbon through the domain boundary.

## What the files in `cmake-build-debug/stoyan` are

`slit.geo`/`slit.msh`/`solver.log` are the **open-ended** slit, run on 0.9.2 (built 2026-10-05):
6614 nodes, 18 steps to 100 µs at 2 A, J/Jc 0.044, success. `single_tape_100mm.geo` is the
**enclosed** slit (0.1 mm gap, 35–65 mm) without a mesh; its default refinement field along the
slit (`Field[2].VIn = 0.025`) gives 178 747 nodes and 1.19 M tets, so the "coarse diagnostic" the
colleague's agent describes was not coarse. The deck feeds `input curves : [1]`,
`output curves : [3]` and defines curves 2 and 4 without any condition.

## What the run showed (devel 337a318, slit refinement 0.2 mm, 21 322 nodes, 146 380 tets)

- Mesh topology from the `.msh` (scratch script): the shell's one-facet edges form the outer rim
  (split into two open chains once `select_edges` drops the end-face edges) and one closed loop of
  304 edges on the slit rim that touches no end-face node — the "isolated" case of
  `CutFactory::refuse_closed_side_curve_sharing_nodes` (`cl_CutFactory.cpp:2013`).
- With the condition corrected to `[1,2]` / `[3,4]` the run passes `create_cuts` and completes.
  `slit.bfm` holds seven curves (four terminal, three side) and two cut sidesets: `cut_1`
  (4314 faces, the terminal sheet from the tape's left edge to the domain boundary over the full
  length) and `cut_2` (1783 faces, under the slit, x ∈ [−1.7, 1.0] mm, z ∈ [31, 66] mm).
- `iv_results.csv`: `I_1 = 2 A` (fixed, the ramp), `I_2 = −1.0 A` at `U_2 ≈ 2e-16` — the second
  abstract dof is the free hole generator, carrying the current linked by a loop through the slit
  around one ribbon. Nodal `Jz` on the wedge block gives 1.0 A per ribbon before, inside and after
  the slit; the hole-interior air nodes are single-valued (149 nodes, φ ≈ −1.0).
- The earlier claim in this session that "the cohomology runs before thin-shell duplication, so
  no cut spans the hole" was **wrong**: `MaxwellFactory::create_cuts_sub_master` calls
  `tFactory.create_thin_shell_cuts()` before `tFactory.run()` (`cl_MaxwellFactory.cpp:1181`), so
  the shells are already opened when `compute_cohomologies` runs and a hole in a shell is an H1
  generator of that complex. Recorded in memory so it is not re-derived.

## The deck error

Every curve whose sideset A or B is a shell sideset is pushed onto that shell's
`terminal_curves()` (`cl_MaxwellFactory.cpp:~3620-3660`), and `close_terminal_loops` duplicates
all of them; only the curves named in a current condition reach `mTerminals` and get a fixed
current. A terminal curve outside every condition is therefore a free cut dof, i.e. the ribbon is
shorted through the boundary. That explains the ~0.6 current density in the "unfed" ribbon of the
open-slit run. Remedy sent to Stoyan: one bracket group per crossing for the enclosed slit
(`[1,2]` in, `[3,4]` out, the rutherford convention); a separate 0 A condition on `[2]`/`[4]` for a
deliberately unfed ribbon (not run here).

## Not verified / open

- The 1 A / 1 A split follows from symmetry; an asymmetric case (unequal ribbon widths or a
  defect in one ribbon) is the real test that the free hole dof redistributes current correctly.
- The inner slit edges show no current peak at 0.2 mm resolution for a 0.1 mm gap; the mesh
  cannot resolve it, so this says nothing.
- The colleague's "merged sideset, one ribbon at zero Jz" result may be the mixed-winding problem
  that Yufan's alignment pass in 9f569c7 repairs; not confirmable without his deck (low
  confidence).
- `slit.bfm` records `/curves/closed` as 1 for all seven curves, including the open outer side
  chains; not pursued (the terminal refusal at association passed, so `is_closed()` was false
  then). Possibly a writer defect or a post-`close_terminal_loops` state; worth one look if the
  bfm curve flags are ever read back.
- Server's note after sending: 337a318 also carries two interface-processor commits (e30bd98,
  411e69a) that are being reverted before 0.9.3; they do not touch this case. Stoyan's 200 mm deck
  is to be run on the released 0.9.3.

## Mail

Sent 2026-10-09 by the Server session on Christian's yes: from belfem@lbl.gov to stoyan@fnal.gov,
reply-to cmesse@lbl.gov, signed Claude, subject "Your slit-tape deck: the BELFEM abort, and one
change your input deck needs"; relay accepted (250 OK). Text swept by Codex luna/medium, no
suspected errors flagged. The final text is not in the repository (scratchpad `mail/final.md`).

## Files Updated

| File | Change |
|---|---|
| `devlog/dl20261009_stoyan_enclosed_slit_tape.md`, `devlog/README.md` | this record |
