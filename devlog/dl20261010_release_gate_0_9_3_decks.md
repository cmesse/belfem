# Devlog 2026-10-10 — Release gate for 0.9.3: cfaf7e8 against 12a8eec on seven decks, and the final tip 6d5c56c against the f8ce949 and 12a8eec references on all 21 shipped decks

**Date:** 2026-10-10 (runs 2026-10-09 11:36 to 2026-10-10 00:21)
**Topic:** The deck comparisons the Server session requested for 0.9.3 on Christian's word: (A) DrClean's interface-copy ordering cfaf7e8 against 12a8eec on the seven decks of the DR-198 gate; (B) an intermediate iron/coil check of 337a318 (first-facet sidesets, later reverted) against f8ce949, kept for the record; (C) the final tip 6d5c56c (cut-copy ordering 6d5c56c on top of the f3841bc stack 7011f44) against f8ce949 references, 12a8eec for the seven of (A), on every shipped deck, with generator and cut counts recorded per deck. Christian ruled tape_quench_usermat out of the result table
**AIs involved:** Claude Fable 5.1 (runs, tables); Server session (brief, sequencing, pipeline); DrClean (cfaf7e8, 337a318, 6d5c56c); tests session (7011f44, the corc_twolayer generator)
**Claude Confidence:** high — every number is executed; medium on the attribution of the one history that moves (tapestack2d_racetrack), see the tolerance rerun
**Literature References:** N/A

## Method

Serial, release build on the MKL flavor (`-O2`), two OpenMP threads, one solver at a time, every
run in a fresh directory from one regenerated mesh per deck, no `.bfm`, no memdump. Simulation
times cut to reach a first save (1 s costheta and gantry, 20 s undulator2d, 100 ms dipole,
corc_twolayer and the 3-D tape decks, 200 ms the 2-D layered decks, 1 ms the corc decks, 10 to
40 ms the circuit and disk decks, 100 s racetrack_usermat under a 90 min cap, 25 min cap
elsewhere). Final-tip runs under `-v 3`, which prints the 1-cochain count; the references of the
first twelve decks were run again for one step under `-v 3` to obtain theirs. Per deck and side:
base nodes and elements; 1-cochains (generators); the edge count before and after the cuts; the
enriched mesh's node, edge and facet totals and the pinned nodes from the `.bfm` (the facet total
stands in for a cut-face count, which the release binary does not print at any level); prescribed
and free dofs; refusals; the terminal history column by column. corc_twolayer's mesh comes from
the 6cb9685 generator (`numTurns = 2`; the f8ce949 tree's 4/3 winds a layer's tapes into one
conductor and aborts in the cohomology), md5 322ea514, 87,760 nodes. Scripts: `final_gate.sh`,
`gate_table.py` in the session scratchpad; the consolidated table `final_table.txt`.

## (A) cfaf7e8 against 12a8eec, seven decks

Counts identical on all seven; histories to round-off: helix 6.5e-11; rlc_circuit 8.9e-12;
tapestack2d_racetrack 4.1e-11; tapestack3d I_1 identical, U_1 8e-7, the 1e-10-level I_2 is
noise; tapestack3d_quench I_1 identical, U_1 1.1e-6; pancake I_1 identical, U_1 1.9e-8 (fourth
attempt: the OOM killer took three at 25 GB resident alongside other runs; alone it peaks at
38 GB); racetrack_usermat reached its full 100 s on both binaries (24 steps, 1286965 free and
1022 prescribed dofs both), the three tape currents identical, the voltages to 1.3e-8 — the
first racetrack result comparison in any round.

## (B) 337a318 against f8ce949, iron and coil decks (superseded by the revert, kept for the record)

costheta 1.0e-9, dipole 0.0, undulator2d 1.07e-5 (its own Picard-1e-6 noise), gantry currents
identical and 913 voltages at a median of 1.0e-8 with one 1.4 mV column at 3.2e-3, tapestack3d
spot check U_1 3.6e-7; all counts identical, no refusal. The first-facet work left 0.9.3 on
Christian's ruling the same evening.

## (C) 6d5c56c against the references, all 21 decks

Every count identical on every deck: generators, edges before and after the cuts, enriched node,
edge and facet totals, pinned nodes, prescribed and free dofs. No refusal. Histories on the
carrying columns: costheta 1.4e-9; gantry currents identical, 696 voltages above a thousandth
of the largest at median 0 and max 2.5e-5 (the 1.4 mV column U_058 at 1.3e-3, as in every pair
of this day); undulator2d median 0, max 2.8e-5; dipole, disk_pulse, disk_ramp identical;
tapestack3d U_1 3.3e-7; corc_twolayer 5.4e-9; helix I_2 9.6e-10, U_1 8e-12; rlc_circuit 5.1e-9
over 188 rows; tapestack2d_layered 2.5e-7; corc_periodic_bc currents to 1e-13, U_1 1e-11;
corc_garber_effect 1.1e-11; rutherford 3.7e-13; tapestack_circuit 1.7e-9 over 27 rows;
tape_with_contact_and_defect U_1 5e-9; pancake 5.7e-7; racetrack_usermat 8.6e-9 over the full
100 s; tapestack3d_quench U_1 1e-6. The secondary cut quantities of helix, corc_periodic_bc,
tapestack3d, tapestack3d_quench and tape_with_contact_and_defect sit at 1e-10 to 1e-19 and
carry no information. tape_quench_usermat: setup identical, one step, not a gate check.

**The one history beyond solver tolerance:** tapestack2d_racetrack. I_1 identical; U_1 (0.5 V)
agrees to 1e-11 through 1.2 ms, then 1.5e-7 at 1.4 ms, 5e-6 at 1.6, 1.5e-4 at 1.8, 8.6e-3 at
2.0, 7.6e-3 at 2.2 ms. The adaptive Δt sequences of the two runs are identical for 32 steps and
part at step 33 (1.7805 against 1.7789 ms); after that the n = 30 power law amplifies the
round-off of the renumbered system. The same deck agreed to 4.1e-11 in (A), where nothing was
renumbered. Rerun at Newton tolerance 1e-10 instead of the deck's 1e-7, same cut and cap, on both binaries:
U_1 agrees between 12a8eec and 6d5c56c at every row to 1e-10 or better (3.6e-12 at 0.2 ms,
1.7e-10 at 1.4 ms, 9.6e-11 at 2.0 ms), and the Δt sequences no longer part. The divergence at
1e-7 is the controller amplifying round-off that the tolerance leaves in the iterate, not a
changed answer; the 1e-7 answer itself differs from the 1e-10 answer by 0.7 % at 2.0 ms on
either binary, the same pattern the corc deck showed on 2026-09-29.

## What this says

The cut-copy ordering and the f3841bc stack change no count on any shipped deck — the cut sets
and dof sets are what they were — and no result beyond solver tolerance, except where an adaptive
controller on a strongly nonlinear deck amplifies the round-off a renumbered system produces.
Equal counts on this toolchain do not show that the new numbering is platform independent; that
is the ordering commit's own claim and its own test.

## Not done, by design

No register row filed (the Server's rule for this gate). No source change. `main` untouched.
Three detours recorded for the next gate: costheta's geometry has no physical groups, so a lone
`Physical Point` makes gmsh export only that point unless `Mesh.SaveAll = 1` precedes it; the
gate script generates a mesh only when it creates the mesh directory itself; `pkill -f` with a
pattern that matches its own command line kills the shell that runs it.
