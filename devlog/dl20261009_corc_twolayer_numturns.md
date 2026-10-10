# corc_twolayer sets up again: its generator winds 2 turns per period, not 4/3; tape_with_contact_and_defect gets its CHANGELOG entry

**Date:** 2026-10-09
**Purpose:** Two tasks from Christian, relayed by the Server session that coordinates 0.9.3.
Task 1: `examples/corc_twolayer` stops in setup on devel at `Cohomology::clean_spfa`; find the
cause by bisecting generator, code, and deck, and make it set up and run again without touching
the closed cohomology core. Task 2: give `examples/tape_with_contact_and_defect`, new since
0.9.2, a CHANGELOG entry that states what the deck does today.
**Module:** `examples/corc_twolayer` (generator parameter, README), `examples/README.md`,
`examples/rutherford/README.md`, `CHANGELOG.md`
**AIs:** Claude Opus 5.5 ("tests" session); Codex `gpt-5.6-luna`/medium (language sweep of the
new prose, slug `corc_twolayer_numturns`)
**Status:** COMPLETE. `numTurns = 2` ruled by Christian (AskUserQuestion). No `src/` change;
no register row touched.
**Artifacts:** scratchpad only (gone with the session): release trees of 0d4d1b5 and v0.9.2
(mkl flavor, serial), run directories `tip_shipped`, `v092_shipped`, `c1_noper`, `c2_solder`,
`tip_res2`, `tip_augmesh`, `tip_turns2`, `twcd`

## Task 1: the cause

The shipped deck on the mesh `python/main.py` writes today aborts after 42 s on 0d4d1b5:
"generator 0 is a thick cut", with a 9-edge loop. The three axes:

| axis | test | result |
|---|---|---|
| generator | today's `corc.msh` vs the one in Christian's 2026-09-02 runs (`belfem_transfer/cmake-build-debug/corc_solder*`) | bit-identical (md5 `e8d24b2e…`); deck files unchanged since v0.9.2 |
| code | v0.9.2 (ef6d47f) build, same mesh, same deck | identical abort, same loop, same node ids: no regression, nothing to bisect |
| deck | `periodic` block removed, nothing else | sets up and steps (t = 800 ms) |

The loop sits nowhere near the conductor: r = 52–56 mm (domain radius 60), fixed φ ≈ 60°,
running axially through the periodic seam twice (z 0 → 25.13 mm, folded, 0 → 25.13 mm again),
and the periodic map is a pure translation. Refining the outer air (`domainResolution` 4 → 2 mm,
149,635 nodes) still aborts, now on a 9-edge loop over a single period, so the class integrates
to at least 10 along one axial period. `domainResolution = 3` failed in the generator itself
(`CapError: no cap node within 1e-6 of dom slot 1`).

Christian answered the first report with "corc_twolayer worked with periodic before". It had,
on a mesh that predates the 0.9.0 history: his runs `corctest` and `dr79` (2026-08-26) and
`g_sign_corc` (2026-08-28) under `belfem_transfer/cmake-build-claude/` use the shipped topology,
periodic block included, on a mesh 37.70 mm long; today's is 25.13 mm. With `pitch = 3` that is
`numTurns = 2` against `4/3`. Christian's standalone generator (`~/codes/corc/main.py`) still
carries `numTurns = 2`; the in-tree `main.py` has carried `4.0 / 3.0` since v0.9.0.

- The August mesh (md5 `1361a088…`) with the shipped deck on 0d4d1b5: 21 steps to t = 1 s,
  no error.
- Today's generator with only `numTurns = 2`: 87,760 nodes, 556,398 elements, 3551 cap pairs;
  the shipped deck sets up and steps. First save at t = 100 ms on 0d4d1b5 (step 6, about
  52 s of wall time per step, serial, on a shared machine): all six tape currents are 2.2 A,
  which is the ramp's 200 A × 0.11 s / 10 s; tape voltages are 3.0249e-6 V (one layer) and
  3.0058e-6 V (the other); the seventh generator, the free axial loop of the periodic model,
  carries 2.195 A at 6.6e-13 V.

With 4/3 turns the twist carries tape 4 → 5 → 6 → 4, so each layer is one conductor that
runs three times through the seam. The deck addresses six conductors, one bracket per tape, so
even past the cut stage it would meet the count refusal that `rutherford` met in September
(dl20260918_rutherford_cable_generator.md, "5 conditions, 2 cohomology generators"). That same
devlog records the outer-air thick-cut refusal for rutherford's chained tapes, which refinement
did not cure either. With 2 turns every tape closes on itself, and topology and deck agree.

Christian ruled: restore `numTurns = 2`. The `main.py` comment now says why the number must be
whole for this deck. The committed generator reproduces the tested mesh (md5 `322ea514…`).

Side observation, not pursued: the topology the generator prints (`corc_topology.conf`: block 2
a solder conductor, terminals 13/14), dropped into the deck, aborts at
`cl_Cohomology.cpp:1088` (`aNumConditions <= n`), a conditions-versus-generators count check.
The shipped deck does not use that topology.

For Gregory, should the chained case matter: the 4/3 mesh with the shipped deck is a 42-s
reproducer of the thick-cut refusal (loop nodes 19286→51210, 41524→51210, 41524→48995,
48995→53116, 23295→53116, 19979→40682, 40682→42742, 40349→42742, 22602→40349; 23295≡19979 and
22602≡19286 are periodic pairs).

## Task 2: tape_with_contact_and_defect

No session had a run on record (Yufan's never ran it), so it was run on 0d4d1b5, serial, with
the plugin built against that tree and `tape.msh` from `gmsh -3 tape.geo`. It sets up and steps
through the sigmoid at Δt = 10 ms, with coupled Newton converging, up to about 100 A. Past the
defect's Ic (about 110 A at t ≈ 0.28 s), two steps are rejected (steps 26 and 29). After that
none are rejected, but the controller shrinks every step: 1.2 ms at 262 ms, 0.1 ms at 282 ms,
about 45 s of wall time per step. The entry says so.

CHANGELOG: rutherford entered the tree first (2026-09-23), so its sentence now reads "Twenty
decks now ship". The new entry closes with twenty-one: 19 in 0.9.2, plus rutherford and
tape_with_contact_and_defect, with `pancake_tape` removed in f6618b9.

## Also changed

- `examples/README.md`: corc_twolayer's periodic column changes from "twisted" to "translated".
  The cap map is a pure translation and was one under 4/3 as well.
- `examples/rutherford/README.md`: the aside "as in `corc_twolayer`" on permuting closures is
  removed.
- `examples/corc_twolayer/README.md`: the deck declares the interlayer bed air, not solder;
  the periodic sentence describes the 2-turn closure.
