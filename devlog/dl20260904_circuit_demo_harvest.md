# Circuit demo harvest: four transient tests in, the demo executable out, DR-39 struck

**Date:** 2026-09-04
**Purpose:** Execute `todo/closed/circuit_demo_harvest_plan.md` (adopted, re-scoped and jury-audited
2026-08-30) end to end: harvest the four circuits hardcoded in `src/executables/electricalCircuit.cpp`
into asserted fixed-Δt transient tests, run them, delete the executable, retarget the documents that
named it, and strike DR-39.
**Module:** `tests/circuit`, `src/executables`, `config`, docs, `todo/`
**AIs involved:** Claude (implementation), Codex (language sweep over the three user-facing rewrites)

## What landed

**Tests — `tests/circuit/test_ElectricalCircuit.cpp`, no new TU, no CMake change (P5 held).**

| case | plan step | gate |
|---|---|---|
| `ResistiveSinePinsCurrentConvention` | R2 | exact, mixed absolute/relative 1e-12; source-branch current asserted **negative**, resistor's positive |
| `FredericCircuitTwoRegimes` | R4 | amplitude 0.998 A over steps 100–499 and 1.506 A over 1201–1600, both from the live component values, ±1 %; switch fire pinned to step 500 via the switch branch current (0 at 499, nonzero at 500) |
| `DiodeBridgeFullWaveRectifies` | R5 | `v(2)−v(3) ≥ −1e-6` every step; exactly two maxima in the second period; peak = iterated fixed point of `v = 10 − 2·Vt·ln(v/(R·Is)+1)` ±1 % |
| `NetlistTwinMatchesFredericCircuit` | R6 | same circuit as a `.cir` through `NgspiceCircuitFactory`; every node voltage and branch current within 1e-12 (absolute floor) of the hand-built run at every one of 1600 steps; node names resolved through `node_index()`, labels compared lowercase |

Helper change: `solve_attempt`/`take_step` take an optional `uint * aNumIterations` so a test can
report Newton health. Existing callers unchanged. `mixed_tol` and a `write_scratch` copy (no HDF5
dependency) added beside them. `create_frederic_circuit()` writes its values in the SPICE parser's
arithmetic form (`200 * 1e-6`, `3.183 * 1e-3`, `24.975 * 1e-3`, `1.0 / 50.0`) so the twin is
bit-for-bit, and puts `t_switch` half a step before the 25 ms boundary so the latch is never a
floating-point tie.

**Deletion (R7).** `src/executables/electricalCircuit.cpp` removed, its `Add_Executable` block in
`src/executables/CMakeLists.txt` and its `BELFEM_INSTALL_EXECUTABLES` entry in `config/globals.cmake`
with it. `Controller::solve_circuit()` is now the only circuit Newton loop in the tree.

**Docs (R8).** The plan's thirteen sites plus two the grep turned up beyond them
(`examples/README.md`, `examples/scripts/Allclean`). The generator `scripts/update_doc_index.py`
was fixed before the two `.dox` files it emits. `CLAUDE.md`, `doc/README.md`,
`src/executables/doc/README.md`, `src/circuit/doc/circuit_usage_guide.md` §11, `Doxyfile.in`,
`examples/scripts/README.md` and `examples/scripts/Allrun` no longer name the binary. Dated records
(`devlog/`, `doc/lessons_learned_evidence.md`, `todo/closed/*`, `todo/debt_register_closed.md`)
kept as written. `scripts/check_doc_claims.py`: 38/38 after the `CLAUDE.md` edit.

**Register (R9).** DR-39 struck and archived to `todo/debt_register_closed.md`; `[P]` count 8 → 7.
(a) discharged by deletion; (b) re-homed: `todo/closed/ngspice_parser_plan.md` Phase 6 now carries
the row's three constraints (not circuit-scoped → Input Contract change in both artifacts; the
7-slot `mValues` schema collides with PULSE and cannot hold PWL at all; PWL storage must be
`Cell< real >`), §7 and Phase 4 of that plan are retargeted so they neither point at the deleted
file nor re-propose the shared helper, and the PULSE-now / PWL-on-demand split is logged in its
§11 as a recommendation, not a decision.

## Evidence

`make check` run by Christian 2026-09-04, 19/19 suites. Each of the four cases appears by name
as `[ OK ]` in `cmake-build-debug/Testing/Temporary/LastTest.log` (17:49 PDT), so they ran rather
than being inferred from the suite total. Runtimes 3 / 20 / 8 / 38 ms — well inside `fast`.

## Pre-registration outcomes (plan §6)

| # | prediction | outcome |
|---|---|---|
| P1 | R5 converges with ω ≡ 1 | **held.** Every step converged; worst case **11** Newton iterations per step. R1′ (relaxation in the helper) therefore **not implemented**, and the 2026-08-30 downgrade of "relaxation is not optional" to "protective" was right |
| P2 | R6 agrees to 1e-12 over all 1600 steps | **held** at 1e-12 with an absolute floor |
| P3 | switch fires at end of step 500 in both runs | **held** |
| P4 | source-branch current negative, resistor's positive | **held** — the MNA convention is now pinned by a test |
| P5 | no new TU, no CMake change | **held** |
| P6 | deleting the file breaks nothing | **held as invocation** — nothing ran it; the residual `electricalCircuit` grep hits are dated records only |
| P7 | any R5 `nan` would come from commutation, not overflow | not exercised — no `nan` occurred |

## Residue

- The `make check` that proves the executable list still links after R7 (`src/executables` and
  `config/globals.cmake` changed) has **not** run on the post-deletion tree: the green run above
  predates the deletion. It is the one gate this devlog cannot claim.
- `examples/scripts/Allclean` still removes `CircuitResults.txt`; nothing under `src/` writes that
  name either. Left alone — not this plan's question.
- `src/circuit/doc/circuit_usage_guide.md` §11's second bullet cites `todo/ngspice_parser_plan.md`,
  which is both a doc-cites-todo violation and a stale path (the file is in `todo/closed/`). Not
  touched: out of scope, flagged for the next doc pass.
- Language sweep over the three rewritten user-facing passages: see the closing note below.

## Language sweep

Codex (`gpt-5.6-luna`, medium; exchange slug `circuit_demo_harvest`) swept the three rewritten
user-facing passages (`src/executables/doc/README.md` overview, `circuit_usage_guide.md` §11 first
bullet, `examples/scripts/README.md` circuit paragraph). Two phrasing improvements applied as
returned. Two factual flags, both accepted and fixed: `db2exo` had been called a "property-table
tool" (pre-existing wording; it converts a property database to Exodus), and "ordinary h-phi runs"
excluded coupled h-ɸ/T runs — now "ordinary `belfem` runs, magnetic-only or coupled as the deck
decides" at all three places that used the phrase (`examples/scripts/README.md`,
`examples/scripts/Allrun`, `examples/README.md`). `git status` after the round: Codex modified
nothing.
