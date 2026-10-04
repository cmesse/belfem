# DR-105 R1 lands: a floor retry answers the reason it was cut — a stall-guard cut at the Δt floor suspends the guard's Picard verdict instead of replaying the attempt

**Date:** 2026-10-03 (afternoon and evening; the same day's campaign devlog is
`dl20261003_dr105_weekend_campaign.md`)
**Purpose:** Christian's "Go" on plan v2. `reset_timestep` now takes the reason for the cut at all
twenty call sites; when a cut demanded by the magnetic stagnation guard finds Δt already at the
minimum, the guard's Picard verdict is suspended until the next accepted step, so the floor retry's
raised iteration budget and watchdog window can bound the attempt instead of the guard ending it at
the same iterate twenty times. Print-, policy- and signature-level change in the controller; no
acceptance criterion moved.
**Module:** `src/fem/kernel` (`cl_FEM_Controller.hpp`: `ResetReason`, `mStallSuspended`, the
signature; `cl_FEM_Controller.cpp`: twenty call sites, the floor block, the guard's Picard arm,
`finalize`); `src/fem/doc/timestepping_strategy.md` §4, `src/fem/kernel/doc/nonlinear_controller_theory.md`
(guard row), `doc/input_file_reference.md` (floor policy prose; no key contract changed, schema
untouched); `todo/dr105_picard_stall_plan.md`, `todo/debt_register.md`
**AIs involved:** Claude Fable 5.1 (implementation, gates, verification); Codex `gpt-6-astra`/xhigh +
Grok `grok-4.7`/xhigh (blind code jury); Codex `gpt-5.6-luna`/medium (language sweep of the three
doc passages)
**Status:** code landed and committed; **reviewed, and gated by execution only in part** — the
mechanism's own gate (an R1 Picard run meeting the floor and recovering) had not been observed when
this was written; see §Gates for what ran and what is still running
**Literature References:** Messe et al. 2023 §4, Eqs. 12–14 (the lagged Picard update the guard
judges; the Eq. 13 hand-over the escalation keeps); the controller theory note's "always-accept
Picard" reading
**Verification:** syntax gate on both TUs; private HEAD tree (`5da57f7` worktree + the diff, MKL,
Blaze, debug): `make belfem`, then `make check` 17/17 with the jury diff and again 17/17 with the two
post-round fixes; R0a (guard switched off by `stall tolerance : 0`, serial, 25 ms): the plateau the
guard kills resolves by itself in every quench step (flat runs up to 13 iterates, 28 steps with a
run ≥ 8), the 32 rejections all by the watchdog at the Δt cap, no floor retry — the measurement R1
is sized by; the serial R1 run to 20 ms: 0 floor retries, but on a different trajectory than the
HEAD abort (see §What the gate taught), so it is not the reproducer; after the commit, `make check`
on the committed HEAD `b01911f` in the private tree (controller sources byte-identical to the
commit), at Christian's request because he was away from the office: 17/17 in 36.75 s
**Exchange:** `tmp/ai_exchange/review_dr105_picard_stall.md` (code pre-registration, two auditor
entries, verification + reconciliation, the post-round fixes note) and `.diff`

## Summary

The diff is 23 header lines and about 70 source lines. A private `enum class ResetReason { Stall,
Watchdog, Budget, LineSearch, SolverFailure, NaN, Circuit }` is the argument of `reset_timestep`;
every one of the twenty call sites names its reason, so a site that forgets does not compile. In
the at-floor branch of `reset_timestep`, after the budgets and watchdog windows are multiplied by
the retry factor, a `Stall` reason sets `mStallSuspended` — identical member math on every rank,
before the rank-0 print, which gains a line naming what was done. In
`magnetic_stagnation_forces_reset` the Picard arm, while suspended, still offers the one Newton
escalation per attempt and then withholds the cut; the Newton→Picard demotion arm is untouched.
`finalize` clears the flag in the block that already restores the deck budgets, so "flag set"
implies "floor regime", and the segregated driver, which never calls `finalize` after a rejection,
needs no extra reset.

The code jury found no defect: both auditors confirmed the lifecycle, the MPI placement, the
preserved fallbacks and the hard bound (the raised iteration ceiling). What they corrected was
prose and two labels: the suspension starts with the first floor retry, one attempt after the
cut that pins Δt (the clamp branch does not arm it) — the doc had said "for the rest of the floor
regime" as if from the pinning cut; the watchdog is a conditional bound, not a hard one; "the
guard returns at the sixth iterate" is the DR-105 trace, not an invariant (window + 1); the box
line blamed the stall guard on later non-stall floor cuts; the segregated thermal watchdog site
was labelled `Budget`. The last two were applied with Christian's approval after the round, the
rest as exempt doc and comment edits, and the input reference gained the cost note for
`floor retries : 0` (a suspended attempt spends the raised ceiling on every retry).

## What the gate taught

The planned reproducer was the serial `p2_np1` deck, which aborted at 15.41 ms on HEAD and is
step-identical across two cold runs of the same binary (567 headers compared). On the R1 binary
the same deck, with `simulation time` shortened to 20 ms for the gate, passed 16 ms with no floor
retry — but on a different trajectory from step 197 on, with the thermal residual of the first
iterate already different at the fourth printed digit. Two controls separated the causes: the HEAD
binary on the shortened deck differs from the HEAD 0.3 s run from the first step's iterate lines,
and the HEAD and R1 binaries on the same shortened deck differ from the first step's second
iterate. So two runs of one binary and deck replay exactly, while any change — a binary rebuilt
from the same sources plus this diff, or the `simulation time` key, which enters nothing but the
driver's loop bound — moves the residuals at the 1e-2 relative level from the first solve. The
mechanism of that sensitivity was not pursued (memory-layout-dependent summation order in the
MKL-backed solves is the candidate; recorded as an observation, not a finding). Consequence: a
cross-binary replay of the HEAD abort is not available, and R1's own gate is an R1 Picard run
that meets the floor on its own trajectory: the warm start from the aborted run's 15 ms save
(`g4_restart_r1`) and the two full-length 0.3 s runs, serial and eight ranks (`g2_serial_r1_full`,
`g3_np8_r1_full`), launched at close. Every pure-Picard run of this deck so far met the floor
between 10.6 and 33 ms.

## R0 diagnostics

- **R0a, guard off, serial, 25 ms** (above): the plateau resolves given patience; the watchdog
  is the bound; Δt never reached the floor. This is the direct support for the suspension.
- **R0c, np = 4 cold on HEAD, 13 ms:** clean — the 11.745 ms abort of the campaign's np = 4 run
  did not reproduce, as expected for a parallel trajectory.
- **R0b, Anderson off, serial, 20 ms:** clean, 1173 accepted, 66 rejected, no floor retry. In the
  quench plain Picard needs more iterates (median 27 against 19 with Anderson 3) and shows the same
  first-body kick (+23.9 against +23.6 dB), but almost no flat plateau: longest flat run median 2,
  max 5, two steps with a run ≥ 5, against 4, 10 and nine steps with Anderson. The map makes the
  kick and a slow descent; Anderson's mixing makes the flat cluster the five-sample guard trips on.

## Open

- The mechanism gate: a suspended floor retry observed to recover (or to end on the watchdog or
  budget) on the R1 binary; R3a's Picard-vs-hybrid temperature comparison from the same runs.
- R2 (the T5 deck on the hybrid) waits for the Picard gate runs, which need the Picard deck.
- O4 (`target iterations : 100`), the table extension (O3's separate task).
- The run-to-run sensitivity of this deck to any binary or deck change (observed, not explained).

## Files Updated

- `src/fem/kernel/cl_FEM_Controller.hpp`, `src/fem/kernel/cl_FEM_Controller.cpp`
- `src/fem/doc/timestepping_strategy.md`, `src/fem/kernel/doc/nonlinear_controller_theory.md`,
  `doc/input_file_reference.md`
- `todo/dr105_picard_stall_plan.md`, `todo/README.md`, `todo/debt_register.md`
- `devlog/README.md`, `devlog/campaigns/release_1.0.md`
