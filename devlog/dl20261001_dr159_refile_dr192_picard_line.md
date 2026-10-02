# DR-159 struck and re-filed as DR-192; the 1e-11 deck makes two builds agree and hits the linear floor; the iteration line names the state — and the 2026-09-29 "lagged Picard residual" mechanism was wrong

**Date:** 2026-10-01 (session 2026-09-30 16:00 to 2026-10-01)
**Purpose:** Christian asked whether DR-159's "PREMISE REFUTED" made the row obsolete. It did not: the
row carried a live residue (marginal exits on a 1e-7 deck, a 12 % two-build divergence, a gate
restriction) and two rulings. He ruled strike-and-re-file and act on both rulings; the deck ruling
produced two 1 ms runs on two binaries of the same source at 1e-11, the print ruling produced a
plan whose premise its own jury refuted, a corrected mechanism, a second plan, and a landed change
**Module:** `src/fem/kernel` (`Controller::iterate_coupled`, `iterate_magnetic`, `iterate_thermal`,
the two print helpers, `cl_FEM_Controller.hpp`), `src/sparse` (`cl_SolverPETSC.cpp`, one comment),
`src/fem/kernel/doc/nonlinear_controller_theory.md`; `todo/debt_register.md`,
`todo/debt_register_closed.md`, `todo/test_gap_backlog.md`, `todo/picard_line_head_residual.md`;
`examples/corc_periodic_bc/input.conf` (`tolerance : 1e-10`, Anderson on)
**AIs involved:** Claude Fable 5.1 (register, runs, plans, implementation, verification, gates); Codex
`gpt-6-astra` (plan v1 high, plan v2 xhigh, code high) + Grok `grok-4.7` (same depths), blind juries
**Status:** code landed, gated in a private tree, code jury clean (no code finding; three doc
sentences corrected); `make check` 17/17 in the private tree and pass in Christian's shared tree, 2026-10-01; the deck ruled and landed 2026-10-02 (`tolerance : 1e-10`, Anderson on), DR-192 struck with its
executable evidence; live register 6, archive 181
**Literature References:** Messe et al. 2023 §4 (always-accept Picard, pre-update residual — the
comment at `cl_FEM_DofMgr_SolverData.cpp:2346-2349` cites it) and §2.7 (the < 1e-11 recommendation
for HTS Picard, the reason the deck was tightened); Ern, Stephansen & Zunino 2009 not touched
**Verification:** executable, on a private HEAD tree in the session scratchpad (`git worktree` at
`0af5c62`, `/opt/scls/mkl`, Blaze, debug, the same 36 options as `cmake-build-debug`, built by
Claude under the 2026-09-29 O5 precedent). DR-192: `corc_periodic_bc` at `tolerance : 1e-11`,
Anderson on, 1 ms, on the shared `cmake-build-debug` binary (pinned copy, `LD_LIBRARY_PATH`) and on
the private rebuild, compared field by field with `scipy.io.netcdf_file` at all ten saved steps.
Print change: syntax gate, private-tree build, three decks (`tapestack2d_layered`;
`tape_quench_usermat` fully coupled and `coupling : segregated`, plugins built in the scratchpad)
read line by line against the certificates; then `make check` in the shared tree, Christian, 2026-10-01: pass, and `make check` in the private tree with the diff, Claude, 2026-10-01: 17/17 (19.95 s) — Christian authorized Claude to run the suite in its own tree that day
**Exchange:** `tmp/ai_exchange/review_picard_line_head_residual.md` (v1 pre-registration, two
refuting entries, verification; v2 pre-registration, two entries, verification; code pre-registration,
two entries, verification; gate evidence) and `.diff` — ephemeral, distilled here and in DR-192

## Summary

"Premise refuted" on DR-159 meant the 2026-09-19 claim was dead — no step of `corc_periodic_bc` is
accepted unconverged — but the row still held the thing that claim had been built on: from 0.55 ms
the head residual sits within 2× of the deck's 1e-7 in a tenth of the steps, two builds of the same
source part by 12 % of the fields there, and any field-level gate on the deck had to stop before
0.6 ms. The row's title stated the dead claim, so it was struck into the closed register with its
refutation evidence intact and the residue re-filed as DR-192 with the gate restriction and the two
rulings; the test-gap backlog's O2 is marked resolved.

**The deck ruling, executed.** With `tolerance : 1e-11` and Anderson on, both binaries converge every
accepted step to 1e-11 (314 steps, max head 9.9e-12, 2.06 Picard iterates on average, no relaxed
exit, no stall exit) and agree to 8.5e-11 on Jz, 2e-12 on phi and 3e-11 on Hx at every saved step
through 1 ms, I_1 identical, U_1 to 4.8e-10, no node above 1e-6 — the divergence is gone, which is
the first executable confirmation that it was the loose tolerance. The cost is new: 40 timestep
rejections, identical on both binaries, each from the stall detector on a residual sitting
bit-identical at −107.7 dB (1.7e-11) for five iterates at Δt = 0.01 ms. That is the linear-solve
floor (STRUMPACK at its default `relative tolerance : 1e-10`) landing just above the target; the
cut halves Δt and the floor drops under it. 354 attempts instead of 108, Δt down to 0.001 ms,
5 h 17 min against ~75 min. Two variants were launched to settle the deck: V1 keeps 1e-11 and sets
the linear tolerance to 1e-12 — nine rejections by 0.56 ms, so the floor is not the linear exit
test alone; V2 uses 1e-10 — 108 steps, no rejection, Δt never cut, 1 h 37 min, max head 9.65e-11,
the shipped deck's shape at the literature's order. V2's fields differ from the 1e-11 run's by up
to 12 % of Jz at 0.9 ms, which measures the two Δt histories and not reproducibility; the 1e-10
deck on the second binary (B2) agrees with V2 to 1e-11 through 0.7 ms and 3e-8 on Jz at 1 ms. Christian ruled 1e-10 with Anderson; the shipped deck
carries it since this session, and DR-192 is struck: the gate restriction is replaced by a 1e-7
relative field tolerance for regression gates on this deck through 1 ms, and the inherent
marginal exit is now four orders below the fields' reproducibility.

**The print ruling, twice.** The first plan proposed printing a "head" residual beside a "lagged"
one on the Picard line, on the mechanism the 2026-09-29 session had written into the register: that
the line prints `r = A_old x_new − b`. Codex and Grok refuted it independently from the source: on
the Picard path `solve_from_residual` backs the residual vector up before the solve and restores it
after (`mRhsBackup`, "Picard reports the PRE-update residual, always-accept semantics, Messe 2023
§4"), so the body line already prints the head residual of the state the solve started from, one
iterate behind the solve count, and the accepted iterate's residual appears only on the certificate
line. Only the Newton branch recomputes a lagged post-update residual — the comment there says
"same lagged-operator semantics as the Picard branch", which is what the 2026-09-29 reading had
copied. The module theory doc had the right statement all along (`nonlinear_controller_theory.md`,
the pre-update paragraph). The refutation of DR-159 stands on the accept test, which reads the
head; the mechanism sentence in the register, the closed row and my memory note were corrected.

Christian chose, from three options, to print the line where the head is measured and index it by
the state. Plan v2 went through a second jury at xhigh (both: print-only, no decision changes;
ten findings on the spec of the resulting log, the doc paragraph and the gate, all folded in), was
approved, and landed: three head prints added under the existing rank-0 / info-level gate, three
body prints removed, `print_line` deleted, the two per-field helpers take the residual as an
argument, the format string untouched. The column now runs 0 … n and its last line equals the
certificate — in the DR-159 case, step 63 reads `… 4, residual 0.004370 ( -23.60 dB )`,
`… 5, residual 0.000000 ( -70.06 dB )`, `timestep succeeded │ magnetic : 9.863e-08`, and the
misread is no longer available. Gated on all three loops, since a magnetic-only deck runs
`iterate_coupled` and the two segregated loops need a thermal deck with `coupling : segregated`
(Codex caught that the first gate missed them).

## DR-159 → DR-192

The closed DR-159 row keeps the 2026-09-29 status text as evidence, prefixed with the strike record
and the correction of its mechanism sentence. DR-192 (P2, `[MIXED]` `[P]`; struck 2026-10-02 once the deck landed) carried: the measured
marginal exits (10 / 8 / 1 of 108 / 108 / 152 steps within 2× of 1e-7 on the three 2026-09-29
variants), the two-build divergence of 2026-09-19 as its consequence, the corrected observability
statement, the 1e-11 result above, the pre-0.6 ms gate restriction until the deck is settled, and
the two rulings' state. `test_gap_backlog.md` O2 → resolved; T-67's stall assertion reads "a stall
exit is a cut, not an acceptance".

## The 1e-11 runs

| run | binary | tolerance | linear rel. tol. | accepted | rejected | Δt min | wall |
|---|---|---|---|---|---|---|---|
| A1 | `cmake-build-debug` at `0af5c62`, pinned | 1e-11 | 1e-10 (default) | 314 | 40 | 0.001 ms | 5 h 17 min |
| B1 | private rebuild of `0af5c62`, same options | 1e-11 | 1e-10 | 314 | 40 (same steps) | same | 5 h 12 min |
| V1 | A | 1e-11 | 1e-12 | 314 | 40 (same as A1) | 0.001 ms | 5 h 10 min |
| V2 | A | 1e-10 | 1e-10 | 108 | 0 | 0.01 ms (never cut) | 1 h 37 min |
| B2 | B | 1e-10 | 1e-10 | 108 | 0 (identical Δt history to V2) | 0.01 ms | 1 h 36 min |

A1 vs B1, worst relative difference per field (scale = max |field| at that step):

| t [ms] | phi | Jz | Jx | Hx |
|---|---|---|---|---|
| 0.1 | 8.3e-13 | 1.6e-11 | 1.2e-11 | 1.3e-11 |
| 0.5 | 1.4e-12 | 2.1e-12 | 2.0e-12 | 2.1e-11 |
| 0.7 | 1.4e-12 | 4.6e-12 | 4.3e-12 | 2.2e-11 |
| 1.0 | 3.2e-13 | 8.5e-11 | 4.8e-11 | 1.5e-11 |

V2 vs B2 at 1e-10, two binaries, worst relative difference per field:

| t [ms] | phi | Jz | Jx | Hx |
|---|---|---|---|---|
| 0.1 | 8.1e-13 | 1.4e-11 | 9.7e-12 | 1.2e-11 |
| 0.7 | 4.1e-13 | 9.7e-12 | 5.5e-12 | 2.2e-12 |
| 0.8 | 4.9e-13 | 6.0e-10 | 5.3e-10 | 1.1e-10 |
| 1.0 | 9.5e-13 | 3.0e-8 | 1.3e-8 | 6.4e-9 |

The two binaries took different Picard paths at 1e-10 (mean 2.52 vs 3.52 iterates per step,
4 vs 36 relaxed exits, identical Δt history) and still met within 1e-8, because each step
converges to the same state at that tolerance; at 1e-7 the same path difference became the 12 %
field split. For comparison, the 2026-09-29 1e-7 pair (shipped vs Anderson, different knobs, same binary) agrees
to 1e-12 through 0.5 ms and parts by 71 % of Jz at 0.6 ms. The rejection mechanism, step 28 of A1:
`Picard 1 … -13.02 dB`, then five lines at exactly `-107.68 dB`, then `timestep rejected │ retrying
step 28 │ Δt = 0.0050 ms`. The deck's `stall tolerance : 0.001` is doing what it is for; the floor
is the linear solve's.

## The Picard line

What the log shows per trip, from the source (`cl_FEM_Controller.cpp` at `0af5c62`): the head is
measured (`compute_residual`, `tEpsilonHead = residual( mIteration )`), the NaN guard, the accept
test — certificate and return — else the body (`solve_from_residual`, Picard restores r, so
`mEpsilon` equals the head), `++mIteration`, print. The line printed after solve k therefore carried
the residual of the state before it, and the state after the last solve was printed only by the
certificate. After the change the print sits between the NaN guard and the accept test in all
three loops and names the state: index = body solves of that field so far, residual = that
state's head, label and `relax` = the algorithm and the initial trial ω the next solve starts with.
Deliberate consequences, all stated in the plan and the theory doc: Newton's post-update lagged
number is no longer printed (the guards still use it); coupled thermal rows are indexed by the
thermal count and appear only on trips that measured a thermal head; the two coupled rows are no
longer adjacent; a head followed by a solver failure or a cut prints its line before the reset
block; a coupled trip whose magnetic head certifies while thermal iterates repeats the magnetic
index; `accepted at stall` prints the last thermal body residual and stays the labeled uncertified
exception.

## Juries

- Plan v1 (astra/high + grok-4.7/high): both refuted the premise from `cl_FEM_DofMgr_SolverData.cpp:2305`
  / `:2350`; Codex also the thermal row on frozen trips, the sentinel width and the −110 dB
  arithmetic; Grok also five mis-placed cites (two inside a NaN branch), the predictor-residual range
  and the `mIteration` / `mIteration2` mismatch of the coupled thermal row. All CONFIRMED; plan
  withdrawn; three options put to Christian.
- Plan v2 (astra/xhigh + grok-4.7/xhigh): print-only confirmed by both; Codex: `relax` is the initial
  trial ω (the coupled Newton line search halves it after the line), the specimen's ω was 0.594 not
  0.5, the gate missed `iterate_magnetic` / `iterate_thermal`, the test rationale was overstated;
  Grok: the label mutations after line-search exhaustion and thermal escalation, the early returns
  that gain a line, the coupled index repeat, the stall exception, the pair split, seven cites, the
  PETSc comment. All CONFIRMED, all folded into the plan; nothing touched the code move.
- Code round (astra/high + grok-4.7/high): no code finding from either; both confirm print-only,
  the four call sites after the NaN return and before the accept test, the gate and the broadcast
  order. Three doc findings on the theory paragraph, all CONFIRMED and corrected in the same
  session: the stall sentence described only the previously-latched path (a newly detected stall
  exits on the trip that printed its head — Codex and Grok independently); "repeats the index on
  the next line" names the thermal row, the repeat is on the next magnetic line (Grok); the
  certificate prints `%10.3e` while the line's fixed field shows `0.000000` below 5e-7, so the dB
  cell is the one to compare (Grok). A NaN head printing no line and `min iterations` lines at
  tolerance were added as clauses.

## Open

- The `git worktree` `treeB` was removed at session end; the pinned binary copy `binA` and the
  six run directories go with the scratchpad.
- The 2026-09-29 devlog's `:103-105` mechanism sentence is wrong as written (records are kept); this
  entry, DR-192 and the closed DR-159 row carry the correction.

## Files Updated

- `src/fem/kernel/cl_FEM_Controller.hpp`, `src/fem/kernel/cl_FEM_Controller.cpp`,
  `src/sparse/cl_SolverPETSC.cpp` (comment), `src/fem/kernel/doc/nonlinear_controller_theory.md`
- `todo/debt_register.md` (DR-192 filed in DR-159's place), `todo/debt_register_closed.md`
  (DR-159), `todo/test_gap_backlog.md` (O2), `todo/picard_line_head_residual.md` (new),
  `todo/README.md`
- `examples/corc_periodic_bc/input.conf`
- `devlog/README.md`, `devlog/campaigns/release_1.0.md`
