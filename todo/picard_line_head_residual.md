# The iteration line names the state: the residual is printed where it is measured

**Date:** 2026-09-30
**Purpose:** The per-iterate log line `Magnetic Picard k, residual x ( dB ), relax ω` is printed after
body solve k and carries the pre-update residual of the state that solve started from, so the
residual column runs one iterate behind the solve count and the accepted iterate's residual appears
only on the `timestep succeeded` line. Read naively, a −25 dB last line before a converged accept
looks like an accepted stall (the 2026-09-19 reading that became DR-159, premise refuted 2026-09-29;
residue DR-192). The change moves the print to the head, where the residual is measured, and
indexes the line by the state: the column then runs 0 … n and ends with the accepted iterate.
Print-only: no controller decision changes. Chosen by Christian 2026-09-30 (option B of the three
put to him after the first plan was withdrawn by its jury, §10).
**Module:** `src/fem/kernel` (`cl_FEM_Controller.hpp`, `cl_FEM_Controller.cpp`)
**AIs involved:** Claude (plan + implementation), Codex astra (audit), Grok 4.7 (third voice)
**Status:** ✅ COMPLETE (2026-10-01) — landed: the iteration line of each field is printed at its head and names the state (three head prints added, three body prints removed, `print_line` deleted, the two helpers take the residual); the theory doc paragraph written and corrected on the code jury's three doc findings. Verified by execution in a private HEAD tree (`treeB`, MKL/Blaze debug): syntax gate, build, `tapestack2d_layered` 97 steps, `tape_quench_usermat` fully coupled 4 steps and segregated 3 steps — every block's last line equals its certificate. Juries: plan v1 withdrawn (astra/high + grok-4.7/high), plan v2 clean on the code (astra/xhigh + grok-4.7/xhigh), code round no code finding (astra/high + grok-4.7/high). Owed: Christian's `make check-fast` on the shared trees. Residual follow-ups: none on the print; DR-192 item 1 (the deck) continues on the register

> **Scope guards:**
> - No change to any controller decision: `mEpsilon` / `mEpsilon2`, `mResidualHistory`, the
>   stagnation guard, the watchdog, the divergence strikes, the line search, Anderson and the ω
>   growth keep reading exactly what they read today.
> - The line format is unchanged: `   │ Magnetic Picard %5u, residual %s (%7.2f dB ), relax %7.5f │`
>   with `residual_string()`; only its position in the loop and the meaning of the index change.
> - The certificate lines (`timestep succeeded`, `magnetic converged`, `substep succeeded`,
>   `accepted at stall`) are untouched.
> - The circuit line (`cl_FEM_Controller.cpp:2540`) is out of scope.

---

## 1. What the log shows today (source trace, confidence high)

Per trip of `iterate_magnetic`: assemble at the current state x_k, `compute_residual()`
(`:2116`), `tEpsilonHead = residual( mIteration )` (`:2118`); NaN guard (`:2121-2125`); accept predicate
(`:2127`), return branch `:2130` → certificate line (`:2140-2146`) and return. Otherwise the body:
`solve_from_residual()` (`:2166`), which on Picard backs r up and restores it after the update
(`cl_FEM_DofMgr_SolverData.cpp:2305`, `:2350`), so `mEpsilon = residual( mIteration )` at `:2178`
is the same number as `tEpsilonHead`; `++mIteration` (`:2191`); `print_line_magnetic( tOmega )`
(`:2193-2196`) prints `mIteration` and `mEpsilon`. So line k carries the residual of x_{k−1},
and x_n (the accepted state) is printed only by the certificate. On Newton the body recomputes a
post-update residual with the lagged operator (`:2291-2298` in the solver data), so there line k
carries a number that is neither x_{k−1}'s nor x_k's head.

`iterate_coupled`: magnetic head `:1098-1101`, NaN guard `:1107-1111`, accept `:1113`; thermal
head `:1601-1604` inside `if( tUpdateThermal )` (`:1489-1494`), NaN guard `:1606-1610`; body stores
`:1225-1226` and `:1644-1645`; `++mIteration` at `:1774-1777` (gated on `tSolveM`), `++mIteration2` at `:1696`; one `print_line(
tOmega, tOmega2 )` at `:1785` prints both rows after the bodies — the thermal row with
`mIteration`, not `mIteration2` (`:3598`), and whenever `mKernel2 != nullptr` (`:3591`), including
trips on which no thermal head was measured.

`iterate_thermal`: head `:2358-2361`, NaN guard `:2363-2374` (two reset flavors), accept `:2377`,
certificate `:2385-2393`; body store `:2420`; `++mIteration2` `:2432`; `print_line_thermal( tOmega )`
`:2434-2437` with `mIteration2`.

## 2. Mechanism

- **R1** `print_line_magnetic( const real aEpsilon, const real aOmega )` and
  `print_line_thermal( const real aEpsilon, const real aOmega )` (`.hpp:802-805`, `.cpp:3606-3633`):
  the residual becomes an argument; the body of each helper is otherwise unchanged (same format,
  `residual_string( aEpsilon )`, `10·log10( aEpsilon )`, `mIteration` / `mIteration2`).
- **R2** `iterate_magnetic`: call `print_line_magnetic( tEpsilonHead, tOmega )` under the
  `mCommRank == 0 && gLog.info_level() > 0` gate right after the NaN guard (`:2125`) and before the
  accept test (`:2127`); delete the body call (`:2193-2196`). `tOmega` is bound at `:2100` and clamped at `:2104`
  before the head.
- **R3** `iterate_coupled`: `print_line_magnetic( tEpsilonHead, tOmega )` after `:1111`, before
  `:1113` (`tOmega` is bound at `:1058-1059`, handed over at `:1087`); `print_line_thermal( tEpsilonTHead, tOmega2 )` after `:1610`,
  before `:1612` (`tOmega2` is copied at `:1588`); delete the `print_line( tOmega, tOmega2 )` call
  and its comment (`:1780-1786`). Delete `Controller::print_line` (`.hpp:798-799`, `.cpp:3577-3601`)
  — it has no other caller (`grep` this session: one call site).
- **R4** `iterate_thermal`: `print_line_thermal( tEpsilonTHead, tOmega )` after the closing brace of the NaN branch at `:2375`, before
  `:2377` (`tOmega` is bound at `:2349-2352`); delete the body call (`:2434-2437`).
- **R5** `src/fem/kernel/doc/nonlinear_controller_theory.md`, a "reading the iteration line"
  paragraph after the pre-update paragraph (`:30-33`) and consistent with the Newton sentences at
  `:41-46` (Newton's stored ε is the updated iterate and leads Picard at the switch — still true
  of the controller's `mEpsilon`, no longer of the printed column): the line names the state,
  k = body solves of that field so far (magnetic and thermal count separately in coupled mode),
  its residual is the head of that state under the fresh assembly, the label and `relax` are the
  algorithm and initial trial ω the next solve starts with (the coupled Newton line search may
  halve ω after the line), the column runs 0 … n on the success path and the certificate repeats
  the last value, a coupled magnetic-skip trip repeats an index, and `accepted at stall` is the
  uncertified exception that prints the last thermal body's residual. Doc text, no jury.
- **R5b** `src/sparse/cl_SolverPETSC.cpp:270`: the comment names `Controller::print_line()`;
  point it at `print_line_magnetic()` in the same diff.
- **R6** devlog, `todo/README.md`, DR-192 status cell (item 2).

## 3. What the log shows after (the DR-159 case, step 63 of the 2026-09-29 run)

```
   │ Magnetic Picard     0, residual 0.018005 ( -17.45 dB ), relax 1.00000 │
   │ Magnetic Picard     1, residual 0.000004 ( -54.45 dB ), relax 1.00000 │
   │ Magnetic Picard     2, residual 0.000002 ( -57.65 dB ), relax 1.00000 │
   │ Magnetic Picard     3, residual 0.006855 ( -21.64 dB ), relax 1.00000 │
   │ Magnetic Picard     4, residual 0.004370 ( -23.60 dB ), relax 0.50000 │
   │ Magnetic Picard     5, residual 0.000000 ( -70.06 dB ), relax 0.59427 │
   ├────────────────────┬─────────────────────────┬────────────────────────┤
   │ timestep succeeded │  magnetic :  9.863e-08  │                        │
```
Same residuals as today, one line more, index shifted by one: the drop from −23.6 dB to −70.1 dB
is in the column. (The block is the 2026-09-29 step-63 log, scratchpad `dr159/shipped/run.log`,
with the certificate's value added as the state-5 line.) The `relax` on line k is the **initial
trial ω** the next solve starts with: the value the growth rule left after solve k (line 5 shows
0.594, the rule's `β + γ·(2/π)·atan((ε0−ε)/ε0)` with β = 1.1, γ = 0.4 (`cl_FEM_Controller.hpp:75-76`) applied to 0.5 — found by
Codex), on a promotion trip the first-flip damped store (`:2070-2071`, `:1016-1017`), and on a
certification-only trip the ω of a solve that does not run. In the coupled Newton line search the
trial ω is halved after the print (`:1449`), so the ω the accepted trial applied is on no line.
Today's line k shows the ω after that search, so the two columns are not a plain one-line shift
on that path (Grok).

## 4. Deliberate changes of log content (not defects)

- **Newton** loses its post-update lagged number from the log; the next head (the true residual
  of that iterate) is printed instead. The line search and the guards still consume the lagged
  value; nothing but the print changes.
- **Coupled thermal rows** are indexed by `mIteration2` (thermal solves) instead of `mIteration`
  (Grok, 2026-09-30, on the v1 round: `:3598` vs `:3629`), and a thermal row appears only on a
  trip that measured a thermal head (`tUpdateThermal`), so a frozen trip prints no stale thermal
  row.
- **A NaN head prints no line**; the reset's own `timestep rejected` block follows. Today a NaN
  body residual prints `nan` on the line first.
- **`min iterations`**: a converged head still prints its line, then the forced body runs and
  the next head prints again — the column shows the forced iterate.
- The algorithm label on line k is the algorithm the next solve starts with (promotion is
  decided at `:2063-2064` / `:995-1058` before the head). Today's label is the live algorithm
  at body-print time, which on two coupled paths is not the one that ran: after a line-search
  exhaustion (`:1393-1402`, message `:1412-1416`) today's line says Picard for Newton trials and
  prints the Newton ω; after a thermal flat escalation (`:1673`) today's thermal row says Newton
  for a Picard body. The head print is the consistent one (Grok).
- **Early returns gain a line, and a `nan` body line is lost.** A measured head followed by a
  `solve_failed` return (`:2169-2175`, `:1216-1222`, `:1635-1641`, `:2403-2417`), a line-search
  cut (`:1430-1434`), a thermal watchdog cut (`:1701-1704`) or a thermal NaN (`:1606-1610`) now
  prints its state line before the reset block; today those paths return without a line. A
  finite head whose body goes NaN prints the finite head; today the `nan` body value is printed
  (`:2193` before `:2207`). The "column ends with the accepted iterate" statement is the success
  path (Grok).
- **Coupled magnetic-skip trips repeat the index.** When the magnetic head certifies but thermal
  still runs (`tSolveM` false, `:1127-1128`), `mIteration` does not advance (`:1774`), so the
  next magnetic head line carries the same index with a different residual (thermal moved the
  state). Today's print shows the same index too; nothing shifts on that trip (Codex, Grok).
- **The coupled rows are no longer a pair.** The magnetic line prints after `:1111`, the thermal
  line after `:1610`, with the magnetic body (and any line-search message) between them (Grok).
- **`accepted at stall` is the one certificate that does not repeat a printed head.** It prints
  `mEpsilon2` of the last thermal body (`:1644`, `:1740-1743`) and the accepting trip measures
  no thermal head (`:1489-1492`); on Newton that stored value is the lagged recompute, so it
  differs from the last thermal line. Uncertified by design, labeled so (Codex, Grok).

## 5. What could break

- Nothing functional (R2-R4 move a print and delete a helper; R1 changes a private signature
  with three callers). Confidence high.
- Width: unchanged format, unchanged widths; `%5u` and `%7.2f` are what shipped.
- Consumers of the line: no tracked parser of these iteration lines (`git grep` over
  `scripts/`, `python/`, `tests/`, `doc/`; `tests/mesh/test_BfmFile.cpp:218` captures stdout for a
  mesh warning, not for these lines); the ephemeral `parse.py` of the DR-192 runs reads the
  certificate line for acceptance and would need its index shifted for the per-iterate trace.
- MPI: `tEpsilonHead` is rank-uniform after the broadcast in `residual()`
  (`cl_FEM_DofMgr_SolverData.cpp:2420-2424`); the print is rank-0 gated as today.

## 6. Gate

- Syntax gate on both files (`-fsyntax-only -std=gnu++17` with the tree's warning flags).
- Executable, three loops: magnetic-only decks run `solve_coupled()` (`belfem.cpp:212-216`),
  so R3's magnetic line is gated by `examples/tapestack2d_layered` (minutes); R3's thermal line
  needs a fully coupled deck and R2 + R4 (`solve_magnetic()` / `solve_thermal()`, `belfem.cpp:272`,
  `:283`) need the same deck with `coupling : segregated` — the cheapest thermal example
  (`examples/corc_twolayer`, `tapestack3d`, `tapestack3d_quench`, `tape_quench_usermat`, all 3D),
  run for two or three timesteps each in this session's private HEAD tree (`treeB`, idle since the
  DR-192 runs) and killed; the lines are read against §3 and §4. Then Christian's `make check-fast`.
  A loop whose deck does not run in this session is reported as not gated, not as gated (Codex).
- No test is added for these lines: none of the existing tests reads them, and a stdout-gold
  assertion on a log line would pin formatting rather than behavior.

## 7. Open questions

- **O1** keep `print_line`'s deletion in this diff, or leave the dead helper for a later sweep?
  Plan: delete (no dead code).

## 8. Steps

- [x] R1 helper signatures — landed 2026-10-01
- [x] R2 `iterate_magnetic` — landed 2026-10-01
- [x] R3 `iterate_coupled` + delete `print_line` — landed 2026-10-01
- [x] R4 `iterate_thermal` — landed 2026-10-01
- [x] R5 theory doc paragraph — landed 2026-10-01
- [x] R6 devlog, index, DR-192 cell — 2026-10-01

## 9. Audit trail

- `tmp/ai_exchange/review_picard_line_head_residual.md` — v1 plan round (Codex `gpt-6-astra`/high,
  Grok `grok-4.7`/high, both refuting v1); v2 plan round (astra/xhigh, grok-4.7/xhigh); code round
  to follow on the own-files diff.


---

## 10. HISTORICAL — plan v1 (withdrawn 2026-09-30 by its plan jury)

> **HISTORICAL (reconciled 2026-09-30).** v1 proposed printing a "head" residual beside a "lagged" one on the body line. Codex astra and Grok 4.7 independently showed the premise wrong: on Picard `solve_from_residual` restores the pre-update residual (`cl_FEM_DofMgr_SolverData.cpp:2305`, `:2350`), so the body line already prints the head of the state the solve started from and the two cells would carry one number; only Newton recomputes a lagged post-update residual (`:2291-2298`). Five v1 line cites were also wrong (two would have put a store inside a NaN path), the width invariants were false on sentinels, and the −110 dB specimen was −70.1 dB. Verification and reconciliation: `tmp/ai_exchange/review_picard_line_head_residual.md`. The v1 text follows unchanged.

> **Scope guards:**
> - No change to any controller decision: `mEpsilon`, `mResidualHistory`, the stagnation guard,
>   the watchdog, the divergence strikes, Anderson and the ω growth all keep reading the lagged
>   residual exactly as today.
> - The circuit line (`cl_FEM_Controller.cpp:2540`, `Step %2u, residual= …`) is a different loop
>   and is out of scope.
> - The 76-column box geometry of the log is kept; every line stays 76 display columns.

---

## 1. The two numbers (source trace, confidence high)

| where | what is measured | who reads it |
|---|---|---|
| top of every trip — `iterate_coupled` `cl_FEM_Controller.cpp:1099-1101`, `iterate_magnetic` `:2107-2110`, thermal `:1603` (coupled) and `:2358-2360` (`iterate_thermal`) | `compute_jacobian_and_rhs()` at the current iterate, `compute_residual()`, `tEpsilonHead = residual( mIteration )` — `‖A_new x − b_new‖ / ‖b_new‖` | the certified exit (`:1113`, `:1612`, `:2125`, `:2377`); printed only on the `timestep succeeded` / `magnetic converged` line (`:2140-2146`, `:1747-1753`) |
| body, after `solve_from_residual()` — `:2180-2181`, `:2420`, and the coupled body | `mEpsilon = residual( mIteration )`: the post-solve recompute `r = A_old x_new − b` with the operator assembled at the previous iterate (`cl_FEM_DofMgr_SolverData.cpp:2280-2300`, "same lagged-operator semantics as the Picard branch") | `print_line*` (`:3578-3620`), `mResidualHistory` (`:2199`), every guard |

For a strongly nonlinear operator (n = 25 power law near J/Jc ≈ 0.87) the lagged residual of the
last body iterate sits at 1e-2 to 1e-3 while the head residual of the same iterate, measured one
trip later, is 1e-7 (`corc_periodic_bc`, 2026-09-29, three 1 ms runs).

## 2. Mechanism

- **R1** `cl_FEM_Controller.hpp` (member block `:289-305`): two members beside `mEpsilon` /
  `mEpsilon2`,
  ```cpp
  //! head residual of the iterate the last body solve started from
  //! ( measured with the freshly assembled operator; print only )
  real mEpsilonHead  = BELFEM_REAL_MAX ;
  real mEpsilonHead2 = BELFEM_REAL_MAX ;
  ```
- **R2** assign them at the four head sites, right after the NaN guard so a NaN head never
  reaches the print: `:1101` → `mEpsilonHead = tEpsilonHead` (after `:1107-1111`),
  `:1603` → `mEpsilonHead2` (after `:1606-1610`), `:2110` → `mEpsilonHead` (after `:2115-2119`),
  `:2360` → `mEpsilonHead2` (after `:2363-2367`). Every print site runs after its head site in the
  same call (`print_line` `:1785` after `:1101`; `print_line_magnetic` `:2195` after `:2110`;
  `print_line_thermal` `:2436` after `:2360`), so the member is never printed unassigned.
- **R3** the three print helpers (`print_line` `:3578-3600`, `print_line_magnetic` `:3606-3618`,
  `print_line_thermal` `:3621-3633`) print
  ```
     │ Magnetic Picard    5, head 4.37e-03 ( -23.6 dB), lag 6.9e-03, ω 0.500 │
     │ Magnetic Picard    5, head 9.86e-08 (-110.0 dB), lag 4.4e-03, ω 1.000 │
     │ Thermal  Newton    5, head 4.37e-03 ( -23.6 dB), lag 6.9e-03, ω 0.500 │
  ```
  format `"   │ %-8s %-6s %4u, head %8.2e (%6.1f dB), lag %7.1e, ω %5.3f │"` with
  `10·log10( mEpsilonHead )` in the dB cell — 76 display columns, measured by script against the
  existing line and the box rules. `head` is always scientific: `residual_string()` (`:630-635`)
  prints `0.000000` for anything below 5e-7, which is exactly the range the head lives in. The
  lagged residual keeps two significant digits; its dB is dropped (it is the number nobody should
  be reading for convergence). ω keeps three decimals (the growth rule moves it in the second).
- **R4** `src/fem/kernel/doc/nonlinear_controller_theory.md`: a short "reading the iteration
  line" paragraph — `head` is the true residual of the state the body solve started from, `lag` the
  linearized post-solve estimate the controller acts on, and the certified line carries the head of
  the last iterate. Doc text, no jury.
- **R5** `devlog`, `todo/README.md`, DR-192 status cell (item 2).

## 3. Semantics the reader must know (stated in the doc, R4)

The head on line `k` belongs to the iterate the body solve `k` started from (`x_{k−1}`), the lag to
the iterate it produced (`x_k`); the head of `x_k` appears on line `k+1` or on the certified line.
This one-line offset is inherent: the head of `x_k` does not exist until the next trip assembles
at `x_k`. Printing the body line one trip late instead (so that head and lag name the same iterate)
would lose the line on every reset path (NaN head, solver failure, cut) and reorder it around the
promotion label — rejected.

On a coupled no-body trip (magnetic certified, thermal body runs) `mEpsilon` already equals the
head (`:1137`), so `head` and `lag` print the same value — correct, not a defect. On a frozen
thermal trip the thermal line prints the stale `mEpsilon2` today and will print the stale
`mEpsilonHead2` too — unchanged behavior.

## 4. What could break

- Nothing functional: no decision reads the new members. Confidence high (R2 adds stores only).
- Line width: `%4u` holds iterations to 9999 (`max iterations` decks use ≤ 200); `%6.1f` holds
  dB down to −999.9; `%5.3f` holds ω ≤ 9.999 (`mOmegaMax` ≤ 1 on every path).
- Consumers of the old format: none tracked — `git grep` over `scripts/`, `python/`, `tests/`,
  `doc/` finds no parser of `Magnetic Picard` / `, residual `; the only parser is the 2026-09-29
  scratchpad `parse.py`, ephemeral.
- `ω` is a two-byte UTF-8 character in a line that already carries box-drawing characters and
  `Δ` (`print_header`), so the terminal accounting is by display column already.

## 5. Gate

- Syntax gate on the two files (`-fsyntax-only -std=gnu++17`, the tree's warning flags).
- Executable: rebuild the private HEAD tree of this session (`treeB`, already configured) with
  the diff applied once the DR-192 comparison run has released it, run `examples/2D_Tapestack`
  for a few steps, and read the new lines against the box. Then Christian's `make check-fast`.
- Nothing in `tests/` asserts on the log line; no test is added — a test that parses stdout
  would be the same false-precision class as T-67's log-driven assertions.

## 6. Open questions

- **O1** the labels: `head` / `lag` (the source's own words, chosen here) vs `true` / `est.` —
  Christian's call, no code consequence.

## 7. Steps

- ~~[ ] R1 header members~~ (v1 step, withdrawn with the plan)
- ~~[ ] R2 four head-site stores~~ (v1 step, withdrawn with the plan)
- ~~[ ] R3 three print helpers~~ (v1 step, withdrawn with the plan)
- ~~[ ] R4 controller theory doc paragraph~~ (v1 step, withdrawn with the plan)
- ~~[ ] R5 devlog, index, DR-192 cell~~ (v1 step, withdrawn with the plan)

## 9. Re-opened design (after the plan jury)

What the log actually shows on Picard, per trip k: the head residual of x_k is measured, and if
not converged the body solves and prints that same number as `Magnetic Picard k+1, residual`.
So the residual column is state-indexed but one behind the solve count, and the accepted
iterate's residual appears only on the `timestep succeeded` line. Options:

- **B — state-indexed line at the head.** Print `Magnetic <alg> k, residual head(x_k)` right
  after the head measurement of every trip (k = solves so far), and drop the body print. The
  column then ends with the accepted iterate. Newton loses its post-update lagged number from
  the log (the guards still use it); the coupled loop prints the magnetic and thermal lines at
  their own heads instead of as a pair.
- **A — relabel only.** Keep the print where it is; index the line by the state it describes
  (`mIteration` before the increment) so `Picard 4, residual 4.4e-03` followed by
  `timestep succeeded 9.9e-08` reads as states 4 and 5.
- **C — doc only.** The theory doc already says pre-update; add the one-behind reading there.

## 8. Audit trail

- `tmp/ai_exchange/review_picard_line_head_residual.md` — plan round (Codex `gpt-6-astra`/high,
  Grok `grok-4.7`/high), code round to follow on the own-files diff.
