# Harvest the electricalCircuit Demo Circuits into Transient Unit Tests

**Date:** 2026-08-28
**Purpose:** The four demo circuits hardcoded in `src/executables/electricalCircuit.cpp` are the only
exercise the circuit module's *transient* path has ever had, and three of the four are commented out
against an API that no longer exists. Harvest all four into `tests/circuit` as fixed-timestep
transient tests with analytic gates, add one netlist-equivalence test tying the ngspice front-end to
those gates, then delete the executable — which removes the second copy of the circuit Newton loop.
~~and the whole of DR-118 with it~~ (struck 2026-08-30: DR-118 was fixed in place 2026-08-29 and
struck on its own gate, so deletion no longer discharges it — see the re-scope note).
**Module:** `tests/circuit` (+ `src/executables`, `config`, docs)
**AIs involved:** Claude (exploration + plan), Codex (audit), Grok (third voice)
**Status:** **COMPLETE 2026-09-04 — R0/R2/R4/R5/R6 landed and RAN green in Christian's `make check`, R7 deleted the executable, R8/R9/R10 done, DR-39 STRUCK; R1′ not needed (P1 held, 11 iterations worst case). Post-deletion link gate (`make check` on the tree without the executable) still owed — see `devlog/dl20260904_circuit_demo_harvest.md`.** Previously: ADOPTED 2026-08-30 (Christian's ruling on DR-39(a)), RE-SCOPED against the tree the
same day, PROMOTED to the DR-39 closing plan (second ruling, below), and AUDITED 2026-08-30 by a
blind Codex + Grok jury round whose three P0s and six P1s are applied throughout — see §8.** Deletion is
the chosen discharge for DR-39(a); the shared-`CircuitSolver` shape is refuted for good. Two steps
(R1, R3) are already satisfied in-tree by work that landed after this plan was written and are
struck; the §1 DR-118 premise expired the same way. What remains is R2/R4/R5/R6 (with R1′
**conditional**, see below), then R7–R10. **This plan is now the whole of what DR-39 needs to be
struck** — leg (b) is re-homed by ruling, not implemented here. No source modified by this plan
yet; implementation is to be dispatched to a fresh session against the audited text.
**2026-09-04: R0 done; R2/R4/R5/R6 WRITTEN into `tests/circuit/test_ElectricalCircuit.cpp`
(`ResistiveSinePinsCurrentConvention`, `FredericCircuitTwoRegimes`,
`DiodeBridgeFullWaveRectifies`, `NetlistTwinMatchesFredericCircuit`) plus the
`aNumIterations` out-parameter on `solve_attempt`/`take_step`; syntax-checked against the debug
tree's flags, NOT yet run — their boxes stay open until Christian's `make check`. R1′ untouched
pending the R5 iteration count. R7–R10 not started.**
**First run, 2026-09-04: R2, R4, R6 GREEN on the first execution (so P2, P3, P4 hold). R5 RED at
step 100 — and the cause is neither of P7's two candidates.** Probe (scratchpad, linked against
the prebuilt library): with phase 0 exactly the four samples where the source is ~1e-13 V fail
(steps 100/200/300/400), every other step converges in ≤ 9 iterations, and halving Δt fails the
same way at its own zero-crossing samples. At those samples the exact solution is ~1e-14 in every
unknown (the bridge has no reactive element) and `ElectricalCircuit::residual()` =
`norm(mRHS)/mRHSnorm` with `mRHSnorm = norm(J·x)` (`cl_ElectricalCircuit.cpp:988-1001`) divides
roundoff by roundoff — the Newton plateaus at ~1e-7 although the solution is converged to machine
precision. **This is a library-level normalization singularity, not a Newton failure; neither a
smaller Δt nor R1′ can touch it. Filed for the register in R9 (new row), not fixed here — the
same `residual()` feeds `Controller::solve_circuit()`.** Test response: a half-step phase offset
on the source so the zero crossings fall between samples; with it R5 converges at every step,
worst case **11 iterations at step 1** (cold start), 9 at each commutation. **P1 holds → R1′ is
NOT implemented.**

> **Re-scope note, 2026-08-30 (tree-checked, static).** Three things changed under this plan
> between 2026-08-28 and today:
> 1. **DR-118 was fixed in place on 2026-08-29, not deleted** — `electricalCircuit.cpp:113-193` is
>    now controller parity (`set_timestep → shift → compute_MNA_matrix` per attempt, reject =
>    `shift_back` + halve + `continue`, dt floor). The §1 row "Duplicate Newton loop, **and a broken
>    one**" is struck below. The deletion case survives on different grounds and is *not* weaker:
>    the binary is uninvokable by a user, and three of its four circuits reference a `SourceType`
>    enum that exists nowhere in the tree.
> 2. **The transient driver and the RLC gate already exist** — `tests/circuit/test_ElectricalCircuit.cpp:56-81`
>    (`solve_attempt`/`take_step`) and `:84-99, 535+` (`create_rlc_circuit`,
>    `RLCRingAcrossRejectedStep`). R1 and R3 are struck; R1′ replaces R1.
> 3. **The netlist front-end shipped 2026-08-27**, so R6 builds on production code rather than on
>    a promise.
>
> **The one coverage hole that makes the ordering non-negotiable:** `create_diode` appears in the
> entire test suite only at `test_CircuitAdjacency.cpp:257`, structurally, never solved. R5 is the
> first transient exercise the diode Newton will ever have, and it must land before the demo — the
> only other artifact that ever ran one — is deleted.
>
> **Freeze consequence, recorded not hidden:** R7/R8 change `BELFEM_INSTALL_EXECUTABLES` and four
> doc sites, i.e. a user-visible contract, on the 2026-08-30 design-freeze day. Ruled anyway, in
> preference to shipping the stub through 1.0 undecided.

> **Closing ruling, 2026-08-30 (Christian) — this plan closes DR-39 outright.**
> Asked what "close DR-39" covers, Christian ruled: **execute (a), re-home (b), then strike the
> row.** Concretely:
> - **(a)** is this plan, R2/R4/R5/R6 then R7–R10.
> - **(b)** `Pulse`/`PWL` is ruled a **feature, not debt**. It already has a home
>   (`todo/closed/ngspice_parser_plan.md` Phase 6) and a user-visible behaviour today (the netlist parser
>   names both and refuses them, `cl_NgspiceCircuitFactory.cpp:430-436`), so nothing is silently
>   lost by taking it off the register. R9 records the re-homing and strikes DR-39.
> - **(c)** was discharged 2026-08-27 when the netlist front-end shipped.
>
> **Recommendation carried into the re-homing, Claude, medium confidence — Phase 6 should split.**
> The two SPICE source types are not equally relevant to what BELFEM is for.
> `PULSE(v1 v2 td tr tf pw per)` is a trapezoid: ramp → plateau → ramp → repeat, which *is* the
> standard HTS magnet duty cycle and the standard AC-loss cycle, and BELFEM can only approximate it
> today with `Ramp` or `Sigmoid`. `PWL` is the general breakpoint table, whose only magnet use is
> replaying a *measured* current profile — for which no BELFEM deck syntax exists either — and it
> is the leg that costs the most: its variable-length payload is exactly what the MPI rule routes
> through `share`/`receive` rather than `broadcast`, and its storage must be `Cell< real >` because
> `cl_SourceFunction.hpp:43-44` bans `cl_Vector.hpp` on plugin-ABI grounds. Proposal: fund `PULSE`,
> gate `PWL` on a real deck that needs it. **This is a recommendation to Phase 6, not a decision
> taken here** — R9 writes it down, it does not act on it.

> **Scope guards (from the task brief):**
> - IN scope: the four demo circuits, a shared fixed-Δt transient driver in the test TU, one netlist
>   equivalence deck, deletion of `electricalCircuit.cpp` and its install entry, the doc and
>   debt-register follow-through.
> - OUT of scope: implementing DR-39(b) `Pulse`/`PWL` — re-homed by the closing ruling above to
>   parser plan Phase 6, and R9 is the only step that touches it; extracting a shared
>   `CircuitSolver` class (DR-39(a)'s original shape — deletion discharges it more cheaply); adaptive
>   timestepping (it belongs to `Controller`, not to `ElectricalCircuit`); any *fix* to DR-121,
>   DR-123 or DR-124 — the tests are laid out to avoid those defects, not to cover them.
> - Compatibility: `electricalCircuit` disappears as a user-facing binary. Nothing in `src/`,
>   `examples/` or any shipped deck invokes it (verified: the only non-doc references are its own
>   two CMake lines and `config/globals.cmake:14`).

---

## 1. Current Behaviour and How It Fails

`src/executables/electricalCircuit.cpp` is a **205**-line driver (it grew with the DR-118 rewrite;
"192" was the 2026-08-28 count) that builds circuits **in C++** and writes `circuitAnalysis.txt`.
One circuit is live (`:46-65`); three are commented out (`:67-88`). It is built (`src/executables/CMakeLists.txt:44-46`) and **installed**
(`config/globals.cmake:14`).

| Failure | Mechanism | Evidence |
|---|---|---|
| Three of four circuits cannot compile | they call `create_voltage_source( SourceType::SINE, amplitude, freq, n1, n2 )`; the live API is `create_voltage_source( SourceFunction *, n1, n2, label )` | `electricalCircuit.cpp:69,76,82` vs `cl_ElectricalCircuit.hpp:209` (both anchors corrected 2026-08-30 by audit) |
| The R-only demo would build a singular circuit | it is written for the shared `ElectricalCircuit(4)` at `:46` but only references nodes 0 and 1, leaving 2 and 3 floating | `electricalCircuit.cpp:75-78` |
| The diode bridge rectifies into nothing | the four diodes form a correct bridge (`0→2, 3→1, 3→0, 1→2`) but there is no load between the DC terminals 2 and 3 | `electricalCircuit.cpp:82-88` |
| ~~Transient behaviour is gated by eyeballing a text file~~ **half stale, corrected 2026-08-30 by audit.** The demo still asserts nothing. But "the only automated circuit coverage solves to a single DC operating point" is **false**: `RLCRingAcrossRejectedStep` (`test_ElectricalCircuit.cpp:535`) is a 3000-step transient, and the switch latch has `SwitchLatchRevertsOnRejectedStep` (`:158`). **The genuine hole is the diode Newton, and only the diode Newton** | `test_ElectricalCircuit.cpp:158, 535`; DC-only cases at `test_NgspiceFactory.cpp:160-171`, `test_HybridFactory.cpp:47-70` |
| ~~Duplicate Newton loop, and a broken one~~ **STRUCK 2026-08-30 — half stale.** The loop is still duplicated (`:134-164` vs `cl_FEM_Controller.cpp:2528-2585`), and that duplication is DR-39(a). It is no longer *broken*: DR-118 was fixed in place 2026-08-29 | ~~its reject path is DR-118~~ — the reject path is now controller parity; the surviving defect is duplication alone, and the two copies are provably equivalent (`residual()` = `norm(mRHS)/mRHSnorm`, `cl_ElectricalCircuit.cpp:997-1001`, holds no state) | `cl_FEM_Controller.cpp:2528-2585`; `electricalCircuit.cpp:134-164` |

**Note on step numbering:** the IDs below are historical and deliberately not renumbered — R1 and
R3 are struck in place, R0 was added 2026-08-30 ahead of R1′, and R2/R4/R5/R6/R7–R10 keep the
numbers the 2026-08-28 plan gave them. Renumbering would break every cross-reference in the devlogs
and in the DR-39 row.

**Bottom line, rewritten 2026-08-30 after the audit round refuted the old one.** The old text said
BDF history, the switch latch and the diode Newton are "exercised by exactly one artifact, whose own
timestep control is a filed defect". Both halves are now false: **BDF has a gate**
(`RLCRingAcrossRejectedStep`), **the latch has a gate** (`SwitchLatchRevertsOnRejectedStep`), and
DR-118 was struck 2026-08-29 so the demo's timestep control is no longer defective. What survives is
narrower and sufficient: **the diode Newton has no gate at all** — `create_diode` is solved nowhere
in the suite — and the demo is the only artifact that has ever solved one. R5 closes that hole; the
other three steps preserve demo content that would otherwise be deleted unrecorded.

## 2. Architecture: Why the Direct API Is the Right Spine

The gates are about MNA assembly, BDF stamping, the switch latch and the diode Newton — the layer
under `create_resistor`/`create_inductor`/… . Building through `NgspiceCircuitFactory` instead would
route every gate through the parser, so a failure could not distinguish a physics defect from a
lexing defect. So: **four gates on the direct API, plus one equivalence test** (R6) proving the
netlist front-end reproduces one of them exactly, which is where parser coverage belongs and which
additionally exercises the O5 name→index lookup, the O8 `order` directive and the `switch` directive.

Rejected: netlist-only (cannot localise a failure); direct-API-only (leaves the front-end tied to
the physics gates by nothing). *("four gates on the direct API" above was written before R3 was
found already in-tree; the direct-API gates are R2, R4 and R5.)*

**Fixed Δt throughout.** The demo's adaptive control (`:180-193`) is the executable's; `Controller`
owns the real policy, and these gates are about MNA/BDF/latch/Newton, not about step control.
~~its reject path is DR-118~~ and ~~keeps every test away from DR-123, whose trigger is a switch
firing inside a *rejected* step~~ — **both struck 2026-08-30: DR-118 and DR-123 were fixed and
struck on 2026-08-29, and the rejected-step switch case is now itself a test**
(`SwitchLatchRevertsOnRejectedStep`, `test_ElectricalCircuit.cpp:158`). Fixed Δt stays the choice
on its own merits: a gate that also varies the step cannot say which of the two it failed on.

## 3. Gap Table

| # | State / behaviour | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| ~~1~~ | fixed-Δt transient driver (shift → Newton) | all four gates | **YES, since 2026-08-29** — `solve_attempt`/`take_step`, `tests/circuit/test_ElectricalCircuit.cpp:56-81` | ~~(c)~~ **satisfied** | residue is the ω update only, which `take_step` does not carry (`set_omega(1.0)`, no relaxation) → **R1′** |
| 2 | resistive circuit, sine source | exact algebraic gate | commented out, wrong API | (a) rebuild from `:75-78` | 2-node circuit, ground = last index (`cl_ElectricalCircuit.cpp:38-46`) |
| ~~3~~ | series RLC at resonance | BDF order > 1 under a dynamic load | **YES, since 2026-08-29** — `create_rlc_circuit` (order 2, not 4) + `RLCRingAcrossRejectedStep`, `test_ElectricalCircuit.cpp:84-99, 535+` | ~~(a)~~ **satisfied** | the in-tree version gates the ring across a *rejected* step rather than the resonance amplitude; judged sufficient — do not rebuild it |
| 4 | switched RL → RL‖C, BDF3 | the `Switch` latch on the real time axis | live demo `:45-64` | (a) rebuild verbatim | `Switch::shift` fires at the first `mTime ≥ mSwitchTime` (`cl_Switch.cpp:47-58`, latch at `:53-57`) |
| 5 | diode bridge + **load resistor** | Newton on four exponentials | bridge present, load absent | (c) — add `create_resistor( 100, 2, 3 )` | decided 2026-08-28, Christian; see §5 O1 |
| 6 | netlist deck for row 4 | front-end ↔ physics tie | no | (c) | `* belfem: order` and `* belfem: switch` directives, `cl_NgspiceCircuitFactory.hpp:62-70` |
| 7 | `electricalCircuit` build + install entries | deletion | present | (c) | `src/executables/CMakeLists.txt:44-46`, `config/globals.cmake:14` |
| 8 | doc references to the executable | deletion | ~~four live references~~ **THIRTEEN, established by the 2026-08-30 audit round** | (c) | the full inventory is in R8 — do not work from this cell |

### 3.1 Cross-cutting findings

- **Ground is the last node index**, not node 0 (`cl_ElectricalCircuit.cpp:38-46`); every harvested
  circuit must be sized so its return node is the highest index. This is what makes the R-only demo
  singular under the shared `ElectricalCircuit(4)`.
- **BDF order ramps by itself.** `compute_coefficients` sizes the Vandermonde from the *filled*
  history count, `n = mH.size() + 1` (`cl_BDF.cpp:33`), so a BDF3/BDF4 element started from rest is
  order-1 on step one and climbs. A cold start at order 4 is therefore safe — and correct, since the
  zero history *is* the true solution for t < 0. (high confidence)
- **Component creation order is load-bearing.** DR-121: a `Switch` created before a voltage source
  shifts every later V-branch MNA stamp. Circuit 4 creates V first, as the demo did; the test must
  say so in a comment or a later reorder will silently trip a known defect.
- ~~**Newton relaxation is not optional for the diode gate.**~~ **DOWNGRADED 2026-08-30 to
  *protective, and therefore conditional* (medium confidence).** The overflow risk is real and now
  quantified: `Diode::compute_current` evaluates `exp( ( v+ − v− ) / Vt )` (`cl_Diode.cpp:27-31`),
  so at Vt = 0.026 V any Newton iterate driving a diode above **≈ 18.4 V** returns `inf` and the
  solve is lost. But the iterate does not start from nothing — every step warm-starts from the
  previous converged `mX`, and at Δt = 1/(60·200) the 10 V / 60 Hz source moves by at most
  10·2π·60·Δt ≈ **0.31 V per step**, so the Newton start is close and a full step is unlikely to
  overshoot by 18 V. Note also that the demo's own ω update is **a no-op on the first iteration**:
  it enters with `epsilon0 = BELFEM_REAL_MAX`, takes the `epsilon < epsilon0` branch, multiplies by
  `beta + gamma·(2/π)·atan(1) ≈ 1.3`, and is clamped back to 1 (`cl_FEM_Controller.cpp:2564-2574`)
  — so relaxation cannot protect the *first* iterate of a step in any case. This is why **R1′ is
  now conditional on R5 actually failing** (§4), with the decision rule pre-registered in §6.

### 3.2 Findings added 2026-08-30 while making the plan implementation-ready

All four were read out of the tree in the same pass; each changes a gate or removes a trap.

- **The voltage-source branch current is the NEGATIVE of the delivered current (high confidence).**
  MNA node rows here are *current leaving the node*: a resistor's conductance stamp is
  `cl_ElectricalCircuit.cpp:645-670` and its **residual** contribution `+G·(v+ − v−)` is `:761-781`;
  a voltage source stamps `+i` into the same residual row (`:816-821`). The residual block is the
  evidence that fixes the sign of the *solved* current — both cites are given because the audit
  round found the single `:761-781` cite ambiguous. So for a source driving a load from `n+` to `n−`, `i < 0`. Cross-checked against
  the in-tree divider test, which passes with v(0) = 5, v(1) = 2.5 through R1 = R2 = 1 kΩ: the node-0
  KCL is `2.5 mA + i = 0`, i.e. **i = −2.5 mA**. **R2 and R4 must gate on this sign explicitly**
  rather than on `std::abs`, so the test pins the convention instead of hiding it; a future sign
  regression would otherwise pass silently.
- **`t_switch` on an exact step boundary is a floating-point tie, not a determinism guarantee
  (high confidence).** `ElectricalCircuit::shift()` advances `mTime += mDeltaTime` and then passes
  the **end-of-step** time to `Switch::shift`, which latches on `aTime >= mSwitchTime`
  (`cl_ElectricalCircuit.cpp:195-221`, `cl_Switch.cpp:46-57`). `mTime` after 500 accumulated
  additions of `50e-6` differs from a directly-computed `1.25/50` at the ~1e-18 level with an
  unpredictable sign, so a `t_switch` sitting *exactly* on a boundary fires at step 500 or step 501
  depending on the compiler and optimisation level — and, worse for R6, possibly *differently* in
  the hand-built and netlist runs. The fix costs nothing: put `t_switch` half a step early. See R4.
- **Resistor and diode currents ARE live after `solve()` (high confidence).** `ElectricalCircuit::solve()`
  calls `compute_current()` on every RESISTOR / INDUCTOR / CAPACITOR / TERMINALPAIR / DIODE /
  SUPERCONDUCTOR after `update_components()` (`:962-986`), so R5 may read the bridge load current
  straight from `get_current_on_component()`. Only the unknown-current components (voltage source,
  switch) take their current from the solution vector; nothing in the plan needs a hand-computed
  Ohm's law.
- **An open switch is not a singularity (high confidence).** The open branch stamps `i_sw = 0` on
  its own row (`cl_ElectricalCircuit.cpp:879-886`), so R4's pre-switch circuit — in which node 2
  hangs off nothing but the capacitor and the open switch — is well posed: the capacitor's
  discretized conductance ties node 2 to ground.

## 4. Ordered Steps

- [x] ~~**R1** — `tests/circuit/test_TransientCircuits.cpp`, plus its entry in
      `tests/circuit/CMakeLists.txt` SOURCES. Anonymous-namespace helper: `advance( circuit, dt )`
      = `set_timestep` → `shift` → `compute_MNA_matrix` → relaxed Newton to 1e-10, returning the
      iteration count so a test can assert convergence.~~ **SUPERSEDED 2026-08-30 — the driver
      landed elsewhere first.** `tests/circuit/test_ElectricalCircuit.cpp:56-81` already has
      `solve_attempt` + `take_step` doing exactly this sequence. A second TU would duplicate the
      helper, which is a poor look on a de-duplication row.
- [x] **R0** *(first, and it is not optional; done 2026-09-04 — all four cases went into the existing TU, no new file, no `CMakeLists.txt` change)* — **Harvest into the existing
      `tests/circuit/test_ElectricalCircuit.cpp`, not a new TU.** R1 was going to create
      `test_TransientCircuits.cpp`; a second TU would duplicate `solve_attempt`/`take_step`, which
      is an indefensible look on a *de-duplication* row. No new file, no new `CMakeLists.txt` entry —
      `test_ElectricalCircuit.cpp` is already in `tests/circuit/CMakeLists.txt:11` and the suite is
      already `TESTLABELS fast` (`:3`). Every step below adds `TEST(...)` cases and, where needed,
      anonymous-namespace builders next to `create_rlc_circuit`. Budget check: R2 (400 steps) +
      R4 (1600) + R5 (400) + R6 (1600) ≈ 4000 solves on systems of order 4–6, which stays inside
      `fast`.
- [x] ~~**R1′**~~ *(replaces R1; **CONDITIONAL — attempt R5 first**; **NOT NEEDED, 2026-09-04: R5 converged at every step with ω ≡ 1, worst case 11 iterations — struck on its own decision rule**)* — `solve_attempt` pins
      `set_omega( 1.0 )` and never relaxes (`test_ElectricalCircuit.cpp:61`). The 2026-08-28 plan
      asserted R5 could not converge without relaxation; §3.1 downgrades that to *protective*, and
      §3.2 shows the demo's own ω update is a **no-op on the first iteration** anyway, so it cannot
      protect the one iterate most at risk. **Decision rule, pre-registered in §6: build R5 against
      the existing ω ≡ 1 helper. Only if R5 fails to converge does R1′ land** — and then it lands as
      the ω update from `cl_FEM_Controller.cpp:2564-2574`, added to `solve_attempt` behind a
      **defaulted** parameter so all existing callers keep ω ≡ 1 and their current results
      bit-for-bit. Adding an unexercised branch to a shared test helper on the strength of a
      prediction that turned out false would be worse than not adding it. If R1′ is skipped, say so
      in the devlog with the observed iteration count — that is the evidence, not the silence.

      **Two corrections from the 2026-08-30 audit round, both adopted:**
      - **R1′ is NOT the response to an overflow.** `ElectricalCircuit::solve()` writes
        `mX -= mOmega * tdX` and returns (`cl_ElectricalCircuit.cpp:948-960`); `solve_attempt` has
        no rollback. Once `mX` holds `inf`, relaxation only scales an update against a poisoned
        state. **If R5 produces `nan`, the response is a smaller fixed test Δt (or explicit
        first-step damping), not ω.** R1′'s trigger is *slow or non-convergence*, not `nan`.
      - **`solve_attempt` returns `bool` only** (`test_ElectricalCircuit.cpp:57-72`), so the
        iteration count R1′'s decision rests on is not observable today. Add an **optional
        out-parameter** (`uint * aNumIterations = nullptr`) as part of R5 — independently of whether
        relaxation is ever added — or the §6 P1 evidence cannot be recorded at all.
- [x] **R2** *(after R0)* — **Resistive.** `ElectricalCircuit( 2 )` — node 0 live, node 1 ground
      (ground is the LAST index, §3.1; its voltage is 0.0 from the `ElectricNode` constructor,
      `cl_ElectricNode.cpp:20-24`, and `update_components()` never writes it — so asserting
      `v( ground ) == 0` is a tautology, not a gate; the real gate is the *difference*).
      `create_voltage_source( sine 1 V / 50 Hz, 0, 1, "V1" )` then
      `create_resistor( 1.0, 0, 1, "R1" )`. Δt = 50 µs, 400 steps (one period). The circuit is purely
      resistive and carries no discretisation error at all, so the gates below are exact — but
      **they must be mixed absolute/relative, not purely relative** (audit round, 2026-08-30):
      the window contains the sine's zeros exactly. The source sample is `sin( π k / 200 )`, which
      is **0 at k = 200 and k = 400**, so a relative comparison there divides by ~1e-16 and means
      nothing. Use `EXPECT_NEAR( got, want, 1e-12 * std::max( 1.0, std::abs( want ) ) )` or an
      equivalent absolute floor. Gates at **every** step:
      - `v(0) − v(1)` equals the source sample. Take the source sample from the component
        (`circuit.component( 0 )->get_value()` — `get_value()` is virtual on `Component`,
        `cl_Component.hpp:71-72`, so **no cast and no extra include are needed**), not from a
        re-evaluated `sin()` in the test — re-deriving it would let a source-function defect and a
        stamping defect cancel. Note `get_value()` returns the value written by the last `shift()`
        (`cl_VoltageSource.cpp:52-56`), i.e. the source at the END of the step, which is the sample
        the stamp used.
      - **`get_current_on_component( 0 ) == −( v(0) − v(1) ) / R`, with the minus sign asserted, not
        absorbed into `std::abs`** (§3.2). This is the step that pins the branch-current convention;
        it is the only gate in the plan that would catch a future sign regression.
      - `get_current_on_component( 1 )` (the resistor, live after `solve()` — §3.2) equals
        `+( v(0) − v(1) ) / R`, i.e. **opposite in sign** to the source branch. The pair is the
        whole point: same magnitude, opposite convention.
- [x] ~~**R3** *(after R1)* — **Series RLC at resonance.** 4 nodes, sine 1 V / 10 kHz, L = C = 15.915e-6
      (both BDF4), R = 1 Ω, Δt = 1 µs (100 steps/period), 4 periods. ζ = 1/(2Q) = 0.5, so the transient
      is down by e⁻¹²·⁶ ≈ 3e-6 well before the last period. Gates over the final period:
      peak \|i\| = 1.000 A ±2 %, and the resistor carries the **whole** source voltage
      (`v(2) − v(3)` tracks the source), which is the resonance condition stated without extracting a phase.~~
      **ALREADY IN TREE 2026-08-30:** `create_rlc_circuit` + `RLCRingAcrossRejectedStep`
      (`test_ElectricalCircuit.cpp:84-99, 535+`) build the same circuit at BDF **order 2** and gate the
      ring across a rejected step instead of the resonance amplitude. Different gate, same coverage
      intent; judged sufficient — do not rebuild it. Anything wanting the order-4 amplitude gate
      specifically should reopen this as its own row, not smuggle it in here.
- [x] **R4** *(after R0)* — **Frederic's circuit, two regimes.** Verbatim topology and values from
      `electricalCircuit.cpp:45-64`, created in **exactly this order** (V before the switch — and say so in a
      comment, **citing the regression test rather than an open defect**: DR-121 was struck
      2026-08-29 and `SwitchBeforeVoltageSourceAgrees` (`test_ElectricalCircuit.cpp:110`) is what
      would now catch a reorder. The ordering is hygiene backed by a gate, not a live landmine —
      corrected 2026-08-30 by audit):
      `ElectricalCircuit( 4 )`, ground = node 3; `create_voltage_source( sine 1 V / 50 Hz, 0, 3, "v1" )`;
      `create_inductor( 200*1e-6, 3, 0, 1, "l1" )`; `create_resistor( 1.0, 1, 3, "r1" )`;
      `create_switch( false, 24.975*1e-3, 1, 2, "s1" )`; `create_capacitor( 3.183*1e-3, 3, 2, 3, "c1" )`.
      **Labels are lowercase deliberately** — the netlist parser case-folds every token
      (`cl_NetlistParser.cpp:274-286`), so R6 can only compare labels if R4 uses the folded form
      (audit round; live precedent `test_NgspiceFactory.cpp:92,115` expects `"is"`/`"sc1"`).
      **Δt = 50 µs** (Christian's call, 2026-08-30), run to 80 ms = 1600 steps.

      **The switch time is 24.975 ms, not the demo's 25 ms, and that is deliberate — do not
      "restore" it.** §3.2: the latch tests the *end-of-step* time against `t_switch`, and 25 ms is
      exactly a 500-step boundary, so the comparison is a floating-point tie that fires at step 500
      or 501 depending on the build. 24.975 ms sits half a step early (499.5 steps), which is
      ~2·10¹³ ulp clear of any boundary — and the **physical firing instant is unchanged**: the
      first end-of-step time ≥ 24.975 ms is still 25.0 ms, the end of step 500. Determinism gained,
      behaviour identical. This also makes R6 possible at all (a 1-ulp difference between
      `24.975e-3` and `spice_number_to_si( "24.975m" )` cannot move the fire step; a 1-ulp
      difference at an exact boundary could).

      Gates on the **source branch current** `get_current_on_component( 0 )`, whose sign is negative
      by §3.2 — gate the amplitude `max| i |` over a window, which is sign-agnostic, and assert the
      sign once at a known-nonzero step:
      - last **pre**-switch period — **accepted steps 100…499**, t = 5.00…24.95 ms: amplitude
        **0.99803 A ±1 %**, from \|1 / ( R + jωL )\| with ωL = 0.0628 Ω.
      - last **post**-switch period — **accepted steps 1201…1600**, t = 60.05…80.00 ms: amplitude
        **1.50562 A ±1 %**, from \|1 / ( jωL + R‖(1/jωC) )\|. **ωRC is 0.9999689, NOT 1.0
        exactly** (audit round — 3.183 mF is a rounded `1/(2π·50)`), so R‖C ≈ 0.500016 − 0.500000j
        rather than exactly ( 1 − j )/2. Computing the expected amplitude from the live component
        values, as this step already requires, absorbs the difference; the ±1 % band is 200× larger
        than it.
      - **Two window facts that are easy to get wrong, and both were wrong in this plan until
        2026-08-30.** *(i)* The period at 50 Hz is 20 ms = **400** steps at Δt = 50 µs, not 200 —
        an amplitude gate taken over a half period can miss the peak entirely depending on phase,
        so the window must be a full 400 samples. *(ii)* **Step 500 is already post-switch.**
        `shift()` fires the latch and `compute_MNA_matrix()` then stamps with the switch CLOSED, so
        the solution recorded at the end of step 500 belongs to the new regime; the pre-switch
        window must stop at step **499**. Settling is comfortable either way: the pre-switch L/R
        constant is 200 µs (4 steps), so step 100 is ~25 τ in.
      - **compute both from the component values in the test**, not from the literals above, so a
        transcription slip cannot be the thing that passes.
      The post-switch amplitude is 51 % higher, so the gate fails loudly if the switch never fires or
      fires in the wrong regime. τ_RC = 3.183 ms, so the 60 ms window opens ~11 τ after the switch
      and the 80 ms end is ~17 τ. Discretisation is not a factor at this tolerance: ωΔt = 0.0157 and
      the BDF order climbs to 3 within three steps (§3.1), so the amplitude error is O(10⁻⁶) against
      a ±1 % band; if either gate misses by more than a few tenths of a percent, suspect the
      circuit, not the tolerance.
- [x] **R5** *(after R0; the first transient exercise the diode Newton has ever had)* — **Diode
      bridge.** `ElectricalCircuit( 4 )`, ground = node 3. The demo's four diodes unchanged
      (Is = 0.1e-3 A, Vt = 0.026 V), sine 10 V / 60 Hz across 0–1, **plus the load the demo lacked**,
      `create_resistor( 100.0, 2, 3 )` (O1). Creation order: **V first** (DR-121), then
      `create_diode( 0.1e-3, 0.026, 0, 2 )`, `( 3, 1 )`, `( 3, 0 )`, `( 1, 2 )`, then the load.
      **Creation-order note corrected 2026-08-30 by audit: R5 has no switch, so DR-121 does not
      apply to it.** Create V first anyway as house order; the rationale belongs to R4, not here.
      Bridge check, done statically so the test does not have to discover it: with v(0) > v(1) the
      path is 0→D1→2→R→3→D2→1; with v(1) > v(0) it is 1→D4→2→R→3→D3→0. Both deliver **+** to node 2,
      which is what makes it full-wave. Δt = 1/(60·200) ≈ 83.3 µs, 2 periods = 400 steps.
      Gates:
      1. **Rectification** — `v(2) − v(3) ≥ −1e-6 V` at every step. **The floor is 1e-6, not 1e-9**
         (audit round): `solve_attempt`'s tolerance is a *relative residual*,
         `norm(mRHS)/mRHSnorm ≤ 1e-9` (`test_ElectricalCircuit.cpp:59`,
         `cl_ElectricalCircuit.cpp:997-1001`), which bounds no node voltage to 1e-9 V. A 1e-9 V
         floor would flake at the commutation instants.
      2. **Full-wave** — two maxima of `v(2) − v(3)` per source period. Fails if any single diode is
         mis-stamped (a half-wave result has one). Count maxima on the second period only, and with
         a threshold at half the observed peak, so numerical wobble near zero cannot manufacture
         extra maxima.
      3. **Peak output = 9.6426 V ±1 %** — the fixed point of `v = 10 − 2·Vt·ln( v/(R·Is) + 1 )`
         (two diodes in series, each dropping `Vt·ln( I/Is + 1 )`, matching `cl_Diode.cpp:27-31`
         exactly). **Iterate it in the test**; the literal above is the expected answer, not the
         gate.
      4. **Newton health** — assert `solve_attempt` returned true at every step and, if the helper is
         extended to report it, record the worst-case iteration count in the devlog. That number is
         the evidence for or against R1′ (§6).
      **Failure-mode watch, rewritten 2026-08-30 after the audit refuted the first version.** The old
      text named IEEE overflow (`exp( ΔV / 0.026 )` above ≈ 18.45 V) as the risk and R1′ as the
      response. Both are wrong. A 10 V source cannot drive any single bridge diode to 18.45 V, so
      that mode is unreachable; and R1′ could not repair it anyway (see R1′). **The real risk is
      commutation.** At each source zero crossing the conducting pair swaps and the newly forward
      diodes start from reverse bias — the classic diode-Newton jump, which the 0.31 V/step warm-start
      argument does not cover. `dIdV = (Is/Vt)·exp( ΔV/Vt )` reaches ≈ 4·10¹⁶⁵ at ΔV = 10 V: finite,
      but the Jacobian is catastrophically ill-conditioned there. **If R5 fails or flakes, look at the
      zero crossings first.** Mitigations in order of preference: a smaller Δt near the crossings is
      not available (fixed Δt is a design choice, §2), so the options are a globally smaller Δt, then
      R1′. Is = 0.1 mA is a deliberately soft exponential, which is why P1 may still hold.
- [x] **R6** *(after R4)* — **Netlist equivalence.** R4's circuit as a `.cir` deck written to
      `::testing::TempDir()`, built through `NgspiceCircuitFactory`, driven by the same helper, and
      compared against a fresh R4 run **step by step** — not just at the end, or a transient
      divergence that re-converges would pass.

      **Do NOT use `FileGuard` — R6 would not compile on a `USE_HDF5=OFF` tree.** Both auditors
      caught this independently: `FileGuard` is defined at `test_ElectricalCircuit.cpp:203`, inside
      the `#ifdef BELFEM_HDF5` block that runs `:189`–`:521`, and HDF5 is optional
      (`CMakeLists.txt:94`). R6 writes a `.cir`, not HDF5. **Copy `write_scratch` from
      `test_HybridFactory.cpp:43-50`** (`::testing::TempDir()` + `ofstream`, no HDF5 dependency) —
      note it does not delete the file, which is acceptable for scratch but should not be mistaken
      for a guard. Take ownership of the built circuit through `NgspiceCircuitFactory::circuit()`,
      which hands over exactly once (`cl_NgspiceCircuitFactory.cpp:55-63`).

      Deck (node names chosen so first-appearance packing reproduces the hand-built indices — see
      the caution below):
      ```
      * R4 equivalence deck -- Frederic's circuit
      V1 n0 0 SIN(0 1 50)
      L1 n0 n1 200u
      R1 n1 0 1
      * belfem: switch S1 n+=n1 n-=n2 state=open t_switch=24.975m
      C1 n2 0 3.183m
      * belfem: order L1 3
      * belfem: order C1 3
      .end
      ```
      **Bit-equality is achievable and the gate should be 1e-12, but only if the test writes the
      hand-built values in the same arithmetic form** (verified 2026-08-30 against
      `fn_spice_number.cpp:124-190` and `cl_NgspiceCircuitFactory.cpp:421-426`): `200u` evaluates to
      `200.0 * 1.0e-6`, `3.183m` to `3.183 * 1.0e-3`, `24.975m` to `24.975 * 1.0e-3`, and
      `SIN(0 1 50)` to `set_periodic( Sine, 1.0, 1.0 / 50.0, 0.0 )`. So R4 must be written as
      `200*1e-6`, `3.183*1e-3`, `24.975e-3` and `1.0 / tFreq` with `tFreq = 50` — **not** as
      pre-multiplied literals such as `2.0e-4`, which differ in the last bit and would drift the
      1600-step traces past 1e-12. If the gate still misses at 1e-12, **report the discrepancy;
      do not loosen the tolerance silently** — the loosening is the finding.

      **Caution on node indices, tree-checked:** the netlist packs non-ground nodes by first
      appearance and puts ground last (O9). With the deck above that yields n0→0, n1→1, n2→2,
      "0"→3, matching the hand-built circuit exactly — but **resolve every node through
      `node_index( name )` anyway** (`cl_NgspiceCircuitFactory.cpp:200`) rather than relying on the
      coincidence. A later edit to the deck's line order would otherwise change the mapping
      silently. Component *indices* line up for the same reason (components are created in merged
      line order), so `get_current_on_component( 0 )` is the voltage source on both sides; assert
      that with `component( 0 )->get_label()` (`cl_Component.hpp:96-97` — the accessor is
      `get_label()`, not `label()`). **Compare against `"v1"`, lowercase**: the parser folds case
      (`cl_NetlistParser.cpp:274-286`), so `"V1" != "v1"` and a label assertion written the obvious
      way fails. R4 is specified with lowercase labels for exactly this reason.

      **The comparison tolerance needs an absolute floor, like R2's.** Node voltages and branch
      currents in this circuit pass through zero every half period, so a purely relative 1e-12
      comparison is undefined at those samples. Use
      `EXPECT_NEAR( b, a, 1e-12 * std::max( 1.0, std::abs( a ) ) )` or equivalent.

      **The `* belfem: order` lines are load-bearing, not decoration:** the factory's default BDF
      order is **1** (`cl_NgspiceCircuitFactory.cpp:452-453`), so without them the deck builds a
      BDF1 circuit against R4's BDF3 and the traces separate immediately. If R6 misses at 1e-12,
      check in this order: (1) the order directives took effect, (2) the arithmetic forms above,
      (3) whether the switch fired on the same step.

      Coverage delivered: `SIN(...)`, `* belfem: order` at BDF3, and
      `* belfem: switch ... state=open t_switch=`. The V card precedes the switch directive, as the
      DR-121 guard requires and as `create_switch` ordering demands.
- [x] **R7** *(after R2, R4, R5, R6 green — and "green" means Christian has RUN them, not that they
      compile)* — Delete `src/executables/electricalCircuit.cpp`; drop the two lines
      `src/executables/CMakeLists.txt:44-45` plus the `include( ... Add_Executable.cmake )` at `:46`
      that belongs to them, and remove `electricalCircuit` from the
      `BELFEM_INSTALL_EXECUTABLES` list at `config/globals.cmake:14`. Both anchors re-verified
      2026-08-30. **Deleting the file before the gates have run is the one irreversible mistake
      available in this plan** — the demo is the only artifact in the tree that has ever solved a
      diode bridge.
- [x] **R8** *(after R7)* — Docs. **The 2026-08-30 audit round found this step's inventory wrong by
      nine sites. It said four; there are thirteen.** Codex found one, Grok found four, and re-running
      the grep myself found four more — including the one that matters most. Work from the table
      below, not from §3 row 8.

      | # | site | what it says / why it matters |
      |---|---|---|
      | 1 | `CLAUDE.md:136` | executables list |
      | 2 | `src/executables/doc/README.md:19` | table row. **Also `:9-11`**, which claims every executable reads `input.conf` — false today *because of* this binary, and true once it is gone |
      | 3 | `doc/README.md:59` | executables blurb |
      | 4 | `src/circuit/doc/circuit_usage_guide.md:585-595` | §11, stale in **three** ways (below) |
      | 5 | **`scripts/update_doc_index.py:124`** | **THE GENERATOR.** It emits the `executables` blurb that produces sites 6 and 7. **Fix this one first** — patching the generated files without it means the dead name reappears on the next doc-index refresh |
      | 6 | `doc/doxygen_nav.dox:52` | generated from site 5 |
      | 7 | `doc/groups.dox:266` | generated from site 5 |
      | 8 | `Doxyfile.in:1166` | comment naming the three "public executables" |
      | 9 | `examples/scripts/README.md:63-72` | ~~**user-facing and already FALSE**: says the binary *reads* `circuitAnalysis.txt`~~ **FACTUAL ERROR FIXED 2026-08-30** (see the note below); the paragraph still names the binary and must go with it |
      | 10 | `examples/scripts/Allrun:120-125` | **a shipped script**, not a doc: a block comment telling maintainers not to add a circuit heuristic because this binary ignores `input.conf`. It does **not** invoke the binary — the compatibility claim survives — but the comment is nonsense afterwards. ~~It carried the same false "opens circuitAnalysis.txt" framing~~ **fixed 2026-08-30** |
      | 11 | `todo/closed/ngspice_parser_plan.md` §7 (`:231+`) and Phase 4 (`:531-534`) | §7 is titled "extend `electricalCircuit.cpp`" and calls it "the correct reference loop"; Phase 4 still says "factor the Newton loop into a shared helper" — **the exact shape this ruling refuted.** Retargeting these is R9's business, listed here so the grep adds up |
      | 12 | `todo/closed/shared_library_and_install_plan.md:32,420` | install-set inventory |
      | 13 | `todo/README.md:365-368` | this plan's own index entry — update with R10 |

      **Deliberately NOT touched, and the implementer must not "helpfully" fix them:** `devlog/*` and
      `doc/lessons_learned_evidence.md` (INC-159, INC-444) are **dated records, kept as written**;
      `todo/closed/*` and `todo/debt_register_closed.md` are archived; and
      `.claude/worktrees/debt-register-sweep/` is **a peer's checkout — editing another worktree is
      forbidden.**

      **`circuit_usage_guide.md` §11 was stale in three ways. TWO ARE ALREADY FIXED (2026-08-30,
      Christian's approval) — only the deletion-coupled one is left for R8:**
      1. **STILL OPEN, and R8's job:** it offers the executable as the "existing FEM-free runner"
         (`:585-592`). That bullet dies with the file. Nothing to correct — it is to be *removed*
         once R7 lands, with the netlist path (already described in the next bullet) as the answer.
      2. ~~it calls the ngspice work "deferred, nothing written" and points at
         `todo/deferred/ngspice_parser_plan.md`~~ **FIXED 2026-08-30.** That was a broken path and a
         false status — v1 shipped 2026-08-27. The bullet now records the importer as shipped, names
         `NetlistParser`/`NgspiceCircuitFactory` and the hybrid `circuit → file` key, and lists
         `circuitrun` and PULSE/PWL as what remains.
      3. ~~it states the executable "calls `shift()` at the *end* of each step (output is written at
         the start, i.e. the previous step's solution)"~~ **FIXED 2026-08-30.** It was false in
         **three** ways, not one: the MNA matrix is recomputed on *every* attempt, not only when Δt
         changes (`electricalCircuit.cpp:120`); `shift()` runs at the **top** of each attempt
         (`:119`), not the end; and the solution is written **after acceptance** (`:182-183`), so
         each line is that step's own result, not the previous step's.

      **Already done, 2026-08-30 (Christian approved the edit) — do not redo:** the *factual errors*
      at sites 9, 10 and 4 (claims 2 and 3) are corrected in the tree. They were stale independently
      of DR-39 — false the moment the DR-118 fix landed on 2026-08-29, and false about the file
      direction since they were written — so they were not left to wait behind R7's gates. R8's
      remaining work at those sites is only the *deletion* churn: remove the bullet and the
      paragraph that name a binary that no longer exists.

      Then run the Codex language sweep over the touched sections of sites 2, 4 and 9 — those are
      user-facing guides, so the sweep is mandatory (CLAUDE.md, "Prose Gets a Language Sweep").
      *(The 2026-08-30 factual corrections were not swept: they replaced wrong statements with right
      ones, which is currentness work, not prose. The sweep is owed on whatever R8 rewrites.)*
      Tell it that file paths, class names, CMake option names and the header block are off limits.
      Sites 5–8 and 11–13 are generators, config and working artifacts: **not swept.**

- [x] **R9** *(after R7)* — Debt register: **strike DR-39 outright** — ID and description struck,
      status cell left unstruck and carrying the closure evidence, per the register's "How to read a
      row". **Three** things must happen, or the strike is not defensible:
      1. **(a) discharged by deletion**, with the four gates named and the fact that they RAN
         (`make check` result and date). If a gate has not run, the row does not get struck —
         "pending run is not a row a static sweep can close".
      2. **(b) re-homed, not implemented**, by Christian's ruling of 2026-08-30: `Pulse`/`PWL` is a
         feature living at `todo/closed/ngspice_parser_plan.md` Phase 6, whose user-visible behaviour today
         is the parser's explicit refusal (`cl_NgspiceCircuitFactory.cpp:430-436`,
         `.hpp:62`). Carry (b)'s three constraints ACROSS into that Phase 6 entry before striking
         — they are the row's most expensive content and they must not die with it: (i) it is not
         circuit-scoped, `boundary_condition_function_type()` feeds Maxwell and thermal BCs too, so
         it is an Input Contract change in **both** artifacts; (ii) `mValues` is a fixed 7-slot
         `Cell< real >` broadcast whole, and **PULSE's seven parameters COLLIDE with that schema
         rather than merely needing more room** — all seven slots are already named
         `BELFEM_BCVAL_AMPLITUDE … FUZZYNESS` (`cl_SourceFunction.hpp:46-52`), so "add slots" is the
         wrong mental model (wording corrected 2026-08-30 by audit); **PWL does not fit at all**,
         its variable length being exactly what the MPI rule routes through `share`/`receive`;
         (iii) PWL storage must be `Cell< real >` (`cl_SourceFunction.hpp:43-44`). Add the split
         recommendation from the closing ruling (fund PULSE, gate PWL on a real deck) as an open
         question on Phase 6, marked as Claude's recommendation at medium confidence — **not** as a
         decision.
      3. **Retarget `todo/closed/ngspice_parser_plan.md` §7 and Phase 4 in the same edit.** §7 (`:231+`) is
         titled "extend `electricalCircuit.cpp`" and calls that file "the correct reference loop";
         Phase 4 (`:531-534`) still says "factor the Newton loop out of `Controller::solve_circuit()`
         into a shared helper — resolving the top-vs-bottom `tEpsilon0` drift". **That is the shape
         Christian refuted.** Leaving them makes the parser plan point at a tombstone and re-propose
         the extraction. The FEM-free-runner *need* is real and should survive as an unfunded item;
         the *mechanism* must not.

      **Anchor trap — do NOT follow Appendix A here.** Appendix A cites the ε₀ assignments as
      `cl_FEM_Controller.cpp:2143` and `electricalCircuit.cpp:155`. Both are wrong and were verified
      wrong on 2026-08-30: `:2143` is `dump_system_if_requested( …, "magnetic", … )` in the
      **magnetic** Newton, and `:155` is the `omega *=` print path. The live assignments are
      controller **`:2547`** and demo **`:162`** (the row's existing `:2546` is off by one). R9
      **must not re-litigate the drift claim at all** — it was downgraded to structural in the row on
      2026-08-30 — and copying Appendix A's numbers would replace correct cites with magnetic-solver
      noise.
      ~~with the 2026-08-13 "behavioural drift" claim corrected as a false positive~~ and
      ~~DR-118 discharged outright~~ — **both already done ahead of this step**: the drift was
      downgraded to structural in the row on 2026-08-30 (with the stale `residual()` anchor
      corrected to `cl_ElectricalCircuit.cpp:997-1001`), and DR-118 was fixed in place and struck on
      2026-08-29, so it is not this plan's to discharge.
      Also: move the row to `todo/debt_register_closed.md`, update the `[P]` count and the
      circuit-module row list in the preamble, and run `scripts/check_doc_claims.py` — it guards the
      register's DR counts.
- [x] **R10** — Gate: **`make check` (Christian runs — not the implementer; the user runs builds).**
      `check-fast` is enough for the circuit suite alone, but R7/R8 touch `src/executables` and
      `config/globals.cmake`, so the full `check` is what proves the executable list still links.
      Confirm the four new cases were actually **compiled into** the built `test_circuit` binary and
      ran — a green suite total is not evidence that a specific case ran (the DR-117 and DR-140
      precedents). Then a `devlog/dlYYYYMMDD_circuit_demo_harvest.md` dated to the session that lands
      it, and the `devlog/README.md` entry. Record in it: whether R1′ was needed and the observed
      Newton iteration counts; whether R6 held at 1e-12; and that DR-39 is struck.
      ~~`dl20260828_…`~~ — the plan did not land on the day it was written; the 2026-08-30 ruling and
      re-scope are recorded in `devlog/dl20260830_dr39_ruling.md`.

## 5. Open Design Questions

- **O1 — Bridge load resistance.** RESOLVED 2026-08-28 → 100 Ω (Christian approved completing the
  bridge). At 100 Ω the peak diode current is 96 mA and each drop is 0.179 V, so the two-drop
  signature is ~3.6 % of the peak: large enough to gate on, small enough that the circuit is still
  recognisably a rectifier.
- **O2 — Tolerances.** ~~1e-12 for R2 and R6, 1-2 % for R3-R5~~ **REVISED 2026-08-30 now that R3 is
  struck and the discretisation is quantified.** Current proposal: **1e-12 relative for R2**
  (purely resistive, no discretisation error exists) and **for R6** (bit-comparable, see the
  arithmetic-form note in that step); **±1 % for R4 and R5**. The ±1 % is deliberately loose
  relative to the real error — at ωΔt = 0.0157 under BDF3 the R4 amplitude error is O(10⁻⁶) — so it
  gates *topology and switching*, not integration accuracy. ~~**Still open for an auditor:**~~ **RESOLVED 2026-08-30 by the audit round: keep ±1 %, do not
  tighten R4 to 0.1 %.** The auditor's reasoning, adopted: the windows open at steps 100 and 1201,
  past the BDF ramp, and discretisation error at ωΔt = 0.0157 under BDF3 is O(10⁻⁶) — so tightening
  would convert a topology/switching gate into a discretisation gate, which §2 explicitly refused.
  ±1 % still fails loudly on a missing switch (0.998 vs 1.506 is 51 %) and, for R5, on a missing
  two-diode drop (0.357 V = 3.7 % against a ±96 mV band). **This is a design call recorded from a
  single auditor voice — Christian may overrule it.**
- ~~**O3 — Do we owe the DR-123 test here?**~~ **RESOLVED 2026-08-30, by events.** The register's
  owed reproducer was "switch crossing `t_switch` inside a rejected step"; DR-123 was fixed and
  struck 2026-08-29 and the case is now a green test in the tree
  (`SwitchLatchRevertsOnRejectedStep`, `test_ElectricalCircuit.cpp:158`). Nothing owed here, and
  nothing to note in R9.

## 6. Pre-registration

Written **before** any code exists, so that a surprise is reportable as a finding rather than
absorbed as a tuning step. Protocol §11: an outcome that contradicts a line below is evidence, and
it goes in the devlog whichever way it falls.

| # | Prediction | If it holds | If it fails |
|---|---|---|---|
| P1 | **R5 converges with ω ≡ 1** — the warm start moves ≤ 0.31 V per step and no diode approaches the ≈ 18.4 V overflow bound | **R1′ is not implemented.** Record the worst-case Newton iteration count in the devlog as the evidence | R1′ lands: the ω update from `cl_FEM_Controller.cpp:2564-2574` behind a defaulted parameter, existing callers bit-for-bit unchanged. Report that the 2026-08-30 downgrade in §3.1 was wrong |
| P2 | **R6 agrees with R4 to 1e-12 over all 1600 steps** — the SPICE number path reproduces `200*1e-6`, `3.183*1e-3`, `24.975*1e-3` and `1.0/50.0` bit-for-bit | the deck↔API tie is exact and stated as such | **Do not loosen the tolerance.** Report where the traces separate and by how much; a 1-ulp component-value difference and a wrong-step switch fire look nothing alike, and the second is a defect |
| P3 | **The switch fires at the end of step 500 in both R4 and R6** (first end-of-step time ≥ 24.975 ms is 25.0 ms) | the off-boundary `t_switch` did its job | if the two runs fire on different steps despite the half-step margin, the tie analysis in §3.2 is wrong and the finding is larger than this plan |
| P4 | **R2's source-branch current is negative** and the resistor's is positive, same magnitude | the MNA convention is pinned by a test for the first time | a positive source-branch current means §3.2 misread the stamps; correct §3.2, do not flip the test to match |
| P5 | **No new test file and no `CMakeLists.txt` change are needed** | the de-duplication row does not itself duplicate a helper | if a new TU turns out unavoidable, say why in the devlog — it is the one outcome that undercuts the row's own argument |
| P7 | **R5's `nan`, if it happens, comes from commutation at the source zero crossings — not from IEEE overflow, which a 10 V source cannot reach** (added 2026-08-30 by audit) | the mitigation is a globally smaller Δt, then R1′ | if a diode genuinely exceeds 18.45 V during an iterate, the warm-start bound in §3.1 is wrong and that is the finding |
| P6 | **Deleting `electricalCircuit.cpp` breaks nothing.** Restated 2026-08-30 after the audit: the claim is true as *invocation* and was false as *documentation* — no script or deck runs the binary (`examples/scripts/Allrun:120-125` explicitly refuses to select it), but thirteen files name it | `make check` stays green, and after R8/R9 `grep -rn electricalCircuit` returns only `devlog/`, `doc/lessons_learned_evidence.md`, `todo/closed/`, `todo/debt_register_closed.md`, this plan, and the peer worktree | any *executable* consumer found is a scope change, not a fix-in-passing |

**Not pre-registered, because it is not this plan's to decide:** whether Phase 6 should build `PWL`
at all. The closing ruling records a recommendation; the decision belongs to whoever funds Phase 6.

## 7. Definition-of-Done Checklist

- [x] Every gap-table row mapped to a step or an open question.
- [x] Each claimed gap backed by a citation, not assumption.
- [x] Ordered steps with dependencies.
- [x] Open questions logged, not silently decided (O1 and O3 resolved; O2 still open for an auditor).
- [x] Predictions pre-registered before implementation (§6).
- [x] R2, R4, R5, R6 land in the **existing** `tests/circuit/test_ElectricalCircuit.cpp`.
- [x] `make check` green with the harvested cases (Christian runs the build), **and** the four new
      cases confirmed present in the built binary and reported as run — not inferred from the suite
      total.
- [x] `grep -rn electricalCircuit` returns only dated records (`devlog/`, `doc/lessons_learned_evidence.md`),
      archived rows (`todo/closed/`, `todo/debt_register_closed.md`), this plan, and the peer
      worktree. **All 13 live sites in R8 closed, generator (site 5) included.**
- [x] `src/circuit/doc/circuit_usage_guide.md` §11 corrected on all three counts, and sites 2, 4
      and 9 swept by Codex.
- [x] `scripts/update_doc_index.py:124` fixed **before** the two files it generates, so they do not
      regress on the next refresh.
- [x] DR-39 struck, archived to `todo/debt_register_closed.md`, and (b)'s three constraints carried
      across into `todo/closed/ngspice_parser_plan.md` Phase 6 **before** the strike.
- [x] `scripts/check_doc_claims.py` run after the register edit.

## 8. Audit Trail

- Exchange thread: `tmp/ai_exchange/circuit_demo_harvest.md` (Codex + Grok, plan round 2026-08-28).
- Exchange thread: `tmp/ai_exchange/review_dr39_close.md` (Codex `gpt-5.6-terra`/high + Grok
  `grok-4.6`/high, jury, blind, DR-39 closing plan round 2026-08-30). **Round complete:
  pre-registration frozen before dispatch, every citation re-verified against the tree,
  reconciliation table appended. Three P0s and six P1s CONFIRMED and applied above.** Both auditors
  stayed read-only (`src/`, `tests/`, `config/` byte-identical to the pre-round snapshot).
  **Status: REVIEWED, NOT VERIFIED — nothing was built and nothing was run.**

---

## Appendix A — Decision: delete the executable rather than extract a `CircuitSolver`

DR-39(a) has been carried since 2026-08-13 as "extract the shared Newton loop into a `CircuitSolver`
class". The premise was three copies; there are two, and the second one is a hardcoded demo with no
input file and no netlist. Extracting a class to be shared between one production caller and one demo
is work spent to keep the demo. Deleting the demo leaves `Controller::solve_circuit()` as the single
implementation and stops the release shipping a testbench as a user-facing binary — provided the
demo's *content*, the four circuits, is preserved, which is what R2/R4/R5 are for.

**Amended 2026-08-30 — one leg of this argument expired, and the conclusion is unchanged.** As
written, the appendix leaned on the demo having "a filed defect in its timestep control (DR-118)"
and on deletion discharging DR-118 with it. Neither holds: DR-118 was **fixed in place** on
2026-08-29 (`electricalCircuit.cpp:113-193`, controller parity) and struck on its own forced-rejection
gate. So the demo's loop is now *correct*, and deletion discharges only DR-39(a).

What survives is the part that was always doing the work, and it got stronger rather than weaker:

- **The binary cannot be used.** It takes no input file and no arguments beyond `-v`, and writes
  `circuitAnalysis.txt` with hardcoded node and component indices. It can solve Frederic's circuit
  and nothing else, yet it ships installed (`config/globals.cmake:14`).
- **Three of its four circuits are not merely stale, they are unbuildable.** They call
  `create_voltage_source( SourceType::SINE, … )`, and `SourceType` exists **nowhere in the tree** —
  `rg -w 'SourceType' src/` returns only those three commented lines (`:69,76,82`).
- **The FEM-free-runner need it was standing in for now has a better answer.** The netlist front-end
  shipped 2026-08-27, so a file-driven circuit runner is a real feature to be specified and funded
  (parser plan Phases 4–5), not a hardcoded demo to be kept alive. Keeping the stub as a placeholder
  for that feature ships the placeholder and not the feature.

**Also refuted, permanently:** the `CircuitSolver` extraction itself. The two copies are the same
~15-line loop, differ only in where their constants come from and which box they print, and are
provably equivalent in the one place it matters (`residual()` = `norm(mRHS)/mRHSnorm`,
`cl_ElectricalCircuit.cpp:997-1001`, no hidden state). A shared helper would buy a seam, not safety.
Christian ruled on 2026-08-30: no class, delete the demo. Do not re-propose the extraction.

**Correction to the register, high confidence.** ~~This paragraph's line numbers are wrong and R9
must not copy them — see the anchor trap in R9.~~ **STRUCK 2026-08-30 by the audit round, in place:
`cl_FEM_Controller.cpp:2143` is `dump_system_if_requested` in the magnetic Newton and
`electricalCircuit.cpp:155` is the `omega *=` print line. The live ε₀ assignments are `:2547` and
`:162`. The argument below is sound; only its anchors are not.** DR-39's row asserts "one real behavioural drift —
the controller updates `tEpsilon0_Circuit` at the TOP of the loop (`:1999`), the executable at the
BOTTOM (`:159`)". The two are the same recurrence: the controller's top-of-loop
`tEpsilon0 = tEpsilon` copies the *previous* residual before the new one is computed, which is
exactly the state the executable's bottom-of-loop assignment leaves behind; both start from
`BELFEM_REAL_MAX`. ~~Current line numbers are `cl_FEM_Controller.cpp:2143` and `electricalCircuit.cpp:155`.~~
**Both refuted 2026-08-30; the live assignments are `cl_FEM_Controller.cpp:2547` and
`electricalCircuit.cpp:162`.** The ngspice-parser round of 2026-08-27 already downgraded this to
"structural-only" (`devlog/README.md:21`) without amending the row; R9 amends it.
