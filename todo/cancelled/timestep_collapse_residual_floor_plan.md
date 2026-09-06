# Timestep Collapse on Highly Nonlinear Decks: Diagnosis and Mitigation

> **CANCELLED 2026-09-03** (todo/ currentness sweep, round 3): superseded — same as its design twin; the M1–M8 mitigations were never approved and the collapse mechanism was identified elsewhere. Status lines and checkboxes below are as they stood at closure and are not maintained.

**Date:** 2026-08-10
**Purpose:** Both `greg5` (2D tapestack, magnetic-only Picard) and `sidecoatings` (3D coated
conductor, coupled Newton) ratchet Δt down to the deck minimum and abort. The mechanism, in one
sentence: **the controller's promote/demote and stall machinery destroys more convergence progress
than the nonlinear iteration can regain, and the guards then read the resulting flat patch as a
stall and cut a timestep that was never the problem.**
**Module:** `src/fem/kernel` (`cl_FEM_Controller.{cpp,hpp}`), secondary
`src/fem/kernel/cl_FEM_DofMgr_SolverData.cpp`, `src/fem/iwg/cl_IWG_Timestep.hpp`
**AIs involved:** Claude (this analysis). Codex/Grok audit **not yet run** — see R0.
**Status:** OPEN — analysis complete against the `sidecoatings` t = 13.8762 ms log (timestep 694),
which is reconstructed **iteration by iteration** below with no unexplained events. **No source
modified, no run performed.** Every mitigation in §4 is unapproved. The two decks have *different*
dominant causes (D1–D4 for `sidecoatings`, D5 for `greg5`); do not treat them as one bug.

> **Scope guards:**
> - Companion file `todo/controller_picard_tapestack_regression.md` holds the `matfix` →
>   `sideconnectors` controller diff audit (jury round, 2026-08-10). This plan **supersedes its
>   §0.1 BDF row** (see D5) and otherwise builds on it rather than repeating it. D-numbers here are
>   independent of that file's D-numbers; cross-references are explicit.
> - Out of scope: `cl_ThinShellFactory`, the side-connector work, the material-database repair
>   (`todo/closed/example_deck_and_material_db_repair.md`) — except where R4 measures them.
> - The `sidecoatings` warm start from `.bfm` is a workaround, not a fix, and is treated as such.

---

## 1. What the t = 13.8762 ms log actually shows

Deck: `cmake-build-debug/sidecoatings/input.conf` — `algorithm : Newton`, `tolerance : 1e-6`,
`max iterations : 200`, `min iterations : 2`, `target iterations : 50`, `method : bdf1`,
`minimum timestep : 100 ns`. Δt at the failing attempt is already **1e-7 s = the floor**.

The whole 133-iteration attempt is reproduced by four rules. Every promotion, demotion, ω jump and
the final cut is accounted for; **nothing in the log is unexplained**. Confidence: high.

| # | rule | citation |
|---|---|---|
| P | promote Picard → Newton when `mEpsilon < mEpsilonSwitch` (**1e-4**) and `mIteration > 1` | `cl_FEM_Controller.cpp:748`; `hpp:112` |
| T | Newton "trust" growth: first-trial Newton that halves ε sets `tGrowth = 2.0` | `cpp:1287-1294` |
| S | stall = mean absolute deviation of the last 5 dB samples `< mStallBand` (**0.2 dB**) | `cpp:641-700`; `hpp:82-83` |
| L | after a Newton stall, `mJustPicard` latches Picard for the rest of the attempt | `cpp:685`, `:785-788` |

**Reconstruction.** ε orbits 1e-4 — *exactly* the promotion threshold — for the first 81 iterations:

| its | event | ε | consequence |
|---|---|---|---|
| 1–2 | Picard | 9.6e-5 → 9.0e-5 | below `mEpsilonSwitch` ⇒ rule P arms |
| 3 | **promote** (first flip, ω = 0.5·ω_P) | 4.1e-4 | +6 dB kick |
| 4 | demote, Picard inherits ω = 0.25 | 7.0e-3 | **+18 dB lost** |
| 5–12 | Picard grind | back to 9.0e-5 | 8 iterations to undo iteration 3 |
| 13 | **promote** again (ω = 0.537, no damping — `mFirstFlip` spent) | 2.1e-5 | good step… |
| 13→14 | **rule T**: ω = 0.537 × 2 = 1.074 → clamp **1.000** | 3.6e-5 | immediate regression |
| 15–16 | ω halves, demote | 1.1e-3 | **+17 dB lost** |
| 17–22 | Picard grind | back to 9.1e-5 | 6 iterations |
| 23–41 | **promote → kick → demote → grind**, third time | 8.2e-4 → 1.1e-4 | 19 iterations |
| 42 | Picard blows up on its own at ω = 0.494 | 1.9e-3 | Picard is only stable here for ω ≲ 0.3 |
| 53→54 | **rule T** again: ω = 0.619 × 2 → clamp **1.000** | 1.5e-5 → 3.7e-5 | second identical overshoot |
| 59, 66 | two more promote/demote cycles | | |
| 71 | **rule S** fires (MAD of −37.68…−37.76 dB = 0.079 dB) | 1.7e-4 | Picard "stalled" |
| 72 | `try_escalate_to_newton()`, ω reset to ω₀ = **1.0** | 3.0e-4 | escalation kick |
| 76–81 | escalated Newton collapses | 7.0e-3 | **+27 dB lost** |
| 82 | **rule L**: Newton stalls, Picard latched for the rest of the attempt | 6.8e-3 | *the machinery finally gets out of the way* |
| 82–127 | **pure Picard, no promotions** | 6.8e-3 → **1.2e-5** | **−27.4 dB in 45 iterations, monotone** |
| 128–133 | ω sawtooth flattens the curve | 2.0e-5 … 4.4e-5 | |
| 133 | **rule S** fires: MAD = **0.199 dB** vs `mStallBand` 0.2 | 4.4e-5 | escalation already spent ⇒ **cut Δt ⇒ abort at the floor** |

**Bottom line.** Once the controller stopped intervening (iteration 82), the nonlinear problem
converged steadily at ≈ 0.6 dB/iterate and was 1.2e-5 against a 1e-6 tolerance — roughly 25 more
Picard iterations from success, with 67 of the deck's 200 still unspent. It was killed by a
stall detector whose threshold it missed by **0.001 dB**. The preceding 81 iterations — 61 % of the
budget — were spent losing and regaining the same 20 dB, four times, to a promote/demote cycle that
triggers because the residual plateau sits on top of `mEpsilonSwitch`.

The Δt floor and the `BELFEM_ERROR` are therefore **symptoms two levels down**: the guards cut a
step that was converging, the PID hold (`cpp:1961-1970`) slows recovery ~3×, and 24 halvings later
the floor abort delivers a diagnosis ("the residual is not responding to Δt") that is *not what
happened here*.

---

## 1.1 Corroborating trace — timestep 922, t = 19.7633 ms (added 2026-08-10, same day)

Christian supplied a second `sidecoatings` log: the warm-started run at Δt = 0.0389 ms, **accepted**
after 88 iterations / 24 s — the same disease, survived. It independently re-confirms three
defects, adds one new one, and delivers the first *measured* residual floor.

| its | event | matches |
|---|---|---|
| 3–4 | promote at ε = 2.3e-5 < 1e-4, Newton kicks +20 dB to 2.1e-3, ω collapses 0.25 → 0.035 | **D1**, third independent occurrence |
| 5–9 | demoted Picard grinds at ω ≈ 0.02, −0.08 dB/it; stall guard fires at it 9 (MAD ≈ 0.10 < 0.2) | D2's guard, here firing on a *genuinely* slow phase |
| 10–12 | escalation enters at ω₀ = 1.0; two good steps, then **+35 dB blow-up** to 3.7e-2 — the largest kick in either trace. The thermal kernel is kicked +18 dB as collateral (it 11) | **D4**, now measured; coupled collateral is new |
| 12–17 | Newton flat at −14.x dB; stall guard demotes with the `mJustPicard` latch | the it-82 rescue mechanism, again |
| 17–36 | latched Picard descends −15 → −49.3 dB | Picard-alone converges, again |
| 37, 41, 45, 47, 51 | **repeated +10…+30 dB blow-ups each time the growth rule pushes ω back above ≈ 0.2–0.43**, each followed by α-halving and re-descent | **D9 — new, see below** |
| 58–87 | magnetic reaches −66 dB (3e-7, converged); thermal descends to **6.36e-6 and floors** (−51.97 dB, bit-flat) against `tolerance : 1e-6` | first measured floor |
| 88 | `mThermalStalled` latches after 5 bit-flat iterates → warning → `run_coupled()` exits (`cpp:2347-2359`) → **step accepted** | the in-tree accept-on-floor escape, thermal-only |

**What is new:**

1. **A genuine residual floor above tolerance now exists as a measurement — but it is the
   *thermal* kernel** (6.36e-6 vs 1e-6), with max T = 77.15 K (essentially isothermal — the
   relative thermal residual at near-zero thermal load is plausibly conditioning-dominated).
   The magnetic kernel demonstrably reaches −66 dB at this state. This *sharpens O1*: the
   floor hypothesis and the convergence hypothesis are **both right, for different kernels**.
   R3's Δt-response measurement must be run per-kernel.
2. **The thermal-stall accept escape (`cpp:1156-1159`, `:2347-2359`) is the existing precedent
   for M9-class acceptance** — magnetic has only the inert opt-in absolute tolerance. The
   asymmetry is now demonstrated load-bearing: this step lived because of it.
3. Deck note: `nonlinear thermal { min relaxation : 0.01 }` pinned ω₂ at 0.01 for ~10 iterations
   while thermal chased the moving magnetic baseline (its 22–35) — deck-side aggravator.

### D9 — The ω growth rule has no memory of the discovered stability ceiling. **P1. Confidence: high.**

- [ ] Resolved
- Both adaptation tails (`cpp:1305-1330` coupled, `:1508-1530` magnetic) grow ω by the arctan rule
  on every improvement and halve on regression — but nothing records *where* regression happened.
  In this trace Picard blew up at ω ≈ 0.43, 0.39, 0.33, 0.21, 0.21 within one attempt (and at
  0.494 in the timestep-694 trace): after each α-collapse the growth rule climbed straight back
  into the same ceiling, costing a +10…+30 dB kick each time — ≥ 5 kicks in 88 iterations. The
  it-128–133 sawtooth that finally killed timestep 694 is the same pattern at small amplitude.
- Mitigation candidate (for the design-doc revision, not yet approved): an attempt-local ω
  ceiling — on a regression exceeding some dB threshold, record `mOmegaCeiling = fraction ×`
  the ω that caused it; the growth rule clamps to the ceiling for the rest of the attempt
  (decaying or resetting per timestep). This is AIMD with memory, and it also gives MIT-3b/O2
  the "largest safe ω" statistic the auditors asked for.

---

## 1.2 The restart-passes-wall series (added 2026-08-27 night)

Four times on 2026-08-27, tapestack3d_coarse walls that Δt reduction could not cross were
crossed by a warm restart — the fourth time with **identical deck config and identical binary**
(gauged Picard + thermal MUMPS: in-run cascade at t = 4250.88 ms burned 15 rejections down to
Δt = 2 µs without moving; the restart from the last save cleared the same t in its first
attempts and stood 18 ms past it with one rejection, `out2.txt`). A restart rewinds five state
machines at once; the two with mechanistic backing here:

1. **BDF5 step-ratio history poisoning (D5's consequence 1, now with a live candidate case).**
   The memdump holds the history of a converged cruising step (healthy ratios); the cascade
   fills `mH` with ratios like 1 ms : 2 µs, and the variable-step BDF5 coefficients then cancel
   catastrophically — a residual floor that does not respond to Δt, self-amplified by every
   further halving. The restart rewinds to pre-poisoned history.
2. **Controller-state accumulation** (collapsed ω, latches, watchdog counters, PID hold) — the
   restart is M9 performed by hand, plus the history flush M9 lacks.

**Discriminator (one deck word):** `scheme : bdf1` — no multi-step ratio history to poison. If
walls become restart-INSENSITIVE under BDF1, mechanism 1 dominates. Ties into R1 and M8.

- [ ] **M10 — Automatic in-run restart (upgrade of M9).** After N consecutive rejections:
  flush the BDF history to a fresh ramp, clear the Anderson history, reset ω and the guard
  latches, retry at the initial-Δt cap — i.e., exactly what a manual warm restart does, without
  the operator. M9's "one conservative attempt" survives as the stage after M10 fails. Needs
  the D5/M8 coefficient verification first, or M10 will mask it.

---

## 2. Defect Tracker

### D1 — Promote/demote chatter when the residual plateau sits on `mEpsilonSwitch`. **P0. Confidence: high.**

- [ ] Resolved
- `cpp:748` promotes on a **single** sample crossing a **fixed absolute** threshold (`hpp:112`,
  default 1e-4) with no hysteresis, no minimum dwell, and no memory of the previous promotion's
  outcome within the same timestep attempt.
- In the log this fires at iterations 3, 13, 23, 53, 59, 66 — every time ε dips under 1e-4 — and
  each promotion cost 6–18 dB and 6–19 recovery iterations. Iteration 63 (ε = 1.00e-4, printed
  `0.000100`) demoted and iteration 65 (9.6e-5) re-promoted: the solver is switching algorithms on
  the fourth significant digit.
- `mFirstFlip` damps only the **first** promotion of an attempt (`cpp:762`, deliberate per the
  Codex+Grok audit note). Promotions 2..n enter at full inherited ω, which is why 13 and 53 were
  the worst.
- **NOT deck-recoverable at zero cost** — `tolerance switch` exists (`cpp:2504-2506`) but any fixed
  value can be landed on by some other step's plateau. Lowering it to e.g. 1e-8 disables Newton
  entirely, which is a legitimate *experiment* (R2) but not a fix.

### D2 — `mStallBand = 0.2 dB` kills slow-but-real convergence. **P0. Confidence: high.**

- [ ] Resolved
- `hpp:83`; guard `cpp:641-700`; call site `cpp:1257-1261`. On `matfix` the band was 0.001 dB.
- **Measured, not inferred:** at iteration 133 the last five dB samples are
  −42.84, −43.07, −43.18, −43.36, −43.53; mean −43.196; MAD = **0.1993 dB** against a 0.200 dB
  threshold. The run died on a 0.35 % margin.
- The band converts any convergence slower than ≈ 3.8 %/iterate into an abandoned timestep
  (arithmetic in `controller_picard_tapestack_regression.md` D1). A Picard iteration on a
  power-law conductor with n = 25 at ω ≈ 0.1 is *routinely* slower than that and still converging.
- A mean-absolute-deviation test cannot distinguish "flat" from "slowly and monotonically
  descending": both have small MAD. The correct discriminator is the **trend** (sign and magnitude
  of a fitted slope), not the scatter. This is the design defect; 0.2 vs 0.001 is only its gain.
- **Deck-recoverable:** `nonlinear { stall tolerance : 0.001 ; }`.

### D3 — Newton "trust" growth doubles ω into the clamp and immediately regresses. **P1. Confidence: high.**

- [ ] Resolved
- `cpp:1287-1294`: a first-trial Newton step that at least halves the residual sets `tGrowth = 2.0`
  instead of the ≤ 1.3 arctan rule.
- Both occurrences in the log behave identically: ω 0.537 → 1.074 → clamped 1.000 (it 14, ε worsens
  2.1e-5 → 3.6e-5) and ω 0.619 → 1.238 → clamped 1.000 (it 54, ε worsens 1.5e-5 → 3.7e-5). **2 of 2
  trust jumps ended in immediate regression**, each followed by an α-halving cascade and a demotion.
- The rule was tuned on `ts34`, where the complaint was an endgame crawl at ω ≈ 0.42. It generalises
  badly to a residual plateau: one good Newton step near a floor is not evidence of a contractive
  tangent, it is a coin flip (which is exactly the reasoning already written down for
  `mOmegaNoiseBand`, `hpp:104-110`).
- **NOT deck-recoverable.**

### D4 — The escalation resets ω to ω₀ = 1.0 on a state known to be unstable above ω ≈ 0.3. **P1. Confidence: high.**

- [ ] Resolved
- `cpp:561` (`mOmegaNewton = mOmega0`) with the rationale that a collapsed Picard ω would freeze
  the Newton step. But the escalation fires *precisely* when Picard has stalled, i.e. when the
  controller has just spent 70 iterations discovering the stable ω range.
- In the log, iteration 42 shows Picard itself diverging at ω = 0.494; the escalation at iteration
  72 then re-entered at ω = 1.0 and cost +27 dB over iterations 72–81 before rule L latched Picard.
- Suggested direction (unapproved): enter the escalated Newton at `min(ω₀, 2·ω_picard_observed)` or
  at the largest ω that produced an accepted step this attempt.

### D5 — `greg5` silently switched from BDF1 to **BDF5** across the branch. **P0 for `greg5`. Confidence: high (code), untested (effect).**

- [ ] Resolved
- **This supersedes `controller_picard_tapestack_regression.md` §0.1, row `mTimeStepping`**, which
  records "No difference". That row was written against a deck variant reading `scheme : bdf1`.
- The deck actually in the tree, `cmake-build-debug/greg5/input.conf:39`, reads **`scheme : bdf5`**.
- `matfix` parsed **only** `method` (`git show matfix:…cl_FEM_Controller.cpp`, line 2146), so
  `scheme : bdf5` was **ignored** and the run was BDF1.
- `sideconnectors` prefers `scheme` over `method` (`cpp:2750-2756`), so the *same deck* now runs
  **BDF5**.
- Consequences to test, in order of likelihood:
  1. After repeated Δt halvings the step-size history `mH` holds ratios like
     1e-7 : 1e-4 — the variable-step BDF5 coefficients (`compute_bdf_coeffs_5`) then involve large
     cancelling terms, which is a textbook route to a residual floor that *does not respond to Δt*
     (i.e. the exact wording of the abort message).
  2. The growth clamp is 1.2 for `order_active() >= 4` (`cpp:1980-1983`), but the classical
     zero-stability step-ratio bound for BDF5 is ≈ 1.13 — the clamp is above it.
  3. BDF5 is only A(51.8°)-stable. Pure diffusion eigenvalues are real-negative and safe, but the
     HTS power-law tangent is not guaranteed to stay there.
- **Deck-recoverable and near-free to test:** change one word to `scheme : bdf1` (R1).
- Separately: **the legacy key silently changing meaning is itself the defect.** Any deck in the
  wild carrying an unread `scheme :` line changed integration order on this branch with no log line.

### D6 — The imposed current in `sidecoatings` is ADC-quantised noise differentiated by the solver. **P1. Confidence: high (data), medium (as a contributor to this failure).**

- [ ] Resolved
- `cmake-build-debug/sidecoatings/lib/current.cpp:90-97` **linearly interpolates** a measured table,
  `I_vs_t_regular_smooth.txt` (6606 rows, uniform 53 µs spacing).
- The waveform ramps to 137 A by ≈ 2.5 ms at 77 585 A/s, then **plateaus for the remaining 350 ms**.
  On the plateau the samples are quantised to a **0.0189412 A** ADC step and dither by ±1–2 counts.
- A piecewise-linear interpolant of a dithering staircase has a derivative that is a **random square
  wave: dI/dt = 0, ±357, ±715 A/s, changing sign at essentially every 53 µs node.** Around 13.8 ms
  the tabulated dI/dt is −357 A/s while the physical current is flat.
- Why the failure appears at ≈ 13.8 ms and not at 3 ms: during the ramp the dither is 0.5 % of
  dI/dt; on the plateau it is **100 % of it**. The solver is being asked to track sign-alternating
  E-fields through an n-value power law, which is the worst possible input for a lagged-conductivity
  Picard map.
- Quantitative correspondence to test: relative dither = 0.0189412 / 137.06 = **1.38e-4**, against
  an observed residual plateau of 1e-5…1e-4 and an output voltage (`iv_results.csv`) fluctuating
  1.42e-5 … 1.98e-5 V (±25 %) between saves. Suggestive, **not yet proven causal** — R4 decides.
- This is a *user-deck* defect, but BELFEM should not be silently defenceless against it (see M6).

### D7 — A timestep cut never says why. **P2. Confidence: high. Blocks every investigation below.**

- [ ] Resolved
- `reset_timestep()` (`cpp:1836-1890`) prints only the closing box rule. Of the six paths that reach
  it, only two announce themselves (`watchdog_magnetic` `cpp:600-606`, and the escalation
  `cpp:573-579`). The stall cut, the divergence cut, the `max iterations` cut, the backtracking-budget
  cut and the solver-failure cut are **indistinguishable in the log**.
- Reconstructing §1 required reading five guards against a hand-transcribed console dump. That cost
  is paid again on every future report.

### D8 — Consequences of D1–D4 are amplified by the recovery policy. **P1. Confidence: high.**

- [ ] Resolved
- Already documented as D4 in `controller_picard_tapestack_regression.md` (PID + `mPostFailureHold`
  ⇒ ≈ 3 accepted steps of zero growth after every cut, versus ≈ 2 growth steps on `matfix`).
- Recorded here only to name the product: **cut rate × recovery cost = the ratchet.** D1–D4 raise
  the numerator; the PID hold lowers the denominator. Neither alone reaches the Δt floor.

### Non-findings — checked, do not re-audit

- **The Δt floor abort itself is correct behaviour** (`cpp:1847-1857`). Grinding below the deck
  minimum would be worse. The message is accurate about *what* it observed and misleading about
  *why*; fix the diagnosis (M5), not the abort.
- **Picard reports the pre-update residual by design** (`cl_FEM_DofMgr_SolverData.cpp:2286-2290`),
  and this is correct — see D10 (retracted) in the companion file. It does mean every guard acts on
  a one-iteration-lagged signal.
- **`mAbsoluteEpsilonTarget` defaults to 0.0** (`hpp:118`) so the absolute-tolerance escape the
  abort message advertises is inert unless the deck sets it. The key is real
  (`cpp:2494-2501`, `absolute tolerance`).
- **Anderson is off in both decks** and all its helpers early-return.

---

## 3. Investigation — ordered, cheapest first

### R0 — Cross-review this analysis before acting on it *(no dependencies)*

- [ ] Run `/cross-review` on §1–§2. The §1 reconstruction is the load-bearing claim and it was
      produced by a single voice from a console transcript; the promote/demote arithmetic and the
      0.1993 dB MAD both need an independent recomputation.
- [ ] Specifically ask the auditors to attack: (a) the claim that iterations 82–127 were genuinely
      converging rather than drifting toward a floor; (b) whether D6's dither can produce a
      Δt-independent *relative* residual at all, or only a Δt-proportional one.

### R1 — `greg5`: one-word BDF test *(no dependencies; ~1 run)*

- [ ] Run `greg5` unchanged, capture the Δt trace.
- [ ] Change `scheme : bdf5` → `scheme : bdf1` (**only that**) and re-run.
- **Decision rule:** recovers ⇒ D5 is the `greg5` regression and the companion file's §0.1 must be
  corrected before anyone tunes a controller default for that deck. Does not recover ⇒ `greg5` is
  the controller story after all and folds into R2.

### R2 — `sidecoatings`: the controller-off experiment *(no dependencies; ~1 run, no recompile)*

The single most informative run available. Every knob below is an existing input key.

- [ ] Warm-start from `tapestack3d.bfm` at 13.5 ms and re-run the failing window with:

      ```
      nonlinear {
          tolerance         : 1e-6 ;
          max iterations    : 400 ;      // the budget was never the binding constraint
          min iterations    : 2 ;
          target iterations : 50 ;
          algorithm         : Newton ;
          tolerance switch  : 1e-8 ;     // D1: never promote — isolates the Picard path
          stall tolerance   : 0.001 ;    // D2: restore the matfix band
          stall window      : 20 ;       // D2: a longer window is a weaker MAD test
          watchdog window   : 0 ;        // companion D2: watchdog off
          max relaxation    : 0.30 ;     // D3/D4: cap omega below the observed instability
          min relaxation    : 0.02 ;
      }
      ```

- **Prediction, stated before the run so it can be wrong:** the step converges in 60–120 Picard
  iterations at Δt = 1e-4 s, with no cut. If it does, D1–D4 are the whole `sidecoatings` story and
  §4 M1–M4 are the fix.
- [ ] Then relax the knobs one at a time (`tolerance switch` back to 1e-4 first, then
      `max relaxation` back to 1.0) to rank D1 against D3/D4.

### R3 — Measure the Δt response the abort message asserts *(after: R2)*

The abort claims the residual does not respond to Δt. **Nobody has measured this.** For BDF1 with a
converged history, the first-iteration residual must scale **linearly in Δt**; anything flatter is a
real, separate defect.

- [ ] From the same warm start, run the failing step at fixed Δt = 1e-4, 1e-5, 1e-6, 1e-7 s with
      `adapt timestep : false` and record ε at iteration 1.
- [ ] Fit slope of log ε₁ vs log Δt. **1.0 = healthy** (the message is a misdiagnosis, close it via
      M5). **0.0 = a Δt-independent residual term exists** — then, and only then, hunt it:
- [ ] If the slope is flat, re-enable `write_residuals_to_mesh()`
      (`cl_FEM_DofMgr_SolverData.cpp:2440`, commented out; implementation
      `cl_FEM_DofManagerBase.cpp:185`) behind an env gate and view r in ParaView. The residual's
      *spatial* home names the culprit: air φ rows, thin-shell layer rows, edge-coating rows, cut λ
      rows, or the bearing constraint.

### R4 — Is the current dither driving it? *(after: R2; independent of R3)*

- [ ] Median- or Savitzky-Golay-filter `I_vs_t_regular_smooth.txt` (the plateau is physically flat;
      the 0.0189 A steps are instrument quantisation) and re-run the failing window.
- [ ] Alternatively swap `MyCurrent` for a constant 137 A over the same window — a cleaner
      falsification.
- **Decision rule:** recovers with the same controller settings ⇒ D6 is causal and M6 is justified.
  Recovers only *together* with R2's knobs ⇒ D6 is an aggravator; record it and move on.

### R5 — Instrument the cut reason *(no dependencies, but do it before any further runs)*

- [ ] Give `reset_timestep()` a reason argument and print it inside the box, e.g.
      `Δt cut: Picard stall ( MAD 0.199 dB < 0.200 ) at eps 4.4e-5, iteration 133`.
      Cheap, once-per-cut (no hot path), and it retires D7 permanently.
- [ ] Also print, on each cut, the attempt summary: `eps_first`, `eps_best`, iteration of the best,
      and the number of promote/demote flips. R3's slope then falls out of a normal run for free.

---

## 4. Mitigation — proposals, none approved

Split deliberately: M1–M2 are **defaults**, M3–M6 are **design**. The companion file's O2 question
(should the retuned guards be conditioned on `mAlgorithm == NewtonRaphson`?) is *not* answered by
this log — `sidecoatings` **is** a Newton deck and still died on them, so algorithm-conditioning
alone would not have saved it. **That is a new argument against option (a) in that file's O2.**

- [ ] **M1 (D2) — Replace the MAD stall test with a trend test.** Fit a least-squares slope over the
      window; declare a stall only when the slope is non-negative *and* the scatter is small. A
      monotone −0.6 dB/iterate descent then survives any band setting. Keep `stall tolerance` as the
      scatter gate so existing decks still parse.
- [ ] **M2 (D2) — Until M1 lands, restore `mStallBand` to a value that cannot kill a converging
      run,** and require the window to be long enough for the MAD to mean something
      (`mStallWindow = 5` is very few samples for a 0.2 dB test).
- [ ] **M3 (D1) — Hysteresis and dwell on the Picard↔Newton handoff.** Promote on
      `ε < mEpsilonSwitch`, demote only on `ε > 10·mEpsilonSwitch`, enforce a minimum dwell of
      2–3 iterations in each algorithm, and **do not re-promote within one attempt after a
      promotion that regressed** — track the outcome, not just the count. This alone would have
      removed iterations 3–81 from the failing log.
- [ ] **M4 (D3/D4) — Retire the ×2 trust jump near a residual plateau** (require two consecutive
      trusted Newton steps, or cap growth at the arctan rule whenever the accepted ε is within a
      decade of the best-ever ε for this attempt), and enter the escalated Newton from the largest ω
      that produced an accepted step this attempt rather than from ω₀.
- [ ] **M5 (D7) — Make the Δt-floor abort report evidence, not a hypothesis.** With R5's
      bookkeeping the controller can state the measured ε₁-vs-Δt slope across the last cuts and say
      either "the residual scales with Δt — this is a nonlinear-iteration failure, not a step-size
      one" or "the residual is flat in Δt — the floor is real". The current text asserts the second
      unconditionally.
- [ ] **M6 (D6) — Defend against non-smooth sources.** At minimum, warn when a user-defined or
      tabulated source has sign-alternating first differences on the scale of the timestep. Better:
      offer monotone C¹ (PCHIP) interpolation for tabulated sources instead of piecewise-linear, so
      dI/dt is at least continuous. Document that measured waveforms must be filtered before use.
- [ ] **M7 (D5) — Never let a previously-ignored key change physics silently.** Log the resolved
      integration scheme and its source key at startup, and emit a loud notice when `scheme` is
      honoured on a deck that also carries `method`, or when the order is > 2 with
      `adapt timestep : true`.
- [ ] **M8 (D5) — Tighten the step-growth clamp for BDF5** from 1.2 to ≤ 1.13 (`cpp:1980-1983`),
      and consider refusing / warning on BDF ≥ 4 with adaptive stepping until the variable-step
      coefficients are verified against a manufactured solution
      (`todo/closed/bdf_nonlinear_mass_verification.md` may already cover part of this — check before
      duplicating).
- [ ] **M9 (safety net) — A last-resort mode before the floor.** When Δt has been cut N times in a
      row, stop tuning and switch to a fixed conservative configuration (Picard, fixed ω, guards
      disabled, large iteration budget) for one attempt. If *that* fails, the abort is honest. This
      is what iteration 82 did by accident, and it worked.

---

## 5. Open Questions

**O1 — Was the iteration 82–127 descent genuine convergence or an approach to a floor?**
*Update 2026-08-10 (§1.1):* at t = 19.76 ms the magnetic kernel reached −66 dB while the *thermal*
kernel floored at 6.36e-6 > 1e-6. Both hypotheses hold, for different kernels — R2/R3 must
separate them; a magnetic-only floor claim is now the less likely reading.
−27.4 dB over 45 iterations is monotone and fast enough to look genuine, but the last six samples
flatten at 1.2e-5…4.4e-5, which is also what a floor looks like from above. **R2 settles it and
nothing else does.** Everything in §4 is contingent on the answer: if it was a floor, M1–M4 buy
iterations but not convergence, and the real work is D6 / R3.

**O2 — Should `mEpsilonSwitch` be absolute at all?** A fixed 1e-4 threshold is a fixed point for
chatter whenever a deck's residual plateau lands on it. A relative criterion (e.g. promote when the
Picard contraction ratio degrades past a bound) has no such preferred value. Larger change; flag it,
do not start it. Related: `todo/closed/nonlinear_iteration_strategy_near_quench.md`.

**O3 — Do the two decks share a cause at all?** This plan asserts they do not: `greg5` is D5
(silent BDF1 → BDF5), `sidecoatings` is D1–D4 (+ D6). R1 and R2 are independent and can run in
parallel; **do not let one result be read as evidence about the other deck.**

**O4 — What does `matfix` do on this exact `sidecoatings` step?** The companion file establishes the
controller diff is large. A `matfix` run of the same warm start would say whether `sidecoatings`
ever worked here, or whether this deck is newly reaching a state that was previously unreachable
(the side-connector and edge-coating work is on this deck's path).

---

## 6. Definition of Done

- [ ] R0 cross-review complete; §1 reconstruction independently confirmed or corrected.
- [ ] R1 run: `greg5` BDF1-vs-BDF5 verdict recorded; companion file §0.1 corrected either way.
- [ ] R2 run: O1 answered.
- [ ] R3 run (or explicitly closed by R2): the ε₁-vs-Δt slope is a **measured number** in this file.
- [ ] R4 run: D6 promoted to causal or demoted to aggravator.
- [ ] R5 landed: no cut in any future log is anonymous.
- [ ] Agreed subset of M1–M9 implemented, with `doc/input_file_reference.md` synced in the same turn
      for every new or changed key, and a Codex prose pass on the touched sections.
- [ ] Regression gate: `sidecoatings` passes 13.8 ms **from a cold start**; `greg5` completes 15 s;
      a quench case is no worse than today; `make check-fast`.
- [ ] Devlog written; this file moved to `todo/closed/` and `todo/README.md` updated.

---

## 7. Audit Trail

- Primary evidence: the `sidecoatings` timestep-694 console log supplied by Christian, 2026-08-10.
  **It is not in the repository** — the run's stdout was not captured to a file. R5 exists partly so
  this cannot happen again; until then, keep the transcript with this file if it is needed for a
  re-audit.
- Companion: `todo/controller_picard_tapestack_regression.md` (jury round, 2026-08-10). This plan
  supersedes its §0.1 BDF row (D5) and adds a counter-argument to its O2.
- Data inspected directly: `cmake-build-debug/sidecoatings/{input.conf,iv_results.csv}`,
  `.../lib/{current.cpp,defect.cpp,I_vs_t_regular_smooth.txt}`, `cmake-build-debug/greg5/input.conf`.
- **Nothing in §1 is a measurement of a run performed by this session.** The reconstruction is
  arithmetic over a supplied log plus source reading; the confidence levels reflect that.
