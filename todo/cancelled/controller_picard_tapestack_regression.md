# FEM Controller: Pure-Picard Decks Regressed by Quench-Tuned Cut Defaults

> **CANCELLED 2026-09-03** (todo/ currentness sweep, round 3): overtaken — the reproducer deck was never run, and the controller has since been rebuilt around it (PID timestep, Anderson, watchdog and penalty opt-in); the analysis no longer describes the tree. Status lines and checkboxes below are as they stood at closure and are not maintained.

**Date:** 2026-08-10
**Purpose:** Gregory's 2D tapestack case (magnetic-only, `algorithm : Picard`) ran fine on
`matfix` and "struggles a lot" on `sideconnectors`, while quench cases improved. Mechanism in
one sentence: quench-tuned timestep-cut guards now assume a Newton rescue path that pure-Picard
decks cannot take, so stall detection becomes a timestep cut instead of algorithm demotion.
**Module:** `src/fem/kernel` (`cl_FEM_Controller.{cpp,hpp}`, secondary: `cl_FEM_DofMgr_SolverData.cpp`)
**AIs involved:** Claude (pre-registration + verification), Codex (jury auditor), Grok (jury auditor)
**Status:** OPEN — jury review complete (2026-08-10, `--jury`); all citations verified; **no source
modified**. Blocked on R1: the deck has not been run. Rankings are control-flow certainty,
**not measured proof**. The D3 severity disagreement is routed to Christian because code reading
cannot resolve it.
**Currentness sweep 2026-08-11:** findings unchanged and nothing here is closed — the anchors were
re-baselined against the current tree (`try_escalate_to_newton` `:542` → `:600`, its
`mAlgorithm == NewtonRaphson` gate at `:612`, the deck parse `:2518` → `:2935`, and the
retracted-D10 residual site `cl_FEM_DofMgr_SolverData.cpp:2286` → `:2788`). The D10 retraction
was re-read at the new anchor and still holds: the comment there states outright that
`mFieldValues` is deliberately not refreshed so `residual()` reports the pre-update force
residual.

> **Scope guards:**
> - This is a **controller** investigation. The `matfix` → `sideconnectors` diff also touches
>   `cl_ThinShellFactory` (+900 lines), `cl_FEM_DofMgr_BlockData`, `cl_MeshChecker` and others.
>   Those are **out of scope** unless R1 fails to reproduce, in which case the assumption that the
>   regression is controller-local must itself be re-tested (see O3).
> - **No fixes have been applied and none are approved.** The review round explicitly excluded fixes.
> - In scope: every control-flow path that is *live* for a magnetic-only, pure-Picard deck.

---

## 0. Reproducer

Mesh: `cmake-build-debug/greg4/2D_tapestack.msh` (also present in `greg3/`).
Deck as supplied by Christian (from Gregory), verbatim in the relevant parts:

```
solver {
    linear    { library : strumpack ; }
    nonlinear { tolerance : 1e-7 ; max iterations : 100 ; max relaxation : 1.0 ;
                target iterations : 50 ; algorithm : Picard ; }
    timestep  { initial timestep : 0.01 s ; maximum timestep : 0.2 s ;
                minimum timestep : 1 ns ; simulation time : 15 s ;
                adapt timestep : true ; save every : 0.2 s ;
                anderson stabilization : false ; scheme : bdf1 ; }
}
materials { ybco { builtin : ybco ; jc : 4e10 A/m^2 ; n : 25 ; ec : 1e-4 V/m ;
                   resistivity type : power-law ; } }
layers : tape { ybco : 0.25 mum ;  (x4) }
homology { algorithm : generalized pellikka ; }
topology { thinshell : tape { sidesets : 5,6 ; }  air { blocks : 1:6 ; }
           curves { 1 : [5,6] ; } }
boundary conditions { current { input curve : [1] ; type : ramp ; amplitude : 160 A ;
                                period : 10 s ; offset : 0.0 s ; }
                      bearing { nodes : 1 ; } }
```

Branches: good = `matfix` (`2371339c`), bad = `sideconnectors` (`1d6ef305` + uncommitted at
the time of writing; that tree is now committed — `bc578b5e`, verified clean 2026-08-10).
Merge base `429913dc`. Controller diff: 869 insertions / 41 deletions.

### 0.1 How this deck resolves — verify before trusting anything below

| quantity | value | consequence | citation |
|---|---|---|---|
| `mAlgorithm` | `Picard` | **`try_escalate_to_newton()` can never fire** — it requires `mAlgorithm == NewtonRaphson` | `cl_FEM_Controller.cpp:600`, gate at `:612`; parsed `:2935-2957` |
| `mKernel2` | `nullptr` | driver is `run_coupled()` / `iterate_coupled()`; every thermal branch dead | `src/executables/hphirun.cpp:104-106` |
| `mAndersonDepth` / `2` | 0 | all `anderson_*` helpers early-return; the Anderson half of the diff is **inert** | `cl_FEM_Controller.cpp:2830`, `:474-475`, `:494-518` |
| `mRelativeEpsilonTarget` | 1e-7 | watchdog skip-threshold `10·tol` = 1e-6 | — |
| `mIterationTarget` | 50 | the deck *expects* ~50-iteration steps | — |
| ~~`mTimeStepping`~~ | ~~BDF1~~ | ~~on `matfix` the `scheme` key was unread (`method` only) → BDF1 by default there too. **No difference.**~~ **SUPERSEDED 2026-08-10 — see below.** | `cl_FEM_Controller.cpp:2750-2756` |

> **CORRECTION (2026-08-10, Claude). The BDF row above is wrong for the deck actually in the tree.**
> `cmake-build-debug/greg5/input.conf:39` reads **`scheme : bdf5`**, not `bdf1`. `matfix` parsed only
> `method` (`git show matfix:src/fem/kernel/cl_FEM_Controller.cpp`, line 2146), so that line was
> **ignored** there and the run was BDF1; `sideconnectors` prefers `scheme` (`cpp:2750-2756`), so the
> same deck now runs **BDF5** with adaptive stepping. This is a candidate *dominant* cause of the
> `greg5` regression that is entirely independent of the controller, and it is testable by changing
> one word. It is tracked as **D5 / R1 in `todo/timestep_collapse_residual_floor_plan.md`** — run
> that before tuning any controller default for this deck. The §0 deck transcript above was taken
> from a `bdf1` variant Christian supplied, which is why the round missed it. Confidence: high.

**Bottom line:** Anderson, thermal and Newton paths are dead for this deck. The live surface is the
ω-adaptation and timestep-cut policy, plus the §1.1 asymmetry that makes those guards harsher for
Picard than for the quench decks they were tuned on.

---

## 1. Current Behaviour and How It Fails

| Failure | Mechanism | Evidence |
|---|---|---|
| Δt cut on ordinary slow Picard convergence | stall band widened 200× | `cl_FEM_Controller.hpp:83`; `cpp:695-698` → `:1257-1261` |
| Δt cut before the deck's own iteration target | new watchdog, window 30 < target 50 | `hpp:91`; `cpp:586-612`, `:1264-1269` |
| Slow residual growth never damped | Eq. 14 rule suspended for 0–5 % regressions | `cpp:1305-1312` vs `messe2023.txt:574-575` |
| Δt cannot recover between cuts | PID + post-failure growth hold | `cpp:1937-1971`, `hpp:250` |
| Run aborts instead of grinding | hard `BELFEM_ERROR` at the Δt floor | `cpp:1847-1857` |

### 1.1 Cross-cutting finding — the algorithm asymmetry (both auditors, independently)

Every new stall guard is written on the assumption that detection leads to a *rescue*: stall →
demote Newton to Picard, or escalate Picard to Newton, and keep the timestep. That path is
unconditionally unavailable when `mAlgorithm == Picard`, because `try_escalate_to_newton()` bails at
`cl_FEM_Controller.cpp:600` ( the `mAlgorithm == NewtonRaphson` gate at `:612` ). `magnetic_stagnation_forces_reset()` then returns
`! try_escalate_to_newton()` = `true` (`cpp:695-698`) and the caller cuts Δt (`cpp:1257-1261`).

**The same retuned default is a soft demotion for a hybrid quench deck and a hard timestep cut for a
pure-Picard deck.** This single asymmetry explains why quench improved while the tapestack regressed,
and it is the most likely place a fix belongs. Confidence: high.

---

## 2. Defect Tracker

Findings from the 2026-08-10 jury round: Claude pre-registered, Codex and Grok audited blind in
parallel, and Claude re-verified every citation before inclusion.

### D1 — `mStallBand` widened from 0.001 to 0.2 dB. **P0. Confidence: high.**

- [ ] Resolved
- **Raised by:** Claude (pre-reg C1), Codex (#1), Grok (P0b) — independently, 3/3.
- `cl_FEM_Controller.hpp:83`; on `matfix` `real mStallBand = 0.001 ; // dB`.
- Guard: `magnetic_stagnation_forces_reset()`, `cl_FEM_Controller.cpp:641-700`; window is 5 dB
  samples pushed at `cpp:1212`; call site `cpp:1257-1261`.
- **Arithmetic:** for constant contraction ratio ρ the 5 samples are arithmetic with spacing
  `d = −10·log10(ρ)` dB and MAD exactly `1.2·d`. Fires when `1.2·d < mStallBand`:

  | band | fires below | i.e. contraction slower than |
  |---|---|---|
  | 0.001 dB (`matfix`) | d < 8.3e-4 dB | ρ > 0.99981 — effectively never |
  | **0.2 dB (`sideconnectors`)** | **d < 0.167 dB** | **ρ > 0.962, i.e. < 3.8 %/iterate** |

- The header rationale (`hpp:76-81`) is drawn from the *sidecoatings* case, which has a residual
  floor **above** tolerance. The tapestack has no such floor, so the widened band converts
  slow-but-real convergence into an abandoned step.
- **Deck-recoverable:** `nonlinear { stall tolerance : 0.001 ; }` (`cpp:2571-2574`).

### D2 — Progress watchdog is new, default-on, and shorter than the deck target. **P0. Confidence: high.**

- [ ] Resolved
- **Raised by:** Claude (C2), Codex (#2), Grok (P0a) — 3/3.
- `watchdog_magnetic()` `cl_FEM_Controller.cpp:586-612`; default `mWatchdogWindow = 30`
  (`hpp:91`); call site `cpp:1264-1269` → `reset_timestep()`. **No analogue exists on `matfix`.**
- Cuts Δt when no **strictly** new residual minimum has been seen for 30 iterations while ε > 10·tol.
  Two properties bite: (a) `mBestEpsilon` updates only on strict `<` (`cpp:588`), so a plateau or a
  sawtooth returning to an earlier value never re-arms the clock; (b) 30 is **below the deck's own
  `target iterations : 50`**, so a step legitimately needing 60–80 Picard iterations is cut at ~30.
- **Deck-recoverable:** `nonlinear { watchdog window : 0 ; }` (`cpp:2585-2588`).

### D3 — `mOmegaNoiseBand = 0.05` suspends the cited damping rule. **P1. Severity DISPUTED.**

- [ ] Resolved
- **Raised by:** Claude (C3, "most likely mechanism"), Codex (#3, medium), Grok (P2, ~30 %).
- `hpp:110`; branches `cpp:1305-1312` (coupled) and `cpp:1516-1522` (magnetic). New here.
- On `matfix` the rule was strict: `ε < ε₀` → grow ω, **anything else** → `mNumIterationsDiv += 1`
  and `ω *= mAlpha` (0.5). Now `ε₀ ≤ ε < 1.05·ε₀` takes an **empty branch**: ω held *and*
  `mNumIterationsDiv` held, so the `mMaxNumIterationsDiv` relaxation reset (`cpp:1315-1320`) never
  arms either.
- **Literature conflict (found by Claude while verifying Grok's citation complaint, not raised by
  either auditor):** Messe 2023 Eq. 14 is `ω_{k+1}/ω_k = α if ε_{k+1} > ε_k`, else the arctan growth
  term, with α = 0.5, β = 1.1, γ = 0.4 — i.e. exactly `mAlpha`/`mBeta`/`mGamma`
  (`literature/papers/fem/messe2023.txt:574-582`). The rule is **strict**. The neutral band suspends
  precisely it. Meanwhile `cpp:862-863` invokes "Eq. 14: damp on regression and continue" to justify
  removing the Picard line search. **Together the two changes leave a 0–5 % regression band with no
  damping mechanism of any kind** — no line search, no `mAlpha`, no divergence counter.
  This is why Claude ranks D3 well above Grok's 30 %; evidence tier is **literature**, not just source
  trace.
- **NOT deck-recoverable** — no input key.

### D4 — PID timestep control is default-on, holds growth after failures, and has no input key. **P1. Confidence: high.**

- [ ] Resolved
- **Raised by:** Claude (C5), Codex (#4), Grok (P1) — 3/3 on the mechanism.
- `adjust_timestep()` `cl_FEM_Controller.cpp:1921-1990`; `mUseLegacyTimestepControl = false`
  (`hpp:250`) with an explicit "no input key" comment (`hpp:246-249`).
- Steady state (e₀ = e₁ = e₂) collapses the multiplicative PID to `φ = (N_target/N)^kI` with
  kI = 0.30 (`hpp:259`), against legacy `(N_target/N)^0.5`:

  | N (target 50) | legacy φ | PID φ |
  |---|---|---|
  | 25 | 1.414 | 1.231 |
  | 40 | 1.118 | 1.069 |

- Plus `reset_timestep()` sets `mPostFailureHold = 2` which caps φ ≤ 1 for two `adjust_timestep`
  calls (`cpp:1864`, `:1967-1970`), and the `mIteration0` gate skips one more → **≈3 accepted steps
  with zero growth after every cut**.
- **Net:** one cut costs ×0.5; undoing it at N = 25 takes ≈3.3 growth steps *plus* ≈3 held steps ≈ 6
  accepted steps, versus ≈2 on `matfix`. D1/D2 raise the cut rate; D4 lowers the recovery rate ~3×.
  The product is a Δt ratchet.
- **Sub-claim REFUTED:** Grok's "PID can grow Δt *more* aggressively after cheap steps". PID exceeds
  legacy only when `(e₁/e₀)^0.15 · e₀^0.20 > 1` — e.g. e₀ = 0.5, e₁ = 1.5 gives 1.448 vs 1.414, a 2 %
  difference the shared 1.5 clamp (`cpp:1984-1987`) erases immediately. Grok's *second* point
  (post-failure hold prolongs recovery) is the real effect and stands.
- **NOT deck-recoverable.**

### D5 — `reset_timestep()` hard-aborts at the Δt floor; symptom, not cause. **P1. Confidence: high.**

- [ ] Resolved
- **Raised by:** Claude (C6), Grok (2/3).
- `cl_FEM_Controller.cpp:1836-1870` (was a bare `mDeltaTime *= 0.5 ;` on `matfix`). Δt clamps to the
  floor for exactly one attempt, then `BELFEM_ERROR`s. Deck minimum is 1 ns against 0.01 s initial —
  24 halvings, so this is reached only after sustained cut thrash. If Gregory reports an **abort**
  rather than slowness, this is the proximate message; the cause is still D1–D4 driving Δt down there.

### D6 — `mAndersonCommits` is written four times, never read, and documented as active. **P1. Confidence: high. Inert for this deck.**

- [ ] Resolved
- **Raised by:** Grok — **single raiser, needs human adjudication.** Confirmed by direct grep.
- Written `cl_FEM_Controller.cpp:482`, `:497`, `:507`, `:515`; declared `hpp:173-174`. No read exists
  in either file.
- `hpp:171-172` claims "Used to suppress the AIMD growth while Anderson supplies the step direction
  ( O2 )", directly contradicted by `cpp:1281-1286` ("the growth branch stays ACTIVE while Anderson
  is on"). Dead bookkeeping plus a false comment; risks a future O2 patch reintroducing growth
  suppression inconsistently.

### D7 — `watchdog window` is parsed into a `uint` without validation. **P2. Confidence: high.**

- [ ] Resolved
- **Raised by:** Codex — **single raiser.**
- `cl_FEM_Controller.cpp:2585-2588` assigns `get_int()` straight into the `uint` member. The *same
  diff* added `BELFEM_ERROR` guards to `target iterations` (`:2546-2550`) and `anderson depth`
  (`:2576-2582`). A negative value wraps to a huge `uint`, silently *disabling* the watchdog.
  Setup-tier code, so the generous-`BELFEM_ERROR` rule applies.

### D8 — A stale comment says Picard still resets on +9 dB. **P2. Confidence: high.**

- [ ] Resolved
- **Raised by:** Grok — single raiser.
- `cl_FEM_Controller.cpp:1002` still says "a non-escalated Picard iterate still resets on the +9 dB
  threshold as before". The block containing that check (`:1004-1008`) became unreachable for Picard
  when `:876-880` short-circuited the accept path.

### D9 — The line-search comment over-claims Messe 2023 paper1 §4, Eq. 14. **P2. Confidence: high.**

- [ ] Resolved
- **Raised by:** Grok — single raiser, literature-verified.
- `cl_FEM_Controller.cpp:862-863`. Eq. 14 is the α/β/γ relaxation update
  (`literature/papers/fem/messe2023.txt:574-582`); it says nothing about a line search. The
  *engineering* argument for the removal is sound (see D10) — only the citation is over-claimed.
  The §2.7 → §4 renumbering elsewhere in the diff is **correct** (§2.7 is solvers; the nonlinear
  strategy is §4).

### D10 — ~~Picard line-search removal is a regression~~ **FALSE POSITIVE (retracted 2026-08-10, Claude).**

- [x] Retracted — do not re-raise.
- Pre-registered by Claude as C4; self-refuted before the auditors returned; Codex (#5) reached the
  same place independently.
- **Why it is correct:** `cl_FEM_DofMgr_SolverData.cpp:2788-2795` computes the Picard residual
  `r = A·x − b` **deliberately without refreshing `mFieldValues`**, so ε is the **pre-update**
  residual and is exactly ω-independent within a trial. ω is consumed only in the dof update
  (`cl_FEM_DofMgr_SolverData.cpp:2222`, `:2266-2270`), never in the assembly. Since
  `compute_jacobian_and_rhs()` (`cl_FEM_Controller.cpp:831`) re-assembles at the *restored* state on
  every backtrack, the old reject loop produced a **bit-identical residual on all 8 halvings** — it
  could never accept, and cost 8 wasted solves plus an ω collapse.
- The Newton branch by contrast **does** refresh and recompute at the updated iterate
  (`cl_FEM_DofMgr_SolverData.cpp:2213-2237`, new on `sideconnectors`), so retaining the line search
  for Newton only is right.
- **Consequence worth carrying forward:** for a pure-Picard deck **every** controller decision runs on
  a one-iteration-lagged residual. D1, D2 and D3 all act on exactly that lagged signal.

### Non-findings — checked and cleared, do not re-audit

- **Anderson machinery is not the cause.** `anderson stabilization : false` zeroes both depths
  (`cpp:2830`); helpers early-return (`cpp:474-475`, `:492-493`). All three voices agree.
- **Divergence rule 1 → 3 strikes** (`cpp:1226-1228`) is strictly **more lenient** than `matfix`;
  it reduces cuts. Second-order only: the two extra tolerated iterates above +10 dB also land in the
  stall window and the watchdog clock.
- **`scheme` vs `method` key** — `scheme : bdf1` was *ignored* on `matfix` and the default was
  already BDF1 (`cpp:2750-2756`). No behavioural difference.
- **Absolute-tolerance escape** in `run_coupled` / `run_magnetic` (`cpp:2355-2366`) is inert at the
  0.0 default; `mAbsoluteResidual` is a genuine norm set unconditionally
  (`cl_FEM_DofMgr_SolverData.cpp:2424`).
- **New setup-tier Δt-window `BELFEM_ERROR`s** (`cpp:2725-2737`): 0.01 ∈ [1e-9, 0.2], passes.
- **`std::max(std::min(...))` → `std::clamp`** on the ω flips: all inside
  `mAlgorithm == NewtonRaphson` branches, dead for this deck.
- **`BELFEM_DUMP_SYSTEM` blocks** (`cpp:832-845`, `:1410-1423`): env-gated, inert. But they are
  duplicated verbatim with a mis-indented copy at `:1410` and a function-local `static int
  tDumpCount` in each — diagnostic probes that should be stripped before release (P2, housekeeping;
  see `todo/` probe-removal tracking).
  **Update 2026-08-10** — the blind jury round on the controller diff (`review_pid_timestep_impl_b2`,
  Codex + Grok independently) found two substantive defects in them, which raises the priority of
  stripping (or fixing) them above pure housekeeping:
  1. Both dumps fire **before** `impose_voltage_bcs()` (`cpp:451-461`), which patches the RHS for
     `Voltage` / `CircuitVoltage` conditions. On any voltage- or circuit-driven deck,
     `sysdump_*.hdf5` is therefore **not** the system handed to `solve()`. Harmless on
     current-driven decks (the BC loop patches nothing), which is why it has not bitten yet.
  2. The two function-local `static int tDumpCount` are independent counters writing the **same**
     `sysdump_0..3.hdf5` names. Only one iterate path runs per process today, so they do not
     collide in practice — but nothing enforces that.
  Minor, same blocks: `std::getenv` is called on **every** nonlinear assembly (the "production
  paths are untouched" claim is not literally true, though the cost is negligible against an
  assembly), and the accompanying `#include <string>` is redundant (`typedefs.hpp:15`).
  Third jury pass (2026-08-10 evening, Codex): `save_system` opens with `FileMode::NEW` =
  `H5F_ACC_TRUNC` (`cl_HDF5.cpp:44-49`), so with the hard-coded `sysdump_0..3` names every
  diagnostic rerun silently destroys the previous run's dumps — worth `NEW` → timestamped
  names or an existence check if the probes stay.

---

## 3. Ordered Steps

### R1 — Run the deck-recoverable A/B. **Do this first; everything else is gated on it.**

- [ ] Re-run the §0 deck on `sideconnectors` unchanged, capture the full log (Δt history, reset
      frequency, and any `Watchdog: no residual improvement` / stagnation prints).
- [ ] Re-run with **only** the deck changed:

      ```
      nonlinear { tolerance : 1e-7 ; max iterations : 100 ; max relaxation : 1.0 ;
                  target iterations : 50 ; algorithm : Picard ;
                  stall tolerance : 0.001 ;   # restores the matfix band  (D1)
                  watchdog window : 0 ; }     # disables the new watchdog (D2)
      ```

- [ ] Re-run the same deck on `matfix` for the baseline Δt trace.
- **Decision rule:** recovers ⇒ D1 + D2 are dominant and the fix is the defaults (R2). Does not
  recover ⇒ D3 / D4 are implicated; they have **no input key**, so R3 is needed to bisect.

### R2 — Fix the algorithm asymmetry rather than the numbers *(after: R1, if it recovers)*

- [ ] Decide with Christian whether the retuned guards should be conditioned on
      `mAlgorithm == NewtonRaphson` (§1.1), or whether the defaults themselves should move.
      Conditioning is the smaller change and addresses the actual root cause; renumbering the
      defaults risks re-regressing the quench cases the retune was for.
- [ ] Whatever lands, update `doc/input_file_reference.md` in the same turn (input-reference sync
      rule) and have Codex do a prose pass on the touched sections.

### R3 — Give `mUseLegacyTimestepControl` and `mOmegaNoiseBand` an input surface *(after: R1, if it does not recover)*

- [ ] Temporary flag (or deck keys) so D3 and D4 can be bisected without a recompile per trial.
      Both auditors flagged the absence of an input surface as a **process risk**, independent of
      whether these are the culprit: a behaviour-changing default that cannot be A/B'd from a deck
      cannot be regression-tested either.

### R4 — Clear the confirmed side defects *(independent of R1; none of these affect the tapestack)*

- [ ] D6 — either wire `mAndersonCommits` into the growth decision or delete it, and fix
      `hpp:171-172` either way.
- [ ] D7 — add the `BELFEM_ERROR` guard to `watchdog window`, matching its siblings.
- [ ] D8 — fix the stale +9 dB comment at `cpp:1002`.
- [ ] D9 — correct the Eq. 14 citation at `cpp:862-863` (the *argument* is right; cite the residual
      semantics of `cl_FEM_DofMgr_SolverData.cpp:2788-2795` instead).
- [ ] Strip the duplicated `BELFEM_DUMP_SYSTEM` probes at `cpp:832-845` / `:1410-1423` — the
      preferred disposition now that DR-02 Gate A has served its purpose. If they are still
      needed, fix in place instead: one member flag read once in the constructor + one member
      counter shared by both paths, and move the call **after** `impose_voltage_bcs()` so the
      dump is the system actually solved (both defects above; proposal recorded in the
      2026-08-10 session).

### R5 — Regression gate *(after: R2 or R3)*

- [ ] Tapestack deck completes 15 s without cut thrash, Δt trace comparable to `matfix`.
- [ ] A quench case still behaves as well as it does today — the retune must not be undone.
- [ ] `make check-fast`.

---

## 4. Open Design Questions

**O1 — How severe is D3 really?** Claude: primary mechanism (literature-backed — Eq. 14's strict rule
is suspended). Codex: medium. Grok: ~30 %. The band was tuned on sidecoatings, which has a residual
floor *above* tolerance; the tapestack has no such floor. **Not resolvable by code reading — needs
the R1 traces, then Christian's physics judgement.** Physics/design questions are never settled by
vote.

**O2 — Should the retuned guards be algorithm-conditional (§1.1)?** Options: (a) gate the new
defaults on `mAlgorithm == NewtonRaphson`; (b) soften the defaults globally; (c) leave the defaults
and document that pure-Picard decks must override them. (a) targets the root cause; (b) risks
re-regressing quench; (c) is the status quo and pushes the cost onto every user.

> **New evidence against (a) (2026-08-10, Claude).** The `sidecoatings` collapse analysed in
> `todo/timestep_collapse_residual_floor_plan.md` §1 happens on a deck with `algorithm : Newton`,
> and it dies on `mStallBand` at a measured MAD of 0.199 dB. So the widened stall band harms Newton
> decks too, and conditioning it on `mAlgorithm == NewtonRaphson` would not have saved that run.
> The asymmetry of §1.1 is still real, but it is not the whole story — prefer fixing the stall test
> itself (that plan's M1). Confidence: high.

**O3 — Is the regression definitely controller-local?** The assumption throughout is yes. It has not
been tested: the `matfix` → `sideconnectors` diff also rewrote `cl_ThinShellFactory` (+900 lines),
which is squarely on this deck's path (`thinshell : tape`, 4 layers). If R1 fails to reproduce a
controller-driven cut pattern, **re-open this question before tuning anything.**

---

## 5. Definition of Done

- [ ] R1 traces captured on both branches and both deck variants.
- [ ] O1 adjudicated by Christian with the traces in hand.
- [ ] O2 decided and implemented.
- [ ] Tapestack no longer thrashes; a quench case is unchanged (R5).
- [ ] D6–D9 closed or explicitly deferred with a reason.
- [ ] `doc/input_file_reference.md` synced for any new/changed key.
- [ ] Devlog entry written; this file moved to `todo/closed/` and `todo/README.md` updated.

---

## 6. Audit Trail

- Exchange thread: `tmp/ai_exchange/review_controller_tapestack_regression.md` — full record of the
  2026-08-10 `--jury` round: Claude's frozen pre-registration, the blind auditor entries, citation
  verification, and reconciliation table. **Ephemeral — distil anything still needed into this file
  before it is swept.**
- Working diff used as the auditors' subject:
  `git diff matfix -- src/fem/kernel/cl_FEM_Controller.cpp src/fem/kernel/cl_FEM_Controller.hpp`
  (regenerate rather than relying on `tmp/ai_exchange/_subject_controller_tapestack.md`).
- **Per-AI contribution:** Codex caught D7 and independently confirmed D1/D2/D4 plus the D10
  retraction. Grok caught D6, D8, D9 and framed the algorithm asymmetry of §1.1 most sharply. Claude
  pre-registered D1–D5, self-retracted D10 before dispatch, and found the D3 literature conflict
  while verifying Grok's D9. Every auditor citation was re-checked against the tree before it was
  recorded here; two auditor sub-claims were refuted (D4 growth aggressiveness; D10).
- **Nothing here is measured.** No run was performed in the review session. All rankings are
  control-flow certainty. Treat R1 as the first real evidence.
