# Mitigation Design: Controller Stall/Handoff Guards and Δt-Cut Diagnosis

> **CANCELLED 2026-09-03** (todo/ currentness sweep, round 3): superseded — designed before the cause was known; the κ ∝ Δt finding of the penalty opt-in work (2026-09-01) explains the collapse the M-mitigations here were guessing at. Status lines and checkboxes below are as they stood at closure and are not maintained.

**Date:** 2026-08-10
**Purpose:** Concrete, code-level mitigation design for the Δt-collapse failure analysed in
`todo/timestep_collapse_residual_floor_plan.md` (the `sidecoatings` t = 13.8762 ms abort and the
`greg5` regression). Mechanism being fixed, in one sentence: the controller's stall detector cannot
tell "flat" from "slowly descending", its Picard↔Newton handoff chatters on a fixed threshold, and
its recovery moves (trust growth, escalation ω) are tuned for contraction, not plateaus — so the
guards cut timesteps that were converging.
**Module:** `src/fem/kernel` (`cl_FEM_Controller.{cpp,hpp}`), `doc/input_file_reference.md`
**AIs involved:** Claude (design), Codex + Grok (jury audit of this document — pending)
**Status:** DESIGN — **not audited, not approved, no source modified.** Line numbers cite the
current `sideconnectors` working tree (`1d6ef305` + uncommitted when written; that tree is
now committed — `bc578b5e`, verified clean 2026-08-10, so re-anchor line numbers there).

> **Scope guards:**
> - This designs **mitigations only**. The diagnosis it rests on is
>   `todo/timestep_collapse_residual_floor_plan.md` §1–§2 (referenced below as **TC-D1…TC-D7**);
>   the run gates R1–R4 there remain the evidence plan and are *not* replaced by this document.
> - **Explicitly out of scope:** `mOmegaNoiseBand` (disputed, awaits R1 traces + Christian —
>   companion file O1), the progress-watchdog window (deck-recoverable), PID gain retune
>   (`todo/closed/pid_timestep_controller_plan.md` owns it), any change to the Picard residual semantics
>   (settled — companion D10 retraction), source-table smoothing beyond a warning (TC-D6 / M6:
>   PCHIP interpolation is deferred).
> - Every default below is chosen under the constraint **"a quench deck that converges today must
>   behave identically or better"**; each MIT carries an explicit quench-regression argument. The
>   final arbiter is the R5 regression gate, not the argument.
> - Implementation only after Christian approves; route through a fresh Fable session per the
>   critical-refactor policy if approved.

---

## 1. What is being changed, at one glance

| ID | Change | Fixes | Risk | Behaviour change for a healthy quench run |
|---|---|---|---|---|
| MIT-4 | Δt cuts state their reason; floor abort reports measured evidence | TC-D7 | none | log lines only |
| MIT-5 | Log resolved BDF scheme; error on `scheme`/`method` conflict; warn on order ≥ 3 + adaptive Δt; growth clamp 1.2 → 1.12 for order ≥ 4 | TC-D5 | minimal | none (BDF1/BDF2 decks) |
| MIT-1 | Stall = flat **trend**, not small scatter | TC-D2 | low | none (a floored residual has slope ≈ 0 and still fires) |
| MIT-2 | Hysteresis + dwell + failed-promotion latch on the Picard↔Newton handoff | TC-D1 | medium | none (one promotion, residual keeps falling → identical path) |
| MIT-3 | Newton trust growth needs a 2-streak; escalation enters at 2·ω_Picard, not ω₀ | TC-D3, TC-D4 | low | trust ×2 arrives one accepted step later |
| MIT-6 | Last-resort grind mode before the Δt floor | TC-D8, symptom of all | medium | none until ≥ 3 consecutive cuts |
| MIT-7 | Warn when a tabulated/user source has sign-alternating dI/dt at quantisation scale | TC-D6 | none | log line only |

Rollout in that order: MIT-4/5 are pure instrumentation and land first so every later A/B run
produces attributable logs; MIT-1/2/3 are the behavioural fixes; MIT-6 is the safety net; MIT-7 is
independent.

---

## 2. MIT-1 — Trend-aware stall detector *(TC-D2; supersedes companion D1's number-tuning)*

### Current code

`magnetic_stagnation_forces_reset()` (`cl_FEM_Controller.cpp:641-700`) declares a stall when the
**mean absolute deviation** of the last `mStallWindow = 5` dB samples (`hpp:82`) falls below
`mStallBand = 0.2 dB` (`hpp:83`). MAD measures scatter; a perfectly monotone descent of
0.17 dB/iterate has MAD ≈ 0.2 and is indistinguishable from a flat floor. That is the measured
kill at `sidecoatings` iteration 133 (MAD 0.1993 vs band 0.200) and the arithmetic behind the
`greg` tapestack cuts (any contraction slower than 3.8 %/iterate reads as a stall).

### Design

Replace the MAD test with a two-part test over the same window:

1. **Trend:** least-squares slope `b` (dB per iteration) of the window samples against their index.
2. **Scatter:** RMS deviation of the samples about the fitted line.

Stall **iff** `b > -mStallSlope` **and** `rms < mStallBand`.

```cpp
// inside magnetic_stagnation_forces_reset(), replacing the mean/MAD block
// x = 0..n-1, y = dB samples; closed-form least squares
real tSxx = 0.0, tSxy = 0.0, tMeanY = 0.0 ; // ( x mean is (n-1)/2 )
...
const real tSlope = tSxy / tSxx ;                 // dB / iteration
const real tRms   = std::sqrt( tSse / tN ) ;      // about the fit
if ( tSlope > -mStallSlope && tRms < mStallBand ) { /* stall path as today */ }
```

- New member `real mStallSlope = 0.05 ; // dB/iteration` (`hpp`, next to `mStallBand`): a residual
  descending faster than ~1.1 %/iterate is never a stall, no matter how quiet.
- `mStallBand` keeps its value (0.2 dB) and its key, but its meaning tightens from "total scatter"
  to "scatter about the trend". For a flat series the two coincide (fit slope ≈ 0, rms ≈ MAD up to
  the usual 1.25 factor), so floored-residual behaviour is preserved.
- New input key `stall slope` (nonlinear section, parsed beside `stall tolerance` at
  `cpp:2571-2574`), validated `>= 0` with a setup-tier `BELFEM_ERROR` (error-tier-by-call-rate
  rule). `stall slope : 0` restores exactly today's fire condition (slope gate always passes).
- The **thermal** side has no stagnation guard (only `watchdog_thermal`), so nothing to mirror.

### Checks against the three known traces

| trace | window slope | fires? | correct? |
|---|---|---|---|
| `sidecoatings` it 129–133 (−42.84…−43.53 dB) | −0.17 dB/it | **no** (today: yes) | yes — it was converging |
| quench residual floor (flat at −55 dB ± noise) | ≈ 0 | **yes**, as today | yes |
| tapestack slow Picard, ρ = 0.97 (−0.13 dB/it) | −0.13 | **no** (today: yes) | yes — also fixes companion D1 without touching the band |
| genuine ω-sawtooth limit cycle | ≈ 0, large rms | no — **by scatter**, as today | yes: that case belongs to the watchdog (`cpp:1264-1269`), which is untouched |

- [ ] Implemented
- [ ] `doc/input_file_reference.md` synced (`stall slope`, revised `stall tolerance` semantics)

**Quench-regression argument:** the guard was retuned for the sidecoatings-class residual *floor*;
a floor has slope ≈ 0 and still trips. Only descending windows change outcome. Confidence: high.

---

## 3. MIT-2 — Handoff hysteresis, dwell, and a failed-promotion latch *(TC-D1)*

### Current code

`cpp:748`: promote to Newton the moment a **single** ε sample crosses the fixed
`mEpsilonSwitch = 1e-4` (`hpp:112`); demote the moment it does not. No dwell, no hysteresis, no
memory of how the last promotion went. On a plateau sitting *on* the threshold this switched
algorithms six times in one attempt, at 6–18 dB per switch (companion plan §1). `mFirstFlip`
(`cpp:762`) damps only the first promotion.

### Design

Three additions to the flip block at `cpp:744-788`, all per-attempt state reset in
`initialize_timestep` / `reset_timestep`:

1. **Dwell.** New members `uint mSwitchDwell = 3 ;` and `uint mLastSwitchIteration = 0 ;`. Neither
   promotion nor demotion may fire within `mSwitchDwell` iterations of the previous flip. The
   existing `mIteration > 1` promotion guard is subsumed.
2. **Hysteresis.** Promotion stays `ε < mEpsilonSwitch`. Demotion becomes
   `ε > mSwitchHysteresis · mEpsilonSwitch` with `real mSwitchHysteresis = 10.0 ;` — Newton is not
   abandoned for wandering a factor of a few above the switch point; it is abandoned when it has
   demonstrably kicked the residual a decade up. Between the two thresholds the current algorithm
   holds.
3. **Failed-promotion latch.** On promotion, record `mPromotionEpsilon = mEpsilon`. If, while
   Newton runs, `ε > 4 · mPromotionEpsilon` (a +6 dB regression against the promotion entry —
   distinct from the hysteresis bound, which is absolute), demote immediately **and** set the
   existing `mJustPicard` latch: this attempt has proven the tangent unhelpful at this state, so
   Picard finishes the timestep alone. This is exactly the mechanism that rescued the failing
   attempt at iteration 82, made deliberate and cheap (it 3–4 would have latched it, saving ~78
   iterations).

Interaction notes (each needs an auditor's eye):
- The **escalation** path (`try_escalate_to_newton`, `cpp:531-587`) bypasses this flip block by
  design (sets the algorithm directly, comment at `cpp:770-775`) — unchanged.
- The stagnation guard's own Newton→Picard demotion (`cpp:668-690`) also sets `mJustPicard` and is
  unchanged; the latch here is additive, not a replacement.
- `mFirstFlip`'s 0.5 damping on the first promotion stays. With the latch, "promotions 2..n enter
  undamped" (TC-D1) loses most of its bite because a promotion that regresses no longer recurs.
- The stall window is already cleared on every algorithm change (`cpp:790-794`) — hysteresis makes
  those clears rarer, which also lengthens the effective evidence window of MIT-1. Beneficial, but
  the auditors should confirm no guard *depends* on frequent clears.
- New input keys `handoff dwell` (uint, ≥ 1) and `handoff hysteresis` (real, ≥ 1.0), setup-tier
  validated. Behaviour-changing defaults must be A/B-able from a deck (companion R3 principle).

- [ ] Implemented
- [ ] `doc/input_file_reference.md` synced

**Quench-regression argument:** a quench deck promotes once and the residual keeps falling — dwell
and hysteresis never engage on a monotone pass through the threshold; the latch requires a +6 dB
regression that today's healthy runs do not exhibit (if one did, it would burn iterations exactly
as sidecoatings did). Confidence: high for dwell/hysteresis, medium for the latch threshold value
(the 4× factor is chosen from one trace; the jury should challenge it).

---

## 4. MIT-3 — Trust growth needs a streak; escalation enters at Picard's earned ω *(TC-D3, TC-D4)*

### 4a. Trust streak

`cl_FEM_Controller.cpp` ( the `tGrowth = 2.0` block in `iterate_coupled`; the old `cpp:1287-1294`
citation is stale ): one first-trial Newton step that halves ε sets `tGrowth = 2.0`.

**Correction, 2026-08-18 (three-vendor round):** the original "2 of 2" reading does not survive
arithmetic. In the tapestack3d step-26 trace only ONE event is this rule — ω 0.483 → 0.966. The
second suspected doubling, 0.580 → 0.693, is a ratio of 1.195, i.e. ordinary arctan growth, and
the following 0.693 → 0.347 is the punishment factor. MIT-3a therefore addresses ONE of the two
ω excursions in that trace; the other belongs to the arctan rule and needs the D9 ω-ceiling.

Change: new member `uint mNewtonTrustStreak = 0 ;` (reset on every reject, demotion, promotion,
and timestep start). The qualifying condition (`tRanNewton && tBacktracks == 0 &&
mEpsilon < 0.5 * mEpsilon0`) increments the streak; anything else zeroes it. `tGrowth = 2.0`
requires `mNewtonTrustStreak >= 2`; the first qualifying step gets the ordinary arctan growth
(≤ 1.3, with 0.5·ε contraction: ≈ 1.29).

Cost to the motivating `ts34` endgame crawl: the geometric recovery arrives one accepted step
later (ω ×1.29 then ×2 instead of ×2 ×2). Two iterations of delay against the 11 the rule was
built to save — acceptable, and the plateau failure mode (one lucky step near a floor) can no
longer fire it.

### 4b. Escalation entry ω

`cpp:561`: `mOmegaNewton = mOmega0 ;` re-enters at full relaxation on a state where the controller
has just spent an entire window learning the stable ω range (the failing log: Picard itself
diverged at ω = 0.494; escalation entered at 1.0 and lost 27 dB).

Change: `mOmegaNewton = std::clamp( 2.0 * mOmegaPicard, mOmegaMin, mOmega0 ) ;` — the tangent gets
double the headroom Picard had earned, never more than the configured start. The stated worry in
the current comment ("a collapsed Picard ω would freeze the Newton step") is bounded by the 2×
factor plus the growth rule; a Newton that is genuinely contractive earns its way up within 2–3
accepted steps (and MIT-3a's streak then doubles it).

- [ ] Implemented (4a)
- [ ] Implemented (4b)

**Quench-regression argument:** 4a delays one growth step on genuinely contracting Newtons; 4b only
affects the escalation path, which exists for Picard-breakdown rescue and fires on struggling runs
by definition. Confidence: high (4a), medium (4b — the 2× factor is a judgement call; the
alternative "largest ω that produced an accepted step this attempt" needs an extra tracked member
and is offered to the jury as O2).

---

## 5. MIT-4 — Every Δt cut states its reason; the floor abort reports evidence *(TC-D7)*

### Design

- New `enum class TimestepCutReason : uint8_t { Divergence, PicardStall, Watchdog,
  ThermalWatchdog, IterationBudget, BacktrackBudget, ThermalDivergence, SolverFailure, NanResidual }`
  (`cl_FEM_Controller.hpp`).
- `reset_timestep()` gains a required `aReason` parameter plus an optional detail string; **all
  nine call sites** are updated (`cpp:866`, `:1257-1261`, `:1264-1269`, the thermal watchdog site
  `:1183-1186`, `:1251-1255` divergence, the backtrack site `:1004-1008` region, `run_coupled` /
  `run_magnetic` failure paths). One line prints inside the box before the closing rule:
  `│ Δt cut: Picard stall ( slope -0.01 dB/it, rms 0.08 < 0.20 ) at eps 4.4e-05, it 133 │`.
- A ring buffer (reuse `ShiftRegister`, depth 8) records per cut: reason, Δt, ε at iteration 1,
  best ε, iterations spent. The **floor abort** (`cpp:1847-1857`) appends this table and, when ≥ 2
  entries share a reason, the measured ratio ε₁(Δt)/ε₁(Δt/2) — the empirical Δt-response slope of
  companion R3, produced for free. The abort text's unconditional claim ("the nonlinear residual is
  not responding to the timestep size") is replaced by whichever of the two diagnoses the recorded
  ratios support, falling back to "insufficient data" wording.
- Cold path only (once per cut / once per abort); no allocation at cut time (ring sized in the
  constructor per the member-scratch rule).

- [ ] Implemented
- [ ] Abort message wording reviewed (it is user-facing)

**Quench-regression argument:** log-only. Confidence: high.

---

## 6. MIT-5 — BDF scheme transparency and the order-≥4 growth clamp *(TC-D5)*

- After the scheme parse (`cpp:2750-2796`): unconditional startup log line naming the resolved
  scheme **and the key it came from** (`scheme` / `method` / default).
- If **both** keys exist and disagree → setup-tier `BELFEM_ERROR` (silent precedence is how the
  `greg5` deck changed integration order without a log line). If both exist and agree → accept.
- If resolved order ≥ 3 **and** `adapt timestep : true` → loud warning citing the step-ratio
  stability bounds (the variable-step high-order path is unverified;
  `todo/closed/bdf_nonlinear_mass_verification.md` tracks the verification).
- Growth clamp (`cpp:1980-1983`): `order_active() >= 4` tier drops 1.2 → **1.12** (classical BDF5
  zero-stability ratio bound ≈ 1.13; sit under it). BDF3's 1.4 and the default 1.5 stay.
- **Not** changing the `scheme`-over-`method` precedence itself: it is now released behaviour on
  this branch, and the conflict error covers the dangerous case.

- [ ] Implemented
- [ ] `doc/input_file_reference.md` synced (scheme/method conflict now an error)

**Quench-regression argument:** quench decks run BDF1/BDF2; the clamp change touches order ≥ 4
only. Confidence: high.

---

## 7. MIT-6 — Last-resort grind mode before the floor *(safety net; the accidental it-82 rescue, made deliberate)*

### Design

- New member `uint mConsecutiveCuts = 0 ;` — incremented in `reset_timestep()`, zeroed on every
  **accepted** timestep.
- When `mConsecutiveCuts >= mGrindThreshold` (default 3), the next attempt runs in **grind mode**,
  announced in the log:
  - algorithm latched to Picard for the attempt (`mJustPicard = true` — same mechanism as today's
    Newton-stall demotion, so no new state semantics);
  - ω started at `min( mOmegaPicard, 0.25 )` and adapted only by the existing arctan/α rules
    (promotion, escalation, trust jump all masked);
  - MIT-1 stall guard and the magnetic watchdog disabled **for this attempt only**;
  - iteration budget for the attempt raised to `2 · mMaxNumIterations`;
  - divergence guard (+10 dB, 3 strikes) and the NaN guard stay armed — grind mode must not
    protect a genuinely diverging run.
- Outcome: converges → counter cleared, normal policy resumes next step. Fails → `reset_timestep`
  as usual; when the floor is then reached, the MIT-4 abort can truthfully say that an
  unsupervised low-ω Picard run at this Δt was tried and failed — the abort's diagnosis becomes
  *earned*.
- Input key `last resort : <n> ;` (0 = off) in the nonlinear section. **Default ON at 3** is the
  designed intent — the failing log is the direct evidence (iterations 82–127 converged under
  exactly this policy) — but the default is flagged to Christian as **O1** below, since it changes
  behaviour on every deck that today reaches 3 consecutive cuts.

- [ ] Implemented
- [ ] O1 decided (default on/off)
- [ ] `doc/input_file_reference.md` synced

**Quench-regression argument:** inert until three consecutive cuts, a state that today leads to
the floor spiral in the observed traces; a quench deck that legitimately needs Δt cut 3× in a row
pays one grind attempt (bounded by 2× the iteration budget) before cutting resumes. That cost is
real and is the reason for the O1 flag. Confidence: medium.

---

## 8. MIT-7 — Non-smooth source warning *(TC-D6; warn-only this round)*

At tabulated/user-defined source load (the `cl_SourceFunction` user-defined path), if the sampled
first differences alternate sign in more than half of any 100-sample stretch while `|ΔI|` stays
within a small multiple of its own median (the quantisation-dither signature), print a one-time
warning naming the file and suggesting pre-filtering. No behaviour change; PCHIP interpolation for
tabulated sources is deferred to its own task (needs monotonicity-preservation design and belongs
with the source-function code, not the controller).

- [ ] Implemented (warning only)

---

## 9. Explicitly not done, and why

| item | why not |
|---|---|
| Revert `mStallBand` to 0.001 | MIT-1 removes the failure mode without renumbering; 0.001 would undo the quench tuning the band was set for |
| Condition guards on `mAlgorithm == NewtonRaphson` (companion O2 option a) | refuted as sufficient: `sidecoatings` **is** a Newton deck and died on the same guards (companion O2 amendment, 2026-08-10) |
| `mOmegaNoiseBand` changes | disputed severity (companion O1); awaits R1 traces + Christian |
| Watchdog window/logic changes | deck-recoverable (`watchdog window : 0`); watchdog semantics (strict new minimum) are orthogonal to the MAD defect |
| PID gains / post-failure hold | owned by `todo/closed/pid_timestep_controller_plan.md`; MIT-6 reduces how often the hold matters |
| `mEpsilonSwitch` becoming relative (contraction-ratio based) | right long-term question, wrong blast radius for this round — carried as O3 |

---

## 10. Open questions for Christian

**O1 — MIT-6 default:** last-resort grind mode ON at 3 consecutive cuts, or OFF (opt-in key)?
Designed intent is ON; the counterargument is the bounded-but-real cost on decks that legitimately
cut repeatedly.

**O2 — MIT-3b entry ω:** `2 · mOmegaPicard` (stateless, proposed) vs "largest ω that produced an
accepted step this attempt" (one more tracked member, tighter). Either bounds the observed +27 dB
escalation kick.

**O3 — carried:** should `mEpsilonSwitch` be relative to the attempt's residual trajectory instead
of absolute? (Companion plan O2; out of scope here.)

---

## 11. Verification plan (after approval, after implementation)

1. **Unit-level:** MIT-1 slope/rms on synthetic windows (flat, monotone 0.17 dB/it, sawtooth);
   MIT-2 flip sequence on a scripted ε trace reproducing §1 of the companion plan (expect: 1
   promotion, 1 latch, 0 further flips).
2. **Case A/B (Christian runs, per policy):** `sidecoatings` cold start through 13.8 ms;
   `greg5` with `scheme : bdf1` (TC-R1) and, separately, with the mitigations on `bdf5`;
   one quench deck unchanged — Δt trace and iteration counts vs today.
3. **`make check-fast`** (USE_TEST toggle per the shared-tree rule).
4. Every landed key documented in `doc/input_file_reference.md` in the same turn + Codex prose
   pass on the touched sections.

---

## 12. Audit Trail

- Diagnosis basis: `todo/timestep_collapse_residual_floor_plan.md` (2026-08-10, this session) and
  `todo/controller_picard_tapestack_regression.md` (2026-08-10 jury round).
- This document's jury round: pre-registration + exchange thread under
  `tmp/ai_exchange/` — see §13 reconciliation once complete.
