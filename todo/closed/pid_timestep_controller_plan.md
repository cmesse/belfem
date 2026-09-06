# PID Timestep Control for the FEM Controller (with Legacy Switch)

**Date:** 2026-08-09
**Purpose:** Replace the memoryless iteration-target sqrt rule and the blind `Δt *= 0.5` cut with
a PID-filtered step-size controller in log space, so that the timestep responds smoothly to the
solver's cost signal instead of oscillating between over-growth and cut cascades ("bumping against
the 1 ns floor"). The current behaviour stays available behind an internal member switch (legacy
mode). Two independent hardening items ride along: the minimum-timestep floor is actually enforced
(today `reset_timestep` halves without bound), and the memdump warm start restores the crashed
run's Δt instead of re-entering the transient at the deck's initial timestep.
**Module:** `src/fem/kernel` (`cl_FEM_Controller.hpp/.cpp`); one touch in `Controller::save_memdump`
/ `load_memdump` (no `Mesh` API change)
**AIs involved:** Claude (plan + draft), Codex + Grok (jury audit of this file before any source
edit), Christian (policy decisions O1–O3)
**Status:** IMPLEMENTED 2026-08-09 (R3+R4 landed) and **COMMITTED in `1d6ef305`** — "added
PID timestepping", 2026-08-09; commit state verified 2026-08-10, working tree clean — on Christian's
execute order; O2/O3 resolved (abort at floor + C5 escape hatch; PID default). Remaining:
round-3 diff jury (R6, in flight), R5 A/B validation (Christian: legacy vs PID, bdf1 vs bdf5),
O1 gain calibration.

> **Scope guards:**
> - NO change to the nonlinear iteration itself (relaxation AIMD, Picard/Newton handoff, watchdog,
>   stall band). Those are separate findings in the running review round (C2–C6) and must not be
>   entangled with the Δt controller.
> - NO local-truncation-error estimator in phase 1 (that needs IWG_Timestep support — deferred,
>   see O4). Phase 1 keeps the iteration count as the control signal, same as legacy.
> - NO new input-file keys in phase 1 — the mode switch is an internal member per Christian's
>   instruction (2026-08-09). Deck exposure comes after validation (and then updates
>   `doc/input_file_reference.md` same turn, per standing rule).
> - The thermal sub-step machinery (`mCouplingFactor`, `reset_thermal`) is untouched.
> - The save-grid trimming block in `adjust_timestep` (`cl_FEM_Controller.cpp:1837-1875`) is
>   untouched; only the φ (growth-ratio) computation is mode-switched.

---

## 1. Current Behaviour and How It Fails

Controller state, working tree `sideconnectors` (2026-08-09):

| Failure | Mechanism | Evidence |
|---|---|---|
| Δt collapses through the configured floor | `reset_timestep()` does `mDeltaTime *= 0.5` with no clamp and no abort; the only `mDeltaTimeMin` clamps sit in `adjust_timestep()`, which runs only for accepted steps (`finalize`, gated on `!mReset`) | `cl_FEM_Controller.cpp:1746` vs `:1826,1861`; sidecoatings log prints `Δ t = 0.0000 ms` against a 100 ns deck minimum |
| Growth is memoryless and re-climbs into the same cliff | accepted-step ratio `φ = sqrt(target/N)` uses only the last step's N; after a cut cascade the very next accepted steps grow at the clamp ceiling again | `cl_FEM_Controller.cpp:1802,1817` |
| Failure response is a fixed halving with no history | every cut is ×0.5 regardless of how the attempt failed or how many cuts preceded it; only SINGULAR soft-fails are counted (`mSolverFailCount`, aborts at 8 — `cl_FEM_Controller.cpp:831,1960`), while non-singular cuts (watchdog, max-iterations, residual floor — the sidecoatings case) grind unbounded | `cl_FEM_Controller.cpp:1746`; log: identical −46 dB floor at every halved Δt |
| Warm start re-enters the transient at the deck Δt | memdump meta stores timestep index / timestamp / checksum only; `set_params` re-installs `initial timestep` | `cl_Mesh.cpp:2992-3017`, `cl_FEM_Controller.cpp:2522-2527` |

**Bottom line:** the step-size logic is a high-gain, memoryless I-controller with an unbounded
failure path; it oscillates between ceiling growth and cut cascades, and nothing distinguishes
"this Δt is slightly too big" from "no Δt will ever converge".

## 2. Architecture: PID in Log Space, Legacy as a Degenerate Case

**Key observation:** the legacy rule is already a pure I-controller. With the cost error
`e = N / N_target`, legacy computes `φ = e^(−1/2)` — an integral controller with gain 0.5 and no
memory. The control literature's diagnosis of exactly this controller class (Söderlind 2002) is
that high-gain I-control on a noisy signal produces step-size oscillation and rejection bursts;
the remedy is PI/PID filtering with modest gains. So the new mode is a strict generalization:
**kP = kD = 0, kI = 0.5 reproduces legacy exactly** (same clamps), which gives a clean A/B and a
trivially safe fallback.

Controller form (Valli, Carey & Coutinho 2002, multiplicative/log-space PID):

```
φ = (e₁/e₀)^kP · (1/e₀)^kI · (e₁²/(e₀·e₂))^kD
```

with `e₀ = N_now/N_target`, `e₁`, `e₂` the two previous accepted-step errors (all initialized
to 1). φ then passes through the existing order-aware ratio clamp (`tPhiMax` = 1.5/1.4/1.2 by
active BDF order, `cl_FEM_Controller.cpp:1809-1817`), which is kept verbatim in both modes — it
encodes the variable-step BDF ratio bounds and is orthogonal to the controller choice.

Failure path (both modes): cut ×0.5 **clamped at `mDeltaTimeMin`**; a cut demanded while already
at the floor aborts loudly with a diagnostic (the residual is not responding to Δt — the disease
is elsewhere). This is the ONLY new terminal condition, and it is philosophy-consistent: a
rejected timestep is an "expected algorithmic failure" that must feed the retry policy, and
aborting is legitimate only where retry is exhausted (`doc/coding_philosophy.md:621-624`) — which
at the floor it is. The originally drafted consecutive-cuts abort was DROPPED after audit
(§6 D3): it could fire while Δt was still above the floor, i.e. while recovery remained possible.
PID mode additionally (SUNDIALS practice, Brenan–Campbell–Petzold §5): resets the error history
to 1 (a rejected attempt's cost says nothing about the next accepted step's, and stale history
would kick φ through the P/D terms) and holds growth (`φ ≤ 1`) via `mPostFailureHoldSteps = 2`,
so the controller cannot re-climb over the same cliff it just fell off — this is the specific
anti-"bump against 1 ns" mechanism, and it also tempers the warm-start re-entry.

**Hold semantics, exactly (audit R3):** the first accepted step after a cut skips
`adjust_timestep` entirely — the existing `mIteration0 > 0` gate sees the zeroed counter from
`reset_timestep` (`cl_FEM_Controller.cpp:1710,144-145,1795`) — so no growth and no hold
decrement happen there. The hold therefore suppresses growth on the *second and third* accepted
steps after the cut; combined with the forced skip, growth resumes on the fourth. This is
documented behaviour, not a defect: the net effect ("no re-growth for ≥3 accepted steps after a
cut") is what the mechanism is for.

**Legacy-equivalence caveat (audit, Codex-2):** the A/B point (kP, kI, kD) = (0, 0.5, 0) matches
the legacy *growth formula* but the PID mode still applies the post-cut hold and history reset.
For a strict controller-formula A/B, also set `mPostFailureHoldSteps = 0`; otherwise the lever
tests "legacy gain + rejection hold".

**Rejected alternative:** controlling on solution change per step (`‖Δq‖/‖q‖`, Valli et al.'s
original signal) or on a predictor–corrector LTE estimate. Both are better signals and both are
deferred (O4): the first introduces a new tolerance with physical-units semantics, the second
needs IWG history access. Phase 1 must be signal-compatible with legacy to be comparable.

### Literature basis

None of these are in `./literature/` (checked 2026-08-09); cited per the fallback rule with DOIs:

- Gustafsson, Lundh & Söderlind 1988, *A PI stepsize control for the numerical solution of
  ordinary differential equations*, BIT 28:270–287. doi:10.1007/BF01934091 — PI replaces the
  elementary (deadbeat-I) controller; kills rejection oscillation.
- Gustafsson 1994, *Control-theoretic techniques for stepsize selection in implicit RK methods*,
  ACM TOMS 20(4):496–517. doi:10.1145/198429.198437 — PI gains for implicit methods; convergence-
  failure restart strategy.
- Söderlind 2002, *Automatic control and adaptive time-stepping*, Numerical Algorithms
  31:281–310. doi:10.1023/A:1021160023092 — the control-theory framing used in §2; why high-gain
  I-control oscillates.
- Söderlind 2003, *Digital filters in adaptive time-stepping*, ACM TOMS 29(1):1–26.
  doi:10.1145/641876.641877 — H211b filter; smooth step sequences improve multistep (BDF)
  stability ("computational stability", also Söderlind & Wang 2006, JCAM 185:225–243,
  doi:10.1016/j.cam.2005.03.008).
- Valli, Carey & Coutinho 2002, *Control strategies for timestep selection in simulation of
  coupled viscous flow and heat transfer*, Commun. Numer. Meth. Engng 18(2):131–139.
  doi:10.1002/cnm.475 — the multiplicative PID form adopted here; demonstrates PID step control
  without an embedded LTE estimator, on coupled nonlinear flow+heat problems (structurally the
  closest match to BELFEM's coupled magnetic+thermal setting).
- Brenan, Campbell & Petzold 1996, *Numerical Solution of IVPs in DAEs* (SIAM), Ch. 5 — BDF
  production practice: cut factor on convergence failure, growth suppression after rejection,
  bounded consecutive-failure count with loud abort.

The gains question (O1): Valli et al. tuned (kP, kI, kD) = (0.075, 0.175, 0.01) for an
accuracy-based signal. Our phase-1 signal (iteration count) is quantized and noisier, and the
legacy-equivalent gain is kI = 0.5. Draft default: **PI with kI = 0.30, kP = 0.15, kD = 0** —
derivative action on a quantized count amplifies noise (Söderlind 2003's argument for
low-order filters), so kD exists as a member but ships 0.

## 3. Draft Implementation (jury: audit this section)

### 3.1 New members (`cl_FEM_Controller.hpp`, next to the existing cut/fail state
`mSolverFailCount` / `mDeltaTimeMin`, ~line 216-232 — audit R6: keep control-mode + cut state
together, not in the omega/Anderson block)

```cpp
//! timestep-control mode: true restores the legacy behaviour
//! ( iteration-target sqrt rule + memoryless halving ). Internal
//! switch by design — no input key until the PID mode is validated
bool mUseLegacyTimestepControl = false ;

//! PID gains on the cost error e = N_iterations / N_target, applied
//! multiplicatively ( log-space PID, Valli et al. 2002, CNM 18:131 ).
//! Legacy rule == pure I with kI = 0.5. Derivative ships disabled:
//! the iteration count is quantized and D-action amplifies its noise
real mCtrlKi = 0.30 ;
real mCtrlKp = 0.15 ;
real mCtrlKd = 0.0 ;

//! cost-error history: e_n, e_{n-1}, e_{n-2} of ACCEPTED steps only;
//! 1.0 = on target. Reset to 1.0 on every timestep cut ( a rejected
//! attempt's cost carries no information about the next accepted one )
real mCtrlErr0 = 1.0 ;
real mCtrlErr1 = 1.0 ;
real mCtrlErr2 = 1.0 ;

//! post-rejection growth hold: number of accepted steps that must
//! pass after a cut before the controller may grow the timestep
//! again ( SUNDIALS-style eta cap; prevents re-climbing the cliff )
uint mPostFailureHold      = 0 ;
uint mPostFailureHoldSteps = 2 ;

```

~~consecutive-cuts abort counter (`mConsecutiveCuts` / `mMaxConsecutiveCuts = 25`)~~ — DROPPED
after audit (D3): it could abort while Δt was still above the floor, i.e. while the retry policy
still had room, violating the expected-algorithmic-failure tier
(`doc/coding_philosophy.md:621-624`). The floor abort is the sole new terminal condition;
persistent singular systems keep their existing `mSolverFailCount` abort at 8.

### 3.2 `adjust_timestep()` — mode-switched φ only (`cl_FEM_Controller.cpp:1792`)

Replace the single line `real tPhi = std::sqrt( ... )` (cpp:1802) with:

```cpp
real tPhi ;
if ( mUseLegacyTimestepControl )
{
    // legacy: memoryless I-control with gain 0.5 ( pre-2026-08 rule )
    tPhi = std::sqrt( static_cast< real >( mIterationTarget ) / mIteration ) ;
}
else
{
    // sanity: both counts must be positive here ( guarded by the
    // mIteration0 gate in practice — audit O-B )
    BELFEM_ASSERT( mIteration > 0 && mIterationTarget > 0,
        "adjust_timestep called with a zero iteration count or target" );

    // shift the cost-error history and load this step's sample;
    // the max() floors a degenerate zero sample so it cannot poison
    // the P/D terms of the following steps ( audit A1 )
    mCtrlErr2 = mCtrlErr1 ;
    mCtrlErr1 = mCtrlErr0 ;
    mCtrlErr0 = static_cast< real >( std::max( mIteration, 1u ) )
              / static_cast< real >( mIterationTarget ) ;

    // multiplicative PID ( Valli et al. 2002 ): phi =
    // (e1/e0)^kP * (1/e0)^kI * (e1^2/(e0*e2))^kD
    tPhi = std::pow( mCtrlErr1 / mCtrlErr0, mCtrlKp )
         * std::pow( 1.0 / mCtrlErr0, mCtrlKi )
         * std::pow( ( mCtrlErr1 * mCtrlErr1 )
                     / ( mCtrlErr0 * mCtrlErr2 ), mCtrlKd );

    // a recent cut vetoes growth until the hold expires
    if ( mPostFailureHold > 0 )
    {
        tPhi = std::min( tPhi, 1.0 ) ;
        --mPostFailureHold ;
    }
}
```

Everything downstream (order-aware `tPhiMax` clamp, coupling-factor adaptation, min/max clamp,
save-grid trimming) is unchanged and shared by both modes.

### 3.3 `reset_timestep()` — floor enforcement, both modes (`cl_FEM_Controller.cpp:1746`)

Replace `mDeltaTime *= 0.5 ;` with (floor test reworked per audit O-F — "this cut would not
change Δt" is the abort condition; the run gets exactly one attempt AT the floor):

```cpp
// the minimum timestep is a hard floor in BOTH control modes. The
// halving lands on the floor once ( one attempt runs AT the minimum );
// a further cut demanded there cannot change Delta t, so the retry
// policy is exhausted and the run must fail loudly with the actual
// diagnosis instead of grinding ( sidecoatings collapse: a
// Delta-t-independent residual floor above the tolerance )
if ( 0.5 * mDeltaTime < mDeltaTimeMin )
{
    BELFEM_ERROR( mDeltaTime > mDeltaTimeMin,
        "timestep cut requested at the configured minimum ( %.3e s ):\n"
        "the nonlinear residual is not responding to the timestep size.\n"
        "Check the achievable residual floor against the tolerance\n"
        "( conditioning ), the time integration scheme, and the tape model.",
        ( double ) mDeltaTimeMin );
    mDeltaTime = mDeltaTimeMin ;
}
else
{
    mDeltaTime *= 0.5 ;
}

if ( ! mUseLegacyTimestepControl )
{
    // the failed attempt's cost history is meaningless for the next
    // accepted step, and stale samples would kick phi through the
    // P/D terms; hold growth for the following accepted steps
    mCtrlErr0 = 1.0 ;
    mCtrlErr1 = 1.0 ;
    mCtrlErr2 = 1.0 ;
    mPostFailureHold = mPostFailureHoldSteps ;
}
```

### 3.4 `finalize()` — no counter to clear

~~clear `mConsecutiveCuts` next to `mSolverFailCount`~~ — obsolete with D3 (counter dropped).
No `finalize` change remains in this plan.

### 3.5 Warm-start Δt persistence (`save_memdump` / `load_memdump`, cpp:2765/2814)

Audit-corrected (R1 compile blocker: `HDF5::active_group()` returns `hid_t` **by value**
(`cl_HDF5.hpp:127`) while the scalar helpers take `hid_t &` (`hdf5_tools.hpp:92,159`) — a
temporary cannot bind; an lvalue is required. Also: read while the `meta` group is still
selected, and declare a local `herr_t`).

In `save_memdump`, inside the `meta` group, after `mMesh->save_meta(...)` and before
`close_active_group()`:

```cpp
hid_t  tGroup  = tFile.active_group() ;
herr_t tStatus = 0 ;
hdf5::save_scalar_to_file( tGroup, "delta_time", mDeltaTime, tStatus );
```

In `load_memdump` on rank 0, after `load_meta` and BEFORE `close_active_group()` of the `meta`
group (old dumps lack the key — backward compatible):

```cpp
hid_t tGroup = tFile.active_group() ;
if ( hdf5::dataset_exists( tGroup, "delta_time" ) )
{
    herr_t tStatus = 0 ;
    hdf5::load_scalar_from_file( tGroup, "delta_time", mDeltaTime, tStatus );
}
```

All-ranks collective shape (audit O-E): inside the existing `if ( tExists )` block that already
broadcasts `mRunningTimeStep` (cpp:2868-2875), executed by every rank:

```cpp
broadcast( mDeltaTime );
mDeltaTime = std::clamp( mDeltaTime, mDeltaTimeMin, mDeltaTimeMax );
```

No HDF5 call on non-root ranks. The post-`finalize` save ordering already persists the
NEXT step's Δt (audit O-D — correct warm-start quantity, keep).

(Helpers: `hdf5::dataset_exists` defined at `hdf5_tools.hpp:35-44`; scalar save/load —
`hdf5_tools.hpp:92,159` — used the same way in `cl_Mesh.cpp:2995-3006`.)

Note the BDF history *spacing* (`mH`, `mStepCount`) is deliberately NOT persisted (O5): on
restart the startup ramp re-enters at BDF1 and rebuilds order, which is exactly the safe
behaviour observed to work. Restoring Δt removes the cliff; restoring `mH` would add IWG API
surface for no stability gain.

## 4. Ordered Steps

- [x] **R1** — Jury audit of this draft (slug `review_pid_timestep_draft`, Codex + Grok,
  2026-08-09); findings verified and reconciled into §6, must-fixes folded into §3.
- [x] **R2** — Christian's execute order 2026-08-09: implement with draft defaults — O2 floor
  abort (the round-1 C5 absolute-tolerance wiring provides the deck-level accept-at-floor
  escape hatch), O3 PID default, O1 draft gains pending R5. Christian set the sidecoatings deck
  to bdf1 / anderson off for the scheme A/B.
- [x] **R3** — Implemented 2026-08-09 (Claude, Christian's execute order) in
  `cl_FEM_Controller.hpp/.cpp`, together with the round-1 P1 package (neutral band, stall band
  0.2 dB, magnetic absolute tolerance live, watchdog promotion re-anchor, setup validation).
  Not compiled — build is Christian's.
- [x] **R4** — Implemented 2026-08-09: `delta_time` scalar in the memdump meta group, lvalue
  `hid_t` + local `herr_t`, read inside the open group, all-ranks broadcast + clamp.
- [ ] **R5** — A/B validation on sidecoatings: legacy vs PID at bdf1/bdf2/bdf5. Acceptance
  criterion per audit O-G: **fewer cut cascades, or a clean diagnostic abort at the floor — not
  "sidecoatings reaches tolerance"** (the Δt-independent residual floor is a separate disease,
  round-1 S4). Christian runs; per standing rule no builds from Claude. *(after: R3, R4)*
- [x] **R6** — Done 2026-08-09: round-3 jury on the controller diff (6 correctness findings,
  all fixed same turn — §6 D7-D9); `src/fem/doc/timestepping_strategy.md` (+ README) written,
  Codex prose polish applied; `doc/input_file_reference.md` + `nonlinear_controller_theory.md`
  synced (Codex caught the cross-doc drift); devlog `dl20260809_pid_timestep_controller.md`.
  Deck exposure of the PID keys stays conditional on R5 and spins off as its own todo then.
  ~~Note for the next anchor sweep: the §4.2/§4.4 Site cells in the input reference are
  table-wide drifted (pre-existing) and were left on the table's old baseline.~~ Done in R8:
  §4.2/§4.3/§4.4 anchors re-baselined against the current `set_params` (they were ~300 lines
  stale).
- [x] **R7** — Batch 2 (Christian 2026-08-10), implemented same day: (a) O2 revision — floor
  abort replaced by bounded iteration-budget escalation with warning box and deck-value
  restore in `finalize`; (b) warm restart made VISIBLE and controllable — opt-OUT per
  Christian's final ruling: `timestep { restart : ... ; }` defaults true (resume, announced
  by a WARM RESTART banner with resume t/timestep/Δt — previously silent,
  `hphirun.cpp:83` / `hphiTrun.cpp:90` load unconditionally), `false` ignores the dump with
  a log notice and starts fresh; (c) new deck key
  `solver { compute conditioning : true ; }` — the member existed but the constructor-time
  arming ran before `set_params`, so it was unreachable from a deck; the arming now happens
  at parse time with the null-solver guard (round-1 Grok L2); (d) the controller box header
  shows the configured scheme (BDF1…BDF5, Explicit, CN, Galerkin) instead of "Timestep".
  Docs synced same turn (input reference §4 direct-keys table + §4.4 rows, strategy doc
  §4/§5). Round-4 jury on the incremental diff: closed — 2 confirmed doc/comment fixes
  applied, the unbounded-floor-retry dissent routed to Christian, and findings in his
  concurrent DR-02 dump code reported untouched. *(not compiled — build is Christian's)*
- [x] **R10** — Conditioning diagnostic reworked (Christian 2026-08-10): κ sampled once per
  timestep at the **first iterate** instead of leaving the MUMPS error analysis armed for
  every solve; thermal κ added; `n/a` where a solver cannot supply it cheaply; footer closing
  row now carries `κM`/`κT` beside the timing in Christian's 22+48 layout (verified to render
  at exactly 71 columns, compact `1.43e8` formatter, no `%-Ns` padding of the two-byte κ).
  Docs synced same turn.
- [x] **R9** — Floor backstop (Christian 2026-08-10, accepting the round-4 jury's middle
  ground): new deck key `timestep { floor retries : n ; }`, default 20, `0` = unlimited
  (the pre-backstop behaviour). The escalation is unchanged; on exceeding the cap the run
  stops with a diagnosis naming the two deck-side remedies. The warning box now counts
  retries down. Docs synced same turn (input reference §4.4 + §4.5 guidance, strategy doc
  §4). Closes the only P1 both round-4 auditors raised.
- [x] **R8** — Input-reference extension (Christian 2026-08-10): the new keys documented in
  full (`solver { compute conditioning }`, `timestep { restart }`), every row this campaign
  touched revised, the internal (non-deck) controller constants listed explicitly, §4.5 field
  guidance for a collapsing Δt with the κ(A) discriminator, §13 restart pitfall, §2 memdump
  note, and the §4.2–§4.4 `Site` anchors re-baselined (~300 lines stale). Codex prose pass on
  the touched sections per the standing sync rule.

## 5. Open Questions

- **O1** — Gains for the quantized iteration-count signal. Draft: PI (0.15, 0.30, 0). The
  legacy-equivalent point (0, 0.5, 0) is one A/B lever. Needs empirical calibration on
  sidecoatings + one healthy case (no floor), R5.
- **O2** — Floor policy. ~~RESOLVED 2026-08-09 → abort~~ ~~REVISED 2026-08-10 → never abort~~
  **FINAL 2026-08-10: escalate, then a generous backstop** (Christian accepted the round-4
  jury's middle ground). A cut demanded at `mDeltaTimeMin` doubles the iteration budgets
  (`max iterations` + `watchdog window`, both fields) per floor retry, capped at 4× the deck
  values, prints a warning box counting the retries down, and retries at the minimum; budgets
  restore on the next accepted step (`finalize`). After `floor retries` consecutive attempts
  (deck key, default 20, `0` = the unlimited behaviour) the run stops with a diagnosis
  pointing at `compute conditioning` and `absolute tolerance`. The `absolute tolerance`
  escape remains the sanctioned accept for a Δt-independent floor. Implemented R7 + R9.
- **O3** — Default mode. RESOLVED 2026-08-09 → PID default, legacy via internal switch
  (Christian's instruction). Gregory's decks inherit PID; the R5 A/B and the legacy switch are
  the rollback path.
- **O4** — Phase 2 signal: predictor–corrector LTE (DASSL-style, needs IWG_Timestep history
  access) or WRMS solution change. Deferred; new todo when phase 1 validates.
- **O5** — Persist BDF `mH`/`mStepCount` in the memdump? Draft says no (startup ramp re-anchors
  safely). Revisit only if restarted runs show first-steps instability at restored Δt.
- **O6** — Warm-start fidelity is partial (audit O-C): `mCouplingFactor` / `mDeltaTime2` are
  not persisted (thermal Δt is rebuilt as `mDeltaTime/mCouplingFactor` at deck value). Scoped
  out; revisit with O4.

## 6. Defects (from the jury audit 2026-08-09, all verified by Claude against source)

- [x] **D1 (CRITICAL, compile)** — §3.5 as first drafted could not compile: `active_group()`
  returns `hid_t` by value (`cl_HDF5.hpp:127`), the scalar helpers take `hid_t &`
  (`hdf5_tools.hpp:92,159`); also missing `herr_t tStatus` and the read had to happen before
  `close_active_group()`. Found by Grok (R1) + Codex (#1), independently, blind. Fixed
  2026-08-09 in §3.5 (lvalue `tGroup`, local status, read inside open group).
- [x] **D2 (HIGH, spec honesty)** — "legacy bit-identical except floor" was false with the
  consecutive-cuts abort shared across modes; and the (0, 0.5, 0) A/B lever silently included
  the PID hold. Found by Codex (#2, #3) + Grok (R2), blind agreement. Fixed 2026-08-09: D3
  drop + §2 caveat + §7 rewrite.
- [x] **D3 (HIGH, policy)** — consecutive-cuts abort (25) could fire while Δt was still above
  the floor — an abort inside a live retry policy, against
  `doc/coding_philosophy.md:621-624`; it also overlapped the existing `mSolverFailCount`=8
  singular-abort (Grok R5). Fixed 2026-08-09: counter dropped entirely; floor abort is the sole
  new terminal condition.
- [x] **D4 (MEDIUM)** — hold semantics mis-stated: the first accepted step after a cut skips
  `adjust_timestep` entirely (`mIteration0` gate), so "2 accepted steps" under-counted the
  forced skip. Found by Grok (R3), anticipated in pre-registration A7. Fixed 2026-08-09: §2
  documents the true sequence (no re-growth for ≥3 accepted steps).
- [x] **D5 (MEDIUM)** — degenerate `e₀ = 0` sample could poison the P/D history (pre-reg A1;
  Grok O-B adds the `BELFEM_ASSERT`). Fixed 2026-08-09 in §3.2 (`max(N,1)` floor + assert).
- [x] **D6 (LOW)** — helper citations pointed at call sites, not the definition
  (`hdf5_tools.hpp:35-44` is the API); floor test's `BELFEM_EPSILON` band was microscopic
  relative to 100 ns and is replaced by the O-F "cut would not change Δt" test; member
  placement moved next to `mSolverFailCount`. Found by Codex (minor) + Grok (R4, O-F, R6).
  Fixed 2026-08-09 in §3.1/§3.3/§3.5.
- **Retained as documented behaviour:** BELFEM_ERROR side-effect style was already flagged in
  pre-registration (A2) and the increment no longer exists after D3.

Round 3 (implementation diff, jury 2026-08-09, thread `review_pid_timestep_impl`):

- [x] **D7 (HIGH on the opt-in path)** — after exhausted Newton backtracking, `mEpsilon` was
  restored to `mEpsilon0` but `mEpsilonAbs` kept the rejected trial's value, so an
  absolute-tolerance deck could decide on a state that was rolled back. Found by Codex
  (anticipated as pre-reg B1). Fixed 2026-08-09: `tEpsilonAbs0` captured beside `mEpsilon0`
  and restored in the exhausted-budget exit (verified by Claude).
- [x] **D8 (MEDIUM)** — consistency gaps of the new absolute criterion and the promotion
  re-anchor: the thermal flat-stall gate tested only the relative magnetic target (Codex),
  and the thermal Picard→Newton promotion did not re-anchor the thermal watchdog (Grok).
  Fixed 2026-08-09 at all three sites.
- [x] **D9 (MEDIUM)** — input/restart hygiene: memdump `delta_time` not NaN-hardened (helper
  asserts are debug-only), negative `absolute tolerance` silently disabling the escape, and
  the zero-iteration divide surviving on the legacy branch. Found by Grok. Fixed 2026-08-09:
  BELFEM_ERROR after the rank-0 load, ≥ 0 checks at both parse sites, BELFEM_ASSERT hoisted
  above the mode branch.
- **Policy findings (no code change, Christian's execute order):** PID default and the 0.2 dB
  stall band flip behaviour for all decks without deck keys (Grok) — the internal legacy
  switch, the (0, 0.5, 0) gain point, and the `stall tolerance` key are the rollback levers;
  both are top R5 watch items.

## 7. Definition of Done

Both modes compile and run the sidecoatings deck. Legacy mode preserves today's controller
behaviour with exactly ONE deliberate exception, shared by both modes: the `minimum timestep`
floor is enforced with a one-attempt-at-floor clamp and a diagnostic abort when a further cut is
demanded there (§3.3 — the round-1 C1 bug fix). A cut cascade therefore terminates in a
diagnostic abort instead of a sub-minimum grind; a warm start resumes at the dumped Δt;
R5 A/B recorded in the devlog.

## 8. Audit Trail

- `tmp/ai_exchange/review_fem_controller_timestep.md` — round 1 (controller + settings, jury,
  closed with reconciliation 2026-08-09).
- `tmp/ai_exchange/review_pid_timestep_draft.md` — round 2 (this draft, jury, closed with
  reconciliation 2026-08-09).
- Round 3 (post-implementation diff audit) — see §4 R6. Distill all before sweep.
