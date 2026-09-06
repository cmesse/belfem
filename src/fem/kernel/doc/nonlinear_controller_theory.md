# The Nonlinear Controller: Hybrid Iteration, Relaxation, and Timestep Control {#fem_kernel_nonlinear_controller_theory}

**Date:** 2026-08-06
**Purpose:** Theory and implementation reference for `fem::Controller` — the hybrid
Picard/Newton strategy, the adaptive relaxation, the safety nets (line search, stall guard,
progress watchdog, escalation, flip latch), and the adaptive timestep loop for coupled
magneto-thermal runs.
**Module:** `src/fem/kernel` (`cl_FEM_Controller`)

---

## 1. Role and structure

`Controller` owns the transient outer loop of the h-φ and coupled thermal solvers. It
starts each timestep, drives the nonlinear iteration of one or two kernels, decides
acceptance, adapts relaxation and step size, and retries with a halved step when an attempt
fails. Each outer iteration of a fully coupled run solves the magnetic system and then the
thermal system in Gauss–Seidel fashion. Since 2026-08-29 the complete nonlinear loop of a
timestep lives in the controller: drivers call `solve_coupled()` (or `solve_magnetic()` /
`solve_thermal()` for segregated runs with thermal sub-stepping) once per timestep and retry
on `reset()`. Each trip of the loop first **head-measures** the committed state — a fresh
assembly plus residual, no solve — and only then decides whether to take a body solve. A
timestep exits when each field's head certificate meets either its relative tolerance or its
optional absolute tolerance, with at least `min iterations` body solves taken: the state a
timestep commits is always the state whose residual was measured (the certified
exit). A thermal field accepted at a latched flat-stall is the one deliberate exception, and
the exit line labels it as uncertified.
Since 2026-08-27, a spent thermal iteration budget can cut the attempt (Sections 4 and 6).

The residual is the relative linear residual ε = ‖A x − b‖ / ‖b‖ of each kernel's system
(Messe et al. 2023, Eqs. 10–11). It is reported per iteration and broadcast so every MPI rank takes
identical control-flow decisions. Both Picard paths, plain and Anderson-mixed, evaluate ε
**pre-update**: it is the honest out-of-balance force residual of the previously accepted
iterate under its own assembly, which is the published convergence criterion. The Anderson
path must not test the updated iterate against the frozen assembly operator. At ω = 1 the
mixed update solves that system, so the "residual" degenerates to solver roundoff and
reports false convergence, observed as thin-shell checkerboarding under an
instantly-"converged" log. Mixing quality is available separately as the fixed-point
residual ‖G(x)−x‖/‖x‖ (`DofManager::fixed_point_residual()`). Bathe §8.4.4 warns that
increment-based measures alone can under-report true error on stiff problems; ε remains the
gate. Only the Newton path reports the freshly updated iterate. Its tangent differs from the
solved operator, so the recompute is non-degenerate; its ε therefore leads Picard ε by one
iterate at the switch thresholds. The certified exit does not change any of this: within the
loop every residual keeps its historical meaning, and always-accept Picard stays exactly as
published — the head certificate is an **additional** measurement of the committed state
that gates only the exit.

## 2. Hybrid Picard → Newton iteration

Each timestep starts on Picard (Messe et al. 2023, Section 4): the relaxed fixed-point step is robust
far from the solution but linearly convergent. Once the residual falls below the switch
tolerance (`tolerance switch`, default 1e-4 magnetic / 1e-3 thermal), the iteration is
promoted to the configured terminal scheme, typically quasi-Newton (Messe et al. 2023, Eq. 13); it is
demoted back if the residual rises above the switch again. Each scheme keeps its **own**
relaxation store (ω_Picard, ω_Newton); a promotion or demotion carries the **active** store's
clamped value across the switch, so the incoming scheme starts from the relaxation level
reached by the outgoing scheme.

The **first** promotion of an attempt enters Newton at half the carried relaxation
(ω_Newton = 0.5 · ω_Picard): the first full step after a promotion routinely overshoots from
a nearly-converged iterate, and the growth rule (Section 3) earns the other half back once
the tangent proves contractive. Only the first flip is damped — re-damping every
promote/demote cycle would ratchet ω down by ×0.5 per cycle. The same first-flip damping
applies to the thermal promotion. A timestep reset re-arms the latch.

Three latches protect against switch pathologies:

- **Stagnation fallback:** if a Newton phase stalls (flat residual window, Section 4), the
  solver falls back to Picard for the rest of the attempt (`mJustPicard`). After a
  promotion-Newton stall, Picard resumes from **its own** last relaxation at half throttle
  (0.5 · ω_Picard): the collapsed ω of a stalled Newton reflects the tangent's failure, not
  Picard's contraction history. After an *escalated*-Newton stall the stale Picard ω predates
  a genuine Picard breakdown, so the conservative min-merge of both stores is kept there.
- **Line-search fallback (added 2026-08-10):** a *promoted* Newton tangent that rejects all
  eight backtracking trials — omega 0.5 down to 0.5/2⁷ — has produced a step direction that
  is unusable at any scale, which is a tangent problem, not a timestep problem (observed
  signature: a net-transport-current deck with a point bearing, where the tangent carries a
  near-null φ-gauge mode; cutting Δt made it worse because promotion re-fired at every step
  size). The exhausted search restores the last accepted iterate and sets the same
  `mJustPicard` latch: Picard, which was contracting when it promoted, finishes the timestep
  at half throttle, with the stall window cleared and the progress watchdog re-anchored to
  the restored state. The *escalated* tangent (`mForceNewton`) is exempt — there Picard has
  already broken down, so the timestep cut remains the only lever. While either latch holds,
  the promotion is not even attempted: attempting it would re-anchor the watchdog's
  best-residual clock on every iterate and disable the guard that has to end a fruitless
  Picard grind.
- **Thermal flip latch:** in the coupled loop a running thermal Newton is kept until the
  residual retreats a full decade above the switch; the segregated thermal path has no such
  band. A residual that still orbits the switch flips the scheme every few iterations while
  each Newton entry degrades the residual; the third demotion within one attempt pins the
  thermal solver to Picard for the
  rest of that attempt (and, like the magnetic latches, suppresses further promotion
  attempts while it holds).

## 3. Adaptive relaxation

After every iteration the relaxation of the active scheme is adapted from the residual
history (Messe et al. 2023, Eq. 14, with α = 0.5, β = 1.1, γ = 0.4):

```
ω ← ω · α                                          if ε_{k+1} ≥ 1.05 ε_k          (decay)
ω held                                             if ε_k ≤ ε_{k+1} < 1.05 ε_k   (noise band, 0.2 dB)
ω ← ω · [ β + γ (2/π) arctan( (ε_k − ε_{k+1}) / ε_k ) ]   if ε_{k+1} < ε_k     (growth, ≤ ×1.3)
```

**Sign note.** The growth branch as implemented uses the argument (ε_k − ε_{k+1})/ε_k ∈ (0, 1]
— faster improvement grows ω more, up to exactly β + γ/2 = 1.3. Equation 14 as printed in
Messe et al. 2023 has the opposite sign in the arctan argument, which would *shrink* ω after fast
improvement (factor down to 0.9); the paper's own prose ("β > 1 and γ … enable rapid
convergence") supports the implemented behavior, and under the ω ≤ 1 clamp the printed
variant would make γ work *against* recovery on precisely the best steps. The implementation
is authoritative; treat the printed equation as an erratum.

Ten consecutive diverging iterations reset ω to its initial value — a deliberate last-resort
reset that either re-enters the convergence basin or reaches the divergence rules below
quickly. The growth branch stays active when Anderson mixing is enabled; see
`anderson_acceleration_theory.md` for why the mixing parameter must not be starved.

**Newton trust growth (coupled magnetic path).** Two *consecutive* Newton steps accepted on
the first trial (no backtracks) that each at least halve the residual have earned trust: on
the second one ω grows by ×2 instead of the ≤ ×1.3 adaptation crawl (a single such step is
not evidence near a residual floor), so a healthily contracting Newton recovers from a damped
entry geometrically. The test uses the algorithm that actually ran the iteration's solves, so a
same-call escalation cannot misclassify an accepted Picard step as Newton trust.

**Backtracking line search (coupled magnetic path, NEWTON solves only).** Picard iterates
always accept and follow the published Eq. 14 strategy: regression is answered by the
adaptation tail's α damping plus the divergence/stagnation/watchdog guards, never by
intra-iterate rejects. This is forced by the residual semantics: the pre-update ε measures
the iterate's **entry** state — the very state the backup holds — and the lagged post-update
residual of a relaxed Picard step is identically (1−ω)·r, pure ω arithmetic. A Picard reject
loop therefore re-solves the same system blind, restores the state it just condemned, and
can only exit by copying the stale residual — the frozen bit-identical residual print with
ω collapsing to the floor that stalled a run at a flux front, cutting Δt to zero without
ever changing the state. For Newton solves the post-update residual is non-degenerate
(J ≠ A), and the line search stands: each trial is accepted if it is within 0.3 decades
(~2×) of the reference — or, while the residual is still large (log₁₀ ε < 0.8,
where early/mid-settling wander is expected), if the regression stays within **one decade**
of the reference. The decade bound is what stops multi-decade kicks from a deep reference: a
full Newton step from a −50 dB iterate that jumps to −30 dB used to sail through the plain
absolute clause and could spiral the whole timestep.

**The reference is the pre-update residual of the CURRENT assembly** (since 2026-08-10), not
the previous iterate's ε. Those differ by one assembly: the previous ε was measured under
A(x_prev), while every trial is measured under A(x_entry) — and the ω → 0 limit of a trial is
exactly ‖A(x_entry)·x_entry − b(x_entry)‖, the pre-update residual. With the stale reference,
any state where the two assemblies disagree by more than the acceptance band rejected all
eight trials *flat in ω*, independent of the step direction (the greg5 net-current trace:
`Newton n, relax 0.00391`, step dead). Against the honest reference a small enough step is
always acceptable, so the search can only fail on genuinely explosive directions — which the
line-search fallback of Section 2 then catches. The pre-update residual is captured in the
Newton solve branch (`SolverData::mPreUpdateResidual`, broadcast with the other residual
scalars) and is identical across the trials of one iteration, since each retry reassembles at
the restored entry state. A rejected trial is rolled back (fields
and dof values restored), ω is halved, and the solve is retried, with a budget of 8
backtracks per outer iteration. An exhausted budget on a *promoted* tangent triggers the
line-search fallback of Section 2 (accept the restored iterate, latch Picard for the rest of
the timestep); on an *escalated* tangent it cuts the timestep past the grace period and
accepts the restored iterate inside it. Acceptance and rejection are also
the commit and discard points of the Anderson history (Picard accepts commit immediately).

**Moved-baseline detector (Newton rejects only, since Picard no longer rejects; its
original motivating case — a Picard reject storm under a moving thermal baseline — can no
longer occur).** A genuine step overshoot vanishes as ω → 0; a regression that
stays flat (within 0.05 decade) across **two** consecutive genuine ω halvings is
ω-independent, meaning the residual *baseline* moved under the magnetic system — near quench
the staggered thermal update shifts ρ(T) enough that the stale reference is a decade off at
the same state. The detector is armed **only when that premise can hold**: a thermal kernel
exists and the thermal update actually ran in the previous iteration (not frozen). On a
magnetic-only run an ω-independent flat regression is a defect signature (e.g. a wrong
constraint), and an ungated detector was observed accepting a 44 dB kick there. When armed,
such a trial is accepted so the iteration proceeds against the honest
current residual instead of burning the backtracking budget into a spurious timestep cut;
the divergence and stagnation guards still terminate real drift. The accepted trial's staged
Anderson pair is discarded and the window cleared (its fixed-point residual mixes two
thermal baselines), and the watchdog's best-residual tracker is restarted so it does not cut
the step the detector just rescued.

## 4. Safety nets

The safety nets run in this order: the divergence/ceiling group, the stall guard, and then
the watchdog. The coupled *thermal* watchdog is the exception. It runs immediately after
the thermal solve and therefore before the ceiling group. Escalation can intercept one
magnetic divergence verdict:

| Net | Trigger | Action |
|---|---|---|
| soft solver failure | the linear solver reports a failed factorization/solve (e.g. singular matrix at a strained iterate) under the soft-fail contract | cut Δt without touching the fields (the matrix was assembled at the *accepted* iterate, so relaxation cannot cure it); 8 consecutive failing attempts abort loudly |
| divergence rules | ε > 10 on **three consecutive iterates** past the minimum iterations, NaN, or `max iterations` exhausted | cut Δt — except a magnetic Picard breakdown with Newton available and untried, which escalates instead (next row). The persistence (strike counter, reset on any iterate at or below the bar and per attempt) gives the Eq.-14 α damping two halvings to catch a flux-front Picard overshoot: an instantaneous cut selected attempts by overshoot height, not trend, and pinned Δt at the front. The magnetic ceiling depends only on the iterate count. In coupled mode, the thermal ceiling counts completed thermal solves and has been enforced since 2026-08-27; it was inert before that date. Thermal solves run in lockstep with the magnetic iterates while the update gate permits, and frozen iterates do not count. The ceiling fires only when the thermal field is the reason the loop continues. The magnetic field must have reached its relative or absolute target at least once in the attempt, recorded by a per-attempt latch since 2026-08-28. Before that date, the magnetic field had to be at target on the same iterate. Continuing thermal updates could therefore perturb an already-converged magnetic field and disarm the ceiling. The thermal field must also have met neither target and must be neither stall-latched nor currently flattening toward the flat-stall exit. The acceptance test still requires both fields under target simultaneously, so the latch can never accept a degraded magnetic iterate — it can only stop paying for thermal updates that keep spoiling one. An unconditional thermal ceiling would cap the entire coupled loop at the thermal budget while both fields were still converging. The ordering above has one side effect. If the thermal watchdog and thermal ceiling are both due on the same iterate, the watchdog wins because it is evaluated first and returns. Across different iterates, whichever predicate first becomes eligible ends the attempt. The ceiling becomes eligible only after the magnetic field first reaches target. |
| escalation | Picard flooring far above the switch tolerance (the lagged-conductivity map lost contraction through the E–J transition), terminal scheme is Newton, once per attempt, **and ε < 10** | force Newton for the rest of the attempt before cutting Δt |
| stall guard | mean absolute deviation of the last `stall window` residuals (in dB) below `stall tolerance` while above tolerance | Newton → Picard fallback; if Picard itself stalls → cut Δt |
| progress watchdog | **no new residual minimum** for `watchdog window` iterations (default 30) while still above 10× tolerance, except when that field's relaxation grew since the previous watchdog evaluation; that growth is evidence of a line search still recovering from an overshoot | cut Δt. In coupled mode, the effective *thermal* window has been clamped at setup since 2026-08-28. The bound is the magnetic `max iterations` minus the thermal `min iterations` minus one, floored at 1. With lockstep counting, the magnetic ceiling is the last iterate that can run the thermal watchdog, so a wider window could never fire first. The clamp accommodates an early thermal best; a late best or the spare clauses can still let the ceiling win. `0` stays the disable, and segregated decks are never clamped. |
| thermal flat-stall exit | magnetic converged, thermal residual bit-flat (< 0.001 decade change) for 5 consecutive coupled iterations while above tolerance | accept the timestep with a warning box — further iterations provably change nothing (typically T pinned at a material table ceiling, which zeroes the dT-derivatives) |

The escalation's ε < 10 guard exists because escalation is a rescue for Picard stalling at a
*low* floor, not for divergence — and since escalation sets the force-Newton latch, which
disables the +10 dB divergence cut, escalating from a blown-up iterate would leave Newton
flailing with no exit.

**Soft-fail contract.** Under the controller, the sparse solver wrappers (MUMPS, STRUMPACK)
record a failed factorization/solve instead of aborting (`set_soft_fail`); the failed LHS is
never written into dofs or fields, the failure verdict is rank-uniform (MUMPS propagates
errors to every rank; the STRUMPACK parallel path reduces the verdict with an
`MPI_Allreduce`, since STRUMPACK itself does not), and the controller treats the event as a
failed trial that cuts the timestep. Standalone solver users keep the loud abort. The
consecutive-failure counter clears on every completed timestep.

The watchdog exists because the stall guard's flat-band criterion cannot see two important
failure shapes: (a) *sawtooth limit cycles*, where the relaxation collapses and periodically
resets so the residual oscillates with locally-good slopes but makes no net progress, and
(b) *slow monotone creep* of the thermal residual at the relaxation floor, which never
violates the ε > 10 divergence rule. "No new minimum in N iterations" is immune to the shape
of the oscillation. The watchdog tracks each attempt's best residual per kernel; an
escalation restarts the clock against the current residual so the forced Newton phase gets a
full window to prove itself. Setting `watchdog window : 0` disables it.

The **relaxation-growth spare** exists because the bare no-new-best test cannot tell a stall
from a *recovery*. The adaptive relaxation halves ω on an overshoot and regrows it at
roughly β per accepted iterate, so climbing back from ω ≈ 0.1 to 1.0 takes ~17–22 iterates,
longer than a tight window. On the coarse tapestack3d mesh the bare test cut mid-recovery steps whose
best residual was within 1 dB of tolerance, and the cuts *cascaded*: a half-Δt retry
reproduced the same overshoot/recovery shape and was cut again (measured 2.10 → 1.05 →
0.53 ms while the best residual *improved*). The spare reads "ω grew since the previous
iterate" as evidence that the line search is still working. A genuine stall has ω pinned or
shrinking (the floor grinds, and every replayed stall keeps firing), while a recovery keeps
growing ω even when the residual climb itself wobbles (in step 378 the cut iterate had ω growing and
the residual regressing, so the spare is an ω statement, not a residual one). Replayed against both campaign logs: the spare kept every omega-pinned
firing, released the mid-recovery cuts, and produced zero false positives on accepted steps.
Accepted cost: when a step's *other* field creeps inside its own 10× band while this
field's ω happens to grow, the step is no longer cut early. It runs to `max iterations`,
which is the pre-watchdog behavior and is still bounded. The test compares the ω that *produced*
the current residual (passed by the caller, post-line-search); it is never re-derived from
the live algorithm, which a same-call escalation or a stagnation demotion may already have
flipped. The comparison history resets only at attempt boundaries; promotions and escalations
restart the best-residual clock, but they do not invalidate the previously recorded ω.
The spare narrows the sawtooth coverage claimed above: an ω-sawtooth is spared on every
up-phase and can only fire on a down-phase or plateau iterate. A short cycle still fires
within a window or two; a *long* up-phase with no new best is exactly the recovery shape the
spare exists to protect. This is a deliberate trade, replayed but not proven. The "zero false positives"
sentence describes those two logs, not a guarantee for the next mesh.

**Retry hygiene.** A cut restores the savepoint state and restarts the attempt cleanly:
both equations return to Picard, the relaxation stores, divergence counters, flip counter,
residual histories, watchdog trackers, and the Anderson windows are all reset. The
relaxation carry across a cut is deliberate policy: after a divergence the under-relaxation
is kept (restoring full ω would re-diverge at the smaller step); after a clean accuracy cut
the initial ω is restored.

**The thermal update gate** (`update gate`, disabled by default) skips the thermal solve
while the magnetic residual is above a threshold. Freezing the thermal update can convert a
mutual divergence spiral into a limit cycle pinned at the gate rather than resolving it;
once the state is outside the convergence basin, the controller needs a timestep cut, not
update reordering. The knob is retained for experiments
only; when it fires, the thermal solver restarts on Picard at thaw, because a cold Newton
start from a state that moved during the freeze diverges.

## 5. Adaptive timestep control

After a converged step, the next step size targets a fixed iteration cost (Messe et al. 2023, Sec. 4):

```
Δt ← Δt · clamp( φ, 0.5, φ_max ),   φ = (e₁/e₀)^kP · (1/e₀)^kI · (e₁²/(e₀e₂))^kD,   e = N_iterations / N_target
```

bounded by `minimum timestep` / `maximum timestep`, with save-point alignment (`save every`)
temporarily overriding and restoring Δt. A lost attempt halves Δt and re-runs from the
savepoint. With a save grid configured, `save()` **verifies** that the converged time
actually lies on the grid instead of trusting the sticky flag armed when the step was
trimmed — a timestep cut skips the adjust phase, and the retry would otherwise store an
off-grid frame. The controller is a log-space PID on the cost error e (kP = 0.15, kI = 0.30,
kD = 0; Valli, Carey & Coutinho 2002); the legacy rule sqrt( target / used ) is its pure-I point
kI = 0.5 and is kept behind an internal switch. φ_max is order-aware — 1.5 at active BDF order
≤ 2, 1.4 at order 3, 1.2 at orders 4–5 — which keeps the step-size ratios inside the BDF
zero-stability bounds (see `bdf_timestepping_theory.md`); growth is also held at φ ≤ 1 for the
two accepted steps after a cut (see `src/fem/doc/timestepping_strategy.md`).

Two configuration warnings are worth calling out:

- A large `target iterations` value (e.g. 100) steers Δt to the edge of the convergence
  basin, so the controller favors a few very large, marginally convergent steps. Values
  near the default 20 usually give several cheap steps instead of one expensive failure.
- The nonlinear `tolerance` interacts with physics quality: Messe et al. 2023 requires ε ≈ 1e-11 in the
  Newton stage to prevent element-wise current-density checkerboarding in HTS simulations.
  Runs at looser tolerance should verify the current-density field visually.

## 6. Input keys (nonlinear sections)

| Key | Default | Meaning |
|---|---|---|
| `tolerance` (alias `relative tolerance`) | 1e-6 | relative convergence target |
| `absolute tolerance` | 0.0 (disabled) | opt-in absolute escape on the raw ‖Ax−b‖, which is dimensional — no universal default is possible. Live in BOTH sections since 2026-08-09: the loop terminates once the absolute residual is at or below the target. This is the sanctioned way to accept a run whose relative residual sits at a solver noise floor above `tolerance` (see `src/fem/doc/timestepping_strategy.md` §4) |
| `tolerance switch` | 1e-4 (magnetic), 1e-3 (thermal) | Picard → Newton promotion threshold |
| `algorithm` | Picard | terminal scheme: `Picard` or `Newton` |
| `min/max iterations` | 2 / 100 | iteration bounds for each attempt. In coupled mode, the thermal maximum counts completed thermal solves. These solves run in lockstep with the magnetic iterates while the update gate permits. A spent thermal budget cuts the attempt only when the magnetic field has reached its target at least once in the attempt, using a per-attempt latch since 2026-08-28, and the thermal field has not reached either of its own targets. Both keys are validated at setup since 2026-08-28 (positive maximum, nonnegative minimum, maximum at least minimum). See the divergence-rules row in Section 4 |
| `target iterations` (magnetic section only) | 20 | iteration cost the Δt adaptation aims for |
| `min/max relaxation` | 1e-3…1 (magnetic), 0.1…1 (thermal) | ω clamp |
| `stall window` / `stall tolerance` (magnetic section only) | 5 / 0.2 dB (default raised from 0.001 dB on 2026-08-09 — the band must sit above the noise of a floored residual or the Newton→Picard demotion never fires) | flat-band stall guard |
| `watchdog window` | 30 | no-new-minimum watchdog with the relaxation-growth spare; 0 disables |
| `anderson depth` | 0 (off) | Anderson mixing window; see `anderson_acceleration_theory.md`. If no explicit depth is given, `timestep { anderson stabilization : true ; }` fills defaults 3 (magnetic) / 1 (thermal). |
| `scheme` (timestep section; legacy alias `method`) | bdf1 | time integration: `bdf1`…`bdf5`, `explicit`, `crc`/`crank-nicolson`, `galerkin`. Crank-Nicolson and Galerkin hard-error with stiffness; see `bdf_timestepping_theory.md`. |
| `update gate` (thermal) | off | experimental freeze gate, see Section 4 |
| `coupling` (thermal) | fully coupled | `fully coupled` or `segregated` (+ `coupling factor`) |

## 7. Literature

- Messe et al. 2023, Section 4 — residual definition (Eqs. 10–11), hybrid
  Picard/Newton strategy (Eqs. 12–13), relaxation rule (Eq. 14, sign note in Section 3
  above), timestep adaptation, and the checkerboarding tolerance requirement.
- Arsenault et al. 2023 — the magnetodynamic h-φ coupling whose convergence
  behavior this controller manages.
- Bathe, Ch. 8 — general convergence and robustness strategies for nonlinear FE iteration.
