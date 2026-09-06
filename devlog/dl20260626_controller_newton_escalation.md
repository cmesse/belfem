# Controller: Picard-breakdown → Newton escalation

**Date:** 2026-06-26
**Purpose:** Stop the adaptive controller from cutting Δt indefinitely when a magnetic
Picard iteration breaks down at a residual far above `mEpsilonSwitch`; escalate to the
configured Newton tangent once before resetting the timestep.
**Module:** fem/kernel (`cl_FEM_Controller`)

## Problem (CORC, Δt_max = 50 ms)

Reported by Christian. At timestep 66 (t = 3050 → 3100 ms) the magnetic Picard iteration
takes one good step (46 → 11.8 dB) then stalls. The controller halves Δt and retries; every
retry stalls identically, the residual floor barely moving (11.83 dB @ 50 ms → 6.61 dB @
0.39 ms — ×128 in Δt for ~5 dB). Three-way analysis (Claude + Codex + Grok,
`tmp/ai_exchange/corc_controller_stall.md`) consensus:

- `residual 9.000000` is a display clamp (`std::min(mEpsilon, 9.0)`); the true value is the
  dB column. Floor ≈ ε 15, five decades above `mEpsilonSwitch` (1e-4).
- Root cause is a **state-local nonlinear-stiffening event** at the t≈3050 ms current-ramp
  advance (E-J power law, n=20, |J|→Jc), NOT a timestep-size problem. Halving Δt is the wrong
  remedy and the controller applied it blindly.
- **Core defect:** Newton never engages. The Picard→Newton handoff is gated behind
  `mEpsilon < mEpsilonSwitch`, which the stalled residual never reaches, so the solver stays
  pure Picard — and the lagged-conductivity Picard map is not contractive in this regime.

## Change

New `Controller::try_escalate_to_newton()` (cl_FEM_Controller.cpp). On a magnetic Picard
breakdown, if the configured terminal algorithm is Newton and we have not escalated yet this
attempt, it pins the solver to Newton (`mForceNewton`), restarts relaxation from `mOmega0`
(not the collapsed line-search omega), clears the stall window, and returns true so the caller
keeps iterating at the **current** Δt instead of cutting it. Capped at one escalation per
attempt (`mNewtonEscalated`) so it cannot ping-pong; if Newton also breaks down the caller
resets exactly as before. State machine mirrors the existing Newton→Picard fallback in
`magnetic_stagnation_forces_reset()`.

Wiring:
- New members `mForceNewton`, `mNewtonEscalated` (hpp); reset in `initialize_timestep()` and
  `initialize_magnetic()`.
- Handoff condition (both `iterate_coupled` and `iterate_magnetic`) now fires Newton when
  `mForceNewton || (mEpsilon < mEpsilonSwitch && mIteration > 1)`.
- The two post-loop reset sites (coupled, magnetic) call `try_escalate_to_newton()` before
  `reset_timestep()`; refactored the compound coupled condition into `tMagneticReset` /
  `tThermalReset` so escalation fires only on a magnetic (non-NaN, non-thermal, in-budget)
  breakdown — behavior-preserving otherwise.
- `magnetic_stagnation_forces_reset()` Picard branch returns `!try_escalate_to_newton()`;
  Newton-stall branch now also clears `mForceNewton`.

Only engages when `algorithm : Newton` is configured (CORC input.conf does). Pure-Picard runs
are unaffected.

### Follow-up: let the escalated Newton actually run (two guards)

First test (escalation only) showed the banner firing but Newton never iterating — the Δt
still halved every round. Two reset sites were guillotining Newton before it could work; both
now gated on `mForceNewton`:

1. **In-loop line-search guard** (`iterate_coupled`, the backtracking loop): a freshly
   escalated Newton no longer resets on its first over-relaxed overshoot
   (`tLogEpsilon > 0.9 && !mForceNewton`); it damps through the backtracking budget (ω halving)
   and only resets once the budget is exhausted (`tBacktracks >= 8`).
2. **Post-loop "+10 dB past grace" guard**: while Newton is escalated, the
   `mEpsilon > 1E1 && mIteration > mMinNumIterations` reset is suppressed
   (`tMagneticDiverged = !mForceNewton && …`) so Newton gets *several* damped iterations to
   pull the residual down. Termination falls back to the stagnation guard (flat residual hands
   Newton back to Picard, a second Picard stall then resets), the `mMaxNumIterations` ceiling,
   and NaN — all still active.

## Status — VERIFIED WORKING (2026-06-26)

- Builds clean: `make hphirun -j20` (mpicxx on PATH), `-Werror -pedantic-errors`, exit 0.
- Christian confirmed on his working CORC tree: at the breakdown timestep, after the
  "escalating to Newton" banner the solver now runs actual damped Newton iterations at the
  same Δt and carries the step instead of cascading Δt halvings. ("we're looking good!")
- (Claude's local repro stayed blocked by a pre-existing periodic-mesh setup crash on this
  branch — `PeriodicityFactory::map_facets`, facet 49349 has no partner — unrelated to the
  controller; the corc.msh in cmake-build-debug/corc/ does not pass periodic matching with the
  current bfmfile working tree. Verification was done on Christian's run.)

## Next (deferred)

- Controller-side robustness from the consensus writeup (`tmp/ai_exchange/corc_controller_stall.md`):
  `mDeltaTimeMin` floor in `reset_timestep()` and a "is halving even helping?" Δt-insensitivity
  detector, so a genuine structural wall can't drive unbounded bisection.
- The escalation in `iterate_magnetic` (run_magnetic path) still fires but isn't damped — that
  path has no line search, so its post-loop +10 dB reset still cuts it. Not exercised by the
  CORC hphirun case; revisit only if run_magnetic needs it.
