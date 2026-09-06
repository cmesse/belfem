# Kernel theory documentation: BDF, Anderson acceleration, nonlinear controller

**Date:** 2026-07-31 (overnight session)
**Purpose:** Write the theory reference documentation for the time integration and nonlinear
solution machinery: three new documents in `src/fem/kernel/doc/`, audited by Grok
(independent re-derivation of the mathematics) and Codex (claim-by-claim verification
against the source), followed by a Codex prose pass.
**Module:** src/fem/kernel (docs), src/fem/iwg (BDF subject matter)

## What landed

- `src/fem/kernel/doc/bdf_timestepping_theory.md` — variable-step BDF1–BDF5: coefficient
  construction from Lagrange differentiation with the fixed-step limits verified by
  independent derivation (BDF2: 3/2, 2, 1/2; BDF3: 11/6, 3, 3/2, 1/3), the assembled system
  shape and the centralized alternating-sign history combination, startup order ramp,
  step-size history with rejection/savepoint handling, and stability/accuracy notes
  (A(α)-stability, step-ratio bounds, and the 1/(1−λh) overestimation of runaway modes that
  makes step-size sensitivity studies mandatory for quench transients).
- `src/fem/kernel/doc/anderson_acceleration_theory.md` — type-II Anderson mixing of the
  Picard branch: formulation and its exact plain-step reduction at depth 0, the
  column-normalized QR least squares with its safeguard ladder, the β-asymmetry argument
  (β scales the residual content but not ΔXγ — hence flush-on-reject rather than a frozen ω),
  the stage/commit handshake and the five flush-rule classes, mixed-state residual
  semantics, cost, and the `anderson depth` configuration.
- `src/fem/kernel/doc/nonlinear_controller_theory.md` — the hybrid Picard→Newton strategy
  with per-scheme relaxation stores and switch latches, the relaxation adaptation rule with
  the Eq. 14 arctan **sign erratum** documented (the implemented sign matches the paper's
  prose and design intent; the printed equation would penalize the best steps), backtracking
  line search, the safety-net roster (divergence rules, escalation, stall guard, progress
  watchdog, retry hygiene, the experimental update gate), adaptive timestep control, and a
  complete input-key table for both nonlinear sections.
- Index entries in `src/fem/kernel/doc/README.md` (new "Theory References" section).

## Audit round

Grok re-derived the BDF coefficients, the stability claims, the Anderson formulation, and
the controller bounds independently: all verified, three wording tightenings applied. Codex
verified every factual claim against the source and found seven documentation defects, all
fixed: a residual-semantics overclaim, safeguard and safety-net ordering, an incomplete
acceptance rule, missing/overbroad input keys, the `reset_fields` history restoration scope,
and a step-ratio overclaim. A 13-item prose pass followed.

## Two code observations — FIXED same day (2026-07-31, approved by Christian)

1. The Newton branch's post-update residual recompute multiplied with the assembly-time
   `mFieldValues`, which the update loop did not refresh — the residual it reported belonged
   to the pre-update iterate, the same one-iteration lag the Anderson path had to fix for
   its promotion logic. **Fixed:** the Newton update loop now refreshes `mFieldValues`
   alongside the dof values, so the recompute (whose `mRhsBackup` machinery existed for
   exactly this purpose) finally reports the updated iterate's residual. The nonlinear
   fixed-point set is unchanged, but finite-tolerance exits, promote/demote timing, and the
   line search can accept different iterates — notably, the pre-fix Newton backtracking was
   residual-blind to ω (every trial re-measured the same pre-update residual), so this fix
   turns the Newton line search into a real residual test (Grok audit).
2. `IWG_Timestep::reset_fields()` restored step-size history slots `mH(0..2)` but not
   `mH(3)`, which BDF5 reads — and the value could not be recovered by un-shifting, because
   the shift rotation destroys it. **Fixed:** `shift_fields` keeps a one-deep backup
   (`mHDropped`) of the value falling off the end; `reset_fields` restores it, and the
   savepoint mechanism snapshots/restores it alongside `mH`.

Both fixes audited by Codex + Grok (exchange `tmp/ai_exchange/newton_lag_mh3_fixes.md`);
theory docs updated in place to describe the fixed behavior.
