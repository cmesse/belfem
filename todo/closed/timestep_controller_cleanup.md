# Time-Step Controller Cleanup and Future Work

**Date:** 2026-06-26
**Purpose:** Record the resolution of the disabled error-controller `#HACK` in `Controller::adjust_timestep()` and park the larger controller/integrator ideas that are out of scope for now.
**Module:** fem/kernel (`cl_FEM_Controller.cpp`)
**Status:** CLOSED 2026-07-03 — the `#HACK` resolution is done and verified in the tree
(`adjust_timestep()`, `cl_FEM_Controller.cpp:1098`: clean `tPhi = sqrt(mIterationTarget/mIteration)`,
clamp `[0.5, 1.5]`, all dead error-proxy machinery removed). The "Out of scope — parked" and
"Future directions" sections below are **not** part of this task and remain on record here for the
historical trail (the two parked items — save-point restore of an unvalidated step, inert ω
line-search — are actionable bugs; revive as their own todo if/when scheduled).

## Background

`adjust_timestep()` carried a textbook error-controlled step-size formula
`dt *= (tol/err)^(1/(p+1))`, but the error estimate it was fed —
`tErr = mEpsilonFirst * mDeltaTime` — is invalid:

- `mEpsilonFirst` is a **normalized nonlinear residual** (`‖r‖/‖b‖`), i.e. a measure
  of *solver difficulty*, not of *temporal truncation error*.
- `tErr` mixes that dimensionless residual with a time in seconds, then compares it
  against `mRelativeEpsilonTarget = 1e-6` (a dimensionless residual tolerance) —
  dimensionally and in scale incoherent.

For realistic residuals and ms-scale steps the resulting factor is ≪ 1, so it would
relentlessly shrink the step toward `mDeltaTimeMin`. That is why it was overruled with
`tErrFactor = 1.0` (`//#HACK`), which returns the controller to (approximately) the
published iteration-count method (Messe et al. 2023, paper1, Section 4).

Worked example (observed step 8: `mEpsilonFirst = 0.0909`, `Δt = 10 ms`):
`tErr = 9.1e-4` → `tErrFactor = sqrt(1e-6/9.1e-4) = 0.033`; combined with the iteration
factor this *halves* a step that converged in 2 iterations — exactly backwards.

## Done

- [x] Remove the dead temporal-error-proxy machinery (`tErr_EM`, `tErr_Th`, `tErr`,
  `tErrFactor`, `tSafety`, `tOrder`, `tIterFactor`, `tIterRatio`) and the `//#HACK`.
- [x] Keep a single, clean **iteration-count** controller:
  `tPhi = sqrt(mIterationTarget / mIteration)`, sizing the step to hold the solver near
  its target iteration count (Messe et al. 2023, paper1, Section 4).
- [x] Tighten the growth clamp from `[0.5, 2.0]` to **`[0.5, 1.5]`** to match the paper's
  ×1.5 growth rule and curb over-aggressive jumps. *(Behavioral change — flagged for
  review; revert to 2.0 if undesired.)*

## Out of scope now — parked

Separate concerns observed while tracing the controller; not addressed in this cleanup.

- [ ] **Save-point restore reinstates an unvalidated step.** When a step is clamped to
  land on a save point, the original (larger) `Δt` is stashed in `mDeltaTimeTemporary`
  and restored verbatim next step (`cl_FEM_Controller.cpp:1040-1045`), bypassing the
  growth clamp — so the first *real* trial of that larger step can be a hard failure
  (observed 10 ms → restored 25 ms → cut to 12.5 ms). After a cut, do not blindly
  reinstate the parked step; re-enter it through the controller.
- [ ] **Inert ω backtracking line-search** in `iterate_coupled()` (`while(tRun)`,
  ~`:489-563`). The residual is measured at the **pre-solve** field state
  (`mFieldValues`, snapshotted at the end of `compute_jacobian_and_rhs()`), while `ω`
  only scales the **post-solve** update — so halving `ω` (1.0 → 0.0078 → 0.001) leaves
  the residual bit-for-bit identical and the search does nothing. The paper's relaxation
  rule (Eq. 14, across iterations) is the real mechanism and already works; either
  evaluate the residual at the relaxed iterate `u_new`, or remove the inner line-search.

## Future directions (design, not yet scheduled)

Beyond the published controller; recorded from the 2026-06-26 design discussion.

- [ ] **Dual controller:** `Δt = min(Δt_accuracy, Δt_convergence)` — the paper's
  iteration rule as a stability *floor*, plus an accuracy *ceiling* it never had.
- [ ] **Convergence signal = contraction rate** `ρ = ε_{k+1}/ε_k` (or Deuflhard's
  affine-invariant monotonicity test for the Newton phase) instead of quantized
  iteration count — smoother and *predictive* (shrink before failure, not after).
- [ ] **Feed-forward pre-shrink** around analytically-known source events (`|dI/dt|`,
  `|dB/dt|`, current reversals) to avoid the reactive failed-step-then-cut waste.
- [ ] **Integrator upgrade → TR-BDF2** (2nd order, L-stable, *free* embedded local-error
  estimate). Makes `Δt_accuracy` cheap. **Open risk:** verify the n≈30 power-law
  checkerboarding (the reason implicit Euler was chosen for maximal damping) still dies
  under the composite step. Coordinate with `bdf_jacobian_scaling_bug.md` (BDF
  coefficients / `mAlpha` are currently non-functional).
- [ ] **Pseudo-transient continuation (SER)** for DC / initialization phases, where time
  accuracy is not required and `Δt` can serve as the continuation parameter.

## Related

- [[nonlinear_iteration_strategy_near_quench.md]] — high iteration counts near quench as
  a temporal-discretization-error symptom → adaptive Δt with quench detection; the
  convergence-signal and feed-forward items above are the natural mechanism.
- [[bdf_jacobian_scaling_bug.md]] — BDF2-5 currently non-functional (`mAlpha` NaN);
  prerequisite for any higher-order integrator work.
