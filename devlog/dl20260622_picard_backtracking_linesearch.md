# Backtracking line search for the magnetic Picard/Newton iteration

**Date:** 2026-06-22
**Purpose:** Tame the highly nonlinear HTS power-law iteration (Picard oscillation,
60–72 iters/timestep with residual blow-ups) by adding a revert-if-worse
backtracking line search to the relaxation loop.
**Module:** fem/kernel (`cl_FEM_Controller.cpp::iterate_coupled`)

## Problem

The h-φ magnetic solve showed pathological convergence: Picard residuals oscillating
for 60–72 iterations, periodic blow-ups (true residual ~5e4–4e7, shown clamped to
`9.000000` in the residual column while the dB column reports the true value), and a
relaxation factor thrashing 0.05–0.5. Christian hypothesised the no-resistivity buffer
layer was the cause.

## Literature consultation (excluding Messe 2023, per Christian — author distrusts own paper)

- **Buffer is a φ-region, not the culprit.** `cl_IWG_Maxwell.cpp:257-258` routes
  `DomainType::Air` and `DomainType::Buffer` to the *same* scalar-potential weak form
  (`phi_tri3`/`phi_tet4`), so the buffer is linear and carries no resistivity *by
  design* — correct h-φ practice (Arsenault et al. 2021, paper7, ~l.169: "in the air
  domains … no resistivity must be specified"). It cannot be the nonlinear driver.
- **The real cause is the power-law non-coercivity.** Dular et al. 2021 (paper0), §V
  (~l.1100-1112): with the power law the differential resistivity → 0 as J → 0, so the
  curl-curl operator loses coercivity → non-contractive fixed-point map → the
  oscillation/blow-ups. Their remedy is a "regularized version with two limiting
  resistivity values."
- **The globalization fix is backtracking.** Bathe §8.4 (line search) and Deuflhard's
  damped-Newton natural-monotonicity test: reject a step that worsens the residual,
  reduce damping, retry from the last good state.

## ρ_min analysis (concluded: floor is well-placed, NOT the lever)

BELFEM already implements Dular's two-limit regularization: floor `mRhoMin`/`gRhoMin`
= 1e-16 (`cl_Material.hpp`, `powerlaws.hpp` `std::max`), physical ceiling via parallel
combination with ρ_n, plus a `std::clamp(·, gRhoMin, gRhoMax=1e10)` in the assembly
(`mt_maxwell_h.hpp:35-36`). With the run's actual values from `cmake-build-debug/input.conf`
(`J_c=1e10 A/m²`, `n=20`, `E_c=1e-4 V/m`): ρ_c=E_c/J_c=1e-14; the floor engages at
**J*/J_c = 0.785**, just below the loss-relevant band, with spurious `E_floor < 1%` of
`E_c`. Only ~0.5 decade of safe headroom before clipping the sub-critical loss tail.
Conclusion: the floor is sized about right; the oscillation is the *intrinsic*
`ρ ∝ J^19` sensitivity (a 5% J error → ~95% ρ swing) over the [0.8, 1.1]·J_c operating
band, which no floor touches. So the high-value change is the iteration safeguard, not
retuning ρ_min.

## Change implemented

`Controller::iterate_coupled()` now wraps assemble → solve → residual in a backtracking
loop:
- **Snapshot before the step:** `mMesh->field_data(label)` for every dof field **and**
  every `dof->value()` (`mBackupFields`, `mBackupDofValues`). Backing up the dof values
  too is required — the relaxed update reads `dof->value()` as its `x_k` baseline
  (`cl_FEM_DofMgr_SolverData.cpp:2157`), which a field-only revert would miss.
- **Accept / reject:** baseline `tLogEpsilonRef = log10(mEpsilon0)` (the last accepted
  residual; `REAL_MAX` on the first iter of a timestep, so trial 1 auto-accepts). Accept
  when `tLogEpsilon < tLogEpsilonRef + 0.3` (tolerate ≤ ~2× worse — Picard is
  legitimately non-monotone); else revert both snapshots, `tOmega → clamp(0.5·tOmega,
  mOmegaMin, mOmegaMax)`, retry.
- **Termination:** self-terminates as ω → ω_min (update vanishes → residual returns to
  `mEpsilon0` → accepted); an 8-trial cap force-exits to the existing `mEpsilon > 1E1`
  → `reset_timestep()` Δt-cut as a backstop.

## Debugging arc (two dead-loop iterations before it worked)

The first two attempts were silent no-ops: `tLogEpsilon` was seeded to `BELFEM_REAL_MAX`
and copied into the comparison each trial, so `tLogEpsilon + 0.3 < REAL_MAX` always
passed → trial 1 always accepted → the revert branch was unreachable. The tell-tale was
in the trace itself: a line reading `residual 9.000000 ( 47.62 dB )` is internally
inconsistent (`10·log10(9)=9.5 dB`), proving `9.000000` is a display clamp
(`std::min(mEpsilon, ~9)` in the print path) while the dB column shows the true residual
(~5.8e4) — and that blown-up step was being *printed and accepted*, which a working
revert (print is after the loop) would have suppressed. Fix = baseline from `mEpsilon0`
+ flip to reject-if-worse + cap. Also fixed along the way: removed a `--mIteration` that
walked the counter backwards on backtracks.

## Status

Implemented; `mpicxx -fsyntax-only -Werror` clean (with/without the MUMPS define is N/A
here — controller is unconditional). Needs `make reset && make <target>` to land in the
static lib, then a runtime check on timestep 80 (Δt=50 ms): success = the 47.62/76.46 dB
excursions no longer print, and ω no longer climbs back to 1.0 to re-detonate. `0.3` is
the tuning knob (raise → ~0.5 if it reverts too eagerly in the −35 dB tail). Uncommitted.

## Related

- Floor/regularization detail and the crossover table → this session's analysis (above).
- Pairs with the [[dl20260622_mumps_error_messages]] work (same session; MUMPS warning
  +8 / convergence improvement traced separately to the staged `ICNTL(14)=30`).
