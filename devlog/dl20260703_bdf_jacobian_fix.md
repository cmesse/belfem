# Devlog 2026-07-03 — BDF Jacobian Scaling Fix (B1-B4) + input hook + thermal savepoint

**Date:** 2026-07-03
**Topic:** Wire the BDF coefficient pipeline into the timestep cycle, fix the Newton Jacobian α scaling, add the Controller-owned input `method:` hook, and add the thermal savepoint that resolves the `initialize_thermal()` reset todo (`todo/closed/bdf_jacobian_scaling_bug.md`, items B1-B4 + B6 + C1/A1/A2; savepoint + hook as same-day follow-ups). **Todo CLOSED**; V1/B7 spun out to `todo/bdf_nonlinear_mass_verification.md`.
**AIs involved:** Claude (implementation), Codex (audit)
**Claude Confidence:** high post-audit (was medium ~80% pre-audit)
**Codex Audit Confidence:** high — ramp/reset/dispatch confirmed; 3 findings raised, 2 fixed same session (C1 gate, asserts), 1 deferred (B7)
**Literature References:** Hairer & Wanner (1996) II.4 (BDF Jacobian structure); Messe et al. 2023 (paper1) §2.7 (Picard-first per-timestep strategy, tolerance)

## Summary

BDF2-5 were non-functional: `compute_bdf_coefficients()` was never called, so `mAlpha`/`mBeta` stayed NaN and every BDF≥2 scheme produced an all-NaN system matrix. BDF1 worked only because it never reads those coefficients. This session wired the coefficient pipeline in (lazy dirty flag), added the startup order ramp, and added the α factor to the Newton correction `assemble_dJdx`. Two additional defects were found and fixed along the way (A1, A2 below). User approved source edits for this task.

## Key Findings

- Confirmed all 2026-06-22/2026-07-01 audit findings against the current tree before editing (coefficients never computed; `_nok` variants and `mBeta` in the blast radius; startup division-by-zero hazard for orders ≥ 3).
- **A1 (new):** the `aHaveStiffness = false` branch of `set_timestepping_method()` never set `mOrder` for BDF2-5 (`cl_IWG_Timestep.cpp:160-179` pre-edit), so `bdf2_nok`…`bdf5_nok` had history depth 1: `mFieldData` sized 2 and `qold(k>1)` falling back to auto-created zero fields (`cl_FEM_Calculator.cpp:1961-1968`).
- **A2 (new):** on a rejected timestep, `reset_fields()` popped `mH(0)` (size of last *completed* step) but the retry's `shift_fields()` re-pushed `mDeltaTime`, which still held the *abandoned* step size — the Controller assigns the halved `delta_time()` only after the shift (`cl_FEM_Controller.cpp:125/134`, halving at `:1062`). The last completed step size was thereby lost from the variable-step history.
- The lazy hook in `compute_jacobian_and_rhs()` covers all real assembly paths: `IWG_Maxwell` and `IWG_TransientHeatConduction` inherit it; `IWG_MaxwellPostproc` (projection-only) and `IWG_StaticHeatConduction` (static) override it but do not use BDF coefficients. Codex independently confirmed no bypass (`cl_FEM_DofManager.cpp:664`) and no problematic `delta_time()` consumer between reset and the next initialize.
- **C1 (Codex, high — fixed):** `EigenValues::compute_matrices()` temporarily switches to `MassOnly`/`StiffnessOnly` without resetting `mOrder` (`cl_FEM_DofMgr_EigenValues.cpp:118/126`); the ramp could then have dispatched `bdf1` instead of the utility scheme. Fixed by gating `compute_bdf_coefficients()` on `mMethod` ∈ BDF2-5 (all other methods: α = 1, `mTimestepActive = mTimestep`, early return). Also initialized the previously uninitialized `mMethod`.
- **B7 (Codex, high — deferred):** the only `dMdX_times_h` producer (`mt_maxwell_phi.cpp:90`, pre-existing `#HACK` marker) contracts with `phi0` alone; exact for BDF1, an unweighted approximation of the required β-weighted history contraction for BDF2-5 → inexact (typically still convergent) Newton. Needs the weighted history exposed to matrix producers; tracked as B7 in the todo.
- **Restart observation (Codex, medium — accepted):** `mStepCount`/`mH` are not in the memdump, so the first resumed BDF-p step runs BDF1. Required (empty `mH` would divide by zero) and consistent with the circuit-restart cold-start philosophy; future option is persisting `{mH, mStepCount}`.

## Changes Made

- `src/fem/iwg/cl_IWG_Timestep.hpp` —
  - `mAlpha` initialized to 1.0 instead of NaN (B3).
  - New members: `mCoeffsDirty` (lazy recompute flag), `mStepCount` (startup ramp counter), `mHaveStiffness`, `mTimestepActive` (per-step dispatch pointer; differs from the configured `mTimestep` only during the ramp).
- `src/fem/iwg/cl_IWG_Timestep.cpp` —
  - B1: `shift_fields()` / `reset_fields()` / `set_timestepping_method()` set `mCoeffsDirty`; `compute_jacobian_and_rhs()` recomputes coefficients once when dirty (first element assembly of a step, after `delta_time()` is set).
  - B2: `compute_bdf_coefficients()` computes the effective order `min(mOrder, mStepCount)` (≥1) and selects `mTimestepActive` (bdf1-4 or `_nok` counterparts during the ramp); `compute_timestep()` dispatches through `mTimestepActive`.
  - B3: `default:` case sets `mAlpha = 1.0` for order ≤ 1 schemes.
  - A1: `mOrder = 2…5` now also set in the no-stiffness branch of `set_timestepping_method()`.
  - A2: `reset_fields()` restores `mDeltaTime = mH(0)` before un-shifting the step-size history.
  - C1: `compute_bdf_coefficients()` early-returns (α = 1, configured dispatch) for every non-BDF2-5 method; `mMethod` initialized to BDF1.
  - Asserts (Codex checklist item 4): non-null `mTimestepActive` in `compute_timestep()`; positive `mDeltaTime` and positive `mH` slots before the coefficient formulas (debug-only).
- `src/fem/iwg/cl_TimestepMatrices.hpp/.cpp` — B4: `assemble_dJdx( adt, aAlpha = 1.0 )`; the `dMdX_times_x` term is scaled by α; call site passes `mAlpha`. Default parameter keeps other (currently nonexistent) callers safe.
- `src/fem/iwg/doc/iwg_usage_guide.md` — B6: new "Newton Jacobian for Nonlinear Time-Stepping" section plus coefficient-lifecycle and startup-ramp notes.
- `todo/bdf_jacobian_scaling_bug.md` — B1-B4/B6 ticked, A1/A2 recorded, stale claims struck through.

Both edited translation units pass `g++ -fsyntax-only` with the project's real build flags (`-Wall -Werror`). No build or regression run was performed (user runs builds).

## Follow-up (same day): input-file method hook

Christian wired the timestepping method into the input file (commit `a2535c0e`): the
Controller parses `solver { timestep { method: bdf1..bdf5 | explicit | crc } }` in
`set_params()` (`cl_FEM_Controller.cpp:1720`), `MaxwellFactory::create_controller()` and
`hphiTrun.cpp:81` (thermal) pull it via `Controller::euler_method()`; the old
`MaxwellFactory::mTimeStepping` member was removed. Claude review found and fixed two
linkage defects (uncommitted on top of `a2535c0e`):

- **`BELFEM_ERROR( true, … )` in the unknown-method branch never fires** (the macro
  triggers on a *false* check) → typos silently fell back to BDF1. Changed to `false`.
- **Method now set after `create_old_dof_fields()`.** The old flow configured the method
  before `create_field()`; the new flow configures it after (`create_controller`,
  `hphiTrun.cpp:81` — ThermalFactory still pre-sets BDF1 at `cl_ThermalFactory.cpp:152`).
  The old-timestep storage (`mFieldData`, sized `mOrder+1`) and the `phi1..phiN` mesh
  fields were therefore created for order 1, and `shift_fields()` would index past the
  end for BDF ≥ 2. Fix: `set_timestepping_method()` re-runs `create_old_dof_fields()`
  when `mField` is already linked (safe: `IWG::set_field` sets `mField`/`mMesh` together,
  and re-running is idempotent; the temporary EigenValues switch keeps `mOrder`
  unchanged, so its re-run is a no-op rebuild). Both drivers set the method before
  `load_memdump`, so restored field data still lands in existing fields.
- **`thermalTest` lost its method entirely.** Christian's removal of
  `ThermalFactory::mTimeStepping` also removed the only `set_timestepping_method` call
  on `thermalTest`'s path → null dispatch pointer at first assembly. Fixed by a
  constructor default in `IWG_Timestep` (`set_timestepping_method( BDF1 )`; only
  `IWG_Timestep` overrides that virtual, so in-ctor dispatch is safe) — the pointer is
  now never null by construction, and later reconfiguration recreates the storage.

## Follow-up (same day): thermal savepoint — resolves the `initialize_thermal()` todo

The old todo comment asked whether the thermal equation must shift at every sub-step and
how a magnetic reset can revert it. Resolution (implemented):

- **Shifting every sub-step is required.** The thermal BDF history must hold the previous
  sub-step at `delta_time2` spacing; shifting only once per magnetic step would make the
  time derivative wrong by up to the coupling factor. Shift frequency unchanged.
- **The defect was the magnetic reset path.** `reset_timestep()` called
  `mEquation2->reset_fields()` once, though the thermal equation had shifted N sub-steps
  within the failed magnetic step — and N un-shifts cannot recover the state anyway
  (each un-shift duplicates the deepest history slot).
- **Fix: savepoint.** New `IWG_Timestep::make_savepoint()` / `restore_savepoint()`: deep
  copies of all dof-field slots plus `mH`, `mDeltaTime`, `mStepCount`; restore marks the
  coefficients dirty. `initialize_timestep()` and `initialize_magnetic()` take the
  thermal savepoint at the top (before the circuit-failure early-return can trigger a
  reset), and `reset_timestep()` restores it instead of calling `reset_fields()`.
  Storage is allocated once and reused (re-sized only if the order changes); restore is
  idempotent, which also fixes the latent double-reset double-un-shift hazard for the
  thermal equation. `reset_thermal()` (single sub-step retry) keeps `reset_fields()` —
  its one-shift/one-un-shift pairing is correct.
- Also fixed in passing: the thermal block of `initialize_timestep()` assigned
  `delta_time()` *before* `shift_fields()` (opposite of the magnetic ordering), so under
  adaptive Δt the thermal `mH(0)` recorded the new step size instead of the completed
  one; swapped. `reset_timestep()` now also rolls `mTime02` back to `mTime0`.
- The magnetic equation keeps single-shift `reset_fields()` (correct pairing); migrating
  it to the savepoint as well would additionally preserve the deepest history slot
  across rejects — noted as a future option.

## Open Questions

- **B5 regression test** still open: BDF2-vs-BDF1 Newton-iteration comparison on the nonlinear heat problem (`thermalTest`), then the HTS benchmark. Needs a build + run.
- **B7** (from the audit): give the matrix producers access to the β-weighted history sum so `dMdX_times_h` satisfies the `assemble_dJdx` contract for BDF ≥ 2.
- Crank-Nicolson/Galerkin scale K by 0.5·Δt / Δt·2/3 in the system matrix, but `assemble_dJdx` scales `dKdX_times_x` by the full Δt for every scheme. Harmless for BDF (weight is exactly Δt) and currently for CN/Galerkin in practice (used with Picard), but the Newton tangent for those two schemes would be inconsistent — noted, not fixed.
- `reset_fields()` still cannot restore `mH(3)` (the deepest history slot is duplicated on un-shift) — pre-existing, only affects BDF5 across a rejected step. **RESOLVED 2026-07-31:** `shift_fields` now backs up the dropped slot (`mHDropped`) and `reset_fields` restores it; see `dl20260731_kernel_theory_docs.md`.

## Files Updated

- src/fem/iwg/cl_IWG_Timestep.hpp
- src/fem/iwg/cl_IWG_Timestep.cpp
- src/fem/iwg/cl_TimestepMatrices.hpp
- src/fem/iwg/cl_TimestepMatrices.cpp
- src/fem/kernel/cl_FEM_Controller.cpp (input-hook fixes + thermal savepoint wiring)
- src/fem/iwg/doc/iwg_usage_guide.md
- todo/closed/bdf_jacobian_scaling_bug.md (moved from todo/, marked CLOSED)
- todo/bdf_nonlinear_mass_verification.md (new — spun-out V1/B7 follow-up)
- todo/README.md (active → closed; added spun-out follow-up)
- devlog/dl20260703_bdf_jacobian_fix.md (this file)

Note: `src/fem/thermal/cl_ThermalFactory.{cpp,hpp}` also changed (Christian removed
`mTimeStepping`, method now flows from the Controller); the `IWG_Timestep` constructor
default covers `thermalTest`'s otherwise-unset path.
