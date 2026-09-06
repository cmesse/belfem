# Devlog 2026-07-14 — BDF History Tangent Fix (B7 / phi_ferro #HACK) + Ferro-Conductor Semantics (O11)

**Date:** 2026-07-14
**Topic:** Consistent-tangent history contraction for nonlinear-mass producers
(`phi_ferro`), closing B7 of `bdf_nonlinear_mass_verification.md`; plus the decided
dispatch semantics for ferro-conductors (maxwell_kernel_collapse O11).
**AIs involved:** Claude/Fable (analysis + implementation), Codex (audit launched),
Grok (unavailable — truncation issue, see feedback memory)
**Claude Confidence:** high (design + edits); pending Codex audit + Christian's build
**Literature References:** consistent/algorithmic tangent — Belytschko §5.9,
Zienkiewicz & Taylor Vol 2 §4.4 (Simo & Taylor 1985)

## Summary

`phi_ferro` contracted its `dMdx_times_h` block with the single previous state
(`|B·phi0|`) — exact for BDF1, an underweighted approximation (factor ≈ Σβᵢ = α:
3/2 BDF2 … 137/60 BDF5, plus O(Δt) direction error) for BDF ≥ 2. Converged results were
never affected (the residual uses exact β weights); only Newton's rate. Fixed by a
single-source-of-truth history helper.

## Changes Made

- `IWG_Timestep::collect_qhist()` (new, `cl_IWG_Timestep.{hpp,cpp}`): builds
  `β₀q₀ − β₁q₁ + β₂q₂ − β₃q₃ + β₄q₄` into `mCalc->qswap()`, truncated to the new
  `mOrderActive` member (set on both `compute_bdf_coefficients` paths → startup ramp and
  variable Δt respected; order ≤ 1 → `q₀` exactly).
- `bdf2`–`bdf5` RHS bodies deduplicated onto the helper (same expressions relocated —
  FP-identical; the per-method sign patterns now exist in exactly one place).
- `phi_ferro` (`mt_maxwell_phi.cpp`): `dMdx_times_h` contracts against `|B·qhist|`;
  `#HACK: TO BE TESTED` removed; the by-value `phi0` copy (hidden per-element allocation)
  removed with it. `dMdx_times_x`/`M` unchanged.
- Design note: Christian directed "use the ShiftRegister"; code check showed the vector
  history lives as mesh fields (`label+n`, shifted by `shift_fields()`, gathered by
  `qold()`), with `ShiftRegister` (even post-upgrade) not the practical home —
  savepoint/restart/MPI machinery already rides the fields. Christian concurred
  ("not must use ShiftRegister if fields are more practical").
- Deliberately deferred to its own commit: replacing the isotropized tangent with the
  exact rank-1 form (`u uᵀ·μ′/H`, `u = Bᵀ(Bφ)`) — it would break the BDF1 bit-identity
  verification oracle of THIS change.

## Verification Gates (Christian)

1. BDF1 run: bit-identical (qhist ≡ q0 — anything else is a bug).
2. BDF2 on the dipole/iron validation case: converged fields unchanged, per-timestep
   Newton iterations same-or-better.

## O11 (decided, logged in maxwell_kernel_collapse_plan.md §5)

"Ferro" = field-dependent μ WITHOUT current (topological role; cuts may penetrate).
Iron-on-Conductor is legal — the dispatcher must route on the material's μ flags at link
time (`is_constant(mu)`), to variable-μ h-kernels mirroring the fixed `phi_ferro`
(exact rank-1 tangent available at HTS-Newton-term cost). Corollary: today
iron-on-Conductor silently runs with μ0 mass. Remaining: validation problem (Christian).

## Files Updated

- src/fem/iwg/cl_IWG_Timestep.hpp / .cpp (collect_qhist, mOrderActive, bdf2-5 dedup)
- src/fem/maxwell/matrices/mt_maxwell_phi.cpp (phi_ferro)
- todo/bdf_nonlinear_mass_verification.md (B7 ticked)
- todo/maxwell_kernel_collapse_plan.md (O11)
