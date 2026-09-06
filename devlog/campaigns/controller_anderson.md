# Campaign: Nonlinear Controller / Anderson Acceleration

**Status seed:** 2026-08-05 from dl20260722/0730 + `todo/anderson_picard_acceleration_plan.md`
— every claim `[seeded — confirm]` until Christian's correction pass.
**Branch:** `devel` — **COMMITTED** (corrected 2026-08-10): the bundle landed in `6b2a2b98`
(2026-07-30), the residual fix in `4f2c11cd` (08-07), the PID controller in `1d6ef305`
(08-09). The "UNCOMMITTED since 2026-07-30" line this page carried for five days was wrong.
**Plan:** `todo/anderson_picard_acceleration_plan.md`

## Current accepted design `[seeded — confirm]`

Hybrid Picard/QN/Newton controller per Messe et al. 2023 (paper1) §2.7. Repair bundle
(audited, amended, pending approval): R1 ω-carry fix (11 sites, not 4), R2 thermal
flip-count latch, R3 retry hygiene; then opt-in Anderson(m) mixing on the per-kernel
Picard branch (`gels(AbortOnError)` + `ShiftRegister`, fixed-point residual =
`mLhsVector − mFieldValues`, accept-commit staging, 5-class flush inventory). O2
resolved: β = live-ω with growth ACTIVE. ~~Anderson path must refresh `mFieldValues`.~~
**REVERSED (D5, fixed 2026-08-07):** refreshing before the residual multiply made ε the
*linear solver's* roundoff (‖A·A⁻¹b − b‖/‖b‖) rather than a nonlinear measure — false
convergence at `mMinNumIterations`. ε is now the pre-update force residual, identical to
depth 0 (DR-52). ~~`gRhoMin = 0` is correct (1e-16 clamp inside the power law below
0.91·Jc).~~ **Superseded 2026-08-10:** `mRhoMin` is now **0.0** by ruling — the 1e-16 floor
bound for all J below ≈0.874·jc on typical tape constants while `drho_powerlaw_dJ`
differentiates the *unfloored* law, i.e. an inconsistent Newton tangent across the whole
subcritical band. Scope
guards: Picard branch only, `anderson depth` key absent = bit-identical AIMD, no
outer-loop Anderson, no history persistence.

## Last passing reproducer `[seeded — confirm]`

T(3 ms) = 90.5 K solver-path-independent (ts17 campaign). No regression test yet —
synthetic R7 case is part of the plan; A/B tape benchmark gate R8 is Christian's.

## Open P0/P1 `[seeded — confirm]`

- ~~The whole bundle is UNCOMMITTED (DR-10): land or shelve before it rots.~~ **Closed — it
  landed.** What remains of DR-52 is the greg3 A/B *run*, not a commit. Anderson is opt-**in**
  (`mAndersonDepth = 0`, `cl_FEM_Controller.hpp:173-174`); the one-day opt-out experiment of
  2026-08-07 was rolled back the same day.
- Δt-sensitivity of the 20 ms/440 K result; ΔT-per-step limiter proposal (DR-11).
- Plan approval + O2 mixing-β sign-off still pending formally.

## Superseded approaches `[seeded — confirm]`

O2 "hold ω during growth" (rejected — the hold was a ratchet, D4); fixed-β Anderson;
outer-loop Anderson (out of scope by decision).

## Dated entries

dl20260722_coupled_newton_hardening · dl20260730_hphi_ts17_controller_strategy ·
todo/anderson_picard_acceleration_plan.md
