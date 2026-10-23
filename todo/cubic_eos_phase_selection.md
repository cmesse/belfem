# Cubic EoS Phase Selection

**Date:** 2026-10-22
**Purpose:** Let a cubic-EoS `Gas` choose the liquid or vapour root by itself and correctly, instead
of always returning the vapour root unless the caller sets a flag.
**Module:** `src/physics/gasmodels`
**AIs involved:** Claude (analysis); jury and literature pending
**Status:** OPEN — design. Not started. Raised by the 2026-10-22 usage-guide sweep (review_gasmodels
F12), discussed with Christian the same day.

> **Scope guards:**
> - The Helmholtz EoS path is unchanged: it decides the phase through its own vapour-pressure
>   ancillary (`Helmholtz::is_liquid`, `cl_GM_Helmholtz.cpp:720-735`).
> - Until this lands, a cubic `Gas` returns the vapour root unless the caller sets
>   `Gas::set_liquid_flag( true )`, and the flag belongs to the whole object, not to each point.

---

## 1. Current Behaviour

- `EoS_Cubic::v( T, p )` solves the cubic in Z by Cardano and returns the smallest root with
  Z > B when `mParent.is_liquid()` is true, otherwise the largest root (`cl_GM_EoS_Cubic.cpp:158-198`).
- `Gas::mLiquidFlag` defaults to `false` (`cl_Gas.hpp:183`). Its only writer is `Helmholtz::v()`
  (`cl_GM_Helmholtz.cpp:254`), so for a pure cubic gas it stays false.

**Bottom line:** the cubic can never return a liquid state on its own.

## 2. Why Gibbs Selection Alone Is Not Enough (discussion with Christian, 2026-10-22)

- **The input pair decides whether there is any ambiguity.** (T, v) gives p explicitly, so it is
  unique. (p, v) → T is unique for a decreasing α(T) (medium confidence, to be checked). Only
  (T, p) → v can have three real roots, between the two spinodals.
- **The vapour curve must be the cubic's own.** The stable root has the lower Gibbs energy; equal
  Gibbs energy means p equals the cubic's saturation pressure, found from equal fugacities. A switch
  on an external ancillary picks the wrong root near the curve.
- **Metastable states need history.** Between the binodal and a spinodal, subcooled vapour and
  superheated liquid are physical states, and the Gibbs rule cannot return them. A hard switch at
  p_sat also makes v jump, so a Newton iteration chatters across it. Hysteresis fixes both: keep the
  previous phase while its root still exists.
- **History must be per point.** `mLiquidFlag` is one flag per `Gas` object. If one `Gas` serves many
  nodes or quadrature points, the flag holds whichever point was evaluated last.

## 3. Steps

- [ ] **R1** — literature: the cubic's saturation curve from equal fugacities, and root selection
  (route through the `literature` skill).
- [ ] **R2** — compute and cache the cubic's own p_sat(T) (spline over [T_triple, T_crit]).
- [ ] **R3** — design per-point phase state owned by the caller, with hysteresis, and the Gibbs
  rule only for the first evaluation, which has no history.
- [ ] **R4** — plan jury, then code jury (`/code-with-jury`).

## 4. Open Questions

- **O1** — Who owns the per-point phase state: the FEM field (an extra nodal flag), or the caller of
  `Gas::v( T, p )`?
- **O2** — What do the current cubic-EoS users need? If they only use the cubic as the Helmholtz
  initial guess and for gas-phase states, R3 may not be needed.
