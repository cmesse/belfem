# Variable-Step BDF Time Integration (BDF1–BDF5) {#fem_kernel_bdf_timestepping_theory}

**Date:** 2026-07-31
**Purpose:** Theory and implementation reference for BELFEM's variable-step backward
differentiation formula (BDF) time integrators, including the startup order ramp, the
step-size history, and the savepoint mechanism used by the coupled controller.
**Module:** `src/fem/iwg` (`cl_IWG_Timestep`), consumed by `src/fem/kernel` (`cl_FEM_Controller`)

---

## 1. Baseline: the implicit Euler step

BELFEM's transient solvers advance a semi-discrete system

```
M dq/dt + K q = f
```

with an implicit (backward) scheme. The first-order baseline is the implicit Euler step
(Messe et al. 2023, Eq. 9):

```
( M + h K ) q_{n+1} = h f + M q_n ,      h = Δt
```

Both `M` and `K` may depend on the solution (nonlinear materials), so each step is iterated
by the nonlinear controller (see `nonlinear_controller_theory.md`). The scheme family is
selected per equation via `IWG_Timestep::set_timestepping_method()`
(`cl_IWG_Timestep.hpp:98`); the θ-schemes (forward, Crank–Nicolson, Galerkin, backward) and
the higher-order BDF members share one dispatch (`src/sparse/en_SolverEnums.hpp:46-61`).

## 2. Variable-step BDF-p

A BDF scheme of order p replaces the time derivative at `t_{n+1}` by the derivative of the
Lagrange polynomial interpolating the p+1 states `q_{n+1}, q_n, …, q_{n+1-p}`. Because the
adaptive controller changes the step size, BELFEM uses the **variable-step** coefficients:
with the upcoming step `h = Δt` and the previous steps `h_1, h_2, …` (stored newest-first in
`mH`, `cl_IWG_Timestep.cpp:305-308`), define the cumulative sums

```
S_2 = h + h_1 ,   S_3 = S_2 + h_2 ,   S_4 = S_3 + h_3 ,   S_5 = S_4 + h_4 .
```

The discrete system solved each step is

```
( α M + h K ) q_{n+1} = h f + M ( β_0 q_n − β_1 q_{n-1} + β_2 q_{n-2} − β_3 q_{n-3} + β_4 q_{n-4} )
```

truncated to the active order. The alternating sign pattern is centralized in
`IWG_Timestep::collect_qhist()` (`cl_IWG_Timestep.cpp:582-629`); `bdf2()` is the
corresponding assembly wrapper for the canonical shape (`cl_IWG_Timestep.cpp:766-784`).

The coefficients follow from differentiating the interpolating polynomial at `t_{n+1}`
(`compute_bdf_coeffs_2()` … `compute_bdf_coeffs_5()`, `cl_IWG_Timestep.cpp:1008-1140`).
For BDF2:

```
α   = ( 2h + h_1 ) / S_2
β_0 = S_2 / h_1
β_1 = h² / ( h_1 S_2 )
```

**Fixed-step sanity check.** For h = h_1 these reduce to α = 3/2, β_0 = 2, β_1 = 1/2 — the
classical BDF2 stencil (3/2 q_{n+1} − 2 q_n + 1/2 q_{n-1}) / h = q̇. The higher orders reduce
to the classical fixed-step coefficients the same way.

For BDF3 the same construction gives

```
α   = 1 + h/S_2 + h/S_3
β_0 = S_2 S_3 / ( h_1 ( h_1 + h_2 ) )
β_1 = h² S_3 / ( S_2 h_1 h_2 )
β_2 = h² S_2 / ( S_3 ( h_1 + h_2 ) h_2 )
```

and orders 4 and 5 extend the pattern with S_4, S_5 (`cl_IWG_Timestep.cpp:1045-1140`).

## 3. Startup ramp

A BDF-p step needs p history states, which do not exist at the beginning of a run. Instead of
seeding with a lower-accuracy bootstrap of the full history, BELFEM ramps the order: step k
runs BDF-min(k, p). A warm restart is different since 2026-08-15: the memdump carries the
field history AND the integrator state (`bdf_h`, `bdf_step_count`, `bdf_last_dt`), so a
resumed run continues at its earned order with the uninterrupted variable-step coefficients.
Since 2026-08-30 the integrator state is mandatory: a dump written by an older binary lacks the
keys and is refused by name rather than re-entering the ramp, so the ramp
(`compute_bdf_coefficients()`) is now reached only on a genuinely cold start. The active order `mOrderActive` also selects the matching
history truncation in `collect_qhist()`, so an unpopulated `mH` slot is never read (debug
builds assert this, `cl_IWG_Timestep.cpp:358-365`).

## 4. Step-size history, rejection, and savepoints

Three mechanisms keep the variable-step history consistent with the adaptive controller:

| Mechanism | What it does | Where |
|---|---|---|
| `shift_fields()` | pushes the accepted state into the history ring and the completed step onto `mH(0)` | `cl_IWG_Timestep.cpp:282-310` |
| `reset_fields()` | reverses one shift for a rejected timestep attempt: un-shifts the fields and the step-size history, including the deepest slot `mH(3)` (read by BDF5), which the shift would otherwise discard — `shift_fields` keeps a one-deep backup (`mHDropped`) for exactly this reversal, and the savepoint transports it | `cl_IWG_Timestep.cpp:461-495` |
| savepoint (`make_savepoint` / `restore_savepoint`) | deep copy of fields, `mH`, `Δt` and the step counter; restores across **any number** of shifts — required because the thermal equation may sub-step several times inside one magnetic step | `cl_IWG_Timestep.hpp:42-49` |

Coefficients are recomputed **lazily**: shifting or resetting marks them dirty
(`mCoeffsDirty`), and `compute_jacobian_and_rhs()` recomputes them on first use, after the
upcoming `Δt` is known (`cl_IWG_Timestep.hpp:120-124`). This ordering matters: the
coefficients depend on the *upcoming* step and the *already shifted* history.

## 5. Stability and accuracy notes

- In their classical fixed-step form, BDF1 and BDF2 are A-stable; BDF3–BDF5 are A(α)-stable with decreasing wedge angle
  (general theory: Bathe, Ch. 9 for implicit time integration in FEM; the classical BDF
  stability results are standard ODE literature). For the strongly damped
  magneto-thermal systems targeted here, the A(α) restriction is rarely binding; caution is
  only warranted for weakly damped oscillatory modes.
- Variable-step zero-stability requires bounded step-size ratios, and the admissible
  ratio bound shrinks with the BDF order (BDF2: 1+√2 ≈ 2.41; roughly 1.48 at order 3
  and near 1 at orders 4–5 for arbitrary sequences — Grigorieff 1983, quoted in
  Hairer & Wanner; smoothly varying sequences remain stable at all orders, which is
  how production BDF codes run order 5). The controller's growth clamp is therefore
  order-aware (`Controller::adjust_timestep`): ≤ 1.5 per step at active order ≤ 2,
  ≤ 1.4 at order 3, ≤ 1.2 at orders 4–5; shrinking (≥ 0.5) is always allowed.
  Savepoint alignment can still override an individual step outside the
  ratio clamp, bounded only by the min/max step limits.
- **Implicit schemes overestimate unstable modes.** For a locally unstable mode
  dy/dt = λy with λ > 0 (thermal runaway during a quench), implicit Euler produces the
  growth factor 1/(1 − λh), which exceeds the exact e^{λh} and diverges as λh → 1.
  Unconditional stability therefore does **not** imply accuracy through a quench transient:
  the step size must resolve the runaway timescale. This is a property of the method, not a
  defect; it is the reason the temporal resolution of quench simulations must be verified by
  a step-size sensitivity study.

## 6. Quick reference

| Item | Value / location |
|---|---|
| Scheme selection | `set_timestepping_method( EulerMethod, aHaveStiffness )`, `cl_IWG_Timestep.hpp:98` |
| Supported schemes | Static, θ-family (θ = 0, ½, ⅔, 1), BDF2…BDF5, derivative/mass/stiffness-only (`src/sparse/en_SolverEnums.hpp:46-61`) |
| Coefficients | `compute_bdf_coefficients()` + `compute_bdf_coeffs_p()`, `cl_IWG_Timestep.cpp:319-456, 1008-1140` |
| History combination | `collect_qhist()` — single source of truth for the sign pattern |
| Startup ramp | step k runs BDF-min(k, p) |
| Rejection handling | `reset_fields()` (one step) / savepoints (multi-step, coupled) |

## 7. Literature

- Messe et al. 2023, Section 4 — implicit Euler baseline (Eq. 9), residual and
  iteration framework the BDF steps are embedded in.
- Bathe, Ch. 9 — implicit direct time integration in finite element analysis.
- Variable-step BDF coefficient construction via Lagrange differentiation and the classical
  BDF stability theory (A-stability of orders 1–2, A(α)-stability of 3–5, step-ratio bounds)
  are standard numerical ODE results; see e.g. Hairer & Wanner, *Solving Ordinary
  Differential Equations II* (not part of the local reference library).
