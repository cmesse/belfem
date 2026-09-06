# Timestep input extension: scheme key + opt-out Anderson/BDF5 defaults

**Date:** 2026-08-07
**Purpose:** Christian's decision after the Newton-extension jury round
(`dl20260807_newton_extension_jury.md`): Anderson stabilization and BDF5
become opt-OUT rather than opt-in — "we can always change our mind later."
New deck syntax in the `timestep` section; uncommitted, not compiled.
**Module:** `src/fem/kernel` (`cl_FEM_Controller`), docs.

## New input contract (solver → timestep section)

```
timestep
{
    initial timestep : 0.1 ms ;
    maximum timestep : 0.1 s ;
    simulation time  : 6 s ;
    scheme           : bdf5 ;   // bdf1..bdf5, explicit, crank-nicolson, galerkin
    adapt timestep   : true ;
    anderson stabilization : true ;
    save every       : 100 ms ;
}
```

- **`scheme`** is the preferred key; **`method` stays as a legacy alias**
  (existing decks unaffected). `galerkin` added to the accepted values;
  crank-nicolson/galerkin still hard-error downstream with stiffness
  (`IWG_Timestep::set_timestepping_method`) — parseable, unusable, as
  intended ("we don't use these at the moment").
- **Default scheme is now BDF5** (`Controller::mTimeStepping`,
  cl_FEM_Controller.hpp) — decks without a scheme/method key (greg, greg2)
  silently move BDF1 → BDF5.
- **`anderson stabilization`** (bool, default **true**) is a master switch:
  absent or true → Anderson ON with default depths (magnetic 3, thermal 1 ≈
  Aitken, matching the doc examples); false → both depths zeroed.
  Contradiction guard: `false` combined with a nonzero explicit
  `anderson depth` in a nonlinear section is a hard error; an explicit
  `anderson depth : 0` still disables a single field.
- Depths are now forwarded to SolverData **unconditionally** at the end of
  `set_params` (previously only inside the `anderson depth` key branch, which
  is insufficient now that defaults are nonzero). `set_anderson_depth` is
  re-callable (deletes/reallocates the ShiftRegisters), and
  `set_thermal_kernel` repeats the thermal forward for late-linked kernels.

## Scope and caveats

- Opt-out applies to input-file-driven runs (`set_params`); programmatic
  drivers that never call `set_params` keep SolverData's own default
  (depth 0) and the IWG constructor's BDF1 placeholder — unchanged.
- The jury round's caveat stands and is now LOAD-BEARING: BDF5 under the
  ×1.5/×0.5 adaptive-Δt ratios is an unproven-stability default (F6,
  `dl20260807_newton_extension_jury.md`); Anderson-on changes the
  reported-residual semantics vs legacy depth-0 Picard (F3). The G-B1 /
  G-P1 A/B gates are the validation path for the new defaults; each default
  can be flipped back with one key (`scheme : bdf1` /
  `anderson stabilization : false`).
- Docs updated: `anderson_acceleration_theory.md` (opt-out contract, master
  switch), `nonlinear_controller_theory.md` input table (depth defaults,
  new `scheme` row).
