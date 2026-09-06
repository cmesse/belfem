# DR-33 — the reproducer was deleted, and the half that needs no deck is now tested

**Date:** 2026-08-11
**Purpose:** Work out what remains of DR-33 and close what can be closed without a run
**Modules:** `fem/iwg`, `tests/fem`

## The finding

DR-33 had been reduced to one item, V1: *"compare BDF2-vs-BDF1 Newton-iteration count on
the nonlinear heat problem (`thermalTest`)"*. Its reproducer column said `thermalTest`.

**`src/executables/thermalTest.cpp` was deleted on 2026-08-06** in `e4c02ac8` ("various
fixes"), alongside `periodictest.cpp`, as dead-code cleanup. `src/executables/` now holds
`hphirun.cpp`, `hphiTrun.cpp` and `electricalCircuit.cpp` and nothing else, and no file in
`src/` or `config/` still names it. The register and
`todo/bdf_nonlinear_mass_verification.md` both went on citing it for five days. The stale
CMakeFiles directories in the two build trees are what made it look present.

**Reviving it verbatim would revive dead code.** Recovered with
`git show e4c02ac8^:src/executables/thermalTest.cpp`, it is a scratch driver: hardcoded
`costheta.msh` path, hardcoded `SolverType::UMFPACK` — opt-in and OFF by default, the same
trap DR-49 was about — initial temperature 0 K, and its own
`// todo: implement cp and lambda in material definition for thermal problem`. So the
ρ(T) nonlinear mass V1 is *named* for was never actually configured in it. The deletion was
right; the register's pointer to it was the defect.

So DR-33 is blocked on a **decision** — which live deck replaces it — not on a run.

## The split

The defect V1 exists to guard is the one the BDF Jacobian work opened on:
`compute_bdf_coefficients()` was never called, `mAlpha`/`mBeta` stayed NaN, and BDF2-5
assembled an all-NaN system. Part of that machinery needs no solve at all.

`IWG_Timestep` is concrete, sits under no pure virtuals, and its constructor allocates
only its own coefficient buffers — no mesh, no kernel, no dof manager. New
`tests/fem/test_BdfTimestepMethod.cpp`, six cases in the `fem`/`fast` suite:

| test | what it pins |
|---|---|
| `DefaultMethodIsBdf1` | the constructor configures BDF1 itself — the "dispatch pointer is never null even if a driver forgets" guarantee |
| `MethodRoundTrips` | BDF1-5 selection survives `set_timestepping_method` |
| `RampHoldsFirstStepAtOrderOne` | every BDF2-5 run executes its first step at order 1 |
| `Bdf1IsNotRamped` | BDF1 is order 1 from the start and recomputation is idempotent |
| `NonBdfSchemesSkipCoefficientWork` | the method gate that protects the eigenvalue `MassOnly`/`StiffnessOnly` switch from inheriting stale BDF order state (C1 of the original plan) |
| `CrankNicolsonAndGalerkinAreRejected` | both are parseable and then hard-error before any dispatch pointer is built |

The ramp test is the load-bearing one and it is not a tautology: **before the first shift
`mH` is all zeros**, so a ramp that failed to hold step one at order 1 would walk straight
into `compute_bdf_coefficients()`'s own *"Invalid BDF history step size"* assert. The test
pins that it does not.

The CN/Galerkin case is worth its line too. That clause of DR-33 was **struck** on
2026-08-09 rather than fixed, on the reasoning that both methods hard-error before any
tangent is built — parseable, unusable, as intended. Nothing had ever tested the rejection,
so the struck clause rested entirely on a code reading. It does not any more.

## The seam, and what running it changed

With no deck available, a seam became the only route to covering
`compute_bdf_coeffs_2..5` at all. Christian approved it, and it is one line in
`cl_IWG_Timestep.hpp`:

```cpp
friend class BdfCoefficientProbe ;
```

Zero runtime cost, no widened public API, and the probe is defined only in the test and
used nowhere in `src/`. It poses a step-size history and reads back α and β. §4 of the
test then has three cases:

| test | what it pins |
|---|---|
| `ConstantStepReducesToTextbook` | α and β against the published BDF2-5 constants — 3/2, 11/6, 25/12, 137/60 and their β sets (Hairer & Wanner 1996 II.4). An independent reference, not a restatement of the implementation |
| `VariableStepIsExactOnPolynomials` | the property that actually *defines* BDF-p, and the one the constant-step check cannot reach: at non-uniform step size the formula must still differentiate every polynomial of degree ≤ p exactly |
| `RampUsesReducedOrderCoefficients` | a BDF5 run two steps in must produce BDF2's coefficients — `collect_qhist` truncates history to `mOrderActive`, so a mismatch would silently stop the Jacobian being the tangent of its residual |

**Running it corrected the test itself, which is the part worth keeping.** The exactness
case was first written with realistic 3e-4 step sizes, and there it was very nearly
vacuous. Exactness is scale-free, but the *failure* it must detect is not: a scheme of the
wrong order misses by O(h^p), which at h = 3e-4 is around 1e-13 — indistinguishable from
the roundoff of a correct answer. Measured on the first run, BDF5 reproduced the derivative
of t⁶ to 2e-13, i.e. **the test would have passed for a scheme it exists to reject.**

Two changes fixed it. The steps were rescaled to O(0.1), where the separation is fifteen
orders of magnitude — ~1e-15 for the degrees that must be exact, 16–73 % for the degree
that must not be. And the negation at degree p+1 is now *asserted* rather than merely
observed: "exact up to degree p" is a fingerprint of BDF-p only together with its own
failure at p+1, otherwise a scheme of too high an order passes just as happily.

A static check could not have found this. The test compiled, and would have gone green,
while testing almost nothing.

## Scope correction

An earlier note in the plan file — written before this work — claimed the coefficients
could be pinned through the **public** interface. They cannot, and that claim was corrected
in place. `compute_bdf_coefficients()` is indeed mesh-free, but its inputs and outputs are
not reachable: `mAlpha`, `mBeta`, `mH` and `mStepCount` are **private** with no accessor,
and the only writers of `mH`/`mStepCount` are `shift_fields`, `reset_fields` and
`restore_savepoint`, every one of which dereferences `mField` and calls
`DofManager::initialize()`. That is what made the seam above necessary rather than
optional, and why it was put to Christian as a decision instead of assumed.

## One guard that was decided by evidence

`CrankNicolsonAndGalerkinAreRejected` is wrapped in `#ifndef NDEBUG`, matching the
convention already used in `test_Spline.cpp`. That is not cosmetic here. `BELFEM_ERROR` is
always *active*, but its failure mode is not constant: it throws under debug and calls
`error_abort()` under release (`assert.hpp`). Since `set_timestepping_method` lives in
`libbelfem_iwg.a`, the behaviour is baked in at library build time. Compiling the probe
with `-DDEBUG` against the release-built library and running it **took the process down via
MPI_ABORT** rather than throwing — so without the guard, a release-mode `make check` would
have lost the test binary rather than reported a failure.

## Evidence

- `tests/fem/test_BdfTimestepMethod.cpp` — **nine cases** — compiles clean under **both**
  backends (Armadillo and Blaze) × NDEBUG and DEBUG, with the tree's
  `-Wall -Werror -pedantic-errors`.
- The header change is inert for production: `cl_IWG_Timestep.cpp` and
  `mt_maxwell_phi.cpp` recompile clean with it, and a `friend` declaration changes no
  layout, so the prebuilt libraries stayed link-compatible.
- **Executed green** as standalone scratchpad probes linked against the prebuilt
  `cmake-build-debug` libraries:
  - §1-§2: `order_active` = 1 for BDF2, 3, 4 and 5; method round-trip; the non-BDF gate.
  - §4 constant step: α = 1.5, 1.833…, 2.083…, 2.283… against textbook, max coefficient
    deviation 4.4e-16.
  - §4 ragged step: exact on t⁰..t^p to ≤ 1.2e-15 for every order, and wrong by
    16.3 %, 16.2 %, 31.9 %, 73.1 % respectively at degree p+1.
  - §4 ramp: BDF5 configured, α = 3/2 after two steps, 25/12 after four, 137/60 after five.
- The §3 debug case is skipped under NDEBUG by its own guard, as designed — the abort
  observation above is why that guard exists.
- Registered in `tests/fem/CMakeLists.txt`, so it runs under `make check-fast`.

## The run half is parked, not pending

Asked which live deck should replace `thermalTest`, Christian's answer was that **we do not
have one at this time.** That is a stronger statement than "the reproducer was deleted": a
nonlinear-mass ρ(T) thermal transient does not exist in the tree in any form, so V1's run
is not schedulable. It is therefore recorded as **parked**, not as an open action — an
action item nobody can act on is noise in a register whose whole purpose is the release
lens.

The unblock condition is written down so it is recognisable if it ever arrives: any deck
with a temperature-dependent ρcp thermal transient, run once at BDF1 and once at BDF2 at
fixed Δt, comparing Newton iteration counts to ±1. No new tooling is needed — `method :`
is already an input key.

## What that leaves standing, stated plainly

**Nothing in the tree exercises BDF2-5 end to end.** That is the honest residual, and it
should not be overstated in either direction:

- It is bounded. BDF2-5 is opt-**in** — the 2026-08-07 change that briefly made BDF5 the
  default was rolled back the same day — so no deck in use depends on the path. The
  original all-NaN defect is fixed, and its startup ramp is now unit-tested.
- It is real. `compute_bdf_coeffs_2..5`, the variable-step coefficient formulas, are
  verified by no test and no run.

**With no deck available, the seam described under "Scope correction" is now the only route
to that coverage** — which raises its value relative to when it was first noted as an
aside. Least-invasive form: one `friend` probe declaration in `cl_IWG_Timestep.hpp` (zero
runtime cost, no widened public API, no test name inside any production function), plus a
test comparing the computed α/β against the textbook constant-step BDF2-5 coefficients at
uniform Δt, which is what the variable-step formulas must reduce to. That remains
Christian's call and no code has been written for it.
