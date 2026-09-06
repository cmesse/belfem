# BELFEM ODE Module Tests — Detailed Plan

**Date:** 2026-03-24
**Purpose:** Method-level test matrix for the ODE module (`src/numerics/ode/`)
**Depends on:** `tests_0_strategy.md`, `tests_1_containers.md` (Cell, Vector, Matrix, ShiftRegister), `tests_2_linalg.md` (gesv)
**Confidence:** High on BDF coefficient verification. High on RK45/DOP853 against analytic solutions. Medium on adaptive step control edge cases.
**Status:** 47 tests implemented and passing.

---

## Module Overview

| Component | Content | Backend |
|---|---|---|
| `BDF` | Backward Differentiation Formula multi-step integrator (orders 1–6) | Vandermonde + LAPACK gesv |
| `ODE` | Abstract base class for ODE right-hand side | Pure virtual `compute()` |
| `Integrator` | High-level adaptive integrator wrapper | Function pointer dispatch |
| `RK45` | Runge-Kutta-Fehlberg 4(5) explicit adaptive stepper | Cash-Karp coefficients |
| `DOP853` | Dormand-Prince 8(5,3) explicit adaptive stepper | Hairer-Norsett-Wanner coefficients |
| `Status` | Enum: OK, TRAPPED, MAXIT | — |
| `Type` | Enum: RK45 (RK78 commented out) | — |

---

## Key Design Observations

### BDF Coefficient Computation

`BDF::compute_coefficients()` builds a Vandermonde matrix from normalized time steps and solves for the α-coefficients using `gesv` (LAPACK). For **uniform step sizes** (all h equal), the coefficients should match the well-known BDF formulas:

| Order | α₀ | α₁ | α₂ | α₃ | α₄ | α₅ | α₆ | β₀ (=h·coefficient for f) |
|---|---|---|---|---|---|---|---|---|
| BDF-1 | 1 | -1 | | | | | | 1 |
| BDF-2 | 3/2 | -2 | 1/2 | | | | | 1 |
| BDF-3 | 11/6 | -3 | 3/2 | -1/3 | | | | 1 |

These are verifiable against published tables. The `mCoefficients` vector stores these after `compute_coefficients()`.

### BDF Uses ShiftRegister

`ShiftRegister<real>` is a fixed-capacity FIFO buffer. `mH` stores the last N step sizes. `aY` stores the last N+1 solution values. The BDF order ramps up as more steps are pushed, starting at BDF-1 and increasing to the capacity limit.

### BDF `eval` vs `deval`

- `eval(aY, aF)` — Given the derivative `aF` at the current time, compute `y_{n+1}` using the BDF formula. This is forward integration.
- `deval(aY)` — Given the history `aY`, compute the derivative approximation `dy/dt ≈ Σ(αₖ·yₖ)/h`. This is backward differentiation.

Both have scalar and vector overloads.

### RK45 Uses Aliased Work Vectors

In `fn_ODE_RK45.cpp`, lines 1216-1220, all intermediate stage solutions `y1` through `y5` alias `aWork(6)`. This is intentional — each stage overwrites the same buffer because the previous stage's y-value is no longer needed. `ya` also aliases `aWork(6)`. The fourth-order solution overwrites all intermediate stages.

### DOP853 Is NOT Wired Into the Integrator

The `Integrator` constructor switch (line 633) only handles `Type::RK45`. `DOP853` exists as a standalone function with the same signature as `RK45`. To test DOP853, call `DOP853_init()` and `DOP853()` directly — do not go through the `Integrator`.

### The Integrator Returns Mutable References

`epsilon()`, `timestep()`, `time()`, `maxtime()`, `max_num_iterations()` all return `real&` or `uint&`. This is how BELFEM configures the integrator — by assigning to the returned reference.

### Status Enum

- `OK` — step converged
- `TRAPPED` — step converged but was clamped to `maxtime`
- `MAXIT` — maximum iterations exceeded without convergence

---

## Test Strategy: Use Known Analytic Solutions

The ODE module tests should define concrete `ODE` subclasses with known analytic solutions and verify the numerical result is within tolerance.

### Test Equations

| Name | ODE | y(0) | Solution | Domain |
|---|---|---|---|---|
| ExpDecay | y' = -y | 1.0 | e^(-t) | [0, 2] |
| SinCos | y' = cos(t) | 0.0 | sin(t) | [0, π] |
| Quadratic | y' = 2t | 0.0 | t² | [0, 1] |

For BDF, simpler verification: push known values from an analytic solution, call `deval`, compare with analytic derivative.

---

## Test File Structure

```
tests/ode/
├── test_BDF.cpp             # BDF coefficients, eval/deval (scalar + vector)
├── test_ODE_Integrator.cpp  # ODE base class, Integrator, RK45, DOP853
```

---

## 1. BDF Coefficient Verification

### 1.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `BDFCoefficientsOrder1` | Push 1 uniform step → coefficients match BDF-1: [1, -1] |
| `BDFCoefficientsOrder2` | Push 2 uniform steps → coefficients match BDF-2: [3/2, -2, 1/2] |
| `BDFCoefficientsOrder3` | Push 3 uniform steps → coefficients match BDF-3: [11/6, -3, 3/2, -1/3] |
| `BDFCoefficientsOrder4` | Push 4 uniform steps → coefficients match BDF-4: [25/12, -4, 3, -4/3, 1/4] |
| `BDFCoefficientsVariableStep` | Push non-uniform steps → verify sum of all coefficients == 0 (consistency property: constant function has zero derivative) |
| `BDFOrderRampsWithHistory` | (idea from ChatGPT) Push 1, 2, 3 steps → after each `compute_coefficients()`, verify `coefficients().length() == tH.size() + 1` |
| `BDFNoUpdateFlagReusesCoefficients` | (idea from ChatGPT) Compute once, mutate history, call `eval(..., false)` → old coefficients reused |

---

## 2. BDF Eval and Deval

### 2.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `BDFDevalConstantSignalReturnsZero` | (idea from ChatGPT) Fill history with constant value → deval returns 0 |
| `BDFDevalScalarLinear` | For y = t, deval should return dy/dt = 1 within tolerance |
| `BDFDevalScalarQuadratic` | For y = t², deval should approximate dy/dt = 2t |
| `BDFEvalScalarLinear` | For y = t (linear), BDF eval should recover exact value at next step |
| `BDFEvalScalarQuadratic` | For y = t², higher-order BDF with enough history should be near-exact |
| `BDFEvalVectorMatchesScalar` | Vector overload of eval produces same result as scalar on each component |
| `BDFDevalVectorMatchesScalar` | Vector overload of deval produces same result as scalar on each component |

### 2.2 Error Paths `[semantic + debug]`

| Test Name | Guard | What It Verifies |
|---|---|---|
| `BDFHighOrderThrows` | `BELFEM_ERROR` (always) | (idea from Gemini) `ShiftRegister` with capacity >= 7 → constructor throws |
| `BDFEvalCapacityMismatchThrows` | `BELFEM_ASSERT` (debug) | (idea from Gemini) `aY.capacity() != mH.capacity() + 1` → assert fires |
| `BDFDevalCapacityMismatchThrows` | `BELFEM_ASSERT` (debug) | Same for deval |

---

## 3. ODE Base Class

### 3.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `ODEDimensionAccessor` | Constructed with dimension N, `dimension()` returns N |
| `ODECheckEventsDefaultReturnsZero` | Default `check_events()` returns 0.0 |

---

## 4. RK45 Integration (via Integrator)

### 4.1 Tests `[semantic]`

All tests use concrete `ODE` subclasses with known solutions.

| Test Name | Equation | What It Verifies |
|---|---|---|
| `RK45ExponentialDecay` | y' = -y, y(0) = 1 | After integrating to t=1, y ≈ e^(-1) within `tTol = 1e-6` |
| `RK45SinCosine` | y' = cos(t), y(0) = 0 | After integrating to t=π/2, y ≈ 1 within `tTol` |
| `RK45TwoDimensional` | y₁' = y₂, y₂' = -y₁ (simple harmonic) | At t=2π, y₁ ≈ y₁(0), y₂ ≈ y₂(0) (orbit returns) |
| `RK45StatusOK` | Any test equation | `step()` returns `Status::OK` on successful convergence |
| `RK45StatusTrapped` | Integrate with small maxtime | `step()` returns `Status::TRAPPED` when clamped |
| `RK45TrappedPreservesStoredTimestep` | (idea from ChatGPT) When trapped, `aStep` is NOT overwritten (line 1321: only on `Status::OK`) |
| `RK45FixedTimestep` | With `set_auto_timestep(false)` | Integration still converges (no adaptive refinement) |
| `RK45MaxIterationsZeroReturnsMaxit` | (idea from ChatGPT) `max_num_iterations() = 0` → immediate `Status::MAXIT` |

---

## 5. DOP853 Integration (Direct Call)

### 5.1 Tests `[semantic]`

DOP853 is NOT wired into the `Integrator`. Tests call `DOP853_init()` and `DOP853()` directly.

| Test Name | Equation | What It Verifies |
|---|---|---|
| `DOP853ExponentialDecay` | y' = -y, y(0) = 1 | y(1) ≈ e^(-1) within tighter tolerance (1e-10) |
| `DOP853SinCosine` | y' = cos(t), y(0) = 0 | y(π/2) ≈ 1 within tight tolerance |
| `DOP853HigherOrderAccuracy` | Smooth problem | DOP853 achieves tighter tolerance than RK45 for same step count |
| `DOP853StatusTrapped` | Integrate with small maxtime | Returns `Status::TRAPPED` |

---

## 6. Integrator Configuration

### 6.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `IntegratorTypeAccessor` | Constructed with `Type::RK45`, `type()` returns `Type::RK45` |
| `IntegratorEpsilonMutableRef` | `epsilon()` returns mutable reference, assignment works |
| `IntegratorTimestepMutableRef` | `timestep()` returns mutable reference, assignment works |
| `IntegratorMaxtimeMutableRef` | `maxtime()` returns mutable reference, assignment works |
| `IntegratorMaxIterationsMutableRef` | `max_num_iterations()` returns mutable reference |
| `IntegratorTimeDoesNotDriveStep` | (idea from ChatGPT) Set `time()` to arbitrary value, pass separate `aT` to `step()` → step uses `aT`, not `time()`. Documents the actual contract. |

### 6.2 Error Paths `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `IntegratorUnknownTypeThrows` | (idea from Gemini) Casting invalid value to `Type` and passing to constructor → `BELFEM_ERROR` fires |

---

## 7. What We Do NOT Test (Deferred)

- Parallel ODE integration (does not exist)
- Event detection (commented out in RK45, `check_events` returns 0 by default)
- RK78 (commented out in Type enum)
- DOP853 via Integrator (not wired up — only RK45 is dispatched)
- Stiff problems (BDF is designed for stiff systems but testing requires stiff-specific test equations)
- BDF order > 6 (constructor asserts capacity < 7)

---

## 8. Implementation Notes for Claude Code

1. **Define concrete ODE subclasses inside the test file.** Each test equation is a small class inheriting from `belfem::ode::ODE` with a known analytic solution.

2. **For BDF tests, use ShiftRegister.** The BDF constructor takes a `ShiftRegister<real>&` for step sizes. Push uniform or variable steps, then call `compute_coefficients()` and verify.

3. **For RK45 tests, use the `Integrator` wrapper.** Create an ODE, construct `Integrator(ode, Type::RK45)`, set `timestep()`, `epsilon()`, `maxtime()`, then loop calling `step()` until convergence or trapped.

4. **For DOP853 tests, call directly.** Use `DOP853_init(ode, work)` then `DOP853(ode, t, y, step, work, epsilon, maxIter, tmax, autoStep)`.

5. **Tolerances for integration tests:** Use `tTol = 1e-6` for RK45 and `tTol = 1e-10` for DOP853 (8th order should be much more accurate).

6. **The BDF example code (bdf.cpp) has undefined variables** (`tZ`, `tG`, `tX`, `tdX` at lines 67-69). This is a bug in the example, not the library. Ignore it.

7. **BDF constructor asserts capacity < 7** via `BELFEM_ERROR` (always active). Do NOT test BDF-7 or higher.

8. **(from ChatGPT) `Integrator::time()` does NOT drive `step()`.** The `time()` property is stored but unused by the step function. The actual integration state comes from the `aT` argument passed to `step()`. Tests should not assume `time()` advances automatically.

9. **(from ChatGPT) Trapped step preserves stored timestep.** At line 1321, `aStep = h` only executes when `aStatus == Status::OK`. When TRAPPED, the original step is preserved. This is a subtle contract worth testing.

10. **(from ChatGPT) The example `ode.cpp` sets `maxtime() = 0.930` but `tTime = 1.0`.** This means `aTmax < t0`, which forces `h = aTmax - t0 < 0`. Either this is deliberate backward-integration support or a design hole. The test plan does NOT test this edge case — it needs investigation first.

11. **BDF eval/deval capacity checks use `BELFEM_ASSERT`, which is compiled out under `NDEBUG`.** Wrap those error path tests in `#if BELFEM_ASSERTIONS_ACTIVE` (the macro is exported by `assert.hpp`). Do not restate the condition as `#ifndef NDEBUG` — that is not equivalent when both `NDEBUG` and `DEBUG` are defined.

12. **BDF high-order check uses `BELFEM_ERROR`, which is compiled in for every build.** Leave it unguarded: `EXPECT_THROW` works in release too, because each test `main` calls `belfem::assert::set_throw_on_error( true )`. Note that "always active" describes the *check*, not the reaction — without that call a release build would `MPI_Abort` instead of throwing, killing the whole test binary.

---

## 9. Source Bugs Found and Fixed

Two source bugs were found and fixed during test implementation.
Full details in `devlog/dl20260324_test_suite_bugs.md` under "ODE Module".

| ID | Severity | Files | Description |
|---|---|---|---|
| BUG-ODE-1 | Medium | `cl_ODE_Integrator.cpp`, `en_ODE_Type.hpp` | DOP853 missing from Integrator dispatcher. Fixed by user: added `Type::DOP853` to enum and dispatch case to constructor. |
| BUG-ODE-2 | Low | `cl_ODE_Integrator.cpp` | `time()` was a dead property — `step()` never updated `mTime`. Fixed: `mTime = aT` after integration call. |
| BUG-ODE-3 | Low | `fn_ODE_RK45.cpp`, `fn_ODE_DOP853.cpp` | Persistent TRAPPED status after step rejection. Fixed by resetting status at the top of each iteration. |

**Design observation (not a bug):**
- BDF vector `eval` mutates `aY(0)` in-place — intentional for efficiency

---

## 10. Codex Audit Checklist

- [ ] Concrete ODE subclasses defined with known analytic solutions
- [ ] BDF coefficient tests verify against published BDF tables
- [ ] BDF eval/deval tested for scalar and vector overloads
- [ ] RK45 tested via Integrator with at least 2 test equations
- [ ] DOP853 tested via direct call with at least 1 test equation
- [ ] Status enum values (OK, TRAPPED, MAXIT) all tested
- [ ] Integrator configuration accessors tested
- [ ] `EXPECT_NEAR` for all floating-point comparisons
- [ ] BELFEM naming conventions (`t` prefix for locals)
- [ ] Deferred items documented in test file header
