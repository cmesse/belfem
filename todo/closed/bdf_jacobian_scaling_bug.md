# BDF Jacobian Scaling Bug in TimestepMatrices

**Date:** 2026-01-21
**Purpose:** Document critical bug in Newton Jacobian assembly for BDF2-BDF5 time-stepping schemes
**Severity:** High (impacts convergence for nonlinear transient problems)
**Affected Modules:** `src/fem/iwg/`
**Refreshed:** 2026-06-22 — Codex-audited; the originally proposed fix is **unsafe** and **insufficient** (see below)
**Implemented:** 2026-07-03 — B1-B4 and B6 completed (Claude, user-approved). See devlog `dl20260703_bdf_jacobian_fix.md`.

## ✅ CLOSED 2026-07-03

BDF2-5 are functional. The core defect (`compute_bdf_coefficients()` never called →
`mAlpha`/`mBeta` permanently NaN → all-NaN system) is fixed, together with the α scaling
of the Newton correction and every hazard found in the two audit rounds. Landed and
Codex-audited: **B1** lazy coefficient wiring, **B2** startup order ramp (step 1 BDF1 →
step p BDF-p, reject-safe), **B3** α = 1 for order ≤ 1, **B4** α scaling in
`assemble_dJdx`, **B6** docs, **C1** `mMethod` gate protecting the temporary
`MassOnly`/`StiffnessOnly` switches, **A1** `mOrder` in the no-stiffness branch, and
**A2** the reset/retry step-size history. The follow-up input-file `method:` hook and the
thermal savepoint (resolving the `initialize_thermal()` todo) landed the same day — see
the devlog. Initial run confirmed working by Christian.

**Two forward-looking items were spun out to the active
[`bdf_nonlinear_mass_verification.md`](bdf_nonlinear_mass_verification.md):** V1 (was
B5) the BDF2-vs-BDF1 regression test, and B7 the `dMdX_times_h` producer contract (the
history-derivative term is contracted with `phi0` only — exact for BDF1, an unweighted
approximation for BDF2-5; affects Newton *rate*, not converged correctness). The checklist
and analysis below are preserved as the historical record.

---

## ⚠️ Status Update (2026-06-22, Codex safety audit)

**Do NOT apply the "3-line fix" below as written — it would break BDF1 (NaN) and would not
help BDF2-5 (which are currently non-functional).** Verified against the current tree
(Claude + Codex):

1. **`mAlpha` is permanently `NaN`.** It is initialized to `BELFEM_QUIET_NAN`
   (`cl_IWG_Timestep.cpp:39`) and assigned only inside `compute_bdf_coeffs_2/3/4/5()`
   (`:772,797,827,866`). Those are called only from `compute_bdf_coefficients()`
   (`:290`), which is **declared, defined, and never called anywhere in the repo**
   (whole-tree grep; Codex concurs). So the BDF coefficients are never computed.

2. **Consequence: BDF2-5 do not currently work.** `bdf2..bdf5` do `aJ *= mAlpha`
   (`:508,534,560,588`) → all-NaN system matrix. The doc's premise ("BDF1 correct,
   BDF2-5 just inaccurate") is wrong: BDF2-5 are non-functional, and **BDF1 is the only
   working scheme** — it works precisely because it never reads `mAlpha`
   (`bdf1()`, `:472-486`, builds `M + Δt·K` directly).

3. **Why the proposed fix is unsafe.** `assemble_dJdx()` is called **once,
   order-agnostically** (`:381`), and `mdJdx` feeds the global Newton Jacobian
   (`cl_FEM_DofManager.cpp:477/489/523/535` → `assemble_newton(…, dJdx())`). Passing
   `mAlpha` (=NaN) would poison BDF1's tangent → NaN. The non-BDF schemes
   (`explicit_euler`, `crank_nicolson`, `galerkin`) all scale M by **1.0**, not `mAlpha`.

The original analysis (below) is correct **about the α term in the Jacobian** but assumes a
working BDF coefficient pipeline that does not currently exist. Keep this task **active**.

## Implementation Checklist (added 2026-07-01, supersedes the 2026-06-22 numbered fix sequence)

Re-verified against the current tree on 2026-07-01 (Claude). Findings 1-3 above still held, with two additions:

- **Blast radius is larger than stated.** In addition to `bdf2..bdf5`, the `_nok` variants
  also scale by `mAlpha` (`bdf2_nok..bdf5_nok`, `cl_IWG_Timestep.cpp:630,656,682,709`),
  and **all** BDF≥2 schemes build their RHS history term from `mBeta`, which is
  NaN-initialized (`cl_IWG_Timestep.cpp:40`). The earlier claim "history term assembly:
  no issues found" (below, "Other Observations") was stale in practice: the β formulas
  are correct but were never evaluated before B1.
- **Startup hazard.** `mH` is zero-initialized (`cl_IWG_Timestep.cpp:41`), so during the
  first p−1 steps of a BDF-p run, `compute_bdf_coeffs_3/4/5` divide by zero. BDF2's
  `tH1` (`:766`) is whatever `mDeltaTime` held at the first `shift_fields()`. Wiring B1
  alone still yields NaN/inf at startup; B2 is required for orders ≥ 2.

Work items in dependency order:

- [x] **B1 — wire `compute_bdf_coefficients()` into the timestep cycle (the real blocker).**
  Implemented with a lazy dirty flag inside `IWG_Timestep`: `shift_fields()` /
  `reset_fields()` set `mCoeffsDirty = true`, and `compute_jacobian_and_rhs()`
  (`cl_IWG_Timestep.cpp:370`) recomputes once when the flag is set. This is
  self-contained: no base-class virtual, no changes to `Controller::initialize_*` call
  sites (`cl_FEM_Controller.cpp:125,224,316`), and the required ordering ("after
  `delta_time` is set, after `shift_fields`", `cl_IWG_Timestep.hpp:90-96`) holds
  because the first element assembly of a step happens after both. This
  also covers standalone drivers (`thermalTest.cpp:206`).
- [x] **B2 — startup order ramp.** Dispatch on `min(mOrder, steps taken)`: step 1 runs
  BDF1, step 2 BDF2, …, so no coefficient formula reads an empty `mH` slot. Requires a
  step counter in `IWG_Timestep`; `reset_fields()` rolls it back.
- [x] **B3 — `mAlpha = 1.0` for order ≤ 1.** Set it in the `default:` case of
  `compute_bdf_coefficients()` (`cl_IWG_Timestep.cpp:316-319`) and/or initialize
  `mAlpha = 1.0` instead of NaN (`cl_IWG_Timestep.hpp:35`, `cl_IWG_Timestep.cpp:39`), so
  every scheme that scales M by 1.0 gets the correct derivative factor in B4.
- [x] **B4 — α scaling in `assemble_dJdx`.** Do this only after B1-B3; with NaN α, it
  poisons BDF1's tangent. Add the `aAlpha = 1.0` parameter, scale `mdMdX_times_x`
  (Option 1 below), and pass `mAlpha` at the call site (`cl_IWG_Timestep.cpp:381`).
- [ ] **B5 — regression test.** Compare BDF2-vs-BDF1 Newton-iteration count on the
  nonlinear heat problem (see "Testing Strategy"); optionally add the HTS benchmark
  comparison.
- [x] **B6 — documentation.** Added the Newton-Jacobian section to
  `src/fem/iwg/doc/iwg_usage_guide.md` (snippet under "Next Steps") and refresh the stale
  claims in this file: the "3 lines / 1-2 hours" estimate and the "history term assembly
  correct" observation.

Codex's 2026-07-03 implementation audit confirmed the ramp/reset design and raised
three findings: two were fixed in the same session, and one remains a follow-up.

- [x] **C1 — ramp hijack of temporary `MassOnly`/`StiffnessOnly` switches (Codex, high).**
  `EigenValues::compute_matrices()` (`cl_FEM_DofMgr_EigenValues.cpp:118/126`) temporarily
  switches the method without resetting `mOrder`; with stale `mOrder > 1` and
  `mStepCount < mOrder`, the ramp could dispatch `bdf1` instead of the configured
  utility scheme, and stale α could leak into `assemble_dJdx`. Fixed by gating
  `compute_bdf_coefficients()` on `mMethod` being BDF2-5: every other scheme gets
  `mAlpha = 1.0` and `mTimestepActive = mTimestep` unconditionally.
- [ ] **B7 — `dMdX_times_h` producer contract for BDF ≥ 2 (Codex, high; follow-up).**
  `assemble_dJdx()` subtracts `dMdX_times_h` assuming it holds
  `(∂M/∂x)·(Σ β_i·x^{n-i})`, but the only producer (`mt_maxwell_phi.cpp:90`, already
  marked `#HACK: TO BE TESTED`) contracts with `phi0` alone. That is exact for BDF1 but an
  unweighted single-state approximation for BDF2-5. Consequence: inexact (but typically
  still convergent) Newton for nonlinear-mass BDF ≥ 2. Fix requires exposing the
  β-weighted history sum to the matrix producers (e.g. precompute it into
  `mCalc->qswap()` before `compute_mkf`); scope beyond this task.
- **Observation — restart ramps from BDF1 (Codex, medium; accepted behavior).**
  `mStepCount` and `mH` are not persisted in the memdump, so the first resumed step of a
  BDF-p run executes BDF1 even though old field slots were saved
  (`cl_Mesh.cpp:3005`). This matches the deliberate cold-start philosophy adopted for the
  circuit restart (`todo/closed/restart_circuit_state.md`: order ramps from 1); without
  persisted `mH`, ramping is *required* to avoid dividing by empty history
  slots. Future option: save `{mH, mStepCount}` alongside the fields.

Two additional defects were found and fixed alongside B1-B4 (2026-07-03):

- [x] **A1 — `mOrder` missing in the no-stiffness branch.** The `aHaveStiffness = false`
  branch of `set_timestepping_method()` never set `mOrder` for BDF2-5, so `bdf2_nok`
  through `bdf5_nok` ran with history depth 1 (`mFieldData` sized 2, wrong `mH` shift depth,
  `qold(k>1)` reading auto-created zero fields). Now mirrors the stiffness branch.
- [x] **A2 — `reset_fields()` corrupted the step-size history on retry.** The un-shift
  popped `mH(0)` (size of the last *completed* step), but the retry's `shift_fields()`
  re-pushed `mDeltaTime`, which still held the size of the *abandoned* step — the
  Controller assigns the halved `delta_time()` only after the shift
  (`cl_FEM_Controller.cpp:125/134`, halving at `:1062`). `reset_fields()` now restores
  `mDeltaTime = mH(0)` before un-shifting.

## Executive Summary

The `TimestepMatrices::assemble_dJdx()` method is missing the BDF coefficient α when assembling the mass matrix derivative contribution `(∂M/∂x)·x^n`. This causes incorrect Newton Jacobian for BDF2-BDF5, leading to slower or failed convergence in nonlinear transient problems.

**Impact:**
- ✅ BDF1 (Backward Euler): **Correct** (α = 1)
- ❌ BDF2-BDF5: **Incorrect** (missing α ≈ 1.5-3.6 scaling factor)

## Theoretical Background

### First-Order Nonlinear Transient System

```
M(x)·ẋ + K(x)·x = f(t)
```

### BDF Discretization (order p)

```
M(x^n)·(α·x^n - Σ β_i·x^{n-i}) + Δt·K(x^n)·x^n = Δt·f^n
```

Where:
- α = BDF coefficient for current time level
- β_i = BDF coefficients for history terms
- Variable time-step formulas in `IWG_Timestep::compute_bdf_coeffs_*`

### Residual

```
r(x^n) = M(x^n)·(α·x^n - Σ β_i·x^{n-i}) + Δt·K(x^n)·x^n - Δt·f(x^n)
```

### Exact Newton Jacobian

```
J = ∂r/∂x^n
  = α·M + α·(∂M/∂x)·x^n - (∂M/∂x)·(Σ β_i·x^{n-i})
    + Δt·K + Δt·(∂K/∂x)·x^n - Δt·(∂f/∂x)
```

**Key point:** The `(∂M/∂x)·x^n` term must be scaled by α.

## Current Implementation

### File Locations

- `src/fem/iwg/cl_TimestepMatrices.cpp:131-155` - `assemble_dJdx()`
- `src/fem/iwg/cl_IWG_Timestep.cpp:362` - Call site
- `src/fem/kernel/cl_FEM_DofMgr_SolverData.cpp:1572-1601` - Newton assembly

### Bug Location

**File:** `src/fem/iwg/cl_TimestepMatrices.cpp:135-137`

```cpp
if ( mFlags.test( static_cast< index_t >( MatrixFlag::dMdX_times_x )) )
{
    mdJdx += mdMdX_times_x ;  // ← BUG: Missing α scaling factor
}
```

**Should be:**
```cpp
mdJdx += aAlpha * mdMdX_times_x ;  // Need α factor for BDF consistency
```

### Current Assembly Logic

```cpp
// IWG_Timestep.cpp:362 (in compute_jacobian_and_rhs)
mTimeStepMatrices->assemble_dJdx( mDeltaTime ) ;  // Only passes Δt, not α
```

```cpp
// TimestepMatrices.cpp:131-154 (assemble_dJdx)
void TimestepMatrices::assemble_dJdx( const real adt )
{
    mdJdx.fill( 0.0 ) ;

    if ( mFlags.test( static_cast< index_t >( MatrixFlag::dMdX_times_x )) )
    {
        mdJdx += mdMdX_times_x ;  // Missing α!
    }

    if ( mFlags.test( static_cast< index_t >( MatrixFlag::dMdX_times_h )) )
    {
        mdJdx -= mdMdX_times_h ;  // Correct (no scaling)
    }

    if ( mFlags.test( static_cast< index_t >( MatrixFlag::dKdX_times_x )) )
    {
        mdJdx += mdKdX_times_x*adt ;  // Correct (Δt scaling)
    }

    if ( mFlags.test( static_cast< index_t >( MatrixFlag::dFdX )) )
    {
        mdJdx -= mdFdX*adt ;  // Correct (Δt scaling)
    }
}
```

## Proposed Fix

### Option 1: Pass α to assemble_dJdx (RECOMMENDED)

**Modify signature:**

```cpp
// cl_TimestepMatrices.hpp:194-195
void
assemble_dJdx( const real adt, const real aAlpha = 1.0 );  // Add α parameter with default
```

**Update implementation:**

```cpp
// cl_TimestepMatrices.cpp:131-154
void
TimestepMatrices::assemble_dJdx( const real adt, const real aAlpha )
{
    mdJdx.fill( 0.0 ) ;

    if ( mFlags.test( static_cast< index_t >( MatrixFlag::dMdX_times_x )) )
    {
        mdJdx += aAlpha * mdMdX_times_x ;  // FIX: Scale by BDF coefficient
    }

    if ( mFlags.test( static_cast< index_t >( MatrixFlag::dMdX_times_h )) )
    {
        mdJdx -= mdMdX_times_h ;
    }

    if ( mFlags.test( static_cast< index_t >( MatrixFlag::dKdX_times_x )) )
    {
        mdJdx += mdKdX_times_x*adt ;
    }

    if ( mFlags.test( static_cast< index_t >( MatrixFlag::dFdX )) )
    {
        mdJdx -= mdFdX*adt ;
    }
}
```

**Update call site:**

```cpp
// cl_IWG_Timestep.cpp:362
mTimeStepMatrices->assemble_dJdx( mDeltaTime, mAlpha ) ;
```

**Changes required:**
1. `src/fem/iwg/cl_TimestepMatrices.hpp:195` - Add parameter
2. `src/fem/iwg/cl_TimestepMatrices.cpp:131` - Update signature and line 137
3. `src/fem/iwg/cl_IWG_Timestep.cpp:362` - Pass mAlpha

**Pros:**
- Minimal code changes
- Default parameter maintains backward compatibility
- Clear mathematical meaning

**Cons:**
- None

### Option 2: Compute dMdX_times_x with α·x^n in IWG

Have each physics IWG contract `∂M/∂x` with `α·x^n` instead of `x^n` when computing the derivative.

**Pros:**
- No change to TimestepMatrices interface

**Cons:**
- Requires changes in all physics IWG implementations
- Less clear separation of concerns
- α is a time-stepping parameter, not a physics parameter

### Recommendation

**Use Option 1.** It's cleaner, requires minimal changes, and keeps the time-stepping logic where it belongs (in the time-stepping classes, not in physics IWGs).

## Testing Strategy

### 1. Unit Test: Simple Nonlinear Heat Equation

**Problem:**
```
ρ(T)·c·∂T/∂t = ∇·(k·∇T) + Q
```

Where ρ(T) = ρ₀·(1 + β·T) (nonlinear mass matrix).

**Test cases:**
- 1D rod with known analytical solution
- Compare BDF1, BDF2, BDF3 convergence rates
- Fixed time-step to isolate Jacobian accuracy effect

**Expected results:**
- **Before fix:** BDF2-BDF5 require more Newton iterations than BDF1
- **After fix:** All BDF orders have similar Newton iteration counts

### 2. Regression Test: Magnetodynamics

**Problem:**
```
σ·∂A/∂t + ∇×(ν(B)·∇×A) = J_source
```

Where ν(B) = ν₀/(1 + |B|/B_sat) (nonlinear reluctivity → nonlinear mass).

**Test cases:**
- Existing HTS benchmarks from `literature/papers/fem/messe2023.txt`
- Monitor Newton iteration counts for BDF2 vs BDF1
- Check convergence history with tight tolerance (ε < 10⁻¹¹ per Messe et al. 2023)

**Expected results:**
- **Before fix:** BDF2 may require Picard pre-iterations or fail to converge
- **After fix:** BDF2 should converge faster, matching theoretical quadratic convergence

### 3. Verification: Compare with Reference Implementation

Compare against GetDP, COMSOL, or FreeFEM++ for same nonlinear transient problem.

## Literature References

### BDF Theory

**Standard references:**
- Gear (1971) "Numerical Initial Value Problems in Ordinary Differential Equations"
- Hairer & Wanner (1996) "Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems"
  - Chapter II.4: Backward Differentiation Formulas
  - Establishes that BDF Jacobian must include ∂M/∂x contributions scaled by α

**From BELFEM literature:**
- **Bathe (2016), Chapter 9** - Implicit time integration for dynamics (focuses on Newmark/Houbolt for second-order systems, but principles apply)
- **Hughes (2000), Chapter 8** - Time-stepping for parabolic problems (first-order systems)
- **Messe et al. (2023) paper1, Section 2.7** - Emphasizes tight convergence tolerance (ε < 10⁻¹¹) for HTS simulations, indicating sensitivity to Jacobian accuracy

### Newton-Raphson for Transient Problems

General FEM principle (Hughes, Bathe, Zienkiewicz):
> For nonlinear M(u) and K(u), the consistent tangent requires **full derivative** of residual w.r.t. unknowns.

Omitting Jacobian terms leads to "modified Newton":
- Still convergent (if radius of convergence is reached)
- **Slower** convergence (no longer quadratic)
- **More iterations** required
- May **fail** for stiff problems (common in HTS with field-dependent μ_r)

## Related Code Review

### Positive Aspects

The overall architecture is **well-designed**:

✅ **Clean separation of concerns:**
- `TimestepMatrices` - Storage and element-level assembly
- `IWG_Timestep` - Time-stepping logic (BDF coefficients, method selection)
- `SolverData` - Global assembly and Newton solve

✅ **Efficient contraction:**
- Storing `(∂X/∂x)·v` as n×n instead of n×n×n tensor is correct

✅ **Flag-based assembly:**
- Flexible for different physics (with/without K, D, derivatives)

✅ **Separate J and A matrices:**
- `mSystemMatrix` (A) for Picard
- `mJacobianMatrix` (J) for Newton
- Allows hybrid strategies (Picard → Quasi-Newton → Newton)

✅ **Variable time-step BDF:**
- Correctly implemented in `compute_bdf_coeffs_*` functions (cl_IWG_Timestep.cpp:742-801)

### Other Observations

**No issues found in:**
- BDF coefficient calculation (correct variable time-step formulas)
- ~~History term assembly (β coefficients applied correctly)~~ *(stale — the β formulas
  are correct, but they were not evaluated before B1; see the 2026-07-01 checklist note)*
- RHS assembly (f terms scaled by Δt correctly)
- Element-to-global assembly (DofMgr_SolverData.cpp)

~~**The bug is isolated to this single location:** Missing α in line 137 of
cl_TimestepMatrices.cpp.~~ *(stale — the blast radius includes the never-called coefficient
pipeline, the `_nok` variants, missing `mOrder` in the no-stiffness branch (A1), and
the reset/retry step-size history (A2); see the Implementation Checklist)*

## Implementation Priority

**Priority:** High

**Rationale:**
1. Affects all nonlinear transient simulations with BDF2-BDF5
2. HTS magnetodynamics (BELFEM's primary application) has nonlinear μ_r(H) → nonlinear mass
3. Paper1 (Messe et al. 2023) emphasizes tight convergence tolerance needed for HTS
4. Fix is simple and low-risk (3 lines changed)

~~**Estimated effort:** 1-2 hours (implementation + testing)~~ *(stale — the actual scope
was B1-B4 plus A1/A2; implementation landed 2026-07-03, and testing (B5) remains open)*

**Risk:** Low
- Default parameter ensures backward compatibility
- BDF1 behavior unchanged (α = 1)
- Existing tests should still pass

## Success Criteria

1. **Code review:** Jacobian assembly matches theoretical formula
2. **Unit test:** BDF2 Newton iterations ≈ BDF1 iterations (±1) for nonlinear heat equation
3. **Regression test:** Existing HTS benchmarks converge faster or with fewer Picard pre-iterations
4. **No regressions:** All existing tests pass

## Next Steps

1. **Implement fix** (Option 1)
2. **Create unit test** for nonlinear heat equation
3. **Run regression tests** on HTS benchmarks
4. **Update documentation** in `src/fem/iwg/doc/iwg_usage_guide.md`:
   ```markdown
   ## Newton Jacobian for Nonlinear Time-Stepping

   For BDF schemes with nonlinear M(x), the Newton Jacobian includes:
   - α·(∂M/∂x)·x^n: Current iterate derivative (scaled by BDF coefficient α)
   - (∂M/∂x)·(Σ β_i·x^{n-i}): History term derivatives (not scaled)
   - Δt·(∂K/∂x)·x^n: Stiffness derivative (scaled by time-step)
   - Δt·(∂f/∂x): Load derivative (scaled by time-step)

   Reference: Hairer & Wanner (1996) II.4 for BDF Jacobian structure.
   ```

## References

- **Code locations:**
  - `src/fem/iwg/cl_TimestepMatrices.hpp:44-199`
  - `src/fem/iwg/cl_TimestepMatrices.cpp:1-160`
  - `src/fem/iwg/cl_IWG_Timestep.hpp:26-293`
  - `src/fem/iwg/cl_IWG_Timestep.cpp:27-805`
  - `src/fem/kernel/cl_FEM_DofMgr_SolverData.cpp:1572-1601`

- **Literature:**
  - Messe et al. (2023) - paper1 (BELFEM core paper)
  - Hairer & Wanner (1996) - BDF theory
  - Bathe (2016) §8, §9 - Implicit time integration
  - Hughes (2000) Ch. 8 - Parabolic time-stepping

---

**Status:** CLOSED 2026-07-03 (B1-B4, B6, A1, A2, C1; Codex-audited). Verification (V1) and
the B7 history-derivative contract spun out to `bdf_nonlinear_mass_verification.md`.
**Assigned:** Claude (implementation, complete)
**Target completion:** Done
