# Nonlinear Iteration Strategy for HTS Problems Near Quench

**Date:** 2026-02-03
**Purpose:** Theory analysis and recommendations for improving nonlinear solver performance near critical material transitions (quench)
**Module:** `src/fem/kernel/cl_FEM_Controller.*pp`
**Problem:** Iteration counts reach ~100 near J → J_c instead of target 5-20
**Status:** THEORY REFERENCE — **partly implemented, partly superseded, one part refuted by
events** (reassessed 2026-08-09). This is the oldest file in `./todo/` (2026-02-03) and
predates the entire controller campaign; keep it for its literature argument, not as a work
plan.

> **2026-08-09 currentness sweep — recommendation by recommendation:**
>
> 1. **"Adaptive time-stepping" (PRIORITY 1) — the Δt-control half is implemented and has
>    moved on.** `Controller::adjust_timestep()` (`cl_FEM_Controller.cpp:1870`) adapts Δt from
>    the iteration count with the `sqrt( mIterationTarget / mIteration )` rule, gated by
>    `adapt timestep` in the `timestep` section. Replacing exactly that rule with a log-space
>    PID is now its own live plan, **`pid_timestep_controller_plan.md`** — that file, not this
>    one, is where Δt-control work belongs.
> 2. **The quench-detection trigger is NOT implemented.** Nothing computes a J/J_c margin in
>    the Controller and shrinks Δt ahead of the transition; the controller only reacts to
>    iteration counts after the fact. This is the one genuinely unclaimed idea in this
>    document, and Q1–Q3 below are still the right open questions for it. Related but
>    distinct: the ΔT-per-step limiter proposal (`debt_register.md` DR-11).
> 3. **"Lower iteration target near singularities" — not implemented;** `mIterationTarget`
>    is a single deck-level constant (`:2488`).
> 4. **"NOT line search" — this recommendation was overtaken by evidence.** A backtracking
>    line search was added to `iterate_coupled` in 2026-06 and is load-bearing today
>    (`:805-994`), and Anderson(m) mixing shipped in `6b2a2b98` with a five-case test suite (it
>    was briefly opt-out, and is opt-in again after the 2026-08-07 rollback — but implemented
>    and wired either way). The claim below that line search / Anderson "won't help near
>    singularity" should be read as a 2026-02 hypothesis that practice did not sustain — do
>    not cite it as settled. The nuance history actually delivered is finer than either
>    position: the Picard line search was *retired* in 2026-08-07's addendum 2 as degenerate
>    for lagged-operator Picard, while the Newton line search stayed. So "not line search"
>    was half right, for a reason this document did not anticipate.
> 5. **"NOT arc-length"** — still the right call; nothing has challenged it.
>
> Net: the 5–10× speedup claim was never measured, and the two mechanisms that would deliver
> it (PID Δt control, quench-margin trigger) are now tracked as (1) a separate plan and (2)
> the single open item of this file.

Original status line: THEORY COMPLETE - Ready for implementation

---

## Executive Summary

**Problem diagnosis:** High iteration counts (~100) near quench are a **symptom, not a disease**. The root cause is temporal discretization error when Δt is too large for the rapidly changing material response near J → J_c.

**Literature consensus (Bathe §8.4.1, Zienkiewicz Vol 2 Ch 3):**
> "The primary procedure for reaching convergence (if convergence difficulties are encountered) is to **decrease the magnitude of the load step**." (Bathe, p. 758)

**Recommended approach:**
1. **Adaptive time-stepping with quench detection** (PRIORITY 1) - Reduce Δt automatically when approaching J_c
2. **Lower iteration target** near singularities (5-10, not 20)
3. **NOT line search** - Wrong tool for material nonlinearity
4. **NOT arc-length** - Applies to structural instabilities, not material singularities

**Expected improvement:** 5-10× reduction in iteration count (100 → 10-20) by matching Δt to physics time scale.

---

## Theory: Why 100 Iterations Near Quench?

### Material Stiffness Singularity

HTS constitutive law (power-law E-J relation):

```
σ(J) = E_c (J/J_c)^n    where n ≈ 35
```

Material tangent modulus:

```
dσ/dJ = (n·E_c/J_c) · (J/J_c)^(n-1)
```

**Critical observation:** As J → J_c, the tangent modulus **explodes exponentially**:

| Margin to J_c | J/J_c | (J/J_c)^34 | Relative stiffness |
|---------------|-------|------------|---------------------|
| 50% margin | 0.50 | 1.78 × 10⁻¹¹ | 1.0× (baseline) |
| 10% margin | 0.90 | 0.0262 | 1.47 × 10⁹ |
| 1% margin | 0.99 | 0.714 | 4.01 × 10¹⁰ |
| 0.1% margin | 0.999 | 0.893 | 5.01 × 10¹⁰ |

**Interpretation:** The effective time scale for material response **shrinks by 9-10 orders of magnitude** as J → J_c. Using the same Δt that worked at 50% margin requires proportionally more iterations to resolve the rapid changes.

### Why This Isn't a Nonlinear Solver Problem

From Bathe §8.4.1 (p. 775):

**Newton-Raphson quadratic convergence properties:**

1. **If** current iterate is "sufficiently close" to solution U*
2. **And** tangent stiffness matrix K satisfies Lipschitz continuity
3. **Then** convergence is quadratic: error^(k+1) = O(error^k)²

**Near quench, condition (2) fails:**
- Lipschitz constant L in ||K(u₁) - K(u₂)|| ≤ L||u₁ - u₂|| becomes **enormous**
- Material stiffness changes by 10¹⁰ over small displacement changes
- Quadratic convergence radius shrinks to ~10⁻¹⁰ × Δt

**Practical consequence:** Newton iteration still works, but the "basin of attraction" becomes microscopic. Each iteration makes progress, but convergence is linear (not quadratic) because we're outside the quadratic convergence region.

**The fundamental issue:** Δt is too large for the physics happening at this time scale. No amount of clever nonlinear iteration will fix discretization error.

---

## Why Line Search Fails (And Why Arc-Length Doesn't Apply)

### Line Search Diagnosis

**Line search theory (Nocedal & Wright Ch 3):** Find ω ∈ (0,1] minimizing:

```
φ(ω) = ||F(x^k + ω Δx)||
```

**Why it returns ω = 1.0 in your implementation:**

**Hypothesis 1: State-dependent operators not updated**

For line search to work, must evaluate F(x + ω Δx) with:
1. Updated DOF values: `x_trial = x + ω Δx`
2. **Recomputed derived fields:** `J = curl(H)` at new state
3. **Updated material properties:** `σ(J_new)` with power-law
4. **Full element loop:** Reassemble residual with new σ

If residual uses stale material properties, φ(ω) doesn't vary → defaults to ω=1.0.

**Hypothesis 2: Non-descent direction (Picard iterations)**

From Bathe §8.4.2 (BFGS method, p. 759):
> "Line search requires a descent direction where ΔU^T · (R - F) > 0"

**Picard iterations:** `x^(k+1) = G(x^k)` where direction `G(x) - x` is **not guaranteed descent**.

**Practical consequence:** Line search applied to Picard will fail or behave erratically.

### Arc-Length Methods (Bathe §8.4.3)

**What arc-length solves (p. 761-763):**
- Structural snap-through (buckling, post-collapse)
- Load-displacement limit points (dλ/du = 0)
- Enables tracking beyond structural instability

**Governing equations:**
```
λR - F(u) = 0           (n equations, n+1 unknowns)
f(Δλ, Δu) = 0           (constraint: spherical arc-length)
```

Where λ is **load multiplier**.

**Why this doesn't apply to quench:**

| Feature | Arc-Length Target | HTS Quench Problem |
|---------|------------------|-------------------|
| Nonlinearity type | Geometric (buckling) | Material (power-law) |
| Critical point | Load maximum (∂λ/∂u = 0) | Material singularity (σ→∞) |
| Load control | λR (proportional loading) | Time-dependent (magnetodynamics) |
| Solution | Continue in λ-u space | No load parameter to vary |
| Jacobian | Becomes singular at limit point | Well-conditioned but stiff |

**Conclusion:** Arc-length is the wrong tool. The problem isn't traversing a limit point - it's resolving rapid material changes.

---

## Recommended Strategy: Adaptive Time-Stepping with Quench Detection

### Theoretical Foundation

**From Zienkiewicz & Taylor Vol 2 (Ch 3, Nonlinear Problems):**
> "For problems with severe material nonlinearity, automatic time step control based on iteration count provides robust convergence without trial-and-error."

**From Bathe §8.4.1 (p. 758):**
> "The primary procedure for reaching convergence is to decrease the magnitude of the load step."

### Algorithm

**Step 1: Compute quench proximity indicator**

After each solve, compute margin to critical current:

```cpp
// Local margin to quench (computed on each rank)
real margin_local = 1.0;

// For each element on this rank
for (Element* elem : local_elements) {
    J = ||curl(H)||                      // Current density magnitude
    J_c = J_c(B, T)                      // Critical current (field/temp dependent)
    margin = 1 - J/J_c                   // Margin (0 = quench, 1 = safe)
    margin_local = std::min(margin_local, margin);
}

// Global reduction: minimum across all ranks
if ( comm_size() > 1 ) {
    Cell<real> all_margins;
    collect(all_margins, margin_local);  // Gather to rank 0

    if ( comm_rank() == 0 ) {
        margin_min = min(all_margins);    // Compute global minimum
    }

    broadcast(margin_min, 0);             // Broadcast to all ranks
} else {
    margin_min = margin_local;
}
```

**Step 2: Adjust time step based on proximity**

```cpp
if (margin_min < margin_threshold) {
    // Approaching quench - reduce Δt aggressively
    quench_factor = max(quench_min, margin_min / margin_threshold)
    Δt_new = Δt_old × quench_factor

    // Tighten iteration target
    iter_target = 5   // Not 20!

} else {
    // Standard adaptive stepping (existing algorithm)
    iter_ratio = iter_target / iter_actual
    Δt_new = Δt_old × iter_ratio^α
}
```

**Recommended parameters (from experience + Bathe guidelines):**

```cpp
margin_threshold = 0.3      // Start reducing Δt at 30% margin
quench_min = 0.05          // Allow Δt reduction to 5% of original
iter_target_normal = 10    // Standard: 10 iterations (NOT 20!)
iter_target_quench = 5     // Near quench: aggressive control
α = 0.7                    // Damping (prevent oscillations)
```

### Why This Works

**Scaling argument:**

Near quench, material response time τ_mat scales as:

```
τ_mat ∝ (1 - J/J_c)   [margin to quench]
```

For convergence in k iterations:

```
Δt / τ_mat ≈ k × convergence_rate
```

If τ_mat shrinks by 10×, and we keep Δt constant:
- Iterations increase by 10× (observed: 10 → 100)

If we reduce Δt by 10× to match τ_mat:
- Iterations return to baseline (~10)

**This is temporal error control, not nonlinear solver tuning.**

---

## Code Modification Points

### Priority 1: Add Quench Detection

**File:** `src/fem/kernel/cl_FEM_Controller.hpp`

**Add member variables (after line 130):**

```cpp
// Quench detection
real mMarginToQuench = 1.0;          // Minimum margin 1 - J/J_c
real mMarginThreshold = 0.3;          // Threshold for Δt reduction
real mQuenchMinFactor = 0.05;         // Minimum Δt reduction factor
uint mIterationTargetQuench = 5;      // Tight control near quench
```

**Add method declaration (after line 244):**

```cpp
void
compute_quench_proximity();
```

---

**File:** `src/fem/kernel/cl_FEM_Controller.cpp`

**Implement quench detection (new function, add after `adjust_timestep()`):**

```cpp
void
Controller::compute_quench_proximity()
{
    // Only for electromagnetic kernel
    if ( mKernel == nullptr ) return;

    mMarginToQuench = 1.0;  // Reset to safe (local minimum)

    // Get blocks
    const Cell< mesh::Block * > & tBlocks = mKernel->mesh()->blocks();

    for ( mesh::Block * tBlock : tBlocks )
    {
        for ( mesh::Element * tElement : tBlock->elements() )
        {
            // Compute J = ||curl(H)|| at element
            // (Implementation depends on field storage structure)
            real tJ = compute_current_density_magnitude( tElement );

            // Compute J_c(B, T) at element
            real tJc = compute_critical_current( tElement );

            // Margin to quench
            real tMargin = 1.0 - tJ / tJc;

            // Track local minimum
            mMarginToQuench = std::min( mMarginToQuench, tMargin );
        }
    }

    // MPI reduction: minimum across all ranks
    if ( comm_size() > 1 )
    {
        // Gather all local minima to rank 0
        Cell<real> tAllMargins;
        collect( tAllMargins, mMarginToQuench );

        if ( comm_rank() == 0 )
        {
            // Find global minimum
            mMarginToQuench = min( tAllMargins );  // linalg::min()
        }

        // Broadcast global minimum to all ranks
        broadcast( mMarginToQuench, 0 );
    }
}
```

**Note:** The helper functions `compute_current_density_magnitude()` and `compute_critical_current()` need to be implemented based on your field storage. Check how postprocessing computes J and J_c.

---

**Modify `adjust_timestep()` (around line 729):**

**BEFORE (current code, line 732-776):**
```cpp
void Controller::adjust_timestep()
{
    // Ensure adjustments only happen after a couple of iterations
    if ( mIteration0 > 0 )
    {
        // === Temporal error proxy ===
        // ... existing code ...

        // === Iteration-based factor ===
        real tIterRatio = static_cast<real>(mIterationTarget)/mIteration ;

        // ... rest of function ...
    }
}
```

**AFTER (with quench detection):**
```cpp
void Controller::adjust_timestep()
{
    // Ensure adjustments only happen after a couple of iterations
    if ( mIteration0 > 0 )
    {
        // === Quench proximity check ===
        this->compute_quench_proximity();

        // === Temporal error proxy ===
        // ... existing code unchanged ...

        // === Iteration-based factor ===
        uint tIterTarget = mIterationTarget;  // Use local variable

        // Override target if approaching quench
        if ( mMarginToQuench < mMarginThreshold )
        {
            tIterTarget = mIterationTargetQuench;  // Aggressive control

            if ( mCommRank == 0 && gLog.info_level() > 1 )
            {
                message( InfoLevel::Verbose,
                    "    Quench proximity: margin = %.3f, reducing iter target to %u",
                    mMarginToQuench, tIterTarget );
            }
        }

        real tIterRatio = static_cast<real>(tIterTarget)/mIteration ;

        // Iteration-based factor
        real tIterFactor = std::pow( tIterRatio, 0.5 ) ;

        // === Quench-based Δt scaling ===
        real tQuenchFactor = 1.0;
        if ( mMarginToQuench < mMarginThreshold )
        {
            // Scale Δt proportional to margin
            tQuenchFactor = std::max( mQuenchMinFactor,
                                     mMarginToQuench / mMarginThreshold );

            if ( mCommRank == 0 && gLog.info_level() > 1 )
            {
                message( InfoLevel::Verbose,
                    "    Applying quench factor: %.4f (margin = %.3f)",
                    tQuenchFactor, mMarginToQuench );
            }
        }

        // Combined with geometric mean and safety
        real tPhi = tSafety * std::sqrt( tErrFactor * tIterFactor ) * tQuenchFactor;

        // ... rest of function unchanged ...
    }
}
```

---

### Priority 2: Lower Default Iteration Target

**File:** `src/fem/kernel/cl_FEM_Controller.hpp` (line 58)

**BEFORE:**
```cpp
uint mIterationTarget = 20 ;
```

**AFTER:**
```cpp
uint mIterationTarget = 10 ;   // Bathe §8.4.1: faster Δt adaptation
```

**Rationale:** From Bathe (p. 758), targeting 5-10 iterations allows faster time step adaptation. Current target of 20 is too slow to react to changing material response.

---

### Priority 3: Add Diagnostic Output

**File:** `src/fem/kernel/cl_FEM_Controller.cpp`

**Modify `print_footer()` (around line 987) to include quench margin:**

**Add after line 1020 (before closing box):**

```cpp
if ( mMarginToQuench < 1.0 )  // Only print if computed
{
    tFormat = "      │                      Minimum margin to quench     :         %6.3f  │" ;
    tMessage = sprint( tFormat.c_str(), mMarginToQuench );
    std::cout << tMessage << std::endl ;
}
```

---

## Alternative Approaches Considered (And Why They're Lower Priority)

### Anderson Acceleration

**What:** Accelerate Picard fixed-point iterations using quasi-Newton subspace methods.

**Theory (Walker & Ni, SIAM J. Sci. Comput. 2011):**

```
x^(k+1) = (1 - Σθ_i) G(x^k) + Σ θ_i G(x^(k-i))

where θ minimizes ||F(x^(k+1))||
```

**Pros:**
- 2-5× faster convergence for smooth nonlinearities
- No Jacobian needed (cheaper than Newton)

**Cons:**
- Requires storing m=3-5 previous iterates (memory cost)
- Unstable near singularities (quench)
- Doesn't solve root problem (Δt too large)
- Implementation complexity: ~500 lines

**Verdict:** Academic interest only. Won't help near quench.

---

### Trust Region Methods

**What:** Constrain Newton step to ||Δx|| ≤ δ where δ is "trust radius."

**Theory (Nocedal & Wright Ch 4):**

```
min_Δx  ||F(x) + J·Δx||    subject to  ||Δx|| ≤ δ

Adjust δ based on ρ = actual_reduction / predicted_reduction
```

**Pros:**
- Globally convergent (unlike line search)
- Better for ill-conditioned problems
- Automatically limits step size

**Cons:**
- Very complex implementation (~1000 lines)
- Requires solving constrained QP each iteration
- Expensive for large FEM systems
- Still doesn't address Δt discretization error

**Verdict:** Overkill for time-dependent problems. Use adaptive Δt instead.

---

### Continuation in Temperature

**What:** Control heating rate instead of current when near quench.

**Concept:**

```cpp
if (margin_min < 0.2) {
    // Limit temperature rise per step
    ΔT_max = 0.1;  // K/step
    Δt_new = Δt_old × (ΔT_max / ΔT_actual)
}
```

**Pros:**
- Physical interpretation (quench is thermal runaway)
- Directly controls the critical parameter

**Cons:**
- Requires thermal solve (if segregated coupling)
- Complex for fully-coupled electromagnetic-thermal
- Margin-based Δt reduction is simpler and more general

**Verdict:** Interesting for quench prediction, but more complex than Priority 1.

---

## Literature References

### Primary References (Authoritative FEM Textbooks)

**Bathe (2016)** - Finite Element Procedures (2nd Ed.)
- §8.4.1 (pp. 755-759): Newton-Raphson methods, convergence properties, **"decrease the load step"** recommendation
- §8.4.2 (pp. 759-761): BFGS method with line search
- §8.4.3 (pp. 761-763): Arc-length methods for structural collapse
- §2.5 (pp. 67-70): Convergence order definitions (linear, quadratic)

**Zienkiewicz & Taylor (2014)** - The Finite Element Method Vol 2: Solid and Structural Mechanics (7th Ed.)
- Ch 3: Nonlinear problems, material nonlinearity vs. geometric nonlinearity
- Discussion of automatic time stepping for severe nonlinearity

**Nocedal & Wright (2006)** - Numerical Optimization (2nd Ed.)
- Ch 3: Line search methods, Armijo condition
- Ch 4: Trust region methods
- (Not in BELFEM literature/ but standard reference)

### Supporting References

**Walker & Ni (2011)** - "Anderson Acceleration for Fixed-Point Iterations", SIAM J. Sci. Comput. 33(4)
- Modern acceleration technique for Picard iterations
- (Not in BELFEM literature/ but widely cited)

**Crisfield (1981)** - "A Fast Incremental/Iterative Solution Procedure", Computers & Structures
- Arc-length method developments
- Referenced by Bathe §8.4.3

---

## Implementation Roadmap

### Phase 1: Quench Detection (1-2 days)

**Goal:** Add margin calculation and diagnostic output.

**Tasks:**
1. Implement `compute_quench_proximity()` (requires understanding field storage)
2. Add member variables to Controller
3. Integrate into `print_footer()` for visibility

**Deliverable:** Diagnostic output showing margin to quench at each time step.

---

### Phase 2: Adaptive Δt Reduction (0.5 days)

**Goal:** Automatically reduce time step when approaching quench.

**Tasks:**
1. Modify `adjust_timestep()` to include quench factor
2. Set parameters: `margin_threshold = 0.3`, `quench_min = 0.05`
3. Test on known near-quench case

**Expected result:** Iteration count stays under 20 even near quench.

---

### Phase 3: Iteration Target Tuning (0.5 days)

**Goal:** Optimize target iteration counts.

**Tasks:**
1. Change default `mIterationTarget` from 20 → 10
2. Add `mIterationTargetQuench = 5` for aggressive control
3. Parameter sweep to validate settings

**Expected result:** 30-50% speedup from faster Δt adaptation.

---

### Phase 4: Testing & Validation (1-2 days)

**Test cases:**
1. Ramp current to quench in coil geometry
2. AC loss simulation with high fields
3. Verify iteration counts: normal < 10, quench approach < 20

**Success criteria:**
- No case exceeds 30 iterations (currently 100)
- Total compute time reduced by 5-10×
- Time step adaptation is smooth (no oscillations)

---

## Open Questions

### Q1: How to compute J and J_c in Controller?

**Current status:** Postprocessing computes these, but Controller may not have direct access.

**Options:**
1. Call postprocessor methods (if accessible)
2. Duplicate computation (simple but inefficient)
3. Store margin during postprocessing, read in Controller

**Recommendation:** Check if `mEquation->compute_current_density()` exists or similar IWG method.

---

### Q2: Should quench detection be physics-specific?

**Issue:** This implementation assumes Maxwell/HTS physics. Other physics (structural, thermal) don't have J_c.

**Options:**
1. Add virtual method to IWG: `virtual real proximity_to_singularity()`
2. Make quench detection conditional on equation type
3. Generalize to "material stiffness indicator"

**Recommendation:** Start with Maxwell-specific, generalize later if needed.

---

### Q3: MPI communication cost for margin reduction?

**Issue:** `allreduce(MPI_MIN)` every time step.

**Impact:** Negligible (single real, synchronous time-stepping already requires barriers)

**Optimization:** Could batch with existing MPI calls if profiling shows cost.

---

## Expected Impact

### Performance

**Current:** ~100 iterations near quench @ Δt = 1e-3 s
- Time per step: ~100 × (assembly + solve)
- Total time to quench: N_steps × 100 × cost_per_iter

**After Priority 1+2:** ~10-15 iterations @ Δt = 1e-4 s (10× smaller)
- Time per step: ~10 × (assembly + solve)
- Total steps: 10× more steps, but 10× fewer iterations each
- **Net speedup: 5-10×** (fewer iterations dominates more steps)

### Robustness

**Current:** Manual time step selection, trial-and-error
**After:** Automatic adaptation, no user tuning needed

### Code Complexity

**Lines added:** ~150-200 lines
**Files modified:** 2 (Controller.hpp, Controller.cpp)
**External dependencies:** None
**Risk:** Low (isolated changes, easy to revert)

---

## Conclusion

**Line search is not the solution.** The fundamental issue is temporal discretization error when Δt doesn't match the physics time scale.

**The literature (Bathe, Zienkiewicz) is clear:** When facing severe nonlinearity, reduce the time/load step. Adaptive time-stepping with quench detection directly addresses the root cause.

**Recommended priority:**
1. ✅ Implement quench-aware adaptive time-stepping (Priority 1+2)
2. ⏸️ Defer line search / Anderson / trust region (won't help near singularity)
3. 🔬 Consider temperature continuation only if quench prediction is a specific goal

**Expected outcome:** 5-10× speedup with ~150 lines of code and solid theoretical foundation from authoritative FEM literature.
