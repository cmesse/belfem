# BDF Nonlinear-Mass: Verification & History-Derivative Contract

**Date:** 2026-07-03
**Purpose:** The two forward-looking items left open after the BDF Jacobian scaling bug was
implemented and closed (`closed/bdf_jacobian_scaling_bug.md`). The α-scaling, coefficient
wiring, startup ramp, and the reset/savepoint paths are all in the tree and Codex-audited;
what remained was (1) a regression test that BDF2+ now converges like BDF1, and (2) making
the Newton *history* derivative exact for nonlinear mass at BDF ≥ 2.
**Status (2026-08-13): CLOSED — DR-33 closed on Christian's ruling.** The parked run
half is retired by execution: the tapestack3d quench campaign ran BDF1→BDF5 end to end,
coupled, with temperature-dependent ρcp and savepoint restores — the fixed-Δt Newton-count
A/B is waived with it. The nine-case unit coverage stays in `check-fast`.

Previous status (2026-08-11): V1 SPLIT — the no-deck half is done and in `check-fast`
(`tests/fem/test_BdfTimestepMethod.cpp`); the run half is PARKED because no deck exists to
run it on (Christian), its driver `thermalTest` having been deleted on 2026-08-06
(`e4c02ac8`). B8 added and fixed the same day.** B7 is done
and in tree: `phi_ferro_newton` contracts `dMdx_times_h` against
`IWG_Timestep::collect_qhist()`, and the `#HACK` marker is gone. B8 (below) found that both
Maxwell producers contracted that history as a **magnitude** where the isotropized tangent
requires a **signed projection** — fixed at both sites, run gate pending. Priority note: the 2026-08-07 opt-out change briefly made **BDF5 the default scheme**
(`dl20260807_timestep_input_optout.md`), which would have put every key-less deck on the
multi-step history path V1 guards — but it was **rolled back to BDF1 the same day** with the
Anderson residual fix (`dl20260807_greg3_false_convergence_jury.md`, addendum 1; confirmed in
tree: `cl_FEM_Controller.hpp:189`). So V1 stays P2: it guards a path decks must now select
explicitly. Tracked as `debt_register.md` DR-33.
**Severity:** Medium — the implemented fixes make BDF2-5 functional and give an exact
current-state tangent; B7 only affects the *rate* of Newton convergence for nonlinear-mass
BDF ≥ 2, not correctness of the converged solution.
**Affected Modules:** `src/fem/iwg/`, `src/fem/maxwell/matrices/`
**Module:** fem/iwg

## Background

`assemble_dJdx()` builds the Newton correction (`cl_TimestepMatrices.cpp:131`):

```
dJ/dx = α·(∂M/∂x)·x^n − (∂M/∂x)·(Σ β_i·x^{n-i}) + Δt·(∂K/∂x)·x^n − Δt·(∂f/∂x)
```

B4 fixed the α factor on the current-state term. The history term subtracts
`mdMdX_times_h`, which is *contracted* by the physics IWG — `assemble_dJdx` assumes it
already holds `(∂M/∂x)·(Σ β_i·x^{n-i})`. See the residual/Jacobian derivation in
`closed/bdf_jacobian_scaling_bug.md` ("Theoretical Background") and Hairer & Wanner (1996)
II.4.

## Work Items

- [ ] **V1 — regression test (was B5). BLOCKED 2026-08-11: the driver it names is gone.**
  `src/executables/thermalTest.cpp` was deleted on 2026-08-06 (`e4c02ac8`, dead-code
  cleanup alongside `periodictest.cpp`); `src/executables/` now holds only `hphirun`,
  `hphiTrun` and `electricalCircuit`. Reviving it verbatim would revive dead code: the
  recovered file hardcodes `costheta.msh` and `SolverType::UMFPACK` (opt-in, OFF by
  default), starts from T = 0 K, and carries its own
  `// todo: implement cp and lambda in material definition` — so the ρ(T) nonlinear mass
  this item is named for was never configured in it. **A deck decision is owed before
  this can be a run.** Original intent, unchanged: compare BDF2-vs-BDF1 Newton-iteration
  count on a nonlinear heat problem with ρ(T) temperature-dependent mass. Success:
  BDF2 iterations ≈ BDF1 (±1) at a fixed Δt, where before the fix BDF2 produced an
  all-NaN system. Optionally add the HTS benchmark (nonlinear μ_r(H)) comparison and
  check the segregated-path magnetic-reset-after-thermal-substeps case exercises the new
  savepoint restore (`Controller::reset_timestep`).
  - **The item splits, and the cheap half is DONE (2026-08-11):** `tests/fem/test_BdfTimestepMethod.cpp`,
    six cases in the `fem`/`fast` suite. `IWG_Timestep` is concrete, sits under no pure
    virtuals, and its constructor allocates only its own coefficient buffers, so it stands
    up with no mesh, kernel or dof manager. Pinned: the constructor's BDF1 default (the
    "dispatch pointer is never null" guarantee), method round-trip over BDF1-5, the
    **startup ramp holding step one at order 1 for every BDF2-5 run** — which matters
    because `mH` is all zeros before the first shift, so a broken ramp walks straight into
    `compute_bdf_coefficients()`'s own *"Invalid BDF history step size"* assert — the
    non-BDF gate that protects the eigenvalue `MassOnly`/`StiffnessOnly` switch from stale
    order state (C1), and the CN/Galerkin hard rejection, which had never been pinned
    despite being the reason that clause of DR-33 was struck rather than fixed.
  - **What the cheap half does NOT cover, and why — correcting an earlier claim in this
    file.** The variable-step coefficient formulas (`compute_bdf_coeffs_2..5`) are *not*
    reachable: `mAlpha`, `mBeta`, `mH` and `mStepCount` are **private** with no accessor,
    and the only writers of `mH`/`mStepCount` are `shift_fields`, `reset_fields` and
    `restore_savepoint`, all of which dereference `mField` and call
    `DofManager::initialize()`. Covering them needs either a dof-manager fixture or a
    deliberate seam (a `friend` probe, or public const accessors plus a way to seed the
    step-size history) — **a production-header change, so Christian's call, not assumed.**
  - **The Newton-iteration-count half is PARKED, not pending (Christian, 2026-08-11):
    there is no deck to run it on, and none is planned.** A nonlinear-mass ρ(T) thermal
    transient does not exist in the tree in any form, so this is not schedulable and is
    not carried as an open action. **Unblock condition, written down so it is
    recognisable if it ever arrives:** any deck with a temperature-dependent ρcp thermal
    transient, run once at BDF1 and once at BDF2 at fixed Δt, comparing Newton iteration
    counts (±1). No new tooling is needed — `method :` is already an input key.
  - **The coverage fact this leaves standing: nothing in the tree exercises BDF2-5 end to
    end.** Bounded, and it should not be overstated — BDF2-5 is opt-**in** (default rolled
    back to BDF1 on 2026-08-07), so no deck in use depends on it, and the original all-NaN
    defect is fixed with its startup ramp now unit-tested. What no test and no run touches
    is `compute_bdf_coeffs_2..5`, the variable-step coefficient formulas.
  - [x] **The seam is BUILT (2026-08-11, Christian's go-ahead).** One
    `friend class BdfCoefficientProbe ;` line in `cl_IWG_Timestep.hpp` — zero runtime
    cost, no widened public API, the probe defined only in the test and used nowhere in
    `src/` — plus §4 of `tests/fem/test_BdfTimestepMethod.cpp`:
    - `ConstantStepReducesToTextbook` — α and β against the published BDF2-5 constants
      (3/2, 11/6, 25/12, 137/60 and their β sets; Hairer & Wanner 1996 II.4, an
      independent reference rather than a restatement of the implementation).
    - `VariableStepIsExactOnPolynomials` — the property that actually *defines* BDF-p,
      and the one the constant-step check cannot reach: at non-uniform step size the
      formula must still differentiate every polynomial of degree ≤ p exactly. This is
      what a mistyped cumulative sum would break while leaving the uniform case intact.
    - `RampUsesReducedOrderCoefficients` — a BDF5 run two steps in must produce BDF2's
      coefficients. `collect_qhist` truncates the history to `mOrderActive`, so an α
      computed at the configured order would be contracted against a shorter history and
      the Jacobian would silently stop being the tangent of its residual.
  - **Running it corrected the test itself, and the reason is worth keeping.** The
    exactness case was first written with realistic 3e-4 steps, where it was very nearly
    vacuous: a scheme of the wrong order misses by O(h^p), which at that size is ~1e-13 —
    indistinguishable from the roundoff of a correct answer. Measured: BDF5 reproduced
    d/dt of t⁶ to 2e-13, so the test would have passed for a scheme it exists to reject.
    Steps rescaled to O(0.1), where the separation is fifteen orders of magnitude (~1e-15
    exact against 16–73 % wrong), and the **negation at degree p+1 is now asserted**
    rather than merely observed — "exact up to degree p" is a fingerprint of BDF-p only
    together with its own failure at p+1.

- [x] **B7 — `dMdX_times_h` producer contract for BDF ≥ 2. IMPLEMENTED 2026-07-14
  (Claude/Fable, Christian's go-ahead via the maxwell_kernel_collapse Q2 discussion;
  Codex audit launched):** new single-source-of-truth `IWG_Timestep::collect_qhist()`
  builds `β₀q₀ − β₁q₁ + β₂q₂ − β₃q₃ + β₄q₄` truncated to `mOrderActive` (a new member set
  by `compute_bdf_coefficients`, so the startup ramp and variable Δt are respected;
  order ≤ 1 → `q₀` exactly, preserving BDF1 bit-identity). The four `bdf2-5` RHS bodies
  now consume the same helper (their inline sign patterns deduplicated — FP-identical,
  same expressions relocated), and `phi_ferro` contracts `dMdx_times_h` against
  `|B·qhist|` instead of `|B·phi0|`; the `#HACK` marker is gone. Note for the audit: the
  dof-ordered qhist against the shape-ordered B rests on the same node↔dof ordering
  invariant the residual assembly itself uses (`aJ·tSwap` in every bdfN body).
  **Codex audit 2026-07-14 (high confidence): no blocking defect.** FP-identity of the
  bdf2-5 dedup verified; `mOrderActive` coverage verified on all dirty-marking paths incl.
  the eigenvalue MassOnly/StiffnessOnly switch; **BDF1 identity PROVEN** — `qold(0)` and
  `node_data("phi0")` read the same field entries because node dofs store `node->index()`
  as their field index (`cl_FEM_Dof.cpp:20`) and element nodal linking is node-ordered
  (`cl_FEM_Element.cpp:223`); `qswap` aliasing safe (no other writer). One medium watch
  item, not blocking: `DofManager::compute_full_matrices()`
  (`cl_FEM_DofManager.cpp:817,860`) calls `compute_mkf` directly, bypassing the
  coefficient refresh — a caller before normal assembly could give `phi_ferro` a stale
  `mOrderActive`; harmless today (in-tree `save_IV` callers run post-solve, and that path
  consumes only M/K, discarding the dMdx blocks), but keep in mind if
  `compute_full_matrices` ever gains a pre-solve caller.
  Original item:** The only producer,
  `phi_ferro()` (`src/fem/maxwell/matrices/mt_maxwell_phi.cpp:90`, marked
  `#HACK: TO BE TESTED`), contracts `∂M/∂x` with `phi0` (the single previous state)
  rather than the β-weighted history sum `Σ β_i·x^{n-i}`. Exact for BDF1 (β_0 = 1, one
  history state), an unweighted single-state approximation for BDF2-5. Consequence:
  inexact (but typically still convergent) Newton for nonlinear-mass BDF ≥ 2.
  - Proposed approach: precompute the β-weighted history sum (the same `tSwap` vector the
    `bdfN` RHS assembly already builds, `cl_IWG_Timestep.cpp:499-509`) into a
    calculator-visible buffer (e.g. `mCalc->qswap()`) before `compute_mkf`, and have the
    matrix producer contract `∂M/∂x` against that buffer for the history term instead of
    `phi0`. Verify the current-state term still uses `x^n`.
  - Cross-check against the residual derivation before implementing; the history term is
    the only place the multi-state weighting enters the tangent.

- [x] **B8 — signed history projection, not a magnitude. FIXED 2026-08-11 (Claude, from Christian's
  task brief; independently confirmed by Codex and Grok).**
  B7 moved the producers to the *right vector* (`collect_qhist`), but both Maxwell
  sites still reduced it to the *wrong scalar*: `phi_ferro_newton` used
  `h0 = |B·qhist|`, and the nonlinear-μ branch of `h_newton_mu` used
  `h0 = |E·qhist|`. The isotropized tangent instead requires the signed projection
  `h0 = ĥ·Hvec_hist = dot(Hvec, Hvec_hist)/H` (slot contract
  `cl_TimestepMatrices.hpp` `(dM_ik/dx_j)·h_k`; the thermal producer's
  `dot(Nvec, qhist)` was already signed). The two agree only for parallel, co-oriented
  history — under field reversal (AC zero-crossing) the magnitude flips the sign of the
  tangent's history term. Both sites now compute the dot-product form with shared scratch
  (`"hcur"`/`"hhist"`, registered in `IWG_Maxwell::create_custom_vectors_and_matrices`,
  with no per-intpoint allocation). At `H ≤ BELFEM_EPSILON`, both sites set the history
  scalar to 0: a zero-subgradient choice that degrades the tangent to the Picard matrix,
  the correct limit. The residual is untouched, so its roots are invariant, although
  iteration paths may differ. The remaining gate is a Newton-iteration-count A/B run on
  a ferro AC deck through a zero-crossing (expect improvement), with a monotonic-ramp
  control (expect ~no change). Derivation and H→0 decision:
  `devlog/dl20260811_ferro_history_scalar_sign.md`. Tracked as `debt_register.md` DR-64.

## Notes / Constraints

- The current-state α factor (B4) is already correct and Codex-confirmed; do not touch it.
- BDF1 must remain exact — the `phi0` contraction is the correct BDF1 history term, so any
  change must reduce to `phi0` when the effective order is 1 (including the startup ramp's
  first step of a BDF-p run).
- `phi_ferro` also feeds `dMdX_times_x` from `H1` (current state) — that path is unaffected.

## Success Criteria

1. **V1:** BDF2 Newton-iteration count within ±1 of BDF1 on `thermalTest`; no NaN.
2. **B7:** for a nonlinear-mass BDF2 step, the assembled `dJdx` history block matches a
   finite-difference tangent to the RHS history term within solver tolerance.
3. **No regression:** BDF1 behaviour unchanged (bit-for-bit where feasible).

## References

- Predecessor (implemented + closed): `closed/bdf_jacobian_scaling_bug.md`
- Session record: `devlog/dl20260703_bdf_jacobian_fix.md`
- Hairer & Wanner (1996) II.4 — BDF Jacobian structure
- Messe et al. 2023 (paper1) §2.7 — HTS convergence tolerance

---

**Status:** OPEN — implementation of the core BDF fix is done and audited; this is the
verification (V1) and the one exactness follow-up (B7).
**Assigned:** TBD (V1 needs a build + run; B7 is a producer-side change)
