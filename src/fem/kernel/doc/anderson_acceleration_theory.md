# Anderson-Accelerated Picard Iteration {#fem_kernel_anderson_acceleration_theory}

**Date:** 2026-07-31
**Purpose:** Theory and implementation reference for the Anderson mixing of the
Picard branch (opt-in: off by default, master switch
`timestep { anderson stabilization : true ; }`): formulation, safeguards, the
stage/commit handshake with the nonlinear controller, and the
history-invalidation rules.
**Module:** `src/fem/kernel` (`fn_FEM_anderson_mixing.hpp`, `cl_FEM_DofMgr_SolverData`,
`cl_FEM_Controller`)

---

## 1. Motivation

The Picard stage of the nonlinear iteration (Messe et al. 2023, Eq. 12) is a relaxed
fixed-point iteration

```
x_{k+1} = ( 1 − ω ) x_k + ω G( x_k ) ,      G( x ) = A(x)⁻¹ b(x)
```

whose convergence rate is bounded by the spectral radius of the relaxed iteration map. The
scalar relaxation ω rescales the step but cannot exploit the *direction* information
contained in successive residuals. For the stiff HTS power law (n ≈ 35) and partitioned
magneto-thermal coupling, this spectral-radius limit is what makes Picard phases slow.

Anderson mixing (Anderson 1965; type-II form as analyzed by Walker & Ni 2011) replaces the
scalar update with a multisecant extrapolation over a short window of previous iterates. Its
depth-1 limit is Aitken-style dynamic relaxation (Irons & Tuck), the standard accelerator for
partitioned coupling. In BELFEM it accelerates each kernel's own Picard iteration; in fully
coupled runs, where the peer field moves between outer iterations, this is *block-Anderson* —
each field is accelerated while the other evolves, the standard partitioned-coupling usage.

## 2. Formulation

Let `r_k = G(x_k) − x_k` be the **fixed-point residual** and let the window hold the m most
recent committed pairs `(x_i, r_i)` plus the current pair. Form the difference columns

```
Δx_j = x_{j+1} − x_j ,   Δr_j = r_{j+1} − r_j        ( j over consecutive chain pairs )
```

solve the small linear least-squares problem

```
γ = argmin ‖ r_k − ΔR γ ‖₂
```

and update

```
x_{k+1} = x_k + β r_k − ( ΔX + β ΔR ) γ .
```

With an empty window (or on any fallback) the update reduces **exactly** to the relaxed
Picard step — this is the regression anchor of the implementation: `anderson depth : 0`
bypasses all Anderson code and runs the plain relaxed update, written in increment form
`x -= ω δ` since the increment-form change, algebraically identical to the historical
`(1−ω) x + ω G` (`cl_FEM_DofMgr_SolverData.cpp`, Picard branch).

Throughout this document, β denotes the Anderson mixing parameter; in BELFEM it is the
controller's **live relaxation ω**. The
controller's backtracking line search therefore damps the Anderson step directly, and its
relaxation adaptation (growth on improvement, decay on worsening) remains active; the
overshoot guard specific to Anderson is the flush-on-reject rule of Section 4, not a frozen ω.
Note the asymmetry that motivates that rule: β scales the `r` and `ΔR` content of the step,
but **not** the `ΔX γ` term — a rejected mixed step cannot be rescued by halving β alone,
which is why rejection empties the window and retries with plain damped Picard.

## 3. The least-squares solve and its safeguards

The columns of ΔR go collinear precisely as the iteration converges, so the small solve is
performed by QR (`lapack::gels`) on **column-normalized** ΔR — normal equations would square
the condition number exactly where it explodes. γ is un-scaled after the solve. The solve
runs on the master rank inside the existing solve path; no additional MPI communication is
introduced (`fn_FEM_anderson_mixing.hpp`).

Safeguards (the illegal-argument abort fires immediately after the `gels` call, before any
soft handling; the soft guards then act in the order listed):

| Guard | Trigger | Action |
|---|---|---|
| window clamp | more columns than free dofs (tiny systems) | window truncated to n; always active, also in release |
| vanishing column | ‖Δr_j‖ below machine tolerance | treat as rank-deficient |
| illegal argument | `gels` reports `info < 0` | programming error — always aborts |
| rank deficiency | `gels` reports a zero diagonal (`info > 0`) | drop the **oldest** column, retry once |
| coefficient guard | non-finite γ or ‖γ‖∞ > 10² | drop the oldest column, retry once |
| final fallback | retry failed too | plain relaxed Picard step; the pair is **not** committed |

The 10² bound on the (un-scaled, unit-column) coefficients is an engineering default; the
mixing coefficients of a healthy window are O(1)–O(10), and larger values indicate the
least-squares geometry has degenerated even though `info == 0`.

## 4. Stage, commit, and flush: history validity

The Anderson window is only valid for a contiguous run of **accepted** Picard iterates from
**one** timestep attempt of **one** kernel. Two mechanisms enforce this.

**The stage/commit handshake.** The solve path cannot know whether the controller's line
search will accept the trial it just produced. `SolverData` therefore *stages* the pair
`(x_k, r_k)` (captured before the update overwrites the dofs: `x_k` = the assembly-time
free-dof values, `G(x_k)` = `x_k − δ`, reconstructed from the increment `δ` the increment-form
solve returns); the controller *commits* it into the window
only when the trial is accepted, and *discards* it on rejection. A pair whose own
least-squares solve fell through to the plain step is not staged either — with one exception:
the bootstrap pair of an empty window always stages, otherwise the window could never fill.

**The flush rules.** The window is emptied whenever the map it was built from changes:

1. at every attempt and sub-step start, and on every retry after a cut,
2. on every Picard↔Newton switch (promotion, demotion, stagnation fallback, escalation),
3. on a rejected trial — together with the discard; the retry then runs plain damped Picard,
4. when the thermal update gate freezes/thaws the thermal solve (the state moves while the
   kernel is frozen),
5. on a memory-dump restore (history is never persisted).

**Residual semantics.** The Anderson path reports the same **pre-update** ε as the depth-0
Picard path: ‖A(x_k)·x_k − b(x_k)‖/‖b‖, the published out-of-balance force criterion
(Messe et al. 2023 §4, Eq. 10–11). It must NOT evaluate the mixed state under the lagged
operator. At β = ω = 1 the bootstrap/fallback step is exactly x = A⁻¹b, so that "residual"
is the direct solver's roundoff (≈ −124 dB), independent of nonlinear consistency. The
controller then declares convergence on an under-iterated state, observed as thin-shell
checkerboarding; removing the mixing restored physical fields. An earlier revision refreshed
the field-values vector here to avoid one-iterate-late controller judgments (ts16
false-promotion cycle). The pre-update gate keeps that lag; it costs at most one extra
Picard iterate and avoids the dishonest lagged-operator residual. Per-iterate mixing quality
is exposed separately as the **fixed-point residual** ‖G(x_k)−x_k‖/‖x_k‖
(`fixed_point_residual()`, broadcast alongside ε), Walker & Ni's monitored quantity and
Bathe's (§8.4.4) increment criterion. That quantity alone can under-report true error on
stiff maps, so it never replaces ε as the stop test. The Newton branch keeps its post-update
recompute: its tangent differs from the solved operator, so the recompute is non-degenerate
there.

## 5. Cost and configuration

Memory (master rank, compact free-dof length n): (m+1) solution snapshots, (m+1) residual
snapshots, one n×m least-squares scratch, and the `gels` workspace — all preallocated members,
grown once. The per-iteration cost is O(n·m) plus one m-column QR, negligible against the
sparse factorization. The history containers are `ShiftRegister< Vector<real> >` (fixed
depth, no reallocation, index 0 = newest).

```
solver {
    timestep          { ... ; anderson stabilization : true ; }  // opt-in master switch, default false
    nonlinear         { ... ; anderson depth : 3 ; }   // magnetic; alias: nonlinear magnetic
    nonlinear thermal { ... ; anderson depth : 1 ; }   // depth 1 ≈ Aitken
}
```

depth m ∈ [0, 8]. Anderson is opt-in: absent keys leave the legacy Picard
update bit-identical. `anderson stabilization : true` fills missing depths
with magnetic 3 and thermal 1. An explicit `anderson depth` always wins,
including explicit 0. `anderson stabilization : false` hard-errors with any
explicit nonzero depth.

| Item | Location |
|---|---|
| mixing step (pure math, unit-tested) | `fn_FEM_anderson_mixing.hpp` — `anderson_mixing_step()` |
| staging, write, fixed-point residual | `SolverData::anderson_update()` |
| commit / discard / clear API | `SolverData::anderson_commit/discard/clear()` |
| controller hooks (accept/reject, flush sites) | `cl_FEM_Controller.cpp`, see `nonlinear_controller_theory.md` |
| unit tests | `tests/fem/test_AndersonMixing.cpp` (plain-step identity, one-step affine convergence, window shrink, fallback, clamp) |

## 6. Literature

- Anderson, D. G. 1965, *Iterative procedures for nonlinear integral equations*, J. ACM 12,
  doi:10.1145/321296.321305.
- Walker, H. F. & Ni, P. 2011, *Anderson acceleration for fixed-point iterations*, SIAM J.
  Numer. Anal. 49, doi:10.1137/10078356X — the type-II form implemented here.
- Irons, B. & Tuck, R. C. 1969, *A version of the Aitken accelerator for computer iteration*,
  Int. J. Numer. Meth. Eng. 1 — the depth-1 limit; standard in partitioned FSI coupling
  (Küttler & Wall 2008, doi:10.1007/s00466-008-0255-5).
- Messe et al. 2023, Section 4 — the relaxed Picard baseline (Eq. 12) and the hybrid
  iteration this accelerator plugs into.

(The Anderson/FSI references are not part of the local reference library; DOIs given for
independent access.)
