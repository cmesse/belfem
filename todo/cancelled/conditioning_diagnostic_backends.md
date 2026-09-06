# Conditioning Diagnostic: Backends Beyond MUMPS

> **CANCELLED 2026-09-03** (todo/ currentness sweep, round 3): low relevance — Christian had already ruled STRUMPACK, PARDISO, UMFPACK and PETSc-without-MUMPS unsupported; the remaining SuperLU and PETSc-via-MUMPS conditioning backends serve a diagnostic nobody is blocked on. Status lines and checkboxes below are as they stood at closure and are not maintained.

**Date:** 2026-08-10
**Purpose:** Give the `compute conditioning` diagnostic a solver-native backend beyond MUMPS,
using the Arioli/Demmel/Duff error analysis those libraries already implement, so it stops
depending on an eigenvalue estimate that cannot resolve a badly conditioned mixed h-phi system.
**Module:** `src/sparse` (+ `src/fem/kernel`)
**Status:** IN PROGRESS ( R1-R2 of R6 ). MUMPS forward/backward-error accessors landed 2026-08-10.
The eigenvalue fallback stays LIVE ( decided 2026-08-10, Christian, after a brief deactivation the
same day ): it is still correct wherever kappa is moderate, so it degrades to a warning and a NaN
instead of being removed. The deck parser now warns that conditioning without MUMPS is not
recommended.
STRUMPACK, PARDISO, UMFPACK and PETSc-without-MUMPS are accepted as unsupported ( Christian,
2026-08-10 ) -- the number cannot be had from them for reasonable effort.
Remaining: SuperLU ( R3 ), PETSc-via-MUMPS ( R4 ), Controller routing ( R5 ), docs ( R6 ).

> **Scope guards:**
> - `EigenValues` is NOT deleted and NOT disabled. It is intact, computes the large end correctly,
>   remains the basis for eigenvalue work and for generalized problems later, and stays the
>   fallback for solvers with no native error analysis. It is simply not a reliable kappa
>   estimator for a badly conditioned mixed formulation.
> - The target quantity is the ADD (Arioli, Demmel & Duff) linear-system condition number and its
>   companions, NOT the spectral ratio.

---

## 1. Why the eigenvalue path is not enough on its own

`compute_conditioning()` returned `|lambda_max| / |lambda_min|`. ARPACK accepts a Ritz value when
its error bound falls below `tol * |lambda|`, and that bound cannot go below the backward error of
the matvec, roughly `eps_mach * |lambda_max|`. So the small end is reachable only while

```
tol  >  eps_mach * kappa        ( ~ 1e-16 * kappa )
```

| kappa | tol needed | reachable? |
|---|---|---|
| 1e9  | 1e-7 | yes |
| 1e11 | 1e-5 | yes, marginal |
| 1e12 | 1e-4 | yes ( current default ) |
| **1e17** | **> 1** | **no — asks for >100 % relative error** |

A mixed h-phi formulation reaches kappa ~ 1e17 by construction — the H and phi blocks and the
constraint rows carry very different scalings. That is characteristic of the formulation, **not**
numerical singularity (Christian, 2026-08-10; cf. Dular et al. 2021 on mixed h-phi stability).

**Measured 2026-08-10:** with tol = 1e-4 and ncv = 20 the run exhausts the restart limit and
converges nothing. Confidence: high — the mechanism is arithmetic, and the measurement agrees.

Raising `maxit`, `ncv` or loosening `tol` further cannot help.

> **RECONCILED 2026-08-28.** The sentence that followed here — "a spectral transform
> (shift-and-invert) would be needed, and that reintroduces the factorization-reuse problem" —
> was **half right and its conclusion was wrong**. A spectral transform is indeed needed; it does
> **not** have to be shift-and-invert, and therefore does not have to reuse a factorization. A
> **spectral fold**, `OP = sigma*I - A` with `sigma` just above `max|lambda|`, is a spectral
> transform costing one extra matvec term per product: it turns the small end of `A` into the
> large end of `OP`, which is the end Arnoldi converges. No factorization, no JOB 5/6 problem, no
> new TPL. Implemented 2026-08-28 — see `todo/cancelled/arpack_small_end_configuration.md`.
>
> Two limits from that work belong here, because they bound what the fallback can ever report:
> `lambda_min` comes out of a cancellation, so its relative error is about
> `max( tol, eps_mach ) * kappa` — the same arithmetic wall as the table above, now reached from
> the other side — and the fold declines outright on a spectrum straddling zero, where the
> smallest-magnitude eigenvalue is interior and no fold reaches it. So the case for a
> solver-native ADD number (R3-R5 below) is unchanged; what changed is that the fallback now
> answers where it can, instead of exhausting the restart limit first.

**The path is nonetheless kept.** Where kappa is moderate the estimate is perfectly good, and it
is the only option for a solver that supplies nothing else. What changed is the failure mode:
`EigenValues` warns and returns `BELFEM_QUIET_NAN` rather than aborting, `arpack::check_naupd`
treats `-9999` and `check_neupd` treats `-14` as expected algorithmic outcomes rather than errors,
and the footer shows `n/a` for that step. This todo is about giving the other solvers something
better, not about removing the fallback.

## 2. What MUMPS already gives — and what BELFEM throws away

`mRInfoG` is sized 20 (`cl_SolverMUMPS.cpp:38`) and MUMPS fills the whole error-analysis block
under `ICNTL(11) = 1`, which `Controller::arm_conditioning_*` already sets. Only two entries have
accessors:

| RINFOG | `mRInfoG` idx | quantity | exposed today |
|---|---|---|---|
| 4 | 3 | `\|\|A\|\|_inf` | no |
| 5 | 4 | `\|\|x\|\|_inf` | no |
| 6 | 5 | scaled residual | no |
| 7 | 6 | omega1, componentwise backward error | no |
| 8 | 7 | omega2 | no |
| **9** | **8** | **forward error bound `\|\|dx\|\|/\|\|x\|\|`** | **no** |
| 10 | 9 | COND1 | `get_cond0()` (`:1269`) |
| 11 | 10 | COND2 | `get_cond1()` (`:1280`) |

**RINFOG(9) is arguably the better diagnostic for "how hard is this to solve":** kappa bounds the
worst case, the forward error bound says how many digits of `x` were actually lost. It costs
nothing — the value is already in memory.

## 3. PETSc — via its MUMPS interface

Verified in the linked headers (`/opt/scls/mkl/include/petscmat.h`):

```
:2401  MatMumpsSetIcntl ( Mat, PetscInt, PetscInt )
:2407  MatMumpsGetInfog ( Mat, PetscInt, PetscInt * )
:2408  MatMumpsGetRinfo ( Mat, PetscInt, PetscReal * )
:2409  MatMumpsGetRinfog( Mat, PetscInt, PetscReal * )
```

So PETSc does not compute ADD itself — it **forwards MUMPS's**. The route is:

1. `PCFactorSetMatSolverType( pc, MATSOLVERMUMPS )` — PETSc here has MUMPS
   (`petscconf.h:106`)
2. `PCFactorGetMatrix( pc, &F )` to reach the factored matrix
3. `MatMumpsSetIcntl( F, 11, 1 )` to arm the error analysis
4. `MatMumpsGetRinfog( F, 9 / 10 / 11, &value )` for the forward error bound and COND1/COND2

Same numbers, same semantics as the direct MUMPS path, so both backends would report a
comparable quantity. Note it requires the PETSc solve to be `PCLU` with MUMPS — it yields nothing
for an iterative PETSc configuration.

*(PETSc's other option, `KSPComputeExtremeSingularValues` (`petscksp.h:325`), is a different
quantity — sigma_max/sigma_min of the PRECONDITIONED operator — and needs a Krylov method, so it
is unavailable in the `KSPPREONLY` + LU configuration BELFEM selects at
`cl_SolverPETSC.cpp:479-484`.)*

## 4. SuperLU — native, and it fits what BELFEM already has

Verified in `slu_ddefs.h`:

```
:110  dgssvx( options, A, perm_c, perm_r, etree, equed, R, C, L, U, work, lwork,
              B, X, recip_pivot_growth, rcond, ferr, berr,
              Glu, mem_usage, stat, info )
:210  dgscon( norm, L, U, anorm, rcond, stat, info )
:214  dgsrfs( trans, A, L, U, perm_c, perm_r, equed, R, C, B, X, ferr, berr, stat, info )
```

Two independent routes, and the cheap one fits BELFEM's current structure:

- **`dgscon`** estimates `rcond` (1-norm or inf-norm) **from the existing L and U factors**.
  BELFEM's wrapper already calls `dgstrf` directly and retains `mL` / `mU`
  (`cl_SolverSUPERLU.cpp:207-232`), so this is one extra call plus a work array — no second
  factorization, no expert driver, no restructuring.
- **`dgsrfs`** performs iterative refinement and returns `ferr` / `berr` — the forward and
  backward error bounds, i.e. the ADD quantities proper, directly comparable to MUMPS RINFOG(9)
  and RINFOG(7)/(8).
- `dgssvx` bundles both but would mean replacing the hand-rolled `dgstrf` path.

`dgscon` gives a norm-based kappa; `dgsrfs` gives the ADD pair. **They are different quantities** —
see O1.

## 5. Ordered Steps

- [x] **R1** — expose the MUMPS block already in memory. DONE 2026-08-10: `get_forward_error()`
      ( RINFOG(9) ) and `get_backward_error()` ( max of RINFOG(7)/(8) ) added as virtuals on
      `solver::Wrapper` and overridden in `MUMPS`. Both return a signalling NaN unless
      `ICNTL(11) = 1`, matching `get_cond0/1`.
- [x] **R2** — O1 RESOLVED 2026-08-10 ( Christian ): report BOTH. The condition number is the
      more descriptive figure for the user-facing footer ; the forward error is the cheaper thing
      to steer on and gets its own getter so a timestepper can read it.
- [ ] **R3** — SuperLU, now that O1 wants both numbers:
      - `dgscon( norm, mL, mU, anorm, &rcond, &stat, &info )` for the condition number. The
        factors are already retained ( `cl_SolverSUPERLU.hpp:66-67` ), so this needs only an
        `anorm` of the matrix and a stat object. Contained.
      - `dgsrfs` for `ferr` / `berr`. Two constraints found while scoping: the wrapper keeps no
        equilibration state, so it must pass `equed = "N"` with unused `R`/`C` ; and `B`/`X` only
        exist as `SuperMatrix` inside `solve( Vector, Vector )`, so the call belongs there.
      - **`dgsrfs` performs iterative refinement, i.e. extra triangular solves on every call.**
        It must be gated behind an error-analysis switch like MUMPS's
        `set_mumps_error_analysis` ; the SuperLU wrapper has no such switch yet. Adding it is
        part of this step.
- [ ] **R4** — PETSc: `MatMumpsGetRinfog` behind the LU-with-MUMPS configuration; report
      unavailable for iterative configurations.
- [ ] **R5** — Controller: route `compute conditioning` to whichever backend the active solver
      supports, and keep the "n/a for this solver" message for the rest. Remove the temporary
      MUMPS-only guard added 2026-08-10.
- [ ] **R6** — `doc/input_file_reference.md`: state which solvers support the key and which
      quantity is reported (standing input-reference sync rule).

## 6. Open Design Questions

- **O1 — RESOLVED 2026-08-10 ( Christian ): expose both.** The condition number ( ADD ) stays the
  user-facing figure because it is the more descriptive one ; the forward error becomes a separate
  getter, useful to the timestepper. Every backend must supply both under the same names or
  declare them unavailable. Original question kept below for the reasoning.

  *(original)* **which quantity does `mConditionNumber` mean?** MUMPS COND1/COND2 (ADD), SuperLU
  `dgscon` rcond (norm-based kappa_1), SuperLU `dgsrfs` ferr/berr (ADD), and the retired spectral
  ratio are four different objects. The footer prints one number under one name; it must mean the
  same thing on every solver or it cannot be compared across runs. **Not decided.**
  My recommendation: the **forward error bound** (MUMPS RINFOG(9), SuperLU `ferr`), because it
  answers the stated purpose — how hard was this system to solve — and both libraries produce it
  natively.
- **O2 — cost.** `ICNTL(11) = 1` runs iterative refinement plus condition estimation on every
  armed solve. `ICNTL(11) = 2` gives backward error only, more cheaply. If O1 lands on
  ferr/RINFOG(9) the full setting is required; if it lands on backward error only, the cheap one
  suffices. **Not decided.**
- **O3 — universal fallback.** The scaled residual
  `\|\|Ax - b\|\|_inf / ( \|\|A\|\|_inf \|\|x\|\|_inf + \|\|b\|\|_inf )` costs ONE matvec, needs no
  factorization and no transpose, and works for every solver including STRUMPACK and PARDISO,
  which offer nothing else. Worth adding as the baseline even after R3/R4? **Not decided.**

## 7. Why not compute ADD ourselves

Verified 2026-08-10: there is **no transpose solve anywhere in the `Solver` API**
(`cl_Solver.hpp`, `cl_SolverWrapper.hpp`). The Hager-Higham estimator behind ADD needs
`\|\| |A^-1| |A| |x| \|\|`, which reduces to estimating `\|\|A^-1 D\|\|_inf` and therefore needs
products with both `A^-1` and `A^-T`. Adding that means a new entry point across every wrapper,
on top of the JOB 5/6 refactorization behaviour in the MUMPS wrapper
(`cl_SolverMUMPS.cpp:340-350`) that already blocked shift-and-invert. Using each library's own
implementation (R3, R4) avoids all of it.

## 8. Definition of Done

- [ ] One documented quantity, reported identically by every backend that supports it.
- [ ] `compute conditioning : true` either produces that number or says plainly it cannot.
- [ ] No solver silently reports a different quantity under the same name.
- [ ] `doc/input_file_reference.md` states the support matrix.
