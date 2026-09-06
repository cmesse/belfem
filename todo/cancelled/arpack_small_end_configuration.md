# ARPACK Small End: Automatic Configuration and Spectral Folding

**Date:** 2026-08-28
**Purpose:** Make the `compute conditioning` eigen fallback actually produce a number for the
thermal system, and stop it burning wall clock when it cannot. Mechanism: size `ncv` / `maxit`
from `( n, nev, which end )` instead of one hardwired floor, replace the mode-1 `'SM'` request with
a **spectral fold** ( `lambda_min = sigma - mu_max`, `mu` the spectrum of `sigma*I - A`,
`sigma = lambda_max` ), and latch the diagnostic off for the rest of the run once an end is proven
unreachable.
**Module:** `src/fem/kernel` (+ `src/sparse` Fortran drivers)
**AIs involved:** Claude (exploration + plan), Codex (audit), Grok (third voice)
**Status:** **IN PROGRESS — R1-R6 landed; 12 defects found and fixed ( D1-D12 ); the CENTRAL
QUESTION IS STILL UNANSWERED.** The second test run ( `out.txt`, 2026-08-28 21:48 ) settled which
end fails: the plain `'LM'` run for `rho` **converges**, the **fold** does not — 572 restarts, 6303
products, nothing. So the remaining question is not a bug but a method question: can a polynomial
( non-inverting ) Krylov transform reach the small end of a discretized diffusion operator at all?
My analysis says no — a shift is affine and preserves relative gaps, so the fold cures the FILTER
and not the SEPARATION — but the jury round did not engage it and **nobody has measured this
spectrum**. Next step is the executable gate ( §6.1 ), not more code.
*(Earlier: IN PROGRESS — R1-R6 landed, FIRST TEST RUN FAILED and its regression is fixed;
awaiting the second test run.)* The 2026-08-28 run
( `cmake-build-debug/tapestack3d/out4.txt` ) proved two things and refuted one: the **latch works**
( eigen-analysis cost goes to 0 ms after it fires, `out4.txt:302`, and the run continued, so the
collective structure holds under 4 ranks ), the **cost per failing step fell** 4747 -> 2519 ms, and
the **derived restart cap broke a converging end** ( D7 ). The fold itself is STILL UNTESTED — the
run never reached it, or if it did, the message could not say so ( D8 ). Plan audited twice and the code audited once by Codex; five plan defects and two code defects
adopted, one code defect ( D2 ) found by Claude before the audit reported it, one ( D1 ) found while
editing and not raised by any audit. Compiles clean: `gfortran -fsyntax-only -cpp -fopenmp` on both
drivers, `g++ -std=gnu++17 -fsyntax-only` with the build's own flags on the C++.
**Nothing is verified — no build has been run and no test has executed.**
*(Earlier: PLAN — revised 2026-08-28 after the Codex plan audit)* ( `tmp/ai_exchange/arpack_small_end.md` ).
Codex returned "revise before implementation" and was right on six counts; the fold survives but its
shift, its guards, its error table and its cost model all changed. Grok was **unavailable** for this
round — its CLI refuses to start because it cannot resolve `/run/podman/podman.sock` to build its
sandbox deny list, so this is a TWO-voice round, not three. No source modified.

> **Scope guards (from the task brief, Christian 2026-08-28):**
> - **Order is fixed:** R1-R6 land first, Christian test-runs `tapestack3d`, and only then does the
>   PETSc link ( `KSPComputeExtremeSingularValues` ) start. That link is **R7+, a separate step**,
>   not part of this batch.
> - **SLEPc stays deferred.** `todo/deferred/slepc_eigensolver_integration.md` is not revived here.
> - **No new `input.conf` keys.** The configuration is automatic, so the two-artifact input
>   contract is not triggered by R1-R6. R7 may need one; that is decided when R7 is planned.
> - Eigen*vectors* remain out of scope. `rvec = .false.` stays.
> - The serial matvec's missing scalar accumulator is **not** in this batch (see O3).

---

## 1. Measured Starting Point

From `cmake-build-debug/tapestack3d/out2.txt`, BDF5 step 308, 2026-08-28:

| quantity | value | source |
|---|---|---|
| thermal dofs `n` | 97095 | `out2.txt:74-76` |
| thermal solver | PETSc, `preconditioner : asm` -> GMRES | `input.conf`, `cl_SolverPETSC.cpp:625-634` |
| small-end outcome | `converged 0 of 1 in 301 restarts ( ncv = 20 )` | `out2.txt:203` |
| cost of that failure | **4747 ms per timestep**, every timestep | `out2.txt:216` |
| timestep wall clock | 16 s | `out2.txt:219` |
| magnetic kappa | 4.66e7 ( MUMPS COND1, not the spectral ratio ) | `out2.txt:217` |

So the diagnostic spends ~30 % of the timestep to print `n/a`. That is the second defect, and it is
independent of whether the small end is reachable at all.

## 2. Root Cause

`job = 0` maps to `which = 'SM'` with `bmat = 'I'` and `iparam( 7 ) = 1`
( `arpacktools.f90:69,149,158` ; `parpacktools.f90:103,269,278` ) — regular-mode Arnoldi on
`OP = A`. IRAM converges from the **exterior** of the spectrum of `OP`; `'SM'` asks for the region
where the filter polynomial has the least resolution, and on a diffusion operator that is also
where the eigenvalues are densest. `ncv = 20` leaves a filter of degree 19 to isolate one
eigenvalue out of that cluster. Confidence: high — this is the case the ARPACK users' guide sends
to shift-invert, and `todo/deferred/slepc_eigensolver_integration.md` §1 records the same mechanism.

**Not** the cause here: the accuracy floor `tol > eps_mach * kappa` documented at
`cl_FEM_DofMgr_EigenValues.hpp:75-88`. At `kappa ~ 1e7` and `tol = 1e-4` that inequality has three
orders of headroom. It bites later (see O1), not now. Confidence: high.

**Also true, and the reason "configure it better" is a fair description:** every tuning knob is
hardwired and unreachable. `set_subspace_size`, `set_tolerance`, `set_num_minvals`,
`set_num_maxvals` have **zero callers** outside the class, and `mNumMaxIter` has no setter at all
( `cl_FEM_DofMgr_EigenValues.hpp:60-90` ). `ncv = 20`, `maxit = 300`, `tol = 1e-4` apply to both
ends of the spectrum although the two ends differ by 20x in matvec count
( measured: `SM` 1231 numop vs `LM` 61, deferred SLEPc plan §1.1 ).

## 3. The Fold

Let `rho` be the magnitude the existing `'LM'` run already returns
( `cl_FEM_DofMgr_EigenValues.cpp:330-345` reduces to `sqrt( re^2 + im^2 )` ), and set

```
sigma = rho * ( 1 + gFoldMargin )        gFoldMargin = 1e-2          B = sigma*I - A
```

**The margin is load-bearing, not cosmetic** ( Codex re-audit §1 ). `rho` is a Ritz estimate and
may land *below* the true spectral radius. When it does, `sigma - lambda_max` goes negative and can
win the modulus reduction, returning `lambda_max` while every guard still passes:

```
spec( A ) = { 0.99995, 1.0 },  rho converged to 1e-4 lands at 0.9999
spec( B ) = { -5e-5, -1e-4 }   'LM' picks -1e-4  ->  sigma - mu = 1.0 = lambda_MAX, positive, accepted
with sigma = rho*1.01:
spec( B ) = { 9.949e-3, 9.899e-3 }  both positive  ->  sigma - mu = 0.99995 = lambda_min   correct
```

Verified numerically 2026-08-28. The formal requirement is `sigma >= ( lambda_min + lambda_max )/2`;
1 % clears it by two orders for any Ritz estimate worth accepting, and it is **free**, because
`sigma` cancels out of the answer entirely ( §3.3 ). A residual consistency check backs it up:
`lambda_min_alg <= rho` must hold, since `lambda_min <= lambda_max <= rho` — R4 rejects the fold if
it does not.

For a **real** spectrum every eigenvalue satisfies `lambda in [ -rho, rho ]`, so with `sigma > rho`
`spec( B ) = { sigma - lambda } subset ( 0, 2*sigma ]` is **strictly positive** — no assumption
about the sign of the spectrum is needed. The largest-magnitude eigenvalue of `B` is
therefore its algebraically largest, and

```
mu_max = sigma - lambda_min_alg          =>          lambda_min_alg = rho - mu_max
```

where `lambda_min_alg` is the **algebraically** smallest eigenvalue. Both runs are exterior `'LM'`
problems, which already converge in ~5 restarts. No factorization, no new TPL, no
reverse-communication solve — the JOB 5/6 refactorization blocker recorded in the deferred SLEPc
plan §1 is never touched.

### 3.1 Why the shift is the MAGNITUDE and not the signed lambda_max

> **CORRECTED 2026-08-28 (Codex, plan audit §1 — a real defect in the first draft, kept here
> because the negative result is what makes the current design safe).**

The first draft folded about the signed `lambda_max` and claimed that `sigma - mu_max <= 0` was a
self-detecting guard for a spectrum straddling zero. **That is false.** Codex's counterexample,
reproduced numerically 2026-08-28:

```
spec( A ) = { -10, 1 }
'LM' on A returns the largest MAGNITUDE, i.e. the signed value -10, so sigma = -10
spec( B ) = { 0, -11 }        'LM' on B returns mu = -11
sigma - mu = 1                positive, real, finite, converged -- every proposed guard passes
                              and the answer is lambda_MAX, not lambda_min
```

Folding about `rho = |−10| = 10` instead gives `spec( rho*I - A ) = { 20, 9 }`, both positive,
`mu_max = 20`, `lambda_min_alg = 10 - 20 = -10` — the true algebraically smallest value, and its
negative sign is a **correct rejection** rather than a silent wrong answer.

### 3.2 Admissibility contract

The sign of `lambda_min_alg` is now a certificate, not a guess:

| `lambda_min_alg = rho - mu_max` | meaning | action |
|---|---|---|
| `> 0` | every eigenvalue lies in `( 0, rho ]`, so the spectrum is **positive definite** and `lambda_max_alg = rho` | report `kappa = rho / lambda_min_alg` |
| `<= 0` | the spectrum straddles zero or touches it. The smallest-**magnitude** eigenvalue is then interior, and no fold can reach it | report unavailable, **once**, and latch |

That second row is the whole reason the D-SM fallback is dropped (O2, resolved below): the case
where the fold declines is exactly the case where `'SM'` returns a smallest *modulus* — a different
quantity from the `lambda_min_alg` the derivation needs — and it is also the interior problem that
does not converge in the first place. Falling back to it would buy a wrong number slowly.

Remaining guards, both runs:

| risk | detection | action |
|---|---|---|
| complex spectrum | `\|imag( mu )\| > tol * \|real( mu )\|` on the winning Ritz value | unavailable + latch |
| `rho` unconverged | the `'LM'` run returned NaN | unavailable, no fold attempted |
| converged but non-extremal Ritz value | **not detectable** — see O4 | accepted risk, documented |

The last row is Codex's §1 residue and is honest rather than solved: `NCONV >= nev` says ARPACK
converged what it was asked for, not that the value is the global extremum. Confidence that this
matters in practice for an exterior `'LM'` request: low. It is logged as O4, not papered over.

### 3.3 Accuracy — what the fold does and does not buy

The fold beats the **convergence rate** problem. It does not beat the **arithmetic floor**.

`sigma` cancels exactly. `B` is defined by whichever floating-point `sigma` we choose, so
`mu_max = sigma - lambda_min_alg` holds for that number, and `lambda_min_alg = sigma - mu_max`
subtracts the same stored value back out. Codex confirms this independently ( audit §2 ). The
consequences: `rho` may be computed at a loose tolerance, and the **same stored double** must be
passed to the matvec and to the final subtraction — not recomputed, not rounded in between.

What does not cancel is the cancellation itself. `mu_max ~ sigma`, so an error of
`max( tol, eps_mach ) * sigma` in `mu_max` lands undivided on `lambda_min_alg`:

```
rel_err( lambda_min_alg )  ~  max( tol, eps_mach ) * kappa
```

| kappa | `eps_mach * kappa` | `tol = 1e-10` term | 1 % on lambda_min? |
|---|---|---|---|
| 1e4 | 2e-12 | 1e-6 | yes |
| 1e7 | 2e-9 | 1e-3 | yes |
| 1e8 | 2e-8 | 1e-2 | marginal |
| 1e10 | 2e-6 | 1 | no at `tol = 1e-10`; the roundoff floor alone would still allow it |

> **CORRECTED 2026-08-28 (Codex, audit §2).** The first draft's table mis-stated the roundoff
> floor: `eps_mach * kappa` at `kappa = 1e13` is ~2e-3, i.e. 0.2 %, so the draft's "no" was wrong,
> and the roundoff-only 1 % threshold is near `kappa ~ 0.01 / eps ~ 4.5e13`. What actually binds at
> `tol = 1e-10` is the **requested tolerance**, not roundoff. Hence `tol = 1e-10` for the folded run.

Two further limits, both stated rather than fixed:

- **Non-normality.** The drivers are the nonsymmetric `dnaupd` / `pdnaupd` in mode 1. For a
  non-normal matrix a small Ritz residual can accompany a much larger eigenvalue error, scaled by
  the eigenvalue condition number. The model above is therefore a **heuristic estimate, not a
  bound** ( Codex audit §2, confidence medium ).
- **The reported ratio is not `kappa_2`.** `rho / lambda_min_alg` equals the 2-norm condition
  number only for a normal matrix. The admissibility contract proves the spectrum is real-positive,
  which is necessary but not sufficient. See O3.

## 4. Heuristic Sizing

Replace the single `mSubspaceSize = 20` floor with a function of `( end, n, nev )`. Both runs are
exterior now, so the two profiles differ only in tolerance:

```
lambda_max run ( 'LM' on A ):      ncv_want = max( 2*nev + 20, 20 ),  tol = 1e-4
folded run    ( 'LM' on sigma*I - A ): ncv_want = max( 2*nev + 20, 20 ),  tol = 1e-10

feasibility, in this order -- the floor is a HARD ARPACK requirement, so it is tested, never applied
as a max() after the caps ( Codex audit §6 ):

    ncv_min  = nev + 2                                    ARPACK: ncv - nev >= 2, hence ncv >= nev+2
                                                          ( NOT 2*nev+2 -- that is remark 4's
                                                            recommendation, not a requirement;
                                                            Codex re-audit §2 )
    ncv_cap  = min( n, gBasisBudget / ( 8 * n ) )         basis is n x ncv doubles. Both 8*n and
                                                          n*ncv are formed in a 64-bit type FIRST --
                                                          int_t may be int32, where 8*n overflows
                                                          above n ~ 2.7e8 ( Codex re-audit )
    if ( ncv_min > ncv_cap )  ->  diagnostic UNAVAILABLE, latched, one message
    ncv      = clamp( ncv_want, ncv_min, ncv_cap )

    maxit    = ceil( budget / max( 1, ncv - nev ) ),  budget = 2000 matvecs
```

> **CORRECTED 2026-08-28 (Codex, audit §5, independently found by Claude before the audit
> returned).** The first draft asserted "one restart costs `ncv - nev` matvecs" and sold the budget
> as a hard ceiling. Both are false. Measured ( deferred SLEPc plan §1.1 ): `SM` 1231 numop over
> 122 restarts = **10.1** per restart at `ncv - nev = 19`; `LM` 61 over 5 = **12.2**. ARPACK boosts
> `nev` internally between restarts so `np = ncv - nev` shrinks, and the initial Arnoldi
> factorization costs a further `ncv`. The two measured ratios differ, so no fixed multiplier
> exists. `maxit` is therefore a **restart cap whose worst case is bounded by the budget** — actual
> spend is typically half. A true operation budget would have to be enforced inside the reverse
> communication loop; that is out of scope here, and instead R4 **prints the achieved `NUMOP`**
> ( `info( gInfoNumOperations )`, already populated ) so the real cost is visible rather than modeled.

The explicit setters stay authoritative: a caller that has called `set_subspace_size` /
`set_tolerance` is not overridden by the heuristic ( R1 adds the "explicitly set" flags ).

## 4.0 THE MEASUREMENT (2026-08-28, reproducer — top of the evidence ladder)

`sysdump_thermal_3.hdf5`, dumped with `BELFEM_DUMP_SYSTEM=1` and analysed by
`tmp/analyse_thermal_spectrum.py`:

| quantity | value |
|---|---|
| n / nnz | 97095 / 1622509, CSR |
| asymmetry `\|\|A - A^T\|\|_max / \|\|A\|\|_max` | **3.549e-16 — the matrix is SYMMETRIC** |
| `lambda_max` | 1.830271 |
| six smallest | 6.1760e-8, 6.3500e-8, 6.3606e-8, 6.3799e-8, 6.4471e-8, 6.4749e-8 |
| `lambda_min` | 6.175959e-08 > 0 — **symmetric positive definite** |
| **kappa_2** | **2.9635e+07** |
| relative gap at the small end | **9.510e-10** |
| the same gap under shift-invert | 2.818e-02 |
| amplification | 2.964e+07 |

**F1 of the pre-registration is CONFIRMED BY MEASUREMENT, not by agreement.** The gap the fold has
to work with is 9.5e-10, five to six orders below the ~1e-4 at which a polynomial Krylov method
stops being able to isolate an eigenvalue. No `ncv`, no `maxit`, no tolerance changes that: the six
smallest eigenvalues agree to three significant figures, which is the dense smooth-mode cluster a
diffusion operator has by construction. **The fold is structurally dead for this operator class**,
and the 572 restarts / 6303 products it spent were the arithmetic being honest.

Three further consequences, each independent of the fold:

1. **The premise was off by two and a half orders.** kappa_2 is 2.96e7, not the ~1e5 the task
   assumed. The ~1e4-1e5 figure came from MUMPS COND1, an Arioli/Demmel/Duff *linear-system*
   condition number — a different object, as this plan has said throughout. The measured value is
   also close to the magnetic system's MUMPS number ( 4.7e7 - 2.7e8 across the run ), so the thermal
   system is not the well-conditioned one it was believed to be.
2. **BELFEM is using the wrong ARPACK driver family.** The matrix is symmetric to machine precision,
   and `src/sparse` contains **no** `dsaupd` / `dseupd` / `pdsaupd` — only the nonsymmetric
   `dnaupd` path. Lanczos would be cheaper ( short recurrence, far less re-orthogonalization ) and
   carries the convergence theory the nonsymmetric driver does not. It would NOT fix the gap.
   Filed as a separate item; do not fold it into this plan.
3. **The fold's ACCURACY design was sound; only its convergence was not.** At kappa = 2.96e7 the
   cancellation model gives `rel_err ~ tol*kappa = 1e-10 * 3e7 = 3e-3`, which is a perfectly good
   diagnostic. Had it converged, the number would have been usable. The tight tolerance was not the
   problem.

## 4.1 Defects Found While Implementing

- [x] **D1 (HIGH) — the serial ARPACK matvec has an OpenMP data race.** `a` and `b` are declared at
      subroutine scope ( `arpacktools.f90:96-97` before the fix ) and are absent from the
      `private` clause of the `!$omp parallel do` that used them, so they defaulted to **shared**:
      every thread in the team wrote the same pair of row bounds, and a row could be summed with
      another row's limits. The distributed driver never had the bug — it inlines `pointers(i)`
      into the loop header ( `parpacktools.f90:327-334` ). `USE_OPENMP` is ON by default, so any
      single-rank run with more than one thread had a racy matvec; the 2026-08-10 serial
      measurements in the deferred SLEPc plan §1.1 were taken at 1 and 2 threads and are suspect
      at 2. **Fixed 2026-08-28** as part of R3: the row bounds moved into the loop header and the
      row sum into a private `acc`, matching the distributed driver. Found by Claude while editing;
      not reported by either audit round, because the plan did not ask about the existing loop.
- [x] **D2 (HIGH) — an infeasible subspace was indistinguishable from a failed iterate.** Both
      drivers returned a bare `BELFEM_QUIET_NAN` when `configure()` found no legal `ncv`, and
      `compute_conditioning()` mapped every NaN to `NotConverged` — so `EigenOutcome::Infeasible`
      could never be produced, and a condition that can never improve spent three strikes before
      latching. **Fixed 2026-08-28** with an `mSubspaceInfeasible` member, set by whichever driver
      declined and cleared at the top of every `run()`. Found by Claude while self-reviewing, and
      independently by the Codex code audit ( which read the tree before the fix landed ) — the
      agreement is worth recording, but it is two readings of the same static text, not evidence
      from a run.
- [x] **D3 (MEDIUM) — an explicitly set `ncv` bypassed every legality check.** `configure()` tested
      feasibility against the automatic floor and cap and then copied `mSubspaceSize` over the
      result unvalidated, so an explicit value could sit below ARPACK's legal minimum ( where the
      Fortran driver would silently raise it, meaning it was never authoritative ) or above the
      basis budget ( where it would allocate past it ). **Fixed 2026-08-28**: the explicit value is
      clamped into the legal range and the caller is told when it had to move. Found by the Codex
      code audit.
- [x] **D4 (MEDIUM) — the feasibility floor was below what the drivers actually allocate.** The
      first revision set the floor at ARPACK's hard bound `nev + 2`, on Codex's correct objection
      that remark 4's `2*nev + 1` is a recommendation. But **both drivers apply that recommendation
      themselves** — `ncv = max( 2*nev + 1, ncvmin )` at `arpacktools.f90:154` and
      `parpacktools.f90:249` — so a subspace checked against `nev + 2` would be approved and then
      allocated larger than the budget it was checked against. **Fixed 2026-08-28**: the floor is
      `max( nev + 2, 2*nev + 1 )`, i.e. what will actually be allocated. Identical for `nev = 1`,
      the only value in use. Found by Claude while applying D3.
- [x] **D5 (LOW) — the derived restart cap silently weakened a public entry point.** The product
      budget of 2000 is calibrated for the exterior requests the diagnostic now makes, and would
      have cut `compute_smallest_eigenvalues()` — an INTERIOR `'SM'` request, still public, though
      no longer reached from the conditioning path — from the old fixed 300 restarts to 105.
      **Fixed 2026-08-28**: the budget is job-aware, 20000 for job 0 ( giving ~1052 restarts ) and
      2000 for the exterior ends. Found by Claude while checking behaviour preservation.
- [x] **D7 (HIGH, REGRESSION, caught by Christian's test run 2026-08-28) — the derived restart cap
      broke an end that had been converging.** First run of the new code: `PARPACK converged 0 of 1
      in 96 restarts ( ncv = 22 )`, thermal kappa still `n/a`
      ( `cmake-build-debug/tapestack3d/out4.txt:209` ). 96 restarts is exactly the derived cap,
      `2000 / ( 22 - 1 ) = 95`. The old code ran `maxit = 300`, and the old output carried exactly
      ONE failure block per timestep although the old `compute_conditioning()` ran BOTH ends and
      printed on every failure.
      > **QUALIFIED 2026-08-28 (Codex follow-up audit).** I wrote that this **established** that the
      > `'LM'` end had been converging within 300 restarts. It does not. One failure block proves
      > that exactly one end failed, **not which one** — Codex confirmed there is no path where a
      > run fails to converge silently, but the old messages were unlabeled, so `'SM'` failed and
      > `'LM'` succeeded ( overwhelmingly likely, and what the 2026-08-10 measurements predict ) and
      > the reverse are both admissible. The 2519 ms cannot arbitrate either: the 0.83 ms/product
      > figure is from a different matrix and excludes distribution and restart overhead.
      > Three readings remain open, and they do not have the same remedy:
      > **(A)** the plain run failed at 95 where it used to converge at 300 — the cap is the cause
      > and the fix works; **(B)** the plain run succeeded and the FOLD failed at 95 — the cap was
      > too small for the tighter tolerance, and the fix still works; **(C)** the plain run fails at
      > 300 too, i.e. `'LM'` was the end failing all along — the cap is innocent and this fix does
      > **nothing**. D8's labels are what discriminate them, which is why both fixes shipped
      > together.
      **Root cause: the product budget was calibrated on the 2026-08-10 measurement at n = 10428 and
      applied unguarded to a matrix nine times larger.** The general lesson is not about the number:
      a heuristic that can WEAKEN a working configuration is a regression with a formula in front of
      it. The same trap was avoided by hand for job 0 ( D5 ) and then walked into on job 1, because
      D5 was treated as a special case instead of a rule.
      **Fixed 2026-08-28:** `mNumMaxIter` ( 300 ) is now a FLOOR that the budget may only raise,
      never replace, so no automatic sizing can land below the historical default. Budgets are
      keyed to how hard the request is — 4000 plain, 12000 folded ( its tolerance is 1e-10 against
      1e-4, six more decades, which costs restarts ), 20000 for the interior job 0. `maxit` is a
      cap, not a target, so headroom is only ever paid for by a run that was going to fail anyway,
      and the latch bounds how often that happens.
      **Cost qualification ( Codex ):** the new per-attempt worst case can exceed BOTH the previous
      change and the original fixed-300 scheme, because a plain-success/fold-failure sequence pays
      300 restarts plus 571 on each of three strikes. The product budgets are heuristics, not
      enforced ceilings — nothing stops a run mid-flight. That is a deliberate
      convergence-versus-cost trade, and it is the reason the achieved product counts are now
      printed rather than modeled.
- [x] **D8 (MEDIUM, same run) — a failure did not say WHICH run failed.** Both messages printed
      `ncv` and the restart count but not which of the two ends was in flight, so `out4.txt` cannot
      distinguish "the plain `'LM'` run died" from "the plain run succeeded and the fold died" —
      the two have completely different remedies. R4 had also promised to print the achieved
      `NUMOP` and did not; the code audit did not catch it because the audit brief did not ask.
      **Fixed 2026-08-28:** an `mRunLabel` member names the run in every message, both failure paths
      print the achieved product count and the tolerance in force, and the success path prints
      `lambda_max`, `lambda_min` and both product counts — the numbers anyone re-tuning the budget
      actually needs.
- [x] **D9 (P1, jury round 2026-08-28, Codex) — a transient indefinite iterate switched the
      diagnostic off permanently.** `NotPositiveDefinite` latched on its FIRST occurrence, justified
      by a comment of mine asserting it is "a property of the problem, not of this timestep". That
      assertion is false: the Jacobian is rebuilt every timestep, the deck is nonlinear and Δt
      varies, so the spectrum moves with it, and nothing says an indefinite iterate at step k
      implies one at step k+1. One transient sign flip would have suppressed every later estimate.
      **Fixed 2026-08-28:** only `Infeasible` latches immediately — it depends on the problem SIZE
      and the basis budget, which cannot come right next step. Everything else, `NotPositiveDefinite`
      included, takes the strike count, so a spectrum that really is indefinite by construction
      still latches, it just has to prove it over `mMaxFailures` steps instead of being assumed.
- [x] **D10 (P2, same round, Codex) — the resolution-limit report ignored an explicit tolerance.**
      `configure()` honours `set_tolerance()` and returns `mEpsilon`, but the resolution test and its
      message both used `mFoldEpsilon` unconditionally, so an explicit caller got the wrong threshold
      and a printed tolerance that was not the one in force. **Fixed 2026-08-28:** `configure()`
      records what it chose in `mEffectiveTolerance` and the report reads that, so the rule lives in
      one place. `configure()` is no longer `const` — it has a side effect now, and claiming
      otherwise would be a lie in the signature.
- [x] **D11 (P2) — the success message was not rank-guarded.** Found by Claude while fixing D10:
      `message()` does not filter by rank ( every other message in the file guards it explicitly ),
      so the new success line would have printed once per rank — four times on this deck — and on a
      non-master rank it would have read `mEffectiveTolerance` as a NaN, since `configure()` runs on
      the master alone under PARPACK. **Fixed 2026-08-28** with an explicit `mCommRank == 0` guard
      around the reporting only, leaving the returned value computed on every rank.
- [x] **D12 (P2, pre-registration F8) — a failing fold said nothing about the end that worked.**
      `rho` was reported only on FULL success, which the 2026-08-28 runs never reached, so no log
      could say what `lambda_max` was or how cheap the first end is. **Fixed 2026-08-28:** `rho` and
      its product count are printed as soon as that run converges, independently of what the fold
      then does.
- [x] **D6 (LOW) — `a`-prefixed locals in new code.** `aResult` in `reduce_extremum()` and `aKappa`
      in `compute_conditioning()`; `a` is the ARGUMENT prefix. **Fixed 2026-08-28** to `tResult` and
      `tKappa`. The pre-existing `aResult` locals in `run()` / `run_arpack()` / `run_parpack()` are
      left alone — the minimal-rename rule applies to code the change is not otherwise touching.
      Found by the Codex code audit.

## 5. Ordered Steps

- [x] **R1** — `EigenValues`: `mSubspaceSizeExplicit` / `mToleranceExplicit` flags set by the
      existing setters, an `mNumMaxIter` setter, and a private `configure( ... )` implementing §4
      including the **feasibility test** ( not a trailing `max()` ) and 64-bit budget arithmetic.
- [x] **R2** — `run_arpack` / `run_parpack`: keep the extremal **magnitude** as the return value,
      and additionally record the signed real and imaginary parts **copied from the exact slot that
      won the modulus reduction** ( Codex audit §6 — the winner's identity, not just its value ).
      (after: R1)
- [x] **R3** — Fortran drivers: `sigma` argument on `arpack_standard_eigen` and
      `parpack_standard_eigen`; the reverse-communication slot computes `y = sigma*x - A*x` when it
      is non-zero. In `parpacktools.f90` this is `y(i) = sigma*x(i) - acc`, where **`x` is the local
      `nloc` slice of `workd`** ( `parpacktools.f90:307` ), *not* `xglobal` — `xglobal(i)` would be
      wrong on every rank above zero, and the global-index form would need
      `xglobal( tDispls(rank+1) + i )`. In `arpacktools.f90` the same line also introduces the
      scalar accumulator the serial driver still lacks, because the shifted form requires it
      anyway. Update both `extern "C"` prototypes. (after: R2)
- [x] **R4** — `compute_conditioning()`: `rho` first, then the fold at `sigma = rho`, then the §3.2
      admissibility test. **The master decides an outcome code — `Ok`, `NotPositiveDefinite`,
      `Complex`, `NotConverged`, `Infeasible` — and broadcasts it before any branch that may call
      `run()`** ( Codex audit §4 ). Print the achieved `NUMOP`. (after: R3)
- [x] **R5** — Latch: `mDiagnosticUnavailable`, set **only from the broadcast outcome code**, never
      from a rank-local floating-point comparison. Not cleared by `reset()` ( which today clears
      only `mMatrixFlag`, `cl_FEM_DofMgr_EigenValues.cpp:789` ), so it persists across Jacobian
      rebuilds — which is the intent. This is what stops the 4747 ms/step bleed. (after: R4)
- [ ] ~~**R5b** — direct-`SM` fallback ladder.~~ **STRUCK 2026-08-28 (O2 resolved).** The fallback
      would only run when the fold declines, i.e. on a spectrum straddling zero, where `'SM'`
      returns the smallest *modulus* — not the `lambda_min_alg` the derivation needs ( Codex audit
      §6 ) — and where the interior problem does not converge anyway. Dropping it also removes the
      17 s worst case.
- [x] **R6** — Devlog + `todo/README.md`; reconcile `todo/conditioning_diagnostic_backends.md` §1,
      which says a spectral transform requires the factorization-reuse fix. The fold is a spectral
      transform that does not. (after: R5)
- [ ] **R7** — **SEPARATE STEP, after Christian's test run.** PETSc `KSPSetComputeSingularValues` +
      `KSPComputeExtremeSingularValues` behind a capability flag ( GMRES-class KSP only ), falling
      back to the eigen path when unavailable. The deferred SLEPc plan dismissed this call, but only
      for the `KSPPREONLY` + LU configuration; the thermal deck runs GMRES + ASM, where it is
      available. It reports `sigma_max/sigma_min` of the **preconditioned** operator — a different
      quantity, which R7's planning must name explicitly.

## 6. Open Design Questions

- **O1 — what should the footer show when the fold hits the resolution floor?** The cell is a fixed
  six characters ( `conditioning_string`, `cl_FEM_Controller.cpp:604-640` ), so a `>` prefix does not
  fit without shifting the column. R4 proposes: return the value, print the resolution caveat as a
  separate message line, leave the footer format alone. Note Codex ( audit §6 ) is right that this
  is a **resolution-limited estimate, not a rigorous lower bound** — the message must not say
  "kappa >= X" as though it were proved. **Not decided.**
- **O2 — RESOLVED 2026-08-28 ( Claude, on Codex audit §6 ): drop the direct-`SM` fallback.** See the
  struck R5b.
- **O3 — the footer calls the number "Conditioning" for two different quantities.** The MUMPS branch
  reports ADD COND1; the eigen branch reports a spectral ratio that equals `kappa_2` only for a
  normal matrix. This predates the plan ( deferred SLEPc plan §1, last row ) and the fold does not
  fix it. **Not decided** — but it is now the largest remaining honesty gap in the diagnostic.
- **O4 — a converged Ritz value is not certified extremal.** `NCONV >= nev` does not prove the
  global extremum was found ( Codex audit §1 ). No cheap certificate exists inside ARPACK's
  interface. Accepted risk for an exterior request; logged so it is not rediscovered as a defect.
- **O5 — the serial matvec's missing scalar accumulator** ( `arpacktools.f90:196-204` versus the
  `acc` form at `parpacktools.f90:365-372`, measured 2.2x at `-O2` ) is now folded into R3, because
  the shifted expression needs the accumulator anyway. Previously O3 in the first draft.

## 4.2 Found by the Jury Round but NOT This Plan's Work

`scripts/cross_review.sh` reviews the whole working-tree diff, which is shared with a concurrent
session. Two confirmed defects in that session's files came back with the round and are recorded
here only so they are not lost — **they belong to whoever owns that work**:

- `hdf5::get_groups()` ( `src/io/hdf5_tools.hpp` ) returned EVERY link at the location, while the
  new `db2exo` ( `src/physics/database/db2exo.cpp`, untracked ) calls `select_group()` on each name
  it returns, so any dataset or named datatype at that location would fail the group open rather
  than being skipped. Fixed 2026-08-28 at Christian's instruction: the callback now checks
  `H5Oget_info_by_name3( ... H5O_INFO_BASIC ... )` and appends only `H5O_TYPE_GROUP`, and skips
  non-hard links rather than following them out of the location.
- Trailing whitespace at `src/io/hdf5_tools.hpp:78` failed `git diff --check HEAD`. Cleared by the
  same edit; the gate is clean.
- **Second-order defect in that fix, raised by the fix audit and also fixed:** the first version
  turned a failed `H5Oget_info_by_name3` into "skip this entry", so a permission failure or a
  corrupt object would have produced a short group list indistinguishable from a complete one, and
  `get_groups()` ignored `H5Literate2`'s return value entirely. The callback now returns `-1` on a
  query failure ( the contract: negative aborts the iteration and becomes the iterate's return
  value ) while a non-group still returns 0, and `get_groups()` checks the result with an
  always-active `BELFEM_ERROR` — it runs once per file, so the tier is affordable and the message
  survives a release build.

## 6.1 What Is NOT Done

**The measurement, which is now the only thing that moves this forward.** `tmp/analyse_thermal_spectrum.py`
is written, compiles, and its dependencies are present ( h5py 3.16, scipy 1.17, numpy 2.4 ). It needs
one `SpMatrix::save()` call in `compute_conditioning()` and one run, and it reports: whether the
thermal Jacobian is symmetric at all ( which decides both whether the spectral ratio IS kappa_2 and
whether `dsaupd` should replace `dnaupd` ), the true `lambda_min` by sparse-LU shift-invert, the
actual kappa, and the small-end relative gap beside what shift-invert would amplify it to. That
single number decides between "tune further", "shift-invert", and "use a different quantity".


- **Nothing has been executed.** Syntax checks are not a build and a build is not a run. The
  conditioning number, the fold's accuracy, the latch, and the OpenMP race fix are all unverified
  ( protocol §11 evidence ladder: this is "reviewed", not "verified" ).
- Grok never saw any of it. Its CLI would not start on this machine
  ( `/run/podman/podman.sock`, permission denied while building the sandbox deny list ), so the
  third voice is missing from the plan round AND the code round. A standalone briefing was written
  for manual submission: `tmp/ai_exchange/grok_briefing_arpack_small_end.md`.
- R6 ( devlog, README, reconciling `conditioning_diagnostic_backends.md` ) and R7 ( the PETSc
  singular-value link ) are open by design — R7 waits on the test run.

## 7. Definition of Done

- [ ] `tapestack3d` reports a thermal conditioning number, or a stated lower bound, or says once
      that it cannot — and in the last case costs nothing on subsequent steps.
- [ ] The fold shifts by `rho * ( 1 + gFoldMargin )`, never by a signed `lambda_max` and never by a
      bare `rho` (§3, §3.1).
- [ ] `lambda_min_alg <= rho` is checked before the ratio is formed.
- [ ] Every branch that may call `run()` is reached from a broadcast outcome code.
- [ ] Collective structure unchanged: every rank still enters and leaves `run_parpack` together.
- [ ] `make check-fast` clean.
