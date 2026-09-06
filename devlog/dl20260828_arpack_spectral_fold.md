# ARPACK Small End: Spectral Fold, Automatic Sizing, and a Latch

**Date:** 2026-08-28
**Purpose:** Record the change that replaces the non-converging `'SM'` request in the conditioning
diagnostic with a spectral fold, sizes the ARPACK knobs automatically, and stops the diagnostic
retrying an unreachable spectrum every timestep.
**Module:** `src/fem/kernel`, `src/sparse`
**AIs involved:** Claude (exploration, plan, code), Codex (three audit rounds), Grok (**unavailable**)

## What prompted it

Christian's `tapestack3d` run reported `n/a` for the thermal conditioning number on every timestep
and spent 4747 ms per step doing so, out of a ~16 s step (`out2.txt:203-219`, n = 97095, PETSc+ASM,
4 ranks). The reasonable question was whether ARPACK was simply configured badly, given MUMPS put
the magnetic system at 4.66e7.

Two things had to be separated first. MUMPS reports an Arioli/Demmel/Duff *linear-system* condition
number (RINFOG(10)); the ARPACK path reports a *spectral ratio*. They are different quantities and
neither predicts the other. And the accuracy-floor argument recorded at
`cl_FEM_DofMgr_EigenValues.hpp` (`tol > eps_mach * kappa`) was not what was biting — at that kappa
it has orders of headroom. The actual cause is structural: `job = 0` maps to `which = 'SM'` with
`iparam(7) = 1`, regular-mode Arnoldi on `OP = A`, which converges from the exterior of the spectrum
and has essentially no resolution at the interior, where a diffusion operator is also densest.

## What landed

`sigma = rho * ( 1 + 1e-2 )` with `rho = max|lambda|` from the existing `'LM'` run; the
reverse-communication slot then applies `OP = sigma*I - A`, whose spectrum is strictly positive for
a real spectrum, so its largest eigenvalue is `sigma - lambda_min` — an exterior request again.
`lambda_min = sigma - mu`, and `sigma` cancels exactly because the same stored double goes in and
comes back out. Both ends are now problems Arnoldi is good at.

Around it: `ncv`, `maxit` and `tol` are sized from `( n, nev, end )` instead of three hardwired
members that no caller in the tree ever reached; the master decides an `EigenOutcome` and broadcasts
it before any branch that may re-enter a collective; and an end proven unreachable latches off for
the rest of the run instead of costing 4.7 s per step forever.

Full reasoning, the corrected error analysis, and the audit trail:
`todo/arpack_small_end_configuration.md`.

## What the audits changed

Codex ran three rounds (plan, plan re-check, code) and was right every time it objected. Two
findings changed the design rather than the wording:

1. **The first draft folded about the signed `lambda_max`.** `'LM'` returns the largest *magnitude*,
   so on `spec(A) = {-10, 1}` the shift is `-10`, the fold returns `1`, and every guard passes while
   the answer is `lambda_max`. Folding about the magnitude instead makes `spec(sigma*I - A)`
   non-negative unconditionally, and the sign of the result becomes a genuine positive-definiteness
   certificate. Reproduced numerically before accepting it.
2. **An underestimated `rho` breaks it too**, since `sigma - lambda_max` can go negative and win the
   modulus reduction. Hence the 1 % margin, which is free because `sigma` cancels.

It also killed the direct-`'SM'` fallback (it answers a different question — smallest *modulus*,
not smallest algebraic value — in exactly the case the fold declines), corrected the cost model
(`maxit` is a restart cap, not a matvec ceiling: measured 10.1 and 12.2 products per restart at
`ncv - nev = 19`, because ARPACK boosts `nev` internally), and caught two implementation defects.

## Found while editing, by neither audit

**The serial matvec had an OpenMP data race.** `a` and `b` were declared at subroutine scope in
`arpacktools.f90` and were absent from the `private` clause of the loop that used them, so they
defaulted to shared and every thread wrote the same row bounds. `USE_OPENMP` is ON by default, so
any single-rank run with more than one thread had a racy matvec. The distributed driver never had
it — it inlines `pointers(i)` into the loop header. Note that the 2026-08-10 serial timings recorded
in `todo/deferred/slepc_eigensolver_integration.md` §1.1 were taken at 1 and 2 threads, and the
2-thread numbers are suspect.

## The measurement that ended it

`BELFEM_DUMP_SYSTEM=1` ( Christian pointed out `save_system` already existed — the dump only needed
routing at the thermal field, since both call sites dumped the magnetic one ) plus
`tmp/analyse_thermal_spectrum.py`:

```
n = 97095, nnz = 1622509
asymmetry ||A-A^T||/||A|| = 3.5e-16          -> SYMMETRIC
lambda_min = 6.175959e-08, lambda_max = 1.830271   -> SPD
kappa_2 = 2.9635e7
small-end relative gap        = 9.51e-10     <- what the fold has
same gap under shift-invert   = 2.82e-02     <- amplification 3.0e7
```

The six smallest eigenvalues agree to three significant figures. That is the dense smooth-mode
cluster a diffusion operator has by construction, and 9.5e-10 is five orders below where a
polynomial Krylov method can still isolate an eigenvalue. **The fold is structurally dead for this
operator class** — confirmed by reproducer, not by agreement. Two side findings: kappa_2 is 2.96e7,
not the ~1e5 the task assumed ( that number was MUMPS COND1, a different object ), and the matrix is
symmetric while BELFEM had no symmetric ARPACK driver at all.

Christian's rulings that followed: report the solver-native metric where available and fall back to
shift-invert, RAM is the user's responsibility, sample once per run rather than per step, and
implement the missing symmetric drivers. The first three are planned in
`todo/conditioning_shift_invert_fallback.md`; the fourth landed the same day —
`arpack_symmetric_eigen` / `parpack_symmetric_eigen` with their own `check_saupd` / `check_seupd`,
since the symmetric error tables are not the nonsymmetric ones.

## Sequel, 2026-08-29: shift-invert lands and the gate passes

The fold was retired and replaced by ARPACK mode 3 with `sigma = 0`, driven from C++ with the
master running ARPACK and every rank servicing each inverse application through a collective MUMPS
solve. A scoped freeze on the wrapper turns those into `JOB 3` solve-only calls against one
factorization.

Verified by execution against a 1D Laplacian with an analytic spectrum — picked because its small
end is dense ( relative gap 1.8e-6 ), the same pathology that defeated the fold:

```
lambda_max  rel 2.11e-15   lambda_min  rel 6.25e-12   kappa_2  rel 6.25e-12   GATE PASSED
MUMPS: 1 x JOB 6, 60 x JOB 3      small end converged in ONE restart
```

`make check` green on top of that ( Christian ), which is the wiring gate the probe cannot give.

**The finding worth carrying forward: running it caught a CRITICAL defect that three audit rounds
could not.** The plan audit recommended `SymmetryMode::GeneralSymmetric` for the dedicated solver,
since the matrix is symmetric and MUMPS could halve the factor memory. But BELFEM stores the FULL
matrix and passes `SymmetryMode` straight into MUMPS `SYM`, and nothing anywhere extracts a
triangle — so `SYM = 2` double-counts every off-diagonal and factorizes a different matrix. It does
not fail. It returned `lambda_min = -2.5e-19` against a true `2.46e-6`, with ARPACK reporting
`info 0` and `nconv 1`. First-solve residual was 7.4e16; with `SYM = 0` it is 1.3e-12. No static
reading finds that, because it is a property of how BELFEM feeds MUMPS rather than of the code on
the page.

Two traps recorded for the next caller: the ARPACK drivers require a Fortran-based matrix and the
base must be restored before MUMPS sees it ( otherwise the matvec reads one element before
`values` and the abort surfaces inside an unrelated `malloc` ); and changing a class layout
invalidates every TU that allocates it, so a scratchpad probe must recompile `cl_Solver.cpp` too,
not just the file it edited.

Detail: `todo/conditioning_shift_invert_fallback.md` §3.5.

## Session close, 2026-08-29

Landed after the gate passed: **R5/R8** — `Controller::finalize_run()` runs the eigen fallback ONCE
after the time loops, wired into all THREE executables that call `finalize()`, with the quantity
NAMED beside each number ( "MUMPS COND1, rhs-dependent" versus "kappa_2, spectral" ) because those
two differ by ~40x on the same operator and would otherwise read as a bug. Arming changed to
"held until a capture actually happens" rather than first-step-only, on the DR-127 session's
finding that a timestep can now legitimately take zero solves.

Four debt-register outcomes, two of them mine:

- **DR-140** ( filed, open ) — the MUMPS symmetric-mode contract. Fixed across two sessions:
  wrapper guard plus three corrected defaults here, `cl_IWG.hpp` base default flipped there on
  Christian's ruling. Left OPEN deliberately: the defect is fixed but the CAPABILITY is not, and a
  struck row would tell a future reader symmetric modes are handled when the truth is they are
  refused.
- **DR-142** ( filed and struck same day ) — `EigenValues` owned four raw pointers and deleted none
  of its copy/move operations. Self-reported: one of those four is a pointer this session ADDED, so
  an existing hazard was extended without noticing, and four audit rounds missed it because it lay
  outside the diff's stated scope. The audits did their job; the scope was set too narrowly.
- DR-141 ( `Solver`, same shape ) stayed with the session that found it — split rather than folded,
  because one row owning edits in two sessions' files is how gates rot.

**The lesson worth carrying: running it found what reading it could not.** Three static audit rounds
read the symmetry code and none could see the defect, because it is a property of how BELFEM feeds
MUMPS rather than of the code on the page. The probe found it in one run. Corollary, learned twice
the hard way the same day: a hand-linked probe fails OPEN — a stale object makes a guard ABSENT, and
absent is indistinguishable from broken from outside. Diff object mtimes against source mtimes, and
pair a guard test with its inverse so a stale binary cannot show green in both directions.

Also retired one recurring failure mode: `todo/debt_register.md`'s header restates a count its own
rows carry, and nothing recomputed it. Three sessions edited that sentence today and it went stale
twice. `scripts/check_doc_claims.py` now recounts the open `[P]` and `[W]` rows and diffs them
against the header — and caught the author of the rule getting it wrong within the hour.

## Status

**Reviewed, not verified.** Both Fortran drivers pass `gfortran -fsyntax-only -cpp -fopenmp` and the
C++ passes `g++ -std=gnu++17 -fsyntax-only` under the build's own flags. Nothing has been built and
nothing has been run. Christian's test run is the next step, and R7 (the PETSc
`KSPComputeExtremeSingularValues` link, available here because the thermal deck runs GMRES+ASM
rather than the `KSPPREONLY`+LU configuration the deferred SLEPc plan assumed) waits on its result.

Grok was absent from every round: its CLI refuses to start on this machine because it cannot resolve
`/run/podman/podman.sock` to build its sandbox deny list. A standalone briefing for manual
submission is at `tmp/ai_exchange/grok_briefing_arpack_small_end.md`.
