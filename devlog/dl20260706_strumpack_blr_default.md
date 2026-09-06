# STRUMPACK BLR Compression Off by Default (Rank-Dependent Newton Stall Root-Caused)

**Date:** 2026-07-06
**Purpose:** Root-cause and fix erratic, rank-count-dependent nonlinear convergence with STRUMPACK; BLR compression now off unless explicitly requested.
**Module:** sparse

## Symptom

Same coupled h-φ timestep (hphirun, N = 72,285) behaved wildly differently
per MPI rank count: 2 procs converged in a few iterations, 1 proc exploded
mid-Newton (−61 dB → −23.5 dB at Newton 9) and ground at 3e-3, 4 procs never
got past a ~4e-4 residual floor with the relaxation controller thrashing
between Newton and Picard. MUMPS on the same case was fine.

## Root cause (verified in STRUMPACK 7.2.0 sources)

- BELFEM's old N-based heuristic (pre-dating this week) enabled BLR
  compression for N ≥ 50k — for serial AND MPI alike.
- STRUMPACK's default BLR relative tolerance is **1e-4**
  (`BLROptions.hpp:46-48`); BELFEM set only the leaf size, never the
  tolerance → the factorization is a lossy preconditioner good to ~4 digits.
- With compression on, `KrylovSolver::AUTO` solves via GMRES(30)
  preconditioned by that lossy factorization
  (`SparseSolverMPIDist.cpp:348-352`) — and `solve_internal` returns
  `SUCCESS` **unconditionally** (`:389`), so GMRES stagnation near the
  1e-4 preconditioner quality is silent. BELFEM's return-code check
  cannot see it.
- Rank dependence: large fronts use distributed BLR (`BLRMatrixMPI`) whose
  tile grid follows the process count → preconditioner quality varies with
  ranks (worst at 4 here: first-solve floor −33.8 dB vs −61.8 dB at 1/2).
- The "2 procs better than 1 proc" paradox was controller path divergence,
  not a solver ranking: solver noise seeds different iterates; 2p hit its
  jump at relax=1.0 and recovered via Picard at 0.5, 1p hit it with relax
  already ground to 0.24 and its fallback ran at 0.12.

**Confirmed by experiment:** `compression scheme : off ;` restored clean,
rank-independent convergence at 1/2/4 procs (Christian, 2026-07-06). Note
this input key only works since the S2 fix (dl20260705) — previously an
explicit `off` was silently overridden back to BLR.

## Change

`strumpacktools.cpp` — `set_strumpack_options` compression switch:

- `AUTOMATIC` (the default) and `OFF` now both disable compression for
  STRUMPACK. The old N ≥ 50k auto-BLR heuristic is gone.
- Explicit `compression scheme : blr` still works and now gets, in addition
  to the N-based leaf size, `BLR_options().set_rel_tol( compression_cutoff() )`
  = 1e-8 (first actual use of `SolverParameters::mCompressionCutoff`) so an
  explicit request doesn't fall into the same 1e-4 trap.
- MUMPS untouched: its `Automatic` BLR mode is conservative and handled the
  case well.

## Verification

- `-fsyntax-only` clean with real flags.make flags.
- Christian's 1/2/4-proc runs with compression off: clean convergence
  (that experiment is exactly what the new default reproduces).

## Follow-up same day: small-system turtle retired

With compression off, 8 procs at 9k DOFs/proc runs beautifully — directly
refuting the `hatch_small_turtle` warning box ("too small for distributed
STRUMPACK", "solver hanging", "restart with no more than 1 proc"). Its
empirical basis is gone: the hang risk was the pre-7.1.2 upstream bug plus
the BLR-noise thrashing fixed above, and the 50k/80k-per-proc thresholds
were calibrated in that same era. Removed the function, declaration, call
site, and now-orphaned `<algorithm>`/`<iostream>` includes; also the
Joan G. Stark art credit in this file's header (the turtle was its only
jgs piece — Joan personally approved BELFEM's use of her art, and the
oversubscription turtle in `Wrapper::hatch_turtle` remains, credited in
its own file). The OpenMP oversubscription warning stays — that failure
mode is real and current.

## Open / follow-ups

- Optional hardening: warn when `Krylov_iterations()` hits `maxit`
  (closes the silent-stagnation hole for users who opt into BLR).
- If BLR is ever wanted for truly large systems, benchmark the crossover
  N and the tolerance/leaf-size trade-off before re-introducing any
  automatic heuristic.
