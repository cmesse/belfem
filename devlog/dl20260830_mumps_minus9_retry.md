# Devlog 2026-08-30 — MUMPS −9 Workspace Retry + DR-144 R12/R13 (Joint Round)

**Date:** 2026-08-30
**Topic:** DR-106's in-wrapper `-9` retry and DR-144's two deferred steps, landed as one round over `cl_SolverMUMPS.cpp`
**AIs involved:** Claude (plan + implementation), Codex (jury, terra high→xhigh), Grok (jury, 4.6 high→xhigh)
**Claude Confidence:** high on the design, anchored to the MUMPS user guide; the retry itself is unexercised until R8
**Codex/Grok Audit Confidence:** both high on the structure; Grok ~75% "safe to land" bounded by the then-unread PDF (since read by Claude) and the missing run gate
**Literature References:** MUMPS 5.9.0 user guide (error −9; JOB 5 definition; ICNTL(14) phase note; INFOG rank-uniformity §"all processors would return with INFOG(1) = −8"); linked library verified 5.9.1 via `dmumps_c.h`
**Verification:** `g++/gfortran -fsyntax-only` clean on all touched TUs under Armadillo, Blaze, and no-MUMPS define sets — below the compile/link rung. **Reviewed, not verified.** R7 (`make check`) and R8 (a deck rerun reproducing −9, tape_quench BDF5/mumps) are owed and recorded in the todo; the tree was held by the live `tape_quench_usermat` run all session.

## Summary

One plan+audit → implement+audit cycle closed DR-106's code work and DR-144's R12/R13, which edit
the same two solve overloads. Process: O1/O3/O4 resolved by Christian up front (ICNTL(23) out;
cap 240 persist-bounded; no deck key), design written into `todo/dr106_dr144_joint_session.md`
§3.1, blind plan jury, amendments, implementation, blind code jury at xhigh (MPI-collective/ABI
safety boundary), four prescribed fixes applied, module guide updated + swept.

## Key Findings

- The shim exported only rank-local `INFO`; MUMPS puts the true code on the failing rank and `−1`
  elsewhere, so any retry keyed on it is a strand/deadlock. Fix: export `INFOG(1:40)`
  (`mumpstools.f90`, new `aInfoG` argument), whose (1:2) are documented rank-uniform.
- Plan jury must-fixes that survived verification: five call sites not four (Blaze `#else` arm);
  do NOT retarget `mumpstools_create_solver` to `INFOG(1)` (unproven after a failed `JOB=−1`;
  wrong rollback key would reopen DR-143's slot leak); `free()` never restored `MemoryRelaxation`
  (only the ctor writes it); `error_message` had to move to `INFOG` with `print_soft_fail`.
- Retry job adjudicated against both auditors' conditional objections with the manual open:
  JOB 5 (= 2+3) is the documented −9 recovery through this shim; JOB 2 alone would return the
  RHS copy as "solution" (the shim pre-copies `aY→aX`); ICNTL(14) is read at BOTH analysis and
  factorization, so the raise needs no re-analysis. A rejected JOB 5 (−3) exits the ladder into
  the soft-fail arms: status quo, no hang.
- Code jury catches, all fixed: the `free()` restore sat outside `#ifdef BELFEM_MUMPS` while
  `mIParameters` is sized only inside (release UB in a no-MUMPS build, unreachable today); the
  −9 box printed raw `INFOG(2)` (millions convention when negative) — shortage dropped from the
  box; `MemoryRelaxation == 0` would loop the collective forever (guard + clamp added); the
  controller string claimed a retry STRUMPACK/PETSc do not have (reworded, 4 sites).
- Split verdicts adjudicated: `Vector<int_t>` kept for `mInfoG` (matches sibling `mInfo`;
  header comment documents the deviation). G-E's "third-category miss" framing softened per Grok:
  the 8-strike abort caps an already-correct retry; the defect was the message.
- Noted, not touched: `error_message`'s −29 arm reads index 22 written for rank-local INFO
  (pre-existing, low); `mumpstools_free_solver` still exports rank-local INFO (teardown);
  matrix-RHS overload does not store SymmetryMode (pre-existing); OOM risk of a persisted 240
  is the accepted O3 design, ICNTL(23) the future lever.

## Changes Made

- `src/sparse/mumpstools.f90` — `aInfoG` export; copy-out comment carries the manual quote
- `src/sparse/mumpstools.hpp` — widened `mumpstools_solve` prototype (shim + wrapper are ONE rebuild: `bind(c)` argument-shift smash if either side is stale)
- `src/sparse/cl_SolverMUMPS.hpp/.cpp` — `mInfoG` member; anonymous-namespace relaxation constants; shared `escalate_workspace()` ladder in both overloads (base restore after the loop); error arms/decodes on `INFOG`; G11 bails hard-error unless soft-fail; `free()` restores the default under `#ifdef BELFEM_MUMPS`; stale `f90:307` comment fixed
- `src/fem/kernel/cl_FEM_Controller.cpp` — 4 abort strings no longer say "persistently singular"
- `src/fem/kernel/cl_FEM_DofMgr_EigenValues.cpp` — G10 unconditional `mK` recapture + null error
- `src/sparse/doc/solver_memory_and_compression.md` — three stale claims fixed, Codex-swept
- Trackers: `todo/dr106_dr144_joint_session.md` (status + boxes + O-closures), `todo/dr144_mumps_lifecycle_plan.md` (R12/R13 ticked, 11/11), `todo/debt_register.md` (DR-106 and DR-144 retagged `[CODE]`→`[RUN]`, resolutions appended); `check_doc_claims.py` 37/37 after

## Open Questions

- R7/R8 gates (Christian, once the tree is free). R8 is the only gate that exercises the retry.
- §7 of the joint todo: the 148 400 vs 86 215 free-dof discrepancy remains undiagnosed, not this row.
- O5: does `−8` join the ladder if ever observed (same ICNTL(14) remedy).

## Files Updated

- src/sparse/{mumpstools.f90, mumpstools.hpp, cl_SolverMUMPS.hpp, cl_SolverMUMPS.cpp}
- src/fem/kernel/{cl_FEM_Controller.cpp, cl_FEM_DofMgr_EigenValues.cpp}
- src/sparse/doc/solver_memory_and_compression.md
- todo/{dr106_dr144_joint_session.md, dr144_mumps_lifecycle_plan.md, debt_register.md}

## Addendum — same day: gates ran, both rows struck

Christian rebuilt ( shim + wrapper in one pass, 14:49–14:50 ) and restarted `tape_quench_usermat`
on the new binary. The restart supplied R8 at the top rung of the evidence ladder:

- **Step 791 ( t = 109.18 ms ), magnetic solve:** `-9` at ICNTL(14) = 30, ladder walked
  60 → 120 → 240 ( raw INFOG(2) shortfall 214 634 → 3 421 817 → 8 327 358 — it grows with the
  relaxation because more pivoting is admitted ), succeeded at 240, and the step **converged**
  ( magnetic 5.95e-08, thermal 1.00e-08, "timestep succeeded" ). This exact class of event was a
  timestep cut before today.
- **Eigen/conditioning probe:** its own MUMPS instance laddered independently to 240 and the
  analysis completed ( 6013 ms, COND1 5.55e9 ). The second ladder starting from the default is
  correct, not a persistence failure — that instance is created and freed per call, and
  persistence is per-instance by design.
- **Zero soft-fail boxes** in the whole log; run healthy per Christian.
- `check-fast` ran 15/15 suites green incl. `sparse` at 14:50 ( full suite rides the nightly CI
  once committed ).

**DR-106 and DR-144 struck and archived** on this evidence ( Christian's ruling );
`check_doc_claims.py` 37/37 after the strike. Log citations:
`cmake-build-debug/tape_quench_usermat/out.txt:34758-34812`,
`cmake-build-debug/Testing/Temporary/LastTest.log`.

## Addendum 2 — same day: the ladder's log hygiene

The restart above was healthy, but the `-9` events printed through the frame of the timestep box.
Two independent leaks, both fixed:

**1. MUMPS' own lines ( `mumpstools.f90` ).** Every soft `-9` emitted

```
 On return from DMUMPS, INFOG(1)=              -9
 On return from DMUMPS, INFOG(2)=         2327741
```

despite `ICNTL(4) = 0` and `ICNTL(1) = -1`. Root cause read out of `tmp/MUMPS_5.9.1/src/`:
the ERRORG section at `dmumps_driver.F:2360` is guarded by

```fortran
      IF (id%MYID.EQ.MASTER.and.MPG.GT.0.and.
     & id%INFOG(1).lt.0) THEN
```

`MPG = id%ICNTL(3)` ( `:749` ), and unlike `PROKG` ( `:752-753` ) this site carries **no**
`ICNTL(4)` test — so `ICNTL(4) = 0` does not reach it and neither does `ICNTL(1)`. Five further
`IF ( MPG .GT. 0 )` warning sites in the same file share the pattern ( `:1529, 1801, 1819, 1834,
1854, 1868` ). Silencing the *stream* is the only lever. `ICNTL(1)`, `ICNTL(2)` and `ICNTL(3)` now
all follow the info level together: real units at `-v 4` and up, `-1` below. The `JOB = -1` create
was never affected — the driver leaves `MPG = 0` for that job ( `:508`, and the `ICNTL` read block
at `:744` excludes `JOB = -1` ).

**2. Our own rung report ( `cl_SolverMUMPS.cpp`, `escalate_workspace` ).** The
`MUMPS out of workspace ( -9 ): retrying with ICNTL(14) = %i` line was a bare `message()` and
punched the same hole. It is now drawn as a box row, `"   │%-71s│"` — the identical inner width
`print_soft_fail()` writes into, verified against the controller's top border ( 76 columns total ).

Not changed, deliberately: the two `MUMPS returned a warning on proc %i ...` sites
( `cl_SolverMUMPS.cpp:397, :1863` ) are still unboxed. Their payload is a variable-length warning
string that a 71-column field would truncate, and a truncated warning is worse than a ragged line.

Verification tier: **reviewed**, not verified — `gfortran -fsyntax-only` on the shim and
`g++ -std=gnu++17 -fsyntax-only` on the wrapper ( real `flags.make` defines and includes ) both
clean, and the box widths checked by rendering the rows. The executable gate is the next run that
takes a `-9`; Christian builds.
