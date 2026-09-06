# MUMPS Solver-Pool Exhaustion: the tapestack3d Step-8 Abort

**Date:** 2026-08-29
**Purpose:** Diagnose and fix the tapestack3d run abort at timestep 8 ("couldn't
initialize MUMPS solver"); harden the diagnostic path that triggered it.
**Module:** `src/sparse` (mumpstools registry, MUMPS wrapper), `src/fem/kernel`
(shift-invert conditioning), `tests/sparse`
**AIs involved:** Claude (diagnosis + fix), three peer Claude sessions (cross-checks),
Codex + Grok (plan audit and code audit rounds)
**Exchange:** `tmp/ai_exchange/mumps_pool_recycling.md` (pre-registration, both audit
rounds, reconciliation), `.../mumps_pool_recycling.diff`

## The failure

`cmake-build-debug/tapestack3d`, 8 ranks, coupled h-φ/T, magnetic on MUMPS, thermal on
PETSc, conditioning diagnostic on for both. Seven timesteps ran clean; at step 8 all
ranks aborted in `MUMPS::initialize` via `Controller::finalize →
compute_conditioning → EigenValues::run_shift_invert → Solver::solve`.

## Root cause (DR-143)

The mumpstools instance registry (default 8 slots) allocated IDs monotonically:
`gNumSolvers + 1`, freed slots never rescanned, counter reset only on full drain. The
production magnetic MUMPS instance pinned slot 1 for the whole run, so the drain never
fired; the per-timestep shift-invert conditioning fallback for the PETSc thermal field
created and freed one MUMPS instance per step, burning IDs 2–8 over steps 1–7. Step 8's
create found the pool "full" (7 of its 8 slots actually free) and returned −1 into a
bare `BELFEM_ERROR`.

The arithmetic pinned the diagnosis: 1 long-lived instance + 7 churned IDs = failure on
exactly the 8th step, and two-long-lived would have failed at step 6. A latent allocator
defect, not a regression in the conditioning fallback — nothing before the per-timestep
conditioning redesign ever freed a MUMPS instance mid-run while another lived. A peer
session confirmed the previous night's binary (660294f6) had no shift-invert path at
all, which is why its serial run showed only a soft ARPACK miss instead.

Secondary defect: the create-site error ignored `soft_fail()` while both solve paths
honour it, so a diagnostic whose own comment says "must not kill the run" killed the
run.

## The fix (three pieces, one round)

1. **Allocator** (`mumpstools.f90`): first-free scan over `gOccupied` is now the
   authority; `gNumSolvers` is an occupancy count recomputed on every free; an
   `allocated()` guard covers the whole free body (closes a latent post-drain
   out-of-bounds the recount would otherwise hit); the freed struct's user pointers
   (`irn`/`jcn`/`A`/`rhs`) are nullified so a recycled slot's `JOB = -1` never sees the
   previous tenant's pointers.
2. **Consumer** (`cl_FEM_DofMgr_EigenValues.cpp`): `run_shift_invert` keeps its MUMPS
   instance across timesteps instead of create+free per call. Bonus: the stable matrix
   pointer makes `select_job` return JOB 5, so the per-step analysis phase is also
   saved. The last step's factors staying resident is the documented, accepted price of
   the opt-in diagnostic.
3. **Soft-fail-consistent create** (`cl_SolverMUMPS.cpp`): the wrapper is marked
   initialized only AFTER a successful create; a failed create under `soft_fail` flags
   and returns with the flag false (so the next call retries — after the allocator fix
   a slot may have freed); both solve overloads bail flagged on `mSolverID <= 0` before
   touching the Fortran registry. This shape is Grok's: the first sketch was refuted
   because `is_initialized()` reads a flag set at the top of `initialize`, so the solve
   assert could never catch a failed create and release builds would have indexed
   `gSolvers(0)`.

Regression test `SolverLifecycle.MUMPSPoolRecyclesFreedSlots` (pinned instance + 12
create/free churns) is the only executable gate that exercises recycling — the consumer
fix removes the production trigger, a verification hole Grok caught in the plan round.

## Method notes

- Three peer Claude sessions cross-checked the diagnosis before the fix round; the
  cold reader corrected the secondary defect's mechanism and found the post-drain
  hazard, and the exchange records which refutations survived (one of theirs did not:
  it had checked HEAD where the claim was about 660294f6).
- Both vendors audited plan and code. The plan round's blocking find (Fix C dual-flag
  unsafety) came from Grok; Codex independently required the both-overloads coverage.
- Pre-existing lifecycle gaps surfaced by the audits (failed `JOB = -1` keeps its slot;
  `tInfo` ignored; init-then-free-before-solve skips the Fortran free) were kept out of
  the diff and filed as DR-144.

## Code round completed (later the same day)

The code round had shipped with only one vendor on it: Codex audited the landed diff,
Grok had audited only the plan. The Grok leg was dispatched afterwards against a diff
verified byte-identical to the one Codex read, and returned the same verdict — no
blocking defect, nothing to edit mid-round. Its tightenings sharpen rather than reverse:

- **Pool-exhaustion create is non-collective.** The `-1000` branch returns without the
  `MPI_BARRIER` the success path takes, so a pool that was mixed-occupancy across ranks
  would deadlock inside Fortran, ahead of every C++ guard Fix C adds. The old high-water
  allocator had the identical shape, so this is inherited, not introduced — but it is
  the concrete form of the deadlock risk both vendors rated only medium-high. Filed to
  DR-144.
- **`report_unavailable` mislabels a solver failure.** Its reason ternary has no
  `EigenOutcome::SolverFailed` arm, so a MUMPS soft-fail — including registry
  exhaustion — is reported as "no legal Krylov subspace fits the basis budget". That is
  the same misdirection class that made DR-143 itself expensive to find. Filed to
  DR-144; it is an executable change and gets its own round.
- **The regression test is a ceiling test, not a recycle assertion.** It proves twelve
  churns survive a pool of eight; it never asserts the ID was reused, and it would pass
  on the OLD allocator if anyone raised `gMaxNumSolvers` to 13 or more. It also does not
  exercise Fix C, since with recycling the soft-fail branch is never reached. Sharpening
  it needs a C++ accessor for pool state — not worth a round today, but recorded so
  nobody cites the gate for more than it covers.

Two claims were corrected in this pass. The devlog's earlier reading that a failure on
the old allocator aborts the process is wrong: `tests/sparse/test_sparse_main.cpp:22`
sets `throw_on_error`, so it surfaces as an ordinary gtest failure. And the comment at
`mumpstools.f90:154-156` overstates its case — `NULLIFY` is legal and correct on the
freed struct, but slot-reuse safety rests on the `JOB = -2` / `JOB = -1` sequence, not on
the four-pointer belt, which touches only the components BELFEM itself associates.

## Status

Fix landed in the working tree, reviewed not verified. Both vendors have now audited
both rounds. Syntax gates: all four TUs pass `-fsyntax-only` with build-tree flags
(below the compile rung). Owed gates: `make check` green including the new lifecycle
test (wiring), then the tapestack3d deck rerun past step 8 with per-step thermal
conditioning numbers (reproducer). Audit verdicts and the citation-verification table
are appended to the exchange file; see the reconciliation entries there.
