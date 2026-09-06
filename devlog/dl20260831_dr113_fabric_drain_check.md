# DR-113: the comm fabric gets its tripwire — comm_drain_check()

**Date:** 2026-08-31
**Topic:** preventive hardening of the point-to-point fabric against the DR-100
defect class ( register option ( a ); option ( b ) tag separation stays deferred ).

## What landed

The two-tags-per-rank-pair fabric relies purely on execution order, so a stray
message never fails where it is born — it poisons an innocent recv an arbitrary
distance later ( DR-100: one month from defect to first crash ). The hardening
makes the DR-100 probe ladder permanent:

- **`comm_drain_check( const char * aLabel )`** in `commtools` ( declaration
  `commtools.hpp`, definition `commtools.cpp` next to `comm_tag()`, whose
  contract it enforces; all vendor calls stay inside the wrapper layer, L-21 ).
  Gate `#if defined( BELFEM_MPI ) && BELFEM_ASSERTIONS_ACTIVE` — compiled out
  in release, no-op in serial, guarded against a pre-init `gComm`.
- Mechanism: `comm_barrier()` closes the exchange; each rank `MPI_Iprobe`s both
  fabric tags of every pair it belongs to ( explicit source and tag, never
  `MPI_ANY_TAG` — MUMPS/PARPACK share `MPI_COMM_WORLD` ); the first-hit verdict
  is folded through the existing `allreduce` ( MPI_MAX ); on a global hit EVERY
  rank raises `BELFEM_ERROR` — detectors name label/source/tag/kind/bytes,
  the others point at the peer's error box.
- The allreduce is the load-bearing trick: it replaces the trailing barrier as
  the fence against the next exchange's sends AND makes the failure symmetric,
  so the debug throw policy ( per-rank lldb, reaction never rank-dependent )
  survives with no rank stranded in a barrier. Codex round 1 wanted a direct
  `comm_abort` here; that was rejected as a policy breach, and Grok's allreduce
  design shipped instead.
- Five chokepoints, all existing collective contexts: `DofManager::solve` and
  `solve_from_residual` ( the solve/distribute boundary ),
  `DofManager::postprocess` ( where DR-100's stray was born ),
  `SolverData::collect_matrices` entry ( where it struck ),
  `Controller::finalize` exit ( the per-timestep boundary — finalize is thereby
  collective in assert-active builds ).

## Audit trail

Full three-stage round in `tmp/ai_exchange/dr113_fabric_hardening.md`, both
auditors at gpt-5.6-terra / grok-4.6, effort xhigh ( safety boundary: MPI
collectives ). Plan round: both accept-with-changes — the shared load-bearing
catch was that my planned "BELFEM_ERROR → MPI_Abort kills peers" story is false
where the check exists ( debug throws, `gThrowOnError = BELFEM_ASSERTIONS_ACTIVE` ).
Code round: Codex accept-with-changes ( label null-guard, loop-variable prefix,
finalize contract note — all applied ), Grok accept. Final hunks archived as
`tmp/ai_exchange/dr113_code.diff`.

Known residuals, recorded in the register row: a `comm_check` throw mid-sweep is
not collective ( only reachable on a broken MPI ); a surplus send born inside
`collect_matrices` before `collect_rhs_vector` consumes it is not covered by the
entry drain. The check is best-effort by nature — a stray still in flight can
escape one sweep; an abort is always real, a clean pass proves nothing.

## Evidence status

**G0 and G1 met; the load-bearing gate is still owed.**

- **G0 ( build )** — MET. `make check` builds the debug tree with the change.
- **G1 ( suite )** — MET. `make check` passes ( Christian, 2026-08-31 ), so the
  change regresses nothing.
- **G2 ( no false positive at a live chokepoint )** — **NOT met, and `make check`
  does not touch it.** Of the whole suite only `tests/comm` registers
  `TESTRANKS` ( `2 4` ); every other suite, `tests/fem` included, runs at one
  rank, where `comm_drain_check` returns at the `size() < 2` guard before any
  collective. The np=2 / np=4 binary exercises `commtools` primitives and never
  constructs a `DofManager`, `SolverData` or `Controller`, so **no test executes
  any of the five chokepoints with more than one rank** ( grep: no
  `comm_drain_check` under `tests/` ). A green suite here is evidence about the
  build and about everything else in the tree — not about this check's live
  path. Still owed: a parallel gate deck ( RLC np=4 or 2D_Tapestack np=2, debug
  binary ) through several saves with zero drain aborts.
- **G3 ( the tripwire actually fires )** — owed. A deliberate unmatched send in
  a scratch build; every rank must error and name the nearest boundary.

Static gates behind the above: strict syntax pass
( `-Wall -Werror -pedantic-errors`, debug-tree flags ) on all four touched
translation units, plus release-shape and serial-shape passes on
`commtools.cpp`.

## Records

- `todo/debt_register.md` DR-113 row updated ( option ( a ) landed, gates owed;
  ( b ) still deferred ); `check_doc_claims.py` 37/37 after the edit.
- `src/comm/doc/comm_usage_guide.md`: new `comm_drain_check` section ( Codex
  language sweep applied ).
