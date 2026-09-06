# PETSc Learns to Fail Softly; BELFEM_ERROR Learns Syslog

**Date:** 2026-08-15
**Purpose:** Record the response to the 21:19 `DIVERGED_ITS` run abort —
soft-fail routing for PETSc, the thermal deck's move to a direct solver,
a linear iteration-budget key, and the syslog hook Christian requested
**Module:** sparse, core, fem/kernel, tests
**Round:** `tmp/ai_exchange/petsc_softfail_syslog_plan.md` +
`petsc_softfail_syslog_impl.patch`

## The abort

The restarted quench run died at 21:19:
`PETSc KSPSolve did not converge ( DIVERGED_ITS after 10000 iterations )`.
Two stacked problems. The thermal ASM/Krylov solve at rtol 1e-10 is
marginal on this system — ten decades of residual reduction from a
one-level additive Schwarz preconditioner — so whether a given solve lands
under the tolerance is roundoff luck. And the convergence-reason check
added in `a8694a0e` (2026-08-12) escalates unconditionally: exactly the
anti-pattern `doc/coding_philosophy.md` §error-handling names, a
`BELFEM_ERROR("did not converge")` inside an algorithm whose controller
has a retry policy. Before that commit the same non-convergence was
*silently accepted* (the DR-45 defect class); after it, a recoverable
stall killed the run. Neither behaviour is the designed one: a failed
trial should cut the timestep.

## What landed (plan + dual audit → code + dual audit)

- **PETSc soft-fail** (`cl_SolverPETSC.cpp`): the fail decision is
  `MPI_Allreduce`d before branching — PETSc documents the reason as a
  collective property, but one rank returning while peers enter
  `collect_lhs`/`KSPView`/the barrier is a deadlock, and STRUMPACK's
  return codes already lost that exact bet. Soft-fail covers the
  iteration-class reasons only (ITS, DTOL, NULL, the two BREAKDOWNs);
  NaN/Inf and preconditioner failures stay hard, because KSP reuse after
  them is not established (would need `KSPReset` + re-set operators —
  future work if ever needed). The next attempt is clean under the
  default: the diverged iterate is never scattered back, and
  `initial guess : false` means PETSc zeros the start vector.
  `KSPGetConvergedReason`'s own failure is folded into the hard verdict.
  The controller's 8-strike abort strings (all four) no longer claim
  "persistently singular" for what may be non-convergence.
- **`max iterations` for linear sections** (`SolverParameters` → PETSc
  `KSPSetTolerances`): deck-tunable Krylov budget, validated positive
  before the uint store (the parser rounds a real), zero never reaches
  PETSc (there `max_it 0` means zero iterations, not "default"),
  `-ksp_max_it` still overrides. Carried through copy ctor and
  `synchronize()` (payload 9 → 10). STRUMPACK's refinement cap stays a
  fixed safety net, deliberately not deck-tunable (DR-76). Both input
  contract artifacts updated, with the linear-vs-nonlinear same-name
  warning spelled out.
- **Thermal deck → MUMPS** (tapestack3d `input.conf`) — **WITHDRAWN the
  same night, see the W2 section below.** The reasoning at the time: the
  thermal system is cheaply factorizable, a direct solve delivers ~1e-15
  — which would retire the rtol×T drift (DR-70/72) more completely than
  any tolerance — and the thousands-of-Krylov-iterations cost per solve
  is gone. Audited clean on symmetry, matrix format and backend
  assumptions; memory named as the one real risk and accepted at 4 ranks.
  It was the memory, and something else besides.
- **BELFEM_ERROR → syslog** (Christian's ask, adapted from a Gemini
  signal-handler template — the signal handler itself is out of scope,
  syslog is not async-signal-safe): `assert::error()` now writes the
  error-box message to syslog (identity `belfem`, LOG_CRIT, one line per
  message line, MPI rank in the payload since LOG_PID only identifies
  local PIDs). The placement survived a genuine reviewer conflict: Grok
  wanted it after the throw branch (quiet EXPECT_THROW tests), Codex
  showed debug builds THROW by default (`gThrowOnError` initializes to
  `BELFEM_ASSERTIONS_ACTIVE`), so that placement would never log in the
  builds developers actually run. Resolution: fires before the branch
  under its own `syslog_on_error()` flag, and all 14 test mains disable
  it beside their existing `set_throw_on_error( true )`. Read it back
  with `journalctl -t belfem` (or `SYSLOG_IDENTIFIER=belfem`).

## Process notes

Both plan audits were substantive: Codex caught the debug-throw
contradiction and the discarded `KSPSetTolerances` status; Grok refuted
the "early return skips only the scatter" claim with the full collective
inventory, corrected the retry-hazard analysis (the opted-in danger is
the STALE guess, not the diverged iterate), and flagged that after the
deck switch a quench run can no longer certify the PETSc soft-fail — that
executable gate needs its own PETSc deck.

## Phase-3 outcome

Both audits returned; all findings adopted same-session. Codex caught the
real hole — the fail decision was collective but the soft/hard
CLASSIFICATION was still rank-local — so both verdicts now ride one
`MPI_Allreduce( int[2], MAX )`: any rank failing fails everyone, any rank
seeing a hard reason (or a failed `KSPGetConvergedReason`, Grok's
addition) makes it hard everywhere. Grok added `KSP_DIVERGED_NULL` to the
soft set (GMRES Hessenberg breakdown, same family), fixed the rank-0
message lying about a peer's divergence, generalized the magnetic
8-strike strings the thermal fix had missed, and got the stale-guess
hazard and the full soft/hard reason partition into both contract
artifacts. Codex also caught the pre-init syslog sentinel (now "rank ?")
and the "iteration-class partition" comment overclaim. Grok explicitly
withdrew its plan-phase syslog-placement position on the facts.

## W2 failed on first contact — reverted the same night (DR-77)

The first run on the rebuilt binary refuted the thermal-direct-solve
decision within ten minutes, on both counts at once.

**Memory:** the MUMPS thermal factorization alongside the magnetic
STRUMPACK drove the four ranks to 51.6 GiB with swap exhausted (8131 of
8191 MiB) and a heavy imbalance (17.8 / 13.4 / 10.5 / 9.8 GiB). The risk
was flagged in the plan and audited as "plausible, not cheap until a run
exists" — the run existed, and it was not affordable. A monitor caught the
deterioration before the OOM killer did.

**Numerics, and this part is unexplained:** from a memdump two earlier
runs had processed identically, magnetic iterate 1 reproduced bit-for-bit
(0.068987, −11.61 dB) and iterate 2 collapsed from −39.93 dB to −3.35 dB,
with the thermal residual exploding to 3.3e2 then 7.1e2 and Δt falling
50 → 25 → 12.5 → 6.25 ms. Only the thermal solve runs between those two
iterates, so the thermal answer was wrong from the first solve. Whether
that follows from the thrashing or is an independent defect in the MUMPS
thermal path is undetermined — the plan audits checked symmetry mode,
matrix format and backend assumptions and all three came back clean, so
the mechanism is one none of the three of us modeled. Registered as DR-77
with the experiment that would separate the two causes (thermal-only peak
RSS, plus a small-deck A/B where memory is not a factor).

Deck reverted to PETSc the same night, with `max iterations : 500` added —
affordable precisely because W1 landed: exhausting the budget now costs a
rejected timestep instead of the run. The block carries a do-not-switch
comment pointing at the register. **W1, W1b and W4 are unaffected by this
reversal**; only W2 is withdrawn. The lesson worth keeping: three clean
audits of the things we thought to check did not make the change safe, and
the memory caveat we all wrote down as acceptable was the one that bit.

Status: **reviewed, not verified** — syntax gates green on every touched
TU (-Wreorder). Executable gates owed: `make check-fast`; `journalctl -t
belfem` after a forced error; the quench restart clearing step 95/96. The
PETSc soft-fail path is now exercised by the production deck again, so a
marginal thermal solve there certifies it directly.
