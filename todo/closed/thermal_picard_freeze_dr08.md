# Thermal Picard Freeze and the Roundoff-Sensitive Picard/Newton Switch (DR-08)

**Date:** 2026-08-12
**Purpose:** Record the first clean reproducer of the DR-08 coupled freeze, and
propose a fix for both the freeze and the mechanism that exposed it
**Module:** `src/fem/kernel` (Controller Picard/Newton handoff), `src/sparse`
(PETSc convergence reporting)
**Status: CLOSED (2026-08-13, Christian's ruling — DR-08 closed in the register).**
R1 verified live (KSP reason + iteration counts printed through the whole tapestack3d
campaign and drove the ASM/GAMG finding); R2/R3 ran in the production binary with the
freeze never recurring over hundreds of coupled steps. R1b and R4 are retired unrun with
the close; R5's scope note and O1 stand below for whoever next touches the handoff.

Original header: PLAN — reproducer captured, root cause of the freeze NOT yet
established. No code changed.
**Register:** `debt_register.md` DR-08 (P1, open since 2026-07-22, blocking,
gate recorded as "freeze case deck" — this is that deck).

---

## 1. The reproducer

Two runs of `cmake-build-debug/tapestack3d`, warm-started from the **same**
`memdump.hdf5` at t = 10.0000 ms, timestep 6, Δt = 3.0240 ms. Identical deck,
identical mesh, identical binary. The **only** difference is the decomposition:

| | ranks × threads | thermal at Picard 2 | switch | outcome |
|---|---|---|---|---|
| A | 8 × 2 | 0.000**100** (−40.01 dB) | ≤ 1e-4 → **Newton** | converges in 6 its, 3 min 20 s |
| B | 4 × 4 | 0.000**102** (−39.93 dB) | > 1e-4 → **Picard** | **frozen**, no progress |

Both runs are **bit-identical** through `Magnetic Picard 2` (−156.54 dB). They
diverge first at the thermal residual, by about 2 % — the ordinary consequence of
an assembly summing in a different partition order.

In run B the thermal residual is then **identical to every printed digit** at
iterations 2, 3, 4 and 5:

```
Thermal Picard 2, residual 0.000102 ( -39.93 dB ), relax 1.00000
Thermal Picard 3, residual 0.000102 ( -39.93 dB ), relax 1.00000
Thermal Picard 4, residual 0.000102 ( -39.93 dB ), relax 1.00000
Thermal Picard 5, residual 0.000102 ( -39.93 dB ), relax 1.00000
```

The magnetic side is promoted to Newton in both runs and sits at machine
precision throughout, so it is not the differentiator.

## 2. Two separate defects

### D1 — the freeze itself (root cause NOT established)

An exactly constant residual across four iterations means the thermal state is
not moving. Three candidate mechanisms, not yet distinguished:

1. **The linear solve returns its input.** `KSPSetTolerances` sets the KSP
   *relative* tolerance to the solver's own `mEpsilon`
   (`cl_SolverPETSC.cpp:393-399`), and `KSPSetInitialGuessNonzero` can be TRUE
   (`:510`). If the initial guess already satisfies rtol, PETSc returns after
   zero iterations with the vector unchanged, the Picard update is zero, and the
   next nonlinear residual is bit-identical. **This matches the observation
   exactly and is the leading hypothesis.**
2. **A genuine Picard fixed point that is not a root** — the lagged operator has
   a fixed point where the full nonlinear residual is 1.02e-4.
3. **A stale print** — the thermal was not updated and the previous residual was
   re-displayed. Judged unlikely here: `tUpdateThermal` gates on the *magnetic*
   residual (`:1200-1202`), which is identical in both runs, so it cannot explain
   why A progressed and B did not.

**Nothing currently distinguishes these**, because neither the KSP iteration
count nor its converged reason is ever read — see D3.

### D2 — the promotion gate can deadlock, and is roundoff-sensitive

Promotion to Newton requires `mEpsilon2 < mEpsilonSwitch2`
(`cl_FEM_Controller.cpp:1215-1219`). But when Picard cannot reduce the residual,
Picard is exactly what has to get below the threshold for Newton to be allowed to
run. **Newton is gated behind progress that only Newton can make.**

Independently: comparing a floating-point residual against a fixed constant makes
the *algorithm* a function of the partition count. Run A and run B executed
different solvers on identical input. That is a reproducibility defect in its own
right, separate from whether either result is correct.

The code already knows this region is fragile — `:1254-1258` states "the switch
has no hysteresis, so a residual orbiting the tolerance switch flips the
algorithm every few iterates", and there is a three-demotion chatter latch. Those
mitigations target *oscillation across* the threshold. This case sits just above
it and never oscillates, so none of them fire.

### D3 — PETSc convergence is never checked

`PETSC::solve` inspects the `PetscErrorCode` from `KSPSolve`
(`cl_SolverPETSC.cpp:168-201`) but never calls `KSPGetConvergedReason` or
`KSPGetIterationNumber` — neither appears anywhere in the file. A solve that hits
`maxits`, or diverges, returns code 0 and is accepted silently. Same class as
DR-45's MUMPS `INFO(1)=+1`: the library reports a problem through a channel
nobody reads.

## 3. Proposed fix

**R1 — instrument before changing behaviour (do this first).**
Add `KSPGetConvergedReason` + `KSPGetIterationNumber` after `KSPSolve`. Hard-error
on negative reasons; report the iteration count at Detailed level. This is
independently justified, cheap, and it settles D1 by showing whether the frozen
solves are returning after zero iterations. Do not guess at D1 before this runs.

**R2 — break the promotion deadlock (stall escape).**
Promote to Newton when Picard *stalls*, not only when it crosses the threshold:
if the thermal residual improves by less than a small factor over N consecutive
Picard iterates, promote regardless of the absolute value. The machinery already
exists — `mBestEpsilon2` / `mBestEpsilonIteration2` (the watchdog re-anchor) and
`mThermalFlatCount` — so this is a new promotion clause, not new bookkeeping.
Contract with R3: an escape promotion must NOT reset the chatter latch.

**R3 — hysteresis on the switch.**
Promote at `eps < epsSwitch`, demote only at `eps > k * epsSwitch` with k of
order 3-10. Removes the roundoff sensitivity of the *demotion* edge and makes the
existing chatter latch a backstop rather than the primary defence. Note this
alone would NOT have fixed run B, which never crossed the threshold at all —
R2 is the load-bearing change; R3 is hygiene.

**R4 — decide what reproducibility across rank counts is worth.**
Even with R2 and R3, residuals differ by roundoff between decompositions and any
threshold test can select differently. If bitwise-comparable behaviour across
rank counts is wanted, that is a larger conversation (deterministic reduction
order); if not, R2's stall escape at least ensures both decompositions reach the
same *algorithm* eventually, and the answer is set by the nonlinear tolerance
rather than by the partitioner.

## 4. What this invalidates

The 8×2 vs 4×4 throughput experiment run the same evening is **confounded** and
must be repeated. The raw step times (3:20 vs 5:26) compare 6 iterations against
9, because one run was promoted to Newton and the other was not. Per iteration
the difference is ~8 % (33.3 s vs 36.2 s), not the ~39 % the step times suggest.
Re-run with `tolerance switch` set well away from where the residual lands, so
both decompositions take the same path.

## 5. Steps

- [x] **R1** instrument `KSPSolve` with converged reason + iteration count —
      **done 2026-08-12** (`cl_SolverPETSC.cpp`, after the error-code
      coordination block): `KSPGetConvergedReason` + `KSPGetIterationNumber`,
      `BELFEM_ERROR` on negative reasons (uniform across ranks — the reason is
      a property of the collective KSP object, so erroring everywhere cannot
      deadlock), iteration count at Detailed level, rank-0 gated because
      `message()` has no rank gate of its own.
- [~] ~~**R1b**~~ retired unrun with the DR-08 close (2026-08-13) — was: re-run the 4×4 warm restart at Detailed level and record whether
      the frozen solves return at zero iterations ( gate: Christian runs )
- [x] **R2** stall-escape promotion — **done 2026-08-12, as a PORT rather than
      a new clause**: the magnetic side already had exactly this mechanism
      (`try_escalate_to_newton` / `mForceNewton`), so the fix is its thermal
      twin `try_escalate_thermal_to_newton` / `mForceNewton2` /
      `mNewtonEscalated2`, called from the existing stagnation detector at
      **two** flat iterates ( the >= 5 stall warning stays as the backstop ).
      **Deviation from §3 as planned, following the magnetic precedent:** the
      escalation CLEARS `mJustPicard2` rather than respecting it — the magnetic
      twin does the same, and the ping-pong the latch guards against is capped
      here by `mNewtonEscalated2` ( once per attempt ) instead. Same lifecycle
      as the magnetic pair: cleared in `initialize_timestep` /
      `initialize_thermal`; same `mEpsilon2 < 1E1` divergence guard.
- [x] **R3** hysteresis on the demotion edge — **done 2026-08-12**: a running
      Newton is kept unless the residual retreats a full decade above the
      switch (`10.0 * mEpsilonSwitch2`, hardcoded, no input key ). The
      fresh-timestep reset of `mEpsilon2` to `BELFEM_REAL_MAX` sits far above
      any band, so first iterates stay Picard; the post-freeze Picard restart
      ( ts12 ) keeps priority via `!mThermalFrozen` in the keep condition.
- [ ] **R5 ( NEW, found in flight 2026-08-13 ):** R2's stagnation gate cannot
      fire in a genuinely COUPLED breakdown. `try_escalate_thermal_to_newton`
      is called from the flat-residual detector, which gates on the MAGNETIC
      side already being converged
      ( `mEpsilon <= mRelativeEpsilonTarget || mEpsilonAbs <= mAbsoluteEpsilonTarget` ).
      Observed on `tapestack3d` at 1.05 Ic, step 577: both physics stuck with
      residuals oscillating at their relaxation floors — magnetic ~-28 dB
      ( relax 0.005 ), thermal ~-6 dB ( relax 0.011 ) — for 22+ iterates.
      `mThermalFlatCount` never increments there, so the thermal escalation is
      unreachable by construction. What carried the step was the pre-existing
      MAGNETIC `try_escalate_to_newton` plus a rejection and a Δt halving;
      the retry then converged in 8 iterates.
      **R2 is therefore correctly scoped to "thermal alone lags a converged
      magnetic side" and does NOT address coupled stalls.** Decide whether a
      coupled-breakdown escalation is wanted at all, or whether
      escalate-then-cut on the magnetic side is the right and sufficient
      answer — the evidence so far is that it is: every breakdown this run
      recovered after one cut.
- [~] ~~**R4**~~ retired with the DR-08 close (2026-08-13) — cross-rank reproducibility accepted as roundoff-level, not owed a ruling
- [ ] repeat the 8×2 vs 4×4 throughput experiment once the path is stable
- [ ] **O1 ( new, out of scope here ):** the staggered `iterate_thermal` path
      has its own copy of the handoff and did NOT receive R2/R3 — the
      reproducer lives in `iterate_coupled` and the controller convention is
      minimal edits. Port when the staggered driver next gets attention.

**Status of the code changes: reviewed, not verified.** All three TUs pass
`-fsyntax-only` under the tree's own `-Wall -Werror` in both the release and
the assert-active configuration. Nothing has been built or run; the gate is
R1b plus a re-run of the 4×4 reproducer, which after R2 must promote to Newton
by the third flat iterate instead of freezing.
