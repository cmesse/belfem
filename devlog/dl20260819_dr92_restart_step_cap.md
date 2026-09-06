# DR-92: Cap the First Warm-Restart Step

**Date:** 2026-08-19
**Purpose:** Session record of the DR-92 hunt and fix: why a warm restart at the
dumped timestep produced pathological solves, what it was not, and the one-line
cap that resolves it.
**Module:** fem/kernel ( Controller )
**AIs involved:** Claude ( investigation, code ), Codex + Grok ( theory, design and
code audits, three rounds each )
**Method:** Christian's five-step structure: theory+audit, design+audit,
code+audit, test, document plus prose sweep
**Verification:** VERIFIED BY EXECUTION. Full acceptance test passed on the
shipped binary ( four criteria, below ). Plan: `todo/dr92_restart_step_ramp.md`.
Measurement record: `tmp/ai_exchange/bnorm_anomaly.md`.

## The defect

A warm restart that re-entered at the dumped timestep produced a first
preconditioned residual four orders too large, every linear solve pinned at its
iteration cap, a Newton grind and usually a rejected step, while a **live** run
at the identical state, time and Δt was healthy.

Measured in a sandbox replica, all from dumps the campaign wrote:

| restore | first Δt | `GMRES it. 0` | Krylov | outcome |
|---|---|---|---|---|
| memdump_2350 | 0.154 ms | 0.0297 | 16 | healthy |
| memdump_2700, clamped | 5 ms | 0.903 | 16 | healthy, 2 iterates, −138.8 dB |
| memdump_2700 | 50 ms | **10851** | 50 = cap | sick, rejects |
| memdump_2800 | 50 ms | **21245** | 50 = cap | sick |
| campaign, live | 50 ms | 0.019–0.043 | ~16 | healthy |

## What Measurements Ruled Out

- **MPI** ( field-collect clobber, partition, aura, ghost ownership, reduction
  order ): a **serial** restart reproduces it, and the collect is
  `mCommSize > 1`-guarded. This also retired DR-93 as the cause.
- **BDF integrator restore**: a probe showed `mH` is exactly the dump's vector,
  correctly shifted.
- **The φ/φ0 hole at the current-injection nodes**: real, and present
  *identically* in healthy restores, so it is inert.
- **The imposed current**: correct, `max|fixed value|` = I(t+Δt).
- **Step-size ratio**: a ratio-1.00 restore at 50 ms is *sicker* than a
  ratio-10.83 one.
- **BDF order**: BDF1 drops `it. 0` from 21245 to 4.17 but still caps out, so it
  is an amplifier and not the cause.

## The mechanism

**In the measured cases the first post-restart Δt decided the outcome, and
reducing Δt afterwards did not recover the process.** The decisive measurement:
inside the *same* process, the sick 50 ms restore rejected and retried at 25 ms,
and `it. 0` stayed at 9888.

Leading hypothesis, inferred from the measurements plus a source trace and **not
independently proven**: solver-entry state derived from the process's first
Jacobian is reused on later value updates ( `reorder_internal` runs once, `update_matrix_values` reapplies the
stored permutation; BELFEM initialises once and never frees in the timestep
loop ). BDF assembles `J = αM + ΔtK`, so the first Δt decides whether that
permutation comes from a mass- or a stiffness-weighted matrix.

**Christian's contribution closed the design question.** Every dump is written at
a save point, so the last completed step is a save-point value rather than a step
the controller chose on merit, which rules `bdf_last_dt` out as the cap. And since
the *live* run picks 50 ms at that state and is healthy, the measurements point
away from dump content as the immediate fix: the shipped mitigation controls the
first solver initialization instead of adding more persisted state.

## The fix

`Controller::load_memdump` caps the re-entry step at the deck's `initial
timestep`, the step the deck already declares it can start cold from, so no new
deck key. `mDeltaTimeInitial` is stored at parse time rather than snapshotted, so
the invariant does not depend on call order. The controller then grows back under
its existing order-clamped limiter.

## Acceptance test: all four criteria passed

1. Banner prints `delta t = 5.0000 ms` instead of the dumped 50, with an
   explicit cap message.
2. First solve `it. 0` = **0.903286** in **16** Krylov iterations, reproducing
   the manually-clamped experiment to six significant figures.
3. Step converges in **2 Picard iterates** to −139.10 dB.
4. **The same process climbs 5 → 50 ms with zero rejections**, `it. 0` rising
   only 0.903 → 1.54. Five save-point snaps ( to as low as 0.12 ms ) are absorbed
   without incident.

Criterion 4 is the one that matters for sufficiency: the earlier grow-back ladder
was two processes, this is the production path in one.

## What this reverses, deliberately

The Δt half of the DR-38 restart-cliff policy: the dumped step is no longer
adopted verbatim. The BDF-history half is kept, and it is what actually prevents
the order cliff. `save_memdump` still writes `delta_time`, which remains useful
for diagnostics and for a future loader that can re-enter safely.

Cost: recovery is `initial × 1.2ⁿ`, ~13 steps on tapestack3d. Cheap against a
rejected grind; not free, and longer on decks whose `initial timestep` sits far
below their working step. That is the knob to raise if restarts feel slow.

## Method notes worth keeping

- **The audits earned their cost three times over.** Codex found four overclaims
  in the theory and supplied the solver-reuse source evidence behind the leading
  hypothesis. Grok caught that my grow-back ladder was a *second process* and
  therefore never tested state-healing, which produced the in-process
  discriminator that refuted my own T3 and T4. Grok also killed design C by
  tracing that `initialize` leaves `reordered_ = false`, making the seed idea a
  silent no-op.
- **Two of my own claims were retracted after measurement**: that the φ clobber
  was the mechanism, and that the step-size ratio was. Both were plausible,
  both were wrong, and each cost only a four-minute sandbox run to kill.
- **A reproducible sandbox turned archaeology into experiment.** Once the
  anomaly could be reproduced on demand from a dump, each hypothesis cost
  minutes instead of days.
