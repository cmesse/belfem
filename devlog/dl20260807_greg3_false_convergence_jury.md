# greg3 false convergence — three-AI jury round + live confirmation

**Date:** 2026-08-07
**Purpose:** Detective jury round on the greg3 tapestack oscillations
(checkerboard blobs along the thin shell, `bfield.png`/`bfield2.png`) and
the suspicious instant convergence (−124.43 dB at Picard 1). Read-only;
Christian ran the decisive experiment same day and confirmed the verdict.
**Module:** fem/kernel (SolverData, Controller, Anderson mixing)
**Thread:** `tmp/ai_exchange/review_greg3_false_convergence.md`
(pre-registration frozen before dispatch; blind Codex + Grok jury on the
neutral brief `tmp/ai_exchange/greg3_false_convergence_brief.md`;
verification + reconciliation appended).

## Verdict (P0, 3/3 independent, source-traced, CONFIRMED by A/B run)

**The reported residual under Anderson mixing is the linear solver's own
roundoff, not a nonlinear convergence measure.** Chain:

1. The 2026-08-07 opt-out defaults (`dl20260807_timestep_input_optout.md`)
   arm Anderson depth 3 + BDF5 on any deck without the keys — greg3 has
   neither (`cl_FEM_Controller.hpp:140,163`, unconditional forward at the
   end of `set_params`).
2. `anderson_update` refreshes `mFieldValues` to the mixed iterate
   (`cl_FEM_DofMgr_SolverData.cpp:2771-2780`), then the residual multiply
   `:2289` uses the SAME lagged matrix that was just solved. On the
   bootstrap iterate (empty window, β = ω = 1) the mixing step returns
   exactly x = A⁻¹b (`fn_FEM_anderson_mixing.hpp:55-56`), so
   ε = ‖A·A⁻¹b − b‖/‖b‖ ≈ 3.6e-13 = −124.43 dB — by algebraic identity,
   independent of nonlinear consistency. The legacy depth-0 path
   deliberately reports the PRE-update residual and is immune; the refresh
   (added for the ts16 promotion-lag fix) silently destroyed the criterion
   for every Anderson run, latent while Anderson was opt-in (jury F3 of
   `dl20260807_newton_extension_jury.md` made concrete).
3. `run_coupled`/`run_magnetic` exit at `mMinNumIterations = 2` with ε
   "converged" → one effective lagged-Picard step per timestep, true
   residual never measured → checkerboarding, the canonical
   loose-tolerance artifact (Messe et al. 2023, paper1, **§4** — NOT §2.7;
   Codex correction, the wrong §2.7 cite also sits in the
   `cl_FEM_Controller.cpp:125` comment).
4. **The log's second line is the tell (Codex, round's best find,
   Claude-verified):** at iterate 2 the acceptance reference IS the fake
   floor, so the coupled line search (hphirun drives `iterate_coupled`
   even magnetic-only) rejects every honest trial, halves ω seven times
   (1 → 2⁻⁷ = 0.0078125 → prints 0.00781), then the `tBacktracks >= 8`
   escape (`cl_FEM_Controller.cpp:960-963`) COPIES `mEpsilon = mEpsilon0`
   and accepts the restored state — reproducing the bit-identical
   −124.43 dB. Corollary: 8 wasted assemble+solve cycles per timestep,
   and the reject storm is in-log proof the honest residual was far above
   tolerance. Claude's pre-registered cross-timestep ω-decay mechanism
   was refuted (ω resets each timestep, `:117`); Grok's "unexplained"
   stance superseded by the verified trace.

**Live confirmation (Christian, same day):** rerunning greg3 with
`anderson stabilization : false` removes the oscillations. That is the
reproducer tier — verdict closed.

## Secondary findings

- `algorithm : Newton` never runs: promotion needs `mIteration > 1` and
  the loop is gone at 2 (P1, 3/3).
- BDF5 remains silently active and is a separate validation risk (jury F6
  zero-stability under ×1.5/×0.5 Δt adaptation); it cannot produce the
  residual signature (3/3). **Open lead:** in the confirmed Anderson-off
  run, Δt collapses toward zero at t ≈ 6.5 s (I ≈ 104 A = 0.65·Ic;
  deck Ic = 4·jc·0.25 µm·4 mm = 160 A = the ramp target). Candidate
  drivers: genuine n=100 stiffness as flux penetration deepens vs the F6
  BDF5/Δt-ratio hazard vs the −124 dB-era escape logic. Decisive gate =
  `scheme : bdf1` A/B (G-B1).
- Grok watchlist (single-raiser, needs adjudication): the Newton branch's
  post-update lagged residual shares the hazard qualitatively (zero only
  where J ≈ A, so not an identity — latent).
- Fix direction routed to Christian (NOT executed): pre-update residual in
  the Anderson path vs fixed-point residual ‖G(x)−x‖ vs reassembled true
  residual; and whether Anderson/BDF5 stay opt-out.

## Attribution

Codex: backtracking-escape trace (round's best find), §4 literature
correction, explicit-over-implicit flag on opt-out defaults. Grok:
BELFEM_EPS-floor exclusion, defaults table incl. `cl_MaxwellFactory.cpp:736`,
Newton-branch watchlist. Claude: pre-registration (C1 chain), verification
pass, ω-arithmetic reconciliation. Christian: the field report, the images,
and the decisive A/B run.

## Addendum (same day): fix executed on Christian's ruling

Literature consult decided the residual semantics (Messe 2023 paper1 §4 =
the validated force criterion; Bathe §8.4.4 = increment criteria alone
under-report on stiff maps; Walker & Ni monitor ‖G−x‖). Christian picked
option (a) + fixed-point signal, plus a defaults rollback. Landed
(working tree, not compiled; focused diff
`tmp/ai_exchange/anderson_residual_fix.diff`, Codex audit thread
`tmp/ai_exchange/anderson_residual_fix_audit.md`):

- **Honest ε on the Anderson path:** the `mFieldValues` refresh in
  `anderson_update` is REMOVED — ε is the pre-update force residual
  ‖A(x_k)x_k − b(x_k)‖/‖b‖, identical to depth-0 (D5 of the Anderson
  campaign reversed; the ts16 promotion-lag costs at most one extra
  Picard iterate). Newton keeps its post-update recompute (J ≠ A,
  non-degenerate).
- **Fixed-point residual exposed:** ‖G(x_k)−x_k‖/‖x_k‖ computed in
  `anderson_update`, broadcast in `residual()`,
  `SolverData/DofManager::fixed_point_residual()`. Diagnostic only — no
  controller decision consumes it yet; deliberately NOT a convergence
  gate (Bathe Eq. 8.119 caveat).
- **Defaults rolled back:** `mTimeStepping` BDF5 → BDF1,
  `mAndersonDepth/2` 3/1 → 0/0. `anderson stabilization` is now opt-IN:
  `true` fills 3/1 where no explicit `anderson depth` key exists
  (explicit key always wins, incl. explicit 0); `false` + explicit
  nonzero depth stays a hard error. The one-day opt-out experiment
  (dl20260807_timestep_input_optout.md) is thereby reverted in substance;
  `scheme` key + `galerkin` parsing stay.
- **Docs synced:** `anderson_acceleration_theory.md` (residual semantics
  rewritten, opt-in config), `nonlinear_controller_theory.md` (§1
  residual paragraph, input table), `doc/input_file_reference.md`
  (§4.2/4.3/4.4 rows + shifted :line cites) — per Christian's new
  standing rule: input-contract changes update the reference doc in the
  same turn + Codex prose pass (dispatched,
  `tmp/ai_exchange/input_reference_prose.md`). §2.7 → §4 cite fixed at
  three controller comment sites.

## Addendum 2 (same day): the 6.4 s Δt→0 hang — Picard line-search retirement

With honest residuals the run froze at t ≈ 6.405 s, Δt → 0, showing the
same fingerprint at −63.94 dB: bit-identical ε over three iterates, relax
1.0 → 2⁻⁷ → ω-floor. Root cause (source-traced, Christian approved both
fixes): the ported matfix backtracking line search is **degenerate for
Picard solves**. The pre-update ε judges the iterate's ENTRY state — which
is exactly what the backup restores — and within one frozen assembly no
measurement can rank a Picard trial (the lagged post-update residual of a
relaxed step is identically (1−ω)·r). So a regression triggers 8 blind
identical re-solves, the escape copies `mEpsilon = mEpsilon0` (the frozen
print), the state never changes, `reset_timestep` restores ω = 1 (ε < 1
branch) and the next attempt repeats — the watchdog cuts Δt forever.
Decisive literature fact: Messe 2023 paper1 §4 Eq. 14 has NO Picard line
search — α = 0.5 damping on regression, keep iterating, Δt cut only after
k_max; the adaptation tail already implements Eq. 14 verbatim (growth
β + γ·(2/π)·atan, α on regression).

Landed (working tree, not compiled; diff
`tmp/ai_exchange/picard_linesearch_retirement.diff`, Codex audit
`tmp/ai_exchange/picard_linesearch_retirement_audit.md`):

- **Picard iterates always accept** in `iterate_coupled` (commit Anderson
  pair, exit trial loop); regression handling = Eq. 14 α branch in the
  adaptation tail + divergence (+10 dB post-loop rule) / stagnation /
  watchdog guards. The line search, backtracking budget, and
  moved-baseline detector remain for NEWTON solves only, where the
  post-update residual is non-degenerate (ts34 protection preserved).
  Side effect: the ts1727 Picard reject storm under a moving thermal
  baseline can no longer occur (no Picard rejects to storm).
- **Probe `BELFEM_PROBE_QHIST`** (opt-in env, DR-25 hygiene — strip
  later): `IWG_Timestep::qhist_deviation()` = max relative L2 deviation
  between the live dof fields (`mFieldData(0)`) and the one-step history
  q° (`mFieldData(1)`), printed at the first iterate of each attempt.
  Tests the second finding: the −63.94 dB first-iterate residual is
  Δt-INDEPENDENT, which is impossible for a consistent history (as
  Δt → 0, ε → 0 since ‖b‖ → ‖Mq°‖/Δt diverges); frozen ε ≈ 4e-7 implies
  ‖x₀ − q°‖/‖q°‖ ≈ 4e-7 — a retry/savepoint hygiene defect if confirmed.
- Theory doc updated (line-search section Newton-only + degeneracy
  rationale, moved-baseline scope). No input-file contract change, so the
  input reference is untouched.

Rejected mitigation: deck-level `max relaxation : 0.5` (Christian tested
the idea's direction — heavier damping was unproductive, many more
timesteps needed; the structural fix above is the answer).

**Audit outcome (Codex, same day):** no P0. P1 confirmed + fixed — with
Picard always accepting, a doomed iterate (+10 dB / exhausted budget)
reached the coupled thermal solve before the post-solve reset;
`tMagneticDoomed` now gates `tUpdateThermal` with the reset's own
conditions (`mIteration + 1`, mForceNewton exemption; NaN caught by the
gate comparison). Doomed iterates take the freeze path, so the
resume-after-freeze Picard restart applies; greg3 (magnetic-only) is
unaffected. P2 (`aMax` naming) refuted by in-repo return-value idiom
(`aResult`, `aResidual`).

**Probe v2 (same day):** Christian's run printed `#qhist deviation : 0`
everywhere — v1 was TRIVIALLY zero by construction: it compared q° against
the live fields right after `shift_fields`, which had just DEFINED q° as a
copy of those fields (it verified the shift/reset round-trip, but corrupt
level-0 state would be copied into q° and still match). v2 measures the
real invariant: `probe_qhist_snapshot()` captures the dof fields when a
timestep is ACCEPTED (`finalize`, `!mReset`, rank 0, env-gated), and
`qhist_deviation()` compares q° at the first iterate of the NEXT attempt
against that snapshot — spanning every retry cycle. Chronology note: the
Δt-independence inference for the first-iterate residual was over-read
from a single attempt box (iterates 2/3 were escape copies, only one fresh
measurement existed); with the escape gone for Picard, the new logs show
the true ε(Δt) trajectory directly.

## Addendum 3 (same day): flux-front stall — pure Picard + persistent divergence rule

Probe v2 delivered: 169/169 zero — BDF history exactly equals the accepted
state across all retry churn; the desync hypothesis is dead. The full
out.txt also showed the first-iterate residual scaling ~linearly with Δt
(−27.9 dB @ 10 ms → −65 dB @ 1.4 µs), retiring the Δt-independence worry.
The stall anatomy instead:

1. **Newton counterproductive on this net-current deck** (bearing gauge
   mode): promotion at ε < 1e-4 REGRESSED Picard's −64.9 dB iterate to
   −60.5, flailed to the ω floor, stagnation-latched, re-crawled — 19-30
   its/step vs Picard's one-step 37 dB crash-through. Remedy applied at
   deck level: `algorithm : Picard` (matches the recorded bearing-mode
   remedy). Healthy steps dropped to 2 iterates.
2. **Flux-front ω=1 overshoot** (t ≈ 6.28 s, n=100 front jump): from
   −64 dB a full Picard step spikes to +9…+14 dB; Eq. 14 recovers it
   within 2-3 α halvings when allowed (logged: +8.91 → +5.05 → −65.9,
   converged in 8 its).
3. **The +10 dB rule was trend-blind** (jury F9 confirmed in the field):
   it cut every attempt whose overshoot peaked above ε=10 at the first
   post-min iterate — selecting by overshoot HEIGHT, not trend — so Δt
   sawtoothed at ~µs forever (100 timesteps for 0.1 ms). The paper cuts
   only on k_max.

Fix (Christian pre-approved the conditional): `mDivergenceStrikes` —
ε > 10 must persist for THREE consecutive iterates before the cut fires
(two α halvings of grace), counter reset per attempt
(initialize_timestep hygiene + reset_timestep) and on any iterate at or
below the bar; the doomed-thermal gate mirrors the pending value;
iterate_magnetic gets the same rule; NaN/k_max/mForceNewton semantics
unchanged. Theory doc divergence-rules row updated. Codex audit:
`tmp/ai_exchange/divergence_strikes_audit.md`.

## Addendum 4 (same day): field guidance + probe swept

- **`doc/input_file_reference.md` §4.5 added** (Codex-polished, one
  opening-paragraph imprecision caught: the literal `algorithm : Picard`
  default vs the paper's Picard→quasi-Newton staging are now stated
  separately): when Newton / Anderson / BDF2-5 help, when to turn them
  off, and the log fingerprints of each pathology; closing rule of thumb
  = return to the baseline (`Picard`, `bdf1`, Anderson off) before
  diagnosing anything. Table rows cross-link the section.
- **`BELFEM_PROBE_QHIST` removed completely** (both hooks, the snapshot
  member, the deviation method, and the includes added for it) after
  refuting the history-desync hypothesis 169/169. New standing policy
  (Christian): probes are stripped at session end, or Claude asks to
  remove one the moment its question is answered.
- Field state at sweep time: run past the 6.4 s barrier; front hiccups =
  honest ω=1 blowups riding Δt* (in-step recovery or clean 3-strike cut)
  plus occasional Newton promotion storms; recommended deck = pure
  Picard, optional strike budget 3→5 pending Christian's call.
