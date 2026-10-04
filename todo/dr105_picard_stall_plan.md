# DR-105: the Picard stall at the quench front — a floor retry that changes something

**Date:** 2026-10-03
**Purpose:** The T5 `tape_quench` deck (pure Picard, Anderson 3) aborts on the floor-retry limit at
10–33 ms into the quench, serially and in parallel, while the same deck on the hybrid completes
0.3 s twice. The weekend campaign and a blind two-auditor survey established the abort mechanism:
the stagnation guard ends every attempt at the sixth iterate, a pure-Picard deck cannot escalate,
the Δt cut walks to the floor, and the floor retry raises budgets the guard never consumes while
resetting ω, so twenty retries replay one attempt. This plan makes the floor retry respond to the
reason the attempt was cut, moves the deck to the hybrid, runs the diagnostics that separate the
plateau's mechanism from the guard, and routes the physics (1770 K, material tables, five pure-Picard
decks) to Christian. No code is edited before the plan jury and his approval.
**Module:** `src/fem/kernel` (`cl_FEM_Controller.hpp/.cpp`: the stagnation guard, `reset_timestep`'s
floor block, `finalize`); `tmp/dr105_weekend/deck` (T5 deck, not shipped); `todo/debt_register.md`
**AIs involved:** Claude Fable 5.1 (campaign, brief, plan); Codex `gpt-6-astra`/xhigh + Grok
`grok-4.7`/xhigh (blind survey on the brief; plan jury to follow at the same depth)
**Status:** ✅ R1 VERIFIED, R3a ANSWERED 2026-10-03 — the np = 8 pure-Picard run completes 0.3 s on R1 (one floor retry, recovered) and matches the hybrid's T_max to 0.1 K throughout; R2 moot; R3b (tolerance sensitivity), R3c/O3 (tables), O4 open; the serial R1 run continues detached — earlier: R1 VERIFIED BY EXECUTION 2026-10-03 19:55 — `g3_np8_r1_full` met the floor at 26.45 ms with the HEAD abort's six-iterate −43 dB signature, the suspended retry converged on its eighth iterate and was accepted, Δt grew back to the cap three steps later; R3a from the full-length runs pending — earlier: R1 LANDED 2026-10-03 (code jury astra/xhigh + grok-4.7/xhigh: no defect; two post-round fixes applied on Christian's approval; make check 17/17 in the private tree; R0a supports the suspension; the mechanism's own gate — a suspended floor retry on the R1 binary — not yet observed: a cross-binary replay of the HEAD abort is impossible, the deck's trajectory moves with any binary or deck change; `g4_restart_r1`, `g2_serial_r1_full`, `g3_np8_r1_full` running; `make check` on HEAD `b01911f` in Claude's private tree at Christian's request: 17/17) — earlier: PLAN v2 APPROVED 2026-10-03 (Christian: plan, O1 unconditional, O2 leave the five decks, O3 report against the 199 K envelope; R0 launched the same hour; O4 open); implementation is a later, separately approved round — v1's R1 (a patience factor) withdrawn on the plan jury's finding that it delays the cut by at most three iterates; R1 is now the suspension of the guard's Picard verdict at the floor with the reason passed to `reset_timestep`; R0c, R3a–c and two F-sentences corrected on the same round (astra/xhigh + grok-4.7/xhigh, all CONFIRMED). Awaiting Christian's approval and the rulings O1–O4; no source edited
**Literature References:** Messe et al. 2023 §4, Eqs. 12–14 (lagged-conductivity Picard, the
quasi-Newton hand-over, the Eq. 14 relaxation rule; the 10⁻¹¹ final tolerance against
checkerboarding) — not §2.7, which is kernel layout; Walker & Ni 2011 (Anderson type II, as
implemented); Anderson 1965
**Exchange:** `tmp/ai_exchange/review_dr105_picard_stall.md` (pre-registration, two survey entries,
verification and reconciliation; the plan round appends to it). Brief: `tmp/dr105_weekend/STALL_BRIEF.md`

> **Scope guards:**
> - No acceptance criterion changes: a step is accepted only on a certified head residual, as today.
>   `absolute tolerance` is not used as a repair.
> - No physics key of the T5 deck is touched by this plan; residue 3 of DR-105 (the defect, the
>   current table, and now the material tables' validity range) is Christian's ruling (O3).
> - The thermal flat-stall latch is a separate mechanism and is out of scope unless the jury
>   shows it shares the defect.
> - Nothing in `src/homology`.

---

## 1. Findings the plan rests on (survey round, all CONFIRMED by source or log)

- **F1 — the abort.** `magnetic_stagnation_forces_reset` (`cl_FEM_Controller.cpp:877-937`): with the
  5-sample ring full and ε above tolerance, a mean absolute deviation < 0.2 dB declares a stall
  (`.hpp:91-102`). On a pure-Picard deck `try_escalate_to_newton` (`:695-706`) cannot fire, so
  the caller cuts Δt. At the floor (`:2679-2717`) the cut becomes a retry that raises
  `max iterations` and the watchdog window ×2^k (cap 4) — which the guard never consumes, since it
  returns at iterate 6 — and the next attempt starts again at ω = `mOmega0` (`initialize_timestep`, `:170-173`,
  unconditional on the coupled driver; `:2596-2599` on the reset path). The twenty retries
  therefore replay one attempt to printed precision (1 distinct trace among 21 in three of
  the four aborts; 7 in the MUMPS run). `BELFEM_ERROR` at `:2693` ends the run (exit 134).
- **F2 — the attempt that dies.** Predictor at −45…−53 dB (5e-6 … 3e-5, 5–30× the tolerance),
  Anderson-Picard with ω = 1 on the predictor and the first body and the α = 0.5 cut from the
  second body on, the first body lifts the residual by +9…+30 dB, then a flat cluster. The
  cluster's level improves ≈ 6 dB per Δt halving (−29 → −36 → −42 dB on step 929; the last step
  to −43.5 dB is the floor clamp, 0.0013 → 0.0010 ms, not a halving), so the floor is where the
  next halving does not exist; within an
  attempt it does not respond to ω (1.0 vs 0.5). Surviving attempts show plateaus of ≤ 4 samples
  by construction (the guard kills 5), so the plateau's own length and mechanism are **not
  established** — Anderson depth 3 is active throughout and plain Picard was never run.
- **F3 — what is not the cause.** MPI (one abort is serial) and STRUMPACK (one is MUMPS) are not
  necessary conditions. Conditioning is neither established nor excluded: the footer's
  `|λ|max/|λ|min` (8.5e12–1.7e14, no trend with T) is not κ₂ (`:3877-3882`), and the plateau
  (4.5e-5) sits ≈ 3 decades below that ratio times ε. The linear solve's backward error was not
  measured. Thermal enters a stalled attempt at −48 dB and converges in one body.
- **F4 — the hybrid.** `algorithm : Newton` completes 0.3 s at np = 1 and np = 8 (15824 / 15826
  steps, 208 / 211 rejections, 0 floor warnings). The same guard verdict fired 414 times and
  escalated to Newton instead of cutting (234 accepted, 180 cut later and retried at a smaller Δt,
  none reached the floor). Median 3 iterates per step in the quench against 19 for pure Picard.
- **F5 — 1770 K.** No one-step Newton overshoot in either log (largest per-step T_max change
  0.4 K; T_max identical to pure Picard to 0.1 K at every common time to 33 ms; serial and np = 8
  agree to 1e-5). But the user plugin's tables end at 282 K (`cp_HTS`), 300 K (`cp_Cu`), 303 K
  (`cp_Hast`), 390 K (`rho_Hast`), 199 K (`k_Hast`); `cu_cp` / `hast_cp` / `rint_cp` hold the last
  knot above (`matlib.cpp:305-310`), copper ρ is clamped at 1358 K (`:94-99`). Above ≈ 300 K the
  hot spot's heat capacity is frozen: the temperature is the implemented model outside its tables,
  with no evidence of a solver artifact and no claim to validated physics. The deck's 1e-6 tolerance is five decades
  looser than Messe et al. 2023 §4's 10⁻¹¹.
- **F6 — deck facts that matter.** `target iterations : 100` (default 20) pushes Δt up; `hts_jc`
  / `hts_n` return 1.0 at T ≥ 90 K inside the plugin, so the deck's `critical temperature` key
  only guards the library's branch; `J/Jc_max` in the footer is |J| above T_crit.
- **F7 — consumers.** The controller default is Picard (`.hpp:272`); five shipped decks run pure
  Picard (`corc_periodic_bc`, `costheta`, `helix`, `tapestack_circuit`, `undulator2d`), the other
  fifteen the hybrid. The pure-Picard path is live for users.

## 2. Mechanism statement

Established: the run dies because the recovery policy does not match the failure. The guard says
"stuck"; the floor retry answers with more iterations and a longer watchdog, which the guard never
lets the attempt reach, and resets ω, so nothing differs between retries. Not established: why
the Anderson-Picard map sits on a flat cluster for 4+ iterates after one full update from a
near-converged predictor, and whether the cluster would resolve by itself at the floor if the
guard let it (every surviving attempt at larger Δt did resolve, 8–27 iterates).

## 3. Steps

### R0 — diagnostics before the code change (runs only, no source edit; ≈ 3 h of machine time)

- [x] **R0a** — **RESULT 2026-10-03:** serial, guard off, 25 ms, exit 0: 1338 accepted (840 in the quench), 32 rejected — every one by the watchdog at Δt ≈ 0.02 ms (31–38 iterates), none at the floor, zero floor retries. In the quench the steps take a median of 21 iterates (max 58); the longest flat run per step is 4 (median), up to 13; 114 steps carry a flat run ≥ 5 and 28 a run ≥ 8, and all of them resolved to tolerance. So the plateau the guard kills does resolve when left alone, within the iteration budget, and the watchdog is the bound that fires on the genuinely slow steps — at the cap, not at the floor. R1's suspension is sized by this. (`tmp/dr105_weekend/r0a_stall0`) Serial, T5 deck, cold, `stall tolerance : 0` (legal: the band is stored unchecked
  and the test is strict `<`, so 0 never fires), to 25 ms: does the plateau resolve when the guard
  lets it, and within how many iterates? **Record the stopping predicate of every cut attempt**
  (watchdog — conditional: a new minimum refreshes it, a residual within 10× tolerance or a rising
  ω suppresses it, `:806-830` — or the iteration budget, `:1818`, or neither): a watchdog or budget
  stop makes the measured length a lower bound. Outcomes: (i) the plateau resolves inside the
  budget — R1 as designed suffices and the measured length is its evidence; (ii) the run still
  walks to the floor on the watchdog or the budget — the plateau does not resolve at any Δt, R1
  still lands (a stuck attempt should then end on the budget, and the abort is right), and the
  mechanism of the plateau becomes the open row.
- [x] **R0b** — **RESULT 2026-10-03:** serial, Anderson off, 20 ms, exit 0: 1173 accepted, 66 rejected, no floor retry. In the quench (9.5–15.3 ms) plain Picard takes a median of 27 iterates per step (max 63) against 19 (max 47) with Anderson 3; the first-body kick is the same (+23.9 vs +23.6 dB median); the flat cluster is mostly Anderson's: longest flat run median 2, max 5, two steps with a run ≥ 5, against 4, 10 and nine steps with Anderson. So the map produces the kick and a slow descent; Anderson's mixing produces the flat plateau the 5-window guard trips on. (`r0b_anderson0`) Serial, T5 deck, cold, `anderson depth : 0`, to 20 ms: does plain Picard show the
  same first-body kick and cluster? (Codex.) Separates the mixing from the map.
- [x] **R0c** — **RESULT 2026-10-03:** np = 4 cold on HEAD ran 13 ms clean (813 accepted, 62 rejected), so the campaign's np = 4 abort at 11.745 ms was trajectory chance; the serial cold run (deterministic for one binary and deck) is the positive control — but no cross-binary replay exists, see the devlog. (`r0c_np4`) The reproducer: `p1_np4` cold at np = 4 aborted at 11.745 ms in 30 min; rerun it on
  HEAD to see whether it reproduces (np = 8 did not). If not, the serial cold run to 15.4 ms
  (2 h, deterministic) is the positive control. A memdump restart is a legitimate warm start —
  the dump carries the BDF history of both fields and the loader restores it at the saved order
  (`:5147-5169`, `:5370-5375`, `:5436-5478`) — but not a reproducer of the failing attempt: it
  drops the Anderson history (`:5431-5434`), caps the re-entry Δt at `initial timestep`
  (`:5393-5401`), and holds the last accepted save, not the rejected attempt (`belfem.cpp:241-248`).

### R1 — the floor retry responds to the reason the attempt was cut (code, `src/fem/kernel`)

v1 proposed a patience factor on the guard; the plan jury showed it delays the cut by at most
`tFactor − 1 ≤ 3` iterates (overlapping 5-windows, `tFactor` capped at 4) and never reaches the
budgets the retry raises. v2 makes the trade §2 names:

- [x] **R1a** `reset_timestep( const ResetReason aReason )`: an `enum class ResetReason { Stall,
  Watchdog, Budget, Divergence, SolverFailure, NaN, LineSearch, Thermal }` argument at all 20 call
  sites (a missed site is a compile error); stored in `mLastReset` as identical member math on
  every rank, outside the rank-0 print, like the floor factor (`:2683-2684`). `Stall` is the only
  value consumed; the others feed the retry line's text.
- [x] **R1b** The floor block (`:2679-2717`): when `mLastReset == Stall`, set `mStallSuspended =
  true` beside the ×`tFactor` on the budgets and the watchdog windows; the retry line names what
  was raised ("stall guard suspended, iteration budget 800, watchdog 60"). Any other reason leaves
  the guard armed, as today.
- [x] **R1c** `magnetic_stagnation_forces_reset`: while `mStallSuspended` the Picard arm keeps
  its one Newton escalation per attempt and withholds the cut (landed 2026-10-03 this way, not
  as the literal skip v2 wrote — a hybrid deck at the floor keeps its rescue; the code jury noted
  the divergence and the doc matches the code); the Newton→Picard demotion arm is unchanged. The attempt is then bounded by the raised watchdog (conditional) and
  the raised iteration budget (hard, `:1818`), which is exactly the budget the retry exists to grant.
- [x] **R1d** Resets: `mStallSuspended = false` and `mLastReset = None` on the accepted-step path
  where the deck budgets are restored (`finalize`, `:3065`), and `mLastReset = None` at the top of
  every attempt; no per-iterate counter exists, so the six `mResidualHistory.clear()` sites
  (`:204`, `:336`, `:722`, `:924`, `:1048`, `:1413`, and `:2646` in `reset_timestep`) need no
  change. `mStallSuspended` survives `reset_timestep` (it is set there) and dies with the accept.
- [x] **R1e** `nonlinear_controller_theory.md`: the floor-retry paragraph states the reason-aware
  policy. No new deck key; if the jury or Christian asks for one (e.g. `floor stall suspension :
  false`), `doc/input_file_reference.md` and `doc/input_schema.yaml` change in the same session.
- What could break: nothing away from the floor (`mStallSuspended` is false there); at the floor
  after a non-stall cut nothing changes; a genuinely stuck attempt at the floor now spends the
  raised budget (×2, ×4) on each retry instead of 6 iterates — the cost the suspension buys the
  chance of recovery with, bounded by `mMaxFloorRetries`. The thermal flat-stall latch is untouched.
- Rejected on the jury's arithmetic: the patience factor (v1).

### R2 — the T5 deck (deck only)

- [x] **R2** — closed 2026-10-03 as moot: the T5 deck is a `tmp/` artifact kept on pure Picard for the gate runs; the shipped quench decks already run the hybrid. `tmp/dr105_weekend/deck/input.conf`: `algorithm : Newton` (hybrid), as the shipped
  quench decks are configured. Measured: 3 vs 19 iterates per step in the quench, two complete
  runs. This does not fix the pure-Picard path and is not a substitute for R1.

### R3 — the 1770 K (runs and a ruling, no code)

- [x] **R3a** — **RESULT 2026-10-03:** the np = 8 pure-Picard R1 run completed 0.3 s (16733 steps, 544 rejections, one floor retry, 3.35 h) and its T_max equals the hybrid's to 0.1 K at every 25 ms from 10 to 300 ms: Newton did not leave the trajectory; the 1770 K is the model's answer under both algorithms. Field-level differences (T up to 55 K at 147 ms, Jz up to 78 %) are the two 1e-6 trajectories parting, as on corc at 1e-7. The serial R1 run continues (floor met at 11.77 ms, recovered on the first floor retry). With R1 landed: serial pure Picard, cold, to 0.3 s (≈ 13 h) — the check for
  "Newton left the trajectory": at every 1 ms save, max over nodes of |T_Picard − T_hybrid| and
  the same for the thermal field's L2 norm, against the hybrid's own serial-vs-np8 scatter (1e-5
  relative) and against the local temperature rise per save; plus I_1 and U_1. The two runs also
  differ in their Δt histories (adaptation reads the iteration count, `:2844-2868`), so a
  difference has two candidate causes; a disagreement therefore opens a measurement (same Δt
  sequence forced by `adapt timestep : false` at the cap) before any Newton row.
- [ ] **R3b** Tolerance sensitivity on the hybrid: `tolerance : 1e-9` and `1e-11` (the paper's
  terminal criterion, Messe et al. 2023 §4) to 50 ms against 1e-6; the 1e-11 run may sit on a
  linear-solve floor as the corc deck did (DR-192) — that is itself a result.
- [ ] **R3c** Ruling (O3): the plugin's active tables end at 199 K (`k_Hast`), 282 K (`cp_HTS`),
  300 K (`cp_Cu`, `k_HTS`), 303 K (`cp_Hast`), 390 K (`rho_Hast`), with the endpoint held above and
  copper ρ clamped at 1358 K. The validity envelope is the minimum over every active property of
  every active layer — 199 K as the tables stand, reached at ≈ 47 ms. Either every table is
  extended to the temperature the study is to cover (NIST cryogenic and high-temperature
  compilations for Cu and Hastelloy C-276, the DOI beside each number), or the deck's validity is
  stated as that envelope and results beyond it are reported as the model outside its tables.

### R4 — register, rulings, close

- [ ] **R4a** DR-105's row: residue 1 (the rerun) discharged twice; residue 2 restated as F1 (a
  controller recovery-policy defect, not a parallel one; the August wall does not reproduce);
  residue 3 widened by F5 (tables). Decide with Christian whether the row stays or R1 gets its
  own row (prefer extending, per the register's policy).
- [ ] **R4b** Devlog of the campaign + this plan; `todo/README.md`; campaign page.

## 4. Gates

- [x] Positive control: the four HEAD aborts (serial, np 4, np 8, np 8 MUMPS) carry the signature; no cross-binary replay exists (the deck's trajectory moves with any binary or deck change).
- [x] R1 binary: `g3_np8_r1_full` at the floor, step 1699, the suspended retry accepted on its eighth iterate (2026-10-03 19:55); the serial deck to 0.3 s pure Picard (R3a) running; `make check` in the private tree; Christian's
  `make check`; the five pure-Picard shipped decks smoke-run for a few steps each (unchanged
  behavior expected: none reaches the floor).
- Mutation check: with `mStallSuspended` forced false, the reproducer must abort again at the
  same step; with `aReason` forced to `Watchdog` at the stall site, likewise.

## 5. Open questions (Christian)

- **O1** R1's shape — **RESOLVED 2026-10-03 → unconditional suspension, no new key** (Christian).
- **O2** The five pure-Picard shipped decks — **RESOLVED 2026-10-03 → leave as they are**; revisit if one floors (Christian).
- **O3** The material tables — **RESOLVED 2026-10-03 → report the 0.3 s runs against the 199 K envelope (≈ 47 ms), beyond it as the model outside its tables; extending the tables is a separate task** (Christian).
- **O4** `target iterations : 100` in the T5 deck against the default 20 (Grok): retune, or leave
  as the recipe measured it?

## 6. Audit trail

- Survey round (blind, astra/xhigh + grok-4.7/xhigh): `tmp/ai_exchange/review_dr105_picard_stall.md`.
- Plan round: same slug, appended.
- Campaign evidence: `tmp/dr105_weekend/` (STATUS.md, ANALYSIS_PLAN.md with amendments, logs).

## 7. Revisions

- **v2 (2026-10-03, plan jury astra/xhigh + grok-4.7/xhigh):** R1 patience factor withdrawn
  (delays the cut by ≤ 3 iterates, never reaches the raised budgets — both auditors); R1 is the
  guard's Picard verdict suspended at the floor with the reason passed to `reset_timestep`; R0a
  records the stopping predicate; R0c's "BDF history not persisted" corrected (it is persisted;
  the dump drops the Anderson history and caps Δt); F2's relaxation sentence and the ω-reset site
  corrected; R3a given a temperature criterion and the Δt confound; R3b at 1e-11; R3c over every
  active table (`k_Hast` 199 K). Verification: `tmp/ai_exchange/review_dr105_picard_stall.md`.
