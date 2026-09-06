# Anderson-Accelerated Picard + Controller State Hygiene

**Date:** 2026-07-30
**Purpose:** Stop ts17-class nonlinear stalls at their mechanism: first repair the controller's
relaxation/switching state machine (the E/D/F consensus bundle from the ts17 campaign), then add
opt-in Anderson(m) mixing to the per-kernel Picard branch — replacing blind AIMD ω-adaptation with
a history-based update built on `gels`/`ShiftRegister`.
**Module:** `src/fem/kernel` (`cl_FEM_Controller`, `cl_FEM_DofMgr_SolverData`); consumes `src/linalg/lapack`
**AIs involved:** Claude (plan + implementation), Codex (audit + prose), Grok (third voice)
**Status:** ✅ **CLOSED 2026-08-11 (Christian's ruling: "that is also done").** Moved to
`closed/`. R1–R7 landed and were diff-audited (Codex+Grok); the Anderson bundle is committed
(`6b2a2b98`), the P0 residual defect found in this design was fixed and committed
(`4f2c11cd`), and Anderson is opt-**in** at depth 0/0. R8's end-to-end gate is closed on
Christian's say-so rather than on a recorded A/B; the register keeps the run evidence question
open under `debt_register.md` DR-52 (greg3 A/B), which is where anyone looking for numbers
should go. The unticked boxes below are left as written — they are the original
definition-of-done, not outstanding work.

*Historic status line:***Status:** IMPLEMENTED + COMMITTED + DIFF-AUDITED — R1–R7 landed; Codex+Grok audit round
complete, findings D1–D3 fixed same session (see §4.0); syntax probe clean. BLOCKED on **R8**:
Christian builds, runs `make check` (5 AndersonMixing cases in `tests/fem/test_AndersonMixing.cpp`),
then the tape A/B (deck: drop `update gate`, `target iterations : 20`, A/B `anderson depth : 3`
incl. the depth-0 spot check FIRST).

**2026-08-09 currentness sweep — two corrections to the line above:**
1. **The bundle is COMMITTED, not uncommitted.** `6b2a2b98` (2026-07-30 23:31) carries
   `fn_FEM_anderson_mixing.hpp`, the Controller/SolverData changes and
   `tests/fem/test_AndersonMixing.cpp`. `debt_register.md` DR-10 ("UNCOMMITTED since
   2026-07-30 — land or shelve") and the `devlog/campaigns/controller_anderson.md` header
   are both stale; DR-10 is closed. Anderson also went **opt-OUT** (default depth 3, master
   switch `anderson stabilization`) on 2026-08-07 — that is in tree and committed too
   (`cl_FEM_Controller.cpp:2627`).
2. **R7 is done, not partial** — the five synthetic cases exist and cover exactly what the
   box asks: `empty_history_is_plain_relaxed_step`, `affine_map_converges_in_one_mixed_step`,
   `degenerate_column_shrinks_window`, `window_is_clamped_to_system_size`,
   `unusable_history_falls_back_to_plain_step`.

3. **A P0 defect was found in this design and has since been FIXED — D4, see §4.4.** The
   Anderson path's residual refresh made ε the linear solver's own roundoff. The fix is
   **committed in `4f2c11cd`** (2026-08-07; commit state verified 2026-08-10), and it
   *reverses* this plan's D5
   decision: `mFieldValues` is deliberately no longer refreshed, and Anderson is opt-**IN**
   again (depth 0/0 default). R8 must be run against the fixed semantics, not the ones
   described in §3 below.

Note that CLAUDE.md's `make tests` is a stale target name; the real one is
`make check` / `check-fast`.

> **Scope guards (decided 2026-07-30, Christian):**
> - Phase 1 (R1–R3) + Anderson in ONE batch; per-kernel **Picard branch only** — the Newton
>   branch and the memoryless line search (measure G of the ts17 ranking) are OUT of scope.
> - Activation is **opt-in via input key**; a deck without the key must reproduce current
>   behavior bit-for-bit (depth = 0 hard-bypasses to the untouched legacy Picard block).
> - In fully coupled mode this is **block-Anderson** (each field accelerated while the peer
>   field moves between outer iterations) — accepted openly as standard partitioned practice
>   (O4); a monolithic outer-loop Anderson stays OUT of scope.
> - No persistence: mixing history is per-attempt ephemeral, never saved/restored;
>   `load_memdump` clears it defensively.
> - Validation builds/runs are Christian's (shared build tree; AI does not invoke make).

---

## 1. Current Behaviour and How It Fails

Consensus findings of the 2026-07-30 tri-AI ts17 campaign (devlog
`dl20260730_hphi_ts17_controller_strategy.md`); all verified with citations by Codex and Grok:

| Failure | Mechanism | Evidence |
|---|---|---|
| ω ratchet to floor | promote/demote carries `min(ω_Newton, ω_Picard)` — one-way decay; poisons every Newton entry after a collapse. **Eleven sites, not four** (Grok C2) | `cl_FEM_Controller.cpp:470, 531, 540, 678, 687, 839, 848, 956, 965` (+ related policy at `:1123`) |
| Thermal Picard↔Newton chatter | promotion at `ε₂ < switch` with zero hysteresis; `mJustPicard2` latch unreachable from coupled path | `:671-689`, `:1013-1016` |
| Retry state leaks | `initialize_timestep` resets only the magnetic equation to Picard; `mThermalFrozen` never cleared; `mNumIterationsDiv` shared magnetic/thermal, never reset | `:112`, `:96-122`, `:789-794`, `:923-936`, `:1027-1040` |
| AIMD blind ω-adaptation | grow ≤1.3 / cut ×0.5 on scalar norms only; discards the optimal-ω information in the residual vectors; pins ω at floor under periodic coupled kicks | `:782-800`, `:810-818` |

**Bottom line:** the Picard relaxation machinery has no estimate of the optimal step and a state
machine that leaks across switches and retries; fix the leaks (R1–R3), then give Picard a
principled update (Anderson) instead of retuning constants.

## 2. Architecture: Why Per-Kernel Anderson in the Picard Branch

Type-II Anderson mixing lives where the Picard update already happens: the master-rank per-dof
loop in `SolverData` (`cl_FEM_DofMgr_SolverData.cpp:2206-2224`, inside `is_master()` from
`:2066`), which has `mLhsVector` (= G(x)) and the pre-update compact free-dof values
`mFieldValues` (= x).

> **Residual identity (Grok C1 — load-bearing).** The Anderson residual is the FIXED-POINT
> residual r_k = G(x_k) − x_k = `mLhsVector − mFieldValues`, formed BEFORE the update
> overwrites the dof values. `mRhsVector` after the post-update multiply holds the LINEAR
> residual Ax − b (`:2223-2224`) — that one belongs to the controller's convergence metric
> (`residual()`, `:2345-2350`) and stays untouched. The two must never be conflated.

All required primitives exist in-tree:

- `ShiftRegister< Vector<real> >` — fixed-depth history, owning-type lifetime audited 2026-07-14
  (BDF precedent, `cl_BDF.cpp`). Capacity must be > 0 (`cl_ShiftRegister.hpp:477-480` errors on
  zero), so depth = 0 skips member construction entirely (Grok C8).
- `lapack::gels( A, B, Work, AbortOnError=false )` returning `info` — QR least squares, scratch
  reused across calls (`src/linalg/lapack/fn_gels.hpp:251-317`).
- `lapack::gesvd` — reserved escalation path (truncated-SVD mixing) if the γ-guard proves
  insufficient at n = 35; NOT wired in this round.

Rejected alternatives: monolithic outer-loop Anderson (needs inter-field scaling h vs T —
deferred); normal-equations LS (squares the condition number exactly when ΔR columns go
collinear); NLOPT retuning of α/β/γ (ranked last by all three AIs in the ts17 campaign).

Update (Walker & Ni 2011 type-II; independently re-derived by Grok, confirmed):

```
r_k  = G(x_k) − x_k                      (compact free-dof space)
γ    = argmin ‖ r_k − ΔR_k γ ‖₂          (gels on column-normalized ΔR, m_k ≤ depth)
x_{k+1} = x_k + β r_k − ( ΔX_k + β ΔR_k ) γ
```

depth = 0 (empty γ) reduces exactly to the current relaxed Picard step with β = ω — the
regression anchor. depth = 1 ≈ Aitken.

## 3. Gap Table

| # | State / behaviour | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| 1 | ω carry across Picard↔Newton switches — ALL 11 sites | R1 | min-of-both (defect) | (c) | `:470,531,540,678,687,839,848,956,965`; policy review at `:1123` |
| 2 | Thermal flip-count latch | R2 | latch exists, no trigger in coupled path | (c) | `:691-694` |
| 3 | Retry reset of thermal algorithm / `mThermalFrozen` / `mNumIterationsDiv` | R3 | not reset | (c) | `:96-122` |
| 4 | (x, r) snapshot history per kernel | R4 | absent | (c) | ShiftRegister members in SolverData; capacity > 0 gating |
| 5 | ΔR/ΔX scratch + destructive `gels` copies + Work + γ | R4 | absent | (c) | preallocated members, grown once |
| 6 | Free-dof space for history vectors | R4 | — | (a) | vectors live in COMPACT free-dof space — fixed dofs are absent, not zeroed; `mFieldValues` sized to `mNumberOfFreeDofs`, filled under `!is_fixed()` (`:2015-2029`); Dirichlet enters via `mDirichletMatrix` (`:2072-2103`) (Codex) |
| 7 | γ-guards: `info != 0`, non-finite γ, ‖γ‖∞ > threshold; column normalization | R4 | absent | (c) | soft failure passes `info == 0`; O1 |
| 8 | Accept-commit handshake Controller ↔ SolverData | R4 | absent | (c) | Grok C3 — largest gap; §3.2 |
| 9 | Input key parsing, per kernel (incl. `nonlinear magnetic` alias, `set_params` `:1614-1616`) | R5 | absent | (c) | `nonlinear` / `nonlinear thermal` sections |
| 10 | History flush triggers | R6 | n/a | (c) | inventory §3.1 — extended by Codex |
| 11 | MPI safety of master-only mixing | R4 | — | **verified** | update branch master-only (`:2066`); non-masters run barrier/solve peer branch `:2306-2332`, never touch dofs; residual broadcast `:2342-2379` unchanged (Codex, high) |

### 3.1 Cross-cutting finding: the flush inventory (audited, extended)

Anderson history is valid only within one contiguous run of ACCEPTED Picard iterates of one
attempt of one kernel. Flush (clear history; next step = plain relaxed Picard) on ALL of:

1. attempt/sub-step start and retry: `initialize_timestep()`, `initialize_magnetic()`,
   `initialize_thermal()` (segregated mode uses the latter two — `hphiTrun.cpp:137-166`,
   `cl_FEM_Controller.cpp:218-360`; thermal sub-steps shift fields at `:317-349`),
   `reset_timestep()` (also covers the circuit-failure path `:180-183`, `:291-294`), and
   `reset_thermal()` (`:999-1010` → `:1172-1183`) [Codex],
2. algorithm switch Picard↔Newton at EVERY switch site, including the stagnation latch `:470`
   and the force-Newton escalation `:420-428` [Grok],
3. rejected/backtracked trial — staged pair discarded, never pushed (§3.2),
4. thermal freeze/thaw via the update gate (`mThermalFrozen` transition),
5. `load_memdump()` restore (`:1999-2073`) — defensive, history is never persisted [Codex].

**O4 (peer-kernel policy) — RESOLVED 2026-07-30 (lean, Claude + Codex option b):** in fully
coupled mode the peer field moves between outer iterations, so each kernel's map is
quasi-static, not fixed. We accept this openly as **block-Anderson** (standard partitioned
practice, cf. FSI) rather than flushing on every peer update — the latter would cap the window
at size 0 in coupled mode, making the feature dead there. The γ-guard is the safety net when
peer movement invalidates the history.

### 3.2 Accept-commit handshake (Grok C3 — designed with R4, not after)

Accept/reject lives in the Controller (backup `:504-512`, solve `:587`, accept test `:593-598`,
reject-restore `:601-609`); `SolverData::solve()` cannot know a trial's fate. Therefore:

- `solve()` STAGES the candidate pair (x_k, r_k) in members; it does not push.
- Controller calls `commit_mixing_history()` on accept, `discard_mixing_candidate()` on
  reject/rollback; both are no-ops at depth = 0.
- All §3.1 flush events call `clear_mixing_history()` (staging included).

## 4. Ordered Steps

### 4.0 Implementation Progress (updated 2026-07-30, late session)

**Implemented, pending diff audit + user build (R8):**

- R1: all nine promote/demote/latch carries now copy the clamped ACTIVE-store ω
  (`std::clamp( <other store>, min, max )`, 9 sites); the `reset_timestep` min-carry at the
  divergence branch is deliberately KEPT (conservative carry across a Δt cut is intended there).
- R2: `mThermalFlipCount` counts thermal demotions detected in the `iterate_coupled` handoff;
  the third latches `mJustPicard2` for the attempt.
- R3: `initialize_timestep` now resets the thermal equation to Picard, `mThermalFrozen`,
  `mNumIterationsDiv`, and the flip counter; `mOmega0` documented as internal at its init.
- R4: mixing math lives in the standalone header `src/fem/kernel/fn_FEM_anderson_mixing.hpp`
  (`anderson_mixing_step`, type-II, column-normalized `gels` with logical-column-count call
  into the full scratch allocation — no per-iteration resize; guards: `info>0`/non-finite/
  ‖γ‖∞ > 1e2 → drop-oldest retry once → plain step; `info<0` always fatal). `SolverData`
  stages r = `mLhsVector − mFieldValues` pre-update, gates the legacy block behind
  `mAndersonDepth > 0` (legacy statements untouched for the depth-0 anchor), writes the mixed
  iterate through the same masked dof loop. Commit/discard/clear API per §3.2.
- Flush wiring (merged R6): `initialize_timestep` / `initialize_magnetic` /
  `initialize_thermal` / `reset_timestep` / `reset_thermal` / `load_memdump`, every
  algorithm-switch site (coupled + segregated + escalation + both stagnation latches), gate
  freeze transition, and reject-restore (discard + full flush — with β = live ω the ΔX·γ term
  is not damped by backtracking, so a rejected mixed step must fall back to plain damped
  Picard; rationale from Grok C4).
- O2: commit counters `mAndersonCommits(2)` live in the CONTROLLER (rank-consistent by
  construction — SolverData history is master-only and must not steer control flow); AIMD
  growth branch held while the counter is non-zero, α-decay untouched.
- R5: `anderson depth : 0..8` parsed in both nonlinear sections; thermal forward repeated in
  `set_thermal_kernel` for decks that link the thermal kernel after `set_params`.
- R7: `tests/fem/test_AndersonMixing.cpp` (4 cases: empty-history ≡ plain step, affine
  one-mixed-step convergence, degenerate-column shrink, unusable-history fallback);
  registered in `tests/fem/CMakeLists.txt`. Standalone syntax check of the mixing header
  against real Armadillo + the unified LAPACK layer: 0 errors.

**Known cosmetic debt:** the legacy Picard block sits un-reindented inside the new `else`
(deliberate zero-churn for the depth-0 review); re-indent in a follow-up formatting pass.

**Diff-audit round (2026-07-30 late, Codex + Grok on the uncommitted diff,
`tmp/ai_exchange/anderson_impl_audit.md`):**

- Grok: clean bill on all five focus items — chain indexing (hand-derived 2-column example, no
  oldest/newest swap), shrink drops the oldest pair, O2 reject trace (retry runs with halved β
  AND empty window), workspace threshold = LAPACK documented minimum, all four gtest
  expectations derived to hold (affine case collapses symbolically to x* = c/(1−a), γ = −19).
- **D1 (Codex, HIGH): commit counters advanced on Newton iterates** — nothing is staged
  outside the Picard branch, but the controller helpers incremented unconditionally → AIMD
  growth wrongly held for Newton ω. **Fixed 2026-07-30** (Claude): commit helpers no-op unless
  the live equation algorithm is Picard (rank-consistent guard). Verified by Claude against
  the call sites.
- **D2 (Codex, HIGH): release-unsafe `gels` when window > n_free** — the `tUse ≤ n` assert
  compiles out in release; an underdetermined call overruns the rhs buffer on tiny systems.
  **Fixed 2026-07-30** (Claude): runtime clamp `tUse = min(size, n)` replaces the assert; new
  gtest `window_is_clamped_to_system_size` (n=2, depth=3 → tUsed=2, hand-checked γ=(−1,3)).
- **D3 (Codex, MEDIUM): O3 non-compliance — failed-solve pair was still staged/committable.**
  **Fixed 2026-07-30** (Claude): staging now requires a successful mix OR an empty window (the
  bootstrap pair must stage or the window never fills). Side effect, documented in code: the
  controller commit counter becomes an upper bound of the true window fill — the O2 growth
  hold errs conservative, identically on every rank.
- Codex scope note, accepted as intended: depth-0 bit-identity holds for the SolverData Picard
  block (statement-identical vs HEAD); a FULL run is deliberately not identical to HEAD because
  the R1–R3 hygiene fixes apply to all decks.
- Grok operational caveats (no code change): γ-guard 1e2 needs R8 calibration; clear-on-reject
  intentionally erases accepted commits too (do not "optimize" away without re-auditing O2);
  pre-existing latent issue inherited, NOT introduced: the Picard-branch controller residual is
  evaluated at pre-update `mFieldValues` (Newton was fixed once, Picard never was) — future
  hygiene candidate, out of scope here. *(Superseded by D5 below: it became load-bearing.)*

**First-run findings (Christian's ts16 log, Δt = 0.1 ms, 201-iteration grind, 2026-07-30):**

- **D4 (HIGH, from log): the O2 growth-hold made ω a one-way ratchet.** α-decay and
  backtracking kept shrinking ω, nothing grew it back while the window stayed non-empty →
  ω pinned at 1e-4 for ~130 its, mixed step starved to β ≈ 0 (pure ΔX extrapolation).
  **Fixed 2026-07-30** (Claude): growth branch re-enabled at all four adaptation sites — O2
  amended to "β = live ω, growth ACTIVE, overshoot guarded by flush-on-reject". The
  `mAndersonCommits` counters remain as diagnostics only.
- **D5 (HIGH, from log): the C3 residual lag became load-bearing under Anderson.** The
  controller read each mixed step's quality one iteration late: stale 2e-6 dips triggered
  Newton promotions whose flips cleared the window right after every good extrapolation
  (period-5 cycle in the log), then the true 5e-4 residual demoted back. **Fixed 2026-07-30**
  (Claude): the Anderson write loop refreshes `mFieldValues`, so the post-update multiply
  reports the mixed state's residual (Newton-branch semantics); depth-0 legacy path untouched.
- Watchdog cross-check: this binary predates measure C — with it, the ts16 attempt cuts at
  ~iteration 47 (no new magnetic minimum after it 17) instead of 201.
- Candidate if chatter persists after rebuild: magnetic flip-count latch mirroring R2
  (mJustPicard exists, needs the counter). NOT implemented — expected unnecessary once D5's
  false promotions are gone.

- [x] **R1** — ω-carry fix at ALL ELEVEN sites (`:470,531,540,678,687,839,848,956,965`): copy the
  clamped ACTIVE-store ω across the switch instead of min-of-both; review (not necessarily
  change) the related `reset_timestep` policy at `:1123`. Testable: ts17 log no longer shows
  floor-pinned Newton entries after Picard recovery, on coupled AND segregated decks.
- [x] **R2** — thermal flip-count latch in `iterate_coupled` (after: R1): count promote/demote
  cycles per attempt; on the 3rd, set `mJustPicard2` for the rest of the attempt. Counter resets
  in `initialize_timestep`/`reset_timestep`.
- [x] **R3** — retry hygiene in `initialize_timestep`: `mEquation2->set_algorithm(Picard)`,
  `mThermalFrozen = false`, `mNumIterationsDiv = 0`, R2's flip counter; document
  `mOmega0 = 1.0` as internal (`:43`), do NOT expose as input (ts17 do-not-do list).
- [x] **R4** — Anderson core in `SolverData` Picard branch **together with the §3.2 commit API
  and the §3.1 flush hooks** (after: R1–R3; O2 locked first): members per gap rows 4–7
  (constructed only when depth > 0, C8); r = `mLhsVector − mFieldValues` staged pre-update
  (C1); ΔR column-normalized before `gels` with γ unscaled after (C6); guards per O1; failure
  policy per O3; `info < 0` is a programming error → `BELFEM_ERROR` regardless of soft mode
  (C5). depth = 0 hard-bypasses around the untouched legacy block `:2206-2224` (Codex).
- [x] **R5** — input key `anderson depth : m ;` in `nonlinear` (alias `nonlinear magnetic`) and
  `nonlinear thermal` sections (after: R4); absent/0 = off. Parsed in `Controller::set_params`,
  forwarded per kernel.
- [x] ~~R6 — flush triggers wired separately~~ merged into R4 per Grok's ordering finding
  (splitting R4/R6 lets rejected trials into history in the interim).
- [x] **R7** — synthetic regression: small linear fixed-point test where Anderson(1) must
  converge in ≤2 iterations and depth = 0 reproduces plain Picard exactly.
  *(Done — `tests/fem/test_AndersonMixing.cpp`, 5 cases, committed in `6b2a2b98`.)*
- [ ] **R8** — end-to-end gate (Christian): rebuild, run `tape_hphiTrun` deck (gate removed,
  target 20). FIRST the depth=0 bit-identical check (load-bearing), then A/B with
  `anderson depth : 3` — ts17-class steps must converge in fewer iterations; log γ-norms to
  calibrate O1.
  **Precondition (D4, §4.4):** run this against the FIXED residual semantics. Before the
  2026-08-07 fix the Anderson arm reported ~ −124 dB at the bootstrap iterate by algebraic
  identity and an A/B iteration-count comparison measured nothing. The fix is **committed
  (`4f2c11cd`)** — so R8 is gated only on the run, and on Anderson being explicitly enabled
  (it is opt-in again, default depth 0).

## 5. Open Design Questions

- **O1** — γ-guards: `info != 0`, non-finite γ, and ‖γ‖∞ > 10² (compile-time constant), all →
  shrink window by one and retry once. Threshold is an engineering default, NOT literature-backed
  for n = 35; only meaningful with ΔR column normalization (Grok C6). R8 logs γ-norms to
  recalibrate. *Provisionally set; revisit after R8.*
- **O2** — mixing β — **RESOLVED 2026-07-30 → β = live ω, AIMD growth suppressed (decided
  Christian, per Grok C4):** the Anderson step is scaled by the live line-search ω
  (`tIWG->omega()`), so the controller's backtracking (`:575-578`, `:648`) damps Anderson
  directly; while the window is non-empty the AIMD *growth* branch is skipped (ω held), the
  α-decay branch stays active. Since β IS the live Picard ω store, Grok's carry concern
  (`mOmegaPicard` ≥ β used) is satisfied by construction. Rejected trials additionally flush
  the window per §3.1/3.2. (Codex's fixed-β=1 alternative recorded as rejected for the
  identified backtracking conflict.)
- **O3** — `gels` soft-failure policy — **RESOLVED 2026-07-30 (both auditors):** per-iterate
  fallback (drop oldest column, retry once; else one plain Picard step), bad pair NOT committed;
  no latch-off unless R8 logs show repeated cycling. `info < 0` always hard-errors.
- ~~O4~~ — RESOLVED, see §3.1 (block-Anderson accepted).

## 6. Interface Design

```
solver { nonlinear         { ... ; anderson depth : 3 ; }   // alias section: nonlinear magnetic
         nonlinear thermal { ... ; anderson depth : 1 ; } }
```

- depth m ∈ [0, 8]; 0/absent = off (default). Memory on master, in COMPACT free-dof length
  n_free: (m+1) x-snapshots + (m+1) r-snapshots + ΔR + ΔX (formed on the fly from snapshots)
  + destructive gels copies ≈ (3m+2)…(4m+2) vectors of n_free (Grok C7) + γ (O(m)) + gels Work
  — negligible vs factorization.
- All members preallocated on first activation, sized once, per the no-hot-path-allocation rule;
  none constructed at depth = 0 (ShiftRegister capacity must be > 0).

## 7. Definition-of-Done Checklist

- [ ] Every gap-table row mapped to a step (1→R1, 2→R2, 3→R3, 4–8→R4, 9→R5, 10→R4, 11→verified).
- [ ] depth=0 bit-identical: legacy block untouched, hard bypass, spot-checked FIRST in R8.
- [ ] Flush inventory §3.1 (5 trigger classes) + commit API §3.2 implemented in one batch.
- [ ] O2 signed off by Christian before R4 coding starts.
- [ ] R7 synthetic test passes; R8 A/B run recorded in a devlog with iteration counts + γ-norms.

## 8. Audit Trail

- Exchange threads: `tmp/ai_exchange/hphi_ts17_controller_strategy.md`,
  `tmp/ai_exchange/controller_measures_menu.md` (ts17 campaign, closed),
  `tmp/ai_exchange/anderson_plan_audit.md` (this plan's audit round, 2026-07-30, closed).
- Plan-audit findings, all folded in above: **Grok** — C1 residual identity (`mRhsVector` is
  Ax−b, not G−x), C2 eleven ω-carry sites, C3 accept-commit API, C4 O2 line-search conflict,
  C5 info-sign policy, C6 column scaling, C7 memory count, C8 ShiftRegister capacity-0;
  independent re-derivation of the type-II formula (confirmed). **Codex** — segregated-mode +
  `reset_thermal` + `load_memdump` flush triggers, O4 peer-kernel policy, gap row 6 compact-space
  correction, row 11 MPI verification, depth=0 bypass requirement, input-alias note, prose fixes.
  Findings were re-verified against the cited code before inclusion.
- Decisions 2026-07-30 (Christian): batch scope Phase 1 + Anderson; opt-in input key;
  per-kernel Picard only.


## 4.4 D4 — Anderson residual was the linear solve's roundoff (found and fixed 2026-08-07)

- [x] **D4 (P0, 3/3 independent, source-traced, confirmed by an A/B run, FIXED).** Discovered
  by the greg3 jury round (`devlog/dl20260807_greg3_false_convergence_jury.md`), made concrete
  from F3 of `devlog/dl20260807_newton_extension_jury.md`. **Not a defect in the Anderson
  mathematics — a defect in what the controller then measured.**

  **Mechanism (as it was).** `anderson_update` refreshed `mFieldValues` to the mixed iterate,
  and the residual multiply that followed reused the same lagged matrix that had just been
  solved. On the bootstrap iterate the window is empty and β = ω = 1, so the mixing step
  returns exactly x = A⁻¹b (`fn_FEM_anderson_mixing.hpp:55-56`) and
  ε = ‖A·A⁻¹b − b‖/‖b‖ ≈ 3.6e-13 = −124.43 dB — **independent of nonlinear consistency**.
  The legacy depth-0 path deliberately reports the PRE-update residual and was immune, which
  is why this stayed latent while Anderson was opt-in; the 2026-08-07 opt-out default armed it
  on every deck. Downstream: `run_coupled`/`run_magnetic` exited at `mMinNumIterations = 2`
  with a "converged" ε → one effective lagged-Picard step per timestep → the canonical
  loose-tolerance checkerboard (Messe et al. 2023, paper1, **§4** — not §2.7, which is the
  wrong cite that also sat in the controller comments).

  **Fix, executed the same day on Christian's ruling after a literature consult**
  (Messe 2023 §4 force criterion; Bathe §8.4.4 on increment criteria under-reporting on stiff
  maps; Walker & Ni monitoring ‖G−x‖) — verified in tree 2026-08-09:
  - The `mFieldValues` refresh is **removed**, with the reasoning recorded inline at
    `cl_FEM_DofMgr_SolverData.cpp:2776-2787`. ε is now the pre-update force residual
    ‖A(x_k)x_k − b(x_k)‖/‖b‖, identical to the depth-0 path. Newton keeps its post-update
    recompute (J ≠ A, so it is not degenerate there).
  - **This reverses D5 of this campaign** — the promotion-lag refresh that D5 introduced. The
    accepted cost is at most one extra Picard iterate.
  - `fixed_point_residual()` = ‖G(x_k) − x_k‖/‖x_k‖ is exposed as a **diagnostic only**
    (`cl_FEM_DofManager.cpp:1062`, `cl_FEM_DofMgr_SolverData.hpp:386`); deliberately not a
    convergence gate (Bathe Eq. 8.119 caveat).
  - **Defaults rolled back**, reverting the one-day opt-out experiment in substance:
    `mTimeStepping` BDF5 → BDF1 and `mAndersonDepth`/`2` 3/1 → 0/0
    (`cl_FEM_Controller.hpp:165-166,189`). `anderson stabilization : true` fills 3/1 only
    where no explicit `anderson depth` key exists; an explicit key always wins, including an
    explicit 0. The `scheme` key and `galerkin` parsing stay.
  - Docs synced the same turn (`anderson_acceleration_theory.md`,
    `nonlinear_controller_theory.md`, `doc/input_file_reference.md`).

  **Status of the fix: COMMITTED** — `4f2c11cd` (2026-08-07), with the Picard line-search
  retirement and the PID controller (`1d6ef305`) in the same series; working tree clean,
  verified 2026-08-10 (`debt_register.md` DR-52). The compile half of that gate was already
  retired by DR-55. **The real gate is now the run alone** — R8 here, and the R9 traces of
  `matfix_controller_port_plan.md`.
