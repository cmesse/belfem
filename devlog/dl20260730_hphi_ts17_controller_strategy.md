# hphiTrun ts17 controller stall — three-AI strategy analysis

**Date:** 2026-07-30
**Purpose:** Diagnose why `tape_hphiTrun` grinds 201 + 62 iterations at timestep 17 (Δt = 0.75 / 0.375 ms) before converging at Δt = 0.1875 ms in 13; derive controller strategy improvements. Read-only session (no source edits); Claude analysis, Codex + Grok independent audits, both with verified `file:line` citations.
**Module:** fem/kernel (`cl_FEM_Controller`)

## Symptom

Coupled magneto-thermal h-φ run. At Δt = 0.75 ms the controller burns all 200 allowed iterations: early phase shows frozen thermal residuals, late phase (~its 70–201) is a period-~12 limit cycle — magnetic Newton pinned at ε ≈ 4e-6 (target 1e-6) with ω collapsed to the 1e-4 floor, thermal cycling Picard→Newton→degrade→Picard across its 1e-4 tolerance switch. The halved attempt detonates when the divergence-counter resets ω to 1.0 (residual 0.13 → 3.89). Only the quartered step converges.

## Root causes (consensus, all citations in `cl_FEM_Controller.cpp` unless noted)

1. **Thermal update gate re-enabled by the input deck.** `update gate : 1e-2` in the run's `input.conf` re-activates the experiment knob that was disabled by default after this exact pathology was observed (header comment `cl_FEM_Controller.hpp:86-93`). While magnetic ε > 1e-2 the thermal solve is skipped (`:660-661`) and stale residuals are printed — including the "Thermal Newton 1 / 3082.55 dB" first line of each retry (sentinel `mEpsilon2 = REAL_MAX` from `:204` plus a stale algorithm label; no thermal solve ran).
2. **Thermal Picard↔Newton chattering, no hysteresis, no coupled-path latch.** Promotion fires the instant ε₂ < switch (`:671-674`), demotion the instant it rises back; the `mJustPicard2` latch is only settable in `iterate_thermal` behind an exact-equality test (`:1013-1016`, `BELFEM_EPSILON` in log10 space) that never fires in practice, and there is no thermal stagnation guard in `iterate_coupled`.
3. **Relaxation carry across the switch is a one-way ratchet.** Promote/demote takes `min(mOmegaNewton2, mOmegaPicard2)` (`:678`, `:687`), so once either store hits the floor, every subsequent Newton entry inherits the floor — each thermal Newton step is doomed to degrade the residual. (Both auditors independently; disproved the "wrong store read" hypothesis — `:697-702` rebinds correctly.)
4. **The magnetic stall guard is effectively dead.** Band 0.001 dB over 5 samples (`cl_FEM_Controller.hpp:67-68`, guard `:445-484`) ≈ 0.023 % residual-ratio flatness; the observed ±0.02–0.1 dB plateau can never trip it, so the step dies only at max iterations. History is also cleared on every algorithm flip (`:549-552`).
5. **Divergence-counter ω reset is a coin flip and its state leaks.** `mNumIterationsDiv ≥ 10 → ω = mOmega0` (`:789-794`) rescued attempt 1 but detonated attempt 2; `mOmega0` is hard-coded 1.0 (`:43`, never parsed from input) and `mNumIterationsDiv` is shared between the magnetic and thermal adaptation paths (`:923-936`, `:1027-1040`) and never reset per timestep.
6. **`target iterations : 100` steers Δt to the edge of the convergence basin.** `adjust_timestep` uses φ = √(target/n) (`:1202-1205`); Messe et al. 2023 (paper1) §4 grows Δt 1.5× only when converged in < k_max/2, and the code default target is 20.
7. **State not reset across retry attempts:** `initialize_timestep` resets only the magnetic equation to Picard (`:112`) — `mEquation2`'s algorithm and `mThermalFrozen` carry into the halved attempt. (`mJustPicard2`/`mEpsilon2`/`mFirstFlip2` are correctly reset at `:105/:204/:213`.)

Verified safe: the magnetic backtracking restore (`:575-651`) completes before the thermal branch (`:653`) and restores the last accepted iterate, so thermal fields need not be in the backup for the current ordering.

## Literature check

Paper1 (Messe et al. 2023) §4, Eq. 12–14: the Picard→Newton handoff follows the paper, but the relaxation rule (α = 0.5, β = 1.1, γ = 0.4, trial-and-error) does **not** match Eq. 14 as printed: the paper's arctan argument (ε_{k+1}−ε_k)/ε_k ∈ (−1, 0] gives a factor in [0.9, 1.1] (ω *shrinks* on fast improvement), while the code (`:785`, thermal `:812`) uses (ε_k−ε_{k+1})/ε_k ∈ (0, 1], factor in (1.1, 1.3] (ω *grows* on fast improvement). The two agree only at stall. Under the ω ≤ 1 clamp the code's sign is the recovery-friendly one (Grok: paper sign is "anti-recovery"; do not flip to paper fidelity without re-fitting β, γ). Which sign is intended is an author question. Deviations: the gate (not in the paper), input tolerance 1e-6 vs the paper's ε_n = 1e-11 (checkerboarding guard — policy choice, flagged), thermal switch 1e-4 vs paper ε_p = 1e-3. The repo library has no coverage of Aitken/Anderson acceleration for partitioned coupling (Bronshtein's Aitken–Neville is unrelated interpolation); those recommendations rest on external FSI literature (Irons & Tuck 1969; Küttler & Wall 2008, doi:10.1007/s00466-008-0255-5).

## Agreed ship order (config now; code changes need approval)

1. **Config:** drop `update gate` from input.conf; `target iterations : 100 → 20`.
2. **Progress watchdog** replacing the stall guard: windowed-slope/projection cut (fire only on: full window, single algorithm, ε > 10× tol, projected iterations exceed remaining budget AND median dB-improvement below threshold, past grace period). Merges the "wider stall band" and "early abandon" proposals while avoiding false cuts on non-monotone coupled histories; also the required safety net once the gate is gone (the +10 dB rules don't catch a sustained 0.5–1.5 spiral).
3. **Thermal flip-count latch** (2–3 promote/demote cycles in one attempt → `mJustPicard2` for the rest) + **fix the ω carry** at `:678/:687` to copy the clamped active-store ω instead of min-of-both.
4. **Retry hygiene in `initialize_timestep`:** `mEquation2->set_algorithm(Picard)`, clear `mThermalFrozen`, reset `mNumIterationsDiv`.
5. **Relaxation strategy (answer to "can α-β-γ be chosen better?"):** retuning the constants is the lowest-leverage knob — the rule's AIMD form discards the optimal-ω information. Newton stage: memoryless line search (restart ω = 1 each iteration, interpolated backtracking; Bathe §8.4) — removes the ω-floor pinning outright. Picard/coupling stage: Aitken Δ² dynamic relaxation on the exchanged coupling fields (parameter-free secant estimate of 1/(1−λ); more realistic than Anderson over full FEM dof vectors). Interim AIMD tweaks if kept: severity-dependent α, remember-last-good-ω re-entry.

## Round 2 — independent measure ranking (same day)

Follow-up prompted by an external (web-API) Claude analysis proposing NLOPT tuning of the Eq. 14 constants. Eleven candidate measures (A–K: config fixes, watchdog, hysteresis/latch, ω-carry fix, retry hygiene, Newton line search, Aitken, Anderson, NLOPT tuning, 1D toy landscape) were put to Codex and Grok as a neutral, unordered menu; Claude's ranking was locked before reading theirs. Result: all three code-aware rankers produced the same top-3 (A: remove `update gate`; B: `target iterations` → ~20; C: slope/projection progress watchdog), the same small-code bundle next (E ω-carry fix — now known to affect four sites, magnetic `:531/:540` included — plus D flip-count latch and F retry hygiene), G afterwards, H before I, and J (NLOPT retune) last or gated behind the structural fixes ("tuning constants around known defects optimizes noise"). The web analysis, which had no code access, independently ranked parameter-free adaptation over constant-retuning — consistent with the rest.

**Eq. 14 sign discrepancy resolved as a probable paper typo:** the code's `(ε_k−ε_{k+1})/ε_k` (grow ω on fast improvement, four sites `:785/:812/:924/:1028`) contradicts the printed equation but matches the paper's own prose ("β > 1 and γ … enable rapid convergence", messe2023.txt:581-584) and is recovery-friendly under the ω ≤ 1 clamp, where the printed sign would shrink ω after good steps. Keep the code; correct the equation in any erratum/follow-up; author confirmation pending.

**NLOPT practicality note (Codex):** `en_Opt_Algorithm.hpp:27-39` wraps only local NLopt methods — a global DIRECT/CRS stage for constant-tuning would need a wrapper extension or external orchestration.

## Attribution

Codex: gate parse trace, `mOmega0`/`mNumIterationsDiv` leak, backtracking-safety proof, early-abandon ranking. Grok: ω min-carry ratchet, F2-vs-F4 causal separation, five false-cut failure modes that shaped the watchdog design. Exchange thread: `tmp/ai_exchange/hphi_ts17_controller_strategy.md` (ephemeral, distilled here).
