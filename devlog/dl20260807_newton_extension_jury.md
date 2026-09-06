# Maxwell Newton-extension convergence regression — three-AI jury round

**Date:** 2026-08-07
**Purpose:** Read-only jury investigation of the convergence regression
reported after the Newton equation for the Maxwell problem was extended
(new material-derivative tangent channels). Hypotheses adjudicated:
H1 = derivative errors, H2 = Anderson + BDF5 strategy, H3 = too-strict
timestep restart criteria (Gregory), H4 = bearing φ-gauge mode (background).
**Module:** fem/maxwell (tangent kernels), fem/kernel (Calculator,
Controller, SolverData), fem/iwg (Timestep), physics/materials (powerlaws,
Kohler metals), numerics BDF.
**Thread:** `tmp/ai_exchange/review_newton_extension_convergence.md`
(pre-registration frozen before dispatch; blind Codex + Grok jury;
verification + reconciliation). Brief:
`tmp/ai_exchange/newton_extension_brief.md`. **Review round was read-only;
F1 + F4 fixes applied afterwards on Christian's go-ahead (see addendum).**

## Confirmed findings (severity, agreement)

1. **F1 — P1, the round's only outright algebra error (Codex, single-raiser,
   Claude-verified by independent derivation):** the piecewise three-regime
   model's flux-flow tangent uses `dt_dc = 1/(a·√(b²+ac))`; the correct
   implicit derivative of a·t² − 2b·t − c = 0 is `1/(2·√(b²+ac))`. Off by
   factor a/2; for a < 0 the tangent's SIGN flips. Present in all eight
   `drho_piecewise_dJ` overloads (`powerlaws.hpp:1315,1400,1476,1547,1623,
   1698,1771,1844`). Newton-only, flux-flow window (j1, j3] only, piecewise
   decks only. At jury time, the fix was one line ×8 and awaited Christian's
   approval (granted and executed same day — see addendum).
2. **F2 — P1 conditional (3/3):** HTS `mFundRhodB`/`mFundRhodBeta` are
   hardwired `return_zero` (`cl_FEM_Calculator.cpp:153-154`) while ρ itself
   evaluates Jc(B,θ,T) via `jc_eval` (`powerlaws.hpp:221`). Inconsistent
   Newton tangent for field-dependent-Jc decks (dρ_PL/djc = −n·ρ_PL/jc,
   n ≈ 30). DOWNGRADED from prime suspect: the todo scope guard records
   that current CORC jc/n are CONSTANT (Christian 2026-07-21), so the
   channel is dead on those decks. Activates the moment a Jc(B,θ) tape law
   is attached (`todo/powerlaw_jc_n_field_derivatives.md`).
3. **F3 — P1 (Claude + Codex):** depth-0 Picard reports the PRE-update
   residual (`cl_FEM_DofMgr_SolverData.cpp:2269-2289` never refreshes
   `mFieldValues`), while Newton (`:2217-2238`) and Anderson (`:2771-2780`)
   report post-update. Documented as intentional legacy, but every
   controller threshold acts on differently-lagged ε depending on strategy;
   A/B convergence comparisons across strategies are systematically skewed.
4. **F4 — P1 latent (3/3):** `compute_drhodT` (and `compute_drhodb`/
   `compute_drhodbeta`) read `mRhoClamped` without refreshing rho
   (`cl_FEM_Calculator.hpp:2570-2580`, `:2867-2888`), unlike
   `compute_drhodj`. Live call orders mitigate (thermal kernel computes rho
   first, `mt_thermal_h.cpp:92`), but the contract is one reorder away from
   a silent wrong tangent.
5. **F6 — P1 hazard, PLAUSIBLE (Claude + Grok, literature tier only):**
   variable-step BDF5 zero-stability requires step ratios near 1. The
   controller applies ×1.5 growth / ×0.5 cuts, while
   `bdf_timestepping_theory.md:110-114` justifies the clamp against the
   BDF2 bound (1+√2), an order-mismatched argument. The published,
   validated strategy (Messe 2023 paper1 §4, Eq. 9-14) is BDF1 +
   Picard/Quasi-Newton solving with the Picard operator A (Eq. 13). The
   consistent-tangent extension and BDF5 are both departures from the
   validated baseline; Anderson likewise has no smoothness guarantee on
   n≈30 power laws (implementation itself verified clean, Walker & Ni
   type-II, hygienic window lifecycle).
6. **F9 — P2, reviewer DISAGREEMENT (Codex vs Grok):** Codex: the +10 dB
   divergence rule is trend-blind after 2 iterations. Grok: the stall band
   (0.001 dB MAD / 5 samples) mathematically cannot trip on healthy
   ≳0.01 dB/it progress — Gregory's "too strict for slow-but-sound" is
   weak as stated; "aggressive once bit-flat" is intentional design.
   Undecidable by review → fixed-Δt gate.

Positive clearances (3/3): variable-step BDF2-5 coefficients are EXACT (two
independent hand derivations); collect_qhist signs are correct;
retry/savepoint bookkeeping round-trips; Anderson mixing is correct + opt-in
(depth 0 default at jury time); dρ/dJ power-law chain rule and the metal
β-chain tangent row algebra are correct; Kohler FD derivative family is
structurally correct. The extension's ALGEBRA is sound; its defects are
omissions plus F1.

**Causal filter (Grok, verified):** with `algorithm : Picard` and Anderson
off, NONE of the extended code executes (`cl_IWG_Maxwell.cpp:318-331` routes
to `h_picard`, which computes only mu/rho). A pure-Picard regression cannot
be caused by the Newton extension — it must come from BDF order/Δt policy,
residual semantics, deck/material changes, or mis-attribution. The bearing
gauge mode (probe-tier evidence, dl20260806) remains the best-supported
Newton-flat driver on net-current decks; F1/F2 are additional, separable
Newton-quality defects.

## Executable gates proposed (pending Christian)

- **G-N1:** piecewise-deck Newton A/B with the F1 one-line fix.
- **G-B1:** `method : BDF1|BDF2` vs `BDF5`, controller untouched — decides F6.
- **G-P1:** `adapt timestep : false` fixed-Δt run — decides F9.
- Deck audit: does any failing deck use the piecewise model (F1) or a
  JcFunction with B/θ dependence (F2)?

## Attribution

Codex: the F1 flux-flow derivative error (round's best find), residual-lag
P1 framing, trend-blind divergence rule. Grok: the architecture filter D6,
C5 impact downgrade via the todo scope guard, H3 refutation math, bearing
V7 re-rank. Claude: pre-registration (C1-C11), BDF hand-derivations,
literature addendum (Messe 2023 §4 baseline, BDF2-bound mismatch in the
theory doc, Anderson smoothness caveat). Gregory: the field report and the
restart-criteria hypothesis. Christian: investigation direction (material
functions + Newton matrices; literature check on Anderson/BDF5).

## Addendum (same day): fixes executed + deck audit

On Christian's go-ahead: **F1 fixed** — `dt_dc = 0.5/√(b²+ac)` in all eight
`drho_piecewise_dJ` overloads (powerlaws.hpp); **F4 fixed** —
`compute_drhodT`/`compute_drhodb`/`compute_drhodbeta` now refresh rho before
reading `mRhoClamped` (cl_FEM_Calculator.hpp). Working tree only, not
compiled. F2/F3/F5/F6 deliberately not touched (separate decisions).

Deck audit: `greg`/`greg2` = Picard, BDF1, constant jc/n, so F1 and F2 are
dormant there; Gregory's regression needs a different driver. `corc` =
Newton + **bdf5** + adaptive Δt + metal layers, so it runs the full suspect
stack; G-B1 (BDF1/2 vs BDF5 A/B) is the decisive next run. (Chronology
note: the audit reflects the defaults at jury time; the same-day opt-out
change — `dl20260807_timestep_input_optout.md` — means decks without a
scheme key, greg/greg2 included, run BDF5 + Anderson from now on.)
