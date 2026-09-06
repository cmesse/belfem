# Coupled Newton Hardening: Honest Residuals, Parallel Fixes, B/β Tangent Channel

**Date:** 2026-07-22
**Topic:** hphiTrun coupled-run campaign — lagged-Newton-residual fix, thermal-kernel
parallel crash, PETSc parallel matrix type, and the ρ(|B|,β) tangent wiring
**AIs involved:** Claude (analysis + fixes), Codex (controller analysis + derivation
verification), Grok (coupling numerics, C3 remedy refutation)
**Thread:** `tmp/ai_exchange/maxwell_h_kernels_audit.md`

## Summary

Continuation of the 2026-07-21 campaign. The coupled magneto-thermal solver went
from 110 all-Picard iterations/timestep to sustained dual-Newton convergence in
~35 iterations, serial and 4-proc parallel both working. Root causes closed:
the never-firing Newton promotion (deck's `tolerance switch : 1e-18` sentinel),
the lagged Newton residual reporting, a thermal-kernel parallel aura crash, and
the PETSc parallel matrix type. The remaining magnetic-Newton endgame freeze
(state-dependent, ~−55 dB) was root-caused to the missing ρ(|B|,β) tangent
channel and the fix implemented (Codex-verified derivation).

## Fixes (chronological)

1. **Newton promotion sentinel (deck):** `tolerance switch : 1e-18` in the
   magnetic `nonlinear` block made `mEpsilon < mEpsilonSwitch` unreachable —
   `algorithm : Newton` armed with an impossible trigger. Removed by Christian.
2. **Lagged Newton residual** (`cl_FEM_DofMgr_SolverData.cpp`): Newton reported
   the PRE-update residual (`residual()` reads `mRhsVector`, last written
   before the solve), so every controller decision ran one iteration behind.
   Fix: new `mRhsBackup` member preserves b across the solve; post-update
   recompute `r = A·x_new − b` mirrors the Picard branch. First attempt
   recomputed against the overwritten vector (`A·x − r_old` → residual
   exactly 1.000000) — the b-preservation form is the correct one.
3. **Thermal-kernel parallel aura crash** ("Key 1 not found",
   `cl_FEM_DofMgr_BlockData.cpp`): the facet-neighbor aura block
   (`:370-393`) enrolled master/slave elements without a block-selection
   guard; on the thermal kernel (air not selected) an air aura element hit
   the unguarded `aBlockIndexMap(block_id)` at `:406`. Fix: `key_exists`
   guard, same idiom as the owned path. Secondary find, NOT yet fixed
   (Christian's call): the owned path's `continue` at `:258` skips the slave
   check when the master's block is unselected — silent under-collection for
   BC sidesets facing air.
4. **PETSc parallel matrix type** (`cl_SolverParameters.cpp`): CSR default +
   no auto-selection → parallel PETSc always asserted. Both ctors now default
   to AIJ when `SolverType::PETSc && comm_size() > 1`; explicit
   `matrix format` key still wins. (Enum spelling: `PETSc`, lowercase c.)

5. **Thermal update gate — tried, measured, DISABLED by default.** Christian's
   ts17 diagnosis (thermal updates from a badly unconverged magnetic state
   poison ρ(T) via the garbage j²ρ source) led to a gate skipping the thermal
   solve while `mEpsilon > mThermalUpdateGate`. The rerun refuted it as a
   rescue: ts17 turned into a 180-iteration limit cycle PINNED at the gate
   (magnetic grinds to ~0.98 at ω≈0.03, dips under, thermal updates,
   knockback to 1.2-1.5, repeat) — a state outside the basin needs a timestep
   cut, not update reordering. Default now `BELFEM_REAL_MAX` (off); the
   `update gate` key remains as an experiment knob.
6. **Thermal Picard→Newton handoff ported into `iterate_coupled`** (Christian's
   "call thermal Newton later"): mirrors `iterate_thermal` — thermal starts
   Picard each timestep, promotes only when `mEpsilon2 < mEpsilonSwitch2`
   with `mIteration2 > 1`, demotes on rise, per-algorithm ω store
   (`mOmegaPicard2`/`mOmegaNewton2`, pointer-bound for the tail adaptation),
   plus a Picard restart after any freeze (inert while the gate is off).
   Kills the sticky-Newton failure (ts12 its 14-18; ts16 it.1 thermal Newton
   at residual 1.0). **DECK PREREQUISITE:** `nonlinear thermal { tolerance
   switch : 1e-10 }` must become ~1e-3, else thermal Newton never engages —
   same sentinel trap as the magnetic 1e-18.
7. **PENDING (proposed, awaiting Christian):** sustained-high timestep cut —
   `mEpsilon > ~0.5` for ~10 consecutive iterations past the grace period →
   `reset_timestep()`. This is what actually saves ts17-class steps (cut at
   ~it.10, rerun at Δt/2) — the existing bail-outs are blind to the failure
   mode (divergence cut needs >10, stagnation guard needs a FLAT history;
   the cycle peaks at 1.5 and swings ±2 dB).

## Controller analysis (Codex + Grok, ranked backlog — NOT yet implemented)

- Grok refuted Claude's C3 remedy direction: the ω collapse IS what converges
  a non-contracting 1-1 staggered map; the defect is the atan GROWTH half
  re-arming the instability (sawtooth visible at ts12 its 5-6/9-10) plus the
  1-1 alternation itself. Ranked fixes: oscillation-aware ω growth freeze,
  magnetic sub-cycling, Aitken on the exchanged fields, damping-policy
  symmetrization, thermal handoff port.
- Thermal sticky-Newton risk DEMONSTRATED (ts12 its 14-18: thermal Newton
  walks to residual 0.87 after a magnetic excursion; no Picard restart in the
  coupled path). The `iterate_thermal` handoff (`:891`) needs porting into
  `iterate_coupled`.
- Newton-promotion bounce: first post-promotion step at carried ω ≈ 1
  regresses 20-30 dB; B/β tangent should soften; a first-step ω cap is a
  two-line mitigation.

## The ρ(|B|,β) tangent channel (implemented, Codex-verified derivation)

Endgame evidence: ts7 froze at −55.27 dB in Newton while Picard ground past to
−60; ts12/13 converged in Newton — the freeze is state-dependent, appearing
when the residual concentrates in the tangent-blind ρ(B,β) directions of the
metal/alloy layers (cryogenic Kohler magnetoresistance).

Derivation (exchange thread, Codex-confirmed incl. the previously unmodeled
**C-channel through β's ĵ dependence**):

```
dKdx += ( Cᵀj ) ⊗ Row,   Row = dρ/d|B|·(µ+|h|·dµ/dh)·ĥᵀE + dρ/dβ·dβ/dq
dβ/dq = −sign(c)/√(1−c²)·[ (1/|h|)(ĵ−c·b̂)ᵀE + (1/|j|)(b̂−c·ĵ)ᵀC ]
```

Implementation:
- `Alloy::drhodB/drhodbeta` overrides (table chain rules mirroring Metal's) —
  Codex's completeness blocker.
- MaxwellData: `drhodb`/`drhodbeta` cache slots, `mFundRhodB/Beta` pointers,
  clamp-zeroed wrappers, metal variants with the standard β-memoization
  guard; ctor defaults `return_zero`, metal branches (bulk + TS) wired.
- `add_rho_field_tangent()` helper in `mt_maxwell_h.cpp`, called from
  `h_newton_mu0` and `h_newton_mu`: variable-µ chain (µ + |h|dµ/dh covers all
  µ variants uniformly), β-term 3D-only (2D β ≡ π/2 by convention), kink
  guards |c| ∈ (1e-8, 1−1e-8), field guards 1e-6 (bj_angle precedent),
  preallocated workspaces `vE`/`vC`/`Rw`, no per-point construction, dt owned
  by `assemble_dJdx`.

**Falsifiable prediction:** the state-dependent −55 dB Newton freeze moves
down substantially or vanishes; ts7-class timesteps converge in Newton
without the Picard limp.

## Session-end verification (ts17/ts18 reruns)

- **Thermal handoff verified in-trace:** Picard early, promotion/demotion
  dancing on the threshold exactly as designed, sustained-Newton finales to
  −112/−124 dB at relax 1.0 throughout.
- **Cut-and-retry economics:** ts17 = 5 wasted its + 43-it/15 s converged
  retry at Δt/2; ts18 = genuine state detonation at t≈2.2 ms, Picard-stall
  escalation fired, cut to Δt=0.025, converged 32 its/11 s (vs yesterday's
  180-iteration limit cycle on the same step).
- **B/β tangent verdict: PARTIAL.** ts17-retry finale converges in Newton
  freeze-free; ts18-retry still froze at −58.74 dB (was −55.27) before the
  stagnation guard demoted and Picard finished. Wall moved, not demolished.
  Remaining suspects: non-local hn(φ) coupling (unmodeled by design),
  ∂K/∂T cross-block (out of scope per thermal plan guard), or the drhodb
  path not engaging — one-element FD probe would settle the last.
- Open question: what cut ts17-attempt-1 after 5 its at residual 1.54
  (below every reset threshold Claude traced) — suspect Christian lowered
  the divergence cut to ~1E0; confirm next session.

## Open (next session: "fix the remaining derivatives" — Christian)

- jc/n lookup-table derivatives: `todo/powerlaw_jc_n_field_derivatives.md`
  (R1: JcFunction deval_dB/dbeta/dT; the plan is Codex-polished and the
  consumer machinery from this session is ready for the HTS branch).
- `drho_powerlaw_dT` jc(T)/n(T) legs → the `dbdT = 0` placeholder in
  `compute_drhodT_hts` (formula recorded in the plan's O3).
- The −58.74 dB freeze investigation (FD probe of drhodb first).
- Secondary aura find (`:258` continue-skips-slave) awaiting decision.
- Controller polish backlog: ω-sawtooth growth freeze, first-post-promotion
  ω cap, sub-cycling, Aitken (Grok's ranked list).
- 4-proc coupled verification with all of today's fixes.
