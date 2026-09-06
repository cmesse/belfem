# Anderson-accelerated Picard + controller hygiene — implementation session

**Date:** 2026-07-30 (evening; follows the ts17 analysis session, same day)
**Purpose:** Implement plan `todo/anderson_picard_acceleration_plan.md` R1–R7: controller
relaxation/switching state-machine repairs plus opt-in per-kernel Anderson(m) mixing of the
Picard branch. Decisions by Christian: batch scope Phase 1 + Anderson, opt-in input key,
Picard-only, O2 = β = live ω with AIMD growth suppressed.
**Module:** fem/kernel; consumes the new unified `src/linalg/lapack` layer (parallel session by a
second Claude instance; `gels` gained the `AbortOnError`/`info` pattern on request from this one).

## What landed (uncommitted; ~440 insertions, 5 files + 2 new)

- **`fn_FEM_anderson_mixing.hpp` (new):** standalone type-II Anderson step (Walker & Ni).
  Fixed-point residual r = G(x)−x; difference columns from the committed ShiftRegister history
  chained with the current pair; column-normalized QR least squares via the typed `lapack::gels`
  backend called with the *logical* column count against the full scratch allocation (leading
  dimension carries the stride — no per-iteration resize, valid for Armadillo and Blaze
  spacing); γ un-scaled after the solve. Guards: `info > 0`, non-finite γ, or ‖γ‖∞ > 1e2 →
  drop the oldest column, retry once, else plain relaxed step; `info < 0` always fatal.
- **`SolverData`:** members + `set_anderson_depth` / `anderson_commit` / `anderson_discard` /
  `anderson_clear`; `anderson_update` stages (x, r) from `mLhsVector`/`mFieldValues` *before*
  the update overwrites the dofs (the fixed-point residual, NOT the linear `mRhsVector` —
  Grok's C1 from the plan audit), then writes the mixed iterate through the same masked dof
  loop. The legacy Picard block is statement-identical inside a `depth == 0` bypass (deliberate
  un-reindent for the bit-identity review; formatting follow-up noted in the plan). History
  registers are heap-allocated only while depth > 0 (ShiftRegister rejects capacity 0).
- **Controller:** all nine ω promote/demote/latch carries now copy the clamped ACTIVE-store ω
  (the min-of-both ratchet is gone; the `reset_timestep` min-carry stays by design); thermal
  flip-count latch (3rd demotion in an attempt pins `mJustPicard2`); retry hygiene in
  `initialize_timestep` (thermal → Picard, `mThermalFrozen`, `mNumIterationsDiv`, flip counter);
  commit-on-accept / discard+flush-on-reject around the backtracking loop; flush at every
  trigger of plan §3.1 (attempt/sub-step inits, resets, every algorithm-switch site, gate
  freeze, `load_memdump`); AIMD *growth* held while the Anderson window is non-empty (O2),
  tracked by controller-side commit counters so every rank branches identically; input keys
  `anderson depth : 0..8` in both nonlinear sections, thermal forward repeated in
  `set_thermal_kernel`.
- **Test (`tests/fem/test_AndersonMixing.cpp`):** empty history ≡ plain relaxed step;
  affine uniform-contraction map converges in ONE mixed step (hand-checked: β = 0.5, γ = −19,
  x₂ = 10c = x* exactly); degenerate oldest column → window shrink; unusable history → plain
  fallback.

## Design points worth remembering

- **Reject ⇒ full flush, not just discard.** With β = live ω, backtracking damps the β r and
  β ΔR γ terms but NOT the ΔX γ term of a mixed step — halving ω cannot rescue a bad
  extrapolation. Clearing the window on reject makes the retry a plain damped Picard step,
  which the line search *can* control. (Grok C4 rationale, extended.)
- **Commit counters live in the Controller, not SolverData.** The history is master-only; any
  controller branch keyed on it would diverge across MPI ranks. The controller counts its own
  commit/clear calls — identical code path on every rank.
- **Verification state:** standalone syntax check of the mixing header against real Armadillo +
  the unified LAPACK layer: 0 errors (re-run after the audit fixes). NOT yet compiled into the
  project, NOT yet run — R8 (build, `make tests`, tape A/B with `anderson depth : 3`, γ-norm
  logging) is Christian's morning gate.

## Diff-audit round (same session)

Grok independently re-derived the implementation from the code alone: clean bill on chain
indexing, shrink-oldest semantics, the O2 reject trace (retry = halved β + empty window),
workspace floor, and all five gtest expectations. Codex confirmed depth-0 statement identity
of the legacy Picard block and complete §3.1 flush coverage, and found three real defects,
all fixed in-session: (D1) commit counters advanced during Newton iterates → helpers now
guard on the live algorithm; (D2) the `tUse ≤ n` assert was release-unsafe on tiny systems →
runtime clamp + new `window_is_clamped_to_system_size` gtest; (D3) a failed-solve pair was
still staged → staging now requires a successful mix or the empty-window bootstrap. Inherited
(pre-existing, untouched): the Picard-branch controller residual is evaluated at pre-update
`mFieldValues` — the Newton branch was fixed for this once, Picard never was; future hygiene
candidate.
