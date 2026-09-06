# Devlog 2026-08-26 — DR-72: Increment-Form Picard Solve Landed

**Date:** 2026-08-26
**Topic:** DR-72 remedy A′ — the Picard branch of `SolverData::solve` now solves
A·δ = A·x − b instead of A·x_new = b
**AIs involved:** Claude (implementation), Codex + Grok (plan audit and code audit, both
rounds, both approve)
**Claude Confidence:** high on the algebra and branch mechanics; the executable gates are owed
**Codex Audit Confidence:** high (plan and code rounds)
**Grok Audit Confidence:** high (~95%) on the constraint scorecard; open risks carried below
**Literature References:** Bathe 2016 §8.4.1 (incremental-iterative scheme, p. 755), §8.4.4
Eq. 8.120 (increment-normalized criterion — cited as motivation, not implemented); Messe et
al. 2023 §4 Eq. 10-11 (pre-update residual criterion, preserved)
**Verification:** VERIFIED BY EXECUTION 2026-08-27 (see §Gate Results below) — adiabatic A/B
on a dr30_bulk-derived deck in an own build tree (`cmake-build-claude`, Debug, np=2):
absolute form drifts 4.4e-2 K in 10 steps at rtol 1e-6, increment form holds 1.4e-9 K
(3e7x); `make check` 14/14, `check-fast` 9/9. Production-deck parity still owed.

## Summary

DR-72's re-scoped remedy is implemented: the Picard branch computes the residual
r = A·x_k − b with a master-local SpMV *before* the solve, backs it up, solves for the
increment δ, updates `value -= ω·δ` (exactly the relaxed Picard update in exact
arithmetic), stages Anderson with R = −δ, and restores r for `residual()`. An iterative
linear solver's relative tolerance is now measured against the increment scale instead of
the absolute 77 K / 300 K field offset — the scale-free criterion DR-70's mitigation
(default rtol 1e-10) only approximated. Under an inexact Krylov solve the increment form is
*tighter* than the absolute form (lagged residual ≲ rtol·ε_nl instead of ≲ rtol at ω = 1);
that asymmetry is the feature, and parity gates must judge physics, not bitwise fields.

## Key Findings

- Only the Picard branch solved the absolute system; Newton was already in increment form
  and served as the in-file template (`cl_FEM_DofMgr_SolverData.cpp`, both arms of the same
  switch).
- Reset safety (Christian's explicit concern) traced clean twice, independently by Claude
  and Codex: every controller solve site checks `solve_failed()` before `residual()`;
  `reset_timestep` (`cl_FEM_Controller.cpp:2164+`) restores fields/savepoints and clears
  Anderson history but never reads `mRhsVector`/`mLhsVector`/`mRhsBackup`; the next
  assembly re-derives all three (`reset_rhs_vector`, `update_field_values`, the solve
  itself). An abandoned iterate leaves no state a consumer reaches.
- The NaN path (`std::isnan(mRhsNorm)`) is unchanged in effect: `mRhsNorm` is still ‖b‖ and
  is the NaN carrier; the controller's NaN guards fall through to `reset_timestep` as
  before (Grok, high ~90%).
- No current vector solver wrapper modifies the caller's RHS in place (Codex checked
  SuperLU, MUMPS, PETSc, STRUMPACK individually); the backup/restore is conservative and
  keeps the branch robust.
- The residual-reporting contract is untouched: `residual()` still returns
  ‖A·x_k − b‖ / ‖b‖ at the pre-update iterate ( `mFieldValues` deliberately stale — the
  D5-reversal rule stands, though the mechanism is now backup-backed rather than
  recompute-backed ).
- `mLhsVector` between solves now holds δ, not the absolute solution. Consumers: the
  opt-in `initial guess : true` deck flag (default false; post-change the stale guess is a
  previous δ — the right scale, except across a Picard→Newton promotion) and the debug
  HDF5 `LHS` dataset (now δ on Picard decks, as on Newton decks since forever).

## Changes Made

- `src/fem/kernel/cl_FEM_DofMgr_SolverData.cpp` — Picard branch restructured to increment
  form; Anderson staging `R = −mLhsVector`; `anderson_update` head and body comments
  updated (the body comment now warns that re-inserting a post-update multiply would
  compute A·x − r = b and pin ε at 1).
- `src/fem/kernel/cl_FEM_DofMgr_SolverData.hpp` — `mRhsBackup` comment broadened (Newton
  preserves b, Picard preserves r, exclusive per call). NOTE: this hunk was swept into
  commit `c65ea8d3` ("added bscco-2223 dataset") by a concurrent session; the `.cpp` change
  is the uncommitted remainder.
- `src/fem/kernel/doc/anderson_acceleration_theory.md` — two sentences updated (depth-0
  form; staged pair reconstructs G = x − δ), per Grok's code-audit finding.
- `todo/debt_register.md` — DR-72 status column updated same-session (implemented,
  reviewed-not-verified, gates listed).

## Gate Results (2026-08-27, own build tree, Christian-authorized, ≤2 procs)

Build: fresh `cmake-build-claude/` (GNU 11, Open MPI, MKL/Blaze, `USE_DEBUG=ON`,
`USE_TEST=ON`), `make -j2`, clean. Suite: **`make check` 14/14, `check-fast` 9/9.**

Adiabatic drift A/B — dr30_bulk-derived tapestack3d deck, adiabatic (no thermal BC
section → empty BC factory), sigmoid drive cut to ~0.2 mA, PETSc/ASM thermal at
deliberately loose `relative tolerance : 1e-6`, nonlinear thermal 1e-3 (below nothing —
above the 1e-4 switch, so every iterate stays Picard), 10 fixed 5 ms steps, cold start,
`mpirun -np 2`. A = same tree relinked with the HEAD (absolute-form) TU compiled
out-of-band (shared source tree never touched); B = increment form. Binary identity
proven by archive-member byte-compare (12,003,680 B = increment object; absolute
compiles to 12,003,272 B).

| run | max\|T−77\| at t=55 ms | growth | thermal residual floor |
|---|---|---|---|
| A (absolute) | **4.373e-2 K** | linear, ~3.7e-3 K/step ≈ 48×rtol·T, all 257,681 nodes off | −49…−52 dB (the linear exit test) |
| B (increment) | **1.364e-9 K** | non-accumulating | **−108.7 dB = 1e-11** |

Three signatures beyond the headline 3×10⁷ drift reduction:

- **Magnetic parity is exact to the printed digit** (−20.38 dB iterate 1, −156.54 dB
  iterate 2, STRUMPACK direct — identical A vs B): behind a direct solve the increment
  form is roundoff-identical, as the audits predicted.
- **B reaches 1e-11 through an ITERATIVE solver at rtol 1e-6** — the Messe et al. 2023
  nonlinear target the headroom rule said only a direct linear solve could reach. The
  increment form's Krylov error scales with the increment, so the printed nonlinear
  residual is no longer floored at the linear tolerance.
  `Controller::check_iterative_solver_headroom` is therefore over-conservative for the
  Picard-increment path (it aborted the first gate deck at nonlinear 1e-7 / rtol 1e-6,
  a combination the new code demonstrably satisfies) — recorded, not changed; relaxing
  it needs its own round.
- **Zero soft-fails, DIVERGED, or Δt resets in either log** — Grok's PETSc noise-floor
  risk did not fire in 10 settled steps; a long fully-settled run remains the fuller test.

By-catch, registered as **DR-111**: the first 2D gate attempt (2D_Tapestack + thermal
sections) aborts in the FIRST thermal assembly on both ranks — `Assertion adV >= 0.0
failed` (Negative Jacobian determinant), `Calculator::dV` → `T_h_picard` →
`IWG_MaxwellThermal::compute_mkf`. That was the first-ever smoke of the 2D thin-shell
thermal path (EF_QUAD4TS fix 4a42d982 smoke runs were pending); the 3D path is clean
under the same debug guard. A gdb probe measured the failing determinant:
`-5.000000000022252e-12` — exactly the healthy wafer Jacobian magnitude
(1 µm/2 × ~20 µm/2), sign flipped: the winding is reversed, the geometry intact.
Root cause: `create_elements_on_blocks_line2` assumes a right-hand pairing between
the LINE2 facet direction and the layer-offset normal that nothing enforces, and the
layer elements are born after the MeshChecker pass, outside its scope (full trace in
the DR-111 row). Debug-only guard + negative mass in release = the DR-06
cold-spot mechanism, now with a reproducer (`cmake-build-claude/gate_dr72/`). Also
tripped en route: `check_iterative_solver_headroom` (working as designed).

## Open Questions / Owed Gates

- ~~Adiabatic mini-deck A/B at loose rtol~~ — done, see §Gate Results (4.4e-2 → 1.4e-9 K).
- ~~`make check` / `check-fast`~~ — done, 14/14 and 9/9.
- Forced-thermal + magnetic parity on a PRODUCTION deck (real drive, real Joule): still
  owed. The gate deck gives magnetic parity at printed precision and thermal physics
  near-zero by construction; a forced run is the remaining evidence.
- PETSc noise-floor soft-fail on a long fully-settled adiabatic run (Grok ~55%): did not
  fire in 10 steps; unresolved for 1e4-step scale.
- DR-111 diagnosis (2D TS thermal negative Jacobian) — separate row, reproducer preserved.
- Exchange thread `tmp/ai_exchange/dr72_increment_picard.md` distilled here; GC-eligible.

## Files Updated

- src/fem/kernel/cl_FEM_DofMgr_SolverData.cpp
- src/fem/kernel/cl_FEM_DofMgr_SolverData.hpp
- src/fem/kernel/doc/anderson_acceleration_theory.md
- todo/debt_register.md
- devlog/README.md
