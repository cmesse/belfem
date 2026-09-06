# Devlog 2026-08-28 — DR-97 latch, DR-119 coupled clamp, DR-120 parse validation

**Date:** 2026-08-28
**Topic:** Three controller-guard fixes landed in one jury round: the thermal-budget
magnetic-hit-target latch (DR-97 extension), the coupled-mode thermal watchdog
window clamp (DR-119), and setup validation of the nonlinear iteration keys plus
siblings (DR-120)
**AIs involved:** Claude (plan + implementation), Codex + Grok (plan audit AND code
audit, both blind), Codex (language sweep)
**Claude Confidence:** high
**Codex Audit Confidence:** high (code round: approved, zero in-scope defects)
**Grok Audit Confidence:** high (code round: "landed as reconciled", zero defects)
**Verification:** compile/link + focused regression — Christian built the tree and
ran the test suite on 2026-08-29, all green, on top of the session's
`g++ -fsyntax-only` probe. The three behaviour-specific run gates (DR-97 A/B,
DR-119 trip with ε2 > 10× tol and non-growing ω2, DR-120 negative-deck setup
abort) are NOT covered by the suite and survive in the register rows.

## Summary

Christian asked to fix DR-97, DR-119 and DR-120 together and ruled the three scope
questions up front: extend DR-97's conservative budget predicate with a per-attempt
latch, fix DR-119 with a coupled-mode clamp (not a default change, not a warning),
and widen DR-120 to the sibling wrap-prone keys. Full plan+audit→code+audit with
both vendors; every accepted tightening applied. All changes in
`src/fem/kernel/cl_FEM_Controller.cpp` / `.hpp`.

## Key Findings

- The step-267 gap of the landed DR-97 predicate is real and closable: `run_coupled`
  requires both fields under target simultaneously, so latching "magnetic reached
  target at least once this attempt" cannot accept a degraded magnetic iterate — it
  only stops paying for thermal updates that keep spoiling a converged one.
- Grok's plan audit sharpened the DR-119 derivation: the last iterate that can run
  `watchdog_thermal` is the magnetic ceiling itself (`tMagneticDoomed` skips the
  thermal block on the ceiling iterate), so the bound
  `magnetic max − thermal min − 1` trips at exactly the last capable iterate in the
  measured case. Both vendors demoted the clamp from structural guarantee to
  conservative policy sized to an early thermal best.
- Two defects found beyond the DR-120 row's text, both landed: the magnetic
  `stall window` clamp never saw negatives (the `std::max<uint>` promoted the
  negative before clamping — found by Claude in the pre-plan read), and
  `coupling factor` unconditionally divides `initial timestep` (found by Codex;
  Grok had argued the key was dead — resolved against the tree in Codex's favor:
  segregated `0` was a division by zero).
- A magnetic-only deck defaults `mIsFullyCoupled = true`, so the clamp's setup
  notice had to be gated on thermal-section presence (Grok C3), not on the flag.

## Changes Made

- `cl_FEM_Controller.hpp`: `bool mMagneticHitTarget` beside `mThermalStalled`.
- `cl_FEM_Controller.cpp`: latch set + predicate substitution in
  `iterate_coupled()`; latch resets in `initialize_timestep()` and
  `reset_timestep()`; validated parses for `max iterations` / `min iterations`
  (both blocks, with a `max ≥ min` cross-check that also fires against the other
  key's default), `watchdog window` (both blocks, ≥ 0, 0 stays the disable),
  `stall window` (≥ 0 before the cast, 0/1 still clamp to 2), `coupling factor`
  (> 0); coupled-mode clamp of `mWatchdogWindow2Deck`/`mWatchdogWindow2` after the
  deck snapshot with a rank-0, thermal-section-gated `message()` notice.
- Input contract, both artifacts: `doc/input_schema.yaml` (constraints on all
  seven touched keys; latch note on thermal `max iterations`; clamp note on
  thermal `watchdog window`; stall-window constraint corrected from "≥ 2" to
  "≥ 0 + clamp"), `doc/input_file_reference.md` §4.2/§4.3.
- `src/fem/kernel/doc/nonlinear_controller_theory.md`: divergence-rules row,
  progress-watchdog row and parameter table updated for latch + clamp (the theory
  doc still carried the pre-latch "at target on the same iterate" claim — Grok C5).
- Codex language sweep applied to the touched reference/theory prose; two
  load-bearing technical sentences the sweep's replacements dropped were kept.
- Register: DR-97/DR-119/DR-120 status cells updated (landed, reviewed not
  verified, gates surviving); DR-137 filed (magnetic watchdog vs its own ceiling,
  Grok's flag, deliberately not folded into DR-119).

## Open Questions

- The three run gates: DR-97 A/B (coupled deck, thermal budget < magnetic budget);
  DR-119 gate must show ε2 > 10× thermal tolerance and non-growing ω2 at the trip
  iterate at InfoLevel::Default or louder; DR-120 negative-deck setup abort.
- DR-137 (new): whether the magnetic watchdog deserves the same clamp.
- The register's separate coupled line-search hypothesis (thermal residual
  re-measured against a moving magnetic solution) remains untested and possibly
  more important than the wall-clock cost this session recovered.

## Files Updated

- src/fem/kernel/cl_FEM_Controller.cpp / .hpp
- doc/input_schema.yaml, doc/input_file_reference.md
- src/fem/kernel/doc/nonlinear_controller_theory.md
- todo/debt_register.md (three rows amended, DR-137 filed)
