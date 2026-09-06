# Devlog 2026-08-24 — penalty input keys land (DR-99 closed)

**Date:** 2026-08-24
**Topic:** The Maxwell IWG's three penalty slots gained an `input.conf`
path: optional sub-blocks `coulomb gauge penalty { chi }` and
`nitsche ghost penalty { eta, k_reg }` under `nonlinear magnetic`,
with deliberately asymmetric absence semantics. Closes DR-99.
**AIs involved:** Claude (design + implementation), Codex + Grok (blind
plan round, blind code round). Design agreed with Christian in-session
(unit-based `k_reg : 1e-3 Ohm ;` was his correction to the draft).
**Claude Confidence:** high
**Codex Audit Confidence:** high (both rounds)
**Grok Audit Confidence:** high ~90% (both rounds; validator source read)
**Verification:** executed — builds clean; six-case parse smoke on a
2D_Tapestack scratch deck (accept/accept/error/convert/error/error, all
exact); 2-rank MPI smoke clean; `make check` 13/14 with the single
failure attributed by both vendors to the pre-existing uncommitted
solver-default change (1e-10 → 1e-9) in `cl_SolverParameters.hpp`;
`check_doc_claims.py` 34/34; schema YAML valid.

## Summary

Parse lives at the end of the `nonlinear magnetic` reads in
`Controller::set_params`. Semantics: gauge block absent = χ 0 (off — the
assembly kernels then skip the G-operator entirely); gauge block present
without `chi` = setup error; ghost block absent = constructor defaults
(4.0, 1e-3 Ohm) stand, because the ghost is load-bearing and must not
switch off by omission. All `set_penalty` calls are collective
(rank-0 set + broadcast) and sit outside the rank guard; only the new
one-line slot log is rank-0. The resolved slots are logged once at setup.

The plan round earned its keep four times over: the original plan's
"include cl_IWG_Maxwell.hpp from the controller" would not have compiled
(kernel CMake has no maxwell include path — both vendors); `get_real` on
the dimensionless keys would have silently SI-scaled `eta : 4 mOhm ;`
into 0.004 (both) — the keys now go through `get_value( key, "-" )` and
reject dimensioned values; `k_reg = 0` would NaN an SC–SC ghost facet
(`km = ρ/h` is exactly zero at `mRhoMin = 0`, so the regularized harmonic
mean evaluates 0/0) — the vendors disagreed here and the source settled
it for the strict `> 0` guard; and my plan facts were a day stale
(Christian had already zeroed the χ slot), which reframed DR-99's actual
debt as "no input path", not "default on".

The code round found only documentation residues (both vendors:
"the parse itself does not need another design round"): a Purpose line
still denying the knobs exist, the schema using `unit:` where the
belfem-conf validator only reads `dimension:` (with `unit:` the checker
would false-pass exactly the deck the C++ rejects), a stale Restriction
paragraph, and the DR-99 reproducer cell. All fixed and gate-verified.
One honest gap is now documented rather than hidden: angle tokens are
scale-only, so `eta : 4 deg ;` slips the dimension check and becomes
4·π/180.

## Changes Made

- src/fem/kernel/cl_FEM_Controller.cpp — parse block + slot log.
- src/fem/maxwell/matrices/mt_maxwell_h.cpp — k_reg comment updated to
  name the deck key (the "expose via psi()" note was about to become
  false); stray `;` removed.
- doc/input_file_reference.md §4.2.1 (new) + doc/input_schema.yaml
  (`sections:` node under `nonlinear`) — the two-artifact rule.
- src/fem/maxwell/doc/ghost_penalty_stabilization.md §3 + Purpose +
  Restriction; src/fem/maxwell/doc/README.md — exposure story updated.
- todo/debt_register.md — DR-99 closed.
- Exchange: tmp/ai_exchange/gauge_input_keys.md (plan, two
  reconciliations) + per-vendor files.

## Open / Recorded (inherited, not this round's)

- belfem-conf `_spec_for` does not resolve section aliases: a deck using
  `nonlinear magnetic` may never get the nested keys schema-checked
  (pre-existing for every key in the section).
- doc/input_file_reference.md still documents the linear-solver default
  tolerance as 1e-10; the header now says 1e-9 — companion edit of
  Christian's uncommitted SolverParameters change, together with
  test_Solver.cpp:301.
- h_ghost insulator-limit comment (`mt_maxwell_h.cpp` ~:437) imprecise.

## Files Updated

- see Changes Made; devlog/README.md index line added; Codex language
  sweep over §4.2.1 and the module-doc rewrite run post-devlog.
