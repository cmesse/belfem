# Devlog 2026-08-27 — DR-97: coupled-mode thermal iteration budget enforced

**Date:** 2026-08-27
**Topic:** `nonlinear thermal { max iterations }` was silently inert in fully-coupled mode; landed the conservative enforcement plus the full doc half
**AIs involved:** Claude (primary), Codex (plan + code audit, prose sweep), Grok (plan + code audit)
**Claude Confidence:** high on the landed logic (three concurring static traces); the regression question is deferred to the A/B by design
**Codex Audit Confidence:** high (code approved, two doc findings)
**Grok Audit Confidence:** high (code approved, three doc findings; refuted the original predicate in the plan round)
**Literature References:** none required — controller policy, not a formulation change (protocol §5)
**Verification:** reviewed, NOT verified — `g++ -fsyntax-only` with the kernel build-tree flags ran green on the edited TU (compile-level only, no link/run); the A/B regression gate named in the register remains owed

## Summary

DR-97's code half is landed. `iterate_coupled()` now consults the thermal iteration budget
through a new `tThermalBudgetSpent` clause in `tThermalReset` (`cl_FEM_Controller.cpp`, search
for the identifier): the budget counts completed thermal solves and cuts the attempt only when
it is spent while the thermal field is the sole keep-alive — magnetic at its relative-or-absolute
target, thermal at neither, not stall-latched, not mid-flat-streak, and only on an iterate that
actually solved thermal (`tUpdateThermal` guard). The doc half closed in the same session:
`doc/input_file_reference.md` §4.3 gained a `max iterations` row, `doc/input_schema.yaml` gained
`default: 100` plus a scope note on the thermal key, and
`src/fem/kernel/doc/nonlinear_controller_theory.md` §1/§4/§6 were corrected (that file carried
the same "the ceiling works" overclaim as the input contract).

## Key Findings

- The plan round split the vendors on the load-bearing point. The pre-registered predicate
  (reject whenever the budget is spent and thermal is unconverged) was refuted by Grok: with the
  update gate off (the default), the thermal solves in lockstep with the magnetic iterates, so
  `mIteration2` tracks the coupled iterate count and an unconditional thermal ceiling caps the
  whole coupled loop at the thermal budget while both fields are still converging. Codex had
  approved the same predicate with the A/B as the gate. Christian chose the conservative form.
- Honest cost of that choice, against the two measured tapestack3d grinds: the conservative
  clause catches the step-276 shape (magnetic under target at the iterate-16/17 check while
  thermal alone kept the loop open) and may miss step 267 if the magnetic sat re-perturbed above
  target on every post-budget check. Part, not all, of the measured ~23/~44 min per grind is
  recovered.
- Both code audits found no logic defect in the landed clause. Their doc findings (watchdog wins
  a same-iterate tie with the ceiling; a flat streak defers the cut rather than accepting the
  step; "lockstep" needs the update-gate qualifier) were applied to the code comment, both input
  artifacts, and the theory doc.
- Grok independently closed the plan round's Newton-escalation worry: the iterate that newly
  sets `tMagneticReset` skips the thermal solve via `tMagneticDoomed`, so the budget clause
  cannot newly coincide with a magnetic reset and never newly suppresses
  `try_escalate_to_newton()`.
- Two spin-off rows filed: DR-119 (the coupled thermal watchdog cannot fire before the magnetic
  ceiling when the magnetic budget ≤ the thermal watchdog window) and DR-120 (both nonlinear
  `max iterations` keys parse into `uint` unvalidated; a negative value silently disables the
  ceiling).

## Changes Made

- `src/fem/kernel/cl_FEM_Controller.cpp` — the `tThermalBudgetSpent` clause and comment
- `doc/input_file_reference.md` §4.3 — new `max iterations` row (Codex-swept)
- `doc/input_schema.yaml` — thermal `max iterations`: `default: 100` + scope note
- `src/fem/kernel/doc/nonlinear_controller_theory.md` — §1 convergence paragraph, §4 overview
  and divergence-rules row, §6 key table (Codex-swept)
- `todo/debt_register.md` — DR-97 settled and retagged `[F]`→`[P]` (A/B gate owed), DR-119 and
  DR-120 added, freeze-lens preamble updated

## Open Questions

- The A/B regression gate: decks with a small thermal budget that currently grind to convergence
  past it while the magnetic is at target would now be rejected; include coincident
  magnetic-escalation cases. Until it runs, the change is reviewed, not verified.
- The DR-97 row's separate untested hypothesis (thermal residual degrading at the relaxation
  floor as a coupled line-search interaction) is untouched by this session and remains open.

## Files Updated

- src/fem/kernel/cl_FEM_Controller.cpp
- doc/input_file_reference.md
- doc/input_schema.yaml
- src/fem/kernel/doc/nonlinear_controller_theory.md
- todo/debt_register.md
- devlog/README.md
