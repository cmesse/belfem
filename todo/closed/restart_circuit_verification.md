# Coupled-Circuit Restart: Verification & Documentation Tail

> **CLOSED 2026-09-03** (todo/ currentness sweep, round 3): the unit restart suite exists (`tests/circuit/test_ElectricalCircuit.cpp`); R6 coupled A/B is a run gate. Status lines and checkboxes below are as they stood at closure and are not maintained.

**Date:** 2026-07-03
**Purpose:** Close out the three verification/documentation items left after the coupled-circuit
dynamic-state restart was **implemented** (`{time, delta_time, x, prev_x}` in the `/circuit` group
of `memdump.hdf5`, scatter on load, R4 length guard, Controller wiring — all landed and tri-AI
audited). The implementation plan is closed at `todo/closed/restart_circuit_state.md`; this file
tracks only what its Definition-of-Done still lists as open.
**Module:** `src/circuit`, `src/fem/kernel` (Controller), `src/fem/iwg` (R5)
**Status:** OPEN — R5 and the coupled R6 remain; the format-publication item closed 2026-08-29.
**Superseding note 2026-08-29 (format-v2 round):** the 2026-08-09 sweep's findings are
historical — `tests/circuit/test_ElectricalCircuit.cpp` now carries the unit-level restart
suite (roundtrip, continuation equivalence full/partial, switch latch, old-format refusal,
pre-first-shift dump), and `circuit_usage_guide.md` §9 documents the `/circuit` v2 group;
the standalone-L/C v1 limitation is superseded by per-component state. R6 (coupled A/B)
is still the item with teeth. Tracked as `debt_register.md` DR-124 (DR-38's successor row).
(Correction, caught by the Codex triage leg: an earlier version of this note said BDF5 had
become the default scheme on 2026-08-07. It did, for one day, and was rolled back to BDF1
the same day with the Anderson residual fix — `cl_FEM_Controller.hpp:189`. R5's multistep
history question stands on its own merits; it is not newly urgent.)

> **DR-124 (filed 2026-08-27, from the DR-115 fix audit):** static trace suggests R6 would
> fail today — the first post-restart `compute_MNA_matrix()` re-runs its one-shot sizing
> block (`mMNA == nullptr`) and re-constructs `mX`/`mPrevX` as zeros, wiping what
> `load_state()` just restored. If the R6 test is ever written, expect it to be red before
> the DR-124 fix; details in the register row and `tmp/ai_exchange/dr115_fix_round.md`.

## Open Items

- [ ] **R5 — Confirm `IWG_Timestep` cold-starts cleanly on a load.** Verify the multistep history
  stages seed from the loaded field / restart the order ramp rather than reading stale stage
  buffers. FEM side, independent of the circuit state. Trace `mFieldData` handling
  (`cl_IWG_Timestep.cpp:203,214,237,267-268,333-334`) across a `load_memdump`.

- [ ] **R6 — End-to-end coupled restart test.** No such test currently exists. CORC with circuit:
  run to T, dump, restart, compare against the uninterrupted run within tolerance (R12
  methodology from the closed meshfile plan). Must include:
  - a time-dependent source crossing the restart point (exercises the restored clock),
  - a switch that fired *before* the dump (v2 restores the latch from the dump — the check is
    that it comes back fired and does NOT re-fire; the v1 "self-correcting re-fire" is gone,
    superseded 2026-08-29 by latch persistence, **pending Christian's ack** of that design change),
  - a first-step rejection after restart (exercises `mPrevX` / `shift_back()`).
  Serial + one MPI width (circuit is rank-0 only; the broadcast pattern covers the rest).
  Note: the format-v2 round (2026-08-29) landed unit-level continuation equivalence
  (`WarmRestartContinuationEquivalence`, full and partial histories, 1e-12 twin traces) —
  R6's remaining value is the coupled path itself, incl. the controller's restart-entry Δt cap.

- [x] **Publish the `/circuit` restart format in the module doc.** Done 2026-08-29:
  `circuit_usage_guide.md` §9 documents the v2 layout (vectors + `n_components` +
  `c<NNN>_`-prefixed per-component state) and marks the v1 standalone-L/C limitation
  superseded — v2 dumps and restores the histories, so the limitation paragraph this bullet
  originally asked for is recorded as historical instead.

## Related

- Closed implementation plan: `todo/closed/restart_circuit_state.md` (full gap table, audit trail,
  §6 format sketch, v1-limitation rationale).
- Parent (closed): `todo/closed/meshfile_refactor_plan.md` (R6b).
