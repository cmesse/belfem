# /circuit Restart Format v2: Per-Component State (DR-124 residue)

**Date:** 2026-08-29
**Purpose:** Remove the accepted v1 limitation of the circuit memdump — L/C histories, component currents, terminal-pair state and the switch latch now dump and restore; continuation is unit-level exact
**Module:** circuit, tests/circuit

## Round

Own plan+audit → code+audit thread: `tmp/ai_exchange/dr124_circuit_restart_v2.md`. Plan audits: both vendors request-changes, convergent; all amendments adopted (see the reconciliation section there). Everything below is reviewed, not verified — the circuit suite rerun and DR-124's coupled R6 remain the gates. Scope held hard: `src/circuit` + `tests/circuit` + the module guide + this bookkeeping; no controller, container, io, or executable edits (parallel sessions hold parts of the tree; peers were coordinated with before editing).

## Landed

- `Component` gains virtual `save_state`/`load_state( hid_t, prefix )` with empty defaults; overridden by `Inductor` (`h`/`i` registers + current), `Capacitor` (`h`/`v` + current), `FEMTwoTerminals` (three registers + `mIn`/`mVn`, and `set_current(mIn)` on load — a terminal pair is not an unknown-current dof, so the coupled Jacobian would otherwise inject 0 instead of I_n on the first post-restart step; Grok's blocker), `Switch` (`is_closed`/`is_switched`; rollback snapshots mirrored, not dumped), `VoltageSource` (`value_time`) and `CurrentSource` (current) so a stamp-without-shift is coherent (Codex's amendment).
- New header-only `fn_circuit_state_io.hpp`: `save_shift_register`/`load_shift_register` over the container's public API (`cl_ShiftRegister.hpp` untouched). Dumps carry each register's capacity; load hard-errors on mismatch in EITHER direction (order-2 dump into order-3 would silently become a partial BDF3 window — Grok's blocker 2). Rebuild is `clear()` + oldest-first `push()`; the full-dump revert flag differs (`CanRevert` vs live `CanRevertFull`) until the next push, so `load_state(); shift_back();` is forbidden and documented — every production reject path shifts first.
- `ElectricalCircuit::save_state` walks components under `c%03u_` creation-index prefixes with a `n_components` count and per-component `type` dataset; `load_state` validates count and types before dispatching (a reorder across types refuses instead of swapping histories — Codex's identity amendment; a same-type reorder of an edited deck is NOT detectable and is out of contract, a restart continues the same deck). Pre-v2 dumps (no `n_components`) hard-error with the delete-the-memdump message.
- L/C loaders rebuild companions via `update_companions()` when the restored registers are non-empty — the round-3 DR-138 helper is what makes the loaded objects coherent for both legal continuations. NO `compute_current()` walk after load (Grok's veto: it would clobber restored inductor current).
- Tests: `WarmRestartContinuationEquivalence` (500-step dump, 500-step continuation vs uninterrupted twin — inductor current, capacitor current, node voltage at 1e-12 per step), `WarmRestartPartialHistory` (dump after ONE order-2 step: partial-fill rebuild, the state the big test never hits), `WarmRestartSwitchLatch` (fired latch restored WITHOUT a shift; no double toggle), `WarmRestartRefusesOldFormat` (hand-written v1 dump → throw).
- Docs: `circuit_usage_guide.md` §9 publishes the v2 layout and the loader's contract; the v1 standalone-L/C limitation is marked historical. `todo/restart_circuit_verification.md`: format-publication bullet ticked; R6's re-fire sub-bullet rewritten for latch persistence, **pending Christian's ack** of that design change (the v1 self-correcting re-fire is superseded).

## Code audit

Grok: APPROVE — all nine reconciliation amendments traced present; it specifically hunted the reverse-loop off-by-one, the `hid_t&` binding, DR-138 assert trips, fixture index errors, and a v1 dump slipping past `n_components`, and cleared each. Codex: request-changes, four items, all applied: (1) the blocker — an empty-history dump violated the stamp-without-shift contract (fresh `mRL`/`mRC` indeterminate); fixed by `load_state` seeding component timesteps/BDF1 companions from the restored `delta_time` before the component walk, with the `WarmRestartBeforeFirstShift` regression; (2) validation-first guard order (format/count/type checks before any state mutation); (3) reorder wording made honest (cross-type refuses; same-type on an edited deck is undetectable — the same-deck rule is a contract); (4) the roundtrip comparisons upgraded to bitwise `EXPECT_EQ`. Grok's closure items applied: `VoltageSource::mValueTime = 0.0` in-class init, paired-register lockstep `BELFEM_ERROR` in the L/C loaders (corrupt-dump OOB guard), and the stale `restart_circuit_verification.md` preamble rewritten. Honest coverage note kept from Grok: TERMINALPAIR state restoration has no unit test — that path belongs to the coupled R6.

## Gate result (2026-08-29, Christian's run)

Full suite 15/15 green. The circuit binary's 85 tests all passed on first run: the 1e-12 continuation equivalence (full and partial registers), the pre-first-shift dump, the switch latch, the old-format refusal — and the DR-138 RLC twin, red this morning, is green under its 1e-3 bound. Predictions P1/P2 of this round and the DR-138 round's pre-registration all confirmed. The v2 restore and the reject-retry semantics are VERIFIED at unit level; the coupled R6 A/B is now the only run gate on DR-124.

## Notes

- P3 correction absorbed from the plan audit: R6 does not become a formality — the controller caps the first post-restart Δt at `mDeltaTimeInitial` regardless of circuit state. This round removes the last KNOWN circuit-state gap, nothing more.
- Syntax gates green on all eight touched circuit TUs and the test file (module flag sets).
- Format change is pre-release-legal per DR-112's ruling; a restart is a continuation of the same deck, never a migration path for an edited one — every mismatch is a hard error.
