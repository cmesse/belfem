# Circuit Cluster: DR-118 + DR-121 + DR-123 + DR-124 Fixed

**Date:** 2026-08-29
**Purpose:** One plan+audit → code+audit round over the four open circuit-module debt rows; all four fixes landed, two by-catch rows filed
**Module:** circuit, executables, tests/circuit

## Round

Pre-registered plan, both-vendor plan audit (Codex + Grok, convergent approve-with-amendments), implementation, then the both-vendor code audit: Codex returned request-changes with five findings (demo banner one dt behind the solved state, dt-floor message vs predicate mismatch, HDF5 test cleanup missing on failure paths, premature audit wording in this devlog, bare `k` test loop indices), of which the first four were applied the same session and the fifth answered with the `tests/io/test_HDF5.cpp` bare-`k` precedent; Grok returned APPROVE with non-blocking follow-ups, all applied: both `fixed =` register cells rewritten to the gates actually landed, `WarmRestartKeepsState` now takes two steps so the dumped `prev_x` is nonzero (its comparison was zeros-vs-zeros), explicit `<algorithm>` include, and the DR-138 mechanism clause corrected to the post-cluster call order. Grok also discharged the plan round's HDF5-wiring risk: `Add_Test.cmake` links the whole `libbelfem` and adds `src/io` includes regardless of `LIBLIST`. Exchange: `tmp/ai_exchange/circuit_cluster_dr118_121_123_124.md`. Everything below is reviewed, not verified — the run gates are owed (`make check-fast`, one demo run, and DR-124's coupled R6).

## Landed

- **DR-121** — `compute_MNA_matrix` gets a `SWITCH` case advancing the shared unknown-current counter (`cl_ElectricalCircuit.cpp`); no MNA stamp (switch stamps stay in the Jacobian walk; both matrices share the `mVertices` pattern incl. the switch diagonal, so the raw-data merge stays aligned). The ngspice factory's v1 switch-before-V refusal is removed; its regression test flipped to acceptance.
- **DR-123** — `Switch::shift_back()` override restores `mIsClosed`/`mIsSwitched` from snapshots taken at the top of every `Switch::shift()` (ctor-initialized, one-deep like the L/C `ShiftRegister`s); `ElectricalCircuit::shift_back()` walk gets the `SWITCH` case (no `compute_current` — the switch current is an unknown-current dof restored by `update_components()`). Live consumer: the coupled `reset_timestep()` path.
- **DR-124** — the one-shot block in `compute_MNA_matrix` zero-initializes `mX`/`mPrevX` only when either length differs from the dof count, so `load_state()`'s restored vectors survive the first stamp. Row stays open on the coupled R6 gate; `load_state` still restores only `{time, delta_time, x, prev_x}` (no L/C histories, no switch latch).
- **DR-118** — `electricalCircuit.cpp` time loop restructured to controller parity: `set_timestep → shift → compute_MNA_matrix` at the top of every attempt, reject = `shift_back` + halve + `continue` (same interval re-solved, first-step rejection legal), grow path no longer restamps, dt-floor `BELFEM_ERROR` ends the previously-possible infinite halving. Output semantics deliberately corrected: t=0 row before the loop, accepted-state rows with aligned time tags, no rows for rejected attempts, step numbering off-by-one gone.

## By-catch

- Fixed in the DR-118 file: `delete tFunction` after `delete tCircuit` was a double free (`VoltageSource` owns and deletes its `SourceFunction`). Structurally gone: the old first-stamp-before-first-shift path that solved against uninitialized `miLh`/`miCh`.
- Filed **DR-138**: suspected inductor current loss through `ElectricalCircuit::shift_back()` — revert then `compute_current()` against un-reverted `mRL`/`miLh` companions; the retry pushes the wrong current into history. Inductor-specific, live on demo and coupled reject paths. Static trace (Grok), endorsed (Codex), no gate yet. Pre-registered reading of a red `RLCRingAcrossRejectedStep`: this defect, not a DR-118 regression.
- Filed **DR-139**: `~ElectricalCircuit` deletes `mJ` but never `mMNA` (leak; deliberately not fixed in this diff on both auditors' scope instruction).

## Tests

New `tests/circuit/test_ElectricalCircuit.cpp` (wired into `SOURCES`): `SwitchBeforeVoltageSourceAgrees` (mixed creation order vs V-first twin, topology chosen so a misplaced stamp falls outside the switch column's sparsity), `SwitchLatchRevertsOnRejectedStep` (fire → revert → re-fire through the circuit walk), `WarmRestartKeepsState` + `WarmRestartControllerParity` (HDF5 roundtrip discriminators for the DR-124 wipe, `#ifdef BELFEM_HDF5`; explicitly NOT the R6 gate), `RLCRingAcrossRejectedStep` (order-2 series RLC at resonance, forced reject-and-retry vs never-rejected twin on the realigned grid, analytic V/R amplitude as loose sanity). `SwitchBeforeVoltageSourceRefused` became `SwitchBeforeVoltageSourceAccepted`.

The plan's P1 ("T1 throws today") was retracted in reconciliation: the throw is topology-dependent; voltage agreement is the discriminator. T4 gates library semantics only — the demo executable is review-only until someone runs it.

## Gate results (2026-08-29, Christian's run of the circuit suite, 79 tests)

- `SwitchBeforeVoltageSourceAgrees`, `SwitchLatchRevertsOnRejectedStep`, `WarmRestartKeepsState`, `WarmRestartControllerParity`: **green** — DR-121, DR-123 and DR-124's wipe fix are verified at unit level (rows updated). The flipped factory test `SwitchBeforeVoltageSourceAccepted` is green too.
- `RLCRingAcrossRejectedStep`: **red**, max post-retry deviation 7.4e-3 A on the 1 A carrier vs the 1e-3 bound — the pre-registered DR-138 outcome, and the magnitude matches the one-step-ahead-companion mechanism (~amplitude·2π·dt/T ≈ 6e-3 with dt = T/1000), not truncation (orders smaller) and not a counter bug (would be O(1)). DR-138 is now runtime-confirmed; its row carries the number. Per the pre-registration this is NOT a license to widen the DR-118 fix; DR-138 takes its own round.
- Operational consequence: the circuit suite stays red until DR-138's fix lands (or the test is temporarily gated with a DR-138 reference) — routed to Christian.

## DR-138 fix (round 3, same thread)

Christian: fix now. Plan round split the vendors, and the tree decided it: Codex approved a Switch-style snapshot of `mRL`/`miLh` at the top of `shift()`; Grok refuted it — every live path (`take_step`, the demo, both controller inits) calls `set_timestep()` before `shift()`, and `Inductor::set_timestep` clobbers `mRL` with the BDF1 seed `L/dt`, so the snapshot would capture the seed, not the accepted BDF2 `1.5·L/dt`, leaving a residual `V_n·dt/(3L) ≈ 2e-3 A` at T4's reject point (which sits, usefully, at the worst case: current zero, inductor voltage peak). Grok's trace checked out; Codex's magnitude derivation (`I_wrong − I_n ≈ h·I′ ≈ 6.3e-3 A`, matching the observed 7.4e-3, three orders above truncation) stands.

Landed (Grok's required shape): `Inductor::shift_back()` reverts the registers, then — when non-empty — recomputes the companions from the reverted state through a new private `update_companions()` shared with `shift()`; the reverted registers are exactly the post-accepted-shift state, so the walk's `compute_current()` reconstructs the accepted current bit-for-bit before the retry pushes it back into history. Empty-register guard keeps the `set_timestep` seed and zero source on a first-step rejection (exact at the zero cold state). `Capacitor` gets the identical recompute copy for rollback coherence (its history is voltages, so its half is latent). In-class `miLh = 0.0`/`miCh = 0.0` provide explicit no-history initialization (defensive; on a real first rejection the rejected shift had already written the companions, so no live path read them uninitialized). No new members. Syntax gates green.

Pre-registered rerun outcome: `RLCRingAcrossRejectedStep` green under 1e-3; a residual red would be a finding about the bound (variable-step BDF residual after the two half-steps), never a license to loosen it. **Rerun RAN the same day (Christian, full suite 15/15): GREEN — the DR-138 fix is verified at unit level.**

Code audit, both vendors: Grok APPROVE (traced T4's rejection as the `CanRevertFull` register mode with the recompute load-bearing for `miLh`; its optional precondition assert inside `update_companions()` applied in both classes). Codex request-changes on two wording overclaims only — capacitor first-step comment ("zero source" only on a zero-voltage cold start) and the "removes indeterminate reads" claim (the inits are defensive; the rejected shift had already written the companions) — both applied; no code defect from either vendor. Grok's "fabricated Codex clause" flag was a snapshot race (it read the exchange as dispatched, before Codex's section landed) and is resolved in the exchange, not by editing the accurate clause.

Shared-checkout coordination: peers pinged before editing; src/circuit confirmed unclaimed by all four active sessions; register and devlog index taken append-only (DR-127 was struck by the parallel session in the meantime — its 24-`[P]` recount stands, unchanged by this round).

## Notes

- DR-118 fixed in place; ngspice parser plan Phase 4 (`CircuitSolver` extraction) stays frozen and inherits the corrected loop.
- DR-124 process: the row demanded a discriminator before the fix; fix and discriminator landed together, and the red/green demonstration (guard reverted vs applied) is available at gate time — Christian's call to exercise or waive.
