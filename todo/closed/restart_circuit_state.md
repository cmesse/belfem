# Coupled-Circuit Dynamic-State Restart Plan

**Date:** 2026-07-01
**Purpose:** Persist the coupled electrical circuit's dynamic state in the `memdump.hdf5` restart
file so a CORC restart resumes with a **warm circuit** instead of re-creating it cold from
`input.conf`. This is the one restart item left open after the `.bfm` mesh save/load refactor
(`todo/closed/meshfile_refactor_plan.md`, R6b). It extends the **restart** state
(`Controller::save_memdump`/`load_memdump`), not the `.bfm` enrichment cache.
**Module:** `src/circuit` (state + serialization), `src/fem/kernel` (Controller wiring)
**AIs involved:** Claude (exploration + plan), Codex (prose pass + technical re-check), Grok
third-voice audit
**Status:** CLOSED 2026-07-03 — implementation landed and tri-AI audited (R1–R4 done). The three
remaining verification/doc items (R5 IWG cold-start check, R6 end-to-end restart test, `/circuit`
module-doc publication) were spun out to the active `todo/restart_circuit_verification.md`. Original
in-progress status preserved below for the historical record. **R2/R3 implemented
(Christian) all on 2026-07-01**: `save_state`/`load_state` of `{time, delta_time, x, prev_x}`,
the `update_components()` scatter, Controller wiring, hphirun reorder (D1 fixed), and the R4
length guard (Claude, user-approved). **Tri-AI audit passed 2026-07-01** (Codex + Grok, §8):
one fix applied (`cl_Circuit.hpp` now includes `hdf5_types.hpp` for `hid_t` self-containment),
one residual deferred (release-mode `tStatus` checks, D6 precedent). R1 verified closed.
**Codex re-check 2026-07-01:** the implementation is present in the current tree; this file is
partly historical. **Still open: R5 (IWG cold-start check), R6 (end-to-end restart test —
Christian), and module-doc publication of the `/circuit` restart format/v1 limitation.**

> **Scope guards (from the task brief + the 2026-06-30/07-01 restart-philosophy decisions):**
> - **Structure is not persisted.** The circuit is always rebuilt from `input.conf`
>   (`ElectricalCircuitFactory`); only the *dynamic* state is stored in the restart file and
>   re-attached to the rebuilt circuit. Input stays the source of truth (same D21 philosophy as
>   the `.bfm`).
> - **Minimal state (decided 2026-07-01, Christian): persist only the solution vectors `mX` and
>   `mPrevX` (+ the clock, see O2) — no ShiftRegister serialization.** Unlike the mesh, the
>   circuit topology never changes, so everything else is reconstructed at runtime: the restored
>   solution is scattered back into nodes/components, and the **first regular `shift()` re-seeds
>   the BDF registers by itself** (it pushes the components' *current* values —
>   `cl_Inductor.cpp:49`, `cl_FEMTwoTerminals.cpp:36-39`), with the order ramping from 1. This
>   matches the FEM side exactly: state restart, cold-started integrator.
> - **Deliberately NOT persisted:** the FEM BDF multistep history (`IWG_Timestep::mFieldData`/
>   `mH`), the controller's adapted Δt (`mDeltaTime`), and — per the decision above — the circuit's
>   L/C/terminal-pair ShiftRegister histories. Δt is re-read from the (possibly edited) input file.
> - Circuit solve is **rank 0 only** (`cl_FEM_Controller.cpp:147`); save/load of circuit state is
>   rank-0 only, matching `save_memdump`'s existing guard (`:1783`).
>
> **Known v1 limitation (accepted):** a standalone **inductor's current is a true state** that the
> Norton companion model hides — `Inductor::compute_current()` derives it from the
> history-dependent `miLh` (`cl_Inductor.cpp:35-38`), so it is *not* recoverable from `mX` (node
> voltages) alone. For CORC-class circuits (sources + resistors + terminal pairs) `mX` carries
> everything that matters; circuits with standalone L (or C) restart with a perturbed L/C state.
> If that ever matters, the fallback is per-component register persistence (the earlier draft's
> `Component::save_state` design) — deferred, not planned.

---

## 1. Historical Cold-Restart Failure Mode

This section describes the **pre-fix behaviour** that motivated the implementation. The core
circuit-state persistence described here is now implemented; the current open work is tracked in
§4/§7.

The restart file already works for the FEM side: `Controller::save_memdump`
(`cl_FEM_Controller.cpp:1781-1822`, rank 0) writes `meta` (timestep/timestamp/running_timestep/
checksum), `fields`, `globals` (+ `fields2`/`globals2` for the thermal mesh);
`Controller::load_memdump` (`:1825-1900`) restores them, broadcasts time to all ranks, and the
run resumes warm. Current `hphirun` ordering attaches the circuit first, then loads the memdump
(`set_circuit` before `load_memdump`), and saves at each saved step.

Historically, the circuit saved nothing: no serialization code existed in `src/circuit/`. On
restart `ElectricalCircuitFactory` rebuilt the circuit from `input.conf` with freshly-constructed
components, and `Controller::set_circuit` attached it after the memdump load. That cold circuit
failed against warm FEM fields in these ways:

| Failure | Mechanism | Evidence |
|---|---|---|
| **Wrong drive for the whole rest of the run — circuit clock replays from zero** | `ElectricalCircuit::mTime = 0.0` is private with **no setter** (`cl_ElectricalCircuit.hpp:92`); only `shift()` advances it (`cl_ElectricalCircuit.cpp:180`). Time-dependent sources are evaluated from *this* clock: `VoltageSource::shift` → `mValueTime = mFunction->compute(aTime)` (`cl_VoltageSource.cpp:51-54`, same for `CurrentSource`). A restart at T therefore replays the source waveform from t≈0 — a sine drive, ramp, or pulse is time-shifted forever, not just transiently. | verified 2026-07-01 |
| **Switch replays** | `Switch::shift` fires when `aTime >= mSwitchTime && !mIsSwitched` (`cl_Switch.cpp:43-50`). Cold restart resets `mIsSwitched=false` **and** the clock, so a switch that fired before the dump re-fires mid-run at the wrong wall-clock time. *(Self-corrects once the clock is restored — see §3 row 6.)* | `cl_Switch.hpp:27-33` |
| **Bad Newton start + broken first-step reject** | `mX` (node voltages + unknown branch currents) and `mPrevX` restart at zero (`cl_ElectricalCircuit.hpp:71-74`). The first MNA Newton solve starts from zero instead of the converged state, and a timestep rejection on the first restarted step (`Controller::reset_timestep` → `mCircuit->shift_back()`, `cl_FEM_Controller.cpp:1074`) "restores" zeros via `mX = mPrevX` (`cl_ElectricalCircuit.cpp:209`). | verified 2026-07-01 |
| **Cold BDF history in L/C/terminal pairs** | `Inductor`/`Capacitor`/`FEMTwoTerminals` registers (`mI`/`mV`/`mH`) start empty; the discretized history source (`miLh`/`miCh`) is wrong in the first steps. | `cl_Inductor.cpp:44-67` |

Purely algebraic components (Resistor, Diode, Superconductor) and the sources' `mValueTime` are
**not** failures because they are recomputed from the restored voltages/time on the next
`shift()`/`compute_current()` pass (see §3 rows (a)). Under the minimal-state decision, the L/C/TP
history row is likewise **accepted** rather than fixed: the registers re-seed from the scattered
solution on the first `shift()`, order ramping from 1 (see scope guards + v1 limitation).

**Bottom line:** the FEM side resumes at time T while the circuit resumes at t=0 with zero state
and a replaying drive; the coupled run silently diverges from an uninterrupted one. The fix is a
minimal state blob — one clock and two vectors — plus a scatter on load; everything else
reconstructs at runtime because the circuit topology is input-defined and never enriched.

---

## 2. Architecture: `{time, x, prev_x}` in a `/circuit` Group, Everything Else Reconstructed

Reuse the existing restart-file pattern: `save_meta`/`save_fields`/`save_globals` each own an HDF5
group, with rank-0 write and rank-0 read + broadcast. Add:

- `ElectricalCircuit::save_state( … )` / `load_state( … )` — writes/reads `mTime`, `mX`, `mPrevX`.
  As class members they can restore the private `mTime` without adding a public setter. On load,
  **scatter `mX`** back into node voltages, unknown branch currents, and algebraic
  `compute_current()` — the exact loop `shift_back()` already runs
  (`cl_ElectricalCircuit.cpp:214-255`, minus the register `revert()`s); factor it into a shared
  helper rather than duplicating it.
- Controller wiring inside the existing `save_memdump`/`load_memdump` (a `circuit` group written
  when `mCircuit != nullptr`, rank 0 only).
- **No per-component serialization.** The first `initialize_timestep` after the load runs the
  normal sequence — `set_timestep` → `shift()` → MNA solve (`cl_FEM_Controller.cpp:149-152`) —
  and `shift()` itself re-seeds the registers from the scattered state: `Inductor::shift` pushes
  `get_current()`, `Capacitor::shift` its voltage, `FEMTwoTerminals::shift` `mIn`/`mVn`
  (`cl_Inductor.cpp:44-67`, `cl_FEMTwoTerminals.cpp:34-40`). Sources recompute from the restored
  clock; a fired switch re-fires immediately (correctly) on that same first `shift()`.

**Re-attachment keying is a non-issue in this design:** no per-component records exist. The only
consistency requirement is that the rebuilt circuit's solution-vector length
(`mNumberOfNodes-1 + mNumberOfUnknownCurrents`) matches the stored `x` — a cheap guard against an
edited `input.conf` (R4, policy O4). *(The earlier draft's creation-index + fingerprint scheme is
superseded by the 2026-07-01 minimal-state decision; see O1.)*

---

## 3. Gap Table

Class: **(a)** rebuilt deterministically on load (from what) / **(b)** open question / **(c)** must
be saved.

| # | State | Needed for | Class | Citation / rationale |
|---|---|---|---|---|
| 1 | Circuit clock `mTime` | source waveforms, switch timing | **(c) save** (or seed from mesh time — O2) | `cl_ElectricalCircuit.hpp:92`; advanced only in `shift()` (`.cpp:180`); no setter exists |
| 2 | Solution vector `mX` (node voltages + unknown branch currents) | Newton start, BC fix values, register re-seeding | **(c) save**, scatter on load | `cl_ElectricalCircuit.hpp:71`; scatter pattern in `shift_back()` (`.cpp:215-224`) |
| 3 | `mPrevX` | `shift_back()` on a first-restarted-step reject | **(c) save** (cheap; else first reject is corrupt) | `cl_ElectricalCircuit.hpp:74`; reject path `cl_FEM_Controller.cpp:1074` |
| 4 | L/C register histories (`mI`/`mV`/`mH`) | BDF companion model (`miLh`/`miCh`) | **(a — accepted approximation)** re-seeded by the first `shift()` from the scattered state, order ramps from 1 (decided 2026-07-01, Christian); exact only where the component current is in `mX` — see the v1 limitation | `cl_Inductor.cpp:44-67` |
| 5 | `FEMTwoTerminals` `mIn`/`mVn` + registers | first coupled iteration | **(a — accepted)** refilled by the first FEM solve's `compute_circuit_current`; one-iteration transient accepted | `cl_FEMTwoTerminals.cpp:34-40,55-60` |
| 6 | Switch `mIsClosed`, `mIsSwitched` | correct topology branch | **(a)** given a restored clock: rebuilt `mIsClosed` is the input-file initial state, and the first `shift()` at T ≥ `mSwitchTime` toggles it immediately — the post-fire state is reproduced. (Only exact for the single-fire semantics the `mIsSwitched` guard implies.) | `cl_Switch.cpp:43-50` |
| 7 | Source `mValueTime` | drive value | **(a)** recomputed by the next `shift()` from the restored clock | `cl_VoltageSource.cpp:51-54` |
| 8 | R/D/SC currents | output/consistency | **(a)** `compute_current()` from restored node voltages — same recompute `shift_back()` does | `cl_ElectricalCircuit.cpp:227-245` |
| 9 | MNA matrix `mMNA`, Jacobian `mJ`, `mRHS` | solve | **(a)** rebuilt every step (`compute_MNA_matrix` in `initialize_timestep`) | `cl_FEM_Controller.cpp:151` |
| 10 | Component identity / keying | re-attachment | **moot** — no per-component records in the minimal design; only the `x`-length guard remains (R4) | superseded, see O1 |
| 11 | Relaxation `mOmega`, `mDeltaTime` | solver behaviour | **(a)** config / set fresh each step (`set_timestep`, `cl_FEM_Controller.cpp:149`) — Δt deliberately not persisted (scope guard) | `cl_ElectricalCircuit.hpp:89,95` |

### 3.1 Cross-cutting findings

- **The clock (row 1) is the highest-value single item.** `x`/`prev_x` alone do not fix the
  source-replay failure — the drive is evaluated from the circuit's own clock. Conversely,
  restoring the clock alone already fixes sources *and* (row 6) the switch. It is one scalar;
  store it (O2).
- **Load-order hazard (closed by O3/D1):** `hphirun.cpp` used to call `load_memdump` before
  `set_circuit`, so the Controller had no circuit to restore into at load time. Current code calls
  `set_circuit` first.
- **The scatter is the load-side workhorse.** Restoring `mX` without scattering leaves nodes and
  component currents at zero, and then the first `shift()` would seed the registers from zeros —
  silently reproducing the cold-history failure. Scatter must run inside `load_state`, and the
  register re-seeding (row 4/5) depends on it.

---

## 4. Ordered Steps

### 4.0 Implementation Progress (updated 2026-07-01)

**Landed (Christian, 2026-07-01):** `Circuit::save_state`/`load_state` pure virtuals
(`cl_Circuit.hpp`), `ElectricalCircuit::save_state`/`load_state` writing/reading
`{time, delta_time, x, prev_x}` (`cl_ElectricalCircuit.cpp:965-987`, `#ifdef BELFEM_HDF5`), and the
Controller wiring: `save_memdump` writes the `circuit` group inside the rank-0 block, `load_memdump`
reads it behind `mCircuit != nullptr && group_exists` (`cl_FEM_Controller.cpp:1818-1823/1869-1874`).

- [x] **D1 — HIGH (load is a silent no-op) — hphirun ordering.** ~~`hphirun.cpp` calls
  `load_memdump` (`:78`) **before** `set_circuit` (`:81`), so at load time `mCircuit == nullptr`
  and the new `circuit`-group guard skips the restore without a message.~~ **Fixed 2026-07-01
  (Christian)** — `set_circuit` now precedes `load_memdump` (O3 option (i); verified in the
  diff). The optional belt-and-suspenders `BELFEM_ERROR` for "`circuit` group present but no
  circuit attached" was not added — low priority, protects only future drivers.
- [x] **Load-side scatter.** ~~`load_state` restores `mX` but does not push it into node
  voltages / component currents.~~ **Fixed 2026-07-01 (Christian)** — new private
  `ElectricalCircuit::update_components()` (node voltages + unknown branch currents from `mX`),
  factored out of the two inline copies in `shift_back()` and `compute_correction()` (pure
  refactor there) and called at the end of `load_state`. Verified adequate (Claude): the
  algebraic R/D/SC `compute_current()` calls are deliberately **not** made at load — nothing
  consumes those currents before the first Newton loop recomputes them, and calling L/C/TP
  `compute_current()` at load would read uninitialized companion values. `Capacitor::shift`'s
  node-voltage read (`cl_Capacitor.cpp:51`) is now correctly seeded.
- [x] **Length guard (R4).** ~~`update_components()` indexes `mX` up to
  `mNumberOfNodes+mNumberOfUnknownCurrents-2` — a stored `x` from an edited (smaller) circuit is
  OOB in release.~~ **Fixed 2026-07-01 (Claude, user-approved)** — `BELFEM_ERROR` in `load_state`
  checks `mX` and `mPrevX` lengths against `mNumberOfNodes - 1 + mNumberOfUnknownCurrents` before
  `update_components()`; message follows the O4 lean ("delete the memdump to restart with a cold
  circuit").
- **`delta_time` added to the format (Christian).** Stored as self-description; on load it is
  overwritten by `set_timestep( controller Δt )` before first use (`cl_FEM_Controller.cpp:149`) —
  deliberate: input stays authoritative for Δt (scope guard). See the O2 resolution.

- [x] **R1 — Verify the runtime-reconstruction chain end to end (read-only).** **Done 2026-07-01
  (tri-AI: Claude + Codex + Grok, all citations independently verified).**
  `BDF::compute_coefficients` sizes from `mH.size() + 1` (`cl_BDF.cpp:33`) — a one-entry
  register yields a well-posed 2×2 Vandermonde solve = exact backward-Euler; downstream loops
  key off `coefficients().length()`, so only `mI(0)`/`mV(0)` participate. `shift()` pushes
  before computing coefficients (`cl_Inductor.cpp:47-50`), so the register is never empty at
  solve time. Sources/switch pick up the restored clock on the first `shift()`; `revert()`
  unavailability after restore is safe (`shift()` always precedes any `shift_back()` on the
  restarted step, `cl_FEM_Controller.cpp:149-150` vs `:1074`).
- [x] **R2 — `ElectricalCircuit::save_state`/`load_state`.** **Done (Christian, 2026-07-01):**
  save/load of `{time, delta_time, x, prev_x}` + the load-side scatter via the factored-out
  `update_components()` helper (see §4.0). The length guard is tracked separately (R4).
- [x] **R3 — Controller + hphirun wiring.** `save_memdump`: add the `circuit` group inside the
  rank-0 block (`cl_FEM_Controller.cpp:1801` vicinity), guarded on `mCircuit != nullptr`. Load
  side per O3's resolution (recommended: reorder `hphirun.cpp` so `set_circuit` precedes
  `load_memdump` — `set_circuit` only stores the pointer and partitions BCs, no dependency on
  loaded fields; re-verify). Loading an old memdump **without** a `circuit` group must stay legal
  (`group_exists` guard → cold circuit; warning optional), so pre-feature dumps keep working.
  **Done (Christian, 2026-07-01):** Controller side + the hphirun reorder (D1). Residue: the
  optional no-circuit hard error (§4.0 D1 note) and the cold-circuit *warning* when an old dump
  has no `circuit` group (currently silent skip — acceptable). (after: R2)
- [x] **R4 — Consistency guard.** On load, `BELFEM_ERROR`-check the stored `x` length against the
  rebuilt circuit's `mNumberOfNodes-1 + mNumberOfUnknownCurrents`. Mismatch policy per O4.
  **Done 2026-07-01** (both `x` and `prev_x` checked in `load_state`). (after: R2)
- [ ] **R5 — Confirm `IWG_Timestep` cold-starts cleanly on a load** (seed history stages from the
  loaded field / restart the order ramp) rather than reading stale stage buffers. *(Carried from
  the parent plan; FEM side, independent of R1-R4.)*
- [ ] **R6 — End-to-end coupled restart test.** CORC with circuit: run to T, dump, restart,
  compare against the uninterrupted run within tolerance (the R12 methodology). Include: a
  time-dependent source crossing the restart point, a switch that fired *before* the dump, and a
  first-step rejection after restart (exercises `mPrevX`). Serial + one MPI width (circuit is
  rank-0 only; the existing broadcast pattern covers the rest). (after: R3, R4)

---

## 5. Open Design Questions (not silently decided)

- **O1 — Re-attachment keying. RESOLVED 2026-07-01 (moot) → minimal-state design (Christian).**
  No per-component state records exist, so no keying is needed; the earlier creation-index +
  fingerprint scheme is superseded. Only the `x`-length guard remains (R4). Revisit only if the
  deferred per-component register persistence (v1-limitation fallback) is ever built.
- **O2 — Persist the circuit clock, or seed it from the mesh time? RESOLVED 2026-07-01 → store it
  (Christian's implementation).** `save_state` writes `time` *and* `delta_time`
  (`cl_ElectricalCircuit.cpp:965-987`). On `delta_time`: the (x, prev_x) pair is only
  interpretable as a rate together with the Δt between them, so storing it keeps the file
  self-describing — but note neither `prev_x` nor `delta_time` is currently *consumed* after a
  restart: the first `initialize_timestep` runs `set_timestep(controller Δt)` (`:149`) and then
  `shift()`, whose first statement `mPrevX = mX` (`cl_ElectricalCircuit.cpp:177`) overwrites both
  before any reader (rejection can only occur after `shift()`, and its revert uses the *new* Δt
  symmetrically). They ride as self-description/future-proofing; the live state is `time` + `x`.
  **Caution kept deliberate:** the loaded `delta_time` must never win over the input-file Δt —
  today `set_timestep` overwrites it before first use, preserving the "edit Δt and continue"
  decision. A debug assert that stored `time` matches the mesh timestamp at save remains worth
  adding.
- **O3 — Load-order fix. RESOLVED 2026-07-01 → option (i), implemented (Christian).**
  `hphirun.cpp` reordered: `set_circuit` before `load_memdump` (closes D1). The suggested
  `BELFEM_ERROR` for a `circuit` group with no circuit attached was not added — optional; it only
  protects future drivers from repeating the mis-ordering.
- **O4 — Length-mismatch policy. RESOLVED 2026-07-01 → hard error (implemented, R4).**
  `BELFEM_ERROR` with the message telling the user to delete `memdump.hdf5` to accept a cold
  circuit deliberately. Note the guard is weaker than the superseded fingerprint — same-length
  edits (e.g. swapping a resistor value) pass silently; that is consistent with the "input is the
  source of truth, state is just the last solution" philosophy.

---

## 6. `/circuit` Group Sketch (restart file, not `.bfm`)

Written by rank 0 into `memdump.hdf5` alongside `meta`/`fields`/`globals`.

```
/circuit                        (opt; only when a circuit is attached)
    time         scalar   real    circuit clock mTime at dump (O2)
    delta_time   scalar   real    Δt of the last completed step — self-description for prev_x;
                                  overwritten by set_timestep(input Δt) before first use (O2)
    x            [n]      real    solution vector (node voltages + unknown branch currents)
    prev_x       [n]      real
```

n = `mNumberOfNodes-1 + mNumberOfUnknownCurrents`; the load-side length check is the R4 guard.
When the format ships, document it in `src/circuit/doc/circuit_usage_guide.md` (or a sibling
format note) and mark this section historical.

---

## 7. Definition-of-Done Checklist

- [x] Every §3 (c)-row lands in a step; (a)-rows demonstrably rebuilt by R1 verification. R6 keeps
  the end-to-end check open.
- [x] Open questions O2-O4 resolved in place, none silently.
- [x] Old memdumps without `/circuit` still load via the `group_exists` guard. The cold-circuit
  warning is now optional/nice-to-have, not a blocker.
- [ ] R6 passes: restarted CORC coupled run matches uninterrupted within tolerance, including the
  source-crossing / fired-switch / first-step-reject variants.
- [ ] The v1 limitation (standalone L/C state) is documented in the module doc alongside the
  format; this plan moved to `todo/closed/` with a DONE Status.

---

## 8. Audit Trail

- **Tri-AI implementation audit done 2026-07-01** (thread: `tmp/ai_exchange/restart_circuit_state.md`,
  distilled here + devlog, GC-eligible). **Codex:** pass on all six audit questions; confirmed
  Claude's pre-finding that `cl_Circuit.hpp` was not self-contained for `hid_t` (fallback typedef
  lives in `hdf5_types.hpp:15-21`) — **fixed** by a direct include; flagged ignored `tStatus` in
  `load_state`/`save_state` as a release-mode gap (assert-only status in `hdf5_tools.hpp:563-567`)
  — **deferred** per the meshfile-plan D6 precedent (low value for a self-produced restart file;
  the length guard catches truncation-shaped corruption). **Grok (third voice):** independently
  confirmed the three highest-risk items — refactor equivalence (~98%), one-entry BDF degradation
  (~95%), HDF5 active-group semantics (~97%); zero refutations; corrected a Claude citation (the
  second `update_components()` call site is `solve()` at `cl_ElectricalCircuit.cpp:862`, not
  "compute_correction"). All cited lines re-verified by Claude before acceptance.
- Codex prose pass on the plan draft: done 2026-07-01.
- Exploration: Claude subagent sweep of `src/circuit` + Controller coupling, 2026-07-01; the
  load-bearing citations (circuit members, `shift`/`shift_back`, `Inductor`/`FEMTwoTerminals`
  `shift` bodies, memdump wiring, hphirun ordering, source/switch `shift`) were re-verified by
  direct read the same day.
- Design simplification to `{time, x, prev_x}`: decided 2026-07-01 by Christian ("the topology of
  the circuit is not altered, unlike with the mesh — we can fully reconstruct it at runtime");
  Claude verified the enabling mechanism (`shift()` re-seeds registers from current component
  state) and identified the standalone-inductor caveat + the clock as the one extra scalar.

---

## Related

- Parent (closed): `todo/closed/meshfile_refactor_plan.md` (R6/R6b, §6.4).
- Restart file: `Controller::save_memdump`/`load_memdump` (`cl_FEM_Controller.cpp:1781/1825`);
  mesh side `Mesh::save_meta/save_fields/save_globals` (`cl_Mesh.cpp:2971+`).
- Circuit module guide: `src/circuit/doc/circuit_usage_guide.md`.
- Also still open in the parent (kept there, not here): the `/meta/format_version` root attribute
  (O9 / §6.4) — a reader version gate; add before the `.bfm` format is frozen.
