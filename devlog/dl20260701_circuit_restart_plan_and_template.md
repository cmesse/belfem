# Devlog 2026-07-01 — Circuit-Restart Plan + Todo-Plan Template

**Date:** 2026-07-01
**Topic:** Expanded `todo/restart_circuit_state.md` into a full implementation plan (minimal-state
design) and extracted a reusable plan template from the closed meshfile refactor plan.
**AIs involved:** Claude (exploration + drafts), Codex (prose pass on both files)
**Claude Confidence:** high on the code findings (all load-bearing citations re-verified by direct
read); the plan itself is unaudited (Codex technical audit pending).

## Summary

Two documentation deliverables, no source modified:

1. **`todo/plan_template.md`** — canonical structure for substantial todo plans, extracted from
   `todo/closed/meshfile_refactor_plan.md`: header with living Status line, scope guards, failure
   table, gap table with (a)/(b)/(c) classes, ordered `Rn` steps + `Dn` defect tracker (severity,
   "Fixed YYYY-MM-DD", FALSE-POSITIVE retractions kept in place), `On` open questions,
   definition-of-done, audit trail, decision appendices, lifecycle. Registered in
   `todo/README.md`; Claude memory updated so future plans use it automatically.

2. **`todo/restart_circuit_state.md`** — rewritten from a stub into a template-form plan for
   persisting the coupled circuit's dynamic state in `memdump.hdf5` (parent R6b follow-up).

## Key Findings (circuit exploration, verified 2026-07-01)

- **The circuit clock is the worst cold-restart defect:** `ElectricalCircuit::mTime` is private
  with no setter (`cl_ElectricalCircuit.hpp:92`), advanced only by `shift()`
  (`cl_ElectricalCircuit.cpp:180`), and time-dependent sources evaluate from it
  (`cl_VoltageSource.cpp:51-54`). A cold restart replays the drive waveform from t≈0 for the
  entire remainder of the run — not a transient.
- **Load-order hazard:** `hphirun.cpp:78` calls `load_memdump` *before* `set_circuit` (`:81`) —
  the Controller has no circuit to restore into at load time (plan O3).
- **`shift()` re-seeds registers from current component state** (`Inductor::shift` pushes
  `get_current()`, `cl_Inductor.cpp:49`; `FEMTwoTerminals::shift` pushes `mIn`/`mVn`,
  `cl_FEMTwoTerminals.cpp:36-39`) — the mechanism that makes the minimal-state design work.
- **Switch self-corrects given a restored clock:** rebuilt `mIsClosed` is the input initial state
  and the first `shift()` at T ≥ `mSwitchTime` toggles it immediately (`cl_Switch.cpp:43-50`).
- **Standalone inductor current is NOT recoverable from `x`:** the Norton companion hides it in
  the history-dependent `miLh` (`cl_Inductor.cpp:35-38`). Accepted v1 limitation; irrelevant for
  CORC-class circuits (sources + resistors + terminal pairs).
- Circuit solve is rank-0 only (`cl_FEM_Controller.cpp:147`); memdump already rank-0-guarded.

## Decisions (Christian, 2026-07-01)

- **Minimal-state design:** persist only `{time, x, prev_x}` in a `/circuit` group — no
  ShiftRegister serialization, no per-component `save_state` virtuals, no fingerprint keying
  (only an `x`-length guard). Rationale: circuit topology is input-defined and never enriched,
  so — unlike the mesh — everything else reconstructs at runtime; registers re-seed order-1 via
  the first `shift()`, matching the FEM-side "state restart, cold-started integrator" philosophy.
  This superseded the initial draft's per-component register-persistence design the same day.

## Changes Made

- `todo/plan_template.md` — new (+ Codex prose pass applied).
- `todo/restart_circuit_state.md` — rewritten as the full plan (R1–R6, O2–O4; O1 resolved moot).
- `todo/README.md` — template convention note in Layout; restart entry updated.
- Claude memory: template convention saved.

## Same-Day Implementation + Tri-AI Audit (updated later on 2026-07-01)

Christian implemented the minimal-state restart the same evening; Claude reviewed each increment
and added the R4 guard; Codex + Grok audited the final diff. All plan steps except R5/R6 closed.

- **Christian:** `Circuit::save_state`/`load_state` pure virtuals; `ElectricalCircuit`
  implementation storing `{time, delta_time, x, prev_x}` in a `/circuit` group of
  `memdump.hdf5`; `update_components()` factored out of `shift_back()`/`solve()` and called from
  `load_state` (the scatter); Controller `save_memdump`/`load_memdump` wiring; hphirun reorder
  (`set_circuit` before `load_memdump` — closed D1, the silent-no-op load).
- **Claude review catches during implementation:** D1 (load was a silent no-op before the
  reorder); the missing scatter (verified load-bearing: `Capacitor::shift` reads node voltages,
  `cl_Capacitor.cpp:51`); `delta_time`/`prev_x` are self-description, not live state (both
  overwritten by `set_timestep` + `shift()` before any consumer — O2 analysis).
- **Claude edits (user-approved):** R4 length guard (`BELFEM_ERROR` in `load_state`, O4 hard-error
  policy); `cl_Circuit.hpp` now includes `hdf5_types.hpp` (audit fix, `hid_t` self-containment).
- **Audit (Codex full 6-question pass + Grok third-voice on the top 3, zero refutations):**
  refactor is a pure extraction; one-entry BDF register → exact backward-Euler
  (`cl_BDF.cpp:33`, closes R1); HDF5 group semantics + rank-0-only restore clean; no collective
  imbalance. Deferred residual: release-mode `tStatus` checks (assert-only in
  `hdf5_tools.hpp:563-567`) per the meshfile D6 precedent. Thread
  `tmp/ai_exchange/restart_circuit_state.md` distilled + GC-eligible.

## Open Questions

- ~~Plan O2/O3/O4~~ — all resolved 2026-07-01 (O2 store the clock + `delta_time`; O3 hphirun
  reorder, option (i); O4 hard error on length mismatch).
- Remaining before the plan closes: **R5** (confirm `IWG_Timestep` cold-starts cleanly on a
  load) and **R6** (end-to-end CORC coupled restart vs uninterrupted run — Christian, after
  `make reset && make hphirun`).
