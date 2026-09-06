# Devlog 2026-07-01 - Restart Circuit State Re-check

**Date:** 2026-07-01
**Topic:** Re-check and refresh `todo/restart_circuit_state.md`
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Re-read `todo/restart_circuit_state.md` against the current circuit restart implementation and refreshed stale todo/index wording. No source code was modified.

## Key Findings

- The core circuit restart implementation is present: `ElectricalCircuit::save_state/load_state` stores `{time, delta_time, x, prev_x}` and scatters `x` with `update_components()`; `Controller::save_memdump/load_memdump` writes/reads the optional `/circuit` group; `hphirun` and `hphiTrun` attach the circuit before loading the memdump.
- The todo file was partly stale: §1 still read as if no circuit serialization existed, R3 was half-open, and `todo/README.md` still described the work as a pending unaudited plan.
- Remaining items still look real: R5 (`IWG_Timestep` restart/old-history check), R6 (end-to-end coupled restart test), and circuit module documentation of the `/circuit` schema plus the standalone L/C v1 limitation. Static evidence suggests mesh old-step fields may already be saved and loaded because `Mesh::save_fields/load_fields` covers all mesh fields, but the behavior has not been validated end to end.

## Changes Made / Proposed

- Updated `todo/restart_circuit_state.md` to mark the old cold-restart failure section as historical, note the current implementation status, close R3, and tick completed definition-of-done items.
- Updated `todo/README.md` to reflect implemented/audited status and the actual `/circuit` contents `{time, delta_time, x, prev_x}`.

## Open Questions

- Does the current all-field memdump path already close R5 for `IWG_Timestep`, or should the old-field/timestep-history state gain explicit restart validation or versioned persistence before the coupled restart test?

## Files Updated

- `todo/restart_circuit_state.md`
- `todo/README.md`
- `devlog/dl20260701_restart_circuit_state_recheck.md`
- `devlog/README.md`
- `tmp/ai_exchange/restart_circuit_state_recheck.md`
