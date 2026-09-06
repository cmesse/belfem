# Devlog 2026-05-15 - Periodic Input Syntax Trace

**Date:** 2026-05-15
**Topic:** Existing input syntax for defining periodic front/back sidesets
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Traced `cmake-build-debug/input.conf`, `MaxwellFactory::create_periodic()`, and `PeriodicityFactory` to identify the currently implemented BELFEM syntax for periodicity. The existing hook is plane-based and lives under `topology -> periodic`; it does not directly accept sideset IDs. Updated the local input file with the working plane-based syntax and drafted a plan for adding explicit sideset filters.

## Key Findings

- `MaxwellFactory` looks for exactly one unlabeled `periodic` subsection under `topology`.
- The periodic block requires `source` and `target` keys, each with three IDs passed to `PeriodicityFactory::set_master_plane()` and `set_slave_plane()`.
- In `cmake-build-debug/corc.msh`, physical sideset 7 lies on `z = 0`, and physical sideset 8 lies on `z = 17.39482682`.
- Nodes/vertices `1,2,3` define the sideset-7 plane; `5,6,7` define the corresponding sideset-8 plane.

## Changes Made / Proposed

- Updated `cmake-build-debug/input.conf` with the currently supported periodic block:

```txt
topology
{
    periodic
    {
        source : 1, 2, 3 ;
        target : 5, 6, 7 ;
    }
}
```

- Added `todo/periodic_input_extension_plan.md` for the proposed sideset-filter extension.
- Updated `todo/periodic_bc_fix_plan.md` to link the refined input plan.

## Open Questions

- Whether BELFEM should add a compact `sidesets : 7 @ 8 ;` parser helper after the safer `source sidesets` / `target sidesets` filter path.

## Files Updated

- `cmake-build-debug/input.conf`
- `todo/periodic_input_extension_plan.md`
- `todo/periodic_bc_fix_plan.md`
- `devlog/dl20260515_periodic_input_syntax_trace.md`
- `devlog/README.md`
