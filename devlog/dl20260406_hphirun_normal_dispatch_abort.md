# Devlog 2026-04-06 — hphirun Normal Dispatch Abort

**Date:** 2026-04-06
**Topic:** Read-only investigation of `hphirun` aborting with `No normal function assigned`
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Investigated the reported `hphirun` crash during thin-shell Nitsche-DG-Ghost setup. The abort is triggered in `Calculator::allocate()` because ghost thin-shell sidesets carry `PENTA6TS` master/slave element types, but the calculator's normal-function dispatch still lacks a `GeometryType::PENTA` branch.

## Key Findings

- `DofManager::init_work()` allocates calculators for all selected sidesets, including ghost sidesets (`src/fem/kernel/cl_FEM_DofManager.cpp:305`).
- Maxwell selects `DomainType::Ghost` sidesets and routes them to `h_ghost()` (`src/fem/maxwell/cl_MaxwellFactory.cpp:1707`, `src/fem/maxwell/cl_IWG_Maxwell.cpp:556`).
- `ThinShellFactory` builds ghost facets with shell master/slave elements, specifically face `1` on the master and face `0` on the slave (`src/mesh/cl_ThinShellFactory.cpp:1696`).
- `SideSet` stores `mMasterType` / `mSlaveType` from those shell elements, so ghost sidesets expose `PENTA6TS` to the calculator (`src/fem/kernel/cl_FEM_SideSet.cpp:41`).
- `Calculator::allocate()` supports `PENTA6TS` for Nedelec data but not for normal dispatch; `mesh::geometry_type( PENTA6TS )` returns `GeometryType::PENTA`, which falls into the `default` error path (`src/fem/kernel/cl_FEM_Calculator.cpp:357`, `src/fem/kernel/cl_FEM_Calculator.cpp:734`, `src/mesh/meshtools.cpp:352`).
- `h_ghost()` itself later asserts `PENTA6TS` master/slave elements, confirming this is the expected active path rather than corrupted setup (`src/fem/maxwell/matrices/mt_maxwell_h.cpp:2112`).

## Changes Made / Proposed

- No source changes made.
- Appended the investigation result to `todo/ai_exchange.md`.
- Proposed next source-level direction: add prism normal support to `Calculator` or special-case the thin-shell prism ghost path.

## Open Questions

- Whether the intended framework behavior is to support prism-master sidesets generally or only thin-shell TS prism facets.
- Whether `PENTA18TS` should share the same normal path once higher-order ghost support is enabled.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260406_hphirun_normal_dispatch_abort.md
