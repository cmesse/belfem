# Devlog 2026-04-23 — InterfaceTsConnector H1-H3 Audit

**Date:** 2026-04-23
**Topic:** Independent Codex audit of `InterfaceTsConnector` assembly hypotheses H1-H3 (plus H5 spot-check)
**AIs involved:** Claude, Codex
**Claude Confidence:** medium (~40% on H1, ~30% on H2/H3)
**Codex Audit Confidence:** high on H1/H3, medium-high on H2-for-this-case, high on H5

## Summary

Read-only source audit confirms the connector penalty path is dispatched and assembled in the live Jacobian path. H1-H3 do not explain the reported `~10^7` gap for the current linear `PENTA6TS` ↔ `HEX8TS` setup.

## Key Findings

- H1 (assembly reachability) is clean: `InterfaceTsConnector -> h_penalty` dispatch and sideset assembly loop are active in the nonlinear solve path (`src/fem/maxwell/cl_IWG_Maxwell.cpp`, `src/fem/kernel/cl_FEM_DofManager.cpp`, `src/fem/iwg/cl_IWG_Timestep.cpp`).
- H2 is a latent design bug but not active for this case: `IWG::number_of_dofs_per_element( SideSet* )` uses node-count logic for `Ghost`/`InterfaceTsConnector` (`src/fem/iwg/cl_IWG.cpp`), yet linear `PENTA6TS`(6) + `HEX8TS`(8) coincidentally matches edge DOF counts and current hardcoded workspace sizes in `IWG_Maxwell` (`src/fem/maxwell/cl_IWG_Maxwell.cpp`).
- H3 (ordering) is clean for current model: sideset DOFs are edge-only and linked in master-then-slave order (`src/fem/maxwell/cl_Maxwell_FieldList.cpp`, `src/fem/kernel/cl_FEM_Element.cpp`).
- H5 spot-check is clean: no `set_sources(...)` assignment found for `OuterEdgesLongitudinal`; source linking is on inner-duplicate longitudinal edges only (`src/mesh/cl_ThinShellFactory.cpp`).

## Changes Made / Proposed

- Reviewed the existing Codex verdict entry in `todo/ai_exchange.md` and cross-validated it against an independent source trace performed in this session.
- No source-code changes proposed in this session (read-only investigation).

## Open Questions

- Highest-value next checks remain runtime/geometry-focused: side-facet orientation correctness under all wrap cases, connector target-edge selection in layer views, and whether enforced trace quantity matches intended physics.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260423_interface_ts_connector_h123_audit.md