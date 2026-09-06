# Devlog 2026-04-06 — Ghost DOF Count Review

**Date:** 2026-04-06
**Topic:** Read-only investigation of the ghost-sideset `12 vs 6` local DOF mismatch
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Investigated the new `compute_element_dof_connectivity()` assert for a ghost sideset element reporting `12` local DOFs while `IWG::number_of_dofs_per_element( SideSet * )` expected `6`. The mismatch is real: the generic `IWG` sideset helper still assumes an H-phi-style mixed interface, while ghost H-H coupling links edge/face DOFs from both master and slave shell elements.

## Key Findings

- Maxwell uses `SideSetDofLinkMode::MasterAndSlave` for sidesets (`src/fem/maxwell/cl_IWG_Maxwell.cpp:42`).
- The ghost DOF table is built from `Conductor` DOFs (`src/fem/maxwell/cl_Maxwell_FieldList.cpp:116`), so ghost sidesets carry edge/face H-formulation DOFs.
- `Element::link_dofs_master_and_slave()` counts edge/face contributions from both master and slave sides (`src/fem/kernel/cl_FEM_Element.cpp:791`).
- `IWG::number_of_dofs_per_element( SideSet * )` still counts edge/face DOFs only on the master side in its generic edge-dof branch (`src/fem/iwg/cl_IWG.cpp:672`).
- `DofManager` uses this undercounted value to size temporary Jacobian/RHS buffers for sideset assembly (`src/fem/kernel/cl_FEM_DofManager.cpp:698`), so the problem is not limited to the debug assert.

## Changes Made / Proposed

- No source changes made.
- Appended the investigation result to `todo/ai_exchange.md`.
- Proposed likely fix directions:
  - special-case `DomainType::Ghost` in `IWG::number_of_dofs_per_element( SideSet * )`, or
  - replace the asymmetric generic count with one that reflects the actual linked master/slave element layout.

## Open Questions

- Whether any non-ghost Maxwell sidesets also use `MasterAndSlave` with edge/face DOFs on both sides and therefore share the same latent undercount.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260406_ghost_dof_count_review.md
