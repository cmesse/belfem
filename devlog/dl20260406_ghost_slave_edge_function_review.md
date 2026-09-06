# Devlog 2026-04-06 — Ghost Slave Edge Function Review

**Date:** 2026-04-06
**Topic:** Read-only investigation of `Slave Edge function has not been assigned` in the ghost thin-shell path
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Investigated the new `Calculator::Es()` abort during `h_ghost()` assembly. The failure is caused by a stale calculator assumption: slave edge basis functions are only allocated for `"edge_a"`, while ghost sidesets carry conductor `"edge_h"` DOFs on both master and slave shell elements.

## Key Findings

- `h_ghost()` explicitly requests both master and slave edge basis functions through `Em()` and `Es()` (`src/fem/maxwell/matrices/mt_maxwell_h.cpp:2124`).
- `Calculator::Es()` aborts when `mEdgeFunctionSlave` is null (`src/fem/kernel/cl_FEM_Calculator.cpp:2103`).
- In `Calculator::allocate_memory()`, `mEdgeFunctionsSlave` is allocated only if `"edge_a"` is present in `all_fields()` (`src/fem/kernel/cl_FEM_Calculator.cpp:219`, `src/fem/kernel/cl_FEM_Calculator.cpp:313`).
- Ghost sidesets inherit their DOF table from `Conductor`, so they carry H-formulation edge/face DOFs rather than `"edge_a"` (`src/fem/maxwell/cl_Maxwell_FieldList.cpp:116`, `src/fem/maxwell/cl_Maxwell_FieldList.cpp:410`).
- `Element::link_dofs_master_and_slave()` already counts and links edge DOFs on both master and slave, confirming the ghost sideset is H-H edge/edge coupling (`src/fem/kernel/cl_FEM_Element.cpp:791`).
- `Calculator::link( Element * )` would assign `mEdgeFunctionSlave` if the slave edge-function pool existed, so the root problem is missing allocation rather than the final link step (`src/fem/kernel/cl_FEM_Calculator.cpp:984`).

## Changes Made / Proposed

- No source changes made.
- Appended the investigation result to `todo/ai_exchange.md`.
- Proposed likely fix directions:
  - allocate slave edge functions based on the actual linked sideset DOF layout rather than the presence of `"edge_a"` in `all_fields()`, or
  - special-case `DomainType::Ghost` so the slave edge-function pool is created for ghost H-H sidesets under the `edge_h` path as well.

## Open Questions

- Whether the same stale `"A is always slave"` assumption also affects slave-side face basis allocation for higher-order ghost elements.
- Whether `SideSetDofLinkMode::MasterAndSlave` still carries other hidden H-phi assumptions that should be split from the new ghost H-H path.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260406_ghost_slave_edge_function_review.md
