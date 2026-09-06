# Devlog 2026-04-02 — h_ghost Resistivity Strategy

**Date:** 2026-04-02
**Topic:** Read-only recommendation on whether `h_ghost()` should recover resistivity from state or store it in a mesh field
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** `todo/thinshell_selective_nitsche_coupling.md`

## Summary

Reviewed whether the ghost interface kernel should use a stored mesh field for resistivity recovery.

Recommendation: do not introduce a mesh field as the primary mechanism yet. The interface kernel already has access to the current local master/slave DOF state, while a stored resistivity field would need continuous updates and would still be a poor proxy for HTS materials whose effective resistivity is state-dependent.

## Key Findings

- `Calculator::q()` already reconstructs the current local DOF vector for the active element from the mesh fields at assembly time ([src/fem/kernel/cl_FEM_Calculator.cpp](/home/christian/codes/belfem/src/fem/kernel/cl_FEM_Calculator.cpp#L1822)).
- Master/slave-coupled local DOFs are linked in master-first, slave-second order, which is compatible with splitting the ghost state into the two shell sides directly inside `h_ghost()` ([src/fem/kernel/cl_FEM_Element.cpp](/home/christian/codes/belfem/src/fem/kernel/cl_FEM_Element.cpp#L829)).
- A persistent resistivity field would need to be refreshed every nonlinear iteration / timestep update and would still collapse quadrature/state dependence into a coarse stored value, which is especially questionable for HTS materials.
- A local on-the-fly coefficient evaluation is better aligned with the current kernel structure and with the fact that ghost stabilization is only assembled on a subset of interfaces.

## Changes Made / Proposed

- No source changes made.
- Logged the recommendation in `todo/ai_exchange.md`.

## Open Questions

- For the first practical on-the-fly implementation, should the interface coefficient use a crude `rho(T)` proxy or a state-dependent HTS proxy based on local current magnitude?
- If caching becomes necessary later, should it be done per element, per interface, or per block/material region?

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260402_h_ghost_resistivity_strategy.md
