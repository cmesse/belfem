# Devlog 2026-04-07 — Patch-Test Trace Step 1 and Step 2

**Date:** 2026-04-07
**Topic:** Read-only follow-up on the layered thin-shell patch-test failure: LINE2 T-matrix ordering and `EF_PENTA6TS::E()` constant-field reproduction
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Completed two of the remaining read-only patch-trace checks from the thin-shell topology devlog.

- The LINE2 node-hanging path does not appear to have a sign bug tied to “lower global node ID at local facet position 1”.
- The shell edge interpolation operator `EF_PENTA6TS::E()` reproduces a uniform in-plane field exactly on a flat unit triangle.

These two checks shift suspicion away from the shell basis itself and toward the thin-shell interface attachment / coupling path.

## Key Findings

- In [cl_MaxwellFactory.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_MaxwellFactory.cpp), `hang_thinshell_edges_on_nodes_bottom()` and `hang_thinshell_edges_on_nodes_top()` populate hanging edge sources by local facet position, not by global node ID.
- In [cl_FEM_DofMgr_DofData.cpp](/home/christian/codes/belfem/src/fem/kernel/cl_FEM_DofMgr_DofData.cpp), LINE2 node-source expansion uses coefficients `[+1, -1]`, so the resulting hanging edge DOF is `phi(local node 0) - phi(local node 1)`.
- Therefore, if the lower global node ID sits at facet position `1`, the sign still follows the local edge orientation. That is not, by itself, a bug.
- In [cl_EF_PENTA6TS.cpp](/home/christian/codes/belfem/src/fem/interpolation/nedelec/cl_EF_PENTA6TS.cpp), the unit-triangle pseudoinverse gradients are:
  - `nabla_xi = (-1, -1, 0)`
  - `nabla_eta = ( 1,  0, 0)`
  - `nabla_zeta = ( 0,  1, 0)`
- With those gradients and the code’s own lower/upper face duplication, a constant field `H = (Hx, Hy, 0)` is reproduced exactly by the shell operator when using the oriented edge integrals `q0 = Hx`, `q1 = -Hx + Hy`, `q2 = -Hy` on both faces.
- So `EF_PENTA6TS::E()` passes the flat constant-field reproduction check and is unlikely to be the source of the `~20 dB` multilayer patch-test failure.

## Changes Made / Proposed

- No source-code changes were made.
- Added an audit note to `todo/ai_exchange.md`.

## Open Questions

- Does the remaining patch-test defect live in shell-air / shell-shell attachment rather than in the shell basis?
- Does a simplified run with uniform `H_inf` and without cuts/current BCs still fail, confirming the issue is internal to the layered coupling path?

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260407_patchtest_trace_step12.md
- devlog/README.md
