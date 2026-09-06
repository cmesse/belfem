# Devlog 2026-04-07 — Thin-Shell `bn` Projection Audit

**Date:** 2026-04-07
**Topic:** Read-only verification of Claude's claim that thin-shell kernels in `mt_maxwell_h.cpp` double-count tangential field by using the full air-side gradient as `bn`
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Verified Claude's claim directly in code. The thin-shell kernels currently compute

`bn = -0.5 * mu0 * ( Bm * phi_m + Bs * phi_s )`

where `Bm` and `Bs` are the full physical-space gradients of the scalar `phi` fields on the master and slave air elements. That means `bn` contains both normal and tangential components of the air-side field, while `bt = -mu0 * E * edge_q` already supplies the shell tangential field. So `b = bt + bn` double-counts the tangential component.

Projecting `bn` onto the facet normal is therefore the correct fix.

This is a real bug for HTS material evaluation, but it is not a convincing explanation for the Cu/Ag patch-test failure if the metal resistivity is effectively independent of `B`.

## Key Findings

- The affected pattern appears at all listed thin-shell sites in [mt_maxwell_h.cpp](../src/fem/maxwell/matrices/mt_maxwell_h.cpp):
  - `bn = tCalc->Bm(0) * phi_m + tCalc->Bs(0) * phi_s;`
  - `bn *= -0.5 * constant::mu0;`

- `Calculator::Bscalar_master()` and `Calculator::Bscalar_slave()` in [cl_FEM_Calculator.hpp](../src/fem/kernel/cl_FEM_Calculator.hpp) compute
  - `mInvJm * dNdXi`
  - `mInvJs * dNdXi`
  respectively, so they return full gradients, not normal projections.

- The shell field `bt = -mu0 * E * edge_q` comes from `EF_PENTA6TS` in [cl_EF_PENTA6TS.cpp](../src/fem/interpolation/nedelec/cl_EF_PENTA6TS.cpp), which is built from the surface Jacobian pseudoinverse and in-surface gradients. That supports the intended split:
  - `bt`: tangential shell field
  - `bn`: normal complement from neighboring `phi`

- Calling `tCalc->normal()` in the metal shell kernels should be safe in the current path:
  - `get_normal_calculator()` returns the visible shell sideset calculator in [cl_FEM_Calculator.cpp](../src/fem/kernel/cl_FEM_Calculator.cpp)
  - `Calculator::allocate()` assigns `mFunNormal` for the supported master geometries in the same file
  - HTS shell kernels already call `tCalc->normal()` on the same calculator path

## Changes Made / Proposed

- No source-code changes were made in this audit.
- Added a Codex audit entry to `todo/ai_exchange.md`.

## Open Questions

- None for the correctness of the projection itself.
- The larger unresolved issue remains the patch-test failure, which this fix is unlikely to resolve for constant-resistivity Cu/Ag stacks.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260407_bn_projection_audit.md
