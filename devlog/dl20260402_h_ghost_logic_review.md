# Devlog 2026-04-02 — h_ghost Logic Review

**Date:** 2026-04-02
**Topic:** Read-only audit of the current `h_ghost()` logic after the latest master/slave convention clarification
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** `todo/thinshell_selective_nitsche_coupling.md`

## Summary

Reviewed `src/fem/maxwell/matrices/mt_maxwell_h.cpp:h_ghost()` using the clarified convention that the master element is the lower shell element and the slave element is the upper shell element.

The current face-slot mapping is coherent with that convention: the interface is the top face of the master element and the bottom face of the slave element. The remaining logic defects are an inverted denominator assertion, a reversed sign in the slave-side derivative operator, and the still-present HTS-to-zero coefficient guard.

## Key Findings

- `BELFEM_ASSERT( rho_m + rho_s < BELFEM_EPSILON, ... )` is inverted and therefore aborts on normal interfaces in debug builds while leaving a divide-by-zero path in release builds ([src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L1835), [src/core/assert.hpp](/home/christian/codes/belfem/src/core/assert.hpp#L150)).
- The master/slave face mapping is coherent: ghost facets use master face `1` and slave face `0`, which means `Em` is evaluated on the master top face and `Es` on the slave bottom face ([src/mesh/cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L1696), [src/fem/interpolation/fn_IF_initialize_integration_points_on_facet.cpp](/home/christian/codes/belfem/src/fem/interpolation/fn_IF_initialize_integration_points_on_facet.cpp#L1181)).
- `EF_PENTA6TS::E()` activates columns `3..5` on the top face and `0..2` on the bottom face, so the current use of `Em(:,3..5)` and `Es(:,0..2)` is consistent with the active interface basis ([src/fem/interpolation/nedelec/cl_EF_PENTA6TS.cpp](/home/christian/codes/belfem/src/fem/interpolation/nedelec/cl_EF_PENTA6TS.cpp#L282)).
- The slave derivative sign is still reversed: with local thin-shell DOF order `[bottom, top]`, a through-thickness derivative should use `[-E_bottom | +E_top]/h` on both elements, but `Ds` is assembled as `[+Es | -Es]/hs` ([src/mesh/cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L1258), [src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L1889)).
- The HTS coefficient path still forces any HTS side below `T_crit` to zero, which collapses the ghost term on metal/HTS interfaces and conflicts with the stated HTS material contract and the existing thin-shell kernels that use `rho_powerlaw()` / `rho_piecewise()` ([src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L1831), [src/physics/materials/cl_Material.hpp](/home/christian/codes/belfem/src/physics/materials/cl_Material.hpp#L517), [src/physics/materials/cl_Material_YBCO.hpp](/home/christian/codes/belfem/src/physics/materials/cl_Material_YBCO.hpp#L43)).

## Changes Made / Proposed

- No source changes made.
- Logged the audit in `todo/ai_exchange.md`.

## Open Questions

- Should `h_ghost()` use a crude positive proxy such as `rho(T)` on both sides, or should it evaluate an effective HTS resistivity consistent with the active material model?
- Is there any intended formulation in which `Ds` is defined with an outward-normal sign instead of the common through-thickness `z` derivative sign?

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260402_h_ghost_logic_review.md
