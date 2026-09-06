# Devlog 2026-04-02 — Thin-Shell Facet Integration Review

**Date:** 2026-04-02
**Topic:** Read-only audit of `fn_IF_initialize_integration_points_on_facet.cpp` for `QUAD*TS` / `PENTA*TS`
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Reviewed the thin-shell facet integration-point logic against the TS element definitions in `cl_Element_QUAD*.hpp` and `cl_Element_PENTA*.hpp`. The current code routes TS element types through the generic `QUAD` / `PENTA` facet logic, which does not match the actual TS face numbering or orientation conventions.

## Key Findings

- `initialize_integration_points_on_facet()` dispatches only by `mesh::geometry_type()` in [src/fem/interpolation/fn_IF_initialize_integration_points_on_facet.cpp](/home/christian/codes/belfem/src/fem/interpolation/fn_IF_initialize_integration_points_on_facet.cpp#L36) and [src/fem/interpolation/fn_IF_initialize_integration_points_on_facet.cpp](/home/christian/codes/belfem/src/fem/interpolation/fn_IF_initialize_integration_points_on_facet.cpp#L94), so `QUAD*TS` and `PENTA*TS` fall through to ordinary `QUAD` / `PENTA` handling.
- `PENTA6TS` face `0` is `[0,1,2]` in [src/mesh/cl_Element_PENTA6TS.hpp](/home/christian/codes/belfem/src/mesh/cl_Element_PENTA6TS.hpp#L56), while the corresponding volume bottom face is `PENTA6` face `3` = `[0,2,1]` in [src/mesh/cl_Element_PENTA6.hpp](/home/christian/codes/belfem/src/mesh/cl_Element_PENTA6.hpp#L76). The geometric face matches, but the orientation is reversed.
- `PENTA6TS` face `1` = `[3,4,5]` in [src/mesh/cl_Element_PENTA6TS.hpp](/home/christian/codes/belfem/src/mesh/cl_Element_PENTA6TS.hpp#L64) matches volume `PENTA6` face `4` = `[3,4,5]` in [src/mesh/cl_Element_PENTA6.hpp](/home/christian/codes/belfem/src/mesh/cl_Element_PENTA6.hpp#L84).
- `QUAD4TS` face `0` = `[0,1]` in [src/mesh/cl_Element_QUAD4TS.hpp](/home/christian/codes/belfem/src/mesh/cl_Element_QUAD4TS.hpp#L53) matches volume `QUAD4` face `0` in [src/mesh/cl_Element_QUAD4.hpp](/home/christian/codes/belfem/src/mesh/cl_Element_QUAD4.hpp#L52), but `QUAD4TS` face `1` = `[3,2]` in [src/mesh/cl_Element_QUAD4TS.hpp](/home/christian/codes/belfem/src/mesh/cl_Element_QUAD4TS.hpp#L59) is the reversed version of volume `QUAD4` face `2` = `[2,3]` in [src/mesh/cl_Element_QUAD4.hpp](/home/christian/codes/belfem/src/mesh/cl_Element_QUAD4.hpp#L64).
- `number_of_facets()` and `number_of_orientations()` are also inconsistent for TS types because they are still derived from generic geometry in [src/mesh/meshtools.cpp](/home/christian/codes/belfem/src/mesh/meshtools.cpp#L563) and [src/mesh/meshtools.cpp](/home/christian/codes/belfem/src/mesh/meshtools.cpp#L791).
- `Calculator::link( mesh::Facet * )` indexes slave integrations as `index_on_slave * 3 + orientation_on_slave - 1` in [src/fem/kernel/cl_FEM_Calculator.cpp](/home/christian/codes/belfem/src/fem/kernel/cl_FEM_Calculator.cpp#L1014), which aligns with triangular-face semantics, not the current generic `PENTA` precompute layout.

## Changes Made / Proposed

- No source changes made.
- Added collaboration notes to `todo/ai_exchange.md`.
- Added this devlog entry.

## Open Questions

- Whether TS facet integration should be handled by dedicated `ElementType::*TS` branches in `initialize_integration_points_on_facet()`, or by introducing TS-aware helpers below the existing geometry dispatch.
- Whether `number_of_facets()` / `number_of_orientations()` should also gain TS-specific overrides to match the actual shell element facet model.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260402_thinshell_facet_integration_review.md
