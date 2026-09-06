# Devlog 2026-03-30 — Facet Test Prototype Review

**Date:** 2026-03-30
**Topic:** Read-only audit of the facet prototype before converting it into gtests
**AIs involved:** Codex
**Codex Audit Confidence:** high

## Summary

Reviewed the prototype in [src/fem/interpolation/test_facets.cpp](/home/christian/codes/belfem/src/fem/interpolation/test_facets.cpp#L39) and ran the built executable `cmake-build-debug/bin/test_facets`.

The prototype is useful and largely aligned with the intended integration-test scope: it checks that facet quadrature weights match surface quadrature weights and that parent-element and facet-element mappings land on the same physical points. The main issue is not missing harness machinery, but that the prototype already exposes concrete higher-order slave-orientation mismatches.

## Key Findings

- The prototype covers master and slave facet mappings for `TET4/10`, `HEX8/20/27`, `PENTA6/15/18`, and `PYRA5/13/14` via [src/fem/interpolation/test_facets.cpp](/home/christian/codes/belfem/src/fem/interpolation/test_facets.cpp#L49).
- Its test strategy is consistent with the current plan split: [todo/tests/tests_11_interpolation.md](/home/christian/codes/belfem/todo/tests/tests_11_interpolation.md#L393) says `initialize_integration_points_on_facet` belongs in integration tests rather than interpolation metadata tests.
- All weight comparisons passed in the executable run.
- Point-mapping failures were observed only for higher-order slave orientations:
  - `HEX20`: facet 1 / printed orientation 2, facet 2 / printed orientation 2, facet 3 / printed orientations 1 and 3
  - `HEX27`: same failing pattern as `HEX20`
  - `PENTA15`: facet 2 / printed orientations 2 and 3
  - `PENTA18`: facet 2 / printed orientations 2 and 3
- The likely defect location is the slave helper logic in [src/fem/interpolation/fn_IF_initialize_integration_points_on_facet.cpp](/home/christian/codes/belfem/src/fem/interpolation/fn_IF_initialize_integration_points_on_facet.cpp#L615) and [src/fem/interpolation/fn_IF_initialize_integration_points_on_facet.cpp](/home/christian/codes/belfem/src/fem/interpolation/fn_IF_initialize_integration_points_on_facet.cpp#L893), not the basic prototype loop.
- The prototype does not cover 2D edge-facet behavior (`TRI*`, `QUAD*`) even though the API and [src/fem/interpolation/cl_IF_IntegrationData.cpp](/home/christian/codes/belfem/src/fem/interpolation/cl_IF_IntegrationData.cpp#L105) support those slave/master paths.

## Changes Made / Proposed

- No source-code changes made.
- Proposed next step: convert the prototype into gtests with explicit assertions and isolate the currently failing higher-order slave orientations as named regressions instead of treating the whole prototype as a blanket pass.

## Open Questions

- Should the immediate gtest file encode only the currently passing cases, or should it also include expected-failure/regression cases for the known higher-order slave mismatches?
- Do you want the new gtest file to cover only 3D volume facets, matching this prototype, or the full public API including 2D edge-facet mappings?

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260330_facet_test_prototype_review.md
