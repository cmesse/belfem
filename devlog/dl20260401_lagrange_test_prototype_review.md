# Devlog 2026-04-01 — Lagrange Test Prototype Review

**Date:** 2026-04-01
**Topic:** Read-only audit of `src/fem/interpolation/test_langrange.cpp`
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Reviewed the standalone Lagrange interpolation test prototype and wrote the findings to `tmp/test_lagrange_review_codex.md`. The current prototype is not trustworthy as a validation harness because the assertion macro is broken and the derivative/Hessian checks use incorrect row indexing.

## Key Findings

- `TEST_NEAR` in `src/fem/interpolation/test_langrange.cpp:45-50` compares against `aEpsilon` instead of `aExpect` and inverts the pass/fail logic.
- The first-derivative check uses `N_xi.row( k )` at `src/fem/interpolation/test_langrange.cpp:184`, where `k` is the quadrature-point index rather than the derivative-direction index.
- The second-derivative checks repeat the same indexing mistake at `src/fem/interpolation/test_langrange.cpp:280` and `src/fem/interpolation/test_langrange.cpp:361`.
- The 3D Hessian branch additionally rescales the wrong temporary and fills the `d2/dzeta2` slot from the wrong row at `src/fem/interpolation/test_langrange.cpp:343-344`.
- The prototype omits `TET20`, `TET35`, and `HEX64`, although the factory supports them.
- `src/fem/interpolation/CMakeLists.txt:40-41` wires `test_langrange.cpp` to the executable name `test_facets`, so there is no dedicated `test_langrange` target.

## Changes Made / Proposed

- Added `tmp/test_lagrange_review_codex.md`
- Appended the audit summary to `todo/ai_exchange.md`

## Open Questions

- Whether the intention is to keep this file as a temporary standalone executable or replace it directly with a GoogleTest target.
- Whether the current `test_facets` executable name in `src/fem/interpolation/CMakeLists.txt` is an accidental leftover or an intentional temporary reuse.

## Files Updated

- tmp/test_lagrange_review_codex.md
- todo/ai_exchange.md
