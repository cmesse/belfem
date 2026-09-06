# Devlog 2026-08-30 — DR-110 Capacity Contract

**Date:** 2026-08-30
**Topic:** Retire the `Basis` source-capacity question by design ruling
**AIs involved:** Codex, Christian
**Codex Audit Confidence:** high
**Verification:** static source trace plus Christian's design ruling; no executable gate required

## Summary

DR-110 is closed as accepted caller-discipline behavior. `Basis` does not retain a source-allocation capacity; callers allocate exactly the source count they append.

## Changes Made

- Documented the allocation/append invariant on the public `Basis` declarations and beside the allocation implementation.
- Moved DR-110 to the closed register with the ruling and future-caller constraint.

## Files Updated

- `src/mesh/cl_Mesh_Basis.hpp`
- `src/mesh/cl_Mesh_Basis.cpp`
- `todo/debt_register.md`
- `todo/debt_register_closed.md`
- `devlog/README.md`
