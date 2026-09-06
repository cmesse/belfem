# Devlog 2026-07-17 - Calculator qold Cache Audit

**Date:** 2026-07-17
**Topic:** Independent audit of working-tree `Calculator::qold` cache acceleration
**AIs involved:** Claude, Codex
**Claude Confidence:** high on F1, medium-high on F2, medium on F3-F5
**Codex Audit Confidence:** high for F1/F5/F6 and added lifecycle blocker; medium-high for F2/F3; medium for F4
**Literature References:** N/A

## Summary

Codex audited Claude's `tmp/ai_exchange/qold_field_cache.md` findings against the working-tree diff in `cl_FEM_Calculator`, `cl_IWG`, and `cl_IWG_Timestep`, plus `Mesh::create_field()` and field ownership. The Codex section in the exchange file was filled with line-cited findings.

## Key Findings

- Blocking: the qold cache is initialized only from `Calculator::link(Group*)` (`src/fem/kernel/cl_FEM_Calculator.cpp:1167-1191`), while ordinary block assembly reaches `qold()` through `IWG::link_to_group()` and `Calculator::link(Element*)` without building `mQold` (`src/fem/kernel/cl_FEM_DofManager.cpp:635-659`, `src/fem/iwg/cl_IWG.cpp:383-388`, `src/fem/iwg/cl_IWG_TransientHeatConduction.cpp:56-62`).
- Confirmed F1: the cache key stride uses the inclusive max field index (`src/fem/kernel/cl_FEM_Calculator.cpp:2265-2297`), so `(s, max)` aliases `(s+1, 0)` when dof fields include both endpoints.
- Refined F2: standard Controller paths call `shift_fields()` before assembly and BDF startup caps `mOrderActive`, so first reads are ordered correctly in normal timestep flows; custom direct assembly without a prior shift remains risky.
- Confirmed F6 safe: cached `Vector<real>*` pointers to `mesh::Field::mData` remain stable across later `Mesh::create_field()` calls because fields are heap objects stored by pointer (`src/mesh/cl_Mesh.cpp:520-547`, `src/mesh/cl_Mesh_Field.hpp:33-43`).

## Changes Made / Proposed

- Updated `tmp/ai_exchange/qold_field_cache.md` under `## CODEX Audit`.
- No source code was modified.
- Proposed direction: rebuild the cache from a deliberate timestep-storage lifecycle hook and use a `mMaxDofFieldIndex + 1` stride or a pre-sized direct table with null checks.

## Open Questions

- Whether `set_timestepping_method()` after a cache build is supported; if so, it must also invalidate/rebuild `Calculator::mQold`.

## Files Updated

- `tmp/ai_exchange/qold_field_cache.md`
- `devlog/dl20260717_qold_cache_audit.md`
- `devlog/README.md`
