# Devlog 2026-05-14 - DofData Periodic Step 2 Audit

**Date:** 2026-05-14
**Topic:** Read-only audit of latest DofData periodic DOF changes
**AIs involved:** Codex
**Codex Audit Confidence:** high for sequencing/build assessment; medium for target-already-hanging runtime invariant until tested on a periodic thin-shell case
**Literature References:** N/A

## Scope

Audited the latest changes in:

- `src/fem/kernel/cl_FEM_DofMgr_DofData.cpp`
- `src/fem/kernel/cl_FEM_DofMgr_DofData.hpp`
- `src/fem/kernel/cl_FEM_Dof.cpp`
- `src/fem/kernel/cl_FEM_Dof.hpp`

## Findings

The new `collect_hanging_dofs()` flow now rebuilds hanging-DOF membership from all `mDOFs` after periodic DOF constraints are applied. This addresses the main Step 2 sequencing gap: periodic DOFs newly marked as hanging can now enter `mHangingDOFs` before serial removal or MPI exchange.

The smaller-ID periodic direction is implemented through `DofData::entangle()`, and `Dof::entangle()` unfolds an already-hanging representative by copying its sources and weights. This matches the intended source-chain behavior.

The remaining correctness risk is the case where the selected periodic target already has sources. `Dof::entangle()` asserts that the target has no sources. That is acceptable only if documented as an invariant; otherwise the implementation needs explicit replacement/merge semantics.

## Clarity Recommendations

- Rename `entangle()` to something more specific, such as `inherit_sources_from()` on `Dof` and `constrain_periodic_pair_by_id()` on `DofData`.
- Add a `type_id()` equality assertion in `DofData::entangle()`.
- Prefer stack `DynamicBitset` objects over `new`/`delete` in `collect_hanging_dofs()`.
- Consider reusing `set_sources()` inside `Dof::entangle()` to keep source validation in one place.

## Verification

- `make -C cmake-build-debug hphirun -j2` passed.

## Files Touched

- `devlog/dl20260514_dofdata_periodic_step2_audit.md`
- `devlog/README.md`

