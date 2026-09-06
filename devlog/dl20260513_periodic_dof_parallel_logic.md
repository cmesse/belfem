# Devlog 2026-05-13 - Periodic DOF Parallel Logic Audit

**Date:** 2026-05-13
**Topic:** Audit of periodic DOF source direction and parallel communication requirements
**AIs involved:** Codex
**Codex Audit Confidence:** high for current DofData communication mechanics; medium-high for no extra mesh entity communication under the existing periodic distributor path
**Literature References:** N/A

## Scope

Read-only audit of the proposed rule: for each periodic DOF pair, hang the larger-ID DOF onto the smaller-ID DOF. Checked whether this requires communicating both periodic DOFs and their mesh entities to worker ranks.

## Findings

The smaller-ID rule is compatible with the existing hanging-DOF communication model, provided periodic DOFs are registered into `mHangingDOFs` after the periodic source relation is created.

Current parallel `DofData::collect_hanging_dofs()` already supports missing source DOFs on worker ranks:

- rank 0 sends the global hanging-DOF ID list,
- each worker selects the hanging DOFs it actually has,
- rank 0 sends the needed source DOF IDs for those local hanging DOFs,
- workers request missing source DOFs,
- rank 0 sends source DOF metadata: global index, type, entity type, index on entity, and mesh-basis ID,
- workers create placeholder source DOFs attached to local mesh entities.

Therefore, the periodic fix should not need a new DOF communication channel. It should reuse the existing hanging-DOF source import path.

The mesh-entity side also appears covered for normal periodic node/edge/face DOFs. `Mesh_Distributor` expands selected periodic nodes, edges, and faces to their periodic counterparts, and `ProtoMesh::create_periodicities()` reconstructs local periodic pointers. That means a worker that has the periodic hanging DOF should also have the source DOF's mesh basis, or the existing source-DOF import will fail visibly when it calls `mMesh->node/edge/face(id)`.

## Required Implementation Shape

The important missing piece remains `mHangingDOFs` registration:

1. Create normal mesh-basis hanging source relations.
2. Apply periodic DOF relations.
3. Rebuild `mHangingDOFs` from every DOF where `Dof::is_hanging()` is true.
4. Continue with the existing serial removal or MPI hanging-DOF communication.

For each periodic pair, use the smaller global DOF ID as representative:

- if the smaller-ID DOF is free, set the larger-ID DOF source to it with weight `1.0`;
- if the smaller-ID DOF is already hanging, copy its non-hanging source list and weights onto the larger-ID DOF;
- do not set a source to a hanging DOF, because `Dof::set_sources()` asserts that source DOFs are non-hanging.

The implementation must also decide how to handle the case where the larger-ID periodic target already has sources. Current `Dof` does not support source replacement; all `set_source(s)` paths assert that sources are not already allocated. A small, explicit replacement helper or an invariant/assertion is needed.

## Caveats

- No end-to-end MPI periodic test was run in this audit.
- The conclusion assumes periodic source bases are node/edge/face entities covered by the existing distributor's periodic counterpart expansion.
- If future periodic constraints target entities outside that expansion path, then mesh-basis communication would need to be revisited.

## Files Touched

- `devlog/dl20260513_periodic_dof_parallel_logic.md`
- `devlog/README.md`

