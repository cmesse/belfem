# Devlog 2026-04-24 — ThinShell MasterNode Originals Audit

**Date:** 2026-04-24
**Topic:** Whether `ThinShellFactory::mMasterNodes` can be assumed to contain originals, and whether `->original()` can be removed from `cl_ThinShellFactory.cpp`
**AIs involved:** Codex
**Codex Audit Confidence:** high

## Summary

Read-only audit of the `CutFactory -> MaxwellFactory -> ThinShellFactory` path.

Verdict:

- In the normal Maxwell thin-shell workflow on rank 0, `ThinShellFactory::mMasterNodes` is populated from the pre-duplication thin-shell sideset nodes and is therefore original with respect to the thin-shell opening step.
- This is **not** a justification for removing every `->original()` in `cl_ThinShellFactory.cpp`.
- Several `original()` calls are still semantically required because they operate on facet/edge/curve nodes or layer duplicates, not directly on `mMasterNodes`.

## Key Findings

- `MaxwellFactory::create_cuts_sub_master()` stores `CutFactory::thin_shell_master_nodes()` before `CutFactory::run()`, and `MaxwellFactory::create_thinshells()` later passes that saved list into `ThinShellFactory`. This means the `ThinShellFactory` list is captured before cohomology-cut interface duplication runs.
  - `src/fem/maxwell/cl_MaxwellFactory.cpp:752-767`
  - `src/fem/maxwell/cl_MaxwellFactory.cpp:849-858`

- `CutFactory::duplicate_nodes_on_face_sidesets()` fills `mThinShellMasterNodes` from `tNodes( tIndex )`, i.e. the nodes already present on the thin-shell sidesets before the new thin-shell duplicate nodes are appended to the mesh. The corresponding slave list may contain the new duplicates.
  - `src/homology/cl_CutFactory.cpp:1738-1828`

- The thin-shell duplicate creation path does **not** call `set_original()` on those new thin-shell duplicate nodes. By contrast, the later generic cut-duplicate path does call `set_original()`.
  - Thin-shell opening path: `src/homology/cl_CutFactory.cpp:1761-1805`
  - Generic cut-duplicate path: `src/homology/cl_CutFactory.cpp:2701-2714`

- After cohomology processing, `CutFactory::restore_thin_shell_sidesets()` restores the original thin-shell facets, and `Mesh::finalize()` then calls `update_facet_nodes()`, which repopulates each facet element from its master element. With the current `fix_facet_masters()` convention, ThinShellFactory sees facet surface elements on the original/master side.
  - `src/homology/cl_CutFactory.cpp:848-880`
  - `src/mesh/cl_Mesh.cpp:683-700`

- Removing all `->original()` from `cl_ThinShellFactory.cpp` is unsafe. Many calls are used to canonicalize facet/edge/curve lookups into the shell-node index space or to compare duplicate-aware identities.
  - Node collection: `src/mesh/cl_ThinShellFactory.cpp:483-522`
  - Higher-order normal accumulation: `src/mesh/cl_ThinShellFactory.cpp:655-678`, `src/mesh/cl_ThinShellFactory.cpp:904-939`
  - Layer element assembly: `src/mesh/cl_ThinShellFactory.cpp:1183-1350`
  - Temporary edge canonicalization: `src/mesh/cl_ThinShellFactory.cpp:1377-1476`
  - Layer edge/face cloning: `src/mesh/cl_ThinShellFactory.cpp:1503-1560`
  - Side-connector facet preprocessing: `src/mesh/cl_ThinShellFactory.cpp:2664-2724`

- A narrower simplification may still be possible for call sites where the operand is known to be taken directly from `mMasterNodes` / `aNodes` selected from it, but that requires a separate per-call-site audit and should be guarded by explicit assertions rather than inferred globally.
  - Candidate region: `src/mesh/cl_ThinShellFactory.cpp:1016-1096`

## Changes Made / Proposed

- No source-code changes made.
- Added this read-only devlog and updated `devlog/README.md`.

## Open Questions

- There is no explicit runtime assertion proving every entry of `ThinShellFactory::mMasterNodes` satisfies `tNode == tNode->original()`. The current guarantee is architectural/control-flow based, not encoded as a local invariant.
- If future workflows load or construct meshes that already carry duplicate-linked sideset nodes before `create_thin_shell_cuts()`, the current assumption should be re-checked.

## Files Updated

- devlog/dl20260424_thinshell_masternode_originals_audit.md
- devlog/README.md
