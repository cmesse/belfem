# Devlog 2026-04-06 — Thin-Shell Patch-Test Topology Audit

**Date:** 2026-04-06
**Topic:** Read-only audit of thin-shell DOF creation and interface topology after metal-only multilayer runs still showed large residuals
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** medium-high (~80%)
**Literature References:** N/A

## Summary

Investigated whether the large residual in metal-only thin-shell stacks is caused by a structural topology / DOF issue rather than by material nonlinearity or solver tuning. The audit traced:

- thin-shell mesh extrusion and ghost-facet creation,
- sideset and block element creation,
- hanging-edge source construction on shell-air and shell-shell interfaces,
- DOF-level T-matrix expansion for hanging edges,
- shell-specific block/facet ownership logic.

I did not yet find a smoking-gun sign error in the edge-to-node or edge-to-edge hanging transfer itself. The strongest remaining conceptual suspect is a possible mismatch between the shell extrusion normal and the facet master/slave convention used when attaching shell bottom/top faces to adjacent volume domains.

## Key Findings

- `MaxwellFactory::hang_thinshell_edges_on_nodes_bottom/top()` and `hang_thinshell_edges_on_edges_bottom/top()` in `src/fem/maxwell/cl_MaxwellFactory.cpp` are internally consistent on local ordering:
  - bottom uses master-face ordering directly,
  - top reorders slave nodes/edges into master orientation,
  - node-based hanging edges map shell edge endpoints to volume nodes by local facet position,
  - edge-based hanging edges compute an orientation sign from matched local node indices.

- `DofData::create_dofwise_t_matrices_master()` in `src/fem/kernel/cl_FEM_DofMgr_DofData.cpp` expands linear edge-to-node hanging constraints with coefficients `[+1, -1]`, which is consistent with `h = -grad(phi)` on an oriented air edge.

- `BlockData::collect_thin_shell_facet_ids()` in `src/fem/kernel/cl_FEM_DofMgr_BlockData.cpp` has a real bug:
  - it counts only selected thin-shell facets,
  - then populates the output vector with facets from all thin shells.
  This can corrupt the facet list when multiple thin shells exist, but it is probably not the root cause for a single-stack Cu/Ag patch-test failure.

- The dedicated thin-shell FEM element constructor in `src/fem/kernel/cl_FEM_Element.cpp` appears unused in the active Maxwell path.
  - `SideSetData::create_sidesets()` builds ordinary `SideSet` objects in `src/fem/kernel/cl_FEM_DofMgr_SideSetData.cpp`.
  - `SideSet::initialize_elements()` then always uses the generic sideset-element constructor in `src/fem/kernel/cl_FEM_SideSet.cpp`.
  - I could not find an active call site for the specialized `Element(..., aLayers, ...)` constructor.
  This means the 2D-only `compute_edge_directions_thinshell()` routine is likely dead or at least non-critical for the present Maxwell run.

- The strongest conceptual suspect after this audit is the shell-orientation convention:
  - `ThinShellFactory` extrudes layers using normals derived from the geometric facet node ordering in `src/mesh/cl_ThinShellFactory.cpp`.
  - `MaxwellFactory::create_hanging_edges_and_facets()` assumes shell face `0` corresponds to the original facet master side and face `1` to the slave side.
  - That assumption is only safe if the facet node ordering implies a normal from master to slave.
  - I have not yet proved that this is guaranteed for the active 3D mesh path.

- `IWG::number_of_dofs_per_element( SideSet* )` in `src/fem/iwg/cl_IWG.cpp` undercounts `ThinShell` and `Ghost` sidesets by ignoring edge/face DOFs.
  - In first-order triangular cases this can match by coincidence because `#nodes == #edges == 3`.
  - Worth fixing, but not yet tied directly to the observed patch-test failure.

## Changes Made / Proposed

- No source-code changes were made in this audit.
- Added a Codex audit entry to `todo/ai_exchange.md` summarizing the current evidence and ranking of suspects.

## Open Questions

- Is the thin-shell extrusion normal guaranteed to point from facet master to facet slave for the actual 3D element types in use?
- If that convention is correct, does the remaining patch-test failure come from the internal ghost transmission algebra in `h_ghost()` rather than from the DOF construction?
- Is the `collect_thin_shell_facet_ids()` bug active in the current case, or only latent because there is just one thin shell?

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260406_thinshell_patchtest_topology_audit.md
