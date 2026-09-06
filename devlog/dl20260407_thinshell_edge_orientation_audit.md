# Devlog 2026-04-07 — Thin-Shell Edge Orientation Audit

**Date:** 2026-04-07
**Topic:** Read-only audit of thin-shell edge orientation versus facet orientation, with emphasis on the ID-based edge rule
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Checked whether the ID-based edge canonicalization can silently flip thin-shell edge orientations relative to the corresponding facet, with focus on the active first-order `TRI3 -> PENTA6TS` path used in the multilayer patch-test investigation.

Conclusion: I do not see an orientation mismatch bug in that path. The shell element preserves the original facet's local node/edge ordering, while the edge objects are canonicalized by IDs and then corrected back to local orientation through explicit edge-direction signs.

## Key Findings

- [cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp) inserts `PENTA6TS` face nodes in the exact original facet order:
  - bottom face `0,1,2`
  - top face `3,4,5` corresponding to the same facet nodes `0,1,2`
- [cl_Element_PENTA6TS.hpp](/home/christian/codes/belfem/src/mesh/cl_Element_PENTA6TS.hpp) preserves the same local triangular edge cycle on both shell faces:
  - bottom: `(0,1), (1,2), (2,0)`
  - top: `(3,4), (4,5), (5,3)`
- [cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp) canonicalizes temporary facet edges by node IDs in `create_temporary_edges()`, and `create_edges_on_layers()` copies that orientation onto shell layers verbatim.
- [cl_Mesh.cpp](/home/christian/codes/belfem/src/mesh/cl_Mesh.cpp) and [cl_FEM_Element.cpp](/home/christian/codes/belfem/src/fem/kernel/cl_FEM_Element.cpp) then recover the local sign by comparing each element's local edge node order with the stored edge object's order.
- [cl_EF_PENTA6TS.cpp](/home/christian/codes/belfem/src/fem/interpolation/nedelec/cl_EF_PENTA6TS.cpp) uses `aElement->edge_directions( mS )`, so the shell interpolation acts on the element-local orientation, not only on the canonical edge orientation.
- [cl_MaxwellFactory.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_MaxwellFactory.cpp) accounts for reversed edge orientation in shell-to-volume hanging-edge coupling by computing `tWeight = +/- 1`, and the slave-side top interface is explicitly reordered into master orientation before the sign test.
- [cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp) creates internal shell-shell ghost facets with master face `1`, slave face `0`, and orientation `1`, which is consistent because both faces were built from the same original facet order.

## Changes Made / Proposed

- No source-code changes made.
- Added a Codex audit note to `todo/ai_exchange.md`.

## Open Questions

- The current conclusion applies to the active first-order `TRI3 -> PENTA6TS` path. Higher-order `TRI6 -> PENTA18TS` deserves a separate audit.
- The patch-test failure is more likely to live in the layered interface coupling path than in thin-shell edge orientation itself.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260407_thinshell_edge_orientation_audit.md
- devlog/README.md
