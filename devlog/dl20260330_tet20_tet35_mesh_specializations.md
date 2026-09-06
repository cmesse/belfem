# Devlog 2026-03-30 — TET20 And TET35 Mesh Specializations

**Date:** 2026-03-30
**Topic:** Add missing mesh element specializations for `TET20` and `TET35`
**AIs involved:** Codex
**Codex Audit Confidence:** high

## Summary

Verified that the cubic and quartic tetrahedron parametric node order already matched the interpolation-side definitions and then added the missing mesh specializations for `TET20` and `TET35`.

The new headers provide `type()`, `dimension()`, `get_nodes_of_facet()`, `get_corner_nodes_of_facet()`, `get_nodes_of_edge()`, and `get_edges_of_facet()`, using the existing `TET10` facet orientation convention while preserving `TRI10` / `TRI15` surface-node order on extracted facets.

## Key Findings

- The parametric node ordering in [src/mesh/cl_Element_Factory.cpp](/home/christian/codes/belfem/src/mesh/cl_Element_Factory.cpp#L509) and [src/mesh/cl_Element_Factory.cpp](/home/christian/codes/belfem/src/mesh/cl_Element_Factory.cpp#L525) already matched [src/fem/interpolation/lagrange/cl_IF_TET20.hpp](/home/christian/codes/belfem/src/fem/interpolation/lagrange/cl_IF_TET20.hpp) and [src/fem/interpolation/lagrange/cl_IF_TET35.hpp](/home/christian/codes/belfem/src/fem/interpolation/lagrange/cl_IF_TET35.hpp).
- Gmsh reference faces in `more/gmsh/tet20.msh` and `more/gmsh/tet35.msh` were used to derive the high-order face-node order, then rotated as needed to match BELFEM's existing tetra facet orientation convention from [src/mesh/cl_Element_TET10.hpp](/home/christian/codes/belfem/src/mesh/cl_Element_TET10.hpp).
- `HEX64` already had a dedicated mesh specialization; only `TET20` and `TET35` were missing.

## Changes Made / Proposed

- Added [src/mesh/cl_Element_TET20.hpp](/home/christian/codes/belfem/src/mesh/cl_Element_TET20.hpp)
- Added [src/mesh/cl_Element_TET35.hpp](/home/christian/codes/belfem/src/mesh/cl_Element_TET35.hpp)
- Updated [src/mesh/cl_Element_Factory.cpp](/home/christian/codes/belfem/src/mesh/cl_Element_Factory.cpp#L29) to include the new headers
- Updated [todo/tests/tests_10_mesh.md](/home/christian/codes/belfem/todo/tests/tests_10_mesh.md#L96) to remove the now-stale warning about missing tetra mesh headers

## Open Questions

- None for the mesh specializations themselves. The next reasonable step is to add or update automated tests that exercise `TET20` and `TET35` in the mesh catalog and facet-extraction paths.

## Files Updated

- src/mesh/cl_Element_TET20.hpp
- src/mesh/cl_Element_TET35.hpp
- src/mesh/cl_Element_Factory.cpp
- todo/tests/tests_10_mesh.md
- todo/ai_exchange.md
- devlog/dl20260330_tet20_tet35_mesh_specializations.md
