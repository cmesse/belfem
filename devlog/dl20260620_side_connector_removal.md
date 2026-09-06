# Devlog 2026-06-20 — Side-Connector Code Removal

**Date:** 2026-06-20
**Topic:** Removed the dead "side connector" thin-shell feature ( deemed unphysical on 2026-06-05 — free binormal-H at the fold; the HEX8TS wrap was deleted then, but the code removal was left pending ). Closes that thread.
**Module:** src/fem/maxwell, src/mesh, src/fem/kernel, src/homology
**AIs involved:** Claude ( scoping + the `h_tb`/`h_tb_t` removal ), a focused subagent ( the mechanical multi-file removal ), Christian ( scope: what stays vs goes ).
**Claude Confidence:** high — 0 remaining references to the removed symbols; all touched translation units compile under `-DDEBUG -Werror`.

---

## What was removed

The side-connector feature was still threaded through ~11 files via three `DomainType` values and the `h_penalty` IWG.

- **`DomainType::LeftCoating` (6), `RightCoating` (7), `InterfaceTsConnector` (31)** — removed from `en_DomainType.hpp` ( gaps left at 6/7/31, no renumbering ) and `en_DomainType.cpp` ( `to_string` ), plus every switch-case across `cl_IWG_Maxwell.cpp`, `cl_IWG.cpp`, `cl_MaxwellFactory.cpp`, `cl_Topology.cpp`, `cl_MaxwellPostprocessor.cpp`, `cl_Maxwell_FieldList.cpp`. Stale connector comments cleaned in `cl_FEM_Element.cpp` and `cl_Mesh_OrderConverter.cpp`.
- **`h_penalty`** ( the `InterfaceTsConnector` IWG ) — removed from `mt_maxwell_h.{hpp,cpp}`.
- **`h_tb` / `h_tb_t`** ( temperature-bulk h-matrices ) — these were dispatched *only* by the now-removed `LeftCoating`/`RightCoating` cases, so they became dead; removed ( decl + def ). `h_ghost` ( the working ghost-facet stabilization ) was left fully intact right after them.

## What was deliberately kept

- **`ThinShellFactory::compute_binomial_vectors`** ( + its "side-connector sign" `BELFEM_ERROR` ) — needed later ( Christian ); self-contained, no dependency on the removed code.
- **Ghost-facet stabilization** — `mCreateGhostFacets`, `create_ghost_facets`, `GhostFacets`, `DomainType::Ghost`, `h_ghost`. This is the *working* stabilization; the side connectors were a separate dead end ( not the ghost facets, which were initially confused — corrected mid-task ).
- **`cl_EF_HEX8TS.cpp`** ( harmless connector special-case ) and **`cl_EdgeCutter` / `corctest.cpp`** ( general edge utility, not the FEM connector feature ).

## Verification

`grep` for `LeftCoating|RightCoating|InterfaceTsConnector|h_penalty|h_tb|h_tb_t` across `src/` ( excluding `doc/` ): **0 matches**. Syntax-checks ( `mpicxx -fsyntax-only`, real `flags.make`, `-DDEBUG`/`-Werror` ): clean for `en_DomainType.cpp`, `cl_IWG_Maxwell.cpp`, `mt_maxwell_h.cpp`, `cl_MaxwellFactory.cpp`, `cl_Topology.cpp`, `cl_Maxwell_FieldList.cpp`, `cl_MaxwellPostprocessor.cpp`, `cl_IWG.cpp`. Full `make` confirmation pending ( Christian runs builds ).

## Files Updated

- src/mesh/en_DomainType.hpp, en_DomainType.cpp
- src/fem/maxwell/matrices/mt_maxwell_h.hpp, mt_maxwell_h.cpp
- src/fem/maxwell/cl_IWG_Maxwell.cpp, cl_MaxwellFactory.cpp, cl_MaxwellPostprocessor.cpp, cl_Maxwell_FieldList.cpp
- src/fem/iwg/cl_IWG.cpp
- src/homology/cl_Topology.cpp
- src/fem/kernel/cl_FEM_Element.cpp, src/mesh/cl_Mesh_OrderConverter.cpp ( comments )
