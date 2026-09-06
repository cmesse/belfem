# Devlog 2026-04-09 — Buffer Cut Plan Review

**Date:** 2026-04-09
**Topic:** Read-only review of `todo/buffer_cut_implementation_plan.md`
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** Messe et al. 2023 (paper1) §1.3, Alves et al. 2022a (paper5) §II-B / Eq. (2), Alves et al. 2024 (paper8) §II, Schnaubelt et al. 2023 (paper9) §IV, Pellikka et al. 2013

## Summary

Reviewed the proposed buffer-cut rerouting plan against the current Maxwell, thin-shell, homology, and thermal implementations. The single-cut insulated-limit idea is physically plausible, but the implementation plan currently targets the wrong ownership points and assumes nodal structures that do not exist in the current thin-shell mesh.

## Key Findings

- Step 1 currently would not preserve usable layer data: `ThinShellFactory::create()` moves `Layer::Nodes/Edges/Faces` into the mesh with `append_move()`, which clears the source containers before the temporary `Layer*` objects are deleted. Persisting those `Layer*` objects in `ThinShell` without redesign would leave empty containers. (`src/mesh/cl_ThinShellFactory.cpp:283-299`, `src/containers/cl_Cell.hpp:502-515`)
- Step 3 misidentifies `set_original()` as a thermal-healing mechanism. It is duplicate bookkeeping, and the current thin-shell implementation does not create interface-side `NodeDuplicates` at all. Adjacent layer blocks still share the same `Node*` objects. (`src/mesh/cl_ThinShellFactory.cpp:1070-1088`, `src/mesh/cl_ThinShellFactory.cpp:1243-1275`, `src/fem/maxwell/doc/thinshell_postprocessor_node_sharing.md:57-107`)
- Step 4 assumes a distinct buffer-node / HTS-partner-node pair can be rewired inside `create_hanging_edges_and_facets()`, but cut support is created earlier by `CutFactory::run()`, while `create_hanging_edges_and_facets()` only handles interface edge/face hanging relations. (`src/fem/maxwell/cl_MaxwellFactory.cpp:681-780`, `src/fem/maxwell/cl_MaxwellFactory.cpp:1102-1210`)
- Step 5 is too late in the thermal lifecycle and also collides with current block selection: thermal hanging constraints are established during DOF-manager initialization, and `ThermalFactory` still omits `DomainType::Buffer` from the selected thermal block list. (`src/fem/kernel/cl_FEM_DofMgr_DofData.cpp:2327-2465`, `src/fem/thermal/cl_ThermalFactory.cpp:63-104`, `src/fem/thermal/cl_IWG_MaxwellThermal.cpp:175-180`)

## Changes Made / Proposed

- No source-code changes.
- Added this devlog entry.
- Appended a Codex audit note to `todo/ai_exchange.md`.

## Open Questions

- Should the implementation pivot to an explicit `NodeDuplicates` extension for thin-shell layer interfaces before attempting any buffer-cut rerouting?
- Should cut rerouting live in `CutFactory` / thin-shell cut restoration rather than `MaxwellFactory::create_hanging_edges_and_facets()`?
- Is the intended thermal treatment to include buffer blocks in the thermal kernel, or to keep thermal continuity through shared nodes and constrain only electromagnetic duplicates?

## Files Updated

- devlog/dl20260409_buffer_cut_plan_review.md
- devlog/README.md
- todo/ai_exchange.md
