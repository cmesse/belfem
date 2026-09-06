# Devlog 2026-04-09 — Buffer Cut Plan v2

**Date:** 2026-04-09
**Topic:** Read-only trace and rewrite of the buffer-cut implementation plan
**AIs involved:** Claude, Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** Messe et al. 2023 (paper1) §1.3, Alves et al. 2022a (paper5) §II-B / Eq. (10), Alves et al. 2024 (paper8) §II, Schnaubelt et al. 2023 (paper9) §IV, Pellikka et al. 2013

## Summary

Traced the real thin-shell, cut, kernel, and thermal data flow behind the buffer-cut proposal and rewrote the implementation plan around the actual prerequisites. The physics direction remains plausible, but the viable hook points differ from the original plan: node duplicates must be introduced at the buffer interface inside `ThinShellFactory`, the cut reroute must happen after `create_thinshells()` but before the Maxwell kernel is constructed, and the thermal continuity hook must be registered before the thermal kernel / DOF manager freeze the hanging-basis state.

## Key Findings

- Thin-shell layer interfaces still share the same `Node*` across adjacent blocks. `ThinShellFactory::Layer` has no `NodeDuplicates`; the current `hasDuplicates` path only duplicates edges/faces. [confidence: high; `src/mesh/cl_ThinShellFactory.hpp:35-44`, `src/mesh/cl_ThinShellFactory.cpp:1167-1180`, `src/mesh/cl_ThinShellFactory.cpp:1258-1275`, `src/mesh/cl_ThinShellFactory.cpp:1473-1661`]
- Persisting full `Layer*` data in `ThinShell` is not viable. `ThinShellFactory::create()` moves the temp containers into the mesh with `append_move()`, which clears the source cells, and then deletes the temp `Layer*`. [confidence: high; `src/mesh/cl_ThinShellFactory.cpp:280-299`, `src/mesh/cl_ThinShellFactory.cpp:342-346`, `src/containers/cl_Cell.hpp:502-515`]
- The cut topology is created before thin-shell layer blocks exist. `create_cuts()` precedes `create_thinshells()`, so the current `CutFactory` timing cannot see the thin-shell buffer block even though `Topology::select_blocks()` now classifies `DomainType::Buffer` as phi. [confidence: high; `src/fem/maxwell/cl_MaxwellFactory.cpp:395-398`, `src/fem/maxwell/cl_MaxwellFactory.cpp:681-724`, `src/homology/cl_Topology.cpp:398-420`]
- The real mechanical meaning of the cut jump is a node source relation that includes abstract cut nodes. The reroute must therefore rewrite the relevant node source/original/duplicate relation before the kernel collects the hanging basis. [confidence: high; `src/homology/cl_CutSet.cpp:42-69`, `src/homology/cl_CutFactory.cpp:2619-2749`, `src/fem/kernel/cl_FEM_Kernel.cpp:103-107`]
- Thermal block selection still omits `DomainType::Buffer`, and thermal continuity cannot be added late because both `Kernel::Kernel()` and `DofManager::initialize()` establish the hanging-basis / hanging-DOF state before solve time. [confidence: high; `src/fem/thermal/cl_ThermalFactory.cpp:63-104`, `src/fem/thermal/cl_ThermalFactory.cpp:133-157`, `src/fem/kernel/cl_FEM_Kernel.cpp:103-107`, `src/fem/kernel/cl_FEM_DofManager.cpp:137-147`, `src/fem/kernel/cl_FEM_DofMgr_DofData.cpp:2327-2429`]
- Saved thin-shell node indices are not a stable post-factory matching key because the mesh later renumbers nodes and may refresh facet nodes from masters. [confidence: high; `src/mesh/cl_ThinShellFactory.cpp:480-520`, `src/mesh/cl_Mesh.cpp:428-434`, `src/mesh/cl_Mesh.cpp:683-701`, `src/fem/maxwell/cl_MaxwellFactory.cpp:421-425`]

## Changes Made / Proposed

- Wrote the revised plan in `todo/buffer_cut_implementation_plan_v2.md`.
- Appended the headline correction to `todo/ai_exchange.md`.
- Added this devlog entry.
- Updated `devlog/README.md`.
- No source-code changes.

## Open Questions

- Should the first node-duplicate implementation split all buffer-adjacent interface nodes or only the cut-intersection subset? The plan is written so either can work, provided the thermal pair export matches the final duplication scope. [confidence: medium (~70%)]
- What is the cleanest pre-kernel thermal constraint API: mesh sidecar consumed by `Kernel`/`DofData`, or a thermal-factory registration path that populates node pair data before `Kernel::Kernel()` runs? [confidence: medium (~75%)]
- The sign/orientation of the rerouted cut still needs one debugger trace on a real buffer case before coding. [confidence: medium (~75%)]

## Files Updated

- todo/buffer_cut_implementation_plan_v2.md
- devlog/dl20260409_buffer_cut_plan_v2.md
- devlog/README.md
- todo/ai_exchange.md

