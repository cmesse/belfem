# Devlog 2026-04-08 — Buffer Phase 4 Audit

**Date:** 2026-04-08
**Topic:** Read-only audit of thin-shell buffer-domain plumbing for the H-phi formulation
**AIs involved:** Codex
**Codex Audit Confidence:** high

## Summary

Reviewed the current Phase 4 buffer plumbing against the six requested audit points. The core magnetic solve path looks internally consistent: `Buffer` dispatches only to scalar-`phi` kernels, the new topology guard fires before sideset consumers, and loading Magnesia into the magnetic kernel does not by itself create a hidden `rho()` dependency.

Two real gaps remain. First, `Topology::block_map()` still drops `DomainType::Buffer`, so the new air/buffer postprocessor selection is not actually reachable. Second, the thermal IWG still has no `Buffer` branch, so coupled thermal runs would still abort on a buffer block.

Go / no-go recommendation: **no-go for a sign-off end-to-end `hphirun cmake-build-debug/input.conf` run**. The solve path itself is close, but the advertised buffer H/B postprocessing path is still inert because `Topology::block_map()` omits `Buffer`. A smoke run is still reasonable if the immediate goal is only "does the magnetic solve start and assemble?", but I would not treat that as milestone sign-off.

## Key Findings

- **Major | confidence: high** `Topology::update_block_map()` still excludes `DomainType::Buffer` in `src/homology/cl_Topology.cpp:514-523`, even though `select_blocks()` already classifies `Buffer` as a phi block in `src/homology/cl_Topology.cpp:409-413`. `MaxwellFactory::create_postprocessors()` passes `mTopology->block_map()` into the air postprocessor at `src/fem/maxwell/cl_MaxwellFactory.cpp:1876-1881`, and `MaxwellPostprocessor::select_blocks_and_materials()` only sees what that map contains at `src/fem/maxwell/cl_MaxwellPostprocessor.cpp:161-180`. Result: the new Buffer-aware postprocessor logic is still unreachable, so the Phase 4 recovery fix is inert on the current path.

- **Moderate | confidence: high** The thermal side is still incomplete. `IWG_MaxwellThermal::link_to_group()` routes `Ferro` and `Air` to `T_phi` at `src/fem/thermal/cl_IWG_MaxwellThermal.cpp:175-179`, but `DomainType::Buffer` still falls into the default `BELFEM_ERROR` path at `src/fem/thermal/cl_IWG_MaxwellThermal.cpp:181-184`. This does not block `hphirun`, but `hphiTrun` or any coupled thermal path with a buffer block would still abort.

- **Minor | confidence: high** The explicit buffer-boundary side paths are still not fully mirrored. `MaxwellFactory::set_sidesets()` admits `AirPeriodic`, `AirSymmetry`, and `AirAntiSymmetry` but not the corresponding `Buffer*` cases in `src/fem/maxwell/cl_MaxwellFactory.cpp:1709-1718`, even though `IWG_Maxwell`, `FieldList`, and `Topology::select_sidesets()` already know about `BufferSymmetry`, `BufferAntiSymmetry`, and `BufferPeriodic`. `CutFactory::unflag_symmetry_sidesets()` also still omits `BufferSymmetry` at `src/homology/cl_CutFactory.cpp:2797-2800`. These do not block the currently supported "internal thin-shell buffer layer" case because external buffer boundaries are intentionally unsupported, but they are still missed `Air`-family switch sites.

## Verification Notes

- **Phi-path / `rho` safety | confidence: high** `IWG_Maxwell` dispatches `DomainType::Buffer` only to `phi`, `phi_tri3`, `phi_tet4`, or `phi_tri6_tet10` in `src/fem/maxwell/cl_IWG_Maxwell.cpp:252-280`. Those kernels in `src/fem/maxwell/matrices/mt_maxwell_phi.cpp:24-140` assemble only `trans(B)*B` or `mu0*trans(B)*B` into `M()`; none of them calls `aCalc->material()` or `rho()`. Only `phi_ferro()` at `src/fem/maxwell/matrices/mt_maxwell_phi.cpp:48-96` touches material data. So removing the old buffer guard in `set_material` is safe for the magnetic phi path.

- **Material-map / Magnesia loading | confidence: high** `Magnesia` is `MaterialType::NonMetal` in `src/physics/materials/cl_Material_Magnesia.cpp:23-35` and does not define `rho`; `Material::rho()` would assert if called in `src/physics/materials/cl_Material.hpp:1429-1456`. I did not find any iteration over the kernel material map that unconditionally calls `rho()`. `Kernel::add_material()` only stores pointers in the map/vector at `src/fem/kernel/cl_FEM_Kernel.cpp:726-734`, and the downstream Maxwell logic queries materials per block type, not by "all loaded materials". So binding Magnesia to the buffer label is safe on the current magnetic path.

- **Topology guard ordering | confidence: high** The new `BELFEM_ERROR` guard in `Topology::detect_sideset_types()` runs before both `select_blocks()` and `select_sidesets()` in `Topology::run()` at `src/homology/cl_Topology.cpp:38-45`. Unsupported buffer-touching-non-buffer default sidesets therefore fail before any consumer classifies or uses them. The supported internal thin-shell buffer case does not hit this guard: `detect_sideset_types()` only rewrites sidesets that start as `DomainType::Default` at `src/homology/cl_Topology.cpp:174-176`, while the internal layer interfaces are handled inside the thin-shell machinery, and `ThinShellFactory::create_ghost_facets()` already skips interfaces adjacent to `Buffer` at `src/mesh/cl_ThinShellFactory.cpp:1777-1783`.

- **Substring match in `create_buffers()` | confidence: medium (~75%)** I did not find any current BELFEM built-in material name, repo-local example input, or thin-shell layer label where `"buffer"` appears as part of a clearly non-buffer material name. For the present tree and the supplied `cmake-build-debug/input.conf`, `find("buffer") != npos` looks operationally safe. The residual risk is future user-defined labels such as `buffered_hastelloy` or `buffer_cap_cu`, which would be misclassified even if they were not intended as true insulators. So the substring rule is acceptable as a stopgap, not as a durable contract.

- **MPI / rank split | confidence: high** The new root-only buffer node-flagging in `create_postprocessors()` mirrors the existing Air path exactly at `src/fem/maxwell/cl_MaxwellFactory.cpp:1767-1780`. The all-rank `init_fields()` addition also mirrors Air and guards with `block_exists()` at `src/fem/maxwell/cl_MaxwellFactory.cpp:2028-2048`. The MPI-specific bug I found is the `block_map()` omission above: once `mTopology->synchronize_maps()` calls `update_block_map()` on every rank at `src/homology/cl_Topology.cpp:121` and `src/fem/maxwell/cl_MaxwellFactory.cpp:694,738`, `Buffer` is dropped consistently on every rank, so the postprocessor omission is parallel-consistent but still wrong.

- **Checked but not escalated | confidence: medium-high** `InterfaceProcessor` still seeds its phi-side bitset from `groups( DomainType::Air )` only in `src/homology/cl_InterfaceProcessor.cpp:241-296`. I did not classify that as a current blocker because the supported buffer case stays inside the thin-shell path, `create_ghost_facets()` skips buffer-adjacent layer interfaces, and unsupported external buffer sidesets are now rejected earlier in `detect_sideset_types()`.

## Changes Made / Proposed

- No source edits.
- Wrote this devlog entry.
- Appended the headline findings to `todo/ai_exchange.md`.
- Added this entry to `devlog/README.md`.

## Open Questions

- If explicit `buffer` topology blocks outside the thin-shell auto-tag path are going to be supported, `MaxwellFactory::collect_material_labels_from_domains()` still treats `DomainType::Buffer` like a default/inactive case rather than a real material-bearing phi block. I did not classify that as a blocker for the current thin-shell milestone because the current input uses thin-shell layer promotion, not explicit buffer domains.

- If explicit buffer symmetry or periodic boundaries are expected soon, the remaining `Buffer*` boundary cases should be completed through the switch tables before that feature is advertised.

## Files Updated

- /home/christian/codes/belfem/devlog/dl20260408_buffer_phase4_audit.md
- /home/christian/codes/belfem/devlog/README.md
- /home/christian/codes/belfem/todo/ai_exchange.md
