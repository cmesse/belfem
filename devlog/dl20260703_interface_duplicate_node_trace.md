# Devlog 2026-07-03 — Interface Duplicate Node Trace

**Date:** 2026-07-03
**Topic:** Locate interface node duplication for ParaView output and coil extension point
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Read-only trace of where BELFEM duplicates nodes at Maxwell material interfaces for postprocessing/visualization. The active duplication path is the homology `InterfaceProcessor`, invoked from `CutFactory::compute_thin_cuts_and_duplicate_interface_nodes()`.

## Key Findings

- `CutFactory::run()` calls `compute_thin_cuts_and_duplicate_interface_nodes()` after cohomology and before mesh refinalization (`src/homology/cl_CutFactory.cpp:155`).
- That helper invokes `InterfaceProcessor` after thin-cut processing and duplicate/original linking (`src/homology/cl_CutFactory.cpp:680-700`).
- `InterfaceProcessor` currently seeds only Air and Ferro block bitsets (`src/homology/cl_InterfaceProcessor.cpp:281-293`) and admits interfaces if at least one side is Air or Ferro while excluding Air-Air (`src/homology/cl_InterfaceProcessor.cpp:336-341`).
- Actual duplicate allocation happens in `InterfaceSet::duplicate_nodes()` (`src/homology/cl_InterfaceProcessor.cpp:147-170`), relinking happens in `InterfaceSet::relink_elements()` (`src/homology/cl_InterfaceProcessor.cpp:220-235`), and mesh insertion/original duplicate unification happens in `InterfaceProcessor::add_duplicate_nodes_to_mesh()` / `unify_duplicates()` (`src/homology/cl_InterfaceProcessor.cpp:441-595`).
- Topology already classifies coil interfaces as `InterfaceAirCoil` and `InterfaceFerroCoil` (`src/homology/cl_Topology.cpp:309-334`) and includes them in `phi_interface_ids()` (`src/homology/cl_Topology.cpp:501-507`), but the current interface duplication filter does not explicitly track `DomainType::Coil`.

## Changes Made / Proposed

- No source changes made.
- Proposed extension point: add a coil membership bitset in `InterfaceProcessor` and include it in the interface admission predicate if coil-side duplicate visualization is desired.

## Open Questions

- Whether coil interface sidesets should remain geometry/postprocessing-only or be activated in Maxwell assembly. `cl_IWG_Maxwell.cpp` currently errors if `InterfaceAirCoil` or `InterfaceFerroCoil` reaches assembly.

## Files Updated

- `devlog/dl20260703_interface_duplicate_node_trace.md`
