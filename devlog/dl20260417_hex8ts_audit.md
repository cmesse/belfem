# Devlog: HEX8TS Side-Connector Audit

**Date:** 2026-04-17
**Author:** Junie (Codex role)

## Context
Investigating why the `HEX8TS` side-connector might fail to carry current or produce residual errors in thin-shell tape models. Claude proposed 5 hypotheses (H1-H5).

## Findings

### 1. Thermal Solver Omission (CRITICAL)
The most significant finding is that `LeftCoating` and `RightCoating` domain types are completely missing from `ThermalFactory::create_thermal_kernel`. 
- **Code:** `src/fem/thermal/cl_ThermalFactory.cpp` lines 92-100 (block selection) and 178-183 (block activation).
- **Impact:** Side-connector elements are not present in the thermal problem. They have no thermal DOFs. Maxwell thermal kernels (`h_tb_t`) sampling these nodes will see uninitialized (zero) temperature data. Joule heating in the wrap is not accounted for in the thermal simulation.

### 2. Thermal Kernel Node Sampling
The `h_tb_t` kernel (Maxwell-Thermal coupling for the thin beam) has a logic flaw in its node sampling.
- **Code:** `src/fem/maxwell/matrices/mt_maxwell_h.cpp` lines 1865-1875.
- **Logic:** It samples nodes `{1,2,5,6}` or `{0,3,4,7}`. In the `HEX8TS` topology, these are the **segment end-faces** (cross-sections of the wrap along the curve), not the **interface faces** (tape-side vs air-side).
- **Impact:** The temperature used for resistivity calculation is an average of two inner and two outer nodes at the segment ends, which is physically less representative than averaging the inner face nodes `{0,1,2,3}`.

### 3. Verification of Core Architecture
Hypotheses **H1** (Wiring) and **H2** (Cut Topology) were investigated and found to be structurally sound:
- **Edge Orientation:** Nedelec assembly in `HEX8TS` is invariant to element handedness because it uses `absDetJ` for the integration volume and the Piola transform handles the sign consistency of the basis functions relative to the mesh edges.
- **Cut Persistence:** `ThinShellFactory` correctly iterates over `NodeDuplicates` during layer extrusion, ensuring that if a node on the initial sideset was on a cohomology cut, all extruded layer nodes are also on the cut.

## Conclusion
The zero-current symptom is likely an interpretation issue (H5, symmetric bypass) or a post-processing artifact (missing block IDs in postprocessor). However, the thermal coupling is definitely broken for the side connectors due to the `ThermalFactory` omission.
