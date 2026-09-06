# Buffer Cut Rerouting: Implementation Plan v2

**Date:** 2026-04-09  
**Author:** Codex, on Claude's rewrite request  
**Purpose:** Replace `todo/buffer_cut_implementation_plan.md` with a version that matches the current BELFEM thin-shell, cut, kernel, and thermal pipelines.  
**Design note:** `src/fem/maxwell/doc/buffer_cut_topology.md`  
**Status:** Read-only rewrite. No source code changes yet.

---

## Scope

The physics rationale is unchanged: keep the single-cut contract and reroute the cut jump through the buffer-side node so that the imposed transport current encloses the HTS loop only. That remains plausible in the insulated-limit topology described in `src/fem/maxwell/doc/buffer_cut_topology.md` and is consistent with the local literature on thin cuts / cohomology cuts and multiply connected current paths. [confidence: high; `literature/papers/fem/alves2024.txt:17-18,174-222`, `literature/papers/fem/alves2022a.txt:279-296`, `literature/papers/fem/schnaubelt2023.txt:343-351`]

This rewrite only corrects the implementation plan. It does not re-derive the physics and it does not reopen the separate contact-impedance milestone.

---

## Fixed Premises

1. **At most one buffer layer and one HTS layer per tape.** Enforce with `BELFEM_ERROR`.  
2. **Linear elements only** for the first implementation. The current shell hanging pass already enforces this at `src/fem/maxwell/cl_MaxwellFactory.cpp:1114-1115`.  
3. **Thermal continuity must remain exact** across the stack, including any buffer-interface node pairs introduced for the magnetic reroute.  
4. **No-buffer tapes must remain bit-identical.** The new node-duplicate logic must not change ordinary conductor-conductor interfaces.  
5. **One cut, one current, no terminal changes.** The reroute is internal to the cut/shell interaction.

---

## Trace Summary

### A.1 Node sharing at layer interfaces

- `ThinShellFactory::Layer` currently has `hasDuplicates`, `Edges`, `EdgeDuplicates`, `Faces`, `FaceDuplicates`, and `GhostFacets`, but **no `NodeDuplicates` container**. [confidence: high; `src/mesh/cl_ThinShellFactory.hpp:35-44`]
- Shared `Node*` ownership across adjacent thin-shell blocks is created when element blocks are assembled directly from `Layer::Nodes`. The same layer array is used as the top row of block `b-1` and the bottom row of block `b`. [confidence: high; `src/mesh/cl_ThinShellFactory.cpp:1167-1180`, `src/mesh/cl_ThinShellFactory.cpp:1258-1275`]
- The existing `hasDuplicates` mechanism is **edge/face-only**. It is set from material changes in `ThinShellFactory::create()` and later drives `EdgeDuplicates` / `FaceDuplicates` plus ghost facets. [confidence: high; `src/mesh/cl_ThinShellFactory.cpp:187-195`, `src/mesh/cl_ThinShellFactory.cpp:1473-1513`, `src/mesh/cl_ThinShellFactory.cpp:1519-1557`, `src/mesh/cl_ThinShellFactory.cpp:1560-1628`, `src/mesh/cl_ThinShellFactory.cpp:1638-1661`, `src/mesh/cl_ThinShellFactory.cpp:1775-1783`]
- `create_buffers()` does not change `hasDuplicates` and does not add any node-level split at the buffer interface. It only tags the block as `DomainType::Buffer` and makes buffer edges/faces source themselves. [confidence: high; `src/mesh/cl_ThinShellFactory.cpp:1665-1751`]
- Therefore the buffer-conductor node pair required by the rerouted cut does not exist today. The safe extension is **not** to overload `hasDuplicates`; it is to add a separate node-duplicate path that fires only on buffer-adjacent interfaces, leaving ordinary conductor-conductor node sharing untouched. [confidence: high; code trace above, plus `src/fem/maxwell/doc/thinshell_postprocessor_node_sharing.md:154-227`]

### A.2 Layer data lifetime

- Temporary `Layer*` objects are allocated in `ThinShellFactory::create()`, their containers are moved into the mesh with `append_move()`, and then the `Layer*` objects are deleted. `append_move()` clears the source container after the move. [confidence: high; `src/mesh/cl_ThinShellFactory.cpp:180-185`, `src/mesh/cl_ThinShellFactory.cpp:280-299`, `src/containers/cl_Cell.hpp:502-515`, `src/mesh/cl_ThinShellFactory.cpp:342-346`]
- `ThinShell` does **not** retain per-layer node/edge/face containers. It stores block pointers, thicknesses, material labels, and a saved `mNodeIndices` list. [confidence: high; `src/mesh/cl_ThinShell.hpp:42-53`, `src/mesh/cl_ThinShell.cpp:55-60`]
- The saved `mNodeIndices` are not a stable post-factory lookup key. They are built from the master shell node selection in `collect_nodes()`, then later the mesh renumbers all nodes with `Mesh::update_node_indices()`. Facet nodes can also be rewritten from their masters by `Mesh::update_facet_nodes()`. [confidence: high; `src/mesh/cl_ThinShellFactory.cpp:480-520`, `src/mesh/cl_ThinShellFactory.cpp:353`, `src/mesh/cl_Mesh.cpp:428-434`, `src/mesh/cl_Mesh.cpp:683-701`, `src/fem/maxwell/cl_MaxwellFactory.cpp:421-425`]
- Minimal persistent data should therefore be **compact sidecar metadata**, not full `Layer` ownership. The current shell block order plus a final list of buffer duplicate/original node pairs is sufficient for the reroute and the thermal follow-up. [confidence: high]

### A.3 CutFactory pipeline and the real reroute hook

- `MaxwellFactory` calls `create_cuts()` before `create_thinshells()`. [confidence: high; `src/fem/maxwell/cl_MaxwellFactory.cpp:395-398`]
- `create_cuts()` reads thin-shell input only into `Protoshell` metadata and marks the original mesh sidesets as `DomainType::ThinShell`. It does **not** create the layered thin-shell blocks or layer nodes. [confidence: high; `src/fem/maxwell/cl_MaxwellFactory.cpp:681-724`, `src/fem/maxwell/cl_MaxwellFactory.cpp:2476-2626`]
- `create_cuts_sub_master()` calls `CutFactory::create_thin_shell_cuts()` before `CutFactory::run()`. That step duplicates nodes on the original thin-shell sideset, stores `mThinShellMasterNodes` / `mThinShellSlaveNodes`, relinks slave elements, and duplicates the shell facets. [confidence: high; `src/fem/maxwell/cl_MaxwellFactory.cpp:747-770`, `src/homology/cl_CutFactory.cpp:1573-1613`, `src/homology/cl_CutFactory.cpp:1616-1832`, `src/homology/cl_CutFactory.cpp:1834-1907`, `src/homology/cl_CutFactory.cpp:1910-1958`]
- The actual cut topology is then created in `CutFactory::run()`: cohomologies are computed, thin cuts are created, interface nodes are duplicated, and the mesh is finalized again. [confidence: high; `src/homology/cl_CutFactory.cpp:107-164`, `src/homology/cl_CutFactory.cpp:817-838`]
- Abstract cut nodes are created in `CutProcessor::create_abstract_nodes()`. Cut duplicate nodes get source lists `[abstract nodes..., original node]` in `CutSet::create_duplicates()`. Those source lists are later converted back into original/duplicate bookkeeping by `CutFactory::link_node_duplicates_and_originals()`. [confidence: high; `src/homology/cl_CutProcessor.cpp:780-787`, `src/homology/cl_CutSet.cpp:42-69`, `src/homology/cl_CutProcessor.cpp:1167-1183`, `src/homology/cl_CutProcessor.cpp:1188-1237`, `src/homology/cl_CutFactory.cpp:2619-2749`]
- The cut coefficient becomes a global `phi` DOF because `IWG_Maxwell` marks abstract nodes with abstract dof type `phi`, `DofData` creates those DOFs, and `IWG_Maxwell::collect_abstract_node_dofs()` gathers them for the controller / BC path. [confidence: high; `src/fem/maxwell/cl_IWG_Maxwell.cpp:186-191`, `src/fem/kernel/cl_FEM_DofMgr_DofData.cpp:1068-1083`, `src/fem/maxwell/cl_IWG_Maxwell.cpp:723-744`, `src/fem/maxwell/cl_IWG_Maxwell.cpp:88-100`, `src/fem/kernel/cl_FEM_Controller.cpp:839-850`]
- **Critical timing correction:** the thin-shell buffer block does not exist when `CutFactory` runs, even though `Topology::select_blocks()` now classifies `DomainType::Buffer` as phi. [confidence: high; `src/homology/cl_Topology.cpp:398-420`, `src/fem/maxwell/cl_MaxwellFactory.cpp:395-398`, `src/fem/maxwell/cl_MaxwellFactory.cpp:845-855`]
- Therefore the viable hook in the current codebase is **not** an internal `CutFactory::run()` modification. It is a **post-CutFactory, post-`create_thinshells()`, pre-kernel MaxwellFactory pass** that rewrites the relevant node source/original/duplicate relation before the kernel freezes the hanging-basis graph. [confidence: high; `src/fem/kernel/cl_FEM_Kernel.cpp:103-107`]

### A.4 Thermal block selection and timing gap

- `ThermalFactory::create_thermal_kernel()` still omits `DomainType::Buffer` when it selects the thermal blocks. [confidence: high; `src/fem/thermal/cl_ThermalFactory.cpp:63-104`]
- The later thermal activation logic already includes `DomainType::Buffer`, and `IWG_MaxwellThermal` already dispatches Buffer blocks to `T_phi`. [confidence: high; `src/fem/thermal/cl_ThermalFactory.cpp:172-190`, `src/fem/thermal/cl_IWG_MaxwellThermal.cpp:175-180`]
- Thermal continuity cannot be added as a late post-processing patch. `Kernel::Kernel()` collects / expands the mesh hanging basis before the field is created, and `DofManager::initialize()` later turns that mesh-basis state into DOF-wise T-matrices. [confidence: high; `src/fem/kernel/cl_FEM_Kernel.cpp:103-107`, `src/fem/thermal/cl_ThermalFactory.cpp:133-157`, `src/fem/kernel/cl_FEM_DofManager.cpp:137-147`, `src/fem/kernel/cl_FEM_DofMgr_DofData.cpp:2327-2429`, `src/fem/kernel/cl_FEM_DofMgr_DofData.cpp:3562-3594`]
- Because the mesh basis state is shared, the thermal fix should use a **thermal-only pre-kernel constraint path** or a dedicated DofManager-side sidecar. It should not overwrite the magnetic cut sources after the Maxwell kernel has been built. [confidence: high]

---

## Revised Prerequisites

### Prerequisite 1: Node duplicates at buffer-adjacent thin-shell interfaces

This is a mesh-factory change. The current shell stack has no node-level split between buffer and conductor layers.

### Prerequisite 2: Post-CutFactory, pre-kernel reroute hook

The reroute still consumes CutFactory output, but in the current build order it must run from `MaxwellFactory` after `create_thinshells()` and before `Kernel::Kernel()` is called.

### Prerequisite 3: Thermal pre-kernel continuity hook

Thermal equality constraints for the buffer duplicate/original node pairs must be registered before the thermal kernel collects the hanging basis and before the thermal DOF manager initializes.

---

## Ordered Implementation Steps

## Step 1: Add a buffer-only node-duplicate path to `ThinShellFactory`

**What changes**

- Extend `ThinShellFactory::Layer` with a dedicated node path, for example `bool hasNodeDuplicates` plus `Cell<Node*> NodeDuplicates`, instead of reusing `hasDuplicates`. Touch points are `src/mesh/cl_ThinShellFactory.hpp:35-44` and the layer-setup logic in `src/mesh/cl_ThinShellFactory.cpp:187-195`.
- Set `hasNodeDuplicates` only on interfaces adjacent to a buffer layer. The existing `hasDuplicates` logic stays responsible for edge/face duplication on ordinary material changes. `create_buffers()` at `src/mesh/cl_ThinShellFactory.cpp:1665-1751` is the current place where buffer blocks are known by material label; this logic needs to feed the new node flag without changing the old edge/face flag semantics.
- In `create_nodes_on_layers()` (`src/mesh/cl_ThinShellFactory.cpp:991-1092`), allocate and wire the new `NodeDuplicates` array for flagged interfaces, including periodicity and duplicate/original bookkeeping analogous to the existing layer-node bookkeeping.
- In the shell element builders (`src/mesh/cl_ThinShellFactory.cpp:1167-1180`, `src/mesh/cl_ThinShellFactory.cpp:1258-1275` and the analogous higher-order variants, even if the first implementation still errors for non-linear cases), route the interface row through `NodeDuplicates` for the block side that needs to separate from the shared layer node, analogous to how `link_elements_with_edges()` already selects `EdgeDuplicates` for the lower interface of the upper block at `src/mesh/cl_ThinShellFactory.cpp:1574-1575`, `src/mesh/cl_ThinShellFactory.cpp:1602-1604`.

**Why**

This satisfies **Prerequisite 1**. Without a node-level split, there is no buffer-side node available to carry the rerouted cut jump.

**Depends on**

None.

**Risk**

High. This is the first place where a wrong ownership choice can silently change both shell topology and the no-buffer baseline.

**Estimated effort**

1-2 days.

**Notes**

- The no-buffer path stays bit-identical if `hasNodeDuplicates` fires only on buffer-adjacent interfaces.
- The current `hasDuplicates` flag should not be repurposed. That would change conductor-conductor node sharing and create regressions outside the buffer feature.

## Step 2: Add shell discovery and compact buffer metadata in `MaxwellFactory`

**What changes**

- Add a new discovery pass immediately after `create_thinshells()` in `src/fem/maxwell/cl_MaxwellFactory.cpp:397-400`, before kernel construction begins at `src/fem/maxwell/cl_MaxwellFactory.cpp:463`.
- Scan `mMesh->thin_shells()` and each shell's block order (`src/mesh/cl_ThinShell.hpp:42-53`, `src/mesh/cl_ThinShell.hpp:149-153`) to locate:
  - exactly one buffer layer,
  - exactly one HTS layer,
  - whether the HTS lies above or below the buffer in the layer stack.
- Enforce the one-buffer / one-HTS contract with `BELFEM_ERROR`.
- Build a transient `BufferCutShellInfo`-style record local to `MaxwellFactory`. Minimal data should be:
  - `ThinShell*`,
  - buffer block id / index,
  - HTS block id / index,
  - HTS-side orientation relative to the buffer,
  - later, the final buffer duplicate/original node pairs identified by the reroute.
- Do **not** try to persist full `Layer` ownership into `ThinShell`. The temp layer containers are gone by design. [See `src/mesh/cl_ThinShellFactory.cpp:280-299`, `src/mesh/cl_ThinShellFactory.cpp:342-346`, `src/containers/cl_Cell.hpp:502-515`]
- Do **not** rely on saved integer node indices as the long-lived matching key. They are renumbered later. Prefer `Node*` identity and element/facet connectivity. [See `src/mesh/cl_ThinShellFactory.cpp:480-520`, `src/mesh/cl_Mesh.cpp:428-434`, `src/mesh/cl_Mesh.cpp:683-701`]

**Why**

This supplies the minimum shell-specific data needed by **Prerequisite 2** without reviving the broken “store the whole `Layer`” design from v1.

**Depends on**

Step 1.

**Risk**

Medium.

**Estimated effort**

0.5-1 day.

**Notes**

- If HTS detection via attached `Material*` is not yet available at this point in the pipeline, use the thin-shell material labels plus an input-material/property probe. Do not delay the reroute until after kernel creation just to get `Block::material()`, because that would be too late for the hanging-basis graph.

## Step 3: Add the actual cut reroute as a post-thinshell, pre-kernel Maxwell pass

**What changes**

- Insert a new pass in `MaxwellFactory` after `create_thinshells()` and before the Maxwell kernel is constructed. Current anchor points are `src/fem/maxwell/cl_MaxwellFactory.cpp:397-400` and `src/fem/maxwell/cl_MaxwellFactory.cpp:463`.
- The new pass should:
  1. identify cut-crossing shell nodes from the CutFactory-opened shell data (`mThinShellMasterNodes`, `mThinShellSlaveNodes`, plus the duplicate/original links produced by `CutFactory`),
  2. find the corresponding buffer-side node created by Step 1 on the HTS-adjacent interface,
  3. recover the current cut source relation on the old node,
  4. release the old node back to a free unknown,
  5. reattach the cut jump to the buffer node with the same abstract-node source and the correct sign.
- Mechanically, “rerouting” means rewriting the node source/original/duplicate relation that the CutFactory created through `CutSet::create_duplicates()` and `CutFactory::link_node_duplicates_and_originals()`. The cut current itself remains the same abstract `phi` DOF path. [confidence: high; `src/homology/cl_CutSet.cpp:42-69`, `src/homology/cl_CutFactory.cpp:2619-2749`, `src/fem/maxwell/cl_IWG_Maxwell.cpp:186-191`, `src/fem/kernel/cl_FEM_DofMgr_DofData.cpp:1068-1083`]
- Run this pass before the Maxwell kernel constructor collects the hanging basis at `src/fem/kernel/cl_FEM_Kernel.cpp:103-107`.
- Keep the existing single-cut contract: no extra abstract nodes, no new terminals, no second cut coefficient.

**Why**

This is **Prerequisite 2**. It is the actual implementation of the rerouted cut in the only hook point that matches the current pipeline.

**Depends on**

Step 1 and Step 2.

**Risk**

High. This is the core topological change.

**Estimated effort**

1-2 days.

**Notes**

- The exact sign/orientation of the partner node remains a required one-case debug check before coding. `buffer_cut_topology.md` gives the intended rule, but the real sign still depends on the current shell/cut orientation convention. [confidence: medium (~75%)]
- `create_hanging_edges_and_facets()` is still the wrong ownership point for this logic. By that stage the cut graph should already be correct.

## Step 4: Export the final buffer duplicate/original node pairs for thermal use

**What changes**

- Once Step 3 has identified and rerouted the relevant buffer cut nodes, store the resulting duplicate/original node pairs in a compact sidecar accessible to the thermal setup.
- The minimal exported object is a list of `Node*` or stable ids for pairs that need `T_dup = T_org`.
- If Step 1 duplicates an entire buffer-adjacent interface layer, this exported pair list should cover every duplicate/original pair that must remain thermally continuous. The magnetic cut rewrite itself still only touches the cut-intersection subset.

**Why**

This provides the handoff required by **Prerequisite 3**.

**Depends on**

Step 3.

**Risk**

Medium.

**Estimated effort**

0.5 day.

**Notes**

- Keep this data compact and explicit. Do not infer it later from stale node-index arrays.

## Step 5: Add the thermal pre-kernel continuity hook and fix block selection

**What changes**

- Fix the easy gap first: include `DomainType::Buffer` in the two thermal block-selection switches at `src/fem/thermal/cl_ThermalFactory.cpp:69-85` and `src/fem/thermal/cl_ThermalFactory.cpp:88-104`.
- Add a thermal-only pre-kernel hook that consumes the node-pair list from Step 4 before `mThermalKernel` is constructed at `src/fem/thermal/cl_ThermalFactory.cpp:135`. This hook must be in place before the thermal kernel collects the mesh hanging basis at `src/fem/kernel/cl_FEM_Kernel.cpp:103-107`.
- Preferred implementation direction: a DofManager/DofData-side registration path that creates thermal node constraints `T_dup = T_org` during thermal initialization without mutating the magnetic mesh basis state.
- Keep the later thermal activation path unchanged; it already treats Buffer blocks as `T_phi` domains. [See `src/fem/thermal/cl_ThermalFactory.cpp:172-190`, `src/fem/thermal/cl_IWG_MaxwellThermal.cpp:175-180`]

**Why**

This satisfies **Prerequisite 3** and preserves thermal continuity across the buffer split introduced for magnetics.

**Depends on**

Step 4.

**Risk**

High. The thermal and magnetic kernels share the same mesh object, so accidental cross-talk through shared basis state is the main danger.

**Estimated effort**

1-2 days.

**Notes**

- A late post-`mThermalField->initialize()` patch is not sufficient. By then both the kernel hanging-basis graph and the thermal DOF-wise hanging matrices have already been established.

## Step 6: Verification and regression matrix

**What changes**

- Add a targeted verification checklist for:
  - no-buffer regression,
  - magnetic cut reroute correctness,
  - thermal continuity.

**Why**

This confirms all three prerequisites without broadening the scope.

**Depends on**

Steps 1-5.

**Risk**

Low.

**Estimated effort**

0.5-1 day.

**Verification targets**

- **No-buffer baseline:** a tape without a buffer must remain bit-identical. This is the main regression gate for Step 1.
- **Single buffer tape:** the number of abstract/current DOFs must remain unchanged at one imposed current coefficient. [confidence: high]
- **Current distribution:** HTS carries the imposed transport current; buffer carries zero current; substrate response is solver-determined and expected to be near zero in the insulated limit. [confidence: high for qualitative behavior; medium (~75%) for exact tolerance]
- **Magnetic field:** compare against the no-buffer reference with a tolerance, not bit equality. [confidence: high]
- **Thermal continuity:** all exported duplicate/original buffer pairs satisfy `T_dup = T_org` in the coupled run. [confidence: high]
- **One debug trace:** inspect one cut/shell crossing and confirm the sign/orientation before hardcoding the reroute rule. [confidence: high]

---

## Remaining Gotchas

- `create_buffers()` still detects buffers by material-name substring (`"buffer"`). That is acceptable for this plan but remains brittle. [confidence: high; `src/mesh/cl_ThinShellFactory.cpp:1689-1699`]
- If the new node-duplicate path duplicates full buffer interfaces, the thermal continuity hook must re-merge the full thermal pair set, not only the one cut-crossing node. [confidence: medium (~80%)]
- If the global Maxwell pipeline is ever reordered so that thin-shell blocks exist before `CutFactory::run()`, the reroute hook could migrate inward. In the current code it cannot. [confidence: high]

