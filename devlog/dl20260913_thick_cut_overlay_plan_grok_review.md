# Devlog 2026-09-13 — Thick-cut overlay plan, independent Grok review

**Date:** 2026-09-13
**Purpose:** Second technical review of `todo/thick_cut_overlay_block.md`, after Codex's
same-day pass and Christian's all-or-nothing ruling. Source re-opened; Codex's F1–F7 were
not taken as given.
**Module:** `src/homology`, `src/mesh`, `src/fem/kernel`, `src/fem/maxwell`
**AIs involved:** Grok grok-4.6
**Verification:** static source review on this worktree. Literature: Pellikka et al. 2013
eq. 4.5 (`literature/papers/topology/pellikka2013.txt:365-372`) and Alves et al. 2022b eq. 10
(`literature/papers/fem/alves2022b.txt:320-323`). No build, no `make check`, no solver run.
Reviewed, not verified.

## Verdict

The overlay block is the right carrier (O17). The 2026-09-09 plan is **not** an implementation
checklist in its present form. Christian's 2026-09-13 ruling (all thin if the attempt
succeeds, all thick if it fails) already kills the hybrid skip surgery that R2 was built
around; Codex F1–F7 are the other reasons not to start from that R2 text.

This pass independently confirms those seven findings and adds nine comments the first review
did not make load-bearing: three dof-count formulas, the material-map mechanism, when overlay
edges actually exist, MaxwellData at two Calculator sites, distributor aura vs `gNoIndex`,
`DomainType` integer priority, one-block sparse Z, the terminal-alignment snapshot, and a
default-change gate. Live checkboxes are C8–C16 in the plan §9.

No implementation authorization is inferred. The closed-core exit remains owner-applied
(protocol §7.1).

## Independent confirmation of Codex F1–F7

| ID | Claim | This pass | Evidence re-opened here |
|---|---|---|---|
| F1 | Decide globally before thin topology commits; R2 vs R3 sequencing is inconsistent | **agree, high** | `CutProcessor` ctor `collect_facets` at `:71`, then `duplicate_nodes`/`relink_elements`/`create_thin_cut_sidesets`/`finalize` at `:84-92`. Factory calls the processor then `InterfaceProcessor` (`cl_CutFactory.cpp:460-473`). Overlay copies taken before `InterfaceProcessor` do not follow a later parent relink (`cl_InterfaceProcessor.cpp` relinks by physical master-block ids). |
| F2 | `collect_elements` runs the thin-only classifier before the non-tight abort | **agree, high** | `determine_cut_case_2d/3d` at `cl_CutData.cpp:333-354`, seven 3D patterns then abort at `:600-608`, coefficients from ±1 bitsets at `:231-248`. The plan's "inadmissible pattern stays an abort" is true for the thin path and fatal for an unrectified thick generator that reaches this function. |
| F3 | The infeasible branch is not a ten-line swap at the `BELFEM_ERROR` | **agree, high** | Flag/level mutation and `error.exo` at `cl_Cohomology.cpp:496-532` with the comment that abort makes corruption acceptable; abort at `:534-544`; normal cleanup at `:629-641`. A return at `:534` skips that cleanup. |
| F4 | Canonical original-id pairs ≠ live edge orientation | **agree, high** | `Mesh::compute_edge_directions` compares element-local endpoints to the existing `Edge` endpoints (`cl_Mesh.cpp:1844-1860`). `EdgeFactory::edge_key` orders by **current node index**, not original id (`cl_EdgeFactory.cpp:327-336`). |
| F5 | `AnnulusForcedRepair` is a repairable perturbation | **agree, high** | `tests/homology/test_Cohomology.cpp:255-260,301-340`: add an exact coboundary, then `clean()` must restore a unit representative. That is the opposite of SPFA infeasibility. |
| F6 | Buffer layers hang their own edges after the overlay would have been built | **agree, high on the omitted producer; medium on the physics** | `ThinShellFactory::create_buffers` (`:2049-2103`) types non-resistive layers `Buffer` and hangs edges on their own nodes. Buffer dispatches with Air (`cl_IWG_Maxwell.cpp:286-310`). `tapestack3d` has a magnesia buffer. |
| F7 | `>= dim` does not drop a closed TET4 cochain; §2 sign conflicts with R6 | **agree, high** | Closed `z_ab = p_b - p_a` on a complete simplex has at least `dim` nonzero edges. `phi()` assembles `μ0 ∫ BᵀB` (`mt_maxwell_phi.cpp:41-45`); for `h = −Bφ + Ψ I` the coupling is `−BᵀΨ`, as R6 writes. Pellikka eq. 4.5 uses `+ φ dni`; Alves eq. 10 uses `−φ ∇w`. BELFEM's `compute_air` is `H = −B φ` (`cl_MaxwellPostprocessor.cpp:471-477`). R1 is the sign gate; do not derive it. |

Disagreement of substance with Codex: none. Additional items below are omissions, not refutations.

## Additional findings (G1–G10)

### G1 — Three dof-count formulas, not one (P1 for R5)

`DofData::num_dofs_per_element` (`cl_FEM_DofMgr_DofData.cpp:3439-3455`) is what R5 names.
`IWG::number_of_dofs_per_element` (`cl_IWG.cpp:509-519`) is a parallel sum over node/edge/face
tables and is the count `SolverData` asserts against in DEBUG (`cl_FEM_DofMgr_SolverData.cpp:1021-1026`).
`Element::link_dofs` computes `mNumberOfDofs` the same way and `malloc`s that many pointers
(`cl_FEM_Element.cpp:169-174`) **before** any trailing abstract-dof link. Updating only the
first formula is a debug abort on the first overlay element.

### G2 — Material map does not see internally created blocks (P1 for R3)

`collect_material_labels_from_domains` walks `mDomains` (`cl_MaxwellFactory.cpp:2758-2814`).
Overlay blocks are not input domains. A new `switch` case there never fires. Assign
`mMaterialBlockAssignment[id] = "air"` after the block exists. The `select_blocks` /
`link_to_group` / activation-mode tables remain a separate, still-required set (gap row 9).
Unknown types in the `:834-848` switch become `Inactive`; unknown types in `select_blocks`
(`:773-791`) are silently omitted.

### G3 — Overlay edges do not survive the factory (P1 for O8/R6)

`CutFactory::run` resets edges and faces after the processor (`cl_CutFactory.cpp:172-174`).
`MaxwellFactory::create_edges_and_faces_on_mesh` rebuilds edges only for Conductor groups
(`:1588-1597`). O8(a) is an addition to **that** call, plus `Block::set_edges_flag()` on the
overlay.

`Calculator` sets `tHaveH` from `tEquation->all_fields()` (`cl_FEM_Calculator.cpp:977-987`),
which is IWG-global. HPhi already lists `edge_h`, so any overlay block with `has_edges()` gets
an `EdgeFunction`. That is convenient for O8(a) and is also why MaxwellData must be skipped
explicitly: `allocate()` constructs it for every non-Air block (`:1233-1247`) and
`link_maxwell()` returns early only for `DomainType::Air` (`:1740-1742`). A new type is not Air.

`TimestepMatrices::assemble_J` is `J = M + K Δt` (`cl_TimestepMatrices.cpp:116-124`). Writing
ΔM into `M()` is the correct magnetodynamic slot; Codex's `:741` citation is `collect_qhist`,
which is the history that multiplies against that M. The claim is right, the line is nearby.

### G5 — Thin-shell owner inheritance is not a complete MPI precedent (P1 for R4)

Thin shells have their own nodes. Overlay copies **share** air node objects, so they appear in
air node-stars.

- Kernel graph construction already skips neighbours with `index() == gNoIndex`
  (`cl_FEM_Kernel.cpp:275-276`). Dropping overlay from the metis graph is safe **there**.
- Distributor aura does not skip (`cl_Mesh_Distributor.cpp:411-420`). `DynamicBitset::set` is
  a debug assert and a release out-of-range write (`cl_DynamicBitset.hpp:524-529`).
  `gNoIndex` is `numeric_limits<index_t>::max()` (`src/core/typedefs.hpp:57`).
- After `update_element_indices` the overlay likely has a valid mesh index, in which case the
  failure mode is extra ghosts / possible double assembly, not a crash. Either way it is L-11.
  Enumerate ignore-vs-ghost for overlay in every node-star walk before R4.

The plan's thin-shell "side authority" citations do not walk `node->element`.
`cl_ThinShellFactory.cpp:2270` is a comment about `is_closed()`. `compute_side_authority`
(`:2633-2653`) reads the facet master/slave. Row 6 should not treat those sites as the blast
radius.

### G7 — Append `DomainType`; do not insert (P1 for R3)

`fix_facet_masters` compares domain-type integers as master/slave priority
(`cl_MaxwellFactory.cpp:1531-1538`). The header comment on Air/Coil/Ferro/Conductor (`en_DomainType.hpp:23-25`)
is that rule. `Topology` sizes `mTypeMap` as `UNDEFINED+1` (`cl_Topology.cpp:168`). Inserting a
value in the middle renumbers every later type. Append before `UNDEFINED`. `to_string`'s
default is `"Undefined"` (`en_DomainType.cpp:93-96`).

### G8 — One overlay block with sparse Z (O14, recommendation)

All-or-nothing means every generator is thick together. Per-block dof counts are uniform
(`num_dofs_per_element` is per block), so the cheap mesh is **one** overlay block whose local
dof list is `[φ, I_1, …, I_{n_k}]` with n_k = number of generators and Z structurally zero
where a generator does not cross the element. That assembles ∫μ ψ_g·ψ_h once. Per-set blocks
(plan row 21 option a) still work; they are extra bookkeeping that v1 no longer needs in
order to mix realizations.

### G4 / C14 — Snapshot after terminal alignment

`updatekGeneratorsFromHomology` builds new cochains as linear combinations, deletes the old
ones, and replaces the container (`cl_Cohomology.cpp:1189-1208`). The factory then calls
`clean()` again (`cl_CutFactory.cpp:364-366`). Overlay Z is that final basis. A status flag
taken inside the constructor `clean()` does not label the generators the factory will hand
to `CutProcessor`. Combined with F3: the enabled policy must be known at construction, and
the snapshot must be after the second clean.

### G10 — Default `auto` is a behaviour change (O2, L-04)

A deck that omits the new key and currently dies on the thin-cut abort would, under default
`auto`, run the overlay and produce a solution. That is a default change. Gate it on a deck
that omits the key. Production default `thin` (today's abort) plus a test-only force-thick
hook keeps the existing contract. Christian's ruling is the dispatch policy, not automatically
the omitted-key default.

## What the 09-09 plan got right and should keep

- Additive assembly with zero φφ on the overlay; parents stay in the air block.
- Abstract node as the shared I dof; cell dofs are the wrong carrier
  (`cl_FEM_DofMgr_DofData.cpp:566-600`; `cl_FEM_Element.cpp:226`).
- `IWG_Maxwell::collect_abstract_node_dofs` does **not** fix the dofs (override at
  `cl_IWG_Maxwell.cpp:647-659`); `set_currents` fixes current-driven ones (`:96-121`);
  voltage rows stay free and `Controller::impose_voltage_bcs` adds `V Δt` (`:514-545`).
  Overlay does not change that chain. Row 24 is correct.
- Hidden-block infrastructure is unfinished: no `Block::hide()`, Exodus filters on count
  only (`cl_Mesh_ExodusWriter.cpp:127-142`), ProtoMesh restores hidden for sidesets only
  (`cl_ProtoMesh.cpp:1425` vs `:1333-1342`).
- Volume coupling does not impose the conductor trace; `create_hanging_edges_and_facets`
  sources only the two air nodes (`cl_MaxwellFactory.cpp:1694-1707`) and the T-matrix edge
  conversion is bounded by `number_of_nodes()` (`cl_FEM_DofMgr_DofData.cpp:3562-3575`).
- `CutFactory::cuts()` is never populated (`:2332`); do not hang registration on it.
- `TimestepMatrices::assemble_J` = `M + K Δt` makes the M slot the right place for ΔM.
- Config tag: the file is `src/fem/maxwell/fn_mesh_config_tag.hpp:339-347`, not under
  `src/mesh/`. The schema warning at `doc/input_schema.yaml:1187-1196` is the right one.

## Literature

Alves et al. 2022b eq. 10 is the BELFEM-shaped statement: `h = Σ h_e w_e − Σ φ_n ∇w_n + Σ I_i ψ_i`.
Pellikka et al. 2013 eq. 4.5 is the same construction with the opposite nodal sign convention
(`+ φ dni`) and with integer cohomology coefficients unrestricted to ±1. Integer `|z_e| ≥ 2` is
in-scope for the thick path; tightness is not required. Neither paper specifies BELFEM's overlay
block; that is a carrier choice (O17), and Christian already chose it.

## Timeline comment

The 11–14 working-day figure was calibrated on mixed-mode skip logic that is now struck.
Removing that surgery helps R2. The remaining surgeries are the core-exit lifecycle (F3),
the thin-classifier bypass (F2), the interface current term (R7, historically defect-prone),
E-on-φ-block plus MaxwellData skip (R6/G3), and distributor aura/owner (R4/G5). Re-estimate
after C1–C16 are folded into the recipe. Devlog throughput is not this plan's critical path.

Serial-first (current-driven TET4, one rank, `|z_e|=2` unit test) is still the right floor.

## Files updated

- this devlog; `devlog/README.md`
- `todo/thick_cut_overlay_block.md` §9 C8–C16 and the scope-clarification pointer
- `todo/README.md` pointer
- `tmp/ai_exchange/thick_cut_plan.md` (ephemeral)

No `src/` modified. R0 remains unticked.
