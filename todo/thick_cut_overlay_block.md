# Thick Cuts as an Overlay Block: Fallback When the Thin Cut Cannot Be Built

**Date:** 2026-09-09
**Purpose:** When the thin-cut attempt fails (non-tight cohomology representative, or a generator
with no unit representative), do not abort: impose the currents with the cocycles themselves,
h = −∇φ + Σ_g I_g ψ_g, ψ_g = Σ_e z_{g,e} w_e (Pellikka et al. 2013 eq. 4.5; Alves et al. 2022b
eq. 10, `literature/papers/fem/alves2022b.txt:320-323`). Mechanism: **one** hidden **overlay block**
of copies of every φ-domain element crossed by any generator (same node objects, new element ids)
whose element dof list is [φ at the nodes, I_1 … I_n of the abstract nodes] and whose IWG assembles
only the φ–I coupling and the I–I Gram block. The air block is untouched; no node is duplicated.
**Module:** `src/homology` (dispatch, generator record), `src/mesh` (block, Exodus, `.bfm`,
distributor), `src/fem/kernel` (dof slot, element linking, partition), `src/fem/maxwell` (kernel,
interface, postprocessing)
**AIs involved:** Claude (blind pass + plan), Codex gpt-6-astra/high (blind pass), Grok grok-4.6/high
(blind pass); prior round 2026-09-06 (`tmp/thickcut/`, Codex terra/high + Grok high) reconciled in
`devlog/dl20260909_thick_cut_overlay_plan.md`; plan review round 2026-09-13 (Codex F1–F7, Grok
G1–G10, both high effort) folded in by Claude the same day.
**Status:** PLAN — drafted 2026-09-09 from three blind passes and the 2026-09-06 handoff. Second
iteration 2026-09-09: difficulty grades and a timeline. **Third iteration 2026-09-13: rebased on
Christian's all-or-nothing ruling and on the Codex and Grok plan reviews** — the per-generator
hybrid recipe is struck, the thin attempt is split from its commit, the core exit is specified as a
lifecycle contract rather than a ten-line swap, the coefficient record is oriented explicitly, the
buffer-layer trace is a third interface site, and the gates test the production dispatcher.
**Pending Christian's approval; no source modified.** Gregory Giard's standing ruling (2026-09-03)
is "eventual target, not now"; Christian is considering doing it now. One change sits in the closed
core (O13: a policy-driven status exit in `clean_spfa`, AI drafts, Christian or Gregory applies).
Not scheduled until Christian says so.

**Scope ruling — Christian, 2026-09-13:** always attempt thin cuts first. If that attempt fails,
use thick cuts for **all generators in the model**. There are no mixed thin/thick cases.
Per-generator diagnostics ("generator 3 was the one that failed, and why") stay; per-generator
realization does not.

> **Scope guards:**
> - **Cohomology core untouched.** `cl_Cohomology`, `cl_Homology`, `cl_SimplicialComplex`,
>   `cl_Chain`, `cl_Cochain`, `fn_Smith` are closed to AI edits; every step below lives outside them.
>   Nothing here changes how a generator is computed, only how it is turned into a constraint.
> - **Two triggers, both in v1.** (i) The non-tight abort in `CutProcessor::collect_facets`.
>   (ii) The non-unit abort, which fires *inside* the core's constructor path (`Cohomology::clean()`
>   from both ctors, `cl_Cohomology.cpp:36-40, :58-61`, `clean_spfa` at `:189-192`, reached from
>   `cl_CutFactory.cpp:344` and again `:364-366`). It fires only after SPFA has proven that no unit
>   representative exists (`:534-544`), i.e. on a valid integer cocycle — exactly the object the thick
>   path consumes. The change is a status exit instead of the abort, taken **before** the
>   diagnostic mutation at `:496-532` and with the normal cleanup at `:629-641` preserved (O13).
>   AI drafts, Christian or Gregory applies. **Non-unit generators are the whole point of the thick
>   cut, not a blocker.**
> - **Thin cut stays the first attempt; the decision is global and is taken before any thin
>   topology is committed.** No duplicate node, relinked element or thin sideset from a partially
>   successful thin attempt survives a fallback (F1).
> - **v1 element scope: TET4 air, 3D.** Second order is excluded by the interface hanging code
>   itself (`create_hanging_edges_and_facets` aborts for order > 1, `cl_MaxwellFactory.cpp:1625`).
>   Ferro in the band is excluded (O4). 2D keeps today's abort in v1 (O16).
> - **Parents are never removed** from the air block. The overlay adds stiffness; it must never
>   re-assemble the φφ block.
> - The 2024 XFEM sub-element prototype (`tmp/thickcut/recovered/`) is **not** on this path
>   (three-voice verdict 2026-09-06: dead end as a solver feature; keep for visualisation only).
> - Handoff step ids R1–R7 in `tmp/thickcut/HANDOFF.md` map onto the steps here as noted in §4.

---

## 1. Current Behaviour and How It Fails

The thin-cut pipeline runs inside the `CutProcessor` **constructor** (`cl_CutProcessor.cpp:21-93`,
no `run()`), called from `CutFactory::compute_thin_cuts_and_duplicate_interface_nodes()`
(`cl_CutFactory.cpp:452-484`) on rank 0 only (`MaxwellFactory::create_cuts_sub_master`,
`cl_MaxwellFactory.cpp:1213-1247`; workers run an empty factory, `:1288-1302`). The constructor
order matters for the fallback: `collect_edges` (`:50`) → `collect_elements` (`:52`, which runs the
**thin-only local classifier** `determine_cut_case_2d/3d`, `cl_CutData.cpp:333-354`) →
`collect_facets` (`:71`, the non-tight aborts) → `collect_nodes` → `determine_cut_sets` →
`create_cut_sets` → bitsets → `create_abstract_nodes` (`:82`) → `duplicate_nodes` (`:84`) →
`relink_elements` (`:86`) → `collect_duplicates` → `create_thin_cut_sidesets` (`:90`) →
`mMesh->finalize()` (`:92`). Everything up to and including `collect_facets` is a **feasibility
attempt** that mutates only flags and indices; everything after it **commits** thin topology.

| Failure | Mechanism | Evidence | Fallback? |
|---|---|---|---|
| Non-tight sheet — 3D, no conjugated face survives | positive-case faces do not close | `cl_CutProcessor.cpp:583-588` | **yes (v1)** |
| Non-tight sheet — 3D, pruning eats > ¼ of the emitted faces | idem (corc_solder class, measured 2026-09-02) | `cl_CutProcessor.cpp:589-594` | **yes (v1)** |
| Non-tight sheet — 2D, no conjugated edge survives | idem | `cl_CutProcessor.cpp:384-392` | status exit, but the 2D thick builder is O16 (not v1): still aborts, with the fallback named in the message |
| No unit representative (`\|z_e\| ≥ 2`) | `clean_spfa` aborts after SPFA proves infeasibility | `cl_Cohomology.cpp:534-544` (closed core), called from `cl_CutFactory.cpp:344, :364-366` | **yes (v1)** — policy-driven status exit, owner-applied (O13) |
| Inadmissible local edge pattern | generator or closedness defect upstream | `cl_CutData.cpp:600-608` (3D), `:413-416` (2D) | **no** — stays an abort on the thin path; the thick path never calls this classifier (F2) |
| Accidental path — non-tet air element | one element type for the whole mesh | `cl_CutProcessor.cpp:33-35`; `BELFEM_ERROR( false, "This should not happen" )` at `:683, :802, :1220` | O6 |

Everything the fallback needs is alive at the non-tight abort and dies right after it: the
`Cohomology` object (deleted with the factory, `cl_CutFactory.cpp:79-81`), per generator the cocycle
edges with periodic mates (`CutData::collect_edges`, `cl_CutData.cpp:189-219`), the crossed elements
(`collect_elements`, `:264-305`). One loss already happens before the abort: `collect_coefficients`
keeps only ±1 as two bitsets (`cl_CutData.cpp:231-248`). Mesh edges and faces exist during the
factory (`CutFactory::run`, `:119-128`) and are **reset afterwards** (`:172-174`); the FEM later
recreates edges for Conductor blocks only (`cl_MaxwellFactory.cpp:1588-1597`). On a `.bfm` reload
the whole cohomology stage is skipped (`cl_MaxwellFactory.cpp:697-720`), so whatever the overlay
needs must be in the file.

**Bottom line:** both aborts sit in constructors, and the data that would let a fallback run is
transient and already lossy. The fallback needs (i) a global verdict "thin feasible / not" taken
before the thin topology is committed, (ii) a persisted generator record with full integer
coefficients, keyed on something that survives the edge reset and oriented explicitly, and (iii) a
carrier for the global current dofs inside volume elements — which is what the overlay block is.

## 2. Architecture: Why the Overlay Block Is the Right Spine

- **Assembly is additive across blocks.** An overlay element with dofs [φ_1..φ_n, I_1..I_m] and a
  local matrix whose φφ sub-block is zero adds exactly the coupling block M_φI = −μ0 ∫ BᵀΨ and the
  Gram block M_II = μ0 ∫ ΨᵀΨ to the global system, with Ψ = E Z (sign of the coupling as in R6;
  the reference-tetrahedron probe in the Codex review gives −1/6 μ0 for the node-1/current entry
  with B the gradient operator of `phi()`, `mt_maxwell_phi.cpp:34-45`; **R1 still owns the sign
  of Z itself**). Symmetry is preserved because both the row and the column are assembled by the
  same element. The air block's `maxwell::phi*` kernels stay untouched.
- **One overlay block, sparse Z (O14 resolved).** Under the all-or-nothing rule every generator is
  thick together, and per-block dof counts are uniform, so the overlay is **one** block whose
  element dof list is [φ, I_1 … I_n] with n = number of generators and Z structurally zero where a
  generator does not cross the element. That assembles the cross term ∫μ ψ_g·ψ_h once, for free.
  Per-generator-set blocks (the 09-09 design) still work; they were bookkeeping for a mixed
  realization that no longer exists.
- **Per-block dof tables are the existing idiom.** `DofData::num_dofs_per_element` is per block
  (`cl_FEM_DofMgr_DofData.cpp:3439-3455`), `FieldList::collect_block_dofs` is a switch on the block's
  `DomainType` (`cl_Maxwell_FieldList.cpp:239-295`), and `IWG_Maxwell::link_to_group` picks the kernel
  by `DomainType` (`cl_IWG_Maxwell.cpp:257-440`). A new block type slots into all three. The
  alternative — augmenting the air kernel in place on tagged elements — needs per-element dof counts
  that do not exist and would touch every air kernel; it is the same integrals (Grok), so the choice
  is a carrier choice, recorded as O17 (Christian chose the overlay).
- **The abstract node already is the global dof.** It carries a `phi`-type dof
  (`cl_IWG_Maxwell.cpp:211`; flags `cl_FEM_DofMgr_DofData.cpp:1067-1094`), is selected on every rank
  (`cl_Mesh_Distributor.cpp:395-398`), is fixed for current-driven runs (`IWG_Maxwell::set_currents`,
  `cl_IWG_Maxwell.cpp:100-122`) and free for voltage-driven ones (`cl_FEM_Controller.cpp:514-545`;
  `collect_abstract_node_dofs` does not fix it, `cl_IWG_Maxwell.cpp:647-659`), is written to the csv
  (`:4228-4274`) and saved in the `.bfm` (`cl_Mesh_BfmFile.cpp:340-348`). Nothing in that chain changes.
- **Rejected:** the geometric sub-element cut (needs a global side labelling, which a non-tight
  band does not have — Grok 2026-09-06, accepted); the multi-sheet thin cut (after condensation it
  is this design with more bookkeeping — `tmp/thickcut/docs/cut_representative_options.md` §2).

## 3. Gap Table

Classes: **(a)** rebuildable from existing data, **(b)** open question, **(c)** must be handled explicitly.

| # | State / behaviour | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| 1 | Global verdict "thin feasible / not" **before** thin topology is committed | dispatch | no — the ctor runs attempt and commit back to back (`cl_CutProcessor.cpp:50-92`); the three non-tight `BELFEM_ERROR`s are inside `collect_facets` (`:384-392, :583-594`) | c | split the ctor: attempt (edges, elements, facets) returns a status; `commit_thin_cuts()` runs `collect_nodes` … `finalize`. On failure the attempt's tentative edge indices (`cl_CutProcessor.cpp:163-171`) and element indices (`cl_CutData.cpp:315-321`) and flags are reset. The verdict *is* the prune outcome turned into a status (no other tightness query exists — Grok). Per-generator diagnostics logged; realization is global (F1) |
| 2 | Deck key selecting the realization policy | user control, `.bfm` invalidation | no | c | read next to `algorithm` in `create_cuts` (`cl_MaxwellFactory.cpp:1163-1180`); **must** be added to `src/fem/maxwell/fn_mesh_config_tag.hpp:338-347` (the schema warns exactly this, `doc/input_schema.yaml:1187-1196`); `doc/input_file_reference.md:678`; not a new `CutAlgorithm` value (`en_CutAlgorithm.hpp:21-28` selects the reduction, `cl_CutFactory.cpp:280-353`). Default is a **behaviour change** if `auto` (O2, G10) |
| 3 | Crossed-element list = Whitney support | overlay membership | partly — `collect_elements` keeps elements with ≥ dim flagged edges (`cl_CutData.cpp:294-304`) **and** runs the thin classifier on each (`:333-354`) | c | the thick builder extracts support itself: every φ-block element with ≥ 1 edge in supp(z), **without** `determine_cut_case_*` (F2: a closed integer cochain such as `[1,1,−2,3,2,1]` on one TET4 is pattern 63 and aborts the classifier). The ≥ dim threshold does *not* lose support for a closed cochain (F7: on a full simplex z_ab = p_b − p_a has ≥ dim nonzero edges); threshold 1 is chosen to decouple from thin assumptions, not to repair a loss. Closedness is validated separately |
| 4 | Signed **integer** z_{g,e} per support edge, **live-edge oriented** at assembly time | kernel, interface, postproc | no — `collect_coefficients` keeps ±1 bitsets only (`cl_CutData.cpp:231-248`); `CutData::weight()` returns −1/0/+1; the raw map is `Cochain::getSimplicesMap()` (`:189`) | c | **Record = edge key + coefficients (Christian, 2026-09-13).** At snapshot time (edges still exist) store per support edge the pair (`edge->node(0)->original()->id()`, `edge->node(1)->original()->id()`) and z_1..z_n relative to that order. The pair is the edge's identity in id space and survives the reset; note `EdgeFactory::edge_key` keys on current node **indices** (`cl_EdgeFactory.cpp:319-336, :477-490`), so the rebuild resolves the pair through the element's `get_nodes_of_edge`, not through the factory key. **Orientation (F4):** the rebuilt edge runs lower index → higher index, which need not match the stored order, so on rebuild compare the stored pair with the live edge's endpoints, the `compute_edge_directions` idiom (`cl_Mesh.cpp:1844-1860`): match keeps z, reversal negates it. `EF_TET4` then applies local → live (`cl_EF_TET4.cpp:94-95`). Two distinct operations, one comparison each, not a double flip. `original()` keeps periodic mates and interface duplicates on their original's key |
| 5 | Overlay block creation, element ids, shared nodes | mesh | precedent only | a | thin-shell pattern: fresh elements, `append( mMesh->elements(), block->elements() )` **and** `mMesh->add_block` (`cl_ThinShellFactory.cpp:232-246`); `add_block` alone only pushes + maps (`cl_Mesh.cpp:2284-2288`); `unfinalize()` clears `mElements` (`:909`) so a `finalize()` after it recollects from all blocks (`:783-786`); ids from `max_element_id()` (`cl_CutProcessor.cpp:875`); `Block` **owns and deletes** its elements (`cl_Block.cpp:34-38`) so copies must be fresh objects; a block's element type is its first element's (`cl_Block.hpp:126-134`) → one block per element type (TET4 only in v1, so one block) |
| 6 | Node→element back-links now include the copies | everything walking `node->element(e)` | side effect | c | **Enumerated (G5):** kernel metis graph (`cl_FEM_Kernel.cpp:270-303`) already skips neighbours with `index() == gNoIndex` (`:275-276`) — dropping the overlay from the graph is safe there; distributor aura (`cl_Mesh_Distributor.cpp:405-420`) does **not** skip, and `DynamicBitset::set` is debug-asserted only (`cl_DynamicBitset.hpp:524-529`; `gNoIndex = numeric_limits<index_t>::max()`), so an overlay with `gNoIndex` is a release OOB write and one with a valid index is an extra ghost. The thin-shell "side authority" sites cited on 09-09 do **not** walk `node->element` (`cl_ThinShellFactory.cpp:2270` is a comment; `compute_side_authority` `:2633-2653` reads facet master/slave) — struck from the blast radius. Decision per walk: ignore vs ghost (R4) |
| 7 | Block hidden from Exodus | output | no — `Block::mIsHidden` exists (`cl_Block.hpp:50`, getter `:167`) but **has no setter** and is documented VTK-only; `ExodusWriter::collect_blocks` filters on element count only (`cl_Mesh_ExodusWriter.cpp:114-140`); sidesets honour `is_hidden()` (`:57, :183, :462`) | c | add `Block::hide()`; filter both passes of `collect_blocks`; header `mNumElements` (`:170-174`) and field loops (`:922-1015`) iterate `mBlocks` |
| 8 | Hidden flag transported | MPI, `.bfm` | half — `.bfm` packs the bit into the group tag (`cl_ProtoMesh.cpp:1270`, `cl_Mesh_BfmFile.cpp:466`, load `:525`) but `ProtoMesh` restores it for sidesets only (`cl_ProtoMesh.cpp:1425`; blocks `:1333-1342` do not); distributor block data is id/etype/dtype/thickness (`cl_Mesh_Distributor.cpp:1777-1833`) | a | derive `hidden` from the new `DomainType` on every path |
| 9 | New `DomainType` for the overlay block | dof table, kernel, postproc, names, activation | no — `Cut = 8` is the thin-cut *sideset* type (`cl_MaxwellFactory.cpp:1242`; `cl_Maxwell_FieldList.cpp:383`); typing the overlay `Air` **doubles the air stiffness** via `link_to_group` | c | **Append the enumerator before `UNDEFINED`, never insert (G7):** `fix_facet_masters` compares the integer values as master/slave priority (`cl_MaxwellFactory.cpp:1531-1538`, the header comment `en_DomainType.hpp:23-25`), and `Topology` sizes `mTypeMap` as `UNDEFINED+1` (`cl_Topology.cpp:168`). Tables: `en_DomainType.hpp` + `to_string` (default `"Undefined"`, `en_DomainType.cpp:93-96`); `FieldList::collect_block_dofs` **default aborts** (`cl_Maxwell_FieldList.cpp:288-293`); `link_to_group` (`cl_IWG_Maxwell.cpp:286-310` is the Air case); `MaxwellPostprocessor` ctor (`cl_MaxwellPostprocessor.cpp:56-132`); **`select_blocks`** (`cl_MaxwellFactory.cpp:773-793`, unknown types silently omitted) and the activation switch at `:829-848` (unknown types → `Inactive`); `mBlockTypes` broadcast is automatic once the rank-0 map holds it (`cl_IWG_Maxwell.cpp:144-160`); labels free from `set_block_and_sideset_names` (`:1252-1275`) |
| 10 | Material for the overlay block | Calculator, kernel | no — `collect_material_labels_from_domains` walks the **input** `Domain` objects (`cl_MaxwellFactory.cpp:2758-2814`); the overlay is created internally and is never a deck domain, so a new `switch` case there never fires (G2) | c | assign `mMaterialBlockAssignment[overlay_id] = "air"` after the block exists |
| 11 | Overlay block excluded from the φ topology | a second cohomology run must not see doubled elements | n/a | c | `Topology::phi_block_ids` (`cl_Topology.hpp:88-91`); `IWG_Maxwell::custom_postprocess` flags Air/Buffer nodes only (`cl_IWG_Maxwell.cpp:590-600`), fine |
| 12 | Generator record persisted | restart, `.bfm` reuse | no cohomology record in the `.bfm` (abstract/pinned/orphaned `:340-370`, duplicates `:1289-1370`); reload skips cohomology (`cl_MaxwellFactory.cpp:697-720`) | c | §6; malformed-record checks are release-active `BELFEM_ERROR`s (sizes, generator bounds, missing parents) |
| 13 | Generator record distributed | MPI | no | c | `share`/`receive` pair (CLAUDE.md: never `broadcast` for variable-size payloads); elements arrive grouped by block index (`cl_Mesh_Distributor.cpp:1099-1113`); symmetric send/receive including ranks with zero overlay elements; the record must reach every rank that runs recovery on a band parent or sources an interface edge, not only overlay owners |
| 14 | Overlay element owner = parent's owner | MPI | precedent exists for thin shells: their blocks are dropped from the metis graph (`cl_FEM_Kernel.cpp:226-232`) and inherit the facet owner (`:448-460`); owners are reset at `:189-193` and finalized by `update_ownerships` at `:516` | c | same pattern: drop overlay from the graph, copy owner from the parent before `:516`; needs a parent link on the copy (O7). The thin-shell precedent covers ownership only, **not** the shared-node aura (row 6): shells have their own nodes, the overlay shares air nodes |
| 15 | Abstract dofs visible on the ranks that assemble the band | MPI | no — dof→proc flags are derived from geometric element nodes of selected blocks (`cl_FEM_DofMgr_DofData.cpp:1176-1190`); receiving the abstract node ≠ receiving its dof (Codex) | c | mark every abstract dof consumed by every rank owning an overlay element (one block, all generators) |
| 16 | +n dofs per overlay element — **three count formulas, one malloc** | dof manager | no — `DofData::num_dofs_per_element` (`cl_FEM_DofMgr_DofData.cpp:3439-3455`), `IWG::number_of_dofs_per_element` (`cl_IWG.cpp:509-519`, the count `SolverData` asserts against in DEBUG, `cl_FEM_DofMgr_SolverData.cpp:1021-1026`) and `Element::link_dofs` (`cl_FEM_Element.cpp:169-174`, sizes `mNumberOfDofs` and `malloc`s `mDOFs`) are three parallel sums of nodes·dpn + edges·dpe + faces·dpf; cell dofs exist but are **per-element unique** (`:566-600`) — the wrong carrier; `link_dofs` has the placeholder "here is where we would link #celldof" (`:226`) | c | add a per-block "abstract dof" count to **all three** (G1) — updating one is a debug abort on the first overlay element; n = number of generators |
| 17 | Linking the abstract dofs into the element | dof manager | precedent only — `Element::link_lambda_dofs` appends a facet's lambda dof (`cl_FEM_Element.cpp:1292-1313`, needs `mFacet`); `link_dofs( aDofManager, blockID )` is the block-element entry (`:143`) | c | new trailing step after the malloc is sized for it; the abstract dofs stay free (non-hanging) unknowns in the fixed order of `iwg->abstract_nodes()`; local Jacobian sized from the new count (`cl_IWG_Maxwell.cpp:280-282` sizes from `number_of_local_dofs()`; `Calculator::qold()/q()` traverse the actual local list, `cl_FEM_Calculator.cpp:2894-2915, :2995-3000`) |
| 18 | `Calculator::E(k)` on a φ block | kernel | no — the edge function is built only if the IWG-global fields hold `edge_h` (`cl_FEM_Calculator.cpp:977-987`, true for HPhi) **and** `block->has_edges()` (`:1041-1049`); it needs `Element::edge_directions` (`cl_EF_TET4.cpp:94-95`); air blocks have no edges after the factory reset | c | **O8 resolved (G3):** overlay edges are created in `MaxwellFactory::create_edges_and_faces_on_mesh` (`cl_MaxwellFactory.cpp:1588-1597`, today Conductor groups only) by adding the overlay block to the `create_edges` call, plus `Block::set_edges_flag()` (`cl_Block.hpp:175`) on the overlay; no edge dofs. `MaxwellData` must be skipped for the new type at **both** `allocate()` (`:1233-1247`) and `link_maxwell()` (`:1740-1742`) — both skip Air only, and the new type is not Air |
| 19 | Overlay kernel | assembly | no | c | new `maxwell::phi_thick_cut` next to `phi` (`mt_maxwell_phi.cpp:26-45`), writes the same `aMatrices->M()` slot: `TimestepMatrices::assemble_J` is `J = M + K Δt` (`cl_TimestepMatrices.cpp:116-124`) and `collect_qhist` (`cl_IWG_Timestep.cpp:741`) multiplies the history against it, so magnetodynamics and the timestep history see ΔM like air; φφ block zero; for a fixed I the column is a RHS contribution, for a free I it is live |
| 20 | Interface edges in the support: h_e = φ_n0 − φ_n1 + Σ_g I_g z_{g,e} — **three sites** | conductors, thin shells, buffers | no — sources are the two air nodes (`create_hanging_edges_and_facets`, `cl_MaxwellFactory.cpp:1690-1710`); the T-matrix edge branch converts node sources in a loop bounded by `edge->number_of_nodes()` (`cl_FEM_DofMgr_DofData.cpp:3562-3575`, weights `:3600-3680`), so simply appending a third source is not enough; thin shells hang in `hang_thinshell_edges_on_nodes_bottom/top` (`cl_MaxwellFactory.cpp:1901, :1972`); **buffer layers** hang their in-plane edges on their own nodes in `ThinShellFactory::create_buffers` (`cl_ThinShellFactory.cpp:2049-2103`), and layer nodes copy the *authority* node's sources and weights (`:2790-2816`) — in the thin path that copies the cut duplicate's abstract-node source into every layer, in the thick path there is nothing to copy (F6, confirmed 09-13) | c | explicit "extra current term" path in the conversion; weight z_{g,e} (precedent: `CutSet::create_duplicates` sources the abstract node with weight 1, `cl_CutSet.cpp:42-80`); once per shared edge; periodic slave chains preserved (`:1672-1674` uses `original()`); buffer in-plane edges carry the z of the surface edge they project onto (O18 for the physics check) |
| 21 | Interface node duplication for thick generators | consistency | `InterfaceProcessor` runs after the processor (`cl_CutFactory.cpp:466-476`), ties duplicates with weight 1 (`cl_InterfaceProcessor.hpp:27-31`) and **relinks** air-side master elements of the interface sidesets to duplicates (`cl_InterfaceProcessor.cpp:95-141, :231-246`) | c | overlay copies are built **after** the `InterfaceProcessor`, so `copy->node(i) == parent->node(i)` holds at partition and serialization time (F1); asserted in the builder (O9) |
| 22 | **Overlapping generators** (two bands share elements) | correctness | n/a | c | one block, sparse Z, dof list [φ, I_1 … I_n]: the kernel assembles the full ΔM with Ψ = E Z, cross terms included (Codex critical → O14 resolved by G8) |
| 23 | H reconstruction on band elements | postproc and the solve | no — `compute_air` H = −Bφ (`cl_MaxwellPostprocessor.cpp:471-479`); `l2_phi` (`mt_maxwell_l2_phi.cpp:41, :69`); `IWG_Maxwell::custom_postprocess` (`cl_IWG_Maxwell.cpp:586-600`); **shell trace** `compute_h_trace` φ-side branch (`cl_FEM_Calculator.hpp:2694-2702, :2742-2746`) feeds normal-field recovery during the solve, not only output | c | add + Σ I_g ψ_g on the **physical** element (nodal fields live on shared nodes; the overlay must not run recovery/projection or it doubles the mass) — O10 |
| 24 | `CutFactory::cuts()` as a registry | — | **never populated** (only `return mCuts`, `cl_CutFactory.cpp:2332`), so the `DomainType::Cut` loop in `create_cuts_sub_master` (`cl_MaxwellFactory.cpp:1240-1243`) is a no-op (Grok) | a | do not hang overlay registration on it; note for a separate cleanup |
| 25 | Current / voltage bookkeeping on the abstract dof | controller | yes, unchanged | a | `cl_FEM_Controller.cpp:514-545, :3042, :4228-4274` |
| 26 | Snapshot timing of Z | record | n/a | c | the basis handed to `CutProcessor` is the **post-alignment** one: `updatekGeneratorsFromHomology` builds new cochains as linear combinations, deletes the old ones and replaces the container (`cl_Cohomology.cpp:1189-1208`), then `clean()` runs again (`cl_CutFactory.cpp:364-366`). Snapshot after the last `clean()`; a verdict taken in the constructor `clean()` does not label the final basis (C14) |
| 27 | Infeasibility fixture | gates | none — `AnnulusForcedRepair` (`test_Cohomology.cpp:255-340`) adds an exact coboundary to a unit generator and expects `clean()` to **repair** it; it is the opposite of SPFA infeasibility (F5) | c | keep it for integer-preservation and direct builder tests; add a distinct fixture whose class has no unit representative (owner-approved, since only the core can certify that) — R9 |
| 28 | Tests | gates | none for `CutProcessor` or the solve; `tests/homology` checks generators + pairing (`test_Cohomology.cpp:214-240`, `fast`); `tests/fem/test_InterfaceOrientation.cpp` has the Ampère telescope idiom (`:270-280`) | c | §4 R9 tests the **production dispatcher**, not only a forced thick |

### 3.1 Cross-cutting findings

- **Edge identity does not survive the pipeline.** Everything keyed on `Edge*` or edge id at the
  abort site is gone after `reset_edges()`; the persisted record is keyed on original node-id pairs
  and **converted to the live edge orientation on rebuild** (row 4).
- **The thin classifier runs before the tightness test** (row 3). A generator that is a valid
  integer cocycle but not a unit one would abort in `determine_cut_case_3d` before `collect_facets`
  ever runs. On the thick path the support is extracted by the builder, never by the classifier.
- **Shared node objects have a blast radius** (row 6), now enumerated: one safe walk (kernel
  graph), one unsafe walk (distributor aura), two false citations struck.
- **Five tables key on the block type and abort or silently misbehave on an unknown one**
  (rows 9, 10, 16): dof table (abort), kernel dispatch (Air = doubled stiffness), `select_blocks`
  (silently not assembled), activation switch (`Inactive`), material map (never reached).
- **Volume coupling alone does not impose the current on conductors** (row 20). The overlay fixes
  Ampère in the air; the conductor trace still reads −∇φ until the interface relation carries
  + I z_e — at three sites, the third being the buffer layers that today inherit I by copying a
  duplicate's sources.
- **The interface path is linear-only already** (`cl_MaxwellFactory.cpp:1625`), fixing the v1 element
  scope to TET4 without a new restriction.

## 4. Ordered Steps

Handoff mapping: R1 ↔ handoff R1, R2 ↔ R2, R3+R4+R5 ↔ R3, R6 ↔ R4, R7 ↔ R5, R8 ↔ R6, R9 ↔ R7.
Struck sub-steps are kept visible per `todo/plan_template.md`.

- [ ] **R0 — Approval and sign-off.** Christian approves the plan; Gregory informed (his ruling
      2026-09-03 stands unless he lifts it). Decided: O3 (model-wide fallback), O8 (edges via
      `create_edges_and_faces_on_mesh`), O14 (one block, sparse Z), O17 (overlay). Still to decide
      before R2: O2 (default), O13 (form and owner of the core exit), O16 (2D out of v1), O18
      (buffer traces).
- [ ] **R1 — Orientation spike (read-only, no new paths).** On tapestack3d (healthy thin cut,
      `examples/tapestack3d/input.conf`): read the stored cochain z, compute ⟨z, γ⟩ for the paired
      homology cycle, compare with the imposed I and with ∮H from the saved solution. Fixes O1 (sign
      of z_e in h_e = φ_n0 − φ_n1 + I z_e). Also counts overlapping generators on a corc-class mesh
      (informational). Gate for everything after.
- [ ] **R2 — Global verdict and dispatch.** (after: R1)
  - [ ] **Split `CutProcessor` into attempt and commit (F1).** The ctor keeps `collect_edges`,
        `collect_elements`, `collect_facets` and stores a status instead of the three non-tight
        `BELFEM_ERROR`s (`cl_CutProcessor.cpp:384-392, :583-594`); a new `commit_thin_cuts()`
        runs `collect_nodes` … `create_thin_cut_sidesets` … `finalize()` (`:73-92`). On failure
        the attempt's flags, the rewritten edge indices (`:163-171`) and element indices
        (`cl_CutData.cpp:315-321`) are reset. Abstract nodes are created on both paths
        (`create_abstract_nodes`, `:846-854`, preserving the existing ordering). The
        inadmissible-pattern abort (`cl_CutData.cpp:600-608`) **stays** on the thin path.
  - [ ] ~~skip thick generators in `duplicate_nodes/relink_elements/create_thin_cut_sidesets`
        and `CutSet` creation~~ — struck 2026-09-13: no mixed realization.
  - [ ] **Core status exit (owner-applied, O13).** Behind a policy known **at construction**
        (a ctor argument, since both ctors call `clean()`; a later setter cannot rescue an abort
        inside the ctor): when the thick fallback is enabled and SPFA proves infeasibility,
        branch **before** the diagnostic mutation at `cl_Cohomology.cpp:496-532` (it unflags
        candidates, flags obstruction edges, rewrites element levels and writes `error.exo` on the
        stated assumption that the process is about to abort), record "no unit representative",
        stop rectifying the basis (model-wide fallback makes further repair pointless), and fall
        through to the normal cleanup at `:629-641` (`unflag_everything`, progress). The factory
        reads the status after the **last** `clean()` (`cl_CutFactory.cpp:364-366`, row 26). AI
        drafts the diff and a test on a **genuinely infeasible** fixture (row 27; the owner
        certifies infeasibility); Christian or Gregory applies (§7.1).
  - [ ] **Thick support extraction (F2).** New `ThickCutBuilder` (or a `CutData` method that the
        thin path never calls): for every generator, snapshot the raw integer map
        (`Cochain::getSimplicesMap()`) including periodic mates with the same expansion rule the
        thin path uses (`cl_CutData.cpp:249-258`, re-checked against endpoint correspondence,
        F4), collect every φ-block element with ≥ 1 edge in the support, and validate
        **closedness** (zero circulation on every face) and the cycle pairing — never
        `determine_cut_case_*`. Gate: the closed six-edge TET4 example `[1,1,−2,3,2,1]` is
        accepted; a deliberately non-closed cochain is rejected.
  - [ ] `CutFactory::compute_thin_cuts_and_duplicate_interface_nodes` (`cl_CutFactory.cpp:452-484`):
        thin feasible → `commit_thin_cuts()` and today's flow; otherwise `link_node_duplicates…`
        (no-op), `InterfaceProcessor` (still needed, row 21), **then** the thick builder makes the
        overlay block, then `collect_orphan_nodes`. Log which generator(s) failed and why.
  - [ ] Deck key `homology { cut : thin | thick | auto }` (name: O2) in `create_cuts`
        (`cl_MaxwellFactory.cpp:1163-1180`), passed into `CutFactory` **and** into the `Cohomology`
        ctor policy; `src/fem/maxwell/fn_mesh_config_tag.hpp:338-347`; `doc/input_schema.yaml:1178+`;
        `doc/input_file_reference.md:678`. `thin` = today's abort. The **omitted-key default is a
        separate decision** (O2): `auto` turns a fatal deck into a running one, so it needs a run
        gate on a deck that omits the key (L-04). `thick` forces the overlay without a thin attempt
        (test and diagnostics; whether it is a production value is O2).
  - [ ] ~~`CutData::collect_elements` threshold 1 for the thick path~~ — superseded by the
        builder's own extraction (F2, F7).
- [ ] **R3 — Overlay block on the mesh (rank 0).** (after: R2)
  - [ ] New `DomainType` **appended before `UNDEFINED`** (`en_DomainType.hpp`, G7) + `to_string`;
        `Block::hide()` setter; hidden derived from the type on the `.bfm` and distributor paths
        (row 8).
  - [ ] Build **one** block: fresh copies of every crossed φ element, same node objects (asserted
        `copy->node(i) == parent->node(i)`), ids from `max_element_id()`, parent id kept on the copy
        (O7); built after the `InterfaceProcessor` and before the factory's
        `unfinalize/reset_edges/reset_faces/finalize` (`cl_CutFactory.cpp:172-175`) — the edges
        still exist then, so Z can be taken per element and local edge; registration via the
        thin-shell pattern (`append` + `add_block`, `cl_ThinShellFactory.cpp:232-246`);
        `Topology::phi_block_ids` excludes the type.
  - [ ] Add the type to `select_blocks` (`cl_MaxwellFactory.cpp:773-793`) and the activation switch
        at `:829-848`; set `mMaterialBlockAssignment[overlay_id] = "air"` after creation (row 10;
        O4 refuses Ferro).
  - [ ] Row 6: overlay dropped from the metis graph (safe, `:275-276`); the distributor aura walk
        (`cl_Mesh_Distributor.cpp:411-420`) gets an explicit rule (ignore, or ghost with the
        parent) — decided here, implemented in R4.
  - [ ] Exodus: `is_hidden()` filter in both passes of `ExodusWriter::collect_blocks`
        (`cl_Mesh_ExodusWriter.cpp:114-140`).
- [ ] **R4 — Generator record: persist and distribute.** (after: R3)
  - [ ] `.bfm`: one support-edge table (n0, n1, z_1..z_n by original node id, row 4) plus the
        parent id per copy (§6), under the block; the same table serves the volume kernel and the
        three interface sites; load path rebuilds without the `Cohomology` object (reload skips
        it, `cl_MaxwellFactory.cpp:697-720`); per-element Z is assembled at link time from the
        element's local edges (`get_nodes_of_edge` → table lookup via a setup-time `Map` on the
        pair → endpoint comparison with the live edge); malformed-record checks release-active
        (unknown pair on a support edge is an error, not a zero); round-trip test.
  - [ ] Distributor: `share`/`receive` of the record, symmetric on every rank including empty local
        sets; copies land with their parents (drop from the metis graph and copy the owner,
        thin-shell pattern `cl_FEM_Kernel.cpp:226-232, :448-460`, before `update_ownerships` `:516`);
        the aura rule from R3 implemented; the record reaches every rank that recovers on a band
        parent or sources a band interface edge.
  - [ ] Abstract dofs consumed on every rank owning an overlay element
        (`cl_FEM_DofMgr_DofData.cpp:1176-1190`).
- [ ] **R5 — Dof slot and linking.** (after: R3)
  - [ ] `FieldList::collect_block_dofs` case for the new type (φ node dofs only,
        `cl_Maxwell_FieldList.cpp:239-295`).
  - [ ] Per-block abstract-dof count n in **all three** formulas (G1): `DofData::num_dofs_per_element`
        (`cl_FEM_DofMgr_DofData.cpp:3439-3455`), `IWG::number_of_dofs_per_element` (`cl_IWG.cpp:509-519`),
        `Element::link_dofs` sizing (`cl_FEM_Element.cpp:169-174`); local Jacobian/timestep buffers
        sized from it; the DEBUG identity in `SolverData` (`:1021-1026`) is the gate.
  - [ ] `Element::link_dofs` trailing step "link abstract dofs" modelled on `link_lambda_dofs`
        (`cl_FEM_Element.cpp:1292-1313`), in `iwg->abstract_nodes()` order; the abstract dofs stay
        free unknowns.
- [ ] **R6 — Coupling kernel.** (after: R4, R5)
  - [ ] Overlay edges: add the overlay block to the `create_edges` call in
        `create_edges_and_faces_on_mesh` (`cl_MaxwellFactory.cpp:1588-1597`) and `set_edges_flag()`
        on it; `MaxwellData` skipped for the type at `allocate()` (`cl_FEM_Calculator.cpp:1233-1247`)
        **and** `link_maxwell()` (`:1740-1742`); confirm `E(k)` is built (`:1041-1049`) with no edge
        dofs and no `nedelec_data_h()` call.
  - [ ] `maxwell::phi_thick_cut`: ΔM = μ0 ∫ [ 0, −BᵀΨ ; −ΨᵀB, ΨᵀΨ ] (sign of Z per R1), Ψ = E Z with
        Z the n_edges × n integer matrix in live orientation, structurally sparse; dispatch case in
        `link_to_group` (`cl_IWG_Maxwell.cpp:257-440`); same `M()` slot; scratch for Ψ and ΨᵀΨ
        sized once at link time, no per-call `Matrix` temporaries.
  - [ ] Probe (L-01/L-09), with a precise observable: zero Z in the **volume kernel only** and
        measure the missing magnetic energy / residual on the band; a prescribed I stays fixed and
        an intact interface term can still impose circulation, so "the current is lost" is not the
        observable. Then restore.
- [ ] **R7 — Interface relation, three sites.** (after: R5)
  - [ ] Conductor edges in the support: extra current terms (abstract nodes, weights z_{g,e}) in
        `create_hanging_edges_and_facets` (`cl_MaxwellFactory.cpp:1690-1710`) and an explicit path in
        the T-matrix edge conversion (`cl_FEM_DofMgr_DofData.cpp:3562-3680`); once per shared edge;
        periodic chains preserved.
  - [ ] Thin-shell top/bottom edges: same rule in `hang_thinshell_edges_on_nodes_bottom/top`
        (`cl_MaxwellFactory.cpp:1901, :1972`).
  - [ ] **Buffer-layer in-plane edges** (F6): same rule for the edges `create_buffers` hangs on
        their own nodes (`cl_ThinShellFactory.cpp:2049-2103`), with z taken from the surface edge
        the layer edge projects onto; the normal-direction edges carry zero. O18 settles whether
        the cocycle's restriction to the shell surface is the right extension (it is what the thin
        path does implicitly through the copied duplicate sources, `:2790-2816`); until then
        tapestack3d's magnesia buffer is **not** covered by the two outer helpers alone.
  - [ ] Guard O5: an interface edge reached from both sides with opposite sign → refuse.
- [ ] **R8 — Postprocessing and output.** (after: R6)
  - [ ] H, B on band elements gain + Σ I_g ψ_g on the physical element: `compute_air`, `l2_phi*`,
        `custom_postprocess`, and the shell trace `compute_h_trace` (`cl_FEM_Calculator.hpp:2694-2702`,
        which is solve-side) — mechanism per O10; overlay excluded from recovery/projection.
  - [ ] φ written continuous; overlay block hidden (R3) with an optional deck switch to show it.
- [ ] **R9 — Gates.** (after: R8) The gates test the **production dispatcher**, not only a forced
      `thick` (F5).
  - [ ] Unit (element algebra): assembled ΔM against direct quadrature of the enriched basis on one
        TET4, signed coefficients incl. ±2, two overlapping generators (exact cross term 1/6 μ0 on
        the reference tet, Codex probe); `fast`.
  - [ ] Unit (orientation, F4): reversed live edge, permuted node storage, periodic/interface
        duplicates, `.bfm` reload — reconstructed Ψ and signed circulation unchanged; `fast`.
  - [ ] Unit (support extraction, F2): closed `[1,1,−2,3,2,1]` accepted, non-closed rejected; `fast`.
  - [ ] Unit (pairing): overlay build on the annulus meshes of `tests/homology`, ∮ I ψ = I via the
        `inner_ring_pairing` idiom (`test_Cohomology.cpp:240`); `fast`.
  - [ ] Dispatcher: healthy model → all thin, no overlay block; non-tight model → all thick with
        evidence the fallback ran; SPFA-infeasible model (row 27 fixture) → all thick; two
        generators where an early thin candidate succeeds and a later one fails → **no** committed
        thin topology, no thin-attempt duplicates, abstract-node ordering unchanged; invalid
        cocycle → rejected, not silently thick.
  - [ ] tapestack3d: forced `thick` vs thin to solver tolerance on H, loop circulation, magnetic
        energy, terminal quantities, transient losses (raw φ need not agree); ∮H = I at every
        station — **after** O18 is settled for its buffer.
  - [ ] **Both drive modes:** current-driven (I fixed, column is RHS) and voltage/circuit-driven
        (I free, φ–I and I–I live) — a prescribed-current match alone hides a missing I–I block (O15).
  - [ ] The previously aborting corc_solder deck (nonfree) runs on `auto` and still aborts on
        `thin`: winding = I along the cell, tape J not rim-only; `examples/corc_periodic_bc`.
  - [ ] MPI: np 2 and 4 on tapestack3d thick, incl. a rank with no overlay elements and a band that
        spans partitions; overlay owner == parent owner; abstract dofs present where assembled;
        serial-equal within solver tolerance.
  - [ ] `.bfm` round trip (coefficients and identity preserved, cache invalidated by the key) and
        restart identical; Exodus has no overlay block and air H includes I ψ.
  - [ ] Omitted-key deck (L-04): whatever O2 decides, one deck without the key runs and behaves as
        documented.
  - [ ] `make check-fast`.
  - [ ] A generator with |z_e| = 2 imposes the right current — the case the thin cut cannot
        represent at all and the sharpest single proof that the overlay does something new. The
        unit test may seed from `AnnulusForcedRepair` (integer preservation through the builder,
        with `clean()` **not** called in between, or it repairs the coefficient away); the
        end-to-end run uses the infeasible fixture / a coarse corc-class deck.
- [ ] **R10 — Docs.** (after: R9) `src/homology/doc/thick_thin_cuts_and_conjugate_edges.md` gains the
      ψ section; `cohomology_theory_and_implementation.md` §FEM Discretization; module README;
      Codex language sweep; devlog; `doc/input_*` already in R2.

### 4.0 Implementation Progress

Not started. This subsection becomes the Dn tracker once R0 is approved.

### 4.1 Difficulty grading — third iteration (2026-09-13)

Grades: **M** = mechanical, one plan+audit → code+audit round, AI-led, Christian reviews the diff.
**G** = guided, AI-led but a Christian decision or a run gate sits inside the loop. **S** = surgery,
Christian leads, AI audits and writes tests. A *round* is one plan+audit → code+audit cycle; in this
repository that is a single session (same-day precedents: 2026-07-03, 2026-07-16, 2026-08-04, §4.2).
"Run" = Christian builds and runs (the user runs builds, CLAUDE.md).

| Step | Item | Grade | Rounds | Why this grade |
|---|---|---|---|---|
| R0 | approval, O2/O13/O16/O18 decided | — | 0 | decisions, not work |
| R1 | orientation spike on tapestack3d | G | 1 + run | AI writes the probe; the answer needs a run and Christian's reading of the sign |
| R2 | `CutProcessor` attempt/commit split, status, tentative-state reset | **S** | 1 | Christian's construction; smaller than the struck hybrid skip (no per-generator threading through `relink_elements`), but the reset of tentative indices/flags is the part to get right once |
| R2 | core status exit in `clean_spfa` (owner-applied) | **S** | 1 + owner | not ten lines: ctor policy, early branch before `:496`, fall-through to `:629-641`, status readable after the last `clean()`; AI drafts + infeasible-fixture test, Christian/Gregory applies |
| R2 | thick support extraction + closedness validation | G | 1 | new code, no thin coupling; the periodic-mate expansion rule is the one thing to check against the thin path |
| R2 | `CutFactory` dispatch to the thick builder | M | (same round) | one call site, `cl_CutFactory.cpp:452-484` |
| R2 | deck key + ctor policy + config tag + schema + reference | M | 1 | precedent: `homology.algorithm`; input-contract rule is mechanical; O2 decides the default |
| R3 | `DomainType` (appended), `Block::hide()`, hidden derived on `.bfm`/distributor paths | M | 1 | enum + setter + three `if`s |
| R3 | build the one overlay block, parent link (O7), node-identity assert, register before the factory reset | G | 1 | thin-shell registration precedent; placement after `InterfaceProcessor` is now fixed |
| R3 | `select_blocks`, activation switch, material assignment | M | (same round) | two switch cases + one map write |
| R3 | aura rule decision (row 6) | G | (same round) | one decision, one walk |
| R3 | Exodus `is_hidden()` filter | M | (same round) | two loops |
| R4 | `.bfm` edge table (§6) + link-time Z assembly with endpoint comparison + load without `Cohomology` + round-trip test | G | 1–2 | precedent: side-connector record 2026-08-04; the comparison is the `compute_edge_directions` idiom and has its own unit test |
| R4 | distributor `share`/`receive` + owner inheritance + aura rule | **S** | 2 + np gates | MPI ownership (L-11); the thin-shell owner pattern exists but the shared-node aura is new |
| R4 | abstract dofs consumed on overlay-owning ranks | G | (same round) | one flag pass |
| R5 | `FieldList` case | M | (with R5 below) | one switch case |
| R5 | abstract-dof count in three formulas + buffer sizing | G | 1 | three parallel sums; the DEBUG assert is the gate |
| R5 | `link_dofs` trailing "link abstract dofs" | G | 1 | precedent `link_lambda_dofs` |
| R6 | overlay edges via `create_edges_and_faces_on_mesh`, `set_edges_flag`, two `MaxwellData` skips | G | 1 + run | was graded S while the E source was open; O8 is resolved and the sites are known; the run confirms `E(k)` exists with no edge dofs |
| R6 | `phi_thick_cut` kernel + dispatch | M | 1 | 30 lines of quadrature once R1 fixed the sign |
| R6 | inert-kernel probe | M | (same round) | zero Z in the volume kernel, measure the energy defect |
| R7 | conductor interface: extra current term in `create_hanging_edges_and_facets` and the T-matrix edge conversion | **S** | 2 | the most defect-prone region in the record (hanging-edge orientation devlogs 2026-04-07, 04-24, 08-25); the conversion loop bound means a representation change, not an append |
| R7 | thin-shell top/bottom edges | **S** | 1 | same rule, second site; after the conductor path is green |
| R7 | buffer-layer in-plane edges | **S** | 1 | third site; needs O18 first |
| R7 | O5 guard | M | (same round) | detect opposite signs on one edge, refuse |
| R8 | + Σ I ψ in `compute_air`, `l2_phi*`, `custom_postprocess`, `compute_h_trace` | G | 1–2 | four sites, one helper; O10 decides the mechanism; the shell trace is solve-side so it needs a run |
| R8 | φ continuous, hidden block, show switch | M | (same round) | already true once R3 lands |
| R9 | four `fast` unit tests (algebra, orientation, extraction, pairing) | M | 1–2 | pure comparisons; `tests/fem/test_EdgeFunctions.cpp` and `inner_ring_pairing` idioms |
| R9 | dispatcher gates incl. the infeasible fixture | G | 1 + owner | the fixture needs the owner to certify infeasibility |
| R9 | tapestack3d thin vs thick, both drive modes | G | runs | Christian runs; AI writes the comparison script |
| R9 | corc_solder (nonfree), periodic corc | G | runs | Christian runs |
| R9 | MPI np 2/4 | G | runs | Christian runs; this is where R4's surgery is proven |
| R9 | `.bfm` round trip, restart, Exodus, omitted-key deck | M | (with R4 test) | |
| R10 | docs + sweep + devlog | M | 1 | |

**Totals:** ~22 rounds, of which **6 are surgeries** (R2 attempt/commit split, R2 core exit,
R4 distributor, R7 × 3 sites); the rest are AI-led rounds Christian reviews. Every run gate is
his. Against the 09-09 count: the hybrid-skip surgery is gone, the E-on-φ-block surgery dropped to
guided, and two surgeries were added (the core exit is no longer "trivial"; the buffer site is
new).

### 4.2 Calibration from the devlog

| Comparable landing | Scope | Elapsed | Source |
|---|---|---|---|
| `.bfm` I/O backend rewrite (ProtoMesh), 25 defects, serial + np 2/4 gates | mesh + persistence + MPI, larger than this plan's R3+R4 | 2026-06-22 → 07-01, **9 days** | `todo/closed/meshfile_refactor_plan.md` |
| Side-connector `.bfm` record | one persistence record + audit + 3 follow-ups | **same day**, 2026-08-04 | `devlog/dl20260804_sideconnector_bfm_persistence.md` |
| Interface node duplication (InterfaceProcessor R3/R5) | homology-side processor change + serial validation | **same day**, 2026-07-03 | `todo/closed/interface_node_duplication_coil_ferro.md` |
| Periodic cap cut emission | CutProcessor change + field check | **same day**, 2026-07-16 | `todo/closed/periodic_cap_cut_emission.md` |
| Thin-shell conductor normal field, R1–R5 + six tests | Calculator/postproc surgery in the shell trace | 2026-09-01 → 09-02, **2 days** | `todo/closed/thinshell_conductor_normal_field.md` |
| ParMETIS/PT-Scotch wiring, D1–D5, np 2/4 green | build + MPI partition path | 2026-08-31 → 09-04, **4 days** | `todo/closed/parmetis_ptscotch_wiring.md` |
| Shared library + install spine, Darwin + Linux gates | build system, cross-platform | 2026-08-25 → 09-03, **9 days** | `todo/closed/shared_library_and_install_plan.md` |

The pattern: a plan of this shape here closes in **days, not weeks**, and the elapsed time is set
by the run gates and the surgeries, not by the mechanical rounds. The 2024 thin cut itself went from
first sub-element prototype to node duplication + abstract nodes in five weeks (07-09 → 08-15,
`tmp/thickcut/REPORT.md` §1) — with none of today's infrastructure, and that is the "uninformed"
reference: a developer new to this code would plan **6–8 weeks** for the surface this plan touches.
Both reviewers noted, correctly, that devlog throughput is not this plan's critical path; the
surgeries and the run turnaround are.

### 4.3 Timeline (informed, third iteration)

Critical path: R1 → R2 → R3 → {R4 ∥ R5} → R6 → R7 → R8 → R9. R10 and the four unit tests of R9 can
run alongside R6–R8. R4 and R5 are independent once R3 has landed. Two items sit outside AI
throughput entirely: the owner-applied core exit (R2) and the owner-certified infeasible fixture
(R9); both can be prepared early and land whenever Christian or Gregory has a slot.

| Phase | Steps | Working days | Christian's part |
|---|---|---|---|
| A — spike and dispatch | R1, R2 | 2–3 | reads the sign (R1); owns the attempt/commit split; applies or forwards the core exit |
| B — mesh and carrier | R3, R4 (serial part), R5 | 2–3 | O7; aura rule; `.bfm` record review |
| C — physics, serial, current-driven | R6, R7 (conductor), R8, R9 unit tests + tapestack3d serial | 3–4 | interface surgery; runs |
| D — parallel, shells, buffers, free current, corc | R4 (distributor), R7 (shells, buffers), R9 remaining gates | 4–5 | distributor surgery; O18; np and nonfree runs |
| E — docs | R10 | 1 | sign-off |
| **Total** | | **12–16 working days** ≈ 2.5–3 calendar weeks | |

Up from 11–14: the hybrid skip is gone (−1), but the core exit is a real lifecycle change (+0.5),
the buffer site is a third interface surgery (+1) and the dispatcher gates need an owner-certified
fixture (+0.5). Still "doable within two weeks" only on the serial-first cut below, and only if
phases C and D do not uncover a second defect behind the first (L-10). Phase D is the part most
likely to slip: every "wrong only in parallel" in the record (L-11) cost a day or more.

**Serial-first cut** (if the two weeks must hold): phases A–C give a working, current-driven,
TET4, conductor-only thick cut on one rank in **8–10 days**, including the |z_e| = 2 unit test and
the dispatcher gates; D and E are the following week. tapestack3d's buffer keeps it out of phase C
unless O18 turns out to be "zero support on the buffer" for that deck.

## 5. Open Design Questions (not silently decided)

- **O1 — Sign convention.** Mesh-oriented z_e entering h_e = φ_n0 − φ_n1 + I z_e with or without
  a flip. Settled by R1, not by derivation: Alves eq. 10 uses −φ ∇w (BELFEM's `H = −Bφ`), Pellikka
  eq. 4.5 uses +φ dn_i (`pellikka2013.txt:365-372`); both are self-consistent, the code is one of
  them. Do not apply the `edge_directions` sign twice (Codex); the stored → live endpoint comparison of
  row 4 is a *different* sign and is applied once.
- **O2 — Key name and default.** `homology { cut : thin | thick | auto }`. **The omitted-key
  default is a behaviour change if `auto`** (G10): a deck that dies today on the thin-cut abort
  would run and produce a solution. Options: (a) default `thin` (today's contract; `auto` is
  opt-in; a `thick` value exists for tests and diagnostics); (b) default `auto` with a run gate on
  an omitted-key deck (L-04) and a loud log line naming the fallback. Christian's ruling is the
  *dispatch policy*, not automatically the omitted-key default. Recommendation: (a) for the first
  landing, flip to (b) once the tapestack3d equivalence gate is green. Christian.
- **O3 — Realization policy.** ~~Keep the thin cut where it builds and overlay only the failed
  generators, or `thick` forces the overlay for all?~~ **Resolved by Christian, 2026-09-13:** always
  try thin first; any failure selects thick for the whole model. No mixed realization. Implemented
  as the attempt/commit split in R2.
- **O4 — Ferro in the band.** μ(h) then depends on I (Newton coupling). v1: refuse if any crossed
  element is Ferro; log.
- **O5 — Thin shells reached from both sides** with opposite signs inject 2I. Detect and refuse
  (v1) or restrict z to one side by construction (needs the cleaner, outside AI scope).
- **O6 — Element coverage.** The Calculator's Nédélec dispatch has no PYRA5
  (`cl_FEM_Calculator.cpp:1095-1150`); HEX8/QUAD4 exist. v1 TET4 only; one block per element type.
- **O7 — Parent link on the copy.** Carrier for "parent element id": a new `Element` member, the
  `.bfm` tags dataset (`cl_Mesh_BfmFile.cpp:551`), or an index into the block record (§6). Decide in R3.
- **O8 — Source of E on a φ block.** ~~(a) create edges on the overlay block … (b) evaluate Whitney
  functions from node coordinates inside the kernel~~ **Resolved (G3, 2026-09-13): (a),** placed
  in `create_edges_and_faces_on_mesh` (`cl_MaxwellFactory.cpp:1588-1597`) because the factory
  resets all edges after the cuts (`cl_CutFactory.cpp:172-174`); `set_edges_flag()` on the
  overlay turns `E(k)` on since `edge_h` is IWG-global; `MaxwellData` skipped at both Calculator
  sites. Verify at the R6 run that a block with edges and zero edge dofs passes the dof manager.
- **O9 — `InterfaceProcessor` with thick generators.** It relinks air-side master elements of the
  interface sidesets to duplicates (`cl_InterfaceProcessor.cpp:95-141, :231-246`); the overlay is
  built after it, and the builder asserts `copy->node(i) == parent->node(i)`. Remaining check:
  nothing in it assumes thin-cut duplicates exist (it takes `mAbstractNodes`, which exist on both
  paths).
- **O10 — Where the ψ correction is added in postprocessing.** Parent air element flagged and
  reading the record (preferred: fields live on shared nodes, and the shell trace is solve-side),
  or an accumulate-mode overlay pass. Decide in R8.
- **O11 — Conditioning of M_II on a large band** (corc: ~7.5k of 190k nodes). Measure in R9; a tight
  representative (option 1, min-cost tension, `cut_representative_options.md` §4) shrinks it and is
  independent of this plan.
- **O12 — Second order.** z lives on the lowest-order Whitney part of TET10/TRI6 only; excluded in
  v1 by the interface path anyway.
- **O13 — Form and owner of the core status exit.** The |z_e| ≥ 2 abort is inside `clean_spfa`
  (`cl_Cohomology.cpp:534-544`), called from both constructors and again from the factory
  (`cl_CutFactory.cpp:344, :364-366`). **What the exit has to be (F3, C3):** a policy passed to the
  `Cohomology` ctor (both ctors call `clean()`, so a setter is too late); a branch **before** the
  diagnostic mutation at `:496-532` (which corrupts flags and levels and writes `error.exo` on the
  assumption that the process is about to abort); a recorded "non-unit basis" status; rectification
  stopped for the whole basis (model-wide fallback), pocket removal skipped; the normal cleanup at
  `:629-641` still executed. Whether Gregory applies it himself, and the exact API shape (ctor
  flag + status getter is the minimum), is the owner's decision. AI drafts, does not apply (§7.1).
  Not a blocker for phases B–C, which can be built and unit-tested on the non-tight trigger.
- **O14 — Overlapping generators.** ~~Options (a) one block per generator set … (b) cross-term
  blocks~~ **Resolved (G8, 2026-09-13): one overlay block, dof list [φ, I_1 … I_n], Z structurally
  sparse.** The cross term is assembled once by construction. Frequency of overlap in practice is
  still worth counting on corc in R1, for O11.
- **O15 — Which drive mode first.** Current-driven only (column is RHS; smaller change) or both from
  the start. R9 requires both before "done"; v1 may land current-driven first.
- **O16 — 2D.** The 2D abort exists (`cl_CutProcessor.cpp:384-392`) and TRI3 Whitney functions exist;
  Gregory scoped option 2 as 3D-first. **Recommendation: out of v1.** The status exit is shared
  code and lands for 2D too, but the thick builder is TET4-only, so a 2D non-tight failure still
  aborts — with a message that names the fallback as 3D-only. Christian.
- **O17 — Overlay block vs in-place augmentation of the air kernel.** Same integrals; the overlay
  costs a `DomainType`, Exodus/partition/link handling, and keeps `phi()` untouched; in-place needs
  per-element dof counts in the air block. Christian's brief chose the overlay; recorded so the
  alternative is visible (L-13).
- **O18 — Cocycle on thin-shell buffer layers (F6, new 2026-09-13).** The cocycle is computed on
  the φ-domain mesh before the shell layers exist; its restriction to the shell surface is what
  the top/bottom helpers see. In the thin path, buffer layers inherit I by copying the duplicate's
  sources (`cl_ThinShellFactory.cpp:2790-2816`), i.e. the surface restriction is implicitly
  extruded through the stack. The thick path has to do that explicitly: in-plane buffer edges carry
  the z of the surface edge they project onto, normal edges carry zero. **Check before R7:** is
  that extrusion the formulation Christian intends for a multi-layer shell (same I through every
  layer's trace), and does tapestack3d's cocycle have support on the tape surface at all? If the
  support is zero there for that deck, the tapestack3d gate can run before the buffer site lands.
  Christian.

## 6. Generator Record (to be audited against the code before R4)

One overlay block (element type T = TET4 in v1), n generators. Revised 2026-09-13 on Christian's
proposal to key the record on the edge's own endpoints: one edge table replaces the per-element
coefficient array and the separate interface-edge list.

| field | type | source | consumer |
|---|---|---|---|
| generator indices | `Cell< index_t >` into `mesh->abstract_nodes()`, in the factory's existing order | R2 | dof linking, controller ordering |
| parent element id per copy | `Cell< id_t >` (n_copies) | R3 (O7) | owner inheritance, postproc, node-identity assert |
| support-edge table | per support edge: `id_t n0, n1` = `edge->node(0/1)->original()->id()` in the edge's own order at snapshot time, then `int z_1 … z_n` relative to that order, **full integers**, structural zeros where a generator does not cross | raw `Cochain::getSimplicesMap()` snapshot after the last `clean()`, taken while the factory's edges still exist; periodic mates expanded with the thin path's rule (`cl_CutData.cpp:249-258`, re-checked against endpoint correspondence, F4) | volume kernel (per-element Z built at link time from `get_nodes_of_edge` + lookup + endpoint comparison with the live edge, row 4); all three R7 interface sites (same lookup on the hanging edge); postproc |

`.bfm`: under the block's group, datasets `generators`, `parents`, `edges`. The record must not
need the `Cohomology` object to load, must survive the edge reset, and is validated on load with
release-active errors (sizes, generator bounds, parent existence, and a support edge of an overlay
element that has no table row). Storage is owned `Cell` for transport and `Matrix` for the
algebra, sized at setup; the pair → row `Map` is built once at link time.

## 7. Definition-of-Done Checklist

- [ ] Every gap-table row mapped to a step or an open question (rows 1–28 ↔ R2–R9, O1–O18).
- [ ] Each claimed gap backed by a citation re-opened in the session that cites it.
- [ ] Ordered steps with dependencies.
- [ ] Open questions logged, not decided.
- [ ] Dispatcher gates (healthy, non-tight, infeasible, late failure, invalid) green; tapestack3d
      thin vs thick in both drive modes; corc_solder winds by I; `.bfm` round trip; np 2/4;
      overlapping-generator, orientation and extraction unit tests; |z_e| = 2 unit test and deck;
      omitted-key deck; `make check-fast`.
- [ ] Docs updated and swept; devlog written.

## 8. Audit Trail

- Exchange thread: `tmp/ai_exchange/thick_cut_plan.md` — Claude blind pass, Codex gpt-6-astra/high
  blind pass, Grok grok-4.6/high blind pass (both told not to open `tmp/thickcut/`; both confirmed
  they did not); 2026-09-13 plan reviews (Codex, Grok) and Claude's reconciliation.
- Handoff from the 2026-09-06 session (another machine): `tmp/thickcut/{HANDOFF,REPORT,FINDINGS}.md`,
  `docs/cut_representative_options.md`, `audits/`. Its R1–R7 and O1–O6 are absorbed above.
- Devlogs: `devlog/dl20260909_thick_cut_overlay_plan.md` (blind passes and first plan),
  `devlog/dl20260913_thick_cut_overlay_review.md` (Codex F1–F7),
  `devlog/dl20260913_thick_cut_overlay_plan_grok_review.md` (Grok G1–G10),
  `devlog/dl20260913_thick_cut_overlay_plan_rebase.md` (this reconciliation).
- **Contamination note:** Claude read `tmp/thickcut/HANDOFF.md` (the one-page summary) *before* its
  own pass, by accident; `REPORT.md`, `FINDINGS.md`, the options document and the audits were read
  only afterwards. The handoff's design matches Christian's brief, so the effect is small, but the
  Claude pass is not blind to the R1–R7 outline.

### 8.1 Blind-pass reconciliation (2026-09-09)

All three passes agree on: the three non-tight abort sites; `en_CutAlgorithm` is the wrong axis;
the key must enter the config tag; no `.bfm`/MPI generator record exists; `num_dofs_per_element`
has no scalar slot; `E(k)` is gated on `has_edges`; typing the overlay `Air` doubles the stiffness;
the interface relation needs + I z_e; `compute_air` needs + I ψ; no solver-level cut test exists.

| Finding | Claude | Codex | Grok | Re-opened |
|---|---|---|---|---|
| Non-unit abort is inside the core's ctor path; needs an owner-applied status exit (O13) | missed | ✔ | ✔ | `cl_Cohomology.cpp:36-40, :58-61, :189-192, :534-544`; `cl_CutFactory.cpp:344, :364-366` |
| Overlapping generators need ∫μ ψ_g·ψ_h (O14, row 22) | missed | ✔ (critical) | — | design argument; no code site |
| `collect_coefficients` keeps ±1 only; snapshot the raw map first (row 4) | missed | ✔ | ✔ | `cl_CutData.cpp:231-248` |
| `collect_elements` ≥ dim threshold under-counts the Whitney support (row 3) | ✔ (from handoff) | — | ✔ | `cl_CutData.cpp:296-302` — **refuted 2026-09-13 (F7)** for closed cochains; the row now says why threshold 1 is still used |
| No `Block::hide()` setter; `.bfm` saves the bit but `ProtoMesh` restores it for sidesets only (rows 7, 8) | partial | ✔ | ✔ | `cl_ProtoMesh.cpp:1270, :1333-1342, :1425`; `cl_Mesh_BfmFile.cpp:466, :525` |
| Thin-shell precedent for block registration and owner inheritance (rows 5, 14) | missed | ✔ | ✔ | `cl_ThinShellFactory.cpp:232-246`; `cl_FEM_Kernel.cpp:226-232, :448-460` |
| `select_blocks` list and the activation switch must admit the new type (row 9) | missed | ✔ | ✔ | `cl_MaxwellFactory.cpp:773-793, :829-848` |
| Abstract dof visibility on workers is derived from element nodes (row 15) | missed | ✔ | — | `cl_FEM_DofMgr_DofData.cpp:1176-1190` |
| Cell dofs are per-element unique, wrong carrier (row 16) | missed | — | ✔ | `cl_FEM_DofMgr_DofData.cpp:566-600`; `cl_FEM_Element.cpp:226` |
| `MaxwellData` not allocated for Air; overlay kernel must not need it (row 18) | missed | — | ✔ | `cl_FEM_Calculator.cpp:1740-1742` |
| Edge-source conversion loop is bounded by `number_of_nodes()` (row 20) | missed | ✔ | — | `cl_FEM_DofMgr_DofData.cpp:3562-3575` |
| Shell trace `compute_h_trace` is solve-side, not only output (row 23) | missed | ✔ | — | `cl_FEM_Calculator.hpp:2694-2702` |
| Current-driven vs free I: both gates needed (O15) | missed | ✔ | ✔ | `cl_IWG_Maxwell.cpp:100-122`; `cl_FEM_Controller.cpp:514-545` |
| `CutFactory::cuts()` never populated (row 24) | missed | — | ✔ | `cl_CutFactory.cpp:2332` |
| Node→element back-link blast radius (row 6) | ✔ | ✔ | ✔ | `cl_FEM_Kernel.cpp:270-303`; `cl_Mesh_Distributor.cpp:405-420` — narrowed 2026-09-13 (G5) |
| `Topology::phi_block_ids` must exclude the overlay (row 11) | ✔ | — | — | `cl_Topology.hpp:88-91` |
| `.bfm` reload skips cohomology (row 12) | missed | — | ✔ | `cl_MaxwellFactory.cpp:697-720` |
| `tests/fem/test_InterfaceOrientation.cpp` as a gate foundation (row 28) | missed | ✔ | — | `:270-280` |

Disagreements: none of substance. Grok leans toward asking whether the overlay is the right carrier
at all (O17); Codex accepts it and pushes the overlapping-generator requirement. Neither auditor
opened the handoff, and both independently reproduced its R1–R7 shape from the tree — which is the
strongest evidence in this session that the handoff's plan is grounded.

## 9. Plan review round — 2026-09-13

Reports: [Codex F1–F7](../devlog/dl20260913_thick_cut_overlay_review.md),
[Grok G1–G10](../devlog/dl20260913_thick_cut_overlay_plan_grok_review.md). Both static, both high
effort, no build or run. Claude re-opened every cited site the same day before folding
(`devlog/dl20260913_thick_cut_overlay_plan_rebase.md`); all seventeen findings held. Each comment
below names where it landed in the body; the body is the recipe, this section is the ledger.

- [x] **C0 — Scope ruling recorded.** All thin on success; all thick on failure. Header, O3.
- [x] **C1 — Rebase R2/R3 around the global decision (F1).** R2 is now an attempt/commit split of
      `CutProcessor` with tentative-state reset; the overlay is built after the `InterfaceProcessor`
      and before the factory's edge reset, with the node-identity assert (rows 1, 21; O9).
- [x] **C2 — Bypass the thin-only classifier on the thick path (F2).** Own support extraction with
      closedness validation; the `[1,1,−2,3,2,1]` example is an R9 unit gate (row 3, R2).
- [x] **C3 — Core exit as a lifecycle contract (F3).** Ctor policy, branch before `:496`, cleanup
      preserved, status after the last `clean()`; "ten-line" struck everywhere (O13, R2, §4.1).
- [x] **C4 — Stored vs live orientation (F4).** Record keyed on the edge's own endpoints by
      original id (Christian's proposal), one endpoint comparison on rebuild, an R9 unit gate
      (row 4, §6, R4, O1).
- [x] **C5 — Gates for the production dispatcher (F5).** Healthy / non-tight / infeasible / late
      failure / invalid; `AnnulusForcedRepair` re-scoped; infeasible fixture is owner-certified
      (rows 27–28, R9).
- [x] **C6 — Buffer traces (F6).** Confirmed from `create_buffers` and the layer-node source copy:
      third interface site in R7, physics question O18, tapestack3d gate sequenced after it.
- [x] **C7 — Corrections and re-estimate (F7).** §2 sign aligned with R6; threshold rationale
      rewritten; 2D scoped out of v1 (O16); inert-kernel observable made precise; §4.1/§4.3 redone.
- [x] **C8 — Three dof-count formulas (G1).** Row 16, R5.
- [x] **C9 — Material map (G2).** Row 10, R3.
- [x] **C10 — Overlay edges after the factory reset; `MaxwellData` at two sites (G3).** Row 18,
      O8 resolved, R6.
- [x] **C11 — Aura vs `gNoIndex`; false side-authority citations (G5).** Row 6 rewritten; R3
      decides the walk rule, R4 implements it.
- [x] **C12 — Append the `DomainType` (G7).** Row 9, R3.
- [x] **C13 — One block, sparse Z (G8).** §2, row 22, O14 resolved.
- [x] **C14 — Snapshot after terminal alignment (G4).** Row 26, R2, §6.
- [x] **C15 — Default is a behaviour change (G10).** O2 rewritten with a recommendation; R9 gate.
- [x] **C16 — Rebase the body, not only §9.** Done: row 1, R2, §4.1 and §4.3 no longer describe a
      per-generator skip; struck sub-steps stay visible.

Open after this round, for Christian: O2 (default), O13 (owner and API of the core exit), O16 (2D
out of v1), O18 (buffer extrusion). Nothing else in the two reviews is unaddressed.
