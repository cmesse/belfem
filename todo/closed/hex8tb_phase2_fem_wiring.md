# HEX8TB Phase 2 — FEM Wiring for the Side Connector Wall Element

**Date:** 2026-07-29
**Purpose:** Make the degenerate HEX8TB wall element (4 longitudinal edge dofs) operational: derive and implement its edge function, wire it into dof counting, IWG and
MaxwellFactory dispatch, and verify orientation and physics on a single-tape case. Phase 1
(mesh-level enum/helper coverage) landed 2026-07-29 and is NOT part of this plan.
**Module:** `src/fem/interpolation/nedelec`, `src/fem/maxwell`, `src/fem/kernel`
**AIs involved:** Claude (plan + implementation support), Codex + Grok (audits)
**Status:** OPEN — construction, edge function, and `.bfm` persistence are complete (R1–R3,
D1–D7, B1–B6; B7 deprioritized). `h_side_connector` kernel scaffolding landed 2026-08-05
(Christian) and went through the frozen jury protocol; review fixes were applied the same day,
including the `h = ht + hb + hn` recovery via the master-block/sideset-hop route. D7
(side-node original() normalization) was fixed 2026-08-06; the fusing-experiment gate
(flat no-cut tape, `mFuseEdges` on, vs Norris) is now unblocked; the §6 orientation
gate and the `reset_source_container()` clobbering audit remain open before any
cut/periodic fused run. R5 plumbing P1–P5 are DONE (2026-08-08/09: facet wiring, dof table +
group activation, dispatch, material derivation, `edge coating : on` input keys;
P2 anchors corrected 2026-08-08: the abort points are `FieldList::collect_block_dofs`
and the group-activation switch, NOT the metis-feeding `tBlocks` selection).
**P6 (WIP stop) and R6 (postprocessor lists) are both DONE (2026-08-09)** — the
`edge coating : on` path runs through to assembly. **Next is the smoke run, then R7**;
P8 (doc sync) and P9 (parallel deferreds, out of scope for iteration 1) remain. R4 (geometry
interpolation) DONE 2026-08-08: HEX8 Lagrange reused via case labels, exact-cuboid
metric keeps the derivatives — rationale in `side_coating_wall_element.md` §4.5.
Gate: Christian's smoke run. ~~Reminder: B7 — the `.bfm` checksum ignores
connector settings; delete any cached mesh before the first connector-enabled run.~~
**No longer necessary since 2026-08-11 (DR-21):** the mesh-configuration stamp covers the
connector flag and coating width, so a settings change now rebuilds by itself.

> **Scope guards:**
> - The **coated** fuse algorithm (`h_in` internal dof) is designed but not implemented; the
>   uncoated wall (all side traces fused to g) is the first target. Never tie `h_in` to g
>   (the one genuine overconditioning risk — 3-AI consensus 2026-07-28).
> - **Closed side curves** are out of scope (guarded in `compute_binomial_vectors`); the
>   facet-based sign anchor + threefold demo are a separate campaign (see O4).
> - **Selective connector layers** (`mConnectorsForAllLayers = false`) are out of scope
>   (loud guard in `create()`).
> - Parallel distribution of connector blocks is out of scope for the first iteration.

## Context (what exists after Phase 1)

- Element template `cl_Element_HEX8TB.hpp` (`ElementTemplate<8,8,4,6,1>`): nodes follow the
  plain HEX8 convention; the 4 dof edges are (0,1), (3,2), (4,5), (7,6) so that **all four
  point in +curve direction** (deliberate, decided 2026-07-29, Christian).
- `ThinShellFactory::create()` builds one connector block per side curve
  (`DomainType::LeftCoating` / `RightCoating`), elements assembled edge-driven with a
  per-station orientation flag (`cl_ThinShellFactory.cpp:407-473`). Outer nodes and edges
  hang on their mid-plane sources; weights come later from
  `DofData::create_dofwise_t_matrices_master()`.
- Mesh-level helper coverage (meshtools, to_string, both MeshCheckers, CurvedElementChecker,
  vtktools) landed 2026-07-29. The MeshCheckers **error out** on a negative-volume HEX8TB
  instead of swapping — a node swap would corrupt the slot-tied edges.

## Ordered steps

- [x] **R1 — Derive the HEX8TB edge function** (theory session, literature-first).
  DONE 2026-07-31 (3-AI consultation, unanimous): `E_k = s_k · F_k(η,ζ) · ∇ξ` with
  F_k = ⅛(1±η)(1±ζ) — the earlier sketch `½·s_k·N_k(ξ,ζ)·∇η` had a frame-label mismatch.
  Exact-cuboid metric (block w/d + midpoint L, nodes orient the frame only); constant
  Nabla; C gives j_t ≡ 0, j in the (b,n) plane; h_t = stream function of the cross-flow,
  transfer current between boundary points = Δh_t (= the r′ stencil — the Ampère/boundary
  routes meet there; the volumetric-vs-r′ calibration is fixed in R5). Deliverable landed
  as `src/fem/maxwell/doc/side_coating_wall_element.md` (comprehensive campaign doc)
  instead of a nedelec/doc page.
- [x] **R2 — `fn_num_nedelec_dofs.hpp`**: HEX8TB case → 4. DONE 2026-07-31.
- [x] **R3 — `cl_EF_HEX8TB`** in `src/fem/interpolation/nedelec/` + factory case in
  `cl_EdgeFunctionFactory.cpp`. DONE 2026-07-31 (link/update_nabla/E/C implemented per
  R1 consensus, Codex+Grok audited; field recovery for the material law is designed but
  NOT implemented — see the doc §5).
- [x] **R4 — Geometry interpolation**: DONE 2026-08-08 — both case labels landed
  (`cl_IF_InterpolationFunctionFactory.cpp:198`, `cl_FEM_Calculator.cpp:632`); compile
  gate = Christian's next build. Decision (Christian, same day): HEX8TB reuses the
  plain HEX8 Lagrange functions for nodal interpolation, while all spatial derivatives
  and integration weights stay on the EF's exact-cuboid metric — the committed
  metric/interpolation mismatch is documented as the lesser evil vs the exact
  isoparametric Jacobian in `side_coating_wall_element.md` §4.5 (the exact J is not a
  drop-in: it needs a Piola-transformed basis and re-imports the (L/w)² conditioning;
  the r′ calibration absorbs metric constants anyway). Implementation = two case
  labels: `cl_IF_InterpolationFunctionFactory.cpp:196` (join the HEX8/HEX8TS case
  group) and the Calculator nedelec-data switch (add HEX8TB to the
  PENTA6TS/QUAD4TS/HEX8TS/HEX8 → `nedelec_data_linear_h` group, `cl_FEM_Calculator.cpp:~634`
  — the plan files the latter under R6, same edit class). The dV path needs NO change:
  `dV_hex` already prefers `mEdgeFunction->det_J()` when an edge function is attached.
  (after: R1)
- [ ] **R5 — IWG + MaxwellFactory wiring** — kernel scaffolding `h_side_connector` landed
  2026-08-05 (Christian; jury round `review_h_side_connector`: pre-registered, F1–F11
  triaged, fixes applied; kernel math = M(µ0) + K(ρ CᵀC) with h = ht + hb + hn recovery
  and the ρ(T,‖B‖,β) clamped material law). Remaining PLUMBING, in dependency order:

  ### R5 plumbing checklist (from the 2026-08-05 review insights)

  - [x] **P1 — connector facet wiring.** DONE 2026-08-08 (Christian + Claude; compile
    gate = next build): connector passes added to BOTH
    `BlockData::link_thin_shell_facets_serial` and `_parallel` — the fem wall element
    gets its RECOVERY facet via the id+1 invariant (`mesh->facet( element->id() + 1 )`;
    the id+1 object is the facet's QUAD4 in the hidden sideset, NOT a member of
    `mMesh->elements()`), supplying `index_on_master` for `side_connector_xi_eta`.
    BONUS: the fem master element is wired via `set_master` (facet master → fem block →
    fem element), which unblocks deleting the kernel's per-link "wire the pointer
    directly" re-derivation (fold into P7). Parallel variant ships the connector id
    list rank-0 → all (`share`, zero-length safe), and derives facet + master locally;
    aura copies without a shipped facet are skipped (assembly runs owned elements
    only — full aura facet shipping stays a P9 item). Jury round
    `review_hex8tb_wiring_batch` (2026-08-08) hardened the failure paths per the
    frequency rule (setup runs once → always-active checks): owned-⇒-facet and
    null-facet-master are `BELFEM_ERROR` on both variants; the guard triad
    gained a thickness/material count lock and a positive-thickness floor.
    Codex's "connector blocks never reach the kernel" claim was REFUTED
    (Maxwell path: `add_equation` → `create_field` → `IWG_Maxwell::initialize`
    selects ALL blocks; `KernelParameters` selection feeds metis only).
  - [x] **P2 — dof table + group activation.** DONE 2026-08-08 (Christian): connector
    cases landed in `FieldList::collect_block_dofs` (share the Conductor edge-h dof
    table, `cl_Maxwell_FieldList.cpp:258-259`) and in the group-activation switch
    (`cl_MaxwellFactory.cpp:586-587`); metis `tBlocks` selection correctly left
    untouched. RESIDUALS: (a) the outer-edge T-matrix reach check
    (`create_dofwise_t_matrices_master` must produce weights for the hanging connector
    outer edges) moves to the smoke-run gate; (b) RESOLVED 2026-08-08: the two
    connector cases briefly added to the SIDESET selection switch were removed by
    Christian (connectors are blocks, not sidesets; the cases were dead code plus
    a latent `collect_sideset_dofs` trap — flagged by the wiring jury, W1/P6).
    (Original correction note 2026-08-08, kept for the record: the `tBlocks`
    domain-type switch feeds ONLY the metis
    the `tBlocks` domain-type switch, `cl_MaxwellFactory.cpp:523-539`, feeds ONLY the metis
    partition graph via `KernelParameters::select_blocks` — do NOT add connector cases there.
    Layer blocks are explicitly reset from the partition bitset (`cl_FEM_Kernel.cpp:226-232`)
    and connectors must stay out of metis the same way; their owners come from the D6 pass.
    The fem blocks are built from the IWG's own selection, which already contains ALL mesh
    blocks via `set_block_types_in_magnetic_equation`, `cl_MaxwellFactory.cpp:1773-1786`.)
    The two required wiring spots: (a) `FieldList::collect_block_dofs`
    (`cl_Maxwell_FieldList.cpp:253-292`) — connector types currently hit
    `default: BELFEM_ERROR "Invalid block type"` inside `IWG_Maxwell::initialize`, the first
    abort on the path; add the connector case with an edge-h dof table like conductors
    (shared inner edges resolve to the SAME dofs as the shell edges via the shared Edge
    objects; outer edges are hanging, so their T-matrix weights come from
    `DofData::create_dofwise_t_matrices_master` — verify connector outer edges are reached
    in that pass). (b) the fem-group domain-type update switch
    (`cl_MaxwellFactory.cpp:577-598`) — connector blocks currently fall through to
    `default:` → `GroupActivationMode::Inactive`; add both connector cases.
  - [x] **P3 — dispatch.** DONE 2026-08-08: `case LeftCoating/RightCoating →
    mFunMKF = & maxwell::h_side_connector` added to `IWG_Maxwell::link_to_group`,
    mirroring the Ghost case (compile gate = next build). The calculator scratch
    (`t`, `hb`, `B`, `Tseam`) was already allocated (`cl_IWG_Maxwell.cpp:477-484`).
    Decision recorded: the term lives in IWG_Maxwell (no dedicated IWG_MaxwellWall).
  - [x] **P4 — material assignment (review F10).** DONE 2026-08-08 under the
    edge-coating policy (jury round `review_edge_coating_policy`, Christian's
    decision): there is NO separate coating-material input key — the wall inherits
    the outer stabilizer layer's material (same plating step). Landed: derivation
    on the create AND reload paths (`assign_thin_shell` loops now cover
    `side_connector_blocks()`, label = `materials().first()`); the
    `assign_materials` domain-switch hole closed (Left/RightCoating cases push
    the label and gate on `MaterialType::PureMetal` — YBCO subclasses Metal but is
    typed HTS, so the class test would pass it — plus a warn-if-not-copper);
    the kernel's material null-check upgraded BELFEM_ASSERT → BELFEM_ERROR
    (release-active, the deref below is unconditional). OPEN (deferred, not P4):
    the r′ calibration dial (effective ρ vs an optional interface-resistance
    term, doc §4.3) — decided at R7 calibration time.
  - [x] **P5 — input keys (policy decided 2026-08-08, wired 2026-08-09).** The flag
    is `edge coating : on` — OPT-IN for the release (unanimous jury + Christian;
    default-ON reconsidered only after R7 passes AND B7 gets its cache stamp).
    Guard triad when the flag is on (implemented as loud `BELFEM_ERROR`s in the
    live `create_side_connectors`, 2026-08-08): ≥ 3 layers, identical top/bottom
    layer material (label equality covers RRR — RRR lives on the materials-section
    entry), identical top/bottom thickness (relative tolerance). Width is DERIVED
    from the outer layer (`mConnectorWidth` default is now 0.0 = derive; the old
    hardcoded 1e-4 = 100 µm vs the ~20 µm plated wall is gone). WIRED: the flag
    and optional `edge coating width` key parse per tape in the `thinshell`
    topology section (`read_thin_shell_data`, `cl_MaxwellFactory.cpp:2584-2610`;
    fatal on 2-D meshes, on width-without-flag, and on non-positive width), ride
    the `Protoshell` (the config carrier), and set
    `mCreateSideConnectors`/`mConnectorWidth` per tape at the top of the 3-D
    branch in `ThinShellFactory::create()` (`cl_ThinShellFactory.cpp:305-306`).
    `doc/input_file_reference.md` §8 updated same turn (Codex prose pass applied:
    bool-token list corrected — only `true`/`on`/`yes`/`1` are true, anything
    else silently reads false — and the WIP stop is now stated).
  - [x] **P5b — SideLayerOld purge (2026-08-09).** The entire legacy side-element
    scaffolding was removed (~1470 lines): the dead `create_side_connectors`
    overload, its nine helpers (`create_outer_nodes`,
    `collect_inner_nodes_and_edges`, `create_tangential_edges`,
    `create_binomial_edges`, `create_extra_sideset`, `create_side_elements`,
    `preprocess_binomial_edges`, `create_facet_table`, `map_extra_nodes`,
    `check_inside`), the uncalled `compute_side_indices`, the `SideLayerOld`/
    `SideFacet` structs, dead-only members (`mSideIndices`, `mIndices`,
    `mYbcoLayer`, `mNodeMap`, `mElementBitset`, `mElementMapper`, `mX`, `mXi`),
    the commented old-signature call blocks in `create()`, the commented
    node-source block in `create_buffers`, and orphaned includes
    (`cl_IF_ElementMapper.hpp`, `fn_posv.hpp`, `fn_trans.hpp`, `cl_Queue.hpp`,
    `fn_sum.hpp`). Kept deliberately: the ctor's Air/Ferro element flag-1 loop
    (mutates global mesh state that CutFactory/EdgeCutter may read downstream)
    and `mConnectorsForAllLayers` (guarded in the live builder). Verified: no
    dead names remain, all live functions present, brace balance intact.
    Compile gate = Christian's next build.
  - [x] **P6 — WIP stop removed (2026-08-09)** after the full batch (P1–P5,
    P7a kernel migration) compiled — the `edge coating : on` path now runs
    through to assembly. First smoke run: delete any cached `.bfm` first (B7 —
    the checksum ignores connector settings).
  - [x] **P7a — MaxwellData dispatch migration (2026-08-09, jury round
    `review_sideconnector_maxwelldata`, all four P0s fixed pre-implementation).**
    `h_side_connector` now runs entirely through MaxwellData: new ctor branch
    (`cl_FEM_Calculator.cpp:~166-215`) with PureMetal + constant-mu gates
    (Christian's decision: copper walls, no dMdx port), full pointer set incl.
    `mFundMudH`, metal ρ family reused verbatim; new field variants
    `compute_h_side_connector` (ht + hb + hn, per-element frame via
    `prepare_side_connector_frame`, lazy after `reset()`),
    `compute_T_side_connector` (seam interpolation, `gTbulk` fallback,
    `compute_T_fem` clamp contract); `side_connector_edge_function` moved to
    `cl_FEM_Calculator.cpp` keeping the EF_PENTA6TS friendship;
    workspaces relabeled `normal`/`tangent`/`binomial` (+`Tseam`), the `"b"`/
    `"hb"` aliasing hacks are gone. Newton part = `add_rho_field_tangent` in
    the new `h_side_connector_newton` (dρ/d‖B‖ + dρ/dβ; dρ/dj ≡ 0 for
    metals — the old broken `drho` block is deleted); Picard/Newton dispatch
    in `link_to_group` mirrors conductors. Thermal-coupling fixes: ctor nulls
    the thermal peer via `block_exists` (never adopts the EmptyBlock
    calculator), `allocate()` gets a connector exception, link dispatcher
    routes connectors to `link_element_maxwell` even on coupled runs.
    `save_resistivity` added like the conductor kernels. Compile gate =
    Christian's next build.
  - [x] **P7b — kernel residuals:** revisit the per-element T-copy, now in
    `prepare_side_connector_frame` (visualization write into the mesh "T"
    field from assembly = MPI hazard; move to a postproc/setup pass or gate
    on serial). Master-calc relink documented safe for sequential assembly
    (jury R2, `cl_FEM_DofManager.cpp` block loops).
    **FIXED 2026-08-09 (committed since — working tree clean at `bc578b5e`, verified
    2026-08-10; run gate still pending):**
    viz copy moved to master-side `MaxwellPostprocessor::copy_seam_temperatures()`
    (SideConnector case of `run()`); the assembly-side Tseam read restored to the
    station-ordered `original()` route; `compute_hanging_dofs` restored into the
    vector solve path (was dead since 22432b05).
    **CONFIRMED LIVE 2026-08-09 (probe tier, 9-rank run, e-s.00004):** the
    save-time field gather collects dof-carrying entities only, so wall-node T
    written on ranks != 0 never reaches the file — connector_27 (zero
    master-owned nodes) is 100% frozen at 77.0 while its coincident tape
    partners read 77.013; connector_29's owner-0 nodes are 100% copied, all
    other owners frozen. TWO stacked fixes needed: (1) move the viz copy to the
    master-side SideConnector postprocessor; (2) restore the station-ordered
    read `T( node(tTapeFace[k])->original()->index() )` — the live facet-order
    read is permuted (3/3 proof: PENTA face 2 = {0,3,5,2} up-leg first, HEX8TB
    face 2 = {2,3,7,6} != station {3,2,6,7}, `compute_orientation()` unused).
    Full record: `tmp/ai_exchange/seam_temperature_copy.md` (swept; distilled in
    `devlog/dl20260809_sideconnector_viz_regression_round.md`). Related new
    find: `compute_hanging_dofs()` dead in the vector solve path since 2025-03
    (`cl_FEM_DofMgr_SolverData.cpp:2299` behind the aborting `default:`) — walls
    are the first heavy hanging-edge consumers, consequence trace pending.
  - [ ] **P8 — doc sync.** `side_coating_wall_element.md` §5 still describes the pre-scaffold recovery
    route (per-gap averages); the implemented route is: master layer-block
    calculator + `get_normal_calculator` hop → tape-level `compute_hn` for h_n, penta EF
    evaluated at the mapped lateral-face point for h_b. Update §5/§6 once the plumbing
    stabilizes (docs state facts inline, no defect IDs).
  - [ ] **P9 — parallel deferreds** (out of scope for iteration 1, keep visible; sharpened
    2026-08-08 by a read-only trace): element/facet shipping is already generic —
    `Distributor::select_entities` selects by owner over the global containers, so HEX8TB
    elements and recovery facets travel once D6 gives them owners, and the source-closure
    fixpoint brings the hanging outer nodes/edges along. The real gaps: (a)
    `send/receive_thinshell_data` (`cl_Mesh_Distributor.cpp:1681-1767`) carries only layer
    blocks + thicknesses + materials — no connector record (the parallel counterpart to the `.bfm`
    `coatings`/`widths`/`seams` datasets from B1/B2/B4), so non-root ranks rebuild
    ThinShells without `side_connector_blocks()`; ~~(b) `send_block_data` does not ship block
    thickness, so `EF_HEX8TB::link()` sees NaN width off-root~~ (b) FIXED 2026-08-09
    after the first parallel abort (EF_HEX8TB width assert): `proto::GroupData`
    gained `mThickness`, `send/receive_block_data` broadcast it, and the
    ProtoMesh block materialization applies it — connector walls AND thin-shell
    layers now carry thickness off-root independently of the shell record;
    (c) the T-copy write in
    parallel (see P7b); (d) aura shipping of the recovery-facet master to neighbor ranks.
    Connector ownership itself is already handled (inherits from recovery-facet master, D6).

  - [x] **B8 — reload hang normalization (2026-08-09, Christian's option-1
    ruling).** The viz-decoupled wall inner edges hang 1:1 on their tape sheet
    edge; on `.bfm` reload, `reconstruct_edge_connectivity`'s last-wins tie
    resolution hands ALL tied element slots to the twin sheet, orphaning the
    tape primary (no element flags it → zero dofs → the T-matrix builder
    aborted with "invalid number of dofs on source edge : 0"). Fix:
    `ProtoMesh::normalize_edge_hangs()` re-points every 1:1 edge-on-edge hang
    to the tie representative (same key/last-wins resolver as the slot
    reconstruction; unique keys self-resolve, weight transfers unchanged —
    layer twins are same-oriented by construction). ORDERING IS LOAD-BEARING:
    it must run AFTER `load_hanging_entities()` — the first attempt lived at
    the end of `reconstruct_edge_connectivity`, which `load_edge_data()` calls
    BEFORE the hangs exist, and was a silent no-op. Called from
    `BfmFile::load()`. KNOWN ASYMMETRY, accepted with
    the ruling: at tied (non-cut) stations a create-run carries two sheet
    traces, a reload-run one — this matches what the slot reconstruction
    already does to the tape elements themselves (pre-existing, first exposed
    by the wall hang). If through-thickness branching must survive reload at
    tied stations, the tie-break needs sheet awareness instead (option 2,
    bigger surgery, fixes the latent tape infidelity too).

  Legacy consumer wiring (`connect_side_connectors_with_facets`, `h_penalty`/`h_tb`)
  stays NOT recovered — clean start confirmed. (after: R2, R3; P-items after the
  review-round fixes of 2026-08-05)
- [x] **R6 — Calculator / Postprocessor (completed 2026-08-09)**: the Nedelec data
  path case in `cl_FEM_Calculator.cpp` landed 2026-08-08 together with R4 (HEX8TB
  joined the `nedelec_data_linear_h` group); the postprocessor side landed
  2026-08-09: HEX8TB added to both thin-shell type lists in
  `cl_FEM_Postprocessor.cpp` (side-local aura-skip — wall recovery reads master
  state, like the layers), `Left/RightCoating` accepted by the metal branch of
  `MaxwellPostprocessor::select_blocks_and_materials` (PureMetal gate guarantees
  no-jc), and `MaxwellFactory::create_postprocessors` now collects
  `side_connector_blocks()` (rank-0 collect + broadcast of ids AND true domain
  types — the types must travel because off-root shells lack the connector
  record, P9a) into a SEPARATE ThinShellConductor instance (s==2 pass with its
  own block map): a Postprocessor instance is single-element-type by
  construction (`compute_node_matrices` sizes the patch matrices once from
  `mElementType` and asserts homogeneity), so PENTA6TS layers and HEX8TB walls
  must never share one — the first merged attempt died in
  `compute_element_coeffs` on exactly that (Blaze size mismatch, caught on
  Christian's first postproc-enabled run). Seam nodes are recovered by BOTH
  instances (side-local claims; wall instance runs second and wins the write) —
  acceptable, both are legitimate side-local estimates. UPDATE 2026-08-09
  (Christian's ParaView finding: E·q misses the normal/binomial components):
  the walls now have their own `MaxwellPostprocessorType::SideConnector` —
  `compute_side_connector` recovers the FULL field via
  `MaxwellData::compute_h` (same ht + hb + hn route as the assembly kernel),
  the type owns its block-selection case (Left/RightCoating, no jc filter —
  PureMetal gate) and the s==2 factory pass creates it; the temporary
  Left/RightCoating acceptance in the ThinShellConductor case is reverted.
  This supersedes the tangential-only caveat below for H and B. Note: the
  postproc H on walls is the wall's own `E·q` (tangential part only, without the
  hb/hn recovery) — same simplification the thin-shell layers already accept;
  J = `C·q` is the physically meaningful wall output. (after: R3)
- [ ] **R7 — Verification gate**:
  - [ ] numeric orientation check on both connector signs (resolves O1);
  - [ ] single-tape open-curve case: wall current `(h_l − h_{l+1})/L` vs analytic r′;
  - [ ] R0.3 probe: interior + buffer-adjacent side-trace values vs g on a corc multilayer;
  - [ ] twisted-helix periodic regression (cross-block partners);
  - [ ] cut-terminates-on-side-curve station test (perpendicular-cut trap detection).

## .bfm persistence tracker (3-AI gap analysis 2026-08-04, all items 3/3 agreement)

The reload path never rebuilds connectors (`MaxwellFactory::create_thinshells` early-returns
on a non-empty `thin_shells()`), so anything the file does not carry is absent for the whole
run. What DOES round-trip: HEX8TB elements + blocks (incl. domain type), recovery facets
(master/slave/index/orientation), facet ids (so `facet id = wall id + 1` holds), hidden
sideset flags, hanging outer nodes/edges with sources+weights, HEX8TB edge slots. Facets are
saved without topology, so the loader re-derives their nodes from the master — which is
exactly what `set_master(..., true)` does at construction.

- [x] **B1 (P0):** RESOLVED 2026-08-04 (Christian + Claude): per-shell connector record in
  `save_thinshell_data`/`load_thinshell_data` — `coatings` (connector block ids) alongside
  the renamed `layers`; the loader repopulates `side_connector_blocks()`, so the ownership
  pass in `Kernel::partition_mesh` sees them again. Audit fix: the save loop indexed the
  container with its own uninitialized output buffer (`( tIDs[ k ] )` → `( k )`).
- [x] **B2 (P0):** RESOLVED 2026-08-04: new `widths` dataset carries the connector block
  thickness, and the loader stamps it back onto the block — this is now the ONLY record of
  it, since `Block::mThickness` is not part of the group data and
  `ThinShell::set_thicknesses` only reaches layer blocks (`cl_ThinShell.cpp:39-52`).
  `EF_HEX8TB::link()`'s width lookup no longer sees NaN after a reload.
- [x] **B3 (P1):** RESOLVED 2026-08-04: re-derived at load, no new dataset. The facet's two
  dof edges ARE the wall's inner-face edge slots, and the slave face index is persisted —
  `Element_HEX8TB::get_edges_of_facet` now supports the lateral faces 0/2 (bottom edge
  first, then top, matching the facet slot convention), and `load_thinshell_data`'s
  coating loop rebuilds each recovery facet's edge container from
  `slave()->get_edges_of_facet( index_on_slave() )`. Runs after
  `reconstruct_edge_connectivity`, so the wall slots exist. (Also softens D2: only faces
  1/3 error out now — they carry no longitudinal dof edges.)
- [x] **B4 (P1):** RESOLVED 2026-08-04: `create_side_connectors` now pushes
  `side_connector_sidesets()` next to the block push (one block + one sideset per curve,
  index-aligned — the bfm `seams` dataset pairs them by position), and the loader restores
  both. This fixed the FRESH path as well, where the container had always been empty.
  **B6 is now unblocked.**
- [x] **B5 (P2, latent):** RESOLVED 2026-08-04 without persistence: the wall's sheet
  convention is SLOT-BASED and deterministic — bottom slots 0/1 take the duplicate sheet
  when it exists, top slots 2/3 always the primary (`create_side_connectors`, tA-tD
  selection). `reconstruct_edge_connectivity` now applies the same rule: HEX8TB slots 2/3
  resolve via a new find-FIRST lookup (first = primary, since primaries are appended to
  the mesh before duplicates and the sort is stable), everything else keeps last-wins.
  No behavior change when no duplicates exist (first == last). Layer blocks' own
  last-wins convention untouched (pre-existing, accepted).
- [x] **B6 (P2):** RESOLVED 2026-08-04 by Christian's design call: instead of teaching the
  topology snapshot about connectors, the factory now refuses to CACHE enriched meshes —
  `BELFEM_ERROR( ! mUseEnrichment, ... )` before `mMesh->save()` in `MaxwellFactory`. Cut
  enrichment is currently dead anyway (the implemented bubble functions are not the right
  space; something different is needed per the literature). The remaining type-map leak of
  connector entities on reload is benign: they land in the `Left/RightCoating` buckets,
  and every consumer switches on specific domain types — `Topology::select_blocks` sends
  them through `default: pass` (verified), so they never enter the phi/non-phi selections.
- [x] **B7 (P2, deprioritized by Christian 2026-08-04)** — ~~nothing invalidates a stale
  `.bfm` when connector settings change — the checksum hashes only the original gmsh
  geometry, so flipping the connector flag or changing the width silently reuses the
  cached mesh.~~ **The processing-options stamp half is DONE (2026-08-11, DR-21).**
  `meta/config` + `meta/config_text` now carry a fingerprint of the mesh-defining settings,
  and cache reuse requires it to match alongside the geometric checksum. The tag walks the
  whole `topology` subtree, so **the connector flag and `edge coating width` are both in it**
  — verified on the `sidecoating` deck, whose canonical text contains
  `edge coating = true`. A mismatch rebuilds and names the changed line; a `.bfm` written
  before the stamp existed rebuilds once with its own message. **So the "delete any cached
  mesh before the first connector-enabled run" advice elsewhere in this file is now
  obsolete.**
  **Not done — the second half of the original ask:** a real `/meta/format_version` gate.
  The stamp answers "were these settings used?", not "can this reader understand this
  file?". The `layers` rename break recorded in `src/mesh/doc/bfm_file_format.md` §6 still
  has no reader gate. That is a separate, smaller item.

B7's stale-settings half is closed; its format-version half is the only open item of this set. Validation for the whole set: a save→reload cycle
with the WIP factory error temporarily gated. Analysis + audit record:
`tmp/ai_exchange/sideconnector_bfm_gaps.md`.

## Defect tracker (carried over from the factory audits)

- [x] **D8 (HIGH, fusing path — THE verified blob generator, 2026-08-06):** at 61 of 151
  rim stations of the Garber CORC (fused run, `BELFEM_PROBE_FUSED_ROWS` v2 dump,
  210,066 rows), the composed edge-on-node constraints of DIFFERENT interface levels of
  the SAME station disagree by the free-cut λ: some levels read `h_e = Δφ`, others
  `h_e = Δφ ± λ` (dof 179126), some hang on different node representatives altogether
  (original vs cut duplicate, e.g. 3533 / 35356 / 157464). Adjacent levels' tangential
  traces then differ by λ (≈ amps, time-oscillating) → inter-level jump → element-scale
  J = curl h blobs strung along the rim; constraint-level, hence algorithm- and
  solver-independent (persists under Picard+MUMPS). Individually every row is
  well-formed: sign coherence 141,901/141,901 clean, zero internally inconsistent
  hang-pair groups, zero-weight λ entries are correct cancellations (+λ−λ across an
  edge on one side). The defect is CROSS-level branch mixing: `connect_side_edges`
  sources the ORIGINAL mid-surface curve nodes (pure Δφ) while the outer-interface
  anchors (`hang_thinshell_edges_on_nodes_bottom/top`) resolve through sheet/cut
  duplicate nodes (Δφ ± λ) at cut-adjacent stations. Fix direction: all levels of one
  station must reference the SAME branch of g. Related: the 4b seam diagnostic's
  sign-coherence defect at in-plane cuts 0/3 (CutData). Minor loose ends from the same
  dump: 6 of 156 edges per interior block absent (flagged or unhanged — identify), and
  425 rows with weight-sums ±2/±3 (multi-λ content, census pending).
  **FIX IMPLEMENTED 2026-08-06 (committed since; closure below):** cut-aware single
  authority landed in `cl_ThinShellFactory` (`compute_side_authority` + rewired
  `connect_side_edges( Layer* )` / `connect_side_nodes`) — every rim station now
  resolves against the cut-composed air volume nodes the anchors read. The 6-missing
  edges (→ D9) and 425 multi-λ rows (→ D10) are classified and tracked in
  `todo/deferred/side_edge_fusing_cut_aware_plan.md` (§4a, §5). **Gate criterion REVISED
  same evening:** first post-fix run kills the oscillating blobs (incoherence
  fixed) but reveals that 90/348 rim stacks have genuinely different λ branches
  above vs below the tape (Christian's insight) — residual cross-level "mixing"
  at exactly those stacks is now the reopened-O1 physics signal, not a bug.
  **CLOSED 2026-08-08 (Christian's disposition):** the defect this tracker names —
  the cross-level incoherence — is fixed (committed) and verified end-to-end
  (oscillating blobs gone on the Garber run); the residual two-valuedness at the
  90 stacks is O1 physics, owned by the fusing plan. The confirmatory station
  census and physics gates live on ONLY as `mFuseEdges` re-enable preconditions
  in `todo/deferred/side_edge_fusing_cut_aware_plan.md` (its R4/R5 — fusing-plan step IDs,
  not this plan's R4). D8 cannot affect this campaign: the connector and fuse
  paths are mutually exclusive branches in `ThinShellFactory::create()`, and
  `mFuseEdges` is off in source (`cl_ThinShellFactory.hpp:280`).

- [x] **D7 (HIGH, fusing path):** RESOLVED 2026-08-06 (Fable session, on Christian's
  go-ahead): `connect_side_nodes` now indexes with `original()->index()`, matching
  every other `Layer::Nodes` consumer (`:1521-1524` etc.) and the D1 pattern.
  Verified worse than "wrong node": `create_node_container` resets ALL master
  nodes to `gNoIndex` and re-indexes only the sorted original surface nodes, so a
  cut-duplicate station indexed the container with the `gNoIndex` sentinel —
  out-of-bounds, not merely mis-aimed. Unblocks re-enabling `mFuseEdges` for the
  no-cut experiment (gate 2 of the fusing campaign); the §6 orientation gate and
  the `reset_source_container()` clobbering audit remain open before any
  cut/periodic fused run is trusted.
  Original finding: `connect_side_nodes` indexes the layer node array with
  `aTargetNodes( tOrg->index() )` — NO `original()` normalization
  (`cl_ThinShellFactory.cpp:3007-3009`) — the exact sibling of the bug D1 fixed in
  `compute_side_edge_indices` on 2026-08-03. On cut-carrying side curves (cuts terminate
  on side curves; Gregory's periodic CORC) the node hang can grab the wrong layer node.
  Found by Grok 2026-08-06 (physics jury, `tmp/ai_exchange/review_side_edge_fusing_physics.md`),
  citation-verified by Claude. Top verified suspect for the gap-edge J oscillations that
  led to the `mFuseEdges` off switch; latent while `mFuseEdges = false`, MUST be fixed
  before the fuse is ever re-enabled. Related jury outcomes: the edge-on-node hang IS the
  g anchor (DofData LINE2 +1/−1, `cl_FEM_DofMgr_DofData.cpp:3560-3645` — "missing anchor"
  refuted); coverage l = 1..N−2 coherent (outer interfaces φ-anchored via
  `hang_thinshell_edges_on_nodes_bottom/top`); doc §6 orientation gate still unrun.
  Discriminating probe agreed by both auditors (extends R7/R0.3): dump per-station
  hanging-source rows (target edge id, direction flag, source node ids
  original-normalized, weights, h − g residual) on a flat 3-layer NO-CUT tape with
  `mFuseEdges` on, plus one cut/periodic station on the corc mesh.

- [x] **D6 (MEDIUM, parallel):** RESOLVED 2026-08-03 (rule by Christian): side connector
  elements inherit ownership from the MASTER of their recovery facet (= the layer block
  element the wall spans), and the recovery facet follows the same owner. Implemented in
  `Kernel::partition_mesh` right after the layer-block ownership pass (masters have final
  owners there), via `mesh->facet( id+1 )->master()->owner()`; debug owner asserts extended
  to connector blocks. Recovery facets are NOT in the metis facet graph (built from
  `tShell->facets()` only), so the min-consistency sweep cannot interfere.
- [x] **D1 (HIGH):** RESOLVED 2026-08-03: `compute_side_edge_indices` now keys with
  `original()->index()` on both ends of the stride walk, matching the original-keyed
  edge map — cut-duplicate stations resolve correctly. (= audit R1c, Claude 2026-07-29,
  confirmed Codex + Grok; fix approved by Christian.)
- [ ] **D2 (LOW):** `Element_HEX8TB::get_edges_of_facet` only supports facets 4/5; side
  facets error out. This is fine for the wall element's purpose — document or extend when facet-based
  BCs touch connector blocks.
- [x] **D3 (LOW):** RESOLVED by lifecycle (found 2026-08-03): `unfinalize()` clears
  `mElements` (`cl_Mesh.cpp:896`), and the re-finalize after the ThinShellFactory calls
  `collect_elements_from_blocks()` over ALL blocks — connector blocks included — so
  HEX8TB elements DO enter `mMesh->elements()` and the element maps.
- [x] **D4 (MEDIUM, R5 gate):** RESOLVED 2026-08-03 (Christian's design): the factory no
  longer writes `physical_tag` on HEX8TB elements — the material machinery owns the tag.
  `EF_HEX8TB::link()` recovers the layer thickness through the recovery facet instead:
  `mesh->facet( element->id() + 1 )->master()->block_id()` (facet id = wall id + 1 by
  construction; the master IS the layer block the wall spans). Lifecycle verified:
  `unfinalize()` clears `mFacets`/`mElements` (`cl_Mesh.cpp:896-897`), the re-finalize
  after the ThinShellFactory re-collects from ALL sidesets and rebuilds the facet map.
- [x] **D5 (MEDIUM, semantics):** RESOLVED 2026-08-03, by the same design: the wall spans
  exactly one shell layer block (sheet j to sheet j+1), and the recovery facet's master
  is that block — d = master block thickness, unambiguous. No layer-vs-mean question.

## Open questions

- **O1 — Element handedness (audit C3):** Grok leans toward the `tSign==1` winding being inverted
  vs. the standard HEX8 Jacobian; Codex leans toward it being fine; node numbering is deliberate
  (Christian, 2026-07-29 — HEX8 convention, edge scheme differs). The MeshChecker guard
  converts a wrong sign into a loud error. PARTIALLY DEFUSED 2026-07-31: the implemented
  EF builds its own right-handed orthonormal frame (b = n×t) and detJ = Lwd/8 > 0 by
  construction, so the operator math never sees the node-coordinate Jacobian sign. The
  numeric R7 check remains to verify the factory geometry is not actually inverted
  (MeshChecker would fire) and that the frame's b agrees in sign with the factory binormal
  where recovery projections will use it.
- **O2 — Coated variant:** one new internal edge dof `h_in` per station, coupled only through the
  wall r′; enumeration of the four virtual dofs per the concept slides. Requires the
  atomicity guard: a coated station must never ship without its wall element.
- **O3 — Outer edge sheet selection (audit N2):** the connector bottom prefers
  `OuterEdgeDuplicates` when present, while the legacy path always took primary outer
  edges. Confirm intent for material-interface levels (probably correct dual-sheet symmetry;
  decide and document).
- **O4 — Closed-loop tapes / threefold demo:** needs (a) facet-based sign anchor in
  `compute_binomial_vectors` (also removes the terminal-topology requirement for open
  curves), (b) periodic tangent stencil at the seam, (c) the perpendicular-cut treatment
  (duplicated outer nodes at the crossing station) — the demo forces the cut-crossing
  design rather than guarding it away.
- ~~**O5 — Edge direction machinery:**~~ RESOLVED 2026-08-03 (analysis, 3-AI check
  pending numeric confirmation): the standard machinery already covers HEX8TB —
  `fem::Element::link_dofs` (`cl_FEM_Element.cpp:180-183`) and the aura ctor (:44-48)
  call `compute_edge_directions()` whenever the mesh element has edges. For HEX8TB the
  exact-match branch always resolves, because the factory inserts edges whose node
  objects ARE the element's node-table entries in both `tIsReversed` branches: forward
  → edge intrinsic n0→n1 == local +ξ → s=+1; reversed → s=−1 — exactly the
  global-dof-to-local-ξ conversion the EF needs. No `set_edge_direction` wiring
  required. Guarded by the MeshChecker no-swap policy (a node swap would break the
  slot-tie this relies on).
- **O6 — Curve closure has NO reliable oracle (found 2026-07-30, tape_hphiTrun; sharpened
  by Christian):** two independent problems. (1) `Curve::is_closed()` is semantically
  overloaded: `CurveFactory::sort_end_nodes` sets it geometrically
  (`cl_CurveFactory.cpp:1070`), but the CutFactory overwrites it on side curves to mean
  "not shared with another shell" (`cl_CutFactory.cpp:1591`, feeding the duplication
  exclusion at `:1608`) — true on every single-shell model regardless of chain topology.
  (2) The node-list convention for genuinely closed loops differs per pipeline:
  `CurveFactory::collect_nodes` repeats the start node at the end
  (`cl_CurveFactory.cpp:1392–1428`), while `close_terminal_loops` stores one node per
  segment with NO repeat (`cl_CutFactory.cpp:2021–2027`) — so `first == last` is not a
  closure test either (confirmed in the tape_hphiTrun bfm). Consequences: the
  `compute_side_edge_indices` walk covers the closing segment only under the repeat
  convention; `compute_binomial_vectors` now relies solely on the terminal-anchor search
  (loud error, convention-independent) instead of any closure guard. Before closed side
  curves become a real use case (threefold demo, O4), the side-curve pipeline must pick
  and enforce ONE node-list convention at construction, and the flag should be split or
  renamed (e.g. `is_shared`). Side notes from the same bfm inspection: two of the four
  terminal rings end on a duplicate node instead of closing (`closed=1`, chains 15→4778 /
  12→4777); and the LINE3 case in `close_terminal_loops` falls through into the
  "invalid element type" error (missing `break`, `cl_CutFactory.cpp:2030–2045`) —
  latent until order-2 terminals exist.

## Audit trail

- `tmp/ai_exchange/sideconnector_wip_audit.md` — rounds 1–3 (SideLayer, refactor, HEX8TB
  creation; ephemeral, distilled here and into the devlogs).
- `devlog/dl20260728_side_connector_revival.md`, `devlog/dl20260729_side_connector_fixes_hex8tb.md`.
- Physics consensus record: fusing theory + trap catalog (2026-07-28, Claude + Codex + Grok).
