# 2D h-φ Thin-Shell Pipeline — Gap Analysis vs. the 3D Reference Path

**Date:** 2026-07-01
**Purpose:** Procedure-by-procedure comparison of the 2D h-φ thin-shell pipeline against the
3D reference path; inventory of what is MISSING, STUBBED, or WRONG in 2D.
**Method:** Read-only audit. Five parallel stage audits (one per pipeline stage), findings
cross-checked; the highest-impact claims were independently re-verified in a second pass
(marked "re-verified" below). Companion checklist: `todo/2d_thinshell_todo.md`.
> **ARCHIVED SNAPSHOT — not a tracker, not maintained (marked 2026-08-11).**
> This is a read-only audit taken on 2026-07-01. Its *reasoning* — the
> procedure-by-procedure comparison against the 3D reference path — is what it is for, and
> that does not expire. Its *currency* does: several gaps it analyses are fixed in code
> (the EF_QUAD4TS basis/curl operator in `4a42d982`, the top-node tie in `43474c9f`, the
> LINE2 case of `to_master_orientation`, the LINE3 fall-through), and every `file:line` in
> it is baselined to 2026-07-01 and has drifted. **The live artifact is
> `todo/2d_thinshell_todo.md`;** check any finding here against that checklist and the tree
> before acting on it.

**Confidence convention:** claims are code-confirmed unless tagged `[INFERENCE]`
(with what would confirm them). Classifications: MISSING / STUB / PRESENT-BUT-WRONG /
PRESENT-BUT-UNTESTED / OK / DEAD.

---

## Currentness note (2026-08-09)

This file is an **analysis snapshot**, not a task tracker — the live checklist is
`todo/2d_thinshell_todo.md`, which carries the 2026-08-04 execution update. Two of the
gaps analysed below have since been closed in code and their descriptions here should be
read as historical: the `EF_QUAD4TS` basis/curl-operator gap (fixed in `4a42d982`,
root-cause record `devlog/dl20260804_2d_thinshell_layer_alternation.md`) and the flipped
top-surface node tie (fixed in `43474c9f`). Every `file:line` in this document refers to
HEAD of 2026-07-01 (or 2026-07-10 in the addendum) and has drifted; re-locate by symbol.

## Second-iteration addendum (2026-07-10)

The analysis below is the 2026-07-01 snapshot; its line numbers refer to that HEAD. A
second-iteration verification pass (Claude subagents, Grok, and Codex; distilled in
`devlog/dl20260710_2d_thinshell_plan_v2.md`) resolved the open inferences. Read these
corrections first:

- **Gap G3 (`compute_element_adjacencies` poisoning) is FIXED on HEAD** (2026-07-02):
  the 2D branch now resets element containers and calls
  `Mesh::reset_connectivity(ElementToElement)` (cl_CutFactory.cpp:2761-2773).
- **E3 answered:** a single-layer 2D shell creates NO Ghost sideset (`nb==1` early
  return, cl_ThinShellFactory.cpp:1837) — the `h_ghost` blockers (E1/E2) gate
  multi-layer shells only, so the "minimum path to M1" at the end of this file is
  overstated: drop `h_ghost`/scratch from it.
- **D2 top-node pairing: CONFIRMED cross-wired** (layer insert order n2=top(F1),
  n3=top(F0), cl_ThinShellFactory.cpp:1218-1221, versus `get_top_nodes` + orientation flip).
- **A4/G4 orientation sign: CONFIRMED mirrored** (2D `cross(shell, tape)` vs 3D
  `cross(boundary, shell)` — anticommuted operands).
- **B3 terminal halving: LARGELY REFUTED** — the parser doubles 2D terminal lists
  (cl_MaxwellBoundaryConditionFactory.cpp:125, :144-151), so the halving is consistent
  for any parsed deck; only defensive validation remains.
- **E4 sizing mix: CONFIRMED silent wrong-physics in debug AND release** — Armadillo
  expression assignment silently resizes (cl_AR_Vector.hpp:375-381); BELFEM's
  asserting `op_VectorPlus` is bypassed by expression templates.
- **E5 escalated:** no live writer of layer-element edge directions was found (the
  tape-roller path is dead) — the default bitset violates the `EF_QUAD4TS::E()`
  invariant `mS[1] = -mS[0]`. Caveat: 3D PENTA6TS works in production, so either
  a writer was missed or 3D tolerates the defaults — see task T1 in the todo.
- The current work plan, sprint scope, and re-baselined line numbers live in
  `todo/2d_thinshell_todo.md` (v2).

---

## Executive summary

The 2D thin-shell path is **not a stub** — every pipeline stage has real, reached 2D code
(dispatched almost everywhere by runtime `mesh->number_of_dimensions()` branches, plus
`ModelDimensionality::TwoD` and `ElementType` switches). The mesh-side machinery
(shell opening, cohomology, thin cuts, QUAD4TS layer creation) is mostly present and
plausibly correct for **linear, serial, non-periodic** cases. The pipeline is broken by a
small number of concrete defects, the hardest of which sit at the downstream end:

1. **Assembly-killer:** `maxwell::h_ghost` is hard-coded for PENTA6TS (6 edge dofs); on 2D
   ghost sidesets it reads/writes out of bounds in release and asserts in debug
   (mt_maxwell_h.cpp:1766-1937). Every multi-layer 2D shell hits this.
2. **Assembly-suspect:** thin-shell work vectors `bt`/`bn`/`n` are hard-coded size 3 while
   `b`/`bm`/`bs` are sized `mesh::dimension(type)` = 2 in 2D; `h_ts_*` mixes them
   (cl_IWG_Maxwell.cpp:594-596 vs mt_maxwell_h.cpp:137/2271).
3. **Edge-hanging breaks:** conductor-above-shell aborts via `to_master_orientation` on a
   1-edge LINE facet (fn_to_master_orientation.cpp:676); the air-above top-node pairing is
   suspected cross-wired (QUAD4TS facet-CCW top vs PENTA position-aligned top).
4. **Silent-topology trap:** 2D terminal curves are user-provided and unvalidated; if the
   `topology/curves` entry is missing or its first id isn't the tape, the tape endpoints are
   silently duplicated (torn tip) with no diagnostic (cl_MaxwellFactory.cpp:2497-2526,
   cl_CutFactory.cpp:1664-1699).
5. **Order-2 is blocked end-to-end in 2D** (five independent blockers), and **2D + periodic**
   is unsupported at two independent points.
6. **Connectivity poisoning:** `CutFactory::compute_element_adjacencies()` is face-based; in
   2D it stamps `Connectivity::ElementToElement` on an empty adjacency, which the mesh never
   resets and the partitioner later trusts (re-verified; see Stage 2/3).

---

## Stage 1 — Spatial cut creation

### 3D reference sequence

1. `MaxwellFactory::MaxwellFactory` — cl_MaxwellFactory.cpp:66-70 — `read_mesh()` (sets
   `mComputeCohomologies`; fresh path refinalizes with `Connectivity::Compute` :304-306),
   constructs `mesh::Topology` (rank 0).
2. `MaxwellFactory::read_domain_types` — cl_MaxwellFactory.cpp:107 (def :313-342) — sets
   block/sideset `DomainType` from the input `topology` section.
3. `MaxwellFactory::create_curves` — call :108-114, def :347-374 — 3D:
   `CurveFactory::intersect(a, b, id)` (cl_CurveFactory.cpp:561-615) builds terminal curves
   at sideset∩sideset intersections.
4. `create_magnetic_kernel()` → `create_cuts()` — cl_MaxwellFactory.cpp:432 → :744-793.
5. `read_thin_shell_data` — :749, def :2394-2547 — builds `Protoshell`s
   (cl_Protoshell.hpp:40) from `topology{thinshell{...}}` + `layers`; associates terminal
   curves (3D: matches `sideset_a` OR `sideset_b` and swaps so `sideset_b` = tape,
   :2510-2523); tags tape sidesets `DomainType::ThinShell` (:2539-2545).
6. `Topology::run` — :752 → cl_Topology.cpp:38-47 — `collect_block_and_sideset_types`
   (:164), `detect_sideset_types` (:213), `select_blocks` (:416), `select_sidesets` (:467).
   No dimension branch anywhere in the file.
7. `MaxwellFactory::fix_facet_masters` — :773, def :1066-1134 — orients facets so the master
   is the higher-priority domain type. Dimension-agnostic.
8. `create_cuts_sub_master` — :774, def :798-833 — constructs `mesh::CutFactory` (:802,
   ctor cl_CutFactory.cpp:42-61), `create_thin_shell_cuts()` (:808), harvests
   master/slave node lists (:813-814), `create_terminal_list()` (:818, def :978-1061),
   `set_terminals` (:820), then `tFactory.run()` (:823).
9. `CutFactory::create_thin_shell_cuts` — cl_CutFactory.cpp:1600-1643 — the mesh surgery:
   - :1605 gate on `mTopology->groups(DomainType::ThinShell)`;
   - :1613 `create_curves_for_thinshells()` (:2380-2388) → 3D
     `create_side_curves_for_thinshells_3d()` (:2430-2465) →
     `CurveFactory::thin_shell_side_curves` (cl_CurveFactory.cpp:74-162) fills
     `Protoshell::side_curves()`; then `mMesh->create_curve_map()`;
   - :1615 `duplicate_nodes_on_face_sidesets()` (:1646-1868) — 3D open/closed side-curve
     classification :1700-1757 (closed side-curve nodes NOT duplicated); duplicates created
     :1792-1809; periodic pairs fixed :1812-1841; fills
     `mThinShellMasterNodes`/`mThinShellSlaveNodes` :1845-1865;
   - :1616 `relink_slave_elements_with_duplicate_nodes()` (:1871-1950);
   - :1619 `duplicate_and_relink_facets()` (:1953-2008) — originals parked in
     `mTemporaryThinShellSidesets`, tape facets replaced by two one-sided copies;
   - :1622 `relink_non_thinshell_facets()` (:2010-2046);
   - :1624 `close_terminal_loops()` (:2048-2192);
   - :1627 duplicates appended to `mMesh->nodes()`;
   - :1630-1633 3D `orient_terminal_curves()` (:2194-2252 + `_sub` :2254-2294);
   - :1639 `update_node_indices()`.
10. `CutFactory::run` — cl_CutFactory.cpp:116-192 — refinalize (:121-122), `create_edges`
    (:123-124), **3D-only** `create_faces`/`finalize_faces` (:125-129), periodicity
    (:132-136), `Curve::assign_edges` (:139-142); then handoff to Stage 2/3
    (`compute_cohomologies` :164, `compute_thin_cuts_and_duplicate_interface_nodes` :171,
    `restore_thin_shell_sidesets` :172, `compute_element_adjacencies` :178, refinalize
    :181-184).

### 2D current state

| # | 3D procedure | 2D counterpart | dispatch | classification |
|---|---|---|---|---|
| 1 | `CurveFactory::intersect` (cl_CurveFactory.cpp:561) | `from_2d_sidesets` (:731-822) | runtime if, cl_MaxwellFactory.cpp:357-367 | PRESENT-BUT-UNTESTED |
| 2 | terminal-curve association 3D a/b swap (cl_MaxwellFactory.cpp:2510-2523) | 2D `sideset_a`-only match (:2501-2507) | runtime if | PRESENT-BUT-UNTESTED, unvalidated |
| 3 | `Topology::run` (cl_Topology.cpp:38) | same code | none | OK |
| 4 | `fix_facet_masters` (cl_MaxwellFactory.cpp:1066) | same code | none | PRESENT-BUT-UNTESTED |
| 5 | `create_side_curves_for_thinshells_3d` (cl_CutFactory.cpp:2430) | none — `if (dim==3)` guard :2382 | runtime if | MISSING by design; endpoint role taken by row 6 |
| 6 | 3D closed-side-curve node exclusion (:1700-1757) | 2D endpoint exclusion (:1664-1699) | runtime if :1664 | PRESENT-BUT-UNTESTED; silently no-op if `terminal_curves()` empty |
| 7 | relink/duplicate facet surgery (:1871, :1953, :2010) | same code | none | OK (generic over Facet/Node) |
| 8 | `close_terminal_loops` (:2048-2192) | same code | none | PRESENT-BUT-WRONG for LINE3 (see below) |
| 9 | `orient_terminal_curves` (:2194) | `orient_terminal_curves_2D` (:2296-2377) | runtime if :1630-1637 | PRESENT-BUT-UNTESTED; sign suspect (see below) |
| 10 | `create_faces` in `run()` (:125-129) | correctly skipped | runtime if | OK |
| 11 | `compute_element_adjacencies` (:2913-2943) | none — same faces-only code | none (silent) | PRESENT-BUT-WRONG in 2D (see Stage 2) |
| 12 | `restore_thin_shell_sidesets` (:876-910) | same code | none | OK |

Key details:

- **2D terminal-curve semantics differ from 3D.** In 2D the "terminal curve" is the tape
  polyline itself: `from_2d_sidesets` builds `Curve(aID, sideset(aSideSets(0)), nullptr)`
  (cl_CurveFactory.cpp:755) — `sideset_b` is null. `read_thin_shell_data` then matches only
  `tSideSetID == tCurve->sideset_a()->id()` (cl_MaxwellFactory.cpp:2503). So a
  `topology/curves` entry whose **first** listed sideset is the tape is *mandatory* in 2D,
  and nothing enforces it. If absent: (a) the endpoint exclusion (:1678-1698) unflags
  nothing, so the tape tip nodes are duplicated and the tape is torn open; (b)
  `close_terminal_loops` no-ops (`tSegments.size()==0` continue, cl_CutFactory.cpp:2060);
  (c) `orient_terminal_curves_2D` and the Poisson terminal fixing
  (cl_CutFactory.cpp:249-258) no-op. No error or warning anywhere. **Re-verified: gap G1.**
- **`close_terminal_loops` LINE3 missing `break`** (re-verified): the node-container switch
  case `ElementType::LINE3` (cl_CutFactory.cpp:2171-2181) falls through into `default:` →
  `BELFEM_ERROR(false, "Invalid element type at terminal loop %lu")` (:2182-2186). Every
  LINE3 terminal curve — i.e. any second-order run with terminal curves, 2D or 3D — aborts
  with a misleading message. One-line fix.
- **`orient_terminal_curves_2D` sign convention** `[INFERENCE]`: 2D uses fixed `tS = +z`
  (:2310) and tape normal `tB = (-py, px)` (:2356-2357); the 3D reference feeds
  `sideset_a` (= boundary, after the :2513-2516 swap) into the operand slot commented as
  "shell" (:2228-2232), while 2D feeds `sideset_a` (= tape, :2331). The operand roles are
  mirrored between the paths, so if 3D is the validated convention the 2D curve orientation
  may be flipped. Confirm numerically (known geometry, compare imposed current direction)
  or against `src/fem/maxwell/doc/thin_shell_facet_orientation.md`.
- **Dead spatial-cut machinery** (no call sites repo-wide, grep-verified by two independent
  auditors): `CutFactory::create_sidesets_2d/3d` (:605/:745), `create_cut_sideset_2d/3d`
  (:1127/:1388), `flag_nodes_and_facets_of_tape_sidesets` (:2470 — contains a *different*
  2D endpoint heuristic :2487-2529 than the live one), `create_facet` (:2571),
  `connect_facet_to_slave` (:2623), `check_element_types` (:915), `save_edges` (:946), all
  of `SideSetFactory`, all of `CutProcessorManual`. Live cut-sideset creation is
  `CutProcessor`/`CutData` (Stage 2/3). This dead layer can mislead audits.
- Cosmetic: `tProtoQUAD4TS` at cl_CutFactory.cpp:1715-1718 is a botched rename of
  `tProtoShell2` inside the 3D branch.

### Gap

- **G1** — no validation of 2D terminal-curve association (cl_MaxwellFactory.cpp:2497-2526);
  silent mesh corruption when missing.
- **G2** — LINE3 fall-through abort in `close_terminal_loops` (cl_CutFactory.cpp:2171-2186).
- **G3** — `compute_element_adjacencies` has no 2D branch (details under Stage 2).
- **G4** — 2D curve-orientation sign unverified vs the 3D convention `[INFERENCE]`.
- **G5** — dead duplicated cut machinery obscures the live path (cleanup, not correctness).

---

## Stage 2 — Cohomology computation

### 3D reference sequence

1. `create_cuts()` resolves `mCutAlgorithm` from the `homology` input section
   (cl_MaxwellFactory.cpp:758-771; `to_cut_algorithm`, en_CutAlgorithm.hpp:48).
2. `CutFactory::run()` — cl_CutFactory.cpp:116 — edges (:123-124), 3D-only faces
   (:125-129), periodicity backup (:132-136), curve `assign_edges` (:139-142).
3. `compute_poisson_problem()` (:197) or `compute_rcm_problem()` (:272) — node/edge ordering
   (collective; the only part non-root ranks join, `create_cuts_sub_slave`
   cl_MaxwellFactory.cpp:868-893); `rearrange_nodes/edges` (:300/:333).
4. `compute_cohomologies()` — cl_CutFactory.cpp:370:
   - `Homology::suggest_Homology()` (cl_Homology.cpp:234) — suggested H1 generators from
     terminal surfaces; 3D flags terminal sidesets + faces (:270-287, :331-399, output
     terminals :419-432);
   - `Homology::reorient_generators()` (:687) — 3D multiplies by −1 (:693-697);
   - flag φ-block elements/nodes/edges/faces (cl_CutFactory.cpp:415-426);
   - `unflag_symmetry_sidesets()` (:429, def :2843);
   - `SimplicialComplex` (:430; cl_SimplicialComplex.cpp:109 — 1-cell coboundary via faces
     :259-278, 2-cells from faces :311-393, top cells :396-466);
   - algorithm switch (:432-505): Pellikka (:1233), CCR (:507), PellikkaGeneralized
     (:1261), BeltedTree (cl_BeltedTree.cpp:24 + `Cohomology` cl_Cohomology.cpp:49). The
     `reduce_*` homology variants are unreachable (`mSuggestHomologies` hardcoded `true`,
     cl_CutFactory.hpp:87);
   - `Cohomology` ctor (cl_Cohomology.cpp:34) → coboundary matrix
     (cl_SimplicialComplex.cpp:1380), `cohomologyGroupOfChainComplex` (:91) → Smith normal
     form (fn_Smith.hpp:336/393/398), `generatorsOfCohomology` (:150), `clean()` (:195);
   - `updatekGeneratorsFromHomology()` (cl_Cohomology.cpp:420, called
     cl_CutFactory.cpp:516) — recombines the basis per terminal condition; 3D orders
     generators [in0,out0,in1,out1,…] (:507-519).
5. Generator→cut conversion and node duplication:
   `compute_thin_cuts_and_duplicate_interface_nodes()` (cl_CutFactory.cpp:841) → the
   `CutProcessor` pipeline (shared with Stage 3, see there), then
   `link_node_duplicates_and_originals()` (:2706), `InterfaceProcessor`
   (cl_InterfaceProcessor.cpp:262 — node-based, no dimension branches), and
   `collect_orphan_nodes()` (:2807).
6. Tail: `restore_thin_shell_sidesets()` (:876), `compute_element_adjacencies()` (:2914),
   refinalize (:181-184).

### 2D current state

| # | 3D procedure | 2D counterpart | dispatch | classification |
|---|---|---|---|---|
| 1 | orchestration `create_cuts`/`sub_master`/`sub_slave` | same code | none | OK |
| 2 | face creation in `run()` (:125-129) | skipped | runtime if | OK |
| 3 | Poisson/RCM + rearrange | same code | none | OK |
| 4 | `suggest_Homology` 3D (cl_Homology.cpp:270-287, 331-399) | 2D block terminals, input-half only (:259-267, :406) | runtime if :262/:331/:419 | PRESENT-BUT-UNTESTED (convention risk, below) |
| 5 | `reorient_generators` 3D ·(−1) (:693-697) | 2D ·(+1) no-op (:698-702) | runtime if | PRESENT-BUT-UNTESTED (empirical sign, below) |
| 6 | φ-flagging incl. faces (cl_CutFactory.cpp:403-427) | 2D twin without faces (:378-401) | runtime if :378 | OK |
| 7 | `unflag_symmetry_sidesets` face sweep (:2870-2892) | none — face loop empty in 2D | none (silent) | PRESENT-BUT-WRONG for 2D symmetry models |
| 8 | `SimplicialComplex::create_complex` 3D | 2D: edge coboundary via elements (:232-256), top cells = elements | `tDim` ifs :232/:311/:402/:433 | OK |
| 9 | coreduction algorithms | same code; empty k=3 maps handled (:1385/:1392 pad) | implicit | OK |
| 10 | `BeltedTree` (cl_BeltedTree.cpp:24) | 2D top-cell branch (:193-200) | runtime if :193 | PRESENT-BUT-UNTESTED (TRI/TET-only edge assumptions :229-231/:277-287) |
| 11 | `Cohomology` Smith/generators/clean | same code | none | OK |
| 12 | `updatekGeneratorsFromHomology` in/out pairing (:513-519) | 2D one-generator-per-condition (:507-508, :521-523) | `tIs3D` :507 | OK (documented) |
| 13 | `compute_element_adjacencies` (cl_CutFactory.cpp:2914) | none — faces-only | none (silent) | PRESENT-BUT-WRONG in 2D |

Key details:

- **2D terminal convention in `suggest_Homology`** `[INFERENCE, medium]`: the 2D branch
  flags only the first half of each BC domain list (`j < aTerminals(i).size()/2`,
  cl_Homology.cpp:260; comment "Only input terminals in 2-D" :259). If a 2D current BC
  lists a single block, the suggested generator is empty → cascades into
  `BELFEM_ERROR("Inconsistency in the definition of the terminals…")`
  (cl_Cohomology.cpp:563-565) or the empty-generator asserts
  (cl_CutData.cpp:266/326). Whether the `[go-blocks | return-blocks]` convention matches
  what the input parser produces was not confirmed — check `create_terminal_list`
  (cl_MaxwellFactory.cpp:979) against a real 2D input deck.
- **Sign handling is empirical, not derived**: `reorient_generators` comments admit "it
  should be +1.0 here, but there must be another sign error elsewhere in the code in 3-D"
  (:696) and "somehow the sign is good in 2-D" (:701). The 2D branch dodges a hack rather
  than resting on a verified convention.
- **`unflag_symmetry_sidesets` has no 2D consistency sweep**: the face loop
  (cl_CutFactory.cpp:2875-2892) is a silent no-op in 2D; the analogous cleanup (unflag top
  *elements* whose edges were all unflagged) does not exist. A 2D φ-element adjacent to a
  symmetry sideset stays flagged with dangling boundary entries in the complex before
  Smith/coreduction runs.
- **`compute_element_adjacencies` poisoning (re-verified end to end)**: the function builds
  element↔element adjacency solely from `mMesh->faces()` (cl_CutFactory.cpp:2917-2938),
  empty in 2D, then stamps `set_connectivity(Connectivity::ElementToElement)` (:2942).
  `Mesh::unfinalize()` deliberately never resets that flag (cl_Mesh.cpp:874-885, reset code
  commented out). `ConnectivityCalculator::connect_elements_to_elements()` then skips the
  correct node-based computation (cl_Mesh_ConnectivityCalculator.cpp:1026-1040; the
  early-out at :1014-1017 also fires once TS connectivity exists). Confirmed consumer:
  `Kernel::distribute_mesh()` calls `connect_elements_to_elements()` at
  cl_FEM_Kernel.cpp:627, :639, :670 ahead of `mesh::Distributor` — so a parallel 2D run
  partitions on an empty element graph `[INFERENCE that the Distributor actually requires
  this graph — confirm in cl_Mesh_Distributor.cpp]`. The comment at
  cl_CutFactory.cpp:2940-2941 claims the thin-shell factory needs the data, but the Stage-4
  audit found no direct adjacency consumer in `ThinShellFactory` (the TS-element adjacency
  is computed separately and correctly at cl_Mesh_ConnectivityCalculator.cpp:1046-1074).
- **Housekeeping (both dims)**: `CutFactory::mCuts` is never populated, so the
  `set_domain_type(DomainType::Cut)` loop at cl_MaxwellFactory.cpp:825-828 is a no-op;
  `[INFERENCE]` cut sidesets presumably get their type via the `cut_NN` label
  (en_DomainType.cpp:122-125) — the call site that applies this was not located. Stray pasted
  block-scope declaration `void collect_thin_cut_edges(...)` sits inside
  `determine_cut_case_3d` (cl_CutData.cpp:459-463; legal C++, declared function does not
  exist — paste accident on the live 3D path).

### Gap

- Symmetry-sideset consistency sweep missing in 2D (cl_CutFactory.cpp:2875).
- 2D terminal convention unvalidated (cl_Homology.cpp:260) `[INFERENCE]`.
- `compute_element_adjacencies` needs an edge/facet-based 2D branch, or must not stamp the
  flag in 2D (cl_CutFactory.cpp:2914-2943).
- Housekeeping: `mCuts` no-op loop; stray declaration cl_CutData.cpp:459-463; unresolved
  3D sign hack (cl_Homology.cpp:696) that 2D merely happens to dodge.

---

## Stage 3 — Thin-cut computation

### 3D reference sequence

Shell-opening entry (`create_thin_shell_cuts`) is enumerated in Stage 1 (items 9.x).
The thin-cut computation proper runs inside
`CutFactory::compute_thin_cuts_and_duplicate_interface_nodes()` (cl_CutFactory.cpp:841):

1. `CutProcessor` ctor — cl_CutProcessor.cpp:21 — per-generator `CutData` (:43-48), then:
   `collect_edges` :156 (→ `CutData::collect_edges` cl_CutData.cpp:186,
   `collect_coefficients` :231), `collect_elements` :249 (→ `CutData::collect_elements`
   cl_CutData.cpp:264; cut case via `determine_cut_case_3d` :443), `collect_facets` :307
   (3D branch :415-582 — flags faces per cut case, removes self-canceling faces, removes
   φ-boundary/periodic faces :420-423/:525, iterative dangling-face pruning :505-553),
   `collect_nodes` :656 (→ `CutData::flag_nodes` :691-715), `determine_cut_sets` :771,
   `create_cut_sets` :885, `compute_edge_bitsets` :210, `compute_node_bitsets` :589
   (TET4/TET10 kernels :968/:1062), `classify_periodic_pairs` :704,
   `create_abstract_nodes` :872, `duplicate_nodes` :1180 (→ `CutSet::create_duplicates`
   cl_CutSet.cpp:42), `relink_elements` :1191, `collect_duplicates` :1281,
   `create_thin_cut_sidesets` :899 (→ `CutData::add_thin_cut_sidesets_to_mesh`
   cl_CutData.cpp:55; 3D branch :143-178 emits **master+slave** facets, with the Step-5a
   slave-only fallback :159-172).
2. `link_node_duplicates_and_originals()` (cl_CutFactory.cpp:2706); `InterfaceProcessor`
   (:861); `collect_orphan_nodes()` (:2807).

### 2D current state

**Central question answered definitively: 2D thin-cut computation is neither a silent
no-op nor a stub.** A genuine 2D implementation exists and is reached (runtime branches at
cl_CutFactory.cpp:1630/1664; cl_CutData.cpp:69/333/666/729/783; cl_CutProcessor.cpp:317).
The shell-opening routines operate generically on `Facet*`/`Node*` and process LINE facets
correctly.

| # | 3D procedure | 2D counterpart | dispatch | classification |
|---|---|---|---|---|
| 1 | `determine_cut_case_3d` (cl_CutData.cpp:443) | `determine_cut_case_2d` (:365; case table 3/5/6 :387-418; coefficient assert :431) | runtime if :333 | OK (TRI only; QUAD → BELFEM_ERROR :415) |
| 2 | `collect_facets` 3D (cl_CutProcessor.cpp:415-582) | 2D branch (:317-414) | runtime if :317 | PRESENT-BUT-WRONG — two omissions (below) |
| 3 | `flag_nodes` 3D (cl_CutData.cpp:691-715) | 2D (:666-690) | runtime if :666 | OK (incl. periodic partners) |
| 4 | node-bitset kernels TET4/TET10 | TRI3/TRI6 (:924/:1015) | `mElementType` switch (:606/:788/:1206) | OK (non-TRI → "This should not happen" :642/:828/:1246) |
| 5 | `add_thin_cut_sidesets_to_mesh` 3D (cl_CutData.cpp:143-178) | 2D (:69-142) | runtime if :69 | PRESENT-BUT-WRONG — asymmetric (below) |
| 6 | `CutSet::create_duplicates`, `classify_periodic_pairs`, `InterfaceProcessor` | same code | none | OK (node-based) |

Key details:

- **2D `collect_facets` omissions** (both confirmed by two auditors):
  (a) it excludes only `mPhiBoundaries` edges (cl_CutProcessor.cpp:376-383), while 3D
  protects `mPhiBoundariesAndPeriodic` (:420-423, :525) — **periodic sidesets are unhandled
  in 2D thin cuts**; (b) there is no 2D counterpart of the 3D iterative dangling-face
  pruning loop (:505-553) — dead-end thin-cut edges pass through unpruned and later trip
  the cl_CutData.cpp:89-90 asserts in debug, or silently mis-select masters in release.
- **2D thin-cut facets are master-only**: the 2D branch sets only
  `tFacet->set_master(...)` (cl_CutData.cpp:139), with side identification via
  coordinate-epsilon matching (:98-134; comment :94 "something is off with the original
  pointers"), and asserts `tEdge->number_of_elements()==2` (:89) with **no boundary
  fallback** — the 3D Step-5a slave-only fallback (:159-172) has no 2D twin. Cuts touching
  a non-φ boundary abort in debug or corrupt state in release. Consumers expecting master+slave
  pairs on cut facets see `slave()==nullptr` on every 2D cut facet `[INFERENCE on whether
  any live consumer dereferences slave() on Cut sidesets in 2D — check the Cut sideset IWG
  path]`.
- Second-order T-matrix weights are unfinished for both dimensions:
  `// todo: fix weights for second order!` (cl_CutProcessor.cpp:97) — relevant to 2D TRI6.
- Element-type inference is dim+order only (cl_CutProcessor.cpp:33-35): a 2D mesh with QUAD
  φ-blocks would run the TRI3 bitset kernels on 4-edge elements and produce wrong bitsets
  silently before any error fires.

### Gap

- Port periodic-edge protection and dangling-edge pruning into the 2D `collect_facets`.
- Symmetrize `add_thin_cut_sidesets_to_mesh` 2D: master+slave facets, topological side
  selection, boundary (Step-5a-style) fallback.
- 2D terminal-curve provisioning (shared with Stage 1, G1).

---

## Stage 4 — Thin-shell insertion

### 3D reference sequence

1. `create_magnetic_kernel()` — cl_MaxwellFactory.cpp:425 — fresh path `create_cuts()`
   (:432); then `create_edges_and_faces_on_mesh()` (:457), `create_thinshells()` (:458),
   and (rank 0, fresh path) `create_hanging_edges_and_facets()` (:483).
2. `create_edges_and_faces_on_mesh()` — :1139 — `mMesh->create_edges(...)` on Conductor
   groups (:1145); `create_faces(...)` only if `max_element_order() > 1` (:1151).
3. `create_thinshells()` — :898 — .bfm-reload short-circuit (:902-923); `ThinShellFactory`
   ctor (:924, cl_ThinShellFactory.cpp:35); per protoshell `tFactory.create()` (:936);
   material↔block registration (:948-959); unfinalize/finalize (:963-964).
4. `ThinShellFactory::create()` — cl_ThinShellFactory.cpp:108:
   `collect_sidesets` :440 / `collect_facets` :462 (protoshell "shell_NN" sidesets are
   emptied — husks by design, :504) / `collect_nodes` :517; normals via ElementType switch
   :131-157 → `process_nodes_tri3` :726 / `process_nodes_tri6` :833;
   `compute_distances` :989; `create_nodes_on_layers` :1027 (periodic backup :1131-1176);
   `create_temporary_edges` :1391; node table updates :2110/:2147; element creation switch
   :204-240 → `create_elements_on_blocks_tri3` :1285 (PENTA6TS :1308) / `_tri6` :1329
   (PENTA18TS :1356); `create_edges_on_layers` :1514 (+ TRI6: faces :1560/:1679);
   `link_elements_with_edges` :1600; blocks tagged `DomainType::ThinShell` :261-271;
   `create_buffers` :1704; `create_ghost_facets` :278/:1825 (Ghost sideset :288-304);
   periodic `create_periodic_sideset` :310-311/:2304; `new ThinShell(...)` :325
   (cl_ThinShell.cpp:17 — hides both sidesets); layer entities moved into the mesh
   :348-367.
5. `create_hanging_edges_and_facets()` — :1163 — PART 1 interface hanging (:1181-1267,
   deactivates ThinShell husk sidesets :1257-1260); PART 2 per-shell edge hanging
   (:1285-1367) via `hang_thinshell_edges_on_nodes_bottom` :1429 / `_on_nodes_top` :1489 /
   `_on_edges_bottom` :1549 / `_on_edges_top` :1638, built on
   `mesh::get_{bottom,top}_{nodes,edges}` (meshtools.cpp:1219/1309/1336/1377); PART 3
   quadratic facet weights, gated `dim==3 && order==2` (:1372-1373).
6. DOF-manager consumption: `BlockData::create_blocks` (cl_FEM_DofMgr_BlockData.cpp:69),
   `collect_thin_shell_facet_ids` :152, `link_thin_shell_facets` :387 (serial :408 /
   parallel :431); hanging-DOF weights `DofData::collect_hanging_dofs`
   (cl_FEM_DofMgr_DofData.cpp:2328) → `create_dofwise_t_matrices_master` :3471
   (edge-on-node :3555-3671; edge-on-edge incl. LINE2/LINE3 cascade :3673-3902).

### 2D current state

| # | 3D procedure | 2D counterpart | dispatch | classification |
|---|---|---|---|---|
| 1 | `process_nodes_tri3` (:726) | `process_nodes_line2` (:561) | ElementType switch :131-157 | PRESENT-BUT-UNTESTED (mirrors tri3, incl. duplicate handling :613-650) |
| 2 | `process_nodes_tri6` (:833) | `process_nodes_line3` (:653) | same switch | PRESENT-BUT-WRONG — transposed indexing (below) |
| 3 | `create_elements_on_blocks_tri3` (:1285) | `_line2` (:1191; QUAD4TS :1216-1221) | ElementType switch :204-240 | PRESENT-BUT-UNTESTED (QUAD4TS constructible: cl_Element_Factory.cpp:169-172/:456) |
| 4 | `_tri6` (:1329) | `_line3` (:1234; QUAD9TS :1262-1272) | same switch | PRESENT-BUT-WRONG upstream (rows 2, and Stage-5 EF gap) |
| 5 | `create_edges_on_layers`/`link_elements_with_edges` | same code, order-1 branch :1609-1635 | order switch | OK |
| 6 | faces on layers (TRI6 only :227-234) | correctly not needed | ElementType switch | OK |
| 7 | `create_buffers` (:1704) | same code | none | PRESENT-BUT-UNTESTED |
| 8 | `create_ghost_facets` (:1825) | same code | none | PRESENT-BUT-UNTESTED (`top/bottom_facet_index` QUAD entries in meshtools not verified) |
| 9 | `create_periodic_sideset` (:2304) | none | none | MISSING — PENTA-only (below) |
| 10 | `hang_thinshell_edges_on_nodes_bottom` (:1429) | same code; meshtools QUAD cases :1223-1229/:1340-1350 | type switches in meshtools | PRESENT-BUT-UNTESTED |
| 11 | `hang_thinshell_edges_on_nodes_top` (:1489) | same code; `get_top_nodes(QUAD)` = facet 2 (meshtools.cpp:1313-1316) | type switch | PRESENT-BUT-WRONG `[INFERENCE]` — cross-wired pairing (below) |
| 12 | `hang_thinshell_edges_on_edges_top` (:1638) | same code | none | PRESENT-BUT-WRONG — hard abort (below) |
| 13 | `hang_thinshell_edges_on_edges_bottom` (:1549) | same code | none | PRESENT-BUT-UNTESTED |
| 14 | PART 3 quadratic facet sources (:1372-1414) | none | `dim==3 && order==2` gate | MISSING (moot while 2D is linear-only) |
| 15 | `mesh::ThinShell`, `BlockData` linking | same code | none | OK |
| 16 | DOF type request: 2D shells use **edge_h** like 3D | same lists (cl_IWG_Maxwell.cpp:54-75; cl_Maxwell_FieldList.cpp:256-263, :388-394, :411-417) | none | OK (confirmed: 2D h-φ shells are edge-DOF based, not node-h) |

Key details:

- **`process_nodes_line3` is guaranteed broken**: allocates
  `aNodeNormals.set_size(3, aNodes.size())` (cl_ThinShellFactory.cpp:671) but accumulates
  `aNodeNormals(i,0) += …` with `i` = node index (rows, :694-696/:703-705/:711-713) —
  transposed relative to the column-per-node layout the read side uses
  (`.col(tNode->index())` into a length-2 vector, :719, :669). Out-of-bounds for node
  index ≥ 3 and a shape mismatch regardless. Blocks LINE3 shells at the first step.
- **2D edge-hanging, conductor above the shell → hard abort**:
  `hang_thinshell_edges_on_edges_top` calls
  `to_master_orientation(Facet*, Cell<Edge*>&, …)` (cl_MaxwellFactory.cpp:1673) on the
  slave facet's edge set; for a 2D TRI3 slave, `get_edges_of_facet` returns 1 edge
  (cl_Element_TRI3.hpp:90-118) and the size switch handles only 3 and 4 →
  `BELFEM_ERROR("Invalid number of edges for facet")`
  (fn_to_master_orientation.cpp:676-681). Needs a size-1 LINE case (identity + sign).
- **2D edge-hanging, air above → suspected cross-wired pairing** `[INFERENCE, medium-high]`:
  `get_top_nodes(QUAD)` returns facet 2 = {node2, node3} = {top(F1), top(F0)}
  (facet-CCW order, cl_Element_QUAD4TS.hpp:74-79), while the PENTA top facet {3,4,5} is
  position-aligned with the bottom (cl_Element_PENTA6TS.hpp:100-106).
  With `to_master_orientation(LINE2)` mapping volume nodes to {F0,F1}
  (fn_to_master_orientation.cpp:32-37), the lookup
  `tNodesOnVolume(tEdge->node(k)->index())` (cl_MaxwellFactory.cpp:1531) pairs top(F1) with
  the volume node at F0. Confirm with a two-element 2D mesh before fixing.
- **`create_periodic_sideset` is PENTA-only**: `BELFEM_ASSERT(geometry_type(...)==PENTA)`
  (cl_ThinShellFactory.cpp:2320) and a `f<3` side-facet loop (:2324). 2D + periodic
  asserts in debug and mis-enumerates QUAD4TS facets in release. Reached whenever
  `mMesh->has_periodicity()` (:306-322).
- **Legacy dead path (conflict resolved between auditors, re-verified)**: the layers-Element
  ctor (cl_FEM_Element.cpp:89-110) → `link_dofs_thin_shell` (:901) →
  `compute_edge_directions_thinshell` (:1381, with its `dims==2` "tape roller" assert
  :1383-1385 and `physical_tag` read :1388) has **no caller**: `new Element(` exists only at
  cl_FEM_Block.cpp:68/:80 and cl_FEM_SideSet.cpp:230, and the SideSet call passes
  (mode, master, slave) — not a layers cell. The "tape roller" mechanism it references does
  not exist in the tree. This is DEAD code from the pre-refactor sideset-based shell
  assembly; its 2D-only assert must not be read as evidence of a working 2D path. The live
  question it leaves behind: how QUAD4TS layer-element `edge_directions()` are set for
  `EF_QUAD4TS::link()` (Stage 5).
- **Silent skip (both dims)**: `if (tFacets.size() == 0) break;` at
  cl_MaxwellFactory.cpp:1288 is a `break`, not `continue` — one empty ThinShell silently
  skips edge-hanging for all remaining shells.
- 3D-side bug noted in passing (out of scope, report-only): `process_nodes_tri6` writes
  `tM(0)` three times instead of `tM(0)/tM(1)/tM(2)` (cl_ThinShellFactory.cpp:913-915).
- Side-connector residue (report-only, removal pending by design): `EF_HEX8TS`
  (cl_EdgeFunctionFactory.cpp:62-64), HEX8TS cases in meshtools.cpp:1263-1272/:1360-1367,
  cl_FEM_Calculator.cpp:369/:675, fn_num_nedelec_dofs.hpp:34, cl_Element_Factory.cpp:48/:185,
  `ThinShellFactory::compute_binomial_vectors` (:1964-2107) and `compute_side_indices`
  (:76) — no callers.

### Gap

- Fix `process_nodes_line3` indexing (prerequisite for any LINE3 shell).
- Add the LINE-facet (size-1) case to `to_master_orientation`; verify/fix the top-node
  pairing convention for QUAD4TS.
- Add a QUAD4TS branch to `create_periodic_sideset`.
- `break`→`continue` at cl_MaxwellFactory.cpp:1288; stale "only linear triangular" error
  text at :1291-1293 (LINE2 is accepted; message claims otherwise).

---

## Stage 5 — Main solver + preprocessors

Architecture note (matters for reading older docs): `cl_IWG_Maxwell.cpp` is now a thin
dispatcher that binds one function pointer `mFunMKF` per group (`link_to_group`,
cl_IWG_Maxwell.cpp:239-567); all element matrices live in `matrices/mt_maxwell_*.cpp`
operating on the `Calculator` abstraction. Thin shells are assembled on the **layer blocks**
(PENTA6TS / QUAD4TS, `DomainType::ThinShell`) via `h_ts_*`; the tape facet sideset is
`GeometryOnly` (normal-field calculator only, `init_activation_maps`
cl_IWG_Maxwell.cpp:109-118); a `DomainType::Ghost` sideset couples adjacent layers by
penalty (`h_ghost`).

### 3D reference sequence

1. `create_equation` — cl_MaxwellFactory.cpp:408-420 — `ModelDimensionality` from mesh dim;
   `IWG_Maxwell` ctor (cl_IWG_Maxwell.cpp:31-87) declares label-based DOFs (Air "phi",
   Conductor "edge_h" (+"face_h" order 2), Ferro "phi").
2. `initialize` — cl_IWG_Maxwell.cpp:123-219 — `FieldList::initialize`
   (cl_Maxwell_FieldList.cpp:44-232); `collect_block_dofs` (:237-294 — ThinShell blocks get
   the Conductor table); `collect_sideset_dofs` (:299-432 — ThinShell :388, Ghost :411,
   Cut :381, GeometryOnly/Inactive skip :418-423).
3. DOF creation: `DofData::create_dofs` (cl_FEM_DofMgr_DofData.cpp:108); `count_edge_dofs`
   (:1267, `mesh::number_of_edges(type)` — type-driven).
4. Hanging-DOF preprocessing: Stage-4 edge hanging → `collect_hanging_dofs` (:2328) →
   `create_dofwise_t_matrices_master` (:3471; node 1:1 :3491-3553; edge-on-node with 2/3-
   node tables :3555-3671; edge-on-edge + cascade, all LINE2/LINE3 combos :3673-3902).
5. Element/DOF linking: `SideSet::initialize_elements` (cl_FEM_SideSet.cpp:142-238) →
   `Element` ctor (cl_FEM_Element.cpp:56-85); DOF-count formulas for Ghost/ThinShell
   sidesets in `IWG::number_of_dofs_per_element` (cl_IWG.cpp:650-690, type-driven).
6. Calculator setup: edge functions (cl_FEM_Calculator.cpp:256-335; per-facet master :299,
   per-(facet,orientation) slave :326); nedelec-data dispatch :341-380
   (PENTA6TS/QUAD4TS/HEX8TS → `nedelec_data_linear_h` :367-375); dV dispatch :618-696
   (3D PENTA6TS → `dV_ts` :663); normals :767-822 (PENTA :815).
7. Assembly loop: `IWG_Timestep::compute_jacobian_and_rhs` (cl_IWG_Timestep.cpp:370-386) →
   `IWG_Maxwell::compute_mkf` (cl_IWG_Maxwell.cpp:224-234) → `(*mFunMKF)(...)`. Binding
   table (`link_to_group`): Air TRI3/TET4/TRI6/TET10 → `phi_*` (:257-285); Conductor →
   `h_*` (:306-397); ThinShell blocks → `h_ts_*` by material (:437-527); Ghost →
   `h_ghost` (:555-559); default → `BELFEM_ERROR("Not implemented")` (:560-564).
8. Thin-shell matrices: `maxwell::h_ts` (mt_maxwell_h.cpp:2216-2502), `h_ts_metal` (:92),
   `h_ts_hts*` (:468-1765); normal flux `compute_bn` (mt_maxwell_h.hpp:124-167) →
   `Calculator::get_normal_calculator` (cl_FEM_Calculator.cpp:1094-1146, looks up the
   GeometryOnly tape sideset).
9. Ghost penalty: `maxwell::h_ghost` (mt_maxwell_h.cpp:1766-1937); scratch matrices from
   `create_custom_vectors_and_matrices` (cl_IWG_Maxwell.cpp:601-615).
10. BCs: `MaxwellBoundaryConditionFactory` (cl_MaxwellBoundaryConditionFactory.cpp:23-411);
    terminal plumbing `create_terminal_list` (cl_MaxwellFactory.cpp:979-1061); currents via
    abstract-node DOFs (`set_currents` cl_IWG_Maxwell.cpp:92-104,
    `collect_abstract_node_dofs` :720-742); Dirichlet/symmetry wiring
    (cl_MaxwellFactory.cpp:626-678).
11. Timestepping: `IWG_Timestep::{explicit_euler,crank_nicolson,galerkin,bdf1}`
    (cl_IWG_Timestep.cpp:409-480) — matrix-size-agnostic.

### 2D current state

| # | 3D procedure | 2D counterpart | dispatch | classification |
|---|---|---|---|---|
| 1 | `create_equation`, IWG ctor, FieldList | same code (TwoD branch :413-415) | runtime | OK |
| 2 | air/interface/symmetry bindings | `phi_tri3` :262-265, `phi_phi_2d` :295-298, `symmetry_phi_2d` :423-425, `h_symmetry_2d` :544-547 | ElementType / facet dim | OK (bindings exist; bodies not audited) |
| 3 | EF_PENTA6TS | EF_QUAD4TS (cl_EdgeFunctionFactory.cpp:54-56; cl_EF_QUAD4TS.cpp:18-209) | ElementType switch | PRESENT-BUT-UNTESTED, two suspicions (below) |
| 4 | Calculator PENTA6TS branches | QUAD4TS: `dV_ts` :631-634, `normal_quad_straight` :785/:1519-1571, nedelec :368, `slave_integration_index_2d` :398-403; `number_of_orientations(QUAD4TS)=2` (meshtools.cpp:884-892) | ModelDimensionality + type switches | PRESENT-BUT-UNTESTED |
| 5 | `h_ts_*` block assembly | same code via Calculator | none | PRESENT-BUT-WRONG (suspected size bug, below) |
| 6 | `compute_bn`/normal calculator | same code; 2D LINE2-on-TRI3 normals implemented (:1403-1449) | GeometryType switch | PRESENT-BUT-UNTESTED |
| 7 | `h_ghost` + scratch matrices | none — 3D hard-coded | none | PRESENT-BUT-WRONG — the hardest 2D blocker (below) |
| 8 | Ghost/ThinShell DOF-count formulas (cl_IWG.cpp:650-690) | same code | type-driven | OK |
| 9 | hanging-DOF T-matrices | same code; LINE2/LINE3 combos exist (:3692/:3716/:3833) | entity-type switches | PRESENT-BUT-UNTESTED |
| 10 | BC factory | 2D = "no output terminal keys" path (:123-138) | convention, not checked | PRESENT-BUT-UNTESTED (no dim cross-check) |
| 11 | currents via abstract nodes; timestepping; postproc gating (Hz/Bz skipped, cl_MaxwellFactory.cpp:1902-1915) | same code | none / runtime | OK |

Key details:

- **`h_ghost` is 3D-hard-coded (hardest single 2D blocker)**: asserts facet `TRI3` (:1771)
  and master/slave `PENTA6TS` (mt_maxwell_h.cpp:1873-1875), then uses hard-coded
  `n=6; m=3; d=3` (:1878-1880, comment "if you want other elements, adapt the numbers
  below"). QUAD4TS edge matrices `Em/Es` are 2×2 (`mE` allocation, cl_EF_QUAD4TS.cpp:21;
  2 Nedelec dofs per fn_num_nedelec_dofs.hpp:28), so in release builds (asserts compiled
  out, unchecked `Matrix::operator()`) the loops read/write out of bounds — memory
  corruption without a diagnostic; debug aborts. The scratch allocation feeding it is
  equally hard-coded: `e=6`, `K±±` 6×6, `D±` 3×6
  (`create_custom_vectors_and_matrices`, cl_IWG_Maxwell.cpp:601-615, re-verified).
  Every multi-layer 2D shell assembles Ghost groups. `[INFERENCE, unverified]`: a
  single-layer shell might have an empty Ghost sideset and dodge this — check
  `mCreateGhostFacets` logic (cl_ThinShellFactory.cpp:175-303).
- **Thin-shell work-vector size mixing (re-verified)**: `create_custom_vectors_and_matrices`
  sizes `"j","b","h","bm","bs"` as `mesh::dimension(tType)` (= 2 for QUAD4TS) but
  `"bt","bn","n"` as hard-coded 3 (cl_IWG_Maxwell.cpp:585-598). `h_ts_*` mixes them
  (`b = bt + bn`, mt_maxwell_h.cpp:2271 and :137). `[INFERENCE]` outcome depends on
  `Vector` assignment/resize semantics — either a debug length assert or a silently
  inconsistent |b|/β; confirm with a 2D debug run or the linalg backend's `operator=`.
- **EF_QUAD4TS is unvalidated** `[INFERENCE]`: (a) `E()` applies `-mS[1]` citing
  "mS[1] = -mS[0] set in link()" (cl_EF_QUAD4TS.cpp:185-187), but `link()` (:71) takes
  signs verbatim from `aElement->edge_directions(mS)` — nothing enforces the opposite-sign
  invariant, and how QUAD4TS layer elements get their edge directions in the live pipeline
  is unresolved (the legacy tape-roller path is dead, Stage 4). (b) curl parameters
  `B(0,0)=B(0,1)=+0.5` (:157-167, "symmetric current") lack a cited derivation — compare
  against the PENTA6TS reduction / Messe et al. 2023 (paper1). Needs an analytic patch
  test.
- **2D is linear-only, enforced at five independent points**:
  `BELFEM_ERROR(max_element_order()==1)` (cl_MaxwellFactory.cpp:1176);
  `EF_QUAD9TS`/`EF_PENTA18TS` missing (cl_EdgeFunctionFactory.cpp:66-82);
  "Higher order thin shells are not implemented!" (cl_FEM_Calculator.cpp:637);
  the TRI6/TET10-only `TMatrix` (cl_Maxwell_TMatrix.cpp:19-112) with the 3D-only PART 3
  (cl_MaxwellFactory.cpp:1372-1414); the 3D-only facet-source branch
  (cl_FEM_DofMgr_DofData.cpp:2339-2391, todos :2354-2355). Plus the Stage-4
  `process_nodes_line3` bug and the Stage-1 LINE3 loop bug. Any order-2 2D effort must
  clear all of these; until then treat 2D as linear-only.
- BC factory stores `mNumDimensions` (cl_MaxwellBoundaryConditionFactory.cpp:24) but never
  branches on it; 2D-ness is inferred from absent "output terminal" keys (:123-138, voltage
  needs `length`). Misconfiguration fails late (at `link_to_group`'s default
  `BELFEM_ERROR`, cl_IWG_Maxwell.cpp:560-564) or scales silently.
- Latent: `EF_LINE3::...` contains `BELFEM_ERROR(false, "function not implemented")`
  (cl_EF_LINE3.cpp:180) — not reached via the factory today.

### Gap

- `h_ghost` 2D variant + type-dependent scratch sizing (the release-mode memory-corruption
  path).
- `bt/bn/n` sizing consistency for TwoD thin-shell groups.
- EF_QUAD4TS validation (edge-direction/sign provenance + curl coefficients) — including
  answering where QUAD4TS `edge_directions()` are set in the live pipeline.
- (Deferred, order-2): EF_QUAD9TS, 2D TMatrix analog, 2D facet sources.

---

## Cross-stage dependency picture

```
Stage 1  G1 terminal-curve validation ──► Stage 3 thin cuts sane ──► currents correct
Stage 1  G2 LINE3 break ────────────────► any order-2 run possible (with all Stage-4/5 order-2 items)
Stage 2  adjacency 2D branch ───────────► parallel 2D runs (Distributor)
Stage 3  collect_facets periodic+prune ─► 2D periodic / robust cuts
Stage 4  process_nodes_line3 ───────────► LINE3 shells (order-2 chain)
Stage 4  to_master_orientation LINE ────► 2D shells with conductor neighbors
Stage 4  top-node pairing ──────────────► correct hanging-DOF T-matrices ─► Stage 5 solve
Stage 5  h_ghost 2D + scratch sizing ───► multi-layer 2D shells assemble at all
Stage 5  bt/bn/n sizing + EF_QUAD4TS ───► 2D thin-shell physics correct
```

Minimum path to "a linear, serial, non-periodic 2D thin-shell case runs": G1 (or a
correctly-authored input deck), Stage-4 `to_master_orientation` LINE case + top-pairing
verification, Stage-5 `h_ghost`/scratch + `bt/bn/n` sizing, EF_QUAD4TS validation.

---

## Appendix — files inspected

Coverage union of the five stage audits plus the synthesis re-verification pass.
(L…) = specific line ranges read rather than the full file.

**fem/maxwell**
- src/fem/maxwell/cl_MaxwellFactory.cpp (L55-450, 700-1800, 1860-2000, 2394-2553), .hpp
- src/fem/maxwell/cl_IWG_Maxwell.cpp (full), cl_IWG_Maxwell.hpp
- src/fem/maxwell/cl_Maxwell_FieldList.cpp (full)
- src/fem/maxwell/cl_Maxwell_TMatrix.cpp (full)
- src/fem/maxwell/cl_MaxwellBoundaryConditionFactory.cpp (full), .hpp
- src/fem/maxwell/matrices/mt_maxwell_h.hpp (full), mt_maxwell_h.cpp (L28-165, 1766-2526;
  remainder structure-scanned)
- src/fem/maxwell/en_Maxwell_Formulations.hpp, doc/*.md (grep-level)

**homology**
- src/homology/cl_CutFactory.cpp (full), .hpp
- src/homology/cl_CutData.cpp (full), .hpp (accessors)
- src/homology/cl_CutProcessor.cpp (full), .hpp (partial)
- src/homology/cl_CutSet.cpp (full)
- src/homology/cl_InterfaceProcessor.cpp (full)
- src/homology/cl_SideSetFactory.cpp (full)
- src/homology/cl_CutProcessorManual.cpp (full)
- src/homology/cl_Homology.cpp (full), .hpp (L1-80)
- src/homology/cl_Cohomology.cpp (full)
- src/homology/cl_SimplicialComplex.cpp (full), .hpp (targeted)
- src/homology/cl_BeltedTree.cpp (full)
- src/homology/cl_Chain.cpp (L1-120), cl_Cochain.cpp (full)
- src/homology/cl_Topology.cpp (full), .hpp
- src/homology/fn_Smith.cpp/.hpp (structure), en_CutAlgorithm.hpp (full)
- src/homology/doc/*.md (grep-level)

**mesh**
- src/mesh/cl_CurveFactory.cpp (full), .hpp; cl_Curve.cpp (L153-184), .hpp
- src/mesh/cl_Protoshell.hpp (full); cl_ThinShell.cpp/.hpp (full)
- src/mesh/cl_Facet.cpp (L36-66), .hpp
- src/mesh/cl_Mesh.cpp (L762-904, 2560-2600)
- src/mesh/cl_Mesh_ConnectivityCalculator.cpp (L1000-1090)
- src/mesh/meshtools.cpp (L860-960, 1190-1406)
- src/mesh/fn_to_master_orientation.cpp (L19-120, 586-685)
- src/mesh/cl_Element_QUAD4TS.hpp, cl_Element_QUAD9TS.hpp (full);
  cl_Element_PENTA6TS.hpp (L55-160); cl_Element_TRI3.hpp (L85-124);
  cl_Element_LINE2.hpp, cl_ElementTemplate.hpp (L895-915 + grep)
- src/mesh/cl_EdgeFactory.cpp (L530-610), cl_FaceFactory.cpp (L60-140)
- src/mesh/cl_Element_Factory.cpp, Mesh_Enums.hpp, cl_Element.hpp, en_DomainType.cpp
  (L100-140), cl_ProtoMesh.cpp (grep)

**fem/kernel + iwg + interpolation**
- src/fem/kernel/cl_ThinShellFactory.cpp (full), .hpp
- src/fem/kernel/cl_FEM_DofMgr_BlockData.cpp (L60-525)
- src/fem/kernel/cl_FEM_DofMgr_DofData.cpp (L1267-1400, 2280-2470, 3470-3910)
- src/fem/kernel/cl_FEM_Element.cpp (L40-240, 850-1450)
- src/fem/kernel/cl_FEM_SideSet.cpp (L100-320)
- src/fem/kernel/cl_FEM_Calculator.cpp (L126-420, 560-1146, 1400-1850, 2030-2106), .hpp
- src/fem/kernel/cl_FEM_Kernel.cpp (L190-310, 600-690)
- src/fem/kernel/cl_FEM_PhysicalBoundaryCondition.cpp (L76-135), .hpp
- src/fem/kernel/cl_FEM_Group.cpp, cl_FEM_KernelParameters.*, cl_MeshChecker.cpp,
  cl_Pipette.cpp (grep-level)
- src/fem/iwg/cl_IWG.cpp (L380-720), cl_IWG_Timestep.cpp (L340-486)
- src/fem/interpolation/cl_EdgeFunctionFactory.cpp (full)
- src/fem/interpolation/nedelec/cl_EF_QUAD4TS.cpp (full); cl_EF_LINE3.cpp, cl_EF_TRI3.cpp,
  cl_EF_PENTA6TS.cpp, fn_num_nedelec_dofs.hpp (targeted)
- src/fem/interpolation/cl_IF_InterpolationFunctionFactory.cpp, cl_IF_IntegrationData.cpp
  (grep-level)

Not audited (declared [INCOMPLETE] by the stage audits): `FaceFactory::create_faces_2d`
internals; `top_facet_index`/`bottom_facet_index` QUAD4TS entries in meshtools;
EF_QUAD4TS shape-function numerics; `DofMgr_SideSetData` thin-shell sideset assembly;
`ConnectivityCalculator::connect_facets_to_facets` 2D behavior; bodies of `phi_tri3`,
`phi_phi_2d`, `symmetry_phi_2d`, `h_symmetry_2d`; the L2 projection matrices
(`mt_maxwell_l2_*`); the hphirun runner/Controller; `mesh::Distributor`'s exact use of the
element graph.
