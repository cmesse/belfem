# Interface Node Duplication for Ferro-Air and Air-Coil Visualization

> **CLOSED 2026-09-03** (todo/ currentness sweep, round 3): R3/R5 landed and serially validated on the production coil/ferro model since 2026-07-03; the remnant is parallel validation, a run gate. Status lines and checkboxes below are as they stood at closure and are not maintained.

**Date:** 2026-07-03
**Purpose:** Extend the `InterfaceProcessor` node-duplication scheme (currently exercised for the
φ/h conductor interfaces) to ferro-air and air-coil interfaces, so ParaView output shows the
physical field jump at the iron surface and a cleanly separated (dead) coil region. The mechanism
is to duplicate the interface nodes: **ferro-air duplicates hang on their originals with weight 1**
(solve unchanged; side-local postproc recovers each side independently); **air-coil duplicates are fully
decoupled** (no sources — the coil is excluded from the computation, its nodal fields stay NaN).
**Module:** `src/homology` (primary), `src/fem/maxwell`, `src/fem/kernel` (postproc)
**AIs involved:** Claude (exploration + plan), Codex (trace `devlog/dl20260703_interface_duplicate_node_trace.md` + audit rounds 1-2)
**Status:** IMPLEMENTED + SERIAL-VALIDATED (2026-07-03) — R3 (coil decouple) and R5 (side-local
recovery with same-PP sibling merge) landed and confirmed on Christian's tape/coil model
(interfaces crisp, cuts continuous, coils dead); per-side ownership fix in `Mesh::update_ownerships()`
applied after Codex audit round 2. Remaining: R6 guard, R7/R8 parallel validation, re-scoped R4
(strip still active, but the serial run is physical — verify T-matrix build order), O4/O5.
See `devlog/dl20260703_interface_decouple_sidelocal_recovery.md`.

**Re-verified 2026-08-09 (currentness sweep): still accurate; nothing has been closed since,
and the open items are unchanged.** R6 (double-duplication guard), R7/R8 (parallel
validation, iron-yoke + coil reproducer on 2/4 procs) and the re-scoped R4 are all still
open — R7/R8 are run gates that only Christian can clear.
Anchor drift: the ferro source strip R4 refers to is now
`cl_MaxwellFactory.cpp:1938-1960` (`reset_source_container()` on flagged hanging ferro
nodes), not `:1882-1900`; the coil hard-error guard is `cl_IWG_Maxwell.cpp:338,343`, not
`:403-411`. `InterfaceProcessor::unify_duplicates()` is `cl_InterfaceProcessor.cpp:579`.
**Cross-campaign note worth flagging:** `reset_source_container()` clobbering is *also* an
open audit item on the side-connector track (`hex8tb_phase2_fem_wiring.md`, §6 gate) — the
same routine, reached from a different direction. Whoever picks up R4 should read that item
first; a change here could move the connector result and vice versa.

> **Scope guards (from the task brief, Christian 2026-07-03):**
> - Ferro-air: duplicate + weight-1 hang. "Does not make sense for the physics itself, but it
>   makes sense for the visualization where we process ferro and air fields independently."
> - Air-coil: duplicate and **decouple completely**. "The coil domains are to be excluded from
>   the calculation. If we want to be super clean, we can even fill the nodal values with NaN."
> - Coil interfaces stay geometry/postprocessing-only — `cl_IWG_Maxwell.cpp:403-411` hard-errors
>   if `InterfaceAirCoil`/`InterfaceFerroCoil` reaches assembly. This plan must not activate them.
> - Keep the conductor (φ/h) interface path handled by `CutProcessor` out of scope; do not regress it.

---

## 1. Current Behaviour and How It Fails

The duplicate → hang → uncouple pattern already existed for conductor interfaces:

| Stage | Mechanism | Citation |
|---|---|---|
| A. Mesh: duplicate + relink + hang | `InterfaceSet::duplicate_nodes()` / `relink_elements()` / `add_duplicate_nodes_to_mesh()` — duplicates of **non-hanging** originals get `set_sources( original, 1.0 )`; duplicates of hanging originals **copy the original's source list** instead *(corrected per Codex audit 2026-07-03)* | `cl_InterfaceProcessor.cpp:147-243, 466-489` (weight-1 branch `:470-474`, copy branch `:156-169, 475-488`) |
| A'. Register pairs | `unify_duplicates()` — **derives `set_original()` from the sources list** | `cl_InterfaceProcessor.cpp:546-558` |
| B. Solve: tie | hanging dof T-matrix built from `mesh_basis()->source(k)/weight(k)`, type-matched (φ→φ); hanging dof excluded from the system, recomputed post-solve as Σ weight·source | `cl_FEM_DofMgr_DofData.cpp:3471-3546`, `cl_FEM_DofMgr_SolverData.cpp:2524-2538` |
| C. Postproc: uncouple | **CORRECTED (Codex audit 2026-07-03, verified by Claude):** `recover_fields()` canonicalizes each element node to `original()` and adds the contribution to the original **and all registered duplicates** — it is NOT side-local. Any prior conductor-interface side separation had to come from per-domain-type postprocessor block selection and run order (resolved in O3) | `cl_FEM_Postprocessor.cpp:985-1042` (canonicalize `:987`, all-duplicates write `:999-1013`) |

Current ferro-air and air-coil blockers:

| Failure | Mechanism | Evidence |
|---|---|---|
| Ferro-side duplicates get severed | The "removes connections at interface" pass strips sources from **every** hanging node flagged in a ferro block — a weight-1-hung ferro-side interface duplicate is exactly such a node → it becomes a free φ dof (wrong physics) or dead (wrong viz) | `cl_MaxwellFactory.cpp:1882-1899` |
| Coil not explicit in admission | `InterfaceProcessor` seeds only Air and Ferro bitsets; the admission predicate (`:337-341`) *textually* passes air-coil (air on one side) but coil membership is never explicit, and buffer-coil is silently skipped | `cl_InterfaceProcessor.cpp:281-293, 336-341`; Codex trace dl20260703 |
| Coil duplicates would be hung, not decoupled | `add_duplicate_nodes_to_mesh()` unconditionally hangs every duplicate on its original (weight 1) — for a dead coil we want **no** sources at all | `cl_InterfaceProcessor.cpp:466-489` |
| Decoupling breaks pair registration | `unify_duplicates()` finds the original **through the sources list**; a sourceless duplicate would never get `set_original()` / `add_duplicate()`, so postproc would not see the pair | `cl_InterfaceProcessor.cpp:546-558` |
| Buffer ≠ Air inconsistency | `Topology::sideset_type()` treats Buffer as Air (`cl_Topology.cpp:288-289`), but `mIsAirBlock` is seeded from `groups(DomainType::Air)` only → buffer-ferro/buffer-coil interfaces classified differently by the two systems | `cl_InterfaceProcessor.cpp:286-289` |
| No side-local recovery | `recover_fields()` writes each element contribution to the original **and every registered duplicate** — a ferro-air split could not emerge from this code as originally written; whichever domain postprocessor writes the shared field last wins at the interface *(Codex audit 2026-07-03, finding 7 — the biggest blocker the draft plan missed)* | `cl_FEM_Postprocessor.cpp:985-1013` |

**Bottom line:** stage B (the solve-side tie) is generic and needs no new machinery; stage C is
**not**, because recovery wrote to original + all duplicates, so the visualization split needs its
own mechanism (O3). The work was (1) explicit per-interface-class treatment in `InterfaceProcessor`
(tie vs decouple), (2) replacing the deliberate MaxwellFactory ferro source-strip with the
weight-1 tie, (3) making pair registration survive the sourceless (decoupled) case, and
(4) establishing side-local recovery for the split.

## 2. Architecture: Reuse the Hanging-Dof Spine

Architecture choice: keep `InterfaceProcessor` as the single owner of interface duplication and give each
`InterfaceSet` a **treatment tag**:

- `TieWeight1` (ferro-air): duplicate, relink, `set_sources( original, 1.0 )`. The generic
  hanging-dof machinery ties φ through the solve — no IWG, no new dof types, no solver change.
  The **viz split does not come for free** (Codex audit): pre-R5 recovery wrote each side's
  contribution to original + all duplicates, so R5/O3 had to make the ferro-air recovery
  side-local (or mask the cross-side writes) for the split to appear in the output.
- `Decouple` (any coil-touching interface): duplicate, relink, **no sources**; if registration
  is wanted, set `set_original()` / `add_duplicate()` explicitly rather than deriving it from
  sources. Coil-side nodes carry no dof — coil
  blocks are excluded from computation (`cl_Topology.cpp:554`, `update_block_map` skips Coil) —
  and `init_fields()` already fills φ with `BELFEM_QUIET_NAN` for nodes outside air/buffer/ferro
  blocks (`cl_MaxwellFactory.cpp:2104-2113`), so the NaN fill is automatic once the duplicates
  live only in coil elements.

Rejected: handling air-coil inside `CutProcessor` (that machinery is cut/cohomology-specific; coil
interfaces carry no cuts) and inventing an edge-H coupling for the coil side (the coil is intentionally dead —
`cl_IWG_Maxwell.cpp:403-411` enforces this).

## 3. Gap Table

| # | State / behaviour | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| 1 | Air-ferro interface admitted to an `InterfaceSet` | ferro-air split | yes by inspection (predicate passes); runtime unverified | (b) → R1 | `cl_InterfaceProcessor.cpp:337-341` |
| 2 | Air-coil interface admitted | coil split | textually yes, but coil not explicit; buffer-coil skipped | (c) | `:286-293, 337-341`; Codex dl20260703 |
| 3 | Which side gets duplicated (master priority) | correctness of relink + strip exemption | **yes — resolved (O1):** `MaxwellFactory::fix_facet_masters()` flips facets so the higher `DomainType` value is master → ferro(3)/coil(4) master over air(1); `InterfaceSet` duplicates the master side | (a) | `cl_MaxwellFactory.cpp:1069-1105` (flip at `:1092-1096`), `cl_InterfaceProcessor.cpp:94-96`; stale ordering comment `cl_MaxwellFactory.hpp:184` (Codex) |
| 4 | Weight-1 hang survives MaxwellFactory | φ continuity across ferro-air | **no** — ferro strip severs it | (c) | `cl_MaxwellFactory.cpp:1892-1898` |
| 5 | HPhi source refilter keeps the air/ferro original | same | yes — air, buffer, ferro nodes are flagged before the refilter | (a) | `cl_MaxwellFactory.cpp:1810-1873` |
| 6 | Decoupled duplicate registered as viz pair | postproc split | **no** — `unify_duplicates()` discovers originals only through the sources list; and `set_original()` alone is insufficient: the original also needs `add_duplicate()` since postproc enumerates `original()->number_of_duplicates()` (Codex) | (c) → also O6 | `cl_InterfaceProcessor.cpp:503, 546-558`, `cl_Node.hpp:295`, `cl_FEM_Postprocessor.cpp:257` |
| 7 | Coil nodal fields = NaN | "super clean" dead coil | φ: yes, automatic; **all other fields default to 0.0** (`Mesh_Field` init), and registered coil duplicates would even be **overwritten by air-side recovery** (all-duplicates write, row 8) | (c) → O4, O6 | `cl_MaxwellFactory.cpp:2104-2113`, `cl_Mesh_Field.cpp:40` (Codex) |
| 8 | Postproc recovers ferro/air sides independently | discontinuous B/H | **no** — `recover_fields()` canonicalizes to `original()` and writes to original + all registered duplicates; `Gradient` likewise writes one value to all copies. The previous conductor-side split mechanism, if present, was unestablished | (c) → O3 | `cl_FEM_Postprocessor.cpp:985-1013` vs `cl_Gradient.cpp:102-149` (Codex, verified by Claude) |
| 9 | No double duplication of a node already split by `CutProcessor` | triple junctions (air-coil-conductor) | `collect_nodes_and_elements()` walks `original()` + duplicates, but no `is_duplicate()` guard on the clone step — a re-clone would trip the duplicate-of-a-duplicate assertion | (c) | `cl_InterfaceProcessor.cpp:52-74, 147-170`, `cl_Node.hpp:258-263` |
| 10 | Periodicity of duplicates | periodic models | handled — periodic partners paired and backed up | (a) | `cl_InterfaceProcessor.cpp:173-216` |
| 11 | MPI ownership of duplicates | parallel postproc | infrastructure exists — `Mesh::update_ownerships()` assigns originals + registered duplicates the min owner over all attached elements, and the distributor serializes duplicates; R7 verifies the update runs **after** the new duplicates exist | (b) → R7 | `cl_Mesh.cpp:1637-1665`, `cl_Mesh_Distributor.cpp:1035`, `cl_FEM_Postprocessor.cpp:246-264` (Codex) |

### 3.1 Cross-cutting findings

- **Recovery is not side-local (Codex, the biggest gap in the draft).** `recover_fields()` writes
  each element contribution to the original and every registered duplicate
  (`cl_FEM_Postprocessor.cpp:987, 999-1013`) — so within one domain postprocessor, originals
  and duplicates receive identical values, and across postprocessors the last writer wins at
  shared field names. The ferro-air split therefore required either side-local recovery or
  masked writes (O3); coil duplicates, if registered, would be overwritten by air recovery (O6).
- **The ferro strip (`:1892-1898`) appears to be the *previous* design, not an accident.** Git
  blame → `582131de` ("end of year push", 2024-12-13): the same change assigned ferro-air
  duplicate sources and then reset ferro-block node sources before the per-domain postprocessors —
  i.e. it deliberately *decoupled* ferro-air duplicates (Codex, confidence medium). Christian's
  weight-1-tie decision supersedes it; R4 replaces rather than patches it.
- **`unify_duplicates()` keys everything off the sources list.** The decouple path must not rely
  on it; set `set_original()` **and** `add_duplicate()` explicitly at creation time (if the pair
  is registered at all — see O6).

## 4. Ordered Steps

- [ ] ~~**R1 — Ground truth.** Add log-only instrumentation (no behaviour change) for which
  interface sets `create_interface_sets()` builds on an iron-yoke + coil model, the duplicates' hanging
  state before/after `MaxwellFactory`, and — most importantly — **how the existing conductor-side
  viz split actually appears in the output** given that recovery writes to all duplicates
  (evidence for O3). Christian runs the build/model. *(read-only + logging)*~~
  **Obsolete 2026-07-03** — ground truth obtained empirically instead (three serial runs +
  ParaView φ/B evidence); O3 was answered by implementation (no split mechanism ever existed).
- [ ] **R2 — Explicit classification** *(after: R1)*. Partially done 2026-07-03: `InterfaceSet`
  got the treatment tag {`TieWeight1`, `Decouple`}, chosen in its constructor from the block
  `domain_type()` (coil on either side ⇒ `Decouple`) — no separate coil bitset needed for that.
  Still open: decide Buffer≡Air membership (O5) and make the admission predicate explicit per
  class instead of implicit.
- [x] **R3 — Decouple path** *(done 2026-07-03, serial-validated)*. Implemented as: decoupled sets
  skip the `tOriginalMap` registration in `add_duplicates()` (⇒ the flat sourcing loop in
  `add_duplicate_nodes_to_mesh()` skips them via single-`find()` miss); coil duplicates stay
  **unregistered** (O6); `unify_duplicates()` final registration loop guarded by
  `is_duplicate()` — also fixes a latent null-container `add_duplicate()` write for any
  sourceless duplicate (`original()` returns self, container never allocated).
- [ ] **R4 — Replace the ferro strip with the tie** *(re-scoped 2026-07-03)*. The strip at
  `cl_MaxwellFactory.cpp:1882-1900` is **still active**, yet the serial run is physical with
  crisp interfaces — likely the dof-level T-matrices are built from mesh-level sources during
  kernel init *before* the strip runs, so it only clears mesh-level bookkeeping and the solve
  tie survives. Verify that ordering (`create_dofwise_t_matrices_master()` call site vs the
  strip); then decide whether to remove or keep it. The original "severs the tie → wrong physics" premise is
  not supported by the serial evidence.
- [x] **R5 — Side-local recovery for the split** *(done 2026-07-03, serial-validated)*.
  Implemented in `cl_FEM_Postprocessor.cpp` as side-local claiming/accumulation (raw
  `element->node(k)`) **plus the same-PP sibling merge**: a node's recovery patch is its own
  elements ∪ elements of registered siblings selected in the same postprocessor. Same-domain
  pairs (cohomology cuts) recover the full disc → continuous; cross-domain pairs (ferro-air)
  stay side-local → crisp. V and B use the same membership predicate (`is_flagged(0)` /
  `mNodeMatrices.key_exists`) so the projection is consistent by construction. `Gradient`
  confirmed **not** in the output path (never instantiated; air B/H comes from
  `MaxwellPostprocessor`, `cl_MaxwellFactory.cpp:1923-1927`). Coil NaN scope beyond φ = O4,
  still open.
- [ ] **R6 — Double-duplication guard** *(after: R3)*. In the clone step, guard against nodes
  already duplicated by `CutProcessor` (triple junctions); assert rather than silently re-clone.
- [ ] **R7 — Parallel checks** *(after: R3, R4)*. Verify duplicate `owner()` assignment through
  the distributor and exercise periodic pair backup (`:173-216`) on a periodic model if available.
  **Expanded 2026-07-03:** Codex audit round 2 confirmed the shared-pair-owner hole (silent 0.0
  rows on duplicates whose own-side elements live on a non-owner rank); fixed by per-side
  ownership in `Mesh::update_ownerships()` (each node owned from its own elements, duplicate
  falls back to the original's owner when its element list is empty). Verify on 2/4 procs:
  no zero stripes at interfaces/cuts, solve residual histories unchanged (edge/facet owners
  shift → different but consistent dof partition — soft change accepted by Christian).
- [ ] **R8 — End-to-end validation** *(after: all)*. Iron-yoke + coil reproducer, serial + 2/4
  procs: (i) solve identical to reference (φ tie is exact — same iteration counts and residuals),
  (ii) exodus output shows split B/H at the ferro-air interface, (iii) coil region NaN/dead,
  (iv) conductor-interface viz unchanged. Christian runs; results are recorded here + devlog.

## 5. Open Design Questions

- **O1 — Which side is duplicated?** **RESOLVED 2026-07-03 (Codex, verified by Claude) →**
  `MaxwellFactory::fix_facet_masters()` (`cl_MaxwellFactory.cpp:1069-1105`) flips facets so the
  higher numeric `DomainType` becomes master: ferro(3) and coil(4) are master over air(1); it runs
  before `CutFactory` spawns `InterfaceProcessor`, and `InterfaceSet` duplicates the master side
  (`cl_InterfaceProcessor.cpp:94-96`). Side note: stale ordering comment at
  `cl_MaxwellFactory.hpp:184`.
- **O2 — Original purpose of the ferro strip** (`cl_MaxwellFactory.cpp:1892-1898`).
  **Mostly resolved 2026-07-03 (Codex git-blame, confidence medium) →** introduced by `582131de`
  ("end of year push", 2024-12-13) together with ferro-air duplicate source assignment; it
  deliberately severed those hangs before the per-domain postprocessors — the previous *decouple*
  design for ferro-air. Christian's weight-1-tie decision supersedes it. Residual: confirm via R1
  the ordering in R4, then decide whether to remove or keep the strip.
- **O3 — How does the viz split actually work?** **RESOLVED 2026-07-03 (by implementation) →**
  it never worked: no split mechanism existed. The old recovery wrote every element contribution
  to original + all registered duplicates, so both copies carried cross-side-blended values and
  last-writer-wins masked the absence (interfaces were in fact blurred). The mechanism is now
  the same-PP sibling merge introduced by R5: postprocessor membership distinguishes gauge cuts
  (same PP → merged full-disc patch, continuous) from material interfaces (different PP →
  side-local, crisp).
- **O4 — NaN scope on the coil side.** Only φ gets NaN outside air/buffer/ferro
  (`cl_MaxwellFactory.cpp:2104-2113`); all other mesh fields default to **0.0**
  (`cl_Mesh_Field.cpp:40`), so H/B/J on coil-only nodes need explicit NaN init if the "super
  clean" option is wanted (Christian may waive it).
- **O5 — Buffer membership.** Should Buffer count as Air in the admission predicate (consistent
  with `cl_Topology.cpp:288-289`), making buffer-ferro/buffer-coil behave like air-ferro/air-coil?
- **O6 — Register coil duplicates as viz pairs at all?** **RESOLVED 2026-07-03 (Christian) →
  unregistered.** Coil duplicates carry no sources and no `set_original()`/`add_duplicate()`
  registration — recovery never touches them; `update_ownerships()` derives their owner from
  their own coil elements, NaN/init fill survives. Implemented in R3; serial-validated.

## 7. Definition-of-Done Checklist

- [ ] Every gap-table row maps to a step (rows 1→R1; 2→R2; 6→R3; 4→R4; 7,8→R5; 9→R6; 10,11→R7) or an open question (row 3 resolved: O1).
- [ ] O1-O6 resolved and annotated in place (O1/O3/O6 done; O2 mostly done).
- [ ] Solve bit-identical (residual history) on a ferro-air model before/after — the weight-1 hang must be physics-neutral.
- [ ] R8 reproducer passes serial + 2/4 procs; exodus shows the split.
- [ ] Devlog written; Codex audit recorded in §8.

## 8. Audit Trail

- Codex read-only trace: `devlog/dl20260703_interface_duplicate_node_trace.md` (high confidence) —
  located the duplication path, the Air/Ferro-only bitsets, and the `cl_IWG_Maxwell.cpp` coil-interface
  hard errors; proposed the coil-bitset extension point adopted in R2.
- Claude exploration (2026-07-03, this plan): stage B/C mechanism mapped through subagents
  (hanging-dof T-matrices, `compute_hanging_dofs`, `recover_fields` vs `Gradient`); ferro-strip
  collision and `unify_duplicates` sources-dependency found by direct read. Confidence high on
  citations, medium on O1/O2 assumptions.
- Decisions by Christian (2026-07-03): ferro-air = weight-1 hang; air-coil = full decouple,
  coil excluded from calculation, optional NaN fill.
- **Codex audit round 1 (2026-07-03):** 7 findings against the draft. Two corrections folded in:
  (i) stage-A weight-1 claim was over-generalized (hanging originals copy sources instead), (ii) the
  stage-C "independent per-side recovery" claim was wrong — recovery writes to original + all
  duplicates (new failure row + O3/O6). O1 resolved (`fix_facet_masters`), O2 traced to
  `582131de` (deliberate decouple design). Key claims (`fix_facet_masters` flip logic,
  `recover_fields` canonicalization) independently re-verified by Claude against the tree before
  inclusion. Codex ran in a read-only sandbox and could not append to the exchange file; findings
  captured from its stdout and appended to
  `tmp/ai_exchange/interface_node_duplication_coil_ferro_audit.md` by Claude.
- **Codex audit round 2 (2026-07-03, parallel ownership):** requested by Christian as second
  opinion on the side-local recovery's parallel behaviour. Verdict (high confidence, citations
  re-verified by Claude): the shared-pair-owner hole is real — distributor ghosts duplicates
  and their elements but as **aura**; FEM blocks iterate owned-only, so a duplicate whose
  own-side elements live on a non-owner rank is claimed by nobody (silent 0.0 row; collectives
  stay protocol-safe, which is what makes it silent). Preferred fix (a), per-side ownership, was
  applied to `Mesh::update_ownerships()` after Claude verified solver-safety
  (`compute_hanging_dofs()` master-only; thin-shell duplicates already per-side). Exchange:
  `tmp/ai_exchange/postproc_sidelocal_parallel_ownership.md`.
- **Implementation session (2026-07-03, Claude Opus → Fable):** R3 + R5 landed and
  serial-validated. The R5 mechanism went through three iterations (side-local → V/B-consistent
  → same-PP sibling merge) driven by Christian's ParaView evidence (`cohomology.png`,
  `FieldPhi.png`, `FieldB.png`). Full record:
  `devlog/dl20260703_interface_decouple_sidelocal_recovery.md`.
