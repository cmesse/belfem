# Periodic Thin-Cut Continuity Fix

**Date:** 2026-06-08
**Purpose:** Implementation plan for the `phi`-continuity error that appears only when cohomology thin cuts are used with periodic boundary conditions.
**Current status:** 2026-06-18 planning refactor. This is a living implementation plan; the backup source is `todo/periodic_thin_cut_continuity_fix.org.md`.

---

## 1. Current Status / Executive Summary

**✅ FUNCTIONAL BLOCKER CLEARED (2026-06-18):** with Step 9 (Option 2) implemented, the corc reproducer now runs the full pipeline — cuts → cohomology → thin-cut duplication/relink → periodic rebuild (`create_facet_map`/`match_edges`) → `set_entity_dependencies` → assembly → nonlinear solve — and produces a physical HTS J/Jc distribution (`cmake-build-debug/corc.png`). See `devlog/dl20260618_periodic_corc_solves.md`. Remaining: Step 7 correctness acceptance (seam φ-continuity up to the cohomology jump) and Step 8 probe cleanup.

**Live forward path:** Step 9. Key periodic facets/edges by an induced cut duplicate's original periodic identity. Do not synthesize image duplicates, do not create dup-to-dup periodic pairs for these induced one-bit duplicates, and do not clobber the original node's single `mPeriodic`. Value periodicity remains in the hanging/source chain.

**Done:** Steps 1-3 are implemented and statically confirmed: periodic sideset routing, deterministic periodic node-pair backup/append/restore, and removal of post-cut geometric node re-discovery. Step 4 diagnosis is complete. Step 5a slave-only seam face emission is done. Step 5f final-universe periodic edge collection is done and validated.

**Superseded:** Step 5b jump-side exclusion was the wrong framing for corc. Step 6 one-sided duplication / verdict-table preclassification is not the main fix. Old Step 4 diagnostic detours that blamed flagging, orientation, or admission are historical evidence only.

**Open:** Step 9 validation and implementation, Step 7 runtime acceptance, Step 8 cleanup, and residual guards listed below.

**Severity:** High. Periodic cases can silently get the wrong `phi` jump across the periodic interface. This is separate from the loud coarse-mesh `|coeff| > 1` thin-cut limitation.

**Affected modules:** `src/homology/`, `src/fem/maxwell/`, `src/mesh/`

**Key evidence and provenance:** `devlog/dl20260608_periodic_thick_to_thin_audit.md`, `devlog/dl20260608_periodic_cut_claim_audit.md`, `devlog/dl20260611_periodic_step3_restore_audit.md`, `devlog/dl20260611_periodic_step3_hardening_reaudit.md`, `devlog/dl20260616_periodic_5f_fix.md`, `devlog/dl20260617_periodic_facetmap_rootcause.md`, `todo/ai_exchange.md`, `tmp/periodic/step5_6_adapted_strategy.md`, `tmp/periodic/periodic_cut_dup_pairing_strategy.md`.

---

## 2. Workflow Map

The current pipeline, in execution order:

1. `MaxwellFactory::create_cuts()`
   - Verified call order: `create_cuts()` runs before `create_edges_and_faces_on_mesh()` and `create_thinshells()` in `MaxwellFactory::create_magnetic_kernel()` (`cl_MaxwellFactory.cpp:436-439`).
   - Enters the cut creation pipeline. The rank-0 periodicity object used for the cut pipeline is created later inside `CutFactory::run()` when `mMesh->set_periodicity(...)` is called (`cl_CutFactory.cpp:134`); the node-pair backup is taken immediately afterward (`:135`).

2. `CutFactory::run()`
   - Finalizes the mesh, creates edges/faces, and creates periodicity with the full mesh.
   - If periodicity exists, calls `mMesh->periodicity()->backup_node_pairs()` immediately after `mMesh->set_periodicity(...)` (`cl_CutFactory.cpp:135`).

3. Cohomology computation
   - `CutFactory::compute_cohomologies()` produces thick-cut generators on the periodic quotient.
   - The thick-cut generator is not the current root problem.

4. Thin-cut conversion and duplication
   - `CutFactory::compute_thin_cuts_and_duplicate_interface_nodes()`
   - Verified `CutProcessor` constructor sequence (`cl_CutProcessor.cpp:71-107`):
     - `collect_facets()`
     - `collect_nodes()`
     - `determine_cut_sets()` / `create_cut_sets()`
     - `compute_edge_bitsets()`
     - `compute_node_bitsets()`
     - historical/WIP `classify_periodic_pairs()` preclassification is not the live fix path (`cl_CutProcessor.cpp:88`, definition `:975`; table declared as `mPairVerdict` in `cl_CutProcessor.hpp:75-81`)
     - `create_abstract_nodes()`
     - `duplicate_nodes()`
     - `relink_elements()`
     - `collect_duplicates()`
     - `create_thin_cut_sidesets()`
   - Timing fact: `CutProcessor::collect_duplicates()` clears the per-`CutSet` duplicate maps before control returns to `CutFactory`.

5. Duplicate-pair registration and deterministic restore setup
   - Registered dup-to-dup pairs are appended while the relevant duplicate maps still exist.
   - Known registration sites: `CutSet`, thin-shell face duplicates, thin-shell layer nodes, and `InterfaceProcessor`.
   - `CutFactory::link_node_duplicates_and_originals()` is not the place where cohomology-cut periodic duplicate propagation is implemented.

6. Edge/facet rebuild before the Maxwell solve
   - `MaxwellFactory::create_edges_and_faces_on_mesh()` rebuilds conductor volume edges/faces.
   - `ThinShellFactory` may generate periodic wrapper sidesets.
   - Step 5f fixed final-universe edge collection by sourcing stale wrapper edges from the master facet when the wrapper container is empty.

7. Maxwell-pass periodic rebuild
   - `mMesh->periodicity()->update()`
   - `mMesh->periodicity()->set_entity_dependencies()`
   - Verified call site: `cl_MaxwellFactory.cpp:477-480`.
   - In `PeriodicityFactory::update_periodicity()`:
     - `map_facets()` selects aligned master/slave facet arrays.
     - Node-pair backup restore replaces post-cut geometric node matching; the live source does not use the old `flag_periodic_entities_12()` / `collect_nodes_from_flags_12()` story.
     - Edges/faces/facets are still re-derived from the selected facets.
     - `create_facet_map()` and `create_edge_map()` build keys from node-list positions.

**Where Step 9 belongs:** inside the Maxwell-pass rebuild key construction. When a periodic facet/edge corner is an induced cut duplicate (`dup = original + cutDOF`), `create_facet_map()`, `create_edge_map()`, and related match paths must key that corner by the duplicate's original periodic identity while leaving the duplicate unpaired.

**Stale documentation note:** `src/mesh/doc/periodicity.md` still documents `flag_periodic_entities_12()` / `collect_nodes_from_flags_12()`. These names do not appear in the live source path; cleanup item 8g tracks the documentation update.

---

## 3. Current Diagnosis

### What Failed

The Maxwell-pass periodic rebuild fails in `create_facet_map()` with `Node N not flagged`. Diagnostics `#DIAG 5j/5k` measured **178 distinct** offending seam nodes on periodic master facets. They are:

- thin-cut duplicates,
- `dup = original + cutDOF`,
- sourced from the abstract cohomology cut DOF(s) plus a periodic original,
- absent from restored `master_nodes` because `restore_node_pairs()` restores only explicitly paired backup entries.

Measured classification: **178 induced / 0 genuine-geometric-jump / 4 abstract cut DOFs**.

### Why It Failed

The 178 nodes are unpaired by design in the current one-bit branch (`cl_CutSet.cpp:144`, no `set_periodic`, no `add_node_pair_to_backup`). They are not in the restored periodic node list, so `create_facet_map()` cannot assign them a valid key. The integrity check is correct; weakening it would hide a real rebuild coverage defect.

### What This Is Not

- It is not a periodic matching bug. Step 4 showed periodic-edge support/sign, node flag closure, element admission, and periodic-edge head selection are symmetric.
- It is not a genuine geometric one-sided jump in corc. The abstract-aware split found **0** genuine jumps.
- It is not the old Step 5b exclusion case. Excluding these facets would drop real periodic continuity.
- It is not the main Step 6 duplication-policy problem. The duplication is correct as-is for the current corc blocker.

### Why Pairing Was Rejected

Synthesizing an image duplicate and pairing it is unsafe:

- the one-bit branch may not create an image duplicate,
- pairing the seam duplicate to the plain periodic original would overwrite the original's single `mPeriodic`,
- `set_entity_dependencies()` resets slave sources from the master, so clobbering `mPeriodic` would corrupt the T-matrix semantics.

### Why Original-Identity Keying Is Chosen

For node DOFs, the duplicate's value is already represented by its source chain:

```text
S = o_S + cutDOF
o_S = o_M through original periodicity
therefore S = o_M + cutDOF
```

The rebuild needs topological correspondence for facets/edges, not a new value constraint. Therefore key the facet/edge corner through `o_S`/`o_M` and leave value periodicity to the hanging/source chain.

### What 9a Must Validate First

Before implementing key substitution broadly, confirm that facet/face DOF pairing needs only geometric/topological facet correspondence from original-identity keying, not the duplicate node object itself. The known 2nd-order face-DOF slave path is still a TODO (`cl_FEM_DofMgr_DofData.cpp:2354`) and remains out of scope.

---

## 4. Completed Work

### Step 1 - Periodic sideset routing / cut-face trim fix

**Status:** [x] implemented / confirmed, with runtime visual check still pending.

Purpose: keep periodic seam faces available to the thin-cut pipeline while still using seam edges as peel guards. Periodic sidesets must not be ordinary `phi` boundaries: their edges protect the peel, but their faces must not be trimmed away as external `phi` boundary faces.

Implemented notes:

- [x] Split `CutProcessor` boundary roles into `mPhiBoundaries` for trimming and `mPhiBoundariesAndPeriodic` for edge flagging.
- [x] Periodic sidesets remain in the boundary-edge flagging list (`cl_CutProcessor.cpp:409`).
- [x] Periodic sidesets are excluded from the cut-face trim list (`cl_CutProcessor.cpp:460`).
- [x] `Topology::select_sidesets()` routes periodic types to `mPhiPeriodicIDs` through `phi_periodic_ids()`, and `CutFactory` passes that list to `CutProcessor`.
- [x] Generic `DomainType::Periodic` routing fixed. `ConductorPeriodic` remains deliberately excluded.
- [x] Geometric periodic-sideset tagging fixed the corc case where periodicity is declared by plane points only (`PeriodicityFactory::tag_periodic_sidesets()`, called at the end of `MaxwellFactory::create_periodic()` before `Topology::run()`).

Historical reconciliation note: the original record contains two successive corc observations that should both remain visible. Step 1h found a generic `DomainType::Periodic` routing gap where seam edges were unguarded and the **peel** ate cut faces; it also records that the trim fix was still correct but corc did not exercise the D3 trim path at that stage. Step 1i later found the actual plane-points corc gap: no sideset carried a periodic domain type, so seam sidesets were ordinary phi boundaries and the **boundary trim** killed all in-plane cut faces (**279/54/27/101 -> 0**) before `PeriodicityFactory::tag_periodic_sidesets()` fixed the classification. Treat peel-vs-trim as a staged-history detail; both implemented fixes remain part of Step 1.

Still open for acceptance:

- [ ] **1g)** Inspect the debug thin-cut mesh and confirm that the cut reaches the periodic interface.

### Step 2 - Duplicate-pair propagation where both sides are duplicated

**Status:** [x] implemented / confirmed for both-bit / registered duplicate pairs.

Purpose: when both originals in a periodic pair are duplicated by the same cut pattern, their duplicates must also become a periodic pair and be registered while the original-to-duplicate map still exists.

Implemented notes:

- [x] Propagation lives in `CutSet::create_duplicates()`, not a later `CutFactory` pass.
- [x] Periodic originals are handled in a dedicated pass using `tOrg->periodic()`.
- [x] Both duplicates are created together and connected by `set_periodic(tDupA, tDupB)`.
- [x] Registered pairs are appended through `Periodicity::add_node_pair_to_backup(tDupA, tDupB)`.
- [x] Originals' periodic pointers and flags remain unchanged; scratch flag 6 gates each periodic pair once.
- [x] Periodic node closure was added in `CutData::flag_nodes()` so the symmetry assert is meaningful.
- [x] Self-periodic periodic duplicate handling now asserts loudly.

Scope correction: this does **not** apply to the 178 Step-9 induced one-bit duplicates. They are intentionally unpaired and rely on the hanging/source chain plus original periodicity.

Residual:

- [ ] **2f RESIDUAL:** Copying master/slave flags 1/2 is likely unnecessary because the live rebuild is backup/restore, not flag-based. Do not add it speculatively; confirm whether any live path consumes flags 1/2 on duplicates.

### Step 3 - Deterministic periodic node-pair backup/append/restore

**Status:** [x] implemented / statically confirmed; runtime acceptance remains open.

Purpose: original periodic pairs and registered duplicate pairs must survive the cut pipeline deterministically. Post-cut geometric node matching is no longer allowed to rediscover coincident nodes.

Implemented notes:

- [x] `backup_node_pairs()` is called immediately after `mMesh->set_periodicity(...)` in `CutFactory::run()` (`cl_CutFactory.cpp:135`).
- [x] Pointer backup lifetime was audited: no backed-up node is deleted during the backup-to-restore window.
- [x] Thin-shell face duplicate pairs are appended.
- [x] Cut duplicate pairs are appended.
- [x] `InterfaceProcessor` duplicate pairs are appended.
- [x] Thin-shell layer node pairs are appended with a geometry-derived tolerance.
- [x] `PeriodicityFactory::update_periodicity()` restores backed-up node pairs and then re-derives only edges/faces/facets.
- [x] A second post-cut `update()` after backup consumption trips `BELFEM_ERROR` through `mNodePairsRestored`.
- [x] Proto periodicity recomputes Hesse forms, validates plane sizes, and is marked restored so it cannot later geometric-match a duplicated mesh.
- [x] `pair_orientation(..., tolerance)` replaced strict on-plane checks; `mMaxTolerance` protects restore-time corruption checks.

Runtime acceptance still open:

Re-scoped 2026-06-11: 3g/3h are by-construction properties of the final design, and 3i is enforced by existing asserts. A debug reproducer run that completes without tripping any named assert discharges 3g-3i automatically; only 3j and 1g require deliberate manual inspection.

- [ ] **3g)** Confirm the restore branch executed with a non-empty backup at update time.
- [ ] **3h)** Confirm no mixed original-to-duplicate pair entered the restored lists; facet coverage is enforced by the rebuild integrity checks.
- [ ] **3i)** Debug reproducer completes silently through edge/face/facet rebuild asserts.
- [ ] **3j)** Final `phi` field is continuous across the periodic seam except for the intended cohomology jump.

### Step 4 - Diagnostics and final diagnosis

**Status:** [x] diagnosis complete.

Final result: the seam cut-bit asymmetries are the legitimate period-direction / z-wrapping case, not a periodic matching bug. The earlier diagnostic forks are retained in [Section 6](#6-superseded--historical-paths) as historical evidence. Original three-way provenance came from Claude's diagnostic pass, Codex's thick-to-thin audit, and Codex's claim confirm/refute audit; the final status here is the post-2026-06-18 resolved view.

Key discharged findings:

- [x] corc generators are unit-coefficient for this reproducer; the coarse-mesh `|coeff| > 1` limitation is separate.
- [x] periodic cochain support/sign is symmetric by construction.
- [x] node flags are periodic-symmetric.
- [x] element admission is symmetric.
- [x] periodic-edge head selection corresponds under the periodic map: **779 correspond / 0 not**.
- [x] the remaining cut-bit asymmetry comes from **non-periodic cohomology edges** used by a period-direction generator, not from direct periodic facet-edge mismatching.
- [x] final classification for the relevant corc seam case: **0 in-plane / 0 vertical-partner edges, all seam-to-interior**.

Residual diagnostics / guards:

- [ ] **4a RESIDUAL:** Keep generator coefficient summaries as optional permanent instrumentation for other meshes; no longer a corc gate.
- [ ] **4d RESIDUAL:** After Step 9, verify no restored pair is original-to-duplicate and every selected periodic-facet node is either keyed by original identity or rejected by a named residual guard.
- [ ] **4e RESIDUAL:** In `PeriodicityFactory::select_sidesets()`, skip/reject empty sidesets; `tag_periodic_sidesets()` already has its own empty skip, but the geometric selector remains unguarded.
- [ ] **4f RESIDUAL:** Log selected source/target sidesets and facet counts in `map_facets()` to expose stale/generated periodic-plane sidesets.

### Step 5a - Slave-only seam face emission

**Status:** [x] done / Codex-audited.

`CutData::add_thin_cut_sidesets_to_mesh()` no longer dereferences `tFace->master()` unconditionally. If a seam face is slave-only (`master()==nullptr`, set by `fix_face_slaves()` on the target plane), it emits a one-sided facet from the slave element through `set_master(tFace->slave(), index_on_slave, true)`. A both-null face is guarded by `BELFEM_ERROR`.

Forward caveat: if `cut_*` sidesets are ever routed into a true two-sided `DomainType::Cut` consumer, revisit this. `Calculator::link()` currently dereferences `slave()` unconditionally (`cl_FEM_Calculator.cpp:1080`).

### Step 5f - Final-universe periodic edge collection

**Status:** [x] done / validated.

Problem: final periodic rebuild reached stale or edge-less periodic wrapper facets. A wrapper `has_edges()` skip was tried and reverted because conductor periodic wrappers also lacked wrapper edge containers; it would silently drop real conductor periodic edge constraints.

Implemented fix: `collect_edges()` uses wrapper edges first, then falls back to the facet master via `get_edges_of_facet(index_on_master())` when the wrapper container is empty. The master fallback is gated on TRI3 periodic seam facets; invalid PENTA lateral QUAD cases fail loudly.

Measured/validated:

- `#DIAG 5f`: **996/996** edge-less wrappers were TRI3 with `master_has_edges==1`, sidesets 7/8, masters on conductor blocks 1/2.
- Periodicity active with **279/54/27/101** faces.
- All **996** pass-2 wrappers were master-sourced.
- TRI3 gate never fired.
- `#DIAG 5h`: `wrapper-vs-master-node-mismatch 0/996`, so master-sourced edges matched wrapper nodes.
- The run advanced past `collect_edges`.

New blocker exposed afterward: post-cut `match_edges` asserted on **565 src-only / 580 tgt-only** seam-edge asymmetry (net +15; nodes symmetric 846/846). This is consistent with the Step 9 induced-duplicate keying problem, not caused by the 5f fix.

---

## 5. Live Forward Plan - Step 9

**Step 9 - Original-identity keying for induced cut duplicates**

**Status:** [ ] open / live forward path.

Decided 2026-06-18 as Option 2. Findings and audit are in `tmp/periodic/periodic_cut_dup_pairing_strategy.md`; Codex + Grok confirmed. Step 9 runs before Step 7 acceptance.

- [x] **9a) Validation guard — DONE 2026-06-18 (Claude trace + Codex high-confidence cross-check).** The duplicate node object is **not** needed in the matched lists. `set_entity_dependencies()`'s node path (`cl_Mesh_Periodicity.cpp:53-75`) iterates only the paired originals; the dup's node-DOF periodicity is carried by its own source container `{o_S, cutDOF}` plus `o_S → o_M`, flattened by `Mesh::expand_hanging_basis_sources()` (`cl_FEM_Kernel.cpp:103-107`). The DOF manager builds periodic DOFs from mesh-basis **source containers**, not from `master_facets()`/`Facet::periodic()` (zero such reads in `DofData`). corc is linear (node-φ + edge-H only), so the facet/face-DOF path is benign (2nd-order slave path is the known TODO `cl_FEM_DofMgr_DofData.cpp:2354`, out of scope). **Plan refinement: 9b and 9c must land together** — `set_entity_dependencies`'s edge block (`:77-106`) runs *before* the facet block (`:107-160`) and asserts `A->node(0)->periodic()==B->node(0)` (`:89-90`), which fails on an unpaired dup endpoint; fixing only 9b would just move the abort to that edge assert.
- [x] **9b) Facet keying by original identity — IMPLEMENTED 2026-06-18 (pending corc run).** The substitution is simply **`node->original()`**: `Node::original()` returns the geometric original for a cut dup (`set_original` at `cl_CutFactory.cpp:2800`) and `this` for thin-shell dups / non-dups, so it is automatically scoped to cut dups. `create_facet_map()` uses `node(k)->original()` in the integrity check (`:833/:949/:956`) and the node-position key (`:979-980` 2D, `:995` 3D); `match_edges` likewise (`:1306-1312`). Codex + Grok CONFIRMED the keying is correct for corc (set_index over `aNodes`, source↔target match preserved, no both-bit regression *provided* a dup and its original are never collected on the same side — the 9c guard catches that).
- [x] **9c) Edge keying + edge-dependency pairing — IMPLEMENTED 2026-06-18 (Codex+Grok pre-audited, Christian-approved).** `create_edge_map()`/`match_edges` already key by `node(i)->original()->index()`. Two further pieces applied (compile clean, gnu++17 debug+release):
  - **(i) `set_entity_dependencies` edge path (`cl_Mesh_Periodicity.cpp:82-90`) — pointer-based pairing.** The master/slave edge lists are collected independently per side, so positional `mMasterEdges(e) ↔ mSlaveEdges(e)` was unreliable (only the *node* and *facet* lists are positionally paired, via restore/`map_facets`; both auditors confirmed `:82-83` was the **only** active positional edge consumer). Now: `Edge * B = A->periodic()` (the authoritative partner `match_edges` set) + `BELFEM_ERROR( B != nullptr, … )` (always-active — `B` is dereferenced immediately, so a debug-only check would null-deref in release) + original-aware orientation asserts `A->node(i)->original()->periodic() == B->node(i)->original()` (backward-compatible: non-dup `original()==self`).
  - **(ii) `create_edge_map()` integrity guard — `BELFEM_ERROR` (release-active, Christian's call for peace of mind).** Degenerate (`A == B` self-loop) and collision (key already maps to a *different* edge) → loud error, index-reported (IDs are arbitrary). Should never fire on a healthy mesh (thin-shell edges are on the thin shell's own nodes; the seam uses the cut dup, not its original). *Note (auditors): this covers the **target** side; a SOURCE-side key collision is not separately checked — a post-match involution check (`E->periodic()->periodic()==E`) would cover both, deferred for now.*
  - Downstream sweep clean (both auditors): `compute_edge_directions` already original-aware; the DOF manager reads source containers, not node `periodic()`; the thin-shell map already uses `->original()`. The ≥4-node higher-order edge swap remains unimplemented (`cl_Mesh_PeriodicityFactory.cpp:~1292`) — out of scope for linear corc. (Doc cleanup: `src/mesh/doc/periodicity.md:84` still documents the abandoned positional edge pairing — 8g.)
- [ ] **9d) 3D faces.** Extend original-identity keying to seam faces if they exist.
- [ ] **9e) Genuine-jump residual guard.** If a seam duplicate's original is non-periodic, original identity cannot resolve it. This is the abstract-aware genuine-geometric-jump bucket, measured **0** in corc. Exclude that entity from periodic matching or fail with a named error; do not silently pass it through Step 9.
- [x] **9f) Re-run corc — DONE 2026-06-18.** `create_facet_map()` / `match_edges` / `set_entity_dependencies` all pass; the run reaches and runs the magnetic solve, producing a physical HTS J/Jc distribution (`cmake-build-debug/corc.png`). Proceed to Step 7 acceptance.

Out of scope for Step 9:

- image duplicate synthesis,
- dup-to-dup periodic pairs for the 178 induced one-bit duplicates,
- clobbering original `mPeriodic`,
- rotational or anti-periodic maps,
- full 2nd-order face-orientation support,
- the separate `|coeff| > 1` limitation.

---

## 6. Superseded / Historical Paths

### Step 5b - Jump-side exclusion strategy

**Status:** [ ] **SUPERSEDED** by Step 9 for corc.

Original idea: filter periodic facets/edges containing one-sided cut duplicates out of the periodic mapping universe. Draft signature:

```cpp
bool facet_is_jump_side( Facet * f )
{
    for ( uint k = 0; k < f->number_of_nodes(); ++k )
    {
        if ( f->node( k )->is_hanging() && ! f->node( k )->is_periodic() ) return true;
    }
    return false;
}
```

Why superseded: the 178 offending corc nodes are induced-periodic through sources (`original + cutDOF`), not genuine jumps. Excluding them would drop real periodic continuity.

Surviving value:

- **RESIDUAL:** for a genuine geometric jump where the duplicate's original is non-periodic, exclusion or a named error is still required. This is now Step 9e, and corc measured **0** such nodes.
- Historical K3 caveat remains valid: exclusion is unsafe for entities carrying edge or 2nd-order face DOFs unless scoped or guarded, because it can drop required periodic constraints.
- Historical K4 correction remains valid if exclusion is ever revived: filter quotient pairs after centroid pairing and compact aligned source/target arrays in lockstep; do not filter per copy.

### Step 6 - One-sided duplication / verdict-table policy

**Status:** [ ] **SUPERSEDED** as the main fix.

The Step 6 plan assumed seam asymmetry required changing duplication policy through `CutProcessor` preclassification plus one-sided/both-sided/fragmentation branching. The 2026-06-18 source analysis shows the corc blocker is an induced-periodic keying problem; the duplication is correct as-is for the live fix. Do not implement 6a-6c/6f as the Step 9 solution.

Historical context retained:

- The sign-blind `abs(tCase)` admission proposal was audited and rejected; positive-only admission is the intended single-imposition mechanism.
- The legitimate period-direction jump was confirmed: every periodic-edge pathway is symmetric, and the asymmetry enters through non-periodic interior cohomology edges used by the z-wrapping generator.
- Classification before Step 9: about **208** cut-set asymmetric pairs and about **38** fragmentation pairs were observed as part of the historical Step 6 framing.
- [x] **5e)** Design review gate done 2026-06-16: the adapted Step 5/6 strategy (`tmp/periodic/step5_6_adapted_strategy.md`) was reviewed by Codex and Grok. That review is provenance now; Step 9 supersedes the Step 5b/6 implementation path.

Step 6 item status:

- [ ] **6a SUPERSEDED:** `CutProcessor` preclassification pass is not the live fix path.
- [ ] **6b SUPERSEDED:** verdict-table-driven one-sided duplication is not the live fix path.
- [ ] **6c SUPERSEDED:** cross-pattern fragmentation wiring is not needed for Step 9.
- [ ] **6d RESIDUAL:** release correctness concern survives. In release, debug-only symmetry handling can fall into the both-dup path for asymmetric pairs; verify/fix this independently of Step 9.
- [ ] **6e RESIDUAL:** promote symmetry-policy violations and relink "Could not find duplicate" (`cl_CutProcessor.cpp:1471`, currently `BELFEM_ASSERT`) to named runtime errors or handled branches once final policy is settled.
- [ ] **6f SUPERSEDED:** one-sided safety guard tied to Step 5b exclusion is superseded; Step 9e replaces it for genuine geometric jumps.

Source-context note from this refactor: Step 6 preclassification / one-sided-branch WIP exists in the current tree (`classify_periodic_pairs()` call at `cl_CutProcessor.cpp:88`, definition at `:975`; `mPairVerdict` at `cl_CutProcessor.hpp:75-81`; verdict types at `cl_CutSet.hpp:26-41`). Because Step 6 is superseded as the main fix, treat that code as parked/WIP or cleanup material unless Step 9 validation explicitly reuses it.

### Historical Step 4 Diagnostic Forks

These probes are important evidence but should not be presented as the active fix path.

- **4b sign-pair classifier:** implemented and Codex-audited. corc result: **911 reports**, cross-cut join -> **448** quotient pairs, **76% (339)** one-sided reports, which explains why Step 6 looked dominant before Step 9 superseded it, and **39** within-cut sign defects (**20** `{+,+}` + **19** `{-,-}`) concentrated in cuts 0 (**19**) and 3 (**18**). These sign defects are a separate D4-2 / CutData sign-coherence issue, not Step 9.
- **4c CutSet symmetry survey:** implemented and audited. It logged asymmetric periodic pairs, skipped them in debug to continue surveying, and eventually stopped at the relink null assert. This proved the old single symmetry assert was too strict, but Step 9 does not require a Step 6 duplication-policy rewrite.
- **4g periodic cochain coherence:** implemented, then found tautological. `CutData::collect_coefficients()` copies the master edge's sign to `periodic()->index()` unconditionally, `collect_edges()` adds the mate to the support, and slave edges are unflagged before the raw 1-cochain is built; therefore `weight(e)==weight(mate)` holds by construction. The cochain layer is periodic-symmetric, so one-sided node bits and 4b face defects originate downstream of the cochain. Do not restore this probe as a root-cause test.
- **4h flagging hypothesis:** historical. Counts were cut 0 -> **1654/51**, cut 1 -> **606/0**, cut 2 -> **436/0**, cut 3 -> **1073/231** SET/BLOCKED, **282** blocked bits, **98/754** seam nodes unflagged. This initially suggested asymmetric flagging, but that conclusion was later refuted.
- **4i correction:** found **0** unflagged cohomology-edge endpoints with a flagged partner; **2037** dropped endpoints were all partner-unflagged (**1935** non-periodic far-ends, **102** periodic-both-unflagged). All **248** asymmetric bit-0 nodes were flagged, so flagging was not the asymmetry source. This temporarily shifted suspicion to orientation/head selection.
- **4j admission/orientation correction:** edge admission was symmetric (`asym-visited 0`), periodic-node map had **536** pairs, and periodic-edge head selection corresponded under T (**779 correspond / 0 not**). This ruled out admission and periodic-edge head selection; the final correction is the non-periodic-edge period-direction case below.
- **Final 4h/4i/4j result:** the **248** cut-bit asymmetries came from non-periodic cohomology edges. **195/248** were explained by z=0 nodes receiving extra cuts through non-periodic edges. Logs showed **128326** non-periodic vs **4051** periodic bit-set lines. The definitive classification for the relevant case was **0 in-plane / 0 vertical-partner, all seam-to-interior**.

### Sign-Defect Material

**Status:** [ ] **RESIDUAL / separate issue.**

The 4b cut-0/3 sign defects are independent of Step 9. They may matter for in-plane generators or generated `cut_*` sidesets but are not the Maxwell-pass induced-duplicate blocker.

### Step 5c / 5d Residuals

- [ ] **5c RESIDUAL:** generated `cut_*` sideset guard remains live for cut-0/3 in-plane generators and sign defects, but is moot for corc's period-direction transverse cut faces. Do not tag `cut_*` sidesets as `DomainType::Cut` until the one-sided facet emission from Step 5a is reconciled with consumers such as `Calculator::link()`.
- [ ] **5d RESIDUAL:** promote periodic-facet node coverage from a debug-only check to a named always-active check, and add equivalent edge coverage. Under Step 9, induced duplicates are keyed by original identity; anything still missing after that is a defect or a Step 9e genuine-jump residual.

---

## 7. Acceptance Checklist

### Already Discharged

- [x] Step 1 static checks: periodic seam edges are flagged before the peel loop, and cut faces on periodic sidesets are not removed by the trim loop.
- [x] Step 2 static checks: both-bit registered periodic cut duplicates pair dup-to-dup, never original-to-duplicate, and originals remain unchanged.
- [x] Step 3 static checks: backup/restore replaced post-cut geometric node matching; registered duplicate-pair append sites are complete for the current rank-0 path.
- [x] Step 4 diagnosis: seam asymmetry is legitimate period-direction / z-wrapping behavior, not periodic matching, flagging, orientation, or admission failure.
- [x] Step 5a slave-only seam face emission compiles and is audited.
- [x] Step 5f final-universe edge collection fix is validated through the periodic rebuild edge collection stage.

### Runtime Checks Still Pending

- [ ] **7a)** Build with `make reset && make hphirun -j 20`.
- [ ] **7b)** Run the translational periodic reproducer.
- [ ] **7c)** Confirm the current Step 9 diagnostics / focused trace show the expected classification: induced duplicate vs genuine geometric jump, source/target ownership, and abstract-source structure. Historical 4b/4c diagnostics may remain useful provenance but should not be required unless still present in the current tree.
- [ ] **7d / 1g)** Confirm in debug VTK output that the cut reaches the periodic seam.
- [ ] **7e)** Confirm duplicate source definitions:
  - both-bit periodic pairs have paired duplicates with matching source definitions after `set_entity_dependencies()`,
  - induced exactly-one-bit seam duplicates have `q = p + sum(I)` and are not periodically tied by accident.
- [ ] **7f)** Confirm the periodic rebuild restores every registered duplicate-to-duplicate node pair and reports no unexpected selected periodic-facet/edge node missing from the restored/keyed universe. Under Step 9, induced dups are keyed by original identity; only a genuine geometric jump is excluded or rejected by Step 9e.
- [ ] **7g / 3j)** Confirm final `phi` field is continuous across the periodic interface up to the intended cohomology jump.
- [ ] **9f)** Re-run corc after Step 9 and confirm `create_facet_map()` / `match_*` pass and the run reaches the solve.

Build note from `CLAUDE.md`: use `make reset && make hphirun -j 20`; plain `make` may not refresh static libraries.

### Manual Inspection Checks

- [ ] **9a)** Validate the original-identity keying assumption before implementing 9b-9d broadly.
- [ ] **3g)** Confirm the restore branch executed with a non-empty backup at update time.
- [ ] **3h)** Confirm mixed original-to-duplicate pairs have no path into restored node lists.
- [ ] **3i)** Confirm debug rebuild asserts complete silently after Step 9.
- [ ] **Q4)** Review Gregory's "extra periodic cut" interpretation. Current confidence medium: the in-plane seam cut and the period-direction generator may be the same quotient object; Step 9 handles the induced duplicates either way, but the generator interpretation still needs theory review.

### Residual Guards

- [ ] **2f RESIDUAL:** confirm no live consumer needs copied flags 1/2 on duplicates.
- [ ] **4e RESIDUAL:** skip/reject empty sidesets in `select_sidesets()`.
- [ ] **4f RESIDUAL:** log selected periodic sidesets and facet counts in `map_facets()`.
- [ ] **5c RESIDUAL:** handle generated `cut_*` sidesets safely without routing one-sided facets into two-sided consumers.
- [ ] **5d RESIDUAL:** add named always-active facet and edge coverage checks.
- [ ] **6d RESIDUAL:** fix or explicitly retire release-only asymmetric-duplication behavior.
- [ ] **6e RESIDUAL:** promote relink/symmetry policy failures to named runtime errors after the final policy lands.
- [ ] **9e RESIDUAL:** reject or explicitly exclude genuine geometric jump duplicates whose originals are non-periodic.

### Cleanup Checks

- [ ] Complete Step 8 cleanup only after Step 9 and Step 7 pass.
- [ ] Remove temporary diagnostics only after preserving useful evidence in devlogs or this plan.
- [ ] Update stale docs so future debugging does not rediscover the old flag-based rebuild story.

---

## 8. Cleanup / Dead Code

Do this after the functional fix lands and passes Step 7, so cleanup does not obscure the correctness diff.

- [ ] **8a)** Remove dead `CutFactory::create_cut_sideset_2d()`.
- [ ] **8b)** Remove dead `CutFactory::create_cut_sideset_3d()`.
- [ ] **8c)** Remove dead `CutFactory::create_sidesets_2d()`.
- [ ] **8d)** Remove dead `CutFactory::create_sidesets_3d()`.
- [ ] **8e)** Remove unused `SideSetFactory`, or document a current consumer if one is intentionally revived.
- [ ] **8f)** Remove the legal but dead block-scope declaration inside `determine_cut_case_3d()`.
- [ ] **8g)** Update `src/mesh/doc/periodicity.md` so it no longer claims cut-duplicate periodicity is propagated by `link_node_duplicates_and_originals()` or by `flag_periodic_entities_12()` / `collect_nodes_from_flags_12()` unless the implementation is moved there.
- [ ] **8h)** Remove or document a consumer for vestigial flag-1/2 copy blocks in both periodic duplicate paths: `cl_CutFactory.cpp:1835-1839` and `cl_ThinShellFactory.cpp:1164-1168` (refs from 2026-06-12). Live master/slave markers are flags 4/5 in `ThinShellFactory::flag_periodic_nodes()`; flag 2 uses are scratch/local.
- [ ] **8i)** Remove temporary diagnostics once Step 7 passes.
  - Newer probes in `cl_Mesh_PeriodicityFactory.cpp`: `#DIAG 5f`, `#DIAG 5g/5i`, `#DIAG 5h/5j/5k`, guarded includes `<iostream>/<set>/<map>`.
  - `check_health()` and calls at `:366/:369`, plus the parked facet-based `match_nodes_and_edges` block in `update_periodicity`, are Christian's WIP; keep or remove per his call, not as part of Step 9.
  - `cl_CutProcessor.cpp`: A/B in-plane periodic-face probes at `:472` and `:603`.
  - `cl_CutProcessor.cpp`: probe 4b `#DIAG 4b` block and `tSeamCase` survey.
  - `cl_CutProcessor.hpp` / `.cpp`: probe 4h, `mSurveyNodeBits`, and guarded debug include.
  - `cl_CutProcessor.cpp`: probe 4i.
  - `cl_CutProcessor.cpp`: probe 4j.
  - `cl_CutSet.cpp`: probe 4c and guarded debug include. Do not blindly restore the old strict symmetry assert; it is too strict for legitimate one-sided membership. Under Step 9, relax/remove it or convert it to the final named policy.
  - Probe 4g was already removed on 2026-06-16. Do not restore it: it was tautological and unsafe because `weight()` can index out of range on non-cohomology periodic edges.
  - The earlier `DIAG seam element` log was already removed.
  - Pre-existing debug couts, independent decision: `cl_CutFactory.cpp:892` (`#add sideset`) and `:2907` (`#SAVING CURVE`).
- [ ] **8j RESIDUAL / source-context cleanup candidate:** if Step 6 WIP remains in source after Step 9 (`classify_periodic_pairs()` at `cl_CutProcessor.cpp:975`, `mPairVerdict` at `cl_CutProcessor.hpp:75-81`, and one-sided duplication branch comments), either remove it or document why it is still needed. It is not the live fix path.

---

## 9. Confirmed Side Findings

### DOF Manager Handles Hanging + Periodic Nodes

`Periodicity::set_entity_dependencies()` (`cl_Mesh_Periodicity.cpp:46-72`) supports explicitly paired hanging + periodic nodes. For a master/slave pair `(A, B)`, it resets the slave's source container and then makes the slave inherit the master's sources if the master is hanging, or hang directly off the master otherwise.

Consequences:

1. For both-bit registered cut duplicates, `set_periodic(tDupA, tDupB)` in `CutSet::create_duplicates()` is necessary and correct.
2. The slave duplicate's own sources are wiped and replaced by the master's sources; this is harmless because the originals are periodically tied.
3. The slave-facet periodic path is unhandled (`cl_FEM_DofMgr_DofData.cpp:2354` TODO); second-order / face-DOF handling is out of scope.
4. The above applies only to explicitly paired duplicates. The 178 Step-9 induced duplicates are hanging but intentionally unpaired; they rely on their own source chain plus the original's periodicity, and Step 9 keys their facets by original identity.

### Translational Periodicity Sign

`CutData` copies the same plus/minus bit to the periodic partner edge. That is correct for the translational periodic maps used here. Rotational and anti-periodic maps are out of scope.

### Coarse-Mesh `|coeff| > 1` Limitation

This is a separate loud failure mode. `CutData::weight()` stores only boolean plus/minus state, so it cannot represent `|coeff| > 1`. If a non-unit coefficient survives `Cohomology::clean()`, `determine_cut_case_3d()` should throw instead of silently producing the continuity error described here. This needs a separate multi-cut or coefficient-splitting design.

### Dead Thin-Cut Sideset Path

The live thin-cut sideset path is:

```text
CutProcessor::create_thin_cut_sidesets()
-> CutData::add_thin_cut_sidesets_to_mesh()
```

Dead cut sideset creation paths are tracked in Step 8.

---

## 10. Open Questions

- [x] **Q1)** RESOLVED: periodic facets carry both originals and duplicates because facet nodes follow the relinked element via `update_facet_nodes()`. With Step 3, geometric matching no longer competes originals against twins post-cut.
- [x] **Q2)** RESOLVED: CutSet dup-to-dup pointers survive the first post-cut `reset_nodes()`; restore re-sets them authoritatively.
- [x] **Q3)** CLOSED: Step 1 alone was insufficient; Steps 1-3 are necessary but not sufficient. The residual failure is seam-specific and handled by Step 9.
- [ ] **Q4)** Is Gregory's "extra periodic cut" a separate zero-jump generator, or does the seam closure handled here make it unnecessary? Current confidence medium. On the periodic quotient, H1 gains a generator for the period-direction loop; its natural thin-cut representative is a cross-section and may be the seam plane itself. The corc A-probe counts (**279/54/27/101** in-plane faces across the four cuts) are consistent with one dominant in-plane generator. If so, the in-plane seam cut and Gregory's extra periodic cut are the same object, and Step 9 handles the induced duplicates from that generator. The theory review with Gregory remains for generator interpretation, not the fix. Literature note: Alves et al. 2022 (paper6) sidestepped PBCs for the Roebel case (`n x h = 0` sufficed; Section VII / Fig. 7 discussion) and notes generator representatives are non-unique.
- [x] **Q5)** RESOLVED: the generic hanging-source cascade is not a clean fallback. The single-source node path asserts if the source basis is itself hanging, and the multi-source path resolves only one level. Explicit pair registration remains load-bearing for registered dup-to-dup pairs. Step 9 relies on ordinary duplicate source expansion plus original periodicity, not the rejected fallback cascade.
- [x] **Q6)** RESOLVED: no code path deletes a backed-up node during the backup-to-restore window.

---

## 11. Out of Scope

- Rotational or anti-periodic maps.
- Full second-order face-orientation handling.
- The separate coarse-mesh `|coeff| > 1` limitation.
