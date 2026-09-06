# Periodic Ties for Conductor Volume Edges on Seam Planes

**Date:** 2026-07-10
**Purpose:** Fix the silently missing periodic constraints on bulk-conductor
edge dofs at periodic planes (helix example: J pattern rotated 90° across the
seam, r = −0.96, while net current and air field are correct).
**Module:** mesh (PeriodicityFactory), fem/maxwell (MaxwellFactory), fem/kernel (DofData)
**Status:** DONE 2026-07-10 — D1–D3 all fixed and verified. R4: cap-pair Jz correlation +0.999 (was −0.96); D3 parallel BFM-reload guard in place. R5 (corc regression watch) and R6 (parallel helix rerun) are non-blocking watch items, not open defects. O3 (λ₁/λ₀ = 0.9899 structural gap) is a separate future investigation.

Evidence and audit thread: `devlog/dl20260710_helix_seam_edge_periodicity.md`,
`tmp/ai_exchange/helix_seam_sign_flip.md`.

## Defects

- [x] **D1 (CRITICAL, found Claude/Grok 2026-07-10):** solve-phase periodic
  rebuild ties zero conductor edges. Chain:
  `create_edges(false, Conductor, {}, false)` gives cap facet wrappers no
  edges (`cl_MaxwellFactory.cpp:1147-1152`) → `collect_edges` reads
  `tFacet->element()->has_edges()` (`cl_Mesh_PeriodicityFactory.cpp:593-603`)
  → empty lists → `match_edges` early-out (`:840`) →
  `set_entity_dependencies` edge loop n=0 (`cl_Mesh_Periodicity.cpp:79`).
  Conductor cap edge dofs stay free at both planes. Confidence: high
  (statically verified + matches all field measurements).
  **Fixed 2026-07-10 (Claude, Christian-approved direction):** master-element
  fallback in `collect_edges` (`get_edges_of_facet(index_on_master)` when
  wrapper has no edges); pairing safety re-verified — `match_edges` is
  key-based on endpoint originals in the aligned backup lists (enumeration-
  order independent) and aligns intrinsic direction by node swap
  (`cl_Mesh_PeriodicityFactory.cpp:863-877`) with `compute_edge_directions()`
  rebuilt after (`:415`). Runtime verification = R4.
- [x] **D2 (LOW, found Grok 2026-07-10):** the empty edge collection is
  silent — no log, no check. **Fixed 2026-07-10:** named always-active
  `BELFEM_ERROR` in `update_periodicity` — zero collected master edges while
  any master-plane facet's master element carries edges is now fatal.

- [x] **D3 (HIGH, found Christian 2026-07-10, parallel run):** BFM reload
  aborts in `finalize_edges` → `connect_edges_to_edges` →
  `Mesh::unflag_all_edges`: the not-finalized branch walks ALL blocks ×
  elements calling `tElement->edge(e)` for the static per-type count —
  air elements (block 5) correctly have no containers →
  `mHaveEdges` assert (`cl_ElementTemplate.hpp:678`). Only the BFM path
  runs `finalize_edges` INSIDE `finalize()` (mIsFinalized still false,
  `cl_Mesh.cpp:796`); the fresh path calls it after finalize completed →
  safe mEdges branch. Helix = first BFM mesh with a partial (conductor-
  only) edge universe. **Fixed 2026-07-10 (Claude):** `has_edges()` /
  `has_faces()` guards in the element-walking branches of
  `unflag_all_edges` / `unflag_all_faces` (`cl_Mesh.cpp`); template
  override returns plain `mHaveEdges` for every concrete type, so the
  guard is total. Rest of the path audited clean:
  `connect_edges_to_elements` and `compute_edge_directions` already
  guarded; edge→element containers hold only edge-carrying elements.
  Syntax-verified with real flags; parallel rerun = R6.

## Open questions

- [x] **O1: RESOLVED 2026-07-10 → Option A (decided Christian):**
  `collect_edges` falls back to
  `tFacet->master()->get_edges_of_facet(tFacet->index_on_master())` when the
  wrapper has no edges (wrapper-first preserved for thin shells). Christian's
  caveats (edge enumeration and intrinsic direction may differ between the
  sides) verified handled: pairing is key-based, direction aligned by node
  swap at match time; second-order Nédélec explicitly unsupported
  (`match_edges` node-swap asserts < 4 nodes).
- [x] **O2: RESOLVED 2026-07-10 → periodic wins on slave rim (decided
  Christian).** Mechanics verified against the cascade code:
  - The `DofData:3495-3499` "source is hanging" assert is NODE-path only;
    the EDGE→EDGE 1:1 path has a real cascade: if the source edge dof is
    itself hanging, the slave dof expands through the source's sources with
    multiplied weights (`cl_FEM_DofMgr_DofData.cpp:3695-3707`). So
    slave-rim (periodic, w=+1) → master-rim (interface-hanging on air φ
    node differences, `:3563-3671`) flattens to φ-node sources. Supported.
  - Required guard: `create_hanging_edges_and_facets` PART 1
    (`cl_MaxwellFactory.cpp:1167+`, runs AFTER `set_entity_dependencies`)
    must SKIP edges that already carry a periodic edge-source, else
    `allocate_source_container` asserts non-empty container
    (`cl_Mesh_Basis.cpp:290-298`). Master-plane rim edges are unaffected
    (set_entity_dependencies only adds sources to SLAVE edges B, never A).
  - Constraint web is consistent: slave rim edge → master rim edge → master
    air φ nodes; slave air φ nodes → master air φ nodes (node branch).
    Cap-interior slave edges are trivial 1:1 on non-hanging master edges.
  - **Ordering hazard (feeds R2):** dof-level `is_hanging()` of the source
    is evaluated DURING the same `mHangingDOFs` loop (`:3486`); if the
    slave-rim dof is processed before the master-rim dof received its
    node-expansion, the cascade check sees a not-yet-hanging source and
    links dof-on-unresolved-dof — flatness is assumed downstream (cf. the
    NODE-path assert). R2 must guarantee master-before-slave order (or
    iterate the loop until flat).
- [ ] **O3:** free generator λ₁ settles at exactly λ₁/λ₀ = 0.9899 —
  measured UNCHANGED by the seam fix and scale-invariant across timesteps
  (0.245/0.2475 pre-fix at t=2.46s; 0.0013742/0.0013882 post-fix at
  t=0.0039s). So it is NOT a seam symptom: the ratio is structural.
  Hypotheses to check: the free generator's loop encircles a slightly
  different surface than the driven one (support/one-element-ring
  difference), or a genuine physical partition (λ₁ measures the strand at
  a different station with some return path). Low priority; revisit with
  the generator sign-debt trace.

## Plan

- [x] **R1:** decide O1 (Christian), implement the edge collection fix.
  Done 2026-07-10: master-element fallback in `collect_edges`.
- [x] **R2 (after R1):** implement the O2 policy. Done 2026-07-10:
  (a) skip-guard in `create_hanging_edges_and_facets` PART 1 — edges whose
  `source(0)` is EDGE-type (periodic tie) are flagged and skipped, so
  `allocate_source_container` never collides; node-sourced interface hangs
  are untouched. (b) cascade flatness via a deferred pass in
  `create_dofwise_t_matrices_master`: 1:1 edge-on-edge dofs whose source
  edge is mesh-hanging but dof-unresolved are collected and flattened after
  the main loop (hard error on unresolvable / deeper-than-one-level chains).
- [x] **R3 (after R2):** D2 coverage check (named, always-active). Done
  2026-07-10 in `update_periodicity`.
- [x] **R4:** rerun helix (Christian); acceptance MET, verified 2026-07-10
  on `hphi_results.e-s.00012` of the post-fix run: Jz correlation
  +0.9985..+0.9991 on all four cap pairs (was −0.96); means equal across
  each pair; max point-wise |ΔJz| ≈ 0.5-0.8% of mean (discretization
  level). NOTE: O3 ratio λ₁/λ₀ = 0.9899 is UNCHANGED by the fix and
  scale-invariant across timesteps → structural, not a seam symptom; O3
  stays open as an independent question.
- [ ] **R5:** corc regression (thin-shell periodic path must be unaffected —
  wrapper-first collection preserved; wrappers with edges never hit the
  fallback).
- [ ] **R6:** parallel helix run past the D3 abort (Christian); watch for
  further partial-edge-universe assumptions downstream of the BFM reload
  (distributor, edge-direction sync `cl_Mesh.cpp:1931` receive side is safe
  by construction — send side filters on `has_edges()`).

## Related but separate (not this plan)

- Generator sign debt: `reorient_generators` 3D −1 fudge
  (`cl_Homology.cpp:687-704`, "must be another sign error elsewhere");
  input≡output-under-periodicity not detected by the literal chain-equality
  check (`cl_Homology.cpp:388-391`). Did NOT cause the helix pattern defect
  (means were correct) but should be traced once seam ties exist.
- Postproc: cap-node Bz/φ correlations ~0.85 (vs +1.0 in-plane) — likely
  one-sided recovery at seam duplicates; re-measure after fix.
