# Helix "Flipped" Current Density: Missing Periodic Ties on Conductor Volume Edges

**Date:** 2026-07-10
**Purpose:** Root-cause analysis of the apparent input/output current-density flip in the helix example. Three-AI session (Claude primary with quantitative field forensics, Codex + Grok auditors). Fix implemented same day after Christian's approval — see "Applied Fix" below and `todo/periodic_volume_edge_seam_fix.md`.
**Module:** mesh (PeriodicityFactory), fem/maxwell (MaxwellFactory), fem/kernel (DofData)

## Symptom

After the cross-set InterfaceSet fix (same day), the helix runs to completion,
but the current density at the end faces "looks flipped" between input and
output in the visualization (`cmake-build-debug/helix/helix.png`).

## Quantitative forensics (from `hphi_results.e-s.00048` + `iv_results.csv`)

1. **Net current is CORRECT everywhere.** Mean Jz ≈ +3.18e5 A/m² at all
   8 conductor caps; Jz·A = 0.250 A ≈ imposed ramp I₀ = 0.2475 A. No sign
   flip of the transport current — input/output are NOT actually swapped.
2. **The J pattern violates periodicity.** Point-wise over each of the four
   identified face pairs (123/123 cap nodes coincide under the pure
   z-translation): Jz fluctuation correlation = −0.96. Dipole analysis: all
   four z=0 caps have identical local-frame fluctuation dipoles (azimuthal),
   all four z=5 caps identical (radial) — the z=5 pattern is the z=0 pattern
   **rotated by ~90°**. The solve is 4-fold symmetric but NOT
   screw-symmetric. This is what reads as "flipped" in the colormap.
3. **Air is perfectly periodic.** In-plane B correlation +1.0000 across the
   planes (1304/1304 matched air nodes); that is also why the net current is
   right — Ampère's circulation through the tied air enforces it.
4. `iv_results.csv`: λ₀ driven (ramp), λ₁ free with V₁ ≈ 0 settling ~2.5 mA
   behind I₀ (re-check this gap after the fix).

## Root Cause (Grok audit V1-V7, citations re-verified against source)

The solve-phase periodic rebuild silently ties ZERO conductor edges:

1. `CutFactory::run` deletes the cut-time edge universe after the cuts
   (`reset_edges`, `cl_CutFactory.cpp:166-167`) — intentional.
2. `MaxwellFactory::create_edges_and_faces_on_mesh` recreates edges only on
   conductor **volume** elements: `create_edges(false, Conductor, {}, false)`
   (`cl_MaxwellFactory.cpp:1147-1152`). The periodic cap **facet wrapper
   elements never receive edge containers.**
3. `PeriodicityFactory::collect_edges` collects via `tFacet->element()`
   guarded by `has_edges()` (`cl_Mesh_PeriodicityFactory.cpp:593-603`) →
   collects nothing on either plane.
4. `match_edges` early-outs on the empty list (`:840`);
   `Periodicity::set_entity_dependencies` edge loop runs zero times
   (`cl_Mesh_Periodicity.cpp:79`). The orientation asserts never execute —
   **the failure is completely silent.**
5. The DofManager is NOT at fault: EDGE→EDGE 1:1 elimination exists and
   would fire if the mesh sources existed
   (`cl_FEM_DofMgr_DofData.cpp:3673-3711`; gate = `basis_is_hanging()`).

Consequence: every conductor cap edge dof is free at both planes; the
conductor H field carries no seam constraint. Corc never exposed this —
its conductors are thin shells (no volume edges on periodic planes). The
helix is the first solid conductor crossing a periodic plane.

Codex's alternative ranking (generator orientation / `reorient_generators`
3D −1 fudge) is refuted as the *pattern* mechanism: with enforced seam ties,
J adjacent to the planes would be point-wise periodic regardless of
generator signs. The generator sign debt is real but separate
(`cl_Homology.cpp:687-704` fudge; one-vs-two generators for periodically
closed strands, `cl_Homology.cpp:388-391` literal-equality check).

## Applied Fix (Christian-approved: Option A; O2 = periodic wins on slave rim)

Christian's caveat — edge enumeration and intrinsic direction from the
master-element fallback may differ between the sides — was verified handled
before implementing: `match_edges` pairs by endpoint-original keys over the
ALIGNED backup-restored node lists (collection order irrelevant) and
physically swaps the target edge's nodes on reversed correspondence
(`cl_Mesh_PeriodicityFactory.cpp:863-877`), with `compute_edge_directions()`
rebuilt afterwards (`:415`). Second-order Nédélec is out of scope (the swap
asserts < 4 nodes).

Changes (all syntax-verified with the real per-target build flags):

1. **`collect_edges` master-element fallback**
   (`cl_Mesh_PeriodicityFactory.cpp`): when the facet wrapper carries no
   edge container (bulk conductor caps), edges are reconstructed via
   `tFacet->master()->get_edges_of_facet(index_on_master)`. Wrapper-first
   preserved — thin-shell wrappers with edges never hit the fallback.
2. **Named coverage check** (`update_periodicity`): zero collected master
   edges while a master-plane facet's master element carries edges is now a
   fatal `BELFEM_ERROR` — the silent-hole class is dead.
3. **PART 1 skip-guard** (`create_hanging_edges_and_facets`,
   `cl_MaxwellFactory.cpp`): edges already carrying a periodic edge→edge
   source (source(0) is EDGE-type) are flagged and skipped — periodic wins
   on the slave rim; interface node-sources land only on master-plane rim
   edges. Avoids the `allocate_source_container` empty-container assert.
4. **Deferred cascade flattening**
   (`DofData::create_dofwise_t_matrices_master`): 1:1 edge-on-edge dofs
   whose source edge is mesh-hanging but whose dof was not yet resolved
   (container-order hazard) are collected and flattened after the main
   loop; unresolvable or deeper-than-one-level chains are hard errors.
   This also removes a latent order-fragility of the pre-existing
   edge-on-edge cascade.

**Verified same day** on the post-fix rerun (`hphi_results.e-s.00012`):
cap-pair Jz correlations **+0.9985..+0.9991** on all four pairs (was
−0.96), means equal per pair, max point-wise |ΔJz| ≈ 0.5–0.8 % of mean.
One surprise: the λ₁/λ₀ ratio (0.9899) is UNCHANGED by the fix and
scale-invariant across timesteps → structural, not a seam symptom — kept
open as O3 in the todo. corc regression (R5) still pending.

## Session Artifacts

- Exchange: `tmp/ai_exchange/helix_seam_sign_flip.md` (measurements, Codex
  + Grok audits, resolution)
- Correction from the same session recorded there: verified input.conf
  source/target are point-entity IDs (2,3,4 = z=0 air circle; 294,282,284 =
  z=5) — the configuration is correct.

## Same-day follow-up: parallel BFM-reload abort (D3)

The parallel run (BFM reload path) aborted in `finalize_edges` →
`connect_edges_to_edges` → `Mesh::unflag_all_edges`: "Edges for element
99328 on block 5 have not been allocated". Christian's read confirmed —
block 5 is air and correctly carries no edge containers; the fallacy is in
the walker. The not-finalized branch of `unflag_all_edges`
(`cl_Mesh.cpp:1110`) iterated all blocks × elements calling
`tElement->edge(e)` for the STATIC per-type edge count (TET4 → 6) with no
`has_edges()` guard. Only the BFM path executes `finalize_edges` inside
`finalize()` (mIsFinalized still false, `cl_Mesh.cpp:796`); the fresh path
calls it after finalize completed and takes the safe mEdges branch. The
helix is the first BFM mesh with a partial (conductor-only) edge universe —
same theme as the seam fix above.

Fix: `has_edges()` / `has_faces()` guards in the element-walking branches
of `unflag_all_edges` / `unflag_all_faces`. The `ElementTemplate` override
returns plain `mHaveEdges` for every concrete type, so the guard is total.
Path audited onward: `connect_edges_to_elements` and
`compute_edge_directions` were already guarded, and edge→element
containers hold only edge-carrying elements, so `connect_edges_to_edges`
is safe by construction after the guard. Syntax-verified with real flags;
parallel rerun (R6) with Christian.
