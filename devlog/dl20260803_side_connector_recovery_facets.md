# Devlog 2026-08-03 — Side Connectors: Recovery Facets, Defect Tracker Emptied

**Date:** 2026-08-03
**Topic:** The recovery-facet design discussion, three Claude+Codex+Grok audit rounds on
Christian's new `ThinShellFactory::create_side_connectors`, and the closure of every
open defect (F1–F6, B3, D1–D6) in the side-connector construction.
**AIs involved:** Claude (primary, incl. real-compiler probes), Codex + Grok
(independent adversarial audits, every round)
**Claude Confidence:** high (every finding file:line-cited and confirmed by at least
two of three AIs; all fixes compile-verified against real build flags)
**Literature References:** none new (physics consensus 2026-07-28/31 carried)

## Summary

Christian proposed reviving the facet idea from the old HEX8TS attempt as the wall's
channel to neighbor fields; the design discussion reframed it as a **recovery facet** —
a read-only pairing object, not a physics carrier (the old wrap failed because its
facet carried penalty terms disciplining a free dof; a facet that only reads cannot
reproduce that). Christian then implemented `create_side_connectors` building one
hidden sideset per side curve: a QUAD4 facet per wall element pairing the shell
layer-block element (master, lateral face) with the HEX8TB (slave, inner face), with
the invariant **facet id = wall element id + 1**. Three audit rounds followed; the
defect tracker of `todo/hex8tb_phase2_fem_wiring.md` is now empty.

## Key results

- **Round 1 (construction audit):** four blockers, all confirmed independently by all
  three AIs — map type mismatch (compile error, verified with mpicxx), map lookup keyed
  on a layer-edge copy whose index is never set (`gNoIndex`), `create_edge_to_face_map`
  recording only edge-slot 0 of every facet (flag-all + break), and the sideset
  registered once per layer gap (m-fold delete). Plus: facet edge insertion went to the
  Vertex base container that `Facet::edge()` never reads (Codex), and the facet winding
  matched the slave face for one connector sign but not the other (Claude + Grok,
  hand-derived independently). Confirmed correct: slave face choice (2/0), master
  lateral-face convention (= base edge slot, PENTA6TS table verified), block-element
  order, ownership chain, and the id+1 invariant (held at both element and Facet level).
- **Fixes:** automated node relinking (`set_master(..., true)` — master face nodes are
  the same node objects, so linking is free and canonicalizes the winding; finalize's
  `update_facet_nodes` rewrites facet nodes from the master anyway, so construction now
  matches post-finalize state); slave orientation via `compute_orientation()` rather
  than a per-sign constant (the master's lateral-face cycle starts at a face-index-
  dependent corner — PENTA6TS faces 0/1 wind bottom-edge-first, face 2 up-leg-first —
  so the offset is mesh-dependent, not two cases); edge container moved onto the
  wrapped QUAD4; map key → temp edge; map redesigned to flag only the current curve's
  side edges with no break, arming Christian's corruption abort (two side edges on one
  facet); sideset registration moved after the gap loop.
- **Round 2/3 (D-closures):**
  - **D4/D5 (Christian's design):** the factory no longer writes `physical_tag` (the
    material machinery owns it); `EF_HEX8TB::link()` recovers the layer thickness via
    `mesh->facet( id+1 )->master()->block_id()`. Since the wall spans exactly one layer
    block and the recovery facet's master IS that block, the gap thickness is
    unambiguous — D5 never existed in this frame.
  - **Lifecycle verified:** `unfinalize()` clears `mFacets`/`mElements`; the
    MaxwellFactory re-finalize after the ThinShellFactory recollects from ALL sidesets
    (hidden included) and all blocks — so the facet map resolves at assembly time,
    and **D3 closed for free** (HEX8TB elements do enter `mMesh->elements()`).
  - **B3 (found independently by both auditors):** `create()` reset the facets' edge
    containers before the connector branch — the reset now sits just before the
    temp-edge deletion, preserving its dangling-pointer purpose.
  - **`add_sideset` hardening (Christian's rule of symmetry):** it now registers the
    sideset's facets in `mFacetMap` (map-only — pushing into `mFacets` would suppress
    finalize's collect guard), so facet lookups work independent of the finalize cycle.
  - **O5 closed by analysis (both auditors failed to refute):** the standard
    `compute_edge_directions()` machinery covers HEX8TB — the factory inserts edges
    whose node objects are the element's node-table entries, so forward stations
    resolve s=+1 and reversed s=−1, exactly the global-dof→local-ξ conversion the edge
    function needs. Guarded by the MeshChecker no-swap policy.
  - **D6 (Christian's rule):** side-connector elements inherit ownership from the
    master of their recovery facet; implemented in `Kernel::partition_mesh` after the
    layer-block ownership pass, facet follows the same owner. The recovery facets are
    not in the metis facet graph (built from `tShell->facets()` only), so the
    min-consistency sweep cannot interfere.
  - **D1:** `compute_side_edge_indices` keys now original-normalized on both ends of
    the stride walk — cut-duplicate side-curve stations resolve.

## Changes made

- `src/fem/kernel/cl_ThinShellFactory.{hpp,cpp}` — recovery-facet block rework
  (relink + orientation + element-level edges), per-curve `create_edge_to_face_map`
  with corruption abort, map key fix, sideset registration once per curve, facet
  edge-container reset moved after the connector branch, D1 key normalization,
  physical-tag write removed (Christian + Claude).
- `src/fem/interpolation/nedelec/cl_EF_HEX8TB.cpp` — thickness recovery through the
  recovery facet; physical tag no longer consumed.
- `src/mesh/cl_Mesh.cpp` — `add_sideset` registers facets in `mFacetMap`.
- `src/fem/kernel/cl_FEM_Kernel.cpp` — connector ownership inheritance in
  `partition_mesh` + debug owner asserts for connector blocks.
- `src/fem/maxwell/doc/side_connector_wall_element.md` — recovery-facet subsection in
  §3.2, thickness-recovery route in §4.1, resolved caveats in §6.
- `todo/hex8tb_phase2_fem_wiring.md` — D1–D6 ticked, O5 struck, status updated.
- `tmp/ai_exchange/side_connector_recovery_facets.md` — full three-round audit record
  (ephemeral).

## Open items

Defect tracker empty. Before the wall-term (r′) theory session: Christian's numeric
smoke run of the construction (`mCreateSideConnectors = true` up to the WIP stop:
MeshChecker volumes on both signs, `compute_orientation` values, facet resolution).
Watch items: lateral `get_edges_of_facet` unimplemented on HEX8TB (faces 4/5 only) and
PENTA6TS (faces 3/4 only); `mConnectorWidth` still hardcoded; R4 geometry-interpolation
factory case pending.
