# Devlog 2026-07-29 — Side Connectors: Fuse Fixes, SideLayer Rework, HEX8TB Element

**Date:** 2026-07-29
**Topic:** Three audit/fix rounds on the new side-connector construction in
`ThinShellFactory::create()`, the periodicity design decision for hanging outer nodes, and
the landing of the degenerated HEX8TB wall element with mesh-level (Phase 1) wiring.
**AIs involved:** Claude (primary), Codex + Grok (independent adversarial audits, every round)
**Claude Confidence:** high (all findings code-cited and 3-AI confirmed; the one open
disagreement — element handedness — is logged as O1 in the phase-2 plan, not decided)
**Literature References:** none new this session (physics consensus from 2026-07-28 carried)

## Summary

The uncoated side-edge fusing and the new SideLayer-based connector construction went
through three Claude+Codex+Grok audit rounds, alternating with Christian's fixes. All
critical defects found were fixed same-day; the degenerated HEX8TB wall element
(`ElementTemplate<8,8,4,6,1>`, 4 longitudinal edge dofs, all oriented +curve) now exists,
and connector blocks are assembled per side curve. Memory ownership of the whole
construction chain was verified clean by all three AIs. Phase 1 helper wiring (enum tables,
checkers, VTK) landed; Phase 2 (edge function, IWG, verification) is planned in
`todo/hex8tb_phase2_fem_wiring.md`.

## Key results

- **Fuse branch + subroutines hardened:** order-aware `compute_side_edge_indices`
  (stride per element order, closed-curve node-repeat assert), guards for degenerate
  curves, `original()`-normalized orientation checks. Remaining known defect: raw-index
  map keys (tracked as D1 in the phase-2 plan).
- **SideLayer rework (Christian):** edge-driven node collection with index
  remap/restore — duplicate-safe by construction because temp and layer edges reference
  original-normalized nodes (verified). Node hanging and outer-edge hanging moved into the
  ctor; outer entities hang on mid-plane sources with weights deferred to
  `DofData::create_dofwise_t_matrices_master()` (claim verified against the DofData code).
- **Periodicity design settled:** outer connector nodes own no dofs (pure hanging), so
  periodic constraints flow through their sources — the ctor periodicity block was deleted
  rather than fixed; no backup registration needed. Cross-block partner corruption risk
  (twisted helix) eliminated by not writing periodicity at all.
- **Cut/perpendicular-crossing theory:** layers (buffers included) already carry cut node
  duplicates (`create_nodes_on_layers` mirrors original/duplicate links per layer);
  `hasDuplicates`/EdgeDuplicates are the material-interface twin sheets, not cut twins.
  A cut surface exiting through a side curve mid-segment would need duplicated outer
  nodes; periodic-seam cuts absorb the problem by construction. Closed side curves are
  guarded (no terminal anchor for the binormal sign); the facet-based anchor + threefold
  demo are logged as O4.
- **HEX8TB landed (Christian):** node numbering follows plain HEX8 (deliberate); the 4 dof
  edges (0,1),(3,2),(4,5),(7,6) all point +curve. Element assembly is edge-driven with a
  per-station orientation flag — the node-container ordering (and the `unique(Cell<T*>)`
  pointer-sort trap) is no longer load-bearing. Block ownership: mesh owns block, block
  owns elements, thin shell holds non-owning refs — no leaks, no double delete (3-AI).
- **Open disagreement (O1):** element handedness for `tSign==1` — Grok leans inverted,
  Codex leans correct; node scheme is deliberate per Christian. The MeshCheckers now
  error out loudly on a negative-volume HEX8TB instead of node-swapping (a swap would
  corrupt the slot-tied edges). Numeric det(J) check gates Phase 2.

## Changes made (all uncommitted until today's commit)

- `src/fem/kernel/cl_ThinShellFactory.{hpp,cpp}` — fuse + SideLayer + connector block
  construction, guards, fixes from all three rounds (Christian + Claude).
- `src/mesh/cl_Element_HEX8TB.hpp` (new), `Mesh_Enums.hpp` (HEX8TB=125 + to_string incl.
  HEX8TS alias), `cl_Element_Factory.cpp` (Christian).
- Phase 1 helper wiring (Claude): `meshtools.cpp` (geometry_type, linear type,
  number_of_faces, number_of_nedelec_dofs, get_bottom_nodes; node/edge counts and orders
  by Christian), `cl_Mesh_CurvedElementChecker.cpp` (linear), `vtktools.cpp`
  (VTK_HEXAHEDRON + node IDs), kernel + postproc `cl_MeshChecker.{hpp,cpp}`
  (`swap_hex8tb` = loud error).
- `todo/hex8tb_phase2_fem_wiring.md` (new plan), `tmp/ai_exchange/sideconnector_wip_audit.md`
  (rounds 1–3, ephemeral).

## Open items

See `todo/hex8tb_phase2_fem_wiring.md` (R1–R7, D1–D3, O1–O5). Blocked on the edge-function
derivation (next theory session).
