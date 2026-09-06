# Devlog 2026-06-08 - Periodic Cut Claim Audit

**Date:** 2026-06-08
**Topic:** Read-only audit of specific periodic thick-to-thin cut claims C1-C6
**AIs involved:** Claude, Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high overall; medium where runtime outcome depends on facet master ownership
**Literature References:** N/A

## Summary

Audited six falsifiable claims about the periodic h-phi cohomology-cut pipeline. The current tree confirms that periodic sidesets are included in `phi_boundary_ids()` and then trimmed as hard boundaries by `CutProcessor`, before duplicate nodes are created. It also confirms that cohomology cut duplicates do not receive periodic pointers or flags 1/2, and that the live post-cut periodicity update is a coordinate-only facet/node rematch.

## Key Findings

- C1 confirmed: `Topology::select_sidesets()` classifies `AirPeriodic`, `BufferPeriodic`, and `FerroPeriodic` as phi boundaries, and `CutProcessor::collect_facets()` unflags faces on all `mBoundaries` without a periodic exception (`src/homology/cl_Topology.cpp:433-445`, `src/homology/cl_CutProcessor.cpp:454-461`).
- C1 caveat: the peeling-loop continuation claim is not supported as stated, because boundary sideset edges are flagged and the peel test skips flagged boundary edges (`src/homology/cl_CutProcessor.cpp:401-407`, `:503-507`).
- C2 confirmed: `CutSet::create_duplicates()` and `CutFactory::link_node_duplicates_and_originals()` set hanging sources and original/duplicate links, but do not propagate `periodic()` or flags 1/2 (`src/homology/cl_CutSet.cpp:59-68`, `src/homology/cl_CutFactory.cpp:2787-2791`).
- C3 confirmed conditionally: the live rebuild is `Periodicity::update()` -> `PeriodicityFactory::update_periodicity()` -> `match_nodes()`, which pairs nodes by transformed-coordinate coincidence and cannot distinguish originals from cut duplicates (`src/fem/maxwell/cl_MaxwellFactory.cpp:476-479`, `src/mesh/cl_Mesh_PeriodicityFactory.cpp:620-642`).
- C4 confirmed: stale `CutFactory::create_sidesets_*`, `create_cut_sideset_*`, and `SideSetFactory` paths have no live call sites; the live path is `CutProcessor::create_thin_cut_sidesets()` -> `CutData::add_thin_cut_sidesets_to_mesh()`.
- C5 uncertain in the absolute form: selected 3D bad coefficient patterns fail loudly, but there is no global non-unit assert after `clean()`, and unsupported coefficients map to weight zero in `CutData::weight()`.
- C6 confirmed for translational periodicity: periodic edge matching enforces same endpoint orientation before `CutData` copies same-sign coefficients to partner edges.

## Changes Made / Proposed

- No source changes.
- Proposed first repair site: `CutProcessor` boundary trimming. Periodic sidesets should not be treated as hard cut-stop faces; keep enough boundary-edge information to prevent peeling, but do not delete cut faces merely because they lie on a periodic sideset.
- Proposed follow-up repair: propagate or rebuild periodic ownership for cohomology cut duplicates so post-cut periodic constraints cannot omit or mispair duplicate nodes.

## Open Questions

- On the failing mesh, do periodic boundary facets reference originals, duplicates, or a mix after `Mesh::finalize()` copies nodes from facet masters?
- Does a non-unit coefficient ever survive `Cohomology::clean()` in a selected cut element on the failing case, or is the coefficient issue only an independent coarse-mesh failure mode?

## Files Updated

- `todo/ai_exchange.md`
- `devlog/dl20260608_periodic_cut_claim_audit.md`
- `devlog/README.md`
