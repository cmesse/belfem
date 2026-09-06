# Devlog 2026-06-08 - Periodic Thick-to-Thin Cut Audit

**Date:** 2026-06-08
**Topic:** Read-only audit of periodic topology preservation from `CutFactory::compute_cohomologies()` into `CutFactory::compute_thin_cuts_and_duplicate_interface_nodes()`
**AIs involved:** Codex
**Codex Audit Confidence:** high for code-path observations; medium for downstream physical impact without a reproducer
**Literature References:** Local homology documentation only (`src/homology/doc/README.md`, `src/homology/doc/thick_thin_cuts_and_conjugate_edges.md`)

## Summary

Traced `CutFactory::run()` through cohomology generation and thick-to-thin conversion. The cohomology path has explicit periodic topology support, but the handoff to `CutProcessor` reduces the data to edge-index coefficients in `Cochain::getSimplicesMap()`. The thin-cut generation reconstructs faces, node bitsets, duplicates, and relinking from current mesh topology rather than from a persistent periodic/topological cut descriptor.

## Key Findings

- `SimplicialComplex::create_complex(..., true)` and `Cohomology::clean()` explicitly use periodic nodes, edges, and faces during thick-cut computation.
- `CutData::collect_edges()` and `collect_coefficients()` expand a periodic cochain edge to its paired edge, but they copy the same plus/minus coefficient without storing an orientation/sign under the periodic map.
- `CutData::determine_cut_case_3d()` uses local edge directions and coefficients to choose the local thin-cut face, so local orientation is not completely discarded.
- `CutProcessor` carries thin-cut duplication as cut-membership bitsets and hanging node sources; it does not carry explicit conjugate edge/face metadata, element-side provenance beyond master/slave faces, or periodic node/edge/face equivalence classes.
- `CutFactory::link_node_duplicates_and_originals()` links cut duplicates to originals, but it does not set `periodic()` pointers for cut duplicates. The current `src/mesh/doc/periodicity.md` describes a flag-based duplicate periodicity propagation API that is not present in the current `Periodicity` class.
- The current working tree has a stray block-scope `collect_thin_cut_edges(...)` declaration embedded inside `determine_cut_case_3d()`. This is legal C++ and not a compile blocker, but it is dead/no-op copy-paste corruption that should be removed during cleanup.

## Changes Made / Proposed

- Added this read-only audit devlog.
- No source-code changes were made.

## Open Questions

- Whether the periodic counterpart coefficient should sometimes be negated depends on the exact periodic map orientation and should be verified with a minimal periodic tetra/prism mesh reproducer.
- Whether periodic boundary facets after the later mesh rebuild expose all relevant cut duplicates depends on which element side owns each boundary facet.

## Files Updated

- `devlog/dl20260608_periodic_thick_to_thin_audit.md`
- `devlog/README.md`
