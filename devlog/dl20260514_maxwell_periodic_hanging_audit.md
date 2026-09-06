# Devlog 2026-05-14 - Maxwell Periodic Hanging Audit

**Date:** 2026-05-14
**Topic:** Audit whether MaxwellFactory can create exactly-one-hanging periodic entities
**AIs involved:** Codex
**Codex Audit Confidence:** medium-high for core Maxwell h-phi/thin-shell edge path; medium for global node/side-connector scope

## Summary

Read-only audit of MaxwellFactory periodic/hanging sequencing and entity source creation.
For the core MaxwellFactory h-phi interface and thin-shell edge paths, a periodic pair where only one entity is hanging appears to indicate non-periodic or inconsistent geometry/topology rather than a valid case.

## Key Findings

- `MaxwellFactory::create_magnetic_kernel()` creates edges/faces, creates thin shells, then updates periodic node/edge/face pairs before `create_hanging_edges_and_facets()` runs.
- `create_hanging_edges_and_facets()` assigns hanging sources to interface/thin-shell edges by processing all relevant sidesets/facets. It does not explicitly mirror to periodic counterparts, so symmetry depends on the corresponding periodic interface/thin-shell facet existing and being classified consistently.
- `Periodicity::match_edges()` asserts matched master/slave counts for selected conductor/thin-shell edges, but it does not verify that the later hanging/source state is symmetric.
- Linear MaxwellFactory face handling does not create hanging periodic faces; the higher-order facet-source block is unreachable after the linear-element guard.
- Node-level hanging can be created in ThinShellFactory/side-connector paths, so the global DofData guard is fully proven only if those paths also preserve periodic symmetry.

## Changes Made / Proposed

- No source changes made.
- Added this devlog entry and updated `devlog/README.md`.

## Open Questions

- If side connectors remain in scope, separately audit whether their node-source paths can create exactly-one-hanging periodic node pairs.
