# Devlog 2026-04-30 — Side Connector Plumbing Audit

**Date:** 2026-04-30
**Topic:** Read-only audit of side-connector partitioning and Maxwell facet plumbing
**AIs involved:** Codex
**Codex Audit Confidence:** high for source-level blockers, medium-high for workflow risk
**Literature References:** N/A

## Summary

Audited the working-tree changes since `HEAD` for side-connector partitioning in `Kernel::partition_mesh()` and Maxwell connector-to-facet wiring in `MaxwellFactory::connect_side_connectors_with_facets()`.

## Key Findings

- `ThinShellFactory::create_side_elements()` still has stale code in the `aSign < 0` branch, referencing removed `mSideConnectorFacetMap` and `aFacetIDs`; this is a compile blocker.
- The connector facet map is indexed by layer loop `j` where it should be indexed by curve-segment loop `k`, so many connector elements would be wired to the wrong generated/shell facets after the compile blocker is fixed.
- Setting `fem::Element::mMaster/mSlave` on block elements is incompatible with current `Calculator::link()`, which treats any non-null master as sideset/interface state and calls `Group::master_integration()` on a `Block`.
- The current single `mFacet` pointer is being asked to represent both a generated connector surface facet and the adjacent shell facet. Those have different contracts, so the plumbing needs either a second facet pointer or a dedicated connector-adjacency structure.
- `h_tb()` has its stiffness contribution commented out for non-thermal connector runs, which suppresses the connector's resistive curl-curl term.
- Hidden connector sidesets avoid Exodus output, but HDF5 writing and any explicit FEM/postprocess sideset selection still assume facets have masters.

## Changes Made / Proposed

- No source changes made.
- Added this devlog entry and updated `devlog/README.md`.

## Open Questions

- Decide whether connector elements should store adjacent shell facets separately from generated connector surface facets, or whether the normal-recovery path should be refactored to avoid overloading `fem::Element::mFacet`.
- Decide whether pseudo connector sidesets should be excluded from HDF5/restart output or promoted to real mesh facets with master connectivity.

## Files Updated

- devlog/dl20260430_side_connector_plumbing_audit.md
- devlog/README.md
