# Devlog 2026-03-30 — High-Order Orientation Table Generator

**Date:** 2026-03-30
**Topic:** Independent brute-force generation of orientation tables for `HEX64`, `TET20`, and `TET35`
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Created a temporary Python generator under `./tmp` to derive high-order facet orientation tables from the existing linear orientation tables plus the reference parametric node coordinates. The goal is to produce an independent comparison artifact before any hardcoded tables are added to `ElementFactory`.

## Key Findings

- The generated tables have the expected sizes for the facet element types:
  - `TET20` -> `TRI10` facets -> `10 x 12`
  - `TET35` -> `TRI15` facets -> `15 x 12`
  - `HEX64` -> `QUAD16` facets -> `16 x 24`
- The generated matrices are in the same format consumed by the existing facet tests: each entry is a reference-element volume-node index for one oriented facet node.
- The generator is independent of mesh facet-connectivity code. It uses only:
  - the linear orientation tables (`TET4`, `HEX8`)
  - linear facet shape functions (`TRI3`, `QUAD4`)
  - high-order facet-node parameter coordinates (`TRI10`, `TRI15`, `QUAD16`)
  - reference parent-element node coordinates (`TET20`, `TET35`, `HEX64`)

## Changes Made / Proposed

- Added `tmp/generate_orientation_tables_highorder.py`
- Generated `tmp/orientation_tables_highorder_generated.md`
- Appended a summary for Claude/Codex coordination to `todo/ai_exchange.md`

## Open Questions

- Whether Claude’s hardcoded `ElementFactory` tables will match these generated matrices exactly, or differ by a facet-node ordering convention.
- Whether `HEX64` should follow the same row grouping style as the existing `HEX27` orientation table comments when it is eventually hardcoded.

## Files Updated

- tmp/generate_orientation_tables_highorder.py
- tmp/orientation_tables_highorder_generated.md
- todo/ai_exchange.md
