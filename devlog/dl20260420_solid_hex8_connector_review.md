# Devlog 2026-04-20 - Solid HEX8 Connector Review

**Date:** 2026-04-20
**Topic:** Read-only review of solid `HEX8` side connector rewrite
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** Local Nedelec implementation convention in `src/fem/interpolation/doc/nedelec.md` Section 5

## Summary

Reviewed the rewritten `ThinShellFactory::create_side_connectors()` path that now builds solid `HEX8` side-connector elements instead of reduced `HEX8TS` elements. The new edge topology is structurally plausible, but the current tree has compile blockers and one important layer-view selection error that would still tie cap elements to buffer edge views instead of outer-conductor edge views.

## Key Findings

- `cl_ThinShellFactory.hpp` no longer declares `collect_nodes(...)`, so `cl_ThinShellFactory.cpp` fails direct syntax checking.
- `cl_EF_HEX8.cpp` has a `ddF` typo and an out-of-bounds `dF(3,8)` derivative entry.
- Current `EF_HEX8` scalar factors now satisfy unit local circulation on all 12 `HEX8` edges; remaining issues are derivative typos.
- The new solid connector still chooses the wrong inner longitudinal edge view on lower and upper cap elements at duplicate layers.
- `h_tb_t()` still samples old connector node sets and must be updated for the solid-HEX8 node ordering.

## Changes Made / Proposed

- Appended the audit to `todo/ai_exchange.md`.
- No source-code edits were made.

## Open Questions

- Whether the side connector should explicitly store role-specific inner edge pointers for the selected layers `{0,1,m-2,m-1}` to avoid reintroducing the `Edges` vs `EdgeDuplicates` ambiguity.

## Files Updated

- `todo/ai_exchange.md`
- `devlog/dl20260420_solid_hex8_connector_review.md`
- `devlog/README.md`
