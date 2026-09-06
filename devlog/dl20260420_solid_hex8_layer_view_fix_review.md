# Devlog 2026-04-20 - Solid HEX8 Layer View Fix Review

**Date:** 2026-04-20
**Topic:** Review of proposed side-connector layer-view fix
**AIs involved:** Claude, Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Reviewed Claude's proposed fix for the solid-`HEX8` side connector layer-view bug. The proposal preselects the conductor-side inner longitudinal edge view once per `SideLayer`, removing the Lo/Hi-dependent duplicate-edge branch in `create_side_elements()`.

## Key Findings

- The proposed per-`SideLayer` view selection is correct for the current four-layer connector selection `{0, 1, m-2, m-1}`.
- `s=1` should select `Edges` because it is the top view of the lower outer wrap block.
- `s=2` should select `EdgeDuplicates` when present because it is the bottom view of the upper outer wrap block.
- The proposal does not address separate compile blockers or stale thermal node sampling from the prior audit.

## Changes Made / Proposed

- Appended the verdict to `todo/ai_exchange.md`.
- No source-code edits were made.

## Open Questions

- Whether to add a runtime diagnostic or assert documenting the four-side-layer assumption and selected view per `SideLayer`.

## Files Updated

- `todo/ai_exchange.md`
- `devlog/dl20260420_solid_hex8_layer_view_fix_review.md`
- `devlog/README.md`
