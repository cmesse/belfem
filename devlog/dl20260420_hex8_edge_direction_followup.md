# Devlog 2026-04-20 - HEX8 Edge Direction Follow-Up

**Date:** 2026-04-20
**Topic:** Follow-up review after flipping `mesh::HEX8` bottom-edge directions
**AIs involved:** Codex
**Codex Audit Confidence:** medium-high (~80%)
**Literature References:** Local Nedelec implementation convention in `src/fem/interpolation/doc/nedelec.md` Section 5

## Summary

Reviewed the user's follow-up change that flips the local bottom-edge directions in `mesh::HEX8::get_nodes_of_edge()`. This makes the current `EF_HEX8` scalar factors satisfy unit local circulation on all 12 edges, but it changes the mesh-level convention for linear hexahedra.

## Key Findings

- `EF_HEX8` now has local unit circulation on all 12 `HEX8` edges with the edited bottom-edge directions.
- `HEX20`, `HEX27`, and `HEX64` still use the old bottom-edge convention, so `HEX8` now differs from the high-order hex family.
- Edge creation itself likely remains safe because edge keys are canonicalized by node index and local element directions are recomputed by comparison.
- HEX8-to-high-order conversion should be checked before accepting the convention change permanently.

## Changes Made / Proposed

- Appended the follow-up audit to `todo/ai_exchange.md`.
- No source-code edits were made.

## Open Questions

- Should the high-order hex element definitions be updated to match the new `HEX8` bottom-edge convention?

## Files Updated

- `todo/ai_exchange.md`
- `devlog/dl20260420_hex8_edge_direction_followup.md`
- `devlog/README.md`
