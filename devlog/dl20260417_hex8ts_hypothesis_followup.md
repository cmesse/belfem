# Devlog 2026-04-17 — HEX8TS Hypothesis Follow-Up

**Date:** 2026-04-17
**Topic:** Read-only follow-up audit of Claude's H1-H5 side-connector hypotheses
**AIs involved:** Codex
**Codex Audit Confidence:** high

## Summary

Performed a read-only follow-up audit of the current `HEX8TS` side-connector source paths after reading the active `todo/ai_exchange.md` thread. The main outcome was negative evidence against the remaining "mechanical break" hypotheses: I did not find a source-level smoking gun in the inner-edge source wiring or in the CutFactory-to-ThinShellFactory relink path. The strongest additional blind spot is instead a modeling/precondition one: side connectors are silently skipped unless the stack has at least three layers and the outer material labels match.

## Key Findings

- `collect_inner_nodes_and_edges()` plus `create_side_elements()` and `Element::compute_edge_directions()` make the inner-edge/source path structurally consistent; no obvious H1 sign bug was found in the magnetic assembly path.
- A wrong `tSign` would mainly affect connector geometry / left-right labeling and therefore thermal-face sampling in `h_tb_t()`, not select a different magnetic kernel path.
- `CutFactory::create_thin_shell_cuts()` records master/slave node pairs before thin-shell creation, and `ThinShellFactory::collect_nodes()` plus `create_nodes_on_layers()` explicitly propagate those duplicate/original relations into the extruded layer nodes; no obvious H2 cut-orphaning path was found.
- The documented buffer-cut reroute is still not implemented in the current `CutFactory` / `MaxwellFactory` code paths. The mutual-exclusion issue with side connectors is therefore a future guard requirement, not an active runtime conflict today.
- Wrap material assignment itself is deterministic (`tMaterials(0)`), but the modeling path is narrower than the documentation story implies: connectors are only created when `materials().size() >= 3` and `first == last` (case-insensitive), and the connector thickness is hard-coded to `0.1e-3` m.

## Changes Made / Proposed

- Added a new Codex audit entry to `todo/ai_exchange.md` with source citations and calibrated confidence.
- No source-code changes were made.

## Open Questions

- If the user still sees a connector-specific magnetic failure after the known postprocessor / Pipette / thermal issues, the next highest-value check is geometric rather than topological: verify that `compute_binomial_vectors()` extrudes the connector to the intended side of the tape for the exact problematic side curve.
- The current code should probably error, not silently return, when the input requests a production wrap configuration that violates the `materials().first() == materials().last()` precondition.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260417_hex8ts_hypothesis_followup.md
- devlog/README.md
