# Devlog 2026-04-07 — Single-Edge Lifecycle Trace for Thin-Shell Edge ID 57008

**Date:** 2026-04-07
**Topic:** Read-only lifecycle trace of the shell edge printed as `edge 57008` in the single-layer thin-shell patch test
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Traced the shell edge printed as `edge 57008` through construction, shell-element linkage, connectivity registration, reindexing, hanging-source assignment, and the final node-based edge T-matrix branch. The static code path looks structurally correct for the air/TET4 – shell/PENTA6TS – air/TET4 case.

Important correction: the debug print's `edge 57008` field is the edge **ID**, not the edge **index**. For this linear case, the printed `field 12840` implies the edge mesh index is `12840`.

The main unresolved branch-differential after the trace is `connect_edges_to_edges()` if this edge participates in curve-facet connectivity. Otherwise, the best next step remains a runtime print of the resolved sources/weights for this exact edge DOF.

## Key Findings

- The actual mesh edge is created in `ThinShellFactory::create_edges_on_layers()`, not in `EdgeFactory` or `connect_nodes_to_edges()`.
- Prototype edge creation and layer-edge cloning both satisfy the new `Vertex` semantics: `allocate_node_container(2)` is followed by exactly two distinct slot writes.
- `PENTA6TS` edge linkage is bottom-face edges first, then top-face edges; this matches both the shell local DOF order and `MaxwellFactory::hang_thinshell_edges_on_nodes_bottom/top()`.
- `connect_nodes_to_edges()` and `connect_edges_to_elements()` are unchanged vs `devel`.
- `connect_edges_to_edges()` differs from `devel`; the extra facet-derived neighbor path is only live if `edge->number_of_facets() > 0` for this edge.
- The active hanging-source path for this edge is `hang_thinshell_edges_on_nodes_bottom()`, and its node-order logic is consistent with the shell bottom-face / master-volume-facet ordering.
- The node-sourced edge branch in `create_dofwise_t_matrices_master()` is logic-identical to `devel`; only a debug print was added downstream of `set_sources()`.

## Changes Made / Proposed

- No source-code changes.
- Added this devlog and a corresponding Codex note to `todo/ai_exchange.md`.

## Open Questions

- Does edge ID `57008` have any attached facets through `CurveFactory` in the actual run?
- Because `collect_nodes()` now sorts by vertex ID, are shell edge IDs numerically stable across `ghost` and `devel` for geometry-based comparisons?
- What are the actual runtime source node IDs and weights for DOF `119663` after `create_dofwise_t_matrices_master()`?

## Files Updated

- `devlog/dl20260407_single_edge_57008_trace.md`
- `todo/ai_exchange.md`
