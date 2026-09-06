# Side-Edge Fusing: Cut-Aware Single Authority (D8) — R1–R3

**Date:** 2026-08-06
**Purpose:** Session log — archaeology, design, Codex audit, and implementation of
the single-authority fix for the fuse-to-g branch mixing (D8) behind the rim J
blobs on the fused CORC runs.
**Module:** `src/fem/kernel` (`cl_ThinShellFactory`)
**Plan:** `todo/side_edge_fusing_cut_aware_plan.md` (kept live; this log is the
session record)
**AIs:** Claude (diagnosis, design, implementation), Codex (design audit, 8
findings). Branch convention signed off by Christian (bottom sheet, fixed);
Prof. Sirous review of the cut-station flattening folds into the R5 physics gates.

## What was established (R1, read-only)

- The outer-interface anchors (`hang_thinshell_edges_on_nodes_bottom/top`,
  `cl_MaxwellFactory.cpp:1469/:1529`) contain no cut logic: their cut-awareness is
  inherited entirely from reading `aFacet->master/slave()->get_nodes_of_facet(...)`
  on volume elements the CutFactory already relinked. Ordering:
  `create_cuts()` (:444) → `create_thinshells()` (:470, fuse inside) → anchors
  (:1310). Consequence: at fuse time the anchor sources do not exist yet, but the
  same INPUT (the relinked volume element→node links) is available.
- The fuse (`connect_side_edges` / `connect_side_nodes`) sourced from the
  mid-surface temp/curve nodes, which are branch-blind: `create_temporary_edges`
  (`cl_ThinShellFactory.cpp:1695`) keys by original-normalized pairs, `unique()`s
  (collapsing twins), and rebuilds nodes from the positional container. One
  station, two branches of g → the 61/151 cross-level λ(t) mismatches (D8).
- Probe-dump census (210,066 rows, script
  `tmp/ai_exchange/side_edge_fusing_r1_census.py`): D9 = all 9 interior blocks
  missing the same six edge offsets {1, 50, 53, 102, 105, 154} (reversal-symmetric
  pairs summing to 155, seam/cut-termination stations) — hypothesis: the twin
  collapse above. D10 = the 425 |weight-sum| ≥ 2 rows are compositions through up
  to three of the FOUR generator dofs (bases 179123–179126) at ~55 seam-zone
  stations; no anomalous weights anywhere; verdict rides on the R4 gate.

## Design (R2, v2 after Codex audit)

Single authority per rim station, resolved once in `ThinShellFactory` by matching
`original()` identity directly on the volume facet nodes (cut duplicates carry
`original()` links, `cl_CutFactory.cpp:2657-2658`) — no orientation conventions, no
second cut implementation. Edges are side-resolved exactly like the anchors (the
λ jump lives between midpoints, never within one station cluster — zero flattening
for edges); only the single per-station node dof picks one branch (first facet in
traversal order — the O1 flattening point, deliberately localized).

Codex audit (8 findings, all verified against cited lines before acceptance;
reconciliation table in `tmp/ai_exchange/side_edge_fusing_r2_audit.md`):

- **D11 (refuted draft, real):** layer nodes have no `original()` links to the
  mid-surface (`new Node` + `set_index` only) — v2 dropped the layer-node hop.
- **D12 (concern, adopted):** authority side = master unless the master block is a
  conductor, slave fallback (identity matching makes the fallback orientation-free).
- **D13 (concern, verified, pre-existing):** periodicity `update()` +
  `set_entity_dependencies()` run after the fuse and reset periodic-slave node
  sources; node-fuse snapshots at seam stations can be superseded. Identical
  exposure in the pre-change code; pointer semantics are unavailable (DofData NODE
  branch asserts non-hanging sources, `cl_FEM_DofMgr_DofData.cpp:3506`).
  Documented; explicit seam-station check added to the R4 gate.
- **D14 (refuted draft, real):** the `connect_side_edges( SideLayer* )` overload is
  never called and the side-connector path hard-aborts — SideLayer left untouched,
  deferred to the HEX8TB campaign.

## Implementation (R3, uncommitted)

`src/fem/kernel/cl_ThinShellFactory.{hpp,cpp}`:

- New `compute_side_authority( aCurve, aFacets, aEdges, aEdgeIndices,
  aEdgeSources, aNodeSources )`: per curve position, the source PAIR (flat,
  aligned to the temp edge's own node order — the DofData ±1 convention is
  untouched); per station, ONE authority node (first-writer-wins over the original
  index). Includes the D9 detector (repeated temp-edge index → loud
  `message( InfoLevel::Default, ... )` warning, NOT an assert — it is expected to
  fire on the greg deck and must not block the R4 probe run) and `BELFEM_ERROR`
  guards for unmatched stations (unlinked duplicates) and both-sides-conductor.
- `connect_side_edges( Layer* )` and `connect_side_nodes` rewired onto the
  authority containers; `original()->index()` container indexing preserved.
- Facet lookup reuses the existing `create_edge_to_face_map`.

## Next

Christian builds; R4 = fused Garber setup with `BELFEM_PROBE_FUSED_ROWS=1`, rerun
the census (target: 0/151 mixed stations, was 61; pairing stays 141,901/141,901;
D9 warning count vs the 6-per-block evidence; seam-station node rows for D13).
Then R5 physics gates (Norris flat tape, fused-vs-free A/B, net-current sum rule).

## Addendum (same evening): O1 reopened — the upper and lower cuts do not match

Christian built and ran the fused Garber deck (`greg/garber2.png`): the
alternating, time-oscillating rim blobs are GONE — the D8 incoherence is fixed —
but the rim shows a coherent edge-parallel current DEPRESSION where Norris-type
edge peaking is expected. Christian's immediate question ("the cuts on the upper
and lower side don't necessarily match") was tested against the pre-fix dump by
comparing the bottom-sheet vs top-sheet ANCHOR rows per rim stack (both anchor
territory, untouched by R3): **90 of 348 full stacks carry different λ branches
above vs below** (almost all ±λ₄ = 179126), in extended stretches along the rim.
So the two-valuedness of g at the rim is segment-wide, real transport-current
physics (the MMF wrapping the tape edge), and single-branch fusing suppresses its
through-thickness transition — the flattening premise behind O1's closure was
wrong. O1 REOPENED with three options (free rims for ship / graded two-branch
authority with complementary weights, needs DofData 4-source edge conversion /
HEX8TB wall element as the proper physics); decision Christian + Prof. Sirous.
The R3 code stays — it removed the incoherence and is the substrate for either
refinement. Christian agreed on option (a) for the release; `mFuseEdges` was
switched OFF in source the same night (`cl_ThinShellFactory.hpp:276` — the tree
had carried `true` since fusing became unconditional in b42044da, which is why
Gregory's production runs were fused; free rims are the ship state again). Data +
options: `todo/side_edge_fusing_cut_aware_plan.md` §4a(e), §9 O1. Campaign
continues next session (O1 decision with Prof. Sirous, post-fix R4 probe dump).

## References

- Plan + evidence: `todo/side_edge_fusing_cut_aware_plan.md` (§4a/§4b), D8 tracker
  in `todo/hex8tb_phase2_fem_wiring.md`.
- Audit thread: `tmp/ai_exchange/side_edge_fusing_r2_audit.md` (ephemeral).
- Related: `dl20260806_matfix_controller_port_executed.md` (Addendum 2),
  `src/fem/kernel/doc/bearing_gauge_eigenmode.md` (the separately-resolved Newton
  paralysis — not this bug).
