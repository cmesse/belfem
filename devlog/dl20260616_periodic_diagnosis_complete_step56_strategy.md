# Devlog 2026-06-16 — Periodic thin-cut diagnosis COMPLETE; Step 5/6 strategy adapted, reviewed, and documented

**Date:** 2026-06-16
**Topic:** Closed the corc periodic thin-cut root-cause investigation (probes 4b–4j + NZ), adapted the Step 5/6 strategy to the confirmed mechanism, had Codex + Grok review it, and refactored the module doc + plan.
**Module:** homology (`cl_CutProcessor`, `cl_CutSet`, `cl_CutData`)
**AIs involved:** Claude (diagnosis, strategy, doc/plan refactor), Codex + Grok (independent hypothesis reviews and strategy review).
**Claude Confidence:** high — the conclusion is by direct, classified data, cross-confirmed by both auditors.

## Summary

The investigation converged after **three reversals**, each refuted by the next probe — exactly why every "root cause" was confirmed with direct data before any fix:
1. **flagging** (4h: 282 blocked bits, head unflagged) → refuted by 4i (flags are periodic-symmetric; the blocks are symmetric drops);
2. **orientation / head-selection** (by elimination) → refuted by the complete-map trace (779 heads correspond / 0 not) and the node-level membership match (248/248 from periodic edges);
3. settled: the asymmetry comes from **non-periodic interior cohomology edges**.

**Confirmed root cause:** the 248 periodic-pair cut-bit asymmetries are the **legitimate period-direction (z-wrapping) cut jump**, not a bug. `match_edges` matches only periodic-*facet* edges, so interior/slab edges are non-periodic by design; a period-wrapping generator's cochain runs on them and the cut genuinely separates periodic partners. Classification (all-node z + all-edge survey): **0 in-plane, 0 vertical-partner, all seam→interior** → no matching bug. Two faces: ~208 cut-set asymmetry + ~38 per-element pattern fragmentation (`flip_node_bitsets` keys CutSet membership by the element-local hex). The `CutSet` symmetry assert is therefore **too strict**.

The mesh is a full 3D slab (17004 nodes, continuous z) — the earlier "one element thick" read was an artifact of only seam-node coords; corrected.

## Strategy (Steps 5/6) — adapted + reviewed

Drafted `tmp/periodic/step5_6_adapted_strategy.md`. **Codex + Grok both** rejected my first two fix options and converged on the design:
- **B1** (global per-node pattern in `flip_node_bitsets`) FLAWED — the per-element pattern is load-bearing (`determine_cut_sets`→`flip_node_bitsets`→`relink_element` use the element-local hex).
- **B2** (CutSet-local cross-pattern pairing) unsafe — corrupts the jump signature via `set_entity_dependencies`.
- **Adopted:** a **CutProcessor preclassification pass** (after `compute_node_bitsets`, before `duplicate_nodes`) producing per-`(pair, CutSet)` verdicts; `CutSet` consumes them. Grok refinements R1–R7 (per-(pair,CutSet) union-OR verdicts, per-pattern source wiring, the **release** mis-pair is a correctness bug not just a missing assert, persisted 4b data, multiplicity, 5c downgrade, 5a stays live).

## Changes Made

- `src/homology/doc/thick_thin_cuts_and_conjugate_edges.md` — refactored the PBC section: confirmed period-direction picture, periodic-facet-vs-interior-edge table, global-vs-per-element-signature note; refuted the old in-plane-seam consequence; cut-case tables + hexagon left untouched; theory only (Codex-verified faithful).
- `todo/periodic_thin_cut_continuity_fix.md` — rewrote Step 6 (6a preclassification … 6f acceptance) to the reviewed design; downgraded 5c (R6), marked 5e done; updated the Status line + intro to the confirmed mechanism (the in-plane/flagging/orientation framings marked refuted). Codex fidelity-checked; two stale tails (lines 5, 20) fixed.
- `tmp/periodic/step5_6_adapted_strategy.md` — the reviewed strategy (scratch).
- Diagnostic probes 4b/4c/4h/4i/4j + NZ remain in source (debug-only, tracked in plan 8i).

## Next

Implement Step 5 (rebuild policy; 5a/5f front-loadable) then Step 6 (preclassification pass + verdict-driven CutSet branch + release fix). Keep separating the legitimate jump from the 4b cuts-0/3 sign defects (D4-2).

## Files Updated

- src/homology/doc/thick_thin_cuts_and_conjugate_edges.md
- todo/periodic_thin_cut_continuity_fix.md, todo/ai_exchange.md
- tmp/periodic/step5_6_adapted_strategy.md
- devlog/dl20260616_periodic_diagnosis_complete_step56_strategy.md, devlog/README.md
