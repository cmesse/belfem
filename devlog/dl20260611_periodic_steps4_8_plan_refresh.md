# Devlog 2026-06-11 — Periodic Steps 4-8 Plan Refresh (overnight session)

**Date:** 2026-06-11 (evening) / 2026-06-12
**Topic:** Audit and refresh of Steps 4-8 in `todo/periodic_thin_cut_continuity_fix.md`; code trace for upcoming roadblocks; Codex cross-audit. No source changes (plan-audit session per task constraints).
**AIs involved:** Claude (primary), Codex (confirm/refute audit via `ask_codex.sh`)
**Claude Confidence:** high for plan verdicts and line corrections; medium-high for the new roadblocks (Codex-confirmed); medium for the Q4 topology hypothesis
**Codex Audit Confidence:** high overall; medium-low on the Q4 enrichment
**Literature References:** Alves et al. 2022 (paper6), §VII / Fig. 7 discussion — PBCs sidestepped for the Roebel case (`n×h=0` sufficed); generator representatives non-unique. No direct precedent in the BELFEM FEM corpus for cohomology cuts lying in periodic quotient boundaries.

## Summary

Verified every Steps 4-8 item and Q3/Q4 against the committed tree (HEAD
`798a421`), applied verdicts directly to the plan, traced the code paths the
remaining steps will touch, and cross-audited with Codex (claims P1-P5,
R1-R6; full thread at the tail of `todo/ai_exchange.md`). Three substantive
outcomes: (1) a **new Step-7 blocker independent of one-sided duplication**
— the final-universe rebuild collects edges over air periodic facets whose
wrapper elements have no edge containers (now plan item 5f); (2) the
one-sided hazard **extends to the edge chain** (`create_edge_map()` has no
integrity check; folded into 5b/5d); (3) one of my claims was **refuted** by
Codex and verified so: generated `cut_*` sidesets are *untyped* —
`CutFactory::mCuts` is never populated, making the `DomainType::Cut` tagging
loop a no-op (latent dead-wiring defect, recorded in 5c).

## Recap of Steps 1-3 as implemented (cross-checked against source)

- **Step 1** (seam faces survive to the cut): trim/edge-flag split
  (`mPhiBoundaries` vs `mPhiBoundariesAndPeriodic`, `cl_CutProcessor.cpp`);
  `Topology::select_sidesets()` routes periodic domain types incl. generic
  `Periodic` (1h, `cl_Topology.cpp:446`); geometric
  `tag_periodic_sidesets()` (1i, `cl_Mesh_PeriodicityFactory.cpp:1304`,
  called from `cl_MaxwellFactory.cpp:401`) — the actual corc fix, since the
  corc input declares periodicity by plane points and no sideset carries a
  periodic type.
- **Step 2** (dup↔dup pairing in `CutSet::create_duplicates()`): periodic
  pass fenced in `has_periodicity()`, flag-6 pair gating, own-original
  sources, `set_periodic(tDupA, tDupB)`, self-periodic assert; periodic node
  closure in `CutData::flag_nodes()` (2k).
- **Step 3** (deterministic pairing): `backup_node_pairs()` at
  `cl_CutFactory.cpp:135`; registration at all four duplication sites
  (CutSet `:119`, thin-shell faces `cl_CutFactory.cpp:1832`, layers
  `cl_ThinShellFactory.cpp:1161` with geometric half-thickness tolerance,
  InterfaceProcessor `:210` with membership assert); restore branch +
  post-consumption guard in `update_periodicity()`; proto hardening
  (`create_hesse()` + degeneracy guard + plane validation +
  `mNodePairsRestored`); `pair_orientation(A,B,tol)` with `mMaxTolerance`
  restore corruption assert. Exactly one post-cut `update()` remains
  (`cl_MaxwellFactory.cpp:476` ff.).
- Diagnostics in tree: only the A/B in-plane face probes
  (`cl_CutProcessor.cpp:472/:603`); the seam-element case log was removed
  after answering the sign-distribution question.

## Key Findings

1. **5f (NEW, blocks Step 7, Codex-confirmed high):** final-universe rebuild
   vs conductor-only edges — `collect_edges()`
   (`cl_Mesh_PeriodicityFactory.cpp:569-603`) iterates `tFacet->edge(e)`
   over air periodic facets; `ElementTemplate::edge()` has only a debug
   `mHaveEdges` assert (`cl_ElementTemplate.hpp:678`) before returning the
   unallocated pointer. The final rebuild has never executed on a cut
   periodic mesh, so this was unreachable until now.
2. **Edge-chain coverage gap (folded into 5b/5d, Codex-confirmed):**
   `create_edge_map()` (`:706-721`) has no integrity check; missing keys
   fail loudly (`cl_Map.hpp:223-239`) but stale-index key *collisions* are
   silent.
3. **R2 refuted/corrected:** generated `cut_*` sidesets are untyped;
   `CutFactory::mCuts` has no producer (`cl_CutFactory.hpp:53`,
   `cl_CutFactory.cpp:2700`), so `cl_MaxwellFactory.cpp:810-812` tags
   nothing. 5c reworded; dead wiring recorded.
4. **4a discharged for corc:** the run reached `cl_CutSet.cpp:103`, past
   `determine_cut_case_3d()`'s non-unit throw — corc generators are
   unit-coefficient.
5. **8h re-verified:** live periodic markers are flags 4/5
   (`ThinShellFactory::flag_periodic_nodes()`, `:2257`); the flag-1/2 copy
   blocks (`cl_CutFactory.cpp:1835-1839`, `cl_ThinShellFactory.cpp:1164-1168`)
   remain vestigial.
6. **Q4 enriched (hypothesis, medium):** the quotient's period-direction
   generator's dual cut is a transverse cross-section — possibly the seam
   plane itself, which would make the in-plane cut and Gregory's "extra
   periodic cut" the same object. Codex: consistent with the A-probe counts
   (279/54/27/101), not decisive. Needs theory review with Gregory.
7. **Q3 proposed closed** as superseded by the reproducer campaign
   (Christian decides).

## Changes Made

- `todo/periodic_thin_cut_continuity_fix.md`: Status header updated; 4a
  downgraded; 4b-4f annotated with current line refs and context; Step 5
  purpose + 5a-5d revised (edges added, severities corrected, 5c reworded
  per R2); new 5f added (5e gate now covers 5a-5d + 5f); 6b design note;
  7c prerequisite note; 8h/8i refreshed with verified refs and the actual
  diagnostics inventory; Q3/Q4 annotated. Ordering 4→5→6→7→8 confirmed.
- `todo/ai_exchange.md`: claims entry (P1-P5, R1-R6), Codex audit (via
  wrapper), resolution entry.
- No source files touched.

## Open Questions

- Step 5b policy choice (exclude cut-side facets/edges vs quotient-cut-aware
  rule) — needs Christian's decision before Step 6.
- 5f remedy choice (allocate periodic-sideset edges vs skip unallocated
  facets) — design decision.
- Q3 closure — recommended, awaiting Christian.
- Q4 — theory review with Gregory.
- Housekeeping: `todo/ai_exchange.md` is ~4650 lines; protocol says archive
  beyond ~500. Suggest archiving the pre-2026-06-11 threads to
  `todo/ai_exchange_archive_20260612.md` at the next natural break.

## Files Updated

- todo/periodic_thin_cut_continuity_fix.md
- todo/ai_exchange.md
- devlog/dl20260611_periodic_steps4_8_plan_refresh.md
- devlog/README.md
