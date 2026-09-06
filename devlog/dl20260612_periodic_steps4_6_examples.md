# Devlog 2026-06-12 — Steps 4-6 Worked Examples (scratch drafts + Codex design review)

**Date:** 2026-06-12
**Topic:** Discussion-ready, patch-shaped examples for Steps 4-6 of `todo/periodic_thin_cut_continuity_fix.md`, written to `./tmp/periodic/` (scratch, gitignored). No source changes. Also: third-AI (Grok) usage policy recorded; conjugate-edge doc verified earlier the same day (separate thread).
**AIs involved:** Claude (drafts), Codex (design review, dry run of the 5e gate). Grok not invoked — no Claude/Codex disagreement arose (standing policy: Codex default auditor; Grok as third voice on disagreement, `GROK_EFFORT=high`, citations verified before trust).
**Claude Confidence:** medium-high for draft mechanics (API-checked); the 5b option-A sufficiency claim was deliberately put to Codex and came back scoped.
**Codex Audit Confidence:** high overall.
**Literature References:** vocabulary anchored to `src/homology/doc/thick_thin_cuts_and_conjugate_edges.md` (verified against source 2026-06-12, V1–V5 thread).

## Summary

Christian asked for concrete examples that make the *purpose* of each
remaining step visible for discussion. Four files now sit in
`./tmp/periodic/`: `step4_seam_diagnostics.cpp` (quotient-face sign-pair
verdicts, one-bit pair logger, post-rebuild coverage scan),
`step5_emission_and_rebuild.cpp` (5a emission variants, 5b policy options
A/B, 5c sideset tagging fix, 5d tier promotion, 5f air-edge guard),
`step6_one_sided_duplication.cpp` (three-way branch + cross-pattern
precheck), and a `README.md` narrative (one-sentence purpose per step, the
4→5→6 dependency logic, and the six open decisions). Codex reviewed the
drafts as a dry run of the 5e gate; its corrections were folded back in.

## Key Findings (from the Codex K-review, reconciled)

- **K2 (FLAWED→fixed):** emitting a slave-only seam face from its slave
  element must normalize node ordering via `to_master_orientation()` — the
  draft's verbatim `set_master(slave, idx, true)` would carry slave-local
  ordering (`cl_Facet.cpp:36`, `cl_Face.cpp:68`).
- **K3 (scoped):** 5b option A (exclude jump-side quotient pairs from the
  mapping universe) is airtight for node-DOF air seams — node pairs are
  restored from the backup independently of facet selection — but would drop
  periodic constraints for edge/2nd-order-face DOFs; the implementation needs
  a named non-node-DOF guard.
- **K4 (FLAWED→redesigned, verified):** in the final update,
  `reset_faces()/reset_facets()` run *before* `map_facets()`
  (`cl_Mesh_PeriodicityFactory.cpp:326` ff.), so the filter cannot key
  quotient pairs by `face->periodic()`; correct shape: centroid-pair first,
  classify, remove both copies, compact both aligned arrays in lockstep.
- **K5 (simplification adopted):** the 6b cross-pattern guard becomes a
  CutProcessor-level precheck before `duplicate_nodes()` using existing
  `mCutSets`/`node_bitset()` state — no plumbing into `CutSet`.
- **K1/K6 (SOUND):** 4b quotient-face classification is implementable
  (per-CutData map scoping noted); the 5f guard correctly asks the wrapper
  element (`Facet::element()->has_edges()`).

## Changes Made

- `tmp/periodic/`: four new files (drafts + README), Codex corrections folded.
- `todo/ai_exchange.md`: query + Codex K-review + resolution thread; also the
  earlier doc-verification thread (V1–V5) from the same day.
- No source files touched; no plan-file changes needed (the drafts implement
  the plan as refreshed on 2026-06-12).

## Open Questions (the decision list for Christian — also in tmp/periodic/README.md)

1. 5a: variant A (normalize, recommended) vs B (skip).
2. 5b: option A (quotient-pair exclusion, recommended, with the K3 scope
   guard) vs option B (jump-aware keying).
3. 5c: tag cut sidesets at emission + type-filter; retire or populate `mCuts`.
4. 5f: skip-guard (recommended for h-φ) vs allocating air seam edges.
5. 6b: precheck shape accepted? (Codex simplification adopted in draft.)
6. One-sided duplicate marker flag vs implicit
   `is_hanging() && !is_periodic()` signature.

## Files Updated

- tmp/periodic/README.md (+3 draft .cpp files)
- todo/ai_exchange.md
- devlog/dl20260612_periodic_steps4_6_examples.md
- devlog/README.md
