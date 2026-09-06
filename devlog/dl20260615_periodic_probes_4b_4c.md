# Devlog 2026-06-15 — Periodic thin-cut probes 4b + 4c (implemented, audited, first results)

**Date:** 2026-06-15
**Topic:** Implemented Step 4 diagnostics (probes 4b and 4c) from `todo/periodic_thin_cut_continuity_fix.md` directly in source (user-approved), Codex-audited both, and read the first corc 4b results.
**Module:** homology (`cl_CutProcessor.cpp`, `cl_CutSet.cpp`)
**AIs involved:** Claude (implementation, analysis), Codex (audit of both probes, ran `mpicxx -fsyntax-only`). Grok not invoked (no disagreement).
**Claude Confidence:** high — both probes API-checked and syntax-verified; the 4b result interpretation is cross-checked (per-cut + cross-cut join).
**Codex Audit Confidence:** high.

## Summary

Two debug-only diagnostic probes landed and were audited:

- **4b** — `CutProcessor::collect_facets()`, after the cut-case admission loop: per-quotient-face sign-pair classifier (`#DIAG 4b`). Codex confirmed it compiles and is read-only; one wording slip (lower- vs higher-id side prints) corrected in a comment, no logic change.
- **4c** — `CutSet::create_duplicates()`, at the symmetry assert: surveys *every* asymmetric periodic node pair (`#DIAG 4c`) instead of dying on the first; skips that pair's duplication and continues. Codex caught a real **release-build blocker** (B1: `tBitA/tBitB` unused once `BELFEM_ASSERT` vanishes under `-DNDEBUG -Werror`) — fixed by scoping the declarations + assert inside the `#if !NDEBUG` guard. Verified clean in debug and strict release.

## Key result — corc 4b survey (out.txt, 911 face reports)

The gdb backtrace confirmed the run still dies at `cl_CutSet.cpp:103` ("Invalid periodic flagging"); 4b ran fully because `collect_facets` precedes `duplicate_nodes`.

Cross-cut join → 448 distinct quotient pairs:
- **76% (339) genuinely one-sided** — mate face conjugated in *no* cut. This is the seam-jump structure Step 6 is designed for, so the data **supports** the one-sided-duplication policy as the dominant case (my first-glance "refutes" read was wrong — UNPAIRED *is* the expected case).
- **39 within-cut sign defects** (20 `{+,+}` double-imposition + 19 `{-,-}` hole), concentrated in **cuts 0 (19) and 3 (18)** — the in-plane generators (cut 1 cleanly one-sided, cut 2 = 2). These are a **separate** CutData sign-coherence problem (D4-2), not fixed by Step 6.

Conclusion: corc has two independent issues — a healthy one-sided seam (Step 6) and a localized sign defect in the in-plane cuts (needs a CutData fix). 4c will confirm the node-level picture and tie the 39 sign-defects to specific node pairs.

## Changes Made (source — debug-only, removal tracked in plan 8i)

- `src/homology/cl_CutProcessor.cpp`: probe 4b `#DIAG 4b` block.
- `src/homology/cl_CutSet.cpp`: probe 4c `#DIAG 4c` block + debug-guarded `#include <iostream>`; original symmetry assert retained (now `tBitA == tBitB`) inside the guard.

## Docs / tracking

- `todo/periodic_thin_cut_continuity_fix.md`: 4b/4c marked done with results; **8i removal inventory updated** with exact probe locations and the assert-restore note.
- `todo/ai_exchange.md`: 4b + 4c query/verdict/resolution threads (Codex's 4c verdict relayed by Claude — its sandbox was read-only).
- Memory: `project_4b_seam_diagnostic_result.md`.

## Next

Build (`make reset && make hphirun -j 20`) and run corc to collect the `#DIAG 4c` survey; correlate one-sided node pairs against the 4b faces, and check the 39 sign-defect faces are *not* one-sided node pairs (confirming they need the CutData sign fix, not Step 6).

## Files Updated

- src/homology/cl_CutProcessor.cpp, src/homology/cl_CutSet.cpp
- todo/periodic_thin_cut_continuity_fix.md, todo/ai_exchange.md
- devlog/dl20260615_periodic_probes_4b_4c.md, devlog/README.md
