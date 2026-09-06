# CORC cap-current anomaly: focused rim-defect audit (A1/A2/B)

**Date:** 2026-08-25
**Purpose:** Static audit of the two jury rim-defect candidates for the
corc_solder periodic-plane |J| anomaly. No source edits. No executable gate
in this session — conclusions are **reviewed**.
**Module:** `src/fem/maxwell`, `src/mesh` periodicity, `src/homology`, `src/fem/kernel`
**Exchange:** `tmp/ai_exchange/corc_rim_defect_audit.md`

## Verdict (one line)

A1 is a latent `original()`-key hazard that a non-degenerate InterfaceCondAir
facet cannot fire; A2 (PART 1 hanging of untied half-cut rim edges onto air φ)
is live in the shared hanging path; B (closed tape-cap cycles injected into
the bulk suggested 1-chain) is live in both homology lineages and mixes
generators globally, with any local cut-jump change inherited from support
overlap.

## What was not done

No `BELFEM_PROBE_FUSED_ROWS` census, no cut-azimuth dump, no source edit.
Whether A2 sits at the measured gap-strip *centers* is a cut-geometry fact
the probe dump in the chat audit is written to decide.

## Pointers

The full Q1–Q5 audit with file:line evidence is the Grok session reply of
this date (not duplicated here).
