# Devlog 2026-08-11 — Workflow Best-Practices Document

**Date:** 2026-08-11
**Topic:** Experience report on the multi-AI development method, written for an external reader
**AIs involved:** Claude
**Claude Confidence:** high (for the dated facts; the selection of "what worked" is judgement)
**Verification:** documentation only — nothing compiled or run

## Summary

Christian asked for a shareable markdown file describing how we approach problems — not the
code — with emphasis on which policies proved successful *after* they were adopted. Written
as `doc/ai_workflow_best_practices.md` and indexed in `doc/README.md`.

The document is deliberately positioned as the third of three: `METHODOLOGY.md` maps our
terms onto published vocabulary, `doc/ai_collaboration_protocol.md` is the normative rule
set, and the new file is the experience report — each practice paired with what changed
after adoption and what it costs. It states its evidence inline (measurements, outcomes)
and cites no devlog, todo file, or tracker ID, per the documentation rule.

## Source material

The devlog corpus 2026-03 → 2026-08 (248 entries), the six campaign pages, the debt
register, the protocol, `METHODOLOGY.md`, and the commit log. Introduction dates confirmed
from git rather than from prose: protocol `2026-03-16`, `ask_codex.sh` `2026-06-11`, the
per-topic sharded exchange `2026-06-18` (protocol §10), `ask_grok.sh` `2026-07-23`,
cross-review tooling + `METHODOLOGY.md` + evidence ladder `2026-08-05`, campaign pages +
debt register + falsification battery `2026-08-06`.

## Practices selected, and the evidence behind each

Grouped as reviews and evidence (blind cross-vendor jury, pre-registration, testimony-not-
verdict, the evidence ladder, the review stop condition, calibrated confidence), division
of labour (read-only default, human adjudication of physics and design, stated auditor
scope, human-owned executable gates), memory (two-tier artifacts with distillation,
recording rejected approaches, the campaign/register compression layer, currentness sweeps,
token anchors, code→document mechanical checking, "point, don't restate"), and domain
practices that generalise (conditional literature-first, make silent failures loud, defaults
need a run gate, standalone probes).

Two negative results are stated as prominently as the positive ones: unanimity across three
models has repeatedly been refuted by a source read, and mechanical checking cannot supply a
premise that was never written down.

## Open questions

- Whether this belongs in `doc/` or should be merged into the root `METHODOLOGY.md`, which
  currently carries the published-vocabulary mapping only.
- A Codex prose pass has not been run; the standing prose rule covers `./todo/` files, but
  this file is intended for an external reader and would benefit from one.

## Files Updated

- doc/ai_workflow_best_practices.md (new)
- doc/README.md (index entry)
