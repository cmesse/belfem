# Devlog 2026-08-30 — Explicit depth selection for the AI auditor wrappers

**Date:** 2026-08-30
**Purpose:** Record the session that made model and reasoning effort explicit, validated and
recorded for every `ask_codex.sh` / `ask_grok.sh` invocation.
**Module:** `.claude/scripts`, `scripts/cross_review.sh`, `doc/ai_collaboration_protocol.md`

## Summary

Christian proposed that Claude must select effort — and, for Codex, model — before calling another
AI. Implemented under the standing plan+audit → code+audit rule with both vendors, as
`todo/ai_wrapper_tier_selection.md` (now complete).

The wrappers previously passed neither knob. Every Codex audit this project ever dispatched
inherited `~/.codex/config.toml`; every Grok audit inherited `~/.grok/config.toml`. Neither
recorded which depth produced a finding, so a thin round-1 result was indistinguishable from a
thin deep one.

## Key Findings

- **`gpt-5.6-sol` defaults to `low` effort.** Escalating by model alone is therefore a downgrade.
  This is the strongest argument for the rule Christian asked for: the two knobs must move together
  or not at all. Observed locally on 2026-08-30; not repository-verifiable.
- **The stale `GROK_EFFORT`-hard-400 belief was memory-only.** `doc/lessons_learned_evidence.md`
  INC-487 already recorded the expiry correctly; only the extra-repo memory file was stale.
  `ask_grok.sh:28` had a different, real defect: it advertised `none|minimal|max`, which the
  current model does not offer.
- **An unexported assignment never reaches the child.** The first draft pinned `--quick` with a
  plain assignment in `cross_review.sh`, which the wrappers would never have seen — the unattended
  post-commit hook would have kept running at the interactive default forever. Found by Grok in the
  plan round, confirmed by experiment, fixed with a prefix assignment onto `run_auditor`.
- **Presence is not validity.** Found by Codex in the code round: the driver checked that the four
  variables were non-empty but not that they were legal, so a typo failed one leg, let the other
  bill, and exited 0 with a synthetic failure entry that read like a vendor outage.

## Changes Made

- `ask_codex.sh`: `CODEX_MODEL` / `CODEX_EFFORT`, allowlist-validated before stdin is read, passed
  as `-m` and `-c model_reasoning_effort=` ahead of the positional `-`.
- `ask_grok.sh`: `GROK_MODEL` / `GROK_EFFORT`, validated, with `--model` inside `run_grok` so the
  `--resume` salvage pass keeps the depth.
- Both: `(model=…, effort=…)` on the exchange header, per-knob `[defaulted]`.
- `cross_review.sh`: per-mode depth, value validation before dispatch, prefix-assignment onto the
  child, depth on the failure entry.
- `doc/ai_collaboration_protocol.md` §9.1 (the depth table), `.claude/commands/cross-review.md`,
  `CLAUDE.md`, and the permission allowlist.

## Verification

The R10 gate ran (a)-(j) green. Two checks carry the weight: an 821 036-byte prompt through the
newly-flagged argv path, 6.4× the historical E2BIG breaking point; and `--quick` stamping
`luna`/`medium` while the invoking shell exported `sol`/`xhigh`, which is what proves the
child-environment fix. `codex exec` accepting `xhigh` and `grok --model` are now measured rather
than read off a cache. `check_doc_claims.py`: 37/37.

Two jury rounds: plan (Codex REVISE, Grok APPROVE-WITH-CORRECTIONS) and code (ten findings, all
ten confirmed against the tree and fixed). In the code round Codex found the finding that costs
money and Grok found the three that make the policy toothless; only two of ten were found by both.

## Open Questions

- The `sol` rung is allowlisted but named by no row of the table, pending an A/B against
  `terra/xhigh` on a real BELFEM diff. Nothing in this session settles which is better.
- The wrapper defaults were chosen to match the vendor configs on this date. They do not track
  drift: if either config moves, the two diverge silently.
- Whether the permission matcher treats a `VAR=value` prefix as part of the command string was not
  verified. Extra rules were added defensively rather than from a spec.

## Files Updated

`.claude/scripts/ask_codex.sh`, `.claude/scripts/ask_grok.sh`, `scripts/cross_review.sh`,
`.claude/commands/cross-review.md`, `doc/ai_collaboration_protocol.md`, `CLAUDE.md`,
`.claude/settings.local.json`, `todo/ai_wrapper_tier_selection.md`, `todo/README.md`.
Exchange threads: `tmp/ai_exchange/ai_wrapper_tier_selection.md`,
`tmp/ai_exchange/review_ai_wrapper_tier_code.md`.
