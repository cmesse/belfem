# Codex depth tiers moved to gpt-6.1-terra and gpt-6-luna

**Date:** 2026-11-15
**Purpose:** Christian asked whether `ask_codex.sh` was current, given that `gpt-6.1-terra`,
`gpt-6.1-sol` and `gpt-6-astra` are now available. They were not: every live tier still named a model
of the 5.6 series. Bring the wrappers, the depth table and the recipes onto the new ids.
**Module:** meta (jury tooling)
**AIs involved:** Claude Opus 5.5 (primary); Codex `gpt-6.1-terra`/high and Grok `grok-4.7`/high
(code jury)
**Status:** COMPLETE — the tier change and the jury's text fixes landed; the checker gap is filed
as DR-205
**Verification:** `check_practice_drift.py` (DRIFT = the four pre-existing P3 index gaps only; no
P2 drift), `check_practice_drift.py --self-test` passed, `check_doc_claims.py` 38/38, `bash -n`
on both scripts, and the wrapper rejects an unknown id with an `Allowed:` line equal to its case
arm. `gpt-6.1-terra` was accepted in two live jury calls on 2026-11-15. `gpt-6-luna` and
`gpt-6.1-sol` are present in `~/.codex/models_cache.json` but have not yet been called live.
**Exchange:** `tmp/ai_exchange/review_codex_model_ids.md`

## Summary

Christian's rulings (2026-11-15, by quiz): the terra rows and the wrapper default move to
`gpt-6.1-terra`; the luna rows (ordinary prose sweep, unattended `--quick`) move to `gpt-6-luna`;
`gpt-6.1-sol` replaces `gpt-5.6-sol` as the escape hatch that no row names; `gpt-6-astra` stays on
the critical-framework row. The `gpt-5.6` ids stay allowlisted so earlier rounds can be
reproduced, and they go on the anchored by-design line of §9.1. Without that line, the recent
devlog headers that name `gpt-5.6-terra` would raise P2a: the checker excuses older groks in
records automatically, but not older Codex ids.

## What landed

- `.claude/scripts/ask_codex.sh`: case arm and `Allowed:` message have the same ten ids
  (P2e). The default is `gpt-6.1-terra`. The header comment now says the model default follows the
  table and only the effort default is the 2026-08-30 config observation.
- `scripts/cross_review.sh`: the same arm; `--quick` is pinned to `gpt-6-luna`/medium; the
  example in the missing-depth message is updated.
- `doc/ai_collaboration_protocol.md` §9.1: table, recipes, by-design line, the defaults
  paragraph. `.claude/commands/cross-review.md` and the `CLAUDE.md` prose-sweep recipe are
  updated too.
- Not changed, deliberately: `ask_copilot.sh`. This account's Copilot accepts only `auto`
  (`ask_copilot.sh:18-22`), so the reserve seat cannot name the new tiers; enabling them is an
  account setting. The `gpt-5.6-terra` header specimen at protocol §6 is also left as is, because
  it is a filled 2026-09-21 example and rewriting the id would falsify it.

## Code jury

Both auditors found no defect that blocks the change. All findings were checked against the
source:

- CONFIRMED (Codex, Grok): the §9.1 defaults paragraph still said the defaults reproduce the
  vendor config. Rewritten.
- CONFIRMED (Codex, Grok): the new comment said P2c "holds" the default to the terra row. P2c
  (`check_practice_drift.py:545-551`) only requires the default to be allowlisted and to appear in
  the Codex column, so `gpt-6-luna` would pass. Narrowed to what the checker does.
- CONFIRMED (Grok): "the same rule that keeps older groks" was false. Groks are excused in
  records only (`:524-527`), and a row naming one is P2c DRIFT (`:537-542`). An id on the
  by-design line is excused everywhere. The sentence now says so and warns against it. The
  missing checker rule is **DR-205** `[W]` (Christian: fix the text, register the rule).
- CONFIRMED (Grok): "the gpt-5.6 ids and everything older" overclaimed, because
  `gpt-5.5`/`5.4`/`5.4-mini` are not on the line. Dropped.
- PARTLY REFUTED (Grok, "vendor ids unchecked"): terra was accepted live. Luna and sol have not
  been called.
- Pre-existing, not changed: `cross_review.sh` duplicates the Codex arm, which P2b/P2e do not
  read.

## Open

- DR-205: retired-id anchor for `check_practice_drift.py`.
- First live `gpt-6-luna` call: the next prose sweep or `--quick` run. A vendor rejection there
  fails closed.
