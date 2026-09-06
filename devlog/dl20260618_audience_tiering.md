# Devlog 2026-06-18 — Audience-Tiered Collaboration Artifacts

**Date:** 2026-06-18
**Topic:** Refactor the AI-collaboration rules so *audience* (reader) is the organizing axis; move the AI-to-AI exchange to ephemeral, sharded, per-task files.
**AIs involved:** Claude (implementation). No Codex/Grok audit this session (the change is to docs + the wrapper scripts themselves; verification was by direct test of the resolver).
**Claude Confidence:** high on the doc/script edits and the verified wrapper behavior; medium on the one open question below (the delta hook).
**Literature References:** N/A (collaboration-process change, not a formulation).

## Summary

Introduced a first-class distinction by **reader**:

- **AI-only artifacts** — read only by AIs (the exchange). Ephemeral, machine-parseable, disposable, never committed.
- **AI+human artifacts** — read by both (`./todo/` progress files, `./devlog/` session logs). Durable, tracked, curated.

Lifetime follows audience: an artifact is ephemeral *iff no human ever needs to return to it*. The tier boundary is a **distillation** step (signal lifted into the devlog before the scratch is swept), which is what makes aggressive GC of the AI-only tier safe.

The concrete structural change: the AI-to-AI exchange moved from a single durable `./todo/ai_exchange.md` (now ~4900 lines, far past its archive threshold) to **sharded, ephemeral, per-task files** at `./tmp/ai_exchange/<slug>.md`, one topic per file.

## Changes Made

**`doc/ai_collaboration_protocol.md`** (authoritative — 7 edits):
1. §2 retitled to "Audience Tiers & the Communication Channel"; added the Audience-Tiers definition + the "format follows audience" rule. (Folded into §2 rather than inserting a new numbered section, to avoid renumbering the many `§N` cross-references.)
2. §2 channel changed to `./tmp/ai_exchange/<slug>.md` (sharded, ephemeral, AI-only); entry format unchanged; Rules updated (one file per topic, no archiving).
3. §6: the devlog is now explicitly the durable AI+human distillation of the ephemeral exchange; conclusions are lifted into the devlog (and `./todo/` progress files) before the scratch is GC-eligible.
4. §7: `./tmp/ai_exchange/` added to the always-allowed write list.
5. §10 File Summary: added an **Audience** column; added the `./tmp/ai_exchange/<slug>.md` row (AI-only, "Ephemeral; safe to delete anytime — recommended sweep of files older than 14 days"); marked the old `./todo/ai_exchange.md` as frozen legacy.
6. Nonfree exception reframed in audience terms: the AI-only shared exchange is a vendor-bound leak surface, so nonfree work suppresses the AI-only tier entirely and keeps only the local human-readable record.
7. §9: documented that the wrappers resolve a per-task `./tmp/ai_exchange/<slug>.md` path; §8 Quick Start and the §9 workflow updated to the new path.

**`AGENTS.md`:** one pointer line under "Collaboration Protocol" naming the audience-tier principle and pointing to §2/§6.

**`CLAUDE.md`:** updated both `todo/ai_exchange.md` references (the "Key points" bullet and the Minimum Session Compliance Checklist) to the new scheme. (Historical references in `./devlog/` and `./todo/` were left untouched — forward-only migration.)

**`.claude/scripts/ask_codex.sh` and `ask_grok.sh`:** see "Wrapper behavioral changes" below.

**Infra:** created `./tmp/ai_exchange/`; added an explicit, commented `tmp/ai_exchange/` entry to `.gitignore` (it was already covered by the existing `tmp/` rule — kept explicit so it survives any future narrowing).

## Wrapper behavioral changes (review before relying on them)

Both `ask_codex.sh` and `ask_grok.sh` now resolve the exchange path instead of hardcoding `todo/ai_exchange.md`:

- **Slug precedence:** `$AI_EXCHANGE_SLUG` (sanitized to `lowercase_underscore`) → else `sess_<first-8-of-$CLAUDE_CODE_SESSION_ID>` → else `scratch`. Resolved path: `./tmp/ai_exchange/<slug>.md`.
- **Verified behavior** (live test of the resolver, no AI invoked):
  - default (session id present): `tmp/ai_exchange/sess_798c6b0c.md`
  - `AI_EXCHANGE_SLUG=periodic_thin_cut`: `tmp/ai_exchange/periodic_thin_cut.md`
  - no session id: `tmp/ai_exchange/scratch.md`
- The wrappers `mkdir -p` the dir and create the file if absent. `ask_grok.sh`'s previous *refuse-to-create* guard was relaxed for this ephemeral tier (it still keeps every safety behavior: read-only sandbox, the "sandbox could not be applied" hard-fail, timestamped `# CODEX`/`# GROK` headers, stdout echo).
- The audit prompt now tells the auditor to **read the existing thread at the resolved path first** before auditing.
- Both pass `bash -n`. To pin a topic, export `AI_EXCHANGE_SLUG=<topic>` before invoking.

## Session-ID availability (the flagged uncertainty — RESOLVED)

The task asked me to verify, not assume, that a session tag reaches the wrapper subprocess. **Verified:** `CLAUDE_CODE_SESSION_ID` (a 36-char UUID, e.g. `798c6b0c-...`) is exported and *propagates to a spawned script* (confirmed by running a standalone probe script, exactly as the wrappers are spawned). So the session-tag fallback is real, not assumed; I used `sess_<first 8 hex>` for a clean, collision-resistant filename. If this var is ever absent, the resolver degrades to `scratch`.

## Delta hook — UPDATED (follow-up approved in a later message)

`.claude/hooks/ai_exchange_delta.sh` was rewritten to track the **newest** file under `tmp/ai_exchange/` (i.e. whichever topic file a wrapper last appended to) instead of the hardcoded legacy `todo/ai_exchange.md`. State is now path-aware (`<path>\t<size>`), so switching topic/session files resets the offset cleanly; the older single-number state format degrades to a safe reset. Tested end-to-end: first-read injects + records state, no-new-content is silent, an appended entry shows only the delta, and a non-trigger prompt exits silently. `bash -n` clean.

## Open Questions / Follow-ups
- `ask_grok.sh` is currently git-ignored/untracked (only `ask_codex.sh` is un-ignored in `.gitignore`). Pre-existing; not changed here. Flagging in case both wrappers should be tracked.
- Migration is forward-only: the existing `todo/ai_exchange.md` (~4900 lines) and its archives were left intact as historical record.

## Files Updated

- doc/ai_collaboration_protocol.md
- AGENTS.md
- CLAUDE.md
- .claude/scripts/ask_codex.sh
- .claude/scripts/ask_grok.sh
- .gitignore
- tmp/ai_exchange/ (created; git-ignored)
- devlog/dl20260618_audience_tiering.md (this file) + devlog/README.md
