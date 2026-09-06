# ask_grok.sh aligned to Grok CLI 0.2.111 headless API

**Date:** 2026-07-23
**Purpose:** Diagnose and fix `.claude/scripts/ask_grok.sh` so Claude can spawn Grok as the secondary auditor the same way `ask_codex.sh` spawns Codex.
**AIs involved:** Grok (diagnosis + fix). Claude as the intended caller of the wrapper.

## Problem

`ask_grok.sh` was written against an older/incorrect mental model of the Grok headless API. Against **grok 0.2.111** (`~/.grok/docs/user-guide/14-headless-mode.md`, `22-permissions-and-safety.md`) several defaults were wrong or weak:

| Old assumption | Actual API (0.2.111) |
|----------------|----------------------|
| `--permission-mode plan` = read-only exploration | `plan` is accepted for compatibility but **does not enable a CLI permission policy**. Plan mode is a separate interactive planning feature. Only `default` and `bypassPermissions` change policy via the CLI flag. |
| Plain stdout is a reliable final-answer channel | Structured capture is better with `--output-format json` and the `.text` field (same role as Codex `--output-last-message`). |
| Large preamble on argv is fine | Prefer `--prompt-file` for large multi-line role preambles (quoting / ARG limits / cleanliness). |
| Permission mode alone enforces read-only | Real read-only contract is `--sandbox read-only` **plus** a tool allowlist (`--tools read_file,grep,list_dir`). Without the allowlist, headless cancels shell tools that would prompt, but config `permission_mode = always-approve` can still let non-read tools fire (sandbox still blocks project writes). |

Secondary issues already partially handled in the prior draft (kept): intermittent empty responses (retry loop), fused `##` headings (perl split), missing `##` warning.

Older note (2026-06): `GROK_EFFORT=high` failed on `grok-build` with a 400 for `reasoningEffort`. **No longer true** for default model `grok-4.5` — `--effort high` works.

## Fix (`.claude/scripts/ask_grok.sh`)

- Invoke via `--prompt-file` + `--output-format json` + `--verbatim`; extract `.text` with a small `python3` helper; surface non-`EndTurn` `stopReason` as a warning (e.g. turn-budget cutoffs).
- Enforce read-only with `--sandbox read-only` and default `--tools read_file,grep,list_dir` (`GROK_TOOLS` env override). Drop default `--permission-mode plan`; optional `GROK_PERMISSION_MODE` only if the caller opts in.
- Keep retry loop, sandbox-not-applied hard fail, exchange `# GROK` attribution, stdout echo (Codex parity).
- `.gitignore` already un-ignores `ask_grok.sh` (local change).

## Verification

End-to-end smoke (`AI_EXCHANGE_SLUG=ask_grok_smoke_test_v2`):

```text
RC=0
## Audit & Verdict  (first line of AGENTS.md = "# AGENTS.md")
sandbox=read-only, tools used: [read_file]
```

Also verified separately: sandbox `read-only` blocks project file writes even if the agent claims to overwrite; `--effort high` returns normal JSON on grok-4.5.

## Caller tips (Claude)

```bash
AI_EXCHANGE_SLUG=<topic> .claude/scripts/ask_grok.sh "focused audit prompt..."
# larger multi-file audits:
GROK_MAX_TURNS=45 GROK_EFFORT=high AI_EXCHANGE_SLUG=<topic> .claude/scripts/ask_grok.sh "..."
# if literature web tools are needed (still no shell):
GROK_TOOLS=read_file,grep,list_dir,web_search,web_fetch .claude/scripts/ask_grok.sh "..."
```

**Do not call bare `grok -p` for audits.** That path recreates the permission-cancel footgun.

---

## Follow-up (same day): narration-only / Cancelled root cause

Claude reported exit 0 with only lead-in narration ("I'll audit…"), sometimes empty; cardano/ferrari 3/3 fail streaks; empty-stdout gate missed the narration-only variant.

### Root cause (session-proven)

When headless grok batches `run_terminal_command` and that tool needs a permission prompt, headless **cancels the tool and ends the entire turn**:

```text
permission_resolved { tool: run_terminal_command, decision: cancelled }
turn_ended { outcome: cancelled, cancellation_category: permission_cancelled }
JSON: stopReason=Cancelled, text="I'll run that Python one-liner…"   # lead-in only
process exit: 0
```

Reproduced directly with `--permission-mode default`. Math audits that "want a quick python check" hit this almost deterministically → **streaky, non-independent failures**. Fresh retries without removing shell keep failing the same way. Successful prompts (trivial file-read, brainstorms that never shell out) pass first try.

Fusion of lead-in onto the first `##` with no newline is a separate text-plumbing artifact (still fixed by the perl split).

### Wrapper hardening (v2)

1. **Prevent:** default `--tools read_file,grep,list_dir` (refuse GROK_TOOLS that include shell); preamble forbids shell and requires final message to BE the audit.
2. **Detect:** quality gate rejects `stopReason=Cancelled`, missing `##`, or body `< GROK_MIN_CHARS` (default 200). Empty **and** narration-only both fail the gate (non-zero exit; nothing appended).
3. **Salvage:** on reject with a `sessionId`, one `--resume` finish pass ("write the full ## audit now, no shell") reuses files already read.
4. **Retry:** default 5 fresh attempts with quadratic backoff (cluster-aware), not 3 independent coin-flips.

### Verification

- Smoke: AGENTS.md first-line audit, attempt 1/2, RC=0, `# GROK` recorded.
- Ferrari Q1 mini (historically shell-hungry): attempt 1/3, RC=0, body ~5 KB, 12 `##` headings, tools allowlist only.
