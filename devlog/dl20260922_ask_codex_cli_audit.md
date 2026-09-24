# ask_codex CLI audit

**Date:** 2026-09-22
**Purpose:** Review `.claude/scripts/ask_codex.sh` against the installed Codex CLI and the collaboration protocol, with particular attention to evolving CLI arguments.
**Module:** `meta` (AI collaboration tooling)
**AIs involved:** Codex `gpt-6-astra`/high
**Status:** REVIEWED — no source changes; one small follow-up is recommended before the next wrapper edit.
**Verification:** `codex-cli 0.156.0` local `codex exec --help`; `bash -n .claude/scripts/ask_codex.sh` clean. No audit request was dispatched.

## Outcome

The active invocation is compatible with the installed CLI: `exec`, `--sandbox read-only`, `-C`, `-m`, `-c model_reasoning_effort=<effort>`, stdin through positional `-`, and `--output-last-message` are all accepted by Codex CLI 0.156.0. The wrapper correctly puts options before the stdin marker and avoids the earlier argv-size/open-stdin failure class.

`codex exec` now also offers `--ephemeral` (run without persistent session files). Add it on the next approved wrapper change: the wrapper is a one-shot peer-audit path and its only intended record is the disposable `tmp/ai_exchange/<slug>.md` file. Without it, Codex can retain a separate local session transcript outside that lifecycle. Confidence: high (installed CLI help; protocol §2 and §10).

## Secondary observations

- The wrapper intentionally lets the user's Codex configuration load. `--ignore-user-config` exists in CLI 0.156.0, but should not be added without a policy decision: it could improve reproducibility while breaking installations that need configured providers or settings. Confidence: medium.
- The model validator still accepts historical `gpt-5.5`, `gpt-5.4`, and `gpt-5.4-mini`, whereas protocol §9.1 only deliberately allowlists `gpt-5.6-sol` outside its named tiers. Reconcile that wording/list when model policy is next updated; no current invocation failure was observed. Confidence: medium.
- Parallel Codex/Grok appenders do not lock the shared exchange file, so their multiple `printf` writes can theoretically interleave. This is pre-existing and distinct from the solved temporary-file collision. Confidence: medium.

## Files Updated

- `devlog/dl20260922_ask_codex_cli_audit.md`
- `devlog/README.md`
