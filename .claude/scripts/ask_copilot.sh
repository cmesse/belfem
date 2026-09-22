#!/bin/bash
# ask_copilot.sh — request an independent audit from GitHub Copilot CLI and record it in the
# AI-only exchange. Third seat: a substitute for a Codex or Grok seat that is at its usage
# limit, never a replacement for both.
#
# Usage:
#   ask_copilot.sh "Audit this specific claim or file section..."
#   ask_copilot.sh -      (reads prompt from stdin)
#   echo "prompt" | ask_copilot.sh
#   AI_EXCHANGE_SLUG=periodic_thin_cut ask_copilot.sh "..."   # pin the per-task topic file
#   COPILOT_BIN=/path/to/copilot ask_copilot.sh "..."         # override the CLI location
#
# Depth selection (set both; the stamp says which knob was left unchosen):
#   COPILOT_MODEL   auto | claude-fable-5.1 | claude-fable-5 | claude-opus-5 | claude-opus-4.8 |
#                   claude-opus-4.7 | claude-sonnet-5 | claude-haiku-4.5 | gpt-6-astra |
#                   gpt-5.6-terra | gpt-5.6-sol | gpt-5.6-luna | gpt-5.5 | gpt-5.4 | gpt-5.4-mini |
#                   gemini-3.8-flash | kimi-k3
#                   (default auto). `auto` is the only name this account accepts on the
#                   command line today: it routes to gpt-5.6-luna and ignores the effort knob.
#                   Every other name is in the Copilot catalog but policy-disabled for this
#                   account, so it is refused before any request is billed; enabling a model
#                   is a toggle on github.com/settings/copilot/features, not a wrapper change.
#   COPILOT_EFFORT  low|medium|high|xhigh (default medium); ignored with a warning by `auto`
#
# Example:
#   COPILOT_MODEL=claude-opus-5 COPILOT_EFFORT=high ask_copilot.sh "..."
#
# The response is appended under a # COPILOT header to the per-task exchange file and echoed
# to stdout, matching ask_codex.sh / ask_grok.sh, so cross_review.sh can collect it the same way.
#
# Read-only is enforced by the CLI's own permission rules, not by trust: every file-writing
# tool is denied, the shell is allowed only for a fixed list of read-only commands, the
# GitHub MCP server is off, and the working directory is never trusted (no hooks, skills or
# plugins of the checkout are loaded). Denial rules beat every allow rule in the CLI.

set -u

# --- Resolve the Copilot CLI ---------------------------------------------------
if [ -n "${COPILOT_BIN:-}" ]; then
    COPILOT="$COPILOT_BIN"
elif command -v copilot >/dev/null 2>&1; then
    COPILOT="$(command -v copilot)"
else
    COPILOT="$HOME/.local/bin/copilot"
fi
if [ ! -x "$COPILOT" ]; then
    echo "ask_copilot.sh: Copilot CLI not found or not executable: $COPILOT" >&2
    echo "Install the GitHub Copilot CLI, put it on \$PATH, or set COPILOT_BIN=/path/to/copilot." >&2
    exit 1
fi

BELFEM="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

# --- Depth selection: model + reasoning effort -------------------------------
COPILOT_MODEL_CHOSEN=1
COPILOT_EFFORT_CHOSEN=1
if [ -z "${COPILOT_MODEL:-}" ];  then COPILOT_MODEL="auto";    COPILOT_MODEL_CHOSEN=0;  fi
if [ -z "${COPILOT_EFFORT:-}" ]; then COPILOT_EFFORT="medium"; COPILOT_EFFORT_CHOSEN=0; fi

case "$COPILOT_MODEL" in
    auto|claude-fable-5.1|claude-fable-5|claude-opus-5|claude-opus-4.8|claude-opus-4.7|claude-sonnet-5|claude-haiku-4.5|gpt-6-astra|gpt-5.6-terra|gpt-5.6-sol|gpt-5.6-luna|gpt-5.5|gpt-5.4|gpt-5.4-mini|gemini-3.8-flash|kimi-k3) ;;
    *)
        echo "ask_copilot.sh: unknown COPILOT_MODEL '$COPILOT_MODEL'" >&2
        echo "Allowed: auto claude-fable-5.1 claude-fable-5 claude-opus-5 claude-opus-4.8 claude-opus-4.7 claude-sonnet-5 claude-haiku-4.5 gpt-6-astra gpt-5.6-terra gpt-5.6-sol gpt-5.6-luna gpt-5.5 gpt-5.4 gpt-5.4-mini gemini-3.8-flash kimi-k3" >&2
        exit 1
        ;;
esac

# Capped at xhigh for the same reason as the Codex wrapper: `max` buys nothing demonstrated.
case "$COPILOT_EFFORT" in
    low|medium|high|xhigh) ;;
    *)
        echo "ask_copilot.sh: unknown COPILOT_EFFORT '$COPILOT_EFFORT'" >&2
        echo "Allowed: low medium high xhigh" >&2
        exit 1
        ;;
esac

COPILOT_MODEL_TAG="$COPILOT_MODEL"
COPILOT_EFFORT_TAG="$COPILOT_EFFORT"
if [ "$COPILOT_MODEL_CHOSEN" -eq 0 ];  then COPILOT_MODEL_TAG="$COPILOT_MODEL [defaulted]";   fi
if [ "$COPILOT_EFFORT_CHOSEN" -eq 0 ]; then COPILOT_EFFORT_TAG="$COPILOT_EFFORT [defaulted]"; fi
if [ "$COPILOT_MODEL_CHOSEN" -eq 0 ] || [ "$COPILOT_EFFORT_CHOSEN" -eq 0 ]; then
    echo "ask_copilot.sh: depth not selected by the caller — using model=$COPILOT_MODEL effort=$COPILOT_EFFORT" >&2
fi
if [ "$COPILOT_MODEL" = "auto" ]; then
    # the router refuses an effort setting; the stamp records what auto resolved to
    EFFORT_ARGS=()
    COPILOT_EFFORT_TAG="n/a (auto routing ignores effort)"
else
    EFFORT_ARGS=( --reasoning-effort "$COPILOT_EFFORT" )
fi

# --- Resolve the per-task AI-only exchange file: ./tmp/ai_exchange/<slug>.md --
sanitize_slug() {
    printf '%s' "$1" | tr 'A-Z' 'a-z' | tr -c 'a-z0-9_' '_' | sed 's/_\{2,\}/_/g; s/^_//; s/_$//'
}
if [ -n "${AI_EXCHANGE_SLUG:-}" ]; then
    SLUG="$(sanitize_slug "$AI_EXCHANGE_SLUG")"
elif [ -n "${CLAUDE_CODE_SESSION_ID:-}" ]; then
    SLUG="sess_$(sanitize_slug "${CLAUDE_CODE_SESSION_ID:0:8}")"
else
    SLUG="scratch"
fi
[ -n "$SLUG" ] || SLUG="scratch"
EXCHANGE_DIR="$BELFEM/tmp/ai_exchange"
EXCHANGE="$EXCHANGE_DIR/$SLUG.md"
mkdir -p "$EXCHANGE_DIR"
[ -f "$EXCHANGE" ] || : > "$EXCHANGE"

# --- Collect the prompt from argument or stdin ---------------------------------
if [ $# -eq 0 ] || [ "$1" = "-" ]; then
    AUDIT_PROMPT=$(cat)
else
    AUDIT_PROMPT="$*"
fi
if [ -z "$AUDIT_PROMPT" ]; then
    echo "Usage: ask_copilot.sh \"<audit prompt>\"  OR  ask_copilot.sh - (reads stdin)" >&2
    exit 1
fi

FULL_PROMPT="You are the Copilot AI conducting a read-only code audit for the BELFEM project.
Before answering, read these two files in the project root:
  - AGENTS.md
  - doc/ai_collaboration_protocol.md
They define the auditor role (precision audit), the audit checklist, and the required format.

This task's AI-only exchange thread is at:
  $EXCHANGE
Read it first (it may be empty or hold prior entries) so you have the conversation context before auditing.

Key rules:
- Prefer independent verification over agreement; challenge assumptions.
- State a confidence tier (high / medium / low) and ground it in evidence.
- Format your response starting with '## Audit & Verdict' then 'Confidence: <tier>'.
- Reference specific file paths and line numbers wherever possible.
- This is a read-only audit: do not modify any files. You have no file-writing tool and only
  read-only shell commands; do not try to work around that.

AUDIT REQUEST:
$AUDIT_PROMPT"

TMPOUT=$(mktemp "${TMPDIR:-/tmp}/copilot_response_XXXXXX")
TMPERR=$(mktemp "${TMPDIR:-/tmp}/copilot_err_XXXXXX")
TMPPROMPT=$(mktemp "${TMPDIR:-/tmp}/copilot_prompt_XXXXXX")
LOGDIR=$(mktemp -d "${TMPDIR:-/tmp}/copilot_logs_XXXXXX")
trap 'rm -rf "$TMPOUT" "$TMPERR" "$TMPPROMPT" "$LOGDIR"' EXIT

# The prompt goes through stdin, never argv: a large review bundle exceeds the per-argument
# limit, stdin does not. `-p -` is not documented, so the prompt is piped and -p is omitted:
# the CLI reads a piped stdin as the non-interactive prompt.
printf '%s\n' "$FULL_PROMPT" > "$TMPPROMPT"

# Read-only shell allow list. Anything not listed is denied in non-interactive mode, and
# `write` is denied outright, so an edit attempt fails instead of landing in the checkout.
READ_TOOLS=(
    --deny-tool 'write'
    --allow-tool 'shell(cat)' --allow-tool 'shell(head)' --allow-tool 'shell(tail)'
    --allow-tool 'shell(sed)' --allow-tool 'shell(awk)' --allow-tool 'shell(grep)'
    --allow-tool 'shell(rg)'  --allow-tool 'shell(ls)'  --allow-tool 'shell(find)'
    --allow-tool 'shell(wc)'  --allow-tool 'shell(nl)'  --allow-tool 'shell(diff)'
    --allow-tool 'shell(git diff)' --allow-tool 'shell(git log)' --allow-tool 'shell(git show)'
    --allow-tool 'shell(git status)' --allow-tool 'shell(git blame)' --allow-tool 'shell(git grep)'
)

"$COPILOT" \
    -C "$BELFEM" \
    --model "$COPILOT_MODEL" \
    "${EFFORT_ARGS[@]}" \
    -s \
    --no-auto-update \
    --disable-builtin-mcps \
    --log-dir "$LOGDIR" \
    "${READ_TOOLS[@]}" \
    < "$TMPPROMPT" \
    > "$TMPOUT" 2>"$TMPERR"

COPILOT_EXIT=$?
if [ $COPILOT_EXIT -ne 0 ]; then
    echo "copilot failed (exit $COPILOT_EXIT):" >&2
    cat "$TMPERR" >&2
    exit $COPILOT_EXIT
fi

RESULT=$(cat "$TMPOUT")
if [ -z "$RESULT" ]; then
    echo "copilot returned an empty response" >&2
    cat "$TMPERR" >&2
    exit 1
fi

# `auto` picks the model per request; the stamp carries what it resolved to, read from the
# session log, so a finding stays attributable to a model and not to a router.
if [ "$COPILOT_MODEL" = "auto" ]; then
    RESOLVED=$(grep -ho '"model":"[^"]*"' "$LOGDIR"/*.log 2>/dev/null | grep -v '"auto"' | sort | uniq -c | sort -rn | head -1 | sed 's/.*"model":"\([^"]*\)".*/\1/')
    [ -n "$RESOLVED" ] && COPILOT_MODEL_TAG="auto -> $RESOLVED"
fi

TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S %Z')
{
    printf '\n---\n\n'
    printf '# COPILOT %s  (model=%s, effort=%s)\n' "$TIMESTAMP" "$COPILOT_MODEL_TAG" "$COPILOT_EFFORT_TAG"
    printf '%s\n' "$RESULT"
} >> "$EXCHANGE"

printf '%s\n' "$RESULT"
