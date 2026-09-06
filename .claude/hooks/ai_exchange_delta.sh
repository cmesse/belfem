#!/usr/bin/env bash
# Injects new content from the AI-only exchange when an AI audit is pasted/mentioned.
# The exchange is now sharded into ephemeral per-task files under tmp/ai_exchange/
# (see doc/ai_collaboration_protocol.md §2). This hook watches the most-recently-
# modified file there — i.e. whichever topic file a wrapper last appended to — and
# shows only the delta since the last read. State is path-aware so switching topic
# files (or sessions) resets the offset cleanly.

INPUT=$(cat)
PROMPT=$(echo "$INPUT" | jq -r '.prompt // ""')

echo "$PROMPT" | grep -qiE \
  '(codex|ai[[:space:]]*(audit|review|findings)|(from|by)[[:space:]]+(codex|gpt|gemini|another[[:space:]]+ai)|## (finding|issue|confidence)|audit[[:space:]]+(report|thread|result))' \
  || exit 0

# repository root = two levels above this script (.claude/hooks/)
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXCHANGE_DIR="$ROOT/tmp/ai_exchange"
STATE="$ROOT/.claude/ai_exchange_pos.txt"

# newest exchange file = the one a wrapper most recently appended to
EXCHANGE=$(ls -t "$EXCHANGE_DIR"/*.md 2>/dev/null | head -1)
[ -n "$EXCHANGE" ] || exit 0

# state format: "<path>\t<size>" (older single-number state resets safely)
STATE_PATH=$(cut -f1 "$STATE" 2>/dev/null)
STATE_SIZE=$(cut -f2 "$STATE" 2>/dev/null)
case "$STATE_SIZE" in (''|*[!0-9]*) STATE_SIZE=0 ;; esac

SIZE=$(stat -c %s "$EXCHANGE" 2>/dev/null || echo 0)

# reset offset if the newest file changed (new topic/session) or was truncated
if [ "$STATE_PATH" = "$EXCHANGE" ]; then LAST=$STATE_SIZE; else LAST=0; fi
[ "$SIZE" -lt "$LAST" ] && LAST=0

[ "$SIZE" -gt "$LAST" ] || exit 0

BASENAME=$(basename "$EXCHANGE")
if [ "$LAST" -eq 0 ] && [ "$SIZE" -gt 15000 ]; then
    NEW=$(tail -c 15000 "$EXCHANGE")
    MSG="Last ~15KB of tmp/ai_exchange/$BASENAME (first read; future triggers show only new entries):"
else
    NEW=$(tail -c +$((LAST + 1)) "$EXCHANGE")
    MSG="New content in tmp/ai_exchange/$BASENAME since last read:"
fi

printf '%s\t%s\n' "$EXCHANGE" "$SIZE" > "$STATE"
printf "%s" "$NEW" | jq -Rs --arg msg "$MSG" \
    '{hookSpecificOutput: {hookEventName: "UserPromptSubmit", additionalContext: ($msg + "\n\n" + .)}}'
