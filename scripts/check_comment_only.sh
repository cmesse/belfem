#!/usr/bin/env bash
# Gate for comment-only commits (doc/commenting_guidelines.md §13, the sweep rules).
#
# For every C/C++ and Fortran source that differs between BASE and NEW, strip the
# comments from both revisions (scripts/strip_comments.py) and diff what is left.
# A non-empty diff means a token outside a comment changed, and the gate fails.
# Line numbers embedded by BELFEM_ERROR / BELFEM_ASSERT through __LINE__ are
# expected to move, which is why this compares stripped source and not binaries.
#
# Usage:
#   scripts/check_comment_only.sh                # HEAD vs the working tree
#   scripts/check_comment_only.sh <ref>          # <ref> vs the working tree
#   scripts/check_comment_only.sh <a>..<b>       # commit range
#
# Exit 0 only when every changed file was stripped on both sides and differs in
# comments only. A file the stripper cannot read is reported and fails the gate;
# the script never stops early, so the summary line is always printed.
set -uo pipefail

RANGE="${1:-HEAD}"
if [[ "$RANGE" == *..* ]]; then
    BASE="${RANGE%%..*}"; NEW="${RANGE##*..}"
else
    BASE="$RANGE"; NEW=""            # empty NEW = working tree
fi

REPO="$(git rev-parse --show-toplevel)" || exit 2
cd "$REPO" || exit 2
STRIP="$REPO/scripts/strip_comments.py"

if [ -n "$NEW" ]; then
    FILES=$(git diff --name-only "$BASE" "$NEW" -- '*.cpp' '*.hpp' '*.h' '*.c' '*.f90')
else
    FILES=$(git diff --name-only "$BASE" -- '*.cpp' '*.hpp' '*.h' '*.c' '*.f90')
fi

if [ -z "$FILES" ]; then
    echo "check_comment_only: no source files differ ($RANGE)"
    exit 0
fi

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

fetch() {     # $1 = ref ("" = working tree), $2 = path, $3 = output; empty file if absent
    if [ -z "$1" ]; then
        if [ -f "$2" ]; then cp "$2" "$3"; else : > "$3"; fi
    else
        git show "$1:$2" > "$3" 2>/dev/null || : > "$3"
    fi
}

FAIL=0
COUNT=0
for f in $FILES; do
    COUNT=$((COUNT + 1))
    case "$f" in
        *.f90) LANG_ARG=fortran ;;
        *)     LANG_ARG=c ;;
    esac
    fetch "$BASE" "$f" "$TMP/old.src"
    fetch "$NEW"  "$f" "$TMP/new.src"
    if ! python3 "$STRIP" "$LANG_ARG" < "$TMP/old.src" > "$TMP/old.txt"; then
        echo "CANNOT STRIP ($BASE): $f"; FAIL=1; continue
    fi
    if ! python3 "$STRIP" "$LANG_ARG" < "$TMP/new.src" > "$TMP/new.txt"; then
        echo "CANNOT STRIP (${NEW:-working tree}): $f"; FAIL=1; continue
    fi
    if diff -q "$TMP/old.txt" "$TMP/new.txt" > /dev/null; then
        echo "comment-only:     $f"
    else
        echo "NOT comment-only: $f"
        diff -u --label "$f (base)" --label "$f (new)" "$TMP/old.txt" "$TMP/new.txt" | head -40
        FAIL=1
    fi
done

if [ "$FAIL" -ne 0 ]; then
    echo "check_comment_only: FAILED ($RANGE, $COUNT files) -- a non-comment token changed or a file could not be stripped"
    exit 1
fi
echo "check_comment_only: OK ($RANGE, $COUNT files) -- every changed source differs in comments only"
