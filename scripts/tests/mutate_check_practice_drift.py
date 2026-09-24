#!/usr/bin/env python3
# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""
Mutation suite for `scripts/check_practice_drift.py`.

The checker's `--self-test` is its gate, and a green gate proves nothing unless it can go red.  Each
entry below edits one piece of the checker — disables a probe, or reintroduces a defect an earlier
version had — and the self-test must then fail.  A mutant that survives is a probe or a fix the
fixtures do not actually pin.

Every anchor is an exact substring of the checker's source and must occur there exactly once, so an
edit to the checker can break an anchor.  That is deliberate: change the checker and this suite
together, and add a mutant for each new probe or fixed defect.

Usage: python3 scripts/tests/mutate_check_practice_drift.py [checker] [workdir]
Defaults: the checker next to this directory, and a fresh temporary directory, removed afterwards.
The unmutated checker's self-test runs first and must pass: if it is already red, every mutant would
look killed, so the suite stops with exit 2, printing the child's exit status (or that it timed out)
and the last 20 lines of its output — which need not include the failing fixture, because the
self-test prints each FAIL where it occurs and only a count at the end.  An anchor that does not occur exactly once
is reported as BROKEN and fails the run.  Otherwise exit 0 only if every mutant is killed.  A mutant
that hangs is stopped after 180 s and counted as killed.

A kill is any non-zero exit of the mutant's self-test.  It shows the edit made the gate go red, not
which fixture did it.
"""

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

target = sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).resolve().parent.parent / "check_practice_drift.py")
created = len(sys.argv) <= 2
work = sys.argv[2] if not created else tempfile.mkdtemp(prefix="mutate_practice_drift_")
src = open(target, encoding="utf-8").read()


def self_test(path):
    """(exit status, combined output) of `path --self-test`; status None if it timed out."""
    try:
        r = subprocess.run([sys.executable, path, "--self-test"], capture_output=True,
                           encoding="utf-8", errors="replace", timeout=180)
        return r.returncode, r.stdout + r.stderr
    except subprocess.TimeoutExpired:
        return None, ""


def finish(code):
    if created:
        shutil.rmtree(work, ignore_errors=True)
    sys.exit(code)


# A kill means something only if the unmutated self-test passes.
status, output = self_test(target)
if status != 0:
    why = "timed out after 180 s" if status is None else f"exited {status}"
    print(f"the unmutated checker's --self-test {why}, so no kill below would mean anything; its last "
          "output:", file=sys.stderr)
    print("\n".join(output.strip().splitlines()[-20:]), file=sys.stderr)
    finish(2)

M = {
 'P2a disabled': ('findings.append((DRIFT, "P2a"', 'None and findings.append((DRIFT, "P2a"'),
 'P2b columns disabled': ('for model in sorted(models):', 'for model in []:'),
 'P2c grok default off': ('findings.append((DRIFT, "P2c", d,\n                         f"grok default', 'None and findings.append((DRIFT, "P2c", d,\n                         f"grok default'),
 'P2c codex rejected off': ('findings.append((DRIFT, "P2c", d, f"codex default {d} is rejected', 'None and findings.append((DRIFT, "P2c", d, f"codex default {d} is rejected'),
 'P2c codex no-row off': ('findings.append((DRIFT, "P2c", d, f"codex default {d} is named by no', 'None and findings.append((DRIFT, "P2c", d, f"codex default {d} is named by no'),
 'P2d disabled': ('findings.append((INFO, "P2d"', 'None and findings.append((INFO, "P2d"'),
 'P2e disabled': ('findings.append((DRIFT, "P2e"', 'None and findings.append((DRIFT, "P2e"'),
 'P3 disabled': ('for name in sorted(files - targets)]', 'for name in sorted(files - targets) if False]'),
 'P1a disabled': ('for f in sorted(declared) if not used.get(f)]', 'for f in sorted(declared) if False]'),
 'P1b disabled': ('for f in sorted(used) if f not in declared and', 'for f in sorted(used) if False and'),
 'P4 disabled': ('    if not old:\n        return []', '    if True:\n        return []'),
 'label capture': ('INDEX_ENTRY = re.compile(r"^(?:- |#{2,4} )\\[[^\\]]*\\]\\(([^)\\s]+)\\)")', 'INDEX_ENTRY = re.compile(r"^(?:- |#{2,4} )\\[([^\\]]+)\\]\\(()")'),
 'occurrence count': ('headers.update({m for l in _header_lines(p) for m in MODEL.findall(l)})', 'headers.update([m for l in _header_lines(p) for m in MODEL.findall(l)])'),
 'guard no-op': ('    if empty:\n        raise', '    if False:\n        raise'),
 'policy grok undeclared': ('declared = rows | exceptions | ({policy_grok} if policy_grok else set())', 'declared = rows | exceptions'),
 'arm -> message': ('allowlists = {s: wrapper_models(root / p)', 'allowlists = {s: wrapper_message_models(root / p)'),
 'no P2c precedence': ('if model in declared_for_use or model in reported:', 'if model in declared_for_use:'),
 'D-c11 exception unchecked': ('findings.append((DRIFT, "P2b", model,\n                             f"the by-design', 'None and findings.append((DRIFT, "P2b", model,\n                             f"the by-design'),
 'D-c12 grok row unchecked': ('findings.append((DRIFT, "P2c", model,\n                             f"the Grok column names {model}; the policy', 'None and findings.append((DRIFT, "P2c", model,\n                             f"the Grok column names {model}; the policy'),
 'D-c13 whole-table codex': ('elif d not in codex_rows:', 'elif d not in rows:'),
 'D-c13 arm not checked': ('if d not in allowlists.get("codex", set()):', 'if False:'),
 'D-c14 basename only': ('if os.path.dirname(joined) == base:', 'if True:'),
 'D-c15 P1c disabled': ('return [(DRIFT, "P1c", name,', 'return [] and [(DRIFT, "P1c", name,'),
 'D-c17 reversed tiebreak': ('key=lambda p: p.name)', 'key=lambda p: p.name, reverse=True)'),
 'D-c19 headed denominator': ('probe_template(declared, used, len(window), threshold, denom)', 'probe_template(declared, used, denom, threshold, denom)'),
 'D-c20 OSError uncaught': ('except (ExtractionError, OSError) as exc:', 'except ExtractionError as exc:'),
 'D-c21 policy unguarded': ('"policy grok": _newest_grok(allowlists["grok"]),', '"policy grok": "x",'),
 'D-c22 history is drift': ('declared_for_use = declared | allowlists.get("grok", set())', 'declared_for_use = declared'),
 'D-c16 recent unvalidated': ('"--recent", type=_at_least(1)', '"--recent", type=int'),
 'D-c16 threshold unvalidated': ('"--threshold", type=_percent', '"--threshold", type=float'),
 'D-c23 codex column by prefix': ('codex_rows = columns.get("codex", set())', 'codex_rows = {r for r in rows if _seat(r) == "codex"}'),
 'D-c23 P2b column by prefix': ('for seat, models in sorted(columns.items()):', 'for seat, models in sorted({s: {m for m in rows if _seat(m) == s} for s in columns}.items()):'),
 'D-c24 resolve() again': ('joined = os.path.normpath(os.path.join(base, target))', 'joined = str((Path(base) / target).resolve())'),
 'D-c25 P1c needs P1a': ('missing_date = devlogs_without_date(window)', 'missing_date = devlogs_without_date(window) if not used.get("Date") else []'),
 'D-c26 errors always silent': ('if not quiet:\n            print(f"extraction failed', 'if False:\n            print(f"extraction failed'),
 'D-c27 headers guarded again': ('"devlog window": window,', '"devlog window": used,'),
 'D-c28 tildes not fences': ('FENCE = re.compile(r"^ {0,3}(`{3,}|~{3,})(.*)$")', 'FENCE = re.compile(r"^ {0,3}(`{3,})(.*)$")'),
 'D-c28 info-string rule off': ('if not m or (m.group(1)[0] == "`" and "`" in m.group(2)):', 'if not m:'),
 'D-c30 width unchecked': ('if len(_cells(row)) != len(header):', 'if False:'),
 'D-c30 duplicate header allowed': ('if header.count(name) != 1:', 'if name not in header:'),
 'D-c30 strip all edge pipes': ('s = row.strip()\n    if s.startswith("|"):\n        s = s[1:]\n    if s.endswith("|"):\n        s = s[:-1]', 's = row.strip().strip("|")'),
 'D-c30 second block included': ('elif rows:\n            break', 'elif rows:\n            pass'),
 'O6 P1a on headless window': ('if headed_n:', 'if True:'),
 'header ignores fences': ('        if fenced:\n            continue\n        if line.startswith("## "):', '        if False:\n            continue\n        if line.startswith("## "):'),
 'table ignores fences': ('if (opened := _fence_open(line)) is not None:\n            if rows:\n                break\n            fence = opened\n            continue', 'if False:\n            pass'),
 'D-c28 closer length ignored': ('run[1] >= fence[1]', 'True'),
 'D-c42 tilde closes backtick': ('run[0] == fence[0] and ', ''),
 'D-c35 closer takes any space': ('if m and not m.group(2).strip(" \\t"):', 'if m and not m.group(2).strip():'),
 'D-c35 closer takes any text': ('if m and not m.group(2).strip(" \\t"):', 'if m:'),
 'D-c40 unclosed opener swallows': ('if opened is not None and longest[opened[0]][i + 1] >= opened[1]:', 'if opened is not None:'),
 'D-c40 closer counted at the opener': ('longest[opened[0]][i + 1] >= opened[1]', 'longest[opened[0]][i] >= opened[1]'),
 'D-c40 strict length': ('longest[opened[0]][i + 1] >= opened[1]', 'longest[opened[0]][i + 1] > opened[1]'),
 'D-c40 suffix not carried': ('longest[c][i] = longest[c][i + 1]', 'longest[c][i] = 0'),
 'D-c39 end ignores fences': ('if not fenced[i] and (lines[i].startswith("### ")', 'if (lines[i].startswith("### ")'),
 '_lines back to splitlines': ('    parts = re.split(r"\\r\\n|\\r|\\n", Path(path).read_text(encoding="utf-8", errors="replace"))', '    parts = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()'),
 'v6 table reader (the regression)': ('def _table_rows(proto_path):\n    """The first contiguous block of `|` lines in the depth-selection span, outside fences.\n\n    A fence opened here that never closes is deliberately NOT demoted to text, unlike in\n    `_header_lines`: demoting it would expose any example table inside it as the depth table.\n    Left open, it hides the rest of the span, no rows are found, and the required-extract guard\n    fails the run closed.\n    """\n    rows, fence = [], None\n    for line in _depth_span(proto_path):\n        if fence is not None:\n            if _fence_closes(line, fence):\n                fence = None\n            continue\n        if (opened := _fence_open(line)) is not None:\n            if rows:\n                break\n            fence = opened\n            continue\n        if line.startswith("|"):\n            rows.append(line)\n        elif rows:\n            break\n    return rows\n\n\n', 'def _table_rows(proto_path):\n    """The first contiguous block of `|` lines in the depth-selection span, outside fences.\n\n    Fences matter only before the table: once rows have started, any non-`|` line — a fence\n    opener included — ends the block.\n    """\n    span = _depth_span(proto_path)\n    rows = []\n    for line, fenced in zip(span, _fenced_flags(span)):\n        if fenced:\n            if rows:\n                break\n            continue\n        if line.startswith("|"):\n            rows.append(line)\n        elif rows:\n            break\n    return rows\n\n\n'),
 'D-c43 ambiguity accepted': ('    if len(starts) > 1:', '    if False:'),
 'D-c43 prefix anchor': ('not fenced[i] and DEPTH_HEADING.match(l)]', 'not fenced[i] and l.startswith(DEPTH_ANCHOR)]'),
 'D-c39 start ignores fences': ('if not fenced[i] and DEPTH_HEADING.match(l)', 'if DEPTH_HEADING.match(l)'),
 'D-c45 any whitespace': ('DEPTH_HEADING = re.compile(r"^### 9\\.1(?:[ \\t]|$)")', 'DEPTH_HEADING = re.compile(r"^### 9\\.1(?:\\s|$)")'),
 'D-c46 enclosing fence ignored': ('if (opener := _open_fence_at(lines, start)) is not None:', 'if (opener := None) is not None:'),
 'D-c46 closers never close': ('            if _fence_closes(lines[i], fence):\n                fence = None\n            continue\n        if (opened := _fence_open(lines[i])) is not None:', '            continue\n        if (opened := _fence_open(lines[i])) is not None:'),
 'D-c47 no line numbers': ('f"fenced blocks, at lines {\', \'.join(str(i + 1) for i in starts)}")', 'f"fenced blocks")'),
 'D-c48 space only': ('DEPTH_HEADING = re.compile(r"^### 9\\.1(?:[ \\t]|$)")', 'DEPTH_HEADING = re.compile(r"^### 9\\.1(?: |$)")'),
 'D-c48 no end-of-line': ('DEPTH_HEADING = re.compile(r"^### 9\\.1(?:[ \\t]|$)")', 'DEPTH_HEADING = re.compile(r"^### 9\\.1(?:[ \\t])")'),
 'D-c49 duplicate lines off by one': ('str(i + 1) for i in starts)', 'str(i) for i in starts)'),
 'D-c49 heading line off by one': ('f"depth anchor at line {start + 1} lies', 'f"depth anchor at line {start} lies'),
 'D-c49 opener line off by one': ('f"{opener + 1} and never closed")', 'f"{opener} and never closed")'),
 'D-m5 commit bodies counted per occurrence': ('commits.update(set(MODEL.findall(body)))', 'commits.update(MODEL.findall(body))'),
 'D-m5 trailing period absorbed': ('MODEL = re.compile(r"\\b((?:gpt|grok)-\\d+(?:\\.\\d+)*(?:-[A-Za-z][A-Za-z0-9]*)*)")', 'MODEL = re.compile(r"\\b((?:gpt|grok)-[\\d.]+(?:-[A-Za-z][A-Za-z0-9]*)*)")'),
 'D-m5 seat names admitted': ('MODEL = re.compile(r"\\b((?:gpt|grok)-\\d+(?:\\.\\d+)*(?:-[A-Za-z][A-Za-z0-9]*)*)")', 'MODEL = re.compile(r"\\b((?:gpt|grok|Opus)[- ]\\d+(?:\\.\\d+)*(?:-[A-Za-z][A-Za-z0-9]*)*)")'),
 'D-m5 window date order dropped': ('sorted(by_name, key=date, reverse=True)', 'sorted(by_name, key=date)'),
 'D-m5 --strict ignored': ('return 1 if (strict and drift) else 0', 'return 0'),
}
killed = broken = 0
for name, (a, b) in M.items():
    found = src.count(a)
    if found != 1:
        print(f"  BROKEN   {name}: anchor occurs {found} times, must be exactly once")
        broken += 1
        continue
    open(f"{work}/mut.py", "w", encoding="utf-8").write(src.replace(a, b))
    status, _ = self_test(f"{work}/mut.py")
    ok = status != 0
    if status is None:
        name += " (timeout)"
    killed += ok
    print(f"  {'killed  ' if ok else 'SURVIVED'} {name}")
print(f"\n{killed}/{len(M)} mutants killed" + (f", {broken} BROKEN anchors" if broken else ""))
finish(0 if killed == len(M) else 1)
