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
Check the convention documents against practice.

`check_doc_claims.py` verifies those documents against the **tree**: build flags, make targets,
symbol existence.  This script covers the other direction — what the devlogs, the commit bodies
and the auditor wrappers actually do — and catches a class of defect that is invisible to both a
reader and a review round, because each side is individually self-consistent:

  * a depth ruling is spoken, lands in a session's notes, and the table auditors read is never
    edited, so the two disagree for days;
  * a wrapper default is correct when pinned and nothing re-reads it when the policy moves;
  * a devlog is written and its index line is forgotten, so the entry is invisible to anyone
    navigating by the index;
  * a header template describes fields practice has left behind, because each entry is written
    from the last entry rather than from the template.

Two rules keep the probes honest.

**Anchors, never line numbers.**  Every extract is located by a searchable token, following
`doc/input_schema.yaml` and `check_doc_claims.py`.  A line number read from a working copy is
already stale by the time it is committed.

**Exceptions are declared structurally, never inferred from prose.**  A model may be deliberately
allowlisted and named by no depth row.  That is declared on one anchored line the protocol owns.
Reading the section's prose instead would be unsound: prose carries escape-hatch declarations *and*
descriptions of live defects, and no token search separates them — documenting a defect would
silently disable the probe that detects it.

Tiers follow **normativity**, not the shape of the comparison.  DRIFT is a broken rule the document
states.  INFO is a convention moving where the document declines to impose one — the header field
list is explicitly "not a schema", so a session adding a field is not a defect.  Only DRIFT gates
under `--strict`.

`--self-test` is the gate, and every probe has at least one fixture that must make it fire and one
that must keep it silent.  The live run reports and asserts no exact set, so repairing a real drift
can never break the check.

The gate is itself tested by `scripts/tests/mutate_check_practice_drift.py`, whose anchors are exact
substrings of this file: change the two together, and add a mutant for every new probe.
"""

import argparse
import contextlib
import io
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

# A version is dot-separated integers; a suffix segment starts with a letter (gpt-6-astra,
# gpt-5.4-mini, grok-4.7-build-fast).  Seat names such as "Claude Opus 5" fall outside the
# alphabet, so they cannot leak into a model set, and a trailing sentence period is not absorbed.
MODEL = re.compile(r"\b((?:gpt|grok)-\d+(?:\.\d+)*(?:-[A-Za-z][A-Za-z0-9]*)*)")
GROK_VERSION = re.compile(r"grok-(\d+(?:\.\d+)*)")
FIELD = re.compile(r"^\*\*([A-Za-z][A-Za-z /]*):\*\*")
# The capture is the link DESTINATION inside (...), not the label inside [...]; the two differ
# whenever an entry is given a short name.
INDEX_ENTRY = re.compile(r"^(?:- |#{2,4} )\[[^\]]*\]\(([^)\s]+)\)")
FENCE = re.compile(r"^ {0,3}(`{3,}|~{3,})(.*)$")
CASE_ARM = re.compile(r"^\s*([A-Za-z0-9.|-]+)\)\s*;;")
OLD_TITLE = re.compile(r"^# Devlog \d{4}-\d{2}-\d{2}")
DEVLOG_DATE = re.compile(r"dl(\d{8})")

TEMPLATE_ANCHOR = "### Required Header"
TEMPLATE_SUCCESSOR = "### When to Write a Devlog"
DEPTH_ANCHOR = "### 9.1"
# The heading itself, not a prefix of it: `### 9.10` must not count as the depth section.
DEPTH_HEADING = re.compile(r"^### 9\.1(?:[ \t]|$)")
EXCEPTION_ANCHOR = "**Allowlisted but named by no row (by design):**"

WRAPPERS = {"codex": (".claude/scripts/ask_codex.sh", "CODEX_MODEL"),
            "grok": (".claude/scripts/ask_grok.sh", "GROK_MODEL")}

DRIFT, INFO = "DRIFT", "INFO"


class ExtractionError(RuntimeError):
    """A required extract came back empty or an anchor was not found."""


# --------------------------------------------------------------------------- extractors
# Each takes its path(s) explicitly so the self-test can point it at a fixture.

def _lines(path):
    """The file's lines, split on CommonMark line endings only.

    `str.splitlines()` also breaks on form feed, vertical tab and several Unicode separators,
    which Markdown treats as ordinary characters — a closing fence followed by a form feed would
    then look bare and close the fence.
    """
    parts = re.split(r"\r\n|\r|\n", Path(path).read_text(encoding="utf-8", errors="replace"))
    if parts and parts[-1] == "":
        parts.pop()
    return parts


def template_fields(proto_path):
    """Field names declared in the header template's fenced example."""
    lines = _lines(proto_path)
    try:
        start = next(i for i, l in enumerate(lines) if l.startswith(TEMPLATE_ANCHOR))
    except StopIteration:
        raise ExtractionError(f"anchor not found: {TEMPLATE_ANCHOR}")
    successor = next((i for i, l in enumerate(lines[start:], start)
                      if l.startswith(TEMPLATE_SUCCESSOR)), len(lines))

    # Fence-aware: the example's own first line begins with '#' (it is a title inside the
    # block), so a scan that treats '^#' as a heading would reject the very fence it must read.
    fence = [i for i, l in enumerate(lines[start:successor], start) if l.startswith("```")]
    if len(fence) < 2:
        raise ExtractionError("header template fence does not open and close before "
                              f"{TEMPLATE_SUCCESSOR!r}")
    return {m.group(1) for l in lines[fence[0] + 1:fence[1]] if (m := FIELD.match(l))}


def devlog_window(devlog_dir, recent):
    """The newest `recent` devlogs, ordered by filename date then filename.

    Several entries routinely share a date, so the date alone does not order the corpus; the
    filename, ascending, breaks the tie deterministically.  Sorting by name first and then stably
    by date keeps the two directions independent — a single reversed key would reverse both.
    """
    def date(p):
        m = DEVLOG_DATE.search(p.name)
        return m.group(1) if m else "0"
    by_name = sorted(Path(devlog_dir).glob("dl*.md"), key=lambda p: p.name)
    return sorted(by_name, key=date, reverse=True)[:recent]


def _fence_open(line):
    """(character, length) if `line` opens a fenced block, else None.

    CommonMark: three or more backticks or tildes, indented at most three spaces; a backtick
    fence's info string may not itself contain a backtick.
    """
    m = FENCE.match(line)
    if not m or (m.group(1)[0] == "`" and "`" in m.group(2)):
        return None
    return m.group(1)[0], len(m.group(1))


def _bare_run(line):
    """(character, length) if `line` is a bare fence run — a possible closer — else None.

    Bare means nothing after the run but spaces or tabs; a non-breaking space is text, as
    CommonMark has it.
    """
    m = FENCE.match(line)
    if m and not m.group(2).strip(" \t"):
        return m.group(1)[0], len(m.group(1))
    return None


def _fence_closes(line, fence):
    """True if `line` closes `fence`: the same character, at least as long, and nothing after it
    but spaces or tabs — a non-breaking space is text, as CommonMark has it.

    A shorter or different fence inside is text, which is how one fence nests inside another.
    """
    run = _bare_run(line)
    return bool(run and run[0] == fence[0] and run[1] >= fence[1])


def _fenced_flags(lines):
    """For each line, whether it lies inside a fenced block (its delimiters included).

    An opener that is never closed is read as ordinary text, so fences that closed earlier stay
    excluded and the lines after the stray opener are read normally.  CommonMark would instead run
    an unclosed fence to the end of the file, swallowing a devlog's body into a code block and
    hiding every field after the opener.

    One pass: `longest[c][i]` is the longest bare run of fence character `c` at or after line
    `i`, so an opener of `c` with length L has a closer somewhere after it iff
    `longest[c][i + 1] >= L`, and an opener without one is text at once.
    """
    n = len(lines)
    longest = {"`": [0] * (n + 1), "~": [0] * (n + 1)}
    for i in range(n - 1, -1, -1):
        for c in longest:
            longest[c][i] = longest[c][i + 1]
        if (run := _bare_run(lines[i])) is not None:
            longest[run[0]][i] = max(longest[run[0]][i], run[1])
    flags, fence = [], None
    for i, line in enumerate(lines):
        if fence is not None:
            flags.append(True)
            if _fence_closes(line, fence):
                fence = None
            continue
        opened = _fence_open(line)
        if opened is not None and longest[opened[0]][i + 1] >= opened[1]:
            fence = opened
            flags.append(True)
            continue
        flags.append(False)
    return flags


def _header_lines(path):
    """The unfenced lines before the first `## ` heading outside a fenced block.

    Fenced lines are excluded, so an example `**Date:**` or model id cannot stand in for a real
    one.  Fences are resolved by `_fenced_flags`.
    """
    lines = _lines(path)
    out = []
    for line, fenced in zip(lines, _fenced_flags(lines)):
        if fenced:
            continue
        if line.startswith("## "):
            return out
        out.append(line)
    return out


def devlog_fields(paths):
    """Counter of header field names — one count per file — over files that have a header."""
    counts, denom = Counter(), 0
    for p in paths:
        fields = {m.group(1) for l in _header_lines(p) if (m := FIELD.match(l))}
        if fields:
            denom += 1
            counts.update(fields)
    return counts, denom


def devlogs_without_date(paths):
    """Devlogs whose header carries no **Date:** — the one field the protocol requires."""
    return sorted(p.name for p in paths
                  if not any((m := FIELD.match(l)) and m.group(1) == "Date" for l in _header_lines(p)))


def _open_fence_at(lines, target):
    """The index of the opener whose fence is still open at `target`, scanning strictly, or None.

    Strict means CommonMark's reading: a fence opened above runs until a closer, and an unclosed one
    runs to the end of the file.  `_fenced_flags` instead reads an unclosed opener as text.
    """
    fence, opener = None, None
    for i in range(target):
        if fence is not None:
            if _fence_closes(lines[i], fence):
                fence = None
            continue
        if (opened := _fence_open(lines[i])) is not None:
            fence, opener = opened, i
    return opener if fence is not None else None


def _depth_span(proto_path):
    """The depth-selection section, from its heading to the next heading, both outside fences.

    Fences are resolved over the whole protocol file, so a heading inside a fenced example is
    text and neither starts nor ends the section.  The heading must occur exactly once outside
    fences: a stray unclosed opener is read as text, which can expose a quoted copy of the heading,
    and with two candidates the protocol is ambiguous, so the run fails closed rather than pick
    one.  A single heading that a strict reading would still place inside an unclosed fence fails
    closed too.
    """
    lines = _lines(proto_path)
    fenced = _fenced_flags(lines)
    starts = [i for i, l in enumerate(lines) if not fenced[i] and DEPTH_HEADING.match(l)]
    if not starts:
        raise ExtractionError(f"anchor not found: {DEPTH_ANCHOR}")
    if len(starts) > 1:
        raise ExtractionError(f"depth anchor {DEPTH_ANCHOR!r} appears {len(starts)} times outside "
                              f"fenced blocks, at lines {', '.join(str(i + 1) for i in starts)}")
    start = starts[0]
    # The resolver reads a stray opener as text; for the depth section that would let a heading
    # inside an unclosed fence stand, so an unclosed fence above the heading fails the run closed.
    if (opener := _open_fence_at(lines, start)) is not None:
        raise ExtractionError(f"depth anchor at line {start + 1} lies inside a fence opened at line "
                              f"{opener + 1} and never closed")
    end = next((i for i in range(start + 1, len(lines))
                if not fenced[i] and (lines[i].startswith("### ") or lines[i].startswith("## "))),
               len(lines))
    return lines[start:end]


SEAT_COLUMNS = {"codex": "Codex", "grok": "Grok"}


def _cells(row):
    """A table row's cells.  Exactly one edge pipe is removed on each side, so an empty first
    cell (`|| a | b |`) stays a cell instead of shifting the rest left."""
    s = row.strip()
    if s.startswith("|"):
        s = s[1:]
    if s.endswith("|"):
        s = s[:-1]
    return [c.strip() for c in s.split("|")]


def _table_rows(proto_path):
    """The first contiguous block of `|` lines in the depth-selection span, outside fences.

    A fence opened here that never closes is deliberately NOT demoted to text, unlike in
    `_header_lines`: demoting it would expose any example table inside it as the depth table.
    Left open, it hides the rest of the span, no rows are found, and the required-extract guard
    fails the run closed.
    """
    rows, fence = [], None
    for line in _depth_span(proto_path):
        if fence is not None:
            if _fence_closes(line, fence):
                fence = None
            continue
        if (opened := _fence_open(line)) is not None:
            if rows:
                break
            fence = opened
            continue
        if line.startswith("|"):
            rows.append(line)
        elif rows:
            break
    return rows


def table_models(proto_path):
    """Models named by the depth table, per seat column: {seat: set of models}.

    The column is found by the table's own header cell, never inferred from a model's prefix —
    a gpt- id written in the Grok cell is a Grok-column entry, and the probes must see it there.
    A header that does not name each seat exactly once, or a row whose width differs from the
    header's, raises ExtractionError.
    """
    rows = _table_rows(proto_path)
    if not rows:
        return {}
    header = _cells(rows[0])
    # A table the probes cannot read unambiguously fails closed rather than being half-read.
    columns = {}
    for seat, name in SEAT_COLUMNS.items():
        if header.count(name) != 1:
            raise ExtractionError(f"depth table header must name {name!r} exactly once, "
                                  f"found {header.count(name)}")
        columns[seat] = header.index(name)
    for k, row in enumerate(rows[1:], 2):
        if len(_cells(row)) != len(header):
            raise ExtractionError(f"depth table row {k} has {len(_cells(row))} cells, "
                                  f"the header has {len(header)}")
    out = {seat: set() for seat in columns}
    for row in rows[1:]:
        cells = _cells(row)
        for seat, i in columns.items():
            out[seat] |= set(MODEL.findall(cells[i]))
    return out


def exception_models(proto_path):
    """Models declared deliberately unnamed, from the anchored line only.

    May legitimately be empty — a protocol with no escape hatch is valid — so this extract is
    never guarded non-empty.
    """
    return {m for l in _depth_span(proto_path)
            if l.startswith(EXCEPTION_ANCHOR) for m in MODEL.findall(l)}


def wrapper_models(sh_path):
    """The allowlist the wrapper enforces, read from its `case` arm.

    The arm is what accepts or rejects a model; the rejection message is only a description of
    it.  The effort arm (low|medium|high|xhigh) contributes nothing, being outside the alphabet.
    """
    found = set()
    for l in _lines(sh_path):
        if m := CASE_ARM.match(l):
            found |= {t for t in m.group(1).split("|") if MODEL.fullmatch(t)}
    return found


def wrapper_message_models(sh_path):
    """The models the wrapper's rejection message claims to accept."""
    return {m for l in _lines(sh_path) if "Allowed:" in l for m in MODEL.findall(l)}


def wrapper_default(sh_path, variable):
    """The model assigned to `variable` when the caller sets none.

    Anchored on the variable name, so an unrelated `*_MODEL` guard elsewhere in the file cannot be
    mistaken for the seat's default.
    """
    guard = re.compile(r'if \[ -z "\$\{' + re.escape(variable) + r':-\}" \];\s*then\s*'
                       + re.escape(variable) + r'="([^"]+)"')
    for l in _lines(sh_path):
        if m := guard.search(l):
            return m.group(1)
    return None


def record_models(devlog_paths, repo_path, n_commits, commit_texts=None):
    """Models named in recent devlog headers and recent commit bodies, one count per record.

    Returns (headers, commits, commit_status).  Headers only on the devlog side: the same token
    in body prose ("grok-4.6 stays allowlisted for reproducing older rounds") is a mention.
    Commit bodies corroborate a model already seen in a header and never trigger a finding alone,
    because a body can name a model while describing a rule rather than what audited that commit.

    A commit read that fails is reported in `commit_status`, never collapsed into an empty set that
    the report would describe as read.
    """
    headers = Counter()
    for p in devlog_paths:
        headers.update({m for l in _header_lines(p) for m in MODEL.findall(l)})

    status = f"{n_commits} commit bodies"
    if commit_texts is None:
        if n_commits <= 0:
            commit_texts, status = [], "no commit bodies requested"
        else:
            try:
                out = subprocess.run(
                    ["git", "-C", str(repo_path), "log", f"-{n_commits}", "--format=%B%x00"],
                    capture_output=True, check=True, encoding="utf-8", errors="replace").stdout
                commit_texts = [t for t in out.split("\0") if t.strip()]
                status = f"{len(commit_texts)} commit bodies"
            except (subprocess.CalledProcessError, OSError) as exc:
                commit_texts = []
                status = f"commit bodies UNAVAILABLE ({exc.__class__.__name__}); corroboration off"
    commits = Counter()
    for body in commit_texts:
        commits.update(set(MODEL.findall(body)))
    return headers, commits, status


def index_targets(readme_path):
    """Filenames that are the link destination of an index entry.

    An index line may mention a second filename in its prose, an entry's label may differ from
    its destination, and a destination may point into another directory; matching any of those
    would count a file as indexed when it is not.
    """
    # Lexical, not Path.resolve(): resolve() follows a symlinked devlog out of the directory and
    # can raise on a malformed destination, and one bad index line must not end the check.
    base = os.path.normpath(os.path.dirname(os.path.abspath(readme_path)))
    out = set()
    for l in _lines(readme_path):
        if not (m := INDEX_ENTRY.match(l)):
            continue
        target = m.group(1).split("#", 1)[0]
        if not target:
            continue
        # A link into another directory, or a URL ending in the same name, lands outside the
        # README's own directory and indexes a different file.
        joined = os.path.normpath(os.path.join(base, target))
        if os.path.dirname(joined) == base:
            out.add(os.path.basename(joined))
    return out


def devlog_files(devlog_dir):
    return {p.name for p in Path(devlog_dir).glob("dl*.md")}


def title_forms(paths):
    old = sum(1 for p in paths if OLD_TITLE.match(next(iter(_lines(p)), "")))
    return old, len(paths)


def head_commit(repo_path):
    try:
        return subprocess.run(["git", "-C", str(repo_path), "rev-parse", "--short", "HEAD"],
                              capture_output=True, check=True, encoding="utf-8",
                              errors="replace").stdout.strip()
    except (subprocess.CalledProcessError, OSError):
        return "unknown"


# --------------------------------------------------------------------------- guard

def guard_required(extracts):
    """Fail on an empty required extract; never skip.

    A probe whose failure mode is silently reporting success is worse than no probe: an empty
    extract makes every comparison that consumes it pass vacuously.  The by-design exception set
    is deliberately absent from this mapping — empty means "no escape hatch", which is a valid
    document and a clean self-test fixture.
    """
    empty = sorted(name for name, value in extracts.items() if not value)
    if empty:
        raise ExtractionError("required extract empty: " + ", ".join(empty))


# --------------------------------------------------------------------------- probes

def _seat(model):
    return "codex" if model.startswith("gpt-") else "grok"


def _newest_grok(allowlist):
    """The greatest plain-versioned grok; suffixed ids such as grok-4.7-build-fast are not seats."""
    versions = {m: tuple(int(x) for x in v.group(1).split("."))
                for m in allowlist if (v := GROK_VERSION.fullmatch(m))}
    return max(versions, key=versions.get) if versions else None


def probe_models(columns, exceptions, allowlists, messages, defaults, headers, commits):
    """P2a-P2e.  `allowlists`, `messages` and `defaults` are keyed by seat."""
    findings = []
    reported = set()

    # The Grok seat is always the latest allowlisted grok.  That model is declared by the policy
    # itself, whether or not a row happens to name it.
    policy_grok = _newest_grok(allowlists.get("grok", set()))
    rows = set().union(*columns.values()) if columns else set()
    codex_rows = columns.get("codex", set())
    grok_rows = columns.get("grok", set())
    declared = rows | exceptions | ({policy_grok} if policy_grok else set())
    # Older groks stay allowlisted so an earlier round can be reproduced, and a devlog recording a
    # round run while one of them was the latest is history, not drift.  The default and the rows
    # are held to the latest below; only use is relaxed.
    declared_for_use = declared | allowlists.get("grok", set())

    # P2c - a default that is not what the policy names.  "Equals the latest", not "appears in
    # some row": a default left on an older row would otherwise pass while the policy moved on.
    if (d := defaults.get("grok")) and policy_grok and d != policy_grok:
        findings.append((DRIFT, "P2c", d,
                         f"grok default is {d}; the policy names the latest allowlisted grok, "
                         f"{policy_grok}"))
        reported.add(d)
    # A non-grok token in the Grok column is misplaced, which P2b reports; here only groks.
    for model in sorted(m for m in grok_rows if m.startswith("grok-")):
        if policy_grok and model != policy_grok:
            findings.append((DRIFT, "P2c", model,
                             f"the Grok column names {model}; the policy names the latest "
                             f"allowlisted grok, {policy_grok}"))
            reported.add(model)
    # The Codex default must be a model the Codex arm accepts and the table's Codex column names;
    # the same id written only in the Grok column satisfies neither.
    if d := defaults.get("codex"):
        if d not in allowlists.get("codex", set()):
            findings.append((DRIFT, "P2c", d, f"codex default {d} is rejected by ask_codex.sh's own arm"))
            reported.add(d)
        elif d not in codex_rows:
            findings.append((DRIFT, "P2c", d, f"codex default {d} is named by no Codex-column row"))
            reported.add(d)

    # P2a - used in recent records, declared nowhere.  Headers trigger; commit bodies only
    # corroborate.  A model P2c already diagnoses is not reported twice.
    for model in sorted(headers):
        if model in declared_for_use or model in reported:
            continue
        also = f", and in {commits[model]} commit bodies" if commits.get(model) else ""
        findings.append((DRIFT, "P2a", model,
                         f"{model} named in {headers[model]} recent devlog headers{also}; "
                         "no depth row names it and it is not declared by design"))
        reported.add(model)

    # P2b - a row naming a model its own seat's wrapper would reject.  Bound to the seat: every
    # gpt row is absent from the grok allowlist and vice versa.
    # Each column is checked against its own seat's arm, so a model in the wrong column is caught.
    for seat, models in sorted(columns.items()):
        for model in sorted(models):
            if seat in allowlists and model not in allowlists[seat]:
                findings.append((DRIFT, "P2b", model,
                                 f"the {SEAT_COLUMNS[seat]} column names {model}, "
                                 f"which ask_{seat}.sh would reject"))
    # The by-design line says "allowlisted"; an id no wrapper accepts makes that line false.  The
    # line has no columns, so its seat is read from the id's prefix.
    for model in sorted(exceptions):
        seat = _seat(model)
        if seat in allowlists and model not in allowlists[seat]:
            findings.append((DRIFT, "P2b", model,
                             f"the by-design line names {model}, which ask_{seat}.sh would reject"))
            reported.add(model)

    # P2e - the rejection message lists a different set from the arm that enforces it, so a
    # caller told which models are allowed is told something false.
    for seat in sorted(allowlists):
        arm, said = allowlists[seat], messages.get(seat, set())
        if said and said != arm:
            diff = ", ".join(sorted(arm ^ said))
            findings.append((DRIFT, "P2e", seat,
                             f"ask_{seat}.sh's 'Allowed:' message disagrees with its case arm "
                             f"on {diff}"))

    # P2d - allowlisted, unruled, not an escape hatch.  Informational: a wrapper may carry a
    # model no row needs.  A model already reported as DRIFT is not repeated here.
    for seat, allow in sorted(allowlists.items()):
        for model in sorted(allow - declared - reported):
            findings.append((INFO, "P2d", model,
                             f"{model} allowlisted in ask_{seat}.sh, named by no depth row"))
    return findings


def probe_index(files, targets):
    """P3 - a devlog outside its own index."""
    return [(DRIFT, "P3", name, f"{name} has no index entry in devlog/README.md")
            for name in sorted(files - targets)]


def probe_template(declared, used, window_n, threshold, headed_n):
    """P1a (DRIFT) and P1b (INFO).

    The directions are not symmetric.  A declared field nobody writes means the template
    describes something dead: a fact about the document, so it fires at zero uses with no
    threshold.  A field practice uses but the template omits is not a defect — the template states
    that its list is not a schema — so it is reported only once it is common (`threshold`, in
    percent of the window), and only as INFO.
    """
    # With no headed entry in the window, "used by none" says nothing about the template; P1c
    # reports each missing Date instead.
    findings = []
    if headed_n:
        findings += [(DRIFT, "P1a", f,
                      f"template declares **{f}:** — used by none of the last {window_n} entries")
                     for f in sorted(declared) if not used.get(f)]
    findings += [(INFO, "P1b", f,
                  f"**{f}:** used by {used[f]}/{window_n} recent entries, not declared by the template")
                 for f in sorted(used) if f not in declared and 100 * used[f] / window_n >= threshold]
    return findings


def probe_dates(missing):
    """P1c - a recent devlog without **Date:**, the one field the protocol requires."""
    return [(DRIFT, "P1c", name, f"{name} has no **Date:** field") for name in missing]


def probe_titles(old, total):
    """P4 - title form.  Informational: the older form is explicitly not to be churned."""
    if not old:
        return []
    return [(INFO, "P4", "title",
             f"{old}/{total} recent entries use the retired '# Devlog <date>' title form")]


# --------------------------------------------------------------------------- live run + report

def run_live(root, recent, commits_n, threshold, commit_texts=None):
    root = Path(root)
    proto = root / "doc/ai_collaboration_protocol.md"
    devlog_dir = root / "devlog"
    readme = devlog_dir / "README.md"

    window = devlog_window(devlog_dir, recent)
    declared = template_fields(proto)
    used, denom = devlog_fields(window)
    columns = table_models(proto)
    rows = set().union(*columns.values()) if columns else set()
    exceptions = exception_models(proto)
    allowlists = {s: wrapper_models(root / p) for s, (p, _) in WRAPPERS.items()}
    messages = {s: wrapper_message_models(root / p) for s, (p, _) in WRAPPERS.items()}
    defaults = {s: wrapper_default(root / p, v) for s, (p, v) in WRAPPERS.items()}
    headers, commit_hits, commit_status = record_models(window, root, commits_n, commit_texts)
    files, targets = devlog_files(devlog_dir), index_targets(readme)
    old, total = title_forms(window)
    missing_date = devlogs_without_date(window)

    guard_required({
        "template fields": declared,
        # The window must be non-empty; its headers need not be — a headerless entry is P1c.
        "devlog window": window,
        "depth table rows": rows,
        "Codex column models": columns.get("codex"),
        "codex allowlist": allowlists["codex"],
        "grok allowlist": allowlists["grok"],
        "codex default": defaults["codex"],
        "grok default": defaults["grok"],
        # With no plain-versioned grok the policy names nothing, and P2c would skip silently.
        "policy grok": _newest_grok(allowlists["grok"]),
        "devlog files": files,
        "index targets": targets,
    })

    findings = probe_models(columns, exceptions, allowlists, messages, defaults, headers, commit_hits)
    findings += probe_index(files, targets)
    findings += probe_template(declared, used, len(window), threshold, denom)
    findings += probe_dates(missing_date)
    findings += probe_titles(old, total)

    window_note = (f"at {head_commit(root)}: {denom} of the newest {len(window)} devlogs "
                   f"(newest date first, then name) carry a header, back to "
                   f"{window[-1].name if window else '-'}; {len(files)} files, "
                   f"{len(targets)} index entries; {commit_status}")
    return findings, window_note


def report(findings, window_note):
    print(f"window: {window_note}\n")
    for tier in (DRIFT, INFO):
        group = [f for f in findings if f[0] == tier]
        if group:
            print(f"{tier} ({len(group)}):")
            for _, probe, _, message in group:
                print(f"  {probe:4s} {message}")
            print()
    if not findings:
        print("no drift found")
    return sum(1 for f in findings if f[0] == DRIFT)


def run(root, recent=40, commits=40, threshold=25.0, strict=False, commit_texts=None, quiet=False):
    """One live run; returns the process exit code."""
    try:
        findings, window_note = run_live(root, recent, commits, threshold, commit_texts)
    except (ExtractionError, OSError) as exc:
        if not quiet:
            print(f"extraction failed: {exc}", file=sys.stderr)
        return 2
    if quiet:
        drift = sum(1 for f in findings if f[0] == DRIFT)
    else:
        drift = report(findings, window_note)
    return 1 if (strict and drift) else 0


# --------------------------------------------------------------------------- self-test
# Every probe has a fixture that must make it fire and one that must keep it silent.  Each case
# is written into a fresh directory, so no file planted for one case is visible to the next.

_ROW_GROK = "| Plan audit | `gpt-5.6-terra`, `high` | `grok-4.7`, `high` |"

_PROTO = """### Required Header

```markdown
# A topic phrase

**Date:** 2026-01-01
**Purpose:** what it was for
**AIs involved:** Claude Opus 5
```

### When to Write a Devlog

### 9.1 Depth Selection

{exceptions}| Subject | Codex | Grok |
|---|---|---|
{row}

### Workflow
"""

_DEVLOG = """{title}

**Date:** 2026-01-01
**Purpose:** what it was for
**AIs involved:** Claude Opus 5, Codex `gpt-5.6-terra`/high, Grok `grok-4.7`/high{extra}

## Summary
"""

_README = "## Entries\n- [dl20260101_clean.md](dl20260101_clean.md) — a clean entry\n"


def _sh(variable, default, arm, message=None):
    return ('if [ -z "${' + variable + ':-}" ]; then ' + variable + '="' + default + '"; fi\n'
            'case "$' + variable + '" in\n'
            '    ' + "|".join(arm) + ') ;;\n'
            '    *) echo "Allowed: ' + " ".join(message or arm) + '" >&2 ;;\n'
            'esac\n'
            'case "$EFFORT" in\n'
            '    low|medium|high|xhigh) ;;\n'
            '    *) echo "Allowed: low medium high xhigh" >&2 ;;\n'
            'esac\n')


def _tree(base, name, *, row=_ROW_GROK, exceptions="", template_extra="",
          title="# A topic phrase", devlog_extra="", readme=_README,
          codex=("gpt-5.6-terra", ["gpt-5.6-terra"], None),
          grok=("grok-4.7", ["grok-4.7"], None), extra_devlogs=()):
    root = Path(base) / name
    for sub in ("doc", "devlog", ".claude/scripts"):
        (root / sub).mkdir(parents=True)
    proto = _PROTO.format(exceptions=exceptions, row=row)
    if template_extra:
        proto = proto.replace("**AIs involved:** Claude Opus 5\n",
                              "**AIs involved:** Claude Opus 5\n" + template_extra + "\n")
    (root / "doc/ai_collaboration_protocol.md").write_text(proto)
    (root / ".claude/scripts/ask_codex.sh").write_text(_sh("CODEX_MODEL", *codex))
    (root / ".claude/scripts/ask_grok.sh").write_text(_sh("GROK_MODEL", *grok))
    (root / "devlog/dl20260101_clean.md").write_text(_DEVLOG.format(title=title, extra=devlog_extra))
    (root / "devlog/README.md").write_text(readme)
    for fname, text in extra_devlogs:
        (root / "devlog" / fname).write_text(text)
    return root


def self_test():
    failures = []

    def check(label, condition):
        print(f"  {'ok  ' if condition else 'FAIL'} {label}")
        if not condition:
            failures.append(label)

    def find(root, commit_texts=None):
        return run_live(root, 40, 0, 25, commit_texts)[0]

    def has(findings, probe, subject=None):
        return any(p == probe and (subject is None or s == subject) for _, p, s, _ in findings)

    plain = _DEVLOG.format(title="# A topic phrase", extra="")

    with tempfile.TemporaryDirectory() as base:
        print("clean tree:")
        f = find(_tree(base, "clean"))
        check("no findings at all", not f)
        try:
            run_live(_tree(base, "no_exceptions"), 40, 0, 25)
            check("an empty by-design list is accepted by the guard", True)
        except ExtractionError as exc:
            check(f"an empty by-design list is accepted by the guard ({exc})", False)

        print("P1 - template fields:")
        f = find(_tree(base, "p1a", template_extra="**Retired:** nobody writes this"))
        check("P1a fires on a declared field nobody uses", has(f, "P1a", "Retired"))
        f = find(_tree(base, "p1b", devlog_extra="\n**Module:** meta"))
        check("P1b reports an undeclared practice field", has(f, "P1b", "Module"))
        check("P1b is INFO, never DRIFT", all(t == INFO for t, p, _, _ in f if p == "P1b"))

        print("P2a - used, declared nowhere:")
        rogue = "\n**Voices:** Codex `gpt-9.9-rogue`, then `gpt-9.9-rogue` again"
        f = find(_tree(base, "p2a", devlog_extra=rogue))
        check("P2a fires on a model named by no row", has(f, "P2a", "gpt-9.9-rogue"))
        msg = next((m for _, p, s, m in f if p == "P2a"), "")
        check("P2a counts headers, not occurrences", "named in 1 recent devlog headers" in msg)
        f = find(_tree(base, "p2a_commit", devlog_extra=rogue),
                 commit_texts=["gpt-9.9-rogue gpt-9.9-rogue gpt-9.9-rogue"])
        msg = next((m for _, p, s, m in f if p == "P2a"), "")
        check("P2a counts commit bodies, not occurrences", "and in 1 commit bodies" in msg)
        f = find(_tree(base, "p2a_nohdr"), commit_texts=["reviewed by gpt-9.9-rogue"])
        check("a commit body alone never triggers P2a", not has(f, "P2a"))
        f = find(_tree(base, "p2a_suffix", devlog_extra="\n**Voices:** `grok-4.7-build-fast`"))
        check("a suffixed id is its own model, not the policy grok",
              has(f, "P2a", "grok-4.7-build-fast"))

        print("P2b - a row its own wrapper would reject:")
        f = find(_tree(base, "p2b", row="| Plan audit | `gpt-5.6-terra`, `gpt-5.6-luna` | `grok-4.7`, `high` |"))
        check("P2b fires on a codex row the codex arm rejects", has(f, "P2b", "gpt-5.6-luna"))
        f = find(_tree(base, "p2b_healthy"))
        check("P2b is silent on a healthy two-seat table", not has(f, "P2b"))

        print("P2c - a default off the policy:")
        f = find(_tree(base, "p2c_grok", grok=("grok-4.6", ["grok-4.7", "grok-4.6"], None)))
        check("P2c fires on a lagging grok default", has(f, "P2c", "grok-4.6"))
        f = find(_tree(base, "p2c_codex",
                       codex=("gpt-5.6-luna", ["gpt-5.6-terra", "gpt-5.6-luna"], None)))
        check("P2c fires on a codex default named by no row", has(f, "P2c", "gpt-5.6-luna"))
        f = find(_tree(base, "p2c_prec", grok=("grok-4.6", ["grok-4.7", "grok-4.6"], None),
                       devlog_extra="\n**Voices:** `grok-4.6`"))
        hits = [p for _, p, s, _ in f if s == "grok-4.6"]
        check("one defect is one finding, and P2c is the one reporting it", hits == ["P2c"])

        print("P2d - allowlisted, unruled:")
        f = find(_tree(base, "p2d", codex=("gpt-5.6-terra", ["gpt-5.6-terra", "gpt-5.5"], None)))
        check("P2d reports an unruled allowlisted model", has(f, "P2d", "gpt-5.5"))
        check("P2d is INFO", all(t == INFO for t, p, _, _ in f if p == "P2d"))
        f = find(_tree(base, "p2d_policy", row="| Plan audit | `gpt-5.6-terra`, `high` | `high` |"))
        check("the policy grok is declared even when no row names it",
              not any(s == "grok-4.7" for _, _, s, _ in f))

        print("by-design exceptions:")
        f = find(_tree(base, "excn", exceptions=EXCEPTION_ANCHOR + " `gpt-5.6-sol`\n\n",
                       codex=("gpt-5.6-terra", ["gpt-5.6-terra", "gpt-5.6-sol"], None),
                       devlog_extra="\n**Voices:** `gpt-5.6-sol`"))
        check("a declared exception is not reported, even when used",
              not any(s == "gpt-5.6-sol" for _, _, s, _ in f))

        print("P2e - message and arm disagree:")
        f = find(_tree(base, "p2e", codex=("gpt-5.6-terra", ["gpt-5.6-terra"],
                                           ["gpt-5.6-terra", "gpt-5.5"])))
        check("P2e fires when the message lists a model the arm rejects", has(f, "P2e", "codex"))
        check("the allowlist is read from the arm, not the message", not has(f, "P2d", "gpt-5.5"))

        print("P3 - index coverage:")
        f = find(_tree(base, "p3", extra_devlogs=[("dl20260102_unindexed.md", plain)]))
        check("P3 fires on an unindexed devlog", has(f, "P3", "dl20260102_unindexed.md"))
        f = find(_tree(base, "p3_mention", extra_devlogs=[("dl20260102_unindexed.md", plain)],
                       readme="## Entries\n- [dl20260101_clean.md](dl20260101_clean.md) — see "
                              "`dl20260102_unindexed.md`\n"))
        check("a filename mentioned in an entry's prose is not indexed",
              has(f, "P3", "dl20260102_unindexed.md"))
        f = find(_tree(base, "p3_label", readme="## Entries\n- [short name](dl20260101_clean.md)\n"))
        check("the link destination counts, not the label", not has(f, "P3"))
        f = find(_tree(base, "p3_heading",
                       readme="## Entries\n### [dl20260101_clean.md](dl20260101_clean.md) — x\n"))
        check("a heading-form entry indexes its file", not has(f, "P3"))

        print("P4 - title form:")
        f = find(_tree(base, "p4", title="# Devlog 2026-01-01 — a topic"))
        check("P4 reports a retired title form", has(f, "P4"))
        check("P4 is INFO", all(t == INFO for t, p, _, _ in f if p == "P4"))

        print("extraction:")
        check("a trailing sentence period is not part of the model",
              MODEL.findall("defaults to grok-4.7.") == ["grok-4.7"])
        check("seat names are outside the model alphabet",
              MODEL.findall("Claude Opus 5 and Claude Fable 5.1") == [])
        check("a suffixed grok id does not break version ordering",
              _newest_grok({"grok-4.7-build-fast", "grok-4.6"}) == "grok-4.6")
        try:
            find(_tree(base, "guard", row=""))
            check("the guard rejects an empty required extract", False)
        except ExtractionError:
            check("the guard rejects an empty required extract", True)
        _, _, status = record_models([], Path(base) / "not-a-repo", 5)
        check("a failed commit read is reported, not described as read", "UNAVAILABLE" in status)

        print("v3 - declarations the wrappers must honour:")
        f = find(_tree(base, "excn_rejected", exceptions=EXCEPTION_ANCHOR + " `gpt-9.9-ghost`\n\n"))
        check("P2b fires on a by-design id no wrapper accepts", has(f, "P2b", "gpt-9.9-ghost"))
        f = find(_tree(base, "grok_row_lag",
                       row="| Plan audit | `gpt-5.6-terra`, `high` | `grok-4.6`, `high` |",
                       grok=("grok-4.7", ["grok-4.7", "grok-4.6"], None)))
        check("P2c fires on a Grok row naming an older grok", has(f, "P2c", "grok-4.6"))
        f = find(_tree(base, "grok_row_ok"))
        check("P2c is silent when the Grok row names the latest", not has(f, "P2c"))
        f = find(_tree(base, "codex_default_rejected",
                       row="| Plan audit | `gpt-5.6-terra`, `gpt-5.6-luna` | `grok-4.7`, `high` |",
                       codex=("gpt-5.6-luna", ["gpt-5.6-terra"], None)))
        check("P2c fires on a Codex default its own arm rejects, even when a row names it",
              has(f, "P2c", "gpt-5.6-luna"))
        f = find(_tree(base, "codex_default_grok_column",
                       codex=("grok-4.7", ["gpt-5.6-terra", "grok-4.7"], None)))
        check("a Grok-column token does not satisfy the Codex default", has(f, "P2c", "grok-4.7"))
        f = find(_tree(base, "codex_prec", codex=("gpt-5.6-luna", ["gpt-5.6-terra", "gpt-5.6-luna"], None),
                       devlog_extra="\n**Voices:** `gpt-5.6-luna`"))
        hits = [p for _, p, s, _ in f if s == "gpt-5.6-luna"]
        check("a Codex default off the table is one P2c finding, not also P2a", hits == ["P2c"])

        print("v3 - history is not drift:")
        f = find(_tree(base, "old_grok_use", grok=("grok-4.7", ["grok-4.7", "grok-4.6"], None),
                       devlog_extra="\n**Voices:** `grok-4.6`"))
        check("an older allowlisted grok recorded in a header is not P2a", not has(f, "P2a", "grok-4.6"))
        f = find(_tree(base, "unlisted_grok", devlog_extra="\n**Voices:** `grok-9.9`"))
        check("a grok no wrapper allowlists is still P2a", has(f, "P2a", "grok-9.9"))

        print("v3 - index destinations:")
        f = find(_tree(base, "p3_elsewhere", extra_devlogs=[("dl20260102_other.md", plain)], readme="## Entries\n- [dl20260102_other.md](dl20260102_other.md)\n- [x](../archive/dl20260101_clean.md)\n"))
        check("a link into another directory does not index the local file",
              has(f, "P3", "dl20260101_clean.md"))
        f = find(_tree(base, "p3_url", extra_devlogs=[("dl20260102_other.md", plain)], readme="## Entries\n- [dl20260102_other.md](dl20260102_other.md)\n- [x](https://example.org/dl20260101_clean.md)\n"))
        check("a URL ending in the same name does not index it", has(f, "P3", "dl20260101_clean.md"))
        f = find(_tree(base, "p3_dotslash", readme="## Entries\n- [x](./dl20260101_clean.md)\n"))
        check("a ./ link to the local file indexes it", not has(f, "P3"))

        print("v3 - the required Date field:")
        nodate = plain.replace("**Date:** 2026-01-01\n", "")
        f = find(_tree(base, "p1c", extra_devlogs=[("dl20260102_nodate.md", nodate)],
                       readme=_README + "- [dl20260102_nodate.md](dl20260102_nodate.md) — x\n"))
        check("P1c fires on one missing Date while other entries carry the field",
              has(f, "P1c", "dl20260102_nodate.md"))
        f = find(_tree(base, "p1c_clean"))
        check("P1c is silent when every entry has a Date", not has(f, "P1c"))

        print("v3 - windows and denominators:")
        wdir = Path(base) / "window" / "devlog"
        wdir.mkdir(parents=True)
        for n in ("dl20260101_z.md", "dl20260101_a.md", "dl20260102_m.md"):
            (wdir / n).write_text(plain)
        check("the window is newest date first, then name ascending",
              [p.name for p in devlog_window(wdir, 3)] ==
              ["dl20260102_m.md", "dl20260101_a.md", "dl20260101_z.md"])
        headless = "# A topic phrase\n\n## Summary\n"
        root = _tree(base, "p1b_window", devlog_extra="\n**Module:** meta", extra_devlogs=[("dl20260102_headless.md", headless)],
                     readme=_README + "- [dl20260102_headless.md](dl20260102_headless.md) — x\n")
        f = run_live(root, 40, 0, 60)[0]
        check("P1b's share is of the window, not of headed files", not has(f, "P1b", "Module"))

        print("v3 - failures are loud:")
        gone = _tree(base, "missing_proto")
        (gone / "doc/ai_collaboration_protocol.md").unlink()
        err = io.StringIO()
        with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
            rc = run(gone, commits=0)
        check("a missing input is exit 2 with a message, not a traceback",
              rc == 2 and "extraction failed" in err.getvalue())
        try:
            find(_tree(base, "suffix_only", grok=("grok-4.7-build-fast", ["grok-4.7-build-fast"], None)))
            check("an allowlist with no plain grok fails the guard", False)
        except ExtractionError:
            check("an allowlist with no plain grok fails the guard", True)
        bad = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--recent", "-1"],
                             capture_output=True, text=True)
        check("--recent -1 is rejected by the parser", bad.returncode == 2 and "must be" in bad.stderr)
        bad = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--threshold", "nan"],
                             capture_output=True, text=True)
        check("--threshold nan is rejected by the parser", bad.returncode == 2 and "must be" in bad.stderr)

        print("v4 - columns are read from the table, not from a prefix:")
        f = find(_tree(base, "col_grok_cell",
                       row="| Plan audit | `gpt-5.6-terra`, `high` | `grok-4.7`, `gpt-5.6-luna` |",
                       codex=("gpt-5.6-luna", ["gpt-5.6-terra", "gpt-5.6-luna"], None)))
        check("a Codex default named only in the Grok column is P2c", has(f, "P2c", "gpt-5.6-luna"))
        check("a gpt id in the Grok column is P2b against the Grok arm", has(f, "P2b", "gpt-5.6-luna"))
        f = find(_tree(base, "col_codex_cell",
                       row="| Plan audit | `gpt-5.6-terra`, `gpt-5.6-luna` | `grok-4.7`, `high` |",
                       codex=("gpt-5.6-luna", ["gpt-5.6-terra", "gpt-5.6-luna"], None)))
        check("the same default named in the Codex column is silent",
              not any(s == "gpt-5.6-luna" for _, _, s, _ in f))
        f = find(_tree(base, "col_grok_in_codex",
                       row="| Plan audit | `gpt-5.6-terra`, `grok-4.7` | `grok-4.7`, `high` |"))
        check("a grok id in the Codex column is P2b against the Codex arm", has(f, "P2b", "grok-4.7"))
        noheader = _PROTO.replace("| Subject | Codex | Grok |", "| Subject | Codex | Other |")
        root = _tree(base, "col_no_grok_header")
        (root / "doc/ai_collaboration_protocol.md").write_text(
            noheader.format(exceptions="", row=_ROW_GROK))
        try:
            run_live(root, 40, 0, 25)
            check("a table without a Grok column header fails the guard", False)
        except ExtractionError:
            check("a table without a Grok column header fails the guard", True)

        print("v4 - index destinations are lexical:")
        root = _tree(base, "symlinked")
        (root / "elsewhere").mkdir()
        (root / "elsewhere/real.md").write_text(plain)
        os.symlink(root / "elsewhere/real.md", root / "devlog/dl20260102_link.md")
        (root / "devlog/README.md").write_text(_README + "- [dl20260102_link.md](dl20260102_link.md)\n")
        f = run_live(root, 40, 0, 25)[0]
        check("a symlinked devlog indexed under its own name is not P3", not has(f, "P3"))
        root = _tree(base, "loop")
        os.symlink("loop", root / "devlog/loop")
        (root / "devlog/README.md").write_text(_README + "- [x](loop/dl20260103_x.md)\n")
        try:
            run_live(root, 40, 0, 25)
            check("a malformed index destination does not end the check", True)
        except Exception as exc:
            check(f"a malformed index destination does not end the check ({type(exc).__name__})", False)

        print("v4 - headers:")
        headless = "# A topic phrase\n\n## Summary\n"
        root = _tree(base, "all_headless", readme=_README)
        (root / "devlog/dl20260101_clean.md").write_text(headless)
        try:
            f = run_live(root, 40, 0, 25)[0]
            check("an all-headless window reports P1c instead of failing the guard",
                  has(f, "P1c", "dl20260101_clean.md"))
            check("an all-headless window does not blame the template (no P1a)", not has(f, "P1a"))
        except ExtractionError as exc:
            check(f"an all-headless window reports P1c instead of failing the guard ({exc})", False)

        print("v5 - fenced blocks in a devlog header:")
        def missing(label, text):
            p = Path(base) / f"fence_{label}.md"
            p.write_text(text)
            return bool(devlogs_without_date([p]))
        p = Path(base) / "fence_stray_ends.md"
        p.write_text("# T\n\n```markdown\n**Purpose:** before\n## ex\n**Date:** after\n")
        hdr = "\n".join(_header_lines(p))
        check("a stray opener is text: what follows it is header until a ## ends it",
              "**Purpose:** before" in hdr and "**Date:** after" not in hdr)
        check("an unclosed fence is a stray line: a Date after it still counts",
              not missing("unclosed_then_date", "# T\n\n```\n**Date:** 2026-01-01\n## Summary\n"))
        check("a later unclosed fence does not re-admit an earlier closed fence's Date",
              missing("closed_then_unclosed", "# T\n\n```\n**Date:** 1999-01-01\n```\n\n**Purpose:** x\n\n"
                                              "```\nnot closed\n## Summary\n"))
        check("a closer followed by a non-breaking space does not close the fence",
              missing("nbsp_close", "# T\n\n```\nfoo\n```\u00a0\n**Date:** 1999-01-01\n```\n\n"
                                    "**Purpose:** x\n\n## Summary\n"))
        check("a fenced **Date:** does not stand in for a missing one",
              missing("fenced_date", "# T\n\n```\n**Date:** 2026-01-01\n```\n\n**Purpose:** x\n\n## Summary\n"))
        check("a shorter fence inside a longer one is text (nesting)",
              not missing("nested", "# T\n\n````markdown\n```markdown\n## ex\n```\n````\n\n"
                                    "**Date:** 2026-01-01\n## Summary\n"))
        check("a shorter closer inside a longer fence does not end it",
              missing("short_closer", "# T\n\n````\n```\n**Date:** 1999-01-01\n````\n\n**Purpose:** x\n\n"
                                      "## Summary\n"))
        check("a tilde line does not close a backtick fence",
              missing("tilde_in_backtick", "# T\n\n```\n~~~\n**Date:** 1999-01-01\n```\n\n**Purpose:** x\n\n"
                                           "## Summary\n"))
        check("a closer followed by a form feed does not close the fence",
              missing("formfeed_close", "# T\n\n```\nfoo\n```\f\n**Date:** 1999-01-01\n```\n\n"
                                        "**Purpose:** x\n\n## Summary\n"))
        # The one-pass resolver against its definition: re-scan, demoting each unclosed opener.
        def reference(lines):
            as_text = set()
            while True:
                flags, fence, opened_at = [], None, None
                for i, line in enumerate(lines):
                    if fence is not None:
                        flags.append(True)
                        if _fence_closes(line, fence):
                            fence = None
                        continue
                    if i not in as_text and (op := _fence_open(line)) is not None:
                        fence, opened_at = op, i
                        flags.append(True)
                        continue
                    flags.append(False)
                if fence is None:
                    return flags
                as_text.add(opened_at)
        import random
        rng = random.Random(20260923)
        vocab = ["```", "````", "~~~", "~~~~", "```py", "``` a`b", "``` ", "```\u00a0", "  ```",
                 "text", "## h", "**Date:** x", ""]
        agree = all(_fenced_flags(d) == reference(d)
                    for d in ([rng.choice(vocab) for _ in range(rng.randint(0, 14))] for _ in range(3000)))
        check("the one-pass resolver agrees with the re-scan definition on 3000 random documents", agree)
        check("a tilde fence hides its ## heading",
              not missing("tilde", "# T\n\n~~~\n## ex\n~~~\n\n**Date:** 2026-01-01\n## Summary\n"))
        check("a longer closer closes the fence",
              not missing("long_close", "# T\n\n```\n## ex\n````\n\n**Date:** 2026-01-01\n## Summary\n"))
        check("a closer with trailing text does not close the fence",
              missing("text_close", "# T\n\n```\n## ex\n``` not a closer\n**Date:** 2026-01-01\n"
                                    "## Summary\n"))
        check("a backtick in an info string means the line opens no fence",
              not missing("info_backtick", "# T\n\n``` a`b\n**Date:** 2026-01-01\n```\n## Summary\n"))
        f = find(_tree(base, "fenced_model", devlog_extra="\n```\nCodex `gpt-9.9-rogue`/high\n```"))
        check("a model id inside a fenced example is not a header use", not has(f, "P2a"))

        print("v5 - the depth table is read strictly:")
        def table_err(label, proto_text):
            root = _tree(base, label)
            (root / "doc/ai_collaboration_protocol.md").write_text(proto_text)
            try:
                run_live(root, 40, 0, 25)
                return None
            except ExtractionError as exc:
                return str(exc)
        clean_proto = _PROTO.format(exceptions="", row=_ROW_GROK)
        err = table_err("dup_grok", clean_proto.replace("| Subject | Codex | Grok |",
                                                        "| Subject | Codex | Grok | Grok |")
                        .replace("|---|---|---|", "|---|---|---|---|")
                        .replace(_ROW_GROK, _ROW_GROK + " `grok-4.6`, `high` |"))
        check("a duplicate Grok header fails closed", err is not None and "exactly once" in err)
        err = table_err("no_codex", clean_proto.replace("| Subject | Codex | Grok |", "| Subject | Other | Grok |"))
        check("a missing Codex header fails closed", err is not None and "Codex" in err)
        err = table_err("uneven", clean_proto.replace(_ROW_GROK, "| Plan audit | `gpt-5.6-terra`, `high` |"))
        check("a row narrower than the header fails closed", err is not None and "cells" in err)
        root = _tree(base, "double_pipe")
        (root / "doc/ai_collaboration_protocol.md").write_text(
            clean_proto.replace(_ROW_GROK, "|| `gpt-5.6-terra`, `high` | `grok-4.7`, `high` |"))
        try:
            cols = table_models(root / "doc/ai_collaboration_protocol.md")
            check("an empty first cell keeps the columns aligned",
                  cols == {"codex": {"gpt-5.6-terra"}, "grok": {"grok-4.7"}})
        except ExtractionError as exc:
            check(f"an empty first cell keeps the columns aligned ({exc})", False)
        spaced = clean_proto.replace(_ROW_GROK + "\n", _ROW_GROK + "\n\nprose between\n\n"
                                     "| Later | `gpt-9.9-y`, `high` | `grok-4.7`, `high` |\n")
        root = _tree(base, "second_table")
        (root / "doc/ai_collaboration_protocol.md").write_text(spaced)
        f = run_live(root, 40, 0, 25)[0]
        check("a second | block after the table is not part of it", not has(f, "P2b", "gpt-9.9-y"))

        try:
            f = find(_tree(base, "fenced_before", exceptions="```\n| a | b | c |\n```\n\n"))
            check("a | line in a fenced block before the table is not its header", not f)
        except ExtractionError as exc:
            check(f"a | line in a fenced block before the table is not its header ({exc})", False)
        example = ("```markdown\n| Subject | Codex | Grok |\n|---|---|---|\n"
                   "| X | `gpt-9.9-example`, `high` | `grok-4.7`, `high` |\n"
                   "## a heading inside the example\n```\n\n")
        try:
            f = run_live(_tree(base, "example_table_in_fence", exceptions=example), 40, 0, 25)[0]
            check("a fenced example straddling a heading is skipped and the real table read", not f)
        except ExtractionError as exc:
            check(f"a fenced example straddling a heading is skipped and the real table read ({exc})", False)
        quoted = _PROTO.format(exceptions="", row=_ROW_GROK).replace(
            "### 9.1 Depth Selection",
            "```\n### 9.1 Depth Selection\n| Subject | Codex | Grok |\n|---|---|---|\n"
            "| X | `gpt-9.9-example`, `high` | `grok-4.7`, `high` |\n```\n\n### 9.1 Depth Selection", 1)
        root = _tree(base, "quoted_heading")
        (root / "doc/ai_collaboration_protocol.md").write_text(quoted)
        try:
            f = run_live(root, 40, 0, 25)[0]
            check("a fenced copy of the section heading is not the section", not f)
        except ExtractionError as exc:
            check(f"a fenced copy of the section heading is not the section ({exc})", False)
        forged = ("```\n### 9.1 Depth Selection\n| Subject | Codex | Grok |\n|---|---|---|\n"
                  "| X | `gpt-9.9-example`, `high` | `grok-4.7`, `high` |\n\n")
        root = _tree(base, "forged_anchor")
        (root / "doc/ai_collaboration_protocol.md").write_text(_PROTO.format(exceptions="", row=_ROW_GROK)
            .replace("### 9.1 Depth Selection", forged + "### 9.1 Depth Selection", 1))
        try:
            run_live(root, 40, 0, 25)
            check("an unclosed fence before a forged section heading fails closed", False)
        except ExtractionError as exc:
            text = (root / "doc/ai_collaboration_protocol.md").read_text().split("\n")
            hits = [n for n, l in enumerate(text, 1) if l == "### 9.1 Depth Selection"]
            check("an unclosed fence before a forged section heading fails closed, naming the lines",
                  "appears 2 times" in str(exc) and f"at lines {hits[0]}, {hits[1]}" in str(exc))
        root = _tree(base, "later_subsection")
        (root / "doc/ai_collaboration_protocol.md").write_text(_PROTO.format(exceptions="", row=_ROW_GROK)
            .replace("### Workflow", "### 9.10 A later subsection\n\n### Workflow"))
        try:
            f = run_live(root, 40, 0, 25)[0]
            check("a ### 9.10 heading is not a second depth section", not f)
        except ExtractionError as exc:
            check(f"a ### 9.10 heading is not a second depth section ({exc})", False)
        for label, heading in (("tab_anchor", "### 9.1\tDepth Selection"), ("bare_anchor", "### 9.1")):
            root = _tree(base, label)
            (root / "doc/ai_collaboration_protocol.md").write_text(
                _PROTO.format(exceptions="", row=_ROW_GROK).replace("### 9.1 Depth Selection", heading, 1))
            try:
                f = run_live(root, 40, 0, 25)[0]
                check(f"{heading!r} is the depth heading", not f)
            except ExtractionError as exc:
                check(f"{heading!r} is the depth heading ({exc})", False)
        root = _tree(base, "formfeed_anchor")
        (root / "doc/ai_collaboration_protocol.md").write_text(_PROTO.format(exceptions="", row=_ROW_GROK)
            .replace("### Workflow", "### 9.1\f0 not a heading\n\n### Workflow"))
        try:
            f = run_live(root, 40, 0, 25)[0]
            check("### 9.1 followed by a form feed is not the depth heading", not f)
        except ExtractionError as exc:
            check(f"### 9.1 followed by a form feed is not the depth heading ({exc})", False)
        root = _tree(base, "anchor_in_unclosed_fence")
        (root / "doc/ai_collaboration_protocol.md").write_text(_PROTO.format(exceptions="", row=_ROW_GROK)
            .replace("### 9.1 Depth Selection", "```\n\n### 9.1 Depth Selection", 1))
        try:
            run_live(root, 40, 0, 25)
            check("a lone depth heading inside an unclosed fence fails closed", False)
        except ExtractionError as exc:
            text = (root / "doc/ai_collaboration_protocol.md").read_text().split("\n")
            heading = text.index("### 9.1 Depth Selection") + 1
            opener = max(n for n, l in enumerate(text[:heading - 1], 1) if l == "```")
            check("a lone depth heading inside an unclosed fence fails closed, naming both lines",
                  f"depth anchor at line {heading} lies inside a fence opened at line {opener} "
                  "and never closed" in str(exc))
        unclosed_example = ("```markdown\n| Subject | Codex | Grok |\n|---|---|---|\n"
                            "| X | `gpt-9.9-example`, `high` | `grok-4.7`, `high` |\n\n")
        try:
            f = run_live(_tree(base, "unclosed_example", exceptions=unclosed_example), 40, 0, 25)[0]
            check("an unclosed fenced example table is never adopted as the depth table",
                  not any(s == "gpt-9.9-example" for _, _, s, _ in f))
        except ExtractionError as exc:
            check("an unclosed fenced example table is never adopted as the depth table "
                  "(fails closed on an empty table)", "depth table rows" in str(exc))

        print("exit policy:")
        drifty = _tree(base, "exit_drift", extra_devlogs=[("dl20260102_unindexed.md", plain)])
        clean = _tree(base, "exit_clean")
        check("--strict exits 1 on DRIFT", run(drifty, commits=0, strict=True, quiet=True) == 1)
        check("--strict exits 0 on a clean tree", run(clean, commits=0, strict=True, quiet=True) == 0)
        check("without --strict DRIFT exits 0", run(drifty, commits=0, quiet=True) == 0)
        info_only = _tree(base, "exit_info", title="# Devlog 2026-01-01 — a topic")
        check("--strict exits 0 on INFO alone",
              run(info_only, commits=0, strict=True, quiet=True) == 0)

    print()
    if failures:
        print(f"self-test FAILED: {len(failures)} case(s)")
        return 1
    print("self-test passed")
    return 0


def _at_least(low):
    def parse(text):
        n = int(text)
        if n < low:
            raise argparse.ArgumentTypeError(f"must be >= {low}")
        return n
    return parse


def _percent(text):
    x = float(text)
    if not 0.0 <= x <= 100.0:     # also rejects nan and inf
        raise argparse.ArgumentTypeError("must be a percentage in [0, 100]")
    return x


def main():
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("--recent", type=_at_least(1), default=40, help="devlogs in the window (default 40)")
    ap.add_argument("--commits", type=_at_least(0), default=40, help="commit bodies to read (default 40)")
    ap.add_argument("--threshold", type=_percent, default=25.0,
                    help="percent of the window a practice field needs before P1b reports it")
    ap.add_argument("--strict", action="store_true", help="exit 1 when a DRIFT finding is reported")
    ap.add_argument("--self-test", action="store_true", help="run the fixtures; this is the gate")
    args = ap.parse_args()

    if args.self_test:
        return self_test()
    return run(Path(__file__).resolve().parent.parent, args.recent, args.commits,
               args.threshold, args.strict)


if __name__ == "__main__":
    sys.exit(main())
