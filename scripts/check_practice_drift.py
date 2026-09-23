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
"""

import argparse
import re
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path, PurePosixPath

# A version is dot-separated integers; a suffix segment starts with a letter (gpt-6-astra,
# gpt-5.4-mini, grok-4.7-build-fast).  Seat names such as "Claude Opus 5" fall outside the
# alphabet, so they cannot leak into a model set, and a trailing sentence period is not absorbed.
MODEL = re.compile(r"\b((?:gpt|grok)-\d+(?:\.\d+)*(?:-[A-Za-z][A-Za-z0-9]*)*)")
GROK_VERSION = re.compile(r"grok-(\d+(?:\.\d+)*)")
FIELD = re.compile(r"^\*\*([A-Za-z][A-Za-z /]*):\*\*")
# The capture is the link DESTINATION inside (...), not the label inside [...]; the two differ
# whenever an entry is given a short name.
INDEX_ENTRY = re.compile(r"^(?:- |#{2,4} )\[[^\]]*\]\(([^)\s]+)\)")
CASE_ARM = re.compile(r"^\s*([A-Za-z0-9.|-]+)\)\s*;;")
OLD_TITLE = re.compile(r"^# Devlog \d{4}-\d{2}-\d{2}")
DEVLOG_DATE = re.compile(r"dl(\d{8})")

TEMPLATE_ANCHOR = "### Required Header"
TEMPLATE_SUCCESSOR = "### When to Write a Devlog"
DEPTH_ANCHOR = "### 9.1"
EXCEPTION_ANCHOR = "**Allowlisted but named by no row (by design):**"

WRAPPERS = {"codex": (".claude/scripts/ask_codex.sh", "CODEX_MODEL"),
            "grok": (".claude/scripts/ask_grok.sh", "GROK_MODEL")}

DRIFT, INFO = "DRIFT", "INFO"


class ExtractionError(RuntimeError):
    """A required extract came back empty or an anchor was not found."""


# --------------------------------------------------------------------------- extractors
# Each takes its path(s) explicitly so the self-test can point it at a fixture.

def _lines(path):
    return Path(path).read_text(encoding="utf-8", errors="replace").splitlines()


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
    filename breaks the tie deterministically.
    """
    def key(p):
        m = DEVLOG_DATE.search(p.name)
        return (m.group(1) if m else "0", p.name)
    return sorted(Path(devlog_dir).glob("dl*.md"), key=key, reverse=True)[:recent]


def _header_lines(path):
    out = []
    for line in _lines(path):
        if line.startswith("## "):
            break
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


def _depth_span(proto_path):
    lines = _lines(proto_path)
    try:
        start = next(i for i, l in enumerate(lines) if l.startswith(DEPTH_ANCHOR))
    except StopIteration:
        raise ExtractionError(f"anchor not found: {DEPTH_ANCHOR}")
    end = next((i for i, l in enumerate(lines[start + 1:], start + 1)
                if l.startswith("### ") or l.startswith("## ")), len(lines))
    return lines[start:end]


def table_models(proto_path):
    """Models named by the depth table's rows.  Rows declare the tiers."""
    return {m for l in _depth_span(proto_path) if l.startswith("|") for m in MODEL.findall(l)}


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

    An index line may mention a second filename in its prose, and an entry's label may differ
    from its destination; matching either would count a file as indexed when it is not.
    """
    out = set()
    for l in _lines(readme_path):
        if m := INDEX_ENTRY.match(l):
            out.add(PurePosixPath(m.group(1).split("#", 1)[0]).name)
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


def probe_models(rows, exceptions, allowlists, messages, defaults, headers, commits):
    """P2a-P2e.  `allowlists`, `messages` and `defaults` are keyed by seat."""
    findings = []
    reported = set()

    # The Grok seat is always the latest allowlisted grok.  That model is declared by the policy
    # itself, whether or not a row happens to name it.
    policy_grok = _newest_grok(allowlists.get("grok", set()))
    declared = rows | exceptions | ({policy_grok} if policy_grok else set())

    # P2c - a default that is not what the policy names.  "Equals the latest", not "appears in
    # some row": a default left on an older row would otherwise pass while the policy moved on.
    if (d := defaults.get("grok")) and policy_grok and d != policy_grok:
        findings.append((DRIFT, "P2c", d,
                         f"grok default is {d}; the policy names the latest allowlisted grok, "
                         f"{policy_grok}"))
        reported.add(d)
    if (d := defaults.get("codex")) and d not in rows:
        findings.append((DRIFT, "P2c", d, f"codex default {d} is named by no depth row"))
        reported.add(d)

    # P2a - used in recent records, declared nowhere.  Headers trigger; commit bodies only
    # corroborate.  A model P2c already diagnoses is not reported twice.
    for model in sorted(headers):
        if model in declared or model in reported:
            continue
        also = f", and in {commits[model]} commit bodies" if commits.get(model) else ""
        findings.append((DRIFT, "P2a", model,
                         f"{model} named in {headers[model]} recent devlog headers{also}; "
                         "no depth row names it and it is not declared by design"))
        reported.add(model)

    # P2b - a row naming a model its own seat's wrapper would reject.  Bound to the seat: every
    # gpt row is absent from the grok allowlist and vice versa.
    for model in sorted(rows):
        seat = _seat(model)
        if seat in allowlists and model not in allowlists[seat]:
            findings.append((DRIFT, "P2b", model,
                             f"depth row names {model}, which ask_{seat}.sh would reject"))

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


def probe_template(declared, used, denom, threshold):
    """P1a (DRIFT) and P1b (INFO).

    The directions are not symmetric.  A declared field nobody writes means the template
    describes something dead: a fact about the document, so it fires at zero uses with no
    threshold.  A field practice uses but the template omits is not a defect — the template states
    that its list is not a schema — so it is reported only once it is common (`threshold`, in
    percent of the window), and only as INFO.
    """
    findings = [(DRIFT, "P1a", f,
                 f"template declares **{f}:** — used by none of the last {denom} entries")
                for f in sorted(declared) if not used.get(f)]
    findings += [(INFO, "P1b", f,
                  f"**{f}:** used by {used[f]}/{denom} recent entries, not declared by the template")
                 for f in sorted(used) if f not in declared and 100 * used[f] / denom >= threshold]
    return findings


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
    rows = table_models(proto)
    exceptions = exception_models(proto)
    allowlists = {s: wrapper_models(root / p) for s, (p, _) in WRAPPERS.items()}
    messages = {s: wrapper_message_models(root / p) for s, (p, _) in WRAPPERS.items()}
    defaults = {s: wrapper_default(root / p, v) for s, (p, v) in WRAPPERS.items()}
    headers, commit_hits, commit_status = record_models(window, root, commits_n, commit_texts)
    files, targets = devlog_files(devlog_dir), index_targets(readme)
    old, total = title_forms(window)

    guard_required({
        "template fields": declared,
        "devlog headers": used,
        "depth table rows": rows,
        "codex allowlist": allowlists["codex"],
        "grok allowlist": allowlists["grok"],
        "codex default": defaults["codex"],
        "grok default": defaults["grok"],
        "devlog files": files,
        "index targets": targets,
    })

    findings = probe_models(rows, exceptions, allowlists, messages, defaults, headers, commit_hits)
    findings += probe_index(files, targets)
    findings += probe_template(declared, used, denom, threshold)
    findings += probe_titles(old, total)

    window_note = (f"at {head_commit(root)}: {denom} of the newest {len(window)} devlogs "
                   f"(ordered by filename date, then name) carry a header; {len(files)} files, "
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
    except ExtractionError as exc:
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
        check("an empty by-design list is accepted by the guard", True)

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


def main():
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("--recent", type=int, default=40, help="devlogs in the window (default 40)")
    ap.add_argument("--commits", type=int, default=40, help="commit bodies to read (default 40)")
    ap.add_argument("--threshold", type=float, default=25.0,
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
