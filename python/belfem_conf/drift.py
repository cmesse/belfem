"""The two drift directions.

schema -> code   every anchor still resolves. Catches renames and removals.
code -> schema   every parse-call literal is known to the schema. Catches the
                 dangerous case: a key added to a factory and to no document.

The second is the one no amount of human discipline covers, because whoever
adds the key is precisely the person unaware the document exists. It is also
the check that would have caught both holes found during development.

Not implemented: doc -> code. A claim can be false while every key and anchor
is correct — the reference once described a shipped feature as unimplemented,
and no grep finds that. It needs a reader.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from . import enums
from .repo import Sources
from .schema import Anchor, Schema, known_keys

# The accessors that read a deck key. `section`/`section_exists` are included
# because a section name is part of the contract too. Checked against
# cl_Input_Section.hpp 2026-08-13 (an earlier draft listed a `get_uint` that
# exists nowhere in src/ — this list is transcribed, keep it reviewed).
PARSE_CALLS = (
    "key_exists", "get_real", "get_bool", "get_int", "get_string",
    "get_value", "get_ids", "get_id_groups", "get_reals", "section_exists",
    "section",
)

_LITERAL = re.compile(
    r"\b(" + "|".join(PARSE_CALLS) + r")\s*\(\s*\"([^\"]+)\""
)


@dataclass
class Finding:
    kind: str          # "unresolved-anchor" | "wrong-consumer" | "unknown-key"
    detail: str
    where: str


@dataclass
class Report:
    findings: list[Finding]
    anchors_declared: int
    anchors_checked: int
    anchors_resolved: int
    keys_seen: int
    keys_known: int
    files_scanned: int
    enums_checked: int = 0

    @property
    def ok(self) -> bool:
        return not any(f.kind != "wrong-consumer" for f in self.findings)


def check_anchors(schema: Schema, sources: Sources) -> list[Finding]:
    """schema -> code."""
    findings: list[Finding] = []
    for anchor in schema.anchors:
        token = anchor.token.strip()
        # a bare key literal is stored with its quotes; search for them too,
        # so `"unit"` cannot match the word `unit` inside a comment
        hits = sources.find(token, anchor.consumer)
        if hits:
            continue

        # widen: the token may be real but attributed to the wrong consumer,
        # which is a much weaker finding than not existing at all
        anywhere = sources.find(token)
        if anywhere:
            findings.append(
                Finding(
                    "wrong-consumer",
                    f"{token!r} not in {anchor.consumer}, found in "
                    f"{anywhere[0][0]}:{anywhere[0][1]}",
                    anchor.path,
                )
            )
        else:
            findings.append(
                Finding("unresolved-anchor", f"{token!r} resolves nowhere", anchor.path)
            )
    return findings


def check_keys(schema: Schema, sources: Sources) -> tuple[list[Finding], int]:
    """code -> schema. Returns (findings, number of distinct literals seen).

    The scan is over the file as ONE string, not line by line. It matters:
    `mesh.unit` is parsed only as `get_string(\\n "unit" )` split across two
    lines in cl_MaxwellFactory.cpp, and a per-line regex never sees it — a
    silent false "clean" in exactly the direction this check exists for.
    """
    known = known_keys(schema)
    seen: dict[str, tuple[str, int]] = {}

    for path in sources.consumer_paths():
        numbered = sources.lines(path)          # comment lines already gone
        text = "\n".join(line for _n, line in numbered)
        linemap = [n for n, _line in numbered]
        for match in _LITERAL.finditer(text):
            n = linemap[text.count("\n", 0, match.start())]
            seen.setdefault(match.group(2).lower(), (path.name, n))

    findings = [
        Finding("unknown-key", f"{key!r} is parsed but absent from the schema",
                f"{where[0]}:{where[1]}")
        for key, where in sorted(seen.items())
        if key not in known
    ]
    return findings, len(seen)


def run(schema: Schema, sources: Sources) -> Report:
    anchor_findings = check_anchors(schema, sources)
    key_findings, keys_seen = check_keys(schema, sources)

    enum_refs = enums.collect_refs(schema.data)
    enum_findings = [
        Finding("enum-mismatch", problem, where)
        for where, problem in enums.check(schema.data, sources.root)
    ]

    resolved = len(schema.anchors) - sum(
        1 for f in anchor_findings if f.kind == "unresolved-anchor"
    )

    return Report(
        enums_checked=len(enum_refs),
        findings=anchor_findings + key_findings + enum_findings,
        # declared vs checked is the coverage self-report: if these disagree the
        # walker is skipping anchors, which is exactly how the first prototype
        # managed to pass a check it had never run
        anchors_declared=schema.anchor_fields_seen,
        anchors_checked=len(schema.anchors),
        anchors_resolved=resolved,
        keys_seen=keys_seen,
        keys_known=len(known_keys(schema)),
        files_scanned=len(sources.consumer_paths()),
    )
