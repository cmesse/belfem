"""Command line entry point."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import drift, schema as schema_mod
from .document import Document, read_verbatim
from .parser import ParseError
from .repo import Sources, find_root

GREEN, RED, YELLOW, DIM, RESET = (
    "\033[32m", "\033[31m", "\033[33m", "\033[2m", "\033[0m"
)


def _paint(enabled: bool):
    if enabled:
        return GREEN, RED, YELLOW, DIM, RESET
    return "", "", "", "", ""


def cmd_drift(args) -> int:
    root = find_root(Path(args.root) if args.root else None)
    sources = Sources(root)
    sch = schema_mod.load(root)
    report = drift.run(sch, sources)

    green, red, yellow, dim, reset = _paint(sys.stdout.isatty() and not args.no_color)

    # Coverage first, deliberately. A checker that does not say how much it
    # examined is indistinguishable from one that passes everything.
    print(f"{dim}schema   {root / 'doc/input_schema.yaml'}{reset}")
    print(
        f"anchors  {report.anchors_resolved}/{report.anchors_checked} resolve"
        f"   (declared in file: {report.anchors_declared})"
    )
    if report.anchors_checked < report.anchors_declared:
        print(
            f"{red}         WARNING: {report.anchors_declared - report.anchors_checked}"
            f" declared anchors were not checked — the walker is skipping some{reset}"
        )
    print(
        f"keys     {report.keys_seen} parsed literals across "
        f"{report.files_scanned} consumers; schema knows {report.keys_known}"
    )
    print(
        f"enums    {report.enums_checked} value list(s) re-derived from "
        f"to_string() in src/"
    )
    print()

    unknown = [f for f in report.findings if f.kind == "unknown-key"]
    unresolved = [f for f in report.findings if f.kind == "unresolved-anchor"]
    misfiled = [f for f in report.findings if f.kind == "wrong-consumer"]
    enum_bad = [f for f in report.findings if f.kind == "enum-mismatch"]

    if enum_bad:
        print(f"{red}ENUM VALUE LISTS OUT OF STEP WITH THE CODE ({len(enum_bad)}){reset}")
        print("  The C++ resolves these by comparing against to_string(), so the")
        print("  switch is the authority and the schema's list is a copy.")
        for f in enum_bad:
            print(f"    {f.detail}   {dim}{f.where}{reset}")
        print()

    if unknown:
        print(f"{red}KEYS PARSED BUT NOT IN THE SCHEMA ({len(unknown)}){reset}")
        print("  A key the code reads and no document describes. This is the")
        print("  failure mode human review does not catch.")
        for f in unknown:
            print(f"    {f.detail}   {dim}{f.where}{reset}")
        print()

    if unresolved:
        print(f"{red}ANCHORS THAT RESOLVE NOWHERE ({len(unresolved)}){reset}")
        for f in unresolved:
            print(f"    {f.detail}   {dim}{f.where}{reset}")
        print()

    if misfiled and args.verbose:
        print(f"{yellow}ANCHORS FILED UNDER THE WRONG CONSUMER ({len(misfiled)}){reset}")
        print(f"  {dim}Resolve, but not where the schema says. Not a failure.{reset}")
        for f in misfiled:
            print(f"    {f.detail}   {dim}{f.where}{reset}")
        print()
    elif misfiled:
        print(f"{dim}({len(misfiled)} anchors resolve under a different consumer; "
              f"-v to list){reset}")
        print()

    if report.ok:
        print(f"{green}drift: clean{reset}")
        return 0

    print(f"{red}drift: {len(unknown) + len(unresolved) + len(enum_bad)} finding(s){reset}")
    print(
        f"{dim}If you added an input key, add it to doc/input_schema.yaml and\n"
        f"doc/input_file_reference.md — see CLAUDE.md, \"The Input Contract:\n"
        f"Two Artifacts, One Rule\".{reset}"
    )
    return 1


def cmd_roundtrip(args) -> int:
    """Parse and re-serialise decks, requiring byte-identical output.

    This is the gate the rest of the tooling stands on. If a deck cannot be
    reproduced byte for byte, the parser is losing something, and no validator
    or editor built on it can be trusted not to mangle a file on save.
    """
    green, red, yellow, dim, reset = _paint(
        sys.stdout.isatty() and not args.no_color)

    if args.decks:
        decks = [Path(d) for d in args.decks]
    else:
        root = find_root(Path(args.root) if args.root else None)
        decks = sorted(root.glob("examples/*/input.conf"))
        if not decks:
            print(f"{red}no decks found under examples/{reset}")
            return 1

    failures = 0
    for deck in decks:
        label = str(deck)
        try:
            original = read_verbatim(deck)
            doc = Document.loads(original, path=deck)
            produced = doc.dumps()
        except ParseError as exc:
            print(f"  {red}PARSE {reset}{label}\n         {exc}")
            failures += 1
            continue
        except OSError as exc:
            print(f"  {red}READ  {reset}{label}\n         {exc}")
            failures += 1
            continue

        if produced != original:
            failures += 1
            where = next(
                (i for i, (a, b) in enumerate(zip(original, produced)) if a != b),
                min(len(original), len(produced)),
            )
            line = original.count("\n", 0, where) + 1
            print(f"  {red}DIFFER{reset} {label}   first difference at line {line}")
            continue

        probed, problem = (0, None) if args.no_probe else _probe_spans(original)
        if problem:
            failures += 1
            print(f"  {red}SPAN  {reset} {label}\n         {problem}")
            continue

        sections = sum(1 for _ in doc.walk())
        stmts = _count_statements(doc.root)
        note = f", {probed} value spans probed" if probed else ""
        print(f"  {green}ok    {reset} {label}   "
              f"{dim}{len(original)} bytes, {sections} sections, "
              f"{stmts} statements{note}{reset}")

    print()
    if failures:
        print(f"{red}round-trip: {failures} of {len(decks)} deck(s) not "
              f"reproduced byte-for-byte{reset}")
        return 1
    print(f"{green}round-trip: {len(decks)} deck(s) reproduced byte-for-byte{reset}")
    return 0


def cmd_check(args) -> int:
    from . import validate as validate_mod
    from .validate import ERROR, WARNING

    root = find_root(Path(args.root) if args.root else None)
    schema = schema_mod.load(root)
    green, red, yellow, dim, reset = _paint(
        sys.stdout.isatty() and not args.no_color)

    decks = ([Path(d) for d in args.decks] if args.decks
             else sorted(root.glob("examples/*/input.conf")))
    if not decks:
        print(f"{red}no decks given and none found under examples/{reset}")
        return 1

    worst = 0
    for deck in decks:
        try:
            doc = Document.load(deck)
        except ParseError as exc:
            print(f"{red}{deck}: {exc}{reset}")
            worst = 1
            continue

        result = validate_mod.validate(doc, schema.data, root)
        errors = result.errors
        warns = sum(1 for d in result.diagnostics if d.severity == WARNING)
        # a note is accepted-and-ignored input, not a warning — a deck whose
        # only diagnostics are notes is ok, and painting notes yellow made
        # sidecoating's inert `/` look like a problem forever
        mark = (f"{red}FAIL{reset}" if errors
                else f"{yellow}warn{reset}" if warns else f"{green} ok {reset}")
        print(f"{mark}  {deck}")

        for d in result.diagnostics:
            colour = (red if d.severity == ERROR
                      else yellow if d.severity == WARNING else dim)
            print(f"      {colour}{d.severity}{reset} {deck.name}:{d.line}  {d.message}")
            print(f"            {dim}{d.where}{reset}")
            if d.hint:
                print(f"            {dim}{d.hint}{reset}")

        if args.verbose or not args.decks:
            done = ", ".join(f"{v} {k}" for k, v in sorted(result.checked.items()))
            print(f"      {dim}checked: {done}{reset}")
        if args.verbose:
            for s in result.skipped:
                print(f"      {dim}not checked: {s}{reset}")
        print()
        worst = max(worst, 1 if errors else 0)

    if not args.verbose:
        print(f"{dim}-v lists what was not checked, and why{reset}")
    return worst


def _count_statements(section) -> int:
    return len(section.statements) + sum(
        _count_statements(s) for s in section.sections)


def _all_statements(section) -> list:
    out = list(section.statements)
    for child in section.sections:
        out.extend(_all_statements(child))
    return out


PROBE = "@@belfem-conf-probe@@"


def _probe_spans(original: str) -> tuple[int, str | None]:
    """Rewrite every value in turn and confirm only that span moved.

    Without this the round-trip check is close to tautological: `dumps()`
    returns the original text when nothing was edited, so it would pass even if
    the parser produced no tree at all. This exercises the spans themselves —
    for each statement, replacing its value must change exactly those bytes and
    nothing else, and re-parsing the result must read the new value back.
    """
    reference = Document.loads(original)
    targets = [s for s in _all_statements(reference.root) if s.value_span]

    for stmt in targets:
        doc = Document.loads(original)
        here = next(
            (t for t in _all_statements(doc.root) if t.span == stmt.span), None)
        if here is None:            # invariant break, not a deck problem
            return len(targets), (
                f"re-parse did not reproduce the span of {stmt.key!r}")
        doc.replace_value(here, PROBE)
        produced = doc.dumps()

        start, end = stmt.value_span
        if produced != original[:start] + PROBE + original[end:]:
            return len(targets), (
                f"editing {stmt.key!r} changed bytes outside its value span")

        # the edit is after the key, so this statement's start does not move
        again = Document.loads(produced)
        back = next(
            (t for t in _all_statements(again.root) if t.span[0] == stmt.span[0]),
            None)
        if back is None or back.value != PROBE:
            return len(targets), (
                f"re-parsing after editing {stmt.key!r} did not read it back")

    return len(targets), None


def main(argv: list[str] | None = None) -> int:
    # common options go on a parent parser so they are accepted on either side
    # of the subcommand: `belfem-conf --root X drift` and `belfem-conf drift
    # --root X` both work, which is what anyone actually types
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--root", help="repository root (default: auto-detect)")
    common.add_argument("--no-color", action="store_true", help="plain output")

    parser = argparse.ArgumentParser(
        prog="belfem-conf",
        parents=[common],
        description="Tools for the BELFEM input.conf contract.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    drift_p = sub.add_parser(
        "drift",
        parents=[common],
        help="check doc/input_schema.yaml against the C++ parse sites",
        description=(
            "Compares the schema with the code in both directions: every anchor "
            "must still resolve, and every parsed key must be described. "
            "Exits non-zero on either."
        ),
    )
    drift_p.add_argument(
        "-v", "--verbose", action="store_true",
        help="also list anchors filed under the wrong consumer",
    )
    drift_p.set_defaults(func=cmd_drift)

    rt_p = sub.add_parser(
        "roundtrip",
        parents=[common],
        help="parse and re-serialise decks, requiring byte-identical output",
        description=(
            "The gate the rest of the tooling stands on. With no arguments, "
            "runs over every examples/*/input.conf in the checkout."
        ),
    )
    rt_p.add_argument("decks", nargs="*", help="decks to check (default: all examples)")
    rt_p.add_argument(
        "--no-probe", action="store_true",
        help="skip the per-value span probe (identity check only, much weaker)",
    )
    rt_p.set_defaults(func=cmd_roundtrip)

    chk_p = sub.add_parser(
        "check",
        parents=[common],
        help="validate deck(s) against doc/input_schema.yaml",
        description=(
            "Checks structure, unknown keys, enum values, unit dimensions and "
            "cross-references. Reports what it did NOT check, and why."
        ),
    )
    chk_p.add_argument("decks", nargs="*", help="decks (default: all examples)")
    chk_p.add_argument("-v", "--verbose", action="store_true",
                       help="list coverage and the checks that were skipped")
    chk_p.set_defaults(func=cmd_check)

    args = parser.parse_args(argv)
    return args.func(args)
