#!/usr/bin/env python3
"""Resolve @@file|regex@@ placeholders in markdown to file:line citations.

Usage: resolve_doc_cites.py [--check] doc.md [doc2.md ...]

A placeholder @@cl_Foo.cpp|Foo::bar\(@@ becomes `cl_Foo.cpp:NNN`, where NNN is the
first line of cl_Foo.cpp (searched under src/) matching the regex. With --check,
existing `file:NNN` citations are listed with the text found at that line instead,
so a reader can judge drift. Line numbers drift with every edit; the placeholder
form is the one to keep in the source of a document that is regenerated.
"""
import re, sys, pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
SRC = ROOT / "src"

def find(fname):
    hits = sorted(SRC.rglob(fname))
    return hits[0] if hits else None

def resolve(text, doc):
    def sub(m):
        fname, pat = m.group(1), m.group(2)
        p = find(fname)
        if p is None:
            print(f"{doc}: FILE NOT FOUND {fname}", file=sys.stderr); return f"{fname}:?"
        rx = re.compile(pat)
        for i, line in enumerate(p.read_text().splitlines(), 1):
            if rx.search(line):
                return f"{fname}:{i}"
        print(f"{doc}: NO MATCH in {fname} for /{pat}/", file=sys.stderr)
        return f"{fname}:?"
    return re.sub(r"@@([A-Za-z0-9_.]+)\|(.+?)@@", sub, text)

def check(text, doc):
    for m in re.finditer(r"`?([A-Za-z_][A-Za-z_0-9]*\.(?:cpp|hpp|f90)):(\d+)", text):
        p = find(m.group(1)); n = int(m.group(2))
        if p is None: print(f"{doc}: {m.group(0)} -> FILE NOT FOUND"); continue
        lines = p.read_text().splitlines()
        txt = lines[n-1].strip()[:80] if n <= len(lines) else "<beyond EOF>"
        print(f"{doc}: {m.group(1)}:{n} -> {txt}")

if __name__ == "__main__":
    args = sys.argv[1:]
    mode_check = "--check" in args
    for a in args:
        if a.startswith("--"): continue
        p = pathlib.Path(a); t = p.read_text()
        if mode_check: check(t, p.name)
        else:
            out = resolve(t, p.name)
            if out != t: p.write_text(out); print(f"resolved {p.name}")
