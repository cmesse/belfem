#!/usr/bin/env python3
"""Put the BELFEM license header on every source file and drop CLion bylines.

Walks the tree (nonfree/, literature/, tmp/, archive/, build trees and .git are
skipped), and for every .cpp/.hpp/.h/.c/.f90/.F90/.py/.m file:

  * removes an IDE byline block: the `// Created by <name> on <date>.` line
    with the bare `//` lines directly around it, wherever it sits in the first
    thirty lines (the block is repeated in a few files);
  * inserts the license block at the top when the file does not carry one,
    after a shebang and an encoding line in python.

Comment-only by construction. `--dry-run` lists what would change; `--apply`
writes. The C++ gate is scripts/check_comment_only.sh; the python gate is the
`ast` comparison this script performs itself on every file it rewrites.
"""
import ast, re, sys
from pathlib import Path

MARK = 'BELFEM -- The Berkeley Lab Finite Element Framework'
BODY = [
    'BELFEM -- The Berkeley Lab Finite Element Framework',
    'Copyright (c) 2026, The Regents of the University of California,',
    'through Lawrence Berkeley National Laboratory (subject to receipt of any required',
    'approvals from the U.S. Dept. of Energy).  All rights reserved.',
    '',
    'Developers: Christian Messe, Gregory Giard',
    '',
    'See the top-level LICENSE file for the complete license and disclaimer.',
]
C_HEADER = '/*\n' + ''.join((' * ' + l).rstrip() + '\n' for l in BODY) + ' */\n'
HASH_HEADER = ''.join(('# ' + l).rstrip() + '\n' for l in BODY)
BANG_HEADER = ''.join(('! ' + l).rstrip() + '\n' for l in BODY)
PERCENT_HEADER = ''.join(('% ' + l).rstrip() + '\n' for l in BODY)

SKIP = {'nonfree', 'literature', 'tmp', 'archive', '.git'}
C_EXT = {'.cpp', '.hpp', '.h', '.c'}
F_EXT = {'.f90', '.F90'}
P_EXT = {'.py'}
M_EXT = {'.m'}
BYLINE = re.compile(r'^\s*//\s*Created by .* on .*$')


def strip_byline(lines):
    """Remove every byline block in the first thirty lines; return (lines, count)."""
    n = 0
    i = 0
    while i < min(len(lines), 30):
        if BYLINE.match(lines[i]):
            lo = i
            while lo > 0 and lines[lo - 1].strip() == '//':
                lo -= 1
            hi = i + 1
            while hi < len(lines) and lines[hi].strip() == '//':
                hi += 1
            # one blank line that separated the block from what follows
            if hi < len(lines) and lines[hi].strip() == '' and (lo == 0 or lines[lo - 1].strip() == ''):
                hi += 1
            del lines[lo:hi]
            n += 1
            i = lo
        else:
            i += 1
    return lines, n


def add_header(lines, ext):
    if ext in C_EXT:
        return C_HEADER.splitlines(True) + ([] if lines and lines[0].strip() == '' else ['\n']) + lines
    if ext in F_EXT:
        return BANG_HEADER.splitlines(True) + ['\n'] + lines
    if ext in M_EXT:
        return PERCENT_HEADER.splitlines(True) + ['\n'] + lines
    # python: keep shebang and encoding line first
    k = 0
    if lines and lines[0].startswith('#!'):
        k = 1
    if k < len(lines) and re.match(r'^#.*coding[:=]', lines[k]):
        k += 1
    return lines[:k] + HASH_HEADER.splitlines(True) + ['\n'] + lines[k:]


def process(path, apply):
    ext = path.suffix
    text = path.read_text(encoding='utf-8')
    lines = text.splitlines(True)
    lines, nb = strip_byline(lines) if ext in C_EXT else (lines, 0)
    added = False
    if MARK not in text:
        lines = add_header(lines, ext)
        added = True
    if not nb and not added:
        return None
    new = ''.join(lines)
    if ext in P_EXT:
        if ast.dump(ast.parse(text)) != ast.dump(ast.parse(new)):
            sys.exit(f'python AST changed: {path}')
    if apply:
        path.write_text(new, encoding='utf-8')
    return nb, added


def main():
    apply = '--apply' in sys.argv
    root = Path('.')
    files = []
    for p in root.rglob('*'):
        if any(part in SKIP or part.startswith('build') or part.startswith('cmake-build') for part in p.parts):
            continue
        if p.is_file() and p.suffix in C_EXT | F_EXT | P_EXT | M_EXT:
            files.append(p)
    nb = na = 0
    for p in sorted(files):
        r = process(p, apply)
        if r:
            b, a = r
            nb += b; na += a
            print(f'{"byline " if b else "       "}{"header " if a else "       "} {p}')
    print(f'{"applied" if apply else "dry run"}: {nb} byline block(s) removed, {na} header(s) added, {len(files)} files scanned')


if __name__ == '__main__':
    main()
