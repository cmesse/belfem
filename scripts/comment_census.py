#!/usr/bin/env python3
# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""Comment census over the BELFEM sources.

Counts, per module and optionally per file, the comment classes that
doc/commenting_guidelines.md asks to remove or to review: banners, sub-banners,
narration candidates, commented-out statements, TODO/FIXME, dated lines, review
provenance, pointers into the working record, and Doxygen @brief blocks that
restate a name. Compiler directives (!$omp, !dir$, #pragma) and string literals
are excluded; the six units of the closed cohomology core are counted in a
separate row and never listed as targets.

Usage:
    scripts/comment_census.py [--files [N]] [--root src] [--closed]

    --files [N]   also print the N worst files (default 25)
    --root DIR    directory to walk (default: src, relative to the repo root)
    --closed      list the per-file counts of the closed core instead

The "noise" column is a first-word pattern on short comment lines and is an
over-approximation: a stage heading or a symbol annotation matches it too.
Every candidate is read before it is removed; the number is for the devlog's
before/after table, not a deletion list.
"""
import argparse
import collections
import os
import re
import sys

NOISE = re.compile(
    r'^\s*(//|!)\s*(increment|initiali[sz]e|reset|allocate|deallocate|free|delete|'
    r'loop over|get my|get the|set the|tidy up|clean up|wait for|wait until|'
    r'create (the )?communicator|return|copy|call|compute|update|write|read|open|'
    r'close|check|add|remove|clear|fill|assign|save|load|send|receive|resize|'
    r'populate|start|end|done|now|finally|first|then|grab|set proc|reset counters?)\b',
    re.I)
BANNER = re.compile(r'^\s*(//|!)\s*[-=]{8,}\s*$')
SUBBAN = re.compile(r'^\s*(//|!)(\s-){6,}\s*$')
DEAD = re.compile(r'^\s*//\s*[A-Za-z_][\w:<>.\-]*\s*(\(|=|\+\+|--).*;\s*$')
TODO = re.compile(r'\b(TODO|FIXME)\b', re.I)
DATED = re.compile(r'(//|^\s*!).*20[0-9]{2}-[0-9]{2}-[0-9]{2}')
PROV = re.compile(r'(//|^\s*!).*\b(codex|grok|claude)\b', re.I)
RECORD = re.compile(r'(//|^\s*!).*(ai_exchange|\bINC-[0-9]+|\bDR-[0-9]+|\bL-[0-9]{2}\b)')
BRIEF = re.compile(
    r'@brief\s+(default |copy )?(constructor|destructor|get|set|compute|check|create|return)\b',
    re.I)
STRLIT = re.compile(r'"[^"\\]*(?:\\.[^"\\]*)*"')
DIRECTIVE = re.compile(r'^\s*!\$|^\s*!dir\$|^\s*#pragma', re.I)

CLOSED_CORE = {'cl_Cohomology', 'cl_Homology', 'cl_SimplicialComplex',
               'cl_Chain', 'cl_Cochain', 'fn_Smith'}
EXTENSIONS = ('.cpp', '.hpp', '.h', '.c', '.f90')
COLUMNS = ['lines', 'comment', 'banner', 'subban', 'noise', 'dead', 'todo',
           'dated', 'prov', 'record', 'brief?']


def census_file(path):
    is_fortran = path.endswith('.f90')
    is_comment = re.compile(r'^\s*!') if is_fortran else re.compile(r'^\s*//')
    counts = [0] * len(COLUMNS)
    with open(path, errors='ignore') as handle:
        lines = handle.read().split('\n')
    counts[0] = len(lines)
    for line in lines:
        if DIRECTIVE.match(line):
            continue
        stripped = STRLIT.sub('""', line)
        commented = bool(is_comment.match(line))
        has_comment = commented or '//' in stripped or (is_fortran and '!' in stripped)
        counts[1] += commented
        counts[2] += bool(BANNER.match(line))
        counts[3] += bool(SUBBAN.match(line))
        counts[4] += bool(NOISE.match(line) and len(line.strip()) < 45)
        counts[5] += bool(DEAD.match(line))
        counts[6] += bool(has_comment and TODO.search(stripped))
        counts[7] += bool(DATED.search(stripped))
        counts[8] += bool(PROV.search(stripped))
        counts[9] += bool(RECORD.search(stripped))
        counts[10] += bool(BRIEF.search(line))
    return counts


def module_of(rel):
    parts = rel.split(os.sep)
    return os.sep.join(parts[:2]) if len(parts) > 2 else parts[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--files', nargs='?', const=25, type=int, default=0)
    parser.add_argument('--root', default='src')
    parser.add_argument('--closed', action='store_true')
    args = parser.parse_args()

    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    root = os.path.join(repo, args.root)
    if not os.path.isdir(root):
        sys.exit(f'comment_census: no such directory {root}')

    rows = []
    for dirpath, _, files in os.walk(root):
        for name in files:
            if not name.endswith(EXTENSIONS):
                continue
            path = os.path.join(dirpath, name)
            rel = os.path.relpath(path, root)
            closed = os.path.splitext(name)[0] in CLOSED_CORE
            rows.append((rel, closed, census_file(path)))

    if args.closed:
        for rel, closed, counts in sorted(rows):
            if closed:
                print(f'{rel:50s}' + ''.join(f'{x:7d}' for x in counts))
        return

    per_module = collections.defaultdict(lambda: [0] * len(COLUMNS))
    total = [0] * len(COLUMNS)
    closed_total = [0] * len(COLUMNS)
    for rel, closed, counts in rows:
        for i, x in enumerate(counts):
            total[i] += x
            per_module[module_of(rel)][i] += x
            if closed:
                closed_total[i] += x

    header = f'{"module":28s}' + ''.join(f'{c:>8s}' for c in COLUMNS)
    print(header)
    print('-' * len(header))

    def weight(counts):
        return counts[4] + counts[5] + counts[7] + counts[8] + counts[9]

    for module, counts in sorted(per_module.items(), key=lambda kv: -weight(kv[1])):
        if counts[0] < 300:
            continue
        print(f'{module:28s}' + ''.join(f'{x:8d}' for x in counts))
    print('-' * len(header))
    print(f'{"TOTAL":28s}' + ''.join(f'{x:8d}' for x in total))
    print(f'{"closed core (6 units)":28s}' + ''.join(f'{x:8d}' for x in closed_total))

    if args.files:
        print(f'\nworst {args.files} files by noise + dead + dated + prov + record '
              '(closed core excluded):')
        worst = sorted((r for r in rows if not r[1]), key=lambda r: -weight(r[2]))
        for rel, _, c in worst[:args.files]:
            print(f'{rel:62s} {c[0]:6d}  noise {c[4]:4d}  dead {c[5]:3d}  '
                  f'todo {c[6]:3d}  dated {c[7]:3d}  prov {c[8]:3d}  banner {c[2]:4d}')


if __name__ == '__main__':
    main()
