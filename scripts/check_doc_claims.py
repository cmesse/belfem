#!/usr/bin/env python3
"""
Check the mechanically-verifiable claims in the convention documents against the tree.

`CLAUDE.md` and `doc/coding_philosophy.md` are read at the start of every AI session and
by every new contributor.  A false statement in them is not inert — it is acted on.  The
2026-08-11 sweep of `CLAUDE.md` found nine claims that were false against the repository,
four of which would have produced wrong code: a denied MPI overload that exists, a hand-rolled
`aligned_alloc` pattern the codebase does not use, `BELFEM_ERROR` prescribed for convergence
failure, and inverted matrix-backend defaults.

Six of those nine needed no judgement to catch.  They were facts with a single source in the
tree: a compiler flag, a make target, an executable name, a module's presence in
`src/CMakeLists.txt`, a CMake default, the existence of a named symbol.  This script checks
that class and only that class.  Prose, rationale, and design advice are out of scope — they
need a reader, and that is what the cross-review round is for.

Following the same rule as `doc/input_schema.yaml`: **every probe anchors on a searchable
token, never a line number.**  A line number read from a working copy is already stale by the
time the copy is committed.

**A probe that finds nothing FAILS; it never skips.**  Three bugs of exactly that shape were
found in this script while it was being written and audited — a probe returned an empty result
and the check consuming it passed vacuously while appearing to run.  A checker whose failure
mode is silently reporting success is worse than no checker, so every fact the checks depend on
is declared in REQUIRED_FACTS and verified non-empty before anything else runs.

    python3 scripts/check_doc_claims.py [--verbose]

Exit code 0 if every claim holds, 1 otherwise.  Failures name the document, the claim, and
the tree fact that contradicts it, so the fix is obvious without re-deriving the probe.
"""

import argparse
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent

CLAUDE = 'CLAUDE.md'
PHILOSOPHY = 'doc/coding_philosophy.md'
CONVENTION_DOCS = [CLAUDE, PHILOSOPHY]

# not a convention document — it is checked for one specific thing, a prose
# count that duplicates what its own rows carry
REGISTER = 'todo/debt_register.md'

# Facts whose emptiness would silently disable a check. Probed, then asserted non-empty.
REQUIRED_FACTS = ['debug_flags', 'release_flags', 'make_targets', 'executables',
                  'library_kinds', 'share_overloads', 'options', 'apple_backend',
                  'other_backend']

# CMake is case- and space-tolerant: `if(`, `if (`, `IF (` are all the same command.
IF_RE = re.compile(r'\s*if\s*\(', re.I)
ELSEIF_RE = re.compile(r'\s*elseif\s*\(', re.I)
ELSE_RE = re.compile(r'\s*else\s*\(', re.I)
ENDIF_RE = re.compile(r'\s*endif\s*\(', re.I)


def read(relpath):
    path = ROOT / relpath
    return path.read_text() if path.exists() else ''


def cmake_branches(text, condition):
    """Return (if-branch, else-branch) of `if( <condition> )`, respecting nesting.

    Both halves of this must be depth-aware, and every form CMake accepts has to count.
    A regex that stops at the first `else()` picks up an unrelated earlier block; one that
    stops at the first `endif()` is cut short by any nested `if()`; and one that matches
    only `if(` misses the `if (` spelling used elsewhere in this repository, so its `endif`
    decrements a depth that was never incremented.  Each mistake silently empties a branch,
    which turns the check reading it into a no-op that still reports success — the reason
    REQUIRED_FACTS exists.

    `elseif` at the top level closes the if-branch, exactly as CMake reads it.
    """
    lines = text.splitlines()
    start = next((i for i, l in enumerate(lines)
                  if re.match(r'\s*if\s*\(\s*' + condition + r'\s*\)', l, re.I)), None)
    if start is None:
        return '', ''
    depth, split, end = 0, None, len(lines)
    for i in range(start, len(lines)):
        line = lines[i]
        if IF_RE.match(line):
            depth += 1
        elif ENDIF_RE.match(line):
            depth -= 1
            if depth == 0:
                end = i
                break
        elif depth == 1 and split is None and (ELSE_RE.match(line) or ELSEIF_RE.match(line)):
            split = i
    if split is None:
        return '\n'.join(lines[start + 1:end]), ''
    return '\n'.join(lines[start + 1:split]), '\n'.join(lines[split + 1:end])


def grep_tree(pattern, subdir, suffixes):
    """(count, first path) of `pattern` across a subtree."""
    rx = re.compile(pattern)
    n, where = 0, None
    for path in sorted((ROOT / subdir).rglob('*')):
        if path.suffix not in suffixes or 'cmake-build' in str(path):
            continue
        try:
            hits = rx.findall(path.read_text(errors='ignore'))
        except OSError:
            continue
        if hits and where is None:
            where = path.relative_to(ROOT)
        n += len(hits)
    return n, where


def probes():
    """Facts read out of the tree, each from its single authoritative source."""
    cmake = read('CMakeLists.txt')
    gcc = read('config/compiler/config_gcc.cmake')
    srccmake = read('src/CMakeLists.txt')
    commtools = read('src/comm/commtools.hpp')

    debug_block, release_block = cmake_branches(gcc, 'USE_DEBUG')
    facts = {
        'debug_flags': re.findall(r'BELFEM_CXXFLAGS\s+"([^"]*)"', debug_block),
        'release_flags': re.findall(r'BELFEM_CFLAGS\s+"([^"]*)"', release_block),
        'make_targets': set(re.findall(r'add_custom_target\(\s*([\w-]+)', cmake)),
        'executables': sorted(p.stem for p in (ROOT / 'src/executables').glob('*.cpp')),
        'unbuilt_modules': set(re.findall(r'^\s*#\s*add_subdirectory\(\s*(\w+)\s*\)', srccmake, re.M)),
        'vtk_gated_modules': set(re.findall(r'add_subdirectory\(\s*(\w+)',
                                            cmake_branches(srccmake, 'USE_VTK')[0])),
        'share_overloads': set(re.findall(r'^\s*share\(\s*(\w+)\s*<', commtools, re.M)),
        'options': dict(re.findall(r'option\(\s*(\w+)\s*"[^"]*"\s*(ON|OFF)\s*\)', cmake)),
    }
    facts['aligned_alloc_sites'], facts['aligned_alloc_where'] = grep_tree(
        r'\b(?:aligned_alloc|posix_memalign)\s*\(', 'src', {'.cpp', '.hpp', '.h'})
    facts['fno_exceptions_sites'], _ = grep_tree(
        r'fno-exceptions', 'config', {'.cmake', '.txt'})
    facts['fno_exceptions_sites'] += len(re.findall(r'fno-exceptions', cmake))

    # Library kinds, across every CMake file rather than the project helper alone, so
    # that a document saying "STATIC only" is checked against every exception in the
    # tree -- the tracked tree: build trees, scratch, archives and parallel checkouts
    # hold copied decks whose kinds say nothing about this project.
    not_the_project = ('cmake-build', '/build/', '/tmp/', '/archive/', '/nonfree/',
                       '/literature/', '/.claude/', '/.git/')
    facts['library_kinds'] = set()
    facts['shared_lib_files'], facts['module_lib_files'] = [], []
    for path in sorted(ROOT.rglob('*.cmake')) + sorted(ROOT.rglob('CMakeLists.txt')):
        rel = '/' + path.relative_to(ROOT).as_posix()
        if any(s in rel for s in not_the_project):
            continue
        try:
            body = path.read_text(errors='ignore')
        except OSError:
            continue
        for kind in re.findall(r'add_library\([^)]*?\b(STATIC|SHARED|MODULE)\b', body, re.S):
            facts['library_kinds'].add(kind)
            if kind == 'SHARED':
                facts['shared_lib_files'].append(rel[1:])
            elif kind == 'MODULE':
                facts['module_lib_files'].append(rel[1:])

    apple = re.search(r'if\s*\(\s*APPLE\s*\)(.*?)else\s*\(\)(.*?)endif\s*\(\)', cmake, re.S)
    if apple:
        on = lambda block, lib: bool(re.search(r'USE_MATRIX_' + lib + r'[^\n]*"\s*ON\s*\)', block))
        facts['apple_backend'] = 'Blaze' if on(apple.group(1), 'BLAZE') else 'Armadillo'
        facts['other_backend'] = 'Blaze' if on(apple.group(2), 'BLAZE') else 'Armadillo'
    return facts


def checks(f):
    """(claim, doc, holds, detail) — one row per checkable statement."""
    rows = []

    # Guard first: a fact that came back empty means a probe broke, and every check reading
    # it would pass for the wrong reason. Fail here rather than silently drop those rows.
    for name in REQUIRED_FACTS:
        rows.append((f'probe: {name}', 'scripts/check_doc_claims.py', bool(f.get(name)),
                     'probe returned nothing — the checks reading it would pass vacuously'))
    if not all(f.get(n) for n in REQUIRED_FACTS):
        return rows

    NEG = r'\b(never|not|no longer|rather than|instead of|nor|no)\b'

    def forbid(doc, pattern, claim, detail):
        """Fail if the document asserts `pattern` — but a *denial* of it is fine.

        "release is -O2 (never -O3)" and "-O3 is not used" are both correct, so the window
        is checked on both sides of the match; an earlier version looked only behind and
        would have failed the second wording.
        """
        text = read(doc)
        hits = [m for m in re.finditer(pattern, text, re.I | re.M)
                if not re.search(NEG + r'[^.\n]{0,40}$', text[max(0, m.start() - 70):m.start()], re.I)
                and not re.match(r'[^.\n]{0,40}' + NEG, text[m.end():m.end() + 70], re.I)]
        rows.append((claim, doc, not hits,
                     detail + (f' — document says {hits[0].group(0)!r}' if hits else '')))

    def require(doc, pattern, claim, detail):
        rows.append((claim, doc, bool(re.search(pattern, read(doc), re.I)), detail))

    release, debug = ' '.join(f['release_flags']), ' '.join(f['debug_flags'])

    for doc in CONVENTION_DOCS:
        # 1. optimisation flags
        forbid(doc, r'-O3\b', 'release optimisation level',
               f'config_gcc.cmake release branch sets {release!r}')
        if '-O2' in release:
            require(doc, r'-O2\b', 'release optimisation level', f'tree sets {release!r}')
        if '-Og' in debug:
            require(doc, r'-Og\b', 'debug optimisation level', f'tree sets {debug!r}')

        # 2. exceptions: probed, not assumed
        if f['fno_exceptions_sites'] == 0:
            forbid(doc, r'`?-fno-exceptions`?[^.\n]{0,40}(?:for maximum speed|release build)',
                   'exception flag', 'no CMake file sets -fno-exceptions')

        # 3. hand-rolled alignment
        if f['aligned_alloc_sites'] == 0:
            forbid(doc, r'^[^\n]*=\s*\(\w+\s*\*\)\s*aligned_alloc', 'aligned allocation',
                   'src/ contains no aligned_alloc or posix_memalign call')

        # 4. MPI overload denials must match the header
        for denied in re.findall(r'no\s+`?share\(?\s*(\w+)', read(doc)):
            rows.append((f'share({denied}) denial', doc, denied not in f['share_overloads'],
                         f'commtools.hpp declares share() for {sorted(f["share_overloads"])}'))

    # 5. test targets
    forbid(CLAUDE, r'\bmake tests\b', 'test target',
           f"CMakeLists defines {sorted(t for t in f['make_targets'] if 'check' in t)}")
    if 'check' in f['make_targets']:
        require(CLAUDE, r'\bmake check\b', 'test target', 'CMakeLists defines a check target')
    # Bidirectional: the docs must state whichever default CMake actually sets.
    # The pre-2026-08-14 form only fired when CMake said OFF, so a flip to ON
    # silenced the probe exactly when the documents became wrong.
    tTestDefault = f['options'].get('USE_TEST')
    if tTestDefault in ('ON', 'OFF'):
        tOther = 'OFF' if tTestDefault == 'ON' else 'ON'
        for doc in (CLAUDE, PHILOSOPHY):
            require(doc, r'USE_TEST[^\n]{0,30}\b' + tTestDefault + r'\b', 'test default',
                    f'CMakeLists sets option(USE_TEST ... {tTestDefault})')
            forbid(doc, r'USE_TEST`?\s+(?:(?:which\s+)?defaults\s+)\**' + tOther + r'\b',
                   'test default (stale)',
                   f'CMakeLists sets option(USE_TEST ... {tTestDefault})')

    # 6. executable inventory
    claude = read(CLAUDE)
    listed = re.search(r'\*\*executables/\*\*:([^\n]*)', claude)
    rows.append(('executable inventory', CLAUDE, bool(listed), 'no executables/ entry found'))
    if listed:
        named = set(re.findall(r'`(\w+)`', listed.group(1)))
        missing, extra = set(f['executables']) - named, named - set(f['executables'])
        rows[-1] = ('executable inventory', CLAUDE, not (missing or extra),
                    f'src/executables holds {f["executables"]}; document adds {sorted(extra)}, '
                    f'omits {sorted(missing)}')

    # 7. modules commented out of the build, and modules behind a default-OFF option
    for module in sorted(f['unbuilt_modules']):
        entry = re.search(r'\*\*' + module + r'/\*\*:([^\n]*)', claude)
        rows.append((f'{module}/ build status', CLAUDE,
                     bool(entry) and bool(re.search(r'not built|commented out|incomplete',
                                                    entry.group(1), re.I)),
                     f'src/CMakeLists.txt keeps add_subdirectory({module}) commented out'))
    if f['options'].get('USE_VTK') == 'OFF':
        for module in sorted(f['vtk_gated_modules']):
            entry = re.search(r'\*\*' + module + r'/\*\*:([^\n]*)', claude)
            rows.append((f'{module}/ USE_VTK gate', CLAUDE,
                         bool(entry) and 'USE_VTK' in entry.group(1),
                         f'{module} is built only under USE_VTK, which defaults OFF'))

    # 8. matrix backend defaults
    rows.append(('matrix backend defaults', CLAUDE,
                 bool(re.search(f['other_backend'] + r'[^\n]*default[^\n]*except[^\n]*Apple',
                                claude, re.I)),
                 f'CMakeLists defaults to {f["other_backend"]}, and to '
                 f'{f["apple_backend"]} on Apple'))

    # 9. library kind — SHARED exists, but only in the user-plugin templates
    if 'SHARED' in f['library_kinds']:
        rows.append(('library kind', CLAUDE,
                     bool(re.search(r'STATIC', claude)) and bool(
                         re.search(r'template', claude, re.I)),
                     f'project libraries are STATIC; SHARED appears only in '
                     f'{f["shared_lib_files"]}, so the document must name that exception'))
    forbid(CLAUDE, r'both static and shared library builds', 'library kind',
           'Add_Library.cmake creates STATIC libraries only')

    # 9b. plugin library kind — the templates and example decks build MODULE (dlopen'ed)
    if 'MODULE' in f['library_kinds']:
        rows.append(('plugin library kind', CLAUDE,
                     bool(re.search(r'\bMODULE\b', claude)),
                     f'plugins build MODULE libraries in {f["module_lib_files"]}, '
                     f'so the document must name that kind'))

    # 10. Doc-vs-doc, not doc-vs-tree: CLAUDE.md must not contradict the philosophy's ruling
    #     that non-convergence returns a status. This one enforces a policy rather than
    #     deriving a fact — src/ does contain BELFEM_ERROR convergence sites, which is a
    #     question about the source, not about the document.
    for doc in CONVENTION_DOCS:
        forbid(doc, r'BELFEM_ERROR[^\n]*\bconvergence failures\b', 'error tier (policy)',
               'coding_philosophy makes non-convergence a status/retry outcome, not an abort')

    # ------------------------------------------------------------------
    # todo/debt_register.md: the header sentence restates a count the rows
    # already carry, so nothing recomputes it and nothing fails when it drifts.
    # Three sessions edited that file on 2026-08-29 alone and the count went
    # stale twice. This is exactly the class this script exists for: a
    # checkable claim in prose that the tree contradicts.
    #
    # Counted on the OPEN rows only — a struck row is `| ~~DR-nn~~` and is no
    # longer part of the open population it is being compared against.
    register = read(REGISTER)
    if register:
        for tag in ('P', 'W'):
            actual = len(re.findall(
                r'^\| DR-\d+\b[^|\n]*\*\*\[' + tag + r'\]\*\*', register, re.M))

            # the header states each group as "( <n> rows" following the tag
            stated = re.search(
                r'`\[' + tag + r'\]`\s*\((\d+)\s+rows', register)

            rows.append((
                f'debt register: [{tag}] row count',
                REGISTER,
                bool(stated) and int(stated.group(1)) == actual,
                (f'header says {stated.group(1)}, tree has {actual} open [{tag}] rows'
                 if stated else
                 f'header states no count for [{tag}] — tree has {actual}')))

        # [F] is stated in a different sentence and in WORDS, not digits:
        #   "**The `[F]` list is one row**" / "**The `[F]` list is EMPTY**"
        # so the numeric parse above cannot see it. It went stale twice on
        # 2026-08-29 — once when a row was added and once when the last one
        # was struck — and both times the [P] check passed while [F] lied.
        # [F] is the freeze lens, so a wrong count there is read on exactly
        # the day it matters most.
        words = {'no': 0, 'zero': 0, 'empty': 0, 'one': 1, 'two': 2, 'three': 3,
                 'four': 4, 'five': 5, 'six': 6, 'seven': 7, 'eight': 8,
                 'nine': 9, 'ten': 10}

        actual_f = len(re.findall(
            r'^\| DR-\d+\b[^|\n]*\*\*\[F\]\*\*', register, re.M))

        stated_f = re.search(
            r'The\s+`\[F\]`\s+list\s+is\s+(?:\*\*)?([A-Za-z]+|\d+)', register)

        if stated_f:
            token = stated_f.group(1).lower()
            value = int(token) if token.isdigit() else words.get(token)
        else:
            value = None

        rows.append((
            'debt register: [F] row count',
            REGISTER,
            value is not None and value == actual_f,
            (f'header says {stated_f.group(1)!r}, tree has {actual_f} open [F] rows'
             if stated_f else
             f'header states no count for [F] — tree has {actual_f}')))

    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument('--verbose', action='store_true', help='list passing claims too')
    args = ap.parse_args()

    rows = checks(probes())
    failed = [r for r in rows if not r[2]]

    for claim, doc, holds, detail in rows:
        if args.verbose or not holds:
            print(f'{"ok  " if holds else "FAIL"} {doc}: {claim}')
            if not holds:
                print(f'       tree: {detail}')

    print(f'\n{len(rows) - len(failed)}/{len(rows)} claims hold.')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
