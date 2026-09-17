#!/usr/bin/env python3
# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""Strip comments from a C/C++ or Fortran source, keeping everything else.

Used by scripts/check_comment_only.sh: two revisions of a file are stripped
and diffed, so a comment-only change leaves nothing behind. The stripper
therefore has to be exact about what is *not* a comment:

  C/C++    - string literals ("...", with escapes), character literals
             ('...'), raw string literals (R"delim(...)delim"), and
             preprocessor lines are copied through untouched, including a
             '#' on a macro continuation line, which is what makes
             `g++ -fpreprocessed` unusable here.
           - `//` to end of line and `/* ... */` (possibly spanning lines)
             are removed.
  Fortran  - `!` to end of line outside '...' and "..." literals is removed;
             `!$` directives (OpenMP, and `!dir$`-style vendor directives)
             are kept because the compiler reads them.

Blank lines and trailing whitespace are dropped from the output so that a
removed comment line does not leave a difference behind.

Usage:  strip_comments.py {c,fortran} < in > out
"""
import sys


def strip_c(text):
    out = []
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        # raw string literal R"delim( ... )delim", with an optional u8/u/U/L
        # prefix; the R must not be the tail of a longer identifier
        if c == 'R' and text[i + 1:i + 2] == '"':
            k = i
            if k > 0 and text[k - 1] in 'uUL':
                k -= 1
                if k > 0 and text[k - 1:k + 1] == 'u8':
                    k -= 1
            if k == 0 or not (text[k - 1].isalnum() or text[k - 1] == '_'):
                j = text.find('(', i + 2)
                if j != -1:
                    delim = text[i + 2:j]
                    end = text.find(')' + delim + '"', j)
                    if end != -1:
                        out.append(text[i:end + len(delim) + 2])
                        i = end + len(delim) + 2
                        continue
        if c == '"' or c == "'":
            # string or character literal; copy with escapes
            j = i + 1
            while j < n and text[j] != c:
                if text[j] == '\\':
                    j += 1
                if text[j:j + 1] == '\n':
                    break          # unterminated on this line: give up on it
                j += 1
            out.append(text[i:j + 1])
            i = j + 1
            continue
        if c == '/' and i + 1 < n and text[i + 1] == '/':
            j = text.find('\n', i)
            i = n if j == -1 else j       # keep the newline
            continue
        if c == '/' and i + 1 < n and text[i + 1] == '*':
            j = text.find('*/', i + 2)
            if j == -1:
                i = n
                continue
            # keep the newlines inside the block so later line-based
            # diffs stay aligned; they are dropped as blank lines anyway
            out.append('\n' * text[i:j + 2].count('\n'))
            i = j + 2
            continue
        out.append(c)
        i += 1
    return ''.join(out)


def strip_fortran(text):
    out = []
    for line in text.split('\n'):
        quote = ''
        res = []
        i, n = 0, len(line)
        while i < n:
            c = line[i]
            if quote:
                res.append(c)
                if c == quote:
                    quote = ''
            elif c in ('"', "'"):
                quote = c
                res.append(c)
            elif c == '!':
                if line[i:i + 2] == '!$' or line[i:i + 5].lower() == '!dir$':
                    res.append(line[i:])
                break
            else:
                res.append(c)
            i += 1
        out.append(''.join(res))
    return '\n'.join(out)


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in ('c', 'fortran'):
        sys.exit(__doc__)
    text = sys.stdin.read()
    stripped = strip_c(text) if sys.argv[1] == 'c' else strip_fortran(text)
    lines = [l.strip() for l in stripped.split('\n')]
    sys.stdout.write('\n'.join(l for l in lines if l) + '\n')


if __name__ == '__main__':
    main()
