#!/usr/bin/env python3
# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""Check the MATLAB transcriptions of the C++ tables against the C++ files.

The zero-tests in nedelec_tet10/ and lagrange/ compare the generator against
tables typed into the MATLAB scripts, not against the C++ source. This script
closes that gap: it parses every `mX( i, k ) = expr;` line of
cl_EF_TET10.cpp and every `ad2NdXi2( i, j ) = expr;` line of
cl_IF_PENTA18.hpp, parses the corresponding MATLAB tables, and reports every
entry whose symbolic difference is not zero. The MATLAB face tables hold the
eight active rows ( act = [1 2 4 5 7 8 10 11] ) of the twelve C++ rows.

Run from anywhere:  python3 compare_tables.py   (exit 0 = no differences)
"""
import re, sys
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
SRC = HERE.parents[1]                      # src/fem/interpolation
xi, eta, zeta = sp.symbols('xi eta zeta', real=True)
tau = 1 - xi - eta - zeta
ns = {'xi': xi, 'eta': eta, 'zeta': zeta, 'tau': tau}
for name, sym in (('xi', xi), ('eta', eta), ('zeta', zeta), ('tau', tau)):
    for f in (2, 4, 8, 16, 32):
        ns[f'{name}{f}'] = f * sym


def parse(expr):
    expr = re.sub(r'(\d)\.(?!\d)', r'\1', expr)      # C++ `1.` -> `1`
    return sp.sympify(expr.replace('^', '**'), locals=ns)


def cpp_tables(path, pattern):
    tables = {}
    for name, i, j, expr in re.findall(pattern, path.read_text()):
        tables.setdefault(name, {})[(int(i), j)] = parse(expr)
    return tables


def matlab_columns(path):
    """`name = [e1\ne2\n...];` blocks -> {name: [expr, ...]}"""
    cols = {}
    for name, body in re.findall(r'^(m\w+)\s*=\s*\[(.*?)\];', path.read_text(), re.S | re.M):
        cols[name] = [parse(l.strip()) for l in body.strip().splitlines() if l.strip()]
    return cols


def main():
    bad = 0
    # TET10 Nedelec tables
    cpp = cpp_tables(SRC / 'nedelec' / 'cl_EF_TET10.cpp',
                     r'(m[GHUVW](?:xi|eta|zeta)?)\(\s*(\d+),\s*(k)\s*\)\s*=\s*([^;]+);')
    ml = {}
    ml.update(matlab_columns(HERE / 'nedelec_tet10' / 'tet10_function.m'))
    ml.update(matlab_columns(HERE / 'nedelec_tet10' / 'tet10_derivatives.m'))
    act = [0, 1, 3, 4, 6, 7, 9, 10]
    for name, col in sorted(ml.items()):
        rows = range(12) if name[1] in 'GH' else act
        if len(col) != len(rows):
            print(f'{name}: MATLAB has {len(col)} rows, expected {len(rows)}'); bad += 1; continue
        for m, c in zip(col, rows):
            d = sp.expand(m - cpp[name][(c, 'k')]) if (c, 'k') in cpp[name] else None
            if d is None:
                print(f'{name}({c}): missing in C++'); bad += 1
            elif d != 0:
                print(f'{name}({c}): MATLAB {m}  C++ {cpp[name][(c, "k")]}  diff {d}'); bad += 1
    n_tet = sum(len(v) for v in ml.values())
    # PENTA18 second derivatives: MATLAB is one-based, C++ zero-based
    hdr = (SRC / 'lagrange' / 'cl_IF_PENTA18.hpp').read_text()
    cpp18 = {(int(i), int(j)): parse(e) for i, j, e in
             re.findall(r'ad2NdXi2\(\s*(\d+),\s*(\d+)\s*\)\s*=\s*([^;]+);', hdr)}
    ml18 = {(int(i) - 1, int(j) - 1): parse(e) for i, j, e in
            re.findall(r'ad2NdXi2\(\s*(\d+),\s*(\d+)\s*\)\s*=\s*([^;]+);',
                       (HERE / 'lagrange' / 'penta18.m').read_text())}
    if set(cpp18) != set(ml18):
        print(f'PENTA18: index sets differ ({len(cpp18)} vs {len(ml18)})'); bad += 1
    for k in sorted(set(cpp18) & set(ml18)):
        d = sp.expand(cpp18[k] - ml18[k])
        if d != 0:
            print(f'ad2NdXi2{k}: MATLAB {ml18[k]}  C++ {cpp18[k]}  diff {d}'); bad += 1
    print(f'compare_tables: {n_tet} TET10 entries and {len(ml18)} PENTA18 entries compared, {bad} difference(s)')
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
