# TET10 Nédélec tables

**Date:** 2026-09-16
**Purpose:** Generator and zero-tests of the edge and face function tables in `cl_EF_TET10.cpp`, plus the Lagrange derivation for `cl_IF_TET10.hpp`.
**Module:** `src/fem/interpolation`

The C++ tables were produced by `tet10_generate.m`; `tet10_function.m` and `tet10_derivatives.m` hold a transcription of them and check it against the generator symbolically. All scripts are MATLAB scripts (shared workspace), not functions, so run order matters.

| Script | Kind | What it does | Needs |
|---|---|---|---|
| `tet10_generate.m` | generator | builds `g, h` (12 edge functions) and `u, v, w` (12 face candidates, three per face) under the pinned node map `lambda = [xi; zeta; eta; tau]` | Symbolic Toolbox |
| `tet10_shortcuts.m` | support | the `xi2 … tau8` shorthand used in the C++ tables | after the generator |
| `tet10_function.m` | zero-test | compares `mG, mH` (12 rows each) and the 8 active rows of `mU, mV, mW` against the generator; prints five `expand(...)` results that must be zero | runs the two above |
| `tet10_derivatives.m` | zero-test | the fifteen derivative tables `mGxi … mWzeta` | runs the two above |
| `check_tet10.m` | driver | runs both zero-tests and asserts every residual; `matlab -batch check_tet10` exits 0 on success | |
| `parse_main.m` | support | prints one table (`mWzeta`) in C++ syntax; the pattern for regenerating another | |
| `defelement.m` | cross-check | the same construction against the published DefElement expressions of the degree-2 Nédélec tetrahedron; every `expand(...)` must print zero; some pairs match with a sign flip and say so | standalone |
| `tet10_lagrange.m` | derivation note | derives `cl_IF_TET10.hpp::N` by inverting a Vandermonde matrix over the Pascal basis in EXODUS node order; the plot at the end is a sketch | standalone |

**The node-map pin.** EXODUS numbers the TET10 so that node 1 carries `zeta` and node 2 carries `eta`. The naive triangle-extension map `[xi; eta; zeta; tau]` is left-handed and yields tables that fail edge conformity; the generator pins `[xi; zeta; eta; tau]` and `defelement.m` confirms it against an independent source. Anyone regenerating a table must keep that line.

**Coverage.** The generator emits three face candidates per face; the C++ stores all twelve rows of `mU, mV, mW`; the zero-tests compare the eight active ones (`act = [1 2 4 5 7 8 10 11]`). The C++ `precompute()` also fills the Lagrange tables `mNxi, mNeta, mNzeta`; those come from `tet10_lagrange.m`, not from the generator. `../compare_tables.py` checks the MATLAB transcriptions against `cl_EF_TET10.cpp` itself.
