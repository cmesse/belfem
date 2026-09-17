# PENTA18 second derivatives

**Date:** 2026-09-16
**Purpose:** Symbolic zero-test of the hand-written 6×18 table `cl_IF_PENTA18.hpp::d2NdXi2`.
**Module:** `src/fem/interpolation`

| Script | Kind | What it does |
|---|---|---|
| `penta18.m` | zero-test | defines the eighteen PENTA18 shape functions, differentiates them twice, and subtracts a transcription of the C++ table; `X` must be all zeros |
| `check_penta18.m` | driver | runs it and asserts; `matlab -batch check_penta18` exits 0 on success |

The row order of the table is `xi xi`, `eta eta`, `zeta zeta`, `eta zeta`, `zeta xi`, `xi eta`, the same as the C++. `../compare_tables.py` checks the transcription against the header. The runtime check of the same table is the central-difference test in `tests/fem/test_LagrangeInterpolation.cpp`. Symbolic Math Toolbox required.
