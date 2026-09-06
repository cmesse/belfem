# DR-42 — gesvd provided-buffer tests and spline column-overload tests

**Date:** 2026-08-11
**Purpose:** Close the two coverage gaps the 2026-08-05 jury round raised as P1s
**Modules:** `linalg/lapack`, `numerics/spline`, `tests/linalg`, `tests/math`

## What the row asked for

DR-42 carried two Grok-only findings from the gesvd/posv/Spline jury round of
2026-08-05, both Claude-verified at the time and both routed to Christian:

1. the gesvd wrapper's **provided-buffer branch is untested**, complex flavors
   especially;
2. the spline **`aCol` overloads are untested**.

Both are now written. Neither has been through `make check-fast` — that is
Christian's, and it is the only thing left on the row.

## gesvd: which branch had never run

`fn_gesvd.hpp` splits on the length of the caller's `Work`:

- **too short** → ask LAPACK for the optimal size with `lwork = -1`, resize,
  then decompose;
- **long enough** → trust the buffer and derive `lwork` from its length.

Both existing tests (`GesvdReconstruct`, `GesvdEconomy`) enter with a default
constructed `Vector`, so both take the first branch every time. The second was
dead in the suite.

That branch is where the arithmetic is delicate. There is one caller-owned
real-valued buffer; LAPACK's `work` lives in its head — for complex `T`
reinterpreted in place, two reals per entry — and the complex-only
`5*min(m,n)` real `rwork` in its tail. An off-by-one in either the length
division or the tail offset does not fail loudly, it overruns.

Four new typed tests, each running over `float`, `double`,
`complex<float>`, `complex<double>`:

| test | what it pins |
|---|---|
| `GesvdProvidedWorkExactMinimum` | `Work` at exactly the length the wrapper calls sufficient. Asserts the length is **unchanged** after the call — the query branch resizes, so an unchanged length is the observable proof the other branch ran |
| `GesvdProvidedWorkOversized` | minimum + 137. The tail moves far from where the minimum puts it, and for complex the split lands at `2*lwork < length - rwork`, which the exact-minimum case cannot reach |
| `GesvdWorkBufferReuse` | one buffer across three calls — query-grown, reused at the same shape, then reused on a **smaller** problem, where `lwork` comes from the old length while the `rwork` tail is sized from the new `min(m,n)` |
| `GesvdSingularValuesOnly` | `jobu = jobvt = 'N'`: `U`/`VT` are never sized, so the wrapper passes `ldu = ldvt = 1` and the `data()` of empty matrices. Cross-checked against the full decomposition |

All four pass `AbortOnError = false` and assert on the returned `info`. That is
deliberate, and it is what makes the exact-minimum test worth having: the test
mirrors the wrapper's own `lwork` formula, so it is not a second opinion on
what LAPACK needs — **LAPACK gives that verdict itself.** If the formula
understates the true minimum, `?gesvd` returns `info = -13` and the test goes
red instead of the process aborting.

## spline: the overload that was being used as its own reference

§3.4 of `test_Spline.cpp` documents the deliberate extrapolation behaviour by
asserting `eval( x ) == eval( x, 0 )` for `x < x_min`. It uses the column
overload as the trusted reference — while nothing tested the column overload.

The contract has two halves, and only the first can be checked without
trusting the class:

1. `eval( x, k )` is the polynomial of column `k` at `x`;
2. `eval( x )` is `eval( x, find_col( x ) )`.

New §8 covers both, for all five overloads (`eval`, `deval`, `ddeval`,
`entropy`, `dentropy`):

- `EvalWithColumnUsesGivenColumn` / `EntropyWithColumnUsesGivenColumn` —
  every valid column at three abscissae: the left knot, the midpoint, and a
  point far outside the interval. The last one is what separates "uses column
  `k`" from "clamps the way `find_col` does". Expected values come from
  coefficients written by hand, so they are known in closed form.
- `ColumnOverloadsAgreeWithFindCol` — `EXPECT_DOUBLE_EQ`, not `EXPECT_NEAR`:
  the two forms are the same arithmetic in the same order, so anything short
  of identity is a defect.
- `SplineColumnDebug` — the `aCol < mNumberOfIntervals` range assert (the
  coefficient table is `n()` columns wide but only `n()-1` of them are
  intervals, so the boundary is one short of the table) and the entropy-mode
  assert on the column forms.

### Why §8 is outside the SuperLU gate

The table is written straight through `matrix_data()` on a spline built with
the empty-container constructor. No solver is in the loop, so §8 compiles and
runs in every configuration. DR-49's whole premise was that the construction
sections were gated out and these tests would have been invisible with them;
building them solver-free removes the dependence rather than betting on the
gate.

While in the file, its header comment was corrected — it still claimed
"Requires UMFPACK", which DR-49 had already made false.

## Evidence

Static, both files:

- clean compile under **Armadillo** and **Blaze**, each in `NDEBUG` and
  `DEBUG`, with the tree's own `-Wall -Werror -pedantic-errors`.

Executed, as standalone scratchpad probes carrying the same assertions
(the suite itself needs a build, which is Christian's):

- **gesvd**, 4 flavors × 4 cases, all green against LAPACK. The complex
  exact-minimum case is the informative one: `lwork = 10` with a 15-entry
  `rwork` tail, buffer length 35, `info = 0` — LAPACK accepted the wrapper's
  minimum. Reconstruction error `1.4e-15` (`complex<double>`), `9.5e-07`
  (`complex<float>`).
- **spline**, all five column overloads exact to `0.000e+00` against the
  closed form, column-free vs column form **bitwise** identical at all seven
  abscissae, and 7/7 asserts fired.

Caveat on the spline probe: it linked the prebuilt (NDEBUG, Blaze)
`libbelfem_spline.a` for the constructor while compiling the header-inline
evaluators with `-DDEBUG` in its own translation unit. That is what let the
debug asserts be exercised at all; the library side contributes only scalar
setup, so the mix does not affect what was measured. It is a probe, not a
substitute for `make check-fast`.

## Left open

`make check-fast` with `USE_TEST=ON`. Nine new test cases, of which the four
gesvd ones instantiate over four types — so 16 + 5 assertions' worth of newly
executing code. A red result there is a finding, not a regression.

## Side observation, not a defect

The stale April build tree (`build/`) is configured with `USE_SUITESPARSE` and
without `USE_SUPERLU`, so `test_Spline.cpp`'s §3–§7 would compile out there.
`USE_SUPERLU` is `ON` by default (`CMakeLists.txt:69`) and the live tree
(`cmake-build-debug/`) defines `BELFEM_SUPERLU`, so DR-49's re-gating holds
for any default configure. Noted only so a future reader who opens the old
tree does not read it as a regression.
