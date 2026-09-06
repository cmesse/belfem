# Nédélec Volume-Element Test Battery, and the Third Defect It Caught on Arrival

**Date:** 2026-08-14
**Purpose:** Session record — the circulation/curl regression battery for the volume Nédélec
elements (Christian's ask, same day as the TET4/TET10 table fixes), and defect D3 that the
battery found within minutes of first running.
**AIs:** Claude; builds on the falsification-tooling design (R12 there) and the day's defect
work (`todo/nedelec_edge_function_defects.md`)

## What landed

- **`tests/fem/support/cl_EF_TestVolume.hpp`** (new fixture): single code-built TRI3, TRI6,
  TET4 or TET10 with caller-provided corner coordinates, exact-midpoint midside nodes,
  canonical edges from `get_nodes_of_edge()` (with a flip option), self-mastered `mesh::Face`
  objects for the TET10, explicit curved-flag control, and the minimal real
  Kernel/DofManagerBase/Block chain of the proven `TS_TestStack` pattern.
- **`tests/fem/test_EdgeFunctions.cpp`** (new, 22 tests, wired into `tests/fem/CMakeLists.txt`):
  per element type, on the reference AND a distorted element (the distorted vertices are the
  exact rational coordinates of the day's sympy probes, so both verification chains meet on
  the same numbers): circulation identity (own-edge value 1, TRI6 pair 1/2, TET10 pair 1,
  zero on foreign edges, zero for all face dofs); `C` against central-difference curl in
  physical space (exact for these polynomial degrees, test-side Jacobian from the pinned node
  map); EXODUS `detJ > 0` (pins the node-map booby trap directly); edge-flip sign isolation;
  TET10 curved/straight path equivalence. The pinned conventions in the test double as
  documentation: the TET parameter table (node 1 carries zeta) is written out with the trap
  warning.
- The pre-existing Lagrange suite (`test_LagrangeInterpolation.cpp`, Kronecker delta,
  partition of unity, FD first/second derivatives, 24 element types) was run alongside:
  Christian's ask for "tests for both" was half-fulfilled already; the battery completes the
  Nédélec half.

**Result: 46/46 green** (24 Lagrange + 22 Nédélec) in a standalone probe against the prebuilt
libraries with the freshly compiled fixed `cl_EF_TET10.cpp`. Negative control: linked against
the pre-fix stale library, the battery reproduces the D2 circulation signature (2/3, −1)
exactly — it demonstrably catches the class it was built for.

## D3 — found by the battery on its first full run

The TET10 curved-path tests returned values of order 1e14. Valgrind was clean (computed
garbage, not uninitialized reads); a private-access probe printing the Jacobian entries showed
the mechanism: **the eta- and zeta-columns of the Jacobian were identical** (b=c, e=f, h=i),
so J was singular and `1/detJ` amplified the nablas. Root cause: rows 8 and 9 of the
`mNeta`/`mNzeta` shape-derivative tables in `EF_TET10::precompute()` (the ζτ and ητ midside
nodes) carried each other's values — the third instance of the naive-node-map trap, four
entries. `mNxi` and the `mD` second-derivative table were checked row-by-row against the
pinned shape set in `tmp/tet10/tet10_lagrange.m` (whose node table is correct) and are clean.
Fixed in place; battery green afterwards.

**The masking lesson, now baked into the battery:** `link()` evaluates its Jacobian at
evaluation point 0. With point 0 on the first edge (eta = 0), the crossed entries coincide
with the correct values, so the straight path looked healthy all along; only the curved path
(per-point Jacobian) exposed the defect. The battery now places an interior point at column 0
so both paths always see a fully generic Jacobian.

## Status

All three defects of the day (D1 TET4 edge term, D2 TET10 scalar/derivative tables, D3 TET10
shape-derivative rows) are fixed in source and covered by the battery. The probe run is the
evidence so far; the formal gate is `make check` in a `USE_TEST=ON` tree after Christian's
rebuild, plus the 3D bulk-conductor run. Tracker: `todo/nedelec_edge_function_defects.md`
(R3, in progress); battery record: `todo/falsification_tooling.md` R12.
