# First Suite Run: Stale Default-Tolerance Expectation in test_sparse

**Date:** 2026-08-13
**Purpose:** Fix the one failure Christian's first-ever `USE_TEST=ON` run surfaced
in `test_sparse`, and explain the missing `test_homology` binary.
**Module:** `tests/sparse`

## The failure is not a debug/release difference

`SolverParameters.DefaultParametersValid` failed in the release-flags tree, but it
would fail identically under `-Og`: the test pins the default relative tolerance
at `1e-8`, while the source default in `cl_SolverParameters.hpp` is `1e-10`. The
history is a deliberate drift the test never followed: `1e-6` → `1e-8`
(`45264244`, "fix PETSc settings") → `1e-10` (`5aae789b`, "fixed newton
matrices"), with the test written in the first-tier batch (`eda97cec`) against
the `1e-8` value. It could never fail visibly before because the suite had never
executed — `USE_TEST` defaults OFF everywhere. This run was the first execution,
which is exactly the class of finding the gate exists to produce.

**Fix:** the expectation now pins `1e-10`, with a comment recording the drift
history and the rule that a deliberate default change updates the test in the
same commit. The neighboring `RelativeToleranceStored` test (`set` → `get`
round-trip of `1e-6`) was checked and is correct as written. Syntax-checked
against the tree's release flags (`-O2 -DNDEBUG -Werror`).

## test_homology "not built"

Not a wiring defect: the target is fully generated in the build tree
(`tests/homology/CMakeFiles/test_homology.dir` rule exists, and
`CTestTestfile.cmake` registers the `homology` test). The binary was simply not
yet compiled when ctest ran — the cmake regeneration that picked up the new
`add_subdirectory( homology )` and the build that consumed it did not complete
in the same pass. A plain `make` (or `make test_homology`) from the build
directory produces it; `make check-fast` then includes the eight homology cases.

## NDEBUG verification of both affected suites (probe tier)

Since every earlier probe of the homology suite ran under `-DDEBUG`, both
binaries were rebuilt standalone with the tree's exact release flags
(`-O2 -DNDEBUG -Wall -Werror -pedantic-errors`, full define set from the
regenerated `flags.make`) against the prebuilt libraries and run:

- `test_sparse` with the tolerance fix: **52/52**.
- `test_homology` (all eight cases): **8/8** — no `BELFEM_ASSERT` dependence in
  the fixtures, as intended.

The wider NDEBUG-cleanliness question resolved empirically: the test authors
guarded assert-dependent cases with `#ifndef NDEBUG` throughout (`test_Matrix`,
`test_DynamicBitset`, `test_Cell`, …), which is consistent with the first
release-flags run failing only on the stale sparse expectation. Note
`assert.hpp:26-28`: `BELFEM_ASSERT` is active when `!NDEBUG || DEBUG`, so the
`#ifndef NDEBUG` test guards are conservative — in a hypothetical
NDEBUG+DEBUG build they skip tests whose asserts would actually fire — which
can under-test but never false-fail.
