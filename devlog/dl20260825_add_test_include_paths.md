# Add_Test.cmake: dead `SSF_SRC_DIR` include paths

**Date:** 2026-08-25
**Purpose:** Why the nonfree test suites never built, and the fix
**Module:** config/scripts (build system)

## Finding

`config/scripts/Add_Test.cmake` took the `core`, `comm`, `containers`, `linalg`,
`math/graph` and `sparse` include paths from `${SSF_SRC_DIR}` — a pre-BELFEM
name that has never been defined in this repository (present since the
initial commit `ba14a92d`). The lines expanded to `-I/core`, `-I/comm`, … and
did nothing.

The twelve open-source suites under `tests/` compiled anyway because the
root `CMakeLists.txt` builds the `banner` executable at directory scope
(line ~283) *before* `add_subdirectory( tests )` (line ~302), and
`Add_Executable.cmake`'s `include_directories` calls leak into every
subdirectory added afterwards. `nonfree/tests` is added at line ~272 —
before the leak — so `test_kepler` and `test_manta` saw only the dead
paths and died on `cl_Communicator.hpp: No such file or directory`.

ctest therefore listed both suites (tests #1 and #2 of 14) but reported
them **"Not Run"** — the executables did not exist. They never *failed*;
the 2026-08-06 devlog/todo line "kepler/manta don't compile (nonfree
include paths)" was an accurate symptom without a cause. Both records
carry an addendum pointing here.

## Fix

`Add_Test.cmake` now mirrors `Add_Executable.cmake`: `${BELFEM_SOURCE_DIR}`
replaces the dead variable, and `numerics/spline`, `numerics/bezier` and
`${CMAKE_BINARY_DIR}/generated` are added, so a test tree is
self-sufficient wherever it is added. The open-source suites see the same
paths as before (they inherited these already); nothing else changed.

## Executable gate

`make test_kepler test_manta` builds both. `ctest -R "kepler|manta"`:
kepler passes; manta 5/6 — `MantaSurface.ScaleRobustness` fails an
inverse-map round-trip by ~1e-3 (`nonfree/tests/manta/test_Surface.cpp:300`).
That is a module finding, recorded in `nonfree/devlog`, not a build issue.
`scripts/check_doc_claims.py`: 34/34 hold.
