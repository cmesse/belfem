# `make install` and the shared `libbelfem` on Linux — R4 and R9 gated

**Date:** 2026-08-27
**Purpose:** Session record for the Linux half of `todo/shared_library_and_install_plan.md`. The
Darwin session the same day left R4, R9(a) and R9(b) owed; all three now have evidence.
**Module:** `CMakeLists.txt`, `config/` (no edits — this session was a gate, not a change)
**Platform:** el9.8 x86_64, GCC 11 (system), CMake 4.3.3, Armadillo backend, MKL/SCLS, Debug,
Open MPI, 20 cores

## What ran

The tree already carried R1–R8b from the Mac. Nothing in `src/` or `config/` was edited here. The
build tree `cmake-build-debug` was reconfigured in place — a second tree was not an option, the
filesystem is at 94% and that tree alone is 56 GB:

```
cmake -DUSE_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=<scratch>/prefix .
make -j20
```

Flipping the library type recompiles nothing: the modules are OBJECT libraries and were already
PIC, so only the aggregate library, the six executables and the test binaries relink. That took
about two minutes.

## Results

| Gate | Result |
|---|---|
| R4 artefacts | `lib/libbelfem.so.0.9.0` (287 MB Debug) → `.so.0` → `.so`, all three present |
| R4 `ldd` | clean on `banner`, `belfem`, `hphirun`, `hphiTrun`, `electricalCircuit`, `material` |
| R4 runtime | `bin/banner` runs with `LD_LIBRARY_PATH` unset |
| R4 `make check` | 12 of 13 — see below |
| R6/R7 install | `bin` 18 MB, `include` 6.1 MB, `lib64` 275 MB, `share` 92 MB; `Allrun` keeps its execute bit |
| R7 churn | a file dropped into `share/material` appears in the prefix after `make install` with no re-configure |
| R9(b) | installed `hphirun` runs the shipped `2D_Tapestack` deck with `BELFEM_DATA` **unset** |
| R9(a) | `DESTDIR` staged, tree moved: `$ORIGIN/../lib64` first in `RUNPATH`, `libbelfem.so.0` resolved out of the moved `lib64/`, deck runs with `BELFEM_DATA` set |
| O8 | the moved tree still fails with `BELFEM_DATA` unset — identical to Darwin |

The executables shrank from ~250–300 MB each to ~4 MB, which is what makes an installed tree
reasonable to ship: `lib64` carries the debug information once instead of six times.

`GNUInstallDirs` resolves to `lib64` here, as expected on el9, and the install destination and the
`$ORIGIN`-relative rpath are both derived from it, so they cannot disagree. That was G8's whole
point and it is now exercised rather than argued.

## The one test failure

`containers/Genome.RandomizeLogScaleBranch` fails:

```
tests/containers/test_Genome.cpp:287: Expected: (tOutput(1)) <= (tMax(1) + BELFEM_EPSILON),
                                     actual: 100.00000000000013 vs 100
```

An absolute epsilon on a quantity of magnitude 100: the log-scale round trip lands 1.3e-13 over the
upper bound. This is not a linkage effect and cannot be one — the static and the shared library are
assembled from the *same* object files, so the machine code behind the test is byte-identical either
way.

**Fixed the same session, on Christian's word.** The defect is in the class, not the test:
`Genome::set_values()` and `Genome::randomize()` both clamp to `[min, max]` on the way in, and
`get_values()` did not clamp on the way out, so `exp( log( tMin ) + 1·( log( tMax ) − log( tMin ) ) )`
handed back a few ulp more than the caller's maximum. `src/containers/cl_Genome.hpp` now clamps the
decoded value against `mMinVals( i )` / `mMaxVals( i )`, which makes the bounds an invariant rather
than a near-miss. Six lines, no interface change, and `Genome` has no consumer outside its own test
yet — `grep` over `src/`, `nonfree/` and `tests/` finds only `tests/containers/test_Genome.cpp`.

The gate was differential and standalone, because the build tree was busy: a scratchpad probe links
`cl_Genome.hpp` against `assert.cpp` and `stringtools.cpp` only, and was run twice — against a copy
of the header with the clamp removed, where it reproduces `100.00000000000013` exactly as ctest
reported it, and against the fixed header, where the same axis returns `100` and ten `randomize()`
trials stay inside their bounds. The existing test is the regression test; nothing new was added.

## Two things worth knowing

**The `::` in the build-tree RUNPATH is CMake's, not ours.** A build-tree binary carries

```
/usr/lib64:/opt/scls/mkl/lib64:/opt/intel/…/lib/intel64:<build>/lib::/opt/scls/mkl/lib
```

and an empty entry means "the current directory" to the loader. It looked like a defect in the
fifteen config files that append to `BELFEM_RPATH`. It is not: `BELFEM_RPATH`, `CMAKE_BUILD_RPATH`
and the directory `LINK_DIRECTORIES` property were each printed and are all free of empty entries,
and a minimal CMake project with the same variables does not reproduce it. The trailing separator is
emitted by CMake itself to fence off the portion of the runpath the compiler driver appends — the
install log names that portion out loud, "Set **non-toolchain portion** of runtime path", and the
`/opt/scls/mkl/lib` after the empty entry is the toolchain's. The installed binaries have no empty
entry, because the install-time rewrite replaces the CMake-controlled portion in one piece. Nothing
to fix; recorded so the next person does not re-open it.

**The plugin templates cannot build against an install.** Recorded as O9 in the plan. The installed
branch of `UserMaterialTemplate.cmake` puts `${BELFEM_DIR}/include` alone on the include path, while
R6 installs headers as `include/belfem/<module>/…`, so the flat `#include "cl_Material.hpp"` the
template itself demonstrates cannot resolve. In the same file, `BELFEM_CACHE_LOCATIONS` is set and
never read, so a plugin compiles with no backend define at all. Both are R11's, which already owns
the templates.

## State left behind

`cmake-build-debug` is configured `USE_SHARED_LIBS=ON` with `CMAKE_INSTALL_PREFIX` pointing at a
session scratch directory that will not survive. Flipping either back is a relink, not a rebuild.
Nothing was written to `/usr/local`, and `share/` and `examples/` in the source tree are as they
were — the churn probe file was removed.

The gate itself is a shell script rather than prose, kept in the session scratchpad; it inspects
and runs, never builds, so it cannot collide with a build in another worktree. The shipped
`2D_Tapestack` deck simulates 20 ms, which is hours; the gate truncates the installed copy to
0.05 ms and runs `Allclean` first, because BELFEM warm-restarts from a leftover `memdump.hdf5` and
would otherwise decide it had already finished.
