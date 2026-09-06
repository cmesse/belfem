# CMake Generator Pinning and Collision Guard

**Date:** 2026-07-30
**Purpose:** Diagnose the recurring "build tree corruption" that forced full cache wipes, and pin the project to the Unix Makefiles generator
**Module:** build system (`CMakeLists.txt`)

## Symptom

A terminal `make hphiTrun -j 20` in `cmake-build-debug/` re-ran the full CMake
configure (the whole TPL discovery banner: SUPERLU, MUMPS, STRUMPACK, PETSC,
METIS, SCOTCH, HDF5, EXODUS, NLOPT), reported success, and then died:

```
-- Build files have been written to: /home/christian/codes/belfem/cmake-build-debug
make[1]: CMakeFiles/Makefile2: No such file or directory
make[1]: *** No rule to make target 'CMakeFiles/Makefile2'.  Stop.
make: *** [Makefile:771: hphiTrun] Error 2
```

Recovery had always been a complete cache wipe, and the failure kept coming
back.

## Root cause

Two generators were sharing one build tree. Forensics on
`cmake-build-debug/` at the time of the failure:

| artifact | mtime | written by |
|---|---|---|
| `Makefile` | 22:14 | Unix Makefiles generator |
| `CMakeCache.txt` (`CMAKE_GENERATOR:INTERNAL=Ninja`) | 22:22 | Ninja generator |
| `build.ninja`, `CMakeFiles/rules.ninja` | 22:23 | Ninja generator |
| `CMakeFiles/Makefile2`, `CMakeFiles/Makefile.cmake` | — | absent |

CLion 2026.2 configures the Debug profile with Ninja explicitly — from
`cmake-build-debug/CMakeFiles/clion-Debug-log.txt`:

```
/opt/scls/gcc/bin/cmake -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_MAKE_PROGRAM=/home/christian/Applications/clion/bin/ninja/linux/x64/ninja \
  -G Ninja -S /home/christian/codes/belfem -B /home/christian/codes/belfem/cmake-build-debug
```

Its cache reset deletes `CMakeCache.txt` and `CMakeFiles/` but leaves the
top-level `Makefile` behind. That orphan is what a terminal `make` then picks
up, and the failure chain is:

1. `make hphiTrun` matches the stale rule at `Makefile:771`, which first
   depends on `cmake_check_build_system` (`Makefile:1170`).
2. That rule runs `cmake --check-build-system CMakeFiles/Makefile.cmake 0`.
   The file is gone, so CMake does a **full re-configure** — the TPL banner.
3. The re-configure regenerates **Ninja** files, because the cache says Ninja
   and no `-G` is passed.
4. Control returns to make, which recurses into `CMakeFiles/Makefile2` — a
   file the Ninja generator never writes. Error.

Every terminal `make` therefore churned the tree between the two generators.
Only wiping the cache cleared it, which is why the problem read as random
corruption rather than a configuration conflict.

## Changes

Both blocks sit in `CMakeLists.txt` immediately after the in-source-build
guard, before `add_subdirectory( src )`.

### 1. Generator pin (`CMakeLists.txt:188-204`)

`option( ALLOW_NINJA ... OFF )` plus a `FATAL_ERROR` for any generator not
matching `Makefiles`. The message names the CLion settings path and the shell
equivalent, and points at `-DALLOW_NINJA=ON` as the escape hatch. Placed
*before* the cleanup block below so a rejected configure aborts without
deleting anything.

### 2. Collision cleanup (`CMakeLists.txt:206-228`)

At configure time, remove the entry points belonging to the other generator:

| active generator | removed |
|---|---|
| `Ninja` | `Makefile`, `CMakeFiles/Makefile2`, `CMakeFiles/Makefile.cmake` |
| `*Makefiles` | `build.ninja`, `.ninja_deps`, `.ninja_log`, `CMakeFiles/rules.ninja` |

No false positives: a healthy Ninja tree contains no `Makefile`, and a healthy
Makefiles tree contains no `build.ninja`. Object files are untouched — both
generators use the same `<subdir>/CMakeFiles/<target>.dir/` layout, so only the
dependency bookkeeping differs.

## Verification

The cleanup block proved itself unprompted: editing `CMakeLists.txt` triggered
a CLion auto-reload, and `cmake-build-debug/CMakeFiles/clion-Debug-log.txt:47`
recorded

```
-- Removing Makefile: leftover from a different generator
```

Both branches were then exercised in a scratch tree against the real source:

- `cmake -G Ninja` → aborts at `CMakeLists.txt:197` with the pin message,
  leaving the tree untouched.
- `cmake -G "Unix Makefiles"` into a directory seeded with dummy `build.ninja`
  and `CMakeFiles/rules.ninja` → logs both removals, configures and generates
  clean, writes `Makefile` and `CMakeFiles/Makefile2`.

## Follow-up

The pin is enforced but the IDE profile still requests Ninja, so CLion will
report the fatal error on its next reload until the generator dropdown is
changed (Settings → Build, Execution, Deployment → CMake → Debug → Generator →
`Unix Makefiles`, then Reset Cache and Reload). This could not be scripted:
CLion holds `.idea/workspace.xml` in memory and rewrites it on save, and it
passes `-G` explicitly, so the `CMAKE_GENERATOR` environment variable has no
effect.

The separate `build/` tree was already a clean Unix Makefiles configuration and
was unaffected throughout.
