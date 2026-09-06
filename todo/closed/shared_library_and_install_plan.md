# Shared Libraries, `make install`, and RPATH (Linux + Darwin)

> **CLOSED 2026-09-03** (todo/ currentness sweep, round 3): the shared-library and install spine is in the tree and gated on Darwin and Linux; R11 (flip `USE_SHARED_LIBS` default) is a release-policy call, and the retired `hphirun`/`hphiTrun` names were removed from `BELFEM_INSTALL_EXECUTABLES` the same day. Status lines and checkboxes below are as they stood at closure and are not maintained.

**Date:** 2026-08-25
**Purpose:** Turn BELFEM from a build-tree-only, static-archive project into something that can be
installed and run from `$CMAKE_INSTALL_PREFIX`. Three coupled changes: build the project as a
shared library, add `install()` rules that also deploy `share/` and `examples/` in a way that
survives frequent content churn, and give both the build tree and the installed tree working
RPATHs on ELF and Mach-O.
**Module:** `CMakeLists.txt`, `config/` (+ one small change in `src/comm`, `src/physics`)
**AIs involved:** Claude (exploration + plan); Codex + Grok audit owed before implementation
**Status:** IN PROGRESS — plan drafted 2026-08-25, amended 2026-08-26 (R8b, G20, O6). **Implemented
and gated on Darwin 2026-08-27** (macOS 15.7 x86_64, GCC 16, Blaze): R1, R2, R3, R5, R6, R7, R8,
R8b and R10 done on branch `claude`; the order was inverted from the draft because the Mac was
available first. **Gated on Linux the same day** (el9.8 x86_64, GCC 11, Armadillo, MKL/SCLS, Debug):
R4 and both halves of R9 pass, so `lib64` vs `lib`, `$ORIGIN`, `DESTDIR` relocation and the
compiled-in data directory are all confirmed on ELF. **Still owed:** R11 and R12. Three findings
the plan did not anticipate are recorded as O7 (`gComm`/`gLog` defined in `main()`), O8 (data
lookup not relocatable) and O9 (the plugin templates cannot compile against an install). The Codex + Grok audit (§8) has not run; the
implementation commits are small enough to audit as a set.

> **Scope guards:**
> - **Out of scope:** Windows; `USE_VTK`/`visualizer` and `nonfree/` install rules beyond making
>   them not break the configure; a CPack/RPM/conda package; changing which TPLs are linked;
>   any change to solver or physics behaviour.
> - **Compatibility promises kept:** the static build must remain available and must remain the
>   default until R9 passes. `make`, `make check`, `make check-fast`, `make reset`, `make doc`
>   keep working unchanged. Sources and per-module `CMakeLists.txt` files are not restructured.
> - **Compatibility promise explicitly dropped (needs sign-off — see O2):** the 24 per-module
>   archives `libbelfem_<module>.a` and, with them, the `make libbelfem_<module>.a` per-module
>   build gate used throughout the devlogs and the scratchpad-probe workflow.
> - **In scope:** the 24 `src/` library targets, the 4 executables (`banner`, `hphirun`,
>   `hphiTrun`, `electricalCircuit`), the test executables, headers, `share/`, `examples/`.

---

## 1. Current Behaviour and How It Fails

The project builds 24 static archives into `${CMAKE_BINARY_DIR}/lib` and links every executable
and test against all of them inside `-Wl,--start-group` / `--end-group`
(`config/scripts/Add_Executable.cmake:56-60`, `config/scripts/Add_Test.cmake:45-52`). There is no
`install()` rule anywhere in the main tree — `grep -rn "install(" --include=CMakeLists.txt
--include=*.cmake .` returns hits only in the two plugin *templates* and in
`examples/2D_Undulator/lib/CMakeLists.txt`, never in `CMakeLists.txt`, `src/` or `config/`.

| Failure | Mechanism | Evidence |
|---|---|---|
| Cannot build shared at all | `STATIC` is hardcoded; the artifact suffix is part of the **CMake target name** | `config/scripts/Add_Library.cmake:27`, `config/globals.cmake:6` |
| Hard link error if switched | `-fPIC` reaches C++ only; Fortran and C never get it, and three `.f90` files are in `libbelfem_sparse` | `config/compiler/config_gcc.cmake:102` vs `:19,32`; `config_icc.cmake:58` vs `:21,27` |
| Undefined symbols when creating a dylib | 5 mutual dependency cycles among the module archives | measured, §1.1 |
| `.so`/`.dylib` land outside `lib/`, survive `make reset` | only `RUNTIME` and `ARCHIVE` output dirs are set | `CMakeLists.txt:262-263`, `:289-295` |
| No shared-library link interface | `target_link_libraries` for the library targets is commented out | `config/scripts/Add_Library.cmake:36` |
| `make install` does nothing | no `install()`, no `GNUInstallDirs` | whole tree |
| Installed binaries find no data files | data lookup is `$BELFEM_DATA` or CWD-relative only | `src/comm/cl_Communicator.cpp:71-76`, `fn_GT_data_path.cpp:56-70`, `fn_material_data_path.cpp:28-32` |
| Plugin authors cannot build against an install | the template probes `${BELFEM_DIR}/include`, which no rule ever creates | `src/physics/materials/UserMaterialTemplate.cmake:58-60` |
| Binaries need `LD_LIBRARY_PATH` even in the build tree | `BELFEM_RPATH` is collected but never emitted as `-Wl,-rpath` | `config/compiler/finalize_compiler.cmake:61-62,74-77` |

### 1.1 The cycles (measured, high confidence)

`nm -g` over the 24 archives in `cmake-build-debug/lib`, resolving each undefined symbol to its
defining archive, gives five mutual pairs:

```
comm <-> core        core <-> io        integration <-> mesh
interpolation <-> kernel                iwg <-> kernel
```

This is why `--start-group` is there. It also decides the architecture: ELF permits undefined
symbols in a `.so`, so naive per-module `.so`s would *appear* to work on Linux, but Mach-O `ld64`
links dylibs under a two-level namespace with `-undefined error` as the default, so a genuine
cycle cannot be built on Darwin without `-flat_namespace -undefined suppress` — deprecated, and it
breaks symbol interposition. A per-module split would therefore be Linux-only *and* would hide the
cycles until the Mac gate.

### 1.2 What the RPATH variable actually does today

Fifteen config files append to `BELFEM_RPATH` (HDF5, Exodus, TinyXML2, MKL, MPI, SuiteSparse,
MUMPS, STRUMPACK, PETSc, SuperLU, METIS, SCOTCH, NLOPT, VTK, SCLS, the compiler's own libdir).
`finalize_compiler.cmake:74-77` feeds the list into the directory `LINK_DIRECTORIES` property and
nothing else. No `-Wl,-rpath` is ever emitted; the single exception is the Intel path,
`config_icc.cmake:72`. The variable is a link-search list wearing an RPATH name.

**Bottom line:** the tree has no shared-library path, no install path, and no RPATH path, and the
five archive cycles make "one shared library per module" the wrong shape for a project that has to
link on macOS.

---

## 2. Architecture: One Aggregate `libbelfem`, Built From Per-Module OBJECT Libraries

`config/scripts/Add_Library.cmake` becomes an `OBJECT` library factory. Each module keeps its own
`CMakeLists.txt`, its own `set( LIBNAME … )`, its own `include_directories()` and its own
`add_compile_options()` (`src/visualizer/CMakeLists.txt:11` needs the last one) — none of that
changes. A single new target `belfem` then links every module's `$<TARGET_OBJECTS:…>` into one
library whose type follows `USE_SHARED_LIBS`.

Why this and not per-module shared libraries:

1. It dissolves all five cycles by construction, on both platforms, with no source refactor.
2. `-Wl,--start-group` and the reversed dependency-ordering loops in `Add_Executable.cmake:41-48`
   and `Add_Test.cmake:34-42` disappear: there is one library to link.
3. The `extern` globals in `src/core/globals.hpp:19-24` (`gTbulk`, `gRhoMin`, `gRhoMax`,
   `gBelfemDataPath`, defined only where `BELFEM_INITIALIZE_GLOBALS` is set) stay intra-library
   rather than becoming exported cross-library data references.
4. One `install(TARGETS)` entry, one SONAME, one rpath to reason about.
5. Link time drops: 27 link steps against 24 archives each become 27 link steps against one.

**Rejected:** per-module `.so`/`.dylib`. It requires breaking five cycles first — an open-ended
refactor of `core`/`io`/`comm` and of `kernel`/`iwg`/`interpolation` — for no benefit BELFEM
consumes, since nothing loads a subset of the framework. Revisit only if a genuine need for
partial loading appears.

**Rejected:** keeping static archives and installing those. It ships a 4-executable install where
each binary carries its own copy of the framework, and it leaves the dlopen'd user-material and
source-function plugins (`cl_Material.cpp:764`, `cl_Material_UserDefined.cpp:33`,
`cl_SourceFunction.cpp:224`) resolving framework symbols only through the executable's `-rdynamic`
export table — which works, but means a plugin can never be tested against anything but a full
executable.

---

## 3. Gap Table

| # | Item | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| G1 | Library type is not selectable | shared build | no — `STATIC` literal | (c) | `Add_Library.cmake:27` |
| G2 | Target name embeds `.a` | `.so` vs `.dylib` naming | no | (c) | `globals.cmake:6`; consumed at `Add_Library.cmake:27,33,39`, `Add_Executable.cmake:43,47`, `Add_Test.cmake:36,40`, `src/core/CMakeLists.txt:48-49`, `src/physics/materials/CMakeLists.txt:50` — **5 live files, small blast radius** |
| G3 | `-fPIC` missing for Fortran/C | linking `.f90` objects into a `.so` | no | (c) | `config_gcc.cmake:102` is CXX-only; `splinalg.f90`, `arpacktools.f90`, `parpacktools.f90` are in `libbelfem_sparse` |
| G4 | Five archive cycles | any dylib on Darwin | worked around with `--start-group` | (c) | §1.1, measured |
| G5 | No link interface on library targets | Darwin dylib creation; correct `DT_NEEDED` | no | (c) | `Add_Library.cmake:36` (commented out) |
| G6 | `CMAKE_LIBRARY_OUTPUT_DIRECTORY` unset | `.so`/`.dylib` in `lib/`; `make reset` cleans them | no | (c) | `CMakeLists.txt:262-263`, `:289-295` |
| G7 | No `VERSION`/`SOVERSION` | upgradeable installed library | no | (c) | `PROJECT_VERSION` exists, `CMakeLists.txt:28` |
| G8 | No `GNUInstallDirs` | `lib` vs `lib64` — **el9 uses lib64**, and the rpath must match | no | (c) | platform is `el9_8` |
| G9 | No `install(TARGETS)` | anything installed at all | no | (c) | — |
| G10 | No header install | the plugin workflow, which already assumes `include/` | no | (c) | `UserMaterialTemplate.cmake:58-60`, `UserLibraryTemplate.cmake` |
| G11 | Generated `belfem_version.hpp` not installed | headers that include it | no | (c) | `src/core/CMakeLists.txt:31` |
| G12 | `share/` not installed | gas tables + material DBs at `$prefix/share/belfem` | no | (a) — `install(DIRECTORY)` globs at install time | §6 |
| G13 | `examples/` not installed | decks at `$prefix/share/belfem/examples` | no | (a) — same | §6, O1 |
| G14 | Installed binaries cannot find `share/` | a usable install | no | (c) | `cl_Communicator.cpp:71-76`; gastables searches CWD-relative only (`fn_GT_data_path.cpp:63-67`), materials returns `""` outright (`fn_material_data_path.cpp:30-31`) |
| G15 | No RPATH emitted anywhere | running without `LD_LIBRARY_PATH` | no | (c) | §1.2 |
| G16 | No `@rpath`/`@loader_path` handling | Darwin | no | (c) | — |
| G17 | Ad-hoc codesign happens pre-install | Darwin: install rewrites load commands and voids the signature | no | (b/c) | `Add_Executable.cmake:68-77`; O3 |
| G18 | Stale `-L` baked into `CMAKE_EXE_LINKER_FLAGS` | points at the build tree | harmless but wrong for an install | (c) | `CMakeLists.txt:264-267` |
| G19 | No CMake package config | downstream `find_package(belfem)` | no | (b) | O4 — nice-to-have, not release-blocking |
| G20 | `make check` cannot find `share/` | gastables / materials tests actually running | no — they self-skip | (c) | §3.2 |

### 3.1 Cross-cutting findings

- **G3 is the first thing that will bite.** Nothing else can be tested until Fortran is compiled
  PIC. Fix it by deleting the hand-rolled `-fPIC` from both compiler configs and setting
  `CMAKE_POSITION_INDEPENDENT_CODE ON` once, which covers C, C++ and Fortran uniformly.
- **G8 and G15 are the same bug on Linux.** An rpath of `$ORIGIN/../lib` on a distro that installs
  to `lib64` produces a binary that links and then fails at exec. Derive both from
  `CMAKE_INSTALL_LIBDIR`; never spell either literally.
- **G14 makes the difference between an install that exists and one that works.** Everything else
  in §4 can be verified by `ldd`/`otool`; G14 is only caught by running an installed `hphirun`
  from a directory that is not the source tree. R9's gate must do exactly that.
- **Keep default visibility.** The `dlopen` plugin path relies on it, and `-rdynamic`
  (`config_gcc.cmake:99`) is what makes it work today. Do **not** add `-fvisibility=hidden` as
  part of this work; it is a separate decision with its own audit.
- **`USE_EXAMPLES` is unrelated to `./examples/`.** It gates small coding demos inside `src/`
  (`src/fem/kernel/CMakeLists.txt:70` and six siblings). The install of the `./examples/` decks
  must not be tied to it.

### 3.2 `make check` has no route to `share/` (added 2026-08-26, Christian)

A test run cannot assume `$BELFEM_DATA` is set in the environment, and in this tree the fallbacks
do not cover for it:

- `config/scripts/Add_Test.cmake:22` calls `add_test` with **no `WORKING_DIRECTORY`**, so ctest runs
  each test with CWD = its own build directory — e.g. `cmake-build-debug/tests/physics/gastables`.
- `fn_GT_data_path.cpp:63-67` climbs at most three levels (`../../../share/fluid`), which from there
  resolves to `<build>/share/fluid`. That does not exist, and cannot, for any out-of-source build.
- `fn_material_data_path.cpp:30-31` has no CWD fallback at all: it returns `""` unless
  `gBelfemDataPath` is set.

So the data-dependent tests reach their files only when the developer happens to have exported
`BELFEM_DATA`, and otherwise skip themselves — by design (`tests/physics/CMakeLists.txt:30-32`), but
it means a green `make check` does not establish that they ran.

**The source tree's own `share/` is the correct first guess for a test**, and CMake knows where it
is (`${CMAKE_SOURCE_DIR}/share`). See R8b. Note this is a *test-harness* fix and deliberately not a
library one: compiling the configuring machine's source path into `libbelfem` would ship
`/home/christian/codes/belfem/share` inside a distributed binary. Only `BELFEM_INSTALL_DATADIR`
(R8) belongs in the library.

*Incidental, fixed 2026-08-27:* the comment at `tests/physics/CMakeLists.txt:30` said the data
lives in `share/fluidprop`; the directory has been `share/fluid` since the 2026-08-12 rename. It was
the only stale `fluidprop` reference left in a live source or config file.

---

## 4. Ordered Steps

- [x] **R1 — PIC everywhere.** Delete `-fPIC` from `config_gcc.cmake:102` and `config_icc.cmake:58`;
      set `set( CMAKE_POSITION_INDEPENDENT_CODE ON )` in the top-level `CMakeLists.txt` before
      `add_subdirectory( src )`. *Gate:* full static rebuild still green; `readelf -d` on a Fortran
      object shows no text relocations. Independently useful, lands first, low risk.

- [x] **R2 — Decouple target name from artifact suffix.** *(after: R1)* Introduce
      `BELFEM_MODULE_TARGET( <name> )` → plain target name `belfem_<name>`, with
      `set_target_properties( … OUTPUT_NAME belfem_<name> )` kept so filenames are unchanged in the
      static build. Touch the 5 live files listed in G2. *Gate:* static build byte-identical
      artifact names; `make check` green.

- [x] **R3 — OBJECT libraries + one aggregate target.** *(after: R2; R2 folded in, 2026-08-27)* `Add_Library.cmake` emits
      `add_library( belfem_<name> OBJECT ${SOURCES} )`; a new `config/scripts/Add_BelfemLibrary.cmake`
      (included after `add_subdirectory( src )`) creates
      `add_library( belfem ${all objects} )` whose type follows `USE_SHARED_LIBS`, attaches the TPL /
      Fortran / OpenMP / MPI link interface (G5), sets `VERSION`/`SOVERSION` from `PROJECT_VERSION`
      (G7), and sets `LIBRARY_OUTPUT_DIRECTORY` alongside the existing two (G6). `Add_Executable.cmake`
      and `Add_Test.cmake` collapse to `target_link_libraries( … belfem … )`; the reversed
      ordering loops and `--start-group` go away. Add `option( USE_SHARED_LIBS … OFF )` — **default
      OFF until R9 passes**, then flipped in R11. *Gate:* static aggregate build green, `make check`
      green, before shared is even tried.

- [x] **R4 — Shared build on Linux.** *(after: R3; done 2026-08-27, evidence in
      `devlog/dl20260827_shared_libs_install_linux.md`)* Configured with `-DUSE_SHARED_LIBS=ON`.
      *Gate:* `libbelfem.so.0.9.0` (287 MB Debug) + `.so.0` + `.so` symlinks in `lib/`; `ldd` clean
      on all six executables; `make check` 12/13 — the one failure,
      `Genome.RandomizeLogScaleBranch`, was a pre-existing defect that cannot depend on the library
      type, since both types are linked from the same object files — `Genome::get_values()` did not
      clamp its decoded value to the caller's bounds the way `set_values()` and `randomize()` do,
      so the log branch overshot `max` by a few ulp; fixed the same session in
      `src/containers/cl_Genome.hpp`; `bin/banner` runs
      with `LD_LIBRARY_PATH` unset. `make reset` removal is reasoned, not executed: the `.so` is
      written into `${CMAKE_BINARY_DIR}/lib`, which the target removes wholesale, and running it
      would have cost a full rebuild of a 56 GB tree.

- [x] **R5 — RPATH block.** *(after: R4)* In the top-level `CMakeLists.txt` after
      `include( GNUInstallDirs )`: `CMAKE_SKIP_BUILD_RPATH OFF`,
      `CMAKE_BUILD_WITH_INSTALL_RPATH OFF`, `CMAKE_INSTALL_RPATH_USE_LINK_PATH ON`,
      `CMAKE_BUILD_RPATH` from `BELFEM_RPATH`, and `CMAKE_INSTALL_RPATH` =
      `<origin>/../${CMAKE_INSTALL_LIBDIR}` followed by `${BELFEM_RPATH}`, where `<origin>` is
      `$ORIGIN` on ELF and `@loader_path` on Mach-O. On Apple also `CMAKE_MACOSX_RPATH ON` and
      `CMAKE_INSTALL_NAME_DIR "@rpath"`. Set these through the CMake *variables*, never as raw
      `-Wl,-rpath` in `CMAKE_EXE_LINKER_FLAGS` — `$ORIGIN` passed as a flag gets eaten by shell
      expansion. Drop the now-redundant `-L` at `CMakeLists.txt:264-267` (G18). *Gate:*
      `readelf -d bin/hphirun | grep RUNPATH` lists the TPL directories; the binary runs with
      `LD_LIBRARY_PATH` unset.

- [x] **R6 — `install(TARGETS)` + headers.** *(after: R5)* `include( GNUInstallDirs )`;
      install `belfem` (`LIBRARY`/`ARCHIVE` → `${CMAKE_INSTALL_LIBDIR}`) and the four executables
      (`RUNTIME` → `${CMAKE_INSTALL_BINDIR}`) under `EXPORT belfemTargets`; install public headers
      to `${CMAKE_INSTALL_INCLUDEDIR}/belfem` preserving the module directory layout, plus the
      generated `belfem_version.hpp` (G11).

- [x] **R7 — `share/` and `examples/`.** *(after: R6)* The two `install(DIRECTORY …)` rules of §6.
      *Gate:* touch a new file into `share/material`, re-run `make install` **without re-running
      cmake**, confirm it appears in the prefix. That gate is the whole point of the requirement.

- [x] **R8 — Installed data lookup.** *(after: R6)* Compile `BELFEM_INSTALL_DATADIR` into the
      `comm` module (`target_compile_definitions`, value `${CMAKE_INSTALL_FULL_DATADIR}/belfem`)
      and consult it in `Communicator::set_globals()` (`cl_Communicator.cpp:71-76`) **after**
      `$BELFEM_DATA` and **before** the existing CWD-relative fallbacks, so no current behaviour
      changes. This is the plan's only `src/` change beyond CMake and needs explicit approval.

- [x] **R8b — Point `make check` at the source `share/`.** *(independent of R1-R8; can land first)*
      In `config/scripts/Add_Test.cmake`, after the `add_test` at line 22:
      `set_tests_properties( ${TESTNAME} PROPERTIES ENVIRONMENT
      "BELFEM_DATA=${CMAKE_SOURCE_DIR}/share" )`. Uses the path CMake already knows, so no test
      depends on the developer's environment. Keep it in the ctest property, **not** in a compiled
      define — see §3.2. *Gate:* with `BELFEM_DATA` unset in the shell, the gastables tests execute
      instead of skipping. Expect this to turn currently-skipped tests into running ones, so budget
      for real failures surfacing on the first run; that is the step working, not a regression.

- [x] **R9 — Linux end-to-end gate.** *(both halves green 2026-08-27)* *(after: R7, R8)* Two runs, because the original single
      gate asked for something R8 cannot deliver (found on Darwin 2026-08-27, see O8):
      **(a) relocation:** `make install DESTDIR=/tmp/stage`, move the staged tree elsewhere, run an
      installed `hphirun` on an installed example deck from a directory outside the source tree
      with `LD_LIBRARY_PATH` unset and `BELFEM_DATA` pointing at the moved `share/belfem` — proves
      the library, rpath and signature survive a move;
      **(b) data lookup:** `make install` to the configured `CMAKE_INSTALL_PREFIX` (no `DESTDIR`),
      same run with `BELFEM_DATA` **unset** — proves R8's compiled-in data dir.
      **Together these are the definition of done for Linux.**
      *Result:* both green. (b) installed to a scratch prefix, `hphirun` ran the shipped
      `2D_Tapestack` deck from `<prefix>/share/belfem/examples` with `BELFEM_DATA` **unset** and both
      `LD_LIBRARY_PATH` and the source tree out of the picture — R8 proven on ELF. (a) staged with
      `DESTDIR`, the tree moved, `readelf -d` shows `$ORIGIN/../lib64` first, `ldd` resolves
      `libbelfem.so.0` out of the moved `lib64/`, and the deck runs with `BELFEM_DATA` pointing at
      the moved `share/belfem`. The same run with `BELFEM_DATA` unset fails, exactly as O8 predicts
      on Darwin — the behaviour is identical on both platforms.

- [x] **R10 — Darwin gate.** *(after: R9 — ran first instead, 2026-08-27, evidence in
      `devlog/dl20260827_shared_libs_darwin.md)* Same sequence on the Mac: `.dylib` with
      `@rpath` install name (`otool -D`), `@loader_path/../lib` in the executables (`otool -l`),
      `make check` green, installed run green. Expect to find at least the codesign issue (O3) and
      possibly the `BELFEM_FORTRANLIBS` absolute-path block at `config_gcc.cmake:116`.

- [ ] **R11 — Flip the default.** *(after: R10 and R9, both now green)* `USE_SHARED_LIBS` default ON; update
      `CLAUDE.md` (the "Project libraries are built STATIC" line), `doc/` build documentation, and
      the two plugin templates so they link `belfem` from an install — see **O9** for the two
      concrete defects in them that the Linux gate turned up.

- [ ] **R12 — Devlog + doc sweep.** *(after: R11)* `devlog/dl2026MMDD_shared_libs_and_install.md`;
      Codex language sweep over the user-facing build documentation only.

---

## 5. Open Design Questions (not silently decided)

- [x] **O1 — Which `examples/` files get installed?** *(decided 2026-08-27: exclude the `Allclean`
      list, keep every `.msh`; installed `examples/` measured at 20 MB on the Mac.)* `examples/` is 76 MB, of which ~42 MB is
      generated `.bfm` and ~1 MB is `belfem.log`. `examples/scripts/Allclean` is the authoritative
      list of what a run produces: `*.exo`, `*.e-s.*`, `*.csv`, `*.bfm`, `memdump.hdf5`,
      `belfem.log`, `out.txt`, `err.txt`, `CircuitResults.txt`, `circuitAnalysis.out`,
      `slurm-*.out`, `slurm-*.err`. Excluding exactly that list takes the install to ~34 MB.
      A second, separable question: three of the four `.msh` files have a `.geo` beside them and are
      therefore regenerable (`inductor` 7.7 M, `costheta` 2.0 M, `2D_tapestack` 6.1 M); `corc.msh`
      (17 M) has no `.geo` and is a source file. Dropping the regenerable meshes would give ~18 MB
      but makes three decks require gmsh before they run.
      **Recommendation:** exclude the `Allclean` list, keep every `.msh`. Derive the exclude
      patterns from `Allclean` and say so in a comment, so the two lists cannot drift.
- [ ] **O2 — Losing `make libbelfem_<module>.a`.** *(implemented 2026-08-27 in R3 — the per-module
      gate is now `make belfem_<module>`; the sign-off below is still owed, the change is one
      commit to revert if refused.)* The per-module archive is a *make target*, and it
      is the standard per-module compile gate in the devlogs (`dl20260822_coulomb_gauge_stepBC.md:14`,
      `dl20260608_superlu_wrapper_audit_fixes.md:87`) and the basis of the scratchpad-probe workflow
      (probes link the prebuilt per-module `.a`s). After R3 the equivalent is `make belfem_<module>`
      — an OBJECT library is still an individually buildable target, so the compile gate survives
      under a new name — and probes link the single `libbelfem.{a,so}` instead, which is simpler.
      **Needs Christian's sign-off** because it invalidates muscle memory and every written
      instruction of that form.
- [x] **O3 — Darwin codesign ordering.** *(settled 2026-08-27: `install(CODE)` re-signs the four
      executables and the dylib after the install-time rpath rewrite; `codesign -vv` on the
      relocated install reports "satisfies its Designated Requirement". `--deep` dropped. Gated on
      x86_64 only — on Intel an invalid signature would not have killed the process, so the
      arm64 failure mode is prevented, not reproduced.)* `Add_Executable.cmake:68-77` ad-hoc signs the build-tree
      binary `POST_BUILD`. `make install` rewrites the load commands to swap build rpath for install
      rpath, which invalidates that signature; on Apple Silicon an invalid signature means the
      kernel refuses to exec. Expected fix: an `install(CODE …)` that re-runs `codesign -f -s -` on
      the installed artifacts (and dropping `--deep`, which Apple has deprecated). *Confidence
      medium — cannot be tested from Linux; settle it in R10.*
- [ ] **O4 — Ship a `belfemConfig.cmake`?** *(2026-08-27: follow-up. The `EXPORT belfemTargets`
      set exists on every installed target, so the config file is one `install(EXPORT)` away.)* Would let plugin authors write `find_package(belfem)`
      instead of the hand-rolled include-path probing in `UserMaterialTemplate.cmake:44-68`. Cheap
      once R6 has an `EXPORT` set. Release-blocking or follow-up?
- [x] **O6 — Should the test harness override an existing `$BELFEM_DATA`?** *(decided 2026-08-27:
      override, as recommended; landed with R8b.)* ctest's `ENVIRONMENT`
      property *sets* the variable for the test process; it does not defer to an inherited value,
      and `ENVIRONMENT_MODIFICATION` has no set-if-unset operation, so true "first guess" semantics
      would need a wrapper script. **Recommendation: override.** A test suite should be hermetic,
      and a stale `$BELFEM_DATA` left pointing at an older install is exactly what makes
      `make check` results irreproducible between machines. Recorded rather than assumed, because
      it does mean a developer cannot point `make check` at alternative data without editing CMake.

- [ ] **O7 — `gComm` / `gLog` are defined in `main()`, not in the library (found 2026-08-27, R4 on
      Darwin).** The first dylib link failed with `Undefined symbols: _gComm, _gLog`. Every
      executable defines both itself — `Communicator gComm; Logger gLog( <level> );` — so each
      main picks its own verbosity; the library only carries `extern` declarations
      (`cl_Communicator.hpp:222`, `cl_Logger.hpp:136`). That is ~50 definition sites across
      `src/`, all 15 test mains, and 9 mains in `nonfree/`, and it is the documented convention
      (`src/core/doc/core_usage_guide.md:189`). ELF tolerates the dangling reference inside a
      `.so`; Mach-O's two-level namespace does not.
      **Interim (landed):** `-Wl,-undefined,dynamic_lookup` on the `belfem` target, Darwin +
      shared only (`Add_BelfemLibrary.cmake`). dyld then resolves the two symbols from the
      executable at load time — the mechanism Python extension modules rely on. Cost: every
      undefined symbol in the dylib now surfaces at load instead of at link, on Darwin.
      **Recommendation:** define both in the library (`cl_Communicator.cpp`, `cl_Logger.cpp`),
      replace the per-main definitions with a call that sets the verbosity, and drop the linker
      option. Correct on both platforms and makes the plugin story honest (a plugin loaded by a
      non-BELFEM host would otherwise have no `gComm`). It is a ~50-file mechanical edit that
      reaches into `nonfree/` and changes a documented convention, so it **needs Christian's
      sign-off** and its own commit; not done on autopilot.

- [ ] **O8 — Data lookup is not relocatable (found 2026-08-27).** R8 compiles
      `CMAKE_INSTALL_FULL_DATADIR` in as an absolute path, so a staged tree that is moved finds
      its library (rpath is `@loader_path`/`$ORIGIN`-relative) but not its `share/` — the first
      Darwin gate run failed on exactly this (`File sp-ap.hdf5 does not exist`) and passed once
      `BELFEM_DATA` was set. Making data lookup relocatable means resolving the executable's own
      path at runtime (`_NSGetExecutablePath` on Darwin, `/proc/self/exe` on Linux, `dladdr` on
      both) and trying `<exe>/../share/belfem` before the compiled prefix. Worth doing — it is what
      makes a tarball install work — but it is platform-specific `src/comm` code, so it is a
      follow-up with its own audit, not part of this plan. Until then, a moved install needs
      `BELFEM_DATA`.

- [x] **O9 — The plugin templates cannot compile against an installed tree (found 2026-08-27,
      Linux; ~~open~~ **fixed and gated 2026-08-28**).** `UserMaterialTemplate.cmake:57-59` and its
      `UserLibraryTemplate.cmake` sibling
      handle the installed case by putting `${BELFEM_DIR}/include` alone on the include path. R6
      installs headers as `include/belfem/<module>/…`, preserving the module layout, so a plugin
      source that writes `#include "cl_Material.hpp"` — the form the templates' own example uses —
      finds nothing. The source-tree branch works because it lists six module directories
      explicitly; the installed branch needs the same six under `include/belfem`.

      **Fixed in both templates 2026-08-28** after the three-AI round on the material template
      (`tmp/ai_exchange/review_user_material_template.md`). Both now derive one
      `BELFEM_HEADER_ROOT` — `<dir>/src` or `<dir>/include/belfem` — and build the module list
      against it, so the two branches cannot drift apart again. Both now **probe for a header**
      (`cl_Material.hpp`, `cl_Vector.hpp`) rather than for a directory name: `EXISTS
      "${BELFEM_DIR}/include"` succeeds after any `make install`, so the old directory probe
      accepted the tree and the failure surfaced later as a missing header instead of at
      configure. *Gate, run:* headers staged into a scratch prefix exactly as `:340-346` installs
      them (626 `.hpp` under `include/belfem/<module>/`); the material template configures against
      it, reports `Header root: <prefix>/include/belfem`, and links `libmyalloy.so`.

      ~~Separately, the `BELFEM_CACHE_LOCATIONS` block below it is dead code: it is set and never
      read, so a plugin is compiled without `BELFEM_ARMADILLO`/`BELFEM_BLAZE` and silently
      disagrees with the framework about the matrix ABI.~~ **Dead code confirmed and removed
      2026-08-28 (DR-37 R13); the ABI consequence was a false positive, retracted here rather
      than deleted.** There is no matrix ABI for a plugin to disagree about: the user-material
      API is backend-free by design (`cl_Material.hpp:24-25` forbids including `cl_Vector.hpp` /
      `cl_Matrix.hpp`, both `set_user_defined_polynomial` overloads take `Cell<real>` /
      `std::vector<real>`, and `BhCurve` appears only as a forward-declared pointer), and
      `tests/physics/backendfree/` pins it. Verified by standalone build: the template plus
      `example_user_material.cpp` configures, compiles and links with no backend define, exports
      `MyAlloy_init`, resolves zero Blaze/Armadillo symbols, and needs exactly one host symbol
      (`belfem::Material::set_constant`). The block is now replaced by a comment stating the
      contract, so it is not re-added.

- [x] **O5 — Default `CMAKE_INSTALL_PREFIX`.** *(decided 2026-08-27: left at CMake's `/usr/local`.
      The Darwin gate installed under `DESTDIR` and under an explicit scratch prefix; nothing was
      written to `/usr/local`.)*

---

## 6. Installed Layout and the Two `install(DIRECTORY)` Rules

```
$CMAKE_INSTALL_PREFIX/
├── bin/                     belfem, hphirun, hphiTrun, electricalCircuit,
│                            material, gas, db2exo
│                            (banner dropped 2026-08-31, superseded by
│                            `belfem --version`)
├── lib64/                   libbelfem.so.0.9.0 → .so.0 → .so     (lib/ on non-multilib)
├── include/belfem/          public headers, module layout preserved, + belfem_version.hpp
└── share/belfem/
    ├── fluid/               from ./share/fluid
    ├── material/            from ./share/material
    └── examples/            from ./examples
```

```cmake
install( DIRECTORY ${CMAKE_SOURCE_DIR}/share/
         DESTINATION ${CMAKE_INSTALL_DATADIR}/belfem
         USE_SOURCE_PERMISSIONS )

install( DIRECTORY ${CMAKE_SOURCE_DIR}/examples/
         DESTINATION ${CMAKE_INSTALL_DATADIR}/belfem/examples
         USE_SOURCE_PERMISSIONS
         # run output, per examples/scripts/Allclean — keep the two lists in step
         PATTERN "*.exo"              EXCLUDE
         PATTERN "*.bfm"              EXCLUDE
         PATTERN "*.csv"              EXCLUDE
         PATTERN "memdump.hdf5"       EXCLUDE
         PATTERN "belfem.log"         EXCLUDE
         PATTERN "out.txt"            EXCLUDE
         PATTERN "err.txt"            EXCLUDE
         PATTERN "CircuitResults.txt" EXCLUDE
         PATTERN "slurm-*"            EXCLUDE )
```

Three details that carry the requirement:

1. **This is already the dynamic form.** `install(DIRECTORY)` walks the tree at *install* time, so
   a file added to `share/` or `examples/` after the tree was configured is picked up by the next
   `make install` with no re-configure. The anti-pattern to avoid is `file(GLOB)` + `install(FILES)`,
   which freezes the file list at configure time.
2. **The trailing slash on the source path is load-bearing.** Without it the rule installs into
   `share/belfem/share/` and `share/belfem/examples/examples/`.
3. **`USE_SOURCE_PERMISSIONS` is not optional here.** `examples/*/Allrun`, `Allclean` and
   `examples/scripts/*` are executable (`-rwxr-xr-x`), and the `DIRECTORY` form's default
   permissions differ from the `FILES` form's. Spell it out rather than relying on the default.
4. `install(DIRECTORY)` adds and overwrites but never removes. A file deleted from `share/` stays
   in a prefix that was installed into before. Acceptable; worth one sentence in the docs.

**DESTDIR, for the record:** `DESTDIR` is a staging prefix that is stripped when the package is
unpacked on the target, so RPATH must encode `CMAKE_INSTALL_PREFIX` and never `DESTDIR`. The
`$ORIGIN` / `@loader_path`-relative first entry in R5 is what makes the tree relocatable, so
`make install DESTDIR=/staging` followed by moving the tree resolves correctly — which R9 tests
explicitly rather than assuming.

---

## 7. Definition-of-Done Checklist

- [ ] Every gap-table row (G1–G19) maps to a step or an open question.
- [ ] Each claimed gap carries a citation; the cycle claim is measured, not assumed.
- [ ] `USE_SHARED_LIBS=OFF` still produces the current static build and green `make check`.
- [x] `USE_SHARED_LIBS=ON` on Linux: `ldd` clean, no `LD_LIBRARY_PATH`, `make check` 12/13 with the
      one failure pre-existing and unrelated to linkage ( 2026-08-27 ).
- [ ] With `BELFEM_DATA` unset in the shell, `make check` runs the data-dependent tests rather than
      skipping them (R8b).
- [x] Adding a file to `share/` and re-running `make install` deploys it with no re-configure (R7)
      — gated on both platforms.
- [~] A relocated, staged install runs an example deck outside the source tree with
      `LD_LIBRARY_PATH` unset (R9(a), green) — but it still needs `BELFEM_DATA`, because the
      compiled-in data directory is absolute. That half is O8 and is deliberately out of this plan.
      The non-relocated install needs neither variable (R9(b), green).
- [x] The same on Darwin, with `otool -D`/`-l` evidence and a valid signature after install (R10).
- [ ] O1–O5 resolved in place with dates, not silently decided.
- [ ] Codex + Grok audit round on this plan before R1 starts.

---

## 8. Audit Trail

- Exchange thread: `tmp/ai_exchange/shared_libs_install.md` — **not yet opened.** Per
  `feedback_three_vendor_audits`, both Codex and Grok audit this plan before implementation, and
  each of R1–R11 goes through plan+audit → code+audit.
- Evidence produced while drafting (Claude, 2026-08-25): the archive dependency graph was computed
  with `nm -g` over `cmake-build-debug/lib/*.a` and the symbol→owner resolution done in a scratch
  script; the `examples/` size and mesh/`.geo` classification were measured with `du`/`find`. Both
  are reproducible from the tree. Everything else in §1 and §3 is *reviewed* — read from the CMake
  files, not executed.
