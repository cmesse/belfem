# The SCLS Toolchain {#doc_scls_toolchain}

**Date:** 2026-09-22
**Purpose:** Name the dependency stack BELFEM is built and tested against, how to
install and activate it, and what its presence changes in the configure step

## The short version

**BELFEM is developed against [SCLS](https://belfem.lbl.gov/scls.html), the Scientific
Core Library Stack.** SCLS packages every third-party library BELFEM links by default —
Open MPI, BLAS/LAPACK, Armadillo and Blaze, STRUMPACK, MUMPS, SuperLU, PETSc, ARPACK,
HDF5, Exodus, METIS and SCOTCH — as one coherent set, built with one compiler and one
integer width. It carries a current CMake as well. SCLS is developed alongside BELFEM at
Lawrence Berkeley National Laboratory, and a machine with SCLS installed is the environment
in which BELFEM is built, tested and run.

```bash
sudo dnf install https://belfem.lbl.gov/scls/el9/x86_64/Packages/scls-release-2026-1.el9.noarch.rpm
sudo dnf install scls-gcc
source /opt/scls/gcc/share/scls/activate
```

Then configure BELFEM as usual. The activation script exports `$SCLS` and `$SCLS_FLAVOR`,
and BELFEM's CMake reads both without further arguments.

## Is it required?

Not by the build system, and that distinction is worth being precise about.

`config/system/find_scls.cmake` runs only `if( DEFINED ENV{SCLS} )`. The package finder
(`config/scripts/belfem_find_package.cmake`) searches in this order: an explicit
`-D<PACKAGE>_DIR`, an environment `<PACKAGE>_DIR`, the SCLS root, any extra roots the caller
passes, and finally `/usr` and `/usr/local`. A build that finds every library under `/usr`
configures and links.

What it does not get you is the thing SCLS exists to provide: **a stack whose parts agree
with each other.** The libraries BELFEM links are not independent. PETSc, STRUMPACK and
MUMPS are each built against a particular MPI and a particular ScaLAPACK. BLAS, LAPACK and
ScaLAPACK share an integer width. The matrix backend, the solvers and BELFEM itself must use
the same compiler and the same Fortran name mangling. Distribution packages are built
independently, so they are not guaranteed to line up. When they do not, the failure rarely
arrives at the link stage — it arrives later, inside a factorization, and is hard to read.

So: assembling the dependencies by hand is possible and is nobody's mistake. It is simply
not the configuration BELFEM is tested in, and a problem reported from one is harder for us
to reproduce. **If you are starting out, install SCLS.**

Two things sharpen the point:

- **Open MPI, and no other MPI.** `USE_MPI` defaults `ON`, so a default build needs one, and
  the configure refuses MPICH and Intel MPI (see [MPI support](@ref doc_mpi_support)). Every
  SCLS flavor ships Open MPI, so activating one settles the question.
- **CMake ≥ 4.0**, which is newer than several distributions ship. SCLS carries its own
  CMake, and activation puts it on `PATH`.

## Installing and activating

The canonical instructions, with the current release numbers, live at
<https://belfem.lbl.gov/scls.html>. In outline:

| Platform | How |
|---|---|
| Enterprise Linux 9 / 10, Amazon Linux 2023 | install the `scls-release` RPM, then `dnf install scls-<flavor>` |
| Ubuntu 24.04 LTS | install the `scls-archive-keyring` package, then `apt install scls-<flavor>` |
| Any other Linux | `git clone https://github.com/cmesse/scls.git`, write a `flavor.conf`, then run `./scls build all` |
| macOS (Intel and Apple Silicon) | the same source build; SCLS's `doc/MACOS_BUILD.md` covers the prerequisites |

**The binary repositories are a convenience, not the portability boundary.** They exist where
they add something — installation and removal through the system package manager. The same
recipes and flavors build the same stack from source into an isolated prefix everywhere else:
an enterprise Linux release older than the repositories cover, a cluster login node you do not
administer, a distribution on neither packaging path, and macOS. A source build requires what
any bootstrap does — a working system toolchain, the system prerequisites each recipe
declares, disk, and hours. The `macos` flavor is the largest because it builds its own GCC to
get Fortran.

Native Apple Silicon (arm64) is a supported platform, not an experiment: both the stack and
BELFEM build there.

Each flavor has one activation script. To switch flavors, source the new one: the paths are
replaced rather than accumulated.

```bash
source <prefix>/share/scls/activate
```

`<prefix>` is the `prefix:` value in the SCLS flavor file `flavors/<flavor>.yaml`. For the
packaged flavors (`gcc`, `mkl`, `debug`), use `/opt/scls/<flavor>`. A source-built stack uses
the prefix set for its flavor. For example, `macos` uses `/opt/scls`, so its activation script is
`/opt/scls/share/scls/activate`.

It sets `PATH`, `PKG_CONFIG_PATH`, `CMAKE_PREFIX_PATH`, `SCLS` and `SCLS_FLAVOR`, and, where
those are installed, `PETSC_DIR` / `SLEPC_DIR` / `MPI_HOME`. It deliberately does **not** touch
`LD_LIBRARY_PATH`: runtime lookup goes through rpaths, and BELFEM's own link line follows the
same policy.

`scripts/scls_env.sh` in this repository is a different thing and is easy to confuse with the
above. It is not an SCLS script: it is a small Linux-only pin that BELFEM's own hooks and
helper scripts source so they cannot run against a stray toolchain. It prepends
`/opt/scls/gcc/bin` and a pinned Python prefix to `PATH`, then checks that `cmake` and
`mpicxx` resolve under `/opt/scls` — warning in an interactive shell, failing in a
non-interactive one. A user building BELFEM does not need it; `activate` is the supported
entry point.

### The three public flavors

| Flavor | Stack | Use it for |
|---|---|---|
| `gcc` | GCC + OpenBLAS, LP64 integers, x86-64-v3 | the default production build |
| `mkl` | GCC + Intel oneAPI MKL for BLAS/LAPACK/ScaLAPACK | sites standardized on MKL; also gives BELFEM PARDISO |
| `debug` | GCC + Netlib reference BLAS/LAPACK, `-Og -g` | Valgrind, sanitizers, a debugger |

The source build additionally supports flavors such as `intel`, `gcc-mkl-cuda` and `macos`.
`mkl` is a Linux-only flavor: BELFEM's `config/linalg/config_mkl.cmake` stops the configure
when `USE_MKL` is set on an Apple host. The `macos` flavor supplies OpenBLAS instead, and
since the preset above keys on the flavor *name*, `USE_MKL` stays `OFF` there without being
asked to.

## What activation changes in the BELFEM configure

Four separate mechanisms read `$SCLS`, and it helps to know they are separate.

**1. Build defaults, from the flavor name.** `config/system/find_scls_flavor.cmake` runs
before the `option()` calls and presets three switches:

| `$SCLS_FLAVOR` | `USE_MKL` | `USE_PARDISO` | `USE_DEBUG` |
|---|---|---|---|
| name contains `mkl` | `ON` | `ON` | `OFF` |
| `debug` | `OFF` | `OFF` | `ON` |
| anything else (`gcc`, …) | `OFF` | `OFF` | `OFF` |
| `$SCLS` unset | `OFF` | `OFF` | `OFF` |

These are *defaults only*. An explicit `-DUSE_DEBUG=…` or a `ccmake` edit always wins. A value
already in a build tree's cache is never overwritten, so repointing `$SCLS` does not change that
tree's knobs. The configure summary prints the flavor it read, which is how a mismatch becomes
visible. BELFEM checks nothing against the host: `$SCLS_FLAVOR` is your statement about which
toolchain you sourced.

**2. Library, include and rpath directories.** `find_scls.cmake` adds `$SCLS/lib64` (or
`$SCLS/lib`) to the link directories and to BELFEM's rpath, and `$SCLS/include` to the include
path. If `$SCLS` is set but points at a tree with neither, the configure aborts — a half-sourced
environment fails loudly instead of silently falling back to `/usr`.

**3. Package discovery.** `$SCLS` enters the search order of `belfem_find_package` as
described above, ahead of the system prefixes.

**4. Compilers, only if the flavor ships them.** If neither `-DCMAKE_CXX_COMPILER` nor `$CXX`
is set and `$SCLS/bin/g++` exists, BELFEM pins the `gcc`/`g++`/`gfortran` triple from that
directory. The block exists because CMake's default C-compiler search looks for `cc` before
`gcc`; on macOS that would pair AppleClang with SCLS's `gfortran`. The Linux binary flavors
ship no compiler binary — they are library stacks built against the system GCC — so on Linux
the block does nothing, and the system compiler is used, which is the intended policy there.

## Building without SCLS

If you are assembling the dependencies yourself, the list is the set of default-`ON` options at
the top of `CMakeLists.txt`. Point CMake at anything outside `/usr` with `-D<PACKAGE>_DIR=…`,
one per library, and keep one rule: **one compiler, one MPI, one BLAS integer width across all
of them.**

Two BELFEM options have no SCLS package behind them — `USE_SUITESPARSE` and `USE_TINYXML2`.
Both default `OFF`, so neither is in the way of a default build.

## Where this is stated elsewhere

`README.md` and [Getting Started](@ref doc_getting_started) name SCLS at the point where a
reader is about to run `cmake`. This page is the long form; the SCLS site itself is the
authority on release numbers, platform coverage and package contents, and is where to look
when the commands above have moved on.

This page does not document SCLS's own build system, its package recipes, or how to add a
package to a source build. Those belong to the SCLS repository.
