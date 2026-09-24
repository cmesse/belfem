# Getting Started {#doc_getting_started}

**Date:** 2026-08-14
**Purpose:** From a fresh checkout to a running simulation — build, first run,
and where to go next. This is the user entry point; the contributor entry
point is the [Coding Philosophy](@ref doc_coding_philosophy).

---

## 1. Build

### 1.1 The dependency stack

Start with **[SCLS](https://belfem.lbl.gov/scls.html)**, the Scientific Core
Library Stack. It packages every third-party library BELFEM links by default —
Open MPI, BLAS/LAPACK, Armadillo and Blaze, STRUMPACK, MUMPS, SuperLU, PETSc,
ARPACK, HDF5, Exodus, METIS and SCOTCH — built together, with one compiler and
one integer width, and it carries a current CMake as well. It is what BELFEM is
built, tested and run against. On Enterprise Linux 9:

```bash
sudo dnf install https://belfem.lbl.gov/scls/el9/x86_64/Packages/scls-release-2026-1.el9.noarch.rpm
sudo dnf install scls-gcc
source /opt/scls/gcc/share/scls/activate
```

Packages also exist for EL10, Amazon Linux 2023 and Ubuntu 24.04. The RPM and DEB
repositories are a convenience, not a limit on where the stack runs: the same
recipes build the same stack from source into an isolated prefix on any other Linux
system, and on macOS, where Intel and Apple Silicon are both supported platforms.
Activation exports `$SCLS` and `$SCLS_FLAVOR`, and BELFEM's CMake reads both
automatically.

You *can* build against libraries you assembled yourself: the search falls back to
`/usr` and `/usr/local`, and `-D<PACKAGE>_DIR=…` points it elsewhere. We do not test
that configuration, and an inconsistent stack fails inside a solver rather than at
the link stage. [The SCLS toolchain](@ref doc_scls_toolchain) covers the flavors,
their configure-time presets, and what to watch for when building without it.

### 1.2 Configure and build

Beyond the stack, BELFEM requires CMake (≥ 4.0), a C++17 compiler and a Fortran
compiler. It also requires **Open MPI**: `USE_MPI` defaults `ON`, so a stock
configure needs one. A distribution's CMake is often older than version 4.0, so
use the one from the toolchain you build against. Open MPI is the only supported
MPI: MPICH and Intel MPI are untested, and the configure step refuses them (see
[MPI support](@ref doc_mpi_support)). The generator is pinned to Unix
Makefiles.

```bash
mkdir build && cd build
cmake ..
make -j8
```

The defaults build the Maxwell modules with STRUMPACK, MUMPS, SuperLU and
PETSc, HDF5 and Exodus I/O, and the test suite (`USE_TEST=ON`). On an ordinary
checkout this is an optimized build: `USE_DEBUG` defaults to `OFF` (`-O2` with
`NDEBUG`), so
`BELFEM_ASSERT` is compiled out. Configure with `-DUSE_DEBUG=ON` for `-Og -g`
and the full assertion set when you need it. The full option list sits at the
top of `CMakeLists.txt`.

The toolchain adds two exceptions. If `$SCLS` names the `debug` flavor,
`USE_DEBUG` presets to `ON`. An `mkl` flavor presets `USE_MKL` and `USE_PARDISO`
to `ON`. An explicit `-DUSE_DEBUG=…` always wins, and when a flavor is in play
the configure summary prints it. The full table is in
[The SCLS toolchain](@ref doc_scls_toolchain).

```bash
make check    # run the test suite
make doc      # generate this documentation (requires Doxygen)
```

## 2. First simulation

The `examples/` directory ships working decks. Most decks carry the
*geometry* (`.geo`) — generate the mesh with [gmsh](https://gmsh.info) first; three
(`corc_twolayer`, `pancake`, `rutherford`) build theirs with `python/main.py` instead, as
their READMEs say.
Starting from the repository root:

```bash
cd examples/helix
gmsh -3 helix.geo -o helix.msh
../../build/bin/belfem               # serial
mpirun -np 4 ../../build/bin/belfem  # MPI
```

The `bin/` path follows your build directory — `build/bin/` after the §1
quick start, `cmake-build-debug/bin/` in a stock IDE setup.

`belfem` is the solver, and the deck decides what it solves: an unlabeled
`linear thermal` or `nonlinear thermal` solver section requests the coupled
h-ɸ/T problem, otherwise the run is magnetic-only. It announces the choice at
startup, so there is nothing to pass on the command line. It reads
`input.conf` from the working directory and writes the field solution to an
Exodus file, the current/voltage history to `iv_results.csv`, and a restart
dump to `memdump.hdf5`. The Exodus file is named after the mesh — `helix.msh`
gives `helix.e-s`, viewable in ParaView — except on the segregated coupled
path (`coupling : segregated`), which still writes a fixed
`hphi_results.e-s`.

Three things first-time users trip over:

- **A rerun resumes.** A directory holding `memdump.hdf5` continues the
  previous run (a WARM RESTART banner names the resume point). Delete the
  dump or set `timestep { restart : false ; }` for a fresh start.
- **The first run builds material databases.** Metals with an `RRR` value
  cache a lookup table (`Copper_RRR50.hdf5`, …) in the run directory; the
  first build takes a while, later runs reuse it.
- **The mesh comes from the `.geo`, not from git.** No example ships a
  `.msh`; generate it with gmsh 4.x (fresh generation keeps the point
  elements BELFEM needs — do not convert old meshes across format versions).

The per-deck prerequisites and cleanup scripts are described in
[Examples](@ref doc_examples), which lists what each deck needs before it will
run; the helper scripts shared between decks are covered in
[Shared example scripts](@ref doc_examples_scripts).

## 3. Writing your own deck

Copy the nearest example and edit it. Every section and key of `input.conf`
is documented in the [Input File Reference](@ref doc_input_file_reference) —
syntax, units, defaults, aliases, and the pitfalls each key carries. The
machine-readable twin (`doc/input_schema.yaml`) drives `python/belfem-conf`,
which validates a deck without running the solver:

```bash
python/belfem-conf check path/to/input.conf   # from the repository root
```

`belfem-conf` needs PyYAML.

## 4. Where to go next

| If you want to | Read |
|----------------|------|
| Understand every `input.conf` key | [Input File Reference](@ref doc_input_file_reference) |
| Understand a module's internals | [Module Documentation](@ref doc_index) |
| Contribute code | [Coding Philosophy](@ref doc_coding_philosophy) — BELFEM's conventions differ deliberately from mainstream C++ |
| Trace a formulation to the literature | [Literature References](@ref doc_literature_references) |
