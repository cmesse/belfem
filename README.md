# BELFEM

**Berkeley Lab Finite Element Framework** — a C++17 finite element framework
built for high-performance computing, with a focus on the electromagnetics of
high-temperature superconductors (HTS).

BELFEM solves coupled magneto-thermal problems on distributed meshes. Its
Maxwell module implements the h-φ formulation with thin-shell elements and
automatically generated cohomology cuts; the surrounding infrastructure
provides MPI-parallel mesh handling, backend-agnostic linear algebra, and
interfaces to the major sparse direct solvers (STRUMPACK, MUMPS, PETSc,
SuperLU, and others). Lumped electrical circuits can be coupled to the field
problem, and thermal quench simulation is built in.

## Building

BELFEM builds with CMake and requires a C++17 compiler, a Fortran compiler,
and — for parallel runs — **Open MPI** (the only supported MPI; MPICH and
Intel MPI are refused at configure time). The generator is pinned to Unix
Makefiles.

```bash
mkdir build && cd build
cmake ..            # defaults: MPI, Maxwell modules, STRUMPACK/MUMPS/PETSc,
                    # HDF5 + Exodus I/O, tests ON, optimized flags
make -j8
make check          # run the test suite
make doc            # generate the documentation (requires Doxygen)
```

Unless `$SCLS` names the `debug` flavor, this is already an optimized build: `USE_DEBUG`
defaults to `OFF` (`-O2 -DNDEBUG`), so the assertions are compiled out. Configure with `-DUSE_DEBUG=ON` for `-Og -g` and the full
assertion set. The main configuration options are listed at the top of `CMakeLists.txt`.

## First simulation

The `examples/` directory ships working input decks. Every example carries the
geometry (`.geo`) — generate the mesh with [gmsh](https://gmsh.info) first.
From the repository root:

```bash
cd examples/helix
gmsh -3 helix.geo -o helix.msh
../../build/bin/belfem         # use your build directory's bin/
```

`examples/README.md` explains the workflow, per-deck prerequisites, and the
material-database caches. The complete `input.conf` reference lives in
`doc/input_file_reference.md`.

## Documentation

`make doc` renders the full documentation — module guides, theory notes, and
the API reference — into `<build>/doc/html/index.html`. The hand-written
sources live in `doc/` and `src/<module>/doc/`; start with
`doc/coding_philosophy.md` before contributing code, as BELFEM deliberately
departs from mainstream modern-C++ conventions in favor of HPC performance.

## Citing BELFEM

If you use BELFEM in academic work, please cite:

> C. Messe et al., *"BELFEM: a special purpose finite element code for the
> magnetodynamic modeling of high-temperature superconductor tapes"*,
> Superconductor Science and Technology, 2023.
> DOI: [10.1088/1361-6668/acf7f9](https://doi.org/10.1088/1361-6668/acf7f9)

Machine-readable citation metadata is in `CITATION.cff`.

## License

See the `LICENSE` file. BELFEM was developed at Lawrence Berkeley National
Laboratory under U.S. Department of Energy funding; the notice below governs
its use and redistribution.

---

Berkeley Lab Finite Element Framework (BELFEM) Copyright (c) 2026,
The Regents of the University of California,
through Lawrence Berkeley National Laboratory
(subject to receipt of any required approvals from the U.S. Dept. of Energy).
All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative 
works, and perform publicly and display publicly, and to permit others to do so.
