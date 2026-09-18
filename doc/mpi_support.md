# Supported MPI Implementation {#doc_mpi_support}

**Date:** 2026-09-18
**Purpose:** State which MPI BELFEM is built and validated against, and why the build
refuses the alternatives

## The requirement

**BELFEM requires Open MPI.** It is the only MPI implementation the framework has been
built and validated against.

**MPICH is not supported, and neither is Intel MPI**, which is MPICH-derived. This is not a
packaging oversight to be worked around — it is a statement about what has been tested:

- BELFEM has never been validated on MPICH.
- On that path, **PETSc crashed** when BELFEM called into it.

A build that reaches the link stage on MPICH is therefore not evidence that it works. The
failure that matters arrives later, inside a solver, and is harder to read than a missing
library.

### Intel compilers

Using the Intel compilers is fine, and is a route that has been used. **Build Open MPI with
the Intel compilers**; do not substitute Intel MPI for it. The MPI implementation and the
compiler are independent choices, and only the Open MPI half is fixed here.

## How the requirement is enforced

### At configure time

`config/system/find_mpi.cmake` preprocesses a small probe against `mpi.h` and reads the
implementation macros the header defines — `OPEN_MPI` for Open MPI, `MPICH` /
`MPICH_VERSION` / `I_MPI_VERSION` for the MPICH family. A MPICH-family header stops the
configure with a message naming the reason.

The check is deliberately asymmetric. An implementation it cannot identify produces a
**warning, not an error**: a probe that fails to run is not evidence of an unsupported MPI,
and an exotic-but-working setup should not be broken by a detection miss.

To configure anyway — accepting that the result is untested and that PETSc is the known
casualty:

```
cmake -DALLOW_UNTESTED_MPI=ON ...
```

This mirrors `-DALLOW_NINJA=ON`, the generator guard in the top-level `CMakeLists.txt`: a
hard default with an explicit, visible override.

`config/system/find_mpi.cmake` then links a two-call MPI program (`MPI_Init`, `MPI_Finalize`)
with the wrapper used by the build. CMake runs its own compiler test at `project()` with the
plain C++ compiler. BELFEM switches to `mpicxx` afterward. Without this probe, nothing links
against `libmpi` until the first executable, after the whole library has compiled. A failed
probe stops the configure with the linker's output. When that output names PMIx, the message
adds the package to install (`scls-<flavor>-pmix` under SCLS) and an `ldd` command that shows
which `libpmix.so.2` the wrapper's `libmpi.so` resolves to. Open MPI 5 needs the PMIx it was
built against, and an older system copy leaves its `PMIx_*` imports unresolved. This is a hard
stop, unlike the detection miss above: a wrapper that cannot link two calls cannot link `belfem`.

On Linux under SCLS, the configure also runs `ldd` on the probe when the prefix ships a PMIx.
It warns if `libpmix.so.2` resolves outside the prefix. The link uses the rpath, with the
prefix first. The loader reads `LD_LIBRARY_PATH` before `RUNPATH`, so an environment
containing another PMIx runs the program against that one. The configure summary prints the
resolved path as `PMIx`.

### At link time

Two places name Open MPI's libraries directly, and both are correct rather than
accidental:

| file | what it names |
|---|---|
| `config/linalg/config_mumps.cmake` | `-lmpi_mpifh`, `-lmpi_usempif08`, `-lmpi_usempi_ignore_tkr` — the Fortran interface MUMPS needs |
| `config/linalg/config_mkl.cmake` | `-lmkl_blacs_openmpi_${BELFEM_MKL_SUFFIX}` — MKL ships one BLACS per MPI |

MKL's BLACS is selected unconditionally, including on the Intel-compiler path. That is the
correct choice given the guidance above: the Intel route runs on Open MPI.

The rpath for every executable in both the build and install trees lists the toolchain prefix,
the compiler's own library directory when the compiler does not live under `/usr`, and the
third-party directories. It never lists the system library directories: `/usr/lib64`,
`/usr/lib`, `/lib64`, `/lib`, and the Debian multiarch directory under `/usr/lib`. Those
directories are searched by default, so adding them to the rpath adds nothing at run time. At
link time, GNU ld resolves a shared library's own dependencies through the rpath directories
first, in order. A system directory ahead of the prefix would let a distro PMIx shadow the
SCLS copy on every link. `config/compiler/detect_gcc.cmake` skips the compiler libdir when the
compiler lives under `/usr`, and `belfem_prune_system_libdirs()`
(`config/scripts/belfem_prune_system_libdirs.cmake`) removes the system directories before the
list becomes the rpath.

**Do not "fix portability" by adding `-lmpifort` or an MKL `blacs_intelmpi` / `blacs_mpich`
variant.** Those flags would let a MPICH build link, which advertises support that does not
exist and replaces an immediate, legible link error with a later crash inside PETSc.
Supporting MPICH means validating it, not renaming libraries.

## Checking what you have

```
mpicxx -showme:version      # Open MPI answers; MPICH-family wrappers do not recognise it
mpirun --version            # "Open MPI" or "HYDRA" / "Intel(R) MPI"
ldd $(mpicxx --showme:libdirs | cut -d' ' -f1)/libmpi.so | grep pmix   # show which PMIx libmpi.so resolves to
```

The configure step reports what it found, so a mismatch between the wrapper on `PATH` and
the one BELFEM used is visible in the CMake output.

## Related

- `config/system/find_mpi.cmake` — `MPI_HOME` handling, the implementation guard, the link probe and the PMIx check
- `config/scripts/belfem_prune_system_libdirs.cmake` — keeps the system library directories off the rpath
- `src/comm/doc/README.md` — the communication layer built on top of this
