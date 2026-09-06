# Supported MPI Implementation {#doc_mpi_support}

**Date:** 2026-08-11
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

### At link time

Two places name Open MPI's libraries directly, and both are correct rather than
accidental:

| file | what it names |
|---|---|
| `config/linalg/config_mumps.cmake` | `-lmpi_mpifh`, `-lmpi_usempif08`, `-lmpi_usempi_ignore_tkr` — the Fortran interface MUMPS needs |
| `config/linalg/config_mkl.cmake` | `-lmkl_blacs_openmpi_${BELFEM_MKL_SUFFIX}` — MKL ships one BLACS per MPI |

MKL's BLACS is selected unconditionally, including on the Intel-compiler path. That is the
correct choice given the guidance above: the Intel route runs on Open MPI.

**Do not "fix portability" by adding `-lmpifort` or an MKL `blacs_intelmpi` / `blacs_mpich`
variant.** Those flags would let a MPICH build link, which advertises support that does not
exist and replaces an immediate, legible link error with a later crash inside PETSc.
Supporting MPICH means validating it, not renaming libraries.

## Checking what you have

```
mpicxx -showme:version      # Open MPI answers; MPICH-family wrappers do not recognise it
mpirun --version            # "Open MPI" or "HYDRA" / "Intel(R) MPI"
```

The configure step reports what it found, so a mismatch between the wrapper on `PATH` and
the one BELFEM used is visible in the CMake output.

## Related

- `config/system/find_mpi.cmake` — `MPI_HOME` handling and the implementation guard
- `src/comm/doc/README.md` — the communication layer built on top of this
