# DR-35 and DR-36 — five seeded build claims checked against the tree

**Date:** 2026-08-11
**Purpose:** Resolve two `[seeded — confirm]` register rows by inspection
**Modules:** `sparse`, build configuration

Both rows came from the 2026-08-05 seeding pass and had never been checked. Between them
they carried five sub-claims, and every one turned out to be mechanically decidable — no
code written, nothing to run beyond linking a binary.

**Four of the five are now settled. One survives, and it is bigger than the row said.**

## DR-35 — STRUMPACK

### Version: 8.0.0, everywhere

The row asked to confirm ≥ 7.1.2 for the small-problem hang fix. All three SCLS toolchains
carry **8.0.0**:

```
/opt/scls/mkl/lib64/libstrumpack.so.8.0.0
/opt/scls/gcc/lib64/libstrumpack.so.8.0.0
/opt/scls/debug/lib64/libstrumpack.so.8.0.0
```

corroborated independently by PETSc's own record of what it was built against
(`petscpkg_version.h:301-303`: MAJOR 8, MINOR 0, SUBMINOR 0), and `ldd` on a binary linked
this session resolves `libstrumpack.so.8.0`. Confirmed.

### The Krylov hardening is not needed — the return-code path already carries it

The row wanted a warning when BLR runs out of Krylov iterations. STRUMPACK 8.0.0 has a
dedicated return code for exactly that:

```cpp
NO_CONVERGENCE,     /*!< The iterative solver did not converge. */   // StrumpackParameters.hpp:55
```

BELFEM already decodes it (`strumpacktools.cpp:27`), and the wrapper treats **every**
non-SUCCESS code through the soft-fail contract or a `BELFEM_ERROR` quoting the decoded
message — on the serial path (`cl_SolverSTRUMPACK.cpp:262-283`) and the distributed one
(`:336-373`) alike, the latter reducing the failure flag across ranks precisely so a
rank-local SUCCESS cannot mask a peer's failure. A separate `Krylov_iterations() == maxit`
probe would be redundant with that.

**Confidence, split, because the two halves are not equally supported:**

- The *handling* is **high** — it is right there at the cited lines.
- Whether STRUMPACK *emits* `NO_CONVERGENCE` at maxit rather than SUCCESS is **medium**.
  The source is not installed (headers and `.so` only), the header's `solve()`
  documentation does not state it, and the only positive evidence is the string
  `"maximum Krylov iterations"` in the binary. Settling it needs the STRUMPACK source or a
  contrived non-converging run.

Exposure is bounded either way: BLR compression is OFF by default.

## DR-36 — three build claims, two dead

### "dynamic MKL link untested" — refuted; it is the daily driver

`config_mkl.cmake:38-39` links dynamically by construction — *"link dynamically;
finalize_compiler.cmake turns this into -L and rpath"* — and puts the libdir on
`BELFEM_RPATH`. `readelf -d` on a binary linked this session lists all five MKL libraries
as `NEEDED` shared objects:

```
libmkl_scalapack_lp64.so.2   libmkl_intel_lp64.so.3   libmkl_gnu_thread.so.3
libmkl_core.so.3             libmkl_blacs_openmpi_lp64.so.2
```

resolving into `/opt/intel/oneapi/redist/lib`. The binary then ran. Not untested — it is
the configuration in daily use.

### "PETSc `-lX11`/`-lmpfr` removal unconfirmed" — confirmed removed

The call is commented out at `config_petsc.cmake:33` with its reason recorded ("only needed
if PETSc was built with X support"), and the full `ldd` closure of a PETSc-linked binary
(`libpetsc.so.3.25` present) contains no libX11 and no libmpfr.

### "MPICH needs `-lmpifort`" — the row's framing was wrong, and so was mine

I first wrote this up as a portability gap: two hardcoded OpenMPI sites, fix by adding the
MPICH variants. **Christian corrected the premise:**

> We have not validated BELFEM for MPICH actually. So far, PETSc always crashed when
> calling it. I always compiled openmpi with intel when I used it.

That changes what the item *is*. OpenMPI is the only validated MPI. MPICH — and Intel MPI,
which is MPICH-derived — is an unvalidated configuration carrying a **known PETSc crash**,
and even the Intel-compiler route goes through a self-built OpenMPI rather than Intel MPI.

So the two hardcoded sites —

| site | what is hardcoded |
|---|---|
| `config_mumps.cmake:26` | `-lmpi_mpifh -lmpi_usempif08 -lmpi_usempi_ignore_tkr` — "assume openmpi (mpich needs -lmpifort)" |
| `config_mkl.cmake:60-62` | `-lmkl_blacs_openmpi_${BELFEM_MKL_SUFFIX}` — "BLACS, assume openmpi" |

— are **consistent with the actual support posture, not defects.** They also corroborate
the correction: `config_mkl.cmake` picks the OpenMPI BLACS *unconditionally*, including on
the Intel-compiler path, which is exactly right for someone who builds OpenMPI with Intel
rather than using Intel MPI.

And it retires the fix I proposed. Adding `-lmpifort` and an MKL `blacs_intelmpi` variant
would advertise support that has never existed, and would move a MPICH user from an
unexplained link error to a crash inside PETSc — strictly worse, because the failure would
arrive later and further from its cause.

### The gap that is real: the requirement is written down nowhere

Checking the correction turned up something neither the row nor I had looked for.
`README.md`, everything under `doc/`, and `CLAUDE.md` contain **no occurrence** of
"OpenMPI", "MPICH", "Intel MPI" or "impi". And `config/system/find_mpi.cmake` inspects only
`MPI_HOME` — whether it is set, and whether it sits inside SCLS — and never which
implementation it points at.

So BELFEM has a hard requirement on OpenMPI whose only trace in the repository is two
comments buried in the linalg configs. A new open-source user on MPICH gets an unexplained
link failure, or links and then crashes inside PETSc, with nothing anywhere to tell them
why.

**Both done the same session, on Christian's go-ahead.** Neither is "add MPICH support".

### 1. `doc/mpi_support.md`

States the requirement, names the PETSc crash as the reason rather than leaving it as
folklore, gives the Intel-compiler guidance (build Open MPI with them; do not substitute
Intel MPI), and records the two hardcoded link sites as **correct** — with an explicit
instruction not to add `-lmpifort` or an MKL `blacs_intelmpi` / `blacs_mpich` variant, so
the fix retired above is not re-proposed by the next reader who finds those flags and sees
a portability bug. Indexed in `doc/README.md`, with pointers from `src/comm/doc/README.md`
and `CLAUDE.md`'s Build Configuration section.

Cited by file and flag string rather than by line number, per the lesson recorded in
`CLAUDE.md`'s input-contract section.

### 2. The configure-time guard

In `config/system/find_mpi.cmake`. It preprocesses a small probe against `mpi.h` and reads
the implementation macros **the header itself defines** — `OPEN_MPI` against
`MPICH` / `MPICH_VERSION` / `I_MPI_VERSION` — so nothing depends on wrapper flags, compiler
names or `mpirun` output. `CMAKE_CXX_COMPILER` is already `mpicxx` at that point
(`detect_compiler.cmake` runs first), so `mpi.h` is on the default include path.

A MPICH-family header aborts the configure with a message naming PETSc and the override;
`-DALLOW_UNTESTED_MPI=ON` proceeds with a warning. That mirrors `-DALLOW_NINJA=ON` in the
top-level `CMakeLists.txt` — a hard default with an explicit, visible escape hatch — rather
than inventing a new convention.

**The asymmetry is deliberate and is the part worth defending: an implementation the probe
cannot identify produces a warning, not an error.** A probe that fails to run is not
evidence of an unsupported MPI. Erring hard there would convert a detection miss into a
broken build for someone whose setup works, which is a worse failure than the one the guard
exists to prevent.

### Verification

All five paths exercised through a standalone `cmake -P` harness that includes the real
`find_mpi.cmake` — no BELFEM configure, so the shared build tree was never touched:

| input | outcome |
|---|---|
| real `mpicxx` | `Open MPI`, silent |
| `g++` (no `mpi.h`, rc ≠ 0) | warning, proceeds |
| synthetic MPICH | `FATAL_ERROR` with the full message |
| synthetic MPICH + `ALLOW_UNTESTED_MPI=ON` | warning, proceeds |
| synthetic unrecognized | warning, proceeds |

The probe was also run directly against the real toolchain: `mpicxx -E` emits
`BELFEM_MPI_FLAVOR OPENMPI`, and `g++` fails with rc = 1 rather than producing a false
MPICH verdict. `scripts/check_doc_claims.py` reports 32/32 after the `CLAUDE.md` edit.

## Method note

Nothing here needed a build or a run of BELFEM proper. Four of five claims fell to reading
the CMake config plus `readelf`/`ldd` on a binary linked against the tree's own libraries,
and the fifth to reading two installed headers. This is the same **code → claim** direction
`scripts/check_doc_claims.py` was built for, applied by hand to a register row rather than
to a convention document — and it found the same class of defect: statements that were true
when written and quietly stopped being true.

## Left open

- DR-35's medium-confidence half: whether STRUMPACK emits `NO_CONVERGENCE` at maxit. Not
  worth chasing unless BLR is turned on.
- MPICH support itself — out of scope **by decision**, not pending test. The guard makes
  that decision visible instead of leaving it to a link error.

Both register rows are closed.

## Postscript on method

Four of the five claims fell to inspection, which is the happy case. The fifth is the
instructive one: it was *checkable* — I confirmed both hardcoded sites and the `NEEDED`
entries in the binary, and every fact I reported about them was correct — and I still drew
the wrong conclusion from it, because the missing premise was not in the tree. "MPICH is
unvalidated and PETSc crashes there" is not visible from CMake files, link lines or `ldd`
output. A mechanical claim → tree check catches stale facts; it cannot supply a support
posture that was never written down. That is precisely the gap the two proposed steps
close.
