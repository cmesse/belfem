# Devlog 2026-07-01 - CMake belfem_find_package Migration

**Date:** 2026-07-01
**Topic:** Tidy top-level CMake includes from L100 onward; migrate library configs to `belfem_find_package()`
**AIs involved:** Claude Code (interactive with Christian)
**Literature References:** N/A

## Summary

Walked the `CMakeLists.txt` include chain from line 100 down, call by call, fixing
latent bugs and migrating the hand-rolled `XXX_DIR → SCLS → /usr` library searches
to the new `belfem_find_package()` helper (introduced with the nlopt integration).
Every migrated finder was verified with a configure-only run against the SCLS
install (`/opt/scls/gcc`); all resolve shared libraries, matching the
dynamic-first linking policy.

## Design Decisions (Christian's calls)

- **Dynamic linking everywhere** unless only a static library exists;
  `find_library`'s default `.so`-before-`.a` order implements this.
- **Link order is solvers → linalg → MPI.** The previously empty
  `BELFEM_SOLVER_LIBS` (linked before `BELFEM_MATRIX_LIBS` in
  `Add_Executable.cmake`) now holds MUMPS/STRUMPACK/PETSc; graph libraries
  (METIS/SCOTCH) stay in the matrix list so they link after the solvers that
  need them; MPI Fortran libs are appended at the tail.
- **No SCLS layout assumptions:** SCLS is the preferred root, not the assumed
  layout. Explicit `-lsuperlu -larpack` stay on the Armadillo line because only
  the SCLS Armadillo build is guaranteed to carry those deps itself.
- **MKL:** dropped macOS support (fail fast, MKL is no longer supported there);
  switched from static `libmkl_*.a` paths to dynamic `-lmkl_*` flags with the
  libdir registered in `BELFEM_RPATH`. Untested — no MKL on this machine.
- **PETSc** (used for thermal problems): supports both prefix installs (via
  helper) and classic `PETSC_ARCH` in-tree builds; `-lX11`/`-lmpfr` commented
  out pending a compile test.

## Changes Made

- `config/system/find_unix.cmake` — reduced to the Unix guard; fixed invalid
  `error()` call (was never a CMake command); removed APPLE deployment-target
  block duplicated from top-level pre-`project()` code.
- `config/system/find_mpi.cmake` — removed broken guard comparing the MPI root
  against the SCLS *lib dir* (hard CMake error when `SCLS` unset); rpath append
  now unconditional (dedup in finalize); warns when `MPI_HOME` points outside
  SCLS.
- `config/system/find_mkl.cmake` — fixed `message( FATAL_ERROR, ...)` stray
  comma that demoted the abort to a printout.
- `config/linalg/config_mkl.cmake` — dynamic MKL (see above); lp64/ilp64
  handled via `BELFEM_MKL_SUFFIX`; oneAPI 2024+ `lib` layout fallback; dead
  commented blocks removed. Netlib branch untouched.
- `config/linalg/config_matrix.cmake` — backend checks (both/neither) moved to
  top with fail-fast; Armadillo link line unchanged.
- `config/linalg/config_suitesparse.cmake` — migrated to helper
  (`umfpack.h`, 9 libs, `INCLUDE_SUFFIXES suitesparse .`); GPL warning kept,
  "partlylicensed" typo fixed.
- `config/linalg/config_superlu.cmake` — migrated (`slu_ddefs.h`, `superlu`);
  `BELFEM_METIS` guard and `ARMA_DONT_USE_SUPERLU` logic unchanged.
- `config/linalg/config_mumps.cmake` — migrated (`dmumps_c.h`; dmumps,
  mumps_common, esmumps, pord) into `BELFEM_SOLVER_LIBS`; identical
  Intel/GNU MPI-Fortran branches collapsed (openmpi assumed, mpich hint kept).
  Note: the MPI Fortran libs sit before METIS/SCOTCH in the matrix list (their
  configs are included later), not at the absolute tail — harmless with shared
  libs.
- `config/linalg/config_strumpack.cmake` — migrated
  (`StrumpackSparseSolver.hpp`; strumpack, 4× butterflypack, slate_lapack_api,
  slate, + zfp when `USE_STRUMPACK_ZFP`) into `BELFEM_SOLVER_LIBS`; CUDA flag
  block kept; redundant `BELFEM_FCFLAGS` include append dropped (finalize adds
  `BELFEM_INCLUDES` to Fortran flags already).
- `config/linalg/config_petsc.cmake` — rewritten for both layouts (see above);
  "Turn of MPI" typo fixed; errors if `PETSC_ARCH` set without `PETSC_DIR`.
- `config/linalg/config_metis.cmake`, `config_scotch.cmake` — migrated
  (metis.h/parmetis.h; scotch.h/ptscotch.h); def dedup guards kept because
  finalize does NOT dedupe `BELFEM_DEFS` and config_superlu also sets
  `BELFEM_METIS`.
- `config/io/config_hdf5.cmake` — migrated (`hdf5.h`/`hdf5_hl.h`; hdf5_hl,
  hdf5; `-lz -ldl` stay as plain flags).
- `config/io/config_exodus.cmake` — migrated (`exodusII.h`; exodus, netcdf).
- `config/io/config_tinyxml2.cmake` — NEW; absorbs the inline tinyxml2 snippet
  from `CMakeLists.txt`.
- `config/scripts/Add_Test.cmake` — links `BELFEM_SOLVER_LIBS` (tests
  previously received solver libs via the matrix list).
- `CMakeLists.txt` — `if( Apple )` → `if( APPLE )` (macOS silently got the
  Armadillo default); duplicate `option( USE_EXODUS )` removed; redundant
  second `find_package(VTK REQUIRED)` removed (find_vtk.cmake already calls it
  with hints); IO section reduced to three symmetric includes.

## License Fact-Check: SCOTCH

The `USE_SCOTCH` option said "(non-commercial use only)" — outdated since
SCOTCH 4.0 (Feb 2006). Current SCOTCH (SCLS ships 7.0.11) is CeCILL-C,
LGPL-like: linking into commercial/proprietary software is permitted;
obligations apply only to modifications of SCOTCH itself. No restriction on
BELFEM's BSD distribution. Option text corrected to "(CeCILL-C, LGPL-like)".
Source: https://www.labri.fr/perso/pelegrin/scotch/

## Verification

Configure-only runs in a scratch build dir, each resolving `.so` libraries
under `/opt/scls/gcc/lib64`:

- defaults (SuperLU, HDF5, Exodus, METIS)
- `-DUSE_SUITESPARSE=ON`
- `-DUSE_MUMPS=ON -DUSE_SCOTCH=ON` — link.txt confirmed order:
  dmumps/mumps_common/esmumps/pord → scalapack/lapack/cblas/blas → mpi_mpifh…
- `-DUSE_STRUMPACK=ON -DUSE_SCOTCH=ON`
- `-DUSE_PETSC=ON`
- `-DUSE_VTK=ON`

## Audit (Codex + Grok, 2026-07-01)

Both auditors reviewed the diff (thread: `tmp/ai_exchange/cmake_find_package_audit.md`,
distilled here). Confirmed: no stale `BELFEM_IO_LIBS`/`find_tplibs` consumers;
esmumps→scotch order correct; helper mechanics sound; macOS changes net-positive.
Post-audit fixes applied:

- **Debian multiarch** (both, confirmed blocker for stock-Debian default builds):
  helper now also searches `lib/${CMAKE_LIBRARY_ARCHITECTURE}`. HDF5's split
  `hdf5/openmpi` Debian layout still needs `-DHDF5_DIR` — documented limitation.
- **`-Wl,--no-as-needed`** prepended to the dynamic MKL block (Codex; needed on
  `--as-needed`-default linkers). Codex's companion claim that `BELFEM_RPATH`
  yields no runtime rpath was REFUTED empirically: `LINK_DIRECTORIES` emits
  `-Wl,-rpath` in link.txt on this host and built binaries carry the RUNPATH.
- **arpack ordering** (both, pre-existing): `config_matrix.cmake` now included
  before `config_arpack.cmake`, so dedup keeps `-larpack` after `-larmadillo`
  (verified in link.txt).
- **`set( BELFEM_SOLVER_LIBS )`** init added (Grok, readability).
- "MPI Fortran at tail" wording corrected (see config_mumps note above).

## Rank-Ordered Link Line (follow-up, same day)

Christian provided the SCLS dependency graph with build groups (rank 0-10).
Rule adopted: **dependents link before dependencies — higher rank first;
same-rank order is free.** Implementation replaces the ad-hoc
`BELFEM_SOLVER_LIBS`/`BELFEM_MATRIX_LIBS`/`BELFEM_OTHER_LIBS`/
`BELFEM_PERFORMANCELIBS` lists entirely:

- `belfem_link_libraries( RANK libs... )` macro (in belfem_find_package.cmake)
  appends to `BELFEM_TPL_RANK_<N>` buckets.
- `finalize_compiler.cmake` assembles `BELFEM_TPL_LIBS` from rank 10 down to 0;
  `Add_Executable.cmake`/`Add_Test.cmake` link
  `${BELFEM_LIBS} ${BELFEM_TPL_LIBS} ${BELFEM_FORTRANLIBS} ${BELFEM_OPENMPLIBS}`.
- Rank map: petsc 9; strumpack stack 8 (atomic block, internal order already
  descends: strumpack → butterflypack → slate → zfp); armadillo, exodus 7;
  arpack, mumps block, netcdf 6; hdf5, parmetis, scalapack, scotch 5;
  superlu, suitesparse (not in SCLS graph; needs blas+metis), MPI-Fortran 4;
  metis 3; nlopt, tinyxml2, profiler 2; blas/lapack group + atomic MKL block 1;
  gomp/pthread/m/dl/z 0.
- METIS_LIBRARIES and EXODUS_LIBRARIES are split (parmetis 5 / metis 3;
  exodus 7 / netcdf 6).
- The earlier arpack/matrix include-order fix became obsolete (ranks handle
  it); no dedup on the TPL list — buckets are curated, and duplicate compat
  entries (armadillo line) are intentional.

Verified: all-solvers configure (`MUMPS+STRUMPACK+SUITESPARSE+PETSC+SCOTCH`)
passes; banner link.txt order matches the rank map exactly, every graph edge
satisfied. Caveat: option cache in existing build dirs may pin old
USE_SCOTCH=OFF — the new ON default only applies to fresh caches.

## make reset Obsoleted (follow-up, same day)

Root cause of the "static libs not updated by plain make" problem:
`Add_Executable.cmake`/`Add_Test.cmake` linked project libraries as free
`-lbelfem_*` strings with only `add_dependencies` for build order — no
file-level dependency, so executables never relinked after an archive was
rebuilt. Fixed by linking against the library TARGETS
(`BELFEM_LIB_TARGETS`); CMake now relinks automatically.

Two latent issues surfaced by the first real link test:

- `config_gcc.cmake` (Clang branch) smuggled `-Wl,-export_dynamic` through
  `BELFEM_LIBS`; moved to `CMAKE_EXE_LINKER_FLAGS`.
- The project archives are genuinely circular (`core` → `io::file_exists`
  while `io` precedes `core` in the reversed order). The old setup linked
  the base libs TWICE (top-level `BELFEM_LIBS` block + Add_Executable
  prepending), which masked the cycle. Now the archives are wrapped in
  `-Wl,--start-group/--end-group` on Linux (ld rescans); macOS ld64 rescans
  archives by default. The redundant top-level `BELFEM_LIBS` block was
  removed.

Verified in a scratch build dir: `make banner` links; `touch
src/core/banner.cpp && make banner` recompiles AND relinks (mtime advances);
`USE_TEST=ON` configures and `test_containers` builds and links through the
new Add_Test.cmake. `make reset` remains available but should no longer be
needed for normal development — pending confirmation on a full build.

## Open Items

- Dynamic MKL link line untested (no MKL install here) — verify on an MKL box.
- PETSc thermal executable must be compiled once to confirm `-lX11`/`-lmpfr`
  are not needed (both commented out).
- MPICH users need `-lmpifort` instead of the openmpi Fortran libs (comment in
  config_mumps.cmake).
- Deferred tidy candidates: numeric `COMPILER_ID` magic values, double
  `globals.cmake` include (finalize + top level), `BELFEM_LIBLIST` vs
  `BELFEM_LIBLIST_BASE` near-duplication, gperftools has no search logic.
- Full `make` not yet run — configure-level verification only.
