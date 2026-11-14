# Allrun: pin MKL_NUM_THREADS and OPENBLAS_NUM_THREADS (and decide MKL_DYNAMIC) next to OMP_NUM_THREADS

**Date:** 2026-09-25
**Purpose:** Make the thread count that `examples/scripts/Allrun` chooses (`NTHREADS`) govern MKL and OpenBLAS as well, even when the user's or site's environment exports `MKL_NUM_THREADS` or `OPENBLAS_NUM_THREADS`, so that the oversubscription guard cannot be silently defeated on the SCLS `mkl` flavor.
**Module:** `examples/scripts` (+ `doc/parallel_execution.md`)
**AIs involved:** Claude (belfem host session, finding); Claude Opus 5.5 (R1/R3, 2026-10-14)
**Status:** OPEN — R1, R3, R4 done 2026-10-14 (`devlog/dl20261014_allrun_blas_thread_pins.md`); R2 waits on O1 (Christian); R5 waits on a run by Christian

> **Scope guards:**
> - SCLS deliberately does NOT pin MKL threads in the packages. Christian ruled (2026-09-25) that thread control belongs in BELFEM's Allrun. No SCLS change here.
> - Only the environment Allrun exports, plus documentation. No solver or kernel change.

## Background

On 2026-09-25 SCLS unified the `mkl` flavor to a single MKL threading layer (gf_lp64 + gnu_thread + core; scls commit 6a4d268). Before that, `libscalapack` linked `libmkl_sequential` and `libarmadillo` linked `libmkl_rt` (SCLS devlog `dl20260925_mkl_threading_uniformity.md`). As a result, **ScaLAPACK's MKL BLAS calls, which used to be single-threaded, now run threaded**, and so do MUMPS, STRUMPACK, SLATE, PETSc and SLEPc through it.

What devel already does (checked at 68a4ba0, confidence high):
- `examples/scripts/Allrun:396` exports `OMP_NUM_THREADS="$NTHREADS"`. `libmkl_gnu_thread` honours `OMP_NUM_THREADS`, so in a clean environment MKL already follows `NTHREADS`.
- Allrun refuses ranks × NTHREADS > cores (`examples/scripts/Allrun:184`, Slurm path `:326`).
- `examples/{corc_garber_effect,costheta,helix}/Allrun` are thin `exec ../scripts/Allrun` wrappers and inherit this.
- `doc/parallel_execution.md` covers `OMP_NUM_THREADS` (`:10`, `:21`, `:30`, `:48`) but never mentions `MKL_NUM_THREADS` or `MKL_DYNAMIC`.

## Gap

**`MKL_NUM_THREADS` takes precedence over `OMP_NUM_THREADS` for MKL.** If it is inherited (Intel `setvars.sh`, module files and site profiles set it), MKL ignores `NTHREADS`, and `MPI_RANKS × MKL threads` can oversubscribe the node with no warning from the guard. Hand-launched runs outside Allrun are already a documented failure mode (INC-345, `devlog/dl20260816_bdf_history_restore.md:90-97`).

## Steps

- [x] **R1** In `examples/scripts/Allrun` next to `:396`, export `MKL_NUM_THREADS="$NTHREADS"`. Print a one-line NOTE if an inherited value differed, so users see that it was overridden. Extended 2026-10-14 to `OPENBLAS_NUM_THREADS`. It matters only for a pthreads OpenBLAS: the SCLS OpenMP build ignores it for the thread count (probe, see the devlog).
- [ ] **R2** Decide `MKL_DYNAMIC`. The default `TRUE` lets MKL use fewer threads than asked, which changes reduction order. That matters for the determinism statement in `devlog/dl20260810_examples_and_rho_database_repair.md:140-148`. Options: leave it, export `FALSE`, or document only. → **O1**
- [x] **R3** Verify that the Slurm path (`NTHREADS="${SLURM_CPUS_PER_TASK:-$NTHREADS}"`, `:312`) reaches the same export before `run_solver`. Reviewed 2026-10-14 (static): the `sbatch` branch exits before the export; the re-entered `--slurm` run passes the export block before `run_solver srun`, and `srun` propagates the environment by default. Not run on a cluster.
- [x] **R4** `doc/parallel_execution.md`: add `MKL_NUM_THREADS` / `MKL_DYNAMIC` to the policy section, and state that hand-launched `mpirun` must set `OMP_NUM_THREADS`, `MKL_NUM_THREADS` and `OPENBLAS_NUM_THREADS`.
- [ ] **R5** Measure once on the `mkl` flavor: run a deck with `MKL_NUM_THREADS=8` pre-exported and check that after R1 the STRUMPACK/MUMPS traces show `NTHREADS`, not 8.

## Open questions

- **O1** `MKL_DYNAMIC`: pin `FALSE` for reproducibility, or leave MKL's default? (Christian)

## Definition of done

Allrun-launched runs use exactly `NTHREADS` MKL threads regardless of the inherited environment (R5 evidence), and the policy doc says so.
