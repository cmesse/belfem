# Allrun pins MKL_NUM_THREADS and OPENBLAS_NUM_THREADS to NTHREADS

**Date:** 2026-10-14
**Purpose:** Christian's question was whether `OPENBLAS_NUM_THREADS`, like `MKL_NUM_THREADS`, is a BELFEM or an SCLS item in the recurring OpenMP+MPI oversubscription trouble. He then asked for the Allrun fix.
**Module:** `examples/scripts` (`Allrun`); `todo/mkl_num_threads_allrun.md`
**AIs involved:** Claude Sonnet 5.5 (the triage), Claude Opus 5.5 (the edit); no auditor seat, because this is a shell export block
**Status:** COMPLETE for R1 (both variables) and R3; R2/O1 (`MKL_DYNAMIC`), R4 (`doc/parallel_execution.md`) and R5 (measurement) open
**Verification:** probe level. `bash -n examples/scripts/Allrun` passes. The new block was extracted and sourced under macOS `/bin/bash` 3.2.57 with `set -euo pipefail` and `NTHREADS=2` in three environments: nothing inherited; `MKL_NUM_THREADS=8`; `MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=16`. All three end with OMP/MKL/OpenBLAS = 2/2/2, and the NOTE is printed exactly for the differing inherited values. No deck was run, and no solver trace was checked (that is R5).

## Summary

**The owner is BELFEM, not SCLS.** On 2026-09-25 SCLS ruled that thread control belongs to the application launcher; see the scope guard in `todo/mkl_num_threads_allrun.md`, and scls `changelogs/scalapack.md`, `todo/backlog.md:243`. Allrun is where `NTHREADS` is chosen and where ranks × threads is guarded (`Allrun:184`, Slurm `:326`). The SCLS side only owes documentation: its records mention only MKL, and OpenBLAS should be added there. That edit is not part of this commit.

**Why OpenBLAS needs its own pin.** SCLS builds OpenBLAS with `USE_OPENMP=1` (scls `recipes/openblas.yaml:38`). That build honours `OMP_NUM_THREADS`, so a clean environment already follows `NTHREADS`. An inherited `OPENBLAS_NUM_THREADS` takes precedence over it, however, exactly as `MKL_NUM_THREADS` does for MKL. Either one can set the BLAS thread count per rank past the guard without any warning (confidence high for MKL, which the todo file documents; medium-high for OpenBLAS, from the OpenBLAS env precedence and not measured here).

## What landed

`examples/scripts/Allrun`, after `export OMP_NUM_THREADS="$NTHREADS"`: a loop over `MKL_NUM_THREADS` and `OPENBLAS_NUM_THREADS` that exports each as `$NTHREADS`. It prints `NOTE: overriding inherited <var>=<old> with <NTHREADS>` when an inherited value was set and differed. It reads the inherited value with `${!var:-}`, which is safe under `set -u` on Bash 3.2. Every non-MKL, non-OpenBLAS flavor gets two harmless extra variables.

`MKL_DYNAMIC` was deliberately not touched. Choosing it is O1 in the todo file, and that ruling is Christian's.

## Open Questions

- O1 / R2: whether `MKL_DYNAMIC` should be pinned to `FALSE` (Christian).
- R4: `doc/parallel_execution.md` should name the two variables and say that a hand-launched `mpirun` must set them. It is user-facing, so it needs the Codex language sweep.
- R5: run a deck on the `mkl` flavor (and on an OpenBLAS flavor) with a large pre-exported value, and check that the solver traces show `NTHREADS`.
- SCLS: extend `todo/backlog.md:243` and `doc/MKL_ABI_POLICY.md` with OpenBLAS (a docs-only change; no pin in `recipes/environment.yaml`).

## Files Updated

- `examples/scripts/Allrun`
- `todo/mkl_num_threads_allrun.md` (R1 and R3 ticked, Status rewritten, OpenBLAS added to the scope)
- `todo/README.md` (the index line)
- `devlog/README.md` (the index line)
- `devlog/dl20261014_allrun_blas_thread_pins.md`
