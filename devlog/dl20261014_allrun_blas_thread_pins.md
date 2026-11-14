# Allrun pins MKL_NUM_THREADS and OPENBLAS_NUM_THREADS to NTHREADS

**Date:** 2026-10-14
**Purpose:** Christian's question was whether `OPENBLAS_NUM_THREADS`, like `MKL_NUM_THREADS`, is a BELFEM or an SCLS item in the recurring OpenMP+MPI oversubscription trouble. He then asked for the Allrun fix.
**Module:** `examples/scripts` (`Allrun`); `todo/mkl_num_threads_allrun.md`
**AIs involved:** Claude Sonnet 5.5 (the triage), Claude Opus 5.5 (the edit and the second round); a retroactive code jury of Codex `gpt-5.6-terra`/high and Grok `grok-4.7`/high on the Allrun hunk and the R4 text, plus a Codex `gpt-5.6-luna`/medium language sweep of the R4 text
**Status:** R1, R3 and R4 done. R2/O1 (`MKL_DYNAMIC`) waits on Christian's ruling, and R5 (measurement) on a run by Christian
**Verification:** probe level. `bash -n examples/scripts/Allrun` passes. The new block was extracted and sourced under macOS `/bin/bash` 3.2.57 with `set -euo pipefail` and `NTHREADS=2` in three environments: nothing inherited; `MKL_NUM_THREADS=8`; `MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=16`. All three end with OMP/MKL/OpenBLAS = 2/2/2, and the NOTE is printed exactly for the differing inherited values. No deck was run, and no solver trace was checked (that is R5).

## Summary

**The owner is BELFEM, not SCLS.** On 2026-09-25 SCLS ruled that thread control belongs to the application launcher; see the scope guard in `todo/mkl_num_threads_allrun.md`, and scls `changelogs/scalapack.md`. Allrun is where `NTHREADS` is chosen and where ranks × threads is guarded (`Allrun:184`, Slurm `:326`). On the SCLS side, `doc/MKL_ABI_POLICY.md` gets one sentence saying the launcher owns thread control.

**OpenBLAS: corrected later the same day.** The first version of this entry said an inherited `OPENBLAS_NUM_THREADS` takes precedence over `OMP_NUM_THREADS` in the SCLS OpenBLAS. **That was wrong.** SCLS builds OpenBLAS with `USE_OPENMP=1`. That build parses `OPENBLAS_NUM_THREADS` (`driver/others/openblas_env.c`) but ignores it for the thread count: `blas_get_cpu_number` uses it only under `#ifndef USE_OPENMP` (`memory.c`), and `num_cpu_avail` takes `omp_get_max_threads()` on every call (`common_thread.h:143-178`). A probe against the installed `/opt/scls` OpenBLAS 0.3.33 confirms it: with `OMP_NUM_THREADS=2`, a DGEMM runs on 2 threads whether `OPENBLAS_NUM_THREADS` is 8, 1 or unset, and `OPENBLAS_NUM_THREADS=1` alone gives 16 (all logical CPUs). So on SCLS only `OMP_NUM_THREADS` limits OpenBLAS. The `OPENBLAS_NUM_THREADS` export stays because it is correct for a pthreads OpenBLAS outside SCLS, such as a distribution's. The MKL half holds: `MKL_NUM_THREADS` does override `OMP_NUM_THREADS`.

## What landed

`examples/scripts/Allrun`, after `export OMP_NUM_THREADS="$NTHREADS"`: a loop over `MKL_NUM_THREADS` and `OPENBLAS_NUM_THREADS` that exports each as `$NTHREADS`. It prints `NOTE: overriding inherited <var>=<old> with <NTHREADS>` when an inherited value was set and differed. It reads the inherited value with `${!var:-}`, which is safe under `set -u` on Bash 3.2. Every non-MKL, non-OpenBLAS flavor gets two harmless extra variables.

`MKL_DYNAMIC` was deliberately not touched. Choosing it is O1 in the todo file, and that ruling is Christian's.

## Open Questions

- O1 / R2: whether `MKL_DYNAMIC` should be pinned to `FALSE` (Christian).
- R5: run a deck on the `mkl` flavor (and on an OpenBLAS flavor) with a large pre-exported value, and check that the solver traces show `NTHREADS`.
- The jury's further findings, listed below, are waiting on Christian's decision.

## Files Updated

- `examples/scripts/Allrun`
- `todo/mkl_num_threads_allrun.md` (R1 and R3 ticked, Status rewritten, OpenBLAS added to the scope)
- `todo/README.md` (the index line)
- `devlog/README.md` (the index line)
- `devlog/dl20261014_allrun_blas_thread_pins.md`

## Jury round (retroactive) and the R4 doc text

cdabf2e went in without the jury that `code-with-jury` requires for a shell script, so a retroactive round ran on the Allrun hunk plus the R4 doc paragraph: Codex `gpt-5.6-terra`/high and Grok `grok-4.7`/high, blind, slug `review_allrun_blas_thread_pins`. The first dispatch was lost to a host reboot and was re-run unchanged.

- **OpenBLAS wording:** Codex confirmed the correction above and narrowed "never reads" to "ignores for the thread count". Grok argued, from the variable's name appearing in the dylib, that the SCLS OpenBLAS honours it. The probe above refutes that. Grok was right that the installed library is 0.3.33, not 0.3.34.
- **`MKL_DOMAIN_NUM_THREADS`** (all three): a per-domain MKL setting outranks `MKL_NUM_THREADS`, so the pin is not complete against it. Documented in `doc/parallel_execution.md`; whether Allrun should act on it is Christian's call.
- **Bash 3.2:** Grok's concern about `${!v:-}` under `set -u` was refuted by running it on `/bin/bash` 3.2.57.
- **Single-raiser, pre-existing:**
  - Multi-node hand `mpirun` needs `-x` (Codex; now in the doc).
  - `SLURM_EXPORT_ENV=NONE` defeats the `srun` path (Codex).
  - `hatch_turtle()` checks only `OMP_NUM_THREADS` (Grok, `cl_SolverWrapper.cpp:440`).
  - `doc/coding_philosophy.md:718` cites `hatch_turtle` at `:195-245`; it is at `:391` (Grok).
  
  These were reported to Christian for a decision and not deferred.
- **`MKL_DYNAMIC` (O1):** the facts from both auditors match the pre-registration. Christian rules on them; they are not settled by vote.

R4 landed: `doc/parallel_execution.md` has the BLAS variables in the recommendation block and a paragraph under "Launching", and was given a Codex `gpt-5.6-luna`/medium language sweep.

**Verification (this round):** probe level. The OpenBLAS thread-count probe and the Bash 3.2 probe were executed. No deck was run.
