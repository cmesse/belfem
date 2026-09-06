# Handoff: SCLS ships STRUMPACK with OpenMP tasking disabled

**Date:** 2026-08-17
**For:** the SCLS-focused session ( this is a BUILD-SYSTEM issue, not a
BELFEM one )
**Status:** finding confirmed, cause not yet identified, no fix attempted
**Impact:** likely leaves parallel factorization performance on the table
in every BELFEM production run

## The finding

`/opt/scls/mkl` ships a STRUMPACK 8.0.0 whose `StrumpackConfig.h` has
**both** OpenMP tasking features compiled OUT:

```
/opt/scls/mkl/include/StrumpackConfig.h
    /* #undef STRUMPACK_USE_OPENMP_TASKLOOP */
    /* #undef STRUMPACK_USE_OPENMP_TASK_DEPEND */
```

These are not recipe choices. STRUMPACK sets them from CMake
`try_compile` probes, i.e. "does this compiler support `omp taskloop` /
`omp task depend`". **Both probes succeed on this host** — verified
directly:

```
$ g++ -fopenmp -fsyntax-only <taskloop test>   # compiles clean
```

So the shipped library was configured in an environment where the probes
FAILED, and the result was silently baked in. That is the thing to
investigate: the RPM/build environment, not STRUMPACK and not BELFEM.

## Why it matters

STRUMPACK's multifrontal factorization uses OpenMP tasking to
parallelise the elimination tree. With `TASKLOOP` and `TASK_DEPEND`
undefined, it falls back to coarser scheduling.

Measured on the reference machine ( i9-10900X, 10 physical cores ),
BELFEM tapestack3d, 2.47 M magnetic dofs, same step from the same
restart state:

| configuration | numeric factorization |
|---|---|
| 1 rank x 10 threads | 21.2 s |
| 4 ranks x 2 threads | 20.5 s |
| 4 ranks x 4 threads ( oversubscribed ) | 23.5 s |

The factorization time is **flat across a 5x change in thread count**,
then degrades under oversubscription. Assembly in the same comparison
scales 3.6-3.8x with rank count, so the machine and the measurement are
fine — it is specifically this phase that does not scale.

That flatness is consistent with a bandwidth-bound kernel. It is ALSO
consistent with tasking being compiled out. **The two hypotheses have
not been separated**, and separating them is the point of this handoff:
if it is the build flags, a rebuild is free performance in the phase
that currently dominates a 4-rank step ( ~21 s of factorization against
~9 s magnetic + ~20 s thermal assembly ).

## What already exists to work with

A local, patchable STRUMPACK build is at
`/home/christian/codes/strumpack-local/` ( built 2026-08-17 ):

- `configure_local.sh` — the full configure line, derived from the SCLS
  recipe ( MKL variant, cloned from github.com/cmesse/scls ). It matches
  the shipped library's configuration, and **deliberately forces the two
  flags OFF** with `-DSTRUMPACK_USE_OPENMP_TASK_DEPEND=FALSE
  -DSTRUMPACK_USE_OPENMP_TASKLOOP=FALSE` so the local library is a
  like-for-like stand-in. **Delete those two lines to get the build with
  tasking ENABLED** — that is the A/B.
- `install/` — the like-for-like library, `StrumpackConfig.h` verified
  byte-identical to the system one.
- `strumpack_pivot_counter.patch` — an unrelated local instrumentation
  patch ( counts replaced tiny pivots ). Keep or drop; it is orthogonal
  to this issue.

BELFEM switches between system and local with no source change:

```bash
cmake <belfem-src> -DSTRUMPACK_DIR=/home/christian/codes/strumpack-local/install
cmake <belfem-src> -USTRUMPACK_DIR      # back to /opt/scls
```

## Suggested line of attack ( for the SCLS session )

1. **Reproduce the probe failure.** Configure STRUMPACK inside the SCLS
   RPM build environment and capture `CMakeFiles/CMakeError.log` /
   `CMakeOutput.log` for the two OpenMP probes. The likely suspects are
   a compiler/flag mismatch in the RPM env ( no `-fopenmp` reaching the
   probe, a different compiler than the one used for the real build, or
   a sandbox where the probe cannot link ).
2. **Check the blast radius.** The same probe pattern may have silently
   disabled features in OTHER SCLS packages built in that environment —
   worth grepping the installed `*Config.h` / `config.h` headers under
   `/opt/scls` for `#undef` of features this host supports.
3. **Fix in the recipe environment**, then rebuild STRUMPACK and confirm
   `StrumpackConfig.h` shows both features DEFINED.
4. **Measure the payoff** with the A/B above: local build WITH tasking
   vs the system library, same BELFEM deck, comparing `factor time`
   from a verbose run. Note the methodological caveat below.

## Methodological caveat ( learned the hard way tonight )

Do **not** compare step wall-times. Changing thread counts or library
scheduling perturbs floating-point reduction orders, which moves the
nonlinear convergence path: the same step from the same dump took 4
iterates in one run and 10 in another, with identical per-iterate cost.
Compare the `factor time = ...` lines from `hphiTrun -v` directly, or
normalise by iterate count. See `doc/parallel_execution.md`.

## What this is NOT

- Not a BELFEM defect. BELFEM's own settings are documented in
  `doc/parallel_execution.md` and were measured tonight.
- Not the cause of the recent slowdown campaign ( that was DR-89, a
  tolerance-contract bug, since fixed ).
- Not urgent for correctness — only for speed.
- **Does NOT reopen the thread policy.** Christian's ruling ( 2026-08-17 )
  stands regardless of this fix: `OMP_NUM_THREADS = real_cores / ranks`,
  no hyperthreading, no oversubscription. Measured on 4 ranks: 4 threads
  each cost ~15 % factorization time AND ~13 GiB of summed footprint
  versus 2 threads. A tasking-enabled STRUMPACK changes how well the
  granted threads are used, not the cost of over-granting them.
