# Parallel Execution: Ranks, Threads and the Allocator {#doc_parallel_execution}

**Date:** 2026-08-17
**Purpose:** How to choose MPI ranks, OpenMP threads and allocator settings
for a production run, and why those are the right choices

## The recommendation

```bash
export OMP_NUM_THREADS=$(( REAL_CORES / MPI_RANKS ))
export MALLOC_ARENA_MAX=2
mpirun --stream-buffering 0 -np $MPI_RANKS stdbuf -oL belfem >> out.txt
```

For the binding flags a hybrid run should add, see "Launching" below.

Three rules, in order of how much they cost you when broken:

1. **Give the assembly ranks, not threads.** BELFEM's element loop carries
   no OpenMP pragmas: it is serial within a rank. On a default build
   `OMP_NUM_THREADS` accelerates only the third-party numerics (STRUMPACK,
   MKL, threaded BLAS, PARDISO) — BELFEM's own three OpenMP kernels are
   gated on `BELFEM_OMP`, which is OFF (see "BELFEM's own OpenMP" below).
   A run with one rank and ten threads leaves the dominant cost entirely
   unparallelised.
2. **Count real cores, not hyperthreads.** Both the assembly (memory-bound
   element loops) and the factorization (BLAS3 kernels already saturating
   the vector units) gain nothing from a second logical thread on the same
   physical core, and the extra context costs cache.
3. **Do not oversubscribe.** `MPI_RANKS x OMP_NUM_THREADS` must not exceed
   the real core count. MPI ranks and OpenMP threads compete for the same
   cores; oversubscription turns a saturated machine into a thrashing one.

Example, a 10-core / 20-thread workstation ( the reference machine ):

| ranks | threads/rank | total | verdict |
|---|---|---|---|
| 4 | 2 | 8 | recommended, measured below |
| 8 | 1 | 8 | more assembly parallelism, factorization unchanged |
| 1 | 10 | 10 | worst case: assembly entirely serial |
| 4 | 5 | 20 | oversubscribed via hyperthreads — do not |

## Launching: bind explicitly, then check placement

The three rules above choose `ranks x threads`. They do not place those threads on
cores, and launcher and runtime defaults are not a placement policy.

**Nothing sets `OMP_NUM_THREADS` for you.** With the variable unset, both GNU
libgomp and LLVM/Intel libomp take the popcount of the process affinity mask,
which counts every hyperthread. Neither runtime inspects SMT, and the OpenMP
standard leaves the initial value of `nthreads-var` implementation-defined. So an
unset variable is not a recommendation; it is whatever the mask happened to
contain.

**Open MPI binds, but not the way a hybrid run needs.** Its documented defaults
are bind-to-core for two processes or fewer, **bind-to-package** above that, and
bind-to-none when oversubscribed. On a single-package workstation "package" is
the whole machine, so every rank's mask is every logical CPU and every OpenMP
runtime independently concludes it owns the node. That is how ten ranks come to
request twenty threads on ten cores without anything having chosen it.

One BELFEM-specific wrinkle: a pre-main constructor in `cl_Communicator.cpp`
exports binding-policy-`none` variables for Open MPI, PRRTE and Hydra (see
`src/comm/doc/comm_usage_guide.md`). Under `mpirun` the launcher has already
applied its binding before the executable starts, so those exports are believed
inert there, but the interaction with explicit `--bind-to` flags has not been
measured. `--report-bindings` settles it per machine; trust that output over
either default.

Bind on purpose instead. Use one launcher family per run; `srun --cpu-bind` and
`mpirun --map-by ...:PE=` can fight if combined:

```bash
# workstation, one package
mpirun -np $NRANKS --map-by slot:PE=$OMP_NUM_THREADS --bind-to core \
       --report-bindings ...

# multi-socket node: map by the locality object so a rank's threads
# stay inside one package
mpirun -np $NRANKS --map-by package:PE=$OMP_NUM_THREADS --bind-to core \
       --report-bindings ...

# under Slurm, let Slurm do it ( this is what examples/scripts/Allrun does )
sbatch --ntasks-per-node=$RANKS_PER_NODE --cpus-per-task=$NTHREADS ...
srun --mpi=pmix -n $NRANKS ...
```

`--report-bindings` prints the placement it actually chose; use it at least once when
bringing up a new machine instead of trusting the flags.

**A Slurm portability trap:** sites differ on whether `--cpus-per-task` counts
logical or physical CPUs when SMT is on, so treating `SLURM_CPUS_PER_TASK`
as `OMP_NUM_THREADS` is not portable as a blind rule. Check with
`--report-bindings` or `lscpu -p=CORE,SOCKET` on the target machine before
trusting it.

## Why: the measured split

One coupled timestep of the tapestack3d study ( 2.47 M magnetic +
0.90 M thermal dofs ), the SAME step re-run from the same restart state,
serial versus four ranks — one variable changed:

| quantity | 1 rank x 10 thr | 4 ranks x 2 thr | 4 ranks x 4 thr |
|---|---|---|---|
| magnetic assembly | 33.3-34.1 s | 9.1-9.4 s | 9.1-9.4 s |
| thermal assembly | 74.8-75.1 s | 19.9 s | 19.9 s |
| numeric factorization | 21.2 s | 20.5 s | **23.5 s** |
| cheap re-solves | 13.8 s | 4.3 s | - |
| per-iterate cost | ~150 s | ~73 s | ~70 s |

The third column is the oversubscription control ( 4 x 4 = 16 threads on
10 physical cores, run 2026-08-17 ):

- **Assembly is bit-for-bit the same as the 2-thread run.** Threads
  cannot touch it — this is the cleanest confirmation that the element
  loop is serial per rank.
- **The factorization gets ~15 % SLOWER** ( 23.5 s vs 20.5 s ). Together
  with the flat 21.2 s at ten serial threads, that says the multifrontal
  phase is bandwidth-bound: past the physical cores, extra threads add
  contention, not throughput. Hyperthreading does not rescue it.
- **Oversubscription also costs MEMORY**: the 4 x 4 run peaked near
  54 GiB summed against ~41 GiB for 4 x 2 at the same steps — the
  per-rank solver workspace scales with the thread count. Ruled
  ( Christian, 2026-08-17 ): **no thread oversubscription, period** —
  and the ruling explicitly survives any future STRUMPACK rebuild with
  OpenMP tasking repaired, because better task scheduling improves what
  the granted threads do; it does not repeal the contention and
  workspace cost of granting too many.

Read the table as three separate facts:

- **Assembly scales almost ideally with ranks.** This is the dominant
  term ( ~75 % of a serial step ) and the entire reason to run in
  parallel. A `perf` profile of the assembly phase is 94 % BELFEM's own
  code and essentially 0 % BLAS — there is nothing there for threads to
  accelerate.
- **The factorization does not care.** It was already saturating the
  machine through OpenMP in the serial run ( 537 / 573 GFlop/s observed
  in the multifrontal phase ), so redistributing it across ranks moves
  the work without shrinking it. Do not expect rank count to buy solver
  time.
### Measure phases, not steps

**Step wall-time is a NOISY metric for comparing settings.** Changing the
thread count changes MKL/STRUMPACK reduction orders, which perturbs the
last bits of the residual, which moves the convergence path. (Until
2026-09-01 BELFEM's own `matvec_csc` reduction reassociated too; with
`BELFEM_OMP` OFF that kernel is now bitwise deterministic across thread
counts, so what remains is third-party.) The same
step from the SAME restart dump took 4 iterates ( 4 min 53 s ) in one
run and 10 iterates ( 11 min 40 s ) in another, with per-iterate cost
identical at 70-73 s. Anyone A/B-ing a setting on step wall-time alone
would have concluded the second configuration was 2.4x worse; it was
not. Compare per-phase timings ( the Verbose `time for computing
Jacobian and residual` and `factor time` lines ) or normalise by the
iterate count.

The cleanest overall number comes from two draws that happened to need
the SAME iterate count: the serial warm start and the 4x4 run both took
10 iterates on the identical step, at 26 min 33 s versus 11 min 38 s —
an honest **2.3x overall speedup at four ranks**. Quote that, not the
5.4x ( a lucky 4-iterate draw ) and not the 3.6-3.8x ( assembly phase
only, diluted by the non-scaling factorization ).

Between-run variance is also LARGER than pure summation-order noise
when the RANK COUNT changes: partitioning, elimination ordering and the
parallel tree all change with it, so the linear solves deliver
structurally different ( all tolerance-legal ) error vectors. Thread
count alone reassociates sums; rank count changes the arithmetic path.

- **Part of any speedup you measure is convergence, not parallelism.**
  In the comparison above the step also needed 4 magnetic iterates
  instead of 10. The honest parallel speedup is the ~3.7x on assembly;
  the rest came from solver-tolerance work ( see the debt register ).

## Memory: ranks redistribute, they do not reduce

The sparse factor is the ceiling and it is **rank-independent**:
21.2 GiB serial, 21.2 GiB across four ranks. What changes with rank
count:

- Per-rank mesh and matrices shrink ( partitioning works ), but the
  TOTAL grows, because every rank carries an aura of ghost entities and
  parallel runs additionally build assembly tables ( ~0.5 GiB on this
  problem ).
- Measured on the reference machine: serial peaked at 47-51 GiB;
  four ranks peaked near 41 GiB summed with 13-17 GiB free.

So more ranks is not a memory remedy. If a run does not fit, the levers
are the problem size, the solver's compression settings ( see
`src/sparse/doc/solver_memory_and_compression.md` ) or out-of-core, not
the rank count.

NOTE on reading the memory printout: `SolverData` SUMS each entry across
ranks before printing. "Mesh : 6 GiB" on a four-rank run means 6 GiB in
total, roughly 1.5 GiB per rank — not 6 GiB per rank.

## The allocator settings

**`MALLOC_ARENA_MAX=2`.** glibc gives concurrent threads independent
heap arenas ( default ceiling 8 x cores ), and memory freed into one
arena is neither returned to the OS nor available to the others. Under a
workload that churns hundreds of MiB per linear solve, that stranding
lets the resident set ratchet upward while nothing is leaked. Capping the
arena count pools the free memory instead.

Honest status: this setting is **recommended but its contribution is not
isolated**. The measurement that demonstrably cured the ratchet was the
`malloc_trim()` call the Controller now makes after every accepted
timestep ( 13.7 GiB released at a step boundary; the resident floor fell
from 46 to 32 GiB ), and that measurement was taken on a run that did
NOT export `MALLOC_ARENA_MAX`. Keep it — the lock-contention cost at a
few threads per rank is negligible — but do not credit it for a healthy
run without an A/B.

**`stdbuf -oL`.** Not performance: BELFEM mixes `std::cout` with
`fprintf(stdout)`, and a redirected run is fully buffered, so `out.txt`
looks empty for minutes. Line-buffering makes progress visible.

**`--stream-buffering 0`.** Also not performance, and not the same thing.
`stdbuf -oL` only helps output that ends in a newline. The progress bar
redraws in place with `\r` and never writes one, so line buffering cannot
release it — and on top of that Open MPI holds the child's stdout in its
I/O forwarder until the job ends. The result is that the bar shows nothing
at all and then a finished 100 % bar. Measured on Open MPI 5.0.10 with a
40-frame bar: without the flag the whole bar is delivered in one chunk at
exit, with it in roughly ninety. It must be the command-line option —
exporting the MCA variable the flag sets in the rank
(`OMPI_MCA_ompi_stream_buffering=0`) reaches the process but does not
stream, so the flag evidently also carries a forwarder-side setting.
`srun` has no equivalent, which is why `examples/scripts/Allrun` adds this
only on its `mpirun` path.

## Related

- `doc/mpi_support.md` — which MPI implementation is supported, and why
- `src/sparse/doc/solver_memory_and_compression.md` — solver-side memory,
  compression and the per-library tolerance contracts


## BELFEM's own OpenMP, and `OMP_STACKSIZE`

**BELFEM's own hand-written OpenMP is off on a default build.** Three Fortran
kernels — `src/sparse/splinalg.f90` (`matvec_csr`, `matvec_csc`),
`arpacktools.f90` and `parpacktools.f90` (the ARPACK reverse communication
matvecs) — carry `!$omp` directives, but they are gated on the `BELFEM_OMP`
define, which comes from `USE_BELFEM_OPENMP` and is **OFF**.

**That does not make BELFEM's objects OpenMP-free.** The matrix backend is
header-only and parallelizes its own expression evaluation whenever OpenMP is
available: Armadillo because BELFEM asks it to
(`src/linalg/armadillo/armadillo.hpp`, `#define ARMA_USE_OPENMP` under
`#ifdef OMP`), Blaze on its own (`blaze/system/SMP.h`, OpenMP mode whenever
`_OPENMP` is defined). Those pragmas compile *into BELFEM translation units*.
Measured on a default Darwin release build: **129 of 339 objects reference
`GOMP_*`.** Both platforms are affected — Linux defaults to Armadillo, Apple to
Blaze. So a `gomp_thread_start` frame with BELFEM code above it does **not**
imply a BELFEM pragma.

This matters for *reading a stack trace*, not for the stack-overflow failure
below: neither backend uses an OpenMP reduction, so neither creates the
per-thread private array that caused it.

This is a separate switch from `USE_OPENMP`, which stays **ON**: the
third-party solvers need `-fopenmp`, and the places that query or set their
thread budget stay on the plain `OMP` define. Those are deliberate and are
commented as such in the source — `hatch_turtle()`
(`src/sparse/cl_SolverWrapper.cpp`), the startup banner's `Threads Used` line
(`src/core/banner.cpp`) and PARDISO's `gParameters(3)`
(`src/sparse/pardisotools.f90`). Do not "unify" them onto `BELFEM_OMP`:
that would delete the oversubscription warning and silently serialize PARDISO.

**Why OFF.** Both matvecs run master-only — the residual
(`cl_FEM_DofMgr_SolverData.cpp`, inside an `is_master()` block) and the
eigenvalue power iteration (`cl_FEM_DofMgr_EigenValues.cpp`) — and at roughly
two flops per nonzero they are noise beside assembly and a direct
factorization. Against that, `matvec_csc`'s `reduction(+:y)` gave each worker
thread a **private copy of the whole result vector on that thread's stack**.
An OpenMP worker gets 2 MiB on Darwin, so the run died with SIGBUS above

    n = 2 MiB / 8 bytes = 262,144 rows

which ordinary 3D problems cross easily. The threading never paid for that
exposure, so it was switched off rather than repaired (DR-155).

**`OMP_STACKSIZE` — read this before turning `BELFEM_OMP` back on.**

- **Turning `BELFEM_OMP=ON` re-arms the defect.** The private copy is still on
  the worker stack. Either raise `OMP_STACKSIZE` for every run, or move the
  scratch to the heap first. `OMP_STACKSIZE=64M` was verified to carry a
  295,315-dof case that crashed at the default.
- **`OMP_STACKSIZE` is a diagnostic knob with `BELFEM_OMP` OFF, not a
  requirement.** The failure this section describes needs an OpenMP *array*
  reduction, which materializes one n-sized private copy per thread on that
  thread's stack. **The matrix backend cannot hit it**: neither Blaze nor
  Armadillo uses an OpenMP reduction of any kind — Blaze partitions the
  destination so each thread writes its own slice of shared memory, Armadillo
  uses element-wise loops with `critical`/`atomic`. So gating `BELFEM_OMP` off
  genuinely removes this failure class from BELFEM's own code; the backend's
  threads are not a second route to it. The third-party solvers have not been
  audited for it and no BELFEM run has faulted that way — so if a kernel ever
  does fault inside a `gomp_thread_start` frame, raising `OMP_STACKSIZE` is
  worth one try before assuming a BELFEM bug, but nothing here says it will be
  needed. If a third-party kernel ever faults inside a
  `gomp_thread_start` frame, raise it before assuming a BELFEM bug.
- The limit is per worker thread and is **not** the 8 MiB main-thread stack.
  `OMP_NUM_THREADS=1` therefore hides the whole class of failure, because the
  work then runs on the master thread — useful as a diagnostic, not as a fix.

**If you want the threading back**, measure it first: build with
`-DUSE_BELFEM_OPENMP=ON`, run the same deck both ways at equal
`OMP_NUM_THREADS` with `OMP_STACKSIZE` raised, and compare **per-phase**
timings rather than step wall-time (see "Measure phases, not steps" above).
No such measurement exists in the repository; the kernels were threaded on the
reasoning that a scatter needs it, never on a number.
