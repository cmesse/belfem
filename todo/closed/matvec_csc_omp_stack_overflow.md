# `matvec_csc` overflows the OpenMP worker stack above 262,144 rows

**Date:** 2026-09-01 (rev. 1 jury round; rev. 2 oversight sweep; rev. 3 O1/O3 resolved;
**rev. 4 — RESOLVED by removing the parallelism**)
**Purpose:** `matvec_csc`'s `!$omp reduction(+:y)` places each worker thread's private copy of the
whole result vector on that thread's stack. libgomp gives a worker 2 MiB on Darwin, so any
multithreaded sparse CSC product with more than 262,144 rows walks off the guard page and the run
dies with SIGBUS/SIGSEGV inside the kernel. Move the per-thread scratch to the heap — **and
reproduce the reduction's semantics while doing it, which the first draft of R1 did not.**
**Module:** `src/sparse` (`splinalg.f90`, `cl_SpMatrix.cpp`)
**Debt row:** DR-155 (`todo/debt_register.md`)
**AIs involved:** Claude (found, isolated, measured; pre-registration), Codex `gpt-5.6-terra`/high
and Grok `grok-4.6`/high (jury audit 2026-09-01, `tmp/ai_exchange/review_matvec_csc_omp_stack_overflow.md`)
**Status:** **FIXED IN THE WORKING TREE, NOT YET COMPILED OR TESTED** (2026-09-01, Christian's
ruling). The defect is closed by **removing the parallelism rather than repairing it**: BELFEM's own
`!$omp` directives are now gated on a new `BELFEM_OMP` define which is **OFF by default**, so
`matvec_csc` no longer creates a per-thread private copy of `y` at all. The heap-buffer design
(R1, O1's communicator buffer, D14's initialisation fork) is **superseded and not implemented** —
kept below as the record of why it was hard, which is part of why it was not chosen.
**Owed:** a build and `make check` (Christian runs builds), plus the R5 deck gate. See §5.

> **Scope guards (CORRECTED 2026-09-01 by the jury round — the first two were wrong):**
> - In scope: the single `reduction(+:y)` site at `splinalg.f90:103`. It is the **only** OpenMP
>   array reduction in the Fortran tree (`grep -rnE "reduction\s*\(" src/ --include=*.f90` returns
>   one hit), so the kernel-side defect is contained to this one subroutine.
> - Out of scope: `matvec_csr`, which gathers into a private **scalar** (`splinalg.f90:35,52-56`)
>   and cannot hit this.
> - ~~Out of scope: MKL builds — the kernel is inside the `#else` of `BELFEM_MKL` at
>   `cl_SpMatrix.cpp:1650,1665`, so `USE_MKL=ON` never reaches it.~~ **FALSE — see D2.** That guard
>   covers only the 5-argument overload. The 2-argument
>   `SpMatrix::multiply( const Vector<real> &, Vector<real> & )` (`cl_SpMatrix.cpp:1737-1774`) has
>   no `BELFEM_MKL` branch at all and calls `matvec_csc` unconditionally at `:1747`. **MKL builds
>   crash on the eigenvalue path.** Widened into this plan per Christian's ruling of 2026-09-01:
>   one DR-155 sweep, not a new row.
> - ~~Not a Darwin-only bug. Darwin hits it first because its worker stack is 2 MiB; Linux's 8 MiB
>   default moves the same cliff to 1,048,576 rows rather than removing it.~~ **UNEARNED — see O2.**
>   The 8 MiB figure was assumed, not measured, and there are two incompatible readings of what
>   replaces it. It is possible the Linux cliff is the *same* 262,144 rows.
> - A second, in-tree reason Darwin hit it first, independent of stack size: MKL is a hard
>   configure error on macOS (`config/linalg/config_mkl.cmake:2-3`), so every Darwin build takes
>   the Fortran `#else` on the 5-argument path.

---

## 1. Current Behaviour and How It Fails

`matvec_csc` (`src/sparse/splinalg.f90:62-115`) parallelises the column loop. CSC is a scatter, so
every thread writes into overlapping entries of `y` (`:110`). The kernel resolves that with an array
reduction (`splinalg.f90:100-104`):

```fortran
!$omp parallel do &
!$omp shared(n,m,nnz,x,pointers,values,indices,off,base) &
!$omp private(a,b,xj) &
!$omp reduction(+:y) &
!$omp schedule(static)
```

The comment at `splinalg.f90:94-99` states the design intent exactly, and it is sound: an array
reduction gives each thread a private `y` and combines once, instead of one atomic per nonzero,
which serialised the inner loop. It even names the cost — *"n * 8 bytes per thread of scratch"*.

What it does not say is **where** those bytes live. `y` is an explicit-shape dummy
(`dimension(n), intent(inout)`, `splinalg.f90:77`), so gfortran materialises each thread's private
copy in that thread's stack frame. The scratch is stack, not heap, and nothing bounds it against
the thread stack limit.

**The correctness of the current kernel rests on two separate zeroings, and only one of them is
visible in the source.** `y = 0.0d0` at `:92` zeroes the caller's buffer; `reduction(+:)` then
identity-initialises each *private* copy to 0 before the region and adds those copies back into the
already-zeroed original. Any hand-written replacement must supply the second one explicitly (D1).

| Failure | Mechanism | Evidence |
|---|---|---|
| Hard crash — SIGBUS / SIGSEGV inside `matvec_csc._omp_fn.0` on a `gomp_thread_start` worker | private `y` of `n * 8` bytes exceeds the 2 MiB worker stack | crash report faulting thread, and the reproducer in §1.1 |
| Same crash on MKL builds, eigen path only | 2-arg overload has no MKL branch (D2) | `cl_SpMatrix.cpp:1747`, callers at `cl_FEM_DofMgr_EigenValues.cpp:1670,1706,1712,1732-1733,1755-1759` |
| Linux behaviour | unmeasured; two incompatible predictions | see O2 |

### 1.1 Measurement

Isolated reproducer: a standalone driver calling `matvec_csc` directly through the
`fspblas.hpp:38-49` binding with a diagonal matrix, so the only variable is `n`. Built against
`/tmp/belfem/lib/libbelfem.dylib`, macOS 15.7.9 x86_64, gcc 15 / libgomp, `OMP_NUM_THREADS=16`:

| n | private `y` per thread | result |
|---|---|---|
| 261,000 | 2039.1 KiB | ok |
| 262,000 | 2046.9 KiB | ok |
| **262,144** | **2048.0 KiB = 2 MiB exactly** | **crash** |
| 263,000 | 2054.7 KiB | crash |
| 295,315 (`examples/tape_quench_usermat`) | 2307.1 KiB | crash |

Two controls, both at n = 295,315, both survive:

- `OMP_STACKSIZE=64M` → completes, correct result.
- `OMP_NUM_THREADS=1` → completes. The reduction then runs on the master thread, which has the
  8 MiB main stack.

`OMP_NUM_THREADS=1` alone would confound thread count with stack size; `OMP_STACKSIZE=64M` at 16
threads varies only the stack. Together the two controls isolate the mechanism.
**Confidence: high** on the mechanism.

**Calibration correction (2026-09-01).** The sweep brackets the boundary to
(262,000 → ok, 262,144 → crash] — a 144-row, 1152-byte window. No point between 262,001 and
262,143 was run, so the earlier phrasing *"lands on 262144 to the byte"* and *"not an estimate"*
claimed a resolution the sweep does not have. That the window's upper edge coincides with
2 MiB / 8 B remains strong evidence for the mechanism. Either bisect the window (8 runs) or say
"within 1.2 kB of 2 MiB". The same overclaim propagated to `todo/debt_register.md` and
`devlog/README.md` and is corrected in the register with this revision.

### 1.2 How it presents in a real run

`examples/tape_quench_usermat` on this machine, meshed with `gmsh -3` and run with the installed
`belfem`, reaches 295,315 free dofs and dies in the first BDF1 step:

```
libbelfem.0.9.0.dylib    matvec_csc._omp_fn.0
libgomp.1.dylib          gomp_thread_start
libsystem_pthread.dylib  _pthread_start
```

(`~/Library/Logs/DiagnosticReports/belfem-2026-09-01-010821.ips`, EXC_BAD_ACCESS / SIGBUS.)

The call that dies is the residual, `mSystemMatrix->multiply( mFieldValues, mRhsVector, 1.0, -1.0 )`
(`cl_FEM_DofMgr_SolverData.cpp:2368`, `:2386`) — once per Newton step, **not** inside a Krylov
loop: MUMPS and STRUMPACK are direct, and ARPACK carries its own kernel
(`arpacktools.f90:228-237`). The genuinely per-iteration caller of `matvec_csc` is the eigenvalue
power iteration through the 2-argument overload (D2).

**Bottom line:** the kernel's scratch allocation strategy is coupled to an OS thread-stack limit
that nothing in BELFEM sets, checks, or documents, and the coupling is invisible until a problem
crosses a dof count that ordinary 3D work crosses easily.

## 2. Resolution: Remove the Parallelism, Do Not Repair It

**Chosen 2026-09-01 (Christian).** The threading in these kernels never earned the exposure:

- **The matvecs are master-only.** `mSystemMatrix` is the full `mNumberOfFreeDofs` square
  (`cl_FEM_DofMgr_SolverData.cpp:396-400`) and the residual at `:2368` sits inside
  `if ( mKernel->is_master() )` at `:2308`; the eigen matvecs are inside
  `if( mParent->parent()->is_master() )` (`cl_FEM_DofMgr_EigenValues.cpp:1704`). So the kernel is
  not MPI-parallel, and OpenMP was its only parallelism.
- **It is noise anyway.** Roughly two flops per nonzero, next to a full assembly and a direct
  factorization per Newton step. Assembly and the solver dominate; the matvec does not.
- **Nobody ever measured it.** The `splinalg.f90:94-99` comment measures *atomics versus
  reduction* — two threaded variants — never threaded versus serial. There is no benchmark,
  devlog or comment in the tree establishing the threading ever paid.

**What landed.** A new `BELFEM_OMP` define, distinct from the existing `OMP`:

| define | meaning | set by | state |
|---|---|---|---|
| `OMP` | OpenMP exists in this build; here is what the **third parties** will do with it | `USE_OPENMP` (`config_gcc.cmake:47`, `config_icc.cmake:76`) | **ON** |
| `BELFEM_OMP` | **BELFEM's own** kernels may thread | `USE_OPENMP AND USE_BELFEM_OPENMP` (`CMakeLists.txt:179-181`) | **OFF** |

All 12 `!$omp` blocks in `splinalg.f90`, `arpacktools.f90` and `parpacktools.f90` are wrapped in
`#ifdef BELFEM_OMP`, and the four `use omp_lib` guards in the ARPACK pair moved to it. Verified by
preprocessor: `cpp -P -DBELFEM_OMP splinalg.f90` yields 11 `!$omp` lines, `cpp -P` alone yields 0.

**Nine `#ifdef OMP` sites deliberately did NOT move**, and three of them carry an explicit
"NOT a mistake" comment in the source, because moving them would break third-party behaviour on
every default build:

| site | what it does | consequence if moved |
|---|---|---|
| `cl_SolverWrapper.cpp` `hatch_turtle()` | the oversubscription warning | the warning silently disappears — and third-party threading is exactly what we cannot influence |
| `pardisotools.f90` `gParameters(3)` | sets **PARDISO's** thread count | PARDISO silently serialised |
| `banner.cpp` `Threads Used` | reports the third-party budget | user loses the number |
| `banner.cpp:21,438`, `cl_SolverWrapper.cpp:16,284`, `armadillo.hpp:26` | `<omp.h>`, helpers, Armadillo's own OpenMP | build/feature loss |

**Two things this does NOT do**, both documented in `doc/parallel_execution.md`:

1. It does not make the process single-threaded. STRUMPACK, MKL and threaded BLAS still spawn
   workers on the same 2 MiB Darwin stacks, so **`OMP_STACKSIZE` remains a live knob** and a fault
   inside a `gomp_thread_start` frame is not automatically a BELFEM bug.
2. It does not repair the kernel. **Turning `BELFEM_OMP=ON` re-arms DR-155**, since the private
   copy is still on the worker stack. The source comment at the reduction says so.

**Side benefit:** with its own reduction gone, `matvec_csc` is now **bitwise deterministic across
thread counts**, removing one contributor to the convergence-path noise documented at
`doc/parallel_execution.md` ("Measure phases, not steps").

---

## 2a. Superseded: Why the Heap Buffer Was the Alternative Spine

Keep the algorithm — the rejection of per-nonzero atomics in the `splinalg.f90:94-99` comment
still stands, and nothing measured here contradicts it. Change only *where the per-thread copy
lives* — **and restore, by hand, the initialisation the reduction was providing for free.**

A Fortran `allocatable` array is heap-allocated regardless of thread. A per-thread column of a
`(n, nthreads)` buffer gives:

- the same one-combine-at-the-end cost model the reduction had,
- an allocation that fails loudly instead of hitting a guard page,
- a memory cost that is now explicit rather than implied.

`ybuf` is column-major, so **each thread's column is contiguous** — right for the scatter. The
combine must be written in the same direction: `y = y + ybuf(:,t)` over `t`, never an outer loop
over rows reading `ybuf(i,1), ybuf(i,2), …` at stride `n`, which is the cache-thrash inversion
`doc/coding_philosophy.md` forbids (D8).

Rejected alternatives:

- **Raise `OMP_STACKSIZE` from inside `main`** (`setenv` before the first parallel region). Works,
  but it is a global side effect set by a library on behalf of its host, it must run before any
  OpenMP use anywhere in the process, and it leaves the real coupling in place — the next larger
  problem re-crashes at whatever number was hard-coded.
- **Per-nonzero atomics.** Already rejected on measured grounds by the original author
  (`splinalg.f90:95-97`). Do not reopen without a measurement.
- **The transposed-CSR duality the kernel comment names itself** (`splinalg.f90:98-99`). Evaluated
  and rejected: a CSC-structured product is a scatter under any labelling, so the duality moves the
  scatter to the other kernel rather than removing it — `cl_SpMatrix.cpp:1678-1681` already selects
  between the kernels on exactly that duality. It is also unavailable in practice, because
  `fn_matrix_type.hpp:25-36` forces CSC for MUMPS, UMFPACK and SuperLU. Recorded here so nobody
  re-derives it or concludes "just store CSR".
- **Row-blocking** — run the column loop once per row window, bounding the scratch to
  `B * 8 * nthreads` with no allocation at all. Rejected on the per-block re-scan of the whole
  column structure (`nnz` reads × `nblocks`), but recorded as the only option that removes the
  allocation question entirely.
- **`-fmax-stack-var-size` / `-fstack-arrays`.** Not a fix: those govern *local* array placement,
  not OpenMP privatisation, and a private copy cannot live in static storage. Stays rejected.

## 3. Gap Table

| # | State | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| 1 | per-thread private `y` | collision-free scatter | yes, but on the stack | (c) explicit | `splinalg.f90:103` |
| 2 | **identity-initialisation of each private copy** | correctness | yes — *implicitly*, by `reduction(+:)` | (c) explicit | **D1.** `allocate` does not zero |
| 3 | team size actually entering the region | **combining** `ybuf` (sizing uses the cap — D14) | no — never queried | (c) explicit | capture `omp_get_num_threads()` *inside* the region; `omp_get_max_threads()` sizes (D1, D14) |
| 4 | allocation failure path | large n × many threads | no — stack overflow is silent | (c) explicit | **D5**: no abort primitive exists in `src/sparse/*.f90` |
| 5 | combine step | correctness | implicit in `reduction` | (c) explicit | must reproduce **reduction semantics**, not just the `+` |
| 6 | combine direction | performance | n/a | (c) explicit | **D8**: contiguous `y + ybuf(:,t)`, not strided over rows |
| 7 | serial (1-thread) path | avoid n×1 scratch for nothing | no | (a) rebuildable | branch on the team size and write `y` directly |
| 8 | `USE_OPENMP=OFF` build | compiling at all | yes today (`!$omp` are comments) | (c) explicit | **D4**: `#ifdef OMP` + `use omp_lib`, per `arpacktools.f90:13-15` |
| 9 | MKL builds | — | ~~not affected~~ **affected** | (c) explicit | **D2**: `cl_SpMatrix.cpp:1737-1774` has no MKL branch |
| 10 | `matvec_csr` | — | not affected | (a) | gathers into a private scalar, no array reduction |
| 11 | in-suite regression at realistic `n` | catching a wrong combine | fixture exists at `tN = 64` | (c) explicit | **D7**: `tests/sparse/test_SpMatrix.cpp:1782` |
| 12 | documentation of BELFEM's own OpenMP | truthfulness | no — two docs deny it exists | (c) explicit | **D9** |

### 3.1 Cross-cutting findings

- **The memory cost is real but does not invert the tuning advice.** `ybuf` is
  `n * 8 * nthreads` bytes — at 16 threads and 1M dofs, 128 MiB per rank *per matrix*; at the
  4M × 16 case R4 mandates, 512 MiB. `doc/parallel_execution.md:19-23` already says give ranks, not
  threads, and `:119-126` already records that oversubscription costs memory; heap scratch
  **reinforces** that advice rather than inverting it. What must change is the *reason* given, and
  the two false sentences in D9.
- **Both call sites pass caller-owned `aY.data()`** (`cl_SpMatrix.cpp:1705-1715` and
  `:1747-1757`), and the kernel zeroes `y` itself (`splinalg.f90:92`). The transposed path reaches
  the same kernel with **swapped dimensions** via the CSR↔CSC duality (`:1678-1681`), so `n` at
  kernel entry is `mNumCols` when `aTransposedFlag` is set. **Size `ybuf` from the kernel's own
  `n`, never from a matrix member.**
- **`mSwap` is NOT the precedent it looked like** (D6). It exists only under `BELFEM_NETLIB`
  (`cl_SpMatrix.hpp:85-88`, `cl_SpMatrix.cpp:1927-1932`), so it is absent in exactly the builds
  D2 just brought into scope; it is one vector of `max(rows,cols)` with no thread multiplier
  (`:714-715`); and it is **live across the kernel call** whenever `aBeta != 0`, which is the
  residual path that crashed (`:1673` → `:1729`, driven from
  `cl_FEM_DofMgr_SolverData.cpp:2368,2386`). Reusing or aliasing it is a silent overwrite.

## 4. Defects Found (jury round, 2026-09-01)

Severities are the reconciliation table's, in `tmp/ai_exchange/review_matvec_csc_omp_stack_overflow.md`.

- **D1 — CRITICAL — the draft R1 produces silently wrong matvecs.** Raised by Grok, confirmed by
  Claude. `allocate( ybuf(n,nthreads) )` does not zero. The current kernel is correct because of
  *two* zeroings (`splinalg.f90:92` for the caller's buffer, and `reduction(+:)`'s identity for
  each private copy); the draft R1 removed the second and did not replace it, so the combine sums
  heap garbage. Second route to the same hole: sizing from `omp_get_max_threads()` — the cap, not
  the team the region receives — leaves columns of threads that never ran uninitialised, and
  "combine across columns" adds them. This is exactly the failure the plan's original audit trail
  named as worse than the crash. (*rev. 2:* the remedy is in the zeroing and the combine, not the
  sizing — sizing from the cap is correct and unavoidable; see D14 for the two closed designs.)
- **D2 — HIGH — the MKL scope guard is false.** Raised independently by Christian, Claude, Codex
  and Grok. `SpMatrix::multiply( const Vector<real> &, Vector<real> & )` (`cl_SpMatrix.cpp:1737-1774`)
  has no `#ifdef BELFEM_MKL`: it switches on `mType` and calls `matvec_csc` (`:1747`) or
  `matvec_csr` (`:1761`) directly, while the 5-argument sibling routes through
  `mkl_sparse_multiply` (`:1553`). Live callers: `cl_FEM_DofMgr_EigenValues.cpp:1670, 1706, 1712,
  1732-1733, 1755-1759` (inverse iteration, per-iteration) and `operator*` at `cl_SpMatrix.hpp:641-648`.
  Three consequences: MKL builds get no MKL where the matvec is innermost; DR-155's crash is live in
  MKL builds; and two implementations of one operation have already drifted — the 5-arg overload
  validates dimensions with two `BELFEM_ASSERT`s (`:1633`, `:1638`), the 2-arg one validates nothing.
- **D3 — HIGH — the D2 repair has a trap.** Forwarding the 2-arg overload to the 5-arg one as
  `this->multiply( aX, aY, 1.0, 0.0, false )` is instruction-for-instruction equivalent on the
  Fortran side (`mSwap` save at `:1673` skipped, `aY *= aAlpha` at `:1725` skipped, `aY += mSwap`
  at `:1729` skipped; `:1678`/`:1683` reduce to the 2-arg switch; `:1680-1681` give
  `&mNumRows`/`&mNumCols`), and gains the two missing asserts. But the Fortran kernel zeroes `y`
  itself, so 2-arg callers are entitled to pass an **uninitialised** `aY` — and the eigen driver
  does: `mZreal.set_size( mN )` at `cl_FEM_DofMgr_EigenValues.cpp:1661` takes no fill value and
  `:1706` is its first use. `mkl_sparse_d_mv` (`:1600-1603`) computes `y := αop(A)x + βy` and the
  Inspector-Executor sparse BLAS is not documented to short-circuit `β == 0`; `0.0 * garbage` is 0
  for finite garbage and **NaN for a NaN/Inf bit pattern**. `operator*` is already safe
  (`cl_SpMatrix.hpp:646` constructs `Vector<real> aY( aA.n_rows(), 0.0 )`); the eigen sites are not.
- **D4 — HIGH — the draft R1 breaks the `USE_OPENMP=OFF` build.** Raised by Codex and Grok.
  `USE_OPENMP` is a real option (`CMakeLists.txt:78`); `-fopenmp` and the `OMP` definition are
  added only under it (`config/compiler/config_gcc.cmake:38-47`, `list( APPEND BELFEM_DEFS "OMP" )`
  at `:47`). Today `splinalg.f90` has no language-level OpenMP dependency. The in-tree pattern to
  copy is `arpacktools.f90:13-15` — `#ifdef OMP` / `use omp_lib` / `#endif` (also
  `arpacktools.f90:330`, `pardisotools.f90:16`).
- **D5 — HIGH — RESOLVED 2026-09-01 by O1 (no Fortran allocation remains, so the C++ `BELFEM_ERROR` covers it). Kept for the record: R3 could not be deferred, and named a primitive that does not exist.**
  `grep -rniE "error stop|BELFEM_ERROR|\bstop\b" src/sparse/*.f90` returns nothing. There is no
  Fortran abort helper in the tree, so "`stat=` + `BELFEM_ERROR` equivalent" has to be invented
  (`error stop`, or a `bind(C)` abort wrapper) **before** R1 ships — shipping `stat=` with no
  branch is another garbage-buffer path. `doc/coding_philosophy.md:603-606` makes failed
  allocation an always-active `BELFEM_ERROR`, not an assert.
- **D6 — HIGH — `mSwap` is a bad precedent twice over.** See §3.1.
- **D7 — HIGH — the gate the plan needs already has a home.** `tests/sparse/test_SpMatrix.cpp:1782`,
  `TEST( SpMatrix, CscMultiplyMatchesReferenceAcrossThreadCounts )`, was written for this kernel
  (its comment: *"the per-nonzero atomic was replaced by an array reduction, so the result must be
  thread-count independent"*). It builds a colliding CSC matrix (`:1789-1795`), compares against a
  dense reference (`:1806-1813`) and loops `{1, 2, 8}` threads (`:1816-1822`) — but at `tN = 64`
  (`:1786`), three orders of magnitude below the cliff. The dense oracle is O(n²) and cannot follow
  `tN` upward, so the large case needs a **sparse** oracle as a sibling test, not a bigger `tN`
  on this one. Meanwhile §7's out-of-repo driver leaves nothing behind.
- **D8 — HIGH — the draft's combine is the strided direction.** "Combine across columns in a second
  parallel loop over rows" reads `ybuf(i,t)` at stride `n`. Contiguous is
  `y = ybuf(:,1); do t = 2, nt; y = y + ybuf(:,t); end do` (or a tree). R6 could fail on the
  combine alone, independent of allocation cost.
- **D9 — HIGH — two documents assert this defect is impossible.**
  `doc/coding_philosophy.md:674`: *"compiled BELFEM code contains **zero** OpenMP pragmas and spawns
  no threads; the only OpenMP touchpoints are read-only queries."* False — 19 `!$omp` directive
  lines across `splinalg.f90` (11), `arpacktools.f90` (4), `parpacktools.f90` (4), and the DR-155
  crash frame is `matvec_csc._omp_fn.0` on `gomp_thread_start`. `CLAUDE.md` already carries the
  correct version. `doc/parallel_execution.md:20-22`: *"`OMP_NUM_THREADS` accelerates only the
  third-party numerics (STRUMPACK, MKL)"* — false for the same reason.
  `scripts/check_doc_claims.py` passes 37/37 over this tree, so this class of claim is outside its
  coverage; worth a checker item in its own right.
- **D10 — MEDIUM — a thread cap must bind the team, not only the buffer.** Raised by Codex. If O3
  caps the count used to size `ybuf` but the region still takes the runtime's team,
  `omp_get_thread_num()+1` indexes past the buffer. Any cap must appear as `num_threads(...)` on
  the construct.
- **D11 — MEDIUM — the reproducer recipe (§7) is under-specified and partly wrong.** `pointers` is
  `dimension( m + 1 )` (`splinalg.f90:75`) and is read at `:107` as `pointers(j+1)`, so
  `pointers(j) = j` for `j = 1..m` leaves the last entry unset. `values` is never assigned, so the
  stated oracle *"the correct result is `y == x`"* does not follow. And `cl_SpMatrix.cpp:1667` is a
  comment: both sites pass `this->indexing_base()` (`:1668`, `:1741`), so base 0 is the C
  convention, not "what both call sites pass". Fixed in §7 below.
- **D12 — MEDIUM — O1's stated rationale is wrong.** The plan argued for the C++ member on
  "called once per Krylov iteration". Refuted by Grok and confirmed: the crashing path is one
  residual `multiply` per Newton step (`cl_FEM_DofMgr_SolverData.cpp:2368`, `:2386`); direct
  solvers do not call this kernel inside the linear solve. The per-iteration caller is the eigen
  path via D2's overload. O1 must be re-argued on **lifetime and the O3 cap**, not call frequency.
- **D13 — FALSE POSITIVE (retracted 2026-09-01) — "inner `i` is not private → data race."**
  Raised by Grok; refuted on two independent grounds. (a) In Fortran, the iteration variable of a
  *sequential* loop inside a parallel construct is predetermined private in the innermost enclosing
  construct; only `j` is the worksharing iterator, so `i` at `splinalg.f90:109` is already private.
  (b) `tests/sparse/test_SpMatrix.cpp:1782-1830` runs a colliding CSC product at 2 and 8 threads
  against a dense reference; a shared scatter index could not survive it (last recorded green:
  `check-fast` 2026-08-30). Recorded so it is not re-audited. **The forward-looking half stands:**
  R1 must add `tid` to the private list, and anything that stops being a worksharing iterator when
  the nest is restructured must be privatised explicitly.
- **D14 — HIGH (rev. 2) — rev. 1's own R1 asked for something unimplementable.** In overcorrecting
  D1, the revised R1 said "size from the team size that actually enters the region, not
  `omp_get_max_threads()`" — but the team size cannot be known before the region is entered
  (`omp_get_num_threads()` returns 1 outside), and the allocation must happen before it. **Sizing
  from the cap is fine and necessary; D1's defect lives in the zeroing and the combine, never in the
  sizing.** Two closed designs, one of which R0 must pick by name:
  **(a)** allocate `max_threads` columns; *inside* the region each thread zeroes **its own column**
  before scattering; capture `nt = omp_get_num_threads()` inside (single/master into a shared);
  combine columns `1..nt` only. **(b)** allocate and zero **all** columns before the region; combine
  all columns — unused ones add zeros, costing bandwidth but never correctness.
  Spec subtlety that belongs in R1 either way: `omp_get_max_threads()` bounds the team only when
  the construct carries no larger `num_threads` clause — we write the construct, so there is none;
  and if O3's cap lands, `num_threads(min(cap, max))` on the construct is what makes the bound real
  (D10 restated from the sizing side). Note (a) is also the NUMA-right choice: the thread that
  zeroes a column is the thread that first-touches its pages; (b) puts every page on the master's
  node. R6 runs on a UMA Mac and structurally cannot see that difference — logged in R6.
- **D15 — MEDIUM (rev. 4, found during the R5 run) — "BELFEM's compiled code has no OpenMP" is
  wrong a SECOND time, and D9's fix did not go deep enough.** The matrix backend is header-only and
  compiles its **own** OpenMP into BELFEM translation units: Armadillo because BELFEM asks it to
  (`src/linalg/armadillo/armadillo.hpp`, `#define ARMA_USE_OPENMP` under `#ifdef OMP`), Blaze on its
  own (`blaze/system/SMP.h`, OpenMP mode whenever `_OPENMP` is defined; `blaze/config/SMP.h` defaults
  shared-memory parallelisation on). Measured on the sandbox Darwin release build:
  **129 of 339 objects reference `GOMP_*`.** Platform-independent — Linux defaults to Armadillo,
  Apple to Blaze. So gating `BELFEM_OMP` off does **not** make BELFEM's objects OpenMP-free, and a
  `gomp_thread_start` frame with BELFEM code above it does not imply a BELFEM pragma. Corrected in
  `doc/coding_philosophy.md` and `doc/parallel_execution.md` (both of which this session had already
  corrected once, for D9 — the first correction was true but incomplete).
  **Scope limit, checked on Christian's challenge and it is the important half:** this is NOT a
  second route to the DR-155 failure. That mechanism needs an OpenMP *array* reduction, and
  **neither backend uses a reduction of any kind** — `grep -rn "reduction *("` returns **zero** hits
  in both `blaze/` and `armadillo_bits/`. Blaze partitions the destination so each thread writes its
  own slice of shared memory (`parallel shared`, `for schedule`, `sections`); Armadillo runs
  element-wise loops with `critical`/`atomic` (58 `parallel for schedule`, 12 `critical`, 5
  `atomic`). So gating `BELFEM_OMP` off genuinely removes this failure class from BELFEM's own code.
  The first draft of this finding implied `OMP_STACKSIZE` mattered *because of* the backend; it does
  not, and the docs were corrected again to say so. The residual truth is narrower and only about
  reading stack traces.
- **D16 — the "nine sites stay on `OMP`" decision is EXECUTED-verified, not just argued.**
  `hatch_turtle()` fired in the R5 run on this machine: `OMP_NUM_THREADS` unset, so the runtime
  returned the affinity-mask popcount — *"threads requested by this rank : 16, PHYSICAL cores
  available here : 8"*, with the banner's own line *"the runtime defaulted to the affinity mask,
  which counts every hyperthread. It is not a recommendation."* Under the blanket
  `OMP` -> `BELFEM_OMP` sweep first proposed, that warning would have vanished on every default
  build. This is the blind-default case the O3 discussion predicted, observed live.

## 5. Ordered Steps

**Rev. 4:** R1-R6, R8 and R9 are struck — they implemented the heap buffer, which was not the
chosen resolution. They are kept, struck through, because they are the record of the design that
was rejected and of the defects (D1-D14) that made it expensive.

- [x] **R0** — **DONE 2026-09-01 (Christian).** Superseded in outcome by R11: the design fork
  closed on "remove the parallelism" rather than on where the buffer lives.
- [ ] ~~**R1** — Replace the array reduction with a heap `ybuf(n, nthreads)` reproducing reduction
  semantics.~~ **STRUCK (rev. 4)** — not implemented; the directives are gated off instead.
- [ ] ~~**R2** — Serial short-circuit.~~ **STRUCK** — the kernel is now serial unconditionally.
- [ ] ~~**R3** — Allocation failure path.~~ Already folded into R1 by D5, and struck with it.
- [ ] ~~**R4** — Re-run the §1.1 sweep at n up to 4M x 16 threads.~~ **STRUCK** — with the
  directives gated off there is no per-thread copy to size. Superseded by R13.
- [ ] ~~**R6** — Benchmark the heap buffer against the current kernel.~~ **STRUCK.** The
  measurement that matters now is R14 (threaded vs. serial), which is the one nobody ever ran.
- [x] **R8** — **DONE 2026-09-01 (Christian: "`SpMatrix::multiply` should use MKL if MKL is
  available").** The 2-argument overload (`cl_SpMatrix.cpp`) no longer switches on `mType` and calls
  the Fortran kernels; it forwards to the 5-argument one as
  `this->multiply( aX, aY, 1.0, 0.0, false )`, which carries the `BELFEM_MKL` branch. Three
  consequences: MKL builds now use `mkl_sparse_d_mv` on the ARPACK inverse-iteration path; the
  overload gains the two dimension `BELFEM_ASSERT`s it never had (D2's drift); and **the only two
  remaining `matvec_csc`/`matvec_csr` call sites are now inside the `#else` of `BELFEM_MKL`**, so an
  MKL build never reaches the Fortran kernels — the plan's originally-*false* scope guard is true as
  an end state.
  The forward is instruction-for-instruction equivalent on the Fortran branch: `aBeta == 0.0` skips
  both the `mSwap` save and `aY += mSwap`, `aAlpha == 1.0` skips `aY *= aAlpha`, and
  `aTransposedFlag == false` reduces `tRunCsr`/`tIndices`/`tNumRows`/`tNumCols` to exactly the old
  `switch`.
  **D3's fill landed as `#ifdef BELFEM_MKL` only, and that asymmetry is deliberate**: the Fortran
  kernels zero `y` themselves, so callers of this overload may pass an uninitialised `aY` and the
  eigen driver does (`set_size()` with no fill value, then straight into here). `mkl_sparse_d_mv` is
  not documented to short-circuit `beta == 0`, and `0.0 * NaN` is NaN. On the Fortran branch the
  kernel's own zeroing already covers it, so an unconditional fill would be a second O(n) pass on a
  hot path. The reasoning is in the source comment, not only here.
  **Gates:** `g++ -fsyntax-only` green on `cl_SpMatrix.cpp` (netlib/Blaze), and `aY.fill( 0.0 )`
  compile-probed on `Vector<real>`. **The MKL branch itself is UNVERIFIABLE on this machine** —
  `config/linalg/config_mkl.cmake:2-3` makes MKL a hard configure error on macOS, so the one path
  this change exists to enable is the one that cannot be compiled here. It needs a Linux MKL tree;
  folded into R13.
- [ ] ~~**R9** — Add a large-`n` collision test with a sparse oracle.~~ **STRUCK as specified.**
  The cliff is unreachable with the directives off, so the large-`n` case has nothing to catch.
  `CscMultiplyMatchesReferenceAcrossThreadCounts` (`tests/sparse/test_SpMatrix.cpp:1782`) **stays**
  — it now passes trivially, and it is the guard that fires the day anyone sets
  `USE_BELFEM_OPENMP=ON`.

**What actually landed and what it still owes:**

- [x] **R11** *(new, rev. 4)* — Introduce `BELFEM_OMP`, distinct from `OMP`. Option
  `USE_BELFEM_OPENMP` (default OFF, `CMakeLists.txt:88`), define appended under
  `USE_OPENMP AND USE_BELFEM_OPENMP` (`:179-181`). All 12 `!$omp` blocks in `splinalg.f90`,
  `arpacktools.f90`, `parpacktools.f90` wrapped in `#ifdef BELFEM_OMP`; the four `use omp_lib`
  guards in the ARPACK pair moved from `OMP`. `splinalg.f90`'s two subroutine headers gained the
  same guarded `use omp_lib` (Christian) — not needed by the directives alone, but the file's
  OpenMP dependency is now explicit and any future `omp_get_*` in the kernel compiles rather than
  failing under `implicit none`. Nine third-party `#ifdef OMP` sites left alone, with
  "NOT a mistake" comments at `hatch_turtle()`, `pardisotools.f90` and `banner.cpp`.
  **Executed gates:** preprocessor — 11 `!$omp` lines with the define, 0 without; and
  `gfortran -fsyntax-only` **green on all three kernels in BOTH token states** (6/6; gcc 16.1.0,
  `parpacktools` needs `-I` at the MPI include). That is a syntax gate, not a build.
- [x] **R12** *(new, rev. 4)* — Documentation. `doc/parallel_execution.md` gains a
  "BELFEM's own OpenMP, and `OMP_STACKSIZE`" section — why OFF, the 262,144 arithmetic, the
  re-arming warning, and that `OMP_STACKSIZE` still matters for third-party workers.
  D9's two false claims corrected (`coding_philosophy.md`, `parallel_execution.md` rule 1) and the
  reassociation sentence re-attributed. `CLAUDE.md`'s build-options bullet updated.
  A `DR-155` note added at the reduction itself (`splinalg.f90`).
- [◐] **R13** *(2 of 3 configurations DONE 2026-09-01)* — **Build gate.** Configure and build clean, then
  `make check`. Three configurations now: `-DUSE_BELFEM_OPENMP=OFF` (default) **and**
  `-DUSE_BELFEM_OPENMP=ON`, the latter to prove the token is a real switch and not a one-way
  deletion; **plus a `-DUSE_MKL=ON` build on Linux for R8**, since the MKL branch cannot be compiled
  on Darwin at all and is the whole point of that change. On the MKL build, an eigenvalue case is
  the one that exercises it — `mkl_sparse_d_mv` reached through the 2-argument overload, with the
  uninitialised-`aY` fill in play.

  **Configuration 1 DONE — default (`USE_BELFEM_OPENMP=OFF`, `USE_DEBUG=ON`, `USE_MKL=OFF`,
  Darwin/Blaze, `cmake-build-debug`):** clean build, then `make check` **15/15, 0 failed**
  (fast 33.03 s over 8 tests, mpi 10.12 s over 2, total 44.81 s).
  **The fix is confirmed PRESENT IN THE BINARY, not inferred from the suite result.**
  `nm` on the built `splinalg.f90.o` (Sep 1 03:40; `lib/libbelfem.a` 03:46) shows **zero** OpenMP
  symbols. Falsification control on the same source: compiled `-DBELFEM_OMP -fopenmp` it yields
  `_matvec_csc._omp_fn.0`, `_matvec_csr._omp_fn.0`, `_GOMP_parallel` and `_GOMP_atomic_*`;
  compiled without the token, zero. `matvec_csc._omp_fn.0` is the literal faulting frame of the
  §1.2 crash report, so its absence is caused by the gating and by nothing else. The control also
  establishes most of what configuration 2 was for at the object level: the source still compiles
  and still outlines correctly with the token ON, so the token is a switch and not a one-way
  deletion.
  **What this green does NOT establish, and the honesty matters here:**
  (a) `CscMultiplyMatchesReferenceAcrossThreadCounts` (`tests/sparse/test_SpMatrix.cpp:1782`)
  passed, but with the kernel serial its `{1, 2, 8}` sweep is **trivially** identical — the suite
  does not exercise the thing that broke, and cannot until the token is ON;
  (b) this was a **debug** build (`USE_DEBUG=ON`), so release is unbuilt;
  (c) the MKL forwarder (R8) is untouched by it and remains unbuildable on Darwin.

  **Configuration 2 DONE — Release + shared (`USE_DEBUG=OFF`, `USE_SHARED_LIBS=ON`,
  `USE_BELFEM_OPENMP=OFF`, `USE_MKL=OFF`, Blaze), sandbox tree `/tmp/build`:** clean configure and
  build, 0 errors, `make install` to `/tmp/belfem-dr131` green. This closes limit (b) above — the
  change is now built in **both** `USE_DEBUG` states, and the shared build is the one that exercises
  the Darwin `dlopen` / `dynamic_lookup` path. Gating re-confirmed in this binary too: `nm` on its
  `splinalg.f90.o` shows zero OpenMP symbols. The deck then ran on it (R5).
  **Configuration 3 — Linux `-DUSE_MKL=ON` — remains owed**, and it is the only one that can compile
  R8's forwarded MKL branch at all. On the ON build the pre-existing DR-155
  crash returns above 262,144 rows unless `OMP_STACKSIZE` is raised — that is expected, and is
  itself the discrimination that the gate is testing the right thing.
- [x] **R5** — **RAN GREEN 2026-09-01.** `examples/tape_quench_usermat`, `gmsh -3 tape.geo`
  (21,303 nodes / 146,265 elements), a **Release + `USE_SHARED_LIBS=ON`** build from the current
  source installed to `/tmp/belfem-dr131`, plugins built against that prefix, **default environment
  with `OMP_STACKSIZE` unset** (verified absent from the environment before launch).
  **The first BDF1 step — the exact step that died — completed and converged:**
  `BDF1 1 | t : 0.5000 ms`, `timestep succeeded | magnetic : 1.033e-10 | thermal : 6.684e-12`,
  58 s. That is this gate's bar, and it passed.
  It then ran on far past the bar: **54 accepted timesteps, 6 rejected, t = 10.98 ms, 6 h 50 m
  wall, ZERO crash signatures** (`SIGBUS|SIGSEGV|Segmentation|Abort trap|Killed|bad_alloc` all
  absent from 2003 log lines). Stopped deliberately, not by failure — see the note below.
  Evidence log: `/tmp/build/tape_quench_r5_evidence.log`.
  **Why it was stopped rather than completed:** the deck is `simulation time : 0.3 s` and reached
  10.98 ms in 6 h 50 m, with per-step cost 5-13 min through the stiff patch. Completion projects to
  **~114 hours** — not viable, and not what R5 asks for. Deck completion belongs to DR-131, which is
  being closed on a cheaper in-scope deck (`disk_pulse`).
  **A live by-product worth recording:** `hatch_turtle()` fired at startup on the blind default —
  *"threads requested by this rank : 16, PHYSICAL cores available here : 8"* with its own line
  *"the runtime defaulted to the affinity mask, which counts every hyperthread"*. See D16.
- [ ] **R10** *(optional now, was O2's gate)* — `OMP_DISPLAY_ENV=TRUE`, or the worker-stack probe,
  on Darwin and on the CI host. **No longer decides severity** — the defect is unreachable on a
  default build — but it is the number that would tell a future implementer whether re-arming
  `BELFEM_OMP` is safe anywhere. Cheap; worth doing before anyone turns the token on.
- [ ] **R14** *(new, rev. 4 — the measurement nobody ever ran)* — Before `USE_BELFEM_OPENMP` is
  ever defaulted back ON, measure threaded vs. serial: same deck, same `OMP_NUM_THREADS`, both
  builds, `OMP_STACKSIZE` raised on the ON side, compared on **per-phase** timings rather than step
  wall-time (`doc/parallel_execution.md`, "Measure phases, not steps"). ARPACK's kernels deserve
  their own number: they sit in the Arnoldi inner loop and have a much better prior for mattering
  than `matvec_csc` did.

## 6. Open Design Questions

- **O1 — Where does the buffer live? — SUPERSEDED (rev. 4): there is no buffer.** The kernel no
  longer threads, so no per-thread scratch exists, and the `CommunicationMutable` /
  `rbuffer()` / `ibuffer()` prototype that O1 had settled on was **reverted the same day** — with
  the kernel serial and the `gesv` pivot use ruled noise before the release, it had no consumer
  left. Everything below is kept as the reasoning that would apply again if `BELFEM_OMP` is ever
  turned on, **not as a description of the current tree**. Was: **a process-level scratch owned by
  the communicator.** Neither of the two shapes the plan first
  posed. Fortran-side `allocatable` allocates on every call and leaves D5 to be invented; a
  per-`SpMatrix` member is resident for the object's lifetime and multiplies by the number of live
  matrices (`mSystemMatrix` + `mDirichletMatrix` in `cl_FEM_DofMgr_SolverData.hpp:101-108`,
  `mK` + `mM` in the eigen driver) — order 512 MiB per rank at 1M × 16 where one buffer costs
  128 MiB. `mSwap` is not the precedent it looked like (D6).
  The chosen shape is `CommunicationMutable`'s `rbuffer()` / `ibuffer()`: one buffer per rank, an
  owner whose lifetime is the process, and allocation failure that is an ordinary `BELFEM_ERROR` on
  the C++ side — **which dissolves D5 entirely.** Access is by call, **not by inheritance**:
  `SpMatrix` needs `gComm.comm_mutable()->rbuffer( ... )`, not a `CommunicationObject` base.
  **Lifetime contract (Christian, 2026-09-01):** the communicator is always closed last — every
  `main()` ends with `gComm->finalize()` — so an `SpMatrix` still working after that is a caller
  contract violation, not a case the buffer defends against. No dead-communicator guard.
  **Two roles, two mechanisms** (they differ in how their demand varies):
  - the **swap** role (`mSwap`'s job) is thread-independent and long-lived → register the demand.
    `allocate_swap()` (`cl_SpMatrix.cpp:712-715`) is already the single choke point and is already
    called from all five sites where dimensions settle (`:105`, `:203`, `:334`, `:1470`, `:1855`),
    and it already computes the right quantity — `max( mNumRows, mNumCols )`, which is what the
    transposed path needs since the kernel enters with `n = mNumCols` (`:1680`) and
    `mDirichletMatrix` is rectangular. Register / deregister with a rescan to the high-water mark
    over live members, so the buffer **shrinks** when a large matrix dies.
  - the **`ybuf`** role is transient and depends on `nthreads`, which changes at runtime — the
    suite's own `tests/sparse/test_SpMatrix.cpp:1816-1822` calls `omp_set_num_threads` in a
    `{1, 2, 8}` loop between `multiply` calls on an already-constructed matrix, so a demand
    registered at construction is stale *low* on the second iteration. → plain grow-on-demand
    `rbuffer( n * nt )` at the call site, never registered.
  **Follow-on work this creates, tracked but out of DR-155's kernel scope:** retiring `mSwap` in
  favour of the registered buffer (a clean win — it removes the `BELFEM_NETLIB` asymmetry that D2
  made load-bearing), which needs (i) two *named* regions rather than one buffer plus an offset
  convention, because `mSwap` is live **across** the kernel call whenever beta != 0 (`:1673` →
  `:1705` → `:1729`) so the two roles coexist within a single call; (ii) re-registration in
  `operator=( SpMatrix && )` (`:1873`), which *steals* `mSwap`/`mSwapSize` (`:1929-1931`) instead
  of going through `allocate_swap()` as copy-assign does at `:1855`; and (iii) shrink hysteresis, so
  a build/destroy loop over temporaries does not realloc on every pair. Do it **after** the kernel
  is green, as its own step.
- **O2 — Is the Linux cliff real, and where? — DEMOTED (rev. 4): no longer a severity question.**
  Unreachable on a default build. It still matters as a precondition for ever re-arming
  `BELFEM_OMP`, so R10 is kept as optional rather than struck.
  The plan's original 1,048,576-row prediction assumed libgomp inherits an 8 MiB worker stack. Two
  incompatible readings emerged and neither is decidable from this tree:
  (a) *Grok:* macOS's pthread default for secondary threads is 512 KiB, yet the measured cliff is
  exactly 2 MiB — so libgomp is **setting** 2 MiB, and would on Linux too, putting the Linux cliff
  at the **same 262,144 rows** and making the plan's figure simply wrong.
  (b) *Claude:* the worker stack derives from `RLIMIT_STACK` at process start, which is site policy —
  schedulers and containers set it lower, `ulimit -s unlimited` changes glibc's fallback rather than
  removing the bound — so the cliff is `stack_limit / 8` and can be **below** Darwin's.
  Both agree the 1,048,576 figure is unearned. **R10 settles it for the price of one command.**
  If (a) holds, this stops being a porting note and becomes a live limit on every platform.
- **O3 — Should the kernel cap its own thread count? — MOOT (rev. 4): the kernel has no threads.**
  The sizing table below is kept because it is the memory a future `BELFEM_OMP=ON` build would
  need, and because it is half of why the token is OFF. Was resolved as: **no cap, document the
  cost instead.** The footprint is `max( rows, cols ) * nthreads * 8` bytes **per
  rank**:

  | case | n | threads | scratch |
  |---|---|---|---|
  | typical | 1 M | 20 | **152 MiB** — noise against what the rest of the solve needs |
  | large | 4 M | 16 | 488 MiB |
  | extreme | 10 M | 128 | **~10 GB** — acceptable on a machine that has 10 M dofs and 128 cores |

  **`T` is what the rank *sees*, and it is NOT bounded (corrected 2026-09-01, Christian).**
  No-oversubscription is a strong recommendation, not something BELFEM can enforce — nothing stops
  a user setting `OMP_NUM_THREADS` as they like. Worse, the *blind default* is the common case and
  the plan must size for it: with `OMP_NUM_THREADS` unset the runtime hands back the affinity-mask
  popcount, **every hyperthread included** (`cl_SolverWrapper.cpp:435-437` says exactly this), so a
  128-thread / 64-core node gives `T = 128` without the user having chosen anything — a 2× buffer
  for zero speed on a bandwidth-bound scatter. Therefore:
  - **Size from `omp_get_max_threads()`, re-read on every call.** Undersizing is a buffer overrun
    (D10, D14), so the value the rank sees is the only safe bound — however unreasonable it is.
  - **The allocation failure path stops being theoretical** and becomes the main guard. Its
    `BELFEM_ERROR` must name `n`, `nthreads` and the byte count, so the user can see what asked
    for the memory. Written into R1.
  - **Surface the cost where the user already gets told.** `Wrapper::hatch_turtle()`
    (`cl_SolverWrapper.cpp:390-470`) already fires on both oversubscription and the blind default,
    already computes `omp_get_max_threads()`, `physical_cores_available()` and `gComm.node_size()`,
    and its banner already reads *"This costs time and memory"* — DR-155 gives that a named number.
    Two limits to note when wiring it: it is called only from MUMPS and STRUMPACK init
    (`cl_SolverMUMPS.cpp:191`, `cl_SolverSTRUMPACK.cpp:175`), so other solvers and the eigen path
    never see it.

  One thing that *does* hold: **it does not compound with rank count.** The buffer is per rank and
  `n` is the rank's *local* row count, so a node running `R` ranks × `T` threads lands near
  `n_global * T * 8` regardless of `R` — the local `n` shrinks as ranks grow.

  Logged, not reopened: the hyperthread case is the one place a kernel-side cap would actually pay
  (hyperthreads do not speed a bandwidth-bound scatter), so if the buffer ever does become the
  binding constraint, capping the *team* at physical cores is where to look first — bound via
  `num_threads(...)` per D10, never by sizing alone.

  What R4 must watch at the extreme end is therefore not only memory but the combine: at
  10 M × 128 the D8-correct contiguous form is ~1.3e9 additions per matvec, and it is serial.
  See R6, which now carries that as a named measurement rather than a footnote.
  Note also that the *stack* cliff depends on `n`, not on thread count (except 1 vs. >1); thread
  count only becomes a memory question once the copy is on the heap. Do not mix the two.
  D10 still stands for any future cap: it would have to bind `num_threads(...)` on the construct,
  not only the buffer size.
- **O4 — Where does R8's zeroing live? — MOVED with R8 to its own row (rev. 4).** In the forwarder (`aY.fill(0.0)`, local, makes the
  contract explicit) or in the callers (`set_size( mN, 0.0 )` at
  `cl_FEM_DofMgr_EigenValues.cpp:1659-1661` and siblings). Recommendation: the forwarder — the
  overload's contract is what changes, so fixing callers is whack-a-mole.

## 7. Reproducer

Kept out of the repo; rebuild from this description. **R9 is what makes this durable — this recipe
is a crash detector, not a regression gate.** A driver defining `gComm` and `gLog` (the library
declares them extern — see `config/scripts/Add_BelfemLibrary.cmake:44-58`), linking `libbelfem` and
`-lgomp`, declaring the `fspblas.hpp:38-49` binding, and calling `matvec_csc` with:

- `n = m = nnz`, `base = 0` (the C convention; note both call sites actually pass
  `this->indexing_base()`, `cl_SpMatrix.cpp:1668` and `:1741` — D11),
- `indices(j) = j` for `j = 1..nnz`,
- `pointers(j) = j` for `j = 1..m+1` — **`pointers` is `dimension(m+1)`** (`splinalg.f90:75`) and
  entry `m+1` is read at `:107`; the earlier recipe left it unset (D11),
- `values(j) = 1.0` — **must be set**, or `y == x` is not an oracle and the driver only detects
  crashes (D11),
- `n` on the command line.

Sweep `n` under `OMP_NUM_THREADS=16`. The transition is sharp and lands in (262,000, 262,144].

## 8. Definition-of-Done Checklist

- [ ] Every gap-table row mapped to a step or an open question.
- [x] **O1 and O3 decided** (2026-09-01): communicator-owned buffer, no cap. D5 dissolved with
  it. **Still owed before R1 starts:** the D14 design — (a) or (b) — picked by name.
- [ ] **D1 closed: `ybuf` initialisation and team-size handling written and reviewed** — this is the
  one whose failure mode is silent.
- [ ] R4 passes at 4M rows × 16 threads with a default environment.
- [ ] R5 passes: the deck runs without `OMP_STACKSIZE`.
- [ ] R6 shows no regression below the old cliff, **and** the combine direction was benchmarked
  before R1 was finalised.
- [ ] R8 landed: MKL builds use MKL on the 2-argument path, with the D3 fill.
- [ ] R9 landed: a large-`n` collision case with a sparse oracle is in the suite and green.
- [ ] R10 run: `OMP_STACKSIZE` recorded on Darwin and Linux; O2 answered and DR-155's severity
  re-checked against it.
- [ ] D9's two false sentences corrected (R7).

## 9. Audit Trail

Found and measured in a single session while gating DR-131's plugin fix on Darwin
(`devlog/dl20260901_darwin_plugin_gate.md`); the crash was incidental to that work, and the
isolation, the boundary and both controls are that session's own executed gates.

**Jury round 2026-09-01** — Claude pre-registration frozen, then Codex `gpt-5.6-terra`/high and
Grok `grok-4.6`/high in parallel and blind, followed by a Claude verification pass over every
citation. Record: `tmp/ai_exchange/review_matvec_csc_omp_stack_overflow.md`. Outcome: **both
auditors declined to approve R1 as specified**, and the round produced one CRITICAL (D1) that no
single reviewer had — Grok found it, Claude confirmed it against the source. D2 was raised
independently four times (Christian, Claude, Codex, Grok). One finding was refuted and is recorded
as D13 rather than deleted. Three of Claude's own pre-registered findings were corrected in the
verification pass — one rested on a malformed `grep` — and those corrections are stated in the
exchange file rather than quietly absorbed.

**Second-pass sweep, same day (rev. 2, Claude/Fable).** A cold re-read of rev. 1 against the
tree, scoped to oversights rather than re-litigation. Found: rev. 1's own R1 demanded an
unimplementable sizing (D14 — the team size is unknowable before the region; the fix is two closed
initialisation designs, R0 now picks one by name); R10's env print may be inconclusive and gained a
definitive worker-stack probe as fallback; R1's requirements were marked branch-invariant with the
C++ guard precedent named (`banner.cpp:21-22`, `cl_SolverWrapper.cpp:16-17`); R8 was annotated with
its end state (post-R8 an MKL build never reaches `matvec_csc`, the two overloads being the
kernels' only callers); and R6 gained two logged blind spots (UMA first-touch, serial-combine
cost) plus R1 an integer-kind note. No step was struck and no severity moved.

**Rev. 3, same day — O1 and O3 resolved by Christian in conversation.** The buffer is a
process-level scratch owned by the communicator (`gComm.comm_mutable()->rbuffer()`), reached by
call rather than by inheritance from `CommunicationObject`; no thread cap, the cost documented
instead. Two of the reviewer's objections were overruled and are recorded as such: a
dead-communicator guard on the destructor is unnecessary, because the communicator closing last is
a **contract** — every `main()` ends with `gComm->finalize()`, so an `SpMatrix` working after that
is a caller violation, not a case to defend against; and the residency concern was answered with
numbers rather than a cap (152 MiB typical, ~10 GB at the 10 M x 128 extreme, per rank, not
compounding with rank count). The exchange also produced the two-roles split — register the
thread-independent swap demand, keep `ybuf` grow-on-demand, because `omp_set_num_threads` changes
the thread count at runtime inside the suite's own test. `mSwap` retirement is now tracked follow-on
work under O1 rather than part of the kernel change.

**Rev. 4, same day — Christian resolved it by removing the parallelism.** The turn came from a
question the reviewers had not asked: *what does this threading actually buy?* Checking the call
paths answered it — both matvecs are master-only, so OpenMP was their only parallelism, but they sit
next to a full assembly and a direct factorization and are noise against those. And no measurement
of threaded-vs-serial exists anywhere in the tree: the kernel comment measures atomics against the
reduction, two threaded variants. A defect whose fix needed a buffer registry, a two-role split,
reduction semantics reproduced by hand and a 512 MiB scratch budget was closed instead by a
default-OFF build token.

Two of the reviewer's positions were corrected in the process. A Darwin-only ban was argued against
and dropped in favour of OFF everywhere — O2 being unresolved, a Darwin carve-out could have left
the same cliff live on Linux, and platform-conditional threading means platform-divergent rounding
in a framework that validates against published benchmarks. And the proposed blanket
`OMP` -> `BELFEM_OMP` sweep was narrowed after an inventory of all 13 guard sites showed **nine must
stay** — moving `hatch_turtle()`, PARDISO's `gParameters(3)` or the banner's thread count would have
silently deleted the oversubscription warning, serialised a third-party solver, and hidden the
number, all on the default build. Those three now carry "NOT a mistake" comments naming the reason.

The `gesv` pivot buffer was ruled noise on the same grounds and deferred past the release, which
leaves `rbuffer`/`ibuffer` without a consumer for now. Also asked and answered: a buffer registry
would **not** have prevented a third-party crash — MUMPS, STRUMPACK and MKL allocate their own
workspace on their own threads, and the existing instruments for that are the in-wrapper MUMPS
ICNTL(14) ladder and `hatch_turtle()`, which is the second reason that guard had to stay on `OMP`.

**Evidence standing (`doc/ai_collaboration_protocol.md` §11):** the crash mechanism is
**verified** (reproducer, executed 2026-09-01). Everything in §4 is **reviewed, not verified** —
source trace plus AI agreement, with no executable gate run.
**Rev. 4's change is REVIEWED, NOT VERIFIED, and not even compiled.** The only executed check is a
preprocessor one: `cpp -P -DBELFEM_OMP src/sparse/splinalg.f90` yields 11 `!$omp` lines and
`cpp -P` alone yields 0, which proves the gating works and nothing else. R13 (build + `make check`,
both token states) and R5 (the deck on a default environment) are owed and are what would make this
verified.
