# DR-155's plan went to a jury and came back unapproved — the fix would have been silently wrong

**Date:** 2026-09-01
**Topic:** Three-AI jury round on `todo/matvec_csc_omp_stack_overflow.md`; DR-155 widened to carry
the missing MKL branch on the 2-argument `SpMatrix::multiply`
**Module:** `src/sparse` (plan and register only — no source changed)
**AIs:** Claude (pre-registration + verification), Codex `gpt-5.6-terra`/high, Grok `grok-4.6`/high
**Exchange:** `tmp/ai_exchange/review_matvec_csc_omp_stack_overflow.md`

---

## Why this round happened

DR-155 was filed the same day, from a measured Darwin crash: `matvec_csc`'s `!$omp reduction(+:y)`
puts each worker's private copy of the whole result vector on the thread stack, and the run dies
once `n * 8` crosses libgomp's 2 MiB worker stack. The mechanism was never in dispute in this round
and is not disputed now. What went to the jury was the **plan** — and the plan's own audit trail had
asked for exactly this, on the grounds that *"a wrong reduction here produces silently wrong
matvecs rather than a crash."*

That turned out to be the right thing to be afraid of.

## The finding that justified the round

**D1, CRITICAL, raised by Grok and confirmed against the source: the draft R1 would have produced
silently wrong matvecs.**

The current kernel is correct because of **two** zeroings, and only one is visible in the source.
`y = 0.0d0` (`splinalg.f90:92`) zeroes the caller's buffer. `reduction(+:)` then
identity-initialises each *private* copy to zero before the region and adds those copies back into
the already-zeroed original. R1 replaced the reduction with `allocate( ybuf(n,nthreads) )` and a
hand-written combine — removing the second zeroing and not replacing it. `allocate` does not zero;
the combine would have summed heap garbage into `y`.

The same hole a second way: R1's gap table sized `ybuf` from `omp_get_max_threads()`, which is the
**cap**, not the team size the region actually receives. Columns belonging to threads that never ran
stay uninitialised, and "combine across columns" adds them.

Neither Claude's pre-registration nor Codex's audit had this. It is the whole return on the round.

## The scope guard that was false, and Christian's ruling

The plan and the register both said `USE_MKL=ON` never reaches the kernel. **False.** The
`#ifdef BELFEM_MKL` block (`cl_SpMatrix.cpp:1650` → `:1665` → `:1731`) bounds only the 5-argument
`SpMatrix::multiply`. The 2-argument overload (`:1737-1774`) has no MKL branch at all and calls
`matvec_csc` at `:1747` unconditionally. Its callers are the ARPACK inverse iteration
(`cl_FEM_DofMgr_EigenValues.cpp:1670, 1706, 1712, 1732-1733, 1755-1759`) and `operator*`
(`cl_SpMatrix.hpp:641-648`). So MKL builds crash on the eigen path, and the plan contradicted
itself — its own §3.1 already cited `:1747-1757` as a call site.

Raised independently **four** times: by Christian in conversation, and by all three reviewers.

Christian's ruling: **one DR-155 sweep, not a new row** — same file, same overload pair, and the
false guard is DR-155's own scope error. Extend-don't-branch, whose carve-out is for a genuinely
separate subsystem; this is not one.

Christian's framing of the same site is the more useful one: *if MKL is set, we always want to use
the MKL routines.* Forwarding the 2-argument overload to the 5-argument one is
instruction-for-instruction equivalent on the Fortran side at α=1, β=0, untransposed (checked branch
by branch), gains the two dimension asserts the 2-arg overload lacks, and deletes one of DR-155's
two call sites. **But it needs `aY.fill( 0.0 )` first** (D3): the Fortran kernel zeroes `y` itself,
so 2-arg callers are entitled to pass an uninitialised `aY` — and the eigen driver does
(`mZreal.set_size( mN )`, `cl_FEM_DofMgr_EigenValues.cpp:1661`). `mkl_sparse_d_mv` computes
`y := αop(A)x + βy` and the IE sparse BLAS is not documented to short-circuit β = 0; `0.0 * NaN` is
NaN. `operator*` is already safe — it constructs `Vector<real> aY( aA.n_rows(), 0.0 )`.

## What else the round changed

- **The Linux cliff figure was unearned.** "8 MiB → 1,048,576 rows" was assumed. Two incompatible
  readings now exist and neither is decidable statically: Grok's — macOS's pthread default is
  512 KiB yet the cliff is exactly 2 MiB, so libgomp is *setting* 2 MiB and Linux would sit at the
  **same 262,144 rows**; and Claude's — the worker stack follows `RLIMIT_STACK`, so it is site
  policy and can be *below* Darwin's. Both agree the plan's number is wrong.
  **This decides DR-155's severity and costs one command:** `OMP_DISPLAY_ENV=TRUE` on any
  OpenMP-linked binary prints the effective `OMP_STACKSIZE`. Filed as R10, ahead of everything else.
- **"Measured to the byte" overstated the sweep.** The boundary is bracketed to
  (262,000, 262,144] — a 1152-byte window; no point between was run. That the upper edge coincides
  with 2 MiB / 8 B is still strong evidence for the mechanism. Corrected in the plan, the register
  and `todo/README.md`.
- **R3 cannot be a later step.** `grep -rniE "error stop|BELFEM_ERROR|\bstop\b" src/sparse/*.f90`
  returns nothing — there is no Fortran abort primitive to call, so "`stat=` + `BELFEM_ERROR`
  equivalent" names something that does not exist. Shipping `stat=` with no branch is another
  garbage-buffer path. Folded into R1.
- **`mSwap` is a bad precedent twice over.** `#ifdef BELFEM_NETLIB` only (`cl_SpMatrix.hpp:85-88`),
  so absent in exactly the builds D2 brought into scope; and **live across the kernel call**
  whenever β ≠ 0, which is the residual path that crashed (`cl_SpMatrix.cpp:1673` → `:1729`, driven
  from `cl_FEM_DofMgr_SolverData.cpp:2368, 2386`).
- **O1's rationale was refuted.** The plan argued for a C++ member on "called once per Krylov
  iteration". The crashing path is one residual `multiply` per *Newton* step; MUMPS and STRUMPACK
  are direct and do not call this kernel inside the solve. O1 must be re-argued on lifetime, and
  O3 (thread cap) is now a precondition rather than a logged aside — with Codex's trap attached:
  a cap must bind `num_threads(...)` on the construct, or `omp_get_thread_num()+1` indexes past the
  buffer.
- **The gate the plan wanted already exists.** Both auditors asked for a collision-heavy CSC test
  against a serial reference. `tests/sparse/test_SpMatrix.cpp:1782` is exactly that, written for
  this kernel, threads {1, 2, 8} — at `tN = 64`. The step is not "write a test" but "add a
  large-`n` sibling", and the dense oracle is O(n²) so the large case needs a sparse one.
- **Two documents assert this defect is impossible.** `doc/coding_philosophy.md:674` — *"compiled
  BELFEM code contains **zero** OpenMP pragmas and spawns no threads"* — and
  `doc/parallel_execution.md:20-22` — *"`OMP_NUM_THREADS` accelerates only the third-party
  numerics"*. Both false against 19 `!$omp` lines in three Fortran files and a `gomp_thread_start`
  crash frame. `CLAUDE.md` already carries the correct version. **`check_doc_claims.py` passes
  37/37**, so this claim class is outside its coverage — worth a checker item.

## What was refuted, and what the reviewers got wrong about themselves

**Grok's "inner `i` is not private → data race" is a FALSE POSITIVE**, kept in the plan as D13 so it
is not re-audited. Two independent grounds: in Fortran the iteration variable of a *sequential* loop
inside a parallel construct is predetermined private in the innermost enclosing construct, so `i` at
`splinalg.f90:109` already is; and `tests/sparse/test_SpMatrix.cpp:1782-1830` runs a colliding CSC
product at 2 and 8 threads against a dense reference, which a shared scatter index could not
survive. The forward-looking half stands: R1 must privatise `tid`.

**Three of Claude's own pre-registered findings were corrected in verification.** The worst: P1-2
claimed the tree had zero `omp_lib` / `omp_get_*` uses and that R1 would introduce the first. That
"zero" came from `grep -F` matching `\|` literally. Redone correctly, the tree already has the
pattern R1 needs — `#ifdef OMP` / `use omp_lib` at `arpacktools.f90:13-15`, with `OMP` defined at
`config/compiler/config_gcc.cmake:47`. Codex got this right and Claude did not; the conclusion
(R1 must guard, or `USE_OPENMP=OFF` breaks) survives with better evidence. Also corrected: the
`mSwap` precedent was cited without its `BELFEM_NETLIB` guard, and the "once per Krylov iteration"
framing was inherited from the plan rather than checked.

## Files touched

- `todo/matvec_csc_omp_stack_overflow.md` — rewritten: corrected scope guards (struck, not deleted),
  D1-D13 defect list, revised R0-R10, O1-O4, corrected reproducer recipe, evidence-standing note
- `todo/debt_register.md` — DR-155 row widened; the `[P]` count line amended (count unchanged, 8)
- `todo/README.md` — entry rewritten
- `tmp/ai_exchange/review_matvec_csc_omp_stack_overflow.md` — the round's record
- `devlog/README.md` — this entry

**No `src/` change.** `check_doc_claims.py` 37/37.

## Evidence standing

The crash mechanism is **verified** — reproducer, executed 2026-09-01 (`dl20260901_darwin_plugin_gate.md`).
Everything this round added is **reviewed, not verified**: source trace plus reviewer agreement, no
executable gate run today. D1 was confirmed against the source rather than by vote, and the one
refutation rests on the OpenMP spec plus an existing regression test whose last recorded green run
is `check-fast` 2026-08-30 — cited, not re-executed here.

**Next step is a gate, not another round** (protocol §11): R10's `OMP_DISPLAY_ENV=TRUE`, which
answers O2 and re-checks DR-155's severity, then O1+O3 as Christian's design call.

## Addendum, same day: rev. 2 of the plan (second-pass oversight sweep)

A cold re-read of the revised plan before any implementation, on Christian's instruction to look
for oversights rather than to code. The sweep found that **rev. 1's own R1 had introduced an
unimplementable requirement** while overcorrecting D1: "size from the team size that actually
enters the region" — but `omp_get_num_threads()` returns 1 outside the region and the allocation
must precede it. Filed as **D14**: sizing from the cap is correct and unavoidable; D1's defect
lives in the zeroing and the combine. R0 now picks one of two closed initialisation designs by
name — (a) each thread zeroes its own column inside the region and the combine runs over the
captured team size, or (b) zero-all before, combine-all after — with (a) recommended, since it is
simultaneously the correctness fix and the NUMA first-touch layout.

Smaller sharpenings, all annotated *(rev. 2)* in place: R10's `OMP_DISPLAY_ENV` may print the
unresolved ICV and gained a definitive fallback (a worker printing its own stack via
`pthread_get_stacksize_np` / `pthread_getattr_np`, which on Darwin also tests Grok's
libgomp-sets-2-MiB hypothesis directly); R1's bullets were declared branch-invariant with the C++
guard precedent named (`#ifdef OMP` + `<omp.h>` at `banner.cpp:21-22`,
`cl_SolverWrapper.cpp:16-17` — the same `OMP` token as Fortran); R8 was annotated with its end
state — the two `multiply` overloads are the kernels' **only** callers, so post-R8 an MKL build
never reaches `matvec_csc` and the struck scope guard becomes true; R6 logged two structural blind
spots (UMA cannot see first-touch; the D8-correct combine is serial, O(n·nt), and must be measured
specifically); and R1 gained an integer-kind note (`omp_get_max_threads()` is a C `int`;
`BELFEM_INT64` and >2³¹-element edge cases). No verdict changed, no step struck, no severity moved.
`check_doc_claims.py` 37/37 after the edit.


## Addendum 2: rev. 3 — O1 and O3 resolved (Christian)

The design fork closed in conversation, on a shape neither the plan nor the jury had proposed:
**a process-level scratch owned by the communicator**, `gComm.comm_mutable()->rbuffer()` /
`ibuffer()`, which Christian had just built. One buffer per rank rather than one per `SpMatrix`,
an owner whose lifetime is the process, and allocation failure that is an ordinary C++
`BELFEM_ERROR` — **which dissolves D5**, the missing Fortran abort primitive, by removing the
Fortran allocation entirely.

Two of the reviewer's objections were overruled, and the corrections are the useful part of the
record:

- **The dead-communicator guard was unnecessary.** The reviewer wanted `deregister()` to no-op if
  the communicator were already gone. Christian's ruling: the communicator closing last is a
  **contract** — every `main()` ends with `gComm->finalize()` — so an `SpMatrix` still working
  after that is a caller violation, not a case the buffer defends against.
- **The residency concern was answered with numbers, not a cap.** `max(rows,cols) * nthreads * 8`
  per rank: 152 MiB at the typical 1 M x 20, ~10 GB at the 10 M x 128 extreme. Acceptable on a
  machine that has 10 M dofs and 128 cores. Two structural reasons it cannot inflate quietly were
  added to the plan: it does not compound with rank count (local `n` shrinks as ranks grow, so a
  node lands near `n_global * T * 8` regardless of `R`), and `T` is already bounded by the standing
  no-oversubscription ruling. **O3 resolved: no cap, document the cost.**

An earlier proposal in the same exchange — derive `SpMatrix` from `CommunicationObject` so it
registers itself, and waive `mSwap` — was refined rather than adopted. Inheritance buys nothing the
buffer needs (access is a call, not a base class), and the existing registry's `free()`-at-finalize
contract suits process-lifetime singletons. Christian's counter-proposal, an explicit
register/deregister with a rescan to the high-water mark over live members, is the right mechanism
for the **swap** role — and `allocate_swap()` (`cl_SpMatrix.cpp:712-715`) is already the single
choke point, already called from all five sites where dimensions settle, already computing the
right quantity (`max(rows,cols)`, which is what the transposed path needs since the kernel enters
with `n = mNumCols`). It does **not** suit `ybuf`, whose demand depends on a thread count that
changes at runtime: `tests/sparse/test_SpMatrix.cpp:1816-1822` calls `omp_set_num_threads` in a
`{1,2,8}` loop between `multiply` calls on an already-constructed matrix, so a demand registered at
construction would be stale *low* on the second iteration. Hence the **two-roles split** now in the
plan: register the swap, keep `ybuf` grow-on-demand.

`mSwap` retirement is tracked as follow-on work under O1, out of the kernel change's scope, with
its three prerequisites named: two *named* buffer regions (the swap is live **across** the kernel
call at β ≠ 0), re-registration in the move-assign that currently steals it, and shrink hysteresis.

One forward-looking item came out of the extreme case: at 10 M x 128 the memory is affordable but
the D8-correct combine is ~1.3e9 **serial** additions per matvec. R6 now carries that as a named
measurement — at high thread counts the combine, not the allocation, is what bites.

**Correction within rev. 3, and it matters more than the resolution it amends:** the reviewer's
second structural reason for needing no cap — *"`T` is already bounded by the standing
no-oversubscription ruling"* — was **wrong**, and Christian said so. No-oversubscription is a strong
recommendation BELFEM cannot enforce; nothing stops a user setting `OMP_NUM_THREADS` as they like.
The sharper case is the *blind default*: with the variable unset the runtime returns the
affinity-mask popcount, **every hyperthread included** — `cl_SolverWrapper.cpp:435-437` already says
exactly this — so a 128-thread / 64-core node yields `T = 128` with the user having chosen nothing,
and a 2× buffer for zero speed on a bandwidth-bound scatter.

Three consequences now in the plan. Sizing must use `omp_get_max_threads()` **re-read on every
call** — undersizing is a buffer overrun, so the value the rank sees is the only safe bound however
unreasonable it is. The allocation-failure path stops being theoretical and becomes the main guard,
so its `BELFEM_ERROR` must name `n`, `nthreads` and the byte count. And the cost gets surfaced where
the user is already being told: `Wrapper::hatch_turtle()` (`cl_SolverWrapper.cpp:390-470`) already
fires on both oversubscription *and* the blind default, already computes every quantity needed, and
its banner already reads *"This costs time and memory"* — DR-155 gives that sentence a named number.
Its two limits are recorded with it: it is called only from MUMPS and STRUMPACK init
(`cl_SolverMUMPS.cpp:191`, `cl_SolverSTRUMPACK.cpp:175`), so other solvers and the eigen path never
see it.

The half of the argument that survived: the footprint does **not** compound with rank count, since
the local `n` shrinks as ranks grow. And one observation logged rather than reopened — the
hyperthread case is the one place a kernel-side cap would actually pay, so that is where to look
first if the buffer ever becomes binding.

**Still open after rev. 3:** D14's (a)/(b) initialisation design (recommendation: (a)), and O2.


## Addendum 3: rev. 4 — resolved by removing the parallelism

Christian closed DR-155 by deleting the parallelism rather than repairing it, after asking the
question none of the three reviewers had: **what does this threading actually buy?**

The call paths answer it. Both matvecs are master-only — `mSystemMatrix` is the full
`mNumberOfFreeDofs` square (`cl_FEM_DofMgr_SolverData.cpp:396-400`) and the residual at `:2368` is
inside `if ( mKernel->is_master() )` at `:2308`; the eigen matvecs sit inside
`if( mParent->parent()->is_master() )` (`cl_FEM_DofMgr_EigenValues.cpp:1704`) — so OpenMP was their
only parallelism, but they are ~2 flops per nonzero beside a full assembly and a direct
factorization per Newton step. And **no measurement of threaded-versus-serial exists anywhere in the
tree**: the kernel comment at `splinalg.f90:94-99` compares atomics against the reduction, two
threaded variants. A defect whose repair needed a buffer registry, a two-role split, reduction
semantics reproduced by hand and a 512 MiB scratch budget was closed with a default-OFF build token.

**What landed.** A new `BELFEM_OMP` define, distinct from `OMP`: `USE_BELFEM_OPENMP` (default OFF,
`CMakeLists.txt:88`) with the define appended under `USE_OPENMP AND USE_BELFEM_OPENMP` (`:179-181`).
All 12 `!$omp` blocks in `splinalg.f90`, `arpacktools.f90` and `parpacktools.f90` wrapped in
`#ifdef BELFEM_OMP`; the four `use omp_lib` guards in the ARPACK pair moved to it, and
`splinalg.f90`'s two subroutine headers gained the same guarded import at Christian's request — the
directives do not need it, but a future `omp_get_*` in that kernel would fail under `implicit none`
without it.

**Two reviewer positions corrected on the way.** A **Darwin-only** ban was argued against and
dropped in favour of OFF everywhere: O2 is unresolved, so a carve-out could have left the identical
cliff live on Linux, and platform-conditional threading means platform-divergent rounding in a
framework that validates against published benchmarks. And the proposed **blanket**
`OMP` -> `BELFEM_OMP` sweep was narrowed after an inventory of all 13 guard sites showed **nine must
stay** — moving `hatch_turtle()` (`cl_SolverWrapper.cpp`), PARDISO's `gParameters(3)`
(`pardisotools.f90`) or the banner's `Threads Used` line would have silently deleted the
oversubscription warning, serialised a third-party solver, and hidden the thread budget, all on the
default build. Those three now carry explicit "NOT a mistake" comments naming the reason, at
Christian's instruction.

**Asked and answered:** would a buffer registry have prevented a *third-party* crash? No. MUMPS,
STRUMPACK and MKL allocate their own workspace on their own threads and never ask BELFEM for a
buffer; the instruments that do work there are the in-wrapper MUMPS ICNTL(14) ladder (DR-106) and
`hatch_turtle()` — which is the second reason that guard had to stay on `OMP`. The `gesv` pivot use
of `ibuffer()` was ruled noise on the same grounds and deferred past the release, leaving
`rbuffer`/`ibuffer` without a consumer for now.

**What this does NOT do**, and both are documented so they are not rediscovered the hard way: the
process is still multi-threaded (STRUMPACK, MKL, threaded BLAS spawn workers on the same 2 MiB
Darwin stacks, so `OMP_STACKSIZE` remains a live knob and a `gomp_thread_start` fault is not
automatically a BELFEM bug); and the kernel is **not repaired** — `USE_BELFEM_OPENMP=ON` re-arms
DR-155, which the source comment at the reduction now says outright.

**Side benefit:** `matvec_csc` is now bitwise deterministic across thread counts, removing one
contributor to the convergence-path noise documented in `doc/parallel_execution.md`.

**Evidence standing.** Executed: the preprocessor check (11 `!$omp` lines with the define, 0
without) and `gfortran -fsyntax-only` green on all three kernels in **both** token states (6/6).
**Owed, and Christian runs builds:** R13, a build plus `make check` in both token states — the ON
build being the proof that the token is a real switch and not a one-way deletion — and R5, the
`tape_quench_usermat` deck on a default environment with no `OMP_STACKSIZE`. Until those run this is
**reviewed, not verified**.

**Files touched in `src/` this session:** `CMakeLists.txt`, `src/sparse/splinalg.f90`,
`arpacktools.f90`, `parpacktools.f90`, `pardisotools.f90`, `cl_SolverWrapper.cpp`,
`src/core/banner.cpp`; docs `doc/parallel_execution.md`, `doc/coding_philosophy.md`, `CLAUDE.md`.


## Addendum 4: the MKL half fixed too

Christian: *"`SpMatrix::multiply` should use MKL if MKL is available."* D2's defect closed the same
session.

The 2-argument `SpMatrix::multiply` no longer switches on `mType` and calls the Fortran kernels; it
forwards to the 5-argument overload as `this->multiply( aX, aY, 1.0, 0.0, false )`, which carries
the `BELFEM_MKL` branch. Three consequences: MKL builds now reach `mkl_sparse_d_mv` on the ARPACK
inverse-iteration path; the overload picks up the two dimension `BELFEM_ASSERT`s it never had, which
was D2's drift between two implementations of one operation; and **both surviving
`matvec_csc`/`matvec_csr` call sites are now inside the `#else` of `BELFEM_MKL`** — so the scope
guard this plan opened by proving *false* ("`USE_MKL=ON` never reaches the kernel") is true as an end
state, for the first time.

The forward is instruction-for-instruction equivalent on the Fortran branch, checked branch by
branch: `aBeta == 0.0` skips both the `mSwap` save and the `aY += mSwap` add-back, `aAlpha == 1.0`
skips `aY *= aAlpha`, and `aTransposedFlag == false` reduces `tRunCsr`, `tIndices`, `tNumRows` and
`tNumCols` to exactly the switch that was deleted.

**The fill is `#ifdef BELFEM_MKL` only, and the asymmetry is the interesting part.** The Fortran
kernels zero `y` themselves, so every caller of this overload is entitled to hand it an
uninitialised `aY` — and the eigen driver does, `set_size()` with no fill value and then straight
into the multiply. `mkl_sparse_d_mv` computes `y := alpha*op(A)*x + beta*y`, and the
Inspector-Executor sparse BLAS is not documented to short-circuit `beta == 0`; `0.0 * y` is 0 for
finite garbage but NaN for a NaN or Inf bit pattern, which uninitialised heap supplies. On the
Fortran branch the kernel's own zeroing already covers it, so an unconditional fill would be a
second pointless O(n) pass on a path that runs once per ARPACK iteration. The reasoning lives in the
source comment, not only in this log.

**Honest limit, and it is the sharp one:** `config/linalg/config_mkl.cmake:2-3` makes MKL a hard
configure error on macOS, so **the one code path this change exists to enable is the one that cannot
be compiled on the machine where it was written.** What ran here: `g++ -fsyntax-only` green on
`cl_SpMatrix.cpp` (netlib/Blaze), and a compile probe confirming `aY.fill( 0.0 )` on
`Vector<real>`. A Linux `-DUSE_MKL=ON` build exercising an eigenvalue case is owed and is folded
into R13.
