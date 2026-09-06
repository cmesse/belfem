# DR-144: MUMPS Instance-Lifecycle Accounting

**Date:** 2026-08-29
**Purpose:** Close the seven lifecycle-accounting gaps filed as DR-144, so that C++ wrapper state
(`MUMPS::mSolverID`, `MUMPS::mInitialized`, `Wrapper::mIsInitialized`) can never disagree with
Fortran pool occupancy (`gOccupied`), so that the exhaustion path is collective, and so that a
solver failure is reported as a solver failure rather than as a basis-budget refusal.
**Module:** `src/sparse` (wrapper + Fortran shim), `src/fem/kernel` (the diagnostic's user message)
**AIs involved:** Claude (exploration + plan), Codex (audit), Grok (third voice)
**Status:** **ROW STRUCK 2026-08-30 ( Christian's ruling ) — plan complete.** 11 of 11 steps done:
R12/R13 landed 2026-08-30 in the DR-106 joint session (`todo/dr106_dr144_joint_session.md`) and the
same-day rebuild compiled them into a `check-fast`-green suite ( 15/15 incl. `sparse` ). Archived row:
`todo/debt_register_closed.md`. The green build was confirmed to have actually exercised the
new test rather than inferred from the suite result (per the DR-117 lesson): `test_sparse` is dated
22:06 against a last source edit of 21:51, all four MUMPS test names are compiled into that binary,
and `sparse` is in `check`'s dependency list (`CMakeLists.txt:422`). **Five steps (R3, R4, R5, R6,
R11) still land reviewed-not-verified** — they have no reachable in-tree trigger, so the green build
says nothing about them, and must not be read as if it did.** Two jury rounds, both vendors, both rounds. The code round found a **P0 that would have
made R1 ship inert** — a parallel session's `Detailed` -> `Verbose` sweep had demoted
`report_unavailable` above the level the debug executables run at, so the new message could not
print. Fixed, along with four smaller reporting defects and a test that did not encode its own
hazard. Two further findings belong to other sessions' uncommitted code and are carried as **G10 and
G11 on this row** (steps R12/R13, deferred by ownership rather than difficulty) — they were briefly
mis-filed as new DR IDs and withdrawn; see §4.2.** Christian
approved the v2 plan and the DR-90 disposition in one ruling. All nine gaps G1-G9 have landed,
including every correction the jury round required. Syntax-checked clean against the tree's
Armadillo config (`cl_SolverMUMPS.cpp`, `cl_SolverWrapper.cpp`, `cl_FEM_DofMgr_EigenValues.cpp`,
`test_Solver.cpp`) and `gfortran -fsyntax-only` clean on `mumpstools.f90`, the latter confirmed real
by a deliberate negative control. **That is a syntax gate, not a correctness gate** — `make check`
(R8) is still owed and is the only thing that can run the new test. **Reviewed, not verified.** Scope approved by Christian 2026-08-29 ("DR-144 in full, one round").
Jury round returned Codex **"revise before implementation"** and Grok **"approve with the required
corrections"**. Both confirmed the gap inventory and F1 (the round's central finding); both rejected
the *step shapes* for R4 and R6 and the off-by-one in R9. v1's O2 recommendation was refuted by a
concrete cross-wrapper destroy sequence. An **eighth gap (G8)** was found by the round and is now in
scope. See §4.0.

> **Scope guards:**
> - Triangle extraction for a symmetric MUMPS mode is OUT of scope. That is DR-140's residue and
>   has its own file (`todo/mumps_symmetric_triangle_extraction.md`).
> - Nothing here changes the soft-fail contract, the frozen-factorization scope, or `select_job()`.
>   The only new coupling is that `mInitialized` acquires one earlier write site.
> - No new deck key, no new public API. `mumpstools_set_num_max_solvers` is either guarded or
>   removed (O1), not extended into a reallocating pool.
> - **G6 is live today** and is the only one that is: it mislabels a real message on the running
>   conditioning diagnostic. G1, G2, G5, G7 and G9 are latent with no in-tree trigger. G3 and G8 are
>   latent but reachable through one public call (`wrapper()->initialize()`), not unreachable — v1
>   oversold that, corrected by Grok. G4 is a real wrapper bug with no matrix-RHS-only caller today.
>   The rest of the round buys correctness under future callers.

---

## 1. Current Behaviour and How It Fails

The pool is eight `DMUMPS_STRUC` slots with an occupancy table (`mumpstools.f90:29-32`). The C++
side tracks the same lifecycle three times over, in members that are set at different moments:

| member | set where | read where |
|---|---|---|
| `MUMPS::mSolverID` | `mumpstools_create_solver` out-arg (`cl_SolverMUMPS.cpp:233-238`) | the `<= 0` guards at `:240`, `:501`, `:661` |
| `MUMPS::mInitialized` | three SOLVE sites only: `:609`, `:643`, `:783` | **exactly one place**, `MUMPS::free():284` |
| `Wrapper::mIsInitialized` | `Wrapper::initialize()` (`cl_SolverWrapper.cpp:69`), called at `cl_SolverMUMPS.cpp:267` | `is_initialized()` (`cl_SolverWrapper.hpp:382-386`) |

**Bottom line:** the three flags are set at three different points in the create/solve/free
sequence, and only the Fortran occupancy table is authoritative. Every gap below is one place where
that divergence becomes observable.

## 2. Architecture: Make the Fortran Table the Single Authority

The fix direction is not to add a fourth flag but to make each existing one mean exactly one thing
and be set at the moment that thing becomes true:

- `gOccupied(k)` = a live `DMUMPS_STRUC` exists in slot `k`. Therefore the occupancy write must be
  rolled back when `DMUMPS( JOB = -1 )` fails, and the create must hand back `aSolverID = -1` so the
  existing C++ `<= 0` gate catches it without a second decision path.
- `MUMPS::mInitialized` = a Fortran instance exists that `free()` must release with `JOB = -2`. Its
  own comment at `:604` already says this. It is therefore established by a successful `initialize()`,
  not by a solve.
- `Wrapper::mIsInitialized` = the wrapper holds a usable solver ID, so `Solver::solve` need not
  retry the create. Already correct after the DR-143 round.

**Alternative considered and — as v1 stated it — refuted.** Keying `free()` on `mSolverID > 0`
instead of on `mInitialized` closes G3 and G4 equally and leaves one handle instead of two, which is
why Codex recommended it. But the ID is never cleared (`MUMPS::free()` at `:275-351` does not write
it), and an ID that outlives its slot names a slot someone else now owns: wrapper A frees slot 5 and
keeps `mSolverID = 5`; wrapper B's create reuses slot 5 by first-free scan; `~A` issues `JOB = -2`
against **B's live instance**. `gOccupied` reads 1, so the shim does not skip it. Today only
`mInitialized` prevents that. The fix is therefore not to choose between the two handles but to
**zero the ID on release** (R10); the consolidation onto a single handle becomes safe only
afterwards, and is deliberately left to a follow-up. See O2.

## 3. Gap Table

Seven gaps as filed, re-verified against the working tree on 2026-08-29 (the tree has moved since
filing: DR-143's first-free scan and the soft-fail create both landed). Class **(c)** = must be
handled explicitly; **(b)** = open question, see §5.

| # | Gap | Still open? | Class | Citation |
|---|---|---|---|---|
| G1 | A failed `DMUMPS( JOB = -1 )` leaves its slot occupied: `gNumSolvers` and `gOccupied` are written BEFORE the call and nothing rolls back on error | yes | (c) | `mumpstools.f90:94-107` |
| G2 | `MUMPS::initialize` ignores the `tInfo` the shim returns; a positive ID with a failing `INFO(1)` passes the `mSolverID > 0` gate | yes | (c) | `cl_SolverMUMPS.cpp:233-240` |
| G3 | `mInitialized` is set only in `solve`, so initialize-then-free-before-any-solve skips `mumpstools_free_solver` and leaks the slot | yes | (c) | set at `:609,:643,:783`; gate at `:284` |
| G4 | The matrix-RHS overload sets `mInitialized` only in its soft-fail branch (`:783`) and never on success; the vector overload sets it at `:643`. An instance used exclusively through successful matrix-RHS solves never receives `JOB = -2` | yes | (c) | `cl_SolverMUMPS.cpp:643` vs `:783` |
| G5 | Pool exhaustion is NON-collective: the `-1000` branch returns without the `MPI_BARRIER` the success path takes | yes | (c) | `mumpstools.f90:92` vs `:103` |
| G6 | `report_unavailable`'s reason ternary has no `EigenOutcome::SolverFailed` arm, so a solver failure is reported as "no legal Krylov subspace fits the basis budget" | yes | (c) | `cl_FEM_DofMgr_EigenValues.cpp:1561-1571`; `SolverFailed = 6` at `.hpp:45`, set at `.cpp:1389` |
| G7 | `mumpstools_free_solver` indexes `gOccupied( aSolverID )` with no bound check, and `mumpstools_set_num_max_solvers` does not reallocate the pool | yes | (b) for the setter, (c) for the bound check | `mumpstools.f90:142`, `:41-48`; declared at `mumpstools.hpp:28` and DEFINED in Fortran, with no C++ caller ( v1 said "declaration-only", corrected by Grok ) |
| **G8** | **NEW, found by the audit round.** A second `MUMPS::initialize()` on the same wrapper creates a slot and overwrites `mSolverID` — losing the first ID — and only THEN calls `Wrapper::initialize()`, which throws on its already-initialized check. Both slots leak, and the first is now unreachable. Not one of the seven filed gaps. Reachable only through the public `wrapper()` back door, since `Solver::solve` guards on `is_initialized()` | yes | (c) | create at `cl_SolverMUMPS.cpp:233-238`, `Wrapper::initialize()` throw at `cl_SolverWrapper.cpp:64-66`, call order at `cl_SolverMUMPS.cpp:267` |
| **G12** | **NEW, raised by Christian 2026-08-29 from reading a footer.** `MUMPS ADD COND2` prints `1.00e0` on most timesteps, and that value is not a measurement: MUMPS initializes both condition numbers to 1.0, partitions the rows, and SKIPS the entire COND2 estimator when no row lands in the second category — `IF (.NOT.LCOND2) GOTO 170`, a return past the estimator and past `ERX = ERX + OMEGA(2)*COND(2)`. The same partition leaves `OMEGA(2)` at exactly 0. So the reported 1.0 is an untouched initializer for a term that does not exist. Traced to the shipped MUMPS source, not to the user guide | yes | (c) | `tmp/MUMPS_5.9.1/src/dsol_aux.F` — `DMUMPS_SOL_OMEGA:886-899`, `DMUMPS_SOL_LCOND:959-1033` |
| **G10** | **NEW, found by the CODE round; folded in under extend-don't-branch rather than branched.** `EigenValues::mK` is captured once and can dangle on a WORKER: `run_shift_invert` fills a null `mK` from `mParent->jacobian()` and never revisits it; `jacobian()` returns `mSystemMatrix`, which `SolverData::reset()` deletes; `DofManager::set_equation()` calls `reset()`; `EigenValues::reset()` clears only `mMatrixFlag`. The master recovers through `link_matrix()` on the next `run()`, a worker does not. The block's own comment ( "`link_matrix()` only ever runs on the master" ) is FALSE — `compute_matrices()` and `run_arpack()` both call it unguarded | yes | (c), DEFERRED | `cl_FEM_DofMgr_EigenValues.cpp` shift-invert entry; `link_matrix()` calls at `:370`, `:456`, `:605` |
| **G11** | **NEW, found by the CODE round; folded in the same way.** The `mSolverID <= 0` bail in BOTH `MUMPS::solve` overloads returns silently even when `soft_fail()` is false, so a caller that solves without a successful `initialize()` gets its solution vector untouched, `failed() == true`, and no abort in debug or release. The `BELFEM_ASSERT( is_initialized() )` below it is unreachable for that case. The comment claims it covers only the soft-failed-create path; the code never tests `soft_fail()` | yes | (c), DEFERRED | `cl_SolverMUMPS.cpp`, both solve overloads |
| **G9** | **NEW, found by the audit round.** `mumpstools_create_solver` passes `aInfo` as `MPI_BARRIER`'s status argument (`:103`). `aInfo` is `integer( int_t )` (`:64`), which under `BELFEM_INT64` is 64-bit, while `MPI_BARRIER` requires a default `integer`. `mumpstools_free_solver` gets this right with its own `integer :: ierr` (`:127`). The value is also immediately overwritten at `:107`, so the barrier status is discarded either way | yes | (c) | `mumpstools.f90:64`, `:103`, `:107` vs `:127`, `:146` |

### 3.1 Cross-cutting findings

**F1. G3 and G4 are one fix, not two, and they are cheaper than the register claims.** The DR-144
row says they "need their own plan round because the `mInitialized` semantics carry the
frozen-factorization and soft-fail machinery." Re-measured: `MUMPS::mInitialized` is **read at
exactly one site**, `MUMPS::free():284`. It touches neither `mFrozen` nor `soft_fail()`; the
frozen scope is carried by `mFrozen.mArmed` and `mMatrix`, and the soft-fail retry is keyed on
`Wrapper::mIsInitialized` through `is_initialized()`. Setting `mInitialized = true` beside the
existing `Wrapper::initialize()` call at `cl_SolverMUMPS.cpp:267` therefore closes both gaps and
makes all three solve-site writes redundant. Confidence: high (grep is exhaustive over `src/`; the
member is private and the header declares it at `cl_SolverMUMPS.hpp:75`).

**F2. G1 subsumes most of G2.** If the shim rolls back and returns `aSolverID = -1` on a failed
`JOB = -1`, the existing `mSolverID <= 0` branch already handles it, and `print_soft_fail` reports
the real MUMPS `INFO(1)` rather than the shim's `-1000`. G2 is then defence in depth rather than the
only line: the two sides must not silently depend on each other. Both are still in scope.

**F3. Two error strings are stale by-catch from DR-140's rename.** `Wrapper::get_cond1()` reports
`"get_cond0() is not implemented"` and `get_cond2()` reports `"get_cond1() ..."`
(`cl_SolverWrapper.cpp:180-199`). Names a function that no longer exists. Picked up here because it
is in the same file set, and it is a message, not behaviour.

## 4. Ordered Steps

> **Revised after the 2026-08-29 jury round.** v1's R4 and R6 were rejected as unsafe by both
> auditors; R9 had an off-by-one; R1's message was too narrow; R2 was missing a sequencing
> constraint. The v1 text is struck rather than deleted so the correction is legible.

- [x] **R1 (G6)** Rewrite the reason ternary at `cl_FEM_DofMgr_EigenValues.cpp:1561-1571`:
      add an `EigenOutcome::SolverFailed` arm, add an **explicit `Infeasible` arm**, and make the
      default print the unmatched numeric outcome. ~~v1: "name the factorization/frozen-solve
      failure"~~ — refuted by Grok: `SolverFailed` is set from `wrapper()->failed()`
      (`.cpp:1383-1389`), which is ALSO raised by a failed or exhausted create
      (`cl_SolverMUMPS.cpp:249-258`), so a factorization-specific sentence would relabel pool
      exhaustion as a factorization failure — better than the basis-budget lie, still the wrong
      subsystem. Word it as a linear-solver/instance failure without pinning it to factors. The
      explicit `Infeasible` arm plus a numeric default is what stops the next enumerator repeating
      G6. Independent of every other step; highest value per line, and the only gap live today.
- [x] **R2 (G3+G4)** Set `mInitialized = true` in `MUMPS::initialize()`, **after** the
      `Wrapper::initialize()` call at `cl_SolverMUMPS.cpp:267` (so a throw there cannot leave
      `mInitialized` true against `mIsInitialized` false), and delete the three now-redundant
      solve-site writes at `:609`, `:643`, `:783`. **Delete, do not leave, the comments at
      `:604-605` and `:778-779`** that justify those writes — after R2 they describe a lie.
      Comment the invariant at the declaration (`cl_SolverMUMPS.hpp:75`): *set by a successful
      create, cleared by `free()`, read only to decide whether `JOB = -2` is owed*.
      **(must land in the same edit as R3 and R5)**
- [x] **R3 (G1)** In `mumpstools_create_solver`, roll back on a failed `JOB = -1`: if
      `tMUMPS%INFO( 1 ) < 0`, clear `gOccupied( aSolverID )`, decrement `gNumSolvers`, `nullify`
      the slot's `irn`/`jcn`/`A`/`rhs` as the free path does (`:157-160`), and set
      `aSolverID = -1` while keeping `aInfo = tMUMPS%INFO( 1 )` so the real cause survives. This
      also makes the subroutine's own header comment true — `:54` already documents
      "-1 on failure". **Open sub-question (Grok, confidence ~60%, MUMPS 5.5.1 guide not in
      `literature/`):** whether a `JOB = -1` that allocated and then failed still owes a
      `JOB = -2`. Occupancy rollback is the right *pool* accounting either way; a MUMPS-side leak
      would be invisible to this round. Record the uncertainty, do not paper over it.
      **(same edit as R2, R5)**
- [x] **R4 (G5 + G9)** ~~Move the `MPI_BARRIER` AHEAD of the first-free scan.~~ **REJECTED by both
      auditors.** After the barrier each rank still scans its own process-local `gOccupied`
      (`:79-82`); on disagreement one rank returns through the exhaustion arm while the other
      enters `DMUMPS`, which is itself collective on `MPI_COMM_WORLD` (`:98`, `:105`) — the hang
      simply moves from the explicit barrier into MUMPS. Corrected shape:
      **(a)** keep a barrier immediately before `JOB = -1` where it is today (`:103`);
      **(b)** give it its own `integer :: ierr`, matching `free_solver:127` — this fixes **G9**,
      the `BELFEM_INT64` kind mismatch, and stops the status being clobbered at `:107`;
      **(c)** add a matching barrier on the exhaustion arm so both *agreeing* outcomes are
      collective. This is hygiene, **not** a mixed-occupancy fix, and the plan must stop claiming
      otherwise. Making mixed occupancy a loud error rather than a hang needs an allreduce of the
      chosen slot / exhaustion flag before `DMUMPS` (Codex's proposal) — a different step, logged
      as O4, not smuggled into a barrier move. **(after: R3, same subroutine)**
- [x] **R5 (G2)** At `cl_SolverMUMPS.cpp:240`, gate on `mSolverID <= 0 || tInfo < 0` rather than on
      the ID alone, and carry `tInfo` into `mInfo( 0 )` on both arms. Note `tInfo > 0` is a MUMPS
      *warning* and must NOT abort the create. **(same edit as R2, R3 — see D6: R2 without R3+R5
      would mark a G1-failed create as owing `JOB = -2`)**
- [x] **R6 (G7)** Bound-check the free, in this order and no other:
      **(1)** `if( .not. allocated( gOccupied ) ) return` — the existing guard at `:138`, kept
      FIRST; **(2)** `if( aSolverID < 1 .or. aSolverID > size( gOccupied ) ) return`;
      **(3)** only then the existing `gOccupied( aSolverID )` test at `:142`.
      ~~v1 put the bound check ahead of the `allocated` guard~~ — that is **Fortran-illegal**:
      `size()` on an unallocated allocatable is undefined, and the guard exists precisely because
      a full drain deallocates (`:174-177`). Both auditors flagged it independently, both HIGH.
      Use `size( gOccupied )`, **not** `gMaxNumSolvers` — after a late setter call the two differ,
      which is the other half of G7. Resolve the setter per O1.
- [x] **R10 (from O2/O3, new)** Zero `mSolverID` in `MUMPS::free()` after the Fortran call returns,
      and again before every create attempt. **This is the correction that makes the whole round
      safe.** Grok's sequence: wrapper A frees slot 5 while keeping `mSolverID = 5`; wrapper B's
      create reuses slot 5 by first-free scan; `~A` runs `free()` again and issues `JOB = -2`
      against **B's live instance**. That is a cross-wrapper destroy, not a benign double free —
      `gOccupied` reads 1, so the shim does not skip it. Today only `mInitialized` prevents it.
      **(same edit as R2)**
- [x] **R11 (G8, new)** Add an always-active `BELFEM_ERROR( ! this->is_initialized(), … )` at the
      **top** of `MUMPS::initialize()`, before `mumpstools_create_solver`. Setup path, runs once,
      so `BELFEM_ERROR` is the right tier (`doc/coding_philosophy.md`, frequency rule).
- [x] **R7 (F3, by-catch — NOT one of the gaps)** Correct the two stale function names in the
      `Wrapper::get_cond1`/`get_cond2` error strings (`cl_SolverWrapper.cpp:180-199`; v1 cited
      `:177-193`, off by three). Codex flagged this as unmapped scope expansion and is
      procedurally right: it is DR-140 rename by-catch, not DR-144. Kept because it is two string
      literals in the same file set, and flagged here for Christian rather than folded in silently.
- [x] **R14 (G12)** Suppress the ADD COND2 footer row when the term is empty. Added
      `Wrapper::get_omega2()` / `MUMPS::get_omega2()` returning `mRInfoG( 7 )` ( RINFOG(8) ),
      captured into a new slot 3 of `mConditionNumbers1/2`, and gated the print on a file-static
      `add_cond2_is_real( omega2 )`. **Keyed on omega2, NOT on `COND2 == 1.0`** — the sentinel value
      is indistinguishable from a genuine estimate near one, while the empty category is the actual
      cause. **Omitted rather than printed as `n/a`**, on Christian's call and because the two say
      different things: `n/a` means "asked, unavailable", and this means "no such quantity for this
      solve". The uncaptured case ( omega2 is NaN ) still prints `n/a` exactly as before, so only the
      one case that would lie is dropped. **Not suppressed unconditionally**: COND2 is a real number
      in a minority of steps — 3 of 14 and 5 of 41 in the tapestack3d logs, ranging 5.4e2 to 7.9e3 —
      and those are the steps where part of the residual went through the padded denominator.
      **COND1 has the identical structure** ( `IF (.NOT.LCOND1) GOTO 130` ) and deliberately did NOT
      get the guard: omega1 is nonzero in every observed step, so the guard could not be exercised,
      and an unexercised guard is worse than a visible asymmetry. Logged as O6.
      **Audited 2026-08-29, both vendors; the JUSTIFICATION was refuted and the gate survived.**
      Both found that `omega2 == 0` does not prove the category was empty — `IW( i, 1 ) = 2` is
      assigned outside the `IF ( TAU .GT. ZERO )` guard ( `dsol_aux.F:891-900` ), so a `TAU == 0`
      row joins the category, sets `LCOND2`, runs the estimator, and leaves omega2 at zero. My
      comments asserted the biconditional. Codex rejected the patch over it; Grok approved the code
      and rejected only the comments. **Grok's reading adopted:** at `omega2 == 0` the term is zero
      either way, so a suppressed COND2 cannot move the bound it exists to explain, and the real
      discriminator ( `IW` ) is internal and unrecoverable from the exposed scalars. The rare
      suppression of a genuinely estimated COND2 is accepted and documented in the helper rather
      than hidden. Comments rewritten to the forward-error-term framing; base-class contract
      restated without MUMPS internals; `compute_conditioning`'s slot-ownership comment extended to
      slot 3; helper renamed `add_cond2_is_real` -> `cond2_row_wanted`, because it returns true for
      NaN and the old name contradicted its own contract.
- [x] **R12 (G10) — DONE 2026-08-30 in the DR-106 joint session** (`todo/dr106_dr144_joint_session.md`,
      step R6; audited by both vendors at xhigh, reviewed-not-verified). Assign `mK = mParent->jacobian()` unconditionally on
      every entry with a `BELFEM_ERROR( mK != nullptr, … )` beside it — the same cost as the null
      test, and what `link_matrix()` already does — and correct the false master-only claim in the
      comment. **Blocked on ownership, not on difficulty:** the block is a different session's
      uncommitted work.
- [x] **R13 (G11) — DONE 2026-08-30 in the DR-106 joint session** (same round as R12, step R5,
      else-shape per the code jury). `if( ! this->soft_fail() ) BELFEM_ERROR( false, … ) ;
      else { this->flag_failure() ; return ; }` in both overloads, so the silent path exists only for
      the caller that asked for it. Same ownership block — pre-existing DR-143 work.
- [x] **R8 — GATE.** *(RAN GREEN 2026-08-29, Christian.)* `make check` green, including `tests/sparse/test_Solver.cpp`. Shared with
      DR-140's owed rerun, so the two rows close together. Process coupling only — DR-144 must not
      block on DR-140's triangle-extraction residue.
- [x] **R9 — GATE.** *(test written and RUN GREEN 2026-08-29 as part of R8.)* A unit test driving the G3/G4/G8 path: `tSolver.wrapper()->initialize( tM,
      SymmetryMode::Unsymmetric, 1 )` then `tSolver.free()`, in a loop. Both `Solver::wrapper()`
      (`cl_Solver.hpp:202-205`) and `Wrapper::initialize( SpMatrix&, … )`
      (`cl_SolverWrapper.hpp:170-178`) are public, so no new API is needed.
      ~~"the eighth cycle exhausts the pool"~~ — **off by one**, caught by both auditors: with
      eight empty slots and no pinned occupant, cycles 1-8 each leak one and **cycle 9** is the
      first exhaustion. Use the existing test's `k < 12` loop, or pin an instance as
      `MUMPSPoolRecyclesFreedSlots` does (`tests/sparse/test_Solver.cpp:205-211`). **Do not
      replace that test** — it gates the first-free allocator on the solve+free path; R9 gates G3
      on the initialize+free path. Different gaps, both needed.

**Not gated executably:** R3, R4, R5, R6 and R11 have no reachable trigger in the tree, so they land
**reviewed, not verified**. Say so in the devlog and the register row rather than letting the shared
`make check` imply coverage it does not give. Optional cheap probe for R6 (Grok): call the C API
with IDs `0`, `-1`, `99` after a live create — must not crash and must not `JOB = -2` the live slot.

### 4.2 Process correction — the two by-catch findings were mis-filed, 2026-08-29

D20 and D21 were first filed as **DR-149 and DR-150**. Christian rejected both. Two things were
wrong, and only one of them was a policy miss:

1. **The extend-don't-branch policy** (`debt_register.md`, Maintenance, set 2026-08-29): *when a
   session's own gate or audit turns up a finding in the same contract or subsystem as the row being
   worked, widen that row's scope cell and put the work in its plan, rather than branching a fresh
   ID. New rows are for genuinely separate subsystems. The register had been inflating one ID per
   finding.* Both findings are the MUMPS wrapper lifecycle and the conditioning diagnostic's own
   reporter — the same subsystem this row already spans, since G6 lives in
   `cl_FEM_DofMgr_EigenValues.cpp` too. Neither is a separate subsystem. They are now **G10 and
   G11** on DR-144, with R12 and R13 as their steps.

2. **DR-149 was an ID collision, which is worse.** A parallel session had filed and withdrawn that
   exact ID hours earlier (its evidence folded into DR-112, its work into
   `todo/closed/restart_timestep_consistency.md`) — and the register's rule is that *IDs are stable and are
   never reused*, withdrawn or not. The scan that missed it matched only lines beginning `| DR-`, so
   it saw rows and not prose, and a withdrawn ID lives in prose. **Scan the whole file for an ID
   before claiming it, not just the row starts.**

Worth recording as its own lesson: the policy that would have prevented this was set the same day,
in the file being edited, by another session — and this session read that file five times without
reading the rules section it sits in.

### 4.1 Audit Round 2 — CODE jury, 2026-08-29 (Codex + Grok, blind/parallel)

Codex: **revise before merge**. Grok: **ship the lifecycle mechanics after fixing D1 and D2 — do not
treat this as a clean bill of health.** Both independently confirmed the lifecycle core (V1-V10 in
Grok's leg): `mInitialized` at successful create, occupancy rollback, id cleared on free, double-init
guard, the barrier kind fix, and the rename's internal consistency. Neither modified source.

**Both auditors' top finding was the same, and it was not in the diff's new code.** The
`InfoLevel::Detailed` -> `Verbose` sweep across this file — a PARALLEL session's uncommitted work,
which the pre-registration explicitly excluded from the round's scope — also demoted
`report_unavailable`. `Detailed` is 3, `Verbose` is 4, `message()` prints on `<=`, and both debug
executables construct `gLog( InfoLevel::Detailed )` (`belfem.cpp:38`). So R1's new `SolverFailed`
sentence would never have printed in the one configuration where somebody is debugging a quiet
diagnostic. **The fix would have shipped inert.** Restored to `Detailed` for those two calls only;
every telemetry demotion around them left alone, which is precisely what both auditors recommended.

| # | Defect | Severity | Found by | Disposition |
|---|---|---|---|---|
| D13 | `report_unavailable` demoted to `Verbose`, above the debug executables' `Detailed` — R1's message invisible where it matters | **P0** | Codex + Grok (independent) | FIXED — those two calls restored to `Detailed`, telemetry demotions untouched |
| D14 | `print_soft_fail` promises "timestep will be cut" on a failed CREATE; the diagnostic consumer arms soft-fail itself and continues the step | HIGH | Grok | FIXED — the `-1000` case gets its own line. Same wrong-subsystem class as G6 |
| D15 | The `NotPositiveDefinite` string still says "no spectral fold can reach it"; the fold was retired for shift-invert, which CAN return a negative lambda_min | HIGH | Grok | FIXED in the same ternary R1 rewrote |
| D16 | `error_message()` has no `-1000` arm, so the HARD-fail create path (soft-fail off) prints a bare `-1000` | LOW | Grok | FIXED — decoded at the `BELFEM_ERROR` |
| D17 | The C++ comment asserts the create-fail branch "cannot desynchronize the ranks"; the Fortran correctly refuses that claim | LOW | Grok | FIXED — C++ comment weakened to match |
| D18 | `MUMPSPoolRecyclesFreedSlots` does not encode the cross-wrapper destroy its own comments describe: it would pass if a churn free destroyed the pin's slot | HIGH (test gap) | Grok | FIXED — post-loop solve through the pin, with `EXPECT_NEAR` |
| D19 | The new test does not exercise the matrix-RHS overload, and nothing here is MPI | MEDIUM (test gap) | Codex + Grok | ACCEPTED, not fixed — no in-tree matrix-RHS-only consumer to drive it. Said out loud in the test's own comment |
| D20 | `EigenValues::mK` captured once; `set_equation()` deletes the matrix and a worker keeps the stale pointer. The block's "link_matrix() only runs on the master" comment is false | MEDIUM (~70% live) | Grok | **NOT FIXED — folded into DR-144 as G10** ( see R12 ). A parallel session's uncommitted block; the rule is not to edit over it. Briefly filed as DR-149 and WITHDRAWN — see §4.2 |
| D21 | `mSolverID <= 0` bail in both solve overloads is silent even with soft-fail OFF | MEDIUM (~60% reachable) | Grok | **NOT FIXED — folded into DR-144 as G11** ( see R13 ). Same reason, pre-existing DR-143 work. Briefly filed as DR-150 and WITHDRAWN |
| D22 | Failed `JOB = -1` slots are recycled without a proven `JOB = -2`; the MUMPS contract is unread | MEDIUM (~70%) | Codex + Grok | ACCEPTED as standing debt — already documented in the shim comment, deliberately not papered over. Occupancy accounting is right either way |
| D23 | Pool is now hard-capped at 8 with the setter deleted; a multi-field deck each holding a diagnostic instance could still reach it | MEDIUM | Grok | LOGGED as O5 — a consequence of resolving O1 to (b), and worth knowing |

**Claude's code-round pre-registration scored:** Q1 (removing the solve-site writes breaks the
soft-fail reuse contract) — NOT found; Grok explicitly verified the opposite (V1). Q2 (`mSolverID`
read somewhere between free and next initialize) — NOT found. Q3 (the exhaustion-arm barrier
deadlocks against a rank on the success arm), which I named as the step I was least sure of — not
raised as a defect; both auditors instead attacked the *comment* that oversold it (D17). Q4
(deleting the setter is out of scope) — not raised as scope; Grok raised its consequence instead
(D23). **Q5 (something only a Fortran reader catches) — MISSED in the direction I expected and hit
in one I did not: the worst finding was in a file I had explicitly told the auditors was not mine.**

### 4.0 Audit Round 1 — plan jury, 2026-08-29 (Codex + Grok, blind/parallel)

Codex: **revise before implementation**. Grok: **approve with the required corrections**. Both
confirmed G1-G7, F1, F2 and the "G3+G4 are one write at successful create" reduction against the
tree. Neither auditor modified source (`git status` checked; the only file that moved during the
round was `src/numerics/sources/cl_SourceFunction.hpp`, a co-running session's work, untouched here).

| # | Defect in the PLAN | Severity | Found by | Disposition |
|---|---|---|---|---|
| D1 | R6's bound check was ordered ahead of the `allocated` guard — `size()` on an unallocated allocatable is undefined | HIGH | Codex + Grok (independent) | ACCEPTED, R6 rewritten with the explicit three-step order |
| D2 | R4's barrier move does not make the exhaustion *decision* collective; it relocates the hang into `DMUMPS` | HIGH | Codex + Grok (independent) | ACCEPTED, R4 rewritten; the allreduce escalation split out as O4 |
| D3 | **G9**: `aInfo` is `integer( int_t )` and is passed as `MPI_BARRIER`'s status; wrong kind under `BELFEM_INT64`, and clobbered at `:107` | HIGH | Grok | ACCEPTED as a new gap, verified against `mumpstools.f90:64` vs `:127` |
| D4 | **G8**: a second `initialize()` creates a slot and overwrites `mSolverID` before `Wrapper::initialize()` throws — both slots leak | HIGH | Grok | ACCEPTED as a new gap, R11 added |
| D5 | v1's O2 recommendation (key `free()` on `mSolverID > 0`) enables a **cross-wrapper destroy** if the ID is not zeroed | HIGH | Grok (refuting), Codex (proposing the same endpoint but WITH zeroing) | ACCEPTED. R10 added. O2 resolved: keep the boolean this round, zero the ID; Codex's single-handle consolidation logged as follow-up |
| D6 | R2 shipped alone would mark a G1-failed create as owing `JOB = -2`; it must land with R3+R5 | MEDIUM | Grok | ACCEPTED, sequencing constraint added to R2/R3/R5 |
| D7 | R1's sentence was factorization-specific and would mislabel pool exhaustion; no explicit `Infeasible` arm, so the next enumerator repeats G6 | MEDIUM | Grok | ACCEPTED, R1 rewritten |
| D8 | R9's "eighth cycle exhausts" is off by one without a pinned occupant | MEDIUM | Codex + Grok (independent) | ACCEPTED, R9 corrected to cycle 9 / `k < 12` |
| D9 | O1's "silent no-op" recommendation conflicts with the explicit-over-implicit error policy | HIGH (Codex) / endorsed-but-louder (Grok) | Codex + Grok | ACCEPTED, O1 resolved to (a)-**loud** or (b) |
| D10 | R7 is unmapped scope expansion relative to "DR-144 in full" | MEDIUM | Codex | PARTIALLY ACCEPTED — kept, but relabelled as by-catch and flagged for Christian rather than presented as a DR-144 step |
| D11 | Citation nits: F3 `:177-193`→`:180-199`; G6 ternary starts `:1561` not `:1563`; `is_initialized()` `:382-386`; G7 "declaration-only" is wrong (defined in Fortran, no C++ caller) | LOW | Grok | ACCEPTED, all four verified and corrected in place |
| D12 | v1 called `gNumSolvers` a counter kept in step with the table; it is **write-only** (`:32,76,94,172`, zero reads), so R3's decrement is bookkeeping, not authority | LOW | Grok | ACCEPTED, noted in R3 |

**Claude's pre-registration scored:** P1 (a live initialize-then-free path exists) — NOT found by
either auditor, and independently narrowed by Claude before the round returned: G3 is unreachable
through `Solver::solve`, which always solves after its lazy init (`cl_Solver.cpp:159-182`), and
reachable only through the public `wrapper()` back door. Grok added the case Claude missed — a
first-solve `BELFEM_ERROR` throw at `:629-632` fires *before* today's `:643` write, so a caught
error after a successful create also leaks. **P2 (the barrier move is not equivalent) HIT, by both
auditors, and was the round's most valuable finding.** P3 (an auditor prefers `mSolverID > 0`) HIT —
Codex proposed exactly that, and Grok refuted the unguarded form; the split is why R10 exists.
P4 (R5's sign / warning handling) partially hit — both noted `tInfo > 0` must not abort. **P5 (a gap
I did not file) HIT TWICE**, G8 and G9.

## 5. Open Design Questions

- **O1 — what happens to `mumpstools_set_num_max_solvers`?** Declared at `mumpstools.hpp:28`,
  defined at `mumpstools.f90:41-48`, called from nowhere in `src/` or `tests/`.
  **RESOLVED 2026-08-29 by the jury round → (a), but LOUD, or (b).** v1 recommended a silent no-op
  once the pool is allocated. Codex rejected that as conflicting with the explicit-over-implicit
  error policy — an initialization-time precondition earns an always-active error, not a silent
  swallow — and Grok independently endorsed (a) only in its loud form, on the same standard the
  DR-140 symmetry guard was held to. Both would block on (c), the reallocating pool, as a feature
  with no requester. **Claude's revised recommendation: (a)-loud** — no-op when the requested size
  matches, hard error when it differs after allocation. (b), deleting the API outright, is equally
  honest and Grok would not block on either; still Christian's call.
- **O2 — key `free()` on `mInitialized` or on `mSolverID > 0`?**
  **RESOLVED 2026-08-29 → keep `mInitialized` this round, and zero `mSolverID` (R10).** The
  auditors split and then converged. Codex recommended collapsing to the single ID handle; Grok
  refuted the plan's unguarded form of exactly that with a concrete cross-wrapper destroy (§4, D5)
  — but Codex's own proposal already included zeroing the ID, so the two agree on the load-bearing
  part and differ only on whether the boolean then goes away. Deleting it is a strictly larger edit
  than this round was scoped for, and is safe *only after* R10 lands.
  **Follow-up, not this round:** once R10 is in, `mInitialized` is redundant and Codex's
  single-handle consolidation becomes the tidier resting state. File it when R8/R9 are green.
- **O3 — should a failed create clear `mSolverID` explicitly?**
  **RESOLVED 2026-08-29 → yes, and widened into R10.** Grok: the pre-call store is insurance (the
  shim always writes the out-arg today at `:83-88`); the **post-free** store is the one that stops
  the cross-wrapper destroy. Both are in R10.
- **O5 — NEW, raised by the code round. The pool is now hard-capped at 8 with no setter.**
  Resolving O1 to (b) removed even the broken lever. Two instances per field (production +
  diagnostic) fits comfortably today, but a coupled deck with more fields each holding a diagnostic
  instance could reach the cap, and there is now nothing to raise it with. If that day comes the
  answer is a reallocating pool, which is the (c) both auditors would have blocked on as premature.
  Logged so the deletion is not mistaken for a claim that 8 is always enough.
- **O6 — NEW. Should COND1 get the same emptiness guard as COND2?** The MUMPS structure is
  identical: `COND(1)` is initialized to 1.0 and its estimator is skipped when no row lands in the
  first category. In practice omega1 is nonzero in every step observed, so the sentinel has never
  been seen — which is exactly why the guard was not added: it would be untestable. Logged rather
  than decided, because "we have never seen it" is a weaker argument than it sounds.
- **O4 — raised by the plan round. Should mixed-occupancy across ranks become a loud error?**
  R4 as corrected makes both *agreeing* create outcomes collective, which is hygiene. It does not
  detect ranks that disagree about which slot is free. Codex's shape: allreduce the chosen slot /
  exhaustion flag and enter `DMUMPS` only on consensus. Grok rates "not constructible from any
  current in-tree caller" at high confidence, and "stays that way under a future caller that skips
  a rank" at **low**. This is a real hardening step and a real cost (a collective on every create).
  Logged, not decided — and deliberately NOT smuggled into R4, which is what v1 effectively did.

### 5.1 Open risks carried forward (not steps, not silently dropped)

- `mumpstools_free_solver` has the **same barrier asymmetry** as create: the unallocated and
  unoccupied paths skip the barrier (`:138-146`). After R2 every rank that created will free, so it
  stays rank-uniform on in-tree paths. Residual, not a new gap.
- `mumpstools_solve` still indexes `gSolvers( tSolverID )` with no bound check (`:289`); the C++
  `mSolverID <= 0` gates (`:501`, `:661`) are the only fence. **Out of scope — R6 does not cover
  it**, and the plan must not imply otherwise.
- `MPI_COMM_WORLD` vs `gComm` throughout the shim is pre-existing and out of scope.

## 6. Definition-of-Done Checklist

- [x] Every gap-table row mapped to a step (G1→R3, G2→R5, G3+G4→R2, G5→R4, G6→R1, G7→R6,
      G8→R11, G9→R4b, G10→R12, G11→R13, G12→R14) or an open question (O1, O4, O5, O6). G10 and G11 are
      DEFERRED by ownership, not unmapped.
- [ ] Each claimed gap backed by a citation re-verified against the working tree, not against the
      register row that filed it.
- [x] Open questions logged, not silently decided. (O1 resolved to (b); O2/O3 resolved into R10;
      O4 and O5 left open.)
- [x] `make check` green (R8) and the initialize-then-free unit test present and passing (R9).
- [◐] DR-144 row updated: gates recorded as RUN, tag kept `[CODE]` `[W]` because G10 and G11 are
      live but latent residue. **NOT struck** — a residual living in another session's artifact is
      still a residual (register rule, "When to strike").
- [x] Devlog written; `todo/README.md` updated; `scripts/check_doc_claims.py` run after every
      register count change — it caught a wrong `[W]` count within a minute of my making it.
- [x] R2, R3, R5 and R10 land in ONE edit (D6 sequencing constraint).
- [x] A code-audit jury round on the resulting diff, with BOTH vendors, before the row is retagged.

## 7. Audit Trail

- Exchange thread: `tmp/ai_exchange/dr144_mumps_lifecycle.md` — Claude's pre-registration (written
  before dispatch), then the blind parallel Codex and Grok legs of 2026-08-29. Every `file:line` in
  the auditors' findings was re-verified against the working tree before it was accepted here; the
  four in D11 were corrections to this plan's own citations and were confirmed one by one.
- Register row: `todo/debt_register.md`, DR-144 (filed 2026-08-29 from the DR-143 audit rounds:
  Codex plan finding 1.4 plus Grok open risk; G4 from Codex's code-audit adjacent findings; G5-G7
  from Grok's code-audit leg).
