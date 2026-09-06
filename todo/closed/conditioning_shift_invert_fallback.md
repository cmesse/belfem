# Conditioning Diagnostic: Solver-Native First, Shift-Invert Fallback

> **CLOSED 2026-09-03** (todo/ currentness sweep, round 3): shift-invert is implemented and `make check` ran green 2026-08-29; the PETSc-native R7 belongs to the (cancelled) PETSc plan. Status lines and checkboxes below are as they stood at closure and are not maintained.

**Date:** 2026-08-28
**Purpose:** Implement Christian's ruling of 2026-08-28: the conditioning diagnostic reports the
metric the active solver can supply natively, and falls back to a **shift-invert** eigenvalue
computation when it cannot. Memory cost of the fallback is the user's responsibility — the
diagnostic is opt-in and off by default.
**Module:** `src/sparse` (wrappers + ARPACK drivers), `src/fem/kernel` (routing)
**AIs involved:** Claude (measurement + plan), Codex (audit), Grok (unavailable — CLI cannot start)
**Status:** **R1-R3 + R6 IMPLEMENTED AND THE PROBE GATE PASSES ( 2026-08-29 ).** Shift-invert
reproduces an analytic spectrum to 6.3e-12 on a matrix with the same dense small end that defeated
the fold, with 1 factorization and 60 frozen reuses. Six code defects found and fixed — four by the
Codex round, **one CRITICAL only by running it** ( §3.5 D1: the audit-recommended
`GeneralSymmetric` silently factorizes the wrong matrix ), one regression I introduced while fixing
another. **`make check` GREEN ( Christian, 2026-08-29 )** against a library rebuilt at 04:14, after the last
source edit at 04:11, with every new symbol present in `libbelfem.a` — so the wiring gate the probe
explicitly could not provide is now passed too. **R5/R8 REDESIGNED 2026-08-29 ( Christian's clarified ruling, second pass ): per TIMESTEP, not per
run.** `finalize_run()` and the once-per-run latches are REMOVED — on runs that get interrupted
( which is how this deck is actually operated ) an end-of-run hook never fires and the user gets the
number NEVER, which is worse than the honest per-step n/a it replaced. Now: a MUMPS field samples on
the FIRST iterate of each timestep ( native estimate, drift across the step negligible next to the
estimate's own error ); a non-native field runs the eigen shift-invert at the END of each timestep,
in `finalize()`, against the converged matrix, which is still intact after the last solve and is
the operator a user expects the number to describe. Same-session fix while wiring: `run_shift_invert`
dereferenced `mK` on every rank but `link_matrix()` is master-only, so a worker's `mK` was NULL —
a guaranteed parallel crash on first use, caught before it ever ran ( workers now use their own
assembly submatrix, exactly as the production solve does ). Force-eigen-under-MUMPS noted as future
work, out of scope. Still owed: the numeric gate against a REAL dumped Jacobian, R7 ( PETSc ), R9
( the input-contract docs ), and a rebuild + `make check` on the R5/R8 delta.
**DR-143 FOUND AND FIXED 2026-08-29 ( separate session, own round — `tmp/ai_exchange/mumps_pool_recycling.md` ):
the per-timestep R5/R8 shape exposed a latent mumpstools allocator defect — the instance registry
never recycled freed slots, so the per-step `free()`/re-init churn plus the pinned production MUMPS
exhausted the 8-slot pool at exactly timestep 8 of the tapestack3d deck ( all-rank abort in
`MUMPS::initialize` ). Fix: first-free scan in the registry, the shift-invert instance now
PERSISTS across timesteps ( the per-call `free()` at the end of `run_shift_invert` is removed —
also buys JOB 5 analysis reuse ), and a soft-failed create no longer aborts ( wrapper flag set
only after a successful create; both solve overloads bail flagged on `mSolverID <= 0` ). The
"diagnostic must not kill the run" promise now holds at the create site too. Gates owed there:
`make check` with the new `SolverLifecycle.MUMPSPoolRecyclesFreedSlots` test, deck rerun past
step 8. Code round closed with both vendors later the same day, no blocking defect. One
finding from it lands in THIS file's code and is filed as DR-144: `report_unavailable`'s
reason ternary has no `EigenOutcome::SolverFailed` arm, so a MUMPS soft-fail is reported
as "no legal Krylov subspace fits the basis budget" — the wrong cause, in the one message
a user reads when the diagnostic goes quiet.**
*(Plan history: v3 ( 2026-08-29 ) — v2 audit returned "revise before implementation" with the
mode-3 semantics verified against **arpack-ng 3.9.1 upstream**, the version actually linked ( §2.2 —
the audit's own citations were to a stale 3.9.0 tree, corrected by Christian ) and
two implementation-blocking defects plus a 25-item guess list; ALL folded into §2.3 below. The
planning session HALTS here per Christian's directive — a separate Fable session implements from
this file. §2.3 is the implementation specification; read it before touching code.)*
R4 ( symmetric drivers ) is DONE and landed. R7 ( PETSc singular values ) is explicitly sequenced
AFTER shift-invert is verified working, per Christian 2026-08-29. Grok unavailable throughout
( CLI sandbox cannot start ).

> **Scope guards:**
> - The diagnostic stays **opt-in** (`compute conditioning : true`) and off by default. Christian's
>   ruling 2026-08-28: "It is the user's responsibility to make sure the system has enough RAM."
>   `DR-77`'s memory finding is therefore a **documented caution, not a blocker**.
> - Track A (solver-native) is owned by `todo/conditioning_diagnostic_backends.md`. This plan owns
>   Track B (shift-invert) and the routing rule between them.
> - SLEPc stays deferred. This uses BELFEM's own ARPACK and `Solver`.

---

## 1. Why a fallback is needed at all — MEASURED, not argued

`sysdump_thermal_3.hdf5`, 2026-08-28, via `BELFEM_DUMP_SYSTEM=1` and
`tmp/analyse_thermal_spectrum.py`:

| quantity | value |
|---|---|
| n / nnz | 97095 / 1622509 |
| symmetry | `\|\|A-A^T\|\|/\|\|A\|\|` = 3.5e-16 — **symmetric** |
| lambda_min / lambda_max | 6.175959e-08 / 1.830271 — **SPD** |
| kappa_2 | **2.9635e7** |
| small-end relative gap | **9.51e-10** |
| gap under shift-invert | 2.82e-2 (amplification 3.0e7) |

The six smallest eigenvalues agree to three significant figures. A polynomial (non-inverting)
Krylov method cannot isolate an eigenvalue out of that cluster at any `ncv`, `maxit` or tolerance —
the spectral fold implemented in `todo/cancelled/arpack_small_end_configuration.md` spent 572 restarts and
6303 products confirming it. **Shift-invert is not an optimization here, it is the only Krylov route
to this end of the spectrum.** Confidence: high — reproducer.

## 2. The blocker, re-verified 2026-08-28

Shift-invert needs `OP*x = A^{-1}x` in the reverse-communication slot, i.e. many solves against ONE
factorization. Neither wrapper can do that today:

| wrapper | behaviour on a repeated `solve()` | citation |
|---|---|---|
| MUMPS | JOB 6 on the first call, **JOB 5 (factorize+solve) on every later call** with the same matrix pointer. No JOB 3 (solve-only) path exists anywhere in the file | `cl_SolverMUMPS.cpp:291-301` |
| PETSc | re-pushes values (`petsctools_update_matrix`) and re-runs `set_preconditioner` / `set_krylovmethod` on **every** `solve()`, so the PC is set up again | `cl_SolverPETSC.cpp:103-131` |

So ~20 Arnoldi products would cost ~20 factorizations. This is the "surgery on the wrappers every
production solve depends on" that deferred `todo/deferred/slepc_eigensolver_integration.md`, and
Christian's ruling now sanctions paying for it — but the change must be **inert for production**:
a flag that defaults to off, leaving JOB 5/6 exactly as they are today.

## 2.1 Architecture ( v2, decided 2026-08-29 — O1 RESOLVED )

**The inverse loop lives in C++, runs ARPACK on the MASTER only, and applies `A^{-1}` through a
collective `Solver::solve`. PARPACK is not involved.**

Why this shape wins ( and why the audit's PARPACK-bridge worry disappears ):

- The Krylov basis at `nev = 1`, `ncv ~ 10` on `n = 97095` is ~8 MB of doubles — there is nothing
  to distribute. The expensive object is the FACTORIZATION, and that is already distributed inside
  MUMPS.
- The collective-solve-with-master-held-matrix pattern is BELFEM's production pattern: workers
  enter `mumpstools_solve` with NULL data pointers ( `cl_SolverMUMPS.cpp:346-352` ), and
  `EigenValues` itself already owns a dedicated `Solver` that every rank calls collectively —
  `compute_lambda_max` does exactly this ( `cl_FEM_DofMgr_EigenValues.cpp:1238,1275` ). The
  implementation extends a precedent instead of inventing a bridge.
- The v1 audit judged C++-owned reverse communication "the safer BELFEM design" and warned that a
  callback-into-Fortran needs an explicit context object to stay thread-sane. The C++ loop needs
  neither callback nor context.

**The collective protocol** ( the part the v2 audit must attack hardest ):

```
every rank:                                  master only:
  loop:                                        dsaupd step ( via bind(c) shim )
    broadcast( tCommand )   <---------------   decide: 1 = apply OP, 0 = done
    if tCommand == 0: break
    mSolver->solve( *mK, tY, tX )   ← COLLECTIVE, all ranks, every OP application
    ( master copies tX/tY into workd slots; workers' vectors are carried but unread )
```

One broadcast per OP application, ~10-20 applications total, once per run. Every rank executes the
same solve sequence, so no rank can starve a collective. The command travels by `broadcast`, which
is collective — never by `share`/`receive`.

**ARPACK entry**: NOT direct `dsaupd_` externs — the tree has zero underscore externs and all
ARPACK access goes through `bind(c)` Fortran ( verified 2026-08-29 ), which keeps name-mangling
assumptions out of C++. Instead, two thin forwarding shims in `arpacktools.f90`:

```
arpack_si_step( ido, n, which, nev, tol, resid, ncv, v, ldv,
                iparam, ipntr, workd, workl, lworkl, info )   ! forwards to dsaupd
arpack_si_extract( ... )                                       ! forwards to dseupd
```

The shims forward verbatim — ALL state lives in the caller's arrays, which C++ owns as members of
`EigenValues`, sized once. Note dsaupd keeps internal SAVE'd locals between reverse-communication
calls of one factorization run: one eigenproblem at a time, process-wide. That matches BELFEM's
deliberate not-thread-safe posture, but it means the shims must never be driven from two
interleaved problems — say so in the shim header comment.

**Mode and transform** ( C1 applied ):

- `iparam(7) = 3` ( shift-invert ), `sigma = 0`, `which = 'LM'`, `bmat = 'I'` with the `ido = 2`
  slot answering `y = x` ( M = I ).
- The OP the solve applies is `A^{-1}` exactly; `dseupd` receives `sigma` and performs the back
  transformation ITSELF — the caller takes the returned values as eigenvalues of A and applies
  **no** `1/nu`. ( C1: applying it manually would invert twice. )
- **Verification item for the audit**: confirm from the dsaupd/dseupd source that mode 3 with
  `bmat = 'I'` is legal for the standard problem ( scipy's eigsh drives it exactly this way, but
  the plan should not rest on scipy's word ), and confirm which slot ( ido = -1 vs 1 ) must use
  `workd(ipntr(3))` in mode 3 — the ARPACK docs distinguish them and getting it wrong corrupts
  the basis silently.
- `A` symmetric ⟹ `A^{-1}` symmetric, so the SYMMETRIC pair ( dsaupd/dseupd, landed in R4 ) drives
  the loop. The nonsymmetric magnetic system is OUT of scope for the fallback until someone asks.

**The factorization** ( C3/C6 applied — scoped, dedicated, frozen ):

- A DEDICATED `Solver( SolverType::MUMPS )` owned by `EigenValues`, created for the diagnostic and
  freed after it — NOT the production thermal solver, which is PETSc in the reference deck and
  useless for this. `#ifndef BELFEM_MUMPS` ⟹ outcome `Infeasible`, honestly reported. ( STRUMPACK
  and UMFPACK reuse are NOT implemented in this round — one wrapper, done properly. )
- Call sequence inside one frozen scope, all owned by `compute_conditioning()`:
  1. first `mSolver->solve( *mK, y, x )` → wrapper runs JOB 6 ( analyze + factorize + solve )
  2. `mSolver->wrapper()->freeze_factorization()` — arms JOB 3
  3. the OP loop ( solves only, JOB 3 )
  4. `unfreeze + free` — always, including on every failure path ( RAII guard object, so an
     exception cannot leak an armed freeze )
- No assembly can intervene: the frozen interval is entirely inside `compute_conditioning()`, which
  must NOT call `compute_jacobian_and_rhs()` between steps 1 and 4 — and therefore `mMatrixFlag`
  handling moves BEFORE the freeze ( C6 ).
- λ_max comes first, from the existing exterior `'LM'` Lanczos run — no solve involved, so it runs
  BEFORE the factorization exists and needs nothing frozen.

**Memory**: one thermal MUMPS factorization at `n = 97095`, once, at end of run — when the
timestep loop is over and the magnetic factors are no longer live. Christian's ruling stands: RAM
is the user's responsibility; `DR-77` is the recorded caution, and its own amendment records a
coarse-deck thermal MUMPS run going for hours without workspace failure.

**Executable gate**: the measured ground truth. `lambda_min = 6.175959e-08`,
`lambda_max = 1.830271`, `kappa_2 = 2.9635e7` on `sysdump_thermal_3.hdf5`. A first implementation
gate that needs NO solver run: load that dump into a scratchpad probe, drive the new loop against
it serially, compare all three numbers. Only then wire it into the Controller.

## 2.2 Mode-3 Ground Rules ( verified against **arpack-ng 3.9.1**, the version actually linked )

> **PROVENANCE CORRECTION, 2026-08-29 ( Christian ).** The v2 audit verified these rules against
> `/opt/tpls_old/build/arpack/SRC/`. That tree **is** arpack-ng — not classic ARPACK — but it is
> **3.9.0** ( `configure.ac`: `AC_INIT([ARPACK-NG],[3.9.0]...)` ), while BELFEM links **3.9.1**:
> `/opt/scls/mkl/lib/libarpack.so.2.1.0` and `libparpack.so.2.1.0` both carry the build path
> `.../rpmbuild/BUILD/arpack-ng-3.9.1/`. The `tpls_old` path is a trap — a stale tree of a
> different version, sitting where a session naturally reaches for sources.
>
> **The rules below are therefore re-verified against upstream `opencollab/arpack-ng`, tag 3.9.1.**
> All five survived unchanged, so no conclusion moved — but the citations now point at the code
> that actually runs. **Methodological note for future sessions: a local source tree is not
> evidence about a linked binary until its version has been matched to the library.**

| rule | evidence, arpack-ng 3.9.1 |
|---|---|
| `bmat='I'` + mode 3 is legal for the standard problem | `SRC/dsaupd.f` validation rejects only `mode .eq. 1 .and. bmat .eq. 'G'` ( `ierr = -11` ); no mode-3/`'I'` rejection exists |
| on `IDO=1` in modes 3/4/5 the operand `B*X` is at `WORKD(IPNTR(3))` | `dsaupd.f` IDO=1 doc: *"In mode 3,4 and 5, the vector B * X is already available in WORKD(ipntr(3)). It does not need to be recomputed in forming OP * X."* |
| `IPNTR(3)` is defined as the `B*X` pointer in shift-and-invert | `dsaupd.f`: *"IPNTR(3): pointer to the vector B * X in WORKD when used in the shift-and-invert mode."* |
| `IDO=-1` uses `IPNTR(1)` ( initialization, forces the start vector into range of OP ) | `dsaupd.f` IDO=-1 doc |
| `IDO=3` is the user-shift request only | `dsaupd.f`: shifts placed via `IPNTR(11)`; unreachable at `IPARAM(1)=1` |
| `LWORKL >= NCV**2 + 8*NCV`, `IPNTR` length 11 | `dsaupd.f` argument docs |
| `dseupd` returns eigenvalues of the ORIGINAL problem | `dseupd.f` D doc: *"D contains the Ritz value approximations to the eigenvalues of A*z = lambda*B*z"*, ascending |
| the back-transform is applied INSIDE dseupd ( C1 ) | `dseupd.f`, `TYPE='SHIFTI'`: `workl(ihd+k-1) = one / workl(ihd+k-1) + sigma` |
| dsaupd state must be untouched between the calls | `dseupd.f`: *"These arguments MUST NOT BE MODIFIED between the the last call to DSAUPD and the call to DSEUPD."* |
| `RVEC=.FALSE.` ⟹ `Z` not referenced | `dseupd.f` |

**Two nuances the re-verification added:**

- dseupd's own `LWORKL` check is **guarded by `rvec`** ( `if (rvec .and. lworkl .lt. ncv**2+8*ncv) ierr = -7` ), so a values-only call would not catch an undersized `workl`. Size it correctly anyway — **dsaupd** requires the same bound unconditionally, and relying on an unchecked path is how a heap corruption gets found in production instead of at the call.
- **`IDO=2` stays a defensive branch, not an expected event** ( unchanged from v2 ): with `bmat='I'` the operand is copied internally rather than requested. Handle it; never rely on seeing it.

**Version-gap finding with a real consequence.** arpack-ng 3.9.1's CHANGES carry
*"[BUG FIX] Ensure that LAPACK RNG state is propagated (regression in 3.9.0)"* and *"[BUG FIX]
Ensure that separate random seeds are used on different parallel thread in D and S versions"*. That
is precisely the random-starting-vector path an `INFO=0` start uses — so the stale 3.9.0 tree and
the linked 3.9.1 differ exactly where this plan would have relied on them being identical.
**Decision for the implementer ( flagged, not pre-empted ):** `INFO=0` ( random start ) makes the
diagnostic non-reproducible run to run, which is a poor property for a number a user compares
across runs. `INFO=1` with a fixed starting vector ( e.g. normalized all-ones ) makes it
deterministic, at the documented risk that a pathological start misses an eigenvalue — low risk for
a well-separated shift-invert problem, and the probe gate would catch it. Choose deliberately and
record the choice.

scipy vendors ARPACK-ng and drives `eigsh( A, sigma=... , M=None )` through exactly the
configuration this plan proposes — read `scipy/sparse/linalg/_eigen/arpack/arpack.py:506-516`
( mode selection ), `:583-612` ( the symmetric reverse-communication loop ), `:560-582` ( sizing ).
What it establishes, and what it cannot:

| rule | evidence |
|---|---|
| `mode = 3` with `bmat = 'I'` is legal for the standard problem | scipy sets exactly this pair when `M is None` ( `:512-516` ) |
| on `ido = 1` in mode 3, the OP input is `workd( ipntr(3) )`, NOT `ipntr(1)` | `:601-603` — THE classic trap; reading `ipntr(1)` corrupts the basis silently |
| `lworkl = ncv*( ncv + 8 )`, `ipntr` dimension 11 — same as our landed symmetric driver | `:569,:581` |
| `iparam(1) = 1` ( exact shifts ) keeps `ido = 3` from ever occurring | scipy raises on `ido = 3` and never sees it |
| `dseupd` receives `sigma` and returns eigenvalues of A — no manual `1/nu` ( C1 ) | scipy applies no post-inversion |

**Robust ido rule for the implementation** ( correct under both classic ARPACK and ARPACK-ng ):

```
ido == -1  ->  input at ipntr(1)     ( startup: B*x not yet available )
ido ==  1  ->  input at ipntr(3)     ( B*x precomputed; equals x for bmat = 'I' )
both       ->  result to ipntr(2), y = A^{-1} * input
ido ==  2  ->  y = x                 ( B = I )
```

**Version-sensitivity flag**: scipy's current loop has NO `ido = -1` branch at all, implying the
ARPACK-ng it vendors may not issue `-1` in this mode — but BELFEM links its own ARPACK-ng build,
and the plan must not assume the two behave identically. The rule above is safe under EITHER
convention, which is why it is stated as the implementation rule rather than the minimal one.
The audit is still asked to confirm against dsaupd's own documentation ( v2 audit question 1 ).

## 2.3 Implementation Specification ( v3 — the v2 audit's blocking defects and guess list, resolved )

**The two implementation-blockers, decided:**

1. **Freeze signature and worker records.** `freeze_factorization( const SpMatrix & aMatrix )` —
   it takes the matrix. Rank split: the MASTER records and validates the full identity ( pointer,
   `n_rows`, `n_cols`, `nnz`, `data()`, pointers/indices pointers ); WORKERS record only a
   rank-uniform frozen flag — their MUMPS-side `mMatrix` stays null by design
   ( `cl_SolverMUMPS.cpp:283` branch is master-only ), so a worker CANNOT build the structural
   record and must not pretend to. Validation runs where the data lives; uniformity comes from the
   freeze being armed/dropped at the same point of the same code path on every rank. **Resolve the
   `mMatrix` shadowing while in there**: `MUMPS` re-declares `Wrapper`'s `mMatrix`
   ( `cl_SolverMUMPS.hpp:44` vs `cl_SolverWrapper.hpp:73` ) — the frozen record goes into its OWN
   struct ( `FrozenFactorization` ) in the MUMPS wrapper, touching neither shadowed member.
2. **The first factorizing solve IS the `IDO=-1` service — there is NO dummy solve.** Sequence:
   dsaupd returns `IDO=-1` → broadcast command 1 → collective solve on `WORKD(IPNTR(1))` → that
   solve runs JOB 6 ( analyze + factorize + solve ) → on success, freeze immediately → result to
   `WORKD(IPNTR(2))` → resume dsaupd. Every later `IDO=±1` request is a frozen JOB 3 solve.

**The RCI switch ( normative ):**

| dsaupd returns | action |
|---|---|
| `IDO = -1` | input `WORKD(IPNTR(1))` → collective solve → output `WORKD(IPNTR(2))` |
| `IDO = 1` | input **`WORKD(IPNTR(3))`** → collective solve → output `WORKD(IPNTR(2))` |
| `IDO = 2` | copy `IPNTR(1)` → `IPNTR(2)` ( defensive; not expected for `bmat='I'` ) |
| `IDO = 99` | broadcast command 0, exit loop |
| `IDO = 3` or anything else | broadcast command 0 so workers exit, then always-active `BELFEM_ERROR` — `IPARAM(1)=1` makes 3 unreachable, so reaching it is a driver bug |

Workers accept commands 0 and 1 ONLY; an unknown command word is an always-active error, not a skip.

**RCI initialization ( normative ):** `IDO=0`, `INFO` per the start-vector decision in §2.2
( `0` = random, `1` = caller-supplied `resid` — decide and record ), `IPARAM(1)=1`,
`IPARAM(3)=maxit`, `IPARAM(7)=3`, `LDV=N`, `LWORKL=NCV*(NCV+8)`, `TOL` from `configure()`,
`WHICH='LM'`, `BMAT='I'`, `sigma=0.0` handed to the extract shim.

**Master-side ARPACK members** ( sized once in `run_shift_invert`, `BELFEM_QUIET_NAN`-filled where
real ): `mSiResid(N)`, `mSiBasis(N*NCV)`, `mSiWorkD(3*N)`, `mSiWorkL(NCV*(NCV+8))`,
`mSiSelect(NCV)`, `mSiD(NCV)`, `iparam(11)`, `ipntr(11)`, dummy `z(1)` with `LDZ=1`.

**The shims are NOT verbatim forwards — they are BELFEM-style `bind(c)` wrappers** exposing ONLY
numeric arguments, exactly like `arpack_standard_eigen` does: `BMAT`, `WHICH`, `HOWMNY`, `RVEC`
and the `SELECT` logical array live INSIDE the Fortran shim as constants/locals, so no character
or Fortran-logical ABI ever crosses the C boundary ( this retires the v2 audit's ABI item 8
wholesale ). `arpack_si_step` forwards the numeric RCI state to dsaupd; `arpack_si_extract` calls
dseupd with `RVEC=.FALSE., HOWMNY='A'` hardcoded and copies `D` out. State between the last step
call and the extract call is untouched by construction, since C++ owns every array.

**Dedicated solver:** a NEW member `mShiftInvertSolver` — `mSolver` stays with
`compute_lambda_max()`, which uses it differently ( the §2.1 citation of `compute_lambda_max` is
precedent for collective ownership, NOT the same code path ). Guard `#ifdef BELFEM_MUMPS` BEFORE
construction on every rank ( `Solver( MUMPS )` throws in a MUMPS-less build,
`cl_Solver.cpp:72` ) — a MUMPS-less build reports `Infeasible` without constructing anything.
Symmetry mode: **`GeneralSymmetric`** via `set_symmetry_mode()` before the first solve — exploits
the declared symmetry without ASSUMING positive definiteness; one measured SPD dump is not an
invariant for every end-of-run Jacobian, and an indefinite matrix then still factors, with the
returned `lambda_min <= 0` mapping to `NotPositiveDefinite` honestly. The dedicated solver
inherits NOTHING from the production solver — fresh `SolverParameters`, no conditioning/error
analysis armed ( `ICNTL(11)` stays 0 in a fresh wrapper, `cl_SolverMUMPS.cpp:33` ). `ICNTL(10)=20`
( iterative refinement, set at init, `cl_SolverMUMPS.cpp:168` ) is DELIBERATELY KEPT: it runs per
inverse application and buys OP accuracy, which is exactly what a 1e-10 eigen tolerance wants.

**Outcome mapping ( extends `EigenOutcome` ):** add `SolverFailed` — a MUMPS factorization or
solve failure is neither `NotConverged` ( ARPACK's fault ) nor `Infeasible` ( size/build );
conflating them would misdirect whoever reads the message. Mapping: JOB 6 fails → `SolverFailed`
( latch: this run is over — end-of-run has no retry ); frozen JOB 3 fails → invalidate freeze,
`SolverFailed`, same latch; dsaupd `INFO=1` / `NCONV<NEV` → `NotConverged`; extract failure →
`NotConverged`; `lambda_min <= 0` → `NotPositiveDefinite`. **Latch semantics change with
once-per-run:** there is no later timestep to spend strikes in, so in the end-of-run path EVERY
non-Ok outcome latches immediately and `mMaxFailures` is inert ( keep the counter for any caller
that still runs per-step; the once-per-run caller just never benefits ).

**Result synchronization ( normative ):** after the loop, master broadcasts — outcome code first,
then `lambda_min`, `NCONV`, `NUMOP`, restart count. Workers read nothing from their stale local
arrays. Same discipline as the fold path ( broadcast-outcome-before-branch ), unchanged.

**DR-127 sequencing:** R1-R3 plus the scratchpad gate are INDEPENDENT of the Controller and may
land while DR-127 is in flight. R5/R8 ( `finalize_run()` wiring ) block until the DR-127 session
lands its restructure — "landed" means its Controller changes are committed and it has released
the file, per the standing cross-session agreement ( §4.1 ). Do not wire the Controller against
the pre-DR-127 shape.

**Scratchpad gate ( normative, before any Controller wiring ):** serial probe per the
scratchpad-probe pattern ( standalone TU linking the prebuilt `.a`s, `-DDEBUG` for header asserts,
ABSOLUTE paths ): load `cmake-build-debug/tapestack3d/sysdump_thermal_3.hdf5` via
`SpMatrix::load`, run the new loop single-rank, require `lambda_max` rel-err `< 1e-4` vs
`1.830271`, `lambda_min` rel-err `< 1e-3` vs `6.175959e-08`, `kappa_2` within a factor 1.01 of
`2.9635e7` ( ground truth measured at scipy tol 1e-10; the looser lambda_min band absorbs the
cancellation-floor arithmetic of §3.3 in the parent plan ). The probe is a probe — it cannot
verify build wiring; `make check` remains the wiring gate.

## 3. Ordered Steps

- [x] **R1 ( v3 — implement per §2.3, which OVERRIDES the v2 wording on the signature and the
      worker record )** — `solver::Wrapper`: scoped freeze API, not a trusted boolean ( C3 ).
      `freeze_factorization()` / `unfreeze_factorization()` virtuals ( default: `BELFEM_ERROR`,
      "this wrapper cannot reuse a factorization" ) plus `supports_factorization_reuse()`
      ( default `false` ). `freeze_factorization()` may only be called after a successful
      factorization and RECORDS the frozen state: matrix pointer, `n_rows`, `n_cols`, `nnz`,
      `data()` pointer, pointers/indices pointers. Every frozen-mode solve re-checks the record
      with **always-active `BELFEM_ERROR`** — a debug assert is the wrong tier for silent
      stale-factor use ( C3 ). No values scan per solve; the contract is scoped exclusivity, with
      an optional `BELFEM_ASSERT`-tier checksum of `data()` at freeze and at unfreeze to catch an
      accidental violation in debug builds. Production paths never call any of it — their JOB
      sequences are bit-identical, and the v2 audit should confirm that claim against the diff.
- [x] **R2 ( v2 )** — MUMPS: frozen ⟹ **JOB 3** ( solve only ). Applies to **both** overloads —
      the vector and the multiple-RHS JOB-selection blocks ( C4, `cl_SolverMUMPS.cpp:281,432` ).
      The freeze is invalidated ( forced unfrozen + record cleared ) in BOTH soft-failure paths
      that already clear `mMatrix` ( C5, `:373,:547` ) and in `free()`. Frozen-state agreement is
      collective by construction: the freeze is armed and dropped at the same point of the same
      code path on every rank, never from a rank-local value. dsaupd's own JOB-3 prerequisite —
      factors surviving in the saved `DMUMPS_STRUC` — was verified by the v1 audit
      ( `mumpstools.f90:29,:261,:125` ). (after: R1)
- [x] **R3 ( v3 — implement per §2.2 + §2.3; the tables there are normative )** — the C++
      shift-invert loop, per §2.1: `bind(c)` shims `arpack_si_step` /
      `arpack_si_extract` in `arpacktools.f90` ( verbatim forwards to dsaupd / dseupd );
      `EigenValues::run_shift_invert()` owning all reverse-communication arrays as members;
      master-only ARPACK, broadcast command word, collective dedicated-MUMPS solves; mode 3,
      `sigma = 0`, `bmat = 'I'`, `which = 'LM'`; **values taken from dseupd as returned — no manual
      `1/nu`** ( C1 ); RAII freeze guard ( C3/C6 ); outcome codes and the latch reused unchanged
      from the fold machinery. Gate BEFORE Controller wiring: the scratchpad probe against
      `sysdump_thermal_3.hdf5` must reproduce `lambda_min = 6.175959e-08` and
      `kappa_2 = 2.9635e7`. (after: R2)
- [x] **R4 — DONE 2026-08-28** ( Christian: "for the thermal problem we can assume the matrix is
      symmetric, but this is not the case for the magnetic problem. You can implement the missing
      functions in the same fashion in which I have implemented the other ones." ) Landed:
      `arpack_symmetric_eigen` ( `arpacktools.f90` ) and `parpack_symmetric_eigen`
      ( `parpacktools.f90` ) driving `dsaupd`/`dseupd` and `pdsaupd`/`pdseupd`, with the same
      argument list, `info` layout and fold argument as the nonsymmetric pair so the call site
      switches on the matrix alone; separate `check_saupd` / `check_seupd` decoders, because the
      symmetric error tables genuinely differ ( `-5` admits `'LA'`/`'SA'`/`'BE'`, `-13` means
      something else, there is no `info = 3`, and `dseupd` carries `-15`/`-16`/`-17` that `dneupd`
      does not ); three real size differences carried through — `workl` is `ncv*( ncv + 8 )` not
      `ncv*( 3*ncv + 8 )`, `ipntr` is 11 not 14, and one real `d` array replaces `dr`/`di` with
      `lambdaimag` zeroed. Selection is **declared, not detected**: `set_symmetric( true )` on the
      thermal field in `Controller::arm_conditioning_thermal()`, magnetic keeps the nonsymmetric
      default, which is the safe default ( the nonsymmetric driver is correct for a symmetric
      matrix, merely slower ). All four TUs syntax-check clean; **nothing built or run**.
      *(Superseded text:)* ~~Use the SYMMETRIC driver where the matrix is symmetric.~~ The measurement says the
      thermal Jacobian is symmetric to 3.5e-16, and `src/sparse` contains **no** `dsaupd` /
      `dseupd` / `pdsaupd` — only the nonsymmetric path. Lanczos is cheaper (short recurrence, far
      less re-orthogonalization) and for a symmetric matrix the spectral ratio IS `kappa_2`, which
      is the number the footer claims to print. Detect symmetry once per matrix rather than
      trusting the physics. (after: R3, and separable from it)
- [ ] **R5** — **WHEN each metric is taken ( Christian's ruling, 2026-08-28, second pass ).** The
      analysis does not belong in the per-solve path at all:
      - **solver-native: sample at the FIRST step of the RUN, then stop.** It is a by-product of a
        solve that happens anyway, but arming it is not free — `ICNTL(11) = 1` runs iterative
        refinement plus condition estimation on every armed solve. Today it is re-armed EVERY
        timestep ( `cl_FEM_Controller.cpp:211` in `initialize_timestep`, `:344` in the magnetic
        path ) and captured every timestep. One sample per run replaces that.
      - **eigen / shift-invert: ONCE, at the END of the run.** "We do not know a priori which step
        will be the last one" — so it is not scheduled, it is hooked. This is the step that makes
        the whole cost problem disappear: the 13223 ms measured per failing step becomes 13 s once.
      Routing between them stays as the ruling states: native where the solver supplies it,
      shift-invert where it does not, an honest statement where neither is possible, and the
      quantity NAMED wherever it is printed ( O2 ). (after: R4)
- [x] **R8 — the end-of-run hook does not exist yet.** `Controller::finalize()`
      ( `cl_FEM_Controller.cpp:2647` ) is per-TIMESTEP, not per-run: it stops the iteration timer
      and prints the step footer, and it is where `compute_conditioning()` is called today
      ( `:2655` ). A new `Controller::finalize_run()` is needed. `belfem.cpp` has TWO time loops —
      coupled ( `:181-220` ) and segregated ( `:222-282` ) — and both fall through to a single
      point before `return gComm.finalize()` ( `:284` ), so ONE call there covers both. Not the
      destructor: this can throw, and a solver failure inside `~Controller` would be unrecoverable.
      (after: R5)
- [x] **R6 ( v2 )** — Retire the fold path inside `compute_conditioning()`: shift-invert replaces
      it as the small-end route wherever MUMPS exists, and where MUMPS does not exist the fold
      cannot converge on this operator class anyway ( measured ). The fold ARGUMENT stays in the
      Fortran drivers — it costs nothing there and the sigma parameter is load-bearing for the
      shims' argument compatibility. O3 thereby resolves to "retire the path, keep the plumbing".
      ~~Retire or cap the spectral fold.~~ It has no remaining niche: wherever a factorization
      exists shift-invert dominates it, and where none exists it cannot converge on this operator
      class anyway. **O3.** (after: R5)
- [ ] **R7 ( re-scoped 2026-08-29, Christian: "when that works, 2) if petsc is available" )** —
      PETSc `KSPSetComputeSingularValues` + `KSPComputeExtremeSingularValues`, **strictly after R3
      is verified against the ground truth**, gated at runtime on `BELFEM_PETSC` and a GMRES-class
      KSP, reported under its own name ( preconditioned-operator ratio — NOT kappa_2, see O2 ).
- [ ] **R9** — `doc/input_file_reference.md` **and** `doc/input_schema.yaml` for any deck-visible
      change, plus the support matrix of which solver yields which quantity. The two-artifact rule
      applies the moment R5 touches a key. (after: R5)

## 3.1 One Caveat on Sampling Once — MEASURED

The ruling rests on "the conditioning will not change significantly". At order-of-magnitude
resolution that holds; more tightly it does not. From the run logs, the magnetic MUMPS number over
consecutive steps:

```
out2.txt   4.66e7 -> 5.10e7 -> 7.50e7 -> 7.95e7 -> 8.47e7 -> 9.85e7 -> 1.17e8
out.txt    8.23e7 ... 2.29e8 -> 2.68e8
```

A factor of about **6 across the observed steps**, drifting monotonically upward as the quench
develops. So a single first-step sample is a fair order-of-magnitude diagnostic and NOT a
description of the run.

**Consequence for the display, which is a real decision and not a detail:** if the first-step value
is printed in EVERY step's footer, a reader sees a constant 4.66e7 for a system that reached 2.68e8
and will conclude it is stable. Either print it once, or label it as the first-step sample. Folded
into O4.

## 3.2 Audit Corrections ( Codex, 2026-08-28 ) — status at v2: C1/C3/C4/C5/C6 are FOLDED INTO
R1-R3 ( v2 ) and §2.1; their boxes tick when the implementing session lands the step that carries
them. C2 is routed to Christian ( serial-PETSc stale-PC question — debt register, not this plan ).
C7 stays an open alternative, unchosen.

- [ ] **C1 (CRITICAL) — the plan's back-transformation is WRONG.** R3 said "`lambda = 1/nu`". With
      ARPACK **mode 3**, `dneupd` / `pdneupd` ( and `dseupd` / `pdseupd` ) **already** transform the
      Ritz values back to eigenvalues of the original problem. Applying the inversion again would
      invert them a second time. Whoever writes R3 takes the values as returned.
- [ ] **C2 (HIGH) — PETSc must be split into serial and distributed, and the serial case may be a
      LIVE BUG.** The plan claimed both paths refactorize. Distributed does — it redistributes and
      calls `petsctools_update_matrix()` per solve, which assembles specifically to announce changed
      values and force a PC rebuild. **Serial does not:** it shares the `SpMatrix` arrays directly,
      later solves update vectors only, and `petsctools_notify_matrix_update()` — which exists
      precisely to announce externally-changed shared values — **is called nowhere in the tree**. So
      serial PETSc may already be solving against a stale preconditioner after the values change.
      That is a production correctness question independent of this plan and belongs in the debt
      register, not here. **Route to Christian.**
- [ ] **C3 (HIGH) — a trusted reuse boolean can silently solve the wrong system.** If values change
      and MUMPS gets JOB 3 it solves with the old factors and returns a valid-LOOKING answer.
      Pointer identity proves only object identity — the wrapper already treats the same pointer as
      "same structure, different values" and deliberately picks JOB 5. A debug-only assert is the
      wrong tier for silent stale-factor use. Wanted instead: an explicitly scoped frozen state
      armed only after a successful factorization, recording pointer, dimensions, nnz, data and
      pattern pointers and factorization state, checked with always-active errors, disarmed after
      the diagnostic. An exact values check cannot be free — `memcmp` is O(nnz) per product — so the
      contract should be scoped exclusivity ( no assembly inside the frozen interval ) with an
      optional debug checksum at entry and exit.
- [ ] **C4 (MEDIUM) — both MUMPS overloads need the policy.** There are two independent JOB-selection
      blocks; the multiple-RHS overload repeats the same policy. Changing only the vector overload
      leaves inconsistent reuse.
- [ ] **C5 (MEDIUM) — factor validity needs invalidating in three places**, not one: both
      soft-failure paths that already clear `mMatrix`, and `free()`.
- [ ] **C6 (MEDIUM) — reuse must be tied to call timing.** `compute_conditioning()` runs in
      `finalize()` while the thermal capture happens right after the first thermal solve; and
      `EigenValues` can itself recompute the Jacobian when `mMatrixFlag` is set, which invalidates
      live factors although the `SpMatrix*` is unchanged. R3 must name which factorization is being
      reused and prove no reassembly intervened.
- [ ] **C7 (LOW, opportunity) — a cheaper SPD estimate exists.** A Hager/Higham 1-norm estimator
      needs only a few solves against existing factors and no Arnoldi basis at all ( for symmetric
      A the transpose solve is the same operation ). It yields an estimate of `kappa_1`, not
      `kappa_2`, so under the naming ruling it is honest — just not the spectral condition number.

## 3.3 Defects Found While Implementing R4

- [x] **D1 (MEDIUM) — the generated driver carried six FABRICATED citations.** The distributed
      symmetric driver was produced by transforming the nonsymmetric one, and the blanket
      `pdnaupd -> pdsaupd` rename silently rewrote verified `file:line` citations inside comments
      into claims about `pdsaupd.f:564`, `pdseupd.f:335`, `pdsaupd.f:522`, `dsaupd.f:504` and
      `pdsaupd.f:526` — line numbers that were checked for the NONSYMMETRIC sources and that
      nobody has ever checked for the symmetric ones. Caught by grepping the generated text for
      `.f:` before compiling. **Fixed 2026-08-28:** the numbers are stripped and the claims restated
      without them. Recorded because the mechanism generalizes — **any** mechanical rename across
      commented code can turn a verified citation into an invented one, and the compiler cannot
      see it.

## 3.4 Defects Found by the Implementation Audit ( Codex, 2026-08-29 ) — ALL FIXED

- [x] **I1 (CRITICAL) — `ido = 2` deadlocked the ranks.** The first cut branched on
      `tIsMaster && tIdo == 2`, but `tIdo` is master-only state: a worker's copy stays 0, so every
      worker fell through to the `else` and blocked inside the collective MUMPS solve while the
      master took its copy branch and went back to ARPACK. A collective decided from rank-local
      state. **Fixed:** the command word is now THREE-valued — Stop / Solve / IdentityCopy — and
      every rank branches on the broadcast, never on `tIdo`. ( Predicted by Claude in the audit
      brief and confirmed by Codex; the prediction is why it was written as a question rather than
      shipped as an assumption. )
- [x] **I2 (CRITICAL) — a second run would have selected JOB 5 against a destroyed instance.**
      `MUMPS::free()` runs JOB -2 and clears `mInitialized`, but did NOT clear `mMatrix`. On the
      next call the wrapper re-initializes, `select_job()` finds the same `SpMatrix` address still
      recorded and returns JOB 5 ( factorize + solve, reusing an analysis ) against a brand-new
      instance that has analysed nothing. **Fixed:** `free()` clears `mMatrix` so the next solve
      selects JOB 6, and clears the soft-failure latch, which is sticky and belongs to the
      instance being torn down.
- [x] **I3 (HIGH) — `SolverFailed` was unreachable, and the diagnostic could kill the run.** The
      dedicated solver never armed soft-fail, so a factorization failure raised an always-active
      error instead of setting `failed()` — the whole job would die to report a conditioning
      number. **Fixed:** `set_soft_fail( true )` plus `clear_failure()` on the dedicated wrapper on
      every rank before any solve.
- [x] **I4 (HIGH) — the frozen scope was not exception-safe.** The ordinary failure break did
      unfreeze, but a hard error inside a frozen solve, or `freeze_factorization()` throwing on its
      own identity check, would leave the scope armed over dead factors. **Fixed:** a
      `FrozenScopeGuard` armed BEFORE the freeze call, released on the normal path.

Confirmed correct by the same audit, and worth recording because they were the parts most likely to
be silently wrong: the `ipntr` operand slots ( `ido = 1` → `IPNTR(3)`, `ido = -1` → `IPNTR(1)`,
result → `IPNTR(2)` ), the whole `iparam` index map, the freeze identity check holding on the master
because `select_job()` sets `mMatrix` on the JOB 6 path, both solve overloads routing through
`select_job()`, both soft-failure paths and `free()` invalidating the frozen record, and the
extraction consuming `dseupd` output without a second reciprocal.

Open, not fixed: `mSiRhs` / `mSiLhs` are sized to the GLOBAL row count on every rank although MUMPS
ignores worker storage — correct but wasteful, two extra global vectors per worker. Left alone
deliberately: the alternative is rank-dependent sizing, which is exactly the kind of asymmetry that
produced I1.

## 3.5 THE GATE PASSED — and it caught a defect no audit round did ( 2026-08-29 )

Run in-session against a **1D Laplacian, n = 2000**, whose spectrum is known in closed form
( `lambda_k = 4 sin^2( k pi / 2( n + 1 ) )` ). Chosen deliberately over a synthetic well-separated
matrix: its small end is dense, relative gap 1.8e-6, so it has the SAME pathology as the thermal
operator and the retired fold would fail on it too.

```
  large end : 1 converged in 82 restarts, 2491 products
  solve  1 : ||Ax-b||/||b|| = 1.343e-12   ( JOB 6, fresh )
  solve  5 : ||Ax-b||/||b|| = 3.322e-12   ( JOB 3, frozen )
  small end : 1 converged in 1 restarts, 61 solves  ( info 0 / 0 )

  lambda_max   got  3.9999975e+00   ref  3.9999975e+00   rel 2.11e-15   PASS
  lambda_min   got  2.4649350e-06   ref  2.4649350e-06   rel 6.25e-12   PASS
  kappa_2      got  1.6227598e+06   ref  1.6227598e+06   rel 6.25e-12   PASS
```

MUMPS JOB census over the run: **1 x JOB 6, 60 x JOB 3.** One factorization, sixty reuses — the
freeze contract doing exactly what it was built for, and the small end converging in ONE restart
where the fold spent 572 and got nothing.

### D1 (CRITICAL) — `SymmetryMode::GeneralSymmetric` silently factorizes the WRONG matrix

The v2 plan audit recommended `GeneralSymmetric` for the dedicated solver, reasoning that the
matrix is symmetric and MUMPS could halve the factor memory. I implemented that. **It is wrong for
this codebase, and nothing static would have caught it.**

BELFEM stores the FULL matrix and passes `SymmetryMode` straight into MUMPS `SYM`
( `cl_SolverMUMPS.cpp:159-160` ); **nothing anywhere extracts a triangle** — verified by grep across
the wrapper and the Fortran shim. MUMPS with `SYM = 1` or `2` expects the lower triangle only, so a
full matrix double-counts every off-diagonal and the factorization is of a different matrix.

It does not fail. It returns a plausible wrong answer:

| | `SYM = 2` ( as recommended ) | `SYM = 0` ( corrected ) |
|---|---|---|
| first solve `\|\|Ax-b\|\|/\|\|b\|\|` | **7.4e+16** | 1.3e-12 |
| converged `lambda_min` | **-2.4967534e-19** | 2.4649350e-06 |
| ARPACK verdict | `info 0 / 0`, `nconv 1` | `info 0 / 0`, `nconv 1` |

ARPACK reported success in both columns. **Fixed 2026-08-29:** the dedicated solver uses
`SymmetryMode::Unsymmetric`, which is also what every production BELFEM solve uses. Do not
"optimize" it without first making the caller supply a triangle.

This is the entry that justifies the whole probe: three audit rounds read this code and none could
see it, because it is a property of how BELFEM feeds MUMPS rather than of the code on the page.

### Two probe-only traps, recorded because the next caller will hit them

- **The ARPACK drivers require a FORTRAN-based matrix** ( `link_matrix()` does this in production ).
  A zero-based matrix makes the matvec's `do j = pointers( 1 ), ...` start at zero and read one
  element before `values`, corrupting the heap; the abort then surfaces inside an unrelated
  `malloc` much later. The base must ALSO be restored before the matrix reaches MUMPS, which does
  its own `create_coo_indices()` and base switch.
- **Changing a class layout invalidates every TU that allocates it.** Adding the
  `FrozenFactorization` member to `MUMPS` changed `sizeof`, and `Solver::create_wrapper()` — which
  `new`s it — lives in another TU. Linking a fresh `cl_SolverMUMPS.o` against a stale
  `libbelfem.a` gave heap corruption inside the MUMPS constructor. A scratchpad probe must
  recompile every TU that sees the changed header, not just the changed one.

## 3.6 R5/R8 Landed 2026-08-29 — and a LATENT PRODUCTION DEFECT found alongside

`Controller::finalize_run()` runs the eigen fallback ONCE, after the time loops, wired into all
THREE executables that call `finalize()` — `belfem.cpp`, `hphirun.cpp`, `hphiTrun.cpp`. Wiring only
`belfem` would have silently removed the diagnostic from the other two, since the eigen call moved
OUT of the per-timestep `finalize()`.

The end-of-run block NAMES the quantity beside each number ( "MUMPS COND1, rhs-dependent" versus
"kappa_2, spectral" ). Without that, a COND1 of ~1e5 printed next to a kappa_2 of 2.96e7 for the
same operator reads as a bug.

**Arming changed from "first step of the run" to "armed until a capture actually happens"**, on the
DR-127 session's finding: `capture_conditioning_*` fires on a SOLVE, not a timestep, and since the
certified-exit restructure a timestep may legitimately take ZERO solves ( `min iterations : 0` with
a predictor already at tolerance — rare on a stock deck, whose default is 2, but a warm restart
lands close to converged by construction, which is where it is most likely ). First-step-only
arming would then capture nothing and the run would report no conditioning at all, silently. A NaN
from `get_cond0()` is now read as "not sampled yet" and the arm is held.

### DR candidate ( ROUTE TO CHRISTIAN ) — `IWG`'s symmetry default is SYM=1 on the production path

Found by the DR-127/DR-128 session while checking my diagnostic-path defect, and **verified
independently here**:

- `cl_FEM_DofMgr_SolverData.cpp:2056` pushes `mParent->iwg()->symmetry_mode()` into the solver on
  creation — so the mode reaches MUMPS on **every production solve**, not just from the eigen path.
- `cl_IWG.hpp:309` defaults that to `SymmetryMode::PositiveDefiniteSymmetric` — **MUMPS SYM = 1**,
  which expects the LOWER TRIANGLE ONLY, exactly as SYM = 2 does.
- BELFEM stores the FULL matrix and never extracts a triangle ( §3.5 D1 ), so any symmetric mode
  factorizes a different matrix, silently.

What saves the tree today is an accident of derivation: every concrete IWG derives from
`IWG_Timestep`, whose own ctor default is `Unsymmetric` ( `cl_IWG_Timestep.hpp:117` ) and shadows
the base. Verified: no site in `src/fem` reaches `IWG`'s default.

**So it is not currently exploitable, and it is one derivation away from being so.** A future IWG
deriving directly from `IWG` gets SYM = 1 and a silently wrong factorization on the production
path — strictly worse than the diagnostic-path instance I hit, which at least only corrupted a
diagnostic. Two candidate fixes, both Christian's call: extract a triangle in the wrapper when a
symmetric mode is set, or make the base default `Unsymmetric` and force symmetric to be explicit.
Not filed in `todo/debt_register.md` because a concurrent session holds that file.

## 3.7 R5/R8 Audit Findings — both fixed

- [x] **A1 (HIGH) — the run-level sample was erased on the next timestep.** `capture_*` stored the
      MUMPS number in `mConditionNumber*`, but those are the PER-STEP display values and are reset
      to NaN at the top of every timestep by design ( "stale-kappa hygiene", `:209`, `:345`,
      `:456` ) so a step that captured nothing prints `n/a` rather than the previous step's figure.
      With one-sample-per-run arming the latch then blocked recapture, so the value was erased on
      the step after the capture and never replaced — only the capturing step would have printed a
      number, and `finalize_run()` would have reported nothing at all. **Fixed:** `mRunConditionNumber*`
      hold the run sample and are never reset; `finalize_run()` reads those for a NATIVE field and
      `mConditionNumber*` for an EIGEN field ( whose value is computed after the last reset ).
      Found by the Codex round — note it is the same failure shape the DR-127 session predicted
      from a different direction, which is why the arming change was made at all.
- [x] **A2 (LOW) — my comment overstated `get_cond0()`'s contract.** It returns NaN when the
      analysis was not armed, but nothing guarantees an armed solve returns a finite estimate. The
      code treated NaN as "not sampled yet" and held the arm, so a degenerate estimate would have
      left `ICNTL(11)` on for the whole run — the per-solve cost this design exists to remove.
      **Fixed:** attempts are capped at `gMaxConditioningAttempts = 8`, after which the field gives
      up and disarms. Confidence that it was reachable: low, per the audit ; the cap costs nothing.

Confirmed clean by the same round: the collective placement in `finalize_run()` ( the rank-0 return
comes after every collective ), the sample latch's rank-consistency ( MUMPS broadcasts RINFOG to all
ranks, so the latches cannot diverge and `set_mumps_error_analysis` stays uniform — my main worry,
refuted ), no-diagnostic being a genuine no-op, and all three executables reaching the call on every
path with no early return after controller construction.

## 3.8 The MUMPS Symmetry Contract — hardened at the wrapper ( 2026-08-29 )

Christian ruled on the adjacent question: "make unsymmetric the default and symmetric explicit."
The DR-127 session took the IWG half ( `cl_IWG.hpp:322`, base default now `Unsymmetric` ). The
wrapper half is here, and it turned out to be broader than the IWG default:

| site | was | now |
|---|---|---|
| `cl_SolverWrapper.hpp:173` `initialize()` default | `GeneralSymmetric` | `Unsymmetric` |
| `cl_SolverMUMPS.hpp` `initialize()` default | `GeneralSymmetric` | `Unsymmetric` |
| `cl_SolverMUMPS.hpp` `mSymmetryMode` member | `GeneralSymmetric` | `Unsymmetric` |

The two `initialize()` defaults are a virtual pair and are kept IDENTICAL deliberately: default
arguments bind statically, so a base/derived mismatch would change the mode with the static type of
the pointer the call goes through.

**Plus an always-active guard**: `MUMPS::initialize()` now rejects any non-`Unsymmetric` mode with a
`BELFEM_ERROR` naming the reason. An error rather than a silent downgrade to SYM = 0, because a
caller that asks for symmetry wants the memory saving and quietly not giving it is its own kind of
lie. This converts the worst failure class ( silently wrong ) into the best ( loud and immediate ),
and it is the honest state of the code until someone implements triangle extraction.

Scope check that makes this safe: **only MUMPS consumes `mSymmetryMode`** — verified across every
wrapper; UMFPACK, PARDISO, PETSc, STRUMPACK and the base all have zero uses. SuperLU takes the enum
but maps symmetric modes to `SymmetricMode = YES`, which is an ORDERING heuristic, not a
triangle-storage contract, so it is unaffected. Nothing in-tree requests a symmetric mode for MUMPS.

## 4. Open Design Questions

- **O1 — RESOLVED 2026-08-29 → C++ master-loop + collective dedicated-MUMPS solves, no PARPACK,
  no callback** ( §2.1 ). The v1 audit's PARPACK local-vector bridge concern is mooted: PARPACK is
  not used for the inverse loop, because the basis is megabytes while the factorization — the only
  big object — is already distributed inside MUMPS. *(original question kept below)*
  ~~where does the reverse-communication loop live?~~ Moving it into C++ is the honest
  structure (the operator is a C++ object) but rewrites both drivers and changes the PARPACK
  collective story, since every rank must call the solve. A C function pointer into Fortran keeps
  the drivers intact but is harder to make MPI-correct. **Not decided — needs the audit's view.**
- **O2 — the footer prints one number under one name, and there are now three.** MUMPS COND1
  (Arioli/Demmel/Duff, componentwise), PETSc `KSPComputeExtremeSingularValues`
  (`sigma_max/sigma_min` of the **preconditioned** operator), and the spectral ratio from
  shift-invert (`kappa_2` only when symmetric). The 2026-08-28 measurement shows how far apart they
  are: MUMPS COND1 was quoted at ~1e4-1e5 for this system while the measured `kappa_2` is 2.96e7.
  **Reporting them under one label would be misleading.** Inherited unresolved from
  `conditioning_diagnostic_backends.md` O1. **Not decided — needs Christian.**
- **O4 — how is a once-per-run number displayed in a per-step footer?** See §3.1: the measured
  drift is ~6x, so a first-step value repeated in every footer misrepresents the run. Options:
  print it only in the step that sampled it; label the column as the first-step sample; or move
  both numbers to an end-of-run summary beside the eigen result. **Not decided — needs Christian**,
  because it is about what the output MEANS, not about mechanism.
- **O3 — retire the fold, or keep it capped?** It cost 12000 products to fail on the measured
  matrix. Keeping it as the no-factorization path means keeping a path that cannot work on the very
  operator class the diagnostic is aimed at. **Recommendation: retire it**, keep the shift-invert
  and native routes, and let "no metric available" be an honest answer. **Not decided.**

## 4.1 Cross-Session Constraints ( agreed 2026-08-28 )

A peer session holds Christian's go for the DR-127 certified-exit restructure, which rewrites
`cl_FEM_Controller.cpp/.hpp`, `cl_FEM_DofMgr_SolverData.*` and `cl_FEM_DofManager.*`. Deconflicted
and now normative in `tmp/ai_exchange/dr127_certified_exit_plan_v4.md`:

- the four `dump_system_if_requested` call sites travel with the **assembly** phase, not the solve —
  the dump's contract is the assembled system BEFORE the solve;
- `finalize()`'s conditioning block is left alone by that session, so **R5/R8's move to
  `finalize_run()` lands after theirs**;
- `mMagneticHitTarget` in the Controller diff belongs to neither session — a third writer is in that
  file, and its ownership is being routed to Christian.

DR-127 is on hold on Christian's word; that session will ping before touching Controller.

## 5. Definition of Done

- [ ] `compute conditioning : true` on the tapestack3d thermal system reports a number, or says
      plainly which metric is unavailable and why.
- [ ] The solver-native metric is armed for exactly ONE solve per run, and the eigen path runs
      exactly ONCE per run. No per-timestep conditioning cost remains.
- [ ] Production solve paths show an unchanged runtime JOB sequence and unchanged numerics while
      reuse is disabled ( the honest form of the v1 "bit-identical" claim — layout/vtable changes
      are unavoidable and irrelevant ).
- [ ] The reported quantity is NAMED wherever it is printed.
- [ ] The shift-invert result is validated against the measured ground truth:
      `lambda_min = 6.175959e-08`, `lambda_max = 1.830271`, `kappa_2 = 2.9635e7`. **This plan has a
      known-good answer to check against, which is rare — use it.**
- [x] `make check` clean — **RUN 2026-08-29 by Christian, green.** This is the WIRING gate, and it
      is the one a scratchpad probe structurally cannot give: the probe hand-links, so a missing
      `LIBLIST` entry or an unbuilt Fortran shim would not surface there. Verified independently
      that the archive was rebuilt AFTER the last edit and contains `arpack_si_step`,
      `arpack_si_extract`, `arpack_symmetric_eigen`, `freeze_factorization` and `select_job`.
      It also covers the regression risk that mattered most: `select_job()` and the `free()`
      changes sit on the PRODUCTION solve path, and the suite exercises solvers with reuse
      disabled.
      **What it does NOT cover:** no test reaches `run_shift_invert()` — it is only entered through
      `compute_conditioning()`, which needs a deck — so the numeric behaviour rests on §3.5's
      analytic gate, not on the suite.
