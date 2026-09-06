# SLEPc Eigensolver Integration

**Date:** 2026-08-10
**Purpose:** Interface SLEPc as a second eigenvalue backend alongside the (now repaired) ARPACK-ng
driver, so BELFEM can resolve the *small* end of the spectrum — which mode-1 Arnoldi structurally
cannot — via a spectral transform backed by a sparse direct solve. Mechanism: wrap the existing
`SpMatrix` CSR arrays in a PETSc `Mat` and drive SLEPc `EPS` with `STSINVERT`.
**Module:** `src/sparse` (+ `src/fem/kernel`, `src/comm`)
**AIs involved:** Claude (exploration + plan)
**Status:** **DEFERRED 2026-08-11** (currentness sweep) — moved to `deferred/`, and the move is
bookkeeping on a decision already taken, not a new one: the deprioritization was measured and
recorded on 2026-08-10, the production line for the conditioning number is
`conditioning_diagnostic_backends.md` (MUMPS accessors landed, SuperLU next), and nothing is
waiting on this file. **The revival condition is stated and cheap:** re-measure the ARPACK path
at CORC scale; if it stops being competitive there, this plan is ready to execute as written.
Nothing here was refuted — R1–R7, the gap table, O1–O3 and Appendix A's four scaffold defects all
stand. *(Superseded status, kept for the record:)* PLAN — drafted, not audited, **deprioritized by measurement on 2026-08-10** (see §1.1):
at the size measured, the repaired ARPACK path is projected to cost ~70–150 ms in a release build,
which a shift-invert factorization would not beat. The case now rests on **large** problems only
and needs a re-measurement at CORC scale before it is worth executing. No source modified.

> **Scope guards:**
> - The repaired ARPACK path stays. This is an *addition*, not a replacement — see O2.
> - Parallel/distributed eigensolves are OUT of scope. `EigenValues::run` is master-rank-only
>   (`cl_FEM_DofMgr_EigenValues.cpp:146`) and stays that way; all PETSc objects are built on
>   `PETSC_COMM_SELF`.
> - Eigen*vectors* are out of scope; only values are extracted.
> - No deck-visible behaviour change without an explicit input key (see R7 and the
>   input-reference sync rule).

---

## 1. Current Behaviour and How It Fails

The repaired ARPACK driver (see `devlog/dl20260810_arpack_driver_configuration_fix.md`) now sizes
its work arrays correctly, uses a workable `ncv`, and decodes both routines' error tables. What it
still cannot do is the small end.

| Failure | Mechanism | Evidence |
|---|---|---|
| `which='SM'` does not converge | Mode-1 Arnoldi builds a Krylov space from `A`; it has essentially no resolution near the smallest-magnitude eigenvalues. Not a tuning parameter — a property of the method. | `arpacktools.f90` sets `bmat='I'`, `iparam(7)=1` |
| Shift-and-invert is blocked | Needs a *solve* in the reverse-communication slot. `Solver::solve` reuses an initialized wrapper (`cl_Solver.cpp:151`) and nothing in the FEM layer calls `free()`, but the MUMPS wrapper drives JOB 5/6 (factorize+solve), not JOB 3. Every Arnoldi matvec would refactorize. | `cl_SolverMUMPS.cpp:340-350` |
| Generalized problems impossible | Driver is hardwired to the standard problem. `Ax = λBx` needs a new driver regardless of backend. | `arpacktools.f90` (`bmat='I'`) |
| The "conditioning" figure is not κ₂ | Returns `\|λ_max\|/\|λ_min\|`, which equals κ₂ only for normal matrices; the h-φ Jacobian is not symmetric. The MUMPS branch meanwhile returns an Arioli/Demmel/Duff *system* estimate, so the two branches of one function report different quantities. | `cl_FEM_DofMgr_EigenValues.cpp` `compute_conditioning`; `cl_FEM_Controller.cpp:2185-2197` |

**Bottom line:** the remaining gap is not speed — it is that the small end of the spectrum requires
a factorization the ARPACK route cannot reuse, and SLEPc brings its own solver stack that can.

## 1.1 Measurement (2026-08-10) — why this is deprioritized

Probe data from the repaired driver, n = 10428, nnz = 70864, nev = 1, ncv = 20, 2 OpenMP threads:

| end | nconv | restarts | numop | time |
|---|---|---|---|---|
| `SM` | 1 | 122 | 1231 | 769 ms |
| `LM` | 1 | 5 | 61 | 35 ms |

Three things fall out, in order of importance:

1. **`SM` converges.** `nconv 1`, `info 0/0` — it is expensive, not broken. The prediction that it
   would hit the restart cap and return nothing was **wrong**. It costs 20× the matvecs of `LM`
   (1231 vs 61), which is the structural cost of the small end, but it does produce a value.
2. **The measurement was taken in a debug build.** `cmake-build-debug/CMakeCache.txt` has
   `USE_DEBUG:BOOL=ON`, and the Fortran flags are `-O0 -fcheck=bounds -fbacktrace`, so every
   `values(j)` / `indices(j)` / `x(indices(j))` in the matvec pays a bounds test. A standalone
   benchmark of the identical kernel measured 0.4879 ms/matvec at those flags versus 0.0865 ms at
   `-O2` — **5.6×**. (Confidence: high on the mechanism; the projection below is an extrapolation,
   **not** a release-build measurement.)
3. **Threading is fine.** 1 thread gives 1424/64 ms against 769/35 ms at 2 threads — 1.85×, near
   linear. An earlier hypothesis that OpenMP overhead dominated at this matrix size is **refuted**.
   Amdahl on that ratio also puts ~92% of the time inside the matvec loop, i.e. ARPACK's internal
   orthogonalization is a minor term.

**Consequence:** projected release cost is roughly 150 ms as-is, or ~70 ms with the scalar-accumulator
fix to the matvec (a further 2.2× at `-O2`, measured in the same benchmark — `y(i)` cannot be held in
a register because `x` and `y` are both pointers into `workd` and may alias). A sparse LU of this
matrix plus a handful of shift-invert solves would land in the same range or worse, so **SLEPc buys
nothing at this size**.

What survives: `SM`'s matvec count grows as the mesh refines and conditioning worsens, while a
factorization scales more gracefully. So the crossover is a real question — it just has to be
demonstrated at CORC scale, with the same probe, in a **release** build, before this plan is worth
executing.

## 2. Architecture: Why SLEPc Is the Right Spine

The cost comparison is **not** "ARPACK vs SLEPc". It is:

- **ARPACK route:** add a solve-only (JOB=3) path across `cl_SolverMUMPS` / `STRUMPACK` /
  `UMFPACK` / `PARDISO` — i.e. surgery on the wrappers every production solve depends on,
  including code under active debugging in the timestep-collapse campaign.
- **SLEPc route:** SLEPc never touches BELFEM's `Solver`. `ST` shift-invert is backed by PETSc's
  own `KSP`/`PC`, and `PCLU` factors **once** and reuses the factorization across all inner solves.

The thing that is expensive and risky in the first route is free in the second. (Confidence: high —
the wrapper JOB behaviour is cited above; PETSc factorization reuse is standard `PCLU` semantics.)

Integration is incremental because PETSc is already wired: `MatCreateSeqAIJWithArrays` is already
used at `cl_SolverPETSC.cpp:416`, `DistMatrixAIJ<PetscInt>` already exists (`petsctools.hpp`), and
`config_petsc.cmake` is a working template for the SCLS prefix-install search.

**Rejected alternative:** porting only the conditioning diagnostic to SLEPc. Christian's ruling is
that conditioning is diagnostic-only and an estimate suffices, so that alone would not justify a
new TPL. The justification is the *capability* (small end, and later generalized problems).

### Verified environment facts (2026-08-10, re-verified against the correct SCLS flavor)

> **Flavor trap — read this first.** Three SCLS trees are installed: `/opt/scls/gcc`,
> `/opt/scls/mkl`, `/opt/scls/debug` (plus `/opt/scls/cea`). BELFEM links **`/opt/scls/mkl`** —
> confirmed by the build's own include flags (`-I/opt/scls/mkl/include`, `flags.make`). The trees
> differ in real ways: the `gcc` flavor's `libarpack.so` links OpenBLAS and has no FEAST, the `mkl`
> flavor links MKL and does. Any capability check run against `/opt/scls/gcc` is about a library
> this build does not load.

| Fact | Value | Source (mkl flavor) |
|---|---|---|
| SLEPc / PETSc | 3.25.1 / 3.25.2 | `/opt/scls/mkl/lib/libslepc.so.3.25.1`, `petscversion.h:8-10` |
| PETSc integer width | **32-bit** (no `PETSC_USE_64BIT_INDICES`) | `petscconf.h` |
| PETSc scalar | real double | `petscconf.h:213` |
| Direct solvers available to `PCLU` | MUMPS, STRUMPACK, SuperLU, SuperLU_dist | `petscconf.h:106,144-147` |
| SLEPc backends | ARPACK, **FEAST**, ScaLAPACK | `slepcconf.h:8,11` |
| ARPACK integer width | 32-bit (`INTERFACE64 0`) | `arpack/arpackdef.h:6` |
| ARPACK's BLAS | MKL LP64 + GNU OpenMP threading | `ldd libarpack.so` → `libmkl_gf_lp64`, `libmkl_gnu_thread`, `libmkl_core` |
| PETSc lifecycle owner | `Communicator` | `cl_Communicator.cpp:114` init, `:279` finalize |
| `PetscInt == int` already asserted at runtime | yes | `cl_Communicator.cpp` (post-`PetscInitialize` type checks) |

The int-width chain is consistent by construction and needs no BELFEM-side guard: MKL LP64 ↔
ARPACK `INTERFACE64 0` ↔ PETSc 32-bit ↔ `USE_MKL_64BIT_API=OFF`. SCLS owns that consistency
(Christian, 2026-08-10); a 64-bit build would rebuild the whole tree as a unit.

`SLEPC_HAVE_FEAST` is worth knowing about but is not the target here: FEAST computes eigenvalues
inside a contour and still needs linear solves, so it carries the same factorization cost as
shift-invert without being a better fit for "the smallest one".

## 3. Gap Table

| # | State | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| 1 | SLEPc discovery + link rank | build | no | (c) | rank **10**, one above PETSc's 9; line assembled 10→0 (`belfem_find_package.cmake:3`) |
| 2 | `SlepcInitialize` / `SlepcFinalize` | EPS type registration | no | (c) | must bracket inside PETSc's; `SlepcInitialize` no-ops the PETSc part if already initialized |
| 3 | `SpMatrix` → `Mat` | operator | partially | (a) | reuse the `cl_SolverPETSC.cpp:416` pattern |
| 4 | Indexing base | `MatCreateSeqAIJWithArrays` needs **0-based** | conflicting | (c) | PETSc path sets `Cpp` (`cl_SolverPETSC.cpp:411`), ARPACK path sets `Fortran` (`cl_FEM_DofMgr_EigenValues.cpp:150`) — whichever backend runs must set its own |
| 5 | CSC matrices | operator | (a) | (a) | CSC fed as CSR yields `Aᵀ`; `σ(A) = σ(Aᵀ)` and eigenvalues match, so this is safe here (unlike the solve path, `cl_SolverPETSC.cpp:263-265`) |
| 6 | Array lifetime | `Mat` validity | no | (c) | `MatCreateSeqAIJWithArrays` **borrows**, does not copy — arrays must outlive the `Mat` and not be reallocated |
| 7 | Spectral transform | the whole point | no | (c) | `STSINVERT` + target 0 + `EPS_TARGET_MAGNITUDE`; `KSPPREONLY` + `PCLU` + MUMPS |
| 8 | Problem type | correctness | no | (c) | `EPS_NHEP` — Jacobian is **not** Hermitian |
| 9 | Convergence reporting | trustworthy result | no | (c) | `EPSGetConverged` before reading values, mirroring the `NCONV` check now in the ARPACK path |
| 10 | κ₂ definition | the diagnostic | no | (b) | see **O1** |

### 3.1 Cross-cutting findings

- **Row 4 is the most likely first bug.** The two backends demand opposite indexing bases on the
  *same* `mK`. Whichever runs must set its own base immediately before use; assuming the incoming
  state is a latent, intermittent defect.
- **Row 6 is the second.** Unlike most PETSc constructors, `MatCreateSeqAIJWithArrays` wraps caller
  memory. If `mK` is rebuilt or resized while a `Mat` is alive, the `Mat` dangles.
- **Row 7 is what makes this worth doing at all.** An `EPS` without `ST` is exactly as slow as the
  ARPACK path was — the port would appear to work and change nothing.

## 4. Ordered Steps

- [ ] **R1** — `config/linalg/config_slepc.cmake` mirroring `config_petsc.cmake` (prefix-install
      branch via `belfem_find_package`, plus the `SLEPC_DIR`/`PETSC_ARCH` in-tree branch);
      `option(USE_SLEPC ...)` defaulting to `USE_PETSC`; hard error if `USE_PETSC` is off;
      `BELFEM_SLEPC` define; `belfem_link_libraries( 10 ... )`; include it after
      `config_petsc.cmake` in `CMakeLists.txt`.
- [ ] **R2** *(after R1)* — `SlepcInitialize` / `SlepcFinalize` in `cl_Communicator.cpp`, nested
      strictly inside the PETSc pair (`:114` / `:279`), under `#ifdef BELFEM_SLEPC`.
- [ ] **R3** *(after R2)* — `src/sparse/slepctools.{hpp,cpp}`: one entry point taking an
      `SpMatrix &`, an end-of-spectrum flag, `nev`, tolerance and iteration cap, filling
      `Vector<real>` real/imag parts and returning `nconv`. Non-SLEPc builds get the stub-typedef
      treatment `petsctools.hpp` already uses. Must set the indexing base itself (row 4), keep the
      borrowed arrays alive across the solve (row 6), and wire `ST` for the small end (row 7).
- [ ] **R4** *(after R3)* — dispatch in `EigenValues::run()`. Keep the ARPACK path intact; see O2
      for which backend is preferred.
- [ ] **R5** *(after R4)* — resolve O1 and make `compute_conditioning()` report a defensible
      quantity, or document precisely what it reports.
- [ ] **R6** *(after R5)* — build gate: configure and link with `USE_SLEPC` both ON and OFF, and
      with `USE_ARPACK=OFF` (the gate repaired 2026-08-10) to confirm the backends are genuinely
      independent.
- [ ] **R7** *(after R6)* — smoke gate on a real deck with `compute conditioning : true`, ARPACK
      vs SLEPc on the same matrix, comparing λ_max and wall time. **If a deck key is added**,
      `doc/input_file_reference.md` is updated in the same turn per the standing input-reference
      sync rule, with a Codex prose pass on the touched section.

## 5. Open Design Questions (not silently decided)

- **O1 — `EPS` ratio or `SVD` for the conditioning number?** `EPS` gives `|λ_max|/|λ_min|`, which
  is not κ₂ for a non-normal matrix; `SVD` (`SVDCROSS` / `SVDTRLANCZOS`) gives the real κ₂ =
  σ_max/σ_min at higher cost. Christian's ruling is that the figure is diagnostic-only and an
  estimate suffices, which argues for the cheap path — but "cheap *and* mislabelled" is the
  current state, and it is currently being weighed against the residual-floor hypothesis in
  `todo/timestep_collapse_residual_floor_plan.md:101-106`. A third option is a Hager–Higham 1-norm
  estimate, which is what MUMPS's `get_cond0()` already reports, giving one consistent quantity
  across all solvers. **Not decided.**
- **O2 — Does SLEPc replace ARPACK, or coexist?** Coexisting keeps a fallback and makes R7's A/B
  possible, at the cost of two code paths. Replacing means one path but drops a working backend
  and makes `USE_SLEPC=OFF` builds lose eigenvalues entirely. Sub-question: if both are present,
  is the preference per-end (SLEPc for small, either for large) or global? **Not decided.**
- **O3 — Shift target for `STSINVERT`.** σ = 0 gives smallest-magnitude directly but requires `A`
  to be factorizable; a singular or near-singular Jacobian would fail in `PCLU`. A small nonzero σ
  avoids that but biases which eigenvalues are found. Needs a decision on the failure path —
  clean diagnostic, or automatic retry at a perturbed σ. **Not decided.**

## 6. Definition-of-Done Checklist

- [ ] Every gap-table row mapped to a step or an open question.
- [ ] `USE_SLEPC` ON and OFF both configure, compile and link.
- [ ] `USE_ARPACK=OFF` with `USE_SLEPC=ON` still produces eigenvalues.
- [ ] SLEPc and ARPACK agree on λ_max for the same matrix (R7).
- [ ] The small end returns a converged value in bounded time — the capability this plan exists for.
- [ ] O1/O2/O3 resolved in place with dates, not silently.
- [ ] `doc/input_file_reference.md` synced if a deck key was added.

---

## Appendix A — Corrections to the circulating scaffold

A SLEPc scaffold (Gemini, 2026-08-10) was reviewed. It is directionally correct but carries four
defects worth recording, because they are the natural things to get wrong here and will otherwise
be re-proposed:

| # | Scaffold | Correct | Why it matters |
|---|---|---|---|
| 1 | `EPSSetWhichEigenpairs(EPS_SMALLEST_MAGNITUDE)` with no `ST` | `STSINVERT` + `EPSSetTarget(0)` + `EPS_TARGET_MAGNITUDE`, with `KSPPREONLY` / `PCLU` / MUMPS | **The whole point.** Without the transform this is exactly the trap the ARPACK driver was in; the port would "work" and be no faster |
| 2 | `EPSSetProblemType(EPS_GHEP)` | `EPS_NHEP` | GHEP asserts Hermitian. The h-φ Jacobian is not symmetric — wrong answers, not just slow ones |
| 3 | Builds an identity `Mat B` for a "generalized" problem | `EPSSetOperators(eps, A, NULL)` | B = I *is* the standard problem; the identity matrix is pure overhead |
| 4 | `Mat` on `PETSC_COMM_SELF`, `EPS` on `PETSC_COMM_WORLD` | one communicator throughout — `PETSC_COMM_SELF`, matching the master-rank-only guard | Mismatched communicators are a parallel-only failure, i.e. invisible in serial testing |

Minor: `SlepInitialize`/`SlepFinalize` are missing their `c`; `EPS_PGN`/`EPS_GN` in its comments are
not SLEPc problem types; and it has no error checking (BELFEM has `petsctools_error_message` for
this).
