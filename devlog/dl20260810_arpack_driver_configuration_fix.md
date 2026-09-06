# ARPACK Driver Configuration Fix — Work Arrays, `ncv`, and Error Decoding

**Date:** 2026-08-10
**Purpose:** Fix the ARPACK-ng eigenvalue driver after "ARPACK is slow" was raised as a reason to
replace it with SLEPc; establish what the slowness actually was
**Files:** `src/sparse/arpacktools.{f90,hpp,cpp}`, `src/sparse/CMakeLists.txt`,
`src/fem/kernel/cl_FEM_DofMgr_EigenValues.{hpp,cpp}`

## How this started

The session opened as a SLEPc question: ARPACK is slow, SLEPc is already in SCLS, should the
conditioning diagnostic move over. The assessment that followed found the slowness was **not
ARPACK's** — the driver was misconfigured in ways that guarantee bad convergence regardless of
which Arnoldi library sits underneath. Christian's ruling: fix the configuration first, and keep
the driver in Fortran (ARPACK is an F77 library; SCLS owns the int-width / MPI / BLAS consistency,
so the 32-vs-64-bit concern raised during the assessment is handled by construction).

The SLEPc port was **not** executed and remains open. Its assessment is recorded below.

## What was actually wrong (confidence: high — all verified against arpack-ng 3.9.0 `SRC/`)

Two of these are live memory-corruption bugs, not merely bad tuning.

1. **`slct` out-of-bounds write, unconditional.** `dneupd.f:330` declares `logical select(ncv)`
   and `dneupd.f:548-551` unconditionally runs `do 10 j = 1,ncv ... select(j) = .false.`. The
   driver declared `slct(nev)` — a **one-element stack array** at `nev = 1` — while passing
   `ncv = 3`. Three logicals written past the end, on every call, for as long as this file has
   existed.
2. **`z` out-of-bounds write, conditional.** dneupd's contract is `N by NEV+1`
   (`dneupd.f:91`); the driver declared `z(n, nev)`. ARPACK declares it `z(ldz,*)` internally, so
   no compiler can catch it, and `dneupd.f:738` copies `nconv` columns. `nconv` reaches `nev+1`
   exactly when a complex-conjugate pair straddles the requested cut — which for a *non-symmetric*
   (`dnaupd`/`dneupd`) driver is the expected case, not an exotic one.
3. **`ncv = 2*nev+1 = 3`** (`arpacktools.f90:90`, old). Remark 4 of `dnaupd.f` gives `2*nev+1` as
   a *floor*; at `nev = 1` that floor is a Krylov subspace of dimension 3, which restarts
   constantly and discards the subspace each cycle.
4. **`ncv` could not be raised in place.** Every work array was sized to an expression that
   equals the required size *only* when `ncv == 2*nev+1` — `workev` declared `6*nev+3` (needs
   `3*ncv`), `workl` declared `(2*nev+1)*(6*nev+11)` against a computed
   `lworkl = ncv*(3*ncv+8)`, `vectors(n, 2*nev+1)` against `ncv` columns. The whole block was
   pinned to one value of `ncv`.
5. **The arrays were on the stack.** `residuals(n)`, `workd(3n)`, `vectors(n,·)`, `z(n,·)` were
   automatic arrays, and `config_gcc.cmake:38` puts `-fopenmp` on the Fortran flags, which keeps
   them there. At `n = 1e5` with a useful `ncv = 20` that is ~18 MB of Arnoldi vectors alone.
   Surviving at `ncv = 3` was the only reason this ever ran.
6. **`mNumMaxIter = mN * mN`** (`cl_FEM_DofMgr_EigenValues.cpp:79`, old) clobbered the sane member
   default of 300. `iparam(3)` counts *implicit restarts*, not matvecs — for `n = 1e5` this asked
   for 1e10 of them, which is why a non-converging run ground instead of failing.
7. **The error decoder used the wrong table.** `arpack::check` decoded results from *both*
   routines against **dneupd's** table. `dnaupd` and `dneupd` do not share one. Consequences:
   `dnaupd info = 3` — *"No shifts could be applied… One possibility is to increase the size of
   NCV relative to NEV"*, i.e. the exact diagnosis of item 3 — had no case at all and fell through
   to `"Unknown Error"`; `dnaupd info = 1` (max iterations, a benign warning) was reported as
   `"Shur reordering failed"`. **ARPACK had been reporting the bug and the decoder was
   swallowing it.**
8. **`NCONV` was never read**, so a partially-converged result was indistinguishable from a
   converged one.

## Changes

`src/sparse/arpacktools.f90`
- Work arrays are `allocatable` and sized from the actual `ncv` per the documented contracts:
  `slct(ncv)`, `z(n,nev+1)`, `workev(3*ncv)`, `vectors(n,ncv)`, `workl(lworkl)`. This fixes items
  1, 2 and 4 and moves the large buffers off the stack (item 5).
- `ncv = min(n, max(2*nev+1, 20))`.
- `ido` is its own variable instead of sharing `info(2)`; the two routines now get separate status
  slots, because they cannot share a decoder.
- `info` extended to 6 entries: dnaupd flag, dneupd flag, `NCONV`, `MXITER`, `NUMOP`, `NUMREO`.
- `ido = 2` services `y = x` (B = I) rather than `cycle`-ing without progress.

`src/sparse/arpacktools.{hpp,cpp}`
- `check` split into `check_naupd` / `check_neupd`, each carrying its routine's real table
  transcribed from the 3.9.0 source. dnaupd `info = 1` and `= 3` are now warnings that let the run
  continue instead of fatal mislabels.
- Named `gInfo*` indices for the status array; added the `typedefs.hpp` include the header always
  needed but got away without (the `.cpp` happened to include it first).

`src/fem/kernel/cl_FEM_DofMgr_EigenValues.{hpp,cpp}`
- Deleted the `mNumMaxIter = mN * mN` clobber (item 6).
- `NCONV` is checked before the result is used, with a debug bound against the λ buffer length.
- **`nev` now follows the setters.** `set_num_minvals` / `set_num_maxvals` previously resized the
  λ buffers while `nev` stayed a hardcoded `const mNumEigenValues = 1`, so asking for five
  eigenvalues computed one. `run()` now takes its count from whichever setter matches the job, and
  both setters reject a non-positive count.
- The returned extremum is picked by scanning the converged set rather than reading slot 0 —
  dneupd promises no ordering the caller can rely on, and with `nev > 1` slot 0 is not necessarily
  the extreme value.
- New `compute_smallest_eigenvalues()` / `compute_largest_eigenvalues()` plus `lambda_real()`,
  `lambda_imag()`, `number_of_converged_values()`, so the computed spectrum is reachable — this is
  the "actual eigenvalues as an asset" Christian asked for. Master-rank only; this solver does not
  run distributed.
- `run()` is guarded by `BELFEM_ARPACK` with a clean error, and `arpacktools.cpp` moved inside the
  `USE_ARPACK` block in `src/sparse/CMakeLists.txt` — it was compiled unconditionally while the
  `.f90` was gated, so `USE_ARPACK=OFF` produced a link error instead of a clean degrade.

## Not fixed — the open decisions

- **`which = 'SM'` is still `'SM'`** (confidence: high). None of the above makes mode-1 Arnoldi
  resolve the small end of the spectrum; that is a structural property of the method, not a
  tuning parameter. The remedy is shift-and-invert (mode 3, `sigma = 0`, `'LM'` on `A⁻¹`), which
  needs a *solve* in the reverse-communication slot. `Solver::solve` does reuse an initialized
  wrapper (`cl_Solver.cpp:151`) and nothing in the FEM layer calls `free()`, but the MUMPS wrapper
  drives JOB 5/6 (factorize+solve) rather than JOB 3 (solve-only), so each shift-invert matvec
  would trigger a full refactorization. Shift-invert therefore depends on exposing a solve-only
  path in the wrappers.
- **Behaviour change to watch:** with the iteration limit now sane, `compute_conditioning()`'s
  `run(0)` may report a clean non-convergence where it previously ground silently. That is the
  intended improvement, but it is a visible change on that path.
- **SLEPc remains open.** SLEPc 3.25.1 and PETSc 3.25.2 are in SCLS, PETSc has MUMPS /
  SuperLU_dist / STRUMPACK, and BELFEM already builds a PETSc `Mat` from its CSR arrays
  (`cl_SolverPETSC.cpp:416`) plus `DistMatrixAIJ<PetscInt>`, so the port is incremental. Note for
  whenever it is taken up: SLEPc would only be *faster* on the small end because `ST`
  shift-invert is a two-line setup, not because its Arnoldi is better — `EPS` with
  `EPS_SMALLEST_MAGNITUDE` and no spectral transform would be exactly as slow. If it happens, the
  target is the **`SVD`** module, not `EPS` (see below). SLEPc in SCLS is itself built
  `SLEPC_HAVE_ARPACK 1`.

## Separate finding: the diagnostic is not a condition number

Confidence: high. `compute_conditioning()` returns `|λ_max| / |λ_min|`, which equals κ₂ only for
normal matrices. The h-φ Jacobian is not symmetric, let alone normal, so the reported figure can
be off by orders of magnitude. The MUMPS branch (`cl_FEM_Controller.cpp:2189`) meanwhile returns
`get_cond0()`, an Arioli/Demmel/Duff error-analysis estimate for the linear *system* — so the two
branches of the same function do not report the same quantity today. A true κ₂ needs singular
values (SLEPc `SVD`, or a Hager–Higham 1-norm estimator against the live factorization, which is
cheap and sufficient given Christian's ruling that the figure is diagnostic-only).

This is worth flagging beyond cosmetics: `todo/timestep_collapse_residual_floor_plan.md:101-106`
carries a live hypothesis that a residual floor is "plausibly conditioning-dominated". That
hypothesis is currently being weighed against a number which is not a condition number, and which
— on the non-MUMPS path — was being produced by the driver fixed above.

## Pending gate

Fortran syntax-checked standalone (`gfortran -fsyntax-only -cpp -fopenmp`, clean); no stale callers
of the removed `arpack::check`. **Not compiled in-tree and not run** — build and smoke gate are
Christian's. The clangd diagnostics on `cl_FEM_DofMgr_EigenValues.cpp` are stale-index noise
(`'armadillo' file not found` — the TU does not parse) and are unchanged from before the edits.
