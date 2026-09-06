# SuperLU Solver Wrapper — Audit and Fixes

**Date:** 2026-06-08
**Purpose:** Record the audit + remediation of the new sequential SuperLU solver wrapper (`cl_SolverSUPERLU.{hpp,cpp}`), unifying the Claude and Codex findings, and document two pre-existing header-integration blockers that prevented the wrapper from ever compiling.
**Module:** `src/sparse`

## Context

New BSD-licensed sequential SuperLU wrapper intended to replace GPL UMFPACK
(single process, no MPI/threading). The matrix is CSC, zero-based; SpMatrix
arrays are wrapped by pointer (no copy); the solve is split into
`symbolic()` / `numeric()` / `solve()` for reuse across Newton iterations
(values change, pattern fixed). Files are new/untracked.

A two-AI read-only audit (Claude broad, Codex precision) was reconciled, then
fixes were applied one step at a time with the user.

## Findings (unified, by severity)

**Tier 1 — correctness/safety**
- `mAC` (SLU_NCP) never allocated before `sp_preorder` writes through it and
  allocates `AC->Store` → UB/crash. (`sp_preorder.c:92`)
- `mGlu` null passed to `dgstrf` (deref in `dLUMemInit`) → crash.
- `mU` (SLU_NC) destroyed with `Destroy_SuperNode_Matrix` (SCformat layout) →
  invalid free / heap corruption. Correct call is `Destroy_CompCol_Matrix`.
- Multi-RHS `solve(Matrix&,Matrix&)`: used `mMatrix->n_cols()` as nrhs, wrapped
  `aLHS.data()` (uninitialised) instead of the packed RHS, and copied the
  untouched RHS back into `aRHS` → wrong results + OOB for nrhs != n.

**Reconciled to NOT a bug (Codex retracted):** the CSR branch is valid —
`CSR(A)` shares arrays with `CSC(Aᵀ)`, and `Trans=TRANS` makes `dgstrs` solve
`(Aᵀ)ᵀx = Ax = b`. Kept; comment corrected (`mAC = Aᵀ·Pc` in CSR mode). A
dedicated CSR test is still TODO.

**Tier 2 — robustness/portability**
- No index-ABI guard: `belfem::int_t` (32/64-bit) handed to SuperLU `::int_t`
  (here `int`). Added `static_assert(sizeof(belfem::int_t)==sizeof(::int_t))`
  plus runtime `nrow/ncol` fit-`int` and `nnz` fit-`::int_t` checks.
- `if(mOptions==nullptr) delete mOptions` was inverted (leak on re-init);
  options were also never freed → fixed (delete in destructor) and `Fact`
  reset to `DOFACT` at the end of `symbolic()`.
- `dgstrf` info: `info<0` (illegal arg) was mislabelled singular; singular
  pivot index is 1-based. Split into illegal-arg / out-of-memory / singular.
- SamePattern value alias (`ACstore->nzval = Astore->nzval`) documented as
  requiring no SpMatrix reassignment/resize between `symbolic()` and `free()`.
- `malloc` results now checked; empty/uninitialised matrix rejected in
  `symbolic()`.

**Tier 3 — clarity**
- Removed dead `error_message()` declaration; corrected `numeric()` error text
  ("called before symbolic()"); corrected the `mGlu` "reused for SamePattern"
  comment (plain SamePattern reuses perm_c/etree, not the LU workspace).

**License:** clean — only `dgstrf`/`dgstrs`/`get_perm_c`/`sp_preorder`/
`dCreate_*`/`Destroy_*` are reachable; no MC64/`dldperm`/ILU path.

## Pre-existing build blockers discovered (the wrapper had never compiled)

No `cl_SolverSUPERLU.cpp.o` existed in any build tree. Two distinct
header-integration failures, surfaced only when building the real target:

1. **Armadillo backend — `superlu_enum_consts.h` guard collision.** Armadillo
   includes that header inside `namespace arma::superlu`, defining the global
   guard `__SUPERLU_ENUM_CONSTS`; the later global `<slu_ddefs.h>` then skips it
   and `fact_t`/`yes_no_t`/… are undefined. Fix: `config/linalg/config_superlu.cmake`
   adds `ARMA_DONT_USE_SUPERLU` under `if(USE_MATRIX_ARMADILLO)` (BELFEM never
   uses `arma::spsolve`). NB: guard var is `USE_MATRIX_ARMADILLO`, not
   `USE_ARMADILLO`.

2. **Blaze backend (the ACTIVE build) — global BLAS prototype clash.**
   `<slu_ddefs.h>` unconditionally declares the global `extern "C"` BLAS
   symbols `dcopy_ daxpy_ dgemm_ dgemv_ dtrsm_ dtrsv_`; Blaze
   (`blaze/math/lapack/clapack/*.h`) declares the same global C symbols with
   incompatible prototypes (`blaze::blas_int_t*` + Fortran char-length args).
   They are the same C-linkage entity, so two prototypes are ill-formed and
   `::` cannot disambiguate. Fix: in `cl_SolverSUPERLU.hpp`, rename SuperLU's
   six (unused, never-called) BLAS prototypes to `belfem_slu_*` for the
   duration of the `<slu_ddefs.h>` include only (`#define` … `#include` …
   `#undef`), leaving the backend's declarations and `libsuperlu.a`'s real
   symbols untouched.

The stale `compile_commands.json` (Armadillo-flavoured) initially masked #2;
the active configuration is Blaze + NETLIB (`flags.make`).

## Verification

`make libbelfem_sparse.a` builds clean under the real Blaze+NETLIB toolchain
(exit 0, no conflicting-declaration errors); `cl_Solver.cpp` — the other TU
pulling both the wrapper header and Blaze — also compiles, confirming the
header-confined macro is include-order independent.

## Remaining / TODO

- Full executable link + a runtime solve to confirm `dgstrf`/`dgstrs` resolve
  and produce correct results (a `make reset && make <target>` is needed for
  downstream consumers to pick up the rebuilt static lib).
- CSR-path unit test (the `CSC(Aᵀ)`+TRANS trick).
- The `ARMA_DONT_USE_SUPERLU` path is correct by construction but not yet
  exercised (current build is Blaze).

## Round 2 — Codex integration-wiring audit (same day)

A follow-up Codex pass found the wrapper was never registered in the solver
enum/dispatch infrastructure. Verified and fixed:

- **Default solver** (`en_SolverEnums.hpp`): added `#elif BELFEM_SUPERLU →
  SolverType::SUPERLU` so a SuperLU-only build no longer defaults to UNDEFINED.
- **`to_string(SolverType)`** (`en_SolverEnums.cpp`): added the `SUPERLU →
  "SUPERLU"` case. This was a latent bug: with SUPERLU hitting the default
  `"UNKNOWN"`, `solver_type("unknown")` actually returned `SUPERLU` (loop match
  at enum slot 1) while `"superlu"` could not be parsed at all.
- **`matrix_type()`** (`fn_matrix_type.hpp`) and **`preferred_matrix_format()`**
  (`fn_preferred_matrix_format.hpp`): added `SUPERLU → CSC` (native SLU_NC, no
  transpose path); previously threw / silently fell through to CSR.
- **RHS preconditions** (`cl_SolverSUPERLU.cpp`): both solve overloads now use
  `BELFEM_ERROR` (was `BELFEM_ASSERT`, compiled out in release — a short RHS
  could read/write past storage).
- **`free()` vs base reset:** Codex's "just call `Wrapper::free()`" would
  **leak** here, because `symbolic()` calls `free()` on matrix-pointer change;
  resetting `mIsInitialized` makes the driver re-`initialize()` next iteration,
  which clears the `mHave*` flags without freeing → orphaned factors. Fixed by
  splitting: `free_factorization()` (SuperLU resources only, called by
  `symbolic()` re-entry) and the public `free()` = `free_factorization()` +
  `Wrapper::free()` (driver/destructor only). Preserves symbolic/numeric reuse.
- **`tData` leak:** matrix solve now frees the packed buffer before the (fatal)
  `dgstrs` error check and scatters only on success.

**Round 3 (Codex residual):** in debug, `BELFEM_ERROR` *throws*
(`assert.hpp:98`) and BELFEM tests catch it, so the `dgstrf`-failure leak of
`mL/mU/mGlu` is test-visible after all. Fixed by setting `mHaveNumeric = true`
right after the `info < 0` (illegal-argument) check and before the OOM/singular
checks: for `info >= 0` the factors exist (valid when singular, partial on OOM)
so `free_factorization()` can release them on a thrown error; `info < 0` returns
before any factor storage is built, so it is deliberately left unowned (the
empty structs are an unreachable wrapper-bug path, and destroying their
uninitialised `Store` would crash).

Open policy note (Codex): `gDefaultSolver` prefers SuiteSparse over SuperLU when
both are enabled (`BELFEM_SUITESPARSE` precedes `BELFEM_SUPERLU`). Left as-is
(explicit SuiteSparse opt-in wins); flip the order if BSD SuperLU should be the
default replacement for UMFPACK.

Rebuild check: `make libbelfem_sparse.a` clean (exit 0), incl. `en_SolverEnums.cpp`.

## Files touched

- `src/sparse/cl_SolverSUPERLU.hpp` (new) — BLAS-prototype rename guard, member
  init, `free_factorization()` decl, comments, dead-decl removal.
- `src/sparse/cl_SolverSUPERLU.cpp` (new) — Tier 1/2/3 fixes; teardown split;
  RHS `BELFEM_ERROR`; tData leak.
- `config/linalg/config_superlu.cmake` (new) — `ARMA_DONT_USE_SUPERLU` guard.
- `src/sparse/en_SolverEnums.hpp` — default-solver `BELFEM_SUPERLU` branch.
- `src/sparse/en_SolverEnums.cpp` — `to_string(SUPERLU)`.
- `src/sparse/fn_matrix_type.hpp` — `matrix_type(SUPERLU) = CSC`.
- `src/sparse/fn_preferred_matrix_format.hpp` — `preferred_matrix_format(SUPERLU) = CSC`.
