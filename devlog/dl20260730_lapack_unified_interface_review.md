# Unified LAPACK Interface — Concept Review (gesv / lapacktools)

**Date:** 2026-07-30
**Purpose:** Three-way review (Claude + Codex + Grok) of Christian's concept
for a backend-agnostic LAPACK interface replacing the per-backend wrappers.
**Module:** linalg
**Session type:** read-only review; no source edits.

## Concept under review

`src/linalg/lapack/lapacktools.hpp` (new) provides a
`lapack::leading_dimension()` helper that answers the one question the two
backends answer differently — the stride between columns (`n_rows()` for
Armadillo, `matrix_data().spacing()` for padded Blaze). `src/linalg/lapack/`
`fn_gesv.hpp` (new) then implements a single generic `gesv` body: owned
extern "C" prototypes for `s/d/c/z gesv_`, a `lapack::gesv<T>` dispatch
template, and shared `belfem::gesv(Matrix, Vector|Matrix, Pivot)` wrappers.
The goal: once leading dimensions and datatypes are handled centrally, the
per-backend `fn_AR_*` / `fn_BZ_*` LAPACK wrappers collapse to one body.

## Verdict — unanimous

**Concept: GO.** All three reviewers agree the architecture is right; it
generalizes the empirically verified gels result (see
dl20260730_gels_least_squares_wrapper.md) and retires real debt: the
Armadillo gesv/posv paths pass `n` as LDA/LDB unconditionally
(`fn_AR_gesv.hpp:193-200`, `fn_AR_posv.hpp:187-195`) — only accidentally
correct because Armadillo never pads. The draft itself is a non-compiling
sketch (expected for a concept file) with one genuine design decision open.

## The design decision: the LAPACK integer — settled

The draft types LAPACK arguments as `int_t`. The auditors initially proposed
a dedicated `lapack_int_t` because the LAPACK integer width is a property of
the *linked LAPACK build* (LP64 vs ILP64) — and because of a macro landmine:
`typedefs.hpp:47` keyed `int_t` off `BELFEM_INT64`, while the only 64-bit
code path in CMake defined the never-consumed `BELFEM_I64`, so `int_t` and
`MKL_INT` would silently disagree under `USE_MKL_64BIT_API`.

Christian resolved it the other way, same session: the macros are now
unified — `USE_MKL_64BIT_API` appends `BELFEM_INT64` in the main
CMakeLists.txt and `BELFEM_I64` is retired. With SCLS guaranteeing that all
third-party libraries in a stack are built with one consistent integer
width (the point of SCLS), `int_t` and the linked LAPACK integer coincide by
construction; a separate `lapack_int_t` would only duplicate that guarantee.
The concept's sizeof `static_assert`s against `MKL_INT` /
`blaze::blas_int_t` remain as the backstop that the invariant holds.
Side effect worth recording: the OLD wrappers hardwire `int` in their
prototypes, so `USE_MKL_64BIT_API` builds were already ABI-broken before
this refactor — the unified interface is what makes ILP64 actually correct.

Consequences: pivot vectors become `Vector<int_t>` end-to-end, which needs
a caller sweep (`cl_Gradient.hpp`, `cl_BDF.hpp`, `fn_circle_from_points.hpp`,
`fn_create_beam_poly.hpp`, tests — all hold `Vector<int>`, identical today,
ill-formed under `BELFEM_INT64`), and the `"%i"` format on `info` must
width-match `int_t`.

## Consolidated blocker list (before replacing per-backend files)

1. Include-guard structure: the `#endif` at `lapacktools.hpp:24` closes the
   guard before the namespace; trailing `#endif` unmatched.
2. Guard/basename collision: `src/linalg/lapack/fn_gesv.hpp` reuses both the
   basename and guard `BELFEM_FN_GESV_HPP` of the public router
   `src/linalg/fn_gesv.hpp`; both directories are always on the include path
   (`config/Add_Library.cmake:7-8`), so a flat include is ambiguous and the
   second header is silently skipped. Found independently by both auditors.
3. Complex types: the `cplx_float_t`/`cplx_double_t` typedefs are reversed
   (`typedef existing new_name;`) and target keywords; the types exist
   nowhere else. Unanimous fix: drop the glue entirely and declare the owned
   prototypes with `std::complex<float/double>*` (existing AR precedent),
   plus layout `static_assert`s. MKL complex glue is unnecessary as long as
   BELFEM owns the prototypes and MKL headers are not in the same TU.
4. `#ifdef MKL` never fires — the build defines `BELFEM_MKL`
   (`config/linalg/config_mkl.cmake:12`). Also found independently twice.
5. `leading_dimension` must be `template<typename T>` (the draft accepts
   only `Matrix<real>` while `gesv<T>` is templated), with
   `std::max<lapack_int_t>` + explicit casts (`n_rows()`/`spacing()` return
   `size_t`), a range guard on the narrowing, `<algorithm>` included, and a
   return on all preprocessor paths.
6. Matrix-RHS overload asserts were copied from the Vector overload and call
   `B.length()` on a `Matrix` (no such member). Correct checks:
   `A.n_rows() == B.n_rows()`, A square, `Pivot.length() >= A.n_rows()`.
   Logical-size guards that Blaze padding can mask in release must be
   `BELFEM_ERROR`, not debug-only `BELFEM_ASSERT` (gels-thread rule).
7. extern "C" migration must be single-wave: the old prototypes in
   `fn_AR_gesv.hpp`, `fn_LAPACK_gesv/getrf/getri/gemm.hpp` use `int*`;
   C-linkage names share one declaration space, so old and new declarations
   cannot coexist in a TU once types diverge. Only live consumer of the old
   raw headers is `fn_invert_symmetric.hpp`.
8. Dispatch primary should fail at compile time
   (`static_assert(dependent_false<T>)`) rather than runtime `BELFEM_ERROR`.
9. Document the contract: owning, contiguous, column-major `Matrix`/`Vector`
   only — no Blaze views/submatrix expressions.

Cosmetics: missing commas between assert format strings and args, a
`complex<float>` printed where `complex<double>` is meant, unqualified
`complex<...>` in messages, `sgesv_` parameter named `l` instead of `ldb`
(inherited from the AR header).

## Migration order (agreed)

1. Fix blockers; land unified `gesv` behind the public router
   `src/linalg/fn_gesv.hpp`; delete AR/BZ gesv prototype blocks same-wave.
2. Re-unify `posv` (fixes the AR LDA/LDB debt) and `gels` (preserve
   `work_size()` and the `BELFEM_ERROR` logical-size rule) via the same
   `leading_dimension` helper.
3. Fold `fn_LAPACK_getrf/getri/gemm.hpp` into the `lapack_int_t` scheme and
   retire the old headers (port `fn_invert_symmetric.hpp`).

Note: raw-LAPACK paths bypass `blaze::gesv`'s internal invariant checks
(`fn_BZ_gesv.hpp:37-59`), so the wrapper guards become the only line of
defense — hence the ERROR-not-ASSERT policy on size checks.

## Implementation (same day, first wave: gesv)

Landed after Christian's go-ahead:

- `lapacktools.hpp` rewritten: fixed guards, `BELFEM_MKL` macro, namespace-
  scope `static_assert`s, templated `leading_dimension<T>()` with
  `std::max<int_t>` + range assert, `dependent_false` helper for the
  dispatch primaries.
- Unified `gesv` lives in `src/linalg/lapack/fn_gesv.hpp`; Christian moved
  the `posv`/`gels` routers there too and deleted the top-level routers, so
  `lapack/` is the single home for LAPACK wrappers and the basename/guard
  collision is resolved by there being only one `fn_gesv.hpp` on the
  include path.
- Guard policy (Christian, overriding the auditor recommendation): size and
  shape preconditions are `BELFEM_ASSERT` — the two-tier system is
  deliberate, debug mode + Valgrind catch logic bugs, release stays lean.
  `BELFEM_ERROR` remains only on LAPACK `info != 0`.
- Deleted `fn_AR_gesv.hpp`, `fn_BZ_gesv.hpp`, `fn_LAPACK_gesv.hpp` (only
  the router included the first two; nothing included the third).
- Pivot sweep `Vector<int>` → `Vector<int_t>` at every gesv call site:
  beam-poly helpers (+ signatures), circle_from_points, `cl_BDF`,
  `cl_Gradient`, materials Copper/Silver/Iron, `test_LinalgSolvers`,
  nonfree `cl_Gas`/`cl_CN_Scheme`, and kepler's `mHermitePivot` (found
  through the changed `create_fifth_order_beam_poly` signature).

**Q2 was resolved the opposite way from the auditor recommendation — by
Christian, correctly.** "Drop the MKL glue, use `std::complex` prototypes"
(unanimous auditor + Claude position) failed on the first Blaze compile:
Blaze's `clapack/gesv.h` puts `float*`/`double*` declarations of
`cgesv_`/`zgesv_` into every TU, so `std::complex` prototypes are a
conflicting C-linkage declaration. Claude's fallback (raw `float*`/
`double*`) fixed Blaze but would have silently broken the MKL path, where
the same symbols are declared with `MKL_Complex8/16` — and since MKL is not
in the current build config, no compile would have caught it. Christian's
original `cplx_float_t`/`cplx_double_t` typedef layer (MKL_Complex8/16
under `BELFEM_MKL`, plain float/double otherwise) is exactly what lets one
prototype line match every vendor's own declarations; the draft's only
defect there was the reversed `typedef` spelling, never the concept.

**Verification** (scratch program, both backends, Netlib LAPACK):
real 5×5 vector-RHS — under Blaze with `spacing() = 6 > n_rows = 5`, so the
padded-stride path was genuinely exercised; real 3×3 with nrhs = 2; complex
3×3 through `zgesv`. All pass, max error ≤ 9e-16, identical across
backends. Syntax-checked with real build flags: `test_LinalgSolvers.cpp`,
`cl_BDF.cpp`, materials Copper/Silver/Iron, helper headers under both
backends. Not checked: `cl_Gradient`/`cl_Surface`/nonfree TUs — the build
tree's per-library `flags.make` are stale (old MKL-era defines; today's
CMake option changes need a cmake re-run); their sweep pattern is identical
to the verified files.

## Second wave: posv + gels (Christian's unified drafts, hardened same day)

Christian moved unified `fn_posv.hpp` / `fn_gels.hpp` into
`src/linalg/lapack/` alongside gesv. Fixes applied on top:

- `fn_gels.hpp`: three missing `reinterpret_cast< cplx_*_t * >` in the
  complex dispatch (`work` in cgels, `b` and `work` in zgels);
  `fn_BZ_polyfit.hpp` re-pointed from the deleted `fn_BZ_gels.hpp` to
  `fn_gels.hpp`.
- `fn_posv.hpp`: the posv family takes a `char * uplo`, and the Fortran ABI
  appends a hidden length argument after the regular list — Blaze's
  `clapack/posv.h` declares it explicitly (`blaze::fortran_charlen_t`), so
  our prototypes must match or every Blaze TU trips a conflicting-C-linkage
  warning (= error under the build's `-Werror`). New `fortran_charlen_t`
  typedef in `lapacktools.hpp` (Blaze's under `BELFEM_BLAZE`, else
  `size_t`), trailing argument on all four prototypes, specializations pass
  `1`. gesv/gels have no char arguments — posv is the only affected family.
- `fn_posv.hpp` latent bugs: Vector overload called
  `leading_dimension( B )` on a `Vector` (no such overload — invisible
  until first instantiation since templates are lazy; now `ldb =
  max(1,n)`, vectors are contiguous under both backends); Matrix overload
  still passed `&n` as LDA/LDB (the exact unpadded-only debt this refactor
  retires) — now `leading_dimension` for A and B, `int` locals widened to
  `int_t`; dispatch primary now `static_assert( dependent_false<T> )`
  consistent with gesv.

**const in the raw prototypes — settled non-const.** Fortran LAPACK never
modifies n/nrhs/lda/ldb/uplo/trans (documented "unchanged on exit"), and C
linkage ignores constness, so `const int_t*` prototypes would be legal and
ABI-identical. But declaration consistency rules them out: Blaze's clapack
headers declare these symbols with non-const `blas_int_t*` in every Blaze
TU, and two extern "C" declarations of one symbol with different parameter
types are ill-formed — the same rule that produced today's complex-type and
charlen conflicts. (MKL's own headers declare const; Blaze defers to MKL
via `#if !defined(INTEL_MKL_VERSION)` when MKL headers are present, so
whichever vendor header is in the TU dictates.) Non-const + `const_cast`
at const call sites (the old `fn_LAPACK_gemm.hpp` pattern) is the only
convention compatible with all vendors; the high-level wrappers own their
locals as non-const, so in practice no casts are needed at all.

Verification after this wave: scratch program extended with SPD posv
(vector + matrix RHS), compiled `-Wall -Werror`, **zero warnings**, all
cases pass on both backends (max err ≤ 9e-16); gesv+posv+gels+polyfit
headers syntax-check together under both backends.

## Third wave: gemm / getrf / getri audit (Christian's refactor, same day)

Christian folded gemm/getrf/getri into `src/linalg/lapack/` with the
const-dispatch + `const_cast`-at-the-boundary pattern and asked for an
audit, specifically the gemm transpose handling. Findings, all fixed and
numerically verified:

**gemm — the transpose suspicion was correct, twice over:**
- `m`/`n`/`k` and all three asserts assumed `'N','N'`. Correct rule: `m` =
  rows of op(A), `n` = cols of op(B), `k` = cols of op(A) = rows of op(B),
  each switching on its trans flag.
- The trans-aware `leading_dimension( A, transa )` rested on a
  misconception: **LDA is a property of the stored array only** — the
  stride between columns as A sits in memory — and never changes with
  trans. Trans only changes which logical dimension LAPACK validates the
  stride against (`lda >= m` for 'N', `>= k` for 'T'), which the physical
  stride satisfies automatically. The Armadillo branch returning
  `n_cols()` under 'T' passed a wrong stride for any non-square matrix;
  the Blaze branch already (correctly) ignored the flag. The trans
  parameter is removed; with it the unused-parameter pragma block became
  unnecessary and was dropped entirely.
- Latent compile errors: `std::complex< float * >` (pointer inside
  complex) in both complex specializations — they could never match, so
  complex gemm fell to the runtime-error primary; `beta` passed by value
  where `const T*` is expected; unqualified `gemm(...)` call could not
  reach `lapack::gemm` from namespace `belfem`.

**getri — wrapper never inverted:** the final call after the workspace
query was `lapack::getrf` (factorization again) instead of
`lapack::getri`, dropping Work/lwork; `A.data` missing parentheses in the
query call; the Pivot auto-resize replaced by an assert (getri *consumes*
getrf's pivot output — a fresh pivot vector is meaningless); square-A
assert added; contract documented (A must hold the getrf factorization).

**getrf:** same `std::complex< float * >` typos; unqualified `getrf` call
→ `lapack::getrf`; sign-compare in the pivot resize; primary's error
message said "getri".

**Cross-cutting:** dispatch primaries unified on
`static_assert( dependent_false<T> )` (gemm/getrf/getri/gels); gels guard
renamed from stale `BELFEM_FN_AR_GELS_HPP`; `std::max( 1, <int_t> )`
ambiguities pinned with `std::max<int_t>` (breaks under BELFEM_INT64);
`%i` format args cast to int; `fn_invert_symmetric.hpp` passed `info` by
value where `int_t*` is expected (now `int_t` + `&info`) and
`cl_Material_Alloy` pivots swept to `Vector<int_t>`.

**Intel pragmas in lapacktools.hpp:** intent was the new LLVM-based Intel
compiler (icx/icpx) — which defines `__clang__`, so the
`#pragma clang diagnostic` branch already covers it; the
`#elif BELFEM_INTEL` branch used classic-icc `#pragma warning` syntax and
was unreachable for icx (and the push/pop guards keyed on different macro
families, risking imbalance). Moot after the trans-parameter removal: the
whole pragma block is gone.

Verification: scratch suite extended to 13 cases — gemm in all four
transpose combinations on non-square operands (2×3 · 3×4) plus an
`alpha=2, beta=3` accumulate case, getrf+getri inversion round-trip,
gels overdetermined — ALL PASS on both backends, `-Wall -Werror`, zero
warnings, max error ≤ 2e-15. `cl_Material_Alloy.cpp` and
`test_LinalgSolvers.cpp` syntax-check clean.

## Fourth wave: gesvd audit

Christian added `fn_gesvd.hpp` for completeness (no consumer yet). The
draft had the two classic gesvd traps, both ABI-level, verified against
Blaze's `clapack/gesvd.h` declarations:

- **The singular values `s` are ALWAYS real** — `float*`/`double*` even in
  the complex flavors. The draft declared them `cplx_*_t` and typed the
  high-level `S` as `Vector<T>`.
- **`cgesvd_`/`zgesvd_` take an extra real `rwork` array** (≥ 5·min(m,n))
  between `lwork` and `info`. Omitting it shifts the argument list — the
  library reads the `info` slot as the rwork pointer. The real flavors do
  NOT have this argument, so the unified dispatch signature carries an
  `rwork` parameter that the s/d specializations ignore.
- Plus two trailing `fortran_charlen_t` for the two job chars (posv
  precedent, Blaze declares them).

New `lapack::real_t<T>` trait in `lapacktools.hpp` (T for real, value_type
for `std::complex<T>`) types `s`/`rwork` in the dispatch and `S` in the
high-level wrapper.

High-level wrapper defects fixed: missing `void` return type; `jubu`/
`jubvt`/`'a '` typos; U rows were `n` instead of `m`, VT dims transposed
for 'S' (correct: U is m×m 'A' / m×mn 'S', VT is n×n 'A' / mn×n 'S');
`ldu`/`ldvt` of 0 for the 'N' jobs (LAPACK requires ≥ 1 always — now
`max(1,·)` semantics via the 1-default); the workspace "query" passed
`&query` (an `int_t*`) as the work array with `lwork = 1` instead of the
`lwork = -1` protocol; `info` never declared; always-query replaced by the
family's grow-only caller-owned Work pattern; internal real `tRWork`
allocated only for complex T; stray `#include <dense/DenseMatrix.hpp>`
removed (it resolved to an unrelated library under /opt/scls).

Note: an editor stale-buffer save restored the draft once mid-audit
(caught by the Armadillo compile picking up the reverted file); the
corrected version was re-applied.

Verification: suite now 16 cases — gesvd real 3×2 'A' with known singular
values (5, 3) and U·diag(S)·VT reconstruction, real 4×3 'S' economy with
shape checks, complex 3×3 through zgesvd exercising the rwork path — ALL
PASS both backends, `-Wall -Werror`, zero warnings.

## Fifth wave: final 3-AI audit, test suite, AbortOnError

**Test suite.** New `tests/linalg/test_lapack.cpp` (registered in the
`test_linalg` target): a gtest typed suite over the four LAPACK datatypes
( float, double, complex<float>, complex<double> ) — 13 tests × 4 types =
52. Coverage: gesv (vector/matrix RHS + singular no-abort), posv
(Hermitian-PD, both RHS forms), gels (overdetermined + underdetermined
min-norm with oversized-B contract), gemm (plain, double-transposed,
α=2/β=3 accumulate), getrf+getri round trip, gesvd ('A' reconstruction
with shape/order checks + 'S' economy). One data set serves all four
flavors via a Scalar<T> factory that drops imaginary parts for real types.
Writing the suite caught a live bug: `beta == 0.0` in gemm does not
compile for complex<float> — now `static_cast<T>(0.0)`.

**Grok final audit:** current Netlib path clean; all prior defects
confirmed fixed in-tree. Actionable residuals applied:
- R1 (HIGH): under Blaze, Blaze's clapack protos stay active even with MKL
  (`mkl_cblas.h` never defines INTEL_MKL_VERSION), so the MKL_Complex glue
  now applies only for `BELFEM_MKL && !BELFEM_BLAZE`; under Blaze the
  real-pair convention always matches Blaze's declarations. Corroborated
  independently by Codex against a real oneAPI 2025.2 install.
- R10: `sizeof(int_t) == sizeof(MKL_INT)` static_assert backstop.
- R6: gemm dimension logic accepts lowercase 'n' job chars.
- R5: gels B-size guards aligned to the settled ASSERT policy.
- R9: `fn_invert_symmetric` info checks upgraded ASSERT → ERROR (runtime
  failures, family rule).
- R3 (partial): range assert on Vector leading_dimension; the per-wrapper
  n/nrhs casts remain unguarded (accepted residual, matters only past
  2^31 rows).
- R4: last `Vector<int>` pivots in `tests/math/test_TensorKernels.cpp` and
  `tests/old/math/tensor/fn_invert_symmetric.cpp` swept to `Vector<int_t>`.
Documented-only residuals (no code change): R2 charlen-vs-MKL-header arity
if `mkl_lapack.h` ever shares a TU; R7 SuperLU proto clash if the
`cl_SolverSUPERLU.hpp` rename shim were dropped; R8 owning/contiguous
contract not restated on every public entry point.

**AbortOnError** (Christian's design): high-level solver wrappers take
`const bool AbortOnError = true` and return `int_t` info, so iterative
callers (planned Anderson stabilization) can suppress the fatal error and
react to `info != 0`. Completed the pattern across the family: the
Vector-RHS overloads of gesv/posv/gels now match their Matrix-RHS
siblings (flag + info return), and both gels workspace-query errors are
gated with an early `return info` on suppressed failure (getri/gesvd
already had this). gemm keeps `void` (BLAS, no info). New regression
test: singular gesv with `AbortOnError = false` must return info != 0
without aborting.

**Codex final audit** (empirical, against the oneAPI 2025.2 install and
Netlib contracts): all prior fixes confirmed in-tree; LDA logic, workspace
queries, gesvd real-S/rwork and job sizing confirmed correct. Its one
actionable blocker is now fixed: under `BELFEM_INT64` + Blaze, nothing set
`BLAZE_BLAS_IS_64BIT` (blaze/config/BLAS.h defaults it to 0), so
`blaze::blas_int_t` stayed 32-bit and the `lapacktools.hpp` static_assert
fired — the backstop catching a real configuration gap, exactly as
designed. `blaze_config.hpp` now defines `BLAZE_BLAS_IS_64BIT` from
`BELFEM_INT64`, keeping Blaze's BLAS integer in sync with `int_t`
(Codex verified the fix via `-fsyntax-only`). Stale `Vector<int>` pivot
mentions in `tests/doc/tests_04_tensor.md`, `tensor_usage_guide.md` and
`linalg_usage_guide.md` updated.

**Accepted constraints** (documented in the exchange thread, not fixable
while matching Blaze's declarations): `mkl_lapack.h` must never share a TU
with these wrappers (const-qualified MKL prototypes and missing char
lengths are conflicting declarations); the same holds for the generic
`lapack.h`/BLAS `fortran.h` under `LAPACK_FORTRAN_STRLEN_END` and for
un-shimmed `slu_ddefs.h` (`cl_SolverSUPERLU.hpp`'s rename shim keeps
SuperLU safe). Char-length policy (superseded 2026-07-31, see tenth wave): now present
on EVERY routine with character arguments.

All 52 tests pass under both backends, `-Wall -Werror`, zero warnings.

## Sixth wave: test_linalg made Blaze-clean (pre-existing defects)

Making the whole `test_linalg` target build and pass under the current
Blaze config surfaced a chain of pre-existing Blaze-side defects, all
fixed:

- `fn_append.hpp` Blaze branch was written with std::vector semantics
  (`clear`/`reserve`/`insert` — none exist on `blaze::DynamicVector`);
  now `resize( old + add, true )` (preserve) + tail `std::copy`, with a
  debug assert against self-append (the pointer would dangle across the
  realloc). Both overloads.
- `test_Matrix.cpp` wrote through row/col views with `operator()`, which
  Blaze's `Row`/`Column` views don't provide — switched to `operator[]`,
  supported by both backends' view types.
- **Real latent backend divergence:** `cl_BZ_Matrix::submat()` forwarded
  `( firstRow, firstCol, lastRow, lastCol )` straight into
  `blaze::submatrix`, whose 3rd/4th arguments are ROW/COLUMN COUNTS — so
  under Blaze, `submat` returned the wrong shape for every call while the
  Armadillo passthrough used inclusive last indices. Now converts with
  `last - first + 1` (Armadillo convention documented as the interface
  contract). No production caller existed yet; the test was the first
  consumer to notice.
- `Matrix.CapacityMatchesProduct` asserted `capacity() == rows*cols`, an
  Armadillo-only invariant — Blaze's capacity includes the SIMD column
  padding (3×5 → 20). Test now checks `>= rows*cols` (physical-allocation
  meaning kept; no production caller uses `Matrix::capacity()`).

Full `test_linalg` (237 tests incl. the 52 LAPACK ones) passes under
BOTH backends.

## Seventh wave: MKL sparse multiply port (link failure under oneMKL 2026)

Turning MKL on broke the `hphiTrun` link: undefined `mkl_dcsrmm_` /
`mkl_dcscmm_` from `cl_SpMatrix::multiply`. Root cause pinned precisely:
these are the NIST-style sparse BLAS routines Intel deprecated years ago —
oneMKL **2025.2 still ships them**, but `/opt/intel/oneapi/mkl/latest`
points to **2026.1**, where they are removed from the libraries (only
`mkl_internal_dcsrmm_` remains in the layered interface lib). Not a
link-order problem — the deprecation cliff.

Fix: ported both `BELFEM_MKL` branches of `SpMatrix::multiply` to the
supported Inspector-Executor sparse BLAS via one file-local helper
(`mkl_sparse_multiply` in `cl_SpMatrix.cpp`): wraps the existing
CSR/CSC arrays in a `sparse_matrix_t` handle without copying
(`mkl_sparse_d_create_csr/csc`, index base read from `mPointers[0]`),
calls `mkl_sparse_d_mv` (alpha/beta/transpose preserved; the simple
overload passes 1/0/no-trans), destroys the handle. Status checks are
`BELFEM_ERROR`. `fspblas.hpp`'s dead `mkl_dcscmm_`/`mkl_dcsrmm_`
prototypes replaced by `#include <mkl_spblas.h>`. `int_t*` → `MKL_INT*`
via `reinterpret_cast` (same width by the lapacktools static_asserts;
under ILP64 `int64_t` = `long` vs `MKL_INT` = `long long` are distinct
types of equal width, so the cast is required and safe).

Verified: `-fsyntax-only -Werror` clean against the current MKL+Blaze
build flags; all four IE symbols (`create_csr/csc`, `d_mv`, `destroy`)
confirmed exported by 2026.1's `libmkl_intel_lp64.so`. Repo swept: no
other deprecated NIST-style MKL calls remain. Known trade-off: the handle
is created per call — if profiling ever shows this in a hot loop, cache
it as a member and invalidate on matrix mutation.

## Eighth wave: sparse multiply hybrid + FSPBLAS retirement

Christian's design, implemented: the plain `y = A·x` overload ALWAYS uses
the portable Fortran kernels (`matvec_csr`/`matvec_csc` in
`splinalg.f90`) on every build — identical numerics everywhere, covers the
ARPACK reverse-communication loop; the extended
`y = α·op(A)·x + β·y` overload uses MKL's Inspector-Executor
`mkl_sparse_multiply` when `BELFEM_MKL` is set (fused alpha/beta, native
transpose) and otherwise falls back to the kernels with the `mSwap`
beta-emulation — `mSwap` stays, it is what makes beta possible since the
kernels overwrite y.

FSPBLAS retired: `USE_FSPBLAS` option removed from CMakeLists.txt and
`config_mkl.cmake` (incl. the `-lfspblas` link line and the MKL/FSPBLAS
mutual-exclusion check); the `dcscmm_`/`dcsrmm_` prototypes deleted from
`fspblas.hpp`, whose own-kernel prototypes are now unconditional.

**Latent bug fixed in the portable extended path:** it silently ignored
`aTransposedFlag` — the existing `MultiplyTransposed` test used a
symmetric matrix (`Aᵀx == A·x`) and could not catch it. Now handled via
the CSR↔CSC duality: `op(A)` keeps the same arrays and calls the sibling
kernel with swapped dimensions (a CSR matrix's arrays read as CSC are
exactly `Aᵀ`; the pointer dimensions line up by construction). Three new
tests lock it in: non-symmetric rectangular (3×5) transpose for CSR and
CSC, plus an `α=2, β=3` accumulate case.

Verified: full sparse suite (40 tests) passes in BOTH flavors — portable
kernels (Netlib defines) and MKL (linked and executed against the real
oneMKL 2026.1 `mkl_intel_lp64/sequential/core`).

## Ninth wave: geev (2026-07-31)

Christian added `fn_geev.hpp`; three-way audit (G1–G8, all confirmed by
Codex + Grok) and rewrite. geev's signature asymmetry is the sharpest in
the family: the real flavors return eigenvalues as separate `wr`/`wi`
arrays with conjugate-pair eigenvectors packed across consecutive VR/VL
columns and have no rwork; the complex flavors take one complex `w` PLUS
a real `rwork(2n)`. Blaze declares all four with two trailing
`fortran_charlen_t`.

Key fixes: `rwork` was declared `int_t*` (ABI break); missing charlens;
the complex<double> specialization called `dgeev_` and the zgeev one
called `cgeev_` (copy-paste); `k < n` with `n` a pointer; non-`inline`
specializations; the high-level wrappers had a dozen compile/contract
errors (return type `info`, undeclared locals, `ldvr` from VL, `jobvl`
lost, W never sized). New `lapack::cplx_t<T>` trait complements
`real_t<T>`.

**Scratch design (Christian's ruling, overrides the gesvd-style internal
tRWork):** ONE caller-owned `Vector<real_t<T>> Work` — LAPACK's work array
in its head (for complex T reinterpreted in place, two reals per entry),
the 2n real scratch in its tail (wr/wi for real flavors — packed into the
complex W after the solve, query-guarded — and LAPACK's rwork for complex
flavors). No allocation inside the wrapper; grow-only with per-flavor
lwork floors. `fn_gesvd.hpp` still has an internal tRWork from before this
ruling — to be aligned when next touched.

Tests: 2 new typed cases × 4 datatypes — upper-triangular real spectrum,
and a rotation-block matrix (eigenvalues +i, −i, 2) verifying the
conjugate-pair eigenvector packing through the A·v = λ·v residual. Suite
now 60 tests, all passing on both backends, `-Wall -Werror` clean.

## Tenth wave: hidden char lengths unified family-wide (2026-07-31)

Christian asked whether gels/gemm/gesvd need the hidden lengths too, and
whether the mechanism is gfortran-specific. Findings:

- gesvd/posv/geev already had them (Blaze forces); gels (1 char) and gemm
  (2 chars) did not. The old rationale for omitting them on gemm — "match
  SuperLU's prototypes" — is void: `slu_ddefs.h` declares `dgemm_` with
  const-qualified pointers (`const char*`, `const int*`), so the
  declarations conflict with ours regardless of char lengths; the
  `#define dgemm_ belfem_slu_dgemm_` rename shim in `cl_SolverSUPERLU.hpp`
  is the actual co-include protection. Nothing else in-tree declares gels.
- gels and gemm prototypes and dispatch calls now carry the trailing
  `fortran_charlen_t` arguments ( `1` per char ), completing the rule:
  every character dummy argument has its hidden length, everywhere.
- Compiler scope: the trailing hidden length is NOT gfortran-specific —
  on Linux x86-64 gfortran, Intel ifort/ifx and LLVM flang all use the
  same convention (length appended by value after the regular argument
  list); the only variation Blaze's `fortran_charlen_t` encodes is the
  WIDTH (int for gfortran ≤ 7, size_t since gfortran 8; Intel uses
  size_t). MKL's C-callable entry points ignore the lengths — passing
  them is harmless under the SysV ABI (extra trailing by-value args).

Verified: full 60-test suite passes on both backends against the
gfortran-built Netlib libraries; `cl_Material_Alloy` (invert_symmetric →
raw gemm dispatch) syntax-checks clean — dispatch signatures are
unchanged, the lengths live only inside the specializations.

## Eleventh wave: gees — Schur decomposition with sort callback (2026-07-31)

Christian scaffolded `fn_gees.hpp` (prototypes and the four per-flavor
SELECT typedefs were correct, including rwork in c/z and both charlens);
the dispatch and wrapper were completed here. gees adds the family's first
callback: SELECT for eigenvalue sorting, with per-flavor arity (two
`const real*` for s/d, one `const complex*` for c/z), plus the LOGICAL
`bwork` array — `Vector<int_t>& BWork` in the wrapper, since Fortran
LOGICAL width tracks the default integer width in consistent builds
(gfortran `-fdefault-integer-8`, `MKL_INT`), referenced and grown only
for `sort='S'`.

Design: user-facing `lapack::gees_select_t<T>` trait keeps callbacks in
`std::complex` terms; the complex dispatch reinterpret_casts the function
pointer to the extern "C" typedef — both auditors rate this practically
universal on the SysV ABI and formally outside the standard (documented
at the cast; a C-linkage trampoline with a file-static callback slot is
the upgrade path if ever needed — viable because BELFEM is
single-threaded per rank). `nullptr` select iff `sort='N'`, asserted.
Single real Work buffer per the geev ruling, with a per-flavor tail:
2n reals (wr/wi) for s/d, n reals (rwork — gees, unlike geev, needs only
n) for c/z. `sdim` via optional out-pointer; info = n+1/n+2 reorder
semantics documented. Blaze declares no gees (only gges), so BELFEM owns
the prototypes.

Audits (Codex + Grok): SHIPPABLE, zero contract or Work-carving defects;
all their substantive residuals were documentation and one test gap, all
applied — factorization stated as `Z·T·Zᵀ` (real) vs `Z·T·Zᴴ` (complex),
and the real-flavor SELECT pair semantics: a conjugate pair is selected
if EITHER member returns true and counts as TWO toward sdim. The new
`GeesSortedConjugatePair` test locks exactly that asymmetry: selecting
only `imag > 0` yields sdim == 2 for real flavors, sdim == 1 for complex.

Armadillo safety (both eigen-wrappers): SCLS Armadillo is an
ARMA_USE_WRAPPER build (`wrapper2_*` symbols), so it never declares raw
`dgeev_`/`dgees_` — no C-linkage conflicts; only a hypothetical
ARMA_DONT_USE_WRAPPER build would need re-checking.

Suite: 72 tests (18 × 4 datatypes), all passing on both backends,
`-Wall -Werror` clean. Family: ten routines.

## Twelfth wave: documentation (2026-07-31)

- All thirteen `belfem::` wrapper functions carry Doxygen blocks
  (`@brief` / `@param[in,out]` / `@return`) documenting in/out semantics,
  destruction of A, scratch contracts and info meanings, in the house
  style sampled from `cl_JcFunction.hpp`.
- New `src/linalg/doc/lapack_usage_guide.md` — maintainer-level guide:
  quick-reference table, shared contracts (error model, pivot policy, the
  three Work models incl. the carved single buffer), the five design
  pillars as copyable one-rule contracts, per-routine notes, the
  add-a-routine checklist, pitfalls table. Scope/structure consulted with
  Codex + Grok beforehand (both: cross-reference `linalg_usage_guide.md`
  rather than duplicate; audience = future wrapper authors; all ten
  routines in the quick reference). Registered in `src/linalg/doc/
  README.md`.
- Consultation flagged two source inconsistencies, fixed before
  documenting: the matrix-RHS `gels` workspace-query error was not gated
  by `AbortOnError` (now gated + early return, like every sibling), and
  the `gesvd` job asserts accepted `'O'` although the wrapper does not
  implement overwrite semantics (restored to `'A'/'S'/'N'`).
- Known deviation documented as-is: `gesvd` still allocates its complex
  rwork internally (predates the single-buffer ruling; align when next
  touched). Stale `gesv`/`posv` sections in `linalg_usage_guide.md` are
  flagged in the new guide's See-Also rather than silently rewritten.
- Suite re-verified after the doc pass and the two fixes: 72/72 both
  backends, `-Wall -Werror` clean.

## Records

- Exchange thread: `tmp/ai_exchange/lapack_unified_interface.md` (Claude
  initial review D1–D11 / Q1–Q5, Grok D12–D19, Codex M1–M5, resolution,
  implementation + Q2 correction).
- Follow-ups: re-unify `posv` (fixes the AR LDA debt) and `gels` (keep
  `work_size()` + `BELFEM_ERROR` rule) via `leading_dimension`; fold
  `fn_LAPACK_getrf/getri/gemm` into the scheme (consumer:
  `fn_invert_symmetric.hpp`).
