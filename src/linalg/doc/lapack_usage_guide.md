# LAPACK Wrapper Usage Guide {#linalg_lapack_usage_guide}

**Date:** 2026-07-31
**Purpose:** How the unified backend-agnostic LAPACK wrappers in
`src/linalg/lapack/` work, the contracts they share, and how to add a new
routine.
**Module:** linalg

This is the maintainer-level guide for the wrapper layer. For the
user-facing `Vector`/`Matrix` API and everyday solver recipes, see
[linalg_usage_guide.md](linalg_usage_guide.md); this document covers what
sits underneath: the ABI glue, the scratch-buffer contracts, and the
conventions every `fn_*.hpp` in `lapack/` follows.

## Overview

Each LAPACK routine is wrapped exactly once, in a single header that works
for both matrix backends (Armadillo and Blaze) and all four LAPACK
datatypes (`float`, `double`, `std::complex<float>`,
`std::complex<double>`). Each header contains three layers:

1. **extern "C" prototypes** for the four flavors (`sgesv_` … `zgesv_`),
   declared to match the vendor headers that may share a translation unit.
2. **A dispatch template** `lapack::<name><T>` with a
   `static_assert( dependent_false<T> )` primary and four `inline`
   specializations that forward to the flavors.
3. **High-level wrappers** in namespace `belfem` taking `Matrix<T>` /
   `Vector<T>` operands. These are the functions user code calls.

The shared machinery — datatype glue, integer policy, leading dimensions,
and type traits — lives in `lapacktools.hpp`.

## Quick Reference

| Routine | Problem | Wrapper | Destroys / overwrites | Scratch |
|---|---|---|---|---|
| `gesv` | A·x = b, general square (LU) | `gesv( A, b\|B, Pivot [, abort] )` | A → LU factors, b → solution | Pivot ≥ n |
| `posv` | A·x = b, SPD / Hermitian PD (Cholesky) | `posv( A, b\|B [, abort] )` | A → Cholesky factor, b → solution | — |
| `getrf` | LU factorization | `getrf( A, Pivot [, abort] )` | A → LU factors | Pivot (grown) |
| `getri` | inverse from LU | `getri( A, Pivot, Work [, abort] )` | A → inverse | Work `Vector<T>` |
| `gels` | min ‖A·x − b‖, full rank | `gels( A, b\|B, Work [, abort] )` | A → QR/LQ, b → solution | Work `Vector<T>` |
| `gemm` | C := α·op(A)·op(B) + β·C (BLAS) | `gemm( A, B, C [, α, β, transa, transb] )` | C | — |
| `gesvd` | A = U·diag(S)·Vᵀ | `gesvd( A, S, U, VT, Work [, jobs, abort] )` | A | Work `Vector<real_t<T>>` (carved) |
| `geev` | A·v = λ·v, general | `geev( A, W, VL, VR, Work [, jobs, abort] )` | A | Work `Vector<real_t<T>>` (carved) |
| `gees` | Schur A = Z·T·Zᵀ/ᴴ | `gees( A, W, VS, Work, BWork [, select, jobs, sdim, abort] )` | A → Schur form T | Work carved + BWork `Vector<int_t>` |

Executable examples for every routine and datatype:
`tests/linalg/test_lapack.cpp` (typed gtest suite, 23 cases × 4 flavors).

## Shared Contracts

### Error model: BELFEM_ASSERT in, info out

- **Preconditions** (square shapes, matching dimensions, pivot lengths,
  valid job flags) are `BELFEM_ASSERT` — active in debug builds only, per
  the two-tier error policy in `doc/coding_philosophy.md`.
- **LAPACK failures** (`info != 0`) raise `BELFEM_ERROR` — unless the
  caller passes `AbortOnError = false`, in which case the wrapper returns
  `info` and the caller decides. Every solver wrapper returns an `int_t`
  `info` value (0 on success). This gives iterative schemes a hook to
  survive a failed solve: call with `false`, test the returned info, and
  react without tearing the process down.
- A suppressed **workspace-query** failure returns early with its info;
  the solve is not attempted with an invalid work size.
- `gemm` is the one exception: it is BLAS, has no `info`, returns `void`,
  and takes no flag.

```cpp
// hard failure ( default ): abort with an error box
gesv( tA, tB, tPivot );

// recoverable failure: react to a singular system
if ( gesv( tA, tB, tPivot, false ) != 0 )
{
    // e.g. damp the update and retry
}
```

### Pivot vectors

Pivots are always `Vector<int_t>` — the element width must match the
LAPACK integer (see Design Pillar 1). Two policies exist:

- `gesv` and `getri` **require** `Pivot.length() >= n` (asserted);
- `getrf` **grows** the pivot vector to `min( m, n )` itself, because it
  produces the pivot data that `getri` later consumes. Never modify a
  pivot vector between `getrf` and `getri`.

### Work buffers: three models

Scratch is always caller-owned and grow-only: the wrapper computes the
documented minimum, asks LAPACK for the optimal size (`lwork = -1` query)
only when the provided buffer is too small, and otherwise uses the full
buffer as-is. Reusing one Work object across calls therefore costs one
allocation total. Three models exist:

**(A) No Work** — `gesv`, `posv`, `getrf`, `gemm`. The routine needs no
workspace.

**(B) `Vector<T>` Work** — `getri`, `gels`. The work array uses
the operand type.

**(C) Single carved real Work** — `geev`, `gees`, `gesvd`. One
`Vector< lapack::real_t<T> >` carries everything, with no internal
allocation:

```
| LAPACK work array            | real scratch tail                    |
| lwork entries of T           | geev: 2n   gees: 2n (real)           |
| ( complex T: reinterpreted   |            n  (complex)              |
|   in place, 2 reals/entry )  | gesvd: 0 (real) / 5·min(m,n) (cplx)  |
```

The head is `reinterpret_cast<T*>( Work.data() )` for complex `T`; a
`tRealsPerT` constant (1 for real `T`, 2 for complex `T`) handles the length bookkeeping.
The tail doubles as `wr`/`wi` for the real flavors — whose eigenvalues
LAPACK returns as separate real and imaginary arrays, packed into the
complex `W` afterwards — and as LAPACK's actual `rwork` for the complex
flavors. Reading the optimal size from `Work(0)` works for both, since
the real part of the first work entry comes first in memory.

### Operands are owning, contiguous, column-major

Only owning `Matrix`/`Vector` objects may be passed — never Blaze views,
submatrices, or expression objects. The wrappers take raw `data()` pointers
plus a leading dimension; a view has neither a stable pointer nor the
stride the helper computes.

## Design Pillars

Each pillar states one contract. Enforcement lives in
`lapacktools.hpp` as `static_assert`s and preprocessor checks.

### 1. `int_t` IS the LAPACK integer

SCLS builds every third-party library in a given stack with a single
integer width, so BELFEM's `int_t` and the linked LAPACK's integer
coincide by construction. `USE_MKL_64BIT_API` switches both:
`BELFEM_INT64` (making `int_t` 64-bit) and `MKL_ILP64`;
`blaze_config.hpp` derives `BLAZE_BLAS_IS_64BIT` from the same macro so
Blaze's `blas_int_t` follows. **Wrong:** introducing a separate LAPACK
integer typedef or declaring prototypes with plain `int`. **Checked
by:** `sizeof( int_t ) == sizeof( blaze::blas_int_t )` /
`sizeof( MKL_INT )` static_asserts and the `MKL_ILP64 ↔ BELFEM_INT64`
paired `#error`.

### 2. Complex glue: `cplx_float_t` / `cplx_double_t`

extern "C" prototypes must match any vendor headers that share the
translation unit, or the compiler rejects the conflicting C-linkage
declarations. Blaze's clapack headers declare the complex flavors with
real-pair pointers (`float*`/`double*`) **in every Blaze TU — also when
MKL is the linked library** (nothing defines `INTEL_MKL_VERSION` there).
Therefore, `cplx_*_t` are `MKL_Complex8/16` only for
`BELFEM_MKL && !BELFEM_BLAZE`, and plain `float`/`double` otherwise. The
dispatch specializations use `reinterpret_cast` from `std::complex<T>*`,
which is layout-compatible by the C++11 array-access guarantee.
**Wrong:** declaring prototypes with `std::complex*` (conflicts with
Blaze) or with `MKL_Complex*` under Blaze.

### 3. `leading_dimension()` is a storage property

It is the stride between columns of the *stored* array: `n_rows()` under
Armadillo, `matrix_data().spacing()` under Blaze, whose columns are
SIMD-padded so the stride may exceed the row count. It **never depends on
a transposition flag** — `trans` only changes which logical dimension
LAPACK validates the stride against, which the physical stride satisfies
automatically. **Wrong:** passing `n_rows()` to LAPACK under Blaze
(silent corruption once padding kicks in), or switching the stride on
`trans` (wrong for every non-square matrix). Vectors are contiguous under
both backends; their overload returns the logical length.

### 4. One hidden Fortran length per character argument

Fortran compilers append a hidden by-value length argument for every
`character` dummy argument, after the regular argument list. All
mainstream Linux compilers (gfortran, ifort/ifx, flang) share this
convention; only the width varies. Under Blaze, `fortran_charlen_t` is
Blaze's own typedef (`int` for gcc ≤ 7, `size_t` since gcc 8); under
Armadillo BELFEM uses `size_t` unconditionally (`lapacktools.hpp:97-101`).
Every prototype carries one
`fortran_charlen_t` per char argument and every call passes `1`. LAPACK
itself never dereferences the lengths, so this is declaration
compatibility (Blaze declares them) and ABI honesty, not a fix for observed
breakage.

### 5. const dispatch, `const_cast` at the boundary

The dispatch templates take `const` pointers for the pure-input
arguments; the extern "C" prototypes stay non-const because they must
match the vendor declarations (Blaze declares non-const, MKL const —
they cannot both be matched, and Blaze is the one that shares BELFEM TUs).
The `const_cast` happens exactly once, inside the specializations.

## Routine-Specific Notes

- **`gesv`** — A must be square and is overwritten with its LU factors.
  `info > 0` means an exactly zero pivot.
- **`posv`** — symmetric positive definite for real `T`, **Hermitian**
  positive definite for complex `T` (not complex symmetric). No
  pivoting; `info > 0` names the first non-positive-definite minor.
- **`getrf` / `getri`** — a pair: `getrf` factorizes (and grows the
  pivot vector), `getri` inverts in place and requires exactly the
  factorization and pivots `getrf` produced.
- **`gels`** — full-rank least squares only (`info > 0` = rank
  deficient; a rank-tolerant `gelsy`/`gelsd` wrapper is the natural
  extension if a rank-deficient consumer appears). The RHS carries both
  input and output: it must be allocated with `max( m, n )` entries/rows
  — the underdetermined branch writes solution rows beyond `m`.
- **`gemm`** — dimension assertions are `op()`-aware: `m`, `n`, `k`
  derive from the transposition flags ('N'/'n', 'T'/'t', 'C'/'c'). With
  `beta == 0`, C is resized to fit; otherwise it must already have the
  shape of `op(A)·op(B)`.
- **`gesvd`** — the singular values `S` are **always real**, including for
  complex `A` (`Vector< lapack::real_t<T> >`). The `'O'` (overwrite) job
  flags are not supported. `U`/`VT` are sized per job flag; both `'N'`
  jobs still pass a valid `ld >= 1`.
- **`geev`** — the eigenvalues `W` are **always complex**
  (`Vector< lapack::cplx_t<T> >`). For real `T`, a complex conjugate
  pair `( W(j), W(j+1) )` has its eigenvectors LAPACK-packed across
  two consecutive columns: `v_j = VR(:,j) + i·VR(:,j+1)`,
  `v_{j+1} = conj( v_j )`. An eigenvalues-only overload skips `VL`/`VR`.
- **`gees`** — real flavors produce a **quasi**-triangular Schur form
  (2×2 blocks for conjugate pairs) with `A = Z·T·Zᵀ`; complex flavors a
  produce a triangular one with `A = Z·T·Zᴴ`. The sorting callback `select` has
  per-flavor arity — `int_t f( const T* wr, const T* wi )` for real,
  `int_t f( const std::complex<R>* w )` for complex — and returns
  `int_t`, not `bool`. **Real-flavor pair semantics:** a conjugate pair
  is selected as a whole if the callback is true for *either* member and
  counts as *two* toward `sdim`; complex flavors select independently.
  The complex callback pointer is bridged to the extern "C" typedef by a
  `reinterpret_cast` — sound on the SysV ABI, formally outside ISO C++
  (documented at the cast). `BWork` is the LOGICAL scratch
  (`Vector<int_t>`; Fortran LOGICAL width tracks the integer width),
  referenced only when sorting.

## Adding a New Routine

1. **Read the vendor declaration first.** If
   `/opt/scls/gcc/include/blaze/math/lapack/clapack/<name>.h` exists,
   your prototypes must match it exactly — argument types, real-valued
   arrays, `rwork` presence, and the number of trailing
   `fortran_charlen_t`. If Blaze does not declare the routine (gels,
   gees), BELFEM owns the prototype; still follow the family rules.
2. Count the character arguments; add one `fortran_charlen_t` per char to
   the prototypes and pass `1` per char in the specializations.
3. Check the Netlib docs for flavor asymmetries: real-valued outputs in
   complex flavors (`gesvd` S), extra arrays (`rwork`, `bwork`),
   eigenvalue splitting (`wr`/`wi`). Use `real_t<T>` / `cplx_t<T>` in the
   dispatch signature where types diverge, and decide the Work model:
   plain `Vector<T>` (model B) or a carved real buffer (model C) when a
   real scratch is needed in any case.
4. Dispatch layer: `static_assert( dependent_false<T> )` primary, four
   `inline` specializations, `const_cast` only there.
5. High-level wrapper: check preconditions with `BELFEM_ASSERT`, size outputs,
   `leading_dimension()` for every matrix operand (`ld >= 1` even for
   unreferenced ones), grow-only Work with the documented `lwork` floor
   and the `lwork = -1` query, an `AbortOnError` flag, and an `int_t` `info` return with
   early return on a suppressed query failure.
6. Add a Doxygen block on the wrapper; add typed tests over all four datatypes
   to `tests/linalg/test_lapack.cpp`, verified against a hand-computed
   reference, using odd matrix dimensions so Blaze's padded stride
   (`spacing() > n_rows`) is actually exercised.

## Pitfalls

| Pitfall | Consequence | Rule |
|---|---|---|
| Blaze view / submatrix passed to a wrapper | wrong stride, silent corruption | owning column-major objects only |
| `mkl_lapack.h`, the generic `lapack.h`, or un-shimmed `slu_ddefs.h` in the same TU | conflicting C-linkage declarations, compile error | never co-include vendor LAPACK/BLAS prototypes; SuperLU goes through the rename shim in `cl_SolverSUPERLU.hpp` |
| Forgetting that solvers destroy A | garbage on reuse | copy first if the original is still needed |
| `gels` RHS sized `m` in the underdetermined case | heap overrun | allocate `max( m, n )` |
| `Vector<int>` pivots | ill-formed under `BELFEM_INT64` | pivots are `Vector<int_t>` |
| Reading `sdim` as “number of true callback returns” (real gees) | off-by-pair | pairs count as two, see routine notes |

## See Also

- [linalg_usage_guide.md](linalg_usage_guide.md) — user-facing
  `Vector`/`Matrix` API and solver recipes (note: its `gesv`/`posv`
  sections predate the unified wrappers; the signatures documented here
  are authoritative).
- `tests/linalg/test_lapack.cpp` — executable examples, every routine ×
  every datatype.
- `doc/coding_philosophy.md` — the two-tier error policy and the
  column-major storage rationale.
- Netlib routine references: <https://www.netlib.org/lapack/explore-html/>
  (gesv, posv, getrf, getri, gels, gesvd, geev, gees; gemm under BLAS).
