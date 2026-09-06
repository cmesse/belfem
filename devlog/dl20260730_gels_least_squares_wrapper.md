# LAPACK `gels` least-squares wrapper for both linalg backends

**Date:** 2026-07-30
**Purpose:** Record the addition of `fn_gels.hpp` (real + complex, Armadillo and Blaze) and the replacement of the Blaze `polyfit` normal-equations solve with a direct least-squares solve
**Module:** `src/linalg`

## The question that started it

Blaze ships no `gels.h` — its LAPACK convenience layer covers `gesv`, `posv`,
`geqrf`, `ormqr` and friends, but not least squares. The apparent options were
to copy the Blaze matrix into a temporary column-major buffer, or to implement
least squares by hand.

Neither is necessary. BELFEM configures Blaze with `columnMajor` storage
(`src/linalg/blaze/blaze_config.hpp:55`, consumed at
`src/linalg/blaze/cl_BZ_Matrix.hpp:26`), so the only difference from Armadillo
is that Blaze pads each column for SIMD alignment — the inter-column stride is
`spacing()`, not `n_rows()`. That is exactly what LAPACK's `LDA` argument is
for. Netlib `DGELS` validates only `LDA >= max(1,M)` and `LDB >= max(1,M,N)`,
and threads the supplied leading dimensions unchanged through the internal
`DGEQRF` / `DORMQR` / `DTRTRS` (and the LQ branch) calls. Blaze's own LAPACK
wrappers are the existence proof: they pass `spacing()` as `lda` throughout.

Both auditors (Codex, Grok) confirmed this independently at high confidence.
The claim is now also confirmed empirically — see the verification section.

## What was added

| File | Contents |
|------|----------|
| `src/linalg/fn_gels.hpp` | backend dispatcher, matching `fn_gesv.hpp` / `fn_posv.hpp` |
| `src/linalg/armadillo/fn_AR_gels.hpp` | `s/d/c/zgels_` declarations, typed dispatch, `Matrix`-level overloads |
| `src/linalg/blaze/fn_BZ_gels.hpp` | same, with `lda`/`ldb` taken from `spacing()` |

Public API, identical on both backends:

```cpp
gels( Matrix< T > & A, Vector< T > & B, Vector< T > & Work );   // single rhs
gels( Matrix< T > & A, Matrix< T > & B, Vector< T > & Work );   // multiple rhs
```

`T` may be `float`, `double`, `std::complex<float>` or `std::complex<double>`.
`A` is destroyed (overwritten by the QR or LQ factorization). `B` carries the
right hand side in and the solution out, so it must be allocated with at least
`max(M,N)` rows. `Work` is caller-owned scratch, grown to the LAPACK-reported
optimal size on first use and reused thereafter — repeated calls from a driver
or optimizer do not reallocate.

Because LAPACK reports the optimal work size in `work[0]`, which stays real
valued for the complex flavors, the size is extracted through a small
`lapack::work_size()` overload set rather than a direct cast.

## Two traps worth remembering

**The leading dimension of `B` must come from the allocation, never from a
formula.** Computing `ldb = max(M,N)` satisfies LAPACK's numeric constraint
while lying about the buffer: in the underdetermined case `DGELS` writes rows
`M+1..N` of `B`, so a `B` allocated with only `M` rows is overrun. Both
overloads now read `ldb` from `B.length()` / `B.n_rows()` (Armadillo) or
`B.matrix_data().spacing()` (Blaze).

**The size check on `B` is `BELFEM_ERROR`, not `BELFEM_ASSERT`.** Under Blaze,
column padding can make the *stride* large enough to hide an under-sized
*logical* `B`, so a debug-only check would leave a release build silently
producing garbage. The check therefore tests logical rows and stays active in
release.

The inverse of the padding question is the real footgun: passing `n_rows()` as
`LDA` on a padded Blaze matrix makes column `j+1` start inside column `j`'s
padding. Padding *values* are harmless as long as the strides are honest.

## `polyfit` on the Blaze backend

The Blaze `polyfit` accumulated `Vᵀ V` and `Vᵀ y` directly and solved with
`blaze::posv` (Cholesky). Forming the normal equations squares the condition
number, which for a Vandermonde system is already large. It now builds the tall
Vandermonde matrix and solves it directly with `gels`.

Measured on a degree-8 fit over 60 samples with `x ∈ [2,5]`, comparing maximum
relative coefficient error against the known exact coefficients:

| method | max relative error |
|--------|--------------------|
| normal equations + Cholesky (previous) | 7.1e+01 |
| `gels` (current) | 1.5e-07 |

The previous path was not merely less accurate at this degree — it returned
values with no correct digits.

The body now uses the template parameter `T` throughout instead of hard-coding
`real`, matching the Armadillo sibling. Cost of the change is memory: the tall
`numSamples × (n+1)` Vandermonde replaces the small `(n+1) × (n+1)` normal
matrix. `polyfit` is a setup-path function, so this is an acceptable trade.

The Armadillo `polyfit` was deliberately left alone — it delegates to
`arma::polyfit`, which is already SVD-based and additionally tolerates rank
deficiency. Switching it to `gels` would have been a robustness regression.

## Choice of routine

`gels` assumes full rank and returns `INFO > 0` (a zero diagonal entry in the
triangular factor) rather than computing a minimum-norm solution when the rank
drops. This is documented at both overloads. If a rank-deficient consumer
appears, `gelsy` (pivoted complete orthogonal factorization) or `gelsd` (SVD)
is the right addition — as a sibling routine, not as a replacement.

## Verification

Both backends were compiled and run against the same test program. All cases
pass on both:

- overdetermined real fit, exact recovery of a linear model
- underdetermined `2×3` system, minimum-norm solution, oversized `B`
- multiple right hand sides
- complex (`zgels`) fit, exact recovery of complex coefficients
- `polyfit` recovery of a known quadratic

The Blaze run exercises the padded-stride path end to end, which is the
empirical confirmation that no copy is required.

Two-translation-unit link check confirms the specializations are emitted as
weak symbols, so the header is safe to include from multiple sources. All five
existing `polyfit` consumers (`cl_Material_Copper/Lead/Silver.cpp`,
`cl_GM_Helmholtz.cpp`, `cl_GM_EoS_Hydrogen.cpp`) compile unchanged on both
backends.

Not covered: `sgels` and `cgels` were compiled but never executed numerically,
and no test was added to the repository test suite.

## Complex `gesv` and `posv`

For symmetry with the new `gels`, the Armadillo `gesv` and `posv` wrappers
gained their complex flavors: `cgesv_`/`zgesv_` and `cposv_`/`zposv_`, with the
matching specializations. The `Matrix`-level overloads are already templated and
needed no change — unlike `gels` there is no workspace array, so the
real-valued-`work[0]` problem does not arise. The pivot vector stays `int` for
the complex `gesv` flavors.

The Blaze side needed no edit at all: `fn_BZ_gesv.hpp` and `fn_BZ_posv.hpp`
delegate to `blaze::gesv` / `blaze::posv`, which are already templated over
`complex<float>` and `complex<double>`. This was confirmed by running the
complex tests against the unmodified Blaze headers.

Note that the complex `posv` flavors solve a **Hermitian** positive definite
system, not a symmetric one — a matrix that is complex symmetric but not
Hermitian is not a valid input. This is noted at the declarations.

Verified numerically on both backends: real `gesv`/`posv` regression, complex
`gesv` (2x2 non-Hermitian), complex multi-rhs `gesv`, and complex `posv` on a
Hermitian positive definite system, all recovering known exact solutions. All
in-tree `gesv`/`posv`/`polyfit` consumers still compile on both backends.

Not extended: `src/linalg/lapack/fn_LAPACK_gesv.hpp` still declares only the
real flavors. It is a separate, currently unused declaration set (only the
`getrf`/`getri`/`gemm` headers in that directory have a consumer), so it was
left alone rather than grown.

## Related change

`src/linalg/lapack/` had been staged for deletion; it was restored, since
`src/math/tensor/fn_invert_symmetric.hpp` still includes three of its headers.
It coexists with the per-backend `namespace lapack` blocks without conflict —
the shared layer declares its wrappers with const-reference parameters and the
per-backend files use pointers, so they are distinct overloads, and the
`extern "C"` declarations are identical.
