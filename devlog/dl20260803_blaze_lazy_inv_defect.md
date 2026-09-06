# Blaze lazy `inv()` defect: silent transposition on `B * inv( A )`

**Date:** 2026-08-03
**Purpose:** Record the root cause of the `cl_EF_QUAD4TS.cpp:50` abort, the backend fix
in `fn_BZ_inv.hpp`, and the QUAD4TS pseudo-inverse simplification
**Module:** linalg (blaze backend), fem/interpolation

## Symptom

Running an example aborted inside `EF_QUAD4TS::link()` at `cl_EF_QUAD4TS.cpp:50`:

```
blaze/math/dense/DynamicMatrix.h:4422:
Assertion `( i<m_ ) || blaze::ASSERT_MESSAGE( "Invalid row access index" )' failed.
  blaze::smpAssign<DMatTransposer<...>>
  blaze::solve1x1<...>
  blaze::solve<...>
  DynamicMatrix::operator=<DMatTransExpr<DMatDMatSolveExpr<...>>>
```

Initially suspected to be fallout from the recent LAPACK/Blaze interface work. It is
not: Blaze is 3.9.0 with headers dated Feb 2023, `fn_BZ_inv.hpp` and `fn_trans.hpp` had
only ever been touched by the copyright commit, the new `src/linalg/lapack/*` wrappers
are not included by that translation unit, and `BLAZE_BLAS_IS_64BIT` is inert in a build
without `BELFEM_INT64`. The code path is backend-specific — it works under Armadillo,
whose `inv()` evaluates eagerly.

## Root cause

`blaze::inv( matrix )` is **lazy**; it returns a `DMatInvExpr` (`DMatInvExpr.h:405`).
Blaze then restructures

```
B * inv( A )   -->   trans( solve( trans( A ), trans( B ) ) )      // MatInvExpr.h:190
```

`solve()` writes the result through a `DMatTransposer` view of the destination and calls
`resize( *X, A.rows(), B.columns() )` (`LSE.h:254`). But

```cpp
DMatTransposer::resize( m, n ) { dm_.resize( m, n ); }             // DMatTransposer.h:353-355
```

forwards the extents **without swapping them**. For a non-square left operand the
destination is therefore resized transposed.

Measured with a probe linked against the real debug build:

| expression       | DEBUG              | NDEBUG             |
|------------------|--------------------|--------------------|
| `2x2 * inv(2x2)` | correct            | correct            |
| `3x2 * inv(2x2)` | dest becomes 2x3   | dest becomes 2x3   |
| `2x1 * inv(1x1)` | abort (assert, OOB)| abort              |

The release column is the dangerous one: no assertion, silently transposed result.

`inv( A ) * B` (inversion on the **left**) takes a different operator
(`MatInvExpr.h:161`) writing into a plain resizable destination, and is unaffected.

Affected call sites:

- `cl_EF_QUAD4TS.cpp:50` — `trans( mJ ) * inv( mGram )`, `2x1 * inv(1x1)` → abort
- `cl_FVM_Factory.cpp:783` — `mJ * inv( trans( mJ ) * mJ )`, `3x2 * inv(2x2)` → silently
  transposed, no diagnostic

`cl_EF_PENTA6TS.cpp:61,64` materialises `mInvGram` into a member before multiplying and
`cl_EF_HEX8TS.cpp:380-384` uses a hand-rolled 2x2 adjugate, so neither was ever exposed.

## Fix 1 — `src/linalg/blaze/fn_BZ_inv.hpp`

Make the backend wrapper eager so the broken restructuring can never be selected,
matching the Armadillo backend:

```cpp
    template < typename T, typename = blaze::EnableIf_t< blaze::IsMatrix_v< T > > >
    blaze::DynamicMatrix< blaze::ElementType_t< T >, BLAZE_DEFAULT_STORAGE_ORDER >
    inv( const T & aExpression )
    {
        return blaze::inv( aExpression );
    }
```

Three details are load bearing, each found by an auditor rather than by inspection:

1. **The parameter must be the raw backend type, not `belfem::Matrix< T >`.**
   `fn_inv.hpp:25-41` calls `inv( aA.matrix_data() )`, so ADL resolves to `blaze::inv`
   and a `Matrix`-typed backend overload is never a candidate — an earlier attempt with
   that signature left the abort in place. This matches the existing pattern in
   `fn_AR_inv.hpp:20-26` and `fn_BZ_det.hpp:21-27`.
2. **The `IsMatrix_v` guard is required.** `blaze::inv()` also has a scalar overload
   gated on `IsScalar_v`, which is true for unknown types, so an unconstrained template
   swallows `belfem::Matrix` and fails inside `blaze/math/shims/Invert.h:78`.
3. **The return type must pin the storage order.** `blaze::evaluate( blaze::inv( expr ) )`
   returned column-major for a plain matrix but **row-major** for a product expression,
   which contradicts the column-major invariant in `CLAUDE.md`.

Cost: the lazy `inv( A ) * B → solve( A, B )` rewrite is no longer available. The one
real (non-debug) site is `cl_Gradient.cpp:206`, inside the integration-point loop.
Measured, best-of-6, 200k reps, 3x3 Jacobian with a 27-node `dNdXi`:

| build          | eager inv | blaze solve | ratio |
|----------------|-----------|-------------|-------|
| `-Og` debug    | 219 ns/it | 157 ns/it   | 1.39x |
| `-O3 -DNDEBUG` | 147 ns/it | 124 ns/it   | 1.19x |

~23 ns per integration point in release on a postprocessing path — accepted.

## Fix 2 — `src/fem/interpolation/nedelec/cl_EF_QUAD4TS.cpp:49-51`

The Gram matrix of a LINE2 facet is 1x1, so `J+ = J^T ( J J^T )^-1` collapses to a
scalar division:

```cpp
            // using the pseudoinverse here. the gram matrix of a line is 1x1,
            // so J+ = J^T ( J J^T )^-1 collapses to a scalar division
            mPseudoInvJ = trans( mJ ) / mGram( 0, 0 );
```

Over 720 geometries (angles 0-180 deg, lengths 1e-3..1e4) the old and new forms are
**bitwise identical**. Over 500 further geometries the shape stays 2x1,
`|J*J+ - 1| <= 1e-13`, and `J+` stays parallel to `J^T` to relative 1e-15.

This is specific to the 1x1 Gram of a 1D-manifold facet. The 2D-manifold elements
(`cl_EF_PENTA6TS.cpp`, `cl_EF_HEX8TS.cpp`, `cl_FVM_Factory.cpp`) need a genuine 2x2
inverse; scalar division would be wrong there except for an isotropic Gram.

Degenerate zero-length facets now give inf/nan rather than a Blaze throw. Accepted:
`cl_EF_QUAD4TS.cpp:61-65` already divides by the square root of the same quantity
without a guard, so no deliberate diagnostic was lost.

## Verification

- 11 translation units including `fn_inv.hpp` compile clean under the real build flags
  (`cl_EF_QUAD4TS`, `cl_EF_PENTA6TS`, `cl_EF_TRI3`, `cl_EF_TRI6`, `cl_EF_TET4`,
  `cl_Cohomology`, `cl_Gradient`, `fn_Mesh_compute_surface`, `mt_maxwell_symmetry`,
  `mt_maxwell_background`, `cl_IWG_StaticHeatConduction`).
- `cl_FVM_Factory.cpp:783` pattern now yields 3x2 with correct values (was 2x3).
- Storage order is column-major for both plain matrices and product expressions.

`cl_FVM_Factory.cpp` itself does not currently compile for unrelated reasons (undeclared
`mNodes` at :894/:913, unused-variable `-Werror` at :154/:522/:524); its object file in
the build tree is stale.

## Fix 3 — the edge functions made allocation free

Eager `inv()` costs exactly **one `posix_memalign` per call** (Blaze `DynamicMatrix`
allocates through `blaze::alignedAllocate`, `Memory.h:94`). Measured with an interposed
allocator: old form 1 per call, every replacement 0. That is a per-element allocation in
`link()`, against CLAUDE.md's rule on temporaries in frequently-called methods, so every
site was converted to write into a pre-sized member:

| site | was | now |
|------|-----|-----|
| `cl_EF_PENTA6TS.cpp:60-61`  | `mInvGram = inv( mGram )` (2x2 Gram) | `inv2( mGram, mInvGram )` |
| `cl_EF_PENTA6TS.cpp:139-140`| `mInvJ = inv( tJ )` (3x3)            | `inv3( tJ, mInvJ )` |
| `cl_EF_QUAD4TS.cpp:112-119` | local `Matrix< real > tJ`            | member `mJ2` |
| `cl_EF_QUAD4TS.cpp:121-122` | `mInvJ = inv( tJ )` (2x2)            | `inv2( mJ2, mInvJ )` |
| `cl_EF_TRI3.cpp:51`         | `mInvJ = inv( mJ )` (2x2)            | `inv2( mJ, mInvJ )` |
| `cl_EF_TRI6.cpp:261`        | `mInvJ = inv( mJ )` (2x2)            | `inv2( mJ, mInvJ )` |
| `cl_EF_HEX8TS.cpp:380`      | hand-rolled 2x2 adjugate (pre-existing) | `inv2( mG, mInvG )` |

`mInvGram` had no `set_size()` in the PENTA6TS constructor — it was sized implicitly by
the old assignment — so `mInvGram.set_size( 2, 2 )` was added.

Routing the thin-shell elements through `inv2`/`inv3` was only possible after Fix 4;
before that their singularity test rejected these matrices outright. The first version of
this fix therefore used inline adjugates, and Fix 4 replaced them with the shared helpers
— including the one in `cl_EF_HEX8TS.cpp` that predates this work.

## Fix 4 — `inv2`/`inv3` test rank instead of scale

Both helpers asserted `std::abs( det ) > BELFEM_EPS`. That is an **absolute** bound on a
dimensional quantity, so it is a scale test wearing a rank test's error message, and it
misclassifies in both directions:

| case | \|det\| | Hadamard ratio | cond∞ | old test |
|---|---|---|---|---|
| equilateral TRI, h=1e-4 | 8.66e-09 | 0.87 | 2.4 | pass |
| equilateral TRI, h=1e-6 | 8.66e-13 | 0.87 | 2.4 | pass, 4000x margin left |
| PENTA6TS Gram, h=1e-4 | 7.50e-17 | 0.60 | 3.0 | **ASSERT** — false alarm |
| thin-shell tJ, L=1e-4, t=1e-9 | 5.00e-14 | 1.00 | 2e5 | pass |
| sliver TRI, 1e-9 rad, h=1 | 1.00e-09 | 1e-09 | 2e9 | **pass** — nearly collapsed |

The same triangle flips verdict purely by moving the mesh scale. Replaced by the
**Hadamard ratio** `|det| / prod( ||row_i|| )`, which lies in [0,1] by Hadamard's
inequality, equals 1 for orthogonal rows, and goes to 0 exactly as the rows become
dependent — dimensionless, scale invariant, and tolerant of the deliberate anisotropy in
thin-shell Jacobians (a condition-number test would penalise a tangent/thickness row pair
for a shape that is entirely intentional). Compared in squared form so no square roots
are needed, and the whole expression lives inside `BELFEM_ASSERT`, so it costs nothing in
release.

`BELFEM_EPSILON` (10x machine epsilon, `typedefs.hpp:90`) rather than `BELFEM_EPS`:
across the tree `BELFEM_EPS` is used as a seed or floor **value** (`Vector< real > tI( n,
BELFEM_EPS )` to keep a matrix non-singular, `mValue = BELFEM_EPS`, `aResidual =
BELFEM_EPS`), while `BELFEM_EPSILON` is the **decision threshold** constant — 187 uses
against 10, and these two asserts were the only places in BELFEM where `BELFEM_EPS` was
used as a comparison tolerance. The deliberate 10x margin is also the right side to err
on for a check whose job is catching collapse. Note the constant appears squared in the
source because `det^2` is compared against `prod( ||row_i||^2 )`; it is not a doubled
margin.

The tolerance is deliberately loose. The check exists to catch a collapsed or
uninitialised matrix — a data or logic bug — not to grade mesh quality; real meshes
contain slivers, and element-quality policing belongs in `MeshChecker`, which already
errors on negative-volume elements. One behaviour change: a matrix at Hadamard ratio
~1e-11 (aspect 1e11, still invertible to ~5 digits) used to assert on scale grounds and
now passes.

Verified against the real criterion in a probe: well-shaped triangle Jacobians and
`h^4`-scaling Gram matrices pass at every scale from h=1 down to h=1e-8; anisotropic
thin-shell Jacobians pass down to a thickness of 1e-12 m against a 1e-4 m tangent
(aspect 1e8); collapsed 2x2, row-dependent 3x3, and a Hadamard ratio of 1e-17 all assert.
`inv2` and `inv3` remain **bitwise identical** to `blaze::inv()` over 20000 random draws.

## Follow-ups (not applied)

- `cl_EF_QUAD4TS.cpp:74-78` writes `mNablaXi`/`mNablaEta`, then `:102-109` overwrites
  both; the first block is dead.
- `src/linalg/fn_inv.hpp:28` has a stray top-level `const` on a by-value return.
- `cl_EF_PENTA6TS.cpp` computes `mDetGram = det( mGram )` right after `inv2`, which
  already returns that determinant. Free to dedupe, left alone to keep this diff focused.

## Collaboration

Codex and Grok both audited. Grok refuted the first version of Fix 1 (wrong parameter
type, does not close `cl_FVM_Factory.cpp:783`) — confirmed independently by runtime
probe. Codex refuted the claim that the lost `solve()` rewrite only affected debug code
(`cl_Gradient.cpp:206`) and caught the row-major return type. Both accepted Fix 2 with
high confidence and both independently reached the same scope conclusion for the
2D-manifold elements.
