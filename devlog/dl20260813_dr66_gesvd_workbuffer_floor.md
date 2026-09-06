# DR-66: gesvd Work-Buffer Floor — Reuse Path Made Reachable Under MKL

**Date:** 2026-08-13
**Purpose:** Rule and fix DR-66 — the caller-provided work-buffer reuse path in
`fn_gesvd.hpp` was unreachable for a buffer sized by the wrapper's own query under MKL.
**Module:** `src/linalg/lapack`

## The ruling

The register posed two options: trust the vendor query as the acceptance floor, or keep
the reference floor and document that query-sized buffers are not reusable. Neither
survives contact:

- **Trusting the query is incoherent as an acceptance test.** The length test at the top
  of `gesvd` must decide *without calling LAPACK* whether the caller's buffer is usable —
  that is the entire point of the reuse path. A query result only exists after a LAPACK
  call, so the acceptance floor can only ever be a formula, and the reference-LAPACK
  minimum is the correct one: any `lwork` at or above it is legal on every conforming
  backend.
- **Documentation alone** would keep the code honest while leaving the overload's reason
  to exist permanently dead under MKL.

The actual defect is a broken invariant, not a wrong floor: the wrapper sized its own
buffer from the MKL query (7 for the 4×3 `'A','A'` double case), which lands *below* its
own acceptance threshold (15) — so a buffer the wrapper itself produced was rejected by
the wrapper on every subsequent call, re-entering the query branch with a shrink-then-grow
`set_size` pair each time, and a shrinking shape sequence (4×3 → 3×2) silently shrank a
shared buffer instead of reusing it.

Christian explicitly allowed backend-conditional handling (LAPACK/OpenBLAS vs MKL) if it
helped; it was considered and deliberately not used — one backend-agnostic path covers
all cases, with no `#ifdef BELFEM_MKL` fork to maintain. The ruling, which is the pattern
for every future wrapper in the unified LAPACK interface:

1. **The reference-LAPACK formula is the acceptance floor** (computable without a call).
2. **The wrapper's own sizing takes `max( query, formula )`** — restoring the invariant
   that a wrapper-sized buffer passes the wrapper's own acceptance test. Oversized
   `lwork` is always legal and, per the gesvd docs, generally faster.
3. **Buffers only grow** — the pre-query resize is now guarded, never shrinking a buffer
   that already covers the query call.

## Source changes (`fn_gesvd.hpp`)

Two edits inside the query branch:

```cpp
// grow-only guard before the workspace query
if ( static_cast< int_t >( Work.length() ) < tRealsPerT + tRWorkSize )
{
    Work.set_size( tRealsPerT + tRWorkSize );
}
...
// floor the vendor optimum at the reference minimum
lwork = std::max( lwork, lapack::work_size( Work( 0 ) ) );
```

On reference LAPACK/OpenBLAS, where the query returns ≥ the formula, behaviour is
bit-identical to before. The reuse branch itself is untouched and still never resizes.

## Test changes (`tests/linalg/test_lapack.cpp`)

`GesvdWorkBufferReuse` previously could not assert buffer lengths because they were
backend-dependent; with the fix they are invariant, so it now asserts (a) the first call
grows the buffer to at least `gesvd_min_work( 4, 3 )`, (b) a same-shape second call and
(c) a smaller-shape (3×2) third call both leave the length untouched. The
`GesvdWorkBufferReuseAtReferenceMinimum` header note no longer describes a live
inefficiency; it pins the exact-minimum acceptance edge, which the query branch cannot
produce when the vendor optimum exceeds the formula.

## Evidence

- **Compile:** `test_lapack.cpp` (which includes `fn_gesvd.hpp`) compiles clean with the
  tree's own flags — Blaze + MKL, `-Wall -Werror -pedantic-errors -std=gnu++17`.
- **Probe (executable, this machine, MKL + Blaze):** a standalone scratchpad probe
  linking the prebuilt `libbelfem_{core,containers,comm}.a`:
  - re-measured the raw workspace query: **optimal lwork = 7 vs reference minimum 15**
    for 4×3 `'A','A'` double — independently confirming the number the 2026-08-11 sweep
    measured;
  - with the fix, one buffer through 4×3 → 4×3 → 3×2 holds length **15 / 15 / 15**,
    `info = 0` throughout, max reconstruction error `‖A − U·S·VT‖∞ ≤ 3.4e-15`.
- **Not run:** the gtest suite itself (`USE_TEST=OFF` in the shared tree). The gate is
  `make check-fast` with `USE_TEST=ON` — the same never-yet-run gate as DR-42/49. The
  probe covers the double flavor only; the fix sits in the flavor-generic template path.

Register: DR-66 → fixed (uncommitted), probe-verified, suite gate outstanding.
