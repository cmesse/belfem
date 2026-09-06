# MUMPS Symmetric Mode: Lower-Triangle Extraction

**Date:** 2026-08-29
**Purpose:** Design — ALGORITHM ONLY — for supplying MUMPS one triangle instead of the full matrix
when `SYM != 0`, so that `SymmetryMode::GeneralSymmetric` / `PositiveDefiniteSymmetric` become
supported rather than refused. Halves factor time and memory on a symmetric system.
**Module:** `src/sparse`
**AIs involved:** Claude (draft), Codex (audit)
**Status:** **DRAFT, AUDITED 2026-08-29 — corrections folded in below. NOT FOR IMPLEMENTATION
TODAY** ( Christian, 2026-08-29:
"worth thinking about implementing tomorrow, but not today" ). No source touched.

---

## 1. Why this exists

`MUMPS::initialize()` currently rejects any symmetric mode with an always-active error. That is
the honest state of the code, not a design: BELFEM stores the FULL matrix, MUMPS with `SYM != 0` wants
exactly ONE triangular representative of each symmetric coordinate, and nothing extracts one — so
a symmetric mode made MUMPS treat `( i, j )` and `( j, i )` as DUPLICATES and SUM them, doubling
every off-diagonal and factorizing a different matrix **without failing**. Measured
2026-08-29 on a matrix with an analytic spectrum: first-solve residual `||Ax-b||/||b|| = 7.4e16`
and a converged eigenvalue of `-2.5e-19` against a true `2.46e-6`, with the eigensolver reporting
success throughout. See `todo/conditioning_shift_invert_fallback.md` §3.5.

This document is what it would take to say yes instead of no.

## 2. Where it hooks — the COO view, master only

MUMPS is fed centrally: the master passes triplets and every worker passes nulls
( `cl_SolverMUMPS.cpp:487-499` ):

```
mumpstools_solve( ..., aMatrix.n_rows(), aMatrix.number_of_nonzeros(),
                  1, aMatrix.rows(), aMatrix.cols(), aMatrix.data(), ... )
```

so the extraction is a MASTER-ONLY transformation of that triplet view, and no MPI enters the
algorithm at all. `create_coo_indices()` has already materialised whichever index array the
storage order lacks ( `rows()` for CSR, `cols()` for CSC ), so **CSR and CSC need no separate
code path** — both present the same `( rows, cols, values )` triple by the time MUMPS is called.
That is the single most useful simplification available here.

> **CORRECTED 2026-08-29 ( Codex audit §3, against the local MUMPS 5.5.1 user guide §5.2.2.1
> p28 ).** This document originally asserted that MUMPS *requires* the lower triangle. **It does
> not.** MUMPS accepts EITHER half, including the diagonal; an entry from the opposite triangle is
> not ignored or rejected, it simply represents that symmetric coordinate. What is fatal is
> supplying BOTH, because duplicate triplets are summed — in symmetric and unsymmetric assembled
> mode alike. So lower extraction is a **BELFEM policy choice**, not a MUMPS requirement, and the
> design should say so. Only indices outside `1..N` are silently ignored.

**Square only.** MUMPS is given `n_rows()` as `N` and never sees `n_cols()`
( `cl_SolverMUMPS.cpp:487` ), so symmetric mode must always-active reject `n_rows() != n_cols()`
before extraction rather than trusting the caller.

**Base invariance, verified 2026-08-29 ( and re-confirmed by the audit, including for CHILD
matrices — a child delegates the base switch to its parent and aliases the shared arrays, so it
cannot hold an independently based copy ):** `set_indexing_base()` shifts `mPointers`, `mRows` and
`mColumns` together ( `cl_SpMatrix.cpp:829-848` ), so `rows[k]` and `cols[k]` are always in the
SAME base and the predicate `rows[k] >= cols[k]` is base-invariant. The algorithm therefore does
not care whether it runs before or after the Fortran switch.

## 3. The algorithm

Split by lifetime, because the pattern is stable across solves while the values are not.

**Build the gather map — once per PATTERN:**

```
m = 0
for k in 0 .. nnz-1 :
    if rows[ k ] >= cols[ k ] :          # lower triangle, diagonal included
        mMap( m )  = k                   # where to fetch the value from
        mIrn( m )  = rows[ k ]
        mJcn( m )  = cols[ k ]
        m = m + 1
mNnzLower = m
```

**Gather the values — once per SOLVE:**

```
for m in 0 .. mNnzLower-1 :
    mVal( m ) = data[ mMap( m ) ]
```

then hand MUMPS `( n, mNnzLower, mIrn, mJcn, mVal )`.

Three properties worth stating explicitly:

- **No allocation per solve.** `mMap`, `mIrn`, `mJcn`, `mVal` are members sized once when the map
  is built, matching the framework rule against temporaries in repeatedly-called methods.
- **One pass, no sort, no search.** Emission order follows the source order, which for CSR is
  row-major and for CSC column-major. MUMPS accepts centralized assembled input in any order.
- **Duplicates are neither created nor removed.** Every kept entry maps to exactly one source
  entry. If the SOURCE matrix contains duplicate `( i, j )` pairs, MUMPS sums them — the same
  behaviour it has today, unchanged by extraction.

**Invalidation — THE DRAFT'S ORIGINAL ANSWER WAS UNSOUND ( audit §7.1 ).** I proposed reusing the
`select_job()` trigger: matrix pointer, dimensions, `nnz`. That is not a pattern fingerprint.
`select_job()` treats the same object as the same pattern by POINTER IDENTITY alone
( `cl_SolverMUMPS.cpp:429` ), and the frozen record adds shape and storage addresses but never
index CONTENTS ( `:405` ). An in-place edit of `mRows` / `mColumns` / `mPointers` that keeps the
allocations and counts intact evades every one of those triggers, and the map would then gather
values from the wrong offsets — silently.

So this needs one of: a structural GENERATION COUNTER owned by `SpMatrix` and bumped by anything
that mutates the pattern ( the honest fix, but a broad `SpMatrix` API change — note the class hands
out mutable `data()` widely, so the same objection that killed a generation counter for the frozen
record applies here ); an explicit contract that the caller must invalidate; or an actual hash of
the index arrays. **Unresolved — this is now the first thing tomorrow has to decide**, because
everything else in the design is contingent on the map being trustworthy.

## 4. The hard part is NOT the extraction

The loop above is ten lines and cannot really go wrong. **The risk is that extraction silently
SYMMETRIZES a matrix that is not symmetric.** Discarding the upper triangle and telling MUMPS
`SYM != 0` asserts `A = A^T`; if that is false, the solve is of `tril(A) + tril(A)^T` and there is
no error anywhere. That is precisely the failure class this whole line of work exists to remove,
so shipping the extraction without an answer here would trade one silent-wrong for another.

Options, in increasing cost:

| option | cost | catches |
|---|---|---|
| document the caller's promise, `BELFEM_ASSERT` only | free in release | nothing in release |
| **structural** check once per pattern: is the sparsity pattern symmetric? | O( nnz ) per pattern | wrong-matrix and wrong-storage errors, not wrong values |
| **numeric** check per solve: `A(i,j) == A(j,i)` within a tolerance | O( nnz ) per solve, plus a lookup per entry | everything |
| numeric check on the FIRST solve of a pattern only | O( nnz ) once | everything except values that lose symmetry later |

> **SUPERSEDED 2026-08-29 by a better design from the audit ( §4, §5 ).** My recommendation was
> "structural once per pattern plus numeric on the FIRST solve only". The audit refuted the
> adequacy of that and it was right: a fixed symmetric pattern can acquire asymmetric VALUES at any
> later nonlinear iterate, so a first-solve check proves one assembly and every later extraction
> silently discards the evidence — recreating the exact failure class this work exists to remove.
> It also made my O2 moot: `SpMatrix::operator()` is not the right tool for the transpose lookup.

**The design to implement instead — validation rides the gather, per solve, for almost nothing:**

```
build, once per pattern:   for each kept lower entry m, also record
                             mPartner( m ) = index of the ( j, i ) entry
                           a MISSING partner is the structural check, detected here

gather, every solve:       mVal( m ) = data[ mMap( m ) ]
                           if | data[ mMap( m ) ] - data[ mPartner( m ) ] | > tol : REFUSE
```

One extra indexed load and one comparison per entry, no search, no allocation, and the structural
check falls out of building the partner map rather than costing a separate pass. That is cheaper
than my per-pattern structural scan AND strictly stronger than my first-solve numeric check,
because it validates the values actually being handed to MUMPS on the solve they are handed over.

Two details the implementer must still settle: the comparison tolerance ( exact equality is wrong
for a matrix assembled by summation, where `A(i,j)` and `A(j,i)` may differ in the last bits ), and
what "refuse" means — an error, or a fall back to `SYM = 0` with a warning. Given §1's history, an
error is the safer default.

**Duplicates break naive pairwise checking ( audit ).** If the source contains duplicate `( i, j )`
triplets, "the" partner of an entry is not well defined and a pairwise comparison is not the right
test — the comparison must be between the SUMS of each coordinate's duplicates. Whether BELFEM's
assembled matrices can contain duplicates at all is a precondition to establish, not assume.

## 5. Cost model — is it even worth it?

Pays for itself only when the factorization dominates:

- **Saved:** MUMPS `LDL^T` against one triangle. INPUT storage halves exactly; factor fill and
  runtime depend on ordering, pivoting, scaling and structure, so the benefit is
  workload-dependent and often large but **not** an algorithmic guarantee of 2x
  ( corrected, audit §7.8 ).
- **Paid ( corrected, audit §7.2 ):** THREE `int_t` arrays — `map`, `irn`, `jcn` — plus one `real`
  value array, i.e. **20 bytes per kept entry** at 32-bit `int_t`, not the 12 this draft first
  claimed; plus another 4 for the partner map of §4 unless it shares storage. One O( nnz ) pass per
  pattern, and one O( nnz/2 ) gather plus comparison per solve.

For the magnetic Jacobian, where the STRUMPACK factor was measured at 2.64e9 nonzeros and 21 GiB,
halving the factor is worth far more than a gather over 2e6 source entries. For a small system it
is a loss. So this should be **opt-in per solver instance and never automatic**, which is also
what keeps it out of every existing production path.

## 6. Preconditions the implementer must handle ( all from the audit, §7 )

- [ ] **Both solve overloads.** The vector path ( `cl_SolverMUMPS.cpp:441` ) and the multiple-RHS
      path ( `:593` ) call MUMPS independently. Touching only the first leaves symmetric multi-RHS
      wrong — the same trap `select_job()` had.
- [ ] **Map construction must follow `create_coo_indices()`**, since both coordinate arrays have to
      exist. Today that call sits immediately before the base switch ( `:465` ).
- [ ] **Derived-buffer lifetime.** MUMPS may retain `IRN`/`JCN`/`A` pointers between JOB calls, so
      the derived arrays must not be reallocated while a factorization is live. Guarding only the
      SOURCE arrays does not prove that — this is the correction to O4, and it interacts with the
      frozen-factorization scope: under JOB 3 no value gather is needed at all.
- [ ] **Degenerate shapes.** `create_coo_indices()` allocates a dummy slot at `nnz == 0`
      ( `cl_SpMatrix.cpp:947` ). Define `N = 0`, `nnz = 0` and `nnz_lower = 0` explicitly rather
      than relying on zero-length Fortran dummy-array behaviour, which is not portable.
- [ ] **Overflow.** The kept count and map offsets are `int_t`, which may be 32-bit. Use checked
      conversions and confirm the lower count is representable in the MUMPS integer ABI.
- [ ] **Square check**, always-active, before extraction — see §2.
- [ ] **Assert both COO pointers are non-null** after creation, so a future `SpMatrix` change
      cannot quietly break the one-path assumption.

## 7. Open questions

> **O2 and O3 are RESOLVED by the audit and struck.** O2 ( is `operator()` affordable for the
> transpose lookup ) is moot: the partner-map design in §4 needs no lookup. O3 ( should a missing
> diagonal be refused ) is answered — for `SYM = 2` a structurally absent diagonal is LEGAL input
> and whether it causes trouble is a numerical property of the whole matrix, not something the
> extraction layer should adjudicate; a blanket rejection would be wrong. For `SYM = 1` it is a
> useful diagnostic and a defensible conservative policy, but stronger than MUMPS's own contract.

- **O1 — `SYM = 1` versus `SYM = 2`. Narrowed by the audit to two defensible answers, still
  Christian's call.** Either map the modes honestly ( `GeneralSymmetric` → `SYM = 2`,
  `PositiveDefiniteSymmetric` → `SYM = 1` with the definiteness precondition documented as the
  caller's, and decode a null/negative-pivot failure into "retry with `SYM = 2`" ), or support only
  `SYM = 2` initially. What the audit rules OUT is silently mapping `PositiveDefiniteSymmetric` to
  `SYM = 2`: that would change pivoting and performance behind the caller's back. Note a successful
  `SYM = 1` factorization does NOT prove positive definiteness, so it cannot be used as a check.
- ~~**O2 — how to test value symmetry affordably.**~~ RESOLVED, see above.
- ~~**O3 — the missing diagonal.**~~ RESOLVED, see above.
- **O6 ( NEW, and now the blocking one ) — how is the extraction map invalidated?** The draft's
  original answer was unsound ( §3 ). A generation counter on `SpMatrix` is the honest fix but is a
  broad API change against a class that hands out mutable `data()` freely — the same objection that
  killed a generation counter for the frozen-factorization record. Nothing else in this design can
  be trusted until this is settled.
- **O5 — is the whole thing worth doing?** The measured driver for symmetric mode was halving a
  large factorization. If the only symmetric system in play is the thermal one — 97095 dofs, small
  next to the magnetic problem — the saving may not repay the risk of a silent symmetrization bug.
  A decision, not a calculation, and Christian's to make.

## 8. Not in this document

No step list, no member names, no file diffs. This is the algorithm and its preconditions only, so
that tomorrow's implementation starts from something audited rather than from a blank file.

**The auditor's overall verdict, kept verbatim because it is the right summary:** retain the
unified extraction design, replace the invalidation mechanism, build a transpose/coordinate-group
map during structural setup, and validate numeric symmetry during every non-frozen value gather.
Until those preconditions are specified, **the current always-active rejection of symmetric modes
remains correct** — i.e. nothing about today's code needs to change for this to sit here safely.
