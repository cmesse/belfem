# SpMatrix: Accessor De-virtualization, Base-Agnostic Matvec, Sentinel Slot

**Date:** 2026-08-18
**Purpose:** Phase A (correctness) and Phase B (accessor performance) of the SpMatrix
             optimization spec. Phase C (assembly restructure) deliberately deferred.
**Module:** `sparse`, plus `splinalg.f90`
**Status:** Implemented. `tests/sparse` builds and passes 74/74 in debug, release and under
            ASan. **The full project build and `make check` have not been run** — the user
            runs those. The prebuilt archives in `cmake-build-debug/lib` are stale relative to
            the working tree, so the test binary was relinked with freshly compiled objects.

## Why

`SolverData::assemble_jacobian` / `assemble_newton` call `SpMatrix::index()` once per element
entry, per nonlinear iteration, per timestep. Each call was an indirect jump through a
member-function pointer into one of four near-identical search functions, none of which could
inline. At the flagship scale that is order 1e9 lookups per assembly pass.

## Phase A — correctness

1. **`index()` truncated under `BELFEM_INT64`.** The chain mixed three types: the search
   functions returned `index_t` (unsigned), the public wrapper returned `int`, and
   `operator()` stored the result in `int_t`. `mNumNonZeros` is `int_t` and `set_nnz` caps
   nnz at the `int_t` maximum (`cl_SpMatrix.cpp:666-681`), so `int_t` is the correct unifier
   and the whole chain now uses it. Only three call sites in the tree consume `index()`
   (`SolverData.cpp:660`, `:1642`, `:1750`), all source-compatible.

2. **Writing an absent entry overflowed the heap in release builds.** `allocate_values()`
   allocated exactly `mNumNonZeros` reals, and the miss sentinel is `mNumNonZeros` itself —
   guarded only by `BELFEM_ASSERT`. It now allocates `nnz + 1` and zeroes the extra slot,
   because assembly accumulates (`+=`) into it and would otherwise read indeterminate memory.
   Three further allocation sites bypassed `allocate_values()` entirely (the external-array
   constructor, the HDF5 load path and the deep-copy path); all now route through it, so the
   invariant has one home. `mNumNonZeros` remains the logical length for `fill()`, HDF5 I/O,
   `std::copy` and MPI — the dump slot is never transferred.

   **Verified to be a real bug, not a theoretical one:** the pre-change code, built with ASan
   and asked to assemble into a structural zero, reports
   `heap-buffer-overflow ... WRITE`. The post-change code does not.

3. **`multiply()` rewrote the whole matrix twice per call.** It converted to Fortran base, ran
   the kernel, and converted back — two O(nnz + pointers) sweeps of the index arrays for a
   conceptually read-only operation. The Fortran kernels now take the base as an argument and
   shift internally; `multiply()` passes `indexing_base()` and converts nothing. The MKL path
   already detected the base from `aPointers[0]`, so only the flips around it were removed.
   `set_indexing_base()` is untouched and remains the solver-boundary operation (MUMPS wants
   1-based, STRUMPACK 0-based).

   Correction to the spec's premise: the old code *did* restore the base correctly (the
   restore is guarded, `if( tOldBase == 0 )`), so there was no base leak. The cost and the
   parent/child hazard were the real motivations.

4. **`matvec_csc` used one `!$omp atomic update` per nonzero**, serialising the inner loop.
   Replaced by an OpenMP array `reduction( + : y )`. The reduction gives each thread a private
   `y` (n * 8 bytes per thread); the fallback if that ever binds is the transposed-CSR duality
   path, which gathers instead of scattering.

   The base shift is `off = 1 - base` for array access and `- base` for the slice upper bound,
   rather than `- 1 + off`, which is the same value with one operation fewer.

## Phase B — accessor

5. **De-virtualized the lookup.** `mIndexFunction` and the four `index_csr/csc_{zero,one}_based`
   functions are gone, replaced by one inline `SpMatrix::position()` in the header that reads
   the base from `mPointers[0]` and branches on `mType`. Both branches are loop-invariant in
   every caller.

   This is also a **correctness** improvement the spec did not claim. `set_indexing_base()` on
   a *child* propagates up to the parent (`cl_SpMatrix.cpp:803-809`), but on a *parent* it did
   not update the child's cached `mIndexFunction`. Since `SolverData` builds the Jacobian as a
   child of System, a parent left in Fortran base would leave the child running a zero-based
   search against one-based arrays. Reading the base live removes the whole class of bug.

6. **The hybrid linear/binary search was measured and rejected.** The spec proposed scanning
   linearly below a threshold of ~32 on the theory that FE rows are short. It is slower at
   every row length measured, so `position()` uses `std::lower_bound` unconditionally:

   | entries per row | binary (ns/lookup) | linear ≤16 | linear ≤32 | linear ≤64 |
   |---:|---:|---:|---:|---:|
   |  8 | **5.22** | 6.52 | 6.88 | 7.65 |
   | 16 | **5.91** | 8.82 | 8.62 | 8.71 |
   | 24 | **6.13** | 7.00 | 11.34 | 11.12 |
   | 32 | **6.60** | 7.49 | 15.33 | 14.82 |
   | 50 | **7.10** | 8.05 | 8.13 | 21.04 |

   A whole slice sits in one or two cache lines, so `lower_bound` costs ~6 well-predicted
   steps, while the scan averages k/2 iterations behind a data-dependent branch. The threshold
   macro was removed rather than set to zero — a knob over dead code is worse than neither.

   **This result is not universal, and the table above should not be quoted as if it were.**
   Fable, re-running the same experiment independently on different hardware, measured linear
   *winning* at 8 entries per row (29 ns vs 43 ns) and losing at 24 and 50 — i.e. the textbook
   crossover does exist there. The conclusion that survives both measurements is narrower than
   the table: the crossover point is microarchitecture- and access-pattern-dependent, the
   difference is small at the row lengths BELFEM actually uses, and neither machine is the
   production cluster. That is the reason the knob was deleted rather than tuned. **If this is
   ever revisited, re-measure on the production nodes — do not resurrect a threshold from
   either of these runs.**

7. **`positions_in_slice()`** added: one merge join over a slice against a strictly ascending
   column list, O(slice + cols) streaming instead of n independent O(log k) probes. Phase C
   builds on it; nothing calls it yet. Christian confirmed (2026-08-18) that an element's dof
   list never repeats a global index, so the list is *strictly* ascending and the merge needs
   no tie handling — asserted in debug so a future change to dof construction fails loudly
   rather than silently dropping a contribution.

## Measurements

Development machine (x86_64), project flags with `-O2 -DNDEBUG`. Accessor benchmark is the
`assemble_jacobian` access shape: 20,000 synthetic elements, tN = 24, tN^2 lookups each.

| entries per row | before | after | |
|---:|---:|---:|---:|
|  8 | 7.37 ns | 5.40 ns | 1.37x |
| 16 | 7.56 ns | 5.79 ns | 1.31x |
| 24 | 8.13 ns | 6.25 ns | 1.30x |
| 32 | 8.33 ns | 6.52 ns | 1.28x |
| 50 | 8.84 ns | 6.90 ns | 1.28x |

Matvec, 500,000 rows x 30 entries, 50 products, 4 threads:

| | per matvec |
|---|---:|
| before (two base sweeps per call) | 30.96 ms |
| after (base passed to kernel) | **16.88 ms** |

**1.83x** — the two deleted sweeps were 45 % of the call.

So Phase B is a ~1.3x constant on the accessor, not the order-of-magnitude that Phase C
targets. That is worth knowing before deciding on Phase C: the assembly restructure has to
earn its numerical-reproducibility cost against a baseline that is now 1.3x faster.

## Production result (2026-08-18, Linux, in BELFEM)

| phase | before | after | |
|---|---:|---:|---:|
| assembly (first sample) | 9,867-9,880 ms | 7,740 ms | **~1.3x** |

That matches the accessor microbenchmark (1.28-1.37x) closely enough to treat the synthetic
measurement as representative of the real access pattern.

**Attribution note.** The same branch also carries the `DynamicBitset` two-level summary
(`dl20260818_dynamicbitset_summary_bitmap.md`), which moved Jacobian initialization from
290 s to 23 s. That win belongs to the bitset change, not to this one: initialization is
dominated by `populate_graph`, and the speedup scales with problem size the way a removed
O(N^2) does, not the way a constant-factor accessor change would. This change's contribution to
initialization is confined to `create_assembly_tables`, whose lookup arithmetic was measured at
0.26 s at production scale and was never a bottleneck.

The 1.83x matvec improvement does not appear in either number above — it lands in solve time,
so it is worth checking there separately.

## Tests

`tests/sparse/test_SpMatrix.cpp` §5, seven new tests (67 -> 74):

`PositionMatchesReferenceInBothBases`, `PositionHandlesLongAndShortSlices`,
`PositionsInSliceMatchesPerEntryLookup`, `MultiplyIsBaseNeutral`,
`MultiplyAlphaBetaIsBaseNeutral`, `WriteToAbsentEntryHitsTheSentinelSlot`,
`CscMultiplyMatchesReferenceAcrossThreadCounts`.

The position tests check CSR and CSC in both bases against a reference derived from the raw
arrays, on a pattern with an empty row, a full row, and first/last-column hits. The multiply
tests assert `indexing_base()` is unchanged across the call, including when the matrix is
deliberately left in Fortran base.

Results: **74/74** in debug (`-Og -g`), release (`-O2 -DNDEBUG`) and under
`-fsanitize=address`. The 67 pre-existing tests — which cover multiply, transpose, COO
indices, HDF5 round trips and parent/child structure sharing — pass unchanged.

## Not done

Phase C (batched assembly in `SolverData`) is deliberately deferred: it changes floating-point
summation order, and in a solver driven to eps < 1e-11 specifically to avoid checkerboarding
(Messe et al. 2023) that deserves its own decision made against a reference solve. Phase A and
B are numerically inert.

`assemble_full_matrices` was left on the per-entry accessor, as the spec directs — it is only
reached when JEDI powers are in use.

## Phase C: recommendation on record

Both reviewers converged on parking it, on stronger grounds than "it changes rounding":

- The accessor is now ~6.5 ns and A/J already share a single lookup, so per-element lookup cost
  is a few microseconds against element kernels (thin-shell Maxwell integration) that very
  likely dominate the assembly phase.
- **Gate it on a profile, not on argument.** Profile one representative transient step of a
  tapestack-class run. If `position()` inside `assemble_jacobian` is under ~10 % of the
  assembly phase, a 2-4x on a non-dominant term does not buy back the cost of losing
  bitwise-stable baselines.
- `positions_in_slice()` is merged and tested, so turning Phase C on later is a small change
  rather than a campaign.

## Follow-up

`CLAUDE.md` claimed BELFEM's own code carries no `omp` pragmas. That is false for the Fortran
side — `splinalg.f90`, `arpacktools.f90` and `parpacktools.f90` all carry them, and item 4
edits one. Corrected in this session.
