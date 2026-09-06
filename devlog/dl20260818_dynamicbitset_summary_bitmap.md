# DynamicBitset: Two-Level Summary Bitmap for O(k) where() and reset()

**Date:** 2026-08-18
**Purpose:** Remove the O(N) full-array scan from `DynamicBitset::where()` and `reset()`, which
             made sparsity-pattern construction quadratic in the DOF count
**Module:** `containers`
**Status:** Implemented. The container's own gtest suite builds and passes (176/176 debug,
            159/159 release, 176/176 under ASan+UBSan), compiled standalone against the SCLS
            gtest at `/opt/scls/lib`. **The full project build and `make check` have not been
            run** — the user runs those.

## Why

`devlog/dl20260818_jacobian_init_quadratic.md` traced a ~200 s "Initialize Jacobian" phase
(2.1e6 DOFs, 5 ranks) to `SolverData::populate_graph`: per row it resets an N-bit workspace and
extracts the marked bits, and both operations scanned all `N/64` words regardless of how few bits
were set. With N rows that is O(N^2).

The design came from an independent sweep by Fable (v1, then v2 after call-site analysis); this
entry records what was implemented and what changed against the spec.

## What changed

`src/containers/cl_DynamicBitset.{hpp,cpp}`

Two summary levels: `mSummary1` holds one bit per data word, `mSummary2` one bit per level-1
word. For N = 2.1e6 that is 32,813 data words, 513 level-1 words, 9 level-2 words — ~1.6 %
memory overhead. `where()` and `reset()` walk the summaries and never touch a zero word.

Two invariants, both documented at the member declarations:

- **Correctness (one-directional):** every nonzero data word has its level-1 bit set, every
  nonzero level-1 word has its level-2 bit set. A summary bit standing over a zero word is
  skipped harmlessly.
- **Tightness (maintained by every mutator, asserted by the tests):** a summary bit is set iff
  the word below it is nonzero. Not needed for correctness, but losing it silently degrades
  `where()` back towards a full scan.

| Method | Change |
|---|---|
| `set( aPos )` | three unconditional ORs; still branch-free |
| `reset( aPos )` | clears the summary bits when the word empties |
| `flip( aPos )` | maintains the summaries in **both** directions — a flip can turn a bit *on* in a word whose summary bit is clear, which `where()` would otherwise miss. This is correctness, not tidiness |
| `reset()` | two-level walk instead of `memset`; moved from the header to the `.cpp` |
| `where_sparse()` | two-level walk; sorted output is structural (ascending L2 → L1 → word → bit) |
| `where_dense()` | delegates to `where_sparse` — with no zero words visited there is nothing left for a count-then-fill variant to win |
| `flip()` | early return on an empty bitset, then `recompute_summaries()` |
| `operator\|`, `\|=` | OR the summaries; OR of two tight summaries is tight |
| `operator^`, `^=`, `&`, `&=` | `recompute_summaries()` — XOR and AND can zero a word, so the result cannot be derived from the operands |
| ctors / dtor / copy+move assignment | allocate, copy, steal and free the summaries alongside `mData` |
| `data()` (mutable) | **removed** — a write through it bypasses the summaries, and after `operator\|` was rewritten to use `aResult.mData` it had no callers left. The const overload stays public |
| `summaries_are_tight()` | new public diagnostic, exact equality in both directions plus the no-bit-past-the-end safety rule |
| `recompute_summaries()` | new private helper, rebuilds both levels in one pass |

`aAssumeSparse` is retained on `where()` for source compatibility and no longer selects an
algorithm; the sole `false` caller is `fn_Graph_spfa.cpp:196`.

### Bug fixes bundled in

1. **`flip()` underflowed on an empty bitset** — `index_t` is unsigned, so the bound
   `mMemorySize - 1` wrapped to `SIZE_MAX` and the loop marched through a null `mData`.
2. **The portable ctz fallback could not compile** — `where_sparse` redeclared `int p = 0;` in a
   scope that already declared `int p;`. Unreachable on GCC/Clang/MSVC, but wrong.
3. **The portable ctz fallback computed wrong indices** in both `where_` variants — it
   right-shifted the block while searching, after which `tBlock &= tBlock - 1` operated on the
   shifted value. Both fallbacks are replaced by one `ctz64()` helper with a hard `#error` when
   no intrinsic exists.
4. **`where()`'s per-bit `if ( tIndex < mNumberOfBits )`** is now a `BELFEM_ASSERT`. Verified
   that tail bits are always zero: the constructor zeroes, `set()`/`flip( aPos )` assert bounds,
   `flip()` masks the last block, `set_from_hex()` guards the position.
5. **`lock()`'s hash sentinel** — a computed hash of exactly 0 would make a locked bitset report
   itself writable; it is now mapped to 1.
6. **`to_int()` dereferenced a null pointer on an empty bitset** (pre-existing, found by Fable in
   the conformance review). `select_to_int_function()` picked `to_int_partial` for any size below
   64 bits, including zero, and that reads `mData[ 0 ]` — null when nothing was allocated.
   Dispatch now selects a new `to_int_zero()` when `mMemorySize == 0`, which keeps the check off
   the hot path in the same way `to_int_fail` already did. Covered by
   `ZeroSizeBitsetSurvivesEveryOperation`.

## Deviations from the v2 spec

**Item 11 was wrong about the blast radius, and so was the analysis that fed it.** The spec (and
the note that informed it) stated that the only mutable `data()` write in the tree is inside
`operator|`. That was based on a search of `src/` only. `tests/containers/test_DynamicBitset.cpp`
had **seven** mutable-overload call sites, one of them a genuine external write:

```cpp
// CountIgnoresTrailingBits, before
tBs.data()[ 1 ] = ~uint64_t( 0 );   // corrupt trailing bits — count() must mask them
```

Resolved by:
- a file-local `data_ptr()` helper in the tests that reads through the **const** overload, which
  preserves the pointer-identity intent of the move tests;
- rewriting `CountIgnoresTrailingBits` to assert the invariant directly — no public mutator can
  leave a bit set at or beyond `size()` — instead of manufacturing a state that is now
  unreachable. That invariant is what licenses fix 4 above, so the two must stay together.

**T4 was strengthened.** The spec's tripwire asserted tightness and `count() == 0`. Tightness
catches a *stale* summary bit, but a mutator that wrongly *clears* one — the dangerous direction,
since `where()` then silently drops entries — is only caught if a reference comparison happens to
hit that word. `PerBitClearingKeepsSummariesTight` now compares `where()` against the naive
reference every round as well.

**The move tests gained summary coverage.** They asserted `mData` pointer identity only, which
would still pass if the summaries were left pointing at the moved-from source.

## Deliberate non-changes

`alloc_words()` zeroes what it allocates, so the copy-assignment path memsets three arrays and
then immediately `memcpy`s over them. Kept: assignment is cold, and allocate-zeroed is the safer
default for a class whose whole correctness argument rests on the summaries never containing
garbage. Flagged by Fable as a wasted memset, which it is.

## Measurements

Standalone benchmark of the real class in the exact `populate_graph` sequence (full `reset()` +
30 `set()` + one `reset( aPos )` + `where()` into a reused `Cell`), `-O2 -DNDEBUG`, on the
development machine — **not** the i9-10900X that produced the 200 s:

| N | before | after |
|---:|---:|---:|
| 2,100,000 | 59.3 s | **0.394 s** |
| 10,000,000 | (quadratic) | **2.128 s** |

**150x** at the production size, and 4.76x the DOFs costs 5.4x the time — near-linear, so the
residual quadratic is gone for any realistic N. Four passes on rank 0: ~237 s -> ~1.6 s.

Prototype figures for the record: a one-level summary reached 5.1 s at N = 2.1e6, and a
one-level design *without* the `reset( aPos )` guard measured 63.0 s — slower than the code it
replaced, because the summaries saturate under per-bit clearing.

## Production result (2026-08-18, Linux, in BELFEM)

First real in-solver measurement, old branch vs new branch:

| phase | before | after | |
|---|---:|---:|---:|
| Jacobian init, magnetic | 290,171 ms (4 min 50 s) | 23,024 ms | **12.6x** |
| Jacobian init, thermal | 30,389 ms | 3,356 ms | 9.1x |

The branch under test carries both this change and the SpMatrix accessor work
(`dl20260818_spmatrix_accessor.md`), so the two are not separated by measurement. The
initialization win is attributed here on the evidence that **the speedup scales with problem
size**: the magnetic system is roughly 3x the thermal one, its *before* time is ~9.5x the
thermal one (consistent with an O(N^2) term dominating), and it gains 12.6x against the smaller
system's 9.1x. A constant-factor accessor improvement would have produced roughly equal
speedups; a removed quadratic produces exactly this pattern. The SpMatrix accessor's own
contribution shows up in the assembly phase instead, at ~1.3x, matching its microbenchmark.

The 290 s starting point is also close to the ~237 s this devlog predicted for the four
full-width `populate_graph` passes, with the balance in the serial graph union and the
per-DOF allocations that this change does not touch.

## Verification

- Standalone harness (`ALL CHECKS PASSED`): sizes `{0, 1, 63, 64, 65, 127, 128, 1024, 4095,
  4096, 4097, 262144, 262145, 1000000}`, randomized sequences over every mutator, all six
  bitwise operators, copy/move construction and assignment, hex round-trip, full-reset cycling,
  and per-bit-clear cycling — each compared against a naive per-bit reference, with
  `summaries_are_tight()` asserted after every mutation. Asserts active.
- Clean under `-fsanitize=address,undefined`. (LeakSanitizer is not active on macOS, so this is
  not a leak check — the class remains a Valgrind candidate.)
- Compiles warning-free in the touched files under `-std=c++20 -Wall -Wextra -Wpedantic`, in
  both `-O1 -g` and `-O2 -DNDEBUG`.
- **`tests/containers/test_DynamicBitset.cpp` built and run for real**, linked against the SCLS
  gtest (`/opt/scls/lib`, 1.17.0) and only `cl_DynamicBitset.cpp` + `assert.cpp` +
  `stringtools.cpp`:
  - debug (`-O1 -g`, asserts live): **176/176 passed** — 104 parameterized + 55 `DynamicBitset`
    + 17 `DynamicBitsetDebug`, matching the count recorded in `tests_01_containers.md`;
  - release (`-O2 -DNDEBUG`): **159/159 passed** — the 17 `DynamicBitsetDebug` cases compile out
    under `NDEBUG` as designed. This is the configuration where the removed per-bit bounds check
    actually matters;
  - `-fsanitize=address,undefined`: **176/176 passed**, no sanitizer diagnostics. (The `ERROR`
    boxes in that log are BELFEM's own assertion output from the `Debug` cases, which trip
    asserts on purpose via `EXPECT_THROW`.)

**Not done:** the full project build, `make check` as a target, and any run of the real solver.
The container suite above was compiled standalone in a scratch directory and never touched the
shared build tree. The 200 s figure has not been re-measured in BELFEM — the performance numbers
are microbenchmarks of the container, on a different machine from the one that reported it.

## Tests added

`tests/containers/test_DynamicBitset.cpp` §2.12, seven tests (169 -> 176 expanded cases;
`tests/doc/tests_01_containers.md` updated):

`SummariesTightAcrossAllMutators`, `SummariesTightAcrossBitwiseOperators`,
`SummariesSurviveCopyAndAssignment`, `FullResetCyclingMatchesReference`,
`PerBitClearingKeepsSummariesTight`, `FlipInterleavedWithSetAndReset`,
`ZeroSizeBitsetSurvivesEveryOperation`.

Note recorded in the suite: every size in the existing `Sizes/DynamicBitset` parameter set is
below 4096, so none of them exercise even a second level-1 word. The new fixed-size tests at
4095/4096/4097 and 262144/262145 are what cover the walk.

## Follow-ups

- `SolverData::populate_graph` still hashes `mDofData->dof( id )` once per adjacency entry
  (measured 5.6 s per pass at 2.1e6 DOFs). Rewriting the neighbor IDs in `aGraphData` to
  `my_index` once would remove it.
- `Kernel::partition_mesh` (`cl_FEM_Kernel.cpp:255-301`) uses the same per-element
  `reset()` + `where()` shape over the element graph and benefits from this change
  automatically; worth measuring the "Partitionig mesh" timer before and after.
- `graph::symrcm`'s disconnected-component fallback (`fn_Graph_symrcm.cpp:139-158`) is O(C*V)
  and untouched by this work.
