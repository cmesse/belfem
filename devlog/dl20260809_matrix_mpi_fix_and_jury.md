# Matrix-MPI Transfer Fix, Rule-of-Five Hardening, Test Gating — Jury-Reviewed

**Date:** 2026-08-09
**Purpose:** Execute the sign-off items from the 2026-08-08 verification pass and jury-review the batch
**Files:** `src/comm/commtools.hpp`, `src/linalg/{armadillo,blaze}/cl_{AR,BZ}_Matrix.hpp`, `src/core/cl_Logger.{hpp,cpp}`, `src/containers/cl_DynamicBitset.{hpp,cpp}`, `src/containers/cl_StringList.cpp`, `CMakeLists.txt`, `tests/fem/CMakeLists.txt`, `src/fem/thermal/CMakeLists.txt`, `tests/comm/test_CommMPI.cpp`, docs

## Fixes (all user-approved)

1. **Matrix MPI overflow (the priority):** all five transfer paths
   (`broadcast`/`send`/`receive`/`distribute`/`collect`) now transmit
   `spacing()*n_cols` — the padded footprint of the current shape — instead of
   `capacity()`, which stays stale-large after a Blaze shrink and overflowed
   the receiver's exact-fit buffer. New `spacing()` accessor on both backends;
   `BELFEM_ERROR` capacity guards on the receive sides; INT_MAX bound on the
   unchunked `broadcast`. Regression test `SendReceiveMatrixAfterShrink`
   added (grow 20×20 → shrink 3×4 → send/receive → value check).
2. **Logger:** copy/move deleted (owns `FILE*`), `fopen` checked. No caller
   copies existed (grep-verified).
3. **DynamicBitset:** copy-assign copies `mIndex`; move ctor `noexcept`;
   allocations checked; **zero-size bitsets skip allocation** (see jury).
   **StringList:** `push` bound + allocations upgraded to `BELFEM_ERROR`.
4. **USE_MAXWELL:** full static audit — OFF-build configures/compiles/links
   cleanly after the tests/fem gate (2026-08-08); kernel's `MaxwellData` is
   kernel-owned; dead `include_directories(fem/maxwell)` removed from thermal.

## Jury round (blind Codex+Grok, `tmp/ai_exchange/review_matrix_mpi_and_test_gating.md`)

- **Caught real (2/2 agreement, fixed post-round):** the new bitset allocation
  guard aborted on `DynamicBitset(0)` where `malloc(0)` may return null —
  reachable in production (`SourceExpander` sizes face bitsets 0 on 2D
  meshes); and doc drift (coding_philosophy ×2, tests_06_comm.md) still
  describing the `capacity()` protocol — synced.
- **Refuted with TPL source:** Grok's P1 "sticky spacing" corruption — Blaze
  `resize` sets `mm_ = addPadding(m)` unconditionally in every branch, so
  spacing is a pure function of the current row count.
- **Routed to Christian:** the pre-session `VERSION 0.1.0 → 0.9.0` bump rides
  in the same tree — land with this batch or split the commit.

## Pending executable gate

MPI `make check` with `-DUSE_TEST=ON` (Christian runs builds), including the
new shrink regression test. All evidence so far is source-trace level.
