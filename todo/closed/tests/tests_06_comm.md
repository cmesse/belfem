# BELFEM Communication Tests — Detailed Plan

**Date:** 2026-03-22
**Purpose:** Method-level test matrix for the comm module (`src/comm/`)
**Depends on:** `tests_0_strategy.md` (conventions), `tests_1_containers.md` (Cell), `tests_2_linalg.md` (Vector/Matrix)
**Confidence:** High on utility functions. Medium on MPI protocol tests (require multi-process execution).

---

## Module Overview

The comm module is BELFEM's MPI communication layer. It provides:

| Component | Content |
|---|---|
| `Communicator` class | Singleton-like global (`gComm`), MPI lifecycle, rank/size, argv, RNG, `CommunicationObject` registry |
| `commtools.hpp/.cpp` | Free functions: `send`, `receive`, `broadcast`, `distribute`, `collect`, `share` for scalars, arrays, Cell, Vector, Matrix, strings. Utilities: `comm_tag`, `comm_split`, `comm_splitcount`, `comm_barrier` |
| `commtypes.hpp` | `comm_type<T>()` → MPI datatype mapping |

**Build modes:** When `BELFEM_MPI` is not defined, all communication functions are no-ops and `Communicator` defaults to rank=0, size=1. The test plan splits accordingly.

---

## Two-Tier Test Architecture

### Tier 1 — Non-MPI Tests (always run, single process)

Test utility functions with non-trivial logic, `Communicator` accessors, type mapping, and `CommunicationObject` lifecycle. These compile and run without MPI and go into the normal `make tests` pipeline.

### Tier 2 — MPI Tests (require `mpirun -np N`)

Test actual data transfer: send/receive round-trips, broadcast correctness, distribute/collect pairs, chunking protocol, string serialization. Require a custom `main()` that calls `gComm.init()` / `gComm.finalize()` instead of `gtest_main`. Run via:

```cmake
add_test( NAME Comm_MPI_2
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
            $<TARGET_FILE:tests_comm_mpi> )

add_test( NAME Comm_MPI_4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
            $<TARGET_FILE:tests_comm_mpi> )
```

Tests that require a minimum number of processes should use `GTEST_SKIP()` if `comm_size() < N`.

---

## Message Protocol Summary

Understanding the send/receive protocol is essential for designing MPI tests:

1. **Scalar send/receive:** Single `MPI_Isend`/`MPI_Irecv` with tag `comm_tag(source, target)`.
2. **Container send/receive (Vector, Cell):** Two-phase protocol:
   - Phase 1: Send size with tag T (where T = `comm_tag(sender, receiver)`, incremented with `++`).
   - Phase 2: Send chunked payload, all chunks use tag T+1.
3. **Matrix send/receive:** Three-value size header (rows, cols, capacity), then chunked payload of `capacity()` elements.
4. **Broadcast (Vector, Cell):** Non-blocking broadcast of size, then non-blocking broadcast of data.
5. **String communication:** Serialized into `Vector<char>`, then sent via the Vector protocol.

**Design note on Matrix capacity:** Matrices are transmitted using `capacity()` (the full allocated buffer size) rather than `n_rows() * n_cols()`. This is deliberate: since all MPI ranks use the same matrix backend, transmitting the full capacity is both faster and more reliable than recomputing layout. Tests should verify this contract.

**Tag protocol detail:** `comm_tag(a, b)` returns an *even* number (`2 * (max*size + min) % max_tag`). The `++` increment for the size header uses the next odd number. This deliberate spacing prevents tag collisions between the size message and data messages. The chunked data messages all share the same tag (T+1), relying on MPI ordering guarantees within a single source-destination-tag triple.

---

## Test File Structure

```
tests/comm/
├── test_CommUtils.cpp           # Tier 1: comm_split, comm_tag, comm_splitcount, type mapping
├── test_Communicator.cpp        # Tier 1: Communicator accessors, CommunicationObject
├── test_CommMPI.cpp             # Tier 2: MPI send/receive/broadcast/distribute/collect
```

The Tier 2 file needs a custom `main()`:

```cpp
// test_CommMPI.cpp
#include <gtest/gtest.h>
#include "commtools.hpp"

int main( int argc, char ** argv )
{
    gComm.init( argc, argv );
    ::testing::InitGoogleTest( &argc, argv );
    int tResult = RUN_ALL_TESTS();
    gComm.finalize();
    return tResult;
}
```

---

## 1. Comm Utilities (Tier 1 — Always Run)

**File:** `test_CommUtils.cpp`

### 1.1 comm_split `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CommSplitZeroLength` | `comm_split(0)` → empty Cell (0 chunks) |
| `CommSplitSmallMessage` | Length < 64K → 1 chunk of exactly that length |
| `CommSplitExactBoundary` | Length = 64K → 1 chunk of 64K |
| `CommSplitBoundaryPlusOne` | Length = 64K+1 → 2 chunks: 64K and 1 |
| `CommSplitLargeMessage` | Length = 3×64K + 100 → 4 chunks: three 64K and one of 100 |
| `CommSplitTotalEqualsInput` | Sum of all chunk sizes equals the input length (for various inputs) |
| `CommSplitNoChunkExceedsMax` | Every chunk size ≤ `gMaxCommChunkLength` |

### 1.2 comm_splitcount `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CommSplitcountConsistentWithSplit` | `comm_splitcount(length)` equals `comm_split(length).size() * (comm_size()-1)` |
| `CommSplitcountVectorConsistentWithSplit` | `comm_splitcount(lengths, root)` equals sum of `comm_split(lengths(p)).size()` for all p ≠ root |
| `CommSplitcountZeroLength` | `comm_splitcount(0)` → 0 |

### 1.3 comm_tag `[semantic]` — **MPI-ONLY** (BUG-C1)

**WARNING (BUG-C1):** `comm_tag()` computes `% gComm.max_tag()`. In non-MPI mode, `max_tag() == 0`, so this is a **modulo-by-zero**. All `comm_tag()` callers in `commtools.hpp` are inside `#ifdef BELFEM_MPI` guards, so this is a latent defect — but it means `comm_tag()` tests MUST be MPI-only. Move to Tier 2.

| Test Name | What It Verifies |
|---|---|
| `CommTagSymmetric` | `comm_tag(a, b) == comm_tag(b, a)` for all tested pairs |
| `CommTagEven` | `comm_tag(a, b) % 2 == 0` — leaves room for size-header increment |
| `CommTagDistinctPairs` | `comm_tag(0, 1) != comm_tag(0, 2)` and `comm_tag(1, 2) != comm_tag(0, 2)` |
| `CommTagNonNegative` | `comm_tag(a, b) >= 0` for all tested pairs |
| `CommTagWithinMaxTag` | `comm_tag(a, b) < gComm.max_tag()` |

### 1.4 comm_type Mapping `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CommTypeInt` | `comm_type<int>()` returns correct MPI type (or 0 without MPI) |
| `CommTypeDouble` | `comm_type<double>()` → `MPI_DOUBLE` (or 0) |
| `CommTypeUnsignedLong` | `comm_type<unsigned long int>()` → `MPI_UNSIGNED_LONG` (or 0) |
| `CommTypeChar` | `comm_type<char>()` → `MPI_CHAR` (or 0) |
| `CommTypeBool` | `comm_type<bool>()` → `MPI_CXX_BOOL` (or 0) |
| `CommTypeFloat` | `comm_type<float>()` → `MPI_FLOAT` (or 0) |

**Note:** In non-MPI builds, `comm_type<T>()` returns 0 for all types. In MPI builds, verify against the known `MPI_*` constants.

---

## 2. Communicator Class (Tier 1 — Always Run)

**File:** `test_Communicator.cpp`

### 2.1 Non-MPI Defaults `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `DefaultRankIsZero` | `gComm.rank() == 0` in non-MPI build |
| `DefaultSizeIsOne` | `gComm.size() == 1` |
| `DefaultNodeSizeIsOne` | `gComm.node_size() == 1` |
| `DefaultMaxTagIsZero` | `gComm.max_tag() == 0` in non-MPI build |

### 2.2 Argument Handling `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `ArgumentsCaptured` | After init with argv `{"prog", "--flag", "val"}`, `arguments().size() == 3` |
| `ArgumentStringFormatted` | `argument_string() == "--flag val"` (space-separated, excludes argv[0]) |
| `ExecPathCaptured` | `exec_path()` matches argv[0] |
| `WorkdirNonEmpty` | `workdir()` returns a non-empty string |
| `SetArgumentsOverwrites` | `set_arguments("custom")` → `argument_string() == "custom"` |

### 2.3 Random Engine `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `RandomEngineUsable` | `gComm.random()()` does not throw, returns a value |

### 2.4 CommunicationObject Lifecycle `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `ObjectRegistersOnConstruction` | Creating a `CommunicationObject` increases `gComm.objects().size()` |
| `ObjectIndexMatchesRegistration` | `obj.index()` matches its position in `gComm.objects()` |
| `ObjectFreeNullsSlot` | After `obj.free()`, `gComm.objects()(obj.index()) == nullptr` |
| `MultipleObjectsRegistered` | Create 3 objects → `gComm.objects().size()` increases by 3, each has unique index |

---

## 3. MPI Point-to-Point (Tier 2 — Requires `mpirun`)

**File:** `test_CommMPI.cpp`
**Minimum ranks:** 2

### 3.1 Scalar Send/Receive `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `SendReceiveInt` | Rank 0 sends `int`, rank 1 receives matching value |
| `SendReceiveReal` | Rank 0 sends `real`, rank 1 receives matching value |
| `SendReceiveIndexT` | Rank 0 sends `index_t`, rank 1 receives matching value |
| `SendToSelfIsNoOp` | `send(val, comm_rank())` does not deadlock or corrupt |

### 3.2 Vector Send/Receive `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `SendReceiveVectorSmall` | 10-element `Vector<real>`: data matches after receive |
| `SendReceiveVectorEmpty` | Empty vector: receive produces length-0 vector, no hang |
| `SendReceiveVectorLargeChunked` | 100K-element vector (exceeds 64K chunk): data matches after receive |
| `SendReceiveVectorExactChunkBoundary` | Vector of exactly 64K elements: correct single-chunk transfer |

### 3.3 Cell Send/Receive `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `SendReceiveCellInt` | `Cell<int>{1,2,3,4,5}`: data matches |
| `SendReceiveCellEmpty` | Empty Cell: no hang, size 0 on receiver |

### 3.4 Matrix Send/Receive `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `SendReceiveMatrixSmall` | 3×4 `Matrix<real>`: dimensions and data match |
| `SendReceiveMatrixCapacity` | Verify that `capacity()` bytes are transmitted (not `n_rows*n_cols` — this is the designed behavior) |
| `SendReceiveMatrixEmpty` | 0×0 matrix: no hang |

### 3.5 String Send/Receive `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `SendReceiveString` | "Hello BELFEM" round-trips correctly |
| `SendReceiveEmptyString` | Empty string round-trips to empty string |
| `SendReceiveStringSpecialChars` | String with spaces, newlines, unicode: data preserved |

### 3.6 Raw Array Send/Receive `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `SendReceiveRawArray` | `real[100]` array: data matches on receiver |
| `SendReceiveRawArrayOverflow` | Receiver allocates less than sender transmits → `BELFEM_ERROR` fires |

---

## 4. MPI Collectives (Tier 2 — Requires `mpirun`)

**File:** `test_CommMPI.cpp`
**Minimum ranks:** 2 (some tests require 4)

### 4.1 Broadcast `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `BroadcastScalar` | Root sets value, all ranks have same value after broadcast |
| `BroadcastVector` | Root has 100-element vector, non-root initializes empty → all match after broadcast |
| `BroadcastVectorEmpty` | Root has empty vector → non-root receives empty, no hang |
| `BroadcastCell` | Root has `Cell<int>{1,2,3}` → all ranks match |
| `BroadcastCellString` | Root has `Cell<string>{"foo", "bar"}` → all ranks match |
| `BroadcastRawArray` | `broadcast(T*, root, length)` → all ranks have same data |

### 4.2 Distribute / Collect Pairs `[semantic]` (minimum 4 ranks)

| Test Name | What It Verifies |
|---|---|
| `DistributeCollectScalar` | Each rank distributes its rank number → collect gathers all rank numbers |
| `DistributeCollectVector` | `distribute(Vector<T>)` followed by matching `collect(Vector<T>)` round-trips correctly |
| `DistributeCollectCellVector` | `distribute(Cell<Vector<T>>)` on root → each rank receives its designated vector |
| `DistributeCollectCellCell` | `distribute(Cell<Cell<T>>)` → each rank receives its Cell |
| `DistributeCollectMatrix` | `distribute(Cell<Matrix<T>>)` paired with `collect(Cell<Matrix<T>>)` → data integrity |
| `DistributeEmptyPayload` | One rank gets an empty vector/cell → no hang or crash |

### 4.3 Barrier `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `BarrierDoesNotDeadlock` | `comm_barrier()` returns on all ranks (trivial but confirms wiring) |

---

## 5. Chunking Protocol Tests (Tier 2)

**File:** `test_CommMPI.cpp`

These specifically validate the two-phase size+data protocol and chunk boundaries.

### 5.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `VectorChunkBoundaryMinus1` | Vector of 64K-1 elements → single chunk, correct data |
| `VectorChunkBoundaryExact` | Vector of exactly 64K elements → single chunk |
| `VectorChunkBoundaryPlus1` | Vector of 64K+1 elements → two chunks, correct data |
| `VectorMultipleFullChunks` | Vector of 3×64K elements → three equal chunks |
| `MultipleConsecutiveTransfers` | Send two different vectors in sequence between same rank pair → tag protocol doesn't collide |
| `BidirectionalTransfer` | Rank 0 sends to rank 1, rank 1 sends to rank 0 simultaneously → no tag collision (comm_tag is symmetric, so both directions use the same base tag — verify this works with MPI ordering) |

---

## 6. Debug-Only Tests

### 6.1 Tier 1 (non-MPI) `[debug]`

| Test Name | What It Verifies |
|---|---|
| `DistributeCellWrongSizeThrows` | `distribute(Cell<T>)` with `size() != comm_size()` → assertion |
| `DistributeVectorWrongSizeThrows` | `distribute(Vector<T>)` with `length() != comm_size()` → assertion |

### 6.2 Tier 2 (MPI) `[debug]`

| Test Name | What It Verifies |
|---|---|
| `RawReceiveOverflowThrows` | Receiver pre-allocates smaller buffer than sender's message → `BELFEM_ERROR` fires |
| `DistributeOffsetWrongLengthThrows` | `distribute(data, offsets)` with `offsets.length() != comm_size()+1` → assertion |

---

## 7. What We Do NOT Test

- MPI library correctness (that's OpenMPI's / MPICH's job)
- Exact RNG seed values (implementation detail)
- Environment variable setup for MPI binding (`OMPI_MCA_*`)
- PETSc initialization paths (separate module concern)
- Performance or timing of barriers/transfers
- `comm_check()` with actual MPI error codes (would require injecting MPI failures)

---

## 8. Implementation Notes for Claude Code

1. **Tier 2 tests need a custom `main()`.** Do not link with `gtest_main`. Instead, write a `main()` that calls `gComm.init(argc, argv)`, then `RUN_ALL_TESTS()`, then `gComm.finalize()`.
2. **Use `GTEST_SKIP()`** for tests that require more ranks than available: `if (comm_size() < 4) GTEST_SKIP() << "Requires 4 ranks";`
3. **Matrix capacity is transmitted deliberately.** The protocol sends `capacity()` (full buffer), not `n_rows()*n_cols()`. This is by design because all ranks use the same backend. Tests should verify `capacity()` bytes arrive, not just the logical dimensions.
4. **Tag symmetry matters.** `comm_tag(a, b) == comm_tag(b, a)`, but send uses `comm_tag(myRank, target)` while receive uses `comm_tag(source, myRank)`. For point-to-point, `source` on the receive side equals `myRank` on the send side, so the tags match. Verify this in the bidirectional test.
5. **The size-header protocol uses tag increment.** Send side: size at tag T (via `tCommTag++`), data at tag T+1. Receive side: size at tag T, data at tag T+1 (via `+ 1`). This must stay in sync. The chunking tests exercise this implicitly.
6. **`CommunicationObject` tests** should use a minimal concrete subclass (the base class is not abstract but `free()` is virtual). A trivial test subclass that counts `free()` calls is valuable.
7. **Non-MPI builds:** Most communication functions are empty. Tier 1 tests should focus on the utility functions and Communicator accessors, not on send/receive (which are no-ops).
8. **Bidirectional transfer test:** This is subtle. Both ranks send and receive with the same `comm_tag` base value (since it's symmetric). MPI guarantees message ordering per (source, dest, tag) triple, and the non-blocking sends/receives are waited on. Verify no deadlock and correct data in both directions.

---

## 9. Codex Audit Checklist

When reviewing Claude Code's test implementation, verify:

- [ ] Tier 2 file has custom `main()` with `gComm.init()` / `gComm.finalize()`
- [ ] Tests use `GTEST_SKIP()` for insufficient rank count, not `#ifdef` or assumptions
- [ ] `comm_split` tests verify sum-of-chunks equals input length
- [ ] `comm_tag` tests verify symmetry and evenness
- [ ] Point-to-point tests cover empty, small, exact-boundary, and over-boundary sizes
- [ ] Matrix tests verify `capacity()` bytes transferred (documented design choice)
- [ ] String tests include empty string and special characters
- [ ] Distribute/collect tested as pairs, not in isolation
- [ ] Debug tests use correct guards: `BELFEM_ERROR` tests always run, `BELFEM_ASSERT` tests wrapped in `#ifndef NDEBUG`
- [ ] No tests for MPI internal correctness or timing
- [ ] BELFEM naming conventions (`t` prefix for locals)
- [ ] `EXPECT_EQ` used for integer/rank comparisons, `EXPECT_NEAR` for floating-point data
