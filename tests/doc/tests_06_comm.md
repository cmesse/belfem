# Communication Module Test Suite

**Date:** 2026-03-23
**Revised:** 2026-09-05 (Tier 2 is now registered with ctest; suite table refreshed)

---

## Two-Tier Architecture

### Tier 1 — Always Run (no MPI required)

**Files:** `test_CommUtils.cpp`, `test_Communicator.cpp`
**Main:** `test_comm_main.cpp`
**Run:** `./test_comm`

| Suite | Tests | Coverage |
|-------|-------|----------|
| `CommSplit` | 7 | Zero, small, boundary, boundary+1, large, sum invariant, max-chunk invariant |
| `CommSplitcount` | 3 | Zero, consistency with split, vector overload |
| `CommType` | 6 | MPI constants (int, double, float, char, unsigned long, bool) |
| `Communicator` | 4 | Rank, size, random engine, arguments |
| `CommunicationObject` | 4 | Registration, index, free, multiple objects |
| **Total** | **24** | |

### Tier 2 — MPI Only (launched through `mpirun` by ctest)

**File:** `test_CommMPI.cpp`
**Main:** `test_commmpi_main.cpp` (owns `gComm`/`gLog`, the launcher sentinel and the verdict fold)
**Binary:** `test_commmpi`, registered by `Add_Test.cmake` as `commmpi_np2` and `commmpi_np4`
(`TESTRANKS 2 4` in `tests/comm/CMakeLists.txt`); both are part of `make check`, not `check-fast`,
and carry the `mpi` label so `ctest -LE mpi` skips them by choice.
**Run by hand:** `mpirun --oversubscribe -np 4 ./test_commmpi` with `BELFEM_TESTRANKS=4` in the
environment. The sentinel in `tests/common/tier2_launcher_sentinel.hpp` fails the run when
`comm_size()` differs from that value or is below 2, so a binary launched without `mpirun` cannot
skip its way to green. Every rank's `RUN_ALL_TESTS()` result is MAX-reduced before `finalize()`, so
one failing rank makes the whole ctest line red.

| Suite | Tests | Min Ranks | Coverage |
|-------|-------|-----------|----------|
| `CommTag` | 5 | 2 (1 needs 3+) | Symmetric, even, distinct pairs, non-negative, within max_tag |
| `CommReduce` | 2 | 2 | `allreduce` MIN: smallest rank wins, zero wins |
| `CommSendRecv` | 16 | 2 (1 needs 4) | Scalar int/real, Vector (small/empty/large chunked/exact chunk boundary), Cell (int/empty), Matrix (small/empty/after shrink), String (normal/empty/special), self-send no-op, ring pass 4-rank |
| `CommBroadcast` | 4 | 2 | Scalar, vector, empty vector, Cell of strings |
| `CommDistributeCollect` | 1 | 2 | Scalar distribute/collect |
| `CommBarrier` | 1 | 2 | Deadlock check |
| `CommChunking` | 3 | 2 | Boundary-1, boundary+1, multiple consecutive transfers |
| `ShareReceive` | 3 | 2 | `share`/`receive` Vector and Cell across a chunk boundary, empty Vector |
| **Total** | **35** | | |

---

## @warning API Gotchas

- `gComm.rank()` and `gComm.size()` return `gNoOwner` (sentinel) BEFORE `init()` — not 0/1.
- `comm_split()` chunk size: `constexpr int gMaxCommChunkLength = 64 * 1024` (65536).
- `send(aData, aTarget)` is a no-op when `target == rank` (self-send).
- Matrix send protocol: 3-value header (rows, cols, transfer length = `spacing()*n_cols`), then that many chunked elements. Never `capacity()` — stale-large after a shrink (fixed 2026-08-09).
- String send: serialized to `Vector<char>`, then sent via Vector protocol.
- `CommunicationObject` destructor does NOT call `free()` — tests must call `free()` explicitly to avoid dangling pointers in `gComm.objects()`.

---

## @note Compilation Flags

Tests must be compiled with `-DBELFEM_MPI` to match the library. Without it, `#ifndef BELFEM_MPI` blocks activate and test wrong code paths.
