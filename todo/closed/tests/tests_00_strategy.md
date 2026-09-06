# BELFEM Test Strategy

**Date:** 2026-03-22
**Purpose:** Master test strategy for the BELFEM framework using Google Test
**Audience:** Claude Code (primary implementer), Codex (audit/review), human developers

---

## 1. Philosophy

BELFEM follows **"Performance First, Safety Through Testing"**. Release builds strip all debug checks (`-O3 -DNDEBUG -fno-exceptions`). Debug builds enforce invariants through `BELFEM_ASSERT` and `BELFEM_ERROR`. The test suite is the primary safety net — it must be thorough enough that release builds can be trusted without runtime guards.

**Core testing principles:**

- Tests validate *contracts*, not implementation details of STL wrappers.
- Every public method gets at least one semantic test. Thin STL pass-throughs (e.g., `Cell::begin()`) get smoke tests; BELFEM-authored logic (e.g., `DynamicBitset::count()`, `ShiftRegister::revert()`) gets full coverage.
- No performance benchmarks in unit tests. Zero-overhead behavior is validated through assembly inspection or dedicated benchmark binaries, not brittle timing checks.
- No thread-safety tests. BELFEM is deliberately single-threaded internally (MPI for parallelism).

---

## 2. Two-Lane Test Architecture

BELFEM's two-tier error system requires two categories of tests:

### 2.1 Semantic Tests `[semantic]`

Always run, in both debug and release builds. These verify correct behavior under valid inputs: construction, mutation, access, algorithms, operators, copy/move semantics.

### 2.2 Debug Invariant Tests `[debug]`

Run **only** in debug builds. These verify that `BELFEM_ASSERT` fires on invalid inputs (out-of-bounds access, locked-bitset mutation, etc.).

**VERIFIED 2026-03-22 (source inspection of `assert.hpp` lines 76–167, `assert.cpp` lines 289–297):**

Both `BELFEM_ASSERT` and `BELFEM_ERROR` call `belfem::assert::error()`, which in debug builds (`!NDEBUG || DEBUG`) executes `throw std::runtime_error`. In release builds, it calls `error_abort()` → `std::abort()` (or `MPI_Abort`). `BELFEM_ASSERT` is compiled out entirely in release builds; `BELFEM_ERROR` is always active.

**Conclusion:** In debug builds (the default, `USE_DEBUG=ON`), use:

```cpp
EXPECT_THROW( expression, std::runtime_error );
```

for both `BELFEM_ASSERT` and `BELFEM_ERROR` tests. `EXPECT_DEATH` is NOT needed.

**Note:** The `error()` function prints a dragon-adorned error box to stderr before throwing. This is cosmetic noise in test output and does not affect correctness.

### 2.3 Memory Correctness Tests `[valgrind]`

Not a separate test category in code, but a CI annotation. Any test exercising manual memory (`DynamicBitset`, `ShiftRegister`, `StringList`) must pass Valgrind / ASan cleanly. Tag these in the test plan so CI knows which test binaries to run under Valgrind.

---

## 3. Conventions

### 3.1 File Layout

One test file per container or module. File naming:

```
tests/containers/test_Cell.cpp
tests/containers/test_DynamicBitset.cpp
tests/core/test_Logger.cpp
...
```

Each test file is self-contained: it includes only `<gtest/gtest.h>` and the BELFEM headers it tests. No cross-test-file dependencies.

### 3.2 Naming

Follow BELFEM naming inside test bodies:

- Use `t` prefix for local variables: `tCell`, `tBitset`, `tResult`
- Use BELFEM types: `index_t`, `real`, `uint`, `string` (not `std::string` unless testing interop)
- Test names describe the contract being verified:

```cpp
TEST( CellSemantic, PushIncreasesSize )
TEST( CellSemantic, MoveConstructorLeavesSourceEmpty )
TEST( CellDebug, OutOfBoundsAccessThrows )
```

**Test group naming convention:**

- `ContainerNameSemantic` — always-run tests
- `ContainerNameDebug` — debug-only tests (wrapped in `#ifndef NDEBUG`)
- `ContainerNameTyped` — typed/parameterized test suites

### 3.3 Typed and Parameterized Tests

For template containers, use typed tests to cover the type matrix:

| Container | Types to test |
|---|---|
| `Cell<T>` | `int`, `real`, `index_t`, `std::string` |
| `Cell<T*>` | `int*` (non-owning pointer semantics) |
| `Bitset<N>` | `N = 8, 64, 128, 256` |
| `ShiftRegister<T>` | `int`, `real` |
| `Map<K,V>` / `OrderedMap<K,V>` | `<string, int>`, `<index_t, real>`, `<int, string>` |
| `Set<Key>` | `int`, `index_t`, `string` |

For `DynamicBitset`, use value-parameterized tests over sizes: `{1, 7, 63, 64, 65, 127, 128, 1024}`.

### 3.4 Debug Guards

Wrap all debug-only tests:

```cpp
#ifndef NDEBUG
TEST( CellDebug, OutOfBoundsAccessThrows )
{
    belfem::Cell< int > tCell( 5, 0 );
    EXPECT_THROW( tCell( 10 ), std::runtime_error );  // or EXPECT_DEATH
}
#endif
```

This ensures release-build CI never runs death/throw tests that would either pass vacuously or crash.

### 3.5 Coverage Targets

| Category | Target |
|---|---|
| BELFEM-authored logic (manual memory, bit arithmetic, algorithms) | Full line + branch coverage |
| Thin STL wrappers (`Cell::push` → `vector::push_back`) | One smoke test per method |
| `print()` / `to_string()` output formatting | Smoke test only (stdout capture) |
| Error paths (`BELFEM_ASSERT`, `BELFEM_ERROR`) | One test per distinct error message |

---

## 4. CMake Integration

Tests integrate with the existing `make tests` target. Recommended CMake additions:

```cmake
# One executable per module
add_executable( tests_containers
    tests/containers/test_Cell.cpp
    tests/containers/test_DynamicBitset.cpp
    tests/containers/test_Bitset.cpp
    tests/containers/test_ShiftRegister.cpp
    tests/containers/test_Map.cpp
    tests/containers/test_Set.cpp
    tests/containers/test_Queue.cpp
    tests/containers/test_StringList.cpp
)
target_link_libraries( tests_containers gtest gtest_main belfem_containers belfem_core )
add_test( NAME Containers COMMAND tests_containers )
```

For Valgrind runs:

```cmake
add_test( NAME Containers_Valgrind
    COMMAND valgrind --leak-check=full --error-exitcode=1 $<TARGET_FILE:tests_containers> )
```

---

## 5. Priority and Phasing

### Phase 1 — Highest Value

| Module | Containers | Rationale |
|---|---|---|
| containers | `Cell<T>`, `DynamicBitset` | Most used, most complex, highest bug risk |

### Phase 2 — Core Infrastructure

| Module | Containers | Rationale |
|---|---|---|
| containers | `ShiftRegister<T>`, `Bitset<N>` | Manual memory (ShiftRegister), used in time-stepping and Genome |
| containers | `Map`, `OrderedMap`, `Set` | Lookup tables throughout framework |

### Phase 3 — Completeness

| Module | Containers | Rationale |
|---|---|---|
| containers | `Queue<T>`, `StringList` | Thin wrappers, low risk |
| containers | `Genome<B,N>` | Application-level, exercises Bitset/Vector |

### Future Modules (separate test plan files)

- `tests_2_core.md` — Logger, Timer, Arguments, assert machinery
- `tests_3_linalg.md` — Vector, Matrix, BLAS/LAPACK wrappers
- `tests_4_comm.md` — MPI communication layer
- `tests_5_mesh.md` — Mesh data structures and I/O
- (and so on per module)

---

## 6. Workflow: Claude Code → Codex

1. **Claude Code** reads this strategy file and the container-specific test plan (`tests_1_containers.md`), then implements the test files.
2. **Claude Code** builds and runs the tests, fixes any compilation or logic errors.
3. **Codex** audits the test files against the test plan: checks for missing tests, incorrect assertions, and naming convention violations.
4. **Codex** reports findings in `./todo/ai_exchange.md` per the collaboration protocol.
5. Human reviews and approves.

---

## 7. Things We Explicitly Do NOT Test

- Thread safety (BELFEM is MPI-based, not threaded internally)
- Performance regressions (use benchmarks, not unit tests)
- Platform-specific behavior (tests assume Linux x86_64)
- Private/internal methods (test through public API only)
- `Genome<B,N>` encode/decode internals (test through `set_values`/`get_values` round-trip)

---

## 8. Reference

- **Coding philosophy:** `doc/coding_philosophy.md`
- **Error handling:** `assert.hpp` (two-tier: `BELFEM_ASSERT` debug-only, `BELFEM_ERROR` always-on)
- **Container headers:** `src/containers/cl_*.hpp`
- **AI collaboration protocol:** `doc/ai_collaboration_protocol.md`
