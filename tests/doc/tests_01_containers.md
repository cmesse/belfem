# Container Module Test Suite

**Date:** 2026-03-23

---

## Test Count by File

| File | Suites | Tests |
|------|--------|-------|
| `test_Cell.cpp` | `Cell` (typed ×4), `CellDebug`, `Cell` (non-typed) | 121 |
| `test_DynamicBitset.cpp` | `Sizes/DynamicBitset` (param ×8), `DynamicBitset`, `DynamicBitsetDebug` | 176 |
| `test_ShiftRegister.cpp` | `ShiftRegister` (typed ×2) | 72 |
| `test_Bitset.cpp` | `Bitset` (typed ×4) | 68 |
| `test_Map.cpp` | `Map` (typed ×4), `MapDebug`, `OrderedMapDebug`, `Map`, `OrderedMap` | 63 |
| `test_Set.cpp` | `Set` | 48 |
| `test_Queue.cpp` | `Queue` | 14 |
| `test_StringList.cpp` | `StringList`, `StringListDebug` | 11 |
| `test_Genome.cpp` | `Genome` | 14 |
| **Total** | **28 suites** | **595** |

---

## Naming Conventions

- All suites use the container class name directly: `Cell`, `DynamicBitset`, `Map`, etc.
- Debug-only tests (inside `#ifndef NDEBUG`) use a `Debug` suffix: `CellDebug`, `DynamicBitsetDebug`.
- No `using namespace belfem;` — all types explicitly qualified.
- `t` prefix for locals per BELFEM convention.

---

## Key Design Decisions

### Typed Tests
- `Cell`: typed over `int, belfem::real, belfem::index_t, std::string` via `CellTestTraits`.
- `ShiftRegister`: typed over `int, belfem::real` only (no `std::string` — UB with malloc/free).
- `Bitset`: typed over `N = 8, 64, 128, 256` using `std::integral_constant`.
- `Map`: typed over 4 tag types covering `{Map, OrderedMap} x {<string,int>, <index_t,real>}`.
- `DynamicBitset`: value-parameterized over sizes `{1, 7, 63, 64, 65, 127, 128, 1024}`.
  §2.12 adds fixed-size tests at 4095/4096/4097 and 262144/262145 as well: `where()` and
  `reset()` walk a two-level summary bitmap, and a bitset needs more than 4096 bits before it
  uses a second level-1 word and more than 262144 before it uses a second level-2 word — every
  size in the parameterized set is below both thresholds and leaves the walk unexercised.

### Move Semantics
- STL-backed wrappers (Cell, Map, Set, Queue): moved-from source checked for validity but NOT asserted empty (STL does not guarantee empty after move).
- Manual-memory containers (DynamicBitset, ShiftRegister): moved-from source asserted null/empty (BELFEM guarantees this).

---

## Source Fixes

Source improvements made during test development are documented in devlog/dl20260324_test_suite_bugs.md.

---

## Codex/Junie Review Items

### Addressed
- [x] DynamicBitset §2.12 summary-bitmap suite added; `flip()` zero-size underflow, the
      non-compiling portable ctz fallback, and the shifted-index fallback bug fixed (2026-08-18)
- [x] DynamicBitset mutable `data()` made private — a write through it bypasses the summary
      bitmaps. `CountIgnoresTrailingBits` no longer corrupts the last block by hand; it asserts
      the invariant (no public mutator can set a bit at or beyond `size()`) instead (2026-08-18)
- [x] Cell `SetSizeGrows` now checks default initialization of appended slots (Codex)
- [x] Cell `InsertMoveVersion` smoke test added (Codex)
- [x] DynamicBitset `ToHexEmptyBitset` test added (Codex)
- [x] DynamicBitset `where()` tests verify actual indices, not just size (Codex)
- [x] Map copy/move assignment tests added (Codex)
- [x] Map `get_entry()` uniqueness check added (Codex)
- [x] Set move-assignment source validity check added (Codex)
- [x] Queue `SelfMoveAssignment` fixed (was asserting non-portable behavior) (Codex)
- [x] Queue `MoveConstructor` fixed (was asserting empty on moved-from std::queue) (Codex)
- [x] Genome `RandomizeLogScaleBranch` test added (Codex)

### Accepted as-is
- Cell `PrintRealSpecialization` "e+" check — adequate for smoke test per strategy
- DynamicBitset move tests assert pointer identity — intentional for manual-memory containers
  (they read the pointer through the const `data()`; the mutable overload is private, see below)

---

## API Gotchas

- `Cell::first()`/`last()` throw `std::out_of_range` (via `vector::at()`), NOT `std::runtime_error`. `operator()` throws `std::runtime_error` (via `BELFEM_ASSERT`).
- `ShiftRegister::fill()` operates on `mSize`, NOT `mCapacity`. Calling on empty register is a no-op.
- `ShiftRegister` uses `malloc`/`free` — do NOT use with `std::string` or other non-trivially-copyable types (except whitelisted `Vector<T>`, `Matrix<T>`).
- `DynamicBitset` has no default constructor — cannot be a fixture member.
- `StringList` is non-copyable, non-movable.
- `Genome` constructor requires `(const Vector<real>&, const Vector<real>&, const Bitset<N>&)` — no default constructor.
- `operator==` on `DynamicBitset` requires both operands locked.
- Bitwise operators (`|`, `&`, `^`) on `DynamicBitset` do NOT require locked state; compound variants (`|=`, `^=`, `&=`) reject locked operands.

---

## Future Work

- [ ] Thread-safety stress tests (not applicable — BELFEM is MPI-based)
- [ ] Valgrind CI annotation for DynamicBitset, ShiftRegister, StringList
- [ ] Performance benchmarks (separate from unit tests)
