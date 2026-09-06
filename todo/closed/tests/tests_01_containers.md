# BELFEM Container Tests — Detailed Plan

**Date:** 2026-03-22
**Purpose:** Method-level test matrix for all containers in `src/containers/`
**Depends on:** `tests_0_strategy.md` (read that first for conventions and architecture)
**Confidence:** High on API surface and semantics (based on full source review). Medium on exact error mechanism (`EXPECT_THROW` vs `EXPECT_DEATH`) — verify per strategy doc instructions.

---

## Container Inventory

| Container | Header | Wraps | Memory | Priority |
|---|---|---|---|---|
| `Cell<T>` | `cl_Cell.hpp` | `std::vector<T>` | Automatic (vector) | Phase 1 |
| `DynamicBitset` | `cl_DynamicBitset.hpp` | — (raw `uint64_t*`) | Manual (`malloc`/`free`) | Phase 1 |
| `ShiftRegister<T>` | `cl_ShiftRegister.hpp` | — (raw `T*`) | Manual (`malloc`/`free`) | Phase 2 |
| `Bitset<N>` | `cl_Bitset.hpp` | `std::bitset<N>` | Automatic (bitset) | Phase 2 |
| `Map<K,V>` | `cl_Map.hpp` | `std::unordered_map` | Automatic | Phase 2 |
| `OrderedMap<K,V>` | `cl_OrderedMap.hpp` | `std::map` | Automatic | Phase 2 |
| `Set<Key>` | `cl_Set.hpp` | `std::unordered_set` | Automatic | Phase 2 |
| `Queue<T>` | `cl_Queue.hpp` | `std::queue` | Automatic | Phase 3 |
| `StringList` | `cl_StringList.hpp/.cpp` | — (raw `char**`) | Manual (`malloc`/`free`) | Phase 3 |
| `Genome<B,N>` | `cl_Genome.hpp` | Uses `Bitset<N*B>`, `Vector<real>` | Automatic | Phase 3 (deferred) |

**Note:** There is no `OrderedSet` class in the codebase. Some external analyses mention it — ignore those references.

---

## 1. Cell\<T\>

**File:** `test_Cell.cpp`
**Approach:** Typed tests over `int`, `real`, `index_t`, `std::string`. Separate non-typed tests for `Cell<T*>` pointer semantics and `Cell<real>::print` specialization.

### 1.1 Construction & Destruction `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `DefaultConstructorIsEmpty` | `size() == 0`, `empty() == true` |
| `SizedConstructorWithValue` | `Cell<int>(5, 42)` → `size() == 5`, all elements are `42` |
| `InitializerListConstructor` | `Cell<int>{1,2,3}` → `size() == 3`, correct element values |
| `CopyConstructorDeepCopies` | Modify copy, verify original unchanged |
| `MoveConstructorTransfersOwnership` | After move, target has original data. Source remains valid (do NOT assert empty — STL-backed wrappers do not guarantee moved-from state) |
| `CopyAssignmentDeepCopies` | Assign to existing Cell, verify independence |
| `MoveAssignmentTransfersOwnership` | Assign via move, target has data. Source remains valid (do NOT assert empty) |
| `SelfCopyAssignment` | `tCell = tCell` does not corrupt state |
| `SelfMoveAssignment` | `tCell = std::move(tCell)` does not corrupt state |

### 1.2 Element Access `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `ParenthesisOperatorReadWrite` | `tCell(i) = x` then `tCell(i) == x` |
| `DataPointerMatchesFirstElement` | `tCell.data() == &tCell(0)` |
| `VectorDataReturnsInternalVector` | `tCell.vector_data().size() == tCell.size()` |
| `FirstAndLast` | `first()` and `last()` return correct elements |
| `ConstAccess` | `const Cell<int>& ref = tCell; ref(0)` compiles and returns correctly |

### 1.3 Element Access `[debug]`

| Test Name | What It Verifies |
|---|---|
| `OutOfBoundsThrows` | `tCell(size)` triggers assertion |
| `FirstOnEmptyThrows` | `Cell<int>().first()` triggers assertion |
| `LastOnEmptyThrows` | `Cell<int>().last()` triggers assertion |

### 1.4 Mutation `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `PushIncreasesSize` | Size goes from `n` to `n+1`, last element matches |
| `PushMoveVersion` | `push(std::move(str))` moves without copying |
| `EmplaceConstructsInPlace` | `emplace(args...)` increases size, element is correct |
| `PopReturnsLastAndShrinks` | `pop()` returns last element, size decreases by 1 |
| `PopFromSingleElement` | Cell with 1 element → `pop()` → `empty() == true` |
| `SetSizeGrows` | `set_size(larger)` pads with default values |
| `SetSizeShrinks` | `set_size(smaller)` truncates |
| `SetSizeWithValue` | `set_size(n, val)` fills new slots with `val` |
| `InsertAtBegin` | Insert at `begin()`, verify element at index 0 |
| `InsertAtMiddle` | Insert at `begin() + k`, verify shift |
| `EraseAtPosition` | Erase at iterator, verify size and remaining elements |
| `EraseRange` | Erase `[first, last)`, verify size and content |
| `ClearMakesEmpty` | `clear()` → `size() == 0`, `empty() == true` |
| `ReserveDoesNotChangeSize` | `reserve(100)` → `size()` unchanged, `capacity() >= 100` |
| `ShrinkToFit` | After `pop` operations, `shrink_to_fit()` reduces capacity |
| `SwapExchangesContents` | `tA.swap(tB)` swaps sizes and data |

### 1.5 Free Functions `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `SortAscending` | `belfem::sort(tCell)` produces ascending order |
| `SortWithComparator` | `belfem::sort(tCell, comp)` uses custom comparator |
| `SortPartial` | `belfem::sort(tCell, comp, n)` only sorts first `n` elements |
| `UniqueRemovesDuplicates` | `{5,2,5,1,2}` → `{1,2,5}` (sorted + unique) |
| `UniqueOnAlreadyUnique` | No-op, size unchanged |
| `ReverseFlipsOrder` | `{1,2,3}` → `{3,2,1}` |
| `AppendConcatenates` | `append(tA, tB)` → `tA.size() == sizeA + sizeB`, `tB` unchanged |
| `AppendMoveClearsSource` | `append_move(tA, tB)` → `tB.size() == 0` |
| `SwapFreeFunction` | `belfem::swap(tA, tB)` exchanges contents |

### 1.6 Pointer Semantics `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CellOfPointersNonOwning` | `Cell<int*>` with raw pointers: `clear()` does not free memory (verify pointed-to values survive) |

### 1.7 Print `[semantic]` (smoke test only)

| Test Name | What It Verifies |
|---|---|
| `PrintIntProducesOutput` | Capture stdout, verify non-empty output |
| `PrintRealSpecialization` | `Cell<real>::print()` uses `%+.15e` format |

### 1.8 Iteration `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `RangeBasedForLoop` | `for (auto& x : tCell)` visits all elements |
| `ConstIteration` | `for (const auto& x : tConstCell)` compiles and works |

---

## 2. DynamicBitset

**File:** `test_DynamicBitset.cpp`
**Approach:** Value-parameterized tests over sizes `{1, 7, 63, 64, 65, 127, 128, 1024}`. Separate fixed tests for operators, lock/hash, conversions. All tests are `[valgrind]` candidates.

### 2.1 Construction & Destruction `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `ConstructorInitializesAllZero` | `count() == 0` for all parameterized sizes |
| `SizeMatchesConstructorArg` | `size() == aNumberOfBits` |
| `MemorySizeCorrect` | `memory() == (aNumberOfBits + 63) / 64` |
| `CopyConstructorDeepCopies` | Modify copy, original unchanged |
| `CopyPreservesLockState` | Copy of locked bitset is also locked with same hash |
| `MoveConstructorTransfers` | Source has `size() == 0` and null `data()` after move |
| `MoveConstructorPreservesHash` | Moved-to bitset retains hash from source |
| `CopyAssignmentDifferentSizes` | Assign larger to smaller, verify reallocation |
| `CopyAssignmentSameSize` | Assign same-size, verify no reallocation (memory size unchanged) |
| `MoveAssignment` | Source is nullified, target has correct data |
| `SelfCopyAssignment` | `tBs = tBs` does not corrupt |
| `SelfMoveAssignment` | `tBs = std::move(tBs)` does not corrupt |

### 2.2 Bit Manipulation `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `SetAndTest` | `set(i)` → `test(i) == true` for various `i` |
| `ResetSingleBit` | `set(i)` then `reset(i)` → `test(i) == false` |
| `ResetAll` | Set several bits, `reset()` → `count() == 0` |
| `ResetAllUnlocks` | Locked bitset → `reset()` → `is_locked() == false` |
| `FlipSingleBit` | `flip(i)` toggles from 0→1 and 1→0 |
| `FlipAll` | `flip()` on known pattern, verify inversion |
| `SetWithBoolValue` | `set(i, true)` and `set(i, false)` |
| `BlockBoundaryBits` | Specifically test bits 62, 63, 64, 65 (straddles uint64 boundary) |
| `LastBitInOddSize` | Size 65: set bit 64, verify `test(64) == true` and `count() == 1` |

### 2.3 Counting `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CountAllZeros` | `count() == 0` on fresh bitset |
| `CountAllOnes` | Set all bits, `count() == size()` |
| `CountSingleBit` | Set one bit, `count() == 1` |
| `CountSparse` | Set bits 0, 32, 64: `count() == 3` (crosses block boundaries) |
| `CountIgnoresTrailingBits` | For size 65, the bits 65–127 in the second block must not be counted |

### 2.4 Lock / Unlock / Hash `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `LockSetsLockedState` | `lock()` → `is_locked() == true` |
| `UnlockClearsLockedState` | `lock()` then `unlock()` → `is_locked() == false` |
| `HashDeterministic` | Lock same bitset twice (unlock between), get same hash |
| `EqualBitsetsEqualHash` | Two bitsets with same content produce same hash after lock |
| `DifferentBitsetsDifferentHash` | Different content → different hash (with high probability) |
| `IndexSurvivesLock` | `set_index(42)` then `lock()` → `index() == 42` |

### 2.5 Lock / Unlock / Hash `[debug]`

| Test Name | What It Verifies |
|---|---|
| `SetOnLockedThrows` | `lock()` then `set(0)` triggers assertion |
| `ResetBitOnLockedThrows` | `lock()` then `reset(i)` triggers assertion |
| `FlipOnLockedThrows` | `lock()` then `flip(i)` triggers assertion |
| `HashOnUnlockedThrows` | `hash()` without `lock()` triggers assertion |
| `CompoundOrOnLockedThrows` | `lock()` then `\|=` triggers `BELFEM_ERROR` |
| `CompoundXorOnLockedThrows` | `lock()` then `^=` triggers `BELFEM_ERROR` |
| `CompoundAndOnLockedThrows` | `lock()` then `&=` triggers `BELFEM_ERROR` |
| `OutOfBoundsSetThrows` | `set(size())` triggers assertion |
| `OutOfBoundsTestThrows` | `test(size())` triggers assertion |

### 2.6 Bitwise Operators `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `OrOperator` | `A \| B` sets union of bits |
| `OrAssignOperator` | `A \|= B` modifies A in place |
| `AndOperator` | `A & B` sets intersection of bits |
| `AndAssignOperator` | `A &= B` modifies A in place |
| `XorOperator` | `A ^ B` sets symmetric difference |
| `XorAssignOperator` | `A ^= B` modifies A in place |
| `OperatorWithSelf` | `A \| A == A`, `A & A == A`, `A ^ A` is all zeros |
| `OperatorWithAllZeros` | `A \| zeros == A`, `A & zeros == zeros` |
| `OperatorWithAllOnes` | `A & ones == A`, `A \| ones == ones` |

### 2.7 Bitwise Operators `[debug]`

| Test Name | What It Verifies |
|---|---|
| `OrMismatchedSizeThrows` | `DynamicBitset(32) \| DynamicBitset(64)` triggers `BELFEM_ERROR` |
| `AndMismatchedSizeThrows` | Same for `&` |
| `XorMismatchedSizeThrows` | Same for `^` |

### 2.8 Conversions `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `ToStringBinaryRepresentation` | 4-bit bitset with bits 0,2 set → `"0101"` (MSB-first) |
| `ToHexRoundTrip` | `set_from_hex(to_hex())` produces identical bitset |
| `ToHexEmptyBitset` | Size 0 → `""` (verify no crash) |
| `SetFromHexKnownPattern` | Known hex string → verify individual bits |
| `ToIntPartial` | Size < 64: `to_int()` returns correct integer |
| `ToIntFull` | Size == 64: `to_int()` returns correct integer |
| `ToIntAllZeros` | `to_int() == 0` |
| `ToRawString` | Verify `to_raw_string()` output format |

### 2.9 Conversions `[debug]`

| Test Name | What It Verifies |
|---|---|
| `ToIntTooLargeThrows` | Size > 64: `to_int()` triggers `BELFEM_ERROR` |

### 2.10 Where (Index Extraction) `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `WhereEmptyBitset` | `where(aBits)` on all-zero → `aBits.size() == 0` |
| `WhereAllSet` | All bits set → `aBits` contains `{0, 1, ..., size()-1}` |
| `WhereSparseMode` | `where(aBits, true)` with 3 bits set → `aBits.size() == 3` |
| `WhereDenseMode` | `where(aBits, false)` with 3 bits set → same result as sparse |
| `WhereSparseVsDenseConsistency` | Both modes produce identical output for same bitset |
| `WhereZeroSizeBitset` | Size 0 bitset → `aBits` is cleared |

### 2.11 Index `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `DefaultIndexIsNoIndex` | `index() == gNoIndex` after construction |
| `SetIndexAndRetrieve` | `set_index(42)` → `index() == 42` |
| `IndexSurvivesCopy` | Copy preserves index |
| `IndexSurvivesMove` | Move preserves index, source resets to `gNoIndex` |

---

## 3. ShiftRegister\<T\>

**File:** `test_ShiftRegister.cpp`
**Approach:** Typed tests over `int`, `real`. All tests are `[valgrind]` candidates due to manual memory.

### 3.1 Construction & Destruction `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CapacityConstructor` | `capacity() == N`, `size() == 0`, `empty() == true` |
| `CapacityWithFillValue` | `capacity() == N`, `size() == N`, all elements equal fill value |
| `InitializerListConstructor` | `{1.0, 2.0, 3.0}` → `size() == 3`, `(0) == 1.0` (newest first) |
| `CopyConstructorDeepCopies` | Modify copy, original unchanged |
| `CopyDoesNotPreserveRevertState` | Copy after push → copy cannot revert |
| `MoveConstructor` | Source nullified, target has data |
| `CopyAssignment` | Deep copy, verify independence |
| `MoveAssignment` | Source nullified |
| `SelfCopyAssignment` | No corruption |
| `SelfMoveAssignment` | No corruption |

### 3.2 Push & Shift `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `PushToEmpty` | Push 1 value → `size() == 1`, `(0) == value` |
| `PushShiftsRight` | Push {A, B, C} → `(0) == C`, `(1) == B`, `(2) == A` |
| `PushBeyondCapacity` | Capacity 3, push 4 values → oldest drops, `size() == 3` |
| `PushConstRef` | `push(const T&)` works correctly |
| `PushLvalueRef` | `push(T&)` works correctly |
| `SizeNeverExceedsCapacity` | Push 100 values into capacity-5 register → `size() == 5` |
| `FullAndEmptyStates` | `full()` and `empty()` reflect actual state |

### 3.3 Revert `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `RevertUndoesPush` | Push A then B, revert → `(0) == A`, size decreases |
| `RevertWhenFull` | Full register, push, revert → recovers dropped value |
| `RevertRestoresSize` | If register was not full before push, revert decreases size |
| `RevertMakesSubsequentRevertFail` | Revert, then revert again → triggers `BELFEM_ERROR` |
| `PushAfterRevertEnablesNewRevert` | Push, revert, push → can revert again |

### 3.4 Revert `[debug]`

| Test Name | What It Verifies |
|---|---|
| `RevertWithoutPushThrows` | Fresh register → `revert()` triggers `BELFEM_ERROR` |
| `DoubleRevertThrows` | Push, revert, revert → triggers `BELFEM_ERROR` |

### 3.5 Access `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `ParenthesisOperatorReadWrite` | `(i) = x` then `(i) == x` |
| `ConstAccess` | `const ShiftRegister& ref` → `ref(i)` works |
| `DataPointerValid` | `data()` is non-null after construction |

### 3.6 Access `[debug]`

| Test Name | What It Verifies |
|---|---|
| `OutOfBoundsThrows` | `(size())` triggers assertion |

### 3.7 Other `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `ClearResetsSize` | `clear()` → `size() == 0`, `empty() == true` |
| `ClearResetsRevertState` | `clear()` after push → cannot revert |
| `FillSetsAllElements` | `fill(42)` → all elements are 42 |
| `FillResetsRevertState` | `fill()` → cannot revert |
| `ReserveDestroysPrevious` | `reserve(newCap)` → `size() == 0`, old data gone |
| `ReserveSameCapacityNoOp` | `reserve(cap)` on same cap → early return, size unchanged |
| `FreeReleasesMemory` | `free()` → `data() == nullptr`, `size() == 0`, `capacity() == 0` |
| `IteratorRange` | `begin()` to `end()` covers exactly `size()` elements |

---

## 4. Bitset\<N\>

**File:** `test_Bitset.cpp`
**Approach:** Typed tests over `N = 8, 64, 128, 256` using `std::integral_constant`.

### 4.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `DefaultConstructorAllZero` | `count() == 0` |
| `SizeEqualsN` | `size() == N` |
| `SetAndTest` | `set(i)` → `test(i) == true` |
| `ResetSingleBit` | `set(i)`, `reset(i)` → `test(i) == false` |
| `ResetAll` | Set several, `reset()` (no arg) — note: `Bitset` has no parameterless `reset()`; only per-index. Verify this. If absent, skip. |
| `FlipBit` | `flip(i)` toggles state |
| `CountMatchesManual` | Set known bits, verify `count()` |
| `CopyConstructorDeepCopies` | Modify copy, original unchanged |
| `MoveConstructor` | Move semantics work (though bitset move is effectively copy) |
| `CopyAssignment` | Self-assignment safe, deep copy |
| `MoveAssignment` | Self-assignment safe |
| `EqualityOperator` | Same bits → `==` true, different → `==` false |
| `InequalityOperator` | Negation of equality |
| `DataExposesStdBitset` | `data()` returns reference to internal `std::bitset` |

---

## 5. Map\<Key, Value\> and OrderedMap\<Key, Value\>

**File:** `test_Map.cpp`
**Approach:** Separate test groups for `Map` and `OrderedMap`. Test with `<string, int>` and `<index_t, real>`.

### 5.1 Shared Tests (apply to both) `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `DefaultConstructorEmpty` | `size() == 0`, `empty() == true` |
| `InsertViaSubscript` | `tMap["key"] = val` → `size() == 1` |
| `LookupViaParenthesis` | `tMap("key") == val` |
| `KeyExists` | `key_exists("key") == true`, `key_exists("missing") == false` |
| `EraseKey` | `erase_key("key")` → `key_exists("key") == false`, size decreases |
| `Clear` | `clear()` → `empty() == true` |
| `CopySemantics` | Copy is independent |
| `MoveSemantics` | Source remains valid after move (do NOT assert empty — STL-backed) |
| `OverwriteExistingKey` | `tMap["k"] = 1; tMap["k"] = 2` → `tMap("k") == 2` |
| `FindReturnsIterator` | `find(key)` returns valid iterator for existing key, `end()` for missing |
| `GetEntryByIndex` | `get_entry(0)` returns a valid key-value pair (**NOTE:** `get_entry()` exists on `Map` only, NOT on `OrderedMap`) |
| `Iteration` | Range-based for loop visits all entries |

### 5.2 Shared Tests `[debug]`

| Test Name | What It Verifies |
|---|---|
| `LookupMissingKeyThrows` | `tMap("nonexistent")` triggers assertion |
| `GetEntryOutOfBoundsThrows` | `get_entry(size())` triggers assertion (**Map only** — OrderedMap does not have `get_entry()`) |

### 5.3 OrderedMap-Specific `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `IterationOrderIsSorted` | Insert keys `{3, 1, 2}`, iterate → keys come out `{1, 2, 3}` |

### 5.4 KeyToString `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `StringKeyToString` | `map::KeyToString(string("hello")) == "hello"` |
| `IntKeyToString` | `map::KeyToString(42) == "42"` |
| `UnsignedKeyToString` | `map::KeyToString(42u) == "42"` |
| `LongUnsignedKeyToString` | `map::KeyToString(42LU) == "42"` |
| `UnknownTypeReturnsUnknown` | `map::KeyToString(3.14)` → `"unknown"` |

---

## 6. Set\<Key\>

**File:** `test_Set.cpp`
**Approach:** Test with `int` and `string` keys.

### 6.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `DefaultConstructorEmpty` | `size() == 0`, `empty() == true` |
| `InitializerListConstructor` | `Set<int>{1,2,3}` → `size() == 3` |
| `InsertNewElement` | `insert(42)` → `second == true`, `size()` increases |
| `InsertDuplicate` | `insert(42)` twice → second returns `false`, `size()` unchanged |
| `InsertMoveSemantics` | `insert(std::move(str))` works |
| `ContainsAndKeyExists` | `contains(x) == key_exists(x)` for present and absent keys |
| `Count` | `count(present) == 1`, `count(absent) == 0` |
| `EraseByKey` | `erase(key)` → returns 1, key no longer found |
| `EraseAbsentKey` | `erase(absent)` → returns 0, size unchanged |
| `EraseByIterator` | `erase(find(key))` removes element |
| `Clear` | `clear()` → `empty() == true` |
| `Reserve` | `reserve(100)` does not change size |
| `Swap` | `tA.swap(tB)` exchanges contents |
| `CopySemantics` | Copy is independent |
| `MoveSemantics` | Source remains valid after move (do NOT assert empty — STL-backed) |
| `Iteration` | Range-based for loop visits all elements |
| `SetDataExposesUnderlying` | `set_data()` returns reference to `unordered_set` |

### 6.2 Set Algebra `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `UnionOperator` | `{1,2} \| {2,3}` → `{1,2,3}` |
| `IntersectionOperator` | `{1,2,3} & {2,3,4}` → `{2,3}` |
| `DifferenceOperator` | `{1,2,3} - {2,4}` → `{1,3}` |
| `SymmetricDifferenceOperator` | `{1,2,3} ^ {2,3,4}` → `{1,4}` |
| `UnionWithEmpty` | `A \| {} == A` |
| `IntersectionWithEmpty` | `A & {} == {}` |
| `UnionWithSelf` | `A \| A == A` |
| `IntersectionWithSelf` | `A & A == A` |
| `DifferenceWithSelf` | `A - A == {}` |
| `SymmetricDifferenceWithSelf` | `A ^ A == {}` |
| `DisjointSetsIntersection` | `{1,2} & {3,4}` → empty |
| `EqualityOperator` | Same elements → `==` true |
| `InequalityOperator` | Different elements → `!=` true |
| `SubsetCheck` | `is_subset_of({1,2}, {1,2,3}) == true` — NOTE: verify this method exists; if not, skip |
| `SupersetCheck` | Same caveat |

---

## 7. Queue\<T\>

**File:** `test_Queue.cpp`
**Approach:** Basic semantic tests. Low priority, thin wrapper.

### 7.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `DefaultConstructorEmpty` | `size() == 0`, `empty() == true` |
| `PushIncreasesSize` | `push(x)` → `size()` increases |
| `PopReturnsFront` | FIFO: push 1,2,3 → pop returns 1, then 2, then 3 |
| `PopDecreasesSize` | After pop, `size()` decreases |
| `ConstructFromCell` | `Cell<int>{1,2,3}` → Queue pops in order 1,2,3 |
| `CopySemantics` | Independent copy |
| `MoveSemantics` | Source empty after move |

---

## 8. StringList

**File:** `test_StringList.cpp`
**Approach:** All tests are `[valgrind]` candidates. Manual `malloc`/`free` of `char*` arrays.
**NOTE:** StringList is explicitly non-copyable and non-movable (copy/move constructors and assignments are `= delete`). Do NOT write copy/move tests.

### 8.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `ConstructAndPush` | Push 3 strings into capacity-5 list, verify `item(i)` matches |
| `ItemReturnsCString` | `strcmp(item(0), "hello") == 0` |
| `DataPointerValid` | `data()` returns non-null, `data()[i]` matches `item(i)` |
| `PushEmptyString` | `push("")` → `item(i)` returns `""` |
| `PushLongString` | 1000-char string survives round-trip |

### 8.2 Tests `[debug]`

| Test Name | What It Verifies |
|---|---|
| `PushBeyondCapacityThrows` | Push `N+1` strings into capacity-N list → triggers assertion |
| `ItemOutOfBoundsThrows` | `item(mCount)` triggers assertion |

---

## 9. Genome\<B, N\> (Deferred)

**File:** `test_Genome.cpp`
**Approach:** Round-trip functional tests only. Not a container in the strict sense, but exercises `Bitset` and `Vector` in a non-trivial way.

### 9.1 Tests `[semantic]` (implement when Phases 1–3 are complete)

| Test Name | What It Verifies |
|---|---|
| `SetValuesGetValuesRoundTrip` | Encode values, decode → values match within tolerance |
| `LinearScaleRoundTrip` | Non-log parameter encodes/decodes correctly |
| `LogScaleRoundTrip` | Log-scale parameter encodes/decodes correctly |
| `ClampingBehavior` | Out-of-range values are clamped to `[min, max]` |
| `EqualMinMaxIsPassThrough` | When `min == max`, output equals min |
| `RandomizeProducesValidValues` | After `randomize()`, `get_values()` all within `[min, max]` |
| `InheritProducesValidValues` | After `inherit(mom, dad)`, values are within `[min, max]` |
| `KillResetsState` | `kill()` → `is_alive() == false`, `fitness() == BELFEM_REAL_MAX` |
| `FitnessSetAndGet` | `set_fitness(x)` → `fitness() == x` |

---

## 10. Implementation Notes for Claude Code

1. **First action:** Determine the assertion mechanism (see `tests_0_strategy.md` Section 2.2) and document it.
2. **Second action:** Implement `test_Cell.cpp` as the template — get the typed test infrastructure, debug guards, and CMake integration working on this file first.
3. **Third action:** Implement `test_DynamicBitset.cpp` — this is the highest-risk container and will exercise parameterized tests.
4. **Then:** Work down the priority list.
5. **Convention:** When a test name says "Throws", use whatever mechanism was determined in step 1. When it says "verify" or "check", use `EXPECT_EQ` / `EXPECT_TRUE`.
6. **Convention:** Tests marked `[valgrind]` don't need special code — they are normal tests that must pass under Valgrind in CI.
7. **Convention:** If a method listed here does not exist in the actual header (e.g., `is_subset_of` for Set), skip the test and note the absence in a comment.

---

## 11. Codex Audit Checklist

When reviewing Claude Code's test implementation, verify:

- [ ] Every row in the test matrices above has a corresponding `TEST` or `TYPED_TEST`
- [ ] Debug tests are wrapped in `#ifndef NDEBUG`
- [ ] BELFEM naming conventions used (`t` prefix for locals, BELFEM types)
- [ ] No cross-file dependencies between test files
- [ ] No performance/timing assertions
- [ ] Move semantics tests verify source is in valid empty/null state
- [ ] Copy semantics tests verify deep independence (modify copy, check original)
- [ ] Free functions tested separately from member functions
- [ ] Parameterized DynamicBitset tests include block-boundary sizes (63, 64, 65)
