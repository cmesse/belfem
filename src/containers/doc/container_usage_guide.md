# BELFEM Container Classes - Usage Guide {#containers_container_usage_guide}

**Date:** 2026-01-16
**Module:** containers
**Purpose:** Comprehensive guide to BELFEM's custom container classes and their usage
**Revision:** 2026-01-16 - Initial version incorporating community feedback

---

## Overview

The `src/containers` module provides custom wrapper classes around STL containers, designed to offer:

1. **Consistent interface** across BELFEM codebase
2. **Enhanced error checking** with bounds validation in debug mode
3. **Domain-specific functionality** (e.g., genetic algorithms, shift registers, bitset operations)
4. **Performance optimizations** for specific use cases

Most containers are header-only template classes. Exceptions: `DynamicBitset` and `StringList` have separate `.cpp` implementation files.

> **Critical: Debug vs. Release Behavior**
> All `BELFEM_ASSERT` checks are **compiled out** in release builds (`-DNDEBUG`).
> - **Debug:** Bounds checking, descriptive error messages
> - **Release:** Zero overhead, but undefined behavior on violations (mirrors STL)

---

## Common Pitfalls

Before diving into the API, here are the most common mistakes:

### 1. **Forgetting Debug/Release Differences**

```cpp
Cell<int> c(5, 0);          // five zeros ( Cell(5) alone only reserves )
int bad = c(10);  // Debug: Clear assertion failure
                  // Release: Undefined behavior (likely crash or garbage)
```

**Solution:** Test in both debug and release builds. Use static analysis tools.

### 2. **DynamicBitset: Modifying After lock()**

```cpp
DynamicBitset bits(100);
bits.set(5);
bits.lock();      // Compute hash, freeze bitset
bits.set(10);     // ERROR in debug, UB in release
```

**Solution:** Only `lock()` when finalized. Use `unlock()` to modify again.

### 3. **ShiftRegister: revert() Only Works Once**

```cpp
ShiftRegister<int> reg(3);
reg.push(1);
reg.push(2);
reg.revert();     // OK: undo push(2)
reg.revert();     // ERROR: can only revert once per push
```

**Solution:** Only call `revert()` immediately after a failed `push()`.

### 4. **Map: operator[] vs operator()**

```cpp
Map<string, int> m;
int x = m["missing"];   // OK: creates entry with default value (0)
int y = m("missing");   // checked access: throws in a debug build,
                        // aborts in release ( assert.hpp:88-93 )
```

**Solution:** Use `[]` for insert-or-access, `()` for checked read.

### 5. **Genome: Insufficient Bit Resolution**

```cpp
// Genome has no default constructor -- it needs the bounds and the scale flags
Genome<4, 3> genome( mins, maxs, scales );   // only 2^4 = 16 levels per parameter
// Better: Genome<12, 3> for 4096 levels
```

**Solution:** Choose B ≥ 10 for continuous parameters. Use log scaling for wide ranges.

### 6. **Thread Safety**

All containers are **not thread-safe**. Concurrent writes require external synchronization (mutexes).

---

## Memory Management Strategies

Different containers use different memory models:

| Container | Strategy | Notes |
|-----------|----------|-------|
| `Cell`, `Map`, `OrderedMap`, `Set`, `Queue` | STL allocators | RAII, exception-safe |
| `Bitset<N>` | Stack (compile-time size) | Minimal overhead |
| `DynamicBitset` | `malloc`/`free` for `uint64_t` blocks | Manual deallocation in destructor |
| `ShiftRegister<T>` | `malloc`/`free` | +1 extra element for revert backup |
| `StringList` | `malloc`/`free` for `char**` | C API compatibility |
| `Genome` | Embeds `Bitset<N*B>` | Stack-based DNA |

---

## Container Reference

### 1. Cell<T> - Dynamic Array

**File:** `cl_Cell.hpp`
**STL Equivalent:** `std::vector<T>`

**Purpose:** Primary dynamic array container, wrapping `std::vector<T>` with additional functionality.

#### Description

`Cell` is BELFEM's most commonly used container, providing:
- Bounds-checked access via `operator()` in debug mode
- Additional utility functions (`sort`, `unique`, `reverse`, `append`)
- Move semantics support
- Direct access to underlying `std::vector` when needed

#### Key Features

```cpp
// Construction
Cell<int> a;                          // Empty cell
Cell<int> b(10, 0);                  // 10 elements initialized to 0
Cell<int> c = {1, 2, 3, 4};          // Initializer list

// Access
int val = a(0);                       // Bounds-checked access (debug mode)
int& first = a.first();               // First element (bounds-checked in debug)
int& last = a.last();                 // Last element (bounds-checked in debug)

// Modification
a.push(5);                            // Add element (copy)
a.push(std::move(temp));              // Add element (move)
a.emplace(args...);                   // Construct in-place
int popped = a.pop();                 // Remove and return last element

// Memory management
a.set_size(100);                      // Resize
a.reserve(1000);                      // Reserve capacity
a.shrink_to_fit();                    // Free unused memory
a.clear();                            // Clear all elements

// Direct STL access
std::vector<T>& vec = a.vector_data();
T* ptr = a.data();

// Iteration
for (auto& elem : a) { /* ... */ }    // Range-based for loop
```

#### Free Functions

```cpp
// Sorting and uniqueness
sort(myCell);                         // Sort in ascending order
sort(myCell, comparator);             // Sort with custom comparator
unique(myCell);                       // Sort + remove duplicates (in-place)
reverse(myCell);                      // Reverse order

// Combining cells
append(cellA, cellB);                 // Append cellB to cellA (copy)
append_move(cellA, cellB);            // Append cellB to cellA (move, clears cellB)
swap(cellA, cellB);                   // Swap contents (noexcept)
```

#### When to Use

- **Primary choice** for dynamic arrays in BELFEM
- Node lists, element connectivity, general collections
- When you need STL vector functionality with bounds checking

#### Performance Notes

- `operator()` has zero overhead in release builds (NDEBUG)
- Debug builds include comprehensive bounds checking
- `append_move()` is more efficient than `append()` when source is temporary

**See:** `Cell` class and free functions in `cl_Cell.hpp`

---

### 2. Map<Key, Value> - Hash Map

**File:** `cl_Map.hpp`
**STL Equivalent:** `std::unordered_map<Key, Value>`

**Purpose:** Wrapper around `std::unordered_map` with enhanced error reporting.

#### Description

Provides hash-based key-value storage with:
- Descriptive error messages for missing keys (debug mode)
- Consistent interface with other BELFEM containers
- Two access patterns: `[]` (insert-or-access) and `()` (checked access)

#### Key Features

```cpp
// Construction
Map<string, int> ages;
Map<int, Vector<real>> data;

// Insertion
ages["Alice"] = 30;                   // Insert or update
ages.map_data().insert({"Bob", 25});  // STL-style insert via the exposed container

// Access
int age = ages("Alice");              // Checked access (asserts if not found)
int& ref = ages["Charlie"];           // Insert-or-access (creates if missing)

// Queries
bool exists = ages.key_exists("Bob");
size_t n = ages.size();
bool empty = ages.empty();

// Removal
ages.erase_key("Alice");
ages.clear();

// Iteration
for (const auto& [name, age] : ages) {
    // Process key-value pairs
}

// Direct STL access
std::unordered_map<Key, Value>& raw = ages.map_data();
```

#### Access Patterns

**Two operators with different semantics:**

- **`operator[](key)`**: Insert-or-access, creates entry if missing (default-constructed)
- **`operator()(key)`**: Checked access, asserts if key not found

```cpp
Map<string, int> m;
m["new_key"] = 5;      // OK: creates entry
int val = m("new_key"); // OK: key exists
int bad = m("missing"); // ERROR: asserts in debug ( names the key ); BELFEM_ERROR in release ( generic "Key not found in map." )
```

#### When to Use

- ID-to-object mappings (e.g., node ID → node pointer)
- Fast lookups by key (O(1) average case)
- When insertion order doesn't matter

#### Performance Notes

- Hash-based: O(1) average lookup/insert, O(n) worst case
- No ordering guarantees
- Use `OrderedMap` if you need sorted iteration

**See:** `Map` class in `cl_Map.hpp`

---

### 3. OrderedMap<Key, Value> - Sorted Map

**File:** `cl_OrderedMap.hpp`
**STL Equivalent:** `std::map<Key, Value>`

**Purpose:** Wrapper around `std::map` for sorted key-value storage.

#### Description

Similar to `Map`, but maintains keys in sorted order using a red-black tree.

#### Key Differences from Map

- **Iteration order:** Keys iterated in sorted order
- **Performance:** O(log n) lookup/insert (slower than `Map`)
- **Use case:** When sorted iteration is required

```cpp
OrderedMap<int, string> sorted;
sorted[3] = "three";
sorted[1] = "one";
sorted[2] = "two";

// Iterates in key order: 1, 2, 3
for (const auto& [key, val] : sorted) {
    cout << key << ": " << val << endl;
}
```

#### When to Use

- Need sorted iteration by key
- Range queries on keys
- Deterministic iteration order required

**See:** `OrderedMap` class in `cl_OrderedMap.hpp`

---

### 4. Set<Key> - Hash Set

**File:** `cl_Set.hpp`
**STL Equivalent:** `std::unordered_set<Key>`

**Purpose:** Wrapper around `std::unordered_set` with set operations.

#### Description

Unordered collection of unique elements with mathematical set operations.

#### Key Features

```cpp
// Construction
Set<int> primes = {2, 3, 5, 7, 11};
Set<string> names(vec.begin(), vec.end());

// Insertion
primes.insert(13);
primes.emplace(17);

// Queries
bool has5 = primes.contains(5);      // or key_exists(5)
size_t count = primes.count(3);      // Returns 0 or 1
bool empty = primes.empty();

// Removal
primes.erase(2);
primes.clear();

// Set operations
Set<int> a = {1, 2, 3, 4};
Set<int> b = {3, 4, 5, 6};

Set<int> union_ab = a | b;           // {1, 2, 3, 4, 5, 6}
Set<int> intersect = a & b;          // {3, 4}
Set<int> diff = a - b;               // {1, 2}
Set<int> sym_diff = a ^ b;           // {1, 2, 5, 6}

// Subset/superset tests
bool is_sub = a.is_subset_of(b);
bool is_super = a.is_superset_of(b);
```

#### Set Operations Summary

| Operation | Operator | Meaning |
|-----------|----------|---------|
| Union | `a \| b` | Elements in a or b |
| Intersection | `a & b` | Elements in both a and b |
| Difference | `a - b` | Elements in a but not b |
| Symmetric Diff | `a ^ b` | Elements in a or b but not both |

#### When to Use

- Unique element collections
- Fast membership testing (O(1) average)
- Mathematical set operations

**See:** `Set` class in `cl_Set.hpp`

---

### 5. Queue<T> - FIFO Queue

**File:** `cl_Queue.hpp`
**STL Equivalent:** `std::queue<T>`

**Purpose:** Wrapper around `std::queue` for FIFO operations.

#### Description

Simple first-in-first-out queue with conversion from `Cell`.

#### Key Features

```cpp
// Construction
Queue<int> q;
Queue<int> q2(myCell);               // Convert Cell to Queue

// Operations
q.push(10);                          // Add to back
int front = q.pop();                 // Remove and return front element

// Queries
size_t n = q.size();
bool empty = q.empty();
```

#### When to Use

- Breadth-first search algorithms
- Task scheduling
- Any FIFO processing

**See:** `Queue` class in `cl_Queue.hpp`

---

### 6. Bitset<N> - Fixed-Size Bitset

**File:** `cl_Bitset.hpp`
**STL Equivalent:** `std::bitset<N>`

**Purpose:** Wrapper around `std::bitset<N>` for compile-time fixed-size bit arrays.

#### Description

Compile-time sized bitset for flags and boolean arrays.

#### Key Features

```cpp
// Construction (N is compile-time constant)
Bitset<64> flags;

// Setting bits
flags.set(5);                        // Set bit 5 to true
flags.reset(5);                      // Set bit 5 to false
flags.flip(5);                       // Toggle bit 5
flags.reset();                       // Clear all bits

// Querying
bool isSet = flags.test(5);          // Test if bit 5 is set
index_t numSet = flags.count();      // Count true bits
index_t size = flags.size();         // Returns N

// Direct access to std::bitset
std::bitset<N>& raw = flags.data();
```

#### When to Use

- Fixed number of boolean flags known at compile time
- Bit manipulation with known size
- Memory-efficient boolean arrays

#### Limitations

- Size **must** be compile-time constant
- Use `DynamicBitset` for runtime-sized bit arrays

**See:** `Bitset` class in `cl_Bitset.hpp`

---

### 7. DynamicBitset - Runtime-Sized Bitset

**File:** `cl_DynamicBitset.hpp`, `cl_DynamicBitset.cpp`
**STL Equivalent:** None (custom implementation)

**Purpose:** Runtime-resizable bitset with advanced features.

#### Description

Advanced bitset implementation with:
- Runtime-determined size
- Locking mechanism for immutability and fast hashing
- Hash-based comparison (FNV-64 algorithm)
- Bitwise operations (|, &, ^)
- Efficient bit extraction
- Serialization support

#### Key Features

```cpp
// Construction
DynamicBitset bits(1000);            // 1000 bits, initially false

// Setting bits
bits.set(42);                        // Set bit 42
bits.reset(42);                      // Clear bit 42
bits.set(10, true);                  // Set bit 10 to value
bits.flip(5);                        // Toggle bit 5
bits.flip();                         // Toggle all bits
bits.reset();                        // Clear all bits (also unlocks)

// Querying
bool isSet = bits.test(42);
index_t numSet = bits.count();
index_t size = bits.size();
index_t mem = bits.memory();         // Number of uint64_t blocks

// Extract set bit indices
Cell<index_t> indices;
bits.where(indices);                 // returns the set bits, ascending
bits.where(indices, true);           // the flag is accepted for source compatibility
bits.where(indices, false);          // ... and ignored: there is only one algorithm now.
                                     // The scan skips zero data words, but still walks
                                     // every level-2 summary word, so its cost grows with
                                     // the size of the bitset as well as the bit count

// Locking mechanism (for hashing/comparison)
bits.lock();                         // Compute FNV-64 hash, make immutable
size_t hash = bits.hash();           // Get hash (requires lock)
bool locked = bits.is_locked();
bits.unlock();                       // Make writable again (resets hash)

// Comparison (requires both locked)
if (bits1 == bits2) {                // O(1) hash check, then O(n) bit-by-bit
    // ...
}

// Bitwise operations
DynamicBitset result = bits1 | bits2;  // OR
DynamicBitset result = bits1 & bits2;  // AND
DynamicBitset result = bits1 ^ bits2;  // XOR (symmetric difference)

bits1 |= bits2;                      // In-place OR
bits1 &= bits2;                      // In-place AND
bits1 ^= bits2;                      // In-place XOR

// Serialization
string binStr = bits.to_string();    // Binary string "01010..."
string hexStr = bits.to_hex();       // Hex string "A3F2..."
bits.set_from_hex("A3F2");           // Initialize from hex

// Integer conversion (only if size ≤ 8·sizeof(index_t) bits — 32 in a default build, 64 under BELFEM_INT64; larger bitsets hit a BELFEM_ERROR)
index_t val = bits.to_int();         // Convert to integer

// Indexing (optional metadata)
bits.set_index(42);                  // Associate external index
index_t idx = bits.index();          // Retrieve index
```

#### Locking Mechanism

The locking mechanism serves two purposes:

1. **Immutability**: Prevents modification of bitset after finalization
2. **Fast Comparison**: Cached FNV-64 hash enables O(1) inequality checks

```cpp
DynamicBitset a(100), b(100);
// ... set bits ...

a.lock();
b.lock();

// Fast comparison using cached hashes
if (a == b) {
    // Hash check first (O(1)), then bit-by-bit if hashes match (O(n))
}

// This would trigger assertion in debug:
// a.set(5);  // ERROR: Can't modify locked bitset
```

> **Note:** `reset()` (clear all bits) automatically unlocks the bitset.

#### Serialization

DynamicBitset supports multiple serialization formats:

```cpp
DynamicBitset bits(128);
// ... set bits ...

// Hex (human-readable)
string hex = bits.to_hex();          // "00A3F2..."
bits.set_from_hex(hex);              // Restore

// Binary string (visual debugging)
string bin = bits.to_string();       // "0010101..."

// Raw bytes (not shown in API, but available via data())
const uint64_t* raw = bits.data();   // Direct access to blocks
```

#### Performance Notes

- **`where()` method**: one algorithm; walks the two-level summary bitmaps, so it skips every zero data word (cost: all level-2 words + touched words + set bits). The `aAssumeSparse` flag is accepted and ignored.
- **Bitwise operations**: Optimized with pointer arithmetic, vectorizable
- **Hash comparison**: O(1) for inequality, O(n) for equality confirmation
- **Storage**: Uses `uint64_t` blocks (8 bytes each), efficient for large bit arrays

#### When to Use

- Runtime-determined boolean arrays
- Cohomology/topology computations (set operations on facets)
- Hash-based bitset comparisons (use in maps/sets)
- Need to extract indices of set bits efficiently

#### Lesson Learned: The where_sparse() Optimization Paradox

During development of the `where_sparse()` routine, six different AI models suggested "obvious" optimizations:
- Hoisting bounds checks out of the inner loop
- Loop splitting (separate fast loop for full blocks vs. partial blocks)
- Aggressive masking to eliminate inner-loop logic

**The silicon reality:** Micro-benchmarking revealed the **original, unoptimized loop was 82% faster** than the theoretically "superior" version.

**Why the simple loop won:**

1. **Branch Predictability**: Because unused bits in the final block are **zeroed by construction** (class invariant), the branch `if (tIndex < mNumberOfBits)` is true ~100% of the time. Modern branch predictors handle this with virtually zero clock-cycle cost.

2. **Instruction Cache Density**: The simple loop resulted in smaller binary code. The "optimized" loop-splitting approach doubled the code size, increasing I-cache pressure and preventing effective compiler auto-vectorization.

3. **Low-Latency Invariants**: Relying on class invariants (high bits always zero) is cheaper than enforcing them at runtime via masking.

> **Silicon Guardrail Rule:** Profile before you perfect. Modern silicon is smarter than most algorithms. Keep hot loops simple and compact—more code is rarely faster code.

**See also:** `doc/lessons_learned.md` for the full case study and additional performance lessons.

**See:** `DynamicBitset` class in `cl_DynamicBitset.hpp` and `.cpp`

---

### 8. ShiftRegister<T> - Fixed-Capacity FIFO with Revert

**File:** `cl_ShiftRegister.hpp`
**STL Equivalent:** None (custom circular buffer)

**Purpose:** Fixed-capacity circular buffer with newest-first ordering and one-step revert capability.

#### Description

A specialized container that maintains the N most recent values, with:
- Newest value at index 0, oldest at index N-1
- Fixed capacity determined at construction
- **One-step revert capability** (undo last push)
- Efficient memory management with `malloc`/`free`

#### Key Features

```cpp
// Construction
ShiftRegister<double> history(5);              // Capacity 5, empty
ShiftRegister<double> temps(5, 20.0);          // Capacity 5, filled with 20.0
ShiftRegister<int> recent = {1, 2, 3, 4, 5};   // From initializer list

// Adding values
history.push(100.5);                           // Add newest value
// Contents: [100.5] (size=1)
history.push(200.3);
// Contents: [200.3, 100.5] (size=2)

// Access (0 = newest, N-1 = oldest)
double newest = history(0);                    // Most recent
double oldest = history(history.size()-1);     // Oldest

// Revert last push (CRITICAL: only works once per push)
history.push(300.0);                           // Now [300.0, 200.3, 100.5]
history.revert();                              // Back to [200.3, 100.5]
// history.revert();                           // ERROR: can't revert twice

// Queries
size_t current = history.size();               // Current number of values
size_t max = history.capacity();               // Maximum capacity
bool isEmpty = history.empty();
bool isFull = history.full();

// Modification
history.clear();                               // Remove all values
history.fill(0.0);                             // Fill all current values with 0.0

// Memory management
history.reserve(10);                           // Change capacity (clears data)
history.free();                                // Release memory explicitly
```

#### Revert Capability

> **Critical:** `revert()` can only be called **once** after each `push()`.

The shift register allocates one extra backup slot internally:

```cpp
ShiftRegister<int> reg(3);          // Capacity 3, allocates 4 slots
reg.push(1); reg.push(2); reg.push(3);
// State: [3, 2, 1] (full)

reg.push(4);
// State: [4, 3, 2]  (1 is in backup slot)

reg.revert();
// State: [3, 2, 1]  (restored from backup)

reg.push(5);                        // Must push before next revert
reg.revert();                       // OK
```

**Use case:** Time integration where you may need to reject a timestep:

```cpp
ShiftRegister<Vector<real>> history(3);
history.push(initialState);

for (size_t step = 0; step < numSteps; ++step) {
    Vector<real> newState = time_integrate(history(0), dt);
    history.push(newState);

    if (!check_convergence(newState)) {
        history.revert();           // Undo failed step
        dt *= 0.5;                  // Reduce timestep
        continue;                   // Retry
    }
}
```

#### Iteration

```cpp
ShiftRegister<double> data(10, 1.0);
for (double val : data) {
    // Iterates from newest (index 0) to oldest
}
```

#### When to Use

- Time-stepping algorithms (storing previous timestep values)
- Moving window computations (moving averages, gradients)
- Iterative solvers with history dependence (BDF methods)
- Undo functionality for single operations

#### Performance Notes

- Uses `std::malloc`/`free` for memory (not `new`/`delete`)
- Efficient `push()`: O(capacity) using `std::move_backward`
- Fixed capacity: no dynamic resizing overhead
- Revert is O(capacity) but only available immediately after `push()`

**See:** `ShiftRegister` class in `cl_ShiftRegister.hpp`

---

### 9. Genome<B, N> - Genetic Algorithm Encoding

**File:** `cl_Genome.hpp`
**STL Equivalent:** None (custom genetic algorithm encoding)

**Purpose:** Encode real-valued parameters as bit strings for genetic algorithms.

#### Description

Template class for genetic algorithm optimization:
- **B**: Bits per parameter (resolution = 2^B levels)
- **N**: Number of parameters
- Supports linear and logarithmic parameter scaling
- Implements crossover, mutation, and fitness tracking

> **Tip:** Use B ≥ 10 for continuous parameters (1024 levels). Use log scaling for parameters spanning orders of magnitude.

#### Key Features

```cpp
// Setup parameter bounds and scaling
Vector<real> minVals = {0.1, 1.0, 0.001};
Vector<real> maxVals = {10.0, 100.0, 1.0};
Bitset<3> logScale;                           // false = linear, true = log
logScale.set(0);                              // First parameter uses log scale

// Create genome (8 bits per parameter = 256 levels, 3 parameters)
Genome<8, 3> individual(minVals, maxVals, logScale);

// Initialize with specific values
Vector<real> params = {0.5, 50.0, 0.1};
individual.set_values(params);

// Or randomize (Gaussian distribution centered in range)
individual.randomize();

// Extract parameter values
Vector<real> decoded(3);
individual.get_values(decoded);

// Genetic operations
Genome<8, 3> mom(minVals, maxVals, logScale);
Genome<8, 3> dad(minVals, maxVals, logScale);
Genome<8, 3> child(minVals, maxVals, logScale);

child.inherit(&mom, &dad);                    // Single-point crossover + 1-bit mutation

// Fitness tracking
individual.set_fitness(0.042);
real fitness = individual.fitness();
bool alive = individual.is_alive();           // fitness != REAL_MAX
individual.kill();                            // Set fitness to REAL_MAX
```

#### Genetic Algorithm Workflow

```cpp
// Typical optimization loop
const size_t BITS = 12;                       // 4096 levels per parameter
const size_t PARAMS = 5;                      // 5 parameters to optimize
const size_t POP_SIZE = 100;

// NOTE: Cell( n ) only RESERVES -- it does not size the container, so a
// range-for over it here would iterate over nothing. Size it with the
// two-argument constructor, or push each element.
Cell<Genome<BITS, PARAMS>*> population(POP_SIZE, nullptr);

// Initialize population
for (auto& genome : population) {
    genome = new Genome<BITS, PARAMS>(mins, maxs, scales);
    genome->randomize();
}

// Evolution loop
for (size_t gen = 0; gen < MAX_GENERATIONS; ++gen) {
    // Evaluate fitness
    for (auto& genome : population) {
        Vector<real> params(PARAMS);
        genome->get_values(params);
        real fitness = objective_function(params);
        genome->set_fitness(fitness);
    }

    // Sort by fitness
    opGenomeSort<BITS, PARAMS> sorter;
    sort(population, sorter);

    // Create next generation (elitism + breeding)
    Cell<Genome<BITS, PARAMS>*> nextGen(POP_SIZE, nullptr);   // size it; Cell(n) only reserves

    // Keep top 10% (elitism)
    for (size_t i = 0; i < POP_SIZE/10; ++i) {
        nextGen(i) = population(i);
    }

    // Breed the rest
    for (size_t i = POP_SIZE/10; i < POP_SIZE; ++i) {
        size_t momIdx = tournament_select(population);
        size_t dadIdx = tournament_select(population);
        nextGen(i)->inherit(population(momIdx), population(dadIdx));
    }

    population = nextGen;
}
```

#### Parameter Encoding

- **Linear scale**: `value = min + (max - min) × (bits / (2^B - 1))`
- **Log scale**: `value = exp(log(min) + (log(max) - log(min)) × (bits / (2^B - 1)))`

Use log scale for parameters spanning orders of magnitude (e.g., conductivity: 0.001 to 100).

#### When to Use

- Parameter optimization via genetic algorithms
- Multi-objective optimization problems
- Non-gradient-based optimization
- Discrete parameter search

**See:** `Genome` class in `cl_Genome.hpp`

---

### 10. StringList - C-String Array for I/O

**File:** `cl_StringList.hpp`, `cl_StringList.cpp`
**STL Equivalent:** None (C-compatible string array)

**Purpose:** Fixed-size C-string array for interfacing with C libraries (e.g., Exodus).

#### Description

Low-level string container using `char**` for compatibility with C APIs.

#### Key Features

```cpp
// Construction
StringList names(10);                // Pre-allocate for 10 strings

// Adding strings
names.push("node_block_1");
names.push("element_block_2");

// Access
const char* first = names.item(0);   // Bounds-checked access

// Raw pointer (for C APIs)
char** raw = names.data();
// Example: exodus_write_names(exoid, raw);
```

#### When to Use

- **Only** when interfacing with C libraries requiring `char**`
- Exodus file I/O
- Otherwise use `Cell<string>`

#### Limitations

- Fixed size at construction
- Manual memory management (RAII via destructor)
- Not recommended for general use

**See:** `StringList` class in `cl_StringList.hpp` and `.cpp`

---

## Container Selection Guide

### Quick Reference

| Need | Use |
|------|-----|
| Dynamic array of values | `Cell<T>` |
| Fast key-value lookup | `Map<Key, Value>` |
| Sorted key-value pairs | `OrderedMap<Key, Value>` |
| Unique element collection | `Set<T>` |
| FIFO queue | `Queue<T>` |
| Fixed-size boolean flags (compile-time) | `Bitset<N>` |
| Dynamic boolean array with operations | `DynamicBitset` |
| Recent value history with rollback | `ShiftRegister<T>` |
| Genetic algorithm parameters | `Genome<B, N>` |
| C API string interface | `StringList` |

### Performance Characteristics

| Container | Access | Insert | Search | Memory Strategy |
|-----------|--------|--------|--------|-----------------|
| `Cell` | O(1) | O(1) amortized | O(n) | STL allocator |
| `Map` | O(1) avg, O(n) worst | O(1) avg | O(1) avg | Hash table |
| `OrderedMap` | O(log n) | O(log n) | O(log n) | Red-black tree |
| `Set` | O(1) avg | O(1) avg | O(1) avg | Hash table |
| `Queue` | O(1) front | O(1) | N/A | STL deque |
| `Bitset<N>` | O(1) | O(1) | O(1) | Stack (compile-time) |
| `DynamicBitset` | O(1) | O(1) | O(1) | `malloc` uint64 blocks |
| `ShiftRegister` | O(1) | O(N) | N/A | `malloc` +1 backup |
| `Genome` | N/A | N/A | N/A | Embedded `Bitset<N×B>` |

---

## Common Patterns

### Pattern 1: Building and Sorting a Cell

```cpp
Cell<index_t> nodeIDs;
nodeIDs.reserve(estimatedSize);      // Avoid reallocations

for (auto* node : nodes) {
    nodeIDs.push(node->id());
}

sort(nodeIDs);                       // Sort in place
unique(nodeIDs);                     // Remove duplicates
```

### Pattern 2: Map with Default Values

```cpp
Map<index_t, real> values;

// Get-or-insert pattern
real& val = values[nodeID];          // Creates entry if missing (default = 0.0)
val += contribution;                 // Safe even if new entry
```

### Pattern 3: Set Operations for Filtering

```cpp
Set<index_t> boundaryNodes = {...};
Set<index_t> activeNodes = {...};

// Find nodes that are boundary AND active
Set<index_t> boundaryActive = boundaryNodes & activeNodes;

// Find active nodes that are NOT on boundary
Set<index_t> interior = activeNodes - boundaryNodes;
```

### Pattern 4: DynamicBitset for Mesh Topology

```cpp
// Track which elements touch a node
DynamicBitset elementFlags(numElements);

for (index_t e = 0; e < numElements; ++e) {
    if (element_touches_node(e, targetNode)) {
        elementFlags.set(e);
    }
}

// Extract element indices efficiently
Cell<index_t> touchingElements;
elementFlags.lock();                 // Finalize bitset (compute hash)
elementFlags.where(touchingElements);        // Extract set indices
```

### Pattern 5: ShiftRegister for Adaptive Time Integration

```cpp
// Store last 3 timestep values for BDF-3
ShiftRegister<Vector<real>> history(3, initialState);   // capacity 3, all slots primed:
                                                        // history(1), history(2) are valid from the first step

for (size_t step = 0; step < numSteps; ++step) {
    Vector<real> newState = bdf3_step(
        history(0),                  // Current state (t^n)
        history(1),                  // Previous state (t^{n-1})
        history(2),                  // Two steps back (t^{n-2})
        dt
    );

    history.push(newState);          // Tentatively accept

    if (!converged(newState, tol)) {
        history.revert();            // Undo failed step
        dt *= 0.5;                   // Reduce timestep
        continue;                    // Retry
    }

    dt = adaptive_dt(newState);      // Update timestep
}
```

---

## Design Philosophy

### Why Wrappers?

BELFEM wraps STL containers for several reasons:

1. **Consistent Interface**: All containers use similar naming (e.g., `set_size()` vs. `resize()`)
2. **Enhanced Debugging**: Bounds checking with descriptive error messages in debug builds
3. **Zero Overhead**: Release builds compile to identical STL performance
4. **Domain Extensions**: Additional functionality for FEM/physics (e.g., `unique()`, set operations)
5. **Future Flexibility**: Can switch implementations without changing client code

### Debug vs. Release Behavior

Most containers have different behavior in debug and release builds:

```cpp
// Debug build (BELFEM_ASSERT active)
Cell<int> c(5, 0);          // five zeros ( Cell(5) alone only reserves )
int bad = c(10);  // ASSERTION FAILURE with clear message:
                  // "Cell index out of bounds: 10 (expect < 5)"

// Release build (NDEBUG defined)
int bad = c(10);  // Undefined behavior (same as std::vector)
```

This provides safety during development while maintaining full performance in production.

> **Recommendation:** Always test in both debug and release modes. Use static analysis tools (e.g., clang-tidy, valgrind) to catch UB in release builds.

### Memory Management

See "Memory Management Strategies" table at the top for details.

---

## Best Practices

### 1. Prefer Cell Over std::vector

```cpp
// Good - BELFEM standard
Cell<Node*> nodes;

// Avoid - use Cell for consistency
std::vector<Node*> nodes;
```

### 2. Reserve Memory When Size is Known

```cpp
Cell<index_t> ids;
ids.reserve(mesh->number_of_nodes());  // Avoid reallocations
for (auto* node : mesh->nodes()) {
    ids.push(node->id());
}
```

### 3. Use Appropriate Map Type

```cpp
// Need fast lookup, don't care about order
Map<id_t, Element*> elementMap;       // O(1) access

// Need sorted iteration
OrderedMap<id_t, Element*> sorted;    // O(log n) access, sorted iteration
```

### 4. Lock DynamicBitsets Before Comparison

```cpp
DynamicBitset a(100), b(100);
// ... set bits ...

a.lock();
b.lock();

if (a == b) {  // Fast hash comparison
    // ...
}
```

> **Caution:** Don't forget to `unlock()` if you need to modify again.

### 5. DynamicBitset where() Has One Strategy

```cpp
DynamicBitset flags(10000);
Cell<index_t> indices;
flags.where(indices);                 // the second argument, if given, is ignored
```

### 6. Understand ShiftRegister Revert Limitation

```cpp
ShiftRegister<double> reg(5);
reg.push(1.0);
reg.push(2.0);
reg.revert();     // OK
reg.push(3.0);    // Must push before next revert
reg.revert();     // OK
```

---

## Examples from BELFEM Codebase

### Example 1: Mesh Node Storage

```cpp
// src/mesh/cl_Mesh.hpp
Cell<Node*> mNodes;                   // All nodes in mesh

// Efficient construction
mNodes.reserve(nodeCount);
for (index_t i = 0; i < nodeCount; ++i) {
    mNodes.push(new Node(i, coords));
}
```

### Example 2: Element-to-Node Connectivity

```cpp
// src/fem/kernel/cl_Element.hpp
Cell<Node*> mNodes;                   // Nodes of this element

// Access element node
Node* node = mNodes(localIndex);      // Bounds-checked
```

### Example 3: Cohomology Facet Operations

```cpp
// src/homology/... (hypothetical)
DynamicBitset facet1(numNodes);
DynamicBitset facet2(numNodes);

// Set bits for nodes in each facet
// ...

// Compute boundary operator (symmetric difference)
DynamicBitset boundary = facet1 ^ facet2;

// Lock and hash for storage in map
boundary.lock();
Map<size_t, DynamicBitset*> facetMap;
facetMap[boundary.hash()] = &boundary;
```

---

## Thread Safety and MPI

**Important:** BELFEM containers are **not thread-safe** by default.

- **Reading**: Multiple threads can safely read from const containers
- **Writing**: Use external synchronization (mutexes, atomics) for concurrent writes
- **MPI**: Each rank has separate container instances (no shared memory)
  - Use `comm` module for synchronization (e.g., all-gather for Cells)

Example with OpenMP:

```cpp
Cell<int> shared(100);

#pragma omp parallel for
for (int i = 0; i < 100; ++i) {
    // Read: safe
    int val = shared(i);

    // Write: NOT safe without synchronization
    #pragma omp critical
    shared(i) = compute(val);
}
```

---

## Related Modules

- **linalg**: `Vector`, `Matrix` containers for numerical linear algebra
- **mesh**: Uses `Cell` extensively for node/element storage
- **fem**: Uses containers throughout for DOF management, assembly
- **comm**: MPI communication helpers for container synchronization
- **homology**: Uses `DynamicBitset` for topology operations

---

## See Also

- [Linalg Module](../../linalg/doc/README.md) - Linear algebra containers (Vector, Matrix)
- [Mesh Module](../../mesh/doc/README.md) - Mesh data structures using containers
- [Homology Module](../../homology/doc/README.md) - Cohomology algorithms using bitsets
- `CLAUDE.md` (repository root) - Documentation organization guidelines
- [Documentation Guidelines](../../../doc/documentation_guidelines.md) - How to document code
- C++ Documentation: `make doc` (Doxygen)

---

**Revision History:**
- 2026-01-16: Initial version incorporating feedback from multiple AI reviewers (Grok, ChatGPT, CBorg, Gemini, Opus)
