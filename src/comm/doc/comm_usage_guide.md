# BELFEM Communication Module - Usage Guide {#comm_comm_usage_guide}

**Date:** 2026-01-16
**Module:** comm
**Purpose:** Comprehensive guide to BELFEM's MPI communication abstraction layer

**Revision History:**
- 2026-01-16: Initial version
- 2026-01-16: Critical corrections from external review (blocking semantics, thread safety, process limits)

---

## Communication Contracts

**Read this first:** The BELFEM communication module has strict semantic contracts:

- **All operations are blocking and synchronous** at the BELFEM API level, even when implemented using non-blocking MPI internally
- **Message matching is deterministic and explicit** (rank + automatically generated tag)
- **No function is thread-safe** unless explicitly stated
- **Non-MPI builds compile the comm calls to no-ops** — `broadcast`/`send`/`receive` are semantically neutral in serial, `collect` is not (see Non-MPI Builds)
- **Collective operations require consistent participation** by all ranks with matching parameters

> **Rule:** Calling any comm function before `gComm.init()` or after `gComm.finalize()` is undefined behavior.

---

## Overview

The `src/comm` module provides a **compile-time abstraction** over MPI (Message Passing Interface) for parallel computing in BELFEM. It wraps MPI operations with type-safe C++ templates and supports both **MPI builds** and **non-MPI builds** (single-process fallback).

### Key Features

1. **Type-Safe Communication**: Template-based wrappers that automatically deduce MPI datatypes
2. **Dual-Mode Operation**: Works with or without MPI (controlled by `BELFEM_MPI` flag)
3. **Container Integration**: Direct support for `Cell<T>`, `Vector<T>`, `Matrix<T>`
4. **Automatic Chunking**: `send`, `receive` and `share` split large messages into chunks of `gMaxCommChunkLength` (65 536 **elements of `T`**, not bytes). `broadcast` does *not* chunk: scalars and raw arrays go through a blocking `MPI_Bcast`, containers through two `MPI_Ibcast` calls — one for the size, one for the whole payload
5. **Global Communicator**: Single `gComm` instance managed by the framework
6. **Zero-Overhead**: Non-MPI builds compile to no-ops with no runtime cost

### When to Use

- Distributed finite element assembly across multiple processes
- Parallel linear algebra operations
- Domain decomposition and load balancing
- Distributed mesh I/O and partitioning

> **Important Limitations:**
> - **Thread Safety:** The comm module is **not** thread-safe. However, when compiled with `BELFEM_STRUMPACK`, MPI is initialized with `MPI_THREAD_MULTIPLE`. Concurrent calls from multiple threads still require external synchronization.
> - **Process Count:** Maximum ~32,768 processes at a typical `MPI_TAG_UB`, set by the `comm_tag()` algorithm (see Process Count Limitation section).
> - **Message Size:** Individual message chunks limited to `int` range (~2GB) even on 64-bit builds.

---

## Common Pitfalls

### 1. **Deadlocks from Send/Receive Mismatches**

```cpp
// WRONG: Both processes waiting to receive
if (comm_rank() == 0) {
    receive(data, 1);  // Waits for proc 1
}
if (comm_rank() == 1) {
    receive(data, 0);  // Waits for proc 0 → DEADLOCK
}

// CORRECT: One sends, one receives
if (comm_rank() == 0) {
    send(data, 1);
} else if (comm_rank() == 1) {
    receive(data, 0);
}
```

### 2. **Forgetting MPI Initialization**

```cpp
// WRONG: Using comm before init
int main(int argc, char** argv) {
    proc_t rank = comm_rank();  // Undefined behavior!
    gComm.init(argc, argv);
}

// CORRECT: Init first
int main(int argc, char** argv) {
    gComm.init(argc, argv);
    proc_t rank = comm_rank();  // OK
}
```

### 3. **Size Mismatches in Distribute/Collect**

```cpp
// WRONG: Cell size doesn't match comm_size()
Cell<int> data(5, 0);  // But we have 4 processes
distribute(data);      // Assertion failure!

// ALSO WRONG, and less obvious: the one-argument Cell constructor only
// RESERVES. size() is still 0, so this trips the same assertion.
Cell<int> data(comm_size());

// CORRECT: size it, do not just reserve
Cell<int> data(comm_size(), 0);
distribute(data);      // size OK ( still needs its matching receiver )
```

### 4. **Broadcasting Uninitialized Data**

`broadcast()` resizes the receiving container itself, so this is **not** something you need to
guard against:

```cpp
Vector<real> v;             // empty on the non-root ranks
if (comm_rank() == 0) {
    v.set_size(100);
    // fill v...
}
broadcast(v, 0);            // non-root ranks are resized and filled by the call
```

The real hazard is the opposite one: reaching for `broadcast` on a payload large enough that its
unchunked container transfer fails. For large or variable-size data use the `share`/`receive`
pair documented below, which chunks.

### 5. **Not Calling comm_barrier() Before Timing**

```cpp
// WRONG: Ranks finish at different times
Timer timer;
expensive_operation();
real time = timer.stop();  // Different on each rank!

// CORRECT: Synchronize before timing
comm_barrier();
Timer timer;
expensive_operation();
comm_barrier();
real time = timer.stop();  // Consistent measurement
```

---

## Global Communicator

### `gComm` - Global Communicator Instance

**File:** Defined in `cl_Communicator.hpp`, declared `extern`

The framework provides a single global `Communicator` instance that manages all MPI state.

#### Initialization

```cpp
int main(int argc, char** argv) {
    // Initialize MPI and parse arguments
    gComm.init(argc, argv);

    // ... application code ...

    // Finalize MPI before exit
    return gComm.finalize();
}
```

#### Core Methods

```cpp
// Get process rank (0 to size-1)
proc_t rank = gComm.rank();

// Get total number of processes
proc_t size = gComm.size();

// Get MPI_Comm handle (MPI builds only)
MPI_Comm& world = gComm.world();

// Get executable path
const std::string& path = gComm.exec_path();

// Get working directory
const std::string& dir = gComm.workdir();

// Access command-line arguments
Cell<std::string>& args = gComm.arguments();

// Random number generator (rank-specific seed)
std::mt19937& rng = gComm.random();

// Maximum MPI tag value
int max_tag = gComm.max_tag();
```

**See:** `cl_Communicator.hpp:53`

---

## Initialization and Lifecycle

### Initialization Sequence

BELFEM's MPI initialization involves several automatic steps:

```
┌─────────────────────────────────────────────────────┐
│ Pre-main: setup_mpi_binding() [GCC constructor]    │
│   ├─ setenv("OMPI_MCA_hwloc_base_binding_policy", "none")│
│   ├─ setenv("PRTE_MCA_hwloc_default_binding_policy", "none")│
│   └─ setenv("HYDRA_BINDING", "none")               │
└─────────────────────────────────────────────────────┘
                         ↓
┌─────────────────────────────────────────────────────┐
│ gComm.init(argc, argv)                              │
│   ├─ MPI_Init() or MPI_Init_thread()               │
│   │  (MPI_THREAD_MULTIPLE if BELFEM_STRUMPACK)     │
│   ├─ PetscInitialize() [if BELFEM_PETSC]           │
│   ├─ Type checks (PetscInt==int, PetscReal==real)  │
│   ├─ MPI_Comm_rank/size()                          │
│   ├─ Seed rank-specific random generator            │
│   └─ set_globals()                                  │
└─────────────────────────────────────────────────────┘
                         ↓
              [ Application code ]
                         ↓
┌─────────────────────────────────────────────────────┐
│ gComm.finalize()                                    │
│   ├─ Free CommunicationObjects                     │
│   ├─ MPI_Barrier()                                  │
│   ├─ PetscFinalize() [if BELFEM_PETSC]             │
│   └─ MPI_Finalize()                                │
└─────────────────────────────────────────────────────┘
```

> **Pre-main Environment Setup:** BELFEM automatically disables MPI process binding at load time (via GCC `__attribute__((constructor))`) to prevent conflicts with internal parallelism. This affects Open MPI, PRRTE, and MPICH Hydra. If you need specific binding, set these environment variables explicitly before launching mpirun.

### Lifecycle Rules

> **Critical:**
> - `gComm.init()` must be called **exactly once**, before any comm function
> - `gComm.finalize()` must be the **last** MPI-related call in `main()`
> - Calling `comm_rank()` or any comm function before `init()` is **undefined behavior**
> - Calling any comm function after `finalize()` is **undefined behavior**

### Process Count Limitation

Due to the `comm_tag()` algorithm, BELFEM has a maximum process limit based on your MPI implementation's `MPI_TAG_UB` value:

```
max_procs ≈ 0.5 × (√(2 × MPI_TAG_UB + 1) + 1)
```

For a typical `MPI_TAG_UB ≈ 2³¹` this gives approximately **32,768 processes**. Exceeding it is a `BELFEM_ERROR` at initialization (`cl_Communicator.cpp:219-224`). (An earlier revision of this page said 46,340, which is `√(2³¹)` — the square root of the tag ceiling rather than the formula above applied to it.)

**Tag Generation Formula:** For source rank `s` and target rank `t`, the tag is computed as:
```
tag = 2 × (max(s,t) × comm_size() + min(s,t)) mod max_tag
```

This ensures unique tags for each process pair but limits the maximum communicator size.

### PETSc Integration

When `BELFEM_PETSC` is enabled:
- MPI is initialized **before** PETSc
- Runtime checks verify `PetscInt` matches `int` and `PetscReal` matches `real`
- Both `PetscFinalize()` and `MPI_Finalize()` are called during cleanup
- Do not manually call `PetscInitialize()` — it's handled automatically

### CommunicationObject Base Class

For objects that manage MPI resources (e.g., persistent communicators, windows), inherit from `CommunicationObject`:

```cpp
class MyMPIResource : public CommunicationObject {
public:
    void free() override {
        // Release MPI resources (windows, comms, etc.)
        MPI_Win_free(&my_window);

        // Unregister from gComm
        CommunicationObject::free();
    }
};
```

Objects are automatically registered with `gComm.objects()` upon construction and cleaned up when `gComm.finalize()` is called. This prevents MPI resource leaks in parallel environments.

---

## Type System

### `comm_type<T>()` - MPI Datatype Mapping

**File:** `commtypes.hpp`

Template function that returns the MPI datatype for C++ type `T`.

#### Supported Types

| C++ Type | MPI Type |
|----------|----------|
| `char` | MPI_CHAR |
| `signed char` | MPI_SIGNED_CHAR |
| `unsigned char` | MPI_UNSIGNED_CHAR |
| `short int` | MPI_SHORT |
| `unsigned short` | MPI_UNSIGNED_SHORT |
| `int` | MPI_INT |
| `unsigned int` | MPI_UNSIGNED |
| `long int` | MPI_LONG |
| `unsigned long` | MPI_UNSIGNED_LONG |
| `long long int` | MPI_LONG_LONG |
| `unsigned long long` | MPI_UNSIGNED_LONG_LONG |
| `float` | MPI_FLOAT |
| `double` | MPI_DOUBLE |
| `long double` | MPI_LONG_DOUBLE |
| `bool` | MPI_CXX_BOOL |
| `std::complex<float>` | MPI_CXX_FLOAT_COMPLEX |
| `std::complex<double>` | MPI_CXX_DOUBLE_COMPLEX |
| `std::complex<long double>` | MPI_CXX_LONG_DOUBLE_COMPLEX |

#### Usage

```cpp
// Automatically deduced in templates
template<typename T>
void send_scalar(const T& value, proc_t target) {
    comm_t type = comm_type<T>();  // MPI_DOUBLE for T=double
    MPI_Send(&value, 1, type, target, 0, MPI_COMM_WORLD);
}
```

**Note:** Non-MPI builds define `comm_t` as `int` and `comm_type<T>()` returns 0.

> **Limitation:** `comm_type<T>()` supports only **trivially copyable, contiguous types** (primitives and std::complex). User-defined structs require manual MPI datatype definitions via `MPI_Type_create_struct()`.

---

## Utility Functions

### `comm_size()` - Number of Processes

```cpp
proc_t size = comm_size();  // Returns gComm.size()
```

**Returns:** Total number of MPI processes (1 in non-MPI builds)

**File:** `commtools.hpp:53`

---

### `comm_rank()` - Current Process Rank

```cpp
proc_t rank = comm_rank();  // Returns gComm.rank()
```

**Returns:** Rank of current process (0 to size-1), always 0 in non-MPI builds

**File:** `commtools.hpp:62`

---

### `comm_barrier()` - Synchronization Barrier

```cpp
comm_barrier();  // All processes wait here until all arrive
```

**Use case:** Synchronize before timing, I/O, or phase transitions

**File:** `commtools.hpp:80`

---

### `comm_tag(source, target)` - Generate Unique Tag

```cpp
int tag = comm_tag(0, 5);  // Unique tag for communication from rank 0 to 5
```

**Use case:** Internal use by send/receive to avoid message confusion

**File:** `commtools.hpp:124`

---

### `comm_drain_check(label)` - Debug Tripwire for the Fabric Ordering Contract

```cpp
comm_drain_check( "DofManager::postprocess" );  // at an exchange boundary
```

Each rank pair shares exactly two tags (see the ordering contract on
`comm_tag()` in `commtools.cpp`). If a send's matching receive is skipped, the
failure does not surface where the send was issued; it silently poisons a later
receive on the same tag. This check probes both fabric tags for every pair that
includes the calling rank. If it finds a queued message, it raises
`BELFEM_ERROR` on **every** rank, and the detecting ranks report the boundary
label, source, tag, and byte count. An allreduce combines the local verdicts, so
no rank is left waiting in a barrier and the debug throw policy stays consistent
across ranks.

**Use case:** call this at exchange boundaries — points that every rank reaches
with no exchange in flight. Current call sites: `DofManager::solve`,
`DofManager::solve_from_residual`, `DofManager::postprocess`,
`SolverData::collect_matrices` (entry), and `Controller::finalize` (exit).

**Contract:**

- Collective wherever assertions are active. A call site reachable by only a
  subset of ranks is itself a deadlock, just like a misplaced `comm_barrier()`.
- No-op in release builds and serial runs; it costs nothing outside debug.
- Best-effort in two ways: a stray still in flight when the sweep runs can
  escape it, and a stray carrying a tag outside the pair's two fabric tags is
  invisible to it. A clean pass is therefore not proof of a clean fabric,
  although a reported stray is always a real one.

**File:** `commtools.hpp:113`

---

### `comm_split(length)` - Split Message into Chunks

```cpp
Cell<int> chunks = comm_split(100000);  // Splits into chunks of ≤ 65 536 elements
```

**Returns:** Cell containing chunk sizes for reliable large message transmission

**Chunk size:** Defined by `gMaxCommChunkLength = 64*1024` elements of `T` (not bytes: 512 KiB for `double`)

> **Rationale:** Chunking avoids implementation-dependent MPI limits on message size and improves robustness on older interconnects and debug builds. The 65 536-element limit balances message overhead with reliability.

**File:** `commtools.hpp:134`

---

## Broadcast Operations

### `broadcast(data, root)` - Broadcast to All Processes

Sends data from `root` process to all other processes. Non-root processes receive and resize automatically.

> **Broadcast Rule:** Only the `root` rank's data is used. All other ranks' input values are **ignored and overwritten**. Non-root ranks do not need to pre-allocate containers — they are automatically resized.

#### Scalar Broadcast

```cpp
int value = 42;
if (comm_rank() == 0) {
    value = 100;  // Only root sets value
}
broadcast(value, 0);  // All ranks now have value = 100
```

**Signature:** `void broadcast(T& aMessage, proc_t aRoot = 0)`

**Constraints:** `T` must be arithmetic type

---

#### Raw Array Broadcast

```cpp
double data[10];
if (comm_rank() == 0) {
    // Initialize data...
}
broadcast(data, 0, 10);  // Broadcast 10 elements
```

**Signature:** `void broadcast(T* aMessage, proc_t aRoot, proc_t aLength)`

---

#### Cell Broadcast

```cpp
Cell<int> cell;
if (comm_rank() == 0) {
    cell = {1, 2, 3, 4, 5};
}
broadcast(cell, 0);  // Non-root ranks automatically resized to 5 elements
```

**Signature:** `void broadcast(Cell<T>& aData, proc_t aRoot = 0)`

**File:** `commtools.hpp:748`

---

#### Vector Broadcast

```cpp
Vector<real> vec;
if (comm_rank() == 0) {
    vec = {1.0, 2.0, 3.0};
}
broadcast(vec, 0);  // All ranks receive 3-element vector
```

**Signature:** `void broadcast(Vector<T>& aData, proc_t aRoot = 0)`

**File:** `commtools.hpp:1076`

---

#### Matrix Broadcast

```cpp
Matrix<real> mat;
if (comm_rank() == 0) {
    mat.set_size(3, 3);
    // fill mat...
}
broadcast(mat, 0);  // All ranks receive 3×3 matrix
```

**Signature:** `void broadcast(Matrix<T>& aData, proc_t aRoot = 0)`

**File:** `commtools.hpp:1879`

---

## Point-to-Point Communication

> **Warning:** Do not mix raw `MPI_Send`/`MPI_Recv` calls with BELFEM comm functions unless you manage tags explicitly. BELFEM assumes exclusive control of tag generation via `comm_tag()`, and mixing can cause extremely hard-to-debug message mismatches.

### `send(data, target)` - Send Data to Specific Process

#### Send Scalar

```cpp
if (comm_rank() == 0) {
    real value = 3.14;
    send(value, 1);  // Send to rank 1
}
```

**Signature:** `void send(const T aData, proc_t aTarget = 0)`

**File:** `commtools.hpp:314`

---

#### Send Raw Array

```cpp
if (comm_rank() == 0) {
    int data[100];
    // fill data...
    send(data, 100, 1);  // Send 100 elements to rank 1
}
```

**Signature:** `void send(T* aData, index_t aLength, proc_t aTarget)`

**File:** `commtools.hpp:401`

---

#### Send Cell

```cpp
if (comm_rank() == 0) {
    Cell<int> cell = {1, 2, 3, 4};
    send(cell, 1);  // Send to rank 1
}
```

**Signature:** `void send(Cell<T>& aData, proc_t aTarget = 0)`

**File:** `commtools.hpp:587`

---

#### Send Vector

```cpp
if (comm_rank() == 0) {
    Vector<real> vec = {1.0, 2.0, 3.0};
    send(vec, 1);  // Send to rank 1
}
```

**Signature:** `void send(Vector<T>& aData, proc_t aTarget = 0)`

**File:** `commtools.hpp:915`

---

#### Send Matrix

```cpp
if (comm_rank() == 0) {
    Matrix<real> mat(3, 3);
    // fill mat...
    send(mat, 1);  // Send to rank 1
}
```

**Signature:** `void send(Matrix<T>& aData, proc_t aTarget = 0)`

**File:** `commtools.hpp:1949`

---

### `receive(data, source)` - Receive Data from Specific Process

#### Receive Scalar

```cpp
if (comm_rank() == 1) {
    real value;
    receive(value, 0);  // Receive from rank 0
}
```

**Signature:** `void receive(T& aData, proc_t aSource = 0)`

**File:** `commtools.hpp:351`

---

#### Receive Raw Array

```cpp
if (comm_rank() == 1) {
    int data[100];
    index_t length = 100;
    receive(data, length, 0);  // Receive from rank 0, length updated
}
```

**Signature:** `void receive(T* aData, index_t& aLength, proc_t aSource)`

**Important:** `aLength` is both input (allocated size) and output (received size)

**File:** `commtools.hpp:491`

---

#### Receive Cell

```cpp
if (comm_rank() == 1) {
    Cell<int> cell;
    receive(cell, 0);  // Automatically resized to received size
}
```

**Signature:** `void receive(Cell<T>& aData, proc_t aSource = 0)`

**File:** `commtools.hpp:672`

---

#### Receive Vector

```cpp
if (comm_rank() == 1) {
    Vector<real> vec;
    receive(vec, 0);  // Automatically resized
}
```

**Signature:** `void receive(Vector<T>& aData, proc_t aSource = 0)`

**File:** `commtools.hpp:1000`

---

#### Receive Matrix

```cpp
if (comm_rank() == 1) {
    Matrix<real> mat;
    receive(mat, 0);  // Automatically resized
}
```

**Signature:** `void receive(Matrix<T>& aData, proc_t aSource = 0)`

**File:** `commtools.hpp:2047`

---

## Collective Operations

> **Symmetry Rule:** `distribute()` is the inverse of `collect()` only if:
> - Communicator size is unchanged between calls
> - Container ordering is consistent across all ranks
> - All ranks participate in both operations

### `distribute(data)` — send one element to every process

`distribute` is a **send-only** operation. It sends `data(p)` to rank `p` for every
`p != comm_rank()` and waits on its own sends; it posts no receives, so it never writes into the
container it is given and delivers nothing unless someone is receiving.

> **It needs a matching receive somewhere.** `distribute` and `collect` are not a fixed pair —
> they are a send-side and a receive-side helper over the same symmetric tag space
> (`comm_tag(s,t) == comm_tag(t,s)`), and either composes with the ordinary point-to-point calls:
>
> 1. **Scatter from one rank** — the root calls `distribute`, every other rank calls the scalar
>    `receive(value, root)`.
> 2. **Gather to one rank** — the workers call `send(data, root)` and the root calls `collect`.
>
> *Pattern 1* below uses both of these, one after the other.
> 3. **All-to-all** — every rank calls `distribute`, then every rank calls `collect`. See the
>    size caveat below before using this one.
>
> An unmatched transfer is not harmless. The send may block, or it may sit queued and be picked
> up by a **later** exchange between the same rank pair, which share a tag. The **`distribute`
> side** must supply a container of size `comm_size()` (`commtools.hpp:813`); a receiver using
> the scalar `receive` supplies only its own variable.

> **Size caveat on the all-to-all form.** `distribute` waits in `MPI_Waitall` before it returns,
> so every rank is still inside `distribute` when the matching receives ought to be posted. That
> only works while the sends complete into eager buffers. For `distribute( Cell<T> )` — one
> element of `T` per peer — every MPI implementation in practice buffers that eagerly, though
> **the standard guarantees nothing**, which is why the shipped test
> (`tests/comm/test_CommMPI.cpp`, `DistributeCollectScalar` around line 581) passes. The `Cell<Vector<T>>`, `Cell<Cell<T>>` and
> `Cell<Matrix<T>>` overloads move arbitrary payloads and can exceed the rendezvous threshold,
> where `distribute(); collect()` on every rank deadlocks. Use a root-centred pairing for those,
> or post receives before sends with raw MPI.

#### Distribute and collect a Cell

```cpp
// every rank sends one value to every other rank ...
Cell<int> send_data(comm_size(), 0);
for (proc_t p = 0; p < comm_size(); ++p) {
    send_data(p) = comm_rank() * 100 + p;   // what I send to rank p
}
distribute(send_data);

// ... and every rank gathers the value each other rank addressed to it
Cell<int> recv_data;
collect(recv_data, comm_rank() * 100 + comm_rank());   // my own slot
// recv_data(p) now holds what rank p sent me: p * 100 + comm_rank()
```

Note that the received values land in the `collect` container, **not** back in the one passed to
`distribute` — `distribute` never modifies its argument.

**Signature:** `void distribute(Cell<T>& aData)`

**Requires:** `aData.size() == comm_size()`

**File:** `commtools.hpp:804`

---

The overloads below take the same shape as the `Cell<T>` case above: each is **send-only** and
needs a matching receiver — `collect` on every rank, or a per-rank `receive` of the corresponding
type (`distribute( Cell<Vector<T>> )` scatters against `receive( Vector<T>, root )`, not only the
scalar form). None of them writes into the container it is given.

#### Distribute Vector

```cpp
Vector<int> data(comm_size());
// Similar to Cell
distribute(data);
```

**Signature:** `void distribute(Vector<T>& aData)`

**File:** `commtools.hpp:1131`

---

#### Distribute Cell of Vectors

```cpp
Cell<Vector<real>> data(comm_size(), Vector<real>());   // size it: Cell(n) only reserves
if (comm_rank() == 0) {
    for (proc_t p = 0; p < comm_size(); ++p) {
        data(p).set_size(p + 1);  // Variable-length vectors
        // fill data(p)...
    }
}
distribute(data);
// the values arrive in the matching collect(), or in the receivers' receive()
```

**Signature:** `void distribute(Cell<Vector<T>>& aData)`

**File:** `commtools.hpp:1239`

---

#### Distribute Cell of Cells

```cpp
Cell<Cell<int>> data(comm_size(), Cell<int>());   // size it: Cell(n) only reserves
// Similar to Cell<Vector<T>>
distribute(data);
```

**Signature:** `void distribute(Cell<Cell<T>>& aData)`

**File:** `commtools.hpp:1325`

---

#### Distribute Cell of Matrices

```cpp
Cell<Matrix<real>> data(comm_size(), Matrix<real>());   // size it: Cell(n) only reserves
if (comm_rank() == 0) {
    for (proc_t p = 0; p < comm_size(); ++p) {
        data(p).set_size(3, 3);  // Each rank gets 3×3 matrix
        // fill data(p)...
    }
}
distribute(data);
```

**Signature:** `void distribute(Cell<Matrix<T>>& aData)`

**File:** `commtools.hpp:2144`

---

#### Distribute Raw Array with Offsets

```cpp
Vector<int> offsets(comm_size() + 1);
// offsets[i] = start index for rank i
const real* data = ...;  // Pre-allocated array
distribute(data, offsets);
```

**Signature:** `void distribute(const T* aData, const Vector<U>& aOffsets)`

**Requires:** `aOffsets.length() == comm_size() + 1`

**File:** `commtools.hpp:1414`

---

### `collect(data, myValue)` — receive one element from every process

`collect` is the **receive-side** helper. It resizes `data` to `comm_size()`, stores `myValue` in
its own slot, and posts a receive against every other rank. It sends nothing, so it blocks until
each of those ranks has sent — either from a matching `distribute`, or from an ordinary
`send(data, root)`, which is how *Pattern 1* gathers worker results. `collect` on its own, with
nobody sending, is a hang rather than a gather.

**Signature:** `void collect(Cell<T>& aData, const T aMyValue = 0)`

**File:** `commtools.hpp:859`

---

The overloads below share the semantics above: each is **receive-only**, and needs every other
rank to be sending — through `distribute`, or through an ordinary `send( data, root )` of the
matching type. The one exception is `collect( T*, offsets )`, which is a narrower shape: it
gathers sizes from everyone but receives payloads only from ranks `1 … N-1`, taking its own
offset from `aOffsets(0)` (`commtools.hpp:1494-1565`).

#### Collect Scalars into Vector

```cpp
Vector<real> all_values;
real my_value = comm_rank() * 1.5;
collect(all_values, my_value);   // requires every other rank in a matching distribute()
```

**Signature:** `void collect(Vector<T>& aData, const T aMyValue = 0)`

**File:** `commtools.hpp:1186`

---

#### Collect Vectors into Cell

```cpp
Vector<real> my_vec(comm_rank() + 1);  // Variable-length per rank
// fill my_vec...

Cell<Vector<real>> all_vecs;
collect(all_vecs, my_vec);
// all_vecs(i) holds the vector rank i sent -- every other rank must be in a
// matching distribute(), or sending to this rank with send()
```

**Signature:** `void collect(Cell<Vector<T>>& aData, Vector<T> aMyData = {})`

**File:** `commtools.hpp:1577`

---

#### Collect Matrices into Cell

```cpp
Cell<Matrix<real>> all_matrices;
collect(all_matrices);
// receive-only: the other ranks must be sending, via distribute() or send()
```

**Signature:** `void collect(Cell<Matrix<T>>& aData)`

**File:** `commtools.hpp:2287`

---

#### Collect Raw Array with Offsets

```cpp
real* data = ...;  // Pre-allocated
Vector<int> offsets(comm_size() + 1);
collect(data, offsets);
```

**Signature:** `void collect(T* aData, const Vector<U>& aOffsets)`

**File:** `commtools.hpp:1494`

---

### `share(data)` — chunked one-to-many send, with `receive` on the other side

`share` sends the caller's vector to every other rank, in chunks. It is **not collective**: it
posts sends only, so the ranks that are meant to get the data must call `receive`. Keep the
`if/else` rank guard below — a rank that calls `share` when it should be receiving leaves every
other rank waiting on data that never comes. (`receive` on the sender is harmless: it returns
immediately when `comm_rank() == aSource`, `commtools.hpp:1007-1008`.)

```cpp
Vector<real> data;
if (comm_rank() == 0)
{
    data = {1, 2, 3};
    share(data);        // root only — chunks the message
}
else
{
    receive(data);      // non-root only — receives the chunks and resizes
}
```

Use this rather than `broadcast` for large or variable-size `Vector<T>` and `Cell<T>` payloads:
container `broadcast` moves the whole payload in a single unchunked `MPI_Ibcast` after a second
one carrying the size (`commtools.hpp:1076-1120`), which this project has observed to fail on
large datasets (see `doc/coding_philosophy.md`). There is no `share` overload for scalars — use `broadcast(T&)` for those.

**Signature:** `void share(Vector<T>& aData)`

---

## String Communication

### String Operations

```cpp
// Broadcast string cell
Cell<string> messages;
if (comm_rank() == 0) {
    messages = {"hello", "world"};
}
broadcast(messages, 0);

// Send string
if (comm_rank() == 0) {
    send("Hello from rank 0", 1);
}

// Receive string
if (comm_rank() == 1) {
    string msg;
    receive(msg, 0);
}
```

**Signatures:**
- `void broadcast(Cell<string>& aData, proc_t aRoot = 0)`
- `void send(const string& aMessage, proc_t aTarget = 0)`
- `void receive(string& aMessage, proc_t aSource = 0)`

**File:** `commtools.hpp:2443-2450` (implementation in `commtools.cpp`)

---

## Common Patterns

### Pattern 1: Master-Worker Communication

```cpp
if (comm_rank() == 0) {
    // Master: distribute work
    Cell<int> work_items(comm_size(), 0);   // size it: Cell(n) only reserves
    for (proc_t p = 0; p < comm_size(); ++p) {
        work_items(p) = p * 100;  // Work ID for each rank
    }
    distribute(work_items);

    // Collect results
    Cell<Vector<real>> results;
    collect(results);
} else {
    // Worker: receive work
    int work_id;
    receive(work_id, 0);

    // Do work
    Vector<real> result = compute(work_id);

    // Send result
    send(result, 0);
}
```

---

### Pattern 2: Parallel Domain Assembly

```cpp
// Each rank assembles local domain
Matrix<real> local_stiffness = assemble_local_domain();

// Gather all local matrices to rank 0. collect() only receives, so the
// workers must be sending -- see Pattern 1.
Cell<Matrix<real>> all_matrices;
if (comm_rank() == 0) {
    collect(all_matrices);
    all_matrices(0) = local_stiffness;   // collect() skips its own slot and
                                         // leaves it an empty matrix
} else {
    send(local_stiffness, 0);
}

if (comm_rank() == 0) {
    // Combine matrices into global system
    Matrix<real> global_K = combine(all_matrices);
}
```

---

### Pattern 3: Load Balancing

```cpp
// Gather work counts from all ranks. collect() only receives, so every rank
// must also be in a distribute() -- here each rank sends its count to everyone.
Vector<int> my_counts(comm_size(), compute_local_work_count());
distribute(my_counts);

Vector<int> work_counts;
int my_work = compute_local_work_count();
collect(work_counts, my_work);

if (comm_rank() == 0) {
    // Compute new distribution
    Vector<int> new_distribution = balance_load(work_counts);
    distribute(new_distribution);
} else {
    int my_new_work;
    receive(my_new_work, 0);
}
```

---

### Pattern 4: Ghost Node Exchange

```cpp
// Define neighbors
Cell<proc_t> neighbors = get_my_neighbors();

// Send boundary data to neighbors
for (proc_t neighbor : neighbors) {
    Vector<real> boundary_data = extract_boundary(neighbor);
    send(boundary_data, neighbor);
}

// Receive ghost data from neighbors
for (proc_t neighbor : neighbors) {
    Vector<real> ghost_data;
    receive(ghost_data, neighbor);
    update_ghosts(neighbor, ghost_data);
}
```

> **Caveat, unresolved.** `send( Vector<T> )` completes its own `MPI_Waitall` before returning
> (`commtools.hpp:983`). Two ranks that are each other's neighbors are therefore both inside
> the send loop before either reaches its receive loop. For payloads small enough to be sent
> eagerly this is fine, and it is what the code does today; for payloads over the MPI
> implementation's rendezvous threshold it can deadlock. Whether any BELFEM ghost exchange
> actually crosses that threshold has **not** been measured — treat the pattern as sound only for
> small ghost layers. For large symmetric exchanges, interleaving is not enough (both peers can
> still reach their send first): post the receives before the sends with non-blocking MPI calls,
> or use `MPI_Sendrecv`.

---

### Pattern 5: Synchronous Timing

```cpp
comm_barrier();
Timer timer;

// Expensive operation
solve_system();

comm_barrier();
real elapsed = timer.stop();

if (comm_rank() == 0) {
    std::cout << "Solve time: " << elapsed << " seconds\n";
}
```

---

## Non-MPI Builds

When compiled without `BELFEM_MPI`, all communication functions compile to no-ops:

- `comm_size()` returns 1
- `comm_rank()` returns 0
- `broadcast()`, `send()`, `receive()` do nothing
- `distribute()` and `collect()` do nothing either — `collect()` does **not** resize the container or store `myValue`, so code that reads `data(comm_rank())` after `collect` must handle the serial case itself
- `comm_barrier()` does nothing

> **Guarantee:** Non-MPI builds compile to no-ops, not to serial equivalents — `collect()` leaves its container untouched. Use non-MPI builds for serial debugging, not for performance benchmarking.

This allows writing parallel code that degrades to serial execution, provided the caller sizes what `collect()` would have filled.

```cpp
// This code works in both MPI and non-MPI builds
Vector<real> data;
if (comm_rank() == 0) {
    data = load_data();
}
broadcast(data, 0);  // No-op in serial build
process(data);       // All ranks have data
```

---

## Compile-Time Configuration

The comm module behavior is controlled by the `USE_*` CMake **options** (the
`BELFEM_*` names are the *compile definitions* those options generate — they
are not cache entries and setting them on the command line does nothing):

| CMake Option | Generated define | Effect | Default |
|--------------|------------------|--------|---------|
| `USE_MPI` | `BELFEM_MPI` | Master switch for MPI functionality | ON |
| `USE_STRUMPACK` | `BELFEM_STRUMPACK` | `MPI_Init_thread()` with `MPI_THREAD_MULTIPLE` (`cl_Communicator.cpp`, search `MPI_Init_thread`) | ON |
| `USE_PETSC` | `BELFEM_PETSC` | PETSc integration (auto-initialize) | ON |

### Build Examples

```bash
# Default build already has MPI + STRUMPACK + PETSc

# Non-MPI build (serial fallback)
cmake -DUSE_MPI=OFF ..

# MPI without STRUMPACK
cmake -DUSE_STRUMPACK=OFF ..
```

---

## Performance Tips

### 1. **Minimize Communication**

```cpp
// BAD: Send each element separately
for (int i = 0; i < n; ++i) {
    send(data[i], target);  // n MPI calls!
}

// GOOD: Send all at once
send(data, n, target);  // 1 MPI call
```

---

### 2. **Let `collect()` post the receives**

```cpp
// SLOWER: the root waits for each rank in turn
if (comm_rank() == 0) {
    for (proc_t p = 1; p < comm_size(); ++p) {
        receive(data[p], p);
    }
} else {
    send(my_data, 0);
}

// BETTER: same shape, but the root posts every receive at once
if (comm_rank() == 0) {
    collect(data, my_data);
} else {
    send(my_data, 0);
}
```

Both are correct; the second overlaps its receives instead of serializing them. Note the rank
guard is still needed — `collect` receives only, so calling it on every rank leaves nobody
sending.

---

### 3. **Overlap Communication and Computation**

```cpp
// Start non-blocking communication
MPI_Request req;
MPI_Isend(..., &req);

// Do independent work while waiting
compute_local_data();

// Wait for completion
MPI_Wait(&req, MPI_STATUS_IGNORE);
```

**Note:** BELFEM's comm layer uses non-blocking internally but waits immediately. For overlap, use raw MPI.

---

### 4. **Balance Load Before Communication**

```cpp
// Ensure all ranks have similar work
balance_workload();

// Then communicate (reduces idle time)
comm_barrier();
distribute(data);
```

---

### 5. **Avoid Unnecessary Barriers**

```cpp
// BAD: Excessive synchronization
for (int iter = 0; iter < niter; ++iter) {
    compute();
    comm_barrier();  // Often unnecessary!
}

// GOOD: Only barrier when needed
for (int iter = 0; iter < niter; ++iter) {
    compute();
}
comm_barrier();  // Once at end
```

---

## Thread Safety

### MPI Thread Support Levels

BELFEM's thread support depends on compile-time configuration:

| Build Configuration | MPI Thread Level | Behavior |
|---------------------|------------------|----------|
| Default (no threading libs) | `MPI_Init()` | Single-threaded MPI (MPI_THREAD_SINGLE) |
| With `BELFEM_STRUMPACK` | `MPI_Init_thread(..., MPI_THREAD_MULTIPLE, ...)` | Requests highest thread support |
| With `BELFEM_PETSC` but without `BELFEM_STRUMPACK` | `MPI_Init()` | Single-threaded; STRUMPACK, when present, decides the level regardless of PETSc |

> **Important:** When BELFEM is compiled with STRUMPACK support, it initializes MPI with `MPI_THREAD_MULTIPLE` (the highest thread support level). If the MPI implementation cannot provide this level, a warning is printed to stderr, but initialization continues with the highest available level.

### Thread Safety Rules

> **Critical:** Do NOT call any BELFEM comm function from multiple threads concurrently, even when `MPI_THREAD_MULTIPLE` is initialized.

**Why?** Every rank pair shares exactly two tags on `gComm.world()` (`comm_tag()` is symmetric in source and target) and MPI matches per (source, tag) in FIFO order, so two threads issuing BELFEM comm calls to the same peer interleave their size/payload messages and each consumes the other's. The per-call request buffers themselves are `malloc`ed per call and private.

**Safe threading approaches:**

1. **MPI_THREAD_FUNNELED** — Only master thread calls MPI:

```cpp
#ifdef OMP   // the define the compiler config sets under USE_OPENMP
#pragma omp parallel
{
    // Compute in parallel
    Vector<real> local_result = compute();

    #pragma omp master
    {
        // Only master thread communicates
        send(local_result, 0);
    }
}
#endif
```

2. **MPI_THREAD_SERIALIZED** — Serialize MPI calls with mutex:

```cpp
std::mutex mpi_mutex;

#pragma omp parallel
{
    // Compute in parallel
    Vector<real> result = compute();

    // Serialize MPI calls
    {
        std::lock_guard<std::mutex> lock(mpi_mutex);
        send(result, 0);
    }
}
```

3. **Thread-local buffers** — Each thread uses separate data:

```cpp
#pragma omp parallel
{
    int tid = omp_get_thread_num();

    // Each thread computes independently
    Vector<real> my_result = compute(tid);

    // Master collects from all threads sequentially
    #pragma omp master
    {
        Cell<Vector<real>> thread_results(omp_get_num_threads(), Vector<real>());   // size it: Cell(n) only reserves
        thread_results(0) = my_result;

        #pragma omp barrier

        for (int t = 1; t < omp_get_num_threads(); ++t) {
            // Collect from other threads (thread-safe via barrier)
        }

        // Now send combined result via MPI
        send_combined(thread_results);
    }
}
```

> **Warning:** Even with `MPI_THREAD_MULTIPLE`, BELFEM's comm functions are **not internally synchronized**. External synchronization is always required for concurrent access.

---

## Debugging MPI Programs

### 1. **Enable MPI Error Checking**

All BELFEM comm functions route the MPI return code through `comm_check(error_code)`, a `BELFEM_ERROR` that reports the implementation's error string — it throws in a debug build and aborts the job in release (see `assert.hpp`, `throw_on_error()`).

---

### 2. **Print Rank in Output**

```cpp
std::cout << "[Rank " << comm_rank() << "] Message\n";
```

---

### 3. **Use Barriers to Isolate Bugs**

```cpp
comm_barrier();
std::cout << "[Rank " << comm_rank() << "] Reached checkpoint 1\n";
comm_barrier();
```

---

### 4. **Check for Deadlocks**

Use `MPI_Wtime()` to detect hangs:

```cpp
double t0 = MPI_Wtime();
receive(data, source);
double elapsed = MPI_Wtime() - t0;
if (elapsed > 10.0) {
    std::cout << "[Rank " << comm_rank() << "] Possible deadlock!\n";
}
```

---

### 5. **Run with Fewer Processes**

Start debugging with 2 processes, then scale up:

```bash
mpirun -n 2 ./my_app  # Easier to debug than -n 128
```

---

### 6. **Enable MPI Verbose Mode**

Enable runtime diagnostics to catch tag mismatches and communication errors:

```bash
# Open MPI
mpirun --mca mpi_verbose 1 -n 4 ./my_app

# MPICH
mpirun -env MPIR_CVAR_DEBUG_SUMMARY 1 -n 4 ./my_app
```

> **Tip:** MPI verbose flags vary by implementation. Check your MPI documentation for specifics.

---

## Related Modules

- **containers**: `Cell<T>` used extensively in comm operations
- **linalg**: `Vector<T>` and `Matrix<T>` have comm support
- **mesh**: Uses comm for distributed mesh operations
- **sparse**: Distributed sparse matrix solvers via PETSc/MUMPS

---

## See Also

- [MPI Standard Documentation](https://www.mpi-forum.org/docs/)
- [Open MPI Documentation](https://www.open-mpi.org/doc/)
- [Containers Module](../../containers/doc/README.md) - BELFEM container classes
- [Linear Algebra Module](../../linalg/doc/README.md) - Vector and Matrix types
- `CLAUDE.md` (repository root) - Documentation guidelines
