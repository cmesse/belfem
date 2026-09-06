# BELFEM Coding Philosophy {#doc_coding_philosophy}

**Date:** 2026-08-08
**Purpose:** Standards for nomenclature, memory management, container selection, and design patterns to ensure high-performance execution and developer clarity

**Revision History:**
- 2026-01-16: Initial documentation synthesizing comm, core, containers, and linalg design principles
- 2026-06-14: B1-scoped source-reconciliation pass (containers, comm, memory model, naming); per-claim ledger in `tmp/whitepaper/coding_philosophy_CHANGES.md`
- 2026-08-08: Full-document revision from the jury-feedback harvest. Every factual claim re-verified against source (V1–V11 ledger in `tmp/whitepaper/coding_philosophy_CHANGES.md`). Memory chapter rewritten around "allocate deliberately, not manually"; build flags, alignment story, smart-pointer census, `-fno-exceptions` status, and testing posture corrected to match the repository; unverifiable performance figures removed or cited to literature; audit annotations moved to the Source Reconciliation appendix.

---

## Table of Contents

1. [Introduction](#introduction)
2. [Core Principles](#core-principles)
3. [Naming Conventions](#naming-conventions)
4. [Container Selection Strategy](#container-selection-strategy)
5. [Memory Management Philosophy](#memory-management-philosophy)
6. [Performance-First Design Patterns](#performance-first-design-patterns)
7. [Error Handling Strategy](#error-handling-strategy)
8. [Thread Safety and MPI](#thread-safety-and-mpi)
9. [Benchmarks We Owe Ourselves](#benchmarks-we-owe-ourselves)
10. [Summary Checklist](#summary-checklist)
11. [Further Reading](#further-reading)
12. [Appendix: Source Reconciliation](#appendix-source-reconciliation)

---

## Introduction

BELFEM (Berkeley Lab Finite Element Framework) is designed for **high-performance finite element simulations** on HPC systems. Every design decision prioritizes:

1. **Zero abstraction penalty** - code should compile to the same machine instructions as hand-written C
2. **Deterministic performance** - no hidden allocations, no unpredictable pauses, no reference counting on hot paths
3. **Cache-friendly data structures** - contiguous memory, predictable access patterns
4. **Developer clarity** - consistent naming, explicit ownership, readable mathematical code

This document explains the "why" behind BELFEM's coding conventions. For contributors joining from modern C++ projects emphasizing safety and smart pointers, understanding these principles is essential.

**Target audience:** Developers contributing to BELFEM, code reviewers, and anyone wondering why the code looks different from typical C++14/17 practices.

**Provenance of claims:** Facts in this document are anchored to the source with `file:line` citations, current as of 2026-08-08. Code blocks are either verbatim excerpts (cited) or explicitly labeled *pseudocode*. Where a claim could not be verified from the repository, it is either cited to external literature and labeled a **literature heuristic**, or it has been removed — see the [Source Reconciliation appendix](#appendix-source-reconciliation).

---

## Core Principles

### 1. Performance First, Safety Through Testing

**Philosophy:** In HPC, a 5% performance regression across a large simulation costs real money and time. BELFEM optimizes for the **release build**, then adds safety checks in **debug builds**.

**In practice:**
- `BELFEM_ASSERT` for bounds checks and logic validation — **compiled out in release** (`NDEBUG`)
- `BELFEM_ERROR` for runtime failures — **always active**
- Debug/release is selected by the CMake option `USE_DEBUG` (`option( USE_DEBUG ...)` in `CMakeLists.txt`), not by `CMAKE_BUILD_TYPE`. With no `$SCLS` it defaults **OFF**, so an ordinary `cmake ..` is already a release build and a debug build requires an explicit `-DUSE_DEBUG=ON`. The default is not hard-coded: it comes from `${BELFEM_DEFAULT_USE_DEBUG}`, which `config/system/find_scls_flavor.cmake:23` initializes to `OFF` and raises to `ON` only when `$SCLS` names the `debug` flavor (`:47-50`). An explicit `-DUSE_DEBUG=…` always wins.
- Actual release flags, per toolchain: GCC/Clang C and C++ `-O2` with `NDEBUG` defined, GCC Fortran `-O2` plus a native flag that is `-march=native` on Linux x86, `-mtune=native` on Apple x86 and empty elsewhere (`config/compiler/config_gcc.cmake`); ICC `-O1 -xHost`, ICC Fortran `-O1` (`config/compiler/config_icc.cmake`). Debug builds define `DEBUG` and use `-Og -g` for GCC/Clang C++, `-O0 -fcheck=bounds -fbacktrace` for GCC Fortran, and `-O0 -g` for ICC.

**Exceptions status (honest statement):** BELFEM's design intent is to be *compatible with* exception-free production builds — no `throw` on hot paths, error macros that abort in release. However, **no build configuration currently sets `-fno-exceptions`**: the flag appears in no CMake file, and `src/core/assert.hpp:191` contains a `throw` that is compiled into every default build, so an actual `-fno-exceptions` build would not compile today. Treat "exception-free production" as a design goal, not a property of the current build system.

**Testing posture (honest statement):** Tests are GoogleTest suites run via `make check`, gated on `USE_TEST` which defaults **ON** (`option( USE_TEST "Build Tests" ON )` in `CMakeLists.txt`; release policy since 2026-08-14). Since 2026-08-30 a GitLab CI server runs the **full** suite nightly, which retires the older statement that nothing ran the tests automatically.

Three qualifications keep this an honest statement rather than a claim of safety. First, the pipeline definition (`.gitlab-ci.yml` at the repository root) runs only on a scheduled or manually started pipeline, on a single runner tagged `belfem-local`, so the tree records what the nightly builds and tests but not the hardware it runs on, and a clone cannot reproduce it. Second, there is still no sanitizer configuration: `-fsanitize` appears in no CMake file, and Valgrind is invoked only by hand. Third, and most important, a passing nightly measures the suite, not the code. This project's own incident catalog records 539 incidents in its locked catalog (addenda run to INC-564) of which 280 were code defects and 38 were closed by adding a test; the suite is not yet the thing that finds our bugs, and a green pipeline must not be read as though it were.

The check-removal bargain — "we may compile the bounds checks out because testing catches the bugs" — is therefore *partly* earned. The automation half now exists. The coverage half does not: the modules carrying the heaviest recorded defect load are also the thinnest covered (`src/mesh` is 52,318 source lines against 915 test lines). Full payment still waits on ASan/UBSan jobs running on every change, and on coverage reaching the subsystems that actually break; wiring both up is a tracked proposal (see appendix).

**Contrast with modern C++:**
Modern C++ emphasizes "safety by default" (smart pointers, exceptions, bounds-checked containers). BELFEM inverts this: **speed by default, safety in development**.

---

### 2. Minimal Dependencies, Maximum Control

**Philosophy:** Every abstraction layer adds overhead. BELFEM stays close to "bare metal" when it matters.

**In practice:**
- Core utilities depend only on the C++ standard library (no Boost, no external frameworks)
- MPI and the linear algebra backends (Armadillo/Blaze) are the **only** framework-level dependencies
- Custom containers (`Cell`, `Vector`, `Matrix`) wrap STL or backend types but expose raw pointers for zero-cost interfacing

**Why:** HPC centers use diverse architectures (x86, ARM, Power). Minimal dependencies ensure portability and give full control over memory layout.

---

### 3. Explicit Over Implicit

**Philosophy:** Hidden behavior (implicit conversions, reference counting, magic destructors) makes performance analysis impossible.

**In practice:**
- Ownership decisions are made once, at setup, and are visible in the owner's destructor
- No automatic type conversions between numerical types
- Manual `malloc`/`free` in the few classes that own raw buffers (`StringList`, `DynamicBitset`, `ShiftRegister`, mesh-entity internal arrays)
- Memory ownership is **crystal clear** from naming and documentation

---

### 4. Mathematical Readability

**Philosophy:** FEM code implements published algorithms. The code should resemble the equations in papers.

**In practice:**
- Strict prefix rules (`a`, `t`, `m`, `g`) for framework code
- **Relaxed rules** for mathematical kernels (allow `A`, `B`, `x`, `y`, `i`, `j`)
- Inline comments cite literature: `// Equation (2.7) from Messe et al. 2023`

**Example:** A stiffness matrix assembly loop should look like the formula in Bathe, not like Java enterprise code.

---

## Naming Conventions

### Prefix-Based Nomenclature

BELFEM uses a **strict prefix system** to eliminate ambiguity and improve code readability at a glance.

| Prefix | Meaning | Scope | Example |
|--------|---------|-------|---------|
| **`a`** | **Argument** | Function parameters (all public APIs) | `aTarget`, `aNumRows`, `aMessage` |
| **`t`** | **Temporary** | Local variables (function scope) | `tCount`, `tStatus`, `tBuffer` |
| **`m`** | **Member** | Class data members (private/protected) | `mMatrix`, `mInfoLevel`, `mCommRank` |
| **`g`** | **Global** | Framework-wide globals (extern singletons) | `gLog`, `gComm`, `gTbulk` |

### Rationale

**1. Readability at a glance** *(pseudocode)*
```cpp
void send(const int aTarget, const Vector<real> & aData)
{
    int tCommTag = comm_tag( comm_rank(), aTarget );  // Clearly a local temporary
    int tLength  = aData.length();                    // Local, computed from argument
    MPI_Send(aData.data(), tLength, ...);
}
```
You know instantly: `aTarget` is a parameter, `tCommTag` is local.

**2. Prevents shadowing bugs** *(pseudocode)*
```cpp
class Communicator
{
    int mRank;  // Member

public:
    void set_rank(int aRank) {  // Parameter
        mRank = aRank;  // No confusion - names are distinct
    }
};
```

**3. Consistency across 300k lines**
All BELFEM modules (core, comm, fem, mesh, etc.) use identical conventions. Switching between modules requires zero mental context switch.

---

### Mathematical Exceptions

When implementing **numerical algorithms** directly from literature, readability trumps strict prefixes.

**Allowed in mathematical kernels:**
- Matrix operations: `A`, `B`, `K` (stiffness), `M` (mass), `J` (Jacobian)
- Vectors: `x`, `y`, `f` (force), `u` (displacement), `b` (RHS)
- Scalars: `alpha`, `beta`, `norm`, `det`
- Indices: `i`, `j`, `k`, `m`, `n`

In `src/math` and `src/physics` kernels the `a`/`t` prefixes may be dropped entirely so the code reads like the derivation (see `CLAUDE.md` for the canonical thermophysical symbol table). The `m` and `g` prefixes always still apply.

**Where the kernel/framework boundary lies:** prefixes resume the moment a variable is no longer part of the mathematical expression. The same function can contain both regimes:

*(pseudocode)*
```cpp
void solve_system()
{
    Timer tTimer;  // Prefix applies - not a math quantity

    // Math kernel - relaxed
    Matrix<real> K = assemble_stiffness();
    Vector<real> f = assemble_force();
    Vector<real> u = K.solve(f);  // A \ b in math notation

    uint64_t tElapsed = tTimer.stop();  // Prefix returns
    message(InfoLevel::Default, "Solve time: %lu ms", tElapsed);
}
```

---

### Additional Naming Rules

**Types:** `CamelCase` (e.g., `DynamicBitset`, `Matrix`, `InfoLevel`)

**Functions:** `lowercase_with_underscores` (e.g., `set_size()`, `comm_tag()`, `string_to_bool()`). Exception: capitalized mathematical symbols such as `T` (temperature) are preserved in the physics modules to avoid colliding with `t` (time).

**Constants:** `BELFEM_UPPERCASE` for macros (e.g., `BELFEM_EPSILON`, `BELFEM_INT64`)

**Files:** the stem follows the prefix, matching what the file defines:
- `cl_` + CamelCase class name for class files: `cl_Communicator.hpp` defines `class Communicator`, `cl_DynamicBitset.hpp` defines `class DynamicBitset`
- `fn_` + snake_case function name for standalone functions: `fn_sprint.hpp` defines `sprint()`
- `op_` for operator overloads: `op_MatrixPlus.hpp`
- plain `lowercase_with_underscores` for module-level utility headers: `commtools.hpp`, `typedefs.hpp`, `assert.hpp`

Documentation files are always `lowercase_with_underscores.md`.

---

## Container Selection Strategy

BELFEM provides **custom container wrappers** around STL for consistency, safety, and FEM-specific features.

### Container Decision Tree

```
Need to store a collection of objects?
│
├─ Numerical vector/matrix?
│  └─ Use Vector<T> or Matrix<T> (linalg backend integration)
│
├─ Dynamic array of one element type T (IDs, pointers, mesh entities)?
│  └─ Use Cell<T> (BELFEM's std::vector wrapper)
│
├─ Key-value lookup?
│  ├─ Unordered? → Use Map<K,V> (wraps std::unordered_map)
│  └─ Sorted iteration needed? → Use OrderedMap<K,V> (wraps std::map)
│
├─ Unique set of values?
│  ├─ Unordered? → Use Set<T> (wraps std::unordered_set)
│  └─ Sorted? → Use OrderedSet<T> (wraps std::set)
│
├─ Large boolean flags (>1000s)?
│  ├─ Fixed size at compile time? → Use Bitset<N>
│  └─ Runtime size? → Use DynamicBitset
│
├─ Fixed-depth history buffer (time-stepping)?
│  └─ Use ShiftRegister<T>
│
└─ FIFO queue?
   └─ Use Queue<T> (wraps std::queue)
```

### Primary Container: `Cell<T>`

**The workhorse** — use for almost all dynamic arrays.

`Cell<T>` is a thin wrapper around `std::vector<T>` (`src/containers/cl_Cell.hpp:44`). The argument for it is **interface and discipline, not allocation performance** — allocation is `std::vector`'s in both cases. What the wrapper buys:

1. **Debug bounds checking with release speed:** `operator()` carries a `BELFEM_ASSERT` (`cl_Cell.hpp:148-168`) that compiles out under `NDEBUG`, leaving raw element access.
2. **A uniform BELFEM surface:** `set_size()` (`cl_Cell.hpp:186,194`), `push()` (`:259`), `data()` (`:117-127`), `first()`/`last()`, and `operator()` everywhere — one API to learn across the framework. One caveat is baked into history: the size accessor is `size()` on `Cell`/`DynamicBitset` but `length()` on `Vector`/`Matrix` (`cl_Cell.hpp:178`; `cl_AR_Vector.hpp:257`). Closing that split with a zero-cost alias is a tracked proposal; until then, use the accessor the container actually has.
3. **FEM utilities as free functions** (`cl_Cell.hpp:428-548`): `sort( aCell )`, `unique( aCell )`, `reverse( aCell )`, `append( aA, aB )`, `append_move( aTarget, aSource )`, and `find_index_in_unique_cell( aCell, aMember )` (binary search; the cell must be sorted and unique).
4. **Well-defined `data()` even when empty:** `data()` returns `mCell.data()`, which is well-defined — *possibly `nullptr`* — for an empty container. That composes safely with MPI, which accepts any buffer pointer when the count is 0. Do **not** null-test `data()` as an emptiness proxy: an empty Armadillo `Vector::data()` is deterministically `nullptr`, while an empty Blaze one may not be. Test `size()`/`length()` instead.

**Example** (real API — note `size()` and the free-function `unique`):
```cpp
Cell< index_t > tNodeIDs( 100, gNoIndex );   // sized + filled via set_size semantics

// Debug: bounds-checked; Release: raw access
for ( index_t k = 0; k < tNodeIDs.size(); ++k )
{
    tNodeIDs( k ) = compute_id( k );
}

unique( tNodeIDs );                          // free function, in-place sort+unique

// Direct MPI access — data() is well-defined even if the cell is empty
MPI_Send( tNodeIDs.data(), tNodeIDs.size(), MPI_UINT64_T, ... );
```

---

### Numerical Containers: `Vector<T>` and `Matrix<T>`

**Use ONLY for linear algebra** — not general-purpose arrays.

**Why separate from `Cell`?**
- Backend integration: wrap Armadillo (`arma::Mat`, `src/linalg/armadillo/cl_AR_Vector.hpp:31,37`) or Blaze (`blaze::DynamicVector`/`DynamicMatrix`, `src/linalg/blaze/cl_BZ_Vector.hpp:38,45`, `cl_BZ_Matrix.hpp:26`) for BLAS/LAPACK acceleration
- Mathematical operators: `+`, `-`, `*`, `/` overloaded for matrix arithmetic
- Layout control: column-major storage for BLAS/LAPACK and Fortran interop (see "Matrix Memory Layout" below)
- Alignment: the **backend** allocates via `posix_memalign` — Blaze's `AlignedAllocator` aligns to the SIMD width of the build (32 bytes under AVX2, 64 only under AVX-512); Armadillo aligns to 32 bytes for allocations ≥ 1 KiB, 16 below, and deliberately never uses 64. BELFEM itself performs no aligned allocation anywhere (see "Alignment Reality" below).

**Example** *(pseudocode)*:
```cpp
// Matrix assembly - use Matrix, not Cell<Cell<real>>
Matrix<real> K(num_dofs, num_dofs);  // Stiffness matrix

for (Element * e : elements) {
    assemble(e->stiffness(), e->dofs(), K);
}

Vector<real> u = K.solve(f);  // Solve via backend (calls LAPACK internally)
```

**Don't use `Vector<T>` for ID lists, node collections, etc.** — use `Cell<T>` instead.

#### Matrix Memory Layout

BELFEM matrices are **column-major** (Fortran convention). Column `j` is contiguous in memory; row `i` is strided.

**Why column-major:**
- Direct interop with BLAS, LAPACK, MKL, MUMPS, STRUMPACK, PARDISO, PETSc — all assume Fortran layout
- Armadillo is natively column-major; the Blaze backend is pinned to it via `BLAZE_DEFAULT_STORAGE_ORDER = blaze::columnMajor` (`src/linalg/blaze/blaze_config.hpp:55`)
- Avoids transposes when calling LAPACK from BELFEM

**Backend caveat — `data()` interpretation:**

Under **Armadillo**, the storage is exactly `A(i, j) == data[j * nrows + i]`, and `capacity()` equals `n_rows * n_cols` (`cl_AR_Matrix.hpp:224-226`).

Under **Blaze**, columns are padded for SIMD alignment (padding is on by default and BELFEM does not disable it), so the inter-column stride — `spacing()`, the "leading dimension" in BLAS terms — may exceed `nrows`, and `capacity()` exceeds `n_rows * n_cols`. The LAPACK wrappers already honor this: "passing n_rows() instead would silently corrupt the result" (`src/linalg/lapack/lapacktools.hpp:150-153`).

**Rule:** always use the `A(i, j)` accessor for element access. Never compute offsets from `data()` directly. Reserve `data()` for bulk operations (BLAS/LAPACK calls that take a leading-dimension argument, MPI transfers that treat the buffer opaquely).

For MPI, the comm layer transmits the padded footprint of the current shape — `spacing() * n_cols` elements from `data()` — so the padded layout round-trips bit-identically between ranks running the same binary (`src/comm/commtools.hpp` documents this design in place; Blaze's spacing is a pure function of the row count, so both sides always agree). It deliberately does *not* transmit `capacity()`: a matrix that shrank keeps its old, larger allocation, and that stale count would overflow the exact-fit buffer on the receiving side (defect found and fixed 2026-08-09; the receive paths now also guard `capacity() >= transfer length` with `BELFEM_ERROR`). A Blaze build's serialized matrix is still NOT byte-compatible with an Armadillo build's — relevant only if matrices are checkpointed across builds with different backends.

**Iteration order matters.** The inner loop must vary the row index for cache-friendly access:

```cpp
Matrix<real> A(nrows, ncols);

// FAST — inner loop over rows (contiguous access)
for (uint j = 0; j < ncols; ++j) {
    for (uint i = 0; i < nrows; ++i) {
        A(i, j) = ...;
    }
}

// SLOW — inner loop over columns (strided by the column spacing; cache thrash)
for (uint i = 0; i < nrows; ++i) {
    for (uint j = 0; j < ncols; ++j) {
        A(i, j) = ...;
    }
}
```

For matrices larger than the cache, the slow form misses on every inner-loop access — with 8-byte reals on 64-byte cache lines, up to 8 loads per line instead of 1. *(Literature heuristic, not a BELFEM benchmark: Drepper, "What Every Programmer Should Know About Memory", §3.3; Agner Fog, "Optimizing software in C++", §9.)*

**Watch out:** many C++ libraries default to row-major (Eigen, NumPy, the default mental model from C-array pedagogy). BELFEM does not.

---

### Specialized Containers

**`DynamicBitset`** — Runtime-sized bit array
- **Use for:** Topology flags (active nodes, boundary faces), boolean masks
- **Memory:** 1 bit per value, packed into `uint64_t` words (a `Cell<bool>` wraps `std::vector<bool>`, whose packing is implementation-defined and exposes no word-level access)
- **Operations:** Fast bitwise AND/OR/XOR, `count()`, sparse/dense conversions
- **Implementation:** Manual `uint64_t*` allocation; `lock()` computes an FNV-1a hash and freezes the bitset so it can serve as a hashable key (`cl_DynamicBitset.hpp:516-526`) — the locking is about immutability for hashing, not thread safety

**`ShiftRegister<T>`** — Fixed-depth history buffer
- **Use for:** Time-stepping with revert capability (adaptive time integration)
- **Memory:** Preallocates all slots once; no reallocation during stepping
- **Operations:** `push()`, `revert()` (rollback failed timestep), indexed access to history
- **Implementation:** `std::malloc` buffer with placement-new slot construction and explicit destructor calls — the one sanctioned pool pattern in the tree (see "The pool pattern, done right" below)

**`Bitset<N>`** — Compile-time sized bit array
- **Use for:** Small fixed flags (element type masks, DOF activity)
- **Memory:** Stack-allocated, zero runtime overhead

**`StringList`** — C-compatible string array
- **Use for:** I/O with C libraries (Exodus, HDF5), MPI string communication
- **Memory:** Manual `char**` allocation, each string a separately owned `char*` (`cl_StringList.cpp:25,49`)
- **Operations:** Conversion to/from `Cell<string>`, direct C API passing

---

### Why Not STL Directly?

**Reason 1: One interface, one discipline.** Every module uses the same containers with the same conventions (`set_size`, `push`, `operator()`, `data()`), so switching modules costs nothing. The wrapper is also the single place where debug checks, FEM utilities, and future instrumentation live. This is the "deep module" argument (Ousterhout): a small interface hiding real functionality beats re-deciding vector idioms at every call site.

**Reason 2: Debug safety with release speed.** `operator()` asserts in debug and compiles to raw access in release — `std::vector::operator[]` gives you only the latter, `at()` only the former (always-on, throwing).

**Reason 3: Backend independence.** `Vector`/`Matrix` code is written once against the BELFEM API and runs on either Armadillo or Blaze; the wrapper absorbs the differences (padding, `memptr()` vs `data()`, resize semantics).

**Reason 4: Errors defined out of existence.** `data()` on an empty container is well-defined and composes with zero-count MPI calls; sized constructors initialize with sentinels; utilities like `unique()` and `find_index_in_unique_cell()` encode the FEM idiom once, correctly.

---

### Entity Flags over Maps (2026-08-09)

Every graph vertex — and therefore every mesh entity (Node, Edge, Face,
Element, Facet) — carries an 8-slot flag bitset
(`graph::Vertex::mFlags`, `uint8_t`; `flag( aIndex = 0 )` /
`unflag( aIndex = 0 )` / `is_flagged( aIndex = 0 )`).

**Prefer flags over `Map` lookups** for marking, dedup, and visited-set
logic on mesh entities: zero allocation, O(1) by construction, and the
entity itself carries the state.

**Slot convention:**

- **Slot 0 (the default argument)** is reserved for flag protocols that
  span *between* functions — one function flags, another consumes
  (e.g. `flag_periodic_nodes()` → periodic sideset creation).
- **Slots 1–7** are for *function-local* algorithms; pick a numbered
  slot to avoid collisions with any slot-0 protocol in flight. Multiple
  slots may be combined when an algorithm needs several independent
  marks per entity.

**Hygiene:** flags are global mesh state. A function-local flag use must
clear the touched range *before* the gather (stale flags from an earlier
user must not leak in) and *after* it (the next user must find a clean
slate). The pre-clear is required even on numbered slots — no convention
guarantees a slot was left clean.

```cpp
// function-local dedup on slot 0 requires both clears
for ( index_t e : tEdgeIndices ) { tEdges( e )->unflag_nodes(); }   // pre
for ( index_t e : tEdgeIndices )
{
    Edge * tEdge = tEdges( e );
    for ( uint k=0; k<tEdge->number_of_nodes(); ++k )
    {
        Node * tNode = tEdge->node( k );
        if ( ! tNode->is_flagged() )
        {
            tNode->flag();
            tUniqueNodes.push( tNode );
        }
    }
}
for ( index_t e : tEdgeIndices ) { tEdges( e )->unflag_nodes(); }   // post
```

---

## Memory Management Philosophy

### Allocate Deliberately, Not Manually

The lead principle is **not** "manual is faster than smart pointers." It is: **allocation, ownership, and layout decisions are made deliberately during setup, so that the solve phase runs on stable, compact, predictable data structures.** Manual allocation is one tool for that; so are the backend allocators, the STL wrapped inside `Cell`, and — in setup-scope code — smart pointers.

Two rules make the deliberateness concrete:

**The entity rule.** Objects instantiated at mesh scale (nodes, elements, facets, edges, faces) contain only fixed-size state, compact handles/indices, and non-owning pointers. Variable-size storage belongs to aggregate owners or flattened tables — never to the entity itself growing ad hoc. (An element's node-pointer array is fixed at its type's node count, allocated once in the constructor: `cl_ElementTemplate.hpp:383,404,426`.)

**The ownership rule.** Ownership is hierarchical and established at setup: the mesh owns its entities, blocks own their elements, the kernel owns its managers. Shared ownership is banned in computational data structures. A raw `T*` is the explicit non-owning borrow.

**What manual/deliberate allocation buys, stated precisely:**
- No per-entity heap traffic *during the solve* — everything is allocated before time-stepping begins
- Direct C interop: `data()` pointers go straight into MPI, BLAS, LAPACK, and Fortran solvers with no marshaling
- Visible cleanup: every owner's destructor lists what it frees; there is no unwind machinery to reason about
- Compatibility with an exception-free discipline: no error path on the hot loop depends on `throw` (see the honest exceptions statement in Core Principles — the `-fno-exceptions` *build* is a goal, not yet a fact)

**What it does NOT buy (claims we no longer make):**
- ~~"`unique_ptr` needs two allocations"~~ — false; a default-deleter `unique_ptr` is a single pointer with a zero-size deleter and exactly one allocation. The control-block/atomic-refcount cost is real but belongs to **`shared_ptr` only**.
- ~~"Manual `free` is deterministic, RAII is not"~~ — false; RAII destruction is fully deterministic. The honest distinction is *visibility*: an explicit `free`/`delete` in the destructor body is easier to audit than destruction implied by scope exit, and it carries no unwind-table baggage. That is a readability argument, and we make it as such.
- A `std::unique_ptr<T[], Deleter>` around an aligned buffer is functionally equivalent to a hand-written RAII wrapper. BELFEM prefers its own named wrappers (`StringList`, `DynamicBitset`, `ShiftRegister`) for uniformity of interface and naming — **a style choice, and we own it as one**.

---

### Ownership Annotations

Ownership is expressed by convention, in the spirit of the C++ Core Guidelines' `owner<T*>` / raw-`T*` split:

- **A plain `T*` is a non-owning borrow, always.** A `Solver` takes `Mesh*` and never deletes it.
- **An owning pointer lives inside a container held by its owner**, and the owner's destructor is the single place that deletes. Where a header documents members, owning pointer members are annotated (`// Owns: deleted in destructor` vs `// Borrows: does not delete`).
- `Cell<T*>` **as a container is always non-owning** — `~Cell()` is `= default` and never touches pointees (`cl_Cell.hpp:113`). Whether the *pointers inside it* are owning is the holder's contract, declared by its destructor.

The sanctioned owning-holder pattern, verbatim from the tree — `SideSet` holds both flavors side by side (`cl_SideSet.hpp:43-44`, `cl_SideSet.cpp:34-41`):

```cpp
Cell< Facet * > mFacets;   // owning — deleted below
Cell< Node *  > mNodes ;   // borrowed from the mesh — never deleted

SideSet::~SideSet()
{
    for ( Facet * tFacet : mFacets )
    {
        delete tFacet ;
    }
    mFacetMap.clear() ;
}
```

Zero-cost `owner<T*>`/`observer<T*>` aliases plus a clang-tidy ownership check would turn this convention into tooling; that is a tracked proposal (appendix), not yet in the tree.

---

### The Canonical Lifecycle: Per-Object Allocation, Hierarchical Ownership

This is how mesh entities actually live and die — per-object `new` at setup, hierarchical `delete` at teardown. (There is no pool or arena allocator for mesh entities in BELFEM; an earlier revision of this document showed one as if it existed.)

**1. Allocate & construct (setup).** The element factory returns one heap object per element (`cl_Element_Factory.cpp:65-193` — a pure creator with no state, `cl_Element_Factory.hpp:32-41`); the Gmsh reader news each node (`cl_Mesh_GmshReader.cpp:277`):

```cpp
return new ElementTemplate< 1, 1, 0, 0, 0 >( aID );          // factory, one per element
mNodes( tCount++ ) = new Node( tID, tX, tY, tZ );            // reader, one per node
```

**2. Use (solve).** Entities are reached through non-owning `Cell<T*>` views and plain `T*` borrows. Nothing allocates or frees on this path.

**3. Destroy (teardown).** Ownership is layered: `Mesh` deletes nodes, edges, faces, vertices, and the containers of aggregates (`cl_Mesh.cpp:205-290`); `Block` deletes its elements (`cl_Block.cpp:34-40`); `SideSet` deletes its facets (`cl_SideSet.cpp:34-41`); a `Facet` deletes its wrapped element (`cl_Facet.cpp:29-32`). Each `delete` appears exactly once, in the destructor of the designated owner — `~Mesh()` explicitly does *not* delete elements and facets, with comments marking whose job it is (`cl_Mesh.cpp:256-266`).

**The pool pattern, done right.** The one placement-new pool in the tree is `ShiftRegister<T>`: `std::malloc` for the buffer (`cl_ShiftRegister.hpp:493`), placement-new per slot (`:109`), explicit `->~T()` per slot (`:126`), single `std::free` (`:250-255`). The rule it demonstrates: **objects constructed by placement new are destroyed by an explicit destructor call plus `free` on the pool — calling `delete` on a pool member is undefined behavior.** If a pooled/SoA element store is ever pursued (a legitimate post-release discussion), this is the discipline it must follow.

**Rule of Five.** Any class owning raw memory declares or deletes all five special members. The tree's raw owners are compliant and serve as the reference examples:
- `ShiftRegister`: all five user-defined, moves `noexcept` (`cl_ShiftRegister.hpp:169,186,203,222,250`)
- `DynamicBitset`: all five user-defined (`cl_DynamicBitset.cpp:43,56,418,457,74`)
- `StringList`: all four copy/move members explicitly deleted (`cl_StringList.hpp:39-43`) — the right call for a C-interop buffer that has no business being copied

---

### Where Smart Pointers Actually Appear

The real census (2026-08-08) is small and specific — and it is *not* "I/O and unit tests" (there are zero smart pointers in `tests/` and zero in `src/io/`):

- **Factory-scope ownership handles:** `MaxwellFactory` and `ThermalFactory` hold and return `std::shared_ptr<Kernel>` and `std::shared_ptr<Controller>` (`cl_MaxwellFactory.hpp:66,128,133`; `cl_ThermalFactory.hpp:39,64`). These are created once at setup; the refcount is touched a handful of times per run, so `shared_ptr`'s atomic overhead is irrelevant here. What matters is the entity rule: the *Kernel* internally manages its members with the hierarchical raw-pointer discipline (`cl_FEM_Kernel.cpp:135-170`); the smart pointer is just the top-level lifetime handle.
- **Scratch with a guaranteed release:** `sprint()` uses a `unique_ptr<char[]>` for its formatting buffer (`fn_sprint.hpp:50`); the startup banner wraps a `popen` handle in a `shared_ptr<FILE>` with a `pclose` deleter (`banner.cpp:50`).
- **Dev drivers:** an earlier edition cited `corctest.cpp:44` here as an example of a throwaway driver using smart pointers freely. That file no longer exists, and the surviving drivers under `src/fem/postproc/` (`normaltest.cpp`, `pentatest.cpp`) use raw `new`/`delete` — so the tree currently has no example of this pattern to point at.

**Rule of thumb (unchanged in spirit, corrected in letter):** if the code runs inside the assembly/solve/communication loops, use the deliberate-allocation patterns above. If it runs once at setup or teardown, a smart pointer as a lifetime handle is acceptable — and for shared top-level objects, already established practice.

---

## Performance-First Design Patterns

### 1. Chunking for Large Messages

**Problem:** MPI counts are `int`-typed in the classic API; large mesh payloads exceed what a single message should carry.

**Solution:** the comm layer transparently splits every large transfer into chunks of `gMaxCommChunkLength = 64 * 1024` **elements of `T`** — not bytes (`src/comm/commtools.hpp:30`). `comm_split()` computes the chunk list from the element count (`commtools.cpp:80-92`), and each chunk goes out as an `MPI_Isend` of `count` elements typed via `comm_type<T>()` (`commtools.hpp:330-333`), followed by a single `MPI_Waitall`.

Element types are primitives (or `std::complex`): `comm_type<T>()` is specialized only for arithmetic types and errors on anything else (`src/comm/commtypes.hpp:38-204`). There is **no generic object serialization** — a `Cell` of arbitrary structs cannot be sent without adding a specialization. Strings are flattened to `Vector<char>` (`commtools.cpp:147-225`).

**Why chunking, when MPI-4 has large counts?** MPI 4.0's `MPI_Send_c` family accepts `MPI_Count` payloads natively. Chunking is retained because production clusters still run older MPI stacks, and the chunk loop costs nothing measurable next to the transfer itself. When MPI-4 is the deployment floor, the chunk layer can collapse to a passthrough without touching any caller.

**Distributing data from root to all ranks: payload size dictates the API.**

| Payload | Function |
|---------|----------|
| Single scalar | `broadcast(T&)` |
| Small fixed-size payload (~≤ 8 entries, e.g., a properties tuple, a 3×3 matrix) | `broadcast(...)` |
| Large or variable-size `Vector<T>` or `Cell<T>` (mesh quantities, ephemeris columns) | `share` / `receive` (`commtools.hpp:1548`, `:1648`) |
| `Matrix<T>` of any size | `broadcast(...)` — no `share( Matrix )` overload exists; for very large payloads use a manual `send` / `receive` loop |

`broadcast` is collective: every rank calls it with the same arguments, and the implementation packs the payload into a single unchunked `MPI_Ibcast`. Fine for small bounded payloads; observed to fail on large ones. For large variable-size `Vector<T>` or `Cell<T>` data, use the asymmetric `share`/`receive` pair, which chunks:

```cpp
if ( comm_rank() == 0 )
{
    // ... fill vec ...
    belfem::share( vec );        // root only — chunks the message
}
else
{
    belfem::receive( vec );      // non-root only — receives chunks
}
```

The `if/else` rank guard is required: unlike `broadcast`, `share` and `receive` are **NOT collective**. Calling either on the wrong rank is a deadlock. There is no `share` overload for scalars — scalar `broadcast(T&)` is always the right choice for single-element messages.

**Matrix transfers ship the raw buffer.** `send(Matrix<T>&)` transmits `n_rows`, `n_cols`, and the transfer length `spacing() * n_cols`, then that many elements of the raw buffer — deliberately including any backend padding, so the padded layout round-trips without disassembly (both ranks run the same binary, so the padding is identical for a given shape). Never `capacity()`: it can be stale-large after a shrink. The receive paths assert the transfer length fits the local buffer.

---

### 2. Alignment Reality

An earlier revision of this document showed BELFEM calling `aligned_alloc(64, ...)` for its vectors. **That is not how it works.** `src/` contains no call to `aligned_alloc` or `posix_memalign` at all; alignment comes entirely from the linear-algebra backends:

- **Blaze** allocates through its `AlignedAllocator`, which calls `posix_memalign` with the alignment the SIMD build requires: 32 bytes under AVX2, 64 bytes only under AVX-512. Its error handling checks the *return code* (the correct idiom — `posix_memalign` does not null the pointer on failure).
- **Armadillo** calls `posix_memalign` with 32-byte alignment for allocations ≥ 1 KiB and 16 bytes below, and deliberately caps below 64.
- BELFEM's own raw buffers (`StringList`, `DynamicBitset`, `ShiftRegister`, mesh-entity arrays) use plain `malloc` — default alignment, which is sufficient for pointer and integer arrays.

So the honest statement is: **numerical data is SIMD-aligned to whatever width the backend's build targets; nothing in BELFEM guarantees 64 bytes.** BLAS kernels handle the residual peel/tail cases either way; the alignment mostly determines whether the vectorized core runs aligned loads. *(The cache-line-split penalty for misaligned streaming access is a literature heuristic — Drepper §6.2, Agner Fog's optimization manuals — not a BELFEM measurement.)*

**If you ever hand-allocate aligned memory** (currently no site does), the rules are:
- C11 `aligned_alloc(align, size)` requires `size` to be a **multiple of `align`** — round up, or it is UB per C11 §7.22.3.1 (glibc happens to tolerate it; strict implementations do not)
- `posix_memalign` reports failure via its **return code** and leaves the pointer indeterminate — check the return value, never the pointer
- Either way the buffer is released with `free()`, never `delete[]`

---

### 3. Contiguous Storage for MPI

**Problem:** `std::vector<std::vector<T>>` is not contiguous — each row is a separate allocation, so it cannot be sent as one buffer.

**Solution:** flatten to a single contiguous object.

```cpp
// Bad: Non-contiguous — cannot MPI_Send as one buffer
std::vector<std::vector<real>> tRows(n_rows);

// Good: Contiguous — one buffer, one (chunked) transfer
Matrix<real> tMatrix(n_rows, n_cols);
send( tMatrix, aTarget );   // comm layer handles sizes, padding, chunking
```

This is also the data-oriented-design argument (Acton): the machine rewards arrays of plain data laid out for the access pattern, not graphs of small heap objects.

---

### 4. Zero-Overhead Abstractions

**Principle:** wrappers should compile to the same machine code as the raw implementation.

**Example:** `Cell<T>`'s `operator()` is an inlined `std::vector` access guarded by `BELFEM_ASSERT` (`cl_Cell.hpp:148-168`). Under `USE_DEBUG=ON` (`-Og -g`, `DEBUG` defined) the bounds check runs; under `USE_DEBUG=OFF` (`NDEBUG`) the assert expands to nothing and the accessor is a raw indexed load.

**Verification:** compare assembly with `objdump -d` or godbolt — the release accessor is indistinguishable from a raw array access.

**No temporary `Vector`/`Matrix` in frequently-called member functions.** Any method expected to run repeatedly at runtime allocates its scratch as members, sized once in the constructor (the caller may loop over it millions of times). Expression-template assignments into an already-sized member evaluate in place and are fine. Temporaries remain acceptable in constructors and one-off setup paths.

---

## Error Handling Strategy

BELFEM uses a **two-tier macro system** plus a third, macro-free category for algorithmic failures.

### The Semantic Rule (primary)

- **`BELFEM_ASSERT` — "if this fires, BELFEM has a bug."** Logic errors, bounds violations, broken invariants, precondition breaches by internal callers. Compiled out under `NDEBUG`.
- **`BELFEM_ERROR` — "this can fail in a correct program."** Bad input files, missing environment, failed I/O, failed allocation, MPI errors. Always active in every build.

```cpp
// Debug-only logic check (compiled out in release)
BELFEM_ASSERT( aIndex < mCell.size(), "Index %lu exceeds size %lu", aIndex, mCell.size() );

// Always-active runtime check
BELFEM_ERROR( tFile.is_open(), "Failed to open: %s", aFilename.c_str() );
```

### The Frequency Rule (placement, 2026-08-08)

The check's runtime cost is what justifies compiling it out — so *where* a requirement is validated follows how often the code runs:

- **Hot paths** (container accessors, per-integration-point math, assembly loops, communication inner loops): `BELFEM_ASSERT` wherever reasonable. These run millions of times; an always-active check has a real cost.
- **Setup and initialization** (factories, mesh enrichment, dof wiring, input parsing — anything that runs once per run): be generous with `BELFEM_ERROR`. The check costs nanoseconds once, and a debug-only assert here means a release build silently continues into a state that crashes — or worse, computes garbage — far from the cause. No per-case deliberation: in once-per-run code, default to the always-active tier.

A runtime requirement too expensive to check per-call is validated **once at construction** with `BELFEM_ERROR` — never demoted to debug-only just because the accessor is hot. The same precondition can legitimately be `BELFEM_ERROR` in a factory and `BELFEM_ASSERT` in the accessor it protects.

### Expected Algorithmic Failure (third category — not an abort)

Solver non-convergence, line-search failure, a rejected timestep: these are **not** errors in either macro's sense. They are anticipated outcomes of numerical algorithms, and they return a status or trigger the retry path — timestep reduction, relaxation, controller fallback — exactly as the nonlinear controller already practices. Reaching for `BELFEM_ERROR("did not converge")` inside an algorithm that has a retry policy above it is a design error: it converts a recoverable state into a run abort.

### What the Macros Actually Do

Both funnel through `belfem::assert::belfem_assert(...)`, which calls `assert::error`. That function prints the error box and then consults `throw_on_error()`: **true throws** the `std::runtime_error`, **false calls `error_abort()`**, which delegates to the comm module's `comm_abort` wrapper (declaration and full contract: `src/comm/cl_Communicator.hpp:235-258`; definition: `src/comm/commtools.cpp:67`). That wrapper aborts the job with `MPI_Abort(MPI_COMM_WORLD, 1)` — `MPI_COMM_WORLD` rather than `gComm.world()`, guarded by `MPI_Initialized`/`MPI_Finalized` — and otherwise falls through to `std::abort()`. `BELFEM_ASSERT` additionally compiles to nothing when assertions are inactive; `BELFEM_ERROR` is always expanded.

The flag is initialized to the build's compile-time behavior — throwing where assertions are active, aborting otherwise — so a debug run still throws and a production run still aborts, at any rank count, with no caller involvement. **A debug run throws in parallel too**, deliberately: that is what keeps a failed check inspectable while the other ranks are still alive. `MPI_Abort` is the production reaction. It is a **test hook, not a configuration knob**: a production run must keep the abort, because a throw that escapes `main` terminates one rank and leaves its peers blocked in a collective, and nothing in `src/` catches BELFEM errors. The state deliberately lives in `assert.cpp` rather than as a header static, because `assert::error` is a function template instantiated in the *calling* translation unit — a per-TU copy would mean a test binary could not change the reaction of an already-compiled library.

The reason the hook exists: `BELFEM_ERROR` covers the failures that survive into release (file I/O, MPI, allocation, unsupported configuration), so those paths are worth unit-testing in the shipped configuration. Each test `main` calls `set_throw_on_error( true )` after `gComm.init`, which makes `EXPECT_THROW` work under `NDEBUG` instead of aborting the whole test binary. Two consequences for test authors. First, guard assertion-tier tests with `#if BELFEM_ASSERTIONS_ACTIVE`, the macro exported by `assert.hpp` — not a restated `#ifndef NDEBUG`, which is not equivalent when both `NDEBUG` and `DEBUG` are defined. Second, the two predicates are distinct: a release binary in throw mode can catch a `BELFEM_ERROR` while `BELFEM_ASSERT` still expands to nothing. What such a test verifies is that the check fires with the right message, **not** that the process aborts; release and test runs take different reactions by construction.

Note precisely what drives the default abort-vs-throw split: the **`NDEBUG` macro**, not the exception model. Exceptions are currently enabled in every build configuration (see Core Principles); an uncaught `std::bad_alloc` from a failed `std::vector` reallocation in a debug build unwinds and terminates with a diagnosable message rather than aborting. Under the aspirational `-fno-exceptions` regime, allocation failure would become an immediate abort inside the library — one more reason that switch is a deliberate future decision, not a checkbox.

---

## Thread Safety and MPI

### Thread Safety Philosophy

**Design decision:** BELFEM is **deliberately not thread-safe** internally.

**Rationale:**
1. BELFEM parallelizes with **MPI** (distributed memory), one single-threaded process per rank
2. Internal mutexes add latency to **every** operation, even single-threaded runs
3. Users needing OpenMP can protect calls externally (`#pragma omp critical`)

**What this means:**
- `Logger`, `Progressbar`, `Profiler` — **not thread-safe**
- `rand()` — **not thread-safe** (use per-thread RNG instances)
- `gLog`, `gComm`, globals — **require external synchronization**

**What IS thread-safe:**
- `Timer`, `Hash` — each thread uses its own instance
- `const` functions on containers

---

### One Parallel Programming Model

**Design decision:** BELFEM deliberately supports a **single parallel programming model** — pure MPI, one rank per core, application code single-threaded. Hybrid MPI+X is a cost we do not pay until a measured problem justifies it.

This is a *scope decision*, not a claim that hybrid parallelism cannot work. What the single model buys:

- One thing to debug, profile, and reason about; no rank×thread 2-D tuning space
- No thread-safety obligations on any BELFEM data structure (see above — and the entire framework is built on that exemption)
- No oversubscription hazards from `n_ranks × n_threads` exceeding the node
- Memory ownership and data locality follow rank boundaries exactly

**The facts, and the second one is easy to miss:**

1. **BELFEM's own hand-written OpenMP is off by default.** Three Fortran kernels
(`src/sparse/splinalg.f90`, `arpacktools.f90`, `parpacktools.f90`) carry `!$omp` directives, gated
on the `BELFEM_OMP` define, which is OFF (`USE_BELFEM_OPENMP`, see `doc/parallel_execution.md`).

2. **BELFEM's compiled objects nevertheless contain OpenMP, from the matrix backend.** Both backends
are header-only and parallelize their own expression evaluation whenever OpenMP is available —
Armadillo because BELFEM asks it to (`src/linalg/armadillo/armadillo.hpp`, `#define ARMA_USE_OPENMP`
under `#ifdef OMP`), Blaze on its own (`blaze/system/SMP.h`: OpenMP mode whenever `_OPENMP` is
defined, and `blaze/config/SMP.h` defaults shared-memory parallelization on). So the pragmas are
compiled *into BELFEM translation units*, not confined to a separate library. Measured on a default
Darwin release build (Blaze, `USE_OPENMP=ON`): **129 of 339 objects reference `GOMP_*`.** This is
platform-independent — the Linux default (Armadillo) enables it too.

The practical consequence is about *reading a stack trace*: **a `gomp_thread_start` frame with
BELFEM code above it does not imply a BELFEM pragma.** It is not a second source of the
stack-overflow class described in `doc/parallel_execution.md` — that needs an OpenMP *array*
reduction, and neither backend uses a reduction of any kind (Blaze partitions the destination,
Armadillo uses `critical`/`atomic`). Third-party solvers (STRUMPACK, MKL, threaded BLAS, PARDISO) also spawn workers.
BELFEM's own OpenMP *touchpoints* are only queries about their budget — `omp_get_max_threads` for the
banner and the oversubscription warning below, and PARDISO's thread parameter — and those stay on the
plain `OMP` define by design.

**The sanctioned exception — solver-internal threading.** Third-party solvers thread internally (STRUMPACK/SLATE, threaded BLAS). BELFEM accommodates them without becoming multithreaded itself:
- STRUMPACK builds initialize MPI with `MPI_Init_thread(..., MPI_THREAD_MULTIPLE, ...)` — requested for STRUMPACK/SLATE compatibility, warning-only on downgrade (`src/comm/cl_Communicator.cpp:100-131,296-309`). All of BELFEM's own MPI traffic still happens from the main thread.
- The solver wrapper warns (`hatch_turtle()`, `src/sparse/cl_SolverWrapper.cpp:195-245`) when the OpenMP threads requested by the ranks sharing a node would oversubscribe its cores.

**The honest counterargument, kept on the table:** pure MPI replicates per-rank state. On high-core-count nodes, rank-level memory replication may eventually force fewer-ranks-with-threads — that decision will be made from **measurements** on the memory-bound configuration, not from folklore in either direction. (Earlier revisions carried specific efficiency percentages for MPI vs. hybrid; they were unmeasured and are gone. See [Benchmarks We Owe Ourselves](#benchmarks-we-owe-ourselves).)

**If you need threading around BELFEM today:** serialize the calls —

```cpp
#pragma omp parallel
{
    #pragma omp critical
    {
        solver->step();   // BELFEM is NOT thread-safe — must serialize
    }
}
```

or better, use finer MPI decomposition (more ranks, one thread each).

---

### MPI Design Patterns

**1. Initialization** — plain `MPI_Init` by default; `MPI_THREAD_MULTIPLE` only under STRUMPACK (`cl_Communicator.cpp:100-131`); the `Communicator` constructor then queries rank, size, and `MPI_TAG_UB`.

**2. The rank ceiling (documented limit).** Message tags encode the communicating rank pair: `comm_tag(aSource, aTarget) = 2*(tMax*gComm.size()+tMin) % gComm.max_tag()` (`src/comm/commtools.cpp:64-73`), with base tag for the size header and base+1 for the payload. From `MPI_TAG_UB` the constructor derives and **enforces** a hard rank ceiling `tMaxRank = 0.5*(sqrt(2*mMaxTag+1)+1)` (`src/comm/cl_Communicator.cpp:169-187`): about **32,768 ranks** at the most generous `MPI_TAG_UB` (`INT_MAX`), less on implementations reporting a smaller tag bound. Until the tag scheme is redesigned, that is BELFEM's scaling limit, and this document makes no core-count claims beyond it. The redesign is tractable: the scheme is fully encapsulated in `commtools.hpp` (no public API takes a tag; tags are recomputed symmetrically on both sides), so per-pair sequence counters or per-phase communicators are a comm-layer-local change — tracked as a proposal in the appendix.

**3. Rank-specific execution:**
```cpp
if ( gComm.rank() == 0 )
{
    print_banner( "MyApp" );
}
```

**4. Global variable synchronization:**
```cpp
if ( gComm.rank() == 0 )
{
    gTbulk = read_from_input_file();
}
broadcast( gTbulk );  // all ranks now have the same value
```

**5. Profiler generates rank-specific files** (`profile.4.2.log` for size=4, rank=2) so each rank's profile is analyzed separately.

---

## Benchmarks We Owe Ourselves

Earlier revisions of this document carried performance numbers (cache-penalty factors, alignment throughput gains, MPI-vs-hybrid efficiency percentages) that no BELFEM measurement supports. They have been removed or replaced with cited literature heuristics. The following measurements would let future revisions state them as facts:

| Claim to substantiate | Measurement needed |
|----------------------|--------------------|
| Cost of wrong-order matrix iteration | Micro-benchmark: row-inner vs column-inner fill on release build, sizes spanning L1/L2/LLC |
| Value of SIMD-aligned numerical buffers | Same BLAS-heavy kernel on Blaze-AVX2 (32 B) vs AVX-512 (64 B) builds |
| Pure-MPI parallel efficiency | Strong+weak scaling study of a representative h-φ problem; report efficiency vs. rank count up to the tag-scheme ceiling |
| `MPI_THREAD_MULTIPLE` init cost | Same run, STRUMPACK build vs plain build, on the production MPI stack |
| Release-vs-debug check overhead | `make check` workload timed under `USE_DEBUG=ON` vs `OFF` |

Until a `benchmarks/` suite exists, treat any performance number outside this document's cited heuristics as unverified.

---

## Summary Checklist

### For New Contributors

| ✅ | What to Remember |
|----|------------------|
| **Naming** | Prefix variables (`a` = argument, `t` = temporary, `m` = member, `g` = global). Relax for pure math. |
| **Containers** | `Cell<T>` for arrays, `Vector`/`Matrix` for linear algebra, `DynamicBitset` for large boolean sets. Size accessor: `size()` on `Cell`/`DynamicBitset`, `length()` on `Vector`/`Matrix`. |
| **Memory** | Allocate deliberately at setup. Plain `T*` = non-owning borrow; owners delete in their destructors; Rule of Five on every raw-owning class. Smart pointers OK as setup-scope lifetime handles. |
| **Performance** | Profile before optimizing. No hidden allocations in loops; scratch objects live as members. |
| **Error Handling** | `BELFEM_ASSERT` = "BELFEM has a bug" (debug-only); `BELFEM_ERROR` = "can fail in a correct program" (always). Once-per-run code defaults to `BELFEM_ERROR`. Convergence failure returns status, never aborts. |
| **Thread Safety** | Core is not thread-safe. One parallel model: pure MPI. Solver-internal threading is the sanctioned exception. |
| **MPI** | Broadcast small payloads; `share`/`receive` (with rank guard!) for large `Vector<T>`. Chunking is in elements of `T` (`gMaxCommChunkLength` = 65 536 elements). |
| **Alignment** | Comes from the backend (`posix_memalign`, SIMD-width dependent). Do not hand-roll aligned allocation without the C11 size-rounding and return-code rules. |
| **Builds** | `USE_DEBUG=OFF` = `-O2` + `NDEBUG`; `USE_DEBUG=ON` = `-Og` + assertions. The default follows the SCLS flavor: `OFF` with no `$SCLS` or the `gcc`/`mkl` flavors, `ON` for the `debug` flavor. Exceptions are currently enabled in all builds; don't rely on `throw` in library code anyway. |

---

### Quick Reference: When to Use What

**Use `Cell<T>` when:**
- Storing collections of mesh entities (nodes, elements, faces)
- Building lists of IDs, indices, pointers
- General-purpose dynamic arrays

**Use `Vector<T>` / `Matrix<T>` when:**
- Performing linear algebra (matrix-vector multiply, solve, eigenvalues)
- Storing DOF values, solution fields, coordinate arrays
- Interfacing with BLAS/LAPACK

**Use `DynamicBitset` when:**
- Tracking large boolean flags (>1000s of values)
- Topology operations (node activity, boundary marking)
- Memory-critical applications (1 bit vs 1 byte per bool)

**Use `ShiftRegister<T>` when:**
- Implementing time-stepping with history (BDF methods)
- Need to revert to previous state (adaptive stepping)
- Fixed-depth buffer with no reallocations

**Use manual memory when:**
- Owning a raw C-interop buffer (`char**` for C libraries, bit storage)
- Building fixed-size per-entity arrays at construction time
- Interfacing with C libraries (MPI, BLAS, Fortran)

**Use smart pointers when:**
- Holding a top-level object's lifetime at setup scope (factory-owned kernels/controllers — established practice)
- Scratch buffers with a guaranteed release on all paths (`sprint`, `popen` handles)
- One-off allocations outside hot paths

---

## Further Reading

### Internal Documentation
- **`CLAUDE.md`** — Project-wide conventions and build system
- **`doc/core_module_overview.md`** — Core utilities (Logger, Timer, typedefs)
- **`src/containers/doc/container_usage_guide.md`** — Container APIs and examples
- **`src/comm/doc/comm_module_overview.md`** — MPI communication patterns

### External References
- **Ousterhout, _A Philosophy of Software Design_** — deep modules (`Cell`, the comm layer, the backend wrappers are exactly this); "define errors out of existence" is the `data()`-valid-when-empty principle
- **Lakos, _Large-Scale C++ Software Design_** — physical design and levelization for a 300k-line codebase
- **Drepper, "What Every Programmer Should Know About Memory"** — the cited source for cache-line and locality heuristics in this document
- **Agner Fog, *Optimization Manuals* (agner.org)** — instruction-level and alignment heuristics
- **C++ Core Guidelines, performance section (Per.\*, ES.84)** — the modern-C++ mainstream *agrees* with no-hidden-allocations and no-needless-heap; these rules are allies, not opponents
- **Acton, "Data-Oriented Design and C++" (CppCon 2014)** — the citation for the contiguity and flat-storage arguments
- **Meyers, *Effective C++* (Items 13-17)** — resource management patterns
- **MPI Standard (Chapter 3)** — communication semantics, buffer requirements

**On *Clean Code*:** BELFEM borrows exactly three things — the boy-scout rule, intention-revealing names, and comments that explain *why* (practiced here as equation citations). Its small-method, OO-heavy decomposition style is deliberately **not** adopted: scattering a numerical kernel across a dozen tiny methods destroys the correspondence between code and published equations that Mathematical Readability exists to protect.

---

## Appendix: Source Reconciliation

This document has been through two source-reconciliation passes (2026-06-14 B1-scoped; 2026-08-08 full). The complete per-claim ledgers — verdict, old text, new text, citation — live in `tmp/whitepaper/coding_philosophy_CHANGES.md`. Summary of what the 2026-08-08 pass changed:

**Corrected against source:**
- Build flags: release is `-O2`+`NDEBUG` via `USE_DEBUG=OFF` (never `-O3`); **no build sets `-fno-exceptions`**, and `assert.hpp:191` would not compile under it — the exception-free build is design intent, not fact
- Alignment: no `aligned_alloc`/`posix_memalign` call exists in `src/`; alignment is backend-provided and SIMD-width dependent (≤32 B under Armadillo and Blaze/AVX2); the former hand-rolled 64-byte `Vector` example was defective (unrounded C11 size; discarded `posix_memalign` return) and is gone
- Mesh entity allocation is per-object `new`/`delete` with hierarchical ownership; the former pool + placement-new example did not describe BELFEM and is gone (`ShiftRegister` is the one real placement-new pool)
- Smart pointers: real census is factory-held `shared_ptr<Kernel>`/`shared_ptr<Controller>`, `sprint()` scratch, `popen` handles, dev drivers — not "I/O and unit tests" (both of which have zero)
- `Cell` examples now use the real API (`size()`, free-function `unique(...)`; there is no `fill()`/`length()`/member-`unique()` on `Cell`)
- `data()`-when-empty: well-defined, *possibly null*; backends disagree on null-ness; never null-test for emptiness
- Testing posture: `USE_TEST` defaults ON since 2026-08-14 (was OFF when this ledger entry was first written). The "no CI" half of this entry expired 2026-08-30, when a GitLab server began running the full suite nightly; no sanitizers and manual-only Valgrind still stand, as does the absence of a committed pipeline definition

**Removed (unverifiable performance figures):** the 8×/2-3×/30-40%/10-30%/>95%/70-85%/95-vs-75 numbers and the ">100,000 cores" claim (superseded by the enforced ~32,768-rank tag ceiling, now documented as a limit). What remains cites Drepper/Agner Fog as literature heuristics.

**Defects surfaced by the 2026-08-08 verification — all fixed 2026-08-09 (details in the CHANGES ledger; three-AI jury round in `tmp/ai_exchange/review_matrix_mpi_and_test_gating.md`):**
1. Blaze-backend Matrix MPI transfer overflow on shrunk matrices → fixed: all five transfer paths send `spacing()*n_cols` instead of `capacity()`, receive paths guard capacity, unchunked `broadcast` bounds its count by `INT_MAX`
2. Copyable `Logger` owning a `FILE*` → fixed: copy/move deleted, `fopen` checked
3. `DynamicBitset`/`StringList` minors → fixed: `mIndex` copied in copy-assign, `noexcept` move ctor, checked allocations (zero-size bitsets skip allocation), `StringList::push` bound always-active
4. `make check` not building `test_core` → fixed: `core` added to both check target dependency lists

Pending executable gate: an MPI `make check` run including the shrink-then-transfer regression test (source-trace evidence only until then).

**Possible future directions (explicitly out of scope before 1.0):** a pooled/SoA element store following the `ShiftRegister` placement-new discipline; the tag-scheme redesign; snippet-extraction CI so this document's code blocks can never drift from source again.

---

**Maintained by Christian Messe, Lawrence Berkeley National Laboratory.**

*Any deviation from these principles should be documented in a GitHub issue with performance justification.*
