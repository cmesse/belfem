# BELFEM Sparse Matrix Module - Usage Guide {#sparse_sparse_usage_guide}

**Date:** 2026-01-16
**Module:** sparse
**Purpose:** Comprehensive guide to BELFEM's sparse matrix and solver API
**Revision:** 2026-01-16 - Initial version

---

## Overview

The `src/sparse` module provides **sparse matrix storage** and **direct/iterative solver** interfaces for BELFEM. It wraps multiple high-performance solver libraries:

- **UMFPACK** (SuiteSparse): Sequential sparse direct solver
- **MUMPS**: Parallel sparse direct solver (MPI)
- **STRUMPACK**: Parallel sparse direct/iterative solver (MPI, low-rank compression)
- **PARDISO**: Shared-memory parallel direct solver (Intel MKL)
- **PETSc**: Scalable iterative solver framework (MPI, Krylov methods)

The abstraction allows switching solver backends at **compile time** without changing user code, enabling portability across different HPC systems.

### Key Features

1. **Multiple Solver Backends**: Compile-time selection via CMake flags
2. **Sparse Matrix Formats**: CSR (Compressed Sparse Row) and CSC (Compressed Sparse Column)
3. **MPI Support**: Distributed sparse matrices and parallel solvers
4. **Flexible Indexing**: C++ (0-based) or Fortran (1-based) indexing
5. **Solver Parameters**: Unified interface for reordering, preconditioners, tolerances
6. **Graph-Based Construction**: Build sparse matrices from connectivity graphs

> **Critical: Solver Selection**
> Solver **availability** is determined at **compile time** via the CMake
> options (the `BELFEM_*` spellings are the generated compile definitions,
> not cache entries):
> - `USE_SUITESPARSE=ON` → UMFPACK available (default OFF)
> - `USE_MUMPS=ON` → MUMPS available (default ON)
> - `USE_STRUMPACK=ON` → STRUMPACK available (preferred, default ON)
> - `USE_PARDISO=ON` → PARDISO available (default OFF; ON under the SCLS `mkl` flavor, in parity with `USE_MKL`)
> - `USE_PETSC=ON` → PETSc available (default ON)
>
> Multiple solvers can be enabled in a single build. **Runtime selection** allows switching
> between available solvers:
> ```cpp
> Solver solver1(SolverType::UMFPACK);  // If BELFEM_SUITESPARSE enabled
> Solver solver2(SolverType::MUMPS);    // If BELFEM_MUMPS enabled
> ```
> Use the default solver (automatically selected at compile time) or specify explicitly.

---

## Common Pitfalls

### 1. **Indexing Base Mismatch**

```cpp
SpMatrix A(graph, SpMatrixType::CSR);
// Default: C++ indexing (0-based)
int base = A.indexing_base();  // Returns 0

// MUMPS requires Fortran indexing (1-based)
A.set_indexing_base(SpMatrixIndexingBase::Fortran);
// Now base is 1
```

**Solution:** Check solver requirements. MUMPS needs Fortran indexing; UMFPACK uses C++ indexing.

### 2. **CSR vs. CSC Format**

```cpp
// SpMatrix defaults to CSC if not specified
SpMatrix A(graph);  // Uses CSC by default

// UMFPACK prefers CSC
SpMatrix A_csc(graph, SpMatrixType::CSC);

// MUMPS, STRUMPACK, PARDISO, and PETSc prefer CSR
SpMatrix A_csr(graph, SpMatrixType::CSR);
```

**Solution:** Always specify format explicitly, or use `preferred_matrix_format()`:
```cpp
SpMatrix A(graph, preferred_matrix_format(SolverType::MUMPS));
```

> **Note:** SpMatrix default constructor uses **CSC** format. For solvers that prefer CSR
> (MUMPS, STRUMPACK, PARDISO, PETSc), always specify `SpMatrixType::CSR` explicitly to avoid performance loss.

### 3. **COO Indices for MUMPS**

MUMPS needs explicit row/column index arrays, but **the solver creates them itself** — 
`SolverMUMPS` calls `create_coo_indices()` on the matrix before it solves
(`cl_SolverMUMPS.cpp:716-725`). A caller does not have to:

```cpp
SpMatrix A(graph, SpMatrixType::CSR);
// no create_coo_indices() needed -- SolverMUMPS does it

// If you built them yourself for another reason, this releases them again:
A.free_coo_indices();
```

### 4. **Overwriting Matrix During Solve**

```cpp
SpMatrix A(graph, SpMatrixType::CSR);
Vector<real> x(n), b(n);
// ... fill A and b ...

Solver solver(SolverType::MUMPS);
solver.solve(A, x, b);  // A's values are untouched; wrappers only toggle the indexing base and may add COO index arrays

// If you need A again:
SpMatrix A_copy;
A_copy = A;  // deep copy via assignment (the copy constructor is deleted)
```

**Solution:** None needed for the values: they survive the solve. Only the indexing base may
have been flipped and restored; copy only if you rely on a specific base afterwards.

### 5. **MPI Distribution**

A whole `SpMatrix` object is not sent over MPI, and there is no `distribute( SpMatrix* )`
overload. The matrix lives on the master rank only; the solver wrappers handle
distribution internally:

- **MUMPS** takes the host-centralized matrix from rank 0 and distributes inside
  the library.
- **STRUMPACK/PETSc** use `sparse::DistMatrix`, which extracts the pointer,
  index, and value arrays on rank 0 and scatters raw array slices to workers.

```cpp
// Master rank builds and assembles; workers pass an empty matrix
SpMatrix A;
if (gComm.rank() == 0) {
    A = SpMatrix(graph, SpMatrixType::CSR);
    // ... assemble A ...
}

// solve() is collective; non-root ranks never read A's content
solver.solve(A, x, b);
```

**Solution:** Build the matrix on rank 0 and let the solver wrapper distribute.

---

## Core Types

### SpMatrix - Sparse Matrix

**Files:** `cl_SpMatrix.hpp`

#### Description

Compressed sparse matrix storage (CSR or CSC format) with element access, matrix-vector multiplication, and I/O capabilities.

#### Sparse Matrix Formats

**CSR (Compressed Sparse Row):**
```
Row-major storage
Pointers[i]: Start of row i in Values array
Indices[j]: Column index of Values[j]
```

**CSC (Compressed Sparse Column):**
```
Column-major storage
Pointers[j]: Start of column j in Values array
Indices[i]: Row index of Values[i]
```

**COO (Coordinate) format (for MUMPS):**
```
Explicit row and column index arrays
Rows[k]: Row index of Values[k]
Cols[k]: Column index of Values[k]
```

#### Construction

```cpp
// From connectivity graph (recommended)
Graph graph;
// ... populate graph with connectivity ...
SpMatrix A(graph, SpMatrixType::CSR);

// From graph with explicit size
SpMatrix A(graph, SpMatrixType::CSC, num_rows, num_cols);

// From existing compressed format
SpMatrix A(SpMatrixType::CSR, n_rows, n_cols, nnz, indices, pointers);

// From dense matrix (for testing)
Matrix<real> dense = {{4, 0, 1},
                      {0, 2, 0},
                      {1, 0, 3}};
SpMatrix A(dense, SpMatrixType::CSR);

// Load from HDF5 file
SpMatrix A("matrix.h5", "StiffnessMatrix");

// Empty matrix
SpMatrix A;
```

#### Element Access

```cpp
SpMatrix A(graph, SpMatrixType::CSR);

// Bounds-checked access (debug mode)
real value = A(0, 0);      // Read element (0,0)
A(2, 3) = 1.5;             // Write element (2,3)

// Read returns zero if element is structurally zero
real zero = A(1, 5);       // Returns 0.0 if not in sparsity pattern

// Write requires element to exist in sparsity pattern
// A(1, 5) = 2.0;          // ERROR if (1,5) is structurally zero

// Raw data access
real* values = A.data();
int_t* indices = A.indices();  // Row indices (CSC) or column indices (CSR)
int_t* pointers = A.pointers();
```

#### Size and Properties

```cpp
index_t rows = A.n_rows();
index_t cols = A.n_cols();
index_t nnz = A.number_of_nonzeros();
SpMatrixType type = A.type();  // CSR or CSC

size_t mem = A.memory();  // Memory in bytes
bool have_coo = A.have_coo_indices();
```

#### Matrix-Vector Multiplication

```cpp
SpMatrix A(graph, SpMatrixType::CSR);
Vector<real> x(n), y(n);

// Simple multiplication: y = A * x
A.multiply(x, y);

// General form: y = alpha * A * x + beta * y
A.multiply(x, y, 2.0, 1.0);  // y = 2*A*x + y

// Transposed multiplication: y = alpha * A^T * x + beta * y
A.multiply(x, y, 1.0, 0.0, true);  // y = A^T * x

// Operator overload
Vector<real> y = A * x;  // Allocates new vector
```

#### Utilities

```cpp
// Fill all values
A.fill(0.0);  // Zero out matrix

// Transpose (in-place)
A.transpose();  // Converts CSR ↔ CSC

// Change indexing base
A.set_indexing_base(SpMatrixIndexingBase::Fortran);  // 1-based
A.set_indexing_base(SpMatrixIndexingBase::Cpp);      // 0-based

// MUMPS compatibility
A.create_coo_indices();  // Create row/col index arrays
A.free_coo_indices();    // Free memory

// Printing (for debugging)
A.print("MyMatrix");     // Print in matrix form
A.print2("MyMatrix");    // Print compressed indices
```

#### Saving and Loading

```cpp
// Save to HDF5
A.save("matrix.h5", "StiffnessMatrix", FileMode::NEW);

// Load from HDF5
SpMatrix B;
B.load("matrix.h5", "StiffnessMatrix");

// Save to existing HDF5 group
hid_t group = ...;
herr_t status;
A.save(group, status);
```

#### Parent/Child Structure Sharing

Two matrices with the same sparsity pattern can share their structure arrays
(pointers, indices, and, when present, the COO index array) to save memory.
The second matrix is constructed as a *child* of the first:

```cpp
// parent owns the structure
SpMatrix* tK = new SpMatrix(graph, SpMatrixType::CSR, n, n, false);

// child aliases the structure and owns only its value array (zero-filled)
SpMatrix* tM = new SpMatrix(tK);
```

This is how `SolverData` pairs `mSystemMatrix`/`mJacobianMatrix` and
`mFullMassMatrix`/`mFullStiffnessMatrix`. Memory saved per pair is roughly
`(nnz + n + 1) * sizeof(int_t)`, plus the COO array when MUMPS is used.

**Ownership and lifecycle rules:**

- While linked, the parent owns the pointer, index, and optional COO index
  arrays; the child owns only its values. `memory()` reports that split.
- One child per parent; chains (a child of a child) are refused.
- Either destruction order is safe: destroying the child unlinks it from the
  parent; destroying the parent first transfers structure ownership to the child.
- `set_indexing_base()` on either side converts the shared arrays once and
  updates both index functions. `create_coo_indices()`/`free_coo_indices()`
  keep both sides in sync in either call order.
- `load()` on a linked matrix (child, or a parent with an attached child, as in
  `SolverData::load_system`) never replaces the structure: it verifies the
  file's pattern against the existing arrays (tolerating an indexing-base
  offset) and loads the values only. A former child that inherited the
  structure after its parent's destruction is an ordinary unlinked matrix
  again, and `load()` rebuilds it from the file as usual. `save()` needs no
  special handling.

**Refused on linked matrices** (`BELFEM_ERROR`): `transpose()`,
`sort_entries()`, `set_type()`, copy-assignment onto a linked matrix, and
move-assignment from or onto a linked matrix. Assigning *from* a linked matrix
into an unlinked matrix is allowed and produces an ordinary deep copy.

#### When to Use

- Primary sparse matrix type for all BELFEM FEM assembly
- Large sparse linear systems (> 1000 DOFs)
- Graph-based connectivity (mesh elements, finite differences)
- Matrix-free iterative methods (via multiply)

**See:** `SpMatrix` class in `cl_SpMatrix.hpp`

---

### Solver - Unified Solver Interface

**Files:** `cl_Solver.hpp`

#### Description

Unified interface to multiple sparse direct and iterative solver libraries. Automatically selects the solver backend compiled into BELFEM.

#### Construction

```cpp
// Use default solver (determined at compile time)
Solver solver;

// Explicitly specify solver type
Solver solver(SolverType::MUMPS);
Solver solver(SolverType::STRUMPACK);
Solver solver(SolverType::PARDISO);
Solver solver(SolverType::UMFPACK);
Solver solver(SolverType::PETSc);

// From solver parameters
SolverParameters params(SolverType::MUMPS);
params.set_reordering_method(ReorderingMethod::METIS);
Solver solver(params);
```

#### Solving Linear Systems

```cpp
SpMatrix A(graph, SpMatrixType::CSR);
Vector<real> x(n), b(n);
// ... assemble A and b ...

Solver solver(SolverType::MUMPS);

// Symmetry mode. MUMPS accepts Unsymmetric ONLY -- see the note below.
solver.set_symmetry_mode(SymmetryMode::Unsymmetric);   // the default

// Solve: A * x = b
solver.solve(A, x, b);  // x is overwritten with solution

// Multiple RHS (solve for multiple vectors)
Matrix<real> X(n, m), B(n, m);
solver.solve(A, X, B);  // Each column is a separate RHS
```

#### Solver-Specific Configuration

**PETSc (Iterative Solvers):**
```cpp
Solver solver(SolverType::PETSc);

// Set preconditioner and Krylov method
solver.set_petsc(
    Preconditioner::GAMG,      // Algebraic multigrid
    KrylovMethod::GMRES,       // Generalized Minimal Residual
    1e-8                       // Relative tolerance
);

// Other preconditioners:
// - Preconditioner::JACOBI (diagonal scaling)
// - Preconditioner::ASM (additive Schwarz)
// - Preconditioner::ILU (incomplete LU)
// - Preconditioner::ICC (incomplete Cholesky)

// Other Krylov methods:
// - KrylovMethod::CG (conjugate gradient, SPD only)
// - KrylovMethod::GMRES (general)
// - KrylovMethod::BCGS (biconjugate gradient stabilized)
```

**MUMPS (Parallel Direct Solver):**
```cpp
Solver solver(SolverType::MUMPS);

// Set reordering for fill reduction
solver.set_mumps_reordering(
    MumpsSerialReodrdering::METIS,
    MumpsParallelReodrdering::PARMETIS
);

// Enable block low-rank compression.
// CAVEAT: only reliable AFTER the wrapper's initialize()
// — a call before it gets its epsilon clobbered from SolverParameters.
// The supported path is the deck: compression scheme : blr with a
// compression cutoff 2-3 decades below the nonlinear tolerance; see
// solver_memory_and_compression.md §4. 1e-4 shown here is the value
// that silently stalled Newton on 2026-07-06 — do not copy it onto an
// HTS deck.
solver.set_mumps_blr(
    MumpsBlockLowRanking::FactorizationAndSolution,
    1e-4  // Compression tolerance (see caveat above)
);

// Error analysis (for debugging)
solver.set_mumps_error_analysis(MumpsErrorAnalysis::Full);
```

#### Command-Line Pass-Through (PETSc and STRUMPACK)

Both libraries read their own options directly from the command line of any
BELFEM executable — no code is required in the executable:

- **PETSc**: `Communicator::init` hands the full `argc`/`argv` to
  `PetscInitialize` (`cl_Communicator.cpp`, search `PetscInitialize`), which
  loads every argument into the PETSc options database. All PETSc options work
  as documented by PETSc (`-ksp_type gmres`, `-pc_type hypre`, `-ksp_monitor`,
  `-options_file petsc.opts`, `-snes_*`, …). Unknown arguments are ignored;
  the warning about unused options is suppressed via `-options_left 0`.
- **STRUMPACK**: the wrapper copies the arguments from the global communicator
  and calls `options().set_from_command_line()` **after** applying the
  settings derived from `input.conf` (`cl_SolverSTRUMPACK.cpp`, search
  `set_from_command_line`; both the serial and the distributed path do this).
  The command line therefore overrides the input file. STRUMPACK flags carry
  the `--sp_` prefix (`--sp_compression NONE`, `--sp_reordering_method metis`,
  …); run any STRUMPACK-enabled executable with `--help` **not** intercepted
  by BELFEM, i.e. consult the STRUMPACK manual for the full list.

STRUMPACK's own progress output is tied to the logger: it is enabled when the
info level is at least 5, e.g. by passing `--verbose` (see the core module's
`Arguments` class).

Example:

```bash
mpirun -np 4 belfem -v 3 --sp_compression NONE -ksp_monitor
```

Each parser skips flags it does not recognize, so BELFEM, PETSc, and
STRUMPACK options can be mixed freely on one command line.

#### Cleanup

```cpp
// Manually free solver memory (also done by destructor)
solver.free();
```

#### When to Use

- All sparse linear system solves in BELFEM
- FEM stiffness matrix solves
- Eigenvalue problems (via iterative methods)
- Nonlinear problems (repeated solves with updated matrices)

**See:** `Solver` class in `cl_Solver.hpp`

---

### SolverParameters - Solver Configuration

**Files:** `cl_SolverParameters.hpp`

#### Description

Configuration object for solver settings, supporting MPI broadcasting and input file parsing.

#### Construction

```cpp
// Default parameters for solver type
SolverParameters params(SolverType::MUMPS);

// From input file section
input::Section* input = ...;
SolverParameters params(input);

// Copy constructor
SolverParameters params2(params);
```

#### Configuration

```cpp
SolverParameters params(SolverType::PETSc);

// Matrix format (for distributed solvers)
params.set_distributed_matrix_type(DistributedMatrixType::AIJ);
// or: DistributedMatrixType::CSR
// or: DistributedMatrixType::CSC

// Reordering for fill reduction
params.set_reordering_method(ReorderingMethod::METIS);
// or: ReorderingMethod::SCOTCH
// or: ReorderingMethod::NATURAL (no reordering)
// or: ReorderingMethod::PARMETIS / PTSCOTCH — explicit parallel library;
//     identical to METIS / SCOTCH on MUMPS and STRUMPACK at > 1 rank, and the
//     BELFEM-side parallel nested dissection on the PETSc path

// Compression (STRUMPACK, MUMPS)
params.set_compression_method(CompressionMethod::BLR);  // Block low-rank
// or: CompressionMethod::OFF

// Preconditioner (PETSc)
params.set_preconditioner(Preconditioner::GAMG);

// Krylov method (PETSc, STRUMPACK)
params.set_krylov_method(KrylovMethod::GMRES);

// Relative tolerance
params.set_relative_tolerance(1e-8);

// Use initial guess (iterative solvers)
params.set_use_initial_guess(true);
```

#### Querying

```cpp
SolverType type = params.type();
DistributedMatrixType mat_type = params.distributed_matrix_type();
Preconditioner prec = params.preconditioner();
KrylovMethod krylov = params.krylov_method();
ReorderingMethod reorder = params.reordering_method();
CompressionMethod compress = params.compression_method();
real tol = params.relative_tolerance();
bool use_guess = params.use_initial_guess();
```

#### MPI Synchronization

```cpp
// On rank 0: set parameters
SolverParameters params(SolverType::MUMPS);
if (gComm.rank() == 0) {
    params.set_reordering_method(ReorderingMethod::METIS);
}

// Broadcast to all ranks
params.synchronize();

// All ranks now have same parameters
```

**See:** `SolverParameters` class in `cl_SolverParameters.hpp`

---

## Solver Backend Comparison

| Backend | Type | MPI | Shared-Mem | Compression | When to Use |
|---------|------|-----|------------|-------------|-------------|
| **UMFPACK** | Direct | No | No | No | Sequential, moderate size (< 100k DOFs) |
| **MUMPS** | Direct | Yes | Yes | BLR | Parallel, general matrices, robust |
| **STRUMPACK** | Direct/Iterative | Yes | Yes | BLR, HSS | Parallel, best performance, modern |
| **PARDISO** | Direct | No | Yes (OpenMP) | No | Shared-memory parallel, Intel CPUs |
| **PETSc** | Iterative | Yes | No | N/A | Very large systems, SPD matrices |

### Performance Characteristics

**UMFPACK (SuiteSparse):**
- Sequential only
- Excellent for unsymmetric matrices
- Mature and stable
- Best for < 100k DOFs on single node

**MUMPS:**
- MPI + OpenMP hybrid parallelism
- Robust for ill-conditioned systems
- Good for general unsymmetric matrices
- Block low-rank compression available
- Reliable fallback when STRUMPACK fails

**STRUMPACK (Recommended):**
- MPI + OpenMP hybrid parallelism
- Best overall performance in BELFEM
- Low-rank compression (BLR, HSS)
- Both direct and iterative solvers
- Can be 2× faster than MUMPS (Messe et al. 2023)
- Preferred for production runs

**PARDISO (Intel MKL):**
- Shared-memory parallelism via OpenMP
- Excellent single-node performance
- Optimized for Intel CPUs
- Limited to shared-memory systems

**PETSc:**
- Scalable to very large systems (millions of DOFs)
- Flexible Krylov method selection
- Requires good preconditioner for fast convergence
- Best for SPD matrices with CG
- GAMG preconditioner works well for elliptic problems

### Default Solver Selection

The default solver is determined at compile time (see `en_SolverEnums.hpp:193-203`):

```cpp
#ifdef BELFEM_STRUMPACK
    gDefaultSolver = SolverType::STRUMPACK;  // Preferred
#elif BELFEM_MUMPS
    gDefaultSolver = SolverType::MUMPS;      // Robust fallback
#elif BELFEM_PARDISO
    gDefaultSolver = SolverType::PARDISO;    // Shared-memory
#elif BELFEM_SUITESPARSE
    gDefaultSolver = SolverType::UMFPACK;    // Sequential
#elif BELFEM_SUPERLU
    gDefaultSolver = SolverType::SUPERLU;    // Sequential fallback
#endif
```

**Recommendation:** Compile with STRUMPACK for production, MUMPS for robustness, PETSc for extreme scale.

---

## Free Functions

### `preferred_matrix_format(solver_type)` - Optimal Matrix Format

Returns the preferred sparse matrix format for a given solver.

```cpp
SolverType solver = SolverType::MUMPS;
SpMatrixType format = preferred_matrix_format(solver);
// Returns SpMatrixType::CSR for MUMPS

// Use when constructing matrix
SpMatrix A(graph, preferred_matrix_format(solver));
```

**Mapping:**
- UMFPACK → CSC
- MUMPS → CSR
- STRUMPACK → CSR
- PARDISO → CSR
- PETSc → CSR (converted to AIJ internally)

**File:** `fn_preferred_matrix_format.hpp`

---

### `matrix_type(solver_type)` - Storage Format for a Solver

Returns `SpMatrixType::CSC` for UMFPACK, SuperLU and MUMPS, `CSR` for PARDISO,
PETSc and STRUMPACK. It never sees a matrix and does not detect symmetry. Note that
`preferred_matrix_format()` answers the same question with CSR for MUMPS; use one
consistently.

```cpp
SpMatrixType tType = matrix_type(SolverType::MUMPS);  // CSC
```

**File:** `fn_matrix_type.hpp`

---

### `compute_permutation(A, graph, forward, backward, index)` - METIS Permutation

Builds a graph from `A`, runs METIS nested dissection on it and fills three
`Cell<int_t>` permutations (row forward/backward, and new→old nonzero positions).
No method argument: METIS only. `graph` must be empty on entry.

```cpp
SpMatrix A(graph_in, SpMatrixType::CSR);
Graph graph;
Cell<int_t> forward, backward, index;

sparse::compute_permutation(A, graph, forward, backward, index);

// Apply permutation (manually or within solver)
// Most solvers do this internally
```

**File:** `fn_compute_permutation.hpp`

---

### `create_graph_from_matrix(A)` - Extract Connectivity Graph

Extracts connectivity graph from sparse matrix structure.

```cpp
SpMatrix A(graph_in, SpMatrixType::CSR);
Graph graph_out;
sparse::create_graph_from_matrix(A, graph_out);  // A must be 0-based; graph_out must be empty

// Use for mesh partitioning, reordering, etc.
```

**File:** `fn_create_graph_from_matrix.hpp`

---

### `rcond(A)` - Reciprocal Condition Number

Estimates the reciprocal condition number via CHOLMOD (requires `USE_SUITESPARSE=ON`,
default OFF, and an OpenMP-enabled build; aborts otherwise).

```cpp
SpMatrix A(graph, SpMatrixType::CSR);
// ... assemble A ...

real rc = rcond(A);
// rc ≈ 1.0: well-conditioned
// rc ≈ 1e-16: ill-conditioned (near-singular)

if (rc < 1e-12) {
    message(InfoLevel::Minimal, "Matrix is ill-conditioned");
}
```

**Note:** Expensive operation (requires factorization). Use for diagnostics only.

**File:** `fn_rcond.hpp`

---

## Common Patterns

### Pattern 1: Solving a Sparse Linear System

```cpp
// Build connectivity graph from FEM mesh
Graph graph = build_dof_graph(mesh);

// Create sparse matrix
SpMatrix K(graph, SpMatrixType::CSR);

// Assemble element contributions
for (Element* e : mesh->elements()) {
    Matrix<real> Ke = e->stiffness();
    Cell<index_t> dofs = e->dof_indices();

    for (index_t i = 0; i < dofs.size(); ++i) {
        for (index_t j = 0; j < dofs.size(); ++j) {
            K(dofs(i), dofs(j)) += Ke(i, j);
        }
    }
}

// Create RHS vector
Vector<real> f(n_dofs, 0.0);
// ... assemble force vector ...

// Solve K * u = f
Vector<real> u(n_dofs);
Solver solver;
solver.solve(K, u, f);

// u now contains solution
```

---

### Pattern 2: Symmetric Positive Definite System

```cpp
SpMatrix K(graph, SpMatrixType::CSR);
// ... assemble symmetric positive definite matrix ...

// A symmetric mode is NOT available with MUMPS: under SYM != 0 it wants ONE
// representative per symmetric coordinate ( either triangle ), while BELFEM
// supplies the full matrix -- see Performance Tip 4, "Exploit Symmetry".
// The matrix being SPD does not change that -- SYM = 1 has the same contract
Solver solver;
solver.set_symmetry_mode(SymmetryMode::Unsymmetric);

// Solve
Vector<real> u(n), f(n);
solver.solve(K, u, f);

// Faster than general solve
```

---

### Pattern 3: Multiple Solves with Same Matrix

```cpp
SpMatrix K(graph, SpMatrixType::CSR);
// ... assemble K (remains constant) ...

Solver solver;

// First solve
Vector<real> u1(n), f1(n);
// ... set f1 ...
solver.solve(K, u1, f1);

// Second solve (matrix already factorized)
Vector<real> u2(n), f2(n);
// ... set f2 ...
solver.solve(K, u2, f2);  // REFACTORIZES: a plain second solve rebuilds the factors on every backend

// To reuse them, arm the frozen scope on a wrapper that supports it (MUMPS):
// if ( solver.wrapper()->supports_factorization_reuse() )
// {
//     solver.wrapper()->freeze_factorization( K );   // after a successful solve
//     solver.solve(K, u2, f2);                       // solve only, against the saved factors
//     solver.wrapper()->unfreeze_factorization();
// }
```

---

### Pattern 4: Iterative Solve with PETSc

```cpp
SpMatrix A(graph, SpMatrixType::CSR);
Vector<real> x(n), b(n);
// ... assemble A and b ...

// Configure PETSc
Solver solver(SolverType::PETSc);
solver.set_petsc(
    Preconditioner::GAMG,  // Algebraic multigrid for elliptic problems
    KrylovMethod::CG,      // Conjugate gradient (for SPD)
    1e-10                  // Tight tolerance
);

// Use initial guess (optional)
x.fill(0.0);  // Zero initial guess

solver.solve(A, x, b);
```

---

### Pattern 5: MPI Parallel Solve

```cpp
// Master rank: build and assemble the matrix; workers hold an empty shell
SpMatrix A;
if (gComm.rank() == 0) {
    Graph graph = build_dof_graph(mesh);
    A = SpMatrix(graph, SpMatrixType::CSR);
    // ... assemble A ...
}

// solve() is collective on all ranks; the wrapper distributes internally
// (MUMPS: host-centralized input; STRUMPACK/PETSc: sparse::DistMatrix
//  scatters pointer/index/value slices from rank 0)
Vector<real> x(n), b(n);
Solver solver(SolverType::MUMPS);  // or STRUMPACK
solver.solve(A, x, b);
```

---

### Pattern 6: Saving and Reusing Matrix

```cpp
// Assemble and save
SpMatrix K(graph, SpMatrixType::CSR);
// ... assemble K ...
K.save("stiffness.h5", "K", FileMode::NEW);

// Later: load and solve
SpMatrix K_loaded("stiffness.h5", "K");
Solver solver;
Vector<real> u(n), f(n);
solver.solve(K_loaded, u, f);
```

---

## Performance Tips

### 1. **Choose the Right Solver**

```cpp
// Sequential (< 100k DOFs)
Solver solver(SolverType::UMFPACK);

// Parallel, general (100k - 10M DOFs)
Solver solver(SolverType::STRUMPACK);  // Fastest
Solver solver(SolverType::MUMPS);      // Most robust

// Parallel, very large (> 10M DOFs, SPD)
Solver solver(SolverType::PETSc);
solver.set_petsc(Preconditioner::GAMG, KrylovMethod::CG, 1e-8);
```

### 2. **Enable Compression for Large Systems**

```cpp
// MUMPS with block low-rank compression
Solver solver(SolverType::MUMPS);
solver.set_mumps_blr(
    MumpsBlockLowRanking::FactorizationAndSolution,
    1e-4  // Compression tolerance (larger = more compression, less accuracy)
);

// Can reduce memory by 50-90% for large 3D problems
```

### 3. **Use Optimal Reordering**

```cpp
// METIS generally best for 3D FEM
Solver solver(SolverType::MUMPS);
solver.set_mumps_reordering(
    MumpsSerialReodrdering::METIS,
    MumpsParallelReodrdering::PARMETIS
);

// Can reduce fill-in by 10-100×
```

### 4. **Exploit Symmetry**

> **Not available with MUMPS.** A symmetric mode would cut factor time and memory substantially,
> but under `SYM != 0` MUMPS wants exactly ONE representative of each symmetric coordinate. Either
> triangle is acceptable — it does not insist on the lower one — and what is fatal is supplying
> BOTH: `(i,j)` and `(j,i)` are then summed as duplicates. BELFEM stores and hands over the FULL
> matrix and nothing extracts a triangle, so every off-diagonal would arrive twice and MUMPS would
> factorize a different matrix -- without failing. Measured on a matrix with a known spectrum:
> a first-solve residual of 7.4e16 and a converged eigenvalue of -2.5e-19 against a true 2.46e-6.
>
> `MUMPS::initialize()` therefore rejects any symmetric mode with an always-active error rather
> than accepting one it cannot honor. To enable it, supply one triangle from the wrapper — see
> `todo/mumps_symmetric_triangle_extraction.md`.
>
> SuperLU is unaffected: its `SymmetricMode` is a pivoting and ordering heuristic, not a
> triangle-storage contract, so it reads the full matrix as intended.

```cpp
// MUMPS: unsymmetric only
solver.set_symmetry_mode(SymmetryMode::Unsymmetric);
```

### 5. **Use Optimal Matrix Format**

```cpp
// Use preferred format for your solver
SolverType solver_type = SolverType::MUMPS;
SpMatrix A(graph, preferred_matrix_format(solver_type));  // Returns CSR

// Or explicit mapping:
// UMFPACK → CSC
// MUMPS, STRUMPACK, PARDISO, PETSc → CSR
SpMatrixType format = (solver_type == SolverType::UMFPACK)
                       ? SpMatrixType::CSC : SpMatrixType::CSR;
SpMatrix A(graph, format);
```

### 6. **Minimize Matrix Copies**

```cpp
// SpMatrix is not copy- or move-constructible: build it in place
SpMatrix A(graph, SpMatrixType::CSR);
// ... assemble A ...
solver.solve(A, x, b);
// A still available for analysis

// a deep copy goes through the assignment operator on an empty matrix
SpMatrix A_backup;
A_backup = A;   // explicit, intentional copy
```

---

## MPI Considerations

### Matrix Distribution

BELFEM uses host-centralized sparse-matrix input for MPI-capable solvers:
rank 0 owns the full `SpMatrix`, and the solver wrapper or backend distributes
the work internally. MUMPS consumes the rank-0 matrix through the library
interface; STRUMPACK and PETSc build `sparse::DistMatrix` slices from the
rank-0 pointer, index, and value arrays.

```cpp
// Host-centralized input for MPI-capable solvers: rank 0 owns the full matrix,
// the solver wrapper distributes internally
SpMatrix A;
if (gComm.rank() == 0) {
    A = SpMatrix(graph, SpMatrixType::CSR);
}
solver.solve(A, x, b);   // collective call
```

### Solver Parallelism

| Solver | MPI | OpenMP | Best Use Case |
|--------|-----|--------|---------------|
| UMFPACK | No | No | Sequential |
| MUMPS | Yes | Yes | Hybrid parallel |
| STRUMPACK | Yes | Yes | Hybrid parallel (best) |
| PARDISO | No | Yes | Shared-memory |
| PETSc | Yes | No | Pure MPI |

**Hybrid Parallel Example (STRUMPACK/MUMPS):**
```bash
# 4 MPI ranks, 8 OpenMP threads each = 32 cores total
export OMP_NUM_THREADS=8
mpirun -np 4 ./my_belfem_app
```

### BELFEM's Own OpenMP Kernels: `USE_BELFEM_OPENMP`

The OpenMP threads listed above belong to the third-party solvers. BELFEM's own sparse
kernels use a separate switch. That switch is **OFF** by default:

| CMake option | Default | Governs |
|---|---|---|
| `USE_OPENMP` | ON | `-fopenmp` for the whole build; the solvers need it, and the thread-budget queries (`hatch_turtle()` in `cl_SolverWrapper.cpp`, the banner's `Threads Used` line, PARDISO's `gParameters(3)`) key off the plain `OMP` define |
| `USE_BELFEM_OPENMP` | **OFF** | only the `!$omp` directives in `splinalg.f90` (`matvec_csr`, `matvec_csc`), `arpacktools.f90` and `parpacktools.f90`, through the `BELFEM_OMP` define (`CMakeLists.txt`, appended under `USE_OPENMP AND USE_BELFEM_OPENMP`) |

These switches are not interchangeable. Moving the thread-budget queries to `BELFEM_OMP`
would remove the oversubscription warning and silently serialize PARDISO on every default
build. The three sites have comments that state this.

**Why the kernels are OFF.** `matvec_csc` parallelizes its column scatter with an OpenMP
array reduction, `!$omp reduction(+:y)` (`splinalg.f90`, `matvec_csc`). `y` is an
explicit-shape dummy, so gfortran creates a private copy of the entire result vector in
each worker thread's stack frame. A libgomp worker gets 2 MiB on Darwin. As a result, a
threaded CSC product dies with SIGBUS once

    n > 2 MiB / 8 bytes = 262,144 rows

Ordinary 3D meshes can exceed that size easily. The threading also provides no measurable
benefit: both call paths for these kernels run master-only (the Newton residual in
`cl_FEM_DofMgr_SolverData.cpp` and the ARPACK inverse iteration in
`cl_FEM_DofMgr_EigenValues.cpp`). The kernels cost about two flops per nonzero and run
alongside a full assembly and a direct factorization at each step. No threaded-versus-serial
timing for them exists anywhere in the tree. The parallelism was therefore switched off
instead of repaired, as stated in the kernel comment above the reduction. One useful side
effect is that, with the switch OFF, `matvec_csc` is bitwise deterministic across thread
counts.

**Under an MKL build the switch is moot.** The two-argument `SpMatrix::multiply` forwards to
the five-argument overload (`cl_SpMatrix.cpp`, `SpMatrix::multiply`). Its `BELFEM_MKL` branch
calls `mkl_sparse_d_mv`; both remaining `matvec_csc`/`matvec_csr` call sites are inside its
`#else`. An MKL build therefore never reaches the Fortran kernels on these paths.

**When to turn it ON.** Use it only in the following two situations, never in a production
or CI build:

1. **To measure whether the threading ever pays.** Build once with
   `-DUSE_BELFEM_OPENMP=ON` and once without it. Run the same deck with the same
   `OMP_NUM_THREADS`, then compare per-phase timings instead of step wall-time. Raise
   `OMP_STACKSIZE` for the ON run (`OMP_STACKSIZE=64M` carries a 295,315-dof case that
   crashes at the default). If the measurement shows a gain, the private copy must still
   move off the worker stack before ON can become a supported configuration.
2. **To prove the switch is a two-way toggle:** make one build with the token ON and run
   `make check`.

**What ON does not do.** It does not change the third-party solvers' threading. STRUMPACK,
MKL and a threaded BLAS spawn their own workers regardless, and a libgomp worker gets the
same 2 MiB default stack on Darwin. Therefore, `OMP_STACKSIZE` remains a useful diagnostic
knob with the switch OFF, and a fault inside a `gomp_thread_start` frame is not
automatically a BELFEM bug. The switch also
has no effect on an MKL build, for the reason given above.

**`OMP_NUM_THREADS=1` hides the failure rather than fixing it** — the work then runs on the
8 MiB main-thread stack. Use this as a diagnostic, not as a configuration.

See `doc/parallel_execution.md` for the wider picture. The matrix backends (Armadillo,
Blaze) parallelize their own expression evaluation whenever OpenMP is available. Neither
uses an OpenMP reduction, so neither can hit this failure.

---

## Thread Safety

**Thread Safety:** Like other BELFEM modules, sparse matrices and solvers are **not thread-safe**.

- **Reading:** Safe from multiple threads (const operations)
- **Writing:** Requires external synchronization
- **Solving:** One solve per solver instance at a time

**Example with OpenMP:**

```cpp
SpMatrix K(graph, SpMatrixType::CSR);

// Parallel assembly (each thread writes different elements)
#pragma omp parallel for
for (int e = 0; e < n_elements; ++e) {
    Matrix<real> Ke = compute_element_stiffness(e);
    Cell<index_t> dofs = get_element_dofs(e);

    // Critical section: writing to sparse matrix
    #pragma omp critical
    {
        for (index_t i = 0; i < dofs.size(); ++i) {
            for (index_t j = 0; j < dofs.size(); ++j) {
                K(dofs(i), dofs(j)) += Ke(i, j);
            }
        }
    }
}

// Serial solve (no threading needed)
Solver solver;
solver.solve(K, u, f);
```

---

## Enumerations Reference

### SolverType

```cpp
enum class SolverType {
    UMFPACK,     // SuiteSparse sequential
    SUPERLU,     // Sequential direct, default-ON (USE_SUPERLU)
    MUMPS,       // Parallel direct
    STRUMPACK,   // Parallel direct/iterative (preferred)
    PARDISO,     // Shared-memory parallel
    PETSc,       // Iterative framework
    UNDEFINED
};
```

### SpMatrixType

```cpp
enum class SpMatrixType {
    CSC,         // Compressed Sparse Column
    CSR,         // Compressed Sparse Row
    UNDEFINED
};
```

### SymmetryMode

```cpp
enum class SymmetryMode {
    Unsymmetric = 0,                // General matrix
    PositiveDefiniteSymmetric = 1,  // SPD (Cholesky)
    GeneralSymmetric = 2,           // Symmetric (LDLT)
    UNDEFINED
};
```

### Preconditioner (PETSc only)

```cpp
enum class Preconditioner {
    NONE,       // No preconditioning
    JACOBI,     // Diagonal scaling
    BJACOBI,    // Block Jacobi
    ASM,        // Additive Schwarz
    GAMG,       // Geometric algebraic multigrid (best for elliptic)
    ILU,        // Incomplete LU
    ICC,        // Incomplete Cholesky (SPD only)
    LU,         // Direct solve as preconditioner
    UNDEFINED
};
```

### KrylovMethod (PETSc, STRUMPACK)

```cpp
enum class KrylovMethod {
    PREONLY,    // Direct solve only (no Krylov)
    CG,         // Conjugate Gradient (SPD only)
    GMRES,      // Generalized Minimal Residual (general)
    BCGS,       // BiConjugate Gradient Stabilized
    CGS,        // Conjugate Gradient Squared
    TFQMR,      // Transpose-Free QMR
    AUTO,       // Let solver choose
    UNDEFINED
};
```

### ReorderingMethod

```cpp
enum class ReorderingMethod {
    NATURAL = 0,  // No reordering (original order)
    METIS = 1,    // METIS graph partitioner (recommended)
    SCOTCH = 2,   // SCOTCH graph partitioner
    AUTOMATIC = 3, // the default: library's own choice
    PARMETIS = 4, // explicit parallel ND (PETSc: BELFEM-side; MUMPS/STRUMPACK: as METIS)
    PTSCOTCH = 5, // explicit parallel ND (PETSc: BELFEM-side; MUMPS/STRUMPACK: as SCOTCH)
    UNDEFINED = 6 // appended, never inserted: the value is synchronized as a uint
};
```

### CompressionMethod

```cpp
enum class CompressionMethod {
    OFF = 0,       // No compression
    BLR = 1,       // Block Low-Rank compression (explicit opt-in only)
    AUTOMATIC = 2, // the deck default — identical to OFF on every
                   // library since 2026-08-16; see
                   // solver_memory_and_compression.md
    UNDEFINED = 3
};
```

---

## Related Modules

- **linalg**: Dense `Matrix<T>` and `Vector<T>` for element-level operations
- **fem/kernel**: Uses sparse module for global system assembly and solve
- **mesh**: Provides connectivity graphs for sparse matrix construction
- **comm**: MPI utilities for distributed sparse matrices

---

## Debugging Tips

### 1. **Check Matrix Symmetry**

```cpp
SpMatrix A(graph, SpMatrixType::CSR);
// ... assemble ...

// Verify structural symmetry
for (index_t i = 0; i < A.n_rows(); ++i) {
    for (index_t j = i+1; j < A.n_cols(); ++j) {
        real Aij = A(i, j);
        real Aji = A(j, i);
        if (std::abs(Aij - Aji) > 1e-12) {
            message(InfoLevel::Minimal,
                    "Asymmetry at (%lu,%lu): %e vs %e", i, j, Aij, Aji);
        }
    }
}
```

### 2. **Visualize Sparsity Pattern**

```cpp
SpMatrix A(graph, SpMatrixType::CSR);
A.print("Stiffness");  // Print matrix in readable form

// For very large matrices:
message(InfoLevel::Default,
        "Matrix: %lu × %lu, nnz = %lu, fill = %.2f%%",
        A.n_rows(), A.n_cols(), A.number_of_nonzeros(),
        100.0 * A.number_of_nonzeros() / (A.n_rows() * A.n_cols()));
```

### 3. **Check Conditioning**

```cpp
real rc = rcond(A);
if (rc < 1e-12) {
    message(InfoLevel::Minimal,
            "Matrix is ill-conditioned (rcond = %e)", rc);
}
```

### 4. **Enable Solver Diagnostics**

```cpp
// MUMPS error analysis
Solver solver(SolverType::MUMPS);
solver.set_mumps_error_analysis(MumpsErrorAnalysis::Full);

// Check residual after solve
Vector<real> residual = A * x - b;
real res_norm = norm(residual);
message(InfoLevel::Default, "Residual norm: %e", res_norm);
```

---

## See Also

- **Sparse Solver Libraries:**
  - [UMFPACK Documentation](https://people.engr.tamu.edu/davis/suitesparse.html)
  - [MUMPS User Guide](http://mumps.enseeiht.fr/)
  - [STRUMPACK Documentation](https://portal.nersc.gov/project/sparse/strumpack/)
  - [PETSc Manual](https://petsc.org/release/manual/)
  - [Intel MKL PARDISO](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2024-0/pardiso.html)

- **BELFEM Modules:**
  - [Linear Algebra Module](../../linalg/doc/README.md) - Dense vector/matrix operations
  - [FEM Module](../../fem/kernel/doc/README.md) - Finite element assembly
  - [Mesh Module](../../mesh/doc/README.md) - Connectivity graphs

- **Build System:**
  - `CLAUDE.md` (repository root) - Build configuration and solver selection
  - CMakeLists.txt - Compile-time solver flags

---

**Revision History:**
- 2026-01-16: Initial comprehensive documentation

**Contributors:**
- Based on BELFEM sparse module by Christian Messe and Gregory Giard
- Documentation synthesized from code analysis and usage patterns

**Prepared by the BELFEM documentation team.**
