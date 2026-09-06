# Sparse Matrix Module Documentation {#sparse_index}

Documentation for BELFEM's sparse matrix storage and solver interfaces.

## Contents

### Guides

- [**sparse_usage_guide.md**](sparse_usage_guide.md) - Comprehensive usage guide for BELFEM sparse matrices and solvers
  - Overview of solver backends (UMFPACK, MUMPS, STRUMPACK, PARDISO, PETSc)
  - SpMatrix, Solver, and SolverParameters API reference
  - Matrix formats (CSR, CSC, COO) and indexing
  - Solver-specific configuration (reordering, compression, preconditioners)
  - BELFEM's own OpenMP kernels (`USE_BELFEM_OPENMP`): why OFF, when ON
  - MPI distribution and parallel solving
  - Performance tips and common patterns
- [solver_memory_and_compression.md](solver_memory_and_compression.md) — choosing STRUMPACK vs MUMPS, per-rank memory model, when BLR compression is safe (unified contract, 2026-08)

## Quick Reference

### Core Types

| Type | Description | File |
|------|-------------|------|
| `SpMatrix` | Sparse matrix (CSR/CSC format) | `cl_SpMatrix.hpp` |
| `Solver` | Unified sparse solver interface | `cl_Solver.hpp` |
| `SolverParameters` | Solver configuration | `cl_SolverParameters.hpp` |

Two matrices with identical sparsity patterns can share their structure arrays
using the explicit child constructor `SpMatrix( SpMatrix * aParent )`.
See "Parent/Child Structure Sharing" in the
[usage guide](sparse_usage_guide.md).

### Sparse Matrix Formats

| Format | Description | Preferred By |
|--------|-------------|--------------|
| **CSR** | Compressed Sparse Row | MUMPS, STRUMPACK, PARDISO, PETSc |
| **CSC** | Compressed Sparse Column | UMFPACK |
| **COO** | Coordinate (row/col indices) | MUMPS (required) |

### Solver Backends

| Solver | Type | Parallelism | Best For |
|--------|------|-------------|----------|
| **UMFPACK** | Direct | Sequential | < 100k DOFs, unsymmetric |
| **MUMPS** | Direct | MPI + OpenMP | Robust, general matrices |
| **STRUMPACK** | Direct/Iterative | MPI + OpenMP | Best performance (recommended) |
| **PARDISO** | Direct | OpenMP | Shared-memory, Intel CPUs |
| **PETSc** | Iterative | MPI | Very large, SPD matrices |
| **SuperLU** | Direct | Sequential | Default-ON backend (`USE_SUPERLU`); the wrapper registers itself as non-MPI (`cl_SolverSUPERLU.cpp:28-30`) |

### Backend Selection

```cmake
# Compile-time selection via CMake flags
cmake -DUSE_STRUMPACK=ON ..    # Recommended (default ON)
cmake -DUSE_MUMPS=ON ..        # Robust fallback (default ON)
cmake -DUSE_PARDISO=ON ..      # Shared-memory (opt-in)
cmake -DUSE_SUITESPARSE=ON ..  # Sequential UMFPACK (opt-in)
cmake -DUSE_PETSC=ON ..        # Iterative (default ON)
```

`USE_OPENMP` (default ON) gives the third-party solvers their threads. BELFEM's **own**
Fortran matvec kernels use threads only when `USE_BELFEM_OPENMP` is enabled. That option is
**OFF** and should remain OFF outside a measurement session: the `matvec_csc` array
reduction creates an n × 8 byte private copy on each 2 MiB worker stack and crashes above
262,144 rows. See the usage guide, "BELFEM's Own OpenMP Kernels".

### Basic Usage Pattern

```cpp
// Build sparse matrix from graph
Graph graph = build_connectivity_graph();
SpMatrix K(graph, SpMatrixType::CSR);

// Assemble. There is no free `assemble()` helper in src/sparse -- element
// matrices go into the sparse matrix through its own accessors.
for (Element* e : elements) {
    Matrix<real> Ke = e->stiffness();
    // ... add Ke into K at the rows/cols given by e->dofs()
}

// Solve K * u = f
Vector<real> u(n), f(n);
Solver solver;
solver.solve(K, u, f);
```

## Free Functions

| Function | Description | File |
|----------|-------------|------|
| `preferred_matrix_format()` | Optimal format for solver | `fn_preferred_matrix_format.hpp` |
| `matrix_type()` | Storage format for a solver type (CSC/CSR) | `fn_matrix_type.hpp` |
| `compute_permutation()` | Fill-reducing reordering | `fn_compute_permutation.hpp` |
| `create_graph_from_matrix()` | Extract connectivity | `fn_create_graph_from_matrix.hpp` |
| `rcond()` | Condition number estimate | `fn_rcond.hpp` |

## See Also

- **Source code** - header files in `src/sparse/`
- [Linear algebra module](../../linalg/doc/README.md) - Dense matrices and vectors
- [FEM module](../../fem/kernel/doc/README.md) - Uses sparse module for assembly/solve
- [Root documentation](../../../doc/README.md) - General BELFEM documentation
