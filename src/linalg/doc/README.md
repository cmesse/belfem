# Linear Algebra Module Documentation {#linalg_index}

Documentation for BELFEM's backend-agnostic linear algebra API.

## Contents

### Guides

- [**linalg_usage_guide.md**](linalg_usage_guide.md) - Comprehensive usage guide for BELFEM linear algebra
  - Overview of backend abstraction (Armadillo vs. Blaze)
  - Vector and Matrix API reference
  - 26 free functions (solvers, decompositions, utilities)
  - Performance tips and common patterns
  - Backend selection and interoperability
- [**lapack_usage_guide.md**](lapack_usage_guide.md) - Maintainer guide for the unified LAPACK wrappers (`src/linalg/lapack/`)
  - gesv, posv, getrf/getri, gels, gemm, gesvd, geev, gees — one body, both backends, four datatypes
  - Shared contracts: AbortOnError/info error model, pivot policy, the three Work-buffer models
  - Design pillars: int_t = LAPACK integer, complex glue, leading dimensions, hidden Fortran char lengths
  - Checklist for adding a new routine, pitfalls table

## Quick Reference

### Core Types

| Type | Description | File |
|------|-------------|------|
| `Vector<T>` | Column vector (n×1) | `cl_Vector.hpp` |
| `Matrix<T>` | Dense matrix (m×n) | `cl_Matrix.hpp` |

### Categories of Free Functions

| Category | Functions |
|----------|-----------|
| **Linear Algebra** | `dot`, `cross`, `crossmat`, `norm`, `trans` |
| **Solvers** | `gesv`, `posv`, `inv`, `inv2`, `inv3`, `det`, `eigen` |
| **Utilities** | `linspace`, `append`, `combine`, `reverse`, `sort`, `unique` |
| **Aggregates** | `sum`, `max`, `min` |
| **Polynomials** | `polyval`, `dpolyval`, `ddpolyval`, `polyfit` |
| **Statistics** | `r2` |

### Backend Selection

```cmake
# Armadillo (default on Linux; Blaze is the default on Apple)
cmake -DUSE_MATRIX_ARMADILLO=ON -DUSE_MATRIX_BLAZE=OFF ..

# Blaze
cmake -DUSE_MATRIX_BLAZE=ON -DUSE_MATRIX_ARMADILLO=OFF ..
```

## See Also

- **Source code** - header files in `src/linalg/`
- `Armadillo backend` - Armadillo-specific implementations
- `Blaze backend` - Blaze-specific implementations
- [Sparse module](../../sparse/doc/README.md) - Sparse matrix solvers
- [Root documentation](../../../doc/README.md) - General BELFEM documentation
