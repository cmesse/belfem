# BELFEM Linear Algebra Module - Usage Guide {#linalg_linalg_usage_guide}

**Date:** 2026-01-16
**Module:** linalg
**Purpose:** Comprehensive guide to BELFEM's backend-agnostic linear algebra API
**Revision:** 2026-01-16 - Initial version

---

## Overview

The `src/linalg` module provides a **backend-agnostic** linear algebra API for BELFEM. It abstracts over two high-performance backends:

- **Armadillo** (default on Linux): LAPACK/BLAS-based, mature ecosystem
- **Blaze**: Expression template library, optimized for small-to-medium matrices

The abstraction allows switching backends at **compile time** without changing user code, enabling performance tuning and portability.

### Key Features

1. **Backend Independence**: Write once, compile with either Armadillo or Blaze
2. **Expression Templates**: Lazy evaluation for complex expressions (e.g., `A*B + C*D`)
3. **Bounds Checking**: Debug-mode assertions (zero overhead in release)
4. **Operator Overloading**: Natural mathematical syntax (`+`, `-`, `*`, `/`)
5. **Direct Backend Access**: Escape hatch via `vector_data()` / `matrix_data()`

> **Critical: Backend Selection**
> Backend is chosen at **compile time**. The CMake cache options are
> `USE_MATRIX_ARMADILLO` and `USE_MATRIX_BLAZE` (`CMakeLists.txt:81-87`); the
> `BELFEM_ARMADILLO` / `BELFEM_BLAZE` names are the **generated compile defines**
> (`config/linalg/config_matrix.cmake`), so `-DBELFEM_ARMADILLO=ON` on the command line does
> nothing. The default is **Armadillo everywhere except Apple, where it is Blaze**.
> You cannot mix backends in a single build.

---

## Common Pitfalls

### 1. **Mixing Vector/Matrix with Raw Backend Types**

```cpp
Vector<real> v(3);
arma::vec armadillo_vec = v;  // ERROR: type mismatch

// Correct:
arma::vec armadillo_vec = v.vector_data();  // Explicit conversion
```

**Solution:** Use `.vector_data()` or `.matrix_data()` for backend interop.

### 2. **Element-wise vs. Matrix Multiplication**

```cpp
Matrix<real> A(3, 3), B(3, 3);
Matrix<real> C = A * B;   // Matrix multiplication
// For element-wise: use backend-specific syntax or implement wrapper
```

**Solution:** `*` is matrix multiplication. For element-wise, access backend directly.

### 3. **Initializer List Ambiguity**

```cpp
Vector<real> v = {3};     // Single element: v = [3]
Vector<real> w(3);        // Three zeros: w = [0, 0, 0]
Vector<real> u = {1,2,3}; // Three elements: u = [1, 2, 3]
```

**Solution:** Use explicit syntax: `Vector<real> v(3, 0.0)` for size+fill.

### 4. **Row/Column Lifetime**

```cpp
Matrix<real> A(3, 3);
auto row = A.row(0);      // Returns backend view (Armadillo: subview_row)
A.set_size(5, 5);         // row is now INVALID (dangling reference)
```

**Solution:** Don't hold views across resizing operations.

### 5. **Debug vs. Release Behavior**

```cpp
Vector<real> v(5);
real x = v(10);  // Debug: Assertion "Index 10 out of bounds"
                 // Release: Undefined behavior
```

**Solution:** Test in both modes. Use static analysis for release builds.

---

## Core Types

### Vector<T> - Column Vector

**Files:** `cl_Vector.hpp`, backend implementations in `armadillo/cl_AR_Vector.hpp` or `blaze/cl_BZ_Vector.hpp`

#### Description

Template class wrapping backend column vector (stored as n×1 matrix in Armadillo).

#### Construction

```cpp
// Empty vector
Vector<real> v;

// Size only (uninitialized)
Vector<real> v(10);

// Size + fill value
Vector<real> v(10, 0.0);

// Initializer list
Vector<real> v = {1.0, 2.0, 3.0};
Vector<int> v = {1, 2, 3, 4, 5};

// From backend expression (advanced)
Vector<real> v = trans(someMatrix.row(0));
```

#### Element Access

```cpp
Vector<real> v(5, 1.0);

// Bounds-checked access (debug mode)
real x = v(0);              // First element
v(4) = 2.0;                 // Last element

// Raw pointer (for C/Fortran/MPI)
real* ptr = v.data();
const real* cptr = v.data();  // const version
```

#### Size Operations

```cpp
size_t n = v.length();      // Number of elements
v.set_size(20);             // Resize (no fill)
v.set_size(20, 0.0);        // Resize + fill
v.fill(1.0);                // Fill with value
```

#### Arithmetic Operations

```cpp
Vector<real> a(3), b(3), c(3);

// Assignment
a = 1.0;                    // Fill with scalar
a = {1, 2, 3};              // From initializer list
a = b;                      // Copy
a = std::move(b);           // Move

// In-place operations
a += b;                     // Element-wise addition
a -= b;                     // Element-wise subtraction
a *= 2.0;                   // Scalar multiplication
a /= 2.0;                   // Scalar division
a %= b;                     // Element-wise multiplication (Armadillo .%, Blaze %)

// Binary operations (return new vector)
c = a + b;
c = a - b;
c = a * 2.0;
c = 2.0 * a;
c = a / 2.0;
```

#### Iteration

```cpp
Vector<real> v(5, 1.0);
for (real& x : v) {
    x *= 2.0;
}

for (const real& x : v) {
    std::cout << x << " ";
}
```

#### Backend Interoperability

```cpp
Vector<real> v(3);

// Get underlying backend type
// Armadillo: arma::Mat<real> (column vectors are n×1 matrices)
// Blaze: blaze::DynamicVector<real>
auto& backend_vec = v.vector_data();

// Use backend-specific features
#ifdef BELFEM_ARMADILLO
    arma::Mat<real>& av = v.vector_data();   // n×1 Mat, not arma::vec
    av.save("vector.dat", arma::raw_ascii);
#elif BELFEM_BLAZE
    blaze::DynamicVector<real>& bv = v.vector_data();
    // Blaze-specific operations
#endif
```

#### Printing

```cpp
Vector<real> v = {1.0, 2.0, 3.0};
v.print("MyVector");
// Output:
// MyVector = [ ...
// +1.000000000000000e+00; ...
// +2.000000000000000e+00; ...
// +3.000000000000000e+00; ];
```

#### When to Use

- Primary vector type for all BELFEM computations
- FEM shape function evaluations
- Nodal coordinates, DOF values
- Right-hand side vectors for linear systems

**See:** `Vector` class in `cl_Vector.hpp` and backend implementations

---

### Matrix<T> - Dense Matrix

**Files:** `cl_Matrix.hpp`, backend implementations in `armadillo/cl_AR_Matrix.hpp` or `blaze/cl_BZ_Matrix.hpp`

#### Description

Template class wrapping backend dense matrix (column-major in both backends).

#### Construction

```cpp
// Empty matrix
Matrix<real> A;

// Size only (uninitialized)
Matrix<real> A(3, 3);

// Size + fill value
Matrix<real> A(3, 3, 0.0);

// Initializer list
Matrix<real> A = {{1, 2, 3},
                  {4, 5, 6},
                  {7, 8, 9}};

// From backend expression
Matrix<real> A = inv(B);
```

#### Element Access

```cpp
Matrix<real> A(3, 3);

// Bounds-checked access (debug mode)
real x = A(0, 0);           // Top-left element
A(2, 2) = 1.0;              // Bottom-right element

// Row/column access (returns backend view)
auto row0 = A.row(0);       // View of row 0
auto col1 = A.col(1);       // View of column 1

// Submatrix (inclusive indices)
auto sub = A.submat(0, 0, 1, 1);  // Top-left 2×2 block

// Raw pointer (column-major in both backends)
real* ptr = A.data();
```

#### Size Operations

```cpp
size_t rows = A.n_rows();
size_t cols = A.n_cols();
size_t total = A.capacity();  // backend capacity -- NOT n_rows * n_cols in general

A.set_size(5, 5);             // Resize (no fill)
A.set_size(5, 5, 0.0);        // Resize + fill
A.fill(1.0);                  // Fill with value
```

> **Do not size a raw transfer from `capacity()`.** Under Blaze the columns are padded for SIMD
> alignment, so `spacing()` — the inter-column stride — can exceed `n_rows()` and `capacity()`
> is the padded total (`cl_BZ_Matrix.hpp:254-269`). For an MPI send of a whole matrix the
> payload length is **`spacing() * n_cols()`**. Element access always goes through `A(i, j)`;
> never compute an offset from `data()`.

#### Setting Rows/Columns

```cpp
Matrix<real> A(3, 3);
Vector<real> v = {1, 2, 3};

A.set_row(0, v);              // Set row 0 to v
A.set_col(1, v);              // Set column 1 to v
```

#### Arithmetic Operations

```cpp
Matrix<real> A(3, 3), B(3, 3), C(3, 3);

// Assignment
A = 1.0;                      // Fill with scalar
A = B;                        // Copy
A = std::move(B);             // Move

// In-place operations
A += B;                       // Element-wise addition
A -= B;                       // Element-wise subtraction
A *= 2.0;                     // Scalar multiplication
A /= 2.0;                     // Scalar division
A *= B;                       // Matrix multiplication (A = A * B)

// Binary operations (return new matrix)
C = A + B;                    // Element-wise addition
C = A - B;                    // Element-wise subtraction
C = A * B;                    // Matrix multiplication
C = A * 2.0;                  // Scalar multiplication
C = A / 2.0;                  // Scalar division
```

#### Matrix-Vector Multiplication

```cpp
Matrix<real> A(3, 3);
Vector<real> x(3), b(3);

// Implicit via operator overload (defined in operators/)
b = A * x;                    // Matrix-vector product
```

#### Printing

```cpp
Matrix<real> A = {{1, 2}, {3, 4}};
A.print("MyMatrix");
// Output:
// MyMatrix = [ ...
// +1.000000000000000e+00, +2.000000000000000e+00; ...
// +3.000000000000000e+00, +4.000000000000000e+00; ];
```

#### When to Use

- Jacobians, stiffness matrices (small, dense)
- Rotation matrices, transformation matrices
- Local element matrices before assembly

**See:** `Matrix` class in `cl_Matrix.hpp` and backend implementations

---

## Free Functions

### Linear Algebra Operations

#### `dot(a, b)` - Dot Product

```cpp
Vector<real> a = {1, 2, 3};
Vector<real> b = {4, 5, 6};
real result = dot(a, b);    // 1*4 + 2*5 + 3*6 = 32
```

**File:** `fn_dot.hpp`

---

#### `cross(a, b)` - Cross Product

```cpp
Vector<real> a = {1, 0, 0};
Vector<real> b = {0, 1, 0};
Vector<real> c = cross(a, b);  // [0, 0, 1]
```

**Requirements:** Both vectors must have length 3.

**File:** `fn_cross.hpp`

---

#### `crossmat(n, A, out)` - Normal-Matrix Cross Product

Computes cross product of a normal vector with each column of a matrix.

```cpp
// 2D version: out[k] = n × A[:,k]
Vector<real> n = {nx, ny};
Matrix<real> A(2, ncols);
Vector<real> nxA(ncols);
crossmat(n, A, nxA);        // Each element is 2D cross product

// 3D version: out[:,k] = n × A[:,k]
Vector<real> n3 = {nx, ny, nz};
Matrix<real> A3(3, ncols);
Matrix<real> nxA3(3, ncols);
crossmat(n3, A3, nxA3);     // Each column is 3D cross product
```

**Signature:**
- 2D: `void crossmat(const Vector<real>& n, const Matrix<real>& A, Vector<real>& out)`
- 3D: `void crossmat(const Vector<real>& n, const Matrix<real>& A, Matrix<real>& out)`

**Use case:** FEM flux computations, boundary integrals

**File:** `fn_crossmat.hpp`

---

#### `norm(v)` - Euclidean Norm

```cpp
Vector<real> v = {3, 4};
real n = norm(v);           // 5.0
```

**File:** `fn_norm.hpp`

---

#### `trans(A)` - Transpose

```cpp
Matrix<real> A = {{1, 2, 3},
                  {4, 5, 6}};
Matrix<real> AT = trans(A);  // 3×2 matrix
```

**File:** `fn_trans.hpp`

---

### Matrix Decompositions and Solvers

#### `inv(A)` - Matrix Inverse

```cpp
Matrix<real> A = {{4, 7}, {2, 6}};
Matrix<real> Ainv = inv(A);
// Ainv ≈ [[0.6, -0.7], [-0.2, 0.4]]
```

**Note:** For 2×2 and 3×3, specialized fast implementations `inv2()` and `inv3()` exist.

**Files:** `fn_inv.hpp`, `fn_inv2.hpp`, `fn_inv3.hpp`

---

#### `det(A)` - Determinant

```cpp
Matrix<real> A = {{1, 2}, {3, 4}};
real d = det(A);            // -2.0
```

**File:** `fn_det.hpp`

---

#### `gesv(A, x, pivot)` - General Linear System Solve

Solves Ax = b via LU decomposition. **Overwrites A with LU factorization and x with solution.**

```cpp
Matrix<real> A = {{3, 2}, {1, 2}};
Vector<real> x = {5, 3};     // Initialize x with RHS b
Vector<int_t> pivot(2);        // Pivot indices (required, size ≥ n)
gesv(A, x, pivot);           // A is overwritten, x now holds solution [1, 1]
```

**Signature:** `int_t gesv(Matrix<T>& A, Vector<T>& x, Vector<int_t>& pivot, bool AbortOnError = true)` — returns LAPACK `info` (see the LAPACK guide)
**Modifies:** Both `A` (LU factorization) and `x` (RHS → solution)
**Requires:** `pivot.length() >= A.n_rows()`
**File:** `fn_gesv.hpp`

---

#### `posv(A, x)` - Positive Definite System Solve

Solves Ax = b assuming A is symmetric positive definite (uses Cholesky). **Overwrites A with Cholesky factorization and x with solution.**

```cpp
Matrix<real> A = {{4, 2}, {2, 3}};  // SPD matrix
Vector<real> x = {6, 5};            // Initialize x with RHS b
posv(A, x);                         // A overwritten, x now holds solution
```

**Signature:** `int_t posv(Matrix<T>& A, Vector<T>& x, bool AbortOnError = true)`
**Returns:** LAPACK `info` — 0 on success; > 0 means the leading minor of that order is not
positive definite. With `AbortOnError = true` (the default) a non-zero `info` raises `BELFEM_ERROR` instead — which throws in a debug build and aborts in release.
**Modifies:** Both `A` (Cholesky factorization) and `x` (RHS → solution)
**Requires:** A must be symmetric positive definite
**Faster than `gesv()` (~2× for SPD matrices)**

**File:** `src/linalg/lapack/fn_posv.hpp`

---

#### `eigen(A, values [, abortOnComplex])` / `eigen_sym(A, values)` - Eigenvalues

Eigenvalues only (no eigenvectors). `eigen` returns the number of complex eigenvalues
(aborts on the first one by default); `eigen_sym` reads the upper triangle and returns
ascending real eigenvalues.

```cpp
Matrix<real> A = {{2, 1}, {1, 2}};
Vector<real> values;
eigen_sym(A, values);  // values = [1, 3]
```

**File:** `fn_eigen.hpp`

---

### Utility Functions

#### `linspace(start, end, n)` - Linearly Spaced Vector

```cpp
Vector<real> v = linspace(0.0, 1.0, 11);
// v = [0.0, 0.1, 0.2, ..., 0.9, 1.0]
```

**File:** `fn_linspace.hpp`

---

#### `append(a, b)` - Vector Concatenation (In-Place)

```cpp
Vector<real> a = {1, 2};
Vector<real> b = {3, 4};
append(a, b);               // a is now [1, 2, 3, 4], b unchanged
```

**Signature:** `void append(Vector<T>& a, Vector<T>& b)`
**Modifies:** First argument `a` in-place
**File:** `fn_append.hpp`

---

#### `combine(a, b, out)` - Concatenate Vectors

Writes `a` followed by `b` into `out` (resized; must not alias an input). Overloads
take three or four inputs. Unlike `append`, the inputs are untouched.

**File:** `fn_combine.hpp`

---

#### `reverse(v)` - Reverse Vector

```cpp
Vector<real> v = {1, 2, 3};
Vector<real> r = reverse(v);  // [3, 2, 1]
```

**File:** `fn_reverse.hpp`

---

#### `sort(v)` - Sort Vector (In-Place)

```cpp
Vector<real> v = {3, 1, 2};
sort(v);                    // v is now [1, 2, 3]
```

**Signature:** `void sort(Vector<T>& v)`
**Modifies:** Argument in-place
**File:** `fn_sort.hpp`

---

#### `unique(v)` - Unique Elements (In-Place)

Sorts vector and removes duplicates in-place.

```cpp
Vector<real> v = {1, 2, 2, 3, 1};
unique(v);                   // v is now [1, 2, 3]
```

**Signature:** `void unique(Vector<T>& v)`
**Modifies:** Argument in-place (sorts then removes duplicates)
**File:** `fn_unique.hpp`

---

### Statistical/Aggregate Functions

#### `sum(v)` - Sum of Elements

```cpp
Vector<real> v = {1, 2, 3, 4};
real s = sum(v);  // 10.0
```

**File:** `fn_sum.hpp`

---

#### `max(v)` - Maximum Element

```cpp
Vector<real> v = {1, 5, 3};
real m = max(v);  // 5.0
```

**File:** `fn_max.hpp`

---

#### `min(v)` - Minimum Element

```cpp
Vector<real> v = {1, 5, 3};
real m = min(v);  // 1.0
```

**File:** `fn_min.hpp`

---

### Polynomial Operations

#### `polyval(p, x)` - Evaluate Polynomial

Evaluates polynomial using **highest-power-first** coefficient ordering via Horner's method.

```cpp
// Coefficient convention: p[0]*x^n + p[1]*x^(n-1) + ... + p[n]
Vector<real> p = {3, 2, 1};  // 3x^2 + 2x + 1
real y = polyval(p, 2.0);    // 3*4 + 2*2 + 1 = 17
```

> **Important:** Coefficients are **highest-power-first**, not lowest-power-first.

**File:** `fn_polyval.hpp`

---

#### `dpolyval(p, x)` - Polynomial Derivative

First derivative: dp/dx.

**File:** `fn_dpolyval.hpp`

---

#### `ddpolyval(p, x)` - Polynomial Second Derivative

Second derivative: d²p/dx².

**File:** `fn_ddpolyval.hpp`

---

#### `polyfit(x, y, degree, coeffs)` - Polynomial Fit

Least-squares polynomial fit of given degree. The coefficient vector is the last
argument; it is resized to degree + 1, highest power first.

```cpp
Vector<real> x = {0, 1, 2, 3};
Vector<real> y = {1, 3, 7, 13};  // Roughly y = 1 + 2x²
Vector<real> p;
polyfit(x, y, 2, p);
```

**File:** `fn_polyfit.hpp`

---

### Specialized Functions

#### `r2(approximated, exact)` - Coefficient of Determination

Statistical measure of fit quality (R² score).

```cpp
Vector<real> y_true = {1, 2, 3, 4};
Vector<real> y_pred = {1.1, 1.9, 3.2, 3.9};
real r2 = r2(y_pred, y_true);  // ~0.98; the exact data is the SECOND argument
```

**File:** `fn_r2.hpp`

---

## Backend Comparison

| Feature | Armadillo | Blaze |
|---------|-----------|-------|
| **Maturity** | Mature (10+ years) | Mature (5+ years) |
| **Dependencies** | LAPACK/BLAS | Header-only |
| **Expression Templates** | Yes | Yes (advanced) |
| **Sparse Matrices** | Yes | Yes |
| **Small Matrix Opt** | Good | Excellent |
| **Large Matrix Opt** | Excellent (LAPACK) | Good |
| **Default in BELFEM** | **Yes** | No |
| **Ease of Installation** | Moderate (needs LAPACK) | Easy (header-only) |

### When to Use Which Backend

**Use Armadillo (Linux default) when:**
- Need sparse matrix support
- Large dense matrices (> 100×100)
- Interfacing with existing LAPACK code
- Need mature, stable ecosystem

**Use Blaze when:**
- Small to medium matrices (< 50×50)
- Header-only build preferred
- Maximum performance for expression templates
- Modern C++ template metaprogramming

**Recommendation:** Stick with **Armadillo** unless you have specific performance requirements.

---

## Backend Selection at Compile Time

### CMake Configuration

```cmake
# Armadillo (default on Linux; Blaze is the default on Apple)
cmake -DUSE_MATRIX_ARMADILLO=ON -DUSE_MATRIX_BLAZE=OFF ..

# Blaze
cmake -DUSE_MATRIX_BLAZE=ON -DUSE_MATRIX_ARMADILLO=OFF ..
```

### Preprocessor Checks

```cpp
#ifdef BELFEM_ARMADILLO
    // Armadillo-specific code
    arma::mat A = ...;
#elif BELFEM_BLAZE
    // Blaze-specific code
    blaze::DynamicMatrix<real> A = ...;
#endif
```

---

## Performance Tips

### 1. **Use Expression Templates**

```cpp
// Bad: creates temporary matrices
Matrix<real> temp1 = A * B;
Matrix<real> temp2 = temp1 + C;
Matrix<real> result = temp2 * D;

// Good: single operation via expression templates
Matrix<real> result = (A * B + C) * D;
```

### 2. **Reserve Size When Known**

```cpp
// Bad: resize during loop
Vector<real> v;
for (int i = 0; i < 1000; ++i) {
    v.set_size(i+1);  // Reallocates every iteration
}

// Good: allocate once
Vector<real> v(1000);
for (int i = 0; i < 1000; ++i) {
    v(i) = compute(i);
}
```

### 3. **Prefer Specialized Inverses**

```cpp
// For 2×2 matrices
Matrix<real> A(2, 2), Ainv(2, 2);
real detA = inv2(A, Ainv);  // Ainv filled in place, determinant returned; faster than inv(A)

// For 3×3 matrices
Matrix<real> B(3, 3), Binv(3, 3);
real detB = inv3(B, Binv);  // Faster than inv(B)
```

### 4. **Use SPD Solvers When Applicable**

```cpp
// If A is symmetric positive definite
// posv solves in place and returns LAPACK info, not the solution
Vector<real> x = b;
posv(A, x);                   // faster than gesv() for SPD; x now holds the solution
```

### 5. **Avoid Row/Column Copies**

```cpp
// Bad: copies row
Vector<real> row = A.row(0);

// Good: use view (but watch lifetime!)
auto row_view = A.row(0);
real sum = 0.0;
for (size_t i = 0; i < A.n_cols(); ++i) {
    sum += row_view(i);
}
```

---

## Common Patterns

### Pattern 1: Solving Linear System

```cpp
Matrix<real> K(n, n);     // Stiffness matrix
Vector<real> u(n);        // Initialize with force vector f
Vector<int_t> pivot(n);     // Pivot indices

// ... assemble K and u (u contains RHS) ...

gesv(K, u, pivot);        // K overwritten, u now contains solution
```

### Pattern 2: Rotation Matrix Construction

```cpp
real theta = M_PI / 4;  // 45 degrees
Matrix<real> R(2, 2);
R(0, 0) = cos(theta);
R(0, 1) = -sin(theta);
R(1, 0) = sin(theta);
R(1, 1) = cos(theta);

Vector<real> v = {1, 0};
Vector<real> v_rot = R * v;  // Rotated vector
```

### Pattern 3: Coordinate Transformation

```cpp
// Transform from local to global coordinates
Matrix<real> T = compute_transformation_matrix();
Vector<real> local_coords = {1, 2, 3};
Vector<real> global_coords = T * local_coords;
```

### Pattern 4: Polynomial Interpolation

```cpp
// Fit polynomial to data
Vector<real> x = {0, 1, 2, 3, 4};
Vector<real> y = {1, 2.1, 3.9, 6.2, 8.8};
Vector<real> p;
polyfit(x, y, 2, p);  // Quadratic fit

// Evaluate at new points
real y_interp = polyval(p, 2.5);
```

### Pattern 5: Jacobian Computation

```cpp
// Numerical Jacobian via finite differences
Matrix<real> jacobian(m, n);
Vector<real> x0(n);
real h = 1e-8;

for (size_t j = 0; j < n; ++j) {
    Vector<real> x_plus = x0;
    x_plus(j) += h;
    Vector<real> f_plus = residual(x_plus);

    Vector<real> x_minus = x0;
    x_minus(j) -= h;
    Vector<real> f_minus = residual(x_minus);

    Vector<real> df = (f_plus - f_minus) / (2*h);
    jacobian.set_col(j, df);
}
```

---

## Thread Safety and MPI

**Thread Safety:** Like containers, Vector and Matrix are **not thread-safe**.

- **Reading:** Safe from multiple threads
- **Writing:** Requires external synchronization

**MPI:** Each rank has independent copies. Use `comm` module for synchronization.

**Example with OpenMP:**

```cpp
Matrix<real> A(100, 100);

#pragma omp parallel for collapse(2)
for (int i = 0; i < 100; ++i) {
    for (int j = 0; j < 100; ++j) {
        // Safe: each thread writes different element
        A(i, j) = compute(i, j);
    }
}
```

---

## Sparse Matrices

BELFEM's linalg module focuses on **dense** matrices. For **sparse** matrices, use:

- **`src/sparse/`** module: Wrappers for MUMPS, STRUMPACK, PETSc, etc.
- **Direct backend access:** `arma::sp_mat` (Armadillo) or `blaze::CompressedMatrix` (Blaze)

---

## Related Modules

- **containers**: `Cell<T>` for dynamic arrays (complements Vector for non-numeric data)
- **sparse**: Sparse matrix solvers (MUMPS, STRUMPACK, PETSc, PARDISO)
- **fem**: Uses Vector/Matrix extensively for element matrices
- **numerics**: Polynomial evaluation, splines, integration (uses linalg)

---

## See Also

- [Armadillo Documentation](http://arma.sourceforge.net/docs.html)
- [Blaze Documentation](https://bitbucket.org/blaze-lib/blaze/wiki/Home)
- [Containers Module](../../containers/doc/README.md) - BELFEM container classes
- [Sparse Module](../../sparse/doc/README.md) - Sparse matrix solvers
- `CLAUDE.md` (repository root) - Documentation guidelines
- C++ Documentation: `make doc` (Doxygen)

---

**Revision History:**
- 2026-01-16: Initial version
- 2026-01-16: Corrected factual errors based on external review (in-place function signatures, polynomial conventions, storage order claims, crossmat semantics)
