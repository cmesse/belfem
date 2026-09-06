# BELFEM Tensor Tests — Detailed Plan

**Date:** 2026-03-22
**Purpose:** Method-level test matrix for the tensor module (`src/math/` or `src/tensor/`)
**Depends on:** `tests_0_strategy.md` (conventions), `tests_2_linalg.md` (Vector/Matrix assumed tested)
**Confidence:** High on API surface. High on physics identity approach. Medium on exact Voigt convention (verify against literature during implementation).

---

## Module Overview

A 4th-order tensor library for continuum mechanics / elasticity. All meaningful operations are restricted to 3×3×3×3 tensors (81 elements, column-major flat storage). The module is split into three layers:

| Layer | Content |
|---|---|
| `Tensor<T>` class | Storage, access, fill, operators, conversion, contraction (member API) |
| `belfem::tensor::*` free functions | `identity`, `fill`, `add`, `subtract`, `contract42`, `contract44`, `rotate42`, `kelvin_christoffel`, `mat_to_ten`, `ten_to_mat`, `equal_equal`, `invert_symmetric` |
| High-level API | `ddot()`, `kelvin_christoffel()`, `rotate()`, `compliance_matrix()`, `fiber_polarization()` |

**Backend note:** `contract42`, `contract44`, `kelvin_christoffel`, and `rotate42` have separate Armadillo and Blaze implementations (selected at compile time). Both are hand-unrolled pointer-arithmetic kernels. The same tests validate whichever backend is active.

**All operations use `BELFEM_ASSERT` (debug only).** No `BELFEM_ERROR` in this module.

---

## Testing Philosophy

This module is best tested through **invariants and reference implementations**, not method counting. The hand-unrolled kernels (thousands of lines of explicit index arithmetic) are the primary risk — a single wrong index produces subtly incorrect material response. The strongest defense is:

1. **Naive reference loops** — slow but obviously correct nested-loop implementations of contractions and rotations, compared against the unrolled versions.
2. **Physics identities** — rotation invariance, identity contraction, round-trip conversions, symmetry preservation.
3. **Known analytical solutions** — isotropic elasticity with textbook Lamé parameters.

---

## Floating-Point Comparison

```cpp
namespace
{
    const belfem::real tEps = 1e-12 ;  // for exact-in-theory results
    const belfem::real tTol = 1e-9 ;   // for results through inversion / many FP ops
}
```

**Important:** `operator==` on `Tensor` is **exact elementwise comparison** (via `tensor::equal_equal`). Use it only for copied/filled tensors. For anything involving arithmetic, use entrywise `EXPECT_NEAR`.

---

## Reference Loop Helpers

Claude Code should implement these as test-local helper functions — they are slow but correct by inspection:

```cpp
// Naive A_ijkl * B_kl = C_ij  (contraction of rank-4 with rank-2)
void ref_contract42( const real * A, const real * B, real * C )
{
    for( int i = 0; i < 3; ++i )
    for( int j = 0; j < 3; ++j )
    {
        real tSum = 0.0;
        for( int k = 0; k < 3; ++k )
        for( int l = 0; l < 3; ++l )
            tSum += A[ l*27 + k*9 + j*3 + i ] * B[ l*3 + k ];
        C[ j*3 + i ] = tSum;
    }
}

// Naive A_ijmn * B_mnkl = C_ijkl  (double contraction of two rank-4)
void ref_contract44( const real * A, const real * B, real * C )
{
    for( int i = 0; i < 3; ++i )
    for( int j = 0; j < 3; ++j )
    for( int k = 0; k < 3; ++k )
    for( int l = 0; l < 3; ++l )
    {
        real tSum = 0.0;
        for( int m = 0; m < 3; ++m )
        for( int n = 0; n < 3; ++n )
            tSum += A[ n*27 + m*9 + j*3 + i ] * B[ l*27 + k*9 + n*3 + m ];
        C[ l*27 + k*9 + j*3 + i ] = tSum;
    }
}

// Naive A_mnop = B_ijkl * R_im * R_jn * R_ko * R_lp  (tensor rotation)
void ref_rotate42( const real * B, const real * R, real * A )
{
    std::fill( A, A + 81, 0.0 );
    for( int m = 0; m < 3; ++m )
    for( int n = 0; n < 3; ++n )
    for( int o = 0; o < 3; ++o )
    for( int p = 0; p < 3; ++p )
    for( int i = 0; i < 3; ++i )
    for( int j = 0; j < 3; ++j )
    for( int k = 0; k < 3; ++k )
    for( int l = 0; l < 3; ++l )
        A[ p*27 + o*9 + n*3 + m ] += B[ l*27 + k*9 + j*3 + i ]
            * R[ m*3 + i ] * R[ n*3 + j ] * R[ o*3 + k ] * R[ p*3 + l ];
            // NOTE: verify column-major layout of R matches Matrix<real>::data()
}

// Naive C_ik = A_ijkl * n_j * n_l  (Kelvin-Christoffel)
void ref_kelvin_christoffel( const real * A, const real * n, real * C )
{
    for( int i = 0; i < 3; ++i )
    for( int k = 0; k < 3; ++k )
    {
        real tSum = 0.0;
        for( int j = 0; j < 3; ++j )
        for( int l = 0; l < 3; ++l )
            tSum += A[ l*27 + k*9 + j*3 + i ] * n[ j ] * n[ l ];
        C[ k*3 + i ] = tSum;
    }
}
```

**CRITICAL:** The index layout in these reference loops must match BELFEM's flat storage order: `data[L*27 + K*9 + J*3 + I]` for `T(I,J,K,L)`. Verify by checking `Tensor::operator()` before writing the first test.

---

## Test File Structure

```
tests/math/
├── test_Tensor.cpp                # Class semantics, conversion, operators, fill
├── test_TensorKernels.cpp         # Contraction, rotation, Kelvin-Christoffel (reference-checked)
```

---

## 1. Tensor\<T\> Class Semantics

**File:** `test_Tensor.cpp`

### 1.1 Construction & Destruction `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `SizedConstructor` | `Tensor(3,3,3,3)` → `capacity() == 81`, `is_3333() == true` |
| `SizedConstructorNon3333` | `Tensor(2,2,2,2)` → `capacity() == 16`, `is_3333() == false` |
| `FillValueConstructor` | `Tensor(3,3,3,3, 7.0)` → all 81 elements are `7.0` |
| `FromElasticityMatrix` | Construct from known 6×6 matrix → verify specific tensor entries |
| `CopyConstructorDeepCopies` | Modify copy, original unchanged |
| `MoveConstructorTransfers` | Source `data() == nullptr` after move |
| `CopyAssignment` | Deep copy, verify via `operator==` |
| `MoveAssignment` | Source `data() == nullptr`, target has data |
| `ScalarAssignment` | `tTen = 5.0` fills all 81 entries with `5.0` |
| `MatrixAssignment` | `tTen = tElasticityMatrix` populates from 6×6 |

### 1.2 Construction `[debug]`

| Test Name | What It Verifies |
|---|---|
| `FromMatrixWrongSizeThrows` | 6×6 constructor with non-6×6 matrix → assertion |
| `CopyAssignmentSizeMismatchThrows` | Assign 3×3×3×3 to 2×2×2×2 → assertion |
| `MoveAssignmentSizeMismatchThrows` | Same for move |
| `MatrixAssignmentNon3333Throws` | `tTen(2,2,2,2) = tMatrix` → assertion |

### 1.3 Access `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `ParenthesisReadWrite` | `tTen(i,j,k,l) = x` then `tTen(i,j,k,l) == x` |
| `FlatLayoutMatchesOperator` | `data()[L*27 + K*9 + J*3 + I] == operator()(I,J,K,L)` for all I,J,K,L |
| `DataPointerNonNull` | `data()` is non-null after construction |
| `CapacityMatchesProduct` | `capacity() == I*J*K*L` |
| `Is3333` | True for (3,3,3,3), false for others |

### 1.4 Access `[debug]`

| Test Name | What It Verifies |
|---|---|
| `IndexIOutOfBoundsThrows` | `tTen(3,0,0,0)` on 3×3×3×3 → assertion |
| `IndexJOutOfBoundsThrows` | `tTen(0,3,0,0)` → assertion |
| `IndexKOutOfBoundsThrows` | `tTen(0,0,3,0)` → assertion |
| `IndexLOutOfBoundsThrows` | `tTen(0,0,0,3)` → assertion |

### 1.5 Operators `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `PlusEqualsScalar` | Each element increased by scalar |
| `MinusEqualsScalar` | Each element decreased by scalar |
| `TimesEqualsScalar` | Each element scaled |
| `DivideEqualsScalar` | Each element divided |
| `PlusEqualsTensor` | Elementwise addition |
| `MinusEqualsTensor` | Elementwise subtraction |
| `BinaryPlus` | `A + B` matches componentwise |
| `BinaryMinus` | `A - B` matches componentwise |
| `EqualityExact` | Copy of tensor is `== original` |
| `EqualityDifferent` | Modified tensor `!= original` |

### 1.6 Operators `[debug]`

| Test Name | What It Verifies |
|---|---|
| `PlusEqualsTensorNon3333Throws` | Non-3×3×3×3 tensors → assertion |
| `BinaryPlusNon3333Throws` | Same for binary `+` |

---

## 2. Constitutive Fill & Conversion

**File:** `test_Tensor.cpp` (same file, separate test group)

### 2.1 Isotropic Fill `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `FillAbProducesCorrectStructure` | `fill(a, b)` → diagonal entries are `a + 4b/3`, off-diagonal coupling is `a - 2b/3`, shear entries are `b` |
| `FillIsotropicElasticity` | `fill_isotropic_elasticity(E, nu)` → derives K, G internally; result matches `fill(K, G)` |
| `FillIsotropicConsistentWithVoigtMatrix` | Convert to 6×6 → verify λ+2μ on normal diagonal, λ on normal off-diagonal, μ on shear diagonal |
| `FillIsotropicKnownSteel` | E=200e9, ν=0.3 → verify C11, C12, C44 match textbook values |

### 2.2 Isotropic Fill `[debug]`

| Test Name | What It Verifies |
|---|---|
| `FillAbNon3333Throws` | `fill(a,b)` on non-3×3×3×3 → assertion |
| `FillIsotropicNon3333Throws` | Same for `fill_isotropic_elasticity` |

### 2.3 Orthotropic Fill `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `FillOrthotropicProducesSymmetricMatrix` | Convert to 6×6 Voigt → matrix is symmetric |
| `FillOrthotropicIsotropicLimit` | Use E1=E2=E3=E, all ν equal, all G = E/(2(1+ν)) → matches isotropic result |

### 2.4 Identity Tensor `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `IdentityTensorValues` | I_ijkl = 0.5(δ_ik δ_jl + δ_il δ_jk) — check all 81 entries |
| `IdentityContractionWithMatrix` | `I % S ≈ S` for symmetric 3×3 matrix S |
| `IdentityContractionWithTensor` | `I % C % I ≈ C` for any C (identity is neutral for symmetric double contraction) |

### 2.5 Compliance Matrix `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `ComplianceMatrixIsotropic` | Isotropic parameters → known S_ij values |
| `ComplianceMatrixSymmetric` | Result is symmetric: `S(i,j) == S(j,i)` |
| `StiffnessTimesComplianceIsIdentity` | `inv(S) * S ≈ I(6×6)` |

### 2.6 Mat ↔ Ten Conversion `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `MatToTenToMatRoundTrip` | `mat_to_ten(C, A)` then `ten_to_mat(A, C2)` → `C2 ≈ C` for a general 6×6 |
| `TenToMatToTenRoundTrip` | `ten_to_mat(A, C)` then `mat_to_ten(C, A2)` → `A2 == A` (exact, no arithmetic) |
| `CountingTensorMapping` | Fill `A[i] = i`, call `ten_to_mat`, verify specific known positions (documents the Voigt convention) |
| `IsotropicTensorToVoigtMatrix` | Known isotropic tensor → verify the 6×6 has correct λ+2μ / λ / μ structure |
| `ConstructorFromMatrixMatchesMatToTen` | `Tensor(C6x6)` matches calling `mat_to_ten` directly |
| `ToMatrixMatchesTenToMat` | `tTen.to_matrix(M)` matches calling `ten_to_mat` directly |

### 2.7 Conversion `[debug]`

| Test Name | What It Verifies |
|---|---|
| `ToMatrixWrongSizeThrows` | `to_matrix()` with non-6×6 target → assertion |
| `ToMatrixNon3333Throws` | On non-3×3×3×3 tensor → assertion |

---

## 3. Contraction Kernels (Reference-Checked)

**File:** `test_TensorKernels.cpp`

### 3.1 contract42 — A_ijkl B_kl = C_ij `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `Contract42VsReferenceLoop` | Random-ish tensor A and matrix B: unrolled result matches `ref_contract42` entrywise |
| `Contract42IdentityTensor` | Identity tensor contracted with B → C ≈ B (for symmetric B) |
| `Contract42ZeroTensor` | Zero A → C is all zeros |
| `Contract42ZeroMatrix` | Zero B → C is all zeros |
| `Contract42IsotropicHydrostatic` | Isotropic C contracted with εI (hydrostatic strain) → σ = 3Kε I |
| `Contract42IsotropicPureShear` | Isotropic C contracted with pure shear strain → σ has only shear components |
| `Contract42MatrixOverloadMatchesPointer` | `contract42(A, Matrix, Matrix)` matches `contract42(A, B.data(), C.data())` |
| `Contract42MemberDdotMatchesFreeFunction` | `tTen.ddot(B, C)` matches `ddot(tTen, B, C)` |
| `Contract42OperatorPercentMatchesDdot` | `tTen % B` matches `ddot(tTen, B, C)` |

### 3.2 contract44 — A_ijmn B_mnkl = C_ijkl `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `Contract44VsReferenceLoop` | Random-ish tensors: unrolled result matches `ref_contract44` entrywise |
| `Contract44IdentityLeft` | `I % B ≈ B` (for B with minor symmetry) |
| `Contract44IdentityRight` | `A % I ≈ A` (for A with minor symmetry) |
| `Contract44ZeroLeft` | Zero A → C is all zeros |
| `Contract44Associativity` | `(A % B) % C ≈ A % (B % C)` for generic tensors |
| `Contract44MemberDdotMatchesFreeFunction` | `tA.ddot(tB, tC)` matches `ddot(tA, tB, tC)` |
| `Contract44OperatorPercentMatchesDdot` | `tA % tB` matches member `ddot` |

### 3.3 contract42 / contract44 `[debug]`

| Test Name | What It Verifies |
|---|---|
| `Contract42Non3333Throws` | Non-3×3×3×3 tensor → assertion |
| `Contract42WrongMatrixSizeThrows` | Non-3×3 matrix → assertion |
| `Contract44Non3333Throws` | Any non-3×3×3×3 argument → assertion |

---

## 4. Rotation Kernel

**File:** `test_TensorKernels.cpp`

### 4.1 rotate42 `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `Rotate42VsReferenceLoop` | Random-ish tensor and rotation matrix: unrolled result matches `ref_rotate42` entrywise within `tTol` |
| `Rotate42IdentityRotation` | R = I(3×3) → rotated tensor equals original |
| `Rotate42IsotropicInvariant` | Isotropic tensor rotated by arbitrary R → result ≈ original |
| `Rotate42TwoRotationsMatchComposition` | `rotate(rotate(B, R1), R2) ≈ rotate(B, R2*R1)` |
| `Rotate42SymmetryPreserved` | Input with minor symmetry → output has minor symmetry |
| `Rotate42HighLevelApiMatchesLowLevel` | `rotate(B, R, A)` matches `tensor::rotate42(B.data(), R, A.data())` |

---

## 5. Kelvin-Christoffel

**File:** `test_TensorKernels.cpp`

### 5.1 kelvin_christoffel `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `KelvinChristoffelVsReferenceLoop` | Random-ish tensor and direction: unrolled matches `ref_kelvin_christoffel` |
| `KelvinChristoffelIsotropicXDirection` | Isotropic C with n={1,0,0} → Γ is diagonal with P-wave and S-wave moduli |
| `KelvinChristoffelSymmetricForSymmetricTensor` | Minor-symmetric input → Γ is symmetric |
| `KelvinChristoffelScalingProperty` | Scaling n by s → Γ scales by s² |
| `KelvinChristoffelHighLevelMatchesLowLevel` | `kelvin_christoffel(Ten, Vec, Mat)` matches `tensor::kelvin_christoffel(...)` |

### 5.2 kelvin_christoffel `[debug]`

| Test Name | What It Verifies |
|---|---|
| `KelvinChristoffelNon3333Throws` | Non-3×3×3×3 tensor → assertion |
| `KelvinChristoffelWrongVectorLengthThrows` | Vector not length 3 → assertion |
| `KelvinChristoffelWrongMatrixSizeThrows` | Target not 3×3 → assertion |

---

## 6. Invert Symmetric

**File:** `test_TensorKernels.cpp`

### 6.1 `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `InvertSymmetricRoundTrip` | `C_inv % C ≈ I_sym` (symmetric identity tensor) |
| `InvertSymmetricIsotropic` | Known isotropic stiffness → compliance tensor matches analytical S |
| `InvertSymmetricPreservesSymmetry` | Result has minor and major symmetries |

### 6.2 `[debug]`

| Test Name | What It Verifies |
|---|---|
| `InvertSymmetricNon3333Throws` | Non-3×3×3×3 → assertion |
| `InvertSymmetricWorkVectorTooShortThrows` | `aWork.length() < 72` → assertion |
| `InvertSymmetricPivotVectorTooShortThrows` | `aPivot.length() < 36` → assertion |

---

## 7. Fiber Polarization (Lower Priority)

**File:** `test_TensorKernels.cpp`

### 7.1 `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `FiberPolarizationProduces3333Tensor` | Output is 3×3×3×3 |
| `FiberPolarizationIsotropicLimit` | Isotropic input → known analytical Eshelby result for fibers |

### 7.2 `[debug]`

| Test Name | What It Verifies |
|---|---|
| `FiberPolarizationWrongMatrixSizeThrows` | Non-6×6 elasticity matrix → assertion |
| `FiberPolarizationNon3333Throws` | Non-3×3×3×3 tensor → assertion |

---

## 8. Implementation Notes for Claude Code

1. **Write the reference loops first.** Before any test assertions, implement `ref_contract42`, `ref_contract44`, `ref_rotate42`, `ref_kelvin_christoffel` as test-local functions. They are the ground truth against which the unrolled kernels are validated.
2. **Verify flat storage layout.** Before writing reference loops, confirm that `Tensor::operator()(I,J,K,L)` indexes `data()[L*27 + K*9 + J*3 + I]` for the 3×3×3×3 case. If the stride pattern differs, adjust all reference loops accordingly.
3. **Rotation matrix layout.** The `rotate42` kernel takes `const T * R`, which is `Matrix<real>::data()` — column-major. The reference loop's index `R[col*3 + row]` must match this.
4. **"Random-ish" test data.** Do not use actual random numbers. Use deterministic small integers or simple patterns (e.g., `A[i] = sin(0.3*i + 0.7)`) so tests are reproducible. Avoid symmetric patterns that might mask index bugs.
5. **Separate exact vs approximate comparison.** Use `operator==` only for copied/filled tensors. Use entrywise `EXPECT_NEAR` for anything involving arithmetic, inversion, or rotation.
6. **All assertions in this module are `BELFEM_ASSERT` (debug only).** Wrap all failure-path tests in `#ifndef NDEBUG`.
7. **The `ten_to_mat` Blaze path uses `operator()(i,j)` while Armadillo uses `data()`.** The round-trip test covers both paths implicitly via the backend-agnostic API.

---

## 9. Codex Audit Checklist

When reviewing Claude Code's test implementation, verify:

- [ ] Reference loop helpers are present and use the correct flat-storage index formula
- [ ] At least one "random-ish" (deterministic, non-trivial) tensor is tested against each reference loop
- [ ] Mat↔ten round-trip test uses a non-trivial asymmetric-looking 6×6 matrix
- [ ] Isotropic rotation invariance test uses a non-trivial rotation matrix (not just axis-aligned)
- [ ] Identity contraction test uses a non-symmetric matrix to verify the symmetric identity formula
- [ ] `EXPECT_NEAR` used for all arithmetic results, `operator==` only for fill/copy
- [ ] Debug tests wrapped in `#ifndef NDEBUG`
- [ ] BELFEM naming conventions (`t` prefix for locals)
- [ ] High-level API tests verify agreement with low-level `tensor::` namespace functions
- [ ] `invert_symmetric` round-trip checks `C_inv % C ≈ I_sym`, not just element comparison
