# BELFEM Quaternion Tests — Detailed Plan

**Date:** 2026-03-22
**Purpose:** Method-level test matrix for the quaternion module
**Depends on:** `tests_0_strategy.md` (conventions), `tests_2_linalg.md` (Vector/Matrix assumed tested)
**Confidence:** High on API surface and mathematical identities. High on ownership semantics (fully reviewed).

---

## Module Overview

The quaternion module consists of one class and three free-function headers:

| File | Content |
|---|---|
| `cl_Quaternion.hpp` | `Quaternion<T>` class + non-member operators, `dot`, `cross`, `slerp` |
| `fn_quaternion_to_rotation_matrix.hpp` | Unit quaternion → 3×3 rotation matrix |
| `fn_quaternion_from_rotation_matrix.hpp` | 3×3 rotation matrix → unit quaternion (Shepperd's method) |
| `fn_quaternion_rotate_vector.hpp` | Optimized in-place vector rotation |

**Key design feature:** Quaternion supports dual ownership. An optional `T * aData` parameter on most constructors lets the quaternion either own its buffer (`malloc`/`free`) or borrow an external one. The `const bool mOwnData` flag is set at construction and cannot change. This affects move semantics, destruction, and the return types of `conj()`/`inv()`/operators (which always return owning quaternions).

---

## Floating-Point Comparison Helpers

```cpp
namespace
{
    const belfem::real tEps = 100.0 * BELFEM_EPSILON ;  // for unit-norm checks, identity comparisons
    const belfem::real tTol = 1e-12 ;                    // for rotation/conversion round-trips
}
```

**Rotation equivalence helper** — since `q` and `-q` represent the same rotation, component-wise `operator==` is not sufficient for rotation comparisons:

```cpp
// Returns true if q1 and q2 represent the same rotation (i.e., q1 ≈ ±q2)
bool same_rotation( const Quaternion<real> & q1, const Quaternion<real> & q2 )
{
    real tDotAbs = std::abs( dot( q1, q2 ) );
    return std::abs( tDotAbs - 1.0 ) < tEps ;
}
```

Use `same_rotation()` for all rotation-equivalence assertions. Use `operator==` only when testing the operator itself or when sign is known.

---

## Test File Structure

```
tests/math/
├── test_Quaternion.cpp          # All quaternion tests in one file (module is compact)
```

Given the module's size (~780 lines total), a single test file with well-named test groups is cleaner than splitting across 4–5 files.

---

## 1. Construction & Ownership

### 1.1 Owned-Mode Construction `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `DefaultConstructorZeroes` | `Quaternion<real>()` → `a()==0, b()==0, c()==0, d()==0`, `data() != nullptr` |
| `FourArgConstructor` | `Quaternion<real>(1,2,3,4)` → components match |
| `FromVector3PureQuaternion` | `Quaternion(Vector{x,y,z})` → `a()==0, b()==x, c()==y, d()==z` |
| `FromAxisAngleProducesUnitQuat` | `Quaternion(axis, angle)` → `norm() ≈ 1.0` |
| `FromAxisAngleNormalizesAxis` | Non-unit axis `{2,0,0}` with angle π/2 → same result as unit axis `{1,0,0}` |
| `IdentityStaticMethod` | `Quaternion<real>::identity()` → `(1,0,0,0)` |
| `InitializerListAssignment` | `tQ = {1,2,3,4}` → components match |

### 1.2 Borrowed-Mode Construction `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `BorrowedBufferUsesExternalPointer` | `real buf[4]; Quaternion<real> q(buf)` → `q.data() == buf` |
| `BorrowedBufferWriteThrough` | Modify `q.a()` → `buf[0]` changes |
| `BorrowedBufferDestructorDoesNotFree` | Construct borrowed quaternion, destroy it, verify buffer is still valid (read back values) |
| `FourArgWithExternalBuffer` | `Quaternion(1,2,3,4, buf)` → `buf[0]==1, buf[1]==2, ...` |
| `FromVectorWithExternalBuffer` | `Quaternion(vec3, buf)` → `buf[0]==0, buf[1..3]` match vector |
| `FromAxisAngleWithExternalBuffer` | `Quaternion(axis, angle, buf)` → result written into buf |
| `CopyConstructorWithExternalBuffer` | `Quaternion(src, buf)` → copies values into buf, `data() == buf` |

### 1.3 Copy & Move Semantics `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CopyConstructorDeepCopies` | Copy an owned quaternion, modify copy, original unchanged |
| `CopyAssignmentCopiesValues` | `q2 = q1` → values match, `q2.data() != q1.data()` (separate buffers) |
| `SelfCopyAssignment` | `q = q` does not corrupt |
| `MoveConstructorOwned` | Move from owned: target gets source's pointer, source `data() == nullptr` |
| `MoveConstructorBorrowed` | Move from borrowed: target gets same external pointer, source pointer **not** nulled |
| `MoveAssignmentCopiesValues` | Move assignment copies component values (does not transfer ownership — `mOwnData` is const) |
| `ConjReturnsOwningQuaternion` | `conj()` on a borrowed quaternion returns an owning quaternion (`data()` differs from original buffer) |
| `InvReturnsOwningQuaternion` | Same for `inv()` |
| `BinaryOperatorReturnsOwning` | `q1 + q2` where both are borrowed → result owns its buffer |

### 1.4 Construction `[debug]`

| Test Name | What It Verifies |
|---|---|
| `FromVectorWrongLengthThrows` | `Quaternion(Vector(2))` triggers assertion |
| `FromAxisAngleWrongLengthThrows` | `Quaternion(Vector(5), angle)` triggers assertion |
| `FromAxisAngleZeroAxisThrows` | `Quaternion(Vector{0,0,0}, angle)` triggers assertion |
| `InitializerListWrongSizeThrows` | `tQ = {1,2,3}` triggers assertion |

### 1.5 Memory `[valgrind]`

All ownership tests are `[valgrind]` candidates. The borrowed-buffer destructor path is the primary Valgrind target — must not double-free or leak.

---

## 2. Accessors & Iterators

### 2.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `AccessorsMatchData` | `a() == data()[0]`, `b() == data()[1]`, `c() == data()[2]`, `d() == data()[3]` |
| `ConstAccessors` | `const Quaternion& ref` → `ref.a()`, `ref.data()` compile and return correct values |
| `IteratorRange` | `end() - begin() == 4` |
| `RangeBasedForLoop` | Range-based for visits 4 elements matching `{a(), b(), c(), d()}` |

---

## 3. Algebraic Operations

### 3.1 Compound Operators `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `PlusEquals` | `(1,2,3,4) += (5,6,7,8)` → `(6,8,10,12)` |
| `MinusEquals` | `(6,8,10,12) -= (5,6,7,8)` → `(1,2,3,4)` |
| `TimesEqualsScalar` | `(1,2,3,4) *= 2.0` → `(2,4,6,8)` |
| `DivideEqualsScalar` | `(2,4,6,8) /= 2.0` → `(1,2,3,4)` |
| `HamiltonProductInPlace` | `q1 *= q2` matches Hamilton product formula (verified against hand computation) |
| `HamiltonProductNotCommutative` | `q1 * q2 != q2 * q1` for general q1, q2 |
| `HamiltonProductAssociative` | `(q1 * q2) * q3 ≈ q1 * (q2 * q3)` |

### 3.2 Binary Operators `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `Addition` | `q1 + q2` matches componentwise sum |
| `Subtraction` | `q1 - q2` matches componentwise difference |
| `ScalarMultiplyRight` | `q * 3.0` scales all components |
| `ScalarMultiplyLeft` | `3.0 * q` same result |
| `ScalarDivision` | `q / 2.0` divides all components |
| `QuaternionMultiplication` | `q1 * q2` matches Hamilton product |

### 3.3 Division `[debug]` — NOTE: `/=` uses `BELFEM_ERROR` (always active)

| Test Name | What It Verifies |
|---|---|
| `DivideByZeroThrows` | `q /= 0.0` triggers `BELFEM_ERROR` (test in both debug and release) |

### 3.4 Equality `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `EqualitySameQuaternion` | `q == q` → true |
| `EqualityWithinEpsilon` | Two quaternions differing by `0.5 * BELFEM_EPSILON` per component → true |
| `InequalityBeyondEpsilon` | Differ by `2 * BELFEM_EPSILON` in one component → false |
| `InequalityOperator` | `q1 != q2` iff `!(q1 == q2)` |
| `EqualityIsNotRotationEquivalence` | `q` and `-q` → `q != -q` (documents intentional behavior) |

---

## 4. Special Functions

### 4.1 Norm `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `NormIdentity` | `identity().norm() == 1.0` |
| `NormZero` | `Quaternion(0,0,0,0).norm() == 0.0` |
| `NormScaling` | `(q * s).norm() ≈ q.norm() * |s|` |
| `NormMultiplicative` | `(q1 * q2).norm() ≈ q1.norm() * q2.norm()` |

### 4.2 Conjugate `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `ConjNegatesImaginary` | `conj(a,b,c,d)` → `(a,-b,-c,-d)` |
| `DoubleConjIsIdentity` | `q.conj().conj() ≈ q` |
| `QTimesConjIsScalar` | `q * q.conj()` → `(norm²,0,0,0)` |

### 4.3 Inverse `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `QTimesInvIsIdentity` | `q * q.inv() ≈ identity()` for non-zero q |
| `InvOfUnitIsConj` | For unit quaternion: `inv() ≈ conj()` |

### 4.4 Inverse and Normalize `[debug]` — NOTE: both use `BELFEM_ERROR` (always active)

| Test Name | What It Verifies |
|---|---|
| `InvZeroThrows` | `Quaternion(0,0,0,0).inv()` triggers `BELFEM_ERROR` |
| `NormalizeZeroThrows` | `Quaternion(0,0,0,0).normalize()` triggers `BELFEM_ERROR` |

### 4.5 Normalize `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `NormalizeProducesUnitNorm` | `q.normalize()` → `q.norm() ≈ 1.0` |
| `NormalizePreservesDirection` | Normalized quaternion is proportional to original |
| `NormalizeReturnsSelf` | `&(q.normalize()) == &q` (returns `*this`) |

### 4.6 Dot and Cross `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `DotProductBasic` | `dot({1,0,0,0}, {0,1,0,0})` → `0.0` (orthogonal) |
| `DotProductSelf` | `dot(q, q) ≈ norm(q)²` |
| `CrossProductIsVectorPartOnly` | `cross(q1, q2).a() == 0` (always pure quaternion) |
| `CrossProductAnticommutative` | `cross(q1, q2) ≈ -cross(q2, q1)` (up to sign of imaginary part) |
| `CrossProductMatchesVectorCross` | `cross(q1, q2)` imaginary part matches 3D cross product of `(b,c,d)` parts |

---

## 5. Rotation Semantics

### 5.1 Member `rotate()` `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `IdentityRotationIsNoOp` | `identity().rotate(v) ≈ v` |
| `Rotate90AboutZ` | `(cos45, 0, 0, sin45).rotate({1,0,0})` → `{0,1,0}` |
| `Rotate180AboutX` | `(0,1,0,0).rotate({0,1,0})` → `{0,-1,0}` |
| `Rotate120AboutDiagonal` | 120° about `{1,1,1}/√3` maps x→y→z→x cyclically |
| `RotationPreservesNorm` | `norm(q.rotate(v)) ≈ norm(v)` for arbitrary v and unit q |
| `RotateAndInverseRotateRecovers` | `q.conj().rotate(q.rotate(v)) ≈ v` |
| `CompositionViaMultiplication` | `(q2*q1).rotate(v) ≈ q2.rotate(q1.rotate(v))` |
| `SignEquivalence` | `q.rotate(v) ≈ (-q).rotate(v)` (q and -q are the same rotation) |

### 5.2 Free Function `quaternion_rotate_vector()` `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `OptimizedMatchesMemberRotate` | `quaternion_rotate_vector(q, v, out)` matches `q.rotate(v)` for multiple q and v |
| `OptimizedPreservesNorm` | `norm(result) ≈ norm(input)` |

### 5.3 Rotation `[debug]`

| Test Name | What It Verifies |
|---|---|
| `RotateNonUnitThrows` | `Quaternion(2,0,0,0).rotate(v)` triggers assertion |
| `RotateWrongVectorLengthThrows` | `q.rotate(Vector(5))` triggers assertion |
| `OptimizedNonUnitThrows` | `quaternion_rotate_vector(non_unit, v, out)` triggers assertion |
| `OptimizedInputLengthThrows` | `quaternion_rotate_vector(q, Vector(2), out)` triggers assertion |
| `OptimizedOutputLengthThrows` | `quaternion_rotate_vector(q, v, Vector(5))` triggers assertion |

---

## 6. Quaternion ↔ Rotation Matrix Conversion

### 6.1 `quaternion_to_rotation_matrix` `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `IdentityQuatGivesIdentityMatrix` | `identity()` → `R ≈ I(3×3)` |
| `ResultIsOrthogonal` | `R * Rᵀ ≈ I` |
| `ResultHasDetPlusOne` | `det(R) ≈ +1` |
| `90DegAboutZMatrix` | Known quaternion → known rotation matrix entries |
| `MatrixRotationMatchesQuaternionRotation` | `R * v ≈ q.rotate(v)` for several test vectors |

### 6.2 `quaternion_from_rotation_matrix` (Shepperd's method) `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `IdentityMatrixGivesIdentityQuat` | `R = I` → `q ≈ identity()` |
| `ResultIsUnit` | `q.norm() ≈ 1.0` |
| `BranchPositiveTrace` | Use rotation with positive trace (e.g., small angle) — exercises first branch |
| `BranchXLargest` | Rotation around X axis (r00 dominant) — exercises second branch |
| `BranchYLargest` | Rotation around Y axis (r11 dominant) — exercises third branch |
| `BranchZLargest` | Rotation around Z axis (r22 dominant) — exercises fourth branch |

### 6.3 Round-Trip Consistency `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `QuatToMatrixToQuatRoundTrip` | `q → R → q2`: `same_rotation(q, q2)` for several angles/axes |
| `MatrixToQuatToMatrixRoundTrip` | `R → q → R2`: `R2 ≈ R` entrywise |
| `AllThreeRotationPathsAgree` | For a given q and v: `q.rotate(v) ≈ quaternion_rotate_vector(q,v,out) ≈ R*v` |

### 6.4 Conversions `[debug]`

| Test Name | What It Verifies |
|---|---|
| `ToMatrixNon3x3Throws` | Output matrix not 3×3 → assertion |
| `ToMatrixNonUnitQuatThrows` | Non-unit quaternion → assertion |
| `FromMatrixNon3x3Throws` | Input matrix not 3×3 → assertion |
| `FromMatrixBadDetThrows` | Matrix with determinant ≠ +1 → assertion |

---

## 7. SLERP (Spherical Linear Interpolation)

### 7.1 Endpoint and Midpoint `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `SlerpAtZeroReturnsFirst` | `slerp(q1, q2, 0.0) ≈ q1` (via `same_rotation`) |
| `SlerpAtOneReturnsSecond` | `slerp(q1, q2, 1.0) ≈ q2` (via `same_rotation`) |
| `SlerpMidpoint90Deg` | Between identity and 90° about Z at t=0.5 → 45° about Z |
| `SlerpSameQuaternion` | `slerp(q, q, t) ≈ q` for any t |
| `SlerpResultIsUnit` | Output has `norm() ≈ 1.0` |

### 7.2 Code Path Coverage `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `SlerpShortestPath` | When `dot(q1,q2) < 0`, slerp negates q2 and still produces correct rotation |
| `SlerpNearlyParallelFallback` | When q1 ≈ q2 (cosθ > 1-ε), linear interpolation fallback produces valid unit quaternion |
| `SlerpStandardFormula` | When angle is moderate (e.g., 60°), standard sin-based formula is used |

### 7.3 SLERP `[debug]`

| Test Name | What It Verifies |
|---|---|
| `SlerpNonUnitFirstThrows` | Non-unit q1 triggers assertion |
| `SlerpNonUnitSecondThrows` | Non-unit q2 triggers assertion |

### 7.4 SLERP Notes

`slerp()` does **not** assert `t ∈ [0, 1]`. Extrapolation (t < 0 or t > 1) is valid mathematically and should not be treated as a failure. Do not add bounds-checking tests for t.

---

## 8. Error Mechanism Summary

This module uses both error macros. The test plan must account for which is which:

| Method | Macro | Active In |
|---|---|---|
| Vector constructor length check | `BELFEM_ASSERT` | Debug only |
| Axis-angle axis length check | `BELFEM_ASSERT` | Debug only |
| Axis-angle zero-axis check | `BELFEM_ASSERT` | Debug only |
| `rotate()` unit-quaternion check | `BELFEM_ASSERT` | Debug only |
| `rotate()` vector length check | `BELFEM_ASSERT` | Debug only |
| `quaternion_rotate_vector()` preconditions | `BELFEM_ASSERT` | Debug only |
| `quaternion_to_rotation_matrix()` preconditions | `BELFEM_ASSERT` | Debug only |
| `quaternion_from_rotation_matrix()` preconditions | `BELFEM_ASSERT` | Debug only |
| `operator/=` divide by zero | `BELFEM_ERROR` | Always |
| `inv()` zero quaternion | `BELFEM_ERROR` | Always |
| `normalize()` zero quaternion | `BELFEM_ERROR` | Always |
| `slerp()` unit-quaternion checks | `BELFEM_ASSERT` | Debug only |
| Initializer list length | `BELFEM_ASSERT` | Debug only |

Tests for `BELFEM_ASSERT` paths must be wrapped in `#ifndef NDEBUG`. Tests for `BELFEM_ERROR` paths run in both debug and release.

---

## 9. Implementation Notes for Claude Code

1. **Use `same_rotation()` helper** for all tests that compare quaternions as rotations. Use `operator==` only when testing the operator itself or when the exact sign is known.
2. **Shepperd branch coverage:** To exercise all four branches of `quaternion_from_rotation_matrix`, use rotations that make each branch dominant. E.g., 170° about X (r00 largest, trace negative), 170° about Y (r11 largest), 170° about Z (r22 largest), 30° about any axis (trace positive).
3. **Ownership tests** are `[valgrind]` candidates. Run under Valgrind to catch double-free, use-after-free, or leak in the borrowed-buffer path.
4. **Do not test with `float` or `complex`.** BELFEM uses `real` = `double`. The template parameter exists for generality but is not exercised with other types in the framework.
5. **Combine into one file.** The module is compact enough that `test_Quaternion.cpp` with well-named `TEST` groups (`QuaternionOwnership`, `QuaternionAlgebra`, `QuaternionRotation`, `QuaternionConversion`, `QuaternionSlerp`, `QuaternionDebug`) is cleaner than 5 separate files.
6. **Expression results:** Binary operators (`+`, `-`, `*`, `/`) return owning `Quaternion<T>` objects. No special handling needed (unlike linalg expression templates).

---

## 10. Codex Audit Checklist

When reviewing Claude Code's test implementation, verify:

- [ ] Every row in the test matrices above has a corresponding `TEST`
- [ ] Debug tests use correct guard: `#ifndef NDEBUG` for `BELFEM_ASSERT`, no guard for `BELFEM_ERROR`
- [ ] BELFEM naming conventions (`t` prefix for locals, BELFEM types)
- [ ] `EXPECT_NEAR` used for all floating-point comparisons
- [ ] `same_rotation()` helper used for rotation-equivalence checks, not `operator==`
- [ ] All four Shepperd branches exercised with appropriate rotation matrices
- [ ] Borrowed-buffer tests verify pointer identity (`data() == buf`)
- [ ] Owned-mode move test verifies source `data() == nullptr`
- [ ] Borrowed-mode move test verifies source pointer **not** nulled
- [ ] `slerp` tests cover all three code paths (standard, shortest-path, near-parallel fallback)
- [ ] No bounds-checking test on slerp's `t` parameter
- [ ] Cross-consistency test comparing `rotate()`, `quaternion_rotate_vector()`, and `R*v`
