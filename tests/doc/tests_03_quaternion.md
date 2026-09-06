# Quaternion Module Test Suite

**Date:** 2026-03-23

---

## Test Count by Suite

| Suite | Tests | Coverage |
|-------|-------|----------|
| `Quaternion` | 89 | Construction (owned/borrowed), copy/move semantics, accessors, iterators, compound/binary operators, equality, norm, conjugate, inverse, normalize, dot/cross, rotation (member + free func), quat-matrix conversion (all 4 Shepperd branches), slerp |
| `QuaternionDebug` | 18 | Vector/axis wrong length, zero axis, init-list wrong size, /=0, inv(0), normalize(0), rotate preconditions, free-func preconditions, matrix size/det, slerp unit checks |
| **Total** | **107** | |

---

## Key Design Decisions

### Rotation Equivalence Helper
`q` and `-q` represent the same rotation. Component-wise `operator==` is NOT sufficient. The `same_rotation()` helper checks `|dot(q1, q2)| ≈ 1.0`.

### Dual Ownership Model
- **Owned mode:** `Quaternion()` mallocs its own 4-element buffer. Destructor frees it.
- **Borrowed mode:** `Quaternion(T* aData)` uses external buffer. Destructor does NOT free.
- Move constructor steals pointer from owned, copies pointer from borrowed (source NOT nulled).
- Move assignment always copies values — does NOT transfer ownership (`mOwnData` is `const bool`).
- `conj()`, `inv()`, binary operators always return owning quaternions.

### Tests Added After Review (Codex findings)
- `MoveConstructorOwned` — verifies pointer transfer via saved pointer
- `MoveConstructorBorrowed` — asserts source pointer NOT nulled
- `InvReturnsOwningQuaternion`, `BinaryOperatorReturnsOwning`
- `FromVectorWithExternalBuffer`, `FromAxisAngleWithExternalBuffer`, `CopyConstructorWithExternalBuffer`
- `Rotate120AboutDiagonal`, `OptimizedPreservesNorm`, `ResultHasDetPlusOne`
- `NinetyDegAboutZMatrix`, `IdentityMatrixGivesIdentityQuat`, `FromRotationMatrixResultIsUnit`
- `MatrixToQuatToMatrixRoundTrip`, `AllThreeRotationPathsAgree`

---

## @warning API Gotchas

- `BELFEM_EPSILON` is in the `belfem` namespace — use `belfem::BELFEM_EPSILON`, not bare `BELFEM_EPSILON`.
- Hamilton product: `mData[0]` = scalar (w), `mData[1..3]` = vector (x, y, z). Convention: `i*j = k` (right-hand rule).
- Slerp takes shortest path (negates q2 when `dot < 0`). Falls back to nlerp when nearly parallel. Does NOT bounds-check `t`.
- Shepperd's method (`quaternion_from_rotation_matrix`) has 4 branches — test exercises all via near-180-degree rotations about each axis.

