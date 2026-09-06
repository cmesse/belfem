# Quaternion Usage Guide {#math_quaternion_quaternion_usage_guide}

**Date:** 2026-06-14
**Purpose:** Comprehensive guide to BELFEM's quaternion math helper
**Module:** src/math/quaternion

---

## Overview

`Quaternion<T>` is a small, header-only value type for 3D rotations. It is used
to compose rotations, interpolate orientations (slerp), rotate vectors, and
convert to/from 3x3 rotation matrices. It is a *plain value type*: four inline
components, trivially copyable, no heap allocation, safe to `memcpy` /
MPI-transfer (`cl_Quaternion.hpp:34-39`, `cl_Quaternion.hpp:87-89`).

The module is header-only; its `CMakeLists.txt` declares no library target
(`src/math/quaternion/CMakeLists.txt:1`).

---

## Convention (read this first)

These three choices are the classic silent-bug sources. BELFEM's are:

| Aspect | BELFEM choice | Evidence |
|--------|---------------|----------|
| **Component storage** | **Scalar-first** `q = (w, x, y, z)` | `mData[4]` with accessors `a()=w`, `b()=x`, `c()=y`, `d()=z` (`cl_Quaternion.hpp:34-39`, `:110-148`) |
| **Multiplication** | **Hamilton** (not JPL) | product carries `+v1 x v2` cross term (`cl_Quaternion.hpp:226-243`); tests pin `i*j=k`, `j*k=i`, `k*i=j` (`tests/math/test_Quaternion.cpp:374-413`) |
| **Vector rotation** | **Active**: `v' = q v q*` | `rotate()` embeds `v` as `(0,v)` and computes `q*v*conj(q)` (`cl_Quaternion.hpp:286-304`) |

The source does not use the words *active/passive* or *alibi/alias*; "active"
here is the behavioral fact -- rotating the vector in a fixed frame -- confirmed
by the tests pinning `+90 deg` about `+z`: `(1,0,0) -> (0,1,0)`, matching `R*v`.

**Not implemented:** Euler-angle conversion and quaternion-to-axis-angle
*extraction*. Axis-angle is only an *input* (the `(axis, angle)` constructor). A
separate DIN-9300 Euler helper exists outside this module
(`src/math/tools/fn_rotation_matrix.hpp`).

---

## Quick Reference

### Class API (`cl_Quaternion.hpp`)

| Operation | Signature | Notes |
|-----------|-----------|-------|
| Default | `Quaternion()` | all zeros |
| Components | `Quaternion(a,b,c,d)` | `(w,x,y,z)` |
| Pure quaternion | `Quaternion(const Vector<T>&)` | length-3 vector -> `(0,v)` |
| Rotation | `Quaternion(const Vector<T>& axis, T angle)` | axis normalized internally; angle in radians |
| Identity | `static Quaternion identity()` | `(1,0,0,0)` |
| Accessors | `a() b() c() d()`, `data()`, `begin()/end()` | `a` is the scalar |
| Norm | `T norm()` | Euclidean 4-norm |
| Conjugate | `Quaternion conj()` | negates vector part |
| Inverse | `Quaternion inv()` | `conj()/|q|^2`; **errors** on zero |
| Normalize | `Quaternion& normalize()` | in place; **errors** on zero |
| Rotate vector | `Vector<T> rotate(const Vector<T>&)` | requires unit quaternion |

### Operators and free functions

| Expression | Meaning |
|------------|---------|
| `q1 * q2` | Hamilton product (also `*=`) |
| `q + q2`, `q - q2`, `s*q`, `q*s`, `q/s` | componentwise / scalar |
| `q1 == q2`, `q1 != q2` | componentwise within `BELFEM_EPSILON` (does **not** treat `q` and `-q` as equal) |
| `dot(q1,q2)`, `cross(q1,q2)` | 4-vector dot; vector-part cross as pure quaternion |
| `slerp(q1,q2,t)` | shortest-arc spherical interpolation |
| `quaternion_to_rotation_matrix(q, R)` | `R` is a pre-allocated 3x3 (`fn_quaternion_to_rotation_matrix.hpp`) |
| `quaternion_from_rotation_matrix(R) -> q` | Shepperd's method (`fn_quaternion_from_rotation_matrix.hpp`) |
| `quaternion_rotate_vector(q, in, out)` | allocation-free `v + 2w(uxv) + 2 ux(uxv)`, equals `q v q*` (`fn_quaternion_rotate_vector.hpp`) |

---

## Common Patterns

### Build a rotation and rotate a vector

```cpp
#include "cl_Quaternion.hpp"

Vector< real > tAxis( 3 );
tAxis( 0 ) = 0.0; tAxis( 1 ) = 0.0; tAxis( 2 ) = 1.0;   // +z

Quaternion< real > tQ( tAxis, constant::pi * 0.5 );      // +90 deg about z

Vector< real > tV( 3 );
tV( 0 ) = 1.0; tV( 1 ) = 0.0; tV( 2 ) = 0.0;             // x-axis

Vector< real > tRotated = tQ.rotate( tV );               // -> (0, 1, 0)
```

### Allocation-free rotation in a hot loop

```cpp
#include "fn_quaternion_rotate_vector.hpp"

Vector< real > tOut( 3 );                                // preallocate once
for ( ... )
{
    quaternion_rotate_vector( tQ, tV, tOut );            // no temporaries
}
```

`rotate()` (the member) constructs intermediate quaternions; prefer
`quaternion_rotate_vector()` on critical paths.

### Compose rotations

```cpp
Quaternion< real > tQ = tQ2 * tQ1;     // apply tQ1 first, then tQ2 (Hamilton)
tQ.normalize();                        // renormalize after repeated products
```

### Interpolate orientation

```cpp
Quaternion< real > tQ = slerp( tStart, tEnd, 0.5 );      // unit inputs expected
```

### Matrix round trip

```cpp
#include "fn_quaternion_to_rotation_matrix.hpp"
#include "fn_quaternion_from_rotation_matrix.hpp"

Matrix< real > tR( 3, 3 );
quaternion_to_rotation_matrix( tQ, tR );
Quaternion< real > tBack = quaternion_from_rotation_matrix( tR );  // == +/- tQ
```

---

## Numerical Aspects

- **Unit-quaternion preconditions** are `BELFEM_ASSERT`s (compiled out in
  release): `rotate()`, `quaternion_rotate_vector()`,
  `quaternion_to_rotation_matrix()`, and `slerp()` require `|q| - 1` within
  `100 * BELFEM_EPSILON` (`cl_Quaternion.hpp:288-290`, `:402-405`;
  `fn_quaternion_rotate_vector.hpp:42-45`;
  `fn_quaternion_to_rotation_matrix.hpp:41-45`). Renormalize after long product
  chains to keep these satisfied.
- **Always-active checks** (`BELFEM_ERROR`) guard the genuinely singular cases:
  division by zero scalar, `inv()` of a zero quaternion, and `normalize()` of a
  zero quaternion (`cl_Quaternion.hpp:218`, `:264`, `:272`).
- **`quaternion_from_rotation_matrix`** uses the trace-first variant of
  Shepperd's method: when `trace > 0` it recovers `w` from the trace, otherwise
  it picks the largest of `{r00, r11, r22}` to avoid dividing by a small number,
  then normalizes the result. It asserts `det == +1` but **does not** check
  orthogonality (`fn_quaternion_from_rotation_matrix.hpp:49-98`).
- **`slerp`** negates the second quaternion on a negative dot product (shortest
  arc) and falls back to normalized linear interpolation when the inputs are
  nearly parallel, avoiding division by `sin(theta) -> 0`
  (`cl_Quaternion.hpp:407-424`).

---

## Status and Pitfalls

- **Tested, not yet consumed in production.** The only call sites in the
  checkout are `tests/math/test_Quaternion.cpp`. No `src/` module currently
  rotates geometry, frames, or material orientation with `Quaternion`. The
  helper is implemented and unit-tested but awaits a first production user.
- **The borrowed-buffer mode is gone, and the suite already says so.**
  `Quaternion` was once able to alias an external buffer; it is now a plain
  inline value type, trivially copyable. `tests/math/test_Quaternion.cpp:117-125`
  records the removal and notes that the old borrowed-mode cases targeted an API
  that no longer exists — so there is nothing stale left to work around here.
- **`q` and `-q` are the same rotation** but compare unequal under `operator==`
  (it compares 4-vectors). Use `slerp`/rotation output if you need
  rotation-equivalence rather than component equality.
- **Multiplication order is Hamilton**: `q2 * q1` applies `q1` first. Mixing in
  a JPL-convention reference will silently transpose your rotations.

---

## Files

| File | Contents |
|------|----------|
| `cl_Quaternion.hpp` | `Quaternion<T>` class, operators, `dot`, `cross`, `slerp` |
| `fn_quaternion_to_rotation_matrix.hpp` | unit quaternion -> 3x3 matrix |
| `fn_quaternion_from_rotation_matrix.hpp` | 3x3 matrix -> unit quaternion (Shepperd) |
| `fn_quaternion_rotate_vector.hpp` | allocation-free active vector rotation |

---

## See Also

- `src/math/tensor/doc/tensor_usage_guide.md` - the sibling math helper (4th-order tensors)
