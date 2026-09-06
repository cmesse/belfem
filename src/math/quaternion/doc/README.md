# Quaternion Module Documentation {#math_quaternion_index}

**Module:** src/math/quaternion
**Purpose:** Index of documentation for BELFEM's quaternion math helper

---

## Overview

The `quaternion` module provides a header-only `Quaternion<T>` value type for 3D
rotations: composition, spherical interpolation (slerp), vector rotation, and
conversion to/from 3x3 rotation matrices.

**Convention at a glance:** scalar-first `(w,x,y,z)`, Hamilton multiplication,
active `v' = q v q*` vector rotation. See the usage guide for the full rationale
and evidence.

---

## Documentation Files

- **[quaternion_usage_guide.md](quaternion_usage_guide.md)** - Comprehensive guide
  - Convention (storage / Hamilton / active rotation) with code evidence
  - Class API, operators, and free functions
  - Common patterns (build rotation, hot-loop rotation, compose, slerp, matrix round trip)
  - Numerical aspects (unit preconditions, Shepperd extraction, slerp fallback)
  - Status (tested, not yet a production user) and pitfalls

---

## Quick Reference

### Key Type

| Type | File | Purpose |
|------|------|---------|
| `Quaternion<T>` | cl_Quaternion.hpp | scalar-first unit/quaternion value type |

### Files

| File | Purpose |
|------|---------|
| cl_Quaternion.hpp | class, operators, `dot`, `cross`, `slerp` |
| fn_quaternion_to_rotation_matrix.hpp | unit quaternion -> 3x3 matrix |
| fn_quaternion_from_rotation_matrix.hpp | 3x3 matrix -> unit quaternion (Shepperd) |
| fn_quaternion_rotate_vector.hpp | allocation-free active vector rotation |

### Convention

```
q = (w, x, y, z)            scalar-first
q1 * q2                     Hamilton product (q2*q1 applies q1 first)
v' = q v q*                 active rotation of a vector in a fixed frame
```

---

## Notes

- Header-only; no library target (`CMakeLists.txt`).
- Not thread-safe by design (BELFEM-wide policy); value semantics make copies cheap and MPI-safe.
- No Euler conversion and no quaternion->axis-angle extraction are implemented.

---

## See Also

- `src/math/tensor/doc/README.md` - sibling math helper (4th-order tensors)
