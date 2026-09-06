# Fourth-Order Tensor Module Documentation {#math_tensor_index}

**Module:** src/math/tensor
**Purpose:** Index of documentation for BELFEM's fourth-order tensor math helper

---

## Overview

The `tensor` module provides `Tensor<T>`, a third-/fourth-order tensor container
whose constitutive helpers target the 3x3x3x3 elasticity case: double contraction (`sigma = C : eps`,
`C : S`), frame rotation, symmetric inversion, and conversion to/from the 6x6
Voigt elasticity matrix.

**Notation at a glance:** full 81-component internal storage; the matrix
boundary is **engineering Voigt** (order `{11,22,33,23,31,12}`, no Mandel
`sqrt(2)` factors); minor symmetry is imposed by the 6x6 conversion and major
symmetry follows when the 6x6 is symmetric. See the usage guide for evidence and
the shear-factor handling.

---

## Documentation Files

- **[tensor_usage_guide.md](tensor_usage_guide.md)** - Comprehensive guide
  - Notation and storage (full tensor, Voigt boundary, shear factors, symmetries)
  - Construction / fill / access / operations quick reference
  - Common patterns (stress from strain, orthotropic, inversion, rotation)
  - Numerical aspects (LAPACK 6x6 inverse, no conditioning check, assumed symmetry)
  - Usage (alloy homogenization, udfiber example) and the materials/elastic seam

---

## Quick Reference

### Key Type

| Type | File | Purpose |
|------|------|---------|
| `Tensor<T>` | cl_Tensor.hpp | third-/fourth-order tensor container; constitutive helpers target 3x3x3x3 |

### Core Operations

| Operation | Entry point | Meaning |
|-----------|-------------|---------|
| Fill isotropic | `fill_isotropic_elasticity(E,nu)` | from `E`, `nu` |
| Fill orthotropic | `fill_orthotropic_elasticity(...)` | from 9 engineering constants |
| Fourth:second | `C.ddot(eps,sigma)` / `C % eps` | `sigma_ij = C_ijkl eps_kl` |
| Fourth:fourth | `A.ddot(B,C)` / `A % B` | `C_ijkl = A_ijmn B_mnkl` |
| Rotate | `rotate(B,R,A)` | `A_mnop = B_ijkl R_im R_jn R_ko R_lp` |
| Identity | `tensor::identity(A)` | `I_ijkl = 0.5(d_ik d_jl + d_il d_jk)` |
| Invert | `tensor::invert_symmetric(A,work,pivot)` | in-place symmetric inverse |
| To matrix | `C.to_matrix(M6x6)` | Voigt stiffness |

### Backend selection

`contract42` and `rotate42` are selected by `BELFEM_ARMADILLO` / `BELFEM_BLAZE`;
`contract44`, `identity`, etc. are backend-agnostic unrolled kernels. Element-wise
add/subtract use shape-agnostic `std::transform` loops over the tensor capacity.

---

## Notes

- Header-only apart from the optional `example_udfiber` executable (`USE_EXAMPLES`).
- Manual memory management (`malloc`/`memcpy`/`free`); pass preallocated work/pivot buffers to `invert_symmetric`.
- Symmetry and conditioning are the caller's responsibility (not validated).
- Primary production consumer: `src/physics/materials/cl_Material_Alloy.cpp` (self-consistent homogenization).

---

## See Also

- `src/math/quaternion/doc/README.md` - sibling math helper (rotations)
