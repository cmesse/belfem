# Fourth-Order Tensor Usage Guide {#math_tensor_tensor_usage_guide}

**Date:** 2026-06-14
**Purpose:** Comprehensive guide to BELFEM's fourth-order tensor math helper
**Module:** src/math/tensor

---

## Overview

`Tensor<T>` is a fourth-order tensor container built for constitutive
(elasticity) work: it stores a full 3x3x3x3 tensor, double-contracts with
second- and fourth-order tensors (`sigma = C : eps`, `C : S`), rotates under a
3x3 frame, inverts a symmetric stiffness/compliance tensor, and converts to/from
the 6x6 Voigt elasticity matrix.

The class is templated and (apart from an optional example executable) used
header-only. Storage is **manually managed** (`malloc`/`memcpy`/`free`) like the
rest of BELFEM's hot-path containers (`cl_Tensor.hpp:82-167`).

The primary production consumer is self-consistent alloy homogenization
(`src/physics/materials/cl_Material_Alloy.cpp`); a worked example lives in
`example_udfiber.cpp`.

---

## Notation and storage (read this first)

| Aspect | BELFEM choice | Evidence |
|--------|---------------|----------|
| **Internal storage** | **Full 3x3x3x3** (81 entries), not Voigt-only/Mandel | capacity = `i*j*k*l`; index `L*27 + K*9 + J*3 + I` (`cl_Tensor.hpp:66-78`, `:441`) |
| **Index order** | `i` fastest (column-major-consistent with BELFEM matrices) | `mOffsetJ=3, mOffsetK=9, mOffsetL=27` (`cl_Tensor.hpp:74-76`) |
| **Matrix boundary** | **Engineering Voigt 6x6**, order `{11, 22, 33, 23, 31, 12}` | `ten_to_mat`/`mat_to_ten` maps; `compliance_matrix` puts `1/G23,1/G31,1/G12` on rows 3,4,5 (`fn_TR_ten_to_mat.hpp:74-109`, `fn_compliance_matrix.hpp:36-50`) |
| **Shear factors** | **No Mandel `sqrt(2)`**; stiffness needs no Voigt factor | direct assignment in conversion tables; tests pin `C44=C55=C66=mu` (`tests/math/test_Tensor.cpp:691-721`) |
| **Symmetries** | Minor `C_ijkl=C_jikl=C_ijlk` imposed by the 6x6 conversion; major `C_ijkl=C_klij` holds when the 6x6 is symmetric | `mat_to_ten` duplicates minor pairs; `compliance_matrix` is symmetric (`fn_TR_mat_to_ten.hpp:120-200`) |

**Key consequence:** the 6x6 boundary is the *stiffness* convention (no factors of
2). If you have a *compliance* matrix (engineering strain, factors of 2 on the
shear), do **not** feed it through `mat_to_ten` directly. The provided paths do
the right thing: `fill_orthotropic_elasticity()` builds compliance, inverts to
stiffness, then converts; `invert_symmetric()` handles the factor-of-2 shear
bookkeeping internally (it doubles the shear-pair terms when forming its 6x6 work
matrix and uses a `0.5`-shear identity, `fn_invert_symmetric.hpp:67-96`,
`:112-159`).

`Tensor<T>` is a *general* container: it holds either third-order (`i,j,k`) or
fourth-order (`i,j,k,l`) tensors of arbitrary dimensions, with the tensor order
reported by `order()`. Construction, element access, element-wise add/subtract,
scalar arithmetic, and equality work for any shape. The constitutive helpers
(contraction, rotation, inversion, Voigt conversion, isotropic/orthotropic
fills), however, assert `is_3333()` -- they are specialized to the 3x3x3x3
elasticity case.

---

## Quick Reference

### Construction and fill (`cl_Tensor.hpp`, `fn_TR_fill.hpp`, `fn_compliance_matrix.hpp`)

| Call | Result |
|------|--------|
| `Tensor<real> A(3,3,3,3)` | uninitialized fourth-order 3x3x3x3 |
| `Tensor<real> A(3,3,3)` | uninitialized third-order 3x3x3 |
| `Tensor<real> A(3,3,3,3, 0.0)` | filled with a scalar |
| `Tensor<real> C(M6x6)` | from a 6x6 Voigt elasticity matrix |
| `A.fill(value)` | uniform fill |
| `A.fill(a, b)` | isotropic: `A_ijkl = a d_ij d_kl + b(d_ik d_jl + d_il d_jk - 2/3 d_ij d_kl)` |
| `C.fill_isotropic_elasticity(E, nu)` | computes `K=E/(3(1-2nu))`, `G=E/(2(1+nu))`, then `fill(K,G)` |
| `C.fill_orthotropic_elasticity(E1,E2,E3, nu23,nu13,nu12, G23,G31,G12)` | builds compliance, inverts, converts |

### Access

| Call | Meaning |
|------|---------|
| `A(I,J,K,L)` / `A(I,J,K)` | element access for fourth- / third-order tensors (bounds-checked in debug) |
| `A.data()` | raw pointer (bulk ops) |
| `A.order()`, `A.size_i()` ... `A.size_l()` | tensor order and per-index sizes |
| `A.capacity()`, `A.is_3333()`, `A.print()` | size / type / dump |

### Operations

| Expression / call | Meaning |
|-------------------|---------|
| `A += B`, `A -= B`, `A + B`, `A - B` | element-wise add/subtract for matching order and shape |
| `A += s`, `A *= s`, `A /= s`, ... | scalar arithmetic |
| `A.ddot(B, C)` / `ddot(A,B,C)` / `C = A % B` | fourth:fourth, `C_ijkl = A_ijmn B_mnkl` |
| `C.ddot(eps, sigma)` / `ddot(C,eps,sigma)` / `sigma = C % eps` | fourth:second, `sigma_ij = C_ijkl eps_kl` |
| `rotate(B, R, A)` | `A_mnop = B_ijkl R_im R_jn R_ko R_lp` |
| `tensor::identity(A)` | symmetric identity `I_ijkl = 0.5(d_ik d_jl + d_il d_jk)` |
| `tensor::invert_symmetric(A, work, pivot)` | in-place inverse of a symmetric tensor |
| `kelvin_christoffel(A, n, Gamma)` | acoustic tensor `Gamma_ik = A_ijkl n_j n_l` |
| `fiber_polarization(C6x6, P)` | Eshelby/Hill polarization tensor for fibers |
| `C.to_matrix(M6x6)` | tensor -> 6x6 Voigt stiffness |
| `A == B` | exact componentwise equality |

`ddot`'s second-order overloads take/return a 3x3 `Matrix<T>` (the
stress/strain), **not** a Voigt 6-vector.

---

## Common Patterns

### Isotropic stress from strain

```cpp
#include "cl_Tensor.hpp"
#include "fn_ddot.hpp"

Tensor< real > C( 3, 3, 3, 3 );
C.fill_isotropic_elasticity( 210.0e9, 0.30 );    // E [Pa], nu

Matrix< real > tEps( 3, 3 );                     // symmetric strain (true tensor strain)
// ... fill tEps ...

Matrix< real > tSigma( 3, 3 );
C.ddot( tEps, tSigma );                          // sigma_ij = C_ijkl eps_kl
// equivalently: tSigma = C % tEps;
```

Note the **true** tensor strain: for pure shear use `eps_12 = eps_21 = gamma/2`
(the tests pin `eps_12 = eps_21 = 0.5 -> sigma_12 = mu`).

### Orthotropic material and Voigt round trip

```cpp
Tensor< real > C( 3, 3, 3, 3 );
C.fill_orthotropic_elasticity( E1,E2,E3, nu23,nu13,nu12, G23,G31,G12 );

Matrix< real > tC6( 6, 6 );
C.to_matrix( tC6 );                              // engineering-Voigt stiffness
Tensor< real > C2( tC6 );                        // back to a tensor (minor-symmetric)
```

### Invert a stiffness tensor to compliance

```cpp
#include "fn_invert_symmetric.hpp"

Vector< real > tWork( 72 );                      // >= 72
Vector< int_t > tPivot( 36 );                     // >= 36

Tensor< real > S = C;                            // copy
tensor::invert_symmetric( S, tWork, tPivot );    // S = C^{-1}, in place
```

`invert_symmetric` reduces to a 6x6, doubles the shear-pair terms, LU-factorizes
(`getrf`), inverts (`getri`), multiplies by a `0.5`-shear identity (`gemm`), and
writes back with symmetric duplication. The caller supplies the work/pivot
buffers so the routine performs no allocation.

### Rotate a tensor into another frame

```cpp
#include "fn_rotate.hpp"

Matrix< real > tR( 3, 3 );                       // rotation matrix (e.g. from a quaternion)
Tensor< real > tRotated( 3, 3, 3, 3 );
rotate( C, tR, tRotated );                        // A_mnop = C_ijkl R_im R_jn R_ko R_lp
```

`rotate` and `ddot(.,Matrix,.)` dispatch to a backend kernel selected by
`BELFEM_ARMADILLO` / `BELFEM_BLAZE`; the remaining kernels (`contract44`,
`identity`, ...) are backend-agnostic hand-unrolled loops. Element-wise
addition and subtraction use shape-agnostic `std::transform` loops over the
tensor capacity, so they work for any order or dimension.

---

## Numerical Aspects

- **Inversion** goes through a 6x6 LAPACK solve (`getrf` + `getri`), with
  `info == 0` enforced by `BELFEM_ERROR` (active in release builds too). There is **no** condition-number or residual
  check; a near-singular elasticity tensor (e.g. `nu -> 0.5`, incompressible)
  will invert without warning. The factor-of-2 shear bookkeeping is what makes a
  tensor inverse (not the naive 6x6 inverse) come out correct
  (`fn_invert_symmetric.hpp:98-159`).
- **Symmetry is assumed, not enforced.** Raw `operator()`/`data()` writes accept
  any values; `invert_symmetric` checks only `is_3333()` and buffer lengths.
  Feeding a non-symmetric tensor produces a silently wrong "inverse."
- **Conversion imposes minor symmetry.** A tensor-to-matrix-to-tensor round trip
  only reproduces the original if it was minor-symmetric to begin with
  (`tests/math/test_Tensor.cpp:667-680`).
- **Performance.** The contraction and rotation kernels are fully unrolled over
  the 81 (or 9) entries, avoiding any inner loop. Element-wise add/subtract
  operations loop over the tensor capacity instead, trading hard-coded unrolling
  for support of arbitrary order and dimension. Pass preallocated work/pivot
  buffers and reuse output tensors across calls.

---

## Usage in the framework

- **Self-consistent alloy homogenization** (`cl_Material_Alloy.cpp`): allocates a
  `Cell<Tensor<real>*>` of elasticity, compliance, Eshelby, Hill-polarization,
  influence, identity and work tensors; fills isotropic component stiffnesses,
  iterates effective `K, G` via `S % M`, `P % B`,
  `tensor::invert_symmetric(A, ...)`, and `Ck % A`
  (`cl_Material_Alloy.cpp:288-373`, `:695-779`).
- **Worked example** (`example_udfiber.cpp`, built only when `USE_EXAMPLES`):
  fiber/matrix orthotropic + isotropic fills, `fiber_polarization`, `ddot`,
  `invert_symmetric`, and conversion back to a 6x6 matrix.

**Seam to materials/elastic work (not documented here):** the constitutive
tensor `C` that consumes `E, G, nu` is the boundary between this helper and the
material-property models. This guide stops at the tensor API; the property
source models are out of scope.

---

## Files

| File | Contents |
|------|----------|
| `cl_Tensor.hpp` | `Tensor<T>` class, operators, member `ddot`, `to_matrix` |
| `fn_TR_mat_to_ten.hpp` / `fn_TR_ten_to_mat.hpp` | 6x6 Voigt <-> full tensor conversion |
| `fn_TR_fill.hpp` | isotropic fill formula |
| `fn_compliance_matrix.hpp` | orthotropic 6x6 compliance |
| `fn_identity.hpp` | symmetric identity tensor |
| `fn_invert_symmetric.hpp` | in-place symmetric inverse via 6x6 LAPACK |
| `fn_ddot.hpp` | free double-contraction wrappers |
| `fn_TR_contract44.hpp` | fourth:fourth kernel (backend-agnostic) |
| `fn_TR_contract42.hpp` + `armadillo/`, `blaze/` | fourth:second kernel (backend-selected) |
| `fn_rotate.hpp`, `fn_TR_rotate42.hpp` + `armadillo/`, `blaze/` | tensor rotation |
| `fn_kelvin_christoffel.hpp` (+ backend) | acoustic / Christoffel tensor |
| `fn_fiber_polarization.hpp` | Hill polarization tensor for fibers |
| `fn_TR_equal_equal.hpp` | componentwise equality over the tensor capacity |
| `example_udfiber.cpp` | optional example executable (`USE_EXAMPLES`) |

---

## References

- Fiber polarization tensor: doi:10.1016/S0020-7683(02)00369-4 (`fn_fiber_polarization.hpp:19-24`).
- General Voigt/elasticity background: see `literature/books/index.md` (Hughes, Bathe).

---

## See Also

- `src/math/quaternion/doc/quaternion_usage_guide.md` - sibling math helper (rotations; produces the `R` for `rotate`)
