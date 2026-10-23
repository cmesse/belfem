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
| **Shear factors** | **No Mandel `sqrt(2)`**; stiffness needs no Voigt factor | direct assignment in conversion tables; tests pin `C44=C55=C66=mu` (`tests/math/test_Tensor.cpp`, `IsotropicTensorToVoigtMatrix`) |
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

"General" refers to the tensor's order and shape, not its element type. `T` must
be a trivially copyable numeric type (every user instantiates `Tensor<real>`).
Storage is obtained with `malloc`, populated with `std::fill_n`/`memcpy`, and
released with `free`, so no constructor or destructor of `T` ever runs
(`cl_Tensor.hpp`, constructors and destructor).

---

## Quick Reference

### Construction and fill (`cl_Tensor.hpp`, `fn_TR_fill.hpp`, `fn_compliance_matrix.hpp`)

| Call | Result |
|------|--------|
| `Tensor<real> A(3,3,3,3)` | uninitialized fourth-order 3x3x3x3 |
| `Tensor<real> A(3,3,3)` | uninitialized third-order 3x3x3 |
| `Tensor<real> A(3,3,3,3, 0.0)` | filled with a scalar |
| `Tensor<real> C(M6x6)` | from a 6x6 Voigt elasticity matrix (`Matrix<real>` only; asserts 6x6) |
| `C = M6x6` | assign from a 6x6 Voigt matrix; asserts only `is_3333()`, not the 6x6 size |
| `A.fill(value)` | uniform fill |
| `A.fill(a, b)` | isotropic: `A_ijkl = a d_ij d_kl + b(d_ik d_jl + d_il d_jk - 2/3 d_ij d_kl)` |
| `C.fill_isotropic_elasticity(E, nu)` | computes `K=E/(3(1-2nu))`, `G=E/(2(1+nu))`, then `fill(K,G)` |
| `C.fill_orthotropic_elasticity(E1,E2,E3, nu23,nu13,nu12, G23,G31,G12)` | builds compliance, inverts, converts; allocates two 6x6 matrices per call. See the Poisson-ratio convention below |

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
| `kelvin_christoffel(A, n, Gamma)` | acoustic tensor `Gamma_ik = A_ijkl n_j n_l`; `n` is a length-3 `Vector<T>`, `Gamma` a pre-sized 3x3 `Matrix<T>` |
| `fiber_polarization(C6x6, P)` | Hill polarization tensor `P` for a unidirectional fiber in a transversely isotropic host (the Eshelby tensor follows as `S = P : C`) |
| `C.to_matrix(M6x6)` | tensor -> 6x6 Voigt stiffness (`Matrix<real>` only) |
| `A == B` | exact componentwise equality |

The second-order `ddot` overloads take and return a 3x3 `Matrix<T>` for stress
or strain, **not** a Voigt 6-vector.

**The output must not alias an input.** The `ddot` and `rotate` forms write
output entries while still reading the inputs, and none checks for aliasing.
Therefore `ddot(A, B, A)`, `ddot(A, B, B)` and `rotate(B, R, B)` silently return
wrong values (`fn_TR_contract44.hpp` and the `armadillo/`, `blaze/` kernels of
`contract42` and `rotate42`). Use a separate output tensor or matrix.

**`%` allocates.** Each call to `A % B` constructs a new `Tensor<T>(3,3,3,3)`,
and each call to `C % eps` a new 3x3 `Matrix<T>` (`cl_Tensor.hpp`, `operator%`).
The results equal those of `ddot`, which writes into a caller-supplied output
and does not allocate. Use `%` in setup code and `ddot` on hot paths.

`fiber_polarization` takes the fiber axis to be axis 1 (Voigt row 0). It reads
only `C(1,1)`, `C(1,2)`, `C(2,2)` and the shear diagonals `C(3,3)`, `C(4,4)`,
`C(5,5)` (zero-based). As a result, a general anisotropic `C` is silently
treated as transversely isotropic about axis 1 (`fn_fiber_polarization.hpp`).

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
// same result, but allocates: tSigma = C % tEps;
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

**Poisson-ratio convention.** In textbook notation, `nu_ij = -eps_j / eps_i`
under uniaxial stress in direction `i`. This gives
`eps_j = -nu_ij sigma_i / E_i`, and therefore `S_ji = -nu_ij / E_i`. Because
`S` is symmetric, `nu_ij / E_i = nu_ji / E_j`; either form gives the same
off-diagonal entry. `compliance_matrix` (called by
`fill_orthotropic_elasticity`) sets these one-based entries:

- `S_21 = -nu12 / E2`, so the argument `nu12` is the textbook `nu_21`;
- `S_31 = -nu13 / E3`, so the argument `nu13` is the textbook `nu_31`;
- `S_32 = -nu23 / E2`, so the argument `nu23` is the textbook `nu_23`
  (`= nu_32 E2 / E3`).

Thus the `nu12` and `nu13` arguments expect the ratio whose load is in the
second index, whereas `nu23` expects the ratio whose load is in the first. For a
fiber with its axis along 1, `nu12` is the minor ratio. To use a datasheet major
ratio `nu_12`, pass `nu_12 * E2 / E1`. `example_udfiber.cpp` passes the minor
ratio `nu12 = 0.013` with `E1 = 235`, `E2 = 20`, and reports the major ratio
`-S(1,0) * E1` of the homogenized result (`fn_compliance_matrix.hpp`,
`example_udfiber.cpp`).

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

`rotate`, `ddot(.,Matrix,.)` and `kelvin_christoffel` use backend kernels
selected by `BELFEM_ARMADILLO` / `BELFEM_BLAZE`. The remaining kernels
(`contract44`, `identity`, ...) are backend-agnostic, hand-unrolled loops.

Element-wise addition and subtraction use shape-agnostic `std::transform` loops
over the tensor capacity, so they work for any order or dimension.

---

## Numerical Aspects

- **Inversion** goes through a 6x6 LAPACK solve (`getrf` + `getri`), with
  `info == 0` enforced by `BELFEM_ERROR` (active in release builds too). There is **no** condition-number or residual
  check; a near-singular elasticity tensor (e.g. `nu -> 0.5`, incompressible)
  will invert without warning. The factor-of-2 shear bookkeeping is what makes a
  tensor inverse (not the naive 6x6 inverse) come out correct
  (`fn_invert_symmetric.hpp:98-159`).
- **Symmetry is assumed, not enforced.** Raw `operator()`/`data()` writes accept
  any values. `invert_symmetric` checks only `is_3333()` and the buffer lengths;
  all three checks are `BELFEM_ASSERT` and are compiled out in release builds.
  In release builds, a work vector shorter than 72 or a pivot vector shorter
  than 36 is overrun without a message. A non-symmetric tensor produces a
  silently wrong "inverse."
- **Conversion imposes minor symmetry.** A tensor-to-matrix-to-tensor round trip
  reproduces the original only when the original is minor-symmetric
  (`tests/math/test_Tensor.cpp`, `SymmetricTensorMappingRoundTrip`).
- **Performance.** The contraction and rotation kernels are fully unrolled over
  the 81 (or 9) entries, avoiding inner loops. Element-wise add/subtract
  operations instead loop over the tensor capacity, trading hard-coded unrolling
  for support of arbitrary order and dimension. Pass preallocated work/pivot
  buffers and reuse output tensors across calls. `operator%` and
  `fill_orthotropic_elasticity` allocate on every call; `ddot`,
  `fill_isotropic_elasticity` and `invert_symmetric` do not.

---

## Usage in the framework

- **Self-consistent alloy homogenization** (`cl_Material_Alloy.cpp`): allocates a
  `Cell<Tensor<real>*>` of elasticity, compliance, Eshelby, Hill-polarization,
  influence, identity and work tensors; fills isotropic component stiffnesses,
  iterates effective `K, G` via `S % M`, `P % B`,
  `tensor::invert_symmetric(A, ...)`, and `Ck % A`
  (`Alloy::create_tables`, `Alloy::self_consistent_KG`).
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
- General Voigt/elasticity background: Hughes 2000 and Bathe 2016 (full citations in `doc/literature_references.md`).

---

## See Also

- `src/math/quaternion/doc/quaternion_usage_guide.md` - sibling math helper (rotations; produces the `R` for `rotate`)
