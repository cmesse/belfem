# Tensor Module Test Suite

**Date:** 2026-03-23

---

## Test Count by Suite

| Suite | Tests | Coverage |
|-------|-------|----------|
| `Tensor` | 27 | Construction, copy/move, self-assignment, scalar/matrix assign, access, layout, operators |
| `TensorDebug` | 14 | Wrong size, mismatch copy/move, index OOB, non-3333 operators |
| `TensorFill` | 7 | `fill(a,b)`, isotropic elasticity, known steel values, orthotropic, isotropic limit |
| `TensorConversion` | 8 | Compliance matrix, stiffness*compliance=I, mat-ten round-trip, Voigt convention pinning |
| `TensorKernel` | 26 | contract42 (vs ref, identity, zero, member, operator%), contract44 (vs ref, identity L/R), rotate42 (vs ref, identity, isotropic invariance, composition, high/low level), KC (vs ref, isotropic, scaling, symmetry), invert (round-trip, minor+major symmetry), physics (hydrostatic, pure shear) |
| `TensorKernelDebug` | 9 | Non-3333 for all kernels, wrong sizes |
| **Total** | **91** | |

---

## Key Design Decisions

### Reference Loop Validation
All unrolled kernels (`contract42`, `contract44`, `rotate42`, `kelvin_christoffel`) are validated against naive loop implementations in anonymous namespace. Index formula: `A[L*27 + K*9 + J*3 + I]` for `A(I,J,K,L)`.

### Rotation Composition Order
`rotate(B, R, A)` computes `A_mnop = R_im R_jn R_ko R_lp B_ijkl`. Sequential `rotate(A, R1, T); rotate(T, R2, result)` matches `rotate(A, R1*R2, result)` (NOT `R2*R1`). Tests use non-commuting rotations (Z + X axes) to catch order bugs.

### Tolerance for Large-Magnitude Tensors
Isotropic steel (E=200 GPa) gives tensor values ~O(1e11). Absolute tolerance `1e-9` fails for zero entries with ~1e-5 noise. Solution: tolerance relative to tensor max magnitude.

---

## @warning API Gotchas

- `fill(a, b)` is a MEMBER function, NOT `tensor::fill(tTen, a, b)`. The namespace function takes a raw pointer.
- `fill_isotropic_elasticity(E, nu)` is a MEMBER function.
- `invert_symmetric` is IN-PLACE: 3 arguments `(Tensor&, Vector<real>&, Vector<int_t>&)`.
- Identity tensor: `I_ijkl = 0.5 * (delta_ik * delta_jl + delta_il * delta_jk)` — symmetric identity, NOT simple Kronecker product.
- `operator%` between two tensors: `contract44`. Between tensor and matrix: `contract42`.
- Sized constructor `Tensor(3,3,3,3)` does NOT initialize data. Use `Tensor(3,3,3,3, 0.0)` for zeroed tensor.

---

## @note Grok Submission Errors (Rejected)

Grok's draft had 6 API errors: wrong function names (`tensor::fill(tTen, a, b)` instead of member), wrong signatures, wrong memory assumptions, code that wouldn't compile. All rejected; source code used as ground truth.

---

## @todo Future Work

- [ ] Fiber polarization tests (planned §7, deferred)
- [ ] Additional Kelvin-Christoffel directions (currently only 2D propagation tested)
