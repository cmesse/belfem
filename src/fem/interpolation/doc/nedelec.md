# Nédélec Edge Elements in BELFEM {#fem_interpolation_nedelec}

**Date:** 2026-04-14
**Module:** `src/fem/interpolation/nedelec`
**Purpose:** Theory, interface, and implementation conventions for BELFEM's Nédélec (edge) finite elements. This document covers the content that applies to *all* edge elements — their role in H(curl) formulations, the `EdgeFunction` base class, the lifecycle, the factory pattern, and the volume-element implementations (TRI3, TRI6, TET4, TET10, LINE3). The thin-shell variants (QUAD4TS, PENTA6TS) are specialized and documented separately in [nedelec_thinshell.md](nedelec_thinshell.md).

---

## 1. Why Nédélec elements

BELFEM's magnetodynamic solver uses the H-phi formulation (Messe et al. 2023, Arsenault et al. 2023). In conducting subdomains, the unknown is the magnetic field strength **H**, which must be curl-conforming (continuous tangentially across interfaces, discontinuous normally) and whose curl is the current density **J = ∇ × H**. Standard Lagrange (nodal) interpolation is H¹-conforming and imposes continuity of every component across element faces; it over-constrains the tangential-only continuity that **H** actually needs and produces spurious curl-free modes in multiply connected domains.

Nédélec's (1980) edge elements solve both problems. Their DOFs are *edge circulations*

```
DOF_k = ∫_{edge_k} H · dℓ
```

which make tangential continuity automatic across element faces (because adjacent elements share the same edge DOF) while leaving normal components free to jump. The interpolation space lies in H(curl) rather than H¹. For a systematic treatment see Monk (2003) and Nédélec (1980).

In BELFEM, Nédélec elements are used for every conducting block (plus cohomology-cut DOFs at the interface with scalar-phi regions). Air, ferromagnetic, and buffer regions use scalar Lagrange elements for φ. The coupling at the H-φ interface is handled through static condensation of edge DOFs onto nodal values (see `src/fem/kernel/doc/`).

---

## 2. The `EdgeFunction` interface

All BELFEM edge elements derive from the abstract base `EdgeFunction` at
[`cl_EF_EdgeFunction.hpp`](../nedelec/cl_EF_EdgeFunction.hpp). The class owns the usual FEM-element state and a small public API:

### 2.1 Protected state

| Member | Type | Meaning |
|---|---|---|
| `mJ` | `Matrix<real>` | Geometry Jacobian, **transposed** (`mJ(i,j) = ∂x_j / ∂ξ_i`) |
| `mInvJ` | `Matrix<real>` | Inverse of the transposed Jacobian. Column *j* of `mInvJ` is the physical-space gradient of reference coordinate *j* (`∇ξ_j`) |
| `mDetJ` | `real` | Determinant of the Jacobian (may be signed) |
| `mAbsDetJ` | `real` | Absolute value of `mDetJ`, used as the integration weight |
| `mE` | `Matrix<real>` | Edge basis values at current integration point, shape `(ndim × nDOFs)` |
| `mC` | `Matrix<real>` | Curl of the edge basis at current integration point, same shape as `mE` |
| `mSumW` | `real` | Sum of integration weights over the reference element |
| `mNumDofs` | `uint` | Number of edge DOFs this element carries |

The transposed-Jacobian convention is the same as the rest of BELFEM's interpolation module. A consequence is that the columns (not rows) of `mInvJ` are physical-space gradients of reference coordinates — useful to remember when reading basis construction code.

### 2.2 Public API

```cpp
virtual void link( Element * aElement ) = 0;
virtual void precompute( const Matrix< real > & aXi ) = 0;
virtual const Matrix< real > & E( const uint aIndex );
virtual const Matrix< real > & C( const uint aIndex = 0 ) = 0;
real det_J() const;
real abs_det_J() const;
real sum_w() const;
uint ndofs() const;
```

- `link()` is called once per element and computes everything that depends on physical node coordinates but not on the integration point: the Jacobian, per-edge gradient data, edge orientation signs.
- `precompute()` is called once per integration rule and caches everything that depends on the integration point but not on the physical element (reference-space shape values and derivatives).
- `E(k)` and `C(k)` return the basis and its curl at integration point *k*. They're called inside the element assembly loop.
- `det_J()` and `abs_det_J()` return the determinant cached by the most recent `E()` or `C()` call, for use as the integration weight.
- `ndofs()` is the DOF count; `sum_w()` is the total reference-element volume (currently consumed only by tests, e.g. `test_InterfaceOrientation`).

### 2.3 Lifecycle in element assembly

A typical assembly loop in an IWG looks like:

```cpp
// Once per element type, outside the element loop
EdgeFunctionFactory factory;
EdgeFunction * edgeFn = factory.create_edge_function( tElementType );
IntegrationData integration( tElementType );
integration.populate( tOrder, IntegrationScheme::GAUSS );
edgeFn->precompute( integration.points() );

// Per element
for ( Element * tElement : tElements )
{
    edgeFn->link( tElement );

    // Per integration point
    for ( uint k = 0; k < integration.number_of_integration_points(); ++k )
    {
        const Matrix<real> & E = edgeFn->E( k );  // basis (ndim × nDOFs)
        const Matrix<real> & C = edgeFn->C( k );  // curl (ndim × nDOFs)
        real wJ = integration.weights()(k) * edgeFn->abs_det_J();

        // Accumulate element matrices
        Ke += trans(C) * rho * C * wJ;   // stiffness from ∫ ρ (∇×E)·(∇×E) dV
        Me += trans(E) * mu0 * E * wJ;   // mass from ∫ μ₀ E·E dV
    }
}
```

The `edgeFn` and `integration` objects are cached by element type outside the loop — see `src/fem/kernel/cl_FEM_Group` for BELFEM's calculator machinery that does this caching at the block level.

---

## 3. The factory pattern

Edge functions are created through `EdgeFunctionFactory::create_edge_function(ElementType)` at
[`cl_EdgeFunctionFactory.hpp`](../cl_EdgeFunctionFactory.hpp). The factory is stateless; each call constructs a fresh edge function. Current dispatch:

| ElementType | Edge function class | Category | Status |
|---|---|---|---|
| `TRI3` | `EF_TRI3` | Volume (2D) | Production |
| `TRI6` | `EF_TRI6` | Volume (2D), 2nd order | ⚠️ Proof of concept (see §6.5) |
| `TET4` | `EF_TET4` | Volume (3D) | Production |
| `TET10` | `EF_TET10` | Volume (3D), 2nd order | ⚠️ Proof of concept (see §6.5) |
| `LINE3` | `EF_LINE3` | Higher-order edge (1D), 2nd order | ⚠️ Proof of concept (see §6.5) |
| `QUAD4TS` | `EF_QUAD4TS` | Thin shell — see [nedelec_thinshell.md](nedelec_thinshell.md) | Production |
| `PENTA6TS` | `EF_PENTA6TS` | Thin shell — see [nedelec_thinshell.md](nedelec_thinshell.md) | Production |
| `HEX8` | `EF_HEX8` | Volume (3D) | ⚠️ Rectangular bricks only (see §6.6) |
| `HEX8TS` | `EF_HEX8TS` | Thin shell — see [nedelec_thinshell.md](nedelec_thinshell.md) | Production |
| `HEX8TB` | `EF_HEX8TB` | Side-connector wall element | Production |

The DOF count per element type is also reported through the helper `num_nedelec_dofs(ElementType)` in [`fn_num_nedelec_dofs.hpp`](../nedelec/fn_num_nedelec_dofs.hpp):

| ElementType | DOF count | Notes |
|---|---|---|
| `LINE2` | 1 | Circulation along the line |
| `LINE3` | 2 | Linear + quadratic |
| `TRI3` | 3 | One per edge |
| `TRI6` | 8 | Six edges + two face bubbles |
| `TET4` | 6 | One per edge |
| `TET10` | 20 | Lowest order + quadratic edge + face |
| `QUAD4TS` | 2 | Tangential-only, per thin-shell convention |
| `PENTA6TS` | 6 | Three bottom + three top |

`num_nedelec_dofs` is consulted by the DOF allocation code when the IWG needs to know how many edge DOFs a given block contributes. When adding a new edge element type, update both the factory switch and this helper.

---

## 4. Edge orientation

BELFEM's meshes use *global* edges (shared by multiple elements), so each Nédélec basis function needs a sign for the element's local edge direction. The sign is retrieved with `aElement->edge_directions(mS)` inside `link()`. `mS` is a small fixed-size array (one real per edge) containing `+1` if the element's local edge direction matches the global edge direction, `-1` otherwise. All volume edge elements use this convention uniformly; see e.g. `EF_TET4::link()` at [cl_EF_TET4.cpp](../nedelec/cl_EF_TET4.cpp):94.

The sign enters the basis construction as an explicit factor. Without it, the circulation integral `∫ E_k · dℓ` would return the wrong sign on half of the shared edges, breaking tangential continuity. The orientation information comes from the mesh's edge construction (ordering of edge endpoints by node ID or similar); users don't need to touch it directly.

---

## 5. Circulation, unit-circulation normalization, and DOF sharing

### 5.1 The "DOF = circulation" convention

Every BELFEM edge element is normalized so that the DOF value equals the edge circulation. Concretely, for basis function `E_k` associated with edge *k*,

```
∫_{edge_j} E_k · dℓ = s_k · δ_{jk}
```

where `s_k` is the edge's orientation sign. With this normalization, the discrete H field is the simple linear combination

> **One documented exception: TRI6.** The normalization above holds for the first-order elements.
> `TRI6`'s two parent-edge functions each carry circulation `s/2` rather than unity — the hand
> integrals `H_0 - G_0 = 2 - 3t` and `H_1 - G_1 = 3t - 1` both give `1/2`, and the implementation
> agrees (`nedelec_derivation.md` §3.3, `nedelec/cl_EF_TRI6.cpp`). **This does not extend to
> TET10**, whose edge polynomials are twice the TRI6 pair and do integrate to unit circulation
> (`nedelec_derivation.md` §4.1). Do not renormalise one without the other.

```
H(x) = Σ_k h_k · E_k(x)
```

with coefficients `h_k` that are directly the edge circulations of **H** along each edge. Interpretation is unambiguous, and post-processing (current computation, line integrals, cut coefficients) doesn't need to carry around basis normalization factors.

### 5.2 Why this matters for shared edges

When two elements share an edge, the shared DOF must mean the same thing from both sides. If element A's basis were normalized so that `∫ E_k · dℓ = 1` and element B's basis were normalized so that `∫ E_k · dℓ = 2`, a DOF value of `1.0` would correspond to different physical H fields on the two sides — the assembled system would be inconsistent.

All BELFEM edge functions adhere to the unit-circulation convention. When adding a new edge element, verify this by hand (compute the circulation of its basis on a reference element) before wiring it into the factory. This is particularly important for reduced or custom bases such as the thin-shell elements in [nedelec_thinshell.md](nedelec_thinshell.md), where the basis normalization must be chosen specifically to reach unit circulation.

### 5.3 Curl as current density

The curl operator `C = ∇ × E` gives the current density basis:

```
J(x) = ∇ × H(x) = Σ_k h_k · (∇ × E_k(x)) = Σ_k h_k · C_k(x)
```

The stiffness matrix entry from `∫ ρ (∇×H)·(∇×H) dV` is therefore `trans(C) * ρ * C * dV` at each integration point, which is the form used in BELFEM's conductor IWG kernels in `src/fem/maxwell/matrices/mt_maxwell_h.cpp`.

---

## 6. Volume edge elements

BELFEM's volume Nédélec elements follow the canonical Whitney construction:

```
E_{ab} = N_a · ∇N_b - N_b · ∇N_a        (edge from vertex a to vertex b)
```

where `N_a`, `N_b` are the linear Lagrange shape functions of the two end vertices. This construction has unit circulation on its own edge, zero on every other edge, and zero curl on degenerate edges. The sign convention matches `s_k` from `edge_directions()`.

### 6.1 TET4

The lowest-order 3D edge element: four vertices, six edges, six DOFs. The Whitney basis expands to closed-form polynomials in `(ξ, η, ζ)` plus a barycentric `τ = 1 − ξ − η − ζ`. Both basis and curl are constant over the element (because `∇N_a` is constant), so the element is exactly integrated by a single-point Gauss rule — `mSumW = 1/6` (the reference tetrahedron volume). See [cl_EF_TET4.cpp](../nedelec/cl_EF_TET4.cpp).

### 6.2 TET10

The quadratic counterpart of TET4: ten nodes (four vertices + six edge midpoints), 20 DOFs. Twelve edge DOFs (two per geometric edge: a lowest-order plus a quadratic refinement) and eight face DOFs. The increased richness of the basis pays off in convergence order when the H field has significant variation inside an element.

### 6.3 TRI3 and TRI6

The 2D analogs. In BELFEM's 2D H-phi formulation, the conductor's H field is in-plane (H_x, H_y) while φ is scalar, and cohomology cuts provide the transport current coupling. TRI3 and TRI6 implement the Whitney and quadratic Nédélec bases on a triangle. `EF_TRI6` has eight DOFs: six edge (two per geometric edge) plus two interior, per the standard quadratic Nédélec count on a triangle.

### 6.4 LINE3

The 1D edge element used for higher-order 1D refinements. Rarely appears directly in the magnetic solve — it's used by the cohomology-cut machinery and in certain interface constructions.

### 6.5 ⚠️ Status of the second-order elements

The second-order Nédélec elements — **`TRI6`**, **`TET10`**, and **`LINE3`** — are present in the codebase as **proof-of-concept implementations**. They are not currently part of BELFEM's tested and supported feature set.

Concretely:

- **They have not been validated** against analytical or independently-computed reference solutions for an H-φ problem of meaningful complexity. The basis and curl expressions in the source files are derived correctly in principle, but no convergence study or patch test confirms that they integrate and assemble cleanly through the rest of the pipeline.
- **The IWG and assembly machinery** in `src/fem/maxwell` and `src/fem/kernel` is not exercised against quadratic edge DOFs in the routine test cases. Higher-order DOF allocation, hanging-node treatment at the H-φ interface, and cohomology-cut wiring may have latent bugs that only surface when 2nd-order elements are actually used in production.
- **Mesh handling** for quadratic faces and edges (orientation, slave-side reordering, midside node bookkeeping) has been written but, again, not stress-tested against the specific failure modes that 2nd-order Nédélec exposes.
- **Thin-shell counterparts** (e.g. a hypothetical quadratic `PENTA15TS`) do not exist; second-order support is volume-only.

In practice, **production runs use lowest-order elements** (`TRI3`, `TET4`) for the conducting volumes plus the dedicated thin-shell elements for the tape stack. If you need higher-order accuracy on the volume side, plan for additional validation work before relying on results — at minimum a manufactured-solution patch test on a simple geometry, plus careful inspection of the assembled Jacobian for the specific problem class. Treat the existing `EF_TRI6`, `EF_TET10`, and `EF_LINE3` implementations as a starting point, not as a guaranteed-correct black box.

### 6.6 ⚠️ QUAD and HEX edge elements: rectangular geometry only

The hexahedral edge element `EF_HEX8` (and any future quadrilateral one) has one
geometry restriction that simplex elements do not: use it only on perfectly rectangular
elements, meaning rectangles in 2D and bricks in 3D, either axis-aligned or rigidly
rotated. Do not use trapezoidal, sheared, or otherwise distorted shapes. This is not a
BELFEM implementation limit. It is a known property of tensor-product H(curl)/H(div)
elements, and the failure is silent.

**The mechanism.** On triangles and tetrahedra, the reference-to-physical map is
affine, so the covariant (Piola) transformation carries the reference Nédélec space onto
a physical polynomial space with the same approximation power. On general quadrilaterals
and hexahedra, the map is bi-/trilinear rather than affine. The transform then stops
preserving the polynomial space: the mapped basis loses completeness, along with the
interpolation estimates used by the convergence theory. Monk 2003 builds the
hexahedral Nédélec theory on the affine case (Monk, §6.1) and warns that non-affine
hexahedral maps can give non-optimal convergence, or even non-convergence, in H(curl; Ω)
(Monk, §8.2–8.3).

**How bad it gets.** Arnold, Boffi & Falk and Falk, Gatto & Monk are the key studies.
They show that standard mapped families on general, non-affine quadrilateral meshes lose
approximation order, with the lowest-order elements, the ones BELFEM uses in production,
hit hardest. In the H(div) counterpart to this problem, Boffi et al. 2013, Remark 2.5.5,
state that the divergence does not converge for k = 0 on general quadrilateral meshes.
Falk et al. 2011 extend the negative results to 3D hexahedral H(curl) elements. Arnold
et al. 2001 and 2002 give the underlying approximation theory, including the scalar
serendipity case, and Arnold et al. 2005 treats H(div).

**Why this is dangerous in practice.** Nothing asserts. A distorted-hex h-φ model can
assemble, solve, and produce plausible fields. The fields are still wrong, or they
converge at a reduced rate that refinement along the same distorted mesh family will not
repair. Treat this as a meshing rule:

- **Conducting volumes on hex meshes: rectangular bricks only.** Use structured,
  axis-aligned bricks, or a rigidly rotated version of the same mesh. Grading brick sizes
  is fine because each element remains rectangular; shear and trapezoids to fit the
  boundary are not.
- **Geometry does not fit rectangular bricks? Use `TET4`.** Simplex edge elements have no
  such restriction. Any non-degenerate tetrahedron is an affine image of the reference
  element.
- Plain-quadrilateral volume edge elements deliberately do not exist in BELFEM; the
  `EdgeFunctionFactory` refuses `QUAD4`. Any future implementation would have to follow
  the same rectangular-geometry rule.
- The thin-shell elements (`QUAD4TS`, `HEX8TS`, `PENTA6TS`) are separate constructions
  with their own geometry handling. See `nedelec_thinshell.md`; this section is about
  volume elements.

Scalar (Lagrange) elements are far more forgiving: bilinear/trilinear scalar elements keep
their first-order energy-norm convergence on arbitrary non-degenerate quads/hexes, which is
why thermal problems on such meshes are acceptable. The restriction here is specific to the
edge (and face) element spaces of the magnetic solve.

**References:** Monk 2003, §6.1 and §8.2–8.3; Boffi et al. 2013, §2.2.4 and §2.5.5
(Remark 2.5.5); Arnold et al. 2001, 2002, 2005; Falk et al. 2011. Full citations with DOIs
in `doc/literature_references.md`.

---

## 7. Implementation checklist for new volume edge elements

When adding a new volume Nédélec element type:

1. **Create the class** in `src/fem/interpolation/nedelec/`, deriving from `EdgeFunction`. Set `mNumDofs` and `mSumW` in the constructor, and size `mE`, `mC` appropriately.
2. **Implement `link()`** to build the Jacobian, inverse Jacobian, gradients of reference coordinates, and edge orientation signs. If the basis is constant over the element (affine mapping), compute the curl here too.
3. **Implement `precompute()`** to cache reference-space shape data at each integration point (the structure of `mG` and `mH` in `EF_TET4::precompute()` is a useful template).
4. **Implement `E(k)`** as a loop over edges, each row computing `s_k · (affine polynomial in reference coordinates) · (physical gradient)`. The polynomials come from the Whitney expansion.
5. **Implement `C(k)`**. On affine elements this can return `mC` unchanged from `link()`. On non-affine elements it needs to be re-evaluated per integration point.
6. **Register the element** in `EdgeFunctionFactory::create_edge_function` and in `num_nedelec_dofs`.
7. **Verify unit circulation** by hand or with a simple test: integrate the basis along each reference edge and check that the result is `δ_{jk}` up to orientation.
8. **Update this documentation** with the DOF count, file location, and any element-specific notes.

For thin-shell variants, see [nedelec_thinshell.md](nedelec_thinshell.md) for the additional considerations around dimensional reduction, pseudo-inverse Jacobians, and shared-edge compatibility.

---

## 8. Common pitfalls

- **Transposed Jacobian.** BELFEM stores `mJ(i,j) = ∂x_j/∂ξ_i`, so `mInvJ.col(j) = ∇ξ_j`. Mis-reading this as the canonical Jacobian leads to swapped rows and columns and wrong curls. The memory is: *columns* of `mInvJ` are physical gradients of reference coordinates.
- **Forgetting `s_k`.** Without the edge orientation sign, the basis has unit circulation on its own edge but with the wrong sign on half of the shared edges. Tangential continuity is broken silently — no assertion fires, results just become nonsensical.
- **Stale `mDetJ`.** `mDetJ` is set by the most recent `E()` or `C()` call; if assembly needs both the basis value and the determinant at the same integration point, call `E(k)` (or `C(k)`) before reading `det_J()`. For elements where the determinant is set in `link()` (affine case) rather than per-integration-point, this is moot.
- **Using a Lagrange shape-function factory.** `EdgeFunctionFactory` is separate from `InterpolationFunctionFactory`. Nédélec bases are not Lagrange and can't be constructed through the Lagrange factory path.
- **Assuming the reference volume equals 1.** For tetrahedra, `mSumW = 1/6`; for triangles, `mSumW = 1/2`. The integration weight inside the assembly loop is `integration.weights()(k) * abs_det_J()` — `abs_det_J` picks up the physical volume, the reference weight is already normalized to the reference-element measure.

---

## 9. References

- **Nédélec, J.-C. (1980)**, "Mixed finite elements in ℝ³." *Numerische Mathematik* **35**, 315–341. The foundational paper that introduced edge elements for **H**(curl).
- **Monk, P. (2003)**, *Finite Element Methods for Maxwell's Equations*. Oxford University Press. Chapter 5 covers the Nédélec construction, the edge DOF interpretation, and convergence theory.
- **Arnold, D. N. (2018)**, *Finite Element Exterior Calculus*. SIAM. Places Nédélec elements in the FEEC framework and derives them systematically from Whitney forms.
- **Messe, C. et al. (2023)**, "BELFEM: a special-purpose finite-element code for the magnetodynamic modeling of high-temperature superconducting tapes." *Supercond. Sci. Technol.* **36**, 114001. Describes BELFEM's overall H-phi architecture and how edge elements fit into it.

For thin-shell element theory specifically, see [nedelec_thinshell.md](nedelec_thinshell.md). For the interpolation module in general — including non-Nédélec shape functions, integration data, and the InterpolationFunction factory — see the [module index](@ref fem_interpolation_index) and [interpolation_usage_guide.md](interpolation_usage_guide.md).
