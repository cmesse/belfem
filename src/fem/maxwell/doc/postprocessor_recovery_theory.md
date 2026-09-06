# Maxwell Postprocessor: Recovery Theory {#fem_maxwell_postprocessor_recovery_theory}

**Date:** 2026-01-20
**Purpose:** Mathematical foundation and implementation details of field recovery in Maxwell postprocessor
**Module:** `src/fem/maxwell`

> **⚠️ AUTHORITATIVE DOCUMENT:** This document defines the mathematical model for Maxwell field recovery. Other documents may reference or summarize, but NOT redefine the recovery method.

**See Also:**
- [Maxwell Module README](@ref fem_maxwell_index) - Overview and quick reference

---

## Glossary of Electromagnetic Symbols

| Symbol | Meaning | Units (SI) | Notes |
|--------|---------|------------|-------|
| **H** | Magnetic field intensity | A/m | Primary variable in h-φ formulation |
| **B** | Magnetic flux density | T (Tesla) | B = μ·H |
| **J** | Current density | A/m² | J = ∇×H for conductors |
| **Jc** | Critical current density | A/m² | Superconductor threshold |
| **φ** | Scalar magnetic potential | A (Ampere) | In air/ferro regions |
| **A** | Vector magnetic potential | Wb/m | B = ∇×A (not used in h-φ) |
| **μ** | Magnetic permeability | H/m | μ = μ₀·μᵣ, nonlinear for ferro/SC |
| **μ₀** | Vacuum permeability | 4π×10⁻⁷ H/m | Constant |
| **σ** | Electrical conductivity | S/m | For Ohmic conductors |
| **T** | Temperature | K (Kelvin) | For thermal coupling |
| **σ\*** | Recovered field | (varies) | SPR-projected quantity |
| **σₕ** | Raw FEM field | (varies) | Discontinuous across elements |
| **p(x)** | Polynomial basis | - | Monomials for recovery |
| **c** | Polynomial coefficients | - | From least-squares fit |
| **V** | Gramian/Vandermonde matrix | - | SPD, size m×m |

---

## Overview

The Maxwell postprocessor implements a **Superconvergent Patch Recovery (SPR)** method to project discontinuous element-wise electromagnetic fields onto smooth nodal values suitable for visualization and output. This document explains the mathematical theory, implementation approach, and domain-specific aspects.

---

## 1. The Recovery Problem

### 1.1 Why Recovery is Needed

Finite element solutions provide continuous primary variables:
- **Lagrange formulation:** Scalar potential φ (nodal)
- **Nédélec formulation:** Edge-based DOFs for H-field

However, derived quantities (gradients, curls, material-dependent fields) are:
- ✅ Accurate at integration points (Gauss quadrature)
- ❌ **Discontinuous** across element boundaries
- ❌ Unsuitable for direct visualization

**Recovery goal:** Construct smooth, continuous fields with higher accuracy by exploiting superconvergence properties at integration points.

### 1.2 Classical Approaches

| Method | Description | Accuracy | Cost |
|--------|-------------|----------|------|
| **Direct averaging** | Average element values at nodes | O(h^p) | Very low |
| **L² projection** | Global least-squares fit | O(h^(p+1)) | High (global solve) |
| **SPR (Zienkiewicz-Zhu)** | Local polynomial fit to patches | O(h^(p+1)) | Low (local solves) |
| **Gradient recovery (ZZ)** | SPR + error estimation | O(h^(p+1)) | Medium |

BELFEM uses **SPR** for optimal balance of accuracy and computational cost.

---

## 2. Superconvergent Patch Recovery (SPR)

### 2.1 Mathematical Formulation

For each node **i**, recover a smooth field **σ\*** from discontinuous element values **σₕ**:

**Step 1: Define polynomial approximation**
```
σ*(x) = p(x)ᵀ · c
```

where:
- **p(x)** = polynomial basis vector (monomials)
- **c** = unknown coefficients to determine

**Step 2: Weighted least-squares on element patch**

Minimize the error over all integration points in the patch:
```
min  Σₖ wₖ · ||σₕ(xₖ) - p(xₖ)ᵀ·c||²
```

This yields the normal equations:
```
V · c = b

where:
  V = Σₖ wₖ · p(xₖ) · p(xₖ)ᵀ    (Gramian/Vandermonde matrix)
  b = Σₖ wₖ · p(xₖ) · σₕ(xₖ)     (weighted right-hand side)
```

**Step 3: Solve and evaluate**
```
V · c = b           (solve via Cholesky: V is SPD)
σ*(x₀) = p(x₀)ᵀ·c   (evaluate at node location)
```

### 2.2 Polynomial Basis Selection

BELFEM uses **complete polynomial spaces** matching the element interpolation order:

#### 2D Polynomial Bases
```
Order 1:  p(x,y) = [1, x, y]ᵀ                              (3 terms)
Order 2:  p(x,y) = [1, x, y, x², xy, y²]ᵀ                  (6 terms)
Order 3:  p(x,y) = [1, x, y, x², xy, y², x³, x²y, xy², y³]ᵀ   (10 terms)
Order 4:  p(x,y) = [...up to x⁴, y⁴]ᵀ                      (15 terms)
```

#### 3D Polynomial Bases
```
Order 1:  p(x,y,z) = [1, x, y, z]ᵀ                         (4 terms)
Order 2:  p(x,y,z) = [1, x, y, z, x², xy, y², yz, z², zx]ᵀ  (10 terms)
Order 3:  [1, ..., x³, x²y, xy², y³, ...]ᵀ                (20 terms)
Order 4:  [... up to x⁴, y⁴, z⁴]ᵀ                          (35 terms)
```

**Implementation:** `cl_FEM_Postprocessor.hpp:226-366` (inline functions `poly1_2d` through `poly4_3d`)

### 2.3 Weight Function

BELFEM uses **element volume weighting**:
```
wₖ = Vₑ    (volume of element containing integration point k)
```

This ensures:
- Proper scaling across different element sizes
- Numerical stability of the Gramian matrix
- Consistent units for different field types

**Alternative weights considered but not used:**
- Inverse distance: `w = 1/||x-x₀||²` (commented out in code)
- Uniform: `w = 1` (less stable)

**Code reference:** `cl_MaxwellPostprocessor.cpp:730`, `cl_FEM_Postprocessor.cpp:285`

### 2.4 Patch Construction

A **patch** consists of all elements surrounding a node:

**Basic patch (depth=0):**
- All elements connected to the node
- Includes duplicates for periodic boundary conditions

**Extended patch (depth>0):**
- Iteratively adds neighbor elements
- Controlled by `mRecoveryDepth` parameter
- Improves conditioning but increases cost exponentially

**Code reference:** `cl_MaxwellPostprocessor.cpp:641-691`

---

## 3. BELFEM Implementation Details

### 3.1 Two Recovery Modes

#### Mode 1: Per-Node Recovery (the SPR pass)
**Used by:** the base `Postprocessor::run()`, which `MaxwellPostprocessor::run()` calls first
(`cl_MaxwellPostprocessor.cpp:388-390`) before its own per-type work. The recovery itself lives in
`Postprocessor::recover_fields()`, not in the Maxwell subclass.

**Algorithm:**
```
For each node i:
  1. Build element patch (local neighborhood)
  2. Initialize V = 0, b = 0
  3. For each element e in patch:
       For each integration point k:
         - Compute x_k (physical coordinates)
         - Evaluate p(x_k) (polynomial basis)
         - Compute σ_h(x_k) (field value via compute())
         - Accumulate: V += w·p·pᵀ, b += w·p·σ_h
  4. Solve: posv(V, b)  → coefficients c
  5. Evaluate: σ*(x_i) = p(x_i)ᵀ·c
```

**Advantages:**
- Flexible patch sizes (handles irregular meshes)
- Can vary recovery depth per region
- Straightforward to understand

**Disadvantages:**
- Visits elements multiple times (once per adjacent node)
- Per-node matrix solves (function call overhead)

**Code reference:** `cl_MaxwellPostprocessor.cpp:636-770`

---

#### Mode 2: Batch Recovery (`recover_fields`)
**Used by:** Parent `Postprocessor` class (alternative implementation)

**Algorithm:**
```
Phase 1: Precompute Gramian matrices (compute_node_matrices, called once)
  For each element e:                          ← first element loop
    For each integration point k:
      - Compute x_k = Nᵀ·coords (physical coordinates)
      - Evaluate p(x_k) (polynomial basis)
      - Compute tA = w·p·pᵀ (outer product)
      - Scatter tA to all owned adjacent nodes

Phase 2: Assemble RHS (recover_fields, called every timestep)
  For each element e:                          ← second element loop (redundant geometry work)
    For each integration point k:
      - Compute x_k = Nᵀ·coords (same as Phase 1)
      - Evaluate p(x_k) (same as Phase 1)
      - Compute σ_h(x_k) via compute()
      - Compute tB = w·p·σ_h
      - Scatter tB to all owned adjacent nodes

Phase 3: Per-node solve (recover_fields, continued)
  For each node:
    - Copy Gramian: tA = mNodeMatrices(k)      ← copy because posv overwrites
    - Solve: posv(tA, b)                       ← full Cholesky re-factorization every call
    - Evaluate σ*(x_i) = p(x_i)ᵀ·c
```

**Advantages:**
- Separates Gramian assembly from RHS assembly
- Can batch node solves

**Disadvantages:**
- Each element is visited **twice** (once per phase), with redundant coordinate
  transforms and polynomial evaluations
- Gramian factorization is repeated every call to `recover_fields()`, even though
  the Gramian matrices never change (they depend only on geometry)
- Requires more memory (store all node Gramian matrices)
- Less flexible for varying patch depths

**Code reference:** `cl_FEM_Postprocessor.cpp:219-487`

---

#### Mode 3: Gradient Recovery (`Gradient` class)
**Used by:** `cl_Gradient.cpp` (independent implementation)

**Algorithm:**
```
For each node i:
  1. Reset V = 0, C = 0
  2. For each element e connected to node i (including duplicates):
       For each integration point k:
         - Compute x_k (physical coordinates)
         - Compute J = dNdXi · X (Jacobian)
         - Compute gradient: g = J⁻¹ · dNdXi · φ_e
         - Evaluate p(x_k) (polynomial basis)
         - Accumulate: V += p·pᵀ, C += p·gᵀ
  3. Solve: gesv(V, C, pivot)
  4. Evaluate: grad*(x_i) = p(x_i)ᵀ · C
```

**Note:** Uses `gesv` (general LU) instead of `posv` (Cholesky), despite V being
SPD. This is functionally correct but suboptimal.

**Disadvantages:**
- Per-node: visits each element once per adjacent node (same redundancy as Mode 1)
- Rebuilds Vandermonde matrix from scratch every call
- Allocates `mPhi` and `mElX` per element via `set_size()` inside the loop

**Code reference:** `cl_Gradient.cpp:100-226`

---

### 3.2 Numerical Properties

**Gramian Matrix V:**
- **Type:** Symmetric positive definite (SPD)
- **Size:** m × m, where m = number of polynomial coefficients (3-35)
- **Condition number:** Typically well-conditioned for volume-weighted patches
- **Solver:** Cholesky factorization via LAPACK `posv` (optimal for SPD)

**Failure modes:**
- Singular V: patch too small or all integration points collinear
- Ill-conditioned V: highly stretched elements or large aspect ratios
- Recovery: code does not explicitly check; `posv` may fail silently

**Potential improvement:** Add condition number check or SVD fallback.

---

## 4. Maxwell-Specific Field Recovery

### 4.1 Field Computation at Integration Points

The `compute()` function computes physical quantities from primary DOFs:

| Domain Type | Primary DOFs | Field Computation | Code Reference |
|-------------|--------------|-------------------|----------------|
| **Air** | φ (nodal) | **H** = -∇φ = -B(ξ)·φ | `compute_air:909-917` |
| **Ferro** | φ (nodal) | **H** = -∇φ<br>**B** = μ(H)·H | `compute_ferro:919-932` |
| **Conductor** | edge_h (Nédélec) | **H** = E(ξ)·edge_h<br>**B** = μ(H)·H<br>**J** = C(ξ)·edge_h = ∇×H | `compute_conductor:934-949` |
| **Superconductor** | edge_h (Nédélec) | H, B, J (as above)<br>**J/Jc** = J / Jc(B,T,θ) | `compute_superconductor:1008-1060` |
| **Thin-shell SC** | edge_h + φₘ, φₛ | H, B, J, J/Jc<br>(includes normal field from master/slave) | `compute_superconductor_ts:1063-1148` |

**Key operations:**
- **B(ξ)** = gradient operator in reference coordinates (Lagrange elements)
- **E(ξ)** = Nédélec edge interpolation (H(curl) space)
- **C(ξ)** = Nédélec curl operator
- **μ(H)** = nonlinear permeability from material law

### 4.1.1 Linearity Analysis of Field Computations

For precomputed weight optimization, it is critical to know which fields are
**linear** in the DOFs (can be expressed as `M · DOFs`) and which are **nonlinear**:

| Domain | Field | Computation | Linear in DOFs? |
|--------|-------|-------------|-----------------|
| Air | **H** | `-B(ξ) · φ` | **Yes** |
| Ferro | **H** | `-B(ξ) · φ` | **Yes** |
| Ferro | **B** | `μ(‖H‖) · H` | **No** (μ depends on ‖H‖) |
| Conductor | **H** | `E(ξ) · edge_h` | **Yes** |
| Conductor | **B** | `μ(‖H‖) · H` | **No** |
| Conductor | **J** | `C(ξ) · edge_h` | **Yes** |
| Superconductor | **H** | `E(ξ) · edge_h` | **Yes** |
| Superconductor | **B** | `μ(‖H‖) · H` | **No** |
| Superconductor | **J** | `C(ξ) · edge_h` | **Yes** |
| Superconductor | **J/Jc** | `J / Jc(‖B‖, θ, T)` | **No** |

**Consequence:** H and J can be recovered using precomputed element-wise weight
matrices (see Section 13). B and J/Jc require integration-point-level evaluation
due to material nonlinearity.

### 4.2 Material Nonlinearity

For ferromagnetic and superconducting materials:
```
B = μ(H, T, x) · H
```

The permeability is evaluated at each integration point:
- **Ferro:** μ = μ(||H||) from B-H curve
- **Superconductor:** μ = μ₀ (linear), but Jc = Jc(||B||, θ, T, x)

**Thermal coupling:**
If thermal kernel is active, temperature T is interpolated from thermal DOFs at the same integration point.

**Code reference:** `cl_MaxwellPostprocessor.cpp` — the `mMaterial->mu( norm( mH ) )` evaluation sites (search `mMaterial->mu`)

### 4.3 Thin-Shell Specifics

For thin-shell superconductors, the recovery includes:

**Normal field contribution:**
```
bn = -0.5·μ₀·(Bm(ξ)·φₘ + Bs(ξ)·φₛ)
```
where φₘ, φₛ are potentials from master/slave elements.

**Total field:**
```
b = B_tangential + bn · n
```

**Angle for critical current:**
```
β = acos(n · b / ||b||)    (angle between normal and total field)
Jc = Jc(||b||, β, T)
```

This accounts for the angular dependence of critical current in anisotropic superconductors.

**Code reference:** `cl_MaxwellPostprocessor.cpp` — the `mJJc /= mMaterial->jc(...)` normalization block (search `mJJc`)

### 4.4 T-Matrices: Constraint Matrices for Hanging Nodes

**What are T-Matrices?**

T-matrices (transformation/constraint matrices) enforce **continuity of tangential H-field** at mesh refinement interfaces where hanging nodes or hanging facets exist. In h-refined meshes (adaptive refinement), constrained DOFs on fine-mesh boundaries must satisfy compatibility with coarse-mesh DOFs.

**Why are they needed?**

Nédélec edge elements naturally ensure **tangential continuity** of H within conforming meshes. However, at hanging-node interfaces (where fine elements meet coarse elements):

```
Coarse element edge:  [--------]
                         |  |
Fine element edges:   [--][--]   ← Hanging node in middle
```

Without constraints, the two fine-element edge DOFs are independent, potentially creating **discontinuity in tangential H** at the interface. This violates Maxwell's equations:

```
n × (H₁ - H₂) = 0    (tangential H must be continuous across interfaces)
```

T-matrices enforce this constraint by expressing fine-mesh DOFs as linear combinations of coarse-mesh DOFs.

**Mathematical Formulation:**

For constrained DOFs (hanging), express them in terms of independent DOFs:

```
h_constrained = T · h_independent

where:
  T = sparse constraint matrix (m_constrained × n_independent)
  h_constrained = dependent edge DOFs
  h_independent = master edge DOFs
```

**System-level projection:**

When assembling stiffness matrix K and RHS f, constrained DOFs are eliminated via:

```
K_reduced = Tᵀ · K · T
f_reduced = Tᵀ · f
```

This reduces the system size and ensures constraints are satisfied automatically.

**Maxwell-Specific Implementation:**

BELFEM's `maxwell::TMatrix` class (`cl_Maxwell_TMatrix.hpp/cpp`) computes constraint coefficients for **thin-shell interfaces** where master and slave elements meet:

**Key operations:**

1. **Facet linking:** Link TRI6 facet (on TET10 master/slave elements)
   ```cpp
   mCalc->link(aFacet);  // cl_Maxwell_TMatrix.cpp:46
   ```

2. **Gradient computation:** Compute nabla operator (inverse Jacobian)
   ```cpp
   compute_nabla(aIndex);  // cl_Maxwell_TMatrix.cpp:116-138
   mNabla = J⁻¹  (3×4 matrix for tetrahedral reference element)
   ```

3. **Nédélec edge functions:** Evaluate two Nédélec basis functions on facet
   ```cpp
   compute_nedelec_function(i, j, k, aIndex);  // cl_Maxwell_TMatrix.cpp:141-175

   E₁ = (16·λⱼ·λₖ)·∇λᵢ + (-8·λᵢ·λₖ)·∇λⱼ + (-8·λᵢ·λⱼ)·∇λₖ
   E₂ = (-8·λⱼ·λₖ)·∇λᵢ + (16·λᵢ·λₖ)·∇λⱼ + (-8·λᵢ·λⱼ)·∇λₖ
   ```
   where λᵢ, λⱼ, λₖ are barycentric coordinates.

4. **Integration over facet:** Compute constraint coefficients via Gauss quadrature
   ```cpp
   // For each integration point p on facet:
   T(i,j) -= w(p) · (E × n) · ∇φⱼ · dS(p)  // cl_Maxwell_TMatrix.cpp:86-94

   where:
     E = Nédélec edge function
     n = outward normal vector
     ∇φⱼ = gradient of Lagrange basis on slave element
     w(p) = integration weight
     dS(p) = surface Jacobian
   ```

5. **Return constraint vector:** Returns 12-element vector (2 edge functions × 6 facet edges)
   ```cpp
   process(mesh::Facet*) → Vector<real>(12)  // cl_Maxwell_TMatrix.cpp:43-112
   ```

**Parent Class: General T-Matrix Storage**

The parent class `fem::Tmatrix` (`cl_FEM_Tmatrix.hpp/cpp`) stores constraint matrices in **CSR (Compressed Sparse Row)** format:

```cpp
// Stored data:
uint* mPointers;  // Row pointers (CSR format)
uint* mIndices;   // Column indices
real* mValues;    // Non-zero constraint coefficients

// Key operations:
project(Vector<real>& A, Vector<real>& B);  // B = T·A
project(Matrix<real>& A, Matrix<real>& B);  // B = Tᵀ·A·T
```

**Usage in Recovery:**

During postprocessing, T-matrices ensure recovered fields respect hanging-node constraints:

1. **DOF retrieval:** When `update_dofs()` gathers edge DOFs, constrained DOFs are automatically computed from master DOFs
2. **Field evaluation:** Nédélec interpolation (`E(ξ)·edge_h`) uses constrained DOF values
3. **Continuity:** Recovered H, B, J fields are automatically continuous across refinement interfaces

**When are T-Matrices Created?**

T-matrices are precomputed during mesh initialization if:
- ✅ Mesh has hanging nodes (from h-adaptive refinement)
- ✅ Nédélec elements are used (edge-based DOFs)
- ✅ HDF5 file contains T-matrix data (`read_node_tmatrix_from_hdf5`)

For uniform meshes without hanging nodes, T-matrices are not needed.

**Debugging and Validation:**

The Maxwell TMatrix implementation includes **assertion checks** (`cl_Maxwell_TMatrix.cpp:70-79`):
```cpp
// Verify master and slave integration points coincide geometrically:
BELFEM_ASSERT(||P_facet - P_master|| < 1e-10, "Master point mismatch");
BELFEM_ASSERT(||P_facet - P_slave|| < 1e-10, "Slave point mismatch");
```

**Code References:**

| File | Purpose | Key Lines |
|------|---------|-----------|
| `cl_Maxwell_TMatrix.hpp` | Maxwell-specific header | 23-54 |
| `cl_Maxwell_TMatrix.cpp` | Constraint computation | 43-175 |
| `cl_FEM_Tmatrix.hpp` | General CSR storage | 22-82 |
| `cl_FEM_Tmatrix.cpp` | T·v and Tᵀ·A·T operations | (implementation) |
| `hangingnodes.cpp` | Example usage | 25, 46 |

**Literature Reference:**

The mathematical foundation for thin-shell interface constraints and cohomology-based cut generation is covered in:

- **Alves et al. 2022b, Section 3**: "Cohomology basis computation and cut generation"
- **Messe et al. 2023, Section 2.4**: "Interface coupling via static condensation"

The general theory of hanging-node constraints in FEM is covered in:

- **Zienkiewicz & Taylor, Vol. 1, §10.5**: "Hierarchical and h-adaptive meshes"
- **Brenner & Scott (2008)**, *The Mathematical Theory of Finite Element Methods*, Ch. 4: "Constraint elimination"

---

## 5. Theoretical Foundation: Literature References

### 5.1 Primary References

**Superconvergent Patch Recovery:**

1. **Zienkiewicz, O.C. & Zhu, J.Z.** (1992)
   "The superconvergent patch recovery and a posteriori error estimates. Part 1: The recovery technique"
   *International Journal for Numerical Methods in Engineering*, **33**(7), 1331-1364
   **DOI:** 10.1002/nme.1620330702

2. **Zienkiewicz, O.C. & Zhu, J.Z.** (1992)
   "The superconvergent patch recovery and a posteriori error estimates. Part 2: Error estimates and adaptivity"
   *International Journal for Numerical Methods in Engineering*, **33**(7), 1365-1382
   **DOI:** 10.1002/nme.1620330703

3. **Zienkiewicz, O.C. & Zhu, J.Z.** (1987)
   "A simple error estimator and adaptive procedure for practical engineering analysis"
   *International Journal for Numerical Methods in Engineering*, **24**(2), 337-357
   **DOI:** 10.1002/nme.1620240206

**Textbook treatment:**

4. **Zienkiewicz & Taylor** (2000)
   *The Finite Element Method*, Vol. 1, 5th ed.
   Chapter 15: "Errors and Recovery" (pp. 493-544)
   **Sections:** §15.2 (Superconvergence), §15.4 (SPR), §15.6 (Error estimates)

### 5.2 Related FEM Theory

**Superconvergence theory:**

5. **Babuška, I. & Strouboulis, T.** (2001)
   *The Finite Element Method and its Reliability*
   Oxford University Press, Chapter 9

**Moving Least Squares (MLS) connection:**

6. **Lancaster, P. & Salkauskas, K.** (1981)
   "Surfaces generated by moving least squares methods"
   *Mathematics of Computation*, **37**(155), 141-158

### 5.3 Maxwell-Specific References

**Nédélec element postprocessing:**

7. **Monk, P.** (2003)
   *Finite Element Methods for Maxwell's Equations*
   Oxford University Press, Chapters 5-6

8. **Bossavit, A.** (1998)
   *Computational Electromagnetism: Variational Formulations, Complementarity, Edge Elements*
   Academic Press

**BELFEM implementation papers:**

9. **Messe et al. (messe2023.txt)** - BELFEM core architecture and h-φ formulation
10. **Arsenault et al. (arsenault2023.txt)** - Magnetodynamic h-φ coupling and postprocessing

---

## 6. Superconvergence Properties

### 6.1 Why Integration Points are Special

**Theorem (Herrmann, 1978):** For polynomial finite elements of degree p, derivatives (gradients, curls) are superconvergent at Gauss integration points:

```
||σ_exact - σ_h||_∞ = O(h^(p+1))    at Gauss points
||σ_exact - σ_h||_∞ = O(h^p)        elsewhere
```

**Consequence:** Sampling at integration points and fitting a polynomial of degree p recovers one extra order of accuracy.

**Reference:** Zienkiewicz Vol. 1, §15.2

### 6.2 Recovery Convergence Rate

**Standard FEM:**
```
||u - u_h|| = O(h^(p+1))           (solution)
||∇u - ∇u_h|| = O(h^p)             (gradient/curl)
```

**With SPR:**
```
||∇u - ∇u_h*|| = O(h^(p+1))        (recovered gradient/curl)
```

**Practical implication:** For linear elements (p=1), gradients improve from O(h) to O(h²).

### 6.3 Conditions for Superconvergence

SPR achieves optimal convergence when:
1. ✅ Mesh is sufficiently regular (not highly distorted)
2. ✅ Polynomial degree matches element order
3. ✅ Patch contains enough sampling points (at least m = number of coefficients)
4. ✅ Weights are properly scaled (volume-weighting ensures this)

**BELFEM satisfies all conditions** for typical electromagnetics meshes.

---

## 7. Comparison: Recovery vs. Filter

BELFEM's postprocessor applies two smoothing mechanisms internally (there are
no user-callable `process_*()` entry points — the class interface is `run()` /
`initialize()`, driven by the controller at every saved timestep):

| Property | Patch recovery (SPR) | Inverse-distance filtering |
|----------|---------------------|-------------------|
| **Method** | Polynomial least-squares (SPR) | Inverse-distance weighted average |
| **Accuracy** | O(h^(p+1)) | O(h^p) |
| **Cost** | Medium (small matrix solve) | Low (simple averaging) |
| **Smoothness** | C⁰ continuous, polynomial | C⁰ continuous, piecewise |
| **Use case** | Production, accurate visualization | Quick/debug visualization |
| **Code** | `cl_MaxwellPostprocessor.cpp` (search `patch recovery`) | `cl_MaxwellPostprocessor.cpp` |

**Filter formula:**
```
σ*(x₀) = Σₖ wₖ·σₕ(xₖ) / Σₖ wₖ
where wₖ = 1/||x_k - x₀||²
```

This is **Shepard's method** (1968), a simple form of Moving Least Squares (MLS) with inverse-distance weighting.

---

## 8. MPI Parallelization Strategy

### 8.1 Domain Decomposition

Each MPI rank processes:
- **Owned nodes:** Nodes where `node->owner() == mCommRank`
- **Aura elements:** Neighboring elements needed to build patches

**Patch assembly:**
- Local (no communication required)
- Each rank independently builds V and b for its owned nodes

### 8.2 Communication Pattern

**Phase 1: Node distribution** (`collect_nodes`)
- Rank 0 identifies which nodes belong to each rank
- Broadcasts node IDs to all ranks
- Each rank creates local `mNodes` list

**Phase 2: Field computation** (`run`)
- Parallel (no communication)
- Each rank recovers fields for owned nodes

**Phase 3: Global assembly** (`collect_fields`)
- All ranks send their nodal field values to rank 0
- Rank 0 assembles into global field vectors
- Uses blocking `send`/`collect` (could be optimized)

**Code reference:** `cl_MaxwellPostprocessor.cpp:426-500, 374-386`

---

## 9. Error Estimation Extension (Future)

The Zienkiewicz-Zhu (ZZ) error estimator extends SPR:

**Estimate local error:**
```
η_e = ||σ_h - σ*||_{L²(Ω_e)}    (element-wise)
```

**Global error estimate:**
```
||∇u - ∇u_h|| ≈ (Σ_e η_e²)^(1/2)
```

This is used for **adaptive mesh refinement** (h-adaptivity).

**Implementation status in BELFEM:** Not currently implemented, but infrastructure is in place:
- Recovered fields are already computed
- Element-wise fields are available (`compute_element_data`)
- Missing: norm computation and refinement logic

**Reference:** Zienkiewicz Vol. 1, §15.6

---

## 10. Implementation Clarification: SPR vs L² Projection

**Terminology note:** The BELFEM implementation is **true SPR** (Zienkiewicz-Zhu), not simple L² projection.

### What BELFEM Does (SPR):
```
1. Sample fields at superconvergent Gauss integration points
2. Fit polynomial to patch samples via weighted least-squares
3. Solve local SPD system: V·c = b (Cholesky)
4. Extrapolate polynomial to node location
```

### What BELFEM Does NOT Do:
```
❌ Global L² projection: ∫ N·σ* dΩ = ∫ N·σ_h dΩ (would require global solve)
❌ Simple nodal averaging: σ* = average of adjacent element values
❌ Moving least squares (MLS) with distance weighting (except in filter mode)
```

### Code Evidence:
- **Patch construction:** `cl_MaxwellPostprocessor.cpp:641-691`
- **Gauss-point sampling:** `cl_MaxwellPostprocessor.cpp:713-751`
- **Local solve:** `cl_MaxwellPostprocessor.cpp:755` (`posv(mVandermonde, mCoefficients)`)
- **Polynomial extrapolation:** `cl_MaxwellPostprocessor.cpp:757-768`

**This IS the classic SPR algorithm from Zienkiewicz & Zhu (1992).**

---

## 11. What BELFEM Recovery Deliberately Does NOT Do

**By design, the following are NOT implemented:**

### 1. Global L² Projection
- ❌ No global mass matrix assembly
- ❌ No global linear solve for recovery
- **Why:** Local SPR is faster and scales better to large meshes

### 2. Feedback into Solver
- ❌ Recovered fields do NOT influence next solve step
- ❌ No error-driven adaptive remeshing (yet)
- **Why:** Recovery is postprocessing only; avoids coupling issues

### 3. Alternative Recovery Methods
- ❌ No gradient averaging (simple nodal averaging)
- ❌ No equilibrated recovery (REP method)
- ❌ No moving least squares (except filter mode)
- **Why:** SPR provides optimal accuracy/cost tradeoff for BELFEM use cases

### 4. Stability Enhancements
- ❌ No SVD fallback for ill-conditioned patches
- ❌ No adaptive polynomial order reduction
- ❌ No condition number checking
- **Why:** Volume-weighting typically ensures well-conditioned systems

### 5. Advanced Features (Not Yet Implemented)
- ❌ No ZZ error estimator (infrastructure exists, not active)
- ❌ No higher-order recovery (p+2, p+3 polynomials)
- ❌ No recovery on curved manifolds
- **Why:** Current implementation satisfies production needs

**Recovery depth > 0:** Use only if Gramian is ill-conditioned at depth 0 (rare). Increases cost exponentially.

**If you need these features, consult development team before implementing.**

---

## 12. Summary

The Maxwell postprocessor implements a **mathematically rigorous** Zienkiewicz-Zhu Superconvergent Patch Recovery:

✅ **Exploits superconvergence** at Gauss integration points
✅ **Polynomial fitting** via weighted least-squares
✅ **Volume-weighted** Gramian ensures stability
✅ **Cholesky solver** for SPD systems (optimal)
✅ **Domain-aware** field computation (H, B, J, Jc)
✅ **MPI-parallel** with local patch assembly

**Convergence rate:** O(h^(p+1)) for recovered gradients/curls (one order better than raw FEM)

**Key references:**
- **Theory:** Zienkiewicz & Zhu (1992), Int. J. Numer. Methods Eng., 33(7)
- **Textbook:** Zienkiewicz & Taylor, Vol. 1, Chapter 15
- **Implementation:** `cl_FEM_Postprocessor.cpp`, `cl_MaxwellPostprocessor.cpp`

For **potential GPU acceleration** discussion (extreme-scale only), see: `../../../../todo/maxwell_postprocessor_gpu_acceleration.md`

---

## 13. Planned Optimization: Element-Wise Precomputed Weight Matrices

> **Status:** Planned. See `todo/postprocessor_recovery_findings.md` for full
> analysis, memory estimates, cache invalidation rules, and validation gates.

### 13.1 Motivation

All three current recovery modes (Sections 3.1) share the same fundamental
inefficiency: geometry-dependent work (Vandermonde assembly, Cholesky
factorization, coordinate transforms, polynomial evaluation) is repeated
every timestep, even though it depends only on static mesh geometry.

Additionally, per-node recovery (Modes 1 and 3) visits each element multiple
times — once per adjacent node — multiplying the element loop cost by ~6 in 2D.

### 13.2 Key Insight

The SPR result at node `n` can be decomposed as:

```
z_n = p(x_n)ᵀ · V_n⁻¹ · b_n                                              (*)
```

Let `q_n = V_n⁻¹ · p(x_n)` (computed via Cholesky factor/solve, not explicit
inversion). Since `q_n` depends only on geometry, it is precomputable.

For **linear** fields (H, J), where `y_{e,ip} = M_e(x_ip) · DOFs_e`:

```
z_n = Σ_e  [Σ_ip  w_e · dot(q_n, p(x_{e,ip})) · M_e(x_ip)] · DOFs_e
    = Σ_e  W_e[n] · DOFs_e
```

where `W_e[n]` is a small precomputed matrix per element per node:

```
W_e[n] = Σ_ip  w_e · dot(q_n, p(x_{e,ip})) · M_e(x_ip)
```

### 13.3 Element-Wise Storage

Store one matrix `W_e` per element of size `(n_nodes × n_fields) × n_dofs`:

```cpp
// Precomputed during initialize() — one-time cost
Cell< Matrix< real > > mWeightMatrices ;  // one W_e per element

// Runtime — single element loop, each element visited exactly once
for each element e:
    gather DOFs_e
    z_local = W_e * DOFs_e           // small matrix-vector product
    scatter-add z_local to node fields
```

**Advantages over current implementation:**
- Each element visited exactly **once** (vs ~6x in per-node modes)
- **Zero** runtime linear algebra (no Vandermonde, no Cholesky, no poly eval)
- Zero index overhead (DOFs and nodes implicit from element topology)
- Cache-friendly contiguous matrix per element
- Natural FEM scatter-add pattern (identical to stiffness assembly)

### 13.4 Nonlinear Fields

For nonlinear fields (B, J/Jc), DOF-based weights cannot be used because
the field depends nonlinearly on DOFs (see Section 4.1.1).

**Chosen approach:** Use precomputed scalar alpha weights at the integration-point
level for nonlinear fields, preserving numerical equivalence with the current
method. See `todo/postprocessor_recovery_findings.md` for details and acceptance
criteria.

### 13.5 Applicability

| Recovery class | Current mode | Optimization applies? |
|---------------|-------------|----------------------|
| `MaxwellPostprocessor` | Per-node (Mode 1) | Yes — highest impact |
| `Postprocessor` base | Batch (Mode 2) | Yes |
| `Gradient` | Per-node (Mode 3) | Yes (M = J⁻¹ · dNdξ) |

---

## Revision History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.3 | 2026-03-06 | C. Messe | Fixed Mode 2 description (was overstating advantages). Added Mode 3 (Gradient). Added Section 4.1.1 (linearity analysis). Added Section 13 (element-wise W_e optimization plan). |
| 1.2 | 2026-01-20 | C. Messe | Added Section 4.4: T-Matrix explanation (constraint matrices for hanging nodes) |
| 1.1 | 2026-01-20 | C. Messe | Added cross-references, glossary, SPR clarification, "What NOT to do" section |
| 1.0 | 2026-01-20 | C. Messe (via Claude Code) | Initial comprehensive documentation |
