# Slave-Side Bubble Enrichment for h-phi Interfaces

**Date:** 2026-02-25
**Purpose:** Analysis of what it takes to enrich both master and slave sides of thin shell sidesets
**Module:** fem/kernel, fem/maxwell, fem/interpolation

## Background

Bubble enrichment stabilizes the phi-phi interface (Dular et al. 2021, paper0) by adding
internal bubble functions to the element touching the interface, then statically condensing
them out. Currently, only the **master** element is enriched. This document describes how
to extend enrichment to the **slave** side.

### Current DomainType Priority Bug

The enum is `Air=1, Ferro=2, Coil=3, Conductor=4`. The swap condition in
`fix_facet_masters()` is `tM < tS`, making **higher values master**:

```
Actual master priority:   Conductor(4) > Coil(3) > Ferro(2) > Air(1)
Comment says:             Conductor > Air > Ferro > Coil
Intended (for enrichment): Air should be master (enrichment on master side)
```

Air is always slave in the current code. Since enrichment has been disabled
(`mEnrichSideSets = false`), this hasn't caused problems. Enabling slave-side enrichment
removes this constraint entirely — the priority ordering becomes a convention rather
than a functional requirement.

## Architecture of the Current Master-Side Enrichment

Enrichment involves three layers of pre-computed data and two types of integrals.

### Data Layer 1: Block Enrichment (Volume Integration Points)

**File:** `cl_FEM_Block.cpp:130-151`

```
Block::mEnrichmentData[f]  — one per facet of the element type
```

- Created by `Block::initialize_lookup_tables()`
- Bubble function: `create_bubble_function(elementType, f)`
- Integration points: standard VOLUME Gauss points via `populate(order, scheme)`
- **No orientation dependence** — volume points are the same regardless of which face
  is the interface

| Element | Facets | Entries |
|---------|--------|---------|
| TRI3    | 3 edges | 3 |
| TRI6    | 3 edges | 3 |
| TET4    | 4 faces | 4 |
| TET10   | 4 faces | 4 |

**Both master and slave blocks already create this data.** No changes needed.

### Data Layer 2: SideSet Enrichment (Facet Integration Points, Master Side)

**File:** `cl_FEM_SideSet.cpp:485-503`

```
SideSet::mEnrichmentData[f]  — one per facet of the master element type
```

- Bubble function: `create_bubble_function(mMasterType, f)`
- Integration points: facet points mapped to master coordinates via
  `populate_for_master(f, order, scheme)`
- Indexed the same as `mMasterIntegration` (by master facet index)

### Data Layer 3: SideSet Slave Enrichment (DOES NOT EXIST)

This is what we need to add. Indexed by `(facet, orientation)` permutation,
matching `mSlaveIntegration`.

### Calculator Wiring

**File:** `cl_FEM_Calculator.cpp:902-909`

On `link()`, the calculator sets up three enrichment pointers:

```cpp
mMasterVolumeIntegration = tMasterBlock->integration();        // Block standard integration
mVolumeEnrichment        = tMasterBlock->enrichment_data(f);   // Block bubble at volume points
mSideSetEnrichment       = mGroup->enrichment_data(f);         // SideSet bubble at facet points
```

Accessors: `volume_integration()`, `volume_enrichment()`, `sideset_enrichment()`

### Where Enrichment Is Used

| File | Integral Type | What It Computes |
|------|--------------|------------------|
| `mt_maxwell_interface.cpp` | Volume (over master element) | M = U W^{-1} V (static condensation of bubbles) |
| `mt_maxwell_aphi.cpp` | Surface (facet) + Volume | Surface: coupling with slave Ns via bubble gradient. Volume: same U W^{-1} V pattern |
| `mt_maxwell_phi_phi.cpp` | Volume | Same pattern (currently has `exit(0)`, unfinished) |

## The U-W-V Static Condensation Technique

### Mathematical Derivation

Consider an element adjacent to the interface with standard shape functions N_i
(i = 1..nn) and enrichment shape functions F_c (c = 1..nf). The enriched
approximation is:

```
phi = N_i * phi_i + F_c * psi_c
```

where phi_i are the standard nodal DOFs (global) and psi_c are the enrichment DOFs
(element-local). In the current implementation, the enrichment functions are
products F_c = N_i * E_j (PUM-style), where E_j are bubble functions. The index c
runs over all (i, j) pairs, so nf = nn * ne where ne is the number of enrichment
functions E_j.

The element stiffness contribution from this element involves the bilinear form
a(phi, phi') = integral( grad(phi) . mu . grad(phi') dV ). With the enriched
approximation, the element system is:

```
[ K_NN   K_NF ] [ phi ] = [ f_N ]
[ K_FN   K_FF ] [ psi ]   [ f_F ]
```

where:
- K_NN = integral( trans(B1) * B1 dV )  — standard stiffness (nn × nn)
- K_NF = integral( trans(B1) * B2 dV )  — coupling (nn × nf)
- K_FN = integral( trans(B2) * B1 dV )  — coupling (nf × nn) = trans(K_NF) for symmetric problems
- K_FF = integral( trans(B2) * B2 dV )  — enrichment self-coupling (nf × nf)

Here B1 = inv(J) * dN/dxi is the physical gradient of the standard shape functions,
and B2 = inv(J) * H is the physical gradient of the enrichment functions.

Since the enrichment DOFs psi are internal to the element (they don't connect to
any neighboring element), the enrichment load f_F = 0. From the second row:

```
psi = -inv(K_FF) * K_FN * phi
```

Substituting into the first row gives the condensed system:

```
( K_NN - K_NF * inv(K_FF) * K_FN ) * phi = f_N
```

The enrichment correction to the element stiffness is the **Schur complement**:

```
Delta_K = -K_NF * inv(K_FF) * K_FN = -U * inv(W) * V
```

with the notation used in the code:
- **U** = K_NF = integral( trans(B1) * B2 dV )   (nn × nf)
- **V** = K_FN = integral( trans(B2) * B1 dV )   (nf × nn)
- **W** = K_FF = integral( trans(B2) * B2 dV )   (nf × nf)

### The H Matrix (Product Rule in Parametric Space)

The enrichment functions F_c = N_i * E_j require the product rule for their
parametric derivatives:

```
dF_c / dxi_d = (dN_i / dxi_d) * E_j + N_i * (dE_j / dxi_d)
```

This is assembled into the matrix H (nd × nf):

```
H(d, c) = dNdXi1(d, i) * E(0, j) + N1(0, i) * dNdXi2(d, j)
```

where c = i * ne + j (columns indexed by node i, then enrichment function j).

The physical gradient is then B2 = inv(J) * H.

### Implementation in `mt_maxwell_interface.cpp`

```cpp
// Parametric derivatives of standard and enrichment shape functions
const Matrix<real> & dNdXi1 = tCalc->integration()->dNdXi(k);     // (nd × nn)
const Matrix<real> & dNdXi2 = aCalc->volume_enrichment()->dNdXi(k); // (nd × ne)
const Matrix<real> & N1     = tCalc->integration()->N(k);           // (1 × nn)
const Matrix<real> & E      = aCalc->volume_enrichment()->N(k);     // (1 × ne)

// Product rule: H(d, c) = dN_i/dxi_d * E_j + N_i * dE_j/dxi_d
for (uint d = 0; d < nd; ++d)
{
    uint c = 0;
    for (uint i = 0; i < nn; ++i)
        for (uint j = 0; j < ne; ++j)
            H(d, c++) = dNdXi1(d, i) * E(0, j)
                      + N1(0, i)     * dNdXi2(d, j);
}

// Physical gradients
const Matrix<real> & B1 = tCalc->B(k);      // inv(J) * dNdXi1
B2 = tCalc->invJ(k) * H;                     // inv(J) * H

// Accumulate over integration points
U += w(k) * trans(B1) * B2 * dV(k);    // (nn × nf)
V += w(k) * trans(B2) * B1 * dV(k);    // (nf × nn)
W += w(k) * trans(B2) * B2 * dV(k);    // (nf × nf)
```

After the integration loop, the enrichment correction is:

```cpp
aMatrices->M()  = U * inv(W) * V;
aMatrices->M() *= -constant::mu0;
```

The factor -mu0 comes from the weak form of the magnetostatic problem
(integral of mu0 * grad(phi) . grad(phi') dV with the sign convention used
in the h-phi formulation).

### Note on Dead Code

Line 83 of `mt_maxwell_interface.cpp` accumulates the standard stiffness
`trans(B1)*B1` into `aMatrices->M()`, but line 90 overwrites M entirely with
`U*inv(W)*V`. The standard stiffness is assembled by a separate function — this
line is dead code left from development.

### Sizes for Current and Future Implementations

| Enrichment | ne | nf = nn × ne | W size | inv(W) cost |
|------------|----|----|--------|-------------|
| Current face bubble (TET4) | 1 | 4 | 4×4 | Cheap (4×4 inverse) |
| PUM with 3 edge bubbles (TET4) | 3 | 12 | 12×12 | Still cheap |
| Current face bubble (TRI3) | 1 | 3 | 3×3 | Trivial |
| PUM with 2 edge bubbles (TRI3) | 2 | 6 | 6×6 | Trivial |

In all cases, W is a small dense matrix. The cost of inv(W) is negligible
compared to the global sparse solve.

## The Bubble Functions

### TRI3 Bubbles (2D)

Each bubble is quadratic and vanishes on its corresponding edge:

| Edge | File | Formula | Vanishes on |
|------|------|---------|-------------|
| 0 | `cl_IFG_TRI3A.hpp` | 4 xi eta | edge 0 (eta=0) |
| 1 | `cl_IFG_TRI3B.hpp` | 4 eta (1-xi-eta) | edge 1 (1-xi-eta=0) |
| 2 | `cl_IFG_TRI3C.hpp` | 4 xi (1-xi-eta) | edge 2 (xi=0) |

### TET4 Bubbles (3D)

Each bubble is cubic and vanishes on its corresponding face:

| Face | File | Formula | Zero coordinate |
|------|------|---------|-----------------|
| 0 | `cl_IFG_TET4A.hpp` | 27 xi eta (1-xi-eta) | zeta=0 (dN/dzeta=0) |
| 1 | `cl_IFG_TET4B.hpp` | 27 eta zeta (1-eta-zeta) | xi=0 (dN/dxi=0) |
| 2 | `cl_IFG_TET4C.hpp` | 27 xi zeta (1-xi-zeta) | eta=0 (dN/deta=0) |
| 3 | `cl_IFG_TET4D.hpp` | 27 xi eta zeta | 1-xi-eta-zeta=0 |

**Key property:** On the interface face, N=0 always but grad(N) is nonzero in the
normal direction. It is the gradient that provides the stabilization.

**Critical finding:** This property is actually the root of the problem — see next section.

## Why Current Face-Vanishing Bubbles Are Insufficient

### The Inf-Sup Stabilization Requirement

Dular et al. 2021 (paper0) showed that the h-phi interface requires **hierarchical
enrichment** to satisfy the inf-sup condition. The enrichment must increase the
polynomial order of the **trace space** on the interface — i.e., the restriction of
the finite element space to the interface face.

### What Dular Actually Did

Dular's enrichment uses **edge bubble functions that are nonzero on the interface**.
For a triangular interface face with three edges, the enrichment functions are
quadratic edge bubbles that:

- Are nonzero along their respective edge of the interface face
- Vanish at vertices (preserving nodal compatibility)
- Enrich the trace space from piecewise-linear to piecewise-quadratic on the interface

This is exactly what TRI6 mid-edge shape functions provide: they are nonzero on
their edge and vanish at endpoints. On a TET10, the mid-edge nodes on the interface
face are precisely the Dular-type enrichment functions.

### What BELFEM Currently Does (and Why It's Wrong)

BELFEM's bubble functions (`BubbleFace0`–`BubbleFace3` for TET4) **vanish on their
own face**. When BubbleFace0 is used to enrich an element at face 0:

```
BubbleFace0 = 27 * xi * eta * (1 - xi - eta)
On face 0 (zeta = 0):  N = 27 * xi * eta * (1 - xi - eta) — this is nonzero!
```

Wait — BubbleFace0 is actually a face bubble for face 0 that is nonzero on face 0
but vanishes on the other three faces. However, the function chosen for enrichment
at face `f` is `BubbleFace_f`, which is the bubble associated with face `f`. Looking
at the code more carefully:

| Face Index | Bubble | Formula | On its own face |
|------------|--------|---------|-----------------|
| 0 | TET4A | 27 xi eta (1-xi-eta-zeta) | zeta=0 → 27 xi eta (1-xi-eta) ≠ 0 |
| 1 | TET4B | 27 eta zeta (1-eta-zeta) | xi=0 → 27 eta zeta (1-eta-zeta) ≠ 0 |
| 2 | TET4C | 27 xi zeta (1-xi-zeta) | eta=0 → 27 xi zeta (1-xi-zeta) ≠ 0 |
| 3 | TET4D | 27 xi eta zeta | 1-xi-eta-zeta=0 → 27 xi eta zeta ≠ 0 |

Correction: BubbleFace_f is nonzero on face f. The nomenclature in the earlier
section was misleading — each bubble does NOT vanish on its own face. Rather, each
`BubbleFace_f` is the product of the three barycentric coordinates of the vertices
of face f, and it IS nonzero on face f.

However, the bubble function has a critical limitation: it is a **single** cubic
function on the face. The trace on the interface is:

```
BubbleFace0|_{face 0} = 27 * xi * eta * (1 - xi - eta)
```

This is a single cubic bubble — it enriches the trace space by exactly one function
that vanishes at all vertices and all edges (it's an interior face mode). It does
NOT enrich the **edge** trace space. The inf-sup condition requires enrichment of the
edge traces, not just the face interior.

### The Gap Between Face Bubbles and Edge Bubbles

| Property | Face Bubble (BELFEM) | Edge Bubbles (Dular) |
|----------|---------------------|---------------------|
| Number of functions | 1 per face | 3 per face (one per edge) |
| Nonzero on | Face interior only | Along edges |
| Vanishes at | All edges and vertices | Opposite vertices only |
| Trace enrichment | p=3 interior mode | p=2 edge modes |
| Inf-sup | Insufficient | Sufficient |

The face bubble is an interior mode that doesn't improve the approximation along the
edges of the interface, where the coupling between adjacent interface elements occurs.
Dular's edge bubbles enrich precisely these critical edge traces.

## Correct Approach: PUM-Style Enrichment with TET10 Mid-Edge Functions

### Formulation

Instead of adding a single face bubble, use a **Partition of Unity Method (PUM)**
style enrichment. The enriched approximation on an element adjacent to the interface
is:

```
phi = N_i * phi_i + N_i * E_j * psi_j
```

where:
- `N_i` are the standard TET4 shape functions (i = 0..3)
- `phi_i` are the standard nodal DOFs (global)
- `E_j` are **enrichment functions** — the three TET10 mid-edge shape functions
  for the edges of the interface face
- `psi_j` are enrichment DOFs (element-local, to be condensed)

For interface face 0 (nodes 0, 1, 2), the enrichment functions E_j are the TET10
mid-edge functions for edges (0-1), (1-2), (0-2). Each E_j:

- Is quadratic: E_j = 4 * L_a * L_b for edge (a, b)
- Equals 1 at the mid-edge point, 0 at all vertices
- Is nonzero along its edge on the interface face
- Provides exactly the Dular-type trace enrichment

### Why N_i * E_j Works

The product `N_i * E_j` has the following properties:

1. **Vanishes on the face opposite node i.** Since N_i = 0 on the face opposite to
   node i, the product N_i * E_j = 0 there, regardless of E_j.

2. **May be nonzero on other non-interface faces.** For a node i that is ON the
   interface face, N_i * E_j can be nonzero on the two other faces adjacent to
   node i (besides the interface face). This is a "leak" of the enrichment onto
   non-interface faces.

3. **The leak doesn't matter.** The enrichment DOFs psi_j are **statically condensed**
   at the element level. The global solution is:

   ```
   phi_global = N_i * phi_i     (always linear, always continuous)
   ```

   The enrichment terms N_i * E_j * psi_j only live inside the element computation.
   They modify the element stiffness matrix through condensation (M = U * W^{-1} * V)
   but never appear in the global system. Since the enrichment never propagates to
   neighboring elements through shared DOFs, the "incompatibility" on non-interface
   faces is invisible to the global solution.

4. **The condensation captures the interface physics.** The coupling matrices U and V
   in the static condensation M = U * W^{-1} * V involve integrals of the form:

   ```
   U_ij = integral( grad(N_i) . mu . grad(N_i * E_j) ) dV     (volume)
   V_ij = integral( N_s . n x grad(N_i * E_j) ) dS             (surface)
   ```

   Because E_j is nonzero on the interface (unlike the face bubble), the surface
   coupling terms V_ij are nonzero. This is what provides the inf-sup stabilization:
   the condensed enrichment modifies the element stiffness in a way that accounts for
   the higher-order interface behavior.

### Comparison of Approaches

| Aspect | Face Bubble (current) | PUM with TET10 mid-edge |
|--------|----------------------|------------------------|
| Enrichment DOFs per element | 1 | 3 per interface-adjacent node × nodes on face |
| Global DOFs added | 0 (condensed) | 0 (condensed) |
| Trace enrichment | Interior mode (insufficient) | Edge modes (sufficient) |
| Compatibility | Trivially compatible | Formally incompatible, but condensation makes it irrelevant |
| T-matrix impact | None | None (element-local DOFs) |
| Static condensation size | 1×1 (scalar division) | Small dense system (still cheap) |

### No Impact on T-Matrices for Hanging Edges

The T-matrices for hanging edges operate on **global DOFs**. They express hanging
node DOFs as linear combinations of master node DOFs:

```
phi_hanging = sum_j ( w_j * phi_master_j )
```

The enrichment DOFs psi_j are element-local and never enter the global DOF numbering.
They are eliminated during element matrix assembly via static condensation, before
the global system is ever formed. Therefore:

- T-matrices remain unchanged
- `DofData::create_dofwise_t_matrices_master()` is unaffected
- `consolidate_dofs()` and `relink_dofs()` are unaffected
- The cascade handling for edges on multiple interfaces is unaffected

### Neighbor Enrichment

For thin shell tapes, the interface separates two **Lagrange-based elements** (both
TET4 or both TET10). Since both elements can support the PUM-style enrichment, both
the master and slave side can be enriched independently.

This is particularly natural for thin shell configurations:

```
     Block A (Air)          Thin Shell          Block B (Air/Conductor)
  ┌─────────────────┐    ┌───┐    ┌─────────────────┐
  │                 │    │   │    │                 │
  │  TET4 + enrich  │────│   │────│  TET4 + enrich  │
  │                 │    │   │    │                 │
  └─────────────────┘    └───┘    └─────────────────┘
       N_i*E_j*psi_j                  N_i*E_j*psi_j
       (condensed)                    (condensed)
```

Each side is enriched independently with its own set of condensed enrichment DOFs.
The enrichment on one side does not need to "match" the enrichment on the other
side because:

1. The enrichment DOFs are condensed out at the element level
2. The global solution (N_i * phi_i) remains linear and continuous across the interface
3. Each side's enrichment captures its own higher-order behavior through condensation

For non-thin-shell interfaces (e.g., two volume blocks meeting at a face), the same
principle applies: any Lagrange-based element adjacent to the interface can be
enriched, regardless of what's on the other side.

## Orientation Handling

### The Magic Orientation Matrix (TET4)

**File:** `fn_IF_initialize_integration_points_on_facet.cpp:576-580`

For slave-side facet integration, `intpoints_tet(slaveIndex, orientation, ...)` maps
2D triangle integration points to 3D tet coordinates via:

```
tCase = slaveIndex * 3 + orientation     (4 faces x 3 orientations = 12 cases)

tIndex = { { 0, 2, 3, 2, 1, 3, 0, 3, 1, 0, 1, 2 },   // triangle coord 0 -> tet row
           { 3, 0, 2, 3, 2, 1, 1, 0, 3, 2, 0, 1 },   // triangle coord 1 -> tet row
           { 2, 3, 0, 1, 3, 2, 3, 1, 0, 1, 2, 0 },   // triangle coord 2 -> tet row
           { 1, 1, 1, 0, 0, 0, 2, 2, 2, 3, 3, 3 } };  // tet row that stays 0 (on-face)

aPoints.set_row( tIndex(0, tCase), tPoints.row(0) );   // map triangle xi
aPoints.set_row( tIndex(1, tCase), tPoints.row(1) );   // map triangle eta
aPoints.set_row( tIndex(2, tCase), tPoints.row(2) );   // map triangle zeta
// row tIndex(3, tCase) stays 0  -- this is the face-defining coordinate
```

Row 3 of `tIndex` encodes which barycentric coordinate is zero:

- Cases 0-2 (face 0): coord 1 (eta) = 0
- Cases 3-5 (face 1): coord 0 (xi) = 0
- Cases 6-8 (face 2): coord 2 (zeta) = 0
- Cases 9-11 (face 3): coord 3 (1-xi-eta-zeta) = 0

### How 12 Cases Collapse

The 3 orientations per face permute which of the 3 remaining barycentric coordinates
receive the 2D triangle integration coordinates. The bubble function itself depends
only on the face index (which coordinate is zero), not on the orientation.

**For the volume integral:** Integration points are standard tet volume Gauss points,
independent of which face is the interface. Block enrichment data already handles this
with just 4 entries (one per face). **12 cases collapse to 4.**

**For the surface integral:** The integration points on the face DO depend on
orientation (different mappings from 2D triangle to 3D tet coordinates). However,
since the bubble function vanishes on its own face (N=0), only the gradient matters.
The non-zero gradient component is normal to the face:

```
BubbleFace0 on face 0 (eta=0):
  dN/dxi  = 27 * 0 * (...)  = 0       (tangential)
  dN/deta = 27 * xi * (1-xi)          (NORMAL - nonzero)
  dN/dzeta = 0                         (tangential)
```

For different orientations, `xi` in the formula above gets mapped to different
triangle coordinates (triXi, triEta, triZeta), giving:

```
  Case 0: dN/deta = 27 * triXi * (1 - triXi)
  Case 1: dN/deta = 27 * triEta * (1 - triEta)
  Case 2: dN/deta = 27 * triZeta * (1 - triZeta)
```

These are different functions evaluated at different integration point positions.
In isolation (if the integrand only involved the bubble), the symmetry of triangle
Gauss quadrature would make them equivalent. But in practice, the bubble gradient
is multiplied by the standard shape functions (N_s, B_m), which are NOT symmetric
under orientation permutations. **So all 12 cases are distinct in general.**

### Permutation Counts by Element Type

| Element | Faces | Orientations/face | Total permutations | Source |
|---------|-------|-------------------|--------------------|--------|
| TRI/QUAD (2D) | n | 1 (flip only) | n | `number_of_orientations()` returns 1 |
| TET | 4 | 3 | 12 | `number_of_orientations()` returns 3 |
| HEX | 6 | 4 | 24 | `number_of_orientations()` returns 4 |
| PENTA | 3 quad + 2 tri | 4 / 3 | 18 | Mixed |
| PYRA | 4 tri + 1 quad | 3 / 4 | 16 | Mixed |

Currently, only TRI and TET elements have bubble functions implemented
(`create_bubble_function` in `cl_IF_InterpolationFunctionFactory.cpp:274`).

## What Needs to Change

### 0. Replace Bubble Functions with PUM Enrichment Functions

**This is the most critical change.** The current `create_bubble_function()` returns
face bubbles (BubbleFace0–3). These need to be replaced (or supplemented) with
TET10 mid-edge functions for the edges on the interface face.

For interface face `f` with vertex nodes (a, b, c), the enrichment functions are:

```
E_0 = 4 * L_a * L_b     (mid-edge function for edge a-b)
E_1 = 4 * L_b * L_c     (mid-edge function for edge b-c)
E_2 = 4 * L_a * L_c     (mid-edge function for edge a-c)
```

The PUM products `N_i * E_j` (for each node i on the interface face and each
enrichment function j) form the actual enrichment basis. The number of enrichment
DOFs per element is up to `num_interface_nodes * num_edge_functions` (e.g., 3 × 3 = 9
for TET4 with one interface face), though some products may be linearly dependent
or zero, reducing the effective count.

**Note:** The exact bubble function formulas will be generated with a MATLAB program
to ensure correctness.

### 1. SideSet: Add Slave Enrichment Data

**File:** `cl_FEM_Group.hpp` (base class for SideSet and Block)

Add new member alongside `mEnrichmentData`:

```cpp
Cell< IntegrationData * > mSlaveEnrichmentData ;
```

Add accessor:

```cpp
const IntegrationData * slave_enrichment_data( const uint aPermutation ) const ;
```

**File:** `cl_FEM_SideSet.cpp:506-529` (inside the `if(tHaveSlave)` block)

After creating `mSlaveIntegration`, add enrichment creation:

```cpp
if ( mParent != nullptr && mParent->iwg()->enrich_sidesets() )
{
    InterpolationFunctionFactory tFactory ;
    uint tNumFacets = mesh::number_of_facets( mSlaveType ) ;

    // Clear old data
    for ( IntegrationData * tData : mSlaveEnrichmentData )
    {
        delete tData ;
    }
    mSlaveEnrichmentData.set_size( tNumPermutations, nullptr ) ;

    uint tCount = 0 ;
    for ( uint f = 0 ; f < tNumFacets ; ++f )
    {
        InterpolationFunction * tBubble = tFactory.create_bubble_function( mSlaveType, f ) ;

        for ( uint o = 0 ; o < mesh::number_of_orientations( mSlaveType, f ) ; ++o )
        {
            // The first IntegrationData for each face owns the bubble function;
            // subsequent orientations of the same face share it (own = false).
            bool tOwn = ( o == 0 ) ;
            mSlaveEnrichmentData( tCount ) = new IntegrationData( mSlaveType, tBubble, tOwn ) ;
            mSlaveEnrichmentData( tCount )->populate_for_slave( f, o, tIntegrationOrder, tScheme ) ;
            ++tCount ;

            // For subsequent orientations, create a new bubble (simpler than sharing)
            if ( o < mesh::number_of_orientations( mSlaveType, f ) - 1 )
            {
                tBubble = tFactory.create_bubble_function( mSlaveType, f ) ;
            }
        }
    }
}
```

**Note on bubble ownership:** Each `IntegrationData` can own its shape function
(`mOwnShapeFunction = true` in constructor). Since `populate_for_slave()` evaluates
and caches N and dNdXi at the integration points, the simplest approach is to
create one bubble per permutation (each IntegrationData owns its own copy). The
bubble objects are tiny (no state beyond the template parameters), so the memory
cost is negligible.

### 2. Calculator: Add Slave Enrichment Members

**File:** `cl_FEM_Calculator.hpp`

Add alongside existing enrichment members (line 150):

```cpp
// for slave-side enrichment
const IntegrationData * mSlaveVolumeIntegration = nullptr ;
const IntegrationData * mSlaveVolumeEnrichment  = nullptr ;
const IntegrationData * mSlaveSetEnrichment     = nullptr ;
```

Add accessors:

```cpp
const IntegrationData * slave_volume_integration() const { return mSlaveVolumeIntegration ; }
const IntegrationData * slave_volume_enrichment() const { return mSlaveVolumeEnrichment ; }
const IntegrationData * slave_sideset_enrichment() const { return mSlaveSetEnrichment ; }
```

**File:** `cl_FEM_Calculator.cpp:912-921` (in `link()`, after slave section)

Add slave enrichment wiring:

```cpp
if ( aElement->slave() != nullptr && mGroup->parent()->iwg()->enrich_sidesets() )
{
    Block * tSlaveBlock = mGroup->parent()->block( aElement->slave()->element()->block_id() ) ;
    uint tSlaveIndex = aElement->facet()->index_on_slave() ;

    mSlaveVolumeIntegration = tSlaveBlock->integration() ;
    mSlaveVolumeEnrichment  = tSlaveBlock->enrichment_data( tSlaveIndex ) ;

    // Use the same function pointer logic as mSlaveIntegration
    // to compute the permutation index for slave enrichment
    mSlaveSetEnrichment = mGroup->slave_enrichment_data(
        /* permutation index matching mSlaveIntegration */ ) ;
}
```

The permutation index for the slave enrichment must match the one used for
`mSlaveIntegration`. The existing `slave_integration_2d/tet/hex()` functions
already compute this:

- 2D: `index_on_slave()`
- TET: `index_on_slave() * 3 + orientation_on_slave() - 1`
- HEX: `index_on_slave() * 4 + orientation_on_slave() - 1`

We can reuse the same function pointer pattern or compute the index directly.

### 3. Matrix Computation: Add Slave-Side Loop

**File:** `mt_maxwell_interface.cpp`

The existing code computes M = U_m W_m^{-1} V_m for the master element.
Add an analogous computation for the slave:

```cpp
// --- Slave-side bubble enrichment ---
Calculator * tCalcS = aCalc->group()->parent()->block(
    aCalc->element()->slave()->element()->block_id() )->calculator() ;
tCalcS->link( aCalc->element()->slave() ) ;

// Same structure as master:
// U_s, V_s, W_s computed over slave volume integration points
// M += U_s * inv(W_s) * V_s

for ( uint k = 0 ; k < tCalcS->num_intpoints() ; ++k )
{
    const Matrix< real > & dNdXi1 = tCalcS->integration()->dNdXi( k ) ;
    const Matrix< real > & dNdXi2 = aCalc->slave_volume_enrichment()->dNdXi( k ) ;
    const Matrix< real > & N1 = tCalcS->integration()->N( k ) ;
    const Matrix< real > & E  = aCalc->slave_volume_enrichment()->N( k ) ;
    // ... same H, B2, U, V, W assembly as master ...
}
```

**File:** `mt_maxwell_aphi.cpp`

The surface integral (lines 52-111) uses `sideset_enrichment()` for the master bubble
gradient at the interface. For the slave:

```cpp
const Matrix< real > & dNdXi2s = aCalc->slave_sideset_enrichment()->dNdXi( k ) ;
const Matrix< real > & Es      = aCalc->slave_sideset_enrichment()->N( k ) ;
```

Then assemble the slave-side coupling terms analogous to the master-side ones.

## Summary of Changes

| File | Change | Lines |
|------|--------|-------|
| `cl_FEM_Group.hpp` | Add `mSlaveEnrichmentData`, accessor, destructor cleanup | ~15 |
| `cl_FEM_SideSet.cpp` | Create slave enrichment data in `initialize_lookup_tables()` | ~25 |
| `cl_FEM_Calculator.hpp` | Add slave enrichment members + accessors | ~15 |
| `cl_FEM_Calculator.cpp` | Wire slave enrichment in `link()` | ~15 |
| `mt_maxwell_interface.cpp` | Add slave-side volume integral loop | ~35 |
| `mt_maxwell_aphi.cpp` | Add slave-side surface + volume integral | ~40 |
| **Total** | | **~145** |

## Implications

1. **DomainType priority becomes irrelevant for enrichment.** With both sides enriched,
   it doesn't matter which is master. The priority ordering only affects normal direction
   convention and thin shell node duplication consistency (separate issue).

2. **The `mEnrichSideSets` flag controls everything.** No new flags needed. When true,
   both master and slave get enriched. When false (default), no enrichment on either side.

3. **Performance impact is modest.** With PUM-style enrichment, the condensation
   involves a small dense system (up to 9×9 per element for TET4) instead of the
   current 1×1 scalar. This is still negligible compared to global solve time.
   Adding the slave side doubles the enrichment work but enrichment remains a small
   fraction of total assembly time.

4. **Memory:** For TET4 with slave enrichment, the SideSet stores 12 additional
   `IntegrationData` objects (4 faces x 3 orientations). Each is small (integration
   points + pre-evaluated N, dNdXi at those points). Negligible compared to mesh storage.

5. **Correctness of inf-sup stabilization.** The PUM-style enrichment with TET10
   mid-edge functions provides the same trace-space enrichment as Dular et al. 2021
   (paper0), satisfying the inf-sup condition. The current face-bubble approach does
   not enrich the edge traces and is insufficient for stabilization.

6. **No impact on hanging edge T-matrices.** Enrichment DOFs are condensed at the
   element level and never enter the global DOF system. T-matrices, `relink_dofs()`,
   and the cascade handling for edges on multiple interfaces are all unaffected.
