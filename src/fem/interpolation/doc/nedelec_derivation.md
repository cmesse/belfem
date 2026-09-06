# Nédélec Element Derivation {#fem_interpolation_nedelec_derivation}

**Date:** 2026-08-14
**Purpose:** The mathematical derivation behind BELFEM's Lagrange and Nédélec (edge) elements:
interpolation operators, barycentric coordinates, the geometry Jacobian, Whitney edge functions
for TRI3/TRI6/TET4/TET10, curl operators, and the edge/face generation concept.
**Module:** src/fem/interpolation
**Origin:** Distilled from Christian Messe's quasi-magnetostatic theory notes (LaTeX, written in
the Lagrange-multiplier era of the code). Transcription errata in the notes were fixed during
extraction. Each equation block names its TeX source for traceability.
**Scope:** This document covers the volume elements TRI3, TRI6, TET4 and TET10. The LINE3 and
HEX8 elements and the thin-shell family (QUAD4TS, PENTA6TS, HEX8TB) are documented in
[nedelec.md](nedelec.md) and [nedelec_thinshell.md](nedelec_thinshell.md). For the class
framework (EdgeFunction lifecycle, factory, DOF counts), see [nedelec.md](nedelec.md); for the
weak forms these operators serve, see `../../maxwell/doc/maxwell_weak_forms.md`.

---

## 1. Interpolation Functions and Operators

*(source: `introduction/notation.tex`)*

In finite element theory, a scalar field `phi(x)` is discretized by linearly combining an
interpolation function with basis values. These basis values represent the degrees of freedom
against which the system of equations is solved. They are often attributed to the values the
field assumes at the supporting nodes of the computational mesh. In this case, we collect the
nodal values in a column vector:

```
phi  ≈  N * phi_hat  =  [ N^1  ...  N^n ] * [ phi_hat^1, ..., phi_hat^n ]^T        (lagrange_fashion)
```

The superscript is the node number; the hat marks a value attributed to a node. The spatially
dependent *interpolation function* `N` is defined by the chosen element type and is also called
*shape function* or *ansatz function*. Elements that connect the basis to nodal values are
called *Lagrangian elements*. That is by far not the only possibility: *Hermitian elements*
(mechanical beams and plates) connect the basis to nodal values and their spatial derivatives,
and *isogeometric elements* use B-spline functions.

Within the context of the Maxwell equations, it is beneficial to associate the supporting basis
with the *edges* of an element. This way it can be guaranteed that Ampère's circuital law is
fulfilled. For the magnetic field `h`, scalar values `h_tilde^k` are associated with the edges,
and the edge-based interpolation function `E` contains the geometric information that translates
the scalar basis into a vector field:

```
h  ≈  E * h_tilde  =  [ E_x^1 ... E_x^m ]   [ h_tilde^1 ]
                      [ E_y^1 ... E_y^m ] * [    ...     ]
                      [ E_z^1 ... E_z^m ]   [ h_tilde^m ]
```

For obvious reasons, these elements are called *edge elements*; the names *Whitney elements*
and *Nédélec elements* (Nédélec 1980) are also common in the literature.

If a vector field `a(x)` is instead discretized with the Lagrangian approach, the interpolation
function carries no directional information; the vector components are interpolated as
independent scalar fields, with the degrees of freedom in node-wise order:

```
a  ≈  N * a_hat,   N = [ N^1  0    0   ...  N^n  0    0  ]      a_hat = [ a_x^1, a_y^1, a_z^1,
                       [ 0    N^1  0   ...  0    N^n  0  ]               ...,
                       [ 0    0    N^1 ...  0    0    N^n ]              a_x^n, a_y^n, a_z^n ]^T
```

The spatial gradient of the scalar field is obtained by differentiating the shape function,
which creates the `B`-operator (derivatives are indicated by a comma, so `N^i_,x` abbreviates
`dN^i/dx`):

```
h  =  grad phi  ≈  B * phi_hat,     B = [ N^1_,x ... N^n_,x ]
                                        [ N^1_,y ... N^n_,y ]                      (lagrange_b_operator)
                                        [ N^1_,z ... N^n_,z ]
```

In analogy to the gradient operator `B`, the curl operator `C` is introduced. For the
node-interpolated vector field `a`, the magnetic flux density is:

```
b  =  curl a  ≈  C * a_hat

     [  0       -N^1_,z   N^1_,y  |     |  0       -N^n_,z   N^n_,y ]
C  = [  N^1_,z   0       -N^1_,x  | ... |  N^n_,z   0       -N^n_,x ]
     [ -N^1_,y   N^1_,x   0       |     | -N^n_,y   N^n_,x   0      ]
```

(The notes carried sign slips in the first column block of this matrix; the form above is the
corrected one, consistent with `(curl a)_y = a_x,z - a_z,x` and `(curl a)_z = a_y,x - a_x,y`.)

For the edge-interpolated field `h`, the current density follows the same pattern with the
curl of the edge functions:

```
j  =  curl h  ≈  C * h_tilde,    C = [ E^1_z,y - E^1_y,z  ...  E^m_z,y - E^m_y,z ]
                                     [ E^1_x,z - E^1_z,x  ...  E^m_x,z - E^m_z,x ]     (curloperator)
                                     [ E^1_y,x - E^1_x,y  ...  E^m_y,x - E^m_x,y ]
```

### 1.1 Two-Dimensional Simplification

In 2D the operators are built the same way. The magnetic vector potential reduces to the scalar
`a_z`:

```
a ≡ a_z ≈ N * a_hat

b ≡ [ b_x ]  =  C * a_hat,    C = [  N^1_,y ...  N^n_,y ]
    [ b_y ]                       [ -N^1_,x ... -N^n_,x ]
```

For the non-conducting parts, interpolated with the scalar potential `phi`:

```
h ≡ [ h_x ]  ≈  B * phi_hat,     j ≡ j_z ≡ curl grad phi ≡ 0
    [ h_y ]
```

A sign-convention note: this rule is unsigned, as in the original notes. BELFEM's physical
convention in the phi-domain is `h = -grad phi` (Arsenault et al. 2023); the phi kernel
(`src/fem/maxwell/matrices/mt_maxwell_phi.cpp`) deliberately drops the minus because the mass
term `integral mu |h|^2` is even in the sign. Keep that in mind before wiring `B phi_hat` into
anything where the sign survives, such as an interface term.

And for the conducting parts, using the edge field:

```
h ≈ E * h_tilde,     j ≡ j_z ≈ C * h_tilde,    C = [ E^1_y,x - E^1_x,y  ...  E^m_y,x - E^m_x,y ]
```

### 1.2 Axisymmetric Simplification

In the axisymmetric case the magnetic *vector* potential exists only in the tangential
direction (the notes mislabeled it "scalar potential" here):

```
a ≡ a_t ≈ N * a_hat

b ≡ [ b_r ]  =  C * a_hat,    C = [ -N^1_,z         ...  -N^n_,z         ]
    [ b_z ]                       [ N^1/r + N^1_,r  ...  N^n/r + N^n_,r  ]
```

The `1/r` term comes from `b_z = (1/r) d(r a_t)/dr`. The `phi`-domain and conducting-domain
rules follow the 2D pattern with `(r, z)` in place of `(x, y)`; the current density is the
tangential component `j_t` in both cases.

---

## 2. Lagrange Triangle

*(source: `triangle/triangle_lagrange.tex`)*

### 2.1 Parameter Space

The parameter coordinates are the barycentric coordinates `(xi, eta, zeta)` with the unity
requirement:

```
xi + eta + zeta = 1                                                     (barycentric_unity)
```

### 2.2 First-Order Interpolation

Due to the unity requirement, the third coordinate can be expressed through the first two, so
the interpolation for a TRI3 reads:

```
phi(xi, eta)  ≈  [ xi   eta   1 - xi - eta ] * [ phi_hat^1, phi_hat^2, phi_hat^3 ]^T    (N_operator)
```

### 2.3 Geometry Jacobian and B-Operator

The `B`-operator must contain the spatial derivatives of `N` in the geometry space. Applying
the chain rule, the derivatives are evaluated in the parameter space and transformed:

```
B = J^-T * [ N^1_,xi   N^2_,xi   N^3_,xi  ]  =  [ grad xi | grad eta ] * [ 1  0  -1 ]      (bop)
           [ N^1_,eta  N^2_,eta  N^3_,eta ]                              [ 0  1  -1 ]
```

where `J` is the *geometry Jacobian*. In 2D:

```
J^T = [ x_,xi   y_,xi  ]  =  [ N_,xi  ] * [ x_hat  y_hat ]                    (geometryjacobian2d)
      [ x_,eta  y_,eta ]     [ N_,eta ]
```

The transposed symbol in this equation points out a common pitfall when implementing finite
elements: a Jacobian matrix is used, for example, when a Newton-Raphson iteration is performed.
Any math book that is not related to finite elements will provide the definition without the
transpose (Bronshtein). When formulating finite elements, however, the geometry Jacobian is
always used in its transposed form. For that reason, most relevant finite element books omit
the transposed symbol. This may improve readability, but can cause confusion and antagonizing
pain when implementing the equations. Reader's discretion is advised.

BELFEM stores exactly this transposed form: in `EF_TRI3::link()`
([cl_EF_TRI3.cpp](../nedelec/cl_EF_TRI3.cpp):39-42), `mJ` rows run over the parameters, which
is `J^T` in the math-book convention.

Both straight-edged triangles and tetrahedra have the convenient feature that the geometry
Jacobian remains constant within the element, regardless of the order of the interpolation
function. This is not the case for quadrilaterals or hexahedra. For the straight triangle:

```
J^T = [ x_hat^1 - x_hat^3   y_hat^1 - y_hat^3 ]                              (tri3_jacobian)
      [ x_hat^2 - x_hat^3   y_hat^2 - y_hat^3 ]
```

### 2.4 Second-Order Interpolation

The TRI6 adds three midside nodes. The construction of `N` and `B` works the same way as for
the linear element. If the triangle has straight edges, the geometry Jacobian collapses to the
constant form above; for curved edges, `J` depends on `xi` and `eta`.

---

## 3. Nédélec Triangle

*(source: `triangle/triangle_nedelec.tex`)*

Ampère's circuital law states that the integral of the magnetic field over an arbitrary closed
loop within the computational domain must equal the electric current in the enclosed area. This
cannot be unconditionally fulfilled if Lagrangian elements are used. Even worse, it cannot be
guaranteed that the computation converges to the physical solution. The remedy is a basis
connected to the edges rather than the nodes (Nédélec 1980).

### 3.1 First-Order Edge Interpolation

Although `h` is a vector field, the directional information is implicitly given by the
direction of the edge in the geometry space. The edge degree of freedom therefore reduces to a
scalar, which significantly lowers the number of unknowns. It also means that the edge function
depends on the actual shape of the triangle. With the local edge numbering `edge k` connecting
the barycentric pair `(lambda_a, lambda_b)` in cyclic order `(xi,eta), (eta,zeta), (zeta,xi)`,
the three edge vectors are the Whitney 1-forms:

```
E^1 = s_1 * ( xi   * grad eta  - eta  * grad xi   )                          (nedelec_theta1..3)
E^2 = s_2 * ( eta  * grad zeta - zeta * grad eta  )
E^3 = s_3 * ( zeta * grad xi   - xi   * grad zeta )
```

The spatial gradients of the parameter coordinates come from the inverse of the geometry
Jacobian, and the unity requirement gives:

```
grad xi + grad eta + grad zeta = 0
```

This is exactly what `EF_TRI3` implements: `precompute()` stores the barycentric pair per edge
([cl_EF_TRI3.cpp](../nedelec/cl_EF_TRI3.cpp):80-95) and `E()` assembles
`s_k (G_k grad(next) - H_k grad(prev))` (`:104-113`).

The scalar `s_k` indicates the orientation of the edge. To determine it, the topology of the
mesh must be known. Each edge is owned by one global direction; the sign is:

```
s_k = +1   if edge k runs in the same direction as its global (owner) direction
      -1   else
```

The signs are fetched in `link()` via `aElement->edge_directions(mS)`. With `s_k` baked in, the
basis has unit circulation on its own global edge: `int_edge_j E_k . dl = delta_jk` in global
orientation (equivalently `s_k delta_jk` in the element-local direction). This "DOF equals
circulation" convention is the load-bearing contract for DOF sharing between elements; see
[nedelec.md](nedelec.md) §5. Note that first-order Whitney forms satisfy it automatically, with
no extra normalization factor.

### 3.2 Curl Operator

The entries of the `C`-operator are found by deriving the edge functions in the parameter space
and transforming with the chain rule. On a straight-edged triangle, `E^i_x,x` and `E^i_y,y`
vanish, and after a few lines of algebra (the MATLAB fragment `matlab_curl/fragment.m` in the
notes repository reproduces it symbolically) the operator collapses to a constant:

```
C = 2 / det(J) * [ s_1  s_2  s_3 ]
```

Implemented verbatim in `EF_TRI3::link()` ([cl_EF_TRI3.cpp](../nedelec/cl_EF_TRI3.cpp):64-68).
The result is plausible: each Whitney form `lambda_a grad lambda_b - lambda_b grad lambda_a`
has `curl = 2 grad lambda_a x grad lambda_b`, which is constant on an affine element.

### 3.3 Second-Order Interpolation (TRI6)

There is a family of at least four element types introduced by Nédélec himself (Nédélec 1980,
1986), and more in the literature (Zaglmayr 2006). The simplest is the *first kind* family,
with two degrees of freedom per edge and two additional degrees of freedom on the face. The
edge functions are:

```
E^1 = s_1 * { 2 xi (2 xi - 1) grad eta  -  eta (4 xi - 1) grad xi  }
E^2 = s_1 * { xi (4 eta - 1) grad eta   -  2 eta (2 eta - 1) grad xi }
E^3 = s_2 * { 2 eta (2 eta - 1) grad zeta - zeta (4 eta - 1) grad eta }
E^4 = s_2 * { eta (4 zeta - 1) grad zeta  - 2 zeta (2 zeta - 1) grad eta }
E^5 = s_3 * { 2 zeta (2 zeta - 1) grad xi - xi (4 zeta - 1) grad zeta }
E^6 = s_3 * { zeta (4 xi - 1) grad xi     - 2 xi (2 xi - 1) grad zeta }
```

The three face functions are:

```
F^1 = 4 eta (eta - xi - 1) grad xi    -  4 xi (eta - xi + 1) grad eta
F^2 = 4 eta (eta - zeta - 1) grad zeta - 4 zeta (eta - zeta + 1) grad eta
F^3 = 4 xi (xi - zeta - 1) grad zeta   - 4 zeta (xi - zeta + 1) grad xi
```

They fulfill the property `F^1 + F^2 + F^3 = 0`, so any two of them span the face space while
the third is redundant and is not part of the interpolation vector. The factor 4 is not
required; it slightly improves the conditioning of the matrix.

When assigning the degree-of-freedom indices, remember that the order of the two edge dofs is
swapped when the sign of the edge is negative. The curl operator is computed with the same
technique as for the linear element. An implementation should distinguish straight-edged from
curved-edged second-order elements; for the latter, the derivatives of `grad xi` and `grad eta`
must be computed as well.

**Convention note (2026-08-14, TRI6 only):** the notes and the implementation
([cl_EF_TRI6.cpp](../nedelec/cl_EF_TRI6.cpp)) agree on these formulas, but the two parent-edge
functions each carry circulation `s/2` on their edge, not `1` and `0` (hand integral:
`H_0 - G_0 = 2 - 3t` and `H_1 - G_1 = 3t - 1`, both integrate to `1/2`). The blanket
unit-circulation statement in [nedelec.md](nedelec.md) §5 holds for the first-order elements;
for TRI6 (marked proof-of-concept there) the second-order edge dofs are moments of the Whitney
pair, and the formulas here are the authoritative record of what is implemented. Do not
renormalize one without the other. This note does **not** extend to TET10: the notes' TET10
edge polynomials are twice the TRI6 pair and integrate to unit circulation, see §4.1.

---

## 4. Nédélec Tetrahedron

*(source: `triangle/triangle_nedelec.tex`, edge table rewritten to the implemented convention)*

Element numbering schemes are not standardized across libraries, so node and edge numbers must
be defined explicitly. BELFEM follows the EXODUS II scheme (Schoof & Yarberry 1994). The
implemented coordinate assignment is (see
[cl_IF_TET4.hpp](../lagrange/cl_IF_TET4.hpp):28-31):

```
N_1 = xi,   N_2 = zeta,   N_3 = eta,   N_4 = tau = 1 - xi - eta - zeta
```

Note the swap: node 2 carries `zeta` and node 3 carries `eta`. This is a mean booby trap. On
the EXODUS triangle, node k naturally carries the k-th barycentric coordinate, so one would
assume the tetrahedron continues the pattern (node 2 with `eta`, node 3 with `zeta`). It must
not: on an EXODUS-ordered tetrahedron that naive map is left-handed, `det J < 0`, so the
coordinates of nodes 2 and 3 have to be exchanged to keep the Jacobian positive. Check the
map against [cl_IF_TET4.hpp](../lagrange/cl_IF_TET4.hpp) before writing any TET formula;
assuming the triangle pattern is exactly how the two edge-function defects of 2026-08-14
entered the code. With the EXODUS edge topology (edge k connects nodes (1,2), (2,3), (3,1),
(1,4), (2,4), (3,4)), the six Whitney pairs in the implemented convention are:

```
edge 1:  xi   -> zeta      E^1 = s_1 * ( xi   grad zeta - zeta grad xi   )
edge 2:  zeta -> eta       E^2 = s_2 * ( zeta grad eta  - eta  grad zeta )
edge 3:  eta  -> xi        E^3 = s_3 * ( eta  grad xi   - xi   grad eta  )
edge 4:  xi   -> tau       E^4 = s_4 * ( xi   grad tau  - tau  grad xi   )
edge 5:  zeta -> tau       E^5 = s_5 * ( zeta grad tau  - tau  grad zeta )
edge 6:  eta  -> tau       E^6 = s_6 * ( eta  grad tau  - tau  grad eta  )
```

(The notes list the pairs in the generic order `(xi,eta), (eta,zeta), (zeta,xi), ...`, which
does not match the implemented node-coordinate assignment; the table above is the one to code
against, cf. [cl_EF_TET4.cpp](../nedelec/cl_EF_TET4.cpp):179-213.)

> **Historical note (2026-08-14):** the edge 3 entry (`eta -> xi`, code index 2) of
> `EF_TET4::E()` used to compute `eta grad xi - xi grad zeta` instead of
> `eta grad xi - xi grad eta`, so `E` disagreed with its own curl operator `mC(:,2)` and the
> defective basis carried circulation 1/2 on its own edge and -1/2 on edge 1. Found in the
> extraction audit, confirmed by an exact symbolic probe cross-checked against DefElement's
> published degree-0 basis, and fixed the same day
> ([cl_EF_TET4.cpp](../nedelec/cl_EF_TET4.cpp):194-196 now matches the table above).

Using the partition of unity, `tau` is substituted with `1 - xi - eta - zeta` and `grad tau`
with `-(grad xi + grad eta + grad zeta)`. Deriving the edge functions and applying the chain
rule yields the curl operator; on the affine TET each column is again constant:

```
C(:,k) = 2 * s_k * grad lambda_a x grad lambda_b        for edge k = (a -> b)
```

which is what `EF_TET4::link()` precomputes ([cl_EF_TET4.cpp](../nedelec/cl_EF_TET4.cpp):113-138).

### 4.1 Second-Order Interpolation (TET10)

The edge interpolation functions follow the TRI6 pattern, two per edge, with `tau` taking the
role of the fourth coordinate on the edges towards node 4. In the notes' generic coordinate
labels (pairs `(xi,eta), (eta,zeta), (zeta,xi), (xi,tau), (eta,tau), (zeta,tau)`; remap the
pairing per the table in §4 before coding against it, cf. E16 in the drift record):

```
E^1  = s_1 * { 4 xi (2 xi - 1) grad eta      -  2 eta (4 xi - 1) grad xi    }
E^2  = s_1 * { 2 xi (4 eta - 1) grad eta     -  4 eta (2 eta - 1) grad xi   }
E^3  = s_2 * { 4 eta (2 eta - 1) grad zeta   -  2 zeta (4 eta - 1) grad eta }
E^4  = s_2 * { 2 eta (4 zeta - 1) grad zeta  -  4 zeta (2 zeta - 1) grad eta }
E^5  = s_3 * { 4 zeta (2 zeta - 1) grad xi   -  2 xi (4 zeta - 1) grad zeta }
E^6  = s_3 * { 2 zeta (4 xi - 1) grad xi     -  4 xi (2 xi - 1) grad zeta   }
E^7  = s_4 * { 4 xi (2 xi - 1) grad tau      -  2 tau (4 xi - 1) grad xi    }
E^8  = s_4 * { 2 xi (4 tau - 1) grad tau     -  4 tau (2 tau - 1) grad xi   }
E^9  = s_5 * { 4 eta (2 eta - 1) grad tau    -  2 tau (4 eta - 1) grad eta  }
E^10 = s_5 * { 2 eta (4 tau - 1) grad tau    -  4 tau (2 tau - 1) grad eta  }
E^11 = s_6 * { 4 zeta (2 zeta - 1) grad tau  -  2 tau (4 zeta - 1) grad zeta }
E^12 = s_6 * { 2 zeta (4 tau - 1) grad tau   -  4 tau (2 tau - 1) grad zeta }
```

Note that these polynomials are twice the TRI6 pair, and each parent-edge function integrates
to *unit* circulation on its edge, unlike the TRI6 convention of §3.3.

The face functions come in triples per face with the same redundancy property as on the
triangle (each triple sums to zero, two of three are active):

```
F^1  =  16 eta tau grad xi   -  8 tau xi grad eta   -  8 eta xi grad tau
F^2  = - 8 eta tau grad xi   + 16 tau xi grad eta   -  8 eta xi grad tau
F^3  = - 8 eta tau grad xi   -  8 tau xi grad eta   + 16 eta xi grad tau
F^4  =  16 tau zeta grad eta -  8 eta tau grad zeta -  8 eta zeta grad tau
F^5  = - 8 tau zeta grad eta + 16 eta tau grad zeta -  8 eta zeta grad tau
F^6  = - 8 tau zeta grad eta -  8 eta tau grad zeta + 16 eta zeta grad tau
F^7  =  16 tau xi grad zeta  -  8 tau zeta grad xi  -  8 xi zeta grad tau
F^8  = - 8 tau xi grad zeta  + 16 tau zeta grad xi  -  8 xi zeta grad tau
F^9  = - 8 tau xi grad zeta  -  8 tau zeta grad xi  + 16 xi zeta grad tau
F^10 =  16 eta zeta grad xi  -  8 eta xi grad zeta  -  8 xi zeta grad eta
F^11 = - 8 eta zeta grad xi  + 16 eta xi grad zeta  -  8 xi zeta grad eta
F^12 = - 8 eta zeta grad xi  -  8 eta xi grad zeta  + 16 xi zeta grad eta
```

(The notes dropped the `+` before the `16` term of `F^11`; fixed here.)

The multiplicities of the edge and face functions do not matter for the span (one might divide
the former by two and the latter by eight); what matters is that the same convention is used
consistently on both sides of a shared entity. As with the second-order triangle, the two dofs
of an edge are swapped when the edge sign is negative; whether that swapping happens in the
topology table or directly in the edge function vector is up to the programmer.

When two elements share a common face, a *face ownership* must be defined: for example, the
element with the lower ID owns the face, and dofs 1 and 2 of each face are active with respect
to the owner. The borrowing element must then find out how it is oriented with respect to the
owner and activate the corresponding face functions in its `E` and `C` operators.

> **Historical note (2026-08-14):** the scalar tables in
> [cl_EF_TET10.cpp](../nedelec/cl_EF_TET10.cpp) `precompute()` used to be transcribed in the
> notes' generic coordinate labels while the gradient assignments followed the implemented
> node map: eta and zeta were exchanged in every scalar factor touching them, the exact trap
> described in §4. The defect was found by an exact symbolic probe (24 edge-dof conformity
> violations, face candidates leaking onto edges), fixed at its source (the table generator
> `tmp/tet10/tet10_generate.m` used the naive node map), regenerated through the re-pinned
> MATLAB pipeline, and ported back; the final gate parsed the edited C++ tables and confirmed
> exact unit circulation, zero face leakage, and curl consistency for all 24 dofs. The same
> day, the new runtime battery (`tests/fem/test_EdgeFunctions.cpp`) caught a third instance of
> the class: rows 8-9 of the `mNeta`/`mNzeta` shape-derivative tables (which build the
> Jacobian) carried each other's values, making the eta and zeta Jacobian columns identical
> (singular J at generic points); fixed the same day, battery 22/22 green. TET10 remains
> proof-of-concept pending a full 3D regression run.

If the element has curved edges, the gradients `grad xi`, `grad eta`, `grad zeta` are no longer
constant and must be derived as well.

---

## 5. Edge and Face Generation

*(source: `triangle/triangle_nedelec.tex`; concept only, the implemented generation lives in `src/mesh`)*

Most mesh generators provide node coordinates and element topology; the edge information is
implicit. It can be generated with a sort-unique pass: allocate an array over all element
edges, assign each edge the key

```
d_e = a_e + b_e * n
```

where `a_e < b_e` are the two node numbers and `n` is the number of nodes. Populate by looping
over all elements, sort ascending, apply `unique`. The node numbers are recovered with
`a_e = d_e mod n` and `b_e = d_e div n`, and the edge direction can be defined as pointing from
`a_e` to `b_e`. The generation of faces in 3D works exactly the same way, except that `a_e` and
`b_e` now refer to the indices of the two elements sharing the face, and `n` to the total
number of elements on the mesh (keying by sorted node triples works just as well).

BELFEM's actual edge and face construction (including ownership and direction rules used by
`edge_directions()`) is implemented in the mesh module; see `src/mesh/doc/`. The formula above
is the concept, not a citation of the implementation.

---

## References

- **Nédélec (1980)**, "Mixed finite elements in R³", Numerische Mathematik 35, 315-341
- **Nédélec (1986)**, "A new family of mixed finite elements in R³", Numerische Mathematik 50
- **Zaglmayr (2006)**, "High Order Finite Element Methods for Electromagnetic Field Computation", PhD thesis, JKU Linz
- **Monk (2003)**, "Finite Element Methods for Maxwell's Equations", Ch. 5-6
- **Schoof & Yarberry (1994)**, "EXODUS II: A finite element data model", SAND92-2137
- Full citations: `doc/literature_references.md`
