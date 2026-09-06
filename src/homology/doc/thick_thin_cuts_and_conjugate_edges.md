# Thick Cuts, Thin Cuts, and Conjugate Edges {#homology_thick_thin_cuts_and_conjugate_edges}

**Date:** 2026-06-08 (cut cases + 2D worked example added 2026-06-12)
**Purpose:** Canonical reference for how BELFEM converts a cohomology thick cut into a FEM-compatible thin cut: the conjugate-edge structure, the cut-case convention as implemented in `CutData::determine_cut_case_2d/3d()`, a reference-tetrahedron example (3D), and a fully worked hexagon example (2D) covering duplication, relinking, and static condensation.
**Module:** src/homology

---

## Overview

BELFEM imposes scalar-potential (`φ`) jumps for the mixed h-φ formulation using cohomology
cuts. The cohomology engine first produces a **thick cut** — a set of directed edges with ±1
coefficients (the first cohomology generator H¹(K)). The `CutProcessor` then replaces it by an
equivalent **thin cut** living on element faces, which is what the FEM assembly actually
consumes. The mechanism is clearest on a single reference tetrahedron, worked out below.

This document is the worked-example companion to:

- [cohomology_theory_and_implementation.md](cohomology_theory_and_implementation.md) — §FEM
  Discretization derives the Poincaré–Lefschetz reduction in prose.
- [homology_usage_guide.md](homology_usage_guide.md) — Glossary entries *Thick Cut*, *Thin
  Cut*, *Poincaré-Lefschetz Duality*, *Static Condensation*.
- [cohomology_algorithms.md](cohomology_algorithms.md) — how the thick cut itself is computed.

---

## Reference tetrahedron

Vertices in barycentric coordinates. The swap between `p1` and `p2` is intentional and follows
the **EXODUS** local node ordering; in the mathematical reference orientation these coordinates
would be flipped back.

| node | coordinates  |
| ---- | ------------ |
| p0   | [1, 0, 0, 0] |
| p1   | [0, 0, 1, 0] |
| p2   | [0, 1, 0, 0] |
| p3   | [0, 0, 0, 1] |

The six local edges fall into two groups:

- **Face edges** (bounding face p0–p1–p2): `e0: 0→1`, `e1: 1→2`, `e2: 2→0`
- **Apex edges** (incident to vertex p3): `e3: 0→3`, `e4: 1→3`, `e5: 2→3`

---

## Thick cut

Assume the cohomology algorithm has found a **thick cut** whose discontinuity intersects the
three apex edges `e3, e4, e5`. The intersection points are their midsides:

```text
m3 = [0.5, 0,   0,   0.5]
m4 = [0,   0,   0.5, 0.5]
m5 = [0,   0.5, 0,   0.5]
```

The thick cut is then the small internal triangular surface **m3–m4–m5**, across which the
scalar potential `φ` carries a Heaviside-type discontinuity; it caps the apex p3 off from the
rest of the element. A discontinuity on this internal surface could be discretized directly with
an XFEM-like enrichment. **BELFEM does not do this.**

---

## Thin cut

Instead, BELFEM takes the Poincaré–Lefschetz viewpoint and replaces the internal thick cut by an
equivalent **thin cut**. Geometrically, the surface m3–m4–m5 is pushed down onto the external
face **p0–p1–p2**, the face opposite the apex p3. The boundary of this pushed thin cut is the
face loop formed by the edges `e0, e1, e2`.

The jump is then realized algebraically rather than by interior enrichment:

1. Duplicate the `φ` DOFs on the pushed cut face: **p0, p1, p2 → q0, q1, q2**.
2. Relink the **cut tetrahedron** to the duplicated DOFs q0, q1, q2 on one side of the thin cut;
   the **uncut element** on the opposite side of face p0–p1–p2 keeps the original DOFs p0, p1,
   p2.
3. Treat the duplicated DOFs as hanging DOFs and eliminate them by static condensation / change
   of basis.

The physical potential jump is thus represented without ever introducing an XFEM discontinuity
across the internal midside surface — which is the reason thin cuts are preferred.

### Where this happens in the code

The whole thick-to-thin pipeline lives in `cl_CutProcessor.cpp`, driven from
the `CutProcessor` **constructor** (`cl_CutProcessor.cpp:21-93`) — there is no `run()` method; the pipeline executes during construction:

| Step (above)            | Function                                  | Location |
| ----------------------- | ----------------------------------------- | -------- |
| Identify cut node sets  | `determine_cut_sets()` / `create_cut_sets()` | `cl_CutProcessor.cpp:75-76` |
| Flag the pushed face    | `compute_edge_bitsets()` / `compute_node_bitsets()` / `flip_node_bitsets_*()` | `cl_CutProcessor.cpp:78-80`, per-type `flip_node_bitsets_tet4/tet10/tri3/tri6` |
| Abstract jump DOFs      | `create_abstract_nodes()`                 | `cl_CutProcessor.cpp:82, 846` |
| Duplicate p0,p1,p2 → q0,q1,q2 | `duplicate_nodes()`                 | `cl_CutProcessor.cpp:84, 1154` |
| Relink the cut element  | `relink_elements()` / `relink_element()`  | `cl_CutProcessor.cpp:86, 1165` |
| Emit FEM cut sidesets   | `create_thin_cut_sidesets()`              | `cl_CutProcessor.cpp:90, 873` |

The hanging-DOF static condensation that eliminates q0, q1, q2 is performed later, in the FEM
DOF manager — see `DofData::create_dofwise_t_matrices_master()` in
`fem/kernel/cl_FEM_DofMgr_DofData.cpp` (and the mesh-level setup in
`fem/maxwell/cl_MaxwellFactory.cpp`). A transformation matrix `T` is built so that
`φ⁺ = T [φ⁻; I]`, and the local stiffness undergoes the change of basis `Tᵀ K_local T`.

---

## Conjugate edges

`e0, e1, e2` are the **conjugate edge set** of the cut edges `e3, e4, e5`: the thick cut
intersects the apex edges, while the corresponding thin cut lives on the opposite face, bounded
by the conjugate edges. (Under Poincaré–Lefschetz duality the conjugate-edge loop is the
boundary of the thin-cut surface, while the dual Ampère loop threads perpendicular to it — the
reason the two edge sets are called conjugate.)

| role          | edges          | geometry                                  |
| ------------- | -------------- | ----------------------------------------- |
| cut edges     | `e3, e4, e5`   | apex edges the **thick** cut intersects   |
| conjugate edges | `e0, e1, e2` | face loop bounding the **thin** cut       |

---

## Cut cases (code convention)

The cut case classifies how the thick cut crosses a single element. The authoritative
definitions are the pattern tables in `CutData::determine_cut_case_2d()` and
`determine_cut_case_3d()` (`cl_CutData.cpp`); everything below is derived from them.

### 2D reference triangle

Reference triangle in barycentric coordinates, **EXODUS** edge convention (not the
symmetric convention of analytical geometry):

| node | coordinates |
| ---- | ----------- |
| Q1   | [1, 0, 0]   |
| Q2   | [0, 1, 0]   |
| Q3   | [0, 0, 1]   |

with edges `R1: Q1→Q2`, `R2: Q2→Q3`, `R3: Q3→Q1`.

The thick cut crosses **exactly two** edges of the triangle; those carry coefficients
`±1` (one *in*, one *out* — the code asserts `tCoeffIn + tCoeffOut == 0`). The third
edge carries `0` and is the **conjugated edge**. The case index names the conjugated
edge; the sign encodes the jump orientation (**positive** case: the in-edge coefficient
is `+1` along the local edge direction):

| case | c(R1) | c(R2) | c(R3) | in | out | conjugated edge |
| ---- | ----- | ----- | ----- | -- | --- | --------------- |
| +1   |  0    | +1    | −1    | R2 | R3  | **R1**          |
| +2   | −1    |  0    | +1    | R3 | R1  | **R2**          |
| +3   | +1    | −1    |  0    | R1 | R2  | **R3**          |

Negative cases flip all coefficients (in↔out) and keep the same conjugated edge. In
the code's 0-based indexing this is the `face(tCase-1)` / `edge(tCase-1)` rule: **case
k ⇔ conjugated facet index k−1**.

### 3D tetrahedron

The same structure with faces (cf. the reference tetrahedron above). Three flagged
edges sharing a vertex cap that vertex off; the cut is pushed onto the **opposite
face**, which is the conjugated face:

| flagged edges (0-based) | pattern | case | conjugated face (0-based) |
| ----------------------- | ------- | ---- | ------------------------- |
| e1, e2, e5              | 38      | ±1   | face 0                    |
| e0, e2, e3              | 13      | ±2   | face 1                    |
| e0, e1, e4              | 19      | ±3   | face 2                    |
| e3, e4, e5              | 56      | ±4   | face 3                    |

Four flagged edges give the **diagonal cases** `±5, ±6, ±7` (patterns 43, 30, 53):
the cut crosses the element as a quadrilateral, no single conjugated face exists, and
the element contributes no thin-cut face of its own. (Topologically its portion of
the cut surface is expected to be covered by the conjugated faces of its face-case
neighbors; the code does not verify this face by face. As a protective heuristic, since 2026-09-02
`collect_facets` aborts when its dangling-face pruning removes more than a quarter of the
emitted faces, or all of them; the 2D branch aborts on an empty cut. The one observed case of such a loss, the corc_solder decks,
was a generator that is not a tight one-sheet representative: the positive-case faces do
not close, and the loop eats the surface inward from every hole.)

---

## Worked example in 2D: hexagon fan

Let `h = √3/2` and consider seven points

```text
P1: ( 0,    0 )      P2: ( 1,    0 )     P3: ( 0.5,  h )
P4: (-0.5,  h )      P5: (-1,    0 )     P6: (-0.5, -h )
P7: ( 0.5, -h )
```

forming six triangles fanning around P1:

```text
E1: 1,2,3    E2: 1,3,4    E3: 1,4,5
E4: 1,5,6    E5: 1,6,7    E6: 1,7,2
```

### Thick cut

The cohomology generator crosses four edges (oriented as written), each carrying
coefficient `−1`, so the imposed discontinuity is `Δφ = −I` (sign convention: `Δφ`
is measured original side minus duplicate side, consistent with `φ′ = φ + I` below):

```text
D1: 6→5     D2: 6→1     D3: 7→1     D4: 7→2
```

Note that `D1` and `D4` are *boundary* edges: the thick cut enters and exits the
domain through their midsides. The cut band consists of the elements crossed twice:
`E4, E5, E6`.

### Cut cases and conjugated facets

Working out each element's local edges against the coefficients above:

| element | local edges                 | cut edges | case | conjugated facet |
| ------- | --------------------------- | --------- | ---- | ---------------- |
| E4      | 1→5, 5→6, 6→1               | R2, R3    | **+1** | `F2 = (1,5)`   |
| E5      | 1→6, 6→7, 7→1               | R1, R3    | **−2** | (6,7) — boundary, not used |
| E6      | 1→7, 7→2, 2→1               | R1, R2    | **+3** | `F1 = (1,2)`   |

(Sign check for E4: its in-edge `R2 = 5→6` runs against `D1 = 6→5`, so the coefficient
along the local direction is `+1` → positive case. The same bookkeeping yields `−2`
for E5 and `+3` for E6 — all three labels follow from the single cochain above.)

The thin cut is formed by the conjugated facets of the **positive-case** elements
only: `F1 = (1,2)` and `F2 = (1,5)`. Their node set is `{1, 2, 5}`.

### Duplication and relinking

Create duplicates `P1′, P2′, P5′` and relink the **cut-band** elements:

```text
E4: 1′, 5′, 6     E5: 1′, 6, 7     E6: 1′, 7, 2′
```

Two distinct memberships govern this step:

- **band membership → relink**: every cut element (including the negative-case `E5`)
  is relinked to the primed nodes it touches;
- **conjugated-facet membership → duplicate**: only nodes of `F1 ∪ F2` are duplicated.
  `P6` and `P7` belong to the band but to no conjugated facet — every element touching
  them lies *inside* the band, so no jump crosses them and no duplicate is needed.

The discontinuity is now imposed **between the cut element and its neighbor across
each conjugated facet**: `E4 (1′,5′,6)` against `E3 (1,4,5)` across `F2`, and
`E6 (1′,7,2′)` against `E1 (1,2,3)` across `F1`. The midside discontinuity along
`D1…D4` has been shifted entirely onto the conjugated facets — this *is* the thin cut.

### Static condensation

Introduce one abstract node `A1` carrying the current `I` and assign the duplicated
DOFs the relationships

```text
φ′₁ = φ₁ + I        φ′₂ = φ₂ + I        φ′₅ = φ₅ + I
```

The primed DOFs are **hanging**: a coefficient matrix `T` expresses the element-local
DOFs in terms of the free DOFs `(φ₁, φ₂, φ₅, …, I)`, and the element matrix transforms
by the change of basis

```text
M_el = Tᵗ · M′_el · T
```

(`M` for the Ampère formulation, `K` for Gauss). This is the same mechanism as the
tetrahedron example above, executed by `DofData::create_dofwise_t_matrices_master()`.

---

## Relevance to periodic boundary conditions

This local picture matters for the PBC problem because the thin-cut construction depends not
only on *which* thick-cut edges were found, but also on correctly identifying the conjugate
face/edge structure, the side of the cut, the orientation, and any periodic equivalences of
nodes, edges, and faces. If any of these is lost or misidentified under periodicity, the pushed
thin cut can fail to reproduce the thick cut it was derived from, which is a natural source of
the continuity errors seen in periodic cases.

The hexagon example sharpens this. There, the thick cut **ends at a true boundary** (`D1`,
`D4` are boundary edges) and nothing special is needed: band nodes off the conjugated facets
(`P6`, `P7`) are simply not duplicated. On a periodic mesh the analogous "boundary" is a
**quotient seam** — the band continues on the other side through the identification. The corc
diagnosis settled which mechanism is at work, and it is worth stating precisely because two
*a priori* plausible pictures turned out differently.

**What is periodic and what is not.** `PeriodicityFactory` matches edges **only on the periodic
facets** — `update_periodicity()` collects edges from the master/slave seam facets and pairs
them (`update_periodicity()`, `cl_Mesh_PeriodicityFactory.cpp:336`; `collect_edges()`, `:600`). Therefore:

| entity | periodic? |
|--------|-----------|
| seam-plane node | matched to its mirror partner |
| edge **lying in** a seam facet | matched |
| **interior / slab-crossing** edge (seam node → interior, or a tet diagonal across the slab) | **not periodic — by design** |

**The confirmed seam structure (period-direction generators).** A cohomology generator that
**wraps the period** carries its cochain along *interior / slab-crossing* edges, which are not
periodic. The thin cut then reaches a seam node from one side of the identification but not the
other, so the seam node is "in the cut" while its periodic partner is not. This **one-sided cut
membership across a periodic pair is correct** — the two seam representatives genuinely lie on
opposite sides of a period-wrapping cut. It is *not* a flagging, orientation, or cochain bug:
all of those were measured symmetric across the seam (cochain support+sign by construction, node
flags by the `flag_nodes` partner closure, and `check_edge` head selection node-for-node). The
asymmetry enters only through the non-periodic interior edges above. Consequently the
former symmetric-membership assert in `CutSet::create_duplicates()` was **too strict** for these
generators and has been replaced by the three-way branch at `cl_CutSet.cpp:107-124`; one-sided
seam duplication is the correct response, no longer a hypothesis.

The earlier guess that the *selected conjugated facet lies in the seam plane itself* (an
"in-plane" cut) was **refuted** for corc: every asymmetry-carrying edge runs from a seam node
into the interior, none lies in the seam plane. (E5's conjugate `(6,7)` in the hexagon remains
only a *topological* analogy: E5 is a negative case, never emitted under the positive-only rule.)

**Global membership vs per-element signature.** A node's CutSet is keyed by its **element-local**
cut pattern: `flip_node_bitsets_*` writes the node into `mCutSets(localHex)` and relinking later
fetches the duplicate by the **same** local hex (`cl_CutProcessor.cpp:1243`). Two
periodic partners can therefore carry **identical global cut membership** yet land in
**different** CutSets, because their cut contributions arrive from different elements (near one
seam plane vs the other). This is a second, subtler face of the same seam structure and is
handled together with the one-sided case.

The duplication vocabulary is unchanged: band membership → relink, conjugated-facet membership →
duplicate, positive-case selection → single imposition. What periodicity adds is that the
*conjugated-facet membership* a node inherits can differ between periodic partners. The
implementation policy for all of this — the classification of seam pairs and the
duplication/rebuild rules — is stated inline at `cl_CutSet.cpp:107-124`, not
here.

---

## See Also

- **Theory:** [cohomology_theory_and_implementation.md](cohomology_theory_and_implementation.md)
- **Algorithms:** [cohomology_algorithms.md](cohomology_algorithms.md)
- **Usage:** [homology_usage_guide.md](homology_usage_guide.md)
- **Index:** [README.md](@ref homology_index)
