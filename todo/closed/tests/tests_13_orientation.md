# BELFEM Element Orientation Tests — Detailed Plan

**Date:** 2026-03-22 (completed from stub)
**Purpose:** Test plan for facet orientation permutations at element interfaces
**Depends on:** `tests_11_interpolation.md`, `tests_12_integration.md`, `tests_10_mesh.md`
**Status:** **COMPLETE** — Source investigation performed 2026-03-22. All permutations enumerated, test matrix designed.
**Confidence:** High on permutation catalog and index formulas (verified against source). High on TET/HEX test design. Medium on PENTA (complex mixed-face logic). Low on end-to-end mesh-based tests (require infrastructure that does not yet exist).

---

## Why This Matters

When two elements share a facet (edge in 2D, face in 3D), the slave element's local node ordering on that facet may be rotated or reflected relative to the master's. If the integration point mapping from master to slave does not account for this orientation correctly, the quadrature points on both sides of the interface evaluate at **different physical locations**. This causes silent, catastrophic errors in:

- h-φ interface coupling integrals
- Nédélec edge function assembly across element boundaries
- Any sideset-based integral (thin shells, boundary conditions)
- Static condensation of bubble enrichment on the slave side

The orientation handling is implemented in `IntegrationData::populate_for_slave()` and the `facetintpoints::intpoints_tet/hex/penta` functions, using the `(aSlaveIndex, aOrientation)` parameter pair. This plan ensures that every valid permutation maps correctly.

---

## Source Investigation Findings

### Orientation Encoding

**Type:** `uint`, **1-based**.

**Computation** (in `cl_Facet.cpp:104–130` and `cl_Face.cpp:113–150`):
1. Get corner nodes of the slave's facet via `get_corner_nodes_of_facet(slaveFaceIndex)`
2. Take the first slave corner node's ID
3. Get corner nodes of the master's facet
4. Find the position (0-based) in the master corner node list where that ID appears
5. Return `position + 1`

**Storage:** `mOrientationOnSlave` in `Facet` (`cl_Facet.hpp:41`) and `Face` (`cl_Face.hpp:39`).

### Valid Orientations by Geometry

From `number_of_orientations()` in `meshtools.cpp:781–827`:

| Parent Element | Facet Type | Facets | Orientations/Facet | Total Permutations |
|---|---|---|---|---|
| TRI/QUAD (2D) | LINE edges | 3/4 | 1 | 3/4 |
| TET | TRI faces | 4 | 3 | 12 |
| HEX | QUAD faces | 6 | 4 | 24 |
| PENTA | QUAD + TRI | 3 + 2 | 4 + 3 | 18 |
| PYRA | TRI + QUAD | 4 + 1 | 3 + 4 | 16 |

### Permutation Index Formulas (VERIFIED)

| Geometry | Formula | Source |
|---|---|---|
| 2D | `index_on_slave()` | `cl_FEM_Calculator.cpp:378` |
| TET | `index_on_slave() * 3 + orientation_on_slave() - 1` | `cl_FEM_Calculator.cpp:388` |
| HEX | `index_on_slave() * 4 + orientation_on_slave() - 1` | `cl_FEM_Calculator.cpp:398` |

The `SideSet` pre-populates a flat array `mSlaveIntegration` at `cl_FEM_SideSet.cpp:506–528`:
```
tCount = 0
for face f = 0..numFacets-1:
    for orientation o = 0..numOrientations(f)-1:
        mSlaveIntegration[tCount++].populate_for_slave(f, o, ...)
```

Note: `populate_for_slave` receives orientation **0-based** from this loop. The `orientation_on_slave()` accessor returns **1-based**. The `-1` in the index formula converts back.

### Implementation Details

**TET slave** (`fn_IF_initialize_integration_points_on_facet.cpp:790–825`):
Uses a 4×12 lookup table `tIndex` mapping 3 TRI barycentric coordinate rows to 4 TET barycentric coordinate rows. Case index = `aSlaveIndex * 3 + aOrientation`. The 4th row (opposite vertex) is implicitly zero.

```
tIndex = { { 0, 2, 3, 2, 1, 3, 0, 3, 1, 0, 1, 2 },
           { 3, 0, 2, 3, 2, 1, 1, 0, 3, 2, 0, 1 },
           { 2, 3, 0, 1, 3, 2, 3, 1, 0, 1, 2, 0 },
           { 1, 1, 1, 0, 0, 0, 2, 2, 2, 3, 3, 3 } };
```

**HEX slave** (lines 829–871):
Uses a 2×24 axis-mapping table `tIndex` (1-based axis with sign for flip) plus a 24-element sign vector `tSign` for the fixed-coordinate value. Case index = `aSlaveIndex * 4 + aOrientation`.

**PENTA slave** (lines 552–785):
Handles quad faces (faces 0–2, case = `face*4+orientation`, 12 cases) and tri faces (faces 3–4, case = `(face-3)*3+orientation`, 6 cases) with separate switch blocks.

**2D slave** (lines 156–188, `populate_for_slave_tri`):
Simply reverses the integration point order on the edge. Orientation parameter is ignored.

---

## Known Issues Found During Source Investigation

### ISSUE-O1: Hardcoded `*3` in `Calculator::link(Facet*)` (potential bug)

**File:** `cl_FEM_Calculator.cpp:1014–1016`

The `link(Facet*)` method uses `index_on_slave() * 3 + orientation_on_slave() - 1` for ALL element types. This is correct for TET but **wrong for HEX** (should use `*4`). Currently this code path is only used via `get_normal_calculator()` which may only be called for TET-based sidesets in practice, making this a latent bug.

### ISSUE-O2: No PENTA support in Calculator dispatch

**File:** `cl_FEM_Calculator.cpp:785–807`

The switch for `mFunSlaveIntegration` only handles TRI/QUAD, TET, and HEX. PENTA hits the default error case. PENTA-as-slave elements are not supported at the kernel level, even though the integration point permutations are fully implemented.

### ISSUE-O3: No PYRA slave implementation

`number_of_orientations()` supports PYRA (16 total), but no `facetintpoints::intpoints_pyra` slave variant exists, and there is no PYRA handling in `populate_for_slave`.

### ISSUE-O4: Existing test cannot compile

**File:** `test/fem/cl_IntegrationData_Interface.cpp`

This file covers TRI3, TRI6, QUAD4, QUAD9, TET4, TET10, HEX8, HEX27 but depends on nonexistent `FEM_geometry.hpp` and `geometry::create_test_mesh()`. The test pattern is correct but the infrastructure must be created.

---

## Permutation Catalog

### Triangular faces (TRI on TET, TRI on PENTA)

A triangular face has **3 vertices**. Given a fixed master ordering (A, B, C), the slave can see the same face as:

| orientation_on_slave | Slave ordering | Rotation | 0-based (in flat array) |
|---|---|---|---|
| 1 | (A, B, C) | Identity | 0 |
| 2 | (B, C, A) | 120° rotation | 1 |
| 3 | (C, A, B) | 240° rotation | 2 |

### Quadrilateral faces (QUAD on HEX, QUAD on PENTA)

A quad face has **4 vertices**. Given master ordering (A, B, C, D), the slave can see:

| orientation_on_slave | Slave ordering | Rotation | 0-based (in flat array) |
|---|---|---|---|
| 1 | (A, B, C, D) | Identity | 0 |
| 2 | (B, C, D, A) | 90° rotation | 1 |
| 3 | (C, D, A, B) | 180° rotation | 2 |
| 4 | (D, A, B, C) | 270° rotation | 3 |

### 2D edges (LINE on TRI, LINE on QUAD)

An edge has **2 vertices**. The slave sees the edge reversed:

| Orientation | Effect |
|---|---|
| (any) | Integration point order reversed |

This is handled by `populate_for_slave_tri` which simply reverses the point order.

---

## Floating-Point Comparison

```cpp
namespace
{
    const belfem::real tEps = 1e-12 ;  // for physical-space coordinate coincidence
}
```

---

## Test File Structure

```
tests/fem/
├── test_FacetOrientation2D.cpp    # TRI-TRI and QUAD-QUAD edge interfaces
├── test_FacetOrientationTet.cpp   # TET-TET face interfaces (12 permutations)
├── test_FacetOrientationHex.cpp   # HEX-HEX face interfaces (24 permutations)
├── test_FacetOrientationPenta.cpp # PENTA face interfaces (18 permutations)
```

---

## Testing Approach

### Layer 1: Direct Function Tests (no mesh required)

Call `facetintpoints::intpoints_tet()` / `intpoints_hex()` / `intpoints_penta()` directly with known `(face_index, orientation)` pairs. These functions populate a `Matrix<real>` of integration point coordinates. Verify that:

1. Master and slave integration points, when mapped to physical space through shape functions with manually-constructed node coordinates, coincide.
2. Weights from master and slave are identical.

**Minimum includes needed:**
- `cl_IF_IntegrationData.hpp`
- `fn_IF_initialize_integration_points_on_facet.hpp`
- `cl_IF_InterpolationFunctionFactory.hpp`
- `fn_intpoints.hpp`
- `Mesh_Enums.hpp`

### Layer 2: Mesh-Based Tests (requires infrastructure)

Construct actual `Element` and `Facet` objects with known node coordinates. Test `compute_orientation()` and the full `populate_for_master` / `populate_for_slave` pipeline. This requires building the `create_test_mesh` infrastructure referenced in the existing (non-compilable) `test/fem/cl_IntegrationData_Interface.cpp`.

**Recommendation:** Implement Layer 1 first. Layer 2 can follow once the mesh test infrastructure is available.

---

## 1. TET-TET Orientation Tests (12 permutations)

**File:** `test_FacetOrientationTet.cpp`

### 1.1 Full Permutation Sweep (Parameterized) `[semantic]`

For each `(face_index, orientation)` pair, verify master-slave integration point coincidence in physical space.

| face_index | orientation (0-based) | tCase | Face Nodes (TET4) | Description |
|---|---|---|---|---|
| 0 | 0 | 0 | (0,1,3) | Face 0, identity |
| 0 | 1 | 1 | (0,1,3) | Face 0, 120° |
| 0 | 2 | 2 | (0,1,3) | Face 0, 240° |
| 1 | 0 | 3 | (1,2,3) | Face 1, identity |
| 1 | 1 | 4 | (1,2,3) | Face 1, 120° |
| 1 | 2 | 5 | (1,2,3) | Face 1, 240° |
| 2 | 0 | 6 | (0,2,3) | Face 2, identity |
| 2 | 1 | 7 | (0,2,3) | Face 2, 120° |
| 2 | 2 | 8 | (0,2,3) | Face 2, 240° |
| 3 | 0 | 9 | (0,2,1) | Face 3, identity |
| 3 | 1 | 10 | (0,2,1) | Face 3, 120° |
| 3 | 2 | 11 | (0,2,1) | Face 3, 240° |

Test both **TET4** (linear) and **TET10** (quadratic).

### 1.2 Implementation Pattern

```
1. Define two TET4 elements sharing face 0 with known node coordinates
   - Master: nodes at (0,0,0), (1,0,0), (0,1,0), (0,0,1)
   - Slave: shares face (0,1,3) with a different 4th node
2. For each (face_index, orientation):
   a. Call populate_for_master(face_index, integration_order)
   b. Call populate_for_slave(face_index, orientation, integration_order)
   c. For each integration point k:
      - Map master point to physical space: x_master = Σ N_i(ξ_master_k) · x_i
      - Map slave point to physical space: x_slave = Σ N_j(ξ_slave_k) · x_j
      - EXPECT_NEAR(x_master, x_slave, tEps)
3. Verify weights match: w_master[k] == w_slave[k]
```

### 1.3 Lookup Table Verification `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `TetLookupTableRowSumIsConstant` | For each case column in `tIndex`, rows 0–3 contain a permutation of {0,1,2,3} |
| `TetLookupRow4IsOppositeVertex` | Row 4 of `tIndex` gives the zero-row (vertex opposite to the face) |
| `TetZeroCoordOnOpposite` | The barycentric coordinate in the "opposite vertex" row is always 0 |

---

## 2. HEX-HEX Orientation Tests (24 permutations)

**File:** `test_FacetOrientationHex.cpp`

### 2.1 Full Permutation Sweep (Parameterized) `[semantic]`

| face_index | orientations | tCase range | Face Description |
|---|---|---|---|
| 0 | 0, 1, 2, 3 | 0–3 | Front face (nodes 0,1,5,4) |
| 1 | 0, 1, 2, 3 | 4–7 | Right face (nodes 1,2,6,5) |
| 2 | 0, 1, 2, 3 | 8–11 | Back face (nodes 2,3,7,6) |
| 3 | 0, 1, 2, 3 | 12–15 | Left face (nodes 0,4,7,3) |
| 4 | 0, 1, 2, 3 | 16–19 | Bottom face (nodes 0,3,2,1) |
| 5 | 0, 1, 2, 3 | 20–23 | Top face (nodes 4,5,6,7) |

Test both **HEX8** (linear) and **HEX27** (quadratic).

### 2.2 Implementation Pattern

Same as TET but with hex shape functions and unit-cube node coordinates:
- Master: standard unit cube [0,1]³
- Slave: second cube sharing one face, with nodes ordered to produce the desired orientation

### 2.3 Axis Mapping Verification `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `HexAxisTableCoversAllAxes` | Each case in `tIndex` maps exactly 2 of the 3 axes (xi, eta, zeta) |
| `HexSignTableInRange` | Each entry in `tSign` is ±1 (the fixed coordinate is at +1 or -1 face) |
| `HexOrthogonalMapping` | The two mapped axes are orthogonal in parametric space |

---

## 3. PENTA Orientation Tests (18 permutations)

**File:** `test_FacetOrientationPenta.cpp`

### 3.1 Quad Faces (faces 0–2, 12 permutations) `[semantic]`

| face_index | orientations | tCase range | Description |
|---|---|---|---|
| 0 | 0, 1, 2, 3 | 0–3 | Quad face 0 |
| 1 | 0, 1, 2, 3 | 4–7 | Quad face 1 |
| 2 | 0, 1, 2, 3 | 8–11 | Quad face 2 |

### 3.2 Tri Faces (faces 3–4, 6 permutations) `[semantic]`

| face_index | orientations | tCase range | Description |
|---|---|---|---|
| 3 | 0, 1, 2 | 0–2 | Bottom tri face |
| 4 | 0, 1, 2 | 3–5 | Top tri face |

Test both **PENTA6** (linear) and **PENTA18** (quadratic).

**Note:** PENTA-as-slave is NOT supported at the Calculator dispatch level (ISSUE-O2). These tests verify the raw integration point permutation functions, not the full kernel path.

---

## 4. 2D Edge Orientation Tests

**File:** `test_FacetOrientation2D.cpp`

### 4.1 Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `Tri3EdgeReversalFace0` | For TRI3, face 0: slave integration points are reversed copy of master |
| `Tri3EdgeReversalFace1` | Same for face 1 |
| `Tri3EdgeReversalFace2` | Same for face 2 |
| `Tri6EdgeReversalAllFaces` | For TRI6, all 3 faces: slave points reversed, mapped to same physical coords |
| `Quad4EdgeReversalAllFaces` | For QUAD4, all 4 faces |
| `Quad9EdgeReversalAllFaces` | For QUAD9, all 4 faces |

### 4.2 Implementation Pattern

```
1. Define two TRI3 elements sharing edge 0
   - Master: (0,0), (1,0), (0,1)
   - Slave: shares edge (0,0)-(1,0) with a different 3rd node at (0.5, -1)
2. Populate master integration points on face 0
3. Populate slave integration points on face 0
4. Map both to physical coordinates on the shared edge
5. Verify coincidence
```

---

## 5. Weight Consistency `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `TetWeightsMatchAcrossOrientations` | For all 12 TET permutations: slave weights == master weights |
| `HexWeightsMatchAcrossOrientations` | For all 24 HEX permutations: slave weights == master weights |
| `PentaWeightsMatchAcrossOrientations` | For all 18 PENTA permutations: slave weights == master weights |
| `2DWeightsMatchAcrossEdges` | For all 2D edge cases: slave weights == master weights |

---

## 6. Index Formula Regression Tests `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `TetIndexFormula` | For face f=0..3, orientation o=1..3: flat index == `f * 3 + (o - 1)` |
| `HexIndexFormula` | For face f=0..5, orientation o=1..4: flat index == `f * 4 + (o - 1)` |
| `PentaQuadIndexFormula` | For quad face f=0..2, orientation o=0..3: flat index == `f * 4 + o` |
| `PentaTriIndexFormula` | For tri face f=3..4, orientation o=0..2: flat index == `(f-3) * 3 + o` |

---

## 7. Regression Tests for Known Issues

### 7.1 ISSUE-O1 Regression `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CalculatorLinkUsesCorrectMultiplierTet` | `Calculator::link(Facet*)` with TET slave uses `*3` → correct flat index |
| `CalculatorLinkUsesCorrectMultiplierHex` | `Calculator::link(Facet*)` with HEX slave uses `*4` → correct flat index |

**Note:** This test will expose BUG ISSUE-O1 if the hardcoded `*3` is still present for HEX. Mark as expected-to-fail until fixed.

### 7.2 ISSUE-O2: PENTA Calculator Dispatch `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `CalculatorPentaSlaveNotSupported` | Attempting to set up PENTA as slave in Calculator → `BELFEM_ERROR` (documents the limitation) |

---

## 8. What We Do NOT Test

- **PYRA orientation** — no slave implementation exists (ISSUE-O3). Document gap only.
- **Mesh-level `compute_orientation()`** — requires full mesh infrastructure (Layer 2). Deferred.
- **Reflected (mirrored) orientations** — BELFEM only handles rotations, not reflections. Consistent with conformal meshes where element handedness is preserved.
- **Performance** — per strategy doc.
- **Nédélec edge function orientation** — deferred to future edge function test plan.

---

## 9. Implementation Notes for Claude Code

1. **Use Layer 1 approach first.** Call `facetintpoints::intpoints_tet()` etc. directly. Avoid mesh construction until infrastructure exists.
2. **Hardcode node coordinates.** Use simple geometries: unit tet, unit cube, unit prism. Document coordinates in test code comments.
3. **Parameterized tests are essential.** A `TEST_P` over `{ElementType, face_index, orientation}` structs covers 12+24+18 = 54 cases without code duplication.
4. **The existing test file `test/fem/cl_IntegrationData_Interface.cpp` has the right pattern** but cannot compile. Use it as a reference for the test structure, but implement with direct function calls.
5. **Shape function creation:** Use `InterpolationFunctionFactory::create_lagrange_function()` for both master and slave elements.
6. **Integration order:** Use `auto_integration_order()` to get the recommended order for each element type.
7. **Physical-space mapping:** For each integration point, compute `x_phys = Σ N_i(ξ) · x_node_i` where `N` comes from the shape function evaluated at the parametric coordinates.
8. **TET barycentric coordinates:** TRI integration points have 3 barycentric rows. The TET slave function maps these to 4 barycentric rows (one per TET vertex). The row that corresponds to the opposite vertex gets coordinate 0.
9. **HEX axis mapping:** The `tIndex` table uses 1-based axis indices with sign for flip. `|tIndex| = 1` means xi, `2` means eta, `3` means zeta. Negative means flip direction.
10. **PENTA has two separate switch blocks** — one for quad faces (12 cases) and one for tri faces (6 cases). The test parameterization must account for this split.

---

## 10. Codex Audit Checklist

- [ ] All 12 TET permutations tested with physical-space coincidence
- [ ] All 24 HEX permutations tested with physical-space coincidence
- [ ] All 18 PENTA permutations tested (12 quad + 6 tri)
- [ ] 2D edge reversal tested for TRI3, TRI6, QUAD4, QUAD9
- [ ] Weight consistency verified for all permutations
- [ ] Index formula regression tests present
- [ ] ISSUE-O1 (hardcoded `*3`) regression test present (expected-to-fail until fixed)
- [ ] ISSUE-O2 (PENTA not in Calculator) documented
- [ ] ISSUE-O3 (PYRA not implemented) documented
- [ ] Node coordinates hardcoded (no mesh I/O dependency)
- [ ] Both linear and quadratic element orders tested (TET4+TET10, HEX8+HEX27, PENTA6+PENTA18)
- [ ] `EXPECT_NEAR` for all floating-point comparisons
- [ ] BELFEM naming conventions (`t` prefix for locals)
