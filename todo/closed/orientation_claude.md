# Facet Integration Point Orientation Tests

**Date:** 2026-03-26
**Purpose:** Validate that facet integration points and orientation tables are correct
**Module:** mesh, fem/interpolation
**Sources:** Consolidated from orientation_claude.md, orientation_codex.md, orientation_junie.md

## Consensus Across All Three Reviews

All three reviews (Claude, Codex, Junie) agree on:

1. The proof-of-concept `test_facets.cpp` has a fundamental error: it uses 2D facet
   shape functions (TRI3/QUAD4) to evaluate points that live in the 3D volume element's
   parametric space. The volume shape function (TET4, HEX8, etc.) must be used instead.
2. The test must be split into two classes: master embedding and slave orientation.
3. The production FEM code is likely correct -- the error is in the test, not in
   `facetintpoints` or the orientation tables (high confidence).
4. Linear elements should be proven first; higher-order follows from the same logic.

## Key Insight

The `facetintpoints` functions (`intpoints_tet`, `intpoints_hex`, `intpoints_penta`)
produce integration points in the **volume element's parametric space**
(e.g., 4 barycentric coordinates for TET, 3 tensor-product coordinates for HEX).
Therefore, to verify them we must evaluate the **volume element's** shape function
(e.g., TET4), not the facet's shape function (e.g., TRI3).

There are two overloads per element type:
- **Master overload:** `intpoints_tet(masterFaceIndex, w, xi, order)` --
  maps standard 2D facet integration points into the master volume element's parametric space.
- **Slave overload:** `intpoints_tet(slaveFaceIndex, orientation, w, xi, order)` --
  maps the same physical points into the slave volume element's parametric space,
  accounting for the orientation mismatch.

## Shape Functions Needed

For each volume element type, we need both the volume shape function and the
facet shape function:

| Volume Element | Volume Shape Function       | Facet Type(s)  | Facet Shape Function |
|----------------|-----------------------------|----------------|----------------------|
| TET4           | N_tet4(xi, eta, zeta, tau)  | TRI3           | N_tri3(xi, eta)      |
| HEX8           | N_hex8(xi, eta, zeta)       | QUAD4          | N_quad4(xi, eta)     |
| PENTA6         | N_penta6(xi, eta, zeta)     | QUAD4 + TRI3   | N_quad4, N_tri3      |
| PYRA5          | N_pyra5(xi, eta, zeta)      | TRI3 + QUAD4   | N_tri3, N_quad4      |

Volume shape functions (linear):
- TET4:   `N = { xi, eta, zeta, 1 - xi - eta - zeta }` (barycentric, 4 coords)
- HEX8:   Standard trilinear `(1 +/- xi)(1 +/- eta)(1 +/- zeta) / 8`
- PENTA6: Triangular base `(xi, eta, 1-xi-eta)` x linear axial `(1 +/- zeta)/2`
- PYRA5:  Rational shape functions (special case at apex)

Facet shape functions (linear):
- TRI3:  `N = { xi, eta, 1 - xi - eta }`
- QUAD4: `N = { (1-xi)(1-eta)/4, (1+xi)(1-eta)/4, (1+xi)(1+eta)/4, (1-xi)(1+eta)/4 }`

## Class 1: Master Facet Integration Points

**Goal:** Prove that the master-mode `facetintpoints` correctly embed 2D facet
integration points into the volume element's parametric space.

**No orientation table needed for this class.**

**Method:**
For each volume element type E (TET4, HEX8, PENTA6, PYRA5), for each face f:

1. Create a `ReferenceElement` for E with known node coordinates X.
2. Get the facet node coordinates X_facet from `get_nodes_of_facet(f)`.
3. Generate standard 2D integration points on the facet geometry:
   `intpoints(GAUSS, facet_geometry, order, w_2d, xi_2d)`
4. Generate master volume integration points:
   `intpoints_<type>(f, w_3d, xi_3d, order)`
5. For each integration point k:
   - `p_k = N_facet(xi_2d(:,k)) * X_facet`  (2D facet shape function, facet node coords)
   - `q_k = N_volume(xi_3d(:,k)) * X`        (3D volume shape function, all node coords)
   - Assert `||p_k - q_k|| < epsilon`

This proves the embedding is geometrically correct: the 3D parametric
coordinates produced by the master facetintpoints, when evaluated with
the volume shape function, land at the same physical locations as the
2D parametric coordinates evaluated with the facet shape function.

## Class 2: Slave Facet Integration Points

**Goal:** Prove that the slave-mode `facetintpoints` with a given orientation
produce volume-parametric integration points that correspond to the same
physical locations as the 2D facet integration points (with permuted nodes).

**This is where the orientation table matters.**

**Method:**
For each volume element type E, for each face f, for each orientation o:

1. Same `ReferenceElement` and orientation table as above.
2. Get the facet nodes from `get_nodes_of_facet(f)` -> X_facet (slave perspective).
3. Use the orientation table to get the permuted facet node coordinates X_facet_perm
   (master perspective for orientation o).
4. Generate standard 2D integration points:
   `intpoints(GAUSS, facet_geometry, order, w_2d, xi_2d)`
5. Generate slave volume integration points:
   `intpoints_<type>(f, o, w_3d, xi_3d, order)`
6. For each integration point k:
   - `p_k = N_facet(xi_2d(:,k)) * X_facet_perm`  (2D shape function, permuted facet coords)
   - `q_k = N_volume(xi_3d(:,k)) * X`              (3D volume shape function, all node coords)
   - Assert `||p_k - q_k|| < epsilon`

This proves the slave integration points, when evaluated with the volume
shape function, land at the same physical locations as the 2D points
evaluated with the orientation-permuted facet node coordinates.

## Implementation Plan

### Step 1: TET4 -- prove linear tet first

- [ ] Implement Class 1 (master) for TET4: use N_tet4 volume shape function
- [ ] Implement Class 2 (slave) for TET4: use N_tet4 + orientation table
- [ ] Verify all 4 faces x 3 orientations pass

### Step 2: Extend to remaining linear 3D elements

- [ ] HEX8 (6 QUAD4 faces, 4 orientations each)
- [ ] PENTA6 (3 QUAD4 faces + 2 TRI3 faces, mixed orientations)
- [ ] PYRA5 (4 TRI3 faces + 1 QUAD4 face) -- **NOTE: slave intpoints_pyra not yet implemented**

### Step 3: Higher-order elements (orientation tables only)

Once linear elements pass, higher-order elements (TET10, HEX20/27, PENTA15/18,
PYRA13/14) share the same integration point logic -- only the orientation tables
add midedge and face-center node permutations. These can be verified by:
- [ ] Checking that higher-order orientation tables are consistent extensions
      of the linear tables (same corner node pattern, plus correct midedge nodes)
- [ ] Optionally: brute-force enumeration in Python to independently generate
      the higher-order permutation tables

### Step 4: Convert to Google Test (later)

- [ ] Convert the validated test_facets.cpp into proper GTest cases
- [ ] One TEST per element type, with subtests per face and orientation

## Known Gaps (from Junie review)

1. **PYRA slave integration points:** `facetintpoints::intpoints_pyra(slave, orient, ...)`
   does not exist yet. Must be implemented before PYRA5 slave tests.
2. **PYRA13 orientation table:** `create_orientation_table` has PYRA5 and PYRA14 but
   not PYRA13.
3. **Shape function approach:** Junie suggests using `InterpolationFunction` objects
   from the factory instead of manual function pointers. This is cleaner and
   scales to higher-order elements automatically.

## Notes

- The 2D integration points `xi_2d` are the "ground truth" reference.
- The orientation table columns are ordered as: face 0 / orient 0, face 0 / orient 1, ...,
  face F / orient O (sequential, `tRow` counter in the proof of concept).
- Quad faces have 4 orientations, tri faces have 3 orientations.
- The failure pattern in the original test (face 3 passing, others failing for TET4)
  is consistent with the coordinate mismatch: face 3 nodes [0,2,1] happen to partly
  align with the first two barycentric coordinates read by the TRI3 function (Codex).
