# BELFEM Integration Module Tests — Detailed Plan

**Date:** 2026-03-22
**Purpose:** Method-level test matrix for the integration module (quadrature rules)
**Depends on:** `tests_0_strategy.md` (conventions), `tests_2_linalg.md` (Vector/Matrix), `tests_11_interpolation.md` (shape functions)
**Confidence:** High on weight-sum and containment invariants (universal). High on dispatcher routing (code is a switch cascade). High on bugs BUG-Q1 and BUG-Q2 (visible in source). Medium on polynomial exactness (depends on correct reference volumes). Low on end-to-end gmsh tests (depends on mesh I/O availability).

---

## Module Overview

The integration module provides Gauss quadrature rules for all reference geometries used in BELFEM. It consists of ~70 hardcoded point tables sourced from published literature and the quadpy library, a dispatcher that selects the correct table, and supporting utilities.

| Layer | Content | Test Priority |
|---|---|---|
| **Dispatcher** | `intpoints(scheme, geometry, order, weights, points)` — routes to correct rule | **Highest** |
| **1D kernels** | `gauss_line`, `gauss_quad`, `gauss_hex` — tensor-product builders using Fortran `intpoints_gauss` | High |
| **Simplex tables** | TRI (1–79 points), TET (1–236 points), PENTA (1–71 points), PYRA (1–125 points) | High |
| **Optimal QUAD/HEX** | Non-tensor-product rules: QUAD (8–37 points), HEX (6–34 points) — Hammer & Stroud | Medium |
| **auto_integration_order** | Maps ElementType → recommended integration order | Medium |
| **IntegrationScheme enum** | `to_string`, `string_to_integration_scheme` with normalization | Medium |
| **Lobatto** | 1D Lobatto rules via Fortran `intpoints_lobatto` | Low |

### Point Data Provenance

Integration points and weights are copy-pasted from original papers and the quadpy library. Typos in the data tables are extremely unlikely. The weights are chosen to be compatible with BELFEM's interpolation functions and reference domain conventions.

---

## Reference Domain Conventions

| Geometry | Coord rows | Reference volume | Domain description |
|---|---|---|---|
| LINE | 1 | 2.0 | ξ ∈ [-1, 1] |
| QUAD | 2 | 4.0 | (ξ, η) ∈ [-1, 1]² |
| HEX | 3 | 8.0 | (ξ, η, ζ) ∈ [-1, 1]³ |
| TRI | 3 (barycentric) | 0.5 | L₁, L₂, L₃ ≥ 0, L₁ + L₂ + L₃ = 1 |
| TET | 4 (barycentric) | 1/6 | L₁, L₂, L₃, L₄ ≥ 0, Σ = 1 |
| PENTA | 3 | 1.0 | Tri-base (ξ, η ≥ 0, ξ+η ≤ 1), prism axis ζ ∈ [-1, 1] |
| PYRA | 3 | 4/3 | Quad-base (ξ, η), apex at ζ = 1, base at ζ = 0 |

**VERIFIED 2026-03-22:** All coordinate row counts, reference volumes, and simplest-rule weights confirmed against source. TRI uses 3-row barycentric (L1,L2,L3), TET uses 4-row barycentric (L1,L2,L3,L4). PENTA uses 2 Cartesian triangle coords + 1 axial coord (NOT barycentric). `auto_integration_order` mapping also verified: CONSTANT→1, LINEAR→4, SERENDIPITY→7, QUADRATIC→7, CUBIC→10, QUARTIC→13, QUINTIC→16.

---

## Known Bugs

| ID | Location | Issue | Regression Test |
|---|---|---|---|
| BUG-Q1 | `fn_intpoints.cpp` ~line 727 | **Missing `break` after PYRA case.** The PYRA geometry block in the GAUSS scheme dispatcher does not end with `break;` before the outer `default:`. Valid PYRA orders will fall through and fire `BELFEM_ERROR("Invalid geometry type")`. | Test: `intpoints(GAUSS, PYRA, 1, w, p)` should succeed, not error. |
| BUG-Q2 | `fn_intpoints.cpp` ~line 585 | **Wrong error message in PENTA default.** Says `"Invalid Order for GeometryType::TRI"` but is inside the PENTA case. Copy-paste error. | Document; low-priority fix. |

---

## Floating-Point Comparison

```cpp
namespace
{
    const belfem::real tEps = 1e-12 ;  // for weight sums, barycentric constraints
    const belfem::real tTol = 1e-9 ;   // for polynomial exactness (accumulated FP error)
}
```

---

## Test File Structure

```
tests/numerics/
├── test_IntegrationScheme.cpp     # Enum/string conversion
├── test_IntegrationDispatch.cpp   # Dispatcher routing, auto_integration_order
├── test_IntegrationInvariants.cpp # Weight sums, containment, polynomial exactness
├── test_IntegrationEndToEnd.cpp   # Shape function × quadrature on unit-volume elements
```

---

## 1. IntegrationScheme Enum

**File:** `test_IntegrationScheme.cpp`

### 1.1 to_string `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `ToStringGauss` | `to_string(GAUSS)` → `"GaussModern"` |
| `ToStringGaussClassic` | `to_string(GAUSSCLASSIC)` → `"GaussClassic"` |
| `ToStringLobatto` | `to_string(LOBATTO)` → `"Lobatto"` |
| `ToStringUndefined` | `to_string(UNDEFINED)` → `"undefined"` |

### 1.2 string_to_integration_scheme `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `ParseGauss` | `"gauss"` → `GAUSS` |
| `ParseGaussModern` | `"gaussmodern"` → `GAUSS` |
| `ParseGaussClassic` | `"gaussclassic"` → `GAUSSCLASSIC` |
| `ParseLobatto` | `"lobatto"` → `LOBATTO` |
| `ParseCaseInsensitive` | `"GAUSS"`, `"Gauss"`, `"GaUsS"` all → `GAUSS` |
| `ParseWithSpaces` | `"gauss classic"` → `GAUSSCLASSIC` (spaces stripped) |
| `ParseWithUnderscores` | `"gauss_classic"` → `GAUSSCLASSIC` (underscores stripped) |
| `ParseWithEszett` | `"gaußclassic"` → `GAUSSCLASSIC` (ß → ss) |

### 1.3 string_to_integration_scheme `[semantic]` (BELFEM_ERROR — always active)

| Test Name | What It Verifies |
|---|---|
| `ParseUnknownStringThrows` | `"foo"` → `BELFEM_ERROR` fires |

---

## 2. Dispatcher Routing

**File:** `test_IntegrationDispatch.cpp`

### 2.1 GAUSS Scheme — Point Count Verification `[semantic]`

For each `(geometry, order)`, verify that `intpoints(GAUSS, geometry, order, w, p)` returns the expected number of points. This tests the dispatcher switch logic without checking point values.

| Geometry | Order | Expected Points | Rule Source |
|---|---|---|---|
| LINE | 2 | 2 | Gauss 1D (ceil(order/2)+1) |
| LINE | 4 | 3 | Gauss 1D |
| LINE | 8 | 5 | Gauss 1D |
| QUAD | 4 | 8 | gauss_quad8 (optimal) |
| QUAD | 6 | 12 | gauss_quad12 (optimal) |
| QUAD | 3 | 4 | gauss_quad (tensor-product 2×2) |
| HEX | 3 | 6 | gauss_hex6 (Hammer & Stroud) |
| HEX | 5 | 14 | gauss_hex14 |
| HEX | 7 | 34 | gauss_hex34 |
| HEX | 8 | 64 | gauss_hex (tensor-product 4³ fallback) |
| TRI | 1 | 1 | gauss_tri1 |
| TRI | 2 | 3 | gauss_tri3 |
| TRI | 5 | 7 | gauss_tri7 |
| TRI | 10 | 25 | gauss_tri25 |
| TRI | 20 | 79 | gauss_tri79 |
| TET | 1 | 1 | gauss_tet1 |
| TET | 2 | 4 | gauss_tet4 |
| TET | 5 | 15 | Keast |
| TET | 7 | 35 | Shunn & Ham |
| TET | 10 | 81 | Witherden & Vincent |
| TET | 14 | 236 | Zhang, Cui & Liu |
| PENTA | 1 | 1 | gauss_penta1 |
| PENTA | 5 | 16 | gauss_penta16 |
| PENTA | 9 | 71 | gauss_penta71 |
| PYRA | 1 | 1 | gauss_pyra1 |
| PYRA | 3 | 8 | gauss_pyra8 |
| PYRA | 5 | 27 | gauss_pyra27 |

### 2.2 GAUSSCLASSIC Scheme `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `GaussClassicQuadIsTensorProduct` | `intpoints(GAUSSCLASSIC, QUAD, 4, ...)` returns 9 points (3×3 tensor product) |
| `GaussClassicHexIsTensorProduct` | `intpoints(GAUSSCLASSIC, HEX, 4, ...)` returns 27 points (3³) |
| `GaussClassicTriFallsBackToGauss` | `intpoints(GAUSSCLASSIC, TRI, 5, ...)` returns same as `intpoints(GAUSS, TRI, 5, ...)` |
| `GaussClassicTetFallsBackToGauss` | Same for TET |

### 2.3 Lobatto Scheme `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `LobattoLineReturnsCorrectCount` | `intpoints(LOBATTO, LINE, 4, ...)` returns correct number of points |
| `LobattoIncludesEndpoints` | First and last points are ξ = -1 and ξ = +1 |

### 2.4 PYRA Regression (BUG-Q1) `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `PyraOrder1DoesNotError` | `intpoints(GAUSS, PYRA, 1, w, p)` completes without error |
| `PyraOrder3DoesNotError` | `intpoints(GAUSS, PYRA, 3, w, p)` completes without error |
| `PyraOrder9DoesNotError` | `intpoints(GAUSS, PYRA, 9, w, p)` completes without error |

### 2.5 Output Dimensions `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `LinePointsHave1Row` | Points matrix is 1 × N for LINE geometry |
| `QuadPointsHave2Rows` | Points matrix is 2 × N |
| `HexPointsHave3Rows` | Points matrix is 3 × N |
| `TriPointsHave3Rows` | Points matrix is 3 × N (barycentric) |
| `TetPointsHave4Rows` | Points matrix is 4 × N (barycentric) |
| `PentaPointsHave3Rows` | Points matrix is 3 × N |
| `PyraPointsHave3Rows` | Points matrix is 3 × N |
| `WeightsLengthMatchesPointsCols` | `weights.length() == points.n_cols()` for all geometries |

### 2.6 Dispatcher Error Paths `[semantic]` (BELFEM_ERROR — always active)

| Test Name | What It Verifies |
|---|---|
| `GaussTriInvalidOrderThrows` | `intpoints(GAUSS, TRI, 21, ...)` → error |
| `GaussTetInvalidOrderThrows` | `intpoints(GAUSS, TET, 15, ...)` → error |
| `GaussPentaInvalidOrderThrows` | `intpoints(GAUSS, PENTA, 10, ...)` → error |
| `LobattoNonLineThrows` | `intpoints(LOBATTO, QUAD, 4, ...)` → error |

---

## 3. auto_integration_order

**File:** `test_IntegrationDispatch.cpp`

### 3.1 Tests `[semantic]`

| Test Name | Input ElementType | Expected Order | Rationale |
|---|---|---|---|
| `AutoOrderConstant` | Any CONSTANT element | 1 | Minimal rule |
| `AutoOrderLinearTri3` | TRI3 | 4 | LINEAR → 4 |
| `AutoOrderLinearHex8` | HEX8 | 4 | LINEAR → 4 |
| `AutoOrderQuadraticTri6` | TRI6 | 7 | QUADRATIC → 7 |
| `AutoOrderQuadraticTet10` | TET10 | 7 | QUADRATIC → 7 |
| `AutoOrderCubicTet20` | TET20 | 10 | CUBIC → 10 |
| `AutoOrderQuarticTri15` | TRI15 | 13 | QUARTIC → 13 |

---

## 4. Quadrature Rule Invariants (Parameterized)

**File:** `test_IntegrationInvariants.cpp`

These tests apply to **every** supported `(geometry, order)` pair. Implement as value-parameterized tests.

### 4.1 Weight Sum = Reference Volume `[semantic]`

The single most important invariant. For every rule:

```
Σ wᵢ = V_ref
```

where V_ref is the reference domain volume (see table above).

| Test Name | What It Verifies |
|---|---|
| `WeightSumEqualsReferenceVolume` | For every `(geometry, order)`: `sum(weights) ≈ V_ref` within `tEps` |

This should sweep over ALL supported orders for each geometry: LINE orders 1–20, TRI orders 0–20, TET orders 0–14, QUAD orders 1–13, HEX orders 1–7, PENTA orders 0–9, PYRA orders 0–9 (pending BUG-Q1 fix).

### 4.2 Point Containment `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `LinePointsInDomain` | All points ξ ∈ [-1, 1] |
| `QuadPointsInDomain` | All points (ξ, η) ∈ [-1, 1]² |
| `HexPointsInDomain` | All points (ξ, η, ζ) ∈ [-1, 1]³ |
| `TriPointsInSimplex` | All barycentric coordinates ≥ 0 and ≤ 1, and L₁ + L₂ + L₃ ≈ 1 |
| `TetPointsInSimplex` | All barycentric coordinates ≥ 0 and ≤ 1, and Σ Lᵢ ≈ 1 |
| `PentaPointsInDomain` | ξ, η ≥ 0, ξ+η ≤ 1, ζ ∈ [-1, 1] |
| `PyraPointsInDomain` | Points within valid pyramid reference domain |

### 4.3 Positive Weights (Most Rules) `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `WeightsNonNegativeForLowOrder` | For orders ≤ 5 in all geometries: all weights > 0 |

**Note:** Some high-order rules (Keast tet5, tet11) have negative center weights. This is mathematically valid but should be documented, not tested as an error. Only low-order rules should be checked for positivity.

### 4.4 1D Gauss Symmetry `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `GaussLineSymmetric` | Points are symmetric about ξ = 0: ξᵢ = -ξₙ₋₁₋ᵢ, wᵢ = wₙ₋₁₋ᵢ |
| `LobattoLineSymmetric` | Same symmetry for Lobatto rules |

---

## 5. Polynomial Exactness

**File:** `test_IntegrationInvariants.cpp`

A rule of stated order p must integrate all monomials of degree ≤ p exactly over the reference domain.

### 5.1 1D Polynomial Exactness `[semantic]`

| Test Name | What It Verifies |
|---|---|
| `GaussLineIntegratesConstant` | Σ wᵢ · 1 = 2 (for any order) |
| `GaussLineIntegratesXExact` | Σ wᵢ · ξᵢ = 0 (odd function on symmetric domain) |
| `GaussLineIntegratesX2` | Σ wᵢ · ξᵢ² = 2/3 (for order ≥ 2) |
| `GaussLineIntegratesX4` | Σ wᵢ · ξᵢ⁴ = 2/5 (for order ≥ 4) |
| `GaussLineIntegratesX2p` | For each line rule of order p: ξ^p integrates exactly, ξ^(p+1) may not |

### 5.2 TRI Polynomial Exactness `[semantic]`

Reference integral formulas over the standard triangle (L₁, L₂, L₃ barycentric, area = 1/2):

∫ L₁ᵃ L₂ᵇ L₃ᶜ dA = a! b! c! / (a+b+c+2)! × 2

| Test Name | What It Verifies |
|---|---|
| `TriIntegratesConstant` | Σ wᵢ = 0.5 |
| `TriIntegratesL1` | Σ wᵢ · L₁ᵢ = 1/6 (for order ≥ 1) |
| `TriIntegratesL1Squared` | Σ wᵢ · L₁ᵢ² = 1/12 (for order ≥ 2) |
| `TriIntegratesL1L2` | Σ wᵢ · L₁ᵢ · L₂ᵢ = 1/24 (for order ≥ 2) |
| `TriIntegratesL1Cubed` | Σ wᵢ · L₁ᵢ³ = 1/20 (for order ≥ 3) |

### 5.3 TET Polynomial Exactness `[semantic]`

Reference integral over the standard tetrahedron (barycentric, volume = 1/6):

∫ L₁ᵃ L₂ᵇ L₃ᶜ L₄ᵈ dV = a! b! c! d! / (a+b+c+d+3)! × 6

| Test Name | What It Verifies |
|---|---|
| `TetIntegratesConstant` | Σ wᵢ = 1/6 |
| `TetIntegratesL1` | Σ wᵢ · L₁ᵢ = 1/24 (for order ≥ 1) |
| `TetIntegratesL1L2` | Σ wᵢ · L₁ᵢ · L₂ᵢ = 1/120 (for order ≥ 2) |
| `TetIntegratesL1Squared` | Σ wᵢ · L₁ᵢ² = 1/60 (for order ≥ 2) |
| `TetIntegratesL1Cubed` | Σ wᵢ · L₁ᵢ³ = 1/120 (for order ≥ 3) |

### 5.4 QUAD/HEX Polynomial Exactness `[semantic]`

Reference integrals over [-1,1]ᴰ:

∫ ξᵃ dξ = 2/(a+1) for even a, 0 for odd a.

| Test Name | What It Verifies |
|---|---|
| `QuadIntegratesConstant` | Σ wᵢ = 4 |
| `QuadIntegratesXiEta` | Σ wᵢ · ξᵢ · ηᵢ = 0 (for order ≥ 2) |
| `QuadIntegratesXi2` | Σ wᵢ · ξᵢ² = 4/3 (for order ≥ 2) |
| `HexIntegratesConstant` | Σ wᵢ = 8 |
| `HexIntegratesXi2Eta2` | Σ wᵢ · ξᵢ² · ηᵢ² = 8/9 (for order ≥ 4) |

### 5.5 PENTA Polynomial Exactness `[semantic]`

Penta reference domain is tri-base × line. Volume = 0.5 × 2 = 1.0.

| Test Name | What It Verifies |
|---|---|
| `PentaIntegratesConstant` | Σ wᵢ = 1.0 |
| `PentaIntegratesZeta` | Σ wᵢ · ζᵢ = 0 (odd function in ζ direction) |
| `PentaIntegratesXi` | Σ wᵢ · ξᵢ = 1/3 (= ∫ξ over tri [1/6] × ∫1 over line [2] = 1/3. **CORRECTED 2026-03-23:** was 1/6, verified with 1-point rule: w=1.0, xi=1/3) |

### 5.6 PYRA Polynomial Exactness `[semantic]` (after BUG-Q1 fix)

| Test Name | What It Verifies |
|---|---|
| `PyraIntegratesConstant` | Σ wᵢ = 4/3 |

---

## 6. End-to-End Integration with Shape Functions

**File:** `test_IntegrationEndToEnd.cpp`

**Entry point:** BELFEM provides gmsh files in `./more/gmsh/` for each element type. The elements are defined so that the center of mass is at the origin and the volume is 1.0. This enables a powerful end-to-end test: load a unit-volume element, set up shape functions + integration rule, compute ∫1 dV via numerical quadrature, and verify the result equals 1.0.

### 6.1 Volume Integration `[semantic]`

For each element type with a unit-volume gmsh mesh:

| Test Name | What It Verifies |
|---|---|
| `UnitVolumeIntegrationTRI3` | ∫1 dV = 1.0 using TRI3 shape functions + auto-order quadrature |
| `UnitVolumeIntegrationTRI6` | Same for TRI6 |
| `UnitVolumeIntegrationQUAD4` | Same for QUAD4 |
| `UnitVolumeIntegrationQUAD9` | Same for QUAD9 |
| `UnitVolumeIntegrationTET4` | Same for TET4 |
| `UnitVolumeIntegrationTET10` | Same for TET10 |
| `UnitVolumeIntegrationHEX8` | Same for HEX8 |
| `UnitVolumeIntegrationHEX27` | Same for HEX27 |
| `UnitVolumeIntegrationPENTA6` | Same for PENTA6 |

**Implementation pattern:**

```
1. Read unit-volume mesh from gmsh file
2. Extract node coordinates for the single element
3. Create InterpolationFunction for that element type
4. Get integration points via intpoints(GAUSS, geometry, auto_order, w, p)
5. At each integration point:
   a. Evaluate shape functions N and dNdXi
   b. Compute Jacobian J from dNdXi and node coordinates
   c. Compute det(J)
   d. Accumulate: volume += w(k) * abs(det(J))
6. EXPECT_NEAR(volume, 1.0, tTol)
```

**Note:** This test depends on mesh I/O (GmshReader). If unavailable, build the node coordinates manually from the gmsh file specifications. For simple elements (TRI3, TET4, HEX8), hardcoding 3–8 node coordinates is straightforward.

### 6.2 Linear Function Integration `[semantic]`

For the same unit-volume elements, integrate f(x,y,z) = 1 + x (since the centroid is at the origin, ∫x dV = 0 for a symmetric element, so ∫f dV = 1.0).

| Test Name | What It Verifies |
|---|---|
| `LinearIntegrationQUAD4` | ∫(1+x) dV = 1.0 for unit-volume QUAD4 centered at origin |
| `LinearIntegrationHEX8` | Same for HEX8 |
| `LinearIntegrationTET4` | Same for TET4 |

---

## 7. What We Do NOT Test

- **Literal point/weight values** — data is copy-pasted from papers/quadpy; testing exact values is low ROI
- **Fortran routine internals** (`intpoints_gauss`, `intpoints_lobatto`) — tested indirectly through invariants
- **Performance benchmarks** — per strategy doc
- **All 70+ individual table headers** — tested through parameterized invariant sweeps
- **PYRA integration beyond smoke level** — pending BUG-Q1 fix, and author notes uncertainty about order mapping

---

## 8. Implementation Notes for Claude Code

1. **Parameterized test sweep is essential.** The weight-sum and containment tests should cover ALL supported `(geometry, order)` pairs in a single parameterized suite. Build a vector of `{GeometryType, order, expectedVolume, expectedPointRows}` structs.
2. **Reference volume constants:** LINE=2.0, QUAD=4.0, HEX=8.0, TRI=0.5, TET=1.0/6.0, PENTA=1.0, PYRA=4.0/3.0. These are exact rational numbers — compute them in test code, don't hardcode decimals.
3. **BUG-Q1 (PYRA break) will cause test failures.** If the test framework catches the `BELFEM_ERROR` as a throw, the PYRA tests will throw. Document this and mark the PYRA regression tests as expected-to-fail until the bug is fixed.
4. **The tet10 Shunn & Ham weight table looks suspicious.** The weight values appear to match coordinate values, not proper quadrature weights. The weight-sum test will immediately reveal if this is a problem (weights should sum to 1/6).
5. **Polynomial exactness formulas:** For simplex integrals, use the analytical formula: ∫ L₁ᵃ L₂ᵇ ... dV = (a! b! ... / (a+b+...+D)!) × D! × V_ref. Implement this as a helper function.
6. **GAUSSCLASSIC fallback:** The dispatcher explicitly resets `GAUSSCLASSIC` to `GAUSS` for non-QUAD/HEX geometries. Test this by verifying identical output.
7. **End-to-end tests may require mesh I/O.** If GmshReader is not available in the test binary, hardcode node coordinates for a few simple unit-volume elements (TRI3, TET4, HEX8). The gmsh files in `./more/gmsh/` define the exact geometry.
8. **1D symmetry test:** For `gauss_line`, the Fortran routine produces symmetric points. Verify ξᵢ = -ξₙ₋₁₋ᵢ and wᵢ = wₙ₋₁₋ᵢ for all 1D rules.
9. **Barycentric 4th coordinate for TET:** All TET rules compute `L₄ = 1 - L₁ - L₂ - L₃` in a loop. The containment test should verify that L₄ ≥ 0 and L₄ ≤ 1 as well.
10. **The tet5 and tet11 rules have negative center weights.** This is standard for higher-order Keast rules and is mathematically correct. Do not flag negative weights as errors for these rules.

---

## 9. Codex Audit Checklist

When reviewing Claude Code's test implementation, verify:

- [ ] Weight-sum test covers every supported `(geometry, order)` pair, not just a sample
- [ ] Reference volumes are computed as exact rational expressions, not hardcoded decimals
- [ ] Point containment tested for all geometry types
- [ ] Barycentric sum ≈ 1 verified for TRI and TET points
- [ ] PYRA regression test (BUG-Q1) is present and documents the missing-break issue
- [ ] GAUSSCLASSIC fallback tested: non-QUAD/HEX produces same result as GAUSS
- [ ] Polynomial exactness tested for at least LINE, TRI, TET with monomial integrals
- [ ] auto_integration_order mapping tested for at least 5 element types
- [ ] String parsing tests include ß normalization, spaces, underscores, case insensitivity
- [ ] Negative weights documented for Keast rules, not flagged as errors
- [ ] End-to-end volume integration tested for at least 3 element types
- [ ] `EXPECT_NEAR` for all floating-point comparisons
- [ ] BELFEM naming conventions (`t` prefix for locals)
