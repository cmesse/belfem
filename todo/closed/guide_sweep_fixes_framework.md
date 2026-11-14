# Guide-Sweep Code Defects — Critical Framework Part

**Date:** 2026-10-22
**Purpose:** Fix the code defects in critical framework code (`src/mesh`, `src/fem/kernel`, the graph
vertex base of every mesh entity) that the 2026-10-22 usage-guide jury sweep found
(`devlog/dl20261022_usage_guide_jury_sweep.md`). One change, and one commit, per defect.
**Module:** `src/mesh`, `src/fem/kernel`, `src/math/graph`
**AIs involved:** Claude (trace + plan), Codex `gpt-6-astra` (audit), Grok `grok-4.7` (audit)
**Status:** ✅ COMPLETE (2026-10-22) — D1–D10 landed in e55aee3, d0af680, 19f4411, 1e23a2b, 5053447 (D10, found by the follow-up jury), 9abf660 (cubic test), 89fe9e9, 4287445, 17eab8f, abe679b, 975f94f, a90d2cf. Reviewed, not verified: plan jury and code jury (Codex gpt-6-astra/high then gpt-5.6-terra/high after a usage-limit failure, Grok grok-4.7/high), `-fsyntax-only` on all 392 debug-tree TUs. **Owed gate (Christian):** full build incl. tests and nonfree (Vertex copy ban), `test_mesh` (TensorMeshBSpline.*), `test_fem` (IntegrationAuto.*), `examples/rlc_circuit` (`L1` column vs FEM current). O2 moved to `todo/tensor_mesh_control_point_coordinates.md`.

> **Scope guards:**
> - Only the defects listed here. The companion plan `todo/closed/guide_sweep_fixes_modules.md` holds the
>   non-framework items; the cubic-EoS phase question is a separate design plan
>   (`todo/cubic_eos_phase_selection.md`), not a fix.
> - No change to any usage guide (`src/**/doc/*usage_guide*.md`), which another session is editing.
> - Never touch the cohomology core (protocol §7.1).
> - Evidence so far: source trace only, plus one standalone METIS probe (scratchpad, outside the build). Nothing is built.

---

## 1. Defects

| ID | Defect | Evidence | Sev |
|---|---|---|---|
| D1 | `interpolation_order_numeric` has no arm for LINE5 or TRI15, so it hits `BELFEM_ERROR( false, "Unknown element type")`. The element factory (`cl_Element_Factory.cpp:79,95`), the Gmsh map (`meshtools.cpp:471,474`) and the interpolation factory (`cl_IF_InterpolationFunctionFactory.cpp:123,139`) all build both types | `meshtools.cpp:266-326` | HIGH |
| D2 | `interpolation_order()` has no arm for PYRA5/13/14, LINE5, TRI15 or TET35, so they return `UNDEFINED`. A linear PYRA5 therefore reads as not LINEAR: `Calculator::mIsLinear` is false (`cl_FEM_Calculator.cpp:1727,1756,1907`), and `OrderConverter::check_input_mesh` refuses it as "not linear" (`cl_Mesh_OrderConverter.cpp:101,113`) | `meshtools.cpp:210-262` | LOW |
| D3 | `geometry_type()` has no LINE5 arm, so it returns `GeometryType::UNDEFINED` for an element type the factory builds. Found while tracing D1, a sibling of the same omission | `meshtools.cpp:336-340` | LOW |
| D4 | `compute_tmatrix_quad16/hex64` evaluate the **quadratic** B-spline (`bspline_quad`), which fills 3 of 4 entries. The 4th stays at its initial 0, so the cubic T-matrix is wrong. `bspline_cub` exists but has no caller, and its output is **reversed**: `aResult(0) = (1+xi)^3/48` is the right-end basis function, while `bspline_quad(0) = (1-xi)^2/8` is the left-end one | `cl_TensorMeshFactory.cpp:995-1027,1086-1129,1134-1149` | HIGH (wrong geometry for a cubic tensor mesh) |
| D5 | `element_ijk` casts `( x - origin ) * invstep` to unsigned `index_t` **before** clamping. A point below the origin is a negative real → float-to-unsigned conversion out of range, which is UB ([conv.fpint]) | `cl_TensorMeshConfig.hpp:207-213` | LOW |
| D6 | `SolverData::assemble_jacobian` binds `SpMatrix & E = *mEnforcementMatrix;` and `I = *mImpositionMatrix;` without a condition. Both are null unless JEDI force was requested before `initialize()` (`cl_FEM_DofMgr_SolverData.cpp:361,373`), and they are used only under `mUseJediForce` (`:1539-1551`). Forming a reference to `*nullptr` is UB even when the reference is never used | `cl_FEM_DofMgr_SolverData.cpp:1480,1482` | LOW (UB, no observed failure) |
| D7 | `Controller::finalize` writes the circuit output row (`mCircuit->save_timestep()`) **before** `compute_circuit_current()` copies this step's FEM terminal current into the circuit (`FEMTwoTerminals::set_IV`, `cl_FEMTwoTerminals.cpp:56-61`). Every `CircuitResults` terminal column therefore reports the previous accepted step's current. Christian ruled 2026-10-22: reorder, with a jury check | `cl_FEM_Controller.cpp:3059-3065` | MED (output mislabelled by one step) |
| D8 | `graph::Vertex` owns `Vertex ** mVertices` (deleted in `~Vertex` via `reset_vertex_container()`) but declares no copy or move members, so an implicit copy double-frees. `Node`, `Edge`, `Face`, `Facet`, `Segment`, `Basis`, `Dof`, `ElectricNode` and `Component` all inherit it. `coding_philosophy.md`: "Any class owning raw memory declares or deletes all five special members" | `cl_Graph_Vertex.hpp:58-68`, `cl_Graph_Vertex.cpp:24-27` | LOW (latent) |

Retracted before planning:
- **`metis_ndp( g, 1 )` reaches `METIS_NodeNDP` — FALSE POSITIVE (probe, 2026-10-22).** A standalone
  probe against `/opt/scls/lib/libmetis` on a 20×20 grid graph: `METIS_NodeNDP` with npes=1 returns
  `METIS_OK`, and its permutation is identical to `METIS_NodeND`'s. Canaries past `sizes[2*npes-1]`
  are intact (also at npes = 2, 3, 4). BELFEM never reads `sizes` (`fn_Graph_METIS.cpp:155-174`).
  What remains is the commented-out redirect and its false comment "requires at least 2 partitions"
  (`:109-115`): that is dead code to delete (comment-only, no jury).

## 2. Steps

- [x] **R1 (D1, D2, D3) — complete the three element-type switches in `meshtools.cpp`.**
  `interpolation_order_numeric`: add `LINE5`, `TRI15` → 4 (beside the `TET35` arm).
  `interpolation_order`: `PYRA5` → LINEAR, `PYRA13` → SERENDIPITY, `PYRA14` → QUADRATIC (this matches
  the numeric function's 1/2/2 and the QUAD8-vs-QUAD9 split), `LINE5`, `TRI15`, `TET35` → QUARTIC.
  `geometry_type`: add `LINE5` → LINE.
  **What could break:** a PYRA5 block now passes `OrderConverter::check_input_mesh` and then fails
  later, in `upgrade_type`, with the vaguer "UNSUPPORTED ELEMENT TYPE" (`cl_Mesh_OrderConverter.cpp:1133-1137`).
  Mitigation inside R1: in `check_input_mesh`, refuse a block whose type `upgrade_type` cannot map,
  naming the type, before any work is done. `mIsLinear` becomes true for PYRA5 blocks; its only
  consumer asserts on it (`cl_FEM_Calculator.hpp:2707`), which is correct for a linear pyramid.
  The normal-vector dispatch (`cl_FEM_Calculator.cpp:1537`) switches on TRI and QUAD only, so PYRA is
  unaffected. **Gate:** compile + `make check-fast`; the mesh tests that read LINE5/TRI15 Gmsh files,
  if any exist (to be listed at the code jury).
- [x] **R2 (D4) — use the cubic B-spline in the cubic T-matrices.** Rewrite `bspline_cub` to return
  the four uniform cubic B-spline functions on ξ ∈ [-1, 1] in left-to-right order, the same order as
  `bspline_quad`: N0 = (1-ξ)³/48, N1 = (3ξ³-3ξ²-15ξ+23)/48 (currently `aResult(2)`),
  N2 = (-3ξ³-3ξ²+15ξ+23)/48 (currently `aResult(1)`), N3 = (1+ξ)³/48. Call it from
  `compute_tmatrix_quad16` and `compute_tmatrix_hex64`. Partition of unity: Σ = 48/48 at every ξ
  (checked at ξ = 0: 1+23+23+1). **Confirm at the jury:** that the control-point index
  `idx(·, j) = 0..3` runs left to right, as the quadratic case's does, so the order above is the
  right one. **Gate:** a unit check that each row of T sums to 1 and that T reproduces a linear field
  exactly; a cubic tensor mesh case if one exists.
- [x] **R3 (D5) — clamp in the real domain, then cast.**
  `const real tS = std::floor( ( x - mOrigin( i ) ) * mInvElementSteps( i ) );` then
  `return ( index_t ) std::clamp( tS, 0.0, ( real ) ( mNumElementsPerDim( i ) - 1 ) );`.
  Behaviour for x ≥ origin is unchanged (floor and truncation agree for non-negative values).
  **Gate:** compile.
- [x] **R4 (D6) — no references to null.** Bind `E` and `I` inside the `else if ( mUseJediForce )`
  branch, or use `SpMatrix *` with a `BELFEM_ASSERT( tE != nullptr && tI != nullptr )` at the top of
  that branch. Assembly results are unchanged. **Gate:** compile + `make check-fast`.
- [x] **R5 (D7) — compute the circuit current, then save the row.** Swap the two calls at
  `cl_FEM_Controller.cpp:3061-3063`. Why this is safe: `compute_circuit_current` only calls
  `set_current_and_voltage` → `FEMTwoTerminals::set_IV` (stores I, V, sets the component current).
  It runs on rank 0 in the same block, after the step's last solve, and it ran unconditionally right
  after `save_timestep` anyway. So the circuit state the next step starts from is the same; only the
  output row changes. **To confirm at the jury:** that `save_timestep` writes nothing else that
  `set_IV` changes, and that the restart dump (`save_state`, `:5236`) is unaffected.
  **Behaviour change:** the terminal columns of `CircuitResults` now report this step's current;
  belfem-85 must update the circuit guide. **Gate:** `examples/rlc_circuit`, comparing the column
  against the solver's terminal current at the same t.
- [x] **R6 (D8) — Rule of Five for `Vertex`.** Declare the copy and move constructors and assignment
  operators `= delete` in `cl_Graph_Vertex.hpp`. A grep for copy construction of every subclass finds
  none (only `Basis( const Basis & ) = delete`). **What could break:** any by-value copy anywhere in
  `src/`, `tests/` or `nonfree/` becomes a compile error. That is the point, but the compile is the
  gate. **Gate:** full build, including tests and nonfree.
- [x] **R7 (comment only, no jury)** — delete the dead NodeNDP redirect and its false comment
  (`fn_Graph_METIS.cpp:109-115`). Landed beside R6.

## 3. Open questions

- **O1** — PYRA13 as SERENDIPITY: it matches the existing enum's QUAD8/HEX20/PENTA15 grouping, but
  the enum's consumers (`cl_FEM_Calculator.cpp:1537` and others) never test SERENDIPITY for a pyramid.
  Any objection?

## 4. Amendment after plan jury round 1 (2026-10-22)

Record: `tmp/ai_exchange/review_guide_sweep_fixes_framework.md` (verification + reconciliation).
Every change below is what the auditors asked for; there was no reviewer disagreement on a
load-bearing point, so no round 2 on the plan.

- **D2 severity raised to HIGH.** `auto_integration_order` aborts with `BELFEM_ERROR` on
  `UNDEFINED` (`fn_intpoints_auto_integration_order.cpp:49-53`), so every pyramid and TET35 dies in
  automatic integration-order selection. A straight PYRA5 also starts with `mCurvedFlag` set
  (`cl_ElementTemplate.hpp:445-447`); after R1 it does not, which is correct for a linear element.
- **R1 gate:** enable `IntegrationAuto.DISABLED_AutoOrderCoversPyramidsAndQuarticTets`
  (`tests/fem/test_IntegrationAuto.cpp:213`), extend it with LINE5 → 13 and TRI15 → 13, and add an
  `interpolation_order_numeric` check for LINE5/TRI15. The early `OrderConverter` refusal is a
  named `BELFEM_ERROR` listing the element type; it does not call `upgrade_type`.
  One commit for `meshtools.cpp` + test; the `OrderConverter` refusal is its own commit.
- **D4 corrected:** the 4th basis entry is uninitialized, not 0 (`Vector( n )` does not fill,
  `cl_AR_Vector.hpp:55-56`).
- **New D9 (HIGH; Codex, Grok, Claude):** `link_bsplines_quad16` puts `(p, q+1)` and `(p, q+2)` into
  slots 4 and 5 (`cl_TensorMeshFactory.cpp:784-785`); the T-matrix expects `(p+1, q)` and
  `(p+2, q)` (`idx`, `:999-1002`). The bottom-edge control points are never linked, so the projector
  matrix has empty rows for them. quad9, hex27 and hex64 match their tables (script check).
  **R2b:** fix the two slots. Own commit.
- **R2 gate replaced:** the linear-geometry gate was wrong, since control-point coordinates are not
  Greville points and never enter the projection (`cl_DatabaseProjector.cpp:216-300`). New gate:
  every row of T sums to 1, every control point is linked by at least one element, and projecting a
  cubic polynomial field reproduces it at the nodes. **Citation** for the basis:
  Zienkiewicz et al. 2013, §17.2.1, Eq. 17.4 (Cox–de Boor), with t = (ξ + 1)/2.
- **R3:** include `<cmath>`; the D5 claim is narrowed: UB arises when the truncated value does not
  fit in `index_t` (negative ≤ −1, or huge positive). Finite in-range values behave as before.
- **R5:** the change affects the **terminal-current** columns only; lumped columns and node voltages
  are unchanged. The gate compares the `L1` column of `examples/rlc_circuit` with the FEM terminal
  current at the same t. Grok confirmed `save_memdump` runs after `finalize` (`belfem.cpp:243-248`).
- **R7** lands as its own comment-only commit.
- **New O2 (Christian):** control points are placed at node spacing (`step()` = node step,
  `cl_TensorMeshConfig.hpp:70-73`) in every tensor factory, not at Greville points. Their
  coordinates are not used by the projector, so this is cosmetic today. Intended?
