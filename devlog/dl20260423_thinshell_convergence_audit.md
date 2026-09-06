# Thin-Shell H-φ Convergence Audit — Read-Only Investigation

**Date:** 2026-04-23
**Purpose:** Exhaustive audit of the thin-shell refactoring changes on branch `ghost` since commit `b27eee6`, triggered by the observation that `test_fem` (64/64) and `test_mesh` (46/46) pass but the full H-φ thin-shell Maxwell simulation does not converge. No source edits made — see `todo/ai_exchange.md` for Codex's parallel observations.
**Mode:** Read-only. Parallel sub-agent audits consolidated here.
**Audience:** Next session continuing the thin-shell Option-B debug.

---

## 1. Scope

Reviewed the full `git diff HEAD -- src/` (≈ 1500 LOC across 57 files). Parallel sub-agents covered:

1. Element topology files (PENTA6/TS, PENTA18/TS, HEX8/TS, HEX20, HEX27, HEX64, QUAD4TS, QUAD9TS, PENTA15).
2. Integration-points dispatch (fn_IF_initialize_integration_points_on_facet + SideSet::initialize + Calculator::allocate).
3. HEX8 Nédélec basis sign fix (F(0..3) negation vs new edge convention).
4. Calculator refactor (normal dispatch, slave index scheme, mEdgeFunctionsSlave sizing, all link() paths).
5. meshtools helpers and MaxwellFactory hang_thinshell_* consumers.

Personal re-verification of the most suspicious leads. Conclusions below.

## 2. What the audit **cleared**

All of these are ruled out as the convergence culprit for the current linear 3D PENTA6TS thin-shell run:

- **HEX8 Nédélec basis** (`cl_EF_HEX8.cpp`). All 12 F-expressions integrate to +1 under the new HEX8 edge convention; all dF rows match the F derivatives; HEX8 and HEX8TS agree on the 8 loop edges they share. Confidence high. The sign flip was correct.
- **`Calculator::normal_penta` dispatch**. Both volume PENTA and PENTA\*TS now route through `normal_penta`; case 3 (bottom) and case 4 (top) are correctly negatives of each other. No stale `normal_penta_ts` references remain.
- **Slave integration index refactor** (`slave_integration_index_{2d,tet,hex,penta}` returning flat `uint`). Flat index is computed once in `Calculator::link(Element*)` and used for both `mSlaveIntegration` and `mEdgeFunctionsSlave`. `mEdgeFunctionsSlave` is sized by total permutations (sum of orientations), matching `SideSet::mSlaveIntegration`. Loop structures at `cl_FEM_Calculator.cpp:314-328` and `cl_FEM_SideSet.cpp:515-534` use the same `for f { for o }` pattern producing the same flat order.
- **Orientation table merge for PENTA6/PENTA6TS, PENTA18/PENTA18TS, PENTA15**. Only consumed by `test_facets.cpp`. No runtime FEM consumer indexes the table directly. Not a convergence-path concern.
- **`intpoints_penta` master + slave** handles all 5 facets. For facets 3/4 (tris) the slave uses `tCase = (aSlaveIndex - 3) * 3 + aOrientation` with 6 cases (3 tri orientations × 2 tri facets). Slot 12 = (facet 3, orient 0) gives `(xi, eta, -1)`. Slot 15 = (facet 4, orient 0) gives `(xi, eta, +1)`. Integration points line up physically with master's top-tri at `(xi, eta, +1)`.
- **Ghost-facet orientation=1 hardcode** (`cl_ThinShellFactory.cpp:1852`). `compute_orientation` would compute exactly this value under natural extrusion: slave's first corner (shell-block-b.node[0]) has the same `original()->id()` as master's facet-4 position 0 (shell-block-(b-1).node[3]), so orient=1. The hardcode skips the work but yields the identical result. For TRI3, `to_master_orientation` case 1 is the reflection `[0→0, 1→2, 2→1]` — exactly what the CW slave bottom → CCW master top needs.
- **`hang_thinshell_edges_on_{nodes,edges}_{bottom,top}` helpers**. After replacing `get_nodes_of_facet(bottom_idx)` / `get_edges_of_facet(bottom_idx)` with `mesh::get_bottom_nodes` / `mesh::get_bottom_edges`, the shell-bottom node/edge lists come back in natural label order (`node(0..2)` / `edge(0..2)` for PENTA6TS). This matches the volume master's top-facet CCW ordering position-wise. Sequential `set_index(k)` works. Top path (`get_top_nodes` / `get_top_edges`) was never broken.
- **`Calculator::link(mesh::Facet*)`** with hardcoded `index_on_slave()*3+orientation-1`. This overload is only called by `Maxwell_TMatrix::process`, which is TET10-only, and `MaxwellFactory::create_hanging_edges_and_facets` at line 1100 has a `BELFEM_ERROR( max_element_order == 1, ... )` gate — the entire function only runs for linear meshes, so the TMatrix block at line 1297 (guarded by `max_element_order == 2`) is dead for our failing case. Not relevant.
- **TOP path `to_master_orientation`** on the slave side. After this reorder, `tNodesOnVolume` is in master's (air-below's top-facet) CCW order. `get_top_nodes(aElement)` returns shell-top-block.node[3..5] in natural CCW order. Matches position-wise.

## 3. Definite bugs (latent, not triggering linear 3D failure)

**None of these will fire on the current failing case**, but they will bite future runs and should be fixed before quadratic or 2D-shell work resumes:

### 3.1 `meshtools.cpp:1284` — `get_bottom_nodes(HEX20)` reads `node(21)`
HEX20 has only nodes 0..19 (`cl_Element_HEX20.hpp` template `<20,8,12,6,6>`). `aElement->node(21)` is out-of-range. Also `aNodes.set_size(9, nullptr)` is wrong: HEX20 bottom face is QUAD8, not QUAD9. No HEX20 thin shell exists in production, so this is latent.

### 3.2 `meshtools.cpp:1237` — `get_bottom_nodes(QUAD9/QUAD9TS)` uses `node(5)`
Per `cl_Element_QUAD9TS.hpp:62-68`, the bottom-facet midside is `mNodes[4]`, not `mNodes[5]` (node 5 is the right-side midside). Latent — affects only quadratic 2D shells.

### 3.3 `cl_Element_QUAD9TS.hpp` — missing `get_nodes_of_edge` specialization
The template `<9,4,3,4,1>` declares E=3 edges but no `get_nodes_of_edge` override. Default falls through to `get_nodes_of_facet`, which silently returns the wrong nodes for anything querying edge topology on QUAD9TS. Latent until quadratic 2D shells are run through the Nédélec pipeline.

### 3.4 `cl_Element_HEX64.hpp` — pre-existing `get_edges_of_facet` inconsistency
Sub-agent flagged 5 of 6 HEX64 facets returning edges that contradict their node lists. Verified by `git diff` that our changes to HEX64 this round were only the cosmetic `case( 0 )` → `case 0` rewrap plus the edge-direction flips on the bottom loop consistent with HEX8/20/27. The `get_edges_of_facet` inconsistency is pre-existing (predates this branch). HEX64 is used only in `TensorMeshConfig` for cubic mesh generation — not in thin-shell Maxwell paths.

### 3.5 Stale comments in `cl_MaxwellFactory.cpp:1362-1368` and `:1483-1486`
The comments describe reordering the shell-bottom face to hide the "CW for outward-down normal" asymmetry. After switching to `mesh::get_bottom_*` helpers, this reordering no longer happens at the call site (the helpers bypass `get_nodes_of_facet(bottom)` entirely). Cosmetic.

## 4. The actual convergence hunt: where the real bug is likely to live

With the obvious leads ruled out, the remaining space is:

### 4.1 MOST LIKELY — Mesh-level coupling in a path we didn't refactor

Things `test_fem` / `test_mesh` do **not** exercise but the full hphi run does:

1. **`Mesh::update_facet_nodes`** (`cl_Mesh.cpp:694`) runs during mesh finalization for every sideset facet, calling `tFacet->master()->get_nodes_of_facet(tFacet->index_on_master(), tNodes)` and copying the result into the facet's own element node container. For ghost facets, master is a PENTA6TS with `index_on_master = top_facet_index = 4`. That returns `[s.3, s.4, s.5]` in CCW — same as before — so ghost-facet nodes are unchanged from pre-Option-B. Likely clean.

2. **`Facet::compute_orientation`** uses `get_corner_nodes_of_facet` for both master and slave. For PENTA6TS slave at ghost facet (facet 3, nodes `[s.0, s.2, s.1]` CW), compared against master's facet 4 `[m.3, m.4, m.5]` CCW. `compute_orientation` finds orient by matching `slave.first.original()->id()` against master's facet list. If the `original()` chain is set up consistently (which it is for natural extrusion), orient = 1. The hardcoded orient at `set_slave(.., 3, 1)` skips `compute_orientation`, so this is not where the bug lives — **but check whether any non-ghost shell-air sideset facet gets `compute_orientation` called with the new CW bottom ordering** and produces a different orient than it did pre-Option-B. If orientation_on_slave changes, downstream `to_master_orientation` behaves differently.

3. **`FaceFactory`** — whether face-level (not facet-level) connectivity for the shell uses the facet ordering. For a thin-shell mesh with PENTA18TS (F=3 faces for quadratic) this matters more; for PENTA6TS (F=1) the single face is the midsurface and less likely broken.

### 4.2 SECOND MOST LIKELY — Cut / homology setup on thin shells

`CutFactory::create_thin_shell_cuts` (`cl_CutFactory.cpp`) duplicates thin-shell interface nodes and relinks slave-side volume elements. It uses `tFacet->master()->get_nodes_of_facet(tFacet->index_on_master(), tNodes)` and may implicitly assume something about master-facet orientation. With CW/CCW for PENTA6TS facets 3 and 4 matching volume PENTA6 exactly (Option-B canonicalization), the CutFactory should see the same nodes in the same order as before. Worth a targeted read-only trace if nothing else pans out.

### 4.3 THIRD — H-φ interface coupling on the actual shell-air sideset

For the shell-air sideset facet (not ghost), the master is typically the volume PENTA6 below the shell, slave is volume PENTA6 above. Before Option B, the shell element's `get_edges_of_facet(0)` gave `[mEdges[0..2]]` in CCW walking order for use by `hang_thinshell_edges_*`. After Option B, `get_bottom_edges(aElement)` gives the same `[mEdges[0..2]]`. **Verified equivalent.** The hanging-edge coupling should reproduce the pre-Option-B behavior. **If it doesn't, the bug is likely outside this path.**

### 4.4 FOURTH — Ghost facet sideset's own integration table dimensions

Pre-Option-B, ghost facet SideSet for PENTA6TS master/slave had 2 master facets × 3 orientations = 6 permutations. After: 5 master facets × {4,4,4,3,3} = 18 permutations. The slot used by ghost (`index_on_slave=3, orient=1`) is now flat slot 12 (offset table `{0,4,8,12,15}`). Pre-Option-B used flat slot 0 or 3 (for old `index=0, orient=1` with `*3` offset). The CODE uses `slave_integration_index_penta` to compute the right slot, but any cached/stored index from an earlier pass could be stale. **Search for anything that caches a slave-integration slot number and compare against the new flat layout.**

## 5. Suggested next-step investigation plan

1. **Before touching any code**, run the failing hphi test with prints added (instrumentation, not code changes) at these checkpoints and compare against a known-good pre-Option-B run if available:
   - `Facet::compute_orientation` output (`mOrientationOnSlave` value) for every thin-shell-sideset facet.
   - `Calculator::link(Element*)` flat slave index at `tSlaveIntIndex` for each ghost-facet element.
   - The first/last few `(xi, eta, zeta)` integration points produced by `intpoints_penta` for master_index=4 and slave (3, 0) — verify they're physically at the midsurface interface.
   - Hanging-edge source list for one shell-bottom edge and one shell-top edge — verify they point at the physically correct air nodes.

2. **Bisect** by temporarily replacing individual Option-B changes with their pre-Option-B equivalents (in a throwaway branch) to identify which single change breaks convergence:
   - (a) PENTA6TS facet renumbering (2 → 5 facets).
   - (b) PENTA6TS facet 3 CW node ordering.
   - (c) Ghost facet `set_slave(.., 3, 1)` vs old `set_slave(.., 0, 1)`.
   - (d) `get_bottom_nodes`/`get_bottom_edges` helper vs raw `get_nodes_of_facet(3)`/`get_edges_of_facet(3)` fallback.
   - (e) Orientation table merge for PENTA6 + PENTA6TS.

3. **If (d) turns out to be the issue**, the permutation the helpers apply may still be wrong for a specific mesh topology we haven't considered (e.g., non-natural extrusion, periodic boundary, cut duplicates). In that case, the robust fix is the physical-node-ID matching pattern from `hang_thinshell_edges_on_edges_top` (ID → position map via `original()->id()`).

4. **Low-cost sanity check**: build and run with `-DBELFEM_DEBUG` and see if any assertion fires that test_fem/test_mesh don't exercise.

## 6. Reference: Codex's parallel observations (from `todo/ai_exchange.md` and in-session review)

- Codex flagged `Calculator::link(mesh::Facet*)` as containing a stale `*3` hardcode; re-verified here and confirmed it's TET10-only via TMatrix, gated behind `max_element_order == 2`. Not the culprit.
- Codex flagged the `orientation=1` hardcode as potentially wrong for dihedral flips; re-verified here and confirmed it's correct for TRI3 under BELFEM's orientation convention (`to_master_orientation` case 1 IS the reflection, not a rotation).
- Codex flagged the HEX20/QUAD9TS helper bugs (§3.1, §3.2); confirmed. Latent.

## 7. What I am confident did NOT break

- HEX8 Nédélec basis (unit circulation verified analytically, dF consistent, HEX8/HEX8TS agree on shared edges).
- PENTA6TS Nédélec basis (`EF_PENTA6TS.cpp`) — not modified in this batch at all.
- PENTA6TS edge numbering and directions (6 edges, unchanged from pre-Option-B).
- `normal_penta` all 5 facet cases (Jacobian-based, independent of node ordering).
- Slave integration indexing and edge-function indexing (single flat index path, verified consistent).
- `hang_thinshell_edges_on_{nodes,edges}_{bottom,top}` position-wise coupling (traced manually for PENTA6TS / TRI3 surface, aligned).
- β convention in HTS kernels (`std::abs(dot(n,b)) / norm_b`, unsigned tape-axis angle, `Ic(θ)` symmetry preserved).

## 8. Open question for the user

The failing simulation — is it single-layer thin shell (elementsBottom = elementsTop) or multi-layer (stacked PENTA6TS blocks with ghost facets in between)?

- If **single-layer**: only `hang_thinshell_edges_*` is exercised; no ghost facets; the bug is in the shell-air interface coupling somewhere we haven't traced.
- If **multi-layer**: ghost facets are exercised too; the bug could be in the new 18-permutation flat-index scheme or in the `h_ghost` kernel consuming Em/Es at slot 12.

The answer narrows §4 to either §4.3 (single-layer) or §4.4 (multi-layer).

---

## Files reviewed (summary)

Topology: all `cl_Element_*.hpp` listed in §1.1.
Integration: `fn_IF_initialize_integration_points_on_facet.{hpp,cpp}`, `cl_IF_IntegrationData.{hpp,cpp}`, `cl_FEM_SideSet.hpp`.
Calculator: `cl_FEM_Calculator.{hpp,cpp}` lines 290-330, 770-900, 950-1100, 1759-1912, 2045-2090.
MaxwellFactory: `cl_MaxwellFactory.cpp` lines 1080-1700.
ThinShellFactory: `cl_ThinShellFactory.cpp` lines 1800-1860.
Kernel: `mt_maxwell_h.cpp` lines 2020-2096 (`h_ghost`).
Mesh helpers: `meshtools.{hpp,cpp}` with the new helpers §3.

All clean except the latent bugs enumerated in §3.

---

## Addendum — Codex's finding, missed in §1-7 above

After sharing this audit with Codex, Codex identified the actual primary regression: **two default-value flips in `src/mesh/cl_ThinShellFactory.hpp:65-66`** that my audit did not surface because I was tracing call paths rather than comparing default states.

### Addendum.1 — `mCreateGhostFacets = true` → `false` (primary suspect, high confidence)

`src/mesh/cl_ThinShellFactory.hpp:66`. With this off, the following cascade disappears:

- `cl_ThinShellFactory.cpp:147-160` — no `tLayers(tCount)->hasDuplicates = true` at material interfaces, so no node/edge duplication.
- `cl_ThinShellFactory.cpp:244` — `create_ghost_facets()` never called.
- `cl_ThinShellFactory.cpp:255-260` — `DomainType::Ghost` sideset never created (`g == 0`).
- `cl_IWG_Maxwell.cpp:572` — `DomainType::Ghost` → `maxwell::h_ghost` dispatch never fires.
- `mt_maxwell_h.cpp:1926` — the `h_ghost` Nitsche interior-penalty kernel that enforces H-field interface continuity across material boundaries never assembles.

For multi-layer tape (Cu/Ag/YBCO/…), adjacent layers of different materials now share the same Nédélec edge DOFs instead of being coupled through the interior penalty — exactly the symptom of "non-physical currents on top and bottom layers." This is the regression.

### Addendum.2 — `mConnectorsForAllLayers = false` → `true` (intentional per inline comment, but fold coupling incomplete)

`src/mesh/cl_ThinShellFactory.hpp:65`. Inline comment says "must always be true for the penalty in binomial direction." Per `src/fem/maxwell/doc/shell_connector_coupling.md` and Codex's review, the connector side still has:

- `cl_ThinShellFactory.cpp:2034` — `create_vertical_edges` commented out.
- The shell-to-connector fold coupling (`h_fold` kernel per the doc) not implemented.

If the failing case has side curves, this is a credible independent contributor to non-convergence even after fixing the ghost flag.

### Why §1-7 missed it

My audit traced **code paths** (what runs when invoked) and checked internal consistency **within** each path. I did not compare **default member-variable values** against the pre-Option-B state. `mCreateGhostFacets` itself is referenced and correctly used inside the code; my checks of "is `create_ghost_facets` called correctly" didn't catch that the entire `if (mCreateGhostFacets)` block is now gated off at construction time.

**Lesson for future audits**: when looking for regressions, diff default initializers and class-level constant flips as a separate explicit step, not just call paths.

### Correction to §4's suspect ranking

§4's candidate list is largely moot. The real suspect ranking is:

1. **`mCreateGhostFacets = false` regression** (this addendum) — primary, high confidence.
2. **`mConnectorsForAllLayers = true` with incomplete fold closure** — secondary, if the case has side curves.
3. §4's other suspects — tertiary, only if re-enabling ghost facets doesn't restore convergence.

### Recommended action (still read-only)

1. Flip `mCreateGhostFacets` back to `true` in `cl_ThinShellFactory.hpp:66`.
2. Rebuild, re-run the failing hphi case.
3. If convergence returns: case closed for (1). The side-connector fold closure (2) remains a known incomplete feature per the existing docs — address separately.
4. If not: proceed with §5's instrumentation plan.

