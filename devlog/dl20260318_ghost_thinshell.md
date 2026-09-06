# Ghost Method for Thin-Shell Layer Interfaces

**Date:** 2026-03-18
**Purpose:** Implementation log for DG-decoupled thin-shell edges and ghost facet infrastructure
**Module:** `src/mesh/`
**Branch:** `ghost`
**Status:** In progress — mesh infrastructure complete, DOF adjacency and assembly next

---

## 1. Problem Recap

In the h-φ thin-shell formulation, adjacent material layers with extreme resistivity contrast (e.g., YBCO next to hastelloy) share DOFs at their layer boundary. This forces C⁰ continuity of H_t across a jump of 5+ orders of magnitude in transport properties, causing severe ill-conditioning.

The solution is to decouple layer interfaces via DG (separate edge DOFs per side) and couple them weakly through Nitsche terms. See `todo/thinshell_selective_nitsche_coupling.md` for the full design.

---

## 2. Changes Made

### 2.1 ThinShellFactory: Decoupled Edge/Face Creation

**Files:** `cl_ThinShellFactory.hpp`, `cl_ThinShellFactory.cpp`

- **Layer struct** now has three edge groups (`EdgesTop`, `EdgesMid`, `EdgesBottom`) and three face groups (`FacesTop`, `FacesMid`, `FacesBottom`) plus `GhostFacets`, replacing the single `Edges`/`Faces` arrays.

- **`create_edges_on_layers(aOrder, ...)`** assigns edge groups based on layer position and element order:
  - First-order: even layers are boundary layers. Layer[0] gets `EdgesBottom` only, Layer[n-1] gets `EdgesTop` only, interior layers get both.
  - Second-order: even layers (block boundaries) get `EdgesTop` + `EdgesBottom` for DG decoupling; odd layers (mid, interior to block) get `EdgesMid` only.
  - Result: adjacent blocks no longer share edge objects at their boundary.

- **`create_faces_on_layers(aOrder, ...)`** mirrors the edge logic for faces (second-order only).

- **`link_elements_with_edges()`** updated for both orders:
  - First-order: bottom from `Layer[l].EdgesBottom`, top from `Layer[l+1].EdgesTop`
  - Second-order: bottom from `EdgesBottom`, mid from `EdgesMid`, top from `EdgesTop`

- **`link_elements_with_faces()`** updated similarly, plus `tCount++` fix.

- **ID management:** `mMaxID` split into `mMaxGroupID` (blocks/sidesets) and `mMaxElementID` (elements/edges/faces), both as member variables. All creation functions use the member counters directly, avoiding stale local copies.

### 2.2 Ghost Facet Creation

**File:** `cl_ThinShellFactory.cpp`

- **`create_ghost_facets(aOrder, ...)`** creates ghost facets at each inter-block boundary:
  - Loops `b = 1..nb-1`, computes boundary layer via `l += aOrder`
  - For each facet position, creates a placeholder element and a ghost `Facet`
  - Sets master = block b-1's element (face index 1 = top), slave = block b's element (face index 0 = bottom)
  - Stores in `tLayer->GhostFacets`
  - Early return for single-block thin shells

- **`create()`** updated:
  - Calls `create_ghost_facets()` after block setup
  - Moves all edge/face groups to mesh (`append_move` with size checks)
  - Creates ghost `SideSet` with `DomainType::Ghost`, populates with ghost facets
  - Pushes ghost sideset to mesh
  - Passes both sidesets to `ThinShell` constructor

### 2.3 ThinShell: Ghost Sideset Support

**Files:** `cl_ThinShell.hpp`, `cl_ThinShell.cpp`

- Constructor takes `(SideSet*, SideSet* aGhostSideSet)`, hides both from mesh output
- `ghost_id()` returns ghost sideset ID or `gNoID` if null
- `ghost_facets()` returns ghost sideset's facets or empty `mNull` cell
- `mNull` member provides safe empty return for no-ghost case

### 2.4 Distributor / ProtoMesh: Ghost Serialization

**Files:** `cl_Mesh_Distributor.cpp`, `cl_ProtoMesh.hpp`, `cl_ProtoMesh.cpp`

- `ThinShellData::mGhostSideSetID` added (defaults to `gNoID`)
- `send_thinshell_data()` serializes `ghost_id()` into broadcast buffer
- `receive_thinshell_data()` deserializes `mGhostSideSetID`
- `ProtoMesh::create_thinshells()` looks up ghost sideset by ID with `key_exists` fallback to `nullptr`

### 2.5 Periodicity: Mode-Based Edge Matching

**File:** `cl_Mesh_Periodicity.cpp`

- `match_edges()` refactored into mode-based dispatch:
  - Mode 0: volume conductor edges (unchanged logic)
  - Modes 1..max_order: thin-shell edge groups, one mode per face plane
- Each mode builds its own `Map` independently → no key collisions between edge groups that share the same nodes at layer boundaries
- Uses `DynamicBitset` to select edges by local index range within each element
- DOF manager handles hanging+periodic overlap correctly: if master DOF is hanging, slave inherits the same sources (verified at `cl_FEM_DofMgr_DofData.cpp:3501`)

---

## 3. Mesh Lifecycle

Ghost facets flow through the standard pipeline:

1. `ThinShellFactory::create()` → creates ghost facets, pushes ghost sideset to `mMesh->sidesets()`
2. `Mesh::finalize()` → `collect_facets_from_sidesets()` picks up ghost facets into `mMesh->facets()`
3. `update_facet_nodes()` → ghost facet placeholder elements receive nodes from master element's top face
4. `Partitioner::fix_facet_related_ownerships()` → assigns `min(master->owner(), slave->owner())`
5. `Distributor` → serializes/deserializes ghost sideset ID alongside thin-shell data
6. `ProtoMesh::create_thinshells()` → reconstructs ThinShell with ghost sideset on worker procs

---

## 4. Known Bugs to Fix

### 4.1 Periodicity edge matching: wrong slice size (High) — **FIXED**
`Periodicity::match_edges()` now uses `number_of_edges( tShell->element_type() )` (per-plane count: TRI3→3) instead of `number_of_edges( tBlock->element_type() )` (total: PENTA6TS→6). Verified 2026-03-19.

### 4.2 Misleading comments in match_edges — **FIXED**
Comments now correctly read "mode 1: bottom plane, mode 2: middle plane (second order only)". Verified 2026-03-19.

---

## 5. What's Next

### 5.1 Edge-to-edge adjacency for ghost facets
`ConnectivityCalculator::connect_edges_to_edges()` (`cl_Mesh_ConnectivityCalculator.cpp:614`) discovers edge neighbors through shared elements. With decoupled layers, EdgesTop (block k) and EdgesBottom (block k+1) share no elements and are invisible to each other.

**Plan:** After the normal element-based edge-to-edge loop, iterate the ghost sideset's facets. For each ghost facet, the master element's edges (including EdgesTop) and the slave element's edges (including EdgesBottom) should be registered as neighbors of each other. This injects the cross-layer adjacency at the mesh level, which then propagates correctly to the DOF manager's sparsity pattern.

### 5.2 Ghost sideset registration in MaxwellFactory / IWG
The ghost sideset ID must be added to `IWG::selected_sidesets()` so that `SolverData::compute_element_dof_connectivity()` includes ghost facet element DOFs in the sparse matrix graph. Entry point: `MaxwellFactory` setup, where thin-shell sidesets are registered.

### 5.3 Ghost facet assembly (Nitsche terms)
Implement the nonsymmetric Nitsche bilinear form (consistency + penalty) using the ghost facets' master/slave element pairs. This fills the 12×12 interface matrix blocks allocated by the sparsity pattern. See `thinshell_selective_nitsche_coupling.md` Section 3 for the concrete pseudo-code.

### 5.4 Future
- Selective DG (only high-contrast interfaces) — Phase 3
- Symmetric Nitsche (add adjoint consistency) — if convergence order matters
- Full Newton linearization of Nitsche terms

---

## 6. Related Documents

- `todo/thinshell_selective_nitsche_coupling.md` — Full design document
- `src/fem/maxwell/doc/maxwell_usage_guide.md` — Maxwell module guide
