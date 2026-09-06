# Ghost Method for Thin-Shell Layer Interfaces (Day 2)

**Date:** 2026-03-19 (updated 2026-03-20)
**Purpose:** DOF adjacency, sparsity pattern, and Nitsche assembly for ghost facets
**Module:** `src/mesh/`, `src/fem/maxwell/`, `src/fem/kernel/`
**Branch:** `ghost`
**Status:** In progress — assembly implemented, ready for validation
**Scope:** Phase 1 — low-order (edge_multiplicity=1), nonsymmetric Nitsche
**Previous:** `devlog/dl20260318_ghost_thinshell.md`

---

## 1. Starting Point

Mesh infrastructure is complete (see yesterday's log). Bugs 4.1 and 4.2 from yesterday were verified as already fixed in the code — the dev log was stale. Today picks up from Section 5 (What's Next).

---

## 2. Key Findings

### 2.1 Edge-to-edge connectivity is NOT on the critical path
The DOF sparsity pattern is built entirely from element DOF lists (`SolverData::compute_element_dof_connectivity()` at `cl_FEM_DofMgr_SolverData.cpp:957–1064`), iterating volume blocks + selected sidesets. Edge-to-edge adjacency is never consulted. The mesh-level `connect_edges_to_ghost_facets()` and `connect_faces_to_ghost_facets()` helpers have been commented out — they may be deleted later.

### 2.2 Ghost DOF linking already supported
`IWG_Maxwell` uses `SideSetDofLinkMode::MasterAndSlave` (`cl_IWG_Maxwell.cpp:42`). `DomainType::Ghost` is NOT in the forced-Inactive list (`cl_FEM_Element.cpp:314`). So the generic `link_dofs_master_and_slave()` path will collect DOFs from both master and slave elements automatically.

### 2.3 Historical note: all-edges connection rationale (commented-out code)
The commented-out `connect_edges_to_ghost_facets()` connected ALL edges of master/slave to each ghost facet (not just interface-face edges). This was correct for the Nitsche coupling graph (consistency term uses all 6 DOFs per element). However, this code is off the critical path (Section 2.1) and may be deleted (see F1).

### 2.4 select_sidesets() overwrites — registration was broken
`IWG::select_sidesets()` (`cl_IWG.cpp:289`) writes `mSideSetIDs`, the same member as `IWG::set_sidesets()` (`cl_IWG.cpp:1715`). The call at `cl_MaxwellFactory.cpp:480` replaced the full Maxwell sideset list (interfaces, symmetry, thin-shell) with ghost-only IDs. Additionally, `tCount` was not reset before the fill loop (line 477), causing out-of-bounds writes. Fixed in A1.

### 2.5 Edge direction invariant — no bitset extension needed
Master and slave edge directions are **identical** for thin-shell ghost facets (see C1 for full rationale). The existing `Bitset<12>` in `mEdgeDirections` (`cl_FEM_Element.hpp:75`) is sufficient — no extension to 16 bits was needed. The single bitset loaded from master is correct for both sides.

### 2.6 DOF linking dispatch — confirmed correct
Thin-shell PENTA6TS blocks are typed as `DomainType::Conductor` (not `ThinShell`) by the input file domain assignment (`cl_MaxwellFactory.cpp:318`) before `select_blocks()` runs. They ARE selected and exist in the FEM DofManager. The ghost facet's master/slave elements are found with valid block IDs, so the dispatch stays in `link_dofs_master_and_slave()` — no fallback to `FacetAndMaster`. Confirmed independently by Codex (confidence: high).

### 2.7 Nitsche formulation verified against literature
The harmonic weighting formula was verified against Ern, Stephansen & Zunino (2009) and Zunino (2009). Key result: δ⁺ρ⁺ = δ⁻ρ⁻ = ρ_harm = ρ⁺ρ⁻/(ρ⁺+ρ⁻), explicitly stated in Ern et al. eq. (2.16). The nonsymmetric variant only requires η > 0 for coercivity (Zunino 2009, p. 122). Starting value η = 4.

### 2.8 Ghost facet orientation
Ghost facets set slave orientation to 1 (`set_slave(element, 0, 1)` at `cl_ThinShellFactory.cpp:1738`). This is always correct because extrusion preserves node winding across layers.

---

## 3. Todo List

### A. Fix bugs in current code

- [x] **A1.** Remove broken ghost sideset collection + `select_sidesets()` call (`cl_MaxwellFactory.cpp:462–480`). Ghost sidesets go through the typed sideset path instead.

### B. Register ghost sidesets in the typed Maxwell path

- [x] **B1.** Add `DomainType::Ghost` to switch in `set_block_types_in_magnetic_equation()` (`cl_MaxwellFactory.cpp:1712`). Ghost sidesets now included in `set_sidesets(tIDs, tTypes)`.
- [x] **B2.** Add `DomainType::Ghost` case to `FieldList::collect_sideset_dofs()` (`cl_Maxwell_FieldList.cpp:415`). Ghost DOF table populated via `Ghost.push(tDof)` in `initialize()` (line 125), which copies the full Conductor DOF list. **Phase 1 note:** for low-order, this is `edge_h` only. For higher-order (`mHigherOrder=true`), Conductor also includes `face_h` — the ghost DOF table will automatically include it, but `link_dofs_master_and_slave()` will hit the existing assert at `cl_FEM_Element.cpp:297` ("can't have face dofs on both master and slave"). Higher-order ghost support requires relaxing that assert (deferred). Also added `unique_and_rearrange(Ghost)`. Removed stale Ghost→NonDof routing.

### C. Edge direction for dual-edge-DOF sidesets — NO FIX NEEDED

- [x] **C1.** Verified: master and slave edge directions are **identical** for thin-shell ghost facets. All layer edges are cloned from the same originals with preserved node ordering (`cl_ThinShellFactory.cpp:1511`), and PENTA6TS/QUAD4TS element templates are constructed so top and bottom faces traverse the interface in the same direction. The existing single `mEdgeDirections` bitset in `link_dofs_master_and_slave()` is correct. Confirmed independently by Codex (confidence: high for linear, moderate for higher-order templates).
- [x] **C2.** Added `#ifdef DEBUG` check in `link_dofs_master_and_slave()` (`cl_FEM_Element.cpp:808–815`): when both master and slave carry edge DOFs, compares full `Bitset<12>` via new `operator==` (`cl_Bitset.hpp`). Also added `operator!=` for completeness.

### D. Assembly — DONE

- [x] **D1.** Add `DomainType::Ghost` to `IWG_Maxwell::link_to_group()` (`cl_IWG_Maxwell.cpp:554`). Routes to `maxwell::h_ghost()`.
- [x] **D2.** Implement `maxwell::h_ghost()` (`mt_maxwell_h.cpp:1794–1928`). Nonsymmetric Nitsche: penalty + consistency terms. Uses harmonic-weighted ρ, constant material evaluation (order-of-magnitude sufficient), 2D shape function Ê for through-thickness derivative. Assembles 12×12 K from four 6×6 blocks (Kmm, Kms, Ksm, Kss). Verified signs against Ern et al. (2009) by three independent reviewers (Claude, Codex, Grok).

**Bugs found and fixed during D2 review:**
- Alpha formula was inverted (divided instead of multiplied by 1/h⁺+1/h⁻)
- Ksm was not reset to zero before accumulation
- Kms consistency sign was flipped by the `−=` operator
- Ksm penalty sign was wrong (`+=` instead of `−=`)
- Kss was transposed in final 12×12 assembly (leftover from penalty-only draft)

### E. Validation

- [ ] **E1.** Compile and run. Confirm ghost FEM element gets 12 DOFs (6 master `edge_h` + 6 slave `edge_h`). Validates the full registration → linking → sparsity → assembly chain.

### F. Cleanup and documentation

- [ ] **F1.** Decide on commented mesh-connectivity helpers (`connect_edges_to_ghost_facets`, `connect_faces_to_ghost_facets`). If the FEM sideset path provides correct sparsity, delete them.
- [ ] **F2.** Update `thinshell_selective_nitsche_coupling.md` Section 3: still says "assembly via faces — no new sidesets needed", but the implementation uses ghost facets + a dedicated ghost sideset. Align with current architecture.

### G. Robustness guards (from Codex review, 2026-03-20)

- [ ] **G1.** ProtoMesh ghost reconstruction: `create_thinshells()` (`cl_ProtoMesh.cpp:785`) silently accepts a missing ghost sideset ID via `key_exists` fallback to nullptr. Add assertion or warning for distributed runs where the ghost sideset must be present.
- [ ] **G2.** Higher-order guard: no mechanism prevents `mHigherOrder=true` from activating ghost sidesets. Ghost DOFs inherit the full Conductor set (including `face_h`), but `link_dofs_master_and_slave()` asserts against face DOFs on both sides (`cl_FEM_Element.cpp:297`). Add an early guard or clear error message in MaxwellFactory when higher-order + ghost is attempted.
- [ ] **G3.** Hanging edge + ghost interaction: untested combination. DOF-source inheritance in `cl_FEM_DofMgr_DofData.cpp:3497` could interact with ghost interface DOFs. Needs a test case with both hanging thin-shell edges and ghost interfaces.

---

## 4. Priority

A, B, C, D done. **E1 is the next step** — compile and validate. F is cleanup. G is robustness for future.

---

## 5. Future: Nitsche at shell/volume interface (Phase 4 idea)

**Concept:** Replace hanging edges at the thin-shell/volume boundary with Nitsche coupling, giving a uniform DG+Nitsche architecture for all thin-shell interfaces (inter-layer AND shell/volume). This would eliminate the T-matrix cascade logic in `cl_FEM_DofMgr_DofData.cpp:3486+` — the most fragile code path in the DOF manager.

**Why it's appealing:**
- One coupling mechanism everywhere — no hanging DOFs, no cascades, no source inheritance chains
- Removes the `tEdge->dof()` vs `tOther->dof()` subtleties and the two-interface volume edge problem
- Architecturally cleaner: all thin-shell boundary conditions become weak enforcement

**Why it's harder than inter-layer Nitsche:**
- Inter-layer is 1D scalar DG (design doc Section 2, "key simplification"): ∂h_t/∂z is a finite difference of E evaluations, C never appears. Brenner & Scott Ch. 10 applies directly.
- Shell/volume interface is a 2D surface in 3D: the "normal derivative" involves the volume element's curl operator at the interface, not a simple through-thickness difference. This is closer to full H(curl) DG (Houston, Perugia & Schotzau 2004–2005).
- The shell/volume interface is NOT a high-contrast material jump — the conditioning motivation doesn't apply. Hanging edges give exact continuity, which is physically correct here.

**Prerequisite:** Validate Phase 1 inter-layer Nitsche first. If that works well, the shell/volume extension becomes a natural follow-up with established infrastructure.

**Literature needed:** Houston, Perugia & Schotzau (2004–2005) for H(curl) DG theory at mixed-formulation interfaces.

---

## 6. Related Documents

- `devlog/dl20260318_ghost_thinshell.md` — Day 1 log (mesh infrastructure)
- `todo/thinshell_selective_nitsche_coupling.md` — Full design document (needs update, see F2)
- `src/fem/maxwell/doc/maxwell_usage_guide.md` — Maxwell module guide
- Ern, Stephansen & Zunino (2009), IMA J Numer Anal 29:235–256 — SWIP method, harmonic weighting
- Zunino (2009), J Sci Comput 38:99–126 — WIP method, nonsymmetric variant coercivity
