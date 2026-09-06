# Corner-Edge Transition Problem (Thin Shell Conductor/Air)

**Date:** 2026-03-13  
**Author:** Codex (third pair)  
**Purpose:** Handoff note for continuing investigation tomorrow

---

## Problem Context

After fixing thin-shell edge key canonicalization in `ThinShellFactory::create_temporary_edges()` (using `original()->index()` + canonical ordering), many previously duplicated thin-shell edges are now correctly shared.

This exposed a potential corner-case at conductor/air transition boundaries on thin-shell surfaces:
- some adjacent thin-shell facets want **edge-to-edge** hanging coupling (conductor side),
- others want **edge-to-node** hanging coupling (air side),
- but they now reference the **same** thin-shell edge object.

Current hanging setup uses first-writer-wins (`is_flagged()` + skip), so only one relation may survive.

---

## What We Confirmed

1. **Face/edge indexing on PENTA6TS is consistent**
- bottom face = nodes `[0,1,2]`, top face = `[3,4,5]`
- bottom edges are separate from top edges (`0..2` vs `3..5`)
- `hang_*` functions use face indices correctly (`0` for bottom, `1` for top)

2. **`to_master_orientation` handling appears correct**
- separate overloads for node and edge containers
- TRI3 edge permutations differ from node permutations (as expected)
- Claude’s manual permutation check is consistent with code

3. **Thin-shell Nedelec orientation path is internally consistent**
- `EF_PENTA6TS` uses per-edge sign `mS[k]`
- signs come from FEM edge direction logic (`compute_edge_directions`)
- no immediate mismatch found between local TS topology and Nedelec sign handling

4. **Main architectural constraint**
- edge hanging sources are currently assigned once per edge object
- data path in DOF T-matrix logic branches by first source type (`NODE` vs `EDGE`)
- mixed source kinds on one edge are not naturally handled by current implementation

---

## Core Insight

The key fix is likely **correct** (deduplicating by original node identity), but it reveals an existing assumption:

> a hanging edge is expected to have one consistent coupling mode.

At transition edges, this assumption can fail.  
Previously, accidental duplicate edges could mask this by letting each side carry its own mode.

---

## Physical/Modeling Principle Raised by Christian

If an edge touches an air element in this H-phi setup, it should not carry an independent edge DOF in the final effective dependency; it should resolve to node (phi) DOFs.

This suggests for transition edges:
- node-based dependency should likely be prioritized over edge-based dependency.

---

## Open Technical Question

Does an edge-to-edge source on a transition edge eventually collapse to node DOFs via cascading (because source edge is itself hanging)?

- Sometimes possible.
- Not guaranteed in general.
- Therefore not safe to assume equivalence of edge-based and node-based paths.

---

## Candidate Direction (No code change yet)

For edges at conductor/air transition on thin-shell interfaces:

1. either ensure they resolve to **node-based** sources (preferred),
2. or split transition edges so edge-based and node-based relations do not compete on one edge object.

Dual mixed source kinds on the same edge object likely requires larger refactor.

---

## Suggested Diagnostics For Tomorrow

1. Identify transition edges (shared by facets requiring different coupling modes).
2. For those edges, print:
   - edge id
   - number/source type(s)
   - which `hang_*` path assigned them first
3. Check whether transition edges that should be node-resolved end up with only `EDGE` source type.
4. Correlate this with stagnation case vs working reference.

---

## Related Files

- `src/mesh/cl_ThinShellFactory.cpp` (`create_temporary_edges`)
- `src/fem/maxwell/cl_MaxwellFactory.cpp` (`create_hanging_edges_and_facets`, `hang_*`)
- `src/fem/kernel/cl_FEM_DofMgr_DofData.cpp` (`create_dofwise_t_matrices_master`)
- `src/fem/interpolation/nedelec/cl_EF_PENTA6TS.cpp`

