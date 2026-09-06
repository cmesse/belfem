# Thin-Shell Overhang Bug: Node Duplication Failure Analysis

> **CANCELLED 2026-09-03** (todo/ currentness sweep, round 3): low relevance — the literal first-facet bug was fixed months ago; the §9 consistent-block relink has been "still unimplemented" since February 2026 with no user case pulling it. Status lines and checkboxes below are as they stood at closure and are not maintained.

**Date:** 2026-02-20
**Refreshed:** 2026-06-22 — re-audited against the current tree (Claude + Codex + Grok)
**Purpose:** Analyze why cohomology cuts fail when thin shells extend beyond a solid conductor block
**Module:** homology (CutFactory), fem/maxwell (MaxwellFactory)

---

## Currentness note (2026-08-09, anchors re-baselined 2026-08-11)

Re-checked against the tree: the verdict below is unchanged — the literal first-facet bug is
fixed, the upstream normalization exists, and the robust §9 consistent-block relink is **still
unimplemented**, so this file stays ACTIVE. Only the anchors have moved:
`CutFactory::relink_slave_elements_with_duplicate_nodes()` is `cl_CutFactory.cpp:1730`
(quoted below as 1871-1949) — unchanged since the 2026-08-09 sweep — while
`MaxwellFactory::fix_facet_masters()` has moved again and is now
`cl_MaxwellFactory.cpp:1243` (declared `cl_MaxwellFactory.hpp:198`, called from `:937`;
quoted below as :760 / body :1027-1085). It drifted by +129 lines in two days, which is the
standing argument for re-locating by symbol rather than trusting any number in this file.

Cross-reference worth carrying: `fix_facet_masters` has since acquired a second known
weakness on all-air 2-D meshes — its propagation queue is seeded only by domain-type-contrast
facets, so the orientation sweep never runs there (`deferred/2d_thinshell_todo.md` B9,
`debt_register.md` DR-16). That is a different failure of the same routine this document
relies on for its "upstream normalization now exists" argument.

## Status Update (2026-06-22, 3-AI re-audit)

**The literal bug this doc describes — sampling the slave block from the *first facet only* —
is fixed.** `relink_slave_elements_with_duplicate_nodes()` (now `cl_CutFactory.cpp:1871-1949`)
collects the slave block of **every** facet (`1902-1905`) and relinks all of them, instead of
`tSlaveID = facets()(0)->slave()->block_id()`.

**Upstream normalization now exists.** `MaxwellFactory::fix_facet_masters()`
(`cl_MaxwellFactory.cpp:760`, body `1027-1085`) runs **before** the cut step and normalizes
master/slave: cross-type facets (Conductor/Air) flip so the higher `DomainType` is master
(→ slave = Air); same-type (Air/Air overhang) facets are oriented by a BFS from those
cross-type seeds via `check_facet_orientation`. For the geometry in this doc (Conductor/Air
facets share the sideset with the overhang), that BFS makes the slave the consistent phi
side, so collect-all-slaves yields a consistent split.

**Why this stays ACTIVE, not closed:**
- The §9 **consistent-block** policy (the block present on *every* facet) is **not**
  implemented in `CutFactory`. Collect-all-slaves is correct *only* when `fix_facet_masters`
  has already made the slave the consistent side; `CutFactory` neither checks nor enforces
  this invariant.
- `fix_facet_masters` gives orientation **consistency**, not a direct common-block selection
  (Codex). A thin-shell sideset with **no cross-type seed** orients arbitrarily, and the
  relink side may then not match the intended duplicate side (both auditors flagged this).
- `fix_facet_masters` runs only on the cohomology path (`mComputeCohomologies`); other
  callers of `create_thin_shell_cuts()` would bypass it (today only `create_cuts_sub_master`
  calls it).

**Auditor split:** Grok rated the documented single-sideset overhang **FIXED** in the
cohomology pipeline; Codex rated it **PARTIALLY-FIXED** (invariant delegated, not enforced).
Neither ran an overhang regression mesh — both verdicts are code-path analysis.

**Related (separate) bug:** the multi-sideset *protoshell common-block* inconsistency
(`tapestack3d`) is tracked in `src/homology/doc/homology_usage_guide.md` (formerly the standalone `doc/cohomology_block_id_bug.md`); its proposed "common block →
slave" second pass is also unimplemented. The §9 fix here and that second pass are the same
idea at two scopes (per-sideset vs per-protoshell).

**Recommended next step:** implement the §9 consistent-block relink (or extend
`fix_facet_masters` with the common-block second pass) and add a debug assert that
`tBlockBitset` never holds blocks from both physical sides of a thin-shell sideset. The §2–§10
analysis below describes the *original* (pre-fix) code and remains the rationale for that work.

---

## 1. Problem Statement

Consider a 3D geometry whose 2D side view is:

```
(11) ----------------------------- (10)
|                                    |
|    (1) --- (2) ---- (3) ---- (7)   |
|    |       |        |        |     |
|    |  air  | condr. |  air   |     |
|    (0) ----(4) ---- (5) ---- (6)   |
|                air                 |
(8)--------------------------------(9)
```

- **Top thin shell:** (1)--(2)--(3)--(7)
- **Bottom thin shell:** (0)--(4)--(5)--(6)
- **Solid conductor:** box (2)--(3)--(5)--(4)
- **Air blocks:**
  - Air_above: region above the top thin shell
  - Air_left: region between the shells, left of conductor
  - Air_right: region between the shells, right of conductor
  - Air_below: region below the bottom thin shell

The thin shells are defined as internal sidesets (each facet has a master and a slave
element from **different** blocks).  The thin shells divide the air domain into separate
blocks; GMSH creates distinct volumes on each side of each internal surface.

The segments (1)--(2) and (5)--(6) are **overhangs**: portions of the thin shell that
extend beyond the conductor into regions where both sides are air (but from different
blocks).

Cohomology cuts work correctly when the shells do not overhang (i.e., if nodes (1), (0),
(7), (6) did not exist and the shells exactly coincided with the conductor surface).

---

## 2. Background

### Why Cuts Exist

The mixed h-phi formulation solves:

- **H-formulation** in conductors (edge DOFs for the magnetic field H)
- **phi-formulation** in air/ferro (nodal DOFs for the scalar potential phi)

| Cut Type | Purpose | Mechanism |
|----------|---------|-----------|
| **Thin-shell cuts** | Embed 2D shell elements into the mesh by splitting nodes along the shell surface | Node duplication along the shell |
| **Cohomology cuts** | Enforce Ampere's law `[phi] = I` across a surface in the phi domain | Node duplication along a topologically-identified surface |

In `MaxwellFactory::create_cuts_sub_master()` (`cl_MaxwellFactory.cpp:719`), thin-shell
cuts are created **first** via `CutFactory::create_thin_shell_cuts()`, then cohomology
cuts via `CutFactory::run()`.  The thin-shell node duplication alters the mesh topology
that the cohomology step subsequently operates on.

### How Master/Slave Is Assigned

A critical implementation detail: BELFEM assigns facet master/slave based on **element
ID order** (`cl_Mesh_ConnectivityCalculator.cpp:222-248`).  The element with the
**lower ID** becomes the master.  This is deterministic per mesh, but **not related to
block type, domain type, or geometric normal direction**.

GMSH does not produce consistent surface orientations.  The code uses the master block
(from the first facet of a sideset) to establish a reference orientation for the thin
shell.  This works when all facets in a sideset sit between the same two blocks, but
fails for the overhanging configuration.

---

## 3. Cut Creation Workflow

```
MaxwellFactory::create_cuts_sub_master()          [cl_MaxwellFactory.cpp:719]
├─ CutFactory::create_thin_shell_cuts()            [cl_CutFactory.cpp:1603]
│   ├─ create_curves_for_thinshells()              [Identify terminal/side curves]
│   ├─ duplicate_nodes_on_face_sidesets()          [Create duplicate Node objects]
│   ├─ relink_slave_elements_with_duplicate_nodes() [**BUG IS HERE**]
│   ├─ duplicate_and_relink_facets()               [Create 2 facets per original]
│   ├─ close_terminal_loops()                      [Duplicate terminal segments]
│   └─ append_move(mMesh->nodes(), duplicates)     [Add duplicates to mesh]
│
├─ create_terminal_list()
│
└─ CutFactory::run()                               [cl_CutFactory.cpp:110]
    ├─ unfinalize / finalize mesh                  [Rebuild topology]
    ├─ create_edges / create_faces                 [Recreate topological entities]
    ├─ compute_cohomologies()                      [Build simplicial complex from phi blocks]
    │   ├─ Flag nodes/edges/faces in phi blocks (Air, Ferro)
    │   ├─ Build SimplicialComplex
    │   ├─ Reduce/coreduce complex
    │   └─ Compute Smith normal form → cohomology generators
    ├─ compute_thin_cuts_and_duplicate_interface_nodes()
    └─ restore_thin_shell_sidesets()
```

---

## 4. Root Cause: Block-ID-Based Slave Identification

### The Faulty Code

In `relink_slave_elements_with_duplicate_nodes()` (`cl_CutFactory.cpp:1837-1892`):

```cpp
void CutFactory::relink_slave_elements_with_duplicate_nodes()
{
    // ...
    for( id_t tID : mTopology->groups( DomainType::ThinShell ) )
    {
        // ...
        // Lines 1849-1850: Determine master/slave from the FIRST facet
        id_t tMasterID = mMesh->sideset( tID )->facets()(0)->master()->block_id() ;
        id_t tSlaveID  = mMesh->sideset( tID )->facets()(0)->slave()->block_id() ;

        BELFEM_ERROR( tMasterID != tSlaveID, ... );  // Line 1852

        // Lines 1859-1868: Flag ALL elements matching tSlaveID
        for ( Node * tNode : tNodes )
        {
            for ( uint e=0; e<tNode->number_of_elements(); ++e )
            {
                Element * tElement = tNode->element( e ) ;
                if ( tElement->block_id() == tSlaveID )   // ← THE BUG
                {
                    tBitset.set( tElement->index() );
                }
            }
        }
        // ... relink flagged elements with duplicates ...
    }
}
```

### The Assumption and Why It Fails

The code assumes **all facets in a thin-shell sideset sit between the same two blocks**.
It samples `tSlaveID` from the first facet and uses it to identify slave elements for the
entire sideset.

For the top thin shell (1)--(2)--(3)--(7), the block pairs are:

| Region | Side A | Side B |
|--------|--------|--------|
| (2)--(3) (over conductor) | Conductor | Air_above |
| (1)--(2) (left overhang) | Air_left | Air_above |
| (3)--(7) (right overhang) | Air_right | Air_above |

Air_above is consistently on one side, but the **other side varies** (Conductor,
Air_left, or Air_right).

Since master/slave is determined by element ID (lower ID = master), the first facet's
`tSlaveID` could be **any** of these blocks.  Two failure modes arise:

---

## 5. Failure Modes

### Mode A: tSlaveID Picks the Wrong Block

If the first facet happens to have the lower-ID element on the Air_above side, then:

```
tMasterID = Air_above
tSlaveID  = Conductor   (or Air_left, or Air_right)
```

The code relinks all elements with `block_id == tSlaveID`.  For example, if
`tSlaveID == Conductor`:

- **Conductor region (2)--(3):** Conductor elements get relinked.  Air_above keeps
  originals.  This is the **reverse** of the intended split.
- **Overhang (1)--(2):** No elements have `block_id == Conductor` here, so **nothing
  is relinked**.  Both Air_left and Air_above keep original nodes.  The thin shell
  fails to create a topological gap in the overhang.

**Result:** The overhang region has no node split at all.  Air_left and Air_above share
the same nodes across the thin shell.

### Mode B: tSlaveID Happens to Be Correct

If the first facet happens to have the lower-ID element on the non-Air_above side:

```
tMasterID = Conductor (or Air_left, or Air_right)
tSlaveID  = Air_above
```

All elements with `block_id == Air_above` get relinked.  Since Air_above is consistently
on one side of the thin shell across all regions, this is correct.

**However, this only works by luck of element ID ordering.**  There is no guarantee
that the first facet's element IDs are ordered this way.

### Why the Non-Overhanging Case Always Works

When the thin shell does **not** overhang (all facets between Conductor and Air_above),
`tSlaveID` is always one of these two blocks, regardless of element ID ordering.
In Mode A (tSlaveID = Conductor), the conductor elements are relinked and air keeps
originals.  In Mode B (tSlaveID = Air_above), air elements are relinked and conductor
keeps originals.  Either way, the split is **consistent across all facets** because the
block pair is uniform.  The FEM coupling handles the reversed normal correctly.

In the overhanging case, the block pair **changes** along the sideset, so the
block-ID-based approach cannot produce a consistent split.

---

## 6. Downstream Impact on Cohomology

After `create_thin_shell_cuts()` returns, `CutFactory::run()` rebuilds the mesh topology
and computes cohomologies.  The simplicial complex is built from **phi blocks only**
(Air and Ferro, per `Topology::select_blocks()` at `cl_Topology.cpp:356-398`).

In the failure case (Mode A), the overhang has no node split.  This means:

1. **Air_left and Air_above share original nodes** at the overhang.  The simplicial
   complex sees them as one connected region through the thin shell.

2. **The thin shell is topologically invisible** in the overhang region.  The
   simplicial complex does not reflect the physical barrier that the shell creates.

3. **The first homology group H_1** of the phi domain is incorrect.  The cohomology
   algorithm either misses the loop encircling the conductor+overhang assembly, or
   produces a generator that does not correctly account for the thin shell's extent.

4. **Cohomology cuts are not created correctly**, leading to failure to enforce the
   jump condition `[phi] = I` for transport current.

In Mode B (where tSlaveID happens to be Air_above), but the conductor region has its
normal reversed, there may be additional issues with how `duplicate_and_relink_facets()`
creates the paired facets for the FEM thin-shell coupling.

---

## 7. Junction Nodes (2), (3), (4), (5)

The conductor-shell interface nodes deserve special attention.  Node (2) is connected to:

- Conductor elements (non-phi block)
- Air_above elements (phi block, one side of the thin shell)
- Air_left elements (phi block, other side of the thin shell at the overhang)

In Mode A (tSlaveID = Conductor):
- Conductor elements at node (2) are relinked with duplicates
- Air_above and Air_left both keep originals
- Result: **no split** between Air_above and Air_left at node (2)

In Mode B (tSlaveID = Air_above):
- Air_above elements at node (2) are relinked with duplicates
- Conductor and Air_left keep originals
- Result: correct split (Air_above uses duplicate, others use original)

The junction nodes are correctly handled only in Mode B.

---

## 8. Mesh Contracts Consideration

Per the mesh contracts documentation (Section 5), node duplication invalidates
Node-to-Element and Element-to-Node connectivity.  The relinking at line 1881
(`tElement->insert_node()`) modifies element connectivity without rebuilding the
node's element list cache.  When multiple thin-shell sidesets are processed sequentially,
subsequent iterations use stale node-to-element connectivity.

The current code resets the bitset between iterations (line 1887) and re-flags nodes
(line 1889).  The stale connectivity causes redundant lookups but is harmless because:
- Already-relinked elements have their nodes replaced with duplicates
- Duplicates have `index == gNoIndex` → condition at line 1879 fails → no double-relinking

However, this stale-connectivity issue could interact with the block-ID bug if a node
appears on multiple thin-shell sidesets with different block pairs.

---

## 9. Proposed Fix Strategy

### Core Idea

Identify the **consistent side** of the thin shell — the block that appears on **every**
facet (as either master or slave).  Use that block's elements as the relink targets.

### Algorithm

```
For each thin-shell sideset:
    1. Iterate over all facets to find the "consistent block":
       - Collect block IDs from both master and slave of each facet
       - The block that appears on EVERY facet is the consistent side
       - (For non-overhanging sidesets, both blocks appear on every facet;
         pick one deterministically, e.g., the one with higher block_id,
         or use domain-type preference: prefer phi blocks.)
    2. Set tSlaveID = consistent_block_id
    3. Flag all nodes on the sideset (unchanged)
    4. Flag elements matching tSlaveID for relinking (unchanged)
    5. Replace flagged+duplicated nodes with their duplicates (unchanged)
    6. Reset bitset for next sideset (unchanged)
```

### Why This Works

- **Non-overhanging case:** Both blocks (Conductor and Air_above) appear on all facets.
  The tie-breaking rule picks one deterministically.  Same behavior as before.

- **Overhanging case:** Air_above appears on every facet.  Conductor, Air_left, and
  Air_right each appear only on their respective region's facets.  The algorithm
  identifies Air_above as the consistent side and relinks it.

- **Junction nodes:** Air_above elements at junction nodes (2), (3), (4), (5) are
  correctly identified and relinked.  Conductor and Air_left/Air_right keep originals.

### Tie-Breaking for Non-Overhanging Case

When both blocks appear on all facets, several tie-breaking strategies are possible:

1. **Domain-type preference:** Prefer phi blocks (Air, Ferro) over non-phi blocks
   (Conductor, Coil).  This matches the physical intent: the phi domain should have
   the split.

2. **Higher block_id:** Deterministic but arbitrary.

3. **First facet's slave:** Falls back to current behavior.  Since the non-overhanging
   case works with either choice, this is safe.

Option 1 is recommended because it expresses the physical intent clearly and handles
edge cases where the conductor has elements on both sides (unusual but possible).

### The Assertion

The assertion `BELFEM_ERROR(tMasterID != tSlaveID, ...)` at line 1852 remains valid:
every facet of a thin-shell sideset has master and slave from different blocks (the
thin shell separates two distinct mesh regions).  It does not need to be changed.

### Edge Cases to Verify

1. **Multiple sidesets per protoshell:** If a thin shell is defined by multiple
   sidesets, each is processed independently.  The consistent-block detection should
   handle this correctly per-sideset.

2. **Terminal curves at the conductor boundary:** `close_terminal_loops()` uses
   `tOrg->node(k)->index()` to decide whether a node has a duplicate.  With the fix,
   nodes at the conductor-shell interface (2), (3), (4), (5) are correctly duplicated
   (they are interior to the thin shell surface).  Terminal loops at these curves
   should close correctly because the side curve endpoints (at the mesh boundary in
   the z-direction) are NOT duplicated.

3. **Facet duplication (`duplicate_and_relink_facets`):** Creates facets tA
   (master side) and tB (slave side) with `set_master()` that links nodes from the
   volume element.  With correct relinking, tA references original or non-slave-side
   nodes, and tB references duplicate nodes from the consistent side.

4. **Cohomology validation:** After the fix, the simplicial complex should see the
   overhang region as two separate sides of the thin shell connected only at the
   non-duplicated side curve nodes.  The first Betti number should reflect the
   correct conductor-encircling loop topology.

---

## 10. Files Involved

| File | Lines | Role |
|------|-------|------|
| `src/homology/cl_CutFactory.cpp` | 1837-1892 | **Bug location**: `relink_slave_elements_with_duplicate_nodes()` |
| `src/homology/cl_CutFactory.cpp` | 1603-1643 | `create_thin_shell_cuts()` orchestration |
| `src/homology/cl_CutFactory.cpp` | 1646-1835 | `duplicate_nodes_on_face_sidesets()` (node creation) |
| `src/homology/cl_CutFactory.cpp` | 1896-1943 | `duplicate_and_relink_facets()` (facet splitting) |
| `src/homology/cl_CutFactory.cpp` | 1947-2084 | `close_terminal_loops()` (terminal segments) |
| `src/homology/cl_CutFactory.cpp` | 380-535 | `compute_cohomologies()` (simplicial complex) |
| `src/homology/cl_CutFactory.cpp` | 847-874 | `compute_thin_cuts_and_duplicate_interface_nodes()` |
| `src/homology/cl_Topology.cpp` | 356-398 | `select_blocks()` (phi vs non-phi classification) |
| `src/mesh/cl_Mesh_ConnectivityCalculator.cpp` | 222-248 | Master/slave assignment (element-ID-based) |
| `src/fem/maxwell/cl_MaxwellFactory.cpp` | 719-753 | `create_cuts_sub_master()` (entry point) |

---

## 11. References

- Messe et al. 2023 (paper1), Section 2.7: Nonlinear solver strategy and convergence
- Alves et al. 2022 (paper6): Thin-shell theory and cohomology cuts
- Arsenault et al. 2023 (paper3), Section II: Magnetodynamic h-phi coupling
- Mesh contracts documentation: `src/mesh/doc/mesh_contracts_and_invariants.md`, Section 5
- Homology theory documentation: `src/homology/doc/homology_cohomology_theory.md`
