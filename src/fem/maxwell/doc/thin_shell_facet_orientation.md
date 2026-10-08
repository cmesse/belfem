# Thin Shell Facet Orientation {#fem_maxwell_thin_shell_facet_orientation}

**Date:** 2026-03-02
**Purpose:** Documents why facet master/slave orientation matters for thin shells, the three-stage algorithm that ensures consistency, and the pitfalls that arise from mesh operations
**Module:** fem/maxwell, homology, mesh

---

## Key Files

| File | Key Function | Purpose |
|------|-------------|---------|
| `cl_MaxwellFactory.cpp` | `fix_facet_masters()` | Corrects master/slave assignment using DomainType; in 3-D, propagates orientation across same-type facets |
| `cl_MaxwellFactory.cpp` | `align_thin_shell_facets()` | 3-D: one winding per connected patch of a thin-shell sideset, one master block per block pair across its patches |
| `cl_MaxwellFactory.cpp` | `check_thin_shell_block_sides()` | 3-D: rejects a thin-shell sideset if a facet has only one side or if one block appears on both sides |
| `cl_CutFactory.cpp` | `relink_slave_elements_with_duplicate_nodes()` (~line 1407) | Replaces slave element nodes with thin shell duplicates |
| `cl_CutFactory.cpp` | `duplicate_nodes_on_face_sidesets()` (~line 1182) | Creates the duplicate nodes |
| `cl_Mesh.cpp` | `update_facet_nodes()` (~line 628) | Copies master element's face nodes into the facet surface element |
| `cl_Mesh_ConnectivityCalculator.cpp` | `connect_facets_to_elements()` (~line 171) | Initial master/slave assignment by element ID |
| `cl_Facet.cpp` | `flip()` (~line 89) | Swaps master and slave, relinks facet nodes |
| `fn_check_facet_orientation.hpp` | `check_facet_orientation()` (~line 39) | Checks if two adjacent facets have consistent outward normals |
| `en_DomainType.hpp` | `DomainType` enum (~line 19) | Defines priority ordering: Air(1) < Buffer(2) < Ferro(3) < Coil(4) < Conductor(5) |

---

## Problem Statement

Thin shells in BELFEM are surface-embedded conductors (superconducting tapes, thin films) modeled as sidesets rather than volume blocks. They require **topological separation**. The mesh must have duplicate nodes on one side of the shell so that the DOFs on both sides are independent. Which side gets the duplicates must be **consistent across all facets** of a given thin shell sideset.

Three mesh operations make this work:

1. **`update_facet_nodes()`** copies nodes from the facet's **master** volume element into the facet surface element. Therefore, the facet master must be the element that keeps the **original** nodes.

   > **It skips thin-shell facets.** `Mesh::update_facet_nodes()` returns early for sidesets whose
   > `domain_type()` is `ThinShell` or `GeometryOnly` (`cl_Mesh.cpp:637-641`), because those are
   > extrusion geometry whose node lists are built by the `ThinShellFactory` rather than copied
   > from a volume element. The orientation rule below still governs the *interface* facets between
   > the shell and the bulk, which is where master/slave actually matters.

2. **`relink_slave_elements_with_duplicate_nodes()`** replaces interface nodes in **slave** block elements with thin shell duplicates. Therefore, the slave side must receive the duplicates.

3. **`ThinShellFactory`** reads facet nodes via `tNode->facet()` to compute surface normals and build layered thin shell elements. It expects the facet surface element to reference original nodes with valid facet connections.

If master/slave assignment is inconsistent, duplicates go to the wrong side, and ThinShellFactory can compute incorrect or NaN normals. (`update_facet_nodes()` is not the failure path for the thin-shell sidesets themselves — it skips them, per the note above — but it is for the interface facets between shell and bulk.)

---

## Three-Stage Algorithm

### Stage 1: `MaxwellFactory::fix_facet_masters()` (cl_MaxwellFactory.cpp:1473-1560)

#### Initial state

For newly connected facets, `connect_facets_to_elements()` assigns master/slave by element ID: the element with the **lower ID** becomes the master. Facets that already have a master, such as those arriving from the distributor or from a `.bfm` file, keep their inherited master/slave links. This is arbitrary and depends on block ordering in the mesh file.

#### DomainType-based correction

`fix_facet_masters()` overrides this assignment based on `DomainType` priority. The DomainType enum values define this ordering (en_DomainType.hpp:19-73):

```
Air(1) < Buffer(2) < Ferro(3) < Coil(4) < Conductor(5) < Cut(8) < ThinShell(9)
```

The **higher** DomainType value becomes master. For a typical Air/Conductor interface, Conductor (5) becomes master and Air (1) becomes slave.

#### Cross-type facets (deterministic seeds)

For each facet where master and slave have **different** DomainTypes:

```cpp
uint tM = static_cast<uint>( block(facet->master()->block_id())->domain_type() );
uint tS = static_cast<uint>( block(facet->slave()->block_id())->domain_type() );

if ( tM < tS )      // master has lower priority → flip
    tFacet->flip();
else if ( tM > tS )  // master has higher priority → keep
    ;                 // (already correct)
```

These cross-type facets enter the BFS queue as correctly oriented **seeds**.

#### Same-type facets (BFS propagation)

Facets whose two sides have the **same** DomainType (for example, Conductor/Conductor) cannot be resolved by DomainType comparison alone. In 3-D, they are flagged for BFS processing.

In 3-D, BFS propagation from the cross-type seeds resolves these facets. **In 2-D, propagation is
disabled** (`tPropagate` in `fix_facet_masters()`). The thin-shell cut pipeline and the
thin-shell normals both read the tape facets' node order, and they need one master side per tape.
In 2-D that uniformity is an **input assumption, not something the code establishes**: the master is
the lower element id, so a tape's master side is uniform when each neighboring block's element ids form
a contiguous range, which Gmsh output gives and nothing here checks. Interleaved ids would leave a
mixed master along the tape; the cut pipeline's terminal-curve lookup then aborts with a named
`BELFEM_ERROR` rather than orienting anything.

The master rule is the same in both dimensions, except for propagation:

1. **Cross-type facets** (different DomainType on the two sides): the higher DomainType becomes the
   master, in 2-D and 3-D alike.
2. **Same-type facets**: in 3-D, BFS from the cross-type seeds; in 2-D, the master that
   `connect_facets_to_elements` assigned, which is the element with the **lower id**. The lower-id
   rule applies only here, never to a cross-type facet.
3. **Thin-shell sidesets** in 3-D are then made consistent by Stage 2 below.
4. **Deck-signed sidesets** are flipped after these steps (`flip_thin_shell_sidesets()`), whole
   sidesets at a time.

The 3-D propagation:

```cpp
while ( ! tQueue.empty() )
{
    Facet * tFacet = tQueue.pop();

    for ( uint f = 0; f < tFacet->number_of_facets(); ++f )
    {
        Facet * tOther = tFacet->facet( f );
        if ( ! tOther->is_flagged() ) continue;

        tQueue.push( tOther );
        tOther->unflag();

        if ( mesh::check_facet_orientation( tFacet, tOther ) ) continue;
        tOther->flip();   // flip if normals are inverted relative to seed
    }
}
```

`check_facet_orientation()` compares the traversal direction of the shared edge between two facets. If both facets traverse the shared edge in the **same** direction, their outward normals oppose each other, so one facet must be flipped.

#### Persistence across finalize cycles

The master/slave assignment set by `fix_facet_masters()` is preserved across `unfinalize()`/`finalize()` cycles because `FacetToElement` connectivity is not reset during those operations.

---

### Stage 2: `MaxwellFactory::align_thin_shell_facets()` and `check_thin_shell_block_sides()` (cl_MaxwellFactory.cpp:1565-1926)

The 3-D walk can leave a thin-shell sideset inconsistent. It orients a same-type facet from whichever
neighbor reaches it first, and the neighbors include faces that meet the shell at a T-junction: a
block face split at the tape, or the top face of a coil. On an edge shared by three facets, the
shell facet is edge-consistent with only one of the two other faces. A shell between two conductor
blocks can therefore be reached from a face above it at one end and from a face below it at the
other end, and come out with two windings. Both functions run at the end of `fix_facet_masters()`,
in 3-D only, while the facet-to-facet connectivity still exists.

`align_thin_shell_facets()` processes each thin-shell sideset of each protoshell in three steps:

1. **Patches.** A breadth-first pass over edge neighbors in the same sideset splits it into
   connected patches and records, for each facet, whether its winding is opposite to the first
   facet of its patch. A parity conflict means that the patch is not orientable, so the function rejects it with a
   `BELFEM_ERROR`.
2. **Vote per patch.** When the patch has cross-type facets, they decide, because their
   orientation follows the domain-type rule; otherwise all facets decide. The minority is
   flipped, and a tie keeps the winding of the first facet. Where the cross-type facets of one
   patch disagree, the log says so: the vote kept the majority winding, so some of them may now
   have the lower domain type as master.
3. **Relation between patches.** Disconnected patches of one sideset have no common edge, so the
   only well-defined relation between them is through their blocks. Patches with the same two
   blocks on their sides are given the same master block. The reference is the first patch whose
   vote had cross-type facets, or else the patch with the smallest facet id. The other patches
   with that pair are flipped whole where needed. A patch whose block pair no other patch shares,
   or whose facets do not share one pair, keeps its voted winding. Patches are grouped by their
   unordered block pair in one pass, and each group's reference is picked once. When the sideset
   has more than one patch, the log reports the number of patches, how many were aligned to a
   reference and how many of those were flipped, how many kept their winding because no other
   patch has their block pair, and how many kept it because their facets do not share one pair.

`check_thin_shell_block_sides()` then gathers, for each thin-shell sideset, the set of master
blocks and the set of slave blocks over its facets. It refuses the sideset with a `BELFEM_ERROR` when
a facet lacks an element on one side, or when one block is in both sets. Stage 3 relinks each
sideset's slave blocks as a whole. A block that is master on some facets and slave on others would
therefore be opened on both sides. This also refuses a uniformly wound sheet with the same block
on both of its sides. The message names the block and states that a thin shell must separate two
different volume blocks: the remedy is to split the volume block, not the sideset, so that each
facet has one block on each side, and to split the sideset as well only where it genuinely lies
between two different block pairs. It does not catch a split vote across two distinct block pairs, such as
masters {H1, L2} with slaves {L1, H2}, where H1 and H2 are blocks of the higher domain type and
L1 and L2 of the lower; that case is accepted, and the split-vote message is its only
diagnostic.

A processed mesh (`.bfm`) reloaded from an earlier run skips this stage and keeps the windings it was
saved with.

---

### Stage 3: `CutFactory::relink_slave_elements_with_duplicate_nodes()` (cl_CutFactory.cpp:1407-1483)

This function runs after `duplicate_nodes_on_face_sidesets()` has created the duplicate nodes. It relinks slave-side elements so they use duplicates instead of originals.

#### Per-facet slave block collection

The function collects slave block IDs from **all** facets of each thin shell sideset, not just the first one:

```cpp
for ( Facet * tFacet : tFacets )
{
    tBlockBitset.set( mMesh->block( tFacet->slave()->block_id() )->index() );
}
```

This matters for overhanging thin shells, whose facets may lie between different block pairs.

#### Element flagging

The function first flags all elements in the identified slave blocks. It then collects only flagged elements that touch sideset nodes:

```cpp
for ( Node * tNode : tSideSet->nodes() )
{
    for ( uint e = 0; e < tNode->number_of_elements(); ++e )
    {
        Element * tElement = tNode->element( e );
        if ( tElement->is_flagged() )
            tElementBitset.set( tElement->index() );
    }
}
```

#### Node replacement

For each collected element, interface nodes belonging to the thin-shell sideset are replaced with duplicates:

```cpp
tSideSet->flag_all_nodes();
for ( index_t e : tIndices )
{
    Element * tElement = tElements( e );
    for ( uint k = 0; k < tElement->number_of_nodes(); ++k )
    {
        Node * tNode = tElement->node( k );
        if ( tNode->is_flagged() && tNode->index() != gNoIndex )
        {
            tElement->insert_node( mThinShellDuplicates( tNode->index() ), k );
        }
    }
}
tSideSet->unflag_all_nodes();
```

Flagging nodes per sideset prevents nodes from different thin-shell sidesets from being replaced in the wrong context.

---

## Why This Works

After the three stages complete, the following invariants hold:

1. **`fix_facet_masters()`** ensures that for every thin shell facet, `master()` returns the volume element that will keep original nodes, and `slave()` returns the volume element that will receive duplicate nodes.

2. **`align_thin_shell_facets()`** and **`check_thin_shell_block_sides()`** ensure, in 3-D, that every thin-shell sideset has one winding per patch, one master block per block pair, and disjoint master-block and slave-block sets.

3. **`relink_slave_elements_with_duplicate_nodes()`** replaces interface nodes in slave-side elements with thin shell duplicates. Master-side elements keep originals.

4. **`update_facet_nodes()`** (called during `Mesh::finalize()`) copies nodes from the master element into the facet surface element of every facet outside the `ThinShell` and `GeometryOnly` sidesets. Since the master kept originals, those facet surface elements get original nodes. The thin-shell sidesets are skipped; their node lists are restored by `restore_thin_shell_sidesets()`.

5. **`ThinShellFactory`** reads facet nodes and computes surface normals. Because facet surface elements reference original nodes (restored for thin-shell sidesets, copied from the master elsewhere), the node-to-facet connections are valid, and normal computation succeeds.

---

## Pitfalls and Edge Cases

### 1. Same-DomainType sidesets

When both sides of a thin-shell interface have the same DomainType (e.g., two Conductor blocks), DomainType comparison alone cannot determine master/slave. The BFS propagation from cross-type seeds is essential in 3-D; in 2-D the lower-element-id master (see the rule above) is what the cut pipeline relies on.

**Risk:** If a connected component of same-type facets has no cross-type seed, BFS never
visits it in 3-D: only cross-type facets enter the queue, same-type facets are merely
flagged, and the final pass just clears their flags. Such a component keeps the orientation it came in
with from `connect_facets_to_elements` (lower element id is master), whose consistency BFS does not
establish. If that convention is wrong relative to other parts of the mesh, results may be incorrect.

### 2. Overhanging thin shells

A thin shell may overhang a conductor boundary. Some facets may then lie between Conductor/Air and others between Conductor/Conductor. The slave block collection must iterate over **all** facets to capture every slave block ID.

**Historical bug:** The original code used only the first facet's slave block ID, which failed for overhanging geometries. See `todo/orientation_fix_block_ordering_bug.md` (Bug 1).

### 3. Stale node-facet connections after relink

After `relink_slave_elements_with_duplicate_nodes()` replaces volume-element nodes with duplicates, and `update_facet_nodes()` copies master-element nodes into facet surface elements, existing `Node::facet()` pointers can become stale. This is an intermediate state: it lasts until `Mesh::finalize()` rebuilds the node-to-facet connectivity (`connect_nodes_to_facets()`, which runs after `update_facet_nodes()`). Code that runs in between sees:

- Some facet surface elements now reference duplicate nodes (if the master element was relinked)
- Original master nodes may have `number_of_facets() == 0` because their facets now reference different nodes
- Calling `tNode->facet(f)` on these nodes returns nothing useful

**Fix:** Build node-to-facet adjacency directly from the facet list rather than relying on `tNode->facet()`. See `todo/orientation_fix_block_ordering_bug.md` (Bug 3).

### 4. NaN normals from zero-length vectors

When ThinShellFactory's `process_nodes_tri3()` / `process_nodes_tri6()` / `process_nodes_line2()` compute surface normals, nodes with no valid facet connections produce a zero normal vector `tN = [0, 0, 0]`. Normalizing it (`tN /= norm(tN)`) produces NaN values, which then propagate to node coordinates in `create_nodes_on_layers()`.

**Symptom:** LAPACK `posv` error in the MaxwellPostprocessor recovery pass due to NaN entries in the Vandermonde matrix.

### 5. Non-orientable surfaces

A non-orientable surface admits no consistent orientation. The Stage 1 BFS does not detect this: it visits each flagged facet once and never re-checks a resolved one. For a thin-shell sideset in 3-D, the patch pass of `align_thin_shell_facets()` detects the resulting parity conflict and refuses the sideset with a `BELFEM_ERROR`; other sidesets are not checked. In practice, physically meaningful thin-shell surfaces are orientable (they represent real material layers).

### 6. Thin shell duplicate nodes and `original()`

Duplicate nodes created by `duplicate_nodes_on_face_sidesets()` do **not** have `set_original()` called on them. This means `tDuplicate->original()` returns `this` (the duplicate itself), not the original node.

`Facet::compute_orientation()` uses `original()->id()` to match nodes between master and slave faces. Since thin shell duplicates have different IDs from their originals, this matching still works correctly when duplicates appear only on the slave side (the master face has originals). Code that assumes `original()` tracks the duplication relationship will therefore fail for thin-shell duplicates.

---

## Execution Order in MaxwellFactory

These operations occur in the following order during `MaxwellFactory` initialization:

```
1. Mesh::finalize()                            — initial connectivity
2. MaxwellFactory::fix_facet_masters()          — correct master/slave by DomainType + BFS (BFS in 3-D only)
3. CutFactory::run()
   a. duplicate_nodes_on_face_sidesets()        — create thin shell duplicate nodes
   b. relink_slave_elements_with_duplicate_nodes() — replace slave element nodes
   c. restore_thin_shell_sidesets()             — restore original facets
   d. Mesh::finalize()                          — update_facet_nodes() overwrites facet nodes from master (skips thin-shell sidesets)
4. ThinShellFactory                             — compute normals, build layered elements
```

The key interaction is among steps 2, 3b, and 3d: `fix_facet_masters()` determines which side is master, `relink_slave_elements_with_duplicate_nodes()` gives duplicates to the slave side, and `update_facet_nodes()` copies from the master (which has originals) into the surface elements of the non-thin-shell interface facets.

---

## References

- `todo/orientation_fix_block_ordering_bug.md` — detailed report on the three bugs fixed in this area
- Messe et al. 2023 — BELFEM core paper, thin shell formulation
- Alves et al. 2022b — thin-shell theory, cohomology cuts, interface conditions
