# Thin Shell Facet Orientation {#fem_maxwell_thin_shell_facet_orientation}

**Date:** 2026-03-02
**Purpose:** Documents why facet master/slave orientation matters for thin shells, the two-stage algorithm that ensures consistency, and the pitfalls that arise from mesh operations
**Module:** fem/maxwell, homology, mesh

---

## Key Files

| File | Key Function | Purpose |
|------|-------------|---------|
| `cl_MaxwellFactory.cpp` | `fix_facet_masters()` (~line 981) | Corrects master/slave based on DomainType + BFS |
| `cl_CutFactory.cpp` | `relink_slave_elements_with_duplicate_nodes()` (~line 1838) | Replaces slave element nodes with thin shell duplicates |
| `cl_CutFactory.cpp` | `duplicate_nodes_on_face_sidesets()` (~line 1646) | Creates the duplicate nodes |
| `cl_Mesh.cpp` | `update_facet_nodes()` (~line 678) | Copies master element's face nodes into the facet surface element |
| `cl_Mesh_ConnectivityCalculator.cpp` | `connect_facets_to_elements()` (~line 171) | Initial master/slave assignment by element ID |
| `cl_Facet.cpp` | `flip()` (~line 89) | Swaps master and slave, relinks facet nodes |
| `fn_check_facet_orientation.hpp` | `check_facet_orientation()` (~line 39) | Checks if two adjacent facets have consistent outward normals |
| `en_DomainType.hpp` | `DomainType` enum (~line 19) | Defines priority ordering: Air(1) < Ferro(2) < Coil(3) < Conductor(4) |

---

## Problem Statement

Thin shells in BELFEM are surface-embedded conductors (superconducting tapes, thin films) modeled as sidesets rather than volume blocks. They require **topological separation**: the mesh must have duplicate nodes on one side of the shell so that DOFs on the two sides are independent. Which side gets the duplicates must be **consistent across all facets** of a given thin shell sideset.

Three mesh operations interact to make this work:

1. **`update_facet_nodes()`** copies nodes from the facet's **master** volume element into the facet surface element. Therefore, the facet master must always be the element that keeps **original** nodes.

   > **It skips thin-shell facets.** `Mesh::update_facet_nodes()` returns early for sidesets whose
   > `domain_type()` is `ThinShell` or `GeometryOnly` (`cl_Mesh.cpp:687-690`), because those are
   > extrusion geometry whose node lists are built by the `ThinShellFactory` rather than copied
   > from a volume element. The orientation rule below still governs the *interface* facets between
   > the shell and the bulk, which is where master/slave actually matters.

2. **`relink_slave_elements_with_duplicate_nodes()`** replaces interface nodes in **slave** block elements with thin shell duplicates. Therefore, the slave side must be the side that receives duplicates.

3. **`ThinShellFactory`** reads facet nodes via `tNode->facet()` to compute surface normals and build layered thin shell elements. It expects the facet surface element to reference original nodes with valid facet connections.

If master/slave assignment is inconsistent, the wrong side gets duplicates and ThinShellFactory computes incorrect normals (or NaN normals from zero-length vectors). (`update_facet_nodes()` is not the failure path for the thin-shell sidesets themselves — it skips them, per the note above — but it is for the interface facets between shell and bulk.)

---

## Two-Stage Algorithm

### Stage 1: `MaxwellFactory::fix_facet_masters()` (cl_MaxwellFactory.cpp:981-1049)

#### Initial state

After mesh loading, `connect_facets_to_elements()` assigns master/slave purely by element ID: the element with the **lower ID** becomes the master. This is arbitrary and depends on block ordering in the mesh file.

#### DomainType-based correction

`fix_facet_masters()` overrides this assignment based on `DomainType` priority. The DomainType enum values control the ordering (en_DomainType.hpp:19-73):

```
Air(1) < Ferro(2) < Coil(3) < Conductor(4) < Cut(5) < ThinShell(6)
```

The **higher** DomainType value becomes master. For a typical Air/Conductor interface, Conductor (4) becomes master and Air (1) becomes slave.

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

These cross-type facets are added to a BFS queue as **seeds** with known correct orientation.

#### Same-type facets (BFS propagation)

Facets where both sides have the **same** DomainType (e.g., Conductor/Conductor) cannot be resolved by DomainType comparison alone. These are flagged for later processing.

BFS propagation from the cross-type seeds resolves them:

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

`check_facet_orientation()` compares the traversal direction of the shared edge between two facets. If both facets traverse the shared edge in the **same** direction, their outward normals point in opposite directions and a flip is needed.

#### Persistence across finalize cycles

The master/slave assignment set by `fix_facet_masters()` is preserved across `unfinalize()`/`finalize()` cycles because `FacetToElement` connectivity is not reset during those operations.

---

### Stage 2: `CutFactory::relink_slave_elements_with_duplicate_nodes()` (cl_CutFactory.cpp:1838-1910)

This function runs after `duplicate_nodes_on_face_sidesets()` has created the duplicate nodes. It relinks elements on the slave side to use duplicates instead of originals.

#### Per-facet slave block collection

The function collects slave block IDs from **all** facets of each thin shell sideset, not just the first one:

```cpp
for ( Facet * tFacet : tFacets )
{
    tBlockBitset.set( mMesh->block( tFacet->slave()->block_id() )->index() );
}
```

This is critical for overhanging thin shells where different facets may sit between different block pairs.

#### Element flagging

All elements in the identified slave blocks are flagged. Then, only the flagged elements that actually touch sideset nodes are collected:

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

For each collected element, interface nodes that belong to the thin shell sideset are replaced with their duplicates:

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

Per-sideset node flagging ensures that nodes from different thin shell sidesets are not accidentally replaced in the wrong context.

---

## Why This Works

After both stages complete, the following invariants hold:

1. **`fix_facet_masters()`** ensures that for every thin shell facet, `master()` returns the volume element that will keep original nodes, and `slave()` returns the volume element that will receive duplicate nodes.

2. **`relink_slave_elements_with_duplicate_nodes()`** replaces interface nodes in slave-side elements with thin shell duplicates. Master-side elements keep originals.

3. **`update_facet_nodes()`** (called during `Mesh::finalize()`) copies nodes from the master element into the facet surface element. Since the master kept originals, the facet surface element gets original nodes.

4. **`ThinShellFactory`** reads facet nodes and computes surface normals. Because facet surface elements reference original nodes (from the master), the node-to-facet connections are valid, and normal computation succeeds.

---

## Pitfalls and Edge Cases

### 1. Same-DomainType sidesets

When both sides of a thin shell interface have the same DomainType (e.g., two Conductor blocks), DomainType comparison alone cannot determine master/slave. The BFS propagation from cross-type seeds is essential.

**Risk:** If the entire thin shell sideset has same-type blocks on both sides with no cross-type seeds anywhere in the connected component, the orientation is arbitrary but still consistent (all facets will be oriented the same way by BFS). However, if this leads to the wrong convention relative to other parts of the mesh, results may be incorrect.

### 2. Overhanging thin shells

A thin shell may overhang a conductor boundary, so that some facets sit between Conductor/Air and others between Conductor/Conductor. The slave block collection must iterate over **all** facets to capture every slave block ID.

**Historical bug:** The original code used only the first facet's slave block ID, which failed for overhanging geometries. See `todo/orientation_fix_block_ordering_bug.md` (Bug 1).

### 3. Stale node-facet connections after relink

After `relink_slave_elements_with_duplicate_nodes()` replaces volume element nodes with duplicates, and `update_facet_nodes()` copies master element nodes into facet surface elements, the original `Node::facet()` pointers can become stale:

- Some facet surface elements now reference duplicate nodes (if the master element was relinked)
- Original master nodes may have `number_of_facets() == 0` because their facets now reference different nodes
- Calling `tNode->facet(f)` on these nodes returns nothing useful

**Fix:** Build node-to-facet adjacency directly from the facet list rather than relying on `tNode->facet()`. See `todo/orientation_fix_block_ordering_bug.md` (Bug 3).

### 4. NaN normals from zero-length vectors

When ThinShellFactory's `process_nodes_tri3()` / `process_nodes_tri6()` / `process_nodes_line2()` compute surface normals, nodes with no valid facet connections produce a zero normal vector `tN = [0, 0, 0]`. Normalizing this (`tN /= norm(tN)`) produces NaN, which propagates to node coordinates in `create_nodes_on_layers()`.

**Symptom:** LAPACK `posv` error in the MaxwellPostprocessor recovery pass due to NaN entries in the Vandermonde matrix.

### 5. Non-orientable surfaces

If BFS propagation encounters a facet that has already been resolved but with an inconsistent orientation, the surface is non-orientable. This should be detected and reported as an error. In practice, all physically meaningful thin shell surfaces are orientable (they represent real material layers).

### 6. Thin shell duplicate nodes and `original()`

Duplicate nodes created by `duplicate_nodes_on_face_sidesets()` do **not** have `set_original()` called on them. This means `tDuplicate->original()` returns `this` (the duplicate itself), not the original node.

`Facet::compute_orientation()` uses `original()->id()` to match nodes between master and slave faces. Since thin shell duplicates have different IDs from their originals, this matching still works correctly when duplicates appear only on the slave side (the master face has originals). But code that assumes `original()` tracks the duplication relationship will not work for thin shell duplicates.

---

## Execution Order in MaxwellFactory

The relevant operations occur in this sequence during `MaxwellFactory` initialization:

```
1. Mesh::finalize()                            — initial connectivity
2. MaxwellFactory::fix_facet_masters()          — correct master/slave by DomainType + BFS
3. CutFactory::run()
   a. duplicate_nodes_on_face_sidesets()        — create thin shell duplicate nodes
   b. relink_slave_elements_with_duplicate_nodes() — replace slave element nodes
   c. restore_thin_shell_sidesets()             — restore original facets
   d. Mesh::finalize()                          — update_facet_nodes() overwrites facet nodes from master
4. ThinShellFactory                             — compute normals, build layered elements
```

The key interaction is between steps 2, 3b, and 3d: `fix_facet_masters()` determines which side is master, `relink_slave_elements_with_duplicate_nodes()` gives duplicates to the slave side, and `update_facet_nodes()` copies from the master (which has originals) into the facet surface elements.

---

## References

- `todo/orientation_fix_block_ordering_bug.md` — detailed report on the three bugs fixed in this area
- Messe et al. 2023 — BELFEM core paper, thin shell formulation
- Alves et al. 2022b — thin-shell theory, cohomology cuts, interface conditions
