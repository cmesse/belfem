# Thin-Shell Postprocessor: Shared Interface Nodes {#fem_maxwell_thinshell_postprocessor_node_sharing}

**Date:** 2026-04-06
**Status:** Known limitation. **Math is correct; visualization only is affected.**
**Module:** `src/fem/maxwell` and `src/mesh`
**Related:**
- [postprocessor_recovery_theory.md](postprocessor_recovery_theory.md)
- [thin_shell_facet_orientation.md](thin_shell_facet_orientation.md)
- `devlog/dl20260406_thinshell_postprocessor_trace.md`

---

## TL;DR

In a thin-shell stack with adjacent layers of *different* materials (e.g.
`hastelloy / ybco`), the H-formulation DOFs are decoupled correctly across the
layer interface — edges and faces are duplicated, ghost facets are inserted,
and the Nitsche-DG-Ghost penalty couples the two layers only weakly. **The
solver is fine.**

However, the *geometric nodes* sitting on that layer interface are still
shared between the two layer blocks. The Maxwell postprocessor performs nodal
SPR recovery and writes its outputs into a global per-node field. As a result,
fields that are only meaningful on one layer (e.g. `J/Jc` on the YBCO layer)
appear to "leak" onto the other layer (e.g. Hastelloy) when viewed in
ParaView, even though no leakage actually exists in the assembled system.

This document explains the trace, why we are deliberately not fixing it right
now, and what a clean fix would look like when we revisit.

---

## 1. What we observed

Setup (from `cmake-build-debug/input.conf`):

```
layers : tape
{
    hastelloy : 50 mum ;
    ybco      : 1.0  mum ;
}
```

Two `tape` thin shells, each with one Hastelloy layer block and one YBCO layer
block. After running, the recovered `J/Jc` field is non-zero on the
Hastelloy-layer side of the layer interface — exactly at the row of nodes that
sits between Hastelloy and YBCO.

The H field, B field, and assembled DOF system show no equivalent
contamination.

---

## 2. Trace

### 2.1 Layer node sharing in `ThinShellFactory`

`ThinShellFactory::create_nodes_on_layers` (`src/mesh/cl_ThinShellFactory.cpp:980-1134`)
allocates one `Cell<Node*>` *per layer*, in `Layer::Nodes`. The `Layer` struct
itself (`src/mesh/cl_ThinShellFactory.hpp:35-44`) only carries `Edges`,
`EdgeDuplicates`, `Faces`, `FaceDuplicates`, and `GhostFacets`. **There is no
`NodeDuplicates`.**

When element blocks are built (e.g.
`create_elements_on_blocks_line2` at `src/mesh/cl_ThinShellFactory.cpp:1140-1177`,
and the analogous `line3` / `tri3` / `tri6` variants), the bottom and top
nodes of block `b` come straight from the layer node arrays:

```cpp
Layer * tBottom = aLayers( b );
Layer * tTop    = aLayers( b + 1 );
...
tElement->insert_node( tBottom->Nodes( ... ), 0 );
tElement->insert_node( tBottom->Nodes( ... ), 1 );
tElement->insert_node( tTop   ->Nodes( ... ), 2 );
tElement->insert_node( tTop   ->Nodes( ... ), 3 );
```

For `[hastelloy, ybco]`:

| Block          | bottom nodes        | top nodes           |
|----------------|---------------------|---------------------|
| 0 (hastelloy)  | `aLayers(0)->Nodes` | `aLayers(1)->Nodes` |
| 1 (ybco)       | `aLayers(1)->Nodes` | `aLayers(2)->Nodes` |

**The two blocks share the literal `Node*` objects in `aLayers(1)->Nodes`.**

### 2.2 What `hasDuplicates` actually duplicates

`Layer::hasDuplicates` is set in
`src/mesh/cl_ThinShellFactory.cpp:186-194` whenever the materials of two
adjacent layer blocks differ. Following the flag:

- **Edges** — `create_edges_on_layers` (`src/mesh/cl_ThinShellFactory.cpp:1463-1504`)
  builds a parallel `EdgeDuplicates` cell only when `hasDuplicates` is true.
  `link_elements_with_edges` (`src/mesh/cl_ThinShellFactory.cpp:1550-1623`)
  hands the upper block the duplicates as its bottom row, so adjacent blocks
  do not share edges.
- **Faces** — `create_faces_on_layers` (`src/mesh/cl_ThinShellFactory.cpp:1509-1547`)
  and `link_elements_with_faces` mirror the same pattern.
- **Ghost facets** — `create_ghost_facets`
  (`src/mesh/cl_ThinShellFactory.cpp:1656-`) inserts the Nitsche/ghost penalty
  facets only on `hasDuplicates` interfaces, since that is where the edge DOFs
  on each side are distinct.
- **Nodes** — there is no equivalent treatment. `create_nodes_on_layers` always
  produces exactly one node array per layer; adjacent blocks reuse it.

### 2.3 Why the H-system is still fine

The Maxwell H-φ formulation puts conductor DOFs on **edges**
(`edge_h`, plus `face_h` for higher-order) — see `IWG_Maxwell` field setup at
`src/fem/maxwell/cl_IWG_Maxwell.cpp:51`. There are no nodal conductor DOFs in
the assembled system, so the shared nodes do not introduce any extra coupling
into K. The Nitsche-DG-Ghost penalty acts on the duplicated edges and behaves
exactly as designed.

### 2.4 Why the postprocessor leaks

`MaxwellPostprocessor::compute_superconductor_ts`
(`src/fem/maxwell/cl_MaxwellPostprocessor.cpp:585-`) computes `J/Jc` per
integration point of each YBCO thin-shell element. Recovery happens through
the generic `Postprocessor` SPR pipeline:

- `Postprocessor::select_elements_and_owned_nodes`
  (`src/fem/kernel/cl_FEM_Postprocessor.cpp:222-333`) walks the postprocessor's
  blocks and registers, for each selected element, every
  `tElement->node(k)->original()` in the bitset of "owned nodes."
- `Postprocessor::recover_fields`
  (`src/fem/kernel/cl_FEM_Postprocessor.cpp:867-1034`) accumulates each
  element's contribution into those nodes via the standard SPR Vandermonde
  solve, and `synch_target_fields` writes the result into the global per-node
  field.
- `select_blocks_and_materials`
  (`src/fem/maxwell/cl_MaxwellPostprocessor.cpp:233-256`) only enrolls blocks
  whose material has `Jc`, so the Hastelloy block is *not* in the
  superconducting postprocessor.

The leak does **not** come from the Hastelloy block being processed. It comes
from the YBCO block writing into nodes that the Hastelloy block also points at.
For the bottom row of YBCO elements those nodes are literally
`aLayers(1)->Nodes`, which is the same memory the Hastelloy block uses for its
top row. So `JJCz(node_index)` ends up non-zero at indices that ParaView
colors on both sides.

The same mechanism applies to any per-node recovered field on a thin shell
where the layer materials differ. `J/Jc` is just the most visually obvious
because it is only meaningful on one side of the interface.

---

## 3. Decision: do not fix this right now

The clean fix is to introduce `Layer::NodeDuplicates` as a sibling of
`EdgeDuplicates` / `FaceDuplicates`, allocate them when `hasDuplicates` is
true, and route the upper block's bottom row through the duplicates. That
would also require:

- mirroring the cross-facet `set_original`/`add_duplicate` linking inside
  `create_nodes_on_layers` for the duplicates,
- handling periodicity for the duplicates (the periodic linking loop at
  `src/mesh/cl_ThinShellFactory.cpp:1084-1126`),
- re-collapsing the duplicates back to a single nodal DOF for the **thermal**
  module, which uses node-based Lagrange DOFs and currently relies on the
  shared nodes to keep the temperature continuous across layer interfaces.

We are deliberately not paying that cost today, because:

1. The math is correct — the assembled K matrix is unchanged. This is purely a
   visualization artifact.
2. The thermal module would need its own Nitsche-style treatment (or an
   explicit static condensation pass) if we duplicate the nodes, and we are
   planning to migrate the thermal coupling to a finite-volume-element method
   anyway. We do not want to invest in a temporary thermal Nitsche path that
   we will throw away.
3. Workaround: when reading the field in ParaView, threshold or block-by-block
   filter `JJCz` (and any other per-node field that should only live on one
   layer) by the block ID before coloring. The Hastelloy block label format
   added in 2026-04 (`<domain>_<id>_<material>`, see `MaxwellFactory::create_and_assign_materials`
   in `src/fem/maxwell/cl_MaxwellFactory.cpp`) makes this filter
   straightforward.

---

## 4. Future fix sketch

When we revisit (likely after the FVM thermal transition is in place), the
work is:

1. **`src/mesh/cl_ThinShellFactory.hpp`** — extend the `Layer` struct:

   ```cpp
   struct Layer
   {
       bool hasDuplicates = false ;
       Cell< Node *  > Nodes ;
       Cell< Node *  > NodeDuplicates ;     // NEW
       Cell< Edge *  > Edges ;
       Cell< Edge *  > EdgeDuplicates ;
       Cell< Face *  > Faces ;
       Cell< Face *  > FaceDuplicates ;
       Cell< Facet * > GhostFacets ;
   };
   ```

2. **`create_nodes_on_layers` in `src/mesh/cl_ThinShellFactory.cpp`** — when
   `tLayer->hasDuplicates`, allocate a parallel `new Node(...)` set with the
   same coordinates and store them in `NodeDuplicates`. Replicate the
   cross-facet `add_duplicate` / `set_original` linking and the periodicity
   wiring for the new nodes.

3. **`create_elements_on_blocks_line2 / line3 / tri3 / tri6`** — pick the
   bottom row the same way `link_elements_with_edges` already picks the bottom
   edges:

   ```cpp
   Cell<Node*> & tBottom = aLayers(b)->hasDuplicates
       ? aLayers(b)->NodeDuplicates
       : aLayers(b)->Nodes ;
   Cell<Node*> & tTop    = aLayers(b+1)->Nodes ;
   ```

4. **Thermal coupling** — at the time the fix lands, the thermal solver must
   either (a) have moved to FVM and reformulated the temperature continuity
   condition across the layer interface itself, or (b) have a static
   condensation / explicit constraint pass that re-merges the duplicate node
   pair at the temperature DOF level only.

5. **Postprocessor** — no changes required. Once the nodes are physically
   distinct, the existing SPR pipeline writes only into the appropriate
   layer's nodes and the visualization is clean automatically.

---

## 5. Verification checklist (for the future fix)

When the duplicate-node path is implemented, confirm:

- [ ] Assembled K matrix and solution are bit-for-bit unchanged on a problem
      with `hasDuplicates == false` everywhere.
- [ ] On a `hastelloy / ybco` stack, `JJCz` is exactly zero on every node that
      belongs only to the Hastelloy block (no SPR leak).
- [ ] Temperature field across the layer interface is still continuous (this
      is the thermal-side requirement that motivates the deferral).
- [ ] Cohomology cuts and ghost facets still resolve correctly.
- [ ] Periodic thin-shell stacks still satisfy the periodicity constraints on
      both the originals and the duplicates.

---

## 6. Provenance

This document records a joint Claude/Codex investigation on 2026-04-06.
Confidence: **high** that this is a postprocessor visualization artifact and
not a solver bug. See `devlog/dl20260406_thinshell_postprocessor_trace.md` and
the corresponding entry in `todo/ai_exchange.md` for the audit trail.
