# Jacobian Init — Remaining Bottlenecks After the Bitset Fix

**Date:** 2026-08-18
**Purpose:** Track the findings from the 2026-08-18 Jacobian-init investigation that are still
             open after the `DynamicBitset` and `SpMatrix` work landed
**Module:** `fem/kernel` (DofManager / SolverData), `math/graph`
**Status:** Open, not scheduled. The headline bottleneck is fixed; what remains is ~23 s on the
            magnetic system and is only worth attacking if that 23 s starts to matter.

## Where things stand

`devlog/dl20260818_jacobian_init_quadratic.md` found five candidate bottlenecks (C1-C5) in
`DofManager::init_matrices()`. Production result after the two changes that landed:

| phase | before | after |
|---|---:|---:|
| Jacobian init, magnetic | 290 s | 23 s |
| Jacobian init, thermal | 30.4 s | 3.4 s |
| assembly (first sample) | 9.87 s | 7.74 s |

**C1 is fixed** (`dl20260818_dynamicbitset_summary_bitmap.md`) and accounts for essentially the
whole init win. The `SpMatrix` accessor work (`dl20260818_spmatrix_accessor.md`) accounts for
the assembly figure, not the init one — it was a separate spec, not C2.

## Still open, in the order they are likely to matter

| id | finding | site | measured / estimated | status |
|---|---|---|---|---|
| C2 | `mDofData->dof( id )` is an `unordered_map` lookup once per DOF and once per adjacency entry, ~6.5e7 hash probes per `populate_graph` pass | `cl_FEM_DofMgr_SolverData.cpp:816`, `:830`; `cl_FEM_DofMgr_DofData.hpp:406-410` | **measured 5.6 s per pass** at 2.1e6 DOFs on a dev machine; ~22 s over four passes | open |
| C3 | `compute_dof_dof_connectivity` and `unite_dofs` each build `Cell< Vector<id_t> >` of size N and `unique()` every row; under Armadillo `unique()` reallocates per row. Rank 0 already holds the complete graph, so `unite_dofs` is a second serial pass over redundant worker copies | `cl_FEM_DofMgr_SolverData.cpp:1297-1342`, `:1424-1473`; `fn_unique.hpp:33` | est. 8-20 s | open |
| C5 | `graph::symrcm`'s disconnected-component fallback rescans all remaining vertices per component, O(C*V). Harmless on the connected free-DOF graph; the **fixed**-DOF graph is built with `aLinkToSelf = false` and can be highly fragmented | `fn_Graph_symrcm.cpp:139-158`, called from `cl_FEM_DofMgr_DofData.cpp:3319` | conditional; depends on n_fixed and boundary-patch connectivity, neither known for the production run | open |
| C4b | `create_assembly_tables` ships the System COO twice — `Jacobian` is a child of `System` sharing indices, but both sit in the `tNumMatrices = 6` loop | `cl_FEM_DofMgr_SolverData.cpp:605-663` | a few seconds of duplicated GB-scale collect | open |

C4's original claim (the `SpMatrix::index()` binary searches) was **refuted by measurement**
(0.26 s for 6.3e7 lookups) and is closed.

## Two things that were never measured

- **`Kernel::partition_mesh`** (`cl_FEM_Kernel.cpp:255-301`) uses the same per-element
  `reset()` + `where()` shape over the element graph that made `populate_graph` quadratic, so it
  **got faster for free** with the bitset change — by how much is unknown. Check the
  "Partitionig mesh" timer old branch vs new. This may now be the largest serial phase in
  startup, in which case it is the cheapest thing on this page to learn about.
- **The 1.83x `SpMatrix::multiply` improvement** lands in solve time, not init or assembly.
  Unmeasured in production.

## Whether to do any of this

The remaining 23 s is a one-time cost against runs that take hours, so the honest answer today
is probably not. The order of attack if it ever changes:

1. Measure `partition_mesh` first — zero implementation cost, and it may be larger than
   everything below.
2. C2 is the largest measured item and the cleanest fix: rewrite the neighbor IDs in
   `aGraphData` to `my_index` once, after which `populate_graph` stops hashing entirely.
3. C3 needs a design decision, not just an edit — the redundant `unite_dofs` pass exists
   because rank 0 keeps the full mesh, and removing it touches the parallel contract.
4. C5 is a latent risk rather than a known cost. Cheap diagnostic: log the component count
   from `symrcm( aFixedDofs )` on a production run.

## Related

- `devlog/dl20260818_jacobian_init_quadratic.md` — the original three-AI investigation
- `devlog/dl20260818_dynamicbitset_summary_bitmap.md` — C1 fix
- `devlog/dl20260818_spmatrix_accessor.md` — SpMatrix Phase A+B, and Phase C deferred behind a
  profile gate
