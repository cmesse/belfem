# Closed thin-shell side curve aborted setup

**Date:** 2026-10-03
**Purpose:** record the fix for a setup abort when a 3D thin shell's outer edge is a closed loop
**Module:** mesh (CurveFactory)

## Trigger

NI pancake deck with copper leads instead of tape leads
(`NI_coi/demo_clean_dummy_contact_test_nocopperseam_copperlead`, mesh from
`geo_thick/ni_coil_flip_nocopperseam_fuseculead.py`): the shell keeps only the spiral and the two
lead-base arcs, so no part of it reaches the domain boundary. Setup aborted:

```
terminate called after throwing an instance of 'std::out_of_range'
  what():  Mat::operator(): index out of bounds
  arma::arma_stop_bounds_error
  mesh::CurveFactory::create_protocurves
  mesh::CurveFactory::thin_shell_side_curves
  mesh::CutFactory::create_side_curves_for_thinshells_3d
  ...
  fem::MaxwellFactory::create_cuts
```
(gdb backtrace on a release binary, run on a scratch copy of the deck.)

## Cause

`thin_shell_side_curves` builds the side curves from the edges with one facet of the combined
shell; `select_edges` drops edges on domain-boundary sidesets. With tape leads, those dropped edges
at the terminals split the outer edge into open chains. Without them it is one closed loop.

`identify_subchains` assigned an edge its chain index only when the edge led to a node not yet in
the chain. In a closed loop, the last edge joins two nodes already in the chain and kept no index;
`create_protocurves` then indexed `tEdgeCount( tEdge->index() - 1 )` out of bounds. The later steps
already handle closed loops (`sort_end_nodes` sets the closed flag and picks the node with the
smallest id as the start; `sort_segments` stops when the loop closes).

## Change

`src/mesh/cl_CurveFactory.cpp`, `CurveFactory::identify_subchains`: every edge reached from a node
of the current chain gets the chain index, unconditionally. Open chains are unchanged (their edges
receive the same index, possibly more than once). Both callers (`thin_shell_side_curves` and
`create_protocurves`) run the same function on the same nodes.

## Files changed

| File | Change |
|---|---|
| `src/mesh/cl_CurveFactory.cpp` | edge index assignment in `identify_subchains` |
| `CHANGELOG.md` | `[Unreleased]` → Fixed entry |
| `devlog/dl20261003_closed_side_curve.md`, `devlog/README.md` | this record |

## Evidence

- Static: `cl_CurveFactory.cpp` passes `-fsyntax-only` (mpicxx from `/opt/scls/mkl`, `-std=c++17
  -DBELFEM_ARMADILLO`, this tree's `src/` directories and the SCLS headers as include paths). This
  tree has no build directory, so the module's own `compile_commands.json` flags were not used.
- **Not run:** no rebuild, no rerun of the deck. Reviewed, not verified. The gate is a build of this
  tree and the copper-lead deck getting past `create_cuts`. A closed side curve has not been
  exercised downstream before, so later steps may still need attention.
