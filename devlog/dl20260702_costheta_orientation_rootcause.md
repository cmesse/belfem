# Cos-theta quadrant: CW-mesh orientation root cause (edge-direction crash chain)

**Date:** 2026-07-02
**Purpose:** Root-cause the crash cascade blocking the cosθ iron quadrant model
(`cmake-build-debug/costheta.*`): facet-container overflow → edge-direction abort →
"discontinuous thin cuts". All traced to one mesh property.
**Module:** mesh, fem/kernel, homology

## Executive summary
`costheta.geo` defines the four air surfaces (72-75) with loop orientations that make
gmsh emit **19,655 clockwise triangles** (all of surfaces 72,73,74,75; cables and iron
are CCW). BELFEM's `MeshChecker` (cl_FEM_Kernel.cpp:94) **silently reorients** them
inside every first Kernel construction — *after* `CutFactory::run()` has already created
and linked edges (cl_CutFactory.cpp:123) and after facet `index_on_master` was assigned
at load. The flip invalidates every flipped element's edge-slot linkage (and any
orientation-derived structure), producing the observed downstream failures.

## Evidence chain (all measured, not inferred)
1. Non-unit cochain hypothesis **refuted** by instrumentation: `clean()` entry histogram
   shows `max|c| = 1` for all 2×68 generators (diagnostic in cl_Cohomology.cpp, #DIAG).
2. Edge-slot verification sweep (`check_element_edge_slots`, cl_CutFactory.cpp, #DIAG):
   0 mismatches after `create_edges()`, 0 at poisson entry, **39,310 after
   `fem::Kernel` ctor** = 19,655 flipped elements × 2 slots.
3. Python scan of costheta.msh: exactly **19,655 CW tri3**, all in surfaces 72-75.
   19,655 ≡ MeshChecker's flip count. Raw mesh has NO coincident nodes; max facet
   valence 3 — mesh is otherwise healthy.
4. With a CCW-corrected mesh (script-flipped, `costheta.msh.orig` = backup):
   0 mismatches at every sweep, edge-direction crash gone, and the cohomology
   generators collapse from 1400-2200-edge serpentine supports to **~300 edges** —
   the "broken/discontinuous thin cuts" were the same root cause.

## Fixes applied this session
1. `cl_Mesh_ConnectivityCalculator.cpp` — `connect_facets_to_facets()` now resets facet
   neighbor containers at entry ("container is already allocated" when two kernels each
   compute connectivity; first exposed by the two-kernel poisson path).
2. `cl_MaxwellFactory.cpp` — restored the `optimize cuts` input parsing (rank 0 at the
   `algorithm` read; un-commented the slave-side read). The key had become dead code, so
   `optimize cuts : true` silently did nothing (log showed RCM instead of Poisson).
3. `cmake-build-debug/costheta.msh` — CW triangles flipped to CCW by script (backup
   `costheta.msh.orig`). NOTE: the .geo should be fixed properly (outer loop first /
   loop orientation for Plane Surfaces 72-75) and the mesh regenerated.
4. Diagnostics (#DIAG, to be removed after the model runs): clean() histogram
   (cl_Cohomology.cpp), edge-slot sweeps (cl_CutFactory.cpp).

## Root cause #2: stale element adjacency after cuts (SEGV in BlockData)
After the orientation fix, the run reached the magnetic kernel and SEGV'd in
`BlockData::collect_element_indices` (cl_FEM_DofMgr_BlockData.cpp:296,
`element(0)` on a null container). Chain, all verified in gdb + source:
- `Element::allocate_element_container(0)` (cl_Element.cpp:59) freed the container
  but reset `mNumberOfElements` only in the `aSize != 0` branch → stale count over
  a null array; `Element::element()` (cl_Element.hpp:765) is assert-free.
- `CutFactory::compute_element_adjacencies()` builds adjacency from `mMesh->faces()`
  — EMPTY in 2D — so in 2D it calls `allocate_element_container(0)` on every element
  (destroying the poisson-round adjacency) and then falsely stamps
  `Connectivity::ElementToElement` as valid.
- `connect_elements_to_elements` early-returns on that flag
  (cl_Mesh_ConnectivityCalculator.cpp:1021) → the magnetic kernel never recomputes.
- Never seen before because (a) the poisson path was dead code (nothing filled the
  containers before the cuts), and (b) in 3D faces exist and the function is correct.
- Related: `Mesh::unfinalize()` deliberately skips the element-connectivity reset
  ("we assume that we don't add elements") — an assumption the cut factory violates.

**Fixes:** `Element::allocate_element_container` now always resets the counter;
`compute_element_adjacencies` in 2D resets containers + `reset_connectivity(
ElementToElement)` instead of pretending to have computed it.

## Root cause #3: parallel fresh run — Poisson cut kernel not block-selected
`Kernel::partition_mesh()` aborts at cl_FEM_Kernel.cpp:219 ("No blocks selected")
on a fresh (no-.bfm) parallel launch. Chain (gdb-confirmed, 2 ranks):
- `CutFactory::compute_poisson_problem()` builds `fem::Kernel tKernel(&tParams)` with
  `tParams` never given a block selection — `select_blocks()` is called on the *IWG*
  afterwards, which the Kernel ctor's `partition_mesh()` cannot see (it reads
  `mParams->selected_blocks()`). The magnetic kernel does it right
  (cl_MaxwellFactory.cpp:531 selects on the params *before* `make_shared<Kernel>`).
- Only fires fresh + parallel + `optimize cuts:true`: with a .bfm the whole cut stage
  (hence the Poisson solve) is skipped; in serial `partition_mesh` is never called.

**Fix:** `compute_poisson_problem` now calls `tParams.select_blocks(phi_block_ids())`
before constructing the kernel (cl_CutFactory.cpp). No-op in serial.

**Still broken beyond that (NOT fixed):** with the fix, the parallel Poisson kernel
partitions, then aborts deeper in `ProtoMesh::create_t_matrices()` (Cell OOB on rank 1
— hanging-basis/T-matrix data mismatch when distributing the temporary Poisson kernel).
A separate, harder parallel-distribution bug in the optimize-cuts path.

**Resolution — Poisson cut-optimization removed entirely (2026-07-02).** RCM is
"almost as good, significantly faster" (Christian) and, unlike the Poisson path, builds
no kernel in the cut stage so it distributes cleanly. Removed: `mComputePoisson` (both
factories), `CutFactory::compute_poisson_problem` + `rearrange_nodes`/`rearrange_edges`
(Poisson-only), the ctor `aComputePoisson` flag, and the `optimize cuts` input key.
`CutFactory::run()` now always calls `compute_rcm_problem()`. The general-purpose
`IwgType::Poisson` IWG (heat/hangingnodes/poisson executables) is untouched. Verified:
serial AND parallel (2 ranks) each run the full 28-timestep sim, zero crashes.

## Open items
- **Framework hardening (important):** element orientation must be fixed ONCE at mesh
  load, before facets/edges/curves/cohomology are derived — not silently inside Kernel
  construction. Options: move MeshChecker's core to src/mesh and invoke right after
  element creation in the readers; make Kernel's MeshChecker verify-only (assert/log,
  no silent mutation, at minimum report its flip count).
- MeshChecker never reports (`mElementCount` incremented, never printed).
- Next crash layer (under investigation): raw SEGV after the 2nd cohomology (thin-cut
  machinery), caught by PETSc's signal handler; gdb backtrace pending.
- T1-T9 SPFA rectification plan (todo/thin_cut_nonunit_rectification_implementation.md):
  premise (non-unit coefficients) does NOT occur in this model; plan remains valid for
  meshes where it does, but is not the fix for costheta.
- uint8_t facet-capacity truncation in `Vertex::allocate_facet_container`
  (cl_Vertex.cpp:214: mallocs full size, stores capacity mod 256): masked by the
  orientation fix in this model but still a silent-wrap landmine; guard or widen.

## Collaboration note
Hypothesis discipline paid off twice: the user demanded measurement before the SPFA
surgery (refuted), and the bit-identical crash under a supposedly changed configuration
exposed the dead `optimize cuts` parsing.
