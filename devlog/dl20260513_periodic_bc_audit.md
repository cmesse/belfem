# Devlog 2026-05-13 - Periodic BC Audit

**Date:** 2026-05-13
**Topic:** Periodic boundary condition audit across cohomology, thin shells, DOF conversion, parallel distribution, and input wiring
**AIs involved:** Codex
**Codex Audit Confidence:** high for DOF/input gaps; medium-high for cohomology and base thin-shell propagation
**Literature References:** N/A

## Correction

Follow-up branch comparison showed that the cohomology conclusion below was incomplete for the active `periodic_new` branch. The audit observed mesh-level periodic flagging/rebuild hooks and CutFactory duplicate propagation on `periodic_new`, but did not compare against the older `periodic` branch. The periodic-aware `Cohomology::clean()`, `SimplicialComplex`, `CutData`, and Maxwell input-periodicity setup exist on `periodic` and are missing from `periodic_new`.

Corrected confidence: high that `periodic_new` still needs the periodic cohomology/input machinery ported before the DOF sequencing and MPI fixes can be validated.

## Scope

Read-only audit requested while pausing side-connector work. The audit checked whether periodic boundary conditions are present in:

- cohomology/cut generation,
- thin-shell layer generation, excluding side connectors,
- conversion from periodic mesh entities to DOF constraints,
- MPI handling,
- user input / mesh-reader paths.

No source code was modified.

## Findings

### Cohomology and Cuts

Periodic mesh entity handling is present before and after cut generation, provided a `mesh::Periodicity` object already exists before `MaxwellFactory::create_cuts()`.

- `MaxwellFactory::create_cuts()` flags periodic entities before constructing cuts.
- `CutFactory` contains explicit propagation of periodic pointers/flags when duplicating nodes for thin-shell cuts.
- After thin-shell generation, `MaxwellFactory::create_magnetic_kernel()` calls `collect_nodes_from_flags_12()` and `update()` to rebuild periodic node/edge/face lists.
- `Topology::select_sidesets()` includes `AirPeriodic`, `BufferPeriodic`, and `FerroPeriodic` in phi boundary sidesets.

Confidence: medium-high. The expected pipeline exists, but it depends on periodicity being initialized upstream, and one duplicate-linking path appears fragile around unflagging before flag-copy logic.

### Thin Shell Factory

Base thin-shell layer periodicity is handled for nodes, with edges/faces derived later by the periodicity rebuild.

- `ThinShellFactory::create_nodes_on_layers()` maps periodic source nodes to corresponding layer nodes and sets reciprocal periodic pointers.
- The same code copies periodic master/slave flags from original nodes to layer nodes.
- `Periodicity::update()` later derives periodic edges and faces from node periodicity; it has explicit thin-shell edge handling by layer.

Confidence: high for base thin-shell layers. Side connectors remain excluded from this conclusion; their current periodic handling is incomplete and ad hoc.

### DOF Conversion

There is a DOF conversion routine, but it is sequenced incompletely.

- `DofData::create_dofwise_periodicities_master()` maps periodic node/edge/face DOFs by setting slave DOF sources to the master DOF or to the master's source chain.
- However, it is called after the initial `mHangingDOFs` collection. The periodic slave DOFs it marks are not registered back into `mHangingDOFs`.
- Serial assembly may partially work because element consolidation tests `Dof::is_hanging()` directly, but counts, lifetime ownership, hanging-DOF field updates, and any code iterating `hanging_dofs()` can miss periodic slaves.

Confidence: high.

### Parallel Handling

Mesh-level periodicity distribution exists, but DOF-level periodic constraints are not distributed correctly.

- `Mesh_Distributor` serializes periodic node/edge/face IDs.
- `ProtoMesh::create_periodicities()` reconstructs periodic pointers on non-root ranks.
- The parallel hanging-DOF communication path sends only DOFs in `mHangingDOFs`. Since periodic slave DOFs are not appended after conversion, non-root ranks do not receive the corresponding periodic DOF source relations.

Confidence: high.

### Input Status

The input path appears unwired for periodic BCs.

- `PeriodicityFactory` exists and is demonstrated only in `src/mesh/main.cpp`.
- No production Maxwell input path was found that creates a `mesh::Periodicity` from an input section.
- `domain_type(string)` does not parse periodic domain strings.
- `Domain::Domain()` does not construct periodic sideset domain types from input sections.
- No Gmsh `$Periodic` reader support was found.
- HDF5/BFM persistence of `mesh::Periodicity` was not found.

Confidence: high.

## Proposed Plan

1. Wire input periodicity first. Define the supported input syntax for periodic master/slave planes or sidesets, create the `mesh::Periodicity` before `MaxwellFactory::create_cuts()`, and add periodic domain-type parsing if periodic sidesets are intended to come from input sections.
2. Fix DOF registration sequencing. Keep mesh-basis hanging conversion first, apply periodic DOF source relations, then rebuild `mHangingDOFs` from all `mDOFs` where `Dof::is_hanging()` before serial removal or MPI distribution.
3. Handle preexisting slave sources explicitly. The periodic conversion currently assumes slave DOFs do not already have sources; the fix should assert a clear invariant or implement a controlled source replacement/merge policy.
4. Add serial regression coverage with periodic thin shells. Verify that original nodes, cut duplicates, and thin-shell layer nodes participate in periodicity; verify derived periodic edges; keep side connectors out of this test.
5. Add MPI regression coverage. Use a small periodic case spanning ranks and assert that non-root periodic slave DOFs receive sources and element consolidation eliminates the slave global DOFs.
6. Add input/end-to-end coverage once input wiring exists. The test should prove that a user-declared periodic configuration reaches cohomology, thin-shell rebuilding, DOF conversion, and MPI distribution.

## Canonization

2026-05-13: Claude (Opus 4.7, 1M ctx) reviewed Codex's findings against live source and confirms agreement on every section. The two audits diverged only on plan ordering — Codex put input wiring first, Claude put DOF sequencing first. Canonical order is **DOF sequencing first**, on the rationale that the procedural `PeriodicityFactory` API can drive tests without input wiring, so the actual bug should be fixed and verified before widening the user-facing surface.

Canonical step-by-step plan lives in `todo/periodic_bc_fix_plan.md`.

## Files Touched

- `devlog/dl20260513_periodic_bc_audit.md`
- `devlog/README.md`
- `todo/periodic_bc_fix_plan.md` (added 2026-05-13)
- `todo/README.md` (index entry added 2026-05-13)
