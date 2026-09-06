# Devlog 2026-06-01 - Periodicity Factory Plan

**Date:** 2026-06-01
**Topic:** Read-only trace of periodicity factory lifecycle and next implementation plan
**AIs involved:** Codex
**Codex Audit Confidence:** medium-high for lifecycle and blockers; medium for final face-orientation details until tested
**Literature References:** N/A

## Summary

Traced the current `periodic_new` state around `PeriodicityFactory::create_periodicity()`, `CutFactory::run()`, `Periodicity::update()`, and the Maxwell periodic rebuild path. No source code was changed.

## Key Findings

- `periodic_new` already contains the Step 0 periodic cohomology/input commits (`3940461`, `2ca7a04`, `cec4fb7`, `75d35a7`) and the later `fd46a49` work-in-progress periodicity module.
- `PeriodicityFactory::create_periodicity()` is currently unfinished: it maps facets/nodes, starts edge collection, then exits at `src/mesh/cl_Mesh_PeriodicityFactory.cpp:107`.
- `Periodicity::match_faces()` is also unfinished/debug-only in the current tree: it uses hard-coded element IDs and exits at `src/mesh/cl_Mesh_Periodicity.cpp:325-459`.
- The intended lifecycle is two-phase: create initial node/edge/face periodicity before cohomology so `SimplicialComplex` can fold the periodic boundary, then rebuild final edge/face periodicity after cohomology cleanup and Maxwell edge/face recreation.
- `CutFactory::run()` creates temporary full-mesh edges/faces, calls `create_periodicity()`, computes cohomology, then resets edges/faces at `src/homology/cl_CutFactory.cpp:180-183`.
- `MaxwellFactory::create_magnetic_kernel()` later recreates conductor-needed edges/faces and calls `collect_nodes_from_flags_12()` plus `update()` at `src/fem/maxwell/cl_MaxwellFactory.cpp:447-450`.
- The current `flag_periodic_entities_12()` call in `MaxwellFactory::create_cuts()` occurs before `CutFactory::run()` creates `mMesh->periodicity()`, so periodic role flags need to be set after successful initial periodicity creation.

## Proposed Plan

1. Finish `PeriodicityFactory::create_periodicity()` as the initial geometry matcher: select source/target facets, pair facets by projected centroids, pair nodes by in-plane coordinates, construct a `Periodicity`, populate master/slave node and sideset lists, call `update()`, set role flags, and attach it to the mesh.
2. Move derived edge/face pairing responsibility into `Periodicity::update()` only. Remove or ignore the incomplete factory-local edge-map path.
3. Clean up `Periodicity::match_faces()` by removing hard-coded debug output/exit and validating it against the existing branch version, keeping the newer conductor/thin-shell edge filtering in `match_edges()`.
4. Ensure flags 1/2 are established immediately after the initial `update()` so `CutFactory` duplicate propagation can preserve master/slave roles.
5. Smoke-test `PeriodicityFactory` through a small executable or the current Maxwell input path before touching input syntax.
6. After mesh-level periodicity works, re-run or extend the existing DOF periodic tests and only then revisit explicit sideset filters in `todo/periodic_input_extension_plan.md`.

## Claude Consolidation (2026-06-01)

**Claude Confidence:** high for the code-verified current state; medium for the pre-cohomology edge-scope question until tested.

Re-verified all six of Codex's lifecycle claims against the `periodic_new` tree; all
confirmed:

- `create_periodicity()` `exit(0)` at `cl_Mesh_PeriodicityFactory.cpp:107-108`; `Periodicity`
  construction commented out :110-144; `map_edges()` defined :388 but never called.
- `match_faces()` debug/WIP: hard-coded elements 26903/25950 + debug `cout`, `exit(0)` at
  `cl_Mesh_Periodicity.cpp:459`. `fix_face_slaves()` :462 looks complete.
- `CutFactory::run()` creates full-mesh edges/faces :123-129, calls `create_periodicity()`
  :134, resets edges/faces :180-183.
- `MaxwellFactory::create_magnetic_kernel()` rebuild path: `collect_nodes_from_flags_12()` +
  `update()` at `cl_MaxwellFactory.cpp:447-450`.
- `flag_periodic_entities_12()` at `create_cuts()` :764-767 runs **before** the periodicity
  exists → currently a no-op.
- `match_edges()` mode 0 is conductor-scoped (`cl_Mesh_Periodicity.cpp:208-222`) → the
  pre-cohomology edge-scope open question is real (cuts traverse air).

**Reconciliation with the 2026-05 Claude notes:** the earlier claim that
`create_periodicity()` is "never called" was wrong — it is invoked from `CutFactory::run()`;
it is unfinished, not unwired. The DOF-side entanglement/sequencing work
(`cl_FEM_DofMgr_DofData.cpp`) is, by contrast, already complete and is the dead-but-correct
downstream of the mesh-level blocker.

Consolidated the Codex 6-step plan and the original Steps 0–8 into a single authoritative
"Consolidated Resume Plan (2026-06-01)" section at the top of
`todo/periodic_bc_fix_plan.md` (steps C1–C7 + the edge-scope open question). No source code
changed.

## Files Updated

- `devlog/dl20260601_periodicity_factory_plan.md`
- `todo/periodic_bc_fix_plan.md` (consolidated resume plan section)
- `devlog/README.md`
