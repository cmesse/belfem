# Periodic DOFs as Hanging DOFs — Parallel Architecture (Design Draft)

**Date:** 2026-06-04
**Purpose:** Capture the design decision that periodic DOFs are a special case of hanging DOFs, and what that implies for parallel (MPI) handling. Draft for review; to be folded into a `doc/` section once agreed.
**Status:** Draft / under discussion
**Related:** `todo/periodic_bc_fix_plan.md` (DOF sequencing, MPI distribution), `doc/` hanging-DOF / static-condensation notes

---

## Core insight

**A periodic DOF equality is an algebraic hanging/entanglement constraint.**

In the simple case, one DOF is constrained to the other with one source and weight ±1.
However, "periodic slave hangs on periodic master" is only the geometric intuition, not
the full algebraic rule:

- the existing DOF-level `entangle()` policy chooses a deterministic representative
  (currently smaller DOF ID) when neither side is hanging;
- if both sides are already hanging, it unfolds both source rows and merges them;
- exactly one side already hanging is treated as a topology error unless we explicitly
  decide to relax that invariant.

So the solve-time object is not a special periodic relation; it is a source row in the
same hanging-DOF graph used by thin shells and other T-matrices. In the trivial row this is
`target = source`, but the implementation must preserve the existing row-composition
semantics rather than blindly overwriting `Basis::sources()`.

Signs and orientations:

- node rows use `+1`;
- first-order edge rows use `+1` after `PeriodicityFactory::match_edges()` has aligned the
  slave edge orientation; without that normalization the orientation sign must be carried
  into the row;
- face ownership matters for cohomology orientation, but FEM face DOFs are a separate
  higher-order Nedelec issue. If/when those are enabled, periodic face DOF constraints must
  be represented algebraically as rows/T-matrices too; they should not be dismissed as
  "cohomology only."

This means periodic constraints can reuse the existing static-condensation stack:

- basis-source/T-matrix serialization,
- row unfolding,
- assembly in fully-condensed submatrix form on each proc,
- the cross-proc `mHangingDOFs` MPI exchange.

There is no need for a periodicity-specific assembly path. There is also no need for worker
procs to recompute geometric periodicity.

Important correction to older notes: on the active branch, `collect_hanging_dofs()` already
rebuilds `mHangingDOFs` after `create_dofwise_periodicities_master()` via the `tAfter` pass
in `cl_FEM_DofMgr_DofData.cpp`. The old claim that periodic slaves are missed merely
because the list was frozen is stale. The remaining distributed-design problem is how to
make workers see the same periodic-derived source rows without requiring a live geometric
`Periodicity` object on each worker.

## The one thing periodic adds over ordinary hanging DOFs

Ordinary hanging DOFs (h-φ interfaces, thin shells) have a source that is **geometrically
adjacent** — same interface, same elements — so the source DOF is already in the partition,
or pulled in by the adjacency-based Aura "for free." Nobody has to think about it.

A **periodic** hanging DOF's eventual algebraic source/representative may sit on the
**opposite face** of the domain. It is never adjacency-local. That is the *entire*
periodic-specific problem:

> **Make the periodic counterpart entity local to any proc that needs the constrained row.**

Once the counterpart basis exists locally, a periodic hanging/entangled DOF row is
indistinguishable from any other hanging row, and the existing condensation + exchange
handle it unchanged.

## Resolution: pairing before Aura, constraints after Aura

The distributor has two separate jobs:

1. **Before / during Aura selection:** use rank-0 periodic pairings as another adjacency
   relation. Opposite-face partners never appear in an ordinary connectivity halo, so the
   periodic pairing must drive ghost selection. The active `Mesh_Distributor::select_entities()`
   already follows periodic node/edge/face links while selecting the aura.
2. **After the final per-proc entity set is known:** serialize the periodic constraint rows
   for the entities that actually exist on that proc.

This is why the clean distribution answer is not "send a `Periodicity` object to every
worker." The worker needs the periodic partner as a mesh basis/source so the condensed
submatrix can be built correctly; it does not need the geometric periodicity machinery.

The preferred implementation is therefore option B from the design discussion:

- keep the root mesh's `Periodicity` as the geometric/pairing source;
- let the distributor use it to include periodic partners in the aura;
- then augment each proc's `proto::TMatrixData` (or an equivalent basis-source payload)
  with periodic-derived source rows;
- workers reconstruct ordinary `Basis::sources()` through `ProtoMesh::create_t_matrices()`;
- `DofData::collect_hanging_dofs()` then sees periodic constraints through the normal
  hanging path.

This keeps periodicity out of the worker solve contract. It also avoids the ambiguity of
rehydrating worker-local `periodic()` links and then running `create_dofwise_periodicities_master()`
again, which risks double-applying constraints or forcing worker meshes to support
geometric `Periodicity::update()`.

## Corollary: keep `PeriodicityData`, narrow the solve contract

Do **not** delete `proto::PeriodicityData`. It still has a legitimate semantic role:

- saving/loading periodic mesh metadata,
- carrying side-set/pairing information through `ProtoMesh`,
- reconstructing the root/full-mesh `Periodicity` when a persisted mesh is loaded.

But distributed assembly should not depend on worker-local geometric periodicity. For the
solve, periodicity should be absorbed into ordinary hanging/T-matrix source rows before the
worker builds its condensed contributions.

This gives two contracts:

- **Persistence contract:** `PeriodicityData` records that the mesh has periodic geometry.
- **Distributed solve contract:** `TMatrixData` / `Basis::sources()` records the algebraic
  constraints each worker needs.

`PeriodicityFactory` and geometric `Periodicity::update()` should remain serial/root-full-mesh
operations. A serialized worker-local `Periodicity` must not be allowed to trigger geometric
plane matching during `Mesh::finalize()`; local worker meshes do not necessarily contain
complete periodic planes.

## Confirmed exchange behavior

`cl_FEM_DofMgr_DofData.cpp`'s hanging-DOF exchange already resolves source DOFs by global
DOF ID. Workers report hanging DOFs, root returns the source DOF metadata, and workers
create missing source `Dof` objects as needed.

This confirms that source DOFs may be ghost-owned. The non-negotiable requirement is that
the **mesh basis** for the source DOF already exists in the worker mesh, because the worker
constructs missing source DOFs from `mMesh->node/edge/face/...` using the source mesh-basis
ID. Periodic aura selection must therefore include the periodic source entity before the
hanging exchange runs.

So periodic distribution needs:

1. periodic partner entity included in the aura;
2. periodic-derived source row serialized to the worker;
3. no worker-side geometric periodicity update.

## Implementation implications

1. Rank 0 builds full geometric periodicity before distribution.
2. The distributor uses periodic links to close the aura over periodic partner entities.
3. After aura selection, the distributor materializes periodic constraints into the same
   per-proc source-row format used by `populate_t_matrices()`.
4. Row materialization must preserve `entangle()` semantics:
   - neither side hanging: choose deterministic representative;
   - both sides hanging: unfold and merge rows;
   - exactly one side hanging: keep the current topology-error policy unless deliberately
     changed.
5. Workers consume only normal `TMatrixData` for assembly.
6. The worker path should not also run `create_dofwise_periodicities_master()` unless the
   implementation can prove the constraints were not already materialized.

## TODO checklist

### A. Distributor: materialize periodic constraints as source rows

- [ ] Add a helper in `src/mesh/cl_Mesh_Distributor.*pp` that appends
  periodic-derived rows to each proc's `proto::TMatrixData` after
  `select_entities(p)` has finished and before `send_t_matrices()`.
  Current location to integrate: the per-proc loop in
  `Distributor::distribute()` around `select_entities(p)`,
  `populate_periodicity_data(p)`, and `populate_t_matrices(p)`.
- [ ] Do not rely on `Basis::sources()` already being set on the root mesh for
  periodic constraints. The distributor should derive these rows from
  `mMesh->periodicity()->master/slave_nodes/edges/faces()` plus the selected
  bitsets for proc `p`.
- [ ] Make row emission per-proc and local: only emit a periodic row when the
  target basis and all source basis entities exist in that proc's selected
  entity set. If a source is absent, that is an aura-closure bug and should be
  a hard error.
- [ ] Decide the deterministic representative at basis-row level. The current
  DOF-level `entangle()` uses smaller `Dof::id()` when neither side is hanging;
  the basis-level materializer needs an equivalent stable rule, or a verified
  proof that basis ID ordering produces identical DOF representatives for all
  DOF types on the paired bases.
- [ ] Preserve existing row composition semantics. If either periodic side
  already has sources, the materializer must unfold/merge rows like
  `DofData::entangle()` does; it must not simply overwrite a row with
  `target -> source`.
- [ ] Avoid duplicate target rows in `proto::TMatrixData`. If a target already
  has a thin-shell/hanging row from `populate_t_matrices()`, merge the periodic
  row before serializing instead of appending a second target entry that
  `ProtoMesh::create_t_matrices()` would process as a later overwrite.
- [ ] Fix the `mNumberOfAllTMatrices` gate. Today
  `Distributor::count_entities()` sets it from
  `mMesh->number_of_hanging_basis()`, and `send_t_matrices()` returns early
  when it is zero. Periodic-derived rows must be counted/sent even on meshes
  with no pre-existing hanging basis.
- [ ] Fix `populate_t_matrices()` sizing if needed. It currently sizes
  `mTargetIDs` from `mNumberOfAllTMatrices`; per-proc augmented rows may exceed
  that count unless the global count is updated or the function is changed to
  size from the actual per-proc row count.

### B. ProtoMesh: separate persistence periodicity from solve constraints

- [ ] Keep `proto::PeriodicityData` for save/load and full-mesh reconstruction.
  Do not delete it just because MPI assembly will use `TMatrixData`.
- [ ] Split the `ProtoMesh` behavior by context: persisted mesh load may call
  `create_periodicities()`, but MPI worker mesh construction should not require
  a worker-local geometric `Periodicity` object for assembly.
- [ ] Audit `Distributor::receive_periodicity_data()` and
  `mProtoMesh->create_periodicities()` in the worker path. If workers no longer
  need periodic geometry, either skip these calls during MPI distribution or
  create a non-updating data-only periodicity mode.
- [ ] Ensure a `Periodicity` reconstructed from `PeriodicityData` cannot trigger
  geometric `PeriodicityFactory::update()` during `Mesh::finalize()` on a
  partial/distributed mesh. Local worker meshes do not necessarily contain
  complete periodic planes.
- [ ] Keep `ProtoMesh::create_t_matrices()` as the worker solve entry point for
  periodic-derived constraints. Add duplicate-target guards if the distributor
  cannot guarantee uniqueness.

### C. Dof manager: remove or guard the periodic solve path

- [ ] Decide the serial path. Either keep
  `DofData::create_dofwise_periodicities_master()` for serial/root-only runs,
  or move periodic row materialization earlier so serial and MPI both consume
  the same basis-source representation.
- [ ] If periodic constraints are delivered through `TMatrixData`, guard
  `create_dofwise_periodicities_master()` so workers do not double-apply
  periodic constraints from both `Basis::sources()` and `mMesh->periodicity()`.
- [ ] Consider extracting the row-unfold/merge policy from
  `DofData::entangle()` into a reusable helper, or document why the
  distributor's basis-level materializer is algebraically equivalent.
- [ ] Keep the existing hanging-DOF MPI exchange. It already resolves source
  DOFs by global DOF ID and can create missing ghost-owned source DOFs, provided
  the source mesh basis exists locally.
- [ ] Add assertions in `collect_hanging_dofs()` or the periodic materializer
  that catch exactly-one-side-hanging periodic pairs unless the topology policy
  is deliberately relaxed.

### D. Thin-shell periodic facets

- [ ] Add generated thin-shell boundary/interface facets to sidesets that the
  periodicity kernel can discover. `PeriodicityFactory::select_sidesets()` scans
  mesh sidesets on the periodic plane; facets that exist only as thin-shell
  internal/generated geometry but are not present in a suitable sideset are
  invisible to periodic matching.
- [ ] Audit `CutFactory::duplicate_and_relink_facets()` and
  `restore_thin_shell_sidesets()`. During cut creation, thin-shell facets are
  duplicated and temporarily moved; after restore, the periodic update must see
  the final facets that should participate in periodic pairing.
- [ ] Audit `ThinShellFactory::create()` and `create_ghost_facets()`. The
  generated shell sideset (`DomainType::GeometryOnly`) and ghost sideset
  (`DomainType::Ghost`) may not cover the side facets created where a thin shell
  intersects a periodic plane. Those missing facets need a hidden geometry-only
  sideset or equivalent discoverable container.
- [ ] Decide ownership/domain type for these extra periodic thin-shell facets.
  They should be visible to `PeriodicityFactory` but must not become physical
  boundary conditions or ordinary Maxwell sidesets.
- [ ] Add a validation check in periodic matching that reports which source or
  target plane facet count came from regular boundary sidesets versus generated
  thin-shell sidesets. This will make missing thin-shell periodic facets obvious
  instead of surfacing only as a generic "number of facets does not match."

### E. Tests / verification

- [ ] Serial periodic case with no thin shell: periodic-derived rows produce the
  same condensed system as the existing DOF-level periodic path.
- [ ] Two-rank MPI case with no other hanging basis: verifies that
  periodic-derived `TMatrixData` is sent even when `number_of_hanging_basis()`
  is initially zero.
- [ ] MPI case where the periodic counterpart is on the opposite boundary and
  owned by another rank: verifies aura closure plus hanging exchange.
- [ ] Thin-shell periodic case: verifies generated thin-shell facets are seen by
  `PeriodicityFactory` and periodic layer/duplicate nodes produce source rows.
- [ ] Save/load round trip: verifies `PeriodicityData` persists semantic
  periodicity without duplicating or double-applying solve-time T-matrix rows.

## Confidence

- "Periodic equality is a hanging/entanglement constraint": **high.**
- "Workers need source rows, not geometric periodicity": **high.**
- "Existing hanging exchange can handle ghost-owned source DOFs if the mesh basis exists":
  **high** from reading `cl_FEM_DofMgr_DofData.cpp`.
- "Best implementation is per-proc TMatrixData augmentation after aura selection":
  **medium-high**. It matches the current distributor shape, but the basis-level row
  composition helper still needs careful implementation.

## Not in scope here (tracked elsewhere)

- Updating `todo/periodic_bc_fix_plan.md`: its older "frozen `mHangingDOFs`" diagnosis is
  stale for the active branch and should be replaced by this distribution/source-row plan.
- Post-cohomology duplicate/thin-shell periodic nodes getting their `periodic()` links set
  before distribution (Codex finding #3) — same plan; the duplicates must land in the
  pairing before aura selection and source-row materialization.
- Thin-shell periodic facet creation (next implementation step).
