# Periodicity save/load in BfmFile — audit + completion

**Date:** 2026-06-27
**Purpose:** Audit the newly-wired periodic-data save/load in the BFM mesh-file stack, drive the
fixes to a correct round-trip, and reconcile the refactor plan (R7).
**Module:** `src/mesh` (`cl_Mesh_BfmFile`, `cl_Mesh_PeriodicityFactory`, `cl_ProtoMesh`), `src/fem/maxwell`
**AIs:** Claude (primary), Grok + Codex (independent auditors). Three voices, consensus on the CRITICALs.

## What this was

Christian wired `save/load_periodicity_data` + per-entity helpers (planes, nodes, edges, faces,
facets) into `BfmFile`, plus `ProtoMesh::create_periodicitiy` → `PeriodicityFactory::from_proto`.
Audited with Grok and Codex (`tmp/ai_exchange/periodic_bfm.md`), fixed, then completed the reload
contract (R7).

## Findings (all three auditors)

Three **CRITICAL** defects in the first cut — all confirmed independently by Grok and Codex:

- **C1 — wrong loader dispatched** (`cl_Mesh_BfmFile.cpp:1669`): `load_hanging_faces()` instead of
  `load_periodic_faces()`. Because `group_exists`/`dataset_exists` are both `H5Lexists`, the guard
  passes on the periodic `"faces"` *dataset*, then `select_group` does `H5Gopen2` on a dataset → HDF5
  failure when periodic faces exist. **Fixed.**
- **C2 — plane matrix transposed** (`:1681`, `:1701`): `Matrix(3,2)` indexed `(0,k)/(1,k)` for `k<3`
  → column OOB at `k=2`. Plane is master/slave × 3 points ⇒ must be `Matrix(2,3)`. **Fixed** (save
  side `2×3`; load resizes via `load_matrix_from_file`).
- **C3 — unsized plane Cells** (`:1705-1711`): wrote `tMaster(k)/tSlave(k)` into default-empty
  `Cell`s with no `set_size(3)`. **Fixed.**

Plus, after the CRITICALs were fixed:

- **Codex HIGH (the substantive one)** — `from_proto` restored the pair *lists* but never set the
  per-entity `periodic()` back-pointers or the entity-type flags, while `set_entity_dependencies`
  hard-requires `A->periodic() != nullptr` (`cl_Mesh_Periodicity.cpp:88-89`) and the Maxwell trigger
  calls `update()` unconditionally (which throws on a restored periodicity, guard `:369-370`).
- **Grok HIGH** — `create_periodicitiy` gated on node counts, but planes are saved unconditionally;
  a planes-only file would load then silently discard. Resolved by gating on the **planes**.

## Resolution (R7 complete)

- **Pairs stored by authoritative partner, not list position.** `save_periodic_*` writes
  `(entity->id(), entity->periodic()->id())` per column, so loaded `mMasterX(k)↔mSlaveX(k)` are true
  partners — sidesteps the "independently-ordered lists" hazard (`periodicity.md:88`).
  `mMaster/SlaveSideSets` are therefore **not** serialized (unused by `from_proto`).
- **Reload via O1 option (i):** `create_periodicitiy(true)` → `from_proto(..., aCrosslinkEntities=true)`
  → `crosslink(Periodicity*)` (`PeriodicityFactory.cpp:1426`) sets bidirectional `set_periodic()` and
  flags NODE/EDGE/FACE/FACET (each gated on `master_*.size()>0`). No geometric re-match.
- **Periodic DOF constraints** (slave `add_source(master)`) are *not* re-derived — they ride back via
  the hanging-entity data (`save_hanging_nodes` saves every `is_hanging()` entity).
- **Maxwell read-side gate done** (`cl_MaxwellFactory.cpp:481-489`): `update()` +
  `set_entity_dependencies()` run only when `! node_pairs_restored()`.

## Invariant to honor when the write side is wired (R10)

The enriched `.bfm` must be written **after** `set_entity_dependencies()` runs — otherwise the
periodic-slave source containers don't yet exist, and a restored mesh would have `periodic()` links
but no enforced periodic constraints (silent wrong solve). The Maxwell `.bfm` write is still
commented out, so this is not yet exercised.

## Also landed this session

Abstract + orphan node save/load (`cl_Mesh_BfmFile.cpp:147-211`), ID-keyed: the node *data* rides in
the main `ids`/`coords` table; separate `abstract`/`orphaned` ID lists restore category membership on
load via `mProto->node(id)`. Plan row 8 marked done.

## Plan delta

`todo/meshfile_refactor_plan.md`: R7 → done (with the build notes above); O1 → resolved (option i);
O2 → resolved (no geometric `update()` needed on restore); inventory rows 8 + 15 → done; "still
missing" list trimmed (periodicity + abstract/orphan removed). Remaining: R6 (fields/globals/time),
R8 (thin shells), R9 (block↔material), R10–R11 (write-side + reload orchestration), R12 (restart test).
