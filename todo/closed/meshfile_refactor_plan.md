# BFM Mesh-File Save/Load Refactor Plan

**Date:** 2026-06-22
**Purpose:** Plan to repair and refactor the `.bfm` (HDF5) mesh save/load around the
`ProtoMesh`/`Protoshell` reconstruction engine so a Maxwell run can stop, change parameters,
and relaunch from the last stored timestep — loading a *fully enriched* mesh and skipping the
cohomology / cut / thin-shell factories on reload.
**Module:** `src/mesh` (I/O) + `src/fem/maxwell` (reload orchestration)
**AIs involved:** Claude (exploration + plan), Codex (audit, high-confidence), Grok (third-voice audit)
**Status:** ✅ **COMPLETE (2026-07-01).** The ProtoMesh-backed `BfmFile` I/O backend (§2) is built,
audited (tri-AI, multiple rounds), and wired end-to-end. The R12 restart test passes on the CORC
reproducer (periodic + thin-shell + cohomology cuts) in **serial and on 2/4/8 MPI procs** after the
last two blockers landed — D24 (aliased thin-shell facets serialized twice) and D25 (`Topology::run()`
on the enriched mesh diverging from fresh). All D1–D25 defects and R2–R12 steps are resolved; the
format is documented in `src/mesh/doc/bfm_file_format.md`. **Residual follow-ups (minor, additive,
do not block the mesh save/load restart):** coupled-circuit dynamic-state restart (R6b, needed for a
true CORC circuit restart — see `todo/restart_circuit_state.md`) and the `/meta/format_version`
attribute (§6.4). *(R6a global variables is done — `Mesh::save_globals`/`load_globals`.)* Audit threads (swept): `tmp/ai_exchange/meshfile_refactor.md`,
`bfmfile_facet_audit.md`, `bfmfile_edge_audit.md`, `bfmfile_face_audit.md`.

> **Scope guards (from the task brief):**
> - Mesh load/read is **serial, master-proc only**. Parallel mesh distribution is **out of scope**.
> - **Backward compatibility is explicitly dropped** — we choose a clean format over legacy support.
> - Entities in scope: nodes, elements, and (when present) edges, faces, facets, control points,
>   plus the enrichment state (duplicates, cuts, thin shells, ghosts, periodicity, curves/terminals,
>   hanging/T-matrix data) and physical fields at the last timestep.

---

## 1. Summary of Current Save/Load Behaviour and How It Fails

> **HISTORICAL (pre-refactor baseline) — reconciled 2026-06-28.** §1.1-1.3 describe the original
> legacy `HDF5Writer`/`HDF5Reader` world, which **no longer exists** (those files are deleted). They
> are kept only for the rationale that motivated the refactor. What still holds today: the enriched
> write is *not yet wired* into `MaxwellFactory`/`Mesh::save` (R10), and periodicity/thin-shell/field
> restore are unimplemented (R6-R8). For the **current** implementation state and live defects, read
> §4.0 — not this section. Citations to `cl_Mesh_HDF5{Writer,Reader}.cpp` below are dead links.

### 1.1 The write path is disabled
The enriched-mesh BFM write is **commented out**. In `MaxwellFactory::create_magnetic_kernel()`
the block that would dump the post-factory mesh is replaced by a `#warning`:

```
src/fem/maxwell/cl_MaxwellFactory.cpp:469-478
    std::cout << "#warning: BFM mesh saving disabled" << std::endl ;
    /*mesh::HDF5Writer tWriter( mMeshPath, mMesh, mSaveElementConnectivitiesToBfm );
      this->save_material_map( tWriter.file() );
      this->save_thin_shell_facets( tWriter.file() );
      tWriter.file().close();*/
```

As a result, **no enriched `.bfm` is produced** today by a Maxwell run. The plain `Mesh::save("x.bfm")`
path (`cl_Mesh.cpp:353-357` → `mesh::HDF5Writer`) still works, but writes only base mesh
topology — see §1.3 for what it omits.

**Provenance:** the write was disabled in commit `6c52293` ("more work on thin shells", 2025-08-26),
a large refactor that also removed `src/fem/maxwell_old/`. No devlog or todo explains it, so it was
most likely defensive ("turn it off while restructuring") rather than a response to a known crash.
The matching `HDF5Reader` was left wired up and still consumes existing `.bfm` files, so the legacy
reader's bugs (Appendix A) are live today even though the writer is dark.

### 1.2 The reload path half-exists and is internally inconsistent
`MaxwellFactory::read_mesh()` (`cl_MaxwellFactory.cpp:227-313`) already *prefers* a `.bfm`
sidecar: if `mesh.bfm` exists and its checksum matches the source `.msh`, it loads the `.bfm`
and sets `mComputeCohomologies = false` (`:291-296`) — i.e. it intends to skip the cut/cohomology
factories. On reload `create_magnetic_kernel()` calls `load_domain_types_and_material_map()`
(`:433`, defined `:2124`) to restore block/sideset `DomainType` and the block→material map.

But because the write side is disabled, the `.bfm` that the reload path expects **never contains
the enriched state** (thin shells, ghosts, periodicity, cuts beyond node duplicates, fields). The
two halves were never completed together. Concretely:

- `save_thin_shell_facets()` (`:2173-2207`) writes a `"tapes"` HDF5 group that **nothing ever reads**.
- `save_material_map()` (`:2078-2119`) is only called from the disabled block, so the materials
  that `load_domain_types_and_material_map()` reads (`:2135`) are never written by the live path.
- Periodicity, `ThinShell` objects, and physical fields are **never serialized at all**.

### 1.3 Failure modes, classified
For a Maxwell problem that uses any enrichment, attempting save→reload fails as follows:

| Failure | Mechanism | Evidence |
|---|---|---|
| **Silent data loss — periodicity** | Reader never builds a `Periodicity` object → `has_periodicity()==false` → the rebuild block `create_magnetic_kernel:481-485` is skipped entirely. Periodic BCs silently vanish. | `cl_Mesh.hpp:1859-1861`; `cl_MaxwellFactory.cpp:481` |
| **Hard error — periodicity (if restore attempted)** | If a future restore fills `PeriodicityData` via `from_proto`, that sets `mNodePairsRestored=true` and never sets `periodic()` pointers (`PeriodicityFactory.cpp:320`). A subsequent `update()` with no node-pair backup hits `BELFEM_ERROR` at `:369-370`. | verified `PeriodicityFactory.cpp:320, 359-370` |
| **Silent data loss — thin shells** | `HDF5Reader` recreates `SideSet`s but never `ThinShell` objects (`mThinShells` stays empty). Kernel code that iterates `thin_shells()` finds nothing → no thin-shell physics on reload. | reader has no thin-shell path; `cl_Mesh.hpp:90, 1637` |
| **Silent data loss — fields/time** | Writer never serializes fields, global variables, time step/stamp → a "restart" starts from zero state. | `HDF5Writer` ctor `:26-57` lists no field/global/time writes |
| **Exception — element type** | If a reloaded mesh contains a 20-node 3D element (TET20/HEX20) or HEX8TS, type reconstruction throws. See §3 Row 3. | `meshtools.cpp:1181-1184, 1031-1058` |
| **Latent OOB — 2D faces** | `read_faces_2d()` writes `tFaces(f)` without `set_size()` first (3D path does set_size at `:855`). Any 2D mesh saved with a `faces` group corrupts memory on reload. | `cl_Mesh_HDF5Reader.cpp:811-836` vs `:855` |

**Bottom line:** the current state is "write disabled + reader that reconstructs base topology only".
The refactor must (a) re-enable an *enriched* write, (b) extend it to all enrichment + field state,
and (c) make the reload deterministic and ID-keyed via `ProtoMesh`.

---

## 2. Architecture: why ProtoMesh is the right spine

The reconstruction engine already exists and is battle-tested by the MPI distributor.
`ProtoMesh` (`cl_ProtoMesh.{hpp,cpp}`) holds **ID-keyed** proto structs and rebuilds a live mesh:

- `proto::NodeData, ElementData, ElementExtra, EdgeData, FaceData, FacetData, FacetExtra,
  ControlPointData, GroupData, ThinShellData, PeriodicityData, TMatrixData` (`cl_ProtoMesh.hpp:36-154`).
- Builders: `create_nodes/elements/edges/faces/facets/vertices/control_points/blocks/sidesets/
  thinshells/periodicitiy/t_matrices` (`cl_ProtoMesh.hpp:258-299`).
- The MPI `Distributor` packs the live mesh into these structs (e.g. `share_thinshell_data`
  `cl_Mesh_Distributor.cpp:1649-1668`; T-matrix receive `:2276-2281`) and `ProtoMesh` rebuilds.

**The refactor adds a second backend for the same proto structs: HDF5 instead of MPI.**
Writer: live `Mesh` → proto structs → HDF5 groups. Reader: HDF5 groups → proto structs →
`ProtoMesh` builders → live `Mesh`. This reuses the proven, ID-keyed reconstruction and replaces
the current index-keyed `HDF5Reader` ad-hoc path. (The legacy `HDF5Writer`/`HDF5Reader` can be
kept for non-Maxwell base meshes, or retired — see Open Question O7.)

---

## 3. Gap Table

Legend for **Class**:
- **(a)** deterministically rebuildable on load (states *from what*),
- **(b)** ambiguous / open question,
- **(c)** must be saved explicitly.

> **HISTORICAL (reconciled 2026-06-28).** The **"Saved today?"** column and every
> `HDF5Writer.cpp`/`HDF5Reader.cpp`/`HDF5_Tools.cpp` citation below refer to the **deleted** legacy
> path and are dead. The table's **value is its (a)/(b)/(c) classification of what must be saved** —
> that still holds and drove the BfmFile design. For the *current* per-entity status (done / missing),
> read §4.0's "Implemented" + "Still unimplemented" lists, not the "Saved today?" column.

"Saved today?" refers to the (now-deleted) legacy `HDF5Writer`. Citations are historical.

| # | State | Needed for | Saved today? (legacy — dead) | Class | Code citation / rationale (historical) |
|---|---|---|---|---|---|
| 1 | Node coords + IDs | base geometry | YES | (a)→keep | `HDF5Writer.cpp:71-122`; reader `:152-206` |
| 2 | Element ID + topology + block_id | base | YES (topology by node **index**) | (c) save, **ID-key it** | writer `:204-281`; reader `:314-414` |
| 3 | **Element TYPE** | element creation | NO (rebuilt from dim+nnodes+TS-flag) | **(c) must-save** | `element_type_from_numnodes` **errors** on 20-node 3D (TET20/HEX20 collision) `meshtools.cpp:1181-1184`; TS switch only 4/9/6/18 nodes `:1031-1058` so **HEX8TS unrecoverable** (`cl_Element_Factory.cpp:185-187`). Not a deterministic function of (dim,nnodes). |
| 4 | Block/sideset ID + label + domain_type | grouping, physics dispatch | YES | (a)→keep | writer `:844-891`; reader `:537-654` |
| 5 | Facets: master/slave + indices + orient_on_slave | interfaces, BCs | YES (master/slave by **index**) | (c) save, **ID-key it** | writer `:284-403` (`:316-317` store element index); reader `:418-532` |
| 6 | Edges + Faces (2D/3D) topology | H(curl)/H(div) DOFs | YES | (a) relink on load | writer `:648-841`; reader `:688-918`. **Bug:** `read_faces_2d` missing `set_size` `:811-836`. |
| 7 | Node duplicates (cut/TS): dup↔original | cohomology cuts, φ⁺≠φ⁻ | YES (by **index**) | (c) save, **ID-key it** | writer `:158-200`; reader `:210-266` |
| 8 | Abstract + orphan nodes | cut boundary conditions | YES (by **index**) | (c) save, **ID-key it** | ✅ **DONE in BfmFile** (ID-keyed): `save/load_node_data` `cl_Mesh_BfmFile.cpp:147-211` (separate `abstract`/`orphaned` ID lists, re-tagged on load via `mProto->node(id)`). |
| 9 | Hanging / static-condensation T-matrices | hanging edges, enrichment | PARTIAL — `"hanging"` group covers node/edge/facet/element only; **by index** | **(c) must-save, use proto `TMatrixData`** | writer `:435-645` (index-keyed `:597-638`); proto adds **face + control-point** targets, ID-keyed (`cl_ProtoMesh.hpp:53-61`, builder `:830-...`). Not equivalent. |
| 10 | Curves + segments + terminal sideset_a/b | terminals | YES | (a)→keep | writer `:894-1002`; reader `:1091-1250` |
| 11 | Vertices | point entities | YES | (a)→keep | writer `:405-429`; reader `:658-684` |
| 12 | **Physical fields @ last step** (node/element/edge/face) | **restart** | **NO** | **(c) must-save** | `HDF5Writer` ctor `:26-57` writes no fields. Field shape: `cl_Mesh_Field.hpp:34-43` (`mLabel,mID,mEntityType,mData,mWriteToFile`). Currently only Exodus serializes fields (`cl_Mesh_ExodusWriter.cpp:535-632`). |
| 13 | **Global variables** (label,id,value) | restart, total current/energy | **NO** | **(c) must-save** | `cl_Mesh.hpp:92`; `cl_Mesh_GlobalVariable.hpp`. Only Exodus writes them (`ExodusWriter.cpp:654-711`). |
| 14 | **Time step + time stamp** | restart cursor | **NO** | **(c) must-save** | `cl_Mesh.hpp:105-106` (`mTimeStamp,mTimeStep`); serialized only to Exodus (`ExodusWriter.cpp:637-651`). |
| 15 | **Periodicity** (planes + matched node/edge/face/facet ID pairs + master/slave sidesets) | periodic BC | **NO** | **(c) must-save + dedicated reload step** | proto `PeriodicityData` `cl_ProtoMesh.hpp:140-153`; `to_proto` `PeriodicityFactory.cpp:136-217`. ✅ **DONE — see R7** (BfmFile `:1639-1906`; partner-ID pairs + `crosslink` on load; sidesets intentionally not serialized). |
| 16 | **ThinShell objects** (sideset, ghost sideset, blocks, thicknesses, materials) | skip TS factory on reload | **NO** (`"tapes"` written, never read) | **(c) must-save, use proto `ThinShellData`** | proto `ThinShellData` `cl_ProtoMesh.hpp:114-121`; distributor pack `cl_Mesh_Distributor.cpp:1649-1668`; builder `cl_ProtoMesh.cpp:745-817`; dead `save_thin_shell_facets` `cl_MaxwellFactory.cpp:2173`. |
| 17 | Ghost sidesets (facets, `DomainType::Ghost`) | ghost stabilization | PARTIAL — facets serialize via generic sideset path | (b)→(c) for the **link** | `ThinShellFactory.cpp:287-303` creates them; writer iterates all sidesets `:290-308` (no `is_hidden` filter). But the **ThinShell→ghost link** is only restored via row 16's `mGhostSideSetID` (`cl_ProtoMesh.cpp:785`). So facets: rebuildable; link: must-save (row 16). |
| 18 | Block↔material assignment | material reassignment on reload | NO in live path (writer disabled) | **(c) must-save** | `save_material_map` `:2078` (disabled) / `load_domain_types_and_material_map` `:2124`. |
| 19 | Node/element owners (`proc_t`) | parallel distribution | NO | (a) / **out of scope** | proto `NodeData.mOwners` exists; serial reload sets owner=master. Not needed in scope. |
| 20 | Interpolation TYPE (Lagrange/Bernstein/Hermite) | shape functions | NO (assumed Lagrange) | **(b) config-owned** | non-Lagrange exists `Mesh_Enums.hpp:94-106`; IWG defaults Lagrange `cl_IWG.hpp:279`; Maxwell does not set it (`cl_MaxwellFactory.cpp:410-420`). Fine for Maxwell; **not** deterministic mesh state — log as open. |
| 21 | Curved-element flag | integration order | NO | (a) recomputed | reader calls `flag_curved_elements()` `:116`; def `cl_Mesh.cpp:1793-1806`. |
| 22 | Edge directions / orientation | H(curl) sign | NO | (a) for **local** signs; (b) under periodic | recomputed in `finalize_edges()` `cl_Mesh.cpp:981-997` from endpoint `original()->id()` `:1821-1854`. Periodic sign correctness depends on row 15, not this row. |
| 23 | Checksum | cache validity vs source `.msh` | YES | (a)→keep | writer `:61-67`; reader `:140-148`; used `cl_MaxwellFactory.cpp:265-273` |

### 3.1 Cross-cutting finding: index vs ID keying (central risk)
Rows 2, 5, 7, 8, 9 round-trip **today only because serial reload preserves array order** — the
writer stores `entity->index()` and the reader relies on load order matching save order
(`HDF5Writer.cpp:183-196, 249-255, 315-326, 597-638`; reader `crate_nodes` + `update_node_indices`
`HDF5Reader.cpp:177-204`). The moment any reordering (RCM, a different build, a partial mesh) is
introduced, duplicates/facets/abstract/hanging all break **together and silently**. The proto
structs are already **ID-keyed**; the refactor must adopt ID keying throughout. This is the single
most important correctness change.

---

## 4. Ordered Refactor Steps

Dependencies noted as `(after: …)`. Each step is independently testable. Steps R0–R5 build the
new I/O backend; R6–R9 wire enrichment + restart; R10–R12 integrate and verify.

### 4.0 Implementation Progress (updated 2026-06-28 — re-audit)

The new ProtoMesh-backed backend (§2) now lives in `src/mesh/cl_Mesh_BfmFile.{hpp,cpp}`,
driving `ProtoMesh` `populate_*` / `create_*` / `reset_*`. The legacy
`cl_Mesh_HDF5Reader`/`HDF5Writer` **have been retired** (O7 → retire; files deleted), so
Appendix A is gone — its vlen-pattern cautions now apply to the `hdf5::Dataset` helper.

**Implemented and audited (Claude + Codex + Grok, multiple rounds + a holistic double-check):**

- [x] **R2/R3/R4 (base topology + groups)** — save/load of meta, blocks, sidesets (`GroupData`),
  nodes, elements. ID-keyed; element **type** is carried **per block** (`blocks/types`, each block
  single-type) and reconstructed on load via block counts — this resolves the row-3 (dim,nnodes)
  ambiguity (it is *not* stored per element in `ElementData::mTypes`). Group tag bits
  (`mHidden`/`mHasEdges`/`mHasFaces`) decode symmetrically for **both** blocks and sidesets.
- [x] **Facets** — `save_facet_data`/`load_facet_data` + `create_facets` (no-topology path:
  re-derive nodes from master). hvl encoding master/slave/index/orientation; assert-at-save
  that every facet has a master.
- [x] **Edges** — `save_edge_data`/`load_edge_data` + `create_edges`; vlen-per-edge ↔
  length-prefixed proto ↔ builder verified.
- [x] **Faces** — `save_face_data`/`load_face_data` + `create_faces`; mirrors the facet scheme.
- [x] **Element→edge/face incidence (D3)** — `reconstruct_edge_connectivity` (node-key attach,
  IDs preserved) + `reconstruct_face_connectivity` (from `FaceData` master/slave + index).
  **Nothing from `ElementExtra` is serialized** — see the revised Appendix B decision below.
- [x] **Post-load finalization (D4)** — single `mMesh->finalize()` at the end of `load()` does
  directions, inverses, and all id→entity maps; `Connectivity::Compute` is set (ctor default),
  `Connectivity::FacetToElement` is pre-set so facet master/slave is not clobbered.
- [x] **Distributor element/facet table sort** — `populate_element_data`/`populate_facet_data`
  sort each per-target comm table by **block/sideset index** (not id value), letting
  `create_blocks`/`create_sidesets` use a single linear sweep; audited safe for the live MPI path.
- [x] **`suint`→`uchar` narrowing** of proto index/type/orientation fields (safe: `ElementType`
  max = 127 < 256).
- [x] Per-round fixes applied (facet F1–F7/P1–P2/G2–G6, edge E1–E4, face FA1–FA3, meta typo,
  reconstruct edge sanity-check `f`→`e`/gate).

**Open defects found in the 2026-06-26 holistic double-check** (must clear before R5 round-trip):

- [x] **D1 — node-duplicate serialization corruption (HIGH).** ~~`populate_node_data` writes the
  duplicate records into `tIDs` instead of `tDup`/`mDuplicateData`
  (`cl_ProtoMesh.cpp:167,168,172`), corrupting node IDs and leaving duplicate metadata as a bare
  header.~~ **Fixed 2026-06-26** — writes routed to `tDup`; bookkeeping (`tMemCount`) and the
  `create_nodes` restore loop (`cl_ProtoMesh.cpp:261-289`) are now symmetric, so node duplicates
  round-trip. (Was dead on the live path anyway: `BfmFile::save` is unwired — `Mesh::save("*.bfm")`
  throws "not implemented", MaxwellFactory write commented out — and the MPI distributor uses its
  own correct `Distributor::populate_node_data`. Fires only for cut/thin-shell meshes once wired.)
- [x] **D2 — `create_sidesets()` asymmetry (HIGH).** ~~Currently in `save_facet_data:622` →
  double-inserts the live mesh's already-owned facets into fresh `SideSet`s (double-free on
  `~Mesh`); `load_facet_data:816` never calls it → loaded mesh has empty `sidesets()` + leaked
  facets.~~ **Fixed 2026-06-26** — moved to `load_facet_data` after `create_facets`
  (`cl_Mesh_BfmFile.cpp:813,815`); removed from `save_facet_data`. (Grok + earlier Claude, verified.)
  Note: `create_blocks`/`create_sidesets` stay deliberately **defensive** (count + ID-index map)
  — they are shared with the MPI distributor, where element/facet order is not guaranteed; do NOT
  specialize them for the "HDF5 is already ordered" case (cold path, would fork shared logic).
- [x] **D3 — element→edge/face reverse-incidence (HIGH).** ~~`create_edges`/`create_faces` build
  entities but don't relink `Element::edge(k)`/`face(k)` or directions.~~ **Fixed 2026-06-27 — by
  reconstruction, storing nothing.** `reconstruct_edge_connectivity` builds element→edge incidence
  by node-key (`min*N+max` of corner-node indices) attaching the **stored** edges (IDs preserved);
  `reconstruct_face_connectivity` builds element→face from `FaceData`'s master/slave + `index_on_*`.
  Directions are **not** stored — they fall out of `compute_edge_directions()` in finalize (D4).
  This revises Appendix B's "store edge/face links" conclusion (that assumed re-derivation meant
  re-*creating* edges with fresh IDs; the node-key attach keeps the stored IDs). Gated on block
  `mHasEdges`/`mHasFaces`, which is why the symmetric group-tag decode was a prerequisite.
- [x] **D4 — Mesh maps + finalize (HIGH).** ~~`BfmFile::load()` rebuilds no maps/finalize.~~
  **Fixed 2026-06-27** — `load()` ends with a single `mMesh->finalize()`, which (with
  `Connectivity::Compute` set by the `Mesh(uint)` ctor default) runs `create_maps`,
  `finalize_edges` (inverse + `compute_edge_directions` + `create_edge_map`), and `finalize_faces`
  (`create_face_map` + face↔edge). `Connectivity::FacetToElement` is pre-set in `load_facet_data`
  so `connect_facets_to_elements` is skipped and restored facet master/slave survives. (The two
  redundant explicit `finalize_edges`/`finalize_faces` calls were removed — `finalize()` already
  does both.)
- [x] **D5 — coords transpose OOB + orientation mismatch (HIGH).** ~~Two coupled defects on the
  coordinate round-trip (node + control-point coords; all four calls pass `transpose=true`):
  **(a) OOB** — `load_matrix_from_file(...,true)` did `set_size(tDims[0],tDims[1])` then
  `aMatrix(i,j)=tData[j][i]` over a `tDims[0]`-row buffer → flattened access `j*tDims[1]+i` ≫
  allocation (test n=5,d=2: max offset 21 vs 10-double buffer; round-trip BROKEN for every `n≫d`).
  **(b) Orientation** — `populate_node_data` filled `mCoords` **(n×d)** `(node,dim)` while
  `create_nodes` reads `(d×n)` `(dim,node)`; `transpose=true` on both ends didn't reconcile them.~~
  **Fixed 2026-06-27** — single convention: canonical in-memory layout is **(d×n)** `(dim,entity)`
  (`create_nodes`/`create_control_points` and both `populate_*` agree).
  - `populate_node_data` fill swapped to `(dim,node)` to match its `set_size(d,n)` (user).
  - `load_matrix_from_file` transpose branch now `set_size(tDims[1],tDims[0])` + `aMatrix(j,i)=tData[i][j]`
    (`hdf5_tools.hpp:738-786`) → exact inverse of the `save_matrix_to_file` transpose, no OOB.
  - Net round-trip: `save(transpose=true)` writes readable **(n×d)** disk (one entity/row);
    `load(transpose=true)` restores **(d×n)** in memory. *(Repro that drove this: `/tmp/d5_test2.cpp`.)*
- [x] **D6 — loader hardening — MOSTLY MOOT (downgraded 2026-06-27).** ~~No `dataset_exists`/extent
  cross-checks before `H5Dopen1`/`H5Dread`.~~ **Premise was wrong:** every scalar/vector/matrix loader
  guards with `BELFEM_ERROR( dataset_exists(...) )` (`hdf5_tools.hpp:168,339,497,704`), and every vlen
  dataset goes through the `Dataset` read ctor (existence + vlen type-integrity guard, both always-on,
  `cl_HDF5_Dataset.hpp:53-112`). **Only residue (optional, deferred):** release-mode cross-dataset
  extent consistency — loaders loop on the `ids` count and index sibling vlen rows whose per-row bound
  is only a debug `BELFEM_ASSERT`; a one-line `BELFEM_ERROR( tElements.size()==tNumFaces, … )` per
  loader would catch a truncated file in release. Low value for a self-produced restart format.
- [x] **D7 — error paths not RAII-safe (LOW, debug-only).** ~~`BELFEM_ERROR` can throw after HDF5
  handles / heap arrays are allocated, leaking them in debug builds — a throwing ctor never runs
  `~Dataset()`.~~ **Fixed 2026-06-27** — the read-mode ctor now closes every handle/buffer already
  opened *before* raising at each of the four validation throw-sites (existence, vlen type-integrity,
  extent, `H5Dread`), so no partial-construction leak (`cl_HDF5_Dataset.hpp`). No `try/catch` — release
  is `-fno-exceptions`; cleanup-before-throw is the only portable pattern. The write branch has no
  `BELFEM_ERROR` throw-after-allocate. Also explicitly `= delete`d the move ctor/assignment (were
  already implicitly suppressed by the user-declared dtor + deleted copy; now documented).

**Re-audit defects — 2026-06-28 (current tree):**

- [x] **D8 — CRITICAL — `save_element_data` row index never advances.** ~~`set_size(tCount,n)` with
  `tCount` never incremented.~~ **Fixed 2026-06-28** — `set_size( tCount++, n )` (`cl_Mesh_BfmFile.cpp:298`).
- [x] **D9 — CRITICAL — `save_facet_data` row index never advances.** ~~elements/indices loop +
  dormant `aSaveTopology` branch missing `++tCount`.~~ **Fixed 2026-06-27** — live loop advances
  (`tElements.set_size(tCount,…)` / `tIndices.set_size(tCount++,…)`) **and** the topology branch now
  `set_size( tCount++, n )` (`:418`) + `tData.save()` (`:428`).
- [x] **D10 — CRITICAL — `save_face_data` two defects.** ~~(a) `tCount` never advanced; (b)
  `close()` instead of `save()` → face elements/indices never reach disk.~~ **Fixed 2026-06-27** —
  all three cases use `set_size(tCount,N)` / `set_size(tCount++,M)` (one advance per face), and
  `:841-842` now call `tElements.save()` / `tIndices.save()`.
- [x] **D11 — FALSE POSITIVE (retracted 2026-06-27).** ~~`load_face_data` opens `elements`/`indices`
  without per-dataset guard.~~ Backwards: the read-mode `Dataset` ctor (`:861-862`) already runs
  `BELFEM_ERROR( dataset_exists(...) )` + a vlen type-integrity guard (`cl_HDF5_Dataset.hpp:53-112`),
  strictly *more* checking than `load_facet_data`'s manual `group_exists`. No defect.
- [x] **D12 — HIGH (data loss) — `save_control_point_data` never calls `tTopo.save()`** (Grok catch,
  Claude-confirmed). ~~The `"topology"` Dataset was filled but never `save()`d → incidence never
  written.~~ **Fixed 2026-06-28** — `tTopo.save()` added (`cl_Mesh_BfmFile.cpp:965`).

> **D8-D10 are the same row-index / close-vs-save class previously recorded fixed on 2026-06-27
> (Dataset-refactor pass); the current tree does not contain those fixes — likely lost in a later
> file restructure. Flag for the user before the fix pass.** By contrast `save_edge_data`
> (flat index `e`, `:679/693`), `save_control_point_data` (`:961`), `save_node_duplicate_data`
> (`:1057`), and all five `save_hanging_*` (`tIDs(tCount++)`) are correct.

> **Master-proc guard (Constraint B):** `BfmFile::save` (`:37`) and `load` (`:72`) have no
> `comm_rank()==master` guard — add at the `Mesh::save/load` call sites or atop these methods in a
> later pass; do not add now.

**MaxwellFactory ↔ BfmFile integration defects — 2026-06-29 (tri-AI: Claude + Grok + Codex).**
Audit of the save/load path as driven by `hphirun` (`MaxwellFactory` ctor → `read_mesh` →
`create_magnetic_kernel`). **Net: the `.bfm` round-trip cannot complete a load+solve today** — it does
not even compile (D13), and beyond that has double-execution, empty-ids, missing-return, a materials
crash, and double-creation of shells/curves. Work these top-down tomorrow.

- [x] **D13 — CRITICAL (compile blocker) — `Cell< string > tLabels( n,   )` was invalid C++.**
  ~~`save_curve_data` empty second ctor argument → TU did not compile.~~ **Fixed 2026-06-30** —
  `tLabels(n)` / `tTypes(n)` / `tClosed(n)` reserve + `push`; `tIDs(n,0)` + index-assign; `++c` advances
  (`cl_Mesh_BfmFile.cpp:2185-2249`). Compiles, internally consistent. (Audited Claude 2026-06-30.)
- [x] **D14 — CRITICAL — `BfmFile` ctor double-load/save + `checksum()` null-deref.** **Fixed 2026-06-30.**
  Ctor made passive (`:23-29` empty) → double-load/double-save gone (`read_mesh`/`Mesh::save` each call
  once); `checksum()` (`:38-50`) is now self-contained (`new HDF5(OPEN_RDONLY)` → read `meta/checksum` →
  `close()` → `delete`), no `mFile`/`mProto` dependency, leak-free. `get()` sets `mOwnMesh=false` (`:228`)
  so returning a stack `BfmFile`'s mesh is safe (no double-free). (Codex found; verified Claude 2026-06-30.)
- [x] **D15 — CRITICAL — `save_group_data` never stored group IDs.** ~~Loop pushed everything but
  `tIDs`.~~ **Fixed 2026-06-30** — `tIDs.push(tData.mID)` re-added (`cl_Mesh_BfmFile.cpp:270`).
  (Claude found, Codex-confirmed; verified fixed.)
- [x] **D16 — CRITICAL — `read_mesh` fell off a non-void function on first run.** **Fixed 2026-06-30** —
  `read_mesh` restructured: the source mesh is loaded + scaled **outside** the `if(exists)` block, so all
  four paths return (`.bfm` input `:243`; sidecar-match `:290`; no-sidecar/mismatch fallthrough `:302`).
  No fall-off. (all three; verified Claude 2026-06-30.)
  - [x] **D16b — cache the base checksum before enrichment. Fixed 2026-06-30** —
    `size_t tChecksum = aMesh->checksum();` moved above `if(exists)` (`cl_MaxwellFactory.cpp:276`), after
    `scale_mesh`, so `mHash` holds the **base** value before enrichment on every non-`.bfm` path (incl.
    first run). `save_meta_data` then stores the base checksum → the sidecar gate matches on later runs.
- [x] **D17 — CRITICAL — regular block materials unassigned on `.bfm` load → crash.** **Fixed 2026-06-30.**
  ~~`create_block_to_material_map` was only reachable inside `create_cuts`'s `mComputeCohomologies` block,
  and its fill loop was itself `mComputeCohomologies`-gated → empty map on load → `assign_materials`
  keyed lookup `Map::operator()` "Key not found" crash.~~ Two changes: (1) the internal
  `if(mComputeCohomologies)` around the fill loop **removed** — `create_block_to_material_map` now
  populates unconditionally from `collect_material_labels_from_domains` (input domains, load-safe); (2)
  `create_magnetic_kernel` now calls it on the load path via `else { create_block_to_material_map(); }`
  (`cl_MaxwellFactory.cpp:427-430`), mutually exclusive with the fresh `create_cuts()` call. Regular
  blocks now in the map on both paths → no crash. (all three; verified Claude 2026-06-30.)
  **Coupled to D18:** thin-shell *layer* block→material entries still come only from `create_thinshells:886`;
  when D18 makes the load path skip recreation, those must be hydrated from the loaded
  `ThinShell::materials()` or the crash returns for thin-shell meshes. D17+D18 land together for TS meshes.
- [x] **D18 — HIGH — `create_thinshells` not load-aware → double-create or materials-missing.**
  **Fixed 2026-06-30.** `create_thinshells` now early-returns when `mMesh->thin_shells().size() != 0`
  (shells came from the `.bfm`) — no `ThinShellFactory::create`, so no duplicate shells/blocks. Before
  returning, it **hydrates** `mMaterialBlockAssignment` from the loaded `ThinShell::materials()` /
  `blocks()` (fill-if-missing, `++b` unconditional, `cl_MaxwellFactory.cpp:882-901`), so the loaded
  shells register + label exactly like freshly-created ones (closes the D17 coupling for TS meshes).
  `comm_barrier` symmetry preserved (barrier before the early return matches the tail barrier). (Codex
  found; fix + hydration verified Claude 2026-06-30.)
  - [x] **D18b — `++b` in the create loop decoupled. Fixed 2026-06-30** — `cl_MaxwellFactory.cpp:928-937`
    now advances `b` for every layer block (assigns only when `! key_exists`), so the material index can
    no longer drift if a layer block is ever pre-populated.
- [x] **D19 — HIGH — curves duplicated on load.** **Fixed 2026-06-30** — ctor now guards the input
  `create_curves` with `if ( mMesh->curves().size() == 0 )` (`cl_MaxwellFactory.cpp:108`): skipped when
  `load_curve_data` already restored curves, still runs on a fresh mesh. **Bonus:** a parallel
  `if ( ! mMesh->has_periodicity() )` guard added around `create_periodic` (`:116`) — closes the same
  load-vs-fresh double-creation for periodicity. (Codex; verified Claude 2026-06-30.) Note: the factory's
  `mCurves` member is write-only (only `:370`), so it staying empty on load is harmless.
- [x] **D20 — FALSE POSITIVE (retracted 2026-06-30; Christian was right).** ~~Checksum reload gate is
  self-defeating — saved checksum is the enriched mesh's.~~ Both auditors missed the **caching**:
  `Mesh::checksum()` computes only `if (mHash.value()==0)` and returns the cached value otherwise
  (`cl_Mesh.cpp`); `read_mesh:270` calls it on the **base** mesh — after `scale_mesh` (`:267`), before
  enrichment — so `mHash` holds the base value, and the only `mHash.reset()` is *inside*
  `compute_checksum` (`:2375`), so enrichment never invalidates it. ⇒ `save_meta_data`'s
  `mMesh->checksum()` (`:150`) returns the cached **base** checksum; the sidecar gate is meaningful.
  **Caveats (tracked elsewhere):** the gate can't run until the D14 `checksum()` null-deref is fixed,
  and the D16 fresh-mesh branch must call `aMesh->checksum()` before enrichment to cache the base on
  first run.
- [x] **D21 — design DECIDED 2026-06-30 (Christian).** **The input file remains the source of truth.**
  The `.bfm` exists only to **cache the expensive enrichment** (cohomology/cuts, thin shells, hanging
  entities, periodicity) so a rerun can skip recomputing it; it does **not** own physics assignment.
  Concretely: regular block→material comes from `input.conf` on every run (D17:
  `create_block_to_material_map` ungated + called on the load path); thin-shell layer materials ride in
  the `.bfm` per-shell and are hydrated into `mMaterialBlockAssignment` on load (D18). `save_group_data`
  stores labels/domain-types only — intentionally no regular-block material strings. (Corrects an
  earlier stale Claude claim.)

**Equivalence audit — 2026-06-30 (tri-AI: Claude + Grok + Codex).** Question: does enrich→save→load
reproduce a bit-identical FE system? **Consensus verdict: YES for geometry/topology/orientation/enrichment,
with exactly ONE genuinely-divergent break (D22).** Both auditors independently confirmed and Claude
verified the call sites.

- [x] **D22 — HIGH (genuinely divergent) — cut φ DOFs dropped on `.bfm` load.** **Fixed 2026-06-30.**
  ~~`set_abstract_nodes`/`set_orphaned_nodes` ran only inside `create_cuts()` (`mComputeCohomologies`-gated,
  skipped on load), so the restored abstract/orphaned lists never reached the IWG → `count_node_dofs` saw
  empty lists → cut φ DOFs never created → smaller/wrong system on a cohomology restart.~~ The load `else`
  branch now installs them into the equation (`cl_MaxwellFactory.cpp:439/441`), mirroring the fresh path
  (`:766/768`). DOF IDs are entity-ID-keyed and the `.bfm` restores the lists in saved order, so the φ DOF
  numbering matches the fresh run. (Grok found; Codex + Claude confirmed; fix verified Claude 2026-06-30.)

- [x] **D23 — HIGH (parallel-only deadlock/crash) — type map never broadcast on `.bfm` load.** **Fixed 2026-06-30.**
  ~~The load `else` branch (`create_magnetic_kernel`) called `mTopology->run()` (internally rank-0-guarded —
  fills `mTypeMap` on rank 0 only, no collectives) but **not** the collective `mTopology->synchronize_maps()`
  that broadcasts the type map to the other ranks. The fresh path pairs them (`run()` :749 + `synchronize_maps()`
  :753/:786) but that lives inside `create_cuts()`, skipped on load. On `comm_size>1` the non-root ranks kept an
  empty `mTypeMap` → `groups(Conductor)` → `Key unknown not found in map` crash (`Map<DomainType,Vector<uint>*>::operator()`).~~
  The load `else` now calls `run()` **then** `synchronize_maps()` on all ranks, with the rank-0-local mesh work
  (`create_block_to_material_map`, `set_abstract/orphaned_nodes`) guarded by `if(mCommRank==0)` (`cl_MaxwellFactory.cpp:442-453`).
  **Post-fix load path verified parallel-safe by inspection through Kernel construction (`:521`):** `create_edges_and_faces_on_mesh`
  is `mCommRank==mMasterProc`-guarded (`cl_Mesh.cpp:1372`) so rank 0 skipping (edges exist from `.bfm`) while
  non-root enters is a no-op, not a deadlock; `create_thinshells` load branch early-returns after its own
  `comm_barrier()` (`:918`) matching non-root's `:970` barrier (1-per-rank balanced); `synch_material_map` is fully
  collective-balanced (both branches `comm_barrier`+2×`broadcast`) with rank-0 data populated by `create_block_to_material_map`
  (`:447`) + thin-shell loop (`:912`); `set_block_types_in_magnetic_equation` is local (no collectives); the rank-0
  periodicity `update()` (`:470`) is executed identically fresh vs load and proven collective-free by the working
  fresh 2-proc run. **Runtime confirmation still pending (R12):** this sandbox cannot launch `mpirun`, so the
  serial-path fix has not been executed on `comm_size>1` here — needs a `mpirun -n 2`/`-n 4` load run on the working
  MPI host. (Root-caused + fixed + inspection-verified Claude 2026-06-30; supersedes the "MPI `comm_size>1` not traced" residual below.)

**Parallel deployment crash — 2026-07-01 (Claude, root-caused; fix pending approval).**

- [ ] **D24 — CRITICAL (parallel `.bfm` load crash) — aliased thin-shell facets serialized twice →
  duplicate facet IDs → worker abort in `create_facet_extra`.** Repro: CORC, `mpirun` load path;
  SIGABRT `"container is already allocated"` at `cl_Vertex.cpp:216` via `ProtoMesh::create_facet_extra`
  (`cl_ProtoMesh.cpp:1234`) ← `Distributor::run` (`cl_Mesh_Distributor.cpp:215`).
  **Root cause (verified in `tmp/examples/CORC/corc.bfm`):** `ThinShellFactory::create` builds the
  aggregate "tape" sideset by **pointer-copying** the shell sidesets' facets
  (`cl_ThinShellFactory.cpp:119-124`, `collect_facets` `:463-514`), so `shell_03/04/05` and `tape`
  share the same `Facet` objects on the fresh mesh. `ProtoMesh::populate_facet_data` iterates
  **sidesets** (`cl_ProtoMesh.cpp:385-427`), so the shared facets are written twice — the file holds
  16761 facet records but only 12285 unique IDs (4476 duplicated; tape block byte-identical to the
  shell block, verified with h5py). On load, `create_facets` materializes **two distinct facets per
  duplicated ID** (`mFacetMap[id]` keeps the last, `cl_ProtoMesh.cpp:999`); positional `tGeo`
  stamping keeps all rank-0 asserts green. The distributor then ships both copies (facet bitset is
  index-based); the worker's ID-keyed `mFacetMap` resolves both records to one object → second
  `allocate_facet_container()` → abort. Serial load is also **silently wrong**: 16761 vs 12285
  facets; shells own clones instead of aliasing tape's facets (fresh mesh: shared facets report
  `sideset_id()==tape` after `set_sideset_ids`, `cl_Mesh.cpp:855/1753` — which is why the fresh MPI
  path works: one copy shipped, worker shells arrive empty and are deleted,
  `cl_ProtoMesh.cpp:1448-1451`).
  **Refined root cause (2026-07-01, Christian's DIAG in `collect_facets_from_sidesets`):** the
  CutFactory doubling/stitch is **NOT** at fault — `duplicate_and_relink_facets`
  (`cl_CutFactory.cpp:1953`, 2N temp facets for the cohomology cut) and
  `restore_thin_shell_sidesets` (`:876`) are symmetric; the DIAG series reconstructs exactly
  (10126 → +7582 doubling → +4454 cuts → −7582 restore). The duplication is a **triple reference
  created by `ThinShellFactory::create`** (all three born in commit `2dd4f090`, 2025-09-11):
  (1) the original `shell_NN` sidesets keep stale pointer-copies (`collect_facets` copies, never
  clears the source, `cl_ThinShellFactory.cpp:463-514`); (2) the aggregate tape sideset (same
  pointers) is pushed into `mMesh->sidesets()` (`:337`); (3) `ThinShell` wraps that *same* sideset
  and `Mesh::collect_facets_from_sidesets` re-adds `tThinShell->facets()` (`cl_Mesh.cpp:671-677`)
  — redundant with (2) from birth. DIAG: B−A = A−C = 7582 (facets 3× in `mFacets`). Fresh MPI
  tolerates this only by accident: `update_facet_indices` writes the shared object's index 3×
  (last position wins) and the distributor's bitset keys on `facet->index()` → one bit per unique
  object. Ghost + periodic sideset facets are freshly created (`:1865`, `:2342`) — no aliasing.
  **Interim band-aid (Christian, 2026-07-01):** `unique(mFacets)` in `collect_facets_from_sidesets`
  — dedupes the container but (a) does NOT fix the `.bfm` (`populate_facet_data` iterates
  *sidesets* → duplicate records → parallel-load crash survives) and (b) `unique(Cell<T>)` SORTS
  (`cl_Cell.hpp:461-470`) → `Facet*` by heap address → facet order/indices become
  **nondeterministic** across runs (breaks the "container order = saved order" equivalence
  assumption). Remove once the structural fix lands.
  **Fix (recommended, at the source — supersedes the loader options; refined 2026-07-01 pm):**
  - [x] `ThinShellFactory::collect_facets`: **move** membership instead of copying — in the copy
    loop, `tSideSet->reset_facet_container()` after each source `shell_NN` is copied (clears
    facets + `mFacetMap` + counter; deletes nothing). **Done (Christian, 2026-07-01).** No later
    reconstruction is needed or possible: the facet objects are repurposed *in place* into the
    tape facets (masters/slaves relinked to the extruded layer elements, edge containers reset
    `:247`, `set_sideset_id(tape)` `:250`) — the pre-extrusion shell facets cease to exist as such
    regardless of the container. A fresh rerun re-reads the `.msh`; a `.bfm` reload wants the
    extruded state. Verified: the source `tSideSets` are never read again in `create()` after
    `:124`; protoshells don't share sidesets.
  - [x] **Keep the empty `SideSet` husks in `mMesh->sidesets()`** (clear, don't delete):
    `read_domain_types` re-applies input domain types by **unguarded** `mMesh->sideset(tID)`
    (`cl_MaxwellFactory.cpp:338`) on every run incl. reload — presence required, facets not.
    Empty-safety verified: `SideSet::element_type()` → `ElementType::EMPTY`
    (`cl_SideSet.hpp:284-294`); `create_hanging_edges_and_facets` only flips their domain type
    (`:1254-1257`); finalize loops + `update_facet_nodes` (`cl_Mesh.cpp:675`) +
    `DofData` face pass (`cl_FEM_DofMgr_DofData.cpp:2347`) all guard `number_of_facets()==0`;
    workers prove full absence is tolerated post-distribution.
  - [x] **`.bfm` round-trips empty sidesets. Done (Claude, 2026-07-01).**
    `ProtoMesh::create_sidesets( const bool aKeepEmpty = false )` — keeps zero-facet sidesets when
    set (`cl_ProtoMesh.cpp` push/delete loop); distributor call unchanged (default `false`);
    `BfmFile::load_facet_data` passes `true`. Save side needed nothing (GroupData count 0, no
    facet records).
  - [x] **ExodusWriter count/populate filter mismatch. Done (Claude, 2026-07-01).** `ex_put_init`'s
    sideset count (`cl_Mesh_ExodusWriter.cpp:132-140`) filtered on `is_hidden` only while
    `populate_sidesets` (`:454`) also requires `number_of_facets() > 0` — visible-but-empty husks
    would have declared more sidesets than get written. Count now uses the same filter.
  - [x] Remove the redundant `mThinShells` loop in `Mesh::collect_facets_from_sidesets`
    (`cl_Mesh.cpp:648-651/671-677`) — the aggregate sideset is in `mSideSets` since `2dd4f090`.
    **Done (Christian, 2026-07-01).**
  - [x] Revert the `unique()` band-aid + DIAG prints. **Done (Christian, 2026-07-01).**
  - [x] Post-fix sanity (partial): fresh serial CORC run ✓ → new `.bfm` verified duplicate-free
    (52490/52490 unique; shells count 0; tape 7582) ✓ → **worker `create_facet_extra` crash gone** ✓.
    Full serial reload + `mpirun -n 2` reload still pending behind D25.
  - Result: fresh rank-0 mesh == loaded rank-0 mesh (empty husks + tape owns facets), workers
    unchanged, `mMesh->facets()` naturally unique, `.bfm` duplicate-free, loader dedupe unnecessary.

- [◐] **D25 — CRITICAL (load-path divergence) — `Topology::run()` on the enriched mesh builds a
  different type map than fresh.** Found 2026-07-01 via serial-reload postprocessor crash
  ("Wrong element type (is penta6, expect tet4)", `cl_FEM_Postprocessor.cpp:635`); root-caused,
  **fix implemented (pending build+test)**. Fresh path: `mTopology->run()` executes inside
  `create_cuts` (`cl_MaxwellFactory.cpp:749`) **before** enrichment → map sees input blocks/sidesets
  only ("the thin shell blocks are not part of the topology", `:1936`); cut sidesets stay
  `DomainType::Default` on the mesh forever (created post-scan, never registered — no consumer of
  `groups(Cut)` exists). Load path (D23) ran plain `run()` on the **enriched** loaded mesh:
  (a) buffer block 23 (`Buffer`) + layer blocks 20-26 entered the map → Air postprocessor
  (`Air∨Buffer` rule, `cl_MaxwellPostprocessor.cpp:167`) mixed PENTA6 into a TET4 list → assert;
  (b) `detect_sideset_types` **mutated** the loaded mesh (cut sidesets `Default→Inactive`);
  (c) `groups(Buffer)` became non-empty → the (fresh-dead) buffer node-flag loop
  (`cl_MaxwellFactory.cpp:1815`) went live → hanging-source relink diverged. Likely also implicated
  in the parallel reload dof-visibility crash (node 19863 φ, cut duplicate) — retest after this fix.
  **Fix (Christian's synthesis design — no persistence needed):** `Topology::run_on_enriched_mesh()`
  (`cl_Topology.{hpp,cpp}`) = `collect_enrichment_ids()` (exclusion sets from
  `mMesh->thin_shells()`: layer+buffer block ids, tape + ghost sideset ids) + the normal
  scan/select with `key_exists` skips + **no `detect_sideset_types`** (types already final on
  reload; keeps cuts `Default`). φ rule confirmed pure (`select_blocks:382-394`:
  φ=Air∨Buffer∨Ferro, non-φ=Conductor∨Coil — Christian's rule). Exclusion sets empty on the fresh
  path → zero behavior change there. Call site: load branch now `run_on_enriched_mesh()`
  (`cl_MaxwellFactory.cpp:442`); `synchronize_maps()` unchanged.
  - [x] **Build + reload tests PASS (Christian, 2026-07-01):** serial reload (postproc ✓),
    `mpirun -n 2` ✓, `-n 4` ✓, `-n 8` running smoothly. The parallel dof-visibility crash
    (node 19863 φ, key 39726) was indeed the same D25 divergence — resolved by this fix; the
    "step-6 sideset visibility hole" hypothesis is retired unless a future partition retriggers it.
    Fresh serial confirmed earlier the same day (wrote the clean `.bfm`); fresh `-n 2`/`-n 4`
    control runs still worth a pass when convenient.
  - [ ] **Parked question (Christian):** `cl_MaxwellFactory.cpp:1813-1819` claims buffer nodes
    "need their nodes flagged for hanging-node detection", but `groups(Buffer)` has always been
    empty on the fresh path (map predates the TSF) — dead code. Should it be live? (Physics
    question; the D25 fix deliberately keeps it dead on both paths.)

> **Equivalence VERIFIED EQUIVALENT (both auditors, conditional on D22):** node/edge/face/facet/element
> container order = saved order (no load-side re-sort; `update_node_indices` reassigns in container order);
> edge reconstruction is Map-**lookup**-only (no iteration-order reorder); DOF numbering is keyed by entity
> **ID** not container index (`node_dof_id = id*nDof+…`), so even index drift wouldn't renumber DOFs; edge
> directions recomputed from `original()->id()`; facet/face orientations saved + `compute_orientation`
> won't overwrite; hanging T-matrix sources saved/restored in container order; periodic skip of
> `set_entity_dependencies` is safe because the post-dependency slave sources ride in `/hanging`; thin-shell
> node indices use the same `original()->index()` pattern fresh and load; scaling/checksum correct.
> **Residual (confirm via R12, not hard breaks):** bit-exact sparse-matrix *entry* ordering (hanging-DOF
> enumeration order may differ even when DOF IDs match); periodic *edge-cascade* T-matrices not traced
> end-to-end; MPI `comm_size>1` not traced (serial rank-0 audit only).

> **Verified CORRECT on the integration path (no action):** enriched-save gating + ordering (save is
> `mComputeCohomologies`-gated and runs *after* periodicity `set_entity_dependencies` + hanging
> creation — R7 invariant honored); periodicity `update()` restored-guard; `read_domain_types` re-applies
> input types on both paths (input authoritative — validate if input may drift, D21); `create_edges_and_faces_on_mesh`
> self-gates on `edges_exist()/faces_exist()`; `create_cuts` `mComputeCohomologies`-gated.
>
> **Resolved earlier this session (for the record):** thin-shell material silent-loss (now per-shell
> `materials` in `/thinshells`), `find()`/map-ownership (removed), vertex id/node swap, `create_group("curves")`,
> curve labels restored, curve-map rebuilt by `finalize`. **WONTFIX:** ghost-sideset `0` sentinel — a
> sideset id of `0` cannot occur (it would crash Exodus), so `0`=="no ghost" is safe.

**Still unimplemented:** `format_version` attr; **global variables** (row 13); call-site wiring of the
fields restart file. **Out of scope by design (D21):** regular block→material (input authoritative) and
physical fields/time (separate `Mesh::save_fields`/`load_fields` file, not the `.bfm`).
**Implemented (no longer missing):** node duplicates; hanging entities (T-matrix / static-condensation,
R-row 9); control points (D12); abstract + orphan nodes (ID-keyed); periodicity (R7); thin shells
(per-shell blocks/thicknesses/materials); vertices; curves + segments; checksum (base, D16b);
physical fields + time cursor (separate restart file). **`Mesh::save/load` is wired**
(`Mesh::save:354`, MaxwellFactory `read_mesh`/`create_magnetic_kernel`), and the **D13–D21 integration
defects are resolved** (2026-06-30) — see §4.0. R10/R11 remain open only for the field-restart call-site
wiring + globals + an end-to-end test (R12).

- [x] ~~**R0 — Clear the legacy writer/reader bug catalogue first.**~~ **OBSOLETE (2026-06-28):**
  the legacy `cl_Mesh_HDF5Reader`/`HDF5Writer` are deleted from the tree, so Appendix A (W1-W6,
  H1-H3) is moot and has been removed. The reused HDF5 vlen patterns now live in the audited
  `hdf5::Dataset` helper (`cl_HDF5_Dataset.hpp`), whose own defects are tracked as D5-D7.

- [◐] **R1 — Decide and freeze the HDF5 schema** (see §6). **DE-FACTO (2026-06-28):** an
  implemented schema already exists in `BfmFile`, but it diverges from §6 (e.g. hanging stored as
  per-entity-type subgroups with vlen `topology`/`types`/`weights`; node duplicates as a vlen
  `nodes/duplicates`). §6 is now **aspirational/partly stale** — reconcile §6 to the BfmFile layout,
  and still add the `format_version` attribute (not yet written).

- [x] ~~**R2 — Build `to_proto` extractors on the live `Mesh`**~~ **DONE** — realized as
  `ProtoMesh::populate_*` driven by `BfmFile` (element **type** + ID-keyed topology included). See §4.0.

- [x] ~~**R3 — New `ProtoMeshHDF5Writer`**~~ **DONE** — realized as `BfmFile::save_*` + the
  `hdf5::Dataset` vlen helper (not a separate writer class). ID-keyed throughout. See §4.0.
  *(Carries the open D8-D10 row-index regressions.)*

- [x] ~~**R4 — New `ProtoMeshHDF5Reader`**~~ **DONE** — realized as `BfmFile::load_*` driving the
  existing `ProtoMesh` builders; stored element type consumed (no `element_type_from_numnodes`
  guessing). See §4.0.

- [ ] **R5 — Round-trip test on a base (non-enriched) mesh**: msh → save → load → save, assert
  byte/topology stable and checksum preserved. Gate before touching enrichment. (after: R4)

- [◐] **R6 — Serialize physical fields + time + global variables** (rows 12,13,14). **Fields + time DONE
  2026-06-30**, as a **separate** restart file (not the `.bfm`), matching the D21 split: `Mesh::save_fields`
  / `Mesh::load_fields` (`cl_Mesh.cpp:2946-3018`) write a standalone HDF5 — `meta`:{`timestep`,`timestamp`,
  `checksum`}, `fields`:{per-label `data` + `labels[]`/`types[]`}; load **verifies `checksum == mesh.checksum()`**
  (ties the field file to a mesh identity — consistent now that the base checksum is canonical, D16b) and
  recreates fields by label+`EntityType` via `create_field`. Now **wired into the controller**
  (`Controller::save_fields`/`load_fields`, `hphirun.cpp:78,113`). The `.bfm` itself carries no fields
  (by design, D21).
  - [ ] **R6a — Global variables (row 13). STILL MISSING.** `Mesh::save_fields` writes time+fields but
    NOT global variables (total current/energy, `cl_Mesh_GlobalVariable`). Add `{ids, labels, values}` to
    the restart file's `meta`/`globals`. (Christian agreed to add — low complexity, mostly-derived data.)
  - [ ] **R6b — Coupled circuit dynamic state. STILL MISSING (needed for CORC restart).** The electrical
    circuit does its own time-stepping (`shift()` at step end, MNA/BDF state). On restart it is re-created
    from `input.conf` but its time-*T* state (inductor currents, node voltages, its history) is NOT
    restored → the coupled FEM↔circuit run resumes with a cold circuit. Persist the circuit's dynamic
    state alongside the fields. (Christian agreed to add.)
  - **Deliberately NOT persisted (decided 2026-06-30, per OpenFOAM/COMSOL practice):** BDF multistep
    history (`IWG_Timestep::mFieldData`/`mH`) and the controller's adapted Δt (`mDeltaTime`). Restart is a
    **state restart** (loaded fields = new IC) with a **cold-started** integrator (order ramps from 1) and
    Δt **re-read from the (possibly edited) input file** — matching the "it didn't converge, lower Δt,
    continue" use case. TODO: confirm `IWG_Timestep` cold-starts cleanly on a load (seed history stages
    from the loaded field / restart the order ramp) rather than reading stale stage buffers.
  - **Sync fix (in progress):** the load path must **distribute-only** (rank-0 loaded field → all ranks),
    not `synchronize_fields` (collect+distribute) — `collect` clobbers the loaded field for
    non-rank-0-owned nodes in parallel (invisible on serial). Two-mesh field save is a non-issue (the
    thermal mesh is extracted from the magnetic one).

- [x] **R7 — Periodicity save + a dedicated reload contract** (row 15). **DONE 2026-06-27** (audited
  Claude + Grok + Codex; three CRITICALs found and fixed; see devlog). Built in `BfmFile`
  `save/load_periodicity_data` (`cl_Mesh_BfmFile.cpp:1639-1906`) + `from_proto` link path.
  - **Plane** stored as a `2×3` `id_t` matrix (master/slave rows × 3 plane-node cols);
    `load_periodic_planes` `set_size(3)`s the proto Cells. (Earlier `3×2`/unsized-Cell bugs fixed.)
  - **Pairs stored by authoritative partner, not by list position.** `save_periodic_*` writes
    `(node->id(), node->periodic()->id())` per column, so the loaded `mMasterX(k)↔mSlaveX(k)` are
    true partners — sidestepping the "lists are independently ordered" hazard
    (`periodicity.md:88`). ⇒ `mMaster/SlaveSideSets` are **not** serialized (unused by `from_proto`).
  - **Reload via O1 option (i):** `create_periodicitiy(true)` → `from_proto(..., aCrosslinkEntities=true)`
    → `crosslink(Periodicity*)` (`PeriodicityFactory.cpp:1426`) sets bidirectional `set_periodic()`
    **and** flags each populated entity type (NODE/EDGE/FACE/FACET). No geometric re-match.
  - **Periodic DOF constraints** (slave `add_source(master)`) are **not** re-derived on load — they
    ride back through the hanging-entity data (`save_hanging_nodes` saves every `is_hanging()`
    entity), provided the enriched `.bfm` is written *after* `set_entity_dependencies()` (R10 note).
  - **Maxwell read-side gate done:** `cl_MaxwellFactory.cpp:481-489` now runs
    `update()`+`set_entity_dependencies()` only when `! node_pairs_restored()`. **Write-side still
    pending** — the enriched `.bfm` save (R10) must be wired and must sit *after* the periodicity
    block so the slave source containers exist to be saved.

- [x] **R8 — ThinShell + ghost serialization** (rows 16,17). **DONE 2026-06-30** — realized **not** via
  proto `ThinShellData` but as a direct `BfmFile` path: `save/load_thinshell_data` writes `/thinshells`
  {`ids`, `ghost` (`gNoID`→0 sentinel), vlen `blockids`/`thicknesses`, flat per-shell `materials`} and
  `load_thinshell_data` builds `ThinShell` objects directly (the `.bfm` stores the *realized* mesh, not
  the extrusion recipe). Node indices recomputed via `original()->index()` (factory-equivalent). Load-side
  double-creation handled in `MaxwellFactory::create_thinshells` (D18). *Cleanup left:* verify the dead
  `save_thin_shell_facets`/`"tapes"` path is gone (low priority).

- [x] **R9 — Block↔material map. CLOSED by D21 decision (2026-06-30) — won't implement.** The `.bfm` does
  **not** persist regular block→material; input file is authoritative and rebuilds the assignment every run
  (`collect_material_labels_from_domains`). The only material that rides in the `.bfm` is **thin-shell layer**
  material (in `/thinshells`, hydrated on load — D18), because the input-domain pass skips
  `DomainType::ThinShell`. `save_material_map` is therefore not folded into the schema.

- [x] **R10 — Enriched mesh write in `MaxwellFactory`. DONE 2026-06-30.** `mMesh->save(mMeshPath)`
  (`cl_MaxwellFactory.cpp:469`), `mComputeCohomologies`-gated, after enrichment + periodicity
  `set_entity_dependencies` + hanging creation (R7 invariant). The D13–D21 blockers (D14 double-save, D15
  empty group ids, D20 checksum) are resolved. *Separate follow-up:* the controller does not yet call
  `Mesh::save_fields` at the restart cursor (field-restart write — R6/field-wiring, not the mesh write).

- [◐] **R11 — Reload orchestration in `read_mesh` + `create_magnetic_kernel`.** **Mesh reload DONE
  2026-06-30** — `read_mesh` (D16) returns on all paths, `mComputeCohomologies=false` on `.bfm` load skips
  `create_cuts`; `create_thinshells`/`create_curves`/`create_periodic` are now load-aware (D18/D19 + the
  periodicity guard); `mMaterialBlockAssignment` hydrated from input + loaded shells (D17/D18). **Still
  open:** wire `Mesh::load_fields` into the controller so it resumes from `time_step`/`time_stamp`, and
  persist global variables (R6). (after: R10)

- [x] **R12 — End-to-end restart test. DONE 2026-07-01 (Christian).** CORC (periodic + thin-shell +
  cohomology cuts): fresh run → `.bfm` save → relaunch from the sidecar — passes in **serial and on
  2/4/8 MPI procs** after D24 (facet aliasing) + D25 (topology-map divergence) landed. Format
  documented in `src/mesh/doc/bfm_file_format.md`.

---

## 5. Open Design Questions (not silently decided)

- **O1 — Periodicity reload API. RESOLVED 2026-06-27 → option (i).** `from_proto(..., aCrosslinkEntities)`
  sets bidirectional `periodic()` pointers + entity-type flags directly via `crosslink`
  (`PeriodicityFactory.cpp:1426`); no node-pair-backup replay. One code path, no geometric re-match.
- **O2 — Does any post-reload consumer still need a geometric `update()`? RESOLVED → no.**
  `set_entity_dependencies` is the consumer of `periodic()`; with pairs serialized by partner ID and
  crosslinked on load, and slave source containers restored via the hanging-entity data, the restored
  path skips `update()` entirely (`cl_MaxwellFactory.cpp:481-489`, guarded by `node_pairs_restored()`).
- **O3 — Field set to persist for restart. PARTLY DECIDED 2026-06-30 → persist ALL mesh fields.**
  `Mesh::save_fields`/`load_fields` (`cl_Mesh.cpp:2946`) serialize **every** `mFields` entry (by label +
  `EntityType`) plus `time_step`/`time_stamp`, in a separate checksum-guarded file. **Still open:** BDF
  history lives in the IWG (`cl_IWG_Timestep.cpp:258-285`), **not** on the mesh, and is **not** serialized —
  so a restart is **first-order** (resumes from in-mesh field state), not bit-exact, unless IWG history
  (`mH`, shifted fields) is also persisted. Whether that matters is the deciding correctness question for
  **R12** — confirm the end-to-end restart matches an uninterrupted run within tolerance, and if not,
  serialize the IWG history too.
- **O4 — Element-type encoding. RESOLVED → per-block type.** The `.bfm` stores the `ElementType` enum
  per **block** (`/blocks/types`, single-type per block), not per element; element topology is just node
  IDs. Since each block is single-type, this fully disambiguates TET20/HEX20 (the row-3 ambiguity) without
  a per-element field. (Verified in §6.1.)
- **O5 — Control points. RESOLVED 2026-06-28/30.** `BfmFile` saves/loads control points (entities +
  per-element incidence); the incidence-save bug **D12 is fixed** (`tTopo.save()`). One small caveat to
  confirm at leisure: `ProtoMesh::create_control_points` reads `tCoords(2,p)` unconditionally → OOB for a
  **2D** control-point mesh (shared with the distributor) — only an issue if 2D B-spline/tensor meshes
  actually occur; guard with `tNumDim==3` if so.
- **O6 — Interpolation type (row 20).** Confirm Maxwell never uses non-Lagrange interpolation; if it
  can, the format must store it (it is not derivable from the mesh). (ref: `Mesh_Enums.hpp:94-106`)
- **O7 — Retire or keep the legacy `HDF5Writer`/`HDF5Reader`?** **RESOLVED (2026-06-26): retire.**
  `cl_Mesh_HDF5Reader`/`HDF5Writer` deleted; the ProtoMesh-backed `BfmFile` is the sole HDF5 mesh
  path. Appendix A's reader bugs are therefore no longer live (implementation cautions only for the
  reused vlen patterns).
- **O8 — Curve coordinate recomputation. MOOT (resolved by the implementation).** The new `BfmFile`
  curve path does **not** recompute coords — `save_curve_data` stores segment topology as **node IDs** (+
  arclength), and `load_curve_data` rebuilds segments by attaching the already-loaded mesh nodes
  (`mProto->node(id)`), which carry their coords. No `CurveFactory` coordinate recomputation, so the
  enriched/cut-duplicate concern doesn't arise.
- **O9 — Re-enable gating. PARTLY DECIDED 2026-06-30.** The enriched write is **on**, gated by
  `mComputeCohomologies` + the base-checksum match (D16b/D20) — i.e. re-dump only when the source mesh
  changed; no user feature flag added (always-on). **Still open (low priority):** whether to expose a
  `mesh{write_bfm}` flag to suppress the dump, and the `format_version` root attr is **not yet written**
  (§6.4) so a future compatibility break isn't yet detectable.

---

## 6. `.bfm` HDF5 Schema

All cross-references are **by entity ID**, never by array index — they resolve via `Mesh::node(id)`,
`edge(id)`, `face(id)`, `facet(id)`, `sideset(id)`, `block(id)`, the same maps `ProtoMesh` uses
(`cl_ProtoMesh.hpp:186-195`). Because IDs are stable across finalize/reorder
(`mesh_contracts_and_invariants.md §4`), the format is robust to any future reordering, unlike the
legacy index-keyed writer. Variable-length topology uses the `hdf5::Dataset` vlen helper.

> **MOVED (2026-07-01):** the format is verified (R12) and now documented as the authoritative
> **`src/mesh/doc/bfm_file_format.md`** — group-by-group dataset reference (reconciled against the
> code AND an h5py dump of a real file), encoding conventions, positional-membership rules, and
> the reload contract (stored vs recomputed vs synthesized). §6.1-6.4 below are kept as the
> historical working notes; on any format change, update the module doc, not this section.

### 6.1 Implemented in `BfmFile` today (verified 2026-06-27)

Groups in `save()` order. `(opt)` = written only when non-empty / guarded on load by `group_exists`
or `dataset_exists`.

```
/meta
    dimensions     scalar         uint       (num spatial dimensions)
    entities       scalar         index_t    (entity count)
    groups         scalar         index_t    (group count)
/blocks
    ids            [B]            id_t
    labels         [B]            string
    elements       [B]            uint       (element count per block)
    types          [B]            uint       (element type per block — single-type per block)
    domains        [B]            uint       (DomainType)
    tags           [B]            ...        (physical tags)
/sidesets                         (same dataset set as /blocks; domains incl. Ghost = 37)
    ids, labels, elements, types, domains, tags
/nodes
    ids            [N]            id_t
    coords         [N, dim]       real       (saved transposed → N×dim on disk, HDFView-readable)
    abstract       [*]            id_t   (opt) abstract-node IDs (category membership)
    orphaned       [*]            id_t   (opt) orphaned-node IDs
    duplicates     [*][var]       id_t   (opt) per original: [original_id, dup_id, dup_id, ...]
/elements
    ids            [E]            id_t
    physical       [E]            ...    (opt) physical tags
    topology       [E][var]       id_t       (node IDs)
/facets
    ids            [Fc]           id_t
    physical       [Fc]           ...    (opt)
    topology       [Fc][var]      id_t   (opt) node IDs — only if aSaveTopology (save passes false)
    elements       [Fc][var]      id_t       (master[/slave] element IDs)
    indices        [Fc][var]      uchar      (idx_on_master, idx_on_slave, orient_on_slave)
/edges
    ids            [Ed]           id_t
    topology       [Ed][var]      id_t       (node IDs)
/faces
    ids            [F]            id_t
    elements       [F][var]       id_t       (master[/slave] element IDs)
    indices        [F][var]       uchar      (idx_on_master[, idx_on_slave, orient_on_slave])
/control_points
    ids            [Cp]           id_t
    coords         [Cp, dim]      real       (saved transposed)
    topology       [Cp][var]      id_t       (per-element control-point incidence, IDs)
/hanging                          (opt; one subgroup per entity type with >0 hanging)
    /nodes | /edges | /faces | /facets | /control_points
        ids        [H]            id_t       (the hanging targets)
        topology   [H][var]       id_t       (source entity IDs)
        types      [H][var]       suint      (EntityType per source)
        weights    [H][var]       real       (T-matrix weights)
/periodic                         (opt; only if has_periodicity)
    planes         [2, 3]         id_t       (row0 master / row1 slave × 3 plane-node IDs)
    nodes          [2, *]         id_t   (opt) per col: (id, periodic()->id) — saved transposed
    edges          [2, *]         id_t   (opt) authoritative partner pairs
    faces          [2, *]         id_t   (opt)
    facets         [2, *]         id_t   (opt)
```

### 6.2 Implemented since 6.1 was written (2026-06-30)

Now in the `.bfm` (move into §6.1 on the next consolidation):
- `/meta/checksum` — written by `save_meta_data`; read self-contained by `BfmFile::checksum()`. The
  stored value is the **base** (pre-enrichment) scaled-mesh checksum (D16b), so the sidecar re-dump gate
  in `MaxwellFactory::read_mesh` is meaningful across runs.
- `/vertices` {`ids`, `nodes`} (row 11).
- `/curves` {`ids`, `labels`, `sidesets`[2,n], `types`, `closed`, vlen `edges`/`segments`/`lengths`/`topology`}
  (row 10) — segments folded in (no separate `/segments` group).
- `/thinshells` {`ids`, `ghost`, vlen `blockids`/`thicknesses`, flat `materials`} (rows 16/17). Per-shell
  **materials** stored as a flat `Cell<string>` (authoritative labels), de-flattened on load by the
  blockids row lengths.

### 6.3 Deliberately NOT in the `.bfm` (D21 decision)

- **Regular block→material** — input file is source of truth; rebuilt every run by
  `collect_material_labels_from_domains`. (Thin-shell *layer* materials are the exception — they ride in
  `/thinshells` because the input domain pass skips `DomainType::ThinShell`.)
- **Physical fields + time cursor** — a **separate** restart file via `Mesh::save_fields`/`load_fields`
  (`cl_Mesh.cpp:2946`): `meta`:{`timestep`,`timestamp`,`checksum`}, `fields`:{per-label `data`,
  `labels[]`, `types[]`}; checksum-guarded to the mesh identity. The `.bfm` is the enrichment cache only.

### 6.4 Still missing (to add before the format is frozen)

- `/meta/format_version` (uint) — version gate for the reader (O9). Not yet written.
- **Global variables** (row 13) — neither the `.bfm` nor `Mesh::save_fields` persists them yet
  (`save_fields` writes time + fields only). Add to the fields restart file: `{ids, labels, values}`.
- **Call-site wiring** — `save_fields`/`load_fields` exist but are not yet invoked by the controller at
  the restart cursor (R10/R11).

---

## 7. Definition-of-Done Checklist

- [x] Every "missing knowledge" item mapped to a concrete reconstruction step or open question:
  periodic pairing → R7+O1/O2; duplication config → rows 7/§3.1 + R3 (ID-keyed); ghost elements →
  rows 16/17 + R8; curves/terminals → row 10 (kept); fields @ last step → row 12 + R6.
- [x] Gap table with (a)/(b)/(c) classification, one code citation per row (§3).
- [x] Each claimed gap backed by writer/reader (or proto/factory) citation, not assumption.
- [x] Ordered todo with dependencies (§4).
- [x] Open design questions logged, not decided (§5).
- [x] Proposed HDF5 schema with per-entity fields and ID-keyed periodic/duplication state (§6).

---

## 8. Audit Trail

- Full thread: `tmp/ai_exchange/meshfile_refactor.md` (Claude draft → Codex audit → Grok third voice
  → Claude resolution). Swept per protocol §10 after this plan + a devlog capture.
- **Codex (high confidence)** corrected Row 3 to must-save (TET20/HEX20 collision, verified
  `meshtools.cpp:1181-1184`), flagged Row 9 non-equivalence and Row 17 ThinShell-object gap, and the
  pervasive index-keying fragility.
- **Grok (third voice)** established the periodicity reload is a **hard error** not silent loss
  (verified `PeriodicityFactory.cpp:320, 359-370`), that `to_proto` omits master/slave sidesets, and
  that `ProtoMesh::create_periodicitiy` never wires the periodic-finalize contract — driving R7/O1/O2.
- All three findings independently re-verified by reading the cited code before inclusion.

---

## Appendix A — REMOVED (2026-06-28)

The legacy `HDF5Writer`/`HDF5Reader` bug catalogue (W1-W6, H1-H3, M1-M2) was deleted because the
legacy `cl_Mesh_HDF5Reader`/`HDF5Writer` files no longer exist in the tree (O7 → retire). Any
remaining references to "Appendix A" elsewhere in this plan are historical. The one item that was a
*shared* utility — W6, `get_array_size` polarity in `src/io/HDF5_Tools.cpp` — should be re-checked
against the current `src/io/hdf5_tools.hpp` if that helper is still used; it is **not** part of the
BfmFile/`hdf5::Dataset` path. The legacy round-trip test ideas fold into R5/R12.

---

## Appendix B — Decision: store edge/face/facet links, do not re-derive from node topology

**Date:** 2026-06-23
**Status:** Audit only — no source changed. Confirms and pins down the keying choices behind
rows 2/3/5/6/7/9/12/15 and §3.1. Confidence: **high** (grounded in the live `ProtoMesh`
reconstruction and the element-template invariance result below).

> **SUPERSEDED IN PART (2026-06-27) — see §4.0 D3.** The conclusion below ("store the element→
> edge/face links via `ElementExtra`") assumed re-derivation meant *re-creating* edges/faces with
> fresh IDs. It does not have to: `reconstruct_edge_connectivity` attaches the **stored** edges by
> node-key and `reconstruct_face_connectivity` rebuilds element→face from the stored `FaceData`
> master/slave + index — both **preserve the stored IDs**, so the ID-stability argument no longer
> forces serialization. Net: `ElementExtra` is **not** serialized; element→edge/face incidence is
> reconstructed on load and the rest (directions, inverses, maps) comes from `Mesh::finalize()`.
> What still must be stored is unchanged: node/element/edge/face/facet **entities + IDs + the
> primaries** (coords, element type+topology, edge node-topology, face/facet master/slave+index).

### B.1 The question

Must the `.bfm` store element↔edge, element↔face, and element↔facet connectivity, or can a reader
rebuild it from node topology alone — apply each element's template to its node list and re-key
edges/faces by their node sets?

### B.2 Answer

**Store them, ID-keyed.** Edge/face *incidence* is re-derivable, but the edge/face/facet **IDs** and
the enrichment-derived attributes are not — and everything else cross-references those IDs. This is
already the call the battle-tested `ProtoMesh` (MPI distributor) backend makes: it ships every entity
whole and re-derives nothing.

### B.3 Why the incidence alone *is* re-derivable (so the temptation is real)

Element edges/faces are **derived state**. The templates (`get_nodes_of_edge`/`get_nodes_of_facet`,
e.g. `cl_Element_TET4.hpp`) read `mNodes[...]` fresh on every call, and the accessors are pure
positional (`cl_ElementTemplate.hpp:581` `insert_node` ⇒ `mNodes[k]=node`; `:588` `node(k)` ⇒
`mNodes[k]`). So "local edge *i* = slots (p,q)" is a compile-time constant per element type. Node
duplication overwrites slots **in place** (`cl_CutFactory.cpp:1944`
`insert_node( mThinShellDuplicates(tNode->index()), k )`), so a duplicated edge keeps the *same local
slot pair* — only the node identities change (a→p, b→q). Edge identity is therefore a deterministic
function of the `original()`-keyed endpoint pair, and `finalize_edges()` (`cl_Mesh.cpp:981-997`,
directions from endpoint `original()->id()` `:1821-1854`) recomputes local H(curl) signs
deterministically. Hence the Gap Table classes edges/faces as **(a) relink on load** (Row 6, Row 22a).

### B.4 Why we store anyway — three reasons re-derivation breaks

1. **ID stability is load-bearing.** T-matrices (Row 9), periodicity pairs (Row 15) and restart
   fields (Row 12) all cross-reference edges/faces/facets **by ID**. Re-deriving assigns *fresh* IDs
   in factory order, dangling every such reference — silently reintroducing the index/order fragility
   §3.1 calls "the single most important correctness change" to remove.
2. **Edge/face-attached restart fields need the same IDs back.** Row 12 fields can be `entity_type`
   edge/face (e.g. `edge_h`), keyed to edge/face ID. Re-derived IDs force you to reproduce the *old*
   assignment to re-key — i.e. you store the IDs anyway, implicitly.
3. **Re-derivation is only valid *after* enrichment is baked into the node slots — which reload
   skips.** Regen is "free" only because the *post-duplication* element→node topology is treated as
   primary, and that topology is the cut/cohomology/thin-shell factory output that reload bypasses
   (`mComputeCohomologies=false`). The irreducible primaries that **must** be stored: element **TYPE**
   (Row 3 — without it no template can be selected; not a function of (dim,nnodes)), **ID-keyed
   post-duplication element→node topology** (Rows 2, 7), and the **duplicate↔original map** (Row 7).

### B.5 Facets — must store for two further (non-ID) reasons

- **Sideset/block membership** is a physical/BC grouping, never derivable from node geometry.
- **master/slave role + `orientation_on_slave`** is re-derivable only by re-running the
  `to_master_orientation` reconciliation (`fn_to_master_orientation.cpp`) across duplicated
  interfaces — re-invoking the exact thin-shell/cut logic reload is meant to skip. `FacetData` already
  carries `mMasterIDs/mSlaveIDs/mIndicesOnMaster/Slave/mOrientationsOnSlave`
  (`st_ProtoMesh.hpp:87-104`); keep them as Row 5 **(c)**.

### B.6 Precedent — the existing ProtoMesh already decided this

The MPI-distributor reconstruction stores and rebuilds by ID, re-deriving nothing:
- `EdgeData` = edge IDs + node topology; `create_edges` (`cl_ProtoMesh.cpp:456-480`) rebuilds by
  stored ID, never re-keying from element node-sets.
- `ElementExtra` carries `mEdgeData`, `mFaceData` **and** `mEdgeDirections`; `create_element_extra`
  (`:719-740`) wires `insert_edge(...,k)` *and* `set_edge_direction(k,...)` from stored data.
- `FaceData`/`FacetData` store master/slave IDs + indices + `orientation_on_slave`.

### B.7 Implication for R2/R4

Everything rests on the **`ElementData` extraction**: per-element `ElementType` enum (O4 — store the
enum, not (dim,nnodes)) plus **post-duplication** node topology **by node ID**. Get that right and
edges/faces *could* be relinked — but per B.4 we still serialize the
`EdgeData`/`ElementExtra`/`FaceData`/`FacetData` proto structs (schema §6) to keep IDs stable for the
T-matrix/periodicity/field cross-references. **Net rule: re-derive nothing on load; drive the existing
ID-keyed `ProtoMesh` builders.**

### B.8 Decision summary

| Entity | Incidence re-derivable? | IDs re-derivable? | Action |
|---|---|---|---|
| Element→node (+ type + dup map) | — (primary) | — | **Store, ID-keyed** (Rows 2/3/7) — irreducible base |
| Edge + element→edge + directions | yes | **no** | **Store** (referenced by T-matrix/periodicity/fields) |
| Face + element→face | yes | **no** | **Store** |
| Facet (master/slave/idx/orient + sideset) | no | no | **Store, ID-keyed** (Row 5) |
