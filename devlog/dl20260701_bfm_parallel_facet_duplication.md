# `.bfm` Parallel Load Crashes — Facet Aliasing (D24) + Topology Map Divergence (D25)

**Date:** 2026-07-01
**Purpose:** Root-cause and fix the `.bfm` reload crashes: (1) `mpirun` SIGABRT "container is
already allocated" (`cl_Vertex.cpp:216`) — D24; (2) parallel dof "Key 39726 not found" + serial
postprocessor "is penta6, expect tet4" — both symptoms of D25. **End state: `.bfm` reload runs
serial, 2, 4, and 8 procs.**
**Module:** `src/mesh` (`BfmFile`/`ProtoMesh`/`Distributor`), `src/fem/kernel` (`ThinShellFactory`),
`src/homology` (`Topology`), `src/fem/maxwell` (`MaxwellFactory`)
**AI:** Claude (Fable), with Christian implementing/reviewing each patch. Tracked as **D24/D25**
in `todo/meshfile_refactor_plan.md` §4.0.

## Symptom

`mpirun` on the CORC case with a `corc.bfm` sidecar present (i.e. the **load path**) aborts on a
worker rank:

```
Vertex::allocate_facet_container()          cl_Vertex.cpp:216   "mFacets == nullptr" assert
← ProtoMesh::create_facet_extra()           cl_ProtoMesh.cpp:1234
← Distributor::run()                        cl_Mesh_Distributor.cpp:215   (rank != 0 branch)
← Kernel::distribute_mesh → MaxwellFactory::create_magnetic_kernel → hphirun
```

`Facet` derives from `Vertex` (`cl_Facet.hpp:23`), so the neighbor-facet container lives on the
facet itself; the assert means the *same facet object* receives two allocation records.

## Root cause (confidence: high — every link verified against the file and the code)

1. **Fresh mesh aliasing (by design).** `ThinShellFactory::create` builds the aggregate "tape"
   sideset by **pointer-copying** the facets of the connected shell sidesets
   (`cl_ThinShellFactory.cpp:119-124`; `collect_facets` `:463-514` copies `Facet*`). So
   `shell_03/04/05` (set `Inactive`) and `tape` (`GeometryOnly`) share the same 4476 `Facet`
   objects. Ownership is unambiguous: the facets are stamped `sideset_id = tape`
   (`cl_ThinShellFactory.cpp:250`), re-confirmed by `Mesh::set_sideset_ids()` in finalize
   (`cl_Mesh.cpp:855/1753` — last containing sideset wins; tape is last).

2. **The save writes aliases as records.** `ProtoMesh::populate_facet_data` iterates **sidesets**
   (`cl_ProtoMesh.cpp:385-427`), so shared facets serialize once per membership.
   Verified in `tmp/examples/CORC/corc.bfm` with h5py: `facets/ids` has **16761 records, 12285
   unique** — 4476 IDs (52111+) appear at positions 0–4475 (shell_03/04/05) and again at
   12285–16760 (tape), with byte-identical `elements`/`indices` rows.

3. **The load materializes aliases as clones.** `ProtoMesh::create_facets` creates one `Facet` per
   record; `mFacetMap[id]` silently keeps the last (`cl_ProtoMesh.cpp:999`). The positional
   sideset stamp (`tGeo` → `geometry_tag` ≡ `block_id` ≡ `sideset_id`, `cl_Element.hpp:683-710`)
   keeps every rank-0 assert green, so the corrupted mesh (16761 facets) reaches the distributor.

4. **The worker collides on ID.** `Distributor` ships both clones (facet bitset is index-based,
   `select_entities` `cl_Mesh_Distributor.cpp:320-339`); the worker's `create_facets` again maps
   both records to one `mFacetMap` slot, and `create_facet_extra` resolves both neighbor records
   to the *same object* → second `allocate_facet_container()` → SIGABRT.

**Why the fresh MPI path never crashed:** on the fresh mesh each shared facet exists once in
`mMesh->facets()` with `sideset_id == tape`; the distributor ships one copy, worker-side
`shell_03/04/05` arrive empty and are deleted (`cl_ProtoMesh.cpp:1448-1451`), and thin shells find
their sideset by map. The `.bfm` round-trip is the only path that turns the aliases into clones.

**Serial load is silently wrong too:** the loaded rank-0 mesh has 16761 facets instead of 12285,
and the shell sidesets own independent clones instead of aliasing the tape facets. The serial CORC
reload "worked" only because the clones sit in `Inactive` sidesets.

## Recommended fix (option i in D24 — alias-preserving load; not yet implemented)

Keep the schema as-is and define **duplicate ID = alias membership**:

- `ProtoMesh::create_facets`: on a known ID, do not create a second facet — reuse the mapped
  object for membership only; only first-occurrence facets enter `mMesh->facets()` (shrink the
  container afterwards).
- Build BFM sideset membership from the saved per-sideset positional stream via `mFacetMap`
  (BFM-only path); leave the distributor's scan-based `create_sidesets` untouched (its stream is
  duplicate-free by construction).
- Finalize's `set_sideset_ids` then restores tape ownership exactly as on the fresh mesh.

Alternative (option ii, riskier): save each facet only under its owning sideset — but the loaded
rank-0 mesh would then lose the populated `Inactive` shell sidesets (divergence from fresh; every
consumer would need auditing). Details and trade-offs in `todo/meshfile_refactor_plan.md` D24.

**Also check:** `ThinShellFactory::create_periodic_sideset` — if it pointer-shares facets the same
way, the same duplication fires for periodic thin-shell meshes.

## Follow-up (same day): stitch-back exonerated, triple reference root-caused

Christian instrumented `Mesh::collect_facets_from_sidesets` (DIAG A = facets from sidesets,
B = + `mThinShells` loop, C = after `unique()`; CORC run, shells total 7582):

| call | when | A | B | C |
|---|---|---|---|---|
| 1 | initial finalize | 10126 | = | = |
| 2 | CutFactory, shells doubled | 17708 (= 10126 + 7582) | = | = |
| 3 | + cut sidesets, still doubled | 22162 (= +4454) | = | = |
| 4 | post ThinShellFactory | 60072 | 67654 | 52490 |

**The CutFactory doubling/stitch is correct:** `duplicate_and_relink_facets`
(`cl_CutFactory.cpp:1953` — originals moved to temp sidesets, 2N fresh facets open the domain for
cohomology) and `restore_thin_shell_sidesets` (`:876` — originals moved back, doubles deleted) are
symmetric; calls 1-4 reconstruct exactly (+7582 doubling, +4454 cuts, −7582 restore).

**The duplication is a triple reference from `ThinShellFactory::create`** (all three introduced
together in `2dd4f090`, 2025-09-11): stale pointer-copies left in the original `shell_NN` sidesets
(`collect_facets` copies, never clears, `cl_ThinShellFactory.cpp:463-514`); the aggregate sideset
pushed into `mMesh->sidesets()` (`:337`); and `ThinShell` wrapping that same sideset while
`collect_facets_from_sidesets` re-adds `tThinShell->facets()` (`cl_Mesh.cpp:671-677`) — hence
B−A = A−C = 7582 (each tape facet 3× in `mFacets`). Fresh MPI tolerates the triplication only by
accident (`update_facet_indices` last-position-wins + index-keyed distributor bitset → one bit per
unique object). Ghost and periodic sideset facets are freshly created (`:1865`, `:2342`) — the
aliasing is exclusively shell_NN ↔ tape ↔ ThinShell.

**`unique()` band-aid caveats:** it does not fix the `.bfm` (save iterates *sidesets*, so the
duplicate records — and the parallel-load crash — survive), and `unique(Cell<T>)` sorts
(`cl_Cell.hpp:461-470`), i.e. `Facet*` by heap address → nondeterministic facet order/indices.

**Recommended structural fix (supersedes the loader-side options above; details in D24):** move,
don't copy — clear the `shell_NN` containers in `collect_facets`; drop the redundant `mThinShells`
loop in `collect_facets_from_sidesets`; revert the band-aid; audit rank-0 consumers of post-TSF
`shell_NN` facets (workers already see those sidesets deleted).

## D24 fix (landed same day)

- `ThinShellFactory::collect_facets`: `tSideSet->reset_facet_container()` after each source
  sideset's pointer-copy (Christian) — membership *moves* to the tape aggregate; `shell_NN`
  remain as empty husks (required by `read_domain_types`' unguarded id lookup,
  `cl_MaxwellFactory.cpp:338`).
- `Mesh::collect_facets_from_sidesets`: redundant `mThinShells` loop removed (Christian).
- `ProtoMesh::create_sidesets( aKeepEmpty = false )`: BFM load passes `true` so the husks survive
  the round-trip; distributor path unchanged (Claude).
- `ExodusWriter`: `ex_put_init` sideset count aligned with `populate_sidesets`' filter
  (`!hidden && >0`) — visible-but-empty husks would have desynced the declared count (Claude).
- Verified: new `.bfm` has 52490/52490 unique facet records; shells count 0; tape owns 7582;
  worker `create_facet_extra` crash gone.

## D25: Topology map divergence on reload (the remaining two crashes)

After D24, the reload died (a) on 2 procs with `Key 39726 not found` in worker `link_node_dofs`
(node 19863 φ — a cut duplicate, hanging on {abstract 17006, original 13231}), and (b) in serial
in the postprocessor: "Wrong element type (is penta6, expect tet4)"
(`cl_FEM_Postprocessor.cpp:635`). Both had ONE root cause:

**The fresh path builds the topology type map BEFORE enrichment** (`mTopology->run()` inside
`create_cuts`, `cl_MaxwellFactory.cpp:749`) — thin-shell layer/buffer blocks and tape/ghost
sidesets are never in it ("the thin shell blocks are not part of the topology", `:1936`), and cut
sidesets stay `DomainType::Default` on the mesh. The D23 load branch ran plain `run()` on the
**enriched** loaded mesh: buffer block 23 (`Buffer`, PENTA6) entered the map → the Air
postprocessor (`Air∨Buffer` block rule, `cl_MaxwellPostprocessor.cpp:167`) mixed PENTA6 into a
TET4 element list; `detect_sideset_types` mutated cut sidesets `Default→Inactive`; `groups(Buffer)`
went live in the (fresh-dead) buffer node-flag loop (`:1815`) → hanging-source relink diverged →
the parallel dof-visibility crash.

**Fix (Christian's synthesis design — nothing persisted):** new
`Topology::run_on_enriched_mesh()` (`cl_Topology.{hpp,cpp}`): `collect_enrichment_ids()` builds
exclusion sets from `mMesh->thin_shells()` (all shell blocks incl. buffer + tape/ghost sideset
ids); the scan/select loops skip them (`key_exists` guards, empty sets on fresh → zero change);
`detect_sideset_types` deliberately NOT called (types are final on a reloaded mesh; cuts keep
`Default` like fresh). φ rule confirmed pure function of block domain types
(`select_blocks:382-394`: φ=Air∨Buffer∨Ferro, non-φ=Conductor∨Coil). Load branch calls the new
variant (`cl_MaxwellFactory.cpp:442`); `synchronize_maps()` unchanged. Verified no consumer of the
`Cut`/`Default`/`Ghost`/`GeometryOnly` buckets exists.

**Result: `.bfm` reload passes serial (incl. postprocessor), `mpirun -n 2`, `-n 4`, `-n 8`.**

Parked for Christian: the buffer node-flagging loop (`cl_MaxwellFactory.cpp:1813-1819`) claims
buffers need flagging for hanging-node detection but has always been dead on fresh (map predates
the TSF) — should it be live? (D25 keeps it dead on both paths for equality.)

## Verification artifacts

- h5py inspection of `tmp/examples/CORC/corc.bfm` (counts, positions, byte-compare of the
  duplicated blocks) — reproducible one-liner in the D24 entry's history (this session).
- Sideset inventory: shells 1458/1458/1560 (= 4476, `Inactive`), antisymmetry 824/726/720,
  cuts 1862/2303/1374, tape 4476 (`GeometryOnly`, hidden).
