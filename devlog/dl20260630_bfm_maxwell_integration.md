# BfmFile ↔ MaxwellFactory integration — audit + fix (D13–D21)

**Date:** 2026-06-30
**Purpose:** Audit the newly-wired `.bfm` mesh save/load as driven by `hphirun`/`MaxwellFactory`, fix the
blocking defects, and lock the restart design decision.
**Module:** `src/mesh` (`cl_Mesh_BfmFile`, `cl_Mesh`), `src/fem/maxwell` (`cl_MaxwellFactory`)
**AIs:** Claude (primary), Grok + Codex (independent auditors). Tri-AI on three threads
(`thinshell_bfm`, `vertex_curve_bfm`, `maxwell_loadsave`), all swept into this entry.

## Scope

Three audit passes, then a fix sprint with Christian applying changes and Claude re-auditing each:
1. **Thin-shell save/load** (`BfmFile::save/load_thinshell_data`) — material persistence + node-index
   provenance.
2. **Vertex + curve/segment save/load** — new `BfmFile` entity paths.
3. **MaxwellFactory integration** (`read_mesh`, `create_magnetic_kernel`, material assignment) — the
   actual restart wiring.

## Restart design decision (D21, Christian)

**The input file remains the source of truth.** The `.bfm` exists *only* to cache the expensive
enrichment (cohomology/cuts, thin shells, hanging entities, periodicity) so a rerun can skip
recomputing it. It does **not** own physics assignment:
- Regular block→material is rebuilt from `input.conf` every run.
- Thin-shell *layer* materials are the one exception — they ride in the `.bfm` per-shell, because the
  input-domain pass skips `DomainType::ThinShell`; they're hydrated into `mMaterialBlockAssignment` on load.
- Physical fields + time cursor live in a **separate** restart file (`Mesh::save_fields`/`load_fields`),
  not the `.bfm`.

## Defects found and fixed (D13–D21)

Tri-AI; each verified against the tree. All resolved 2026-06-30 unless noted.

- **D13 (compile blocker)** `Cell<string> tLabels( n,  )` in `save_curve_data` — empty ctor arg. → reserve
  `(n)` + `push` (with `tTypes`/`tClosed`); compiles.
- **D14** `BfmFile` ctor auto-`load()`/`save()` collided with explicit call-site calls → double-execution
  + `checksum()` null-deref. → ctor made passive; `checksum()` made self-contained (own `HDF5` open/read/
  close); `get()` sets `mOwnMesh=false` so returning a stack `BfmFile`'s mesh is safe.
- **D15** `save_group_data` dropped `tIDs.push(tData.mID)` (replaced by an id-guard) → empty `"ids"` → zero
  blocks/sidesets on load. → re-added.
- **D16** `read_mesh` fell off a non-void function when input is non-`.bfm` and no sidecar exists. →
  restructured: source mesh loaded+scaled outside `if(exists)`; all four paths return.
  **D16b:** base checksum cached *before* enrichment (`aMesh->checksum()` moved above `if(exists)`), so the
  sidecar gate matches on later runs.
- **D17** regular block materials unassigned on load → `assign_materials` `Map::operator()` "Key not found"
  crash. → removed the internal `mComputeCohomologies` gate in `create_block_to_material_map` and call it on
  the load path (`else` in `create_magnetic_kernel`); labels come from input domains (load-safe).
- **D18** `create_thinshells` not load-aware → duplicate shells/blocks, or (if input omits shells)
  thin-shell materials missing. → early-return when `thin_shells().size()!=0`, and **hydrate**
  `mMaterialBlockAssignment` from the loaded `ThinShell::materials()` first (fill-if-missing). `comm_barrier`
  symmetry preserved. **D18b:** decoupled `++b` in the create loop (off-by-one if a layer block were
  pre-populated).
- **D19** curves duplicated on load (ctor `create_curves` after `load_curve_data`). → guarded with
  `if ( curves().size()==0 )`; bonus parallel `if ( ! has_periodicity() )` guard on `create_periodic`.
- **D20** *retracted (false positive).* The auditors claimed the saved checksum was the enriched mesh's;
  they missed that `Mesh::checksum()` **caches** (`mHash.value()==0` guard) and is called on the base mesh
  before enrichment, and the only `mHash.reset()` is inside `compute_checksum`. With D16b the base value is
  canonical. Christian caught it.
- **D21** *design decided* (above).

## Fields restart mechanism (R6, partial)

`Mesh::save_fields`/`load_fields` (`cl_Mesh.cpp:2946-3018`): a standalone HDF5 — `meta`:{`timestep`,
`timestamp`,`checksum`}, `fields`:{per-label `data`, `labels[]`, `types[]`}. Load verifies
`checksum == mesh.checksum()` (mesh-identity guard; consistent now that the base checksum is canonical).
Rank-0 only. **Still missing:** global variables, and the controller call-site wiring at the restart cursor.

## Cross-check notes

- Two stale-state corrections went *to* Claude: (1) Claude relayed a removed `BfmFile` block-material map
  from session memory — Grok refuted, verified removed; (2) D20 — Christian was right, both auditors wrong.
  Re-read beats memory.
- Codex's standout catches (missed by Claude + Grok): the `Cell(n)`+index-assign reserve misuse in
  `save_curve_data`, the ctor double-load/save (D14), and the curve duplication (D19).

## State

D13–D21 resolved; the `.bfm` save→load path is structurally complete and the materials/periodicity/
thin-shell/curve load-vs-fresh divergences are closed. **Not yet exercised by a build + round-trip** —
that (R12: save → relaunch → solve, compare to uninterrupted) is the next step, plus globals + the
field-restart call-site wiring (R10/R11). Watch `cl_ProtoMesh.cpp:1474/1480` (`create_periodicitiy`
signature) on the next build — flagged by clangd, likely noise but the kind that can be real.
