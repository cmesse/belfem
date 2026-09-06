# Devlog 2026-06-22 — BFM Mesh-File Save/Load Refactor (planning)

**Date:** 2026-06-22
**Topic:** Planning-only investigation + plan for the `.bfm` (HDF5) mesh save/load refactor around `ProtoMesh`/`Protoshell`, enabling restart-from-last-timestep with a fully enriched mesh.
**AIs involved:** Claude (exploration + plan), Codex (audit, high), Grok (third-voice audit, high)
**Claude Confidence:** medium→high (gap classification verified against code)
**Codex Audit Confidence:** high
**Literature References:** N/A (I/O / state-serialization task; not algorithmic)

## Summary
Autonomous overnight planning run. Produced `todo/meshfile_refactor_plan.md`: a self-contained plan
to repair and refactor the broken `.bfm` mesh save/load so a Maxwell run can stop, change parameters,
and relaunch from the last stored timestep, loading the fully-enriched post-factory mesh and skipping
the cohomology/cut/thin-shell factories. No source modified (planning only).

## Key Findings (all citations re-verified before inclusion)
- **The enriched BFM write is disabled**, not merely buggy: `cl_MaxwellFactory.cpp:469-478` comments
  out the `HDF5Writer + save_material_map + save_thin_shell_facets` block and prints a `#warning`.
- **The reload path is half-built and inconsistent:** `read_mesh()` already prefers a checksum-matched
  `.bfm` and sets `mComputeCohomologies=false` (`:227-313`), and `load_domain_types_and_material_map`
  (`:2124`) runs on reload, but the materials it reads are only written by the disabled path, and
  `save_thin_shell_facets` (`:2173`) writes a `tapes` group nothing ever reads.
- **The reconstruction engine already exists:** `ProtoMesh` + ID-keyed proto structs (NodeData,
  ElementData, ThinShellData, PeriodicityData, TMatrixData, …) rebuild a live mesh; today only the MPI
  `Distributor` drives them. The refactor adds an HDF5 backend for the same structs.
- **Element type is NOT deterministically rebuildable** (Codex, verified): `element_type_from_numnodes`
  hard-errors on 20-node 3D elements (TET20/HEX20 collision, `meshtools.cpp:1181-1184`) and the
  thin-shell switch only handles 4/9/6/18 nodes, so HEX8TS is unrecoverable. → must-save element type.
- **Periodicity reload is a HARD ERROR, not silent loss** (Grok, verified): `from_proto` sets
  `mNodePairsRestored=true` and never sets `periodic()` pointers (`PeriodicityFactory.cpp:320`); a later
  `update()` with no node-pair backup hits `BELFEM_ERROR` at `:359-370`. `to_proto` also never fills
  `mMaster/SlaveSideSets`. The reload needs an explicit periodic-finalize contract (plan step R7).
- **Index-vs-ID keying is the central risk:** duplicates, facets (master/slave by element index),
  edges, hanging, abstract/orphan all round-trip only because serial reload preserves array order. The
  proto structs are ID-keyed; the refactor must adopt ID keying throughout.
- Two pre-existing reader bugs noted for early fix: `read_faces_2d` missing `tFaces.set_size()`
  (`HDF5Reader.cpp:811-836`); leftover debug `std::cout` in `link_elements_to_edges` (`:774,:796`).

## Changes Made / Proposed
- **Created** `todo/meshfile_refactor_plan.md` — failure analysis, 23-row gap table (a/b/c + one
  citation each), 13-step ordered todo (R0-R12 with dependencies), 9 open design questions, and a
  proposed ID-keyed HDF5 schema. Updated `todo/README.md` index.
- **Merged + deleted `todo/hdf5_writer_repair.md`** (the older 2026-05-01 Claude+Codex audit of the
  current `HDF5Writer`/`HDF5Reader`). It overlapped with this task but was complementary: a low-level
  bug catalogue, not the architectural refactor. Folded into the plan as **Appendix A** (W1-W6, H1-H3,
  M1-M2) with all line numbers **re-verified against the current tree** (the 2026-05-01 numbers had
  drifted), the disable-commit provenance (`6c52293`, 2025-08-26) added to §1.1, and the re-enable
  feature-flag question added as O9. R0 now consumes the catalogue, split into must-fix-regardless
  (W6, shared `HDF5_Tools.cpp`), fix-now legacy-reader bugs (W4/W5/H1/H2), and implementation cautions
  for the new writer (W1/W2/W3). Verified W1-W6 still present before adopting.
- No source code changed (planning only, per task scope).

## Open Questions (carried into the plan)
- O1/O2: exact periodic reload API — set `periodic()` pointers in `from_proto` vs repopulate node-pair
  backup; confirm no consumer still needs a geometric `update()`.
- O3: does restart need IWG BDF history (`cl_IWG_Timestep.cpp:258-285`), or are in-mesh fields enough
  for first-order restart? (Likely the deciding correctness question.)
- O5/O8: control-point meshes; curve-coordinate recomputation on cut/thin-shell-bounded curves.
- O7: retire vs keep the legacy `HDF5Writer`/`HDF5Reader`.

## Files Updated
- todo/meshfile_refactor_plan.md (new; Appendix A merged in from hdf5_writer_repair.md)
- todo/hdf5_writer_repair.md (deleted — merged into the plan)
- todo/README.md
- devlog/README.md
- (exchange: tmp/ai_exchange/meshfile_refactor.md — AI-only, ephemeral; distilled here before sweep)
