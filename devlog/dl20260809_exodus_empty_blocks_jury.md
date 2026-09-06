# Devlog 2026-08-09 — Exodus Empty-Block / ThinShellFactory ID-Clash Jury Investigation

**Date:** 2026-08-09
**Topic:** Reported empty `shell_*_ybco` blocks in hphiTrun exodus output; suspected ThinShellFactory block-ID clash
**AIs involved:** Claude (pre-registration + verification), Codex + Grok (blind jury)
**Claude Confidence:** high (mechanisms), medium ~65% (which event produced the observation)
**Codex Audit Confidence:** high
**Grok Audit Confidence:** high (code paths), medium (event reconstruction)
**Literature References:** none required (I/O bookkeeping, not formulation)
**Verification:** numeric probe — ncdump of every `hphi_results.e-s.*` / `.bfm` in `cmake-build-debug/{sidecoatings,greg2,greg3,greg4,tapestack,corc}` at sideconnectors @ 77011047, plus source trace; no executable gate run ("reviewed + probed", not end-to-end verified)

## Summary

Christian reported hphiTrun writing more thin-shell blocks than the input defines
(`shell_18_ybco`…`shell_21_ybco` present, only block 23 carrying elements) and
suspected an ID clash in the ThinShellFactory. A frozen jury round
(pre-registration → blind Codex+Grok → citation verification → reconciliation) found
**no ID clash and no artifact on disk matching the report** — every existing exodus
file is exactly consistent with its input deck — but confirmed one latent
ExodusWriter defect and three hardening gaps, and converged on a stale-file /
conflated-run explanation for the observation.

## Key Findings

- **ID-clash hypothesis REFUTED (3/3):** ThinShellFactory draws block AND sideset IDs
  from one counter seeded with `Mesh::max_block_and_sideset_id()` (blocks + sidesets +
  curves + thin shells, cl_Mesh.cpp:2233-2265; seed at cl_ThinShellFactory.cpp:40).
  The "missing" IDs in exodus output (e.g. 18, 26, 28, 30 in sidecoatings) are the
  hidden geometry/ghost/wall sidesets sharing the counter — expected holes, not lost
  or clashing blocks.
- **Multiple `shell_NN_ybco` blocks are legitimate:** one layer block per `layers`
  entry; greg2/greg3 decks define `ybco : 0.25 mum` four times → four populated ybco
  blocks (200 elements each, ncdump-verified). Renaming to `<domain>_<id>_<material>`
  happens in `assign_materials()` (cl_MaxwellFactory.cpp:2493-2497).
- **Most plausible cause of the observation (mechanism confirmed, event unverifiable):**
  the `.e-s` series counter restarts at `.00001` every run (cl_Mesh.cpp:378-398) and
  `Mesh::save()` never removes stale higher-numbered siblings; a rerun after a
  layer-config change leaves old-layout files in the series and a viewer shows their
  blocks as empty ghosts next to the new run's populated blocks. The reported IDs mix
  two known layouts (18–21 = greg 4×ybco; 23 = sidecoatings single ybco). The
  offending file was deleted before the session (dir mtime evidence).
- **P1 latent writer defect (3/3):** `ex_put_init` declares ALL mesh blocks
  (cl_Mesh_ExodusWriter.cpp:129,148) while `populate_blocks()` defines/names only
  non-empty ones (:278, :334-337) and connectivity + field truth table iterate the
  unfiltered list (:399-431, :890-904) — a capacity-0 block on the mesh yields an
  internally inconsistent exodus file. Not triggered by any current production path.
- **P2 hardening gaps (single-raiser, Claude-verified, need Christian's adjudication):**
  `Block::element_type()` empty guard is ASSERT-only in the once-per-run writer path
  (cl_Block.hpp:133; frequency rule says BELFEM_ERROR);
  `Block::number_of_elements()` returns pre-sized capacity, not fill count
  (cl_Block.cpp:22-24 `set_size(aNumElements, nullptr)`) — an under-filled block passes
  the writer's >0 filter and nullptr-derefs; `Mesh::add_block`/`update_block_map`
  silently overwrite on duplicate block IDs (cl_Mesh.cpp:2270-2274, :2433-2441);
  `.bfm` reload after a layer-config change zips new materials against old blocks by
  position (cl_MaxwellFactory.cpp:936-944) — mislabeling hazard.

## Changes Made / Proposed

- **UPDATE same day: cause confirmed and P1 fix executed on Christian's approval.**
  Christian confirmed the observation was ParaView unioning a stale `.e-s` series —
  the stale-file explanation holds; no block-generation bug existed.
- P1 writer fix applied (`cl_Mesh_ExodusWriter.{hpp,cpp}`): new `collect_blocks()`
  gathers the nonempty blocks once per `save()`, id-sorted, into member `mBlocks`;
  header count (`ex_put_init`), element-id map, `ex_put_block`/`ex_put_names`,
  connectivity, and the element-field truth table + `ex_put_var` loops all use that
  one list (previously header/connectivity/truth-table used the unfiltered mesh
  list while definitions were filtered). `mNumElements` is now the sum over the
  written blocks, as exodus requires. Side effect: the element-id map now iterates
  in the same id-sorted order as the block definitions (was container order — a
  latent mismatch on meshes whose block container is not id-sorted).
  Verified: `g++ -fsyntax-only` with the real `libbelfem_mesh.a.dir` flags — green;
  full build/run is Christian's.
- Still proposed only: BELFEM_ERROR guards for `element_type()` on empty blocks and
  duplicate `add_block` IDs; `number_of_elements()` capacity semantics; stale
  `.e-s` sibling cleanup/warning on save (P2s pending adjudication).

## Open Questions

- ~~Which run/viewing session produced the observed mix?~~ RESOLVED same day:
  Christian confirmed it was ParaView reading a stale series.
- Adjudication of the four P2 hardening items above.
- Executable gate for the writer fix: next hphiTrun save producing a valid,
  ParaView-readable `.e-s` file (fix is syntax-checked only).

## Files Updated

- tmp/ai_exchange/review_exodus_empty_blocks.md (full record: pre-registration,
  Codex + Grok audits, verification, reconciliation table)
- tmp/ai_exchange/exodus_empty_blocks_brief.md (neutral jury brief)
- devlog/dl20260809_exodus_empty_blocks_jury.md (this file)
