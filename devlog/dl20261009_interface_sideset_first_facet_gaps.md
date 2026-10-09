# Interface sidesets: nothing the processor does depends on the order of facets within a sideset — the two first-facet gaps closed by refusal and by keying on all facets, the numbering claim narrowed to the interface copies

**Date:** 2026-10-09
**Purpose:** Christian's rulings A–C of 2026-10-09 on the blind audit of cfaf7e8, relayed by
the Server session. (C) Close the two first-facet gaps reported in the previous plan's §5
before 0.9.3: a ferro- or buffer-master sideset spanning several blocks relinked only the
first facet's block, and a coil that was only the slave of a later facet left the set
tie-weight. (A) The partly-interface sideset stays refused; the docs say so. (B) The
numbering claim is narrowed to the interface copies, since the cut processor still appends
cut copies in hash order. Plus the reconciled docs pass. One commit on devel.
**Module:** `src/homology` (`cl_InterfaceProcessor`), `tests/homology`, `src/homology/doc`,
`doc/input_file_reference.md`, `CHANGELOG.md`
**AIs:** Claude Fable 5.1 (DrClean session: plan, implementation, gates), Codex
`gpt-6-astra`/high and Grok `grok-4.7`/high (plan jury and code jury), Codex
`gpt-5.6-luna`/medium (language sweep); the Server session (brief, blind xhigh delta audit
and pipeline on the pushed SHA); the Yufan session (iron/coil decks on the SHA)
**Status:** COMPLETE — pushed to devel. Gates: private debug tree (SCLS `debug`, gcc,
Armadillo, non-MKL): build zero warnings, `test_homology` 62/62 (57 before), `make check`
17/17, on the diff as audited and after the fix round; private release tree (SCLS `gcc`,
non-MKL): all 21 shipped decks through setup with no `.bfm` present on the final binary — the
same 20 of 21 as on the f8ce949 baseline (costheta, dipole, gantry and undulator2d, the iron
and coil decks, among them); `corc_twolayer` fails on both in `Cohomology::clean_spfa`
before the processor runs. Positive controls red on f8ce949: the two gap refusals, and the
two mixed-facet cases in the swapped order (dropped at admission there). No deck trips a new
refusal. Yufan's iron/coil comparison on the SHA is hers to report
**Exchange:** `tmp/ai_exchange/review_interface_first_facet_gaps.md`,
`tmp/ai_exchange/sweep_interface_gaps_docs.md`
**Plan:** `todo/closed/interface_sideset_first_facet_gaps.md`

## Summary

- **Principle.** After this commit no outcome of the interface processor depends on the order
  of facets within a sideset. `create_interface_sets` admits a sideset when any of its facets
  is duplicated at, and keys it by the master blocks and the slave blocks of all its
  duplicated facets (`hex(masters)|hex(slaves)`); the master block handed to the set for the
  relink comes from a duplicated facet, never from `facets()(0)`, and `collect_nodes_and_elements`
  no longer reads the first facet (`mMasterBlocks`, a `Cell< id_t >`, `unique`). The set still
  receives the union of both sides for its coil test.
- **Gap (1), refused.** Gathering every facet's master block and relinking them all was the
  plan; the plan jury showed it is not clearly right — `relink_elements` replaces only flagged
  originals, so with a merged ferro|air sideset and a ferro|ferro sideset on the same corners
  the first set in sorted-key order wins, and whether the material split survives would
  follow the block registration order. Per Christian's instruction, a sideset whose duplicated
  facets have non-conductor masters (private copies, no sweep) must lie on one master block;
  otherwise `BELFEM_ERROR` names the sideset, a facet with its element on each of two blocks,
  and the remedy (one sideset per block). Conductor masters stay exempt: their copies are
  shared and the sweep installs them in every conductor block.
- **Gap (2), refused.** A sideset whose slave-bearing facets do not all agree on touching a
  coil block (master or slave is Coil) is refused with a facet of each kind and the remedy
  (the facets touching the coil form a sideset of their own). With the coil facet first the
  topology's "a conductor and a coil must not touch" stops the run before the processor; with
  the air facet first the processor's message does; either way it stops.
- **Interface facets only.** The code jury showed that admitting a sideset that mixes
  duplicated facets with a boundary facet or a conductor|conductor face is unsafe downstream:
  the topology classifies the sideset by its first facet (`ConductorAntiSymmetry` or
  `Inactive`), the conductor|air coupling branch dereferences `slave()` on every facet, and
  the processor duplicates the nodes of the non-interface facet. Such a sideset is now
  refused, in either order, with a facet of each kind. The air|air mixture keeps its
  master-type message (ruling A).
- **Tests** (`test_InterfaceProcessor.cpp`, five new layouts, a facet-order flag, a
  `boundary_facet` helper): `MultiBlockFerroSideSetIsRefused`, `CoilSlaveOnLaterFacetIsRefused`
  (both orders, both messages), `ConductorFacetInInterfaceSideSetIsRefused` and
  `BoundaryFacetInInterfaceSideSetIsRefused` (both orders), `OverlappingFerroSideSetsStayApart`
  (two ferro sidesets with swapped masters on their common face are keyed apart and the face
  stays split — the regression for the master-aware key). The numbering test's comment says
  that reversing the node order changes insertion order only.
- **Docs.** Guide item 4: the rules as they stand (interface facets only, one master type,
  one coil treatment, one block for private copies, admission and key from all facets), the
  real periodic backup order (within-set pairs set by set, then the cross-set pairs), "every
  sideset present before cut processing", and the cut-copy residue; the `topology` section
  of the input reference: boundary facets, the partly-interface surface, the `.bfm` caveat;
  CHANGELOG: both facets and master elements named, the refusal bullet extended, the
  `### Changed` note narrowed to the interface copies.

## Key Findings

- **Ruling B, recorded, no register row.** `CutProcessor::collect_duplicates`
  (`cl_CutProcessor.cpp:1184-1230`) still appends each cut set's copies in the hash order of
  the set's duplicate map; those nodes precede the interface copies, and edge keys use node
  indices, so the edge ids and the dof numbering of a model with cuts still depend on the
  platform's hash table. Remaining work after 0.9.3.
- **Plan jury.** Both auditors: the gather-all relink is order-dependent through the sorted
  set keys (withdrawn, refused instead); admission by `facets()(0)` still dropped a sideset
  whose first facet is a boundary or conductor|conductor facet (admission from all facets);
  Grok: nine stale line ranges, `Vector` vs `Cell`, the coil-first order aborts in Topology,
  the docs in the working tree were ahead of the code.
- **Code jury.** Both: the newly accepted mixtures are unsafe downstream (refused); Codex: the
  union key merged two private sidesets with different masters (key = masters|slaves, with a
  regression); both: the count invariant was a debug-only assert (`BELFEM_ERROR`); Grok: the
  union key separates overlapping private sidesets of one master, leaving a second copy
  unused — pre-existing in kind (`FerroFerro` already leaves the master's air copy unused),
  now order-free, reported.
- **What still reads `facets()(0)`.** Nothing in `cl_InterfaceProcessor.cpp`.
  `Topology::detect_sideset_types` (`cl_Topology.cpp:229-237`) still classifies a sideset by
  its first facet: a boundary-first sideset with a later conductor|coil facet is typed
  `ConductorAntiSymmetry` and refused by nobody (it has no duplicated facet, so this
  processor leaves it alone). Pre-existing, other file, for Christian's list. Also
  pre-existing: `duplicates_at` is false for buffer|coil while the topology maps buffer to
  air, so a buffer|coil sideset is never seen by the processor's checks. For a conductor
  sideset over several conductor blocks the set-local element list follows the first
  duplicated facet's block; the sweep makes every conductor element identical, so no outcome
  changes.
- **Shared checkout.** While this session worked, another session edited the closed
  cohomology core in place (`cl_Chain`, `cl_Cochain`, `cl_Cohomology`, `cl_SimplicialComplex`,
  `cl_Homology`, `cl_BeltedTree`, `cl_CutData`, `cl_CutFactory`); none of it is in this
  commit (explicit staging), and the gates ran on git worktrees at f8ce949 plus this diff.

## Changes Made

- `src/homology/cl_InterfaceProcessor.hpp` — `mMasterBlocks`, `add_master_block`,
  `mMasterBlockOfSideSet`, `touches_coil`, `refuse_mixed_sidesets` (renamed from
  `refuse_mixed_master_types`), contract comments
- `src/homology/cl_InterfaceProcessor.cpp` — admission and key from all duplicated facets,
  the handed master block, the four-part refusal, `BELFEM_ERROR` on the count invariant
- `tests/homology/test_InterfaceProcessor.cpp` — five layouts, the facet-order flag,
  `boundary_facet`, five tests, one comment
- `src/homology/doc/homology_usage_guide.md`, `doc/input_file_reference.md`, `CHANGELOG.md`
- `todo/closed/interface_sideset_first_facet_gaps.md`, `todo/README.md`, this devlog,
  `devlog/README.md`

## Open Questions

- Christian's list: the Topology first-facet classification; buffer|coil unseen by the
  processor; overlapping private sidesets of one master leave an unused copy; the cut-copy
  hash order (ruling B residue); `corc_twolayer` failing setup in `Cohomology::clean_spfa`.
