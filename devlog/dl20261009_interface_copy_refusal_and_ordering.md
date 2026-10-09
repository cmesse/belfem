# Interface copies: a sideset with masters of more than one domain type is refused, and copy ids and positions follow the mesh, not the hash table

**Date:** 2026-10-09
**Purpose:** Christian's two rulings of 2026-10-09 on the open DR-198 audit findings, relayed
by the Server session: refuse an interface sideset whose facets do not all have masters of one
domain type, before any copy is allocated; and number the interface copies by node id so that
their ids and mesh positions, and with them the edge ids and the dof numbering, are identical
on every platform. One commit on devel.
**Module:** `src/homology` (`cl_InterfaceProcessor`), `tests/homology`, `src/homology/doc`,
`doc/input_file_reference.md`, `CHANGELOG.md`
**AIs:** Claude Fable 5.1 (DrClean session: plan, implementation, gates), Codex
`gpt-6-astra`/high and Grok `grok-4.7`/high (plan jury and code jury), Codex
`gpt-5.6-luna`/medium (language sweep); the Server session (brief, blind xhigh delta audit
and pipeline on the pushed SHA); the Yufan session (eight-deck comparison against 12a8eec)
**Status:** COMPLETE — pushed to devel. Gates: private debug tree (SCLS `debug`, gcc,
Armadillo, non-MKL): build zero warnings, `test_homology` 57/57 (53 before), `make check`
17/17, on the diff as audited and again after the fix round; private release tree (SCLS
`gcc`, non-MKL): all 21 shipped decks run through setup with no `.bfm` present, 20 of 21
reach the dof banner on the final binary — the same 20 as on the c31ea43 baseline binary;
`corc_twolayer` fails on both at `Cohomology::clean_spfa` (`cl_Cohomology.cpp:534`, the
closed core) before the interface processor runs, so it says nothing about this change.
Positive control: the three new tests red against the c31ea43 source. Yufan's deck comparison
on the pushed SHA is hers to report
**Exchange:** `tmp/ai_exchange/review_interface_copy_ordering.md`,
`tmp/ai_exchange/sweep_interface_copy_docs.md`
**Plan:** `todo/closed/interface_copy_refusal_and_ordering.md`

## Summary

- **Ruling 1, the refusal.** `InterfaceProcessor::refuse_mixed_master_types`, the first
  statement of `create_interface_sets`, walks every original sideset. A sideset is an
  interface sideset when at least one of its facets is one this processor duplicates at
  (`duplicates_at`: a slave exists, the two sides are not both air, one side is air or ferro —
  the admission predicate, now per facet); every slave-bearing facet of such a sideset must
  then have a master block of one domain type, or an always-active `BELFEM_ERROR` names the
  sideset, the two facets with their master elements, both domain types and the remedy. The
  domain-type label is `to_string` plus the enumerator value, so the two coatings, which
  `to_string` both print as "Coating", never read as one type; `en_DomainType.cpp` is
  untouched. A sideset without a duplicated facet (an air|air plus a conductor|conductor facet)
  is left alone, as before. The check precedes the first-facet admission, so a sideset whose
  first facet is air|air and whose second is conductor|air is refused instead of silently
  dropped, in either order. Nothing is allocated before it; the two block bitsets became value
  members so a throwing refusal leaks nothing.
- **Ruling 2, the numbering.** Copies are allocated set by set, in the sets' order (the sorted
  hex strings of the block-index bitsets, a fixed function of the mesh's block order), and
  within a set in ascending id of the original, a shared conductor copy by the first set that
  reaches the original; they are appended to the mesh in ascending copy id, which is the
  allocation order. The three loops that reached mesh state through a `Map` — the allocation
  and the within-set periodic loop in `InterfaceSet::duplicate_nodes`, and the cross-set
  periodic pass — walk `sorted_keys` (a `Cell< id_t >` sorted with `belfem::sort`); the
  collected copies are sorted by id in `add_duplicate_nodes_to_mesh` before the sourcing loop,
  `unify_duplicates` and the append. Left as they were, and why: `relink_elements` and the
  conductor sweep (map lookups inside Cell walks), `collect_nodes_and_elements` (fills a lookup
  map from the mesh's node order), `create_interface_sets` (`unique` already sorts the set
  keys; the other maps are lookups), `add_duplicates` (fills the Cell that is sorted
  afterwards and a lookup map), `unify_duplicates` (walks Cells).
- **Tests.** Three new layouts (`MixedSideSet`, `MixedAirFirst`, `MixedAirLast`,
  `InactiveSideSet`) and a reverse-node option on the fixture; four new tests: the mixed
  sideset refused by name with no node appended; the air|air-plus-conductor|air sideset
  refused in either facet order; the inactive sideset left alone; and the numbering, forward
  and reversed node order, with the absolute oracle (SideSetPerBlock: n 9, m 10, a1 11, a2 12
  at positions 8–11; Ferro: 9–15 in set order "3", "5", "6"), every held copy at position
  id − 1, and the periodic backup pinned to masters 9, 11 and slaves 10, 12 on both builds.
  The existing multi-block same-type sideset test (`OneSideSetRelinksEveryMasterBlock`) is the
  accepted case.
- **Docs.** Guide item 4 of the InterfaceProcessor pipeline (the rule and the numbering), the
  `topology` section of the input reference (one master domain type per interface sideset),
  CHANGELOG `### Refusals` and `### Changed` (ids, positions, edge ids and dof numbering
  differ from 0.9.2 and are now the same on every platform; a `.bfm` from an earlier build is
  still reused, because its key, the base-mesh checksum plus the topology tag, does not
  change, and then replays the old numbering and skips the refusal — delete it).

## Key Findings

- **Plan jury.** Both auditors: a check inside the interface sets misses a sideset whose first
  facet is air|air or a boundary facet, because admission drops it by its first facet — the
  check moved to a processor-level pass over all sidesets. Both: "a function of node ids
  alone" overstated — the contract now names the block order and the set order. Both:
  `Cell< id_t >`, not `Vector`, for the key list. Grok: nine stale line ranges in the plan
  (from a read before 00b7972 moved the file), `%lu` casts, the Ferro oracle (sets "3", "5",
  "6"), the positive control as a run; Codex: compare the restored periodic pairs between the
  two builds; a reused `.bfm` bypasses the change.
- **Code jury.** Both: the first version of the guard also refused sidesets the processor
  never duplicates at — fixed by the per-facet `duplicates_at` predicate (Christian's words
  "the facets of an interface sideset"). Both: a throwing refusal leaked the two bitsets —
  value members now. Grok: the header gloss "i.e. the mesh's block order" was false (the
  sorted hex of masks {1,2} and {0,3} is "6" before "9"), one position assertion was vacuous,
  and the periodic section pinned nothing — all fixed.
- **Positive control.** On the c31ea43 source the forward build numbers the SideSetPerBlock
  copies 11, 9, 10, 12 (hash order), and the reversed build differs from it; both refusal tests
  red (four copies appended).
- **Not refused, for Christian's ruling (plan §5).** Facets of one sideset that share the master
  type but differ on the slave side or in block pair pass the check. Today the set is bucketed
  by the first facet's block pair and relinks only the first facet's master block: a
  multi-block conductor sideset works because the DR-198 conductor sweep installs the shared
  copy everywhere; a multi-block ferro or buffer sideset (private copies, no sweep) leaves the
  other blocks' elements on the original; a coil slave on a later facet does not decouple the
  set; a boundary facet first still drops the sideset at admission. None of this is reached by
  a shipped deck as far as the setup gate shows.
- **Gate baseline.** Running the decks required building the six plugins
  (`cmake -DBELFEM_DIR=… src && make`) and regenerating the `pancake` and `rutherford` meshes
  with their shipped generators (the copies in the shared build trees belonged to other
  decks). `corc_twolayer`, freshly meshed by its generator, fails setup on the unchanged tip
  in `Cohomology::clean_spfa`; not a finding of this change, reported to the Server.

## Changes Made

- `src/homology/cl_InterfaceProcessor.hpp` — value-member bitsets, `duplicates_at`,
  `refuse_mixed_master_types`, the numbering contract
- `src/homology/cl_InterfaceProcessor.cpp` — the refusal, `domain_type_label`, `sorted_keys`,
  the three sorted loops, the sort before append
- `tests/homology/test_InterfaceProcessor.cpp` — four layouts, the reverse-node option,
  `index_in_mesh`, four tests
- `src/homology/doc/homology_usage_guide.md`, `doc/input_file_reference.md`, `CHANGELOG.md`
- `todo/closed/interface_copy_refusal_and_ordering.md`, `todo/README.md`, this devlog,
  `devlog/README.md`

## Open Questions

- Plan §5: the same-master mixed slave/block-pair sideset, the coil-on-a-later-facet case and
  the boundary-first drop — Christian's call whether any becomes a refusal.
- `corc_twolayer` fails setup on devel at `Cohomology::clean_spfa` with the mesh its shipped
  generator writes today; Gregory's module.
