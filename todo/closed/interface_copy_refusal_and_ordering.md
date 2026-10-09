# Interface copies: one master domain type per sideset, and id-ordered numbering

**Date:** 2026-10-09
**Purpose:** Two rulings of Christian's quiz of 2026-10-09 on the open DR-198 audit findings. (1) An interface sideset whose facets do not all have masters of one domain type is refused during setup, before any copy is allocated. (2) The ids and the mesh positions of all interface copies become a function of node ids alone, never of hash-map iteration order, so they are identical on every platform. One commit on devel, `src/homology/cl_InterfaceProcessor.{hpp,cpp}` plus tests and docs.
**Module:** `src/homology` (`cl_InterfaceProcessor`), `tests/homology`, `src/homology/doc`, `doc/`, `CHANGELOG.md`
**AIs involved:** Claude (plan, implementation, gates), Codex + Grok (plan jury, code jury); the Server session runs the blind xhigh delta audit and the pipeline on the pushed SHA; the Yufan session reruns the eight DR-198 decks against 12a8eec
**Status:** ✅ COMPLETE 2026-10-09. Landed on devel in one commit: `refuse_mixed_master_types` over every sideset with a duplicated facet (`duplicates_at`, the admission predicate per facet), before admission and before any copy; `sorted_keys` on the allocation, the within-set and the cross-set periodic loops; the copies sorted by id before sourcing, unify and append; the two block bitsets as value members; the numbering contract in the header. Four new tests (two refusals, one acceptance, the numbering forward and reversed with the absolute oracle and the pinned periodic backup), red on c31ea43 where they should be. Guide, input reference §8, CHANGELOG Refusals + Changed (luna sweep). Gates: debug tree zero warnings, `test_homology` 57/57, `make check` 17/17 (twice); release tree, 21 decks through setup with no `.bfm`: 20 of 21 on the final binary, the same 20 as the c31ea43 baseline, `corc_twolayer` failing on both in `Cohomology::clean_spfa` before the processor runs. Juries: plan astra/high + grok-4.7/high (refusal relocated, contract narrowed, `Cell< id_t >`, oracle locked); code astra/high + grok-4.7/high with one fix round (qualifying-facet predicate, bitset leak, header gloss, two test pins). Residue for Christian: §5 (same-master mixed slave/block pair, coil on a later facet, boundary-first drop); Yufan's eight-deck comparison on the pushed SHA.

> **Scope guards (from the brief):**
> - Ruling 1 covers MASTER domain types only. A sideset whose facets share the master type but differ on the slave side or in block pair is reported (§5), not refused.
> - `en_DomainType.cpp` is not touched (the Coating label collision is handled in the message, §2.1).
> - The periodic-refusal message of 00b7972 keeps its text; only new text is added.
> - devel only (tip `c31ea43`); never main. No register rows.

---

## 1. Current Behaviour and How It Fails

| Failure | Mechanism | Evidence |
|---|---|---|
| A mixed sideset is classified by its first facet | `InterfaceSet::collect_nodes_and_elements` takes `mMasterType` from `mSideSets(0)->facets()(0)` and one master block per sideset from its first facet (`cl_InterfaceProcessor.cpp:116-128`); `create_interface_sets` admits or drops a sideset, and buckets it, by its first facet's block pair (`:363-389`). `fix_facet_masters` flips a facet whose master has the lower domain type (`cl_MaxwellFactory.cpp:1471-1482`); Conductor and Ferro both outrank Air, so a sideset grouping a conductor\|air and a ferro\|air facet keeps both masters — a user grouping the walk never splits. The whole set is then shared or private, paired or refused across a periodic plane, by facet order; a sideset whose first facet is air\|air or a boundary facet is dropped without a look at the rest (`:365-379`). | source trace |
| Copy ids follow hash order | `InterfaceSet::duplicate_nodes` allocates `new Node( ++aMaxNodeID, … )` (`:192`) while iterating `mOriginals` (`:183-209`), a `belfem::Map` = `std::unordered_map` (`cl_Map.hpp:15-17`); libstdc++ iteration order depends on insertion history, so the id a copy gets depends on the platform's hash table, not on the mesh. | source trace |
| Copy positions follow hash order | `add_duplicate_nodes_to_mesh` fills `tDuplicates` by iterating `mSharedDuplicates` and then each set's `mDuplicates` through `add_duplicates` (`:615-627`), sources them (`:634-666`), and `append`s that order to `mMesh->nodes()` (`:672`); `unify_duplicates` (`:676-779`) adds the new copies to their originals in that order (`:771-776`). | source trace |
| Periodic backup order follows hash order | the within-set loop (`:211-246`, backup at `:244`) and the cross-set pass (`:493-566`, backup at `:565`) iterate `mOriginals` / `original_map()` and call `add_node_pair_to_backup` in that order. | source trace |

Consequence: edge ids and the dof numbering downstream differ between platforms (the Server's audit caveat on the DR-198 gate: "equal results on this toolchain do not show that the shared-copy append order is platform independent").

**Bottom line:** classify a sideset by all of its facets and refuse a mixed one; hand out copy ids and positions by node id.

## 2. Mechanism

### 2.1 Refusal (R1)

New private `InterfaceProcessor::refuse_mixed_master_types()`, the first statement of `create_interface_sets` — before any sideset is admitted or dropped by its first facet (Codex/Grok plan round 1: a sideset whose first facet is air|air or a boundary facet never reaches a set), and before `mOriginals`, `mElements` or any copy exists (copies are allocated only in `duplicate_nodes`, the sole allocator, which runs after). For every original sideset (`s < mNumOriginalSidesets`), over its facets that have a slave: the master domain type of the first such facet is the reference; every facet whose master block has another domain type raises an always-active

```
BELFEM_ERROR( false,
  "interface sideset %lu mixes master domain types: facet %lu ( element %lu ) has a %s master, facet %lu ( element %lu ) a %s master; split the sideset so that each sideset has one master domain type",
  ( long unsigned int ) sideset id, first facet id, its master element id, label( type0 ), offending facet id, its master element id, label( type ) )
```

Every id cast to `long unsigned int` (`id_t` is `unsigned int`). Facets without a slave (boundary) take no part in the check: they are not interface facets.

`label( DomainType )` is file-local and prints `to_string( type )` followed by the enumerator's integer in parentheses, e.g. `Conductor ( 5 )`, `Coating ( 6 )` against `Coating ( 7 )`: two different types never print as one word, and `en_DomainType.cpp` stays as it is. Setup-time code, so `BELFEM_ERROR` survives a release build.

The ruling covers master types only; §5 records what the same scan sees on the slave side and the block pair.

### 2.2 Ordering (R2) — the rule, stated in the header

> **Copy numbering.** Interface copies are allocated set by set, in the sets' order (sorted hex strings of the block-index bitsets, i.e. the mesh's block order), and within a set in ascending id of the original; a shared conductor copy is allocated by the first set that reaches the original. Copies are appended to the mesh in ascending copy id, which is the allocation order. Both are deterministic from the block order, the set order and the original node ids, and independent of hash-map iteration — hence identical on every platform for the same mesh file.

Implementation, file-local `sorted_keys( const Map< id_t, Node * > &, Cell< id_t > & )` using `belfem::sort( Cell< T > & )` (`cl_Cell.hpp:430-436`); gathered once per set and reused by the allocation and the within-set periodic loop:

| loop | today | change | why |
|---|---|---|---|
| `InterfaceSet::duplicate_nodes`, allocation (`:183-209`, `new Node` at `:192`) | `for ( auto tPair : mOriginals )` | walk `sorted_keys( mOriginals )` | hands out `++aMaxNodeID` |
| `InterfaceSet::duplicate_nodes`, within-set periodic loop (`:211-246`) | same map | the same sorted keys | order of `add_node_pair_to_backup` (`:244`) |
| `InterfaceProcessor::pair_cross_set_periodic_duplicates` (`:493-566`) | `for ( auto tPair : tSet->original_map() )` | sorted keys | order of `add_node_pair_to_backup` (`:565`) |
| `InterfaceProcessor::add_duplicate_nodes_to_mesh` (`:615-672`) | shared map, then each set's map | collect as today, then sort `tDuplicates` by id before the sourcing loop (`:634`), `unify_duplicates` and `append` (`:672`) | mesh position, the `Node::add_duplicate` order in `unify_duplicates` (`:771-776`), source assignment order |
| `InterfaceSet::add_duplicates` (`:280-296`) | iterates `mDuplicates` | unchanged | fills a Cell that is sorted afterwards and a lookup map (`aOriginalMap`) whose iteration order nothing reads |
| `InterfaceSet::relink_elements` (`:250-273`), the conductor sweep (`:577-597`) | map lookups inside element loops | unchanged | `mElements`, blocks and elements are Cells; the maps are only looked up |
| `collect_nodes_and_elements` (`:62-163`) | fills `mOriginals` from `mMesh->nodes()`, `mElements` from block order | unchanged | the map is a lookup; the Cell is deterministic |
| `create_interface_sets` (`:355-410`) | `unique( tStrings )` sorts (`:392`) | unchanged apart from the refusal call | `mSets` already deterministic given the block order; `tMapA`/`tMapB`/`mSetsMap` are lookups |
| `unify_duplicates` (`:676-779`) | iterates `aDuplicates` and `mMesh->nodes()` | unchanged | both Cells; `aDuplicates` arrives sorted |
| sourcing loop in `add_duplicate_nodes_to_mesh` (`:634-666`) | iterates `tDuplicates` | runs after the sort | per-copy, order-insensitive, but now deterministic too |

The set order itself (`unique` on hex strings of the block-index bitset, `mBlockIndices` from `aMesh->blocks()` order, `:313-317`) depends on the mesh's block order, which is file order.

### 2.3 Behaviour change

Copies keep the same count and the same originals; only their ids and their positions in `mMesh->nodes()` change (0.9.2 and every build since hand them out in hash order). Edge ids follow (`edge_key` uses `Node::index()`, `cl_EdgeFactory.cpp:295-305`, ids `++tMaxID` in key order, `:544-569`) and the dof numbering with them; results agreeing to round-off is the hypothesis Yufan's eight-deck comparison on this SHA against 12a8eec tests, not a property shown here. A `.bfm` written before this commit is still reused — the cache key is the base-mesh checksum (`cl_Mesh.cpp:2284-2318`) plus the config tag (`fn_mesh_config_tag.hpp:325-348`), neither of which this commit changes — and a reused cache sets `mComputeCohomologies = false`, which skips `create_cuts` and with it both the refusal and the numbering (`cl_MaxwellFactory.cpp:442-447`, `:673`). Every gate run and Yufan's comparison start without a `.bfm`; the CHANGELOG says so.

## 3. Ordered Steps

- [x] **R1** `cl_InterfaceProcessor.{hpp,cpp}`: private `InterfaceProcessor::refuse_mixed_master_types()` + file-local `label( const DomainType aType )`, first statement of `create_interface_sets`.
- [x] **R2** the ordering: `sorted_keys` into a `Cell< id_t >`, the three sorted loops, the sort of `tDuplicates` in `add_duplicate_nodes_to_mesh` before sourcing, unify and append, the rule in the header (contract wording of §2.2).
- [x] **R3** `tests/homology/test_InterfaceProcessor.cpp`:
  - `Layout::MixedSideSet`: block 2 Ferro, one sideset 11 over `tTop1` (conductor master) and `tTop2` (ferro master). `MixedMasterSideSetIsRefused`: `std::runtime_error` whose message holds `sideset 11`, `Conductor ( 5 )`, `Ferro ( 3 )`, `facet 102`, `element 2` and `one master domain type`; `number_of_new_nodes() == 0` — this measures appended nodes; that nothing was allocated either follows from the check being the first statement of `create_interface_sets`, before `duplicate_nodes`, the sole allocator (stated in the header comment). The id counter is the constructor's by-value argument and cannot be observed.
  - `Layout::MixedAirFirst` / `MixedAirLast`: sideset 11 holds an air|air facet (A1|A2 on the face n, m, t) and the conductor|air facet `tTop1`, in both orders; both refused with `Air ( 1 )` against `Conductor ( 5 )` — the first-facet admission at `:365-379` must not hide the mixed sideset (Codex/Grok).
  - the accepted multi-block same-type sideset is the existing `OneSideSetRelinksEveryMasterBlock` (one sideset over two conductor blocks); it stays green and is named in the test comment.
  - fixture flag `aReverseNodes`: the eight nodes are pushed into `mMesh->nodes()` in reverse slot order (ids unchanged). `CopyNumberingIsIndependentOfNodeOrder`: for `SideSetPerBlock`, `Ferro` and the periodic `SideSetPerBlock`, build forward and reversed; for every slot and both C elements, the held copy's id and its index in `mMesh->nodes()` agree between the two builds; the absolute oracle on the forward build — `SideSetPerBlock`: sets "5" (11: n, m, a1) then "6" (12: a2) → n 9, m 10, a1 11, a2 12, appended at positions 8–11; `Ferro`: sets "3" (13: n, m, q shared conductor), "5" (11: n, m, a1 shared), "6" (12: ferro private n, m, a2) → 9, 10, 11 for n, m, q; 12 for a1; 13, 14, 15 for the ferro copies of n, m, a2; and on the periodic layout the restored `master_nodes()` / `slave_nodes()` id sequences agree between the two builds (Codex). Positive control (G3): the new test TU linked against the `c31ea43` source — the reversed build is expected to differ on libstdc++ (insertion-dependent iteration; a header-level walk gives 3, 1, 2 vs 1, 3, 2 for keys 1–3, Grok), shown by the run, not assumed; the absolute oracle fails there either way if the old walk is stable but unsorted.
  - every existing InterfaceProcessor test stays green (the periodic ones count backup entries, not order).
- [x] **R4** docs: `src/homology/doc/homology_usage_guide.md` (the refusal and the numbering rule in the InterfaceProcessor section); `doc/input_file_reference.md` §8 `topology` (one master domain type per interface sideset); `CHANGELOG.md` `### Refusals` bullet and `### Changed` note (ids and positions of interface copies, and with them edge ids and dof numbering, differ from 0.9.2 and are now independent of the hash table, hence of the platform; results agree to round-off — Yufan's comparison; a sibling `.bfm` replays the old numbering and skips the refusal because its key, base-mesh checksum plus config tag, does not change — delete it). Codex language sweep on the two user-facing passages.
- [x] **R5** gates: (G1) `test_homology` and `make check` in a private debug tree; (G2) all 21 shipped decks through setup in a private release tree with no `.bfm` present (reach the dof banner, then stop the run), against a baseline pass of the same decks on the `c31ea43` binary — a deck that trips the refusal is a finding for Christian, not something to work around; (G3) the numbering test against `c31ea43` (run); (G4) Yufan's eight-deck comparison on the pushed SHA against 12a8eec, fresh `.bfm` (her session).
- [x] **R6** close out: devlog, `devlog/README.md`, this file to `todo/closed/`, push devel, SHA + file list to the Server and to Yufan.

## 4. What could break

- A shipped deck with a mixed sideset: unknown until G2 runs. `pancake_tape` (gone) had a single-type SolderAir sideset over 43 conductor blocks — the multi-block same-type case, which stays accepted. [confidence: medium]
- Periodic decks (helix, corc, rutherford, tapestack_circuit): the backup order changes, the pairs do not; `restore_node_pairs` consumes every entry regardless of order. [high]
- A `.bfm` written before this commit is reused (checksum + config tag unchanged) and then replays the old numbering AND skips the refusal. Delete it before any run that should exercise this commit. Noted in the CHANGELOG. [high]
- Key list: `Cell< id_t >` with `belfem::sort( Cell< T > & )` (`cl_Cell.hpp:430-436`); never `Vector` for an id list. [high]

## 5. Not refused — for Christian's ruling (the slave side and the block pair)

Per the brief, the same scan could see, per sideset, facets that share the master type but have (a) different slave blocks or (b) different master blocks; this plan refuses neither and logs nothing. Today: `create_interface_sets` buckets the sideset by its FIRST facet's block pair (`:381-385`), and `collect_nodes_and_elements` collects the master block of each sideset's FIRST facet only (`:116-126`), so (b) relinks only the first master block's elements through the set — a multi-block conductor sideset works because the DR-198 conductor sweep installs the shared copy in every conductor element (`:577-597`, Conductor only, `:584`); a multi-block FERRO or BUFFER sideset (private copies, `shares_duplicates` is TieWeight1 ∧ Conductor, no sweep) would leave the other blocks' elements on the original, i.e. unsplit. For (a), a facet whose slave is a coil block while the first facet's slave is air would not make the set `Decouple` (`Decouple` is derived from the first pair's bitset, `:47-54`), and there is no coil sweep. Neither is reached by a shipped deck as far as the setup gate shows (G2 says nothing about them either way: it only shows no deck trips the master-type refusal).

## 6. Open Design Questions

- O1 — append order of copies: ascending copy id (= allocation order, shared and private interleaved by set) rather than "shared first". Chosen because it is one rule for both kinds and purely id-based. Reporting choice, stated in the header.

## 7. Audit Trail

- `tmp/ai_exchange/review_interface_copy_ordering.md` — plan and code juries.
