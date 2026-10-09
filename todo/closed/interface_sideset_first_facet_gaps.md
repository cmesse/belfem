# Interface sidesets: no outcome may depend on facet order within a sideset

**Date:** 2026-10-09
**Purpose:** Christian's rulings A–C of 2026-10-09 on the blind audit of cfaf7e8. (C) The two first-facet gaps reported in `todo/closed/interface_copy_refusal_and_ordering.md` §5 are closed before 0.9.3: a ferro- or buffer-master sideset spanning several blocks relinks every master block, not only the first facet's; a sideset whose facets do not all agree on touching a coil is refused. (A) The partly-interface sideset stays refused and the docs say so with the remedy. (B) The numbering claim is narrowed to the interface copies: the cut processor still appends cut copies in hash order, remaining work after 0.9.3. Plus the reconciled docs pass and two cheap acceptance tests. One commit on devel.
**Module:** `src/homology` (`cl_InterfaceProcessor`), `tests/homology`, `src/homology/doc`, `doc/`, `CHANGELOG.md`
**AIs involved:** Claude (plan, implementation, gates), Codex + Grok (plan jury, code jury); the Server session (blind xhigh delta audit and pipeline on the pushed SHA); the Yufan session (iron/coil decks on the SHA)
**Status:** ✅ COMPLETE 2026-10-09. Landed on devel in one commit: admission and set key from all duplicated facets of a sideset (`hex(masters)|hex(slaves)`), the master block handed to the set from a duplicated facet (no `facets()(0)` read left in the file); `refuse_mixed_sidesets` with four checks — one master domain type, one coil treatment, interface facets only, one master block for private copies; the count invariant as `BELFEM_ERROR`. Gap (1) refused rather than handled (plan jury: gather-all relinking is order-dependent through the sorted set keys); gap (2) refused; the mixed-facet sidesets refused too (code jury: unsafe downstream). Five tests, the refusals red on f8ce949. Docs: rulings A and B and the reconciled pass (luna sweep). Gates: debug zero warnings, `test_homology` 62/62, `make check` 17/17; 21 decks through setup on the final release binary, 20/21, the same 20 as the f8ce949 baseline (`corc_twolayer` fails in `clean_spfa` on both, before the processor). Juries: plan astra/high + grok-4.7/high; code astra/high + grok-4.7/high with one fix round. Residue for Christian: Topology's first-facet classification (other file), buffer|coil unseen by the processor, overlapping private sidesets of one master leave an unused copy, the cut-copy hash order (ruling B).

> **Scope guards (from the brief):**
> - Principle: after this commit no outcome of the interface processor may depend on the order of facets within a sideset. Anything that still reads `facets()(0)` is reported (§5), not widened into this commit.
> - `cl_CutProcessor.cpp` is not touched (ruling B: docs only).
> - Ruling A: the partly-interface sideset (conductor|air with air|air) stays refused.
> - devel only (tip `f8ce949`); never main. No register rows.

---

> **REVERTED 2026-10-09 (Christian's ruling):** e30bd98 left devel the same day; the processor is back at cfaf7e8. The work is on `archive/first-facet-sidesets`. See `devlog/dl20261009_first_facet_campaign_reverted.md`.


## 1. Current Behaviour and How It Fails

| Gap | Mechanism | Evidence |
|---|---|---|
| (1) multi-block private-copy sideset | `InterfaceSet::collect_nodes_and_elements` collects the nodes of EVERY facet (`cl_InterfaceProcessor.cpp:110-127`) but the master block of each sideset's FIRST facet only (`:139-151`), so `mElements` holds the elements of that one block; `relink_elements` rewrites only those (`:277-300`). A conductor master is rescued by the DR-198 sweep (`:600-620`, Conductor only); a ferro or buffer master is not: `Ferro F1 \| Air` first and `Ferro F2 \| Air` later leaves F2's elements on the original, unsplit from the air. | source trace, both auditors |
| (2) coil only on a later facet | `InterfaceSet`'s treatment is decided from the bitset of the set, which is the FIRST facet's block pair (`:64-78`, `:404-412`): `Conductor C \| Air A` first and `Conductor C \| Coil K` later stays tie-weight and shared, so the coil face gets a sourced copy instead of a decoupled one. `Topology::detect_sideset_types` classifies by the first facet as well (`cl_Topology.cpp:229-237`), so its "a conductor and a coil must not touch" refusal (`:329`, `:359`) is skipped for the later facet. | source trace, both auditors |

**Bottom line:** refuse a private-copy sideset whose duplicated facets lie on more than one master block (gap (1), since relinking every block is not clearly right — jury round 1); refuse a sideset whose slave-bearing facets do not all agree on touching a coil (gap (2)); and admit and key a sideset from all of its duplicated facets, not the first.

## 2. Mechanism

### 2.1 Gap (1), refused

> **Withdrawn after plan round 1:** gathering the master block of every facet and relinking them all is NOT clearly right. `InterfaceSet::relink_elements` replaces a node only while it is still a flagged ORIGINAL (`:280-293`), so with a merged ferro|air sideset and a ferro|ferro sideset on the same corners the first set in sorted-key order wins, and whether the ferro|ferro split survives depends on the block registration order (Codex and Grok, independently). Christian's instruction for that case: refuse.

A sideset whose duplicated facets have non-conductor masters (ferro, buffer, coil: private copies, no sweep) must have all its duplicated facets on ONE master block; otherwise an always-active `BELFEM_ERROR` names the sideset, one facet with its element on each of two blocks, and the remedy: one sideset per block. Conductor masters are exempt: their copies are shared and the sweep installs them in every conductor block. The set's master block for `collect_nodes_and_elements` is handed over by the processor from the first duplicated facet of each admitted sideset (one block for private sets by this refusal; any one block for conductor sets, the sweep covers the rest), replacing the `facets()(0)` read there; the block list becomes a `Cell< id_t >`.

### 2.1b Admission and set identity from all duplicated facets

`create_interface_sets` admits a sideset when ANY of its facets is duplicated at (today: its first facet, `:388-400`, so a boundary or conductor|conductor facet first dropped the sideset). The set key is the union of the (master, slave) block-index pairs of ALL duplicated facets of the sideset, so neither admission, nor which sidesets share a set, nor the set order (and with it the copy ids) depends on facet order. `InterfaceSet`'s coil test over the bitset is unchanged and, with the coil-agreement refusal, order-independent.

### 2.2 Gap (2), refused

A second check in the same walk as `refuse_mixed_master_types` (renamed `refuse_mixed_sidesets`, since it now checks three things: master type, master block for private copies, coil agreement): for an interface sideset (one with a duplicated facet), over its slave-bearing facets, "touches a coil" = master or slave block is `DomainType::Coil`; the first such facet is the reference; a facet that disagrees raises an always-active

```
BELFEM_ERROR( false,
  "interface sideset %lu mixes facets that touch a coil block with facets that do not: facet %lu ( element %lu ) %s, facet %lu ( element %lu ) %s; split the sideset so that the facets touching a coil form a sideset of their own",
  sideset id, first facet, its master element, "touches coil block N" / "touches no coil", offending facet, its master element, the other label )
```

Ids cast to `long unsigned int`. No attempt to split a set by treatment. After this, the first facet's pair contains a coil block if and only if every facet touches one, so the treatment cannot depend on order.

### 2.3 Rulings A and B, docs

- A: the guide, the input reference and the CHANGELOG say that a sideset mixing duplicated facets with air|air facets is refused, in either order, and that the remedy is to split the surface.
- B: the CHANGELOG `### Changed` note and the guide promise deterministic numbering of the INTERFACE copies only; the cut copies that `CutProcessor::collect_duplicates` appends (`cl_CutProcessor.cpp:1214`) still come in hash order and precede the interface copies, so edge and dof numbering is not yet platform-independent; remaining work after 0.9.3, noted in the devlog, no register row.
- the reconciled pass: the backup order (all within-set pairs, set by set, then the cross-set pairs); "every sideset" → the sidesets present before cut processing; the input reference gains the boundary-facet, partly-interface and `.bfm` sentences; the CHANGELOG refusal bullet names both facets and both master elements.

## 3. Ordered Steps

- [x] **R1** `create_interface_sets`: admission on any duplicated facet, set key = union of all duplicated facets' pairs, master block handed to the set (`InterfaceSet::add_master_block`); `collect_nodes_and_elements` reads the handed blocks (`Cell< id_t >`, `unique`), never `facets()(0)`.
- [x] **R2** `refuse_mixed_sidesets`: master type (as before), private-copy sidesets on one master block, coil agreement; header comments.
- [x] **R3** tests: `MultiBlockFerroSideSetIsRefused` (block 1 and 2 Ferro, ONE sideset over both air facets): refused with the message, nothing appended — red on f8ce949 (no refusal; C2 keeps the original); `CoilSlaveOnLaterFacetIsRefused` (block 2 Coil; sideset 11 = C1|A1 and C1|C2 in both orders: air-first reaches the processor's message, coil-first aborts in `Topology::sideset_type` with "a conductor and a coil must not touch"); accepted cases `ConductorFacetInInterfaceSideSet` (C1|A1 with C1|C2, both orders) and `BoundaryFacetInInterfaceSideSet` (C1|A1 with a slave-less facet of C1, both orders): no throw, n and m get one conductor copy, the shared face conforms — the reversed orders red on f8ce949 (dropped at admission); the numbering test's comment says reversal changes insertion order only.
- [x] **R4** docs per §2.3; Codex language sweep on the user-facing passages.
- [x] **R5** gates: debug tree `test_homology` + `make check`; release tree, all 21 decks through setup with no `.bfm` (costheta, dipole, gantry, undulator2d carry iron or coils and matter most); `corc_twolayer`'s `clean_spfa` failure is known and separate.
- [x] **R6** close out: devlog, index, this file to `todo/closed/`, push devel, SHA + file list to the Server.

## 4. What could break

- A shipped deck with a multi-block ferro or buffer sideset now trips the new refusal: a finding for Christian, not a workaround. Which decks: unknown until R5 (costheta, gantry, undulator2d carry iron). [medium]
- A shipped deck with a coil on a later facet trips the new refusal: a finding for Christian. [medium]
- `mMasterType` from the first collected block: all slave-bearing masters have one type by R1, so no outcome depends on which is first. [high]

## 5. What still reads `facets()(0)`, and whether it can change an outcome

- `create_interface_sets` (`:388-412`): after this commit, admission, the set key and the handed master block come from the set of duplicated facets, not from `facets()(0)`. No read of `facets()(0)` remains in this file.
- `mMasterType` (`:151`): from the handed master block; one type per sideset after cfaf7e8.
- Pre-existing and out of scope: `duplicates_at` is false for buffer|coil (neither side air or ferro) while `Topology::sideset_type` maps buffer to air, so a buffer|coil sideset is typed `InterfaceAirCoil` by the topology and never seen by this processor's checks.
- `Topology::detect_sideset_types` (`cl_Topology.cpp:229-237`) classifies a sideset by its first facet; its coil refusal is now backed by this processor's, but its domain-type label of the sideset (used downstream of the topology) still follows the first facet. Other file; reported.

## 6. Audit Trail

- `tmp/ai_exchange/review_interface_first_facet_gaps.md` — plan and code juries.
