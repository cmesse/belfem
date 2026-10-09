# Interface sidesets, round 2: message class and Topology type from the facet set

**Date:** 2026-10-09
**Purpose:** Christian's rulings of 2026-10-09 on the blind audit of 337a318: keep the "interface facets only" refusal; fix the two small defects of `refuse_mixed_sidesets` with the docs (commit 1); make `Topology::detect_sideset_types` classify a sideset from all of its facets and refuse disagreement (commit 2). Two commits on devel so each audit is attributable.
**Module:** `src/homology` (`cl_InterfaceProcessor`, commit 1; `cl_Topology`, commit 2 — outside the closed core), `tests/homology`, `src/homology/doc`, `doc/`, `CHANGELOG.md`
**AIs involved:** Claude (plan, implementation, gates), Codex + Grok (plan jury; one code jury per commit); the Server session (blind xhigh audit per src/homology diff, pipeline on the final tip); the Yufan session (iron/coil decks plus tapestack3d on the final SHA)
**Status:** ❌ CLOSED BY RULING 2026-10-09 (Christian): "we do not include this in 0.9.3 because it tries to fix something that is not broken, and I would even remove these fixes from devel and archive them." Commit 1 (411e69a) and e30bd98 reverted on devel; the Topology commit never pushed; all of it on `archive/first-facet-sidesets` (93d8318). Reason: a native Gmsh mesh builds one sideset per geometric entity, so no shipped deck has a sideset whose classification depends on facet order (21/21 decks, identical sideset types and dof banners on the strict binary), and the Topology change was not inert (cut-pruning exemptions, dof tables, the Poisson cut-side setup read the auto-typed antisymmetry types). See `devlog/dl20261009_first_facet_campaign_reverted.md`.

> **Scope guards (from the brief):**
> - The "interface facets only" refusal stays.
> - Commit 1 touches `cl_InterfaceProcessor.{hpp,cpp}`, its tests and the docs; commit 2 touches `cl_Topology.{hpp,cpp}` and tests. Nothing else is widened; anything else found is reported.
> - The isolated buffer|coil sideset: report the right behaviour, change nothing without asking.
> - devel only (tip `6cb9685`); never main. No register rows. The shared checkout's local devel carries another session's unpushed `f3841bc` and has diverged from origin; this work is done in a private worktree at `6cb9685` and pushed from there.

---

## 1. Current Behaviour and How It Fails

| Defect | Mechanism | Evidence |
|---|---|---|
| (a) message class depends on facet order | the master-type and coil checks run interleaved in one facet walk (`cl_InterfaceProcessor.cpp:515-548`): [C\|A, C\|K, F\|A] reports the coil mix, [C\|A, F\|A, C\|K] the master-type mix | source trace, Codex, confirmed by the Server |
| (b) wrong reason text | every non-duplicated facet with a slave is labelled "neither side is air or ferro" (`:575`), but an air\|air facet fails `duplicates_at` because both sides are air (`:474`) | both auditors |
| (c) stale wording | header `:148-155` and the numbering test comment say set order is "the sorted hex strings of the block-index bitsets"; the key is `hex(masters)\|hex(slaves)` | both auditors |
| (d) untested claims | a multi-block CONDUCTOR sideset in both facet orders giving the same mesh is untested; `OverlappingFerroSideSetsStayApart` also passes on `0d4d1b5` (first-facet keys already separated it); "every facet touches a coil → accepted and decoupled" is untested | both auditors |
| (e) doc errors | "admitted when any facet is an interface" (false since the interface-only rule); "before any copy exists" (cut copies exist; the check runs after `CutProcessor`, on the sidesets counted before it); "different master blocks never share a set" (it is the `(master set, slave set)` key); the `.bfm` bypass missing from the guide's section; the numbering promise too wide — corrected after the jury: `mCutSets` is an `OrderedMap` and `CutSet` walks its originals sorted by id (`cl_CutProcessor.cpp:700`), so cut-copy IDS are deterministic, but `collect_duplicates` appends them in the hash order of each set's duplicate map (`:1214-1230`), so their POSITIONS, and the edge keys built from positions, still depend on the platform; the interface copies are deterministic for a given incoming node order | both auditors, Grok's correction |
| (f) Topology first facet | `detect_sideset_types` (`cl_Topology.cpp:219-246`) types a `Default` sideset from `facets()(0)`: a boundary facet first → `ConductorAntiSymmetry` (`:264-266`); a later conductor\|coil facet is never read, so "a conductor and a coil must not touch" (`:329`, `:359-360`) never fires, and the processor skips the sideset (no duplicated facet) | both auditors |

## 2. Mechanism

### 2.1 Commit 1 — `refuse_mixed_sidesets` in fixed precedence

Four complete passes per interface sideset, each over the whole facet list before the next starts, so the class reported is a function of the facet set: (1) master domain type over slave-bearing facets; (2) coil contact over slave-bearing facets; (3) interface facets only — the non-duplicated facet's reason is one of three: "a boundary facet without a neighbor" (no slave), "both sides are air", "neither side is air or ferro" (the second is unreachable after pass 1, since only an air|air facet has an Air master; the label function stays honest); (4) one master block for non-conductor masters. Within a pass the two facets named are the first of each kind in facet order — the ids named may vary with order, the class and the remedy do not. A test permutes [C|A, C|K, F|A] over all six orders and pins the master-type class on every order.

**Isolation from Topology (jury round 1):** `Junction::run` runs `Topology::run` before the processor, as production does, so after commit 2 a mixed sideset dies in Topology first. Topology skips explicitly typed sidesets, so the fixture gets `aTypedSideSets`, which pre-types every sideset `Inactive`; the processor never reads sideset types and `Topology::mTypeMap` holds every type. Every processor refusal test (existing and new) runs pre-typed; full-pipeline tests keep `Default` and, in commit 2, pin Topology's message.

Header and test comments: set order is "the sorted keys `hex(master blocks)|hex(slave blocks)`, a fixed function of the mesh's block order".

Tests (d): `OneSideSet` in both facet orders gives the same copy ids, mesh positions and element connectivity (the one positive claim); a sibling of the overlapping-ferro test red on the old union key — two ferro sidesets with the SAME block pair and swapped masters on one face (F1|F2 and F2|F1 only, no air facets): the union key merged them, the master-aware key keeps them apart; a coil|air sideset whose every facet touches the coil is accepted and its set decoupled (copies without sources, `InterfaceTreatment::Decouple`).

Docs (e): the admission sentence ("admitted when any facet passes `duplicates_at`, after the refusal has run, which has made all of them interface facets"); "before any interface copy exists; after `CutProcessor`, on the sidesets present before it"; "two sidesets with different `(master set, slave set)` keys never share a set"; the `.bfm` bypass in the guide's section; the numbering promise: interface copies get the same ids and positions on every platform for the same incoming node order; cut copies get deterministic ids but are appended in hash order, so their positions and the edge keys built from positions still depend on the platform; CHANGELOG: the interface-only rule stated as a new restriction against 0.9.2 (it is already a Refusals bullet; the sentence is sharpened), and the `### Changed` note states the 0.9.3 rule that set order follows the `masters|slaves` key — not attributed to commit 1, which does not change the key.

### 2.2 Commit 2 — `Topology::detect_sideset_types` from all facets

For a `Default` sideset with facets: FIRST compute every facet's type — `sideset_type( master, slave )` for a facet with a slave, `sideset_type( master )` for a boundary facet — over ALL facets (the conductor|coil `BELFEM_ERROR` inside `sideset_type` fires on any such facet wherever it sits, so the contact class wins in every order); THEN, if the types differ, an always-active `BELFEM_ERROR` names the sideset, the first facet and the first facet of another type with both type names, and the remedy ("split the sideset so that each piece has one type"); the sideset takes the common type. `sideset_type`'s Ferro arm returns `InterfaceFerroCoil` for a Coil slave, as the Coil arm does for a Ferro slave: Topology runs before `fix_facet_masters`, so a ferro|coil sideset can carry masters on both sides, and the asymmetric table (`FerroAntiSymmetry` for a ferro master) would make the new check refuse a single-pair sideset — a deliberate one-case fix, named here and in the commit. The six permutations of [C|bnd, C|A, C|K] are tested: the contact error in every order. Explicitly typed sidesets (deck-declared symmetry, antisymmetry, thin shells, cuts) keep `domain_type() != Default` and are skipped as today; nothing is checked for them (reported, not changed). Empty sidesets stay `Inactive`.

Neighbours, reported: the isolated buffer|coil sideset — Topology maps Buffer to Air and types it `InterfaceAirCoil`, the processor's `duplicates_at` is false for it (neither side air or ferro), so it is typed as a coil interface but never duplicated; the right behaviour is Christian's call (either map Buffer to Air in `duplicates_at` too, or type it `Inactive`); unchanged here. The 2-D ferro|ferro curve: `fix_facet_masters` does not propagate same-type orientation in 2-D (`cl_MaxwellFactory.cpp:1459-1464`), so masters may sit on both blocks; the Topology type of such a facet is `Inactive` whichever side is master (Ferro|Ferro → `Inactive`), so the new check cannot fire on it; the processor's one-master-block refusal did not fire on costheta or undulator2d in the 21-deck gate on 337a318 — §5 says why from the meshes.

## 3. Ordered Steps

- [x] **R1** (commit 1) `refuse_mixed_sidesets` as four complete passes; the three reasons; header and test comments.
- [x] **R2** (commit 1) tests, all pre-typed (processor-only): `MessageClassIsAFunctionOfTheFacetSet` (six permutations), `OneSideSetInBothFacetOrdersGivesOneMesh`, `SwappedMasterSideSetsStayApart` (red on the union key), `CoilSideSetIsAcceptedAndDecoupled` (ordinary originals; the inherited sources of a hanging original are the DR-198 residue, reported); the existing processor refusal tests switched to pre-typed sidesets so they survive commit 2.
- [x] **R3** (commit 1) docs per §2.1; luna sweep on the user-facing passages.
- [x] **R4** (commit 1) gates: debug `test_homology` + `make check`; the InterfaceProcessor suite under `USE_DEBUG=OFF` (release `test_homology`); code jury on the diff; push; SHA to the Server.
- [ ] **R5** (commit 2) `Topology::detect_sideset_types` from all facets with the disagreement refusal, the symmetric Ferro|Coil case; full-pipeline tests (`Default` sidesets) on the junction: the six permutations of [C|bnd, C|A, C|K] → the contact refusal; conductor|air with a boundary facet, with conductor|conductor, with air|air → the disagreement refusal in both orders; conductor|air with conductor|buffer → one Topology type (the processor then refuses it, buffer not duplicated at — a Topology-level acceptance only); a ferro|coil sideset with masters on both sides → one type.
- [ ] **R6** (commit 2) gates: debug `test_homology` + `make check`; the Topology suite under `USE_DEBUG=OFF`; release, all 21 decks through setup with no `.bfm` (the gate script deletes them; corc_twolayer's copy regenerated with 6cb9685's generator) and, per deck, the sideset types against 6cb9685 through a scratchpad-only probe printing `sideset id → domain type` after `detect_sideset_types` on both binaries, diffed; code jury; push; SHA to the Server.
- [ ] **R7** close out: devlog, index lines, this file to `todo/closed/`.

## 4. What could break

- Commit 2 refuses a shipped deck whose sideset mixes types (boundary + interface facets, or conductor|conductor + conductor|air): a finding for Christian. The 21-deck gate on 337a318 passed the processor's interface-only rule, which already refuses those mixtures when a duplicated facet is present; Topology's rule also covers sidesets with no duplicated facet (e.g. boundary + conductor|conductor), so a deck can still trip it. [medium]
- The sideset-type comparison per deck needs the types to be visible; a scratchpad probe if they are not printed. [high]

## 5. The 2-D iron decks

From the meshes the shipped  files produce ( over the 4.1 files, ids = entity tags, every line element mapped to the blocks of the two faces that share its edge): costheta has 338 sidesets and undulator2d 953, and in both decks every sideset borders exactly ONE block pair. costheta's iron blocks 69 and 70 meet only the air block 72 (sidesets 273–286 and their siblings); no sideset borders (69, 70). undulator2d's ferro blocks 2 and 3 meet the air blocks (62, 63, 122, 123, …) and the domain boundary (sidesets 10 and 11 are one-sided); no sideset borders (2, 3). No sideset in either deck holds both a one-sided and a two-sided element (a one-sided sideset's pair is `(block,)` for every element, a two-sided one's is `(block, block)` for every element). So neither deck has a ferro|ferro curve under a sideset, the 2-D orientation caveat has no sideset to act on, every sideset has one master block, and both the processor's one-master-block refusal and the new Topology disagreement check see one facet type per sideset. That is why the 21-deck gate on 337a318 passed them and why commit 2 cannot change their sideset types.

## 6. Audit Trail

- `tmp/ai_exchange/review_interface_order_independence_2.md` — plan jury and the two code juries.
