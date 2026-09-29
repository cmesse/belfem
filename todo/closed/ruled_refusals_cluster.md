# Ruled Refusals Cluster: DR-162, DR-171, DR-172 as One Plan Round

**Date:** 2026-09-24
**Purpose:** Three debt rows whose rulings are already taken land through one plan jury and one
code jury: the serial side-connector guard (DR-162, R1 of `todo/test_gap_followups.md`), the
material-deck contradiction refusals with both input artifacts (DR-171, R9), and the narrow
Dirichlet refusal on a non-duplicate multi-source hanging dof (DR-172, R10). Each row keeps its
own commit; the cluster shares the audit rounds because the diffs are small and disjoint.
**Module:** `src/fem/kernel` (DR-162, DR-172), `src/physics/materials` (DR-171), `tests/fem`,
`tests/physics`, `doc/input_file_reference.md`, `doc/input_schema.yaml`
**AIs involved:** Claude (plan, orchestration, commits), Codex `gpt-6-astra` + Grok `grok-4.7`
(plan and code juries), Opus (implementation)
**Status:** ✅ COMPLETE (2026-09-29). All three rows landed on `devel` in three commits behind one
plan jury (Codex `gpt-6-astra`/high + Grok `grok-4.7`/high, 20 findings, none refuted), one code
jury (same seats: one P1, the linker entry gate, fixed with a `ConnectorOnly` case; two P2), and a
quick terra/medium check on the fix delta. Verified: `make check` green in the release and the
debug configuration of `cmake-build-debug` (17/17 each; the twelve new or renamed cases listed in
the logs); the R11a probe run on `costheta` (3 Dirichlet sidesets) and `gantry` (13) to the first
solve counted zero multi-source hanging nodes, duplicate or not, on every one, then the probe was
stripped. Residue: the pre-existing per-rank half of the linker entry gate is a new `[W]` register
row; DR-175's scope-owner round follows (the fixture case (c) leaks a `DofManager` on its refusal)

> **Scope guards:**
> - Three rows, three commits, one plan jury, one code jury on the own-files diff of all three.
> - No formulation change and no change to what any shipped deck computes. DR-171 changes what a
>   deck may say (two refusals) and adds two programmatic contracts (R5: `set_constant( jc )`
>   defaults `ec`; R6: `Material::n` falls back to the constant instead of asserting), so both input
>   artifacts move in the same session; nothing else touches `input.conf`.
> - The cohomology core is not touched. DR-172's duplicate branch (cut duplicates) is documented,
>   not changed.
> - DR-175 (the factory's leaks on refusal) is NOT in scope. Both new refusals in this plan run
>   before any allocation, so they add no leak of that class.

---

## 1. Current Behaviour and How It Fails

| Row | Failure | Mechanism | Evidence |
|---|---|---|---|
| DR-162 | a partial thin-shell selection with side connectors aborts at setup (throw in the test binaries and in debug, `error_abort` in release; `Map::operator()` is `BELFEM_ERROR` in every build, `cl_Map.hpp:224-235`) | the serial linker guards layer blocks (`cl_FEM_DofMgr_BlockData.cpp:465`) but looks the connector block up unguarded (`:481`) and the recovery facet's master block unguarded (`:500-502`); the parallel twin guards both (`:617`, `:660-661`) | source trace |
| DR-171 K-1 | `builtin : pb40sn60rrr5 ; RRR : 100 ;` runs at RRR 5 silently | the factory reads `RRR` (`cl_MaterialFactory.cpp:118-121`) and forwards it, but `Alloy::Alloy` keeps the in-name token and drops the argument when the name carries one (`cl_Material_Alloy.cpp:56-71`) | source trace |
| DR-171 K-2 | `builtin : copper ; curve : RoxieIron ;` yields an Iron under the key `copper` | the curve branch tolerates the key `builtin` with any value (`cl_MaterialFactory.cpp:59-87`, allow-list by name only); the parser stores a valueless key as the string `"true"` (`cl_Input_Section.cpp:133-140`), so the two forms are distinguishable | source trace |
| DR-171 Q3 | a programmatic `set_constant( jc )` material evaluates `rho_powerlaw` to NaN in release | `set_jc_function` defaults `ec` (`cl_Material.cpp:691-695`), `set_constant` does not (`:257`); `constant_property( ec )` is `BELFEM_QUIET_NAN` until set (`:162`), the accessor asserts only in debug (`cl_Material.hpp:1677`) | source trace |
| DR-171 Q4 | `n( B, angle )` on a constant-n material asserts in debug, dereferences null in release | `Material::n` has no constants fallback (`powerlaws.hpp:994-1013`, dox block `:982-992`) while `Material::jc` does (`:1027-1055`) and the assembly path `n_eval_raw` does (`:93-101`) | source trace |
| DR-172 | a Dirichlet on a non-duplicate hanging dof with more than one source is silently lost, on BOTH imposition paths | `pin_dirichlet_dof` handles only the single-source unit-weight shape (`cl_FEM_SideSet.cpp:47-91`); the per-timestep fallback fixes the alias (`:344-350`) and so does the factory pre-fix loop (`cl_MaxwellFactory.cpp:909-936`); the alias has left the regular container, and `compute_hanging_dofs` overwrites it after every solve (`cl_FEM_DofMgr_SolverData.cpp:2595-2609`) | source trace; `tests/fem/test_DirichletPin.cpp:393-419` documents the loss |

**Bottom line:** three contracts that the code enforces on one path and not on its twin. The fixes
are the missing guard, the missing refusal, and the missing fallback, each mirrored from the path
that already has it.

## 2. Rulings (all Christian, quiz)

- **DR-162, 2026-09-22:** partial selection is SUPPORTED. Guard the serial connector loop as the
  parallel twin does; add a side-connector fixture.
- **DR-171, 2026-09-22:** a contradiction is refused, agreeing redundancy is allowed.
  `builtin : pb40sn60rrr5 ; RRR : 5 ;` passes; `RRR : 100` beside it is refused naming both
  values; a `curve` section may not carry a `builtin` key that names a different material.
  **2026-09-24, O1 → rule C:** the curve check refuses only a value that names a KNOWN non-iron
  builtin or an alloy formula; valueless, the iron aliases, and unknown names (a user's nameplate
  such as `aisi305`) pass.
- **DR-172, 2026-09-24:** the NARROW refusal. `BELFEM_ERROR` on
  `is_hanging() && number_of_sources() > 1 && ! is_duplicate()`, naming sideset, node and sources
  as `Bearing::impose_dirichlet` does (`cl_FEM_Bearing.cpp:151-157`); the duplicate branch stays
  inert and is documented as deliberate.

## 3. Consumer Inventory

| Row | Producer / deck population | Impact of the change |
|---|---|---|
| DR-162 | `set_block_types_in_magnetic_equation` hands EVERY non-Coil block, coatings included, to the equation (`cl_MaxwellFactory.cpp:2176-2212`; `create_field` uses that list, `:786`); the kernel-parameter filter at `:745-765` omits coatings on purpose and is not the dof-manager selection. The thermal IWG leaves the connector loop before it runs (`cl_FEM_DofMgr_BlockData.cpp:476-477`). A partial selection reaches the connector loop only from a programmatic caller (a test) | none on shipped decks; the guard is inert where the factory selects everything (high) |
| DR-171 K-1 | the three shipped alloy sections (`examples/pancake:76-77`, `tapestack3d:109-110`, `tapestack3d_quench:145-146`) write `builtin : Pb38Sn62 ; RRR : 1.5 ;` with no in-name token; no deck writes both | none (high) |
| DR-171 K-2 | two curve sections carry the valueless `builtin ;` (`examples/gantry:60-62`, `examples/costheta:44-46`); `undulator2d:52` carries `curve` alone; no deck carries a valued `builtin` beside `curve` | none; the tolerated forms stay tolerated, and a user nameplate stays legal under rule C (high) |
| DR-171 Q3 | programmatic `set_constant( jc )`: `tests/physics/test_MaterialDerivatives.cpp:157-161`, `:680-687` (both set `ec` first, so the default never engages); the factory sets `ec` after `jc` and `n` (`cl_MaterialFactory.cpp:141-153`), and a later `set_constant( ec )` overwrites the default; the only `have( ec )` reader in the tree is `set_jc_function` | none (high) |
| DR-171 Q4 | `Material::n( B, angle )` callers outside `powerlaws.hpp`: `tests/physics/test_JcCustomGuard.cpp:103` (function-backed) | none; the fallback engages only where the code asserted before (high) |
| DR-172 | `tmp/ai_exchange/consumer_impact_dr172.md`: every multi-source hanging node a shipped deck builds is a cut or interface DUPLICATE (P1-P3, P5 by inheritance); the non-duplicate shape has no producer that reaches a Dirichlet sideset (no refiner exists). Three decks impose a Dirichlet at all (`gantry`, `costheta`, `racetrack_usermat`), all through symmetry domain types, which go through `impose_dirichlet` (`cl_MaxwellFactory.cpp:873-878`) | the narrow form fires on nothing reachable today (medium; the R11 run with the R11a count settles it) |

## 4. Ordered Steps

- [x] **R1 (DR-162) — the serial guard and the selection refusal.** In `link_thin_shell_facets_serial`
  (`cl_FEM_DofMgr_BlockData.cpp:479-505`): `if ( ! mBlockMap.key_exists( tMeshBlock->id() ) ) continue ;`
  before the connector block lookup at `:481`, the parallel twin's rule (`:617`). For the master
  block at `:500-502` the twin's silent skip is NOT mirrored: a selected wall whose layer block is
  unselected would keep `mReferenceElement == nullptr` (`cl_FEM_Element.hpp:55`) and
  `prepare_side_connector_frame` dereferences it unconditionally (`cl_FEM_Calculator.cpp:487-489`).
  Setup must not succeed into that crash (`doc/coding_philosophy.md:639`), so the serial linker
  refuses: `BELFEM_ERROR( mBlockMap.key_exists( tMaster->block_id() ), "side connector block %lu is selected but the layer block %lu of its recovery facet %lu is not - select the shell's layer blocks with its connectors" )`.
  Per O4 (ruled), the parallel twin's guard at `:660-661` is split the same way, with one
  correction found at implementation: the block map is NOT rank-uniform (`create_blocks`,
  `:113-129`, registers a block only when it has owned or aura elements on this rank), so the
  parallel refusal tests the IWG's selection (`iwg()->selected_blocks()`, broadcast by the IWG),
  which every rank holds alike; the `key_exists` and `element_exists` halves both stay a
  `continue` (an aura copy whose master did not ship). In serial the block-map test is equivalent,
  because rank 0 owns every element. The refusal is raised only on ranks that hold the wall; in
  release that aborts the job, in throw mode the other ranks stop at the next collective, the same
  class as the existing recovery-facet error in that loop. `set_facet` at `:497` stays before the master check in both. Comment above the loop: the
  connector block and the master block are guarded independently because a programmatic caller
  may select one without the other; the factory selects both.
- [x] **R2 (DR-162) — the fixture.** `tests/fem/test_ThinShellSelection.cpp`: give
  `TwoShellStacks` an optional side connector on stack 0, the production shape reduced to one
  element: a `HEX8TB` wall element in a coating block (`DomainType::LeftCoating`, pushed to
  `side_connector_blocks()` and to the mesh blocks), and a recovery facet whose `Facet` wraps a
  `QUAD4` with id `wall id + 1` (`cl_ThinShellFactory.cpp:555-557`: wall `++mMaxElementID`, face
  `++mMaxElementID`), master = the lower layer prism (`set_master( ..., true )`), slave = the wall,
  its edge container holding the two longitudinal edges in slots 0 and 1 (`:630-636`), in a sideset
  that is BOTH pushed to `side_connector_sidesets()` AND added to the mesh (`:643`), because
  `Mesh::facet` reads `mFacetMap`, which `finalize()` fills from the mesh sidesets
  (`cl_Mesh.hpp:1443-1451`, `cl_Mesh.cpp:743-746`). Three cases: (a)
  `ConnectorOfUnselectedShellIsSkipped`: select stack 1 only, `create_field()` does not throw;
  (b) `SelectedConnectorLinksFacetAndReference`: select all blocks including the coating; the
  wall's fem element carries the recovery facet and its reference is the layer element; (c)
  `ConnectorWithoutItsLayerBlockIsRefused`: select the coating block and stack 1, not blocks 1-2;
  `expect_refusal` on the R1 message. Cases (b) and (c) need `IWG_Maxwell` to size a
  `LeftCoating` element (`cl_IWG_Maxwell.cpp:431-438`), i.e. the wall's edges populated as the
  factory does. If that cannot be built without reproducing the factory's side-layer bookkeeping,
  the fork is O3.
- [x] **R3 (DR-171 K-1) — RRR contradiction, on the NAME, before construction.** `cl_MaterialFactory.cpp`,
  inside `if ( tIsBuiltin )` BEFORE `create_material` (`:121`): decompose `tMatType` with `to_pair`
  (`src/core/stringtools.hpp:313-331`) exactly as `create_material` will (`:520-525`), find a
  component whose lowercased label is `rrr` (the test `Alloy::Alloy` applies, `cl_Material_Alloy.cpp:61`);
  if one exists and the section carries `RRR` and `real( token ) != tRRR`, refuse naming the
  material label, the deck value and the token value. Placed here it never sees YBCO (which
  hardcodes `set_RRR( 50.0 )`, `cl_Material_YBCO.cpp:118`, and would have tripped an
  `is_constant( RRR )` test) and it allocates nothing, so no DR-175-class leak. Exact comparison:
  the token is an INTEGER (`std::stoi` on `\d+`, `stringtools.hpp:319-328`; `rrr1.5` parses as 1),
  the deck value is `strtod` (`stringtools.cpp:317-321`); an agreeing integer pair is bit-identical,
  and `pb40sn60rrr1.5` beside `RRR : 1.5` is refused as `1 != 1.5`, which is the truth about the
  parser and belongs in the message and in R8.
- [x] **R4 (DR-171 K-2) — curve beside a KNOWN non-iron builtin (rule C).** In the curve branch,
  BEFORE `new material::Iron` (`:64`) and before any file is opened: read the `builtin` value,
  LOWERCASE it (`get_string` returns the stored text unchanged, `cl_Input_Section.cpp:256-267`;
  the valueless marker is stored lowercase `"true"`, `:138`, so `builtin : TRUE` also reads as the
  marker), and refuse only when it names a builtin that `create_material` would build as something
  other than `Iron`. Tolerated: `"true"`; the iron aliases `iron`/`ferro`/`fe` (`:498`); any name
  the factory does not know (`aisi305`). Refused: the pure-metal aliases (`:493-512`),
  `hastelloy`/`hastelloyc276`, `ybco`, `magnesia`/`mgo`/`buffer`, and an alloy FORMULA, recognised
  as `to_pair` yielding components whose labels are ALL pure-metal aliases or `rrr`
  (`to_pair( "aisi305" )` yields `aisi`, no metal, so unknown, so tolerated). The name predicate
  is one private static `MaterialFactory::builtin_kind( const string & )` (returns Iron / OtherKnown /
  Unknown) that `create_material` and the curve branch both use, so the alias tables cannot
  drift. The message names the label, the builtin value and the curve. The allow-list call at
  `:77-80` stays.
- [x] **R5 (DR-171 Q3) — `ec` default in `set_constant( jc )`.** `Material::set_constant`
  (`cl_Material.cpp:257`): after storing, if `aProperty == MaterialProperty::jc` and
  `! have( ec )`, `set_constant( ec, gDefaultEcrit )`, the same guard as `set_jc_function`
  (`:691-695`). No recursion: the `ec` branch has no jc default, and the `switch` at `:272-300`
  only cross-fills M/R/beta/q. The literal `1e-4` becomes one named constant used by all three
  sites (`:694`, the new one, `cl_MaterialFactory.cpp:152`).
- [x] **R6 (DR-171 Q4) — constants fallback in `Material::n`.** `powerlaws.hpp:994-1013`: both
  overloads take the `jc` shape (`:1027-1055`): `mNFunction != nullptr ? eval : constant_property( n )`.
  The dox block at `:982-992` loses the "use `constant_property` instead" sentence and the `\pre`;
  it gains the fallback sentence `jc` carries.
- [x] **R7 (DR-171) — tests.** `tests/physics/test_MaterialFactoryDeck.cpp`: in `AlloyRRR`, the
  agreeing pair `builtin : pb40sn60rrr5 ; RRR : 5 ;` constructs with RRR 5, and `RRR : 100` beside
  the same name is refused with a message carrying both `5` and `100`. New suite `CurveShape`:
  the name predicate is tested DIRECTLY on strings, no deck and no database (`"true"` → tolerated,
  `"iron"`/`"FERRO"` → Iron, `"aisi305"` → Unknown, `"copper"`/`"YBCO"`/`"Pb38Sn62"` → OtherKnown);
  plus one deck-level case, `builtin : copper ; curve : RoxieIron ;`, refused with a message naming
  `copper` and `RoxieIron`, which fires before any file is opened and so is independent of whether
  a `bhdata.hdf5` is reachable through `material::data_file`. In `ConstantLaws`, a programmatic
  `create_material( "ybco" )` with `set_constant( jc )` and `set_constant( n )` and no `ec` has
  `is_constant( ec )` true with the named default and a finite `rho_powerlaw`, and `n( 0.7, 0.4 )`
  and `n( 0.7, 0.4, 77.0 )` return the constant.
- [x] **R8 (DR-171) — input contract.** `doc/input_file_reference.md:484-486` (the exemption
  paragraph), `:491` and the `RRR` row at `:545`; `doc/input_schema.yaml:896-899`, `:915`, `:919`,
  the `rrr_token` note at `:1001-1002`, the `RRR` note at `:1028-1033`. New facts: (1) `RRR`
  beside an alloy name that carries an `rrr<N>` token must agree with it, a differing value is
  refused naming both; the token is an integer, so a fractional RRR belongs in the key, not the
  name; (2) a `curve` section may carry `builtin` valueless, naming iron (`iron`/`ferro`/`fe`), or
  naming something the factory does not know; a value that names a known non-iron builtin or an
  alloy formula is refused; (3) `builtin : true ;` (any case) is indistinguishable from the
  valueless key at the parser and is therefore tolerated. Then the Codex language sweep over the
  touched `.md` paragraphs at **terra/medium** (the input reference is the dense-document row,
  `doc/ai_collaboration_protocol.md:461`); no sweep over the yaml.
- [x] **R9 (DR-172) — the narrow refusal, in the shared helper.** `fem::pin_dirichlet_dof`
  (`cl_FEM_SideSet.cpp:47-49`): after the existing shape test returns false for a non-hanging dof,
  refuse a hanging dof with more than one source:
  `BELFEM_ERROR( ! ( aDof->is_hanging() && aDof->number_of_sources() > 1 ), "Dirichlet group %lu ( node %lu ): the dof is condensed onto %u sources; a Dirichlet value on such an alias is not enforced ( the reconstruction overwrites it after every solve ) - only a single-source unit-weight condensation or a cut duplicate is supported on a Dirichlet group", ... )`.
  Both callers exclude duplicates before calling (`cl_FEM_SideSet.cpp:344`, `cl_MaxwellFactory.cpp:913`),
  which becomes part of the helper's contract in its header comment; so one edit covers the
  per-timestep path and the factory pre-fix loop, and the factory comment at `:926-929` stays
  true. The message states what the predicate tests (duplicate, hanging, source count) and
  nothing about where the sources lie. `BELFEM_ERROR`, not `BELFEM_ASSERT`, because a silent
  release build is the defect itself (`doc/coding_philosophy.md:621-624`); the cost is a few
  branches per sideset node per imposition. The comment block at `cl_FEM_SideSet.cpp:337-343` is
  rewritten to state the branches: plain dofs take the direct fix; a single-source unit-weight
  condensation is rerouted to its source; a DUPLICATE keeps the inert alias fix on purpose (a cut
  duplicate reconstructs as its original plus the cut currents at unit weight, `cl_CutSet.cpp:50-53`,
  and fixing the alias would force those currents to zero); the remaining hanging shape is refused
  in the helper.
- [x] **R10 (DR-172) — tests.** `tests/fem/test_DirichletPin.cpp`: `TwoSourcesFallBackToTheAliasFix`
  becomes `TwoSourcesAreRefused` (`expect_refusal` on the new message; the existing
  `number_of_sources() == 2` asserts kept; neither source dof fixed afterwards); new
  `DuplicateWithTwoSourcesKeepsTheAliasFix`: the same `GluedStrip( 1.0, true )` with node 6 marked
  `set_original( node 3 )` BEFORE `new Kernel` (`cl_Node.hpp:297-304`; the Kernel redirects element
  adjacency through `original()`, `cl_FEM_Kernel.cpp:257-276`, harmless on one rank, and
  `set_original` registers no reverse link, so node 3 walks no duplicates); imposition does not
  throw, the alias carries the fixed flag, both sources stay free. The fixture comment says why the
  duplicate case is inert by design and that the flag check before any solve is the same limit the
  replaced test had.
- [x] **R11 (all) — gates.** Per row, in both trees: `test_fem --gtest_filter='ThinShellSelection.*'`
  (R1-R2); `test_physics` whole binary, and `--gtest_filter='AlloyRRR.*:CurveShape.*:ConstantLaws.*'`
  must list every new case (a filter that matches nothing still exits 0, so the run log is read,
  not just the exit code; `test_MaterialDerivatives` is a source of `test_physics`,
  `tests/physics/CMakeLists.txt:14-21`) (R3-R7); `test_fem --gtest_filter='DirichletPin.*'`
  (R9-R10); then `make check-fast`. Christian runs all of these.
- [x] **R11a (DR-172) — the external gate with its count.** A diagnostic probe (IN PLACE
  2026-09-28, `SideSet::impose_dirichlet`, factory-time only; exempt from the round, stripped after
  the run) in `SideSet::impose_dirichlet` counts, per sideset and once, the
  hanging nodes with more than one source, split into duplicates and non-duplicates, and prints
  the two numbers. `gantry` and `costheta` run to the first solve. Oracle: the non-duplicate count
  is zero on every Dirichlet sideset (the producer claim of `consumer_impact_dr172.md` §1 and §3
  holds) and the refusal does not fire; a nonzero non-duplicate count or a firing run REFUTES that
  claim and stops the commit. The duplicate count is recorded in the devlog as the population the
  inert branch serves; a clean run does not by itself prove the anchors are duplicates, the count
  does.
- [x] **R12 — close out.** Three commits (`fem/kernel: ... DR-162`, `materials + input contract:
  ... DR-171`, `fem/kernel + maxwell: ... DR-172`), each after its gate; devlog; strike the three
  rows in `todo/debt_register.md`; tick R1, R9, R10 in `todo/test_gap_followups.md`; the register
  preamble recount: the `[F]` list drops DR-171 and the `[P]` count drops by two (`:19-21`).

## 5. Open Design Questions

- **O1** K-2 width. **RESOLVED 2026-09-24 (Christian, quiz) → rule C** (see §2).
- **O2** K-1 comparison: exact inequality. **RESOLVED by the jury's finding on the parser:** the
  token is an integer, the key a double; exact is the only honest compare, and the integer-token
  fact goes into the message and the reference (R3, R8).
- **O3** R2 cases (b) and (c) need a wall the IWG can size. If the fixture cannot populate the
  wall's edges without reproducing the factory's side-layer bookkeeping: land case (a) alone and
  keep DR-162 OPEN with the fixture clause outstanding (recommended: the master-block refusal
  would then be untested, and the row's ruling asked for a connector fixture), or strike the row on
  (a) plus the serial/parallel symmetry argument. **RESOLVED 2026-09-24 (Christian, quiz): land
  the guard, keep DR-162 OPEN with the fixture clause outstanding.**
- **O4** The master-block hole in the PARALLEL twin (`:660-661`). Refuse there too, so both twins
  hold the same contract (recommended: the block map is rank-uniform, so the refusal is
  collective-safe, and the parallel `element_exists` skip for aura copies stays), or leave the
  parallel twin as it is and note the asymmetry in the linker comment. **RESOLVED 2026-09-24
  (Christian, quiz): refuse in both twins.** Implemented with the IWG selection as the parallel
  predicate, see R1: the plan's "block map identical on every rank" was false (found by the
  implementation agent, confirmed against `create_blocks`).

## 6. What Could Break

- R1: nothing on the factory path (both blocks always present). A programmatic caller that
  selects a coating block without its layer block now gets a setup refusal instead of a later null
  dereference in the calculator. Under O4, the same on every rank that holds the wall.
- R3: a deck that wrote `RRR` beside an `rrr<N>` name and relied on the name winning now stops.
  No shipped deck does. A deck that wrote a fractional RRR into the NAME (`rrr1.5`) beside a
  matching key is refused as `1 != 1.5`; that deck was already running at RRR 1.
- R4: a deck with `builtin : <known non-iron builtin or alloy formula> ;` beside `curve` now
  stops. No shipped deck does. A typo such as `builtin : irom ;` still passes, harmlessly: the
  curve decides the physics either way.
- R5: a caller that set `jc` as a constant, then read `have( ec )` to decide something, sees
  `true` now. The only reader is `set_jc_function` itself.
- R9: a deck that puts a non-duplicate multi-source hanging node on a Dirichlet sideset, or on an
  explicit `dirichlet` boundary condition, now stops at factory time instead of silently losing
  the condition. The consumer sweep found no producer; R11a is the check.

## 7. Definition of Done

- [x] the three code changes land as three commits, each behind its gate in both trees
- [x] the code jury's reconciliation table is appended to `tmp/ai_exchange/review_ruled_refusals_cluster.md`
- [x] both input artifacts updated, `.md` paragraphs swept at terra/medium
- [x] R11a's two counts per Dirichlet sideset recorded in the devlog, probe stripped
- [x] register rows DR-162, DR-171, DR-172 struck (DR-162 conditional on O3); `[F]` line and `[P]` count updated; `test_gap_followups.md` R1, R9, R10 ticked
- [x] devlog written and indexed

## 8. Audit Trail

Exchange: `tmp/ai_exchange/review_ruled_refusals_cluster.md` (plan jury 2026-09-24: Codex
`gpt-6-astra`/high, Grok `grok-4.7`/high; 20 findings, all CONFIRMED, none refuted; the O1
ruling appended mid-round; code jury to follow in the same thread). Consumer sweep for DR-172:
`tmp/ai_exchange/consumer_impact_dr172.md` (2026-09-22). Rulings: `todo/debt_register.md` rows
DR-162, DR-171, DR-172; `todo/test_gap_followups.md` R1, R9, R10.
