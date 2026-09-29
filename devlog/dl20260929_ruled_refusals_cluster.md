# Ruled refusals cluster: DR-162, DR-171, DR-172 landed as three commits behind one plan jury and one code jury

**Date:** 2026-09-29 (plan and juries 2026-09-24, fix step 2026-09-28, gates and close 2026-09-29)
**Purpose:** clean up the debt register by the cheapest closable cluster: three rows whose rulings
were already taken. The serial side-connector guard (DR-162), the material-deck contradiction
refusals with both input artifacts (DR-171), and the narrow Dirichlet refusal on a non-duplicate
multi-source hanging dof (DR-172).
**Module:** `src/fem/kernel` (`cl_FEM_DofMgr_BlockData`, `cl_FEM_SideSet`), `src/physics/materials`
(`cl_MaterialFactory`, `cl_Material`, `powerlaws.hpp`), `tests/fem`, `tests/physics`, the input contract
**AIs involved:** Claude Fable 5.1 (plan, orchestration, verification, commits); Codex `gpt-6-astra`/high
and Grok `grok-4.7`/high (plan jury and code jury, blind); Codex `gpt-5.6-terra`/medium (quick check on
the fix delta, and the language sweep of the reference paragraphs); Claude Opus (five implementation
agents, every diff read by Fable before dispatch)
**Status:** COMPLETE — three commits on `devel` (`bbbe34c` DR-162, `7a16fcd` DR-171, the DR-172 commit
carries this entry); rows struck; DR-190 filed
**Verification:** `make check` green in the release and the debug configuration of `cmake-build-debug`
(17/17 each, 2026-09-28 and 2026-09-29; the twelve new or renamed cases listed in `LastTest.log`); the
R11a probe run in the sandbox on `costheta` and `gantry` to the first solve (debug binary, meshes from the
decks' own `.geo`): zero multi-source hanging nodes, duplicate or not, on all 16 Dirichlet sidesets
**Exchange:** `tmp/ai_exchange/review_ruled_refusals_cluster.md` (plan round, code round, fix step, quick
check, all in one thread), `review_ruled_refusals_cluster{,_v2,_fix}.diff`, `sweep_input_reference_dr171`

## Summary

Christian asked which register cluster could close next. Of the nine rows opened since 0.9.0, three
had rulings and disjoint small diffs: DR-162, DR-171, DR-172. They shared one plan
(`todo/closed/ruled_refusals_cluster.md`), one plan jury, one code jury and a quick check on the fix
delta, and landed as three commits. The error-path-hygiene trio DR-165/167/175 is the natural next
cluster, ordered after this one so DR-175's scope-owner pass covers the new refusal paths once.

Two rulings were taken in flight. **O1, the width of the curve exemption:** Christian raised the
nameplate case, a user's own steel `aisi305` loaded with the shipped `RoxieIron` curve, which the draft
would have refused; the ruling is rule C, refuse only a `builtin` value that names a KNOWN non-iron
builtin or a formula of known pure metals. **O4, the parallel twin:** refuse in both linkers.

## What landed

**DR-162 (`bbbe34c`).** The serial connector loop skips an unselected connector block by `key_exists`,
as it already did for layer blocks, and REFUSES a selected connector whose layer block is unselected:
a wall without a reference element is dereferenced unconditionally in
`MaxwellData::prepare_side_connector_frame` (`cl_FEM_Calculator.cpp:487-489`), so mirroring the parallel
twin's silent skip would have let setup succeed into a crash (Codex 1, Grok 14, own F2 of the plan
round). The parallel twin refuses on the same shape through the IWG's broadcast block selection. The
plan's premise that `mBlockMap` is rank-uniform was FALSE: `create_blocks` registers a block on a rank
only when that rank has owned or aura elements of it (`cl_FEM_DofMgr_BlockData.cpp:113-129`); the
implementation agent caught it. The code jury's P1: `create_blocks` entered the linker only when a
shell had a selected LAYER block (collector `:141-170`), so a connector-only selection never reached
either refusal; the gate now also opens when any selected block has a coating type
(`IWG::block_type`, map broadcast in `IWG_Maxwell::initialize`). `TwoShellStacks` carries a HEX8TB wall
with its recovery facet at wall id + 1 on a mesh sideset (Grok: `Mesh::facet` reads the facet map that
`finalize()` fills from sidesets only) and four cases; the O3 fallback was not needed.

**DR-171 (`7a16fcd`).** K-1: the RRR check moved onto the NAME before `create_material`, because the
draft's post-construction `is_constant( RRR )` test sat in the HTS arm and would have fired on YBCO's
hardcoded `set_RRR( 50.0 )` (Grok 3); placed there it allocates nothing. The in-name token is
`std::stoi` on `\d+`, the key is `strtod`, so the compare is exact and `rrr1.5` means 1 (Codex P3,
Grok 4; my "both through stod" rationale was wrong). K-2: rule C through one alias table shared by
`create_material` and `builtin_kind`, and a static `curve_tolerates_builtin` predicate tested directly
on strings after Grok noted the tolerated set was untested; the classifier tokenizes letter runs
without converting digits after Codex showed `to_pair`'s `stoi` throws on a long nameplate. Q3:
`set_constant( jc )` defaults `ec` through `gDefaultEcrit`, used at all three sites. Q4: both
`Material::n` overloads fall back to the constant. Both input artifacts updated; the three reference
passages swept at terra/medium (dense-document row, protocol §9.1), no suspected errors.

**DR-172 (this commit).** The refusal lives in `pin_dirichlet_dof`, not in the sideset: Grok found
that the Maxwell factory's pre-fix loop (`cl_MaxwellFactory.cpp:909-936`) has the same alias-fix
fallback, so a sideset-only refusal would have left an explicit `dirichlet` boundary condition
alias-fixing the shape at factory time. Both callers exclude duplicates before calling, now a stated
contract. I simplified the agent's helper: `is_hanging()` is `mNumberOfSources > 0`, so its
zero-source guard described an impossible state; the helper reads as the Bearing template does. The
cut-algebra comment was corrected after Grok 3: a fixed alias constrains the SUM of the cut currents,
to zero only when the original is held at the same value. `TwoSourcesFallBackToTheAliasFix` became
`TwoSourcesAreRefused`; `DuplicateWithTwoSourcesKeepsTheAliasFix` marks node 6 a one-way duplicate
before the kernel is built.

## Juries

Plan round, 20 findings, none refuted. Beyond the P1 above: the agreeing-iron test could not rely on
a missing `bhdata.hdf5` because `material::data_file` searches the data root (Codex 3, Grok 12), so the
predicate is tested directly; the DR-172 message must not claim the sources lie outside the sideset,
which the predicate does not test (Codex 4, Grok 8); the external gate must count per sideset, not
just reach a solve (Codex 5, Grok 10); the sweep depth for the input reference is terra/medium, not
luna (Codex P3, `doc/ai_collaboration_protocol.md:461`).

Code round: Codex two P2, Grok three defects and five risks. REFUTED: Grok 1, "`is_pure_metal` does
not compile as C++17" (scoped enums support the built-in relational operators; a probe under
`-std=c++17 -Wall -Werror -pedantic-errors` compiles, and the file had passed the build tree's flags).
Grok 2 objected that the narrow refusal also catches a refinement node whose sources the same sideset
fixes; that is the ruled behaviour, and neither producer Grok cited adds one (the Maxwell tie rewrite at
`:2298-2335` filters existing sources; the `.bfm` reload reproduces what was serialized). Quick check on
the fix delta: the entry gate's PRE-EXISTING half is per rank (`cl_ProtoMesh.cpp:1400-1425` builds a
`ThinShell` record only when the rank holds its sideset and blocks) while the linker is collective;
the half added here is uniform. Filed as **DR-190 `[W]`**, Christian's ruling.

## Gates

`make check` in `cmake-build-debug` as release (2026-09-28 17:28) and, after Christian reconfigured
it, as debug (2026-09-29 14:01): 17/17 both, every new case OK. The tree named `build` was not run;
its last log is from 2026-09-23. R11a in the sandbox with a factory-time probe in
`SideSet::impose_dirichlet`, since stripped:

| deck | Dirichlet sidesets | nodes | duplicate multi-source | non-duplicate multi-source | first timestep |
|---|---|---|---|---|---|
| `costheta` | 330, 331, 332 | 16, 58, 33 | 0, 0, 0 | 0 | succeeded, 27 steps before the kill |
| `gantry` | 3 to 15 (13) | 3 to 53 | 0 on all | 0 | succeeded |

So the consumer sweep's claim holds on both decks, and the cut anchors of `gantry` do not land on its
symmetry line at all, which is what `unflag_symmetry_sidesets` is for.

## Open

- **DR-190** (new, `[W]`): the per-rank half of the linker entry gate. Fix sketch in the row.
- Fixture case (c) and the connector-only case leak a `DofManager` on their refusal (thrown inside
  `Kernel::create_field` before the kernel owns it): DR-175's class, listed there by the leak checker.
- Whether the parallel refusal, raised only on ranks holding the wall, leaves the others at a barrier
  in throw mode: accepted by the debug-must-throw rule; same class as the recovery-facet error beside it.
- `test_gap_followups.md`: R1, R9, R10 ticked; seven steps remain (R2 to R7, R12).

## Files

`src/fem/kernel/cl_FEM_DofMgr_BlockData.cpp`, `cl_FEM_SideSet.cpp/.hpp`;
`src/physics/materials/cl_MaterialFactory.cpp/.hpp`, `cl_Material.cpp/.hpp`, `powerlaws.hpp`;
`tests/fem/test_ThinShellSelection.cpp`, `test_DirichletPin.cpp`, `tests/physics/test_MaterialFactoryDeck.cpp`;
`doc/input_file_reference.md`, `doc/input_schema.yaml`; `todo/debt_register.md` (live 12, archive 173),
`todo/debt_register_closed.md`, `todo/test_gap_followups.md`, `todo/closed/ruled_refusals_cluster.md`,
`todo/README.md`.
