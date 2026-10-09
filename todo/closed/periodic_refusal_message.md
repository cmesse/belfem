# The cross-set periodic refusal names its cause

**Date:** 2026-10-09
**Purpose:** `InterfaceProcessor::pair_cross_set_periodic_duplicates` refuses a periodic pair whose two copies lie on interface sets of different master domain type or different treatment, but reports it as "partner was not duplicated in any interface set". Make the refusal say the real cause, pin all three messages with tests, and document the rule. The filter itself does not change.
**Module:** `src/homology` (tests in `tests/homology`; docs in `src/homology/doc`, `doc/`, `CHANGELOG.md`)
**AIs involved:** Claude (plan + implementation), Codex + Grok (plan jury, code jury); the Server session runs a blind delta audit after the push
**Status:** ✅ COMPLETE 2026-10-09. Landed on devel: `pair_cross_set_periodic_duplicates` remembers the excluded sibling sets that hold the partner's copy, and `refuse_unpaired_periodic_copy` reports the domain-type case, the treatment (coil) case and the genuine "not duplicated" case with both node ids and both labels; what is refused is unchanged. `test_InterfaceProcessor.cpp`: `Layout::Coil`, a `refusal()` helper, four new tests; the two message tests red against the 53b67ec production source (old text), green on the diff. Guide item 4, CHANGELOG `### Refusals` bullet, the `periodic` section of the input reference (Codex luna sweep applied). Gates in a private debug tree (SCLS `debug`, gcc, Armadillo, non-MKL): build zero warnings, `test_homology` 53/53, `make check` 17/17, twice (diff as audited, and after the fix round). Juries: plan terra/high + grok-4.7/high (five confirmed amendments: no "only", within-set branch stays first, casts, if/else, guide rule stated as the predicate); code terra/high + grok-4.7/high (Grok's coil sentence wider than the predicate: fixed; two test pins; `treatment_label` as a switch). No deck run: no shipped periodic deck has ferro, coil or buffer blocks. Residue: none; the wider rule "a coil cannot end on the periodic plane" is Christian's statement, not the code's, and is left to him.

> **Scope guards (from Christian's brief, relayed by the Server session 2026-10-09):**
> - Exactly (a) message, (b) test, (c) documentation. **Do not change what is refused**, only how it is reported.
> - **Do not touch** the first-facet classification (`mMasterType`, `cl_InterfaceProcessor.cpp:100-107`) or the append order of shared copies — not ruled on.
> - No edits in the six closed cohomology-core units. No new register rows.
> - devel only; main stays pinned at v0.9.2.

---

## 1. Current Behaviour and How It Fails

Christian's ruling of 2026-10-08: a periodic pair whose two copies lie on interface sets of
different master **domain type** (not material) is refused. The code refuses it through
`InterfaceProcessor::pair_within` (`cl_InterfaceProcessor.cpp:747-754`): a sibling set is
admitted only when the treatment matches and, for tie-weight sets, the master type matches. The
sibling scan at `:489-503` skips an excluded set with `continue`, so `D` stays `nullptr` and the
error at `:506-509` reads

```
periodic partner %lu of interface node %lu was not duplicated in any interface set
```

| Failure | Mechanism | Evidence |
|---|---|---|
| Misleading message — domain type | partner's copy exists in a tie-weight set of another master type; the scan drops it silently | `:491` (`pair_within` false → `continue`), `:506` |
| Misleading message — treatment | partner's only copy is on a `Decouple` (coil-touching) set, or the node's own set is the decoupled one | `:751-752`, `:506` |
| Genuine message | no set at all holds the partner | fall-through from the `else` at `:487` to `:506` — must stay as is |

A user reading the message looks for a missing sideset. The real cause is geometric: the periodic
plane coincides with a change of domain type (conductor surface on one side, ferro / buffer /
air-master surface on the other) or a coil ends on the plane. Both are avoidable by moving the
plane into a homogeneous section.

**Bottom line:** the filter is right and the diagnosis is wrong; the fix records *why* a sibling
set was excluded and reports that.

## 2. Mechanism

In `pair_cross_set_periodic_duplicates`, the within-set branch at `:483-486` stays first and
untouched. Inside the existing `else` (`:487-504`), the sibling loop records the excluded sets
that hold a copy of `B`, by cause, instead of dropping them; with several, the last in `mSets`
order is the one named (deterministic):

```cpp
InterfaceSet * tOtherType = nullptr ;       // same treatment, other master type
InterfaceSet * tOtherTreatment = nullptr ;  // other treatment
if ( tDuplicates.key_exists( B->id() ) ) { D = tDuplicates( B->id() ); }   // unchanged
else
{
    for ( InterfaceSet * tOther : mSets )
    {
        if ( tOther == tSet || ! tOther->duplicate_map().key_exists( B->id() ) ) continue;
        if ( ! this->pair_within( tSet, tOther ) )
        {
            if ( tSet->treatment() == tOther->treatment() ) tOtherType = tOther ;
            else tOtherTreatment = tOther ;
            continue;
        }
        ... ambiguity check and D = candidate, unchanged ...
    }
}
if ( D == nullptr ) this->refuse_unpaired_periodic_copy( tSet, A, B, tOtherType, tOtherTreatment );
```

Format arguments: ids cast to `long unsigned int` as at `:508-509`; `to_string( DomainType )`
and the treatment label passed through `.c_str()`. Parameters carry the `a` prefix.

`refuse_unpaired_periodic_copy` (new private const member, declared in the header) holds three
`BELFEM_ERROR`s in this order — setup code, so the message survives a release build:

1. **domain type** (`tOtherType != nullptr`):
   `periodic partner %lu of interface node %lu is duplicated on an interface of another domain type ( node %lu: %s, node %lu: %s ): the periodic plane coincides with a change of domain type, which the periodic tie cannot represent; move the periodic plane into a section of one domain type`
   — `%s` are `to_string( master_type() )` of the node's own set and of the excluded set.
2. **treatment** (`tOtherTreatment != nullptr`):
   `interface node %lu and its periodic partner %lu are duplicated on interfaces of different treatment ( node %lu: %s, node %lu: %s ): the decoupled one touches a coil, and a coil cannot end on the periodic plane; move the periodic plane away from the coil`
   — `%s` are `tie-weight` / `decoupled`, from a file-local label function on `InterfaceTreatment`. Symmetric wording, because the loop reaches the pair from either side first.
3. **genuine** (neither): the existing text, unchanged.

When both excluded kinds hold the partner, the code reports the domain-type message first. No
fixture builds that case and the documentation makes no promise about it: without "only", each
message is truthful on its own, so the order carries no correctness claim (Codex, plan round 1).

Nothing else moves: `pair_within`, the within-set branch, the ambiguity check, the identity checks
at `:513-525`, `mMasterType`, the copy order.

## 3. Ordered Steps

- [x] **R1** `cl_InterfaceProcessor.{hpp,cpp}`: the recording loop and `refuse_unpaired_periodic_copy` as in §2. Update the comment at `:480-481`.
- [x] **R2** `tests/homology/test_InterfaceProcessor.cpp`:
  - new `Layout::Coil` — the `set_domain_type` expression for block 2 (`:114-115`) grows a Coil branch; sideset 12 is then a coil|air set, `Decouple`, master type Coil (`:26-33`: a coil on either side decouples; `:357-358`: one side air admits the set). The fixture hardcodes C2 as master of `tTop2` (`:136`), which is what `fix_facet_masters` would leave for a coil over air (`cl_MaxwellFactory.cpp:1501-1508` flips only when the master's type is the lower; `en_DomainType.hpp:23-31`).
  - `PeriodicPairAcrossDomainTypesIsRefused`: `Junction j( Layout::Buffer, true )`, links `{0,1},{1,0},{2,3},{3,2}`; a1 (id 3) in the buffer|air set only, a2 (id 4) in the conductor|air set only; expect `std::runtime_error` containing "change of domain type", both ids, "Buffer" and "Conductor".
  - `PeriodicPairAcrossTreatmentsIsRefused`: `Layout::Coil`, same links; expect "coil cannot end on the periodic plane", both ids, "tie-weight" and "decoupled".
  - `PeriodicPartnerOutsideEveryInterfaceSetIsRefused`: `Layout::SideSetPerBlock`, link `{2,4}` (a1 → q; q is in no set); expect the unchanged "was not duplicated in any interface set".
  - `PeriodicPairsAcrossFerroSetsAreTied`: `Layout::FerroFerro`, links as above; both ferro|air sets have master type Ferro and private copies, so a1–a2 pairs across sets; expect no throw and the two copies tied; backup holds one entry per copy pair (sets 11, 12, 13 each tie n–m, plus a1–a2 → 4).
  - `PeriodicPairsAreTiedOnce` (existing, `:459-464`) stays as the conductor same-type cross-set success case.
  - header comment: add the layout and the refusal cases.
  - Positive control: the new test TU linked against the 53b67ec production source (the gate worktree before the source diff is applied): the two message tests red (old text), green on the diff. The other two pass on both; they are guards.
- [x] **R3** documentation:
  - `src/homology/doc/homology_usage_guide.md`, item 4 of the pipeline list (`:1577-1586`, the shared-conductor-copy paragraph): the refusal stated as the code has it — a periodic pair is tied only between copies of interface sets of the same treatment and, for tie-weight sets, the same master domain type (`pair_within`); the motivating case in Christian's words (the tie identifies like unknowns, h on conductor edges and φ elsewhere, and cannot represent a formulation change across the plane) given as the reason, not as the predicate — the rule also refuses ferro against buffer, both φ blocks (`cl_Topology.cpp:428-441`); the remedy. Any formulation citation only if confirmed from `doc/literature_references.md`. Codex language sweep.
  - `CHANGELOG.md` `[Unreleased]` → `### Refusals`: one bullet for the refusal itself, which 0.9.2 did not have (its cross-set pass tied any sibling copy; `git show v0.9.2:src/homology/cl_InterfaceProcessor.cpp`), naming the message a user now sees.
  - `doc/input_file_reference.md` §`periodic` (`:838`): one sentence, worded so it reads as a property of any periodic model (`Mesh::has_periodicity()`, `:457-459`), not of this key. Not an input-key contract change, so `doc/input_schema.yaml` is untouched. The Maxwell usage guide does not describe `periodic` (grep: no hit), so nothing there.
- [x] **R4** gate: `test_homology` and `make check` in a private debug tree (SCLS `debug` flavor, gcc, Armadillo, non-MKL), the two refusal tests also run against the 53b67ec binary for the positive control. No deck run: no shipped periodic deck has ferro, coil or buffer blocks (the Server checked all `input.conf`), so the refusal is unreachable from the shipped decks, and the diff touches nothing beyond the message path.
- [x] **R5** close out: devlog, `devlog/README.md`, this file to `todo/closed/`, push devel, report SHA + file list to the Server session.

## 4. What could break

- A legal state the new guard must still accept: the within-set branch (`:483`) and admitted
  siblings are untouched, so every pair that paired before still pairs. The only behavior change
  is the text of an abort that already fired. [confidence: high — the three `BELFEM_ERROR`s sit
  exactly where the one did]
- `to_string( DomainType )` is declared in `en_DomainType.hpp` and returns the type name
  (`en_DomainType.cpp:21-49`); `cl_InterfaceProcessor.cpp` already names `DomainType`. [high]
- The fixture with a Coil block: `Topology::run` groups blocks by type and `create_interface_sets`
  admits coil|air (one side is air, `:347-351`); the conductor sweep in `relink_elements` skips
  the coil block. Should the fixture not build the coil case, the brief allows dropping that test
  — it would be reported, not silently skipped. [medium — to be shown by the run]

## 5. Open Design Questions

- O1 — order of the two mismatch messages when both excluded kinds hold the partner. RESOLVED
  2026-10-09 → domain type first in the code, no documented promise, no fixture; neither message
  over-claims once "only" is removed (Claude, after Codex plan round 1).

## 6. Audit Trail

- `tmp/ai_exchange/review_periodic_refusal_message.md` — plan jury and code jury, same slug.
