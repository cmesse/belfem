# The cross-set periodic refusal in InterfaceProcessor names its cause: a change of master domain type, or a tie-weight copy against a decoupled coil copy, with both node ids

**Date:** 2026-10-09
**Purpose:** Christian's task, relayed by the Server session: the refusal that `pair_within`
produces since c1d2bf7 (DR-198) for a periodic pair whose two interface copies lie on sets of
different master domain type or treatment said "partner was not duplicated in any interface
set", which sends a user looking for a missing sideset. Make it say the real cause, pin all
three messages with tests, document the rule, push to devel. What is refused does not change.
**Module:** `src/homology` (`cl_InterfaceProcessor`), `tests/homology`, `src/homology/doc`,
`doc/input_file_reference.md`, `CHANGELOG.md`
**AIs:** Claude Fable 5.1 (DrClean session: plan, implementation, juries, gates), Codex
`gpt-5.6-terra`/high and Grok `grok-4.7`/high (plan jury and code jury), Codex
`gpt-5.6-luna`/medium (language sweep of the two user-facing passages); the Server session
(brief, blind xhigh delta audit and the 4/4 pipeline on the pushed SHA)
**Status:** COMPLETE — pushed to devel. Gates in a private debug tree (SCLS `debug` flavor, gcc,
Armadillo, non-MKL), run twice (the diff as audited, and after the fix round): build with zero
warnings, `test_homology` 53/53 (49 before), `make check` 17/17. Positive control: the two
message tests red against the 53b67ec production source with the old text. No deck run — no
shipped periodic deck has a ferro, coil or buffer block (the Server checked every `input.conf`),
so the refusal is unreachable from the shipped decks and the diff touches nothing beyond the
message path
**Exchange:** `tmp/ai_exchange/review_periodic_refusal_message.md`,
`tmp/ai_exchange/sweep_periodic_refusal_docs.md`
**Plan:** `todo/closed/periodic_refusal_message.md`

## Summary

- **Mechanism.** The sibling scan in `pair_cross_set_periodic_duplicates` used to `continue`
  past a set that `pair_within` excludes, so `D` stayed null and the one error fired. The scan
  now skips sets that do not hold the partner's copy first, then diverts an excluded set into
  `tOtherType` (same treatment, other master type) or `tOtherTreatment` (other treatment) and
  continues; admitted siblings keep the ambiguity check. A null `D` calls the new private
  `refuse_unpaired_periodic_copy`, three `BELFEM_ERROR`s in order: domain type (both ids and
  both `to_string( master_type() )`), treatment (both ids and `tie-weight` / `decoupled`), and
  the unchanged "not duplicated in any interface set" for a partner on no set at all. The
  within-set branch, `pair_within`, the identity checks, `mMasterType` and the copy order are
  untouched. With several excluded sets the last in `mSets` order is named; with both kinds
  the domain-type message wins — deterministic, undocumented, and harmless because neither
  message claims to be the only cause.
- **Tests** (`tests/homology/test_InterfaceProcessor.cpp`): `Layout::Coil` (block 2 Coil, so
  sideset 12 is a decoupled coil|air set); `refusal()` helper; four new tests —
  `PeriodicPairAcrossDomainTypesIsRefused` (`Layout::Buffer`, pins `node 3: Buffer`,
  `node 4: Conductor`), `PeriodicPairAcrossTreatmentsIsRefused` (`Layout::Coil`, pins
  `node 3: tie-weight`, `node 4: decoupled`), `PeriodicPartnerOutsideEveryInterfaceSetIsRefused`
  (a1 → q, the old sentence verbatim), `PeriodicPairsAcrossFerroSetsAreTied` (two ferro|air
  sets of the same master type pair a1–a2; four backup entries). The existing
  `PeriodicPairsAreTiedOnce` remains the conductor same-type cross-set success case.
- **Documentation.** Guide item 4 of the InterfaceProcessor pipeline states the rule as the
  code has it (same treatment and, for tie-weight sets, the same master domain type), gives
  the h-on-edges / φ-elsewhere identification as the reason (Messe et al. 2023, §1), and says
  the rule also refuses ferro against buffer, two φ blocks. CHANGELOG `### Refusals` carries
  the refusal itself, which 0.9.2 lacked (its cross-set pass tied any sibling copy). The
  `periodic` section of the input reference carries one paragraph, worded for any periodic
  model, not for the key.

## Key Findings

- **Plan jury (Codex terra/high, Grok 4.7/high), five confirmed amendments.** "only" in the
  domain-type sentence was false when both excluded kinds hold the partner — dropped, which
  also retired the untested ordering promise; the plan's snippet had omitted the within-set
  branch that must run first (Grok's FerroFerro trace: a scan from the top would hit the
  ambiguity error); casts and `.c_str()`; if/else over a ternary assignment; the guide's
  h-vs-φ reason is wider than the predicate (`cl_Topology.cpp:428-441` puts Buffer and Ferro
  both on the φ side, and `pair_within` still refuses them against each other) — the guide now
  states the predicate and gives the reason as the motivating case, with the citation.
- **Compile gate caught `[[noreturn]]`.** `belfem_assert` is not `noreturn` (the reaction may
  throw), so GCC under `-Werror` rejects a `noreturn` function that ends in `BELFEM_ERROR`.
  Found by compiling the TU with the module's real `flags.make`; the attribute was removed.
- **Code jury (Codex terra/high, Grok 4.7/high).** Codex: clean, one minor on test
  orientation — refuted, the format binds each id to its own set's label, so the fragments
  hold from either visit order (Grok traced the same). Grok, confirmed and fixed: the coil
  sentence "a coil cannot end on the periodic plane" stated a rule wider than the predicate —
  `pair_within` admits any two decoupled sets, and a coil|air pair inside one set is tied with
  no treatment check; the abort is the mixed case only. The message now says the decoupled
  copy carries no unknown so the tie-weight copy has nothing to tie to, remedy "move the
  periodic plane away from the coil"; the input reference says a coil must not touch the
  interface on one side of the plane only. Also fixed: the domain-type test now pins the type
  beside its id, and `treatment_label` is a `switch` without `default`, so a third enumerator
  is a compile error.
- **The wider rule is Christian's, not the code's.** "A coil cannot end on the periodic plane"
  came with the brief. The code accepts a coil that touches the plane on both sides. Nothing
  here changes what is refused; whether it should is his ruling.

## Changes Made

- `src/homology/cl_InterfaceProcessor.hpp` — declaration of `refuse_unpaired_periodic_copy`
- `src/homology/cl_InterfaceProcessor.cpp` — the recording scan, the method, `treatment_label`
- `tests/homology/test_InterfaceProcessor.cpp` — `Layout::Coil`, `refusal()`, four tests
- `src/homology/doc/homology_usage_guide.md` — item 4 of the pipeline list
- `CHANGELOG.md` — `### Refusals` bullet
- `doc/input_file_reference.md` — `periodic` section
- `todo/closed/periodic_refusal_message.md`, `todo/README.md`, this devlog, `devlog/README.md`

## Open Questions

- None in scope. The wider coil rule above is Christian's to rule on; a periodic deck that
  reaches the refusal does not exist among the shipped decks, so the hand-built fixture is the
  only gate for the message path.
