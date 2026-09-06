# Jury audit: is either clean_spfa() run on the gantry path skippable?

**Date:** 2026-08-25
**Topic:** `src/homology`: `Cohomology::clean_spfa()` executes twice per `CutFactory::compute_cohomologies()` (constructor clean + post-`updatekGeneratorsFromHomology` clean). Question: intentional, and is either run skippable without regression?
**Round:** three-AI jury (Claude pre-registered, Codex + Grok blind parallel), static only.
**Record:** `tmp/ai_exchange/review_spfa_double_run.md` (pre-registration, both audits, verification, reconciliation).

## Verdict (unanimous, source-trace evidence, no execution)

**Keep both runs.**

- **Run 2** (`cl_CutFactory.cpp:357`) is load-bearing. Step 7 of
  `updatekGeneratorsFromHomology` (`cl_Cohomology.cpp:1341-1350`) forms integer
  linear combinations whose coefficients are generically non-unit even from unit
  inputs; `CutData` encodes coefficients as two ±1 bitsets and drops everything
  else (`cl_CutData.cpp:241-255`, `cl_CutData.hpp:143-155`), so skipping run 2
  aborts 3D decks in `determine_cut_case_3d` and can silently mis-orient 2D
  release decks (debug-only assert at `cl_CutData.cpp:431-437` passes on 0+0).
  Run 2 is also the only `remove_cut_pockets(true)` pass over the final
  (post-recombination) generators.
- **Run 1** (constructor clean) is *conditionally* redundant for the final
  cohomology classes, but the condition is unproven: the invariance argument
  requires every suggested homology generator to be an exact cycle, and
  `Homology::suggest_Homology` builds them with boundary tracking suppressed
  (`Chain(1,mMesh,true,false)` → `mBoundary = nullptr`), never checking
  closedness; thin-shell generators are raw curve chains. Even under
  invariance, skipping run 1 changes the concrete cut representatives
  (path-dependent greedy + pocket removal) and moves the "no thin cut exists /
  mesh too coarse" certificate from raw to combined classes. Constructor clean
  is separately load-bearing for `tests/homology/test_CohomologyPeriodic.cpp`
  and the (currently dead) `!mSuggestHomologies` factory branches.

Whether masking early raw-class SPFA aborts would be acceptable is a product
decision, explicitly not settled by the jury.

## Incidental P2 findings (pre-existing, no action taken)

- `Cohomology(Mesh*, Mesh*)` field constructor always aborts in `clean_spfa`
  (`mSimplicialComplex` stays null); no in-tree caller.
- `clean_greedy` is dead: `mFunClean` is set to `clean_spfa` in all three
  constructors and nowhere else.
- `tTransInv` is `Matrix<real>` holding integers, implicitly narrowed to `int`
  at `cl_Cohomology.cpp:1347`.
- `CutFactory::create_sidesets_2d/3d` and `SideSetFactory` are unreferenced;
  the live cut consumer is `CutProcessor`/`CutData`.

If a run-1 skip is ever pursued: defer-flag in the constructor (suggest-homologies
path only), then an A/B on gantry plus one thin-shell deck comparing cut
sidesets, CutData cases, and solve behavior. Static review cannot close it.

## Follow-up in the same session (Christian's ask: document + todo currency)

- **Module docs:** the double-run verdict is now stated inline in
  `src/homology/doc/cohomology_algorithms.md` (new subsection "Cleaning runs
  twice on the factory path", plus the stale §7 line citations replaced with
  symbol anchors) and in the `clean()` entry of
  `src/homology/doc/homology_usage_guide.md`. Codex language sweep applied.
- **Todo currency audit of the cohomology module:**
  `thin_cut_nonunit_rectification_implementation.md` was stale against the
  debt register (its "Phase 3 only / T8 part done / T9 open" state predates the
  layer suite of 2026-08-13, the ctest gate of 2026-08-23, and the
  forced-repair fixture + T9 corc discharge of 2026-08-24 that struck DR-23);
  a dated currency note now marks Phase 3 CLOSED, Phase 4 (min-cost L1) the
  only open item; the matching `todo/README.md` entry, frozen at its 2026-07-01
  "Status: PLAN" text, got the same dated update. `closed/cut_pocket_removal_rules.md`
  and the register rows (DR-23 struck, DR-71) were already current.
- **DR-107 registered** (P3 bundle): the round's four hygiene residues:
  dead `clean_greedy`, always-aborting field constructor, unreferenced
  `create_sidesets_2d/3d`/`SideSetFactory`, `tTransInv` real-to-int narrowing.
