# DR-107 homology hygiene cleanup: five dead-code items, double jury round

**Date:** 2026-08-25
**Topic:** `src/homology`: execute DR-107 (residues of the same day's double-`clean_spfa`
jury) on Christian's condition "ask Grok and Codex if they agree, and if so, we can
clean up". Full record: `tmp/ai_exchange/dr107_cleanup_plan.md` + `dr107_cleanup_impl.diff`.

## What landed (uncommitted; net -1201/+135 lines)

- **D1:** `Cohomology::clean_greedy()` (123 lines) deleted along with the whole
  `mFunClean` function-pointer dispatch; `clean()` calls `clean_spfa()` directly.
  The greedy ALGORITHM stays untouched: production rectification is
  `rectify_greedy_sweeps` inside `clean_spfa()`. Six provenance comments reworded
  to state the quotient-fold / plateau-seed semantics without naming the deleted
  member. Both auditors: the pointer is a 2026-07-15 temporary switch, not a
  Phase-4 extension point.
- **D2:** field constructors `Cohomology(Mesh*, Mesh*)` (always aborted: null
  complex meets `clean_spfa`'s guard) and `Homology(Mesh*, Mesh*)` (left
  `mSimplicialComplex` indeterminate, no initializer) deleted with both
  `generatorsFromField` members. Zero callers tree-wide.
- **D3:** `CutFactory::create_sidesets_2d/3d` (236 lines), the write-only
  `mCohomologies` member, and the entire unreferenced `SideSetFactory` TU
  (542 lines) deleted; CMakeLists entry dropped.
- **D4:** `tTransInv` in `updatekGeneratorsFromHomology` is now `Matrix<int>`,
  ending the silent real-to-int narrowing at `addCochainToCochain`. Grok's
  sharpening: `rowExchange` already truncated through `int` under
  `Matrix<real>`, so the int container matches the ring the kernel always
  computed in; `kernelImage<int>` is the production path for `mD`.
- **D5:** dead `updatekGeneratorsFromHomology_old` (101 lines) deleted; its
  free-cut basis-choice caveat lives on as a comment at the live `kernelImage`
  call.
- **Docs (jury-mandated):** usage-guide constructor sections and three examples
  rewritten to the live `SimplicialComplex` path; the false "CutFactory uses
  SideSetFactory during run()" section replaced by the real
  CutProcessor/CutData path with a retirement note; file-tree line dropped;
  "call exactly once" note corrected (the constructor computes the groups).

## Round mechanics

Plan jury: Codex leg died on a vendor usage limit (reported immediately per the
standing rule; Christian added credits; identical blind rerun). Both vendors
then agreed on all five items, D2/D3 conditional on the same-session doc
rewrites, plus `mCohomologies`. Code jury on the implementation diff caught
four defects in the round's own doc edits, all fixed and re-verified in-round:
an orphaned `delete tFlagged;`, `get_Generators()[k]` (9 sites) and
`gen->data()` (2 sites) API rot (Cell has no `operator[]`), a stale
`cl_CutFactory.cpp:440-511` citation in `cohomology_algorithms.md` §8, and a
BeltedTree-constructor error in the rewritten method blurb.

## Gates

- Scoped grep gate (src/ + tests/): zero references to every deleted symbol
  outside the deliberate retirement note. Dated records (devlogs, register,
  lessons evidence, closed todos) keep historical names by design.
- `-fsyntax-only` with the production flags (`-Og -Wall -Werror
  -pedantic-errors`, full define/include set from `flags.make`): clean on
  `cl_Cohomology.cpp`, `cl_Homology.cpp`, `cl_CutFactory.cpp`.
- **OWED (Christian): `make check` / `check-fast`, homology suite 9 fixtures**,
  the real executable gate, especially for D4; then commit + DR-107 strike.

## Routed to Christian (his pre-existing uncommitted hunks, NOT this round)

The blind code jury reviewed the whole `git diff HEAD` and flagged two real
defects in the in-flight Progressbar work that predates this session:
1. `clean_spfa`'s `Progressbar` hides the cursor (`\033[?25l` in `reset()`)
   and only `finish()` restores it; on the coarse-mesh `BELFEM_ERROR` abort
   (throw in debug, abort in release) the cursor stays hidden in the
   operator's terminal, and on a constructor-body throw the bar leaks.
2. The dense-fallback "greedy rectification stalled" WARNING moved from
   `InfoLevel::Default` to `Verbose`, so default runs would never see the
   never-observed-fallback fire.

## Residue (recorded, deliberately untouched)

`create_cut_sideset_2d/3d`, `save_edges`, `coefficient_TMatrix`: equally dead
siblings, next hygiene bundle if wanted. Example 3's `create_flagged_submesh`
names an API that does not exist in `src/` (pre-existing guide fiction).
