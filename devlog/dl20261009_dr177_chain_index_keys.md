# DR-177: Chain keys are simplex indices, the closed core retyped by Christian

**Date:** 2026-10-09
**Purpose:** Settle DR-177 — `Chain::getSimplicesMap` declared with `id_t` keys over an
`index_t` member — which had waited since 2026-09-22 for a ruling on what the map keys are.
Christian asked whether he could make that call himself instead of routing it to Gregory.
**Module:** `src/homology` (`cl_Chain`, `cl_Cochain` — closed core; `cl_BeltedTree`;
`tests/homology`; `doc/homology_usage_guide.md`)
**AIs involved:** Claude Fable 5.1 (primary, static review and the open-file edits); Codex
`gpt-5.6-luna`/medium (language sweep of the guide's Chain and Cochain sections). No jury: the
core edit is Christian's own hand, the test edit is a mechanical rename, the guide is prose.
**Status:** COMPLETE — DR-177 struck; live 3, archive 191
**Verification:** Christian's build and the homology tests green on the default build (his
report, 2026-10-09). The INT64 configure the row named as its gate did not run: no SCLS flavor
exists for it. Claude's part is reviewed, not verified — every claim below is from reading the
headers and bodies.
**Exchange:** `tmp/ai_exchange/dr177_chain_guide_sweep.md` — the sweep only

## Summary

The ruling needed no judgment call: the code already answers it. All 60 sites that fill a
`Chain` or `Cochain` pass an entity `->index()`, none passes `->id()`, and `Cochain` was
`index_t` end to end. `Chain` was the odd one out — its member was `index_t` but its getter and
all six setters took `id_t`. Christian retyped the two Chain files himself (the core is closed
to AI edits, protocol §7.1) and in the same pass renamed the camelCase chain and cochain
accessors to snake_case. Claude renamed the 32 test call sites, rewrote the guide's Chain and
Cochain sections from the headers, and closed the register row.

## The ruling and its evidence

- **Keys are simplex indices.** Fill sites: `cl_SimplicialComplex.cpp` (nodes, edges, faces,
  elements, all `->index()`), `cl_Homology.cpp`, `cl_Cohomology.cpp`, `cl_BeltedTree.cpp:252,284`.
  A grep for `id()` fed into any chain or cochain setter across `src/` and `tests/` returns
  nothing. Confidence high.
- **The row understated the defect.** It said "a one-word signature change". The two Chain
  files held 15 `id_t` tokens: the getter's declaration and inline body, six setters with their
  declarations and inline or out-of-line bodies, and one local in the chain-pairing operator.
  Under `BELFEM_INT64` only the getter fails to compile; the six setters would have narrowed a
  64-bit index into a 32-bit parameter silently, which is the worse failure. Confidence high.
- **The gate in the row cannot run here.** `-DUSE_MKL_64BIT_API=ON` needs a 64-bit SCLS flavor
  and none exists. The honest gate is the default build plus the homology tests, which ran on
  Christian's tree.

## What landed

**Christian, closed core.** Every `Chain` signature takes `index_t`, matching `Cochain`. The
first pass left one `id_t` behind in the inline body of `get_coefficient` (`cl_Chain.hpp:166`),
DR-177 in a new spot; the sweep caught it and he fixed it. Renames, all in the core files:
`get_simplices_map`, `get_coefficient`, `set_coefficient`, `get_dim`, `get_boundary`,
`get_coboundary`, `add_simplex_to_chain`, `remove_simplex_from_chain`, `add_chain_to_chain`,
`remove_chain_from_chain`, and the cochain counterparts `add_simplex_to_cochain`,
`remove_simplex_from_cochain`, `add_cochain_to_cochain`, `remove_cochain_from_cochain`.
In `cl_BeltedTree`, `create_TreeField` and `create_cohomologyField` became `create_tree_field`
and `create_cohomology_field` with `a`-prefixed parameters.

**Claude, open files.** 32 call sites in `test_CutData.cpp` (19), `test_Cohomology.cpp` (11)
and `test_CohomologyPeriodic.cpp` (2) follow the rename, plus one comment. The guide's
"### Chain" and "### Cochain" sections described an API that never existed — a member `mData`,
accessors `data()` and `add()`, an `inner_product`, an `evaluate`, value-returning operators.
Both sections now list the real public API and state what a caller needs: keys are indices,
every operator works in place, the call operator is the pairing or the evaluation, the
constructor flags decide whether the owned boundary and coboundary chains are allocated,
`get_coefficient` returns 0 for an absent index and a zero coefficient erases the entry,
`add_chain_to_chain` refuses a dimension mismatch and propagates into the owned boundary and
coboundary while `remove_chain_from_chain` does not, `Cochain` has no scalar multiply. Each
claim was read against the bodies before it went in. Codex returned eight prose-only edits and
no suspected errors; all applied.

**Left camelCase, by choice.** The ten Pellikka operators in `SimplicialComplex` (`pReduce`,
`pCoreduce`, `pCombine`, `pCocombine`, the two `pGeneralized…`, `reduceOmit`, `coreduceOmit`,
the two `createMatrixFrom…Map`), the six group methods in `Homology` and `Cohomology`
(`quotientGroup`, `…GroupOfChainComplex`, `generatorsOf…`, `updatekGeneratorsFromHomology`),
`Cochain::isBound`/`isCobound` (no callers anywhere), and `fn_Smith::smallestNonzero`
(file-local). The `p`-prefixed names are the paper's vocabulary and protocol §7.1 pins its
attested case to `coreduceOmit` and `pCoreduce` by name; four documents cite them.

## The shared-checkout collision

While this session's guide rewrite sat uncommitted, a peer session committed the interface
processor work as e30bd98 and `git add` took the whole guide file, rewrite included. The peer
noticed and restored those sections in 337a318 with a note that the rewrite stays in the
working tree for its author. It lands in this commit with its source. No content was lost;
the sequence is recorded because it is the third instance of the staged-window hazard the
session-close skill describes.

## For Gregory

An FYI rather than a question, since he owns the module: DR-177 is closed on the evidence
above, by Christian's hand. The 2026-09-22 report
(`dl20260922_homology_report_for_gregory.md` §3) called it one word; it was fifteen, all
type-only, no algorithmic content touched. The renames leave the Pellikka operator names as
they were.

## Files Updated

- `src/homology/cl_Chain.{hpp,cpp}`, `cl_Cochain.{hpp,cpp}`, `cl_SimplicialComplex.{hpp,cpp}`,
  `cl_Homology.cpp`, `cl_Cohomology.cpp`, `cl_CutData.cpp`, `cl_CutFactory.cpp`,
  `cl_BeltedTree.{hpp,cpp}` — Christian
- `tests/homology/test_CutData.cpp`, `test_Cohomology.cpp`, `test_CohomologyPeriodic.cpp`
- `src/homology/doc/homology_usage_guide.md` — Chain and Cochain sections, three method names
- `todo/debt_register.md`, `todo/debt_register_closed.md` — DR-177 struck
- `devlog/campaigns/release_1.0.md`, `devlog/README.md`
