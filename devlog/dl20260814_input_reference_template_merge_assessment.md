# Devlog 2026-08-14 — BELFEM_Templates README vs. Input Reference: Merge Assessment

**Date:** 2026-08-14
**Topic:** Three-AI comparison of Gregory's `tmp/BELFEM_Templates/README.md` (913-line user guide) against `doc/input_file_reference.md` / `doc/input_schema.yaml`, to decide what to merge into the two contract documents. Scope narrowed by Christian mid-session: the example decks are Gregory's to merge into `./examples`; this task concerns only the two doc files.
**AIs involved:** Claude (pre-registration + verification), Codex (audit), Grok (refutation)
**Claude Confidence:** high (all auditor citations re-verified against the tree)
**Codex Audit Confidence:** high
**Grok Audit Confidence:** high (~90%) on parser facts, medium (~75%) on completeness of the pedagogical inventory
**Literature References:** N/A (documentation contract; Messe et al. 2023 touched only where the README's solver tip contradicts it)
**Verification:** parser facts by static source trace (reviewed); the applied doc edits gated executably — `python/belfem-conf drift` clean (125/125 anchors, including the new `duplicate_sections` anchor) and `python/belfem-conf check` ok on all six example decks, both re-run after the Codex prose pass. No solver deck was executed.

## Summary

The README contributes **no new input key, value, alias, unit, or enum** — every token is already in the schema, usually with a larger legal set (3/3 independent walks). Its value is pedagogical, plus **one genuine contract fact missing from both artifacts**: unlabeled same-type sub-sections are all live. The README also contains a set of errors that must not be merged; every long snippet in it carries at least one, so the merge rule is **rewrite from parser truth, never paste**.

## Key Findings

- **Duplicate-section contract (undocumented anywhere):** `cl_Input_Section.cpp:76-89` stores every child section in the ordered `Cell mData` (duplicates preserved) and additionally in a name `Map` (key `type` or `type:label`, later unlabeled duplicate overwrites). The Maxwell/thermal BC factories, topology, materials, layers, and circuit readers all walk `mData` by index — so each repeated unlabeled `current { }` becomes an independent BC (`cl_MaxwellBoundaryConditionFactory.cpp:27-35,150-169`). `layers : tapeN` multiplicity works by a **different** mechanism (label-keyed map, `cl_MaxwellFactory.cpp:2826-2872`). Open risk recorded: a future consumer using `section_exists()/section("type")` sees only the *last* unlabeled duplicate — the schema should state which access pattern a parent uses.
- **README errors confirmed (none where the README is right and the reference wrong, 2/2):** units "mm, m" (any length unit of the right dimension works); solver libraries "mumps/strumpack/petsc" (actual: +umfpack, superlu, pardiso); builtins list of four (actual: +In, Pb, Sn, iron, magnesia, alloy formulas); "only generalized pellikka" (four algorithms); curves `@` glossed "edge_id @ sideset_id" (actually sideset ∩ sideset, `cl_CurveFactory.hpp:58`; the 2-D parser silently skips `@`-tokens, `cl_Input_Section.cpp:608-612`; the 3-arg `get_ids(key, tape, ids)` overload at `:520-553` is dead — no caller); inert `label :` keys inside topology domains; `coil { material }` inert (`cl_FEM_Domain.cpp:29-36`, Codex+Grok); "homology required for all EM simulations" (conditional, `cl_MaxwellFactory.cpp:922-935`, Grok); the range-notation paragraph teaches `1:4` as one BC's id list — re-teaching exactly the missing-bracket failure mode the reference warns about (Grok R1, the highest-hazard item); "mumps small / strumpack large" tip contradicts the published order.

## Changes Made

Merge set approved by Christian ("make it so") and applied, then reworded per the standing Codex prose pass (all five suggested replacements taken):

**`doc/input_file_reference.md`:**
1. §1: replace "two sections of the same type are distinguished by label" with the full two-store contract: ordered duplicates all live for index-walking consumers; the name map returns only the last unlabeled duplicate.
2. §9: state that repeated same-type BC sub-sections are legal and each becomes its own condition (with the generator-count interaction); one sentence for the 2-D rule (absent output list → input reused, `cl_MaxwellBoundaryConditionFactory.cpp:123-128`); one sentence of bearing rationale (fixes the potential to remove the null space).
3. §1 lists/groups: add the intent sentence — brackets group ids that are electrically one terminal.
4. §6: one-line cross-link that multiple labeled `layers : tapeN` blocks resolve per-thinshell by label.

**`doc/input_schema.yaml`:**
5. `grammar`: add `duplicate_sections` field carrying the two-store contract and the per-parent access pattern (index walk vs name lookup).

**Not merged:** all README example-catalog/tips/workflow/contact content (out of scope, Gregory's `./examples` merge); every erroneous item above. The ten confirmed corrections were written up for Gregory as `tmp/BELFEM_Templates/README_errata.md` (in his repository clone, so it travels with his merge; decks unaffected, prose only).

**Gates:** `belfem-conf drift` clean, 125/125 anchors (the two wrong-consumer notes pre-exist this session); `belfem-conf check` ok on all six example decks — the validator already consumes the new `duplicate_sections` field (informational note on the doubled `terminal pair` in `examples/circuit`).

## Open Questions

- None for this task. Gregory owns applying the errata to his README before/while merging the examples.

## Files Updated

- doc/input_file_reference.md (§1 two additions, §6 one, §9 two)
- doc/input_schema.yaml (`grammar.duplicate_sections`)
- tmp/BELFEM_Templates/README_errata.md (new, for Gregory)
- tmp/ai_exchange/input_reference_template_merge.md (ephemeral exchange; distilled here)
- devlog/dl20260814_input_reference_template_merge_assessment.md (this file)
- devlog/README.md (index)
