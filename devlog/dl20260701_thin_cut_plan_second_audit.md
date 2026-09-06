# Devlog 2026-07-01 — Thin-Cut Rectification Plan: Second Audit Round

**Date:** 2026-07-01
**Topic:** Second Codex + Grok audit round for the rectify-or-certify plan
(`todo/thin_cut_nonunit_rectification_implementation.md`); resolution of the periodic-quotient
(QC) and null-complex open questions; plan updated.
**AIs involved:** Claude (primary), Codex (Q1–Q4 audit + prose polish), Grok (Q2–Q3 audit)
**Claude Confidence:** high on all adopted items (every load-bearing citation re-verified)
**Codex Audit Confidence:** medium overall, high on Q1/Q2/F1/F2
**Grok Confidence:** high on Q2 (~90%), medium on Q3 no-feasible-cycle (~70–75%)
**Literature References:** none new this session (theory basis unchanged:
`src/homology/doc/thin_cut_nonunit_rectification.md` and its reference table)

## Summary

Session goal: fold the 2026-07-01 discussion findings into the implementation plan, with
Codex and Grok double-checking the medium/low-confidence claims. Exchange thread:
`tmp/ai_exchange/thin_cut_plan_update_audit.md`. Outcome: two of the plan's open risks
(periodic quotient QC, null `mSimplicialComplex`) are resolved by static analysis; two new
settled design decisions (shared solver core; `aPeriodicity == true` precondition); one
reproducer correction (corc is unit-coefficient); one question deliberately left open
(greedy termination on feasible inputs).

## Key Findings

- **Quotient adjacency equivalence confirmed (3-way, high).** The plan's SPFA rule
  (master∪slave incident edges, `original_edges()` filter, `rep()` fold) reproduces exactly
  the quotient coboundary of `SimplicialComplex::create_complex()`
  (`cl_SimplicialComplex.cpp:144-188`) and the edge set `clean()` fires
  (`cl_Cohomology.cpp:244-301`); sign conventions match. The `!is_flagged()` skip at
  `cl_Cohomology.cpp:277` is redundant given the domain filter — slave seam edges never enter
  `original_edges()` because slave entities are unflagged before complex construction
  (`cl_SimplicialComplex.cpp:116-131`, bitset populated at `:468-471`).
- **No mirroring pass needed in the `c − dθ` write-back (high).** Downstream never reads the
  cochain by slave edge keys: `CutData::collect_coefficients()` (`cl_CutData.cpp:230-261`) and
  `collect_edges()` propagate ±1 bits / edge pointers to `edge->periodic()` themselves; all
  later reads go through the `weight()` bitsets. Residual scope limit (pre-existing): the
  same-sign propagation is documented correct for translational periodic maps only.
- **`aPeriodicity == true` is a hard precondition (both auditors, adopted).** A `false`-built
  complex on a periodic mesh would let slave seam edges into `original_edges()`, desyncing the
  domain from the solver's quotient fold. Production is safe (`cl_CutFactory.cpp:430`).
- **Null-complex risk closed (high).** The field constructor `Cohomology(Mesh*, Mesh*)`
  (`cl_Cohomology.cpp:23-30`) has zero call sites in `src/` and `nonfree/` and would crash
  today (`mSimplicialComplex == nullptr` dereferenced at `clean():199`). Decision:
  `BELFEM_ERROR` guard, no fallback edge domain.
- **corc reproducer correction (Codex, verified).** corc's generators are unit-coefficient
  (`todo/closed/periodic_thin_cut_continuity_fix.md`, discharged findings) — it is a periodic
  regression check, not a non-unit test case; the Phase 0 verdict needs a real coarse-CCT
  input.
- **Greedy termination on feasible inputs — open (auditor split).** Grok: no proof exists;
  exact-rule simulation found 0 feasible cycles in ~1700 instances (transient coefficient
  growth to |c| = 7 — independently supports int64 distances). Codex: monotone-relaxation
  heuristic argues termination. Agreed consequence either way: a `clean()` hang does **not**
  certify Regime 2 — only the SPFA verdict does. Moot for production once Phase 2 lands.

## Doc extension (same session, user-approved)

`src/homology/doc/thin_cut_nonunit_rectification.md` gained §1.1–1.3 answering Christian's
theory questions on the geometry of the non-unit pushed object: (1.1) the push cannot emit a
second facet copy — both sheets land on the same conjugated face, giving one facet with
multiplicity 2 and a non-uniform jump (`2I` in the throat, `I` elsewhere), which is the
structural reason the pipeline is unit-only (node instances multiply, not facet instances);
(1.2) the multi-sheet cut is watertight as an integer 2-chain (`dc = 0` *is* the
watertightness condition, with multiplicities), ending along a branch line where sheets
separate — not an embedded 2-manifold; (1.3) resolvability per regime — Regime 1 rectifies to
a unit cut (nothing multi-sheet survives), Regime 2 resolves only as a multiplicity chain via
`q_k = p + k·I` (topologically exact, locally blind in the throat — the accuracy argument for
refinement-first). Codex prose+sanity pass applied; its technical flags adopted: diagonal-quad
qualifier in §1.1, explicit `dc = 0` balance at the branch line, Haken condition scoped
per-tet, "accuracy argument" wording.

## Decision + step tracker (same session)

Christian decided to **focus on the feasible-regime fix first** — replacing greedy `clean()`
with the global solve (multi-level duplication stays declined; refinement handles Regime 2
manually via the certificate). The plan gained an **ordered step tracker T1–T9**
(`meshfile_refactor_plan.md` §4 style): T1 shared `feasibility_solve()` core → T2 Phase-0
read-only diagnostic → T3 decision gate on the QB verdict → T4 `clean()` body replacement
(write-back) → T5 global post-rectification guards → T6 infeasible-branch messaging → T7
`determine_cut_case_3d()` error-message split → T8 unit tests (feasible / infeasible 3-cycle /
periodic) → T9 end-to-end validation (coarse-CCT + corc regression). Codex prose pass applied;
its scope catches adopted (T3 now reads "T4–T6 still land"; intro no longer promises a
per-step phase anchor that T7 lacked).

## Risk-list review (same session)

Christian asked whether any remaining Risks / open questions could be resolved. Claude
proposed four closures; Codex audited all code-grounded claims. Outcome — **three closed,
one narrowed after a refutation**:

- **Domain scope — closed by design** (T1 enqueues all representatives and relaxes every
  `original_edges()` arc; T4 writes back over the same domain). Reopens only with the Phase-4
  support-localization optimization.
- **Higher-order elements — closed for supported orders.** Cochain lives on the corner-edge
  skeleton only (`flag_corner_nodes` `cl_CutFactory.cpp:397,422`; `edge_key()` from corner
  endpoints `cl_EdgeFactory.cpp:323`; midnode in slot 2, `grab_nodes` `:477-508`). Codex
  caveat adopted: `EdgeFactory` rejects order ≥ 3 (`cl_EdgeFactory.cpp:272`), so the claim is
  scoped to orders 1–2 rather than "any order".
- **Numeric width — closed by design** (int64 distances in T1; feasible write-back fits
  `int`; Grok's transient `|c| = 7` growth confirms wide storage is required).
- **Cochain invariant consistency — Claude's closure REFUTED by Codex, verified, item
  narrowed instead.** `updatekGeneratorsFromHomology()` runs after the constructor's
  `clean()` on the `mSuggestHomologies` path (`cl_CutFactory.cpp:516`) and recombines the
  cleaned generators as *sources* of `addCochainToCochain` (`cl_Cohomology.cpp:586`), which
  reads their boundary/coboundary side-state (`cl_Cochain.hpp:294,299` — the
  `!mIsBound`/`!mIsCobound` guards are live for generators, checked down to the constructor
  flags). Consequence promoted into T4: side-state maintenance is now **required**, and the
  risk item closes at T4 review.

Large-mesh cost stays open (measured at T9); QC stays open (runtime PBC verification, T8c/T9);
greedy-termination theory question unchanged.

## Changes Made / Proposed

- `todo/thin_cut_nonunit_rectification_implementation.md` — updated (user-approved): two new
  design decisions (shared `feasibility_solve()` core; periodicity precondition), Phase 0
  reproducer + hang-inference caveats, Phase 1 equivalence note, Phase 2 infeasible-branch
  wording (release-active error, edge-ID certificate, manual-refinement first iteration),
  Phase 3 validation split (CCT vs corc), QC risk rewritten as largely-resolved, null-complex
  risk ticked resolved, new open-theory risk item, second audit-provenance block. Codex prose
  polish applied per standing rule.
- No source edits (plan remains PLAN status; edits require explicit approval).

## Open Questions

- QB (Phase 0): feasible-but-greedy-failed vs. genuinely infeasible on a real coarse-CCT
  generator — still the pivotal empirical unknown.
- Greedy termination on feasible inputs (theory; practically moot after Phase 2).
- Runtime QC verification on a PBC mesh (Phase 3 unit test).
- Higher-order elements: cochain assumed corner-edge-only (unchanged from first round).

## Files Updated

- todo/thin_cut_nonunit_rectification_implementation.md (audit adoptions + T1–T9 tracker)
- todo/README.md (entry updated)
- src/homology/doc/thin_cut_nonunit_rectification.md (new §1.1–1.3)
- devlog/dl20260701_thin_cut_plan_second_audit.md (new, this file)
- devlog/README.md
- tmp/ai_exchange/thin_cut_plan_update_audit.md (ephemeral thread, sweep-eligible)
