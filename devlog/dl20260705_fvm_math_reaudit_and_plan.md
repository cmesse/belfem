# Devlog 2026-07-05 — FVM MPFA-O Math Re-Audit + Implementation Plan

**Date:** 2026-07-05
**Topic:** Independent re-check of the `src/fvm` MPFA-O math, audit of the 2026-06-26 defect
record (`todo/fvm_module_next_steps.md`), and drafting of the execution plan for the FVM thermal
thin-shell implementation week (`todo/fvm_implementation_plan.md`).
**AIs involved:** Claude (re-audit + plan), Codex (independent verification)
**Claude Confidence:** high (all findings carry counterexamples or direct literature quotes)
**Codex Audit Confidence:** high (all four claims confirmed)
**Literature References:** Agélas & Masson 2008 (F1) §5; Aavatsmark 2002 (F0); Klausen & Winther
2006 (F2); Alves et al. 2022 (paper6) Table 1

## Summary

Christian plans to replace/augment the nodal FEM thermal thin shell with a cell-centered FVM
(MPFA-O) solver — one constant `T` per cell mirroring the constant current density per element of
the h-φ thin shell. Task: (1) re-check the math of the current `src/fvm` skeleton, (2) audit the
2026-06-26 three-AI defect record, (3) write the implementation plan for next week.

The re-audit **retracts one prior three-AI defect, adds two new ones, and proves the module's
headline conservation test blind to its worst bug**. Net defect count unchanged (D1, D3–D7 real);
the fixes differ substantially from the 2026-06-26 plan.

## Key Findings

- **D2 retracted (FALSE POSITIVE).** The 1/3–2/3 continuity points (`cl_FVM_Factory.cpp:449,457`)
  + `0.5·L` sub-edge measure (`:454,464`) + centroid cell point (`:716`) are exactly the
  **Agélas & Masson 2008 (F1) §5** triangle variant (`agelas2008.txt:305-315`): continuity point
  = barycenter with weights 2/3 at vertex `s`, 1/3 at the other edge endpoint;
  `m^s_σ = m_σ/card(V_σ) = L/2` is a face-measure split, not a segment length. For this placement
  `B^s_K = I` ⇒ the MPFA-O scheme is **symmetric and unconditionally coercive on triangles** —
  a strictly better choice than the canonical midpoint the prior audit demanded. Only residual:
  `fvm_mpfao_theory.md:91-96` says "midpoint" and must be rewritten with the citation.
- **D6 (new).** `Factory::key()` (`cl_FVM_Factory.hpp:137-143`) computes `tA*(N+tB)` — not
  injective (N=6: pairs (4,0) and (3,2) both → 24) ⇒ silent facet overwrite in `mFacets` and
  wrong-facet retrieval in `create_subcells`. Fix: `tA*N + tB`, as `cl_FVM_MeshExtractor.cpp:336`
  already does. Likely a misplaced parenthesis.
- **D7 (new).** Distributed branch of `compute_facet_normals` indexes the node-pair-keyed
  `mFacets` map with **facet IDs** (`cl_FVM_Factory.cpp:291,358`; rank-0 path `:259` is correct).
  `Map::operator()` aborts on missing keys (`cl_Map.hpp:224`) ⇒ runtime abort once `commSize > 1`.
- **Row-sum blindness.** The historical "row sums of C = 1, conservation verified" claim can never
  detect the D1 sign bug: every `value` added to `A(k,l)` is also added to `B(k,i)`
  (`:941-942,956-957`), so `A·1 = B·1` and row sums of `C = A⁻¹B` are 1 for *any* normal signs.
  The real gate is a linear-field patch test (plan step R6).
- **Re-confirmed unchanged:** D1 (parallel master/slave normals + `+`/`+` assembly imposes a flux
  jump, not continuity — independently re-derived), D3 (no CMake wiring; `Tensor` 3-arg ctor does
  not exist, `cl_Tensor.hpp:61-121`; malformed `for` at `:986`), D4 (missing `mC = mQ − mP` and
  `set_size` in the MPI branch), D5 (stub solve path).

## Changes Made / Proposed

- **New:** `todo/fvm_implementation_plan.md` — execution plan per `todo/plan_template.md`:
  gap table (19 rows), ordered steps R1–R16 (Phase 0 build/hygiene → Phase 1 local math + patch
  tests → Phase 2 global assembly/BC/transient → Phase 3 PENTA6 extrusion mirroring
  `cl_ThinShellFactory::create` (`cl_ThinShellFactory.cpp:108-280`) → Phase 4 wrapper + Maxwell
  coupling), open questions O1–O5. No source code modified (read-only session per protocol).
- **Updated:** `todo/fvm_module_next_steps.md` — now the defect ledger: D2 struck through as
  FALSE POSITIVE with the Agélas citation, D6/D7 appended, verdict reconciled, the old "Fix D2"
  step rescoped to doc-only, the continuity-point open question RESOLVED (keep 1/3–2/3).
- **Updated:** `todo/README.md` — new "FVM Thermal Thin-Shell" section registering both files.

## Open Questions

- O1 — `commSize > 2` interaction-region ghosting policy (deferred, not next week).
- O2 — Dirichlet mechanics: eliminate boundary `Tc` strongly vs weak flux-form injection (decide
  during R8).
- O4 — centroid/surface computation is a side effect of `compute_normal` (`:709-720`); proposed
  explicit `compute_element_geometry()` pass during R3.
- O5 — keep the fixed Agélas continuity placement vs configurable η (default keep; revisit only
  if the R6 maximum-principle test fails).

## Files Updated

- todo/fvm_implementation_plan.md (new)
- todo/fvm_module_next_steps.md
- todo/README.md
- devlog/dl20260705_fvm_math_reaudit_and_plan.md (this file)
- devlog/README.md
- tmp/ai_exchange/fvm_module_math_audit.md (ephemeral thread, distilled here, GC-eligible)
