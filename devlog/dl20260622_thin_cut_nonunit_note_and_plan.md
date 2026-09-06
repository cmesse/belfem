# Devlog 2026-06-22 — Non-Unit Thin-Cut: Peer Note, Literature, Audit, Implementation Plan

**Date:** 2026-06-22
**Topic:** Turn the 2026-06-12 analysis of non-unit thick-cut coefficients (`|c(e)| ≥ 2`,
coarse CCT meshes) into a shareable peer-facing technical note, reconcile it with newly-added
literature, audit it with Codex + Grok, and produce an implementation plan.
**AIs:** Claude (primary), Codex (audit + prose), Grok (audit, beta).
**Docs first / read-only except the docs and plan listed below; no source code modified.**

## Summary

Christian asked for the 2026-06-12 findings to be tidied for his peers (ultimately Prof.
Sirois). Output is a standalone theory note plus a phased implementation plan. The note's
mathematics survived a two-AI audit with **one real correction** (a feasibility-vs-objective
conflation in the optimal-representative section, caught by Codex). A key design conclusion for
implementation: the rectification is a **graph** algorithm on the mesh 1-skeleton, and its
state belongs in **external index-keyed arrays**, not on `Vertex`.

## What was done

1. **New peer-facing note** `src/homology/doc/thin_cut_nonunit_rectification.md` — lifted the
   substance out of `todo/thin_cut_nonunit_coefficient_algorithm.md` (working notes) into a
   self-contained note: problem → geometric meaning of `|c|≥2` → manifold-vs-fixed-mesh
   obstruction → difference-constraints/Bellman–Ford → **two regimes** → source-field
   foundation → normal-surface identification → proposed algorithm → references. Stripped all
   AI scaffolding (confidence tags, QA/QB codes, provenance). Indexed in the module README.

2. **Literature reconciliation.** Christian added a topology cluster (`papers/topology/`:
   Dey 2011, Chen & Freedman, Dunfield & Hirani, Costantini), the Gross & Kotiuga book, CLRS,
   Haken 1961, and the Dular 1997/1999 h-formulation papers. Verified all citations/DOIs in the
   note against the library; corrected the note's Chen venue/ℤ₂ scope; **fixed a wrong DOI in
   the library** (`dular1999` `10.1109/20.767289` → `10.1109/20.767308`, 8 occurrences across
   `papers/fem/index.md` + `dular1999.md`).

3. **Dular ↔ note cross-check + XFEM question.** Confirmed the Dular papers are the continuous
   foundation: `c − dθ` ↔ curl-kernel freedom `h_s + grad χ` (Dular 1997 §III.A); minimal
   support ↔ optimal representative (§II.B); the cut jump = transport current with binary `q_i`
   (Dular 1999 Eq. 9). Added a "source-field foundation" section. Settled the XFEM question:
   Dular cuts are **not** XFEM (conforming, cohomology lineage — Bossavit/Kotiuga, cited in the
   papers); XFEM/phantom-node (Hansbo & Hansbo 2004) is the in-element alternative BELFEM
   declines by design — kept as a one-line note under option C, not a section.

4. **Two-AI audit** (`tmp/ai_exchange/thin_cut_nonunit_audit.md`). All code citations
   re-verified against source before editing. Applied:
   - **Codex MUST-FIX (Grok missed):** §3 conflated feasibility (L∞ `‖c−dθ‖∞ ≤ 1`) with the
     weighted-L1 objective; deleted the false "feasibility is the `w ≡ 1` special case" and
     separated the optional optimization from the feasibility test.
   - **Both + user:** defined the `{0,13,19,38,56,30,43,53}` 6-bit edge bitmask
     (`pattern = Σ 2ᵏ[edge k]`, `cl_CutData.cpp:447-475`) with a full decode table and a
     one-line enumeration justification.
   - **Both:** corrected the Summary's non-unit failure mechanism (magnitude-blind
     `determine_cut_case_3d`; non-unit edges dropped to `weight()=0` by `collect_coefficients`,
     then rejected as invalid pattern/coefficient).
   - Precision: normal-surface "parallel copies" → "several disk levels"; defined option (A);
     softened refinement termination + "localizes the throat"; Chen "2010/2011"; Gross Ch. 6.
   - Added a **Notation** table after the Summary (resolves θ-vs-Φ — kept `θ` as the integer
     node potential, distinct from physical `φ` — and defines `d` as the coboundary, killing
     the "infinitesimal" misreading).
   - Grok ran clean (all citations accurate); Codex caught the load-bearing math — consistent
     with the Codex-primary / Grok-tiebreak division.

5. **Implementation design + plan** `todo/thin_cut_nonunit_rectification_implementation.md`
   (Codex prose-polished). Settled with Christian:
   - The solve is a **graph** algorithm on the **1-skeleton only** (the complex/faces are used
     upstream to produce `c` and guarantee closedness `dc=0`; the solver consumes that
     invariant). Run it on the mesh `Node`/`Edge` adjacency + `Cochain` directly, **not** a
     rebuilt `graph::Vertex` (unweighted/undirected, no weighted SSSP; archived Tarjan is
     reference-only).
   - **Solver state in external arrays keyed by vertex index** (`θ` signed, `parent`, `inQueue`,
     `count`) — `Vertex` untouched. Rejected putting state on `Vertex` (option 1: `mLevel` is
     unsigned, can't hold negative `θ`, and is insufficient; option 2: `mData` swap forces
     per-vertex allocation + indirection on the hot loop and a base-class refactor with BFS/DFS/
     RCM blast radius). Arrays are zero-permanent-footprint, cache-friendly, MPI-safe.
   - Periodic **quotient** via slave→master `periodic()` fold (mirrors `clean()`).
   - Phased: (0) read-only feasibility diagnostic answering **QB** first, (1) SPFA + cycle
     extraction, (2) rectify-or-certify replacing `clean()`'s body + global `dc=0`/`|c|≤1`
     asserts, (3) tests, (4) optional min-cost-flow / certificate-guided refinement.

## Files changed

- `src/homology/doc/thin_cut_nonunit_rectification.md` (new), `src/homology/doc/README.md` (index)
- `todo/thin_cut_nonunit_rectification_implementation.md` (new), `todo/README.md` (index)
- `literature/papers/fem/dular1999.md`, `literature/papers/fem/index.md` (DOI fix — separate repo)
- `tmp/ai_exchange/thin_cut_nonunit_audit.md` (audit thread; `.out`/`.log` scratch swept)
- this devlog + `devlog/README.md`

## Open / next

- **QB diagnostic** (Phase 0) is the recommended next action: cheap, conclusive, changes nothing.
- Open questions carried in the plan: periodic-quotient correctness (QC), domain scope of the
  index map, higher-order (non-corner) edges, large-mesh cost.
- The note is shareable with Prof. Sirois as-is.

## Addendum — plan design iteration + two-AI audit (later 2026-06-22)

After the plan was drafted, Christian iterated the data-structure design and the plan was
audited by Codex + Grok before any coding.

**Design refinements (with Christian):**
- The solve is a **graph** algorithm on the **1-skeleton only**; the complex/faces are upstream
  (they produce `c` and guarantee closedness `dc=0` — the solver consumes that invariant).
- **No compaction map and no quotient map:** key the scratch directly by `node->index()` (the
  vertex index *is* the offset); fold periodic slaves to masters *inline* via the typed
  `periodic()` pointer. Rejected `reinterpret_cast` for the fold (skips base-subobject offset
  adjustment, UB-prone in the `Node`/`Edge`/`Vertex` hierarchy — implicit upcast / `static_cast`
  only, and the fold needs none).

**Audit (`tmp/ai_exchange/thin_cut_impl_plan_audit.md`); all code-grounded claims re-verified
against `clean()`/`Cochain` before adoption.** Core design held; gaps were in the BELFEM
realization details:
- **Grok** (strong round, F1/F3/F4/F7 verified real): loop over `mGenerators(1)` with per-gen
  reset; relax only the `original_edges()` domain; quotient adjacency needs master∪slave incident
  edges (not endpoint-fold alone — clean() `:269-301`); `rep()` via flag/`is_flagged()` (no
  `is_slave()`).
- **Codex** (corroborated + new BELFEM-fit): `Vector` is LA-only and `containers/Queue` wraps
  `std::queue` → preallocated `Cell`/ring-buffer; `c−dθ` must keep Cochain boundary/coboundary
  side-state consistent; guards must be `BELFEM_ERROR` not `BELFEM_ASSERT`; int64 distances;
  null-`mSimplicialComplex` fallback; and **corrected Grok** — a folded slave→own-master edge
  with `|c|≥2` is a **length-1 infeasibility certificate**, not "always satisfied" (= the theory's
  `|⟨c,z⟩| > length(z)`).
- Grok miss: claimed the archived Tarjan "not in repo" (globbed `src/` only; it's in `archive/`).

All adopted into `todo/thin_cut_nonunit_rectification_implementation.md` (State model, Phases 0–2,
Risks + audit provenance). `todo/README.md` marked "audited, pending approval." Division of labor
held: Grok found real omissions (citations verified), Codex caught the load-bearing BELFEM-fit
issues. Scratch `.out` swept; exchange thread retained. No source modified.

## Addendum (2026-06-22, later) — retired the superseded analysis file

Deleted `todo/thin_cut_nonunit_coefficient_algorithm.md` (the 2026-06-12 working notes). It was
fully superseded — verified before deleting that every substantive item lives in a durable home:
- geometry / existence obstruction / Bellman–Ford / cut cases / **QA seven-case enumeration**
  ({0,13,19,38,56,30,43,53}) / verified-DOI literature → `src/homology/doc/thin_cut_nonunit_rectification.md`;
- algorithm / **QB** (now Phase 0 diagnostic) / **QC** (periodic quotient design) / missing
  post-clean `|c|≤1` guard → `todo/thin_cut_nonunit_rectification_implementation.md`;
- the 3-AI analysis record + Grok's code findings (**`clean()` iterator-mutation fragility**;
  **silent ±2 weight-0 drop in `CutData::collect_coefficients()`**; no global post-clean assert)
  → `devlog/dl20260612_thin_cut_nonunit_algorithm.md:83-101`.

Nothing unique was lost (pure dedup — unlike the same-day `hdf5_writer_repair.md` merge, which
carried unique bug content). Updated `todo/README.md`: removed the analysis-file row, folded a
pointer to the theory note + 0612 devlog as the "analysis/theory home" into the implementation-plan
entry. The two devlog references to the old filename are left intact (immutable historical records).
No source modified.
