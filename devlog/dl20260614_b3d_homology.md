# Devlog: B3d Homology / Cohomology Cut Layer Whitepaper

**Date:** 2026-06-14
**Purpose:** Read-only documentation pass for Tier B3d - the homology/cohomology cut layer
(`src/homology/`) that makes the h-phi mixed formulation work on multiply-connected and periodic
domains.
**Module:** src/homology

## Context

Research-grade content: the topological precomputation (H^1 cohomology cuts) that lets a single-valued
magnetic scalar potential be used. Scope: the cut/(co)homology machinery; mesh topology (B3.0), graph
(B4b), and the Maxwell h-phi IWG (B3d-EM) are seams only. Read-only; `./archive` and `./nonfree`
positively not accessed.

## Method

Three parallel read-only subagents (algorithm; thick/thin + conjugate edges; periodic + consumption +
tests + provenance) plus direct reads of the core files. **Independent cross-check by GROK, not Codex**
(deliberate one-tier substitution - Codex budget low, per brief). Grok's first run hit the 20-turn
budget ("max turns reached"); re-run with GROK_MAX_TURNS=45 succeeded. GROK_EFFORT was NOT set (memory
caveat). Grok agreed on both high-risk facts and is recorded in todo/ai_exchange.md.

## Key findings

- **What/how:** computes the first cohomology H^1 of the phi-region (Air+Buffer+Ferro blocks;
  cl_Topology.cpp:382-388) as integer 1-cochains over edges, via integer Smith Normal Form on the
  coboundary matrices (fn_Smith, cited to Kaczynski-Mischaikow-Mrozek "Computational Homology"). A
  relative Homology from conductor terminals supplies "suggested homology" to orient/re-basis generators
  so each cut links a prescribed terminal current loop.
- **Four CutAlgorithms** (Pellikka, CCR, PellikkaGeneralized, BeltedTree) all implemented; the first
  three use SNF after a KMM/Pellikka coreduction, BeltedTree is a separate spanning-tree/cotree path that
  bypasses SNF.
- **Thick vs thin cuts:** both implemented; thick (cohomology +/-1 edge cochain) is converted to thin
  (pushed facet surface via node duplication + static condensation - Poincare-Lefschetz, NOT XFEM).
- **Conjugate edge:** the word appears ONLY in docs, never in source; realized in code as the signed
  "cut case" -> opposite/zero-coefficient facet face(tCase-1)/edge(tCase-1) (determine_cut_case_2d/3d).
- **phi-jump:** abstract node per cut carries the cut current I; CutSet duplicates conjugate-facet nodes
  with sources {abstract, original} weight +1 (phi' = phi + I).
- **Periodic:** slave-collapse in SimplicialComplex (identify opposite faces -> torus); periodic
  boundaries act as phi-boundaries; duplicate-node periodicity repaired; periodic-on-cut path is WIP.
- **Consumption (seam):** cuts() -> Cell<SideSet*> tagged DomainType::Cut; abstract_nodes() flagged for
  DOFs; Cut -> DOF table in Maxwell_FieldList.
- **Order-independent:** cut-finding is connectivity-level; FE order only picks the cut-facet element
  type. Not constrained by the order-2 edge cap.
- **Tests: NONE** - no topological-correctness test (k holes -> k generators; torus/genus), no SNF unit
  test. Validated only via end-to-end Maxwell apps. Biggest assurance gap.
- **Provenance:** SNF -> Kaczynski-Mischaikow-Mrozek [verified]; Pellikka -> Pellikka et al. 2013 SIAM
  DOI [verified, docs]; PellikkaGeneralized -> unpublished 2025 Giard preprint [partial]; BeltedTree ->
  NO citation [unknown]. Kotiuga/Gross/Dlotko/Specogna/Bossavit/Whitney cited nowhere.
- **Fragile/WIP/dead:** clean() greedy, no termination guarantee, mutates-during-iteration; |coeff|>1
  unsupported (determine_cut_case_3d hard-errors "THIS IS A BUG!"); clean_bfs() dead declaration; stray
  collect_thin_cut_edges splice (cl_CutData.cpp:441-445); SolveInt cout-failure path.

## Output

- Added `tmp/whitepaper/B3d_homology.md` (in-scope files first; confidence tags; step-by-step algorithm;
  thick/thin/conjugate section; provenance table; tests; IS-vs-WILL-BE fragile-path summary; Grok
  cross-check; open questions; assumptions). Pure ASCII.

## Verification

- Confirmed `./archive` and `./nonfree` not accessed.
- Confirmed whitepaper is ASCII-only.
- Cross-check by Grok (not Codex) recorded explicitly in the doc and in todo/ai_exchange.md.
- No source edited, nothing compiled, no tests run (read-only documentation task).
