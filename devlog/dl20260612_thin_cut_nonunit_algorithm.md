# Devlog 2026-06-12 — Thin-Cut Generation for Non-Unit Coefficients: 3-AI Investigation

**Date:** 2026-06-12
**Topic:** Algorithm for generating thin cuts when the cohomology generator carries
`|coeff| > 1` (conductor looping under itself, e.g. CCT magnets, on coarse meshes).
Independent investigations by Claude, Codex, and Grok; comparison and synthesis.
**AIs involved:** Claude (primary), Codex (independent + comparison), Grok (independent + comparison; first substantive investigation run)
**Claude Confidence:** high (core algorithm and obstruction), medium (literature attributions pending DOI check)
**Codex Audit Confidence:** high
**Grok Confidence:** high on code-grounded claims (citations verified by Claude per protocol)
**Literature References:** internal — `src/homology/doc/thick_thin_cuts_and_conjugate_edges.md`; external (NOT in ./literature, unverified) — Costantini 1998 (IEEE TGRS, MCF phase unwrapping); Dey–Hirani–Krishnamoorthy 2011 (SIAM J. Comput., TU/LP optimal homologous cycles); Chen–Freedman 2010/11 (SODA, localization hardness); Dunfield–Hirani 2011 (SoCG, least-area Seifert surfaces); Kotiuga 1987/1989 + Gross & Kotiuga 2004 (harmonic level-set cuts); Haken 1961 (normal surfaces)

## Summary

The question (Christian): coarse CCT-type meshes produce thick-cut generators with
edge coefficients `|c(e)| ≥ 2`, which the thin-cut pipeline cannot consume. Since
Poincaré–Lefschetz suggests a thin cut "should always exist," what algorithm generates
it? Christian's hunch: a documented-but-obscure graph algorithm.

The hunch is confirmed. All three AIs — Codex and Grok working **without access to
Claude's analysis** (independence enforced in the prompts, then compared) —
independently derived the same answer:

1. `|c(e)| ≥ 2` is **sheet multiplicity**: the Poincaré–Lefschetz dual cut surface
   passes more than once through the same dual 2-cell (in a CCT: the spanning surface
   squeezed repeatedly through a one-element-thick inter-turn throat).
2. Poincaré–Lefschetz guarantees the class (and an embedded representative in the
   *manifold*), **not** a unit-coefficient representative in a *fixed* mesh. The exact
   fixed-mesh obstruction is a closed edge-loop `z` with `|⟨c, z⟩| > length(z)` (e.g.
   a bore loop of an N-turn winding meshed with fewer than N edges).
3. The rectification `c′ = c − dθ`, `θ` integer on nodes, with `|c′| ≤ 1`, is a
   **system of difference constraints** — the obscure-but-classical graph problem:
   `θ(j) − θ(i) ≤ 1 + c(e)`, `θ(i) − θ(j) ≤ 1 − c(e)`. **Bellman–Ford** decides it:
   feasible → integer `θ` from shortest-path distances (total unimodularity of the
   node–edge incidence matrix gives integrality for free); infeasible → a **negative
   cycle**, which is precisely a loop with `|⟨c, z⟩| > length(z)` — a constructive
   certificate that *no* cohomologous unit representative exists, localizing the
   too-coarse region.
4. The current greedy `Cohomology::clean()` (`cl_Cohomology.cpp:195`) is a local
   node-firing relaxation with **no feasibility test, no termination guarantee, and no
   objective**; it must be replaced or augmented by the global solve.

## Key Findings

### Mathematics (3-way independent agreement)

- **Obstruction proof:** adding `dθ` telescopes to zero around any closed loop, and a
  unit cochain contributes at most ±1 per traversed edge, so `|⟨c, z⟩| ≤ length(z)` is
  necessary; Bellman–Ford feasibility shows it is also sufficient. The negative-cycle
  certificate is exactly the violating loop (sum of constraint weights around the cycle
  is `length(z) ∓ ⟨c, z⟩`). [3/3 agree, high]
- **Optimal representative:** minimizing `Σ_e w_e |c(e) − (dθ)(e)|` over integer `θ` is
  an exact **min-cost circulation** problem: the LP dual maximizes `⟨c, y⟩` over
  circulations `Ay = 0` with capacities `|y_e| ≤ w_e`; TU + integral data ⇒ integral
  optima on both sides. Feasibility is the `w ≡ 1` special case. (Codex initially held
  this at medium pending the explicit primal/dual; discharged by Claude during
  comparison.) [high]
- **QA proven (was medium ~75%):** unit coefficients + cocycle closedness (zero
  oriented circulation on every tet face) force the per-tet support to be exactly one
  of {0, 13, 19, 38, 56, 30, 43, 53}, each nonzero pattern with exactly the two ± sign
  assignments the code admits. Verified three ways: Claude brute-force enumeration of
  all 3⁶ assignments under the 3 independent face constraints; Codex's independent
  enumeration (matching result, TET4 corner-edge convention); Grok's case analysis.
  **Consequence: after rectification, the existing thin-cut pipeline (cut cases,
  duplication, relinking, condensation) applies unchanged.** Caveat: recheck if
  higher-order elements ever contribute non-corner edges to the cochain. [high]
- **Normal-surface identification:** the seven cut cases of `determine_cut_case_3d()`
  (`cl_CutData.cpp:425-639`) are exactly the normal-disk inventory of a tetrahedron
  (4 vertex-cap triangles ±1…±4 + 3 vertex-pair-separating quads ±5…±7, patterns
  13/19/38/56 and 30/43/53 verified against the vertex stars / separators).
  `|coeff| ≥ 2` = normal coordinates with multiplicities (parallel disk copies) —
  standard in Haken theory, merely unrepresentable by facet-pushing. Grok's framing
  adopted: an *identification of the disk inventory*, useful intuition, not a full
  embedding into normal-surface machinery. [high]
- **Refinement termination (Codex softening adopted):** *certificate-guided* splitting
  restores feasibility in finitely many steps (`⟨c, z⟩` is invariant; each guided split
  lengthens the violating loop); arbitrary refinement elsewhere proves nothing; a
  global-subdivision fallback guarantees termination since the continuous embedded
  representative exists. [high]

### Code findings (Grok net-new, citations verified by Claude per protocol)

- **`clean()` iterator-mutation hazard — verified, refined:** the range-for at
  `cl_Cohomology.cpp:229` iterates `getSimplicesMap()` while `addSimplexToCochain()`
  (`cl_Cochain.hpp:251-275`) erases/inserts in the same map, with no break. **Not UB
  today**: `OrderedMap` wraps `std::map` (`cl_OrderedMap.hpp:28`) — insertion never
  invalidates iterators, and the current element (an edge with `|coeff| ≥ 2`) cannot be
  zeroed by a single ±1 firing, so only non-current elements get erased (safe for
  `std::map`). Latent UB if the backing container ever changes (flat/unordered map);
  plus a benign pass-order wrinkle (entries inserted before the current key are only
  seen on the outer `while` rescan). Classification: **fragility, not active bug**.
- **`clean()` non-termination input (Codex):** an infeasible directed 3-cycle with all
  coefficients 2 returns to the same state after one full firing pass; the
  `while(tIllegalCoeff)` loop (`cl_Cohomology.cpp:224`) has no iteration cap — contrast
  the downstream 3D peel's `BELFEM_ERROR(k<1000, ...)` guard (verified). On feasible
  inputs greedy can still terminate with non-minimal support (no objective).
- **Silent weight-0 drop (Grok, verified):** `CutData::collect_coefficients()`
  (`cl_CutData.cpp:212-243`) tests `== ±1` only — a surviving ±2 edge sets *neither*
  bitset. It fails loudly only if a selected cut element exposes it through the
  "Invalid cut coefficients" path; an edge outside selected elements is silently
  dropped. There is **no global post-clean `|c| ≤ 1` assertion** before `CutData` runs.
- **Diagnosability (Grok):** the `BELFEM_ERROR` messages in `determine_cut_case_3d()`
  cannot distinguish "non-unit survived clean()" from "unit but inadmissible pattern"
  (the latter would now indicate a generator/closedness bug, per QA).
- **Periodic quotient (Grok, extends open QC):** `clean()` is already periodic-aware
  (slave flagging + partner firing); any replacement Bellman–Ford/MCF solver must run
  on the **quotient graph** (identify θ across periodic node pairs) or certificates and
  reduced cochains will be wrong on PBC meshes.
- **Corrected Codex over-claim:** the stray block-scope declaration at
  `cl_CutData.cpp:441-445` is *legal* (dead) C++, not a compile blocker — already
  tracked as Step 8f cleanup in `todo/periodic_thin_cut_continuity_fix.md`.

### Recommended algorithm (consensus)

```text
rectify_to_unit_or_certify( c ):
  1. Build the difference-constraint graph on the (quotient) mesh node set.
  2. Bellman–Ford (or min-cost circulation with weights w_e for an optimal-support cut).
  3. Feasible  -> apply c := c − dθ once; assert closedness and |c| <= 1 globally;
                  run the existing thin-cut pipeline UNCHANGED (justified by QA).
  4. Infeasible -> extract the negative cycle (parent pointers); fail loudly with the
                  certificate loop (edge IDs), recommending targeted refinement of the
                  throat — or feed an automatic certificate-guided local refinement
                  pass and re-run. Multi-level duplication (q_k = p + k·I, one
                  duplicate per branch level; covering-space element labels) remains
                  the heavy general fallback if refinement is ever unacceptable.
```

The infeasible case thereby turns from a crash (or silent hang in `clean()`) into an
actionable, physics-aligned adaptivity instruction — the throat is where field
gradients are large anyway.

## Changes Made / Proposed

No source edits (read-only session per protocol). Files written:

- `todo/thin_cut_nonunit_coefficient_algorithm.md` — problem analysis + proposal (new)
- `todo/README.md` — index entry (updated)
- `todo/ai_exchange.md` — CLAUDE query, CODEX response, GROK response, CLAUDE resolution
- this devlog + `devlog/README.md` index entry

Proposed implementation order (future session, user approval required for source edits):

1. **QB diagnostic first:** prototype the Bellman–Ford feasibility check as a pure
   diagnostic on a real coarse-CCT generator — decides whether observed failures are
   "feasible but greedy `clean()` failed" or "genuinely infeasible".
2. Implement `rectify_to_unit_or_certify()` (quotient-aware), replacing the greedy
   `clean()` body; keep `clean()`'s signature/call sites.
3. Add the global post-rectification `|c| ≤ 1` + closedness assertion before `CutData`.
4. Split the `determine_cut_case_3d()` error message into the two distinguishable
   failure modes.
5. Certificate-guided refinement integration (separate design, interacts with meshing).

## Open Questions

- **QB:** feasible-but-greedy-failed vs. genuinely infeasible on the real CCT
  reproducers (cheap to answer with the step-1 diagnostic).
- **QC:** quotient-graph formulation details under periodicity (and interaction with
  the Step 5/6 seam work in `todo/periodic_thin_cut_continuity_fix.md`).
- **QD:** verify the external literature attributions (Grok judged the concepts sound
  and "no outright misattributions visible," but none were checked against the actual
  papers; all three AIs flag them unverified).
- Higher-order elements: confirm the cochain only ever lives on corner edges (QA was
  proven for the TET4 six-edge convention).

## Files Updated

- todo/thin_cut_nonunit_coefficient_algorithm.md (new)
- todo/README.md
- todo/ai_exchange.md
- devlog/dl20260612_thin_cut_nonunit_algorithm.md (new, this file)
- devlog/README.md

## Process notes

- Independence protocol worked: both secondary AIs were instructed to write phase-1
  analysis before reading Claude's file or the trailing exchange entries; both reported
  honoring it, and both reproduced the core results from the code + internal doc alone.
- Grok (beta) practical notes: `GROK_EFFORT=high` currently fails hard — the backing
  model `grok-build` rejects the `reasoningEffort` parameter (400); ran at default
  effort. All of Grok's load-bearing file:line citations were verified accurate
  (one finding refined: iterator hazard is fragility, not UB). Citation checking is
  cheap, as predicted.
- `todo/ai_exchange.md` is ~4900 lines, far past the ~500-line archive threshold;
  left un-archived because active periodic threads and the Codex delta-position file
  point into it — archive deliberately in a quiet moment.
