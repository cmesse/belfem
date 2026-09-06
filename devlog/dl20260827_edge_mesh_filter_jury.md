# Devlog 2026-08-27 — Jury: tape-edge mesh refinement / checkerboarding hypothesis

**Date:** 2026-08-27
**Topic:** Three-AI literature jury on the hypothesis "too-fine tape-edge meshes are counterproductive; high gradients cause checkerboard-like oscillations; coarser meshes act as a natural filter for what machine precision struggles with"
**AIs involved:** Claude (pre-registered), Codex + Grok (blind jury)
**Claude Confidence:** high
**Codex Audit Confidence:** high
**Grok Audit Confidence:** high (~90%)
**Literature References:** Messe et al. 2023 §4, §5.1–5.2; Dular et al. 2021 §III, §V; Alves et al. 2022a/2022b; Alves et al. 2024 §III-A; Arsenault et al. 2021 §I/§III, 2023; Schnaubelt et al. 2023 §V; Bathe 2016 §7.4.3, §8.2; Hughes 2000 §4.3, Ch. 9; Belytschko et al. 2014 §8.3; Boffi et al. 2013 §8.10; Monk 2003 §5.3; Brenner & Scott Thm. 14.3.1; Zienkiewicz & Taylor Vol. 1 §7.7, Vol. 2 §3.2.7
**Verification:** reviewed only — literature-consistency level throughout; no executable gate ran. Exchange: `tmp/ai_exchange/review_edge_mesh_filter_hypothesis.md` (subject brief `_subject_edge_mesh_filter_hypothesis.md`).

## Summary

Hypothesis REFUTED as stated, 3/3 agreement under a pre-registered decision rule, every
auditor citation re-verified against the source extracts. One kernel of truth (round-off
floor vs. conditioning) and one legitimate neighboring context (structural dynamics)
were carved out explicitly. Physics adjudication remains Christian's.

## Key Findings

- **Edge refinement is prescribed, not warned against.** Messe et al. 2023 §5.1–5.2 use
  tip-refined meshes (5 µm / 20 µm at tips) as the smooth references; coarseness produces
  the staircase loss curve. Alves et al. 2022b calls high-aspect elements near tape
  extremities "necessary"; Arsenault 2021/2023 keep the superconductor mesh 5×/3× finer.
- **Checkerboarding (j ↔ 0/±2jc) is a convergence-tolerance artifact** (Messe et al. 2023
  §4): appears in BELFEM, COMSOL, and GetDP when ε is loose; cure is ε ≈ 1e-11, mesh size
  is not in the causal paragraph.
- **Tape-crossing oscillations are inf-sup/function-space instabilities** (Dular et al.
  2021): present with linear materials, worsen under refinement *only for unstable equal-
  order pairs* (β_δ ~ δ), cured by hierarchical enrichment — never by coarsening. The h-φ
  thin-shell coupling is reported immune (Alves et al. 2022a).
- **The canonical sharp-gradient oscillation runs opposite to the claim:** Galerkin
  oscillates when the mesh is too *coarse* (element Peclet > 2, Bathe §7.4.3).
- **Machine-precision kernel of truth:** digit loss s ≥ t − log₁₀ cond(K) (Bathe Eq. 8.62,
  aggravated by unequal element sizes and stiffness contrast); Zienkiewicz Vol. 2 §3.2.7
  puts the full-Newton tolerance guideline at half machine precision (~1e-8), so BELFEM's
  1e-11 target deliberately sits near the round-off floor. But no source ties the
  attainable floor to tape-edge h; power-law degeneracy lives at small |j| (tape core),
  and formulation choices (dummy ρ_air) dominate conditioning (Arsenault 2021: 2.6e15 vs
  3.1e7). Remedies are solver/space/tolerance-side, never deliberate coarsening.
- **Where the "coarse mesh filters high frequencies" recollection IS true:** structural
  dynamics / wave propagation — the upper discrete spectrum is spurious and Newmark/HHT
  integrators are designed to filter it (Hughes Ch. 9; Bathe mode-superposition cutoff).
  Transfers weakly to parabolic eddy-current problems, which damp high modes physically
  and use L-stable BDF integration.

## Changes Made / Proposed

- No source changes (read-only round). Subject brief + full jury record in
  `tmp/ai_exchange/` (ephemeral; this entry is the distillation).

## Open Questions

- Grok proposed the right executable gate if this ever needs settling beyond literature:
  measure attainable ε vs. tape-edge h at fixed n, formulation, mesh quality, and solver.
  Not scheduled; protocol stop-condition says gate-before-another-round.
- Whether a *specific* BELFEM deck's edge resolution is past diminishing returns for its
  loss functional — Christian's call, per the physics-adjudication rule.

## Files Updated

- devlog/dl20260827_edge_mesh_filter_jury.md (this file)
- devlog/README.md (index line)
