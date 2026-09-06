# Devlog 2026-08-22 — Coulomb Gauge Penalty, Step A (theory note + blind audit)

**Date:** 2026-08-22 (overnight session)
**Topic:** Independent derivation and audit of the proposed Coulomb-gauge penalty
for the h-φ formulation; Q1–Q8 of the task brief answered; G-operator
prerequisites established. Steps B/C not started (stop point).
**AIs involved:** Claude (derivation), Codex + Grok (blind, independent audit)
**Claude Confidence:** high on the math results, medium on conditioning attribution
**Codex Audit Confidence:** high (probes re-executed, output matched)
**Grok Audit Confidence:** high (~85%) (probes re-derived, no shell)
**Literature References:** Monk 2003 §5.5.1, §7.2.1, §7.4, §3.8; Boffi et al. 2013
§11.4; Messe et al. 2023 Eqs. (6)–(8); Dular et al. 1997 §III; Arsenault et al.
2021; Denis et al. 2026
**Verification:** numeric probes P1–P5 executed this session (math-level, output
quoted in the note §9); assembly contract source-traced and independently
re-traced by both auditors. No compiled-class gate ran — that is Step C.4's
deliverable. No source code was modified.

## Summary

Wrote `todo/coulomb_gauge_penalty_theory.md`: independent first-principles
derivation of a Coulomb gauge penalty for the h-formulation, then a term-by-term
audit of Christian's draft formula, then a blind two-auditor round, then
revision. Central results: (1) the draft formula is dimensionally inconsistent
with a dimensionless χ; a consistent K-channel form γ = χ·ρ*/μ² exists and needs
no new plumbing; (2) on TET4/TRI3 and the thin-shell kernels, every
element-local penalty built from ∇h is identically zero on the discrete curl
null space — including cohomology/cut fields, which remain element-wise
constant in Whitney-1 — so the penalty cannot address the conditioning problem
on those meshes at all; (3) on TET4 the full-gradient Gram equals exactly half
the curl-curl Gram, so a G′G penalty there is a pure artificial-resistivity
perturbation; (4) the penalty has genuine (but partly mesh-artifact) content
only on hex and higher-order families. The G-operator itself remains justified
(diagnostics, postprocessing, hex-mesh variants, Step C).

## Key Findings

- Assembly contract for any penalty: K-channel terms are auto-Δt-scaled and
  residual-consistent; frozen-coefficient penalties need no dJdx and no manual
  RHS (`cl_IWG_Timestep.cpp` bdf2, `cl_TimestepMatrices.cpp` assemble_dJdx,
  `cl_FEM_DofMgr_SolverData.cpp` NewtonRaphson branch: r = A·x − b,
  solve (A+dJdx)Δ = r, x ← x − ωΔ).
- Monk §7.4 on a strong div penalty for edge elements: "the straightforward
  answer is 'no'". Edge-element-compatible alternatives: weak (lumped-mass)
  divergence (Monk Eq. 7.43–7.44), tree–cotree gauging (Dular 1997 §III,
  Denis 2026), jump penalties.
- Audit refutations of my draft (both corrected): `mRhoMin` is 0.0 (floor
  removed by ruling 2026-08-10, `cl_Material.hpp`), not 1e-16 — I had recalled
  instead of re-read (L-08); "HEX8TS = legacy wrap" is wrong (only the wrap
  construction was removed; `EF_HEX8TS` remains thin-shell machinery).
- Grok: P5 probe was structurally unable to fail (L-01) — reclassified as a
  construction check; the real falsifiers are Step C.4's compiled-class tests.
- Grok: full G′G is non-vacuous even on box HEX8 (∇(xy) kernel mode) — the
  no-go is a simplex/Whitney property, not universal.
- Seven EF subclasses declare a private `mG` that collides with Step B's
  planned base-class member name — flagged for sign-off before Step B.
- Drive-by (Grok, no action taken): stale comment "A = 3*M + 2*h*K" in
  `cl_IWG_Timestep.cpp` bdf1.

## Changes Made / Proposed

- NEW `todo/coulomb_gauge_penalty_theory.md` (theory note, audited + revised).
- Exchange: `tmp/ai_exchange/coulomb_gauge_stepA.md` (pre-registration +
  reconciliation), `coulomb_gauge_stepA_codex.md`, `coulomb_gauge_stepA_grok.md`.
- No source edits (read-only session, per task rules).

## Open Questions

Listed in the note §14 for Christian: provenance of the 1e17 condition number
(L-04 — measure before remedy); target regime (normal-zone vs global); the ρ*
averaging fork; mesh family of the quench decks; §12 Step-B decisions (LINE3
G contract, HEX8TS stub, mG rename).

## Stop Point

Per the task brief, work STOPS here awaiting Christian's sign-off on Step A.
Steps B (EdgeFunction base class) and C (per-element G) are planned but not
started; task-list items exist and are blocked.

## Files Updated

- todo/coulomb_gauge_penalty_theory.md (new)
- devlog/dl20260822_coulomb_gauge_stepA.md (this file)
- devlog/README.md (index line)
- todo/README.md (registration)
