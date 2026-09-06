# Side-Connector Bridge — Rev. 4 Q&A (conditioning, half partition, discretization)

**Date:** 2026-07-09
**Purpose:** Session log for rev. 4 of `todo/side_connector_effective_resistivity.md` — answers to
Christian's three questions, Codex-audited same day.
**Module:** fem/maxwell (analysis only — no source modified)

## What was asked

1. Stronger reasoning for why adding extra edge DOFs causes ill-conditioning (§2 table entry #2
   and the p-enrichment verdict were assertions, not arguments).
2. Current-partition physics of the two conducting halves around the φ-only buffer:
   (a) distribution if both halves stay isolated; (b) whether they act as one common conductor
   once the side-connector problem is solved.
3. How `a_R = R·I(h)·I(v)` is actually discretized — the stiffness-matrix build in BELFEM notation.

## What was added (all in `todo/side_connector_effective_resistivity.md`)

- **§2 "Why extra edge DOFs ill-condition the system — the strengthened argument."** Dichotomy:
  an extra edge DOF either cannot represent the 20 µm-governed jump (plain p-enrichment, Gibbs-
  limited) or can (Heaviside/sub-cell) and is then subject to five stacked mechanisms:
  (1) O(η)/O(1/η) two-sided spectral spread, order η⁻² ≈ 400 per perimeter element at
  η = t_w/L = 1/20, with cut-FEM κ ~ η^(−(2p+1))-type divergence flagged external [medium];
  (2) moving material contrast (quench-front dρ swings within one Newton step);
  (3) the digit budget — Bathe §8.2.6 Eq. 8.62 (`bathe.txt:27603–27611`) vs ε < 10⁻¹¹
  (Messe 2023, paper1) leaves ~5 decades of κ for the whole system (heuristic);
  (4) the bad subspace is global along the edge and unpivotable, and the only known controlling
  operator (ghost penalty) has no healthy same-phase neighbor to draw from;
  (5) ad-hoc unpaired H(curl) enrichment breaks the discrete de Rham complex
  (Monk `monk.txt:2416–2429`; Arnold `arnold.txt:222–230, 4288–4302, 667–688`) and collides with
  the φ-condensed boundary trace. Capstone: the wall's entire information content is one number
  per unit edge length — couple, don't enrich.
- **§5.3.2 "Discretization: building K from `a_R = R·I(h)·I(v)`."** Lumped form K += R·s·sᵀ
  (rank-one; the transport cut's symmetric sibling) and the distributed wall form
  `K_e = Σ_k w_k·detJ·r′·B(k)ᵀB(k)`, `B = [Em, −Es]` with t̂-projected traces — Kmm/Kms/Ksm/Kss
  sign pattern identical to `contact_impedance_theory.md` §7 with dS → dy. Into `K()` not `M()`;
  no h-Newton term (dr′/dT is a thermal cross-Jacobian if O4 goes monolithic); Joule heat
  `P′ = r′(Bq)²` at the same Gauss points; O1(a)/(b) changes only the DOF map (TᵀBᵀBT via the
  T-matrix path), not the algebra — O1 stays blocking.
- **§5.6 "The two halves: one conductor or two?"** Framing correction: one *cut* ≠ one
  *conductor* — the rerouted cut imposes I around the HTS loop only; the substrate circulation is
  free (Schnaubelt IC1/IC2 in one cut surface). Isolated (r′ → ∞): free coefficient =
  zero-applied-EMF loop equation ⇒ equivalent to an ideal terminal short outside the domain
  (right for soldered ends, but a modeling choice); steady partition ~10⁻⁶ by branch resistance,
  shorted-secondary transients during ramps, and **no local quench escape** — transfer happens
  over the whole tape length or not at all. Floating half ⇒ impose I_sub = 0 (input option, added
  to O2). With the bridge: physically one conductor, discretely two generators + weak-form
  impedance, deliberately — "one common conductor" is the r′ → 0 *limit behavior*, and the
  partition becomes the λ = sqrt(r′/(R′₁+R′₂)) transmission-line solution.

## Audit

Codex round 4 (same day, thread: `tmp/ai_exchange/side_connector_effective_resistivity.md`):
all three sections confirmed. Fixes applied: four→five mechanisms; "at least"→"order η⁻²";
digit budget marked heuristic; mechanism 5 re-scoped to ad-hoc *unpaired* enrichment (compatible
p-enrichment fails by the other branch); PSD phrasing ("adds no indefinite mode" — a huge r′ can
still stiffen the spectrum); dr′/dT cross-Jacobian caveat; terminal-closure claim reframed as a
boundary-condition assumption. §5.3.2 verified against `contact_impedance_theory.md:163–186` and
the T-matrix collapse path (`cl_MaxwellFactory.cpp:1319–1324,1463–1477,1524–1537`;
`cl_FEM_DofMgr_DofData.cpp:3555–3637`); §5.6 verified against `buffer_cut_topology.md` §3.2–3.5
and `schnaubelt2023.txt:337–351`.

## Status

No source modified. O1 (trace bookkeeping) remains the blocking design-note item; O2 gains the
floating-vs-terminated closure as an explicit input decision.
