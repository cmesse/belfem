# Devlog 2026-07-31 — HEX8TB Edge Function: Theory Consultation, Implementation, Campaign Doc

**Date:** 2026-07-31
**Topic:** Three-AI theory consultation (Q1–Q5) for the HEX8TB side-connector edge
function, implementation of `EF_HEX8TB` (metric, Jacobian, E/C operators), and the
comprehensive side-connector documentation page.
**AIs involved:** Claude (primary), Codex + Grok (independent theory answers + adversarial
audits of code and documentation, prose pass by Codex)
**Claude Confidence:** high (all formulas triple-derived independently and identical;
implementation Codex+Grok verified at high confidence)
**Literature References:** Monk 2003 (H(curl) conformity, first-family hex Nédélec);
Alves et al. 2022a (paper5); Schnaubelt et al. 2023 (paperA); Messe et al. 2023 (paper1)

## Summary

Christian's four theory questions (integration metric, Jacobian/Nabla, E/C operators,
h_b/h_n recovery) plus a fifth ([-1,1] vs [0,1] edge parametrization) went through the
established protocol: cleaned prompt → independent answers from Claude, Codex, and Grok →
merge. Q1–Q3 came back unanimous with identical formulas from all three; Q4 split 2:1 and
was resolved with two hard rules; Q5 confirmed Christian's intuition (edge dofs are
ampere-valued circulations — parametrization-independent). Q1–Q3 were implemented in
`EF_HEX8TB` the same day and audited clean; Q4 is documented as a design, deliberately not
implemented. The campaign reference doc landed in
`src/fem/maxwell/doc/side_connector_wall_element.md`.

## Theory results (full record: `tmp/ai_exchange/hex8tb_edge_function_theory.md`)

- **Q1 — metric (unanimous):** exact cuboid. Directions from geometry (tangent from face
  midpoints, frame orthonormalized), magnitudes exactly from the blocks. Decisive: the
  wall's outer-node offsets are a drawing proxy; a node-based metric would couple the wall
  conductance r′ to the visualization and carry Gram-matrix conditioning ~(L/w)².
- **Q2 — Jacobian (unanimous):** J rows = (L/2)tᵀ, (w/2)bᵀ, (d/2)nᵀ; Nabla columns =
  (2/L)t, (2/w)b, (2/d)n; detJ = Lwd/8. Constant per element, computed once in `link()`;
  no Gram matrix, no LU; `update_nabla()` degenerates to a no-op (kept because
  `Calculator::dV_hex` calls it).
- **Q3 — operators (unanimous):** E(:,k) = s_k F_k (2/L) t;
  C(:,k) = s_k[(4/(Ld)) ∂F_k/∂ζ b − (4/(Lw)) ∂F_k/∂η n]. **j_t ≡ 0 exactly** — correct by
  design (transport current lives in the shell). Unifying picture: **h_t is the stream
  function of the cross-flow** (j = ∇h_t × t); the current transferred between two
  boundary points equals their h_t difference — this is the r′ stencil, and it resolved
  the Grok-vs-Claude "j_b or j_n is the transfer current" split (same current, different
  control surfaces).
- **Q4 — recovery (2:1 + refinement, DOCUMENTED ONLY):** h_n from the two adjacent
  layers' shell recovery (compute_hn analog, per-gap average); h_b primarily from the
  shells' in-plane field at the edge station (Codex sided with Claude; Grok preferred
  air-∇φ — rejected: tape edge is a corner singularity, and reading a solved field adds
  no free dof). Hard rules kept from Grok's warning: never add binormal/normal edge dofs;
  never feed recovered components into the curl/stiffness operator. Blend donors with
  N_k = 2F_k. Recovery must route through the adjacent shells' side-curve facets (the
  wall's faces have no volume neighbors).
- **Q5 — parametrization (Christian's question, confirmed):** hex [-1,1] vs simplex [0,1]
  is immaterial: every BELFEM edge function is normalized to unit circulation
  (∫ E_k·dl = s_k δ_jk, dof unit = ampere), and at lowest order the tangential trace along
  a shared edge is constant in both families — conforming pointwise, not just in the
  integral (verified against `cl_EF_TRI3.cpp` Whitney functions).

## Changes made

- `src/fem/interpolation/nedelec/cl_EF_HEX8TB.{hpp,cpp}` — rewritten: exact-cuboid
  metric, orthonormal right-handed frame (b = n×t), constant mJ/mInvJ/mDetJ filled in
  `link()`, direct C assembly via constant curl factors mP = −4/(wL)·n, mQ = +4/(dL)·b,
  no-op `update_nabla`. Dropped dead members (Gram matrices, QUAD4 interpolator,
  Ex/Ey/Ez scratch). Fixed two latent WIP API errors: `get_coordinates` → `get_coords`,
  fem-Element `block_id()` → `element()->block_id()`. Syntax-verified against real build
  flags.
- `src/fem/interpolation/cl_EdgeFunctionFactory.cpp` — HEX8TB case (+ include).
- `src/fem/interpolation/nedelec/fn_num_nedelec_dofs.hpp` — HEX8TB → 4.
- `src/fem/maxwell/doc/side_connector_wall_element.md` (new) — campaign reference:
  physics, failure history, r′ theory, HEX8TB design, EF derivation (Q1–Q3, Q5),
  recovery design (Q4), status. Registered in `src/fem/maxwell/doc/README.md`.
- `src/fem/interpolation/doc/nedelec_thinshell.md` — historical note updated (HEX8TB is
  back as a different element; old HEX8TS wrap references remain stale).
- `todo/hex8tb_phase2_fem_wiring.md` — R1/R2/R3 ticked; R1 sketch frame-label corrected
  (∇η → ∇ξ); O1 partially defused (EF builds its own frame, detJ > 0 by construction);
  new D4 (physical_tag clobbered by `MaxwellFactory::set_physical_tags_for_elements` →
  thickness lookup fragile, R5 gate) and D5 (ζ-gap semantics undecided: layer j vs mean).

## Audit record

- **Implementation audit (Codex + Grok, same day):** all math checks PASS at high
  confidence on both sides — E/C signs exact against the dF tables, frame right-handed,
  unit circulation holds, cuboid metric fills consistent, `dV_hex` safety, no hot-path
  allocations, exactly 4 edge signs. No defects in the EF math. Open findings are wiring
  concerns (D4/D5 above) plus the known WIP stop in the factory
  (`BELFEM_ERROR "not implemented yet"` — deliberate, side connectors cannot reach
  assembly yet).
- **Documentation audit (Codex + Grok):** every EF/mesh/template/recovery claim verified
  EXACT against the code by both. Thirteen findings, all applied — highlights: the
  `nedelec_thinshell.md` §4 summary still declared HEX8TB/HEX8TS removed (contradiction
  with the revised top note, fixed); "duplicate-safe by construction" overclaimed while
  the raw-index station-lookup miss (D1) is open (caveat added); "Ag/Cu interface"
  mislabeled the cited R_ct band (it is Ag/YBCO, inherited from the source todo, fixed);
  the WIP factory abort, the default-off `mCreateSideConnectors` flag, and the missing
  geometry-interpolation factory case are now stated; Schnaubelt tag normalized to
  paper9; pre-existing wrong `src/mesh` links to ThinShellFactory fixed. Codex prose
  pass applied as the final step.

## Open items

Phase 2 continues per `todo/hex8tb_phase2_fem_wiring.md`: R4 (geometry interpolation
decision), R5 (IWG + MaxwellFactory wiring; D4/D5 block-data contract first), R6
(calculator/postprocessor), R7 (verification gates incl. numeric handedness check),
Q4 recovery implementation, coated variant (`h_in`), closed-loop tapes (O4/O6).
