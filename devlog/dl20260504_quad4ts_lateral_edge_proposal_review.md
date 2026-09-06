# Devlog 2026-05-04 — QUAD4TS Lateral-Edge Side-Connector Proposal Review

**Date:** 2026-05-04
**Topic:** Physical assessment of replacing HEX8TS side-connector elements with QUAD4TS-on-lateral-edge in the T̂-N̂ plane, with virtual b̂ thickness
**AIs involved:** Claude
**Claude Confidence:** high (concerns are basis-completeness arguments, not weight/wiring details)
**Literature References:** Messe et al. 2023 (paper1) §3; Schnaubelt 2023 (paper9) §III, Eqs. 6-8; Alves et al. 2022a (paper5) §II.B-C; Alves et al. 2022b (paper6) Table 1

## Summary

Reviewed a proposal floated by the user's professor: replace the failing HEX8TS side-connector elements with QUAD4TS elements, placed on the lateral edges of the thin shell in the T̂-N̂ plane, with virtual thickness in b̂. The advertised benefit was "no new DOFs or BCs."

After initial enthusiasm based on the Schnaubelt 2023 T2TCL TSA analogue, **Prof. Sirous raised two basis-level objections** that we now believe are fatal to the proposal as stated. The QUAD4TS-on-lateral-edge construction cannot represent the b̂ jump that drives the cross-connector current and, even if forced, can only generate current in a non-physical direction.

## Geometric Setup of the Proposal

- QUAD4TS plane: T̂ (along the tape) × N̂ (across the stack), placed on the lateral side surface of the thin shell.
- Virtual thickness: b̂ (the wrap thickness, ≈ `mConnectorWidth`).
- Nédélec edges: only two, both along T̂ (`cl_Element_QUAD4TS.hpp:113-135`):
  - Edge 0 at N̂_min (bottom curve, nodes 0→1)
  - Edge 1 at N̂_max (top curve, nodes 3→2)
- Adjacency: both T̂ edges sit on the H↔φ interface, where φ wraps continuously around the closed lateral boundary of the tape.

## Prof. Sirous's Concerns — Why the Construction Fails

### Concern 1 — DOFs are slaved, no jump representable

Strong h-φ coupling enforces `H_T = -∂φ/∂T̂` on every edge of the H↔φ interface (Messe et al. 2023, paper1, §3 lines 506-510, static condensation as preferred BELFEM choice). Because both Nédélec edges of the QUAD4TS lie on that interface and φ is single-valued around the closed tape boundary, both edge DOFs collapse to the same φ trace.

Consequence: no QUAD4TS basis function can represent an H_T discontinuity across the virtual b̂ thickness. This is a **basis-completeness problem** — not fixable by tuning weights, gram-solves, or facet wiring (which were the suspect items in `devlog/dl20260429_connector_weight_physics_audit.md`).

### Concern 2 — Allowed current is in the wrong direction

The full QUAD4TS basis can only produce `H = H_T(T,N) T̂`. Then

```
curl H = (∂_b H_T) N̂ + (-∂_N H_T) b̂
```

With no b̂-dependence in a 2D element with virtual thickness, the N̂ component drops, leaving `J = -∂_N H_T b̂`. That is current punching laterally out of the tape side. The φ domain closes around the tape with no sink to receive this current — it is non-physical.

### Why the two concerns reinforce each other

- Concern 1 says the model cannot generate **any** non-trivial J from the QUAD4TS DOFs.
- Concern 2 says even if a jump were forced (e.g., by introducing a topological cut or breaking φ continuity), the only allowed J direction is the one that doesn't make physical sense.

Either alone would be serious; together they say the QUAD4TS-on-lateral-edge construction cannot represent the connector physics it was meant to capture.

## Why Schnaubelt 2023 T2TCL TSA Escapes This

Schnaubelt 2023 (paper9) §III, Eqs. 6-8 provides H DOFs on **both** b̂-faces of the collapsed contact layer, so the b̂ jump is explicitly representable and the cross-layer current is `(H⁺_T - H⁻_T)/t_virt`. QUAD4TS has no such two-face DOF structure — both T̂ edges are on a single face.

Reusing Schnaubelt's pattern would require a different element topology: one that carries DOFs on both b̂-faces of the connector (i.e., four independent T̂ Nédélec edges, two at b̂_min and two at b̂_max). At that point we are no longer using QUAD4TS, and the "no new DOFs" benefit that motivated the proposal is gone.

## Open Questions and Follow-Up Directions

1. **Geometric clarification with Prof. Sirous.** Is the QUAD4TS intended to sit somewhere other than the lateral-edge surface where φ is single-valued? If both Nédélec edges land in the H volume (interior of the tape stack), concern 1 partially relaxes. Worth a round-trip before writing the proposal off entirely.
2. **Could topological cuts re-enable a jump?** Adding a cohomology cut along the connector could let φ become multi-valued there (Alves et al. 2022, paper6, Eq. 10). But this re-introduces the cut machinery the proposal was trying to avoid, and concern 2 (wrong current direction) still applies.
3. **What element topology actually represents the connector?** If we want to keep TSA-style virtual thickness and avoid HEX8TS curl-metric blowup (`todo/hex8ts_thin_shell_curl_metric.md`), the natural answer looks like a Schnaubelt-style two-face element with 4 T̂ Nédélec edges. This is a new element class, not QUAD4TS.

## Recommendation

Do not pursue the QUAD4TS-on-lateral-edge construction in its current form. Take Prof. Sirous's two concerns back as the formal counter-argument before sinking implementation effort. If keeping the "no new DOFs" goal is essential, we should look at a two-face TSA element rather than reusing QUAD4TS.

## Files Updated

- `devlog/dl20260504_quad4ts_lateral_edge_proposal_review.md`
- `devlog/README.md`
