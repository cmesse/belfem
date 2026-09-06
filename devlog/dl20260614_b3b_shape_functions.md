# Devlog: B3b FEM Shape Functions / Interpolation Whitepaper

**Date:** 2026-06-14
**Purpose:** Read-only documentation pass for Tier B3b - the interpolation layer
(`src/fem/interpolation/`): nodal Lagrange + edge Nedelec bases, reference->physical mapping, and the
integration-data assembly that evaluates shape functions at quadrature points.
**Module:** src/fem/interpolation

## Context

Source-grounded whitepaper for the basis-function layer beneath FEM assembly. Scope: the interpolation
classes; mesh edge topology (B3.0), quadrature data (B4c, consumed), assembly (B3a), and formulations/
cuts (B3d) are seams only. Read-only; `./archive` and `./nonfree` positively not accessed.

## Method

Parallel read-only subagents (nodal inventory; Nedelec/edge + Heaviside probe; tests; independent
cross-check of the two highest-risk facts), convention-critical files read directly. **Codex was not
available in this harness** - per the brief, an independent general-purpose subagent reader was
substituted (no priming) rather than fabricating a Codex verdict; it agreed on both high-risk facts.

**Tooling note:** the Bash/rg text channel mangled a few words (rendered "curved"->"ni",
"enrichment"->"liment", "CurvedElementChecker"->"lnChecker"). The Read tool gives true text; all quoted
identifiers were Read-confirmed.

## Key findings

- **MISSION-CRITICAL edge order:** end-to-end max usable Nedelec/h-phi order = **2**, capped
  symmetrically on both sides - interpolation has no edge basis above order 2 (TRI6/TET10), and the mesh
  hard-errors at order >= 3 (`cl_EdgeFactory.cpp:273-274`). No dead higher-order edge basis. Honesty
  nuance: order-2 edge elements are documented as **unvalidated proof-of-concept** (`nedelec.md:201-212`)
  and have zero tests, so production-trusted order is effectively 1.
- **Nodal order:** Lagrange covers the full stored range 1-4 with no gap. Correction: **HEX64 is cubic**,
  not quartic; the true 4th-order types are LINE5/TRI15/TET35 only.
- **EF_LINE3 is dead code** (compiled + in num_nedelec_dofs but not in the factory switch, never
  instantiated). **EF_HEX8TS is live-but-deprecated** (removal decided dl20260605, wiring remains).
- **Heaviside enrichment does not exist** anywhere in src/fem (only a phi-jump mention in homology docs);
  multi-material is handled by cohomology cuts + thin-shell reduction. A separate hierarchical interface
  enrichment (Dular 2021) lives in the kernel/maxwell layer and consumes interpolation - a seam.
- **Mapping:** edge functions carry J/invJ/detJ; nodal IntegrationData stores only parametric dNdXi -
  the physical dN/dx = J^-1 dNdXi is applied in the kernel B/C operator (B3a seam). ElementMapper does
  the inverse (physical->natural) map. Curved geometry handled via `Element::is_curved()` function-pointer
  dispatch (`cl_EF_TET10.cpp:110-124`).
- **B4c seam:** `IntegrationData::evaluate_function()` passes each quadrature point (column of mPoints)
  directly to N/dNdXi with no coordinate transformation; volume + facet(master) + slave(orientation-aware)
  paths. Simplex barycentric-vs-cartesian consistency is a B4c-side guarantee, not checked here.
- **Tests:** Lagrange well covered (Kronecker, partition-of-unity, dNdXi/d2NdXi2 FD); facet integration
  tested with orientations. Nedelec/edge family, ElementMapper, physical dNdx, and IntegrationData class
  API have NO active tests (only-dead tests/old/fem test).
- **Other families:** Bernstein (LINE3/TRI6), Hermite C1 (BEAM/PLATE), bubble facet-enrichment (IFG); no
  Raviart-Thomas / H(div) family.

## Output

- Added `tmp/whitepaper/B3b_shape_functions.md` (in-scope files first; confidence tags; family x geometry
  x order x tested inventory; provenance table; mission-critical edge-order + nodal-order verdicts;
  independent cross-check; open questions; assumptions). Pure ASCII.

## Codex audit incorporated (2026-06-14)

Codex ran the drafted audit prompt (recorded in `todo/ai_exchange.md`) and confirmed the core facts
(max edge order = 2; HEX64 cubic; C1-C5; no Heaviside) while tightening three claims. All three were
re-verified against source by Claude and folded into the whitepaper:
- **A3:** "no dead higher-order edge basis" -> "no dead order-3+ basis" (dead/stub wiring exists at
  order <= 2: EF_LINE3; QUAD9TS/PENTA18TS counts).
- **B1 (most important):** the "full nodal range, no gap" claim was too broad. Quartic shape-function
  classes exist, but `mesh::interpolation_order()` has no QUARTIC branch (returns UNDEFINED for
  LINE5/TRI15/TET35), so `auto_integration_order()` hard-errors for all quartics;
  `interpolation_order_numeric()` hard-errors for LINE5/TRI15; and the enum declares unimplemented
  quintic LINE6/TRI21. Usable nodal order via the standard auto path is effectively 3; only TET35 is
  fully wired at order 4. The "Nodal-order reconciliation" section was rewritten.
- **C6:** quartic LINE5/TRI15 are not in the active derivative test; "well covered" tightened.

These are tightenings, not contradictions - no Claude<->Codex disagreement remained, so no Grok tiebreak
was invoked on these points (Grok's earlier pre-emptive third-voice pass already AGREED on the four
original high-risk claims).

## Verification

- Confirmed `./archive` and `./nonfree` not accessed.
- Confirmed whitepaper is ASCII-only (after Codex-driven edits).
- No source edited, nothing compiled, no tests run (read-only documentation task).
