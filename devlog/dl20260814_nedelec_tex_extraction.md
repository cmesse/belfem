# Nédélec LaTeX Notes: Extraction, Drift Record, and the TET4 Edge-Function By-Catch

**Date:** 2026-08-14
**Purpose:** Session record — extraction of the weak-form and Nédélec-derivation chapters of
Christian's pre-BELFEM LaTeX notes (`tmp/nedelec/`) into the module documentation, a drift
record in the notes repository, and the discovery of a confirmed edge-function defect.
**AIs:** Claude (plan, extraction, reconciliation), Codex (two audit rounds + prose pass),
Grok (two audit rounds, refutation role)
**Plan:** `todo/closed/nedelec_tex_extraction_plan.md`

## What landed

Documentation only; no C++ source, no CMake, no `input.conf` keys touched.

- **`src/fem/interpolation/doc/nedelec_derivation.md`** (new): interpolation operators
  N/E/B/C in 3D, 2D and axisymmetric form; barycentric coordinates and Lagrange TRI3/TRI6
  with the transposed-Jacobian pitfall paragraph kept verbatim; Whitney edge functions for
  TRI3/TRI6/TET4/TET10 with the TET table rewritten to the implemented node map (node1↔ξ,
  node2↔ζ, node3↔η); the full TET10 edge/face polynomial set (with the notes' missing `+` in
  `F^11` fixed); the TRI6 circulation-1/2 convention note; the edge/face generation concept.
- **`src/fem/maxwell/doc/maxwell_weak_forms.md`** (new): fundamental lemma, the Gauss and
  Stokes corollaries, least-squares projection (anchored to the live `mt_maxwell_l2_*`
  kernels), the thermal warm-up example, Maxwell equations + MQS assumption, b-conform (a and
  a-v, both marked derived-not-implemented) and the implemented h-conform weak form, with the
  φ-sign convention note (`h ≈ B φ̂` unsigned in the notes; the kernel deliberately drops the
  Arsenault minus on the even mass term).
- **`tmp/nedelec/drift.md`** (new, in the notes repository): what was extracted where; the
  two structural gaps (no cohomology — manual thin cuts vs the automatic thick-then-thin
  pipeline, full treatment deferred to Gregory's forthcoming paper and thesis; no
  hanging-node/hanging-edge concept — λ-interfaces vs condensation); the smaller drifts
  (symmetry penalty vs λ-saddle, antisymmetry removed not replaced, h-a unimplemented,
  Faraday-vs-Gauss + Arsenault 2026 erratum, power-law superset, controller and thin-shell
  evolution); the 18-item errata table.
- Both module-doc README indexes updated; the stale `mt_maxwell_interface`/`mt_maxwell_aphi`
  rows in the maxwell README replaced with the real `matrices/` inventory.

## Method

Plan first (`todo/` per template), then a blind Codex + Grok plan-stage audit, reconciliation,
extraction, then a second blind file-stage audit of the produced markdown, every equation
re-derived by hand against the TeX and against `EF_TRI3/TRI6/TET4`. Prose was polished under
two hard rules from Christian: no em-dashes, and the voice (dry humor, practitioner asides)
stays. All claims are **reviewed, not verified** — static reading and hand algebra, no
executable gate ran.

Notable audit corrections to Claude's initial claims: antisymmetry BCs were *removed*
(hard-errored off), not converted to the penalty form; the cut pipeline is thick-cut-then-thin
(cohomology generates, `CutProcessor` pushes to a thin cut), so automatic generation is the
real drift; TRI3 matches the notes' Whitney forms exactly, no normalization gap; and Grok
refuted the claim that TET10 shares the TRI6 circulation-1/2 convention (the notes' TET10
polynomials are twice the TRI6 pair and integrate to 1).

## The by-catch: two edge-function defects

The E16 convention check (notes' TET edge table vs the code) led Grok to an algebraic
mismatch, which Claude independently confirmed by hand:

- **D1, CONFIRMED (HIGH):** `EF_TET4::E()` edge 2 (`cl_EF_TET4.cpp:193-195`) computes
  `η∇ξ − ξ∇ζ` where the Whitney form is `η∇ξ − ξ∇η`. The code's own comment and the
  separately-coded curl `mC(:,2)` both carry the correct form. Hand-derived consequence:
  circulation 1/2 instead of 1 on its own edge, −1/2 instead of 0 on edge 0 — broken
  tangential conformity wherever the `E` operator is consumed on TET4 (mass matrices, L2
  projections); the stiffness path through `C` is correct, which is presumably why 3D TET4
  results stayed plausible. No test covers edge-function circulation; the planned
  falsification battery (`todo/falsification_tooling.md`) would have caught this in minutes.
- **D2, SUSPECTED (MEDIUM):** the implemented `EF_TET10` edge 0 pair looks like the same
  relabeling slip (Grok hand integral: circulations 2/3 and −1). TET10 is proof-of-concept
  and off the production path.

**Same-day resolution (later in the session):** Christian ruled "I agree with that D1" and
asked for a scratchpad math check first. An exact sympy probe
(`tet4_circulation_probe.py`, session scratchpad) validated its own tooling against
DefElement's published degree-0 N1 tetrahedron basis (all six functions are the Whitney
forms, circulation matrix identity), then reproduced D1 exactly on a random rational tet
(coded: identity except column 2 with +1/2 own / −1/2 foreign; corrected: exact identity;
`mC` equals the curl of the corrected basis in all six columns). Christian applied the
three-line fix in-tree himself (`cl_EF_TET4.cpp:194-196`); D1 is **fixed,
formula-verified** (binary regression run still owed). The TET10 probe
(`tet10_circulation_probe.py`) then **confirmed D2 and upgraded it**: the scalar tables
`mG/mH/mU/mV/mW` in `EF_TET10::precompute()` are transcribed in the notes' generic labels
while the gradient assembly uses the Exodus map, giving 24 edge-dof conformity violations
plus face candidates leaking ±4/3, ±8/3 onto edges; the η↔ζ swap in the scalars alone
restores the exact block-identity circulation matrix and reproduces the notes' polynomials
symbol for symbol. Fix is table-wide (~200 lines incl. nine hand-written derivative tables)
and is pending in `todo/nedelec_edge_function_defects.md` with a probe-first workflow
proposed.

## Exchange

`tmp/ai_exchange/nedelec_tex_extraction.md` (pre-registration, four audit reports, two
reconciliations) — distilled into the plan and this entry; sweepable.
