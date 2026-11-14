# Devlog 2026-09-09 — Thick-cut overlay block: three blind passes and an action plan

**Date:** 2026-09-09
**Purpose:** Planning session (no code). Locate every insertion point for the thick-cut fallback
Christian described — instead of aborting when the thin cut cannot be built, add a hidden block of
copies of the crossed air elements carrying [φ, I] and assembling the Whitney coupling — and deliver
a tickable action plan.
**Module:** `src/homology`, `src/mesh`, `src/fem/kernel`, `src/fem/maxwell`
**AIs involved:** Claude (blind pass, plan), Codex gpt-6-astra/high (blind pass), Grok grok-4.6/high
(blind pass); handoff from the 2026-09-06 session on another machine (`tmp/thickcut/`)

## Summary

Deliverable: `todo/thick_cut_overlay_block.md` — 25-row gap table, steps R0–R10 with live
checkboxes, 17 open questions, a generator-record sketch, and a reconciliation table of what each
pass caught. All `file:line` citations were re-opened in this session. Status: PLAN, pending
Christian's approval; Gregory's "eventual target, not now" ruling of 2026-09-03 stands until lifted.

Process: Claude traced the pipeline from the abort sites outward; Codex and Grok were dispatched
with the same brief and told not to open `tmp/thickcut/` (both confirmed they did not); the handoff
was read afterwards and its R1–R7 absorbed. Contamination note: Claude read the one-page
`HANDOFF.md` before its own pass by accident; the rest of the handoff only afterwards.

## Key Findings

- **The fallback has a hard boundary at the core.** The non-tight abort (three sites in
  `CutProcessor::collect_facets`, `cl_CutProcessor.cpp:384-392, :583-594`) is outside the closed core
  and is the v1 trigger. The non-unit abort is *inside* `Cohomology::clean()`, which both
  constructors call (`cl_Cohomology.cpp:36-40, :58-61, :534-544`; factory `cl_CutFactory.cpp:344,
  :364-366`). The overlay does not need unit coefficients, but it cannot receive such a generator
  without an owner-applied exit. Codex and Grok both found this; Claude missed it. **Corrected
  later the same session (Christian):** this is not a blocker — the abort fires only after SPFA
  has proven that no unit representative exists, on a valid integer cocycle, which is exactly what
  the thick path consumes; the fix is a ten-line status exit behind a flag, AI-drafted,
  owner-applied (O13 reworded), and the |z_e| = 2 gate is unconditional in v1.
- **Two data losses precede the abort:** `collect_coefficients` keeps only ±1 bitsets
  (`cl_CutData.cpp:231-248`) and `collect_elements` keeps only elements with ≥ dim flagged edges
  (`:296-302`); the raw integer map must be snapshotted first.
- **Hidden blocks are unfinished infrastructure:** `Block::mIsHidden` has no setter, the Exodus
  writer never checks it for blocks (`cl_Mesh_ExodusWriter.cpp:114-140`), and the `.bfm` saves the
  bit but `ProtoMesh` restores it for sidesets only (`cl_ProtoMesh.cpp:1425` vs `:1333-1342`).
- **Four block-type tables must all learn the new `DomainType`**, or the run aborts / silently
  skips / doubles the air stiffness: `FieldList::collect_block_dofs` (`cl_Maxwell_FieldList.cpp:288-293`),
  `link_to_group` (`cl_IWG_Maxwell.cpp:286-310`), `select_blocks` (`cl_MaxwellFactory.cpp:773-793`),
  and the switch at `:826-832`.
- **Overlapping generator bands need the I_g–I_h Gram term** (Codex, critical). Per-generator
  blocks miss it; the plan resolves it as one block per generator *set* (O14).
- **Volume coupling alone leaves the conductor trace at −∇φ.** The interface relation
  h_e = φ_n0 − φ_n1 + I z_e is a separate change in `create_hanging_edges_and_facets` and the
  T-matrix edge conversion, whose loop is bounded by `number_of_nodes()`
  (`cl_FEM_DofMgr_DofData.cpp:3562-3575`) — a third source cannot simply be appended.
- **Precedents exist for the two hardest mesh-side pieces:** thin-shell blocks are registered with
  `append( mMesh->elements(), … )` + `add_block` (`cl_ThinShellFactory.cpp:232-246`), dropped from the
  metis graph and given their parent facet's owner (`cl_FEM_Kernel.cpp:226-232, :448-460`).
- **`CutFactory::cuts()` is never populated** (`cl_CutFactory.cpp:2332`), so the `DomainType::Cut`
  loop in `create_cuts_sub_master` (`cl_MaxwellFactory.cpp:1240-1243`) is a no-op (Grok). Noted,
  not acted on.

## Changes Made / Proposed

- New: `todo/thick_cut_overlay_block.md` (plan), this devlog, index entries.
- No source modified. Proposed work is R0–R10 in the plan.

## Open Questions

O1–O17 in the plan. The three that gate the start: O13 (core exit for non-unit generators,
Gregory), O14 (overlapping generators, Christian), O17 (overlay block vs in-place augmentation —
Christian's brief chose the overlay; recorded so the alternative stays visible).

## Files Updated

- `todo/thick_cut_overlay_block.md` (new)
- `todo/README.md` (entry)
- `devlog/dl20260909_thick_cut_overlay_plan.md` (new), `devlog/README.md` (entry)
- `tmp/ai_exchange/thick_cut_plan.md` (ephemeral; Claude query, Codex and Grok entries, resolution)

## Addendum — second iteration: difficulty and timeline

Christian asked for a difficulty grade per checkbox and a timeline calibrated on the devlog rather
than on an uninformed estimate. Added as §4.1–4.3 of the plan. Result: ~20 rounds, of which four
are surgeries Christian leads (R2 mixed-mode skip logic in the thin-cut pipeline, R4 distributor
transport and owner inheritance, R6 Nédélec E on a φ block without edge dofs, R7 the interface
current term in the hanging-edge conversion at both the conductor and the thin-shell site); the
remaining rounds are AI-led with his review; every run gate is his. Calibration rows: the `.bfm`
backend rewrite closed in 9 days, three comparable homology/persistence items closed same-day, the
shell-trace surgery in 2 days, ParMETIS wiring in 4. Informed timeline 11–14 working days (2–3
calendar weeks), serial current-driven TET4 cut in 7–9 days; the uninformed reference is 6–8 weeks.
The one thing the calibration cannot price is a second defect behind the first in phases C–D (L-10).
Correction after Christian's pushback: the non-unit case is v1's whole point, not a blocked gate;
the core change is a ten-line owner-applied status exit (O13 reworded), and the |z_e| = 2 test is
unconditional.
