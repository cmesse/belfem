# Side-Conductivity Bridge: Code Survey, Baseline Correction, Implementation Plan

**Date:** 2026-07-10
**Purpose:** Session record — survey of the side-edge fuse mechanism (main + legacy
`sideconnectors` branch), the plumbing-gap inventory for the edge-wall impedance bridge, a
baseline discrepancy found against the design note, and the phased implementation plan.
Read-only session; no source modified.
**Module:** fem/maxwell (thin-shell PENTA6TS, h-φ)
**Artifacts:** `todo/sideconnector_bridge_plan.md` (the plan, D1/Rn/On tracker),
`tmp/ai_exchange/sideconnector_bridge_survey.md` (full evidence log, ephemeral).

---

## Context

Follow-up to `dl20260706_side_connector_bridge_formulation.md` / `dl20260709_side_connector_rev4_qa.md`.
O1 (trace bookkeeping, BLOCKING) was resolved in discussion (2026-07) to option (c): the minimal
degenerate wall element — corners stay fused to the air trace g = φ₀ − φ₁, interior side edges
fuse into one inner unknown h_in per station, kernel `r′·BᵀB` on [h_in; g] into K. This session
surveyed the code against that design and produced the phased todo.

## Key findings

- **D1 — baseline discrepancy (the headline) [high].** The design note's picture of *today's*
  model ("all side-edge DOFs through the stack fused into one value tied to the air trace") does
  not match the code for multilayer stacks. Verified: only the outermost layers' edges are hung
  on air φ nodes (`cl_MaxwellFactory.cpp:1316-1317,1432-1548`; complete `add_source` call-site
  sweep `:1250,1476,1536,1620,1715`); at the closed side curve both corners resolve to the SAME
  φ pair (un-split curve nodes, `cl_CutFactory.cpp:1546-1598`) — fusing **by common target**,
  not DOF identity. Interior-level side edges are free `edge_h` unknowns, coupled to the corners
  only through each element's ρ/t stiffness and the ghost-facet Nitsche chain, which is
  deliberately broken at the buffer (`cl_ThinShellFactory.cpp:185,1843,1849-1850`).
  Consequences: §2.1 of the design is exactly true only for single-conductor-layer shells;
  today's multilayer model is NOT the clean r′→∞ limit (buffer-adjacent levels are free chain
  ends — possible lateral leak channel, numeric probe R0.3); and Phase 1 inverts — there is no
  interior fuse to release, the work is to CREATE the h_in fuse (existing edge-on-edge T-matrix
  path with cascade, `cl_FEM_DofMgr_DofData.cpp:3673-3713`, carries it). Reported per the §4
  rule; plan is gated on Christian's confirmation.
- **Stack bookkeeping (needed for the fuse):** ThinShellFactory builds N+1 node/edge levels with
  unconditional per-level copies (`:1028-1129,1514-1555`); adjacent blocks share the
  interface-level edge object at non-duplicate levels and hold ghost-coupled `Edges`/`EdgeDuplicates`
  twins at duplicate ones (`:1601-1674`).
- **Recycling inventory:** the two tools the bridge needs most are already on main with zero
  callers — `compute_binomial_vectors` + terminal-curve sign assert
  (`cl_ThinShellFactory.cpp:1965-2027,2016-2021`, the O5 vehicle) and
  `CurveFactory::thin_shell_side_curves` (consumer `cl_CutFactory.cpp:2285`, demo
  `corctest.cpp:50`). The legacy branch contributes only patterns (`create_tangential_edges`);
  its HEX8TS volume geometry is not needed. No input.conf side-connector parsing precedent
  exists (branch enums were assigned programmatically from the sign). No neighbor-curl trap
  code on the branch (`h_penalty` was a normal-jump penalty).
- **Plumbing gaps (contact_impedance_theory §9.1 analogue):** wall needs a new domain type
  (Curve=34 is thermal-only; InterfaceTsCond=30 dead), parser + FEM_Domain + IWG dispatch +
  FieldList cases (`en_DomainType.cpp:97-179`, `cl_FEM_Domain.cpp:28-73`,
  `cl_IWG_Maxwell.cpp:260-565,601-615`, `cl_Maxwell_FieldList.cpp:299-432`); no per-sideset
  scalar coefficient path exists (materials block-indexed `cl_MaxwellFactory.cpp:2292-2329`,
  only global IWG penalty slots) — r′ path is an explicit decision item; **no in-code energy/
  Joule postprocessing exists anywhere** (MaxwellPostprocessor is pointwise fields only) — the
  wall line item P′ = r′·I′² creates the first one.
- **Cut interaction (O2 scope):** cut constraints attach to φ nodes, the h_in fuse to edge
  DOFs — disjoint entity sets, so the fuse cannot rewrite the cut graph [high]; the
  cut-terminates-on-side-curve station case still gets a numeric probe (plan R2.5;
  `cl_CutFactory.cpp:708` carries an unresolved hanging-edges TODO).

## Decisions and status

- Plan written: `todo/sideconnector_bridge_plan.md` — Phases 0 (note correction + probes),
  1 (h_in fuse), 2 (wall entity plumbing), 3 (kernel `h_edge_wall` + energy line item),
  4 (verification ladder: stencil test → energy audit → transmission-line λ sweep → partition
  sweep). O3/O4 kept as tracked deferrals; O5 resolved as-vehicle. Decision items O-A1/O-B1…B5
  logged, not silently chosen.
- **Blocked on:** Christian's confirmation of the D1-corrected baseline (plan §7 Q1) and the
  five questions in plan §7. Codex audit of the survey + plan pending.
- No source modified. Registered in `todo/README.md` and this devlog in `devlog/README.md`.
