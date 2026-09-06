# Devlog 2026-07-28 — Side-Connector Revival: Code Recovery + Fusing Theory Settled

**Date:** 2026-07-28
**Topic:** Restart of the side-connector campaign: branch surgery + recovery of the legacy
factory code, the degenerate 4-virtual-dof wall element concept, and the audited physics
verdict on side-curve trace fusing (two-algorithm design).
**AIs involved:** Claude (primary), Codex + Grok (independent audits)
**Claude Confidence:** high (recovery, fusing physics); medium (buffer-adjacent trace value at
the side curve — probe pending)
**Codex Audit Confidence:** medium-high overall
**Grok Audit Confidence:** high on fusing-essential verdict
**Literature References:** Alves et al. 2022a (paper5), `alves2022a.txt:213–228, 301–315,
414–438, 466–473, 507–512`; Dular et al. 2021 (paper0), `dular2021.txt:386–407`; Schnaubelt
et al. 2023 (paperA), `schnaubelt2023.txt:226–327, 337–351`

## Summary

The side connectors get another go, this time as a **degenerate edge element with four
virtual dofs** (a return to the HEX8TB idea): the red outer dofs alias the fused side trace,
one new internal edge dof `h_in` (blue) is the only new unknown, and the element's stiffness
is the collapsed wall resistance r′ — derivable equivalently from the natural boundary term
`∮ δhᵀ(n×e) dS` and from Ampère's law (concept slides: `tmp/sindeconnectors.pdf`).

Christian confirmed the D1 finding of `todo/sideconnector_bridge_plan.md`: interior-level
side-curve edges are NOT fused today. This is a *lucky accident* — nothing needs untangling;
fuses are pure additions. Decision: **two fusing algorithms**, selected per side-curve
station: (i) uncoated → fuse all station side edges to the air trace g = φ₀ − φ₁;
(ii) coated → fuse interior edges to one new internal dof `h_in`, coupled to g only through
the wall element's r′ stiffness. Both audits confirmed the physics and the software shape
(one shared station enumerator, two fuse-target policies).

## Key Findings

- **Free interior side traces are not "insulation."** The do-nothing condition at the
  unmodeled lateral wall is n×e = 0 — a zero-impedance lateral redistribution path between
  layer edges within each conducting half (the r′ → 0 limit), because no separate wall
  impedance is charged. The uncoated slit edge (per-layer J·n = 0 + tangential-H continuity
  with the single-valued air trace) is an **essential** condition: h_0 = … = h_N = g. Fusing
  is therefore *not* overconditioning (test: a constraint is redundant only if fused and
  unfused solutions coincide — Lagrange-multiplier-on-a-natural-BC analogy). Refinements
  accepted from audit: the trace jump h_l − h_{l+1} is the wall-current functional, not
  per-layer sheet outflow; only the *exterior* path is cost-free — mass and both TSA blocks
  tax the same dofs inside the adjacent virtual layers.
- **All-traces-to-g is the correct uncoated fuse** (unanimous): pairwise h_l = h_{l+1}
  without the g anchor floats a common mode against the air or secretly collapses to g via
  the already-hung corners.
- **Buffer-adjacent traces are φ-hung, not free:** `ThinShellFactory::create_buffers` hangs
  every buffer-element edge on its own nodes (`cl_ThinShellFactory.cpp:1731–1784`), and
  conductor/buffer interface levels are shared edge objects, so each half is pinned at both
  ends. Whether the buffer-side value resolves exactly to g at the side curve depends on the
  node-copy source cascade — needs a debug probe (the side node-tying block at `:1789–1844`
  is commented out).
- **Implementation traps catalogued** (exchange file `tmp/ai_exchange/sideconnector_fusing.md`,
  to be distilled into the plan): multi-edge T-matrix sources are LINE2-only; the source
  branch is chosen from `source(0)` — never mix node and edge sources on one target edge;
  hanging cascades flatten one level only (`cl_FEM_DofMgr_DofData.cpp:4094–4096`);
  allocate-then-add discipline; flag-once (don't re-touch air-hung corners); EdgeDuplicates
  twins must be enumerated; buffer φ-only edges excluded from carrier enumeration; `h_in`
  must be a real dof owned by an active element (the wall element must reference it); never
  hang `h_in` on g (the one genuine overconditioning risk); periodic side curves and
  cut-terminating stations need dedicated regressions.
- **Literature:** no paper states the multilayer lateral-edge fuse explicitly. Alves 2022a is
  closest (lateral entities not duplicated; exterior traces tied to φ; ignoring lateral faces
  of the virtual prism = today's free/short baseline). Dular 2021 treats lateral no-current
  as an essential constraint. Schnaubelt 2023 supports the coated analogy (collapsed contact
  layer + free current coefficient).

## Changes Made

- **Branch surgery:** `sideconnectors` moved from the stale 2026-05-13 backup commit
  (418c0235, preserved in `sc_backup` and `sc_bakup`) onto the `devel` tip (429913dc).
- **Recovered verbatim from `sc_backup` (uncommitted):**
  - `src/fem/kernel/cl_ThinShellFactory.{hpp,cpp}` — `create_side_connectors` + 9
    subroutines (`create_outer_nodes`, `collect_inner_nodes_and_edges`,
    `create_tangential_edges`, `create_binomial_edges`, `create_extra_sideset`,
    `create_side_elements`, `preprocess_binomial_edges`, `create_facet_table`,
    `map_extra_nodes`), `SideLayer`/`SideFacet` structs, members, `create()` call site.
  - `src/mesh/cl_ThinShell.hpp` — side-connector containers/accessors (+82 lines).
  - `src/mesh/en_DomainType.{hpp,cpp}` — `LeftCoating=6`, `RightCoating=7`,
    `InterfaceTsConnector=31`, `to_string` case.
  - Kept devel's newer shared utilities and `mCreateGhostFacets = true`.
- **Not recovered (deliberate):** MaxwellFactory consumer wiring
  (`connect_side_connectors_with_facets`, IWG dispatch, `h_penalty`/`h_tb`/`h_tb_t`,
  FieldList/thermal entries). Known defects in the recovered code are untouched — defect
  analysis deferred by Christian's instruction.

## Open Questions

- R0.3 numeric probe (Christian runs): interior and buffer-adjacent side-trace values vs g
  on a corc multilayer case; size of the proxy wall current (h_l − h_{l+1})/L.
- Cut-terminates-on-side-curve station test; twisted-helix periodic regression for the fuse.
- Exact enumeration of the degenerate element's four virtual dofs (element template design —
  next theory item with Christian).
- Whether the recovered geometric machinery (outer nodes, binomial edges) is kept for the
  new element or slimmed to the collapsed-wall variant.

## Files Updated

- src/fem/kernel/cl_ThinShellFactory.hpp / .cpp (recovery)
- src/mesh/cl_ThinShell.hpp (recovery)
- src/mesh/en_DomainType.hpp / .cpp (recovery)
- tmp/ai_exchange/sideconnector_fusing.md (audit thread: Claude + Grok + Codex + resolution)
- devlog/dl20260728_side_connector_revival.md (this file)
