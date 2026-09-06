# Side-Conductivity Edge-Wall Bridge — Phased Implementation Plan

**Date:** 2026-07-10
**Purpose:** Implement the collapsed edge-wall impedance bridge of
`side_connector_effective_resistivity.md` (rev. 4) with O1 resolved to option (c) — the
minimal degenerate wall element (N=1): corner side edges stay anchored to the air trace
g = φ₀ − φ₁, interior side edges are fused into one inner unknown `h_in` per station, and a
curve kernel assembles `r′·BᵀB` (B = [1, −1] on [h_in; g]) into K on wall-active segments.
**Module:** `src/fem/maxwell` (+ `src/fem/kernel`, `src/mesh`, `src/homology` touchpoints)
**AIs involved:** Claude (survey + plan, 2026-07-10), Codex (audit pending), Christian (decisions)
**Status:** SUPERSEDED (2026-08-08) — never executed; it was blocked on the Phase 0 design-note
correction. The side-connector implementation route became the HEX8TB wall element campaign
(`todo/hex8tb_phase2_fem_wiring.md`): the degenerate wall is now a real degenerate volume element
with its own edge function, not a curve kernel on [h_in; g]. The coated `h_in` variant survives
there as open question O2. Still-valid material from this plan: the D1 baseline finding
(only outermost layers hang on air φ; interior side edges are free unknowns) and the
recyclables inventory (§1). Original status: PLAN — survey complete, blocked on Phase 0
design-note correction; no source modified.

> **Scope guards (from the task brief):**
> - The bridge is **weak-form only**: no cohomology generator, no cut change, no new mesh scale.
>   t_w never appears as quadrature-time geometry — it lives only inside r′.
> - Corner edges stay fused to the air trace exactly as today (I′ = 0 at outer surfaces is the
>   wall's physical BC). Edge-curve φ nodes stay single-valued (closed-side-curve rule survives).
> - Drop the in-plane curl block and the mass term (design note §5.3.1/§2.3); kernel goes into
>   `K()`, no magnetic Newton term. PSD for r′ > 0.
> - Fuse-toggle and kernel are ONE atomic operation per segment — never ship an active (fused)
>   segment without its kernel.
> - Do NOT revive HEX8TS, `h_penalty`, or any neighbor-curl current recovery (trap: neighbor
>   C-matrices yield the plate current, not the wall current — design note §5.3.1).
> - O3 (consistency terms) and O4 (thermal feedback) are tracked deferrals, not scope.

Survey evidence: `tmp/ai_exchange/sideconnector_bridge_survey.md` (2026-07-10, ephemeral —
distilled here and into `devlog/dl20260710_sideconnector_bridge_survey.md`).

---

## 1. Current Behaviour (survey result) and the Baseline Correction

### 1.1 What the code actually does at the side curve [high confidence]

- **Corners:** only the outermost layers' edges are hung on air φ nodes — the thin-shell hang
  loop processes exclusively `blocks().first()` / `blocks().last()`
  (`cl_MaxwellFactory.cpp:1316-1317`, `:1330`, `:1359`; bodies `:1432-1490`, `:1492-1548`;
  `add_source` at `:1476`, `:1536`). Weights are built in
  `DofData::create_dofwise_t_matrices_master()`: LINE2 coefficients [+1, −1] →
  `h_edge = φ₀ − φ₁` (`cl_FEM_DofMgr_DofData.cpp:3628-3638`; LINE3 `:3639-3660`; cascades
  through hanging source nodes `:3569-3579`, `:3615-3621`).
- **At the closed side curve** the air nodes above and below are the SAME φ nodes (closed
  side-curve nodes deliberately not duplicated, `cl_CutFactory.cpp:1546-1598`) → both corner
  edges resolve to the same value g. Fusing **by common constraint target**, not DOF identity.
- **Interior levels are free.** The complete `add_source` call-site list in
  `cl_MaxwellFactory.cpp` is `:1250, :1476, :1536, :1620, :1715` — none touches interior-level
  edges. They are ordinary free `edge_h` DOFs, coupled to the corners only by each PENTA6TS's
  through-thickness (ρ/t) stiffness and, at duplicated levels, ghost-facet Nitsche coupling
  (`hasDuplicates = (tA != tB) && both have rho`, `cl_ThinShellFactory.cpp:185`; ghost creation
  skips non-duplicate levels and buffer blocks, `:1843`, `:1849-1850` — the chain is
  deliberately broken at the buffer).
- **Stack bookkeeping:** N+1 node/edge levels; per-level `new Node`/`new Edge` for every surface
  entity (`cl_ThinShellFactory.cpp:1028-1129`, `:1514-1555`); block b takes bottom edges from
  level l (`EdgeDuplicates` if flagged) and top edges from level l+1 `Edges` (`:1601-1674`) —
  adjacent blocks share the interface-level edge object at non-duplicate levels.

### 1.2 Discrepancy vs the design baseline (task-brief §2.1) — D1

| Assumed (§2.1) | Found | Evidence |
|---|---|---|
| ALL side-edge DOFs through the stack fused into one value | Only the two corners are anchored (to the same g); interior levels are free unknowns | `cl_MaxwellFactory.cpp:1316-1317` + call-site sweep |
| Fusing = the hardwired r′ → ∞ limit | True only for a single-conductor-layer shell; for N_layers > 1 the buffer-adjacent levels are free chain ends, so today's model is NOT the clean r′ → ∞ limit | §1.1 above |
| Phase 1 = "unfuse" interior edges | Inverted: there is no interior fuse to release — Phase 1 must CREATE the fuse (interior edges of both halves → one `h_in`) | `cl_FEM_DofMgr_DofData.cpp:3673-3713` (edge-on-edge T-matrix path exists) |

**Bottom line:** the design's target picture (corners anchored, one interior plateau `h_in`,
jump h_in − g penalized by r′) is *constructible with existing machinery*, but the design note's
description of *today's* model is wrong for multilayer stacks, and the Phase-1 exit test
("all fused ⇒ bit-identical to today") is unsatisfiable as stated. Per the task's §4 rule this
blocks R1+ until Christian confirms the corrected baseline (R0.2/Q1). Physics side-question:
whether today's free interior side traces constitute a real lateral leak channel in production
runs, or the metals' ρ/t stiffness pins them in practice — needs a numeric probe (R0.3).

## 2. Architecture: Why the T-Matrix Fuse + Curve Kernel Is the Right Spine

The two halves of the bridge ride battle-tested engines:

1. **`h_in` via the hanging-DOF/T-matrix path** — the same mechanism that already condenses
   boundary edges onto φ pairs and cascades through multiply-hung sources
   (`cl_FEM_DofMgr_DofData.cpp:3471-3770`). No new solver structure, no saddle rows; assembly
   sees `TᵀBᵀBT` exactly like every thin-shell kernel (design note §5.3.2).
2. **The kernel as an `h_ghost`-shaped K-only sideset kernel** minus the SIPG consistency terms
   and with the physical coefficient r′ instead of a stabilization α
   (`mt_maxwell_h.cpp:1767-1937` is the pattern; `contact_impedance_theory.md` §7 the block
   structure). It is a physical Robin term, NOT a penalty — do not conflate with `h_ghost`
   (normal-direction, between stacked layers; wrong plane, wrong role).

Rejected alternatives (recorded in the design note): resolving/enriching the wall (§2, three
failures + five-mechanism conditioning argument), reviving HEX8TS (unwetted binormal), riding
the cohomology cut (no conductance), h-activating the buffer (contrast disease).

Recycling inventory (survey S2): side-curve identification
(`CurveFactory::thin_shell_side_curves`, consumed at `cl_CutFactory.cpp:2285`, demo in
`corctest.cpp:50`) and the binormal frame + side sign (`compute_binomial_vectors`,
`cl_ThinShellFactory.cpp:1965-2027`, sign assert `:2016-2021`, currently zero callers) are
**already on main**. The legacy `sideconnectors` branch contributes patterns only
(`create_tangential_edges` for curve-element creation); its outer-node/binomial-edge volume
geometry served the wrap element and is not needed. No input.conf side-connector parsing
precedent exists (branch enums were assigned programmatically from the sign).

## 3. Gap Table

Classes: **(a)** missing plumbing (mechanical), **(b)** open design decision (needs Christian),
**(c)** verification/physics probe.

| # | Gap | Class | Evidence |
|---|---|---|---|
| G1 | Design note describes wrong baseline (D1) | (c)→note fix | §1.2 |
| G2 | No fuse of interior side edges into `h_in` | (a) after O-A1 | `cl_MaxwellFactory.cpp` call-site sweep |
| G3 | No wall domain type wired (Curve=34 is thermal-only; InterfaceTsCond=30 dead) | (a) after O-B1 | `en_DomainType.hpp:19-78`, `.cpp:97-179`; `cl_FEM_Domain.cpp:64-67` |
| G4 | No 1D-curve kernel precedent in Maxwell; group/calculator matrices are Ghost-hardcoded | (a) | `cl_IWG_Maxwell.cpp:398-401,555-558,601-615` |
| G5 | No per-sideset scalar coefficient path (materials block-indexed; only global IWG penalty slots) | (b) O-B2 | `cl_MaxwellFactory.cpp:2292-2329`; `cl_IWG.cpp:1073`; `mt_maxwell_h.cpp:1809-1810` |
| G6 | No wall curve elements / sideset creation along side curves | (a) | recyclables in §2 |
| G7 | No in-code energy/Joule postprocessing at all (wall line item has nothing to plug into) | (b) O-B4 | `cl_MaxwellPostprocessor.hpp:22-48`; module-wide grep |
| G8 | Cut-terminates-on-side-curve interaction unverified; `cl_CutFactory.cpp:708` TODO on hanging edges | (c) | survey S3.3 |
| G9 | Wall sideset DOF table entry missing | (a) | `cl_Maxwell_FieldList.cpp:299-432` |

## 4. Ordered Steps

### Phase 0 — Survey write-up, design-note correction, baseline probes (gate for everything)

- [x] **R0.1** Survey written to `tmp/ai_exchange/sideconnector_bridge_survey.md`; this plan
      created; devlog `dl20260710_sideconnector_bridge_survey.md`. *(Claude, 2026-07-10)* — S
- [ ] **R0.2** Update `side_connector_effective_resistivity.md`: record the O1 resolution as
      option (c) (§5.4), and correct §2.1-equivalent text per D1 (corners-anchored /
      interiors-free baseline; r′ → ∞ limit statement; ghost-chain and buffer-break facts).
      Exit: note re-audited by Codex; Christian signs off the corrected baseline. — S
- [ ] **R0.3** Numeric baseline probe (Christian runs): on the corc multilayer case, report the
      interior side-station traces vs g and the per-layer lateral Ampère circulation
      (h_l − h_{l+1} around boundary QUAD4s). Exit: quantified answer to "does today's model
      leak laterally at interior levels?" feeding the corrected note. — M
- [ ] **R0.4** Micro-verification: do buffer blocks' elements allocate edge containers in
      `link_elements_with_edges` (`cl_ThinShellFactory.cpp:1601-1674` iterates all blocks)
      despite the edge-allocation rule (`thin_shell_virtual_domains.md:21-30`)? Read-only
      check; affects how the h_in fuse enumerates "interior side edges of the two halves". — S

### Phase 1 — The `h_in` fuse (after: R0.2 sign-off; atomic with Phase 3 at activation)

- [ ] **R1.1** Decide O-A1 (h_in carrier). Recommended: reuse an existing interior side edge
      (buffer-adjacent level of the lower half) as the carrier DOF and hang the other interior
      side edges on it with weight ±1 via `Edge::add_source(Edge*, weight)` + the existing
      edge-on-edge T-matrix path (`cl_FEM_DofMgr_DofData.cpp:3673-3713`) — no
      `cl_Maxwell_FieldList` change needed. Exit: decision recorded here. — S
- [ ] **R1.2** Segment-wise wall-activation bookkeeping: mark side-curve stations (LINE
      segments) wall-active; default inactive. Data lives with the wall entity of Phase 2
      (creation order: identify side curve → per-station edge stacks → activation flags).
      Files: `cl_ThinShellFactory` (side-curve station → per-level edge lookup),
      `cl_MaxwellFactory` (activation from input). Exit: with no active segments the built
      system is **bit-identical to today** (D1-corrected sense: nothing touched). — M
- [ ] **R1.3** For active segments, install the fuse: per station, all interior-level side
      edges (both halves, both `Edges`/`EdgeDuplicates` twins at duplicated levels — R0.4
      informs the enumeration) become hanging DOFs sourced on the h_in carrier with weight
      matching the binormal-frame edge orientation (O5/R2.4). Sign rule: node-id-based edge
      direction vs curve tangent. Exit: debug assembly shows one free h_in per station;
      T-matrices verified on a 2-station toy mesh; solve with fuse but WITHOUT kernel is
      gated off (atomicity: activation requires r′, enforced at input parsing). — L

### Phase 2 — Wall entity plumbing (parallel to Phase 1 where independent)

- [ ] **R2.1** input.conf: `wall { ... }` subsection inside the `thinshell` section (pattern:
      the `layers` parser, `cl_MaxwellFactory.cpp:2474-2498`): r′ value (Ω·m, with unit
      parsing), segment selection (O-B3), optional `off` default. Exit: parsed struct
      round-trips; rejects wall-without-r′. — M
- [ ] **R2.2** Domain type (O-B1): new enum `ThinShellWall` (recommended over reusing
      Curve=34, which is a thermal-BC concept, or resurrecting dead InterfaceTsCond=30) +
      string entry (`en_DomainType.cpp:97-179`) + `cl_FEM_Domain.cpp:28-73` dispatch. Topology:
      assign explicitly from MaxwellFactory like Ghost (bypasses `cl_Topology.cpp:251-414`
      production — Ghost precedent). Exit: wall sidesets typed end-to-end. — M
- [ ] **R2.3** Wall curve elements: recycle `CurveFactory::thin_shell_side_curves`
      (`cl_CutFactory.cpp:2285` consumer as reference; standalone demo `corctest.cpp:50`) and
      create LINE2 wall elements + sideset along active segments, master-linked to the
      per-station edge stack (pattern: legacy `create_tangential_edges`, branch-only). Exit:
      wall sideset visible in output mesh; station↔element map asserted. — L
- [ ] **R2.4** Orientation (O5): wire `compute_binomial_vectors`
      (`cl_ThinShellFactory.cpp:1965-2027`, first caller since the June removal — keep its
      terminal-curve sign assert `:2016-2021`) to fix the Ampère-loop orientation per wall
      element; store the sign. BᵀB is orientation-proof, but off-diagonal block signs and the
      reported direction of I′ are not. Exit: sign flips mesh-orientation-invariantly on a
      mirrored mesh test. — M
- [ ] **R2.5** Cut-interaction probe (O2 scope, G8): corc case with a cut terminating on the
      side curve + wall segments crossing it; verify the φ-node cut constraint graph is
      untouched (constraints attach to disjoint entity sets — φ nodes vs edge DOFs, survey
      S3.3 [high]) and the T-matrix cascade folds the cut jump correctly through g. Exit:
      constraint-graph dump identical modulo the new h_in rows; solve converges. — M
- [ ] **R2.6** DOF table: wall-sideset case in `collect_sideset_dofs`
      (`cl_Maxwell_FieldList.cpp:299-432`) exposing the station's edge_h (h_in carrier) + the
      φ pair; group matrix allocation for the wall blocks in `cl_IWG_Maxwell.cpp`
      (Ghost precedent `:601-615`, note its first-order PENTA6TS assert `:603`). Exit: DOF
      manager builds wall groups without touching other tables. — M
- [ ] **R2.7** Coefficient path (O-B2): implement the decided r′ storage
      (recommended: `Map<id_t /*wall sideset*/, real>` on the Maxwell IWG, set by the factory;
      per-station r′(y) deferred until O-B5 says otherwise). Exit: kernel reads r′ per group;
      no global-penalty-slot abuse. — S

### Phase 3 — Kernel (after: R1.3 + R2.6; atomic with Phase 1 activation)

- [ ] **R3.1** `maxwell::h_edge_wall()` in `mt_maxwell_h.cpp`: per integration point on the
      wall LINE element, `K_station += r′ · w_k · detJ · BᵀB` with the t̂-projected traces,
      B = [E_in, −E_g] on [h_in; g] — structurally `h_contact_impedance`
      (`contact_impedance_theory.md` §7) with dS → dy, assembled through the T-matrices
      (g resolves to the φ pair automatically). K() only; no M; no Newton term; NO SIPG
      consistency terms (O3 deferred); no neighbor-curl anywhere. Pattern:
      `h_ghost` block loop (`mt_maxwell_h.cpp:1852-1911`) minus the Dm/Ds terms. Exit:
      element PSD check for r′ > 0; assembles on the toy mesh. — M
- [ ] **R3.2** Dispatch: `case DomainType::ThinShellWall → h_edge_wall` in
      `cl_IWG_Maxwell::link_to_group` (`cl_IWG_Maxwell.cpp:260-565`). Exit: end-to-end solve
      on the toy case; r′ → large reproduces h_in → g. — S
- [ ] **R3.3** Energy line item: integrate P′ = r′·I′² over wall elements per timestep and
      report it as an explicit output (O-B4 decides the vehicle — no in-code energy
      postprocessing exists today, survey S3.4, so this creates the first line item; design
      it so bulk terms can join later). Exit: P′ column appears; equals kernel energy
      `qᵀK_wall q` to machine precision. — M

### Phase 4 — Verification ladder (after: Phase 3)

- [ ] **R4.1** Stencil unit test: degenerate wall element vs a resolved-sliver 2D reference;
      agreement to TSA truncation error. No physics driver needed. Exit: automated test in the
      suite. — M
- [ ] **R4.2** Energy-audit regression: input power = stored + bulk Joule + wall Joule to
      solver tolerance; r′ → ∞ reproduces the (corrected-baseline) insulated balance. Needs
      Q2 answered (where the existing balance check lives). Exit: regression test wired. — M
- [ ] **R4.3** Transmission-line benchmark (design note §8.1): λ = sqrt(r′/(R′₁+R′₂)); sweep
      r′ = 1e-8…1e-4 Ω·m, check λ ≈ 0.5 mm…10 cm and convergence order (doubles as the O3
      empirical check); flag the tangential-mesh-must-resolve-λ guideline in the docs. Exit:
      measured λ within tolerance across the sweep; guideline documented. — L
- [ ] **R4.4** Partition sweep + O2 closure: solver-determined substrate current behaves per
      `buffer_cut_topology.md` §3.4, moving smoothly from the insulated limit as r′ decreases;
      exercise the floating-half caveat (terminal-model input option) at least as an xfail
      placeholder. Exit: monotone partition curve; O2 annotated RESOLVED or split out. — L

## 5. Defect Tracker

- [ ] **D1 (HIGH, design-doc):** design baseline mismatch — §2.1 of the task brief (and the
      corresponding language in `side_connector_effective_resistivity.md`) claims all side-edge
      DOFs through the stack are fused to the air trace; code fuses only the corners; interior
      levels are free (evidence §1.2). Found by Claude 2026-07-10 (survey). Fix = R0.2 +
      Christian's confirmation (Q1). Not a code defect — the code may even be *leaky* for
      multilayer stacks (R0.3 decides), in which case a code-side D2 would be opened.

## 6. Open Design Questions

- [ ] **O-A1** h_in carrier: reuse existing interior edge DOF + T-matrix fuse (recommended,
      zero new field plumbing) vs a dedicated new DOF/field label (cleaner postproc name,
      more plumbing). Decide at R1.1.
- [ ] **O-B1** Domain type: new `ThinShellWall` enum (recommended) vs reuse `Curve`(34) vs
      resurrect `InterfaceTsCond`(30). Decide at R2.2.
- [ ] **O-B2** r′ coefficient path: IWG-owned sideset map (recommended) vs mini wall-material
      object vs per-element field. Same design question as the rint interface property
      (`contact_impedance_theory.md` §9.1). Decide at R2.7.
- [ ] **O-B3** Wall segment specification in input.conf: whole side curve per tape (simplest)
      vs gmsh-tagged sub-segments vs coordinate windows. Needs Christian (Q3).
- [ ] **O-B4** Energy line-item vehicle: log table vs hdf5 timeseries vs both — depends on
      where the existing (external?) energy-conservation check lives (Q2).
- [ ] **O-B5** r′(y) spatial variability in v1: constant per wall section (recommended) vs
      per-station values. Needs Christian (Q5).
- [ ] **O2 (inherited, tracked):** single rerouted cut + free substrate coefficient yields the
      solver-determined partition, or Schnaubelt's explicit IC2 needed; floating-half closure
      as input option. Verified empirically at R2.5/R4.4.
- [ ] **O3 (deferred):** consistency terms for the H(curl) Robin trace pair (Juntunen &
      Stenberg). Rationale: Robin conditions are natural; the R4.3 convergence study is the
      empirical gate — add terms only if order degrades.
- [ ] **O4 (deferred):** thermal feedback — route P′ = r′·I′² into the thermal problem.
      Rationale: needs the thermal-coupling milestone; R3.3 already exposes P′ at the Gauss
      points, and if r′ = r′(T) later, dr′/dT·I′(h)·I′(v) belongs to the thermal
      cross-Jacobian, not the magnetic stiffness (design note §5.3.2).
- [x] **O5** Orientation: vehicle confirmed — `compute_binomial_vectors` + terminal-curve sign
      assert survive on main with zero callers (`cl_ThinShellFactory.cpp:1965-2027`,
      `:2016-2021`). Wiring is R2.4. *(resolved as-vehicle 2026-07-10, Claude survey)*

## 7. Questions for Christian

1. **D1 baseline:** confirm the corrected picture (corners anchored to g, interior side-edge
   traces free, ghost/element-stiffness chain broken at the buffer). Was the "all fused"
   description based on single-conductor-layer runs? Do production multilayer runs show any
   lateral-leak symptom (R0.3 will quantify)?
2. **Energy check:** where does the current energy-conservation check live (external script /
   notebook)? Nothing in-repo integrates power (survey S3.4) — the wall line item needs a home.
3. **Segment granularity (O-B3):** is whole-side-curve-per-tape activation enough for the CORC
   cases, or do you need gmsh-tagged sub-segments from day one?
4. **Production stacks:** how many conductor layers per half do current runs use? (Sets the
   urgency of R0.3 and the size of the h_in fuse per station.)
5. **r′ variability (O-B5):** constant r′ per wall section in v1, or spatially varying r′(y)
   (e.g. slit-edge damage profiles) from the start?

## 8. Definition of Done / Audit Trail

Done when R4.1–R4.4 pass in the test suite, the design note carries the corrected baseline and
the O1(c) record, and the λ mesh guideline is documented. Audit trail:
`tmp/ai_exchange/sideconnector_bridge_survey.md` (survey, 2026-07-10, to be swept after Codex
audit), `devlog/dl20260710_sideconnector_bridge_survey.md` (session record),
`devlog/dl20260706_side_connector_bridge_formulation.md` + `dl20260709_side_connector_rev4_qa.md`
(formulation history). Codex audit of this plan: **pending** — priority items S0/D1 and the
R1.3 edge-enumeration claim.
