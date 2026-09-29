# Test-Gap Follow-Ups: The Rulings of 2026-09-22 Turned Into Work

**Date:** 2026-09-22
**Purpose:** Christian's rulings from the Tier 1 test-gap run (the "quiz", `todo/closed/test_gap_tier1_execution.md` §4)
produced a queue of small source changes, one fixture and two physics questions. This plan holds
them so the next session can land each through the standing loop (plan → jury → code → jury) with a
consumer inventory before anything that changes a value or refuses an input.
**Module:** `src/mesh`, `src/fem/kernel`, `src/physics/materials`, `src/io`, `tests/`
**AIs involved:** Claude (plan + orchestration), Codex and Grok (juries), Opus (implementation)
**Status:** OPEN — 5 of 12 steps done (R1, R8, R9, R10, R11); R1, R9, R10 landed 2026-09-29 via `todo/closed/ruled_refusals_cluster.md`; every step ruled by Christian on 2026-09-22 unless marked; R7 and R10 re-ruled 2026-09-24

> **Scope guards:**
> - Each item is its own row: own-files diff, jury, gate in both trees, one commit.
> - Consumer inventory (`memory: feedback_consumer_impact_before_behavior_change`) before D19,
>   D21, D22, D23 and D24 land — they change values or refuse inputs.
> - The cohomology core stays closed; nothing here touches it.

## 1. Steps

- [x] **R1 (D12, DR-162)** — landed 2026-09-29 (`ruled_refusals_cluster.md` R1-R2; guard in both linkers, refusal on a connector without its layer block, entry gate widened, four fixture cases). Original text: Guard the serial side-connector loop in `cl_FEM_DofMgr_BlockData.cpp`
  (`:481`, `:501`) with `key_exists` as the parallel twin does (`:617`, `:660-661`); add a
  side-connector fixture to `tests/fem/test_ThinShellSelection.cpp`. Ruled: partial selection is
  SUPPORTED. Depth astra/high.
- [ ] **R2 (D13, O1 + O14)** One shared supported-type list for `get_top_nodes` /
  `get_bottom_nodes` (`src/mesh/meshtools.cpp:1389-1411`): refuse TRI3, TET4, HEX20, the
  higher-order hexes and quads alike; the T-01 HEX20 test case then pins the refusal. Ruled: HEX20
  joins the refused set. Depth astra/high.
- [ ] **R3 (D13a, O13)** Comment at the QUAD9TS / PENTA18TS extrusion branches of
  `cl_ThinShellFactory.cpp` (`:1507`, `:1600`): built but not assembled — `cl_MaxwellFactory.cpp:1711`
  refuses them. Comment-only, no round.
- [ ] **R4 (D14, O16)** Delete the callerless `CurveFactory::intersect( id, id, DynamicBitset & )`
  (`cl_CurveFactory.cpp:643-690`) — it never unflags before sorting. Depth astra/high.
- [ ] **R5 (D18, O11)** `BELFEM_ASSERT` (debug tier) on the master-facet node count in
  `Mesh::update_facet_nodes` (`cl_Mesh.cpp:648`) plus a regression case in
  `tests/mesh/test_FacetNodes.cpp` that runs in debug. Depth astra/high.
- [ ] **R6 (D19, O7)** Ferromagnets: `set_RRR` includes `rho_mag( T )` so that
  RRR = ρ(293 K)/ρ(0 K) as the header states (`cl_Material_Metal.hpp:117-123`,
  `Ferromagnetic::rho_custom`). Changes iron/nickel resistivity values → consumer inventory first
  (every deck with `iron`/`ferro`/`nickel`), physics check against the RRR data source, then the
  T-10 RRR identity test is extended to the ferromagnets. Depth terra/high.
- [ ] **R7 (D20, DR-168)** **Ruled 2026-09-24 (quiz): a Codex + Grok physics round on Ern et al. 2009 vs the coded operator comes FIRST (the DR-163 shape, magnesia φ-layer case worked explicitly); the fixture is built to what the round says separates the forms.** The unequal-thickness / unequal-ρ two-shell fixture that separates the
  ghost-interface JUMP of outward slopes from the SWIP weighted average; Christian reads the jump as
  intended, Claude's assessment (04:15 entry in the tracker) says the two may coincide up to the
  factor 2 depending on the slave slope's sign convention. The fixture decides; no kernel edit before
  it. Depth astra/high.
- [x] **R8 (D21, DR-170)** — landed `9373e44` 2026-09-23 via `todo/materials_tangent_cluster.md` R6 (exact endpoint limit `mNff·ρ_n/J` in all eight overloads, not a frozen-knot fallback). Original text: `drho_piecewise_dJ`: replace the `BELFEM_ERROR` at the vanishing
  discriminant (`powerlaws.hpp:1839-1842`) by the frozen-knot fallback `drho_piecewise_dT` uses
  (`:2305-2307`). Ruled: `n = 4` is legal; the abort tier in a Newton tangent is the defect.
  Depth terra/high.
- [x] **R9 (D22, DR-171)** — landed 2026-09-29 (`ruled_refusals_cluster.md` R3-R8; rule C for the curve exemption, both input artifacts). Original text: `MaterialFactory`: refuse contradictions, allow agreeing redundancy
  (`builtin : pb40sn60rrr5 ; RRR : 5 ;` passes, `RRR : 100` beside it is refused naming both;
  `builtin` beside `curve` naming a different material is refused). Update
  `doc/input_file_reference.md` and `doc/input_schema.yaml` in the same session (input contract
  rule). Depth terra/high.
- [x] **R10 (D23, DR-172)** — landed 2026-09-29 (`ruled_refusals_cluster.md` R9-R10; refusal in the shared helper, gantry + costheta probe counts zero). **Confirmed 2026-09-24 (quiz).** NARROW refusal only: `BELFEM_ERROR` on
  `is_hanging() && number_of_sources() > 1 && ! is_duplicate()` in the `impose_dirichlet` fallback
  (`cl_FEM_SideSet.cpp:344-350`), naming sideset, node and sources as `Bearing::impose_dirichlet`
  does (`cl_FEM_Bearing.cpp:151-157`); the duplicate branch stays inert and is documented as
  deliberate (a cut duplicate's alias `phi_dup = phi_org + Σ I_k` must not be over-constrained,
  `consumer_impact_dr172.md`). A blanket refusal would fire on `gantry`'s cut anchors. Christian:
  "refusing sounds reasonable, but it might cause a regression" — the narrow form fires on nothing
  in the tree; the gate that settles the medium-confidence rows is a `gantry` or `costheta` run
  counting multi-source hanging nodes per Dirichlet sideset. Depth astra/high.
- [x] **R11 (D24, DR-163)** — ruled A and landed `31737f5` 2026-09-24; physics round RUN 2026-09-23 (`tmp/ai_exchange/review_dr163_physics.md`): both voices recommend the tangent follows the raw value (the Bézier joins the raw law C1 at j1); the literature pair is Duron 2004 + the piecewise header, Messe §2.6 is the Kim pointer; ruling owed. Original text: Physics round (Codex + Grok; Messe et al. 2023 §2.6, Rhyner 1993) on
  which side of the piecewise resistivity carries the parallel combination below `j1`; Christian
  leans to "parallel is right" and asked for a jury. Consumer inventory before any fix: changing
  `rho_piecewise` changes every piecewise residual.
- [ ] **R12 (D15, O3)** Optional: a setter for `mRhoMin` as an EXPERT setting (zero default; HTS
  values of order 1e-17, 1e-16 already too big) with a two-sided floor test. Low priority; if it
  becomes a deck key the input contract rule fires.

## 2. Register links

DR-162, DR-163, DR-168, DR-170, DR-171, DR-172 (ruled rows), plus the by-catch rows the consumer
sweep filed: DR-180 (null sideset B on second-order 2-D curves), DR-181 (`background_phi` reads
`dS` at point 0), DR-182 (`Periodic` stamp removes a kernel-bearing sideset). Not in this plan:
DR-177 and O10 / O18 go to Gregory by `devlog/dl20260922_homology_report_for_gregory.md`.

## 3. Open questions

- **O0** Tier 2 of the triage (48 rows, 13 Astra-only, N-20 … N-30) is the next run; the Astra crosswalk is complete, nothing is untiered (`todo/test_gap_backlog.md` §3a).

- **O1** R6 changes shipped iron/nickel results; which decks are reference runs Christian wants
  re-baselined afterwards?
- **O2** R7's fixture geometry: two layers of unequal thickness and ρ ratio 10, or the production
  magnesia-on-conductor pair (ρ_s → ∞)?

## 4. Audit trail

Rulings: `todo/closed/test_gap_tier1_execution.md` §4 (quiz entries 03:35–04:15, 2026-09-22). Consumer
inventories: `tmp/ai_exchange/consumer_impact_*.md` (ephemeral; distilled into the tracker table).
