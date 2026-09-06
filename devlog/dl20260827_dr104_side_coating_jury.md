# DR-104 jury round: physics adjudication of the side-coating thermal exclusion

**Date:** 2026-08-27
**Purpose:** Three-AI jury audit (blind, Codex + Grok) of Christian's engineering decision to
exclude the copper side-coating walls from the thermal problem. Not a diff review — the
subject was a physics-interpretation brief (`tmp/ai_exchange/dr104_physics_brief.md`);
findings and reconciliation in `tmp/ai_exchange/review_dr104_side_coating_thermal.md`.
**Module:** fem/thermal, fem/maxwell (thin-shell side connectors)

## The claim under audit

Christian (in chat, citing "DR-111" — the register row is DR-104; DR-111 is the unrelated 2D
thin-shell negative-Jacobian abort): the side coatings exist because the MgO buffer does not
conduct current; all layers conduct heat, and the tape cross-section is orders of magnitude
larger than the side-coating cross-section, so the coatings are thermally irrelevant — an
engineering decision.

## Outcome (3/3 agreement on direction; scope wording is Christian's call)

- **Electrical necessity: confirmed** in-model and physically. The rho-less buffer is
  retagged `DomainType::Buffer` (`cl_ThinShellFactory.cpp:2055-2064`) and runs φ-only, so
  the wall is the model's only inter-half current path — matching surround-electroplated
  REBCO tape.
- **Thermal transport and mass: exclusion is sound**, but for a stronger reason than the
  quoted cross-section ratio (which is the *longitudinal* comparison, ~100× both walls,
  +2.4% of tape copper). The load-bearing fact is that inter-half heat crosses the broad
  faces: the MgO face conductance beats the copper rim path by ≥3e2–1e3× (and the layer
  wafers share interface node sets — block b uses `aLayers(b)`/`aLayers(b+1)`,
  `cl_ThinShellFactory.cpp:1624-1641` — so that path is actually assembled). Wall heat
  capacity is a ~1% correction. Grok's refinement: through-thickness Cu↔Cu the rim
  conductance is comparable to the Hastelloy path, but irrelevant — the stack thermally
  lumps through-thickness in ~µs, far below quench timescales.
- **The argument does not cover the dropped Joule source, which is the real DR-104.**
  Wall dissipation is assembled in the magnetic stiffness (`mt_maxwell_h.cpp:318-319`) and
  consumed by no thermal kernel (`cl_ThermalFactory.cpp:71-96` selects only
  Conductor|ThinShell|Ferro|Buffer) — coupled energy is not conserved and loss accounting
  under-counts by the wall share. Magnitude is regime-dependent: with the live pure-copper
  wall (r′ ~ 1e-8 Ω·m) it is a transfer-front-local correction (~0.2 W/front vs ~kW/m
  stabilizer, ≪1% globally); with the design doc's edge-interface band
  (1e-6..1e-4 Ω·m, `side_connector_wall_element.md:42-48`) it can be locally first-order.
  All numbers are scaling estimates — **reviewed, not verified**; the executable gate is
  the row's with/without-wall-share Joule A/B.
- **Recommended disposition (route to Christian):** keep DR-104 open at P2, narrowed to
  source-term deposition — lump the wall's `∫ρ|j|²dV` onto the existing seam/tape thermal
  nodes; do **not** grow it into thermally meshing the walls. The seam-temperature
  machinery for r′(T) already exists (`MaxwellData::compute_T_side_connector`,
  `cl_FEM_Calculator.cpp:685-722`, `mHaveSeamT`).

## Corrections the jury made to the brief (both auditors, independently; verified)

- The brief conflated the broad Ag/YBCO face contact resistance (1e-12..1e-10 Ω·m²,
  `contact_impedance_theory.md` §10.1 — a theory note, no implementation) with the r′
  formula's *edge-plating* terms (1e-6..1e-4 Ω·m of edge). Neither is wired into live
  code; the shipped wall resistivity is pure Cu ρ(T,B,β).
- "Walls have no temperature" was overstated: no thermal *dofs*, but seam-interpolated T
  feeds the wall's metallic ρ(T) already; only localized wall self-heating is
  unrepresented.
- Cite fixes: width fallback is `cl_ThinShellFactory.cpp:477`; node sharing is the block
  builder, not `create_nodes_on_layers`; the live geo is
  `cmake-build-debug/tapestack3d/tapestack3d.geo` (4 mm — `more/geo/` variant has 12 mm).

## By-catch

- Stale paragraph in `side_connector_wall_element.md:160-164`: still claims the factory
  ends in the WIP `BELFEM_ERROR`, contradicted by §6 (`:396-398`). Doc-fix candidate,
  not applied (no fixes mid-round).
- Grok challenges the DR-104 row clause "the error grows exactly in the scenarios the
  connectors exist to model": developed current sharing has wall transfer current → 0, so
  the missing heat is transfer-front-local, and the current tapestack3d campaign (periodic
  10 mm slice, soldered stack, slow 2 kA sigmoid, metallic r′) is a weak discriminator.
  Single-raiser, scaling-only — flagged for Christian, not auto-resolved.
- Read-only check after the Grok round: `doc/input_file_reference.md` / `doc/input_schema.yaml`
  changed mid-window (16:33), consistent with Christian's concurrent session (input-contract
  sync), not with either auditor's brief; source-file edits predate dispatch. No breach
  attributed.

No P0/P1 findings. No fixes applied. Next step per protocol §11: the executable A/B gate,
not another review round.

## Follow-up the same day: adiabatic comparison, ruling, doc rename

Christian's follow-up question — how does the neglected wall energy compare with the
error of the adiabatic-tape assumption — settled the disposition. The deck runs with no
thermal boundary condition at all, and no convection BC type even exists in
`cl_ThermalBoundaryConditionFactory.cpp` (Dirichlet/Neumann only), so all deposited heat
is retained. Against LN₂ pool boiling (h ~ 1e2..1e4 W/m²K on ~1.1e-2 m²/m of stack
perimeter, stack heat capacity ~8 J/(m·K), thermal time constant 0.1..1 s vs the 10 s
ramp) the adiabatic assumption mis-states the energy balance by of order the total
dissipated energy — roughly two orders of magnitude more than the wall's ≲1% share, and
conservative in sign where the wall omission is unconservative.

**Ruling (Christian): implementing the coatings into the thermal model has no practical
benefit; the modeling choice is documented and justified in the design doc so it can be
caught and changed whenever the need arises.**

Executed:
- `side_connector_wall_element.md` **renamed** `side_coating_wall_element.md` (git mv;
  Doxygen anchor updated in `doc/doxygen_nav.dox`; campaign page
  `devlog/campaigns/side_connector_wall_element.md` renamed alongside). Live references
  updated (`src/fem/maxwell/doc/README.md`, `nedelec_thinshell.md`, `todo/README.md`,
  `hex8tb_phase2_fem_wiring.md`, `side_edge_fusing_cut_aware_plan.md`,
  `falsification_tooling.md`); dated devlogs and `todo/closed/` archives keep the old
  name as written, and the doc header carries a "formerly" note so those records stay
  decodable. Code identifiers (`h_side_connector` etc.) deliberately untouched.
- New **§6.1** in the doc records the full ruling inline: the four dropped terms with
  their justifications (face-path dominance, ~1% mass, seam-T feedback, front-local
  metallic Joule under the ~100x adiabatic margin), the cosmetic loss-accounting gap,
  the two conjunctive revisit triggers (cooling BC AND edge-interface r′ calibration),
  and the permanent refusal of thermal wall meshing. §6's pending list no longer names
  wall Joule as pending; the stale §3.2 WIP-stop paragraph (the round's by-catch) fixed
  in the same pass.
- **DR-104 re-scoped**: [CODE] → [RULING], P2 → P3, "missing physics" → documented
  conditional deferral pointing at §6.1; strike is Christian's call.
- Codex language sweep dispatched over the revised sections
  (`tmp/ai_exchange/dr104_doc_rename_sweep.md`).

## Closure

**DR-104 struck and archived the same day (Christian: close it — the modeling choice is
documented and justified).** Row moved to `debt_register_closed.md` with the ruling; the
live register's `[P]` count updated (21 → 20). The row's substance survives in
`side_coating_wall_element.md` §6.1: dropped terms with justifications, the accepted
loss-accounting gap, and the two conjunctive revisit triggers (cooling BC AND
edge-interface r′ calibration). Struck is not verified — the magnitudes are scaling
estimates, acceptable because nothing is owed while the triggers do not hold.
