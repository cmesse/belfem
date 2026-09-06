# Side Conductivity: Bridge Formulation and Problem Statement

**Date:** 2026-07-06
**Purpose:** Session record — problem statement for the HTS-tape side conductivity (current escape
around the φ-only buffer at the plated tape edge), three AI review rounds, and the worked-out
"resistor in an H-based formulation" bridge kernel. Read-only session; no source modified.
**Module:** fem/maxwell (thin-shell PENTA6TS, h-φ)
**Artifacts:** `todo/side_connector_effective_resistivity.md` (rev. 3 — the living document),
`tmp/ai_exchange/side_connector_effective_resistivity.md` (review thread).

---

## Context

Christian's email (2026-07-06, meeting with Mathias Schmidt, LLNL/MFEM): all cut-resolving routes
for the electroplated side copper fail on the 20 µm : 400 µm scale separation — explicit meshing
(conditioning), extra DOFs (intrinsic), Heaviside enrichment (oscillations; ghost-penalty/Nitsche
needs a healthy same-phase neighbor that a uniformly-cut edge row does not have). Original proposal:
fold the copper into an effective anisotropic ρ⊥ on perimeter elements. Physical target: does
quench current escape around the insulating buffer to the substrate side (CORC current sharing)?
Prof. Sirois' extra-edge-DOF idea and Christian's own XFEM+Nitsche idea both doubted by Christian;
doubts confirmed this session.

## Review rounds (all 2026-07-06, same thread)

1. **Codex + Grok on the ρ⊥ proposal:** sound as a lumped shunt only with (a) parallel-*addition*
   conductance with the true σ·A/L geometry + series contact resistance (not the naive 20/400),
   (b) tensor `Cᵀ R C`, `R = ρ∥ I + (ρ⊥−ρ∥) n nᵀ` — never a scalar ρ on `Cᵀ C` — plus a matching
   nonlinear Newton term, (c) calibrated coefficients.
2. **Codex + Grok on Nitsche + buffer:** no-healthy-neighbor argument confirmed; BELFEM ghost
   facets verified inter-layer/normal-only (`cl_ThinShellFactory.hpp:256-262`, no in-plane path);
   the φ-only buffer has **no curl-curl channel** (`cl_IWG_Maxwell.cpp:255-269`,
   `mt_maxwell_phi.cpp`) so ρ⊥ there is a no-op, and a one-layer-thick ρ⊥ cannot bridge the buffer.
   Both independently recommended a lateral **collapsed contact-impedance bridge** between the two
   conducting halves; cut/IC2 interaction flagged for verification.
3. **Fable formulation session (this log):** findings below.

## Key findings (round 3)

- **Physics reframe.** Slitting exposes all layer cross-sections; surround plating forms a U-jacket.
  Quench current enters the jacket through the broad top Ag/Cu interface (R_ct ≈ 10⁻¹²–10⁻¹⁰ Ω·m²,
  distributed over the transfer length), NOT through the 1.5 µm REBCO edge — no 20 µm field
  structure exists that a mesh must resolve. The wall metal contributes ≈10⁻⁸ Ω·m of edge; slit-edge
  *interfaces* plausibly 10⁻⁶–10⁻⁴ Ω·m — interface-controlled, i.e. calibratable parameters, not
  fields. Second candidate path: **buffer pinholes** (NI-magnet transverse-resistance anomaly;
  external literature, unverified). Model must expose both paths as independent dials.
- **Resistor in the h-formulation** (answers "how do we implement a 1D resistor with H-based
  DOFs?"): a resistor is `a_R(h,v) = R·I(h)·I(v)` added to K, where the current `I(·)` is a
  **jump or circulation of H** by Ampère — for a collapsed sheet `K = n̂×[[H]]`; for a collapsed
  wire `I = ∮H·dl = Σ±h_e` (Nédélec DOFs ARE circulations). No discrete curl evaluation, no new
  DOFs. Cuts (ideal circulation source), h_ghost (jump penalty), and impedance kernels are one
  family: `coefficient × functional × functional`.
- **Edge-wall bridge kernel:** `a_wall = ∫_Γedge r′(y)·[[h_∥]]·[[v_∥]] dy` on the side curve, with
  `I′ = [[H_∥]]` = wall transfer current per unit edge length (horizontal Ampère loop at buffer
  height) and `r′ = ρ_Cu·h_path/t_w + ΣR_ct/t_w` (Ω·m). Limits correct (r′→∞ insulated, r′→0
  short); dissipation `r′|I′|²` is the thermal source; no free parameter → the HEX8TS unwetted-DOF
  failure cannot recur (the constitutive law supplies what the wrap's basis could not).
- **TSA role-reversal + correction flag for `contact_impedance_theory.md` [medium-high, Codex
  audit pending]:** collapsing a conducting sheet yields `(ρ/t)[[H_t]]²` (in-sheet current) and
  `ρt|curl_t H_t|²` (through-sheet current, R_ct = ρt). For rint/T2TCL (transfer ⊥ sheet) the
  contact-resistance physics is in the ρt term — the one that note's §3.2/§3.4 drops while calling
  the retained ρ/h term "the contact-impedance term". Schnaubelt 2023 §III keeps both matrices.
  Harmless for the internal Ag/YBCO rint only because its transfer length (~2 µm) is sub-mesh; a
  §3.4-style kernel on a genuine large-R_ct T2TCL would silently lose the radial resistance. Feeds
  into Milestone C regardless of this task. For the edge wall the roles reverse: transfer runs
  *within* the wall sheet → the jump term IS the physics.
- **Inert-bridge trap (key discrete design question):** the kernel only acts if the side-curve
  bookkeeping lets the inside/outside traces differ — wherever the wall exists, the shell-edge↔air
  tangential continuity must become the impedance condition. Candidate vehicles: distinct
  conductor-layer edge traces (rint §7 pattern) or a φ-jump form riding the duplicated-node /
  cut-termination machinery already present at side curves. O1 in the todo.
- **Pinhole variant cannot live on buffer φ DOFs** (`∇×∇φ ≡ 0`); it must couple the adjacent
  conductor traces (Hastelloy↔YBCO) with the ρt-type term.
- **1D analytic benchmark:** transmission-line current-transfer length `λ = sqrt(r′/(R′₁+R′₂))`
  ≈ 0.5 mm (pristine) / 1 cm (r′ = 5×10⁻⁶) / 10 cm (r′ = 5×10⁻⁴) — validates kernel, cut
  bookkeeping, and convergence order at once, and λ is itself the physical answer.

## Decisions and status

- Effective-ρ⊥ **demoted to optional refinement**; the bridge is the primary model. Rejected:
  cut-riding (no conductance), h-activating the buffer (contrast disease; fallback only),
  p-enrichment, XFEM+Nitsche.
- Open items tracked with checkboxes in `todo/side_connector_effective_resistivity.md` §5.4/§8
  (O1 trace bookkeeping, O2 cut/IC2, O3 consistency terms, O4 thermal feedback, O5 orientation;
  V1–V5 validation).
- **Codex audit (round 3, same day): all four claims confirmed.** (1) TSA split verified against
  the printed 1/h + h/6 matrices in Alves 2022a (`alves2022a.txt:397–436`) [high]; (2) the §3.4
  correction confirmed [high], but the "~2 µm harmless" subclaim downgraded to *conditional* —
  `λ_ct = sqrt(R_ct/(R_sq,1+R_sq,2))` spans µm to hundreds of µm against metal sheets [medium];
  (3) wall functional, units, limits verified, sign-convention note added [high]; (4) the
  inert-bridge trap confirmed with hard evidence — closed side curves are deliberately not
  node-split (`cl_CutFactory.cpp:1546–1598`), air-adjacent shell boundary edges are hung on air φ
  node differences (`cl_MaxwellFactory.cpp:1319–1324,1432–1537`;
  `cl_FEM_DofMgr_DofData.cpp:3555–3637`), and the conductor-air interface kernel is disabled
  (`cl_IWG_Maxwell.cpp:398–401`) → **O1 elevated to BLOCKING**: an explicit side-wall trace pair
  (or a duplicated-φ-node jump form relaxing the closed-side-curve rule locally) must exist before
  any kernel is written. λ formula + 0.5 mm–10 cm range verified [high].
- **Q&A addendum (same day):** Christian asked whether the wall enters as (a) an added stiffness
  `∫(∇×δh)'ρ(∇×h) dV` or (b) the natural boundary term `∮ δh'(n×e) dS` with Ohm's law. Answer
  (todo §5.3.1): **same term** — the boundary integral is the variational residue of the volume
  stiffness; derive via (a) (`C_wall → (E⁻−E⁺)/t_w`, symmetric PSD jump–jump blocks), justify via
  (b) (`[[e_t]]=0` + Ampère `K=n×[[h]]` + sheet Ohm). Trap flagged: `e = ρ·C·h` with the
  *neighbors'* C measures the tape's current, not the wall's (only valid for Nitsche consistency
  terms). Project onto the curve tangent only — the other jump channel is the HEX8TS binormal
  ghost. Cohesive-zone-element analogy; buffer = `r′→∞` limit via the function space.
- **Still pending:** external literature verification of the transverse-resistance/pinhole
  measurements; Christian's go-ahead for the design note. No source modified this session.

Attribution: side-connector physics input from the meeting with Mathias Schmidt (LLNL); the
extra-edge-DOF alternative is Prof. Sirois' suggestion (assessed §2 of the todo).
