# Greg CORC convergence regression — three-AI jury investigation

**Date:** 2026-08-06
**Purpose:** Read-only detective work on Gregory's report: convergence on the
Garber CORC model is significantly worse on `sideconnectors` than on `matfix`
(+ his manual uint16→uint32 hack), Newton "does nothing, always reverts to
Picard", and NEW current-density oscillations appear at the gaps between the
tapes. Claude investigation + pre-registration, Codex + Grok blind jury round,
verification + reconciliation.
**Module:** fem/kernel (Controller, ThinShellFactory), deck `cmake-build-debug/greg/`
**Thread:** `tmp/ai_exchange/review_greg_corc_regression.md` (ephemeral; this
file is the distillation). Brief: `tmp/ai_exchange/greg_corc_regression_brief.md`.
**No source modified.**

## Deck facts (hphirun, `iterate_coupled` path)

`algorithm : Picard`, tolerance 1e-10 (relative-only — see dead knob below),
`max iterations : 100`, `target iterations : 50`, Δt₀ = 1 µs, Δt_max = 10 µs,
no `method` key ⇒ BDF1, no `anderson depth` ⇒ Anderson off, magnetic-only.
3 tapes, 1 layer, twist-extruded periodic cell, thin shells with 10 ybco
layers, gaps 0.01 mm, 90 A / 50 Hz transport current, STRUMPACK.

## Root finding: two controller evolution lines that never merged

merge-base `matfix`↔`sideconnectors` = `429913dc` (≈ Jul 27).

- **matfix-only** (Jul 27–28): `8c52161f`, `7e5df34c`, `2371339c` — the
  trace-driven line-search hardening (dl20260727_controller_line_search.md,
  exists only on matfix): damped first Newton entry ×0.5, acceptance-hole fix
  (absolute accept bounded to +1.0 decade of reference), trust growth,
  stagnation-latch ω handling, moved-baseline detector; the D11–D18/N3
  Newton-tangent fix session; CN/Galerkin hard-error; MUMPS/STRUMPACK
  soft-fail contract.
- **sideconnectors-only** (Jul 28–Aug 5): ts17 controller work (watchdog,
  escalation, ω-carry copy, retry hygiene), Anderson infrastructure, Newton
  residual refresh + BDF5 mH(3) restore, and the side-connector campaign
  (`b42044da` … `d3231b7a`).

Neither branch contains the other's fixes. STRUMPACK BLR-off (`f0d4a5a3`) is
shared; the maxwell/thermal matrix producers are near-identical between
branches (assert-level deltas only).

## Confirmed findings (3/3 unless noted; full table in the exchange thread)

1. **Newton is entirely disabled in Greg's deck.** All three Newton entry
   points gate on the CONFIGURED algorithm (`cl_FEM_Controller.cpp:507`,
   `:681`, `:1050`); with `algorithm : Picard` promotion AND the
   Picard-breakdown escalation are dead, so the controller's only answer to
   the HTS E-J contraction loss is cutting Δt.
2. **"Newton reverts to Picard" (his Newton decks):** the ε-switch promotion
   enters Newton at Picard's full ω (`:683-687`, no ×0.5 damped first entry —
   matfix `:584-600` has it) and the handoff is re-evaluated every iteration
   (`:690-697`): one Newton kick above mEpsilonSwitch = 1e-4 demotes straight
   back; the flat-band stall guard additionally latches Picard for the rest of
   the timestep. Scope: escalated Newton (`mForceNewton`) is exempt from the
   ε-demotion (Codex). Exactly the ts34 pathology matfix fixed.
3. **Line-search acceptance hole** at `:750`: the absolute clause
   (`tLogEpsilon < 0.8`) accepts any trial below ~+8 dB regardless of a
   multi-decade regression; matfix `:686-696` bounds it to within 1.0 decade.
4. **`target iterations : 50`** steers Δt to the convergence-basin edge
   (φ = √(target/n), `:1475-1498`; default and recommendation: 20).
5. **Accepted timesteps always satisfy the RELATIVE 1e-10** — no
   unconverged-accept path (`:941`, `run_coupled` `:1770-1784`). The J
   oscillation is therefore NOT classic adjacent-element checkerboarding
   (Messe 2023 paper1, messe2023.txt:548-560 — that pattern is 0/±2jc in
   adjacent elements at loose ε_n; paper wants 1e-11, deck has 1e-10).
6. **Prime structural suspect for the gap oscillations — unconditional 3D
   side-edge fusing (sideconnectors-only, `b42044da`):** with side connectors
   OFF (deck has none; `mCreateSideConnectors = false`), the else-branch at
   `cl_ThinShellFactory.cpp:360-385` still ties the interior layers'
   (l = 1..N−2, i.e. 8 of Greg's 10) duplicated side-edge AND side-node dofs
   to the base layer as hanging sources (`connect_side_edges` `:2928`,
   `connect_side_nodes` `:2999`), for every auto-detected tape side curve
   (`cl_CutFactory.cpp:2300` → `CurveFactory::thin_shell_side_curves` — the
   detection exists on both branches, the FUSING only on sideconnectors:
   `git grep connect_side_edges matfix -- src/fem/kernel/` is empty). This is
   a discretization change at exactly the gap-facing tape edges, on exactly
   the branch where the oscillations appeared. Whether the tie is the correct
   side-wall discretization is the open r′ wall-term question (campaign R5) —
   physics ruling belongs to Christian/Prof. Sirous, not a jury vote.
7. **Free λ₄** (axial-MMF generator, dl20260804): shared code on both
   branches — a conditioning contributor at most, not branch-differential.
8. **EF_QUAD4TS layer-alternation fix (`4a42d982`) is inert here** — 2D-only;
   the 3D PENTA6TS basis has no −mS hack.

New single-raiser findings (Grok; confirmed by verification, need human
adjudication on priority): `mAbsoluteEpsilonTarget` is parsed (`:1912`) but
consumed nowhere — dead knob, termination is relative-only; the watchdog
window 30 can cut legitimately slow Picard crawls.

## Open disagreement → executable gate

Grok ranks controller-driven Δt-thrash first for the oscillations; Codex ranks
the side-edge fusing first. Both agree on the discriminating experiment:

- **G0 (cheapest, decisive):** same sideconnectors binary, deck with
  `adapt timestep : false`, fixed Δt = 1 µs, `target iterations : 20`,
  `tolerance : 1e-11`. Oscillations vanish ⇒ controller/Δt-history; persist ⇒
  formulation (fusing first, then λ₄ probe).
- **G1:** the same fixed-Δt deck on matfix vs sideconnectors (branch physics
  A/B).

## Recommendations (pending Christian's approval, NOT applied)

1. **Port the matfix controller hardening onto sideconnectors** (acceptance
   fix, damped first Newton entry, trust growth, stagnation-ω, latch re-arm,
   solver soft-fail) and reconcile with the ts17/Anderson work — the two lines
   are complementary. Highest leverage for both of Greg's convergence
   complaints; prerequisite before any Newton-knob advice to him.
2. **Interim deck for Greg on sideconnectors:** `target iterations : 20`,
   `tolerance : 1e-11`; keep `algorithm : Picard` and Anderson off for the
   G0/G1 baseline runs (isolate variables first — Anderson is a separate,
   later A/B).
3. **Decide the side-wall physics** (keep/condition/revert the unconditional
   interior-layer side-edge fusing for plain thin shells) after G0/G1.

## Attribution

Codex: escalation-exemption scope on the demotion claim, fusing re-rank with
CutFactory trace, λ₄ dense-row secondary framing, port line-cites. Grok:
"H2 is not causal for the Picard deck" refutation, dead `mAbsoluteEpsilonTarget`
knob, watchdog false-cut risk, iterate_coupled-vs-iterate_magnetic path nit,
Δt-thrash ranking + G0 gate design. Gregory: the field report and the matfix
A/B observation that triggered the archaeology.

---

## Round 2 (same day): side-edge fusing physics — jury round

Christian added an experimental `mFuseEdges` switch (default false,
`cl_ThinShellFactory.hpp:276`); with the fusing OFF the gap-edge oscillations
are GONE. Question put to the jury: is leaving the rim edges free required or
forbidden? Thread: `tmp/ai_exchange/review_side_edge_fusing_physics.md`.

**Physics consensus (3/3, verified against literature):**
- The collapsed t→0 limit demands ONE single-valued, air-anchored rim trace:
  Alves 2022a never duplicates lateral-boundary entities
  (alves2022a.txt:215-217); Alves 2022b duplicates surface entities "while the
  nodes at the extremities are not" (alves2022b.txt:409-411) and φ-anchors the
  outer interfaces (h_N = −∇φ⁺, h_0 = −∇φ⁻, :451-468); Dular 2021 imposes the
  lateral j·n = 0 essentially (dular2021.txt:393-399). No published thin-shell
  formulation keeps free duplicated multi-level rim traces.
- Free rims are therefore NOT "correct physics" — they are BELFEM's
  historically validated high-aspect-ratio approximation (all Norris-class
  validation ran free); fuse-to-g is the continuum target. "Oscillation gone"
  ≠ "free is right": removing a constraint always smooths.

**Implementation findings (citation-verified):**
- The b42044da tie IS attempted fuse-to-g, not pairwise fusing: edge-on-node
  hanging runs through `DofData::create_dofwise_t_matrices_master`
  (cl_FEM_DofMgr_DofData.cpp:3560-3645) with LINE2 coefficients (+1,−1) ⇒
  h_e = φ(n0) − φ(n1) = g. Claude's "missing g anchor" was refuted.
- Coverage l = 1..N−2 is coherent by design: the outer interfaces are
  φ-anchored surface-wide by `hang_thinshell_edges_on_nodes_bottom/top`
  (cl_MaxwellFactory.cpp:1469/:1529). Grok's "internally inconsistent partial
  coverage" weakened to a rim-station probe question.
- **P1 defect candidate (Grok, confirmed as code fact):
  `connect_side_nodes` indexes `aTargetNodes( tOrg->index() )` WITHOUT
  original()-normalization (cl_ThinShellFactory.cpp:3007-3009), while
  `compute_side_edge_indices` normalizes (:2915-2927). Gregory's side curves
  carry cut duplicates (cuts terminate on side curves; periodic cell), so the
  node hang can grab the wrong layer node at cut stations — top suspect for
  the observed element-scale alternation.**
- Systematic dup-vs-org sign flips unlikely in the clean case:
  `Mesh::compute_edge_directions` original-normalizes both sides
  (cl_Mesh.cpp:1821-1850). The doc §6 orientation gate remains unrun.
- Claude's "free rim = wrap short around the buffer" was REFUTED as stated
  (the June 2026 wrap short needed SHARED free edges; independent free rims
  create no metallic term). The residual air-crack coupling question and the
  buffer guard-rail severity are routed to Christian.

**Agreed next executable gate:** hanging-source probe — dump per-station
source rows (target edge, sign, source node ids, weights, h − g residual) for
a 3-layer homogeneous tape and one cut/periodic station with fusing ON.

**Recommendations (ruling = Christian + Prof. Sirous):** keep
`mFuseEdges = false` as default now (validated approximation; Gregory
unblocked); re-land the uncoated fuse as either the complete dof-manager
all-to-g or — structurally safer, both auditors — Alves-style non-duplication
of rim entities in the factory, behind the §6 gates (orientation check,
single-tape analytic benchmark, cut-station test); guard buffer-split stacks
against silent free-rim runs before any quench/current-sharing claims.
