# Consistent Newton Tangent for the Coulomb Gauge Penalty

**Date:** 2026-08-27
**Purpose:** The gauge penalty's Newton cross-term was deliberately skipped when the penalty was
believed vacuous. At current sharing the skipped term is ~(n−1)× larger than the kept one, and the
2026-08-27 gauged tapestack3d_coarse run walls at t = 4.17 s ~~on exactly its signature~~
(attribution RETRACTED post-audit — see §2; the wall was thermal-first). Add the
consistent linearization as rate/consistency hygiene.
**Module:** `src/fem/maxwell/matrices/mt_maxwell_h.cpp`
**Status:** IMPLEMENTED 2026-08-27 late evening, audit-corrected. Christian ordered coding before
the plan audit returned ("while you code…"); the audit landed mid-implementation and its
corrections were applied: q() bound once per element (was per-Gauss-point), scratch `Gq` (vector,
resized by first assignment) + `GtGq` (e×1 matrix) registered in
`create_custom_vectors_and_matrices`, Grok's frozen C7 shape used verbatim. Both TUs compile
`-Wall -Werror` at production flags; in-tree `make` currently blocked by the unrelated libbelfem
CMake refactor mid-flight. Theory doc §7 updated same session (the omission is no longer policy).
OWED: G-FD gate, `make check-fast` once the tree configures.

## 1. Evidence chain (all 2026-08-27)

1. Gauged warm restart helped at a state where chi-off cascaded (Christian's observation).
2. Mechanism: penalty ∝ ρ (campaign form γ = χρ*/μ²) — vacuous sub-critically (campaign verdict,
   correct in its regime), alive at current sharing. Corroborated by the DR-105 four-run
   comparison: `BDF1+gauge/strumpack` was the only zero-solver-failure tape_quench variant.
3. Clean gauged run (`tapestack3d_coarse/out.txt`) walls at t ≈ 4.17 s with magnetic **Newton**
   blowups (−60.6 → −19.9 dB within one iterate) — Picard iterates healthy.
4. Source (`mt_maxwell_h.cpp`): the power-law channel is consistently linearized
   (`dKdx += (Cᵀj)(Cᵀj)ᵀ·drho/|j|`), the gauge channel is not — `K += GᵀG·(χρ)` with **no**
   corresponding `dKdx` block, in both HTS Newton kernels.

## 2. The math

Residual gauge contribution per integration point, q = element dof vector, J = C q:

    R_g = GᵀG q · χ ρ(|J|) · w dV

Consistent tangent:

    ∂R_g/∂q = χ w dV · [ GᵀG ρ   +   (GᵀG q) ⊗ ( dρ/d|J| · Jᵀ C / |J| ) ]
                        ^ kept        ^ SKIPPED — this plan adds it

Magnitude, AUDIT-CORRECTED (both auditors, verified): the (n−1) ratio is the skipped block's
action **on q relative to the kept gauge block only**. Against the PHYSICAL power-law tangent
the skipped term is O(χ·‖Gq‖²/‖Cq‖²) — with this deck's chi = 1e-4 on TET4/PENTA6TS, an O(1e-4)
rate effect, exactly as `coulomb_gauge_penalty_theory.md` §7 already stated ("costs Newton rate
at O(χ), never correctness"). The plan's original "~95 % wrong" sentence over-claimed; the fix
is consistency hygiene and becomes load-bearing only at large chi or on element modes where
Cq → 0 with Gq ≠ 0 (HEX8 kernel modes). **The 4.17 s wall attribution is RETRACTED** (Grok C2,
log-verified): the first cut at that t was the thermal watchdog / PETSc DIVERGED_BREAKDOWN, the
magnetic +40 dB jump came one iterate AFTER the thermal explosion in fully-coupled mode, and the
run recovered at Δt = 0.5 ms with chi unchanged — a Δt-invariant O(χ) Jacobian hole cannot
explain a Δt-dependent thermal-first cascade. Consistent-linearization
requirement: Bathe §8.4 (quadratic convergence needs the full derivative; an inconsistent tangent
degrades to non-convergent iteration when the omitted term dominates).

**The added block is nonsymmetric** ( (GᵀGq)(Cᵀĵ)ᵀ ). MUMPS unsymmetric path handles it — the
deck already runs `linear magnetic { library : mumps }`. Note for the audit: confirm no consumer
assumes K-sector symmetry (STRUMPACK config, conditioning estimator).

## 3. Site inventory (verified on source, 2026-08-27)

| kernel | chi site | gauge K term | gauge tangent | verdict |
|---|---|---|---|---|
| `h_picard` (:117) | :125 | yes | n/a (Picard) | correct as-is |
| `h_newton_mu0` (:156) | :163 | yes | **MISSING j-channel** | **fix here** |
| `h_newton_mu` (:206) | :211 | yes | **MISSING j-channel** | **fix here** |
| `h_side_connector` (:295) | :300 | yes | n/a (Picard-form) | correct as-is |
| `h_side_connector_newton` (:332) | :337 | yes | j-channel ≡ 0 (metal, dρ/dJ = 0) | correct for j; see O2 |

Gauge is ASSEMBLED only in `mt_maxwell_h.cpp` (6 greppable `penalty( 2` hits total — the sixth
is the Controller's rank-0 print; Grok C6). ~~Thin-shell and phi kernels carry no gauge term; out
of scope.~~ CORRECTED (Grok C5): ThinShell groups bind these same kernels, so thin-shell IS
covered by the fix; only the phi kernels genuinely carry no gauge assembly.

## 4. Steps

- [x] **D1** — Derive and freeze the exact update formula per kernel (this file §2), audit-checked
  — and found to MATCH the formula already written in `coulomb_gauge_penalty_theory.md` §7
  (Grok V2): independent re-derivation, same result.
- [ ] **D2** — Implementation shape: inside the existing `if ( tUseGauging )` blocks of
  `h_newton_mu0` / `h_newton_mu`, guarded by `mx->norm_j( k ) > BELFEM_EPSILON` (same guard as
  the power-law tangent), reusing `j = mx->compute_j( k )`:
  `Gq = G * q` (dof vector access: same source `compute_j` uses — resolve the accessor during
  implementation, member scratch per the no-hot-path-allocation rule), then
  `dKdx_times_x() += Gq_col * trans( Ctj ) * ( chi * drho / norm_j * wdV )` where
  `Gq_col = trans( G ) * ( G * q )`. Scratch vectors as members/Calculator matrices, NOT locals
  (hot path). Naming per framework convention.
- [ ] **D3** — Decide O1 (drho reuse) and O2 (metal B-channel) per audit verdicts.
- [x] **R1** — Plan audit round: Codex + Grok jury, returned 2026-08-27 22:14/22:22
  (`tmp/ai_exchange/review_gauge_newton_tangent.md`). Both endorse the fix; neither blocks the
  patch; both refute the causal story and the original D2 shape.
- [x] **R2** — Implemented (order inverted on Christian's instruction; audit corrections folded
  in same session — Grok's C7 shape verbatim). Formal diff audit: Grok's C7 effectively
  pre-audited the exact shipped shape; a separate diff round is Christian's call.
- [ ] **R3 / gates:**
  - **G-FD (THE verification gate — Codex finding 1, adopted):** directional-derivative probe:
    [R_g(q+εv) − R_g(q−εv)]/2ε vs the assembled gauge tangent applied to v, nonzero chi, Gq,
    Cq. Neither G1 nor G2 can verify the linearization (both auditors); this gate can. OWED.
  - **G1 (RECAST per Grok C3/R2):** Christian's evening run is gauged Picard + thermal MUMPS —
    two changes, so it is a machine-health run, NOT a tangent gate. A thermal-watchdog failure
    in it is pre-registered as INCONCLUSIVE for the magnetic tangent.
  - **G2 (RECAST per Grok verdict):** a post-fix Newton traversal of 4.6–5.5 s is progress
    evidence but does NOT verify this linearization (the wall was thermal-first); equally, a
    G2 failure does not indict it.
  - **G3 (no-regression, premise updated by O3):** `chi : 0` decks bit-identical
    (`tUseGauging` gates the block; chi-absent decks now get chi = 1e-4 by design);
    `make check-fast` green. The "sub-critical tangent ~0" claim is dropped (Grok R5: rho is
    small, not zero, sub-critically — a quantitative bound needs the G-FD probe, not prose).
- [x] **O1 — RESOLVED (both auditors + verified):** `compute_drhodj` caches per intpoint
  behind `is_current` (`cl_FEM_Calculator.hpp:2714+`); the kernels' local `drho` is reused,
  no second call, no ordering hazard.
- [ ] **O2 — WIDENED (Grok C8 + Codex 5):** the missing gauge field-channel is
  χ·(GᵀGq)⊗Rwᵀ with the COMPLETE ρ(|B|,β) derivative row Rw, for **HTS and metals alike**
  (HTS binds dρ/dB through the jc(T,|B|,θ) family). Still O(χ); follow-up, deliberately not in
  this patch (theory doc §7 records it as the remaining omission).
> **2026-09-03 currentness sweep — O3's ruling was REVERSED and this file did not record it.**
> `chi` went back to **opt-in on 2026-09-01** (`todo/penalty_opt_in_plan.md`): the Controller's
> absent-block branch now writes `set_penalty( 0.0, 2 )`, not 1e-4
> (`src/fem/kernel/cl_FEM_Controller.cpp:4454-4461`, whose comment states the reversal and its
> reason — on the tape decks the 1e-4 default was invisible to the conditioning estimate and
> kappa tracked the timestep instead). The IWG constructor default **is** still 1e-4
> (`src/fem/maxwell/cl_IWG_Maxwell.cpp:54`), which is exactly why the Controller writes the 0
> explicitly rather than skipping it. **Consequence for this plan: O3's note below that
> "chi-absent decks now get chi = 1e-4 by design" is false, and G3's premise flips back — a
> chi-absent deck is once again an ungauged deck.** Read O3 as history.

- [x] **O3 — RULED 2026-08-27 late night (Christian: "Let's set gauge chi to 1e-4 by default.
  I feel better when it's on").** Landed with the tangent fix in the binaries: ctor default
  `mPenalty(2) = 1e-4`, Controller absent-block branch 1e-4, opt-out `chi : 0`. Input contract
  updated in the same session (`input_file_reference.md` + `input_schema.yaml`, Codex sweep
  dispatched). NOTE this flips G3's premise: bit-identity is now owed for `chi : 0` decks
  (chi-ABSENT decks get the gauge by design).

## 5. Gap table

| gap | status |
|---|---|
| dof-vector accessor in kernel scope | open — resolve in D2 (read `MaxwellData::compute_j`) |
| scratch allocation policy for Gq | open — member scratch, sized once (coding standard) |
| nonsymmetry consumers | open — audit question (§2) |
| 2D h kernels | ~~none found~~ CORRECTED (Codex 9, Grok C5): the same `h_newton_mu0`/`h_newton_mu` serve 2D, 3D, AND ThinShell groups — the fix covers all of them by construction; G is 4×n in 2D |
