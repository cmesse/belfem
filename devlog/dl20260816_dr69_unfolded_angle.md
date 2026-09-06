# DR-69 Implemented: the Unfolded β Reaches the jc Tables

**Date:** 2026-08-16
**Purpose:** Record the plan+audit → code+audit round that landed the
unfolded field-to-normal angle, and the executed evidence
**Module:** fem/kernel, physics/materials, maxwell
**Round:** `tmp/ai_exchange/dr69_signed_angle_plan.md`

## Christian's refined design

Instead of the ruling's original mechanism (material advertises whether
its Jc(θ) is even), `Calculator::bn_angle` now always returns the
UNFOLDED angle θ ∈ [0, π] and each law decides internally. The survey
made this the minimal design: ModifiedKim is even by construction
(cos²/sin² only — no fold code needed, pinned by two new unit tests),
and JcFunctionDatabase already wraps into its angular window. The whole
defect was the one `abs(dot)` upstream.

## What the audits caught (the round earned its cost again)

- **Both refuted F3:** the Database wrap mapped θ = π to 0 — in eval AND
  all three derivative paths. Fixed with `<=` at all four sites.
- **Grok killed C4's width test:** a [−π/2, π/2] window has width π and
  still mis-wraps; and below-range angles EXTRAPOLATE (the "automatic
  clamp" covers T and B only). The guard became a coverage test in the
  JcFunctionDatabase ctor: `mAngleMin ≤ tol && mAngleMax ≥ π − tol`,
  BELFEM_ERROR at setup tier — a folded export now fails loudly at load.
- **Codex found the missed consumer:** user-defined material plugins
  receive the angle raw; the contract changes from [0, π/2] to [0, π]
  and a stale .so gets the new angle silently. Documented as an ANGLE
  CONTRACT block on the class.
- **Grok killed G2-as-deck:** no input.conf path constructs ModifiedKim
  (`file` → Database, constants otherwise) — the gate became a unit
  test, which now exists and passes.
- **Both corrected the old comment's normal story:** n is the outward
  normal of the mid-surface facet's MASTER volume element (master →
  slave), which is also the layer-stack direction — the minus prefix on
  the thinshell `sidesets` key flips both together via `Facet::flip()`.
  The "n points away from the conductor" sentence was a Kim-era
  leftover; the conductor IS the mid-surface.

## The divergence settled by execution

Codex measured a 35.5% jc endpoint difference in sp-ap and demanded
θ=π regression checks; the raw Robinson data says the poles agree to
≤2.6%. Both were right about their layer: the 35% was CONTROL-POINT
level (the register's standing caveat about order-2 spline storage).
Evaluated through `JcFunctionDatabase`, the endpoints agree —
jc(0)/jc(180) = 1.002 at 0.5 T/84 K and 1.027 at 0.3 T/80 K. The
smoothed export did NOT break the physical 180°-periodicity.

## Executed evidence (Christian's build/test authorization, study paused)

- Full rebuild green; **ctest 14/14** including the two new
  `JcFunctionModifiedKim` cases; `check_doc_claims.py` 34/34.
- Scratchpad probe evaluating shipped `sp-ap.hdf5` through
  `JcFunctionDatabase` with unfolded angles: coverage guard passes, and
  **jc(70°)/jc(110°) = 1.460 at 0.5 T/84 K** — the unfolded path now
  expresses the measured asymmetry (matches the original DR-69 audit's
  1.461 spot check and the raw figshare data verified in
  `dl20260816_dr69_gate2_provenance.md`). Before tonight both angles
  read the same folded value.

## Contract and docs

Input-contract pair updated in the same session: the thinshell
`sidesets` sign is now load-bearing physics (jc-lobe selection), stated
in both `doc/input_file_reference.md` and `doc/input_schema.yaml`
(anchor `"bn_angle"`). Deck exposure named by Grok:
`examples/2D_Tapestack` and `examples/2D_Undulator` carry unsigned
sidesets with sp-ap — their normals must be aligned to the measurement
convention before asymmetry-sensitive results are trusted; same for
tapestack3d itself.

## Remaining gates (Christian's runs)

- G1: constant-jc control deck — the J/Jc asymmetry must vanish
  (doubles as the DR-07 tangent regression).
- G3: tapestack3d frame A/B — expect asymmetric jc/ρ changes between
  the ±x tape halves; the physics payoff gate.
- The β Newton tangent stays deliberately unbound (needs the signed
  ∂θ/∂q chain; DR-07 residual, unchanged).

## Closure

Code-verify signed off (Codex, V1–V6, high confidence). Its one generic
edge — the ctor guard tolerates endpoints within 1e-6 while the runtime
wrap compared exactly, so a third-party table with endpoints an ulp
inside [0, π] could swap poles at exactly θ = 0/π — was closed the same
hour: all four wrap sites now delegate to a shared `wrap_angle()` helper
that snaps boundary values within the guard's tolerance and keeps the
±π periodicity beyond it. By-catch fixed: `materials_usage_guide.md`'s
`angleNxB` row described the metals' current-to-field angle and listed a
phantom `angleNdotB` — rewritten to the real (now unfolded) convention.
Re-gated after both changes: rebuild clean, ctest 14/14, probe identical
(θ = π still reads its own node through the helper).
