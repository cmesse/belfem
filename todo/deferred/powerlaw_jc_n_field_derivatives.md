# Power-Law Field Derivatives for Lookup-Table jc and n

> **DEFERRED 2026-09-03** (todo/ currentness sweep, round 3): the β channel is genuinely missing code (`drho_powerlaw_dbeta` has no hits) but is only a Newton-tangent refinement; unscheduled since its DR-69 gate was struck. Status lines and checkboxes below are as they stood at closure and are not maintained.

**Date:** 2026-07-21
**Purpose:** Extend the Newton tangent of the collapsed h-kernels when critical
current density `jc` and power-law exponent `n` come from lookup tables or other
field-dependent `JcFunction`s: `jc = jc(|B|, β[, T])` and
`n = n(|B|, β[, T])`. This adds the analytic `drho_powerlaw_dB` /
`drho_powerlaw_dbeta` family and the required `JcFunction` derivative hooks —
including `deval_dT`, which the THERMAL solver's consistent tangent needs for
the jc(T)/n(T) legs of `drho_powerlaw_dT` (consumer lands with the thermal
collapse; formula recorded in O3).
**Module:** `src/physics/materials` (powerlaws, JcFunction), `src/fem/kernel`
(MaxwellData dispatch), `src/fem/maxwell` (h_newton kernels)
**Status: |B| channel and T-leg both LANDED and COMPLETE ( T-leg finished
2026-08-14 with the Bézier knot motion ); β channel is the only open item.
The gate it was waiting on is gone — see the note directly below.**

> **2026-09-03 currentness sweep (Grok, verified):** this status said the β channel was
> "gated on DR-69's signed-angle work". **DR-69 is struck and archived** — it appears zero
> times in `todo/debt_register.md` and is struck in `todo/debt_register_closed.md`. So β is
> **unblocked but unscheduled**: nothing is holding it any more, and equally nothing is
> pulling it. It needs Christian to either schedule R1–R5 or park this file in `deferred/`;
> a dangling gate is the one state it should not stay in. The knot motion was staged out when the
T-leg first shipped and removed the next day, after the live run stalled
exactly where the staging predicted: a finite-difference check showed the
frozen-knot form alone reproduced only 4-33 % of dρ/dT across the blend, so
the staged term was the dominant one. Round: `tmp/ai_exchange/drhodt_tleg_plan.md`
( phases 1-3 plus a K4 defect round ); record: `devlog/dl20260813_drhodt_tleg.md`. The T-leg shipped via a
second three-AI round the same day ( `tmp/ai_exchange/drhodt_tleg_plan.md` ):
`drho_{powerlaw,piecewise}_dT` + helpers + 8 wrappers + per-law dispatch,
`compute_drhodT_hts` retired ( its reconstruction was sign-flipped and
divergent at the flux-flow crossover ), and `compute_drhodT` now zeroes on a
clamped T iterate like its cp/λ siblings. Knot MOTION in the piecewise blend
LANDED 2026-08-14 (see the header line above and the `Bézier blend: FULL
derivative` comment in `powerlaws.hpp`) — the staged-remainder caveat that
stood here is closed; the derivative is exact to roundoff in the blend. Originally drafted 2026-07-21 (Claude,
from Christian's direction); no code written yet. Held as conditional future work
behind a scope guard until **2026-08-12, when that guard was crossed.** The `tapestack3d`
deck now attaches the measured SuperPower AP characteristic to its HTS material
(`ybco { file : sp-ap.hdf5 ; resistivity type : piecewise ; }`), replacing the
constant `jc`/`n` that made this plan future work. `jc` and `n` are now
`jc(T, log10|B|, θ)` splines, so `mFundRhodB` / `mFundRhodBeta` returning
`return_zero` is no longer harmless: the Newton tangent is inconsistent with the
residual it is supposed to differentiate. This is no longer conditional work.

**What the activation changes, and what it does not.** The converged answer is
unaffected — the residual is the true one, only the tangent is approximate — so
this is a convergence-rate defect, not a wrong-physics defect. The cost is
concentrated exactly where `d jc / d|B|` is steepest, which for the `tapestack3d`
deck is the Ic crossing. Expect Newton to grind or the timestep controller to cut
hard there; that symptom now has a known cause and should not be chased as a new
bug.

**The underlying spline derivatives already exist; everything above them does
not.** `Database` differentiates its own quadratic spline analytically
(`cl_Database.hpp:65-77`): `evaluate_derivx` / `_derivy` / `_derivz` map onto
`sp-ap.hdf5`'s three axes T, log10(B), θ. That removes the need for finite
differencing at the bottom of the chain — but the `JcFunction` hooks, the
`drho_powerlaw_dB`/`_dbeta` family with its defect and piecewise overloads, the
MaxwellData dispatch, and the thermal jc(T)/n(T) leg are all still to be written,
as steps R1-R5 below set out. Two chain-rule factors must not be forgotten when
the hooks are written — the table stores **log10(Jc)** and axis 1 is
**log10(B)**. **CORRECTED 2026-08-13 (both audit voices):** for the B-leg the
ln10 CANCELS against d log10(B)/dB:
`d jc / d|B| = jc * ( d log10 jc / d log10 B ) / |B|` — no ln10. The ln10
survives only on the linearly-stored axes:
`d jc / dT = jc * ln10 * d log10 jc / dT` and likewise for the angle. The
2026-08-12 version of this paragraph carried a spurious ln10 on the B-leg.

**Related, same activation:** `deval_dT` (below) is the thermal leg of the same
gap. On 2026-08-12 `mFundRhodT` was found **unassigned entirely** in the bulk
branch of the MaxwellData dispatch — a null member pointer that the thermal
Newton jumped through (SIGSEGV at address 0). That is fixed separately by binding
`compute_drhodT_metal/_hts/_bulk` in all three bulk cases, but binding the
pointer does not make the handler exact: if `compute_drhodT_hts` differentiates
only the ρ(T) legs and treats `jc(T)` / `n(T)` as constant, the thermal tangent
carries the same inconsistency as the B/β channel. Verify before reading Newton
behaviour at the crossing.

**HISTORICAL — superseded by the 2026-08-12 activation above.** The paragraph below
was written on 2026-08-09, when the guard still held; its "still holds", "harmless
only because current decks use constant jc/n" and "trigger to activate is
unchanged" are **no longer true** and are kept only because they record how the
channel was confirmed dead in tree, and that two independent reviews reached the
same conclusion under the then-correct premise.

> **Re-verified 2026-08-09 (currentness sweep): the scope guard still holds, and the dead
> channel is confirmed in tree.** For HTS materials `mFundRhodB` and `mFundRhodBeta` are still
> hardwired to `MaxwellData::return_zero` (`cl_FEM_Calculator.cpp:163-164`), while ρ itself
> evaluates `jc_eval`/`n_eval` — i.e. exactly the inconsistency this plan exists to remove,
> harmless only because current decks use constant jc/n. Independent confirmation: the
> 2026-08-07 Newton-extension jury raised the same thing as finding **F2** (P1 conditional,
> 3/3 agreement) and downgraded it from prime suspect **on the strength of this file's scope
> guard** (`devlog/dl20260807_newton_extension_jury.md`). Trigger to activate is unchanged:
> the first Jc(B, θ) or Jc(T) tape characteristic attached to an HTS material.
> `debt_register.md` DR-07. Prerequisite machinery in tree: Metal/Database dB/dβ derivatives
> (3-AI-audited 2026-07-21, incl. the HEX27 N[23] fix); consumers
> (`compute_drhodb/dbeta` + E-channel kernel block) shared with the metal/alloy
> tangent work and tracked there.

> ~~**Scope guard — why this is future work:** in the current corc material law,
> jc and n are CONSTANT (Christian, 2026-07-21). Then ρ_PL depends on |J| only,
> `jc_eval`/`n_eval` return constants, and the existing `drho_powerlaw_dJ`
> tangent is already exact — nothing here is needed. This plan activates the
> moment a JcFunction with normB/angle (or T) dependence is attached to an HTS
> material, e.g. a measured Jc(B, θ) tape characteristic.~~
>
> **Guard crossed 2026-08-12** by exactly the case it named: a measured Jc(B, θ, T)
> tape characteristic (`sp-ap.hdf5`) attached to the `ybco` material of the
> `tapestack3d` deck. The premise "jc and n are CONSTANT" no longer holds, so
> `drho_powerlaw_dJ` alone is no longer the exact tangent.

---

## 1. The math

Notation follows `powerlaws.hpp`: Ec critical field, jc critical current density,
n power-law exponent, ρ_n = ρ(T) normal-state resistivity, J = |j|.

```
ρ_PL  = (Ec/jc) · (J/jc)^(n−1)                       ( flux-flow branch )
ρ_eff = 1 / ( 1/ρ_n + 1/ρ_PL )                       ( parallel combination )
```

Partial derivatives of the unclamped flux-flow branch, from
ln ρ_PL = ln Ec − n·ln jc + (n−1)·ln J:

```
∂ρ_PL/∂J   =  (n−1)/J        · ρ_PL        ( existing, drho_powerlaw_dJ )
∂ρ_PL/∂jc  = −(n/jc)         · ρ_PL
∂ρ_PL/∂n   =  ln( J/jc )     · ρ_PL
```

Chain rule to the field variables, with jc = jc(|B|, β) and n = n(|B|, β) from
the JcFunction pair:

```
dρ_PL/d|B| = ρ_PL · [ −(n/jc)·∂jc/∂|B| + ln(J/jc)·∂n/∂|B| ]
dρ_PL/dβ   = ρ_PL · [ −(n/jc)·∂jc/∂β   + ln(J/jc)·∂n/∂β   ]
```

Parallel-combination chain factor. This matches the `_dJ` family because ρ_n
has no B or β dependence:

```
dρ_eff/dX = dρ_PL/dX / ( 1 + ρ_PL/ρ_n )²        for X ∈ { |B|, β }
```

(General form for reference: dρ_eff/dX = (ρ_eff/ρ_n)²·∂ρ_n/∂X
+ (ρ_eff/ρ_PL)²·∂ρ_PL/∂X. The T-derivative — not in scope here — needs BOTH
terms because ρ_n(T), jc(T), and n(T) all move.)

### Guards and consistency rules

- **J → 0:** return 0 for `normJ < BELFEM_EPSILON` (same guard as `_dJ`); this
  also protects the ln(J/jc) singularity.
- **Clamp consistency:** when ρ_PL sits on the `mRhoMin` clamp, the consistent
  derivative is 0. The existing `_dJ` family shares this gap (benign for n > 2
  because ρ_PL ~ J^(n−1) is astronomically small in the clamped regime) — decide
  ONCE for the whole family and document it.
- **Input consistency:** jc, n, ρ_PL, ρ_n, and the chain factor must come from
  the SAME (J, T, |B|, β) tuple as the matching `rho_powerlaw` call. In
  MaxwellData, keep this under the same memoization guard as the
  `compute_rho_powerlaw_* / compute_drhodj_powerlaw_*` pairing
  (`cl_FEM_Calculator.hpp:2806-2850`).

## 2. JcFunction API extension

`JcFunction` (`cl_JcFunction.hpp:134`) today exposes `eval(B, angle)` /
`eval(B, angle, T)` plus the `depends_on(JcParameter)` dependency bitset. The
derivative API should mirror that shape:

- `deval_dB(B, angle[, T])`, `deval_dbeta(B, angle[, T])`, and
  `deval_dT(B, angle, T)`, with default implementations returning 0. A function
  that does not depend on a parameter has zero derivative, matching the
  dependency-bitset policy. The dB/dβ pair feeds the MAGNETIC Newton tangent
  (this plan); **`deval_dT` feeds the THERMAL solver's consistent tangent**
  (∂ρ/∂T through jc(T) and n(T), consumed by `drho_powerlaw_dT` — see O3) and
  belongs to the same API extension so all three derivative hooks land
  together, per subclass, in one audit round.
- **Per-subclass support:**
  - `JcFunction_ModifiedKim` — analytic closed-form derivatives; differentiate
    the Kim expression, no finite differences.
  - `JcFunction_Database` — table-backed derivatives via
    `Database::evaluate_derivx/y/z` (shape-function derivative machinery
    verified 2026-07-21, incl. the HEX27 N[23] deta/dzeta swap fix; check each
    table's axis-to-argument mapping against the chain factors, cf. the
    `inv_element_step` dimension bug B1).
  - `JcFunction_UserDefined` — optional user-supplied derivative callbacks or a
    central-difference fallback with documented step sizes (follow the
    `gFinDiffDeltaB` / `gFinDiffDeltaAngle` pattern from
    `Metal::drhodB_kohler`).
- Angle-domain care: β from `bn_angle` lives in [0, π/2]. Table/Kim derivative
  evaluation at β = 0 and β = π/2 must not fold through the reflection; use
  one-sided differences or analytic edge formulas.

## 3. Material API extension

Mirror the `drho_powerlaw_dJ` overload matrix (`powerlaws.hpp:1027-1212`) for
`drho_powerlaw_dB` and `drho_powerlaw_dbeta`:

- 4-arg `(normJ, T, normB, angleNxB)` primary overload with `jc_eval`/`n_eval`
  dependency routing (O1 policy). For variable X ∈ {normB, angle}, zero the
  independent derivative legs; short-circuit the total derivative only when
  neither `jc` nor `n` depends on X.
- 8-arg defect twin `(…, x, y, z, t)`: the defect modulation
  jc → jc·d(x,y,z,t) also scales ∂jc/∂B by d because d is field-independent.
  ∂ρ_PL/∂jc uses the modulated jc; verify the −n/jc factor uses that MODULATED
  value.
- Piecewise family: the three-regime `rho_piecewise` needs branch-aware B/β
  derivatives, matching `drho_piecewise_dJ` (power-law regime = this plan's
  formula; normal regime = 0; Bezier blend regime = differentiate the blend in
  log-log space with respect to the branch endpoints' B/β dependence). This is
  the hardest piece — own step, own audit.

## 4. Consumers (shared with the metal/alloy tangent work)

- MaxwellData: the HTS branch of the `compute_drhodb/dbeta` dispatcher pair
  selects the powerlaw variants when the `jc` or `n` supplier depends on
  normB/angle; it uses `return_zero` only when both are constant in that
  variable (the current corc case).
- Kernels: the E-channel tangent block in `h_newton_mu0`/`h_newton_mu`
  (∂|B|/∂q and ∂β/∂q rows) is identical for the metal and HTS channels — build
  once, feed either derivative. Tracked in
  `maxwell_kernel_collapse_plan.md` (R12 follow-up note).

## 5. Steps

- [ ] **R1 — JcFunction derivative API**: `deval_dB`/`deval_dbeta`/`deval_dT`
  base + dependency-routed zero defaults; ModifiedKim analytic; Database
  table-backed; UserDefined fallback. `deval_dT` is built here even though its
  consumer (`drho_powerlaw_dT`, thermal Newton — O3) lands later. Unit probe:
  finite-difference vs analytic at random states (mattest-style).
- [ ] **R2 — `drho_powerlaw_dB`/`_dbeta`** 4-arg + defect overloads per §1/§3;
  Codex audit of the derivation (chain factor, modulated-jc factor, guards).
  *(after: R1)*
- [ ] **R3 — piecewise branch-aware B/β derivatives** (own audit round).
  *(after: R2)*
- [ ] **R4 — MaxwellData HTS-branch dispatch** into `compute_drhodb/dbeta`
  (prerequisite: the dispatcher pair from the metal/alloy work exists).
  *(after: R2)*
- [ ] **R5 — verification**: (a) element-level finite-difference Jacobian check
  against the assembled tangent on a one-element probe; (b) convergence A/B on
  a Jc(B,θ)-table corc variant at Δt = 100 ms — the acceptance signature is
  quadratic Newton contraction at relax 1.0 and no Picard→Newton switch
  regression (cf. the 2026-07-21 timestep-36/38 pathology). *(after: R4)*

## 6. Open questions

- [ ] **O1 — isotropized vs exact E-channel row:** the ∂β/∂q row mixes E- and
  C-channel terms; decide exact rank-1 form vs a `phi_ferro`-style isotropized
  approximation (project convention for the mass term; Messe et al. 2023
  (paper1) Sec. 2.7 for the iteration strategy). Needs a short derivation note
  + Codex check BEFORE implementation.
- [ ] **O2 — clamp-consistent derivatives:** zero the whole `drho_*` family on
  the `mRhoMin` clamp, or keep the current benign inconsistency? One decision
  for `_dJ`, `_dB`, `_dbeta` together.
- [x] **O3 — T-derivative for the thermal solver — LANDED 2026-08-13
  ( three-AI round, `tmp/ai_exchange/drhodt_tleg_plan.md` ):**
  `drho_powerlaw_dT` + `drho_piecewise_dT` ( 4-arg + defect overloads ),
  `djc_eval_dT`/`dn_eval_dT` helpers, 8 kernel wrappers, per-law dispatch;
  `compute_drhodT_hts` retired. The two-term formula below is exactly what
  shipped; the *consumer-site paragraph* below it is retained as history but
  its recipe is WRONG — see the strike-through note. Piecewise follows its
  own residual branch-for-branch ( raw dp0dT in PL, dadT in normal/T_crit,
  frozen-knot partial in the blend incl. the (1−t)²·dlnρ1dT term; knot
  MOTION landed 2026-08-14 — the blend tangent is now exact, the old
  weakest-tangent-in-flux-flow caveat no longer applies ).
  Original text: `drho_powerlaw_dT`
  (thermal-coupled Newton) needs the two-term general chain rule — BOTH legs
  of §1's general form, because ρ_n(T), jc(T), and n(T) all move:

  ```
  dρ_PL/dT  = ρ_PL · [ −(n/jc)·∂jc/∂T + ln(J/jc)·∂n/∂T ]
  dρ_eff/dT = (ρ_eff/ρ_n)²·dρ_n/dT + (ρ_eff/ρ_PL)²·dρ_PL/dT
  ```

  The jc/n legs consume R1's `deval_dT`; the ρ_n leg is the existing
  `drhodT` machinery. The consumer belongs to the thermal-tree collapse
  (R11 of the kernel-collapse plan /
  `thermal_matrices_cleanup_and_newton_plan.md`), noted here so the formula
  does not get re-derived wrong and so R1 ships the hook it will need.

  ~~**Concrete consumer site:** `compute_drhodT_hts` … `b = a·c/(c−a)`
  recovers ρ_PL from the memoized ρ_eff — fill `dbdT` there.~~
  **RECIPE REFUTED AND RETIRED ( 2026-08-13 round, Grok #2 + algebra ):**
  the reconstruction was wrong twice — sign-flipped ( the parallel identity
  gives `b = a·c/(a−c)`, so the coded form returned b < 0 and its dadT
  weighting `c²/(a−2c)²` DIVERGED at the flux-flow crossover ρ_PL = ρ_n ),
  and structurally invalid for piecewise ( raw ρ_PL / Bézier ρ_FF are not
  parallel combinations, so no reconstruction exists ). Do not resurrect it,
  including a "fixed" `a−c` variant: the shipped implementation recomputes
  ρ_PL directly and uses the closed factors
  `∂c/∂a = (c/a)²`, `∂c/∂b = 1/(1+ρ_PL/ρ_n)²`. `compute_drhodT_hts` is
  deleted; the eight `compute_drhodT_{powerlaw,piecewise}_{ts,bulk}[_defect]`
  wrappers replace it.

## 7. Literature

- Messe et al. 2023 (paper1), Sec. 2.6-2.7 — power-law formulation, hybrid
  Picard/Newton strategy, tolerance policy.
- Arsenault et al. 2023 (paper3) — anisotropic Jc(B, θ) characteristics in the
  h-φ context (the use case that activates this plan).
- Bathe 2016, Ch. 8 — consistent-tangent requirements for Newton robustness.

---

## Implementation record (2026-08-13, |B| channel landed)

Audit-gated implementation the same night the run exposed the cost: design
pre-registered in `tmp/ai_exchange/jc_derivative_plumbing.md`, audited by Codex
and Grok independently BEFORE any code, findings applied, then implemented.

**Landed — the |B| channel end to end:**

- [x] `JcFunction::deval_dB/_dbeta/_dT` virtuals, default 0 ( exact for
      constants, conservative for ModifiedKim/UserDefined; unlike eval() the
      defaults do not abort — the Newton consumer early-outs on zero ).
- [x] `JcFunctionDatabase` overrides for all three legs, clamp-consistent
      ( a clamped leg returns exactly 0 — the clamped value is constant )
      with the CORRECTED chain factors: B-leg no ln10, T/θ legs with ln10.
- [x] `Material::djc_eval_dB` / `dn_eval_dB` routing helpers ( null-function
      → 0, mirroring jc_eval — the audit's null-check requirement ).
- [x] `Material::drho_powerlaw_dB` ( + defect ): dp0 with the unfloored law,
      parallel factor (1+rhoPL/rhon)^-2 with floored rhoPL — the exact
      drho_powerlaw_dJ convention, per audit C4.
- [x] `Material::drho_piecewise_dB` ( + defect ): consistent with
      rho_piecewise's OWN residual — RAW dp0 in the power-law regime ( no
      parallel factor: rho_piecewise returns raw rhoPL there, Codex C5
      sharpening ), exactly 0 in the normal regime and above T_crit, staged
      0 in the Bézier blend ( boundaries also move with jc; tangent jump at
      j1 documented ).
- [x] 8 MaxwellData wrappers mirroring the compute_rho_* preambles; dispatch
      bound in both HTS branches ( ts + bulk × powerlaw/piecewise ×
      defect/plain ).
- [x] plan-file chain-rule correction ( the 2026-08-12 ln10 error ).

**Deliberately NOT landed, with reasons ( the audits' findings ):**

- [ ] **β channel — BLOCKED, do not bind.** Both auditors independently
      refuted C1 for β: `add_rho_field_tangent` differentiates bj_angle
      ( field–current, the metal Kohler variable ) while HTS uses bn_angle
      ( field–tape-normal, with n fixed ). Binding `mFundRhodBeta` would
      apply the wrong ∂β/∂q rows. Needs a bn_angle derivative path in the
      kernel first: dβ/dq through E only, no C row, n constant per element.
- [x] **T-leg — LANDED 2026-08-13** ( the `b = a·c/(c−a)` sign question is
  ANSWERED: sign-flipped AND divergent at c = a/2, resolved by eliminating
  the reconstruction — see O3 ). ~~deferred pending a formula
      check. Codex flagged `compute_drhodT_hts`'s parallel reconstruction
      `b = a*c/(c-a)` as the NEGATIVE of the positive solution of
      1/c = 1/a + 1/b. Verify or fix
      that before adding the jc(T)/n(T) legs — or replace the reconstruction~~
      with a Material-level drho_*_dT. Not rushed unaudited.
- [ ] **Existing inconsistency found by the round, not fixed here:**
      `drho_piecewise_dJ` applies the parallel factor in its power-law
      branch ( `powerlaws.hpp:1745-1750` ) although `rho_piecewise` returns
      the RAW power law there — the dJ tangent is mildly inconsistent with
      its own residual today. Factor is (1+rhoPL/rhon)^-2 ≈ 1 in the deep
      power-law regime, so the effect is small; recorded for its own fix.

**Status: |B| channel VERIFIED BY EXECUTION ( 2026-08-13, DR-07 closed on
Christian's ruling ).** The recompiled 8-proc tapestack3d run warm-restarted
into the exact regime that defined the gate and crossed 1.30 → 1.44 Ic plus
the full magnetic penetration transient in mostly 2-iterate / 28 s steps —
against 16–50 iterates and 3–14 min per step on the old binary ( ~3x
throughput through the knee, magnetic escalations 5 → 0 ). Still owed: the
constant-jc bit-identity control run ( doubles as the DR-69 asymmetry
discriminator ). The β channel and the T-leg remain open here as R-items,
with DR-69's convention ruling now a prerequisite for the β design.
