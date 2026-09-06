# dρ/dT T-Leg: the Quench-Feedback Tangent Lands

**Date:** 2026-08-13
**Purpose:** Record the three-AI round that implemented the temperature
derivative of the HTS resistivity laws — the term whose absence stalled the
thermal Newton the moment the quench study produced real heating
**Module:** physics/materials, fem/kernel

## The trigger

The tapestack3d run hit its first thermal-Newton stall at t = 3.875 s —
the very step Joule heating became nonzero. Signature: thermal residual
frozen at −49…−51 dB against a −70 dB target, relaxation collapsing
0.69 → 0.043, magnetic side converged and idling. Three rejections followed
within 25 steps, Δt driven from 50 ms to 3.5 ms. The cause was marked in the
source since the thermal collapse:

    // todo: neglecting dn/dT and djc/dT terms for now
    real dbdT = 0.0 ;

`T_h_newton`'s feedback block `dfdx += Nᵀ·N·(dρ/dT·|j|²)` therefore carried
no superconductor contribution at all — and near 77 K the jc(T) chain gives
dρ/dT ≈ 2ρ per kelvin, which IS the quench feedback.

## Two more defects found during planning

**The retired reconstruction was wrong twice.** `compute_drhodT_hts`
recovered the power-law branch from the memoized parallel value as
`b = a·c/(c−a)`; the parallel identity gives `b = a·c/(a−c)`, so the coded b
was negative, and the dadT weighting it fed became `c²/(a−2c)²` — divergent
at the flux-flow crossover ρ_PL = ρ_n — instead of the correct `(c/a)²`.
Grok additionally showed the reconstruction is structurally invalid for the
piecewise law (raw ρ_PL and Bézier ρ_FF are not parallel combinations, so no
reconstruction exists). Both auditors verified the algebra independently.
Grok also correctly bounded the claim: the live stall is the missing dp0/dT
in the power-law regime; the pole is a second bug that would have fired
later, at the crossover.

**The clamp guard was inconsistent** (Codex): `compute_drhodT` zeroed on
`mRhoClamped` only, while its siblings `compute_dcpdT`/`compute_dlambdadT`
also honor `mTClamped`. A Newton must not push on a flat clamp.

## What shipped (round: `tmp/ai_exchange/drhodt_tleg_plan.md` + `_diff.patch`)

Exact mirror of the |B| channel from the same morning, plus the T-specific
deltas:

- `djc_eval_dT` / `dn_eval_dT` helpers (powerlaws.hpp; conservative-zero
  caveat for 3-arg UserDefined without a `deval_dT` override).
- `Material::drho_powerlaw_dT` (4-arg + defect): closed parallel factors
  `∂c/∂a = (c/a)²`, `∂c/∂b = 1/(1+ρ_PL/ρ_n)²`, floored value in factors,
  unfloored law differentiated, and — unlike the dB channel — NO
  `(djc==0 && dn==0)` early-out, because the ρ_n term is a constant-jc
  material's only correct T-dependence.
- `Material::drho_piecewise_dT` (4-arg + defect): branch-for-branch
  consistent with `rho_piecewise`'s own residual — dadT above T_crit and in
  the normal regime, raw dp0/dT in the power-law regime, and in the Bézier
  blend the frozen-knot partial
  `ρ_FF·[(2t−t²)·dadT/ρ_n + (1−t)²·dlnρ1/dT]` with
  `dlnρ1/dT = −djcdT/jc + 2.5·ln10/n²·dndT` (proposed by Grok in review;
  re-derived independently by Codex in the audit). Knot MOTION stays staged
  out — the tangent is weakest in flux-flow until that lands, and a future
  stall at J ≈ j1 is the staged remainder, not a regression.
- 8 kernel wrappers `compute_drhodT_{powerlaw,piecewise}_{ts,bulk}[_defect]`
  with preambles copied from the dB family; per-law dispatch in both HTS
  branches; `compute_drhodT_hts` deleted.
- D6b: `compute_drhodT` guard is now `(mTClamped || mRhoClamped)`.

Deliberately unchanged: the β channel (bn_angle vs bj_angle, DR-69 gate),
the folded angle convention (DR-69's signed ruling is separately gated), and
the known `drho_piecewise_dJ` PL-branch inconsistency.

## Process and verification

Full two-phase protocol: plan pre-registered with seven questions → Codex
and Grok reviewed in parallel (both approved; Codex added the clamp guard
and the UserDefined caveat, Grok corrected the stall attribution, hardened
the staging language, killed the todo's O3 reconstruction recipe, and
contributed the ρ1 blend term) → implementation → both audited the diff.
Codex phase-3: PASS on all seven audit items at high confidence, including
an independent re-derivation of `dlnρ1/dT` and a stale-index hunt on the
new clamp guard. Grok phase-3 in the same thread.

Status: **reviewed, not verified** — all six syntax gates pass (three TUs ×
asserts/NDEBUG under the tree's own flags), nothing built or run. The live
gate is the running quench study itself: after the rebuild, the thermal
Newton should hold 2–4 iterates through the regime that stalled at 10+
tonight, and the rejection cadence at the knee should collapse the way the
magnetic side's did when its |B| tangent landed.

Consumer note for whoever reads the Jacobian later: the T channel feeds
`T_h_newton` only. Magnetic kernels consume dJ/dB/dβ and are bit-unchanged;
`T_h_picard` uses no tangent and is bit-unchanged; converged states are
identical everywhere (residuals untouched) — this changes Newton's rate,
not its answer.

---

## Addendum (2026-08-14): the knot motion lands — and it was the dominant term

The staging documented above ("expect the tangent to be weakest in
flux-flow") was confirmed by the run within hours: the thermal Newton
descended cleanly in the power-law regime and then bounced and froze at
−60 dB once elements crossed j1. This addendum removes the staging.

**What was missing.** `rho_piecewise`'s Bézier blend is built from knots that
are themselves functions of temperature — `j1 = jc·10^{2.5/n}`,
`rho1 = (Ec/jc)·10^{2.5(n−1)/n}`, `j3 = j1·(ρ_n/ρ1)^{1/nff}` and
`j2 = j1·(ρ_n/ρ1)^{1/(n−1)}`. As jc(T) falls the whole transition slides
toward lower current, so an element at fixed J moves *deeper into the
blend*: tParam carries its own T-derivative. The first implementation froze
the knots and differentiated only the control-point values.

**The magnitude, which justified doing this immediately.** Against a central
finite difference at 77 K REBCO constants, the frozen-knot form alone
reproduces only **4–33 %** of dρ/dT across the blend — 4 % just above j1,
rising to 33 % near j3. The staged term was not a refinement; it carried most
of the derivative. That is why Newton stalled precisely and only in
flux-flow.

**What shipped.** With `dlnj1dT = −dlnρ1dT` exactly (both follow from
`j1/jc = 10^{2.5/n}`), the knot chain is

    dlnj3dT = dlnj1dT + ( dlnρ_ndT − dlnρ1dT ) / nff
    dlnj2dT = dlnj1dT + ( dlnρ_ndT − dlnρ1dT ) / (n−1)
              − dndT/(n−1)² · ln( ρ_n/ρ1 )        ← j2 is the only knot whose
                                                    EXPONENT carries n
    dA,dB,dC = ( the corresponding dln combinations ) / ln10
    dS/dT    = ( 2b·dB + dA·c + a·dC ) / ( 2·tSqrt )
    dt/dT    = ( ( dB + dS/dT ) − t·dA ) / a
    knot     = rhoFF · 2(1−t) · ln( ρ_n/ρ1 ) · dt/dT

the last line being a pleasant collapse: the two weight derivatives
(−2(1−t) on log ρ1, (2−2t) on log ρ_n) combine and the ln10 factors cancel.

**A defect the audit caught, and the correction to my own reasoning.**
Codex refuted the guard I first wrote. At `J = j3` the discriminant reduces
to `R²(1/p − 1/q)²` with `p = n−1`, `q = nff`, so it vanishes when `p == q` —
**n = 4 at the default nff = 3**, legal under the `n > 1` precondition, and
exactly where j2 == j3. My `BELFEM_ASSERT( tSqrt > eps )` would have aborted
a legitimate material, and release would have divided by zero. (Noted in
passing: the pre-existing `|a| > eps` assert has its own degeneracy at n = 7,
where a = 0 exactly. Out of scope, same class of trap.)

Grok then refuted my *justification* for the replacement guard. I wrote that
dt/dT is a genuine branch point and the derivative unbounded. It is not:
dt/dT carries 1/tSqrt, but the weight derivative it multiplies is ∝ (1−t),
and at p == q one has a == b, hence 1 − t = −tSqrt/a. **The product is
finite — the singularity is removable.** Verified numerically: the product
converges to 6.9858e-3 as J → j3. So the guard prevents an indeterminate 0/0
in floating point, not an infinite physical derivative, and the comment now
says so.

The fallback also turns out to be exact rather than merely safe: at
tSqrt = 0 the root gives t = 1, where rhoFF = ρ_n and the frozen part reduces
to dρ_n/dT — the normal-branch derivative — so the tangent stays continuous
across j3.

**One more hardening from the re-audit.** Codex noted that `tSqrt` is
computed before the guard, so a discriminant that roundoff pushes slightly
NEGATIVE yields a NaN — and NaN fails every ordinary comparison, so a
`tSqrt <= tol` test would be silently bypassed and the NaN would propagate
into the tangent. The guard now tests the discriminant itself with a negated
comparison, `!( tDisc > tol )`, which takes the fallback for NaN as well.
Demonstrated: with a discriminant of −1e-18 the old form is bypassed and the
new form fires.

**Verification.** A standalone probe `#include`s the live blend branch
verbatim from the edited header (not a transcription) and compares against a
central finite difference: agreement to 1e-9…1e-12 relative for 77 K REBCO
and n = 12 across the blend; all values finite for n = 4 and n = 3.99; the
guard engages at the exact degenerate endpoint and returns 2.0e-9. An h-sweep
confirmed that the only apparent outliers — const-jc points at 1e-17
magnitude — are finite-difference cancellation noise, converging to the
analytic value at h = 1e-4…1e-2.

Both auditors PASS the live tree — Grok on all eight items ("a go for the
weekend rebuild"), Codex on the K4 fix with the NaN caveat now closed.
Status remains **reviewed, not verified**: six syntax gates green, the probe
is a numerical check of the formula rather than a run of the solver. The live
gate is the weekend quench study, where the thermal Newton should now hold
through flux-flow instead of stalling at −60 dB.
