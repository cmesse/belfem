# Thermal Expansion Rebuilt on Heat Capacity

**Date:** 2026-08-19 / 2026-08-20
**Purpose:** Replace the low-temperature half of α(T) with a Grüneisen branch anchored on
c_p(T), across nine materials; fix the defects that surfaced while doing it; add Aluminum and
Chromium.
**Methodology:** `src/physics/materials/doc/thermal_expansion_from_heat_capacity.md`

## What started it

A request to make the elastic properties thermodynamically consistent with α and c_p, via a
fitted Grüneisen parameter. The proposal as first written was degenerate: a single shared γ(T)
multiplying both K and G cancels identically out of ν, so the six Bézier parameters could not
touch the one quantity they were meant to produce, and were redundant with the Wachtman
coefficients on E. The discussion turned instead to using the Grüneisen identity in the
direction where it is predictive — deriving α from c_p rather than fitting γ — and the elastic
work was deferred.

That turned out to matter, because α was badly wrong.

## The defect

Every metal's ΔL/L was a cubic Bézier with its first two control points at equal height, a
constraint imposed to satisfy the third law at T = 0. Near the origin the parametric map is
linear and the curve quadratic, so that constraint forces **α ∝ T**, where physics requires
a₁T + a₃T³. Copper's α(20 K) came out at 3.156e-6 against a reference 0.30e-6.

The fit was not at fault. Between 4 K and 20 K copper's true expansion is 0.00014 % while the
inter-dataset scatter is 0.0099 % — the signal is two orders of magnitude below the noise, so
least squares had no information there and the constraint filled the vacuum. The 0 K datum in
the Touloukian compilation (−0.3298486 %, seven significant figures) is itself inconsistent
with the 4–12 K cluster (−0.335 to −0.336) in a physically impossible direction.

No refit fixes this. With x locally linear in the parameter and y at best cubic, α ∝ T² is the
steepest a cubic Bézier can reach; T³ is unreachable at any control-point placement.

## What replaced it

α = C(T)·c_p(T) below a split temperature, with C from a low-order polynomial in T matched to
the Bézier at the split. c_p supplies the shape, so the Debye crossover, the T³ region and the
correct non-zero electronic slope at the origin all come for free rather than being fitted.

Two properties made this land cleanly. C varies by under 1 % across 0–200 K, so the polynomial
does very little work — it is nearly α = const·c_p. And the split-temperature condition
`dC/dT = 0` is exactly the condition that makes the join C¹ automatically, since at a stationary
point of C the slope match reduces to the value match.

The user proposed both the split rule `min(0.618·θ_D, cap)` and the turning-point guard that
falls back from cubic to quadratic when the cubic would put a stationary point inside the
interval. Both are in use: five materials take the cubic, four the quadratic.

Selected results, α in 1e-6/K:

| | 20 K before | 20 K after | reference |
|---|---|---|---|
| Copper | 3.156 | 0.288 | 0.30 |
| Silver | 7.164 | 1.168 | ~1.2 |
| Aluminum | 4.463 | 0.213 | 0.24 |
| Chromium (50 K) | 6.10 | 0.503 | — |

## The cap moved twice

Set at 200 K on the reasoning that a lower anchor limits how far the branch extrapolates. Then
measured: the `dln(C)/dT` sign change sits at 0.43–0.51·θ_D where it exists, so `0.618·θ_D`
clears it by construction, while a flat cap breaks that relationship above θ_D ≈ 425 K. At
200 K Aluminum cleared the guard by 5.2 K and Chromium failed outright. Raising the cap to
273.15 K costs at most 5 % in cryogenic α — measured, not assumed — and fixes both.

## What the γ diagnostic found

γ = 3αK/(ρ c_v) at the split temperature is logged, not enforced. It is an independent cross
check on four separately fitted quantities and reads O(2) for every metal. It found four
defects that nothing else did:

- **Lead** — `mYoungPoly`'s T³, T² and T terms sum to ~45 kPa across 0–600 K against a
  19.888 GPa constant. E(T) is frozen at its 0 K value. Not fixed; elastic work is deferred.
- **Indium** — α(293 K) and γ both 27 % low. The expansion curve was refitted against the
  23-point TPRC dataset: rms 0.0846 % → 0.00525 %, α(293) 23.4e-6 → 29.5e-6, matching the
  finite differences of the data itself. The guard flipped from −1.41e-3 to +2.40e-4.
- **Nickel** — `ref_density` was 890 kg/m³ instead of 8900, giving γ = 20.76. Fixed; the
  corrected value gives 1.73. This affected the thermal mass matrix, not just the diagnostic.
- **Iron** — `create_alpha()` was defined and declared but never called. Iron had no thermal
  expansion property at all. Wired.

## Code defects fixed along the way

The first implementation of `create_low_temperature_alpha` had ten. Four were fatal: a missing
`1 +` on the length ratio (α came out negative and 830× too large, so every log was NaN), the
polynomial fitted to ln C but consumed linearly, and two constructor-ordering faults — the
cryo polynomial built after the spline that samples it, and c_p not yet constructed when the
fit read it. The rest were a sign error in the central second difference, `d2cpdT2_spline`
returning the first derivative, a wrong `d²α/dT²` chain rule, and a guard that was an assert
rather than an error and demanded `d²c_p/dT² > 0` — which is false for copper at 212 K
(−0.00508), so it would have fired on the reference material.

`Bezier::d3ydx3` was also wrong: it implemented the x-by-y chain rule with `dxdXi` where
`d3xdXi3` belongs, and divided by ẏ⁵ instead of ẋ⁵. I additionally "fixed" `dddpoint` to
{−6, 18, −18, 6} and reverted it — the basis is on ξ ∈ [−1,1], so {−0.75, 2.25, −2.25, 0.75}
was correct. A comment now says so.

Ordering is no longer left to convention: `create_low_temperature_alpha` errors if c_p is
absent. The subtlety is that c_p must carry its *final routing* — Copper and Silver flip it to
custom in `create_debye()`, Lead leaves it spline-routed — so fitting against one and
evaluating against another would be silently inconsistent.

## Chromium: an exclusion imposed, then lifted

Chromium was initially excluded from the Grüneisen coupling because its Néel transition
(311 K) and spin flip (123 K) both sit in range and would corrupt a single-γ model. Once its
c_p arrived the numbers contradicted the reasoning: `dln(C)/dT` runs smoothly and monotonically
through both transitions with no feature at either, so the guard failure was the Bézier
artifact, not magnetism. And the exclusion was choosing the larger error — 12.4× at 50 K
against a few-percent anomaly.

Lifted, with the anchor at 273.15 K, above the sign change at 261.8 K and below the Néel point,
inside one coherent antiferromagnetic phase. The margin is 11.4 K, the thinnest in the roster.
γ reads 0.91 against a literature 1.3–1.5 — possibly real, since magnetostriction partly
cancels lattice expansion in Cr, but it is the one cross check that does not corroborate.

## Honest limitations

The construction fixes α, not integrated cool-down strain. For copper the total 293→4 K
contraction goes from 0.009 % too much to about 0.014 % too little; the previous model was
accidentally better on the integral because its α error compensated in the opposite direction.
Anyone reading this as a thermal-strain improvement will be disappointed.

α = C·c_p assumes one Grüneisen channel. Where c_p carries a magnetic contribution — Iron,
Nickel, Chromium — that contribution is scaled by the lattice γ, which is not correct.

## State

Nine materials carry the branch: Copper, Silver, Lead, WhiteTin, Indium, Nickel, Iron,
Aluminum, Chromium. Magnesia does not — its own c_p machinery, and no derived `debye0K`
because MgO has `q ≠ 1`.

Aluminum and Chromium are new. Constants, c_p and α are fitted; `create_mech()` is still a
stub in both, so neither constructs. Neither is registered in `cl_MaterialFactory` by request.
Iron's existing factory registration was left alone — it serves the BH path.

`gTAlphaSwitchMax`, `mTAlphaSwitch`, `alpha_switch_temperature()` and `return_zero` (which
replaced three byte-identical clones) are new on `Material`.

**Nothing here has been compiled or run.** All of it is static review. The
`create_alpha`/`create_alpha_table` split moved calls across constructor boundaries in eight
files and is the part most worth a build check.

## Open

- Elastic data: Lead's frozen E(T); the ν quadratics that make K non-monotonic for copper
  (140.6 GPa at 0 K, 137.9 at 150 K, 155.0 at 800 K, where K must fall); `create_mech()` for
  Aluminum and Chromium.
- Whether K and G should become the primary stored pair with E and ν derived. Argued for on
  three independent grounds: ν is ill-conditioned near 0.5 (`dK/K = 2dν/(1−2ν)`, so a 0.01
  error in ν costs 6 % in K at room temperature and 33 % near melting), `Alloy::create_tables`
  already converts to K and G immediately and back at the end, and the current E/ν pairs
  produce unphysical K(T).
- Magnesia.
