# Thermal Expansion from Heat Capacity {#physics_materials_thermal_expansion_from_heat_capacity}

**Date:** 2026-08-25 (first version 2026-08-20)
**Purpose:** How BELFEM represents the linear thermal expansion coefficient α(T) and the
specific heat c_p(T), why the two are coupled through the Grüneisen equation of state below a
split temperature, and what to do when adding or refitting a material.
**Module:** `src/physics/materials`

> Every number in this document was checked against the curves in the tree. Reference values
> for α and c_p come from literature and carry their own uncertainty: typically a few percent,
> and more below 50 K. When a conclusion relies on a ratio rather than an absolute value, the
> text says so explicitly.

---

## 1. The problem this solves

Thermal expansion data is usually reported as ΔL/L, an integral quantity measured relative to
a fixed reference temperature. BELFEM uses 293.15 K as that reference, but the solver uses the
**tangent** coefficient

\f[\alpha(T) = \frac{1}{L}\frac{dL}{dT}\f]

`cl_Material.hpp` defines this explicitly as `(1/l)·∂l/∂T`, *not* `(1/l)·Δl/ΔT`. The mean
(secant) coefficient from a reference temperature is a different function. For copper at 20 K,
the two differ by about a factor of 40, so the distinction is not cosmetic.

Fitting a curve to ΔL/L and then differentiating it fits one quantity but asks the result for
another. At low temperature the data does not constrain the derivative at all. For copper,
between 4 K and 20 K, the true expansion is about 0.00014 %, while the scatter between
contributed datasets over the same interval is about 0.0099 %. The signal is roughly two
orders of magnitude below the noise, so a least-squares fit can do almost anything there and
still look excellent.

It did. A cubic Bézier through ΔL/L with its first two control points at equal height is
linear in T near the origin, because the parametric map is locally linear and the curve is
locally quadratic. Physics requires α = a₁T + a₃T³ with the lattice term dominant above a few
kelvin. The resulting error, measured against reference values:

| material | α(20 K) as fitted | reference | factor |
|---|---|---|---|
| Copper | 3.16e-6 | 0.30e-6 | 10.5× |
| Aluminum | 4.46e-6 | 0.24e-6 | 18.6× |
| Silver | 7.16e-6 | 1.2e-6 | 6.0× |
| Chromium (at 50 K) | 6.10e-6 | — | 12.4× vs c_p |

Refitting a cubic Bézier cannot fix this. With x locally linear in the parameter and y at best
cubic, the steepest attainable behavior is α ∝ T². The α ∝ T³ behavior is unreachable for
any placement of the control points. The representation has to change, not just the fit.

## 2. The physics

Integrating the Grüneisen relation gives the equation-of-state form used here: thermal
expansion is proportional to thermal energy, so

\f[\alpha(T) = C(T)\,c_p(T), \qquad C = \frac{\gamma\,\rho_{\rm ref}}{3K}\f]

The reference density belongs in C, not the current density. Using ρ(T) makes α implicitly
self-referential, because ρ is computed from the integral of α (`Material::density_custom`,
`cl_Material.hpp:1484`). Using ρ_ref removes that self-reference; the remaining residual is
second order in strain.

This representation is used because three useful properties follow directly:

- **α ∝ c_p.** The shape is inherited from a well-determined curve instead of being fitted
  directly. c_p is measured directly, not recovered as the derivative of an integral.
- **Third-law behavior is structural.** c_p(0) = 0 gives α(0) = 0, and c_p′(0) = γ_Sommerfeld
  gives the correct *non-zero* electronic slope. Imposing dα/dT(0) = 0 would be wrong.
- **The Debye crossover comes along.** It needs no separate low-temperature branch, knots, or
  matching.

C varies by less than 1 % across 0–200 K, because γ and K each vary by less than 1 % over that
range.

## 3. How c_p is represented

`Metal::create_cp` (`cl_Material_Metal.cpp:79`) builds a five-segment curve, or six if a
material passes a third control net:

| range | form |
|---|---|
| T < 0.02·θ_D | `mCpPolys(0) = { beta, 0., gamma, 0. }`, i.e. βT³ + γT (`cl_Material_Metal.cpp:183`) |
| 0.02·θ_D … T₁ | fifth-order beam polynomial, C² in log-log |
| T₁ … T₂ | `mCpBezierLow`, cubic Bézier in ln(c_p) over ln(T) |
| T₂ … T₃ | `mCpBezierMedium`, same |
| T₃ … T₄ | `mCpBezierHigh`, same — optional (`Rx`, `Ry` arguments); Iron and Chromium use it to carry the curve close to the melting point |
| above the last net | linear extrapolation in T |

Log-log coordinates fit this problem well. c_p spans four decades between 4 K and 300 K, but
only about eight units in log space. The T³ region is a straight line of slope 3, and the
electronic term is a line of slope 1, so the asymptotics become linear constraints instead of
curvature conditions at a singular point. The construction in §4 also centers on the log-log
slope, which the Bézier provides through its analytic tangent.

θ_D is **derived from β** whenever β, the molar mass and the atom count q per formula unit are
set: β = 12π⁴ q R/(5 θ_D³), so `debye0K` is recomputed whenever any of M, R, β or q is written.
Until 2026-08-25 the derivation ran only for q = 1; it now serves the compounds as well (MgO,
q = 2: 939 K; YBCO, q = 13: 456 K), which is what lets them use the cryogenic branch.

## 4. The construction

This is implemented once in `SplineLookupTable::create_low_temperature_alpha` (`cl_Material_SplineLookupTable.cpp:161,170`)
and used by all nine pure-metal classes. `HastelloyC276`, `YBCO` and `Magnesia` keep their own α(T) curve above the split and enter
the same construction through the anchored variant described in §7b.

**Split temperature.** `mTAlphaSwitch = min( 0.618·θ_D, gTAlphaSwitchMax )`
(`cl_Material_SplineLookupTable.cpp:253`), with the cap at 273.15 K (`cl_Material.hpp:47`). The Debye factor
keeps the anchor above the steep part of c_p; the cap keeps it inside the range where expansion
measurements carry signal. The 0.618 is a convention, not a derivation.

**Matching.** Define h = ln C = ln α − ln c_p. At T\*, take α, dα/dT and d²α/dT² from the
expansion Bézier, and take c_p and its two derivatives from §3. A cubic in T with no linear
term,

\f[p(T) = aT^3 + bT^2 + d, \qquad p'(0) = 0\f]

is then determined by matching h, h′ and h″ at T\*. Below T\*, α = exp(p(T))·c_p(T)
(`Material::alpha_custom`, `cl_Material.hpp:1639`, routed at `cl_Material.cpp:447`); above it, the Bézier is used unchanged. The
absence of a linear term is not an assumption about α. It follows from C = γρ_ref/3K, where
both γ and K are even in T near the origin.

**Two guards.**

1. `dln(C)/dT > 0` at T\*. Failing this check triggers `BELFEM_ERROR`. Where C decreases with
   temperature, the fitted expansion curve has collapsed relative to c_p, and extrapolating it
   downward inflates α.
2. A turning-point test. The cubic's second stationary point sits at
   T_e/T\* = (h″T\* − 2h′)/(h″T\* − h′). It lands inside (0, T\*) exactly when h″T\* > 2h′. In
   that case the construction falls back to a quadratic that matches only h and h′. Both
   branches are in use: Copper, Iron and Nickel take the cubic, Silver, Lead, WhiteTin, Indium,
   Aluminum and Chromium the quadratic (table in §6).

**Boundary condition.** The α spline's left-hand derivative is
exp(p(0))·dc_p/dT(0) = C(0)·γ_Sommerfeld, not the Bézier's slope. For copper the two differ by
a factor of 367.

## 5. The Grüneisen diagnostic

`create_low_temperature_alpha` also computes

\f[\gamma = \frac{3\alpha K}{\rho_{\rm ref}\,c_v}, \qquad c_v = c_p - \frac{9\alpha^2 T K}{\rho_{\rm ref}}\f]

at the split temperature. The model does not use it. The diagnostic is an independent
cross-check on four separately fitted quantities. It reads O(2) for every metal, and it depends
only on the anchor, not on the extrapolation below it. Below T\*, α = C·c_p makes γ nearly
constant by construction, so a flat γ there is not confirmation of anything.

The K it uses comes from the same E(T) and ν(T) the solver sees: since 2026-08-24 these are
the Wachtman modulus and the Grüneisen-derived Poisson ratio of `Metal::create_mech`
(`cl_Material_Metal.cpp:1105`), fitted against Blanke 1989
for most metals. It is checked rather than enforced because the elastic data of several
materials is still provisional. The assertion uses a deliberately wide 0.2–10 band: that is a scale check, not a
precision check. It is not printed during construction — the natural place to surface it is the
material report alongside density, molar mass and the Sommerfeld coefficient.

This diagnostic has already found four defects that nothing else surfaced:

- **Lead** — the former `mYoungPoly` had no working temperature dependence; its T³, T² and T
  terms summed to about 45 kPa across 0–600 K against a 19.888 GPa constant, so E was frozen at
  its 0 K value. Retired on 2026-08-24, when every metal's E(T) became a Wachtman curve.
- **Indium** — expansion curve 27 % low, in both α(293 K) and γ. Refitted 2026-08-19.
- **Nickel** — `ref_density` was 890 kg/m³ rather than 8900, giving γ = 20.76. Fixed.
- **Iron** — `create_alpha()` was defined and declared but never called, so the material had no
  expansion property at all. Wired 2026-08-19.

## 6. Results

| material | θ_D [K] | T\* [K] | branch | γ(T\*) | literature γ |
|---|---|---|---|---|---|
| Lead | 104.1 | 64.3 | quadratic | — | 2.7–2.8 |
| Indium | 108.8 | 67.2 | quadratic | inconclusive | ~2.4 |
| WhiteTin | 199.2 | 123.1 | quadratic | 2.06 | 2.1–2.3 |
| Silver | 226.5 | 140.0 | quadratic | 2.30 | 2.3–2.5 |
| Copper | 343.8 | 212.5 | cubic | 1.98 | 1.96–2.00 |
| Aluminum | 417.9 | 258.2 | quadratic | 2.19 | 2.1–2.2 |
| Nickel | 456.0 | 273.15 | cubic | 1.73 | ~1.9 |
| Iron | 470.0 | 273.15 | cubic | 1.75 | ~1.7 |
| Chromium | 592.7 | 273.15 | quadratic | 0.91 | 1.3–1.5 |

Lead's and Indium's γ values are not diagnostic: both have ν near 0.45, so `1 − 2ν ≈ 0.1` and
K is dominated by the uncertainty in ν. Chromium's low value may be physical — magnetostriction
partly cancels lattice expansion there — but it is the one entry that does not corroborate.

Copper after correction, against reference values:

| T [K] | 20 | 50 | 77 | 100 | 150 |
|---|---|---|---|---|---|
| model [1e-6/K] | 0.288 | 3.844 | 7.668 | 10.050 | 13.085 |
| reference | 0.30 | 3.80 | 7.90 | 10.30 | 13.50 |

## 7. Adding or refitting a material

1. **Fit ΔL/L, not α.** Use 293.15 K as the reference and make the value exactly zero there.
   Enforce that by eliminating one control point analytically rather than by penalizing it.
2. **Keep the low-temperature data in the fit.** Those ΔL/L values are reliable even where the
   derivative is noise, and they stop the curve's left end from swinging. Excluding data below
   200 K made Indium's α *worse*, not better, despite improving the residual.
3. **Keep the control polygon monotone in y.** Unconstrained fits reach slightly lower
   residuals by letting α go negative in gaps where nothing constrains it. For Indium this cost
   0.0002 % of residual and removed an α of −10.3e-6 near 10 K.
4. **Set β after M** so `debye0K` is derived, and check the derived value against literature.
5. **Fit c_p** and check `0.02·θ_D < T₁`, which `Metal::create_cp` asserts.
6. **Check the split.** T\* must clear the `dln(C)/dT` sign change. Measured, that sits at
   0.43–0.51·θ_D where it exists at all, so 0.618·θ_D clears it by construction — but the flat
   cap breaks that relationship once θ_D exceeds about 425 K. Aluminum clears it by 63 K,
   Chromium by 11.4 K.
7. **Read the logged γ.** It should be O(2). An order-of-magnitude miss is a unit or scale
   error, not imprecision.

Call order is fixed: **constants, then c_p, then α — and only then mechanics, ρ and the
Debye curve.** `create_alpha()` builds the expansion Bézier and finishes by calling
`create_cryo_expansion()` (`cl_Material_SplineLookupTable.cpp`), which fits the branch,
sets the boundary condition and builds the spline and its integral. Two dependencies pin the
order:

- **c_p before α**, because the branch is fitted against c_p and its first two derivatives.
  `create_low_temperature_alpha` enforces this with a `BELFEM_ERROR` on `have(cp)`.
- **α before `create_debye()`** was a third constraint until 2026-08-25: the Debye inversion
  from c_p, `Metal::cp_from_debye` (`cl_Material_Metal.cpp:612`),
  carries the dilation term 9α²K T/ρ and therefore needs α, the elastic data and `density`.
  Copper, Silver and Lead now invert through `Metal::compute_debye_from_cv`
  (`cl_Material_Metal.cpp:527`), which drops that term — at their
  12–14 K anchors it is of order 10⁻⁵ of c_p, and the routine refuses above roughly 100 K, where
  it would not be. The elastic moduli are now consumed during construction only by
  `create_mech()` itself, and α only by that routine's validity guard, so `create_debye()` may
  run anywhere after `create_cp()` and `create_alpha()`. A Debye construction that reads no c_p
  at all — Iron's and Nickel's inversion from measured resistivity — is free to run even earlier.

c_p must also carry its **final routing** before the fit, or the branch is fitted against one
curve and evaluated against another. Copper and Silver therefore set c_p custom at the end of
their own `create_cp()` rather than leaving it to `create_debye()`; Lead keeps it
spline-routed throughout, which is equally consistent. `Copper::create_alpha` is the
reference implementation.

## 7b. The non-metals: an anchored branch

Magnesia, HastelloyC276 and YBCO do not fit their expansion data as a dL/L Bézier: Magnesia's
Bézier is too flat below ~200 K, Hastelloy's and YBCO's α(T) are polynomials in T that are
*linear* near 0 K. In all three, α/c_p decreases with temperature everywhere, so the guard of
§4 rejects any split — the curves' derivatives are not trustworthy at any temperature below
room temperature, while their room-temperature values are data. Since 2026-08-25 they use the
**anchored** form of the branch (`SplineLookupTable::create_cryo_expansion_anchored`,
`cl_Material_SplineLookupTable.cpp:175`):
only the *value* of α at the split is taken from the curve; its first two derivatives **at the
split** come from the Grüneisen relation C = α/c_p = C(T\*)·K(T\*)/K(T), with K from the
material's own E and ν, so dln C/dT = −K′/K there. Below the split the construction of §4 then
applies unchanged: ln C is the fitted cubic (or quadratic) with zero slope at 0 K, not 1/K(T).
The guard of §4 is skipped on this path, and that is load-bearing, not cosmetic: Magnesia's K
rises imperceptibly at the split (−K′/K = −5×10⁻⁵), which the strict check would refuse; the
resulting drift of C across the branch is of order 1 %.

Two ordering facts follow. The branch needs c_p and K, so it runs after `create_cp()` and
`create_mech()`; Hastelloy's `create_cp()` in turn reads α, E and ν through `cp_from_debye`, so
its constructor builds the plain polynomial α first (split at 0 K, no branch), then E, ν and
c_p, and attaches the branch last (`HastelloyC276::create_alpha_cryo`). YBCO does the same
after `set_young` / `set_poisson` have scaled its moduli and before its λ spline is sampled —
that spline's Callaway input `grueneisen( T )` reads α, and with the old linear α it evaluated
to 828 at 2 K instead of ≈ 2.1.

Values before and after (α in 10⁻⁶/K):

| material | 4 K | 20 K | 50 K | 77 K | 100 K | 150 K | 293 K |
|---|---|---|---|---|---|---|---|
| Magnesia, before | 1.26 | 3.59 | 5.37 | 6.31 | 6.93 | 8.00 | 10.32 |
| Magnesia, after | 0.0001 | 0.011 | 0.24 | 1.03 | 2.23 | 5.39 | 10.32 |
| Hastelloy, before (polynomial) | 0.45 | 2.26 | 5.6 | 8.6 | 11.2 | 16.7 | 17.0 |
| Hastelloy, after | 0.055 | 0.32 | 3.29 | 6.96 | 9.76 | 13.9 | 17.0 |
| YBCO, before (polynomial) | 0.31 | 1.51 | 3.58 | 5.1 | 6.51 | 8.78 | 12.0 |
| YBCO, after | 0.0016 | 0.10 | 2.10 | 4.12 | 5.59 | 8.03 | 11.95 |

The room-temperature values are unchanged by construction. Hastelloy's plateau value of
17×10⁻⁶ is the polynomial's own and predates this work; it is high against handbook values
near 12×10⁻⁶ and is worth a look at the source data.

The 0 K Debye temperature these three need is derived from β with the atom count q per
formula unit (§3): Magnesia sets M and q = 2, YBCO q = 13, Hastelloy runs with q = 1.

## 8. Pitfalls

| pitfall | symptom | cause |
|---|---|---|
| Fitting α instead of ΔL/L | good residual, wrong derivative | the data constrains the integral, not the slope |
| Matching at the fit's own boundary | α at T\* badly off | a Bézier endpoint has data on one side only, so its slope is nearly free |
| Trusting the Bézier below T\* | α off by 3–19× | it is never evaluated there; only its value at T\* matters |
| Shortening the fit range | worse, not better | fewer constraints let the curve swing |
| Reading γ below T\* as validation | false confidence | α = C·c_p makes γ flat by construction there |
| Reversing the sign of h | α negative or enormous | h = ln α − ln c_p, not the sum |
| Building α before c_p | `BELFEM_ERROR` at construction | the branch is fitted against c_p and its first two derivatives |

## 9. What is not represented

- **Magnetic anomalies.** Neither the log-log Bézier form for c_p nor a smooth expansion Bézier
  can represent a λ peak. Chromium's Néel transition at 311 K and spin flip at 123 K, along
  with the magnetic contributions in Iron and Nickel, are smoothed over. For Chromium the
  anchor is deliberately placed between the two transitions, but the branch still extrapolates
  through the lower one.
- **A magnetic Grüneisen channel.** α = C·c_p uses a single γ. If c_p carries a magnetic
  contribution, that contribution is scaled by the lattice γ. That is not correct and can even
  have the wrong sign.
- **Integrated cool-down strain is not improved.** The construction fixes α, not the integral.
  For copper the total 293→4 K contraction changes from 0.009 % too much to about 0.014 % too
  little. The previous model was accidentally closer for the integral because its α error
  compensated. The residual is set by the quality of the 0 K anchor in the expansion data.
- **Magnesia** does not use this path — it has its own c_p machinery. It does, however, set
  `q = 2` and `beta` (`cl_Material_Magnesia.cpp:56,59`), from which the base setter derives and
  stores `debye0K`; the earlier claim that it has none was wrong.
