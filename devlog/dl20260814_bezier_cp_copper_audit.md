# Bezier-Based cp for Copper, Silver and Indium: Audit and Fixes

**Date:** 2026-08-14
**Purpose:** Session record — audit of the new `Metal::create_cp()` Bezier construction for the
metal specific heat, the fixes that make it live for Copper, and the `Bezier::xi_by_x` defect that
porting Silver and Indium exposed.
**Module:** `src/physics/materials`, `src/numerics/bezier`
**AIs:** Claude (audit, numeric replication, fixes), Christian (implementation, ruling to apply
all findings)

## What the change does

`Copper::create_cp()` used to build cp from nine polynomials, eight switch temperatures and four
hand-rolled transition polynomials, partly in `T` space and partly in `x = ln(T)` space, with a
`k == 0` / `k < 3` / `else` ladder in the middle to handle the space flip. It is replaced by a
generic `Metal::create_cp( Px, Py, Qx, Qy )` that takes two cubic Bezier curves in log-log space
(`ln cp` over `ln T`), joins them at a shared control point, connects the low end to the
Sommerfeld-Debye cubic `cp = gamma*T + beta*T^3`, and extrapolates linearly above the last
control point.

Copper's whole cp definition is now four vectors of four numbers.

The connection to the Sommerfeld-Debye cubic went through two iterations in this session. The
first was a cubic beam polynomial in `T` space, C1, which needed a `max( exp( Px(0) ), T0 + 5 )`
clamp and a `T1 += 1.0` monotonicity walk to stay well behaved. Christian replaced it with a
fifth-order beam polynomial built entirely in log-log space, matching value, slope **and
curvature** at both ends, which lets the clamp and the walk go away. Because `x = ln T` and
`y = ln cp` is a diffeomorphism, matching `y, y', y''` at `x0` is equivalent to matching
`cp, cp', cp''` at `T0`, so the transition is now C2 on both joints. Measured: the old version
left a 25% curvature jump at `T0` (0.0310 below, 0.0386 above); the new one matches to
finite-difference precision (0.0310 vs 0.0313). The curve itself barely moves, at most
0.004 J/(kg K) over the window, so this buys a smooth `dcp/dT` for the Newton path rather than
a better fit.

The transformations are worth writing down, since they are what the log-log space costs:

```
dy/dx   = T cp' / cp
d2y/dx2 = T^2 cp'' / cp - (dy/dx)^2 + dy/dx
```

Both were checked against finite differences of `ln( cp( exp( x ) ) )` and agree.

## Quality of the fit

The construction was replicated numerically (Bernstein basis, the bisection inversion from
`Bezier::xi_by_x`, the beam polynomial, both segments) and compared against Touloukian/NIST
copper:

| T [K] | 30 | 50 | 100 | 300 | 600 | 1000 |
|---|---|---|---|---|---|---|
| error | +0.08% | -0.55% | -0.63% | -0.01% | +0.33% | -0.29% |

Sub-percent from 30 K to 1000 K, C1-continuous at all four breakpoints, strictly monotone over
the full sampling range. Below 30 K the deviation grows to a few percent, which is where the
Sommerfeld-Debye connection takes over and where the 4 K spline sampling grid is coarse anyway.
The method is sound; every defect below was in the plumbing around it.

## Defects found and fixed

### P0-1: the old override shadowed the new implementation

`cl_Material_Copper.hpp` still declared `cp_custom` plus its own `mXCpSwitch`, `mTCpSwitch` and
`mCpPolys`. Those members shadowed the new ones in `Metal`, so `Metal::create_cp` filled
`Metal::mCpPolys` (3 entries) and `Metal::mTCpSwitch` (4) while virtual dispatch sent every cp
call to `Copper::cp_custom`, which read `Copper::mCpPolys` at length zero. Debug build: bounds
assert. Release: out-of-bounds read.

Copper never declared `dcpdT_custom`, so `Metal::dcpdT_custom` was used and did read the new
data. cp and its derivative would have come from two different piecewise definitions.

Fixed by deleting the three shadowing members, the `cp_custom` declaration and its definition.

### P0-2: beta, gamma and debye0K were never set for Copper

The old `create_cp` produced them (`set_constant( beta, ... )`, `set_constant( gamma, ... )`);
the new one consumes them. Nothing filled the gap, and `Copper::set_constants()` only set
`T_max`, `ref_density`, `T_ref_density` and `M`. All three asserts at the top of
`Metal::create_cp` fired, and `create_cp()` runs before `create_debye()` in the constructor
anyway, so there was no ordering that would have rescued it.

Nickel already showed the intended pattern. Copper now sets both coefficients in
`set_constants()`, after `M`, using the values from the old fit:

```cpp
this->set_constant( MaterialProperty::gamma, 0.010969801407435 );
this->set_constant( MaterialProperty::beta,  0.000753034913427 );
```

The order is load-bearing: `Material::set_constant` derives
`debye0K = (2.4 pi^4 R / beta)^(1/3)` on the `beta` write, and that needs `M` (hence `R`)
already present. This reproduces theta0 = 343.75 K.

### P0-3: missing exp() on the last control point

`real y3 = Qy( 3 )` took the logarithm as if it were the value. The low-side twin one screen up
does `std::exp( Py( 0 ) )` correctly. The linear extrapolation above `T3` therefore returned
6.26 instead of 522 J/(kg K), a factor of 83, and `dydT3` was wrong by the same factor.

Not academic: `SplineLookupTable::create_spline` samples in 4 K steps up to
`ceil( T_max / 4 ) * 4 = 1360 K` while `T3 = 1357.0 K`, so the last spline knot landed in that
branch and wrecked the top of the table. After the fix the knot sequence reads
521.05 / 521.92 / 522.80 across 1352 / 1356 / 1360 K.

### P1: an assertion that failed on correct data

`Py( 3 ) == Qy( 3 )` should be `Py( 3 ) == Qy( 0 )`: the two curves share the *last* P and the
*first* Q. For Copper, `Py(3) = Qy(0) = 5.627793` while `Qy(3) = 6.257940`, so the consistency
check aborted the debug build on data that was in fact consistent. The check is an exact `==`
on doubles, which is fine while both lists are literals copied between the two argument
vectors, and brittle the day they are computed.

### P2: robustness

Three copies of `mCpBezierLow == nullptr`, two of which were meant to be `mCpBezierHigh`, so the
high curve was never checked and one line was a duplicate. Added asserts that the four bases
have four control points each (a short vector would silently resize `mX` while `mWork` stayed at
4) and that `Px` and `Qx` are monotonic, which the bisection in `xi_by_x` assumes without saying
so.

Two further items existed only in the first (cubic, C1) iteration and were retired with it: the
beam-poly anchor was read at `xi = -1` while `T1` could be clamped to `T0 + 5`, and the
`T1 += 1.0` monotonicity walk had no upper bound. The quintic version needs neither.

### Defects in the C2 rewrite

The rewrite introduced four of its own, all in the same seam — the boundary between log-log
space and linear `T` space, which is exactly where the original defects lived too:

- `mCpPolys( 2 ) = { dydx3, y3 - dydx3 * Tx }` referenced an undefined `Tx`. Compile error.
- `dydx3 = mCpBezierHigh->dydx( x3 )` is `d(ln cp)/d(ln T)`, but `mCpPolys( 2 )` is evaluated by
  `cp_custom` as `polyval( ..., T )`, i.e. linear in `T`. The slope needs the same carry-back as
  the value: `dcp/dT = cp * dy/dx / T`. For copper the raw log-log slope is 0.567 while the true
  `dcp/dT3` is 0.218, so the extrapolation above melting was steeper by a factor 2.6. The value
  at `T3` was right either way, which is why it would not have been obvious.
- `dcpdT_custom`'s last branch still read `std::exp( dpolyval( mCpPolys( 2 ), T ) ) / T`, the
  log-space pattern applied to a polynomial that is already linear in `T`. Now a plain
  `dpolyval`.
- Dropping the clamp means `T0 < T1` is no longer guaranteed by construction. If a metal's Debye
  temperature puts `theta/50` above the first control point, the interval inverts, the quintic
  Vandermonde degenerates and the `T < T0` branch swallows the range the Bezier was meant to
  cover. Now a `BELFEM_ERROR` naming both temperatures. Copper is comfortable: 6.88 K vs
  12.09 K.

## Porting Silver and Indium: the bracket in `xi_by_x`

Silver and Indium were ported to the new `Metal::create_cp` next. Silver came over clean; Indium
still carried the Copper pattern in its header (`mXCpSwitch`, `mTCpSwitch`, `mCpPolys` plus a
`Cell< Bezier * > mCpBeziers` that nothing fills any more, and a destructor loop over it), all
now removed, along with a stray `#` null directive between two method declarations.

The interesting part is what the two new fits exposed. Both were badly wrong on first evaluation:
Silver had an 18% discontinuity at `T3` and a *negative* extrapolation slope above it, Indium was
+119% at the bottom of its Bezier range. The cause is neither the fits nor the new cp code.

**`Bezier::xi_by_x` bisects on `[-1.25, 1.25]`** — it extrapolates the curve by 12.5% past each
end before searching. A Bezier only interpolates its control polygon for `-1 <= xi <= 1`; outside
that, if the outer control-point spans are lopsided, the cubic folds back on itself. The
bisection then sees **two** roots of `x( xi ) = X` inside its bracket and can converge on the one
from the extrapolated branch, silently, with no iteration-count warning.

Auditing every Bezier basis in the tree:

| curve | outer spans | monotone on [-1,1] | monotone on [-1.25,1.25] | inversion at the endpoints |
|---|---|---|---|---|
| Silver cp high | 0.766 / 1.880 / **0.038** | yes | **no** | last point returns `xi = +1.25` |
| Indium cp low | **0.009** / 0.308 / 0.136 | yes | **no** | first point returns `xi = +1.25` |
| Iron reduced magnetization | 0.999 / 0.001 / 0.000 | yes | **no** | last point returns `xi = +1.008` |
| Lead / Magnesia / WhiteTin alpha | 50 / 243 / 307 | yes | **no** | endpoints happen to be safe |
| the other twelve | — | yes | yes | ok |

**Every curve in the tree is monotone on `[-1, 1]`; five are not on `[-1.25, 1.25]`.** So the
bracket is the whole defect, and the fix is to search the domain the curve actually has. Queries
outside it now saturate at the end rather than folding over. Copper never triggered it because
all six of its spans are balanced, which is why the first audit pass came back clean.

### Why in-domain queries are not protected by the piecewise switch

Christian's (correct) objection: `cp_custom` switches curves at the breakpoints, so the Bezier
is only ever *queried* inside its own `[X0, X3]`. The answer is that the defect sits in the
solver bracket, not in the query value. Quantified on the Silver high curve
(spans 0.766 / 1.880 / **0.038**):

- The outer control span sets the parametric speed at the end: `dx/dxi( +1 ) = 1.5 * 0.038 =
  0.058`, while the acceleration there is `1.5 * ( X1 - 2 X2 + X3 ) = -2.76`. The extrapolated
  cubic therefore stops almost immediately: x rises only until `xi = 1.020` (overshoot 6e-4),
  then plunges to `x( 1.25 ) = X3 - 0.078`.
- The descending branch re-crosses every `x` in `[X3 - 0.078, X3]` — in temperature,
  **every T between 1132 K and 1223 K, all of them in-domain**, acquires a second root inside
  the `[-1.25, 1.25]` bracket.
- For those queries the bracket is not even a bracket: at `x = X3`, `f(-1.25) = -2.91` and
  `f(+1.25) = -0.078` — same sign. The bisection's invariant is void, it wanders, and it
  terminates *successfully* (`|f| < 1e-12`) on the extrapolated branch.

The second, sharper hit point: **`Metal::create_cp` itself queries the exact domain endpoints**
— `y( x1 )`, `dydx( x1 )`, `d2ydx2( x1 )` anchor the quintic at the left end of the low curve,
and `dydx( x3 )` anchors the extrapolation at the right end of the high curve. Endpoint queries
are maximally fold-vulnerable, which is how Indium's entire quintic window got poisoned (wrong
anchor values, +119% at 10 K) and Silver's `dcpdT3` came out negative, even before any runtime
cp evaluation.

Why `[-1, 1]` is provably safe: the derivative of a cubic Bezier is a Bernstein-weighted
combination of the three control spans, and Bernstein weights are non-negative on the domain —
so increasing control abscissae give a strictly monotone `x( xi )` on `[-1, 1]`, hence a unique
root. That is exactly the precondition the new `Px`/`Qx` monotonicity asserts pin down. The
bracket fix is therefore exact, not a workaround; the saturation branch covers the endpoint
queries and the few spline knots past `T_max`.

Blast radius, since this function also serves thermal expansion, Kohler, Debye and the Iron
magnetization curve: for every in-domain query the root is unchanged and the bisection converges
in fewer steps. Only out-of-domain queries change behaviour, and the sole ones in the tree are
the thermal-expansion curves at the very top of their spline grid, where `create_spline` rounds
`T_max` up to a multiple of 4 K. Copper queries 2.0 K past the end of its alpha curve, Silver
0.92 K; the resulting change in alpha is 1.3e-4 relative, at a knot above the melting point.
Everything else lands inside.

With that fixed, all three materials are continuous at every breakpoint, monotone both densely
and on the spline knots, and free of overshoot in the quintic window:

| | theta [K] | T0 / T1 / T2 / T3 [K] | cp(T1) vs `gamma*T + beta*T^3` | extrapolation above T3 |
|---|---|---|---|---|
| Copper | 343.8 | 6.9 / 12.1 / 113.7 / 1357.0 | 1.047 | 1.0 K (0.1%) |
| Silver | 226.5 | 4.5 / 16.1 / 83.5 / 1223.2 | 1.200 | 11.9 K (1.0%) |
| Indium | 108.8 | 2.2 / 13.3 / 20.9 / 272.7 | 1.038 | 157.0 K (36.5%) |

One observation on the fit data rather than the code: Indium's last control point is at 272.7 K
while its `T_max` is 429.7 K, so 36.5% of the temperature range is served by the linear
extrapolation. It extrapolates to 241.6 J/(kg K), which is physically sensible, but pushing the
control point to `ln( 429.75 ) = 6.063` would put the range under the curve.

An earlier flag on Silver — its first control point sitting 20% above `gamma*T + beta*T^3` at
`T1` — is withdrawn as a concern: the ratio `cp( T1 ) / cubic` measures each metal's low-T
Debye-temperature dip blended with the Debye-function flattening, and the observed ordering
(Cu 1.05 minimal dip, In 1.04 dip cancelled by flattening at `T1 = theta/8`, Ag 1.20 mild dip,
Pb 1.84, Sn 2.45 strong soft-phonon dips) is exactly the physically expected one. The quintic
absorbing that gap is the design working, not a fit wobble; see the Lead/White Tin section.

## Lead and White Tin ports

Both came over with clean headers on the cp members except White Tin, which still carried the
shadow `mCpPolys`/`mTCpSwitch` pair and a stale `mTCpSwitch = { 1, 20 }` write ahead of the
`Metal::create_cp` call (dead but confusing, same disease as Copper/Indium; removed). Lead's
`create_debye` calls `this->cp_custom( 12.0 )`, which now dispatches to `Metal::cp_custom`; the
constructor runs `create_cp` first and 12 K lands mid-domain on the low Bezier, so the handoff
is sound.

Numerically both pass everything Copper/Silver/Indium pass: C1 at all four joints (relative
jumps <= 2e-6), densely monotone, and near-zero extrapolation share (Lead 7.6 K = 1.3% of
range, slope +0.066; White Tin 0.08 K = 0.0%). The Qy control polygons of both high curves are
non-monotone (Lead dips -0.021, White Tin -0.32 in log space), but the curves themselves never
dip: min `d(ln cp)/d(ln T)` is +0.071 (Pb) and +0.126 (Sn) — the convex-hull property smooths
the polygon dip into a plateau, so cp stays strictly increasing.

Physics checks, all consistent:

- **Constants match literature.** Pb: `gamma = 1.4718e-2 J/(kg K^2)` is 3.05 mJ/(mol K^2)
  (literature 2.98-3.1), derived `theta0 = 104.1 K` (literature 105). Sn: `gamma` and `beta`
  taken directly from O'Neil 1964, derived `theta0 = 199.2 K` (literature ~200).
- **The theta_D dips are absorbed where they belong.** The quintic's log-log slope peaks at
  3.64 (Pb) and 4.45 (Sn) inside the transition window — a super-Debye rise (slope > 3) is the
  signature of a dipping Debye temperature, and Pb/Sn are the strong soft-phonon cases. Spot
  values track estimate-grade references: Pb cp(5 K) = 1.70 vs ~1.7 for `theta_eff ~ 88-90 K`,
  cp(10 K) = 14.6, cp(300 K) = 128.2 vs 128.5 (CRC), cp(600 K) = 143.9 (ratio to Dulong-Petit
  120.4 = 1.19, plausible with dilation + electronic terms). Sn cp(4 K) = 0.193 inside O'Neil's
  measured range, cp(300 K) = 228.1 vs the CRC 228.
- **White Tin mid-range (50-200 K) sits 5-9% above my reference estimates** — but those
  estimates are low-to-medium confidence Debye-function arithmetic, not tables, and the
  Touloukian curve set the fit was made against is the authority. Worth one glance at the
  dataset, not flagged as a defect.
- **Superconductivity of the metals is out of scope by design** (Christian's ruling, this
  session): the database serves HTS modeling, where Cu, Ag, Pb, Sn and In enter as stabilizer,
  solder and matrix materials. Their own LTS transitions (Pb Tc = 7.19 K, Sn 3.72 K,
  In 3.41 K) are deliberately not modeled — normal-state fits are the intended model, not an
  approximation to be tightened. YBCO is the only class carrying jc, n and the E-J power law.
  (Magnitude check, kept for the record: even below Tc the lattice dominates the electronic
  term 9:1 already at 4 K for Pb, so the distinction would be invisible in enthalpy anyway.)

No physics defects found in either port.

### Pre-existing, off this path: `Bezier::xi_by_y` never converged

Two bugs in the same bisection loop: the midpoint was `0.5 * ( tXi1 - tXi0 )` instead of
`0.5 * ( tXi0 + tXi1 )`, and the bracket update wrote `tF = tF0` instead of `tF0 = tF`. That
breaks `x(y)`, `dxdy` and `d2xdy2`. Replicating both variants on the copper thermal-expansion
curve: the old loop fails to converge for every sampled point, so it was a guaranteed
"Too many iterations" abort rather than a silent inaccuracy, which is why nothing had noticed.
The fixed loop recovers xi exactly. The cp path only uses `xi_by_x`, which was correct.

## Evidence posture

The cp construction is **verified numerically** in replication, in both iterations: the six
quintic constraints are reproduced to 5e-13, value and slope are continuous at all four
breakpoints, the analytic `dcpdT_custom` agrees with finite differences of `cp_custom` to 1e-8
relative across every branch including the two above `T3`, the curve is monotone over the whole
sampling range, and the reference comparison is the table above. The C++ itself is **reviewed,
not compiled** — see below.

One junction is C1 and not C2, by construction rather than by defect: the two Bezier curves meet
at `T2` sharing a control point and a tangent direction (the copper data is G1 there, both sides
give slope 0.674899), but nothing constrains their curvatures. That is a property of joining two
independent cubics, not something the connection code could fix.

## Iron and Nickel ports: the ferromagnets

Both use the same four-vector call; both are structurally clean (no shadow members, connection
points exact, `T0 < T1`, C1 at all interior joints, densely monotone to `T_max`). The
ferromagnet-specific design point: **both fits stop below their Curie peaks by construction** —
Nickel's `T_max = 600 K` sits under Tc = 631 K, Iron's `T_max = 860 K` (a pre-existing
resistivity-validity cap) under Tc = 1043 K — and both have `T3 > T_max`, so the linear
extrapolation branch is unreachable and every spline knot lands inside the high Bezier. The
critical approach to Tc is visibly captured: `dcp/dT` at the top is 0.78 (Ni, 600 K) and 0.72
(Fe, 860 K) against ~0.45 at room temperature, both still rising.

Spot checks (estimate-grade references): Ni cp(300) = 442.6 vs 444 CRC, cp(500) = 530.1,
cp(600) = 595 approaching the lambda peak; low-T consistency with `gamma*T + beta*T^3` at 10 K
to +0.2%. Fe cp(300) = 449.4 vs 449 CRC, cp(600) = 573.0 vs JANAF 574, cp(800) = 686.4 vs ~690.
Derived theta0: Ni 456.0 K (literature 450-477), Fe 470 exactly (set explicitly). Sommerfeld
gammas in molar units: Ni 7.46, Fe ~5.0 mJ/(mol K^2) — both literature-consistent, Ni's large
band-ferromagnet gamma included.

Two fixes applied in Iron:

- **beta/gamma were double-set.** `set_constants()` writes the new Touloukian pair
  (beta = 3.35249e-4, gamma = 0.0893559), then `create_debye()` — which runs before
  `create_cp()` — overwrote both with its older values (beta re-derived from theta = 470,
  numerically identical to 5 digits; gamma = 9.0896e-2, 1.7% off). The overwrite shifted the
  quintic anchor zone by 1.3-1.6% below 12 K relative to the values the fit was made with.
  `create_debye` now sets only `debye0K = 470` and defers beta/gamma to `set_constants`.
- **Leftover probe removed:** `mTDebyeSwitch.print( "T" )` fired on every Iron construction
  (probe-removal policy).

Two report-only observations:

- **Nickel's `Tcurie = 672 K` refuted and corrected to 631 K** — a three-stage finding that
  reversed twice, worth recording in full. Usage traced on request: exactly one consumer chain
  — `Nickel::compute_mred` scales the reduced-magnetization Bezier by it, feeding only
  `Ferromagnetic::rho_mag`'s blend `m^2 A_magnon T^2 + (1 - m^2)^2 A_spin_disorder`, hence
  `rho_custom` and `compute_debye_from_rho`; nothing else (cp, mu, BH) touches it, and the
  blend is sensitive to the choice (at 600 K the weights nearly swap: 0.38/0.39 vs 0.22/0.61).
  Claude first identified 672 as the paramagnetic Weiss temperature, and a
  `MaterialProperty::Tweiss` split was built on that reading (Iron checked in the same pass:
  its 1043 K is the genuine Curie point in both roles, unchanged throughout). Christian then
  challenged the identification (external sources quote ~654 K), and the check against the
  code's own data settled it: **the reduced-magnetization Bezier is calibrated to the true
  Tc**, tracking the canonical Crangle-Goodman shape on `t = T/Tc` and matching measured
  m(600 K) ~ 0.45-0.50 when scaled by 631 (model 0.47) while 672 gives 0.61 — so 672
  de-calibrated the very fit it scaled. Reported Weiss temperatures for Ni scatter over
  650-675 K (fit-window dependent), so the identification was never pinnable; and Nickel's
  transport is not wired yet (`set_rho_i_ref` commented out), so 672 could not have been a
  resistivity calibration either. Final state: `Tweiss` removed again, enum numbering
  restored (the renumber had been audited safe — no integer casts, no serialization, no
  string map — with the user-material .so ABI as the one rebuild caveat, now moot),
  `Tcurie = 631` with the calibration evidence in the comment, `compute_mred` back on
  `Tcurie`. Lesson kept: the earlier "load-bearing, do not correct without refitting" warning
  protected a calibration that never existed. By-catch, report-only: the
  `A_spin_disorder = 3e-8` line's own comment quotes a typical range of 80-120 nOhm m, i.e.
  8e-8 to 1.2e-7 — the set value is below its own quoted range.
- `Iron::create_alpha()` exists (Bezier to 1185 K) but the constructor never calls it, so
  `mThermalExpansion` stays null and alpha is simply not enabled for Iron — pre-existing, and
  harmless because `alpha_custom` is only routed via `set_custom( alpha )` inside the uncalled
  function. Nickel's port wired its `create_alpha` into the ctor; Iron's looks like the same
  step not yet taken.

Gate: `cl_Material_Iron.cpp`, `cl_Material_Nickel.cpp`, `cl_Material_Ferromagnetic.cpp` pass
`g++ -fsyntax-only` with the build tree's flags; the cp constructions verified numerically in
replication as above.

## Pre-commit review: Fe/Ni resistivity and thermal conductivity

Christian added transport to the ferromagnets: `create_debye_and_rho()` in both (a Bloch-
Grueneisen anchoring that subtracts `rho_mag` from the measured total, solves for the BG
amplitude at 295 K and re-anchors `rho_i_ref` at 273.15 K), a new even-polynomial
`debye_custom` for Nickel, and Hust-form lambda coefficient sets for both. Kohler and mech
deferred. Reviewed numerically before commit:

**Iron verifies end to end.** theta(295) = 400.0 as documented; the anchoring reproduces the
total rho(295) = 9.8e-8 exactly, and the derived `rho_i_ref(273.15) = 7.68384e-8` matches the
old commented-out calibration value 7.68395e-8 to five digits — the new machinery lands on the
old number independently. Lambda through `Metal::lambda_custom` with the NBSIR 84-3007
coefficients: 129 (100 K), 83.2 (300 K), 50.4 (600 K) W/(m K) — all inside reference bands;
the large low-T values are the ideal-crystal envelope, correct for `rho_0 = 0` semantics.

**Nickel: structure verifies, the lambda coefficient set does not.** theta(295) = 390.00 (the
White & Woods target), theta positive over the whole range (min 209 at 631 K); anchoring
reproduces rho(295) = 7.0e-8 exactly, `rho_i_ref = 5.089e-8`. But evaluating the committed
lambda chain (validated on Iron through the identical code path) gives **48.5 W/(m K) at
300 K against the handbook 90.7** — low by 1.5-2.6x over 100-600 K (100 K: 62.7 vs ~165;
600 K: 42.7 vs ~66). The high-T plateau of the fit is analytically `T^0.181 / 136.4`, which
cannot reach the Touloukian values the comment cites, so this is the wrong coefficient set
(possibly an earlier fit iteration), not an evaluation subtlety. Flagged for refit;
coefficients not touched.

**Two structural defects fixed:**

- **`rho_0` was never set for either metal → NaN transport.** The old Iron was safe by
  omission (its `create_rho` never called `set_custom( rho )`, so `Ferromagnetic::rho_custom`
  was dead); the new wiring makes it live, and both `rho_custom` and `Metal::lambda_custom`
  read `rho_0` first (`beta = rho_0 / L0`). Default-constructed Fe/Ni returned NaN for rho
  AND lambda. Both now set `rho_0 = 0` (ideal-crystal residual, the Indium/Lead/WhiteTin
  pattern; `set_RRR` overrides for real samples). Benign for both lambda forms: Fe has
  `p0 = p7 = 0`, and for Ni `pow( 0, -1.46 ) = inf` puts `C = p7/inf = 0`.
- **Nickel called `set_RRR` before `create_debye_and_rho`**, i.e. before `rho_i_ref` exists.
  The old "will cause an error" comment is wrong on both counts: transport exists now, and
  with NaN `rho_i_ref` the regula falsi exits silently on its first NaN comparison, leaving
  `rho_0 = NaN` with no diagnostic. The RRR block moved to the end of the ctor (Iron already
  had it there); stale comments removed in both files.

**Report-only, conscious-decision items:** Nickel's `T_max` moved 600 → 631 K, which makes
the cp linear tail reachable for the spline knots at 604-632 K (cp 598 → 620, smoothing over
the lambda-anomaly peak the Bezier never fitted) and pushes the alpha spline past its Bezier's
600 K endpoint (saturation clamp returns the endpoint value — alpha flat over 600-632 K).
Both defensible, both worth knowing. And `Metal::set_RRR` defines RRR from
`( rho_i + rho_0 ) / rho_0` only — for a ferromagnet the experimental RRR includes
`rho_mag( 273 )` (for Ni ~1.3e-8 against rho_i ~5.1e-8, a ~25% definitional difference), so
stated RRR values from measurements are not directly this parameter.

Gate: Iron, Nickel, Ferromagnetic TUs pass `g++ -fsyntax-only` with build flags; all transport
numbers above verified in replication (BG integral J_4.5 by quadrature).

## Compile blocker resolved: six missing powerlaw declarations

The blocker was larger than the one symbol the first error showed. A mechanical diff of every
out-of-line `Material::` definition in `powerlaws.hpp` against the declarations in
`cl_Material.hpp` (matching name AND arity) found six missing, all from the DR-07 T-leg work:
`djc_eval_dT`, `dn_eval_dT`, and both overloads (4-arg and 8-arg defect variant) of
`drho_powerlaw_dT` and `drho_piecewise_dT`. The dB siblings of all six were declared; only the
T-leg twins were not — the definitions landed with their design comment but the header half of
the interface was never written. All six declared now, with doc comments distilled from the
DR-07 design block (parallel-factor and floor conventions, the no-early-out rule for the rho_n
term, the branch-by-branch contract of the piecewise variant).

**Gate:** all ten materials TUs (`cl_Material`, `Metal`, `Copper`, `Silver`, `Indium`, `Lead`,
`WhiteTin`, `YBCO`, `HastelloyC276`, `Nickel`) plus `cl_Bezier.cpp` now pass `g++ -fsyntax-only`
with the build tree's exact `flags.make` flags. Syntax-verified, not yet built or run.

## Shadow-member sweep, concluded

Final census of `mCpPolys` / `mTCpSwitch` / `mXCpSwitch` / `mCpBezier*` across the materials
headers: the five ported metals (Copper, Silver, Indium, Lead, WhiteTin) are clean; `Metal`
holds the one real set; `Magnesia` keeps its own legitimately (it derives from
`SplineLookupTable`, not `Metal` — no shadowing). **`YBCO` and `HastelloyC276` still shadow**
`Metal::mCpPolys`/`mTCpSwitch`, but theirs are NOT removable: both keep their own `cp_custom`
(YBCO's cp carries the superconducting transition, Hastelloy is an alloy — neither fits the
two-Bezier scheme) and those overrides consume the members.

That shadowing exposed a latent regression the sweep was worth doing for: `Metal` now overrides
`dcpdT_custom`, and YBCO/Hastelloy inherit it WITHOUT overriding it — with Metal's cp members
empty and the Bezier pointers null for them. The path is not live today (`dcpdT_custom` is only
routed via `set_custom( cp )`, which of the Metal family only Copper and Silver call; YBCO and
Hastelloy run cp through splines), but it was one `set_custom( cp )` away from a bounds assert.
`Metal::dcpdT_custom` now null-guards on `mCpBezierLow` and falls back to
`Material::dcpdT_custom` (the finite difference), restoring the base-class contract for
subclasses that manage cp themselves.

## Follow-ups

- Build and run the Copper/Silver/Indium/Lead/WhiteTin paths (`make` + a cp table sanity look).
- The other metals still call the old cp scheme; Silver in particular still carries the nine
  polynomial ladder that Copper just shed, and is the obvious next port.
- `Metal::create_cp` reads `debye0K` only to place `T0 = theta * 0.02`. Worth a comment on why
  theta/50 is the validity bound of the Sommerfeld-Debye cubic, since the constant is otherwise
  unexplained.
