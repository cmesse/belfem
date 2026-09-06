# Jury Review of the Wachtman Elasticity Commit + Plumbing Fixes

**Date:** 2026-08-24
**Topic:** Three-AI review of `bcc47bcd` ("add elasticity data for pure metals") — Wachtman E(T)
plus Grüneisen-derived ν(T) for Al, Cr, Cu, In, Fe, Pb, Ni, Ag, Sn — and the plumbing fixes
applied on approval the same night.
**Participants:** Christian, Claude (Fable), Codex, Grok
**Exchange record:** `tmp/ai_exchange/review_wachtman_elastic.md`
**Related:** `dl20260824_iron_kohler_review.md` (earlier the same day; its §5 "elastic data
remains open" is superseded by this entry)

## 1. What the commit adds

`Metal::create_mech( E0, b, T0, T2, nu2 )`: Young's modulus as a Wachtman curve
E = E0 − b·T·exp(−T0/T) (E0, b in GPa), sampled into the E spline; the Grüneisen parameter from
the room-temperature anchor (K_T = E/(3(1−2ν₂)) → K_S → γ = α_V K_S/(ρ c_p)), then held constant
to derive ν(T) = ½ − E/(6K_T(T)) with K_S(T) = γρc_p/α_V on the spline grid, ν(0 K) extrapolated
with zero slope. Nickel carries two Wachtman branches joined by a C¹ Bezier bridge across the
Curie region (the ΔE dip to 124 GPa near 487 K is real — Christian's Blanke table confirms it).

## 2. Review outcome

All three reviewers converged; every citation verified. By source trace **the commit had not
been built**: three compile errors (Silver's leftover zero-arg `create_mech()`, Indium's orphaned
`create_mech()` definition, a public→protected move in `Metal.hpp` that broke `Alloy` and
`mattest`), and every material failed at construction — seven still overrode `E_custom` with
the old `polyval( mYoungPoly )` whose fill routines were deleted (virtual dispatch never reached
the Wachtman code), Nickel never constructed its bridge Bezier, Copper/Silver/Lead called
`create_mech` after `create_debye` although `cp_from_debye` evaluates K(T), and the E spline was
clamped with a NaN slope because `create_spline( E )` defaults its Tangent slope to NaN (Grok's
find — the one defect that also hit Iron, the only metal without an override). Inside the ν
loop: `theta( n )` out-of-bounds write, first grid point at T = 0 where α = 0, `cp` frozen at
room temperature. The thermoelastic algebra itself was correct throughout.

Refuted by data: my concern about the nickel magnitude (Blanke: 183.7 GPa at 293 K, minimum
124 at 487 K, 195.7 at 637 K — the fit reproduces it). Left as Christian's decisions: the
Wachtman curves are nearly flat below room temperature (copper softens 1.4 GPa from 4 K to
293 K where the deleted copper.org polynomial softened 8.6 GPa; Blanke's tables start near
93 K, so the cryo range is the form's own extrapolation); `set_kohler_dependencies` no longer
`set_custom( lambda )`; dynamic (adiabatic) moduli labelled isothermal (sub-percent on E);
two Grüneisen parameters now coexist (constant property vs `Metal::grueneisen( T )`).

Record note: Al/Cr/Ni are being registered in `cl_MaterialFactory.cpp` in Christian's
uncommitted working tree; at HEAD they are not.

## 3. Fixes applied (on approval)

- `Metal::create_mech` rewritten: preconditions as `BELFEM_ERROR`; Wachtman block guarded on
  `have( E )` so Nickel keeps its own E; **E spline built here for everyone with
  `create_spline( E, dEdT_custom( 0 ) )`**; null-spline check; ν loop with `theta( k )`, first
  point at dT, `cp( T )` updated, `K_S → K_T → ν`; 0 K extrapolation kept; ν range guard on
  (−1, 0.5); ν spline installed through the public `set_spline` (the new `friend class Metal`
  removed). Comment "isotropic young's modulus" → adiabatic bulk modulus.
- `Metal.hpp`: access layout restored (the `set_bh_curve … set_table_flags` block is public
  again, `protected:` reinstated before `create_J_spline`).
- Dead elastic overrides removed: `E_custom`/`nu_custom` + `mYoungPoly*`/`mPoissonPoly*`/
  `mTYoungSwitch` in Aluminum, Chromium, Silver, Lead, WhiteTin, Indium, Copper (Copper's
  `E_custom` body and copper.org ν formula, Indium's Kim & Ledbetter polynomials retired to git
  history — both now come from `create_mech`).
- Debye/elastic interdependence removed (Christian's design, superseding an interim reorder):
  the three Debye inversions sit at T₁ = 14 K (Cu), 12 K (Ag), 12 K (Pb), where the dilation
  term (c_p − c_v)/c_p = α_V·γ_G·T is of order 1e-5. `Metal::cp_from_debye` is split into
  `cv_from_debye` + dilation; a new `compute_debye_from_cv` inverts the c_v form through a
  shared `invert_debye` regula falsi (member-function pointer, no duplicated solver) and refuses
  with `BELFEM_ERROR` where 6·α(T)·T > 1e-3 — a guard that needs only α, which every metal has
  before `create_debye`. Cu/Ag/Pb call it; their constructors keep the committed order.
  `compute_debye_from_cp` stays with an always-active E/ν precondition for any future high-T
  caller. Silver's zero-arg `create_mech()` call became the five-arg call, duplicates removed.
- Nickel: `mYoungBezier = new Bezier( x, y )`; `dEdT_custom( 0 ) = 0`.
- Copper/Silver: the leftover `set_custom( E )` in `create_debye()` (a relic of the polynomial E)
  raised `have( E )` before `create_mech` ran and silently skipped the Wachtman assignment —
  Christian caught it after the c_v redesign restored the constructor order. Removed; and
  `Metal::create_mech` now decides on `E0 > 0` rather than on the have-flag (Nickel passes 0
  and must have declared its own E, enforced by `BELFEM_ERROR`).
- Docs: Al/Cr class comments rewritten (registered, elastic data present); Cr comment
  "270 GPa" corrected to the 279 GPa the parameters actually give.

## 4. Status

Reviewed, not verified — no build in this session. Executable gate at the next rebuild:
`make` (the three compile errors are gone by inspection), then construct each of the nine
metals and print E, ν at 4 / 77 / 293 K; ν must sit in (0.2, 0.5) and rise slightly with T,
E must reproduce the Wachtman values (Cu 134.4 / 134.4 / 133.0 GPa). Copper's `create_debye`
now runs on the Wachtman K(T) instead of the old polynomial — the inverted Debye curve will
shift slightly; worth a glance.
