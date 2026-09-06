# Cryogenic α for the Non-Metals, and the Callaway Optical-Channel Fix

**Date:** 2026-08-25 (overnight)
**Topic:** The two defects the documentation overhaul surfaced, fixed with executable gates:
(1) the Callaway kernel's dead optical phonon–electron channel and the YBCO refit that followed;
(2) the metals' old low-temperature thermal-expansion defect in Magnesia, HastelloyC276 and YBCO.
Plus two by-catch defects.
**Participants:** Christian (brief), Claude (Fable); Codex + Grok jury audit of the code diff
(`tmp/ai_exchange/review_alpha_callaway_fixes.md`)

## 1. Callaway kernel and the YBCO refit

`debye.f90` formed `U = hbar·omega_opt/(kB T)` on an `omega_opt` that is already an energy, so the
optical channel evaluated to zero. Fixed (`U = omega_opt/(kB T)`).

The refit was run outside the build: the kernel compiled standalone with gfortran and driven from
Python through `ctypes`, with YBCO's temperature-dependent inputs (θ, v_g, ρ, G, γ, ρ_i, ρ₀)
dumped by the class itself via the `mattest` scratch driver, against Christian's Sommerfeld et al.
2003 table (69 points, 3.8–131 K). Harness and data:
`scratchpad …/ybco/{fit_callaway.py, sommerfeld2003.dat, ybco_inputs.csv}`.

**Result: the fit is degenerate.** The electronic term the class adds, λ_e = L₀T/(ρ_i(T)+ρ₀) —
a normal-state Wiedemann–Franz estimate carried below T_c — reproduces the measured curve on its own
to 10.1 % rms, peak included (7.5/19.5/35.9/43.9/32.9/25.9/17.4/14.9 W/m·K against
8.6/25/36.4/40.8/34.5/27.7/16.5/12.8). Every fit — four parameters, wide bounds, ε free, Γ free —
converges to that same floor by suppressing λ_ph to ≲ 1 W/(m·K), with parameters at their bounds.
Two pre-existing choices cause it: YBCO's resistivity anchor (60 µΩ·cm at 273 K, marked "[fitted]")
is ~4× below literature, which inflates λ_e; and the old Umklapp strength b = 8.17 had been
calibrated against a Grüneisen input of 828 at 2 K — the linear-α artifact of §2 — so with the sane
γ ≈ 2.1 the phonon channel is wide open. The committed parameters with the *old* kernel had the same
10.1 %: the phonon channel never carried information in this model.

**Interim set in the class:** b = 26.14, d = 5.0, Δ₀/(k_BT_c) = 2.1 and λ_opt = 1.06 unchanged
(literature-like), and an *effective* Γ = 1.0 replacing the isotopic 1.1×10⁻⁵ (a coated conductor's
point-defect scattering). Gate: the rebuilt class gives 10.0 % rms, max 22.7 % at 10 K. Documented
as interim in `YBCO::lambda_custom` and in the Callaway guide §8.

**For Christian:** making the phonon term identifiable needs a physical electronic term first — a
literature ρ_ab(T) anchor and a superconducting-state λ_e — which is a modelling decision.

## 2. Thermal expansion of Magnesia, HastelloyC276, YBCO

All three had the metals' old defect (α linear in T near 0 K: Magnesia from a flat-start ΔL/L Bézier,
the other two from α polynomials). Scans of dln(α/c_p)/dT showed it *negative at every temperature*
for Magnesia and YBCO — the curves are too flat relative to c_p everywhere, so the metals' branch
(which matches the curve's derivatives at the split) cannot be attached at any split.

Design chosen: an **anchored** branch. Only the value of α at the split is taken from the curve
(data-backed at room temperature); the derivatives come from the Grüneisen relation itself,
C(T) = C*·K*/K(T) with K from the material's own E and ν, so dln C/dT = −K′/K by construction and
α ≈ C·c_p below the split. `SplineLookupTable::create_cryo_expansion_anchored`; the fit core was
split out of the Bézier wrapper (`create_low_temperature_alpha( alpha, α', α'', poly, aCheckGuard )`),
the split is settable (`set_alpha_switch_temperature`, 0 K = "no branch yet"), and θ_D(0 K) is now
derived with the atom count q (MgO 939 K, YBCO 456 K). The guard on dln C/dT is skipped on the
anchored path — there is no fitted curve to reject (Magnesia's −K′/K is −5×10⁻⁵).

Ordering: Hastelloy's `create_cp()` reads α, E, ν through `cp_from_debye`, so it builds the plain
polynomial α first, then ρ, E, ν, μ, c_p, and attaches the branch last; YBCO attaches it after
`set_young/set_poisson` and before its λ spline, whose Callaway input `grueneisen( T )` reads α.

Gate (debug build, α in 10⁻⁶/K at 4 / 20 / 50 / 77 / 100 / 293 K):

| | before | after |
|---|---|---|
| Magnesia | 1.26 / 3.59 / 5.37 / 6.31 / 6.93 / 10.32 | 0.0001 / 0.011 / 0.24 / 1.03 / 2.23 / 10.32 |
| Hastelloy | 0.45 / 2.26 / 5.6 / 8.6 / 11.2 / 17.0 | 0.055 / 0.32 / 3.29 / 6.96 / 9.76 / 17.0 |
| YBCO | 0.31 / 1.51 / 3.58 / 5.1 / 6.51 / 12.0 | 0.0016 / 0.10 / 2.10 / 4.12 / 5.59 / 11.95 |

YBCO's `grueneisen( 2 K )` went from 828 to 2.08 as a consequence. Hastelloy's plateau of 17×10⁻⁶
is the pre-existing polynomial's value, high against handbook ~12 — flagged, not changed.

## 3. By-catch

- **Hastelloy did not construct in debug builds since `acd4bbf0`**: that commit moved `create_cp()`
  ahead of `create_alpha()`, and `cp_from_debye` asserts on the missing α. Three reviewers missed it in
  the `acd4bbf0` round (the Cu/Ag/Pb dependence was caught, Hastelloy's not). Order restored.
- **`HastelloyC276::create_mu()` registered `MaterialProperty::nu` instead of `mu`** (Codex, from the
  language sweep of all places): the susceptibility was never routed, μ stayed μ₀. My first
  reaction — route it — was **refuted by Grok in the jury round (P0)**: `mu_custom()` evaluates the
  low-T constant on its high-T branch (χ = 0.347 at 300 K, μ_r = 1.347 on the substrate),
  `set_custom( mu )` leaves `dmudH` on the constant path (a debug assert once μ is not constant),
  and the Maxwell calculator would switch the material to its non-constant-μ path. Reverted to
  μ = μ₀ with a comment listing the three things that must be settled before χ(T) goes live.
  **Decision for Christian: is Hastelloy's χ(T) meant to reach the solver?**

## 4. Jury audit and status

Codex and Grok audited the code diff (`tmp/ai_exchange/review_alpha_callaway_fixes.md`). Both
confirmed the kernel fix, the q-aware Debye derivation, the anchored 2-jet algebra (Grok re-derived
it) and the constructor ordering. Applied from the audit: the Hastelloy μ revert above; always-active
preconditions in the anchored branch (K, c_p, α finite and positive, split above the finite-difference
step) and in the fit core (split set); comments that had overclaimed "C(T) = C\*K\*/K(T) below the
split" corrected to what is imposed (the 2-jet at the split; ln C is the fitted polynomial below);
the guard skip documented as load-bearing for Magnesia; `test_callaway` aligned to `lambda_custom`'s
L₀; a scratch `Vector` → `Cell`. Not taken: Codex's preference for keeping `mSplines` private behind a
friend (Grok: protected is the right replacement for a member the subclass reads) — recorded.

Both defects fixed, gates run in a debug build of `mattest` before and after the audit fixes.
`make check` (Christian, morning): 13/14 passed; the one failure was
`YBCOThermalConductivity.RepresentativeTemperaturesAreFinitePositive`, a golden-value test pinning
λ(20/77/100/300 K) to 1e-9 against the old model — the new curve sits 0.4–1.8 % away (the small
phonon share on the dominant electronic term). Baseline re-baked in `tests/physics/test_YBCO.cpp`
with the previous values kept in a comment; the physics test passes again. No commit. `mattest.cpp` now carries the YBCO dump, the α probe and
the split scan as modes — probes, not tests; a `check-fast` test locking α(20 K) of every material to
within a factor two of its c_p-derived value would have caught all of tonight's α defects at once.
