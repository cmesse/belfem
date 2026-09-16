# Quasi-Harmonic Elastic Moduli: ν(T) and E(T) from the Thermal Strain

**Date:** 2026-09-15 (rev. 10: ceilings set from Blanke's digitized curves; O7 lead level pending)
**Purpose:** Replace the one-anchor constant-γ closure in `Metal::create_mech` — which makes copper's
Poisson ratio fall with temperature (0.363 → 0.344) where measurement rises (0.3401 → 0.3471,
Ledbetter 1981, 10.1002/pssa.2210660209) — by a closure in which the **isothermal** bulk modulus and the
shear modulus follow the logarithmic volumetric thermal strain ε(T) = 3∫₀ᵀ α dT with one softening
constant each, K_T = K₀·exp(−δ_K ε), G = G₀·exp(−δ_G ε); E and ν follow from K_T and G. The moduli are
then consistent with the fitted expansion curve by construction and with the measured elastic data
through δ_K, δ_G. c_p, ρ and α enter once more, in R1, through the exact identity that converts the
measured adiabatic bulk moduli to isothermal ones before the fit; the served curves are pure exponentials.
**Module:** `src/physics/materials` (`cl_Material_Metal`, nine pure-metal classes, tests, docs)
**AIs involved:** Claude (exploration + plan), Codex (audit: terra/high on rev. 1, terra/xhigh on rev. 5),
Grok (audit: 4.6/high on rev. 1, 4.6/xhigh on rev. 5)
**Status:** ✅ COMPLETE 2026-09-15, pending Christian's commit. Landed: `Metal::create_mech( E2, nu2, T2, deltaK, deltaG )`
with quasi-harmonic isothermal K and G on the logarithmic thermal strain, nine constructors with fitted constants and
DOIs, Nickel's ΔE machinery removed with a documented rationale, Chromium constant ν, Lead at the static level (crystal K, Blanke E; handbook
values reproduced; O7 closed 2026-09-15), validity ceilings from the shape comparison with
Blanke (Cu 1000, Ag 900, Cr 570, Sn 400 K), the γ diagnostic printed at 298.15 K, `tests/physics/test_MetalElastic.cpp`
(4 tests), docs rewritten and swept. Verified: `test_physics` 14/14 and `make check` 17/17 on `cmake-build-debug`
(debug + MKL) after the final rebuild; three plan-jury rounds and one code jury, every finding verified and fixed.
No residual. Per Christian, no follow-up todo is spun out; this file moves to `todo/closed/` with the commit.

> **Scope guards:**
> - IN: `Metal::create_mech`, `Metal::E_custom`, a **new** `Metal::nu_custom` (the base
>   `Material::nu_custom` is a `BELFEM_ERROR` stub, `cl_Material.cpp:532-535`; ν is today a hand-built
>   vector + `set_spline`, `cl_Material_Metal.cpp:1182-1229`), the nine call sites (Al, Cr, Cu, Ag, In, Pb,
>   Sn, Fe, Ni), a regression test, the module docs and class header comments.
> - OUT: `YBCO::create_mech`, `Magnesia::create_mech`, `HastelloyC276::create_mech` (own constructions,
>   `cl_Material_YBCO.cpp:145`, `cl_Material_Magnesia.cpp:90`, `cl_Material_HastelloyC276.cpp:90`).
>   Magnesia builds mech **before** α because its anchored cryogenic α branch needs K
>   (`cl_Material_Magnesia.cpp:32-34`) — a K-from-ε closure would be circular there; never "unify" it onto
>   this model. `UserDefined` (user `E_custom`, `cl_Material_UserDefined.hpp:332-335`). `Alloy` needs no
>   edit: it forms K, G from each component's `E( T )`, `nu( T )` per grid point
>   (`cl_Material_Alloy.cpp:330-341`) and mixes those, and its Schapery α uses K (`:391-393`); its E, ν
>   **and α** will move silently, so R3 probes one alloy and R4 smoke-tests one.
> - Kept: the public accessors `E( T )`, `nu( T )`, `K( T )`, `G( T )` and their units; `K` and `G` stay
>   the isotropic identities on the served E and ν (`cl_Material.hpp:2228-2238`); the `material` report
>   layout; the `grueneisen` constant property the report prints (`main.cpp:357-359`).
> - Dropped: the Wachtman parameters (E₀, b, T₀) for all nine metals, including Nickel's two-branch
>   Wachtman + Bézier E(T) with the ΔE dip (`Nickel::create_young`, `E_custom`, `dEdT_custom`,
>   `mYoungData`, `mYoungBezier`; ruling O3). `Metal::dEdT_custom` goes too: nothing calls it
>   (grep 2026-09-15: only the declarations and the two overrides); the E spline's 0 K slope is the literal
>   `0.0` passed at `cl_Material_Metal.cpp:1150`.
> - **User-visible changes to call out (served isothermal values):** copper at room temperature moves
>   from E = 133.0 GPa, ν = 0.344 (Wachtman fit, `cl_Material_Copper.cpp:64`) to E ≈ 127.8 GPa,
>   ν ≈ 0.343 (Ledbetter 1981, 10.1002/pssa.2210660209, Table 2, converted to isothermal — the measured
>   adiabatic values are 128.17 GPa and 0.3471; the all-measurement average is 123.5 ± 0.7 GPa,
>   Ledbetter & Naimon 1974, 10.1063/1.3253150, Table 10). Nickel moves from E = 184.0 GPa, ν = 0.31
>   (demagnetized Blanke fit) to ≈ 223 GPa, ≈ 0.297 adiabatic at 300 K (Alers, Neighbours & Sato 1960,
>   10.1016/0022-3697(60)90125-6, Table 2, Hill), slightly lower served — **+21 %**. Every other metal
>   gets the same before/after column in R1 (§6).

---

## 1. Current Behaviour and How It Fails

`Metal::create_mech( E0, b, T1, T2, nu2 )` (`cl_Material_Metal.cpp:1122-1230`):

1. E(T) = E₀ − b·T·exp(−T₀/T) (Wachtman), `cl_Material_Metal.hpp:845-850`, sampled into a 4 K
   spline (`cl_Material_SplineLookupTable.cpp:131-150`) with the literal 0 K slope `0.0` (`:1150`).
2. At the anchor T₂: K_T = E/(3(1−2ν₂)), K_S from the thermoelastic correction, γ = α_V K_S/(ρ c_p)
   (`:1169-1172`). γ is stored as the `grueneisen` constant and **held constant in T**.
3. On the grid: K_S(T) = γ ρ c_p / α_V, K_T from K_S, ν = ½ − E/(6 K_T) (`:1197-1200`); 0 K point
   extrapolated with zero slope (`:1205-1210`); range guard (−1, ½) (`:1213-1218`).

| Failure | Mechanism | Evidence |
|---|---|---|
| ν(T) has the wrong sign of slope for copper | Step 3 assigns the entire T-variation of α_V/c_p to K. Printed 3α/c_p rises 12 % from 10 K to 293 K → model K falls 12 % (implied K_T 163.6 → 142.1 GPa); measured adiabatic B falls 3.3 % (144.46 → 139.74 GPa, Ledbetter 1981, 10.1002/pssa.2210660209, Table 2), isothermal ≈ 6.1 % (probe, §2) | probe on the `material copper` table, 2026-09-15; `devlog/dl20260915_copper_poisson_trend.md` |
| The sign is set at the split temperature | The cryogenic α branch requires dln(α/c_p)/dT > 0 **at T\*** (one-point check, `cl_Material_SplineLookupTable.cpp:335-338`); with constant γ, K_S ∝ c_p/α_V must then fall through T\* | source trace |
| Copper's E(T) nearly flat below room temperature | T₀ = 922 K (`cl_Material_Copper.cpp:64`) suppresses softening below ~300 K: −1.1 % from 0 to 293 K vs −7.5 % measured (138.62 → 128.17 GPa, Ledbetter 1981 Table 2). A property of that fit, not of the form: Indium's T₀ = 44.6 K (`cl_Material_Indium.cpp:32`) is active in the cryogenic range. Per-metal ν slopes on the current tree are **not yet probed** (R1) | `devlog/dl20260824_wachtman_elastic_review.md:36-39` (its "Blanke starts near 93 K" is unverified against Blanke) |
| Copper's room-temperature E is high | Wachtman E(293) = 133.0 GPa vs 128.2 measured / 123.5 ± 0.7 average | Ledbetter 1981; Ledbetter & Naimon 1974 Table 10 |
| Same closure for every metal | all nine constructors call the routine | `cl_Material_{Aluminum:46,Chromium:61,Copper:64,Silver:82,Indium:32,Lead:51,WhiteTin:47,Iron:41,Nickel:39}.cpp` |
| Adiabatic anchor treated as isothermal | ν₂ is a dynamic literature value; `:1169` labels the derived K isothermal; `grueneisen( T )` (`:684-703`) and the `cp_from_debye` dilation term (`:619-620`) treat the served K as K_T. The usage guide says the tabulated moduli are adiabatic and the effect on E is 0.3–0.5 % (`materials_usage_guide.md:596-597`); on **K** it is 3.0 % for Cu at 295 K, hence 0.0041 on ν — comparable to ν's whole 5–295 K change (0.0070) | O2; probe §2 |

**Bottom line:** inverting the Grüneisen relation for K with a constant γ copies the (unknown,
≈ 10 %) temperature variation of γ one-to-one into a modulus that in reality varies by a few percent;
combined with a flat E(T) that inverts dν/dT. The closure has to change, not its constants.

## 2. Architecture: Why the Thermal Strain Is the Right Spine

**Quasi-harmonic picture.** An elastic modulus depends on volume, not on temperature directly, and
(∂ln M/∂ln V)_P = −δ_M is a slowly varying constant. For the **isothermal** bulk modulus this is the
Anderson–Grüneisen parameter δ_T = −(1/(α_V B_T))(∂B_T/∂T)_P, and integrating at 1 bar gives
B_T(T) = B_T(0)·exp(−∫₀ᵀ α_V δ_T dT) — Garai & Laugier 2007, Eq. (8)–(10) (J. Appl. Phys. 101, 023514,
10.1063/1.2424535; arXiv:physics/0601101), who note δ_T is not constant in general. With constant δ and
l(T) = exp(∫_{T_ref}^T α dT) (`Metal::l`, `cl_Material_Metal.cpp:368-375`; integral installed by
`finish_cryo_expansion`, `cl_Material_SplineLookupTable.cpp:232-243`):

    ε(T) ≡ 3 ∫₀ᵀ α dT = 3 ln( l(T) / l(0) ) = ln( V/V₀ )          (logarithmic volumetric strain, O4)
    K_T(T) = K₀ · exp( −δ_K ε(T) ),   G(T) = G₀ · exp( −δ_G ε(T) )
    E = 9 K_T G / ( 3 K_T + G ),      ν = ( 3 K_T − 2 G ) / ( 6 K_T + 2 G )

**The law is written for K_T and the fit is made on K_T (rev. 6).** Every dataset in O1 is ultrasonic, i.e.
adiabatic; the shear modulus is the same in both descriptions (a pure shear has no volume change), the
bulk modulus is not. R1 therefore converts every measured K_S datum to K_T with the exact identity
K_T = K_S / (1 + α_V² T K_S / (ρ c_p)), using BELFEM's own α(T), ρ(T), c_p(T) at the datum's temperature,
**before** the regression. At 0 K the two coincide, so K₀ is state-independent; only δ_K changes. The
served curves are then the pure exponentials above — no serve-time conversion, no 0/0 at T = 0, no
c_p feature leaking into E or ν, and the citation to the isothermal law is honest. The previous
revision's serve-time conversion was the source of both round-2 P0s (§8).

**Structural properties, each exact for the served curves:**

- consistency with the expansion data by construction — the moduli inherit the fitted α curve and its
  Debye crossover; no T₀ to fit. Below the lowest elastic datum (4–5 K) the curves follow the cryogenic
  α branch, a **controlled extrapolation** (α = C·c_p with the Sommerfeld slope), not free-hand;
- dE/dT = dν/dT = 0 at 0 K, because dε/dT = 3α → 0 — **provided the cryogenic α branch is installed**
  (R2 asserts `alpha_switch_temperature() > 0`; a provisional α spline with a Bézier slope at 0 K would
  break the claim);
- ν(T) is non-decreasing wherever α ≥ 0 **iff δ_G > δ_K**: d ln(G/K_T)/dT = −(δ_G − δ_K)·3α. With
  δ_K > 0 both moduli soften monotonically. These are theorems about the served curves, so the
  construction-time check in 3.1 runs over the **whole** served range and is a logic guard;
- ν ∈ (−1, ½) automatically while K_T, G > 0 (the exponential keeps them positive; ν → ½ is approached,
  never reached, for Indium ν₂ ≈ 0.45 and Lead ν₂ ≈ 0.44 — the sampled-grid range guard stays as a
  belt-and-braces check on spline overshoot);
- γ(T) = α_V K_S/(ρ c_p) = α_V K_T/(ρ c_v) stays a **diagnostic** (`Metal::grueneisen( T )` exists,
  `cl_Material_Metal.cpp:684-703`, and already treats the served K as K_T) instead of a closure — a
  wrong δ shows up as an unphysical γ(T), the role γ already plays for α
  (`thermal_expansion_from_heat_capacity.md` §5).

**Pilot fit, copper (2026-09-15, probe).** Ledbetter 1981 Table 2 (polycrystalline Cu, 5–295 K, 16
points; 10.1002/pssa.2210660209) against ε from the trapezoid integral of BELFEM's printed α(T), with the
adiabatic B converted to K_T using BELFEM's printed α, c_p and ρ = ρ_ref/l³:

| modulus | δ | M₀ [GPa] | rms | max dev |
|---|---|---|---|---|
| K_T (isothermal, **the served one**) | 6.52 | 144.26 | 0.10 % | 0.16 % |
| K_S (adiabatic, for the record) | 3.46 | 144.42 | 0.04 % | 0.10 % |
| G (identical in both) | 8.67 | 51.53 | 0.23 % | 0.36 % |

ε(295 K) = 0.947 %. K_S/K_T − 1 = 3.00 % at 295 K, 0.58 % at 100 K, 0 at 5 K. Served copper then reads
E_T = 127.78 GPa and ν_T = 0.3430 at 295 K, 138.62 GPa and 0.3401 at 5 K (measured adiabatic: 128.17 /
0.3471 and 138.62 / 0.3401). δ_G > δ_K holds in both descriptions. The δ values must be **refit in R1
against the code's own ε (spline integral, 4 K grid)** — the trapezoid on the coarse report table is only
accurate to a few percent in ε. Ledbetter's own Varshni fits (Table 4) give B(0) = 144.49, G(0) = 51.72,
E(0) = 138.61 GPa, ν(0) = 0.3401 — cross-check for K₀, G₀.

**Rejected alternatives.**
- *γ(T) model inside the present closure* — γ is the poorly known quantity; the inversion
  K = γρc_p/α_V amplifies its error into K one-to-one.
- *Tabulate ν(T) directly per metal* — full-range data exist for few metals (Christian, 2026-09-15); a
  table carries no structural guarantee of dν/dT(0) = 0 or monotonicity. Measured ν(T) is instead the
  **fitting target** for δ_K, δ_G.
- *Varshni form C(T) = C(0) − s/(exp(t/T) − 1)* (Ledbetter 1981, Eq. 4) — fits copper's five constants
  well, but its parameters carry no link to α(T), so ν would again be free of the expansion data.
- *Keep Wachtman E, add a K(T) model* — two unrelated representations for one thermodynamic object;
  the mismatch is exactly what produced the present sign error.
- *Fit the adiabatic K_S and convert at serve time* (rev. 3–5) — puts c_p(T) into every served modulus
  (Nickel's magnetic c_p rise toward T_C would re-appear as a wiggle in E, ν), divides 0/0 at T = 0 where
  `create_spline` samples, and misattributes Garai's isothermal law to an adiabatic fit. Rejected in
  round 2 (Codex P0, Grok "spec holes").

## 3. Gap Table

| # | State | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| 1 | ε(T) = 3 ln(l(T)/l(0)) | K_T(T), G(T) | `l( T )` exists, referenced to T_ref = 293.15 K | (a) from `l( T )` and `l( 0 )` cached once at construction | `cl_Material_Metal.cpp:368-375` |
| 2 | anchor (T₂, E₂, ν₂) — **isothermal**, from R1 — → K₂, G₂ → K₀, G₀ | model constants | adiabatic E₀ from Wachtman today | (a) K₂ = E₂/(3(1−2ν₂)), G₂ = E₂/(2(1+ν₂)), M₀ = M₂·exp(δ_M ε(T₂)); E₂ enters in GPa and is converted to Pa **at the `create_mech` entry**, as today (`:1138`); served E(T₂) = E₂ exactly | isotropic identities, `cl_Material.hpp:2228-2238` |
| 3 | δ_K, δ_G per metal | the T-dependence | absent | (c) fitted offline (R1) on the code's ε(T): ln K_T and ln G vs ε, K_T from the measured K_S by the exact identity | O1 data; pilot above |
| 4 | E(T), ν(T) served to callers | all consumers | 4 K splines (`create_spline`; `set_spline` routes `E_spline`/`nu_spline`, `cl_Material_SplineLookupTable.cpp:75-83`) | (a) keep the spline plumbing; sample the new `E_custom`, `nu_custom`; **`set_have( nu )` before `create_spline( nu, … )`** — it routes to `nu_custom` (`cl_Material.cpp:1083-1087`) and debug-asserts `have( nu )` (`:1040`) | today's ν bypasses `create_spline` |
| 5 | 0 K boundary slope of the E and ν splines | correct low-T shape | literal `0.0` for E (`:1150`); hand-built ν with a Tangent BC | (a) `create_spline( E, 0.0 )` and `create_spline( nu, 0.0 )` — both slopes are exactly 0 (α(0) = 0). No `dEdT_custom`: nothing calls it; `Metal`'s override is removed with the Wachtman members. End BC becomes Parabolic (`cl_Material.cpp:1045-1046`) where today's ν uses NoCurvature (`:1221-1229`): affects the last 4 K interval only; accepted | `create_spline` signature `:1036-1046` |
| 6 | Nickel: E(T) with the ΔE dip across T_C = 631 K | ν for Ni | `Nickel::E_custom` (two Wachtman branches + Bézier bridge from 465.4 K, `cl_Material_Nickel.hpp:103-135`, `cl_Material_Nickel.cpp:295-341`), `create_mech( 0., 0., 0., 293.15, 0.31 )` anchored to a **demagnetized** E(293) = 184.0 GPa (Blanke fit) | (c) **RESOLVED O3 → removed.** Nickel is treated like the other metals, fitted to the **saturated-state** constants of Alers, Neighbours & Sato 1960 (10 kOe, 10.1016/0022-3697(60)90125-6, Table 2, 20 K steps 0–760 K); the header documents why (§6a) | the dip is the ΔE effect of the demagnetized state, absent at saturation (Ledbetter & Reed 1973 §15, 10.1063/1.3253127); the coded curve is already −18 % at 400 K |
| 7 | Chromium: Néel point 311 K, spin-flip 123 K | E, ν for Cr | the class smooths every anomaly by design (`cl_Material_Chromium.hpp:28-72`) | (a) smooth curve, **no** exception; R1 anchors Cr **away** from 311 K (the header's ±20 K halo) — e.g. the 250 K row of Palmer & Lee — and fits on the windows outside both anomalies | O3 |
| 8 | Iron near T_C | — | T_max = 860 K < T_C = 1043 K (`cl_Material_Iron.cpp:198-201`) | (a) never served | — |
| 9 | above the fitted window to T_max — **RESOLVED 2026-09-15 by ceilings**: T_max lowered where the quasi-harmonic E departs > 5 % in shape from Blanke's digitized curve (`tmp/blanke_wachtman/`): Cu 1000, Ag 900, Cr 570, Sn 400 K; Fe, Al unchanged | high-T E for thermal-stress runs | Wachtman fit against Blanke (293 K … T_max) | (c) exp(−δ ε) continues from the α Bézier; R3 tabulates the deviation from the retired Wachtman curve per metal **into §6**; no gate, by design. Near T_m the real G collapses (Nadal & Le Poac 2003, 10.1063/1.1539913, Eq. 5–6); the quasi-harmonic form does not — accepted, mechanics near melting is out of BELFEM's use | Cu at 1000 K: measured E ≈ 0.6 E₀; the model gives ≈ exp(−8·0.05) ≈ 0.67 E₀ |
| 10 | α range vs T_max (ε above the Bézier's last knot runs on a clamped α, `Bezier::xi_by_x`, comment `cl_Bezier.hpp:104`, clamp `cl_Bezier.cpp:97-98`) | ε(T) near T_max | Nickel: α Bézier to 600 K, T_max 631 (`cl_Material_Nickel.cpp:121, 85`); **Chromium: to 1700 K, T_max 2180** (`cl_Material_Chromium.cpp:132, 87`, documented in its header) | (a) inherited; R1 prints the α range vs T_max for all nine; R5 documents Ni and Cr | Grok finding |
| 11 | adiabatic vs isothermal | labelling; `cp_from_debye` dilation term (`:619-620`); `grueneisen( T )` (`:684-703`) | mislabelled; both consumers treat the served K as K_T | (a) **O2 RESOLVED → isothermal fit (rev. 6)**: served K is K_T by construction; both consumers are **correct as written**; ρ in the R1 conversion is BELFEM's ρ(T) = ρ_ref/l³, c_p and α the served curves | Δν = 0.0041 at 295 K for Cu (§2 probe) |
| 12 | `grueneisen` constant for the report | `main.cpp:357-359` | set from the anchor at `:1173` | (a) γ = α_V K_S/(ρ c_p) at the anchor with K_S = K_T(T₂)·(1 + α_V² T₂ K_T/(ρ c_p)) computed **directly from the constants**, not via `grueneisen( T₂ )` (which calls `K( T )` → `nu( T )` before the ν spline exists) | Codex finding, round 1 |
| 13 | `compute_debye_from_cp` precondition on E and ν | Debye curve above ~30 K | `BELFEM_ERROR( have( E ) && have( nu ) )` `:530-533`; no metal constructor calls it today | (a) unchanged: `create_mech` still sets both flags | constructor order unchanged |
| 14 | Wachtman members `mWachtmanYoung`, `Metal::E_custom`, `Metal::dEdT_custom` | — | live (`cl_Material_Metal.hpp:93, 839-866`) | (c) replaced by named scalars `mK0`, `mG0`, `mDeltaK`, `mDeltaG`, `mL0`; `E_custom` body rewritten; `dEdT_custom` override deleted (uncalled); no subclass override points remain | `Vector` is for linear algebra (`doc/coding_philosophy.md`) |
| 15 | regression test | gate | none touches E or ν (`tests/physics/CMakeLists.txt:14-19`); the physics directory sets no `TESTLABELS`, so `check-fast` skips it (`Add_Test.cmake:43-45`) | (c) R4 with numeric oracles, run under `make check` | — |
| 16 | docs and comments | user-facing | `materials_usage_guide.md:119, 254, 271-272, 585-597`; `material_property_sources.md:25-33`; `thermal_expansion_from_heat_capacity.md:148-150` (cites `:1105`, stale); class headers `cl_Material_Aluminum.hpp:33-35`, `cl_Material_Chromium.hpp:30-32`; `doc/literature_references.md` has none of the elastic sources | (c) **R5**, Codex language sweep | Christian 2026-09-15: state that the roster is built from **polycrystalline** (isotropic) data |
| 17 | Alloy | inherits | mixes component K, G (`cl_Material_Alloy.cpp:336-341`), Schapery α uses K (`:391-393`) | (a) no edit; R3 probes one formula alloy before/after, R4 smoke-tests one | numbers move silently otherwise |

### 3.1 Cross-cutting findings

- **Sign of dν/dT is a property of two numbers.** R1 prints δ_G − δ_K per metal. `create_mech` refuses
  (BELFEM_ERROR, once-per-run) δ_K ≤ 0 and δ_G ≤ δ_K, and additionally checks on the sampled 4 K grid
  over the **whole served range [0, T_max]**: α ≥ 0, ν non-decreasing, ν ∈ (−1, ½). With pure exponentials
  these are theorems given the two δ conditions and α ≥ 0, so a failure is a logic error (bad constants
  or a bad α curve), never an expected algorithmic state — the always-active tier is right. R4 checks the
  **same** window. (Rev. 5 had `min(T₂, T_max)` here and `T_max` in R4; round 2 caught the mismatch.)
- **Anchor consistency.** E₂ and ν₂ are the **isothermal-converted** values of one dataset at its own
  room-temperature row (Cu 295 K, Ni/Al/Pb/Fe/Sn 300 K, In 295 K, Ag per its table, Cr 250 K), and T₂ is
  the same number in the R1 fit and the constructor; the `create_mech` docstring says "isothermal, as
  converted in R1" so nobody pastes a handbook adiabatic value. Mixing a Blanke E with a Wolfram ν at
  "room temperature" is the present practice (`cl_Material_Copper.cpp:62-63`).
- **Cryogenic branch required.** `l( 0 )` is finite for any α spline, so `isfinite` is not the guard.
  R2 asserts `alpha_switch_temperature() > 0` (a cryogenic branch is installed, `cl_Material_SplineLookupTable.cpp:264-265`)
  and `l( 0 ) > 0`, and caches `l( 0 )`.
- **Hill averaging done right in R1.** For cubic crystals K is unique, (C₁₁+2C₁₂)/3; only G is averaged,
  G_H = (G_V + G_R)/2 with G_V = (C₁₁−C₁₂+3C₄₄)/5 and G_R = 5(C₁₁−C₁₂)C₄₄/(4C₄₄+3(C₁₁−C₁₂)). White tin is
  tetragonal: six constants, Voigt and Reuss bounds from the full stiffness/compliance tensors, then the
  Hill mean. Three temperatures for two exponents leaves one residual — Sn is the weakest dataset and R1
  says so in its table.

## 4. Ordered Steps

- [◐] **R1 — Data and fits (offline, Python, `tmp/`). DONE 2026-09-15 except O5–O7** (`tmp/r1_elastic/datasets.py`,
  `fit_r1.py`, `r1_results.json`, `<metal>_0_300.txt` from `material <metal> -t 0 300 1` on the debug binary of 16:29;
  the three materials sources newer than that binary are the peer's warning fixes, not curve changes). Needs from Christian: a reconfigured tree and a
  fine-grid `material <metal>` table (or a one-line probe printing `l( T )`, `alpha( T )`, `cp( T )`,
  `density( T )` at 1 K steps) for all nine metals. Per metal: collect the dataset (O1), Hill-average
  single-crystal sets (3.1), convert every K_S datum to K_T with the identity and BELFEM's α, ρ, c_p at the
  datum's temperature, fit δ_K on ln K_T and δ_G on ln G vs the code's ε over 0–300 K, record residuals,
  δ_G − δ_K, the 0 K limits, the isothermal anchor (T₂, E₂, ν₂), the room-temperature E and ν change vs
  the current tree, and the α range vs T_max. Also probe ν(4/77/293 K) on the **current** tree for all
  nine, so the "inverted for every metal" hypothesis is measured. Output: §6 table.
- [x] **R2 — `Metal::create_mech` rewrite + nine call sites, one compilable step** (after: R1). **DONE 2026-09-15**:
  implemented, code-jury findings fixed (§8), built clean in the debug tree. New
  signature `create_mech( E2, nu2, T2, deltaK, deltaG )`; E₂ GPa → Pa at entry. Preconditions, all
  `BELFEM_ERROR`: `have( alpha )`, `have( cp )`, `have( ref_density )`, `alpha_switch_temperature() > 0`,
  finite E₂ > 0, −1 < ν₂ < ½, T₂ finite and inside the spline range, finite δ_K > 0, δ_G > δ_K, `l( 0 ) > 0`,
  finite K₀, G₀ > 0. Cache `l( 0 )`. New `Metal::eps`, `E_custom`, `nu_custom`; `set_have( E )`,
  `set_have( nu )` **before** `create_spline( E, 0.0 )` and `create_spline( nu, 0.0 )`; then the sampled
  whole-range checks of 3.1; `grueneisen` constant per row 12. Remove `mWachtmanYoung`, `Metal::dEdT_custom`,
  the γ-inversion loop, and Nickel's `create_young`/`E_custom`/`dEdT_custom`/`mYoungData`/`mYoungBezier`.
  `cp_from_debye` and `grueneisen( T )` stay as they are. Update the nine constructors with R1's numbers
  and a source + DOI comment each; Nickel gets the §6a header block.
- [x] **R3 — Probe** (after: R2). **DONE 2026-09-15** (`tmp/r1_elastic/probe_r3.sh`, table below). `material <metal>` for all nine: E, ν at 4 / 77 / 293 K,
  ν monotone over the grid, E(293 K … T_max) against the retired Wachtman curve — deviations tabulated
  into §6. One formula alloy (e.g. `Sn60Pb40`) before/after: E, ν, α at 77 and 293 K.
- [x] **R4 — Regression test `tests/physics/test_MetalElastic.cpp`** (after: R1, R2). **DONE 2026-09-15**, 4 tests, all
  green in the debug tree; oracles in SI, served (isothermal) convention:
  - every pure metal: ν(T_{k+1}) ≥ ν(T_k) − 1e-9 on the 4 K grid over [4 K, T_max]; ν(8 K) − ν(4 K) < 2e-5;
    ν ∈ (−1, ½); E(4 K) > E(293 K); K(4 K)/K(293 K) within 2 % of the source's isothermal-converted ratio
    recorded per metal in §6 (rev. 5's universal 1.10 bound was already false for Indium's adiabatic data,
    B(0)/B(295) = 1.110, Kim & Ledbetter 1998, 10.1016/S0921-5093(98)00490-0);
  - copper (Ledbetter 1981, 10.1002/pssa.2210660209, Table 2, converted): G(5 K) = 51.72 ± 0.5 GPa,
    G(295 K) = 47.57 ± 0.5 GPa (G is state-independent), K(5 K) = 144.46 ± 1.5 GPa (conversion vanishes),
    K(295 K) = 135.7 ± 1.5 GPa, ν(295 K) = 0.343 ± 0.002, E(295 K) = 127.8 ± 1.0 GPa;
  - one formula alloy constructs and gives ν ∈ (−1, ½) at 77 and 293 K.
  Register in `tests/physics/CMakeLists.txt`; runs under `make check` (label `fast` only if wall time < 10 s).
- [x] **R5 — Docs** (after: R2). **DONE 2026-09-15**: drafted (usage guide §9.7 + three passages, property sources table, bibliography
  and gaps, thermal-expansion §5, `doc/literature_references.md` Tier 6, class headers Al/Cr/Ni/Pb); Codex language sweep pending. `materials_usage_guide.md` §9.7 rewritten for the new closure ("isothermal"
  replaces "dynamic (adiabatic) … not corrected"), lines 119/254/271-272 updated; `material_property_sources.md`:
  E and ν columns name the fit source with DOI per metal, plus the sentence that the pure-metal roster
  represents **isotropic polycrystals** (Hill averages where the source is a single-crystal set), and the
  before/after room-temperature table from §6; `thermal_expansion_from_heat_capacity.md` §5 updated (stale
  `:1105` anchor replaced by a greppable token); class headers Al, Cr (α clamp 1700–2180 K next to the
  Néel text), and the Nickel block (§6a); `doc/literature_references.md` gains the O1 sources with DOIs.
  Codex language sweep on the touched sections (luna/medium).
- [◐] **R6 — Gate** (after: R4, R5). `test_physics` 14/14 and **`make check` 17/17 passed** 2026-09-15 (debug tree);
  `check_doc_claims.py` 38/38 after the doc edits. Open: Christian's commit, then the move to `todo/closed/` with the Lead
  follow-up (O7 polycrystal dataset) split into its own note. `make check` green with the new test; `scripts/check_doc_claims.py`
  after the doc edits. Then close: devlog, register line, this file to `todo/closed/`.

## 5. Design Questions — all resolved by Christian, 2026-09-15

- **O1 — Data source per metal. RESOLVED (papers delivered to `tmp/papers/`, `tmp/papers2/`; DOIs read
  from the files or Crossref-resolved where the PDF lacks one).**

  | metal | source | DOI | data | note |
  |---|---|---|---|---|
  | Cu | Ledbetter 1981, phys. stat. sol. (a) 66, 477 | 10.1002/pssa.2210660209 | polycrystal B, G, E, ν, 5–295 K, 16 pts (Table 2) | complete; Varshni fits in Table 4 |
  | Cu (check) | Ledbetter & Naimon 1974, JPCRD 3, 897 | 10.1063/1.3253150 | RT averages (Table 10), single-crystal C_ij(T) compilation | RT anchor cross-check |
  | Al | Kamm & Alers 1964, J. Appl. Phys. 35, 327 | 10.1063/1.1713309 | single crystal 0–300 K, 20 K steps, incl. B (Table I) | Hill |
  | Ag | Neighbours & Alers 1958, Phys. Rev. 111, 707 | 10.1103/PhysRev.111.707 | C₄₄, C′, C_L 4.2–300 K (figures/tables) | Hill; digitize if no table — quality risk, not a blocker |
  | Pb | Waldorf & Alers 1962, J. Appl. Phys. 33, 3266 | 10.1063/1.1931149 | single crystal 0–300 K, 20 K steps (Table I) | Hill |
  | Sn | Rayne & Chandrasekhar 1960, Phys. Rev. 120, 1658 | 10.1103/PhysRev.120.1658 | six tetragonal constants at 300, 77, 4.2 K (Table IV) | three points, two exponents, one residual — weakest set |
  | In | Kim & Ledbetter 1998, Mater. Sci. Eng. A 252, 139 | 10.1016/S0921-5093(98)00490-0 | polycrystal C_l, G, B, E, ν as Varshni fits 5–300 K (Table 1, Eq. 6) | complete via the fit: B 46.99 → 42.33 GPa, G 6.84 → 4.39, ν 0.4306 → 0.4498 (probe) |
  | Cr | Palmer & Lee 1971, Phil. Mag. 24, 311 | 10.1080/14786437108227390 | single crystal 4.2–345 K, 5 K steps, zero field (table p. 317) | anomalies at 123 K and 311 K in the data; fit outside them; anchor at 250 K |
  | Cr (check) | Bolef & de Klerk 1963, Phys. Rev. 129, 1063 | 10.1103/PhysRev.129.1063 | 77, 298, 500 K (Table II) | high-T check |
  | Fe | Rayne & Chandrasekhar 1961, Phys. Rev. 122, 1714 | 10.1103/PhysRev.122.1714 | smoothed single crystal 4.2–300 K, 20 K steps (Table I) | Hill |
  | Ni | Alers, Neighbours & Sato 1960, J. Phys. Chem. Solids 13, 40 | 10.1016/0022-3697(60)90125-6 | single crystal at 10 kOe (saturated), Table 2, 20 K steps 0–760 K | fit 0–300 K, anchor 300 K; the paper notes a residual intrinsic magnetic contribution below T_C even at saturation — absorbed by the fit, stated in §6a |
  | Fe, Ni (check) | Ledbetter & Reed 1973, JPCRD 2, 531 | 10.1063/1.3253127 | compilations by source (Tables 5–19), §15 on the ΔE effect | RT cross-check only; Table 6 holds RT best values, not a T-series |
  | theory | Garai & Laugier 2007, J. Appl. Phys. 101, 023514 | 10.1063/1.2424535 | B_T = B_T(0)·exp(−∫α_V δ_T dT), Eq. 8–10 | the spine — for the **isothermal** K, hence the R1 conversion |
  | context | Köster & Franz 1961, Metall. Rev. 6, 1 | 10.1179/mtlr.1961.6.1.1 | ν(T) of metals, §VI | direction only |
  | context | Nadal & Le Poac 2003, J. Appl. Phys. 93, 2472 | 10.1063/1.1539913 | G(T) to the melting point, G(T_m) = 0 | row 9; not adopted |

  Anderson 1966 was not delivered; the "Wachtman from Grüneisen" attribution is not relied on anywhere.

- **O2 — Adiabatic or isothermal. RESOLVED → isothermal** (Christian: serve the moduli relevant for
  mechanical elastic calculations; BELFEM's mechanics is quasi-static, so isothermal — Landau & Lifshitz,
  *Theory of Elasticity*, §6, textbook, recalled). **Where the conversion lives (rev. 6, implementation
  decision):** in R1, on the data, with the exact identity K_T = K_S/(1 + α_V² T K_S/(ρ c_p)) and
  BELFEM's α, ρ(T) = ρ_ref/l³, c_p at the datum's temperature. G needs no conversion. The served curves are
  pure exponentials; `grueneisen( T )` and `cp_from_debye`, which already treat K as K_T, are correct as
  written. Copper at 295 K: K_T = 135.67 GPa (K_S 139.74), E 127.78, ν 0.3430 (adiabatic 128.17 / 0.3471);
  at 5 K identical in both. The constructor's E₂, ν₂ are therefore isothermal numbers, and the R4 oracles
  compare G directly and K, E, ν against the converted Ledbetter row.
- **O3 — Anomalous metals. RESOLVED** (Christian: "agreed", after `tmp/ai_exchange/review_nickel_elastic.md`;
  jury 3/3 for the smooth model). Nickel does **not** keep its own E(T): (1) the dip is the ΔE effect of the
  demagnetized state and vanishes at saturation (Ledbetter & Reed 1973 §15, 10.1063/1.3253127); BELFEM's
  nickel sits in tesla fields; (2) the anchor is state-dependent by ~18 % (190.5 GPa at H = 0 vs 225.6 GPa
  at 6.2 kOe on one annealed sample, Giebe & Blechschmidt 1931 in L&R Table 9; compiled RT values
  170–231 GPa); (3) the coded curve is already −18 % at 400 K and has a sign change of dE/dT at 487 K — a
  design hazard for a future thermal-stress tangent (none is wired today: `IwgType::PlaneStress` /
  `LinearElasticity` are enum values only, `en_IWGs.hpp:71-72`; `src/fem/postproc/fn_FEM_mises_planestress.cpp:64`
  is a `BELFEM_ERROR` stub); (4) option A would have glued an anelastic demagnetized E to a saturated K,
  producing a ν of neither state. Saturation removes the domain-wall contribution; a smaller intrinsic
  magnetic contribution remains below T_C (Alers, Neighbours & Sato 1960, results section) and is absorbed
  by the fit. What is lost: an unmagnetized nickel part near 500 K is served an E about 40 % too stiff —
  documented in the Nickel header (§6a), Chromium precedent (`cl_Material_Chromium.hpp:28-72`). Chromium:
  smoothed by class design. Iron: never served near T_C.
- **O4 — Strain and functional form. RESOLVED** (Christian: "sounds good"). ε = 3∫₀ᵀ α dT = ln(V/V₀) and
  M₀·exp(−δ ε) = M₀(V/V₀)^(−δ), the exact integral of "d ln M / d ln V = −δ constant". Rejected: the
  engineering strain V/V₀ − 1 (differs by ε²/2; the first draft mixed the two) and the linear form
  M₀(1 − δε) (δ no longer a constant log-derivative, can cross zero, and for indium δ_G·ε ≈ 0.45 at room
  temperature gives 0.64 G₀ exponential vs 0.55 G₀ linear). R1's fit is a straight-line regression of
  ln M on ε.

- **O5 — Iron's isothermal ν is flat, not rising. RESOLVED 2026-09-15 (Christian: option (i)).** Rayne & Chandrasekhar's zero-field data give adiabatic
  ν 0.2849 → 0.2884 (4 → 300 K), but after the isothermal conversion ν_T(300) = 0.2852: the whole rise was the
  adiabatic K_S. The fit returns δ_K = 7.91 > δ_G = 7.84, so the served ν_T drifts from 0.2845 to 0.2844 —
  a change of 1e-4, below any measurement. The strict guard δ_G > δ_K of 3.1 would refuse Iron at
  construction. Options: (i) relax the guard to a tolerance on the **sampled** ν — no decrease larger than
  1e-3 anywhere on the grid, and ν(T_max) ≥ ν(0) − 1e-3 — and drop the strict δ ordering (keep δ_K, δ_G > 0);
  (ii) nudge Iron's δ_G to δ_K + 0.1 by hand (a fudge, and it would be visible in the R1 residuals).
  Recommendation: (i) — "monotone within data accuracy" is what the physics supports for a metal whose
  isothermal ν is genuinely flat, and the guard still catches a sign error (a real inversion is 1e-2, not 1e-4).
- **O6 — Chromium's data do not support the quasi-harmonic model between 120 and 330 K. RESOLVED 2026-09-15 (Christian: option (i), constant ν).** Palmer & Lee's
  zero-field constants (C₁₂ derived as C₁₁ − 2C′) give K rising from 188 GPa at 120 K to 192.5 at 175 K, then
  collapsing to 151 GPa at the Néel point (310 K) and recovering to 171 GPa in the paramagnetic phase; G drops
  from 121 to 115 GPa across the spin-flip region (120–175 K) and is flat above; measured ν falls from 0.237
  (0–120 K) to 0.196 at 310 K, 0.225 at 330 K. Every fit window gives a served ν that decreases somewhere,
  and because ε(300 K) is only 0.32 % the δ values are ill-conditioned (δ_K 38 on 0–120 K, 75 on 175–300 K).
  The class already declares the anomalies unrepresented (`cl_Material_Chromium.hpp:28-72`). Options:
  (i) **constant ν**: δ_G = δ_K = 15.8 (the G fit over 0–120 K plus the paramagnetic 315–330 K points), K₀ =
  190.06, G₀ = 121.19 GPa from the cryogenic plateau, hence ν = 0.2367 (isothermal ≈ 0.236) at every T,
  anchor 100 K; smooth, monotone-within-tolerance under O5(i), documented as "SDW anomalies in K and ν not
  represented; ν held at its cryogenic value"; (ii) fit the paramagnetic phase (315–345 K) and accept a
  10–15 % error in K below the Néel point — wrong where BELFEM uses chromium (cryogenic); (iii) keep the
  present constant-γ closure for Cr alone — no, it inverts ν there too (0.216 → 0.210). Recommendation: (i).
  Note the room-temperature value the tree serves for Cr, ν = 0.21 ("Wolfram Cloud"), is neither the
  cryogenic 0.237 nor the paramagnetic 0.225.
- **O7 — Lead: the Hill average of a soft, strongly anisotropic crystal over-estimates the polycrystal. RESOLVED 2026-09-15, final ruling: static level** (crystal K from Waldorf & Alers, E from Blanke's static-type curve; E 16.25, G 5.66 GPa, ν 0.435 at 300 K; constants 16.25 / 0.4347 / 300 / 7.31 / 11.20). Earlier the same day: (i) with (ii) as interim.
  Waldorf & Alers' single-crystal constants (Zener anisotropy 2C₄₄/(C₁₁−C₁₂) = 4.1) Hill-average to
  G(300 K) = 8.51 GPa and E = 24.0 GPa, against the tree's 16.5 GPa and the handbook's ≈ 16 GPa; the Reuss
  bound alone is 6.6 GPa. Kim & Ledbetter 1998 (10.1016/S0921-5093(98)00490-0) measured the same effect on
  indium (anisotropy 4.7): the real polycrystal's G and E were 8 % below the VRH prediction. For lead the
  gap is larger, and Waldorf & Alers' crystal carried 0.02 % Bi to pin dislocations, so pure lead is softer
  still. Options: (i) find a polycrystalline ultrasonic dataset for lead (Ledbetter's group measured several
  soft metals in the 1980s — a paper request); (ii) use the Hill numbers and document the +45 % change;
  (iii) scale the Hill G to the tree's room-temperature E (a fudge). Recommendation: (i), with (ii) as the
  interim so the plan is not blocked. The same caution applies, more mildly, to tin (anisotropic tetragonal,
  three temperatures only).

## 6. Interface Design and Constants

```cpp
// Metal (protected) — quasi-harmonic elastic constants, all at 0 K, Pa and dimensionless
real mK0     = BELFEM_QUIET_NAN ;   // isothermal bulk modulus at 0 K [Pa] ( equals the adiabatic one there )
real mG0     = BELFEM_QUIET_NAN ;   // shear modulus at 0 K [Pa]
real mDeltaK = BELFEM_QUIET_NAN ;   // Anderson-Grueneisen softening constant of the isothermal bulk modulus
real mDeltaG = BELFEM_QUIET_NAN ;   // softening constant of the shear modulus ( > mDeltaK )
real mL0     = BELFEM_QUIET_NAN ;   // l( 0 K ), cached: eps( T ) = 3 ln( l( T ) / mL0 )

/**
 * Quasi-harmonic elastic moduli ( Garai & Laugier 2007, 10.1063/1.2424535, Eq. 10, constant delta ):
 *   eps( T ) = 3 int_0^T alpha dT = 3 ln( l( T ) / l( 0 ) )
 *   K( T ) = K0 exp( -deltaK eps ),  G( T ) = G0 exp( -deltaG eps )      ( K isothermal )
 *   E = 9 K G / ( 3 K + G ),  nu = ( 3 K - 2 G ) / ( 6 K + 2 G )
 * Poisson's ratio rises monotonically with T iff deltaG > deltaK ( enforced ); dnu/dT( 0 ) = 0 exactly.
 * The anchor values are ISOTHERMAL, as converted from the ultrasonic ( adiabatic ) data in the offline fit
 * that produced deltaK and deltaG; the served E( T2 ) equals E2.
 * @param E2      isothermal Young's modulus at the anchor temperature [GPa], converted to Pa here
 * @param nu2     isothermal Poisson's ratio at the anchor
 * @param T2      anchor temperature [K] — the same T2 the offline fit used
 * @param deltaK  softening constant of the isothermal bulk modulus ( > 0 )
 * @param deltaG  softening constant of the shear modulus ( > deltaK )
 */
void create_mech( const real E2, const real nu2, const real T2, const real deltaK, const real deltaG );

real eps( const real T ) const ;                    // 3 ln( l( T ) / mL0 )
real E_custom( const real T ) const override ;      // 9KG/(3K+G) from the four constants
real nu_custom( const real T ) const override ;     // (3K-2G)/(6K+2G)
// no dEdT_custom: the 0 K spline slope is the literal 0.0 passed to create_spline( E, 0.0 ), exact since alpha( 0 ) = 0
```

Naming: the math/physics exemption applies (`E2, nu2, T2, deltaK, deltaG`, `eps`); members keep the `m`
prefix. `eps` is the strain, not `BELFEM_EPSILON`.

Per-metal constants table — filled by R1 (pilot rows from the coarse-ε probe; refit against the code's ε).
E₂, ν₂ are the isothermal-converted anchor; "tree → new" is the served room-temperature E:

| metal | source (DOI) | data | T₂ [K] | E₂ [GPa] | ν₂ | δ_K (isoth.) | δ_K (adiab.) | δ_G | K₀ [GPa] | G₀ [GPa] | rms K / G | K(4)/K(293) model / src | E(293) tree → new | ν(4/77/293) tree → new | K_S/K_T−1 at 293 K |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Cu | Ledbetter 1981, 10.1002/pssa.2210660209 | 16 pts, 5–295 K, polycrystal | 295 | 127.45 | 0.3433 | 6.47 | 3.43 | 8.60 | 144.24 | 51.52 | 0.11 % / 0.24 % | 1.064 / 1.064 | 133.0 → 127.6 | 0.363/0.361/0.344 → 0.3404/0.3406/0.3433 | 3.0 % |
| In | Kim & Ledbetter 1998, 10.1016/S0921-5093(98)00490-0 | Varshni fits 5–295 K, polycrystal | 295 | 12.85 | 0.4464 | 7.61 | 5.01 | 20.56 | 47.01 | 6.87 | 0.04 % / 0.43 % | 1.174 / 1.173 | 12.8 → 12.9 | 0.433/0.435/0.450 → 0.4303/0.4330/0.4463 | 5.7 % |
| Pb | Waldorf & Alers 1962, 10.1063/1.1931149 | 16 pts, 0–300 K (cubic Hill) | 300 | 23.84 | 0.4043 | 7.31 | 4.05 | 13.23 | 48.65 | 11.30 | 0.15 % / 0.34 % | 1.166 / 1.167 | **16.5 → 24.0 (+45 %, O7)** | 0.448/0.445/0.440 → 0.3922/0.3942/0.4040 | 7.2 % |
| Ag | Neighbours & Alers 1958, 10.1103/PhysRev.111.707 | 14 pts, 0–300 K (cubic Hill) | 300 | 80.40 | 0.3652 | 6.97 | 3.73 | 8.58 | 108.49 | 32.81 | 0.15 % / 0.30 % | 1.089 / 1.089 | 84.2 → 80.7 | 0.384/0.382/0.370 → 0.3627/0.3629/0.3651 | 4.1 % |
| Al | Kamm & Alers 1964, 10.1063/1.1713309 | 16 pts, 0–300 K (cubic Hill) | 300 | 69.51 | 0.3409 | 6.88 | 3.50 | 9.72 | 79.22 | 29.20 | 0.15 % / 0.25 % | 1.084 / 1.085 | 73.4 → 69.8 | 0.364/0.363/0.350 → 0.3359/0.3361/0.3407 | 4.1 % |
| Fe | Rayne & Chandrasekhar 1961, 10.1103/PhysRev.122.1714 | 16 pts, 4.2–300 K (cubic Hill, zero field) | 300 | 212.20 | 0.2844 | 7.91 | 4.81 | 7.84 | 172.26 | 86.72 | 0.34 % / 0.15 % | 1.048 / 1.049 | 213.1 → 212.6 | 0.330/0.325/0.290 → 0.2845/0.2845/0.2844 (**flat, O5**) | 1.9 % |
| Ni | Alers, Neighbours & Sato 1960, 10.1016/0022-3697(60)90125-6, Table 2 (10 kOe) | 16 pts, 0–300 K (cubic Hill) | 300 | 222.35 | 0.2935 | 6.14 | 3.18 | 10.88 | 187.22 | 92.63 | 0.14 % / 0.20 % | 1.042 / 1.042 | 184.0 → 223.0 (+21 %) | 0.304/0.302/0.310 → 0.2876/0.2879/0.2933 | 2.0 % |
| Cr | Palmer & Lee 1971, 10.1080/14786437108227390 (zero field; C₁₂ = C₁₁ − 2C′) | 62 pts 0–330 K; **no window fits (O6)** | (100) | (297.5) | (0.2367) | 38.0 (0–120 K) | 34.9 | 30.8 (0–120 K); 15.8 (0–120 + 315–330 K) | 190.06 | 121.19 | 0.14 % / 0.05 % (0–120 K) | 1.125 / — | 278.7 → ≈ 272–282 | 0.216/0.215/0.210 → ν **falls** in every window | 0.5 % |
| Sn | Rayne & Chandrasekhar 1960, 10.1103/PhysRev.120.1658 | 3 pts (tetragonal VRH) | 300 | 48.29 | 0.3471 | 7.12 | 4.08 | 21.54 | 58.26 | 24.34 | 0.54 % / 0.18 % | 1.103 / 1.097 | 55.3 → 48.7 | 0.377/0.374/0.360 → 0.3167/0.3205/0.3462 | 4.3 % |

ε(300 K) [%]: Cu 0.99, In 2.17, Pb 2.16, Ag 1.26, Al 1.23, Fe 0.62, Ni 0.69, Cr 0.32, Sn 1.42. The isothermal
δ_K exceeds the adiabatic one by 2.9–3.3 for every metal (the conversion factor grows ∝ T above θ_D, so it
looks like extra softening in ε). α range vs T_max: Pb 600 vs 600.61 K and Al 933.45 vs 933.47 K are sub-kelvin
overhangs (fine); Ni 600 vs 631 and Cr 1700 vs 2180 are the two clamps (row 10).

**Current-tree probe (ν at 4 / 77 / 293 K, column "tree"):** the constant-γ closure inverts ν for Cu, Pb, Ag,
Al, Fe (0.330 → 0.290), Sn and Cr; Indium (0.433 → 0.450) and Nickel (0.304 → 0.310) happened to come out
rising. Seven of nine, measured — the rev. 1 hypothesis "all nine" was wrong for two.

**R3 probe (2026-09-15, new debug binary, `material <metal> -t 0 T_max 1`).** All nine construct; served ν
non-decreasing everywhere (max drop 0 on the 1 K grid). Served values and the deviation of the new E(T) from the
retired Wachtman curve above room temperature (row 9: reported, not gated):

| metal | E(4) | E(77) | E(293) [GPa] | ν(4) | ν(77) | ν(293) | ν(T_max) | new E vs Wachtman, 293 K … T_max |
|---|---|---|---|---|---|---|---|---|
| Cu | 138.14 | 137.40 | 127.56 | 0.3400 | 0.3410 | 0.3430 | 0.362 (1358 K) | −7.5 % … +34 % |
| Al | 78.00 | 77.58 | 69.82 | 0.3360 | 0.3360 | 0.3410 | 0.362 (933 K) | −8.9 % … −4.3 % |
| Ag | 89.41 | 88.45 | 80.67 | 0.3630 | 0.3630 | 0.3650 | 0.377 (1235 K) | −5.7 % … +24 % |
| In | 19.66 | 18.44 | 12.90 | 0.4300 | 0.4330 | 0.4460 | 0.454 (430 K) | +0.7 % … +10 % |
| Pb | 31.46 | 30.10 | 24.03 | 0.3920 | 0.3940 | 0.4040 | 0.418 (601 K) | +46 % … +85 % (O7 interim) |
| Sn | 64.09 | 61.99 | 48.74 | 0.3170 | 0.3200 | 0.3460 | 0.375 (505 K) | −12 % … +18 % |
| Fe | 222.76 | 222.36 | 212.63 | 0.2840 | 0.2840 | 0.2840 | 0.284 (860 K) | −1.6 % … +6.6 % |
| Cr | 299.84 | 299.36 | 285.58 | 0.2370 | 0.2370 | 0.2370 | 0.237 (2180 K) | −49 % … +2.5 % |
| Ni | 238.56 | 237.75 | 222.96 | 0.2880 | 0.2880 | 0.2930 | 0.306 (631 K) | — (no Wachtman before) |
| Sn60Pb40 | 52.33 | 50.45 | 39.79 | 0.3420 | 0.3450 | 0.3650 | 0.376 (400 K) | — |

Copper's printed Grüneisen parameter is now 2.03 (was 2.13; literature 1.96–2.00). The high-temperature
deviations are the expected consequence of a 0–300 K fit extrapolated along ε: copper's new E at the melting point
is 34 % above the Wachtman value, chromium's common δ = 15.8 softens it to half the old value at 2180 K. Neither
range is in BELFEM's mechanical use; both are documented, not gated.

**R1 findings that need a ruling (O5–O7, §5).**

### 6a. Nickel header block (drafted for R2/R5; wording to be kept with the code)

```
 * ELASTIC MODULI - WHAT IS AND IS NOT REPRESENTED
 *
 * The moduli follow the same quasi-harmonic closure as the other metals
 * ( Metal::create_mech ): K and G soften with the volumetric thermal strain,
 * fitted to the single-crystal constants of Alers, Neighbours and Sato 1960
 * ( 10.1016/0022-3697(60)90125-6 ), measured in a saturating field of 10 kOe,
 * Hill-averaged to the isotropic polycrystal; room-temperature cross-check
 * against Ledbetter and Reed 1973 ( 10.1063/1.3253127, Table 6 ).
 *
 * Nickel's Young's modulus in the DEMAGNETIZED state is lower and shows a deep
 * minimum below the Curie point ( Blanke 1989: 183.7 GPa at 293 K, 124 GPa at
 * 487 K, 195.7 GPa at 637 K; the fit carried here until 2026-09 was already
 * 18 % below its room-temperature value at 400 K ). That is the Delta-E
 * effect: domain walls move under stress and add strain; at magnetic
 * saturation the walls are immobile ( Ledbetter and Reed 1973, section 15 ).
 * The effect depends on field, stress and annealing state ( 190.5 GPa at H = 0
 * against 225.6 GPa at 6.2 kOe on one sample, Giebe and Blechschmidt 1931,
 * ibid. Table 9 ), so it is a property of the magnetic state, not of the
 * lattice. BELFEM's nickel parts sit in tesla-level fields and are saturated;
 * the saturated moduli are therefore the ones served, and the Delta-E dip is
 * deliberately not represented. A smaller intrinsic magnetic contribution
 * remains below the Curie point even at saturation ( Alers et al. 1960 ) and
 * is absorbed by the fit. Until 2026-09 the class carried Blanke's demagnetized
 * curve as a two-branch Wachtman fit with a Bezier bridge; it was retired
 * because it describes a state the solver does not simulate, its magnitude is
 * uncertain by tens of percent, and its sign change of dE/dT inside 465-631 K
 * would be a hazard for the Newton tangent of a thermal-stress run.
 * Consequence: an unmagnetized nickel part near 500 K is served an E about
 * 40 % too stiff, and the room-temperature E rose from 184 to about 223 GPa.
 * The thermal expansion Bezier ends at 600 K and is clamped up to T_max.
 * Decision: Christian Messe, 2026-09-15.
```

Present anchors, for the record (all `create_mech( E0, b, T0, T2, nu2 )`): Al 75.706 / 0.07071 / 649.62 /
293.15 / 0.35; Cr 279.837 / 0.08741 / 900.83 / 293.15 / 0.21; Cu 134.447 / 0.11492 / 922.14 / 293.15 /
0.344; Ag 86.647 / 0.06298 / 587.35 / 293.15 / 0.37; In 19.603 / 0.02701 / 44.60 / 295 / 0.4498;
Pb 19.011 / 0.03196 / 383.6 / 293.15 / 0.44; Sn 63.532 / 0.21789 / 600.08 / 293.15 / 0.36;
Fe 214.886 / 0.18823 / 998.94 / 293.15 / 0.29; Ni 0 / 0 / 0 / 293.15 / 0.31 (+ `create_young`).

## 7. Definition-of-Done Checklist

- [ ] Every gap-table row mapped to a step or a resolved question.
- [ ] Each claimed gap backed by a citation; rows 9 and 10 carry marked assumptions.
- [ ] Ordered steps with dependencies; R2 is one compilable step.
- [x] Open questions logged, not silently decided (O1–O4 all resolved by Christian, 2026-09-15).
- [ ] End-to-end gate: `make check` with `test_MetalElastic` green on the R4 oracles (served isothermal
  convention, stated tolerances); `material copper` shows ν rising monotonically 4 K → T_max and
  E(295 K) = 127.8 ± 1.0 GPa, ν(295 K) = 0.343 ± 0.002.

## 8. Audit Trail

- Rev. 6 (2026-09-15): round 2 at xhigh on rev. 5 (`tmp/ai_exchange/review_elastic_quasiharmonic.md`,
  Codex terra/xhigh, Grok 4.6/xhigh, blind) — "not R2-ready", architecture endorsed by both. Folded in:
  - **P0 (both):** the isothermal law was cited but the adiabatic K_S was fitted and converted at serve time
    → convert the data in R1, fit and serve K_T (§2). Codex named the fix; Grok added the 0/0 at T = 0 and
    the c_p-wiggle consequence, both now moot.
  - **P0 (both):** scope, pilot and DoD quoted adiabatic numbers as the served result → all served numbers
    are isothermal (Cu 127.8 / 0.343 at 295 K); R4 oracles converted, G compared directly.
  - **P1 (Codex):** ∂E/∂K = 9G²/(3K+G)², not 27G² — and no analytic derivative is needed: `dEdT_custom` is
    never called, the spline gets the literal 0.0 (Grok) → `dEdT_custom` removed (row 5, 14).
  - **P1 (both):** monotonicity window `min(T₂, T_max)` vs `T_max` → whole served range in both (3.1).
  - **P1 (Grok):** R4's universal K(4)/K(293) < 1.10 false for Indium (1.110 adiabatic) → per-metal
    data-driven bound; "stiffens" wording fixed.
  - **P2 (Codex):** R4 had no numeric tolerances → oracles in SI with DOIs; "nothing to extrapolate" →
    controlled cryogenic extrapolation (§2).
  - **Grok:** stale rows (Status "rev. 3", gap 6 "L&R Table 6", gap 3 "O1 fallback", "R7", "R6 alloy smoke",
    Ni window 0–520 vs ANS 0–760, T₂ 293 vs 300) → fixed; Chromium α clamp 1700 vs 2180 K and T₂ inside the
    Néel halo → row 10, row 7 (anchor 250 K); Hill "arithmetic mean" wording → 3.1 formulas; ANS's residual
    magnetic contribution at saturation → O3, §6a; `cl_Bezier` citation split (comment `.hpp:104`, clamp
    `.cpp:97-98`); `fn_FEM_mises_planestress.cpp` is in `src/fem/postproc/`.
  - Not adopted: Grok's "cap the construction check at min(T₂, T_max) and exempt Nickel's c_p wiggle" —
    moot once the fit is isothermal and the served curves are pure exponentials.
- Rev. 5 (2026-09-15): O4 confirmed; Alers, Neighbours & Sato 1960 delivered.
- Rev. 4 (2026-09-15): O3 resolved by Christian after the narrow-claim Nickel brief
  (`tmp/ai_exchange/review_nickel_elastic.md`, Codex terra/medium + Grok 4.6/high, 3/3 for the smooth
  model; the jury returned after the ruling and its refinements are in O3).
- Rev. 3 (2026-09-15): O1 and O2 resolved with the second paper batch.
- Rev. 2 (2026-09-15): round 1 (Codex terra/high, Grok 4.6/high) on rev. 1 — 20 findings, all confirmed:
  two ε definitions, chain rule, `Vector` bundle, O2 blocks R2, `nu_custom` absent, `set_have( nu )`,
  step dependencies, Ni test window, `check-fast` labels, copper RT anchor, docs list, `l( 0 )` guard, Alloy.
- Diagnosis round: `tmp/ai_exchange/review_copper_poisson.md` (Codex confirmed; Grok leg failed — CLI not
  signed in). Distilled into `devlog/dl20260915_copper_poisson_trend.md`.
- R2 implemented 2026-09-15 (Claude, on approval). Constants as served: Cu (127.45, 0.3433, 295, 6.47, 8.60);
  In (12.85, 0.4464, 295, 7.61, 20.56); Pb (23.84, 0.4043, 300, 7.31, 13.23, INTERIM); Ag (80.40, 0.3652, 300, 6.97, 8.58);
  Al (69.51, 0.3409, 300, 6.88, 9.72); Fe (212.20, 0.2844, 300, 7.91, 7.84); Ni (222.35, 0.2935, 300, 6.14, 10.88);
  Cr (298.67, 0.2371, 100, 15.8, 15.8 — constant ν); Sn (48.29, 0.3471, 300, 7.12, 21.54). Guard: sampled ν may
  drop at most 1e-3 (O5); α ≥ 0; ν ∈ (−1, ½). Code jury: `tmp/ai_exchange/review_elastic_r2.md`.
- R1 executed 2026-09-15 (Claude, sandbox; `tmp/r1_elastic/`): nine datasets transcribed with DOIs, Hill
  averages (cubic formulas; tetragonal VRH for Sn), isothermal conversion with the tree's α, c_p, ρ_ref/l³ from
  `material <metal> -t 0 300 1`, exponential fits. Findings O5–O7 filed. Seven of nine current-tree ν curves
  are inverted; Indium and Nickel were not.
- **Round 3?** Protocol §11: a new round needs a material change — rev. 6 has one (where the conversion
  lives). It is a specification change within an endorsed architecture, so a round 3 at xhigh is
  Christian's cost call; the alternative is to go to R1 and let R3/R4 be the gate.
