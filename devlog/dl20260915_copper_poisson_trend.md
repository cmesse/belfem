# Devlog 2026-09-15 — Copper Poisson ratio falls with temperature: the constant-γ closure

**Date:** 2026-09-15
**Topic:** Why `material copper` prints ν = 0.363 at 0 K falling to 0.344 at 293 K, and what a thermophysically consistent ν(T) would need
**AIs involved:** Claude (Fable 5.1), Codex (gpt-5.6-terra, high); Grok leg failed (CLI not signed in)
**Claude Confidence:** high (mechanism), medium (literature magnitudes)
**Codex Audit Confidence:** high
**Literature References:** Ledbetter 1981, phys. stat. sol. (a) 66, 477, 10.1002/pssa.2210660209 (polycrystalline Cu 5–295 K, Table 2 read); Ledbetter & Naimon 1974, JPCRD 3, 897, 10.1063/1.3253150 (Table 10); Garai & Laugier 2007, J. Appl. Phys. 101, 023514, 10.1063/1.2424535 (Eq. 8–10); Köster & Franz 1961, Metall. Rev. 6, 1, 10.1179/mtlr.1961.6.1.1 (§VI); single-crystal C_ij of Cu at 300 K (arXiv:1605.09237 Table 1, its ref. 45) and the low-T set (NIST atomman page) used for the first Hill estimate before the papers arrived
**Verification:** focused regression — `test_physics` 14/14 green (4 new `MetalElastic` tests) on `cmake-build-debug` @ devel + this diff, 2026-09-15; R3 probe on all nine metals; `make check` 17/17 passed (same tree, 2026-09-15); `check_doc_claims.py` 38/38. Read-only session; a peer session was editing warning-related code in the same checkout.

## Summary

Christian noticed that the copper report shows Poisson's ratio decreasing with temperature and
expected the opposite. He is right. The trend is not a coding defect but a consequence of the
ν closure in `Metal::create_mech` (`cl_Material_Metal.cpp:1122-1230`): one room-temperature
anchor (T₂, ν₂) fixes a Grüneisen parameter γ, γ is held constant, and the bulk modulus at every
other temperature is recovered as K_S = γρc_p/α_V. That inversion puts the entire temperature
variation of α_V/c_p into K. For copper the printed 3α/c_p rises 12 % between 10 K and 293 K, so
the model's bulk modulus falls 12 % (implied K_T 163.6 → 142.1 GPa), whereas real copper's K
falls about 3.5 %. The Wachtman E(T) is nearly flat in the same range (−1.1 % against −7.5 %
in the literature), so with ν = ½ − E/(6K_T) the over-stiff K at low temperature drives ν toward
½ at 0 K. Codex confirmed the mechanism independently and sharpened the wording: it is a
one-anchor thermoelastic closure with no independent low-temperature modulus information, not
"constant γ" alone; the flat E opposes the sign rather than causing it.

Literature check (Hill average of the single-crystal constants): ν = 0.338 at low T, 0.345 at
300 K; K 142.0 → 137.1 GPa (−3.5 %), G 51.5 → 47.3 GPa (−8 %), E 137.8 → 127.4 GPa (−7.5 %).
ν rises because the shear modulus softens faster than the bulk modulus.

## Key Findings

- **F1 — unphysical trend** (Claude + Codex, high on direction). BELFEM 0.363 → 0.344 vs
  literature 0.338 → 0.345. Ledbetter 1981 (polycrystalline Cu, 4–295 K, tabulates K, G, E, ν) is
  the direct discriminator; its table has not been read in this session.
- **F2 — mechanism** (Claude + Codex, high). `cl_Material_Metal.cpp:1172` (γ from the anchor),
  `:1197` (K_S = γρc_p/α_V), `:1200` (ν = ½ − E/6K_T). Below the split temperature T* the
  expansion branch enforces dln(α/c_p)/dT > 0 (`cl_Material_SplineLookupTable.cpp:335`), so with
  constant γ the model's K is forced to soften with T by construction.
- **F3 — Wachtman flatness** (Claude, fact high; role: Codex says it opposes the sign). T₀ = 922 K
  at `cl_Material_Copper.cpp:64` suppresses all softening below ~300 K. Already an open decision
  in `dl20260824_wachtman_elastic_review.md:36-39`. The "Blanke starts near 93 K" statement there
  was never checked against Blanke.
- **F4 — adiabatic anchor treated as isothermal** (both, P2). Magnitude only; shifts γ (printed
  2.13 vs literature ≈ 2.0), does not flip the sign.
- **F5 — all nine `create_mech` metals share the closure** (both). Per-metal sign not probed.
- **Codex P2:** `nu(1)`, `nu(2)` read without an active size check (`:1205-1210`); the E spline
  grid is dT = 4 K to T_max, hundreds of points for every built-in metal — latent only.
- **Grok:** wrapper failed with "Not signed in" before any API call; `grok login --device-code`
  restores the leg. One-voice round; Codex agreement is the weakest evidence tier anyway.

## Christian's design constraints (stated 2026-09-15)

- ν must be consistent with the thermophysical properties (c_p, α, ρ) and the thermal-expansion
  data — that was the original intent of the Grüneisen construction.
- For some metals ν is known at a few temperatures, never over the whole range.
- Known structure: dν/dT = 0 at 0 K, and ν increases monotonically with T.
- The metal properties are built from **polycrystalline** (isotropic) data; this should be
  stated in the module documentation (`materials_usage_guide.md` §9.7 and
  `material_property_sources.md` are the natural places). Not done in this read-only session.

## Proposal (not implemented — for Christian's adjudication)

Keep the coupling to thermophysics, but invert it the stable way round. In the quasi-harmonic
picture an elastic modulus depends on volume, not on temperature directly, so

    M(T) = M₀ · exp( −δ_M · ε_V(T) ),   ε_V(T) = 3 ∫₀ᵀ α dT   (volumetric thermal strain)

for M ∈ {K, G}, with one Anderson–Grüneisen-type constant δ_M each (the isothermal δ_T for K is
tabulated for many elements; for Cu the Hill numbers above give δ_K ≈ 3.4, δ_G ≈ 8.0 with
ε_V(293 K) ≈ 1.0 %). E and ν follow from K and G. Properties:

- **Consistency with expansion data by construction** — ε_V is the integral BELFEM already forms
  for ρ(T), so the moduli inherit the fitted α curve and its Debye crossover; no fitted T₀.
- **dν/dT = 0 at 0 K structurally** (dε_V/dT = α_V → 0), and E flat at low T as it should be.
- **Monotone ν(T)** exactly when δ_G > δ_K, i.e. shear softens faster than compression, which is
  the normal case for metals and can be checked per material from the two constants.
- **Wachtman is the special case**: Anderson 1966 derives E₀ − bT·exp(−T₀/T) as an
  approximation to the thermal-energy integral; using ε_V (≈ γρE_th/K) directly removes the
  approximation and the ill-constrained T₀.
- **Data requirements**: two constants per metal. With ν at several temperatures (Cu, Ledbetter
  1981) fit δ_K, δ_G by least squares. With only room-temperature (E, ν) and an E(T) curve
  (Blanke), δ_E comes from the E curve and one of δ_K or δ_G needs a literature value (δ_T or
  the pressure derivative K′ for the bulk modulus). Which fallback to accept is a physics
  decision.
- **γ stays as a diagnostic**, computed from the new K(T) as α_V K_S/(ρc_p), so a wrong δ shows
  up as an unphysical γ(T) — the same role γ already plays for α.

Alternative considered: keep the current closure but replace constant γ by a γ(T) model. Rejected
as a recommendation because γ(T) is exactly the poorly known quantity, and the inversion K = γρc_p/α_V
amplifies its error into K one-to-one.

## Plan and plan audit (same day)

Christian asked for a plan+audit and delivered eleven primary sources to `tmp/papers/` (DOIs in the
plan, O1). Plan: `todo/elastic_moduli_quasiharmonic.md`, registered in `todo/README.md`. Blind jury
(Codex gpt-5.6-terra/high, Grok grok-4.6/high after `grok login`) on rev. 1 returned 20 findings; all
citations re-opened, all confirmed, folded into rev. 2 — the load-bearing ones:

- **Both:** rev. 1 mixed two strain definitions, 3∫α dT (logarithmic) and (l/l₀)³−1 (engineering).
  Adopted the logarithmic strain ε = 3 ln(l(T)/l(0)), which is Garai & Laugier 2007 Eq. (10)
  (10.1063/1.2424535) with constant δ and what the pilot fit used.
- **Codex:** the adiabatic/isothermal question (O2) blocks implementation, because `grueneisen( T )`
  and the `cp_from_debye` dilation term label the served K as isothermal; for copper the choice is
  worth 2.9 % on K and 0.0045 on ν at 293 K — comparable to ν's whole cryogenic change (0.0070).
- **Grok:** `Metal::nu_custom` does not exist (new work, and `set_have( nu )` must precede the spline);
  copper's room-temperature E drops from 133.0 to 128.2 GPa (Ledbetter 1981, 10.1002/pssa.2210660209)
  — a user-visible data correction, not just a slope change; Nickel's test window is 4–293 K.

**Pilot fit (probe):** Ledbetter 1981 Table 2 (16 points, 5–295 K) against ε from BELFEM's own α:
δ_K = 3.46 (rms 0.04 %), δ_G = 8.67 (rms 0.23 %); the derived ν matches the measured ν within 0.0005
over the whole range (0.3406 vs 0.3401 at 5 K; 0.3475 vs 0.3471 at 295 K). To be refit on the code's
spline-integrated ε in R1.

Measured copper, 5 → 295 K (Ledbetter 1981, 10.1002/pssa.2210660209, Table 2): B 144.46 → 139.74 GPa
(−3.3 %), G 51.72 → 47.57 (−8.0 %), E 138.62 → 128.17 (−7.5 %), ν 0.3401 → 0.3471.

**Rev. 3 (same evening).** A second paper batch (`tmp/papers2/`) resolved O1: indium via Kim &
Ledbetter 1998 (10.1016/S0921-5093(98)00490-0, Varshni fits 5–300 K: ν 0.4306 → 0.4498), chromium via
Palmer & Lee 1971 (10.1080/14786437108227390, 4.2–345 K). Christian resolved O2 by rule — serve what
mechanical elastic calculations need — which for BELFEM's quasi-static mechanics means the isothermal
moduli: G is identical, K_T = K_S/(1 + α_V² T K_S/(ρ c_p)) is applied at serve time from the material's
own α, ρ, c_p, and the two existing consumers that treat K as isothermal need no edit. Copper at 293 K
then serves K_T = 135.7 GPa, E = 127.8 GPa, ν = 0.343 (from Ledbetter 1981, 10.1002/pssa.2210660209).
Nadal & Le Poac 2003 (10.1063/1.1539913) records what the quasi-harmonic form does not do: G → 0 at
melting; accepted as out of BELFEM's use.

**Rev. 4 — Nickel (same evening).** Christian asked whether Nickel should keep its own E(T) with the
ΔE dip at the Curie point or take the smooth model for a cleaner Newton iteration. Narrow-claim jury
brief dispatched (`tmp/ai_exchange/review_nickel_elastic.md`). My argument, which he accepted: the dip is
the ΔE effect of the **demagnetized** state and vanishes at saturation (Ledbetter & Reed 1973 §15,
10.1063/1.3253127); BELFEM's nickel sits in tesla fields, so Blanke's curve describes the wrong state;
the anchor itself is state-dependent by ~18 % (190.5 vs 225.6 GPa at 6.2 kOe, Giebe & Blechschmidt 1931
in L&R Table 9); saturated single-crystal constants (Alers, Neighbours & Sato 1960, 10 kOe, L&R Table 6)
Hill-average to E ≈ 223 GPa, ν ≈ 0.30 vs the tree's 183.7 / 0.31. Ruling: Nickel is treated like the other
eight metals; `create_young`, the E override and the Bézier bridge go; the rationale is documented in the
Nickel header (block drafted in the plan, §6a). O3 closed; only O4 remains to confirm. The Nickel jury returned after the ruling and concurred 3/3;
its refinements (the coded curve is already −18 % at 400 K; the T-series lives in Alers, Neighbours &
Sato 1960, 10.1016/0022-3697(60)90125-6, not in Ledbetter & Reed's Table 6; +21 % room-temperature E
change; α clamped 600–631 K; no elasticity weak form is wired yet) are in the plan.

**Rev. 6 — round 2 at xhigh (same evening).** Codex (terra/xhigh) and Grok (4.6/xhigh) endorsed the
architecture and sent rev. 5 back on two P0s with one root: I had placed the adiabatic → isothermal
conversion at serve time. That (a) cited Garai & Laugier's isothermal law (10.1063/1.2424535) for an
adiabatic fit, (b) divided 0/0 at T = 0 where the spline samples, (c) would have leaked c_p features
(Nickel's magnetic c_p rise) into E and ν, and (d) made the scope, pilot and DoD quote adiabatic numbers as
the served result. Fix: convert the measured K_S to K_T in R1 with the exact identity and BELFEM's own α,
ρ, c_p, then fit and serve pure exponentials. Probe on Ledbetter 1981 (10.1002/pssa.2210660209): the
isothermal copper bulk modulus is still an exponential in ε to 0.10 % rms, δ_K = 6.52 (adiabatic 3.46),
K₀ = 144.26 GPa; served copper at 295 K: K 135.7, E 127.8 GPa, ν 0.343. Also fixed: ∂E/∂K (9G², not
27G²) and the discovery that `dEdT_custom` is never called (the spline gets a literal 0.0) so it is
removed; the monotonicity window (whole served range in construction and test); a universal K-ratio bound
that Indium's own data violate (B(0)/B(295) = 1.110, 10.1016/S0921-5093(98)00490-0); numeric R4 oracles;
Chromium's α clamp (1700 vs 2180 K) and anchor away from the Néel point; several stale rows. Whether a
round 3 runs on rev. 6 is Christian's cost call.

**R1 executed (same evening, Christian: "can you do this in your sandbox?").** Ran the existing debug
`material` binary (no build) for the nine metals at 1 K steps, transcribed nine datasets with DOIs
(silver and lead from page images, chromium from Palmer & Lee's 5 K table with C₁₂ = C₁₁ − 2C′), Hill-averaged
(cubic; tetragonal VRH for tin), converted K_S → K_T with the tree's own α, c_p, ρ, and fitted ln K_T and
ln G against the tree's ε. All nine fit to ≤ 0.5 % rms; isothermal δ_K is 6.1–7.9 for the well-behaved metals,
δ_G 7.8–21.5. Scripts and outputs in `tmp/r1_elastic/`. Current-tree probe: the constant-γ closure inverts ν
for seven of nine metals (not for In and Ni). Three data-driven questions filed for Christian: **O5** iron's
isothermal ν is flat (δ_K 7.91 vs δ_G 7.84) so the strict monotonicity guard must become a tolerance;
**O6** chromium's spin-density-wave anomalies (K 192 → 151 GPa toward the Néel point, Palmer & Lee 1971,
10.1080/14786437108227390) defeat the model — recommend constant ν = 0.237 with δ_G = δ_K; **O7** lead's
Hill average gives E = 24.0 GPa against the tree's 16.5 (anisotropy 4.1; Kim & Ledbetter saw the polycrystal
8 % below VRH for indium) — recommend a polycrystalline lead dataset.

**O5–O7 ruled, R2 implemented (late evening).** Christian: "Agreed on all three, go ahead with R2." Source
edits (first of the session, on approval): `Metal::create_mech( E2, nu2, T2, deltaK, deltaG )` rewritten with the
quasi-harmonic K_T, G on ε = 3 ln(l/l₀), isothermal anchor, `BELFEM_ERROR` preconditions, `set_have` + `create_spline`
for E and ν, sampled guards over the whole served range (α ≥ 0, ν ∈ (−1, ½), ν drop ≤ 1e-3); named scalar members
replace the `Vector` bundle; `Metal::dEdT_custom` removed (uncalled); new `Metal::eps`, `nu_custom`. Nickel's
`create_young`/`E_custom`/`dEdT_custom`/`mYoungData`/`mYoungBezier` removed, header block added. Nine constructors
carry the R1 constants with source and DOI; Chromium serves a constant ν = 0.2371; Lead is marked INTERIM. Al and Cr
headers updated. 14 files. Syntax-checked with the debug tree's exact flags on every materials source (clean);
guard logic probed on the R1 tables (all nine pass); code jury (Codex terra/high, Grok 4.6/high) dispatched on the
frozen diff `tmp/ai_exchange/review_elastic_r2.diff`. **Not built** — reviewed, not verified; R3 (probe) and R4
(test) next.

**Code jury, build, gates (night).** The R2 jury (Codex terra/high, Grok 4.6/high) returned P1/P2 only: the
report-constant recovery of K_S was linearized (fixed: exact K_T/(1−q)); "served E(T₂) = E₂" was an overclaim
off-grid (reworded); the guard checked only the analytic curve at knots (now also the served spline between knots);
no test (written); docs stale (rewritten); Chromium's constructor still cited Armstrong & Brown above the Palmer & Lee
call (removed); Copper/Silver file headers stale (fixed with DOIs); Nickel "absorbed by the fit" vs the 0–300 K window
(reworded); Lead's warning only in the .cpp (class header added). Then, on Christian's request, built in
`cmake-build-debug` (single core, `make material test_physics`, clean); `test_physics` 14/14 green including the four
new `MetalElastic` tests (monotone ν for eight metals, copper pinned to Ledbetter 1981 converted, chromium constant ν,
`Sn60Pb40` bounded); R3 probe: all nine construct, ν non-decreasing everywhere, copper's γ now 2.03. High-temperature
deviations from the retired Wachtman curves are tabulated in the plan (Cu +34 % at 1358 K, Cr −49 % at 2180 K — outside
BELFEM's mechanical use, documented, not gated). R5 docs rewritten (usage guide §9.7, property sources with DOIs,
thermal-expansion §5, literature references Tier 6); Codex language sweep (luna/medium) applied, including its flag that
Garai & Laugier's result is for the bulk modulus and the extension to G is BELFEM's own.

**Blanke's digitized curves and the ceilings (late night).** Christian delivered his Wachtman fitting tool and the
digitized Blanke 1989 E(T) curves (`tmp/blanke_wachtman/`), which start at 75–105 K — so the retired fits had a
cryogenic anchor after all, and the earlier "unconstrained below room temperature" statement was wrong; what is true
is that Blanke's copper softens 4.3 % between 100 and 284 K where Ledbetter's polycrystal softens 6.1 %. In the
overlap below 300 K Blanke lies 2–4 % *above* the ultrasonic literature for Cu, Al, Ag (so it is not a static
compilation either; static moduli lie below dynamic ones), within 1 % for Fe and Cr, 33 % below the dynamic Hill
value for Pb (a static-type 16 GPa), 5–12 % above the tetragonal Hill for Sn. Ruling (Christian): measured sources with
data points set the levels; Blanke, uncited and without data points, sets only the validity ceiling, at the
temperature where the quasi-harmonic E(T) departs from Blanke's curve by more than 5 % in shape (both normalized at
room temperature): **Cu 1358 → 1000 K, Ag 1235 → 900 K, Cr 2180 → 570 K, Sn 505 → 400 K**; Fe (855 K) and Al (never
beyond 10 %) unchanged. Rebuilt: `test_physics` 14/14, `make check` 17/17; the `material` tool clamps its tables at
the new ceilings. Also tested and not adopted: (i) Christian's alternative of keeping the Wachtman E with an
optimized γ(T) Bézier — the γ*(T) that Ledbetter's copper implies is smooth (1.89–2.03), so the design is feasible,
but with the existing Wachtman E even a perfect γ inverts ν again (0.345 → 0.337), because that E is 3 % low at 5 K
and 4 % high at 295 K; a refitted three-parameter Wachtman cannot serve 5–1300 K better than ±1 %/15 % either;
(ii) a joint fit of δ_K, δ_G with Blanke's high-temperature points — improves the high range to 4–7 % rms at the
cost of 1.4–2.4 % rms on the cryogenic G, and fails for iron (Blanke's curve carries the approach to the Curie point,
δ_G → 15.6). Static versus dynamic: BELFEM has no mechanical module yet; Christian holds that static moduli would be
the relevant ones for magnet design; the served moduli are dynamic made isothermal, documented as such.

**Closing check — is the resulting γ(T) plausible? (Christian's condition for resolving the session.)** γ is not
prescribed anywhere in the new closure; computed from the served K and the tree's α, ρ, c_p it reads (20 K → 293 K):
Cu 1.88 → 2.03 (lit. 1.96–2.00), Al 2.03 → 2.11 (2.1–2.2), Ag 2.30 → 2.43 (2.3–2.5), In 2.26 → 2.20 (≈ 2.4),
Pb 2.42 → 2.74 (2.7–2.8), Sn 1.88 → 2.21 (2.1–2.3), Fe 1.46 → 1.74 (1.6–1.7), Ni 1.65 → 1.78 (1.8–1.9), Cr 1.07 → 1.04
(low, spin-density-wave metal). A gentle monotone rise to the room-temperature literature value for all metals
except indium (flat, 8 % low) and chromium (flat, as its constant-ν construction implies). γ depends on K_S, α, ρ
and c_p only, so the lead level question (O7) does not touch it. The `material` report now prints γ at 298.15 K,
the temperature of the printed density, labelled as a diagnostic (Christian's request); the constant is evaluated
there in `create_mech` rather than at the anchor. Final rebuild: `test_physics` 14/14, `make check` 17/17; the report
reads γ @ 298.15 K = 2.0346 (Cu), 2.7407 (Pb), 1.0407 (Cr). Session resolved (Christian's condition met).

**Lead, final (Christian: "the values closer to the handbook are more reasonable for our purposes").** Lead is the
one metal served at the static level: bulk modulus and δ_K from the Waldorf & Alers crystal (10.1063/1.1931149; a bulk
modulus is not relaxed by anelasticity), E level and shape from Blanke's static-type curve 95–300 K, which with that K
gives E 16.25, G 5.66 GPa, ν 0.435 at 300 K (handbook 16 / 5.6 / 0.44; Köster & Franz 1961 Table IX 0.45), δ_G = 11.20,
ν monotone from 0.43 at 95 K. Constants (16.25, 0.4347, 300, 7.31, 11.20). O7 closed; the "interim" header replaced by
the static-level rationale; `material_property_sources.md` and the usage guide updated. The γ report label lost its
"( diagnostic )" suffix at Christian's request; the constant is still the diagnostic value at 298.15 K.

## Changes Made / Proposed

- `src/physics/materials/main.cpp` — γ report line at 298.15 K, labelled diagnostic
- `src/physics/materials/cl_Material_{Copper,Silver,Chromium,WhiteTin}.cpp` — `T_max` ceilings with reason
- `tests/physics/test_MetalElastic.cpp` (new) + `tests/physics/CMakeLists.txt`
- `src/physics/materials/doc/{materials_usage_guide,material_property_sources,thermal_expansion_from_heat_capacity}.md`, `doc/literature_references.md`
- `src/physics/materials/cl_Material_Lead.hpp` — interim-data header
- `src/physics/materials/cl_Material_Metal.{hpp,cpp}` — new closure, guards, members
- `src/physics/materials/cl_Material_Nickel.{hpp,cpp}` — ΔE machinery removed, rationale header
- `src/physics/materials/cl_Material_{Aluminum,Chromium}.hpp` — header comments
- `src/physics/materials/cl_Material_{Aluminum,Chromium,Copper,Silver,Indium,Lead,WhiteTin,Iron}.cpp` — constants + DOIs

## Changes Made / Proposed (earlier: none)

None to source. Exchange thread: `tmp/ai_exchange/review_copper_poisson.md` (pre-registration,
Codex audit, verification, reconciliation table); brief: `review_copper_poisson_brief.md`.

## Open Questions

- Christian commits (source, test, docs, plan, devlog); then the plan moves to `todo/closed/`.
- None on the elastic closure. Christian commits.
- Lead: find a polycrystalline ultrasonic dataset (O7 interim in place).
- Check per metal whether the current ν(T) is inverted (probe: print ν at 4/77/293 K).
- Document the polycrystalline/isotropic basis of the metal roster.
- Re-run the Grok leg after `grok login --device-code` if a second voice is wanted.

## Files Updated

- devlog/dl20260915_copper_poisson_trend.md (new)
- devlog/README.md (index line)
- todo/elastic_moduli_quasiharmonic.md (new, rev. 2)
- todo/README.md (registration)
