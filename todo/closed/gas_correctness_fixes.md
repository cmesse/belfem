# Gas Modules: Physical-Correctness Fixes from the Literature Verification

**Date:** 2026-08-06
**Purpose:** Checklist of the defects found by the two-round literature verification of
`src/physics/{gastables,gasmodels}` against the primary sources (devlog
`dl20260806_gas_literature_verification.md`). Every constant and formula named here was
verified against the paper by page-image or text extraction; confidence high unless stated.
**Module:** `../../src/physics/gasmodels`, `../../src/physics/gastables`
**AIs involved:** Claude (five + three verification subagents, adjudication)
**Status:** FIXES APPLIED 2026-08-06 (Claude, on Christian's go-ahead) — D1-D12, H1-H3, O2
comment, O3 rows landed; D2b left open (Christian has not ruled on the MLSG region).
Numerically validated (see devlog): p(Tc,rhoc) = 4.599200 MPa, w = 205.95 m/s, NBP
restored, lambda matches NIST far-field, enhancement diverges at Tc. Codex + Grok jury
review pending; build/test gate is Christian's.
**2026-08-09: D16 added and fixed** — a sign typo in methane residual term 34, found from the
`make check` abort in `Helmholtz::v()`; it broke the critical region only, which is why the
D1 campaign's far-field validations all passed.
**Re-verified 2026-08-09 (currentness sweep): accurate as written; two open boxes remain,
both awaiting Christian, neither a code defect.** D2b (the MLSG scaled-EoS window) and D15
(the deuterium reference state) are both scope/physics decisions. The two gates named in the
status line — the Codex+Grok jury review and Christian's build/test run — have not been
recorded anywhere since, so treat them as still outstanding. Sibling files:
`../gasmodels_open_source_migration.md` (the move itself) and `../nitrogen_eos_completion.md`
(explicitly out of scope here).

> **Scope guards:**
> - `cl_GM_EoS_Nitrogen.*` is WIP and excluded.
> - The departure-at-gPref convention, k_ij support, and the design questions from
>   `dl20260805_gas_modules_jury_audit.md` are NOT re-opened here.
> - Everything verified clean is listed at the bottom to avoid re-audit.

---

## A — Live defects (wrong numbers in reachable paths)

- [x] **D1 (CRITICAL) — Methane Gaussian centers swapped.**
      `cl_GM_EoS_Methane.cpp:108-114` (+ `.hpp:215-220`): code uses ψ = {1.07, 1.11, 1.11, 1.11}
      as the δ-center and γ = 1 as the τ-center. Setzmann & Wagner Tables 35/36: Δᵢ = 1.0
      pairs with δ, γᵢ = 1.07/1.11 with τ. **Fix: swap the roles** (mPsi ↔ mGamma contents, or
      swap their use in `update_g`). Verified consequence of the swap: p(Tc,ρc) = 4.567 MPa
      instead of 4.599200, ∂p/∂ρ < 0 at the critical point (mechanically unstable), cv/cp/w
      off ~2× on the critical isochore at 195 K; exact far from Tc, hence invisible at the
      1 MPa test sweep. After fixing, enable `DISABLED_Methane_Cp_FD_NearCritical`.
      **Fixed 2026-08-06 (Claude): mPsi/mGamma contents swapped to the paper pairing; verified p(Tc,rhoc) = 4.599200 MPa, w = 205.95 m/s via numerical mirror.**

- [x] **D2 (CRITICAL) — Methane λ_cr missing both critical factors.**
      `cl_GM_HelmholtzTransport_Methane.cpp:260-267`: Friend 1989 Eq. (18) is
      λ_cr = [Λ\*·k·ρc²/(6π·η·pc)] · [(T/ρ)(∂p/∂T)_ρ]² · **χ_T\*^((γ−ν)/γ) · F(T\*,ρ\*)** —
      the code computes only the first two factors. `chi()` (:304-313) and `f()` (:317-326)
      are correctly transcribed (Eqs. 19a / 20, F_T = 2.646, F_ρ = 2.678, F_A = −0.637) and
      `mConstExpLambdaCr` = (γ−ν)/γ = 0.468 (:178) is correct — all three are **dead code,
      never called**. Measured consequence: total λ +60 % at 300 K/1 atm, +90 % at 600 K,
      +7 % NBP liquid, and no critical enhancement at Tc. **Fix (minimum):** multiply the
      `lambda_cr()` return by `pow( max(chi(),0), mConstExpLambdaCr ) * f()`; compute the
      signed T\* = 1 − 1/τ once in `lambda_cr()` (the `f()`-internal |T\*| is fine for F
      itself; the MLSG branch needs the signed value). Validation: λ_cr must vanish far from
      critical (restoring the already-verified λ0+λex NIST agreement) and peak near (Tc, ρc).
      **Fixed 2026-08-06 (Claude): chi^0.468 and f() wired into lambda_cr with a chi <= 0 early-out; far-field lambda back to NIST (0.0344 at 300 K/1 atm), enhancement diverges toward Tc.**

- [ ] **D2b (LOW, after D2) — MLSG scaled-EoS region for χ_T\* not implemented.**
      Friend Eqs. (23)-(26): inside |T\*| < 0.03 AND |ρ\*| < 0.25 (185-196 K,
      7.6-12.7 mol/dm³) χ_T\* must come from the scaled EoS (Table 10: Q = 0.1133,
      S = −6.098, W = −1.401, Γ = 0.0801, a = 3.352, b = 0.732, E = 0.287, R = 0.535,
      β = 0.355; Eq. 26 Γ·|T\*|^(−γ) on the critical isochore), else the enhancement carries
      the mean-field divergence instead of the fitted γ = 1.190. Only matters in that small
      window; decide whether BELFEM needs it (cryo methane use sits below, combustion above).
      Full fix spec is in the devlog's round-2 section.

- [x] **D3 (HIGH) — Methane pc repeats the paper's own front-matter typo.**
      `cl_GM_EoS_Methane.cpp:47`: 4.5922e6. Setzmann & Wagner print **both** values:
      4.5922 in the front-matter summary and §2.2 "(4.5922 ± 0.002)", 4.5992 in the vapor-
      pressure block (Eq. 3.2) and in the Eq. (5.3) constraint statement ("constrained to the
      critical parameters given in Eq. (2.3): pc = 4.5992 MPa"). Adjudicated 4.5992 correct:
      the EoS itself evaluates to p(Tc,ρc) = 4.599200 MPa, the vapor-pressure ancillary
      reproduces NBP 101 325 Pa only with 4.5992, Friend Table 1 and NIST agree. **Fix:**
      `mPcrit = 4.5992e6` (the transport header comment `.hpp:46` already says so).
      Consequence of leaving it: p_vap −0.15 % everywhere, wrong p_crit for all data-object
      consumers (SRK helper EoS, Z_crit, transport reducing constants).
      **Fixed 2026-08-06 (Claude): mPcrit = 4.5992e6, adjudication comment added in code; NBP p_vap = 101 324 Pa restored.**

- [x] **D4 (HIGH, tied to D3) — Λ\* rescale algebra.**
      `cl_GM_HelmholtzTransport_Methane.cpp:172-175`: coded 2.231293e9 = 2.235e9 ×
      (Tc_F/Tc_SW)² × (pc_SW/pc_F), but Tc cancels in λ_cr (Tc²/τ² = T²) — the correct
      compensation is Λ\*′ = 2.235e9 × (pc_EoS/4.5992e6) × (162.660/ρc_EoS)². **With D3
      fixed and ρc identical, Λ\*′ = 2.235e9 exactly — drop the rescale.** (Standalone effect
      only 0.014 %; this is about defensible algebra, and it disappears into D3.)
      **Fixed 2026-08-06 (Claude): rescale dropped, Lambda* = 2.235e9 with a comment deriving why it applies unchanged.**

- [x] **D5 (HIGH) — MC-SRK generalized c1 order-of-magnitude typo.**
      `cl_GM_EoS_AlphaFunctionFactory.cpp:427`: 0.16054 → **1.60539** (Coquelet 2004
      Eq. 14, verified against the HAL preprint hal-04180435: c1 = 1.60539ω − 0.10935ω²
      + 0.51780; the code's −0.1094 and 0.5178 are fine). Live for any species outside the
      22 hardcoded CAS entries; water-class ω currently gets 0.560, but ~1.06 is correct.
      **Fixed 2026-08-06 (Claude): 1.6054 (4-decimal rounding matching the sibling coefficients).**

- [x] **D6 (MEDIUM) — Chung β ω²-coefficient transcription error.**
      `fn_GT_idgas_lambda.cpp:46`: 1.368 → **1.3168** (Poling 5th ed. Eq. 10-3.14, five
      consistent occurrences incl. two worked examples). Effect ≤ ~0.4 % on λ (H₂O-class ω);
      applies only to species with critical data but no trans.inp record. While there:
      replace the unidentified "(104 b)"-style comment citations with Poling equation
      numbers (the unnamed secondary source is probably the typo origin).
      **Fixed 2026-08-06 (Claude): 1.3168; all comment citations in idgas_lambda/idgas_mu now name Poling equation numbers.**

- [x] **D7 (LOW live effect) — Oxygen `compute_phi0_tt` power slip.**
      `cl_GM_EoS_Oxygen.cpp:148`: `6·mK(1)/τ³` → must be `6·mK(1)/τ⁴` (Schmidt & Wagner
      Eq. 15: the k₂ term is k₂τ⁻²; the first derivative at :133 is already correct).
      In-range effect up to ~0.03 % on ideal-gas cv at 300 K.
      **Fixed 2026-08-06 (Claude): denominator now tau^2 * tau^2.**

- [x] **D16 (CRITICAL) — Methane residual term 34 has the wrong sign.**
      `cl_GM_EoS_Methane.cpp:90`: `mN(33)` = `+1.387292044e-2`; Setzmann & Wagner Table 35,
      term 34 (d = 4, t = 22, c = 4) is **−0.1387292044·10⁻¹**. The term is the only one of
      the 40 that disagrees with the published set. Consequence: Eq. (5.3) violates its own
      critical constraints — ∂p/∂ρ = −16129 Pa·m³/kg and ∂²p/∂ρ² = −496 at (Tc, ρc) instead
      of zero — so every isotherm up to ≈ 197 K carries a spurious van der Waals loop
      (p ranges over 3.94…4.82 MPa across 135 < ρ < 206 kg/m³ at Tc). Symptom:
      `GASMODELS.Methane_Cp_FD_NearCritical` aborts in `Helmholtz::v()` with
      "Too many iterations for T = 190.595 K, p = 50.000 bar" — the density Newton has no
      unique root to find. Saturation line is off by −4.5 % in p at 180 K and −16 % at 190 K
      against the paper's own ancillary equations (3.2)/(3.4); exact below ~170 K, hence
      invisible to every other test.
      **Fixed 2026-08-09 (Claude): sign restored. Verified against a numerical mirror:
      ∂p/∂ρ and ∂²p/∂ρ² at (Tc, ρc) now vanish, isotherms above Tc are monotone, and p(T,ρ)
      matches CoolProp's Setzmann & Wagner implementation to 5.7e-6 uniformly — that residual
      is BELFEM's Rm/M = 518.2653 J/(kg·K) against the paper's rounded 518.2705.**

- [x] **D17 (CRITICAL) — the cubic entropy reference: `mSref` was compensation, not a
      double count.** The departure spline is anchored on the component's STANDARD entropy,
      `mDepartureSpline.update_data( ..., gTref, data(g)->Sref() / data(g)->M() )`
      (`cl_GM_EoS_Cubic.cpp:628`), so `sdep0( T )` — an entropy *departure* — carried
      +Σ yᵢ sᵢ°/Mᵢ ( air: 6702 J/kg/K ). `Gas::realgas_s` subtracts `sdep0`, and the retired
      `+ mSref` = `idgas_s( gTref, gPref )` = 6702 + 162 ( mixture term ) cancelled it
      exactly. Commit `1c04933f` removed `mSref` and added the mixture term explicitly, which
      left the 6702 uncancelled: every cubic gas returned an absolute entropy one standard
      state too low. Measured as `GASMODELS.Cubic_Entropy` failing with r2 = −56.4 for both PR
      and SRK — an RMS residual of 6648 J/kg/K against VDI D2.2 Table 7 at 250 bar, which is
      6702 minus the ~50 J/kg/K by which the cubic models genuinely miss that table
      ( CoolProp's air EoS reproduces the same table to +50…65 J/kg/K, so the test's
      construction — VDI values plus `s( 298.15, 1 bar )` — is sound ).
      **Fixed 2026-08-09 (Claude): the departure spline's entropy is now anchored at 0.0, so
      `sdep0` means what its name says and the cleaned-up `realgas_s` is exactly the
      pre-`1c04933f` expression, term for term** ( `spline_entropy(Tref)` = Σ yᵢ sᵢ°/Mᵢ is the
      same mass-weighted sum that `sdep0` carried, so the cancellation was exact, not
      approximate ). The commit message's "13.6 instead of 6.86 kJ/kg/K" was reasoned from the
      code on the assumption that `sdep0` is a pure departure; it was not measured.
      Follow-up worth considering: `DISABLED_Gas_Entropy_RefState` is the natural regression
      gate for this, but at `gPref` the cubic still differs from the ideal gas by
      sdep( gTref, gPref ) ≈ −0.5 J/kg/K, i.e. ~7e-5 relative against its 1e-4 tolerance —
      tight enough that it should be re-tolerated before being enabled.

## D13 — transport glue search never terminates on a smooth junction (FIXED)

- [x] **D13 (CRITICAL, runtime abort) — `create_glue_polys_transport` had no
      smooth-junction guard.** Reported by Christian as a regression: the `gas`
      executable aborts in `RefGas::create_glue_polys_transport`
      (`cl_GT_RefGas.cpp:544`, `Assertion tDeltaT <= tDeltaTmax failed`) while
      building the default air mixture.
      **Not caused by the 2026-08-06 work** — `cl_GT_RefGas.cpp` changed that day
      only through the `aT`->`T` parameter rename, and the two functions that did
      change (`fn_GT_idgas_mu`, `fn_GT_idgas_lambda`) feed only *synthesized*
      transport, whereas this loop glues *tabulated* intervals. It is the
      transport twin of the heat-glue defect fixed on 2026-08-05: the search
      widens its window hunting a curvature sign change, the condition
      `( min < 0 && max < 0 ) || ( min > 0 && max > 0 )` excludes zero by
      construction, and across an already smooth junction there is nothing to
      find, so dT runs to the 50 K cap and the BELFEM_ERROR fires.
      `create_glue_polys_heat` received a smooth-junction skip that day;
      `create_glue_polys_transport` did not.
      **Why it surfaces now:** the rebuilt `trans.inp` merges fitted low
      temperature intervals that are constrained to meet the tabulated data in
      value and slope. 11 of the 14 species in the default air composition carry
      one (all but N2O, NO2, O3), so the very first species trips it.
      **First fix 2026-08-06 (Claude): INSUFFICIENT.** Ported the heat guard
      (1e-7 relative smoothness skip). Christian's rebuild still aborted - same
      assert, now reached THROUGH the guard, because the regression has a second
      mechanism the heat case never had: at some junctions the two branches
      genuinely curve in **opposite directions** (an inflection of ln lambda near
      the junction), and a glue whose second derivative never changes sign cannot
      exist there at ANY window width - it matches both end curvatures, so a sign
      change is forced. No smoothness threshold can separate these from
      searchable junctions (the aborting CO junction has a SMALLER kink than a
      neighboring junction that succeeds).
      **Second fix 2026-08-06 (Claude, Fable):** two-phase acceptance. Phase 0 is
      the historical zero-sign-change search, unchanged; phase 1 runs only when
      phase 0 exhausts and accepts exactly ONE curvature sign change - the
      inflection the branches demand, still rejecting wiggly glues. Junctions the
      old condition accepted are therefore glued byte-identically. Verified by
      simulating the exact C++ logic over ALL 34 species of trans.inp
      (scratchpad `glue_final.py`): 0 aborts; phase 2 rescues exactly four
      junctions - CO C@162.5 K, CO2 C@250 K, Xe C@201.2 K (the air aborters, CO
      first in processing order matching Christian's trace) and D2 C@106.3 K,
      which would have hit the first fusion-fuel run - all at the minimal
      dT = 5 window. Still needs Christian's rebuild for the executable gate.

## D14 — element references missing from the rebuilt tables (FIXED, build pending)

- [x] **D14 (CRITICAL, runtime abort) — formation-table references absent.** With
      D13 fixed, the `gas` executable reached the next first-run failure of the
      rebuilt tables: `Gas::create_reference_gases` (`cl_Gas.cpp:342`) anchors
      formation enthalpies on the standard state of each element, and for carbon
      that is `C(gr)` — a condensed species the gas-only table rebuild never
      included. Full closure computed from the shipped compositions: the
      reachable missing references were **C(gr)** (any carbon species, hits
      default air), **Cl2** (CH2CL2, CF2CLBr) and the bromine reference
      (BrF3, CF2CLBr).
      **Fixed 2026-08-06 (Claude, Fable):** the three species added to
      `../../scripts/fluidprop/species.txt` and both tables regenerated with the
      README's canonical invocation. Verified: zero pre-existing lines altered
      in either table, exactly the three thermo records (+ Cl2 transport)
      appended. `e-` (for the ions) was already shipped.
      **Also fixed in code:** `Gas::reference_element` mapped bromine to
      `Br2(cr)` — but bromine is LIQUID at 298.15 K, CEA anchors Br species to
      `Br2(L)`, and the crystal record ends below 298.15 K, so the formation
      table would have anchored on an extrapolation. Now `Br2(L)`, with the
      rationale in a comment.
      **Process note for the record:** the first two regeneration runs used
      default flags and silently DROPPED the NIST-sourced low temperature
      transport intervals of CO, D2, Kr, Ne and Xe — the canonical invocation
      in `../../scripts/fluidprop/README.md` requires `--nist-transport`. Caught by
      diffing against the committed table before shipping; the tool is
      reproducible when invoked as documented.

- [ ] **D15 (physics question -> Christian) — deuterium reference state.**
      `element_to_molecule` has no `D` entry, so element D falls through the
      default and the formation reference for deuterated species is monatomic
      `D` gas rather than `D2`. No abort (D is shipped), but formation
      anchoring of D2O/HD/OD in mixtures differs from the D2 standard state by
      the D2 dissociation enthalpy. Decide: add `D -> D2` (and D2 stays the
      shipped reference), or keep as is deliberately.

## B — Latent defects (dormant in the valid range, fix for consistency)

- [x] **D8 — Oxygen φ0 electronic term: sign and omission.** Two coupled edits make all
      three routines consistent with Schmidt & Wagner Eq. (15), term k₆·ln(1 + ⅔e^(−k₈τ)):
      uncomment the mK(5) term in `compute_phi0` (`cl_GM_EoS_Oxygen.cpp:119` — the commented
      line is already correct), and flip the sign at `:137` to `−mK(5)·mK(7)·…` (phi0_tt at
      :151 is already correct). Dormant below mTmax = 300 K (term ≤ 1e-16) but currently
      s/g omit a term whose second derivative feeds cv — internally inconsistent.
      **Fixed 2026-08-06 (Claude): mK(5) term restored in phi0, sign flipped in phi0_t; all three routines now consistent with Eq. (15).**

- [x] **D9 — Oxygen triple point.** `cl_GM_EoS_Oxygen.cpp:43`: 54.33 → **54.361 K**
      (paper p. 181; the abstract's "54 K" is rounding). Only gates the admission window.
      **Fixed 2026-08-06 (Claude): 54.361 K.**

- [x] **D10 — Dead-code digit errors in the classic alpha factories.** Neither path is wired
      (`init_srk`/`init_pr` use only PM/CCR paths) — fix or delete (see O1), don't leave
      landmines:
      - `create_srk` (`cl_GM_EoS_AlphaFunctionFactory.cpp:34`): Soave ω² coefficient
        0.175 → **0.176** (Soave 1972 Eq. 15).
      - `create_pr78` (`:53,59`): 0.374642 → **0.379642**, 1.487503 is correct (Young et al.
        Eq. 7), and the branch threshold 0.4823 → **ω > 0.491** as published (the coded
        threshold is the continuity crossover induced by the wrong constant; the published
        pair has a small κ jump at 0.491 — if fixing, keep the published form and note the
        jump in a comment).
      **Fixed 2026-08-06 (Claude, per O1 keep-and-fix): 0.176; PR-78 on published constants (0.379642, threshold 0.491) with the kappa jump noted in a comment.**

- [x] **D11 — Lucas viscosity: quantum correction F_Q° omitted.** `fn_GT_idgas_mu.cpp`
      implements F_P° only; Poling Eq. 9-4.19 adds F_Q° = 1.22·Q^0.15·{…} for He (Q = 1.38),
      H₂ (0.76), D₂ (0.52). Dormant (He/H₂ carry trans.inp fits) but silently activates if a
      quantum species without transport data is ever added — relevant to the new D/HD/T2-class
      records. Minimum: a comment; better: implement the correction keyed on the three CAS.
      **Fixed 2026-08-06 (Claude): F_Q implemented per Eq. (9-4.19), keyed on the
      He/H2/D2 CAS numbers, 1/M exponent converted to g/mol. Extended same day per
      the jury round: pH2 (1333-74-0p) reuses Q = 0.76 (same molecule); He-3 has no
      published Lucas Q and deliberately falls through to Fq = 1 (commented).**

- [x] **D12 — Chung β applied to polar species.** Poling restricts the ω-correlation to
      nonpolar fluids; polar default β = (1.32)⁻¹ = 0.758 (~3 % λ effect for H₂O-class).
      The dipole is already in `GasData`; switch on it. Same dormant path as D6.
      **Fixed 2026-08-06 (Claude): beta switches to 1/1.32 when the Lucas reduced dipole (Eq. 9-4.17) exceeds the 0.022 threshold of Eq. (9-4.18).**

## C — Hardening (no wrong numbers today)

- [x] **H1 — Cubic liquid-root physicality filter.** `cl_GM_EoS_Cubic.cpp:162-180` picks the
      smallest root with Z ≥ 1e-6 (paper-conformant); a sub-covolume root (v < b) would
      silently make `chi()` return NaN. Filter liquid roots by Z > B instead. No triggering state was
      constructed (low confidence that one exists in practice).
- [x] **H2 — Helmholtz pressure-range guard.** Leachman's validity range ends at 2000 MPa,
      Schmidt & Wagner's at 81.8 MPa; the code guards only T and silently extrapolates in p.
- [x] **H3 — Comment/label corrections (cosmetic, batched):**
      `cl_GM_HelmholtzTransport_Methane.cpp:323` "Eq. (22)" → Eq. (20);
      `.hpp:45-48` claims Friend's Tc/pc, but the members hold EoS values at runtime;
      `chi()`/`f()` doc comments swap "damping"/"crossover" naming (χ is the symmetrized
      compressibility, F the crossover/damping function);
      `cl_GM_Helmholtz.hpp:317-318` "aTau : T / T_crit" contradicts the τ = Tc/T convention
      used everywhere;
      `cl_GM_EoS_AlphaFunctionFactory.cpp:82` "Table V" → Table IV (Coquelet 2004);
      `cl_GM_EoS_Hydrogen.hpp:31` `mNab` doc ("parameter n from Table 4" → a/b term count;
      stored as `real`, compared against `uint`);
      `cl_GM_Helmholtz.hpp:551` assert format string lacks `%u` for its argument.
      **Fixed 2026-08-06 (Claude): liquid root filter is now Z > B with B = p*b/(R*T).**
      **Fixed 2026-08-06 (Claude): mPmax member set to 2000 MPa (H2), 1000 MPa (CH4),
      81.8 MPa (O2). UPGRADED same day on Christian's ruling after the jury round
      (Codex: debug assert missed `Helmholtz::v` and the transport's own `update_Tp`):
      now an always-active BELFEM_ERROR in `Helmholtz::v()`, the choke point of every
      (T,p) path — "if a gas fails, the model fails too". Open remainder: mTmax has
      never been guarded anywhere (pre-existing; a hard T error would change behavior
      for models like O2 with Tmax = 300 K — needs its own ruling).**
      **Fixed 2026-08-06 (Claude): all seven comment/label sites corrected.**

## O — Notes / decisions (no action unless Christian says so)

- [x] **O1** — Delete vs fix the dead classic factories (D10): fixing keeps Soave/PR-78
      available for benchmarking; deleting follows the D5-precedent from the migration plan
      (dead branches removed). Christian's call.
      **RESOLVED 2026-08-06 -> keep and fix (Christian); executed with D10. Not a default choice.**
- [x] **O2** — R vintage: papers fitted with their own R (H₂ 8.314472, CH₄ 8.31451,
      O₂ 8.31434) vs BELFEM's CODATA-2018 Rm — ≤ 1.5e-5 relative, far below EoS uncertainty.
      At most, note this in code comments; REFPROP-style per-EoS R is NOT recommended here.
      **RESOLVED 2026-08-06 -> comment R (Christian); paper-R notes added at the coefficient tables.**
- [x] **O3** — `mahmoodi2016.pdf` appendix has ~200 components with full PM sets;
      para-hydrogen and helium-3 would be natural cryo additions to `cubicalpha.inp` if ever
      wanted. The C3H6 row is cyclopropane (CAS-consistent, matches the paper's cyclopropane
      row) — the bare formula could be misread as propylene; consider a header note.
      **RESOLVED 2026-08-06 -> add para-hydrogen and helium-3 only, not all ~200 (Christian); rows landed in cubicalpha.inp with a header note (inert until gasdata.inp/thermo.inp carry the species). C3H6 cyclopropane note added to the header too.**
- [x] **O4** — `helium.pdf` (McCarty & Arp) was supplied, but no `EoS_Helium` exists — presumed
      planned work, no defect.
      **RESOLVED 2026-08-06 -> helium is planned, but not now (Christian).**

## Verified clean in this campaign (do not re-audit)

Hydrogen is digit-exact vs Leachman (all three variants, all derivatives, ancillary; the M9
double-pv suspicion was REFUTED analytically). Methane EoS otherwise faithful (40 residual
terms, φ0, and all five derivatives, including the 2026-08-05 phir_tt fix — independently confirmed);
methane η0/ηex/λ0/λex are coefficient-exact vs Friend Tables 8/9, including the two-phase δσ
branch, NIST-validated end to end. Oxygen: all 73 published numbers are exact; phir machinery is correct;
mTmax = 300 K IS the paper's validity limit. Cubic SRK/PR machinery is formula-exact vs
Soave/PR (the code's analytic Ω constants are better than the papers' rounded ones); departure
functions and Jacobian derivatives are exact; mixing rules are exact. Both Coquelet 22-species tables
are digit-exact; CCR-PR generalized Eqs. 26-28 are exact; PM form and parser transform
{2C1, −C2², ⅔C3³} confirmed; `cubicalpha.inp` digit-exact vs the Mahmoodi & Sedigh appendix
(one 6th-digit rounding in a column the parser never reads; CAS numbers are BELFEM-added and
all correct). NASA-9 cp/H/S and derivatives are exact vs RP-1311; CEA transport form and unit
conversions are correct (CEA unit is micropoise → 1e-7); `idgas_mu` is the **Lucas** method
(not Chung), with all constants exact, including the (5/88)·10^(1/6) identity.
