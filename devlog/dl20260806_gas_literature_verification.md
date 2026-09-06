# Gas Modules: Literature Verification Round (physical correctness)

**Date:** 2026-08-06
**Purpose:** Read-only physical-correctness scan of `src/physics/gastables` and
`src/physics/gasmodels` against the source papers Christian supplied in `tmp/eos/`
(session-local PDFs, not part of the repo). Completes the "literature checks" item
left open by the 2026-08-05 jury audit (`dl20260805_gas_modules_jury_audit.md`).
**Module:** physics/gastables, physics/gasmodels
**AIs involved:** Claude (five parallel verification subagents + one direct pass)
**Edit state:** No source files touched. Documentation only (this devlog, index,
tracker tick, memory).

## Method

Five parallel Claude subagents, one per model family, each comparing the code
coefficient-by-coefficient and formula-by-formula against its primary source, with
independent hand differentiation of every derivative chain and (where possible)
numerical reproduction of the papers' own constraints. One direct pass on the oxygen
EoS (its paper is not in `tmp/eos/`; structure verified by hand differentiation
alone). `cl_GM_EoS_Nitrogen.*` excluded (WIP, per Christian).

Literature used: Soave 1972 (`srk.pdf`), Peng & Robinson 1976 (`pr.pdf`),
Mahmoodi & Sedigh 2017 (`cubic.pdf`), Young et al. IECR (`alpha.pdf`),
Leachman et al. 2009 (`hydrogen.pdf`), Setzmann & Wagner 1991 (`methane.pdf`),
Span et al. 2000 (`nitrogen.pdf`, used only for the standard property-relations
table), NASA RP-1311 Parts 1+2 (`cea1.pdf` scanned / `cea2.pdf`), Poling 5th ed.
(`Properties_of_Gases_and_Liquids_P.pdf`), textbook real-gas chapter
(`realgas.pdf`). Additionally the Coquelet–Chapoy–Richon 2004 preprint was fetched
from HAL (hal-04180435) to settle the CCR/MC questions from the primary source.
`helium.pdf` (McCarty & Arp) has no corresponding code — noted as presumed future
work.

## P0 — live wrong numbers, production-exposed

1. **Methane Gaussian centers swapped** (`cl_GM_EoS_Methane.cpp:108-114` with
   `.hpp:215-220`): the code uses ψ = {1.07, 1.11, 1.11, 1.11} as the δ-center and
   γ = 1 as the τ-center; Setzmann & Wagner Tables 35/36 have it the other way
   (Δᵢ = 1.0 pairs with δ, γᵢ = 1.07/1.11 with τ). Numerically proven: with paper
   centers the EoS reproduces its published critical constraints exactly
   (p(Tc,ρc) = 4.599200 MPa, ∂p/∂ρ ≈ 0, w = 205.95 m/s); with the code's swap
   p(Tc,ρc) = 4.567 MPa and **∂p/∂ρ < 0 — mechanically unstable at its own
   critical point** — and on the critical isochore at 195 K cv is 2× too high,
   cp 2× too low, w 123.6 vs 264.1 m/s. Far from critical the Gaussians vanish and
   both forms agree to machine precision, which is why the p = 1 MPa tests never
   saw it. Also destabilizes the v(T,p) Newton near Tc and feeds the transport
   class wrong δ-derivatives. Confidence: high. The existing
   `DISABLED_Methane_Cp_FD_NearCritical` gate probes exactly this region.
2. **Methane λ_cr assembled without susceptibility and crossover**
   (`cl_GM_HelmholtzTransport_Methane.cpp:260-267`): the coded critical
   enhancement is Λ·[1+δ(φr_δ−τφr_δτ)]²/(τ²μ) with **no χ factor and no crossover
   F**; the helpers `chi()` (:304-313) and `f()` (:317-326) exist but have no call
   sites, and the exponent ratio `mConstExpLambdaCr` ≈ 0.468 (:178) is computed and
   never used. Measured: no divergence at Tc (total λ ≈ 0.07 vs NIST peak
   ≳ 0.2 W/m/K) **and** a spurious background at all conditions — total λ +60% at
   300 K/1 atm (0.0548 vs NIST 0.0344), ~+90% at 600 K, +7% NBP liquid — while
   λ0+λex alone match NIST almost exactly. Every λ this class returns is corrupted.
   Confidence: high that the coded form is physically wrong (χ → 0 must kill the
   background at low density and diverge at Tc); the exact intended Friend form
   stays medium pending the Friend 1989 paper (still absent).
3. **Methane pc transposition typo** (`cl_GM_EoS_Methane.cpp:47`): 4.5922e6 vs the
   paper's 4.5992 MPa (stated twice; the transport header's own comment at
   `cl_GM_HelmholtzTransport_Methane.hpp:46` carries the correct value). p_vap
   uniformly −0.15% (NBP 101 170 Pa instead of the ancillary's built-in
   101 324 Pa constraint), T_vap ~0.02 K, and the data-object p_crit is wrong for
   every consumer (SRK helper EoS, Z_crit, transport reducing constants).
   Confidence: high.
4. **MC-SRK generalized c1 order-of-magnitude typo confirmed**
   (`cl_GM_EoS_AlphaFunctionFactory.cpp:427`): code 0.16054, published (Coquelet
   2004 Eq. 14, verified against the HAL preprint) **1.60539**. The prior audit's
   suspicion is settled. Any species outside the 22 hardcoded CAS entries that
   takes this fallback gets roughly half the correct α temperature slope at
   moderate-to-high ω (water-class ω: 0.560 vs ~1.06) → large vapor-pressure and
   departure errors. Confidence: high.

## P1

5. **Chung β ω²-coefficient confirmed a transcription error**
   (`fn_GT_idgas_lambda.cpp:46`): 1.368 vs Poling Eq. 10-3.14's **1.3168** (stated
   five times incl. two worked examples; no 1.368 variant exists — the constant-β
   polar default 0.758 is a different object). Effect ≤ ~0.4% on λ (H2O-class ω);
   fires only for species with critical data but no trans.inp record. The
   unidentified "(104b)/(88)" equation numbers in the comments match neither
   Poling nor RP-1311 — the unnamed secondary source is the likely typo origin.
6. **Oxygen `compute_phi0_tt` power slip** (`cl_GM_EoS_Oxygen.cpp:148`): the τ⁻²
   term's second derivative is 6k₁/τ⁴; code computes 6k₁/τ³ (first derivative at
   :133 is correct). k₁ = −6.65e-5, so the cv/cp error is O(1e-4·R) — negligible
   numerically, but a genuine formula defect. Hand-derived; no paper needed.
7. **Oxygen `compute_phi0_t` electronic-term sign** (`cl_GM_EoS_Oxygen.cpp:137`):
   with E₂ = (2/3)e^(−k₈τ) (hpp:215), d/dτ of k₆·ln(1+E₂) is **−**k₆k₈E₂/(1+E₂);
   code adds it with +. phi0_tt (:151) is consistent with the correct sign, so
   phi0_t disagrees with both neighbors. Latent: below mTmax = 300 K the term is
   ≤ 1e-17; it only wakes above ~1000 K. (The known commented-out mK(5) term in
   phi0 itself (:119) — already routed to Christian — is equally dead below 300 K.)

## P2 — latent / dead-code / dormant-path

8. **Dead-code digit errors in the classic factories** (nothing wires them;
   `init_srk`/`init_pr` use only the PM and CCR paths):
   `create_srk` (:34) has Soave's ω² coefficient as 0.175 vs published 0.176;
   `create_pr78` (:53,59) has 0.374642 vs published 0.379642 and a switch at
   ω = 0.4823 vs published 0.491 — the threshold is exactly the continuity
   crossover of PR-76 with the altered constant, so the deviation is
   self-consistent but non-standard (~0.5% low κ for ω > 0.48). Landmines if ever
   wired.
9. **Chung β applied to polar species** (`fn_GT_idgas_lambda.cpp`): Poling
   restricts the ω-correlation to nonpolar fluids and recommends β = 0.758 as the
   polar default (~3% λ effect for H2O-class). Dormant (fallback-only path); the
   dipole is already in GasData if a switch is ever wanted.
10. **Lucas viscosity omits the quantum correction F_Q°** (`fn_GT_idgas_mu.cpp`;
    Poling Eq. 9-4.19): −14% (H2) / −24% (He) if the fallback ever ran for a
    quantum gas. Currently dormant (He/H2 carry trans.inp fits) — but the new
    D2/HD/T2-class species make this worth a guard or comment.
11. **Liquid-root selection lacks a Z > B physicality filter**
    (`cl_GM_EoS_Cubic.cpp:162-180`): smallest-positive-root per the papers, but a
    sub-covolume root (v < b) would NaN `chi()` silently rather than fail loudly.
    No triggering state constructed; low confidence it can occur in practice.
12. **No pressure ceiling in the Helmholtz range guards**: Leachman validity ends
    at 2000 MPa; code enforces T-range only and silently extrapolates in p.
13. **R vintage**: BELFEM's CODATA-2018 Rm differs from the papers' fitted
    R = 8.314472 (H2) / 8.31451 (CH4) by ~1 ppm — far below EoS uncertainty;
    note only.

## Verified correct (the strong results)

- **Hydrogen vs Leachman 2009: fully clean.** All three variants (para/normal/
  ortho) implemented; every coefficient of Tables 3–8 digit-exact (φ0 a/b pairs
  cross-validated through Table 3's dimensional v_k = −b_k·Tc); all seven
  derivative routines hand-verified, including the Gaussian phir_tt squared
  bracket; vapor ancillary exact. **The suspected double-pv subtraction in
  `set_reference_point` is refuted analytically** — mU0 = mH0 = −h_raw(ref) is the
  unique choice preserving u = h − pv identically; the M9 probe
  (`DISABLED_Helmholtz_u_Consistency`) should pass.
- **Property assembly (`cl_GM_Helmholtz.cpp`)** matches the standard relations
  (Span 2000 Eqs. 56–64) everywhere, including dpdv/dpdT, w², and the Maxwell
  relation in dsdp; v(T,p) Newton judged physically sound.
- **Methane EoS otherwise faithful**: all 40 residual terms, ideal-gas part, and
  all five derivatives (including the 2026-08-05 phir_tt fix, independently
  re-derived and FD-verified — **fix confirmed correct**) match Setzmann & Wagner;
  η, λ0, λex transport pieces reproduce NIST end-to-end (η within ~1% gas and
  liquid).
- **Cubic machinery exact against Soave/PR**: generalized-form mapping, exact
  analytic Ω constants (better than the papers' rounded values), cubic-in-Z
  coefficients, all first/second derivatives, hdep/sdep/cpdep + their p/T
  derivatives incl. the d²a/dT² mixture chain rule, vdW mixing rules, root
  prescriptions, and the departure-spline entropy integral.
- **Alpha functions**: both 22-species Coquelet tables digit-exact vs the primary
  source (incl. catching the paper's own "sulfur dioxide" row that is actually
  H2S — code maps it correctly); CCR-PR generalized Eqs. 26–28 exact; MC c2/c3
  exact; PM form and the parser transform {2C1, −C2², ⅔C3³} confirmed by four
  independent arguments; all four α/α'/α'' chains FD-verified to the noise floor
  on both sides of Tc; branch continuity as published.
- **Gastables layer exact**: NASA-9 cp/H/S and derivatives ≡ RP-1311 Eqs.
  4.9–4.11; CEA transport fit ≡ Eq. 5.1 with the correct **micropoise** → Pa·s
  (1e-7) and µW/cm·K → W/m·K (1e-4) conversions; the viscosity synthesis is the
  **Lucas method** (not Chung as previously assumed) and every Lucas constant is
  algebraically exact in BELFEM's SI units, including the (5/88)·10^(1/6)
  inverse-ξ identity.
- **Oxygen**: derivative structure of all six phir routines matches the 32-term
  Schmidt & Wagner layout (13/11/8 groups, exp(−δ²)/exp(−δ⁴) families) by hand
  differentiation; Tc/pc/ρc match published values.

## Literature still missing (blocks the remainder)

- Friend, Ely & Ingham 1989 (methane transport) — all coefficient tables and the
  exact λ_cr form (P0 #2's intended shape) unverifiable without it.
- Schmidt & Wagner 1985 (O2) — the 32 mN coefficients unverified (structure OK).
- Coquelet 2004 journal version (audit used the HAL preprint) — last-digit risk
  only (cyclohexane C3 1.375 vs preprint's "1.3751").
- Mahmoodi & Sedigh supplementary material — per-species digits of
  `share/fluidprop/cubicalpha.inp` (transform + Tc/pc columns + constraint
  saturation all verified; individual C-digits not).
- McCarty & Arp helium: paper supplied, no code — presumed planned.

## Round 2 (same day): the requested papers arrived

Christian supplied `friend1989.pdf`, `schmidt_wagner_1985.pdf`, and `mahmoodi2016.pdf`
(the M&S paper including its appendix parameter tables). Three more verification
agents closed out everything the morning round left blocked. Tick-off list for the
fix round: **`../todo/closed/gas_correctness_fixes.md`** (registered in `todo/README.md`).

### Methane transport vs Friend 1989 — closed

- Every coefficient verified digit-exact from page renders: collision-integral fit
  (Table 8, incl. the k/3−1 exponent scheme), σ/ε_k (Table 1), F_int pair, all 11
  excess-viscosity g/r/s and all 7 excess-conductivity j/r/s (Table 9), the
  two-phase δσ branch (Eq. 16), and both SI prefactors (each reproduces the paper's
  evaluated mixed-unit constant to ≤0.01 %).
- **The dead helpers are themselves correct**: `chi()` ≡ Eq. 19a (the 0.28631 in the
  paper IS Z_c — using the EoS-consistent value is the right call), `f()` ≡ Eq. 20
  with Table 10's F_T/F_ρ/F_A, and 0.468067 ≡ (γ−ν)/γ. The defect is purely that
  `lambda_cr()` never calls them (P0 #2 confirmed at equation level).
- **λ_cr fix spec** (D2 in the todo): multiply the existing return by
  `pow( max(chi(),0), mConstExpLambdaCr ) * f()`; compute signed T* = 1 − 1/τ once
  in `lambda_cr()` (Eq. 21 is signed; `f()`'s internal |T*| is fine for F itself).
  Optional completeness (D2b): Friend mandates the MLSG scaled EoS for χ inside
  |T*| < 0.03 ∧ |ρ*| < 0.25 — Eq. 23 χ = Q·|ρ*|^(−a)·θ^b/[θ + Ω(θ+R)] with
  Ω = W·T*·|ρ*|^(−1/β) (Eq. 25), θ = 1 + E(1 + S·T*·|ρ*|^(−1/β))^{2β} for
  T* < −|ρ*|^{1/β}/S else 1 (Eq. 24), and χ = Γ·|T*|^(−γ) on the critical isochore
  (Eq. 26); Table 10: Q = 0.1133, S = −6.098, W = −1.401, Γ = 0.0801, a = 3.352,
  b = 0.732, E = 0.287, R = 0.535, β = 0.355, γ = 1.190, ν = 0.633.
- **Λ* rescale** (:172-175): the coded 2.231293e9 = 2.235e9·(Tc_F/Tc_SW)²·(pc_SW/pc_F)
  keeps the SI constant invariant, but Tc cancels in λ_cr — correct compensation is
  pc/ρc² only. With the pc fix below, Λ*′ = 2.235e9 exactly; drop the rescale.

### Methane pc — adjudicated, and the code faithfully copied a paper typo

The two agents "disagreed" because **Setzmann & Wagner print both values**: 4.5922
in the front-matter property summary and §2.2 "(4.5922 ± 0.002) MPa"; 4.5992 in the
vapor-pressure block (Eq. 3.2) and in the Eq. (5.3) constraint statement
("constrained to the critical parameters given in Eq. (2.3): pc = 4.5992 MPa").
Verified directly from the PDF (pages 5/9 vs 27/33). 4.5992 wins on every physical
anchor: the EoS itself evaluates to p(Tc,ρc) = 4.599200 MPa, the ancillary
reproduces NBP 101 325 Pa only with 4.5992, Friend Table 1 and NIST concur. So
`mPcrit = 4.5922e6` is the paper's own front-matter typo, transcribed faithfully.

### Oxygen vs Schmidt & Wagner 1985 — closed, fully clean data

All 73 published numbers digit-exact from page renders: 32 mantissas+exponents of
n, all 32 r/s exponent pairs, all 9 ideal-gas k (term pairing incl. the 2/3
electronic degeneracy factor), Tc/ρc/pc (Weber values = reducing parameters). The
morning round's two hand-derived φ0 findings are **both confirmed against Eq. 15**
(phi0_tt power slip — the only one with in-range effect, ~0.03 % on cv at 300 K;
phi0_t electronic sign flip — dormant), and the commented-out mK(5) term belongs in
φ0 (uncomment + sign fix = D8). New: mTtriple 54.33 vs the paper's 54.361 K, and
one correction to the morning report — **mTmax = 300 K IS the paper's validity
limit** (data ends there; only the ideal-gas part extends to 3000 K), not a
conservative cap.

### cubicalpha.inp vs the M&S appendix — clean

All 26 species digit-exact in all 12 coefficient columns + Tc/pc (extraction was
layout-clean, no ambiguous digits). Sole deviation: H2S PR-PM2 C2 rounded
0.098518 → 0.09852 in a column the parser never reads. CAS numbers are a BELFEM
addition (paper has none) and all 26 are correct — including C3H6 = cyclopropane
(75-19-4), consistent with its coefficients matching the paper's cyclopropane row.
The appendix carries ~200 components; para-hydrogen and helium-3 are ready
transcriptions if the cryo set ever wants them.

## Fix round (same day, approved by Christian)

Christian's rulings: O1 keep-and-fix the classic factories (not a default choice),
O2 comment the R vintage, O3 add para-hydrogen and helium-3 only, O4 helium planned
but not now. All defects fixed except D2b (MLSG region — not ruled on, left open):

- D1 methane Gaussian centers swapped to the paper pairing; D3 pc = 4.5992e6 with
  the adjudication in a code comment; D2 chi^0.468 and f() wired into lambda_cr
  (chi <= 0 early-out for metastable states); D4 rescale dropped, Λ* = 2.235e9
  with the derivation commented.
- D5 c1 = 1.6054; D10 Soave 0.176 + PR-78 published constants with the 0.491
  threshold and the κ-jump note.
- D6 Chung β = 1.3168 + Poling citations throughout idgas_lambda/idgas_mu;
  D12 polar species now get β = 1/1.32 (switch on the Lucas reduced dipole,
  same 0.022 threshold the file already uses); D11 Lucas F_Q quantum correction
  implemented for He/H2/D2 keyed on CAS.
- D7/D8/D9 oxygen: phi0_tt τ⁻⁴ power, electronic term restored in phi0 with the
  sign fixed in phi0_t, T_triple = 54.361 K; paper-R comment (O2).
- H1 cubic liquid root filter Z > B; H2 mPmax debug guard (2000/1000/81.8 MPa);
  H3 all seven comment sites.
- O3: pH2 (CAS 1333-74-0p, CoolProp convention) and He3 (14762-55-1) rows added
  to `cubicalpha.inp` from the M&S appendix — inert until gasdata/thermo carry
  the species; header notes added (incl. C3H6 = cyclopropane).

**Numerical validation** (scratchpad `fixed_code_check.py`, mirroring the fixed
code against the agents' paper reimplementation): p(Tc,ρc) = 4.599200 MPa with
∂p/∂ρ ≈ 0 and w = 205.95 m/s (the paper's critical constraints, previously
4.567 MPa / unstable); NBP p_vap = 101 324 Pa (previously 101 170); λ matches
NIST far-field again — 0.0344 at 300 K/1 atm, 0.0218 at 200 K, 0.184 NBP liquid
(previously +60 %/+90 %) — and λ_cr now diverges toward Tc on the critical
isochore (0.52 W/m/K at 190.6 K → 0.004 at 250 K). Not compiled here — build and
`make check` are Christian's; the near-critical gate
`DISABLED_Methane_Cp_FD_NearCritical` should be enabled with the build.

## Nitrogen WIP review + fix round (same day, Christian's request)

Christian's in-progress `cl_GM_EoS_Nitrogen` (first six functions) reviewed against
Span et al. 2000 (Table 17/18 read from a 200 dpi page render; Eq. 53 from text).
Structure verified correct: mE basis ≡ Eq. (53), 6/26/4 term grouping, all 36
i/j/l exponents digit-exact, β/γ digit-exact, 33/36 N digit-exact. Eleven defects
found and fixed on approval (hpp/cpp):

- Data: mA(5) sign + stray digit (paper a6 = +6.678326e-8); mN k=17 sign
  (+0.635466899859e-3); mN k=34 sign + stray digit (−20.9115600730); mPhi(3)
  20 → 25 (Table 18, k=36).
- Code: `mA(8)` OOB in update_e (a8 is mA(7)); phi0_tt missing factor 2 on the
  a4 term and missing a8² on the Einstein term; update_f Gaussian fill off-by-one
  (k=33 → 32, mF(32) was NaN); update_f now calls the two power-table updaters
  (were never refreshed); phir_d exp-group double mF + wrong δ-power
  (l·δ^(l+1) → l·δ^(i+l)); phir_d Gaussian group spurious outer i factor +
  mDeltaPowL misuse — both groups rewritten as δ^i·(bracket)/δ.
- Compile blockers in the reviewed region: `override` dropped from
  compute_phi0_d (the base declares no such virtual — siblings let the base
  handle ln δ; noted in a comment), `aTau` → mTau in update_tau_pow_j.
- Left as WIP (Christian's active tail): phir_dd/t/tt bodies (syntax + math
  noted in session), ctor/base wiring, critical constants, doi typo.

## Enthalpy reference: adjudicated, CLI changed

Christian observed h = 0 at 298.15 K in the `gastable` executable and asked if
the campaigns changed the reference convention. **They did not** (08-05 fixed
only `h_ref()` units). Traced mechanics: construction discards the data file's
b1/b2 (`cl_GT_RefGas.cpp:252-267`) and re-anchors
H(298.15) := ΔfH°(298.15) + [H°(298.15)−H°(0)] (Hf = record line 2 assigned
enthalpy; Href = interval-line H298−H0 field), which puts h(0 K) = ΔfH°(298.15)
— 0 at 0 K for reference elements, i.e. Christian's preferred scale was already
the internal one. The observed h = 0 at 298.15 was purely the CLI display
subtracting `H_ref()` (≡ H(298.15) post-anchoring). Per his ruling the
subtraction is removed from both branches (`main.cpp`); tables now print the
internal 0 K-anchored scale (N2 shows h(298.15) ≈ +8.67 kJ/mol, h → 0 toward
0 K).

**Physics flag raised, then RESOLVED against the consumers (Christian's
instruction to check combustion and the channel model).** The concern was that
the internal scale attaches ΔfH°(298.15) at 0 K — a hybrid of CEA's convention
(H(298.15) = ΔfH(298.15)) and a true 0 K scale — so that Σν·H_i would differ
from the reaction enthalpy by Σν·[H°(298)−H°(0)]_i, ≈ −2.9 kJ/mol for
H2 + ½O2 → H2O. **It does not propagate.** Both consumers are immune, for
different reasons:

- **Combustion re-references explicitly.** `cl_CN_Scheme.cpp:88-91` builds
  `mH0(k) = Hf(k)/M(k) − h_k(gTref)` and `:315` returns
  `mH(k) = h(k,T,p) + mH0(k)`, i.e. ΔfH°(298.15) + [h(T) − h(298.15)] — the
  CEA convention, assembled from the difference only, so the table's own zero
  cancels identically. This is what feeds the energy Jacobian (`:433`), the
  temperature rate (`:472`) and the heat release (`:539`), which are exactly
  the Σ h·rate sums the flag was about.
- **The channel model only takes differences.** Wall heat is
  `h(T,p) − h(T_wall,p)` (`cl_CH_ChannelODE.cpp:506-507`) and the total-state
  residual is `h(T,p) + ½u² − h_total` (`cl_Channel.cpp:189`); a constant
  offset cancels in both.

So the convention is unobservable to everything in the tree, which makes it a
documentation matter rather than a physics defect. It is now stated precisely
in the `RefGas` class doc — including the correction that H(298.15) is the
formation enthalpy **plus** the sensible term, not the assigned enthalpy alone,
which an earlier draft of that block got wrong. It would only matter if an
absolute enthalpy were ever compared against an outside table.

Equilibrium is a separate path (`Gibbs`/`dGibbsdT`) whose b1/b2 offsets were
already established as load-bearing; unchanged and not re-opened here.

## Status

- [x] Chung β settled: transcription error, 1.368 → 1.3168 (P1 #5)
- [x] MC-SRK c1 settled: order-of-magnitude typo, 0.16054 → 1.60539 (P0 #4)
- [x] PR-78 κ settled: wrong digits but dead code (P2 #8)
- [x] Friend λ_cr closed: helpers correct but never called; full fix spec above
      (round 2; was "coefficient tables blocked" in round 1)
- [x] Hydrogen M9 double-pv suspicion refuted analytically
- [x] Methane phir_tt fix independently confirmed
- [x] ~~Acquire Friend 1989 and Schmidt & Wagner 1985~~ — supplied same day;
      round 2 closed methane transport, oxygen, and cubicalpha.inp
- [x] Methane pc adjudicated: 4.5992 MPa; the paper itself carries the 4.5922 typo
- [x] Defect tick-off list written: `../todo/closed/gas_correctness_fixes.md`
- [x] Fix round executed on Christian's approval (same day) — all items except
      D2b; numerically validated, see the fix-round section
- [x] Codex + Grok jury review of the fix diff — no blocking defect in the
      applied algebra (3-way agreement on D1-D9/H1); PR-78 1.487503 challenge
      REFUTED against alpha.txt:440; confirmed follow-ups: mPmax guard misses
      `Helmholtz::v` and the transport's own `update_Tp` (design → Christian),
      F_Q blind to the new pH2/He3 CAS keys (dormant), mTmax never asserted
      (pre-existing). Thread `tmp/ai_exchange/review_gas_literature_fixes.md`
- [x] Jury follow-ups executed on Christian's rulings: mPmax guard upgraded to
      an always-active `BELFEM_ERROR` in `Helmholtz::v()` (the choke point every
      (T,p) path funnels through, incl. the transport classes' own `update_Tp`;
      the redundant debug assert in `update_Tp` removed); F_Q gains the pH2 key
      (Q = 0.76, same molecule) and a no-published-Q comment for He-3
- [ ] Build + `make check` (Christian); enable
      `DISABLED_Methane_Cp_FD_NearCritical` (D1 gate)
- [ ] D2b — MLSG scaled-EoS region: Christian to rule
- [ ] mTmax guard: still unguarded everywhere (pre-existing); a hard T error
      would change behavior for models like O2 (Tmax = 300 K) — needs a ruling
