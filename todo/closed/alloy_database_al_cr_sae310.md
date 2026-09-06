# Alloy Database Expansion: Aluminum, Chromium, and the Synthesized SAE 310 Gate

**Date:** 2026-08-15
**Purpose:** Extend the elemental material database with Aluminum and Chromium so that the
`Alloy` homogenization tool can synthesize HTS-relevant alloys (first target: SAE 310
stainless), with lambda and rho data completed wherever possible. Mechanism in one sentence:
elemental Bloch-Grueneisen `rho_i` + Hust `lambda` feed `Alloy::create_tables()`, whose
alloy-level `rho_0` carries the disorder scattering that linear mixing cannot produce.
**Module:** `../../src/physics/materials`
**AIs involved:** Claude (plan + numeric verification), Christian (fits, rulings), Codex
(prose pass, applied 2026-08-15)
**Status:** DRAFT; plan written 2026-08-15, not yet audited; no source modified

> **Scope guards:**
> - Superconductivity of the metals is out of scope by design (database serves HTS; only YBCO
>   carries jc/n/power law — standing ruling 2026-08-14).
> - Kohler magnetoresistance and mechanics data (E, nu) for Fe/Ni/Al/Cr are follow-on work, not part of
>   this plan (Christian, 2026-08-14/15); the Alloy mechanics path tolerates their absence
>   only if the synthesized alloy skips E/nu — see O3.
> - Chromium's spin-density-wave antiferromagnetism (Neel point 311 K) is deliberately
>   smoothed, not modeled — Cr enters only as alloy feedstock; see Appendix A.
> - High-temperature use of synthesized Fe-based alloys is capped: elemental Fe/Ni cp carry
>   ferromagnetic contributions a paramagnetic austenite does not have — see R7/O2.

---

## 1. Current Behaviour and Where It Falls Short

The elemental side is largely ready: Fe and Ni carry the new two-Bezier cp
(`Metal::create_cp`), Bloch-Grueneisen `rho_i` anchored against White & Woods
(`cl_Material_Iron.cpp` / `cl_Material_Nickel.cpp`, `create_debye_and_rho()`), and Hust lambda
coefficient sets (Ni set verified against handbook values 100-631 K on 2026-08-15; Fe against
NBSIR 84-3007). Cu, Ag, Pb, Sn, In are complete for cp; Al and Cr do not exist yet.

`Alloy::create_tables()` (`cl_Material_Alloy.cpp:271-420`) synthesizes the following quantities at each 4-K knot:
self-consistent Hill K/G with full tensors, Schapery alpha, Neumann-Kopp cp
(mass-fraction-weighted), mass-weighted Sommerfeld gamma, and two-channel transport:

| Channel | Model | Citation |
|---|---|---|
| rho | `rho_0 + dot( v, rho_i_components )` — alloy-level residual + volume-weighted elemental phonon resistivity | `cl_Material_Alloy.cpp:409-417` |
| lambda | `1 / ( w0 + w_i )`, `w0 = rho_0 / ( L0 T )`, `w_i = dot( v, hust( components ) )` — pure electronic Wiedemann-Franz on the alloy residual + intrinsic component resistivities | `cl_Material_Alloy.cpp:414-419` |

The disorder resistivity is a **parameter, not a mixture**; that design decision makes alloy
synthesis viable at all (SAE 310: measured rho(RT) ~ 87e-8 Ohm m vs ~10e-8
mass-weighted from elements; the difference is the fitted `rho_0`).

| Failure / gap | Mechanism | Evidence |
|---|---|---|
| No Al, no Cr | classes do not exist | `../../src/physics/materials` |
| Alloy lambda underestimates stainless ~40% at RT, ~2x at 4 K | lattice (phonon) heat conduction channel missing — model is purely electronic | model `1/(w0+w_i)` ~ 8-9 W/(m K) at 300 K vs NIST SS310 ~ 13-14; ~0.12 vs ~0.25 at 4 K (numeric replication 2026-08-15, medium-high confidence) |
| Synthesized Fe-alloys inherit spurious magnetic cp above ~500 K | elemental Fe cp fit (T_max 860 K) contains the Curie approach; austenite is paramagnetic | `cl_Material_Iron.cpp` cp fit; dcp/dT rising to 0.72 at 860 K |
| Cr Neel anomalies unrepresentable in `Ferromagnetic` | Brillouin/ferro machinery does not model SDW antiferromagnetism | `cl_Material_Ferromagnetic.cpp` |

**Bottom line:** the elemental pipeline and the Alloy architecture are sound; what stands
between "synthesized 310" and "matches NIST" is two elements and one missing lattice-conduction
term, gated by a measured-alloy comparison that does not exist yet.

## 2. Architecture: Why Element Feedstock + Alloy-Level Residuals Is the Right Spine

Per-alloy direct fits (the HastelloyC276 pattern) are the rejected alternative for the general
case: they scale linearly in effort with the number of compositions and extrapolate poorly
below the measured range. Element feedstock + mixing rules + one fitted `rho_0` (+ one lattice lambda
term, R6) scales to arbitrary compositions and keeps the low-temperature physics (gamma, beta,
Bloch-Grueneisen shape) anchored in element data that exists down to 2 K. Direct fits remain
the fallback for alloys where the synthesis gate (R8) fails; HastelloyC276 stays as-is.

## 3. Gap Table

| # | State | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| 1 | Al: cp (Bezier), rho_i (BG), lambda (Hust), constants | stabilizer + Al-alloy feedstock | no | (c) | R1-R2; theta_D ~ 428 K, gamma 1.35 mJ/(mol K^2); Woodcraft 2005 for lambda |
| 2 | Cr: cp, rho_i, lambda, constants | stainless feedstock | no | (c) | R3-R4; theta_D ~ 630 K; smooth through Neel 311 K (Appendix A) |
| 3 | Lattice lambda channel in Alloy | stainless lambda within ~10% of NIST | no — purely electronic | (c) | R6; quantified gap in §1 |
| 4 | Alloy-level rho_0 for SAE 310 | rho synthesis | mechanism exists, value unfitted | (a) — from measured rho(4 K) | `cl_Material_Alloy.cpp:309` reads `constant_property( rho_0 )` |
| 5 | NIST SS310 reference curves in-repo | R8 gate | no | (c) | R7; NIST cryogenic fits 4-300 K |
| 6 | T_max policy for synthesized Fe-alloys | avoid spurious Curie cp | no | (b) | O2 |
| 7 | E/nu for Al/Cr/Fe/Ni (Alloy mechanics path) | only if synthesized alloys need mechanics | partial | (b) | O3; mech deferred by scope guard |
| 8 | Literature registration | traceability | no | (c) | R9; Woodcraft 2005, NIST cryo, Touloukian Cr volumes into `../../doc/literature_references.md` |

### 3.1 Cross-cutting findings

- **Every new element must follow the post-2026-08-14 conventions:** beta/gamma/rho_0 set once
  in `set_constants()` (single source, avoiding the Iron double-set defect), `set_RRR` only after
  `rho_i_ref` exists (the Nickel ordering defect), no shadow cp members, and lambda coefficient
  sets with `p0 = p7 = 0` unless the metal genuinely carries the RRR-scaled `w_i0` correction.
- **Verify lambda fits through the same evaluation chain before committing.** The Ni "rubbish
  dataset" was caught only because the replication evaluated `Metal::lambda_custom` + `hust`
  end to end against handbook anchors. Make that the standard pre-commit step for every new
  coefficient set (the scratchpad battery from 2026-08-15 is the reference procedure).

## 4. Ordered Steps

- [ ] **R1 — Aluminum: constants + cp.** `Metal` subclass; M from Abundance, gamma/beta from
  literature (theta0 ~ 428 K), T_max at melting 933.47 K; two-Bezier cp against Touloukian.
  Verification: the standard cp battery (continuity, monotonicity, spot values; cp(300) ~ 897
  J/(kg K), Dulong-Petit 3R/M = 924).
- [ ] **R2 — Aluminum: transport.** (after: R1) BG rho_i anchored at a documented reference
  (rho(273) ~ 2.42e-8 Ohm m); Hust lambda per Woodcraft 2005 (RRR-scaled, because Al is the
  RRR-critical stabilizer, so the `w_i0` correction with nonzero p0/p7 is likely wanted here,
  unlike Fe/Cr); `rho_0 = 0` default + set_RRR at ctor end.
- [ ] **R3 — Chromium: constants + cp.** Plain `Metal` (NOT `Ferromagnetic`); smooth fit
  through the 311 K Neel anomaly with the Appendix-A comment; theta0 ~ 630 K, T_max at the
  melting point or a documented validity cap.
- [ ] **R4 — Chromium: transport.** (after: R3) Smoothed BG rho_i (rho(295) ~ 12.5e-8) and
  Hust lambda (lambda(300) ~ 94 W/(m K)), both explicitly documented as SDW-smoothed
  feedstock curves.
- [ ] **R5 — Elemental verification sweep.** (after: R2, R4) Run the full numeric battery on
  Al and Cr (cp + transport, handbook anchors, WF cross-checks at 77 K) and record in devlog.
- [ ] **R6 — Lattice conduction channel in `Alloy::create_tables()`.** Additive
  `w_lattice^-1 = a T^b` (or equivalent) with per-alloy coefficients; default off (a = 0) so
  existing Alloy users (solder) remain bit-identical. Design in §6.
- [ ] **R7 — NIST SS310 reference fits in-repo.** cp and lambda 4-300 K as the comparison
  target (data arrays or fitted polynomials in the test, not a new material class);
  rho(4 K)/rho(300 K) anchors for the `rho_0` fit.
- [ ] **R8 — The synthesis gate.** (after: R5, R6, R7) Synthesize SAE 310 (Fe-25Cr-20Ni by
  mass; Mn/Si/C ignored, since they are invisible in cp/density and absorbed in rho_0), fit `rho_0` and the
  lattice term, then diff against R7 across 4-300 K. Acceptance: cp within ~5%, lambda within
  ~10%, rho within ~5% over the range. This is the falsification gate for the whole element
  pipeline; failure modes route back to the gap table rather than into silent per-alloy tuning.
- [ ] **R9 — Documentation.** Register Woodcraft 2005, the NIST cryogenic fits, and the
  Touloukian Cr volumes in `../../doc/literature_references.md`; update
  `../../src/physics/materials/doc/README.md` (element inventory + Alloy transport model); devlog.

## 5. Open Design Questions (not silently decided)

- **O1 — Lattice term shape.** Simple `a T^b` per alloy, or the two-parameter
  Callaway-lite form used in cryogenic steel fits? `a T^b` is likely sufficient for 4-300 K;
  decide against the R7 data.
- **O2 — T_max policy for synthesized Fe-alloys.** Cap at ~500 K (below the elemental Curie
  contamination, see §1) or subtract a magnetic-cp estimate from the Fe/Ni feedstock above
  400 K? Capping is simpler and covers all HTS use; ruling needed only if a user asks for
  hot-structure temperatures.
- **O3 — Mechanics feedstock.** `Alloy::create_tables()` evaluates `E(T)`/`nu(T)` of every
  component unconditionally (`cl_Material_Alloy.cpp:340-342`). Fe/Ni/Al/Cr mech is deferred by
  scope guard. Does the Alloy path need a "skip mechanics" mode until then, or do the four
  metals get E/nu fits first? Blocks R8 unless resolved.
- **O4 — SAE 310 vs 310S composition.** 310 (0.25 C) vs 310S (0.08 C) differ in carbon only;
  NIST fits are typically for 310S. Decide which composition the gate targets and say so in
  the comparison.

## 6. Lattice Channel Design (R6 artifact)

```cpp
// Alloy::create_tables(), lambda assembly — current:
real w0 = rho_0 / ( constant::L0 * T );
tLambda( i ) = 1.0 / ( w0 + w_i );
// proposed: parallel lattice channel, default-off
// lambda_total = 1/( w0 + w_i ) + a * pow( T, b )    ( a = 0 -> today's behaviour )
```

Coefficients `a`, `b` as alloy-level constants (candidates for two new `MaterialProperty`
entries — note the enum-renumber audit of 2026-08-14: appending at the tail before `UNDEFINED`
is safe in-tree, user-material .so plugins must be rebuilt). Exact storage TBD in R6.

## 7. Definition-of-Done Checklist

- [ ] Every gap-table row mapped to a step or an On.
- [ ] Al and Cr pass the elemental battery (R5) with devlog record.
- [ ] Existing Alloy users (solder) remain bit-identical with the lattice term defaulted off.
- [ ] R8 gate passes with the acceptance bands, or its failure is decomposed into new gap rows.
- [ ] Literature registered; materials doc README updated; input contract untouched (materials
  are C++-configured; verified no input.conf coupling 2026-08-14).

## 8. Audit Trail

- Numeric groundwork (Fe/Ni transport verification, Alloy mixing-rule analysis, stainless
  lambda gap quantification): devlog `dl20260814_bezier_cp_copper_audit.md`, session of
  2026-08-14/15.
- Codex prose pass applied 2026-08-15 (15 wording corrections; exchange thread
  `../../tmp/ai_exchange/alloy_database_plan.md`, ephemeral). Content audit still open.
