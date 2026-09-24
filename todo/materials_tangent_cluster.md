# Debt Cluster A: Materials Tangents and Low-T Fits (DR-163, DR-164, DR-170, DR-173, DR-174)

**Date:** 2026-09-23
**Purpose:** Close five open debt-register rows that share one module (`src/physics/materials`), one test pair and one gate shape — four rows carry a `DISABLED_` pin that is red today and green after the fix; DR-170's test is enabled and goes red only once its `dJ` call is added.
**Module:** `src/physics/materials` (+ `debye.f90`)
**AIs involved:** Claude Fable (plan, juries, close-out), Codex + Grok (plan and code juries, terra/high + grok-4.7/high), Opus (implementation), the `belfem-3c` reviewer session (coordination)
**Status:** CODE COMPLETE 2026-09-24 01:45 — all five rows LANDED on `devel` (DR-173 `fcb8536`, DR-163 `31737f5` added to the three below), each behind both juries and green gates in both trees; DR-164/170/174 struck on Christian's ruling, DR-173/DR-163 strikes owed; the debug-assert leg of every gate owed; residue on the rows. Earlier state: three of five rows LANDED 2026-09-23 23:15 on `devel` (DR-174 `e3a9b6b`, DR-164 `b5c42a9`, DR-170 `9373e44`), each behind the plan jury, the code jury (Codex terra/high + Grok 4.7/high, reconciled) and `test_physics` 147 + 1 designed skip plus `make check-fast` 8/8 in BOTH trees (release only — the debug-assert leg is owed). Remaining: R4 (DR-173, ruling (c) vs (b), probe numbers in R4) and R8/R9 (DR-163, both physics voices recommend A; ruling). **Quiz answered 2026-09-24 00:05:** DR-173 → (c); DR-163 → A; DR-164/170/174 STRUCK and archived; the DR-170 guard band accepted as a bounded approximation. R4 and R9 in implementation (stage 3)

## Context

Christian asked which open cluster to approach next and chose this one. The peer reviewer of the
v0.9.2..HEAD delta (`belfem-3c`, three vendor juries reading the tree) has
**`powerlaws.hpp` in its read set with the physics jury not yet dispatched** and asked us to hold
edits there until it says the round is back; `cl_JcFunction_ModifiedKim.*`, `cl_Material_Lead.cpp`,
`debye.f90`, `cl_Material_YBCO.cpp` and the tests are outside its read set. The plan is therefore
**staged**: stage 1 touches only the free files; stage 2 (`powerlaws.hpp`) waits for the
reviewer's release. The reviewer also handed us one check (commit `99ac9d0` made the constant-n
value overloads read `constant_property( n )` instead of dereferencing a null `mNFunction`;
derivative overloads sharing the `is_constant( n )` precondition must read n the same way).

The register rows are the source of the defects, but two of their fix hypotheses do not survive
a source read (DR-173, DR-174 below) — the plan follows the tree, not the row text.

## Scope guards

- OUT: every other register row; `rho_riva`; `deval_dbeta` (the angle hook must stay unwired,
  `cl_JcFunction.hpp:220-230`); `mRhoMin` floor semantics (follow-up R12); any `input.conf` key.
- OUT: no refit of lead's Touloukian Bézier control points unless the probe in R3 shows the
  Bézier region itself is wrong (it is `exp(...)`, so it cannot go negative).
- The cohomology core is untouched. No file in the reviewer's read set is edited before its release.

## 1. Defects and mechanism, per row

| row | defect (source) | fix shape | stage |
|---|---|---|---|
| **DR-164** P2 | `JcFunctionModifiedKim` has no `deval_dB`; base returns 0 (`cl_JcFunction.cpp:52-58`), so `drho_powerlaw_dB` early-outs to 0 at `powerlaws.hpp:1976` (and the defect overload at `:2006`) while the value depends on B | analytic derivative of `jc = Jc0 (1 + s·B/B0)^(−α)`, `s = √(k² sin²θ + cos²θ)`: `∂jc/∂B = −α · jc · s / (B0 + s·B)`; override in `cl_JcFunction_ModifiedKim.hpp` (header-only class, no `.cpp` exists) | 1 |
| **DR-174** P3 | **Located by instrumentation (probe 5b):** the optical electron–phonon prefactor `A = …·U²V/(V−1)` with `V = exp(U)`, `U = ω_opt/(k_B T) = 720.8/T` (`debye.f90:326-329`) — `exp` overflows at U > 709.8 (the module's `quadruple` is `selected_real_kind(r=60)` = **binary64**, `intpoints.f90:16`, verified), so T < 1.016 K gives `Inf/Inf = NaN` in `K(5)`. A second, masked overflow of the same class sits in `scatter_function` (`:195-197`, `(e^z)²` overflows at z ≈ 355, T < 0.77 K) and in `debye_function` (`:39-48`) | one numerical contract, exact in double: `U²/(1 − e^{−U})` at `:329`; `z^n e^{−z}/(1 − e^{−z})²` in both integrands. Probe 5c: finite to 0.05 K, identical to the library to 1e-15 from 1.5 K up. The `T <= 0` guard stays | 1 |
| **DR-173** P2 | public `cp()` of lead negative at 0.5–2 K. **The row's diagnosis ("first fit region wrong below its switch") is refuted by reading:** `cp_custom` (`cl_Material_Metal.hpp:764-787`) returns `γT + βT³` below `T0 = 0.02 θ ≈ 2.1 K` and `exp( poly5( ln T ) )` up to the first Bézier point (`e^2.026 ≈ 7.6 K`) — both strictly positive. The public path is the cubic spline built by `SplineLookupTable::create_spline` (`cl_Material_SplineLookupTable.cpp:131-148`): knots at 0, 4, 8, … K, start BC = tangent γ at 0 K, value 0 at 0 K. A cubic Hermite segment from (0, 0, slope γ = 0.0147) to a value at 4 K that a T³ term has already grown large (lead: β = 8.3e-3, θ ≈ 105 K; indium's β = 0.0141 is larger yet its public cp stays positive — the sign depends on the whole first segment, not on β alone) dips below zero inside [0, 4 K]. Confidence **medium-high** (source trace, not executed; the probe in R3 settles it) | probe first (R3); then the fix that the probe supports — the leading candidate is to make the spline's first segment T³-consistent: evaluate `cp_custom` directly below the first knot (`cp_spline` dispatch), or change the knot grid (the start BC is already `Tangent` with the correct slope γ = cp'(0)). Decided after the probe, not before | 1 |
| **DR-170** P3 | 4-argument `drho_piecewise_dJ` `BELFEM_ERROR`s on a vanishing discriminant (`powerlaws.hpp:1839-1842`) — that overload carries three always-on errors (`:1835` on `a`, `:1840` on `d`, `:1860` on `imag`); the seven other `dJ` overloads assert on `a` only and still evaluate `1/√(b²+ac)` (`:1385` pattern), so ALL EIGHT are non-finite at the exact degenerate point in every build; `drho_piecewise_dT` returns a guarded fallback (`:2305-2307`). Ruled 2026-09-22: `n = 4` is legal, the abort tier is the defect | in all eight `drho_piecewise_dJ` overloads (`:1324`, `:1409`, `:1495`, `:1566`, `:1637`, `:1714`, `:1789`, `:1873`) guard the discriminant BEFORE any root, `!( d > BELFEM_EPSILON·(b²+1) )`, and return the exact endpoint limit `mNff·ρ_n/J` (probe-verified, see R6); the `a` and `imag` errors of the 4-argument overload are kept (n = 7 makes `a` vanish and the code divides by it — pre-existing window, residue) | 2 |
| **DR-163** P2 | all eight `drho_piecewise_dJ` overloads divide by the parallel factor `(1 + ρ_PL/ρ_n)²` below `j1` (`:1354 … :1899`, eight sites) while `rho_piecewise` returns the raw power law there (`:726-728`); `dB`/`dT` are raw (`:2044-2049`, `:2189-2200`). **Which side is right is physics, not a vote** (`todo/test_gap_followups.md` R11: Christian leans "parallel is right" and asked for a literature-backed jury) | (a) physics round first — Duron et al. 2004 (the parallel combination) and the piecewise header `powerlaws.hpp:1291-1295`; Messe 2023 §2.6 is the Kim-model pointer, not this; neither Duron nor Rhyner is in `literature/` — Codex + Grok single-claim round says what it could not read; (b) consumer inventory: if the VALUE changes, every piecewise deck's residual changes (`cl_FEM_Calculator.hpp:3828-4074` reaches `drho_piecewise_dJ`; `rho_piecewise` sites in the same file); if the TANGENT changes, only Newton iteration counts move; (c) Christian's ruling; (d) apply on both sides consistently — either drop the factor from the eight `dJ` sites, or add the parallel combination to `rho_piecewise` AND its `dB`/`dT` twins | 2 (after ruling) |

Bottom line: four rows are mechanical once ruled; DR-163 carries the one physics decision and is
sequenced last so the rest of the cluster does not wait on it.

## 2. Steps

**Stage 1 — free files (may start once the plan is approved):**

- [x] **R1** (landed `b5c42a9`, 2026-09-23 23:15) `cl_JcFunction_ModifiedKim.hpp`: add `deval_dB( normB, angle, T ) const override`
  with the closed form above (**probe-confirmed 2026-09-23 21:05**: worst 1.8e-10 vs central FD
  over 25 (θ, B) points on the pin's constants); `deval_dT` stays the base zero (the law has no T). Update the
  Doxygen contract note in `cl_JcFunction.hpp:200-207` that names ModifiedKim as "not implemented".
  Consumers that wake up: `drho_powerlaw_dB` (`powerlaws.hpp:1976`, `:2006`; calculator
  `cl_FEM_Calculator.hpp:3670-3764`) AND `drho_piecewise_dB` (`:2034-2037`, `:2068-2071`;
  calculator `:3717, :3741, :3775, :3787`) — a Kim material on either law gains a field tangent.
  The pin compares against the full `rho_powerlaw` chain (parallel factor included), so what must
  match is the whole tangent, not `∂jc/∂B` alone.
- [x] **R2** (landed `e3a9b6b`, 2026-09-23 23:15) `debye.f90` — **mechanism LOCATED 2026-09-23 21:14, fix probe-verified 21:16**.
  Not `scatter_function`: the NaN is born in the optical electron–phonon block —
  `U = omega_opt/(kB*T)` = 720.8/T (`:326`), `V = exp( U )` overflows to Infinity (`:327`), and
  `A = … * ( U*U*V )/( V − 1 )` (`:329`) is Inf/Inf = NaN, which poisons `K(5)` and the whole
  integral. Overflow at exp(709.8) ⇒ T < 1.0155 K, matching the row's "below about 1 K"
  (instrumented scratch copy: inputs finite, `A` = NaN at 0.05–1.0 K, finite at 1.5 K). Fix:
  `( U*U )/( 1 − exp( −U ) )`, algebraically identical, finite for all U > 0 — on the scratch copy
  the rewritten integrator is finite down to 0.05 K (k_ph 9.0e-7 … 5.9e-3 W/(m K)) and identical
  to the library to 1e-15 from 1.5 K to 300 K. Also apply the `e^{−z}` form to `scatter_function`
  and `debye_function` (a second, masked overflow of the same class: `(e^z)²` overflows double at z ≈ 355,
  i.e. T < 0.77 K for θ = 275 K); `debye_derivatives` (`:52-79`, `debye_table` path) shares the form and gets the
  same treatment. Answered 21:20 (Grok's reading, my `kind.f90` check under the module's flags):
  `quadruple = selected_real_kind( r = 60 )` resolves to **kind 8 = binary64** (range 307, huge
  1.8e308) — the Fortran rule picks the smallest precision whose RANGE ≥ 60, and that is double.
  The whole `debye` module therefore runs in double despite its naming; the fix above is exact in
  double, and no kind change is proposed here. Residue for the DR-174 row at close-out: the
  misleading kind name (`intpoints.f90:16`), to be widened onto the row rather than filed anew.
- [x] **R3** Probe (no edit) — **RAN 2026-09-23 20:58** (scratchpad `probe_dr173/`, linked
  read-only against `build/lib/libbelfem.a` of 19:24, Armadillo release). Result: `cp_custom` is
  positive on the whole 0.1–8 K grid; the public spline equals it at the 4 K and 8 K knots
  (to 1e-14) and is **negative from 0.3 K to 2.1 K, minimum −3.23e-2 J/(kg·K) at 1.5 K**
  (19 of 80 points ≤ 0). The other eight builtin metals (copper, silver, aluminum, nickel, iron,
  tin, indium, chromium) stay positive down to 0.1 K. Diagnosis CONFIRMED by execution: the
  defect is the 4 K public spline on a T³ start, not the Touloukian fit. By-catch: between knots
  the spline sits up to +0.148 J/(kg·K) (≈ 5 %) above `cp_custom` at 5.8 K — the knot spacing is
  coarse for lead's low θ everywhere below ~10 K, not only in the first segment.
- [x] **R4** DR-173 fix — **RULED 2026-09-24 00:05 (Christian, quiz): option (c), 1 K knots for cp only; LANDED `fcb8536` 01:45** (step parameter guarded finite and positive after the code jury). The two candidates as they stood before the ruling:
  - (a) *dispatch below the first knot*: `cp_spline`/`dcpdT_spline`/`d2cpdT2_spline` return
    `cp_custom`/`dcpdT_custom`/… for `T < first interior knot` (4 K). Local, every metal keeps its
    values above 4 K; below 4 K every metal moves to its own analytic law (small for the eight,
    the whole defect for lead). Cost: a slope jump at 4 K unless the spline's start segment is
    also matched — the probe's second run measures it.
  - (b) *finer knot grid* (`dT = 1 K` in `SplineLookupTable::create_spline`, `:133`): removes
    the undershoot and the 5 % mid-knot error together, but changes every metal's public cp, λ,
    ρ, … spline by up to a few percent between the old knots → consumer impact on every thermal
    deck; needs its own inventory and A/B before it lands.
  **Probe 2 (20:59, all nine metals):** candidate (a) leaves a C1 break at the 4 K knot —
  public slope vs `dcpdT_custom` at 4 K: lead +18 %, **tin −35 %**, nickel −2 %, the rest < 1 %.
  By-catch: the 4 K grid misrepresents every soft metal below 4 K, not only lead — max relative
  deviation of the public spline from `cp_custom` on 0.1–12 K: lead 201 % (the negative window),
  **tin 48 % at 1.4 K, indium 16 % at 0.8 K**, aluminum/nickel/iron/chromium ≤ 1.6 % (and those
  sit at 10 K, the Bézier region), copper/silver ≈ 0. So (a) fixes the sign for lead but leaves
  tin and indium wrong below 4 K and breaks C1 for all three; (b) fixes all three at the cost of a
  global value change. A third shape, (c): (b) restricted to cp only. The API for it does not exist yet (Codex):
  `SplineLookupTable::create_spline` hard-codes `dT = 4.0` (`:131`); (c) adds a trailing
  `aStep = 4.0` parameter to `Material::create_spline( aProperty, adYdX0, adYdX1 )` and its
  `SplineLookupTable` overload, every existing caller unchanged, and only the cp call at
  `cl_Material_Metal.cpp:218` passes 1.0 — consumer impact limited to cp.
  **Probe 6 (21:30) — the knot grid, measured.** The public cp spline rebuilt from `cp_custom`
  with the library's own `Spline` + `create_helpmatrix` (Tangent start at γ, Parabolic end):

  | dT | lead min / max rel dev | tin max rel dev | indium max rel dev | copper max rel dev | lead pin |
  |---|---|---|---|---|---|
  | 4 K (tree) | **−3.23e-2** / 201 % | 48 % | 16 % | 1.2 % | FAILS (sign and envelope +0.34) |
  | 2 K | +1.5e-3 / 9.8 % (0.7 K) | 0.6 % | 4.1 % | 0.09 % | ok |
  | **1 K** | +1.5e-3 / **0.6 %** (0.4 K) | **0.03 %** | **0.09 %** | **0.01 %** | ok |
  | 0.5 K | +1.5e-3 / 0.06 % | 0.01 % | 0.00 % | 0.00 % | ok |

  Knots at 1 K: lead 602, tin 401, indium 431, copper 1001 — a few kB each, evaluation cost
  unchanged (`find_col` + one cubic). Every candidate spline moves TOWARD its generating function,
  so the "value change" of (b)/(c) is an accuracy gain of ≤ 1.2 % for the hard metals and the
  whole defect for the soft ones. **Recommendation for the ruling: (c) — `dT = 1 K` for cp**
  (`create_spline` gets a per-property step; λ, ρ, … keep 4 K, so the consumer impact is cp only),
  with (b) — 1 K for every property — as the simpler one-constant alternative if a global change
  is acceptable. (a) is withdrawn: it leaves tin/indium wrong below 4 K and breaks C1 at 4 K.
  **No default is taken overnight**: DR-173 goes to the morning quiz as (c) vs (b); the pin stays
  `DISABLED_` until it lands. Probe sources: scratchpad `probe_dr173/probe{,2,6}.cpp`.
- [x] **R5** (Kim and YBCO pins enabled in `b5c42a9`/`e3a9b6b`; the lead pin in `fcb8536`; the raw-residual pin in `31737f5`) Enable pins: rename `DISABLED_KimFieldDerivativeMatchesFiniteDifference`,
  `DISABLED_LeadCpIsNegativeBelowThreeKelvin`, `DISABLED_YbcoRawLambdaIsFiniteBelowOneKelvin`;
  remove lead from the known-defect table in `test_MaterialValues.cpp` (~line 126).

**Stage 2 — `powerlaws.hpp` (only after `belfem-3c` messages that its physics round is back):**

- [x] **R6** (landed `9373e44`, 2026-09-23 23:15; the `a`/`imag` errors kept per the plan jury) DR-170: replace the `BELFEM_ERROR` on the discriminant at `:1840-1842` by a guard
  in the `dT` form (`!( d > BELFEM_EPSILON·(b²+1) )`, tested on the discriminant BEFORE any root
  is taken) returning the **exact endpoint limit `mNff * rhon / normJ`** — not a frozen value.
  Derivation: `dlogrho_dt = 2(1−t)(log ρ_n − log ρ_1)` and at `n = 4` (`a == b`) `1−t = −√d/a`,
  so the `1/√d` in `dt_dc` cancels; the limit is `ρ_n·n_ff/J`. **Confirmed by probe 2026-09-23
  21:02** (`probe3`, constant-n YBCO fixture, `n = 4`, T = 77 K): analytic `dJ`, central FD and
  `n_ff ρ_n/J` agree to 1e-6 at `1 − J/j3 = 1e-7`, and `an/limit = (1−ε)³` exactly — the Bézier's
  upper end is the flux-flow law `ρ_n (J/j3)^{n_ff}` there. `dT` stays finite at the same points.
  Apply to ALL EIGHT `dJ` overloads — the seven siblings assert on `a` only and divide by
  `√(b²+ac)` unguarded (`:1385` pattern), so they are non-finite at the point too (Grok).
  Add the `dJ` call at `J = j3·{1, 1−1e-7, 1−1e-3}` to `PiecewiseDegenerateBlendStaysFinite`
  (`test_MaterialDerivatives.cpp:551`, an ENABLED test — DR-170 has no `DISABLED_` pin) with
  the limit as the expected value at the endpoint, on the 4-argument and the 1-argument overloads.
  The `BELFEM_ERROR`s on `a` (`:1835-1837`) and on `imag` (`:1860-1862`) are KEPT (Codex):
  `a = R(1/n_ff − 2/(n−1))` vanishes at the legal `n = 7`, and the code divides by `a`, so a
  demotion would add a release divide-by-zero; the sibling overloads and the nine value overloads
  already only assert there — a pre-existing window, out of this row's scope (the closed
  `fix_piecewise_degenerate_window.md` declined its R8 on 2026-09-03), recorded as residue on the
  DR-170 row at close-out. R6 changes ONLY the discriminant guard, in all eight overloads.
  Tests: cover every overload family the constant-n fixture can legally call (1-arg, `( J, T )`,
  4-arg, and their `( x, y, z, t )` defect twins); the `( J, |B|, θ )` forms need a `JcFunction`
  and are exercised on the Kim fixture; any overload the fixtures cannot reach is named in the test.
- [x] **R7** (done read-only by the implementation agent, 2026-09-23 22:00: no `Material::n` call in any derivative overload; constant-preconditioned ones read `constant_property( n )`, routed ones `n_eval`) Reviewer's check, corrected (Grok): `:982`/`:995` are `Material::n`, the asserting
  accessor, not `n_eval`; `n_eval_raw` (`:77-83`) already falls back to `constant_property( n )`.
  The check is read-only: no derivative overload calls `Material::n( … )` (today: none, `this->n(`
  is absent from `powerlaws.hpp`), and every constant-preconditioned derivative reads
  `constant_property( n )` (1-argument `dJ` `:1336-1338`) or routes through `n_eval` (4-argument
  `dJ` `:1801`). Record the result; no edit expected.
- [x] **R8** (run 2026-09-23 21:18/21:24, reconciled 22:15; both voices A) DR-163 physics round (literature + Codex/Grok, single-claim depth). Start from
  **Duron et al. 2004** (the parallel combination, `doc/literature_references.md:475-478`) and
  the piecewise header `powerlaws.hpp:1291-1295`, NOT Messe 2023 §2.6 (that section is the
  routing table's pointer for the Kim model, `:128`). Neither Duron 2004 nor Rhyner 1993 is in
  `literature/`; the round says so. Consumer inventory to `tmp/ai_exchange/consumer_impact_dr163.md`
  — nine `rho_piecewise` overloads, both `dB` and both `dT` twins, `drho_piecewise_dT` consumers
  at `cl_FEM_Calculator.hpp:3200-3270`, `rho_piecewise` at `:3808…`; ruling via question form.
- [x] **R9** (ruled A 2026-09-24; LANDED `31737f5` 01:45, with the sibling-overload checks the code jury asked for) DR-163 fix per ruling on all sites of the chosen side; enable
  `DISABLED_PowerLawRegimeTangentMatchesRawResidual` (rewritten to pin whichever side was ruled).

**Close-out (text, no jury):** devlog `devlog/dl20260923_materials_tangent_cluster.md` +
`devlog/README.md`; strike/retag the five rows in `todo/debt_register.md` (6 cells, rebuilt cell by
cell; `[P]`/`[W]` counts refreshed, `scripts/check_doc_claims.py` green); tick the FOLLOWUPS file's R8 (= DR-170) and R11 (= DR-163) in `todo/test_gap_followups.md`
(their IDs collide with this plan's own R6/R8 — different files); this plan is already
registered in `todo/README.md` (§ Added 2026-09-23).

## 3. What could break

- R1: a Newton that previously ran with a zero dρ/dB block on Kim materials now carries the full
  field tangent — convergence should improve, but the sign of `s/(B0 + sB)` must match the FD pin
  at B = 0.2, 1.3, 4.1 T. Consumers: `cl_FEM_Calculator.hpp:3670-3764` (`drho_powerlaw_dB`).
- R2: quadruple-precision rewrite changes nothing above 1 K to roundoff; the YBCO λ spline
  (sampled at 0, 4, 8 … K) is unchanged unless a knot fell in the NaN region — it does not (4 K).
- R4: the cp spline serves every metal; a start-condition change must be gated on the whole
  `MaterialValues` sweep, not on lead alone.
- R6: the guard returns the exact limit `n_ff ρ_n/J` (probe 3), not an approximation. The
  tangent still jumps at `j3` — from `n_ff ρ_n/j3` below to 0 on the normal branch — but that is
  the piecewise law's own shape (ρ = ρ_n constant above `j3`, value C0), unchanged by the fix.
- R9 (if "parallel" is ruled): a value change on every piecewise deck → consumer inventory and a
  deck A/B (one shipped piecewise deck, globals compared) before it lands.

## 4. Gates (Christian runs builds; the plan names them now)

1. **Red first:** `test_physics --gtest_also_run_disabled_tests --gtest_filter='*KimFieldDerivative*:*LeadCp*:*YbcoRawLambda*:*PowerLawRegimeTangent*'`
   on the unfixed tree — record the four reds in the exchange.
2. **Green after, per stage.** Stage 1: `*KimFieldDerivative*`, `*YbcoRawLambda*` enabled
   (and `*LeadCp*` only once R4 is ruled and landed). Stage 2: `*Degenerate*` with the new `dJ`
   call (until R6 lands that call ABORTS the process on `BELFEM_ERROR` rather than failing an
   expectation — run it only after R6); `*PowerLawRegimeTangent*` enabled only after the DR-163
   ruling, rewritten to pin the ruled side. Then full `test_physics` in the debug tree and
   `ctest -R '^physics$'` in release.
3. `make check-fast` on both backends (Armadillo `build/`, Blaze `cmake-build-debug/`) — the
   shared trees are also used by `belfem2-a9`; one `make` per tree.
4. Say **reviewed** until gate 2 has run; **verified** only after.

## 5. Juries and depth

- Plan jury: slug `review_materials_tangent_cluster`, `CODEX_MODEL=gpt-5.6-terra CODEX_EFFORT=high
  GROK_MODEL=grok-4.7 GROK_EFFORT=high` (physics module, not the astra row). Pre-registration
  first; physics claims of DR-163 are routed to Christian, never settled by the jury.
- Code jury: same slug; **own-files-only diff** to `tmp/ai_exchange/<slug>.diff` (shared checkout);
  `git status` after the Grok round.
- DR-163 physics round: `literature` skill first, then a single-claim round (terra/high, grok/high).

## 6. Coordination with the reviewer (`belfem-3c`)

- Stage 2 does not start until its message "physics round is back" arrives (idle-notice
  subscription is set). Any finding it forwards on `powerlaws.hpp` is folded into R6–R9 before the
  code jury.
- Stage 1 files are outside its read set (confirmed by its reply); we still announce the start of
  edits and the commit to it.
- `belfem2-a9` (thickcuts clone) is unaffected except for the shared build trees.

## 7. Open questions for Christian

- **O1** DR-163 direction is his ruling after the physics round (R8) — not asked now.
  **Physics round, Codex leg (21:18, terra/high): recommends A** — the tangent follows the raw
  value. Load-bearing and CONFIRMED by reading (`powerlaws.hpp:741-758`): the Bézier blend's first
  control leg has log-slope exactly `n − 1`, i.e. the blend is built to join the RAW law C1 at j1;
  option B would leave a C0 jump of relative size `ρ1/ρn` and a slope kink at j1 unless the knot
  ordinate, `j2` and the dB/dT legs are re-derived — a blend redesign. **Grok leg (21:24, grok-4.7/high): also A**, same continuity derivation, plus corrections
  folded in here: the error is ≈ 4·10⁻⁶ relative at j1 for the 77 K constants (not ≲ 10⁻⁶) and
  becomes O(1) where jc collapses (defect, near T_crit) — the coded tangent → 0 while the value
  is the raw law; the piecewise raw `dB`/`dT` branches are the ones AFTER the `drho_powerlaw_*`
  legs (my earlier citations named the powerlaw legs); eight value overloads. By-catch for the
  register: `mD = 2.5` (`cl_Material.hpp:286`) is dead and the header prose "2.5 decades above
  Jc" (`powerlaws.hpp:383-384`) contradicts the formula `jc·10^{2.5/n}`. Reconciled in
  `tmp/ai_exchange/review_dr163_physics.md`. **Quiz item, ready:** rule A (tangent follows the
  documented raw value; decks keep their residual; eight one-line edits + enable the pin as
  written + fix the header note) or B (parallel value below j1 — requires re-deriving ρ1, j2, j3
  and the dB/dT legs, i.e. a third law beside `riva`). Recommended: A. **RULED A, 2026-09-24 00:05 (Christian, quiz); R9 in implementation.**
  **Literature status (read 2026-09-23):** Messe et al. 2023 §2.6 says only that HTS resistivity
  is modeled "using the parameters of the power law or variations of it" — it does not define the
  piecewise law. Rhyner 1993 and Duron et al. 2004 are **not in `literature/`**, so neither can be
  read here; the claim that the raw Rhyner law below `j1` is by design rests on
  `src/physics/materials/doc/resistivity_laws.md` §3 (which states it) and §6 ("each law's tangent
  legs differentiate that law's own residual, branch by branch"). By the module's own stated
  convention the *tangent* must follow the *value*; whether the value should be the Duron parallel
  form below `j1` is the physics call. Magnitude of either change: the doc gives ρ_PL/ρ_n ≈ 4·10⁻⁸
  at 1.05 jc for 77 K REBCO constants, so the value difference below `j1` (≈ 1.35 jc at n = 19) is
  ≲ 10⁻⁶ relative — invisible in a deck, visible only in the FD pin and in Newton counts where
  ρ_n is small. Riva 2021 (thesis, EPFL 8754) is cited by the doc as showing measured overcritical
  resistivity rising *more slowly* than the raw law — also not in `literature/`.
- **O2** Commit granularity: one commit per row (the Tier 2 convention) or one for the stage —
  default one per row.
- **O3** Whether the three `todo/` register edits from earlier this session commit on their own
  first (keeps this cluster's diff clean) — recommended yes.
