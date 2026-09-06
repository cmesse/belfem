# Gas Modules Jury Audit (gastables + gasmodels)

**Date:** 2026-08-05
**Purpose:** Overnight three-AI jury audit of the freshly migrated
`src/physics/gastables` and `src/physics/gasmodels` modules — defects,
speedup opportunities, test-suite improvements. `cl_GM_EoS_Nitrogen.*`
excluded (work in progress).
**Module:** physics/gastables, physics/gasmodels

## Method

Two full rounds of the frozen protocol (`doc/ai_collaboration_protocol.md`),
one per module: Claude pre-registration (frozen before dispatch, built from a
four-way parallel deep-read plus personal spot-verification of every
load-bearing claim), then blind parallel Codex + Grok audits via
`scripts/cross_review.sh --jury` against a scope briefing file, then a
citation-by-citation verification pass and reconciliation table. No
compile/run by any party — all verdicts are source-trace tier at best, with
shipped-data inspection (probe tier) where noted.

Full records, all findings, verification labels, and reconciliation tables:

- `tmp/ai_exchange/review_gastables_audit.md` (34 Claude findings, 12 Codex,
  22 Grok)
- `tmp/ai_exchange/review_gasmodels_audit.md` (37+ Claude findings, 12
  Codex, 20 Grok)

These are ephemeral; this devlog is the persistent distillation.

## P0 — silently wrong numbers in production-exposed paths

1. **Remix without cache invalidation** (`cl_Gas.cpp:606-691` +
   `cl_GM_EoS_Cubic.cpp:69-112,661`, 3/3): neither `Gas::remix` nor
   `EoS_Cubic::remix` resets the Statevals cache / the T-keyed a(T) cache;
   post-remix queries at an already-seen (T,p) return the old mixture's
   properties.
2. **`alpha()` gates on the KAPPA statebit** (`cl_GM_EoS_Cubic.cpp:267`,
   `cl_GM_EoS_Idgas.cpp:99`, 3/3): `kappa()` before `alpha()` at the same
   state makes alpha return the never-written slot default 0 → downstream
   `realgas_c`/`realgas_cv`/`realgas_gamma` divide by it.
3. **`realgas_c` missing cp factor** (`cl_Gas.cpp:1088-1094`, Claude only;
   Grok's "identities OK" refuted at verification with a hand derivation):
   `sqrt(p·v·β/(α·cv))` reduces via pβ = α/κ to `sqrt(v/(κ·cv))` — units
   Kelvin; correct is c² = (cp/cv)·v/κ. Probe test added (below).

**Downgraded after Christian's consumer-trace request (2026-08-06
follow-up entry in the exchange file):**

- **`HeatPoly::dSdT` operator precedence** (`cl_GT_HeatPoly.cpp:84-93`,
  3/3) — defect confirmed and the contract (dSdT = Cp/T, fully Rm-scaled)
  is pinned by the Glue/Custom/spline siblings and the `dGibbsdT` chain
  rule; no consumer adds the missing terms. **P0 → P1 latent**: the
  factory leaves every RefGas in SPLINE mode (`create_splines` →
  `set_mode(SPLINE)`), so production/combustion consumed the correct
  spline path; the broken base-poly path needs an explicit
  `set_mode(POLY)`.
- **`RefGas::h_ref()` unit inversion** (`cl_GT_RefGas.cpp:866`, 3/3) —
  defect confirmed (`Href·M` is J·kg/mol² under any reading), but its only
  consumer in the tree is the mass-specific table branch of the `gastable`
  CLI (`main.cpp:102`); all gasmodels/combustion consumers use the molar
  `H_ref()` family consistently. **P0 → P2.**

## P1 highlights

- `TransportPolyEmpty::eval/deval/ddeval` take `const real&` → silently do
  NOT override the base virtuals; gases without transport data abort in
  `poly_Mu` instead of returning 0 (`cl_GT_TransportPolyEmpty.hpp:52-64`).
- `HeatPolyEmpty` doesn't override `dSdT` → base reads 7 coefficients from
  a length-1 vector (OOB / UB in release).
- Leftover `// delete me` block in `cl_GT_RefGas.cpp:526-535` overwrites the
  transport-glue curvature test array → acceptance check passes
  unconditionally + stdout spam.
- 24 unchecked `std::stod/stoi` sites across the parsers = `std::terminate`
  without context under release `-fno-exceptions`; in-house non-throwing
  `to_real` (core stringtools) exists as replacement.
- `cl_GT_InputTransport.cpp:143,146` parse windows clip fractional
  temperatures in shipped trans.inp (63.7→63.0, 126.1→126.0; 20+ lines).
- Unknown/miscased species labels produce a silent zero/NaN RefGas (factory
  never errors); interaction-viscosity existence is ASSERT-only (gone in
  release).
- `expand()` fills the d(energy)/dP Jacobian slot with the `hdep` *value*
  where `total()` uses `dhdp` (`cl_Gas.cpp:2542,2661`).
- **Methane `compute_phir_tt` Gaussian block is linearized** — the Hydrogen
  sibling carries the correct squared-bracket form; methane cv/cp/w wrong
  near the critical region. Existing tests only sweep p = 1 MPa and cannot
  see it (`cl_GM_EoS_Methane.cpp:316-324`).
- Cubic liquid-root selection indexes `mWorkZ(1)/(2)` without checking how
  many real roots cardano returned; `EoS_Cubic::T(p,v)` silently returns
  NaN which `Gas::T` caches; `T_vap()` Newton has no iteration cap.
- Scrambled cache stamps defeat the exp() caches in the Hydrogen/Methane/
  Oxygen Helmholtz classes (recompute every call + latent stale-skip).
- `Helmholtz::set_reference_point` subtracts pv twice for mU0 (u offset by
  ≈ −R·Tref; probe test added).
- Suspected entropy reference double-count in `realgas_s` (+`mSref` on top
  of an absolute-anchored spline entropy; probe test added).
- Literature-blocked constants (papers not in the repo): Chung β 1.368 vs
  1.3168 (`fn_GT_idgas_lambda.cpp:46`), MC-SRK generalized c1 0.16054 vs
  ~1.6054 (`cl_GM_EoS_AlphaFunctionFactory.cpp:427`), undamped methane
  λ_cr with dead crossover helpers (`cl_GM_HelmholtzTransport_Methane.cpp`).
  **Need Christian / the cited papers.**

## Speedup opportunities (P2, confirmed)

- Flow routines (`total`, `expand`, `prandtl_meyer`, `shock`) construct
  Vector/Matrix/pivot temporaries per call; `compute_equilibrium` heap-copies
  its ΔX scratch; equilibrium scratch guard tests the wrong size so it can
  reallocate every call.
- `Helmholtz::v()` evaluates `is_liquid` twice per call (4 pow + exp each)
  and uses constant under-relaxation → linear convergence in the hottest
  function of the class.
- Defeated exp() caches (stamp bugs above) — fixing them is both a
  correctness and a speed win.
- Parsers: range-for by value copies every buffer line; factory opens all
  four data files per MPI rank (`OPEN_RDONLY` instead of the existing
  rank-0+broadcast `OPEN_RDONLY_PARALLEL`).
- `std::pow` with small constant/integer exponents in per-call polynomial
  paths; POLY-mode property calls do a linear scan over polynomials
  (mitigated: factory leaves gases in SPLINE mode).

## Design questions routed to Christian (not settled by vote)

- Reference-pressure convention of the real-gas caloric assembly — already
  documented as unsettled in the big comment above `Gas::realgas_cp`
  (`cl_Gas.cpp:1008-1045`) and in `todo/gasmodels_open_source_migration.md`;
  the `mSref` probe belongs to this complex.
- MPI semantics of `remix_transport` (rank-0-only spline update vs
  `remix_heat`'s all-rank write).
- Missing k_ij binary-interaction support in the cubic mixing rule.
- Whether the silent unknown-species factory behavior is intended API.

## Test-suite changes (committed to working tree tonight)

New **active** tests (pass on current code, lock verified behavior):

- `tests/physics/gastables/cl_GT_HeatPoly.cpp` — FD consistency of the
  NASA-9 poly: dH/dT=Cp, dS/dT=Cp/T, dCp/dT, d²Cp/dT².
- `tests/physics/gastables/cl_GT_RefGas_Consistency.cpp` — POLY vs SPLINE
  crosscheck for Cp/H/S/µ/λ.
- `tests/physics/gasmodels/cl_GM_EoS_Cubic_CallOrder.cpp` — α = pβκ
  identity.
- `tests/physics/gasmodels/cl_GM_Gas_Consistency.cpp` — real-gas dh/dT=cp,
  ds/dT=cp/T (SRK CH4).
- `tests/physics/gasmodels/cl_GM_EoS_Helmholtz_CpConsistency.cpp` —
  supercritical dh/dT=cp for Hydrogen (control for the methane defect).

New **DISABLED_** regression gates / probes (each cites its audit finding;
enable as the corresponding defect is fixed):

- `DISABLED_HeatPoly_dSdT`, `DISABLED_RefGas_dSdT` (C1),
  `DISABLED_RefGas_href_units` (C2),
  `DISABLED_input_transport_fractional_range` (C12),
  `DISABLED_Cubic_AlphaAfterKappa` (M3),
  `DISABLED_Gas_Remix_Statevals` (M1/M2),
  `DISABLED_Gas_SoundSpeed_IdgasLimit` (M5 probe),
  `DISABLED_Gas_Entropy_RefState` (M8 probe),
  `DISABLED_Helmholtz_u_Consistency` (M9 probe),
  `DISABLED_Methane_Cp_FD_NearCritical` (M7).

Both `tests/physics/*/CMakeLists.txt` updated. **Not compiled here** — next
step is `make tests` (the executable gate); the four pre-existing DISABLED
tests (AlphaFunction, Cubic_State, Cubic_Departure, Gibbs) were left as
found.

## Fix round (2026-08-05 evening, approved by Christian)

All confirmed defects fixed except the literature-blocked constants, the
design questions, and the judgment/perf items — full fix list in
`tmp/ai_exchange/review_gas_fixes_scope.md`; self-review, jury audit and
reconciliation in `review_gas_fixes_audit.md`. Highlights:

- P0 set closed: dSdT Rm scaling, h_ref units, remix/Statevals + cubic
  a(T) cache invalidation, alpha/ALPHA statebit, realgas_c cp factor.
- Methane `compute_phir_tt` Gaussian block rederived from the module's own
  `phir_t` convention — independently re-derived and confirmed by Grok in
  the post-fix jury; algebraically identical to the Hydrogen form.
- Parser hardening: new `fn_GT_parse.hpp` (std::from_chars — locale-safe,
  non-throwing, BELFEM_ERROR with context); all 24 stod/stoi sites
  replaced; correct A2+F6.2 composition windows (validated against the
  O2+ record with its `E -1.00` electron pair); 9-wide trans.inp
  temperature windows; buffer-bounds and line-length guards throughout;
  unknown species / empty data path / missing interaction pair are now
  hard errors.
- Consumer sweep before fixing: combustion holds `Gas&` references (the
  post-fix Codex claim of a copy-ctor compile break was REFUTED against
  cl_CN_Injector.hpp:21); archive is not in the build; no combustion or
  archive caller of Gas::T(p,v).
- Post-fix jury residuals fixed same round: Helmholtz-transport downcast
  now checked BEFORE the cast; four more short-line guards; T_from_h
  out-of-range bracket error.
- Routed to Christian (unchanged): shock() failure-path output semantics
  (total state vs last Newton iterate), remix_mass(aRemixHeat=false)
  covolume staleness for realgas consumers, cryo Cp-minimum bisection,
  glue acceptance short-circuits (>= 500 / >= 4000), is_liquid below
  T_triple, Oxygen phi0 commented term, flow-routine scratch members +
  Helmholtz v() relaxation (perf pass), literature constants.
- 9 regression gates enabled (were DISABLED); only
  `DISABLED_Gas_Entropy_RefState` (mSref probe, M8) stays disabled.

## itlr_combustor momentum question (archive, analysis only)

Christian's suspicion: the channel model's momentum RHS is wrong for
dA/dx ≠ 0 (Shapiro path works, channel path did not, in the thesis).
Finding: the archive channel path (`archive/channel/cl_CH_ChannelODE.cpp`,
`compute_channel_ode` + `compute_jacobi_idgas/realgas`) is **correct as it
stands** — symbolic elimination reproduces Shapiro's influence
coefficients exactly (ξ = −(b_mass + b_momentum) term by term, friction
4τw/(Dh·p) = 2k·cf·Ma²/Dh), the frictionless limit yields
(Ma²−1)·u'/u = A'/A, and the realgas Jacobian satisfies
cp − cv = T·v·α²/κ. Key physics: in the NON-conservative (differential)
momentum form the wall-pressure force p·dA/dx cancels identically —
dA/dx belongs ONLY in the continuity row. The classic way to get it wrong
is a mixed derivation: conservative flux d(pA + ρu²A)/dx on the left,
differential matrix on the right, dropping (or double-counting) p·dA/dx —
which produces errors proportional to dA/dx and none in cylindrical
sections, exactly the thesis symptom. Two genuine archive defects found
(NOT fixed — archive out of scope): `cl_CH_ChannelODE.cpp:304`
`cf = 2·v·v·τw/u²` has an extra v (dimensional; cf = 2·v·τw/u²) — this
makes the SHAPIRO path's friction wrong by a factor v (≈ 7 at the
scramjet state 0.43 bar / 1100 K), and the wall-heat term enters the
Shapiro path's η only inside `if( mCombust )` (:325-333), so a
non-combusting heated run drops wall heat entirely.

## Status

- [x] Pre-registration (both modules, frozen before dispatch)
- [x] Jury rounds (Codex + Grok, blind, parallel)
- [x] Verification pass + reconciliation tables
- [x] Test-suite additions (5 active, 10 DISABLED gates/probes)
- [x] Fix round (approved): all confirmed P0/P1 + mechanical P2 fixed;
      post-fix jury round run, residuals fixed, reconciliation appended
- [x] 9 regression gates enabled
- [x] itlr_combustor momentum analysis (archive; report only)
- [ ] Executable gate: build + `make tests` (Christian; needs GCC ≥ 11
      for std::from_chars<double>)
- [ ] Run the remaining probe (entropy ref state, M8) via
      `--gtest_also_run_disabled_tests`
- [x] Literature checks: Chung β (confirmed typo, 1.3168), MC-SRK c1 (confirmed
      typo, 1.60539), PR-78 κ (wrong digits, dead code), Friend λ_cr (χ and
      crossover F never applied — live defect; coefficient tables still blocked
      on the Friend 1989 paper) — see `dl20260806_gas_literature_verification.md`
- [ ] Design calls: shock failure path, remix_mass covolume, k_ij,
      reference-pressure convention, MPI remix_transport
