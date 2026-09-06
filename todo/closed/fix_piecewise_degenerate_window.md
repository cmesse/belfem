# Fix the Piecewise Power-Law NaN in the Degenerate n-Window (tapestack3d DR)

> **CLOSED 2026-09-03** (todo/ currentness sweep, round 3): the Riva law is landed and wired; R8 (harden the legacy `piecewise` window) declined as low value; R7 deck gate owed. Status lines and checkboxes below are as they stood at closure and are not maintained.

**Date:** 2026-08-27
**Purpose:** Eliminate the NaN that killed `tapestack3d_coarse` at timestep 122: for table-driven
jc/n materials, the flux-flow Bézier machinery in the `rho_piecewise` family produces NaN when the
table's n(T,B,θ) approaches 1 (measured transition softening near T_crit), and the degeneracy guard
in the dT-leg returns a value that was already poisoned before the guard runs. Mechanism in one
sentence: hoist and harden the degenerate-window fallbacks so every member of the piecewise family
returns finite values for every (J, T, B, θ) the solver can produce mid-iterate.
**Module:** `src/physics/materials` (powerlaws.hpp); consumer context `src/fem/kernel` (Calculator)
**AIs involved:** Claude (diagnosis + plan), Codex (audit), Grok (third voice), Christian (model decision)
**Status:** IN PROGRESS — R2-R5 implemented 2026-08-27 (riva law + n-floor + wiring + input
contract), three touched TUs syntax-clean, Python FD/finiteness verification passed (0/97,825
nonfinite; dJ 2e-5, dT 1e-6 vs FD). Awaiting Christian's rebuild + R7 deck gate. Open: R5 Codex
language sweep, R6 in-build probe, R8 piecewise hardening, R9 devlog/citations.

> **Scope guards:**
> - The BDF5 rejection-cascade behaviour (predictor blowup over a 16:1 step-ratio spread) is OUT of
>   scope here — it is the *trigger*, not the defect, and belongs to a separate timestep-controller task.
> - `piecewise` keeps its current behaviour for existing decks. The Duron/Riva parallel model enters
>   as a NEW opt-in law `resistivity type : riva` beside `powerlaw` and `piecewise`
>   (decided 2026-08-27, Christian) — not as a change to `piecewise`.

## 1. Diagnosis (evidence: core dump of PID 2640332, 2026-08-27)

Verified from the post-mortem core (`coredumpctl` id 2640332, rank 3), confidence high:

- Crash site: `SolverData::assemble_jacobian`, `cl_FEM_DofMgr_SolverData.cpp:1615`, element 778962,
  block 72 (`PENTA6TS`, `DomainType::ThinShell`, material YBCO), magnetic assembly, iteration 3 of
  the 5th retry of timestep 122 (Δt = 1.1584 ms).
- The full 6×6 element jacobian is NaN. The Calculator's `MaxwellData` cache at the crash:
  `mT = 89.08 K`, `mNormB = 0.1477 T`, `mNormJ = 8.88e10`, `mBeta = 1.2995 rad`, `mRho`, `mdRhodJ`,
  `mdRhodB` all finite, **`mdRhodT = NaN`**. The NaN factory is
  `Material::drho_piecewise_dT(J, T, |B|, θ)` (`powerlaws.hpp:2148-2297`).
- Reproduction (Python transcription against the real `sp-ap.hdf5` grid, trilinear stand-in for the
  quadratic B-spline): at this field/angle the table gives n = 8.6 @ 88 K falling to n = 1.014 @
  90.7 K (values are stored as log10; n crosses 1 near ~90.7 K, jc collapses in the same window).
  For T ≈ 90.7 K and J ≳ 1.8e11 (blend branch), `j2 = j1·(ρn/ρ1)^{1/(n−1)}` overflows to inf,
  `logj2 = inf`, a and b infinite, `tDisc = b² + ac = inf − inf = NaN`, and the function returns NaN
  **through its own degeneracy guard** — branch `GUARD` confirmed in the reproduction.
- Root defect (code): the guard at `powerlaws.hpp:2284` tests `!(tDisc > tol)` (correctly NaN-aware)
  but its fallback `tFrozen` is computed at `:2213-2225`, upstream of the guard, through
  `tSqrt = pow(tDisc, 0.5)` and `tParam = (b + tSqrt)/a` — both NaN when tDisc is NaN or negative.
  The guard catches the *condition* and returns the *poisoned value*.
- Trigger (context, not defect): the timestep-122 rejection cascade degraded the BDF5 predictor
  (first thermal Picard residual per attempt: 3e-4 → 2e-2 → 1e-2 → 1.2e-2 → 1.32e+2 as Δt was cut
  18.5 → 1.16 ms); the wild iterate pushed layer T from 77 K to ~89-91 K at the tape edge. The last
  converged state (t = 4.55 s) is clean: whole T-field within 1.3 mK of 77 K, |J/Jc|max = 1.14.
  Non-reproducibility under lldb / warm restart is expected: the NaN needs a gauss point at
  T ≈ 90-91 K with J > j1 inside one wild rejected-step iterate; any trajectory perturbation
  misses the window.

Secondary findings (confidence high, verified by reading + Python):

- The same n→1⁺ overflow exists in `rho_piecewise` (`:735`, the residual) and
  `drho_piecewise_dJ` (`:1815`): their `n > 1` asserts pass at n = 1.014 while
  `(ρn/ρ1)^{1/(n−1)}` still overflows. The residual itself can NaN the same way.
- `drho_piecewise_dT` has no jc/n asserts at all, unlike its siblings.
- Pre-existing residual/tangent inconsistency: `drho_piecewise_dJ`'s power-law branch applies a
  parallel-resistor factor `1/(1+ρPL/ρn)²` (`:1791`) while the residual `rho_piecewise` returns the
  RAW power law below j1, and the dT-leg documents raw-consistency (`:2166-2168`). One of the three
  is wrong; the tangents must be consistent with the residual (see O3).
- A 1.4M-sample fuzz of the dJ-leg blend over n ∈ [3.2, 5.0] with self-consistent knots produced
  no NaN — the discriminant cancellation alone (n−1 ≈ mNff, D≈0) does not fail in floating point;
  the realized failure needs the n→1 overflow route. Negative result recorded to prevent re-audit.

## 2. Literature

- **Duron, Grilli, Dutoit, Stavrev 2004** (Physica C 401, "Modelling the E–J relation of high-Tc
  superconductors in an arbitrary current range"): the parallel combination
  ρ = ρ_SC·ρ_NS/(ρ_SC+ρ_NS) is the published standard for making the power law valid for arbitrary
  J. Needs no knots, no Bézier, no n > 1 precondition; exact ohmic limit ρ → ec/jc as n → 1.
- **Riva 2021** (EPFL thesis 8754, `./tmp/EPFL_TH8754.pdf`): uses Duron's parallel form for both his
  models (Eq. 5.4); measures the overcritical regime 77-90 K and finds the resistivity slope
  *softens* markedly vs the power law — the physics our sp-ap table encodes as n(T,B,θ) → 1 near
  T_crit. His ρ_ηβ model (Eq. 5.1) is an alternative heuristic with fitted η(T), β(T); §5.1 carries
  the continuity analysis at I → 0⁺, Ic → 0⁺.
- Quantitative A/B (Python, healthy 77 K constants n = 19.25, jc = 4.1e10, ρn = 1.4e-7):
  parallel and current Bézier agree to 4e-8 relative for J ≤ 1.05 jc; differ up to 209× in
  J/jc ∈ [1.4, 55] (Bézier stretches the transition to ~76 jc, parallel saturates by ~3 jc).

## 3. Plan steps

The main line is the new law (R2-R7); the latent piecewise guard defect is hardened separately (R8)
because piecewise remains reachable by every existing deck.

- [x] **R1** — Pre-register this plan in `tmp/ai_exchange/riva_resistivity_law.md` and
  dispatch the Codex + Grok plan audits (three-vendor rule). Done 2026-08-27; both audits landed,
  reconciled below (§3a). Post-Grok `git status` clean w.r.t. Grok (five thin-shell/mesh files
  modified 15:30-15:58 predate the Grok run — parallel session, not a breach).

### §3a Audit reconciliation (2026-08-27, both voices in)

Verdicts on the pre-registered claims — full threads in the exchange file:
- Claim 2 (parallel weights) **CONFIRMED** by both; already implemented in the powerlaw legs.
- Claim 5 first anchor (4e-8 sub-critical agreement) **CONFIRMED** independently by both;
  209× retained as probe-derived only.
- Claim 1 (totality for free) **REFUTED** (both): defect callback `D·jc` unrestricted, `pow(10,x)`
  underflow, constants/custom fallbacks. Totalization for `jc_eff ≤ 0`: return ρn / dρn (fully
  normal), as an expected mid-iterate state — finite fallback, NOT `BELFEM_ERROR` (error-tier
  rule).
- Claim 3 (n-floor) **CONFIRMED as decision, refuted as written**: `dn_eval_dT/dB` must return 0
  while the floor binds; `n_eff == 1` needs a closed form (ρPL = ec/jc, no pow — log-space
  `0·(−∞)` trap at J = 0); J→0 limit at n_eff = 1 is (ec/jc)∥ρn, which near T_crit IS ρn; C0 not
  C1 at the kink (accepted).
- Claim 4 (T_crit continuity automatic) **REFUTED** (both): table Tmax (92) vs T_crit (92.5)
  mismatch means the `T > T_crit` cut can be a jump if jc(Tmax) is not small enough — a property
  of each hdf5, to be VERIFIED in R6 per table; `T == T_crit` boundary convention documented.
- My original R8 fallback **REFUTED** (both, and re-derived by me): the n→1⁺ overflow limit of the
  Bézier is t → 0 (ohmic ρ1), not t → 1 (ρn). R8 rewritten below.
- **Architectural finding (both voices): the existing `powerlaw` law IS Duron's parallel form**
  (`powerlaws.hpp:137-144, :252, :272`; its own docs cite Rhyner 1993, not Duron). `riva` as
  specified = hardened powerlaw. Kept as a separate enum per Christian's decision, but R5 must
  state how riva differs from powerlaw (guards + n-floor + total derivatives), or operators will
  assume a different E-J model. Existing powerlaw keeps its own derivative-leg overflow NaN
  (`drho_powerlaw_dT:2113` infs at ρPL = inf; reproduced in the Python sweep at a table-corner
  n = 86, needs B ≥ 5 T — unreachable in the tapestack deck) — harden or document (O4).
- **Diagnosis nuance (Grok, ~65%):** the core signature (mRho, mdRhodJ finite; mdRhodT NaN) fits
  the n−1 = mNff discriminant route at T ≈ 89 K at least as well as the n→1 overflow at 90.7 K —
  the per-index value caches mean the three cached scalars need not come from the same gauss
  point, so the core cannot arbitrate. Both mechanisms are real, reachable defects; R8 covers
  both regardless, so the fix does not depend on which fired.
- [x] **R2** (after R1, audits reconciled) — Material-side plumbing: replace the `mUsePiecewise`
  bool with a `ResistivityLaw` enum { PowerLaw, Piecewise, Riva } + `set_resistivity_law()` /
  `resistivity_law()`; keep `use_piecewise()` as a shim during the transition, **bound to
  `law == Piecewise` ONLY** (a `law != PowerLaw` shim would silently wire riva to the Bézier
  wrappers — Grok). Parse `resistivity type : riva` at BOTH `cl_MaterialFactory.cpp` sites
  (`:147` builtin, `:204` custom; the `else` is a hard `BELFEM_ERROR`, so a missed site is a
  loud fail, not a silent one). O1 floor: `std::max(n, 1.0)` in `n_eval` (`powerlaws.hpp:73`),
  `dn_eval_dT` (`:117`) and `dn_eval_dB` (`:96`) return 0 while the floor binds.
- [x] **R3** (after R2; implemented 2026-08-27, assembly-path signatures only per Grok) — The `rho_riva` family in `powerlaws.hpp`:
  ρ_riva = ρPL·ρn/(ρPL+ρn) (Duron et al. 2004 Eq. as used in Riva 2021 Eq. 5.4), with
  `rho_riva`, `drho_riva_dJ`, `drho_riva_dT`, `drho_riva_dB` and the same overload ladder as the
  existing families (constant, T-only, (B,θ), (T,B,θ), each ± defect). Requirements:
  - overflow-stable evaluation (branch on ρPL ≷ ρn: `ρn/(1+ρn/ρPL)` vs `ρPL/(1+ρPL/ρn)`, ρPL
    exponent computed in log space with an early ρn return above a log-threshold);
  - analytic derivatives through the parallel weights w = ρn/(ρPL+ρn):
    dρ/dx = w²·dρPL/dx + (1−w)²·dρn/dx, with dρPL/dJ = ρPL(n−1)/J,
    dlnρPL/dT = −(n/jc)·djcdT + ln(J/jc)·dndT, dlnρPL/dB analogous via djcdB/dndB;
  - guards (audit-revised): `jc_eff ≤ 0` or nonfinite (defect callback, underflow) → return
    ρn / dρn·dT / 0 / 0 (fully normal, finite — expected mid-iterate state, never
    `BELFEM_ERROR`); residual overflow early-out (ρPL above a log-threshold → ρn) with MATCHING
    derivative early-outs (dρ/dJ → 0, dρ/dB → 0, dρ/dT → dρn/dT) — without these the sweep
    fails in exactly the window the law exists for; `n_eff == 1` closed form ρPL = ec/jc (no
    pow — the log-space route is 0·(−∞) at J = 0); `J < BELFEM_EPSILON` guard on the dJ leg only
    (residual needs none at n_eff ≥ 1);
  - `T > T_crit` → exactly ρn(T) and dρn/dT; continuity at T_crit is NOT automatic (table Tmax
    vs T_crit mismatch) — verified per table in R6, convention at `T == T_crit` documented;
  - n-floor lives in R2 (central, `n_eval`), not here;
  - NO `n > 1` precondition; document in-code how riva differs from `powerlaw` (which is already
    the same parallel form: guards + total derivatives + n-floor are the delta) and cite Duron
    et al. 2004 where powerlaw currently cites only Rhyner 1993.
- [x] **R4** (after R3; 16 wrappers + 4 three-way dispatch sites, 2026-08-27) — Calculator wiring: `MaxwellData::compute_{rho,drhodj,drhodb,drhodT}_riva_{ts,bulk}{,_defect}`
  wrappers in `cl_FEM_Calculator.hpp` mirroring the piecewise pattern, and the third dispatch branch
  at the four selection sites (`cl_FEM_Calculator.cpp:286,303,382,399`). Sweep for any other
  `use_piecewise()` consumers (postprocessor EJ output included).
- [x] **R5** (with R2-R4, same session; both artifacts updated 2026-08-27; Codex language sweep applied same day — two flagged claims reconciled: the additive-vs-floor desync attribution and the T_crit C0 wording) — Input contract: `resistivity type` gains enum value
  `riva` in BOTH `doc/input_file_reference.md` (`:482` table) and `doc/input_schema.yaml`
  (`:817-822`; keep `case_sensitive: true` and default `power-law` unchanged), with the
  Duron/Riva citations, the "differs from piecewise above ~1.4 jc" caveat, and an explicit
  statement of how riva differs from `powerlaw` (same parallel model; riva adds totalization).
  Must NOT claim exact Riva Eq. 5.1/5.4 compliance: BELFEM keeps `mRhoMin` floor semantics
  (0 by default), not Riva's additive 1e-17 (which would reopen the 2026-08-10 tangent-desync,
  per both audits — deviation documented instead). Codex language sweep over the touched
  `input_file_reference.md` sections before the step closes.
- [◐] **R6** (after R3; Python transcription FD + finiteness sweep passed 2026-08-27, in-build probe against the compiled lib still owed) — Scratchpad probe: link the rebuilt material lib; FD-check all three
  derivative legs against the residual; finiteness sweep over the §1 reproduction grid
  (T ∈ [77, 93], J up to 1e13, sp-ap table values) plus the fuzz corpus; verify riva ≈ raw power
  law below 1.05 jc and riva → ρn above. (Probe policy: no separate audit round.)
- [ ] **R7** (after R6, Christian runs build + deck) — Deck gate: `tapestack3d_coarse` with
  `resistivity type : riva`, cold, 8 procs. Gate: must pass t = 4.55 s without aborting; rejections
  on physics are acceptable. Outputs needed per §5.
- [ ] **R8** (independent of R2-R7, after R1; audit-rewritten) — Harden the `piecewise` degenerate
  window in place (the shipped NaN remains reachable for piecewise decks). TWO distinct
  degeneracies, TWO distinct fallbacks (both audits; t→0 limit re-derived by Claude):
  - finite knots, `!(tDisc > tol)` (n−1 = mNff discriminant, J ≈ j3): the existing documented
    limit t = 1 → ρn / dadT — correct as-is, but HOIST above `tSqrt`/`tParam`/`rhoFF`/`tFrozen`
    so the fallback is computed clean, not through NaN;
  - `j2` nonfinite (n → 1⁺ overflow): the true Bézier limit is t → 0, NOT t = 1 — continue the
    power-law-regime expressions of the respective leg for J ≤ j3, ρn/dadT beyond;
  Coverage: the `j2` line exists in all 8 `rho_piecewise` overloads, all 8 `drho_piecewise_dJ`
  overloads, and both `drho_piecewise_dT` overloads (`:2148`, `:2300` defect twin with the same
  poison-then-guard order); `drho_piecewise_dB` skips the blend and needs no gate. Also disarm
  the `BELFEM_ERROR` on `|a|` in the dJ blend (`:1811-1813`) — a release-mode abort on NaN `a`;
  relax `n > 1.0` asserts to n ≥ 1 (the central floor makes exactly 1.0 reachable).
  Behaviour-preserving outside these two windows. O3 is explicitly NOT part of R8 (it is a
  behaviour change on healthy decks).
- [◐] **R9** — Devlog + citations. Citations DONE 2026-08-27: Rhyner 1993, Plummer & Evetts 1987,
  Duron et al. 2004 (DOI verified against the publisher record) and Riva 2021 added to
  `doc/literature_references.md` (new Tier 5: HTS Material Models); the new module doc
  `src/physics/materials/doc/resistivity_laws.md` covers all three laws and is registered in the
  module README (Codex language sweep dispatched). Devlog + exchange distillation still pending
  at session close.

## 3b. Defect tracker (code-audit round, 2026-08-27)

- [x] **D1** (Codex, HIGH) — subnormal jc passes the positivity gate, `ec/jc` overflows, dB/dT legs
  form `0·inf = NaN`. Fixed 2026-08-27: n ≤ 1 closed form returns `isfinite(rhoPL)`; infinite
  channel resistance takes the fully-normal branch. Verified by Python pathological suite.
- [x] **D2** (Codex, HIGH) — NaN n bypasses `std::max` and the `lg > 250` cap; inf n gives `inf·0`
  at J == jc. Fixed 2026-08-27: cap negated to the NaN-aware `!(lg <= 250.0)`.
- [x] **D3** (Claude, HIGH, found during D1/D2 verification) — `n/jc` and `J/jc` overflow alone for
  tiny finite jc → `inf·0 = NaN` even with zero derivatives. Fixed 2026-08-27: ratios grouped as
  `n·(djc/jc)`, weights as `(w·ρPL)·(w·dln)`, `log(J/jc)` → `log(J) − log(jc)`, all six legs.
- [x] **D4** (Codex, MEDIUM, docs) — "requires file" claim contradicted the factory, which
  correctly accepts constant jc/n. Fixed 2026-08-27 in both input-contract artifacts + code comment.

- [x] **D5** (Grok, HIGH) — dB residual/tangent split at n_eff = 1, J < eps: ohmic ρPL is
  J-independent but B-dependent; dB early-outed to 0 with no n-test. Fixed 2026-08-27 (n-gated
  early-out + ln guard, both overloads); FD-verified at the peak-effect point (jc ~ 1e3).
- [x] **D6** (Grok, MEDIUM) — ec ≤ 0 unguarded: NaN on the n > 1 path, corrupt parallel combo at
  n ≤ 1. Fixed 2026-08-27 (NaN-aware cap catches the log-of-negative; closed-form gate requires
  rhoPL ≥ 0). Follow-up: setup-time BELFEM_ERROR on ec ≤ 0 would fail loudly instead of silently
  going fully normal.
- [x] **D7** (Grok, docs) — both contract artifacts over-claimed "n → 1 falls back to fully
  normal"; actual behaviour is the finite ohmic closed form (ec/jc)∥ρn. Fixed 2026-08-27.
- **Incident:** D1/D2 were clobbered by a concurrent editor save of powerlaws.hpp between
  application and re-verification (parallel session active in the tree); re-applied and verified
  against the file text. Watch for this until the tree quiesces.

## 4. Open design questions

- **O1 — n < 1 handling.** The sp-ap table legitimately produces n < 1 above ~91 K (measured
  softening); for n < 1 the raw ρPL diverges as J → 0.
  **RESOLVED 2026-08-27 → decision (Christian): ohmic floor `std::max(n, 1.0)` enforced centrally
  in `Material::n_eval` (`powerlaws.hpp:73`), NOT per-law; `mRhoMin` stays untouched (no additive
  Riva 1e-17 — document the deviation from Riva Eq. 5.1/5.4 in R5).**
  Implementation corollaries (from the Codex audit, Claude concurs):
  - `dn_eval_dB` (`:96`) and `dn_eval_dT` (`:117`) must return exactly 0 while the floor binds
    (raw n ≤ 1) — differentiating the unfloored law would desynchronize the Newton tangent from
    the residual. Preferred: the gate lives inside `dn_eval_*` (single source of truth; costs one
    extra spline eval on the tangent path). Kink convention: floor binds at raw n ≤ 1, one-sided.
  - The central floor reaches `piecewise` too: at floored n == 1.0 exactly its
    `BELFEM_ASSERT( n > 1.0 )` aborts in debug and its `1/(n−1)` arithmetic still overflows in
    release — R8 must relax those asserts to n ≥ 1 and cover n == 1 in the overflow guard.
  - Behaviour change is confined to raw n < 1, which today is assert-abort or NaN territory for
    every law — strictly an improvement.
- **O2 — Model choice.** ~~Adopt the Duron parallel model as the piecewise default?~~
  **RESOLVED 2026-08-27 → decision (Christian): new opt-in law `resistivity type : riva` beside
  `powerlaw` and `piecewise`; piecewise behaviour untouched for existing decks.**
- **O3 — Which side of the PL-branch inconsistency is right?** dJ-leg applies the parallel factor
  (`:1791`), residual and dT-leg are raw. Consistency demands one convention; raw matches the
  residual today. For `riva` the question disappears (the parallel factor IS the model); for the
  legacy families, decide and align — as its OWN step, NOT inside R8 (behaviour change on healthy
  J ≤ j1 piecewise decks; both audits).
- **O4 — Existing `powerlaw` derivative legs are not total.** `drho_powerlaw_dT:2113-2114` (and
  the dJ/dB legs) go inf/inf → NaN when ρPL overflows; reproduced in the Python sweep at a
  table-corner artifact n = 86 (needs B ≥ 5 T and J ≳ 2e12 — unreachable in the tapestack deck's
  self-field, but real). Either share riva's hardened kernel with powerlaw (Grok's suggestion:
  one parallel-kernel with injected jc/n, fewer duplicated inlines) or document that powerlaw is
  not total. Christian to decide during R3 review. The n-spline's corner blowup (control values
  to 10^14 in log10 space) is a table-quality issue worth a separate look.

## 5. Run bookkeeping (announce-before-run rule)

For R6 I will need, until the gate is adjudicated: `out.txt`, `iv_results.csv`, the exodus frames
covering t ∈ [4.4, 4.7] s, and any core dump if it still aborts. The cold-run crash artifacts of
2026-08-27 are preserved in the session scratchpad (`tapestack3d_crash122/`, plus
`core.2640332`); the systemd journal holds the original core under id 2640332.
