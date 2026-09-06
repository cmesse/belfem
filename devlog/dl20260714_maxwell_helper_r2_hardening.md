# Devlog 2026-07-14 — MaxwellData R2 Hardening (R2a-R2h) + T10 Material-Flag Audit

**Date:** 2026-07-14
**Topic:** Implementation of plan steps R2a-R2h of `todo/maxwell_kernel_collapse_plan.md`
(helper hardening to equivalence-grade), executed by Claude/Fable with Christian's approval;
Codex + Grok audits launched same session.
**AIs involved:** Claude/Fable (implementation), Codex + Grok (audit, pending)
**Claude Confidence:** high (all edits verified by grep/read-back; build handed to Christian)
**Literature References:** consistent/algorithmic tangent — Belytschko §5.9, Zienkiewicz &
Taylor Vol 2 §4.4 (Simo & Taylor 1985); power law — Messe et al. 2023 (paper1) §2.6

## Summary

All remaining R2 items implemented (R2d partially — see below). One NEW latent bug found and
fixed while implementing R2e: `gRhoMin`/`gRhoMax` were never initialized in
`Communicator::set_globals()` (zero-initialized globals) — any executable that did not set
them manually would have clamped every resistivity into `[0,0]` under the new
`compute_rho` clamp. Defaults `0 / BELFEM_REAL_MAX` added (clamp = no-op unless an
executable narrows the window, as hphirun/hphiTrun do).

## Changes Made

- **R2a (D3)** — `cl_FEM_Calculator.cpp` ctor: β-convention comment + loud guard
  `BELFEM_ERROR( !(tIsHTS && tLambdaFieldDependent), … )`; `beta_dummy()` private accessor
  with `BELFEM_ASSERT( mBeta == π/2 )` now used at all 8 bulk-HTS material-call sites in
  `cl_FEM_Calculator.hpp` (replaces the bare `mBeta` dummy arguments).
- **R2b (D6)** — "deliberate self-assignment" comment at `compute_b`.
- **R2c (D4/O1)** — `Material::jc_eval()` / `n_eval()` protected helpers (declared
  `cl_Material.hpp`, defined `powerlaws.hpp`): route on `depends_on(JcParameter::T)` to the
  overridden `eval` arity, constants-fallback when no function. All 8 T-bearing
  powerlaw/piecewise rho+drho overloads now delegate; the four full-dependency
  "wrong powerlaw" assert blocks and the nullptr asserts they guarded are gone (the 4
  remaining asserts are the reduced 3-arg overloads' `!depends_on(T)` — intentionally kept).
  **Verification (per subclass):** ModifiedKim declares normB+angleNxB, overrides
  `eval(B,angle)` ✓; Database declares +T, overrides `eval(B,angle,T)` ✓; UserDefined
  overrides BOTH arities ✓ → the routing always lands on an overridden function.
- **R2d (O2)** — `compute_T_fem` clamps into `[gTmin, mTmax]` (mTmax = material `T_max`
  when constant, else unbounded; cached at ctor) and sets `mTClamped`;
  `compute_dcpdT`/`compute_dlambdadT` return 0 while clamped (consistent tangent).
  `T_clamped()` accessor exposed. **Deferred:** the once-per-step converged-at-clamp
  diagnostic needs Controller-side convergence context — accessors are in place, wiring is
  a follow-up (plan R2d stays [◐]).
- **R2e (O3)** — `compute_rho` clamps into `[gRhoMin, gRhoMax]`, sets `mRhoClamped`;
  `compute_drhodj` forces rho current first and returns 0 while clamped. `rho_clamped()`
  accessor. Plus the `set_globals()` default fix described above (`cl_Communicator.cpp`).
- **R2f (O4)** — `mDensity` cached at ctor (`density(gTroom)` → `ref_density()` fallback →
  NaN + asserting `density()` getter, so magnetics-only models with density-less materials
  do not fail at setup); physics-trap comment at the member.
- **R2g** — two new sections in `src/physics/materials/doc/materials_usage_guide.md`
  (+ revision row 1.1): "Assembly Contract: the Full-Signature Policy for Jc / n" and
  "Density and the Undeformed Mesh (Physics Trap!)".
- Also: `reset()` clears both flags; `<algorithm>`/`<cmath>` includes added to
  `cl_FEM_Calculator.hpp`.

## R2h — T10 Material-Flag Audit Table

Classification the helper uses: `tIsHTS = have(jc)`, `tIsMetal = !jc && depends(rho,normB)`,
`tLambdaFieldDependent = depends(lambda,normB)` (`cl_FEM_Calculator.cpp:85-97`).

| Material (file) | MaterialType | have(jc) | depends(rho,normB) | depends(lambda,normB) | Helper route | Legacy dispatch | Verdict |
|---|---|---|---|---|---|---|---|
| Copper | PureMetal | no | **yes** (Kohler, `cl_Material_Copper.cpp:602-608`) | **yes** | rho_metal + lambda_metal | h_metal / T_h_metal | ✓ match |
| Silver | PureMetal | no | **yes** (`cl_Material_Silver.cpp:487-493`) | **yes** | rho_metal + lambda_metal | h_metal / T_h_metal | ✓ match |
| Indium, Lead, WhiteTin, Iron | PureMetal | no | **verify per material** — base `Metal` sets only T-deps (`cl_Material_Metal.cpp:314-321`); Kohler deps are set where magnetoresistance data exists | ditto | rho_bulk if no normB dep | h_metal → `rho(T,B,β)` via `mFunctionRhoKohler` (`cl_Material_Metal.hpp:452-455`) | ⚠ helper calls `rho(T)`, legacy calls `rho(T,B,β)`→`rho_kohler`; equivalent ONLY if the Kohler pointer degrades to the T-only value without Kohler data — **shadow-compare will decide (R7)** |
| HastelloyC276 | LookupAlloy | no | no | no | rho_bulk + lambda_bulk | h_alloy / T_h_alloy | ✓ match |
| YBCO | HTS | **input-dependent**: jc set by MaterialFactory from input.conf (constant or lookup, `cl_MaterialFactory.cpp:117,186`), NOT by the class | no | no (λ is T-only, `cl_Material_YBCO.cpp:107-108`) | HTS powerlaw route | h_hts | ✓ match when jc configured (always, in practice); an HTS block without jc in the input would fall to rho_bulk — flag at setup? (note for R5) |
| Magnesia | NonMetal (Buffer) | no | no | no | rho_bulk — but `compute_rho` is never called on Buffer blocks (phi formulation) | phi kernels | ✓ inert; ctor now touches only density/T_max — both defined |
| Alloy (user, `cl_Material_Alloy.cpp:38,53`) | **PureMetal (!)** | no | from lookup data | from lookup data | data-dependent | h_metal (type-based!) | ⚠ user alloys are TYPED PureMetal — legacy sends them to h_metal (`rho(T,B,β)`) regardless of their actual deps; helper follows the deps. Same class of mismatch as the plain-metal row; shadow compare decides. **Question for Christian: is `MaterialType::PureMetal` for user Alloys intentional?** |
| UserDefined | UserDefined | config-dependent (`cl_Material_UserDefined.cpp:101,195`) | config | config | dependency-routed (O1) | generic h/h_ts branch on same deps | ✓ by construction after R2c |

**Bottom line:** built-in HTS and lookup-alloy families match exactly; Cu/Ag match exactly;
the risk rows are non-Kohler PureMetals and user Alloys (typed PureMetal), where legacy
dispatches on the *type* and the helper on the *dependency flags*. Both rows are exactly
what the R3 shadow harness exists to adjudicate — no action until it reports.

## Audit Outcome (same session)

- **Codex: "no blocking R2a-R2h defect found"** (high confidence) — verified overload
  delegation + const-correctness, clamp memoization coherence, `beta_dummy()` site
  coverage, gRho-default safety against all in-tree usage, density fallback soundness,
  naming/style. Two caveats, both fixed on the spot: `compute_T_const` now clamps + sets
  `mTClamped` like the FEM path (O2 truly unified), and the reduced 3-arg derivative
  overloads gained the symmetric `!depends_on(T)` asserts (pre-existing asymmetry).
- **Grok: unavailable** — two truncated runs (narration only, exit 0); its four audit
  questions were answered by direct code check instead (no mixed derivative consumers
  exist; table layer already clamped T internally, e.g. `Metal::rho_table`; jc/n routing
  independent per function; UB covered by Codex). Feedback memory updated: stop retrying
  after two truncations.

## R1 Smoke Outcome (same day)

First smoke run SIGSEGV'd (hphirun, jump to 0x0 from `Calculator::link:1300`) → **D12**:
the MaxwellData construction + `mFunLinkElement` selection sat below `allocate()`'s
"done if this is a block" early return (`cl_FEM_Calculator.cpp:833,940`) — dead code for
every BLOCK calculator since the wiring landed, missed by all three read-only audit
passes. Fixed by hoisting the section to right after `mIsAllocated = true` (work vectors
not a precondition — `link_vector()` creates on demand) plus a default initializer on
`mFunLinkElement`. **Rebuild → smoke PASSES (Christian). R1 closed.**

## R3/R4 Cut + R5 Implementation (same day, after the smoke pass)

Christian cut R3 (shadow harness) and R4 (frozen baselines): the git commit is the
baseline and rollback — verification is run-based per §4.1 (revised). R5 then implemented
(Claude): `maxwell::h_calc` / `h_newton_calc` (`mt_maxwell_h.{hpp,cpp}`) and `fem::T_h_calc`
(`mt_thermal_h.{hpp,cpp}`) — the three collapsed kernels, delegating all material math to
`aCalc->maxwell()`, bodies mirroring the legacy expression forms literally, zero call
sites (dispatch still 100% legacy). `T_h_newton` deliberately not written (new physics —
`thermal_matrices_cleanup_and_newton_plan.md` owns it). Pending: build, then the R6
LookupAlloy flip (one family per commit).

## Open Questions

- R2d diagnostic wiring (Controller-side, deferred).
- The two ⚠ rows above (non-Kohler PureMetal, user-Alloy typing) — for Christian.
- Christian: rebuild (R2 landed after the R1 compile) + run the R1 smoke set.

## Files Updated

- src/fem/kernel/cl_FEM_Calculator.hpp (R2a/b/d/e/f + includes + accessors)
- src/fem/kernel/cl_FEM_Calculator.cpp (R2a guard, R2d T_max window, R2f density — ctor)
- src/physics/materials/cl_Material.hpp (jc_eval/n_eval declarations)
- src/physics/materials/powerlaws.hpp (R2c dependency routing, 8 overloads + 2 helpers)
- src/comm/cl_Communicator.cpp (gRhoMin/gRhoMax defaults — latent-bug fix)
- src/physics/materials/doc/materials_usage_guide.md (R2g, two sections + revision row)
- todo/maxwell_kernel_collapse_plan.md (boxes ticked, R1 build noted)
