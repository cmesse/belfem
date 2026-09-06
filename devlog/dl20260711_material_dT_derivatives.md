# Temperature Derivatives for Material Properties (cp, lambda, rho)

**Date:** 2026-07-11
**Purpose:** Audit and fix Christian's new dT-derivative layer in the material framework, plus the underlying Database derivative shape functions
**Module:** physics/database, physics/materials
**AIs:** Claude (primary), Codex + Grok (independent audits, threads in `tmp/ai_exchange/database_deriv_shapefn.md` and `tmp/ai_exchange/material_dT_derivatives.md`)

## Part 1: Database::evaluate_derivx + deval_*dxi shape functions

Christian added `evaluate_derivx()` (2D/3D) and xi-derivative shape functions
`deval_quad4dxi` … `deval_hex64dxi` to `cl_Database.hpp`, copied from
`src/fem/interpolation/lagrange/`. Verification criterion:
`deval_*dxi[n] == d/dxi(eval_*[n])`, checked by symbolic differentiation (sympy),
independently confirmed by Codex and Grok (three-way exact agreement).

**Result:** QUAD4/QUAD16/HEX8/HEX27/HEX64 all correct. `deval_quad9dxi` had 3 bugs
(confidence: high, compiler-independent symbolic proof):
- N[4], N[6]: factor-4 too small (reused `c = xi*eta*0.25` from the value
  function; the lagrange reference uses `c = xi*eta`)
- N[8]: was the eta-derivative `2*eta*(xi2-1)`; correct is `2*xi*(eta2-1)`

Fixed by Christian. Claude removed a leftover unused `xi2` in `deval_hex27dxi`
(warning-as-error in some configs).

Chain rule in `evaluate_derivx` (`* 2 * inv_element_step(0)`) verified correct.
Note: `evaluate()`/`evaluate_derivx()` now share a malloc'd member buffer `mN`
and are no longer re-entrant — accepted, since BELFEM does not use OpenMP here.

## Part 2: dcpdT / dlambdadT / drhodT in the Material framework

Design: derivative function pointers (`mFunctiondCpdT`, `mFunctiondLambdadT`,
`mFunctiondRhodT`) parallel to the value pointers; spline properties use
`Spline::deval`, custom properties fall back to central finite differences
(`gFinDiffDeltaT`, `typedefs.hpp:66`).

Three-way audit (Claude/Codex/Grok, full agreement) found and Claude fixed:

| # | Issue | Fix |
|---|-------|-----|
| 1 | Name hiding: Metal/Alloy 3-arg overrides hid base 1-arg `drhodT`/`dlambdadT` (compile error, confirmed by focused build) | added `using Material::drhodT; using Material::dlambdadT;` to Metal and Alloy |
| 2 | `Metal::drhodT(B,beta,T)` declared `override`, never defined (vtable link error) | defined; dispatches via new twin pointer `mFunctiondRhoKohlerdT` |
| 3 | `drhodT_finite_difference` missing `inline` (ODR) and subtracted `cp(T-dT)` instead of `rho(T-dT)` | both fixed |
| 4 | No derivative counterpart for the `mFunctionRhoKohler` dual dispatch (kohler vs table) | `mFunctiondRhoKohlerdT` wired at both sites (`cl_Material_Metal.cpp` ctor + `populate_rho_database`). Christian upgraded the Kohler branch to a semi-analytic product rule `drhodT_kohler` (`cl_Material_Metal.cpp:582`): f'·g + f·g' with FD only on `rho_i_custom` and on `kohler(B,S,beta)` through S, plus a T<eps → 0 guard matching `rho_kohler`. Accuracy is acceptable since Kohler is always tabulated before FEM runtime. Final names: `drhodT_kohler` / `drhodT_table` (renamed from drho_kohlerdT/drho_tabledT; Claude completed the rename at the two pointer-wiring sites and the stub messages). Claude also fixed a typo in the new function (`g` evaluated at T+dT instead of T, biasing S and f) and removed a formally ODR-violating `inline` on the .cpp definition of the virtual. |
| 5 | Alloy overrode 3-arg `rho`/`lambda` but had no derivative overrides → silently dropped field dependence | added `Alloy::drhodT(B,beta,T)` (table chain rule, mirrors `Metal::drho_tabledT`) and `Alloy::dlambdadT(B,beta,T)` (Wiedemann-Franz quotient rule, mirrors Metal) |
| 6 | 4-arg HTS `dlambdadT` asserted `PureMetal`; value overload asserts `HTS` (debug abort for YBCO) | assert corrected to `HTS` |
| 7 | Base 3-arg `dlambdadT(T,B,beta)` argument order contradicted Metal's `(B,beta,T)` override (same C++ signature!) | **Resolved the opposite way after review (Christian):** the lambda family is uniformly `(T, B, beta)`; Metal and Alloy overrides (value AND derivative) were reordered to match the base. This fixed a REAL pre-existing bug: all call sites (`mt_thermal_h.cpp:83,244,401,505`) pass `(T, norm_b, beta)`, so the old `Metal::lambda(B,beta,T)` override was silently receiving T in the B slot. The rho family stays `(B, beta, T)`. |
| 8 | `dcpdT(T)`/`dlambdadT(T)` called null-able member pointer without the `have(property)` asserts their value twins have | asserts added |

Verified correct as-written (no change needed): pointer wiring at all three
assignment sites (const/custom in `cl_Material.cpp`, spline in
`cl_Material_SplineLookupTable.cpp`); quotient rule in `Metal::dlambdadT`;
`exp(y)*dy/dT` chain rule in `drho_tabledT` (rho stored as ln(rho));
FD helpers for cp/lambda; YBCO/UserDefined `*_custom` paths flow through the
FD fallback correctly.

## Known remaining gaps (documented, intentionally not fixed)

- `drho_tabledT` / `Alloy::drhodT` return the boundary slope for T outside
  `[Tmin,Tmax]` although the clamped value function is constant there
  (true derivative 0). Same for B outside the clamp. Low priority; may even
  help Newton.
- Pre-existing (broader than this patch): UserDefinedMaterial's 3-arg custom
  rho/lambda hooks are not reachable through the public value API; derivative
  API inherits that gap.
- Base value function `lambda(T,B,beta)` still has the (T,...) order while its
  Metal override reads (B,beta,T); only the *derivative* was reordered.

## Status

All fixes applied to working tree, not yet committed. Build/test handoff to
Christian (no build was run by Claude; Codex ran one focused materials build at
~18:46 to confirm the name-hiding error).
