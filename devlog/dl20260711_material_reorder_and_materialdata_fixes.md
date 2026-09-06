# Material Arg-Order Unification (T,B,beta) + MaterialData Mechanical Fixes

**Date:** 2026-07-11
**Purpose:** Complete Christian's reorder of the magneto material accessors to
T-first `(T, B, beta)` and fix the mechanical defects it (and the WIP
MaterialData class) left behind.
**Module:** physics/materials, fem/kernel (MaterialData), fem/thermal, fem/maxwell
**AIs:** Claude (primary + edits), Codex + Grok (independent third audit; thread
`tmp/ai_exchange/material_reorder_and_materialdata_audit.md`)

## Context

Christian unified the magneto accessors on **T-first** `(T, B, beta)` (option (c)
of `todo/rho_lambda_argument_convention.md`) — previously `rho` was T-last
`(normB, angle, T)`, `lambda` already T-first. All args are `real`, so stale
call sites compile silently and swap T↔|B|.

## Audit (3 voices, all high-confidence, cross-verified at file:line)

The reorder was applied INCONSISTENTLY. Confirmed stale sites:

- **mt_thermal_h.cpp rho callers 129/329/409/513** — only the variable was
  renamed (`normB`→`norm_b`); order left T-last while lambda callers in the same
  file were T-first. Core regression (thermal resistivity at T=|B|).
- **mt_maxwell_h.cpp:75** — `rho(gTbulk, beta, norm_b)`, B/beta swapped vs the
  correct sibling :145; **:2172, :2456** — untouched old order.
- **Metal table path** — `rho_table`/`drhodT_table` were reordered to `(T,B,beta)`
  but are reached via `mFunctionRhoKohler` whose contract is `(normB, angle, T)`
  and which `rho_kohler` matches → the DB path (Cu/Ag) got theta=clamp(|B|).
- **Metal & Alloy lambda/dlambdadT internals** — quotient-rule bodies call
  `rho/drhodT(B, beta, T)` against the new T-first wrapper → scrambled.
- Base `Material::rho/drhodT` param NAMES swapped vs decl (harmless; body uses T).
- main.cpp driver + usage-guide docs stale.

## Fixes applied

- Reordered the stale `rho`/`drhodT` callers to `(T, B, beta)`: mt_thermal_h
  (129/329/409/513), mt_maxwell_h (75/2172/2456), Metal.hpp + Alloy.hpp internal
  lambda/dlambdadT calls, base `Material::rho/drhodT` param names, main.cpp 88/379.
- **Chose the pointer-contract fix for the table path:** reverted
  `Metal::rho_table`/`drhodT_table` parameter lists to `(B, beta, T)` (bodies
  unchanged) so they match `rho_kohler` and the `mFunctionRhoKohler` contract.
  Only the PUBLIC `rho`/`lambda` are T-first; the private Kohler helpers stay
  `(normB, angle, T)`.
- **MaterialData (WIP, unwired) mechanical fixes:** `compute_lambda`/
  `compute_dlambdadT` were passing `compute_T(aIndex)` (a real) to the
  index-taking `mFunLambda`/`mFundLambdadT` → now pass `aIndex` (the pointee
  computes T itself). Fixed `mX, mY, mY` → `mX, mY, mZ` in the power-law bulk
  defect rho/drho-dJ.

## Deferred (design-level; NOT touched — need Christian + a build)

MaterialData is not constructed anywhere yet. Left for a supervised pass
(recorded in `todo/thermal_matrices_cleanup_and_newton_plan.md` and the audit
thread): HTS `mFunLambda`/`mFundLambdadT` unassigned + no `compute_lambda_hts`;
`drho/dJ` has no public dispatcher for Newton; per-element peer thermal relink;
memoization keyed only on `aIndex` (no invalidation on `Calculator::link`); ctor
`constant_property(mu)` now assert-exposed by the is_constant() fix. Plus the
usage-guide doc pass (old `(B, angle, T)` order).

## Build

Not built here (Christian runs builds). All edits are argument-order swaps and
type-neutral; the changed MaterialData paths are unwired.
