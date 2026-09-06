# Unify the Kohler resistivity pointer family + user callback on (T, B, beta)

**Date:** 2026-07-13
**Purpose:** Eliminate the argument-order split left by the 2026-07-11 minimal
fix — bring the whole `mFunctionRhoKohler` family AND the UserDefined user
callback onto one T-first `(T, B, beta)` order.
**Module:** physics/materials
**AIs:** Claude (edits), Grok + Codex (read-only verify; thread
`tmp/ai_exchange/kohler_pointer_unify_audit.md`)

## Why

The 2026-07-11 pass made the PUBLIC `rho`/`lambda` T-first but, as a minimal fix
for the broken Metal table path, left the private Kohler helpers +
`mFunctionRhoKohler` pointer T-last `(normB, angle, T)`. Christian judged that
split dangerous (all-`real` args → silent swaps; this area had already produced
three such bugs). Decision: one order everywhere.

## Change (all in src/physics/materials/)

Because every body references its params by NAME, each parameter LIST was
reordered to put `T` first and the bodies were left untouched; only the two
wrapper call-sites actually change behaviour.

- Pointer decls `mFunctionRhoKohler`/`mFunctiondRhoKohlerdT` (cl_Material.hpp:287-288)
  → `(T, normB, angle)`. Pointer TYPE unchanged (arity-only), so the
  `static_cast` assignments (Metal.cpp:40/41/695/696) are untouched.
- Wrapper calls (Metal.hpp:449/467): `mFunctionRhoKohler( B, beta, T )` → `( T, B, beta )`.
- `rho_kohler`/`drhodT_kohler` defs (Metal.cpp:566/582), `rho_table`/`drhodT_table`
  defs (Metal.hpp:453/471), their override decls (Metal.hpp:368/371), and the base
  virtual decls (cl_Material.hpp:1199-1208) → `(T, B, beta)`.
- User callback (Christian follow-up): `mUserRhoFunction`/`mUserLambdaFunction`
  now `(this, T, normB, angle)` (UserDefined.hpp:309/323); registration
  `set_user_defined_function` validation enforces `(T, normB, angleBxJ)`
  (UserDefined.cpp:141-156); 3-arg `lambda_custom` signature + base decls
  `lambda_custom`/`lambda_table` aligned; Doxygen example + member-doc comments
  updated. No in-tree 3-arg user rho/lambda fn exists (example uses T-only
  MatFunc1), so nothing in-tree breaks.

## Verification

3-AI unanimous PASS (high). Grok + Codex + Claude all traced
`Metal::rho(T,B,beta)` → pointer → `rho_kohler`/`rho_table` (and drhodT twins),
each reading pos1=T, pos2=field, pos3=angle; static_casts type-valid;
registration validation matches the callback invocation; no remaining old-order
Kohler/table caller in `src/` or `nonfree/`. Codex also confirmed the DB
*population* order matches (mesh axis x=T, y=log-field, z=beta —
`fn_create_database_mesh.hpp:22`, `Metal.cpp:708/845`), so the table path agrees
with how the database was built.

Codex caveat (pre-existing, NOT caused by this change): the UserDefined 3-arg
rho/lambda path is currently unreachable — public `rho/lambda(T,B,beta)` assert
PureMetal and 3-arg registration only calls `set_custom` (which wires one-arg
functions). So the field-dependent user callback is dead today; the flip makes it
correct for when it is wired. Also flagged: stale old-order `(B,angle,T)` still in
comments/docs (UserDefined.hpp:51, cl_Material.hpp:944, example_user_material.cpp:149,
materials_usage_guide.md:67, dof_manager_usage_guide.md:3228) — doc-lag, deferred.

## Deliberately NOT changed (separate subsystems — Christian to decide)

- Jc user callback `mUserJcFunction3` (cl_JcFunction_UserDefined.hpp) — `(normB, angle, T)`.
- `jc`/`n` registration case (UserDefined.cpp:173-181) — `(normB, angleNxB, T)`.
- `rho_powerlaw` / `drho_powerlaw_dJ` (HTS E-J law) — `(normJ, normB, angle, T)`.

## Build / commit

Not built (Christian runs builds). Not committed — awaiting Christian's go.
Diagnostics seen during editing were the language-server `armadillo`-not-found
cascade, not real errors.
